// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// qd_hmatrix_operator.hpp — QD speed item 4: the Green's function G of
// qd_green_operator.hpp stored as an H-matrix (opt-in, --qd-green-compress
// EPS).  Design: qd-hmatrix-design-rev2 (Andy approved 20 Oct 2026).
//
// The stage product reads G once per stage at the memory-bandwidth limit
// (Vista jobs 1081359-1081361), so its time follows the bytes of G.  Blocks
// of G that couple well-separated patches of the fault have low numerical
// rank and are stored as U·Vᵀ; blocks between nearby patches stay dense.
// On the production G (p = 1, 1000 m) a tolerance of 1e-10 leaves 12.2 of
// 34.9 GiB, with a product error of 2.5e-15 (job 1082503).
//
//   1. Cluster tree: the fault DOFs sorted by position key (so the tree does
//      not depend on the partition), then bisected at the middle of the
//      longest side of the node positions' bounding box down to <= leaf
//      DOFs (as jobs/vista_qd_speed/qd_green_compress.py).
//   2. Block tree: a pair of clusters (t, s) is admissible when
//      max(diam t, diam s) <= eta·dist(t, s); other pairs are split, and a
//      pair with a leaf is a dense block.  Block (t, s) holds the ncomp·|t|
//      rows of t and the 2·|s| columns of s.  An admissible pair whose
//      smaller side (rows or columns) exceeds SEAS_QD_HMAT_SVD_SIDE_MAX is
//      split too, into its four sub-pairs: dgesdd's 32-bit LAPACK needs
//      about 4·min(m, n)² workspace entries, and above 2³¹ its own workspace
//      query overflows and returns a far too small size without an error
//      (p = 2, 1000 m: 41,904 × 27,936 blocks, Vista job 1084328).  No block
//      of the p = 1 production G reaches the cap (largest 20,952 × 13,968).
//   3. Compression: each block is gathered on one worker rank; an admissible
//      block is cut by an exact SVD (LAPACK dgesdd) to the smallest rank k
//      with ||A - A_k||_F <= eps·||A||_F, and kept dense when k·(m + n) >=
//      m·n.  Hence ||G - G_H||_F <= eps·||G||_F.
//   4. Layout: the rows stay where QDGreenOperator puts them (an even split
//      of the fault DOFs, in its global order g, over all ranks).  A rank
//      holds the rows of every block that fall in its range: dense rows,
//      or rows of U·Σ.  Vᵀ of a low-rank block is held by one rank on each
//      node that holds rows of the block, and the short vector w = Vᵀ·x
//      goes from there only to the ranks of that node (no w crosses nodes).
//   5. Product: x = δ_all in tree order; w = Vᵀ·x for the held Vᵀ; w to the
//      node's row ranks; y = dense rows·x + U rows·w.
//
// Cache file (SEAS-QD-HMAT-1), in tree order, so any rank count reads it:
// a 512 B header, the keys, the cluster ranges, T_D, the block table, then
// each block's data from a 4 KiB boundary (dense: m×n; low rank: U·Σ m×k,
// then Vᵀ k×n; row-major).  Blocks are written by their worker ranks, even
// ranks then odd ranks (see QDGreenOperator::WriteFile on Open MPI 5's
// OMPIO), and read back and compared before the run goes on.
//
// Linear-only, like item 3: it needs constant K and loading g(t) = t·ĝ, and
// stays off unless asked for.

#ifndef MFEM_SEAS_QD_HMATRIX_OPERATOR_HPP
#define MFEM_SEAS_QD_HMATRIX_OPERATOR_HPP

#include "mfem.hpp"

#ifdef MFEM_USE_MPI

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

// Largest smaller side of an admissible block cut by one SVD (see 2. above):
// 4·16384² = 2³⁰, half of what 32-bit LAPACK can address.  A test build may
// lower it to exercise the split on a small mesh.
#ifndef SEAS_QD_HMAT_SVD_SIDE_MAX
#define SEAS_QD_HMAT_SVD_SIDE_MAX 16384
#endif

#ifdef MFEM_USE_LAPACK
extern "C" void dgesdd_(char *jobz, int *m, int *n, double *a, int *lda,
                        double *s, double *u, int *ldu, double *vt,
                        int *ldvt, double *work, int *lwork, int *iwork,
                        int *info);
#endif

namespace mfem
{
namespace seas
{

class QDHMatrix
{
public:
   static constexpr int kKey = 6;

   /// @a nf fault DOFs with @a ncomp traction rows each; rank r holds the
   /// rows of fault DOFs [row_lo[r], row_lo[r + 1]) (QDGreenOperator's
   /// layout).  @a keys_all: 6 position keys per fault DOF in global order
   /// g, on every rank.  Builds the cluster and block trees.  Collective.
   QDHMatrix(MPI_Comm comm, long long nf, int ncomp,
             const std::vector<long long> &row_lo,
             const std::vector<real_t> &keys_all, int leaf, double eta)
      : comm_(comm), nf_(nf), ncomp_(ncomp), row_lo_(row_lo), leaf_(leaf),
        eta_(eta)
   {
      static_assert(std::is_same<real_t, double>::value,
                    "--qd-green-compress needs real_t = double (LAPACK "
                    "dgesdd)");
      MPI_Comm_size(comm_, &nranks_);
      MPI_Comm_rank(comm_, &rank_);
      MFEM_VERIFY(static_cast<int>(row_lo_.size()) == nranks_ + 1 &&
                  row_lo_[nranks_] == nf_, "QD H-matrix: bad row layout");
      MFEM_VERIFY(static_cast<long long>(keys_all.size()) == kKey * nf_,
                  "QD H-matrix: " << keys_all.size() << " key values, "
                  "expected " << kKey * nf_);
      MFEM_VERIFY(leaf_ >= 1 && eta_ > 0.0, "QD H-matrix: leaf " << leaf_
                  << ", eta " << eta_ << " (need leaf >= 1, eta > 0)");
      nb_ = static_cast<int>(row_lo_[rank_ + 1] - row_lo_[rank_]);
      SetupNodes();
      BuildTree(keys_all);
      BuildBlocks();
   }

   ~QDHMatrix()
   {
      int done = 0;
      MPI_Finalized(&done);
      if (!done && node_comm_ != MPI_COMM_NULL) { MPI_Comm_free(&node_comm_); }
   }
   QDHMatrix(const QDHMatrix &) = delete;
   QDHMatrix &operator=(const QDHMatrix &) = delete;

   /// Compress G from its dense rows on this rank (QDGreenOperator's layout:
   /// nb·ncomp rows × 2·nf columns, row-major, global columns 2·g + c).
   /// @a mem_node: bytes a node may spend at once on gathered blocks and
   /// their SVDs.  With @a keep, each worker keeps its blocks whole for
   /// WriteFile.  Collective.
   void Compress(const real_t *rows, double eps, double mem_node, bool keep)
   {
#ifndef MFEM_USE_LAPACK
      MFEM_ABORT("--qd-green-compress needs MFEM built with LAPACK "
                 "(MFEM_USE_LAPACK = YES): it calls dgesdd");
#endif
      MFEM_VERIFY(eps > 0.0 && eps < 1.0, "--qd-green-compress: eps "
                  << eps << " outside (0, 1)");
      eps_ = eps;
      kept_ = keep;
      const double t0 = MPI_Wtime();
      std::vector<int> order(blocks_.size());
      std::iota(order.begin(), order.end(), 0);
      // Biggest SVDs first, so the rounds end on small blocks
      std::stable_sort(order.begin(), order.end(), [&](int a, int b)
      {
         return SvdCost(a) > SvdCost(b);
      });
      std::vector<std::vector<int>> rounds;
      PlanRounds(order, mem_node, true, rounds);
      for (const auto &round : rounds)
      {
         RunRound(round, rows, nullptr, keep);
      }
      compress_seconds_ = MPI_Wtime() - t0;
      nrounds_ = static_cast<int>(rounds.size());
      FinishLayout();
   }

   /// y (this rank's nb·ncomp rows) = G_H·δ_all, δ_all = all 2·nf slips in
   /// global order.  Collective over the ranks of each node.
   void Mult(const real_t *delta_all, real_t *y) const
   {
      const size_t nrow = static_cast<size_t>(nb_) * ncomp_;
      for (long long p = 0; p < nf_; p++)
      {
         const long long g = perm_[p];
         xt_[2 * p] = delta_all[2 * g];
         xt_[2 * p + 1] = delta_all[2 * g + 1];
      }
      // w = Vᵀ·x for the Vᵀ held here
      for (const VPart &v : vparts_)
      {
         const real_t *x = xt_.data() + 2 * clusters_[blocks_[v.b].s].lo;
         real_t *w = wbuf_.data() + v.w;
         for (int l = 0; l < v.k; l++)
         {
            w[l] = Dot(v.vt.data() + static_cast<size_t>(l) * v.n, x, v.n);
         }
      }
      // w to the ranks of this node that hold rows of the block
      for (size_t i = 0; i < wsend_idx_.size(); i++)
      {
         wsend_[i] = wbuf_[wsend_idx_[i]];
      }
      MPI_Alltoallv(wsend_.data(), wscount_.data(), wsdispl_.data(),
                    MPITypeMap<real_t>::mpi_type, wrecv_.data(),
                    wrcount_.data(), wrdispl_.data(),
                    MPITypeMap<real_t>::mpi_type, node_comm_);
      for (size_t i = 0; i < wrecv_idx_.size(); i++)
      {
         wbuf_[wrecv_idx_[i]] = wrecv_[i];
      }
      std::fill(y, y + nrow, 0.0);
      for (const RowPart &r : rparts_)
      {
         const Block &bk = blocks_[r.b];
         const int nr = static_cast<int>(r.rows.size());
         if (r.k < 0)
         {
            const real_t *x = xt_.data() + 2 * clusters_[bk.s].lo;
            for (int i = 0; i < nr; i++)
            {
               y[r.rows[i]] += Dot(r.a.data() + static_cast<size_t>(i) * r.n,
                                   x, r.n);
            }
         }
         else
         {
            const real_t *w = wbuf_.data() + r.w;
            for (int i = 0; i < nr; i++)
            {
               y[r.rows[i]] += Dot(r.a.data() + static_cast<size_t>(i) * r.k,
                                   w, r.k);
            }
         }
      }
   }

   /// Write the H-matrix and T_D (this rank's rows, @a td_rows) to @a path.
   /// Needs Compress(..., keep = true).  Collective.
   void WriteFile(const std::string &path, const std::string &info,
                  const real_t *td_rows)
   {
      MFEM_VERIFY(kept_, "QD H-matrix: WriteFile needs Compress(keep)");
      std::vector<real_t> td_tree;
      GatherTD(td_rows, td_tree);
      std::vector<long long> off;
      const MPI_Offset end = DataOffsets(off);
      MPI_File fh;
      int err = MPI_File_open(comm_, path.c_str(),
                              MPI_MODE_CREATE | MPI_MODE_RDWR,
                              MPI_INFO_NULL, &fh);
      MFEM_VERIFY(err == MPI_SUCCESS, "QD H-matrix: cannot create '" << path
                  << "'");
      MPI_File_set_size(fh, end);
      std::vector<real_t> head;
      if (rank_ == 0) { HeadBytes(info, td_tree, head); }
      for (int phase = 0; phase < 3; phase++)
      {
         if (phase < 2 && rank_ % 2 == phase)
         {
            for (const Kept &kb : kept_blocks_)
            {
               WriteChunked(fh, off[kb.b], kb.data.data(), kb.data.size(),
                            path);
            }
         }
         if (phase == 2 && rank_ == 0)
         {
            WriteChunked(fh, 0, head.data(), head.size(), path);
         }
         MPI_File_sync(fh);
         MPI_Barrier(comm_);
         MPI_File_sync(fh);
      }
      // Read back and compare
      int bad = 0;
      std::vector<real_t> back;
      for (const Kept &kb : kept_blocks_)
      {
         back.resize(kb.data.size());
         ReadChunked(fh, off[kb.b], back.data(), back.size(), path);
         if (std::memcmp(back.data(), kb.data.data(),
                         back.size() * sizeof(real_t)) != 0) { bad = 1; }
      }
      if (rank_ == 0)
      {
         back.resize(head.size());
         ReadChunked(fh, 0, back.data(), back.size(), path);
         if (std::memcmp(back.data(), head.data(),
                         back.size() * sizeof(real_t)) != 0) { bad = 1; }
      }
      MPI_File_close(&fh);
      int bad_all = 0;
      MPI_Allreduce(&bad, &bad_all, 1, MPI_INT, MPI_MAX, comm_);
      if (bad_all && rank_ == 0) { MPI_File_delete(path.c_str(), MPI_INFO_NULL); }
      MFEM_VERIFY(!bad_all, "QD H-matrix: '" << path << "' read back "
                  "differently from what was written; the file was deleted");
      kept_blocks_.clear();
      kept_blocks_.shrink_to_fit();
      kept_ = false;
   }

   /// Read the H-matrix from @a path and this rank's rows of T_D into
   /// @a td_rows (nb·ncomp).  Returns false when the file does not exist;
   /// aborts when it was made for another mesh or other settings.
   /// Collective.
   bool ReadFile(const std::string &path, double eps, double mem_node,
                 real_t *td_rows, std::string &info)
   {
      MPI_File fh;
      int err = MPI_File_open(comm_, path.c_str(), MPI_MODE_RDONLY,
                              MPI_INFO_NULL, &fh);
      if (err != MPI_SUCCESS) { return false; }
      eps_ = eps;
      Header h;
      std::memset(&h, 0, sizeof(Header));
      if (rank_ == 0)
      {
         CheckIO(MPI_File_read_at(fh, 0, &h, sizeof(Header), MPI_BYTE,
                                  MPI_STATUS_IGNORE), path, "header");
      }
      MPI_Bcast(&h, sizeof(Header), MPI_BYTE, 0, comm_);
      const Header want = MakeHeader("");
      MFEM_VERIFY(std::memcmp(h.magic, want.magic, sizeof(h.magic)) == 0,
                  "QD H-matrix: '" << path << "' is not an H-matrix cache "
                  "file (SEAS-QD-HMAT-1); a dense --qd-green-cache file "
                  "belongs to --qd-green-cache, not --qd-green-hcache");
      MFEM_VERIFY(h.nf == want.nf && h.ncomp == want.ncomp &&
                  h.nclusters == want.nclusters &&
                  h.nblocks == want.nblocks && h.leaf == want.leaf &&
                  h.eta == want.eta && h.eps == want.eps,
                  "QD H-matrix: '" << path << "' holds N_f " << h.nf << ", "
                  << h.ncomp << " components, leaf " << h.leaf << ", eta "
                  << h.eta << ", eps " << h.eps << ", " << h.nclusters
                  << " clusters, " << h.nblocks << " blocks; this run has "
                  "N_f " << want.nf << ", " << want.ncomp << ", leaf "
                  << want.leaf << ", eta " << want.eta << ", eps " << want.eps
                  << ", " << want.nclusters << ", " << want.nblocks
                  << ".  Remove the file or pass another --qd-green-hcache.");
      h.info[sizeof(h.info) - 1] = '\0';
      info = h.info;
      // Keys, clusters, T_D and block table, checked against this run's
      std::vector<real_t> head(HeadValues());
      if (rank_ == 0)
      {
         ReadChunked(fh, 0, head.data(), head.size(), path);
      }
      for (size_t k = 0; k < head.size(); k += kChunk)
      {
         MPI_Bcast(head.data() + k,
                   static_cast<int>(std::min(kChunk, head.size() - k)),
                   MPITypeMap<real_t>::mpi_type, 0, comm_);
      }
      const real_t *fk = head.data() + kHeadVals;
      real_t ext = 1.0, dmax = 0.0;
      for (long long p = 0; p < nf_; p++)
      {
         for (int c = 0; c < kKey; c++)
         {
            const real_t a = keys_tree_[kKey * p + c];
            ext = std::max(ext, std::abs(a));
            dmax = std::max(dmax, std::abs(a - fk[kKey * p + c]));
         }
      }
      MFEM_VERIFY(dmax <= 1e-9 * ext, "QD H-matrix: '" << path << "' was "
                  "made for another mesh (fault DOF positions differ by "
                  << dmax << ").  Remove the file or pass another "
                  "--qd-green-hcache.");
      const long long *ft = reinterpret_cast<const long long *>(
                               fk + kKey * nf_);
      for (size_t c = 0; c < clusters_.size(); c++)
      {
         MFEM_VERIFY(ft[2 * c] == clusters_[c].lo &&
                     ft[2 * c + 1] == clusters_[c].hi, "QD H-matrix: '"
                     << path << "' has another cluster tree (cluster " << c
                     << ").  Remove the file or pass another "
                     "--qd-green-hcache.");
      }
      const real_t *ftd = fk + kKey * nf_ + 2 * clusters_.size();
      for (int i = 0; i < nb_; i++)
      {
         const long long p = tree_pos_[row_lo_[rank_] + i];
         for (int c = 0; c < ncomp_; c++)
         {
            td_rows[ncomp_ * i + c] = ftd[ncomp_ * p + c];
         }
      }
      const long long *fb = reinterpret_cast<const long long *>(
                               ftd + ncomp_ * nf_);
      for (size_t b = 0; b < blocks_.size(); b++)
      {
         MFEM_VERIFY(fb[4 * b] == blocks_[b].t && fb[4 * b + 1] ==
                     blocks_[b].s && (fb[4 * b + 2] >= -1) &&
                     (blocks_[b].adm || fb[4 * b + 2] == -1),
                     "QD H-matrix: '" << path << "' has another block "
                     "tree (block " << b << ").  Remove the file or pass "
                     "another --qd-green-hcache.");
         blocks_[b].k = static_cast<int>(fb[4 * b + 2]);
      }
      std::vector<long long> off;
      DataOffsets(off);
      for (size_t b = 0; b < blocks_.size(); b++)
      {
         MFEM_VERIFY(fb[4 * b + 3] == off[b], "QD H-matrix: '" << path
                     << "': block " << b << " at byte " << fb[4 * b + 3]
                     << ", expected " << off[b]);
      }
      MPI_Offset size = 0;
      MPI_File_get_size(fh, &size);
      MFEM_VERIFY(size == DataOffsets(off), "QD H-matrix: '" << path
                  << "' has " << size << " bytes, expected "
                  << DataOffsets(off) << " (incomplete write?)");
      // Each worker reads its blocks whole and hands them out
      const double t0 = MPI_Wtime();
      std::vector<int> order(blocks_.size());
      std::iota(order.begin(), order.end(), 0);
      std::stable_sort(order.begin(), order.end(), [&](int a, int b)
      {
         return StoredValues(a) > StoredValues(b);
      });
      std::vector<std::vector<int>> rounds;
      PlanRounds(order, mem_node, false, rounds);
      for (const auto &round : rounds)
      {
         RunRound(round, nullptr, &fh, false, &off, &path);
      }
      MPI_File_close(&fh);
      compress_seconds_ = MPI_Wtime() - t0;
      nrounds_ = static_cast<int>(rounds.size());
      FinishLayout();
      return true;
   }

   /// Bytes of the operator on this rank (rows, Vᵀ, w buffers).
   double LocalBytes() const { return local_bytes_; }

   /// Rank-0 summary lines (sizes, ranks, layout, time).  Collective.
   std::string Summary(const std::string &what) const
   {
      double loc[4] = {near_bytes_, u_bytes_, v_bytes_, double(wrecv_.size())};
      double sum[4], mx[4];
      MPI_Reduce(loc, sum, 4, MPI_DOUBLE, MPI_SUM, 0, comm_);
      MPI_Reduce(loc, mx, 4, MPI_DOUBLE, MPI_MAX, 0, comm_);
      double tot_max = 0.0;
      MPI_Reduce(&local_bytes_, &tot_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm_);
      std::ostringstream o;
      if (rank_ != 0) { return o.str(); }
      long long nadm = 0, nlow = 0, ndense = 0;
      std::vector<int> ks;
      double stored = 0.0;
      for (size_t b = 0; b < blocks_.size(); b++)
      {
         const Block &bk = blocks_[b];
         if (bk.adm) { nadm++; }
         if (bk.k >= 0) { nlow++; ks.push_back(bk.k); }
         else { ndense++; }
         stored += 8.0 * StoredValues(static_cast<int>(b));
      }
      std::sort(ks.begin(), ks.end());
      auto pct = [&](double q)
      {
         return ks.empty() ? 0 : ks[std::min(ks.size() - 1,
                                             size_t(q * (ks.size() - 1)))];
      };
      const double dense = 8.0 * ncomp_ * double(nf_) * 2.0 * double(nf_);
      const double gib = 1073741824.0;
      o << std::fixed << std::setprecision(3)
        << "[QD-HMAT] " << what << " at eps " << std::scientific
        << std::setprecision(1) << eps_ << std::fixed << std::setprecision(3)
        << ", leaf " << leaf_ << ", eta " << eta_ << ": " << clusters_.size()
        << " clusters, " << blocks_.size() << " blocks (" << nadm
        << " admissible, " << nlow << " low rank, " << ndense
        << " dense)\n"
        << "[QD-HMAT] stored " << stored / gib << " GiB vs dense "
        << dense / gib << " GiB (x" << std::setprecision(2) << dense / stored
        << std::setprecision(3) << " smaller); low-rank ranks median/p90/max "
        << pct(0.5) << "/" << pct(0.9) << "/" << (ks.empty() ? 0 : ks.back())
        << "\n"
        << "[QD-HMAT] on the ranks: dense rows " << sum[0] / gib
        << " GiB, U rows " << sum[1] / gib << " GiB, V " << sum[2] / gib
        << " GiB (V copies on " << nnodes_ << " node(s): +"
        << std::setprecision(1) << 100.0 * (sum[0] + sum[1] + sum[2] - stored)
           / stored << "%); largest rank " << std::setprecision(3)
        << tot_max / gib << " GiB, mean " << (sum[0] + sum[1] + sum[2]) /
           nranks_ / gib << " GiB; w received per stage max "
        << static_cast<long long>(mx[3]) << " values\n"
        << "[QD-HMAT] " << nrounds_ << " rounds, " << std::setprecision(1)
        << compress_seconds_ << " s\n";
      return o.str();
   }

private:
   struct Cluster
   {
      long long lo, hi;     ///< range in tree order
      double bmin[3], bmax[3];
      int kid[2];
   };
   struct Block
   {
      int t, s;
      bool adm;
      int k = -1;           ///< rank of a low-rank block; -1 dense
   };
   /// This rank's rows of a block: dense rows (k < 0, n values each) or
   /// rows of U·Σ (k values each); w = offset of the block's w in wbuf_.
   struct RowPart
   {
      int b, k, n;
      long long w;
      std::vector<int> rows;   ///< local row index in y
      std::vector<real_t> a;
   };
   /// Vᵀ (k × n, row-major) of a block whose w this rank computes.
   struct VPart
   {
      int b, k, n;
      long long w;
      std::vector<real_t> vt;
   };
   struct Kept
   {
      int b;
      std::vector<real_t> data;
   };
   struct Header
   {
      char magic[16];
      long long nf, ncomp, nclusters, nblocks, leaf, reserved[3];
      double eps, eta;
      char info[416];
   };
   static_assert(sizeof(Header) == 512, "QD H-matrix header is 512 B");
   static constexpr int kHeadVals = sizeof(Header) / sizeof(real_t);
   static constexpr long long kAlign = 4096;
   static constexpr size_t kChunk = size_t(1) << 27;

   MPI_Comm comm_;
   int nranks_ = 1, rank_ = 0;
   long long nf_;
   int ncomp_;
   std::vector<long long> row_lo_;
   int leaf_;
   double eta_;
   double eps_ = 0.0;
   int nb_ = 0;
   // nodes
   MPI_Comm node_comm_ = MPI_COMM_NULL;
   int nnodes_ = 1, node_rank_ = 0;
   std::vector<int> node_of_;            ///< node index of each rank
   std::vector<int> node_rank_of_;       ///< rank within its node
   std::vector<std::vector<int>> node_ranks_;
   // tree
   std::vector<long long> perm_;         ///< tree position -> g
   std::vector<long long> tree_pos_;     ///< g -> tree position
   std::vector<real_t> keys_tree_;
   std::vector<Cluster> clusters_;
   std::vector<Block> blocks_;
   std::vector<std::vector<int>> block_ranks_;   ///< ranks with rows of t
   // operator
   std::vector<RowPart> rparts_;
   std::vector<VPart> vparts_;
   std::vector<std::vector<int>> vown_;  ///< per block, V rank per node
   bool kept_ = false;
   std::vector<Kept> kept_blocks_;
   std::vector<double> v_load_;          ///< V values per rank so far
   // product buffers
   mutable std::vector<real_t> xt_, wbuf_, wsend_, wrecv_;
   std::vector<long long> wsend_idx_, wrecv_idx_;
   std::vector<int> wscount_, wsdispl_, wrcount_, wrdispl_;
   double near_bytes_ = 0.0, u_bytes_ = 0.0, v_bytes_ = 0.0,
          local_bytes_ = 0.0, compress_seconds_ = 0.0;
   int nrounds_ = 0;

   static real_t Dot(const real_t *a, const real_t *x, int n)
   {
      // Four partial sums: independent chains, so the loop streams at
      // memory speed (as QDGreenOperator::Apply).
      real_t s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
      int j = 0;
      for (; j + 4 <= n; j += 4)
      {
         s0 += a[j] * x[j];
         s1 += a[j + 1] * x[j + 1];
         s2 += a[j + 2] * x[j + 2];
         s3 += a[j + 3] * x[j + 3];
      }
      for (; j < n; j++) { s0 += a[j] * x[j]; }
      return (s0 + s1) + (s2 + s3);
   }

   long long Rows(int b) const
   {
      const Cluster &t = clusters_[blocks_[b].t];
      return ncomp_ * (t.hi - t.lo);
   }
   long long Cols(int b) const
   {
      const Cluster &s = clusters_[blocks_[b].s];
      return 2 * (s.hi - s.lo);
   }
   double SvdCost(int b) const
   {
      const double m = double(Rows(b)), n = double(Cols(b));
      return blocks_[b].adm ? m * n * std::min(m, n) : m * n;
   }
   /// Values stored for block @a b (its k must be known).
   long long StoredValues(int b) const
   {
      const long long m = Rows(b), n = Cols(b);
      return blocks_[b].k < 0 ? m * n : blocks_[b].k * (m + n);
   }
   /// Bytes a worker needs at once for block @a b.
   double WorkBytes(int b, bool svd) const
   {
      const double m = double(Rows(b)), n = double(Cols(b));
      if (!svd || !blocks_[b].adm) { return 8.0 * m * n; }
      const double mn = std::min(m, n);
      // A, Uᵀ and V of dgesdd (jobz = 'S'), its work array (~ 4 mn² + 7 mn
      // for the divide and conquer), and the copy kept for the file.
      return 8.0 * (2.0 * m * n + mn * (m + n) + 4.0 * mn * mn + 7.0 * mn);
   }

   static void CheckIO(int err, const std::string &path, const char *what)
   {
      MFEM_VERIFY(err == MPI_SUCCESS, "QD H-matrix: I/O error " << err
                  << " on " << what << " of '" << path << "'");
   }
   static void WriteChunked(MPI_File fh, MPI_Offset off, const real_t *p,
                            size_t n, const std::string &path)
   {
      for (size_t k = 0; k < n; k += kChunk)
      {
         const int c = static_cast<int>(std::min(kChunk, n - k));
         CheckIO(MPI_File_write_at(fh, off + MPI_Offset(k * sizeof(real_t)),
                                   p + k, c, MPITypeMap<real_t>::mpi_type,
                                   MPI_STATUS_IGNORE), path, "write");
      }
   }
   static void ReadChunked(MPI_File fh, MPI_Offset off, real_t *p, size_t n,
                           const std::string &path)
   {
      for (size_t k = 0; k < n; k += kChunk)
      {
         const int c = static_cast<int>(std::min(kChunk, n - k));
         CheckIO(MPI_File_read_at(fh, off + MPI_Offset(k * sizeof(real_t)),
                                  p + k, c, MPITypeMap<real_t>::mpi_type,
                                  MPI_STATUS_IGNORE), path, "read");
      }
   }

   int OwnerOf(long long g) const
   {
      return static_cast<int>(std::upper_bound(row_lo_.begin(),
                                               row_lo_.end(), g) -
                              row_lo_.begin()) - 1;
   }

   void SetupNodes()
   {
      // Test hook: SEAS_QD_HMAT_NODE_SIZE=n treats each n consecutive ranks
      // as one node, so a single machine runs the multi-node layout (V
      // copies, w kept inside a node).
      const char *fake = std::getenv("SEAS_QD_HMAT_NODE_SIZE");
      if (fake && std::atoi(fake) > 0)
      {
         MPI_Comm_split(comm_, rank_ / std::atoi(fake), rank_, &node_comm_);
      }
      else
      {
         MPI_Comm_split_type(comm_, MPI_COMM_TYPE_SHARED, rank_,
                             MPI_INFO_NULL, &node_comm_);
      }
      MPI_Comm_rank(node_comm_, &node_rank_);
      int leader = rank_;
      MPI_Bcast(&leader, 1, MPI_INT, 0, node_comm_);
      std::vector<int> leaders(nranks_);
      MPI_Allgather(&leader, 1, MPI_INT, leaders.data(), 1, MPI_INT, comm_);
      std::vector<int> ids(leaders);
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      nnodes_ = static_cast<int>(ids.size());
      node_of_.resize(nranks_);
      node_rank_of_.resize(nranks_);
      node_ranks_.assign(nnodes_, {});
      for (int r = 0; r < nranks_; r++)
      {
         node_of_[r] = static_cast<int>(std::lower_bound(ids.begin(),
                                                         ids.end(),
                                                         leaders[r]) -
                                        ids.begin());
         node_rank_of_[r] = static_cast<int>(node_ranks_[node_of_[r]].size());
         node_ranks_[node_of_[r]].push_back(r);
      }
      // MPI_Comm_split_type keys by world rank, so node_rank_of_ is the
      // rank in node_comm_; check it.
      MFEM_VERIFY(node_rank_of_[rank_] == node_rank_, "QD H-matrix: node "
                  "rank mismatch");
   }

   void BuildTree(const std::vector<real_t> &keys_all)
   {
      // Tree order starts from the DOFs sorted by key, which does not
      // depend on the partition.
      std::vector<long long> ids(nf_);
      std::iota(ids.begin(), ids.end(), 0LL);
      auto key = [&](long long g) { return &keys_all[kKey * g]; };
      std::sort(ids.begin(), ids.end(), [&](long long a, long long b)
      {
         return std::lexicographical_compare(key(a), key(a) + kKey, key(b),
                                             key(b) + kKey);
      });
      for (long long i = 1; i < nf_; i++)
      {
         MFEM_VERIFY(!std::equal(key(ids[i - 1]), key(ids[i - 1]) + kKey,
                                 key(ids[i])), "QD H-matrix: two fault DOFs "
                     "have the same position key; cannot order them");
      }
      clusters_.clear();
      clusters_.reserve(4 * (nf_ / std::max(1, leaf_)) + 4);
      Split(ids, 0, nf_, keys_all);
      perm_ = ids;
      tree_pos_.assign(nf_, -1);
      keys_tree_.resize(kKey * nf_);
      for (long long p = 0; p < nf_; p++)
      {
         tree_pos_[perm_[p]] = p;
         std::copy(key(perm_[p]), key(perm_[p]) + kKey,
                   keys_tree_.begin() + kKey * p);
      }
      xt_.assign(2 * nf_, 0.0);
   }

   /// Cluster of ids[lo, hi); returns its index.
   int Split(std::vector<long long> &ids, long long lo, long long hi,
             const std::vector<real_t> &keys_all)
   {
      const int c = static_cast<int>(clusters_.size());
      clusters_.push_back(Cluster());
      Cluster cl;
      cl.lo = lo;
      cl.hi = hi;
      cl.kid[0] = cl.kid[1] = -1;
      for (int d = 0; d < 3; d++)
      {
         cl.bmin[d] = INFINITY;
         cl.bmax[d] = -INFINITY;
      }
      for (long long i = lo; i < hi; i++)
      {
         for (int d = 0; d < 3; d++)
         {
            const double x = keys_all[kKey * ids[i] + d];
            cl.bmin[d] = std::min(cl.bmin[d], x);
            cl.bmax[d] = std::max(cl.bmax[d], x);
         }
      }
      if (hi - lo > leaf_)
      {
         // First longest side, as numpy.argmax
         int ax = 0;
         for (int d = 1; d < 3; d++)
         {
            if (cl.bmax[d] - cl.bmin[d] > cl.bmax[ax] - cl.bmin[ax]) { ax = d; }
         }
         std::stable_sort(ids.begin() + lo, ids.begin() + hi,
                          [&](long long a, long long b)
         {
            return keys_all[kKey * a + ax] < keys_all[kKey * b + ax];
         });
         const long long mid = lo + (hi - lo) / 2;
         cl.kid[0] = Split(ids, lo, mid, keys_all);
         cl.kid[1] = Split(ids, mid, hi, keys_all);
      }
      clusters_[c] = cl;
      return c;
   }

   double Diam(const Cluster &c) const
   {
      double s = 0.0;
      for (int d = 0; d < 3; d++)
      {
         s += (c.bmax[d] - c.bmin[d]) * (c.bmax[d] - c.bmin[d]);
      }
      return std::sqrt(s);
   }
   double Dist(const Cluster &a, const Cluster &b) const
   {
      double s = 0.0;
      for (int d = 0; d < 3; d++)
      {
         const double g = std::max(0.0, std::max(a.bmin[d] - b.bmax[d],
                                                 b.bmin[d] - a.bmax[d]));
         s += g * g;
      }
      return std::sqrt(s);
   }

   void BuildBlocks()
   {
      blocks_.clear();
      Walk(0, 0);
      // Ranks that hold rows of each block's row cluster
      std::vector<std::vector<int>> cl_ranks(clusters_.size());
      for (size_t c = 0; c < clusters_.size(); c++)
      {
         std::vector<int> &rs = cl_ranks[c];
         for (long long p = clusters_[c].lo; p < clusters_[c].hi; p++)
         {
            rs.push_back(OwnerOf(perm_[p]));
         }
         std::sort(rs.begin(), rs.end());
         rs.erase(std::unique(rs.begin(), rs.end()), rs.end());
      }
      block_ranks_.resize(blocks_.size());
      for (size_t b = 0; b < blocks_.size(); b++)
      {
         block_ranks_[b] = cl_ranks[blocks_[b].t];
      }
   }

   void Walk(int t, int s)
   {
      const Cluster &ct = clusters_[t], &cs = clusters_[s];
      const double d = Dist(ct, cs);
      const bool leaf = ct.kid[0] < 0 || cs.kid[0] < 0;
      // A far pair too big for one SVD goes on to its sub-pairs, which are
      // as far apart and smaller, so they are admissible in turn.
      const bool svd_fits = std::min(ncomp_ * (ct.hi - ct.lo),
                                     2 * (cs.hi - cs.lo)) <=
                            SEAS_QD_HMAT_SVD_SIDE_MAX;
      if (d > 0.0 && std::max(Diam(ct), Diam(cs)) <= eta_ * d &&
          (svd_fits || leaf))
      {
         blocks_.push_back({t, s, true, -1});
         return;
      }
      if (leaf)
      {
         blocks_.push_back({t, s, false, -1});
         return;
      }
      for (int a = 0; a < 2; a++)
      {
         for (int b = 0; b < 2; b++) { Walk(ct.kid[a], cs.kid[b]); }
      }
   }

   /// Rounds of at most one block per rank, each node spending at most
   /// @a mem_node bytes at once (a block bigger than that runs alone).
   void PlanRounds(const std::vector<int> &order, double mem_node, bool svd,
                   std::vector<std::vector<int>> &rounds) const
   {
      // Worker of each block in a round: spread over the nodes, then over
      // the ranks of a node.
      std::vector<int> left(order);
      while (!left.empty())
      {
         std::vector<int> round(nranks_, -1), next;
         std::vector<double> mem(nnodes_, 0.0);
         std::vector<int> used(nnodes_, 0);
         int nround = 0;
         for (int b : left)
         {
            const double need = WorkBytes(b, svd);
            int best = -1;
            for (int nd = 0; nd < nnodes_; nd++)
            {
               if (used[nd] >= static_cast<int>(node_ranks_[nd].size()))
               {
                  continue;
               }
               if (mem[nd] + need > mem_node && !(nround == 0)) { continue; }
               if (best < 0 || mem[nd] < mem[best]) { best = nd; }
            }
            if (best < 0) { next.push_back(b); continue; }
            round[node_ranks_[best][used[best]++]] = b;
            mem[best] += need;
            nround++;
         }
         rounds.push_back(round);
         left.swap(next);
      }
   }

   /// One round: every worker gets its block whole (gathered from the
   /// dense rows and cut by SVD, or read from the file), then sends each
   /// rank its rows and each node's V rank its Vᵀ.
   void RunRound(const std::vector<int> &round, const real_t *rows,
                 MPI_File *fh, bool keep,
                 const std::vector<long long> *off = nullptr,
                 const std::string *path = nullptr)
   {
      const int mine = round[rank_];
      std::vector<real_t> full;   // m×n dense, or U·Σ (m×k) then Vᵀ (k×n)
      int k_mine = -2;
      if (rows)
      {
         // Gather: each rank sends the worker its rows of t, columns of s
         std::vector<int> scount(nranks_, 0), sdispl(nranks_, 0),
             rcount(nranks_, 0), rdispl(nranks_, 0);
         std::vector<real_t> send;
         for (int r = 0; r < nranks_; r++)
         {
            sdispl[r] = static_cast<int>(send.size());
            const int b = round[r];
            if (b >= 0) { PackRows(b, rows, send); }
            scount[r] = static_cast<int>(send.size()) - sdispl[r];
         }
         MFEM_VERIFY(send.size() <= size_t(INT_MAX), "QD H-matrix: "
                     << send.size() << " values to send in one round, more "
                     "than MPI's int counts hold");
         if (mine >= 0)
         {
            MFEM_VERIFY(Rows(mine) * Cols(mine) <= INT_MAX, "QD H-matrix: a "
                        << Rows(mine) << " x " << Cols(mine) << " block has "
                        "more values than MPI's int counts hold");
            const Cluster &t = clusters_[blocks_[mine].t];
            const long long n = Cols(mine);
            for (long long p = t.lo; p < t.hi; p++)
            {
               rcount[OwnerOf(perm_[p])] += static_cast<int>(ncomp_ * n);
            }
            int tot = 0;
            for (int r = 0; r < nranks_; r++) { rdispl[r] = tot; tot += rcount[r]; }
            full.resize(static_cast<size_t>(tot));
         }
         MPI_Alltoallv(send.data(), scount.data(), sdispl.data(),
                       MPITypeMap<real_t>::mpi_type, full.data(),
                       rcount.data(), rdispl.data(),
                       MPITypeMap<real_t>::mpi_type, comm_);
         send.clear();
         send.shrink_to_fit();
         if (mine >= 0)
         {
            // Rows arrive by source rank; put them in tree order.
            const Cluster &t = clusters_[blocks_[mine].t];
            const long long n = Cols(mine), m = Rows(mine);
            std::vector<real_t> a(static_cast<size_t>(m * n));
            std::vector<int> pos(rdispl);
            for (long long p = t.lo; p < t.hi; p++)
            {
               const int r = OwnerOf(perm_[p]);
               std::copy(full.begin() + pos[r],
                         full.begin() + pos[r] + ncomp_ * n,
                         a.begin() + ncomp_ * (p - t.lo) * n);
               pos[r] += static_cast<int>(ncomp_ * n);
            }
            full.swap(a);
            k_mine = blocks_[mine].adm ? Truncate(mine, full) : -1;
         }
      }
      else if (mine >= 0)
      {
         full.resize(static_cast<size_t>(StoredValues(mine)));
         ReadChunked(*fh, (*off)[mine], full.data(), full.size(), *path);
         k_mine = blocks_[mine].k;
      }
      // Every rank learns the ranks found in this round
      std::vector<int> kr(nranks_);
      MPI_Allgather(&k_mine, 1, MPI_INT, kr.data(), 1, MPI_INT, comm_);
      for (int r = 0; r < nranks_; r++)
      {
         if (round[r] >= 0) { blocks_[round[r]].k = kr[r]; }
      }
      AssignV(round);
      // Results: rows to their ranks, Vᵀ to each node's V rank
      std::vector<int> scount(nranks_, 0), sdispl(nranks_, 0),
          rcount(nranks_, 0), rdispl(nranks_, 0);
      std::vector<real_t> send;
      if (mine >= 0)
      {
         const Block &bk = blocks_[mine];
         const Cluster &t = clusters_[bk.t];
         const long long m = Rows(mine), n = Cols(mine);
         const long long w = bk.k < 0 ? n : bk.k;   // values per row
         std::vector<std::vector<long long>> rows_of(nranks_);
         for (long long p = t.lo; p < t.hi; p++)
         {
            rows_of[OwnerOf(perm_[p])].push_back(p);
         }
         for (int r = 0; r < nranks_; r++)
         {
            sdispl[r] = static_cast<int>(send.size());
            for (long long p : rows_of[r])
            {
               const real_t *src = full.data() + ncomp_ * (p - t.lo) * w;
               send.insert(send.end(), src, src + ncomp_ * w);
            }
            if (bk.k >= 0 && IsVRank(mine, r))
            {
               const real_t *vt = full.data() + m * bk.k;
               send.insert(send.end(), vt, vt + bk.k * n);
            }
            scount[r] = static_cast<int>(send.size()) - sdispl[r];
         }
         MFEM_VERIFY(send.size() <= size_t(INT_MAX), "QD H-matrix: block "
                     << mine << " sends " << send.size() << " values, more "
                     "than MPI's int counts hold");
         if (keep) { kept_blocks_.push_back({mine, std::move(full)}); }
         full.clear();
         full.shrink_to_fit();
      }
      int tot = 0;
      for (int r = 0; r < nranks_; r++)
      {
         rdispl[r] = tot;
         const int b = round[r];
         if (b >= 0)
         {
            const Block &bk = blocks_[b];
            const long long n = Cols(b);
            const long long w = bk.k < 0 ? n : bk.k;
            rcount[r] = static_cast<int>(ncomp_ * MyRowsOf(bk.t) * w);
            if (bk.k >= 0 && IsVRank(b, rank_))
            {
               rcount[r] += static_cast<int>(bk.k * n);
            }
         }
         tot += rcount[r];
      }
      std::vector<real_t> recv(static_cast<size_t>(tot));
      MPI_Alltoallv(send.data(), scount.data(), sdispl.data(),
                    MPITypeMap<real_t>::mpi_type, recv.data(), rcount.data(),
                    rdispl.data(), MPITypeMap<real_t>::mpi_type, comm_);
      for (int r = 0; r < nranks_; r++)
      {
         const int b = round[r];
         if (b < 0 || rcount[r] == 0) { continue; }
         const Block &bk = blocks_[b];
         const Cluster &t = clusters_[bk.t];
         const int n = static_cast<int>(Cols(b));
         const real_t *src = recv.data() + rdispl[r];
         const long long nmy = MyRowsOf(bk.t);
         if (nmy > 0)
         {
            RowPart rp;
            rp.b = b;
            rp.k = bk.k;
            rp.n = n;
            rp.w = -1;
            const long long w = bk.k < 0 ? n : bk.k;
            for (long long p = t.lo; p < t.hi; p++)
            {
               const long long g = perm_[p];
               if (OwnerOf(g) != rank_) { continue; }
               for (int c = 0; c < ncomp_; c++)
               {
                  rp.rows.push_back(static_cast<int>(
                                       ncomp_ * (g - row_lo_[rank_]) + c));
               }
            }
            rp.a.assign(src, src + ncomp_ * nmy * w);
            src += ncomp_ * nmy * w;
            (bk.k < 0 ? near_bytes_ : u_bytes_) += 8.0 * rp.a.size();
            rparts_.push_back(std::move(rp));
         }
         if (bk.k >= 0 && IsVRank(b, rank_))
         {
            VPart vp;
            vp.b = b;
            vp.k = bk.k;
            vp.n = n;
            vp.w = -1;
            vp.vt.assign(src, src + static_cast<size_t>(bk.k) * n);
            v_bytes_ += 8.0 * vp.vt.size();
            vparts_.push_back(std::move(vp));
         }
      }
   }

   /// Rows of cluster @a c on this rank.
   long long MyRowsOf(int c) const
   {
      long long n = 0;
      for (long long p = clusters_[c].lo; p < clusters_[c].hi; p++)
      {
         if (OwnerOf(perm_[p]) == rank_) { n++; }
      }
      return n;
   }

   /// This rank's rows of block @a b (tree order), columns of s, appended.
   void PackRows(int b, const real_t *rows, std::vector<real_t> &send) const
   {
      const Cluster &t = clusters_[blocks_[b].t], &s = clusters_[blocks_[b].s];
      const size_t ncols = static_cast<size_t>(2 * nf_);
      for (long long p = t.lo; p < t.hi; p++)
      {
         const long long g = perm_[p];
         if (OwnerOf(g) != rank_) { continue; }
         for (int c = 0; c < ncomp_; c++)
         {
            const real_t *row = rows + (static_cast<size_t>(g - row_lo_[rank_])
                                        * ncomp_ + c) * ncols;
            for (long long q = s.lo; q < s.hi; q++)
            {
               const long long h = perm_[q];
               send.push_back(row[2 * h]);
               send.push_back(row[2 * h + 1]);
            }
         }
      }
   }

   /// Cut block @a b (m×n row-major in @a a) by SVD.  Returns k (and
   /// leaves U·Σ, m×k, then Vᵀ, k×n, in @a a) or -1 (a unchanged).
   int Truncate(int b, std::vector<real_t> &a) const
   {
#ifdef MFEM_USE_LAPACK
      // dgesdd (jobz = 'S') needs at most 4·mn² + 7·mn + 64·(m + n)
      // workspace entries (above every query of reference LAPACK and
      // OpenBLAS up to mn = 23,000).  Above 2³¹ its own query overflows and
      // returns a too small size with info 0, so check here, in 64 bits.
      const long long mn64 = std::min(Rows(b), Cols(b));
      const long long need = 4 * mn64 * mn64 + 7 * mn64 +
                             64 * (Rows(b) + Cols(b));
      MFEM_VERIFY(need < INT_MAX, "QD H-matrix: a " << Rows(b) << " x "
                  << Cols(b) << " block needs about " << need << " dgesdd "
                  "workspace entries, more than 32-bit LAPACK can address; "
                  "lower SEAS_QD_HMAT_SVD_SIDE_MAX");
      int m = static_cast<int>(Rows(b)), n = static_cast<int>(Cols(b));
      // Row-major m×n is column-major Aᵀ (n×m): Aᵀ = V Σ Uᵀ, so dgesdd's
      // "U" is V (n×mn) and its "Vᵀ" is Uᵀ (mn×m).
      const int mn = std::min(m, n);
      std::vector<real_t> at(a), sv(mn), vv(static_cast<size_t>(n) * mn),
          ut(static_cast<size_t>(mn) * m);
      std::vector<int> iwork(8 * static_cast<size_t>(mn));
      char jobz = 'S';
      int lda = n, ldu = n, ldvt = mn, lwork = -1, info = 0;
      real_t wq = 0.0;
      dgesdd_(&jobz, &n, &m, at.data(), &lda, sv.data(), vv.data(), &ldu,
              ut.data(), &ldvt, &wq, &lwork, iwork.data(), &info);
      MFEM_VERIFY(info == 0 && wq < double(INT_MAX), "QD H-matrix: dgesdd "
                  "workspace query failed (info " << info << ", lwork " << wq
                  << ") for a " << m << " x " << n << " block");
      lwork = static_cast<int>(wq) + 1;
      std::vector<real_t> work(static_cast<size_t>(lwork));
      dgesdd_(&jobz, &n, &m, at.data(), &lda, sv.data(), vv.data(), &ldu,
              ut.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);
      MFEM_VERIFY(info == 0, "QD H-matrix: dgesdd failed (info " << info
                  << ") for a " << m << " x " << n << " block");
      work.clear();
      work.shrink_to_fit();
      at.clear();
      at.shrink_to_fit();
      // Smallest k with sqrt(sum_{i >= k} s_i²) <= eps·||A||_F
      std::vector<double> tail(mn + 1, 0.0);
      for (int i = mn - 1; i >= 0; i--) { tail[i] = tail[i + 1] + sv[i] * sv[i]; }
      const double fro = std::sqrt(tail[0]);
      int k = 0;
      while (k < mn && std::sqrt(tail[k]) > eps_ * fro) { k++; }
      if (static_cast<long long>(k) * (m + n) >= static_cast<long long>(m) * n)
      {
         return -1;
      }
      a.assign(static_cast<size_t>(k) * (m + n), 0.0);
      for (int i = 0; i < m; i++)
      {
         for (int l = 0; l < k; l++)
         {
            a[static_cast<size_t>(i) * k + l] = ut[static_cast<size_t>(i) * mn
                                                   + l] * sv[l];
         }
      }
      std::copy(vv.begin(), vv.begin() + static_cast<size_t>(k) * n,
                a.begin() + static_cast<size_t>(m) * k);
      return k;
#else
      MFEM_ABORT("QD H-matrix: needs LAPACK");
      return -1;
#endif
   }

   /// The rank on each node that holds Vᵀ of the low-rank blocks of
   /// @a round: among the node's ranks with rows of the block, the one
   /// with the fewest V values so far.  Same on every rank.
   void AssignV(const std::vector<int> &round)
   {
      if (vown_.empty())
      {
         vown_.assign(blocks_.size(), {});
         v_load_.assign(nranks_, 0.0);
      }
      for (int r = 0; r < nranks_; r++)
      {
         const int b = round[r];
         if (b < 0 || blocks_[b].k < 0) { continue; }
         std::vector<int> &vo = vown_[b];
         vo.assign(nnodes_, -1);
         for (int q : block_ranks_[b])
         {
            const int nd = node_of_[q];
            if (vo[nd] < 0 || v_load_[q] < v_load_[vo[nd]]) { vo[nd] = q; }
         }
         for (int nd = 0; nd < nnodes_; nd++)
         {
            if (vo[nd] >= 0)
            {
               v_load_[vo[nd]] += double(blocks_[b].k) * Cols(b);
            }
         }
      }
   }

   bool IsVRank(int b, int r) const
   {
      return !vown_.empty() && !vown_[b].empty() &&
             vown_[b][node_of_[r]] == r;
   }

   /// w offsets and the node exchange of w, after all rounds.
   void FinishLayout()
   {
      // wbuf_: one slot per low-rank block this rank needs (rows or V)
      std::vector<long long> wpos(blocks_.size(), -1);
      long long nw = 0;
      auto slot = [&](int b)
      {
         if (wpos[b] < 0) { wpos[b] = nw; nw += blocks_[b].k; }
         return wpos[b];
      };
      std::sort(rparts_.begin(), rparts_.end(),
                [](const RowPart &x, const RowPart &y) { return x.b < y.b; });
      std::sort(vparts_.begin(), vparts_.end(),
                [](const VPart &x, const VPart &y) { return x.b < y.b; });
      for (VPart &v : vparts_) { v.w = slot(v.b); }
      for (RowPart &r : rparts_) { if (r.k >= 0) { r.w = slot(r.b); } }
      wbuf_.assign(static_cast<size_t>(nw), 0.0);
      // Send: for each node rank d (≠ me), the blocks whose V I hold and
      // d holds rows of, by block; receive likewise from each V rank.
      const std::vector<int> &mates = node_ranks_[node_of_[rank_]];
      const int nn = static_cast<int>(mates.size());
      wscount_.assign(nn, 0);
      wsdispl_.assign(nn, 0);
      wrcount_.assign(nn, 0);
      wrdispl_.assign(nn, 0);
      wsend_idx_.clear();
      wrecv_idx_.clear();
      for (int d = 0; d < nn; d++)
      {
         wsdispl_[d] = static_cast<int>(wsend_idx_.size());
         const int rd = mates[d];
         if (rd != rank_)
         {
            for (const VPart &v : vparts_)
            {
               const auto &br = block_ranks_[v.b];
               if (std::binary_search(br.begin(), br.end(), rd))
               {
                  for (int l = 0; l < v.k; l++) { wsend_idx_.push_back(v.w + l); }
               }
            }
         }
         wscount_[d] = static_cast<int>(wsend_idx_.size()) - wsdispl_[d];
      }
      for (int d = 0; d < nn; d++)
      {
         wrdispl_[d] = static_cast<int>(wrecv_idx_.size());
         const int rd = mates[d];
         if (rd != rank_)
         {
            for (const RowPart &r : rparts_)
            {
               if (r.k >= 0 && vown_[r.b][node_of_[rank_]] == rd)
               {
                  for (int l = 0; l < r.k; l++) { wrecv_idx_.push_back(r.w + l); }
               }
            }
         }
         wrcount_[d] = static_cast<int>(wrecv_idx_.size()) - wrdispl_[d];
      }
      wsend_.assign(wsend_idx_.size(), 0.0);
      wrecv_.assign(wrecv_idx_.size(), 0.0);
      // Every low-rank block with rows here has its V on this node
      for (const RowPart &r : rparts_)
      {
         MFEM_VERIFY(r.k < 0 || vown_[r.b][node_of_[rank_]] >= 0,
                     "QD H-matrix: block " << r.b << " has no V on node "
                     << node_of_[rank_]);
      }
      local_bytes_ = near_bytes_ + u_bytes_ + v_bytes_ +
                     8.0 * (wbuf_.size() + xt_.size());
   }

   void GatherTD(const real_t *td_rows, std::vector<real_t> &td_tree) const
   {
      std::vector<int> cnt(nranks_), dsp(nranks_);
      for (int r = 0; r < nranks_; r++)
      {
         cnt[r] = static_cast<int>(ncomp_ * (row_lo_[r + 1] - row_lo_[r]));
         dsp[r] = static_cast<int>(ncomp_ * row_lo_[r]);
      }
      std::vector<real_t> all(rank_ == 0 ? ncomp_ * nf_ : 0);
      MPI_Gatherv(td_rows, cnt[rank_], MPITypeMap<real_t>::mpi_type,
                  all.data(), cnt.data(), dsp.data(),
                  MPITypeMap<real_t>::mpi_type, 0, comm_);
      td_tree.assign(all.size(), 0.0);
      for (long long p = 0; rank_ == 0 && p < nf_; p++)
      {
         for (int c = 0; c < ncomp_; c++)
         {
            td_tree[ncomp_ * p + c] = all[ncomp_ * perm_[p] + c];
         }
      }
   }

   Header MakeHeader(const std::string &info) const
   {
      Header h;
      std::memset(&h, 0, sizeof(Header));
      std::memcpy(h.magic, "SEAS-QD-HMAT-1", 14);
      h.nf = nf_;
      h.ncomp = ncomp_;
      h.nclusters = static_cast<long long>(clusters_.size());
      h.nblocks = static_cast<long long>(blocks_.size());
      h.leaf = leaf_;
      h.eps = eps_;
      h.eta = eta_;
      std::strncpy(h.info, info.c_str(), sizeof(h.info) - 1);
      return h;
   }

   /// Values (8 B each) of header, keys, clusters, T_D and block table.
   size_t HeadValues() const
   {
      return kHeadVals + kKey * nf_ + 2 * clusters_.size() + ncomp_ * nf_ +
             4 * blocks_.size();
   }

   /// Byte offset of each block's data (from 4 KiB boundaries after the
   /// head); returns the file size.  Needs every block's k.
   MPI_Offset DataOffsets(std::vector<long long> &off) const
   {
      auto up = [](long long x) { return (x + kAlign - 1) / kAlign * kAlign; };
      long long pos = up(static_cast<long long>(HeadValues()) * 8);
      off.resize(blocks_.size());
      for (size_t b = 0; b < blocks_.size(); b++)
      {
         off[b] = pos;
         pos = up(pos + 8 * StoredValues(static_cast<int>(b)));
      }
      return pos;
   }

   void HeadBytes(const std::string &info, const std::vector<real_t> &td_tree,
                  std::vector<real_t> &head) const
   {
      head.assign(HeadValues(), 0.0);
      const Header h = MakeHeader(info);
      std::memcpy(head.data(), &h, sizeof(Header));
      real_t *p = head.data() + kHeadVals;
      std::copy(keys_tree_.begin(), keys_tree_.end(), p);
      p += kKey * nf_;
      long long *cl = reinterpret_cast<long long *>(p);
      for (size_t c = 0; c < clusters_.size(); c++)
      {
         cl[2 * c] = clusters_[c].lo;
         cl[2 * c + 1] = clusters_[c].hi;
      }
      p += 2 * clusters_.size();
      std::copy(td_tree.begin(), td_tree.end(), p);
      p += ncomp_ * nf_;
      std::vector<long long> off;
      DataOffsets(off);
      long long *bt = reinterpret_cast<long long *>(p);
      for (size_t b = 0; b < blocks_.size(); b++)
      {
         bt[4 * b] = blocks_[b].t;
         bt[4 * b + 1] = blocks_[b].s;
         bt[4 * b + 2] = blocks_[b].k;
         bt[4 * b + 3] = off[b];
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_SEAS_QD_HMATRIX_OPERATOR_HPP
