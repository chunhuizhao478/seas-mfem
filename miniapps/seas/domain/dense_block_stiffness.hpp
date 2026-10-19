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

#ifndef MFEM_SEAS_DENSE_BLOCK_STIFFNESS_HPP
#define MFEM_SEAS_DENSE_BLOCK_STIFFNESS_HPP

// y = K u, or y = K u + E v in one sweep, from dense blocks: profiling plan
// item 8a (lts-kernels-design-rev1 §4, option a2; Chunhui Zhao, 2026-10-18).
//
// K (elastic) and E = K_eta (Kelvin-Voigt) are SIPG operators on a DG space,
// so a row of element e couples only to e itself and to the elements that
// share a face with it.  Init() reads the assembled hypre matrices into
//   - one diagonal block per element, upper triangle only,
//   - one block per interior face, the rows of its lower-numbered element;
//     the rows of the other element are the transpose (K and E are symmetric:
//     DGElasticityIPCombinedIntegrator symmetrises each face matrix), so one
//     read of the block serves both products,
//   - one block per shared face: this rank's rows against the DOFs of the
//     neighbour element (the halo),
// stored row by row, the row of K next to the row of E, so one pass over the
// blocks gives both products.  The entries are the assembled ones; only the
// order of the sums changes.  At p = 1 the blocks of both operators take
// about 6 kB per tet against 17 kB for the two CSR matrices.
//
// Every entry of each matrix must land in a block or, for the lower triangle
// and the transposed face rows, match the stored mirror entry to 1e-12 of the
// largest entry (Asymmetry(m)).  Otherwise Init returns false and the caller
// keeps the CSR products.  The halo is the offd column map of K (or of E when E's contains
// K's), exchanged with K's hypre communication package, u and v together.
// After Init the caller may free both HypreParMatrix objects.

#include "mfem.hpp"
#include "../common/stage_timers.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace mfem
{
namespace seas
{

class DenseBlockStiffness
{
public:
   /// Time the halo wait of Mult as StageTimers::DynHaloWait (null: off).
   void SetStageTimers(StageTimers *timers) { timers_ = timers; }

   /// Read the blocks of K (and of E, when not null).  False (nothing
   /// usable) when the space is not DG with one block size, a matrix has an
   /// entry outside the element and face blocks, or the halo cannot be
   /// shared; the caller reduces the flag over all ranks.  Collective
   /// through the face-neighbour exchange of the mesh and space.
   bool Init(ParFiniteElementSpace &fes, HypreParMatrix &K, HypreParMatrix *E)
   {
      ok_ = false;
      nm_ = E ? 2 : 1;
      ParMesh &pm = *fes.GetParMesh();
      pm.ExchangeFaceNbrData();
      fes.ExchangeFaceNbrData();
      n_ = fes.GetVSize();
      hypre_ParCSRMatrix *mat[2] = {K, E ? static_cast<hypre_ParCSRMatrix *>(*E) : nullptr};
      // the halo pattern is hypre's matvec communication package; making it
      // is collective, so it is made here, before any rank can return
      for (int m = 0; m < nm_; m++)
      {
         if (!hypre_ParCSRMatrixCommPkg(mat[m])) { hypre_MatvecCommPkgCreate(mat[m]); }
      }
      bool good = (fes.GetTrueVSize() == n_ && fes.GetMyDofOffset() == fes.GetMyTDofOffset());
      for (int m = 0; m < nm_; m++)
      {
         hypre_CSRMatrix *d = hypre_ParCSRMatrixDiag(mat[m]);
         good = good && hypre_CSRMatrixNumRows(d) == n_ && hypre_CSRMatrixNumCols(d) == n_;
      }
      // element DOFs: every local DOF in exactly one element, one block size
      ne_ = pm.GetNE();
      Array<int> vd;
      nb_ = 0;
      if (ne_ > 0) { fes.GetElementVDofs(0, vd); nb_ = vd.Size(); }
      std::vector<int> dof_elem(n_, -1), dof_pos(n_, -1);
      idx_.assign(static_cast<std::size_t>(ne_) * nb_, 0);
      for (int e = 0; good && e < ne_; e++)
      {
         fes.GetElementVDofs(e, vd);
         if (vd.Size() != nb_) { good = false; break; }
         for (int a = 0; a < nb_; a++)
         {
            const int i = vd[a];
            if (i < 0 || i >= n_ || dof_elem[i] >= 0) { good = false; break; }
            dof_elem[i] = e;
            dof_pos[i] = a;
            idx_[static_cast<std::size_t>(e) * nb_ + a] = i;
         }
      }
      for (int i = 0; good && i < n_; i++) { if (dof_elem[i] < 0) { good = false; } }
      if (!good || nb_ > kMaxNB) { return false; }
      tri_ = nb_ * (nb_ + 1) / 2;

      // interior faces, grouped by their lower-numbered element (the
      // element-order sweep of Mult); adjacency for the lookups below
      std::vector<std::pair<int, int>> fl;
      for (int f = 0; f < pm.GetNumFaces(); f++)
      {
         int e1, e2;
         pm.GetFaceElements(f, &e1, &e2);
         if (e1 >= 0 && e2 >= 0 && e1 != e2) { fl.emplace_back(std::min(e1, e2), std::max(e1, e2)); }
      }
      std::sort(fl.begin(), fl.end());
      fl.erase(std::unique(fl.begin(), fl.end()), fl.end());
      nf_ = static_cast<int>(fl.size());
      own_.assign(ne_ + 1, 0);
      fhi_.resize(nf_);
      for (int f = 0; f < nf_; f++) { own_[fl[f].first + 1]++; fhi_[f] = fl[f].second; }
      for (int e = 0; e < ne_; e++) { own_[e + 1] += own_[e]; }
      // per element: (neighbour, face) for both ends of each face
      std::vector<int> aptr(ne_ + 1, 0);
      for (auto &p : fl) { aptr[p.first + 1]++; aptr[p.second + 1]++; }
      for (int e = 0; e < ne_; e++) { aptr[e + 1] += aptr[e]; }
      std::vector<int> anbr(aptr[ne_]), aface(aptr[ne_]), fill(aptr.begin(), aptr.end() - 1);
      for (int f = 0; f < nf_; f++)
      {
         anbr[fill[fl[f].first]] = fl[f].second;
         aface[fill[fl[f].first]++] = f;
         anbr[fill[fl[f].second]] = fl[f].first;
         aface[fill[fl[f].second]++] = f;
      }
      auto find_face = [&](int e, int g) {
         for (int q = aptr[e]; q < aptr[e + 1]; q++) { if (anbr[q] == g) { return aface[q]; } }
         return -1;
      };

      // halo host: the matrix whose offd column map contains the other's
      hypre_ParCSRMatrix *host = mat[0];
      if (nm_ == 2 && !Contains_(mat[0], mat[1]))
      {
         if (!Contains_(mat[1], mat[0])) { return false; }
         host = mat[1];
      }
      const HYPRE_Int nh = hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(host));
      const HYPRE_BigInt *hmap = hypre_ParCSRMatrixColMapOffd(host);
      // shared faces: (local element, neighbour element); each halo column
      // belongs to one neighbour element at one position
      const int nsf = pm.GetNSharedFaces();
      std::vector<std::pair<int, int>> sl;
      sl.reserve(nsf);
      for (int sf = 0; sf < nsf; sf++)
      {
         int e1, e2;
         pm.GetFaceElements(pm.GetSharedFace(sf), &e1, &e2);
         if (e1 < 0 || e2 >= 0) { return false; }
         sl.emplace_back(e1, -1 - e2);
      }
      std::sort(sl.begin(), sl.end());
      sl.erase(std::unique(sl.begin(), sl.end()), sl.end());
      ns_ = static_cast<int>(sl.size());
      se_.resize(ns_);
      hcol_.assign(static_cast<std::size_t>(ns_) * nb_, static_cast<int>(nh));   // nh = a zero slot
      std::vector<int> h_nbr(nh, -1), h_pos(nh, -1);
      const HYPRE_BigInt *gmap = fes.GetFaceNbrGlobalDofMap();
      for (int s = 0; s < ns_; s++)
      {
         se_[s] = sl[s].first;
         fes.GetFaceNbrElementVDofs(sl[s].second, vd);
         if (vd.Size() != nb_) { return false; }
         for (int k = 0; k < nb_; k++)
         {
            if (vd[k] < 0) { return false; }
            const HYPRE_BigInt g = gmap[vd[k]];
            const HYPRE_BigInt *it = std::lower_bound(hmap, hmap + nh, g);
            if (it == hmap + nh || *it != g) { continue; }   // no entry of this column: zero slot
            const int hp = static_cast<int>(it - hmap);
            if (h_nbr[hp] >= 0 && (h_nbr[hp] != sl[s].second || h_pos[hp] != k)) { return false; }
            h_nbr[hp] = sl[s].second;
            h_pos[hp] = k;
            hcol_[static_cast<std::size_t>(s) * nb_ + k] = hp;
         }
      }
      // per local element: its shared blocks (sl is sorted by element)
      sptr_.assign(ne_ + 1, 0);
      for (int s = 0; s < ns_; s++) { sptr_[se_[s] + 1]++; }
      for (int e = 0; e < ne_; e++) { sptr_[e + 1] += sptr_[e]; }
      auto find_shared = [&](int e, int g) {
         for (int s = sptr_[e]; s < sptr_[e + 1]; s++) { if (sl[s].second == g) { return s; } }
         return -1;
      };

      // the blocks, values interleaved by matrix
      const std::size_t nb2 = static_cast<std::size_t>(nb_) * nb_;
      D_.assign(static_cast<std::size_t>(ne_) * tri_ * nm_, 0.0);
      B_.assign(static_cast<std::size_t>(nf_) * nb2 * nm_, 0.0);
      S_.assign(static_cast<std::size_t>(ns_) * nb2 * nm_, 0.0);
      // layout: row by row, the row of K then the row of E.  Diagonal
      // block: row a holds columns a..nb-1 (packed upper triangle); face
      // and shared blocks: full rows.
      auto dslot = [&](int e, int a, int b, int m) {   // a <= b
         return static_cast<std::size_t>(e) * tri_ * nm_
                + static_cast<std::size_t>(nm_) * (a * nb_ - a * (a - 1) / 2) + m * (nb_ - a) + (b - a);
      };
      auto fslot = [&](int f, int a, int b, int m) {
         return (static_cast<std::size_t>(f) * nb_ * nm_ + static_cast<std::size_t>(a) * nm_ + m) * nb_ + b;
      };
      asym_[0] = asym_[1] = amax_[0] = amax_[1] = 0.0;
      long long bad = 0;
      for (int m = 0; m < nm_; m++)
      {
         // pass 1 stores the upper triangles and the rows of the lower
         // element of each face; pass 2 compares the mirror entries
         for (int pass = 0; pass < 2; pass++)
         {
            hypre_CSRMatrix *d = hypre_ParCSRMatrixDiag(mat[m]);
            const HYPRE_Int *I = hypre_CSRMatrixI(d), *J = hypre_CSRMatrixJ(d);
            const HYPRE_Real *A = hypre_CSRMatrixData(d);
            for (int i = 0; i < n_; i++)
            {
               const int e = dof_elem[i], a = dof_pos[i];
               for (HYPRE_Int p = I[i]; p < I[i + 1]; p++)
               {
                  const int j = J[p], g = dof_elem[j], b = dof_pos[j];
                  const real_t val = A[p];
                  if (g == e)
                  {
                     const bool upper = (a <= b);
                     real_t &slot = D_[dslot(e, std::min(a, b), std::max(a, b), m)];
                     if (pass == 0 && upper) { slot = val; amax_[m] = std::max(amax_[m], std::abs(val)); }
                     if (pass == 1 && !upper) { asym_[m] = std::max(asym_[m], std::abs(val - slot)); }
                     continue;
                  }
                  const int f = find_face(e, g);
                  if (f < 0) { if (pass == 0) { bad++; } continue; }
                  const bool low = (e < g);
                  const std::size_t q = low ? fslot(f, a, b, m) : fslot(f, b, a, m);
                  if (pass == 0 && low) { B_[q] = val; amax_[m] = std::max(amax_[m], std::abs(val)); }
                  if (pass == 1 && !low) { asym_[m] = std::max(asym_[m], std::abs(val - B_[q])); }
               }
            }
            if (pass == 1) { continue; }
            hypre_CSRMatrix *o = hypre_ParCSRMatrixOffd(mat[m]);
            const HYPRE_Int *Io = hypre_CSRMatrixI(o), *Jo = hypre_CSRMatrixJ(o);
            const HYPRE_Real *Ao = hypre_CSRMatrixData(o);
            const HYPRE_Int nc = hypre_CSRMatrixNumCols(o);
            const HYPRE_BigInt *cmap = hypre_ParCSRMatrixColMapOffd(mat[m]);
            std::vector<int> to_host(nc);
            for (HYPRE_Int c = 0; c < nc; c++)
            {
               const HYPRE_BigInt *it = std::lower_bound(hmap, hmap + nh, cmap[c]);
               to_host[c] = (it != hmap + nh && *it == cmap[c]) ? static_cast<int>(it - hmap) : nh;
            }
            for (int i = 0; i < n_; i++)
            {
               const int e = dof_elem[i], a = dof_pos[i];
               for (HYPRE_Int p = Io[i]; p < Io[i + 1]; p++)
               {
                  const int hp = to_host[Jo[p]];
                  const int s = (hp < nh && h_nbr[hp] >= 0) ? find_shared(e, h_nbr[hp]) : -1;
                  if (s < 0) { bad++; continue; }
                  S_[fslot(s, a, h_pos[hp], m)] = Ao[p];
                  amax_[m] = std::max(amax_[m], std::abs(static_cast<real_t>(Ao[p])));
               }
            }
         }
      }
      if (bad > 0) { return false; }
      // the transposed rows are read from the stored ones, so both matrices
      // must be symmetric to round-off (measured: below 1e-16 of the largest
      // entry for the SIPG K and K_eta)
      for (int m = 0; m < nm_; m++) { if (asym_[m] > 1e-12 * amax_[m]) { return false; } }

      // halo pattern of the host
      hypre_ParCSRCommPkg *pkg = hypre_ParCSRMatrixCommPkg(host);
      comm_ = hypre_ParCSRCommPkgComm(pkg);
      const HYPRE_Int nsend = hypre_ParCSRCommPkgNumSends(pkg);
      const HYPRE_Int nrecv = hypre_ParCSRCommPkgNumRecvs(pkg);
      send_procs_.assign(hypre_ParCSRCommPkgSendProcs(pkg), hypre_ParCSRCommPkgSendProcs(pkg) + nsend);
      send_starts_.assign(hypre_ParCSRCommPkgSendMapStarts(pkg),
                          hypre_ParCSRCommPkgSendMapStarts(pkg) + nsend + 1);
      send_map_.assign(hypre_ParCSRCommPkgSendMapElmts(pkg),
                       hypre_ParCSRCommPkgSendMapElmts(pkg) + send_starts_[nsend]);
      recv_procs_.assign(hypre_ParCSRCommPkgRecvProcs(pkg), hypre_ParCSRCommPkgRecvProcs(pkg) + nrecv);
      recv_starts_.assign(hypre_ParCSRCommPkgRecvVecStarts(pkg),
                          hypre_ParCSRCommPkgRecvVecStarts(pkg) + nrecv + 1);
      MFEM_VERIFY(recv_starts_[nrecv] == nh, "DenseBlockStiffness: halo size " << recv_starts_[nrecv]
                  << " != offd columns " << nh);
      send_buf_.resize(static_cast<std::size_t>(nm_) * send_starts_[nsend]);
      recv_buf_.assign(static_cast<std::size_t>(nm_) * (nh + 1), 0.0);   // + the zero slot
      reqs_.resize(nsend + nrecv);
      ok_ = true;
      return true;
   }

   bool Ok() const { return ok_; }
   /// 1 after Init with K only, 2 with K and E.
   int NumMatrices() const { return nm_; }
   /// Block size (3 x nodes per tet), and the counts of this rank.
   int BlockSize() const { return nb_; }
   int NumElements() const { return ne_; }
   int NumFaces() const { return nf_; }
   int NumSharedBlocks() const { return ns_; }
   /// Bytes of the stored blocks and index arrays on this rank.
   long long Bytes() const
   {
      return static_cast<long long>(sizeof(real_t)) * (D_.size() + B_.size() + S_.size())
             + static_cast<long long>(sizeof(int)) * (idx_.size() + hcol_.size() + own_.size()
                                                     + fhi_.size() + sptr_.size() + se_.size());
   }
   /// Largest |A_ij - A_ji| over the mirror pairs of this rank, and the
   /// largest |A_ij|, of matrix m (0 = K, 1 = E).
   real_t Asymmetry(int m) const { return asym_[m]; }
   real_t MaxEntry(int m) const { return amax_[m]; }

   /// y = K u (Init with E == nullptr).
   void Mult(const Vector &u, Vector &y) const
   {
      MFEM_ASSERT(ok_ && nm_ == 1, "DenseBlockStiffness::Mult(u, y): Init with K only");
      Dispatch_<1>(u.HostRead(), nullptr, y.HostWrite());
   }
   /// y = K u + E v (Init with E).
   void Mult(const Vector &u, const Vector &v, Vector &y) const
   {
      MFEM_ASSERT(ok_ && nm_ == 2, "DenseBlockStiffness::Mult(u, v, y): Init with K and E");
      Dispatch_<2>(u.HostRead(), v.HostRead(), y.HostWrite());
   }

private:
   static constexpr int kMaxNB = 168;   // 3 x 56 nodes: p = 5 tets

   /// True when every offd column of sub is an offd column of host.
   static bool Contains_(hypre_ParCSRMatrix *host, hypre_ParCSRMatrix *sub)
   {
      const HYPRE_Int ns = hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(sub));
      const HYPRE_Int nh = hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(host));
      const HYPRE_BigInt *cs = hypre_ParCSRMatrixColMapOffd(sub);
      const HYPRE_BigInt *ch = hypre_ParCSRMatrixColMapOffd(host);
      for (HYPRE_Int j = 0; j < ns; j++)
      {
         const HYPRE_BigInt *it = std::lower_bound(ch, ch + nh, cs[j]);
         if (it == ch + nh || *it != cs[j]) { return false; }
      }
      return true;
   }

   // Two-wide vectors (GCC/Clang vector extension; one NEON register on
   // Vista's Grace cores, half an AVX register on x86).  The row loops below
   // keep two independent partial sums: one chain of multiply-adds is
   // latency bound.
   typedef real_t V2 __attribute__((vector_size(2 * sizeof(real_t))));
   static V2 Load2_(const real_t *p) { V2 v; __builtin_memcpy(&v, p, sizeof(V2)); return v; }
   static void Store2_(real_t *p, V2 v) { __builtin_memcpy(p, &v, sizeof(V2)); }

   /// Returns sum_b r[b] x[b] and adds r[b] xa to y[b] (one block row: the
   /// product with the neighbour's DOFs and the transposed product).
   template <int NBC>
   static real_t RowDotAxpy_(int nb, const real_t *r, const real_t *x, real_t xa, real_t *y)
   {
      const int n = NBC > 0 ? NBC : nb;
      const V2 xav = {xa, xa};
      V2 s0 = {0.0, 0.0}, s1 = {0.0, 0.0};
      int b = 0;
      for (; b + 4 <= n; b += 4)
      {
         const V2 r0 = Load2_(r + b), r1 = Load2_(r + b + 2);
         s0 += r0 * Load2_(x + b);
         s1 += r1 * Load2_(x + b + 2);
         Store2_(y + b, Load2_(y + b) + r0 * xav);
         Store2_(y + b + 2, Load2_(y + b + 2) + r1 * xav);
      }
      for (; b + 2 <= n; b += 2)
      {
         const V2 r0 = Load2_(r + b);
         s0 += r0 * Load2_(x + b);
         Store2_(y + b, Load2_(y + b) + r0 * xav);
      }
      real_t s = 0.0;
      for (; b < n; b++)
      {
         s += r[b] * x[b];
         y[b] += r[b] * xa;
      }
      const V2 t = s0 + s1;
      return s + (t[0] + t[1]);
   }

   /// Returns sum_b r[b] x[b].
   template <int NBC>
   static real_t RowDot_(int nb, const real_t *r, const real_t *x)
   {
      const int n = NBC > 0 ? NBC : nb;
      V2 s0 = {0.0, 0.0}, s1 = {0.0, 0.0};
      int b = 0;
      for (; b + 4 <= n; b += 4)
      {
         s0 += Load2_(r + b) * Load2_(x + b);
         s1 += Load2_(r + b + 2) * Load2_(x + b + 2);
      }
      for (; b + 2 <= n; b += 2) { s0 += Load2_(r + b) * Load2_(x + b); }
      real_t s = 0.0;
      for (; b < n; b++) { s += r[b] * x[b]; }
      const V2 t = s0 + s1;
      return s + (t[0] + t[1]);
   }

   /// Fixed-size loops for the block sizes of p = 1 and p = 2 (12, 30).
   template <int NM>
   void Dispatch_(const real_t *u, const real_t *v, real_t *y) const
   {
      if (nb_ == 12) { Apply_<NM, 12>(u, v, y); }
      else if (nb_ == 30) { Apply_<NM, 30>(u, v, y); }
      else { Apply_<NM, 0>(u, v, y); }
   }

   template <int NM, int NBC>
   void Apply_(const real_t *u, const real_t *v, real_t *y) const
   {
      const int nb = NBC > 0 ? NBC : nb_;
      constexpr int L = NBC > 0 ? NBC : kMaxNB;
      const real_t *x[2] = {u, v};
      // halo: post the exchange of u (and v), as FusedStiffnessPair
      const int nsend = static_cast<int>(send_procs_.size());
      const int nrecv = static_cast<int>(recv_procs_.size());
      const int tag = 7322;
      int q = 0;
      for (int r = 0; r < nrecv; r++)
      {
         MPI_Irecv(recv_buf_.data() + NM * recv_starts_[r], NM * (recv_starts_[r + 1] - recv_starts_[r]),
                   MPITypeMap<real_t>::mpi_type, recv_procs_[r], tag, comm_, &reqs_[q++]);
      }
      for (int s = 0; s < nsend; s++)
      {
         for (int i = send_starts_[s]; i < send_starts_[s + 1]; i++)
         {
            for (int m = 0; m < NM; m++) { send_buf_[NM * i + m] = x[m][send_map_[i]]; }
         }
         MPI_Isend(send_buf_.data() + NM * send_starts_[s], NM * (send_starts_[s + 1] - send_starts_[s]),
                   MPITypeMap<real_t>::mpi_type, send_procs_[s], tag, comm_, &reqs_[q++]);
      }
      std::fill(y, y + n_, 0.0);
      // element-order sweep: when element e is reached, y_e holds the
      // transposed face terms of the lower elements; e adds its diagonal
      // block and its own faces (B x_h into y_e, B^T x_e into y_h)
      real_t xe[NM][L], xh[NM][L], acc[L], yh[L];
      const int *ix = idx_.data();
      const real_t *Dp = D_.data(), *Bp = B_.data();
      const std::size_t nbb = static_cast<std::size_t>(nb) * nb * NM;
      for (int e = 0; e < ne_; e++)
      {
         const int *ie = ix + static_cast<std::size_t>(e) * nb;
         for (int a = 0; a < nb; a++)
         {
            for (int m = 0; m < NM; m++) { xe[m][a] = x[m][ie[a]]; }
            acc[a] = y[ie[a]];
         }
         const real_t *d = Dp + static_cast<std::size_t>(e) * tri_ * NM;
         for (int a = 0; a < nb; a++)
         {
            const int len = nb - a;
            V2 sv = {0.0, 0.0};
            real_t s0 = 0.0;
            for (int m = 0; m < NM; m++)
            {
               const real_t *dr = d + m * len;
               const real_t xa = xe[m][a], *xm = xe[m] + a;
               const V2 xav = {xa, xa};
               real_t *am = acc + a;
               s0 += dr[0] * xa;   // the diagonal entry
               int b = 1;
               for (; b + 1 < len; b += 2)
               {
                  const V2 dv = Load2_(dr + b);
                  Store2_(am + b, Load2_(am + b) + dv * xav);
                  sv += dv * Load2_(xm + b);
               }
               if (b < len)
               {
                  am[b] += dr[b] * xa;
                  s0 += dr[b] * xm[b];
               }
            }
            acc[a] += s0 + (sv[0] + sv[1]);
            d += NM * len;
         }
         for (int f = own_[e]; f < own_[e + 1]; f++)
         {
            const int *ih = ix + static_cast<std::size_t>(fhi_[f]) * nb;
            for (int b = 0; b < nb; b++)
            {
               for (int m = 0; m < NM; m++) { xh[m][b] = x[m][ih[b]]; }
               yh[b] = 0.0;
            }
            const real_t *br = Bp + static_cast<std::size_t>(f) * nbb;
            for (int a = 0; a < nb; a++)
            {
               real_t sa = 0.0;
               for (int m = 0; m < NM; m++, br += nb) { sa += RowDotAxpy_<NBC>(nb, br, xh[m], xe[m][a], yh); }
               acc[a] += sa;
            }
            for (int b = 0; b < nb; b++) { y[ih[b]] += yh[b]; }
         }
         for (int a = 0; a < nb; a++) { y[ie[a]] = acc[a]; }
      }
      {
         ScopedStageTimer wait(timers_, StageTimers::DynHaloWait);
         MPI_Waitall(q, reqs_.data(), MPI_STATUSES_IGNORE);
      }
      const real_t *h = recv_buf_.data(), *Sp = S_.data();
      for (int s = 0; s < ns_; s++)
      {
         const int *ie = ix + static_cast<std::size_t>(se_[s]) * nb;
         const int *hc = hcol_.data() + static_cast<std::size_t>(s) * nb;
         for (int k = 0; k < nb; k++)
         {
            for (int m = 0; m < NM; m++) { xh[m][k] = h[NM * hc[k] + m]; }
         }
         const real_t *br = Sp + static_cast<std::size_t>(s) * nbb;
         for (int a = 0; a < nb; a++)
         {
            real_t sa = 0.0;
            for (int m = 0; m < NM; m++, br += nb) { sa += RowDot_<NBC>(nb, br, xh[m]); }
            y[ie[a]] += sa;
         }
      }
   }

   bool ok_ = false;
   int nm_ = 1, n_ = 0, ne_ = 0, nb_ = 0, tri_ = 0, nf_ = 0, ns_ = 0;
   real_t asym_[2] = {0.0, 0.0}, amax_[2] = {0.0, 0.0};
   std::vector<int> idx_;             ///< element DOFs, nb per element
   std::vector<int> own_, fhi_;       ///< faces of each lower element; the upper element
   std::vector<int> sptr_, se_, hcol_; ///< shared blocks per element; their element; halo slots
   std::vector<real_t> D_, B_, S_;    ///< diagonal (upper), face and shared blocks
   MPI_Comm comm_ = MPI_COMM_NULL;
   std::vector<int> send_procs_, send_starts_, send_map_, recv_procs_, recv_starts_;
   mutable std::vector<real_t> send_buf_, recv_buf_;
   mutable std::vector<MPI_Request> reqs_;
   StageTimers *timers_ = nullptr;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DENSE_BLOCK_STIFFNESS_HPP
