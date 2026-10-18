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

#ifndef MFEM_SEAS_FUSED_STIFFNESS_PAIR_HPP
#define MFEM_SEAS_FUSED_STIFFNESS_PAIR_HPP

// y = K u + E v in one sweep, for the elastic stiffness K and the
// Kelvin-Voigt operator E = K_eta of NitscheWaveOperator.  Profiling plan
// item 4 (Chunhui Zhao, 2026-10-17).
//
// K and E come from the same integrators on the same space, so their
// structural patterns agree; AssembleStiffness finalizes with skip_zeros = 1,
// which drops the few entries that cancel to an exact zero in one matrix and
// not the other.  Init() therefore stores the UNION pattern once (row
// pointers + column indices) with the two value arrays interleaved
// (K_ij, E_ij), a missing entry being an explicit 0.  One sweep reads each
// index once for both products, and u and v cross the partition boundary
// together in one message per neighbour (K's hypre communication package,
// copied).  After Init the caller may free both HypreParMatrix objects: the
// union arrays take 20 bytes per entry against 24 for K and E.
//
// The halo is the offd column map of K (or of E when E's contains K's);
// Init returns false when neither contains the other and the caller keeps
// two hypre SpMVs.

#include "mfem.hpp"

#include <algorithm>
#include <vector>

namespace mfem
{
namespace seas
{

class FusedStiffnessPair
{
public:
   /// Build the union arrays.  False (nothing usable) when neither offd
   /// column map contains the other on this rank; the caller reduces the
   /// flag over all ranks.
   bool Init(HypreParMatrix &K, HypreParMatrix &E)
   {
      ok_ = false;
      hypre_ParCSRMatrix *k = K, *e = E;
      n_ = hypre_CSRMatrixNumRows(hypre_ParCSRMatrixDiag(k));
      if (n_ != hypre_CSRMatrixNumRows(hypre_ParCSRMatrixDiag(e))) { return false; }
      // The halo comes from the matrix whose offd column map contains the
      // other's (normally both are equal); the other's offd columns are
      // renumbered into it.
      hypre_CSRMatrix *ko = hypre_ParCSRMatrixOffd(k), *eo = hypre_ParCSRMatrixOffd(e);
      std::vector<int> e_to_k, k_to_e;
      const bool e_in_k = SubsetMap_(e, k, e_to_k);
      const bool k_in_e = !e_in_k && SubsetMap_(k, e, k_to_e);
      if (!e_in_k && !k_in_e) { return false; }
      hypre_ParCSRMatrix *host = e_in_k ? k : e;
      extra_ = Union_(hypre_ParCSRMatrixDiag(k), nullptr, hypre_ParCSRMatrixDiag(e), nullptr,
                      n_, hypre_CSRMatrixNumCols(hypre_ParCSRMatrixDiag(k)), dI_, dJ_, dV_);
      extra_ += Union_(ko, e_in_k ? nullptr : k_to_e.data(), eo, e_in_k ? e_to_k.data() : nullptr,
                       n_, hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(host)), oI_, oJ_, oV_);
      nnz_k_ = hypre_CSRMatrixNumNonzeros(hypre_ParCSRMatrixDiag(k)) +
               hypre_CSRMatrixNumNonzeros(ko);
      // halo pattern of the host (hypre's matvec communication package)
      if (!hypre_ParCSRMatrixCommPkg(host)) { hypre_MatvecCommPkgCreate(host); }
      hypre_ParCSRCommPkg *pkg = hypre_ParCSRMatrixCommPkg(host);
      comm_ = hypre_ParCSRCommPkgComm(pkg);
      const HYPRE_Int ns = hypre_ParCSRCommPkgNumSends(pkg);
      const HYPRE_Int nr = hypre_ParCSRCommPkgNumRecvs(pkg);
      send_procs_.assign(hypre_ParCSRCommPkgSendProcs(pkg), hypre_ParCSRCommPkgSendProcs(pkg) + ns);
      send_starts_.assign(hypre_ParCSRCommPkgSendMapStarts(pkg),
                          hypre_ParCSRCommPkgSendMapStarts(pkg) + ns + 1);
      send_map_.assign(hypre_ParCSRCommPkgSendMapElmts(pkg),
                       hypre_ParCSRCommPkgSendMapElmts(pkg) + send_starts_[ns]);
      recv_procs_.assign(hypre_ParCSRCommPkgRecvProcs(pkg), hypre_ParCSRCommPkgRecvProcs(pkg) + nr);
      recv_starts_.assign(hypre_ParCSRCommPkgRecvVecStarts(pkg),
                          hypre_ParCSRCommPkgRecvVecStarts(pkg) + nr + 1);
      const HYPRE_Int nc = hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(host));
      MFEM_VERIFY(recv_starts_[nr] == nc, "FusedStiffnessPair: halo size " << recv_starts_[nr]
                  << " != offd columns " << nc);
      w_.resize(2 * static_cast<std::size_t>(n_));
      send_buf_.resize(2 * send_starts_[ns]);
      recv_buf_.resize(2 * nc);
      reqs_.resize(ns + nr);
      ok_ = true;
      return true;
   }

   bool Ok() const { return ok_; }
   /// Entries added by the union (present in one matrix only), this rank.
   long long ExtraEntries() const { return extra_; }
   long long EntriesK() const { return nnz_k_; }

   /// y = K u + E v (u, v, y: local true vectors).
   void Mult(const Vector &u, const Vector &v, Vector &y) const
   {
      MFEM_ASSERT(ok_, "FusedStiffnessPair::Mult before a successful Init");
      const real_t *pu = u.HostRead(), *pv = v.HostRead();
      real_t *py = y.HostWrite();
      const int ns = static_cast<int>(send_procs_.size());
      const int nr = static_cast<int>(recv_procs_.size());
      const int tag = 7321;
      int q = 0;
      for (int r = 0; r < nr; r++)
      {
         MPI_Irecv(recv_buf_.data() + 2 * recv_starts_[r], 2 * (recv_starts_[r + 1] - recv_starts_[r]),
                   MPITypeMap<real_t>::mpi_type, recv_procs_[r], tag, comm_, &reqs_[q++]);
      }
      for (int s = 0; s < ns; s++)
      {
         for (int i = send_starts_[s]; i < send_starts_[s + 1]; i++)
         {
            send_buf_[2 * i] = pu[send_map_[i]];
            send_buf_[2 * i + 1] = pv[send_map_[i]];
         }
         MPI_Isend(send_buf_.data() + 2 * send_starts_[s], 2 * (send_starts_[s + 1] - send_starts_[s]),
                   MPITypeMap<real_t>::mpi_type, send_procs_[s], tag, comm_, &reqs_[q++]);
      }
      // local block while the halo travels; u and v interleaved so that one
      // gather brings both u_j and v_j (one cache line instead of two)
      real_t *w = w_.data();
      for (int j = 0; j < n_; j++) { w[2 * j] = pu[j]; w[2 * j + 1] = pv[j]; }
      const int *I = dI_.data(), *J = dJ_.data();
      const real_t *V = dV_.data();
      for (int i = 0; i < n_; i++)
      {
         // separate sums for K u and E v: two independent add chains (one
         // chain of 2 multiply-adds per entry is latency bound)
         real_t sk = 0.0, se = 0.0;
         for (int p = I[i]; p < I[i + 1]; p++)
         {
            const int j = 2 * J[p];
            sk += V[2 * p] * w[j];
            se += V[2 * p + 1] * w[j + 1];
         }
         py[i] = sk + se;
      }
      MPI_Waitall(q, reqs_.data(), MPI_STATUSES_IGNORE);
      if (!oJ_.empty())
      {
         const int *Io = oI_.data(), *Jo = oJ_.data();
         const real_t *Vo = oV_.data(), *h = recv_buf_.data();
         for (int i = 0; i < n_; i++)
         {
            real_t sk = 0.0, se = 0.0;
            for (int p = Io[i]; p < Io[i + 1]; p++)
            {
               const int j = 2 * Jo[p];
               sk += Vo[2 * p] * h[j];
               se += Vo[2 * p + 1] * h[j + 1];
            }
            py[i] += sk + se;
         }
      }
   }

private:
   /// map[j] = position of sub's offd column j in host's (sorted) offd
   /// column map; false when some column is missing.
   static bool SubsetMap_(hypre_ParCSRMatrix *sub, hypre_ParCSRMatrix *host,
                          std::vector<int> &map)
   {
      const HYPRE_Int ns = hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(sub));
      const HYPRE_Int nh = hypre_CSRMatrixNumCols(hypre_ParCSRMatrixOffd(host));
      const HYPRE_BigInt *cs = hypre_ParCSRMatrixColMapOffd(sub);
      const HYPRE_BigInt *ch = hypre_ParCSRMatrixColMapOffd(host);
      map.resize(ns);
      for (HYPRE_Int j = 0; j < ns; j++)
      {
         const HYPRE_BigInt *it = std::lower_bound(ch, ch + nh, cs[j]);
         if (it == ch + nh || *it != cs[j]) { return false; }
         map[j] = static_cast<int>(it - ch);
      }
      return true;
   }

   /// Union of the patterns of a and b, row by row (a's columns first, then
   /// b's extra ones); values interleaved (a, b), a missing value = 0.
   /// amap / bmap renumber the columns of a / b (nullptr = identity).  Two
   /// passes (count, then fill) so the arrays are allocated once at their
   /// final size.  Returns the number of entries present in only one of the
   /// two.
   static long long Union_(hypre_CSRMatrix *a, const int *amap, hypre_CSRMatrix *b,
                           const int *bmap, int n, int ncol, std::vector<int> &I,
                           std::vector<int> &J, std::vector<real_t> &V)
   {
      const HYPRE_Int *Ia = hypre_CSRMatrixI(a), *Ja = hypre_CSRMatrixJ(a);
      const HYPRE_Int *Ib = hypre_CSRMatrixI(b), *Jb = hypre_CSRMatrixJ(b);
      const HYPRE_Real *Aa = hypre_CSRMatrixData(a), *Ab = hypre_CSRMatrixData(b);
      std::vector<int> pos(ncol, -1);
      // pass 1: row sizes of the union
      I.assign(n + 1, 0);
      for (int i = 0; i < n; i++)
      {
         int cnt = Ia[i + 1] - Ia[i];
         for (HYPRE_Int p = Ia[i]; p < Ia[i + 1]; p++) { pos[amap ? amap[Ja[p]] : Ja[p]] = i; }
         for (HYPRE_Int p = Ib[i]; p < Ib[i + 1]; p++)
         {
            const int c = bmap ? bmap[Jb[p]] : Jb[p];
            if (pos[c] != i) { pos[c] = i; cnt++; }
         }
         I[i + 1] = I[i] + cnt;
      }
      // pass 2: fill
      J.assign(I[n], 0);
      V.assign(2 * static_cast<std::size_t>(I[n]), 0.0);
      std::fill(pos.begin(), pos.end(), -1);
      std::vector<char> inb;
      long long extra = 0;
      for (int i = 0; i < n; i++)
      {
         const int start = I[i];
         int q = start;
         for (HYPRE_Int p = Ia[i]; p < Ia[i + 1]; p++)
         {
            const int c = amap ? amap[Ja[p]] : Ja[p];
            pos[c] = q;
            J[q] = c;
            V[2 * q] = Aa[p];
            q++;
         }
         const int na = q - start;
         inb.assign(na, 0);
         for (HYPRE_Int p = Ib[i]; p < Ib[i + 1]; p++)
         {
            const int c = bmap ? bmap[Jb[p]] : Jb[p];
            if (pos[c] >= start)
            {
               V[2 * pos[c] + 1] = Ab[p];
               if (pos[c] < start + na) { inb[pos[c] - start] = 1; }
            }
            else
            {
               pos[c] = q;
               J[q] = c;
               V[2 * q + 1] = Ab[p];
               q++;
               extra++;
            }
         }
         for (int k = 0; k < na; k++) { if (!inb[k]) { extra++; } }
         MFEM_ASSERT(q == I[i + 1], "FusedStiffnessPair: union row size");
         for (int p = start; p < q; p++) { pos[J[p]] = -1; }
      }
      return extra;
   }

   bool ok_ = false;
   int n_ = 0;
   long long extra_ = 0, nnz_k_ = 0;
   std::vector<int> dI_, dJ_, oI_, oJ_;
   std::vector<real_t> dV_, oV_;
   MPI_Comm comm_ = MPI_COMM_NULL;
   std::vector<int> send_procs_, send_starts_, send_map_, recv_procs_, recv_starts_;
   mutable std::vector<real_t> w_, send_buf_, recv_buf_;
   mutable std::vector<MPI_Request> reqs_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FUSED_STIFFNESS_PAIR_HPP
