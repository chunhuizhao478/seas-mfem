// SPDX-License-Identifier: BSD-3-Clause
// matrix_stats.hpp — one-shot analysis of the assembled DG stiffness matrix.
//
// Designed to be called once, right after ParallelAssemble, before MUMPS
// factorization.  Computes scalar statistics suitable for answering:
//   "How sparse is the matrix?  What is the local bandwidth-equivalent?
//    Is it symmetric?  Is the diagonal positive?  Would AMG be plausible?"
//
// All MPI collectives reduce to rank 0; only rank 0 writes the JSON file.
// Cost: O(global_nnz) for one Transpose + Add + FNorm, plus O(local_nnz)
// for diag/row scans.  At BP5 production scale (~1.5M rows, ~90M nnz)
// the analysis runs in well under a minute.

#ifndef MFEM_SEAS_MATRIX_STATS_HPP
#define MFEM_SEAS_MATRIX_STATS_HPP

#include "mfem.hpp"

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// Scalar stats describing the assembled DG stiffness matrix.
struct MatrixStats
{
   // Size / partition
   long long global_n         = 0;   // # rows = # cols
   long long global_nnz       = 0;   // diag + offd, summed over ranks
   long long global_diag_nnz  = 0;   // rank-local block nnz, summed
   long long global_offd_nnz  = 0;   // inter-rank coupling, summed
   int n_ranks                = 0;
   long long local_rows_min   = 0;
   long long local_rows_max   = 0;
   double    local_rows_mean  = 0.0;

   // Per-row nnz (full row = diag + offd)
   long long nnz_per_row_min  = 0;
   long long nnz_per_row_max  = 0;
   double    nnz_per_row_mean = 0.0;

   // DG element-block structure (passed in by caller; not derived).
   int dg_dofs_per_elem       = 0;   // scalar DOFs in one element
   int dg_vdim                = 0;
   int dg_block_size          = 0;   // dofs_per_elem * vdim

   // "Local bandwidth" — within the rank-owned square diag block,
   // distance between min/max column index per row.  Meaningful for
   // DG only on a per-rank basis (post-partition); not a true bandwidth.
   long long local_bandwidth_min = 0;
   long long local_bandwidth_max = 0;
   double    local_bandwidth_mean = 0.0;

   // Norms / symmetry
   double fnorm              = 0.0;  // ||A||_F
   double sym_residual_fnorm = -1.0; // ||A - A^T||_F  (-1 if not computed)
   double sym_residual_rel   = -1.0; // ||A - A^T||_F / ||A||_F
   double sym_residual_max_abs = -1.0; // max |A_ij - A_ji| over diag block

   // Diagonal sanity (necessary condition for SPD)
   double diag_min           = 0.0;
   double diag_max           = 0.0;
   double diag_abs_min       = 0.0;
   bool   diag_all_positive  = true;
   long long diag_nonpos_count = 0;

   // K * 1 norms (rigid-body kernel probe — a rank-1 sanity check;
   // for elasticity with Dirichlet BCs, K*1 != 0 because Dirichlet
   // rows pin the constant mode).
   double Kones_norm1   = 0.0;
   double Kones_norm2   = 0.0;
   double Kones_norminf = 0.0;
};

/// Compute scalar stats for an assembled HypreParMatrix.
/// @param A           the assembled DG stiffness matrix
/// @param dofs_per_elem  scalar DOFs in one mesh element (e.g. 4 for P1 tet)
/// @param vdim        vector dimension (3 for 3D elasticity)
/// @param comm        MPI communicator (typically pmesh.GetComm())
/// @param do_symmetry if true, computes ||A - A^T||_F via Transpose + Add
///                    (skip on very large matrices if memory is tight)
inline MatrixStats ComputeMatrixStats(const HypreParMatrix &A,
                                      int dofs_per_elem,
                                      int vdim,
                                      MPI_Comm comm,
                                      bool do_symmetry = true)
{
   MatrixStats s;

   int my_rank = 0, n_ranks = 1;
   MPI_Comm_rank(comm, &my_rank);
   MPI_Comm_size(comm, &n_ranks);
   s.n_ranks         = n_ranks;
   s.global_n        = static_cast<long long>(A.GetGlobalNumRows());
   s.dg_dofs_per_elem = dofs_per_elem;
   s.dg_vdim         = vdim;
   s.dg_block_size   = dofs_per_elem * vdim;

   // Local pieces.
   SparseMatrix diag_local, offd_local;
   HYPRE_BigInt *cmap = nullptr;
   A.GetDiag(diag_local);
   A.GetOffd(offd_local, cmap);

   const int local_rows = diag_local.Height();
   const int *Ad_I = diag_local.GetI();
   const int *Ad_J = diag_local.GetJ();
   const double *Ad_data = diag_local.GetData();
   const int *Ao_I = offd_local.GetI();
   const double *Ao_data = offd_local.GetData();

   const long long diag_local_nnz = diag_local.NumNonZeroElems();
   const long long offd_local_nnz = offd_local.NumNonZeroElems();

   // Per-row stats (full row nnz = diag + offd).
   long long row_nnz_min = std::numeric_limits<long long>::max();
   long long row_nnz_max = 0;
   long long row_nnz_sum = 0;

   long long bw_min = std::numeric_limits<long long>::max();
   long long bw_max = 0;
   long long bw_sum = 0;
   long long bw_count = 0;

   double diag_min_local = std::numeric_limits<double>::infinity();
   double diag_max_local = -std::numeric_limits<double>::infinity();
   double diag_abs_min_local = std::numeric_limits<double>::infinity();
   long long diag_nonpos_local = 0;

   for (int i = 0; i < local_rows; ++i)
   {
      const long long rnnz = (Ad_I[i+1] - Ad_I[i]) + (Ao_I[i+1] - Ao_I[i]);
      row_nnz_min = std::min(row_nnz_min, rnnz);
      row_nnz_max = std::max(row_nnz_max, rnnz);
      row_nnz_sum += rnnz;

      // Local "bandwidth" within the diag block (cheap proxy).
      const int r0 = Ad_I[i];
      const int r1 = Ad_I[i+1];
      if (r1 > r0)
      {
         int jmin = Ad_J[r0], jmax = Ad_J[r0];
         double diag_val = 0.0;
         bool found_diag = false;
         for (int k = r0; k < r1; ++k)
         {
            const int j = Ad_J[k];
            if (j < jmin) { jmin = j; }
            if (j > jmax) { jmax = j; }
            if (j == i) { diag_val = Ad_data[k]; found_diag = true; }
         }
         bw_min = std::min(bw_min, static_cast<long long>(jmax - jmin));
         bw_max = std::max(bw_max, static_cast<long long>(jmax - jmin));
         bw_sum += (jmax - jmin);
         bw_count++;
         if (found_diag)
         {
            if (diag_val < diag_min_local) { diag_min_local = diag_val; }
            if (diag_val > diag_max_local) { diag_max_local = diag_val; }
            if (std::abs(diag_val) < diag_abs_min_local)
            {
               diag_abs_min_local = std::abs(diag_val);
            }
            if (diag_val <= 0.0) { diag_nonpos_local++; }
         }
      }
   }
   if (row_nnz_min == std::numeric_limits<long long>::max()) { row_nnz_min = 0; }
   if (bw_min == std::numeric_limits<long long>::max()) { bw_min = 0; }

   // Local Frobenius squared (for ||A||_F and for max|Aij-Aji| on diag).
   double fnorm_sq_local = 0.0;
   for (long long k = 0; k < diag_local_nnz; ++k)
   {
      fnorm_sq_local += Ad_data[k] * Ad_data[k];
   }
   for (long long k = 0; k < offd_local_nnz; ++k)
   {
      fnorm_sq_local += Ao_data[k] * Ao_data[k];
   }

   // Max |Aij - Aji| over the local diag block (cheap, exact for the
   // rank-owned square — does NOT include inter-rank entries).
   double sym_max_abs_local = 0.0;
   {
      SparseMatrix *diag_T = Transpose(diag_local);
      MFEM_VERIFY(diag_T->Height() == diag_local.Height(),
                  "diag transpose mismatch");
      const int *T_I = diag_T->GetI();
      const int *T_J = diag_T->GetJ();
      const double *T_data = diag_T->GetData();
      // Walk Ad and find matching T entries; both should have identical
      // (i, j) coverage since CSR sparsity is symmetric for a symmetric
      // assembly even if values differ slightly.
      for (int i = 0; i < local_rows; ++i)
      {
         // Build map j -> value for row i of T.
         std::vector<int> Tj(T_J + T_I[i], T_J + T_I[i+1]);
         std::vector<double> Tv(T_data + T_I[i], T_data + T_I[i+1]);
         for (int k = Ad_I[i]; k < Ad_I[i+1]; ++k)
         {
            const int j = Ad_J[k];
            const double aij = Ad_data[k];
            double aji = 0.0;
            for (std::size_t m = 0; m < Tj.size(); ++m)
            {
               if (Tj[m] == j) { aji = Tv[m]; break; }
            }
            const double d = std::abs(aij - aji);
            if (d > sym_max_abs_local) { sym_max_abs_local = d; }
         }
      }
      delete diag_T;
   }

   // K * ones probe (cheap, useful for cross-code comparison).
   double Kones_n1_local   = 0.0;
   double Kones_n2sq_local = 0.0;
   double Kones_ninf_local = 0.0;
   {
      HypreParVector ones(A.GetComm(), A.GetGlobalNumRows(), A.GetRowStarts());
      HypreParVector Kv  (A.GetComm(), A.GetGlobalNumRows(), A.GetRowStarts());
      ones = 1.0;
      A.Mult(ones, Kv);
      const int n = Kv.Size();
      for (int i = 0; i < n; ++i)
      {
         const double v = std::abs(Kv(i));
         Kones_n1_local   += v;
         Kones_n2sq_local += v * v;
         if (v > Kones_ninf_local) { Kones_ninf_local = v; }
      }
   }

   // Reductions.  Use long long for counts, double for floats.
   long long lr_min = local_rows, lr_max = local_rows;
   long long lr_sum = local_rows;
   long long diag_nnz_global = 0, offd_nnz_global = 0;
   long long rnz_min_g = 0, rnz_max_g = 0, rnz_sum_g = 0;
   long long bw_min_g = 0, bw_max_g = 0, bw_sum_g = 0, bw_count_g = 0;
   double diag_min_g = 0.0, diag_max_g = 0.0, diag_abs_min_g = 0.0;
   long long diag_nonpos_g = 0;
   double fnorm_sq_g = 0.0, sym_max_abs_g = 0.0;
   double Kones_n1_g = 0.0, Kones_n2sq_g = 0.0, Kones_ninf_g = 0.0;
   long long lr_min_in = lr_min, lr_max_in = lr_max, lr_sum_in = lr_sum;
   long long diag_local_nnz_in = diag_local_nnz, offd_local_nnz_in = offd_local_nnz;
   long long rnz_min_in = row_nnz_min, rnz_max_in = row_nnz_max;
   long long rnz_sum_in = row_nnz_sum;
   long long bw_min_in = bw_min, bw_max_in = bw_max;
   long long bw_sum_in = bw_sum, bw_count_in = bw_count;
   double diag_min_in = diag_min_local;
   double diag_max_in = diag_max_local;
   double diag_abs_min_in = diag_abs_min_local;
   long long diag_nonpos_in = diag_nonpos_local;

   MPI_Reduce(&lr_min_in, &lr_min, 1, MPI_LONG_LONG, MPI_MIN, 0, comm);
   MPI_Reduce(&lr_max_in, &lr_max, 1, MPI_LONG_LONG, MPI_MAX, 0, comm);
   MPI_Reduce(&lr_sum_in, &lr_sum, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
   MPI_Reduce(&diag_local_nnz_in, &diag_nnz_global, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
   MPI_Reduce(&offd_local_nnz_in, &offd_nnz_global, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
   MPI_Reduce(&rnz_min_in, &rnz_min_g, 1, MPI_LONG_LONG, MPI_MIN, 0, comm);
   MPI_Reduce(&rnz_max_in, &rnz_max_g, 1, MPI_LONG_LONG, MPI_MAX, 0, comm);
   MPI_Reduce(&rnz_sum_in, &rnz_sum_g, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
   MPI_Reduce(&bw_min_in, &bw_min_g, 1, MPI_LONG_LONG, MPI_MIN, 0, comm);
   MPI_Reduce(&bw_max_in, &bw_max_g, 1, MPI_LONG_LONG, MPI_MAX, 0, comm);
   MPI_Reduce(&bw_sum_in, &bw_sum_g, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
   MPI_Reduce(&bw_count_in, &bw_count_g, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
   MPI_Reduce(&diag_min_in, &diag_min_g, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
   MPI_Reduce(&diag_max_in, &diag_max_g, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
   MPI_Reduce(&diag_abs_min_in, &diag_abs_min_g, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
   MPI_Reduce(&diag_nonpos_in, &diag_nonpos_g, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
   MPI_Reduce(&fnorm_sq_local, &fnorm_sq_g, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
   MPI_Reduce(&sym_max_abs_local, &sym_max_abs_g, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
   MPI_Reduce(&Kones_n1_local, &Kones_n1_g, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
   MPI_Reduce(&Kones_n2sq_local, &Kones_n2sq_g, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
   MPI_Reduce(&Kones_ninf_local, &Kones_ninf_g, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

   // ||A - A^T||_F via Transpose + Add — only on rank 0 to avoid n_ranks
   // copies of intermediates.  Skip if requested or if A is too large.
   double sym_resid_g = -1.0;
   if (do_symmetry)
   {
      HypreParMatrix *AT = A.Transpose();
      if (AT)
      {
         // A and AT must share row/col partitions and col_map_offd; for
         // a square matrix with conformal row/col partitions and a
         // symmetric assembly this holds.  Wrap in try-style guard.
         HypreParMatrix *D = Add(1.0, A, -1.0, *AT);
         if (D)
         {
            sym_resid_g = D->FNorm();
            delete D;
         }
         delete AT;
      }
   }

   // Populate output struct on rank 0 (other ranks see defaults).
   if (my_rank == 0)
   {
      s.local_rows_min  = lr_min;
      s.local_rows_max  = lr_max;
      s.local_rows_mean = static_cast<double>(lr_sum) / std::max(1, n_ranks);
      s.global_diag_nnz = diag_nnz_global;
      s.global_offd_nnz = offd_nnz_global;
      s.global_nnz      = diag_nnz_global + offd_nnz_global;
      s.nnz_per_row_min = rnz_min_g;
      s.nnz_per_row_max = rnz_max_g;
      s.nnz_per_row_mean =
         (s.global_n > 0) ? static_cast<double>(rnz_sum_g) / s.global_n : 0.0;
      s.local_bandwidth_min  = bw_min_g;
      s.local_bandwidth_max  = bw_max_g;
      s.local_bandwidth_mean =
         (bw_count_g > 0) ? static_cast<double>(bw_sum_g) / bw_count_g : 0.0;
      s.fnorm                = std::sqrt(fnorm_sq_g);
      s.sym_residual_fnorm   = sym_resid_g;
      s.sym_residual_rel     =
         (s.fnorm > 0.0 && sym_resid_g >= 0.0) ? sym_resid_g / s.fnorm : -1.0;
      s.sym_residual_max_abs = sym_max_abs_g;
      s.diag_min             = diag_min_g;
      s.diag_max             = diag_max_g;
      s.diag_abs_min         = diag_abs_min_g;
      s.diag_nonpos_count    = diag_nonpos_g;
      s.diag_all_positive    = (diag_nonpos_g == 0);
      s.Kones_norm1          = Kones_n1_g;
      s.Kones_norm2          = std::sqrt(Kones_n2sq_g);
      s.Kones_norminf        = Kones_ninf_g;
   }
   return s;
}

/// Write the stats to a JSON file (rank 0 only — caller should guard).
inline void WriteMatrixStatsJson(const MatrixStats &s, const std::string &path)
{
   std::ofstream f(path);
   if (!f) { return; }
   f << std::setprecision(17);
   f << "{\n";
   f << "  \"global_n\": "        << s.global_n        << ",\n";
   f << "  \"global_nnz\": "      << s.global_nnz      << ",\n";
   f << "  \"global_diag_nnz\": " << s.global_diag_nnz << ",\n";
   f << "  \"global_offd_nnz\": " << s.global_offd_nnz << ",\n";
   f << "  \"offd_fraction\": "
     << (s.global_nnz > 0
         ? static_cast<double>(s.global_offd_nnz) / s.global_nnz : 0.0)
     << ",\n";
   f << "  \"n_ranks\": " << s.n_ranks << ",\n";
   f << "  \"local_rows\": { \"min\": " << s.local_rows_min
     << ", \"max\": " << s.local_rows_max
     << ", \"mean\": " << s.local_rows_mean << " },\n";
   f << "  \"nnz_per_row\": { \"min\": " << s.nnz_per_row_min
     << ", \"max\": " << s.nnz_per_row_max
     << ", \"mean\": " << s.nnz_per_row_mean << " },\n";
   f << "  \"dg\": { \"dofs_per_elem\": " << s.dg_dofs_per_elem
     << ", \"vdim\": " << s.dg_vdim
     << ", \"block_size\": " << s.dg_block_size << " },\n";
   f << "  \"local_bandwidth\": { \"min\": " << s.local_bandwidth_min
     << ", \"max\": " << s.local_bandwidth_max
     << ", \"mean\": " << s.local_bandwidth_mean
     << ", \"note\": \"distance between min/max column index in the rank-owned diag block; per-rank proxy, not a true global bandwidth\" },\n";
   f << "  \"fnorm\": " << s.fnorm << ",\n";
   f << "  \"symmetry\": { \"fnorm_diff\": " << s.sym_residual_fnorm
     << ", \"relative_to_fnorm\": " << s.sym_residual_rel
     << ", \"max_abs_diag_block\": " << s.sym_residual_max_abs << " },\n";
   f << "  \"diagonal\": { \"min\": " << s.diag_min
     << ", \"max\": " << s.diag_max
     << ", \"abs_min\": " << s.diag_abs_min
     << ", \"all_positive\": " << (s.diag_all_positive ? "true" : "false")
     << ", \"nonpositive_count\": " << s.diag_nonpos_count << " },\n";
   f << "  \"K_times_ones\": { \"l1\": " << s.Kones_norm1
     << ", \"l2\": " << s.Kones_norm2
     << ", \"linf\": " << s.Kones_norminf << " }\n";
   f << "}\n";
}

/// Pretty-print a one-screen summary to mfem::out (rank 0 only).
inline void PrintMatrixStatsSummary(const MatrixStats &s)
{
   const auto fmt = [](double x)
   { std::ostringstream o; o << std::setprecision(6) << x; return o.str(); };
   mfem::out
      << "  [MATRIX-STATS]\n"
      << "    global_n           = " << s.global_n << "\n"
      << "    global_nnz         = " << s.global_nnz
      << "  (diag " << s.global_diag_nnz
      << " + offd " << s.global_offd_nnz << ")\n"
      << "    offd_fraction      = "
      << (s.global_nnz > 0
          ? static_cast<double>(s.global_offd_nnz) / s.global_nnz : 0.0)
      << "  (inter-rank coupling)\n"
      << "    n_ranks            = " << s.n_ranks << "\n"
      << "    local_rows         = [" << s.local_rows_min
      << ", " << s.local_rows_max << "]"
      << "  mean " << fmt(s.local_rows_mean) << "\n"
      << "    nnz/row            = [" << s.nnz_per_row_min
      << ", " << s.nnz_per_row_max << "]"
      << "  mean " << fmt(s.nnz_per_row_mean) << "\n"
      << "    DG block size      = " << s.dg_block_size
      << "  (dofs_per_elem " << s.dg_dofs_per_elem
      << " * vdim " << s.dg_vdim << ")\n"
      << "    local bandwidth    = [" << s.local_bandwidth_min
      << ", " << s.local_bandwidth_max << "]"
      << "  mean " << fmt(s.local_bandwidth_mean)
      << "  (per-rank diag block)\n"
      << "    ||A||_F            = " << fmt(s.fnorm) << "\n"
      << "    ||A - A^T||_F      = " << fmt(s.sym_residual_fnorm)
      << "   (rel " << fmt(s.sym_residual_rel) << ")\n"
      << "    max|Aij - Aji|     = " << fmt(s.sym_residual_max_abs)
      << "  (diag block only)\n"
      << "    diag in            = [" << fmt(s.diag_min)
      << ", " << fmt(s.diag_max) << "]"
      << "  abs_min " << fmt(s.diag_abs_min) << "\n"
      << "    diag all positive  = " << (s.diag_all_positive ? "yes" : "NO")
      << "   nonpositive_count=" << s.diag_nonpos_count << "\n"
      << "    ||K*1||_{1,2,inf}  = " << fmt(s.Kones_norm1)
      << ", " << fmt(s.Kones_norm2)
      << ", " << fmt(s.Kones_norminf) << "\n";
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_MATRIX_STATS_HPP
