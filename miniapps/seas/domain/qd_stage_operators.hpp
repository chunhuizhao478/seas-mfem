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

// qd_stage_operators.hpp — the per-stage face work of a quasi-dynamic stage
// as sparse operators (GPU plan G2).  With δ the owned fault slip ([dip,
// strike] per owned fault DOF), u the displacement (true DOFs) and g(t) the
// Dirichlet boundary values at stored face quadrature points, one stage is
//
//   b   = B·δ + D·g(t)          right-hand side of K u = b
//   T   = C_u·u + C_δ·δ         shear traction, owned fault DOFs [dip, strike]
//   T_n = N_u·u + N_δ·δ         normal traction, owned fault DOFs
//                               (compression positive, as ComputeTraction)
//
// in place of the face loops AssembleSlipContributionIP(+Shared),
// AssembleDirichletLoading and ComputeTraction + RestrictToOwnedFault.
// ElasticityDomainOperator builds the matrices once by probing its
// face-local functions with unit inputs (elasticity_operator_qd.inl).  This
// file holds what it builds them into: row-owned triplets, their conversion
// to HypreParMatrix, the Dirichlet point set, and the products.  hypre does
// the halo exchange that ExpandOwnedToLocalFault and ExchangeFaceNbrData do
// by hand in the loops.

#ifndef MFEM_SEAS_QD_STAGE_OPERATORS_HPP
#define MFEM_SEAS_QD_STAGE_OPERATORS_HPP

#include "mfem.hpp"

#ifdef MFEM_USE_MPI

#include <algorithm>
#include <iomanip>
#include <memory>
#include <numeric>
#include <ostream>
#include <sstream>
#include <vector>

namespace mfem
{
namespace seas
{

/// Partition of a distributed index space with @a n_local entries on this
/// rank, in the layout the HypreParMatrix constructors take: [first, end)
/// with HYPRE_AssumedPartitionCheck(), otherwise every rank's offset.
/// Returns the global size; @a first is this rank's first global index.
inline HYPRE_BigInt MakeQDPartition(MPI_Comm comm, HYPRE_BigInt n_local,
                                    Array<HYPRE_BigInt> &part,
                                    HYPRE_BigInt &first)
{
   int nranks, rank;
   MPI_Comm_size(comm, &nranks);
   MPI_Comm_rank(comm, &rank);
   std::vector<HYPRE_BigInt> all(nranks);
   MPI_Allgather(&n_local, 1, HYPRE_MPI_BIG_INT, all.data(), 1,
                 HYPRE_MPI_BIG_INT, comm);
   std::vector<HYPRE_BigInt> offsets(nranks + 1, 0);
   for (int r = 0; r < nranks; r++) { offsets[r + 1] = offsets[r] + all[r]; }
   first = offsets[rank];
   if (HYPRE_AssumedPartitionCheck())
   {
      part.SetSize(2);
      part[0] = offsets[rank];
      part[1] = offsets[rank + 1];
   }
   else
   {
      part.SetSize(nranks + 1);
      for (int r = 0; r <= nranks; r++) { part[r] = offsets[r]; }
   }
   return offsets[nranks];
}

/// Row-owned sparse entries (local row, global column, value), collected
/// face by face while probing and turned into a HypreParMatrix once.
class QDTriplets
{
public:
   /// Exact zeros (no coupling) are not stored.
   void Add(int row, HYPRE_BigInt col, real_t val)
   {
      if (val == 0.0) { return; }
      row_.push_back(row);
      col_.push_back(col);
      val_.push_back(val);
   }

   /// Sum duplicate (row, col) entries and build the matrix: @a nrows local
   /// rows, row and column partitions as MakeQDPartition returns them.
   std::unique_ptr<HypreParMatrix> Build(MPI_Comm comm, int nrows,
                                         HYPRE_BigInt glob_nrows,
                                         const Array<HYPRE_BigInt> &row_part,
                                         HYPRE_BigInt glob_ncols,
                                         const Array<HYPRE_BigInt> &col_part)
   const
   {
      const size_t n = row_.size();
      std::vector<size_t> perm(n);
      std::iota(perm.begin(), perm.end(), size_t(0));
      std::sort(perm.begin(), perm.end(), [this](size_t a, size_t b)
      {
         if (row_[a] != row_[b]) { return row_[a] < row_[b]; }
         if (col_[a] != col_[b]) { return col_[a] < col_[b]; }
         return a < b;  // equal entries summed in insertion order
      });
      Array<int> I(nrows + 1);
      I = 0;
      std::vector<HYPRE_BigInt> J;
      std::vector<real_t> V;
      J.reserve(n);
      V.reserve(n);
      for (size_t k = 0; k < n; k++)
      {
         const size_t e = perm[k];
         MFEM_VERIFY(row_[e] >= 0 && row_[e] < nrows,
                     "QDTriplets: row " << row_[e] << " outside [0, "
                     << nrows << ")");
         MFEM_VERIFY(col_[e] >= 0 && col_[e] < glob_ncols,
                     "QDTriplets: column " << col_[e] << " outside [0, "
                     << glob_ncols << ")");
         if (k > 0 && row_[e] == row_[perm[k - 1]] &&
             col_[e] == col_[perm[k - 1]])
         {
            V.back() += val_[e];
            continue;
         }
         J.push_back(col_[e]);
         V.push_back(val_[e]);
         I[row_[e] + 1]++;
      }
      for (int r = 0; r < nrows; r++) { I[r + 1] += I[r]; }
      // An empty rank still passes valid pointers.
      if (J.empty()) { J.push_back(0); V.push_back(0.0); }
      return std::unique_ptr<HypreParMatrix>(new HypreParMatrix(
         comm, nrows, glob_nrows, glob_ncols, I.GetData(), J.data(),
         V.data(), row_part.GetData(), col_part.GetData()));
   }

   /// Shift every column by @a offset (columns collected as local indices).
   void OffsetColumns(HYPRE_BigInt offset)
   {
      for (HYPRE_BigInt &c : col_) { c += offset; }
   }

   size_t Size() const { return row_.size(); }

private:
   std::vector<int> row_;
   std::vector<HYPRE_BigInt> col_;
   std::vector<real_t> val_;
};

/// The stage operators and the Dirichlet point set of one domain operator.
struct QDStageOperators
{
   bool with_normal = false;       ///< N_u, N_δ built (elastic sigma_n)
   std::unique_ptr<HypreParMatrix> B;    ///< volume DOFs × owned slip
   std::unique_ptr<HypreParMatrix> D;    ///< volume DOFs × g(t) values
   std::unique_ptr<HypreParMatrix> Cu;   ///< owned shear traction × volume
   std::unique_ptr<HypreParMatrix> Cd;   ///< owned shear traction × slip
   std::unique_ptr<HypreParMatrix> Nu;   ///< owned normal traction × volume
   std::unique_ptr<HypreParMatrix> Nd;   ///< owned normal traction × slip

   /// Dirichlet points on this rank: boundary attribute (selects the
   /// loading function) and physical coordinates; g(t) holds 3 values per
   /// point, point p at entries 3p, 3p+1, 3p+2.
   Array<int> g_attr;
   Vector g_xyz;                   ///< [3 * points]
   Vector g;                       ///< g(t) of the current stage
   int NumPoints() const { return g_attr.Size(); }

   double build_seconds = 0.0;     ///< probing + matrix construction
   long long g_evals = 0;          ///< g(t) evaluations (one per stage)
   double g_eval_seconds = 0.0;    ///< time in the loading callback

   /// Sizes, nonzeros and memory of every matrix, summed over ranks.
   /// Collective; prints on rank 0.
   void PrintSummary(MPI_Comm comm, std::ostream &os) const
   {
      struct Entry { const char *name; const HypreParMatrix *A; };
      const Entry mats[] = {{"B", B.get()}, {"D", D.get()},
                            {"C_u", Cu.get()}, {"C_delta", Cd.get()},
                            {"N_u", Nu.get()}, {"N_delta", Nd.get()}};
      int rank;
      MPI_Comm_rank(comm, &rank);
      double local_pts = NumPoints(), glob_pts = 0.0, bmax = 0.0;
      MPI_Reduce(&local_pts, &glob_pts, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
      MPI_Reduce(&build_seconds, &bmax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
      double total_bytes = 0.0;
      std::ostringstream o;
      for (const Entry &m : mats)
      {
         if (m.A == nullptr) { continue; }
         hypre_ParCSRMatrix *h = *m.A;
         hypre_CSRMatrix *dg = hypre_ParCSRMatrixDiag(h);
         hypre_CSRMatrix *od = hypre_ParCSRMatrixOffd(h);
         // nnz, row pointers, offd column map: bytes as hypre stores them
         double loc[2];
         loc[0] = double(hypre_CSRMatrixNumNonzeros(dg)) +
                  double(hypre_CSRMatrixNumNonzeros(od));
         loc[1] = loc[0] * (sizeof(real_t) + sizeof(HYPRE_Int)) +
                  2.0 * (m.A->Height() + 1) * sizeof(HYPRE_Int) +
                  double(hypre_CSRMatrixNumCols(od)) * sizeof(HYPRE_BigInt);
         double glob[2];
         MPI_Reduce(loc, glob, 2, MPI_DOUBLE, MPI_SUM, 0, comm);
         total_bytes += glob[1];
         if (rank == 0)
         {
            o << "[QD-OPERATORS] " << std::left << std::setw(8) << m.name
              << std::right << std::setw(12) << m.A->GetGlobalNumRows()
              << " x " << std::setw(12) << m.A->GetGlobalNumCols()
              << "  nnz " << std::setw(12) << static_cast<long long>(glob[0])
              << "  " << std::fixed << std::setprecision(1)
              << glob[1] / 1048576.0 << " MiB\n";
         }
      }
      if (rank == 0)
      {
         o << "[QD-OPERATORS] Dirichlet points " << static_cast<long long>(
              glob_pts) << " (g(t) size " << static_cast<long long>(
              3 * glob_pts) << "); all matrices " << std::fixed
           << std::setprecision(1) << total_bytes / 1048576.0
           << " MiB; build " << std::setprecision(3) << bmax
           << " s (slowest rank)\n";
         os << o.str();
         os.flush();
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_SEAS_QD_STAGE_OPERATORS_HPP
