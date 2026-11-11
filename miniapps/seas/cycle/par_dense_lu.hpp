// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// cycle/par_dense_lu.hpp — a dense LU shared by all ranks (ScaLAPACK).
//
// Why (Chunhui Zhao, 2026-11-10: "I would like to emphasize profiling and
// speedup for the current workflow"): the implicit QD fault solve
// (cycle/qd_implicit.hpp, option I) factors its reduced Newton matrix, size
// M = N + m (N fault points, m open points), with a serial LU on every rank:
// every rank does the same 2/3 M^3 flops on one core.  On the UCERF model
// (M = 3,290, 144 ranks) that was 5 factors of about 0.8 s in each 4.8 s QD
// step (Vista job 1133708).  ParDenseLU factors the same matrix once, with
// the ranks sharing the work (ScaLAPACK pdgetrf on a 2-D block-cyclic
// layout), and solves with pdgetrs.
//
// What goes in and out:
//   Factor(M, entry)   entry(i, j) returns the matrix entry (0-based, the
//                      same function on every rank: the callers hold the
//                      replicated Green's function).  Each rank evaluates
//                      only the entries it stores.  Returns false when the
//                      matrix is singular to working precision.  Collective.
//   Solve(b, x)        b, x: the full right-hand side and solution, length
//                      M, the same values on every rank (the callers gather
//                      them already).  Collective.
//
// The layout: nprow x npcol process grid, nprow = floor(sqrt(P)), npcol =
// P / nprow (ranks beyond nprow * npcol hold nothing and only receive the
// solution); square blocks of kBlock.  The pivoting is partial pivoting as
// in the serial LU, so the solution agrees with it to round-off (not bit
// for bit: the order of the floating-point operations differs and depends
// on the rank count).
//
// ScaLAPACK comes with MUMPS (MFEM_USE_MUMPS); without it Available() is
// false and Factor stops with a message.

#ifndef MFEM_SEAS_CYCLE_PAR_DENSE_LU_HPP
#define MFEM_SEAS_CYCLE_PAR_DENSE_LU_HPP

#include "mfem.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

#if defined(MFEM_USE_MPI) && defined(MFEM_USE_MUMPS) && defined(MFEM_USE_DOUBLE)
#define SEAS_HAVE_SCALAPACK
extern "C"
{
   int Csys2blacs_handle(MPI_Comm comm);
   void Cblacs_gridinit(int *ctxt, const char *order, int nprow, int npcol);
   void Cblacs_gridinfo(int ctxt, int *nprow, int *npcol, int *myrow, int *mycol);
   void Cblacs_gridexit(int ctxt);
   int Cfree_blacs_system_handle(int handle);
   int numroc_(const int *n, const int *nb, const int *iproc, const int *isrcproc,
               const int *nprocs);
   void descinit_(int *desc, const int *m, const int *n, const int *mb, const int *nb,
                  const int *irsrc, const int *icsrc, const int *ictxt, const int *lld,
                  int *info);
   void pdgetrf_(const int *m, const int *n, double *a, const int *ia, const int *ja,
                 const int *desca, int *ipiv, int *info);
   void pdgetrs_(const char *trans, const int *n, const int *nrhs, const double *a,
                 const int *ia, const int *ja, const int *desca, const int *ipiv, double *b,
                 const int *ib, const int *jb, const int *descb, int *info);
}
#endif

namespace mfem
{
namespace seas
{

class ParDenseLU
{
public:
   static constexpr int kBlock = 64;   ///< block size of the block-cyclic layout

   static bool Available()
   {
#ifdef SEAS_HAVE_SCALAPACK
      return true;
#else
      return false;
#endif
   }

   explicit ParDenseLU(MPI_Comm comm) : comm_(comm)
   {
#ifdef SEAS_HAVE_SCALAPACK
      int P = 1;
      MPI_Comm_size(comm_, &P);
      nprow_ = std::max(1, static_cast<int>(std::floor(std::sqrt(static_cast<double>(P)))));
      npcol_ = P / nprow_;
      handle_ = Csys2blacs_handle(comm_);
      ctxt_ = handle_;
      Cblacs_gridinit(&ctxt_, "R", nprow_, npcol_);
      int pr = 0, pc = 0;
      Cblacs_gridinfo(ctxt_, &pr, &pc, &myrow_, &mycol_);
      in_grid_ = (myrow_ >= 0 && mycol_ >= 0);
#endif
   }

   ~ParDenseLU()
   {
#ifdef SEAS_HAVE_SCALAPACK
      if (in_grid_) { Cblacs_gridexit(ctxt_); }
      Cfree_blacs_system_handle(handle_);
#endif
   }

   ParDenseLU(const ParDenseLU &) = delete;
   ParDenseLU &operator=(const ParDenseLU &) = delete;

   int GridRows() const { return nprow_; }
   int GridCols() const { return npcol_; }

   /// Assemble this rank's blocks from entry(i, j) and factor.  Collective.
   template <class Entry>
   bool Factor(int M, Entry &&entry)
   {
#ifdef SEAS_HAVE_SCALAPACK
      M_ = M;
      int info = 0;
      if (in_grid_)
      {
         const int nb = kBlock, zero = 0, one = 1;
         locr_ = numroc_(&M, &nb, &myrow_, &zero, &nprow_);
         locc_ = numroc_(&M, &nb, &mycol_, &zero, &npcol_);
         lld_ = std::max(1, locr_);
         descinit_(desca_, &M, &M, &nb, &nb, &zero, &zero, &ctxt_, &lld_, &info);
         MFEM_VERIFY(info == 0, "ParDenseLU: descinit (matrix) info " << info);
         descinit_(descb_, &M, &one, &nb, &nb, &zero, &zero, &ctxt_, &lld_, &info);
         MFEM_VERIFY(info == 0, "ParDenseLU: descinit (right-hand side) info " << info);
         a_.resize(size_t(lld_) * std::max(1, locc_));
         ipiv_.resize(locr_ + nb);
         rows_.resize(locr_);
         for (int li = 0; li < locr_; li++) { rows_[li] = Global_(li, myrow_, nprow_); }
         for (int lj = 0; lj < locc_; lj++)
         {
            const int gj = Global_(lj, mycol_, npcol_);
            double *col = a_.data() + size_t(lj) * lld_;
            for (int li = 0; li < locr_; li++) { col[li] = entry(rows_[li], gj); }
         }
         pdgetrf_(&M, &M, a_.data(), &one, &one, desca_, ipiv_.data(), &info);
         MFEM_VERIFY(info >= 0, "ParDenseLU: pdgetrf argument " << -info);
      }
      // every rank learns the outcome (ranks outside the grid too)
      int bad = (info > 0) ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, comm_);
      ok_ = !bad;
      return ok_;
#else
      (void)M; (void)entry;
      MFEM_ABORT("ParDenseLU: this build has no ScaLAPACK (it comes with MFEM_USE_MUMPS)");
      return false;
#endif
   }

   /// x = A^-1 b; b and x hold all M entries on every rank.  Collective.
   void Solve(const double *b, double *x) const
   {
#ifdef SEAS_HAVE_SCALAPACK
      MFEM_VERIFY(ok_, "ParDenseLU::Solve: no factor");
      std::vector<double> part(M_, 0.0);
      if (in_grid_ && mycol_ == 0)
      {
         std::vector<double> rhs(std::max(1, locr_));
         for (int li = 0; li < locr_; li++) { rhs[li] = b[rows_[li]]; }
         const int one = 1;
         int info = 0;
         pdgetrs_("N", &M_, &one, a_.data(), &one, &one, desca_, ipiv_.data(), rhs.data(),
                  &one, &one, descb_, &info);
         MFEM_VERIFY(info == 0, "ParDenseLU: pdgetrs info " << info);
         for (int li = 0; li < locr_; li++) { part[rows_[li]] = rhs[li]; }
      }
      else if (in_grid_)
      {
         // the right-hand side lives in process column 0; the others take
         // part in the solve with an empty piece
         std::vector<double> rhs(1, 0.0);
         const int one = 1;
         int info = 0;
         pdgetrs_("N", &M_, &one, a_.data(), &one, &one, desca_, ipiv_.data(), rhs.data(),
                  &one, &one, descb_, &info);
         MFEM_VERIFY(info == 0, "ParDenseLU: pdgetrs info " << info);
      }
      // each entry comes from exactly one rank; the others add 0.0
      MPI_Allreduce(part.data(), x, M_, MPI_DOUBLE, MPI_SUM, comm_);
#else
      (void)b; (void)x;
      MFEM_ABORT("ParDenseLU: this build has no ScaLAPACK (it comes with MFEM_USE_MUMPS)");
#endif
   }

private:
   /// Global index of local index l on process p of np (source process 0).
   static int Global_(int l, int p, int np)
   {
      return ((l / kBlock) * np + p) * kBlock + l % kBlock;
   }

   MPI_Comm comm_;
   int nprow_ = 1, npcol_ = 1, myrow_ = -1, mycol_ = -1, handle_ = -1, ctxt_ = -1;
   bool in_grid_ = false, ok_ = false;
   int M_ = 0, locr_ = 0, locc_ = 0, lld_ = 1;
   int desca_[9] = {0}, descb_[9] = {0};
   std::vector<double> a_;
   std::vector<int> ipiv_, rows_;
};

} // namespace seas
} // namespace mfem

#endif
