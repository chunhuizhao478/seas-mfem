// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// bound_constrained_qp.hpp — min 1/2 x'Ax - b'x subject to lo <= x <= hi for
// a parallel SPD matrix A: the damage problem of the staggered scheme when
// its residual is linear in d (p = 2 degradation, AT1 or AT2).
//
// Reduced-space active-set method, the same idea as PETSc vinewtonrsls
// (which MOOSE/RACCOON use through `-snes_type vinewtonrsls`):
//
//   r = A x - b
//   active   = { i : (x_i <= lo_i and r_i > 0) or (x_i >= hi_i and r_i < 0) }
//   inactive = the rest; converged when ||r_inactive|| <= max(rtol ||r_0||, atol)
//   with the active set unchanged from the previous iterate, or on two
//   consecutive iterates
//   step: solve A_II dx_I = -r_I (dx_A = 0) with PCG + BoomerAMG,
//         x <- P(x + s dx), s halved until the energy decreases.
//
// The norms are of the inactive residual (the reduced-space residual PETSc
// monitors).  Errors are reported, never silently accepted: Solve returns
// false when it does not converge within max_it.

#ifndef MFEM_PFF_BOUND_CONSTRAINED_QP_HPP
#define MFEM_PFF_BOUND_CONSTRAINED_QP_HPP

#include "mfem.hpp"

#include <memory>

namespace mfem
{
namespace pff
{

class BoundConstrainedQP
{
public:
   struct Stats
   {
      int iterations = 0;
      int linear_iterations = 0;
      int active = 0;
      real_t r0 = 0.0, r = 0.0;
   };

   real_t rel_tol = 1e-8;
   real_t abs_tol = 1e-10;
   real_t linear_rel_tol = 1e-12;
   int max_it = 50;
   int print_level = 0;

   /// x holds the initial guess on entry (projected onto the bounds).
   bool Solve(const HypreParMatrix &A, const Vector &b, const Vector &lo, const Vector &hi,
              Vector &x, Stats &st) const
   {
      MPI_Comm comm = A.GetComm();
      const int n = x.Size();
      for (int i = 0; i < n; i++) { x(i) = std::min(std::max(x(i), lo(i)), hi(i)); }
      Vector r(n), dx(n), xt(n), rhs(n);
      Array<int> active, prev_active;
      st = Stats();
      bool prev_tol_met = false;
      for (int it = 0; it <= max_it; it++)
      {
         Residual(A, b, x, r);
         active.SetSize(0);
         for (int i = 0; i < n; i++)
         {
            if ((x(i) <= lo(i) && r(i) > 0.0) || (x(i) >= hi(i) && r(i) < 0.0))
            {
               active.Append(i);
            }
         }
         const real_t rn = InactiveNorm(comm, r, active);
         if (it == 0) { st.r0 = rn; }
         st.r = rn;
         st.iterations = it;
         st.active = GlobalCount(comm, active.Size());
         if (print_level > 0 && Mpi::Root())
         {
            mfem::out << "    damage QP it " << it << " |r_I| = " << rn
                      << " active = " << st.active << "\n";
         }
         // Converged on a repeated active set, or when the tolerance is met on
         // two consecutive iterates: a node on a bound with r_i at round-off
         // can enter and leave the set forever (cmame hydromech, step 644:
         // |r_I| = 1.5e-16 with the count alternating 14237/14238).
         const bool tol_met = rn <= std::max(rel_tol * st.r0, abs_tol);
         if (tol_met && (SameSet(comm, active, prev_active) || prev_tol_met))
         {
            return true;
         }
         prev_tol_met = tol_met;
         if (it == max_it) { break; }
         prev_active = active;

         // Reduced Newton step: A_II dx_I = -r_I, dx_A = 0.
         std::unique_ptr<HypreParMatrix> Ar(new HypreParMatrix(A));
         rhs = r;
         rhs.Neg();
         for (int i : active) { rhs(i) = 0.0; }
         std::unique_ptr<HypreParMatrix> Ae(Ar->EliminateRowsCols(active));
         HypreBoomerAMG amg(*Ar);
         amg.SetPrintLevel(0);
         CGSolver cg(comm);
         cg.SetOperator(*Ar);
         cg.SetPreconditioner(amg);
         cg.SetRelTol(linear_rel_tol);
         cg.SetAbsTol(0.0);
         cg.SetMaxIter(2000);
         cg.SetPrintLevel(0);
         dx = 0.0;
         cg.Mult(rhs, dx);
         if (!cg.GetConverged())
         {
            MFEM_ABORT("BoundConstrainedQP: PCG did not converge (final norm "
                       << cg.GetFinalNorm() << " after " << cg.GetNumIterations()
                       << " iterations)");
         }
         st.linear_iterations += cg.GetNumIterations();

         // Projected step with energy backtracking.
         const real_t e0 = Energy(A, b, x);
         real_t s = 1.0;
         bool accepted = false;
         for (int ls = 0; ls < 30; ls++)
         {
            for (int i = 0; i < n; i++)
            {
               xt(i) = std::min(std::max(x(i) + s * dx(i), lo(i)), hi(i));
            }
            if (Energy(A, b, xt) <= e0 + 1e-14 * std::abs(e0)) { accepted = true; break; }
            s *= 0.5;
         }
         if (!accepted)
         {
            MFEM_ABORT("BoundConstrainedQP: line search failed at iteration " << it);
         }
         x = xt;
      }
      return false;
   }

private:
   static void Residual(const HypreParMatrix &A, const Vector &b, const Vector &x, Vector &r)
   {
      A.Mult(x, r);
      r -= b;
   }

   static real_t Energy(const HypreParMatrix &A, const Vector &b, const Vector &x)
   {
      Vector Ax(x.Size());
      A.Mult(x, Ax);
      return 0.5 * InnerProduct(A.GetComm(), x, Ax) - InnerProduct(A.GetComm(), b, x);
   }

   static real_t InactiveNorm(MPI_Comm comm, const Vector &r, const Array<int> &active)
   {
      real_t s = 0.0;
      int k = 0;
      for (int i = 0; i < r.Size(); i++)
      {
         if (k < active.Size() && active[k] == i) { k++; continue; }
         s += r(i) * r(i);
      }
      real_t g = 0.0;
      MPI_Allreduce(&s, &g, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
      return std::sqrt(g);
   }

   static int GlobalCount(MPI_Comm comm, int n)
   {
      int g = 0;
      MPI_Allreduce(&n, &g, 1, MPI_INT, MPI_SUM, comm);
      return g;
   }

   static bool SameSet(MPI_Comm comm, const Array<int> &a, const Array<int> &b)
   {
      int same = (a.Size() == b.Size()) ? 1 : 0;
      for (int i = 0; same && i < a.Size(); i++) { if (a[i] != b[i]) { same = 0; } }
      int g = 0;
      MPI_Allreduce(&same, &g, 1, MPI_INT, MPI_MIN, comm);
      return g == 1;
   }
};

} // namespace pff
} // namespace mfem

#endif
