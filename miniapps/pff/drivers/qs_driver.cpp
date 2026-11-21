// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// qs_driver.cpp — quasi-static staggered phase-field fracture (phase R, the
// RACCOON tutorials).
//
//   mpirun -np N ./qs_driver -c raccoon/mode1_brittle_fracture/configs/mode1.toml
//                            [-ns NUM_STEPS] [-o OUTPUT_DIR]
//
// One load step, as the MOOSE MultiApp pair elasticity.i + fracture.i
// (TransientMultiApp at TIMESTEP_END, fixed-point iterations):
//
//   initial_norm = |R_u(u_old, d)| with the new Dirichlet values (MOOSE
//                  computeResidualL2Norm before iteration 0)
//   for k = 0 .. fixed_point_max_its-1:
//      u   = solve R_u(u; d) = 0                     (Newton, PCG + BoomerAMG)
//      D   = element averages of psie_active(u)       (the CONSTANT MONOMIAL transfer)
//      d   = bound-constrained damage solve, d_n <= d <= 1 (the sub-app is
//            restored to the start of the step at each iteration)
//      end = |R_u(u; d)|                              (norm after TIMESTEP_END MultiApps)
//      converged when end < fixed_point_abs_tol or end/initial_norm < fixed_point_rel_tol
//      (FixedPointSolve::examineFixedPointConvergence; the begin norm is 0 without
//      TIMESTEP_BEGIN MultiApps)
//
// Automatic scaling ([solver] automatic_scaling, MOOSE Executioner
// automatic_scaling with the defaults compute_scaling_once = true,
// resid_vs_jac_scaling_param = 0, off_diagonals_in_auto_scaling = false;
// NonlinearSystemBase::computeScaling at moose 94b6528): one factor per
// displacement component, 1 / max_i |J_ii| over the component's dofs, where J
// is the volumetric Jacobian assembled without constraints (no Dirichlet
// rows, hanging-node constraints not applied) at the first solve (d = 0).
// MOOSE multiplies the kernel residual by the factor; the Dirichlet rows
// (NodalBC set()) hold u - g unscaled.  The fixed-point norms (initial and
// end) are taken of this scaled residual (FEProblemBase::computeResidualL2Norm).
// The displacement Newton test stays unscaled: it is tighter than the scaled
// one, and MOOSE solves the linear system with LU.
// The reaction is the NodalSum of the save_in residual: the sum over the
// boundary's true dofs of the unconstrained internal force at the converged u
// of the last pass, with the d that solve used (the postprocessor runs at
// TIMESTEP_END before the MultiApp).
//
// Output: <dir>/<name>_reaction.csv (time, reaction, fixed-point iterations,
// damage QP iterations), and for the steps in [output] nodal_steps,
// <dir>/<name>_nodal_step<N>.txt with x y u_x u_y d per vertex (np = 1).

#include "mfem.hpp"
#include "../config/qs_config_parser.hpp"
#include "../domain/qs_mesh.hpp"
#include "../domain/qs_damage_operator.hpp"
#include "../integrator/split_elasticity_integrator.hpp"
#include "../solver/bound_constrained_qp.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <iomanip>
#include <memory>
#include <filesystem>

using namespace mfem;
using namespace mfem::pff;

namespace
{

struct EssentialBC
{
   Array<int> tdofs;
   Array<int> spec;   // index into QsCaseConfig::dirichlet per tdof
};

real_t GlobalNorm(MPI_Comm comm, const Vector &v)
{
   return std::sqrt(InnerProduct(comm, v, v));
}

/// Newton for R_u(u; d) = 0 on the true dofs, with X already holding the
/// Dirichlet values.  Converged when |R| <= max(nl_rel_tol |R0|, nl_abs_tol), or
/// when the update is at round-off (|dX| <= 1e-13 |X|).  Returns the iteration
/// count; aborts on failure.
/// Newton solve of the displacement as MOOSE/PETSc SNES runs it (solve_type
/// NEWTON, no line-search change): the residual has the kernel rows scaled by
/// the automatic scaling factors W and the Dirichlet rows u - g (NodalBC,
/// unscaled); a BC with preset = false starts from the old value, so the
/// first step also moves the boundary. Converged when
/// |F| <= max(nl_rel_tol |F_0|, nl_abs_tol). Row scaling does not change the
/// Newton step. Returns the number of Newton iterations.
int SolveDisplacement(ParNonlinearForm &Nfree, Vector &X, const Array<int> &ess,
                      const Vector &G, const Vector &W, const QsCaseConfig &c, int dim,
                      MPI_Comm comm)
{
   Vector R(X.Size()), dX(X.Size()), rhs(X.Size()), Xb(X.Size());
   auto residual = [&]()
   {
      Nfree.Mult(X, R);
      for (int t : ess) { R(t) = X(t) - G(t); }
      Vector Rw(R);
      for (int i = 0; i < Rw.Size(); i++) { Rw(i) *= W(i); }
      return GlobalNorm(comm, Rw);
   };
   const real_t r0 = residual();
   real_t r = r0;
   for (int it = 0; it < 50; it++)
   {
      if (r <= std::max(c.nl_rel_tol * r0, c.nl_abs_tol)) { return it; }
      HypreParMatrix J(dynamic_cast<HypreParMatrix &>(Nfree.GetGradient(X)));
      std::unique_ptr<HypreParMatrix> Je(J.EliminateRowsCols(ess));
      // Rows: J_ff dX_f + J_fb dX_b = R_f, dX_b = R_b = X_b - g_b.
      Xb = 0.0;
      for (int t : ess) { Xb(t) = R(t); }
      rhs = R;
      J.EliminateBC(*Je, ess, Xb, rhs);
      HypreBoomerAMG amg(J);
      amg.SetSystemsOptions(dim, true);
      amg.SetPrintLevel(0);
      CGSolver cg(comm);
      cg.SetOperator(J);
      cg.SetPreconditioner(amg);
      cg.SetRelTol(c.linear_rel_tol);
      cg.SetAbsTol(0.0);
      cg.SetMaxIter(5000);
      cg.SetPrintLevel(0);
      dX = 0.0;
      cg.Mult(rhs, dX);
      if (!cg.GetConverged())
      {
         MFEM_ABORT("displacement PCG did not converge: final norm " << cg.GetFinalNorm()
                    << " after " << cg.GetNumIterations() << " iterations");
      }
      X -= dX;
      r = residual();
      // Round-off floor: a residual tolerance below machine precision (the
      // RACCOON gold tests use nl_abs_tol = 1e-16) cannot be met; stop once
      // the Newton update is at round-off relative to the solution.
      if (GlobalNorm(comm, dX) <= 1e-13 * GlobalNorm(comm, X)) { return it + 1; }
   }
   MFEM_ABORT("displacement Newton did not converge: |R| = " << r << ", |R0| = " << r0);
   return -1;
}

/// MOOSE automatic scaling factors (see the header), from the serial mesh: its
/// element-wise assembled diagonal over all local dofs is libMesh's
/// unconstrained scaling Jacobian.  Identical on every rank.
std::vector<real_t> AutomaticScalingFactors(Mesh &serial, const QsCaseConfig &c,
                                            const IntegrationRule &ir)
{
   // The diagonal of the constrained (true-dof) Jacobian: MOOSE assembles the
   // scaling Jacobian through libMesh's constraint handling, so a node next
   // to a hanging node collects the constrained rows. M-R1b (verbose run,
   // 2026-11-21) prints 1.54762e-06 for mode I, = 1/max diag of the
   // constrained matrix here; the unconstrained diagonal gives 2.0635e-06.
   const int dim = serial.Dimension();
   H1_FECollection h1(1, dim);
   FiniteElementSpace sV(&serial, &h1, dim, Ordering::byNODES);
   FiniteElementSpace sD(&serial, &h1);
   GridFunction d0(&sD);
   d0 = 0.0;
   NonlinearForm nf(&sV);
   auto *integ = new SplitElasticityIntegrator(c.K(), c.G(), c.decomposition, c.pf.deg, d0);
   integ->SetIntRule(&ir);
   nf.AddDomainIntegrator(integ);
   Vector z(sV.GetTrueVSize());
   z = 0.0;
   Vector diag;
   dynamic_cast<SparseMatrix &>(nf.GetGradient(z)).GetDiag(diag);
   MFEM_VERIFY(diag.Size() == sV.GetTrueVSize(), "unexpected Jacobian size");
   const int nt = diag.Size() / dim;   // true vdofs byNODES: (comp, i) -> comp nt + i
   std::vector<real_t> factor(dim, 0.0);
   for (int comp = 0; comp < dim; comp++)
   {
      real_t m = 0.0;
      for (int i = 0; i < nt; i++) { m = std::max(m, std::abs(diag(comp * nt + i))); }
      factor[comp] = (m == 0.0) ? 1.0 : 1.0 / m;
   }
   return factor;
}

void WriteNodal(const std::string &path, ParMesh &mesh, ParGridFunction &u,
                ParGridFunction &d)
{
   Vector ux, uy, dv;
   u.GetNodalValues(ux, 1);
   u.GetNodalValues(uy, 2);
   d.GetNodalValues(dv, 1);
   std::ofstream o(path);
   MFEM_VERIFY(o.good(), "cannot write " << path);
   o << std::setprecision(17) << "# x y disp_x disp_y d\n";
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      const real_t *x = mesh.GetVertex(i);
      o << x[0] << ' ' << x[1] << ' ' << ux(i) << ' ' << uy(i) << ' ' << dv(i) << '\n';
   }
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;

   const char *config_file = "";
   int num_steps_override = -1;
   const char *output_override = "";
   bool fp_verbose = false;
   OptionsParser args(argc, argv);
   args.AddOption(&config_file, "-c", "--config", "Case TOML file.");
   args.AddOption(&num_steps_override, "-ns", "--num-steps", "Override [time] num_steps.");
   args.AddOption(&output_override, "-o", "--output-dir", "Override [output] dir.");
   args.AddOption(&fp_verbose, "-v", "--fp-verbose", "-nv", "--no-fp-verbose",
                  "Print the fixed-point norms of every pass.");
   args.Parse();
   if (!args.Good() || std::string(config_file).empty())
   {
      if (Mpi::Root()) { args.PrintUsage(mfem::out); }
      return 1;
   }

   QsCaseConfig c = QsCaseConfigParser::ParseFile(config_file);
   if (num_steps_override >= 0) { c.num_steps = num_steps_override; }
   if (std::string(output_override).size()) { c.output_dir = output_override; }
   int nsteps = c.num_steps;
   if (nsteps <= 0) { nsteps = static_cast<int>(std::lround(c.end_time / c.dt)); }

   // Mesh and spaces.
   BoundaryNames names;
   Mesh serial = BuildQsMesh(c, names);
   const int max_attr = serial.bdr_attributes.Max();
   const int dim = serial.Dimension();
   const IntegrationRule &ir = IntRules.Get(serial.GetElementBaseGeometry(0),
                                            c.quadrature_order);
   std::vector<real_t> scaling(dim, 1.0);
   if (c.automatic_scaling) { scaling = AutomaticScalingFactors(serial, c, ir); }
   ParMesh pmesh(comm, serial);
   serial.Clear();
   const long long ne = pmesh.GetGlobalNE();
   H1_FECollection h1(1, dim);
   L2_FECollection l2(0, dim);
   ParFiniteElementSpace Vfes(&pmesh, &h1, dim, Ordering::byNODES);
   ParFiniteElementSpace Dfes(&pmesh, &h1);
   ParFiniteElementSpace Lfes(&pmesh, &l2);
   if (Mpi::Root())
   {
      mfem::out << "pff qs_driver: " << c.name << "\n  source: " << c.source
                << "\n  elements " << ne << ", displacement dofs " << Vfes.GlobalTrueVSize()
                << ", damage dofs " << Dfes.GlobalTrueVSize() << ", steps " << nsteps << "\n";
      if (c.automatic_scaling)
      {
         mfem::out << "  automatic scaling factors:" << std::setprecision(10);
         for (int comp = 0; comp < dim; comp++) { mfem::out << ' ' << scaling[comp]; }
         mfem::out << "\n";
      }
   }

   ParGridFunction u(&Vfes), d(&Dfes), D(&Lfes);
   u = 0.0;
   d = 0.0;
   D = 0.0;

   auto make_integrator = [&]()
   {
      auto *integ = new SplitElasticityIntegrator(c.K(), c.G(), c.decomposition, c.pf.deg, d);
      integ->SetIntRule(&ir);
      return integ;
   };
   ParNonlinearForm N(&Vfes);      // with the essential rows removed
   N.AddDomainIntegrator(make_integrator());
   ParNonlinearForm Nfree(&Vfes);  // unconstrained internal force, for the reaction
   SplitElasticityIntegrator *psi_integ = make_integrator();
   Nfree.AddDomainIntegrator(psi_integ);

   // Dirichlet true dofs.
   EssentialBC ess;
   for (size_t s = 0; s < c.dirichlet.size(); s++)
   {
      Array<int> marker = BoundaryMarker(names, c.dirichlet[s].boundaries, max_attr);
      Array<int> list;
      Vfes.GetEssentialTrueDofs(marker, list, c.dirichlet[s].component);
      for (int t : list) { ess.tdofs.Append(t); ess.spec.Append(static_cast<int>(s)); }
   }
   {
      Array<int> sorted(ess.tdofs);
      sorted.Sort();
      sorted.Unique();
      MFEM_VERIFY(sorted.Size() == ess.tdofs.Size(), "a dof has two Dirichlet conditions");
   }
   N.SetEssentialTrueDofs(ess.tdofs);
   // Weights of the fixed-point residual norm: the component's scaling factor
   // on the kernel rows, 1 on the Dirichlet rows (true dofs are byNODES).
   Vector W(Vfes.GetTrueVSize());
   {
      const int ntd = Vfes.GetTrueVSize() / dim;
      for (int comp = 0; comp < dim; comp++)
         for (int i = 0; i < ntd; i++) { W(comp * ntd + i) = scaling[comp]; }
      for (int t : ess.tdofs) { W(t) = 1.0; }
   }
   auto scaled_norm = [&](Vector &v)
   {
      for (int i = 0; i < v.Size(); i++) { v(i) *= W(i); }
      return GlobalNorm(comm, v);
   };
   Array<int> reaction_tdofs;
   {
      Array<int> marker = BoundaryMarker(names, c.reaction_boundaries, max_attr);
      Vfes.GetEssentialTrueDofs(marker, reaction_tdofs, c.reaction_component);
   }

   QsDamageOperator damage(Dfes, c.pf, &ir);
   BoundConstrainedQP qp;
   qp.rel_tol = c.nl_rel_tol;
   qp.abs_tol = c.nl_abs_tol;
   qp.linear_rel_tol = c.linear_rel_tol;

   Vector X(Vfes.GetTrueVSize()), R(Vfes.GetTrueVSize());
   Vector Dd(Dfes.GetTrueVSize()), Dold(Dfes.GetTrueVSize());
   Vector lo(Dfes.GetTrueVSize()), hi(Dfes.GetTrueVSize());
   u.GetTrueDofs(X);
   d.GetTrueDofs(Dd);
   hi = 1.0;

   std::ofstream csv;
   if (Mpi::Root())
   {
      std::filesystem::create_directories(c.output_dir);
      csv.open(c.output_dir + "/" + c.name + "_reaction.csv");
      MFEM_VERIFY(csv.good(), "cannot write " << c.output_dir << "/" << c.name << "_reaction.csv");
      csv << "time," << c.reaction_name << ",fixed_point_its,damage_qp_its\n";
      csv << std::setprecision(12) << 0.0 << ',' << 0.0 << ",0,0\n";
   }

   for (int n = 1; n <= nsteps; n++)
   {
      const real_t t = n * c.dt;
      Dold = Dd;
      lo = Dold;

      // Initial fixed-point norm: old u, new Dirichlet values.
      Vector Xbc(X);
      for (int k = 0; k < ess.tdofs.Size(); k++)
      {
         const DirichletSpec &s = c.dirichlet[ess.spec[k]];
         Xbc(ess.tdofs[k]) = s.value + s.rate * t;
      }
      Nfree.Mult(X, R);
      for (int k = 0; k < ess.tdofs.Size(); k++)
      {
         R(ess.tdofs[k]) = X(ess.tdofs[k]) - Xbc(ess.tdofs[k]);
      }
      const real_t initial_norm = scaled_norm(R);
      // preset BCs (MOOSE default) are applied to the initial guess; the
      // others (preset = false) are reached by the first Newton step.
      for (int k = 0; k < ess.tdofs.Size(); k++)
      {
         if (c.dirichlet[ess.spec[k]].preset) { X(ess.tdofs[k]) = Xbc(ess.tdofs[k]); }
      }

      int k = 0, qp_its = 0, newton_its = 0;
      real_t F_last = 0.0;
      bool converged = false;
      real_t end_norm = 0.0;
      for (k = 0; k < c.fixed_point_max_its; k++)
      {
         // The grid function d holds the damage of the previous pass (or of
         // the last step): the displacement solve uses it.
         newton_its += SolveDisplacement(Nfree, X, ess.tdofs, Xbc, W, c, dim, comm);
         u.SetFromTrueDofs(X);

         // Reaction (NodalSum of save_in): MOOSE evaluates the postprocessor
         // at TIMESTEP_END before the MultiApp runs, so it sees this d.
         Nfree.Mult(X, R);
         F_last = 0.0;
         for (int tdof : reaction_tdofs) { F_last += R(tdof); }

         // Driving energy: element averages of psie_active(u).
         for (int e = 0; e < pmesh.GetNE(); e++)
         {
            Array<int> vdofs;
            Vfes.GetElementVDofs(e, vdofs);
            Vector el;
            u.GetSubVector(vdofs, el);
            real_t avg;
            psi_integ->ElementAverageActiveEnergy(*Vfes.GetFE(e),
                                                  *Vfes.GetElementTransformation(e), el, avg);
            Array<int> ld;
            Lfes.GetElementDofs(e, ld);
            D(ld[0]) = avg;
         }

         // Damage.
         damage.Assemble(D);
         BoundConstrainedQP::Stats st;
         // The fracture sub-app is restored to the start of the step at each
         // pass, so its solve starts from d_n.
         Vector dguess(Dold);
         if (!qp.Solve(damage.A(), damage.b(), lo, hi, dguess, st))
         {
            MFEM_ABORT("step " << n << " fixed-point pass " << k
                       << ": damage QP did not converge, |r_I| = " << st.r
                       << " (|r_I,0| = " << st.r0 << ")");
         }
         qp_its += st.iterations;
         Dd = dguess;
         d.SetFromTrueDofs(Dd);

         N.Mult(X, R);
         end_norm = scaled_norm(R);
         if (fp_verbose && Mpi::Root())
         {
            mfem::out << "  step " << n << " pass " << k << ": end norm " << std::setprecision(6)
                      << end_norm << ", initial " << initial_norm << ", relative "
                      << end_norm / initial_norm << ", Newton its (step so far) " << newton_its
                      << ", damage QP its " << st.iterations
                      << ", |r_I| " << st.r << " / " << st.r0 << "\n";
         }
         if (end_norm < c.fixed_point_abs_tol ||
             end_norm / initial_norm < c.fixed_point_rel_tol)
         {
            converged = true;
            k++;
            break;
         }
      }
      if (!converged)
      {
         if (!c.accept_on_max_fixed_point_iteration)
         {
            MFEM_ABORT("step " << n << ": fixed point not converged in "
                       << c.fixed_point_max_its << " iterations, end norm " << end_norm
                       << ", initial " << initial_norm);
         }
         if (Mpi::Root())
         {
            mfem::out << "  step " << n << ": accepted at the maximum of "
                      << c.fixed_point_max_its << " fixed-point iterations (end norm "
                      << end_norm << ", relative " << end_norm / initial_norm << ")\n";
         }
      }

      // The reaction of the last pass, with the d of its displacement solve.
      real_t Fg = 0.0;
      MPI_Allreduce(&F_last, &Fg, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
      real_t dmax_loc = d.Max(), dmax = 0.0;
      MPI_Allreduce(&dmax_loc, &dmax, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);

      if (Mpi::Root())
      {
         mfem::out << "step " << n << " t = " << t << "  " << c.reaction_name << " = "
                   << std::setprecision(10) << Fg << "  fixed-point its " << k
                   << "  damage QP its " << qp_its << "  max d = " << dmax << "\n";
         csv << std::setprecision(12) << t << ',' << std::setprecision(12) << Fg << ','
             << k << ',' << qp_its << '\n';
      }
      for (int s : c.nodal_output_steps)
      {
         if (s == n)
         {
            MFEM_VERIFY(Mpi::WorldSize() == 1, "nodal output is serial (np = 1)");
            WriteNodal(c.output_dir + "/" + c.name + "_nodal_step" + std::to_string(n) + ".txt",
                       pmesh, u, d);
         }
      }
   }
   return 0;
}
