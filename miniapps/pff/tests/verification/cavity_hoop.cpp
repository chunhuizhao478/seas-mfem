// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// cavity_hoop.cpp — C1 gate 2: the hoop stress around the pressurised
// borehole matches the Lame solution (farms@cdms pulsepower/spe_hftc/PLAN.md,
// Appendix A.8, which uses the infinite-medium limit sigma_tt = p rb^2/r^2).
//
//   mpirun -np 2 ./cavity_hoop [-m tests/verification/meshes/annulus_rb4p15_R16p6_h0p3.msh]
//
// An acrylic annulus r_b <= r <= R (r_b = 4.15 mm, R = 4 r_b, free outer
// edge, plane strain) under a slowly rising borehole pressure: a smooth ramp
// sin^2(pi t / 2 T) to P = 100 MPa over T = 300 us >> R / c_s = 18 us, then
// held, with a small Kelvin-Voigt damping (zeta = 1e-8 s).
// Explicit pff step (cfl 0.4). At t = 400 us the element stresses with
// centroids in r in [1.1 r_b, 0.9 R] are compared with the thick-cylinder
// (Lame) solution for p(t_r):
//   sigma_rr = A - B / r^2, sigma_tt = A + B / r^2,
//   A = p r_b^2 / (R^2 - r_b^2), B = p r_b^2 R^2 / (R^2 - r_b^2),
// as the RMS of |sigma - sigma_Lame| / sigma_tt,Lame(r). The P1 element
// stresses converge at first order (static FE on this mesh: 3.5 % at h = 0.3 mm,
// 1.74 % at h = 0.15 mm), so the gates are: explicit against the static FE
// solution of the same mesh within 1 % (the explicit scheme reaches the
// quasi-static state), and the static FE against Lame within 5 % at h = 0.3 mm.
// -s skips the explicit run.

#include "mfem.hpp"
#include "../../io/gmsh_physical_names.hpp"
#include "../../domain/explicit_elasticity_operator.hpp"
#include "../../solver/central_difference.hpp"
#include "../../solver/cfl.hpp"

#include <cmath>
#include <memory>
#include <iostream>

using namespace mfem;
using namespace mfem::pff;

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const char *mesh_file = "tests/verification/meshes/annulus_rb4p15_R16p6_h0p3.msh";
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Annulus mesh (gmsh 2.2, borehole/outer).");
   bool static_solve = false;
   args.AddOption(&static_solve, "-s", "--static", "-ns", "--no-static",
                  "Diagnostic: solve K u = P b statically instead of the explicit run.");
   args.Parse();
   if (!args.Good()) { args.PrintUsage(mfem::out); return 1; }

   PhysicalNameMap names = ReadPhysicalNames(mesh_file);
   Mesh serial(mesh_file, 1, 1);
   const int max_attr = serial.bdr_attributes.Max();
   ParMesh mesh(MPI_COMM_WORLD, serial);
   serial.Clear();
   H1_FECollection h1(1, 2);
   ParFiniteElementSpace V(&mesh, &h1, 2, Ordering::byNODES), S(&mesh, &h1);

   ExplicitElasticityParams p;
   const real_t E = 2.8e9, nu = 0.35, rb = 4.15e-3, R = 16.6e-3, zeta = 1e-8;
   const real_t T_ramp = 3e-4, T_end = 4e-4;
   p.K = E / 3.0 / (1.0 - 2.0 * nu);
   p.G = E / 2.0 / (1.0 + nu);
   p.density = 1190.0;
   p.decomposition = Decomposition::Spectral;
   p.pressure_marker.SetSize(max_attr); p.pressure_marker = 0;
   p.absorbing_marker.SetSize(max_attr); p.absorbing_marker = 0;
   for (int a : ResolveAttributes({"borehole"}, names, 1)) { p.pressure_marker[a - 1] = 1; }
   PulseLoadParams pp;
   pp.peak_pressure = 1e8;
   PulseLoad pulse(pp);   // only Peak(x) = P is used; the time history is the ramp below

   // Rigid modes: u_y = 0 at (+-R, 0) and u_x = 0 at (0, R), which the
   // axisymmetric solution satisfies (u = u_r(r) e_r).
   Array<int> fixed;
   {
      const real_t pts[3][2] = {{R, 0.0}, {-R, 0.0}, {0.0, R}};
      const int comps[3] = {1, 1, 0};
      for (int k = 0; k < 3; k++)
      {
         for (int v = 0; v < mesh.GetNV(); v++)
         {
            const real_t *x = mesh.GetVertex(v);
            if (std::hypot(x[0] - pts[k][0], x[1] - pts[k][1]) > 1e-9) { continue; }
            Array<int> vd;
            V.GetVertexVDofs(v, vd);
            const int t = V.GetLocalTDofNumber(vd[comps[k]]);
            if (t >= 0) { fixed.Append(t); }
         }
      }
      int nloc = fixed.Size(), nglob = 0;
      MPI_Allreduce(&nloc, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MFEM_VERIFY(nglob == 3, "expected 3 pinned dofs at (+-R,0), (0,R), got " << nglob);
   }
   ExplicitElasticityOperator op(V, S, p, pulse, fixed);
   const int n = V.GetTrueVSize(), nt = S.GetTrueVSize();
   const real_t omega = std::sqrt(op.MaxEigenvalue(60));
   const real_t dt0 = 0.4 * CriticalStep(omega, zeta);
   const int steps = static_cast<int>(std::ceil(T_end / dt0));
   const real_t dt = T_end / steps;
   // Smooth ramp s(t) = sin^2(pi t / (2 T_ramp)) to 1, then held: zero load
   // rate at both ends, so the ring modes are barely excited.
   auto ramp = [&](real_t t)
   {
      if (t >= T_ramp) { return 1.0; }
      const real_t q = std::sin(0.5 * M_PI * t / T_ramp);
      return q * q;
   };

   CentralDifferenceState st;
   st.u.SetSize(n); st.u = 0.0;
   st.v_half.SetSize(n); st.v_half = 0.0;
   st.F_prev.SetSize(n); st.F_prev = 0.0;
   Vector D(nt);
   D = 0.0;
   CentralDifference cd(op);
   CentralDifferenceStepInfo info;
   // Explicit run to t_r (skipped with -s).
   real_t P = pp.peak_pressure;
   if (!static_solve)
   {
      for (int k = 0; k < steps; k++) { cd.Step(st, D, ramp(st.t), zeta, dt, info); }
      P = pp.peak_pressure * ramp(st.t);
   }

   // Static solution of the same discretisation for the same p: K u = p/P_peak b.
   Vector us(n);
   {
      ParBilinearForm k(&V);
      ConstantCoefficient lam(p.K - 2.0 * p.G / 3.0), mu(p.G);
      k.AddDomainIntegrator(new ElasticityIntegrator(lam, mu));
      k.Assemble();
      k.Finalize();
      std::unique_ptr<HypreParMatrix> Km(k.ParallelAssemble());
      std::unique_ptr<HypreParMatrix> Ke(Km->EliminateRowsCols(fixed));
      HypreBoomerAMG amg(*Km);
      amg.SetSystemsOptions(2, true);
      amg.SetPrintLevel(0);
      CGSolver cg(MPI_COMM_WORLD);
      cg.SetOperator(*Km);
      cg.SetPreconditioner(amg);
      cg.SetRelTol(1e-12);
      cg.SetMaxIter(2000);
      Vector bP(op.PressureVector());
      bP *= P / pp.peak_pressure;
      for (int t : fixed) { bP(t) = 0.0; }
      us = 0.0;
      cg.Mult(bP, us);
      MFEM_VERIFY(cg.GetConverged(), "static CG did not converge");
   }

   // Element-centroid polar stresses of a displacement field.
   const real_t lambda = p.K - 2.0 * p.G / 3.0;
   auto stresses = [&](const Vector &U, std::vector<real_t> &rr_out, std::vector<real_t> &tt_out,
                       std::vector<real_t> &r_out)
   {
      ParGridFunction u(&V);
      u.SetFromTrueDofs(U);
      DenseMatrix grad(2);
      IntegrationPoint c;
      c.Set2(1.0 / 3.0, 1.0 / 3.0);
      rr_out.clear(); tt_out.clear(); r_out.clear();
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         ElementTransformation &T = *mesh.GetElementTransformation(e);
         T.SetIntPoint(&c);
         Vector x(2);
         T.Transform(c, x);
         const real_t r = x.Norml2();
         if (r < 1.1 * rb || r > 0.9 * R) { continue; }
         u.GetVectorGradient(T, grad);
         const real_t exx = grad(0, 0), eyy = grad(1, 1), exy = 0.5 * (grad(0, 1) + grad(1, 0));
         const real_t tr = exx + eyy;   // plane strain
         const real_t sxx = lambda * tr + 2.0 * p.G * exx, syy = lambda * tr + 2.0 * p.G * eyy;
         const real_t sxy = 2.0 * p.G * exy;
         const real_t cs = x(0) / r, sn = x(1) / r;
         rr_out.push_back(sxx * cs * cs + syy * sn * sn + 2.0 * sxy * sn * cs);
         tt_out.push_back(sxx * sn * sn + syy * cs * cs - 2.0 * sxy * sn * cs);
         r_out.push_back(r);
      }
   };
   std::vector<real_t> drr, dtt, srr, stt, rad;
   stresses(st.u, drr, dtt, rad);
   stresses(us, srr, stt, rad);
   // sums: [dyn-Lame tt, dyn-Lame rr, static-Lame tt, static-Lame rr, dyn-static tt, dyn-static rr, count]
   real_t loc[7] = {0, 0, 0, 0, 0, 0, 0}, glob[7];
   for (size_t k = 0; k < rad.size(); k++)
   {
      const real_t r = rad[k];
      const real_t A = P * rb * rb / (R * R - rb * rb), B = A * R * R;
      const real_t tt = A + B / (r * r), rr = A - B / (r * r), w = 1.0 / (tt * tt);
      loc[0] += w * (dtt[k] - tt) * (dtt[k] - tt);
      loc[1] += w * (drr[k] - rr) * (drr[k] - rr);
      loc[2] += w * (stt[k] - tt) * (stt[k] - tt);
      loc[3] += w * (srr[k] - rr) * (srr[k] - rr);
      loc[4] += w * (dtt[k] - stt[k]) * (dtt[k] - stt[k]);
      loc[5] += w * (drr[k] - srr[k]) * (drr[k] - srr[k]);
      loc[6] += 1.0;
   }
   MPI_Allreduce(loc, glob, 7, MPITypeMap<real_t>::mpi_type, MPI_SUM, MPI_COMM_WORLD);
   real_t rms[6];
   for (int k = 0; k < 6; k++) { rms[k] = std::sqrt(glob[k] / glob[6]); }
   const bool ok = static_solve ? (rms[2] < 0.05 && rms[3] < 0.05)
                   : (rms[4] < 0.01 && rms[5] < 0.01 && rms[2] < 0.05 && rms[3] < 0.05);
   if (Mpi::Root())
   {
      mfem::out << "cavity_hoop: p = " << P << " Pa";
      if (!static_solve) { mfem::out << " at t = " << T_end << " after " << steps << " steps (dt = " << dt << ")"; }
      mfem::out << ", " << glob[6] << " elements in r in [1.1 r_b, 0.9 R]\n"
                << "  RMS error / sigma_tt,Lame (tt, rr):\n";
      if (!static_solve)
      {
         mfem::out << "    explicit vs Lame   " << rms[0] << ", " << rms[1] << "\n"
                   << "    explicit vs static " << rms[4] << ", " << rms[5] << "  (gate 1 %)\n";
      }
      mfem::out << "    static vs Lame     " << rms[2] << ", " << rms[3]
                << "  (gate 5 %; P1 centroid stresses, first order in h)\n"
                << (ok ? "PASS" : "FAIL") << "\n";
   }
   return ok ? 0 : 1;
}
