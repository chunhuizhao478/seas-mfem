// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_pore_pressure.cpp — checks of multiphysics/poroelastic/pore_pressure_operator.hpp
// against closed forms, with the cmame_revision 2d_hydromech/benchmark
// parameters (farms@cdms elasticity_E1d25.i: E 50 GPa, nu 0.3, K_s 50 GPa,
// K_f 2.24 GPa, porosity 0.008, viscosity 1e-3), np = 1, triangles:
//   - coupling force (ElkPorousFlowEffectiveStressCoupling): for uniform p and
//     d, F_p . u = -alpha(d) p int div u, exact for a linear u;
//   - storage step without flow (k = 0): a uniform volumetric strain increment
//     e gives p = p_old - alpha M e (the undrained response of
//     ElkPorousFlowFullySaturatedMassTimeDerivative), fluid energy 1/2 p^2/M |Omega|;
//   - 1D pressure diffusion with a drained end (PorousFlowFullySaturatedDarcyBase):
//     p(x, t) = sum_m 4/((2m+1) pi) sin((2m+1) pi x / 2L) exp(-(2m+1)^2 pi^2 c t / 4L^2),
//     c = (k/mu) M (Terzaghi's series for the pressure, one-way drainage);
//     the energy lost equals the time-integrated dissipation to first order;
//   - the undrained stiffness ratio of the CFL bound at d = 0:
//     1 + alpha^2 M / (lambda + G) in 2D.

#include "mfem.hpp"
#include "../test_framework.hpp"
#include "../../multiphysics/poroelastic/pore_pressure_operator.hpp"

using namespace mfem;
using namespace mfem::pff;

static const real_t E_ = 50e9, nu_ = 0.3;

static PoroDamageParams Params()
{
   PoroDamageParams P;
   P.K = E_ / 3.0 / (1.0 - 2.0 * nu_);
   P.K_s = 50e9;
   P.K_f = 2.24e9;
   P.phi0 = 0.008;
   P.phi_lb = 0.008;
   P.phi_ub = 0.999;
   P.g_min = 1e-6;
   P.viscosity = 1e-3;
   P.perm_model = PermeabilityModel::None;
   P.k0 = 0.0;
   return P;
}

static void TestCouplingAndStorage()
{
   std::cout << "\n[coupling force and storage step]\n";
   const real_t Lx = 2.0, Ly = 1.0;
   Mesh serial = Mesh::MakeCartesian2D(6, 5, Element::TRIANGLE, true, Lx, Ly);
   ParMesh mesh(MPI_COMM_WORLD, serial);
   H1_FECollection h1(1, 2);
   ParFiniteElementSpace V(&mesh, &h1, 2, Ordering::byNODES), S(&mesh, &h1);
   const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2);
   const PoroDamageParams P = Params();
   Array<int> none;
   PorePressureOperator op(V, S, P, none, ir, PoroPreconditioner::Jacobi);

   // linear u = A x
   const real_t A[4] = {1e-4, 3e-5, -2e-5, 4e-4};   // row major
   VectorFunctionCoefficient uf(2, [&](const Vector &x, Vector &u)
   {
      u(0) = A[0] * x(0) + A[1] * x(1);
      u(1) = A[2] * x(0) + A[3] * x(1);
   });
   ParGridFunction ug(&V);
   ug.ProjectCoefficient(uf);
   Vector U(V.GetTrueVSize()), U0(V.GetTrueVSize());
   ug.GetTrueDofs(U);
   U0 = 0.0;
   const real_t e = A[0] + A[3], area = Lx * Ly;

   for (real_t d : {0.0, 0.5, 1.0})
   {
      Vector D(S.GetTrueVSize()), Pp(S.GetTrueVSize()), F(V.GetTrueVSize());
      D = d;
      const real_t p0 = 3e6;
      Pp = p0;
      F = 0.0;
      op.AddCouplingForce(D, Pp, F);
      const real_t alpha = BiotCoefficient(P, d);
      TEST_REL(F * U, -alpha * p0 * e * area, 1e-12, "F_p . u = -alpha p int div u");

      const PoroPoint pt = EvalPoroPoint(P, d);
      op.Reset();
      PorePressureOperator::StepInfo info;
      op.Step(D, U, U0, 1e-8, Pp, info);
      const real_t pex = p0 - pt.alpha / pt.inv_M * e;
      TEST_REL(Pp.Min(), pex, 1e-9, "storage step: min p = p_old - alpha M e");
      TEST_REL(Pp.Max(), pex, 1e-9, "storage step: max p = p_old - alpha M e");
      TEST_REL(info.fluid_energy, 0.5 * pt.inv_M * pex * pex * area, 1e-9,
               "fluid energy = 1/2 p^2/M |Omega|");
      TEST_NEAR(info.dissipation_rate, 0.0, 1e-30, "no dissipation without flow");
   }
}

static void TestDiffusion(PoroPreconditioner pc, const char *label)
{
   std::cout << "\n[1D pressure diffusion, drained end, " << label << "]\n";
   const real_t L = 1e-2, W = 5e-4;
   const int nx = 80;
   Mesh serial = Mesh::MakeCartesian2D(nx, 2, Element::TRIANGLE, true, L, W);
   ParMesh mesh(MPI_COMM_WORLD, serial);
   H1_FECollection h1(1, 2);
   ParFiniteElementSpace V(&mesh, &h1, 2, Ordering::byNODES), S(&mesh, &h1);
   const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2);
   PoroDamageParams P = Params();
   P.k0 = 1e-15;
   const PoroPoint pt = EvalPoroPoint(P, 0.0);
   const real_t c = pt.mobility / pt.inv_M;
   // drained at x = 0 (attribute 4 of MakeCartesian2D)
   Array<int> marker(serial.bdr_attributes.Max()), drained;
   marker = 0;
   marker[3] = 1;
   S.GetEssentialTrueDofs(marker, drained);
   PorePressureOperator op(V, S, P, drained, ir, pc);

   Vector D(S.GetTrueVSize()), U(V.GetTrueVSize()), Pp(S.GetTrueVSize());
   D = 0.0;
   U = 0.0;
   const real_t p0 = 1e6;
   Pp = p0;
   for (int t : drained) { Pp(t) = 0.0; }
   real_t E0, rate0;
   // energy of the uniform start (before the boundary value is imposed)
   {
      Vector Pu(Pp.Size());
      Pu = p0;
      op.Energies(D, Pu, E0, rate0);
   }
   const real_t T = 0.2 * L * L / c;   // T_v = c t / L^2 = 0.2
   const int steps = 400;
   const real_t dt = T / steps;
   PorePressureOperator::StepInfo info;
   real_t diss = 0.0;
   for (int k = 0; k < steps; k++)
   {
      op.Step(D, U, U, dt, Pp, info);
      diss += dt * info.dissipation_rate;
   }
   ParGridFunction pg(&S), xg(&S);
   pg.SetFromTrueDofs(Pp);
   FunctionCoefficient xf([](const Vector &x) { return x(0); });
   xg.ProjectCoefficient(xf);
   const real_t pi = M_PI;
   real_t err = 0.0;
   for (int i = 0; i < pg.Size(); i++)
   {
      real_t s = 0.0;
      for (int m = 0; m < 200; m++)
      {
         const real_t a = (2 * m + 1) * pi / (2.0 * L);
         s += 4.0 / ((2 * m + 1) * pi) * std::sin(a * xg(i)) * std::exp(-a * a * c * T);
      }
      err = std::max(err, std::abs(pg(i) / p0 - s));
   }
   std::cout << "  c = " << c << " m^2/s, T = " << T << " s, max |p - series|/p0 = " << err
             << ", PCG iterations (last step) " << info.iterations << "\n";
   TEST_ASSERT(err < 5e-3, "pressure profile matches the Terzaghi series at T_v = 0.2 (0.5 % of p0)");
   // Energy: E(T) + dissipated + outflow work = E0. The drained boundary is at
   // p = 0, so it does no work; backward Euler dissipates the remainder.
   TEST_NEAR(info.boundary_power, 0.0, 1e-12 * E0 / dt, "a p = 0 boundary does no work");
   TEST_ASSERT(info.fluid_energy + diss <= E0 * (1.0 + 1e-12), "fluid energy + dissipation <= initial");
   TEST_REL(info.fluid_energy + diss, E0, 0.08, "fluid energy + dissipation = initial to 8 % (BE, jump start)");
}

static void TestCflRatio()
{
   std::cout << "\n[undrained stiffness ratio]\n";
   const PoroDamageParams P = Params();
   const real_t G = E_ / 2.0 / (1.0 + nu_);
   const PoroPoint pt = EvalPoroPoint(P, 0.0);
   const real_t lam = P.K - 2.0 * G / 3.0;
   TEST_REL(UndrainedStiffnessRatio(P, G, 2), 1.0 + pt.alpha * pt.alpha / pt.inv_M / (lam + G), 1e-12,
            "2D ratio = 1 + alpha^2 M/(lambda + G) (largest at d = 0 for cmame)");
   TEST_REL(UndrainedStiffnessRatio(P, G, 3), 1.0 + pt.alpha * pt.alpha / pt.inv_M / P.K, 1e-12,
            "3D ratio = 1 + alpha^2 M/K");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   TestCouplingAndStorage();
   TestDiffusion(PoroPreconditioner::Jacobi, "jacobi");
   TestDiffusion(PoroPreconditioner::AMG, "amg");
   TestCflRatio();
   return TestSummary("test_pore_pressure");
}
