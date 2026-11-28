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

// normal_strain permeability (cmame 3d_hydromech/em_5J/elasticity_mesh2x.i:
// permeability_model = normal_strain, perm_exponent = 10, damage_gradient
// normal, element_size h_c, anisotropic, threshold 0.5, f_c = 1, k0 = 5e-19;
// NDSmallDeformationIsotropicElasticity::updatePermeabilityForCracking):
// k = k0 I + d^b (h_c (1 + eps_nn))^2 / 12 (I - n n). With d linear in x the
// normal is e_x, with u = (e x, 0, 0) eps_nn = e; a pressure gradient along the
// normal sees k0 only, one in the crack plane the full fracture term.
static void TestNormalStrainPermeability()
{
   std::cout << "\n[normal_strain permeability tensor]\n";
   const real_t L = 1e-3;
   Mesh serial = Mesh::MakeCartesian3D(3, 2, 2, Element::TETRAHEDRON, 1.5 * L, L, L);
   ParMesh mesh(MPI_COMM_WORLD, serial);
   H1_FECollection h1(1, 3);
   ParFiniteElementSpace V(&mesh, &h1, 3, Ordering::byNODES), S(&mesh, &h1);
   const IntegrationRule &ir = IntRules.Get(Geometry::TETRAHEDRON, 2);
   PoroDamageParams P = Params();
   P.perm_model = PermeabilityModel::NormalStrain;
   P.k0 = 5e-19;
   P.perm_exponent = 10.0;
   P.d_threshold = 0.5;
   P.fc = 1.0;
   P.anisotropic = true;
   P.normal_source = CrackNormalSource::DamageGradient;
   Array<int> none;
   PorePressureOperator op(V, S, P, none, ir, PoroPreconditioner::Jacobi);

   // Elem::hmax of the tets of a 0.5 mm cube cell: the body diagonal or a face
   // diagonal or... at most the body diagonal, at least the face diagonal.
   const Vector &h = op.ElementHmax();
   real_t hloc[2] = {h.Size() ? h.Max() : 0.0, h.Size() ? -h.Min() : -1e30}, hg[2];
   MPI_Allreduce(hloc, hg, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, MPI_COMM_WORLD);
   TEST_REL(hg[0], 0.5 * L * std::sqrt(3.0), 1e-12, "largest hmax = cell body diagonal");
   TEST_ASSERT(-hg[1] >= 0.5 * L * std::sqrt(2.0) * (1.0 - 1e-12), "hmax >= cell face diagonal");

   const real_t e = 2e-3, d0 = 0.7, dslope = 0.1 / (1.5 * L);
   ParGridFunction dg(&S), ug(&V), pg(&S);
   FunctionCoefficient df([&](const Vector &x) { return d0 + dslope * x(0); });
   VectorFunctionCoefficient uf(3, [&](const Vector &x, Vector &u) { u = 0.0; u(0) = e * x(0); });
   dg.ProjectCoefficient(df);
   ug.ProjectCoefficient(uf);
   Vector D(S.GetTrueVSize()), U(V.GetTrueVSize()), Pp(S.GetTrueVSize());
   dg.GetTrueDofs(D);
   ug.GetTrueDofs(U);
   op.SetState(D, U);

   // reference integrals of d^b h^2 over the mesh at the same points
   real_t loc[2] = {0.0, 0.0}, glob[2];
   for (int el = 0; el < mesh.GetNE(); el++)
   {
      ElementTransformation &T = *mesh.GetElementTransformation(el);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         T.SetIntPoint(&ip);
         const real_t w = ip.weight * T.Weight();
         loc[0] += w;
         loc[1] += w * std::pow(dg.GetValue(T, ip), 10.0) * h(el) * h(el);
      }
   }
   MPI_Allreduce(loc, glob, 2, MPITypeMap<real_t>::mpi_type, MPI_SUM, MPI_COMM_WORLD);
   const real_t vol = glob[0], kfrac = glob[1] * (1.0 + e) * (1.0 + e) / 12.0;

   const real_t gp = 1e9;
   real_t Ef, rate;
   FunctionCoefficient px([&](const Vector &x) { return gp * x(0); });
   pg.ProjectCoefficient(px);
   pg.GetTrueDofs(Pp);
   op.Energies(D, Pp, Ef, rate);
   TEST_REL(rate, P.k0 / P.viscosity * gp * gp * vol, 1e-10, "gradient along the normal: matrix k0 only");
   FunctionCoefficient py([&](const Vector &x) { return gp * x(1); });
   pg.ProjectCoefficient(py);
   pg.GetTrueDofs(Pp);
   op.Energies(D, Pp, Ef, rate);
   TEST_REL(rate, (P.k0 * vol + kfrac) / P.viscosity * gp * gp, 1e-10,
            "gradient in the crack plane: k0 + d^b (h (1 + eps_nn))^2 / 12");
   TEST_ASSERT(kfrac > 1e6 * P.k0 * vol, "fracture term dominates in the plane");

   // below the threshold: matrix permeability in every direction
   dg = 0.4;
   dg.GetTrueDofs(D);
   op.SetState(D, U);
   op.Energies(D, Pp, Ef, rate);
   TEST_REL(rate, P.k0 / P.viscosity * gp * gp * vol, 1e-10, "d < 0.5: k0 I");

   // one step with the tensor in the operator: uniform p, no strain change, no
   // drained boundary -> p stays uniform (the diffusion block has zero row sums)
   dg.ProjectCoefficient(df);
   dg.GetTrueDofs(D);
   Pp = 2.5e6;
   PorePressureOperator::StepInfo info;
   op.Step(D, U, U, 1e-8, Pp, info);
   TEST_REL(Pp.Min(), 2.5e6, 1e-9, "tensor step keeps a uniform pressure (min)");
   TEST_REL(Pp.Max(), 2.5e6, 1e-9, "tensor step keeps a uniform pressure (max)");
}

// The P1-simplex loops against the generic MFEM integrators (same quadrature,
// same formulas): force, one step with a drained boundary, energies; isotropic
// darcy_poiseuille in 2D (cmame 2d_hydromech values) and the normal_strain
// tensor in 3D (cmame 3d_hydromech values).
static void TestFastPath(int dim)
{
   std::cout << "\n[fast path against the generic integrators, " << dim << "D]\n";
   const real_t L = 1e-3;
   Mesh serial = dim == 2 ? Mesh::MakeCartesian2D(7, 5, Element::TRIANGLE, true, 1.4 * L, L)
                 : Mesh::MakeCartesian3D(4, 3, 3, Element::TETRAHEDRON, 1.4 * L, L, L);
   ParMesh mesh(MPI_COMM_WORLD, serial);
   H1_FECollection h1(1, dim);
   ParFiniteElementSpace V(&mesh, &h1, dim, Ordering::byNODES), S(&mesh, &h1);
   const IntegrationRule &ir = IntRules.Get(dim == 2 ? Geometry::TRIANGLE : Geometry::TETRAHEDRON, 2);
   PoroDamageParams P = Params();
   P.k0 = 5e-19;
   P.perm_exponent = 10.0;
   if (dim == 2)
   {
      P.perm_model = PermeabilityModel::DarcyPoiseuille;
      P.wc = 4.13e-7;
   }
   else
   {
      P.perm_model = PermeabilityModel::NormalStrain;
   }
   Array<int> ess_bdr(mesh.bdr_attributes.Max()), drained;
   ess_bdr = 0;
   ess_bdr[0] = 1;
   S.GetEssentialTrueDofs(ess_bdr, drained);
   PorePressureOperator fast(V, S, P, drained, ir, PoroPreconditioner::Jacobi);
   PorePressureOperator gen(V, S, P, drained, ir, PoroPreconditioner::Jacobi);
   TEST_ASSERT(fast.FastPath(), "simplex mesh takes the fast path");
   gen.SetFastPath(false);

   ParGridFunction dg(&S), pg(&S), ug(&V), u0(&V);
   FunctionCoefficient df([&](const Vector &x)
   { return std::min(1.0, std::max(0.0, 0.2 + 0.7 * x(0) / (1.4 * L) + 0.2 * std::sin(9e3 * x(1)))); });
   FunctionCoefficient pf([&](const Vector &x) { return 3e6 * std::cos(4e3 * x(0)) + 2e9 * x(1); });
   VectorFunctionCoefficient uf(dim, [&](const Vector &x, Vector &u)
   { for (int c = 0; c < x.Size(); c++) { u(c) = 1e-3 * L * std::sin(3e3 * x(c) + c) + 2e-3 * x((c + 1) % x.Size()); } });
   VectorFunctionCoefficient u0f(dim, [&](const Vector &x, Vector &u)
   { for (int c = 0; c < x.Size(); c++) { u(c) = 5e-4 * L * std::cos(2e3 * x(c)); } });
   dg.ProjectCoefficient(df);
   pg.ProjectCoefficient(pf);
   ug.ProjectCoefficient(uf);
   u0.ProjectCoefficient(u0f);
   Vector D(S.GetTrueVSize()), Pp(S.GetTrueVSize()), U(V.GetTrueVSize()), U0(V.GetTrueVSize());
   dg.GetTrueDofs(D);
   pg.GetTrueDofs(Pp);
   ug.GetTrueDofs(U);
   u0.GetTrueDofs(U0);

   Vector Ff(V.GetTrueVSize()), Fg(V.GetTrueVSize());
   Ff = 0.0;
   Fg = 0.0;
   fast.AddCouplingForce(D, Pp, Ff);
   gen.AddCouplingForce(D, Pp, Fg);
   Fg -= Ff;
   const real_t fn = std::sqrt(InnerProduct(MPI_COMM_WORLD, Ff, Ff));
   TEST_ASSERT(std::sqrt(InnerProduct(MPI_COMM_WORLD, Fg, Fg)) <= 1e-11 * fn && fn > 0.0,
               "coupling force agrees to 1e-11");

   PorePressureOperator::StepInfo fi, gi;
   Vector Pf(Pp), Pg(Pp);
   for (int t : drained) { Pf(t) = Pg(t) = 1e6; }
   // a large step so that the diffusion block matters next to the storage block
   const real_t dt = dim == 2 ? 1e-6 : 1e-4;
   fast.Step(D, U, U0, dt, Pf, fi);
   gen.Step(D, U, U0, dt, Pg, gi);
   Vector diff(Pf);
   diff -= Pg;
   const real_t pn = std::sqrt(InnerProduct(MPI_COMM_WORLD, Pg, Pg));
   TEST_ASSERT(std::sqrt(InnerProduct(MPI_COMM_WORLD, diff, diff)) <= 1e-8 * pn, "step agrees to 1e-8");
   TEST_REL(fi.fluid_energy, gi.fluid_energy, 1e-8, "fluid energy agrees");
   TEST_REL(fi.dissipation_rate, gi.dissipation_rate, 1e-8, "dissipation rate agrees");
   TEST_REL(fi.boundary_power, gi.boundary_power, 1e-6, "drained boundary power agrees");
   TEST_ASSERT(gi.dissipation_rate > 0.0 && std::abs(gi.boundary_power) > 0.0, "the step has flow");
   // a second step reuses the stored sparsity
   fast.Step(D, U0, U, dt, Pf, fi);
   gen.Step(D, U0, U, dt, Pg, gi);
   diff = Pf;
   diff -= Pg;
   TEST_ASSERT(std::sqrt(InnerProduct(MPI_COMM_WORLD, diff, diff)) <= 1e-8 * pn, "second step agrees to 1e-8");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   TestFastPath(2);
   TestFastPath(3);
   TestNormalStrainPermeability();
   TestCouplingAndStorage();
   TestDiffusion(PoroPreconditioner::Jacobi, "jacobi");
   TestDiffusion(PoroPreconditioner::AMG, "amg");
   TestCflRatio();
   return TestSummary("test_pore_pressure");
}
