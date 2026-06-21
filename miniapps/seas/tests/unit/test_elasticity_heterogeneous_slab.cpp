// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_elasticity_heterogeneous_slab.cpp — Phase 2b Stage 3 of
// document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md
//
// Acceptance test: "a 2-layer depth profile yields the analytically expected
// traction on a slab test."  Builds the heterogeneous (Mode::Coefficient,
// depth-profile) IP ElasticityDomainOperator, imposes a UNIFORM-SHEAR
// displacement u_x = gamma * y (constant strain eps_xy = gamma/2 everywhere),
// recovers the fault traction, and asserts it equals the ANALYTIC value at
// each fault DOF:
//
//   For pure shear, sigma = C : eps has only sigma_xy = 2*mu*eps_xy = mu*gamma
//   (tr(eps)=0 so the lambda term vanishes).  The fault normal is y, so the
//   fault traction is sigma.n = (sigma_xy, 0, 0) and its STRIKE (x) component
//   is mu(depth) * gamma.  Because mu varies with depth (2 layers), the strike
//   traction VARIES with depth — a constant-mu operator would give one value
//   everywhere and fail this test.  This is the heterogeneity proof that
//   Stage 1b's constant-external bit-for-bit test cannot provide.
//
// The expected mu(depth) is taken from the SAME depth-profile evaluator the
// operator's coefficient uses (wrapper->eval_at_xyz), so a per-qp evaluation
// bug (wrong/constant mu) is caught by the per-DOF identity.  The layer
// boundary is aligned with an element-band boundary; fault DOFs exactly on the
// boundary (mixed-layer L2 projection) are skipped.
//
// Serial (Mesh) fixture, no MPI — mirrors test_elasticity_operator.cpp.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "fault_mesh_fixture.hpp"   // test::CreateTestMesh3D (R-503)

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_fail = 0;
int g_checks = 0;

void Check(bool ok, const std::string &msg)
{
   g_checks++;
   if (!ok) { g_fail++; std::cerr << "FAILED: " << msg << "\n"; }
}

using test::CreateTestMesh3D;   // shared fault-mesh fixture (R-503)
}  // namespace

int main()
{
   // Mesh: depth (-z) spans [0, 12]; element-band depth boundaries at 0,3,6,9,12.
   const real_t Lx = 10.0, Ly = 10.0, Lz = 12.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 4, Lx, Ly, Lz);

   // 2-layer depth profile.  Boundary at depth 6 (an element-band boundary):
   //   layer 0 (0–6 km): mu_0 = rho_0 * vs_0^2 = 2600 * 3000^2 = 2.340e10 Pa
   //   layer 1 (6+  km): mu_1 = rho_1 * vs_1^2 = 3000 * 4000^2 = 4.800e10 Pa
   // (depths in metres for the profile; the mesh's Lz units double as metres.)
   const real_t depth_boundary = 6.0;
   std::vector<DepthProfileLayer> layers(2);
   layers[0] = {0.0,            depth_boundary, 6000.0, 3000.0, 2600.0, "constant"};
   layers[1] = {depth_boundary, 1000.0,         8000.0, 4000.0, 3000.0, "constant"};
   const real_t mu0_expect = 2600.0 * 3000.0 * 3000.0;
   const real_t mu1_expect = 3000.0 * 4000.0 * 4000.0;

   std::unique_ptr<DepthProfile1DMaterial> wrapper =
      MakeDepthProfile1DMaterial(layers, 'z');
   MaterialField material = wrapper->field;
   Check(material.mode == MaterialField::Mode::Coefficient,
         "depth-profile MaterialField is Mode::Coefficient");
   Check(material.lambda_coef != nullptr && material.mu_coef != nullptr,
         "depth-profile lambda_coef/mu_coef non-null");

   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   const real_t Vp = 1.0e-9;
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   DomainConfig dcfg;
   ElasticityDomainOperator<Mesh> op(
      mesh, /*order=*/1, *material.lambda_coef, *material.mu_coef,
      Vp, /*Wf=*/Lz, /*lf=*/2.0 * Lx, bc,
      DGMethod::IP, SolverType::CG_AMG, dcfg);

   const int nfd = op.GetNumFaultDOFs();
   Check(nfd > 0, "operator has fault DOFs");

   // Impose pure shear u_x = gamma * y (constant eps_xy); P1-exact projection.
   const real_t gamma = 1.0e-4;
   GridFunction u(&op.GetFESpace());
   VectorFunctionCoefficient ushear(3, [gamma](const Vector &x, Vector &out)
   {
      out = 0.0;
      out(0) = gamma * x(1);   // u_x = gamma * y
   });
   u.ProjectCoefficient(ushear);

   Vector slip_bc(2 * nfd);
   slip_bc = 0.0;   // zero slip => recovered traction is the pure stress part

   Vector trac;
   op.ComputeTraction(u, slip_bc, trac);
   Check(trac.Size() == 2 * nfd, "traction size == 2*nfd");

   Vector x2, x3;   // x3 = depth (= -Z), per GetFaultCoords2D
   op.GetFaultCoords2D(x2, x3);
   Check(x3.Size() == nfd, "fault-coord size == nfd");

   const real_t rel_tol = 1.0e-4;   // patch test is exact for constant strain
   int n_layer0 = 0, n_layer1 = 0;
   real_t sum0 = 0.0, sum1 = 0.0;
   real_t max_dip_ratio = 0.0;

   if (nfd > 0 && trac.Size() == 2 * nfd && x3.Size() == nfd)
   {
      for (int i = 0; i < nfd; i++)
      {
         const real_t depth = x3(i);
         // Skip the layer-boundary DOF row: its L2 projection mixes qps from
         // both layers, so no single-layer analytic value applies.
         if (std::abs(depth - depth_boundary) < 0.5) { continue; }

         // Expected mu from the SAME depth-profile evaluator the operator uses.
         real_t lam_e = 0.0, mu_e = 0.0, rho_e = 0.0;
         wrapper->eval_at_xyz(0.0, 0.0, -depth, lam_e, mu_e, rho_e);
         const real_t exp_strike = mu_e * gamma;

         const real_t dip    = trac(2 * i);
         const real_t strike = trac(2 * i + 1);

         Check(std::abs(std::abs(strike) - exp_strike) <= rel_tol * exp_strike,
               "strike traction == mu(depth)*gamma (analytic patch test)");
         max_dip_ratio = std::max(max_dip_ratio, std::abs(dip) / exp_strike);

         if (depth < depth_boundary) { n_layer0++; sum0 += std::abs(strike); }
         else                        { n_layer1++; sum1 += std::abs(strike); }
      }

      // Pure shear => no dip (normal-plane) traction.
      Check(max_dip_ratio < rel_tol,
            "dip traction negligible vs strike (pure shear)");

      // Heterogeneity actually exercised: DOFs in BOTH layers, with the
      // analytic modulus ratio between them.
      Check(n_layer0 > 0, "fault DOFs sampled in layer 0 (shallow)");
      Check(n_layer1 > 0, "fault DOFs sampled in layer 1 (deep)");
      if (n_layer0 > 0 && n_layer1 > 0)
      {
         const real_t mean0 = sum0 / n_layer0;
         const real_t mean1 = sum1 / n_layer1;
         const real_t ratio_obs = mean1 / mean0;
         const real_t ratio_exp = mu1_expect / mu0_expect;
         Check(std::abs(ratio_obs - ratio_exp) <= 1.0e-3 * ratio_exp,
               "deep/shallow strike-traction ratio == mu_1/mu_0 (heterogeneity)");
         std::cout << "  layer0: " << n_layer0 << " DOFs, mean|strike|=" << mean0
                   << " (expect mu_0*gamma=" << mu0_expect * gamma << ")\n"
                   << "  layer1: " << n_layer1 << " DOFs, mean|strike|=" << mean1
                   << " (expect mu_1*gamma=" << mu1_expect * gamma << ")\n"
                   << "  ratio obs=" << ratio_obs << " expect=" << ratio_exp << "\n";
      }
   }

   std::cout << "=== test_elasticity_heterogeneous_slab ===\n"
             << "  fault DOFs: " << nfd << "\n"
             << "  checks: " << g_checks << ", failures: " << g_fail << "\n"
             << (g_fail == 0 ? "PASS\n" : "FAIL\n");
   return (g_fail == 0) ? 0 : 1;
}
