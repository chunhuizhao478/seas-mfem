// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// Unit tests for seas_config_bridge.hpp (Phase 7a).
// Tests BuildBP5Params round-trip, enum parsing, BoundaryConfig/DomainConfig.
// Run: ./seas_test_config_bridge

#include "mfem.hpp"
#include "../../config/seas_config_bridge.hpp"
#include "../../config/seas_config_parser.hpp"
#include "test_macros.hpp"

#include <cmath>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

void TestBuildBP5ParamsRoundTrip()
{
   std::cout << "\n=== Test: BuildBP5Params round-trip (BP5 preset) ===\n";

   SEASConfig config;
   // Apply BP5 preset (same as benchmark="bp5")
   config.benchmark = "bp5";
   // Manually apply defaults since we're not parsing TOML
   BP5Params bp5_ref;  // Reference: BP5Params default constructor
   config.material.density = bp5_ref.rho;
   config.material.cs = bp5_ref.cs;
   config.material.nu = bp5_ref.nu;
   config.friction.V0 = bp5_ref.V0;
   config.friction.f0 = bp5_ref.f0;
   config.friction.b = bp5_ref.b;
   config.friction.L0 = bp5_ref.L0;
   config.friction.L_nuc = bp5_ref.L_nuc;
   config.friction.a0 = bp5_ref.a0;
   config.friction.amax = bp5_ref.amax;
   config.friction.sigma_n = bp5_ref.sigma_n;
   config.loading.Vp = bp5_ref.Vp;
   config.loading.V_init = bp5_ref.V_init;
   config.loading.V_nuc = bp5_ref.V_nuc;
   config.loading.delta_tau_factor = bp5_ref.delta_tau_factor;
   config.loading.nucleation_eps = bp5_ref.nucleation_eps;
   config.loading.smooth_nucleation = bp5_ref.smooth_nucleation;
   config.fault_geom.Wf = bp5_ref.Wf;
   config.fault_geom.lf = bp5_ref.lf;
   config.fault_geom.hs = bp5_ref.hs;
   config.fault_geom.ht = bp5_ref.ht;
   config.fault_geom.H = bp5_ref.H;
   config.fault_geom.l_vw = bp5_ref.l_vw;
   config.fault_geom.w_nuc = bp5_ref.w_nuc;
   config.time.t_final = bp5_ref.t_final;

   BP5Params p = BuildBP5Params(config);

   // Verify all mapped fields
   TEST_ASSERT(p.rho == bp5_ref.rho, "round-trip: rho");
   TEST_ASSERT(p.cs == bp5_ref.cs, "round-trip: cs");
   TEST_ASSERT(p.nu == bp5_ref.nu, "round-trip: nu");
   TEST_ASSERT(p.V0 == bp5_ref.V0, "round-trip: V0");
   TEST_ASSERT(p.f0 == bp5_ref.f0, "round-trip: f0");
   TEST_ASSERT(p.b == bp5_ref.b, "round-trip: b");
   TEST_ASSERT(p.L0 == bp5_ref.L0, "round-trip: L0");
   TEST_ASSERT(p.L_nuc == bp5_ref.L_nuc, "round-trip: L_nuc");
   TEST_ASSERT(p.a0 == bp5_ref.a0, "round-trip: a0");
   TEST_ASSERT(p.amax == bp5_ref.amax, "round-trip: amax");
   TEST_ASSERT(p.sigma_n == bp5_ref.sigma_n, "round-trip: sigma_n");
   TEST_ASSERT(p.Vp == bp5_ref.Vp, "round-trip: Vp");
   TEST_ASSERT(p.V_init == bp5_ref.V_init, "round-trip: V_init");
   TEST_ASSERT(p.V_nuc == bp5_ref.V_nuc, "round-trip: V_nuc");
   TEST_ASSERT(p.delta_tau_factor == bp5_ref.delta_tau_factor, "round-trip: delta_tau_factor");
   TEST_ASSERT(p.nucleation_eps == bp5_ref.nucleation_eps, "round-trip: nucleation_eps");
   TEST_ASSERT(p.smooth_nucleation == bp5_ref.smooth_nucleation, "round-trip: smooth_nucleation");
   TEST_ASSERT(p.Wf == bp5_ref.Wf, "round-trip: Wf");
   TEST_ASSERT(p.lf == bp5_ref.lf, "round-trip: lf");
   TEST_ASSERT(p.hs == bp5_ref.hs, "round-trip: hs");
   TEST_ASSERT(p.ht == bp5_ref.ht, "round-trip: ht");
   TEST_ASSERT(p.H == bp5_ref.H, "round-trip: H");
   TEST_ASSERT(p.l_vw == bp5_ref.l_vw, "round-trip: l_vw");
   TEST_ASSERT(p.w_nuc == bp5_ref.w_nuc, "round-trip: w_nuc");
   TEST_ASSERT(p.t_final == bp5_ref.t_final, "round-trip: t_final");

   // Verify unmapped fields retain defaults
   TEST_ASSERT(p.V_zero == bp5_ref.V_zero, "unmapped: V_zero");
}

void TestBuildBP5ParamsNonPreset()
{
   std::cout << "\n=== Test: BuildBP5Params with non-preset values ===\n";

   SEASConfig config;
   config.material.density = 3000.0;
   config.friction.b = 0.05;
   config.loading.Vp = 2e-9;

   BP5Params p = BuildBP5Params(config);
   TEST_NEAR(p.rho, 3000.0, 1e-10, "override: density");
   TEST_NEAR(p.b, 0.05, 1e-10, "override: b");
   TEST_NEAR(p.Vp, 2e-9, 1e-20, "override: Vp");

   // Unmapped fields keep BP5Params defaults
   BP5Params def;
   TEST_ASSERT(p.V_zero == def.V_zero, "default: V_zero unchanged");
}

void TestParseDGMethod()
{
   std::cout << "\n=== Test: ParseDGMethod ===\n";
   TEST_ASSERT(ParseDGMethod("IP") == DGMethod::IP, "IP → DGMethod::IP");
   TEST_ASSERT(ParseDGMethod("ip") == DGMethod::IP, "ip → DGMethod::IP");
   TEST_ASSERT(ParseDGMethod("BR2") == DGMethod::BR2, "BR2 → DGMethod::BR2");
   TEST_ASSERT(ParseDGMethod("br2") == DGMethod::BR2, "br2 → DGMethod::BR2");
}

void TestParseSolverType()
{
   std::cout << "\n=== Test: ParseSolverType ===\n";
   TEST_ASSERT(ParseSolverType("cg") == SolverType::CG_AMG, "cg → CG_AMG");
   TEST_ASSERT(ParseSolverType("mumps") == SolverType::MUMPS, "mumps → MUMPS");
   TEST_ASSERT(ParseSolverType("mumps-blr") == SolverType::MUMPS_BLR, "mumps-blr → MUMPS_BLR");
   TEST_ASSERT(ParseSolverType("superlu") == SolverType::SUPERLU, "superlu → SUPERLU");
   TEST_ASSERT(ParseSolverType("strumpack") == SolverType::STRUMPACK, "strumpack → STRUMPACK");
   TEST_ASSERT(ParseSolverType("gmres") == SolverType::GMRES_AMG, "gmres → GMRES_AMG");
}

void TestBuildBoundaryConfig()
{
   std::cout << "\n=== Test: BuildBoundaryConfig ===\n";

   BoundaryTomlConfig bdr;
   bdr.dirichlet_attrs = {5};
   bdr.natural_attrs = {1};
   bdr.fault_attr = 3;
   real_t Vp = 1e-9;

   BoundaryConfig bc = BuildBoundaryConfig(bdr, Vp);

   TEST_ASSERT(bc.dirichlet_attrs.count(5) == 1, "dir attr 5");
   TEST_ASSERT(bc.natural_attrs.count(1) == 1, "nat attr 1");
   TEST_ASSERT(bc.fault_attr == 3, "fault attr 3");
   TEST_ASSERT(bc.default_dirichlet_func != nullptr, "DirichletFunc non-null");

   // Verify DirichletFunc output at y=2000, t=1e9
   Vector x(3); x(0) = 0; x(1) = 2000; x(2) = -10000;
   Vector u(3);
   bc.default_dirichlet_func(x, 1e9, u);
   real_t expected = Vp * 1e9 * 0.5;  // y > 1000 → Vp*t/2
   TEST_NEAR(u(0), expected, 1e-20, "DirichletFunc at y=2000");
   TEST_NEAR(u(1), 0.0, 1e-30, "DirichletFunc u_y=0");
   TEST_NEAR(u(2), 0.0, 1e-30, "DirichletFunc u_z=0");
}

void TestBuildDomainConfig()
{
   std::cout << "\n=== Test: BuildDomainConfig ===\n";

   SolverConfig solver;
   solver.face_basis_type = 4;  // ClosedUniform
   solver.penalty_factor = 2.5;
   solver.blr_tol = 1e-14;
   solver.check_residual = true;

   DomainConfig dc = BuildDomainConfig(solver);

   TEST_ASSERT(dc.face_basis_type == 4, "face_basis_type mapped");
   TEST_NEAR(dc.penalty_factor, 2.5, 1e-10, "penalty_factor mapped");
   TEST_NEAR(dc.blr_tol, 1e-14, 1e-20, "blr_tol mapped");
   TEST_ASSERT(dc.check_residual == true, "check_residual mapped");
}

int main(int argc, char *argv[])
{
   TestBuildBP5ParamsRoundTrip();
   TestBuildBP5ParamsNonPreset();
   TestParseDGMethod();
   TestParseSolverType();
   TestBuildBoundaryConfig();
   TestBuildDomainConfig();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
