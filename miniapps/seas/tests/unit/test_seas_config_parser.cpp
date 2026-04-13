// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// Unit tests for SEASConfig and SEASConfigParser.
// Run: ./seas_test_config_parser
//
// Requires SEAS_USE_TOML for parser tests. Without it, only struct tests run.

#include "mfem.hpp"
#include "../../common/mpi_context.hpp"
#include "../../config/seas_config.hpp"
// Note: TOML parsing tests require a separate binary compiled with
// SEAS_USE_TOML. This test only exercises the SEASConfig struct,
// validation, and CLI overrides (no TOML parsing).
#include "../../config/seas_config_parser.hpp"
#include "../../config/bp5_params.hpp"

#include <cmath>
#include <iostream>
#include <sstream>

using namespace mfem;
using namespace mfem::seas;

// Local test counters (avoid test_macros.hpp static init issues)
static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "FAILED: " << message << "\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_NEAR(value, expected, tol, message) \
   do { \
      num_tests++; \
      real_t _val = (value); \
      real_t _exp = (expected); \
      if (std::abs(_val - _exp) > (tol)) { \
         std::cerr << "FAILED: " << message << " (got " << _val \
                   << ", expected " << _exp << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_PRINT_RESULTS() \
   std::cout << "\n========================================\n" \
             << "  Results: " << num_passed << " passed, " \
             << num_failed << " failed out of " << num_tests << " tests\n" \
             << "========================================\n"

void TestSEASConfigDefaults()
{
   std::cout << "\n=== Test: SEASConfig default values ===\n";
   SEASConfig config;

   TEST_ASSERT(config.mesh.file.empty(), "default mesh file empty");
   TEST_NEAR(config.mesh.scale, 1000.0, 1e-10, "default mesh scale");
   TEST_ASSERT(config.mesh.order == 1, "default mesh order");
   TEST_NEAR(config.material.density, 2670.0, 1e-10, "default density");
   TEST_NEAR(config.material.cs, 3464.0, 1e-10, "default cs");
   TEST_NEAR(config.material.nu, 0.25, 1e-10, "default nu");
   TEST_NEAR(config.friction.V0, 1e-6, 1e-20, "default V0");
   TEST_NEAR(config.friction.f0, 0.6, 1e-10, "default f0");
   TEST_NEAR(config.loading.Vp, 1e-9, 1e-20, "default Vp");
   TEST_ASSERT(config.solver.dg_method == "IP", "default DG method is IP");
   TEST_ASSERT(config.simulation.mode == "qd", "default mode is qd");
   TEST_ASSERT(config.boundary.dirichlet_attrs.empty(), "default dir attrs empty (no preset)");
   TEST_ASSERT(config.boundary.fault_attr == 0, "default fault attr 0 (not set)");
}

void TestSolverConfigDefaults()
{
   std::cout << "\n=== Test: SolverConfig defaults ===\n";
   SolverConfig sc;
   TEST_ASSERT(sc.face_basis_type == BasisType::GaussLobatto,
               "face_basis_type defaults to GaussLobatto");
   TEST_ASSERT(sc.dg_method == "IP", "dg_method defaults to IP");
   TEST_NEAR(sc.blr_tol, 1e-12, 1e-20, "blr_tol defaults to 1e-12");
}

void TestMaterialDerivedProperties()
{
   std::cout << "\n=== Test: MaterialConfig derived properties ===\n";
   MaterialConfig mat;
   mat.density = 2670.0;
   mat.cs = 3464.0;
   mat.nu = 0.25;

   real_t expected_mu = 2670.0 * 3464.0 * 3464.0;
   TEST_NEAR(mat.mu(), expected_mu, 1.0, "mu = rho * cs^2");

   real_t expected_lambda = 2.0 * 0.25 * expected_mu / (1.0 - 2.0 * 0.25);
   TEST_NEAR(mat.lambda(), expected_lambda, 1.0, "lambda = 2*nu*mu/(1-2*nu)");

   // For nu=0.25, lambda should equal mu
   TEST_NEAR(mat.lambda(), mat.mu(), 1.0, "nu=0.25 → lambda == mu");
}

void TestCLIOverrides()
{
   std::cout << "\n=== Test: CLI overrides ===\n";
   SEASConfig config;

   std::vector<std::string> overrides = {
      "material.density=3000",
      "friction.b=0.05",
      "mesh.file=test.msh",
      "solver.dg_method=BR2",
      "time.max_steps=100"
   };

   SEASConfigParser::ApplyCLIOverrides(config, overrides);

   TEST_NEAR(config.material.density, 3000.0, 1e-10, "override density");
   TEST_NEAR(config.friction.b, 0.05, 1e-10, "override friction.b");
   TEST_ASSERT(config.mesh.file == "test.msh", "override mesh file");
   TEST_ASSERT(config.solver.dg_method == "BR2", "override dg_method");
   TEST_ASSERT(config.time.max_steps == 100, "override max_steps");

   // Unchanged fields
   TEST_NEAR(config.material.cs, 3464.0, 1e-10, "cs unchanged");
   TEST_NEAR(config.friction.V0, 1e-6, 1e-20, "V0 unchanged");
}

void TestValidationSuccess()
{
   std::cout << "\n=== Test: Validation passes for valid config ===\n";

   SEASConfig config;
   config.mesh.file = "test.msh";
   config.boundary.dirichlet_attrs = {5};
   config.boundary.natural_attrs = {1};
   config.boundary.fault_attr = 3;
   SEASConfigParser::Validate(config);
   TEST_ASSERT(true, "validation: valid config passes");
}

#ifdef SEAS_USE_TOML
void TestParseValidTOML()
{
   std::cout << "\n=== Test: Parse valid TOML ===\n";

   std::string toml = R"(
benchmark = "bp5"

[mesh]
file = "bp5/mesh/bp5_1000m.msh"
scale = 1000
order = 2

[friction]
b = 0.05

[solver]
dg_method = "IP"
solver_type = "mumps-blr"

[time]
max_steps = 500

[output]
output_dir = "results"
output_prefix = "test_run"
)";

   auto config = SEASConfigParser::ParseString(toml);

   // Preset values from BP5
   BP5Params bp5;
   TEST_NEAR(config.material.density, bp5.rho, 1e-10, "bp5 preset: density");
   TEST_NEAR(config.material.cs, bp5.cs, 1e-10, "bp5 preset: cs");
   TEST_NEAR(config.friction.V0, bp5.V0, 1e-20, "bp5 preset: V0");
   TEST_NEAR(config.friction.f0, bp5.f0, 1e-10, "bp5 preset: f0");
   TEST_NEAR(config.loading.Vp, bp5.Vp, 1e-20, "bp5 preset: Vp");

   // Overridden values
   TEST_ASSERT(config.mesh.file == "bp5/mesh/bp5_1000m.msh", "parsed mesh file");
   TEST_NEAR(config.mesh.scale, 1000.0, 1e-10, "parsed mesh scale");
   TEST_ASSERT(config.mesh.order == 2, "parsed order");
   TEST_NEAR(config.friction.b, 0.05, 1e-10, "override b=0.05 over preset");
   TEST_ASSERT(config.solver.dg_method == "IP", "parsed dg_method");
   TEST_ASSERT(config.solver.solver_type == "mumps-blr", "parsed solver_type");
   TEST_ASSERT(config.time.max_steps == 500, "parsed max_steps");
   TEST_ASSERT(config.output.output_dir == "results", "parsed output_dir");
   TEST_ASSERT(config.output.output_prefix == "test_run", "parsed prefix");
}

void TestParseMinimalTOML()
{
   std::cout << "\n=== Test: Parse minimal TOML (no preset) ===\n";

   std::string toml = R"(
[mesh]
file = "my_mesh.msh"
)";

   auto config = SEASConfigParser::ParseString(toml);
   TEST_ASSERT(config.mesh.file == "my_mesh.msh", "minimal: mesh file");
   TEST_ASSERT(config.benchmark.empty(), "minimal: no benchmark");
   // All other fields should have struct defaults
   TEST_NEAR(config.material.density, 2670.0, 1e-10, "minimal: default density");
}

void TestParseBoundarySection()
{
   std::cout << "\n=== Test: Parse boundary section ===\n";

   std::string toml = R"(
[mesh]
file = "test.msh"

[boundary]
dirichlet = [1, 2, 3, 4]
natural = [5, 6]
fault = 7
)";

   auto config = SEASConfigParser::ParseString(toml);
   TEST_ASSERT(config.boundary.dirichlet_attrs.size() == 4, "4 dirichlet attrs");
   TEST_ASSERT(config.boundary.dirichlet_attrs.count(3) == 1, "has dir attr 3");
   TEST_ASSERT(config.boundary.natural_attrs.size() == 2, "2 natural attrs");
   TEST_ASSERT(config.boundary.fault_attr == 7, "fault attr 7");
}

void TestBP5PresetBitExact()
{
   std::cout << "\n=== Test: BP5 preset matches BP5Params bit-for-bit ===\n";

   std::string toml = R"(
benchmark = "bp5"
[mesh]
file = "dummy.msh"
)";

   auto config = SEASConfigParser::ParseString(toml);
   BP5Params bp5;

   // Bit-level comparison of all doubles
   TEST_ASSERT(config.material.density == bp5.rho, "bit-exact: density");
   TEST_ASSERT(config.material.cs == bp5.cs, "bit-exact: cs");
   TEST_ASSERT(config.material.nu == bp5.nu, "bit-exact: nu");
   TEST_ASSERT(config.friction.V0 == bp5.V0, "bit-exact: V0");
   TEST_ASSERT(config.friction.f0 == bp5.f0, "bit-exact: f0");
   TEST_ASSERT(config.friction.b == bp5.b, "bit-exact: b");
   TEST_ASSERT(config.friction.L0 == bp5.L0, "bit-exact: L0");
   TEST_ASSERT(config.friction.L_nuc == bp5.L_nuc, "bit-exact: L_nuc");
   TEST_ASSERT(config.friction.a0 == bp5.a0, "bit-exact: a0");
   TEST_ASSERT(config.friction.amax == bp5.amax, "bit-exact: amax");
   TEST_ASSERT(config.friction.sigma_n == bp5.sigma_n, "bit-exact: sigma_n");
   TEST_ASSERT(config.loading.Vp == bp5.Vp, "bit-exact: Vp");
   TEST_ASSERT(config.loading.V_init == bp5.V_init, "bit-exact: V_init");
   TEST_ASSERT(config.loading.V_nuc == bp5.V_nuc, "bit-exact: V_nuc");
   TEST_ASSERT(config.loading.delta_tau_factor == bp5.delta_tau_factor, "bit-exact: delta_tau_factor");
   TEST_ASSERT(config.fault_geom.Wf == bp5.Wf, "bit-exact: Wf");
   TEST_ASSERT(config.fault_geom.lf == bp5.lf, "bit-exact: lf");
   TEST_ASSERT(config.fault_geom.hs == bp5.hs, "bit-exact: hs");
   TEST_ASSERT(config.fault_geom.H == bp5.H, "bit-exact: H");
   TEST_ASSERT(config.time.t_final == bp5.t_final, "bit-exact: t_final");
}
#endif // SEAS_USE_TOML

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);

   TestSEASConfigDefaults();
   TestSolverConfigDefaults();
   TestMaterialDerivedProperties();
   TestCLIOverrides();
   TestValidationSuccess();

#ifdef SEAS_USE_TOML
   TestParseValidTOML();
   TestParseMinimalTOML();
   TestParseBoundarySection();
   TestBP5PresetBitExact();
#else
   std::cout << "\n  (TOML parser tests skipped — SEAS_USE_TOML not defined)\n";
#endif

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
