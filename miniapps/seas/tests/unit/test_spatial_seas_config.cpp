// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_seas_config.cpp — Phase 1 of
// document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md
//
// Covers the QD (quasi-dynamic) config-schema additions:
//   * the [solver] table (SolverSpec) — every field round-trips,
//   * the QD [time] knobs (rk45_atol/rk45_rtol/dt_init/dt_max_years/
//     plate_rate_vp/use_petsc_ts) — every field round-trips,
//   * defaults when [solver] / QD [time] keys are absent,
//   * unknown [solver].type is a hard error at parse (Edge Case),
//   * ParseQDSolverType maps every accepted string to the right
//     mfem::seas::SolverType,
//   * a dynamic-only [numerics] key does NOT abort the parser (the
//     "ignored" WARNING is the QD driver's job — plan §C / R-001 — not the
//     shared parser's; here we only confirm the parse succeeds).
//
// As in test_spatial_friction_config.cpp: MFEM is built WITHOUT
// MFEM_USE_EXCEPTIONS here, so MFEM_ABORT calls std::abort().  Abort-expected
// cases run in a forked child; success-path cases run inline.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"
#include "../../domain/elasticity_operator.hpp"  // mfem::seas::SolverType

#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas::spatial;
using mfem::seas::SolverType;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

namespace
{

// Run a parse inside a forked child.  Returns true iff the child aborted
// (non-zero exit, signal, or core dump).  Child stderr is silenced.
bool ParseAbortsInChild(const std::string& toml_text)
{
   ::fflush(stdout);
   ::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0)
   {
      std::cerr << "fork() failed: " << ::strerror(errno) << "\n";
      return false;
   }
   if (pid == 0)
   {
      ::freopen("/dev/null", "w", stderr);
      try
      {
         (void)ParseSpatialFrictionConfigString(toml_text);
      }
      catch (...)
      {
         ::_exit(1);
      }
      ::_exit(0);
   }
   int status = 0;
   while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) { /* retry */ }
   if (WIFEXITED(status))   { return WEXITSTATUS(status) != 0; }
   if (WIFSIGNALED(status)) { return true; }
   return true;
}

// A minimal valid SLW config.  `solver_block` and `extra_time` are spliced
// verbatim into the [solver] table region and the [time] table respectively
// (pass "" to omit).  All other blocks are the required minimum the shared
// parser demands ([meta]/[material_constant_fallback]/[pore_pressure]/[mesh]/
// [velocity]/[stress]/[numerics]/[time]/[output] + a friction block).
std::string QDConfig(const std::string& solver_block,
                     const std::string& extra_time)
{
   std::ostringstream oss;
   oss << "[meta]\n"
       << "schema_version = 1\n"
       << "law = \"slip_weakening\"\n"
       << "[material_constant_fallback]\n"
       << "lambda=32.0e9\nmu=32.0e9\nrho=2670.0\n"
       << "[pore_pressure]\nP_p_pa=0.0\nP_p_grad_pa_per_m=0.0\nmin_sigma_n_pa=0.0\n"
       << "[mesh]\npath=\"/dev/null\"\norder=1\n"
       << "[velocity]\nmodel=\"cvmh\"\ndataset_root=\"/tmp/x\"\noverride_path=\"\"\n"
       << "[stress]\nkind=\"constant_tensor\"\n"
       << "sigma_xx_pa=0.0\nsigma_yy_pa=0.0\nsigma_zz_pa=0.0\n"
       << "sigma_xy_pa=0.0\nsigma_yz_pa=0.0\nsigma_xz_pa=0.0\n"
       << "[numerics]\nader_order=2\nmixed_flux=\"none\"\ncfl=0.5\nuse_pml=false\n"
       << "[time]\ntfinal=\"12s\"\nt_initial=0.0\ndt_initial=\"auto\"\ndt_max=\"0.1s\"\n"
       << extra_time
       << solver_block
       << "[output]\noutput_dir=\"out\"\n"
       << "[friction.slip_weakening]\n"
       << "mu_s_default = 1.1\nmu_d_default = 0.5\nd_c_default = 0.5\n"
       << "cohesion_default = 0.0\n";
   return oss.str();
}

}  // namespace

// T-1  full [solver] + QD [time] block: every field round-trips.
static void T_1_solver_and_qd_time_fields()
{
   std::cout << "\n[T-1] full [solver] + QD [time] parse\n";
   const std::string solver =
      "[solver]\n"
      "type=\"mumps\"\n"
      "ksp_rtol=1.0e-9\n"
      "ksp_atol=1.0e-12\n"
      "ksp_maxit=250\n"
      "amg_elasticity_options=false\n"
      "amg_print_level=2\n"
      "blr_tol=1.0e-8\n"
      "residual_check=false\n";
   const std::string qd_time =
      "rk45_atol=2.5e-8\n"
      "rk45_rtol=1.0e-40\n"
      "dt_init=0.25\n"
      "dt_max_years=0.05\n"
      "plate_rate_vp=1.0e-9\n"
      "use_petsc_ts=true\n";
   const auto cfg = ParseSpatialFrictionConfigString(QDConfig(solver, qd_time));

   TEST_ASSERT(cfg.solver.type == "mumps", "solver.type round-trip");
   TEST_ASSERT(cfg.solver.ksp_rtol == 1.0e-9, "solver.ksp_rtol");
   TEST_ASSERT(cfg.solver.ksp_atol == 1.0e-12, "solver.ksp_atol");
   TEST_ASSERT(cfg.solver.ksp_maxit == 250, "solver.ksp_maxit");
   TEST_ASSERT(cfg.solver.amg_elasticity_options == false,
               "solver.amg_elasticity_options");
   TEST_ASSERT(cfg.solver.amg_print_level == 2, "solver.amg_print_level");
   TEST_ASSERT(cfg.solver.blr_tol == 1.0e-8, "solver.blr_tol");
   TEST_ASSERT(cfg.solver.residual_check == false, "solver.residual_check");

   TEST_ASSERT(cfg.time.rk45_atol == 2.5e-8, "time.rk45_atol");
   TEST_ASSERT(cfg.time.rk45_rtol == 1.0e-40, "time.rk45_rtol");
   TEST_ASSERT(cfg.time.dt_init == 0.25, "time.dt_init");
   TEST_ASSERT(cfg.time.dt_max_years == 0.05, "time.dt_max_years");
   TEST_ASSERT(cfg.time.plate_rate_vp == 1.0e-9, "time.plate_rate_vp");
   TEST_ASSERT(cfg.time.use_petsc_ts == true, "time.use_petsc_ts");
}

// T-2  defaults when [solver] block and QD [time] keys are ABSENT.
static void T_2_defaults_when_absent()
{
   std::cout << "\n[T-2] [solver] absent + QD [time] keys absent ⇒ defaults\n";
   const auto cfg = ParseSpatialFrictionConfigString(QDConfig("", ""));

   TEST_ASSERT(cfg.solver.type == "cg_amg", "default solver.type=cg_amg");
   TEST_ASSERT(cfg.solver.ksp_rtol == 1e-8, "default ksp_rtol");
   TEST_ASSERT(cfg.solver.ksp_atol == 0.0, "default ksp_atol");
   TEST_ASSERT(cfg.solver.ksp_maxit == 2000, "default ksp_maxit");
   TEST_ASSERT(cfg.solver.amg_elasticity_options == true,
               "default amg_elasticity_options=true");
   TEST_ASSERT(cfg.solver.amg_relax_type == 8, "default amg_relax_type=8 (l1-sym-GS)");
   TEST_ASSERT(cfg.solver.amg_aggressive_levels == 0,
               "default amg_aggressive_levels=0 (off; stagnates DG elasticity)");
   TEST_ASSERT(cfg.solver.amg_print_level == 0, "default amg_print_level");
   TEST_ASSERT(cfg.solver.blr_tol == 1e-10, "default blr_tol");
   TEST_ASSERT(cfg.solver.residual_check == true,
               "default residual_check=true (R-006)");

   TEST_ASSERT(cfg.time.rk45_atol == 1e-7, "default rk45_atol");
   TEST_ASSERT(cfg.time.rk45_rtol == 1e-50, "default rk45_rtol");
   TEST_ASSERT(cfg.time.dt_init == -1.0, "default dt_init=-1 (derive)");
   TEST_ASSERT(cfg.time.dt_max_years == 0.1, "default dt_max_years=0.1");
   TEST_ASSERT(cfg.time.plate_rate_vp == -1.0, "default plate_rate_vp=-1");
   TEST_ASSERT(cfg.time.use_petsc_ts == false, "default use_petsc_ts=false");
}

// T-3  unknown [solver].type is a hard error at parse (Edge Case).
static void T_3_unknown_solver_type_aborts()
{
   std::cout << "\n[T-3] unknown [solver].type aborts at parse\n";
   const std::string bad = QDConfig("[solver]\ntype=\"bogus_solver\"\n", "");
   TEST_ASSERT(ParseAbortsInChild(bad),
               "unknown solver.type must abort (no silent default)");
}

// T-4  ParseQDSolverType maps every accepted string to the right enum.
static void T_4_parse_qd_solver_type_mapping()
{
   std::cout << "\n[T-4] ParseQDSolverType string→enum mapping\n";
   TEST_ASSERT(ParseQDSolverType("cg_amg")    == SolverType::CG_AMG,    "cg_amg");
   TEST_ASSERT(ParseQDSolverType("gmres_amg") == SolverType::GMRES_AMG, "gmres_amg");
   TEST_ASSERT(ParseQDSolverType("gmres_ilu") == SolverType::GMRES_BlockILU, "gmres_ilu");
   TEST_ASSERT(ParseQDSolverType("mumps")     == SolverType::MUMPS,     "mumps");
   TEST_ASSERT(ParseQDSolverType("mumps_blr") == SolverType::MUMPS_BLR, "mumps_blr");
   TEST_ASSERT(ParseQDSolverType("superlu")   == SolverType::SUPERLU,   "superlu");
   TEST_ASSERT(ParseQDSolverType("strumpack") == SolverType::STRUMPACK, "strumpack");
}

// T-5  a dynamic-only [numerics] key does NOT abort the parser (the QD
//      driver, not the shared parser, owns the "ignored" warning — §C/R-001).
static void T_5_dynamic_key_does_not_abort_parser()
{
   std::cout << "\n[T-5] dynamic-only [numerics] key parses (no throw)\n";
   // The base config's [numerics] already carries the dynamic-only keys
   // ader_order / mixed_flux / use_pml.  The shared parser must accept them
   // (it is the QD DRIVER, not the parser, that warns they are ignored —
   // §C/R-001).  If the parser aborted on any of them the test binary would
   // die before the assertion below.
   const std::string toml = QDConfig("", "");
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.schema_version == 1,
               "config with dynamic-only [numerics] keys parses successfully");
   TEST_ASSERT(cfg.numerics.ader_order == 2,
               "dynamic key value is still parsed into NumericsSpec");
}

// T-6  (R-101) integer-typed QD [time] knobs parse via toml_double's int
//      branch.  The double QD fields are read with toml_double (not the
//      review's as_floating()-only suggestion, which aborts on TOML
//      integers); this guards that an integer literal like `dt_init = 3`
//      becomes 3.0 rather than aborting.
static void T_6_integer_qd_time_knob_parses()
{
   std::cout << "\n[T-6] integer-valued QD [time] knob parses (toml_double int branch)\n";
   const std::string qd_time =
      "dt_init=3\n"        // integer literal
      "dt_max_years=1\n";  // integer literal
   const auto cfg = ParseSpatialFrictionConfigString(QDConfig("", qd_time));
   TEST_ASSERT(cfg.time.dt_init == 3.0, "integer dt_init -> 3.0");
   TEST_ASSERT(cfg.time.dt_max_years == 1.0, "integer dt_max_years -> 1.0");
}

int main()
{
   std::cout << "=== test_spatial_seas_config (Phase 1 QD schema) ===\n";
   T_1_solver_and_qd_time_fields();
   T_2_defaults_when_absent();
   T_3_unknown_solver_type_aborts();
   T_4_parse_qd_solver_type_mapping();
   T_5_dynamic_key_does_not_abort_parser();
   T_6_integer_qd_time_knob_parses();

   std::cout << "\n=== SUMMARY: " << num_passed << "/" << num_tests
             << " passed, " << num_failed << " failed ===\n";
   return (num_failed == 0) ? 0 : 1;
}
