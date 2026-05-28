// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_friction_config.cpp — Phase 1 of
// spatial_dynamic_rupture_plan.md (rev-3).
//
// Plan §Phase 1 §Acceptance: 12 tests covering every validation rule
// (schema_version mismatch, missing law, both laws present, unknown
// keys, the D-3 alias `d_o`, barrier-sentinel guard R-114, stress-mode
// conflict, time parse, material fallback bounds, mesh.path presence).
//
// MFEM is built WITHOUT MFEM_USE_EXCEPTIONS in this environment, so
// MFEM_ABORT calls std::abort().  We isolate abort-expected tests in a
// child process via fork() + waitpid() so the parent test driver
// survives.  Success-path tests run inline.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas::spatial;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

namespace
{

// Run a parse inside a forked child.  Returns true iff the child
// process aborted (any non-zero exit, signal, or core dump).  The
// child's stderr is silenced to keep test output clean.
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
      // Child: silence stderr, attempt the parse, exit 0 on success.
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

std::string MinimalLSWHeader(int schema_version = 1,
                             const std::string& law = "slip_weakening")
{
   std::ostringstream oss;
   oss << "[meta]\n";
   oss << "schema_version = " << schema_version << "\n";
   if (!law.empty())
   {
      oss << "law = \"" << law << "\"\n";
   }
   oss << "[material_constant_fallback]\n"
       << "lambda=32.0e9\nmu=32.0e9\nrho=2670.0\n"
       << "[pore_pressure]\nP_p_pa=0.0\nP_p_grad_pa_per_m=0.0\nmin_sigma_n_pa=0.0\n"
       << "[mesh]\npath=\"/dev/null\"\norder=1\n"
       << "[velocity]\nmodel=\"cvmh\"\ndataset_root=\"/tmp/x\"\noverride_path=\"\"\n"
       << "[stress]\nkind=\"constant_tensor\"\n"
       << "sigma_xx_pa=0.0\nsigma_yy_pa=0.0\nsigma_zz_pa=0.0\n"
       << "sigma_xy_pa=0.0\nsigma_yz_pa=0.0\nsigma_xz_pa=0.0\n"
       << "[numerics]\nader_order=2\nmixed_flux=\"none\"\ncfl=0.5\nuse_pml=false\n"
       << "[time]\ntfinal=\"12s\"\nt_initial=0.0\ndt_initial=\"auto\"\ndt_max=\"0.1s\"\n"
       << "[output]\noutput_dir=\"out\"\nrestart_prefix=\"cp\"\n"
       << "paraview_volume=\"hdf5\"\nparaview_bulk=\"hdf5\"\nparaview_fault=\"hdf5\"\n"
       << "paraview_volume_dt=\"0.05s\"\nparaview_bulk_dt=\"0.05s\"\nparaview_fault_dt=\"0.001s\"\n"
       << "paraview_volume_zfp_tol=1e-3\nparaview_bulk_zfp_tol=1e-3\nparaview_fault_zfp_tol=1e-12\n"
       << "max_snapshots=5000\ncheckpoint_every_steps=10000\n";
   return oss.str();
}

std::string MinimalLSWBlock(real_t mu_s = 1.1, real_t mu_d = 0.5,
                            real_t d_c = 0.5,
                            const std::string& d_key = "d_c_default")
{
   std::ostringstream oss;
   oss << "[friction.slip_weakening]\n"
       << "mu_s_default = " << mu_s << "\n"
       << "mu_d_default = " << mu_d << "\n"
       << d_key << " = "    << d_c  << "\n"
       << "cohesion_default = 0.0\n";
   return oss.str();
}

std::string MinimalRSBlock()
{
   return std::string(
      "[friction.rate_state]\n"
      "f_0_default = 0.6\n"
      "V_0_default = 1.0e-6\n"
      "eta = \"auto\"\n"
      "a_default = 0.010\n"
      "b_default = 0.015\n"
      "Dc_default = 0.004\n"
      "V_init_default = 1.0e-9\n"
      "sigma_n_default = 50.0e6\n");
}

}  // namespace

// T-1  minimal valid LSW config parses (geoffrey2010 values survive D-3)
static void T_1_minimal_lsw_parses_geoffrey2010()
{
   std::cout << "\n[T-1] minimal LSW config parses with mu_s=1.1\n";
   const std::string toml = MinimalLSWHeader() + MinimalLSWBlock(1.1, 0.5, 0.5);
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.schema_version == 1, "schema_version round-trip");
   TEST_ASSERT(cfg.law == FrictionLawKind::SlipWeakening, "law round-trip");
   TEST_ASSERT(cfg.slip_weakening.has_value(), "slip_weakening present");
   TEST_ASSERT(!cfg.rate_state.has_value(), "rate_state absent");
   TEST_ASSERT(cfg.slip_weakening->mu_s_default == 1.1, "mu_s_default = 1.1");
}

// T-2  schema_version mismatch aborts
static void T_2_schema_version_mismatch_aborts()
{
   std::cout << "\n[T-2] schema_version != 1 aborts\n";
   const std::string toml = MinimalLSWHeader(/*schema_version=*/2)
                            + MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml), "schema_version=2 must abort");
}

// T-3  missing law aborts
static void T_3_missing_law_aborts()
{
   std::cout << "\n[T-3] missing [meta].law aborts\n";
   const std::string toml = MinimalLSWHeader(1, /*law=*/"")
                            + MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml), "missing law must abort");
}

// T-4  both friction blocks present aborts
static void T_4_both_friction_blocks_aborts()
{
   std::cout << "\n[T-4] both [friction.slip_weakening] AND "
                "[friction.rate_state] aborts\n";
   const std::string toml = MinimalLSWHeader()
                            + MinimalLSWBlock()
                            + MinimalRSBlock();
   TEST_ASSERT(ParseAbortsInChild(toml), "both LSW and RS blocks must abort");
}

// T-5  D-3 alias d_o_default works
static void T_5_d_o_alias_works()
{
   std::cout << "\n[T-5] D-3 d_o_default alias resolves to d_c_default\n";
   const std::string toml = MinimalLSWHeader()
                            + MinimalLSWBlock(1.1, 0.5, 0.5, "d_o_default");
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.slip_weakening->d_c_default == 0.5,
               "d_c_default resolved from d_o_default alias");
}

// T-6  D-3 alias collision (both d_c and d_o) aborts
static void T_6_d_o_collision_aborts()
{
   std::cout << "\n[T-6] supplying both d_c_default AND d_o_default aborts\n";
   std::string toml = MinimalLSWHeader();
   toml +=
      "[friction.slip_weakening]\n"
      "mu_s_default = 1.1\n"
      "mu_d_default = 0.5\n"
      "d_c_default = 0.5\n"
      "d_o_default = 0.4\n"
      "cohesion_default = 0.0\n";
   TEST_ASSERT(ParseAbortsInChild(toml),
               "both d_c_default and d_o_default must abort");
}

// T-7  R-114 barrier-sentinel guard on default
static void T_7_r114_default_guard()
{
   std::cout << "\n[T-7] mu_s_default > 1.0e5 aborts (R-114)\n";
   const std::string toml = MinimalLSWHeader()
                            + MinimalLSWBlock(1.0e6, 0.5, 0.5);
   TEST_ASSERT(ParseAbortsInChild(toml),
               "mu_s_default = 1e6 must abort");
}

// T-8  R-114 barrier-sentinel guard on spatial rule
static void T_8_r114_rule_guard()
{
   std::cout << "\n[T-8] [[spatial]] rule with mu_s > 1.0e5 aborts (R-114)\n";
   std::string toml = MinimalLSWHeader() + MinimalLSWBlock();
   toml +=
      "[[friction.slip_weakening.spatial]]\n"
      "kind=\"box\"\n"
      "x_min_m=-1e6\nx_max_m=1e6\n"
      "y_min_m=-1e6\ny_max_m=1e6\n"
      "z_min_m=-1e6\nz_max_m=1e6\n"
      "mu_s = 1e6\nmu_d = 0.5\nd_c = 0.5\n";
   TEST_ASSERT(ParseAbortsInChild(toml),
               "spatial rule with mu_s = 1e6 must abort");
}

// T-9  stress mode conflict: kind=constant_tensor with sidecar_path
static void T_9_stress_mode_conflict_const_with_path()
{
   std::cout << "\n[T-9] [stress] kind=constant_tensor + sidecar_path aborts\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=32e9
rho=2670
[pore_pressure]
P_p_pa=0
P_p_grad_pa_per_m=0
min_sigma_n_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
override_path=""
[stress]
kind = "constant_tensor"
sigma_xx_pa=0
sigma_yy_pa=0
sigma_zz_pa=0
sigma_xy_pa=0
sigma_yz_pa=0
sigma_xz_pa=0
sidecar_path = "/tmp/conflict.h5"
[numerics]
ader_order=2
mixed_flux="none"
cfl=0.5
use_pml=false
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   TEST_ASSERT(ParseAbortsInChild(toml),
               "constant_tensor + sidecar_path must abort");
}

// T-10 stress mode conflict: kind=sidecar_hdf5 with sigma_*_pa
static void T_10_stress_mode_conflict_sidecar_with_sigma()
{
   std::cout << "\n[T-10] [stress] kind=sidecar_hdf5 + sigma_xx_pa aborts\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=32e9
rho=2670
[pore_pressure]
P_p_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
[stress]
kind = "sidecar_hdf5"
sidecar_path = "/tmp/s.h5"
sigma_xx_pa = 1e6
[numerics]
ader_order=2
mixed_flux="none"
cfl=0.5
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   TEST_ASSERT(ParseAbortsInChild(toml),
               "sidecar_hdf5 + sigma_xx_pa must abort");
}

// T-11 time parser round-trips
static void T_11_time_parser()
{
   std::cout << "\n[T-11] SpatialTimeParseSeconds round-trips\n";
   TEST_ASSERT(SpatialTimeParseSeconds("0.05s") == 0.05, "0.05s");
   TEST_ASSERT(SpatialTimeParseSeconds("12s") == 12.0, "12s");
   TEST_ASSERT(SpatialTimeParseSeconds("1.5e-2s") == 1.5e-2, "1.5e-2s");
   TEST_ASSERT(SpatialTimeParseSeconds("3.14") == 3.14, "no suffix");
   TEST_ASSERT(SpatialTimeParseSeconds("auto") == -1.0, "auto sentinel");
}

// T-13  R-004 regression: missing top-level [stress] block aborts.
static void T_13_missing_stress_block_aborts()
{
   std::cout << "\n[T-13] missing [stress] block aborts (R-004)\n";
   // Build header WITHOUT the [stress] block by hand.
   std::ostringstream oss;
   oss << "[meta]\nschema_version = 1\nlaw = \"slip_weakening\"\n"
       << "[material_constant_fallback]\nlambda=32e9\nmu=32e9\nrho=2670\n"
       << "[pore_pressure]\nP_p_pa=0\n"
       << "[mesh]\npath=\"/dev/null\"\norder=1\n"
       << "[velocity]\nmodel=\"cvmh\"\ndataset_root=\"/tmp/x\"\noverride_path=\"\"\n"
       // (no [stress] block)
       << "[numerics]\nader_order=2\nmixed_flux=\"none\"\ncfl=0.5\n"
       << "[time]\ntfinal=\"12s\"\n"
       << "[output]\noutput_dir=\"out\"\n"
       << "[friction.slip_weakening]\n"
       << "mu_s_default=1.1\nmu_d_default=0.5\nd_c_default=0.5\n"
       << "cohesion_default=0\n";
   TEST_ASSERT(ParseAbortsInChild(oss.str()),
               "missing [stress] block must abort");
}

// T-14  R-004 regression: missing top-level [time] block aborts.
static void T_14_missing_time_block_aborts()
{
   std::cout << "\n[T-14] missing [time] block aborts (R-004)\n";
   std::ostringstream oss;
   oss << "[meta]\nschema_version = 1\nlaw = \"slip_weakening\"\n"
       << "[material_constant_fallback]\nlambda=32e9\nmu=32e9\nrho=2670\n"
       << "[pore_pressure]\nP_p_pa=0\n"
       << "[mesh]\npath=\"/dev/null\"\norder=1\n"
       << "[velocity]\nmodel=\"cvmh\"\ndataset_root=\"/tmp/x\"\noverride_path=\"\"\n"
       << "[stress]\nkind=\"constant_tensor\"\n"
       << "sigma_xx_pa=0\nsigma_yy_pa=0\nsigma_zz_pa=0\n"
       << "sigma_xy_pa=0\nsigma_yz_pa=0\nsigma_xz_pa=0\n"
       << "[numerics]\nader_order=2\nmixed_flux=\"none\"\ncfl=0.5\n"
       // (no [time] block)
       << "[output]\noutput_dir=\"out\"\n"
       << "[friction.slip_weakening]\n"
       << "mu_s_default=1.1\nmu_d_default=0.5\nd_c_default=0.5\n"
       << "cohesion_default=0\n";
   TEST_ASSERT(ParseAbortsInChild(oss.str()),
               "missing [time] block must abort");
}

// T-15  R-006 regression: dt_initial = 0 aborts (must be > 0 or "auto").
static void T_15_dt_initial_zero_aborts()
{
   std::cout << "\n[T-15] [time].dt_initial = 0 aborts (R-006)\n";
   std::string toml = MinimalLSWHeader();
   const std::string before = "dt_initial=\"auto\"";
   const std::string after  = "dt_initial=0.0";
   const auto pos = toml.find(before);
   MFEM_VERIFY(pos != std::string::npos,
               "MinimalLSWHeader changed shape — re-grep dt_initial");
   toml.replace(pos, before.size(), after);
   toml += MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml),
               "dt_initial = 0 must abort (must be > 0 or \"auto\")");
}

// T-16  R-007 regression: use_pml as a quoted string aborts with a
// clean toml_bool message (not a toml11 internal abort).
static void T_16_use_pml_string_aborts()
{
   std::cout << "\n[T-16] [numerics].use_pml = \"false\" (quoted) aborts (R-007)\n";
   std::string toml = MinimalLSWHeader();
   const std::string before = "use_pml=false";
   const std::string after  = "use_pml=\"false\"";
   const auto pos = toml.find(before);
   MFEM_VERIFY(pos != std::string::npos,
               "MinimalLSWHeader changed shape — re-grep use_pml");
   toml.replace(pos, before.size(), after);
   toml += MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml),
               "use_pml=\"false\" (quoted) must abort");
}

// T-17  R-201 regression: missing [pore_pressure] block aborts.
static void T_17_missing_pore_pressure_block_aborts()
{
   std::cout << "\n[T-17] missing [pore_pressure] block aborts (R-201)\n";
   std::ostringstream oss;
   oss << "[meta]\nschema_version = 1\nlaw = \"slip_weakening\"\n"
       << "[material_constant_fallback]\nlambda=32e9\nmu=32e9\nrho=2670\n"
       // (no [pore_pressure] block)
       << "[mesh]\npath=\"/dev/null\"\norder=1\n"
       << "[velocity]\nmodel=\"cvmh\"\ndataset_root=\"/tmp/x\"\noverride_path=\"\"\n"
       << "[stress]\nkind=\"constant_tensor\"\n"
       << "sigma_xx_pa=0\nsigma_yy_pa=0\nsigma_zz_pa=0\n"
       << "sigma_xy_pa=0\nsigma_yz_pa=0\nsigma_xz_pa=0\n"
       << "[numerics]\nader_order=2\nmixed_flux=\"none\"\ncfl=0.5\n"
       << "[time]\ntfinal=\"12s\"\n"
       << "[output]\noutput_dir=\"out\"\n"
       << "[friction.slip_weakening]\n"
       << "mu_s_default=1.1\nmu_d_default=0.5\nd_c_default=0.5\n"
       << "cohesion_default=0\n";
   TEST_ASSERT(ParseAbortsInChild(oss.str()),
               "missing [pore_pressure] block must abort");
}

// T-18  R-202 regression: literal numeric dt_initial = -1.0 aborts
// (must not be silently treated as the "auto" sentinel).
static void T_18_dt_initial_literal_negative_aborts()
{
   std::cout << "\n[T-18] dt_initial = -1.0 (literal) aborts (R-202)\n";
   std::string toml = MinimalLSWHeader();
   const std::string before = "dt_initial=\"auto\"";
   const std::string after  = "dt_initial=-1.0";
   const auto pos = toml.find(before);
   MFEM_VERIFY(pos != std::string::npos,
               "MinimalLSWHeader changed shape — re-grep dt_initial");
   toml.replace(pos, before.size(), after);
   toml += MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml),
               "literal dt_initial = -1.0 must abort (not silently \"auto\")");
}

// T-19  R-203 regression: paraview_fault_dt = 0 aborts.
static void T_19_paraview_fault_dt_zero_aborts()
{
   std::cout << "\n[T-19] paraview_fault_dt = 0 aborts (R-203)\n";
   std::string toml = MinimalLSWHeader();
   const std::string before = "paraview_fault_dt=\"0.001s\"";
   const std::string after  = "paraview_fault_dt=0.0";
   const auto pos = toml.find(before);
   MFEM_VERIFY(pos != std::string::npos,
               "MinimalLSWHeader changed shape — re-grep paraview_fault_dt");
   toml.replace(pos, before.size(), after);
   toml += MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml),
               "paraview_fault_dt = 0 must abort");
}

// T-20  R-205 regression: [stress] block without `kind` key aborts.
static void T_20_stress_kind_missing_aborts()
{
   std::cout << "\n[T-20] [stress] without 'kind' key aborts (R-205)\n";
   std::ostringstream oss;
   oss << "[meta]\nschema_version = 1\nlaw = \"slip_weakening\"\n"
       << "[material_constant_fallback]\nlambda=32e9\nmu=32e9\nrho=2670\n"
       << "[pore_pressure]\nP_p_pa=0\n"
       << "[mesh]\npath=\"/dev/null\"\norder=1\n"
       << "[velocity]\nmodel=\"cvmh\"\ndataset_root=\"/tmp/x\"\noverride_path=\"\"\n"
       // [stress] present but missing `kind` key
       << "[stress]\nsigma_xx_pa=0\nsigma_yy_pa=0\nsigma_zz_pa=0\n"
       << "sigma_xy_pa=0\nsigma_yz_pa=0\nsigma_xz_pa=0\n"
       << "[numerics]\nader_order=2\nmixed_flux=\"none\"\ncfl=0.5\n"
       << "[time]\ntfinal=\"12s\"\n"
       << "[output]\noutput_dir=\"out\"\n"
       << "[friction.slip_weakening]\n"
       << "mu_s_default=1.1\nmu_d_default=0.5\nd_c_default=0.5\n"
       << "cohesion_default=0\n";
   TEST_ASSERT(ParseAbortsInChild(oss.str()),
               "[stress] without 'kind' must abort");
}

// T-21  Phase N: the single nucleation kind is "gradual_overstress";
//       absent-kind, "strength_reduction", and "overstress" all abort
//       cleanly (validator rejects anything but the one supported kind).
static void T_21_nucleation_kind_must_be_gradual_overstress()
{
   std::cout << "\n[T-21] [nucleation] kind must be gradual_overstress "
                "(Phase N)\n";
   // Absent `kind` (the round-7 default of "strength_reduction") aborts.
   {
      std::string toml = MinimalLSWHeader() + MinimalLSWBlock();
      toml += "[nucleation]\n"
              "[nucleation.gradual_overstress]\n"
              "center_x_m=0\ncenter_y_m=0\ncenter_z_m=0\n"
              "radius_dip_m=3000\nradius_strike_m=3000\n"
              "delta_tau_dip_pa=0\ndelta_tau_strike_pa=25e6\n"
              "T_nuc_s=\"1.0s\"\n";
      TEST_ASSERT(ParseAbortsInChild(toml),
                  "[nucleation] without `kind` aborts (no implicit default)");
   }
   // Old "strength_reduction" alias is no longer accepted.
   {
      std::string toml = MinimalLSWHeader() + MinimalLSWBlock();
      toml += "[nucleation]\n"
              "kind=\"strength_reduction\"\n";
      TEST_ASSERT(ParseAbortsInChild(toml),
                  "kind=\"strength_reduction\" aborts (Phase N removal)");
   }
}

// T-22  Phase N: a minimal valid `[nucleation]` block with
//       kind=\"gradual_overstress\" parses round-trip into the new
//       GradualOverstressSpec fields.
static void T_22_nucleation_kind_gradual_overstress_parses()
{
   std::cout << "\n[T-22] [nucleation].kind = gradual_overstress parses "
                "(Phase N)\n";
   std::string toml = MinimalLSWHeader() + MinimalLSWBlock();
   toml += "[nucleation]\n"
           "kind=\"gradual_overstress\"\n"
           "[nucleation.gradual_overstress]\n"
           "center_x_m=10.0\ncenter_y_m=20.0\ncenter_z_m=-30.0\n"
           "radius_dip_m=1500.0\nradius_strike_m=2500.0\n"
           "delta_tau_dip_pa=1.0e6\ndelta_tau_strike_pa=2.0e6\n"
           "T_nuc_s=\"1.0s\"\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.nucleation.enabled,
               "[nucleation] block sets enabled = true");
   TEST_ASSERT(cfg.nucleation.kind == NucleationKind::GradualOverstress,
               "kind round-trips to GradualOverstress");
   const auto& g = cfg.nucleation.gradual_overstress;
   TEST_ASSERT(g.center_x_m           ==   10.0, "center_x_m parses");
   TEST_ASSERT(g.center_y_m           ==   20.0, "center_y_m parses");
   TEST_ASSERT(g.center_z_m           ==  -30.0, "center_z_m parses");
   TEST_ASSERT(g.radius_dip_m         == 1500.0, "radius_dip_m parses");
   TEST_ASSERT(g.radius_strike_m      == 2500.0, "radius_strike_m parses");
   TEST_ASSERT(g.delta_tau_dip_pa     ==  1.0e6, "delta_tau_dip_pa parses");
   TEST_ASSERT(g.delta_tau_strike_pa  ==  2.0e6, "delta_tau_strike_pa parses");
   TEST_ASSERT(g.T_nuc_s              ==    1.0, "T_nuc_s parses");
}

// T-23  Phase N: typo'd kind aborts cleanly.
static void T_23_nucleation_kind_typo_aborts()
{
   std::cout << "\n[T-23] [nucleation].kind typo aborts (Phase N)\n";
   std::string toml = MinimalLSWHeader() + MinimalLSWBlock();
   toml += "[nucleation]\n"
           "kind=\"gradual_overstres\"\n"   // intentional typo
           "[nucleation.gradual_overstress]\n"
           "radius_dip_m=3000\nradius_strike_m=3000\nT_nuc_s=\"1.0s\"\n";
   TEST_ASSERT(ParseAbortsInChild(toml),
               "[nucleation].kind = typo must abort");
}

// T-24  R-001 PARSER-LAYER coverage — verifies the Parity-Phase-2
//       default flip (paraview_volume / paraview_bulk default to
//       "off"; paraview_fault stays "hdf5"; paraview_enabled defaults
//       to false).  Does NOT exercise the driver's master-gate logic
//       at drivers/spatial_dyn_driver.cpp:1232-1237 — that lives in
//       the C++ driver and is covered end-to-end by the smoke sbatch
//       jobs/safs/spatial_dyn_smoke_8N_400r_dev_2hr_safs.sbatch
//       (which asserts no PV files appear when paraview_enabled
//       stays false).  R-004 round-3: closing this coverage gap with
//       a free-function extraction is documented as a follow-up.
static void T_24_paraview_default_flip()
{
   std::cout << "\n[T-24] R-001: Parity Phase 2 default flip (volume/bulk = "
                "\"off\", fault stays \"hdf5\")\n";
   // Build a TOML with NO paraview_* keys (deliberately stripping the
   // explicit `paraview_volume = "hdf5"` from MinimalLSWHeader).
   std::ostringstream oss;
   oss << "[meta]\nschema_version = 1\nlaw = \"slip_weakening\"\n"
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
       << "[output]\noutput_dir=\"out\"\nrestart_prefix=\"cp\"\n"
       // NO paraview_volume, paraview_bulk, paraview_fault keys —
       // exercise the parser defaults.
       << "[friction.slip_weakening]\n"
       << "mu_s_default=1.1\nmu_d_default=0.5\nd_c_default=0.5\n"
       << "cohesion_default=0\n";
   const auto cfg = ParseSpatialFrictionConfigString(oss.str());
   TEST_ASSERT(cfg.output.paraview_volume == "off",
               "Parity-Phase-2 default paraview_volume = \"off\"");
   TEST_ASSERT(cfg.output.paraview_bulk == "off",
               "Parity-Phase-2 default paraview_bulk = \"off\"");
   TEST_ASSERT(cfg.output.paraview_fault == "hdf5",
               "paraview_fault default stays \"hdf5\" (documented)");
   TEST_ASSERT(cfg.output.paraview_enabled == false,
               "paraview_enabled default = false (master gate off)");
}

// T-25  PARSER-LAYER coverage of the Parity-Phase-2 OutputSpec
//       extension: the 9 new fields (paraview_every_steps,
//       paraview_fault_legacy_ascii, the three deflate levels, the
//       three regime-adaptive cadences) parse to their documented
//       defaults when omitted from the TOML.  As with T-24, this does
//       NOT exercise the driver's master-gate logic.
static void T_25_paraview_extended_defaults()
{
   std::cout << "\n[T-25] R-001: extended OutputSpec defaults round-trip\n";
   // Same minimal TOML as T-24.
   std::ostringstream oss;
   oss << "[meta]\nschema_version = 1\nlaw = \"slip_weakening\"\n"
       << "[material_constant_fallback]\nlambda=32e9\nmu=32e9\nrho=2670\n"
       << "[pore_pressure]\nP_p_pa=0\nP_p_grad_pa_per_m=0\nmin_sigma_n_pa=0\n"
       << "[mesh]\npath=\"/dev/null\"\norder=1\n"
       << "[velocity]\nmodel=\"cvmh\"\ndataset_root=\"/tmp/x\"\noverride_path=\"\"\n"
       << "[stress]\nkind=\"constant_tensor\"\n"
       << "sigma_xx_pa=0\nsigma_yy_pa=0\nsigma_zz_pa=0\n"
       << "sigma_xy_pa=0\nsigma_yz_pa=0\nsigma_xz_pa=0\n"
       << "[numerics]\nader_order=2\nmixed_flux=\"none\"\ncfl=0.5\nuse_pml=false\n"
       << "[time]\ntfinal=\"12s\"\nt_initial=0\ndt_initial=\"auto\"\ndt_max=\"0.1s\"\n"
       << "[output]\noutput_dir=\"out\"\nrestart_prefix=\"cp\"\n"
       << "[friction.slip_weakening]\n"
       << "mu_s_default=1.1\nmu_d_default=0.5\nd_c_default=0.5\ncohesion_default=0\n";
   const auto cfg = ParseSpatialFrictionConfigString(oss.str());
   TEST_ASSERT(cfg.output.paraview_every_steps == 0,
               "paraview_every_steps default 0");
   TEST_ASSERT(cfg.output.paraview_fault_legacy_ascii == false,
               "paraview_fault_legacy_ascii default false");
   TEST_ASSERT(cfg.output.paraview_volume_deflate_level == -1,
               "paraview_volume_deflate_level default -1");
   TEST_ASSERT(cfg.output.paraview_bulk_deflate_level   == -1,
               "paraview_bulk_deflate_level default -1");
   TEST_ASSERT(cfg.output.paraview_fault_deflate_level  == -1,
               "paraview_fault_deflate_level default -1");
   TEST_ASSERT(cfg.output.paraview_coseismic_dt    < 0.0,
               "paraview_coseismic_dt default unset (< 0)");
   TEST_ASSERT(cfg.output.paraview_nucleation_dt   < 0.0,
               "paraview_nucleation_dt default unset (< 0)");
   TEST_ASSERT(cfg.output.paraview_interseismic_dt < 0.0,
               "paraview_interseismic_dt default unset (< 0)");
}

// T-26 [velocity].use_sidecar = false relaxes the dataset_root /
// override_path non-empty check (production constant-material posture).
static void T_26_velocity_use_sidecar_false_relaxes_dataset_root()
{
   std::cout << "\n[T-26] [velocity].use_sidecar = false parses with "
             "empty dataset_root / override_path\n";
   std::ostringstream oss;
   oss << "[meta]\nschema_version=1\nlaw=\"slip_weakening\"\n"
       << "[material_constant_fallback]\nlambda=32e9\nmu=32e9\nrho=2670\n"
       << "[pore_pressure]\nP_p_pa=0\nP_p_grad_pa_per_m=0\nmin_sigma_n_pa=0\n"
       << "[mesh]\npath=\"/dev/null\"\norder=1\n"
       << "[velocity]\nuse_sidecar=false\nmodel=\"cvmh\"\n"
       << "dataset_root=\"\"\noverride_path=\"\"\n"
       << "[stress]\nkind=\"constant_tensor\"\n"
       << "sigma_xx_pa=0\nsigma_yy_pa=0\nsigma_zz_pa=0\n"
       << "sigma_xy_pa=0\nsigma_yz_pa=0\nsigma_xz_pa=0\n"
       << "[numerics]\nader_order=2\nmixed_flux=\"none\"\ncfl=0.5\nuse_pml=false\n"
       << "[time]\ntfinal=\"12s\"\nt_initial=0\ndt_initial=\"auto\"\ndt_max=\"0.1s\"\n"
       << "[output]\noutput_dir=\"out\"\nrestart_prefix=\"cp\"\n"
       << "[friction.slip_weakening]\n"
       << "mu_s_default=1.1\nmu_d_default=0.5\nd_c_default=0.5\ncohesion_default=0\n";
   const auto cfg = ParseSpatialFrictionConfigString(oss.str());
   TEST_ASSERT(cfg.velocity.use_sidecar == false,
               "use_sidecar=false round-trip");
   TEST_ASSERT(cfg.velocity.dataset_root.empty(),
               "dataset_root empty when use_sidecar=false");
}

// T-27 default use_sidecar=true with empty dataset_root + override_path
// must still abort (back-compat with pre-toggle behaviour).
static void T_27_velocity_use_sidecar_true_default_still_requires_root()
{
   std::cout << "\n[T-27] [velocity] (default use_sidecar=true) with empty "
             "dataset_root / override_path aborts\n";
   std::ostringstream oss;
   oss << "[meta]\nschema_version=1\nlaw=\"slip_weakening\"\n"
       << "[material_constant_fallback]\nlambda=32e9\nmu=32e9\nrho=2670\n"
       << "[pore_pressure]\nP_p_pa=0\nP_p_grad_pa_per_m=0\nmin_sigma_n_pa=0\n"
       << "[mesh]\npath=\"/dev/null\"\norder=1\n"
       // use_sidecar omitted ⇒ default true; both paths empty ⇒ abort.
       << "[velocity]\nmodel=\"cvmh\"\n"
       << "dataset_root=\"\"\noverride_path=\"\"\n"
       << "[stress]\nkind=\"constant_tensor\"\n"
       << "sigma_xx_pa=0\nsigma_yy_pa=0\nsigma_zz_pa=0\n"
       << "sigma_xy_pa=0\nsigma_yz_pa=0\nsigma_xz_pa=0\n"
       << "[numerics]\nader_order=2\nmixed_flux=\"none\"\ncfl=0.5\nuse_pml=false\n"
       << "[time]\ntfinal=\"12s\"\nt_initial=0\ndt_initial=\"auto\"\ndt_max=\"0.1s\"\n"
       << "[output]\noutput_dir=\"out\"\nrestart_prefix=\"cp\"\n"
       << "[friction.slip_weakening]\n"
       << "mu_s_default=1.1\nmu_d_default=0.5\nd_c_default=0.5\ncohesion_default=0\n";
   TEST_ASSERT(ParseAbortsInChild(oss.str()),
               "default use_sidecar=true with empty paths must abort");
}

// T-28  σ_n strength floor: `[friction].sigma_n_strength_floor_pa = 10.0e6`
//        parses round-trip into SpatialFrictionConfig (sliver-blowup plan
//        2026-05-26 §Phase 1).  The `[friction]` table is declared BEFORE
//        the law sub-block.
static void T_28_sigma_n_strength_floor_parses()
{
   std::cout << "\n[T-28] [friction].sigma_n_strength_floor_pa = 10.0e6 "
                "parses\n";
   const std::string toml = MinimalLSWHeader()
                            + "[friction]\n"
                              "sigma_n_strength_floor_pa = 10.0e6\n"
                            + MinimalLSWBlock();
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.sigma_n_strength_floor_pa == 10.0e6,
               "sigma_n_strength_floor_pa round-trips to 10e6");
}

// T-29  σ_n strength floor: absent key ⇒ disabled sentinel -1.0.
static void T_29_sigma_n_strength_floor_absent_disabled()
{
   std::cout << "\n[T-29] absent sigma_n_strength_floor_pa ⇒ -1.0 "
                "(disabled)\n";
   const std::string toml = MinimalLSWHeader() + MinimalLSWBlock();
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.sigma_n_strength_floor_pa == -1.0,
               "absent key defaults to -1.0 (disabled)");
}

// T-30  σ_n strength floor: an EXPLICIT negative value aborts (a negative
//        floor is the disabled sentinel and must be expressed by omitting
//        the key, not by setting it negative).
static void T_30_sigma_n_strength_floor_negative_aborts()
{
   std::cout << "\n[T-30] explicit negative sigma_n_strength_floor_pa "
                "aborts\n";
   const std::string toml = MinimalLSWHeader()
                            + "[friction]\n"
                              "sigma_n_strength_floor_pa = -5.0\n"
                            + MinimalLSWBlock();
   TEST_ASSERT(ParseAbortsInChild(toml),
               "explicit negative sigma_n_strength_floor_pa must abort");
}

// T-31  σ_n strength floor: 0.0 is a VALID value (≡ LSW max(σ_n,0); for RS
//        switches |σ_n|→max(σ_n,0)).  Must parse, not abort.
static void T_31_sigma_n_strength_floor_zero_allowed()
{
   std::cout << "\n[T-31] sigma_n_strength_floor_pa = 0.0 is allowed\n";
   const std::string toml = MinimalLSWHeader()
                            + "[friction]\n"
                              "sigma_n_strength_floor_pa = 0.0\n"
                            + MinimalLSWBlock();
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.sigma_n_strength_floor_pa == 0.0,
               "sigma_n_strength_floor_pa = 0.0 parses (not aborts)");
}

// T-32  R-001: σ_n strength floor key mis-nested under a law sub-block
//        (`[friction.slip_weakening]`) must ABORT, not silently parse to the
//        disabled sentinel.  The floor is a TOP-LEVEL `[friction]` key.
static void T_32_sigma_n_floor_misnested_aborts()
{
   std::cout << "\n[T-32] mis-nested sigma_n_strength_floor_pa "
                "(under [friction.slip_weakening]) aborts (R-001)\n";
   std::string toml = MinimalLSWHeader();
   toml += "[friction.slip_weakening]\n"
           "mu_s_default = 1.1\n"
           "mu_d_default = 0.5\n"
           "d_c_default = 0.5\n"
           "cohesion_default = 0.0\n"
           "sigma_n_strength_floor_pa = 10.0e6\n";   // WRONG: belongs under [friction]
   TEST_ASSERT(ParseAbortsInChild(toml),
               "floor key under [friction.slip_weakening] must abort (R-001)");
}

// T-12 material fallback bounds
static void T_12_material_fallback_bounds()
{
   std::cout << "\n[T-12] material fallback validator (mu <= 0 aborts)\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=-1.0
rho=2670
[pore_pressure]
P_p_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
[stress]
kind = "constant_tensor"
sigma_xx_pa=0
sigma_yy_pa=0
sigma_zz_pa=0
sigma_xy_pa=0
sigma_yz_pa=0
sigma_xz_pa=0
[numerics]
ader_order=2
mixed_flux="none"
cfl=0.5
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   TEST_ASSERT(ParseAbortsInChild(toml), "mu <= 0 must abort");
}

// T-20 (Phase 6 / D3.2): kind="fault_local_prestress" parses; the
// right-lateral / compression-POSITIVE fields land in StressSpec unchanged
// (seeded directly — no Cauchy projection; see T-66 in
// test_compute_safs_params for the per-DOF sign verification).
static void T_20_fault_local_prestress_parses()
{
   std::cout << "\n[FLP-1] [stress] kind=fault_local_prestress parses\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=32e9
rho=2670
[pore_pressure]
P_p_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
[stress]
kind = "fault_local_prestress"
tau_strike_pa = 75.0e6
tau_dip_pa = 0.0
sigma_n_pa = 120.0e6
[numerics]
ader_order=2
mixed_flux="none"
cfl=0.5
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.stress.kind == StressSourceKind::FaultLocalPrestress,
               "kind parses to FaultLocalPrestress");
   TEST_ASSERT(cfg.stress.tau_strike_pa == 75.0e6,
               "tau_strike_pa == +75 MPa (right-lateral positive)");
   TEST_ASSERT(cfg.stress.tau_dip_pa == 0.0, "tau_dip_pa == 0 (pure strike-slip)");
   TEST_ASSERT(cfg.stress.sigma_n_pa == 120.0e6,
               "sigma_n_pa == 120 MPa (compression positive)");
}

// T-21 (Phase 6 / D3.2): fault_local_prestress + a Cauchy sigma key aborts
// (mutual exclusion with kind=constant_tensor).
static void T_21_fault_local_with_sigma_aborts()
{
   std::cout << "\n[FLP-2] [stress] kind=fault_local_prestress + sigma_xx_pa aborts\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=32e9
rho=2670
[pore_pressure]
P_p_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
[stress]
kind = "fault_local_prestress"
sigma_n_pa = 120.0e6
sigma_xx_pa = 1e6
[numerics]
ader_order=2
mixed_flux="none"
cfl=0.5
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   TEST_ASSERT(ParseAbortsInChild(toml),
               "fault_local_prestress + Cauchy sigma_xx_pa must abort");
}

// FLP-3 (Phase 6 / D3.2 R-002): [[stress.patch]] array parses into
// fault_local_patches; FaultLocalPatch::inside() + last-match-wins behave.
static void T_22_fault_local_patches_parse()
{
   std::cout << "\n[FLP-3] [[stress.patch]] parses + inside()/last-match-wins\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=32e9
rho=2670
[pore_pressure]
P_p_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
[stress]
kind = "fault_local_prestress"
tau_strike_pa = 70.0e6
tau_dip_pa = 0.0
sigma_n_pa = 120.0e6
[[stress.patch]]
center_x_m = 0.0
half_x_m = 1500.0
tau_strike_pa = 81.6e6
[[stress.patch]]
center_x_m = 0.0
half_x_m = 3000.0
tau_strike_pa = 78.0e6
[numerics]
ader_order=2
mixed_flux="none"
cfl=0.5
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   const auto &patches = cfg.stress.fault_local_patches;
   TEST_ASSERT(patches.size() == 2, "two [[stress.patch]] parsed");
   TEST_ASSERT(patches.size() == 2 && patches[0].tau_strike_pa == 81.6e6,
               "patch[0] tau_strike = 81.6 MPa");
   TEST_ASSERT(patches.size() == 2 && patches[0].center_x_m == 0.0
               && patches[0].half_x_m == 1500.0,
               "patch[0] box x in [-1500,1500]");
   TEST_ASSERT(patches.size() == 2 && patches[1].tau_strike_pa == 78.0e6,
               "patch[1] tau_strike = 78 MPa");
   // inside(): unconstrained y/z (default NaN center / +inf half) always match.
   TEST_ASSERT(patches.size() == 2 && patches[0].inside(0.0, 9999.0, -9999.0),
               "patch[0] inside at x=0 (y/z unconstrained)");
   TEST_ASSERT(patches.size() == 2 && !patches[0].inside(5000.0, 0.0, 0.0),
               "patch[0] NOT inside at x=5000 (> half_x)");
   // last-match-wins: x=0 is inside BOTH patches; the later patch wins.
   real_t tau = cfg.stress.tau_strike_pa, out = 0.0;
   bool hit = false;
   for (const auto &p : patches)
   { if (p.inside(0.0, 0.0, 0.0)) { out = p.tau_strike_pa; hit = true; } }
   if (hit) { tau = out; }
   TEST_ASSERT(hit && tau == 78.0e6,
               "last-match-wins: x=0 in both patches -> later (78 MPa) wins");
}

// NUM-1 (Phase 6 req 3): [numerics] cfl_safety / fault_iterator / interior_flux
// parse into the enums.  interior_flux="matrix" + mixed_flux="none" is allowed
// (the matrix+mixed_flux mutual-exclusion only fires when mixed_flux != none).
static void T_23_numerics_selectors_parse()
{
   std::cout << "\n[NUM-1] [numerics] cfl_safety/fault_iterator/interior_flux parse\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=32e9
rho=2670
[pore_pressure]
P_p_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
[stress]
kind = "constant_tensor"
sigma_xx_pa=0
sigma_yy_pa=0
sigma_zz_pa=0
sigma_xy_pa=0
sigma_yz_pa=0
sigma_xz_pa=0
[numerics]
ader_order=2
mixed_flux="none"
cfl=0.5
cfl_safety="dg"
fault_iterator="substep"
interior_flux="matrix"
[material]
kind="depth_profile_1d"
profile_csv="/tmp/profile.csv"
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   // NOTE: interior_flux="matrix" now requires a non-Constant [material]
   // (Phase 6 req-3 guard, completed with req 1); the depth_profile_1d block
   // above satisfies it.
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.numerics.cfl_safety == CflSafety::Dg, "cfl_safety=dg");
   TEST_ASSERT(cfg.numerics.fault_iterator == FaultIteratorKind::Substep,
               "fault_iterator=substep");
   TEST_ASSERT(cfg.numerics.interior_flux == InteriorFlux::Matrix,
               "interior_flux=matrix");
   // Defaults (raw / one-shot / scalar) are exercised by every other config
   // test in this file (none of which set these keys) — those still pass,
   // proving the new selectors default to current behaviour (no regression).
}

// NUM-2 (Phase 6 req 3): interior_flux="matrix" + mixed_flux != "none" aborts.
static void T_24_matrix_mixed_flux_aborts()
{
   std::cout << "\n[NUM-2] interior_flux=matrix + mixed_flux=adjacent aborts\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=32e9
rho=2670
[pore_pressure]
P_p_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
[stress]
kind = "constant_tensor"
sigma_xx_pa=0
sigma_yy_pa=0
sigma_zz_pa=0
sigma_xy_pa=0
sigma_yz_pa=0
sigma_xz_pa=0
[numerics]
ader_order=2
mixed_flux="adjacent"
cfl=0.5
interior_flux="matrix"
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   TEST_ASSERT(ParseAbortsInChild(toml),
               "interior_flux=matrix + mixed_flux=adjacent must abort");
}

// NUC-1 (Phase 6 req 6): [nucleation] kind=gradual_overstress_compact_circular
// (TPV102/104) parses into compact_circular.
static void T_25_nucleation_compact_circular_parses()
{
   std::cout << "\n[NUC-1] [nucleation] gradual_overstress_compact_circular parses\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=32e9
rho=2670
[pore_pressure]
P_p_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
[stress]
kind = "constant_tensor"
sigma_xx_pa=0
sigma_yy_pa=0
sigma_zz_pa=0
sigma_xy_pa=0
sigma_yz_pa=0
sigma_xz_pa=0
[numerics]
ader_order=2
mixed_flux="none"
cfl=0.5
[nucleation]
kind = "gradual_overstress_compact_circular"
[nucleation.gradual_overstress_compact_circular]
center_x_m = 0.0
center_y_m = 0.0
center_z_m = -7500.0
radius_m = 3000.0
delta_tau_pa = 45.0e6
T_nuc_s = "1.0s"
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.nucleation.enabled, "nucleation enabled");
   TEST_ASSERT(cfg.nucleation.kind
               == NucleationKind::GradualOverstressCompactCircular,
               "kind == GradualOverstressCompactCircular");
   TEST_ASSERT(cfg.nucleation.compact_circular.radius_m == 3000.0,
               "compact_circular.radius_m == 3000");
   TEST_ASSERT(cfg.nucleation.compact_circular.delta_tau_pa == 45.0e6,
               "compact_circular.delta_tau_pa == 45 MPa");
   TEST_ASSERT(cfg.nucleation.compact_circular.T_nuc_s == 1.0,
               "compact_circular.T_nuc_s == 1 s");
}

// NUC-2 (Phase 6 req 6): [nucleation] kind=instantaneous_overstress_circular
// (TPV31) parses into instantaneous_circular.
static void T_26_nucleation_instantaneous_circular_parses()
{
   std::cout << "\n[NUC-2] [nucleation] instantaneous_overstress_circular parses\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=32e9
rho=2670
[pore_pressure]
P_p_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
[stress]
kind = "constant_tensor"
sigma_xx_pa=0
sigma_yy_pa=0
sigma_zz_pa=0
sigma_xy_pa=0
sigma_yz_pa=0
sigma_xz_pa=0
[numerics]
ader_order=2
mixed_flux="none"
cfl=0.5
[nucleation]
kind = "instantaneous_overstress_circular"
[nucleation.instantaneous_overstress_circular]
center_x_m = 0.0
center_y_m = 0.0
center_z_m = -10000.0
radius_m = 1400.0
taper_m = 200.0
delta_tau_pa = 11.6e6
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.nucleation.kind
               == NucleationKind::InstantaneousOverstressCircular,
               "kind == InstantaneousOverstressCircular");
   TEST_ASSERT(cfg.nucleation.instantaneous_circular.radius_m == 1400.0,
               "instantaneous_circular.radius_m == 1400");
   TEST_ASSERT(cfg.nucleation.instantaneous_circular.taper_m == 200.0,
               "instantaneous_circular.taper_m == 200");
   TEST_ASSERT(cfg.nucleation.instantaneous_circular.delta_tau_pa == 11.6e6,
               "instantaneous_circular.delta_tau_pa == 11.6 MPa");
}

// =====================================================================
//  Phase 6 req 1: [problem]/[boundary]/[fault_geometry]/[hypocenter]/
//  [material] config blocks + their req-7 guards.
// =====================================================================

// CFG1-1: all five new blocks parse; values round-trip.
static void T_40_new_config_blocks_parse()
{
   std::cout << "\n[CFG1-1] [problem]/[boundary]/[fault_geometry]/"
                "[hypocenter]/[material] parse\n";
   const std::string toml = MinimalLSWHeader() + MinimalLSWBlock()
      + "[problem]\ntag=\"tpv205\"\n"
        "[boundary]\nfault_attr=3\nnatural_attrs=[1,2]\nabsorbing_attrs=[4,5,6]\n"
        "[fault_geometry]\nref_normal=[0.0,-1.0,0.0]\nup=[0.0,0.0,1.0]\nkind=\"planar\"\n"
        "[hypocenter]\nx_m=0.0\ny_m=0.0\nz_m=-7500.0\n"
        "nucleation_radius_m=1500.0\nnucleation_taper_m=0.0\n"
        "[material]\nkind=\"constant\"\n";
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.problem.tag == "tpv205", "problem.tag round-trip");
   TEST_ASSERT(cfg.boundary.fault_attr == 3, "boundary.fault_attr == 3");
   TEST_ASSERT(cfg.boundary.natural_attrs.size() == 2
               && cfg.boundary.natural_attrs[0] == 1
               && cfg.boundary.natural_attrs[1] == 2, "natural_attrs == [1,2]");
   TEST_ASSERT(cfg.boundary.absorbing_attrs.size() == 3, "absorbing_attrs size 3");
   TEST_ASSERT(cfg.fault_geometry.ref_normal[1] == -1.0, "ref_normal[1] == -1");
   TEST_ASSERT(cfg.fault_geometry.up[2] == 1.0, "up[2] == 1");
   TEST_ASSERT(cfg.fault_geometry.kind == "planar", "fault_geometry.kind round-trip");
   TEST_ASSERT(cfg.hypocenter.z_m == -7500.0, "hypocenter.z_m == -7500");
   TEST_ASSERT(cfg.hypocenter.nucleation_radius_m == 1500.0,
               "nucleation_radius_m == 1500");
   TEST_ASSERT(cfg.material.kind == MaterialKind::Constant,
               "material.kind == Constant");
}

// CFG1-2: defaults preserved when the new blocks are absent (no SAFS
// regression — existing TOMLs set none of these).
static void T_41_new_blocks_default_when_absent()
{
   std::cout << "\n[CFG1-2] absent new blocks -> struct defaults\n";
   const std::string toml = MinimalLSWHeader() + MinimalLSWBlock();
   const auto cfg = ParseSpatialFrictionConfigString(toml);
   TEST_ASSERT(cfg.problem.tag.empty(), "problem.tag default empty");
   TEST_ASSERT(cfg.boundary.fault_attr == -1, "boundary.fault_attr default -1");
   TEST_ASSERT(cfg.boundary.natural_attrs.empty(), "natural_attrs default empty");
   TEST_ASSERT(cfg.fault_geometry.ref_normal[0] == 0.0
               && cfg.fault_geometry.ref_normal[1] == -1.0
               && cfg.fault_geometry.ref_normal[2] == 0.0,
               "ref_normal default (0,-1,0)");
   TEST_ASSERT(cfg.fault_geometry.up[0] == 0.0
               && cfg.fault_geometry.up[1] == 0.0
               && cfg.fault_geometry.up[2] == 1.0,
               "up default (0,0,1)");
   TEST_ASSERT(cfg.material.kind == MaterialKind::Constant,
               "material default Constant");
}

// CFG1-3 (req 7): overlapping [boundary] attribute sets / bad fault_attr abort.
static void T_42_boundary_overlap_aborts()
{
   std::cout << "\n[CFG1-3] overlapping/invalid boundary attributes abort\n";
   const std::string fault_in_natural = MinimalLSWHeader() + MinimalLSWBlock()
      + "[boundary]\nfault_attr=3\nnatural_attrs=[1,3]\n";
   TEST_ASSERT(ParseAbortsInChild(fault_in_natural),
               "fault_attr in natural_attrs must abort");
   const std::string nat_abs_overlap = MinimalLSWHeader() + MinimalLSWBlock()
      + "[boundary]\nfault_attr=3\nnatural_attrs=[1,2]\nabsorbing_attrs=[2,4]\n";
   TEST_ASSERT(ParseAbortsInChild(nat_abs_overlap),
               "natural/absorbing overlap must abort");
   const std::string bad_fault_attr = MinimalLSWHeader() + MinimalLSWBlock()
      + "[boundary]\nfault_attr=0\n";
   TEST_ASSERT(ParseAbortsInChild(bad_fault_attr),
               "fault_attr=0 (non-positive) must abort");
}

// CFG1-4 (req 7): non-unit-norm or parallel [fault_geometry] axes abort.
static void T_43_fault_geometry_guards_abort()
{
   std::cout << "\n[CFG1-4] non-unit / parallel fault_geometry aborts\n";
   const std::string non_unit = MinimalLSWHeader() + MinimalLSWBlock()
      + "[fault_geometry]\nref_normal=[0.0,-2.0,0.0]\nup=[0.0,0.0,1.0]\n";
   TEST_ASSERT(ParseAbortsInChild(non_unit),
               "non-unit ref_normal must abort");
   const std::string parallel = MinimalLSWHeader() + MinimalLSWBlock()
      + "[fault_geometry]\nref_normal=[0.0,0.0,1.0]\nup=[0.0,0.0,1.0]\n";
   TEST_ASSERT(ParseAbortsInChild(parallel),
               "up parallel to ref_normal must abort");
}

// CFG1-5 (req 7, R-008): hypocenter z > 0 with up[2] > 0 aborts.
static void T_44_hypocenter_positive_z_aborts()
{
   std::cout << "\n[CFG1-5] hypocenter z>0 with up[2]>0 aborts\n";
   const std::string toml = MinimalLSWHeader() + MinimalLSWBlock()
      + "[fault_geometry]\nref_normal=[0.0,-1.0,0.0]\nup=[0.0,0.0,1.0]\n"
        "[hypocenter]\nx_m=0.0\ny_m=0.0\nz_m=500.0\n";
   TEST_ASSERT(ParseAbortsInChild(toml),
               "z_m > 0 with up[2] > 0 must abort");
}

// CFG1-6: [material] kind selectors parse; bad kind / missing path abort.
static void T_45_material_kinds()
{
   std::cout << "\n[CFG1-6] material kind selectors\n";
   const std::string dp = MinimalLSWHeader() + MinimalLSWBlock()
      + "[material]\nkind=\"depth_profile_1d\"\nprofile_csv=\"/tmp/p.csv\"\n";
   const auto cfg = ParseSpatialFrictionConfigString(dp);
   TEST_ASSERT(cfg.material.kind == MaterialKind::DepthProfile1D,
               "kind=depth_profile_1d");
   TEST_ASSERT(cfg.material.profile_csv == "/tmp/p.csv", "profile_csv round-trip");

   const std::string bad = MinimalLSWHeader() + MinimalLSWBlock()
      + "[material]\nkind=\"bogus\"\n";
   TEST_ASSERT(ParseAbortsInChild(bad), "bad material.kind must abort");

   const std::string dp_no_csv = MinimalLSWHeader() + MinimalLSWBlock()
      + "[material]\nkind=\"depth_profile_1d\"\n";
   TEST_ASSERT(ParseAbortsInChild(dp_no_csv),
               "depth_profile_1d without profile_csv must abort");

   const std::string sc_no_path = MinimalLSWHeader() + MinimalLSWBlock()
      + "[material]\nkind=\"sidecar_hdf5\"\n";
   TEST_ASSERT(ParseAbortsInChild(sc_no_path),
               "sidecar_hdf5 without sidecar_path must abort");
}

// CFG1-7 (deferred req-3 guard): interior_flux=matrix + Constant material aborts.
static void T_46_matrix_requires_nonconstant_material()
{
   std::cout << "\n[CFG1-7] interior_flux=matrix + Constant material aborts\n";
   const std::string toml = R"TOML(
[meta]
schema_version = 1
law = "slip_weakening"
[material_constant_fallback]
lambda=32e9
mu=32e9
rho=2670
[pore_pressure]
P_p_pa=0
[mesh]
path="/dev/null"
order=1
[velocity]
model="cvmh"
dataset_root="/tmp/x"
[stress]
kind = "constant_tensor"
sigma_xx_pa=0
sigma_yy_pa=0
sigma_zz_pa=0
sigma_xy_pa=0
sigma_yz_pa=0
sigma_xz_pa=0
[numerics]
ader_order=2
mixed_flux="none"
cfl=0.5
interior_flux="matrix"
[material]
kind="constant"
[time]
tfinal="12s"
[output]
output_dir="out"
[friction.slip_weakening]
mu_s_default=1.1
mu_d_default=0.5
d_c_default=0.5
cohesion_default=0
)TOML";
   TEST_ASSERT(ParseAbortsInChild(toml),
               "interior_flux=matrix with Constant material must abort");
}

int main(int, char**)
{
#ifndef SEAS_USE_TOML
   std::cout << "SEAS_USE_TOML not defined — skipping test.\n";
   return 0;
#else
   std::cout << "Running Phase 1 test_spatial_friction_config\n";
   T_1_minimal_lsw_parses_geoffrey2010();
   T_2_schema_version_mismatch_aborts();
   T_3_missing_law_aborts();
   T_4_both_friction_blocks_aborts();
   T_5_d_o_alias_works();
   T_6_d_o_collision_aborts();
   T_7_r114_default_guard();
   T_8_r114_rule_guard();
   T_9_stress_mode_conflict_const_with_path();
   T_10_stress_mode_conflict_sidecar_with_sigma();
   T_20_fault_local_prestress_parses();
   T_21_fault_local_with_sigma_aborts();
   T_22_fault_local_patches_parse();
   T_23_numerics_selectors_parse();
   T_24_matrix_mixed_flux_aborts();
   T_25_nucleation_compact_circular_parses();
   T_26_nucleation_instantaneous_circular_parses();
   T_11_time_parser();
   T_12_material_fallback_bounds();
   T_13_missing_stress_block_aborts();
   T_14_missing_time_block_aborts();
   T_15_dt_initial_zero_aborts();
   T_16_use_pml_string_aborts();
   T_17_missing_pore_pressure_block_aborts();
   T_18_dt_initial_literal_negative_aborts();
   T_19_paraview_fault_dt_zero_aborts();
   T_20_stress_kind_missing_aborts();
   T_21_nucleation_kind_must_be_gradual_overstress();
   T_22_nucleation_kind_gradual_overstress_parses();
   T_23_nucleation_kind_typo_aborts();
   T_24_paraview_default_flip();
   T_25_paraview_extended_defaults();
   T_26_velocity_use_sidecar_false_relaxes_dataset_root();
   T_27_velocity_use_sidecar_true_default_still_requires_root();
   T_28_sigma_n_strength_floor_parses();
   T_29_sigma_n_strength_floor_absent_disabled();
   T_30_sigma_n_strength_floor_negative_aborts();
   T_31_sigma_n_strength_floor_zero_allowed();
   T_32_sigma_n_floor_misnested_aborts();
   // Phase 6 req 1 + req-7 guards on the new TPV config blocks.
   T_40_new_config_blocks_parse();
   T_41_new_blocks_default_when_absent();
   T_42_boundary_overlap_aborts();
   T_43_fault_geometry_guards_abort();
   T_44_hypocenter_positive_z_aborts();
   T_45_material_kinds();
   T_46_matrix_requires_nonconstant_material();
   std::cout << "\n========================================\n";
   std::cout << "Phase 1 test_spatial_friction_config: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
#endif
}
