// Poroelastic plan rev 3, phase P0 (step P0.1): the option.
//
//   U0-4  ParsePoroConfig: a valid table reads back exactly; every required
//         key missing, misspelt or of the wrong type, unknown keys, the
//         source and starting-pressure tables, --poroelastic 1 without the
//         table, the option off with a table present (warning), and the bulk
//         density check.  The plan's [multiphysics] poroelastic = true next
//         to a [multiphysics.poroelastic] table is not valid TOML (checked).
//
// The driver-level checks of the abort texts arrive with the driver parse in
// step P0.2 (Makefile recipe test-poro-p0).
//
// Usage: mpirun -np 1|2 ./seas_test_poro_p0

#include "mfem.hpp"
#include "../poro_config.hpp"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

// ------------------------------------------------------------------- U0-4
static const char *kValid =
   "[multiphysics.poroelastic]\n"
   "enabled = true\n"
   "porosity = 0.4\n"
   "permeability_m2 = 6.0e-13\n"
   "viscosity_pa_s = 1.0e-3\n"
   "tortuosity = 2.0\n"
   "rho_fluid = 1040.0\n"
   "rho_solid = 2080.0\n"
   "bulk_fluid_pa = 2.5e9\n"
   "bulk_solid_pa = 2.0e10\n"
   "fault_hydraulics = \"permeable\"\n"
   "background_p_pa = 0\n"
   "fluid_penalty = 1.0\n"
   "drag_gamma = \"fitted\"\n"
   "velocity_coupling = \"momentum\"\n";

static const char *kInitTable =
   "[multiphysics.poroelastic.initial_pressure]\n"
   "origin_m = [0.0, -1000.0, -7500.0]\n"
   "magnitude_pa = 5.0e7\n"
   "width_per_m2 = 1.0e-7\n";

static const char *kKeys[] = {"enabled", "porosity", "permeability_m2", "viscosity_pa_s",
                              "tortuosity", "rho_fluid", "rho_solid", "bulk_fluid_pa",
                              "bulk_solid_pa", "fault_hydraulics", "background_p_pa",
                              "fluid_penalty", "drag_gamma", "velocity_coupling"};

static poro::PoroParseResult Parse(const std::string &deck, int cli = -1)
{
   std::istringstream is(deck);
   const toml::value root = toml::parse(is, "u04.toml");
   return poro::ParsePoroConfig(root, cli);
}

/// kValid with the line of `key` removed or replaced by `line`.
static std::string Edit(const std::string &key, const std::string &line)
{
   std::istringstream is(kValid);
   std::ostringstream os;
   for (std::string l; std::getline(is, l);)
   {
      if (l.rfind(key + " =", 0) == 0) { if (!line.empty()) { os << line << "\n"; } }
      else { os << l << "\n"; }
   }
   return os.str();
}

static bool Has(const std::string &s, const std::string &sub) { return s.find(sub) != std::string::npos; }

static void TestU04()
{
   if (g_rank == 0) { std::cout << "--- U0-4: the parser\n"; }
   // a valid table reads back exactly
   {
      const poro::PoroParseResult r = Parse(std::string(kValid) + kInitTable);
      const poro::PoroConfig &c = r.cfg;
      TEST_ASSERT(r.error.empty() && c.enabled, "U0-4 valid tpv5-poro table parses (" << r.error << ")");
      TEST_ASSERT(c.porosity == 0.4 && c.permeability_m2 == 6.0e-13 && c.viscosity_pa_s == 1.0e-3
                  && c.tortuosity == 2.0 && c.rho_fluid == 1040.0 && c.rho_solid == 2080.0
                  && c.bulk_fluid_pa == 2.5e9 && c.bulk_solid_pa == 2.0e10
                  && c.fault_hydraulics == "permeable" && c.background_p_pa == 0.0
                  && c.fluid_penalty == 1.0 && c.drag_gamma == "fitted"
                  && c.velocity_coupling == "momentum",
                  "U0-4 every value reads back exactly (integer 0 accepted as a number)");
      TEST_ASSERT(c.initial_pressure.on && c.initial_pressure.origin_m[1] == -1000.0
                  && c.initial_pressure.origin_m[2] == -7500.0
                  && c.initial_pressure.magnitude_pa == 5.0e7
                  && c.initial_pressure.width_per_m2 == 1.0e-7 && !c.source.on,
                  "U0-4 the starting-pressure table reads back; no source");
      TEST_ASSERT(std::abs(c.BulkDensity() - 1664.0) <= 1e-12 * 1664.0,
                  "U0-4 bulk density (1 - phi) rho_s + phi rho_f = 1664 kg/m^3");
      TEST_ASSERT(poro::CheckBulkDensity(c, 1664.0).empty(), "U0-4 rho 1664 accepted");
      const std::string e = poro::CheckBulkDensity(c, 2670.0);
      TEST_ASSERT(Has(e, "rho = 2670") && Has(e, "1664"), "U0-4 rho mismatch refused: " << e);
   }
   // every required key missing, then of the wrong type
   int miss_ok = 0, type_ok = 0;
   const int nk = sizeof(kKeys) / sizeof(kKeys[0]);
   for (int k = 0; k < nk; k++)
   {
      const std::string key = kKeys[k];
      const poro::PoroParseResult r = Parse(Edit(key, ""));
      if (Has(r.error, key + ": missing")) { miss_ok++; }
      else if (g_rank == 0) { std::cout << "    missing " << key << " gave: '" << r.error << "'\n"; }
      const bool is_string = key == "fault_hydraulics" || key == "drag_gamma"
                             || key == "velocity_coupling";
      const std::string bad = key + (key == "enabled" ? " = 1" : (is_string ? " = 1.0" : " = \"x\""));
      const poro::PoroParseResult r2 = Parse(Edit(key, bad));
      if (Has(r2.error, key + ": must be")) { type_ok++; }
      else if (g_rank == 0) { std::cout << "    wrong type " << key << " gave: '" << r2.error << "'\n"; }
   }
   TEST_ASSERT(miss_ok == nk, "U0-4 each of the " << nk << " keys missing: error names it ("
               << miss_ok << "/" << nk << ")");
   TEST_ASSERT(type_ok == nk, "U0-4 each of the " << nk << " keys of the wrong type: error names it ("
               << type_ok << "/" << nk << ")");
   // misspelt and unknown keys
   {
      const poro::PoroParseResult r = Parse(Edit("porosity", "porosty = 0.4"));
      TEST_ASSERT(Has(r.error, "porosty: unknown key"), "U0-4 misspelt key named: " << r.error);
      const poro::PoroParseResult r2 = Parse(std::string(kValid) + "[multiphysics]\nthermal = 1\n");
      TEST_ASSERT(Has(r2.error, "[multiphysics] thermal: unknown key"),
                  "U0-4 unknown key in [multiphysics] named: " << r2.error);
      const poro::PoroParseResult r3 = Parse(std::string(kValid) + kInitTable + "radius_m = 3.0\n");
      TEST_ASSERT(Has(r3.error, "initial_pressure] radius_m: unknown key"),
                  "U0-4 unknown key in the starting-pressure table named: " << r3.error);
   }
   // values out of range and the choices
   {
      const poro::PoroParseResult r = Parse(Edit("porosity", "porosity = 1.0"));
      TEST_ASSERT(Has(r.error, "porosity: must be in (0, 1)"), "U0-4 porosity 1 refused: " << r.error);
      const poro::PoroParseResult r2 = Parse(Edit("fault_hydraulics", "fault_hydraulics = \"sealed\""));
      TEST_ASSERT(Has(r2.error, "arrives in step P6"), "U0-4 sealed fault refused until P6: " << r2.error);
      const poro::PoroParseResult r3 = Parse(Edit("drag_gamma", "drag_gamma = \"exact\""));
      TEST_ASSERT(Has(r3.error, "drag_gamma: \"exact\" is not one of"), "U0-4 unknown choice named: " << r3.error);
      const poro::PoroParseResult r4 = Parse(Edit("tortuosity", "tortuosity = 0.5"));
      TEST_ASSERT(Has(r4.error, "tortuosity: must be >= 1"), "U0-4 tortuosity < 1 refused: " << r4.error);
   }
   // the source table
   {
      const std::string g = "[multiphysics.poroelastic.source]\nshape = \"gaussian\"\n"
                            "origin_m = [0.0, 0.0, -5000.0]\nrate = 1.0e-6\n";
      const poro::PoroParseResult r1 = Parse(std::string(kValid) + g + "width_per_m2 = 1.0e-7\n"
                                             "time_function = \"smooth_step\"\nrise_time_s = 0.1\n");
      TEST_ASSERT(r1.error.empty() && r1.cfg.source.on && r1.cfg.source.rise_time_s == 0.1,
                  "U0-4 gaussian smooth-step source parses (" << r1.error << ")");
      const poro::PoroParseResult r2 = Parse(std::string(kValid) + g + "time_function = \"step\"\n");
      TEST_ASSERT(Has(r2.error, "width_per_m2: missing"), "U0-4 gaussian without width refused: " << r2.error);
      const poro::PoroParseResult r3 = Parse(std::string(kValid) + g + "width_per_m2 = 1.0e-7\n"
                                             "time_function = \"step\"\nrise_time_s = 0.1\n");
      TEST_ASSERT(Has(r3.error, "rise_time_s: applies to"), "U0-4 rise time with a step refused: " << r3.error);
      const poro::PoroParseResult r4 = Parse(std::string(kValid) + g + "width_per_m2 = 1.0e-7\n"
                                             "time_function = \"ramp\"\n");
      TEST_ASSERT(Has(r4.error, "rise_time_s: missing"), "U0-4 ramp without rise time refused: " << r4.error);
      const std::string w = "[multiphysics.poroelastic.source]\nshape = \"well\"\n"
                            "origin_m = [0.0, 0.0, -5000.0]\nrate = 1.0e-3\ntime_function = \"step\"\n";
      const poro::PoroParseResult r5 = Parse(std::string(kValid) + w);
      TEST_ASSERT(r5.error.empty() && r5.cfg.source.shape == "well", "U0-4 well source parses ("
                  << r5.error << ")");
      const poro::PoroParseResult r6 = Parse(std::string(kValid) + w + "width_per_m2 = 1.0e-7\n");
      TEST_ASSERT(Has(r6.error, "width_per_m2: applies to shape = \"gaussian\" only"),
                  "U0-4 width with a well refused: " << r6.error);
   }
   // on and off
   {
      const poro::PoroParseResult r1 = Parse("[numerics]\ncfl = 0.5\n", 1);
      TEST_ASSERT(Has(r1.error, "--poroelastic 1 needs the [multiphysics.poroelastic] table"),
                  "U0-4 --poroelastic 1 without the table refused: " << r1.error);
      const poro::PoroParseResult r2 = Parse("[numerics]\ncfl = 0.5\n");
      TEST_ASSERT(r2.error.empty() && !r2.cfg.enabled && r2.cfg.warning.empty(),
                  "U0-4 no table: off, no warning");
      const poro::PoroParseResult r3 = Parse(Edit("enabled", "enabled = false"));
      TEST_ASSERT(r3.error.empty() && !r3.cfg.enabled && Has(r3.cfg.warning, "is ignored"),
                  "U0-4 enabled = false with a table: off, warning '" << r3.cfg.warning << "'");
      const poro::PoroParseResult r4 = Parse(kValid, 0);
      TEST_ASSERT(r4.error.empty() && !r4.cfg.enabled && Has(r4.cfg.warning, "--poroelastic 0"),
                  "U0-4 --poroelastic 0 overrides enabled = true, with a warning");
      const poro::PoroParseResult r5 = Parse(Edit("enabled", "enabled = false"), 1);
      TEST_ASSERT(r5.error.empty() && r5.cfg.enabled, "U0-4 --poroelastic 1 overrides enabled = false");
      const poro::PoroParseResult r6 = Parse(Edit("enabled", "enabled = false")
                                             + "porosty = 3\n");
      TEST_ASSERT(r6.error.empty() && !r6.cfg.enabled,
                  "U0-4 with the option off the table is not checked (ignored, as warned)");
      bool threw = false;
      try
      {
         std::istringstream is("[multiphysics]\nporoelastic = true\n[multiphysics.poroelastic]\n"
                               "porosity = 0.4\n");
         (void)toml::parse(is, "plan_layout.toml");
      }
      catch (const std::exception &) { threw = true; }
      TEST_ASSERT(threw, "U0-4 the plan's layout ([multiphysics] poroelastic = true and a "
                  "[multiphysics.poroelastic] table) is invalid TOML, hence the key 'enabled'");
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   if (g_rank == 0)
   {
      std::cout << "=== Poroelastic P0 (plan rev 3) on " << Mpi::WorldSize() << " rank(s) ===\n";
   }
   if (g_rank == 0) { TestU04(); }   // pure parsing, rank 0
   int nf = num_failed;
   MPI_Allreduce(MPI_IN_PLACE, &nf, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "=== Summary: " << num_passed << "/" << num_tests << " passed, "
                << num_failed << " failed ===\n";
   }
   return nf == 0 ? 0 : 1;
}
