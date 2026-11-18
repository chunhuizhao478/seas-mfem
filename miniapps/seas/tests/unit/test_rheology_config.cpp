// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Rheology plan rev 1, unit test U10 (no mesh): the deck tables [rheology]
// and [initial_stress] (rheology/rheology_config.hpp) and the allowed
// combinations of rheology/rheology_stack.hpp.
//
//   absent table or enabled = false gives an empty stack; an enabled model
//   is built with the deck's values; unknown keys, missing or mistyped
//   values, refused combinations and the R7/R8 placeholders each stop with
//   their own message.
//
// Run: ./seas_test_rheology_config

#include "mfem.hpp"
#include "../../rheology/rheology_config.hpp"
#include "test_macros.hpp"

#include <iostream>
#include <sstream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{

toml::value Deck(const std::string &text)
{
   std::istringstream is(text);
   return toml::parse(is, "test_deck.toml");
}

std::string RheologyError(const std::string &text)
{
   RheologyConfig cfg;
   return ParseRheologyConfig(Deck(text), cfg);
}

std::string InitialStressError(const std::string &text)
{
   const real_t up[3] = {0.0, 0.0, 1.0};
   bool present = false;
   InitialStressSpec spec;
   return ParseInitialStressConfig(Deck(text), up, present, spec);
}

bool Has(const std::string &msg, const std::string &part)
{
   const bool ok = msg.find(part) != std::string::npos;
   if (!ok) { std::cerr << "  message was: \"" << msg << "\"\n"; }
   return ok;
}

const char *kDp =
   "[rheology]\n"
   "representation = \"nodal_pm1\"\n"
   "[rheology.drucker_prager]\n"
   "enabled = true\n"
   "cohesion_pa = 1.36e6\n"
   "bulk_friction = 0.1934\n"
   "relaxation_time_s = 0.03\n";

const char *kLitho =
   "[initial_stress]\n"
   "kind = \"lithostatic_ratio\"\n"
   "rho_rock_kg_m3 = 2670.0\n"
   "rho_fluid_kg_m3 = 1000.0\n"
   "gravity_m_s2 = 9.8\n"
   "surface_elevation_m = 0.0\n"
   "b_xx = 0.926793\n"
   "b_yy = 1.073206\n"
   "b_zz = 0.0\n"
   "b_xy = -0.169029\n"
   "b_yz = 0.0\n"
   "b_xz = 0.0\n"
   "taper_top_depth_m = 15000.0\n"
   "taper_bottom_depth_m = 20000.0\n";

/// `text` with the line starting with `key` replaced by `line` ("" drops it).
std::string Edit(const std::string &text, const std::string &key, const std::string &line)
{
   std::istringstream is(text);
   std::string out, l;
   while (std::getline(is, l))
   {
      if (l.compare(0, key.size() + 1, key + " ") == 0)
      {
         if (!line.empty()) { out += line + "\n"; }
      }
      else { out += l + "\n"; }
   }
   return out;
}

} // namespace

void TestOffAndOn()
{
   std::cout << "\n=== U10: off, disabled, enabled ===\n";
   {
      RheologyConfig cfg;
      const std::string e = ParseRheologyConfig(Deck("[mesh]\norder = 2\n"), cfg);
      TEST_ASSERT(e.empty() && !cfg.present && !cfg.Any(), "U10 no [rheology]: off");
      TEST_ASSERT(BuildRheologyStack(cfg)->Empty(), "U10 no [rheology]: empty stack");
   }
   {
      RheologyConfig cfg;
      const std::string off = Edit(kDp, "enabled", "enabled = false");
      const std::string e = ParseRheologyConfig(Deck(off), cfg);
      TEST_ASSERT(e.empty() && cfg.present && !cfg.Any(), "U10 enabled = false: off");
      TEST_ASSERT(BuildRheologyStack(cfg)->Empty(), "U10 enabled = false: empty stack");
   }
   {
      RheologyConfig cfg;
      const std::string e = ParseRheologyConfig(Deck(std::string(kDp)
                                                     + "[rheology.viscoelastic]\nenabled = false\n"
                                                     + "[rheology.damage_breakage]\nenabled = false\n"),
                                                cfg);
      TEST_ASSERT(e.empty() && cfg.drucker_prager, "U10 TPV27 table parses");
      TEST_ASSERT(cfg.dp.cohesion == 1.36e6 && cfg.dp.bulk_friction == 0.1934
                  && cfg.dp.relaxation_time == 0.03, "U10 Drucker-Prager values read");
      std::unique_ptr<RheologyStack> st = BuildRheologyStack(cfg);
      TEST_ASSERT(st->Size() == 1 && st->Describe() == "drucker_prager"
                  && st->NumStateVars() == DruckerPrager::kNumState
                  && !st->IsLinearInStrain(), "U10 stack = {drucker_prager}, 7 state variables");
   }
}

void TestRefusals()
{
   std::cout << "\n=== U10: refusals name the key ===\n";
   const std::string dp = kDp;
   TEST_ASSERT(Has(RheologyError(dp + "foo = 1\n"), "[rheology.drucker_prager] unknown key 'foo'"),
               "U10 unknown key in [rheology.drucker_prager]");
   TEST_ASSERT(Has(RheologyError("[rheology]\nplastic = true\n"), "[rheology] unknown key 'plastic'"),
               "U10 unknown key in [rheology]");
   TEST_ASSERT(Has(RheologyError(Edit(dp, "cohesion_pa", "")), "needs cohesion_pa"),
               "U10 missing cohesion_pa");
   TEST_ASSERT(Has(RheologyError(Edit(dp, "relaxation_time_s", "")), "needs relaxation_time_s"),
               "U10 missing relaxation_time_s");
   TEST_ASSERT(Has(RheologyError(Edit(dp, "enabled", "")), "needs enabled = true or false"),
               "U10 missing enabled");
   TEST_ASSERT(Has(RheologyError(Edit(dp, "bulk_friction", "bulk_friction = \"0.19\"")),
                   "bulk_friction must be a number"), "U10 mistyped bulk_friction");
   TEST_ASSERT(Has(RheologyError(Edit(dp, "enabled", "enabled = 1")),
                   "enabled must be true or false"), "U10 mistyped enabled");
   TEST_ASSERT(Has(RheologyError(Edit(dp, "relaxation_time_s", "relaxation_time_s = -0.03")),
                   "relaxation_time_s must be >= 0"), "U10 negative relaxation time");
   TEST_ASSERT(Has(RheologyError(Edit(dp, "cohesion_pa", "cohesion_pa = -1.0")),
                   "cohesion_pa must be >= 0"), "U10 negative cohesion");
   TEST_ASSERT(Has(RheologyError(Edit(dp, "representation", "representation = \"quadrature\"")),
                   "only \"nodal_pm1\" is implemented"), "U10 representation other than nodal_pm1");
   TEST_ASSERT(Has(RheologyError("[rheology]\n[rheology.damage_breakage]\nenabled = true\n"),
                   "continuum damage-breakage is not implemented yet (rheology plan phase R8)"),
               "U10 damage-breakage refused until R8");
}

// Attenuation plan rev6, section 5.10: the [rheology.viscoelastic] keys.
const char *kVe =
   "[rheology]\n"
   "[rheology.viscoelastic]\n"
   "enabled = true\n"
   "mechanisms = 3\n"
   "f_central_hz = 0.1\n"
   "f_ratio = 100.0\n"
   "Q_s_over_vs = 0.05\n"
   "Q_p_over_Q_s = 1.5\n"
   "moduli = \"relaxed\"\n"
   "width_memory = true\n"
   "in_pml_band = false\n"
   "handback_tol_pa = 0.0\n";

void TestViscoelastic()
{
   std::cout << "\n=== U10: [rheology.viscoelastic] (attenuation plan rev6, 5.10) ===\n";
   const std::string ve = kVe;
   {
      RheologyConfig cfg;
      const std::string e = ParseRheologyConfig(Deck(ve), cfg);
      TEST_ASSERT(e.empty() && cfg.viscoelastic && cfg.Any(), "U10 viscoelastic table parses");
      TEST_ASSERT(cfg.gmb_band.L == 3 && cfg.gmb_band.f_c == 0.1 && cfg.gmb_band.R == 100.0
                  && cfg.gmb_q.from_vs && cfg.gmb_q.q_s_over_vs == 0.05
                  && cfg.gmb_q.q_p_over_q_s == 1.5 && cfg.width_memory && !cfg.in_pml_band
                  && cfg.handback_tol == 0.0, "U10 viscoelastic values read");
      std::unique_ptr<RheologyStack> st = BuildRheologyStack(cfg);
      TEST_ASSERT(st->Size() == 1 && st->Describe() == "viscoelastic"
                  && st->NumStateVars() == 8 * 3 + 6 && st->IsLinearInStrain(),
                  "U10 stack = {viscoelastic}, 8 L + 6 state variables, linear");
   }
   {
      RheologyConfig cfg;
      const std::string fixed = Edit(Edit(ve, "Q_s_over_vs", "Q_s = 50.0"), "Q_p_over_Q_s",
                                     "Q_p = 75.0");
      const std::string e = ParseRheologyConfig(Deck(fixed), cfg);
      TEST_ASSERT(e.empty() && !cfg.gmb_q.from_vs && cfg.gmb_q.Q_s == 50.0
                  && cfg.gmb_q.Q_p == 75.0, "U10 fixed Q_s, Q_p read");
   }
   {
      RheologyConfig cfg;
      const std::string e = ParseRheologyConfig(Deck("[rheology]\n[rheology.viscoelastic]\n"
                                                     "enabled = false\n"), cfg);
      TEST_ASSERT(e.empty() && !cfg.Any() && BuildRheologyStack(cfg)->Empty(),
                  "U10 viscoelastic enabled = false: empty stack");
   }
   const char *each[] = {"mechanisms", "f_central_hz", "f_ratio", "moduli", "width_memory",
                         "in_pml_band", "handback_tol_pa"
                        };
   for (const char *k : each)
   {
      TEST_ASSERT(Has(RheologyError(Edit(ve, k, "")), std::string("needs ") + k),
                  std::string("U10 missing ") + k);
   }
   TEST_ASSERT(Has(RheologyError(Edit(ve, "Q_s_over_vs", "")),
                   "either Q_s_over_vs and Q_p_over_Q_s, or Q_s and Q_p"),
               "U10 half a Q pair");
   TEST_ASSERT(Has(RheologyError(ve + "Q_s = 50.0\n"),
                   "either Q_s_over_vs and Q_p_over_Q_s, or Q_s and Q_p"),
               "U10 both Q pairs");
   TEST_ASSERT(Has(RheologyError(ve + "q_s = 50\n"), "[rheology.viscoelastic] unknown key 'q_s'"),
               "U10 unknown viscoelastic key");
   TEST_ASSERT(Has(RheologyError(Edit(ve, "mechanisms", "mechanisms = 0")),
                   "mechanisms must be 1 .. 10"), "U10 mechanisms = 0");
   TEST_ASSERT(Has(RheologyError(Edit(ve, "mechanisms", "mechanisms = 3.0")),
                   "mechanisms must be an integer"), "U10 mistyped mechanisms");
   TEST_ASSERT(Has(RheologyError(Edit(ve, "f_ratio", "f_ratio = 1.0")),
                   "f_ratio must be > 1"), "U10 f_ratio = 1 with L = 3");
   TEST_ASSERT(Has(RheologyError(Edit(ve, "f_central_hz", "f_central_hz = 0.0")),
                   "f_central_hz must be > 0"), "U10 f_central_hz = 0");
   TEST_ASSERT(Has(RheologyError(Edit(ve, "moduli", "moduli = \"reference\"")),
                   "moduli = \"reference\" needs the relaxed moduli"), "U10 reference refused");
   TEST_ASSERT(Has(RheologyError(Edit(ve, "moduli", "moduli = \"unrelaxed\"")),
                   "\"relaxed\" or \"reference\""), "U10 unknown moduli");
   TEST_ASSERT(Has(RheologyError(Edit(ve, "in_pml_band", "in_pml_band = true")),
                   "decision D4"), "U10 in_pml_band = true refused");
   TEST_ASSERT(Has(RheologyError(Edit(ve, "handback_tol_pa", "handback_tol_pa = -1.0")),
                   "handback_tol_pa must be >= 0"), "U10 negative handback_tol_pa");
   TEST_ASSERT(Has(RheologyError(Edit(ve, "Q_s_over_vs", "Q_s_over_vs = 0.0")),
                   "Q_s_over_vs must be > 0"), "U10 Q_s_over_vs = 0");
   TEST_ASSERT(Has(RheologyError(ve + std::string(kDp).substr(std::string(kDp).find("[rheology.drucker_prager]"))),
                   "attenuation plan section 2.14"), "U10 viscoelastic + drucker_prager refused");
}

void TestCombinations()
{
   std::cout << "\n=== U10: allowed combinations ===\n";
   using K = RheologyKind;
   TEST_ASSERT(RheologyStack::CheckCombination({K::DruckerPrager}).empty(), "U10 {dp} allowed");
   TEST_ASSERT(RheologyStack::CheckCombination({K::Viscoelastic, K::DruckerPrager}).empty(),
               "U10 {ve, dp} allowed");
   TEST_ASSERT(RheologyStack::CheckCombination({K::DamageBreakage, K::Viscoelastic}).empty(),
               "U10 {ve, db} allowed");
   TEST_ASSERT(Has(RheologyStack::CheckCombination({K::DruckerPrager, K::DamageBreakage}),
                   "cannot both be enabled"), "U10 {dp, db} refused");
   TEST_ASSERT(Has(RheologyStack::CheckCombination({K::DruckerPrager, K::DruckerPrager}),
                   "appears twice"), "U10 a model twice refused");
}

void TestInitialStress()
{
   std::cout << "\n=== U10: [initial_stress] ===\n";
   {
      const real_t up[3] = {0.0, 0.0, 1.0};
      bool present = true;
      InitialStressSpec spec;
      const std::string e = ParseInitialStressConfig(Deck("[mesh]\norder = 1\n"), up, present, spec);
      TEST_ASSERT(e.empty() && !present, "U10 no [initial_stress]: absent");
   }
   {
      const real_t up[3] = {0.0, 0.0, 1.0};
      bool present = false;
      InitialStressSpec spec;
      const std::string e = ParseInitialStressConfig(Deck(kLitho), up, present, spec);
      TEST_ASSERT(e.empty() && present
                  && spec.kind == InitialStressSpec::Kind::LithostaticRatio,
                  "U10 lithostatic_ratio parses");
      TEST_ASSERT(spec.B[0] == 0.926793 && spec.B[4] == 1.073206 && spec.B[1] == -0.169029
                  && spec.B[3] == -0.169029 && spec.B[8] == 0.0, "U10 B read into the code frame");
      InitialStress is(spec);
      const real_t x[3] = {0.0, 0.0, -5000.0};
      real_t s[9];
      is.EvalEffective(x, s);
      TEST_ASSERT(std::abs(s[1] - 13.832e6) < 1e3, "U10 parsed deck gives 13.832 MPa at 5 km");
   }
   const std::string L = kLitho;
   TEST_ASSERT(Has(InitialStressError(Edit(L, "b_zz", "")), "needs b_zz (no default)"),
               "U10 missing b_zz");
   TEST_ASSERT(Has(InitialStressError(Edit(L, "b_zz", "b_zz = 0.1")), "must be horizontal"),
               "U10 vertical B refused");
   TEST_ASSERT(Has(InitialStressError(Edit(L, "taper_top_depth_m", "taper_top_depth_m = 25000.0")),
                   "taper_top_depth_m < taper_bottom_depth_m"), "U10 reversed taper refused");
   TEST_ASSERT(Has(InitialStressError(Edit(L, "kind", "kind = \"gradient\"")),
                   "kind = \"gradient\""), "U10 unknown kind refused");
   TEST_ASSERT(Has(InitialStressError(L + "sigma_xx_pa = 1.0\n"), "unknown key 'sigma_xx_pa'"),
               "U10 tensor key in a lithostatic table refused");
   const std::string ct = "[initial_stress]\nkind = \"constant_tensor\"\n"
                          "sigma_xx_pa = 0.0\nsigma_yy_pa = -120.0e6\nsigma_zz_pa = 0.0\n"
                          "sigma_xy_pa = 70.0e6\nsigma_yz_pa = 0.0\nsigma_xz_pa = 0.0\n";
   TEST_ASSERT(InitialStressError(ct).empty(), "U10 constant_tensor parses");
   TEST_ASSERT(Has(InitialStressError(Edit(ct, "sigma_xz_pa", "")), "needs sigma_xz_pa"),
               "U10 constant_tensor needs all six components");
}

int main()
{
   std::cout << "Rheology deck tests (plan rev 1, U10)\n";
   TestOffAndOn();
   TestRefusals();
   TestViscoelastic();
   TestCombinations();
   TestInitialStress();
   std::cout << "\n" << num_passed << " / " << num_tests << " passed, " << num_failed
             << " failed\n";
   return num_failed == 0 ? 0 : 1;
}
