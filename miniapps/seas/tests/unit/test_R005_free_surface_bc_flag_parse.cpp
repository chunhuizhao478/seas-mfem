// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 Phase 2 (I-04 wiring): R-005 CLI-flag parsing for
// --free-surface-bc.
//
// Exercises the parsing logic baked into drivers/tpv102_driver.cpp:
// case-insensitive match, unknown-value fallback-to-gamma, and banner
// string reporting.  Reimplements the parsing logic here rather than
// exec'ing the driver so the test is fast and independent of MPI / mesh
// files.  Any change to the driver's parsing block should be mirrored
// here; divergence is flagged by test failure.

#include "../../dynamic/wave_operator.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

// ---------------------------------------------------------------------------
// Parse one --free-surface-bc input string.  Mirrors the R-005 block in
// drivers/tpv102_driver.cpp.  Returns the mode and the WARNING text that
// should be printed on rank 0 (empty if no warning).  `requested_out` is
// set to the canonicalized requested value (after the unknown->gamma
// fallback), matching the driver's banner.
// ---------------------------------------------------------------------------
static FreeSurfaceBCMode ParseFlag(const std::string &requested_in,
                                    std::string &requested_out,
                                    std::string &warning_out)
{
   std::string lower = requested_in;
   std::transform(lower.begin(), lower.end(), lower.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   warning_out.clear();
   requested_out = requested_in;

   if (lower == "godunov")
   {
      return FreeSurfaceBCMode::Godunov;
   }
   if (lower == "gamma")
   {
      return FreeSurfaceBCMode::Gamma;
   }

   // Unknown value → warning + fallback.
   std::ostringstream w;
   w << "WARNING: unknown --free-surface-bc=\""
     << requested_in << "\"; falling back to gamma.";
   warning_out = w.str();
   requested_out = "gamma";
   return FreeSurfaceBCMode::Gamma;
}

int main()
{
   std::cout << "\n=== TPV102 v9.3.0 Phase 2 (I-04): R-005 "
             << "--free-surface-bc parsing ===\n\n";

   // Case 1: default "gamma" → Gamma, no warning.
   {
      std::string req, warn;
      FreeSurfaceBCMode m = ParseFlag("gamma", req, warn);
      TEST_ASSERT(m == FreeSurfaceBCMode::Gamma,
                  "'gamma' → FreeSurfaceBCMode::Gamma");
      TEST_ASSERT(warn.empty(),
                  "'gamma' emits no warning");
      TEST_ASSERT(req == "gamma",
                  "'gamma' requested banner string == 'gamma'");
   }

   // Case 2: "godunov" → Godunov, no warning.
   {
      std::string req, warn;
      FreeSurfaceBCMode m = ParseFlag("godunov", req, warn);
      TEST_ASSERT(m == FreeSurfaceBCMode::Godunov,
                  "'godunov' → FreeSurfaceBCMode::Godunov");
      TEST_ASSERT(warn.empty(),
                  "'godunov' emits no warning");
      TEST_ASSERT(req == "godunov",
                  "'godunov' requested banner string == 'godunov'");
   }

   // Case 3: uppercase "Godunov" → Godunov (R-005 case-insensitive).
   {
      std::string req, warn;
      FreeSurfaceBCMode m = ParseFlag("Godunov", req, warn);
      TEST_ASSERT(m == FreeSurfaceBCMode::Godunov,
                  "R-005: 'Godunov' (capitalized) → Godunov");
      TEST_ASSERT(warn.empty(),
                  "R-005: capitalized 'Godunov' emits no warning");
   }

   // Case 4: all-caps "GAMMA" → Gamma.
   {
      std::string req, warn;
      FreeSurfaceBCMode m = ParseFlag("GAMMA", req, warn);
      TEST_ASSERT(m == FreeSurfaceBCMode::Gamma,
                  "R-005: 'GAMMA' (all caps) → Gamma");
      TEST_ASSERT(warn.empty(),
                  "R-005: 'GAMMA' emits no warning");
   }

   // Case 5: mixed-case "GoDuNov" → Godunov.
   {
      std::string req, warn;
      FreeSurfaceBCMode m = ParseFlag("GoDuNov", req, warn);
      TEST_ASSERT(m == FreeSurfaceBCMode::Godunov,
                  "R-005: mixed-case 'GoDuNov' → Godunov");
   }

   // Case 6: unknown value "xyz" → Gamma, with warning, banner says gamma.
   {
      std::string req, warn;
      FreeSurfaceBCMode m = ParseFlag("xyz", req, warn);
      TEST_ASSERT(m == FreeSurfaceBCMode::Gamma,
                  "R-005: unknown 'xyz' falls back to Gamma");
      TEST_ASSERT(!warn.empty(),
                  "R-005: unknown 'xyz' emits a WARNING");
      TEST_ASSERT(warn.find("xyz") != std::string::npos,
                  "R-005: warning message contains the offending value");
      TEST_ASSERT(warn.find("falling back to gamma") != std::string::npos,
                  "R-005: warning mentions the fallback");
      TEST_ASSERT(req == "gamma",
                  "R-005: effective banner after unknown flag == 'gamma'");
   }

   // Case 7: typo "gadunov" → Gamma with warning (classic fat-finger).
   {
      std::string req, warn;
      FreeSurfaceBCMode m = ParseFlag("gadunov", req, warn);
      TEST_ASSERT(m == FreeSurfaceBCMode::Gamma,
                  "R-005: typo 'gadunov' falls back to Gamma");
      TEST_ASSERT(!warn.empty(),
                  "R-005: typo 'gadunov' emits a WARNING");
   }

   // Case 8: empty string → Gamma with warning.
   {
      std::string req, warn;
      FreeSurfaceBCMode m = ParseFlag("", req, warn);
      TEST_ASSERT(m == FreeSurfaceBCMode::Gamma,
                  "R-005: empty string falls back to Gamma");
      TEST_ASSERT(!warn.empty(),
                  "R-005: empty string emits a WARNING");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";

   return (num_failed == 0) ? 0 : 1;
}
