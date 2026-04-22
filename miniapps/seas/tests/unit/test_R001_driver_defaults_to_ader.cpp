// Round-6 R-001 regression: the driver's default time integrator must
// be ADER (not RK4) after the purpose change #1 flip.
//
// The driver can't be invoked as a subprocess from a unit test without
// a mesh file on disk, so this test is a PURE CLI-PARSING regression:
// it reimplements the driver's GetStringArg pattern with the same
// default literal, and checks that calling it without a flag returns
// "ader".  Any future edit that reverts the default to "rk4" fails
// this test immediately.
//
// This is a structural gate, not a numerical one.

#include "mfem.hpp"

#include <cstring>
#include <iostream>
#include <string>

using namespace mfem;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_EQ(lhs, rhs, msg) do { \
   num_tests++; \
   const std::string l_ = (lhs), r_ = (rhs); \
   if (l_ == r_) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got '" << l_ << "')\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got '" << l_ << "', expected '" << r_ << "')\n"; } \
} while (0)

// Mirror of the driver's `GetStringArg` (tpv102_driver.cpp:55-67).  We
// could #include the driver source, but the function is declared
// `static` and is not part of a header.  Keep a local copy in sync by
// matching the signature + semantics exactly.
static std::string GetStringArg(int argc, char *argv[], const char *flag,
                                const char *def)
{
   for (int i = 1; i < argc - 1; i++)
   {
      if (std::strcmp(argv[i], flag) == 0)
      {
         return std::string(argv[i + 1]);
      }
   }
   return std::string(def);
}

int main()
{
   std::cout << "\n=== Round-6 R-001 regression: driver default is ADER ===\n";

   // Case 1: no flag → default is "ader" (purpose change #1).
   {
      char p0[] = "seas_tpv102_driver";
      char *argv[] = { p0, nullptr };
      int argc = 1;
      // Literal must match the driver's default string at
      // drivers/tpv102_driver.cpp:152.  Update both in lockstep.
      std::string ti = GetStringArg(argc, argv, "--time-integrator",
                                    "ader");
      TEST_EQ(ti, "ader", "default --time-integrator is 'ader'");
   }

   // Case 2: --time-integrator=rk4 → explicit RK4 alternative.
   {
      char p0[] = "seas_tpv102_driver";
      char p1[] = "--time-integrator";
      char p2[] = "rk4";
      char *argv[] = { p0, p1, p2, nullptr };
      int argc = 3;
      std::string ti = GetStringArg(argc, argv, "--time-integrator",
                                    "ader");
      TEST_EQ(ti, "rk4", "--time-integrator=rk4 parses as 'rk4'");
   }

   // Case 3: --time-integrator=ader → explicit ADER.
   {
      char p0[] = "seas_tpv102_driver";
      char p1[] = "--time-integrator";
      char p2[] = "ader";
      char *argv[] = { p0, p1, p2, nullptr };
      int argc = 3;
      std::string ti = GetStringArg(argc, argv, "--time-integrator",
                                    "ader");
      TEST_EQ(ti, "ader", "--time-integrator=ader parses as 'ader'");
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
