// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_framework.hpp — the minimal assert/report macros of the seas unit tests
// (miniapps/seas/tests/unit/*.cpp), shared by the pff tests.

#ifndef MFEM_PFF_TEST_FRAMEWORK_HPP
#define MFEM_PFF_TEST_FRAMEWORK_HPP

#include <cmath>
#include <iostream>

static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "FAILED: " << (message) << " (line " << __LINE__ << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << (message) << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_NEAR(value, expected, tol, message) \
   do { \
      num_tests++; \
      const double _val = (value); \
      const double _exp = (expected); \
      const double _tol = (tol); \
      if (!(std::abs(_val - _exp) <= _tol)) { \
         std::cerr << "FAILED: " << (message) << ": got " << _val << ", expected " \
                   << _exp << " +- " << _tol << " (line " << __LINE__ << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << (message) << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_REL(value, expected, rtol, message) \
   TEST_NEAR(value, expected, (rtol) * std::abs((double)(expected)), message)

inline int TestSummary(const char *name)
{
   std::cout << "\n" << name << ": " << num_passed << "/" << num_tests
             << " passed, " << num_failed << " failed\n";
   return num_failed == 0 ? 0 : 1;
}

#endif
