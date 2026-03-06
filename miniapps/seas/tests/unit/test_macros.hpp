// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Shared test macros for SEAS miniapp unit tests.
// Include this header instead of duplicating the macros in each test file.

#ifndef MFEM_SEAS_TEST_MACROS_HPP
#define MFEM_SEAS_TEST_MACROS_HPP

#include "mfem.hpp"
#include <iostream>
#include <cmath>

using mfem::real_t;

static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
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
      real_t _tol = (tol); \
      if (std::abs(_val - _exp) > _tol) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         std::cerr << "  Expected: " << _exp << ", Got: " << _val \
                   << ", Diff: " << std::abs(_val - _exp) << "\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_REL_NEAR(value, expected, rel_tol, message) \
   do { \
      num_tests++; \
      real_t _val = (value); \
      real_t _exp = (expected); \
      real_t _tol = (rel_tol) * std::abs(_exp); \
      if (_exp == 0.0) { _tol = (rel_tol); } \
      if (std::abs(_val - _exp) > _tol) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
         std::cerr << "  Expected: " << _exp << ", Got: " << _val \
                   << ", RelErr: " << ((_exp != 0.0) ? \
                      std::abs(_val - _exp) / std::abs(_exp) : std::abs(_val)) \
                   << "\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_PRINT_RESULTS() \
   do { \
      std::cout << "\n========================================\n"; \
      std::cout << "  Results: " << num_passed << " passed, " \
                << num_failed << " failed out of " << num_tests << " tests\n"; \
      std::cout << "========================================\n"; \
   } while (0)

#endif // MFEM_SEAS_TEST_MACROS_HPP
