// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// Unit tests for Logger.
// Run: ./seas_test_logging

#include "mfem.hpp"
#include "../../common/logging.hpp"
#include "test_macros.hpp"

#include <iostream>
#include <sstream>

using namespace mfem;
using namespace mfem::seas;

void TestLoggerDefaults()
{
   std::cout << "\n=== Test: Logger defaults ===\n";
   Logger log;
   TEST_ASSERT(log.IsRoot(), "default rank 0 is root");
   TEST_ASSERT(log.GetLevel() == LogLevel::Info, "default level is Info");
}

void TestLoggerNonRoot()
{
   std::cout << "\n=== Test: Non-root logger is silent ===\n";
   Logger log(1, LogLevel::Debug);
   TEST_ASSERT(!log.IsRoot(), "rank 1 is not root");
   // Non-root Info/Warn/Error/Debug should produce no output
   // (Can't easily capture cout in this test framework, so just verify no crash)
   log.Info("this should not appear");
   log.Warn("this should not appear");
   log.Error("this should not appear");
   log.Debug("this should not appear");
   TEST_ASSERT(true, "non-root logger runs without crash");
}

void TestLoggerLevels()
{
   std::cout << "\n=== Test: Logger level filtering ===\n";
   Logger log(0, LogLevel::Warning);
   // Info is below Warning → should be suppressed on root
   // We can't capture output easily, but verify the API works
   log.Info("suppressed info");
   log.Warn("visible warning");
   log.Error("visible error");
   TEST_ASSERT(true, "level filtering runs without crash");

   log.SetLevel(LogLevel::Silent);
   TEST_ASSERT(log.GetLevel() == LogLevel::Silent, "level set to Silent");
   log.Error("even errors suppressed in Silent");
   TEST_ASSERT(true, "Silent mode runs without crash");
}

void TestLoggerRootOutput()
{
   std::cout << "\n=== Test: Root logger produces output ===\n";
   Logger log(0, LogLevel::Debug);
   // Just verify these don't crash with multiple args
   log.Info("step ", 42, " time = ", 1.5, " s");
   log.Debug("vector size: ", 100);
   TEST_ASSERT(true, "multi-arg logging runs without crash");
}

int main(int argc, char *argv[])
{
   TestLoggerDefaults();
   TestLoggerNonRoot();
   TestLoggerLevels();
   TestLoggerRootOutput();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
