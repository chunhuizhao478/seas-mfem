// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for the TPV104 driver (§4.10 Step 9 gates T_TPV104_SMOKE_1..3).
//
// Scope:
//   - T_TPV104_SMOKE_1 / SMOKE_2 — 100-step mesh-coupled runs are
//     DEFERRED; CLAUDE.md memory `feedback_no_local_reproducer` blocks
//     local mesh-coupled TPV104 runs (no 1000m mesh).  Tracked by the
//     Phase-2 acceptance matrix under Frontera.
//   - T_TPV104_SMOKE_3 — driver banner check.  This test binary launches
//     the `seas_tpv104_driver` binary via popen and scrapes stdout for
//     the four required banner lines.  If the driver binary is not
//     present at the expected path the test skips with a non-counting
//     message.

#include "test_macros.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;

// --------------------------------------------------------------------------
// popen-based driver-binary runner.  Returns empty string on failure.
// --------------------------------------------------------------------------
static std::string RunDriver(const std::string &binary_path,
                             const std::string &extra_args)
{
   const std::string cmd = binary_path + " " + extra_args + " 2>&1";
   FILE *fp = popen(cmd.c_str(), "r");
   if (!fp) { return ""; }
   std::ostringstream oss;
   char buf[256];
   while (std::fgets(buf, sizeof(buf), fp) != nullptr)
   {
      oss << buf;
   }
   pclose(fp);
   return oss.str();
}

static std::string LocateDriver()
{
   const std::vector<std::string> candidates = {
      "./seas_tpv104_driver",
      "../seas_tpv104_driver",
      "../../seas_tpv104_driver",
      "/tmp/seas_tpv104_driver"
   };
   for (const auto &c : candidates)
   {
      std::ifstream f(c);
      if (f.good()) { return c; }
   }
   return "";
}

// --------------------------------------------------------------------------
// T_TPV104_SMOKE_3 — banner content on default flags (R7-001 honest).
// --------------------------------------------------------------------------
void TestBannerDefaults()
{
   std::cout << "\n[T_TPV104_SMOKE_3] driver banner on default flags\n";

   const std::string binary = LocateDriver();
   if (binary.empty())
   {
      std::cout << "  SKIPPED: seas_tpv104_driver binary not found. "
                << "Build it first (e.g. via `make seas_tpv104_driver` "
                << "or the Step-12 target).\n";
      return;
   }
   std::cout << "  using driver binary: " << binary << "\n";

   // `--dry-run` so the driver exits quickly without expecting a mesh.
   const std::string out = RunDriver(binary, "--dry-run");
   TEST_ASSERT(!out.empty(), "driver produced non-empty output");

   // Four mandatory banner lines, R7-001-honest: the banner now
   // describes the runtime dispatch truthfully — Brent hard-coded via
   // EvaluateADERTotal, one-shot wave.AdvanceADER (sub-step iterator
   // not wired), slip-SRW ψ-space with macro-step analytic cadence.
   const std::vector<std::pair<std::string, std::string>> must_have = {
      {"Time integrator: ADER-O2 (one-shot via wave.AdvanceADER)",
       "ADER-O2 + one-shot disclosure"},
      {"Fault iterator: one-shot (sub-step iterator NOT wired",
       "sub-step iterator NOT wired disclosure (R7-001)"},
      {"Friction solver: Brent (hard-coded via EvaluateADERTotal",
       "Brent hard-coded disclosure (R7-001)"},
      {"Friction law: slip-SRW (ψ-space, macro-step analytic",
       "slip-SRW macro-step cadence disclosure (R7-007)"},
   };
   for (const auto &kv : must_have)
   {
      const bool found = out.find(kv.first) != std::string::npos;
      TEST_ASSERT(found,
                  ("banner contains '" + kv.first + "' ("
                   + kv.second + ")").c_str());
      if (!found)
      {
         std::cerr << "  banner text was:\n" << out << "\n";
      }
   }

   // The banner must NOT advertise any misleading "sub-step iterator"
   // or "stable-asinh Newton" as the runtime dispatch on default flags
   // (those would reintroduce R7-001's banner-vs-code divergence).
   const std::vector<std::pair<std::string, std::string>> forbidden = {
      {"Fault iterator: sub-step\n",
       "(no bare 'sub-step' as dispatched iterator)"},
      {"Friction solver: Newton-Raphson (stable-asinh",
       "(no stable-asinh Newton as dispatched solver — R7-001)"},
      {"Friction solver: Newton-Raphson (legacy",
       "(no legacy Newton as dispatched solver)"},
      {"Friction solver: Hybrid",
       "(no Hybrid as dispatched solver)"},
      {"Friction law: aging",
       "(no aging law as default)"},
   };
   for (const auto &kv : forbidden)
   {
      const bool has = out.find(kv.first) != std::string::npos;
      TEST_ASSERT(!has,
                  ("banner does NOT contain '" + kv.first + "' "
                   + kv.second).c_str());
   }
}

// --------------------------------------------------------------------------
// T_TPV104_SMOKE_3b — CLI echo line on non-default flags.
// Under R7-001 option (b), --friction-solver and --fault-iterator do
// NOT change the dispatched solver/iterator.  They show up only in the
// "CLI parsed (banner-only, not dispatched)" echo line, and in
// --ader-order's ADER-O{N} disclosure.
// --------------------------------------------------------------------------
void TestBannerNonDefault()
{
   std::cout << "\n[T_TPV104_SMOKE_3b] driver banner on non-default flags\n";

   const std::string binary = LocateDriver();
   if (binary.empty())
   {
      std::cout << "  SKIPPED (no driver binary).\n";
      return;
   }

   const std::string out = RunDriver(
      binary,
      "--dry-run --friction-solver newton-legacy --fault-iterator oneshot "
      "--ader-order 5");

   TEST_ASSERT(
      out.find("Time integrator: ADER-O5 (one-shot via wave.AdvanceADER)")
      != std::string::npos,
      "banner reflects --ader-order 5");

   // Dispatch lines are always the same regardless of CLI values
   // (R7-001 disclosure invariant): Brent + one-shot.
   TEST_ASSERT(
      out.find("Fault iterator: one-shot (sub-step iterator NOT wired")
      != std::string::npos,
      "fault iterator disclosure is constant on non-default flags");
   TEST_ASSERT(
      out.find("Friction solver: Brent (hard-coded via EvaluateADERTotal")
      != std::string::npos,
      "friction solver disclosure is constant on non-default flags");

   // CLI echo line captures what the user requested.
   TEST_ASSERT(
      out.find("friction_solver=newton-legacy") != std::string::npos,
      "CLI echo shows --friction-solver=newton-legacy");
   TEST_ASSERT(
      out.find("fault_iterator=oneshot") != std::string::npos,
      "CLI echo shows --fault-iterator=oneshot");

   // Explicit 'brent' does NOT change dispatch (it already was brent)
   // but also must not trip the CLI validator.
   const std::string out_brent = RunDriver(
      binary, "--dry-run --friction-solver brent");
   TEST_ASSERT(
      out_brent.find("Friction solver: Brent (hard-coded via EvaluateADERTotal")
      != std::string::npos,
      "--friction-solver brent → still shows hard-coded Brent disclosure");
}

// --------------------------------------------------------------------------
// T_TPV104_SMOKE_3c — R7-004 tri-consistency: banner ↔ [dispatch] lines.
// On `--verify-dispatch`, the driver emits machine-readable dispatch
// tags.  Assert that on ALL CLI values (default, newton-stable,
// newton-legacy, brent, substep, oneshot) the dispatch tag is
// "brent" / "oneshot" / "slip-srw" — the runtime is invariant.
// --------------------------------------------------------------------------
void TestDispatchMatchesBanner()
{
   std::cout << "\n[T_TPV104_SMOKE_3c] R7-004 banner vs dispatch "
             << "tri-consistency check\n";

   const std::string binary = LocateDriver();
   if (binary.empty())
   {
      std::cout << "  SKIPPED (no driver binary).\n";
      return;
   }

   const std::vector<std::string> cli_variants = {
      "--dry-run --verify-dispatch",
      "--dry-run --verify-dispatch --friction-solver newton-stable "
      "--fault-iterator substep",
      "--dry-run --verify-dispatch --friction-solver newton-legacy "
      "--fault-iterator oneshot",
      "--dry-run --verify-dispatch --friction-solver brent",
      "--dry-run --verify-dispatch --friction-solver newton",
   };

   for (const auto &cli : cli_variants)
   {
      const std::string out = RunDriver(binary, cli);
      TEST_ASSERT(
         out.find("[dispatch] rank=0 friction_solver_actual=brent")
         != std::string::npos,
         ("dispatch shows friction_solver=brent under CLI: " + cli).c_str());
      TEST_ASSERT(
         out.find("[dispatch] rank=0 fault_iterator_actual=oneshot")
         != std::string::npos,
         ("dispatch shows fault_iterator=oneshot under CLI: " + cli).c_str());
      TEST_ASSERT(
         out.find("[dispatch] rank=0 friction_law_actual=slip-srw")
         != std::string::npos,
         ("dispatch shows friction_law=slip-srw under CLI: " + cli).c_str());
      // Banner matches the dispatch:
      TEST_ASSERT(
         out.find("Friction solver: Brent (hard-coded via EvaluateADERTotal")
         != std::string::npos,
         ("banner shows Brent hard-coded under CLI: " + cli).c_str());
      TEST_ASSERT(
         out.find("Fault iterator: one-shot (sub-step iterator NOT wired")
         != std::string::npos,
         ("banner shows one-shot under CLI: " + cli).c_str());
   }
}

// --------------------------------------------------------------------------
// T_TPV104_SMOKE_3d — R7-005: unknown --friction-solver aborts.
// --------------------------------------------------------------------------
void TestUnknownSolverAborts()
{
   std::cout << "\n[T_TPV104_SMOKE_3d] R7-005 unknown solver aborts\n";

   const std::string binary = LocateDriver();
   if (binary.empty())
   {
      std::cout << "  SKIPPED (no driver binary).\n";
      return;
   }

   // Unknown value must abort loudly, not silently default.  We cannot
   // portably verify the exit code through popen, but MFEM_ABORT prints
   // a recognisable error string to stderr (captured by `2>&1`).
   const std::string out = RunDriver(
      binary, "--dry-run --friction-solver nonesuch-typo");
   const bool aborted =
      out.find("--friction-solver: unknown value") != std::string::npos
      || out.find("MFEM abort") != std::string::npos
      || out.find("MFEM_ABORT") != std::string::npos;
   TEST_ASSERT(aborted,
               "typo in --friction-solver aborts loudly (R7-005)");
   // And the dry-run OK line must NOT appear (abort happens before).
   TEST_ASSERT(out.find("[dry-run] OK.") == std::string::npos,
               "unknown-solver abort prevents [dry-run] OK.");
}

// --------------------------------------------------------------------------
// Dry-run no-NaN sanity — the driver's built-in `--dry-run` loop must
// complete without reporting NaN.  This is a compact stand-in for
// SMOKE_1 (which requires a real mesh).
// --------------------------------------------------------------------------
void TestDryRunNoNaN()
{
   std::cout << "\n[T_TPV104_SMOKE_1_dryrun] dry-run no-NaN sanity\n";

   const std::string binary = LocateDriver();
   if (binary.empty())
   {
      std::cout << "  SKIPPED (no driver binary).\n";
      return;
   }
   const std::string out = RunDriver(binary, "--dry-run");
   TEST_ASSERT(out.find("[dry-run] OK.") != std::string::npos,
               "dry-run reports OK");
   TEST_ASSERT(out.find("NaN") == std::string::npos
               && out.find("Inf") == std::string::npos,
               "dry-run output contains no NaN / Inf report");
}

// --------------------------------------------------------------------------
// T_TPV104_SMOKE_1 / SMOKE_2 — deferred placeholders.
// --------------------------------------------------------------------------
void DeclareDeferredGates()
{
   std::cout << "\n[T_TPV104_SMOKE_1] mesh-coupled 100-step no-NaN — "
             << "DEFERRED to Frontera (feedback_no_local_reproducer).\n";
   std::cout << "\n[T_TPV104_SMOKE_2] nucleation-triggered rupture "
             << "slip_rate > 1e-3 — DEFERRED to Frontera.\n";
}

int main(int argc, char *argv[])
{
   TestBannerDefaults();
   TestBannerNonDefault();
   TestDispatchMatchesBanner();
   TestUnknownSolverAborts();
   TestDryRunNoNaN();
   DeclareDeferredGates();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
