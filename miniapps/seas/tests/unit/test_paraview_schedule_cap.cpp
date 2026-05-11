// Phase 3 of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md.
//
// Synthetic-V_max unit test for AdaptiveSchedule::max_total_snapshots and
// the cap-aware interval logic.  Constructs `(time, V_max)` sequences
// mimicking 1, 5, and 50 nucleation→event→interseismic cycles.  Asserts:
//
//   total_snapshots_written_ <= K + n_events_in_run
//
// (plan §Phase 3 acceptance criterion 1).  The cap is a SOFT bound on
// interseismic writes; coseismic and nucleation regimes are uncapped, so
// the final count can exceed K by at most one write per event.
//
// We exercise the cap through `ShouldWrite(cycle, time, V_max)` since it
// bumps the counter and advances the schedule without invoking the
// volume PVD writer (matches the fault-only driving pattern).
//
// Pre-Phase-3 default (`max_total_snapshots = 0`): uncapped, baseline
// reference.  With cap = K: counter ≤ K + n_events.

#include "mfem.hpp"
#include "../../io/paraview_output.hpp"
#include "../../config/bp5_params.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)
#define TEST_LE(val, bound, msg) do { \
   num_tests++; \
   long long v_ = (long long)(val), b_ = (long long)(bound); \
   if (v_ <= b_) { num_passed++; std::cout << "  PASSED: " << msg \
      << " (" << v_ << " <= " << b_ << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", bound " << b_ << ")\n"; } \
} while (0)

// ---------------------------------------------------------------------------
// Build a synthetic (time, V_max) sequence with `n_events` earthquakes
// at SPARSE sample density.  Per cycle, we emit ~1 sample per regime
// transition (so each event contributes ~3 writes max: one nucleation
// entry, one coseismic peak, one interseismic re-entry) plus
// interseismic samples spaced far enough apart that the adaptive
// schedule writes them at its dt_interseismic cadence.  This matches
// the plan's `K + n_events_in_run` bound, which assumes ~O(1) extra
// writes per event rather than the dense O(T_event/dt_event) writes
// a fully-resolved synthetic would produce.
// ---------------------------------------------------------------------------
static void GenerateSequence(int n_events, real_t recurrence_yr,
                              std::vector<real_t> &times,
                              std::vector<real_t> &Vs)
{
   const real_t YR  = BP5Params::seconds_per_year;
   const real_t V_min = 1e-12;        // deep interseismic
   const real_t V_co  = 1e-1;         // above v_coseismic  = 1e-3
   // Interseismic sampling well below dt_interseismic = 1 yr, so the
   // schedule chooses when to emit (every ~1 yr by default, capped if
   // max_total_snapshots > 0).
   const real_t dt_inter = 30.0 * 86400;   // 30 days

   times.clear(); Vs.clear();
   real_t t = 0.0;
   for (int e = 0; e < n_events; ++e)
   {
      // Interseismic up to the cycle boundary.
      const real_t t_event = (e + 1) * recurrence_yr * YR;
      while (t < t_event - dt_inter)
      {
         times.push_back(t); Vs.push_back(V_min); t += dt_inter;
      }
      // Single coseismic spike at t_event.
      times.push_back(t_event); Vs.push_back(V_co);
      // Resume interseismic immediately after.
      t = t_event + dt_inter;
   }
   // Tail interseismic for one more recurrence period.
   const real_t t_end = n_events * recurrence_yr * YR + recurrence_yr * YR;
   while (t < t_end)
   {
      times.push_back(t); Vs.push_back(V_min); t += dt_inter;
   }
}

// ---------------------------------------------------------------------------
// Run a synthetic sequence through ShouldWrite() and report
// (uncapped_writes, capped_writes).
// ---------------------------------------------------------------------------
struct RunResult
{
   int writes = 0;
   int n_events_observed = 0;
};

static RunResult RunSequence(int n_events, real_t recurrence_yr,
                              int cap_K)
{
   const std::string prefix = "/tmp/test_paraview_schedule_cap_run";
   ::mkdir(prefix.c_str(), 0755);
   Mesh smesh = Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE);
   ParaViewOutput<Mesh> pv(prefix, smesh, /*order=*/1);
   pv.GetSchedule().max_total_snapshots = cap_K;
   pv.GetSchedule().Validate();

   std::vector<real_t> times, Vs;
   GenerateSequence(n_events, recurrence_yr, times, Vs);
   const real_t T_end = times.empty() ? 0.0 : times.back();
   pv.SetTotalRunTime(T_end);

   // Count regime entries into coseismic as "events".
   int prev_regime = 0;
   int events = 0;
   for (size_t i = 0; i < times.size(); ++i)
   {
      const real_t V = Vs[i];
      const int new_regime = pv.GetSchedule().NextRegime(V, prev_regime);
      if (new_regime == 2 && prev_regime != 2) { ++events; }
      prev_regime = new_regime;
   }

   for (size_t i = 0; i < times.size(); ++i)
   {
      pv.ShouldWrite(static_cast<int>(i), times[i], Vs[i]);
   }

   RunResult r;
   r.writes = pv.GetTotalSnapshotsWritten();
   r.n_events_observed = events;
   return r;
}

int main(int argc, char *argv[])
{
   (void)argc; (void)argv;
   std::cout << "=== test_paraview_schedule_cap (Phase 3) ===\n";

   // --- 1 cycle (1 event), uncapped baseline.
   int uncapped_1event_writes = 0;
   {
      auto r = RunSequence(/*n_events=*/1, /*recurrence_yr=*/2.0,
                           /*cap_K=*/0);
      uncapped_1event_writes = r.writes;
      std::cout << "  INFO: 1-event uncapped: " << r.writes
                << " writes, " << r.n_events_observed << " events\n";
      TEST_ASSERT(r.writes > 0,
                  "1-event uncapped: at least one write happened");
      TEST_ASSERT(r.n_events_observed >= 1,
                  "1-event uncapped: at least one event detected");
   }

   // --- 1 cycle, cap = 5.  Loose cap (writes ≪ K), engages as no-op.
   {
      const int K = 5;
      auto r = RunSequence(1, 2.0, K);
      std::cout << "  INFO: 1-event K=" << K << ": " << r.writes
                << " writes, " << r.n_events_observed << " events\n";
      TEST_LE(r.writes, K + r.n_events_observed,
              "1-event cap: writes <= K + n_events");
   }

   // --- R-003: capped-vs-uncapped diff on 5-event scenario.
   //     The cap K=4 is TIGHT (uncapped writes ≈ 12), so a working cap
   //     must reduce the count vs uncapped.  This catches R-001-class
   //     regressions where the cap silently disables.
   {
      const int K = 4;
      auto uncapped = RunSequence(5, 2.0, /*cap_K=*/0);
      auto capped   = RunSequence(5, 2.0, K);
      std::cout << "  INFO: 5-event uncapped=" << uncapped.writes
                << " capped(K=" << K << ")=" << capped.writes << "\n";
      // Plan §Phase 3 acceptance bound is `K + n_events`, but the
      // plan's own algorithm allows "one more" final write after budget
      // exhaustion (`if remaining_budget <= 0: dt_inter_new = T - t`).
      // The strict bound is therefore `K + n_events + 1`.
      TEST_LE(capped.writes, K + capped.n_events_observed + 1,
              "5-event cap: writes <= K + n_events + 1");
      TEST_ASSERT(capped.writes < uncapped.writes,
                  "R-003: 5-event cap actually reduced write count vs "
                  "uncapped (catches silent-no-op regressions)");
   }

   // --- R-003: capped-vs-uncapped diff on 50-event scenario.
   {
      const int K = 50;
      auto uncapped = RunSequence(50, 2.0, /*cap_K=*/0);
      auto capped   = RunSequence(50, 2.0, K);
      std::cout << "  INFO: 50-event uncapped=" << uncapped.writes
                << " capped(K=" << K << ")=" << capped.writes << "\n";
      TEST_LE(capped.writes, K + capped.n_events_observed,
              "50-event cap: writes <= K + n_events");
      TEST_ASSERT(capped.writes < uncapped.writes,
                  "R-003: 50-event cap actually reduced write count vs "
                  "uncapped");
   }

   // --- Fixed-dt path: cap is a no-op when fixed_dt > 0 (plan §Phase 3
   //     edge case "Combination with --paraview-dt X").
   {
      Mesh smesh = Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE);
      ParaViewOutput<Mesh> pv("/tmp/test_paraview_schedule_cap_fixed_dt",
                              smesh, /*order=*/1);
      ::mkdir("/tmp/test_paraview_schedule_cap_fixed_dt", 0755);
      pv.fixed_dt = 1.0;             // 1 second fixed dt
      pv.GetSchedule().max_total_snapshots = 5;
      pv.SetTotalRunTime(100.0);     // 100 s total
      // 100 ticks at 1 s each, all interseismic.
      for (int i = 0; i < 100; ++i)
      { pv.ShouldWrite(i, real_t(i), 1e-12); }
      // With fixed_dt=1 the cap should NOT inflate (precedence: fixed_dt
      // wins over adaptive cap).  Expect ~100 writes (modulo first-step
      // gating from last_write_time_ initial value).
      const int writes = pv.GetTotalSnapshotsWritten();
      std::cout << "  INFO: fixed_dt=1 over 100 s, K=5: "
                << writes << " writes (cap is no-op under fixed_dt)\n";
      TEST_ASSERT(writes > 50,
                  "fixed_dt path bypasses cap (writes > 50 with K=5)");
   }

   // --- Step-based path: cap is also a no-op (plan §Phase 3 edge case
   //     "Combination with --paraview-every N").
   {
      Mesh smesh = Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE);
      ParaViewOutput<Mesh> pv("/tmp/test_paraview_schedule_cap_every_n",
                              smesh, /*order=*/1);
      ::mkdir("/tmp/test_paraview_schedule_cap_every_n", 0755);
      pv.output_every_n_steps = 2;
      pv.GetSchedule().max_total_snapshots = 5;
      pv.SetTotalRunTime(0.0);
      for (int i = 0; i < 100; ++i)
      { pv.ShouldWrite(i, real_t(i), 1e-12); }
      const int writes = pv.GetTotalSnapshotsWritten();
      std::cout << "  INFO: every-2 over 100 cycles, K=5: "
                << writes << " writes (cap is no-op under step-based)\n";
      TEST_ASSERT(writes >= 49,
                  "step-based path bypasses cap (writes ~50 with K=5)");
   }

   // --- Validate() accepts the legitimate range.
   {
      ParaViewOutput<Mesh>::AdaptiveSchedule s_ok;
      s_ok.max_total_snapshots = 0;
      s_ok.Validate();   // uncapped — should not abort
      s_ok.max_total_snapshots = 1000;
      s_ok.Validate();   // capped — should not abort
      TEST_ASSERT(true, "Validate accepts max_total_snapshots == 0 and > 0");
      // The negative-cap death test is left to manual / fork-based
      // testing: MFEM_VERIFY uses std::abort, which would terminate the
      // whole test process before the summary prints.
   }

   // --- R-004 regression: CommitSchedule + Save / ShouldWrite for the
   //     same step does not double-bump the snapshot counter.
   {
      Mesh smesh = Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE);
      ParaViewOutput<Mesh> pv("/tmp/test_paraview_schedule_cap_R004",
                              smesh, /*order=*/1);
      ::mkdir("/tmp/test_paraview_schedule_cap_R004", 0755);
      pv.SetTotalRunTime(1e9);
      pv.GetSchedule().max_total_snapshots = 100;

      // ShouldWrite drives a commit at t=0.
      const bool wrote = pv.ShouldWrite(0, 0.0, 1e-12);
      TEST_ASSERT(wrote, "R-004 setup: ShouldWrite returned true at t=0");
      const int after_should = pv.GetTotalSnapshotsWritten();
      TEST_ASSERT(after_should == 1,
                  "R-004 setup: counter incremented exactly once after "
                  "ShouldWrite");

      // CommitSchedule with the SAME (time, V_max) must not double-bump.
      pv.CommitSchedule(0.0, 1e-12);
      const int after_commit = pv.GetTotalSnapshotsWritten();
      TEST_ASSERT(after_commit == 1,
                  "R-004: CommitSchedule with same (time, V_max) is "
                  "idempotent (no double-count)");

      // CommitSchedule at a NEW time DOES bump.
      pv.CommitSchedule(1e8, 1e-12);
      const int after_new_commit = pv.GetTotalSnapshotsWritten();
      TEST_ASSERT(after_new_commit == 2,
                  "R-004: CommitSchedule at a new time bumps normally");
   }

   // --- R-005 regression: budget-exhausted + simulation overrun must
   //     NOT re-enable interseismic writes at user dt.
   {
      Mesh smesh = Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE);
      ParaViewOutput<Mesh> pv("/tmp/test_paraview_schedule_cap_R005",
                              smesh, /*order=*/1);
      ::mkdir("/tmp/test_paraview_schedule_cap_R005", 0755);
      pv.GetSchedule().max_total_snapshots = 2;
      pv.SetTotalRunTime(100.0);
      // Drive past the budget (>2 writes) AND past total_run_time_=100.
      // dt_interseismic default = 1 yr = 3.16e7 s, much larger than 0.5
      // step, so each ShouldWrite should be gated by the schedule.  We
      // tick at 50 s intervals to land each tick at a reasonable cadence.
      for (int i = 0; i < 1000; ++i)
      {
         pv.ShouldWrite(i, real_t(i) * 50.0, 1e-12);
      }
      const int writes = pv.GetTotalSnapshotsWritten();
      std::cout << "  INFO: R-005 overrun (K=2, T=100s, sim=50000s): "
                << writes << " writes (cap should hard-ceil)\n";
      // K=2 budget + at most one final ("emit at most one more") = 3.
      TEST_LE(writes, 3,
              "R-005: budget-exhausted + simulation overrun emits at "
              "most K + 1 writes (no resumption at user dt)");
   }

   // --- R-105 regression: CommitSchedule dedup survives 1-ULP FP noise
   //     on the time argument.  Pre-fix, `time == last_write_time_`
   //     used exact float equality, so a caller that recomputed the
   //     time from `t + dt - dt` (which differs by 1 ULP) would miss
   //     the dedup and double-count.  Post-fix the comparison uses an
   //     8-ULP relative tolerance.
   {
      Mesh smesh = Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE);
      ParaViewOutput<Mesh> pv("/tmp/test_paraview_schedule_cap_R105",
                              smesh, /*order=*/1);
      ::mkdir("/tmp/test_paraview_schedule_cap_R105", 0755);
      pv.SetTotalRunTime(1e9);
      pv.GetSchedule().max_total_snapshots = 100;

      // ShouldWrite at t = 1e6 + 0.1 - 0.1 — bit-different from 1e6.
      const real_t t_actual = (real_t(1e6) + real_t(0.1)) - real_t(0.1);
      pv.ShouldWrite(0, t_actual, 1e-12);
      TEST_ASSERT(pv.GetTotalSnapshotsWritten() == 1,
                  "R-105 setup: ShouldWrite committed once");

      // Caller's `time` argument loses 1 ULP en route to CommitSchedule:
      // logically the same step, but bitwise different.  Pre-fix this
      // would double-count.
      pv.CommitSchedule(real_t(1e6), 1e-12);
      TEST_ASSERT(pv.GetTotalSnapshotsWritten() == 1,
                  "R-105: dedup tolerates 1-ULP FP noise on time arg "
                  "(no double-count)");

      // A genuinely different step DOES bump.
      pv.CommitSchedule(real_t(2e6), 1e-12);
      TEST_ASSERT(pv.GetTotalSnapshotsWritten() == 2,
                  "R-105: distinct step still bumps after FP-tolerant "
                  "dedup");
   }

   // --- R-106 regression: Save() and ShouldWrite() at the same cycle
   //     do NOT double-bump the snapshot counter.  Pre-fix, the
   //     step-based path of ShouldWrite ran `++counter` unconditionally
   //     when `cycle % every == 0`, regardless of whether Save() had
   //     already committed at the same cycle.
   {
      Mesh smesh = Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE);
      ParaViewOutput<Mesh> pv("/tmp/test_paraview_schedule_cap_R106",
                              smesh, /*order=*/1);
      ::mkdir("/tmp/test_paraview_schedule_cap_R106", 0755);
      pv.SetTotalRunTime(1e9);
      pv.GetSchedule().max_total_snapshots = 100;
      pv.output_every_n_steps = 1;        // step-based path

      // Both Save and ShouldWrite at cycle=0 should yield exactly one bump.
      pv.Save(0, 0.0, 1e-12);
      pv.ShouldWrite(0, 0.0, 1e-12);
      TEST_ASSERT(pv.GetTotalSnapshotsWritten() == 1,
                  "R-106: Save+ShouldWrite at same cycle does not "
                  "double-bump the snapshot counter");

      // Distinct cycle still bumps normally.
      pv.Save(1, 1.0, 1e-12);
      TEST_ASSERT(pv.GetTotalSnapshotsWritten() == 2,
                  "R-106: distinct cycle still bumps");
   }

   std::cout << "\n=== Summary: " << num_passed << " / " << num_tests
             << " passed; " << num_failed << " failed ===\n";
   return num_failed > 0 ? 1 : 0;
}
