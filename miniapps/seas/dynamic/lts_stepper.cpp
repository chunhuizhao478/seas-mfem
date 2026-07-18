// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// lts_stepper.cpp — the clustered-LTS tick loop (Normative Scheduling of
// PLAN_clustered_lts_ader_2026-07-18.md).
//
// PHASE 2 delivers the pure SEQUENCING here: RunSyncInterval walks a precomputed
// tick table and dispatches Predict/Correct to an ILtsClusterStepper.  The
// concrete stepper (wave-operator per-cluster predictor/corrector) plugs in via
// that interface; this file has zero dependence on the wave operator, so the
// schedule is unit-testable against a mock (Appendix B.4).

#include "lts_stepper.hpp"

namespace mfem
{
namespace seas
{

void RunSyncInterval(const std::vector<LtsTick>& table, ILtsClusterStepper& stepper)
{
   for (const LtsTick& tk : table)
   {
      // Predict opens each due cluster's step at this tick (element-local; any
      // order — no cross-cluster dependence in the predictor).
      for (int c : tk.predict_clusters) { stepper.Predict(c, tk.dt_step[c]); }

      // Correct closes each due cluster's step FINE→COARSE, so a fine cluster's
      // final sub-interval contribution reaches the coarse accumulate buffer
      // before the coarse cluster consumes it in this same tick.  The tick table
      // already sorts correct_clusters ascending (== FINE→COARSE).
      for (int c : tk.correct_clusters) { stepper.Correct(c, tk.dt_step[c]); }
   }
}

} // namespace seas
} // namespace mfem
