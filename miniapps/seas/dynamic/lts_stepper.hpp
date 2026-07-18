// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// lts_stepper.hpp — clustered-LTS tick schedule (Appendix A.3 of
// PLAN_clustered_lts_ader_2026-07-18.md).
//
// PHASE 1 delivers only the TICK-TABLE GENERATOR here (pure arithmetic,
// identical on every rank — the "NORMATIVE: Scheduling" section of the plan).
// The actual multi-cluster stepping loop (lts_stepper.cpp, AdvanceADERCluster
// driving) is a PHASE 2 deliverable and is intentionally absent.
//
// Scheduling contract (plan "NORMATIVE: Scheduling", P-001/P-003/P-004):
//   dt_base           : cluster c steps at dt_c = dt_base * 2^c.
//   T_actual          : length of THIS sync interval (< the nominal
//                       dt_base*2^(Nc-1) only on the final interval before
//                       tfinal, where truncation lives).
//   ticks_per_sync    : ceil(T_actual / dt_base)  (one dt_base step = one tick).
//   predict_due(c,t)  : t % 2^c == 0                 (opens the step at t).
//   correct_due(c,t)  : (t+1) % 2^c == 0  OR  t+1 == ticks_per_sync
//                       (closes the step ending at t+1; the OR handles the
//                        ragged final step).  correct order is FINE -> COARSE.
//   dt_step(c,t)      : min(dt_c, T_actual - dt_base*(t - t%2^c))
//                       (the current-step length; < dt_c only when truncated).

#ifndef MFEM_SEAS_LTS_STEPPER_HPP
#define MFEM_SEAS_LTS_STEPPER_HPP

#include "mfem.hpp"   // mfem::real_t, MFEM_VERIFY

#include <algorithm>
#include <cmath>
#include <vector>

namespace mfem
{
namespace seas
{

/// Global (all-rank-consistent) metadata that gates the per-tick collective
/// count.  Built once in Phase 1 (broadcast per-cluster global counts).
struct LtsGlobalMeta
{
   std::vector<long long> global_elems;        ///< per cluster id
   std::vector<long long> global_fault_faces;  ///< per cluster id
   int  num_state = 9;                          ///< corrector exchanges / correcting cluster
   /// D-2 default: fault faces are rank-interior, so the fault predictor-substep
   /// exchange is dropped GLOBALLY.  Set false only for the (Phase-5) cross-fault
   /// experiment.  When false, a predicting cluster with fault faces contributes
   /// ader_order predictor exchanges.
   bool drop_fault_predictor_exchange = true;
};

/// One tick of a sync interval (Appendix A.3).
struct LtsTick
{
   std::vector<int>    predict_clusters;  ///< due this tick (any order)
   std::vector<int>    correct_clusters;  ///< due this tick, sorted FINE -> COARSE
   std::vector<mfem::real_t> dt_step;     ///< [num_clusters]; current-step length per cluster id
   int                 n_collectives = 0; ///< asserted by the Phase-4 per-tick debug counter
};

namespace detail
{
// ticks_per_sync = ceil(T_actual / dt_base) with a small relative tolerance so
// an exact dyadic interval does NOT gain a spurious 1-ulp final tick, and a
// sub-1e-10 residual step is merged away (plan residual-step-merge rule).
inline long long ticks_per_sync(mfem::real_t T_actual, mfem::real_t dt_base)
{
   const double q = static_cast<double>(T_actual) / static_cast<double>(dt_base);
   long long n = static_cast<long long>(std::floor(q));
   if (q - static_cast<double>(n) > 1e-9) { ++n; }   // ceil, tolerant of FP noise
   return std::max<long long>(n, 1);
}

// 2^c as a real_t (exact for c up to the mantissa width; c <= 31 in practice).
inline mfem::real_t pow2(int c) { return std::ldexp(mfem::real_t(1), c); }
} // namespace detail

/// Build the tick table for ONE sync interval.  Pure arithmetic; identical on
/// every rank (this is what preserves the matched-collective contract).
inline std::vector<LtsTick> BuildTickTable(int num_clusters,
                                           mfem::real_t dt_base,
                                           mfem::real_t T_actual,
                                           int ader_order,
                                           const LtsGlobalMeta& meta)
{
   MFEM_VERIFY(num_clusters >= 1, "BuildTickTable: num_clusters must be >= 1.");
   // REVIEW L-3: the step period is 1LL<<c with c up to num_clusters-1; guard
   // against signed-shift UB.  (Nc is capped at max_clusters=32 upstream, but
   // this makes the tick-loop's own precondition explicit.)
   MFEM_VERIFY(num_clusters <= 62,
               "BuildTickTable: num_clusters " << num_clusters
               << " > 62 would overflow the 1LL<<c step period.");
   MFEM_VERIFY(dt_base > mfem::real_t(0), "BuildTickTable: dt_base must be > 0.");
   MFEM_VERIFY(T_actual > mfem::real_t(0), "BuildTickTable: T_actual must be > 0.");
   MFEM_VERIFY(static_cast<int>(meta.global_fault_faces.size()) >= num_clusters,
               "BuildTickTable: meta.global_fault_faces shorter than num_clusters.");

   const long long nticks = detail::ticks_per_sync(T_actual, dt_base);
   std::vector<LtsTick> table(static_cast<std::size_t>(nticks));

   for (long long t = 0; t < nticks; ++t)
   {
      LtsTick& tk = table[static_cast<std::size_t>(t)];
      tk.dt_step.assign(num_clusters, mfem::real_t(0));

      for (int c = 0; c < num_clusters; ++c)
      {
         const long long period = 1LL << c;              // 2^c ticks per step of cluster c
         const long long phase  = t % period;            // ticks since this step started
         const long long start  = t - phase;             // tick at which the current step began
         const mfem::real_t dt_c = dt_base * detail::pow2(c);
         const mfem::real_t rem  = T_actual - dt_base * static_cast<mfem::real_t>(start);
         tk.dt_step[c] = std::min(dt_c, rem);             // truncated only on the last step

         const bool predict_due = (phase == 0);
         const bool correct_due = ((t + 1) % period == 0) || (t + 1 == nticks);
         if (predict_due) { tk.predict_clusters.push_back(c); }
         if (correct_due) { tk.correct_clusters.push_back(c); }
      }

      // Corrects run FINE -> COARSE (ascending cluster id): the fine side's
      // final sub-interval scatter must reach the accumulate buffer before the
      // coarse side consumes it in the same tick.
      std::sort(tk.correct_clusters.begin(), tk.correct_clusters.end());

      // Per-tick collective count (plan Matched-collectives constraint, refined
      // for the split predict/correct predicates): num_state exchanges per
      // correcting cluster, plus ader_order predictor exchanges per predicting
      // cluster that has GLOBAL fault faces (dropped entirely under D-2).
      int nx = meta.num_state * static_cast<int>(tk.correct_clusters.size());
      if (!meta.drop_fault_predictor_exchange)
      {
         for (int c : tk.predict_clusters)
         {
            if (meta.global_fault_faces[c] > 0) { nx += ader_order; }
         }
      }
      tk.n_collectives = nx;
   }
   return table;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LTS_STEPPER_HPP
