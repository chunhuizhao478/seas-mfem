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
   /// LTS Phase 4b (flux-PREMULTIPLIED forecast): when true, the diff-1 rank-seam
   /// coarse forecast is exchanged at the FINE cluster's CORRECT — one batched
   /// (all-NUM_STATE) collective per correcting cluster c with a coarser neighbour
   /// (c < num_clusters-1).  The coarse side integrates its retained D(k) into the
   /// forecast locally and sends the PROJECTED block; predict is retain-only.  (In
   /// 4a/4b-batched this instead exchanged the raw D(k) at each predicting c>=1 —
   /// ader_order collectives each.)  Set by the driver/test for np>1 multi-cluster
   /// runs; false keeps Stage-1 (diff-0) + single-rank counts unchanged.
   bool exchange_bulk_provider_dk = false;
   /// LTS Track-A A1 (per-tick exchange merge): when true, the seam I-exchange and
   /// the seam forecast exchange are each fired ONCE PER TICK for ALL correcting
   /// clusters together, instead of once per correcting cluster.  At Nc=6 that is
   /// 2 rounds/tick = 64/sync, down from 125.  MUST be set identically to the
   /// stepper's SetMergeTickExchanges(): this field is what the tick table predicts
   /// as `n_collectives`, and the driver ABORTS on a mismatch (P-007 matched
   /// collectives).  Does NOT affect the fault predictor-substep term below, which
   /// is on the predict path and untouched by A1.
   bool merge_tick_seam_exchanges = false;
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

      // Per-tick collective count (plan Matched-collectives constraint).
      // LTS Phase 4b: the seam I-exchange is now ONE batched (all-NUM_STATE)
      // collective per correcting cluster (was num_state per-component in 4a), so
      // the I term is 1*|correct|.  Plus ader_order predictor exchanges per
      // predicting cluster with GLOBAL fault faces (dropped under D-2).
      // A1: the seam I-exchange is ONE batched collective per correcting cluster
      // (4b), or -- under the per-tick merge -- ONE for the whole tick regardless
      // of how many clusters correct.
      int nx = meta.merge_tick_seam_exchanges
               ? (tk.correct_clusters.empty() ? 0 : 1)
               : static_cast<int>(tk.correct_clusters.size());
      if (!meta.drop_fault_predictor_exchange)
      {
         for (int c : tk.predict_clusters)
         {
            if (meta.global_fault_faces[c] > 0) { nx += ader_order; }
         }
      }
      // LTS Phase 4b (flux-PREMULTIPLIED forecast): the diff-1 coarse-forecast
      // exchange moved from PREDICT to CORRECT.  Instead of each predicting cluster
      // c>=1 exchanging its raw D(k) (ader_order batched levels), the COARSE side now
      // INTEGRATES its retained D(k) into a forecast and exchanges THAT once, at the
      // FINE cluster's correct — so every correcting cluster c that has a coarser
      // neighbour (c < num_clusters-1) fires exactly ONE batched (all-NUM_STATE)
      // collective.  Predict is now retain-only (no collective).  c==num_clusters-1
      // (coarsest, no coarser neighbour) is skipped — a rank-uniform gate, matched.
      // (ader_order is unused by this term now; retained in the signature for the
      // fault-predictor term above.)
      if (meta.exchange_bulk_provider_dk)
      {
         if (meta.merge_tick_seam_exchanges)
         {
            // A1: one merged forecast round iff ANY correcting cluster has a
            // coarser neighbour (matches PrepareClusterSeamExchangeTick, which
            // fires the forecast round iff any forecast_per_cluster is set).
            bool any_forecast = false;
            for (int c : tk.correct_clusters)
            { if (c < num_clusters - 1) { any_forecast = true; break; } }
            if (any_forecast) { nx += 1; }
         }
         else
         {
            for (int c : tk.correct_clusters)
            {
               if (c < num_clusters - 1) { nx += 1; }
            }
         }
      }
      tk.n_collectives = nx;
   }
   return table;
}

// ---------------------------------------------------------------------------
// Per-cluster storage (Appendix A.5 / A.6).  PHASE 2: the data structures the
// wave-operator predictor/corrector fill/consume; the tick loop below is what
// sequences those calls.  All indices are DENSE SLOTS from the Phase-1 layout
// (LtsLayout::provider_slot_of_elem / buffer_slot_of_elem).
// ---------------------------------------------------------------------------

/// Raw D(k) derivative stacks retained for PROVIDER (coarse-side) elements so a
/// finer neighbour can time-integrate them (Appendix A.5).  Storage is
/// `[slot][k][block]`, block = NUM_STATE * ndof_per_el.  One epoch counter per
/// slot, bumped every time a predict refreshes that slot; a consumer asserts it
/// reads the CURRENT epoch (GAP-A3 — no stale D(k) reads).
struct LtsDkStore
{
   int n_slots = 0, order = 0, block = 0;
   std::vector<mfem::real_t> data;      ///< n_slots * order * block, RAW/unscaled
   std::vector<long long>    epoch;     ///< n_slots

   void Resize(int n_slots_, int order_, int block_)
   {
      n_slots = n_slots_; order = order_; block = block_;
      data.assign(static_cast<std::size_t>(n_slots) * order * block, mfem::real_t(0));
      epoch.assign(static_cast<std::size_t>(n_slots), 0);
   }
   std::size_t Offset(int slot, int k) const
   { return (static_cast<std::size_t>(slot) * order + k) * block; }
   /// Pointer to slot's raw stack (order*block contiguous, D(0) first).
   mfem::real_t*       Stack(int slot)       { return data.data() + Offset(slot, 0); }
   const mfem::real_t* Stack(int slot) const { return data.data() + Offset(slot, 0); }
   void      BumpEpoch(int slot) { ++epoch[static_cast<std::size_t>(slot)]; }
   long long Epoch(int slot) const { return epoch[static_cast<std::size_t>(slot)]; }
   std::size_t BytesPerRank() const { return data.size() * sizeof(mfem::real_t); }
};

/// Pre-M^-1 flux-contribution accumulate buffers, one per consumer-owning
/// (coarse) element (Appendix A.6 / the Buffer design).  Storage `[slot][block]`.
/// Fine sub-corrects ADD into a slot (++fill); the coarse correct CONSUMES it
/// (reads then zeros, fill→0).  Invariant: every buffer is exactly zero at each
/// sync point (AllZero()).
struct LtsAccumulateBuffers
{
   int n_slots = 0, block = 0;
   std::vector<mfem::real_t> data;   ///< n_slots * block, pre-M^-1 residual units
   std::vector<int>          fill;   ///< n_slots; #fine sub-steps accumulated

   void Resize(int n_slots_, int block_)
   {
      n_slots = n_slots_; block = block_;
      data.assign(static_cast<std::size_t>(n_slots) * block, mfem::real_t(0));
      fill.assign(static_cast<std::size_t>(n_slots), 0);
   }
   mfem::real_t*       Buf(int slot)       { return data.data() + static_cast<std::size_t>(slot) * block; }
   const mfem::real_t* Buf(int slot) const { return data.data() + static_cast<std::size_t>(slot) * block; }

   /// Add a length-`block` contribution into `slot` and count one fine sub-step.
   void AddInto(int slot, const mfem::real_t* contrib)
   {
      mfem::real_t* b = Buf(slot);
      for (int j = 0; j < block; ++j) { b[j] += contrib[j]; }
      ++fill[static_cast<std::size_t>(slot)];
   }
   /// Copy `slot` into `out` (length `block`), zero the slot, reset its fill, and
   /// RETURN the fill count that was consumed (for the fill-count invariant).
   int ConsumeZero(int slot, mfem::real_t* out)
   {
      mfem::real_t* b = Buf(slot);
      for (int j = 0; j < block; ++j) { out[j] = b[j]; b[j] = mfem::real_t(0); }
      const int f = fill[static_cast<std::size_t>(slot)];
      fill[static_cast<std::size_t>(slot)] = 0;
      return f;
   }
   bool AllZero() const
   {
      for (mfem::real_t v : data) { if (v != mfem::real_t(0)) { return false; } }
      return true;
   }
};

/// Abstract per-cluster stepper.  Phase 2 implements this over the wave operator
/// (predict = per-cluster CK predictor; correct = per-cluster corrector); tests
/// implement a mock.  `RunSyncInterval` is the pure, rank-identical tick loop.
struct ILtsClusterStepper
{
   virtual ~ILtsClusterStepper() = default;
   /// Called by RunSyncInterval at the START of each tick (default no-op).  The
   /// multi-rate bulk stepper uses `tick` to derive the consumer-face
   /// sub-interval [a,b] = [t(tick)-t_origin, ...] from the closed-form schedule.
   virtual void BeginTick(int /*tick*/) {}
   virtual void Predict(int cluster, mfem::real_t dt_step) = 0;
   virtual void Correct(int cluster, mfem::real_t dt_step) = 0;

   /// LTS Track-A A1 (per-tick exchange merge).  Called ONCE per tick, AFTER every
   /// Predict and BEFORE the first Correct — the only window in which all of this
   /// tick's pack inputs are final (I_[c] is produced by Predict; seam_coarse_dk_
   /// is written only by PrepareSeamCoarseForecast, also on the predict path) and
   /// no Correct has run yet.  A stepper that merges its seam exchanges packs all
   /// correcting clusters here and fires ONE collective per buffer.
   /// Default no-op => steppers that do not merge are unaffected.
   virtual void BeforeCorrects(const std::vector<int>& /*correct_clusters*/,
                               const std::vector<mfem::real_t>& /*dt_step*/) {}
   /// Counterpart to BeforeCorrects: drop any per-tick buffers so a later tick can
   /// never read stale data.  Default no-op.
   virtual void AfterCorrects() {}
};

/// Drive ONE sync interval from a precomputed tick table (Normative Scheduling):
/// per tick, predict every due cluster (any order), then correct every due
/// cluster FINE→COARSE (the table already sorts `correct_clusters` ascending).
void RunSyncInterval(const std::vector<LtsTick>& table, ILtsClusterStepper& stepper);

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LTS_STEPPER_HPP
