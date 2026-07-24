// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// lts_bulk_stepper.hpp — the multi-rate FAULT-FREE bulk stepper (Phase 2).
//
// Wires the wave operator's per-cluster predictor + corrector behind the pure
// tick loop (RunSyncInterval), managing the per-cluster time-integral storage,
// the provider D(k) store, the accumulate buffers, and — the piece the
// single-cluster path never exercises — the CLOSED-FORM consumer-face
// sub-interval [a,b] derived from the tick schedule (Normative Scheduling):
//
//   t(tick)          = t_s + tick*dt_base
//   t_origin(coarse) = t_s + dt_base*2^c_coarse * floor(tick / 2^c_coarse)
//   a = t(tick) - t_origin(coarse)
//   b = min(t(tick) + dt_fine, t_s + T_actual) - t_origin(coarse)
//
// This class IS the core of the driver's Phase-2 sync-interval loop; the driver
// constructs one, then per sync interval calls SetSyncInterval + RunSyncInterval
// with the tick table.  It is unit-validated by test_lts_multirate (conservation)
// and the degenerate 2-cluster == GTS check.

#ifndef MFEM_SEAS_LTS_BULK_STEPPER_HPP
#define MFEM_SEAS_LTS_BULK_STEPPER_HPP

#include "wave_operator.hpp"
#include "lts_layout.hpp"
#include "lts_stepper.hpp"

#include <cmath>
#include <array>
#include <vector>

namespace mfem
{
namespace seas
{

template <typename MeshType>
class LtsBulkSyncStepper : public ILtsClusterStepper
{
public:
   /// `cluster_id` is the per-element cluster id (used to find a consumer face's
   /// coarse neighbour's level for t_origin).  `Q` is advanced in place.
   LtsBulkSyncStepper(WaveOperator<MeshType>& wave, const LtsLayout& layout,
                      const std::vector<int>& cluster_id, mfem::Vector& Q,
                      int order, mfem::real_t dt_base)
      : wave_(wave), layout_(layout), cluster_id_(cluster_id), Q_(Q),
        order_(order), dt_base_(dt_base)
   {
      const int nc = layout_.num_clusters;
      block_ = NUM_STATE * wave_.GetNDof();
      I_.resize(nc);
      dk_.Resize(static_cast<int>(layout_.provider_elems.size()), order_, block_);
      buf_.Resize(static_cast<int>(layout_.consumer_owner_elems.size()), block_);

      // Precompute per-cluster face arrays (ids / int roles / neighbours).
      fids_.resize(nc); froles_.resize(nc); fnbr_.resize(nc);
      for (int c = 0; c < nc; ++c)
      {
         const LtsCluster& cl = layout_.clusters[c];
         fids_[c] = cl.face_ids;
         fnbr_[c] = cl.face_nbr;
         froles_[c].resize(cl.faces.size());
         for (std::size_t i = 0; i < cl.faces.size(); ++i)
         { froles_[c][i] = static_cast<int>(cl.faces[i]); }
      }
   }

   /// Set the current sync interval [t_s, t_s+T_actual] before RunSyncInterval.
   void SetSyncInterval(mfem::real_t t_s, mfem::real_t T_actual)
   { t_s_ = t_s; T_actual_ = T_actual; }

   void BeginTick(int tick) override { tick_ = tick; }

   /// LTS Phase 4a Stage 2: enable the coarse provider-D(k) ghost exchange (for
   /// diff-1 rank seams).  Must be set consistently with the tick table's
   /// `LtsGlobalMeta.exchange_bulk_provider_dk` (else the matched-collective count
   /// diverges).  Default off (Stage-1 diff-0 / np=1 behaviour).
   void SetExchangeProviderDk(bool v) { exchange_dk_ = v; }

   /// LTS Track-A A1: merge this tick's per-correct seam exchanges into ONE
   /// collective per buffer (125 -> 64 rounds/sync at Nc=6).  Opt-in; default OFF
   /// keeps the per-correct path byte-for-byte.  MUST be set identically on every
   /// rank (it changes the collective count, which the tick table predicts).
   void SetMergeTickExchanges(bool v) { merge_tick_exchanges_ = v; }
   bool MergeTickExchanges() const { return merge_tick_exchanges_; }

   /// A1 merge window (see ILtsClusterStepper::BeforeCorrects).  Packs every
   /// correcting cluster's seam forecast + time-integral and fires exactly two
   /// collectives; each Correct below then reads the merged buffers instead of
   /// exchanging for itself.
   void BeforeCorrects(const std::vector<int>& correct_clusters,
                       const std::vector<mfem::real_t>& dt_step) override
   {
      if (!merge_tick_exchanges_ || correct_clusters.empty()) { return; }
      const int n = static_cast<int>(correct_clusters.size());
      const int nc = static_cast<int>(layout_.clusters.size());
      MFEM_VERIFY(n <= static_cast<int>(fc_pc_.size()),
                  "BeforeCorrects: more correcting clusters than the A1 scratch holds.");
      I_ptrs_.resize(n); dt_pc_.resize(n);
      for (int k = 0; k < n; ++k)
      {
         const int c = correct_clusters[k];
         MFEM_VERIFY(c >= 0 && c < nc && c < static_cast<int>(I_.size()),
                     "BeforeCorrects: cluster id out of range.");
         I_ptrs_[k] = &I_[c];
         dt_pc_[k]  = dt_step[c];
         // SAME predicate the per-correct path passes as `exchange_forecast`
         // (see Correct below): this cluster has a coarser neighbour.
         fc_pc_[k]  = (exchange_dk_ && c < nc - 1);
      }
      wave_.PrepareClusterSeamExchangeTick(
         correct_clusters.data(), n, I_ptrs_.data(), dt_pc_.data(), fc_pc_.data(),
         t_s_, dt_base_, tick_, T_actual_, order_);
   }

   void AfterCorrects() override
   {
      if (merge_tick_exchanges_) { wave_.ClearClusterSeamExchangeTick(); }
   }

   void Predict(int c, mfem::real_t dt_step) override
   {
      std::vector<mfem::real_t> tau(order_);
      for (int o = 0; o < order_; ++o) { tau[o] = dt_step * (o + 0.5) / order_; }
      wave_.ComputeADERSubStepStatesAndIntegralCluster(
         layout_.clusters[c].elems.data(),
         static_cast<int>(layout_.clusters[c].elems.size()),
         Q_, dt_step, order_, tau, Qn_scratch_, I_[c],
         dk_.data.data(), layout_.provider_slot_of_elem.data());
      // Stage 2: retain + exchange this (coarse) cluster's seam-provider D(k) so a
      // finer cross-rank neighbour can integrate its forecast.  c==0 is finest
      // (never a provider); the c>=1 gate is rank-uniform (matched).
      if (exchange_dk_ && c >= 1)
      {
         wave_.PrepareSeamCoarseForecast(c, Q_, dt_step, order_, tau, tick_);
      }
   }

   void Correct(int c, mfem::real_t dt_step) override
   {
      const int nf = static_cast<int>(fids_[c].size());
      std::vector<mfem::real_t> sa(nf, 0.0), sb(nf, 0.0);
      // Closed-form consumer-face sub-interval [a,b] (relative to the coarse
      // neighbour's expansion point).  dt_step here is the FINE cluster's step.
      for (int i = 0; i < nf; ++i)
      {
         if (froles_[c][i] != static_cast<int>(FaceRole::ConsumerFine)) { continue; }
         const int coarse = fnbr_[c][i];
         const int c_coarse = cluster_id_[coarse];
         const long long period = 1LL << c_coarse;
         const mfem::real_t t_origin = t_s_ + dt_base_ * static_cast<mfem::real_t>(period)
                                       * std::floor(static_cast<double>(tick_) / period);
         const mfem::real_t t_tick = t_s_ + dt_base_ * tick_;
         sa[i] = t_tick - t_origin;
         sb[i] = std::min(t_tick + dt_step, t_s_ + T_actual_) - t_origin;
      }
      wave_.AdvanceADERClusterBulk(
         layout_.clusters[c].elems.data(),
         static_cast<int>(layout_.clusters[c].elems.size()),
         fids_[c].data(), froles_[c].data(), fnbr_[c].data(), nf,
         dt_step, order_, I_[c], Q_,
         dk_.data.data(), order_, layout_.provider_slot_of_elem.data(),
         sa.data(), sb.data(), &buf_, layout_.buffer_slot_of_elem.data(),
         // LTS Phase 4a: no fault faces on the bulk path (fault-free); pass the
         // cluster id so the corrector applies this cluster's rank-seam flux
         // (no-op at np=1 / no shared faces).  Stage 2: the schedule gives the
         // diff-1 coarse forecast's [a,b].
         /*fault_face_ids=*/nullptr, /*n_fault_faces=*/0,
         /*cluster_id_for_seam=*/c,
         /*t_s=*/t_s_, /*dt_base=*/dt_base_, /*tick=*/tick_, /*T_actual=*/T_actual_,
         // LTS Phase 4b: run the premultiplied coarse-forecast exchange iff this
         // correcting cluster has a coarser neighbour across a rank seam.  Matches
         // the correct-side +1 collective per c<num_clusters-1 in BuildTickTable.
         /*exchange_forecast=*/exchange_dk_ && (c < layout_.num_clusters - 1));
   }

   /// Invariant (ii): every accumulate buffer is zero at each sync point.
   bool BuffersZero() const { return buf_.AllZero(); }

private:
   WaveOperator<MeshType>& wave_;
   const LtsLayout& layout_;
   const std::vector<int>& cluster_id_;
   mfem::Vector& Q_;
   int order_, block_ = 0;
   mfem::real_t dt_base_ = 0.0, t_s_ = 0.0, T_actual_ = 0.0;
   int tick_ = 0;
   bool exchange_dk_ = false;   // Stage 2: coarse provider-D(k) ghost exchange

   bool merge_tick_exchanges_ = false;    // A1: per-tick merged seam exchanges
   std::vector<const mfem::Vector *> I_ptrs_;   // A1 scratch (BeforeCorrects)
   std::vector<mfem::real_t>         dt_pc_;    // A1 scratch
   // std::vector<bool> is a bitfield with no bool* data(); clusters are capped at
   // 62 by BuildTickTable, so a fixed array needs no allocation and yields bool*.
   std::array<bool, 64>              fc_pc_{};   // A1 scratch
   std::vector<mfem::Vector> I_;          // per-cluster time integral
   std::vector<mfem::Vector> Qn_scratch_; // predictor sub-step scratch (unused output)
   LtsDkStore dk_;
   LtsAccumulateBuffers buf_;
   std::vector<std::vector<int>> fids_, froles_, fnbr_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LTS_BULK_STEPPER_HPP
