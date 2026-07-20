// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// lts_fault_stepper.hpp — the multi-rate FAULT-AWARE stepper (Phase 3 / P3-4).
//
// Extends the fault-free LtsBulkSyncStepper (lts_bulk_stepper.hpp) with the
// per-cluster FAULT HALF: when a due cluster is predicted, evaluate the bulk
// traction traces at THAT cluster's fault QPs, run the friction range Advance
// over the cluster's cluster-contiguous global-QP range at the cluster's step
// dt (writing the imposed states into a persistent global buffer), and — at the
// cluster's correct — apply the cluster's fault-face flux (via the corrector's
// fault_face_ids hook, which consumes the imposed states through the EXACT GTS
// fault flux, ProcessADERFaceToRHS_).
//
// Reduction to GTS: at Nc=1 the single cluster owns every element + every fault
// face, so Predict(0)/Correct(0) reproduce AdvanceADERWithSubStep_Spatial (the
// GTS ADER macro-step with fault friction) to machine-eps (the corrector's
// per-cluster face sweep reassociates the face sum vs the whole-mesh GTS sum;
// the nucleation is absolute vs GTS-incremental — same value, different rounding).
//
// np>=1 under D-2: fault faces are rank-INTERIOR (fault-locality partition,
// P-017), so the fault half is purely local (the interior-only per-cluster
// fault-QP eval + friction range Advance).  At np>1 the BULK half additionally
// has rank seams — handled by the SAME bulk-seam machinery as LtsBulkSyncStepper
// (Correct forwards cluster_id_for_seam + schedule to AdvanceADERClusterBulk;
// Predict runs PrepareSeamCoarseForecast for the diff-1 coarse forecast when
// SetExchangeProviderDk is on).  There are NO cross-rank FAULT faces (D-2,
// asserted in the ctor).

#ifndef MFEM_SEAS_LTS_FAULT_STEPPER_HPP
#define MFEM_SEAS_LTS_FAULT_STEPPER_HPP

#include "wave_operator.hpp"
#include "lts_layout.hpp"
#include "lts_stepper.hpp"
#include "friction_iterator.hpp"     // IFrictionIterator
#include "nucleation_method.hpp"     // INucleationMethod
#include "fault_face_flux.hpp"       // DOFData

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace mfem
{
namespace seas
{

template <typename MeshType>
class LtsFaultSyncStepper : public ILtsClusterStepper
{
public:
   /// `iterator` drives the per-cluster friction range Advance; `dof_data` /
   /// `fault_coords` are the (reordered, cluster-contiguous) per-QP fault state;
   /// `nuc` supplies the ABSOLUTE range nucleation (D-3).  All references must
   /// outlive the stepper.
   LtsFaultSyncStepper(WaveOperator<MeshType> &wave, const LtsLayout &layout,
                       const std::vector<int> &cluster_id, mfem::Vector &Q,
                       int order, mfem::real_t dt_base,
                       IFrictionIterator &iterator,
                       std::vector<DOFData> &dof_data,
                       const std::vector<mfem::Vector> &fault_coords,
                       INucleationMethod *nuc)
      : wave_(wave), layout_(layout), cluster_id_(cluster_id), Q_(Q),
        order_(order), dt_base_(dt_base), iterator_(iterator),
        dof_data_(dof_data), fault_coords_(fault_coords), nuc_(nuc)
   {
      const int nc = layout_.num_clusters;
      block_ = NUM_STATE * wave_.GetNDof();
      I_.resize(nc);
      dk_.Resize(static_cast<int>(layout_.provider_elems.size()), order_, block_);
      buf_.Resize(static_cast<int>(layout_.consumer_owner_elems.size()), block_);

      fids_.resize(nc); froles_.resize(nc); fnbr_.resize(nc);
      for (int c = 0; c < nc; ++c)
      {
         const LtsCluster &cl = layout_.clusters[c];
         fids_[c] = cl.face_ids;
         fnbr_[c] = cl.face_nbr;
         froles_[c].resize(cl.faces.size());
         for (std::size_t i = 0; i < cl.faces.size(); ++i)
         { froles_[c][i] = static_cast<int>(cl.faces[i]); }
      }

      // --- per-cluster fault-QP range (cluster-contiguous by the P-006 reorder) ---
      // The wave operator's interior fault-QP layout is: face at POSITION i in the
      // (reordered) fault_interior_faces_ list owns global QP range [i*nbf,(i+1)*nbf).
      // Under the P-006 reorder that list is sorted by (cluster, mesh-face-id), so a
      // cluster's faces occupy a contiguous block of positions ⇒ a contiguous QP range.
      n_total_fault_qps_ = wave_.GetNumTotalFaultQPs();
      const int nbf = wave_.GetNbfPerFace();
      nbf_ = nbf;
      // REVIEW P3-4 A2 / p4-fault: D-2 invariant — the fault interleave has NO
      // shared (cross-rank) fault faces AT ANY np (fault-locality partitioning,
      // P-017); the per-cluster QP ranges are built only from the interior
      // fault-face list, so any shared QP slot would stay zero in the imposed-state
      // buffer and the corrector would silently consume zeroed states.  Fail loud
      // (cross-rank shared-fault stepping is not supported).
      MFEM_VERIFY(wave_.GetNumSharedFaultQPs() == 0,
                  "LtsFaultSyncStepper: expected zero shared fault QPs (D-2, "
                  "fault-locality), got " << wave_.GetNumSharedFaultQPs()
                  << " — cross-rank shared-fault stepping is not supported.");
      // REVIEW P3-4 A4: on a FAULT-BEARING rank the fault-QP layout must exist
      // (SetFaultDOFData ran); otherwise nbf=0 collapses every range to [0,0), the
      // imposed-state buffer stays unset, and Correct() silently falls to the inline
      // friction re-solve.  LTS Phase 4a: at np>1 a rank may own NO fault faces
      // (fault-free subdomain) — that is legitimate (the fault half is then a no-op,
      // the rank does bulk + bulk seams only).  So require nbf>0 only when this rank
      // actually has fault QPs.
      MFEM_VERIFY(n_total_fault_qps_ == 0 || nbf > 0,
                  "LtsFaultSyncStepper: n_total_fault_qps=" << n_total_fault_qps_
                  << " > 0 but nbf_per_face=" << nbf
                  << " — the fault DOF layout must be built (SetFaultDOFData) before "
                  "constructing the fault stepper.");
      const mfem::Array<int> &fif = wave_.GetFaultInteriorFaces();
      std::map<int, int> face_to_pos;   // mesh face id -> position in fif
      for (int i = 0; i < fif.Size(); ++i) { face_to_pos[fif[i]] = i; }
      qp_begin_.assign(nc, 0);
      qp_end_.assign(nc, 0);
      cluster_fault_faces_.resize(nc);
      for (int c = 0; c < nc; ++c)
      {
         const std::vector<int> &ff = layout_.clusters[c].fault_faces;
         cluster_fault_faces_[c].assign(ff.begin(), ff.end());
         if (ff.empty()) { continue; }
         // REVIEW P3-4 R3 (corrected): BuildLtsLayout tracks fault faces ONLY in
         // cl.fault_faces and never adds them to cl.face_ids (the bulk sweep), so the
         // corrector's role-driven sweep can never see them and the fault loop is
         // their SOLE processor.  Assert that invariant defensively: no cluster fault
         // face may appear in the owned-face sweep (any role) — a future layout that
         // put it there would double- or mis-process it.
         {
            const LtsCluster &clc = layout_.clusters[c];
            std::set<int> sweep_faces(clc.face_ids.begin(), clc.face_ids.end());
            for (int f : ff)
            {
               MFEM_VERIFY(sweep_faces.count(f) == 0,
                           "LtsFaultSyncStepper: cluster " << c << " fault face "
                           << f << " also appears in the bulk face sweep (face_ids) "
                           "— it would be double-processed (layout inconsistency; "
                           "fault faces must live only in cl.fault_faces).");
            }
         }
         int pos_lo = std::numeric_limits<int>::max(), pos_hi = -1;
         for (int f : ff)
         {
            auto it = face_to_pos.find(f);
            MFEM_VERIFY(it != face_to_pos.end(),
                        "LtsFaultSyncStepper: cluster " << c << " fault face "
                        << f << " is not in the wave operator's interior fault-face "
                        "list (layout / operator fault-face sets disagree).");
            if (it->second < pos_lo) { pos_lo = it->second; }
            if (it->second > pos_hi) { pos_hi = it->second; }
         }
         qp_begin_[c] = static_cast<std::size_t>(pos_lo) * nbf;
         qp_end_[c]   = static_cast<std::size_t>(pos_hi + 1) * nbf;
         // The reorder makes each cluster's fault faces contiguous — assert it,
         // else the range sweep would silently span another cluster's QPs.
         MFEM_VERIFY(qp_end_[c] - qp_begin_[c]
                     == ff.size() * static_cast<std::size_t>(nbf),
                     "LtsFaultSyncStepper: cluster " << c << " fault QPs are not "
                     "contiguous (positions [" << pos_lo << "," << pos_hi
                     << "] span " << (pos_hi - pos_lo + 1) << " faces but the cluster "
                     "owns " << ff.size() << ").  The P-006 reorder / cluster ids "
                     "disagree.");
      }

      // --- persistent global imposed-state buffer (disjoint per-cluster ranges) ---
      const std::size_t n_words =
         static_cast<std::size_t>(NUM_STATE) *
         static_cast<std::size_t>(n_total_fault_qps_);
      I_imp_plus_.assign(n_words, 0.0);
      I_imp_minus_.assign(n_words, 0.0);
      // Point the wave operator at the persistent buffer ONCE; Predict(c) writes
      // cluster c's range, Correct(c) reads it through the fault-face flux.
      if (n_total_fault_qps_ > 0)
      {
         wave_.SetSubStepFaultImposedStates(I_imp_plus_.data(),
                                            I_imp_minus_.data(),
                                            n_total_fault_qps_);
      }

      // Capture the ORIGINAL configured friction sub-steps ONCE.  Predict()
      // rescales the iterator to each cluster's dt_step and thereby MUTATES the
      // shared iterator; rescaling always from these originals (not from the
      // mutated GetDeltaT()) keeps Predict() order-independent and obviously
      // correct rather than relying on proportional-rescale invariance.
      cfg_dT0_ = iterator_.GetDeltaT();
      cfg_w0_  = iterator_.GetTimeWeights();

      // Q_pointwise scratch: sized to EXACTLY O = the number of friction sub-steps
      // (REVIEW P3-4 A5) — the range Advance's RunSubSteps_ requires
      // Q_pointwise_*.size() == O exactly.  O is fixed once here (cfg_dT0_.size()).
      Qpw_plus_.resize(cfg_dT0_.size());
      Qpw_minus_.resize(cfg_dT0_.size());
   }

   ~LtsFaultSyncStepper() override
   {
      // Detach the imposed-state buffer so a later GTS path cannot read freed data.
      if (n_total_fault_qps_ > 0) { wave_.ResetSubStepFaultImposedStates(); }
   }

   void SetSyncInterval(mfem::real_t t_s, mfem::real_t T_actual)
   { t_s_ = t_s; T_actual_ = T_actual; }

   /// LTS Phase 4a Stage 2: enable the coarse provider-D(k) ghost exchange for
   /// the BULK diff-1 rank seams (np>1).  Must match the tick table's
   /// LtsGlobalMeta.exchange_bulk_provider_dk.  Default off (np=1 / single cluster).
   void SetExchangeProviderDk(bool v) { exchange_dk_ = v; }

   void BeginTick(int tick) override { tick_ = tick; }

   void Predict(int c, mfem::real_t dt_step) override
   {
      // Rescale the friction iterator's sub-steps to THIS cluster's step dt_step
      // (mirrors AdvanceADERWithSubStep_Spatial): the iterator holds Σ deltaT for
      // some reference dt, and its range Advance asserts Σ deltaT == dt_step.  The
      // predictor eval nodes MUST be the midpoints of those SAME scaled sub-steps
      // (the fault eval reads Q_per_node at them; the corrector's I is node-
      // independent, so the bulk half is unaffected).
      const int O = static_cast<int>(cfg_dT0_.size());
      MFEM_VERIFY(O >= 1, "LtsFaultSyncStepper: iterator sub-steps not set.");
      mfem::real_t cfg_sum = 0.0;
      for (int o = 0; o < O; ++o) { cfg_sum += cfg_dT0_[o]; }
      MFEM_VERIFY(cfg_sum > 0.0, "LtsFaultSyncStepper: Σ deltaT ≤ 0.");
      const mfem::real_t scale = dt_step / cfg_sum;
      std::vector<mfem::real_t> dT(O);
      for (int o = 0; o < O; ++o) { dT[o] = cfg_dT0_[o] * scale; }
      iterator_.SetSubSteps(dT, cfg_w0_);
      const std::vector<mfem::real_t> &deltaT = iterator_.GetDeltaT();
      std::vector<mfem::real_t> tau(O);
      mfem::real_t acc = 0.0;
      for (int o = 0; o < O; ++o) { tau[o] = acc + 0.5 * deltaT[o]; acc += deltaT[o]; }

      // 1. Bulk predictor (retains the cluster's D(k) for its providers).
      wave_.ComputeADERSubStepStatesAndIntegralCluster(
         layout_.clusters[c].elems.data(),
         static_cast<int>(layout_.clusters[c].elems.size()),
         Q_, dt_step, order_, tau, Qn_scratch_, I_[c],
         dk_.data.data(), layout_.provider_slot_of_elem.data());

      // 2. Fault half: evaluate the bulk traction traces at THIS cluster's fault
      //    QPs and advance its friction over [qp_begin,qp_end) at dt_step.  The
      //    imposed states land in the persistent global buffer's cluster range;
      //    Correct(c) consumes them.  Nucleation is the ABSOLUTE range form (D-3),
      //    idempotent so each cluster forcing its own range at its own stage times
      //    does not double-count.
      if (!cluster_fault_faces_[c].empty())
      {
         MFEM_VERIFY(static_cast<int>(Qn_scratch_.size()) == O,
                     "LtsFaultSyncStepper: predictor returned "
                     << Qn_scratch_.size() << " nodes, expected " << O);
         if (static_cast<int>(Qpw_plus_.size()) < O)
         { Qpw_plus_.resize(O); Qpw_minus_.resize(O); }
         // A1: evaluate ONLY this cluster's contiguous interior-fault-face block
         // [fi_begin, fi_end) — writes exactly the [qp_begin,qp_end) slots the
         // range friction Advance reads.  No whole-mesh cost, no uninitialized
         // read, no per-cluster collective.  fi = qp / nbf (ranges are pos*nbf).
         const int fi_begin = static_cast<int>(qp_begin_[c] / nbf_);
         const int fi_end   = static_cast<int>(qp_end_[c]   / nbf_);
         for (int o = 0; o < O; ++o)
         {
            wave_.EvaluateBulkAtFaultQPsCanonicalRange(
               Qn_scratch_[o], fi_begin, fi_end, Qpw_plus_[o], Qpw_minus_[o]);
         }
         const mfem::real_t t_step_start =
            t_s_ + dt_base_ * static_cast<mfem::real_t>(tick_);
         INucleationMethod *nuc = nuc_;
         std::vector<DOFData> &dd = dof_data_;
         const std::size_t qb = qp_begin_[c], qe = qp_end_[c];
         auto nuc_cb = [nuc, &dd, qb, qe](mfem::real_t t, mfem::real_t /*dt*/)
         {
            if (nuc) { nuc->ApplyAbsolute(dd, t, qb, qe); }
         };
         iterator_.Advance(qb, qe, dof_data_, fault_coords_,
                           Qpw_plus_, Qpw_minus_, dt_step, t_step_start,
                           I_imp_plus_.data(), I_imp_minus_.data(), nuc_cb);
      }

      // 3. BULK diff-1 rank seams (np>1): retain + exchange this (coarse) cluster's
      //    seam-provider D(k) so a finer cross-rank neighbour can integrate its
      //    forecast.  Fault faces are rank-interior (D-2), so the fault half above
      //    needs no exchange; only the bulk half has rank seams.  c==0 is finest
      //    (never a provider); the c>=1 gate is rank-uniform (matched collectives).
      if (exchange_dk_ && c >= 1)
      {
         wave_.PrepareSeamCoarseForecast(c, Q_, dt_step, order_, tau, tick_);
      }
   }

   void Correct(int c, mfem::real_t dt_step) override
   {
      const int nf = static_cast<int>(fids_[c].size());
      std::vector<mfem::real_t> sa(nf, 0.0), sb(nf, 0.0);
      for (int i = 0; i < nf; ++i)
      {
         if (froles_[c][i] != static_cast<int>(FaceRole::ConsumerFine)) { continue; }
         const int coarse = fnbr_[c][i];
         const int c_coarse = cluster_id_[coarse];
         const long long period = 1LL << c_coarse;
         const mfem::real_t t_origin =
            t_s_ + dt_base_ * static_cast<mfem::real_t>(period)
            * std::floor(static_cast<double>(tick_) / period);
         const mfem::real_t t_tick = t_s_ + dt_base_ * tick_;
         sa[i] = t_tick - t_origin;
         sb[i] = std::min(t_tick + dt_step, t_s_ + T_actual_) - t_origin;
      }

      // Bulk + fault corrector.  The imposed-state buffer (set once in the ctor)
      // still points here; ProcessADERFaceToRHS_ reads cluster c's range for the
      // cluster's fault faces.
      const int *ff = cluster_fault_faces_[c].empty()
                      ? nullptr : cluster_fault_faces_[c].data();
      const int nff = static_cast<int>(cluster_fault_faces_[c].size());
      wave_.AdvanceADERClusterBulk(
         layout_.clusters[c].elems.data(),
         static_cast<int>(layout_.clusters[c].elems.size()),
         fids_[c].data(), froles_[c].data(), fnbr_[c].data(), nf,
         dt_step, order_, I_[c], Q_,
         dk_.data.data(), order_, layout_.provider_slot_of_elem.data(),
         sa.data(), sb.data(), &buf_, layout_.buffer_slot_of_elem.data(),
         ff, nff,
         // LTS Phase 4a: process this cluster's BULK rank-seam faces (diff-0 +
         // diff-1) via the seam corrector; the schedule gives the diff-1 [a,b].
         // np=1 / no shared faces => no-op.
         /*cluster_id_for_seam=*/c,
         /*t_s=*/t_s_, /*dt_base=*/dt_base_, /*tick=*/tick_, /*T_actual=*/T_actual_,
         // LTS Phase 4b: premultiplied coarse-forecast exchange iff this correcting
         // cluster has a coarser neighbour across a rank seam (c < num_clusters-1).
         /*exchange_forecast=*/exchange_dk_ && (c < layout_.num_clusters - 1));
   }

   /// Invariant (ii): every accumulate buffer is zero at each sync point.
   bool BuffersZero() const { return buf_.AllZero(); }

private:
   WaveOperator<MeshType> &wave_;
   const LtsLayout &layout_;
   const std::vector<int> &cluster_id_;
   mfem::Vector &Q_;
   int order_, block_ = 0;
   mfem::real_t dt_base_ = 0.0, t_s_ = 0.0, T_actual_ = 0.0;
   int tick_ = 0;
   bool exchange_dk_ = false;   // Stage 2: bulk coarse provider-D(k) ghost exchange

   IFrictionIterator &iterator_;
   std::vector<DOFData> &dof_data_;
   const std::vector<mfem::Vector> &fault_coords_;
   INucleationMethod *nuc_ = nullptr;
   std::vector<mfem::real_t> cfg_dT0_, cfg_w0_;  // ORIGINAL configured sub-steps

   std::vector<mfem::Vector> I_;           // per-cluster time integral
   std::vector<mfem::Vector> Qn_scratch_;  // predictor sub-step scratch (reused)
   LtsDkStore dk_;
   LtsAccumulateBuffers buf_;
   std::vector<std::vector<int>> fids_, froles_, fnbr_;

   // Fault half state.
   int n_total_fault_qps_ = 0;
   int nbf_ = 0;                                       // QPs per fault face
   std::vector<std::size_t> qp_begin_, qp_end_;        // per cluster
   std::vector<std::vector<int>> cluster_fault_faces_; // per cluster (mesh face ids)
   std::vector<mfem::real_t> I_imp_plus_, I_imp_minus_;
   std::vector<std::vector<mfem::real_t>> Qpw_plus_, Qpw_minus_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_LTS_FAULT_STEPPER_HPP
