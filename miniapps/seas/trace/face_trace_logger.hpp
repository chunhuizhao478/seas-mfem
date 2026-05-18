// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Face-level runtime tracer for SEAS fault diagnostics.
//
// Rank-local, append-only, barrier-free. No MPI calls inside any method.
// Each rank writes only its own files; non-traced ranks produce nothing.

#ifndef MFEM_SEAS_FACE_TRACE_LOGGER_HPP
#define MFEM_SEAS_FACE_TRACE_LOGGER_HPP

#include "mfem.hpp"
#include "../domain/elasticity_operator.hpp"
#include "../fault/fault_geometry.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <set>

namespace mfem
{
namespace seas
{

/// Per-DOF metadata for a traced fault face DOF.
struct TraceFaceSpec
{
   int fi;              // owned fault face index
   int dof;             // DOF index within face (0..nbf-1)
   int owned_dof_idx;   // global owned DOF index (fi * nbf + dof)
   std::string role;    // "target" or "control"
   real_t cx, cz;       // along-strike (x2), depth (x3) of this DOF
   std::string face_type;   // "interior" or "shared"

   // BP5 metadata (written once at startup)
   real_t depth;
   real_t a, Dc;
   real_t tau_pre_dip, tau_pre_strike;
   real_t V_init_dip, V_init_strike;
};

struct TraceConfig
{
   // Explicit face list (rank -> {fi...})
   int explicit_rank = -1;
   std::vector<int> explicit_fi;

   // Coordinate window (auto-discovery, any rank)
   bool use_coord_window = false;
   real_t x2_min = 0, x2_max = 0, x3_min = 0, x3_max = 0;

   // Control faces
   int num_control_faces = 2;

   // Max traced faces per rank (guard against huge CSV files)
   int max_traced_faces = 50;

   // Sampling
   int warmup_steps = 32;
   int periodic_interval = 500;
   int flush_interval = 64;

   // Event thresholds
   real_t sigma_n_eff_warn = 1e6;      // 1 MPa
   real_t V_mag_warn = 1e-2;           // m/s
   real_t corr_ratio_warn = 0.5;
   real_t jump_factor_warn = 10.0;     // sudden-jump factor

   // Output directory
   std::string output_dir = ".";
};

template <typename MeshType>
class FaceTraceLogger
{
public:
   FaceTraceLogger(const TraceConfig &cfg, int rank)
      : cfg_(cfg), rank_(rank) {}

   ~FaceTraceLogger()
   {
      if (active_ && !finalized_) { Finalize(); }
   }

   void SelectFaces(
      const ElasticityDomainOperator<MeshType> &domain_op,
      const FaultGeometry<MeshType> &geom)
   {
      if (disabled_) { return; }

      nbf_per_face_ = domain_op.GetNbfPerFace();
      const Vector &x2 = geom.GetCoordsX2();
      const Vector &x3 = geom.GetCoordsX3();
      int num_owned_dofs = x2.Size();
      int num_owned_faces = num_owned_dofs / nbf_per_face_;

      // Compute per-face centroids for selection
      std::vector<real_t> face_cx(num_owned_faces, 0.0);
      std::vector<real_t> face_cz(num_owned_faces, 0.0);
      for (int f = 0; f < num_owned_faces; f++)
      {
         for (int k = 0; k < nbf_per_face_; k++)
         {
            int d = f * nbf_per_face_ + k;
            face_cx[f] += x2(d);
            face_cz[f] += x3(d);
         }
         face_cx[f] /= nbf_per_face_;
         face_cz[f] /= nbf_per_face_;
      }

      std::set<int> selected;

      // Mode 1: explicit face list
      if (rank_ == cfg_.explicit_rank)
      {
         for (int fi : cfg_.explicit_fi)
         {
            if (fi >= 0 && fi < num_owned_faces)
            {
               selected.insert(fi);
            }
         }
      }

      // Mode 2: coordinate window
      if (cfg_.use_coord_window)
      {
         for (int f = 0; f < num_owned_faces; f++)
         {
            if (face_cx[f] >= cfg_.x2_min && face_cx[f] <= cfg_.x2_max &&
                face_cz[f] >= cfg_.x3_min && face_cz[f] <= cfg_.x3_max)
            {
               selected.insert(f);
            }
         }
      }

      // Enforce max_traced_faces limit on targets
      if (static_cast<int>(selected.size()) > cfg_.max_traced_faces)
      {
         std::set<int> trimmed;
         int count = 0;
         for (int fi : selected)
         {
            if (count >= cfg_.max_traced_faces) { break; }
            trimmed.insert(fi);
            count++;
         }
         selected = std::move(trimmed);
      }

      // Gather face metadata from domain operator
      const int num_interior = domain_op.GetFaultInteriorFaces().Size();
      const auto &owned_map = domain_op.GetOwnedFaultFaceMap();

      auto get_face_type = [&](int fi) -> std::string
      {
         if (fi < static_cast<int>(owned_map.Size()))
         {
            int local_fi = owned_map[fi];
            return (local_fi < num_interior) ? "interior" : "shared";
         }
         return "unknown";
      };

      // Get BP5 per-DOF parameters from geometry
      const Vector &geom_a = geom.GetAValues();
      const Vector &geom_dc = geom.GetDcValues();
      const Vector &geom_depths = geom.GetDepths();
      const Vector &geom_tau_pre = geom.GetTauPre();
      const Vector &geom_v_init = geom.GetVInit();

      // Build traced DOF specs from selected faces
      for (int fi : selected)
      {
         std::string ft = get_face_type(fi);
         for (int kk = 0; kk < nbf_per_face_; kk++)
         {
            int owned_dof = fi * nbf_per_face_ + kk;
            TraceFaceSpec spec;
            spec.fi = fi;
            spec.dof = kk;
            spec.owned_dof_idx = owned_dof;
            spec.role = "target";
            spec.cx = x2(owned_dof);
            spec.cz = x3(owned_dof);
            spec.face_type = ft;

            // BP5 metadata
            spec.depth = (owned_dof < geom_depths.Size())
                         ? geom_depths(owned_dof) : 0.0;
            spec.a = (owned_dof < geom_a.Size())
                      ? geom_a(owned_dof) : 0.0;
            spec.Dc = (owned_dof < geom_dc.Size())
                       ? geom_dc(owned_dof) : 0.0;
            spec.tau_pre_dip = (2 * owned_dof < geom_tau_pre.Size())
                               ? geom_tau_pre(2 * owned_dof) : 0.0;
            spec.tau_pre_strike = (2 * owned_dof + 1 < geom_tau_pre.Size())
                                  ? geom_tau_pre(2 * owned_dof + 1) : 0.0;
            spec.V_init_dip = (2 * owned_dof < geom_v_init.Size())
                              ? geom_v_init(2 * owned_dof) : 0.0;
            spec.V_init_strike = (2 * owned_dof + 1 < geom_v_init.Size())
                                 ? geom_v_init(2 * owned_dof + 1) : 0.0;

            traced_dofs_.push_back(spec);
         }
      }

      // Mode 3: control faces (nearest to cluster centroid, outside selection)
      if (cfg_.num_control_faces > 0 && !traced_dofs_.empty())
      {
         real_t cluster_cx = 0, cluster_cz = 0;
         for (auto &s : traced_dofs_)
         {
            cluster_cx += s.cx;
            cluster_cz += s.cz;
         }
         cluster_cx /= traced_dofs_.size();
         cluster_cz /= traced_dofs_.size();

         std::vector<std::pair<real_t, int>> dist_fi;
         for (int f = 0; f < num_owned_faces; f++)
         {
            if (selected.count(f)) { continue; }
            real_t dx = face_cx[f] - cluster_cx;
            real_t dz = face_cz[f] - cluster_cz;
            dist_fi.push_back({dx * dx + dz * dz, f});
         }
         std::sort(dist_fi.begin(), dist_fi.end());

         int n_ctrl = std::min(cfg_.num_control_faces,
                               static_cast<int>(dist_fi.size()));
         for (int i = 0; i < n_ctrl; i++)
         {
            int fi = dist_fi[i].second;
            std::string ft = get_face_type(fi);
            for (int kk = 0; kk < nbf_per_face_; kk++)
            {
               int owned_dof = fi * nbf_per_face_ + kk;
               TraceFaceSpec spec;
               spec.fi = fi;
               spec.dof = kk;
               spec.owned_dof_idx = owned_dof;
               spec.role = "control";
               spec.cx = x2(owned_dof);
               spec.cz = x3(owned_dof);
               spec.face_type = ft;
               spec.depth = (owned_dof < geom_depths.Size())
                            ? geom_depths(owned_dof) : 0.0;
               spec.a = (owned_dof < geom_a.Size())
                         ? geom_a(owned_dof) : 0.0;
               spec.Dc = (owned_dof < geom_dc.Size())
                          ? geom_dc(owned_dof) : 0.0;
               spec.tau_pre_dip = (2 * owned_dof < geom_tau_pre.Size())
                                  ? geom_tau_pre(2 * owned_dof) : 0.0;
               spec.tau_pre_strike = (2 * owned_dof + 1 < geom_tau_pre.Size())
                                     ? geom_tau_pre(2 * owned_dof + 1) : 0.0;
               spec.V_init_dip = (2 * owned_dof < geom_v_init.Size())
                                 ? geom_v_init(2 * owned_dof) : 0.0;
               spec.V_init_strike = (2 * owned_dof + 1 < geom_v_init.Size())
                                    ? geom_v_init(2 * owned_dof + 1) : 0.0;
               traced_dofs_.push_back(spec);
            }
         }
      }

      active_ = !traced_dofs_.empty();
      if (!active_) { return; }

      // Build face-level index for summaries (unique fi values)
      std::set<int> unique_faces;
      for (auto &s : traced_dofs_) { unique_faces.insert(s.fi); }
      num_traced_faces_ = static_cast<int>(unique_faces.size());

      // Initialize per-DOF summaries and previous values
      summaries_.resize(traced_dofs_.size());
      prev_V_mag_.resize(traced_dofs_.size(), -1.0);

      // Open output files (fail-soft)
      std::string suffix = "r" + std::to_string(rank_);
      std::string dir = cfg_.output_dir;

      f_meta_.open(dir + "/trace_faces_" + suffix + ".csv");
      f_ts_.open(dir + "/trace_timeseries_" + suffix + ".csv");
      f_events_.open(dir + "/trace_events_" + suffix + ".csv");
      summary_path_ = dir + "/trace_summary_" + suffix + ".csv";

      if (!f_meta_.is_open() || !f_ts_.is_open() || !f_events_.is_open())
      {
         Disable();
         return;
      }

      // Write per-DOF metadata CSV
      f_meta_ << "rank,fi,dof,owned_dof_idx,role,cx,cz,nbf,face_type,"
              << "depth,a,Dc,tau_pre_dip,tau_pre_strike,"
              << "V_init_dip,V_init_strike\n";
      for (auto &s : traced_dofs_)
      {
         f_meta_ << rank_ << "," << s.fi << "," << s.dof << ","
                 << s.owned_dof_idx << "," << s.role << ","
                 << std::scientific << std::setprecision(6)
                 << s.cx << "," << s.cz << "," << nbf_per_face_ << ","
                 << s.face_type << ","
                 << s.depth << "," << s.a << "," << s.Dc << ","
                 << s.tau_pre_dip << "," << s.tau_pre_strike << ","
                 << s.V_init_dip << "," << s.V_init_strike << "\n";
      }
      f_meta_.flush();

      // Write CSV headers
      f_ts_ << "step,time_s,time_yr,rank,fi,dof,owned_dof_idx,role,cx,cz,"
            << "tau_dip,tau_strike,"
            << "tau_stress_dip,tau_stress_strike,"
            << "tau_corr_dip,tau_corr_strike,"
            << "jump_res_dip,jump_res_strike,"
            << "normal_traction,normal_stress,normal_corr,sigma_n_eff,"
            << "V_dip,V_strike,V_mag,psi\n";
      f_ts_.flush();

      f_events_ << "step,time_s,time_yr,rank,fi,dof,owned_dof_idx,"
                << "event,value,threshold,"
                << "sigma_n_eff,normal_traction,normal_stress,normal_corr,"
                << "V_mag\n";
      f_events_.flush();
   }

   bool IsActive() const { return active_ && !disabled_; }

   /// Stage decomposition data from Mult() (overwritten by rejected ODE stages).
   void RecordMult(
      const Vector &traction,
      const Vector &traction_stress,
      const Vector &traction_correction,
      const Vector &jump_residual,
      const Vector &normal_traction,
      const Vector &normal_stress,
      const Vector &normal_correction,
      const Vector &slip_rate,
      real_t sigma_n_base)
   {
      if (!IsActive()) { return; }

      // R-202: the tracer's CommitStep indexes staged_traction_ /
      // staged_stress_ / staged_corr_ / staged_jump_res_ with BP5 2-comp
      // strides (2*dof and 2*dof+1).  If a future caller wires Antiplane
      // (NumSlipComponents == 1) into here, the round-1 R-003 base
      // default sizes those vectors to traction.Size() — which is half
      // the stride — and CommitStep would OOB-read the second half on
      // every iteration.  Refuse loudly rather than corrupt silently.
      const int n_owned = static_cast<int>(traced_dofs_.size());
      MFEM_VERIFY(n_owned == 0 || traction.Size() >= 2 * n_owned,
                  "FaceTraceLogger::RecordMult: traction.Size()="
                  << traction.Size() << " < 2 * n_owned=" << 2*n_owned
                  << ".  Tracer requires BP5 2-component layout; "
                  "Antiplane (1 comp) callers must template-specialise "
                  "the tracer or skip tracer wiring.");

      staged_traction_ = traction;
      staged_stress_ = traction_stress;
      staged_corr_ = traction_correction;
      staged_jump_res_ = jump_residual;
      if (normal_traction.Size() > 0)
      {
         staged_normal_trac_ = normal_traction;
      }
      if (normal_stress.Size() > 0)
      {
         staged_normal_stress_ = normal_stress;
      }
      if (normal_correction.Size() > 0)
      {
         staged_normal_corr_ = normal_correction;
      }
      staged_slip_rate_ = slip_rate;
      staged_sigma_n_base_ = sigma_n_base;
      has_staged_data_ = true;
   }

   void CommitStep(int step, real_t time_s, const Vector &state, bool is_final)
   {
      if (!IsActive() || !has_staged_data_) { return; }

      // Suppress duplicate commit of the same (step, time)
      if (step == last_committed_step_ && time_s == last_committed_time_)
      {
         return;
      }
      last_committed_step_ = step;
      last_committed_time_ = time_s;

      bool sample = ShouldSample(step, is_final);
      bool any_event = false;

      const real_t secs_per_yr = 365.25 * 24.0 * 3600.0;
      real_t time_yr = time_s / secs_per_yr;

      // State layout for BP5: [slip_dip, slip_strike, psi] per node
      constexpr int spn = 3;
      constexpr int psi_idx = 2;

      for (size_t ti = 0; ti < traced_dofs_.size(); ti++)
      {
         auto &spec = traced_dofs_[ti];
         auto &summary = summaries_[ti];
         summary.last_time_s = time_s;

         int dof = spec.owned_dof_idx;

         real_t tau_d = staged_traction_(2 * dof);
         real_t tau_s = staged_traction_(2 * dof + 1);
         real_t ts_d = staged_stress_(2 * dof);
         real_t ts_s = staged_stress_(2 * dof + 1);
         real_t tc_d = staged_corr_(2 * dof);
         real_t tc_s = staged_corr_(2 * dof + 1);
         real_t jr_d = staged_jump_res_(2 * dof);
         real_t jr_s = staged_jump_res_(2 * dof + 1);
         real_t nt = (staged_normal_trac_.Size() > dof)
                     ? staged_normal_trac_(dof) : 0.0;
         real_t ns = (staged_normal_stress_.Size() > dof)
                     ? staged_normal_stress_(dof) : 0.0;
         real_t nc = (staged_normal_corr_.Size() > dof)
                     ? staged_normal_corr_(dof) : 0.0;
         real_t sigma_n_eff = staged_sigma_n_base_ + nt;
         real_t V_d = staged_slip_rate_(2 * dof);
         real_t V_s = staged_slip_rate_(2 * dof + 1);
         real_t V_mag = std::sqrt(V_d * V_d + V_s * V_s);
         real_t psi = (dof < state.Size() / spn)
                      ? state(dof * spn + psi_idx) : 0.0;

         // Update summary
         real_t tau_mag = std::sqrt(tau_d * tau_d + tau_s * tau_s);
         summary.tau_max = std::max(summary.tau_max, tau_mag);
         summary.V_max = std::max(summary.V_max, V_mag);
         summary.psi_min = std::min(summary.psi_min, psi);
         summary.sigma_n_eff_min = std::min(summary.sigma_n_eff_min,
                                            sigma_n_eff);
         summary.max_normal_corr = std::max(summary.max_normal_corr,
                                            std::abs(nc));

         // Event detection
         bool force_row = false;

         auto emit_event = [&](const char *name, real_t val, real_t thresh)
         {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(8);
            oss << step << "," << time_s << "," << time_yr << ","
                << rank_ << "," << spec.fi << "," << spec.dof << ","
                << spec.owned_dof_idx << ","
                << name << "," << val << "," << thresh << ","
                << sigma_n_eff << "," << nt << "," << ns << "," << nc << ","
                << V_mag << "\n";
            event_buffer_.push_back(oss.str());
            force_row = true;
            any_event = true;
            if (summary.first_event.empty())
            {
               summary.first_event = name;
            }
         };

         // Non-finite check
         if (!std::isfinite(tau_d) || !std::isfinite(tau_s) ||
             !std::isfinite(V_mag) || !std::isfinite(psi) ||
             !std::isfinite(sigma_n_eff))
         {
            emit_event("non_finite", 0.0, 0.0);
         }

         // Normal-collapse graduated thresholds
         if (sigma_n_eff <= 0.0)
         {
            emit_event("sigma_n_eff_nonpos", sigma_n_eff, 0.0);
         }
         else if (sigma_n_eff < 1e6)
         {
            emit_event("sigma_n_eff_lt_1mpa", sigma_n_eff, 1e6);
         }
         else if (sigma_n_eff < 5e6)
         {
            emit_event("sigma_n_eff_lt_5mpa", sigma_n_eff, 5e6);
         }

         // Normal correction dominates: |normal_corr| > |normal_stress|
         if (std::abs(ns) > 1e-30 && std::abs(nc) > std::abs(ns))
         {
            emit_event("normal_corr_dominates", std::abs(nc) / std::abs(ns),
                       1.0);
         }

         // V_mag high
         if (V_mag > cfg_.V_mag_warn)
         {
            emit_event("V_mag_high", V_mag, cfg_.V_mag_warn);
         }

         // Correction ratio (tangential)
         real_t tc_mag = std::sqrt(tc_d * tc_d + tc_s * tc_s);
         real_t corr_ratio = tc_mag / std::max(tau_mag, 1e-30);
         if (corr_ratio > cfg_.corr_ratio_warn)
         {
            emit_event("corr_ratio_high", corr_ratio,
                       cfg_.corr_ratio_warn);
         }

         // Sudden jump in V_mag
         if (prev_V_mag_[ti] > 0.0 && V_mag > 0.0)
         {
            real_t ratio = V_mag / prev_V_mag_[ti];
            if (ratio > cfg_.jump_factor_warn)
            {
               emit_event("V_jump", ratio, cfg_.jump_factor_warn);
            }
         }
         prev_V_mag_[ti] = V_mag;

         // Write timeseries row if sampled or forced by event
         if (sample || force_row)
         {
            std::ostringstream oss;
            oss << std::scientific << std::setprecision(8);
            oss << step << "," << time_s << "," << time_yr << ","
                << rank_ << "," << spec.fi << "," << spec.dof << ","
                << spec.owned_dof_idx << "," << spec.role << ","
                << spec.cx << "," << spec.cz << ","
                << tau_d << "," << tau_s << ","
                << ts_d << "," << ts_s << ","
                << tc_d << "," << tc_s << ","
                << jr_d << "," << jr_s << ","
                << nt << "," << ns << "," << nc << "," << sigma_n_eff << ","
                << V_d << "," << V_s << "," << V_mag << ","
                << psi << "\n";
            ts_buffer_.push_back(oss.str());
         }
      }

      rows_since_flush_++;
      if (rows_since_flush_ >= cfg_.flush_interval || any_event)
      {
         FlushBuffers();
      }
   }

   void Finalize()
   {
      if (!active_ || finalized_) { return; }
      finalized_ = true;

      FlushBuffers();
      WriteSummary();  // final definitive write

      if (f_meta_.is_open()) { f_meta_.close(); }
      if (f_ts_.is_open()) { f_ts_.close(); }
      if (f_events_.is_open()) { f_events_.close(); }
   }

   void Disable()
   {
      disabled_ = true;
      active_ = false;
   }

private:
   bool ShouldSample(int step, bool is_final) const
   {
      if (is_final) { return true; }
      if (step <= cfg_.warmup_steps) { return true; }
      if (step > 0 && (step & (step - 1)) == 0) { return true; }
      if (step % cfg_.periodic_interval == 0) { return true; }
      return false;
   }

   void FlushBuffers()
   {
      if (f_ts_.is_open())
      {
         for (auto &row : ts_buffer_) { f_ts_ << row; }
         f_ts_.flush();
      }
      ts_buffer_.clear();

      if (f_events_.is_open())
      {
         for (auto &row : event_buffer_) { f_events_ << row; }
         f_events_.flush();
      }
      event_buffer_.clear();

      rows_since_flush_ = 0;

      // Incremental summary: rewrite on every flush so the file is
      // always up-to-date even if the run aborts.
      WriteSummary();
   }

   void WriteSummary()
   {
      if (summary_path_.empty()) { return; }

      std::ofstream f(summary_path_);
      if (!f.is_open()) { return; }

      f << "rank,fi,dof,owned_dof_idx,role,"
        << "tau_max,V_max,psi_min,sigma_n_eff_min,max_normal_corr,"
        << "first_event,last_time_s\n";

      for (size_t ti = 0; ti < traced_dofs_.size(); ti++)
      {
         auto &spec = traced_dofs_[ti];
         auto &s = summaries_[ti];
         f << std::scientific << std::setprecision(8);
         f << rank_ << "," << spec.fi << "," << spec.dof << ","
           << spec.owned_dof_idx << "," << spec.role << ","
           << s.tau_max << "," << s.V_max << ","
           << s.psi_min << "," << s.sigma_n_eff_min << ","
           << s.max_normal_corr << ","
           << (s.first_event.empty() ? "none" : s.first_event) << ","
           << s.last_time_s << "\n";
      }
      f.flush();
      f.close();
   }

   TraceConfig cfg_;
   int rank_;
   bool active_ = false;
   bool disabled_ = false;
   bool finalized_ = false;
   bool has_staged_data_ = false;
   int rows_since_flush_ = 0;
   int last_committed_step_ = -1;
   real_t last_committed_time_ = -1.0;

   std::vector<TraceFaceSpec> traced_dofs_;
   int nbf_per_face_ = 1;
   int num_traced_faces_ = 0;

   // Staging area (overwritten each Mult, committed on accepted step)
   Vector staged_traction_, staged_stress_, staged_corr_;
   Vector staged_jump_res_, staged_normal_trac_;
   Vector staged_normal_stress_, staged_normal_corr_;
   Vector staged_slip_rate_;
   real_t staged_sigma_n_base_ = 0.0;

   // Buffered CSV rows
   std::vector<std::string> ts_buffer_, event_buffer_;

   // Per-DOF summary accumulators
   struct DofSummary
   {
      real_t tau_max = 0, V_max = 0, psi_min = 1e30;
      real_t sigma_n_eff_min = 1e30;
      real_t max_normal_corr = 0;
      std::string first_event;
      real_t last_time_s = 0;
   };
   std::vector<DofSummary> summaries_;

   // Previous V_mag per traced DOF for jump detection
   std::vector<real_t> prev_V_mag_;

   // File handles
   std::ofstream f_meta_, f_ts_, f_events_;
   std::string summary_path_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FACE_TRACE_LOGGER_HPP
