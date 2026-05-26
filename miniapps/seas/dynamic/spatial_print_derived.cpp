// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/spatial_print_derived.cpp — Phase D implementation.

#include "spatial_print_derived.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

namespace mfem
{
namespace seas
{
namespace spatial
{

namespace
{

constexpr real_t kBarrierSentinel        = 1.0e6;
constexpr real_t kBarrierHalfThreshold   = 0.5 * kBarrierSentinel;

bool EnvSkipEquilibriumGate()
{
   const char* v = std::getenv("SEAS_SKIP_EQUILIBRIUM_GATE");
   return v != nullptr && std::string(v) == "1";
}

#ifdef MFEM_USE_MPI
real_t AllreduceMin(real_t local, MPI_Comm comm)
{
   real_t global = local;
   MPI_Allreduce(&local, &global, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
   return global;
}
real_t AllreduceMax(real_t local, MPI_Comm comm)
{
   real_t global = local;
   MPI_Allreduce(&local, &global, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   return global;
}
real_t AllreduceSum(real_t local, MPI_Comm comm)
{
   real_t global = local;
   MPI_Allreduce(&local, &global, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   return global;
}
int AllreduceSumInt(int local, MPI_Comm comm)
{
   int global = local;
   MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, comm);
   return global;
}
#endif

/// Compute σ_1 eigen-direction + azimuth/plunge from the 6-component
/// symmetric stress tensor in EAST-NORTH-UP frame (x = E, y = N).
/// Azimuth is clockwise from north in degrees (R-004 fix: atan2(x, y)
/// ordering, NOT atan2(y, x)).
void ComputeSigma1AzimuthPlungeMagnitude(const StressSpec& s,
                                         real_t& azimuth_deg,
                                         real_t& plunge_deg,
                                         real_t& magnitude_pa)
{
   DenseMatrix S(3, 3);
   S(0,0) = s.sigma_xx_pa;
   S(1,1) = s.sigma_yy_pa;
   S(2,2) = s.sigma_zz_pa;
   S(0,1) = S(1,0) = s.sigma_xy_pa;
   S(1,2) = S(2,1) = s.sigma_yz_pa;
   S(0,2) = S(2,0) = s.sigma_xz_pa;

   Vector eigvals(3);
   DenseMatrix eigvecs(3, 3);
   S.Eigensystem(eigvals, eigvecs);   // ascending order

   const int i1 = 2;                  // σ_1 = largest (most compressive)
   magnitude_pa  = eigvals(i1);
   const real_t v1x = eigvecs(0, i1);
   const real_t v1y = eigvecs(1, i1);
   const real_t v1z = eigvecs(2, i1);

   // R-004 fix: (x, y) atan2 ordering for clockwise-from-N azimuth.
   azimuth_deg = std::atan2(v1x, v1y) * 180.0 / M_PI;
   if (azimuth_deg < 0.0) { azimuth_deg += 360.0; }

   const real_t v1h = std::sqrt(v1x * v1x + v1y * v1y);
   plunge_deg = std::atan2(-v1z, v1h) * 180.0 / M_PI;
}

/// Regularized rate-and-state steady-state friction coefficient at slip
/// rate V, grounded in this repo's own kernels (no Tandem):
///   ψ_ss = f_0 + b · ln(V_0 / V)               (state_evolution.hpp:196-200)
///   f_ss = a · asinh[(V / 2V_0) · exp(ψ_ss/a)]  (dieterich_ruina.hpp:272-295)
/// with the identical large-(ψ/a) guard the kernel uses
/// (asinh(x) ≈ ln(2x) ⇒ f ≈ a·ln(V/V_0) + ψ) to avoid exp overflow.
real_t SteadyStateFrictionRS(real_t a, real_t b, real_t f0,
                             real_t V0, real_t V)
{
   if (V <= 0.0 || V0 <= 0.0 || a <= 0.0) { return 0.0; }
   const real_t psi_ss     = f0 + b * std::log(V0 / V);
   const real_t psi_over_a = psi_ss / a;
   if (psi_over_a > 700.0)
   {
      return std::max(static_cast<real_t>(0.0),
                      a * std::log(V / V0) + psi_ss);
   }
   const real_t sinh_arg = (V / (2.0 * V0)) * std::exp(psi_over_a);
   return a * std::asinh(sinh_arg);
}

}  // namespace

real_t PrintDerivedAndCheck(
   const PrintDerivedConfig&                cfg,
   const SlipWeakeningPerDOFParams&         lsw,
   const Vector&                            tau_pre_per_dof,
   const Vector&                            sigma_n_eff_per_dof,
   const Vector&                            dof_coords_3d,
   const NucleationSpec&                    nuc,
   const GradualOverstressPerDOFParams&     nuc_params,
   const StressSpec&                        stress,
   real_t                                   mu_bulk,
   real_t                                   cp,
   real_t                                   cs,
   real_t                                   h_min_global,
   real_t                                   dt_cfl,
   real_t                                   tfinal,
   int                                      num_fault_global,
   int                                      num_zero_normal_fallbacks
#ifdef MFEM_USE_MPI
   , MPI_Comm comm
#endif
   , int                                    rank,
   std::ostream&                            out)
{
   if (!cfg.enabled) { return 0.0; }

   const int N = sigma_n_eff_per_dof.Size();
   MFEM_VERIFY(lsw.mu_s.Size() == N, "PrintDerivedAndCheck: lsw size mismatch");
   MFEM_VERIFY(tau_pre_per_dof.Size() == 2 * N,
               "PrintDerivedAndCheck: tau_pre layout mismatch (expected 2N)");
   MFEM_VERIFY(dof_coords_3d.Size() == 3 * N,
               "PrintDerivedAndCheck: dof_coords_3d layout mismatch");

   const bool skip_gate = EnvSkipEquilibriumGate();
   const bool warn_only = skip_gate || !cfg.abort_on_failure;

   // -----------------------------------------------------------------
   // 1.  Per-DOF L_nuc on the non-barrier DOFs.
   // -----------------------------------------------------------------
   real_t Lnuc_local_min = std::numeric_limits<real_t>::infinity();
   real_t Lnuc_local_max = -std::numeric_limits<real_t>::infinity();
   real_t Lnuc_local_sum = 0.0;
   int    Lnuc_local_cnt = 0;
   for (int i = 0; i < N; ++i)
   {
      if (lsw.mu_s(i) >= kBarrierHalfThreshold) { continue; }
      const real_t denom = (lsw.mu_s(i) - lsw.mu_d(i)) * sigma_n_eff_per_dof(i);
      if (denom <= 0.0) { continue; }
      const real_t Lnuc_i = mu_bulk * lsw.d_c(i) / denom;
      Lnuc_local_min = std::min(Lnuc_local_min, Lnuc_i);
      Lnuc_local_max = std::max(Lnuc_local_max, Lnuc_i);
      Lnuc_local_sum += Lnuc_i;
      Lnuc_local_cnt += 1;
   }

#ifdef MFEM_USE_MPI
   const real_t Lnuc_min  = AllreduceMin(Lnuc_local_min, comm);
   const real_t Lnuc_max  = AllreduceMax(Lnuc_local_max, comm);
   const real_t Lnuc_sum  = AllreduceSum(Lnuc_local_sum, comm);
   const int    Lnuc_cnt  = AllreduceSumInt(Lnuc_local_cnt, comm);
#else
   const real_t Lnuc_min = Lnuc_local_min;
   const real_t Lnuc_max = Lnuc_local_max;
   const real_t Lnuc_sum = Lnuc_local_sum;
   const int    Lnuc_cnt = Lnuc_local_cnt;
#endif
   const real_t Lnuc_mean = (Lnuc_cnt > 0) ? (Lnuc_sum / Lnuc_cnt) : 0.0;
   const bool   all_barriers = (Lnuc_cnt == 0);

   // -----------------------------------------------------------------
   // 2.  Outside-asperity equilibrium ratio.
   //     Use Euclidean distance from the nucleation centre (consistent
   //     with the plan's spec).  When nucleation is disabled, every DOF
   //     counts as "outside" — the gate then asks whether ANY DOF
   //     exceeds the strength.
   // -----------------------------------------------------------------
   real_t outside_max_local = 0.0;
   int    outside_argmax_local = -1;
   real_t r_threshold = 0.0;
   if (nuc.enabled)
   {
      r_threshold = cfg.outside_safety_factor *
                    std::max(nuc.gradual_overstress.radius_dip_m,
                             nuc.gradual_overstress.radius_strike_m);
   }
   const real_t cx = nuc.gradual_overstress.center_x_m;
   const real_t cy = nuc.gradual_overstress.center_y_m;
   const real_t cz = nuc.gradual_overstress.center_z_m;
   for (int i = 0; i < N; ++i)
   {
      if (lsw.mu_s(i) >= kBarrierHalfThreshold) { continue; }
      const real_t dx = dof_coords_3d(3*i + 0) - cx;
      const real_t dy = dof_coords_3d(3*i + 1) - cy;
      const real_t dz = dof_coords_3d(3*i + 2) - cz;
      const real_t r = std::sqrt(dx*dx + dy*dy + dz*dz);
      const bool outside = !nuc.enabled || (r > r_threshold);
      if (!outside) { continue; }
      const real_t tau1 = tau_pre_per_dof(2*i + 0);
      const real_t tau2 = tau_pre_per_dof(2*i + 1);
      const real_t tau_mag = std::sqrt(tau1*tau1 + tau2*tau2);
      const real_t strength = lsw.mu_s(i) * sigma_n_eff_per_dof(i);
      if (strength <= 0.0) { continue; }
      const real_t ratio = tau_mag / strength;
      if (ratio > outside_max_local)
      {
         outside_max_local    = ratio;
         outside_argmax_local = i;
      }
   }

#ifdef MFEM_USE_MPI
   struct { double value; int rank_; } in{outside_max_local, rank}, out_maxloc;
   MPI_Allreduce(&in, &out_maxloc, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
   const real_t outside_max = static_cast<real_t>(out_maxloc.value);
   const int    owner       = out_maxloc.rank_;
   real_t outside_coord[3]  = { 0.0, 0.0, 0.0 };
   if (owner == rank && outside_argmax_local >= 0)
   {
      outside_coord[0] = dof_coords_3d(3 * outside_argmax_local + 0);
      outside_coord[1] = dof_coords_3d(3 * outside_argmax_local + 1);
      outside_coord[2] = dof_coords_3d(3 * outside_argmax_local + 2);
   }
   MPI_Bcast(outside_coord, 3, MPITypeMap<real_t>::mpi_type, owner, comm);
#else
   const real_t outside_max = outside_max_local;
   const int    owner       = 0;
   real_t outside_coord[3]  = { 0.0, 0.0, 0.0 };
   if (outside_argmax_local >= 0)
   {
      outside_coord[0] = dof_coords_3d(3 * outside_argmax_local + 0);
      outside_coord[1] = dof_coords_3d(3 * outside_argmax_local + 1);
      outside_coord[2] = dof_coords_3d(3 * outside_argmax_local + 2);
   }
#endif

   // -----------------------------------------------------------------
   // 3.  Nucleation well-posedness, over DOFs INSIDE the nucleation
   //     region (r <= outside_safety_factor · max(radius_*)).  Three
   //     LSW forced-nucleation conditions:
   //       (1) TRIGGER     : |tau_pre + delta_tau| >= mu_s · sigma_n_eff
   //                         at the peak DOF (forcing reaches static yield).
   //                         overshoot := |tau_pre + delta_tau| - mu_s·sigma_n.
   //       (2) LOCKED      : |tau_pre| <  mu_s · sigma_n_eff everywhere
   //                         in-patch (not already failing without forcing).
   //       (3) STRESS-DROP : |tau_pre| >  mu_d · sigma_n_eff everywhere
   //                         in-patch (positive dynamic stress drop, so the
   //                         rupture is self-sustaining and does not re-lock
   //                         when the forcing ends).
   //     delta_tau is the VECTOR overstress (dip, strike); the trigger
   //     compares the magnitude of (pre-stress + overstress), so an
   //     overstress that is not parallel to tau_pre is handled correctly.
   // -----------------------------------------------------------------
   real_t overshoot_local_best   = -std::numeric_limits<real_t>::infinity();
   int    overshoot_local_arg    = -1;
   real_t nuc_peak_local         = 0.0;
   real_t static_ratio_local_max = 0.0;                                     // (2)
   real_t dyn_ratio_local_min    = std::numeric_limits<real_t>::infinity(); // (3)
   if (nuc.enabled && nuc_params.amplitude_dip.Size() == N)
   {
      for (int i = 0; i < N; ++i)
      {
         if (lsw.mu_s(i) >= kBarrierHalfThreshold) { continue; }
         const real_t dx = dof_coords_3d(3*i + 0) - cx;
         const real_t dy = dof_coords_3d(3*i + 1) - cy;
         const real_t dz = dof_coords_3d(3*i + 2) - cz;
         const real_t r  = std::sqrt(dx*dx + dy*dy + dz*dz);
         if (r > r_threshold) { continue; }
         const real_t tau1 = tau_pre_per_dof(2*i + 0);
         const real_t tau2 = tau_pre_per_dof(2*i + 1);
         const real_t ad   = nuc_params.amplitude_dip(i);
         const real_t as   = nuc_params.amplitude_strike(i);
         const real_t A    = std::sqrt(ad*ad + as*as);
         const real_t tau_pre_mag = std::sqrt(tau1*tau1 + tau2*tau2);
         const real_t strength_s  = lsw.mu_s(i) * sigma_n_eff_per_dof(i);
         const real_t strength_d  = lsw.mu_d(i) * sigma_n_eff_per_dof(i);
         // (1) nucleated traction = pre-stress + overstress (vector sum).
         const real_t tau_nuc1    = tau1 + ad;
         const real_t tau_nuc2    = tau2 + as;
         const real_t tau_nuc_mag = std::sqrt(tau_nuc1*tau_nuc1 + tau_nuc2*tau_nuc2);
         const real_t overshoot   = tau_nuc_mag - strength_s;
         if (overshoot > overshoot_local_best)
         {
            overshoot_local_best = overshoot;
            overshoot_local_arg  = i;
         }
         nuc_peak_local = std::max(nuc_peak_local, A);
         // (2) pre-stress static ratio: want < 1 everywhere in-patch.
         if (strength_s > 0.0)
         {
            static_ratio_local_max =
               std::max(static_ratio_local_max, tau_pre_mag / strength_s);
         }
         // (3) pre-stress dynamic ratio: want > 1 everywhere in-patch.
         //     mu_d == 0 ⇒ zero residual strength ⇒ always positive drop.
         const real_t dyn_ratio = (strength_d > 0.0)
                                  ? (tau_pre_mag / strength_d)
                                  : std::numeric_limits<real_t>::infinity();
         dyn_ratio_local_min = std::min(dyn_ratio_local_min, dyn_ratio);
      }
   }
#ifdef MFEM_USE_MPI
   const real_t nuc_peak = nuc.enabled ? AllreduceMax(nuc_peak_local, comm)
                                       : 0.0;
   struct { double value; int rank_; } pin{overshoot_local_best, rank}, pout;
   MPI_Allreduce(&pin, &pout, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
   const real_t overshoot_best  = static_cast<real_t>(pout.value);
   const int    overshoot_owner = pout.rank_;
   const real_t inpatch_static_ratio_max =
      nuc.enabled ? AllreduceMax(static_ratio_local_max, comm) : 0.0;
   const real_t inpatch_dyn_ratio_min =
      nuc.enabled ? AllreduceMin(dyn_ratio_local_min, comm)
                  : std::numeric_limits<real_t>::infinity();
#else
   const real_t nuc_peak                 = nuc_peak_local;
   const real_t overshoot_best           = overshoot_local_best;
   const int    overshoot_owner          = 0;
   const real_t inpatch_static_ratio_max = static_ratio_local_max;
   const real_t inpatch_dyn_ratio_min    = dyn_ratio_local_min;
#endif

   // Compute the amplitude + coords at the overshoot-argmax (peak) DOF for
   // the print.  Also detect the "no DOF inside nucleation patch on any
   // rank" case so we don't silently PASS a config that can't nucleate.
   real_t A_at_argmax      = 0.0;
   real_t coord_at_argmax[3] = { 0.0, 0.0, 0.0 };
   const int local_in_patch_count = (nuc.enabled && overshoot_local_arg >= 0)
                                     ? 1 : 0;
   if (nuc.enabled && overshoot_owner == rank && overshoot_local_arg >= 0)
   {
      const int i = overshoot_local_arg;
      const real_t ad = nuc_params.amplitude_dip(i);
      const real_t as = nuc_params.amplitude_strike(i);
      A_at_argmax = std::sqrt(ad*ad + as*as);
      coord_at_argmax[0] = dof_coords_3d(3*i + 0);
      coord_at_argmax[1] = dof_coords_3d(3*i + 1);
      coord_at_argmax[2] = dof_coords_3d(3*i + 2);
   }
#ifdef MFEM_USE_MPI
   int patch_count_global = local_in_patch_count;
   if (nuc.enabled)
   {
      MPI_Allreduce(&local_in_patch_count, &patch_count_global, 1, MPI_INT,
                    MPI_SUM, comm);
      MPI_Bcast(&A_at_argmax,      1, MPITypeMap<real_t>::mpi_type,
                overshoot_owner, comm);
      MPI_Bcast(coord_at_argmax,   3, MPITypeMap<real_t>::mpi_type,
                overshoot_owner, comm);
   }
#else
   const int patch_count_global = local_in_patch_count;
#endif
   const bool nuc_patch_empty = nuc.enabled && (patch_count_global == 0);

   // -----------------------------------------------------------------
   // 4.  σ_1 azimuth / plunge / magnitude (constant_tensor only).
   // -----------------------------------------------------------------
   real_t s1_azimuth = 0.0, s1_plunge = 0.0, s1_mag = 0.0;
   const bool stress_const =
      (stress.kind == StressSourceKind::ConstantTensor);
   if (stress_const)
   {
      ComputeSigma1AzimuthPlungeMagnitude(stress, s1_azimuth, s1_plunge, s1_mag);
   }

   // -----------------------------------------------------------------
   // 5.  Rank-0 print.
   // -----------------------------------------------------------------
   const int nsteps = (dt_cfl > 0.0)
                       ? static_cast<int>(std::ceil(tfinal / dt_cfl)) : 0;
   const real_t Lnuc_h_min = (h_min_global > 0.0) ? (Lnuc_min / h_min_global) : 0.0;
   const real_t Lnuc_h_max = (h_min_global > 0.0) ? (Lnuc_max / h_min_global) : 0.0;

   if (rank == 0)
   {
      out << "[derived] num_fault_global = " << num_fault_global << "\n";
      out << "[derived] cp = " << cp << " m/s, cs = " << cs << " m/s\n";
      out << "[derived] dt_cfl = " << dt_cfl << " s, nsteps = "
          << nsteps << "\n";
      out << "[derived] mesh h_min = " << h_min_global << " m (scalar cp)\n";
      if (all_barriers)
      {
         out << "[derived] all fault DOFs are barriers; nothing will rupture\n";
      }
      else
      {
         out << "[derived] L_nuc min/max/mean = "
             << Lnuc_min << " / " << Lnuc_max << " / " << Lnuc_mean << " m\n";
         out << "[derived] L_nuc / h_min min/max = "
             << Lnuc_h_min << " / " << Lnuc_h_max
             << " (NB: should be >= 10 in nucleation patch)\n";
      }
      if (nuc.enabled)
      {
         out << "[derived] |tau_pre| / (mu_s * sigma_n_eff) max OUTSIDE "
                "asperity = " << outside_max << "\n";
      }
      else
      {
         out << "[derived] note: no [nucleation] block; rupture must be "
                "driven entirely by tau_pre exceeding strength somewhere.\n";
         out << "[derived] |tau_pre| / (mu_s * sigma_n_eff) max anywhere = "
             << outside_max << "\n";
      }
      out << "[derived] nucleation: "
          << (nuc.enabled ? "enabled" : "disabled") << "\n";
      if (nuc.enabled)
      {
         out << "[derived] nucleation peak |F(r) * delta_tau| = "
             << nuc_peak << " Pa\n";
         if (nuc_patch_empty)
         {
            out << "[derived] WARNING: zero fault DOFs are inside the "
                   "nucleation patch (r < 3 · max(radius_*)).  Check that "
                   "the nucleation centre lies on the fault and that the "
                   "mesh resolves it.\n";
         }
         else
         {
            out << "[derived] most-overstressed DOF at (x="
                << coord_at_argmax[0] << ", y=" << coord_at_argmax[1]
                << ", z=" << coord_at_argmax[2] << ")\n";
            out << "[derived] amplitude |delta_tau| at that DOF = "
                << A_at_argmax << " Pa\n";
            out << "[derived] nucleation overshoot "
                   "(|tau_pre + delta_tau| - mu_s*sigma_n_eff) = "
                << (overshoot_best / 1.0e6) << " MPa "
                << ((overshoot_best >= 0.0) ? "(sufficient)" : "(INSUFFICIENT)")
                << "\n";
            out << "[derived] in-patch max |tau_pre| / (mu_s*sigma_n_eff) = "
                << inpatch_static_ratio_max
                << " (should be < 1: patch initially locked)\n";
            out << "[derived] in-patch min |tau_pre| / (mu_d*sigma_n_eff) = "
                << inpatch_dyn_ratio_min
                << " (should be > 1: positive dynamic stress drop)\n";
         }
      }
      if (stress_const)
      {
         out << "[derived] sigma_1 azimuth = " << s1_azimuth << " deg, "
                "plunge = " << s1_plunge << " deg, magnitude = "
             << s1_mag << " Pa\n";
      }
      else
      {
         out << "[derived] sigma_1 azimuth from sidecar: TBD "
                "(depth-dependent; see stress sidecar README) — R-002\n";
      }
      out << "[derived] fault DOFs with zero-normal fallback: "
          << num_zero_normal_fallbacks << "\n";
   }

   // -----------------------------------------------------------------
   // 6.  Gating.  Aborts (or warns) on:
   //    (a) all barriers,
   //    (b) outside_max >= 1.0 (background supercritical),
   //    (c) nuc enabled and TRIGGER fails: overshoot < 0, i.e. the peak
   //        |tau_pre + delta_tau| never reaches mu_s·sigma_n_eff,
   //    (d) nuc enabled and STRESS-DROP fails: in-patch min
   //        |tau_pre|/(mu_d·sigma_n_eff) <= 1 (negative dynamic stress drop),
   //    (e) nuc disabled and outside_max < 1.0 (cannot nucleate at all).
   //  A WARNING (never an abort) also fires if the patch is already at
   //  static yield before forcing (in-patch max |tau_pre|/(mu_s·sigma_n) >= 1).
   // -----------------------------------------------------------------
   auto fail = [&](const std::string& msg)
   {
      if (warn_only)
      {
         if (rank == 0)
         {
            // R-003 fix: attribute the correct reason — env var vs
            // API caller's abort_on_failure choice — so operator
            // logs don't mislead future readers.
            const char* reason = skip_gate
               ? "SEAS_SKIP_EQUILIBRIUM_GATE=1"
               : "cfg.abort_on_failure=false";
            out << "[derived] WARNING (gate downgraded via " << reason
                << "): " << msg << "\n";
         }
      }
      else
      {
         MFEM_ABORT(msg);
      }
   };

   if (all_barriers)
   {
      fail("[derived] all fault DOFs are barriers; nothing will rupture.");
   }
   else if (outside_max >= 1.0)
   {
      std::ostringstream m;
      m << "[derived] FAIL: max OUTSIDE asperity = " << outside_max
        << " >= 1.0 — pre-stress already exceeds friction strength outside "
           "the nucleation zone at DOF (x=" << outside_coord[0]
        << ", y=" << outside_coord[1] << ", z=" << outside_coord[2]
        << "; owner rank " << owner << ").  "
           "Reduce |tau_pre| or increase mu_s.";
      fail(m.str());
   }
   else if (outside_max >= 0.9)
   {
      if (rank == 0)
      {
         out << "[derived] WARNING: max OUTSIDE ratio = " << outside_max
             << " (close to the 1.0 failure threshold).\n";
      }
   }

   if (nuc.enabled && !all_barriers)
   {
      if (nuc_patch_empty)
      {
         fail("[derived] FAIL: nucleation enabled but ZERO fault DOFs are "
              "inside the patch (3·max(radius_*)).  Check that the "
              "nucleation centre is on the fault.");
      }
      else
      {
         // (c) TRIGGER: the peak in-patch DOF must reach static yield.
         if (overshoot_best < 0.0)
         {
            std::ostringstream m;
            m << "[derived] FAIL: nucleation overshoot = "
              << (overshoot_best / 1.0e6) << " MPa (INSUFFICIENT) — "
                 "|tau_pre + delta_tau| < mu_s*sigma_n_eff at the peak DOF, so "
                 "the forcing never reaches static yield.  Increase "
                 "delta_tau_*_pa or reduce mu_s in the nucleation patch.";
            fail(m.str());
         }
         // (d) STRESS-DROP: positive dynamic stress drop everywhere in-patch.
         if (inpatch_dyn_ratio_min <= 1.0)
         {
            std::ostringstream m;
            m << "[derived] FAIL: in-patch min |tau_pre| / (mu_d*sigma_n_eff) = "
              << inpatch_dyn_ratio_min << " <= 1 — NEGATIVE dynamic stress drop "
                 "(tau_pre <= mu_d*sigma_n_eff somewhere in the nucleation "
                 "patch), so the rupture re-locks once the forcing ends.  "
                 "Lower mu_d or raise tau_pre in the nucleation patch.";
            fail(m.str());
         }
         // LOCKED: warn (never abort) if the patch is already at static yield.
         if (inpatch_static_ratio_max >= 1.0 && rank == 0)
         {
            out << "[derived] WARNING: in-patch max |tau_pre| / "
                   "(mu_s*sigma_n_eff) = " << inpatch_static_ratio_max
                << " >= 1 — the nucleation patch is already at static yield "
                   "before forcing; it ruptures spontaneously (the "
                   "gradual_overstress ramp is moot).\n";
         }
      }
   }
   else if (!nuc.enabled)
   {
      if (outside_max < 1.0)
      {
         std::ostringstream m;
         m << "[derived] FAIL: no [nucleation] block AND max |tau_pre| / "
              "(mu_s * sigma_n_eff) = " << outside_max
           << " < 1.0 — this configuration will not produce any rupture.";
         fail(m.str());
      }
   }

   if (num_zero_normal_fallbacks > 0 && rank == 0)
   {
      out << "[derived] WARNING: " << num_zero_normal_fallbacks
          << " fault DOFs hit the zero-normal fallback; "
             "gradual_overstress F(r) at those DOFs may be miscomputed.\n";
   }

   // PASS iff the background is subcritical AND (no nuc, OR all barriers,
   // OR the nucleation is well-posed: triggers AND has a positive dynamic
   // stress drop).  The LOCKED check is warn-only, so it does not gate PASS.
   const bool nuc_well_posed = !nuc_patch_empty
                               && overshoot_best >= 0.0
                               && inpatch_dyn_ratio_min > 1.0;
   if (rank == 0 && !warn_only && outside_max < 1.0
       && (!nuc.enabled || all_barriers || nuc_well_posed))
   {
      out << "[derived] PASS: initial conditions are well-posed.\n";
   }

   return outside_max;
}

// =====================================================================
// Rate-and-state (aging-law) overload.  PLAN DEVIATION (documented in the
// header): the Phase-3 plan does not specify an RS --print-derived path;
// this is added so the SAFS-RS sbatch (which passes --print-derived) runs.
// All RS physics is grounded in this repo's own kernels (no Tandem):
//   L_nuc = μ·Dc / ((b − a)·σ_n)              [RS analog of the LSW formula]
//   f_ss  = a·asinh[(V/2V_0)·exp(ψ_ss/a)],  ψ_ss = f_0 + b·ln(V_0/V)
//                                              [dieterich_ruina + state_evo]
// =====================================================================
real_t PrintDerivedAndCheckRS(
   const PrintDerivedConfig&                cfg,
   const RateStatePerDOFParams&             rs,
   const Vector&                            tau_pre_per_dof,
   const Vector&                            sigma_n_eff_per_dof,
   const Vector&                            dof_coords_3d,
   const NucleationSpec&                    nuc,
   const GradualOverstressPerDOFParams&     nuc_params,
   const StressSpec&                        stress,
   real_t                                   mu_bulk,
   real_t                                   cp,
   real_t                                   cs,
   real_t                                   h_min_global,
   real_t                                   dt_cfl,
   real_t                                   tfinal,
   int                                      num_fault_global,
   int                                      num_zero_normal_fallbacks
#ifdef MFEM_USE_MPI
   , MPI_Comm comm
#endif
   , int                                    rank,
   std::ostream&                            out)
{
   if (!cfg.enabled) { return 0.0; }

   const int N = sigma_n_eff_per_dof.Size();
   MFEM_VERIFY(rs.a.Size() == N && rs.b.Size() == N && rs.Dc.Size() == N &&
               rs.V_init.Size() == N && rs.f_0.Size() == N &&
               rs.V_0.Size() == N,
               "PrintDerivedAndCheckRS: rs per-DOF vector size mismatch "
               "(expected " << N << ")");
   MFEM_VERIFY(tau_pre_per_dof.Size() == 2 * N,
               "PrintDerivedAndCheckRS: tau_pre layout mismatch (expected 2N)");
   MFEM_VERIFY(dof_coords_3d.Size() == 3 * N,
               "PrintDerivedAndCheckRS: dof_coords_3d layout mismatch");

   const bool skip_gate = EnvSkipEquilibriumGate();
   const bool warn_only = skip_gate || !cfg.abort_on_failure;

   // -----------------------------------------------------------------
   // 1.  Per-DOF RS nucleation length (velocity-WEAKENING DOFs only) and
   //     steady-state friction f_ss(V_init).  A DOF is velocity-weakening
   //     iff b > a; b <= a (velocity-strengthening) cannot host an
   //     instability and is the RS analog of an LSW barrier — skipped for
   //     L_nuc.
   // -----------------------------------------------------------------
   real_t Lnuc_local_min = std::numeric_limits<real_t>::infinity();
   real_t Lnuc_local_max = -std::numeric_limits<real_t>::infinity();
   real_t Lnuc_local_sum = 0.0;
   int    Lnuc_local_cnt = 0;          // == velocity-weakening DOF count
   real_t fss_local_min  = std::numeric_limits<real_t>::infinity();
   real_t fss_local_max  = -std::numeric_limits<real_t>::infinity();
   for (int i = 0; i < N; ++i)
   {
      const real_t fss = SteadyStateFrictionRS(rs.a(i), rs.b(i), rs.f_0(i),
                                               rs.V_0(i), rs.V_init(i));
      fss_local_min = std::min(fss_local_min, fss);
      fss_local_max = std::max(fss_local_max, fss);

      const real_t ab = rs.b(i) - rs.a(i);
      if (ab <= 0.0) { continue; }                 // velocity-strengthening
      const real_t denom = ab * sigma_n_eff_per_dof(i);
      if (denom <= 0.0) { continue; }
      const real_t Lnuc_i = mu_bulk * rs.Dc(i) / denom;
      Lnuc_local_min = std::min(Lnuc_local_min, Lnuc_i);
      Lnuc_local_max = std::max(Lnuc_local_max, Lnuc_i);
      Lnuc_local_sum += Lnuc_i;
      Lnuc_local_cnt += 1;
   }

#ifdef MFEM_USE_MPI
   const real_t Lnuc_min = AllreduceMin(Lnuc_local_min, comm);
   const real_t Lnuc_max = AllreduceMax(Lnuc_local_max, comm);
   const real_t Lnuc_sum = AllreduceSum(Lnuc_local_sum, comm);
   const int    vw_cnt   = AllreduceSumInt(Lnuc_local_cnt, comm);
   const real_t fss_min  = AllreduceMin(fss_local_min, comm);
   const real_t fss_max  = AllreduceMax(fss_local_max, comm);
#else
   const real_t Lnuc_min = Lnuc_local_min;
   const real_t Lnuc_max = Lnuc_local_max;
   const real_t Lnuc_sum = Lnuc_local_sum;
   const int    vw_cnt   = Lnuc_local_cnt;
   const real_t fss_min  = fss_local_min;
   const real_t fss_max  = fss_local_max;
#endif
   const real_t Lnuc_mean = (vw_cnt > 0) ? (Lnuc_sum / vw_cnt) : 0.0;
   const bool   all_vs    = (vw_cnt == 0);   // RS analog of all_barriers

   // -----------------------------------------------------------------
   // 2.  Outside-asperity seed stress ratio |τ_pre| / (f_ss·σ_n) on the
   //     velocity-WEAKENING DOFs (the only ones that can spontaneously
   //     accelerate).  The fault is seeded at V_init by SeedEquilibriumPsi_RS,
   //     so ratio ~ 1 means "near steady state" (slow creep); ratio > 1 means
   //     the VW region is stressed above its steady-state strength and would
   //     accelerate even without nucleation — undesirable OUTSIDE the patch.
   // -----------------------------------------------------------------
   real_t outside_max_local    = 0.0;
   int    outside_argmax_local = -1;
   real_t r_threshold = 0.0;
   if (nuc.enabled)
   {
      r_threshold = cfg.outside_safety_factor *
                    std::max(nuc.gradual_overstress.radius_dip_m,
                             nuc.gradual_overstress.radius_strike_m);
   }
   const real_t cx = nuc.gradual_overstress.center_x_m;
   const real_t cy = nuc.gradual_overstress.center_y_m;
   const real_t cz = nuc.gradual_overstress.center_z_m;
   for (int i = 0; i < N; ++i)
   {
      if (rs.b(i) - rs.a(i) <= 0.0) { continue; }   // skip velocity-strengthening
      const real_t dx = dof_coords_3d(3*i + 0) - cx;
      const real_t dy = dof_coords_3d(3*i + 1) - cy;
      const real_t dz = dof_coords_3d(3*i + 2) - cz;
      const real_t r = std::sqrt(dx*dx + dy*dy + dz*dz);
      const bool outside = !nuc.enabled || (r > r_threshold);
      if (!outside) { continue; }
      const real_t tau1 = tau_pre_per_dof(2*i + 0);
      const real_t tau2 = tau_pre_per_dof(2*i + 1);
      const real_t tau_mag = std::sqrt(tau1*tau1 + tau2*tau2);
      const real_t fss = SteadyStateFrictionRS(rs.a(i), rs.b(i), rs.f_0(i),
                                               rs.V_0(i), rs.V_init(i));
      const real_t strength = fss * sigma_n_eff_per_dof(i);
      if (strength <= 0.0) { continue; }
      const real_t ratio = tau_mag / strength;
      if (ratio > outside_max_local)
      {
         outside_max_local    = ratio;
         outside_argmax_local = i;
      }
   }

#ifdef MFEM_USE_MPI
   struct { double value; int rank_; } in{outside_max_local, rank}, out_maxloc;
   MPI_Allreduce(&in, &out_maxloc, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
   const real_t outside_max = static_cast<real_t>(out_maxloc.value);
   const int    owner       = out_maxloc.rank_;
   real_t outside_coord[3]  = { 0.0, 0.0, 0.0 };
   if (owner == rank && outside_argmax_local >= 0)
   {
      outside_coord[0] = dof_coords_3d(3 * outside_argmax_local + 0);
      outside_coord[1] = dof_coords_3d(3 * outside_argmax_local + 1);
      outside_coord[2] = dof_coords_3d(3 * outside_argmax_local + 2);
   }
   MPI_Bcast(outside_coord, 3, MPITypeMap<real_t>::mpi_type, owner, comm);
#else
   const real_t outside_max = outside_max_local;
   const int    owner       = 0;
   real_t outside_coord[3]  = { 0.0, 0.0, 0.0 };
   if (outside_argmax_local >= 0)
   {
      outside_coord[0] = dof_coords_3d(3 * outside_argmax_local + 0);
      outside_coord[1] = dof_coords_3d(3 * outside_argmax_local + 1);
      outside_coord[2] = dof_coords_3d(3 * outside_argmax_local + 2);
   }
#endif

   // -----------------------------------------------------------------
   // 3.  Nucleation well-posedness, over velocity-WEAKENING DOFs INSIDE the
   //     patch (r <= outside_safety_factor · max(radius_*)).  The RS trigger
   //     analog: the overstressed traction must exceed the steady-state
   //     strength so the patch accelerates,
   //       overshoot := |τ_pre + δτ| − f_ss(V_init)·σ_n_eff   (want >= 0).
   //     Also report the overstress "drive" δτ / ((b − a)·σ_n) — the number
   //     of e-folds of slip rate the overstress spans (larger ⇒ stronger
   //     nucleation).  A patch that is entirely velocity-strengthening cannot
   //     nucleate and is gated as a FAIL.
   // -----------------------------------------------------------------
   real_t overshoot_local_best = -std::numeric_limits<real_t>::infinity();
   int    overshoot_local_arg  = -1;
   real_t nuc_peak_local       = 0.0;
   real_t drive_local_max      = 0.0;
   int    vw_in_patch_local    = 0;
   if (nuc.enabled && nuc_params.amplitude_dip.Size() == N)
   {
      for (int i = 0; i < N; ++i)
      {
         const real_t dx = dof_coords_3d(3*i + 0) - cx;
         const real_t dy = dof_coords_3d(3*i + 1) - cy;
         const real_t dz = dof_coords_3d(3*i + 2) - cz;
         const real_t r  = std::sqrt(dx*dx + dy*dy + dz*dz);
         if (r > r_threshold) { continue; }
         const real_t ad = nuc_params.amplitude_dip(i);
         const real_t as = nuc_params.amplitude_strike(i);
         const real_t A  = std::sqrt(ad*ad + as*as);
         nuc_peak_local  = std::max(nuc_peak_local, A);

         const real_t ab = rs.b(i) - rs.a(i);
         if (ab <= 0.0) { continue; }              // VS: cannot nucleate here
         vw_in_patch_local += 1;

         const real_t tau1 = tau_pre_per_dof(2*i + 0);
         const real_t tau2 = tau_pre_per_dof(2*i + 1);
         const real_t fss  = SteadyStateFrictionRS(rs.a(i), rs.b(i), rs.f_0(i),
                                                   rs.V_0(i), rs.V_init(i));
         const real_t strength = fss * sigma_n_eff_per_dof(i);
         const real_t tau_nuc1 = tau1 + ad;
         const real_t tau_nuc2 = tau2 + as;
         const real_t tau_nuc_mag = std::sqrt(tau_nuc1*tau_nuc1 +
                                              tau_nuc2*tau_nuc2);
         const real_t overshoot = tau_nuc_mag - strength;
         if (overshoot > overshoot_local_best)
         {
            overshoot_local_best = overshoot;
            overshoot_local_arg  = i;
         }
         const real_t drive_denom = ab * sigma_n_eff_per_dof(i);
         if (drive_denom > 0.0)
         {
            drive_local_max = std::max(drive_local_max, A / drive_denom);
         }
      }
   }
#ifdef MFEM_USE_MPI
   const real_t nuc_peak = nuc.enabled ? AllreduceMax(nuc_peak_local, comm)
                                       : 0.0;
   struct { double value; int rank_; } pin{overshoot_local_best, rank}, pout;
   MPI_Allreduce(&pin, &pout, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
   const real_t overshoot_best  = static_cast<real_t>(pout.value);
   const int    overshoot_owner = pout.rank_;
   const real_t drive_max =
      nuc.enabled ? AllreduceMax(drive_local_max, comm) : 0.0;
   const int    vw_in_patch =
      nuc.enabled ? AllreduceSumInt(vw_in_patch_local, comm) : 0;
#else
   const real_t nuc_peak       = nuc_peak_local;
   const real_t overshoot_best = overshoot_local_best;
   const int    overshoot_owner = 0;
   const real_t drive_max      = drive_local_max;
   const int    vw_in_patch    = vw_in_patch_local;
#endif

   // Coords + amplitude at the (VW, in-patch) overshoot-argmax DOF, plus the
   // "any DOF inside the patch on any rank" detector (geometry, friction-
   // agnostic — counts ALL in-patch DOFs, not only VW).
   real_t A_at_argmax        = 0.0;
   real_t coord_at_argmax[3] = { 0.0, 0.0, 0.0 };
   int    local_in_patch_count = 0;
   if (nuc.enabled)
   {
      for (int i = 0; i < N; ++i)
      {
         const real_t dx = dof_coords_3d(3*i + 0) - cx;
         const real_t dy = dof_coords_3d(3*i + 1) - cy;
         const real_t dz = dof_coords_3d(3*i + 2) - cz;
         if (std::sqrt(dx*dx + dy*dy + dz*dz) <= r_threshold)
         {
            local_in_patch_count += 1;
         }
      }
   }
   if (nuc.enabled && overshoot_owner == rank && overshoot_local_arg >= 0)
   {
      const int i = overshoot_local_arg;
      const real_t ad = nuc_params.amplitude_dip(i);
      const real_t as = nuc_params.amplitude_strike(i);
      A_at_argmax = std::sqrt(ad*ad + as*as);
      coord_at_argmax[0] = dof_coords_3d(3*i + 0);
      coord_at_argmax[1] = dof_coords_3d(3*i + 1);
      coord_at_argmax[2] = dof_coords_3d(3*i + 2);
   }
#ifdef MFEM_USE_MPI
   int patch_count_global = local_in_patch_count;
   if (nuc.enabled)
   {
      MPI_Allreduce(&local_in_patch_count, &patch_count_global, 1, MPI_INT,
                    MPI_SUM, comm);
      MPI_Bcast(&A_at_argmax,    1, MPITypeMap<real_t>::mpi_type,
                overshoot_owner, comm);
      MPI_Bcast(coord_at_argmax, 3, MPITypeMap<real_t>::mpi_type,
                overshoot_owner, comm);
   }
#else
   const int patch_count_global = local_in_patch_count;
#endif
   const bool nuc_patch_empty    = nuc.enabled && (patch_count_global == 0);
   const bool nuc_patch_all_vs   = nuc.enabled && !nuc_patch_empty &&
                                   (vw_in_patch == 0);

   // -----------------------------------------------------------------
   // 4.  σ_1 azimuth / plunge / magnitude (constant_tensor only).
   // -----------------------------------------------------------------
   real_t s1_azimuth = 0.0, s1_plunge = 0.0, s1_mag = 0.0;
   const bool stress_const =
      (stress.kind == StressSourceKind::ConstantTensor);
   if (stress_const)
   {
      ComputeSigma1AzimuthPlungeMagnitude(stress, s1_azimuth, s1_plunge, s1_mag);
   }

   // -----------------------------------------------------------------
   // 5.  Rank-0 print.
   // -----------------------------------------------------------------
   const int nsteps = (dt_cfl > 0.0)
                       ? static_cast<int>(std::ceil(tfinal / dt_cfl)) : 0;
   const real_t Lnuc_h_min = (h_min_global > 0.0) ? (Lnuc_min / h_min_global) : 0.0;
   const real_t Lnuc_h_max = (h_min_global > 0.0) ? (Lnuc_max / h_min_global) : 0.0;

   if (rank == 0)
   {
      out << "[derived] friction law = rate-and-state (aging)\n";
      out << "[derived] num_fault_global = " << num_fault_global << "\n";
      out << "[derived] cp = " << cp << " m/s, cs = " << cs << " m/s\n";
      out << "[derived] dt_cfl = " << dt_cfl << " s, nsteps = "
          << nsteps << "\n";
      out << "[derived] mesh h_min = " << h_min_global << " m (scalar cp)\n";
      out << "[derived] velocity-weakening (b>a) DOFs = " << vw_cnt
          << " of " << num_fault_global << "\n";
      out << "[derived] steady-state friction f_ss(V_init) min/max = "
          << fss_min << " / " << fss_max << "\n";
      if (all_vs)
      {
         out << "[derived] all fault DOFs are velocity-strengthening "
                "(b<=a); nothing will nucleate\n";
      }
      else
      {
         out << "[derived] L_nuc(RS)=mu*Dc/((b-a)*sigma_n) min/max/mean = "
             << Lnuc_min << " / " << Lnuc_max << " / " << Lnuc_mean << " m\n";
         out << "[derived] L_nuc / h_min min/max = "
             << Lnuc_h_min << " / " << Lnuc_h_max
             << " (NB: should be >= 10 in nucleation patch)\n";
      }
      if (nuc.enabled)
      {
         out << "[derived] |tau_pre| / (f_ss*sigma_n) max OUTSIDE "
                "asperity (VW DOFs) = " << outside_max << "\n";
      }
      else
      {
         out << "[derived] note: no [nucleation] block; an RS rupture must "
                "grow from tau_pre exceeding the steady-state strength on a "
                "velocity-weakening patch.\n";
         out << "[derived] |tau_pre| / (f_ss*sigma_n) max anywhere (VW DOFs) = "
             << outside_max << "\n";
      }
      out << "[derived] nucleation: "
          << (nuc.enabled ? "enabled" : "disabled") << "\n";
      if (nuc.enabled)
      {
         out << "[derived] nucleation peak |F(r) * delta_tau| = "
             << nuc_peak << " Pa\n";
         out << "[derived] velocity-weakening DOFs inside patch = "
             << vw_in_patch << "\n";
         if (nuc_patch_empty)
         {
            out << "[derived] WARNING: zero fault DOFs are inside the "
                   "nucleation patch (r < 3 · max(radius_*)).  Check that "
                   "the nucleation centre lies on the fault and that the "
                   "mesh resolves it.\n";
         }
         else if (nuc_patch_all_vs)
         {
            out << "[derived] WARNING: the nucleation patch is entirely "
                   "velocity-strengthening (b<=a); it cannot host an "
                   "instability.\n";
         }
         else
         {
            out << "[derived] most-overstressed VW DOF at (x="
                << coord_at_argmax[0] << ", y=" << coord_at_argmax[1]
                << ", z=" << coord_at_argmax[2] << ")\n";
            out << "[derived] amplitude |delta_tau| at that DOF = "
                << A_at_argmax << " Pa\n";
            out << "[derived] nucleation overshoot "
                   "(|tau_pre + delta_tau| - f_ss*sigma_n) = "
                << (overshoot_best / 1.0e6) << " MPa "
                << ((overshoot_best >= 0.0) ? "(drives acceleration)"
                                            : "(below steady-state strength)")
                << "\n";
            out << "[derived] overstress drive |delta_tau|/((b-a)*sigma_n) "
                   "max = " << drive_max
                << " e-folds of slip rate (larger => stronger nucleation)\n";
         }
      }
      if (stress_const)
      {
         out << "[derived] sigma_1 azimuth = " << s1_azimuth << " deg, "
                "plunge = " << s1_plunge << " deg, magnitude = "
             << s1_mag << " Pa\n";
      }
      else
      {
         out << "[derived] sigma_1 azimuth from sidecar: TBD "
                "(depth-dependent; see stress sidecar README) — R-002\n";
      }
      out << "[derived] fault DOFs with zero-normal fallback: "
          << num_zero_normal_fallbacks << "\n";
   }

   // -----------------------------------------------------------------
   // 6.  Gating.  Hard FAILs are the physically unambiguous RS cases:
   //       (a) every DOF velocity-strengthening (b<=a) — nothing nucleates,
   //       (b) nucleation enabled but ZERO fault DOFs in the patch,
   //       (c) nucleation enabled but the patch is entirely
   //           velocity-strengthening.
   //     Resolution (L_nuc/h_min), the overstress trigger, and the
   //     outside-asperity ratio are WARN-only (RS nucleation is more
   //     forgiving than the LSW static-yield gate, and these thresholds are
   //     heuristics rather than hard equilibrium violations).
   // -----------------------------------------------------------------
   auto fail = [&](const std::string& msg)
   {
      if (warn_only)
      {
         if (rank == 0)
         {
            const char* reason = skip_gate
               ? "SEAS_SKIP_EQUILIBRIUM_GATE=1"
               : "cfg.abort_on_failure=false";
            out << "[derived] WARNING (gate downgraded via " << reason
                << "): " << msg << "\n";
         }
      }
      else
      {
         MFEM_ABORT(msg);
      }
   };

   if (all_vs)
   {
      fail("[derived] all fault DOFs are velocity-strengthening (b<=a); "
           "nothing will nucleate.");
   }

   if (nuc.enabled && !all_vs)
   {
      if (nuc_patch_empty)
      {
         fail("[derived] FAIL: nucleation enabled but ZERO fault DOFs are "
              "inside the patch (3·max(radius_*)).  Check that the "
              "nucleation centre is on the fault.");
      }
      else if (nuc_patch_all_vs)
      {
         fail("[derived] FAIL: the nucleation patch is entirely "
              "velocity-strengthening (b<=a) — it cannot host a "
              "rate-and-state instability.  Lower a or raise b in the patch.");
      }
      else
      {
         // WARN-only: overstress that does not exceed steady-state strength.
         if (overshoot_best < 0.0 && rank == 0)
         {
            out << "[derived] WARNING: nucleation overshoot = "
                << (overshoot_best / 1.0e6) << " MPa < 0 — |tau_pre+delta_tau| "
                   "stays below f_ss*sigma_n at the peak VW DOF; the patch "
                   "only creeps faster and may not nucleate within tfinal.  "
                   "Increase delta_tau_*_pa.\n";
         }
      }
   }

   // WARN-only: outside-asperity VW region already above steady-state strength.
   if (outside_max >= 1.0 && rank == 0)
   {
      out << "[derived] WARNING: max OUTSIDE |tau_pre|/(f_ss*sigma_n) = "
          << outside_max << " >= 1 on a velocity-weakening DOF (x="
          << outside_coord[0] << ", y=" << outside_coord[1] << ", z="
          << outside_coord[2] << "; owner rank " << owner << ") — that region "
             "is stressed above its steady-state strength and may accelerate "
             "outside the nucleation zone.\n";
   }

   // WARN-only: under-resolved RS process zone.
   if (!all_vs && Lnuc_h_min > 0.0 && Lnuc_h_min < 10.0 && rank == 0)
   {
      out << "[derived] WARNING: L_nuc / h_min min = " << Lnuc_h_min
          << " < 10 — the rate-and-state process zone is under-resolved; "
             "refine the mesh or increase Dc.\n";
   }

   // R-021 WARN-only: sub-critical patch.  Aging-law nucleation feasibility
   // needs the patch radius >~ L_nuc; if the smallest in-patch L_nuc already
   // exceeds the patch radius, the patch cannot grow a self-sustaining rupture
   // (mirrors the L_nuc/h_min resolution WARN, the opposite — dominant — case).
   if (nuc.enabled && !all_vs && !nuc_patch_empty && !nuc_patch_all_vs
       && rank == 0)
   {
      const real_t patch_radius =
         std::max(nuc.gradual_overstress.radius_dip_m,
                  nuc.gradual_overstress.radius_strike_m);
      if (Lnuc_min > patch_radius)
      {
         out << "[derived] WARNING: in-patch L_nuc (min " << Lnuc_min
             << " m) exceeds the nucleation patch radius (" << patch_radius
             << " m) — the patch is sub-critical and may not nucleate a "
                "self-sustaining rupture; lower Dc or enlarge the patch.\n";
      }
   }

   if (num_zero_normal_fallbacks > 0 && rank == 0)
   {
      out << "[derived] WARNING: " << num_zero_normal_fallbacks
          << " fault DOFs hit the zero-normal fallback; "
             "gradual_overstress F(r) at those DOFs may be miscomputed.\n";
   }

   // PASS iff there is a velocity-weakening region AND (no nucleation, OR the
   // patch is non-empty and not entirely velocity-strengthening).
   const bool rs_well_posed = !all_vs
                              && (!nuc.enabled
                                  || (!nuc_patch_empty && !nuc_patch_all_vs));
   if (rank == 0 && !warn_only && rs_well_posed)
   {
      out << "[derived] PASS: rate-and-state initial conditions are "
             "well-posed.\n";
   }

   return outside_max;
}

}  // namespace spatial
}  // namespace seas
}  // namespace mfem
