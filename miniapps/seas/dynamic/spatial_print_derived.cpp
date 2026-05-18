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
constexpr real_t kNucleationOvershootGap = 0.9;   // 90% of budget

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
   // 3.  Nucleation-budget check (only when nuc enabled).
   // -----------------------------------------------------------------
   real_t overshoot_local_best = -std::numeric_limits<real_t>::infinity();
   int    overshoot_local_arg  = -1;
   real_t nuc_peak_local       = 0.0;
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
         const real_t ad  = nuc_params.amplitude_dip(i);
         const real_t as  = nuc_params.amplitude_strike(i);
         const real_t A   = std::sqrt(ad*ad + as*as);
         const real_t budget = (lsw.mu_s(i) - lsw.mu_d(i))
                              * sigma_n_eff_per_dof(i);
         const real_t overshoot = A - budget;
         if (overshoot > overshoot_local_best)
         {
            overshoot_local_best = overshoot;
            overshoot_local_arg  = i;
         }
         nuc_peak_local = std::max(nuc_peak_local, A);
      }
   }
#ifdef MFEM_USE_MPI
   const real_t nuc_peak = nuc.enabled ? AllreduceMax(nuc_peak_local, comm)
                                       : 0.0;
   struct { double value; int rank_; } pin{overshoot_local_best, rank}, pout;
   MPI_Allreduce(&pin, &pout, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
   const real_t overshoot_best = static_cast<real_t>(pout.value);
   const int    overshoot_owner = pout.rank_;
#else
   const real_t nuc_peak       = nuc_peak_local;
   const real_t overshoot_best = overshoot_local_best;
   const int    overshoot_owner = 0;
#endif

   // Compute the budget at the overshoot-argmax DOF (for the print).
   // Also detect the "no DOF inside nucleation patch on any rank" case
   // so we don't silently PASS a configuration that can't possibly
   // nucleate.
   real_t budget_at_argmax = 0.0;
   real_t A_at_argmax      = 0.0;
   real_t coord_at_argmax[3] = { 0.0, 0.0, 0.0 };
   const int local_in_patch_count = (nuc.enabled && overshoot_local_arg >= 0)
                                     ? 1 : 0;
   if (nuc.enabled && overshoot_owner == rank && overshoot_local_arg >= 0)
   {
      const int i = overshoot_local_arg;
      budget_at_argmax = (lsw.mu_s(i) - lsw.mu_d(i)) * sigma_n_eff_per_dof(i);
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
      MPI_Bcast(&budget_at_argmax, 1, MPITypeMap<real_t>::mpi_type,
                overshoot_owner, comm);
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
            out << "[derived] amplitude at that DOF = " << A_at_argmax << " Pa\n";
            out << "[derived] (mu_s - mu_d)*sigma_n_eff at that DOF = "
                << budget_at_argmax << " Pa\n";
            out << "[derived] nucleation overshoot = "
                << (overshoot_best / 1.0e6) << " MPa "
                << ((overshoot_best >= 0.0) ? "(sufficient)" : "(INSUFFICIENT)")
                << "\n";
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
   //    (b) outside_max >= 1.0 (supercritical),
   //    (c) nuc enabled and overshoot < (1 - 0.9) * budget == budget * (1 - 0.9)
   //        [equivalent to A < 0.9 * budget],
   //    (d) nuc disabled and outside_max < 1.0 (cannot nucleate at all).
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
         // A < 0.9 * budget ⇒ overshoot = A - budget < budget * (-0.1).
         const real_t budget_threshold =
            budget_at_argmax * kNucleationOvershootGap;
         const real_t safety_gap = budget_at_argmax - budget_threshold;
         if (A_at_argmax < budget_threshold)
         {
            std::ostringstream m;
            m << "[derived] FAIL: nucleation overshoot = "
              << (overshoot_best / 1.0e6)
              << " MPa (INSUFFICIENT, gap = " << safety_gap
              << " Pa).  Increase delta_tau_*_pa or reduce mu_s in the "
                 "nucleation patch.";
            fail(m.str());
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

   if (rank == 0 && !warn_only && outside_max < 1.0
       && (!nuc.enabled
           || all_barriers
           || (!nuc_patch_empty
               && A_at_argmax >= budget_at_argmax * kNucleationOvershootGap)))
   {
      out << "[derived] PASS: initial conditions are well-posed.\n";
   }

   return outside_max;
}

}  // namespace spatial
}  // namespace seas
}  // namespace mfem
