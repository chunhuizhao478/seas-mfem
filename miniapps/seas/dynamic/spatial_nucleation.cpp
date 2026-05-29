// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/spatial_nucleation.cpp — Phase N implementation.

#include "spatial_nucleation.hpp"

#include <cmath>

namespace mfem
{
namespace seas
{
namespace spatial
{

// =====================================================================
// Temporal helpers
// =====================================================================

real_t SmoothStep(real_t current_time, real_t t0)
{
   if (current_time <= 0.0) { return 0.0; }
   if (current_time < t0)
   {
      const real_t tau = current_time - t0;
      return std::exp(tau * tau / (current_time * (current_time - 2.0 * t0)));
   }
   return 1.0;
}

real_t SmoothStepIncrement(real_t current_time, real_t dt, real_t t0)
{
   MFEM_ASSERT(std::isfinite(dt) && dt > 0.0,
               "SmoothStepIncrement: dt must be finite and positive; got dt = "
               << dt);
   return SmoothStep(current_time, t0) - SmoothStep(current_time - dt, t0);
}

// =====================================================================
// Spatial helper
// =====================================================================

real_t GaussianFactorFaceLocal(const real_t* dof_xyz,
                               const real_t* dof_basis_dip,
                               const real_t* dof_basis_strike,
                               const real_t* center_xyz,
                               real_t radius_dip_m,
                               real_t radius_strike_m)
{
   MFEM_ASSERT(radius_dip_m    > 0.0,
               "GaussianFactorFaceLocal: radius_dip_m must be > 0; got "
               << radius_dip_m);
   MFEM_ASSERT(radius_strike_m > 0.0,
               "GaussianFactorFaceLocal: radius_strike_m must be > 0; got "
               << radius_strike_m);

   const real_t v[3] = { dof_xyz[0] - center_xyz[0],
                         dof_xyz[1] - center_xyz[1],
                         dof_xyz[2] - center_xyz[2] };

   const real_t dx = v[0] * dof_basis_dip[0]
                   + v[1] * dof_basis_dip[1]
                   + v[2] * dof_basis_dip[2];
   const real_t ds = v[0] * dof_basis_strike[0]
                   + v[1] * dof_basis_strike[1]
                   + v[2] * dof_basis_strike[2];

   const real_t ax = dx / radius_dip_m;
   const real_t as = ds / radius_strike_m;
   const real_t exponent = -(ax * ax + as * as);
   // Avoid denormal-underflow noise on a far-field DOF.
   if (exponent < -700.0) { return 0.0; }
   return std::exp(exponent);
}

// =====================================================================
// Resolver
// =====================================================================

GradualOverstressPerDOFParams ResolveGradualOverstress(
   const GradualOverstressSpec& spec,
   bool                         enabled,
   const Vector&                dof_coords_3d,
   const DenseMatrix&           dof_basis)
{
   GradualOverstressPerDOFParams p;  // zero-sized by default
   if (!enabled) { return p; }

   const int N = dof_coords_3d.Size() / 3;
   MFEM_VERIFY(dof_coords_3d.Size() == 3 * N,
               "ResolveGradualOverstress: dof_coords_3d.Size() must be 3 * N; "
               "got " << dof_coords_3d.Size());
   MFEM_VERIFY(dof_basis.Height() == 9 && dof_basis.Width() == N,
               "ResolveGradualOverstress: dof_basis must have shape (9, N); "
               "got (" << dof_basis.Height() << ", " << dof_basis.Width()
               << ")");
   MFEM_VERIFY(spec.radius_dip_m    > 0.0,
               "ResolveGradualOverstress: spec.radius_dip_m must be > 0");
   MFEM_VERIFY(spec.radius_strike_m > 0.0,
               "ResolveGradualOverstress: spec.radius_strike_m must be > 0");
   MFEM_VERIFY(spec.T_nuc_s         > 0.0,
               "ResolveGradualOverstress: spec.T_nuc_s must be > 0");

   p.amplitude_dip.SetSize(N);
   p.amplitude_strike.SetSize(N);
   p.radial.SetSize(N);

   const real_t center[3] = { spec.center_x_m,
                              spec.center_y_m,
                              spec.center_z_m };

   for (int i = 0; i < N; ++i)
   {
      // `dof_basis` is column-major (9, N); &(dof_basis(r, i)) is a
      // contiguous 3-vector pointer at rows r..r+2 of column i.
      const real_t* dof_xyz       = &dof_coords_3d(3 * i);
      const real_t* basis_dip     = &dof_basis(3, i);
      const real_t* basis_strike  = &dof_basis(6, i);

      const real_t F = GaussianFactorFaceLocal(dof_xyz, basis_dip,
                                               basis_strike, center,
                                               spec.radius_dip_m,
                                               spec.radius_strike_m);
      p.amplitude_dip(i)    = F * spec.delta_tau_dip_pa;
      p.amplitude_strike(i) = F * spec.delta_tau_strike_pa;
      p.radial(i)           = F;
   }
   return p;
}

// =====================================================================
// Per-sub-step accumulator
// =====================================================================

void ApplyGradualOverstressIncrement(
   std::vector<DOFData>&                  dof_data,
   const GradualOverstressPerDOFParams&   params,
   real_t                                 T_nuc_s,
   real_t                                 t_substep_end,
   real_t                                 dt_substep)
{
   // Resolver returned zero-sized arrays — either nucleation is disabled
   // (cfg.nucleation.enabled = false) OR this rank has zero local fault
   // DOFs.  Either way the per-DOF accumulator has nothing to do.
   if (params.amplitude_dip.Size() == 0)
   {
      return;
   }

   // Telescoped to full target: even a sub-ULP increment wouldn't move
   // the accumulator.  Saves the per-DOF loop for every post-nucleation
   // sub-step.
   if (t_substep_end >= T_nuc_s + dt_substep)
   {
      return;
   }

   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(params.amplitude_dip.Size()    == n
               && params.amplitude_strike.Size() == n,
               "ApplyGradualOverstressIncrement: params.amplitude_* size ("
               << params.amplitude_dip.Size()
               << ") != dof_data.size() (" << n << ")");

   const real_t dS = SmoothStepIncrement(t_substep_end, dt_substep, T_nuc_s);
   // R3-004 (mirror of tpv104_nucleation.hpp:124): clamp non-positive
   // increments to 0 — smoothStep is monotone non-decreasing but two
   // close-to-1 values subtracted via different code paths can produce
   // a sub-ULP negative drift; `dS <= 0` guards against it.
   if (dS <= 0.0) { return; }

   for (int i = 0; i < n; ++i)
   {
      dof_data[i].tau1_nuc += dS * params.amplitude_dip(i);
      dof_data[i].tau2_nuc += dS * params.amplitude_strike(i);
      // sigma_n_nuc is intentionally not updated.
   }
}

void ApplyGradualOverstressAbsolute(
   std::vector<DOFData>&                  dof_data,
   const GradualOverstressPerDOFParams&   params,
   real_t                                 T_nuc_s,
   real_t                                 t)
{
   // Disabled (resolver returned zero-sized) or this rank has no fault DOFs.
   if (params.amplitude_dip.Size() == 0) { return; }

   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(params.amplitude_dip.Size()    == n
               && params.amplitude_strike.Size() == n,
               "ApplyGradualOverstressAbsolute: params.amplitude_* size ("
               << params.amplitude_dip.Size()
               << ") != dof_data.size() (" << n << ")");

   // SET (not accumulate) the absolute SCEC ramp value at stage time t.
   // SmoothStep is 0 for t<=0, 1 for t>=T_nuc_s, monotone in between — so this
   // is the no-op-equivalent steady target after the ramp and is idempotent
   // under repeated RK-stage evaluation (no double-apply).
   const real_t S = SmoothStep(t, T_nuc_s);
   for (int i = 0; i < n; ++i)
   {
      dof_data[i].tau1_nuc = S * params.amplitude_dip(i);
      dof_data[i].tau2_nuc = S * params.amplitude_strike(i);
      // sigma_n_nuc is intentionally not updated.
   }
}

// =====================================================================
// Phase 7 — shared in-fault-plane radial distance helper
// =====================================================================

namespace
{
/// In-fault-plane distance from `center` to `dof_xyz`, measured via the
/// per-DOF (dip, strike) basis: r = sqrt(((dof−c)·dip)² + ((dof−c)·strike)²).
/// For the planar y=0 fault (dip=(0,0,∓1), strike=(±1,0,0)) this is exactly
/// the SCEC sqrt(Δalong_strike² + Δdown_dip²).
real_t InFaultPlaneRadius(const real_t* dof_xyz,
                          const real_t* center_xyz,
                          const real_t* dip,
                          const real_t* strike)
{
   const real_t v[3] = { dof_xyz[0] - center_xyz[0],
                         dof_xyz[1] - center_xyz[1],
                         dof_xyz[2] - center_xyz[2] };
   const real_t pd = v[0]*dip[0]    + v[1]*dip[1]    + v[2]*dip[2];
   const real_t ps = v[0]*strike[0] + v[1]*strike[1] + v[2]*strike[2];
   return std::sqrt(pd * pd + ps * ps);
}
}  // namespace

// =====================================================================
// Phase 7 — compact-circular gradual overstress (TPV102/104)
// =====================================================================

real_t CompactBellFactor(real_t r, real_t R)
{
   MFEM_ASSERT(R > 0.0, "CompactBellFactor: R must be > 0; got " << R);
   if (r >= R) { return 0.0; }
   const real_t r2 = r * r;
   const real_t R2 = R * R;
   return std::exp(r2 / (r2 - R2));
}

CompactCircularPerDOFParams ResolveGradualOverstressCompactCircular(
   const GradualOverstressCompactCircularSpec& spec,
   bool                                        enabled,
   const Vector&                               dof_coords_3d,
   const DenseMatrix&                          dof_basis)
{
   CompactCircularPerDOFParams p;   // zero-sized by default
   if (!enabled) { return p; }

   const int N = dof_coords_3d.Size() / 3;
   MFEM_VERIFY(dof_coords_3d.Size() == 3 * N,
               "ResolveGradualOverstressCompactCircular: dof_coords_3d.Size() "
               "must be 3 * N; got " << dof_coords_3d.Size());
   MFEM_VERIFY(dof_basis.Height() == 9 && dof_basis.Width() == N,
               "ResolveGradualOverstressCompactCircular: dof_basis must have "
               "shape (9, N); got (" << dof_basis.Height() << ", "
               << dof_basis.Width() << ")");
   MFEM_VERIFY(spec.radius_m > 0.0,
               "ResolveGradualOverstressCompactCircular: spec.radius_m must be "
               "> 0; got " << spec.radius_m);
   MFEM_VERIFY(spec.T_nuc_s > 0.0,
               "ResolveGradualOverstressCompactCircular: spec.T_nuc_s must be "
               "> 0; got " << spec.T_nuc_s);

   p.amplitude_strike.SetSize(N);
   p.radial.SetSize(N);

   const real_t center[3] = { spec.center_x_m, spec.center_y_m, spec.center_z_m };
   for (int i = 0; i < N; ++i)
   {
      const real_t* dof_xyz      = &dof_coords_3d(3 * i);
      const real_t* basis_dip    = &dof_basis(3, i);
      const real_t* basis_strike = &dof_basis(6, i);
      const real_t  r = InFaultPlaneRadius(dof_xyz, center, basis_dip,
                                           basis_strike);
      const real_t  F = CompactBellFactor(r, spec.radius_m);
      p.amplitude_strike(i) = F * spec.delta_tau_pa;
      p.radial(i)           = F;
   }
   return p;
}

void ApplyGradualOverstressCompactCircularIncrement(
   std::vector<DOFData>&                  dof_data,
   const CompactCircularPerDOFParams&     params,
   real_t                                 T_nuc_s,
   real_t                                 t_substep_end,
   real_t                                 dt_substep)
{
   // Disabled (resolver returned zero-sized) or this rank has no fault DOFs.
   if (params.amplitude_strike.Size() == 0) { return; }

   // Telescoped to the full target — even a sub-ULP increment is a no-op.
   if (t_substep_end >= T_nuc_s + dt_substep) { return; }

   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(params.amplitude_strike.Size() == n,
               "ApplyGradualOverstressCompactCircularIncrement: "
               "params.amplitude_strike size ("
               << params.amplitude_strike.Size()
               << ") != dof_data.size() (" << n << ")");

   const real_t dS = SmoothStepIncrement(t_substep_end, dt_substep, T_nuc_s);
   // Same monotone-non-decreasing round-off guard as the Gaussian path.
   if (dS <= 0.0) { return; }

   for (int i = 0; i < n; ++i)
   {
      dof_data[i].tau2_nuc += dS * params.amplitude_strike(i);
      // tau1_nuc / sigma_n_nuc intentionally not updated (pure strike-slip).
   }
}

void ApplyGradualOverstressCompactCircularAbsolute(
   std::vector<DOFData>&                  dof_data,
   const CompactCircularPerDOFParams&     params,
   real_t                                 T_nuc_s,
   real_t                                 t)
{
   if (params.amplitude_strike.Size() == 0) { return; }

   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(params.amplitude_strike.Size() == n,
               "ApplyGradualOverstressCompactCircularAbsolute: "
               "params.amplitude_strike size ("
               << params.amplitude_strike.Size()
               << ") != dof_data.size() (" << n << ")");

   const real_t S = SmoothStep(t, T_nuc_s);
   for (int i = 0; i < n; ++i)
   {
      dof_data[i].tau2_nuc = S * params.amplitude_strike(i);
      // tau1_nuc / sigma_n_nuc intentionally not updated (pure strike-slip).
   }
}

// =====================================================================
// Phase 7 — instantaneous circular overstress (TPV31, one-shot)
// =====================================================================

real_t CosineTaperFactor(real_t r, real_t R, real_t taper)
{
   MFEM_ASSERT(R > 0.0, "CosineTaperFactor: R must be > 0; got " << R);
   if (r <= R) { return 1.0; }
   if (taper <= 0.0 || r >= R + taper) { return 0.0; }
   // r in (R, R+taper): cosine ramp from 1 (at R) to 0 (at R+taper).
   return 0.5 * (1.0 + std::cos(M_PI * (r - R) / taper));
}

InstantaneousOverstressPerDOFParams ResolveInstantaneousOverstressCircular(
   const InstantaneousOverstressCircularSpec& spec,
   bool                                       enabled,
   const Vector&                              dof_coords_3d,
   const DenseMatrix&                         dof_basis,
   const std::function<real_t(real_t, real_t, real_t)>& mu_at_xyz)
{
   InstantaneousOverstressPerDOFParams p;   // zero-sized by default
   if (!enabled) { return p; }

   const int N = dof_coords_3d.Size() / 3;
   MFEM_VERIFY(dof_coords_3d.Size() == 3 * N,
               "ResolveInstantaneousOverstressCircular: dof_coords_3d.Size() "
               "must be 3 * N; got " << dof_coords_3d.Size());
   MFEM_VERIFY(dof_basis.Height() == 9 && dof_basis.Width() == N,
               "ResolveInstantaneousOverstressCircular: dof_basis must have "
               "shape (9, N); got (" << dof_basis.Height() << ", "
               << dof_basis.Width() << ")");
   MFEM_VERIFY(spec.radius_m > 0.0,
               "ResolveInstantaneousOverstressCircular: spec.radius_m must be "
               "> 0; got " << spec.radius_m);
   MFEM_VERIFY(spec.taper_m >= 0.0,
               "ResolveInstantaneousOverstressCircular: spec.taper_m must be "
               ">= 0; got " << spec.taper_m);

   // Phase 10 (TPV31 spec p. 7): per-DOF mu(depth)/mu_ref amplitude scaling,
   // active only when a reference modulus is set AND a mu lookup is supplied.
   // Otherwise the scale is a uniform 1.0 (byte-identical to the pre-Phase-10
   // path, and to every non-TPV31 config that omits mu_ref_pa).
   const bool mu_scaled =
      (spec.mu_ref_pa > 0.0) && static_cast<bool>(mu_at_xyz);

   p.amplitude_strike.SetSize(N);
   const real_t center[3] = { spec.center_x_m, spec.center_y_m, spec.center_z_m };
   for (int i = 0; i < N; ++i)
   {
      const real_t* dof_xyz      = &dof_coords_3d(3 * i);
      const real_t* basis_dip    = &dof_basis(3, i);
      const real_t* basis_strike = &dof_basis(6, i);
      const real_t  r = InFaultPlaneRadius(dof_xyz, center, basis_dip,
                                           basis_strike);
      real_t mu_scale = 1.0;
      if (mu_scaled)
      {
         const real_t mu_local = mu_at_xyz(dof_xyz[0], dof_xyz[1], dof_xyz[2]);
         MFEM_VERIFY(mu_local > 0.0,
                     "ResolveInstantaneousOverstressCircular: mu_at_xyz("
                     << dof_xyz[0] << ", " << dof_xyz[1] << ", " << dof_xyz[2]
                     << ") returned " << mu_local << " (must be > 0).");
         mu_scale = mu_local / spec.mu_ref_pa;
      }
      p.amplitude_strike(i) =
         CosineTaperFactor(r, spec.radius_m, spec.taper_m)
         * spec.delta_tau_pa * mu_scale;
   }
   return p;
}

}  // namespace spatial
}  // namespace seas
}  // namespace mfem
