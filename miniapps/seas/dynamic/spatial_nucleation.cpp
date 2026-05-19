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

// =====================================================================
// SquareOverstress resolver + applicator
// =====================================================================

SquareOverstressPerDOFParams ResolveSquareOverstress(
   const SquareOverstressSpec& spec,
   bool                        enabled,
   const Vector&               dof_coords_3d)
{
   SquareOverstressPerDOFParams p;
   if (!enabled) { return p; }

   const int N = dof_coords_3d.Size() / 3;
   MFEM_VERIFY(dof_coords_3d.Size() == 3 * N,
               "ResolveSquareOverstress: dof_coords_3d.Size() must be 3 * N; "
               "got " << dof_coords_3d.Size());
   MFEM_VERIFY(spec.T_nuc_s > 0.0,
               "ResolveSquareOverstress: spec.T_nuc_s must be > 0");
   MFEM_VERIFY(!spec.patches.empty(),
               "ResolveSquareOverstress: enabled but patches list is empty");

   p.amplitude_dip.SetSize(N);    p.amplitude_dip    = 0.0;
   p.amplitude_strike.SetSize(N); p.amplitude_strike = 0.0;
   p.radial.SetSize(N);           p.radial           = 0.0;

   for (int i = 0; i < N; ++i)
   {
      const real_t x = dof_coords_3d(3 * i + 0);
      const real_t y = dof_coords_3d(3 * i + 1);
      const real_t z = dof_coords_3d(3 * i + 2);
      // Last-match-wins: iterate forward, overwrite when the DOF is in
      // the patch.  This mirrors the StressPatch override convention
      // documented in spatial_friction.hpp.
      for (const SquareOverstressPatch& pp : spec.patches)
      {
         if (pp.contains(x, y, z))
         {
            p.amplitude_dip(i)    = pp.delta_tau_dip_pa;
            p.amplitude_strike(i) = pp.delta_tau_strike_pa;
            p.radial(i)           = 1.0;
         }
      }
   }
   return p;
}

void ApplySquareOverstressIncrement(
   std::vector<DOFData>&                 dof_data,
   const SquareOverstressPerDOFParams&   params,
   real_t                                T_nuc_s,
   real_t                                t_substep_end,
   real_t                                dt_substep)
{
   if (params.amplitude_dip.Size() == 0) { return; }
   if (t_substep_end >= T_nuc_s + dt_substep) { return; }

   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(params.amplitude_dip.Size()    == n
               && params.amplitude_strike.Size() == n,
               "ApplySquareOverstressIncrement: params.amplitude_* size ("
               << params.amplitude_dip.Size()
               << ") != dof_data.size() (" << n << ")");

   const real_t dS = SmoothStepIncrement(t_substep_end, dt_substep, T_nuc_s);
   if (dS <= 0.0) { return; }

   for (int i = 0; i < n; ++i)
   {
      dof_data[i].tau1_nuc += dS * params.amplitude_dip(i);
      dof_data[i].tau2_nuc += dS * params.amplitude_strike(i);
   }
}

// =====================================================================
// InstantaneousOverstressCircular: resolver + one-shot applicator.
// =====================================================================

namespace
{

// Spatial factor f(r) ∈ [0, 1] for the circular cosine-tapered patch.
real_t circular_cosine_factor(real_t r, real_t ri, real_t ro)
{
   if (r <= ri) { return 1.0; }
   if (r >= ro) { return 0.0; }
   // SCEC TPV31 §7: 0.5 · (1 + cos(π(r-ri)/(ro-ri))) interpolates
   // smoothly from 1 at r=ri to 0 at r=ro.  At r=ri: 0.5·2 = 1.  At
   // r=ro: 0.5·(1 + cos π) = 0.5·0 = 0.
   const real_t arg = M_PI * (r - ri) / (ro - ri);
   return 0.5 * (1.0 + std::cos(arg));
}

}  // namespace

InstantaneousOverstressCircularPerDOFParams
ResolveInstantaneousOverstressCircular(
   const InstantaneousOverstressCircularSpec& spec,
   bool                                       enabled,
   const Vector&                              dof_coords_3d,
   const Vector&                              mu_per_dof)
{
   InstantaneousOverstressCircularPerDOFParams p;
   if (!enabled) { return p; }

   const int N = dof_coords_3d.Size() / 3;
   MFEM_VERIFY(dof_coords_3d.Size() == 3 * N,
               "ResolveInstantaneousOverstressCircular: dof_coords_3d.Size() "
               "must be 3 * N; got " << dof_coords_3d.Size());
   MFEM_VERIFY(mu_per_dof.Size() == N,
               "ResolveInstantaneousOverstressCircular: mu_per_dof.Size() ("
               << mu_per_dof.Size() << ") must equal N (" << N << ")");
   MFEM_VERIFY(spec.radius_inner_m > 0.0,
               "ResolveInstantaneousOverstressCircular: radius_inner_m must "
               "be > 0; got " << spec.radius_inner_m);
   MFEM_VERIFY(spec.radius_outer_m >= spec.radius_inner_m,
               "ResolveInstantaneousOverstressCircular: radius_outer_m ("
               << spec.radius_outer_m << ") must be >= radius_inner_m ("
               << spec.radius_inner_m << ")");
   MFEM_VERIFY(spec.mu_ref_pa > 0.0,
               "ResolveInstantaneousOverstressCircular: mu_ref_pa must be > 0; "
               "got " << spec.mu_ref_pa);

   p.amplitude_dip.SetSize(N);    p.amplitude_dip    = 0.0;
   p.amplitude_strike.SetSize(N); p.amplitude_strike = 0.0;
   p.radial.SetSize(N);           p.radial           = 0.0;

   for (int i = 0; i < N; ++i)
   {
      const real_t dx = dof_coords_3d(3 * i + 0) - spec.center_x_m;
      const real_t dy = dof_coords_3d(3 * i + 1) - spec.center_y_m;
      const real_t dz = dof_coords_3d(3 * i + 2) - spec.center_z_m;
      const real_t r  = std::sqrt(dx*dx + dy*dy + dz*dz);

      const real_t f       = circular_cosine_factor(r, spec.radius_inner_m,
                                                    spec.radius_outer_m);
      const real_t mu_scale = mu_per_dof(i) / spec.mu_ref_pa;
      const real_t tau_nuke = spec.delta_tau_peak_pa * f * mu_scale;
      p.amplitude_dip(i)    = spec.dip_fraction    * tau_nuke;
      p.amplitude_strike(i) = spec.strike_fraction * tau_nuke;
      p.radial(i)           = f;
   }
   return p;
}

void ApplyInstantaneousOverstressCircular(
   std::vector<DOFData>&                              dof_data,
   const InstantaneousOverstressCircularPerDOFParams& params)
{
   if (params.amplitude_dip.Size() == 0) { return; }
   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(params.amplitude_dip.Size()    == n
               && params.amplitude_strike.Size() == n,
               "ApplyInstantaneousOverstressCircular: params.amplitude_* "
               "size (" << params.amplitude_dip.Size()
               << ") != dof_data.size() (" << n << ")");
   for (int i = 0; i < n; ++i)
   {
      dof_data[i].tau1_nuc += params.amplitude_dip(i);
      dof_data[i].tau2_nuc += params.amplitude_strike(i);
   }
}

// =====================================================================
// GradualOverstressCompactCircular: resolver + applicator.
// =====================================================================

GradualOverstressCompactCircularPerDOFParams
ResolveGradualOverstressCompactCircular(
   const GradualOverstressCompactCircularSpec& spec,
   bool                                        enabled,
   const Vector&                               dof_coords_3d,
   const DenseMatrix&                          dof_basis)
{
   GradualOverstressCompactCircularPerDOFParams p;
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
               "ResolveGradualOverstressCompactCircular: spec.radius_m must "
               "be > 0; got " << spec.radius_m);
   MFEM_VERIFY(spec.T_nuc_s > 0.0,
               "ResolveGradualOverstressCompactCircular: spec.T_nuc_s must "
               "be > 0; got " << spec.T_nuc_s);

   p.amplitude_dip.SetSize(N);    p.amplitude_dip    = 0.0;
   p.amplitude_strike.SetSize(N); p.amplitude_strike = 0.0;
   p.radial.SetSize(N);           p.radial           = 0.0;

   const real_t R  = spec.radius_m;
   const real_t R2 = R * R;
   const real_t center[3] = { spec.center_x_m,
                              spec.center_y_m,
                              spec.center_z_m };

   for (int i = 0; i < N; ++i)
   {
      const real_t* dof_xyz       = &dof_coords_3d(3 * i);
      const real_t* basis_dip     = &dof_basis(3, i);
      const real_t* basis_strike  = &dof_basis(6, i);

      const real_t v[3] = { dof_xyz[0] - center[0],
                            dof_xyz[1] - center[1],
                            dof_xyz[2] - center[2] };
      const real_t dx = v[0] * basis_dip[0]
                      + v[1] * basis_dip[1]
                      + v[2] * basis_dip[2];
      const real_t ds = v[0] * basis_strike[0]
                      + v[1] * basis_strike[1]
                      + v[2] * basis_strike[2];
      const real_t r2 = dx * dx + ds * ds;

      // SCEC TPV101/102/104 spec Eq. (13): compact-support bell.
      // F(r) = 0 for r ≥ R; F(0) = 1; C∞ at r = R.
      real_t F = 0.0;
      if (r2 < R2)
      {
         F = std::exp(r2 / (r2 - R2));
      }
      p.amplitude_dip(i)    = F * spec.delta_tau_dip_pa;
      p.amplitude_strike(i) = F * spec.delta_tau_strike_pa;
      p.radial(i)           = F;
   }
   return p;
}

void ApplyGradualOverstressCompactCircularIncrement(
   std::vector<DOFData>&                                 dof_data,
   const GradualOverstressCompactCircularPerDOFParams&   params,
   real_t                                                T_nuc_s,
   real_t                                                t_substep_end,
   real_t                                                dt_substep)
{
   if (params.amplitude_dip.Size() == 0) { return; }
   if (t_substep_end >= T_nuc_s + dt_substep) { return; }

   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(params.amplitude_dip.Size()    == n
               && params.amplitude_strike.Size() == n,
               "ApplyGradualOverstressCompactCircularIncrement: "
               "params.amplitude_* size (" << params.amplitude_dip.Size()
               << ") != dof_data.size() (" << n << ")");

   const real_t dS = SmoothStepIncrement(t_substep_end, dt_substep, T_nuc_s);
   if (dS <= 0.0) { return; }

   for (int i = 0; i < n; ++i)
   {
      dof_data[i].tau1_nuc += dS * params.amplitude_dip(i);
      dof_data[i].tau2_nuc += dS * params.amplitude_strike(i);
   }
}

}  // namespace spatial
}  // namespace seas
}  // namespace mfem
