// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// initial_stress.hpp — the effective initial stress sigma_0'(x) of the
// Nitsche [initial_stress] table (rheology plan rev 1, sections 3.3, 6.3).
// One function feeds both the rheology nodes and the fault prestress
// (NitscheFaultSetup::prestress), so the fault and the bulk cannot disagree
// about the initial state.
//
// kind = "lithostatic_ratio" (TPV26/27 description v13, p. 6-7, written for
// a general up vector e_u and a general horizontal ratio tensor B):
//
//    d(x)     = z_s - x . e_u                    depth below the surface
//    sigma_v  = -rho g d,   P_f = rho_f g d
//    Omega(d) = 1 above d_top, (d_bot - d)/(d_bot - d_top) between, 0 below
//    H        = I - e_u e_u
//    sigma_0  = sigma_v e_u e_u + Omega ((sigma_v + P_f) B - P_f H)
//               + (1 - Omega) sigma_v H
//    sigma_0' = sigma_0 + P_f I
//
// B is symmetric with B e_u = 0.  For TPV26/27 in the code frame (X strike,
// Y fault normal, Z up): B_XX = b11, B_YY = b33, B_XY = b13.
//
// kind = "constant_tensor": sigma_0' is a given effective tensor.
//
// Stresses are tension positive, row-major 3x3, global code frame.
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_INITIAL_STRESS_HPP
#define MFEM_SEAS_INITIAL_STRESS_HPP

#include "mfem.hpp"
#include "rheology_model.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

namespace mfem
{
namespace seas
{

struct InitialStressSpec
{
   enum class Kind { ConstantTensor, LithostaticRatio };
   Kind kind = Kind::ConstantTensor;

   /// ConstantTensor: sigma_0' (Pa), row-major.
   real_t sigma[9] = {};

   /// LithostaticRatio.
   real_t rho_rock = 0.0;            ///< kg/m^3
   real_t rho_fluid = 0.0;           ///< kg/m^3
   real_t gravity = 0.0;             ///< m/s^2
   real_t surface_elevation = 0.0;   ///< z_s along e_u (m)
   real_t B[9] = {};                 ///< horizontal ratio tensor, row-major
   real_t taper_top_depth = 0.0;     ///< m
   real_t taper_bottom_depth = 0.0;  ///< m
   real_t up[3] = {0.0, 0.0, 1.0};   ///< e_u ([fault_geometry] up)
};

class InitialStress
{
public:
   /// Points above the surface by more than this (m) abort; closer ones
   /// (mesh round-off on the free surface) are evaluated at depth 0.
   static constexpr real_t kAboveSurfaceTol = 1e-6;
   /// Asymmetry and |B up| allowed, relative to the largest entry (a tensor
   /// typed in a deck for a tilted up vector is horizontal only to round-off).
   static constexpr real_t kTensorTol = 1e-12;

   /// max |m_ij - m_ji|
   static real_t Asymmetry(const real_t m[9])
   {
      return std::max({std::abs(m[1] - m[3]), std::abs(m[2] - m[6]), std::abs(m[5] - m[7])});
   }
   static real_t MaxEntry(const real_t m[9])
   {
      real_t a = 0.0;
      for (int i = 0; i < 9; i++) { a = std::max(a, std::abs(m[i])); }
      return a;
   }

   /// "" when the spec is usable, else the first problem (deck key names).
   static std::string Validate(const InitialStressSpec &s)
   {
      std::ostringstream e;
      if (s.kind == InitialStressSpec::Kind::ConstantTensor)
      {
         if (Asymmetry(s.sigma) > kTensorTol * MaxEntry(s.sigma))
         {
            e << "[initial_stress] constant_tensor must be symmetric";
         }
         return e.str();
      }
      const real_t *u = s.up;
      const real_t nu = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
      const real_t *B = s.B;
      const real_t bmax = MaxEntry(B);
      if (!(s.rho_rock > 0.0)) { e << "[initial_stress] rho_rock_kg_m3 must be > 0"; }
      else if (!(s.rho_fluid >= 0.0)) { e << "[initial_stress] rho_fluid_kg_m3 must be >= 0"; }
      else if (!(s.gravity > 0.0)) { e << "[initial_stress] gravity_m_s2 must be > 0"; }
      else if (!(s.taper_top_depth >= 0.0 && s.taper_bottom_depth > s.taper_top_depth))
      {
         e << "[initial_stress] need 0 <= taper_top_depth_m < taper_bottom_depth_m, got "
           << s.taper_top_depth << " and " << s.taper_bottom_depth;
      }
      else if (!(nu > 0.0)) { e << "[initial_stress] the up vector is zero"; }
      else if (Asymmetry(B) > kTensorTol * bmax)
      {
         e << "[initial_stress] ratio tensor B must be symmetric";
      }
      else
      {
         for (int i = 0; i < 3; i++)
         {
            const real_t Be = (B[3 * i] * u[0] + B[3 * i + 1] * u[1] + B[3 * i + 2] * u[2]) / nu;
            if (std::abs(Be) > kTensorTol * bmax)
            {
               e << "[initial_stress] ratio tensor B must be horizontal (B up = 0): row "
                 << i << " gives " << Be;
               break;
            }
         }
      }
      return e.str();
   }

   explicit InitialStress(const InitialStressSpec &spec) : s_(spec)
   {
      const std::string err = Validate(spec);
      MFEM_VERIFY(err.empty(), err);
      Symmetrize(s_.sigma);
      Symmetrize(s_.B);
      if (s_.kind == InitialStressSpec::Kind::ConstantTensor) { return; }
      const real_t *u = s_.up;
      const real_t nu = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
      for (int i = 0; i < 3; i++) { eu_[i] = u[i] / nu; }
   }

   const InitialStressSpec &Spec() const { return s_; }

   /// Depth below the surface (m), clamped at 0 within kAboveSurfaceTol.
   real_t Depth(const real_t x[3]) const
   {
      const real_t d = s_.surface_elevation - (x[0] * eu_[0] + x[1] * eu_[1] + x[2] * eu_[2]);
      MFEM_VERIFY(d >= -kAboveSurfaceTol, "[initial_stress] point (" << x[0] << ", " << x[1]
                  << ", " << x[2] << ") is " << -d << " m above surface_elevation_m");
      return std::max(real_t(0.0), d);
   }

   /// Omega(d)
   real_t Taper(real_t d) const
   {
      if (d <= s_.taper_top_depth) { return 1.0; }
      if (d >= s_.taper_bottom_depth) { return 0.0; }
      return (s_.taper_bottom_depth - d) / (s_.taper_bottom_depth - s_.taper_top_depth);
   }

   /// sigma_0'(x), row-major 3x3.
   void EvalEffective(const real_t x[3], real_t sig[9]) const
   {
      if (s_.kind == InitialStressSpec::Kind::ConstantTensor)
      {
         for (int i = 0; i < 9; i++) { sig[i] = s_.sigma[i]; }
         return;
      }
      const real_t d = Depth(x);
      const real_t sv = -s_.rho_rock * s_.gravity * d;
      const real_t pf = s_.rho_fluid * s_.gravity * d;
      const real_t om = Taper(d);
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            const real_t uu = eu_[i] * eu_[j];
            const real_t H = (i == j ? 1.0 : 0.0) - uu;
            sig[3 * i + j] = sv * uu + om * ((sv + pf) * s_.B[3 * i + j] - pf * H)
                             + (1.0 - om) * sv * H + (i == j ? pf : 0.0);
         }
      }
   }

   /// sigma_0'(x) in the rheology layout [xx, yy, zz, xy, yz, xz].
   void EvalEffectiveSym6(const real_t x[3], real_t sig[6]) const
   {
      real_t m[9];
      EvalEffective(x, m);
      sym6::FromRowMajor(m, sig);
   }

private:
   static void Symmetrize(real_t m[9])
   {
      const int ij[3][2] = {{1, 3}, {2, 6}, {5, 7}};
      for (const auto &p : ij) { m[p[0]] = m[p[1]] = 0.5 * (m[p[0]] + m[p[1]]); }
   }

   InitialStressSpec s_;
   real_t eu_[3] = {0.0, 0.0, 1.0};
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_INITIAL_STRESS_HPP
