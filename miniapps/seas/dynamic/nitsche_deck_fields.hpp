// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_deck_fields.hpp — the deck's material and fault stress fields for
// the Nitsche dynamic-rupture driver (TPV6/TPV7/TPV31 plan rev 1, items A1-A3;
// Chunhui Zhao, 2026-10-19).
//
//   NitscheDeckMaterial / BuildNitscheDeckMaterial
//       the MaterialField of [material] kind constant, depth_profile_1d or
//       halfspace_across_fault, the wrappers that own its coefficients, and
//       mu_at(x), the shear modulus at a point where it is single valued
//       (empty for halfspace_across_fault: mu jumps across the fault).
//   MaterialJumps / CountElementsAcrossJumps
//       the surfaces where the material is discontinuous (the halfspace
//       plane; the depths where the layer table jumps) and the elements whose
//       vertices lie strictly on both sides of one.  The operator samples the
//       material once per element at its centroid (NitscheMaterialSample), so
//       such an element would carry one side's material across the jump; the
//       driver refuses the mesh.
//   ScaleByShearModulus
//       sigma(x) -> (mu(x) / mu_ref) sigma(x)   ([stress_scale_by_mu], TPV31).
//   InstantaneousCircularOverstress
//       the [nucleation] kind "instantaneous_overstress_circular" (TPV31) as a
//       global tensor, present from t = 0.
//   MaxFaultNormalDeviation
//       the largest angle (as its sine) between a fault face's normal and the
//       reference normal.  Both nucleation kinds are one global tensor built
//       in the pinned frame (shat, nhat), which gives the intended traction
//       only on a planar fault; the driver refuses them on any other.

#ifndef MFEM_SEAS_NITSCHE_DECK_FIELDS_HPP
#define MFEM_SEAS_NITSCHE_DECK_FIELDS_HPP

#include "mfem.hpp"
#include "heterogeneous_material.hpp"
#include "spatial_nucleation.hpp"
#include "../fault/fault_orientation.hpp"
#include "../spatial/code/spatial_friction.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

using StressFn = std::function<void(const real_t x[3], real_t sigma[9])>;
using ShearModulusFn = std::function<real_t(const real_t x[3])>;

/// Discontinuities of a deck material: one plane (halfspace_across_fault) or
/// a list of depths along an axis (depth_profile_1d).
struct MaterialJumps
{
   bool has_plane = false;
   real_t x0[3] = {0.0, 0.0, 0.0}, n[3] = {0.0, 0.0, 0.0};   // n unit
   int depth_axis = 2;            // 0, 1, 2 = x, y, z
   std::vector<real_t> depths;    // positive depths of the jumps
};

struct NitscheDeckMaterial
{
   std::unique_ptr<DepthProfile1DMaterial> depth;
   std::unique_ptr<HalfspaceAcrossFaultMaterial> halfspace;
   MaterialField field;
   ShearModulusFn mu_at;   ///< empty when mu is not single valued on the fault
   MaterialJumps jumps;
   std::string log;        ///< what was built, for rank 0's log
};

/// Depth of a point as MakeDepthProfile1DMaterial measures it: max(0, -z)
/// for axis z, the raw coordinate for x and y.
inline real_t DeckDepth(const real_t x[3], int axis)
{
   if (axis == 2) { return std::max(real_t(0.0), -x[2]); }
   return x[axis];
}

/// The point at depth d on the axis (the other coordinates zero).
inline void DeckPointAtDepth(real_t d, int axis, real_t x[3])
{
   x[0] = x[1] = x[2] = 0.0;
   x[axis] = (axis == 2) ? -d : d;
}

/// The material of the deck.  Refuses the sidecar kind (not wired in the
/// Nitsche driver).  Each layer top below the surface is a jump when the
/// material just above it differs from the material just below it, judged
/// with the profile's own evaluator.
inline NitscheDeckMaterial BuildNitscheDeckMaterial(const spatial::SpatialFrictionConfig &cfg)
{
   NitscheDeckMaterial m;
   std::ostringstream log;
   log << std::setprecision(10);
   MFEM_VERIFY(!cfg.velocity.use_sidecar,
               "Nitsche driver: the sidecar material ([velocity] use_sidecar = true) is not wired");
   switch (cfg.material.kind)
   {
      case spatial::MaterialKind::Constant:
      {
         const real_t mu0 = cfg.material_fallback.mu;
         m.field = MaterialField::MakeConstant(cfg.material_fallback.lambda, mu0,
                                               cfg.material_fallback.rho);
         m.mu_at = [mu0](const real_t *) { return mu0; };
         log << "  material: constant (lambda " << cfg.material_fallback.lambda << ", mu " << mu0
             << ", rho " << cfg.material_fallback.rho << ")\n";
         break;
      }
      case spatial::MaterialKind::DepthProfile1D:
      {
         const auto &layers = cfg.material.profile_layers;
         const int axis = cfg.material.depth_axis - 'x';
         m.depth = MakeDepthProfile1DMaterial(layers, cfg.material.depth_axis);
         m.field = m.depth->field;
         const DepthProfile1DMaterial *dm = m.depth.get();
         m.mu_at = [dm](const real_t *x) {
            real_t l, mu, r;
            dm->eval_at_xyz(x[0], x[1], x[2], l, mu, r);
            return mu;
         };
         m.jumps.depth_axis = axis;
         log << "  material: depth_profile_1d, depth axis " << cfg.material.depth_axis << ", "
             << layers.size() << " layers:\n";
         for (size_t i = 0; i < layers.size(); i++)
         {
            const auto &L = layers[i];
            log << "    " << i << ": depth [" << L.depth_top_m << ", " << L.depth_bot_m
                << "] m, vp " << L.vp_ms << ", vs " << L.vs_ms << ", rho " << L.rho_kgm3
                << ", " << L.interp << "\n";
         }
         for (size_t i = 1; i < layers.size(); i++)
         {
            const real_t d = layers[i].depth_top_m;
            if (d <= 0.0) { continue; }
            if (!m.jumps.depths.empty() && m.jumps.depths.back() == d) { continue; }
            const real_t eps = 1e-6 * std::max(real_t(1.0), d);
            real_t xa[3], xb[3], a[3], b[3];
            DeckPointAtDepth(d - eps, axis, xa);
            DeckPointAtDepth(d + eps, axis, xb);
            dm->eval_at_xyz(xa[0], xa[1], xa[2], a[0], a[1], a[2]);
            dm->eval_at_xyz(xb[0], xb[1], xb[2], b[0], b[1], b[2]);
            bool jump = false;
            for (int c = 0; c < 3; c++)
            {
               jump = jump || std::abs(a[c] - b[c]) > 1e-6 * std::max(std::abs(a[c]), std::abs(b[c]));
            }
            if (jump) { m.jumps.depths.push_back(d); }
         }
         log << "    jumps at depth";
         for (real_t d : m.jumps.depths) { log << " " << d; }
         log << " m (each must be a mesh surface)\n";
         break;
      }
      case spatial::MaterialKind::HalfspaceAcrossFault:
      {
         const auto &hs = cfg.material.halfspace;
         m.halfspace = MakeHalfspaceAcrossFaultMaterial(hs.vp_near, hs.vs_near, hs.rho_near,
                                                        hs.vp_far, hs.vs_far, hs.rho_far,
                                                        hs.x0, hs.normal);
         m.field = m.halfspace->field;
         // mu jumps across the fault plane: no point value there (mu_at empty)
         const real_t nl = std::sqrt(hs.normal[0] * hs.normal[0] + hs.normal[1] * hs.normal[1]
                                     + hs.normal[2] * hs.normal[2]);
         m.jumps.has_plane = true;
         for (int c = 0; c < 3; c++)
         {
            m.jumps.x0[c] = hs.x0[c];
            m.jumps.n[c] = hs.normal[c] / nl;
         }
         log << "  material: halfspace_across_fault, plane x0 (" << hs.x0[0] << ", " << hs.x0[1]
             << ", " << hs.x0[2] << ") n (" << hs.normal[0] << ", " << hs.normal[1] << ", "
             << hs.normal[2] << ") (must be a mesh surface)\n"
             << "    near, sign((x - x0).n) >= 0: vp " << hs.vp_near << ", vs " << hs.vs_near
             << ", rho " << hs.rho_near << ", mu " << hs.rho_near * hs.vs_near * hs.vs_near << "\n"
             << "    far:                         vp " << hs.vp_far << ", vs " << hs.vs_far
             << ", rho " << hs.rho_far << ", mu " << hs.rho_far * hs.vs_far * hs.vs_far << "\n";
         break;
      }
      default:
         MFEM_ABORT("Nitsche driver: [material] kind sidecar_hdf5 is not wired");
   }
   m.log = log.str();
   return m;
}

/// Elements of this rank with vertices strictly on both sides (more than
/// `tol` m) of a jump: [0] across the plane, then one count per jump depth.
inline std::vector<long long> CountElementsAcrossJumps(const Mesh &mesh, const MaterialJumps &J,
                                                       real_t tol)
{
   std::vector<long long> cnt(1 + J.depths.size(), 0);
   Array<int> ev;
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      mesh.GetElementVertices(e, ev);
      real_t smin = 1e300, smax = -1e300, dmin = 1e300, dmax = -1e300;
      for (int a = 0; a < ev.Size(); a++)
      {
         const real_t *x = mesh.GetVertex(ev[a]);
         if (J.has_plane)
         {
            const real_t s = (x[0] - J.x0[0]) * J.n[0] + (x[1] - J.x0[1]) * J.n[1]
                             + (x[2] - J.x0[2]) * J.n[2];
            smin = std::min(smin, s);
            smax = std::max(smax, s);
         }
         const real_t d = DeckDepth(x, J.depth_axis);
         dmin = std::min(dmin, d);
         dmax = std::max(dmax, d);
      }
      if (J.has_plane && smin < -tol && smax > tol) { cnt[0]++; }
      for (size_t k = 0; k < J.depths.size(); k++)
      {
         if (dmin < J.depths[k] - tol && dmax > J.depths[k] + tol) { cnt[1 + k]++; }
      }
   }
   return cnt;
}

/// sigma(x) -> (mu(x) / mu_ref) sigma(x).
inline StressFn ScaleByShearModulus(StressFn sig, ShearModulusFn mu_at, real_t mu_ref)
{
   MFEM_VERIFY(mu_ref > 0.0, "ScaleByShearModulus: mu_ref = " << mu_ref << " must be > 0");
   return [sig, mu_at, mu_ref](const real_t x[3], real_t s[9]) {
      sig(x, s);
      const real_t f = mu_at(x) / mu_ref;
      for (int i = 0; i < 9; i++) { s[i] *= f; }
   };
}

/// TPV31 nucleation, spatial part: the strike overstress
///   tau(x) = delta_tau CosineTaperFactor(r, R, taper) [mu(x) / mu_ref]
/// (the bracket only when spec.mu_ref_pa > 0) as the global tensor
///   dsigma = -tau (shat (x) nhat + nhat (x) shat),
/// so that -(dsigma nhat).shat = +tau (SCEC h-shear-stress).  r is the
/// in-plane distance from the centre measured with (dhat, shat), as for the
/// compact-circular kind.
inline StressFn InstantaneousCircularOverstress(
   const spatial::InstantaneousOverstressCircularSpec &spec, const real_t shat[3],
   const real_t nhat[3], const real_t dhat[3], ShearModulusFn mu_at)
{
   MFEM_VERIFY(spec.radius_m > 0.0 && spec.taper_m >= 0.0, "instantaneous circular nucleation: "
               "radius_m must be > 0 and taper_m >= 0");
   MFEM_VERIFY(spec.mu_ref_pa <= 0.0 || mu_at, "instantaneous circular nucleation: mu_ref_pa > 0 "
               "needs a shear modulus that is single valued on the fault");
   const std::array<real_t, 3> c = {spec.center_x_m, spec.center_y_m, spec.center_z_m};
   const std::array<real_t, 3> sv = {shat[0], shat[1], shat[2]}, nv = {nhat[0], nhat[1], nhat[2]},
                               dv = {dhat[0], dhat[1], dhat[2]};
   return [spec, c, sv, nv, dv, mu_at](const real_t x[3], real_t s[9]) {
      const real_t e[3] = {x[0] - c[0], x[1] - c[1], x[2] - c[2]};
      const real_t ed = e[0] * dv[0] + e[1] * dv[1] + e[2] * dv[2];
      const real_t es = e[0] * sv[0] + e[1] * sv[1] + e[2] * sv[2];
      real_t tau = spec.delta_tau_pa
                   * spatial::CosineTaperFactor(std::sqrt(ed * ed + es * es), spec.radius_m,
                                                spec.taper_m);
      if (spec.mu_ref_pa > 0.0) { tau *= mu_at(x) / spec.mu_ref_pa; }
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++) { s[3 * i + j] = -tau * (sv[i] * nv[j] + nv[i] * sv[j]); }
   };
}

/// Largest |n x n_ref| over this rank's fault faces: the sine of the largest
/// angle between a face's oriented normal and the unit reference normal.
/// Round-off (about 1e-15) on a planar fault.
inline real_t MaxFaultNormalDeviation(const FaultOrientation &orient, const real_t n_ref[3])
{
   const real_t l = std::sqrt(n_ref[0] * n_ref[0] + n_ref[1] * n_ref[1] + n_ref[2] * n_ref[2]);
   const real_t r[3] = {n_ref[0] / l, n_ref[1] / l, n_ref[2] / l};
   real_t m = 0.0;
   for (int fi = 0; fi < orient.NumLocalFaces(); fi++)
   {
      const real_t *n = orient.Normal(fi);
      const real_t c[3] = {n[1] * r[2] - n[2] * r[1], n[2] * r[0] - n[0] * r[2],
                           n[0] * r[1] - n[1] * r[0]};
      m = std::max(m, std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]));
   }
   return m;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_DECK_FIELDS_HPP
