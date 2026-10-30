// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// finite_width.hpp — the 2.5D finite-width model of a map-view cycle
// (Chunhui Zhao, 2026-10-30: "yes send this paper and equation to 2D thread
// and start implementation and tests").  Derivation, calibration and the
// consistency rules C1-C6:
// /mnt/project-files/2d-inplane/nonplanar-literature/inplane-2p5d-derivation-rev8.pdf
// (Weng & Ampuero 2019, 2020; Lehner, Li & Rice 1981).
//
//   rho u'' = div sigma(u) - beta (u - u_pl),   beta = 4 c^2 mu / W^2
//
// u the in-plane displacement averaged over the seismogenic width W.  On a
// long straight strike-slip fault the static field decays over l = W / (2 c)
// and the stress drop is dtau = c mu D / W (D the width-averaged slip); a
// locked straight fault is loaded at c mu Vp / W.
//
// [cycle.finite_width] (cycle_config.hpp, FiniteWidthConfig):
//   W_m          the seismogenic width W [m]                      [required]
//   calibration  c by name, buried faults only (the note's Table 1;
//                Chunhui Zhao, 2026-10-30: "we will for 2d in plane, only
//                focus on buried fault for now")                 ["wa2020"]
//                  "wa2020"            (1 - nu) pi / 2: Weng & Ampuero 2020
//                                      eq. 21, fitted to 3D dynamic runs;
//                                      sem2dpack's beta = mu (pi (1 - nu) / W)^2
//                  "lehner"            4 / pi: the exact average slip of a
//                                      buried long fault (Lehner, Li & Rice
//                                      1981, gamma = pi / 8)
//   c            the factor as a number instead of calibration
//   substrate    u_pl(x, t) = t u_pl(x, 1)                        [required]
//                  "loading_plane"     +-(Vp t / 2) s on the two sides of
//                                      [cycle.loading]'s plane (kind "bp5":
//                                      y = 0, s = x): the note's eq. 9 with
//                                      the plane as the fault system, for
//                                      straight faults on the plane
//                  "locked_reference"  the QD operator's static displacement
//                                      under the Dirichlet loading at t = 1
//                                      with zero slip and beta = 0: while
//                                      the faults are locked the run is
//                                      loaded as without the term, and the
//                                      term acts on slip alone.  A test
//                                      setting (bend and hand-back checks
//                                      on the same events as the plane-
//                                      strain runs): over many cycles it
//                                      pulls accumulated slip back.
// beta is isotropic, beta(x) = (4 c^2 / W^2) mu(x) on both in-plane
// components (Chunhui Zhao, 2026-10-30: "we will start with beta = beta
// (normal direction), however this could be important"; the note's Step 7:
// with "wa2020" the fault-normal stiffness sqrt(M beta) / 2 is about 1.20 x
// the 3D mode I value, a known model error).  Map view only (a
// 2-D mesh with [fault_geometry] up = (0, 0, +-1)); not with [cycle.lts]
// rate2.
//
// Both regimes take the term through ElasticityDomainOperator::
// SetFiniteWidth with one U_pl1 (the interpolant of u_pl(x, 1)), built once
// by QDPhase and handed to DynPhase: one M_beta and one u_pl(t) in the QD
// solve, G, the dynamic step and the hand-back (rules C1, C2; the dynamic
// regime reads u_pl at t0 + tau).  Kelvin-Voigt damps K alone (C3; the one
// choice, stated in the run log).  The explicit step is cut to dt / sqrt(1 +
// beta dt^2 / (4 rho)) (C6).  Energy bookkeeping (C5) is not part of the
// cycle.

#ifndef MFEM_SEAS_CYCLE_FINITE_WIDTH_HPP
#define MFEM_SEAS_CYCLE_FINITE_WIDTH_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "../domain/boundary_config.hpp"

#include <array>
#include <cmath>
#include <sstream>
#include <string>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// The "loading_plane" substrate: u = sgn(d) (Vp t / 2) s, d = (x - p) . m,
/// with [cycle.loading]'s p, m, s ("bp5": p = 0, m = e_y, s = e_x), no band:
/// the far-field values of MakePlateLoading on both sides, up to the plane.
inline DirichletFunc MakePlateSubstrate(const LoadingConfig &L, real_t Vp)
{
   const bool plate = L.Plate();
   const std::array<real_t, 3> p = plate ? L.plane_point : std::array<real_t, 3> {{0.0, 0.0, 0.0}};
   const std::array<real_t, 3> m = plate ? L.plane_normal : std::array<real_t, 3> {{0.0, 1.0, 0.0}};
   const std::array<real_t, 3> s = plate ? L.slip_dir : std::array<real_t, 3> {{1.0, 0.0, 0.0}};
   return [Vp, p, m, s](const Vector &x, real_t t, Vector &u)
   {
      const int dim = x.Size();
      u.SetSize(dim);
      u = 0.0;
      real_t d = 0.0;
      for (int i = 0; i < dim; i++) { d += (x(i) - p[i]) * m[i]; }
      const real_t h = (d > 0.0) ? 0.5 * Vp * t : (d < 0.0) ? -0.5 * Vp * t : 0.0;
      for (int i = 0; i < dim; i++)
      {
         if (s[i] != 0.0) { u(i) = h * s[i]; }
      }
   };
}

/// The interpolant of u(x, 1) on the vector DG space `fes`, each element
/// evaluated from its inside: a node x is evaluated at x + 1e-6 (x_c - x),
/// x_c the element's centre, so a node on a jump of u (a fault on the
/// loading plane) takes its own element's side.
inline void InterpolateSubstrate(ParFiniteElementSpace &fes, const DirichletFunc &u, Vector &U)
{
   class Inside : public VectorCoefficient
   {
   public:
      Inside(int vd, const DirichletFunc &f) : VectorCoefficient(vd), f_(f) {}
      void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
      {
         Vector x, xc;
         T.Transform(ip, x);
         T.Transform(Geometries.GetCenter(T.GetGeometryType()), xc);
         for (int i = 0; i < x.Size(); i++) { x(i) += 1e-6 * (xc(i) - x(i)); }
         f_(x, 1.0, V);
      }
   private:
      const DirichletFunc &f_;
   };
   ParGridFunction g(&fes);
   Inside c(fes.GetVDim(), u);
   g.ProjectCoefficient(c);
   U = g;
}

/// One paragraph for the run log: W, c, beta, l, the wave cutoff, the
/// locked-fault stressing rate, the substrate.
inline std::string FiniteWidthReport(const CycleConfig &cfg)
{
   const FiniteWidthConfig &F = cfg.finite_width;
   std::ostringstream o;
   if (!F.On()) { o << "off"; return o.str(); }
   const auto &mf = cfg.shared.material_fallback;
   const real_t beta = F.BetaOverMu() * mf.mu;
   o << "W = " << F.W << " m, c = " << F.c
     << (F.calibration.empty() ? std::string(" (deck)") : " (" + F.calibration + ")")
     << ", beta = " << beta << " Pa/m^2 (isotropic), decay length " << F.DecayLength()
     << " m, wave cutoff " << std::sqrt(beta / mf.rho) / (2.0 * M_PI) << " Hz; substrate \""
     << F.substrate << "\"";
   if (F.substrate == "loading_plane")
   {
      o << " (a locked straight fault on the plane is loaded at "
        << F.c * mf.mu * cfg.plate_rate / F.W * 3.15576e7 / 1e3 << " kPa/yr)";
   }
   else { o << " (the loading of the run without the term while the faults are locked)"; }
   return o.str();
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_FINITE_WIDTH_HPP
