// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// thermal_pressurization.hpp — 2D in-plane design rev 3, phase D (section
// 5.9, derivation 10.5): thermal pressurization as a module of its own, for
// SCEC TPV105 (description eqs 11-12).  At every fault point the shear zone
// is a 1-D problem across the fault (z the fault-normal distance),
//
//    dT/dt = alpha_th d2T/dz2 + tau V / (rho_c h sqrt(2 pi)) exp(-z^2 / (2 h^2)),
//    dp/dt = alpha_hy d2p/dz2 + Lambda dT/dt,
//
// with dT/dz = dp/dz = 0 at z = 0 (symmetry), T -> T_ini and p -> p_ini far
// from the fault, and T = T_ini, p = p_ini at t = 0.  The friction law then
// sees the effective normal stress less the pore-pressure change on the
// fault, p(z = 0) - p_ini (the hook is phase D part 2; this file knows
// nothing of the fault operator).
//
// Discretization (derivation 10.5):
//   grid   z_0 = 0, z_{k+1} = z_k + dz0 r^k up to the first node z_M >= z_max
//          (r = growth >= 1), where T = T_ini and p = p_ini (Dirichlet);
//          unknowns at z_0 .. z_{M-1}, the deviations theta = T - T_ini and
//          pi = p - p_ini
//   space  vertex-centred finite volumes: node k owns [z_{k-1/2}, z_{k+1/2}]
//          (z_{-1/2} = 0, the half cell at the symmetry plane; z_{k+1/2} the
//          midpoint), flux alpha (u_{k+1} - u_k) / (z_{k+1} - z_k); the
//          Gaussian source is integrated exactly over each cell,
//          w_k = (erf(z_{k+1/2} / (sqrt(2) h)) - erf(z_{k-1/2} / (sqrt(2) h))) / 2,
//          so the heat that enters a cell is tau V w_k / rho_c per unit time
//   time   the theta scheme, backward Euler (theta = 1, the default) or
//          Crank-Nicolson (theta = 1/2), with tau V of the step held over the
//          step (one-step lag, as the Prakash-Clifton and Kelvin-Voigt layers
//          of the fault path); T first, then p with the source
//          Lambda (theta^{n+1} - theta^n) / dt, so the undrained limit
//          (alpha_hy = 0) gives pi = Lambda theta exactly; two tridiagonal
//          (Thomas) solves per fault point per step, unconditionally stable
//
// Verification: tests/2d_inplane/test_thermal_pressurization.cpp (the
// constant-source solution of Rice 2006 and its finite-h form, Noda and
// Lapusta 2010; grid and step convergence; heat balance).
//
// Scope (Chunhui Zhao, 2026-10-26): uniform properties along z; one z grid
// for every fault point; per-point parameters (uniform in TPV105, stored per
// point like a and V_w).

#ifndef MFEM_SEAS_THERMAL_PRESSURIZATION_HPP
#define MFEM_SEAS_THERMAL_PRESSURIZATION_HPP

#include "mfem.hpp"

#include <cmath>
#include <vector>

namespace mfem
{
namespace seas
{

/// The TPV105 parameters of one fault point (description, eqs 11-12).
struct TPParameters
{
   real_t alpha_th = 0.0;   ///< thermal diffusivity, m^2/s (> 0)
   real_t alpha_hy = 0.0;   ///< hydraulic diffusivity, m^2/s (> 0)
   real_t rho_c = 0.0;      ///< volumetric heat capacity, J/(m^3 K) (> 0)
   real_t Lambda = 0.0;     ///< pore-pressure rise per kelvin, Pa/K (>= 0)
   real_t h = 0.0;          ///< shear-zone width: the source's standard deviation, m (> 0)
   real_t T_ini = 0.0;      ///< initial and far-field temperature, K
   real_t p_ini = 0.0;      ///< initial and far-field pore pressure, Pa
};

/// The z grid: z_0 = 0, z_{k+1} - z_k = dz0 growth^k, up to the first node at
/// or beyond z_max (the far-field node).
struct TPGridSpec
{
   real_t dz0 = 0.0;      ///< first spacing, m (> 0)
   real_t growth = 1.0;   ///< spacing ratio (>= 1)
   real_t z_max = 0.0;    ///< m (> dz0)
};

enum class TPScheme { BackwardEuler, CrankNicolson };

class ThermalPressurization
{
public:
   /// npts fault points (0 is allowed: a rank without fault points), every
   /// one at T_ini and p_ini of default (zero) parameters until
   /// SetParameters.
   ThermalPressurization(int npts, const TPGridSpec &grid,
                         TPScheme scheme = TPScheme::BackwardEuler)
      : n_(npts), theta_(scheme == TPScheme::BackwardEuler ? 1.0 : 0.5)
   {
      MFEM_VERIFY(npts >= 0, "ThermalPressurization: npts = " << npts);
      MFEM_VERIFY(grid.dz0 > 0.0 && grid.growth >= 1.0 && grid.z_max > grid.dz0,
                  "ThermalPressurization: grid dz0 = " << grid.dz0 << " m, growth = "
                  << grid.growth << ", z_max = " << grid.z_max << " m; need dz0 > 0, "
                  "growth >= 1 and z_max > dz0");
      const int kMaxNodes = 100000;
      z_.push_back(0.0);
      real_t dz = grid.dz0;
      while (z_.back() < grid.z_max)
      {
         z_.push_back(z_.back() + dz);
         dz *= grid.growth;
         MFEM_VERIFY((int) z_.size() <= kMaxNodes, "ThermalPressurization: more than "
                     << kMaxNodes << " nodes to reach z_max = " << grid.z_max
                     << " m from dz0 = " << grid.dz0 << " m with growth " << grid.growth);
      }
      m_ = (int) z_.size() - 1;   // unknowns 0 .. m_-1; z_[m_] is the far-field node
      vol_.resize(m_);
      for (int k = 0; k < m_; k++)
      {
         const real_t lo = (k == 0) ? 0.0 : 0.5 * (z_[k - 1] + z_[k]);
         vol_[k] = 0.5 * (z_[k] + z_[k + 1]) - lo;
      }
      prm_.resize(n_);
      w_.assign((size_t) n_ * m_, 0.0);
      T_.assign((size_t) n_ * m_, 0.0);
      p_.assign((size_t) n_ * m_, 0.0);
      dT_.resize(m_);
      a_.resize(m_); b_.resize(m_); c_.resize(m_); r_.resize(m_);
   }

   /// Parameters of point j; its T and p go back to T_ini and p_ini.
   void SetParameters(int j, const TPParameters &prm)
   {
      MFEM_VERIFY(j >= 0 && j < n_, "ThermalPressurization::SetParameters: point "
                  << j << " of " << n_);
      MFEM_VERIFY(prm.alpha_th > 0.0 && prm.alpha_hy > 0.0 && prm.rho_c > 0.0
                  && prm.Lambda >= 0.0 && prm.h > 0.0 && std::isfinite(prm.T_ini)
                  && std::isfinite(prm.p_ini),
                  "ThermalPressurization: point " << j << ": alpha_th = " << prm.alpha_th
                  << ", alpha_hy = " << prm.alpha_hy << ", rho_c = " << prm.rho_c
                  << ", Lambda = " << prm.Lambda << ", h = " << prm.h << ", T_ini = "
                  << prm.T_ini << ", p_ini = " << prm.p_ini << "; need alpha_th, alpha_hy, "
                  "rho_c, h > 0, Lambda >= 0, T_ini and p_ini finite");
      prm_[j] = prm;
      const real_t s = 1.0 / (std::sqrt(2.0) * prm.h);
      real_t lo = 0.0;   // erf at z_{k-1/2}
      for (int k = 0; k < m_; k++)
      {
         const real_t hi = std::erf(0.5 * (z_[k] + z_[k + 1]) * s);
         w_[(size_t) j * m_ + k] = 0.5 * (hi - lo);
         lo = hi;
      }
      for (int k = 0; k < m_; k++)
      {
         T_[(size_t) j * m_ + k] = 0.0;
         p_[(size_t) j * m_ + k] = 0.0;
      }
   }

   /// The same parameters at every point.
   void SetParameters(const TPParameters &prm)
   {
      for (int j = 0; j < n_; j++) { SetParameters(j, prm); }
   }

   /// One step of length dt: point j heated by tau[j] V[j] (shear traction
   /// magnitude in Pa times slip-rate magnitude in m/s, >= 0) over the step.
   void Advance(real_t dt, const real_t *tau, const real_t *V)
   {
      MFEM_VERIFY(dt > 0.0 && std::isfinite(dt), "ThermalPressurization::Advance: dt = "
                  << dt);
      for (int j = 0; j < n_; j++)
      {
         const real_t q = tau[j] * V[j];
         MFEM_VERIFY(q >= 0.0 && std::isfinite(q), "ThermalPressurization::Advance: point "
                     << j << ": tau V = " << tau[j] << " Pa x " << V[j] << " m/s = " << q
                     << "; the shear heating must be finite and >= 0");
         const TPParameters &P = prm_[j];
         MFEM_VERIFY(P.rho_c > 0.0, "ThermalPressurization::Advance: point " << j
                     << " has no parameters (SetParameters)");
         real_t *T = &T_[(size_t) j * m_];
         real_t *p = &p_[(size_t) j * m_];
         const real_t *w = &w_[(size_t) j * m_];
         // temperature: the heat q w_k / rho_c enters cell k per unit time
         const real_t qc = q / P.rho_c;
         for (int k = 0; k < m_; k++) { dT_[k] = T[k]; }
         Assemble_(P.alpha_th, dt, T);
         for (int k = 0; k < m_; k++) { r_[k] += qc * w[k]; }
         Solve_(T);
         // pressure: the source Lambda dT/dt over the step
         for (int k = 0; k < m_; k++) { dT_[k] = T[k] - dT_[k]; }
         Assemble_(P.alpha_hy, dt, p);
         for (int k = 0; k < m_; k++) { r_[k] += P.Lambda * vol_[k] * dT_[k] / dt; }
         Solve_(p);
      }
   }

   int NumPoints() const { return n_; }
   /// Unknowns per point (the far-field node z_M not counted).
   int NumZ() const { return m_; }
   /// The nodes z_0 = 0 .. z_M (M + 1 values; the last is the far-field node).
   const std::vector<real_t> &Z() const { return z_; }
   /// Control-volume lengths of the M unknowns.
   const std::vector<real_t> &Volumes() const { return vol_; }
   const TPParameters &Parameters(int j) const { return prm_[j]; }

   /// T and p on the fault (z = 0) and their changes since t = 0.
   real_t Temperature(int j) const { return prm_[j].T_ini + T_[(size_t) j * m_]; }
   real_t PorePressure(int j) const { return prm_[j].p_ini + p_[(size_t) j * m_]; }
   real_t TemperatureChange(int j) const { return T_[(size_t) j * m_]; }
   real_t PorePressureChange(int j) const { return p_[(size_t) j * m_]; }
   /// The changes at node k (k < NumZ()).
   real_t TemperatureChange(int j, int k) const { return T_[(size_t) j * m_ + k]; }
   real_t PorePressureChange(int j, int k) const { return p_[(size_t) j * m_ + k]; }
   /// rho_c times the integral of T - T_ini over the half line (J/m^2).
   real_t HeatContent(int j) const
   {
      real_t s = 0.0;
      for (int k = 0; k < m_; k++) { s += vol_[k] * T_[(size_t) j * m_ + k]; }
      return prm_[j].rho_c * s;
   }
   /// The fraction of a point's heating that falls inside the grid
   /// (sum of w_k = erf(z_{M-1/2} / (sqrt(2) h)) / 2; 1/2 is the half line).
   real_t SourceFraction(int j) const
   {
      real_t s = 0.0;
      for (int k = 0; k < m_; k++) { s += w_[(size_t) j * m_ + k]; }
      return s;
   }

   /// The state for a checkpoint: the 2 N M deviations T - T_ini, p - p_ini.
   void GetState(Vector &s) const
   {
      s.SetSize(2 * n_ * m_);
      for (size_t i = 0; i < T_.size(); i++) { s(i) = T_[i]; s(T_.size() + i) = p_[i]; }
   }
   void SetState(const Vector &s)
   {
      MFEM_VERIFY(s.Size() == 2 * n_ * m_, "ThermalPressurization::SetState: size "
                  << s.Size() << " != 2 x " << n_ << " points x " << m_ << " nodes");
      for (size_t i = 0; i < T_.size(); i++) { T_[i] = s(i); p_[i] = s(T_.size() + i); }
   }

private:
   /// The theta-scheme system for diffusivity alpha, step dt and the old
   /// values u: (vol/dt + theta A) u^{n+1} = (vol/dt - (1 - theta) A) u^n,
   /// A the finite-volume diffusion operator with u_M = 0 (the source is
   /// added to r_ by the caller).
   void Assemble_(real_t alpha, real_t dt, const real_t *u)
   {
      const real_t th = theta_;
      for (int k = 0; k < m_; k++)
      {
         const real_t lo = (k == 0) ? 0.0 : alpha / (z_[k] - z_[k - 1]);
         const real_t up = alpha / (z_[k + 1] - z_[k]);
         a_[k] = -th * lo;
         b_[k] = vol_[k] / dt + th * (lo + up);
         c_[k] = -th * up;
         real_t Au = (lo + up) * u[k];
         if (k > 0) { Au -= lo * u[k - 1]; }
         if (k + 1 < m_) { Au -= up * u[k + 1]; }
         r_[k] = vol_[k] / dt * u[k] - (1.0 - th) * Au;
      }
   }

   /// Thomas algorithm on (a_, b_, c_) with right side r_; the solution in u.
   void Solve_(real_t *u)
   {
      for (int k = 1; k < m_; k++)
      {
         const real_t f = a_[k] / b_[k - 1];
         b_[k] -= f * c_[k - 1];
         r_[k] -= f * r_[k - 1];
      }
      u[m_ - 1] = r_[m_ - 1] / b_[m_ - 1];
      for (int k = m_ - 2; k >= 0; k--) { u[k] = (r_[k] - c_[k] * u[k + 1]) / b_[k]; }
   }

   int n_ = 0;                  ///< fault points
   int m_ = 0;                  ///< unknowns per point
   real_t theta_ = 1.0;
   std::vector<real_t> z_;      ///< m_ + 1 nodes
   std::vector<real_t> vol_;    ///< m_ control volumes
   std::vector<TPParameters> prm_;
   std::vector<real_t> w_;      ///< n_ x m_ source weights
   std::vector<real_t> T_, p_;  ///< n_ x m_ deviations from T_ini, p_ini
   std::vector<real_t> dT_, a_, b_, c_, r_;   ///< work, m_ each
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_THERMAL_PRESSURIZATION_HPP
