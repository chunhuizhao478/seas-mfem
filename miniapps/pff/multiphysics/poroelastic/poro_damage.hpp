// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_damage.hpp — damage-dependent poroelastic properties of the
// cmame_revision hydromechanical model (u-p Biot, fully saturated, one phase).
//
// MOOSE/FARMS sources (farms branch cdms, src/materials/porousflowmatprops and
// src/materials/phasefield_smalldeform_nonAD):
//   Biot coefficient  ElkPorousFlowDamagedBiotCoefficient
//                     g2 = max((1-d)^2, g_min); alpha = clamp(1 - g2 K / K_s, 0, 1)
//   porosity          ElkPorousFlowDamagedPorosity
//                     damage: phi = phi0 + (1 - phi0)(1 - (1-d)^2)
//                     strain: phi = phi0 + eps_1 (max principal mechanical strain)
//                     then clamp(phi, phi_lb, phi_ub)
//   Biot modulus      ElkPorousFlowDamagedBiotModulus
//                     1/M = max(phi/K_f + (alpha - phi)/K_s, 1e-20)
//   history energy    ElkPorousFlowHistoryEnergyEnhanced
//                     c_p = (phi0 - 1)(1/K_f - 1/K_s) - K/K_s^2  (K intact)
//                     H = max(H_old, <2 psi_e^+ + c_p p^2>_+)
//   permeability      NDSmallDeformationIsotropicElasticity::updatePermeabilityForCracking
//                     exponential      k0 exp(b d) I
//                     darcy_poiseuille k0 I + d^n (w^2/12 I - k0 I), w = d w_c
//                     normal_strain    k0 I + d^n chi_d k_w (I - n n) (anisotropic)
//                                      chi_d = H(d - d_th), w_c = h_c max(1 + n.eps.n, 0),
//                                      w_h = max(f_c w_c, f_c w_r), k_w = w_h^2 / 12
//                     (the rotation applied to the first two is a no-op on an
//                     isotropic tensor and is omitted)
//
// Field equations these feed (ElkPorousFlowEffectiveStressCoupling,
// ElkPorousFlowFullySaturatedMassTimeDerivative, PorousFlowFullySaturatedDarcyBase,
// multiply_by_density = false, relative permeability 1, no gravity):
//   div(sigma_eff(u, d) - alpha(d) p I) = rho u_tt
//   (p - p_old)/(dt M) + alpha (eps_v - eps_v_old)/dt - div(k/mu grad p) = 0

#ifndef MFEM_PFF_PORO_DAMAGE_HPP
#define MFEM_PFF_PORO_DAMAGE_HPP

#include "mfem.hpp"
#include "../../../../linalg/kernels.hpp" // mfem::kernels::CalcEigenvalues<3>
#include "../../common/pff_types.hpp"

#include <algorithm>
#include <cmath>

namespace mfem
{
namespace pff
{

struct PoroDamageParams
{
   real_t K = 0.0;              ///< drained bulk modulus of the intact skeleton [Pa]
   real_t K_s = 0.0;            ///< grain bulk modulus [Pa]
   real_t K_f = 0.0;            ///< fluid bulk modulus [Pa]
   real_t phi0 = 0.0;           ///< initial porosity
   real_t phi_lb = 0.0;         ///< porosity lower bound
   real_t phi_ub = 1.0;         ///< porosity upper bound
   real_t g_min = 1e-6;         ///< floor of (1-d)^2 in the Biot coefficient
   real_t viscosity = 1e-3;     ///< fluid dynamic viscosity [Pa s]
   PorosityModel porosity_model = PorosityModel::Damage;

   // permeability
   PermeabilityModel perm_model = PermeabilityModel::None;
   real_t k0 = 0.0;             ///< intrinsic permeability [m^2]
   real_t perm_exponent = 1.0;  ///< n in d^n (darcy_poiseuille, normal_strain)
   real_t coeff_b = 0.0;        ///< b of the exponential model
   real_t wc = 0.0;             ///< ultimate width of darcy_poiseuille [m]
   real_t d_threshold = 0.5;    ///< chi_d threshold (normal_strain)
   real_t fc = 1.0;             ///< roughness correction f_c
   real_t w_res = 0.0;          ///< residual aperture w_r [m]
   bool anisotropic = true;     ///< K_frac = k_w (I - n n) or k_w I
   CrackNormalSource normal_source = CrackNormalSource::DamageGradient;
   bool regularize_normal = false;
   real_t normal_reg_eps = 1e-8;
   real_t grad_d_tol = 1e-30;
};

inline real_t Clamp(real_t x, real_t lo, real_t hi) { return std::min(std::max(x, lo), hi); }

/// Damaged Biot coefficient alpha(d).
inline real_t BiotCoefficient(const PoroDamageParams &P, real_t d)
{
   const real_t dc = Clamp(d, 0.0, 1.0);
   const real_t g2 = std::max((1.0 - dc) * (1.0 - dc), P.g_min);
   return Clamp(1.0 - g2 * P.K / P.K_s, 0.0, 1.0);
}

/// Damaged porosity; eps_max is the max principal mechanical strain (strain model).
inline real_t Porosity(const PoroDamageParams &P, real_t d, real_t eps_max = 0.0)
{
   real_t phi;
   if (P.porosity_model == PorosityModel::Strain)
   {
      phi = P.phi0 + eps_max;
   }
   else
   {
      const real_t dc = Clamp(d, 0.0, 1.0);
      const real_t g = (1.0 - dc) * (1.0 - dc);
      phi = P.phi0 + (1.0 - P.phi0) * (1.0 - g);
   }
   return Clamp(phi, P.phi_lb, P.phi_ub);
}

/// Biot modulus M(alpha, phi).
inline real_t BiotModulus(const PoroDamageParams &P, real_t alpha, real_t phi)
{
   const real_t denom = phi / P.K_f + (alpha - phi) / P.K_s;
   return 1.0 / std::max(denom, real_t(1e-20));
}

/// Pressure coefficient c_p of the enhanced history energy.
inline real_t HistoryPressureCoefficient(const PoroDamageParams &P)
{
   return (P.phi0 - 1.0) * (1.0 / P.K_f - 1.0 / P.K_s) - P.K / (P.K_s * P.K_s);
}

/// One history update H = max(H_old, <2 psi+ + c_p p^2>_+).
inline real_t EnhancedHistory(const PoroDamageParams &P, real_t psi_active,
                              real_t p, real_t H_old)
{
   const real_t Hc = 2.0 * psi_active + p * p * HistoryPressureCoefficient(P);
   return std::max(std::max(Hc, real_t(0.0)), H_old);
}

/// Effective permeability tensor (3x3 column-major) at one point.
///  d       damage,  grad_d  damage gradient (3),  eps  total strain (3x3),
///  h_c     characteristic length (element size for the cmame inputs).
inline void EffectivePermeability(const PoroDamageParams &P, real_t d,
                                  const real_t *grad_d, const real_t *eps,
                                  real_t h_c, real_t *k)
{
   for (int i = 0; i < 9; i++) { k[i] = 0.0; }
   auto set_iso = [&](real_t s) { k[0] = k[4] = k[8] = s; };

   switch (P.perm_model)
   {
      case PermeabilityModel::None:
         set_iso(P.k0);
         return;
      case PermeabilityModel::Exponential:
         set_iso(P.k0 * std::exp(d * P.coeff_b));
         return;
      case PermeabilityModel::DarcyPoiseuille:
      {
         const real_t w = d * P.wc;
         const real_t kf = w * w / 12.0;
         set_iso(P.k0 + std::pow(d, P.perm_exponent) * (kf - P.k0));
         return;
      }
      case PermeabilityModel::NormalStrain:
      {
         real_t n[3] = {0.0, 0.0, 0.0};
         bool have_normal = false;
         if (P.normal_source == CrackNormalSource::DamageGradient)
         {
            const real_t gn = std::sqrt(grad_d[0]*grad_d[0] + grad_d[1]*grad_d[1] +
                                        grad_d[2]*grad_d[2]);
            if (P.regularize_normal)
            {
               for (int i = 0; i < 3; i++) { n[i] = grad_d[i] / (gn + P.normal_reg_eps); }
               have_normal = true;
            }
            else if (gn > P.grad_d_tol)
            {
               for (int i = 0; i < 3; i++) { n[i] = grad_d[i] / gn; }
               have_normal = true;
            }
         }
         else
         {
            real_t lam[3], vec[9];
            kernels::CalcEigenvalues<3>(eps, lam, vec);
            int imax = 0;
            for (int i = 1; i < 3; i++) { if (lam[i] > lam[imax]) { imax = i; } }
            for (int i = 0; i < 3; i++) { n[i] = vec[3*imax + i]; }
            have_normal = lam[imax] > 0.0;
         }
         const real_t chi = (d >= P.d_threshold) ? 1.0 : 0.0;
         if (!have_normal || chi == 0.0 || d <= 0.0 || h_c <= 0.0)
         {
            set_iso(P.k0);
            return;
         }
         real_t enn = 0.0;
         for (int j = 0; j < 3; j++)
            for (int i = 0; i < 3; i++) { enn += n[i] * eps[i + 3*j] * n[j]; }
         const real_t wcap = h_c * std::max(1.0 + enn, real_t(0.0));
         const real_t wh = std::max(P.fc * wcap * chi, P.fc * P.w_res * chi);
         const real_t kw = wh * wh / 12.0;
         const real_t wt = std::pow(d, P.perm_exponent);
         for (int j = 0; j < 3; j++)
            for (int i = 0; i < 3; i++)
            {
               const real_t dij = (i == j) ? 1.0 : 0.0;
               const real_t frac = P.anisotropic ? kw * (dij - n[i] * n[j]) : kw * dij;
               k[i + 3*j] = P.k0 * dij + wt * frac;
            }
         return;
      }
   }
}

} // namespace pff
} // namespace mfem

#endif
