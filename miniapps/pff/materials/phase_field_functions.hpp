// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// phase_field_functions.hpp — degradation g(d), crack geometric function
// alpha(d), normalisation c0 and the pointwise damage source.
//
// MOOSE/FARMS sources (farms branch cdms):
//   g      PowerDegradationFunction, expression (1-d)^p*(1-eta)+eta, p = 2,
//          eta = 1e-6 in every spe_hftc and cmame_revision input.
//          NDSmallDeformationIsotropicElasticity (model_type AT1/AT2) hard-codes
//          the same p = 2 form.
//   alpha  CrackGeometricFunction: AT1 alpha = d, AT2 alpha = d^2;
//          c0 = 4 int_0^1 sqrt(alpha(s)) ds = 8/3 (AT1), 2 (AT2).
//   source ADPFFSource returns d psi / d d with
//          psi = alpha Gc/(c0 l) + g D            (spe_hftc, D = psie_active)
//          psi = alpha Gc/(c0 l) + 0.5 g H        (cmame_revision, D = H/2)
//   diffusion ADPFFDiffusion: 2 Gc l / c0 grad d . grad q.
//
// All functions are MFEM_HOST_DEVICE so the quadrature-point kernels can run
// on a device later.

#ifndef MFEM_PFF_PHASE_FIELD_FUNCTIONS_HPP
#define MFEM_PFF_PHASE_FIELD_FUNCTIONS_HPP

#include "mfem.hpp"
#include "../common/pff_types.hpp"

#include <cmath>

namespace mfem
{
namespace pff
{

/// Power degradation g(d) = (1-d)^p (1-eta) + eta (PowerDegradationFunction).
struct PowerDegradation
{
   real_t p = 2.0;
   real_t eta = 1e-6;

   MFEM_HOST_DEVICE real_t g(real_t d) const
   { return std::pow(1.0 - d, p) * (1.0 - eta) + eta; }

   MFEM_HOST_DEVICE real_t dg(real_t d) const
   { return -p * std::pow(1.0 - d, p - 1.0) * (1.0 - eta); }

   MFEM_HOST_DEVICE real_t d2g(real_t d) const
   { return p * (p - 1.0) * std::pow(1.0 - d, p - 2.0) * (1.0 - eta); }
};

/// Crack geometric function alpha(d) and its normalisation c0.
struct CrackGeometric
{
   CrackModel model = CrackModel::AT1;

   MFEM_HOST_DEVICE real_t alpha(real_t d) const
   { return model == CrackModel::AT1 ? d : d * d; }

   MFEM_HOST_DEVICE real_t dalpha(real_t d) const
   { return model == CrackModel::AT1 ? 1.0 : 2.0 * d; }

   MFEM_HOST_DEVICE real_t d2alpha(real_t /*d*/) const
   { return model == CrackModel::AT1 ? 0.0 : 2.0; }

   /// c0 = 4 int_0^1 sqrt(alpha(s)) ds.
   MFEM_HOST_DEVICE real_t c0() const
   { return model == CrackModel::AT1 ? 8.0 / 3.0 : 2.0; }
};

/// The phase-field model of one case: toughness, length, g, alpha and the
/// convention for the driving energy.
struct PhaseFieldModel
{
   real_t Gc = 0.0;          ///< fracture toughness [J/m^2]
   real_t l = 0.0;           ///< regularisation length [m]
   PowerDegradation deg;
   CrackGeometric crack;
   DrivingEnergy driving = DrivingEnergy::Active;

   /// D of the damage source from the quantity the elasticity side provides:
   /// psi_e^+ (Active) or the enhanced history H (HalfHistoryEnhanced).
   MFEM_HOST_DEVICE real_t DrivingFrom(real_t q) const
   { return driving == DrivingEnergy::Active ? q : 0.5 * q; }

   /// Pointwise source d psi/d d = Gc/(c0 l) alpha'(d) + g'(d) D (ADPFFSource).
   MFEM_HOST_DEVICE real_t Source(real_t d, real_t D) const
   { return Gc / (crack.c0() * l) * crack.dalpha(d) + deg.dg(d) * D; }

   /// d(Source)/d d, for a Newton linearisation (zero D-derivative).
   MFEM_HOST_DEVICE real_t SourceDerivative(real_t d, real_t D) const
   { return Gc / (crack.c0() * l) * crack.d2alpha(d) + deg.d2g(d) * D; }

   /// Coefficient of grad d . grad q in the damage residual (ADPFFDiffusion).
   MFEM_HOST_DEVICE real_t DiffusionCoefficient() const
   { return 2.0 * Gc * l / crack.c0(); }

   /// Crack energy density Gc/c0 (alpha/l + l |grad d|^2)
   /// (CrackDissipatedEnergyDensity).
   MFEM_HOST_DEVICE real_t CrackEnergyDensity(real_t d, real_t grad_d_sq) const
   { return Gc / crack.c0() * (crack.alpha(d) / l + l * grad_d_sq); }

   /// AT1 elastic threshold: the D at which the source vanishes at d = 0,
   /// D_c = Gc / (c0 l |g'(0)|).  Zero for AT2 (no threshold).
   real_t OnsetDrivingEnergy() const
   {
      if (crack.model == CrackModel::AT2) { return 0.0; }
      return Gc / (crack.c0() * l * std::abs(deg.dg(0.0)));
   }

   /// AT1 uniaxial strength sqrt(3 E Gc / (8 l)) for the Active convention
   /// (the relation used to set Gc from ft in spe_hftc).
   static real_t AT1Strength(real_t E, real_t Gc, real_t l)
   { return std::sqrt(3.0 * E * Gc / (8.0 * l)); }

   /// Inverse of AT1Strength: Gc = 8 l ft^2 / (3 E).
   static real_t AT1ToughnessFromStrength(real_t E, real_t ft, real_t l)
   { return 8.0 * l * ft * ft / (3.0 * E); }
};

} // namespace pff
} // namespace mfem

#endif
