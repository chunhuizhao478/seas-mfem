// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// pff_types.hpp — shared enums of the dynamic phase-field fracture miniapp.
//
// The names follow the MOOSE/FARMS options they reproduce (farms branch cdms),
// so that a case file can be read side by side with the MOOSE input:
//   Decomposition        SmallDeformationIsotropicElasticity::decomposition
//   CrackModel           CrackGeometricFunction (alpha(d), c0)
//   DrivingEnergy        which energy enters g'(d) in the damage equation
//   PermeabilityModel    NDSmallDeformationIsotropicElasticity::permeability_model
//   PorosityModel        ElkPorousFlowDamagedPorosity::porosity_update_model

#ifndef MFEM_PFF_TYPES_HPP
#define MFEM_PFF_TYPES_HPP

#include <string>
#include <stdexcept>

namespace mfem
{
namespace pff
{

/// Tension/compression split of the elastic energy.
enum class Decomposition { None, Spectral, VolDev };

/// Crack geometric function: AT1 alpha = d (c0 = 8/3), AT2 alpha = d^2 (c0 = 2).
enum class CrackModel { AT1, AT2 };

/// Energy that drives damage, D in  Gc/(c0 l) alpha'(d) + g'(d) D.
///  - Active:          D = psi_e^+ (spe_hftc fracture.i: psi = alpha Gc/c0/l + g psie_active)
///  - HalfHistoryEnhanced: D = 1/2 H, H = max_t <2 psi_e^+ + c_p p^2>_+
///    (cmame_revision hydromech: psi = alpha Gc/c0/l + 0.5 g psie_active, the
///    transferred psie_active being ElkPorousFlowHistoryEnergyEnhanced's H)
enum class DrivingEnergy { Active, HalfHistoryEnhanced };

enum class PermeabilityModel { None, Exponential, DarcyPoiseuille, NormalStrain };

enum class PorosityModel { Damage, Strain };

/// Source of the crack normal in the normal-strain permeability model.
enum class CrackNormalSource { DamageGradient, PrincipalStrain };

inline Decomposition ParseDecomposition(const std::string &s)
{
   if (s == "none") { return Decomposition::None; }
   if (s == "spectral") { return Decomposition::Spectral; }
   if (s == "voldev") { return Decomposition::VolDev; }
   throw std::invalid_argument("unknown decomposition '" + s + "'");
}

inline CrackModel ParseCrackModel(const std::string &s)
{
   if (s == "AT1") { return CrackModel::AT1; }
   if (s == "AT2") { return CrackModel::AT2; }
   throw std::invalid_argument("unknown crack model '" + s + "'");
}

inline DrivingEnergy ParseDrivingEnergy(const std::string &s)
{
   if (s == "active") { return DrivingEnergy::Active; }
   if (s == "half_history_enhanced") { return DrivingEnergy::HalfHistoryEnhanced; }
   throw std::invalid_argument("unknown driving energy '" + s + "'");
}

inline PermeabilityModel ParsePermeabilityModel(const std::string &s)
{
   if (s == "none") { return PermeabilityModel::None; }
   if (s == "exponential") { return PermeabilityModel::Exponential; }
   if (s == "darcy_poiseuille") { return PermeabilityModel::DarcyPoiseuille; }
   if (s == "normal_strain") { return PermeabilityModel::NormalStrain; }
   throw std::invalid_argument("unknown permeability model '" + s + "'");
}

inline PorosityModel ParsePorosityModel(const std::string &s)
{
   if (s == "damage") { return PorosityModel::Damage; }
   if (s == "strain") { return PorosityModel::Strain; }
   throw std::invalid_argument("unknown porosity model '" + s + "'");
}

inline CrackNormalSource ParseCrackNormalSource(const std::string &s)
{
   if (s == "damage_gradient") { return CrackNormalSource::DamageGradient; }
   if (s == "principal_strain") { return CrackNormalSource::PrincipalStrain; }
   throw std::invalid_argument("unknown crack normal source '" + s + "'");
}

} // namespace pff
} // namespace mfem

#endif
