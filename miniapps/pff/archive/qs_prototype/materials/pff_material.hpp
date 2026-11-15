// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_PFF_MATERIAL_HPP
#define MFEM_PFF_MATERIAL_HPP

#include "mfem.hpp"
#include "degradation_function.hpp"
#include <cmath>
#include <iostream>

namespace mfem
{

namespace pff
{

/** @brief Material parameters container for phase field fracture.
 *
 * This struct holds all material parameters needed for the phase field
 * fracture simulation, including:
 *   - Elastic properties (E, nu, lambda, mu, K)
 *   - Fracture properties (Gc, l, c0)
 *   - Degradation parameters (eta, p)
 *
 * It also provides utility functions to compute derived quantities.
 */
struct PFFMaterialParameters
{
   // =========================================================================
   // Elastic properties
   // =========================================================================

   /// Young's modulus [Pa]
   real_t E;

   /// Poisson's ratio [-]
   real_t nu;

   /// First Lame parameter: lambda = E*nu / ((1+nu)*(1-2*nu)) [Pa]
   real_t lambda;

   /// Second Lame parameter (shear modulus): mu = E / (2*(1+nu)) [Pa]
   real_t mu;

   /// Bulk modulus: K = E / (3*(1-2*nu)) [Pa]
   real_t K;

   // =========================================================================
   // Fracture properties
   // =========================================================================

   /// Fracture toughness (critical energy release rate) [J/m^2]
   real_t Gc;

   /// Regularization length (phase field width) [m]
   real_t l;

   /// Normalization constant (c0 = 2 for AT2 model)
   real_t c0;

   /// Tensile strength (used to derive Gc if not given directly) [Pa]
   real_t sigma_t;

   // =========================================================================
   // Degradation parameters
   // =========================================================================

   /// Residual stiffness (small value to prevent singularity)
   real_t eta;

   /// Degradation exponent (p = 2 for AT2 model)
   real_t p;

   // =========================================================================
   // Constructors
   // =========================================================================

   /** @brief Default constructor with typical values.
    */
   PFFMaterialParameters()
      : E(40e6), nu(0.25), lambda(0.0), mu(0.0), K(0.0),
        Gc(0.0), l(5e-5), c0(2.0), sigma_t(6.43e6),
        eta(1e-6), p(2.0)
   {
      ComputeLameParameters();
      ComputeFractureToughness(sigma_t);
   }

   /** @brief Constructor with basic elastic and fracture parameters.
    *
    * Computes derived quantities (lambda, mu, K, Gc) automatically.
    *
    * @param[in] E_ Young's modulus [Pa]
    * @param[in] nu_ Poisson's ratio [-]
    * @param[in] l_ Regularization length [m]
    * @param[in] sigma_t_ Tensile strength [Pa] (used to derive Gc)
    * @param[in] eta_ Residual stiffness (default: 1e-6)
    * @param[in] p_ Degradation exponent (default: 2.0)
    */
   PFFMaterialParameters(real_t E_, real_t nu_, real_t l_,
                         real_t sigma_t_, real_t eta_ = 1e-6, real_t p_ = 2.0)
      : E(E_), nu(nu_), lambda(0.0), mu(0.0), K(0.0),
        Gc(0.0), l(l_), c0(2.0), sigma_t(sigma_t_),
        eta(eta_), p(p_)
   {
      ComputeLameParameters();
      ComputeFractureToughness(sigma_t);
   }

   /** @brief Constructor with explicit Gc (no strength-based derivation).
    *
    * @param[in] E_ Young's modulus [Pa]
    * @param[in] nu_ Poisson's ratio [-]
    * @param[in] l_ Regularization length [m]
    * @param[in] Gc_ Fracture toughness [J/m^2]
    * @param[in] eta_ Residual stiffness (default: 1e-6)
    * @param[in] p_ Degradation exponent (default: 2.0)
    * @param[in] use_explicit_Gc Dummy parameter to distinguish constructors
    */
   PFFMaterialParameters(real_t E_, real_t nu_, real_t l_,
                         real_t Gc_, real_t eta_, real_t p_, bool use_explicit_Gc)
      : E(E_), nu(nu_), lambda(0.0), mu(0.0), K(0.0),
        Gc(Gc_), l(l_), c0(2.0), sigma_t(0.0),
        eta(eta_), p(p_)
   {
      (void)use_explicit_Gc; // Suppress unused parameter warning
      ComputeLameParameters();
      // Gc is already set, compute sigma_t for reference
      if (Gc > 0.0 && l > 0.0 && E > 0.0)
      {
         sigma_t = std::sqrt(3.0 * E * Gc / (8.0 * l));
      }
   }

   // =========================================================================
   // Utility methods
   // =========================================================================

   /** @brief Compute Lame parameters from Young's modulus and Poisson's ratio.
    *
    * lambda = E * nu / ((1 + nu) * (1 - 2*nu))
    * mu = E / (2 * (1 + nu))
    * K = E / (3 * (1 - 2*nu))
    */
   void ComputeLameParameters()
   {
      MFEM_VERIFY(E > 0.0, "Young's modulus must be positive");
      MFEM_VERIFY(nu > -1.0 && nu < 0.5,
                  "Poisson's ratio must be in (-1, 0.5)");

      lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
      mu = E / (2.0 * (1.0 + nu));
      K = E / (3.0 * (1.0 - 2.0 * nu));
   }

   /** @brief Compute fracture toughness from tensile strength.
    *
    * For AT2 model: Gc = 8 * l * sigma_t^2 / (3 * E)
    *
    * This relation comes from the analytical solution for a 1D bar
    * under tension, ensuring that the peak stress equals sigma_t.
    *
    * @param[in] sigma_t_ Tensile strength [Pa]
    */
   void ComputeFractureToughness(real_t sigma_t_)
   {
      sigma_t = sigma_t_;
      MFEM_VERIFY(sigma_t > 0.0, "Tensile strength must be positive");
      MFEM_VERIFY(l > 0.0, "Regularization length must be positive");
      MFEM_VERIFY(E > 0.0, "Young's modulus must be positive");

      Gc = 8.0 * l * sigma_t * sigma_t / (3.0 * E);
   }

   /** @brief Set fracture toughness directly.
    *
    * @param[in] Gc_ Fracture toughness [J/m^2]
    */
   void SetFractureToughness(real_t Gc_)
   {
      Gc = Gc_;
      MFEM_VERIFY(Gc > 0.0, "Fracture toughness must be positive");

      // Update sigma_t for consistency
      if (l > 0.0 && E > 0.0)
      {
         sigma_t = std::sqrt(3.0 * E * Gc / (8.0 * l));
      }
   }

   /** @brief Create a DegradationFunction object from the parameters.
    *
    * @return DegradationFunction with eta and p from this struct
    */
   DegradationFunction GetDegradationFunction() const
   {
      return DegradationFunction(eta, p);
   }

   /** @brief Get the coefficient for the diffusion term in damage equation.
    *
    * Coefficient: Gc * l / c0
    *
    * @return Diffusion coefficient [J/m]
    */
   real_t GetDamageDiffusionCoeff() const
   {
      return Gc * l / c0;
   }

   /** @brief Get the coefficient for the reaction term in damage equation.
    *
    * Coefficient: 2 * Gc / (c0 * l)
    *
    * @return Reaction coefficient [J/m^3]
    */
   real_t GetDamageReactionCoeff() const
   {
      return 2.0 * Gc / (c0 * l);
   }

   /** @brief Get the critical strain energy density for damage initiation.
    *
    * psi_c = Gc / (c0 * l)
    *
    * @return Critical strain energy density [J/m^3]
    */
   real_t GetCriticalStrainEnergy() const
   {
      return Gc / (c0 * l);
   }

   /** @brief Print material parameters to output stream.
    *
    * @param[out] os Output stream
    */
   void Print(std::ostream &os = std::cout) const
   {
      os << "Phase Field Fracture Material Parameters:\n"
         << "=========================================\n"
         << "Elastic Properties:\n"
         << "  E (Young's modulus)      = " << E << " Pa\n"
         << "  nu (Poisson's ratio)     = " << nu << "\n"
         << "  lambda (1st Lame param)  = " << lambda << " Pa\n"
         << "  mu (shear modulus)       = " << mu << " Pa\n"
         << "  K (bulk modulus)         = " << K << " Pa\n"
         << "\nFracture Properties:\n"
         << "  Gc (fracture toughness)  = " << Gc << " J/m^2\n"
         << "  l (regularization length)= " << l << " m\n"
         << "  c0 (normalization const) = " << c0 << "\n"
         << "  sigma_t (tensile strength)= " << sigma_t << " Pa\n"
         << "\nDegradation Parameters:\n"
         << "  eta (residual stiffness) = " << eta << "\n"
         << "  p (degradation exponent) = " << p << "\n"
         << "\nDerived Quantities:\n"
         << "  Diffusion coeff (Gc*l/c0)= " << GetDamageDiffusionCoeff() << " J/m\n"
         << "  Reaction coeff (2Gc/c0/l)= " << GetDamageReactionCoeff() << " J/m^3\n"
         << "  Critical psi (Gc/c0/l)   = " << GetCriticalStrainEnergy() << " J/m^3\n"
         << std::endl;
   }

   /** @brief Validate that all parameters are physically meaningful.
    *
    * @return true if all parameters are valid
    */
   bool Validate() const
   {
      bool valid = true;

      if (E <= 0.0)
      {
         std::cerr << "Error: Young's modulus must be positive.\n";
         valid = false;
      }
      if (nu <= -1.0 || nu >= 0.5)
      {
         std::cerr << "Error: Poisson's ratio must be in (-1, 0.5).\n";
         valid = false;
      }
      if (Gc <= 0.0)
      {
         std::cerr << "Error: Fracture toughness must be positive.\n";
         valid = false;
      }
      if (l <= 0.0)
      {
         std::cerr << "Error: Regularization length must be positive.\n";
         valid = false;
      }
      if (eta <= 0.0 || eta >= 1.0)
      {
         std::cerr << "Warning: Residual stiffness eta should be in (0, 1).\n";
         // Not a fatal error, just a warning
      }
      if (p <= 0.0)
      {
         std::cerr << "Error: Degradation exponent must be positive.\n";
         valid = false;
      }

      return valid;
   }
};

} // namespace pff

} // namespace mfem

#endif // MFEM_PFF_MATERIAL_HPP
