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

#ifndef MFEM_PFF_DEGRADATION_FUNCTION_HPP
#define MFEM_PFF_DEGRADATION_FUNCTION_HPP

#include "mfem.hpp"
#include <cmath>

namespace mfem
{

namespace pff
{

/** @brief Degradation function for phase field fracture.
 *
 * Implements the AT2 (Ambrosio-Tortorelli type 2) degradation function:
 *   g(d) = (1-d)^p * (1-eta) + eta
 *
 * where:
 *   - d is the damage variable (0 = intact, 1 = fully damaged)
 *   - p is the exponent (typically 2)
 *   - eta is the residual stiffness (small positive value to avoid singularity)
 *
 * Properties:
 *   - g(0) = 1 (intact material has full stiffness)
 *   - g(1) = eta (fully damaged material has residual stiffness)
 *   - g'(d) < 0 for d in (0,1) (monotonically decreasing)
 */
class DegradationFunction
{
public:
   /** @brief Construct degradation function with given parameters.
    *
    * @param[in] eta Residual stiffness (default: 1e-6)
    * @param[in] p Exponent (default: 2.0 for AT2 model)
    */
   DegradationFunction(real_t eta = 1e-6, real_t p = 2.0)
      : eta_(eta), p_(p), one_minus_eta_(1.0 - eta) {}

   /** @brief Evaluate the degradation function g(d).
    *
    * g(d) = (1-d)^p * (1-eta) + eta
    *
    * @param[in] d Damage variable in [0, 1]
    * @return Degradation factor in [eta, 1]
    */
   MFEM_HOST_DEVICE real_t Eval(real_t d) const
   {
      const real_t one_minus_d = 1.0 - d;
      return std::pow(one_minus_d, p_) * one_minus_eta_ + eta_;
   }

   /** @brief Evaluate the first derivative g'(d).
    *
    * g'(d) = -p * (1-d)^(p-1) * (1-eta)
    *
    * @param[in] d Damage variable in [0, 1]
    * @return First derivative of degradation function
    */
   MFEM_HOST_DEVICE real_t EvalDerivative(real_t d) const
   {
      const real_t one_minus_d = 1.0 - d;
      return -p_ * std::pow(one_minus_d, p_ - 1.0) * one_minus_eta_;
   }

   /** @brief Evaluate the second derivative g''(d).
    *
    * g''(d) = p * (p-1) * (1-d)^(p-2) * (1-eta)
    *
    * @param[in] d Damage variable in [0, 1]
    * @return Second derivative of degradation function
    */
   MFEM_HOST_DEVICE real_t EvalSecondDerivative(real_t d) const
   {
      const real_t one_minus_d = 1.0 - d;
      return p_ * (p_ - 1.0) * std::pow(one_minus_d, p_ - 2.0) * one_minus_eta_;
   }

   /// Get the residual stiffness parameter
   real_t GetEta() const { return eta_; }

   /// Get the exponent parameter
   real_t GetP() const { return p_; }

private:
   real_t eta_;           ///< Residual stiffness
   real_t p_;             ///< Exponent
   real_t one_minus_eta_; ///< Cached 1 - eta
};

/** @brief Crack geometric function for phase field fracture.
 *
 * Implements the crack geometric function for AT2 model:
 *   alpha(d) = d^2
 *
 * The normalization constant for AT2 is c0 = 2.
 */
class CrackGeometricFunction
{
public:
   /** @brief Evaluate the crack geometric function alpha(d).
    *
    * alpha(d) = d^2 (for AT2 model)
    *
    * @param[in] d Damage variable in [0, 1]
    * @return Crack geometric function value
    */
   MFEM_HOST_DEVICE static real_t Eval(real_t d)
   {
      return d * d;
   }

   /** @brief Evaluate the first derivative alpha'(d).
    *
    * alpha'(d) = 2*d (for AT2 model)
    *
    * @param[in] d Damage variable in [0, 1]
    * @return First derivative
    */
   MFEM_HOST_DEVICE static real_t EvalDerivative(real_t d)
   {
      return 2.0 * d;
   }

   /** @brief Evaluate the second derivative alpha''(d).
    *
    * alpha''(d) = 2 (for AT2 model)
    *
    * @param[in] d Damage variable in [0, 1]
    * @return Second derivative
    */
   MFEM_HOST_DEVICE static real_t EvalSecondDerivative(real_t /* d */)
   {
      return 2.0;
   }

   /** @brief Get the normalization constant c0.
    *
    * For AT2 model: c0 = 2
    *
    * @return Normalization constant
    */
   static real_t GetNormalizationConstant()
   {
      return 2.0;
   }
};

} // namespace pff

} // namespace mfem

#endif // MFEM_PFF_DEGRADATION_FUNCTION_HPP
