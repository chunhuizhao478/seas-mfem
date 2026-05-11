// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// material_coefficients.hpp — three `mfem::Coefficient` subclasses that
// derive elastic moduli (λ, μ) and density (ρ) from CVM-H sidecar
// fields (Vp, Vs, density).  These are the only place in the SEAS code
// where the conversion
//
//     ρ(x,y,z) = rho_field.Evaluate(x,y,z)
//     μ(x,y,z) = ρ · Vs(x,y,z)^2
//     λ(x,y,z) = ρ · Vp(x,y,z)^2 − 2 μ
//
// lives.  The rest of the operator code only ever touches λ, μ, ρ.
//
// Each Coefficient holds non-owning references to the `DataField3D`
// objects it needs.  The driver owns the `DataField3D`s and must
// outlive these Coefficients.
//
// Implementing these as derived `Coefficient` classes (rather than
// `FunctionCoefficient` lambdas) keeps the `Eval` signature stable
// and avoids the `std::function` indirection that `FunctionCoefficient`
// adds at every quadrature point.

#ifndef MFEM_SEAS_MATERIAL_COEFFICIENTS_HPP
#define MFEM_SEAS_MATERIAL_COEFFICIENTS_HPP

#include "mfem.hpp"

#include "data_field_3d.hpp"

namespace mfem
{
namespace seas
{

/// Density Coefficient: returns sidecar `density` field trilinearly
/// interpolated to the supplied (T, ip) physical location.
class RhoFromSidecar : public mfem::Coefficient
{
public:
   /// Caller retains ownership of `rho_field`; it must outlive this
   /// Coefficient.
   explicit RhoFromSidecar(const DataField3D& rho_field)
      : rho_field_(rho_field)
   { }

   real_t Eval(mfem::ElementTransformation& T,
               const mfem::IntegrationPoint& ip) override;

private:
   const DataField3D& rho_field_;
};


/// Shear modulus Coefficient: μ = ρ · Vs^2.
/// Both `rho_field` and `vs_field` must outlive this Coefficient.
class MuFromSidecar : public mfem::Coefficient
{
public:
   MuFromSidecar(const DataField3D& vs_field,
                 const DataField3D& rho_field)
      : vs_field_(vs_field), rho_field_(rho_field)
   { }

   real_t Eval(mfem::ElementTransformation& T,
               const mfem::IntegrationPoint& ip) override;

private:
   const DataField3D& vs_field_;
   const DataField3D& rho_field_;
};


/// First Lamé parameter Coefficient: λ = ρ · Vp^2 − 2 μ
///                                      = ρ · (Vp^2 − 2 Vs^2).
/// All three fields must outlive this Coefficient.  The Coefficient
/// does NOT clamp negative λ values; if `Vp^2 < 2 Vs^2` at any sidecar
/// voxel, λ < 0 propagates through — the sidecar builder
/// (`build_velocity_cvmh.py`) is responsible for the per-field range
/// guard.
class LambdaFromSidecar : public mfem::Coefficient
{
public:
   LambdaFromSidecar(const DataField3D& vp_field,
                     const DataField3D& vs_field,
                     const DataField3D& rho_field)
      : vp_field_(vp_field),
        vs_field_(vs_field),
        rho_field_(rho_field)
   { }

   real_t Eval(mfem::ElementTransformation& T,
               const mfem::IntegrationPoint& ip) override;

private:
   const DataField3D& vp_field_;
   const DataField3D& vs_field_;
   const DataField3D& rho_field_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_MATERIAL_COEFFICIENTS_HPP
