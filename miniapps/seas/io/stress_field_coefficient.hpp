// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// stress_field_coefficient.hpp — six-component VectorCoefficient
// backed by a StressField3D (Phase 6 §2 of PLAN_onfaultstress.md).
//
// Returns the six independent components of the symmetric Cauchy
// stress tensor in schema-v1 canonical order:
//
//   v[0] = sigma_xx
//   v[1] = sigma_yy
//   v[2] = sigma_zz
//   v[3] = sigma_xy
//   v[4] = sigma_yz
//   v[5] = sigma_xz
//
// IMPORTANT: this is NOT Voigt notation (which is
// (xx, yy, zz, yz, xz, xy)).  The schema-v1 order is the contract
// shared with the static-equilibration plan §1 (sigma_csm
// VectorCoefficient) and the Phase 4 fault-VTU emission.
//
// Sign convention: pure pass-through, compression POSITIVE SEAS
// (R-501 / R-502).  No sign manipulation here — the source-site
// flip happens once in Phase 3's `bulk_stress_tensor_field`.

#ifndef MFEM_SEAS_STRESS_FIELD_COEFFICIENT_HPP
#define MFEM_SEAS_STRESS_FIELD_COEFFICIENT_HPP

#include "mfem.hpp"

#include "stress_field_3d.hpp"

namespace mfem
{
namespace seas
{

/// 6-component VectorCoefficient that samples a StressField3D at
/// the physical coordinate of the evaluation point.
class StressFieldCoefficient : public mfem::VectorCoefficient
{
public:
   /// Returned value = scale * field(x, y, z) + offset, applied
   /// uniformly to all six components.  Default scale=1, offset=0
   /// (raw pass-through).
   StressFieldCoefficient(const StressField3D& field,
                          real_t scale = 1.0,
                          real_t offset = 0.0)
      : VectorCoefficient(StressField3D::NumComponents())
      , field_(field)
      , scale_(scale)
      , offset_(offset)
   { }

   void Eval(mfem::Vector& v, mfem::ElementTransformation& T,
             const mfem::IntegrationPoint& ip) override
   {
      v.SetSize(StressField3D::NumComponents());
      mfem::Vector x(3);
      T.Transform(ip, x);
      const real_t xc = x[0];
      const real_t yc = x[1];
      const real_t zc = x[2];

      // Schema-v1 canonical order — NOT Voigt.
      v(0) = scale_ * field_.Field(0).Evaluate(xc, yc, zc) + offset_;  // xx
      v(1) = scale_ * field_.Field(1).Evaluate(xc, yc, zc) + offset_;  // yy
      v(2) = scale_ * field_.Field(2).Evaluate(xc, yc, zc) + offset_;  // zz
      v(3) = scale_ * field_.Field(3).Evaluate(xc, yc, zc) + offset_;  // xy
      v(4) = scale_ * field_.Field(4).Evaluate(xc, yc, zc) + offset_;  // yz
      v(5) = scale_ * field_.Field(5).Evaluate(xc, yc, zc) + offset_;  // xz
   }

   const StressField3D& Field() const { return field_; }
   real_t Scale()  const { return scale_; }
   real_t Offset() const { return offset_; }

private:
   const StressField3D& field_;
   real_t scale_;
   real_t offset_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_STRESS_FIELD_COEFFICIENT_HPP
