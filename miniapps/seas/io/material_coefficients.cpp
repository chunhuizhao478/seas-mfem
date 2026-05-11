// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// material_coefficients.cpp — `Eval` definitions for the three
// sidecar-backed elastic Coefficient subclasses.  Implementations are
// out-of-line because each one constructs a small `Vector` on the
// stack for `T.Transform`, which is awkward (though not wrong) to
// inline.

#include "material_coefficients.hpp"

namespace mfem
{
namespace seas
{

real_t RhoFromSidecar::Eval(mfem::ElementTransformation& T,
                            const mfem::IntegrationPoint& ip)
{
   mfem::Vector x(3);
   T.Transform(ip, x);
   return rho_field_.Evaluate(x[0], x[1], x[2]);
}


real_t MuFromSidecar::Eval(mfem::ElementTransformation& T,
                           const mfem::IntegrationPoint& ip)
{
   mfem::Vector x(3);
   T.Transform(ip, x);
   const real_t rho = rho_field_.Evaluate(x[0], x[1], x[2]);
   const real_t vs  = vs_field_ .Evaluate(x[0], x[1], x[2]);
   return rho * vs * vs;
}


real_t LambdaFromSidecar::Eval(mfem::ElementTransformation& T,
                               const mfem::IntegrationPoint& ip)
{
   mfem::Vector x(3);
   T.Transform(ip, x);
   const real_t rho = rho_field_.Evaluate(x[0], x[1], x[2]);
   const real_t vp  = vp_field_ .Evaluate(x[0], x[1], x[2]);
   const real_t vs  = vs_field_ .Evaluate(x[0], x[1], x[2]);
   // λ = ρ (Vp² − 2 Vs²).  No clamp — see header docstring.
   return rho * (vp * vp - 2.0 * vs * vs);
}

} // namespace seas
} // namespace mfem
