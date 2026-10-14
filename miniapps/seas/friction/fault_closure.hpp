// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_closure.hpp — Nitsche dynamic rupture, plan item N-03.
//
// Pointwise friction closure for the Nitsche fault: given the trial traction
// X_tan (the shear traction the bulk would carry with the current slip held
// fixed) and the coupling impedance gamma, return the slip rate V_tan that
// satisfies the projection identity
//
//    tau = [ X_tan - gamma V_tan ]_{tau_s},   V_tan anti-parallel to X_tan,
//
// where [.]_{tau_s} is the projection onto the friction disc |tau| <= tau_s.
// This is the projection form, not a scalar root-find, so a law with a stick
// branch (slip weakening, P9) fits behind the same interface.  gamma is an
// argument on every call; it is never stored as a damping coefficient.
//
// State evolution is exposed as a stateless rate; the time integration of
// the state belongs to the Nitsche operator (P6).
//
// Scope (Chunhui Zhao, 2026-10-14).

#ifndef MFEM_SEAS_FAULT_CLOSURE_HPP
#define MFEM_SEAS_FAULT_CLOSURE_HPP

#include "mfem.hpp"

namespace mfem
{
namespace seas
{

/// State and local parameters of one fault point.
struct FaultPointState
{
   real_t psi = 0.0;   ///< logarithmic state variable [-]
   real_t a = 0.0;     ///< direct-effect parameter at this point [-]
   real_t Dc = 0.0;    ///< critical slip distance at this point [m]
};

class IFaultClosure
{
public:
   virtual ~IFaultClosure() = default;

   /// Projection identity at one point; V_tan anti-parallel to X_tan.
   /// Returns |V|.
   virtual real_t SolveSlipRate(const real_t X_tan[2], real_t gamma,
                                real_t sigma_n, const FaultPointState &s,
                                real_t V_tan[2]) const = 0;

   /// d(state)/dt at (V, s).  Stateless; time integration belongs to the
   /// Nitsche operator (P6).
   virtual real_t StateRate(real_t V, const FaultPointState &s) const = 0;

   /// Seed the state from equilibrium at t = 0 (the continuous condition:
   /// no gamma, no radiation damping).
   virtual void SeedEquilibrium(real_t tau0, real_t V_ini, real_t sigma_n,
                                FaultPointState &s) const = 0;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_CLOSURE_HPP
