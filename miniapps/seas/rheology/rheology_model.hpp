// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// rheology_model.hpp — the pointwise interface every off-fault rheology of
// the Nitsche solver implements (rheology plan rev 1, sections 4 and 6.1).
//
// The solver keeps the elastic stress it has today and adds one correction
// Delta sigma from the rheology (plan eq. 1.1):
//
//    sigma = sigma_0'(x) + C_0 : eps(u) + Delta sigma(eps, q, dt)
//
// sigma_0' is the effective initial stress (initial stress plus fluid
// pressure, plan eq. 3.3), C_0 the isotropic elastic tensor of the mesh, q a
// per-node internal state.  A model sees one rheology node at a time and
// knows nothing about meshes, so each model is unit-tested alone
// (tests/unit/test_rheology_pointwise.cpp).
//
// Symmetric tensors are 6-vectors in the order of the SCEC sample code
// (TPV26/27 description, p. 16): [xx, yy, zz, xy, yz, xz], global code
// frame, with TENSOR shear components (eps_xy, not the engineering
// 2 eps_xy used by constitutive/linear_elastic.hpp).  Stresses are tension
// positive, in Pa.
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_RHEOLOGY_MODEL_HPP
#define MFEM_SEAS_RHEOLOGY_MODEL_HPP

#include "mfem.hpp"

namespace mfem
{
namespace seas
{

/// Symmetric-tensor helpers for the [xx, yy, zz, xy, yz, xz] layout.
namespace sym6
{

constexpr int kSize = 6;

inline real_t Trace(const real_t a[6]) { return a[0] + a[1] + a[2]; }

/// a : b
inline real_t DoubleDot(const real_t a[6], const real_t b[6])
{
   return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
          + 2.0 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5]);
}

/// sig = C_0 : eps = lambda tr(eps) I + 2 mu eps
inline void ApplyIsotropic(real_t lambda, real_t mu, const real_t eps[6], real_t sig[6])
{
   const real_t le = lambda * Trace(eps);
   sig[0] = le + 2.0 * mu * eps[0];
   sig[1] = le + 2.0 * mu * eps[1];
   sig[2] = le + 2.0 * mu * eps[2];
   sig[3] = 2.0 * mu * eps[3];
   sig[4] = 2.0 * mu * eps[4];
   sig[5] = 2.0 * mu * eps[5];
}

/// Row-major symmetric 3x3 -> 6-vector (off-diagonals averaged).
inline void FromRowMajor(const real_t m[9], real_t a[6])
{
   a[0] = m[0];
   a[1] = m[4];
   a[2] = m[8];
   a[3] = 0.5 * (m[1] + m[3]);
   a[4] = 0.5 * (m[5] + m[7]);
   a[5] = 0.5 * (m[2] + m[6]);
}

/// 6-vector -> row-major symmetric 3x3.
inline void ToRowMajor(const real_t a[6], real_t m[9])
{
   m[0] = a[0]; m[1] = a[3]; m[2] = a[5];
   m[3] = a[3]; m[4] = a[1]; m[5] = a[4];
   m[6] = a[5]; m[7] = a[4]; m[8] = a[2];
}

} // namespace sym6

/// What a model may know about one rheology node.
struct RheologyPointContext
{
   real_t x[3] = {0.0, 0.0, 0.0};   ///< node position (m)
   real_t sigma0_eff[6] = {};       ///< sigma_0'(x) (Pa), plan eq. 3.3
   real_t lambda0 = 0.0;            ///< C_0 Lame parameters at the node (Pa)
   real_t mu0 = 0.0;
   real_t t = 0.0;                  ///< time at the end of the step (s)
};

/// One pointwise rheology.  Models are combined by RheologyStack in a fixed
/// order; each one receives the correction of the models before it.
class IRheologyModel
{
public:
   virtual ~IRheologyModel() = default;

   virtual const char *Name() const = 0;

   /// State variables per rheology node.
   virtual int NumStateVars() const = 0;

   /// True when Delta sigma is linear in the strain history (attenuation);
   /// false for plasticity and damage.  Linear-only shortcuts read this.
   virtual bool IsLinearInStrain() const = 0;

   /// State at t = 0.
   virtual void InitState(const RheologyPointContext &ctx, real_t *q) const = 0;

   /// One explicit step to t_{n+1}.
   ///   in:  eps  = total strain at t_{n+1}
   ///        dt   = t_{n+1} - t_n
   ///        q    = state at t_n
   ///        dsig = correction of the models before this one in the stack
   ///               (zero for the first)
   ///   out: q    = state at t_{n+1}
   ///        dsig = correction after this model, so that
   ///               sigma^{n+1} = sigma_0' + C_0 : eps + dsig
   /// Returns true when the step was inelastic at this node (the state
   /// changed), for the per-step yield count.
   virtual bool Update(const RheologyPointContext &ctx, const real_t eps[6], real_t dt,
                       real_t *q, real_t dsig[6]) const = 0;

   /// Largest ratio of the model's tangent stiffness to C_0; the explicit
   /// CFL limit holds only when this is <= 1.
   virtual real_t StiffnessBound(const RheologyPointContext &ctx, const real_t *q) const = 0;

   /// Output fields derived from the state (e.g. plastic strain, eta).
   virtual int NumDiagnostics() const = 0;
   virtual const char *DiagnosticName(int i) const = 0;
   virtual void Diagnostics(const RheologyPointContext &ctx, const real_t *q,
                            real_t *out) const = 0;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RHEOLOGY_MODEL_HPP
