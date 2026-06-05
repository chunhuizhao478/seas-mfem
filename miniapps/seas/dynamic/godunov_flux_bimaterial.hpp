// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// godunov_flux_bimaterial.hpp — Phase R.1 of
// PLAN_phase_R_exact_bimaterial_riemann_rev3.md.
//
// Exact linearised bi-material Riemann solver in the SeisSol
// formulation (Pelties et al. 2012, §2.3; SeisSol
// src/Equations/elastic/Model/ElasticSetup.h::getTransposedGodunovState).
//
// Replaces the average-flux approximation `F = ½(F⁺ + F⁻)` at
// heterogeneous interior faces with a single interface Riemann state
// `Q*` built by characteristic projection.  Per-side flux contributions
// then come from applying each side's OWN normal Jacobian to that
// shared `Q*`.
//
// MATRIX CONVENTION (SeisSol verbatim)
// ------------------------------------
// matR follows SeisSol's getTransposedGodunovState sign convention
// EXACTLY (ElasticSetup.h:102-139, non-acoustic elastic branch).
// The LOCAL cols 0..2 have POSITIVE stress entries (e.g. matR(SXX,0)
// = +(lam_L+2*mu_L)).  This differs from MFEM's GodunovFlux R, which
// uses NEGATIVE stress entries for the same columns.  Both are valid
// eigenvectors of the elastic Jacobian (eigenvectors are defined up
// to a sign), but the SIGN CHOICE MATTERS for the bi-material
// projector because matR mixes LOCAL and NEIGHBOUR materials —
// flipping cols 0..2's signs CHANGES the subspace structure of
// `matR · chi · matR^{-1}` and produces a different bi-material
// formula.  Empirical R.1 testing showed that the MFEM sign
// convention reduces the bi-material formula to the AVERAGE-flux
// approximation (no transmission/reflection effect), whereas
// SeisSol's sign convention gives the exact linearised Riemann
// solver (analytic T_p and R_p reproduced).
//
// FORMULA (SeisSol verbatim)
// --------------------------
// With SeisSol's matR sign convention,
//
//   chi          = diag(1, 1, 1, 0, 0, 0, 0, 0, 0)
//   godunov      = matR · chi · matR^{-1}
//   qGodLocal    = I - godunov     (multiplier of Q_L)
//   qGodNeighbor =     godunov     (multiplier of Q_R)
//   Q* = qGodLocal · Q_L + qGodNeighbor · Q_R
//
// Why this gives standard upwind in the homogeneous limit: SeisSol's
// matR cols 0..2 are eigenvectors of MFEM's A_x with eigenvalue
// −cp_L / −cs_L (not +cp_L) under MFEM's `dQ/dt + A · ∂Q/∂x = 0`
// convention.  The chi-projection therefore selects MFEM-NEGATIVE-
// eigenvalue modes, giving A · godunov = A^-.  Hence
// A · (I-godunov) = A - A^- = A^+, and
//   F = A · Q* = A · ((I-godunov) Q_L + godunov Q_R) = A^+ Q_L + A^- Q_R
// — the standard upwind flux that MFEM's `GodunovFlux::Interior`
// computes.  The "labelling vs. physics" inversion is a notational
// artefact of SeisSol's eigenvector sign choice; the resulting math
// is the standard upwind flux + correct bi-material physics.
//
// R.1.T-1 byte-exact gate: in the homogeneous limit, the bi-material
// formula inverts matR (cols 0..2 LOCAL, 6..8 NEIGHBOUR, 3..5 passive)
// per-call, whereas `GodunovFlux` inverts its R (rightgoing in cols
// 0..2, leftgoing in cols 6..8, passive 3..5) once at construction
// time and caches Ax_plus_/Ax_minus_.  Both paths use
// mfem::DenseMatrixInverse but on different matrices, so the FP
// rounding sequences differ; agreement is expected at the
// 1e-10..1e-12 relative level rather than bit-exact.

#ifndef MFEM_SEAS_GODUNOV_FLUX_BIMATERIAL_HPP
#define MFEM_SEAS_GODUNOV_FLUX_BIMATERIAL_HPP

#include "mfem.hpp"

#include "godunov_flux.hpp"
#include "wave_state.hpp"

namespace mfem
{
namespace seas
{

/// Bi-material exact linearised Riemann solver, SeisSol formulation
/// (Pelties et al. 2012 §2.3).  See file header for sign conventions
/// and the deviation note vs the rev-3 plan text.
class BimaterialFlux
{
public:
   /// @brief Build the per-face Godunov projector matrices `qGodLocal`
   /// and `qGodNeighbor` in the FACE-LOCAL frame
   /// (n̂ → +x̂, t̂₁ → +ŷ, t̂₂ → +ẑ).
   ///
   /// The bi-material interface state from the LOCAL element's POV is
   ///
   ///     Q*_facelocal = qGodLocal · Q_L_facelocal
   ///                  + qGodNeighbor · Q_R_facelocal
   ///
   /// where `Q_L`, `Q_R` are the local- and neighbour-side states
   /// already rotated into the face-local frame.  Each side then
   /// applies its own `A^side · Q*` to compute its own face flux.
   ///
   /// `qGodLocal + qGodNeighbor == I` (verified by R.1.T-7).  In the
   /// homogeneous limit `(lam_L,mu_L,rho_L) == (lam_R,mu_R,rho_R)`,
   /// the bi-material formulation collapses to MFEM's existing
   /// `GodunovFlux::Interior(...)` byte-equivalently up to LU rounding
   /// (verified by R.1.T-1).
   ///
   /// Pre-condition: neither material is acoustic (`mu > epsilon`).
   /// Aborts otherwise — SAFS does not have acoustic regions; adding
   /// the acoustic branch from
   /// SeisSol/src/Equations/elastic/Model/ElasticSetup.h:92-139
   /// is a documented follow-up.
   ///
   /// @param[in]  lam_self, mu_self, rho_self  Local-side material.
   /// @param[in]  lam_nbr,  mu_nbr,  rho_nbr   Neighbour-side material.
   /// @param[out] qGodLocal                     9×9 multiplier of Q_L.
   /// @param[out] qGodNeighbor                  9×9 multiplier of Q_R.
   static void BuildGodunovStateFaceLocal(real_t lam_self, real_t mu_self,
                                          real_t rho_self,
                                          real_t lam_nbr, real_t mu_nbr,
                                          real_t rho_nbr,
                                          mfem::DenseMatrix& qGodLocal,
                                          mfem::DenseMatrix& qGodNeighbor);

   /// @brief Convenience wrapper: build qGodLocal / qGodNeighbor in
   /// the face-local frame AND apply the face rotation to produce the
   /// per-side flux matrices in the GLOBAL frame.
   ///
   /// On output:
   ///
   ///   fluxLocal    = T · A_self_facelocal · qGodLocal    · T^{-1}
   ///   fluxNeighbor = T · A_self_facelocal · qGodNeighbor · T^{-1}
   ///
   /// At runtime the local element's face flux contribution (in the
   /// global frame) is
   ///
   ///   F_h_self_global = fluxLocal · Q_self_global
   ///                   + fluxNeighbor · Q_nbr_global
   ///
   /// `flux_self` supplies the LOCAL-side material (so `A_self` is
   /// built from `flux_self.GetAx()`).  For the neighbour element's
   /// contribution, the caller must invoke this function AGAIN with
   /// self/nbr swapped — the resulting `fluxLocal`/`fluxNeighbor`
   /// then correspond to the neighbour cell's POV.
   ///
   /// Each matrix is sized 9×9; the function calls `SetSize` so the
   /// caller need not pre-size.
   static void BuildPerFaceFluxMatricesGlobal(
      const real_t* nor,
      const GodunovFlux& flux_self,
      const GodunovFlux& flux_nbr,
      mfem::DenseMatrix& fluxLocal,        ///< 9×9 output
      mfem::DenseMatrix& fluxNeighbor);    ///< 9×9 output

   /// @brief Build the per-face CENTRAL-flux matrices in the GLOBAL frame
   /// — the mixed-flux companion to the upwind `BuildPerFaceFluxMatrices-
   /// Global`, used on fault-adjacent faces.  On output:
   ///
   ///   centralSelf = ½ · T · (AxPlus_self + AxMinus_self) · T^{-1}  (= ½·A_self)
   ///   centralNbr  = ½ · T · (AxPlus_nbr  + AxMinus_nbr ) · T^{-1}  (= ½·A_nbr)
   ///
   /// where `A_side` is that side's face-normal Jacobian in the GLOBAL
   /// frame.  At runtime the single-valued central flux is computed via
   /// the SAME `ApplyPerFaceFlux` used by the upwind path:
   ///
   ///   F* = centralSelf · Q_self + centralNbr · Q_nbr
   ///      = ½ ( A_self · Q_self + A_nbr · Q_nbr ).
   ///
   /// SINGLE-VALUEDNESS / CONSERVATION (the central-flux pitfall): the
   /// central flux has NO shared interface state — it is the literal
   /// average of the two PHYSICAL fluxes and is IDENTICAL on both sides of
   /// the face.  Each state is multiplied by ITS OWN side's Jacobian
   /// (`½ A_self·Q_self + ½ A_nbr·Q_nbr`).  Do **NOT** use
   /// `½ A_self·(Q_self + Q_nbr)` (one Jacobian for both states — e.g.
   /// `GodunovFlux::Central` applied per side) under heterogeneity: that is
   /// a DIFFERENT, non-single-valued flux when `A_self ≠ A_nbr` and
   /// silently breaks conservation.  In the homogeneous limit
   /// `A_self == A_nbr == A` both forms reduce to `½ A·(Q_self + Q_nbr)`
   /// == `GodunovFlux::Central` (verified to LU rounding by
   /// tests/unit/test_bimaterial_central_flux.cpp Test 1.1).
   ///
   /// The face-local full Jacobian per side is `GetAxPlus() + GetAxMinus()`
   /// — the SAME construction `GodunovFlux::Central` uses
   /// (godunov_flux.cpp) — so the homogeneous limit reproduces `Central`
   /// to ≤1e-11 relative, rather than via the separately-built `GetAx()`.
   ///
   /// Pre-condition: neither material is acoustic (`mu > 1e-12`).  Aborts
   /// otherwise (SAFS has no acoustic regions — same contract as
   /// `BuildGodunovStateFaceLocal`; the acoustic branch is a documented
   /// follow-up).  `nor` must be a unit vector.  Each matrix is sized 9×9
   /// via `SetSize`; the caller need not pre-size.
   ///
   /// @param[in]  nor         Unit face normal (global frame).
   /// @param[in]  flux_self   LOCAL-side material.
   /// @param[in]  flux_nbr    NEIGHBOUR-side material.
   /// @param[out] centralSelf 9×9 output: ½·A_self (global frame).
   /// @param[out] centralNbr  9×9 output: ½·A_nbr  (global frame).
   static void BuildPerFaceCentralMatricesGlobal(
      const real_t* nor,
      const GodunovFlux& flux_self,
      const GodunovFlux& flux_nbr,
      mfem::DenseMatrix& centralSelf,      ///< 9×9 output: ½·A_self (global)
      mfem::DenseMatrix& centralNbr);      ///< 9×9 output: ½·A_nbr  (global)

   /// @brief Runtime per-QP apply:
   /// `F_h_self = fluxLocal · Q_self + fluxNeighbor · Q_nbr`.
   ///
   /// `F_h_self` is a caller-allocated buffer of size NUM_STATE.
   /// `Q_self` and `Q_nbr` are NUM_STATE-length input buffers.
   /// `fluxLocal` and `fluxNeighbor` must be 9×9.  Aborts otherwise.
   static void ApplyPerFaceFlux(const mfem::DenseMatrix& fluxLocal,
                                const mfem::DenseMatrix& fluxNeighbor,
                                const real_t* Q_self,
                                const real_t* Q_nbr,
                                real_t* F_h_self);
};

}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_GODUNOV_FLUX_BIMATERIAL_HPP
