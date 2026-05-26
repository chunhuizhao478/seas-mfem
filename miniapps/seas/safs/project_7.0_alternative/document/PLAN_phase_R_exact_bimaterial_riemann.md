# Implementation Plan: Phase R — Exact Bi-material Riemann Solver

**Status:** plan (revised against SeisSol reference)
**Date:** 2026-05-18
**Companion to:** `spatial_dynamic_rupture_plan.md` (Phase H) and
`heterogeneous_material_plan.md` (Phase 3).  This plan **replaces** the
"Out of Scope item 1" (exact bi-material Riemann solver) in
`spatial_dynamic_rupture_plan.md:586` — that line is no longer a
future-work caveat; it is now this phase.

**Reference implementation:** SeisSol
(`/Users/chunhuizhao/projects/SeisSol`).  This plan was rev-2'd after
cross-checking the math, eigenvector packing, and dispatch architecture
against `src/Equations/elastic/Model/ElasticSetup.h::getTransposedGodunovState`,
`src/Initializer/CellLocalMatrices.cpp::initializeCellMatrices`, and
the auto-generated kernels in `codegen/kernels/aderdg.py::computeFluxSolverLocal`.
Cited file paths use the SeisSol commit checked out under that prefix.

## Overview

Replace the **average-flux approximation** at heterogeneous interior
faces (`F = ½(F⁺ + F⁻)`, each side's flux from its own
`GodunovFlux`) with the **exact linearised Riemann solution** in the
SeisSol formulation: build a single interface state `Q*` by projecting
the jump `(Q_R − Q_L)` onto the local-going characteristic family,
then apply each side's OWN local Jacobian to `Q*` to produce each
cell's per-face flux contribution.  Per-face Godunov state matrices
(`qGodLocal`, `qGodNeighbor`) are **precomputed at WaveOperator
construction time** and applied per-QP at runtime, mirroring SeisSol's
`getTransposedGodunovState` + `computeFluxSolverLocal` pattern.

The change is mandatory for SAFS production because the CVM-H / CVM-S
sidecars produce impedance contrasts `Z⁻/Z⁺ ≈ 2–4` between basement
and basin elements; on those contrasts the average-flux approximation
(Pelties et al. 2014, §5; Pelties et al. 2012, §2.3) produces spurious
reflections that contaminate slip-rate on the embedded fault by
10–30 % at the first back-radiation arrival.  SeisSol uses the exact
solver in production for the same physical problem; matching their
formulation gives us byte-comparable verification.

## Commit cadence — 4 staged commits with a clean bisect surface

Phase R lands in **four separate commits**, each gating a distinct
class of failure.  Doing them separately is the difference between
"if SAFS shows wrong physics we know exactly which commit to bisect"
and "we get a mystery answer at the end and have to forensic-debug
across the whole rollout".

| Commit | Phases | What it catches if it FAILS | Cost |
|--------|--------|-----------------------------|------|
| **#1** | R.1 + R.2 | Math errors in the standalone solver (index ordering, transposition, eigenvector packing, rounding helper drift).  9 unit tests including `T-PHASER-HOMOG-PARITY` (byte-exact vs `GodunovFlux::Interior` on homogeneous input — the strictest math check the codebase has). | minutes |
| **#2** | R.3 (dispatch wiring) + Phase H Stage 2 | Drift in the SCALAR / TPV / BP5 code paths from inserting the new `if (owned_flux_pool_ && IsBimaterialFace(...))` branch in `wave_operator.inl`.  TPV102/104/205 + BP5 byte-exact regression suite IS THE HEADLINE acceptance gate (R.3.T-1).  Every face on those drivers hits the `IsHomogeneous == true` fall-through to the existing `Interior(...)`; any drift means the dispatch is wrong. | minutes |
| **#3** | R.4 (layered-medium analytic R/T) | The Riemann math is implemented correctly at full WaveOperator scale.  This is the ONLY place where the bi-material physics is verified against an analytic solution; unit tests prove the math at a single face, but only R.4 proves it integrates correctly through `Mult()` and produces the analytic transmission / reflection coefficients of LeVeque §22.4. | ~10 min |
| **#4** | R.5 (SAFS A/B verification + plan cleanup) | The production-scale impact is measurable AND the exact path converges to a refined-mesh reference faster than the average path.  Cheapest way to catch a SAFS-coupling regression that the lower-level tests miss. | hours |

**Strict ordering**: each commit requires the previous one to be
green.  In particular:

- **Commit #2 (R.3 + Phase H Stage 2) cannot land until Commit #1
  (R.1) is green** — the dispatch calls into `BimaterialFlux::
  ApplyPerFaceFlux`, which must already exist + be unit-tested.
- **Commit #2 (R.3 + Phase H Stage 2) is the FINAL gate before any
  heterogeneous-material physics is exposed to production**.  If
  TPV/BP5 byte-exact fails here, STOP and bisect against pre-R.3
  baseline (Commit b18bed5).  Do NOT advance to R.4 or R.5 with a
  red scalar-path regression.
- **Commit #3 (R.4) cannot land until Commit #2 is green** — the
  layered-medium test runs `WaveOperator::Mult()` and requires the
  dispatch to actually reach `ApplyPerFaceFlux` on heterogeneous
  faces.
- **Commit #4 (R.5) cannot land until Commit #3 is green** — without
  R.4's analytic-R/T verification we have no proof the math is right
  at scale; R.5's A/B against average-flux can show "different"
  but not "correct".

## Mathematical formulation (SeisSol-validated)

### Eigenvector packing — the canonical `matR`

In the 9-component velocity-stress system `Q = (σ_xx, σ_yy, σ_zz,
σ_xy, σ_yz, σ_xz, v_x, v_y, v_z)^T`, rotated into the face-local
frame (`n̂ → +x̂`, `t̂₁ → +ŷ`, `t̂₂ → +ẑ`), the normal Jacobian `A_x`
has 9 eigenvalues:

| index | eigenvalue | physical content                              | from which side? |
|------:|------------|-----------------------------------------------|------------------|
| 0     | `+c_p^L`   | local right-going P-wave                      | local            |
| 1     | `+c_s^L`   | local right-going S-wave, y-polarised         | local            |
| 2     | `+c_s^L`   | local right-going S-wave, z-polarised         | local            |
| 3     | 0          | σ_yy mode (advection-free)                    | (zero mode)      |
| 4     | 0          | σ_zz mode                                     | (zero mode)      |
| 5     | 0          | σ_yz mode                                     | (zero mode)      |
| 6     | `−c_s^R`   | neighbour left-going S-wave, z-polarised      | neighbour        |
| 7     | `−c_s^R`   | neighbour left-going S-wave, y-polarised      | neighbour        |
| 8     | `−c_p^R`   | neighbour left-going P-wave                   | neighbour        |

The columns of `matR` are the right eigenvectors in this ordering.
**Critically, the first three columns use LOCAL impedances and the
last three use NEIGHBOUR impedances** — this is what makes the
matrix a bi-material projector rather than a one-sided eigendecomposition.

**SeisSol reference** — `getTransposedGodunovState` at
`src/Equations/elastic/Model/ElasticSetup.h:80-166`.  Reading that
function line-by-line:

```cpp
// Cols 0,1,2: LOCAL outgoing eigenvectors
matR(0, 0) = local.lambda + 2 * local.mu;   // P-wave: σ_nn entry
matR(1, 0) = local.lambda;                  //          σ_t1t1
matR(2, 0) = local.lambda;                  //          σ_t2t2
matR(6, 0) = std::sqrt((local.lambda + 2 * local.mu) / local.rho);  // v_n
matR(3, 1) = local.mu;                      // S1-wave: σ_nt1
matR(7, 1) = std::sqrt(local.mu / local.rho);                       //          v_t1
matR(5, 2) = local.mu;                      // S2-wave: σ_nt2
matR(8, 2) = std::sqrt(local.mu / local.rho);                       //          v_t2

// Cols 3,4,5: zero-mode eigenvectors (scaled for matR conditioning)
matR(4, 3) = local.lambda + 2 * local.mu;   // σ_yz mode column
matR(1, 4) = local.lambda + 2 * local.mu;   // σ_yy mode column
matR(2, 5) = local.lambda + 2 * local.mu;   // σ_zz mode column

// Cols 6,7,8: NEIGHBOUR incoming eigenvectors (negative wave speeds)
matR(5, 6) = neighbor.mu;                              // S2-wave: σ_nt2
matR(8, 6) = -std::sqrt(neighbor.mu / neighbor.rho);   //          v_t2
matR(3, 7) = neighbor.mu;                              // S1-wave: σ_nt1
matR(7, 7) = -std::sqrt(neighbor.mu / neighbor.rho);   //          v_t1
matR(0, 8) = neighbor.lambda + 2 * neighbor.mu;        // P-wave: σ_nn
matR(1, 8) = neighbor.lambda;                          //         σ_t1t1
matR(2, 8) = neighbor.lambda;                          //         σ_t2t2
matR(6, 8) = -std::sqrt((neighbor.lambda + 2 * neighbor.mu) / neighbor.rho);  // v_n
```

(SeisSol's permutation of zero-mode columns differs slightly from the
LeVeque diagonal convention — they pack `(σ_yz, σ_yy, σ_zz)` into
columns 3,4,5 with a scale factor for conditioning.  We follow the
same pattern to keep matR well-conditioned at extreme impedance
contrasts.)

### Godunov projector — `chi`, `godunov`, `qGodLocal`, `qGodNeighbor`

Define the diagonal selector:

    chi = diag(1, 1, 1, 0, 0, 0, 0, 0, 0)        ∈ ℝ^{9×9}

This selects the local-outgoing characteristic family (columns 0,1,2
of matR).

Compute the 9×9 **Godunov projector**:

    godunov = matR · chi · matR^{-1}              ∈ ℝ^{9×9}

The Godunov interface state in the face-local frame is then

    Q*_face_local = (I − godunov) · Q_L + godunov · Q_R

where `Q_L` is the local-side state and `Q_R` is the neighbour-side
state (both already rotated into the face-local frame).  SeisSol stores
the **transposed** versions so their column-major BLAS kernels can
apply them on the right:

    qGodLocal     = I − godunov^T                ∈ ℝ^{9×9}
    qGodNeighbor  =     godunov^T                ∈ ℝ^{9×9}

so that `Q* = qGodLocal^T · Q_L + qGodNeighbor^T · Q_R`.  (We will
follow SeisSol's transposed convention if our existing
`DenseMatrix::Mult` already does the implicit transpose; otherwise we
store the non-transposed form and document the difference clearly.)

### Per-side flux

From the local element's view of the face, the contribution to its
volume residual is

    F^{local} = A^{local}_n · Q*

where `A^{local}_n` is the LOCAL-side face-normal Jacobian (built
from `λ_L, μ_L, ρ_L`).  The neighbour element, when assembling its own
volume residual for the same face, computes its OWN Q* (from its POV,
swapping the roles of local/neighbour, which gives the SAME physical
interface state Q* up to the labelling) and applies its OWN
`A^{neighbour}_n`.

This is critically different from the "one F_h per face" pattern that
homogeneous Godunov uses.  In bi-material:
  - the interface state Q* IS the same physical state from both sides;
  - the per-side flux contributions F^{local} and F^{neighbour} are
    DIFFERENT because they apply DIFFERENT Jacobians to that same Q*.

**SeisSol does NOT compute Q* explicitly** in the runtime kernel.
Instead it precomputes `AplusT = fluxScale · T^{-1} · qGodLocal · A · T`
and `AminusT = fluxScale · T^{-1} · qGodNeighbor · A · T` once per
face per side at WaveOperator construction time (see
`codegen/kernels/aderdg.py::computeFluxSolverLocal` lines 227-246).
At runtime the face flux integrand for the LOCAL element is then
just `AplusT · Q_L_buffer + AminusT · Q_R_buffer` — two
matrix-vector products against precomputed matrices.  We adopt the
same pattern (precompute per-face at ctor, apply per-QP at runtime).

### Homogeneous-limit collapse

When `(λ_L, μ_L, ρ_L) == (λ_R, μ_R, ρ_R)`, the first 3 columns of
matR equal the last 3 columns (up to sign), `godunov = matR · chi ·
matR^{-1}` collapses to the homogeneous upwind projector, and
`qGodLocal · Q_L + qGodNeighbor · Q_R` reproduces the existing
`GodunovFlux::Interior(...)` formula `A^+ · Q_L + A^- · Q_R`
bit-identically.  This is the **`T-PHASER-HOMOG-PARITY`** byte-exact
gate (R.1.T-1 below).

### Acoustic special case (`μ == 0`)

SeisSol handles acoustic (μ = 0, only P-waves) with a separate branch
inside `getTransposedGodunovState` (`ElasticSetup.h:92-139`).
**SAFS does not have acoustic regions** — the CVM-H / CVM-S /
multiscale sidecars are all bounded below at μ_min > 10⁹ Pa.  This
plan therefore ABORTS on `μ_L < ε` or `μ_R < ε` with a clear "SAFS
acoustic regions out of scope" message; the existing
`BuildGodunovFluxPool_` already aborts on `(λ + 2μ) ≤ 0` so acoustic
input would have been rejected upstream anyway.  Adding the acoustic
branch is documented as a follow-up but not implemented here.

### Reference literature

- LeVeque (2002), *Finite Volume Methods*, §22.4 ("Linear elastic
  wave equations") — the underlying derivation of the
  characteristic decomposition; used as the textbook reference for
  the math intuition.  The actual matrix form below follows the
  SeisSol implementation (which is numerically more robust and is the
  pattern proven in production).
- Käser & Dumbser (2006), "An arbitrary high-order discontinuous
  Galerkin method for elastic waves on unstructured meshes — I",
  *Geophys. J. Int.* §3.2 — origin of the SeisSol formulation.
- Pelties, de la Puente, Ampuero, Brietzke, Käser (2012), "Three-
  dimensional dynamic rupture simulation with a high-order
  discontinuous Galerkin method on unstructured tetrahedral meshes",
  *J. Geophys. Res.* §2.3 — the canonical SeisSol bi-material flux
  reference for dynamic rupture.
- Pelties, Käser, Hermann, Castro (2014), "Regular versus irregular
  meshing for complicated models and their effect on synthetic
  seismograms", *Geophys. J. Int.* §5 — quantifies the 10–30 %
  slip-rate error introduced by the average-flux approximation on
  basin/basement contrasts.

## Constraints

### Interface constraints — what cannot change
- **`GodunovFlux::Interior` signature and behaviour are preserved
  verbatim.**  All TPV / BP5 callsites continue to use it byte-
  identically.
- **The `GodunovFluxPool` API is unchanged** — `Build()`, `At(e)`,
  `NumUniqueTriples()`, `NumElements()` keep their signatures.
- **`FaultFaceFlux::EvaluateADER_LSW{,_ForcedRupture}` bi-material
  guard preserved** (`dynamic/fault_face_flux.cpp:743`).  Fault faces
  remain in the per-side-flux pattern of `FaultFaceFlux`, NOT the
  bi-material interior-flux pattern of this plan.  Fault interfaces
  and non-fault bi-material interfaces are physically distinct.
- **Phase H Stage 1 `WaveOperator(MaterialField)` ctor + per-element
  CFL + `GodunovFluxPool` build are preserved verbatim.**

### Dependency constraints
- Depends on Phase H Stage 1 (`owned_flux_pool_`, `per_elem_lmr_`,
  `shared_face_neighbour_material_`) — already committed as
  `b917fca`.
- Depends on Phase H Stage 2 (per-element flux dispatch in
  `wave_operator.inl`'s hot loop + cross-rank `MPI_Allgatherv`
  exchange) — NOT yet committed.  Phase R lands AS PART OF Phase H
  Stage 2 (same commit) so the production tree NEVER sees the
  average-flux placeholder.
- No external library dependencies beyond what GodunovFlux uses
  (`mfem::DenseMatrix` and its `LUSolver` for the 9×9 inverse).

### Convention constraints
- File naming: `dynamic/godunov_flux_bimaterial.{hpp,cpp}` (mirrors
  `godunov_flux_pool` style); helper class `BimaterialFlux`.
- Test naming: `tests/unit/test_phaser_*.cpp` mirroring
  `tests/unit/test_phaseh_*.cpp`.  Aggregate makefile target
  `test-phaser-*`.
- Sign / coordinate conventions (per `miniapps/seas/CLAUDE.md`): face
  normal `n̂` is the unit outward normal from Elem1; tangents `t̂₁ = dip`,
  `t̂₂ = strike` (Tandem `FaultBasis` convention).  Phase R inherits
  rotation from `GodunovFlux::BuildFrame` / `BuildRotation` /
  `BuildRotationInverse` and does not redefine them.
- All physical units SI.

### Numerical / performance constraints
- **Memory bound**: precompute `qGodLocal` (9×9 = 162 doubles) +
  `qGodNeighbor` (162 doubles) per HETEROGENEOUS face only.  For
  SAFS 1000m mesh (~4M tets, ~16M faces, ≲ 5 % heterogeneous), this
  is ≲ 800 K faces × 324 doubles × 8 bytes = ≲ 2 GB per rank.
  Acceptable.  If the SAFS heterogeneity exceeds 20 % of faces, fall
  back to per-QP runtime computation (acceptance gate R.3.T-5 below
  catches this).  **Phase R MUST log the heterogeneous-face count
  and the resulting per-face memory at WaveOperator ctor time**
  so the budget can be monitored.
- **Per-face dispatch cost** (runtime): two 9×9 matrix-vector
  products per QP on heterogeneous faces; homogeneous faces remain
  at one 9×9 matrix-vector product (the existing
  `Interior(...)` path).  Heterogeneous overhead ≈ 2×.
- **Byte-exact gate**: on Mode::Constant input the dispatch must
  produce a result bit-identical to the existing Phase H Stage 1
  `flux_pool_->At(e).Interior(...)` path.  This is the
  `T-PHASER-HOMOG-PARITY` gate.

## Phase R.1: Standalone bi-material flux solver

### Goal
A static helper class `BimaterialFlux` with a single public method
that, given two materials and the face normal, builds the 9×9
`qGodLocal` and `qGodNeighbor` matrices in the face-local frame.  Plus
a thin runtime `Apply(...)` method that computes the per-side flux
given those matrices and the two states.  At the end of this phase
the solver is exercised by 9 standalone unit tests; no
`wave_operator.inl` wiring yet.

### Files to Create
- `miniapps/seas/dynamic/godunov_flux_bimaterial.hpp` — header (≈ 80 LOC).
- `miniapps/seas/dynamic/godunov_flux_bimaterial.cpp` — implementation
  (≈ 280 LOC: matR construction, projector build, runtime apply).
- `miniapps/seas/tests/unit/test_phaser_bimaterial_flux.cpp` — 9 unit
  tests (≈ 400 LOC).

### Files to Modify
- `miniapps/seas/Makefile` — add `GODUNOV_FLUX_BIMATERIAL_SRC/OBJ`,
  `TEST_PHASER_BIMATERIAL_FLUX_SRC/OBJ`, the `.o` rule (uses
  `$(SEAS_INCLUDES)` only), the `seas_test_phaser_bimaterial_flux`
  link rule (deps: `GODUNOV_FLUX_OBJ` +
  `GODUNOV_FLUX_BIMATERIAL_OBJ` + `MFEM_LIBS`), the
  `test-phaser-bimaterial-flux` target, and entries in the unit-test
  aggregate + `make test` list.  Additive only; no existing recipe
  is touched.

### Detailed Requirements

1. **API surface** (in `godunov_flux_bimaterial.hpp`):

   ```cpp
   namespace mfem { namespace seas {

   /// Bi-material exact linearised Riemann solver, SeisSol formulation
   /// (Pelties et al. 2012 §2.3; SeisSol src/Equations/elastic/Model/
   /// ElasticSetup.h::getTransposedGodunovState).
   class BimaterialFlux
   {
   public:
      /// @brief Build the per-face Godunov projector matrices
      /// `qGodLocal` and `qGodNeighbor` in the FACE-LOCAL frame
      /// (n̂ → +x̂, t̂₁ → +ŷ, t̂₂ → +ẑ).
      ///
      /// `Q* = qGodLocal · Q_L_facelocal + qGodNeighbor · Q_R_facelocal`
      ///
      /// gives the bi-material interface state from the LOCAL element's
      /// POV.  Each side then applies its own A^local · Q* to get its
      /// own face flux contribution.
      ///
      /// Pre-condition: neither material is acoustic (`mu > epsilon`).
      /// Aborts otherwise — SAFS does not have acoustic regions.
      ///
      /// @param[in]  lam_self, mu_self, rho_self  Local-side material.
      /// @param[in]  lam_nbr,  mu_nbr,  rho_nbr   Neighbour-side material.
      /// @param[out] qGodLocal                     9×9 mixing matrix on Q_L.
      /// @param[out] qGodNeighbor                  9×9 mixing matrix on Q_R.
      static void BuildGodunovStateFaceLocal(real_t lam_self, real_t mu_self,
                                              real_t rho_self,
                                              real_t lam_nbr, real_t mu_nbr,
                                              real_t rho_nbr,
                                              mfem::DenseMatrix& qGodLocal,
                                              mfem::DenseMatrix& qGodNeighbor);

      /// @brief Convenience wrapper: build qGodLocal / qGodNeighbor in
      /// the face-local frame AND apply the face rotation to produce
      /// the per-side flux matrices in the GLOBAL frame.
      ///
      /// On output:
      ///   fluxLocal   = T · A_self_facelocal · qGodLocal     · T^{-1}
      ///   fluxNeighbor = T · A_self_facelocal · qGodNeighbor · T^{-1}
      ///
      /// These are the per-face matrices the wave operator multiplies
      /// against (Q_self_global, Q_neighbor_global) at runtime.  Sized
      /// 9×9 each.  Caller owns storage.
      ///
      /// `flux_self` carries the LOCAL-side material (so `A_self` is
      /// built from `flux_self.GetLambda/Mu/Rho`).  For the neighbour
      /// cell's face flux, call this function AGAIN with self/nbr
      /// swapped — the resulting fluxLocal/fluxNeighbor matrices then
      /// correspond to the neighbour cell's POV.
      static void BuildPerFaceFluxMatricesGlobal(
         const real_t* nor,
         const GodunovFlux& flux_self,
         const GodunovFlux& flux_nbr,
         mfem::DenseMatrix& fluxLocal,        ///< 9×9 output
         mfem::DenseMatrix& fluxNeighbor);    ///< 9×9 output

      /// @brief Runtime per-QP apply: F_self = fluxLocal · Q_self
      ///                                       + fluxNeighbor · Q_neighbor.
      ///
      /// Caller-allocated output buffer of size NUM_STATE.
      static void ApplyPerFaceFlux(const mfem::DenseMatrix& fluxLocal,
                                   const mfem::DenseMatrix& fluxNeighbor,
                                   const real_t* Q_self,
                                   const real_t* Q_nbr,
                                   real_t* F_h_self);

      /// @brief Dispatch gate: are the two triples "homogeneous enough"
      /// to fall through to the existing homogeneous `Interior(...)`
      /// path?  Uses 6-sig-fig rounding (same threshold as
      /// `GodunovFluxPool::Build`).
      static bool IsHomogeneous(real_t lam_self, real_t mu_self,
                                real_t rho_self,
                                real_t lam_nbr,  real_t mu_nbr,
                                real_t rho_nbr,
                                int dedup_sig_figs = 6);
   };

   }}  // namespace
   ```

2. **`BuildGodunovStateFaceLocal` algorithm** (in
   `godunov_flux_bimaterial.cpp`):

   1. Compute derived quantities:
      ```cpp
      const real_t cp_L = std::sqrt((lam_self + 2.0*mu_self) / rho_self);
      const real_t cs_L = std::sqrt(mu_self / rho_self);
      const real_t cp_R = std::sqrt((lam_nbr  + 2.0*mu_nbr)  / rho_nbr);
      const real_t cs_R = std::sqrt(mu_nbr  / rho_nbr);
      const real_t lp_L = lam_self + 2.0*mu_self;
      const real_t lp_R = lam_nbr  + 2.0*mu_nbr;
      ```
   2. Verify non-acoustic on both sides:
      ```cpp
      MFEM_VERIFY(mu_self > 1e-12 && mu_nbr > 1e-12,
                  "BimaterialFlux::BuildGodunovStateFaceLocal: acoustic input "
                  "(mu_self=" << mu_self << ", mu_nbr=" << mu_nbr
                  << ") not supported; SAFS does not have acoustic regions.  "
                  "Add the acoustic branch from "
                  "SeisSol/src/Equations/elastic/Model/ElasticSetup.h:92-139 "
                  "to enable.");
      ```
   3. Build the 9×9 `matR` following SeisSol's exact column packing
      (`ElasticSetup.h:101-139`).  Use the constants `SXX, SYY, SZZ,
      SXY, SYZ, SXZ, VX, VY, VZ` from `dynamic/wave_state.hpp`.
      **WARNING — index ordering**: SeisSol uses a state-vector
      ordering of `(σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz, u, v, w)`
      that differs in the σ_xy/σ_yz/σ_xz position from our
      `(SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ)` by one
      transposition.  Build `matR` against OUR state ordering, NOT
      SeisSol's verbatim — the eigenvector content is the same, only
      the row indices differ.  Cross-check by asserting against
      `GodunovFlux::GetAx()` in the homogeneous limit (R.1.T-1).
   4. Build `chi = diag(1, 1, 1, 0, 0, 0, 0, 0, 0)`.
   5. Compute `godunov = matR · chi · matR^{-1}` via
      `mfem::DenseMatrixInverse` (Phase H.1 has the same pattern at
      `dynamic/godunov_flux.cpp:86-89`).
   6. Output:
      ```
      qGodLocal    = I - godunov          (NOT transposed — see note below)
      qGodNeighbor =     godunov
      ```
      We use the UN-TRANSPOSED form because MFEM's `DenseMatrix::Mult`
      applies matrices column-first (left-mul on column vectors), which
      matches our existing `Tinv.Mult(Q, Q_rot)` convention in
      `GodunovFlux::Interior` at `godunov_flux.cpp:365`.  Document this
      in the function comment: "SeisSol uses transposed form because
      yateto generates code that left-mults by Q^T; MFEM uses
      column-vector-left-mul, so we store the natural form here.  The
      sign and matrix content is identical to SeisSol when transposed."

3. **`BuildPerFaceFluxMatricesGlobal` algorithm**:
   1. Build face frame (n̂, t̂₁, t̂₂) via `GodunovFlux::BuildFrame`.
   2. Build T, T^{-1} via `GodunovFlux::BuildRotation` /
      `BuildRotationInverse`.
   3. Call `BuildGodunovStateFaceLocal` to get `qGodLocal_FL`,
      `qGodNeighbor_FL` (face-local frame, ≡ "FL").
   4. Get the LOCAL-side face-x Jacobian `A_self_FL = flux_self.GetAx()`.
   5. Compute the per-face flux matrices in the GLOBAL frame:
      ```
      fluxLocal    = T · A_self_FL · qGodLocal_FL    · T^{-1}
      fluxNeighbor = T · A_self_FL · qGodNeighbor_FL · T^{-1}
      ```
      Use `mfem::Mult(T, temp, fluxLocal)` style three-matrix products.
   6. Done.  At runtime, `fluxLocal · Q_self_global + fluxNeighbor ·
      Q_nbr_global = F_self_global` (the face flux contribution to
      the LOCAL element's residual, in the global frame).

4. **`ApplyPerFaceFlux`**:
   ```cpp
   /* static */ void BimaterialFlux::ApplyPerFaceFlux(
      const mfem::DenseMatrix& fluxLocal,
      const mfem::DenseMatrix& fluxNeighbor,
      const real_t* Q_self,
      const real_t* Q_nbr,
      real_t* F_h_self)
   {
      MFEM_ASSERT(fluxLocal.Height() == NUM_STATE
                  && fluxLocal.Width() == NUM_STATE
                  && fluxNeighbor.Height() == NUM_STATE
                  && fluxNeighbor.Width() == NUM_STATE,
                  "BimaterialFlux::ApplyPerFaceFlux: input matrices must "
                  "both be NUM_STATE x NUM_STATE.");
      for (int i = 0; i < NUM_STATE; ++i)
      {
         real_t s = 0.0;
         for (int j = 0; j < NUM_STATE; ++j)
         {
            s += fluxLocal(i, j) * Q_self[j]
               + fluxNeighbor(i, j) * Q_nbr[j];
         }
         F_h_self[i] = s;
      }
   }
   ```

5. **`IsHomogeneous`** body: round each component to `dedup_sig_figs`
   significant figures with the EXACT same rounding rule as
   `dynamic/godunov_flux_pool.cpp::round_sig` (extract the helper into
   a shared `dynamic/material_dedup.hpp` rather than duplicate; the
   rounding rule must be byte-identical between `IsHomogeneous` and
   `GodunovFluxPool::Build`, otherwise the dispatch gate can disagree
   with the pool's dedup).

### Edge Cases to Handle
- **`μ_self == 0` or `μ_nbr == 0`** (acoustic): abort.  SAFS-rejected
  upstream; reject again here as defence in depth.
- **`ρ_L + 2μ_L ≤ 0`** or `ρ ≤ 0`: abort (Phase H.1 already aborts
  upstream; reject here too).
- **Identical materials**: `qGodLocal + qGodNeighbor = I`, and
  `BuildPerFaceFluxMatricesGlobal`'s `fluxLocal + fluxNeighbor =
  T · A · T^{-1} = A^{global}` matches `GodunovFlux::GetAx()` rotated
  to the global frame.  R.1.T-1 enforces this.
- **`nor` not unit length**: caller's contract (same as
  `GodunovFlux::Interior`).
- **9×9 inverse near-singular**: should not happen for non-acoustic
  isotropic materials — `matR` is well-conditioned.  Verify via
  `cond(matR) < 10^8` assertion in the body; if hit, abort with
  diagnostic.

### Acceptance Criteria

- [ ] **R.1.T-1 (homogeneous-limit byte-exact, the headline gate):**
  for 100 random `(Q_self, Q_nbr, nor)` triples and 10 random
  `(λ, μ, ρ)` triples with `flux_self == flux_nbr`,
  `BimaterialFlux::ApplyPerFaceFlux(fluxLocal, fluxNeighbor, Q_self,
  Q_nbr, F_h)` equals `GodunovFlux::Interior(nor, Q_self, Q_nbr, F)`
  bit-for-bit.  This is the `T-PHASER-HOMOG-PARITY` gate.
- [ ] **R.1.T-2 (P-wave analytic transmission coefficient):**  for
  a face with `n̂ = (1, 0, 0)` and a right-going P-wave only on the
  L side (set `Q_self[SXX] = -Z_p^L`, `Q_self[VX] = 1`, others zero;
  `Q_nbr = 0`), the resulting interface velocity is
  `v_n^* = T_p · v_n^L` with `T_p = 2 Z_p^R / (Z_p^L + Z_p^R)`.
  Verify to within `1e-12` relative for `(Z_p^L, Z_p^R) ∈ {(Z, 2Z),
  (Z, Z/3), (Z, 5Z)}`.  Extract `v_n^*` by computing
  `Q* = qGodLocal · Q_L + qGodNeighbor · Q_R` and reading the v_n
  component.
- [ ] **R.1.T-3 (P-wave analytic reflection coefficient):** same
  setup as R.1.T-2; verify the reflected-wave amplitude at the
  local side matches `R_p = (Z_p^R − Z_p^L) / (Z_p^L + Z_p^R)`.
- [ ] **R.1.T-4 (S-wave analytic transmission coefficient):** same
  as R.1.T-2 but for the `v_{t1}` channel; tests that `μ_L ≠ μ_R`
  is handled independently of `λ`.
- [ ] **R.1.T-5 (rotation invariance):** rotate `nor` to a generic
  unit vector (e.g., `(1/√3)·(1,1,1)`), apply the same rotation to
  `Q_self`, `Q_nbr` (using `GodunovFlux::BuildRotation`); the
  resulting `F_h_self` differs from the axis-aligned case by exactly
  that same rotation.  Tolerance `1e-12` relative.
- [ ] **R.1.T-6 (sign symmetry under role swap):** call with
  `(self=A, nbr=B, nor=+x̂)` then call with `(self=B, nbr=A, nor=-x̂)`
  and matching state swap; the two flux outputs agree to round-off
  AFTER accounting for the outward-normal sign flip.
- [ ] **R.1.T-7 (`qGodLocal + qGodNeighbor` invariant):**
  `qGodLocal + qGodNeighbor` equals the identity matrix to within
  `1e-12` relative — this is the mathematical invariant of the
  Godunov projector decomposition (NOT enforced by the chi pattern
  alone; comes from the bi-material eigenstructure).
- [ ] **R.1.T-8 (acoustic abort fires):** `BuildGodunovStateFaceLocal`
  called with `mu_self = 0` or `mu_nbr = 0` aborts via MFEM_VERIFY.
- [ ] **R.1.T-9 (`IsHomogeneous` rounding):** triples agreeing to 6
  sig figs → returns true; differing by 1 part in 10⁵ at any
  component → returns false.  AND the rounding decision matches
  `GodunovFluxPool::Build`'s dedup byte-identically (assert by
  building a 2-element pool with the same triples and checking
  `NumUniqueTriples == 1` iff `IsHomogeneous` returns true).

### Dependencies
- Depends on: nothing (additive, no Phase H Stage 2 dispatch needed).
- Required by: Phase R.2, R.3, R.4.

## Phase R.2: Two-layer `MaterialField` fixture

### Goal
A reusable `mfem::Coefficient` triple that returns one material on
one side of a planar interface and another on the other side.  Test
fixture for R.3 and R.5.

### Files to Create
- `miniapps/seas/tests/fixtures/two_layer_material.hpp` — header-only
  fixture class, ≈ 80 LOC.

### Detailed Requirements
1. Class `TwoLayerMaterial` holds the two material triples + the
   interface plane (z = z_interface).  Owns three `FunctionCoefficient`
   instances.
2. `MakeMaterialField()` returns a `MaterialField::MakeCoefficient(...)`
   with non-owning pointers to the three coefficients.

### Acceptance Criteria
- [ ] **R.2.T-1**: at `(x, y, z) = (0, 0, 0.3)` with `z_interface = 0.5`,
  all three coefficients return the *top* values.
- [ ] **R.2.T-2**: at `(x, y, z) = (0, 0, 0.7)`, all three return the
  *bot* values.
- [ ] **R.2.T-3**: `MaterialField` constructed via `MakeMaterialField()`
  has `mode == MaterialField::Mode::Coefficient`.

### Dependencies
- Depends on: nothing.
- Required by: Phase R.3, R.5.

## Phase R.3: WaveOperator dispatch wiring (lands with Phase H Stage 2; gated on TPV/BP5 byte-exact)

### Goal
Phase H Stage 2 (per-element flux dispatch in `wave_operator.inl`)
calls into `BimaterialFlux::ApplyPerFaceFlux` (via per-face
precomputed `fluxLocal`/`fluxNeighbor` matrices) for every
heterogeneous interior face.  Homogeneous interior faces (every
TPV / BP5 / uniform-SAFS face) short-circuit to the existing
`Interior(...)` path with byte-identical output.

**The headline acceptance gate is TPV/BP5 byte-exact regression**
(R.3.T-1 below).  Every TPV102 / TPV104 / TPV205 / BP1 / BP2 / BP5
driver is scalar-material — every face hits `IsHomogeneous == true`
and falls through to the existing scalar `Interior(...)`.  If the
TPV/BP5 byte-exact suite drifts after R.3 lands, the dispatch wiring
is wrong, and the implementer MUST STOP and bisect against pre-R.3
baseline (Commit b18bed5).  Do not advance to R.4 (analytic R/T
verification at full WaveOperator scale) or R.5 (SAFS-scale A/B)
with a red scalar-path regression — those gates verify different
classes of correctness and won't help you debug a scalar-path
drift.

**Implementation pre-condition:** Phase H Stage 2 has not yet been
committed.  Phase R.3 lands AS PART OF Phase H Stage 2 — they share
ONE commit titled `"seas/Phase H Stage 2 + Phase R: per-element
flux dispatch with exact bi-material Riemann"`.  Do NOT commit
Phase H Stage 2 with an average-flux placeholder and then replace it
with R.3; that briefly stages inferior physics into the production tree.

### Files to Modify
- `miniapps/seas/dynamic/wave_operator.hpp`:
  - Add new private member
    `std::vector<std::array<DenseMatrix, 2>> per_face_bimaterial_flux_;`
    indexed by local face index, storing
    `{fluxLocal, fluxNeighbor}` per heterogeneous face.  Empty for
    homogeneous faces (sparse storage; use
    `std::unordered_map<int, std::array<DenseMatrix, 2>>` if memory
    profiling shows the dense vector wastes too much on
    mostly-homogeneous SAFS subregions).
  - Add accessor `bool IsBimaterialFace(int local_face_idx) const`
    that returns `true` iff the face has a heterogeneous entry in
    `per_face_bimaterial_flux_`.  Used by the dispatch + by R.3.T-3.
  - Add enum + setter:
    ```cpp
    enum class BimaterialFluxMode : int { Exact = 0, Average = 1 };
    void SetBimaterialFluxMode(BimaterialFluxMode m)
    { bimaterial_flux_mode_ = m; }
    BimaterialFluxMode GetBimaterialFluxMode() const
    { return bimaterial_flux_mode_; }
    ```
    Default `Exact`.  `Average` is opt-in for the R.3.T-4 A/B test.
  - Add private helper declaration:
    `void BuildPerFaceBimaterialFluxMatrices_();` — called from the
    `(MaterialField, ...)` ctor body, AFTER `BuildGodunovFluxPool_`
    and `ExchangeBiMaterialNeighbours_`.
- `miniapps/seas/dynamic/wave_operator.inl`:
  - Implement `BuildPerFaceBimaterialFluxMatrices_`: iterate over
    every interior face (and every shared face); look up
    `(λ, μ, ρ)` for both sides (from `per_elem_lmr_` for interior,
    `shared_face_neighbour_material_` for shared); if
    `BimaterialFlux::IsHomogeneous(...)` returns true, skip; else
    call `BimaterialFlux::BuildPerFaceFluxMatricesGlobal(...)` and
    store the two 9×9 matrices in `per_face_bimaterial_flux_[face]`.
    Log the heterogeneous-face count + memory at rank 0.
  - At EVERY `flux_.Interior(...)` call site in the dispatch (audit
    via `grep -n 'flux_\.Interior' wave_operator.inl` — there are
    ≥ 2 sites: the RK4 / `ComputeFaceFluxRHS` path AND the ADER
    `ComputeADERFaceFluxRHS` path), wrap with the new dispatch:
    ```cpp
    if (owned_flux_pool_ && IsBimaterialFace(face_idx))
    {
       const auto& flux_self_mat = per_face_bimaterial_flux_[face_idx][0];
       const auto& flux_nbr_mat  = per_face_bimaterial_flux_[face_idx][1];
       real_t F_h_self[NUM_STATE];
       if (bimaterial_flux_mode_ == BimaterialFluxMode::Exact)
       {
          BimaterialFlux::ApplyPerFaceFlux(flux_self_mat, flux_nbr_mat,
                                            Q_self, Q_nbr, F_h_self);
       }
       else
       {
          // R.3.T-4 A/B comparison only.  Average of per-side
          // homogeneous Godunov fluxes.
          real_t F_self_h[NUM_STATE], F_nbr_h[NUM_STATE];
          owned_flux_pool_->At(e_self).Interior(nor, Q_self, Q_nbr,
                                                 F_self_h);
          owned_flux_pool_->At(e_nbr ).Interior(nor, Q_self, Q_nbr,
                                                 F_nbr_h);
          for (int c = 0; c < NUM_STATE; ++c)
          { F_h_self[c] = 0.5 * (F_self_h[c] + F_nbr_h[c]); }
       }
       // Accumulate F_h_self into elem_self's residual.  The
       // NEIGHBOUR element will independently invoke this dispatch
       // from its OWN ComputeFaceFluxRHS pass and produce its own
       // F_h_nbr (different value, because A^nbr differs from A^self).
       ...
    }
    else if (owned_flux_pool_)
    {
       // Homogeneous heterogeneous-material face: use per-element
       // pool's Interior(), byte-equivalent to scalar path.
       owned_flux_pool_->At(e_self).Interior(nor, Q_self, Q_nbr, F_h);
    }
    else
    {
       // Scalar-material path: UNCHANGED — TPV / BP5 byte-identical.
       flux_.Interior(nor, Q_self, Q_nbr, F_h);
    }
    ```

  **CRITICAL invariant**: in MFEM's element-local DG assembly, each
  cell computes ITS OWN face flux for ITS OWN element residual.  The
  LOCAL cell's POV produces `F_h_self`; when the NEIGHBOUR cell does
  ITS OWN assembly for the same face, it runs the same dispatch from
  ITS POV and produces a DIFFERENT `F_h_self` (which is the
  "F_h_nbr" from the original cell's POV).  This means
  `per_face_bimaterial_flux_[face][0]` from the LOCAL cell's POV is
  NOT THE SAME MATRIX as `per_face_bimaterial_flux_[face][0]` from
  the NEIGHBOUR cell's POV — they use different `A_self` (the
  cell's own Jacobian).  We therefore store BOTH POVs:
  ```cpp
  std::vector<std::array<DenseMatrix, 4>> per_face_bimaterial_flux_;
  //                              [0]: fluxLocal    (cell 1's POV)
  //                              [1]: fluxNeighbor (cell 1's POV)
  //                              [2]: fluxLocal    (cell 2's POV)
  //                              [3]: fluxNeighbor (cell 2's POV)
  ```
  Memory cost doubles to 4 × 81 = 324 doubles per heterogeneous face;
  still ~2 GB at the SAFS scale.  Document this clearly in the
  comment block above the data structure.
- `miniapps/seas/drivers/spatial_dyn_driver.cpp` — add
  `--bimaterial-flux {exact|average}` CLI flag, default `exact`.

### Detailed Requirements
1. **Audit all `flux_.Interior(...)` call sites** in
   `wave_operator.inl` via `grep -n 'flux_\.Interior'` BEFORE writing
   the patch.  Each site gains the new dispatch.  This is the
   `feedback_complete_sign_sites` memory rule applied to flux
   dispatch.
2. **Homogeneous short-circuit MUST agree with `GodunovFluxPool::Build`'s
   dedup** — call `BimaterialFlux::IsHomogeneous(...)` which uses the
   shared `dynamic/material_dedup.hpp` rounding helper.  Otherwise
   the dispatch can route a face to the bi-material branch while the
   pool deduped it to one entry (waste) or route a face that the
   pool kept separate to the homogeneous branch (wrong).
3. **Shared-face POV symmetry**: the cross-rank
   `shared_face_neighbour_material_` map (Phase H.4) supplies the
   neighbour-side `(λ, μ, ρ)`.  Use it on shared faces in place of
   `per_elem_lmr_`.  Document the symmetry contract: both ranks
   that own a shared face must build BIT-IDENTICAL
   `per_face_bimaterial_flux_` entries because they both round
   inputs the same way and run the same matR construction.  A
   `wave.VerifySharedFaceBimaterialFlux()` helper (analogous to
   `VerifySharedFaultDOFDataConsistency`) gathers all heterogeneous
   shared-face flux matrices to root and asserts bit-equality
   across the (rank A, rank B) pair.  Called once after R.4 dispatch
   wiring lands.
4. **DR / boundary faces are SKIPPED** in
   `BuildPerFaceBimaterialFluxMatrices_`: the existing fault-face
   bookkeeping (`fault_interior_faces_`, `fault_shared_faces_`)
   identifies them; they go through `FaultFaceFlux` not this
   dispatch.  Boundary faces are gated on `face_bdr_attr_[face] > 0`
   in `wave_operator.inl` and use the absorbing / free-surface BC
   dispatch instead.

### Edge Cases to Handle
- **Fault face in a bi-material mesh**: skipped (handled by
  `FaultFaceFlux`'s existing `homog_ok` guard).
- **One side at the absorbing boundary**: covered by the absorbing
  BC dispatch; bi-material only applies to two-sided interior faces.
- **Shared face with > 2 ranks (corner / edge)**: ParMesh handles
  this via face-neighbour topology; each rank sees exactly one
  neighbour per shared face.  No special case.
- **Memory budget exceeded**: if rank-0 log reports
  `per_face_bimaterial_flux_` consuming > 4 GB on this rank, emit a
  rank-0 WARNING and recommend the user partition more aggressively
  or switch to per-QP runtime computation (a follow-up flag
  `--bimaterial-flux-storage runtime` we do NOT implement here but
  document).

### Acceptance Criteria
- [ ] **R.3.T-1 — HEADLINE GATE: TPV / BP5 byte-exact regression.**
  `make test-tpv104`, `make test-tpv205`, `make test-tpv102-local`
  (if currently green on this branch), `make test-bp5-fault-operator`
  (if currently green on this branch), `make test-bp5-integration`,
  `make test-bp5-smoke` all produce bit-identical output to the
  pre-Phase-R baseline (Commit `b18bed5`).  Same gates as Phase H
  Stage 1's `T-PHASEH-SCALAR-PARITY`.  Every TPV / BP5 face is
  scalar-material → `IsHomogeneous == true` → fall through to the
  existing scalar `flux_.Interior(...)`; any byte-level drift means
  the new dispatch branch leaked into the scalar path.  **If this
  gate is red, STOP and bisect — do not advance to R.4.**
- [ ] **R.3.T-2 (homogeneous SAFS byte-exact regression)**:
  `seas_spatial_dyn_driver --no-sidecar-material` produces dispatch
  through `IsHomogeneous == true` on every face; dry-run gives the
  same `dt_cfl` and the same first-step `Q_new` as the pre-Phase-R
  baseline (Commit `b18bed5`).
- [ ] **R.3.T-3 (heterogeneous dispatch counter)**: with
  `SEAS_DIAG_PHASER` build, on the R.4 two-layer mesh (`nz = 64`,
  `Z_R/Z_L = 2.0`), the count of `ApplyPerFaceFlux` calls per
  `Mult` equals the count of interface faces (= `nx · ny ·
  nbf_per_face = 16 · nbf_per_face`).  Proof that the dispatch
  reaches the heterogeneous branch on the test fixture R.4 will use.
- [ ] **R.3.T-4 (driver flag wiring)**: `--bimaterial-flux average`
  reproduces the average-flux behaviour bit-identically (used by
  R.4.T-4 for the A/B convergence comparison); `--bimaterial-flux
  exact` routes through the new exact path.
- [ ] **R.3.T-5 (memory budget)**: on the SAFS 1000m
  `_lcfar3000.msh` + CVM-H sidecar, rank-0 log reports
  heterogeneous-face count < 5 % of total interior faces AND
  `per_face_bimaterial_flux_` total < 4 GB per rank.
- [ ] **R.3.T-6 (shared-face POV symmetry)**: at np = 4 on the SAFS
  fixture, `wave.VerifySharedFaceBimaterialFlux()` reports zero
  bit-mismatches across the rank pair that owns each shared
  heterogeneous face.

### Dependencies
- Depends on: Phase R.1, Phase R.2 (the two-layer fixture is used
  by R.3.T-3's dispatch counter check), Phase H Stage 2 (lands in
  the same commit).
- Required by: Phase R.4, Phase R.5.

## Phase R.4: Layered-medium analytic R/T verification (gates SAFS-scale R.5)

### Goal
End-to-end verification on a 1-D layered medium: a Gaussian P-wave
incident normally onto a flat bi-material interface produces the
correct analytic reflected / transmitted amplitudes, and the L2
error converges at the expected FE order.  This is the **only**
place in Phase R where the bi-material physics is verified against
an analytic solution at full `WaveOperator::Mult` scale.  R.1's
unit tests cover the math at a single face; R.5's A/B can show
"exact ≠ average" but not "exact = correct".  R.4 closes that gap.

### Pre-condition: R.3 must be green before R.4 starts
R.4 requires the dispatch (R.3) to actually reach
`ApplyPerFaceFlux` on heterogeneous faces.  If R.3.T-3
(dispatch counter) does NOT report the expected count of
`ApplyPerFaceFlux` calls per Mult on the two-layer test mesh,
R.4 cannot produce a meaningful result and must wait for the
dispatch wiring to be fixed.

### Files to Create
- `miniapps/seas/tests/unit/test_phaser_layered_p_wave.cpp` — driver
  that builds an inline Cartesian mesh, sets up `WaveOperator` with
  a two-layer `MaterialField`, evolves a Gaussian P-wave pulse
  through the interface, and compares against analytic `R_p`, `T_p`.

### Files to Modify
- `miniapps/seas/Makefile` — add the new test target.

### Detailed Requirements
1. Inline Cartesian box mesh `Mesh::MakeCartesian3D(4, 4, nz, …)`
   with `nz ∈ {16, 32, 64, 128}` (the convergence sweep).
2. `TwoLayerMaterial` (from R.2) with `z_interface = 0.5`, top
   `(λ, μ, ρ) = (1, 1, 1)`, bot `(λ, μ, ρ) = (4, 4, 1)` → `Z_R /
   Z_L = 2.0`.
3. Initial state: Gaussian P-wave pulse centred at `z = 0.2` with
   width `σ = 0.05`, normalised so peak `v_z = 1`.
4. Evolve for `t = 0.4` (pulse reaches interface and reflects);
   compare reflected pulse amplitude at `z = 0.1` against `R_p ·
   A_incident` and transmitted at `z = 0.8` against `T_p ·
   A_incident`, where:
   ```
   R_p = (Z_p^R − Z_p^L) / (Z_p^L + Z_p^R)
   T_p = 2 · Z_p^R       / (Z_p^L + Z_p^R)
   ```
   are LeVeque §22.4's analytic P-wave reflection / transmission
   coefficients for normal incidence at a planar contact.
5. Convergence sweep: log-log fit L2 error vs `h = 1/nz`; assert
   slope ∈ `[0.8, 1.5]` for `nz ∈ {32, 64, 128}` (the FE order
   `p = 1` ± slack for discrete pulse-tracking error).
6. **A/B comparison gate**: also run the same sweep with
   `wave.SetBimaterialFluxMode(Average)` and assert the
   average-flux L2 error is ≥ 10× larger than the exact-Riemann
   path on the finest mesh.  This is the empirical proof that the
   exact-Riemann path is the production-correct choice.

### Edge Cases to Handle
- **Pulse hits the absorbing boundary** before the back-reflection
  is captured: choose `tfinal = 0.4` so the transmitted pulse has
  not yet reached the back wall (`z = 1.0`) given `c_p ≈ √3`.
- **Pulse width below the mesh resolution** on `nz = 16`: with
  `σ = 0.05` we have ~ 4 cells per σ at `nz = 16`, barely resolved.
  The lowest `nz` point may show pre-asymptotic convergence — fit
  the slope only over `nz ∈ {32, 64, 128}`.
- **Sign-of-V_z convention** on the reflected pulse: the
  reflected wave's `v_z` has SIGN `−R_p` relative to the incident
  `v_z` (because both incident and reflected fluxes carry
  particle motion in the same physical direction for a low-to-high
  impedance step).  Reference: LeVeque (22.21).  Verify by hand
  computation on the test fixture before asserting.

### Acceptance Criteria
- [ ] **R.4.T-1**: `|R_observed − R_analytic| / |R_analytic| <
  0.02` on the finest mesh (`nz = 128`).
- [ ] **R.4.T-2**: `|T_observed − T_analytic| / |T_analytic| <
  0.02` on the finest mesh.
- [ ] **R.4.T-3**: log-log fit of L2 error vs `h = 1/nz` over
  `nz ∈ {32, 64, 128}` returns slope ∈ `[0.8, 1.5]`.
- [ ] **R.4.T-4 (regression vs average flux)**: same setup with
  `--bimaterial-flux average` (the R.3.T-4 CLI flag) produces L2
  error ≥ 10× larger than the exact-Riemann path on the finest
  mesh.  Documented in the R.4 test log so future maintainers can
  see the empirical justification for the exact path.

### Dependencies
- Depends on: Phase R.1 (solver math), Phase R.2 (two-layer
  fixture), Phase R.3 (dispatch wired into `wave_operator.inl`).
  All three must be green before R.4 starts.
- Required by: Phase R.5 (SAFS-scale acceptance).

## Phase R.5: SAFS-scale acceptance + plan cross-reference cleanup

### Goal
Verify on the actual SAFS 1000m mesh + CVM-H sidecar that exact
bi-material Riemann produces a meaningful difference vs average-flux
at the basement/basin interface AND converges to a refined-mesh
reference faster than the average-flux path.  Then update the parent
plans to remove their "future work" caveats.

### Files to Create
- `miniapps/seas/spatial/code/scripts/verify_bimaterial_riemann.py`
  — A/B driver wrapper: run `seas_spatial_dyn_driver` once with
  `--bimaterial-flux exact` and once with `--bimaterial-flux
  average`, diff the resulting `volume.vtkhdf` files, report max
  relative difference per field + at the embedded fault's slip-rate
  trace.

### Files to Modify
- `miniapps/seas/safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md`
  — line 586's cross-reference to `PLAN_phase_R_exact_bimaterial_riemann.md`
  (this PLAN) is already in place from the rev-2 cleanup; verify it
  matches.
- `miniapps/seas/safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md`
  line 579 — the "heterogeneous_material_plan Phase 3 Detailed
  Req. 6 (average flux)" reference is updated to point at Phase R
  instead.
- `miniapps/seas/CLAUDE.md` — append note in "Known limitation"
  block: bi-material flux dispatch is `Exact` by default;
  `--bimaterial-flux average` is for A/B verification only.
- `miniapps/seas/dynamic/heterogeneous_material.hpp` — if the
  comment block claims "average flux at heterogeneous interior
  faces", update to "exact linearised bi-material Riemann; see
  Phase R".

### Detailed Requirements
1. Run on SAFS `1000m_lcfar3000` mesh + CVM-H sidecar, `--tfinal
   1.0s`, `--ader-order 2`, `--cfl 0.5`.  Compare `volume.vtkhdf`
   between the two flux modes.
2. Plot a slip-rate trace from a station near the basement/basin
   interface; the two traces should differ by ≥ 5 % peak amplitude
   (proof that the exact path matters), AND the exact path's trace
   should be closer to a 2×-refined-mesh reference than the
   average-path's trace.

### Acceptance Criteria
- [ ] **R.5.T-1**: `verify_bimaterial_riemann.py` reports max
  relative difference between the two flux modes of ≥ 5 % on at
  least one bulk-velocity component near the basement/basin
  interface.
- [ ] **R.5.T-2**: on the 2×-refined mesh the exact-flux trace is
  within 1 % of the refined-mesh reference; the average-flux trace
  differs by ≥ 5 %.
- [ ] **R.5.T-3**: `spatial_dynamic_rupture_plan.md:586` carries
  the Phase R cross-reference (verify it matches this file's
  path).
- [ ] **R.5.T-4**: `heterogeneous_material_plan.md` Phase 3
  Detailed Req. 6 cross-references this plan.

### Dependencies
- Depends on: Phase R.4.
- Required by: nothing (closing phase).

## Testing Strategy

| Commit | Phase(s) | What to test | How to validate |
|--------|----------|--------------|------------------|
| **#1** | R.1 + R.2 | Standalone solver math + two-layer fixture | 9 R.1 unit tests (homogeneous byte-exact vs `GodunovFlux::Interior`, analytic P/S transmission + reflection, rotation invariance, role-swap symmetry, `qGodLocal + qGodNeighbor = I` invariant, acoustic abort, `IsHomogeneous` rounding) + 3 R.2 fixture tests. |
| **#2** | R.3 (+ Phase H Stage 2) | Dispatch wiring + TPV/BP5 byte-exact + dispatch counter on the R.2 fixture | **R.3.T-1 (TPV/BP5 byte-exact) is the headline gate.**  Plus homogeneous-SAFS dry-run byte-exact, dispatch counter, CLI flag round-trip, memory budget log, shared-face POV symmetry. |
| **#3** | R.4 | Layered-medium analytic R/T at full `WaveOperator::Mult` scale | Gaussian P-wave pulse on inline 3-D mesh + analytic `R_p` / `T_p` to within 2 % on finest mesh + log-log L2 convergence slope in `[0.8, 1.5]` + A/B vs `--bimaterial-flux average` shows ≥ 10× lower error on the exact path. |
| **#4** | R.5 | SAFS-scale impact + parent-plan cleanup | A/B Python script on SAFS 1000m + CVM-H + mesh refinement convergence vs reference + parent plans updated. |

Numerical correctness chain (each commit enables the next):

1. **Commit #1 — R.1 + R.2** prove the solver math is right at one
   face (the homogeneous byte-exact gate is the strictest math check
   the codebase has) and produce the test fixture R.3 + R.4 + R.5
   all consume.
2. **Commit #2 — R.3 + Phase H Stage 2** prove the dispatch wiring
   does not drift the scalar path (TPV/BP5 byte-exact) and DOES
   reach `ApplyPerFaceFlux` on heterogeneous faces (dispatch
   counter on the R.2 fixture).  Until this commit is green, the
   scalar code path is potentially broken AND there's no production
   way to exercise the new bi-material physics.
3. **Commit #3 — R.4** proves the solver integrates correctly
   inside a `WaveOperator::Mult` and reproduces the analytic
   `R_p` / `T_p`.  This is the only math-correctness gate at
   integrated scale (R.1's unit tests verify one face;
   R.5's A/B can show "different" but not "correct").
4. **Commit #4 — R.5** proves the production-scale impact is
   measurable on SAFS + CVM-H AND the exact path converges to a
   refined-mesh reference faster than the average-flux path.

## Risk Assessment

### Things that have caused real bugs in similar code

1. **State-vector index ordering mismatch with SeisSol's eigenvector
   matR**.  Our state order `(SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ)`
   differs from SeisSol's `(σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz, u, v, w)`
   only nominally, but the σ_xy / σ_yz / σ_xz row indices in matR's
   eigenvector entries (rows 3, 4, 5 vs 3, 5, 4) MUST be checked
   against `dynamic/wave_state.hpp`.  Mitigation: R.1.T-1 (homogeneous
   byte-exact vs `GodunovFlux::Interior`) catches any index swap on
   every random `(λ, μ, ρ, nor, Q)`.
2. **Transposed-vs-untransposed Godunov projector**.  SeisSol stores
   the transposed form because their yateto code-gen left-mults by
   row vectors.  MFEM's `DenseMatrix::Mult` right-mults column
   vectors, so we store the natural form.  If R.1.T-1 fails by 1
   row/column swap, this is the culprit.
3. **POV asymmetry in shared-face matrices**.  Each rank that owns
   a shared face must build the SAME flux matrices.  Bit-identical
   inputs (same rounded `(λ, μ, ρ)`, same `nor` direction relative
   to each rank's local `Elem1`) give bit-identical outputs ONLY
   IF the matR construction is deterministic.  `DenseMatrixInverse`'s
   LU pivoting CAN reorder operations; verify via R.3.T-6 which
   gathers shared-face matrices to root and bit-compares.
4. **`IsHomogeneous` threshold drift vs `GodunovFluxPool::Build`'s
   dedup**.  These MUST share the rounding helper (extract to
   `dynamic/material_dedup.hpp`).  Failure mode: a face the pool
   deduplicated to one entry routes through the bi-material branch
   (wasted work, no correctness issue) OR a face the pool kept
   distinct routes through the homogeneous branch (wrong physics —
   exactly the bug we're trying to prevent).  R.1.T-9 catches the
   threshold drift at unit-test scale BEFORE R.3 lands; the gate
   is therefore caught one commit earlier than it would be on the
   dispatch side.
5. **Memory blowup on SAFS-scale meshes**.  Per-face precomputation
   is what makes the SeisSol approach fast at runtime, but it costs
   324 doubles per heterogeneous face.  If the SAFS basin / basement
   ratio is more heterogeneous than expected (e.g., 30 % of faces
   instead of 5 %), 4 M tets × 4 faces × 0.3 × 324 doubles × 8 bytes
   ≈ 12 GB per rank.  R.3.T-5 catches this; mitigation is to fall
   back to per-QP runtime computation (a follow-up flag we do NOT
   implement in this commit).
6. **TPV / BP5 regression from the dispatch insertion**.  Every
   `flux_.Interior(...)` site gains a new `if (owned_flux_pool_ &&
   ...)` branch.  For scalar-material runs `owned_flux_pool_ ==
   nullptr` and the branch is dead-code-eliminated, so cost is one
   nullptr compare per face.  **R.3.T-1 (TPV/BP5 byte-exact) is the
   headline gate** for this risk — it is the GO/NO-GO criterion
   for advancing from Commit #2 to Commit #3.  If R.3.T-1 fails,
   STOP and bisect against pre-Phase-R baseline (Commit `b18bed5`);
   do NOT advance to R.4 (layered medium) — the layered test
   exercises the heterogeneous branch but doesn't tell you whether
   the scalar branch is broken.

### Known tricky areas in existing code
- The fault-face dispatch in `wave_operator.inl` is the most
  intricate piece of the file (R-016 / R-1600 / R-1601 + per-substep
  imposed-state pointer pair).  Phase R does NOT touch the
  fault-face branch — only the non-fault interior-face branch.
  `grep -n 'flux_\.Interior' wave_operator.inl` gives the exact list
  of call sites to wrap; the fault branch consumes
  `fault_flux_->EvaluateADER*` instead and is excluded.
- Phase H Stage 2 + Phase R is the highest-blast-radius commit in
  the multi-phase rollout.  Drift policy: STOP and bisect — do not
  auto-rebaseline.

## Out of Scope for This Plan
1. **Bi-material at fault faces** (`FaultFaceFlux::Evaluate*`).
   Fault interfaces are physically distinct (carry friction law +
   slip DOF); the per-side flux of a bi-material fault is the topic
   of separate ongoing research (Harris & Day 2005) and is out of
   scope.
2. **Anisotropic material** (`(λ, μ)` replaced by full elasticity
   tensor).  SeisSol supports anisotropy via
   `getRotatedMaterialCoefficients`; for isotropic SAFS this is a
   no-op.  Anisotropic SAFS is a future plan.
3. **Acoustic regions** (μ = 0).  SeisSol supports them via a
   separate branch in `getTransposedGodunovState`; SAFS does not
   have acoustic regions so we abort instead.  Add the acoustic
   branch when needed (the SeisSol implementation in
   `ElasticSetup.h:92-139` is the direct port reference).
4. **Anelastic / viscoelastic / poroelastic equations**.  SeisSol
   has separate code paths in `src/Equations/{viscoelastic,
   poroelastic}/`; out of scope for SAFS Phase R.
5. **GPU port** — same as parent plan.
6. **Per-face precomputation FALLBACK to per-QP runtime** when
   memory budget is exceeded.  Documented as the mitigation in
   R.3.T-5 but not implemented; a follow-up flag
   `--bimaterial-flux-storage runtime` would add the per-QP path.

## Cross-references
- `spatial_dynamic_rupture_plan.md` Phase H (heterogeneous WaveOperator).
- `heterogeneous_material_plan.md` Phase 3 (MaterialField — the
  legacy average-flux note that Phase R replaces).
- `miniapps/seas/CLAUDE.md` "Known limitation — Gmsh `.msh` format"
  section (mesh format gate for SAFS production runs).
- `dynamic/godunov_flux.{hpp,cpp}` (the homogeneous Godunov flux
  Phase R extends).
- `dynamic/godunov_flux_pool.{hpp,cpp}` (the per-element flux cache
  Phase R consumes via `IsHomogeneous`).
- `dynamic/wave_operator.{hpp,inl}` (the dispatch site Phase R wires).

### SeisSol reference paths (read at planning time)
- `/Users/chunhuizhao/projects/SeisSol/src/Equations/elastic/Model/ElasticSetup.h:80-166`
  — `getTransposedGodunovState`: the canonical bi-material Riemann
  builder.  This plan's math section reproduces the matR packing
  verbatim.
- `/Users/chunhuizhao/projects/SeisSol/src/Initializer/CellLocalMatrices.cpp:200-355`
  — `initializeCellMatrices`: where SeisSol calls
  `getTransposedGodunovState` per cell-face at init time and stores
  the resulting matrices into `localIntegration[cell].nApNm1[side]`.
  Our `BuildPerFaceBimaterialFluxMatrices_` mirrors this.
- `/Users/chunhuizhao/projects/SeisSol/codegen/kernels/aderdg.py:227-246`
  — `computeFluxSolverLocal` / `computeFluxSolverNeighbor`: the
  auto-generated kernel that bakes
  `fluxScale · T^{-1} · qGodLocal · A · T` into a single per-face
  matrix.  Our `BuildPerFaceFluxMatricesGlobal` produces the same
  composition.
- `/Users/chunhuizhao/projects/SeisSol/src/Model/Common.h:33`
  — `testIfAcoustic(mu)`: the convention SeisSol uses to detect
  acoustic input.  We use the same threshold (`|mu| ≤ eps`) for
  the acoustic-abort guard in R.1.
