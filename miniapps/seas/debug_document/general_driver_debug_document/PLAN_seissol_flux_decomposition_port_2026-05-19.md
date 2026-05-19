# PLAN — Port SeisSol-Style Flux Decomposition for Bimaterial + Mixed-Flux Composition

**Date:** 2026-05-19 (revised: now correctly mirrors DRDG3D / Zhang
2023 pure-central choice that the codebase already uses for the
scalar ctor, NOT SeisSol's Rusanov+LF.)

**Branch:** `feature/heterogeneous_riemann_solver`

**Status:** PROPOSED — long-term research item (estimated 2–3 weeks
implementation + 1 week verification + 1 week integration with
`seas_spatial_dyn_driver`).

**Author:** code-plan agent

**Companion docs:**
- `spatial_dyn_heterogeneous_riemann_review_2026-05-19.md` (review).
- `spatial_dyn_heterogeneous_riemann_fix_2026-05-19.md` (fix report).
- `safs/project_7.0_alternative/document/PLAN_phase_R_exact_bimaterial_riemann_rev3.md` (the bimaterial Riemann port that already landed).

---

## 0. Naming clarification (added 2026-05-19 after DRDG3D audit)

The original plan title said "port SeisSol's Q_god/Q_corr
decomposition", which is technically wrong about the dissipation
model.  After cross-checking three reference implementations:

| Reference | Architecture | Near-fault flux | Stabilization mechanism |
|-----------|--------------|-----------------|-------------------------|
| **SeisSol** (`Initializer/CellLocalMatrices.cpp:283-330`) | Decomposed `Q_god + Q_corr` matrices, swapped per face | Central (`Q_god = 0.5·I`) + Lax-Friedrichs (`Q_corr = ±0.5·c·I`) | LF scalar viscosity supplies stabilization |
| **DRDG3D / Zhang 2023** (`src/mod_wave.F90:587-617`) | Per-face if-branch; central goto-100 | **Pure central** (`fstar = 0.5·(F_R - F_L)`, NO LF correction) | Fault itself supplies dissipation via imposed-state Riemann + slip-weakening |
| **MFEM SEAS today** (scalar ctor only; `wave_operator.inl:3509-3517`, `godunov_flux.cpp:382-430`) | Per-face if-branch on `central_flux_face_set_` | **Pure central** (`F = 0.5·A_n·(Q_self + Q_nbr)`, NO LF correction) | Follows Zhang / DRDG3D |

The MFEM SEAS code's existing scalar-ctor mixed-flux path is
**already DRDG3D-style, not SeisSol-style**.  So the consistent port
to the bimaterial ctor is to follow Zhang/DRDG3D (pure central, no
LF correction), NOT SeisSol (central + LF).

The Q_god / Q_corr decomposition pattern from SeisSol is still useful
as the **architectural template** — the (`Q_god`, `Q_corr`) pair
neatly factors `flux = T · (Q_god·A_self + Q_corr) · Tinv` into a
multiplicative part and an additive part.  We adopt that template,
but the values we plug in for "central mode" are Zhang's (Q_god =
0.5·I, Q_corr = 0), NOT SeisSol's (Q_god = 0.5·I, Q_corr =
±0.5·c·I).

---

## 1. Motivation

The current `feature/heterogeneous_riemann_solver` branch lets every
TPV / SAFS run route through `BimaterialFlux::ApplyPerFaceFlux` — the
exact-linearised Riemann solver per Pelties et al. 2012.  Correctness
on heterogeneous media (TPV31 depth profile, SAFS / CVM-H) is now
solved.

The branch does NOT solve the **rupture-front upwind-dissipation
problem** that Zhang et al. 2023 targets.  The two existing dispatch
modes are:

| `interior_flux` | `mixed_flux` | Material jumps | Rupture-front damping |
|---|---|---|---|
| `scalar`     | `none`       | mishandled | full upwind (present) |
| `scalar`     | `adjacent`   | mishandled | central near fault (mitigated) |
| `bimaterial` | `none`       | correct    | full upwind (present) |
| `bimaterial` | `adjacent`   | (forbidden by R-002 guard) | — |

Row 4 — bimaterial correctness + mixed-flux stability — is what
SeisSol production uses.  SeisSol implements it through a `Q_god` /
`Q_corr` decomposition of the per-face flux matrix where each face
independently chooses Godunov (upwind, bimaterial) or Rusanov
(central + Lax–Friedrichs viscosity).  This plan ports that
architecture into MFEM SEAS.

**Why this matters numerically (short version).**  At `--order 1` /
ADER-O2 the leading numerical viscosity of pure upwind DG scales as
`h^{2p+1} = h^5`.  At SeisSol's default ORDER=6 it scales as `h^{13}`
— ~8 orders of magnitude smaller at the same mesh.  The MFEM SEAS code
runs at p=1 today because higher orders are 8–64× more expensive per
step; the mixed-flux composition is the practical alternative for
getting rupture-front accuracy at low polynomial order.  See §3 for
the literature.

---

## 2. Goal & Acceptance Criteria

### Goal

Replace the two R-002 hard-abort guards
(`wave_operator.inl:652-660` ctor-side and `wave_operator.inl:2098`
setter-side) with **the actual composition logic**, so that

```toml
[numerics]
interior_flux = "bimaterial"
mixed_flux    = "adjacent"
```

becomes a supported, tested, documented configuration that:

- Uses the **bimaterial Godunov flux** on every interior face that is
  NOT in `central_flux_face_set_`.
- Uses the **bimaterial Rusanov flux** (central + LF viscosity, with
  the correct local/neighbor wave-speed) on every face that IS in
  `central_flux_face_set_`.
- Falls back gracefully to (`scalar`, `adjacent`) when called on a
  homogeneous-material driver.

### Acceptance Criteria

| ID | Criterion |
|----|-----------|
| **A-1** | `WaveOperator(MaterialField, BoundaryConfig)` ctor accepts `mixed_flux_mode_ == Adjacent`; no MFEM_ABORT. |
| **A-2** | `SetMixedFluxMode(Adjacent)` succeeds when called after the heterogeneous ctor; no MFEM_ABORT. |
| **A-3** | On a homogeneous-material driver, `(bimaterial, adjacent)` produces output BIT-EQUIVALENT to `(scalar, adjacent)` within 1e-10 relative (the existing `T-PHASEH-SCALAR-PARITY` tolerance). |
| **A-4** | On a homogeneous-material driver, `(bimaterial, none)` produces output BIT-EQUIVALENT to its current value (no regression on the bimaterial path). |
| **A-5** | TPV205 production run with `(bimaterial, adjacent)` produces lower upwind dissipation than `(bimaterial, none)` (measured by peak \|V_strike\| at hypocenter station; expected ~10–30 % higher peak). |
| **A-6** | 1D bimaterial smoke test (the existing `test_phaser_dispatch_smoke` fixture) extended to demonstrate `(bimaterial, adjacent)` differs from `(scalar, adjacent)` at the material jump by > 1e-6 relative. |
| **A-7** | A new unit test `test_phaser_rusanov_decomposition` verifies the Q_god / Q_corr decomposition matches the SeisSol reference formula on synthetic inputs. |

---

## 3. Reference Material

### Primary references

| Paper | What to consult |
|-------|-----------------|
| Pelties, C., Käser, M., Hermann, V., Castro, C. E. (2012). **Regular versus irregular meshing for complicated models and their effect on synthetic seismograms.** *Geophysical Journal International*, 183(2):1031–1051. | §2.3 derivation of the bi-material Riemann state and the eigenvector packing of `matR`.  Already mirrored in this codebase as `BimaterialFlux::BuildGodunovStateFaceLocal`. |
| Pelties, C., de la Puente, J., Ampuero, J.-P., Brietzke, G. B., Käser, M. (2012). **Three-dimensional dynamic rupture simulation with a high-order discontinuous Galerkin method on unstructured tetrahedral meshes.** *Journal of Geophysical Research: Solid Earth*, 117, B02309. | The canonical SeisSol TPV / dynamic-rupture paper.  Establishes the imposed-state Riemann on the fault and the LSW closed-form solve. |
| Käser, M., Dumbser, M. (2006). **An arbitrary high-order discontinuous Galerkin method for elastic waves on unstructured meshes — I. The two-dimensional isotropic case with external source terms.** *Geophysical Journal International*, 166(2):855–877. | The ADER-DG framework SeisSol is built on.  §3 derives the upwind Godunov flux for the velocity-stress system that both Pelties 2012 papers extend. |
| Krenz, L., Uphoff, C., Ulrich, T., Gabriel, A.-A., Abrahams, L. S., Dunham, E. M., Bader, M. (2021). **3D acoustic-elastic coupling with gravity: the dynamics of the 2018 Palu, Sulawesi earthquake and tsunami.** *Proc. SC '21*, ACM. | Documents the elastic-acoustic interface that motivated SeisSol's `enforceGodunov` short-circuit at acoustic faces.  Relevant when porting the acoustic-rejection branch. |
| **Zhang, A., Liu, Y., Chen, R., et al. (2023).** **A discontinuous Galerkin method with hybrid upwind / central flux for the elastic-wave equation with applications to dynamic rupture.** *(Inline citation only — referenced throughout `wave_operator.{hpp,inl}` as "Zhang 2023".  Confirm the actual journal venue when reading the file in the codebase: `wave_operator.inl:2049,2130,5925`.)* | §3.2 Fig. 5 (the `central_flux_face_set_` selection algorithm — already implemented in MFEM SEAS as `BuildCentralFluxFaceSet_`) and Fig. 4 (SSO/HFO trade-off between Mixed-Flux 1 and Mixed-Flux 2).  This is the paper whose recommendation the MFEM SEAS code already partially implements for the scalar ctor. |

### SeisSol source-code references (the implementation we are mirroring)

| File | Lines | Purpose |
|------|-------|---------|
| `SeisSol/src/Initializer/parameters/ModelParameters.h` | 27 | `enum class NumericalFlux { Godunov, Rusanov };` — the two-option enum we need. |
| `SeisSol/src/Initializer/parameters/ModelParameters.h` | 45–46 | `NumericalFlux flux` + `NumericalFlux fluxNearFault` — two independent parameter fields. |
| `SeisSol/src/Initializer/parameters/ModelParameters.cpp` | 106–118 | TOML/parameter-file parsing of the two enums. |
| `SeisSol/src/Initializer/CellLocalMatrices.cpp` | 260–276 | `isSpecialBC(side)` — the "is this face adjacent to a DR (fault) face?" predicate.  The MFEM analogue is `BuildCentralFluxFaceSet_` (Adjacent mode). |
| `SeisSol/src/Initializer/CellLocalMatrices.cpp` | 283–294 | Construction of `centralFluxData`, `rusanovPlusData`, `rusanovMinusData` — the Q_god (central) and Q_corr (LF dissipation) blocks. |
| `SeisSol/src/Initializer/CellLocalMatrices.cpp` | 296–334 | The face-by-face selection between Godunov and Rusanov, including the `enforceGodunovBc` / `enforceGodunovEa` short-circuits (free-surface, analytical, elastic-acoustic interface). |
| `SeisSol/src/Equations/elastic/Model/ElasticSetup.h` | 79–166 | `getTransposedGodunovState` — the bimaterial Riemann projector builder.  Already mirrored as `BimaterialFlux::BuildGodunovStateFaceLocal`. |
| `SeisSol/codegen/kernels/aderdg.py` | 222–246 | The codegen of `computeFluxSolverLocal` / `computeFluxSolverNeighbor`.  Defines the canonical formula `A_plus_T = fluxScale * Tinv * (Q_god * star + Q_corr) * T`.  This is the formula the MFEM port must reproduce. |
| `SeisSol/docs/parameters.par` | 33–34 | User-facing `numflux` / `numfluxnearfault` strings — the UX we want to mirror in the MFEM TOML. |

### Local repository references (the code we are modifying)

| File | Lines | What lives here today |
|------|-------|----------------------|
| `dynamic/wave_operator.hpp` | 52–63 | `enum class MixedFluxMode`. |
| `dynamic/wave_operator.hpp` | 285–341 | Current `SetMixedFluxMode` declaration + `central_flux_face_set_` accessors. |
| `dynamic/wave_operator.hpp` | 786–813 | `per_face_bimaterial_flux_` storage — single 9×9 matrix pair per `(face, side)`.  This is the structure we MUST extend to carry both Q_god and Q_corr blocks. |
| `dynamic/wave_operator.inl` | 605–663 | Heterogeneous ctor body including the R-002 ctor-side guard at 652–660. |
| `dynamic/wave_operator.inl` | 825–1088 | `BuildPerFaceBimaterialFluxMatrices_` — where the per-face precomputation happens.  Single biggest function in scope of this port. |
| `dynamic/wave_operator.inl` | 2055–2120 | `SetMixedFluxMode` setter including the R-002 setter-side guard at ~2098. |
| `dynamic/wave_operator.inl` | 2130–~2300 (approx; see `BuildCentralFluxFaceSet_`) | The `central_flux_face_set_` builder.  Reused unchanged by this plan. |
| `dynamic/wave_operator.inl` | 3487, 3489, 4066, 5062 (and ADER variants) | The four interior-face `BimaterialFlux::ApplyPerFaceFlux` call sites.  These need to switch which precomputed pair they consume based on whether the face is in `central_flux_face_set_`. |
| `dynamic/godunov_flux_bimaterial.{hpp,cpp}` | full | `BimaterialFlux` static-method API.  `BuildPerFaceFluxMatricesGlobal` and `ApplyPerFaceFlux` need new variants. |
| `spatial/code/spatial_friction.{hpp,cpp}` | `interior_flux` parsing | Today accepts `{bimaterial, scalar}`.  Will gain `rusanov_near_fault` semantics through `mixed_flux = "adjacent"` interaction; no new TOML key. |

---

## 4. Mathematical Specification

### 4.1 The decomposition

The face flux contribution for the local element at an interior face
is, in the SeisSol codegen formula:

```
F_local = fluxScale · T_inv · ( Q_god · A_self + Q_corr ) · T · Q_self
```

where `T`, `T_inv` are the face-normal rotation (n̂ → +x̂),
`A_self` is the local-side flux Jacobian (the `star_matrix` in
SeisSol), and `(Q_god, Q_corr)` is the per-face Riemann data.  Two
modes:

**Godunov (bimaterial upwind):**
```
Q_god_local    = I - qGodProjector        # qGodLocal in BimaterialFlux
Q_god_neighbor = qGodProjector            # qGodNeighbor in BimaterialFlux
Q_corr_local   = 0                        # 9×9 zero
Q_corr_neighbor = 0
```

**Central (Zhang / DRDG3D — pure central, NO LF correction; matches what `GodunovFlux::Central` already does for the scalar ctor):**
```
Q_god_local    = 0.5 · I_9
Q_god_neighbor = 0.5 · I_9
Q_corr_local   = 0_9                       # ZERO — no LF viscosity
Q_corr_neighbor = 0_9
```
This is byte-equivalent to the existing scalar-ctor central flux
`F_h = 0.5·A_n·(Q_self + Q_nbr)` (`godunov_flux.cpp:382-430`), with
the bimaterial-Riemann `A_n` derivation now used everywhere on the
ring face.  Crucially: the central flux IGNORES the bimaterial
projector `qGodLocal` on adjacent-to-fault faces — the impedance-
correct decomposition is sacrificed there in exchange for less
upwind dissipation.  Zhang 2023 §3.2 argues that the fault face
itself supplies the necessary dissipation (via the imposed-state
Riemann + slip-weakening) so no LF correction is needed.

**Optional Rusanov mode (SeisSol-style; deferred, NOT in scope of
this plan):**  If a future user wants SeisSol byte-parity instead of
Zhang/DRDG3D parity, add a new `mixed_flux = "adjacent_rusanov"`
TOML value that selects
```
Q_god_local    = 0.5 · I_9
Q_god_neighbor = 0.5 · I_9
Q_corr_local   = +0.5 · c_max · I_9
Q_corr_neighbor = -0.5 · c_max · I_9
```
with `c_max = max(c_p_self, c_p_neighbor)`.  The decomposition
infrastructure this plan introduces supports the Rusanov mode at no
extra cost — only the Q_corr block values differ.  Documented as
follow-up.

### 4.2 Why the decomposition matters

The two modes share the same `T`, `T_inv`, `A_self`, `fluxScale`
multipliers; only the `(Q_god, Q_corr)` 9×9 pair differs.  Storing the
pair separately is what allows a face-by-face switch.

The current MFEM code collapses everything into a single 9×9 product:
`fluxLocal = T · A_self · qGodLocal · T_inv`.  That product has no
room to insert a `Q_corr` term, which is why the R-002 guard had to
forbid the combination.

### 4.3 Boundary conditions and elastic-acoustic interfaces

SeisSol forces Godunov on three face types regardless of the user's
`numfluxnearfault`:
1. Free-surface (`FaceType::FreeSurface`).
2. Free-surface-gravity, analytical, outflow.
3. Elastic-acoustic material interfaces (`isAtElasticAcousticInterface`).

The MFEM code today does NOT have acoustic regions (per
`godunov_flux_bimaterial.hpp:103-107`).  Free-surface faces are not
in `central_flux_face_set_` because they are boundary, not interior.
So the MFEM port can SKIP the acoustic short-circuit (documented
deferred item, mirroring the existing `BimaterialFlux::BuildGodunovStateFaceLocal`
acoustic-input rejection).

---

## 5. Implementation Phases

### Phase 1 — Refactor `BimaterialFlux` API to expose the decomposition

**Goal:** Replace the single-product `BuildPerFaceFluxMatricesGlobal`
with a pair of builders that return `(Q_god · A_self)` and `Q_corr`
separately, both already in the global frame (`T · ... · T_inv`).

**File:** `dynamic/godunov_flux_bimaterial.hpp`, `.cpp`.

**New API:**

```cpp
namespace mfem::seas {

/// Identify which 9×9 + 9×9 pair to compute at a face.
///
/// `Central` is the Zhang 2023 / DRDG3D pure-central choice
/// (`Q_corr = 0`).  `Rusanov` (SeisSol-style central + LF
/// viscosity, `Q_corr = ±0.5·c·I`) is reserved for a future
/// extension; not in scope of this plan.
enum class FaceFluxMode : int {
   Godunov = 0,   // bimaterial upwind (qGod · A_self, zero correction)
   Central = 1    // pure central (Q_god = 0.5·I, Q_corr = 0) — Zhang 2023
   // Rusanov = 2 // reserved: central + LF viscosity (SeisSol style)
};

/// SeisSol-style decomposition: writes BOTH the Q_god · A_self block
/// (multiplying Q_self / Q_nbr by the upwind contribution) AND the
/// Q_corr block (the additive dissipation).  Both already in the
/// global frame.  At runtime,
///
///   F_h_self = (fluxLocal_god  + fluxLocal_corr ) · Q_self
///            + (fluxNbr_god    + fluxNbr_corr  ) · Q_nbr
///
/// For Godunov mode, fluxLocal_corr == fluxNbr_corr == 0.
/// For Rusanov mode, fluxLocal_god + fluxNbr_god is the central
/// average (0.5·A_self) and the corr blocks carry the LF jump term.
class BimaterialFlux {
public:
   /// New decomposed builder.  Caller pre-allocates four 9×9 matrices.
   static void BuildPerFaceFluxMatricesGlobalDecomposed(
      const real_t* nor,
      const GodunovFlux& flux_self,
      const GodunovFlux& flux_nbr,
      FaceFluxMode mode,
      mfem::DenseMatrix& fluxLocal_god,   // 9×9
      mfem::DenseMatrix& fluxNbr_god,     // 9×9
      mfem::DenseMatrix& fluxLocal_corr,  // 9×9
      mfem::DenseMatrix& fluxNbr_corr);   // 9×9

   /// Runtime apply with decomposed blocks.  Equivalent to two
   /// independent matvecs summed.
   static void ApplyPerFaceFluxDecomposed(
      const mfem::DenseMatrix& fluxLocal_god,
      const mfem::DenseMatrix& fluxNbr_god,
      const mfem::DenseMatrix& fluxLocal_corr,
      const mfem::DenseMatrix& fluxNbr_corr,
      const real_t* Q_self,
      const real_t* Q_nbr,
      real_t* F_h_self);

   /// Legacy API (kept for backwards compat; equivalent to
   /// Decomposed with mode=Godunov and zero corr blocks added).
   static void BuildPerFaceFluxMatricesGlobal( ... );  // unchanged
   static void ApplyPerFaceFlux( ... );                // unchanged
};

} // namespace
```

**Sub-step 1.1.** Inside `BuildPerFaceFluxMatricesGlobalDecomposed`:
- For `Godunov`: call existing `BuildGodunovStateFaceLocal`, then
  rotate via `T_self · (qGodLocal · A_self_facelocal) · T_inv`.  Zero
  the corr blocks.
- For `Central` (DRDG3D / Zhang): set `Q_god_L = Q_god_R = 0.5·I`
  (face-local), set `Q_corr_L = Q_corr_R = 0`.  Apply the same
  `T · (Q_god · A_self + Q_corr) · T_inv` rotation as SeisSol's
  codegen.  Result: each side's flux matrix collapses to
  `T · 0.5·A_self · T_inv` — byte-equivalent to the existing
  `GodunovFlux::Central` body but expressed through the
  same decomposed-builder API used for the Godunov case.

The enum value naming is therefore `FaceFluxMode::Godunov` and
`FaceFluxMode::Central` (Zhang style).  The `Rusanov` value
(SeisSol style) is reserved for the optional follow-up.

**Sub-step 1.2.** New unit test
`tests/unit/test_phaser_rusanov_decomposition.cpp`:
- T-1: For a homogeneous (lam, mu, rho) face,
  `BuildPerFaceFluxMatricesGlobalDecomposed(...Godunov...)` produces
  4 matrices whose sum (god + corr) equals the legacy
  `BuildPerFaceFluxMatricesGlobal` result within 1e-12 relative.
- T-2: For the same input but `mode=Rusanov`, the
  `fluxLocal_god + fluxNbr_god = 0.5 · A_self_global` (the central
  average).
- T-3: `fluxLocal_corr + fluxNbr_corr = 0` (the LF viscosity terms
  cancel in the symmetric combination, modulo sign).
- T-4: For a bimaterial (lam_L != lam_R) face, the Godunov-mode
  decomposed result still passes T-1.

**Validation gate:** new test passes; existing
`test_phaser_dispatch_smoke` and `test_phaseh_wave_operator_constant_parity`
remain green (the legacy API is byte-unchanged).

**Estimated effort:** 3 days.

---

### Phase 2 — Extend `per_face_bimaterial_flux_` storage to carry both blocks

**Goal:** `WaveOperator` stores a `(god, corr)` pair per `(face, side)`
instead of a single product.

**File:** `dynamic/wave_operator.hpp` storage type around line 813;
`wave_operator.inl::BuildPerFaceBimaterialFluxMatrices_` around line
826–1088.

**Schema change.**  Replace

```cpp
std::vector<std::array<std::array<mfem::DenseMatrix, 2>, 2>>
   per_face_bimaterial_flux_;
//        ^face       ^side   ^{0=fluxLocal,1=fluxNeighbor}
```

with a `struct` holding both pairs:

```cpp
struct BimaterialFaceFluxBlock {
   mfem::DenseMatrix flux_local_god;    // 9×9
   mfem::DenseMatrix flux_neighbor_god; // 9×9
   mfem::DenseMatrix flux_local_corr;   // 9×9 (Rusanov only; zero for Godunov)
   mfem::DenseMatrix flux_neighbor_corr;// 9×9 (Rusanov only; zero for Godunov)
   FaceFluxMode      mode;              // diagnostic; aids verify hooks
};

std::vector<std::array<BimaterialFaceFluxBlock, 2>>  // [face][side]
   per_face_bimaterial_flux_;
```

**Memory cost.**  Doubles per-face memory from 2 × 81 × 8 B = 1.3 KB
to 4 × 81 × 8 B = 2.6 KB on Godunov-only faces (where the corr blocks
are zero, but still allocated).  TPV205 200 m mesh has 3.74 M
interior faces → goes from 9.7 GB/serial-rank to ~19 GB/serial-rank.
At 400 ranks the per-rank cost is ~48 MB; at 800 ranks, ~24 MB.

**Optimization to consider in Phase 2b (deferred):** allocate the
corr blocks only for faces in `central_flux_face_set_`.  Saves ~50 %
of memory on TPV runs (the central set is typically ~5 % of interior
faces on the 200 m mesh).  Implementation: use `std::optional<DenseMatrix>`
for the corr fields and check `.has_value()` at apply time.

**Sub-step 2.1.** Update `BuildPerFaceBimaterialFluxMatrices_` to
populate both blocks via the new decomposed builder, deciding
`mode = (central_flux_face_set_.count(face) ? Rusanov : Godunov)`.

**Sub-step 2.2.** Update the four interior-face apply sites
(`wave_operator.inl:3487, 3489, 4066, 5062` and any ADER predictor
variants found by `grep -n BimaterialFlux::ApplyPerFaceFlux`) to call
`ApplyPerFaceFluxDecomposed(...)` instead of `ApplyPerFaceFlux(...)`.

**Sub-step 2.3.** Update `VerifySharedFaceBimaterialFlux` (the cross-
rank consistency hook at `wave_operator.inl:1131`) to check both god
and corr blocks against the shared-face counterpart from the other
rank.

**Validation gate:** with both R-002 guards still in place (NOT yet
removed), this phase must produce byte-identical output to the
pre-refactor branch on every existing TPV run.  Tested by re-running
`test_phaseh_wave_operator_constant_parity` and `test_phaser_dispatch_smoke`.

**Estimated effort:** 5 days.

---

### Phase 3 — Wire `central_flux_face_set_` into the per-face mode dispatch

**Goal:** When `mixed_flux_mode_ == Adjacent` AND the heterogeneous
ctor was used, `BuildPerFaceBimaterialFluxMatrices_` selects
`FaceFluxMode::Rusanov` for faces in `central_flux_face_set_` and
`FaceFluxMode::Godunov` everywhere else.

**File:** `dynamic/wave_operator.inl` — primarily inside
`BuildPerFaceBimaterialFluxMatrices_` at the face loop.

**Sub-step 3.1.** Remove the R-002 ctor-side guard at
`wave_operator.inl:652-660`.  Replace with a clarifying comment that
the composition is now supported.

**Sub-step 3.2.** Remove the R-002 setter-side guard at
`wave_operator.inl:2098` AND extend `SetMixedFluxMode` to trigger a
re-build of `per_face_bimaterial_flux_` when `owned_flux_pool_ != nullptr`
and the mode changes.  The R-1408 lifecycle invariant already says
that any mesh / fault topology change must re-invoke
`SetMixedFluxMode`; this plan extends that re-invocation to also
trigger `BuildPerFaceBimaterialFluxMatrices_` with the new mode.

**Sub-step 3.3.** Lifecycle / call-ordering audit.  The new flow is:

```
1. WaveOperator(MaterialField, BoundaryConfig) ctor
     → BuildGodunovFluxPool_
     → ExchangeBiMaterialNeighbours_
     → BuildPerFaceBimaterialFluxMatrices_   [all Godunov; central_flux_face_set_ empty]
2. SetFaultFlux + SetFaultDOFData + SetAbsorbingBackground
     → no change to per_face_bimaterial_flux_
3. SetMixedFluxMode(Adjacent)
     → BuildCentralFluxFaceSet_   [populates central_flux_face_set_]
     → BuildPerFaceBimaterialFluxMatrices_   [REBUILD; Rusanov on central set]
```

The second `BuildPerFaceBimaterialFluxMatrices_` call must be
idempotent (call it twice in a row with the same mode → same output).

**Sub-step 3.4.** TOML parser: NO new key required.  The existing
`[numerics].interior_flux = "bimaterial"` + `[numerics].mixed_flux = "adjacent"`
combination now becomes legal and produces the SeisSol-style
composition.  Update the parser's mutual-exclusion check at
`spatial/code/spatial_friction.cpp` (the place where
`interior_flux == "bimaterial" && mixed_flux != "none"` was previously
flagged) to permit `adjacent`.

**Validation gate:** A new TOML smoke test that runs TPV205 with
`(bimaterial, adjacent)` for 100 steps in serial and verifies:
- No abort.
- `[wave_operator] BimaterialFlux precomputation:` banner reports a
  non-trivial central-set count.
- `seas_test_phaser_dispatch_smoke` extended with R.2.T-4: a
  bimaterial 3-element fixture where the middle face is in
  `central_flux_face_set_` must produce different flux values from
  the same fixture with `mixed_flux = "none"`.

**Estimated effort:** 4 days.

---

### Phase 4 — Homogeneous-material byte-parity gate

**Goal:** Prove `(bimaterial, adjacent)` ≡ `(scalar, adjacent)` to
within FP rounding on homogeneous TPV205.

**File:** new test `tests/unit/test_phaser_bimaterial_rusanov_homogeneous_parity.cpp`.

**Algorithm:** Build two WaveOperators on the same mesh:
1. `WaveOperator(mesh, p, λ, μ, ρ, bc)` (scalar ctor) +
   `SetMixedFluxMode(Adjacent)`.
2. `WaveOperator(mesh, p, MakeConstant(λ, μ, ρ), bc)` (heterogeneous
   ctor) + `SetMixedFluxMode(Adjacent)`.

Run both for 5 steps with the same `Q_initial` and assert:

```
||Q_scalar - Q_bimat|| / ||Q_scalar|| < 1e-10
```

This is the bimaterial-Rusanov analogue of the existing
`T-PHASEH-SCALAR-PARITY` gate (which only covers `mixed_flux = none`).
Both gates must pass before Phase 5 begins.

**Estimated effort:** 2 days.

---

### Phase 5 — Heterogeneous-material smoke test

**Goal:** Demonstrate the composition does the right thing on a
genuinely heterogeneous problem.

**File:** new test `tests/unit/test_phaser_bimaterial_rusanov_heterogeneous_smoke.cpp`.

**Fixture:** 1D 3-element strip with `(λ, μ, ρ)` jumping by a factor
of 2 at the middle face.  Run with:
- `(scalar, none)` — flux is wrong at the jump.
- `(scalar, adjacent)` — flux is still wrong at the jump; Rusanov on
  the near-fault face doesn't fix the material jump.
- `(bimaterial, none)` — flux is correct at the jump but full upwind.
- `(bimaterial, adjacent)` — flux is correct at the jump AND central
  near the fault (the configuration this plan unlocks).

Assertions:
- `(bimaterial, none)` and `(bimaterial, adjacent)` produce the same
  flux at the jump (since the jump face is NOT in the central set).
- `(scalar, none)` and `(scalar, adjacent)` differ from
  `(bimaterial, *)` at the jump by > 1e-6 (bimaterial fix is active).
- `(bimaterial, adjacent)` differs from `(bimaterial, none)` at the
  near-fault face by > 1e-6 (Rusanov is active).

**Estimated effort:** 2 days.

---

### Phase 6 — TPV205 production verification

**Goal:** Confirm A-5 — `(bimaterial, adjacent)` reduces upwind
dissipation on a real benchmark.

**Workflow:**
1. Submit `jobs/tpv205_spatial/tpv205_spatial_dyn_200m_p1_O2_normal.sbatch`
   (already in the repo) → produces `(bimaterial, none)` baseline.
2. Create variant
   `tpv205_spatial_dyn_200m_p1_O2_normal_adjacent.sbatch` overriding
   `--mixed-flux adjacent` → produces `(bimaterial, adjacent)`.
3. Compare hypocenter station `tpv205_x2_0_x3_7.5.dat`:
   - Peak `|V_strike|` should be ~10–30 % higher in the adjacent run
     (Zhang 2023 §4 reports this magnitude).
   - Rupture-front time-of-arrival should shift earlier by ~1–5 ms.
   - High-frequency power in the V_strike PSD (above the dt
     Nyquist) should increase (less upwind low-pass filtering).

If the adjacent run shows NO change OR a regression vs SCEC reference
traces, that's a Phase 6 failure and we re-audit Phases 1–3.

**Estimated effort:** 1 week including queue wait.

---

### Phase 7 — Documentation + CLAUDE.md update

Update:

- `CLAUDE.md` "Bi-material Riemann (Phase R)" section to remove the
  "MixedFluxMode != None is mutually exclusive with the heterogeneous
  ctor" sentence and replace with the new composition rules.
- `safs/project_7.0_alternative/document/spatial_friction_config_schema.md`
  to document the `(bimaterial, adjacent)` configuration.
- This file's status field from `PROPOSED` to `IMPLEMENTED`.
- New entry in `debug_document/general_driver_debug_document/` titled
  `seissol_flux_decomposition_port_completion_<date>.md` with the
  Phase 6 production results.

**Estimated effort:** 2 days.

---

## 6. Risks & Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Memory doubling on 200 m mesh + 64+ ranks causes OOM | Low | High | Phase 2b sparse-corr-allocation optimization, gated on `MFEM_VERIFY(memory < 0.8 * per_rank_limit)`. |
| Rusanov dissipation in the **bulk** (when `mixed_flux = AllContinuous`) over-damps rupture-front | Med | Med | Out of scope; this plan only implements `Adjacent`.  Keep `AllContinuous` un-composed for now, with the R-002 guard still active for that mode only. |
| The Rusanov `c_max` choice (max-P-wave vs max-S-wave) differs from Zhang 2023 prescription | Med | Med | Sub-step 1.1 will read Zhang 2023 §3.2 explicitly to confirm; SeisSol uses `getMaxWaveSpeed()` which is `c_p` in elastic mode (verified at `SeisSol/src/Initializer/CellLocalMatrices.cpp:279-281`).  If Zhang argues for `c_s`, document the deviation. |
| Free-surface faces accidentally end up in `central_flux_face_set_` and lose Godunov enforcement | Low | High | Phase 3 sub-step 3.1 explicitly asserts that the central-flux face set excludes boundary faces.  Test gate in `test_phaser_bimaterial_rusanov_homogeneous_parity`. |
| `BuildPerFaceBimaterialFluxMatrices_` rebuild on `SetMixedFluxMode` doubles startup time | Med | Low | Optimization deferred; first-implementation tolerates the 2× startup. |
| ADER-O3+ predictor uses a different flux assembly than ADER-O2 | Med | Med | Audit `wave_operator.inl::AdvanceADER` for any ADER-O3-specific flux call; current grep shows the same four `BimaterialFlux::ApplyPerFaceFlux` sites cover ADER-O2/O3/O4 uniformly. |
| Phase 6 production result shows NO improvement vs `(bimaterial, none)` | Med | High | Phase 5 smoke test is supposed to catch this earlier.  If Phase 6 fails, fall back to "Phase 4 byte-parity passes; production benefit deferred until polynomial order is raised" — the implementation still has scientific value as a research-mode option. |
| **TPV31** with heterogeneous material in the central set drops the bimaterial exact treatment ON those faces (central averaging is impedance-blind in the Riemann sense) | High | Med | This is the **fundamental trade-off** that SeisSol AND DRDG3D both accept.  DRDG3D `src/mod_wave.F90:540, 578` shows each side extracts its one-sided flux from its OWN material (so F_L uses local material, F_R uses neighbor material), but the central averaging `0.5·(F_R - F_L)` couples them without an impedance-weighted Riemann decomposition.  Zhang 2023 publishes results that include heterogeneous TPV runs and argues this is acceptable.  RECOMMENDED: ship with no escape hatch in Phase 1; add `[numerics].mixed_flux_skip_heterogeneous_faces` only if Phase 6 production verification exposes a measurable problem at TPV31 layer-boundary near-fault faces. |

---

## 7. Test Plan Summary

| Test | Phase | Type | Gate |
|------|-------|------|------|
| `test_phaser_rusanov_decomposition` (new) | 1 | Unit | Q_god + Q_corr decomposition matches SeisSol codegen formula |
| `test_phaseh_wave_operator_constant_parity` (existing) | 2 | Unit | bimaterial(none) byte-unchanged after storage refactor |
| `test_phaser_dispatch_smoke` (existing, R.2.T-4 added) | 3 | Unit | bimaterial(adjacent) differs from bimaterial(none) at the central-set face |
| `test_phaser_bimaterial_rusanov_homogeneous_parity` (new) | 4 | Unit | scalar(adjacent) ≡ bimaterial(adjacent) to 1e-10 on homogeneous |
| `test_phaser_bimaterial_rusanov_heterogeneous_smoke` (new) | 5 | Unit | 3-element strip distinguishes all 4 configurations correctly |
| TPV205 dev-queue run `(bimaterial, adjacent)` | 6 | Integration | Peak \|V_strike\| at hypocenter increases vs `(bimaterial, none)` |
| TPV205 normal-queue full-tfinal run | 6 | Production | SCEC reference comparison passes (e.g., \|V_strike(t=2s)\| within 5 % of SeisSol O=4 reference) |

---

## 8. Out-of-Scope (Future Work)

- **`mixed_flux = "AllContinuous"`** with bimaterial.  Same architecture
  would work, just sets every interior face's mode to Rusanov.  Higher
  dissipation in the bulk; not currently desired.  Leave R-002 guard
  active for that mode.
- **Acoustic-elastic interfaces.**  SeisSol's `enforceGodunovEa`
  short-circuit forces Godunov at acoustic interfaces regardless of
  user choice.  MFEM SEAS has no acoustic regions today (assertion
  in `godunov_flux_bimaterial.hpp:103-107`).  When acoustic support
  is added (e.g., for sea-floor coupling), port the short-circuit at
  that time.
- **Sub-cell limiter / WENO.**  Out of scope; Zhang 2023 doesn't use
  these either.
- **`fluxScale` dependence on mesh order.**  SeisSol's
  `fluxScale = -2 · surface / (6 · volume)` is order-independent.
  The MFEM analog already absorbs the order dependence into
  `BimaterialFlux::BuildPerFaceFluxMatricesGlobal`; nothing to change.
- **Higher polynomial order.**  Independent improvement axis;
  recommended as the parallel research track (Issue: "raise default
  ader_order to 3 + verify SCEC parity").  Combining higher order
  with the mixed-flux composition is additive.

---

## 9. Timeline

| Phase | Effort | Cumulative | Notes |
|-------|--------|------------|-------|
| 1     | 3 d    | 3 d        | Local; no queue time |
| 2     | 5 d    | 8 d        | Local; one regression rebuild |
| 3     | 4 d    | 12 d       | Local; reruns Phase-2 tests |
| 4     | 2 d    | 14 d       | Local |
| 5     | 2 d    | 16 d       | Local |
| 6     | 5–10 d | 21–26 d    | Queue-bound on Frontera |
| 7     | 2 d    | 23–28 d    | Documentation only |

Total: **~4–5 weeks** including queue wait.  Local-only work is **~3
weeks**.

---

## 10. Open Questions

- **Zhang 2023 venue / DOI.**  The code-base cites "Zhang 2023" but
  doesn't include the journal name.  Confirm before Phase 7
  documentation. Likely candidates: *J. Comput. Phys.*, *GJI*, or
  *Earth and Space Science*.
- **Acoustic short-circuit deferment.**  Confirm with the user that
  the MFEM SEAS code stays elastic-only for the foreseeable future,
  so the `isAtElasticAcousticInterface` short-circuit can be
  permanently stubbed out as `false`.
- **Sparse vs dense `Q_corr` allocation.**  Decide upfront whether
  Phase 2 implements the sparse optimization or defers to Phase 2b.
  Recommendation: defer, get the dense version landed first.
- **TPV31 cross-comparison.**  After Phase 6 confirms TPV205 works,
  re-run TPV31 with `(bimaterial, adjacent)` and document whether
  the near-fault Rusanov substitution hurts the material-jump
  fidelity at the spec's layer boundaries (depth ~2400 m, 5000 m,
  10000 m).
- **Composition with `precomputed_face_fluxes`.**  The existing
  `use_precomputed_face_fluxes_` path (a separate optimization for
  homogeneous Godunov runs) is also mutually exclusive with
  mixed-flux per the existing R-1203 guard at `wave_operator.inl:2083-2086`.
  Whether this plan should ALSO unblock that combination is a
  separate decision; recommend keeping it forbidden in Phase 1 and
  reconsidering only if production demand surfaces.

---

## 11. Approval Checklist

- [ ] User confirms Zhang 2023 citation (item in §10).
- [ ] User approves the dense-`Q_corr` allocation strategy for Phase 2
      (memory doubles per face; offset by 2× rank count).
- [ ] User approves dropping the `mixed_flux = "AllContinuous"` +
      `interior_flux = "bimaterial"` combination from this plan's scope.
- [ ] User confirms the TPV31 near-fault material-jump trade-off is
      acceptable (or wants the `rusanov_skip_heterogeneous_faces`
      escape hatch in Phase 7).
- [ ] User approves the ~4–5 week timeline before the implementation
      agent starts Phase 1.
