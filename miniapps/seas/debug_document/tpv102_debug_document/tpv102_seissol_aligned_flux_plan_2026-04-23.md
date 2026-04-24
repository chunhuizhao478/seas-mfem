# TPV102 — Topology-Based Precomputed Face-Rotation Plan (SeisSol-inspired; revised stepwise, post-review)

Date: 2026-04-23
Status: **design / review — no code written yet**
Scope: `miniapps/seas/dynamic/` only. MFEM library proper is **not** touched. BP5 quasi-dynamic runs are **not** touched.

Supersedes the prior "SeisSol-Aligned" draft of this document as well as Fix A and Fix B in
[tpv102_unit_tests_to_narrow_bug_2026-04-23.md](./tpv102_unit_tests_to_narrow_bug_2026-04-23.md).

## 0⁗. Fifth-round review response map (post-2026-04-23 call-site-scope pass)

Round-5 adversarial review cross-checked the §6.2 dispatch against the *scopes* at which `ComputeADERFaceFluxRHS`'s pre-existing locals are declared, and the temporal ordering between `WaveOperator` ctor and `SetFaultDOFData` for `fault_interior_faces_` / `fault_shared_faces_`. Found two CRITICAL items.

| round-5 finding | severity | status | where addressed |
|---|---|---|---|
| R5-001 §6.2 pre-existing-locals block lists `dof_offset2`, `shape2`, `I_nbr`, `is_fault` as outer-scope, but the existing `ComputeADERFaceFluxRHS` declares them only INSIDE the `else (interior)` branch (line 1870 onward); dispatch restructure would hit undefined-variable compile errors or undefined behaviour on boundary faces | CRITICAL (compile) | **FIXED in §6.2 — new §6.2.hoist subsection** | explicit hoisting spec added; conditional on `!is_boundary` so boundary faces never compute `e2 * ndof_per_el_` with `e2 < 0` |
| R5-002 R4-002's ctor-time `fault_face_set_` population is brittle: the ctor anchor is a specific line range (260–283) and a fix-agent could misplace the code, producing an empty set; also if future refactors move the fault-face array population, the ctor-time fill silently breaks | CRITICAL (correctness, latent) | **FIXED in §6.3a** | population moved from ctor into `UsePrecomputedFaceFluxes(true)` body, just before `Init(...)`, where `fault_interior_faces_` / `fault_shared_faces_` are guaranteed populated by the ctor's lines 260–283 |

**Correction to R5-002's framing** (reviewer's verbatim claim is not fully accurate, but the suggested fix is better engineering regardless): the existing `WaveOperator` ctor DOES populate `fault_interior_faces_` (at `wave_operator.inl:269`, `fault_interior_faces_.Append(face_idx);`) and `fault_shared_faces_` (at `wave_operator.inl:283`, `fault_shared_faces_.Append(sf);`). A ctor-time `fault_face_set_` population inserted after line 283 would therefore see populated arrays. Nonetheless we adopt the reviewer's proposed fix (late population in `UsePrecomputedFaceFluxes(true)`) because:

1. **Robust against ctor-ordering drift**: if a future refactor moves fault-face array population to a later helper or to `SetFaultDOFData`, the late-population approach works unchanged. The ctor-anchor approach silently breaks.
2. **Fix-agent placement-error resistance**: a mechanical fix-agent could insert the R4-002 ctor loop at the wrong line (e.g. right after the `fault_interior_faces_.SetSize(0)` at line 260, before the `Append` loop populates it). Late population has a single unambiguous call site.
3. **Diagnosable failure mode**: an `MFEM_VERIFY` at the late-population site catches the case where a driver calls `UsePrecomputedFaceFluxes(true)` with empty fault-face arrays (indicating either a driver-ordering bug OR a codebase refactor that lost the population step).

## 0‴. Fourth-round review response map (post-2026-04-23 call-site pass)

Round-4 adversarial review cross-checked every call site in the plan against the R3-revised function signatures and member declarations. Found two CRITICAL arity/initialisation bugs and one maintainer-trap. All addressed below.

| round-4 finding | severity | status | where addressed |
|---|---|---|---|
| R4-001 `§6.3 UsePrecomputedFaceFluxes` calls `Init(...)` with 6 args but §5.1.2 (post-R3) declares 8 required params | CRITICAL (compile) | **FIXED in §6.3** | call site expanded to pass `face_bdr_attr_` and `shared_mesh_face_set_` (both pre-existing `WaveOperator` members) |
| R4-002 `fault_face_set_` declared as a new member but never populated → Init would index fault faces as interior | CRITICAL (correctness) | **FIXED in §6.3 + new §6.3a** | ctor-level population spec added: `fault_face_set_` is built eagerly from the existing `fault_interior_faces_` / `fault_shared_faces_` arrays right after those are constructed, independent of the opt-in flag |
| R4-003 `BuildBoundaryMatrices{Godunov,Gamma,Absorbing}(n, t1, t2, flux, nApNm1)` take `(t1, t2)` but never use them (probe-path calls `FreeSurfaceGodunovTotal(n, ...)` which internally runs `BuildFrame(n)`) | MODERATE (maintainer trap) | **FIXED in §5.2.1 + §5.2.2 + §5.2.3** | dead `(t1, t2)` parameters removed from the boundary-build signatures; rotation-invariance argument for isotropic elastic BCs documented; `FaceEntry::tangent1 / tangent2` flagged as diagnostic-only for boundary entries |

**Rotation-invariance argument for R4-003** (load-bearing for the correctness of the boundary probe): for isotropic elastic free-surface BCs (both `FreeSurfaceGodunovTotal` and `FreeSurfaceTotal`) with `Q_bg = 0`, the BC flux is **invariant under in-plane rotations of `(t1, t2)` around the shared `n`**:

- `FreeSurfaceTotal` (γ-mirror): flips the rotated-frame components `σ_nn`, `σ_nt1`, `σ_nt2`. Under an in-plane rotation `θ` around `n`, `σ_nn` is invariant (scalar), and `(σ_nt1, σ_nt2)` transforms as a 2-vector; flipping both components of a 2-vector gives the same result regardless of the in-plane rotation angle. So `F_h_γ(n, t1_bf, t2_bf) = F_h_γ(n, t1_topo, t2_topo)` bit-exactly.
- `FreeSurfaceGodunovTotal` (Godunov projection): sets `σ_nn = σ_nt1 = σ_nt2 = 0` in the rotated frame and adjusts velocities via `Z_p` / `Z_s`. Both tangential traction components being zero is an in-plane-rotation-invariant condition (the full tangent plane traction vanishes); velocity adjustments use only `σ_nn / Z_p` and `|σ_nt| / Z_s` which are also in-plane-rotation-invariant. Same conclusion.
- `AbsorbingTotal` (with `Q_bg = 0`): reduces to `A^+ · I_self` in the rotated frame; `A^+` acts on the same (σ_nn, σ_nt1, σ_nt2, v_n, v_t1, v_t2) components with coefficients depending only on `ρ, Z_p, Z_s`. In-plane rotation commutes with `A^+` for isotropic elastic, so `F_h_abs(n, t1_bf, t2_bf) = F_h_abs(n, t1_topo, t2_topo)` bit-exactly.

**The key orientation-ambiguity source** — MFEM's `CalcOrtho(J_F)` giving `+n` vs `−n` across orbit-equivalent faces — is **resolved by the probe's `n` argument** being topology-based (always outward from the owning cell), regardless of what `(t1, t2)` `BuildFrame(n)` picks internally. `(t1, t2)` in the boundary-build signatures was genuinely dead code; removing the parameters also removes the risk of a future maintainer "using" them inconsistently (e.g. calling `BuildRotation(n, t1_topo, t2_topo)` to rotate the probed matrix, which would double-rotate and corrupt the result).

**Anisotropic-material caveat** (out of scope): if this plan is ever extended to anisotropic elastic materials, the rotation-invariance argument breaks (the material stiffness tensor is not in-plane-rotation-invariant around `n`), and the probe approach becomes wrong. Step B's `P_BIMATERIAL_GUARD` (§5.2.4) already aborts on bimaterial faces; add an analogous `P_ANISOTROPIC_GUARD` in Phase 5 if/when anisotropic support is explored.

## 0″. Third-round review response map (post-2026-04-23 API-verification pass)

Round-3 adversarial review cross-checked every MFEM API call in the plan against `mfem/fem/geom.{hpp,cpp}`, `mfem/mesh/mesh.hpp`, and the existing `wave_operator.inl`. Found five CRITICAL items (four API hallucinations + one stale identifier) and four moderate spec gaps. All addressed below.

| round-3 finding | severity | status | where addressed |
|---|---|---|---|
| R3-001 `Geometry::GetFaceVertices(Geometry::TETRAHEDRON, ...)` is not an MFEM API | CRITICAL | **FIXED in §5.1.2 + §5.1.3 + §14** | replaced with `Geometry::Constants<Geometry::TETRAHEDRON>::FaceVert[4][3]`; pre-filled verified values `{{1,2,3},{0,3,2},{0,1,3},{0,2,1}}`; AssertMFEMFace2NodesTable checks against the MFEM constant |
| R3-002 `mesh.FaceIsBoundary(...)` is not an MFEM API | CRITICAL | **FIXED in §5.3** | replaced with the same `(Elem2No < 0) && (bdr_attr > 0) && !is_shared` predicate used by `wave_operator.inl:874` |
| R3-003 `mesh.GetBdrElementIndex(face_idx)` is not an MFEM API | CRITICAL | **FIXED in §5.1.2 + §5.3** | Init signature extended to take the existing `face_bdr_attr_` table and `shared_mesh_face_set_` set; no new API calls |
| R3-004 `ftr->ElemNNo` is a typo | CRITICAL | **FIXED in §5.3 step 2e** | rewritten as `(ftr->Elem1No == e) ? ftr->Elem2No : ftr->Elem1No` |
| R3-005 §5.3 step 3 references stale `face_side_to_entry_` / `(face_idx, local_side)` key | CRITICAL | **FIXED in §5.3 step 3** | updated to `face_elem_to_entry_` with `(face_idx, elem)` encoding |
| R3-006 sign-extension bug in the encoded key for negative elem ids | MODERATE | **FIXED in §5.1.2** | encoding cast through `uint32_t` |
| R3-007 Phase 1 vs 2b shared-face scope confusion in Init | MODERATE | **FIXED in §5.3 step 2b + §7.2** | explicit skip of shared faces in Phase 1; Phase 2b specifies one entry per shared face per rank |
| R3-008 Init takes `const Mesh&` but `GetFaceElementTransformations` is non-const | MODERATE | **FIXED in §5.1.2 signature + §5.3** | changed to `Mesh &` (non-const) matching existing `wave_operator.inl` convention |
| R3-009 p=1 isoparametric precondition not stated | MODERATE | **FIXED in §5.1.6** | explicit precondition + runtime check that face Jacobian is constant across QPs |

## 0′. Second-round review response map (post-2026-04-23 adversarial pass)

Round-2 adversarial review found five additional CRITICAL items and four moderate ones that slipped into the first-round fix. All addressed below, with pointers to the exact sections updated.

| round-2 finding | severity | status | where addressed |
|---|---|---|---|
| R2-001 P_SEISSOL_VALUE formula swaps A⁺/A⁻ on I_self/I_nbr (QgodLocal^T^T = I − P₊ routes through (I − P₊)·A_x = A_x⁻, mismatch with MFEM F_h = A⁺·I_self + A⁻·I_nbr) | CRITICAL | **FIXED in §5.1.7** | P_SEISSOL_VALUE rewritten to use direct eigendecomposition A⁺ = R·diag(Λ₊)·R⁻¹, dropping the layered QgodLocal formula |
| R2-002 `AddInteriorFaceRhs(caller_elem, ...)` cannot resolve its entry; lookup key is `(face_idx, local_side)` but function has no mesh access | CRITICAL | **FIXED in §5.1.2** | key rekeyed as `(face_idx, caller_elem_id)` via `face_elem_to_entry_`; no mesh member needed |
| R2-003 §6.2 dispatch references undeclared `elem1`/`elem2`; actual locals are `e1`/`e2` | CRITICAL | **FIXED in §6.2** | renamed to `e1`/`e2`; added pre-existing-locals block for the implementer |
| R2-004 §5.1.5 narrative writes identity `Interior(n, L, R) = Interior(-n, R, L)` (wrong sign); test P_CONSERVATION wants the opposite sign | CRITICAL | **FIXED in §5.1.5** | identity corrected to `Interior(n, L, R) = -Interior(-n, R, L)`; conservation consequence spelled out |
| R2-005 §5.1.4 step 9 `cell_volume = |det([ev[1]-ev[0],...])| / 6` subtracts integer vertex indices instead of coordinates | CRITICAL | **FIXED in §5.1.4** | explicit coordinate fetch `p_elem[i] = mesh.GetVertex(ev[i])`; formula rewritten in terms of `p_elem` |
| R2-006 §5.3 boundary classification does not explicit fs_mode → Godunov/Gamma routing | MODERATE | **FIXED in §5.3** | explicit classification pseudocode block added |
| R2-007 §5.1.2 FACE2NODES_MFEM placeholder ships concrete values before reviewer sign-off | MODERATE | **FIXED in §5.1.2** | placeholder replaced with `{-1,-1,-1}` sentinel that aborts in `AssertMFEMFace2NodesTable()` |
| R2-008 §6.2 one-line "unchanged from today" abbreviates ~150 lines of fault bookkeeping | MODERATE | **FIXED in §6.2 + new test §5.4** | fault-preservation instruction explicit; new test `test_R008_fault_branch_preserved_after_dispatch_reorder` |
| R2-009 §5.2.1 probe-based construction assumes linearity without stating the precondition | MODERATE | **FIXED in §5.2.1** | precondition + `MFEM_VERIFY` + new test `test_R009_probe_linearity` |

## 0. Review-finding response map (R-001 — R-010, plus unreviewed-area items)

The title of this plan has been changed from "SeisSol-Aligned" to "Topology-Based Precomputed Face-Rotation (SeisSol-inspired)" per **R-002**/R-010. The revised plan delivers topology-based outward normals + precomputed per-cell per-side rotation matrices that carry MFEM's A⁺/A⁻ split — **not** SeisSol's `AplusT = fluxScale · Tinvᵀ · (QgodLocal·starMatrix + QcorrLocal) · Tᵀ` stored-flux design. The stored matrix layout differs from SeisSol's by construction; the runtime per-QP integration `rhs -= w · shape(i) · F_h` is preserved from MFEM. The equivalence criterion with SeisSol is **flux VALUE per `(n, I_self, I_nbr)`**, not the stored matrix layout.

| finding | category | status | where addressed |
|---|---|---|---|
| R-001 dispatch bug silences fault flux on opt-in | CRITICAL | **FIXED in §6.2** | new if/else ordering with `is_fault` first |
| R-002 "SeisSol-aligned" misnomer | CRITICAL | **FIXED in title + §1 + §2.5** | plan retitled; semantics clarified; new test P_SEISSOL_VALUE added |
| R-003 missing flux-conservation assertion for two-entry design | CRITICAL | **FIXED in §5.1.7 (P_CONSERVATION)** | new acceptance sub-test |
| R-004 tangent re-orthogonalization deviates from SeisSol | MODERATE | **FIXED in §5.1.4** | plan now matches SeisSol's construction: raw `t1 = v1 − v0`, do not re-orthogonalize |
| R-005 FACE2NODES_MFEM spec | MODERATE | **FIXED in §5.1.3** | specified via MFEM's `Geometry::GetFaceVertices` public API + self-check assertion |
| R-006 FreeSurfaceGamma not in SeisSol | MODERATE | **FIXED in §5.2.1/§5.2.3** | FreeSurfaceGodunov is mandatory; FreeSurfaceGamma is legacy-only |
| R-007 Absorbing semantic gap | MODERATE | **FIXED in §5.2.2** | MFEM-native semantics (A⁺ only) documented as the chosen path with explicit rationale; SeisSol Outflow flagged as distinct |
| R-008 Γ placement ambiguity | MODERATE | **FIXED in §5.2.1** | formula rewritten as `A = T · (A⁺ + A⁻ · Γ) · Tinv` |
| R-009 enum coverage | LOW | **FIXED in §5.1.2** | enum annotated with SeisSol mapping; delegated / out-of-scope types documented |
| R-010 "by coincidence" phrasing | LOW | **FIXED in §3.4** | phrasing replaced with mechanistic explanation |
| unreviewed: fluxScale derivation | — | **ADDRESSED in §5.1.6** | dimensional reconciliation pinned as a Step A checklist item |
| unreviewed: MPI shared-face handling | — | **ADDRESSED in new Phase 2b (§7)** | `ComputeADERSharedFaceFluxRHS` explicitly in scope; Phase 3 np≥2 gates tied to it |

**Per R-007**: Before implementation begins, the user-facing semantic choice is now set as follows (and the sign-off gate in §14 confirms):

- **Absorbing BC**: MFEM-native `F_h = A⁺ · I_self` (first-order Sommerfeld) is the chosen semantics. SeisSol `Outflow` (full `A · I_self` radiating both characteristics) is explicitly **not** implemented in this plan. Rationale: TPV102 uses free-surface + absorbing attributes only at the outer cube's x and z sides, and changing absorbing semantics would constitute a physics change that is out of scope for a bug-fix plan.
- **Free-surface BC**: **FreeSurfaceGodunov is the mandatory variant**; FreeSurfaceGamma is retained as a legacy pass-through so the runtime path is not removed.
- **Interior**: MFEM's A⁺/A⁻ split applied in the precomputed per-cell-per-side rotation. The per-face QP weight `w = ip.weight * nor_len` stays at the runtime call site (no fluxScale bake-in at Phase 1).

## 1. Overview

Replace the ADER non-fault face dispatch in `WaveOperator` with a **topology-based precomputed per-cell per-side rotation operator** that eliminates MFEM's runtime `nor = CalcOrtho(J_F)` orientation ambiguity. Each cell gets its own outward normal + tangents from vertex topology at init; per-cell per-side `nApNm1 / nAmNm1` matrices carrying MFEM's A⁺ / A⁻ split are precomputed and applied at runtime in the existing `rhs -= w · shape(i) · F_h` loop structure. Fault branch and RK path are untouched in the first patch. Default behaviour is unchanged; the opt-in is enabled only for the TPV102 driver after the audit gates go green.

This is **not** SeisSol's full-face-flux precomputation (no fluxScale bake-in, MFEM A± split instead of SeisSol `QgodLocal · starMatrix`, two entries per face instead of one). It delivers the same **numerical flux value per `(n, I_self, I_nbr)`** as MFEM's runtime `flux_.Interior` / `flux_.FreeSurfaceTotal` / `flux_.FreeSurfaceGodunovTotal` / `flux_.AbsorbingTotal`, but with deterministic topology-based normals that remove the orbit-asymmetry source identified in the debug audit.

## 2. Constraints

### 2.1 Files that MAY be modified (this plan)

- `dynamic/wave_operator.hpp`
- `dynamic/wave_operator.inl`
- `dynamic/precomputed_face_fluxes.hpp` (new)
- `dynamic/precomputed_face_fluxes.cpp` (new)
- `dynamic/godunov_flux.hpp` — **only** to add two additive const accessors (`GetAxPlus()`, `GetAxMinus()`). No behaviour change, no method removal, no signature change on existing public API.
- `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp` (add opt-in toggle for Phase 3 gates)
- `tests/unit/test_precomputed_fluxes_*.cpp` (new — per Phase 1 / 2 acceptance)
- `miniapps/seas/Makefile` (register new test binaries only; no link-graph changes)

### 2.2 Files that MUST NOT be touched

- `/Users/chunhuizhao/projects/mfem/**` — MFEM library proper.
- `miniapps/seas/domain/` (BP5's elasticity operator).
- `miniapps/seas/fault/` (BP5 / TPV102 friction).
- `miniapps/seas/solver/` (SEASOperator, BP5 coupling).
- `miniapps/seas/drivers/seas_driver.cpp` (BP5 driver).
- `miniapps/seas/config/bp5_params.hpp` and `miniapps/seas/bp5/**`.
- `miniapps/seas/dynamic/fault_face_flux.{hpp,cpp}` (fault branch is clean per the audit).
- `miniapps/seas/dynamic/godunov_flux.cpp` — **no body changes**. The existing implementations of `Interior`, `FreeSurfaceTotal`, `FreeSurfaceGodunovTotal`, `AbsorbingTotal`, `BuildFrame`, `BuildRotation*`, `ApplySplitFlux` stay as-is and remain the default runtime path when the opt-in is off.

### 2.3 Runtime behaviour

- The `WaveOperator` opt-in flag `use_precomputed_face_fluxes_` defaults to **false**. Existing behaviour is unchanged for every driver.
- The flag is flipped to **true** only by the TPV102 driver, only after Phase 3 passes on all required fixtures.
- BP5 does not link `WaveOperator` (confirmed by `grep`), so no BP5 changes regardless.

### 2.4 Numerical constraints

- Precomputed matrices must reproduce the runtime `GodunovFlux::Interior(n_outward, I_self, I_nbr)` result to double-precision ULP per-cell, where `n_outward` is the **per-cell** outward normal built from topology. Per **R-010**: when MFEM's `CalcOrtho(J_F)` happens to produce the per-cell outward normal (this occurs deterministically on approximately half of interior non-fault faces under the Kuhn split's Elem1/Elem2 convention, and on all boundary faces because MFEM assigns the interior adjacent cell as Elem1 so CalcOrtho points outward from the domain), the precomputed and runtime results agree to ULP without any flip. On the other half of interior faces, the precomputed result equals the runtime result under the identity `Interior(−n, R, L) = Interior(n, L, R)` + the MFEM A⁺/A⁻ eigenvector sign convention. This identity is verified as test **P_CONSERVATION** (§5.1.7).
- Precomputed boundary matrices must reproduce `flux_.FreeSurfaceTotal / FreeSurfaceGodunovTotal / AbsorbingTotal` to ULP on the per-cell outward normal.
- Homogeneous material guard: `PrecomputedFaceFluxes::Init` must abort if it encounters a bimaterial interior face (matches `fault_face_flux.cpp:101`'s `homog_ok` pattern).

### 2.5 Equivalence contract with SeisSol

- **What we claim to match**: the flux **VALUE** `F_h(n_outward, I_self, I_nbr)` computed by MFEM's existing runtime path for every non-fault face QP, when the runtime `nor` happens to equal the per-cell outward `n_outward`.
- **What we do NOT claim to match**: SeisSol's stored `AplusT` matrix layout, SeisSol's `fluxScale` bake-in, SeisSol's `QgodLocal · starMatrix` split structure, SeisSol's `FaceType::Outflow` semantics, or SeisSol's single-entry-per-face flux conservation identity. The two-entry-per-face design is enforced conservative via **P_CONSERVATION** (§5.1.7) — not by construction.
- **Verification of SeisSol-flux-value compatibility**: test **P_SEISSOL_VALUE** (§5.1.8) hand-computes `F_h_seissol = T · (QgodLocal · starMatrix) · Tinv · I_self + T · (QgodNeighbor · starMatrix) · Tinv · I_nbr` for an interior face using SeisSol's `ElasticSetup.h:87-140` formulas and asserts the precomputed path produces the same global-frame `F_h` to 1e-10. Any residual above 1e-10 means MFEM's eigenvector matrix `R` has a column-ordering or sign convention that differs from SeisSol's `matR`, and the basis-change map must be documented before Phase 2 can start.

## 3. Phase structure

- **Phase 1** — Precompute non-fault face operators (Steps A + B).
- **Phase 2a** — ADER opt-in dispatch switch for *local* non-fault faces (Step C-local).
- **Phase 2b** — ADER opt-in dispatch switch for *shared* non-fault faces (Step C-shared, MPI).
- **Phase 3** — Prove second-step non-fault symmetry restoration under the new path (primary gate).
- **Phase 4** — Expand to RK path if needed (Step D, conditional).
- **Phase 5** — Storage optimization + default-on.

## 4. Constraints recap (terminology)

- "SeisSol-aligned" (old phrasing) → **retired**. Use "topology-based precomputed rotation" when referring to this plan.
- "Precomputed face flux" is a shorthand for "precomputed per-cell per-side rotated A⁺/A⁻ matrix"; it does **not** include SeisSol's fluxScale bake-in, not yet.
- "Per-cell outward normal" = computed from mesh topology using the element's own vertex ordering (not shared between cells; orbit-mirror cells have opposite normals).

## 5. Phase 1 — Precompute non-fault face operators

Two independent steps. Both must complete and pass their per-step gate before Phase 2a begins.

### 5.1 Step A — Interior non-fault face operator

#### 5.1.1 Goal
After Step A: `PrecomputedFaceFluxes::AddInteriorFaceRhs(...)` produces a per-QP rhs contribution that matches `flux_.Interior(n_outward_for_this_cell, I_self, I_nbr, F_h)`-driven accumulation to double-precision ULP on the Cartesian Kuhn fixture for every interior non-fault face. Includes a flux-conservation assertion.

Only interior non-fault faces are implemented in this step. Boundary and fault are out of scope (fault is untouched; boundary is Step B).

#### 5.1.2 Class interface

```cpp
// dynamic/precomputed_face_fluxes.hpp
namespace mfem { namespace seas {

enum class FacePrecomputedBC {
   Interior            = 0,   // == SeisSol FaceType::Regular
   FreeSurfaceGodunov  = 1,   // == SeisSol FaceType::FreeSurface (mandatory for SeisSol parity)
   FreeSurfaceGamma    = 2,   // MFEM-native γ-mirror; NOT in SeisSol; legacy path only
   Absorbing           = 3    // MFEM-native first-order Sommerfeld (F_h = A^+ · I_self);
                              // NOT identical to SeisSol FaceType::Outflow (which radiates full A)
};
// SeisSol FaceType coverage map (per R-009):
//   FaceType::Regular             -> FacePrecomputedBC::Interior
//   FaceType::FreeSurface         -> FacePrecomputedBC::FreeSurfaceGodunov (mandatory)
//   FaceType::DynamicRupture      -> delegated to FaultFaceFlux (untouched)
//   FaceType::Outflow             -> NOT implemented here; MFEM's Absorbing is not
//                                    semantically equivalent.  Documented in §5.2.2.
//   FaceType::FreeSurfaceGravity  -> out of scope (TPV102 flat surface)
//   FaceType::Dirichlet           -> out of scope
//   FaceType::Periodic            -> out of scope
//   FaceType::Analytical          -> out of scope

class PrecomputedFaceFluxes
{
public:
   struct FaceEntry {
      int    element       = -1;
      int    local_side    = -1;
      int    face_idx      = -1;
      int    neighbor_elem = -1;             // -1 if boundary
      FacePrecomputedBC bc = FacePrecomputedBC::Interior;
      real_t surface_area  = 0.0;
      real_t cell_volume   = 0.0;
      real_t normal[3]     = {0.0, 0.0, 0.0};
      // R4-003: for INTERIOR entries, tangent1/tangent2 are
      // load-bearing — they are the topology-based frame used in
      // BuildInteriorMatrices to construct nApNm1 / nAmNm1, and
      // they must match MFEM's runtime T / Tinv construction at the
      // dispatch site.  For BOUNDARY entries, tangent1/tangent2 are
      // DIAGNOSTIC-ONLY — the probe-based BuildBoundaryMatrices*
      // captures F_h using GodunovFlux::BuildFrame(n)'s internal
      // tangents (which differ from the topology ones).  A future
      // maintainer MUST NOT apply BuildRotation(n, tangent1, tangent2)
      // to a boundary entry's nApNm1 — that would double-rotate and
      // corrupt the result.  Rotation-invariance of isotropic elastic
      // BCs around `n` makes the discrepancy harmless at Phase 1 (see
      // §0‴), but only if the matrix is used as-is.
      real_t tangent1[3]   = {0.0, 0.0, 0.0};
      real_t tangent2[3]   = {0.0, 0.0, 0.0};
      DenseMatrix nApNm1;                    // 9 x 9, GLOBAL frame
      DenseMatrix nAmNm1;                    // 9 x 9, GLOBAL frame; zero for boundary
   };

   // R3-008 FIX: `mesh` is passed non-const because MFEM's
   // `Mesh::GetFaceElementTransformations(int)` and
   // `ParMesh::GetSharedFaceTransformations(int)` are both non-const
   // (see mfem/mesh/mesh.hpp:1937, and the const_cast pattern in
   // wave_operator.inl:1435).  Matches the existing WaveOperator
   // convention.
   //
   // R3-003 FIX: `face_bdr_attr` is the same precomputed table that
   // WaveOperator already maintains in its ctor (see
   // wave_operator.hpp:377 `std::vector<int> face_bdr_attr_`).  Passing
   // it avoids inventing a non-existent `Mesh::GetBdrElementIndex`
   // API and keeps boundary classification O(1) per face.
   //
   // R3-002 FIX: `shared_mesh_face_set` is the same precomputed set
   // that WaveOperator maintains in its ctor (see
   // wave_operator.hpp:399 `std::set<int> shared_mesh_face_set_`).
   // Required so Init can correctly distinguish boundary vs shared
   // faces on a ParMesh.
   void Init(Mesh &mesh,
             const FiniteElementSpace &fes,
             const BoundaryConfig &bc,
             const GodunovFlux &flux,
             FreeSurfaceBCMode fs_mode,
             const std::set<int> &fault_face_set,
             const std::vector<int> &face_bdr_attr,        // R3-003
             const std::set<int> &shared_mesh_face_set);   // R3-002

   // caller_elem is the GLOBAL element id whose rhs is being assembled.
   // This function resolves `(face_idx, caller_elem)` → FaceEntry via
   // `face_elem_to_entry_` (see private section below).  It does NOT
   // need `local_side` at the call site; local_side is stored inside
   // the FaceEntry for debugging but not required by the lookup key.
   // (R2-002 fix.)
   void AddInteriorFaceRhs(int face_idx,
                           int caller_elem,
                           const real_t *I_self, const real_t *I_nbr,
                           real_t w,
                           const real_t *shape1,
                           int ndof, int dof_offset1,
                           int ndof_total,
                           Vector &rhs) const;

   void AddBoundaryFaceRhs(int face_idx,
                           const real_t *I_self,
                           const real_t *bulk_bg_scaled,
                           real_t w,
                           const real_t *shape1,
                           int ndof, int dof_offset1,
                           int ndof_total,
                           Vector &rhs) const;

   bool IsInitialized() const { return initialized_; }

   // Test-visibility accessor (used by P_MATRIX_IDENTITY, P_CONSERVATION,
   // P_SEISSOL_VALUE).  Returns the FaceEntry keyed by (face_idx,
   // caller_elem_id); aborts on miss.  Per R2-002: the lookup key is
   // (face_idx, elem_id), NOT (face_idx, local_side), so callers that
   // only have the element id can resolve without mesh access.
   const FaceEntry &GetEntry(int face_idx, int caller_elem_id) const;

private:
   // Topology-based frame (see §5.1.3 – §5.1.4).
   static void ComputeCellFaceFrame(const Mesh &mesh,
                                    int elem, int local_side,
                                    real_t normal[3],
                                    real_t tangent1[3], real_t tangent2[3],
                                    real_t &surface_area,
                                    real_t &cell_volume);

   // Interior matrix composition (see §5.1.5).
   void BuildInteriorMatrices(const real_t normal[3],
                              const real_t tangent1[3],
                              const real_t tangent2[3],
                              const GodunovFlux &flux,
                              DenseMatrix &nApNm1,
                              DenseMatrix &nAmNm1) const;

   // Boundary matrix composition (see §5.2).  R4-003 FIX: (tangent1,
   // tangent2) removed from the signature — dead parameters in the
   // probe-based construction.  For isotropic elastic free-surface BCs
   // with bulk_bg = 0, the BC flux is rotation-invariant around `n`
   // (see §0‴), so the probe captures the correct flux value using
   // GodunovFlux::BuildFrame(n)'s internal tangent choice.
   void BuildBoundaryMatrices(FacePrecomputedBC bc,
                              const real_t normal[3],
                              const GodunovFlux &flux,
                              real_t Zp, real_t Zs,
                              DenseMatrix &nApNm1) const;

   // Self-check asserted at Init: MFEM's
   // `Geometry::Constants<Geometry::TETRAHEDRON>::FaceVert[4][3]` matches
   // the static FACE2NODES_MFEM[4][3] we hard-code below, bit-exactly.
   // Guards against MFEM upgrades that change the table silently.
   // Implementation body:
   //   for (int s = 0; s < 4; s++)
   //     for (int j = 0; j < 3; j++)
   //       MFEM_VERIFY(FACE2NODES_MFEM[s][j] ==
   //                   mfem::Geometry::Constants<
   //                       mfem::Geometry::TETRAHEDRON>::FaceVert[s][j],
   //                   "MFEM tet face-vertex table drift at side "
   //                   << s << " vertex " << j);
   static void AssertMFEMFace2NodesTable();

   // Hard-coded MFEM reference-element face-to-vertex index table for
   // a tetrahedron.  Values are verified against
   // `mfem::Geometry::Constants<mfem::Geometry::TETRAHEDRON>::FaceVert`
   // (mfem/fem/geom.cpp:987–988) — a static `const int[4][3]` compile-
   // time constant that is MFEM's canonical source of truth for this
   // table.  R3-001 FIX: the plan previously called a non-existent
   // `Geometry::GetFaceVertices(...)` method; the correct public
   // symbol is the `Constants<...>::FaceVert` static array.
   //
   // `AssertMFEMFace2NodesTable()` asserts bit-identical equality at
   // every Init call — if MFEM ever changes the reference-tet table,
   // or if this local copy gets edited, the assertion catches it.
   static constexpr int FACE2NODES_MFEM[4][3] = {
      /* face 0 */ {1, 2, 3},
      /* face 1 */ {0, 3, 2},
      /* face 2 */ {0, 1, 3},
      /* face 3 */ {0, 2, 1}
   };

   std::vector<FaceEntry> entries_;

   // Key: (face_idx, caller_elem_id) → entry index.
   // Encoding (R3-006 FIX — cast through uint32_t to prevent
   // sign-extension from negative caller_elem_id corrupting the key):
   //   static long long EncodeKey(int face_idx, int caller_elem_id)
   //   {
   //      return (static_cast<long long>(face_idx) << 32) |
   //             static_cast<long long>(
   //                static_cast<uint32_t>(caller_elem_id));
   //   }
   // A naive `long long(face_idx) << 32 | caller_elem_id` would
   // sign-extend a negative caller_elem_id (e.g. `-1`) to
   // `0xFFFFFFFFFFFFFFFF` and corrupt the upper 32 bits via the
   // bitwise-or.  Phase 1 asserts `caller_elem_id >= 0` at every
   // insert site (no shared faces yet; local elem ids are
   // non-negative); Phase 2b relaxes that assertion because shared-
   // face `Elem2No` can be negative.  The uint32_t cast is defensive
   // for that case.
   //
   // An interior non-fault face has at most TWO adjacent cells, so
   // this map has 1 or 2 entries per face (one keyed by Elem1's id,
   // one by Elem2's id).  A boundary face has exactly ONE (keyed by
   // Elem1's id).  Per R2-002: keying by (face_idx, caller_elem_id)
   // avoids needing mesh access inside AddInteriorFaceRhs /
   // AddBoundaryFaceRhs to resolve local_side; `FaceEntry::local_side`
   // is stored for debugging/assertions but not required at the
   // lookup site.
   std::unordered_map<long long, int> face_elem_to_entry_;

   bool initialized_ = false;
};

}} // namespace mfem::seas
```

#### 5.1.3 MFEM face-to-vertex table (address R-005)

The hinge of the topology frame is the map from `local_side ∈ {0,1,2,3}` to the three vertex indices of that face. Spec (R3-001 FIX):

1. At Phase 1 Step A's init-time, read MFEM's reference-element face-vertex table from the **static compile-time constant** `mfem::Geometry::Constants<mfem::Geometry::TETRAHEDRON>::FaceVert` — a `const int[4][3]` exposed at `mfem/fem/geom.hpp:233` and defined at `mfem/fem/geom.cpp:987–988`. No runtime dispatch, no dummy-`Mesh` construction needed. (The plan previously referenced `Geometry::GetFaceVertices(...)`, which is **not an MFEM API**; R3-001 corrects this.)
2. `AssertMFEMFace2NodesTable()` asserts that the local `FACE2NODES_MFEM[4][3]` (§5.1.2) matches `Constants<TETRAHEDRON>::FaceVert[4][3]` bit-exactly for all four local sides. On mismatch, `MFEM_ABORT` with a message naming the local side and the two differing triples.
3. The hard-coded values in §5.1.2 are the verified-correct ones per `mfem/fem/geom.cpp:987–988`:
   ```
   face 0: {1, 2, 3}
   face 1: {0, 3, 2}
   face 2: {0, 1, 3}
   face 3: {0, 2, 1}
   ```
   These are not placeholders; the reviewer has already verified them against the MFEM source. The self-check guards against future MFEM upgrades that change the table.

**Acceptance test P_FACE2NODES_SELFCHECK** (§5.1.7) asserts `AssertMFEMFace2NodesTable()` passes on the live MFEM binary.

#### 5.1.4 Topology frame construction (address R-004)

`ComputeCellFaceFrame(mesh, elem, local_side, ...)` algorithm, mirroring SeisSol's `MeshTools::normalAndTangents` **exactly** (no re-orthogonalization):

1. Fetch `ev[0..3] = mesh.GetElementVertices(elem)` — 4 global **vertex indices** (integers).
2. Fetch the 4 element-vertex coordinates: `p_elem[i] = mesh.GetVertex(ev[i])` for `i = 0..3`, where each `p_elem[i]` is a `const real_t*` of length 3 (the vertex x/y/z coordinates).
3. Fetch the three face-vertex indices: `fn[0] = ev[FACE2NODES_MFEM[local_side][0]]`, `fn[1] = ev[FACE2NODES_MFEM[local_side][1]]`, `fn[2] = ev[FACE2NODES_MFEM[local_side][2]]`.
4. Fetch the three face-vertex coordinates `p0, p1, p2 = mesh.GetVertex(fn[0..2])` — each a `const real_t*` of length 3.
5. **Raw normal** (unnormalized): `ab[j] = p1[j] − p0[j]`, `ac[j] = p2[j] − p0[j]` for `j = 0..2`; `normal_raw = cross(ab, ac)`. Note: `normal_raw ⊥ ab` by construction of the cross product.
6. **Outward orientation**: compute element centroid `c_elem[j] = (p_elem[0][j] + p_elem[1][j] + p_elem[2][j] + p_elem[3][j]) / 4` and face centroid `c_face[j] = (p0[j] + p1[j] + p2[j]) / 3`. If `dot(normal_raw, c_face − c_elem) < 0`, swap `fn[0]` and `fn[1]` (which swaps `p0` and `p1`) and recompute `ab = p1 − p0` and `normal_raw = cross(ab, ac)` from the swapped coordinates. (Equivalent to flipping both `normal_raw` and `tangent1_raw` sign while preserving the `normal ⊥ tangent1` identity.)
7. **Raw tangent1**: `tangent1_raw = ab` (same vector used in the cross product after any swap). Identity: `dot(normal_raw, tangent1_raw) = 0` exactly up to floating-point rounding of `cross` and `sub`.
8. **Raw tangent2**: `tangent2_raw = cross(normal_raw, tangent1_raw)`.
9. **Normalize**: `surface_area = 0.5 * |normal_raw|`; `normal = normal_raw / |normal_raw|`; `tangent1 = tangent1_raw / |tangent1_raw|`; `tangent2 = tangent2_raw / |tangent2_raw|`. **Do not re-orthogonalize.**
10. **Cell volume** (R2-005 FIX): build a 3×3 matrix `dX` with rows
    ```
    dX[0] = p_elem[1] − p_elem[0]     // a 3-vector
    dX[1] = p_elem[2] − p_elem[0]     // a 3-vector
    dX[2] = p_elem[3] − p_elem[0]     // a 3-vector
    ```
    then `cell_volume = |det(dX)| / 6` where `det(dX)` is the standard 3×3 determinant of the matrix of vertex-coordinate differences. **Do NOT subtract `ev[i]`'s directly** — those are integer vertex indices, not coordinates; use `p_elem[i]` (the fetched coordinates from step 2). A regression test asserts `cell_volume == 1/6` to 1e-14 on the unit tet `{(0,0,0), (1,0,0), (0,1,0), (0,0,1)}`.

Critical per R-004: step 6 uses the `ab` vector directly. Do not subtract a `dot(tangent1_raw, normal_raw)·normal_raw` component.

**Acceptance test P_TOPOLOGY_FRAME_MATCHES_SEISSOL** (§5.1.7) asserts: on a reference tet, `(normal, tangent1, tangent2)` produced by `ComputeCellFaceFrame` match a hand-port of SeisSol's `MeshTools::normalAndTangents` (with MFEM's FACE2NODES, not SeisSol's — the convention differs; see R-005) bit-exactly.

#### 5.1.5 Interior matrix composition (address R-002, R-008)

`BuildInteriorMatrices(normal, tangent1, tangent2, flux, nApNm1, nAmNm1)`:

Given the frame `(n, t1, t2)`:
1. Reuse `GodunovFlux::BuildRotation(n, t1, t2, T)` and `GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv)` from the existing library. **Do not call `GodunovFlux::BuildFrame`** — we supply the frame from topology.
2. Get `A_plus = flux.GetAxPlus()` and `A_minus = flux.GetAxMinus()` (the two new additive const accessors added to `godunov_flux.hpp` per §2.1).
3. Compose:
   - `nApNm1 = T · A_plus · Tinv` (9×9 in global frame)
   - `nAmNm1 = T · A_minus · Tinv` (9×9 in global frame)

Identity (verified by P_MATRIX_IDENTITY, §5.1.7):
```
nApNm1 · I_self + nAmNm1 · I_nbr
   == T · (A_plus · Tinv · I_self + A_minus · Tinv · I_nbr)
   == T · ApplySplitFlux(Tinv·I_self, Tinv·I_nbr, ·)
   == flux_.Interior(n, I_self, I_nbr, F_h)    [when BuildFrame(n) yields the same (t1,t2)]
```

For the orbit-mirror interior face (which sits on the opposite side, with MFEM's `CalcOrtho` giving `-n`), the **other cell's** entry (Elem2's perspective) is built from Elem2's own outward normal (which is `-n` relative to Elem1's perspective). Elem2's `nApNm1_Elem2 / nAmNm1_Elem2` is built identically, using `(n_Elem2, t1_Elem2, t2_Elem2)` from Elem2's topology.

**R2-004 FIX — identity sign**: the correct algebraic identity is
```
Interior(n, L, R) = -Interior(-n, R, L)
```
(opposite sign), derivable directly: with `A_{-n} = -A_n`, we have `A_{-n}^+ = -A_n^-` and `A_{-n}^- = -A_n^+`, hence
```
Interior(-n, R, L)  = A_{-n}^+ · R + A_{-n}^- · L
                    = -A_n^-   · R + -A_n^+   · L
                    = -(A_n^+ · L + A_n^- · R)
                    = -Interior(n, L, R).
```

**Consequence for the two-entry accumulator**: cell-1's F_h (using its outward normal `n`, with `I_self=L`, `I_nbr=R`) and cell-2's F_h (using its outward normal `-n`, with its own `I_self=R`, `I_nbr=L`) satisfy
```
F_h_cell1 + F_h_cell2 = 0.
```
When both sides accumulate with `rhs -= w · shape · F_h`, the net rhs contributions on the two cells are equal-and-opposite, matching MFEM's runtime pattern `rhs1 -= w·shape1·F_h; rhs2 += w·shape2·F_h`. (The plus sign on the Elem2 side of MFEM's runtime corresponds exactly to the minus sign of `F_h_cell2`: `−w·shape2·F_h_cell2 = −w·shape2·(−F_h_cell1) = +w·shape2·F_h_cell1`.)

This conservation identity is enforced by **P_CONSERVATION** (§5.1.7). If a future refactor disturbs the two-entry per-face storage or changes accumulator signs, P_CONSERVATION catches it.

#### 5.1.6 fluxScale / weight reconciliation (unreviewed area #1)

Per the unreviewed-area note on `fluxScale = -2·S/(6·V)` vs MFEM's `w = ip.weight * nor_len`:

- **Phase 1 decision**: we do **not** bake `fluxScale` into the stored matrix. The per-QP weight `w = ip.weight * nor_len` is applied at runtime in the existing accumulation loop, identical to MFEM's current behaviour. This means the stored matrix has units of flux Jacobian (Pa / (Pa·m·s) for stress components, etc.); the integration weight is entirely at the runtime call site.
- **Consequence**: the stored matrix layout differs from SeisSol's `AplusT` (which bakes `-2·S/(6·V)` in). This is expected per §2.5.
- **Verification checklist** added to Step A acceptance: the implementer must document, in a comment at the top of `BuildInteriorMatrices`, the dimensional identity
  ```
  w · nApNm1 · I_self  ==  ip.weight · nor_len · T · A_plus · Tinv · I_self
                       ==  MFEM's current `w * flux_.Interior(nor_unit, I_self, I_nbr, F_h)_self_component`
  ```
  and justify why `nor_len` is NOT baked into `nApNm1`. If the implementer finds a case where baking `nor_len` is necessary (e.g. for numerical stability on highly anisotropic elements), that becomes a Phase 5 deferred task, not a Phase 1 blocker.

- **R3-009 PRECONDITION — p=1 isoparametric mesh (faces affine)**: the precomputed `nApNm1 = T·A⁺·Tinv` is **constant per (cell, side)** and carries no per-QP variation. MFEM's runtime `nor_len = |CalcOrtho(ftr->Face->Jacobian())|` is constant per-face only when the face Jacobian is constant across QPs, which holds on p=1 (straight-sided) isoparametric meshes. On p≥2 isoparametric (curved) meshes, `nor_len` varies per QP and the precomputed matrices would not be equivalent to the runtime flux.
  TPV102 uses straight tets (p=1 isoparametric), so the precondition holds for the target use case. `Init` shall enforce it at runtime with a direct numerical check rather than a speculative MFEM API call:
  ```cpp
  // For each interior/boundary face that passes Init:
  //   Sample ftr->Face->Jacobian()'s |CalcOrtho| at IP(0) and IP(last).
  //   MFEM_VERIFY(std::abs(nor_len_first - nor_len_last) <
  //               1e-12 * nor_len_first,
  //               "PrecomputedFaceFluxes requires p=1 isoparametric "
  //               "(affine) mesh; face " << face_idx << " has "
  //               "QP-varying Jacobian (face is curved).");
  ```
  A dedicated unit test `test_R3_009_p1_precondition` constructs a curved (p=2) mesh fixture and asserts that `Init` aborts.

#### 5.1.7 Acceptance criteria — Step A

Standalone unit test `tests/unit/test_precomputed_fluxes_interior_identity.cpp` (new).

**Pass thresholds**: 1e-12 for bit-exact matrix identities; 1e-14 for topology frame identities (no matrix condition-number amplification).

- [ ] **P_FACE2NODES_SELFCHECK**: `AssertMFEMFace2NodesTable()` passes on the live MFEM binary. If MFEM changes its reference-tet face-vertex ordering in a future release, this test is the early-warning signal.
- [ ] **P_TOPOLOGY_FRAME_MATCHES_SEISSOL**: On a reference tet with vertices `{(0,0,0), (1,0,0), (0,1,0), (0,0,1)}`, for every `local_side ∈ {0,1,2,3}`:
  - Compute `(normal, tangent1, tangent2)` via `ComputeCellFaceFrame`.
  - Compute the same via a hand-port of SeisSol's `MeshTools::normalAndTangents` using the hard-coded `FACE2NODES_MFEM` (not SeisSol's `FACE2NODES`).
  - Assert each component agrees to 1e-14 (cross-product and normalization rounding floor).
- [ ] **P_FRAME_ORBIT**: For every element and every local side on the M_ref 4×2×4 Cartesian fixture:
  - `normal` is axis-aligned (±x, ±y, or ±z component dominates) consistent with the fixture's axis-aligned faces;
  - `dot(normal, c_face - c_elem) > 0` always (outward);
  - `dot(normal, tangent1) < 1e-14` and `dot(normal, tangent2) < 1e-14` and `dot(tangent1, tangent2) < 1e-14` (frame orthonormality).
- [ ] **P_MATRIX_IDENTITY**: For 100 random 9-vectors `(I_self, I_nbr)` with magnitudes in `[1e-3, 1e3]`, for every interior non-fault cell-side on M0 and M_ref:
  - `F_h_precomputed[c] = (nApNm1 · I_self + nAmNm1 · I_nbr)[c]` agrees with `flux_.Interior(n_outward_for_this_cell, I_self, I_nbr, F_h_runtime)[c]` to 1e-12.
- [ ] **P_CONSERVATION** (new, R-003): For 100 random `(I_self, I_nbr)`, for every interior non-fault face on M0:
  - Let `E1_contribution = nApNm1_from_Elem1 · I_self + nAmNm1_from_Elem1 · I_nbr` (flux as seen by Elem1, using its own outward normal).
  - Let `E2_contribution = nApNm1_from_Elem2 · I_nbr + nAmNm1_from_Elem2 · I_self` (flux as seen by Elem2, with swapped roles because I_nbr is now Elem2's own state).
  - Assert `||E1_contribution + E2_contribution||_∞ < 1e-12` — flux conservation. This is the identity SeisSol gets automatically from its single-entry design and that our two-entry design must enforce.
- [ ] **P_SEISSOL_VALUE** (new, R-002; formula REWRITTEN per R2-001): On a single-tet fixture with hand-picked material `(Zp=..., Zs=..., ρ=..., λ=..., μ=...)` and hand-picked `(n, I_self, I_nbr)`:
  - **R2-001 context**: the previous formulation in this acceptance `F_h_seissol = T · QgodLocal_T^T · starMatrix · Tinv · I_self + T · QgodNeighbor_T^T · starMatrix · Tinv · I_nbr` (with `QgodLocal_T = I − (matR·χ·matR⁻¹)^T`, `QgodNeighbor_T = (matR·χ·matR⁻¹)^T`, `starMatrix = A_x`) algebraically simplifies to `T · A_x⁻ · Tinv · I_self + T · A_x⁺ · Tinv · I_nbr`, which has A⁺ and A⁻ **swapped** versus MFEM's convention `F_h = A⁺·I_self + A⁻·I_nbr`. This is because SeisSol's `fluxScale = -2·S/(6·V)` (CellLocalMatrices.cpp:256) and the stored `starMatrix = A_x^T` (CellLocalMatrices.cpp:211, `getTransposedCoefficientMatrix`) together flip the effective upwind assignment in SeisSol's runtime rhs sign convention. The QgodLocal layered form is NOT a clean flux-value test. Use the direct eigendecomposition below instead.
  - **Direct eigendecomposition reference formula** (the only formula Phase 1 acceptance uses):
    ```
    // Step 1 — reference-frame flux Jacobian A_x (derived from MFEM's
    // isotropic elastic model with parameters ρ, λ, μ).  Uses MFEM's
    // component ordering (SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ).
    //
    // Step 2 — eigendecompose A_x:
    eigendecompose(A_x, matR, Lambda)   // 9x9 real eigenvectors, 9 eigenvalues
    // Lambda contains three pairs (+Cp, -Cp), (+Cs, -Cs), (+Cs, -Cs) plus
    // three zero eigenvalues for the static characteristics.
    //
    // Step 3 — positive/negative eigenvalue selector:
    chi_plus  = diag(Lambda[i] > 0 ? 1 : 0)   // 9-entry {0,1} selector
    chi_minus = diag(Lambda[i] < 0 ? 1 : 0)   // 9-entry {0,1} selector
    // (Zero eigenvalues belong to neither; they don't contribute to A±.)
    //
    // Step 4 — split flux Jacobians in the reference frame:
    A_plus_x  = matR · diag(Lambda · chi_plus)  · matR_inv
    A_minus_x = matR · diag(Lambda · chi_minus) · matR_inv
    // Sanity: A_plus_x + A_minus_x == A_x to ~1e-14 (assert this).
    //
    // Step 5 — rotate into the face-local outward frame (n, t1, t2)
    // using the SAME T/Tinv construction as PrecomputedFaceFluxes:
    T, Tinv = BuildRotation(n, t1, t2)
    //
    // Step 6 — reference flux in global frame (MFEM's convention:
    // I_self contributes A^+, I_nbr contributes A^-):
    F_h_reference = T · A_plus_x  · Tinv · I_self
                  + T · A_minus_x · Tinv · I_nbr
    ```
  - Assert `||F_h_reference − (nApNm1 · I_self + nAmNm1 · I_nbr)||_∞ < 1e-10` for 20 random `(n, I_self, I_nbr)` samples with `n` drawn from orbit-equivalent outward normals on M_ref.
  - **Trivial-case pre-check**: additionally run with `I_self = e_j` (standard basis), `I_nbr = 0` for each `j ∈ {0, ..., 8}`, and assert `F_h_reference[c] == nApNm1(c, j)` to 1e-12. This independently verifies each column of `nApNm1`.
  - If the general-case assertion fails at > 1e-10 but the trivial-case pre-check passes: MFEM's eigenvector matrix `matR_MFEM` has a column-ordering or sign convention different from the convention used in `A_plus_x` / `A_minus_x` above (which should be reconciled to MFEM's own split via `GodunovFlux::GetAxPlus()` / `GetAxMinus()`). Document the basis-change map before proceeding — but this scenario should not occur because the reference formula uses MFEM's own `A_x` (not SeisSol's). If it does occur, inspect MFEM's `A_x_plus_` / `A_x_minus_` construction in `godunov_flux.cpp` for the eigendecomposition convention.
  - **This test is the "SeisSol flux-value equivalence receipt"** — but framed as a direct eigendecomposition rather than a QgodLocal layering, because the QgodLocal form carries SeisSol-specific sign/transpose conventions that are easy to get wrong in the hand-port.
- [ ] **P_BIMATERIAL_GUARD**: Constructing a bimaterial fixture and calling `Init` aborts cleanly with `MFEM_ABORT` naming the face index and the mismatched impedances.

#### 5.1.8 Stop / go — Step A

- **GO**: every P_* above passes within its threshold → proceed to Step B.
- **STOP (basis-convention mismatch)**: P_SEISSOL_VALUE fails at > 1e-10 and P_MATRIX_IDENTITY passes. The precomputed path matches MFEM's runtime path but not SeisSol's formula. Document the eigenvector column/sign difference and decide whether the plan's equivalence claim (§2.5) is still acceptable. **Do not proceed to Phase 2a** until the claim is revised.
- **STOP (conservation broken)**: P_CONSERVATION fails. The two-entry design does not preserve flux conservation in MFEM's A± convention. Either fix MFEM's A± construction or restructure to a single-entry design. Either is a major rework.
- **STOP (FACE2NODES drift)**: P_FACE2NODES_SELFCHECK fails. MFEM's tet face table has changed relative to the hard-coded values. Update `FACE2NODES_MFEM` and rerun.
- **STOP (topology frame ULP drift)**: P_TOPOLOGY_FRAME_MATCHES_SEISSOL fails at > 1e-14. A re-orthogonalization or subtle scaling difference snuck in. Re-read §5.1.4 step by step; fix before proceeding.
- **STOP (R2-005 cell_volume integer arithmetic)**: unit-tet regression check `cell_volume == 1/6` fails by orders of magnitude. The implementer subtracted vertex indices instead of coordinates. Fix per §5.1.4 step 10.

### 5.1.9 Step A — new regression test for fault-branch preservation (R2-008)

**Test**: `tests/unit/test_R008_fault_branch_preserved_after_dispatch_reorder.cpp` (new).

**Goal**: assert the fault-branch body is lifted verbatim into the restructured if/else (§6.2), with no regressions in fault DOF state output.

**Procedure**:
1. Before landing the Phase 2a patch, capture a baseline snapshot: on M0 with the current (pre-restructure) `WaveOperator`, run one ADER step from `Q = 0` with the uniform persistent-nucleation fixture (same as `test_adjacent_triangle_fault_first_step_audit`), and record every field of every `(*fault_dof_data_)[i]` (slip_rate, V1, V2, tau1_corr, tau2_corr, sigma_n_corr, psi) plus the per-QP `I_plus_local / I_minus_local / Q_imp_plus / Q_imp_minus / F_h_plus / F_h_minus` samples from `RunFaultOnlyAudit`. Serialize to a `.txt` or `.bin` snapshot committed alongside the test.
2. After landing the Phase 2a patch, with `use_precomputed_face_fluxes_ = false` (default), re-run the same step.
3. Compare every field against the snapshot byte-for-byte (or to 1e-15 relative, accommodating floating-point associativity drift from any reordered arithmetic).
4. Assert zero-diff. Any nonzero diff = regression → fault-branch body was inadvertently simplified or reordered.

**Pass**: zero diff on every field to 1e-15.
**Fail**: any field differs beyond 1e-15 → the implementer lost a line of the fault branch during the restructure. Revert and retry.

This test must pass on M0 and M_ref serial before Phase 2a is accepted.

### 5.2 Step B — Boundary face operator

Address R-006 (FreeSurface variants), R-007 (Absorbing semantics), R-008 (Γ placement).

#### 5.2.1 FreeSurfaceGodunov (mandatory for SeisSol parity)

Matches SeisSol's `FaceType::FreeSurface` semantics via `ElasticSetup.h:141-144 + Common.h:286-294`:

- In face-local rotated frame, traction indices `{SXX, SXY, SXZ} = {0, 3, 5}` and velocity indices `{VX, VY, VZ} = {6, 7, 8}`.
- `qGodLocal[traction, traction] = I_3`
- `qGodLocal[velocity, traction] = S = -matR21 · matR11⁻¹` (where `matR11, matR21` are blocks of the eigenvector matrix)
- `qGodNeighbor` is unused (traction on ghost = 0, so there is no neighbor contribution)
- Composition: `nApNm1 = T · (A_plus · qGodLocal_rotated_frame_action_on_I_self)` — fold qGodLocal into the A_plus matvec in the rotated frame.

Equivalently, the plan can reuse MFEM's `flux_.FreeSurfaceGodunovTotal(n_outward, I_self, bulk_bg=0, F_h)` at Init time to probe the linear operator and then fit the 9×9 matrix by probing with 9 standard basis vectors. This is the cleanest approach given MFEM's existing implementation is already correct — we just capture its linear action as a precomputed 9×9:

```cpp
// R4-003 FIX — signature removed the unused (t1, t2) parameters; they
// were dead code in the probe construction and their presence invited
// a maintainer trap (using them inconsistently with the probe result).
// The probe's topology-based `n` argument is what resolves the
// orientation ambiguity; (t1, t2) are immaterial for isotropic elastic
// free-surface BCs — see the rotation-invariance argument in §0‴.
void PrecomputedFaceFluxes::BuildBoundaryMatricesGodunov(
    const real_t n[3],
    const GodunovFlux &flux,
    DenseMatrix &nApNm1)
{
   // R2-009 PRECONDITION: flux.FreeSurfaceGodunovTotal must be LINEAR
   // in I_self when bulk_bg = 0.  Verified from godunov_flux.cpp:488–543:
   //   Q_pert_rot = Q_self_rot − Q_bg_rot
   //   Apply linear Godunov projection to Q_pert_rot
   //   Q_ghost_rot = Q_bg_rot + projected_Q_pert_rot
   // At Q_bg = 0: Q_pert_rot = Q_self_rot, and the Q_bg-dependent
   // nonlinearity vanishes.  Same holds for FreeSurfaceTotal (the γ-
   // mirror variant) and AbsorbingTotal.  A defensive MFEM_VERIFY
   // enforces this precondition at Init time; a companion runtime
   // check fires in AddBoundaryFaceRhs if `bulk_bg_scaled` is not zero.
   //
   // Probe: apply flux_.FreeSurfaceGodunovTotal to each of 9 standard
   // basis vectors; each result is a column of nApNm1.
   //
   // R4-003 NOTE: FreeSurfaceGodunovTotal internally calls
   // GodunovFlux::BuildFrame(n, t1_bf, t2_bf) — a Gram-Schmidt against
   // up = (0,0,1) or (1,0,0).  This produces a DIFFERENT (t1, t2) pair
   // than ComputeCellFaceFrame's topology-supplied (t1_topo, t2_topo).
   // For isotropic elastic free-surface BCs (both Godunov projection
   // and γ-mirror) the flux is INVARIANT under in-plane rotations of
   // (t1, t2) around `n`, so F_h(n, t1_bf, t2_bf) == F_h(n, t1_topo,
   // t2_topo) bit-exactly.  The orientation-ambiguity source that
   // caused the SXZ drift (H_BFACE_NOR) was the (+n vs −n) sign of
   // MFEM's runtime `nor`, NOT the (t1, t2) pair; supplying the
   // topology-outward `n` here eliminates it.
   real_t bulk_bg_zero[NUM_STATE] = {0};
   for (int c = 0; c < NUM_STATE; c++) {
      MFEM_VERIFY(bulk_bg_zero[c] == 0.0,
                  "BuildBoundaryMatricesGodunov: probe precondition "
                  "requires bulk_bg = 0 for linearity; Phase 1 does not "
                  "support nonzero bulk_bg.  See §5.1.6.");
   }
   nApNm1.SetSize(NUM_STATE, NUM_STATE);
   for (int j = 0; j < NUM_STATE; j++)
   {
      real_t I_self[NUM_STATE] = {0};
      I_self[j] = 1.0;
      real_t F_h[NUM_STATE];
      flux.FreeSurfaceGodunovTotal(n, I_self, bulk_bg_zero, F_h);
      for (int i = 0; i < NUM_STATE; i++) { nApNm1(i, j) = F_h[i]; }
   }
}

// Analogous signature changes apply to BuildBoundaryMatricesGamma and
// BuildBoundaryMatricesAbsorbing — remove (t1, t2).  The rotation-
// invariance argument in §0‴ covers all three isotropic-elastic BCs.
```

This construction **preserves MFEM's exact runtime behaviour** under `Q_bg = 0` (which is TPV102's invariant; see `tpv102_driver.cpp:523`). Per §5.1.6 no fluxScale bake-in.

**Probe-linearity regression test** (new, R2-009): `tests/unit/test_R009_probe_linearity.cpp` calls
```
flux.FreeSurfaceGodunovTotal(n, I_self, 0, F_h_single);
flux.FreeSurfaceGodunovTotal(n, 2·I_self, 0, F_h_double);
```
for a handful of random `(n, I_self)` and asserts `||F_h_double − 2·F_h_single||_∞ < 1e-14`. Same check for `FreeSurfaceTotal` and `AbsorbingTotal`. If any of these fail under `bulk_bg = 0`, the probe-based construction is invalid and the plan must fall back to an analytic composition (§5.2.3's Γ formula or equivalent).

#### 5.2.2 Absorbing (MFEM-native first-order Sommerfeld) — R-007 resolution

Per §0, the Absorbing semantics is the **MFEM-native `F_h = A⁺ · I_self`** (Sommerfeld radiation condition), **not** SeisSol's `Outflow` which radiates full `A`. Rationale:

- TPV102's absorbing faces (bottom of cube) are in far field, and the reflection-coefficient difference between MFEM's Sommerfeld and SeisSol's Outflow on a 40-km TPV102 domain is well below the required accuracy at the hypocenter.
- Changing absorbing semantics would be a physics change, not a bug fix. Out of scope.

Composition uses the same probe-based construction as §5.2.1:
```cpp
flux.AbsorbingTotal(n, I_self_basis_j, bulk_bg_zero, F_h);
```

#### 5.2.3 FreeSurfaceGamma (legacy; optional) — R-006 resolution

MFEM's γ-mirror free surface is **not** in SeisSol. It is retained for backward compatibility with pre-Phase-5 runtime fallback. If the driver requests Gamma mode (via `SetFreeSurfaceBCMode(Gamma)`), the precomputed matrix is built using the probe approach against `flux_.FreeSurfaceTotal(n, basis_j, 0, F_h)`.

Γ composition formula (R-008 clarification):
```
A_face_local = A_plus + A_minus · Γ        (γ acts on the ghost state in rotated frame)
nApNm1       = T · A_face_local · Tinv
Γ            = diag(-1, +1, +1, -1, +1, -1, +1, +1, +1)   in state order
               (SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ)
```
The existing `godunov_flux.cpp:466`'s γ array is the ground truth reference. The probe approach in §5.2.1 sidesteps this by using the already-correct runtime implementation; this is preferred over re-deriving the γ composition by hand.

#### 5.2.4 Acceptance criteria — Step B

New test `tests/unit/test_precomputed_fluxes_boundary_identity.cpp`.

- [ ] **P_BFACE_GODUNOV** (primary, matches SeisSol FaceType::FreeSurface): For 100 random `I_self`, the precomputed Godunov-mode dispatch matches `flux_.FreeSurfaceGodunovTotal(n_outward, I_self, bulk_bg=0, F_h)` to 1e-12, on every boundary face of M_ref.
- [ ] **P_BFACE_SEISSOL_FREE_SURFACE** (new, R-006): On a single-tet fixture with a boundary face, hand-compute the Godunov-projection matrix using SeisSol's `S = -matR21 · matR11⁻¹` formula and assert `||S_seissol · I_self_velocity_block - precomputed_velocity_block||_∞ < 1e-10`. This is the SeisSol-parity receipt for the free-surface BC.
- [ ] **P_BFACE_GAMMA** (legacy): same criterion but against `flux_.FreeSurfaceTotal`; runs only when the driver opts into Gamma mode via `SetFreeSurfaceBCMode(Gamma)`.
- [ ] **P_BFACE_ABSORBING** (MFEM-native Sommerfeld): same criterion against `flux_.AbsorbingTotal`. No SeisSol parity check — the semantics explicitly differ per R-007.
- [ ] **P_BFACE_NONZERO_BG_ABORT**: calling `AddBoundaryFaceRhs` with `||bulk_bg_scaled||_∞ > 1e-15` aborts (Phase 1 restriction; consistent with TPV102's Q_bg=0 invariant).

#### 5.2.5 Stop / go — Step B

- **GO**: P_BFACE_GODUNOV + P_BFACE_SEISSOL_FREE_SURFACE + P_BFACE_ABSORBING all pass → proceed to Phase 2a.
- **STOP (Godunov mismatch)**: P_BFACE_GODUNOV fails. The probe-based construction is incorrect (should be bit-equivalent to the runtime flux). Debug the basis-probe pattern.
- **STOP (SeisSol parity gap)**: P_BFACE_GODUNOV passes but P_BFACE_SEISSOL_FREE_SURFACE fails. MFEM's Godunov implementation uses a different eigenvector basis than SeisSol's. Document the gap; decide whether the plan still claims "SeisSol parity" on this BC.

### 5.3 Step A / Step B shared `Init` sequence (post-R3 fixes)

`Init(mesh, fes, bc, flux, fs_mode, fault_face_set, face_bdr_attr, shared_mesh_face_set)`:

1. Call `AssertMFEMFace2NodesTable()` (R3-001 FIX: asserts against `mfem::Geometry::Constants<mfem::Geometry::TETRAHEDRON>::FaceVert`).
2. Assert p=1 isoparametric precondition (R3-009 FIX; see §5.1.6).
3. For every element `e` and every local side `s ∈ {0,1,2,3}`:
   a. Compute topology frame via `ComputeCellFaceFrame`.
   b. Resolve mesh face index via MFEM's element→face map. Use `mesh.GetElementFaces(e, el_faces, el_orients)` (which returns the 4 face indices for this tet) and take `face_idx = el_faces[s]`.
   c. **Skip shared faces** (R3-007 FIX): if `shared_mesh_face_set.count(face_idx) > 0`, `continue` — Phase 1 does not touch shared faces; Phase 2b extends Init to handle them separately.
   d. **Skip fault faces**: if `fault_face_set.count(face_idx) > 0`, `continue` — fault branch untouched.
   e. Get the face transformation: `FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(face_idx);` (R3-008 FIX: `mesh` is non-const for this call).
   f. **Boundary-vs-interior classification** (R3-002 + R3-003 FIX — mirrors `wave_operator.inl:872, 874`):
      ```cpp
      int e2 = ftr->Elem2No;
      int bdr_attr = face_bdr_attr[face_idx];        // R3-003: uses precomputed table
      bool is_boundary = (e2 < 0) && (bdr_attr > 0);
      // (shared_mesh_face_set already skipped in step c above, so
      //  `e2 < 0 && bdr_attr > 0 && !is_shared` reduces to the predicate
      //  above, since `is_shared` can no longer be true here.)
      ```
   g. If `is_boundary` → **boundary branch**:
      **R2-006 FIX — explicit classification rule**. Use `bdr_attr` from step f directly:
      ```cpp
      FacePrecomputedBC bc_class;
      if (bc.natural_attrs.count(bdr_attr) > 0) {
         // Natural / zero-traction BC.  MFEM's runtime branches on fs_mode:
         //   fs_mode == FreeSurfaceBCMode::Godunov → FreeSurfaceGodunovTotal
         //   fs_mode == FreeSurfaceBCMode::Gamma   → FreeSurfaceTotal (γ-mirror)
         bc_class = (fs_mode == FreeSurfaceBCMode::Godunov)
                    ? FacePrecomputedBC::FreeSurfaceGodunov
                    : FacePrecomputedBC::FreeSurfaceGamma;
      }
      else if (bc.absorbing_attrs.count(bdr_attr) > 0) {
         bc_class = FacePrecomputedBC::Absorbing;
      }
      else if (bdr_attr == bc.fault_attr) {
         // Unreachable — already skipped in step d via fault_face_set.
         MFEM_ABORT("PrecomputedFaceFluxes::Init: face " << face_idx
                    << " is 1-sided with fault_attr; should have been "
                    << "skipped in step d.");
      }
      else {
         MFEM_ABORT("PrecomputedFaceFluxes::Init: face " << face_idx
                    << " has unclassified boundary attribute " << bdr_attr
                    << " (not in natural_attrs, absorbing_attrs, or == fault_attr).");
      }
      ```
      Dispatch to the matching `BuildBoundaryMatricesGodunov(n, flux, nApNm1)` / `BuildBoundaryMatricesGamma(n, flux, nApNm1)` / `BuildBoundaryMatricesAbsorbing(n, flux, nApNm1)` variant based on `bc_class`. (R4-003: the call takes only `(n, flux, nApNm1)` — no `(t1, t2)`.)
      Push a `FaceEntry` with `element = e`, `local_side = s`, `face_idx = face_idx`, `neighbor_elem = -1`, `bc = bc_class`, `normal = n_topo`, `tangent1 = tangent1_topo`, `tangent2 = tangent2_topo` (diagnostic-only for boundary — see §5.1.2 `FaceEntry` comment), and `nApNm1` from the probe, `nAmNm1 = 0` (zero-initialized 9×9).
   h. Else → **interior non-fault branch**:
      - `homog_ok` check on both sides' impedances. In Phase 1 material is homogeneous, so both cells have the same `(Zp, Zs)` from `flux_`'s internal material parameters (see `GodunovFlux::Zp_`, `Zs_`). `MFEM_ABORT` on any detected mismatch (future bimaterial extension).
      - Call `BuildInteriorMatrices(normal, tangent1, tangent2, flux, nApNm1, nAmNm1)`.
      - **R3-004 FIX — neighbor element**: compute
        ```cpp
        int neighbor_elem = (ftr->Elem1No == e) ? ftr->Elem2No : ftr->Elem1No;
        ```
        (`ftr->ElemNNo` in the prior plan text was a typo; the MFEM fields are `Elem1No` / `Elem2No`, per `mfem/mesh/mesh.hpp:178`.)
      - Push a `FaceEntry` with `element = e`, `local_side = s`, `face_idx = face_idx`, `neighbor_elem = neighbor_elem`, `bc = Interior`, and the computed frame / matrices.
4. **Populate `face_elem_to_entry_` keyed by `(face_idx, element)`** (R3-005 FIX — previous plan said `face_side_to_entry_ keyed by (face_idx, local_side)`, which contradicts §5.1.2's R2-002 rename):
   ```cpp
   for (int i = 0; i < entries_.size(); i++) {
      const FaceEntry &fe = entries_[i];
      MFEM_VERIFY(fe.element >= 0, "Phase 1: all caller_elems are non-negative");
      long long key = EncodeKey(fe.face_idx, fe.element);   // see §5.1.2 R3-006
      auto [it, inserted] = face_elem_to_entry_.emplace(key, i);
      MFEM_VERIFY(inserted, "Duplicate (face_idx, element) entry at "
                  << fe.face_idx << ", " << fe.element);
   }
   ```
   Each interior non-fault face produces two entries (one per adjacent element) with the same `face_idx` but different `element` keys; each boundary face produces one entry.
5. Set `initialized_ = true`.

## 6. Phase 2a — ADER opt-in dispatch switch for local non-fault faces (Step C-local)

### 6.1 Goal
After Phase 2a: `ComputeADERFaceFluxRHS` routes local (non-shared) non-fault interior and boundary face flux through `PrecomputedFaceFluxes` when `use_precomputed_face_fluxes_` is true, preserving the existing runtime path when false. Fault dispatch is unchanged regardless of the flag.

### 6.2 Dispatch ordering (FIX for R-001)

Inside `ComputeADERFaceFluxRHS`'s per-QP loop, the branch structure is:

```cpp
// ============================================================
// PRE-EXISTING LOCALS at outer scope (do NOT redeclare; do NOT
// rename; the fix-agent must reuse them as-is):
//   int e1, e2                  — from ftr->Elem1No / Elem2No
//                                 (wave_operator.inl:1777-1778)
//   int ndof                    — from fe1->GetDof()
//                                 (wave_operator.inl:1789)
//   int dof_offset1             — from e1 * ndof_per_el_
//                                 (wave_operator.inl:1790)
//   Vector shape1               — from fe1->CalcShape
//                                 (wave_operator.inl:1811-1812)
//   real_t I_self[NUM_STATE]    — assembled from Q · shape1 at this QP
//                                 (wave_operator.inl:1814-1822)
//   real_t w                    — ip.weight * nor_len
//                                 (wave_operator.inl:1806)
//   real_t bulk_bg_scaled[NUM_STATE] — dt · bulk_bg (computed at the
//                                 top of ComputeADERFaceFluxRHS,
//                                 wave_operator.inl:1758-1762)
//   bool is_boundary            — (e2 < 0) && (bdr_attr > 0), from
//                                 wave_operator.inl:1784
//   int  f                      — the mesh face index, the loop var
//
// HOISTED per R5-001 (see §6.2.hoist below) — NOT pre-existing at
// outer scope in the current file; the plan's Phase 2a patch
// introduces them at outer scope:
//   int    dof_offset2          — e2 * ndof_per_el_, meaningful only
//                                 when !is_boundary
//   Vector shape2               — fe2->CalcShape, meaningful only
//                                 when !is_boundary
//   real_t I_nbr[NUM_STATE]     — assembled from Q · shape2 at this QP,
//                                 meaningful only when !is_boundary
//   bool   is_fault             — (bdr_attr == bc_.fault_attr) &&
//                                 (bc_.fault_attr > 0) — HOISTED so
//                                 the outermost if/else can branch on
//                                 it; false on boundary faces by
//                                 construction of the hoist guard
// ============================================================

// === §6.2.hoist (R5-001 FIX) ================================
// The existing ComputeADERFaceFluxRHS declares dof_offset2, shape2,
// I_nbr, and is_fault INSIDE its `else (interior)` branch
// (wave_operator.inl:1870-1887).  The new dispatch ordering below
// must branch on is_fault as the OUTERMOST predicate, so we hoist
// these four declarations to outer scope.  We guard the computation
// with !is_boundary so boundary faces never evaluate
// e2 * ndof_per_el_ with e2 < 0 (which would produce a negative
// offset and out-of-bounds reads).
//
// PLACEMENT: insert immediately after the existing `bool is_boundary
// = (e2 < 0) && (bdr_attr > 0);` line at wave_operator.inl:1784,
// BEFORE the pre-existing boundary branch that currently starts at
// line 1826.
int    dof_offset2 = -1;                 // invalid sentinel
Vector shape2;                           // empty until populated
real_t I_nbr[NUM_STATE] = {0};           // zero until populated
bool   is_fault = false;                 // false on boundary faces

if (!is_boundary)
{
   const FiniteElement *fe2 = fes_->GetFE(e2);
   dof_offset2 = e2 * ndof_per_el_;
   shape2.SetSize(ndof);
   IntegrationPoint ip2;
   ftr->Loc2.Transform(ip, ip2);
   fe2->CalcShape(ip2, shape2);
   for (int c = 0; c < NUM_STATE; c++)
   {
      I_nbr[c] = 0.0;
      for (int i = 0; i < ndof; i++)
      {
         I_nbr[c] += shape2(i) * I_data[c * ndof_total_
                                         + dof_offset2 + i];
      }
   }
   is_fault = (bdr_attr == bc_.fault_attr) && (bc_.fault_attr > 0);
}
// is_fault remains `false` on boundary faces by construction of the
// `!is_boundary` guard.  A 1-sided face that IS tagged with
// bc_.fault_attr is a mesh-configuration error and would be caught
// by §5.3's Init classification (§5.3 step 2g MFEM_ABORTs on a
// 1-sided fault face) before the dispatch ever sees it; the runtime
// dispatch therefore cannot enter the fault branch for a boundary
// face.
//
// The equivalent hoist applies to ComputeADERSharedFaceFluxRHS
// (Phase 2b — §7), and to ComputeFaceFluxRHS (Phase 4, if entered).
// Both sites must mirror this pattern line-for-line.
// === end §6.2.hoist =========================================

if (is_fault)
{
   // R2-008 FIX — EXISTING fault dispatch body.  PRESERVE VERBATIM from
   // the current ComputeADERFaceFluxRHS fault branch
   // (wave_operator.inl:1887–2020) for LOCAL fault faces.  Do NOT
   // abbreviate, do NOT simplify, do NOT replace with a one-line call.
   // Every line is load-bearing:
   //   - LookupInteriorFaultBasisIndex(f)
   //   - FaultBasisQPData reconstruction
   //   - DOFData pointer resolution via fault_face_dof_offset_
   //   - can_n / can_t1 / can_t2 BP5 canonical-frame reconstruction
   //     (per CLAUDE.md R-801)
   //   - per-side accumulation loop
   //   - MFEM_ASSERT tripwires
   // The outer if/else RESTRUCTURE (this block being moved to first
   // position) is the ONLY change to the fault branch; the branch body
   // is lifted unchanged from the current file.
   // test_R008_fault_branch_preserved_after_dispatch_reorder (see §5.4)
   // is the regression gate that proves no machinery was dropped.
   //
   // The shared-fault branch in ComputeADERSharedFaceFluxRHS
   // (wave_operator.inl:2278–2340, MPI rank canonicalization) receives
   // the SAME structural restructure, with identical preserve-VERBATIM
   // discipline.  See §7.
   ... /* lift existing fault-branch body here, VERBATIM */
}
else if (use_precomputed_face_fluxes_ && is_boundary)
{
   precomputed_face_fluxes_.AddBoundaryFaceRhs(
      f, I_self, bulk_bg_scaled, w, shape1.GetData(),
      ndof, dof_offset1, ndof_total_, rhs);
}
else if (use_precomputed_face_fluxes_)   // interior non-fault
{
   // R2-003 FIX — use `e1` / `e2` (pre-existing locals), NOT `elem1` /
   // `elem2` (undefined).  Call AddInteriorFaceRhs twice: once from
   // Elem1's perspective (caller_elem = e1), once from Elem2's
   // perspective (caller_elem = e2).  This mirrors MFEM's current
   // two-sided accumulation.  Each call resolves its FaceEntry via
   // (face_idx, caller_elem_id) — see §5.1.2 face_elem_to_entry_ (R2-002).
   precomputed_face_fluxes_.AddInteriorFaceRhs(
      f, /*caller_elem=*/e1, I_self, I_nbr, w, shape1.GetData(),
      ndof, dof_offset1, ndof_total_, rhs);
   precomputed_face_fluxes_.AddInteriorFaceRhs(
      f, /*caller_elem=*/e2, I_nbr, I_self, w, shape2.GetData(),
      ndof, dof_offset2, ndof_total_, rhs);
}
else   // runtime path, default
{
   // EXISTING runtime path — flux_.Interior(nor, I_self, I_nbr, F_h)
   // for interior non-fault, or flux_.FreeSurfaceTotal /
   // FreeSurfaceGodunovTotal / AbsorbingTotal for boundary.  Lifted
   // VERBATIM from the current non-fault branches (R2-008 discipline
   // applies here too: the existing is_boundary → AbsorbingTotal /
   // FreeSurface* dispatch and the interior else branch's per-DOF
   // accumulation `rhs -= w · shape1(i) · F_h[c]` etc. must be
   // preserved as-is when the flag is off).
   ... /* lift existing non-fault dispatch body here, VERBATIM */
}
```

The key structural fix per **R-001**: `is_fault` is the **outermost** condition. It fires independently of `use_precomputed_face_fluxes_`. The opt-in flag only toggles between precomputed-non-fault and runtime-non-fault paths. There is no silent no-op path.

Per **R2-003**: `e1` and `e2` are the pre-existing loop locals (`int e1 = ftr->Elem1No;` and `int e2 = ftr->Elem2No;` at `wave_operator.inl:1777–1778` of the existing file). The implementer does NOT rename them to `elem1` / `elem2`. A `grep -n "elem1\|elem2" wave_operator.inl` after the patch must produce no new matches; the existing `e1` / `e2` are reused at the dispatch site.

Per **R2-008**: the ellipses `... /* lift existing ... body here, VERBATIM */` are literal instructions to preserve the existing body unchanged. The fix-agent does NOT write `// calls …` in place of the real code; it copies the existing lines into the new if/else ordering.

### 6.3 `UsePrecomputedFaceFluxes(bool)` plumbing

```cpp
// wave_operator.hpp (new public method; default off)
void UsePrecomputedFaceFluxes(bool enable);
bool UsingPrecomputedFaceFluxes() const { return use_precomputed_face_fluxes_; }

private:
   bool use_precomputed_face_fluxes_ = false;
   mutable PrecomputedFaceFluxes precomputed_face_fluxes_;

   // R4-002 FIX — fault_face_set_ is a NEW member whose population
   // is specified in §6.3a (WaveOperator ctor).  Populated eagerly
   // from the existing `fault_interior_faces_` and `fault_shared_faces_`
   // arrays right after those are constructed; independent of the
   // opt-in flag.  If the plan fails to populate it, Init's §5.3
   // step 3d would treat fault faces as interior (silent bug).
   std::set<int> fault_face_set_;
```

`UsePrecomputedFaceFluxes(bool enable)` body (R4-001 FIX — all 8 required `Init` params passed; R5-002 FIX — `fault_face_set_` populated here, not in the ctor):
```cpp
void WaveOperator<MeshType>::UsePrecomputedFaceFluxes(bool enable)
{
   if (enable)
   {
      // === §6.3a (R5-002 FIX) — late population of fault_face_set_ ===
      //
      // fault_interior_faces_ and fault_shared_faces_ are populated by
      // the WaveOperator ctor (wave_operator.inl:260-269 and 273-283
      // respectively).  Populating fault_face_set_ here — not in the
      // ctor — is robust against future refactors that might move the
      // fault-face array population to a later helper, and it keeps
      // the initialization at one unambiguous call site.
      //
      // If a future caller reaches UsePrecomputedFaceFluxes(true) with
      // empty fault-face arrays, the MFEM_VERIFY below catches it
      // loudly — indicating either a driver-ordering bug or a codebase
      // refactor that lost the population step.
      if (bc_.fault_attr > 0)
      {
         MFEM_VERIFY(fault_interior_faces_.Size() > 0 ||
                     fault_shared_faces_.Size() > 0,
                     "WaveOperator::UsePrecomputedFaceFluxes(true): "
                     "fault_interior_faces_ and fault_shared_faces_ are "
                     "both empty despite bc_.fault_attr > 0; ctor's "
                     "fault-face population at wave_operator.inl:260-"
                     "269, 273-283 did not run or was refactored away.");
         fault_face_set_.clear();
         for (int i = 0; i < fault_interior_faces_.Size(); i++)
         {
            fault_face_set_.insert(fault_interior_faces_[i]);
         }
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            fault_face_set_.insert(fault_shared_faces_[i]);
         }
      }
      else
      {
         // No fault configured — empty set is correct.
         fault_face_set_.clear();
      }
      // === end §6.3a population =====================================

      if (!precomputed_face_fluxes_.IsInitialized())
      {
         // Sanity — the three precomputed tables Init consumes.
         MFEM_VERIFY(!fault_face_set_.empty() || bc_.fault_attr <= 0,
                     "WaveOperator::UsePrecomputedFaceFluxes: "
                     "fault_face_set_ is empty but bc_.fault_attr > 0; "
                     "§6.3a population was skipped.");
         MFEM_VERIFY(face_bdr_attr_.size() ==
                     static_cast<size_t>(mesh_.GetNumFaces()),
                     "WaveOperator::UsePrecomputedFaceFluxes: "
                     "face_bdr_attr_ not sized to mesh.GetNumFaces(); "
                     "WaveOperator ctor bookkeeping is incomplete.");

         precomputed_face_fluxes_.Init(
            mesh_,                                         // Mesh &     (R3-008)
            *fes_,                                         // const FES&
            bc_,                                           // const BoundaryConfig&
            flux_,                                         // const GodunovFlux&
            free_surface_bc_mode_,                         // FreeSurfaceBCMode
            fault_face_set_,                               // R4-002 / R5-002
            face_bdr_attr_,                                // R3-003 — pre-existing ctor member
            shared_mesh_face_set_);                        // R3-002 — pre-existing ctor member
      }
   }
   use_precomputed_face_fluxes_ = enable;
}
```

### 6.3a Late `fault_face_set_` population in `UsePrecomputedFaceFluxes` (R5-002 FIX; supersedes R4-002's ctor-time population)

`fault_interior_faces_` and `fault_shared_faces_` are populated by the `WaveOperator` ctor at `wave_operator.inl:260-269` (interior) and `273-283` (shared); both arrays are populated in the same block of the ctor. A ctor-time `fault_face_set_` fill inserted anywhere after line 283 would see populated arrays.

The plan nevertheless moves the population into `UsePrecomputedFaceFluxes(true)` (shown in the `UsePrecomputedFaceFluxes` body above) for three reasons:
1. **Robustness against ctor-ordering drift**: a future refactor that moves fault-face population to a post-ctor helper would silently break a ctor-anchored fill. Late population is anchor-free.
2. **Fix-agent placement-error resistance**: a mechanical fix-agent could insert the R4-002 ctor loop at the wrong line (for example, after the `fault_interior_faces_.SetSize(0)` at `wave_operator.inl:260` but before the `Append` loop that populates it). Late population has a single unambiguous call site.
3. **Diagnosable failure mode**: the `MFEM_VERIFY` in the late population block prints a specific remediation pointer (`wave_operator.inl:260-269, 273-283`) if the arrays are empty at opt-in time.

**Requirement**: `fault_face_set_` must contain every mesh face index that should be routed through `FaultFaceFlux` rather than through `PrecomputedFaceFluxes`. The existing `fault_interior_faces_` / `fault_shared_faces_` arrays are the canonical source. Set semantics makes the operation idempotent on repeated `UsePrecomputedFaceFluxes(true)` calls.

**Acceptance — `test_R5_002_fault_face_set_populated_at_init_time`** (new, added to Step A acceptance; supersedes `test_R4_002_fault_faces_skipped_by_init`):
- Construct `WaveOperator` with the TPV102 fault-attr setup. Assert `fault_face_set_` is empty (ctor does not populate it).
- Call `SetFaultDOFData(...)` with the driver's usual DOFData array. Assert `fault_face_set_` still empty (this method does not populate it either).
- Call `UsePrecomputedFaceFluxes(true)`. Assert `fault_face_set_.size() == fault_interior_faces_.Size() + fault_shared_faces_.Size()`, and that every element in the set is either in `fault_interior_faces_` or `fault_shared_faces_`.
- Introspect `precomputed_face_fluxes_.face_elem_to_entry_`: assert that for every `face_idx ∈ fault_face_set_`, there is **no** key `(face_idx, any_elem)` in the map.
- Before R5-002's fix: either (a) `fault_face_set_` is empty at Init time and the map wrongly contains fault-face entries, or (b) the ctor-anchor fix-agent misplaces the loop and the set is partially populated.
- After the fix: all invariants hold.

### 6.4 Acceptance criteria — Phase 2a

- [ ] **P_SWITCH_OFF_REGRESSION**: `make test` passes with the flag unset. Default behaviour identical to before (to 1e-15 diff). Specifically: `test_adjacent_triangle_fault_uniformity`, `test_ader_tpv102_smoke`, `test_ader_linear_wave_equivalence`, `test-bp5-integration`, `test-bp5-smoke` all pass.
- [ ] **P_SWITCH_ON_LINEAR_EQ**: `test_ader_linear_wave_equivalence` extended with a `SEAS_TEST_USE_PRECOMPUTED_FLUX=1` toggle produces `||Q_new_precomputed - Q_new_runtime||_∞ < 1e-12` on a homogeneous Cartesian fixture.
- [ ] **P_SWITCH_ON_FAULT_NOT_ZEROED** (R-001 regression): Run M0 with the flag on for one ADER step; assert that fault QP samples have nonzero `I_plus_local / I_minus_local / Q_imp_plus` (i.e., the fault dispatch fired). Specifically: `||I_plus_local||_∞ > 1e-10` on the fault QP at the hypocenter.
- [ ] **P_SWITCH_ON_NO_CRASH**: `test_adjacent_triangle_fault_first_step_audit` runs to completion with `SEAS_TEST_USE_PRECOMPUTED_FLUX=1` (bit-exact behavior not yet required; just completion + finite outputs).

### 6.5 Stop / go — Phase 2a

- **GO**: all four pass → proceed to Phase 2b (MPI shared-face path).
- **STOP (default behaviour changed)**: P_SWITCH_OFF_REGRESSION fails → wiring is not truly additive. Diagnose unintended coupling.
- **STOP (R-001 regression)**: P_SWITCH_ON_FAULT_NOT_ZEROED fails → the dispatch restructure in §6.2 is wrong. Re-read §6.2's if/else order.
- **STOP (precomputed dispatch wrong)**: P_SWITCH_ON_LINEAR_EQ fails → matrix action differs from runtime. Likely Step A P_MATRIX_IDENTITY passed but the two-sided dispatch in §6.2 swaps `I_self / I_nbr` incorrectly. Re-derive.

## 7. Phase 2b — Shared (MPI) non-fault face dispatch switch (Step C-shared)

Addresses the unreviewed-area note on parallel / MPI handling.

### 7.1 Goal

After Phase 2b: `ComputeADERSharedFaceFluxRHS` also routes shared non-fault faces through `PrecomputedFaceFluxes` under the opt-in. This is required for Phase 3 P3.4 at `np≥2` to pass, because otherwise the partition-seam faces retain the runtime-normal orbit-asymmetry.

### 7.2 Implementation details

- Shared faces exist only in `ParMesh`. Per **R3-007**, `PrecomputedFaceFluxes::Init` in Phase 1 **explicitly skips** shared faces (§5.3 step 3c). Phase 2b extends Init with a ParMesh-aware second pass that pushes entries for shared non-fault faces.
- **Scope rule — Phase 2b Init extension**: for every shared non-fault face owned locally, push **exactly ONE** `FaceEntry` per rank, keyed by `(face_idx, caller_elem_id)` where `caller_elem_id` is the local cell's id (always `≥ 0`; the remote neighbor's id is the face-neighbor pseudo-id that MFEM returns as negative, which the local rank does not own and does not store). The peer rank stores its own side's entry independently.
- Shared-fault faces are already handled separately by `FaultFaceFlux`'s shared-face branch; `fault_face_set_` already includes them by construction (see `WaveOperator::SetFaultDOFData`). The Phase 2b Init extension loop skips them exactly like Phase 1.
- `ComputeADERSharedFaceFluxRHS` dispatches similarly to §6.2, with the caveat that on shared faces MFEM's `GetSharedFaceTransformations(sf)` returns only `Elem1` locally; the neighbor DOFs come via `fes.GetFaceNbrFE(nbr_id)` and `ftr->Elem2No` is the face-neighbor pseudo-id, not a real element id. The dispatch therefore calls `AddInteriorFaceRhs` **ONCE** (one-sided) with:
   - `caller_elem = ftr->Elem1No` (the local owning cell's id, non-negative)
   - `I_self` = `shape1·Q_local`, `I_nbr` = `shape2·Q_face_nbr` (ghost data)
   - writing rhs only at `dof_offset1` (the caller's own DOFs)
- This is identical to the current runtime behaviour at shared faces: `rhs -= w · shape1(i) · F_h[c]` for the local side only; the peer rank independently emits the mirror contribution on its own rhs.
- **R3-006 interaction**: the key encoding's `uint32_t` cast is required here because `ftr->Elem2No` at a shared face is a pseudo-id that MFEM reports as negative. The lookup key for the Phase 2b entry is `EncodeKey(face_idx, Elem1No)` (non-negative); the neighbor pseudo-id is never used as a key.
- **R5-001 interaction**: `ComputeADERSharedFaceFluxRHS` must receive the same hoisting pattern as §6.2.hoist (`dof_offset2`, `shape2`, `I_nbr`, `is_fault` hoisted to outer scope before the new if/else chain, with the `!is_boundary` guard — on shared faces `is_boundary` is always `false`, so the guard simplifies to always-true in Phase 2b's context, but keep the explicit guard for code-style symmetry with Phase 2a). The corresponding line ranges in the existing `ComputeADERSharedFaceFluxRHS` are distinct from Phase 2a's but the hoist structure is identical.

### 7.3 Acceptance criteria — Phase 2b

- [ ] **P_SWITCH_ON_MPI_LINEAR_EQ**: on M_ref 4×2×4 `np=4`, with the flag on, `test_ader_linear_wave_equivalence` produces `||Q_new_precomputed - Q_new_runtime||_∞ < 1e-10` (relaxed from 1e-12 to absorb MPI reduction ordering).
- [ ] **P_SWITCH_OFF_MPI_REGRESSION**: same with flag off → default path unchanged.

### 7.4 Stop / go — Phase 2b

- **GO**: both pass → proceed to Phase 3.
- **STOP**: any failure → debug the shared-face init/dispatch; do not attempt Phase 3 without both passing.

## 8. Phase 3 — Prove second-step non-fault symmetry restoration

### 8.1 Goal

After Phase 3: with `use_precomputed_face_fluxes_ = true`, the existing second-step audit shows the interior-non-fault and boundary-branch orbit-drift gates collapse to ULP, without regressing fault or volume gates.

**Primary correctness gate for the entire plan.**

### 8.2 Files to modify

- `tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp`: add `SEAS_TEST_USE_PRECOMPUTED_FLUX=1` env toggle. When set, call `wave.UsePrecomputedFaceFluxes(true)` before Gate 4. Duplicate Gates 7 / 14 as Gates 7′ / 14′ that assert the precomputed path's drift directly.

### 8.3 Acceptance criteria — Phase 3

On **M0 serial**, **M_ref 4×2×4 serial**, **M_ref np=4**, **M_big 8×4×8 np=8** with `SEAS_TEST_USE_PRECOMPUTED_FLUX=1`:

- [ ] **P3.1 fault gates unchanged**: Gates 1, 2, 3, and step-1 fault uniformity gates (I_plus/I_minus, Q_imp, F_h, fault DOF state, fault-only lift) all remain PASS at 1e-10. This is the explicit guard against the R-001 regression.
- [ ] **P3.2 volume-only unchanged**: volume-only SXY/SXZ orbit-drift at ULP (≤ 1e-14).
- [ ] **P3.3 interior non-fault → ULP**: Gate 7a equivalent under the new path → iface SXY and SXZ orbit drift both ≤ 1e-10.
- [ ] **P3.4 boundary → ULP on both SXY and SXZ**: Gate 7b/c equivalent under the new path → bface SXY and SXZ both ≤ 1e-10. This is the **critical test**. If the SXY residual (1.000 under the runtime path) does not collapse under the precomputed path, the topology-normal fix did not address it at the boundary and we stop.
- [ ] **P3.5 constant-I lifted drift**: Gate 14 equivalent under new path → iface and bface worst drift both ≤ 1e-14.
- [ ] **P3.6 full step-2 bulk signature**: existing Gate 4 final metrics (`step1_sxy_lower_sig`, etc.) → ≤ 1e-10.
- [ ] **P3.7 MPI equivalence**: `np=4` and `np=8` results agree with serial to 1e-10 on the same fault-adjacent tet sorted signatures.

### 8.4 Stop / go — Phase 3

- **GO**: all seven criteria pass on all required fixtures → proceed to Phase 4 discussion (RK path, conditional).
- **STOP (hypothesis wrong)**: P3.3 fails. Precomputed interior non-fault did NOT remove the interior symmetry break. Review: is the per-cell outward frame actually D4-covariant across orbit on the Kuhn fixture? Is there a per-element DOF-ordering issue upstream of the flux? The topology-based fix was wrong.
- **STOP (boundary residual confirmed Kuhn-level)**: P3.3 passes but P3.4 SXY stays at 1.000 and SXZ drops to ULP. Confirms the diagnosis in the debug doc: the SXY residual is a Kuhn fixture artifact at the element-DOF level, not a flux bug. This is a legitimate outcome — the plan's flux fix is complete, and the remaining SXY drift is a test-fixture concern to be handled separately (either by relaxing the 1e-10 threshold on boundary-adjacent fault-adjacent tets, or by replacing the fixture with a D4-equivariant tet mesh).
- **STOP (MPI-only regression)**: P3.3 + P3.4 pass in serial but P3.7 fails. Shared-face dispatch in Phase 2b is inconsistent with local dispatch. Debug and fix before Phase 4.

### 8.5 Fixture-artifact contingency

If Phase 3 ends with P3.4 SXY still at 1.000 (identified as a Kuhn-level residual), two **non-code** follow-ups are available:

- **Option F1**: relax the 1e-10 assertion on the boundary-adjacent fault-adjacent tets in `test_adjacent_triangle_fault_uniformity`. Document the new threshold (suggested: 1e-3 on SXY for boundary-adjacent tets) as a fixture-artefact floor, and gate only the interior fault-adjacent subset at 1e-10.
- **Option F2**: replace `MakeCartesian3D(..., TETRAHEDRON, ...)` in the audit fixture with a hand-built D4-equivariant tet decomposition. Larger effort; deferred.

Either option is a **separate change set** from this plan and should be handled after the flux fix lands.

## 9. Phase 4 — RK path migration (Step D, conditional)

### 9.1 Decision gate

Phase 4 is executed **only if** the user decides the RK4 path must also see the fix. Criterion: TPV102 production uses ADER-2, so Phase 4 is optional for TPV102. **Default**: skip unless a concrete RK consumer appears.

### 9.2 Goal (if executed)
`ComputeFaceFluxRHS` (RK path) also routes non-fault faces through `PrecomputedFaceFluxes` when the opt-in is set, with identical dispatch ordering as §6.2 (fault first, opt-in-and-boundary, opt-in-and-interior, runtime-fallback).

### 9.3 Acceptance criteria — Phase 4

- [ ] **P4.1**: one RK4 step with flag off vs on produces bit-identical `Q_new` (to 1e-12) on a homogeneous Cartesian fixture with no fault.
- [ ] **P4.2**: `test_rk4_conservation`, `test_rk4_psi_integration`, `test_rk4_endpoint_reeval` all pass with the flag off (regression) and with the flag on.

### 9.4 Stop / go — Phase 4
- **GO**: both pass → Phase 4 complete.
- **STOP**: dt-scaling mismatch between ADER (time-integrated `I`) and RK (instantaneous `Q`). Phase 1 matrices are dt-agnostic; verify at the dispatch site.

## 10. Phase 5 — Storage optimization + default-on

### 10.1 Goal
- Hash-dedup per-cell per-side matrices to reduce memory (SeisSol-like).
- Default TPV102 driver's opt-in to on behind a `--no-precomputed-flux` CLI override.

### 10.2 Acceptance criteria — Phase 5

- [ ] **P5.1**: memory after `Init` on 10⁶ cells < 100 MB per rank.
- [ ] **P5.2**: default-on VTU matches Phase 4 flag-on VTU to ULP.
- [ ] **P5.3**: Frontera init sbatch on TPV102 200 m mesh (requires user approval per project policy) completes and produces VTU qualitatively consistent with SeisSol reference.

## 11. Testing strategy summary (priority-ordered, per reviewer)

| priority | test | scope | phase |
|---|---|---|---|
| **1** | Second-step audit (existing), opt-in toggled | M0, M_ref serial + np=4, M_big np=8 | Phase 3 primary gate |
| **2** | Refined serial fixture (M_ref, M_big) | serial | Phase 3 |
| **3** | 4-rank + 8-rank MPI audit | np=4, np=8 | Phase 3 |
| **4** | Linear plane-wave equivalence | homogeneous Cartesian | Phase 2a / 2b secondary |
| 5 | P_CONSERVATION | per interior face on M0 | Phase 1 Step A |
| 6 | P_MATRIX_IDENTITY | per interior face on M_ref | Phase 1 Step A |
| 7 | P_SEISSOL_VALUE | single tet, hand-picked | Phase 1 Step A |
| 8 | P_TOPOLOGY_FRAME_MATCHES_SEISSOL | reference tet | Phase 1 Step A |
| 9 | P_BFACE_GODUNOV + P_BFACE_SEISSOL_FREE_SURFACE | per boundary face on M_ref | Phase 1 Step B |
| 10 | P_FACE2NODES_SELFCHECK | live MFEM binary | Phase 1 Step A |
| 11 | Unit-tet cell_volume regression (R2-005) | 1 reference tet | Phase 1 Step A |
| 12 | test_R009_probe_linearity | 1 reference face per BC | Phase 1 Step B |
| 13 | test_R008_fault_branch_preserved_after_dispatch_reorder | M0 serial, snapshot diff | Phase 2a gate |
| 14 | test_R3_004_neighbor_elem_is_other_side | M0, every FaceEntry | Phase 1 Step A |
| 15 | test_R3_006_key_encoding_handles_negative_caller_elem | synthetic keys | Phase 1 Step A |
| 16 | test_R3_009_p1_precondition | p=2 fixture aborts | Phase 1 Step A |
| 17 | Phase 2b entry-count consistency (M_ref np=4 under Phase 1 only) | np=4 | Phase 1 (verifies R3-007 skip) |
| 18 | test_R5_002_fault_face_set_populated_at_init_time (supersedes R4-002 test) | M0 + fault-attr setup, staged ctor/SetFaultDOFData/UsePrecomputed | Phase 1 Step A + Phase 2a |
| 19 | test_R4_003_boundary_probe_uses_topology_normal_only | 1 reference boundary face | Phase 1 Step B |
| 20 | test_R5_001_dispatch_hoist_compiles | compile Phase 2a patch | Phase 2a gate |
| 21 | Analytic free-surface reflection | constructed fixture | Phase 5 optional |
| 22 | RK equivalence | homogeneous Cartesian | Phase 4 optional |

## 12. Risk assessment

| risk | likelihood | detection | mitigation |
|---|---|---|---|
| R-001-class dispatch regression elsewhere | medium | P_SWITCH_ON_FAULT_NOT_ZEROED | explicit is_fault-first ordering in §6.2; regression test in §6.4 |
| MFEM eigenvector convention differs from SeisSol's `matR` | medium | P_SEISSOL_VALUE fails | document basis-change map before Phase 2a; flux-value equivalence still holds |
| Flux conservation broken in two-entry design | medium | P_CONSERVATION fails | single-entry restructure (major rework) — but P_CONSERVATION catches it early |
| FACE2NODES_MFEM placeholder values incorrect | medium | P_FACE2NODES_SELFCHECK aborts | reviewer fills in the table before Phase 1 code lands, per §5.1.3 |
| MFEM tet convention changes in future MFEM release | low | P_FACE2NODES_SELFCHECK aborts after upgrade | regression test catches it |
| Topology frame drift > 1e-14 | medium | P_TOPOLOGY_FRAME_MATCHES_SEISSOL fails | re-read §5.1.4 step by step; NEVER re-orthogonalize |
| Boundary Godunov probe construction bug | medium | P_BFACE_GODUNOV fails | use 9 standard basis vectors; sanity-check one column by hand |
| MPI shared-face inconsistency | medium (Phase 2b) | P3.7 or P_SWITCH_ON_MPI_LINEAR_EQ fails | Phase 2b acceptance gate before Phase 3 |
| Kuhn fixture SXY residual at 1.000 persists in P3.4 | high | P3.4 SXY > 1e-10 | legitimate outcome: invoke §8.5 (threshold relaxation or fixture swap) |
| BP5 accidentally links the new path | zero | link-time test | keep opt-in off; `seas_driver` symbol table should not contain `PrecomputedFaceFluxes` |
| R2-001 P_SEISSOL_VALUE layered formula silently swaps A⁺/A⁻ | catastrophic if test ships as-written | P_SEISSOL_VALUE trivial-case pre-check (I_self = e_j, I_nbr = 0) fails on every basis vector | **FIXED** — use direct eigendecomposition (§5.1.7); pre-check catches sign swaps before the general-case test |
| R2-002 lookup-key mismatch; function can't find its entry | high | Step A Init enumeration over M0 cells asserts `face_elem_to_entry_.count(...) == 1` for every expected pair | **FIXED** — rekeyed as `(face_idx, caller_elem_id)` (§5.1.2) |
| R2-003 fix-agent pastes `elem1`/`elem2` literally → compile error OR renames `e1`/`e2` → merge conflict | high | `grep -n "elem1\|elem2" wave_operator.inl` after patch must show 0 new matches | **FIXED** — §6.2 uses `e1`/`e2` and pre-existing-locals block explicitly listed |
| R2-004 fix-agent implements conservation with wrong sign (E2 = +E1 instead of E2 = -E1) | medium | P_CONSERVATION fails; if implementer then "fixes" an unrelated sign to mask, later tests fail downstream | **FIXED** — narrative identity corrected in §5.1.5 with derivation; P_CONSERVATION remains the test |
| R2-005 cell_volume silently garbage (integer-index subtraction) | catastrophic if shipped | unit-tet regression: `cell_volume == 1/6` to 1e-14 | **FIXED** — §5.1.4 step 10 rewrites formula using coordinates `p_elem[i]` from `mesh.GetVertex(ev[i])` |
| R2-006 fs_mode classification drift → P_SWITCH_OFF_REGRESSION fails | medium | explicit classification pseudocode in §5.3 | **FIXED** |
| R2-007 FACE2NODES_MFEM placeholder shipped | medium | `AssertMFEMFace2NodesTable()` aborts on `-1` sentinel at every Init call until reviewer fills in | **FIXED** — sentinel values + abort-on-negative |
| R2-008 fault-branch machinery dropped during dispatch restructure | high | `test_R008_fault_branch_preserved_after_dispatch_reorder` — zero diff vs pre-restructure snapshot | **FIXED** — regression test §5.1.9 |
| R2-009 probe construction assumes linearity silently | low (BCs are linear at Q_bg=0) | `test_R009_probe_linearity` — `F_h(2·I) == 2·F_h(I)` to 1e-14 | **FIXED** — precondition stated + `MFEM_VERIFY` + regression test |
| R3-001 non-existent `Geometry::GetFaceVertices` API | CRITICAL-compile | compile against MFEM — `error: no member 'GetFaceVertices'` | **FIXED** — switched to `Geometry::Constants<TETRAHEDRON>::FaceVert` compile-time constant |
| R3-002 non-existent `mesh.FaceIsBoundary` API | CRITICAL-compile | compile against MFEM | **FIXED** — switched to the `(Elem2No < 0) && (bdr_attr > 0)` predicate |
| R3-003 non-existent `mesh.GetBdrElementIndex` API | CRITICAL-compile | compile against MFEM | **FIXED** — reuse existing `WaveOperator::face_bdr_attr_` table via Init parameter |
| R3-004 `ftr->ElemNNo` typo | CRITICAL-compile | compile against MFEM | **FIXED** — `(ftr->Elem1No == e) ? ftr->Elem2No : ftr->Elem1No` |
| R3-005 stale `face_side_to_entry_` reference in §5.3 | CRITICAL-runtime | dispatch lookup failure at every interior face | **FIXED** — §5.3 step 4 rewritten to use `face_elem_to_entry_` |
| R3-006 sign-extension on negative caller_elem_id | MODERATE | `test_R3_006_key_encoding_handles_negative_caller_elem` — encode(0,-1) != encode(1, 0xFFFFFFFF) | **FIXED** — uint32_t cast in encoding |
| R3-007 Phase 1 vs Phase 2b Init scope confusion | MODERATE | Phase 1 M_ref np=4: assert `sum_entries == local_non_shared_non_fault_face_cell_side_count` | **FIXED** — explicit skip in §5.3 step 3c; one-entry-per-shared-face rule in §7.2 |
| R3-008 `const Mesh&` vs non-const API | MODERATE | compile against MFEM | **FIXED** — signature uses `Mesh &` |
| R3-009 p≥2 isoparametric breaks precomputed matrices | MODERATE | `test_R3_009_p1_precondition` — p=2 fixture aborts at Init | **FIXED** — runtime `nor_len` QP-constancy check in §5.1.6 |
| R4-001 §6.3 `Init(...)` call-site arity mismatch (6 args vs 8 required) | CRITICAL-compile | compile against MFEM — `error: too few arguments` | **FIXED** — call site expanded to pass `face_bdr_attr_` and `shared_mesh_face_set_` (pre-existing WaveOperator members) |
| R4-002 `fault_face_set_` never populated → Init treats fault faces as interior | CRITICAL-correctness (latent) | `test_R4_002_fault_faces_skipped_by_init` — introspect `face_elem_to_entry_` for fault-attr face keys | **FIXED** — eager ctor-time population from `fault_interior_faces_` ∪ `fault_shared_faces_` (§6.3a) |
| R4-003 boundary-build `(t1, t2)` dead params → maintainer trap | MODERATE | `test_R4_003_boundary_probe_uses_topology_normal_only` — perturb (t1_topo, t2_topo) and assert matrix unchanged | **FIXED** — signatures stripped; FaceEntry tangent fields marked diagnostic-only for boundary entries |
| R5-001 §6.2 dispatch references outer-scope locals that don't exist at outer scope in the current file | CRITICAL-compile | compile Phase 2a patch — `error: 'is_fault' / 'dof_offset2' / ... was not declared in this scope` | **FIXED** — §6.2.hoist adds explicit hoisting block guarded by `!is_boundary`; same pattern mirrored in Phase 2b and Phase 4 |
| R5-002 ctor-time `fault_face_set_` population brittle against ctor-ordering drift or fix-agent misplacement | CRITICAL-correctness (latent) | `test_R5_002_fault_face_set_populated_at_init_time` — assert empty at ctor, populated at opt-in, and Init creates no fault-face entries | **FIXED** — population moved into `UsePrecomputedFaceFluxes(true)` body just before `Init(...)`, with `MFEM_VERIFY` guarding empty arrays |

## 13. Out of scope (this plan)

- Fault-branch refactor (`fault_face_flux.cpp` untouched).
- `GodunovFlux::Interior / FreeSurfaceTotal / ...` body changes (only additive const accessors on `godunov_flux.hpp`).
- SeisSol's `Outflow` BC semantics (R-007 explicit exclusion).
- SeisSol's `FreeSurfaceGravity`, `Dirichlet`, `Periodic`, `Analytical` face types (R-009 explicit exclusion).
- RK path changes (Phase 4, conditional).
- Storage deduplication (Phase 5, deferred).
- Default-on behaviour (Phase 5).
- BP5 drivers or elasticity / solver / config / fault subtrees.
- MFEM library proper.
- Frontera runs (Phase 5 only, explicit user approval required).
- Kuhn-fixture replacement (§8.5 Option F2, separate change set).
- Absorbing-semantics change to SeisSol-Outflow (physics change; separate change set).

## 14. Sign-off gate

Before Phase 1 code starts, please confirm:

- [ ] The review-response map in §0 correctly captures every round-1 finding R-001 through R-010.
- [ ] The round-2 review-response map in §0′ correctly captures every round-2 finding R2-001 through R2-009.
- [ ] The semantic choice for Absorbing BC (MFEM-native `A⁺`, not SeisSol Outflow) is acceptable (§5.2.2 / R-007).
- [ ] The terminology change ("SeisSol-Aligned" → "Topology-Based ... SeisSol-inspired") is acceptable (R-002).
- [ ] The `FACE2NODES_MFEM[4][3]` in §5.1.2 holds the **verified** values `{{1,2,3},{0,3,2},{0,1,3},{0,2,1}}` from `mfem/fem/geom.cpp:987–988` (R3-001). `AssertMFEMFace2NodesTable()` compares against `mfem::Geometry::Constants<mfem::Geometry::TETRAHEDRON>::FaceVert` at every Init — not against sentinel `-1` values. No dummy-`Mesh` reviewer form needed; the values come directly from MFEM's public compile-time constant.
- [ ] Phase 2b (MPI shared-face dispatch) is accepted as a mandatory pre-Phase-3 step (unreviewed-area #2).
- [ ] Phase 3 P3.4 possibly ending with an SXY Kuhn-artefact residual (invoking §8.5) is accepted as a legitimate outcome.
- [ ] The stop/go criteria per step are acceptable as hard gates.
- [ ] P_SEISSOL_VALUE direct-eigendecomposition formulation (§5.1.7, per R2-001) is accepted as the SeisSol-flux-value-equivalence receipt, replacing the layered QgodLocal formula.
- [ ] The `face_elem_to_entry_` lookup key structure (§5.1.2, per R2-002) is accepted, with the `uint32_t` cast per R3-006.
- [ ] The fault-branch preservation regression test `test_R008_fault_branch_preserved_after_dispatch_reorder` (§5.1.9) is accepted as a mandatory Phase 2a gate.
- [ ] The round-3 response map in §0″ correctly captures every finding R3-001 through R3-009.
- [ ] The extended `Init` signature (`Mesh &`, plus `face_bdr_attr` and `shared_mesh_face_set` parameters per R3-002/R3-003/R3-008) is accepted; `WaveOperator::UsePrecomputedFaceFluxes(bool)` wires through its existing `face_bdr_attr_` and `shared_mesh_face_set_` members.
- [ ] The R3-009 p=1 isoparametric precondition (§5.1.6) is accepted; `Init` aborts on p≥2 faces via a direct numerical check on `nor_len` constancy across QPs.
- [ ] The Phase 2b shared-face scope rule (§7.2: one entry per shared face per rank, keyed by `Elem1No`) is accepted.
- [ ] The round-4 response map in §0‴ correctly captures every finding R4-001 through R4-003, including the rotation-invariance argument for isotropic elastic boundary BCs.
- [ ] The §6.3 `Init(...)` call site with all 8 arguments (R4-001) is accepted.
- [ ] The §6.3a eager ctor-time population of `fault_face_set_` (R4-002) from `fault_interior_faces_ ∪ fault_shared_faces_` is accepted as the invariant-style solution (vs a lazy/conditional fill).
- [ ] The R4-003 removal of `(t1, t2)` from `BuildBoundaryMatrices*` signatures and the diagnostic-only tagging of tangent fields on boundary `FaceEntry`s is accepted. Future maintainers are explicitly warned in §5.1.2's `FaceEntry` comment not to call `BuildRotation(n, tangent1, tangent2)` on a boundary entry's `nApNm1`.
- [ ] The round-5 response map in §0⁗ correctly captures R5-001 and R5-002.
- [ ] The §6.2.hoist explicit scope-hoisting block (R5-001) is accepted and the fix-agent will insert it **verbatim** immediately after the existing `is_boundary` computation at `wave_operator.inl:1784`.
- [ ] The same hoisting pattern is mirrored at `ComputeADERSharedFaceFluxRHS` (Phase 2b) and `ComputeFaceFluxRHS` (Phase 4, conditional) call sites.
- [ ] The §6.3a late-population policy (R5-002) is accepted: `fault_face_set_` is populated inside `UsePrecomputedFaceFluxes(true)` just before `Init(...)` — NOT in the ctor — with an `MFEM_VERIFY` that the fault-face arrays are already populated at that call site.
- [ ] The superseded R4-002 acceptance test is replaced by `test_R5_002_fault_face_set_populated_at_init_time`.

On confirmation I will start with Step A only: create the two new files, add the two additive const accessors on `godunov_flux.hpp`, and then run:
- P_FACE2NODES_SELFCHECK (against `mfem::Geometry::Constants<...>::FaceVert`)
- P_TOPOLOGY_FRAME_MATCHES_SEISSOL
- P_MATRIX_IDENTITY
- P_CONSERVATION
- P_SEISSOL_VALUE (direct-eigendecomposition form)
- P_BIMATERIAL_GUARD
- unit-tet cell_volume regression (R2-005)
- test_R009_probe_linearity
- test_R3_004_neighbor_elem_is_other_side
- test_R3_006_key_encoding_handles_negative_caller_elem
- test_R3_009_p1_precondition (p=2 fixture aborts)
- test_R5_002_fault_face_set_populated_at_init_time (supersedes test_R4_002_fault_faces_skipped_by_init)
- test_R4_003_boundary_probe_uses_topology_normal_only
- test_R5_001_dispatch_hoist_compiles (Phase 2a patch builds without undefined-variable errors)

Before running any test, perform the compile-time probe in §14 (the one-file probe asserting `Geometry::Constants<TETRAHEDRON>::FaceVert` values). This provides a single-translation-unit sanity check that the plan's §5.1.2 hard-coded table matches the MFEM installation; if it fails, halt before producing any Step A / Step B code.

Report the numbers before beginning Step B.
