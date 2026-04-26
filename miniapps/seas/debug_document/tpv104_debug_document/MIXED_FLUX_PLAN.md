# Implementation Plan: Mixed-Flux DG (Zhang et al. 2023, JGR Solid Earth)

## Overview

Zhang, Liu & Chen (2023) document **exactly the problem we have been hunting in
rounds 1–10**: the universal upwind (Godunov) flux on tetrahedral meshes
adjacent to a fault produces **Spatial Spike Oscillations (SSOs)** in the
near-fault stress field whenever the mesh is not perfectly mirror-symmetric
across the fault. Their mesh-quality metric is `r_v = max(V_A/V_B, V_B/V_A)`
between adjacent fault-face tet pairs; SSOs appear once `r_v > ~1.5`.
TPV104 on `tpv104_repro.msh` (frontal-Delaunay tet mesh) routinely has
`r_v > 3` near the rupture front — this is the source of the dip-channel
pollution that we currently mitigate by raising polynomial order
(P=2 → 174× reduction, P=3 → further).

Their fix: **mixed flux** — replace upwind with central flux at a small,
topology-defined subset of faces, dropping the unbalanced upwind dissipation
that produces the SSOs. Two variants:
- **Mixed-Flux 1**: central on every continuous boundary (interior non-fault
  + non-boundary). Eliminates SSOs but reintroduces some HFOs since 97% of
  faces are now central.
- **Mixed-Flux 2** (recommended): central ONLY on the small subset of faces
  IMMEDIATELY ADJACENT to fault elements (~5–10% of faces). Eliminates SSOs
  AND keeps HFOs negligible. This is the production option in Zhang 2023.

Goal: implement Mixed-Flux 2 (with Mixed-Flux 1 as a side-option) gated by
`--mixed-flux={none,adjacent,all-continuous}` CLI flag, default `none`
(byte-identical to current upwind path). Targets the TPV104 dynamic-rupture
dispatch only; BP5 and TPV102 production code paths are untouched.

## Constraints

- **[C-1] BP5 / TPV102 path byte-identical when flag = `none`.** This is the
  same default-OFF-with-CLI-routing pattern as round-7 R-602/R-603 substep
  iteration. No regression on `seas_test_bp5_*`, `seas_test_tpv102_*`,
  `seas_test_ader_*`, or any current `make test` target.
- **[C-2] Editable surfaces (per `feedback_dynamic_folder_editable_for_tpv104`):**
  - **Editable**: `dynamic/godunov_flux.{hpp,cpp}`,
    `dynamic/wave_operator.{hpp,inl}`, `drivers/tpv104_driver.cpp`,
    `tests/unit/test_*.cpp`, `Makefile`.
  - **No-touch**: `bp5/`, `bp1/`, `bp2/`, `domain/`, `fault/`, `solver/`,
    `friction/dieterich_ruina.hpp`, and `dynamic/fault_face_flux.cpp` should
    not be modified by this plan (they don't dispatch the bulk face flux).
- **[C-3] Existing `GodunovFlux::Interior(nor, Q_self, Q_nbr, F_h)` signature
  must remain unchanged** — it is called from many sites (RK4 `Mult`, ADER
  predictor + corrector, fault face flux assembly via `flux_.Interior(can_n,
  Q_imp_*, Q_imp_*, F_h_*)` per-side construction). New `Central(...)` is
  added alongside, not replacing.
- **[C-4] Fault-face dispatch unchanged.** Mixed-Flux 2 leaves the fault
  Riemann (`fault_flux_->Evaluate / EvaluateADER` and the per-side
  `flux_.Interior(can_n, Q_imp_*, Q_imp_*, ...)` flux re-construction)
  exactly as today. The mixed-flux change is ONLY at the
  **interior non-fault face flux** branch in `ComputeFaceFluxRHS` /
  `ComputeADERFaceFluxRHS` (and the shared-face equivalents).
- **[C-5] Numerical: central flux is `F_h = 0.5 (A_n Q⁻ + A_n Q⁺)`** —
  upwind without the dissipation term. Algebraically, since
  `Interior_upwind = 0.5 (A_n Q⁻ + A_n Q⁺) − 0.5 |A_n| (Q⁺ − Q⁻)`,
  the central flux is computed via the same `A_n`-rotation machinery,
  just dropping the `|A_n|`-weighted jump.
- **[C-6] No SeisSol-specific ANTI-feature creep.** This plan implements the
  Zhang 2023 method as published. Other variants (artificial damping,
  Drucker-Prager regularization) are out of scope.

## Mathematical reference (from Zhang 2023)

Velocity-strain form of elastodynamics (Eq. 3 in paper):
```
∂q/∂t = ∂f(q)/∂x + ∂g(q)/∂y + ∂h(q)/∂z
q = (ρv_x, ρv_y, ρv_z, ε_xx, ε_yy, ε_zz, γ_yz, γ_xz, γ_xy)ᵀ
```

Our codebase uses velocity-stress (NUM_STATE = 9 with components
SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ). The mathematics translates
directly: replace ε-based fluxes with σ-based fluxes; the central flux
formula `F_h = 0.5 (A_n Q⁻ + A_n Q⁺)` is form-invariant.

DG weak form (Eq. 7 in paper):
```
∫_Ω ψ_k [ ∂_t q − ∇·(f, g, h) ] dV = ∫_∂Ω ψ_k T⁻¹ (F̂ − f) dS
```

Numerical flux F̂ on face with normal n:
- **Upwind** (current code, GodunovFlux::Interior): characteristic Riemann
  solution at face = `0.5(A_n Q⁻ + A_n Q⁺) − 0.5|A_n|(Q⁺ − Q⁻)`
- **Central** (Mixed-Flux 1, Eq 19): `F̂ = 0.5(f(q) + f(q⁺)) =
  0.5 A_n (Q⁻ + Q⁺)` (same as upwind without dissipation term)
- **Mixed-Flux 2**: central ONLY on `e ∈ E_fault_adj`, upwind elsewhere.

Topological definition of `E_fault_adj` (Figure 5 in paper):
```
E_fault_adj = { face f : f is interior non-fault AND
                          at least one element of {Elem1(f), Elem2(f)}
                          contains another face that IS a fault face }
```

For TPV104's `tpv104_repro.msh` with ~3,840 tets and ~10⁴ fault QPs,
`|E_fault_adj| ≈ 6 × |fault_faces|` (each fault tet has 3 non-fault faces).
The set is small (< 10% of all interior faces).

---

## Phase 1: Add `GodunovFlux::Central` primitive

### Goal
A new public method on `GodunovFlux` that computes the central (arithmetic-
average) numerical flux for the elastodynamic velocity-stress system,
producing the same output buffer layout as `Interior(...)` but without the
`|A_n|` dissipation term. Existing `Interior` callers untouched.

### Files to Modify
- `dynamic/godunov_flux.hpp` — declare `Central(...)` with the same
  signature as `Interior(...)`.
- `dynamic/godunov_flux.cpp` — implement `Central(...)`.

### Detailed Requirements

1. **Function signature** (mirror of `Interior`):
   ```cpp
   void Central(const real_t nor[3],
                const real_t *Q_self,
                const real_t *Q_nbr,
                real_t *F_h) const;
   ```
   - `nor`: unit face normal in global coordinates (3 components).
   - `Q_self`, `Q_nbr`: NUM_STATE = 9 components in global frame.
   - `F_h`: NUM_STATE output buffer, written (not added).

2. **Implementation strategy** (R-1201 corrected). Mirror `Interior`'s
   5-step pattern (`godunov_flux.cpp:350–374`) using the SAME
   eigenvalue-split private members `Ax_plus_` and `Ax_minus_` already
   accessible to `GodunovFlux::Central` as a member function:

   - Step 1: `BuildFrame(nor, t1, t2)` — same as Interior.
   - Step 2: `BuildRotationInverse(nor, t1, t2, Tinv)`,
             `BuildRotation(nor, t1, t2, T)` — same as Interior.
   - Step 3: `Tinv.Mult(Q_self, Q_self_rot)`,
             `Tinv.Mult(Q_nbr,  Q_nbr_rot)` — same as Interior.
   - Step 4 (NEW for Central): central flux in face-rotated frame.
     ```cpp
     // F_central_rot = 0.5 * (Ax_plus_+Ax_minus_) * (Q_self_rot+Q_nbr_rot)
     // Equivalently:
     //   F_central_rot = 0.5 * (Ax_plus_·Q_self_rot + Ax_plus_·Q_nbr_rot
     //                        + Ax_minus_·Q_self_rot + Ax_minus_·Q_nbr_rot)
     real_t Q_sum[NUM_STATE];
     for (int i = 0; i < NUM_STATE; i++)
        { Q_sum[i] = Q_self_rot[i] + Q_nbr_rot[i]; }
     real_t F_rot[NUM_STATE];
     for (int i = 0; i < NUM_STATE; i++) {
        real_t s = 0.0;
        for (int j = 0; j < NUM_STATE; j++)
           s += (Ax_plus_(i, j) + Ax_minus_(i, j)) * Q_sum[j];
        F_rot[i] = 0.5 * s;
     }
     ```
   - Step 5: `T.Mult(F_rot, F_h)` — same as Interior.

   `Ax_plus_` and `Ax_minus_` are private members of `GodunovFlux`
   (declared at `godunov_flux.hpp:211–212`); the new `Central` method
   has access since it is a member.  **DO NOT call
   `GetReferenceStarMatrix(0)`** — that returns the GLOBAL-frame
   Jacobian `A_x` (used by the ADER predictor's CK recursion), which is
   a numerically distinct matrix from the face-rotated `Ax_plus_ +
   Ax_minus_`.  Mixing them silently produces wrong flux values.

   Sign convention for the algebraic identity (used by Phase 5 R-1101
   test): in the face-rotated frame,
   ```
   |A_n|_rot = Ax_plus_ - Ax_minus_                  (positive semi-def.)
   F_upwind_rot - F_central_rot = 0.5 · |A_n|_rot · (Q_self_rot − Q_nbr_rot)
   ```
   where `Q_self` is the "self" (Elem1-side) state passed by the caller.

3. **Verification by algebraic relation**: at every QP, must have
   `Interior(nor, Q_self, Q_nbr, F_up) - Central(nor, Q_self, Q_nbr, F_ce) =
   −0.5 |A_n|(Q_nbr − Q_self)`. The unit test in Phase 5 verifies this
   numerically by comparing both calls on a controlled (Q_self, Q_nbr) pair
   with non-trivial jump.

4. **Edge cases**:
   - `Q_self == Q_nbr`: jump is zero, so `Central == Interior` exactly.
     Test must verify bit-equality in this case.
   - `nor` not unit length: callers normalize before calling. `Central` does
     NOT re-normalize; it assumes the same convention as `Interior`.

### Interfaces
```cpp
// godunov_flux.hpp, after the existing Interior decl
void Central(const real_t nor[3],
             const real_t *Q_self,
             const real_t *Q_nbr,
             real_t *F_h) const;
```

### Edge Cases to Handle
- **Identical states** (`Q_self == Q_nbr`): `Central == Interior` bit-equal.
- **Coordinate-axis-aligned `nor`** (e.g., `(1, 0, 0)`): `Tinv = T = I_block`
  for the appropriate channel permutation; central reduces to
  `F_h[c] = 0.5 · (A_x Q⁻ + A_x Q⁺)[c]` directly.
- **`nor = (0, −1, 0)`** (TPV104 fault-normal direction in global frame):
  routine returns the fault-perpendicular central flux. For non-fault
  faces (where mixed-flux is applied), `nor` is whatever CalcOrtho gave;
  the routine is normal-direction-agnostic.

### Acceptance Criteria
- [ ] `Central(nor, Q, Q, F_h)` produces the same output as
      `Interior(nor, Q, Q, F_h)` to FP precision (tested with a sweep of
      9 channel-isolated Q states and 8 representative `nor` directions).
- [ ] **R-1204 sign convention**: `Interior(nor, Q_self, Q_nbr, F_up)
      − Central(nor, Q_self, Q_nbr, F_ce)` equals
      `+0.5·|A_n|·(Q_self − Q_nbr)` for at least 3 distinct
      `(Q_self, Q_nbr)` pairs and the same 8 normals. Tolerance: 1e−12
      absolute.  `|A_n|` here is `T · (Ax_plus_ − Ax_minus_) · Tinv` — the
      face-rotated absolute-Jacobian, rotated back to global; the test
      computes it explicitly and rotates rather than reading any private
      member that doesn't exist.  `Q_self` is the state of the element
      whose RHS contribution the flux serves (Elem1's side).
- [ ] Existing `seas_test_godunov_flux*` (if any) and ADER tests pass
      bit-identically.

### Dependencies
- Depends on: nothing (uses only existing `BuildRotation`,
  `BuildRotationInverse`, `GetReferenceStarMatrix` machinery).
- Required by: Phase 4.

---

## Phase 2: CLI flag + `WaveOperator` config

### Goal
Driver and wave-operator are aware of a runtime mixed-flux mode. When mode
is `none`, behavior is byte-identical to current path.

### Files to Modify
- `dynamic/wave_operator.hpp` — add `MixedFluxMode` enum, setter,
  getter, and a private member.
- `dynamic/wave_operator.inl` — implement setter; default-constructor
  initializes mode to `None`.
- `drivers/tpv104_driver.cpp` — parse `--mixed-flux <value>` CLI flag,
  validate, call `wave.SetMixedFluxMode(...)`, echo to banner.

### Detailed Requirements

1. **Enum** in `wave_operator.hpp` after the existing `FreeSurfaceBCMode`
   declaration:
   ```cpp
   /// Zhang et al. 2023 mixed-flux mode for interior non-fault face flux.
   /// - None:           upwind everywhere (DEFAULT, byte-identical to
   ///                   pre-Mixed-Flux behavior).
   /// - Adjacent:       central flux on faces immediately adjacent to a
   ///                   fault element, upwind elsewhere (Mixed-Flux 2,
   ///                   recommended).
   /// - AllContinuous:  central flux on every interior non-fault face
   ///                   (Mixed-Flux 1, eliminates SSOs but allows minor
   ///                   HFOs).
   enum class MixedFluxMode : int { None = 0, Adjacent = 1, AllContinuous = 2 };
   ```

2. **Setter / getter**:
   ```cpp
   void SetMixedFluxMode(MixedFluxMode m);
   MixedFluxMode GetMixedFluxMode() const { return mixed_flux_mode_; }
   ```

3. **Private member** in the existing private block:
   ```cpp
   MixedFluxMode mixed_flux_mode_ = MixedFluxMode::None;
   ```

4. **Setter implementation** in `wave_operator.inl`. When mode flips from
   `None` to `Adjacent` or `AllContinuous`, populate
   `central_flux_face_set_` (Phase 3). When mode flips back to `None`,
   clear the set. The setter is non-const (modifies state).

5. **R-1203 cross-check with precomputed-flux state.** The wave operator
   has TWO setters that affect the same dispatch decision —
   `SetMixedFluxMode` and `UsePrecomputedFaceFluxes`.  Both must abort if
   the OTHER is non-default, regardless of call order.  At the start of
   `SetMixedFluxMode(MixedFluxMode m)`:
   ```cpp
   if (m != MixedFluxMode::None && use_precomputed_face_fluxes_) {
      MFEM_ABORT("SetMixedFluxMode("
                 << (m == MixedFluxMode::Adjacent ? "Adjacent"
                                                  : "AllContinuous")
                 << "): mutually exclusive with precomputed face flux "
                 "(UsePrecomputedFaceFluxes is currently true).  "
                 "Disable precomputed flux first.");
   }
   ```
   And SYMMETRICALLY at the start of
   `UsePrecomputedFaceFluxes(bool enable)` (existing method, modified
   under [C-2] editable scope):
   ```cpp
   if (enable && mixed_flux_mode_ != MixedFluxMode::None) {
      MFEM_ABORT("UsePrecomputedFaceFluxes(true): mutually exclusive "
                 "with mixed-flux mode (currently "
                 << static_cast<int>(mixed_flux_mode_)
                 << ").  Disable mixed flux first via "
                 "SetMixedFluxMode(MixedFluxMode::None).");
   }
   ```
   The driver-level abort in Phase 6 becomes a defense-in-depth layer
   that fails fast BEFORE mesh construction; the wave-operator setters
   are the canonical guards.

5. **Driver CLI**:
   ```cpp
   std::string mixed_flux_str =
      GetStringArg(argc, argv, "--mixed-flux", "none");
   WaveOperator<MeshT>::MixedFluxMode mixed_flux_mode;
   if      (mixed_flux_str == "none")           { mixed_flux_mode = ...::None; }
   else if (mixed_flux_str == "adjacent")       { mixed_flux_mode = ...::Adjacent; }
   else if (mixed_flux_str == "all-continuous") { mixed_flux_mode = ...::AllContinuous; }
   else { MFEM_ABORT("--mixed-flux: unknown value '" << mixed_flux_str
                     << "'.  Accepted: none | adjacent | all-continuous."); }
   ```
   - Banner echoes the chosen mode, e.g.,
     `Mixed flux: adjacent (Mixed-Flux 2 per Zhang et al. 2023, central on fault-adjacent non-fault interior faces)`.
   - The `wave.SetMixedFluxMode(mixed_flux_mode)` call goes AFTER
     `wave.SetFaultFlux(...)` (so the fault face geometry is populated
     before the central-flux face set is built).

### Interfaces
```cpp
// wave_operator.hpp public:
enum class MixedFluxMode : int { None = 0, Adjacent = 1, AllContinuous = 2 };
void SetMixedFluxMode(MixedFluxMode m);
MixedFluxMode GetMixedFluxMode() const;
const std::set<int> &GetCentralFluxFaceSet() const;  // test-only, see Phase 3
```

### Edge Cases to Handle
- **Setter called before fault attribute is set**:
  `Adjacent` mode requires `bc_.fault_attr > 0` and a populated
  `fault_interior_faces_`. If either is missing,
  `MFEM_ABORT("SetMixedFluxMode(Adjacent): fault attribute / interior fault faces not yet populated; ...");`.
- **Setter called twice with same value**: idempotent — no-op.
- **CLI value with whitespace / wrong case**: only the exact strings
  `none | adjacent | all-continuous` accepted; anything else aborts.

### Acceptance Criteria
- [ ] Default `--mixed-flux none` (or omitted): driver banner shows
      `Mixed flux: none (upwind everywhere, default)`.
- [ ] `--mixed-flux adjacent`: banner shows the adjacent description; no
      regression on `make test` (which doesn't pass the flag).
- [ ] `--mixed-flux foobar`: aborts with the exact error message above.
- [ ] `wave.GetMixedFluxMode()` reports the configured value.

### Dependencies
- Depends on: nothing (parsing + storage only).
- Required by: Phase 3, Phase 4, Phase 6.

---

## Phase 3: Build the central-flux face set

### Goal
At the moment `SetMixedFluxMode(Adjacent)` is called, identify the set of
interior non-fault faces that touch a fault-adjacent element. Store as
`std::set<int> central_flux_face_set_` keyed by mesh face index.

### Files to Modify
- `dynamic/wave_operator.hpp` — declare private member
  `std::unordered_set<int> central_flux_face_set_;` (R-1206: O(1)
  lookup, not O(log N) — load-bearing for performance at production
  scale where the lookup is in the hot per-face per-step inner loop).
- `dynamic/wave_operator.inl` — implement
  `BuildCentralFluxFaceSet_()` private helper, called from the setter.

### Detailed Requirements

1. **Set construction algorithm** (matches Zhang Figure 5 / Section 3.2,
   with R-1207 explicit MFEM membership conditions):
   ```
   Step 1: Build E_fault_adj = set of element indices e such that some face
           of e is in fault_interior_faces_ ∪ fault_shared_faces_
           (i.e., e is on either side of a fault face).

   Step 2: For each element e ∈ E_fault_adj:
              Array<int> faces, ori;
              mesh.GetElementFaces(e, faces, ori);
              For each face f in faces:
                 // Two-pronged "is interior" check (R-1207):
                 // (a) regular two-sided: GetFaceElementTransformations
                 //     returns ftr with Elem2No >= 0
                 // (b) MPI shared seam: f appears in shared_mesh_face_set_
                 //     even though Elem2No < 0 (neighbor on another rank)
                 FaceElementTransformations *ftr =
                    mesh.GetFaceElementTransformations(f);
                 const bool two_sided_interior =
                    (ftr != nullptr && ftr->Elem2No >= 0);
                 const bool shared_seam =
                    shared_mesh_face_set_.count(f) > 0;
                 const bool is_interior = two_sided_interior || shared_seam;

                 // "Is fault" check: bdr-attribute table populated at ctor.
                 const bool is_fault =
                    (bc_.fault_attr > 0 &&
                     face_bdr_attr_[f] == bc_.fault_attr);

                 if (is_interior && !is_fault) {
                    central_flux_face_set_.insert(f);
                 }
   ```
   For TPV104's tet mesh: each fault element has 1 fault face and 3
   non-fault faces. The total |E_fault_adj| ≈ 2 × N_fault_faces; the total
   `|central_flux_face_set_|` ≈ 6 × N_fault_faces (each face shared by 2
   tets, so unique-set is ~3 × N_fault_faces).

2. **Mode `AllContinuous`**: instead of the algorithm above, populate the
   set with EVERY interior non-fault face:
   ```
   for f in [0, mesh.GetNumFaces()):
      ftr = mesh.GetFaceElementTransformations(f)
      if ftr and ftr.Elem2No >= 0 and face_bdr_attr_[f] != bc_.fault_attr
         and shared_mesh_face_set_.count(f) == 0:
         central_flux_face_set_.insert(f)
   ```

3. **MPI / shared faces**: include shared interior faces (those tracked by
   `shared_mesh_face_set_`) in the candidate iteration for `AllContinuous`.
   For `Adjacent`, a shared face is included iff one of its (local +
   ghost) elements is in `E_fault_adj`. The shared-face per-rank
   identification uses `pmesh.GetSharedFace(sf)` and the existing
   `shared_face_bdr_attr_[sf]` / `fault_shared_faces_` infrastructure.

4. **Logging**: at end of build, if rank 0, print:
   ```
   [mixed-flux] mode=adjacent: |fault_interior_faces|=N_F_int,
                |fault_shared_faces|=N_F_shr,
                |E_fault_adj|=N_E,
                |central_flux_face_set|=N_C
                ({pct}% of total interior non-fault faces).
   ```

5. **Test hook**: public `GetCentralFluxFaceSet() const` returns
   `const std::set<int>&` for unit-test assertions.

### Interfaces
```cpp
// wave_operator.hpp private:
std::set<int> central_flux_face_set_;
void BuildCentralFluxFaceSet_();   // called from SetMixedFluxMode

// wave_operator.hpp public:
const std::set<int> &GetCentralFluxFaceSet() const
{ return central_flux_face_set_; }
```

### Edge Cases to Handle
- **No fault attribute** (`bc_.fault_attr == 0`): `Adjacent` mode aborts at
  the setter (Phase 2 edge case). `AllContinuous` mode just iterates all
  interior faces and produces a populated set.
- **Mesh with zero fault elements**: same as above; `Adjacent` aborts.
- **Re-build after toggling mode**: setter clears the set first, then
  rebuilds. Idempotent on repeated identical mode setting (clears and
  rebuilds an identical set; bit-equivalent to no-op for std::set).

### Acceptance Criteria
- [ ] On `tpv104_repro.msh` (3,840 tets) with `--mixed-flux adjacent`:
      `|central_flux_face_set_|` is between 3 × N_fault_faces and 7 ×
      N_fault_faces (sanity bounds; exact count depends on per-tet face
      sharing).
- [ ] On a 2-tet fault fixture (test mesh): `|central_flux_face_set_| = 6`
      (each tet has 4 faces, 1 is fault; remaining 3 per tet × 2 tets = 6
      non-fault faces, no shared faces).
- [ ] `AllContinuous` populates strictly more faces than `Adjacent` on the
      same mesh.
- [ ] No fault face is ever in the central-flux set (test asserts
      `central_flux_face_set_ ∩ fault_interior_faces_ == ∅`).

### Dependencies
- Depends on: Phase 2.
- Required by: Phase 4.

---

## Phase 4: Wire mixed-flux dispatch into face-flux loops

### Goal
At every interior non-fault face flux site in `wave_operator.inl`, check
whether the face index is in `central_flux_face_set_`. If yes, call
`flux_.Central(...)` instead of `flux_.Interior(...)`. Bit-identical to
current path when `central_flux_face_set_` is empty (mode = `None`).

### Files to Modify
- `dynamic/wave_operator.inl` — four sites:
  - `ComputeFaceFluxRHS` interior non-fault else-branch (RK4 path).
  - `ComputeADERFaceFluxRHS` interior non-fault else-branch (ADER path).
  - `ComputeSharedFaceFluxRHS` interior non-fault branch (RK4 MPI).
  - `ComputeADERSharedFaceFluxRHS` interior non-fault branch (ADER MPI).

### Detailed Requirements

1. **Per-site dispatch pattern**. Find each existing call:
   ```cpp
   flux_.Interior(nor, Q_self, Q_nbr, F_h);
   ```
   Replace with:
   ```cpp
   if (central_flux_face_set_.count(f) > 0)
   {
      flux_.Central(nor, Q_self, Q_nbr, F_h);
   }
   else
   {
      flux_.Interior(nor, Q_self, Q_nbr, F_h);
   }
   ```
   The `f` (mesh face index) is in scope at every site (it's the face
   loop variable). For shared-face sites, `f` is obtained via
   `pmesh.GetSharedFace(sf)`.

2. **Bit-equivalence guard + R-1206 short-circuit**. The dispatch must
   short-circuit on the mode flag so that default-mode runs incur ZERO
   per-face overhead (not just zero behavior change).  Hoist a boolean
   ONCE per `ComputeFaceFluxRHS` / `ComputeADERFaceFluxRHS` call (NOT
   per face):
   ```cpp
   const bool mf_on = (mixed_flux_mode_ != MixedFluxMode::None);
   ```
   Then at every dispatch site:
   ```cpp
   if (use_precomputed_face_fluxes_) { /* precomputed path */ }
   else if (mf_on && central_flux_face_set_.count(f) > 0)
   { flux_.Central(nor, Q_self, Q_nbr, F_h); }
   else
   { flux_.Interior(nor, Q_self, Q_nbr, F_h); }
   ```
   With `mf_on == false`, the `&&` short-circuits and `count(f)` never
   executes — bit-identical AND zero-cost vs pre-change.  At production
   scale (~1.5M faces × 7000 steps), an unconditional set lookup would
   add ~6 hours of wall-clock per simulation even at empty set; the
   short-circuit eliminates this.

3. **Precomputed face flux compatibility**. If
   `use_precomputed_face_fluxes_` is true (TPV102 Phase 2a path), the
   precomputed flux assumes UPWIND. Do NOT route through Mixed-Flux when
   precomputed path is active. Add a guard:
   ```cpp
   if (use_precomputed_face_fluxes_)
   {
      precomputed_face_fluxes_.AddInteriorFaceRhs(...);   // unchanged
   }
   else if (central_flux_face_set_.count(f) > 0)
   {
      flux_.Central(nor, Q_self, Q_nbr, F_h);
      // existing F_h → rhs accumulation continues unchanged
   }
   else
   {
      flux_.Interior(nor, Q_self, Q_nbr, F_h);
      // existing F_h → rhs accumulation
   }
   ```
   Document at the site that `--mixed-flux` and the precomputed-flux flag
   are MUTUALLY EXCLUSIVE; abort at driver level if both are requested
   simultaneously.

4. **PML interaction**. PML damping is applied AFTER the face-flux loop and
   is independent of the flux choice. No change to PML.

5. **Free-surface / absorbing BC branches**. UNTOUCHED. Mixed-flux applies
   only to INTERIOR NON-FAULT faces per Zhang Section 3.2. Boundary fluxes
   (`bdr_attr > 0` cases) keep their existing dispatch.

### Interfaces
No new public interfaces. Modifications are entirely within existing
function bodies.

### Edge Cases to Handle
- **Empty face set + flag != none**: setter aborts (Phase 3 edge case).
- **Face appears in BOTH central set AND fault set**: cannot happen by the
  Phase 3 construction (`central_flux_face_set_` excludes fault faces);
  add an `MFEM_ASSERT` in the dispatch sites for defense-in-depth.
- **Mutual exclusion with precomputed fluxes**: driver MUST abort if
  `--mixed-flux != none` AND precomputed-flux flag is set.

### Acceptance Criteria
- [ ] With `--mixed-flux none`: bit-identical Q_new after one
      `AdvanceADER` step on the 2-tet fault fixture vs current code.
- [ ] With `--mixed-flux adjacent` on the 2-tet fixture: Q_new differs
      from `--mixed-flux none` only at the 6 non-fault faces; the
      difference is `0.5 |A_n|(Q⁺ − Q⁻)` per face per state component
      (algebraic relation per Phase 1 acceptance criterion).
- [ ] All existing TPV102 / ADER / BP5 / friction tests pass (default
      flag = `none`).

### Dependencies
- Depends on: Phase 1 (`Central`), Phase 3 (set built).
- Required by: Phase 5 (test), Phase 6 (driver wiring).

---

## Phase 5: Unit tests

### Goal
Three new unit tests guard the contract:

1. **R-1101**: Central-flux primitive correctness
   (`test_godunov_central_flux.cpp`).
2. **R-1102**: `central_flux_face_set_` topology correctness
   (`test_mixed_flux_face_set.cpp`).
3. **R-1103**: Mixed-flux dispatch bit-identity at flag = `none`
   (`test_mixed_flux_dispatch_none.cpp`).

### Files to Create
- `tests/unit/test_godunov_central_flux.cpp`
- `tests/unit/test_mixed_flux_face_set.cpp`
- `tests/unit/test_mixed_flux_dispatch_none.cpp`

### Files to Modify
- `Makefile` — add `TEST_GODUNOV_CENTRAL_FLUX_*`,
  `TEST_MIXED_FLUX_FACE_SET_*`, `TEST_MIXED_FLUX_DISPATCH_NONE_*` macros
  + 3 build targets following the existing pattern (e.g.,
  `seas_test_tpv104_substep_predictor` recipe).

### Detailed Requirements

#### R-1101: `test_godunov_central_flux.cpp`
Constructs a `GodunovFlux` with TPV104 material parameters
(`rho=2670, cp=6000, cs=3464`). For 8 face normals
(±x, ±y, ±z, plus 2 oblique) and a sweep of (Q_self, Q_nbr) pairs:

| gate | check |
|---|---|
| 1 | `Central(n, Q, Q, F_c) == Interior(n, Q, Q, F_u)` bit-equal |
| 2 | `Interior(n, Q⁻, Q⁺, F_u) − Central(n, Q⁻, Q⁺, F_c) = 0.5 |A_n|(Q⁺ − Q⁻)` to 1e−12 |
| 3 | `Central(n, Q⁻, Q⁺, F_c)` is invariant under simultaneous swap of
     `Q⁻ ↔ Q⁺` and `n ↔ −n` (central flux's anti-symmetry) |

#### R-1102: `test_mixed_flux_face_set.cpp`
**R-1202 corrected**: the 2-tet fixture from `test_ader_tpv102_smoke` has
ZERO interior non-fault faces (every non-fault face on those 2 tets is a
boundary face), so `|central_flux_face_set_| == 0` on it — useful as a
sub-gate but cannot exercise the main contract. The main test must use a
multi-tet fixture with interior non-fault faces.

Primary fixture: `Mesh::MakeCartesian3D(2, 2, 2, Element::TETRAHEDRON,
Lx, Ly, Lz)` — 24 tets (8 hexes × 6-tet split), with interior faces
between cells. Tag the y = Ly/2 plane (interior) as a fault face by
post-processing the boundary-attribute table OR by constructing the mesh
manually with the desired bdr triangle attributes.

Configure `bc_.fault_attr = 3`. Set `Adjacent` mode. Assert:
- `|central_flux_face_set_|` is a strictly POSITIVE integer.
- For every face f in `central_flux_face_set_`:
   - `mesh.GetFaceElementTransformations(f)->Elem2No >= 0` OR
     `shared_mesh_face_set_.count(f) > 0` (i.e., f is interior or a
     shared seam — see R-1207 below).
   - `face_bdr_attr_[f] != bc_.fault_attr` (f is not a fault face).
- For every face f in `central_flux_face_set_`, at least one of
  `(Elem1No, Elem2No)` is in `E_fault_adj` (the construction's defining
  property). This guards against the off-by-one in Phase 3.
- `central_flux_face_set_ ∩ fault_interior_faces_ == ∅`.
- Re-set `None` mode → set is empty (size 0).
- Re-set `AllContinuous` mode → set strictly larger than `Adjacent` (every
  interior non-fault face is included).

Secondary sub-gate (sanity check on the 2-tet fixture): assert
`|central_flux_face_set_| == 0` on the existing 2-tet
`BuildTwoTetFaultMesh` (because all non-fault faces are boundary faces).
Catches the implementer accidentally including boundary faces in the set.

#### R-1103: `test_mixed_flux_dispatch_none.cpp`
2-tet TPV102 fixture. Initialize Q with TPV102 pre-stress.
- Path A: `wave_a.AdvanceADER(Q, dt, 2, Q_new_a)` with mode = None.
- Path B: `wave_b.AdvanceADER(Q, dt, 2, Q_new_b)` after `wave_b.SetMixedFluxMode(None)` (idempotent set).
Assert `Q_new_a == Q_new_b` bit-equal.

### Interfaces
Tests follow the existing `TEST_LE` / `TEST_ASSERT` macro pattern from
`test_ader_linear_wave_equivalence.cpp`. No new infrastructure.

### Edge Cases to Handle
- Test 1 covers `Q_self == Q_nbr` (bit-equality of upwind and central
  when jump is zero).
- Test 2 covers the `bc_.fault_attr = 0` abort (call `SetMixedFluxMode
  (Adjacent)` on a no-fault mesh; expect `MFEM_VERIFY` abort).

### Acceptance Criteria
- [ ] All three tests build via `make seas_test_godunov_central_flux
      seas_test_mixed_flux_face_set seas_test_mixed_flux_dispatch_none`.
- [ ] All three tests run and pass (each `Results: N passed, 0 failed`).
- [ ] Existing test suite (`make test` for TPV104/ADER/TPV102/BP5/friction)
      passes unchanged.

### Dependencies
- Depends on: Phases 1, 2, 3, 4.
- Required by: Phase 6 (verifies before driver wiring).

---

## Phase 6: Driver wiring + smoke test

### Goal
TPV104 driver accepts `--mixed-flux=adjacent`, dispatches the wave operator
through the mixed-flux interior face flux, and produces a complete
production run (no crash, finite slip, sensible peak slip rate).

### Files to Modify
- `drivers/tpv104_driver.cpp` — CLI parse + setter + banner update +
  precomputed-flux mutual-exclusion guard + smoke output.
- `tests/unit/test_tpv104_smoke.cpp` — extend the banner-text gate to
  cover the `Mixed flux:` line.

### Detailed Requirements

0. **R-1205 required call order in the driver** (must be documented in
   `drivers/tpv104_driver.cpp` near the wave-operator construction and
   the iterator-setup region):
   ```
   1. WaveOperator<MeshT> wave(pmesh, order, lambda, mu, rho, bc);
      // bc.fault_attr already set; ctor populates fault_interior_faces_.
   2. wave.SetFaultFlux(&fault_flux);
   3. wave.SetFaultDOFData(&dof_data, nqp_per_face);
   4. wave.SetAbsorbingBackground(Q_bg);
   5. wave.SetMixedFluxMode(mixed_flux_mode);   // <-- AFTER #1; safe
                                                //     because the ctor
                                                //     has populated
                                                //     fault_interior_faces_
                                                //     and bc_.fault_attr
                                                //     is in scope.
   ```
   The setter cross-checks (`bc_.fault_attr > 0` for `Adjacent` /
   `AllContinuous` modes; precomputed-flux exclusion per R-1203) and
   aborts on failure.  For drivers that don't have a fault (e.g.,
   pure-bulk wave-equation tests), only `MixedFluxMode::None` is
   allowed — calling Adjacent on a no-fault mesh aborts.

1. **CLI parse** (Phase 2 already has the parse; this phase wires it into
   the actual `wave.SetMixedFluxMode(...)` call):
   ```cpp
   // After wave.SetFaultFlux(&fault_flux); wave.SetFaultDOFData(...);
   wave.SetMixedFluxMode(mixed_flux_mode);
   if (rank == 0 && mixed_flux_mode != WaveOperator<MeshT>::MixedFluxMode::None)
   {
      std::cout << "[mixed-flux] mode=" << mixed_flux_str
                << "  |central_set|="
                << wave.GetCentralFluxFaceSet().size()
                << "  (Zhang et al. 2023 mixed-flux dispatch)\n";
   }
   ```

2. **Mutual-exclusion guard with precomputed-flux**:
   ```cpp
   const bool use_precomputed_flux =
      HasFlag(argc, argv, "--use-precomputed-face-fluxes");
   if (use_precomputed_flux && mixed_flux_mode !=
       WaveOperator<MeshT>::MixedFluxMode::None)
   {
      MFEM_ABORT(
         "--use-precomputed-face-fluxes is mutually exclusive with "
         "--mixed-flux != none.  The precomputed-flux path assumes "
         "upwind dispatch.");
   }
   ```

3. **Banner update**:
   ```cpp
   const char *mixed_flux_banner =
      (mixed_flux_mode == ::None)          ? "none (upwind everywhere, default)"
    : (mixed_flux_mode == ::Adjacent)      ? "adjacent (Mixed-Flux 2 per Zhang et al. 2023, central on fault-adjacent non-fault interior faces)"
    : (mixed_flux_mode == ::AllContinuous) ? "all-continuous (Mixed-Flux 1, central on every interior non-fault face)"
    : "?";
   std::cout << "Mixed flux: " << mixed_flux_banner << "\n";
   ```

4. **Smoke test extension**:
   `tests/unit/test_tpv104_smoke.cpp` adds:
   - default flag (no `--mixed-flux`): banner contains `Mixed flux: none (upwind everywhere, default)`.
   - `--mixed-flux adjacent`: banner contains
     `Mixed flux: adjacent (Mixed-Flux 2 per Zhang et al. 2023`.
   - `--mixed-flux foobar`: driver aborts.

### Interfaces
No new functions. Driver updates only.

### Edge Cases to Handle
- **--dry-run + --mixed-flux adjacent**: the dry-run path returns before
  populating the central-flux set. Either (a) call
  `wave.SetMixedFluxMode(...)` BEFORE the dry-run check, OR (b) do not
  exercise mixed-flux on dry-run (skip the setter). Recommend (b) — the
  banner echoes the requested mode, but the set is not built; explicit
  comment in the driver.
- **--fault-iterator substep + --mixed-flux adjacent**: should be valid
  (substep dispatch + mixed flux are independent paths). No mutual-
  exclusion required.

### Acceptance Criteria
- [ ] `./seas_tpv104_driver --dry-run` (no flag): banner shows `Mixed flux: none`.
- [ ] `./seas_tpv104_driver --dry-run --mixed-flux adjacent`: banner shows
      adjacent line.
- [ ] `./seas_tpv104_driver --dry-run --mixed-flux foobar`: aborts.
- [ ] `./seas_tpv104_driver --tfinal 0.1 --mesh tpv104_repro.msh
      --mixed-flux adjacent`: completes without crash, slip_strike at
      hypocenter is finite and **within 2% of the no-flag run at the
      same tfinal on the same asymmetric mesh** (Zhang 2023 Fig. 6
      shows rupture-front position deviation < 1% between upwind and
      mixed-flux on the same mesh; 2% is a generous gate). If the
      deviation exceeds 5%, the implementation is almost certainly wrong
      (most likely R-1201 wrong matrix, or R-1207 wrong membership).
- [ ] `make test` green.

### Dependencies
- Depends on: Phases 1–5.
- Required by: nothing (terminal phase).

---

## Testing Strategy

| level | what is tested | tools |
|---|---|---|
| unit (R-1101) | central flux primitive correctness | `seas_test_godunov_central_flux` |
| unit (R-1102) | face set topology | `seas_test_mixed_flux_face_set` |
| unit (R-1103) | dispatch bit-identity at flag=none | `seas_test_mixed_flux_dispatch_none` |
| smoke | driver banner + abort | `seas_test_tpv104_smoke` (extended) |
| regression | BP5 / TPV102 / TPV104 default-flag-off path | `make test` |
| convergence (manual) | Zhang Figure 3 reproduction at our mesh | local run + plot |
| dip-direction (target) | slip_dip on `tpv104_repro.msh` at np=1, dx=500 m | compare against P=1 baseline (-1.0078e-02) and P=2 result (-5.80e-05) |

The dip-direction target is the actual goal of this work. Expected outcome:
mixed-flux at P=1 produces slip_dip in the same regime as P=2 (~10⁻⁴ to
10⁻⁵ m at hypocenter), consistent with Zhang's Figure 6/7 result that
mixed-flux on asymmetric mesh ≈ upwind on symmetric mesh.

## Risk Assessment

| risk | likelihood | mitigation |
|---|---|---|
| **R1: HFOs reintroduced when central is used** | LOW for `Adjacent` (Zhang Figure 8 shows minimal HFO impact when only ~5–10% of faces are central); MEDIUM for `AllContinuous` | Default mode = `None`; users must opt in. Document the trade-off in the banner. |
| **R2: Numerical instability under aggressive central-flux usage** | LOW | Zhang reports stable runs at O=4 with mixed-flux 2 on 11k+ asymmetric tets. Our O ≤ 4 + smaller meshes. |
| **R3: Loss of conservation at the fault face** | NIL | Mixed-flux only changes NON-fault faces. Fault Riemann (current `EvaluateADER`) unchanged → fault conservation is unchanged. |
| **R4: Precomputed-flux path breakage** | MEDIUM | Mutual-exclusion abort in driver. Existing tests use precomputed-flux flag = OFF, so default test path is unaffected. |
| **R5: MPI shared-face inconsistency** | MEDIUM | Same as round-7 R-602/R-603 risk: shared-face dispatch must apply mixed-flux uniformly across rank seams. Phase 4 handles both `ComputeSharedFaceFluxRHS` and `ComputeADERSharedFaceFluxRHS` with the same dispatch logic. R-1003-style abort guard if `nprocs > 1 AND use_substep AND mixed_flux != none` may be needed if the testing reveals inconsistency. |
| **R6: Performance** | LOW | One `std::set::count(f)` per face per step. For TPV104 production (~3M faces, ~7000 macro steps): ~2 × 10¹⁰ set lookups, negligible vs the ~2 × 10¹³ flux op cost. |

## Out of scope
- Drucker-Prager plasticity (Zhang Sec 4 application).
- Thermal pressurization (Zhang Sec 4 application).
- Branch / non-planar faults (Zhang's Wenchuan application).
- Adaptive mesh refinement.
- The Wenchuan-scale earthquake demonstration.

These are downstream Zhang-2023 applications that depend on mixed-flux
working at the base level. This plan delivers ONLY the base mixed-flux
machinery; the applications would each be separate plans.

## Sequencing summary

```
Phase 1 (Central primitive) ──┐
                              ├── Phase 4 (dispatch)
Phase 2 (CLI + config) ───────┤
                              ├── Phase 5 (unit tests)
Phase 3 (face set) ───────────┘                    │
                                                   │
                                    Phase 6 (driver smoke + extended test)
```

Phases 1–3 land independently behind the default `None` flag (zero
behavior change). Phase 4 wires the dispatch but is invisible until the
flag is non-default. Phase 5 + 6 verify and expose the feature.

Estimated total: ~600 LOC + 3 unit tests + Makefile edits + driver
update. Comparable scope to round-7 substep-iteration adoption.
