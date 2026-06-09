# Code Review — Cross-rank material exchange, Phase 1 (adversarial)

**Date:** 2026-06-07
**Reviewer:** code-reviewer agent (adversarial; assume >=3 bugs)
**Plan:** `document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md` (v4, Phase 1)
**Prior review:** `REVIEW_plan_cross_rank_round4_2026-06-06.md` (R-401..R-405)

## Scope (files reviewed)
- `dynamic/bimaterial_wave_operator.hpp` (modified) — accessor decls, counter, getter, relaxed ctor doc.
- `dynamic/bimaterial_wave_operator.inl` (modified) — ctor, `SetupMaterialFaceNbrExchange_`,
  `MaterialAtLocal_`, `MaterialAtNbr_`, `BuildGodunovFluxPool_` routing.
- `dynamic/heterogeneous_material.hpp` (read) — `MaterialField` modes, `EvalAt`, `MakeGridFunction`.
- `tests/parallel/test_bimaterial_seam_material_np2.cpp` (new) — the np=2 gate.
- `Makefile` (modified) — build/target/run rule.
- Cross-checked MFEM: `fem/pgridfunc.cpp::GetValue`, `fem/gridfunc.hpp` GetValue contract,
  `mesh/pmesh.cpp::GetSharedFaceTransformations*`, `mesh/mesh.cpp::GetElementTransformation`.

## Verification independently reproduced
- New gate: **48/48 pass** at np=2 (re-ran `mpirun -np 2 ./seas_test_bimaterial_seam_material_np2`).
- Byte-exact regressions (force-rebuilt the .o to defeat the no-header-dep stale-.o hazard):
  - `seas_test_bimaterial_wave_operator_parity` 9/9 (Coefficient Mult/ADER/MaxDt vs scalar to 1e-9).
  - `seas_test_bimaterial_fault_perside_material` 11/11.
  - `seas_test_bimaterial_central_flux` 6/6.
- The DOCUMENTED-DEVIATION byte-exactness claim (flux-pool routing + GridFunction ctor acceptance
  pulled forward) is **confirmed**: the local Constant/Coefficient path is `MaterialField::EvalAt`
  at the same element transform and centroid IP, with NO added `SetIntPoint` — identical to the
  pre-change `material.EvalAt(e, *T, ip, ...)`. Parity 9/9 corroborates no numeric change.

---

## Findings

### P1-001 — MODERATE (latent CRITICAL for Phase 2/4) — `MaterialAtLocal_` aliases the mesh's shared `Transformation` scratch, which IS `ftr->Elem1`

**Category:** UB / aliasing / API-composability

**Description.**
`MaterialAtLocal_` (Constant/Coefficient branch) obtains the element transform via
```cpp
ElementTransformation *T = const_cast<MeshType &>(mesh_).GetElementTransformation(elem);
```
`Mesh::GetElementTransformation(int)` (mesh.cpp:417) returns `&Transformation` — the mesh's single
shared scratch `IsoparametricTransformation`. `ParMesh::GetSharedFaceTransformations(sf)`
(pmesh.cpp:2936, 2982-2983) sets `FaceElemTr.Elem1 = &Transformation` — **the same object**. So any
caller that holds a live `ftr` from `GetSharedFaceTransformations` and then calls `MaterialAtLocal_`
(on ANY element, including a different one) silently mutates `ftr->Elem1`: its `SetIntPoint` state is
overwritten, and if a different element index is passed, `ftr->Elem1` now describes the WRONG element.

This is exactly the hazard the prompt flagged ("aliasing of the mesh's internal
`GetElementTransformation` scratch with `ftr->Elem1/Elem2`"). The test author is clearly aware: in
`RunModeCase` the PEER read is done first, the LOCAL read last (comment inl-test lines 224-225), and
the purity loop **re-fetches `ftr2` after every local read** (lines 248-250) precisely because
`MaterialAtLocal_` clobbers the scratch. That work-around is invisible from the accessor's signature.

**Why it is only MODERATE for Phase 1, but a latent CRITICAL.**
Phase 1's only consumer is `BuildGodunovFluxPool_`, which does NOT hold an `ftr` while looping — so
Phase 1 is currently safe (parity 9/9 confirms). BUT the accessors are introduced specifically for
the Phase 2 shared-neighbour loop and the Phase 4 shared-fault loop, both of which DO hold `ftr`
across `MaterialAtLocal_`/`MaterialAtNbr_` calls in the same iteration. The plan's own Phase-4 sketch
(plan lines 199-201) calls `MaterialAtLocal_(ftr->Elem1No, ip1, ...)` and then
`MaterialAtNbr_(ftr, ip2, ...)` in one iteration. With the current `MaterialAtLocal_`, the first call
reconfigures `Transformation` (= `ftr->Elem1`) to `Elem1No` at `ip1` — geometrically the same element
so the per-side read is OK, but any *subsequent* re-use of `ftr->Elem1` (normal recompute, upwind
operand, `offset_ip(*ftr->Elem1, ...)`) in that loop reads a scratch left at `ip1`, not the face QP.
The bug will surface as a hard-to-localize wrong result the moment a Phase-2/4 loop interleaves a
`MaterialAtLocal_` between fetching `ftr` and re-using `ftr->Elem1`.

**Trigger.** Phase 2/4: any loop that fetches `ftr = GetSharedFaceTransformations(sf)`, calls
`MaterialAtLocal_(...)`, then re-uses `ftr->Elem1`. Also any future caller that holds `ftr` and calls
`MaterialAtLocal_` on a *different* element (then re-uses `ftr->Elem1`).

**Suggested fix (make the accessor scratch-free — pass the transform in, or use the T-overload):**
```diff
-   void MaterialAtLocal_(int elem, const mfem::IntegrationPoint &ip,
-                         real_t &lam, real_t &mu, real_t &rho) const;
+   // Caller supplies the ElementTransformation (its OWN, not the mesh scratch),
+   // so the accessor never touches Mesh::Transformation and is composable with a
+   // held FaceElementTransformations.  `elem` retained only for GridFunction
+   // GetValue(elem,ip) (which is local-scratch-free for L2 VALUE FEs).
+   void MaterialAtLocal_(int elem, mfem::ElementTransformation &T,
+                         const mfem::IntegrationPoint &ip,
+                         real_t &lam, real_t &mu, real_t &rho) const;
```
and in the body drop the `GetElementTransformation` call (use the passed `T`). The flux-pool loop
already creates its own `T` per element; the Phase-4 loop passes `*ftr->Elem1`. Minimal alternative:
keep the signature but document loudly that `MaterialAtLocal_` invalidates any held `ftr->Elem1`, and
require callers to use it AFTER all `ftr->Elem1` uses in the iteration — but a self-documenting
signature is safer.

**Test case:** `test_P1001_accessor_does_not_clobber_ftr` (np=2): fetch `ftr`; record
`ftr->Elem1->Transform(centroid)` -> `x0`; call `MaterialAtLocal_(other_local_elem, ...)`; assert
`ftr->Elem1->Transform(centroid)` still == `x0` (currently FAILS — scratch overwritten).

---

### P1-002 — MODERATE — Three material GFs sharing ONE space is assumed, not asserted; `MaterialAtNbr_` indexes `mu_gf`/`rho_gf` with `lambda_gf`'s vdofs

**Category:** ASSUMPTION / silently-wrong-material

**Description.**
`MaterialAtNbr_` (GridFunction branch) computes `vd` (face-nbr vdofs) and `fe`/`ndof` ONCE from
`lambda_gf->ParFESpace()`, then the `interp` lambda reuses that same `vd`/`shape`/`ndof` to read
`mu_gf.FaceNbrData()` and `rho_gf.FaceNbrData()` (inl:165-196). This is correct ONLY IF all three GFs
share the SAME `ParFiniteElementSpace` (identical face-nbr dof layout). The `MaterialField` contract
docstring asserts this ("All three GFs share the SAME ParFiniteElementSpace", hpp:81-82), but
**nothing enforces it**: `MakeGridFunction` (heterogeneous_material.hpp:110-124) only checks non-null,
and `SetupMaterialFaceNbrExchange_` only checks each `GetVDim()==1`. If a future CVM build hands in
`mu_gf` on a different space (different order, or a separately-partitioned face-nbr layout), the peer
read returns `mu`/`rho` indexed by `lambda`'s dof map — silently wrong material on every seam.

`MaterialAtLocal_`'s GridFunction path is robust to this (each `GetValue(elem,ip)` uses its own GF's
space), so the asymmetry is `MaterialAtNbr_`-only and easy to miss.

**Trigger.** A GridFunction material whose 3 GFs do not share one ParFESpace (the stated CVM goal).

**Suggested fix:** assert shared space in `MakeGridFunction` and/or `SetupMaterialFaceNbrExchange_`:
```cpp
MFEM_VERIFY(material_->lambda_gf->ParFESpace() == material_->mu_gf->ParFESpace()
            && material_->lambda_gf->ParFESpace() == material_->rho_gf->ParFESpace(),
            "material GridFunctions must share ONE ParFiniteElementSpace "
            "(MaterialAtNbr_ reuses lambda's face-nbr vdof layout for mu/rho).");
```

**Test case:** `test_P1002_gridfunction_disjoint_spaces_aborts` — build with `mu_gf` on a 2nd
(distinct) ParFESpace of the same order; assert the ctor aborts (currently passes silently and would
mis-read at a seam where the two spaces' face-nbr orderings differ).

---

### P1-003 — MODERATE — R-402 only half-enforced: `vdim==1` checked, but space-order == state-order (and L2 map type) is NOT

**Category:** ASSUMPTION / contract gap (plan edge-case unimplemented)

**Description.**
The plan's Phase-1 Edge Cases (plan line 129) require: "GridFunction source space != state space
order: abort with a clear message." The round-4 R-402 fix asked for both `vdim==1` AND "FE order ==
state order". The implementation enforces only `GetVDim()==1` (inl:96-100); there is NO check that the
material GF order equals the state/wave order, and no check that the FE is L2/VALUE-mapped.

Consequences if violated:
- Order mismatch: `BuildGodunovFluxPool_` samples the material at the element CENTROID via the GF; a
  higher/lower-order material GF still returns a self-consistent centroid value, so it would NOT
  abort, it would just silently use a material sampled at a different polynomial resolution than the
  state assumes — divergent from the intended "scalar L2(order p)" contract and from the Coefficient
  path's analytic value.
- Map-type: `MaterialAtNbr_`'s manual `interp` uses `fe->CalcShape(ip_peer, shape)` and a naive vdof
  sign decode, omitting `DofTransformation::InvTransformPrimal` (which MFEM's own `GetValue` applies,
  pgridfunc.cpp:299). For scalar L2/GLL this is identity (correct), but for an H1/ND/RT scalar GF the
  manual decode would diverge from MFEM. The vdim==1 assert does not exclude these.

**Trigger.** A material GF whose order != state order, or a non-L2 scalar space (future CVM).

**Suggested fix:** in `SetupMaterialFaceNbrExchange_`, additionally assert
```cpp
MFEM_VERIFY(material_->lambda_gf->ParFESpace()->GetMaxElementOrder() == /*state order*/ order_,
            "material GridFunction order must equal the state order.");
// and that the FE collection is L2 / map type VALUE, e.g. via FEColl()->GetMapType or a name check.
```
(`order` is available to the ctor; thread it to a member or pass it in.) Alternatively, replace the
manual `interp` with `gf.GetValue(ftr->Elem2No, ip_peer)` (MFEM's face-nbr path) so the
DofTransformation is handled by the library and the only remaining contract is order==state-order.

**Test case:** `test_P1003_gridfunction_order_mismatch_aborts` — material GF at order 2 with a
state/wave order 1 operator; assert ctor aborts (currently passes silently).

---

### P1-004 — LOW — `MaterialAtNbr_` re-implements MFEM's existing face-nbr `GetValue`; the hand-rolled path is redundant risk surface

**Category:** QUALITY / maintainability

**Description.**
`ParGridFunction::GetValue(int i, ip)` (pgridfunc.cpp:270-312) ALREADY handles `i >= GetNE()` as a
face-neighbour read: it does the identical `GetFaceNbrElementVDofs` + `GetFaceNbrFE` +
`face_nbr_data` shape interpolation, plus `DofTransformation::InvTransformPrimal` (which the manual
version omits — see P1-003). The manual `interp` lambda in `MaterialAtNbr_` (inl:178-193)
reimplements this by hand. For scalar L2 the two are equivalent (verified: gate peer reads match to
rel 0.0), but the hand-rolled version is extra surface for the P1-002/P1-003 bugs and silently
diverges for non-L2 spaces.

**Suggested fix:** replace the manual block with
```cpp
lam = material_->lambda_gf->GetValue(ftr->Elem2No, ip_peer);
mu  = material_->mu_gf->GetValue(ftr->Elem2No, ip_peer);
rho = material_->rho_gf->GetValue(ftr->Elem2No, ip_peer);
```
which routes each GF through its OWN space (also fixes P1-002 for free) and applies the proper
DofTransformation. Keep the `nbr_idx>=0` guard for the clear error message.

**Test case:** covered by the existing peer-read gate (would still pass); add a guard test that
`MaterialAtNbr_` equals `gf.GetValue(Elem2No, ip_peer)` to 0 ulp for the GridFunction mode.

---

### P1-005 — LOW — Gate's "PEER mirror" + cref reuse passes correctly but only on a seam-symmetric mesh; the R-301 purity claim is overstated

**Category:** TEST_QUALITY (not a code bug)

**Description.**
The gate's independent peer expectation uses `peer_cx = 2*seam_x - local_cx` (test line 207) and feeds
the **Elem1** reference centroid `cref` as `ip_peer` to `MaterialAtNbr_` (test line 216). Both are
valid ONLY because the mesh is a row of unit hexes symmetric about the seam, so (a) the mirror equals
the peer physical centroid and (b) the Elem1 reference centroid equals the Elem2 reference centroid
(both `(0.5,0.5,0.5)`). On a non-axis-aligned or non-unit seam mesh, `peer_cx` would not equal the
peer centroid and `cref` would not be the peer's correct reference point — the assertions would be
checking the wrong thing. The check is genuinely independent (not tautological) on THIS mesh, but its
correctness is geometry-dependent and undocumented.

Separately, the file/inline comments claim the halfspace-jump case "Covers R-301-class pollution"
(test line 16). R-301 is node-on-jump pollution from PROJECTING a jump material onto GLL nodes; the
Coefficient mode evaluates analytically (no nodes) and the GridFunction mode here is per-element
constant (jump is across the seam, never within an element). So the gate confirms pure per-side reads
at near-seam IPs but does not exercise the node-on-jump projection scenario R-301 describes — it
cannot, because v4 never projects. The claim is harmless but slightly overstated.

**Suggested fix:** add a comment that the mirror/`cref` reuse relies on the seam-symmetric unit-hex
mesh; soften the R-301 comment to "no-projection => no node-on-jump pollution by construction."

**Test case:** none required (test passes); optional: pass `ftr->GetElement2IntPoint()` as `ip_peer`
to make the peer reference frame explicit and mesh-independent for future non-symmetric gates.

---

### P1-006 — LOW — `BuildGodunovFluxPool_` keeps an unused `const MaterialField& material` parameter

**Category:** DEAD_CODE / clarity

**Description.**
After routing through `MaterialAtLocal_` (which reads `material_`), the `material` parameter is unused
(`(void)material;`, inl:228). It is also redundant with the member `material_` (set before the call).
Harmless, but a dead parameter invites a future caller to pass a DIFFERENT `MaterialField` than
`material_`, which the body would then silently ignore.

**Suggested fix:** drop the parameter from `BuildGodunovFluxPool_()` (it already uses `material_`), or
add `MFEM_ASSERT(&material == material_, ...)`.

---

## Items from round-4 — verification

| ID | Round-4 ask | Status in code |
|----|-------------|----------------|
| R-401 | peer `offset_ip` on `ftr->Elem2` Jacobian | N/A to Phase 1 (Phase-4 path not yet implemented); not exercised. |
| R-402 | assert material GF scalar (vdim==1) AND order==state | **HALF done** — vdim==1 asserted (inl:96-100); order/L2 NOT (see P1-003). |
| R-403 | `const_cast` for non-const mesh accessor | DONE (inl:135-136) — but the const_cast is itself the P1-001 aliasing hazard. |
| R-404 | assert `material_ != nullptr` in accessors | DONE (inl:120, 149). |
| R-405 | 3 separate GF exchanges (micro-opt) | As-is (inl:103-105); acceptable, collective-symmetric on all ranks. |

## Collective / ctor-order / byte-exactness checks (all PASS)
- **Collective symmetry (R-5/R-209):** `SetupMaterialFaceNbrExchange_` runs `pmesh.ExchangeFaceNbrData()`
  and (mode==GridFunction, identical on all ranks) the 3 GF exchanges unconditionally — no
  rank-conditional collective. A rank with no shared faces still calls them. OK.
- **Ctor order (R-6):** `SetupMaterialFaceNbrExchange_` (54) precedes `BuildGodunovFluxPool_` (56,
  local-only), `ExchangeBiMaterialNeighbours_` (57, local stub — reads no peer material),
  `BuildPerFaceBimaterialFluxMatrices_` (58). No peer material is read before the exchange. OK.
- **`material_` set before use:** `material_ = &material` (49) precedes all reads. `(void)material;`
  is correct. OK.
- **`nbr_idx = Elem2No - ne_`:** matches MFEM (`Elem2No = NumOfElements + nbr_idx`, pmesh.cpp:2995;
  `ParGridFunction::GetValue` uses `i - GetNE()`, pgridfunc.cpp:275). Correct.
- **Byte-exact local path:** legacy `material.EvalAt(e,*T,ip)` == new `MaterialAtLocal_` Const/Coeff
  branch (same `GetElementTransformation`, same IP, no added `SetIntPoint`). Parity 9/9 confirms.
- **vdof sign decode in `interp`:** harmless no-op for scalar L2/GLL (non-negative vdofs); not wrong
  for Phase-1 use, but see P1-003 for the H1/ND/RT divergence (omits DofTransformation).

---

## Summary

| ID | Severity | Category | One-line |
|----|----------|----------|----------|
| P1-001 | MODERATE (latent CRITICAL Ph2/4) | aliasing/UB | `MaterialAtLocal_` clobbers `Transformation` == `ftr->Elem1` scratch |
| P1-002 | MODERATE | silently-wrong-material | 3-GF-shared-space assumed, not asserted; `MaterialAtNbr_` indexes mu/rho with lambda's vdofs |
| P1-003 | MODERATE | contract gap | R-402 order==state (and L2 map type) not enforced; plan edge-case unimplemented |
| P1-004 | LOW | quality | manual face-nbr interp re-implements `GetValue`; redundant risk surface |
| P1-005 | LOW | test quality | peer mirror/cref reuse + R-301 claim mesh-specific/overstated |
| P1-006 | LOW | dead code | unused `material` parameter in `BuildGodunovFluxPool_` |

- **Critical issues: 0** (no defect breaks Phase-1 as-shipped: the only Phase-1 consumer is the
  scratch-safe flux-pool loop; 48/48 gate + 9/9+11/11+6/6 regressions reproduced).
- Moderate issues: 3 (P1-001 is the one to fix before Phase 2/4 lands; P1-002 and P1-003 close
  unenforced GridFunction preconditions that the v4 goal — GridFunction/CVM first-class — will hit).
- Low issues: 3.

**Verdict: PASS WITH FIXES.**
Phase 1 as delivered is correct and byte-exact for Constant/Coefficient and the documented
deviation (flux-pool routing + GridFunction ctor acceptance) is verified safe. The DEVIATION does not
change any numeric result. However, the new accessors carry a real aliasing hazard (P1-001) that is
inert in Phase 1 only because its sole consumer does not hold an `ftr`; it must be fixed (or its
contract hard-documented and the Phase-2/4 call sites ordered around it) BEFORE the Phase-2 neighbour
loop and Phase-4 fault loop consume the accessors, or those phases will read corrupted `ftr->Elem1`
state. P1-002/P1-003 should land alongside Phase 1 since they protect the GridFunction mode that
Phase 1 just made constructible.
