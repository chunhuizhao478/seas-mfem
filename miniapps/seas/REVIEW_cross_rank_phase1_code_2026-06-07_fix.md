# Fix Report — Cross-rank Phase 1 (against REVIEW_cross_rank_phase1_code_2026-06-07.md)

**Date:** 2026-06-07. Review verdict was PASS WITH FIXES, **0 critical**, 3 moderate, 3 low.
All 6 findings addressed in a single pass. Per the loop policy (no critical ⇒ proceed after fix),
Phase 1 is complete; no re-review mandated, but the gate + 4 regressions were re-run green below.

## Changes made
- **P1-001 (MODERATE, latent CRITICAL):** `MaterialAtLocal_` now takes a caller-supplied
  `ElementTransformation &T` (`MaterialAtLocal_(int elem, ElementTransformation &T, const
  IntegrationPoint &ip, ...)`). The Constant/Coefficient branch calls `material_->EvalAt(elem, T,
  ip, ...)` — NO internal `Mesh::GetElementTransformation`, so the accessor never touches the mesh
  shared scratch that `ParMesh::GetSharedFaceTransformations` aliases as `ftr->Elem1`. This matches
  the existing fault-loop pattern and makes the accessor composable with a held `ftr` (Phase 2/4).
  `BuildGodunovFluxPool_` passes its own per-element `T`. Added the np=2 guard
  `MaterialAtLocal_ leaves a held ftr->Elem1 intact` (records ftr->Elem1 centroid, calls the
  accessor on a DIFFERENT local element via its own transform, asserts ftr->Elem1 unchanged).
- **P1-002 (MODERATE) + P1-004 (LOW):** `MaterialAtNbr_` GridFunction peer read now uses the MFEM
  library path `gf->GetValue(ftr->Elem2No, ip_peer)` per GF (it treats `i >= GetNE()` as a
  face-neighbour read). Each GF reads through its OWN `ParFiniteElementSpace`, so the three GFs need
  NOT share one space (the prior hand-rolled interp reused lambda's vdof map for mu/rho — the bug).
  Removed the manual `GetFaceNbrElementVDofs`/`FaceNbrData`/shape interpolation.
- **P1-003 (MODERATE):** `SetupMaterialFaceNbrExchange_` now asserts, per GF, `vdim==1` AND
  `GetMaxElementOrder() == state order` (`this->GetFESpace().GetMaxElementOrder()`) — the plan's
  "scalar L2(order p)" edge case. Map-type/DofTransformation is now handled by the library `GetValue`
  (so no separate L2 assert is needed for correctness; the order assert pins the contract).
- **P1-005 (LOW):** test comments document that the `peer_cx = 2*seam - local_cx` mirror and the
  `cref`-as-peer-IP reuse are valid only on the seam-symmetric unit-hex mesh; the "R-301" wording
  softened to "no-projection ⇒ no node-on-jump pollution by construction". Local reads use an
  independent `IsoparametricTransformation` (filled overload — scratch-free), so `ftr` stays valid
  without re-fetching.
- **P1-006 (LOW):** dropped the unused `const MaterialField&` parameter from `BuildGodunovFluxPool_()`
  (it reads `material_`, set in the ctor first).

## Verification (re-run after fixes)
- Gate `seas_test_bimaterial_seam_material_np2` (np=2): **56/56** (both np=2 gates + the new P1-001
  guard, all 3 modes, both ranks).
- Byte-exact regressions: parity **9/9**, central-flux **6/6**, perside-material **11/11**,
  mixed-flux-shared (np=2) **14/14**.

## Ready for next phase: YES (0 critical; all findings resolved; tests green).
