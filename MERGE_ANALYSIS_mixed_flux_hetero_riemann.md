# Merge Analysis: `feature/mixed-flux-hetero-riemann` → `system/spatial_dyn_driver`

**Date:** 2026-06-05
**Verdict:** ✅ **Clean merge, no textual conflicts, semantically sound.** Rebuild + full test suite on the merged tree is the only required validation. Do **not** treat the clean `merge-tree` as a license to skip building — two features touch the *same* scalar `WaveOperator` infrastructure.

---

## 1. Branch topology

| Ref | Commit | Position |
|---|---|---|
| merge-base | `04674fc` | common ancestor |
| `feature/mixed-flux-hetero-riemann` (HEAD) | `e8a0db7` | 1 commit past base |
| **local** `system/spatial_dyn_driver` | `3507b78` | **30+ commits past base** (incl. `Merge feat/fault-overint-resample`) |
| `origin/system/spatial_dyn_driver` | `f127ef9` | ancestor of `3507b78` — **origin is BEHIND local** |

⚠️ **Two different targets exist and they are out of sync.** `origin/` (`f127ef9`) is 30 commits behind the local branch (`3507b78`). The merge that matters is into the **local** `3507b78`. The local branch is checked out in the **main repo worktree** (`/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver`), not here.

---

## 2. Scenario A — merge into `origin/system/spatial_dyn_driver` (`f127ef9`)

- `git merge-tree` → **exit 0 (clean)**.
- **Zero overlapping changed files.** The `f127ef9` commit only refreshes **gold data** (`tpv102/gold`, `tpv104/gold`, `tpv205/gold` — 283 files: `.dat`/`.png`/`.sbatch`), touching **no source, no Makefile, no tpv31**.
- The feature touches source + tests + tpv31; the two are fully disjoint.
- Note: `f127ef9`'s "RK45 mixed-flux" gold used the **scalar** mixed-flux path (at base the bi-material operator *aborted* on any mixed mode). The feature keeps that scalar path **byte-exact** (see §4), so this gold stays regenerable.

## 3. Scenario B — merge into **local** `system/spatial_dyn_driver` (`3507b78`) ← the real target

- `git merge-tree --write-tree 3507b78 HEAD` → **exit 0 (clean)**, merged tree `fb50b23`.
- **4 source files modified by BOTH branches** (auto-merged, non-overlapping hunks):

| File | base | feature Δ | target Δ | merged | exact? |
|---|---|---|---|---|---|
| `dynamic/wave_operator.inl` | 6829 | −11 | +77 | 6895 | ✅ 6829−11+77 |
| `drivers/spatial_dyn_driver.cpp` | 3242 | +35 | +165 | 3442 | ✅ 3242+35+165 |
| `dynamic/wave_operator.hpp` | 1125 | +17 | +39 | 1181 | ✅ 1125+17+39 |
| `Makefile` | 5333 | +165 | +88 | 5586 | ✅ 5333+165+88 |

Every merged blob = base + both deltas **exactly** → clean union, no lost/duplicated lines.

### Why the shared-file merge is genuinely clean (not a dangerous splice)
- **`wave_operator.inl`** — feature edits only `SetMixedFluxMode` (~1620, comment) and `ComputeMaxDt`/new `MixedFluxCflFactor_` (~5834). Target edits the ctor, `EvaluateBulkAtFaultQPsCanonical`, the flux-RHS functions, and adds a ~95-line over-integration region (`RebuildFaultQuadrature_`/`SetFaultOverint`). **Disjoint regions.** Merged blob has exactly **1** `ComputeMaxDt`, **1** `SetMixedFluxMode`, **1** `MixedFluxCflFactor_` def.
- **`spatial_dyn_driver.cpp`** — feature hunks at base 803/1134/1244/1896; target hunks at 56–561/1273–1444/2615. **Disjoint.** Merged blob has exactly **1** `main()`.
- **`wave_operator.hpp`** — feature adds `MixedFluxCflFactor_` decl (1, no dup); target adds `SetFaultOverint`/`GetFaultOverint`/`FaultFaceQuadDegree`/`RebuildFaultQuadrature_`/`fault_overint_k_`. No name collisions.
- **`Makefile`** — target's edits (lines 397/703/1844/3588/4357) **never touch the `test:` aggregate rule**. Merged `test:` = target's rule **+ feature's 4 tests** (`test-bimaterial-central-flux/-mixed-flux-{dispatch,cfl,homog-equivalence}`), **zero prerequisites lost**, single `test:` rule, no duplicate target heads.

### Feature-owned files with NO target overlap (merge in verbatim)
`dynamic/bimaterial_wave_operator.{hpp,inl}`, `dynamic/godunov_flux_bimaterial.{cpp,hpp}`, `spatial/code/spatial_friction.{cpp,hpp}` (incl. `MatrixMixedFluxUnderAder`, `seam_continuous`), all new tests, `tpv31/*`, both `.sbatch` jobs. Target did not touch any of these.

---

## 4. Semantic / behavioral assessment (the part `merge-tree` cannot check)

1. **Scalar `ComputeMaxDt` is byte-exact.** Feature refactors the CFL factor into `MixedFluxCflFactor_()` but returns identical values (None→1.0, Adjacent→`rk_aware?0.6:0.9`, AllContinuous→`rk_aware?0.7:0.4`) and the same formula. The *only* new behavior is a fail-loud guard `MFEM_VERIFY(mixed_flux_mode_==None || cfl_rk_aware_)`. ⇒ existing TPV/BP5/scalar-SAFS runs and `f127ef9`'s RK45 gold are unaffected.
2. **Bi-material `mixed_flux=none` is byte-exact.** All 17 new central-flux branches are gated by `mf_on_`/`central_flux_face_set_`; `MixedFluxCflFactor_()` returns 1.0 for None; the new guard passes for None. ⇒ existing heterogeneous runs unchanged.
3. **★ Cross-feature interaction — over-integration × mixed-flux — is already guarded and composes correctly.** The target's `SetFaultOverint(k)` contains:
   ```cpp
   MFEM_VERIFY(k == 0 || (mixed_flux_mode_ == MixedFluxMode::None &&
                          !use_precomputed_face_fluxes_),
               "SetFaultOverint: fault over-integration (k>0) is not compatible "
               "with the mixed-flux / precomputed-face-flux paths ...");
   ```
   At base this never had to consider the bi-material operator (it couldn't do mixed-flux). Post-merge it can — and the guard **still fires correctly** because both operators share the base `mixed_flux_mode_` member. The merged driver call order makes this sound:
   ```
   SetSeamContinuous (1167) → SetMixedFluxMode (1276) → SetFaultOverint (1332) → SetFaultDOFData (1859)
   ```
   `SetFaultOverint` runs **after** `SetMixedFluxMode`, so its guard reads the real mode. The driver comment (line 1320–1321) was already authored: *"Placed AFTER SetMixedFluxMode so SetFaultOverint's guard sees the real mixed-flux mode."* ⇒ a config requesting `interior_flux=matrix + mixed_flux=adjacent + fault_overint>0` aborts cleanly; the two features are mutually excluded exactly where required.

### Residual risk (build/test only, not a conflict)
Both features edit shared scalar `WaveOperator` infrastructure (target rewrote flux-RHS bodies + fault-quadrature; the bi-material path inherits/calls those). No signature changes were made to the flux-RHS overrides, but the *combination* (bi-material mixed-flux present alongside over-integration code) has **never been compiled or tested together**. This is the one thing the analysis cannot certify on paper.

---

## 5. Recommended merge procedure (when you approve execution)

1. The target `3507b78` is checked out in the **main worktree**, and `origin` is 30 commits behind it — run the merge **from the main worktree**, or fast-forward `origin` first. Don't merge into `origin/f127ef9` (you'd lose the 30 local commits' integration surface).
2. `git merge feature/mixed-flux-hetero-riemann` → expected clean (matches `merge-tree` tree `fb50b23`).
3. **Rebuild the full `seas` test suite on the merged tree** (`MFEM_*` overrides per the worktree-build memo). This is mandatory — it is the only check for the §4 residual risk.
4. Run: the feature tests (central 6, dispatch 10, cfl 17, homog 18, wave_operator 25, shared 14 @np2) **and** the target's `test-fault-overint`/`test-fault-resample`/`test-fault-resample-apply`/`test-fault-planar-serial` **and** `test-tpv-config-parse` (2 pre-existing VS-border fails are expected/known).
5. Optional but recommended: add one config-parse assertion that `interior_flux=matrix + mixed_flux=adjacent + fault_overint>0` aborts (locks the §4.3 composition).

## 6. Housekeeping
- Feature worktree has uncommitted noise (`reports/.report-config`, `miniapps/seas/reports/`) — **not** in commit `e8a0db7`, will not merge. Safe to ignore.
- Staged Frontera jobs (tpv31 p1 12N/600c, p2 16N/800c) are untouched by the target and require explicit user approval before submission.
