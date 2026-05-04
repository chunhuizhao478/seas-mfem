# STATUS — feature/safs-quasi-dynamic

**Stream:** [1] SAFS mesh build + quasi-dynamic smoke driver.
**Plan:** `miniapps/seas/safs/smoke_test/PLAN_safs_test.md` (PDF rendered alongside).
**Branch:** `feature/safs-quasi-dynamic` (rooted at `1649cc3`).
**Goal:** end-to-end SAFS smoke driver passing acceptance on Frontera dev queue (8 nodes, 400 cores, 2 hr).

---

## Current state

**Phase:** not started
**Last commit on this branch:** `1649cc3` (inherited from `feature/elasticity-inertia`)
**Local build status:** unverified
**Latest Frontera submission:** none

---

## Iteration log

| iter | date | commit | sbatch_id | result | notes / next action |
|------|------|--------|-----------|--------|---------------------|
| 0    | 2026-05-04 | 1649cc3 | n/a | INIT | branch + worktree created; STATUS.md skeleton placed |

---

## Open questions / known risks

- Plan §Risk Assessment R-733: per-fault `ref_normal = (0,-1,0)` is geometrically inconsistent across the 6 SAFS faults. Documented as acceptable for shakedown; revisit if `--check-residual` reports ill-conditioning.
- Plan §R-006 (moderate, applied): dt_init formula must use `max(V_init, Vp)`, not just `Vp`. If a future override sets `V_nuc != V_init`, re-validate this.
- Pipeline mesh `safs_newset_6_cavity.msh` (1.1 M tets) cannot factor on a laptop — local cycle is build + `make test` ONLY.

---

## Next concrete action

1. `cd ~/projects/seas-mfem-safs/miniapps/seas`
2. `make test` — confirm baseline (157 SAFS-mesh tests + existing seas tests all pass on this branch).
3. Implement Phase 1 of `PLAN_safs_test.md`:
   - Create `safs/code/safs_test_params.hpp`
   - Create `safs/code/safs_boundary_config.hpp`
   - Build `seas_safs_test_driver` Makefile target (additive only).
4. Local acceptance: `make seas_safs_test_driver` no warnings + `make test` still green.
5. Commit on `feature/safs-quasi-dynamic`.

After Phase 1 lands locally, move to Phase 2 (BP5Params override).
