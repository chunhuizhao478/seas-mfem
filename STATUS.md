# STATUS — feature/paraview-compaction

**Stream:** [2] ParaView output compaction for quasi-dynamic SEAS runs.
**Plan:** `miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md`.
**Branch:** `feature/paraview-compaction` (rooted at `1649cc3`).
**Goal:** replace 4M-tiny-files-per-run output with binary VTU (Phase 1) → VTKHDF single-file (Phase 2) → optional ZFP lossy (Phase 2d).

---

## Current state

**Phase:** Phase 0 (recovered untracked WIP files into the branch; not yet committed).
**Last commit on this branch:** `1649cc3` (inherited from `feature/elasticity-inertia`).
**Local build status:** unverified — WIP headers + tests are present but not yet compiled.
**Latest Frontera submission:** none.

---

## Files in flight on this branch (untracked, recovered from main worktree)

```
miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md
miniapps/seas/io/fault_vtkhdf_writer.hpp                        ← Phase 2b
miniapps/seas/io/fault_vtu_binary.hpp                            ← Phase 1
miniapps/seas/tests/unit/test_fault_surface_vtkhdf.cpp           ← Phase 2b
miniapps/seas/tests/unit/test_fault_surface_vtkhdf_mpi.cpp       ← Phase 2b
miniapps/seas/tests/unit/test_fault_surface_vtu_binary.cpp       ← Phase 1
miniapps/seas/tests/unit/test_fault_surface_vtu_gather_mpi.cpp   ← Phase 1
miniapps/seas/tests/unit/test_paraview_hdf_smoke.cpp             ← Phase 2b
```

Note: the previous round of edits to `paraview_output.hpp`, drivers, and existing fault-VTU tests were **reverted** on the main checkout; this branch's `1649cc3` does NOT carry them. Re-apply incrementally per the plan, NOT in one batch.

---

## Iteration log

| iter | date | commit | sbatch_id | result | notes / next action |
|------|------|--------|-----------|--------|---------------------|
| 0    | 2026-05-04 | 1649cc3 | n/a | INIT | branch + worktree created; WIP files recovered; STATUS.md skeleton placed |

---

## Phase ordering

| Phase | Deliverable | Acceptance gate |
|---|---|---|
| 1   | Binary VTU per rank with ZLib level 9 | unit tests pass on np=1,4,8; ≥5× size reduction on TPV102 100-cycle reference run |
| 2   | VTKHDF single-file (one .vtkhdf per simulation) | ParaView opens it cleanly; existing TPV102 visualize_results pipeline still works |
| 2d  | ZFP lossy compression (opt-in) | Round-trip error ≤ 1e-6 on traction/slip; ≥10× over Phase 2 |

Land Phase 1 fully (commit + Frontera reference run) before starting Phase 2. Do NOT batch all phases into one commit.

---

## No-touch boundary

- `friction/`, `fault/` — physics; this branch is I/O only.
- `bp5/`, `bp1/`, `bp2/` driver entry points beyond the documented call-sites.
- The mixed-flux dispatch path (that's stream [4]).

---

## Next concrete action

1. `cd ~/projects/seas-mfem-paraview/miniapps/seas`
2. Inspect the recovered untracked files; confirm they compile (`fault_vtu_binary.hpp` may need including in `paraview_output.hpp`).
3. Implement Phase 1 binary VTU writer + its 2 unit tests; verify they pass.
4. Commit Phase 1 only (do not include Phase 2 / 2d files in this commit).
5. Push branch to origin.
6. Submit a TPV102 100-cycle reference Frontera run to confirm output size reduction.

Phase 2 (VTKHDF) starts after Phase 1 lands.
