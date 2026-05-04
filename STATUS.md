# STATUS — fix/dynamic-rupture-mixed-flux-p2

**Stream:** [4] Debug P>=2 mixed-flux dispatch failure in dynamic-rupture path.
**Reference docs:**
- `miniapps/seas/debug_document/mixed_flux_p2_debug.md` (running debug log; primary).
- `miniapps/seas/document/system_dev/drdg3d_mixed_flux_comparison_2026-04-28.md` (comparison spec).
**Branch:** `fix/dynamic-rupture-mixed-flux-p2` (rooted at `1649cc3`).
**Goal:** identify and fix the P>=2 mixed-flux numerical defect; ship a regression test that catches it.

---

## Current state

**Phase:** Step 1 (understand symptom) — about to begin.
**Last commit on this branch:** `1649cc3` (inherited from `feature/elasticity-inertia`).
**Local build status:** unverified.
**Latest Frontera submission:** none on this branch.

---

## Reproducer (smallest known case)

`miniapps/seas/jobs/tpv104/tpv104_mixed_flux_adjacent_200m_p2_O2_normal.sbatch`
(committed at `1649cc3`; runs P=2 adjacent mixed-flux on 200m TPV104 mesh).

For dev-queue iteration use the same sbatch but adjust `-p development` and shorten `tfinal` if needed (reference: existing `jobs/tpv102/tpv102_200m_p1_1.5s_50rank_dev.sbatch`).

---

## Iteration log

| iter | date | commit | sbatch_id | result | notes / next action |
|------|------|--------|-----------|--------|---------------------|
| 0    | 2026-05-04 | 1649cc3 | n/a | INIT | branch + worktree created; STATUS.md skeleton placed |

---

## Debug skill workflow (per iteration)

Following `code-debug` skill structure. For each cluster iteration:

1. **Symptom**: exact error / wrong-output / NaN observed. Cite log line.
2. **Reproducer**: sbatch + parameters that triggered it.
3. **Hypothesis**: ranked list of possible causes (most likely first).
4. **Investigation**: which hypothesis; what diagnostic was added; what the diagnostic showed; CONFIRMED / ELIMINATED / INCONCLUSIVE.
5. **Root cause** (only after a hypothesis is CONFIRMED): exact file:line + WHY the code is wrong.
6. **Fix** (minimal change): diff applied.
7. **Verify**: re-run sbatch; confirm symptom gone; full test suite green; remove diagnostics.

Each iteration appends a new section to `debug_document/mixed_flux_p2_debug.md`.

---

## No-touch boundary

- `friction/`, `config/bp5_params.hpp` — keep physics constants frozen.
- DG paths unrelated to mixed-flux dispatch (Godunov, central, etc.).
- Quasi-dynamic SEAS solver (`solver/seas_operator.hpp`).

---

## Cross-stream priority

This stream **unblocks stream [3]** (`refactor/dynamic-rupture-unified`).
Land [4] cleanly before starting [3].

---

## Next concrete action

1. `cd ~/projects/seas-mfem-mixedflux/miniapps/seas`
2. `make seas_tpv104_driver` — confirm clean build on this branch.
3. Read existing `debug_document/mixed_flux_p2_debug.md` end-to-end; pull the most-recent Symptom + Hypotheses into Step 1 below.
4. Submit `jobs/tpv104/tpv104_mixed_flux_adjacent_200m_p2_O2_normal.sbatch` from Frontera if a fresh reproducer log is needed.
5. Begin debug skill Step 3 (form hypotheses) → Step 4 (instrument & investigate). Keep diagnostics SEPARATE from any fix.
