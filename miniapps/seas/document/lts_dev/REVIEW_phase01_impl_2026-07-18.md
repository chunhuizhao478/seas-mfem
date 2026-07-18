# LTS Phase 0/1 — implementation review & fix record (2026-07-18)

Branch: `safs-v4_0_0-alt-case1-mfem-speed` (local only — not pushed).
Plan: `document/lts_dev/PLAN_clustered_lts_ader_2026-07-18.md` (rev 4).

## What was implemented (strict to the plan)

| Deliverable | Plan ref | Files | Local validation |
|---|---|---|---|
| Rate-2 clustering (integer-loop binning, maxdiff fixpoint, λ-scan, Nc-cap auto-merge, GAP-A1 assert, speedup stats) | A.1 | `dynamic/lts_clustering.{hpp,cpp}` | `test_lts_clustering` **4043/4043** |
| `--lts-report` driver hook (RAW + PRODUCTION histograms, harmonic/arithmetic speedups) + per-element CFL accessors | Phase 0 | `drivers/spatial_dyn_driver.cpp`, `dynamic/wave_operator.{hpp,inl}`, `dynamic/bimaterial_wave_operator.hpp` | driver links; byte gates green |
| Run-side layout (per-cluster elem/face lists, FaceRole tagging, provider/consumer sets + dense slot maps, local metadata + reduction) | A.2 | `dynamic/lts_layout.{hpp,cpp}` | `test_lts_layout` **87/87** |
| Tick-table generator (predict/correct predicates, FINE→COARSE order, truncated dt_step, matched-collective count) | A.3 | `dynamic/lts_stepper.hpp` | covered by `test_lts_layout` |
| Config parse `[numerics].lts*` + validators + guards (lts+mixed_flux, lts+rk) | Phase 1 | `spatial/code/spatial_friction.{hpp,cpp}`, `drivers/spatial_dyn_driver.cpp` | `test_spatial_friction_config` **314/314** |
| Run-path clustering + layout wiring (gated, still GTS) | Phase 1 | `drivers/spatial_dyn_driver.cpp` | full driver links |

**Byte-exact-off contract:** every LTS branch is gated on `lts != "off"`; the
wave-operator change only *stores* an already-computed `h_e` (`h_min_` still its
min). Confirmed byte-neutral: `test-ader-tpv102-smoke` 4/4, `test-wave-operator`
25/25, `test-phaseh-wave-operator-constant-parity` 52/52.

## Explicitly staged to Phase 1b (documented in-code, not yet implemented)

- **Deterministic serial-mesh rank-0 clustering** (rank-count-independent ids).
- **LTS-aware multi-constraint METIS partition** (companion `.cpp`).
- **Fault-QP reorder** (single-source, P-006).

Reason: these need the material constructed *before* the ParMesh (it is built
after) plus a serial↔local element map, and their acceptance gate (`lts="rate2"`
stations byte-identical under the QP permutation) is only validatable on the
Frontera/Expanse production meshes (no-local-reproducer constraint). Until then
the run-path clustering is built at **np==1 only** (a rank's local mesh can be
disconnected under partition, and the maxdiff fixpoint is per-rank); np>1 logs
the deferral and runs GTS (byte-identical to `lts="off"`).

## Adversarial review (3 reviewers × per-finding verification)

14 raw findings → **1 confirmed**, 13 refuted (each verified against source).

| ID | Sev | Status | Resolution |
|---|---|---|---|
| D-1 | LOW | **CONFIRMED** | `--lts-report` at np>1 could `MFEM_ABORT` in the contiguity assert on a disconnected local submesh. **Fixed:** gate the report's clustering on `nprocs==1` (mirrors the run-path guard); np>1 skips with guidance. |
| C-1 / L-1 | MOD / CRIT | refuted | Same contiguity concern, but unreachable given the np==1 caller gating (single SAFS body is connected); D-1's fix closes the last reachable path. |
| L-3 | LOW | refuted | Applied anyway (defensive): `BuildTickTable` asserts `num_clusters<=62` so `1LL<<c` cannot hit signed-shift UB. |
| L-6 | LOW | refuted | Applied anyway (defensive): `BuildLtsLayout` fails loud on a fault face with `elem2<0` instead of silently demoting it to Boundary (Phase-4 shared-fault concern). |
| C-2 | MOD | refuted | Test strengthened: new **T8** asserts the cost gate binds before nc_cap and is cumulative-vs-baseline (not per-step). |
| L-5 | LOW | refuted | `n_collectives` split into correct/predict sets is a *documented, correct* refinement of the plan's one-line formula for the P-001 split predicates. |
| D-3 | LOW | refuted | New virtuals reorder the vtable — the known clean-rebuild ABI hygiene item, not a code defect. |
| C-3, C-4, C-5, D-2, L-2, L-4 | LOW/MOD | refuted | Test-coverage/perf notes on unreachable or already-invariant-checked paths; verified not defects. |

Post-fix validation: `test_lts_clustering` 4043/4043, `test_lts_layout` 87/87,
`test_spatial_friction_config` 314/314, `test-ader-tpv102-smoke` 4/4, full driver
relinks.

## Commits (local branch `safs-v4_0_0-alt-case1-mfem-speed`)

```
9567808 lts review fixes: D-1 (confirmed) + L-3/L-6/C-2 (defensive)
bfb36b7 lts Phase 1b: run-path clustering + layout wiring (gated, GTS)
fe11ad4 lts Phase 1a: run-side layout + tick table + config parse (A.2/A.3)
3300c80 lts Phase 0: --lts-report driver hook + per-element CFL accessors
f5d007a lts Phase 0: BuildLtsClustering + test_lts_clustering (Appendix A.1/B.1)
```
