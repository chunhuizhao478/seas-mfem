# The "np≥4 attractor" root cause — a station-readout tie-break artifact, not a dynamics defect (2026-07-10)

**Status: ROOT-CAUSED, VERIFIED.** The "binary np≥4 trajectory switch" that gated unify-plan
Phase 6 (PLAN_unify_interior_shared_fault_substep_2026-07-09.md, "Remaining" banner;
R1601_root_cause_2026-07-10.md §"NEW OPEN FINDING") **is not a trajectory switch at all**.
The simulated fields are **bit-exactly partition-invariant** across every partition tested
(np ∈ {1..10} plus a fault-locality partition — 11 partitions, 9 stations, all agree to
≤ 1e-14 except the artifact below). What flips at "np≥4" is *which fault QP the station
writer samples*: stations at along-strike x2 = 0 sit exactly equidistant between two
quadrature points 400 m apart, and `TPV104StationWriter`'s nearest-DOF tie-breaking is
partition-dependent. The "4% attractor gap" is the true spatial difference of the (single,
correct, partition-invariant) solution between those two physical points.

## Symptom

- Observed: hypocenter-station finals (t = 2 s, `seas_tpv104_driver`, symmirror 1000 m mesh,
  P=1/O=2, `--fault-iterator substep --fric-law slip-srw`) split into two "branches":
  strike-slip 5.0021709132 / dip −3.638e-02 at np ∈ {1,2}, vs 4.8122895985 / +1.169e-01 at
  np ∈ {4,10}; each branch internally reproducible to ≥10 digits; V_max identical (12.7387).
- Expected: partition-independent station traces.
- Conditions: previously believed "np ≥ 4"; this doc shows the np-count framing was wrong.
- Recent changes: none causal — the retired R-1601 fallback showed the same "jump" (4.8146);
  the unified path (Phases 2–5) reproduces it *more* cleanly. Pre-existing behavior.

## Reproduction

```
mpirun -np <N> ./seas_tpv104_driver --mesh tpv104/mesh/tpv104_symmirror_1000m.msh \
  --order 1 --tfinal 2.0 --ader-order 2 --fault-iterator substep --fric-law slip-srw \
  --output-dir <dir> --output-prefix sym
# metric: last line of sym_station_x2_0_x3_7.5.dat, cols 2 (h-slip) and 5 (v-slip)
```
Deterministic for every N (re-runs byte-identical).

## Hypotheses and their fates

| # | Hypothesis | Fate | Killing evidence |
|---|---|---|---|
| H3a | np=2 lands on np=1 only because METIS produces no fault seam at np≤3 | TRUE but not the mechanism | np=2/3 have 0 cut fault faces (E1) — but see H3b |
| H3b | The switch = existence of ≥1 shared fault face | **ELIMINATED** | np=6 (9 cut faces) and np=7 (12 cut faces) reproduce np=1 **to 8e-29** (E1); np=4 **with `--partition-fault-locality`, i.e. 0 cut faces, lands on "B"** (E5) |
| H1 | Interior-vs-shared frame-bit flip (`!elem1_on_plus` vs `sign_flipped`) alters seam physics | **ELIMINATED** | same E1/E5 facts: partitions with 27–36 seam QP pairs are bit-exact vs np=1; a seam with 0 cut faces flips the "branch" |
| H2 | Ranks owning zero fault faces (R1601-doc suspect) | **ELIMINATED** | partition dump (E3): every rank owns ≥1 fault face at every np ∈ {2..10} |
| — | Dust-seeded bistability of the y-mirror dip mode | **ELIMINATED** | 1 mm node perturbation (seed ~1e-6 ≫ 1e-16 dust) leaves the partition→"branch" map exactly intact and moves finals only in the 9th digit (E4) |
| — | Partition y-mirror-symmetry selects the branch | **ELIMINATED** | no METIS partition is mirror-symmetric; violation counts (18…311) uncorrelated with branch (E3) |
| ✔ | **Station nearest-QP tie-break is partition-dependent** | **CONFIRMED** | E2 (station × partition matrix), E6 (geometric tie), code reading (E7) |

## Evidence chain

### E1 — np sweep with seam census
`np ∈ {1..10}`: cut-fault-face QP pairs = {0,0,3,15,27,36,69,60,63} for np=2..10.
Hypocenter finals: np ∈ {1,2,3,6,7} → 5.0021709132/−3.638e-02; np ∈ {4,5,8,9,10} →
4.8122895985/+1.169e-01. **No correlation with seam count** (np=7: 36 pairs → "A";
np=4: 3 pairs → "B").

### E2 — the station × partition difference matrix (the decisive measurement)
Max |trace(np) − trace(np=1)| over the full 522-step series, per station:

| station | np2 | np3 | np4 | np4+G1 | np5 | np6 | np7 | np8 | np9 | np10 |
|---|---|---|---|---|---|---|---|---|---|---|
| x2_−12_x3_12 | 0 | 0 | 0 | 0 | 1e−14 | 1e−14 | 1e−14 | 1e−14 | 0 | 1e−14 |
| x2_−12_x3_3 | 0 | 0 | 0 | 0 | 6e−30 | 0 | 0 | 0 | 0 | 1e−29 |
| x2_−9_x3_7.5 | 0 | 0 | 0 | 0 | 1e−17 | 0 | 0 | 0 | 0 | 0 |
| x2_0_x3_12 | 0 | 1e−14 | 6e−30 | 6e−30 | 0 | 0 | 6e−30 | 6e−30 | 0 | 0 |
| **x2_0_x3_3** | **6.7e5** | **6.7e5** | 0 | 1e−32 | 1e−32 | 1e−32 | 8e−29 | 8e−29 | **3.6e6** | 7e−29 |
| **x2_0_x3_7.5** | 1e−12 | 0 | **3.4e6** | **3.4e6** | **3.4e6** | 1e−12 | 0 | **3.4e6** | **3.4e6** | **3.4e6** |
| x2_12_x3_12 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 5e−30 | 0 | 0 |
| x2_12_x3_3 | 8e−29 | 3e−30 | 8e−29 | 8e−29 | 8e−29 | 2e−30 | 0 | 0 | 0 | 0 |
| x2_9_x3_7.5 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Readings:
1. **Seven of nine stations agree across all 11 partitions to ≤1e-14** (the 1e-14/e-29/e-30
   entries are denormal dust in the pre-rupture dip columns and one print-quantum). The
   dynamics are bit-exactly partition-invariant — on a mesh whose partitions cut the fault
   with up to 69 shared QP pairs, through the full unified substep/reconcile path.
2. Only the two **x2 = 0** stations differ macroscopically — and *each has its own partition
   membership set* (x2_0_x3_3 flips for {2,3,9}; x2_0_x3_7.5 flips for {4,4+G1,5,8,9,10}).
   A genuine trajectory switch would move all stations coherently. Per-station flips with
   unrelated membership sets are only consistent with per-station readout selection.
3. The x2_0_x3_3 "flip" is a persistent 1.7% shear-stress offset (6.7e5 Pa at t=2 s) while
   the same run's hypocenter matches np=1 to 11 digits — same artifact, other station.

### E3 — partition structure dump (test-only tool, `$JOB_TMP/dump_part.cpp`)
`Mesh::GeneratePartitioning(np, 1)` (the exact ParMesh-ctor path) reproduces the observed
cut-QP counts for every np (validates the tool). Per-rank fault-face census: **no rank owns
zero fault faces at any np** (kills H2). No partition is y-mirror-symmetric under any rank
relabeling (kills the partition-symmetry hypothesis).

### E4 — mesh-perturbation control
A copy of the mesh with ONE off-fault node moved 1 mm in y (node 32 at (−1000,−1000,−1000)):
partition→branch map **unchanged** ({1,7}→A, {4,5,10}→B), finals move only ~1e-9.
Kills every FP-dust/bistability story: a 1e-6-scale physical asymmetry cannot flip what the
partition trivially flips.

### E5 — fault-locality control (the reframing experiment)
`mpirun -np 4 ... --partition-fault-locality` (G1: every fault face's element pair
co-resident; **zero** shared fault QPs, no reconcile, no seam physics at all) → lands on
"B" (4.8122895985/+1.1691568e-01) bit-exactly. The "branch" is decoupled from fault seams
entirely.

### E6 — the geometric tie
Fault QPs (3-pt triangle rule) nearest to station (x2=0, down-dip 7500 m):
`(x=−166.67 m, z=−7666.7 m)` and `(x=+166.67 m, z=−7333.3 m)`, both at
dist² = 55555.55… m² — equidistant by the mesh's x-symmetry. (Same structure at
(0, 3000): the other flipping station. The x2=±9, ±12 km stations have unique nearest QPs
— exactly the seven stations that never flip.)

### E7 — the code
`dynamic/tpv104_stations.hpp`:
- `FindNearestDOF_TPV104` (:72-92): linear scan over the **rank-local** `fault_coords` in
  local dof order with strict `dist2 < best_dist` — among (near-)tied QPs the **first in
  local enumeration** wins. Local enumeration order is partition-dependent.
- MPI `Open(... comm)` (:130-198): per-rank best distance → `MPI_Allreduce(MIN)` → every
  rank within `tie_tol = max(1e-10, 1e-9·max(dist,1))` (≈2.4e-7 here, ≫ any FP split of the
  exact tie) is a candidate → **winning rank = lowest candidate rank** (second
  `Allreduce(MIN)` on rank id). Net: the station reports whichever tied QP happens to live
  on the lowest-numbered rank owning one, as picked by that rank's local scan order.

At np=1 the serial scan picks one tie member (the 5.0021709132 point); other partitions pick
either member depending on rank layout. Both values are correct samples of the same field —
5.0021709132 and 4.8122895985 are the strike-slip of the SAME solution at two points 400 m
apart (and −3.638e-02 vs +1.169e-01 are the known P=1 dip-pollution values at those points;
the dip channel varies strongly near the hypocenter, which is why the "dip drift" looked so
alarming).

Same pattern (copy-derived) in: `tpv102_stations.hpp`, `tpv205_stations.hpp` (used by the
spatial dyn driver's fault stations), `tpv31_stations.hpp`, and the TPV104 surface-station
writer. Any station coordinate landing on a symmetry line / element boundary of any mesh is
exposed.

## Root Cause

- **What:** partition-dependent tie-breaking in station→fault-QP assignment: local-scan-order
  pick inside `FindNearestDOF_*` + lowest-candidate-rank pick in the MPI owner resolution.
  For stations equidistant from ≥2 fault QPs (x2=0 stations on x-symmetric meshes), different
  partitions report different — individually correct — sample points.
- **Where:** `dynamic/tpv104_stations.hpp:72-92` and `:158-196` (and the three sibling
  writers listed above).
- **Why:** the writer was designed assuming a unique nearest QP; `tie_tol` (R4-007) was added
  to make *ownership* robust, which correctly keeps ties from dropping stations but leaves
  *which tied QP* unspecified — it inherits rank numbering and local ordering.
- **Impact:** every np-vs-np or run-vs-run comparison keyed on an x2=0 station of a symmetric
  mesh is confounded: the R1601 doc's "NEW OPEN FINDING" table, and the hypocenter columns of
  `ADER_ITERATOR_SYM1000_RESULTS.md`, mix the two sample points across rows. **No dynamics
  defect exists**: the unify-plan's np≥4 blocker dissolves, and Phase 6's "metric must
  decrease" acceptance is measurable once the readout is deterministic.

## Fix proposal (for the follow-up /code-implement leg — NOT applied in this diagnosis)

Deterministic, partition-invariant tie-breaking by physical coordinates:
1. In `FindNearestDOF_*`: replace the strict `dist2 < best_dist` scan with a lexicographic
   key `(dist2, x, z, y)` — among distance-ties (within a small relative epsilon or exact),
   prefer smaller x, then smaller z (deterministic geometry, independent of enumeration
   order).
2. In the MPI `Open`: after `global_min_dist`, candidates reduce their best QP's `(x, z)`
   lexicographically (two-phase Allreduce: MIN on x among distance-candidates, then MIN on z
   among x-winners) and ownership goes to the rank holding the winning COORDINATE (rank id
   only as the final tie among bit-equal coordinates, i.e. the duplicated shared-QP copies).
3. Apply identically to the four sibling writers; np=1 traces must be checked — if the serial
   scan already picks the lexicographic winner (likely: smaller x first in structured
   enumeration), np=1 output is byte-unchanged; if not, np=1 gold traces for symmirror change
   ONCE, with this doc as justification, and then all np match np=1 bit-exactly.
4. Regression guard: a station placed at an exact tie on the 2-tet/tet2x2 fixture (or
   symmirror at np ∈ {1,4}) asserting identical trace across np.

### Post-fix trace changes (recorded 2026-07-10; fix = `dynamic/station_nearest_tiebreak.hpp` + all four writers, guard = `tests/parallel/test_station_tiebreak_np2.cpp`)

- **tpv104 symmirror: np=1 traces UNCHANGED** (the serial pick already was the lexicographic
  winner at both tied stations); np ∈ {2,4,7,10} now match np=1 at all nine stations
  (worst 1e-12, one print-quantum value).
- **TPV26 smoke (np=1): two station files changed** — stations (±7.5 km, 7.5 km) are exact
  ties between the smoke fault's two corner QPs (d = 8594.2939 m; the stations sit 8.6 km
  off the tiny scaled fault, extent ±1.45 km).  New picks are the lexicographic winners
  (x = −1454.25 m at the −7.5 station, x = +1391.75 m at the +7.5 station — NOT mirror
  images: lex-min favors smaller x on both sides, by design).  The n-stress column changes
  2.554e7 ↔ 2.444e7 (a 4% apparent "swap").  These stations sample a far-field corner and
  carry no benchmark meaning for the smoke; the byte-exact V_max gates are unchanged
  (TPV26 0.185543, forced-rupture 0.199649).  Post-fix, the forced-rupture smoke's np=4
  station files byte-match np=1 at 15/16 stations (the residual is the pre-existing
  1-row/last-printed-digit FP dust at x2_12_x3_7.5 — not a tie).

## Verification of this diagnosis

- Reproduction: deterministic, all runs; matrix E2 regenerated from scratch this session.
- No production code was modified: all instrumentation was offline (partition dumper +
  python analysis in the job tmp dir); the runs used the stock driver.
- Test suite: untouched (no code changes made during diagnosis).

## Implications for unify-plan Phase 6

1. The **np≥4 attractor blocker is dissolved** — there is no partition-dependent dynamics
   defect on this problem. The unified substep path is bit-exactly partition-invariant at
   every non-tied station across 11 partitions (the strongest multi-rank validation this
   code has ever had).
2. Phase 6 step 1's before/after metric at x2=0 stations must be computed **after** the
   tie-break fix (or on non-tied stations); historical hypocenter tables need a confounding
   caveat, not rewriting.
3. np>1 gold regeneration may become unnecessary: with the readout fixed, np>1 traces are
   expected to equal np=1 bit-exactly on this problem. Sign-off question changes from
   "accept new np>1 numbers" to "accept the one-time np=1 hypocenter-trace change IF the
   serial pick differs from the lexicographic winner".
