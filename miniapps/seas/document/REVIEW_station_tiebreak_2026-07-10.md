# Code Review: station tie-break fix (np≥4 attractor dissolution), 2026-07-10

## Review Scope
- Plan: `debug_document/tpv104_debug_document/np4_attractor_root_cause_2026-07-10.md` §"Fix proposal"
- Files reviewed: `dynamic/station_nearest_tiebreak.hpp` (new),
  `dynamic/tpv{104,205,31,102}_stations.hpp`, `tests/parallel/test_station_tiebreak_np2.cpp`
  (new), `Makefile` (DYNAMIC_HEADERS / obj rule / targets / test-parallel),
  `tests/unit/test_shared_fault_reconcile_cross_rank.cpp` (tolerance override),
  doc edits (R1601 RESOLVED banner, ADER_ITERATOR caveat, unify-plan banner).
- Domain context: `miniapps/seas/CLAUDE.md`; unify plan; the np4 root-cause doc;
  R4-007 history in the writers; memory notes (stale-ABI SEAS_HEADERS trap,
  MUMPS local-env limitation).

## Verified non-issues (measurements/proofs done by this review — do not redo)

1. **Masked lexicographic reduction correctness.** Lexicographic (x,z,y) order is total;
   each rank submits its local minimum, and MIN over per-rank minima = global minimum.
   The sequential masked rounds compare z only among bit-equal x-winners and y only among
   bit-equal (x,z)-winners — correct. RMAX masking cannot collide with a physical
   coordinate; `real_t=float` builds use `numeric_limits<real_t>::max()` consistently.
2. **Within-gate ⊆ local-window containment.** For any candidate rank (gate:
   `d_best_local − G < tol(G)`), every QP q with `d_q − G < tol(G)` satisfies
   `d_q − d_best_local < tol(G) ≤ StationTieTol(d_best_local)` (monotone tolerance,
   `d_best_local ≥ G`), so the local pass-2 window always CONTAINS the rank's true
   global-gate tie subset — the submitted best is never lexicographically larger than the
   rank's true candidate minimum. (The converse overshoot is R-201.)
3. **Collective consistency.** All four writers call `ResolveStationOwnerLex`
   unconditionally for every station (no early `continue` on any rank-local condition);
   4 Allreduce per station on every rank.
4. **Pre-fix failure of the new test.** Traced on the pre-fix code: both ranks tie
   (bit-equal distances), `winning_rank = MIN(0,1) = 0`, rank 0's local pick is its own
   tie member B ⇒ phase 1 reads 222 ≠ 111 ⇒ exit 1. The test genuinely guards the fix.
5. **ofstream reuse / NaN Bcast.** Each scenario uses a distinct prefix; `ofstream::open`
   truncates on re-runs; `MPI_Bcast` of a NaN payload is bit-transport, and NaN
   propagates to a failed `==` check (fails loudly, not silently).
6. **Reconcile-test tolerance 10.0 provably disables the assertion there.** For the
   check `|a−b| ≤ tol·max(|a|,|b|,1)`: `|a−b| ≤ |a|+|b| ≤ 2·max(|a|,|b|,1)` ⇒ the
   relative measure is ≤ 2 for ANY finite a,b ⇒ any tol > 2 disables. No other assertion
   in that test reads the tolerance; the R-101 verifier (its actual oracle) is
   independent, and its 4/4 result includes the negative "guard trips" phase.
7. **TPV26-smoke changed picks are the lexicographic winners.** Station (−7.5, 7.5):
   tie candidates (x,z) = (−1391.75, −1454.25) and (−1454.25, −1391.75), both at
   d = 8594.2939 m ⇒ winner x = −1454.25. Station (+7.5, 7.5): candidates
   (1391.75, −1454.25) / (1454.25, −1391.75) ⇒ winner x = +1391.75. Note the two
   winners are NOT mirror images (lex-min favors smaller x on both sides) — the
   apparent "value swap" (2.554e7 ↔ 2.444e7 n-stress) is the expected consequence, not
   a symmetry bug.
8. **All-ranks-zero-fault corner is production-unreachable** in the native tpv104 driver
   (`MFEM_VERIFY(n_fault_g > 0)` at `tpv104_driver.cpp:1296` precedes writer Open); the
   changed `is_candidate` gate therefore cannot abort a run that previously "worked".
   (Message clarity folded into R-205.)
9. **Stale-reference sweep**: memory file and the unify-plan CURRENT banner updated; the
   remaining "attractor open" texts live in dated historical sections (Phase-2+3 banner,
   two dated REVIEW files) — see R-203 for the one that needs a pointer.

## Findings

### [R-201] [MODERATE] [POSSIBLE] [dynamic/station_nearest_tiebreak.hpp:FindNearestFaultDOFLex + ResolveStationOwnerLex] — local tie window is anchored at the rank-local best, leaving a ~[tol, 2·tol) osculating band where the pick is still partition-dependent

**Category:** EDGE_CASE

**Description:**
The pass-2 window is `d − d_best_local < StationTieTol(d_best_local)`. A candidate
rank's window can therefore admit QPs whose distance exceeds the global gate: up to
`G + 2·tol` when `d_best_local = G + tol⁻`. A QP q* in the band `(G + tol, G + 2·tol)`
with lexicographically smaller coordinates is submitted (and can win) ONLY when it is
co-resident with a `d ≈ G` QP on the same rank; on a partition where q* sits alone on a
rank, its rank fails the gate and q* is never considered. The winner then differs across
partitions — the exact defect class this fix retires, surviving in a distance band of
width `tol ≈ max(1e-10, 1e-9·G)` (sub-micrometer for typical G). Real meshes produce tie
sets that are bit-equal (symmetry) or separated by mesh-scale distances, so this is not
observable in any current configuration — hence POSSIBLE — but it is cheaply closable
and the closure also simplifies the correctness argument to one line ("candidates ≡ the
true global tie set").

**Trigger:**
Synthetic: station S; QP A at distance G (lex-larger coords); QP q* at distance
`G + 1.5·StationTieTol(G)` (lex-smaller). Layout 1: {A, q*} co-resident → q* reported.
Layout 2: A and q* on different ranks → A reported.

**Actual behavior:** pick depends on co-residency inside the osculating band.

**Expected behavior:** pick = lexicographic minimum over the EXACT global tie set
`{q : d_q − G < tol(G)}`, independent of layout.

**Suggested fix:** anchored re-scan protocol — the MPI Open already has
`global_min_dist[s]` in hand before ownership is resolved:
```diff
 // station_nearest_tiebreak.hpp
-inline int FindNearestFaultDOFLex(const std::vector<Vector> &fault_coords,
-                                  int ndof,
-                                  real_t along_strike,
-                                  real_t down_dip)
+/// anchor_dist < 0 (default): window anchored at the local best (serial path).
+/// anchor_dist >= 0: window anchored at |d - anchor_dist| < window_tol (MPI
+/// re-scan with the GLOBAL min); returns -1 when no local QP is in the window.
+inline int FindNearestFaultDOFLex(const std::vector<Vector> &fault_coords,
+                                  int ndof,
+                                  real_t along_strike,
+                                  real_t down_dip,
+                                  real_t anchor_dist = static_cast<real_t>(-1.0),
+                                  real_t window_tol  = static_cast<real_t>(-1.0))
```
with pass 2 using `(anchor_dist >= 0 ? anchor_dist : d_best)` as the window anchor
(`window_tol` likewise, defaulting to `StationTieTol(anchor)`), and each MPI `Open`
replacing its candidate logic with:
```diff
-         const bool is_candidate = std::abs(local_dist[s]
-                                            - global_min_dist[s]) < tie_tol
-                                   && station_dof_[s] >= 0
-                                   && station_dof_[s] < ndof;
+         // Anchored re-scan: candidates ≡ ranks owning a QP in the TRUE
+         // global tie set; the submitted pick is the local lex-min over
+         // exactly that set (R-201).
+         const int lex_dof = FindNearestFaultDOFLex(
+            fault_coords, ndof, stations[s].along_strike,
+            stations[s].down_dip, global_min_dist[s], tie_tol);
+         if (lex_dof >= 0) { station_dof_[s] = lex_dof; }
+         const bool is_candidate = lex_dof >= 0;
```
(tpv102 passes its historical absolute `1e-10` as `window_tol`; the other three pass
their R4-007 `tie_tol`.) The winning rank's `station_dof_[s]` is then the re-scan pick
by construction.

**Test case:**
```cpp
// extend test_station_tiebreak_np2.cpp — phase 5 (osculating band):
// G = sqrt(2)*128 from station "tie"; place q* so that
//   d(q*) = G + 1.5*StationTieTol(G)   (compute z from the closed form)
// with q*.x < A.x.  Layout A: rank 0 = {A, q*}, rank 1 = {far}.
// Layout B: rank 0 = {A}, rank 1 = {q*, far}.
// Assert the reported slip2 is IDENTICAL in both layouts (post-fix: A in
// both, because q* is outside the tol window; or q* in both if inside —
// either way EQUAL).  Pre-R-201-fix, layout A reports q*, layout B reports A.
```

---

### [R-202] [MODERATE] [Makefile] — `*_stations.hpp` are missing from the driver header dependency lists (stale-ABI trap now armed on freshly-edited files)

**Category:** BUG (build system)

**Description:**
`TPV104_HEADERS` (Makefile:5316) lists `tpv104_setup.hpp` but NOT
`dynamic/tpv104_stations.hpp`; make does not chase nested includes, so an edit to a
station header does NOT rebuild the native tpv104 driver object (same for the
tpv205/tpv31/tpv102 driver lists — the station headers appear only in the spatial-driver
obj rule at :2887-2896 and the new test rule). This round dodged the trap only because
`station_nearest_tiebreak.hpp` (in `DYNAMIC_HEADERS`) changed simultaneously. The next
edit that touches ONLY a `*_stations.hpp` ships a stale driver — the exact
"SEAS_HEADERS stale-ABI" failure mode already in the project memory. Pre-existing gap,
but this change made the station headers hot files.

**Trigger:** `touch dynamic/tpv104_stations.hpp && make seas_tpv104_driver` → "up to
date" (no rebuild) → binary runs old writer code.

**Actual behavior:** stale driver binaries after station-writer edits.

**Expected behavior:** any driver embedding a station writer rebuilds when its writer
header changes.

**Suggested fix:**
```diff
 DYNAMIC_HEADERS = dynamic/wave_state.hpp \
                   ...
                   dynamic/station_nearest_tiebreak.hpp \
+                  dynamic/tpv102_stations.hpp \
+                  dynamic/tpv104_stations.hpp \
+                  dynamic/tpv205_stations.hpp \
+                  dynamic/tpv31_stations.hpp \
                   dynamic/tpv102_setup.hpp
```
(One list, covers every `$(SEAS_HEADERS)`-dependent object; slightly over-broad
rebuilds are the accepted cost of this Makefile's explicit-list scheme.)

**Test case:**
```bash
make seas_tpv104_driver ... && touch dynamic/tpv104_stations.hpp && \
  make -q seas_tpv104_driver ...; test $? -ne 0   # -q must report "needs rebuild"
```

---

### [R-203] [LOW] [document/PLAN_unify_interior_shared_fault_substep_2026-07-09.md:132] — the historical Phase-2+3 banner still says the attractor is "Still open" with the falsified suspect

**Category:** QUALITY (doc)

**Description:**
The (dated, historical) "PHASE 2 + PHASE 3 IMPLEMENTED" banner retains: "Still open:
the np≥4 attractor … Prime suspect is now sharpened … This is the Phase 5 collapse's
problem to solve". Both the openness and the suspect are now falsified. The CURRENT
banner (:101) is correct, but a reader scanning chronologically can stop at :132.
Do not rewrite history — add a supersession pointer.

**Suggested fix:**
```diff
-- **Still open:** the **np≥4 attractor** (station strike-slip 4.8123 vs 5.0022 at np≥4; now
+- **Still open** *(superseded 2026-07-10: RESOLVED as a station-READOUT tie-break artifact —
+  the frame-orientation suspect below was falsified; see np4_attractor_root_cause_2026-07-10.md
+  and the Phase-4+5 banner above)*: the **np≥4 attractor** (station strike-slip 4.8123 vs 5.0022 at np≥4; now
```

---

### [R-204] [LOW] [debug_document/tpv104_debug_document/np4_attractor_root_cause_2026-07-10.md] — the TPV26-smoke np=1 trace change is not recorded anywhere durable

**Category:** QUALITY (doc)

**Description:**
The fix changed two np=1 TPV26-smoke station files (`tpv205_station_x2_±7.5_x3_7.5.dat`,
n-stress 2.554e7 ↔ 2.444e7 apparent swap). The justification (exact corner-QP ties at
d = 8594.2939 m; new picks = lex winners x = −1454.25 / +1391.75; winners are not mirror
images by design) exists only in the session conversation. A future `git diff` reader of
those output files (or of a regression harness pinning them) finds no explanation.

**Suggested fix:** append to the np4 doc's "Fix proposal" section:
```markdown
### Post-fix trace changes (recorded 2026-07-10)
- tpv104 symmirror: np=1 traces UNCHANGED (the serial pick already was the lex winner);
  np ∈ {2,4,7,10} now byte-match np=1 at all nine stations (worst 1e-12).
- TPV26 smoke (np=1): stations (±7.5 km, 7.5 km) are exact ties between the smoke
  fault's two corner QPs (d = 8594.2939 m; the stations sit 8.6 km off the tiny scaled
  fault).  New picks are the lexicographic winners (x = −1454.25 m and x = +1391.75 m);
  n-stress column changes 2.554e7 ↔ 2.444e7.  These stations sample a far-field corner
  and carry no benchmark meaning for the smoke; V_max gates unchanged
  (0.185543 / 0.199649).
```

---

### [R-205] [LOW] [dynamic/station_nearest_tiebreak.hpp + tests/unit/test_shared_fault_reconcile_cross_rank.cpp] — two comment gaps that will cost a future debugger time

**Category:** QUALITY

**Description:**
(a) `ResolveStationOwnerLex`'s VERIFY message ("no distance candidate submitted
coordinates; check the candidate gate / global_min_dist") does not name the most likely
real cause — zero fault QPs on every rank (config/mesh error) — nor the NaN-coordinate
corner (NaN fails every `==`, aborting here rather than upstream).
(b) The reconcile test's tolerance comment says "raise the tolerance above the injected
seed" — the actual guarantee is stronger and worth stating: the relative measure is
bounded by 2 for any finite pair, so tol = 10 provably disables the assertion for ALL
fields, present and future.

**Suggested fix:**
```diff
-               "coordinates; check the candidate gate / global_min_dist).");
+               "coordinates).  Most likely: no rank owns any fault QP "
+               "(fault-less mesh / wrong fault attribute), or a NaN fault-QP "
+               "coordinate (NaN fails the coordinate-equality rounds).");
```
```diff
-   // above the injected seed for this test only — the assertion-fires guard
+   // above the injected seed for this test only (the check's relative measure
+   // |a-b|/max(|a|,|b|,1) is bounded by 2 for any finite pair, so tol=10
+   // provably disables it for every payload field) — the assertion-fires guard
```

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-201 POSSIBLE osculating-band partition dependence — closable via
  anchored re-scan; R-202 build-dependency gap on the now-hot station headers)
- Low issues: 3 (R-203 superseded-banner pointer, R-204 durable record of the TPV26-smoke
  trace change, R-205 comment gaps)
- Plan compliance: FULL (fix proposal items 1–4 implemented; the shared-helper deviation
  and the tpv102-gate preservation are documented and sound; acceptance evidence
  regenerated and stronger than required — np=1 symmirror unchanged)
- Verdict: **PASS WITH FIXES**

## Fix Status (applied 2026-07-10, same session)

| ID | Severity | Status | What was done |
|---|---|---|---|
| R-201 | MODERATE (POSSIBLE) | **FIXED** | `FindNearestFaultDOFLex` gained `anchor_dist`/`window_tol` (anchor<0 = historical serial path, verified unchanged); all four MPI `Open`s use the global-min-anchored re-scan (tpv102 passes its absolute 1e-10 gate); test phase 5 (osculating band, two layouts) added — 8/8. Byte-neutral on real runs: symmirror np1+np4 (18 files) and TPV26 smoke (16 files) identical to the post-implement references. |
| R-202 | MODERATE | **FIXED** | All four `*_stations.hpp` added to `DYNAMIC_HEADERS`; probe: `touch tpv104_stations.hpp && make -q seas_tpv104_driver` exits non-zero (rebuild scheduled). |
| R-203 | LOW | **FIXED** | Supersession pointer inserted at the Phase-2+3 banner's "Still open" bullet. |
| R-204 | LOW | **FIXED** | "Post-fix trace changes" record appended to the np4 root-cause doc (symmirror unchanged at np=1; TPV26-smoke corner-tie changes with values, winners, and the not-mirror-images note). |
| R-205 | LOW | **FIXED** | VERIFY message names the fault-less-mesh / NaN causes; reconcile-test comment states the rel≤2 bound that makes tol=10 a provable disable. |

**Post-fix verification (all green):** tie-break 8/8 (incl. phase 5) at np=2; reconcile
cross-rank 4/4; FR parity 68/68; substep parity 39/39; TPV26 smoke V_max 0.185543;
byte-neutrality per R-201 row above.

## Unreviewed Areas
- The three MUMPS-failing QD parallel targets (`test-parallel-domain`,
  `test-mms-parallel`, `test-parallel-fault`) — crash inside `libdmumps` on this
  laptop's mfem-dev env (known local limitation, memory-documented); no changed code in
  their stack. Not attributable to this change; needs the `mumps-test` env or a cluster
  run to go green locally.
- `TPV104SurfaceStationWriter` — different mechanism (`Mesh::FindPoints`); no nearest-QP
  tie-break. It has its own latent multi-rank double-claim risk (every rank that finds
  the point opens the same file), pre-existing and out of scope; noted for a future
  review.
- The symmirror np-matrix and smoke byte-comparisons — verified in the implement leg
  this session; not re-run per the review instructions.
