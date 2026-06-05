# Mixed-Flux + Heterogeneous (Bi-material) Riemann Solver — Plan Review, Bug Report & Proposed Unit Tests

## Summary

Pre-implementation review of `PLAN_mixed_flux_hetero_riemann.md` (companion:
`EXPLORE_mixed_flux_hetero_riemann.md`) against the existing source under
`miniapps/seas/`. **No source has been written yet** (`git status`: only the
markdown docs + the moved PDF are new/modified) — this review verifies the plan's
factual claims about the current code, its internal consistency, mathematical
correctness, and completeness.

Files reviewed (read-only, as-is):
- `dynamic/godunov_flux_bimaterial.hpp` / `.cpp` (Phase 1 primitive site)
- `dynamic/bimaterial_wave_operator.hpp` / `.inl` (Phases 2,3,4)
- `dynamic/wave_operator.hpp` / `.inl` (base hooks, CFL switch, abort)
- `dynamic/godunov_flux.hpp` / `.cpp` (`Central`, accessors, frame helpers)
- `drivers/spatial_dyn_driver.cpp` (Phase 5 wiring + residual guards; Phase 7 CLI)
- `tests/unit/test_bimaterial_wave_operator_parity.cpp` (equivalence reference)
- `makefile` (test target group)
- **(Rev 2)** `tpv31/configs/tpv31.toml`, `tpv31/visualize_results.py`,
  `tpv31/benchmark_data/README.md`, `tpv31/benchmark_data/scec_seisol/*.txt`,
  `jobs/tpv104_spatial/tpv104_p2_rk45_mixedflux_12N_600r_normal.sbatch`,
  `jobs/tpv31_spatial/tpv31_p2_aderO3_noflux_50m_normal.sbatch`,
  the tpv104/tpv31 config + sbatch corpus (Phase 7 premise).

**Review history:**
- **Rev 1** (2026-06-04): 9 plan defects found (3 CRITICAL, 4 WARNING/MEDIUM, 2 LOW).
  Of the user's 10 verification items: 6 fully checked out, 4 revealed plan defects.
- **Rev 2** (2026-06-04): Re-review of the revised plan, which folds the BUG-1…BUG-9
  fixes in and adds **Phase 7 (TPV31 Frontera end-to-end validation)**.
  **BUG-1, BUG-2, BUG-3, BUG-4, BUG-5, BUG-7, BUG-8, BUG-9 are RESOLVED in the plan
  text** (verified against live source — all cited file:line are accurate, code
  unchanged since Rev 1). **BUG-6 is NOT correctly resolved** — the Rev-2 mitigation
  is *unimplementable as written* (its abort condition can never fire); re-classified
  and carried forward, with a new sub-finding. **5 NEW findings** from the Phase-7
  review: 2 CRITICAL (BUG-6 abort dead + the SCEC-tolerance tooling gap), 2 MEDIUM,
  1 LOW. Live code is byte-for-byte unchanged from Rev 1, so all prior code citations
  remain valid.
- **Rev 3** (2026-06-04): Re-review of plan Rev 3, which fixes the five Rev-2 open
  findings (BUG-10..BUG-14). **Verdict: BUG-12, BUG-13, BUG-14 are RESOLVED; BUG-11
  is PARTIALLY RESOLVED (the flags are now in scope, but the metric is under-specified);
  BUG-10's *root cause* is correctly diagnosed (the dead predicate is removed) but its
  *replacement* — the "material-field probe at x_loc/x_nbr" — is NOT implementable as
  written and would false-abort on TPV31's own depth profile.** 6 NEW findings: 2 CRITICAL
  (BUG-15 probe is unimplementable against the real `MaterialField` API; BUG-16 the probe
  false-aborts on legitimate depth variation), 1 HIGH (BUG-19 RMS/peak needs interpolation
  onto a common time base — absent from script and plan), 2 MEDIUM (BUG-17 relative-error
  denominator undefined; BUG-18 gated station column unspecified), 1 LOW (BUG-20 reflection
  geometry/tol under-defined). Live source is STILL byte-for-byte unchanged (`git status`:
  only markdown + the moved PDF), so all prior code citations remain valid; every NEW
  finding is re-verified against the same live `miniapps/seas/` source.
- **Rev 4** (2026-06-04): Re-review of plan Rev 4, which set out to fix the six Rev-3
  open findings (BUG-15/16 CRITICAL, BUG-19 HIGH, BUG-17/18 MEDIUM, BUG-20 LOW) by
  (a) DELETING the material-field probe everywhere and mandating a declarative
  `material.seam_continuous = true` flag as the sole lateral-het guard, and (b) fully
  specifying the SeisSol-overlay metric (interpolate via `np.interp` → peak-normalize →
  named channels). **Verdict on the six: BUG-15, BUG-16, BUG-17, BUG-18, BUG-19, BUG-20
  are ALL RESOLVED** — the probe is gone as a mechanism (referenced only as the rejected
  option), the flag is threaded consistently, and the metric (interp + peak-norm + named
  channels excluding `mu_eff`) is implementable and verified against the real
  `visualize_results.py` column map (the plan's channel names `V_strike`/`slip_strike`/
  `tau_strike`/`sigma_n` MATCH the script's dict keys `:145-155,179-192,205-214` exactly).
  **BUT the Rev-4 fix introduced / left standing THREE NEW defects, one of them a hard
  CRITICAL functional blocker the prior three rounds missed:**
  **BUG-21 (CRITICAL)** — the new `seam_continuous` flag (and the `mixed_flux="adjacent"` TOML
  key it accompanies) trips a SECOND, config-PARSE-level mutual-exclusion guard at
  `spatial/code/spatial_friction.cpp:1162-1167` that the plan NEVER mentions: that validator
  MFEM_ABORTs `interior_flux="matrix"` + `mixed_flux != "none"` inside `LoadSpatialFrictionConfig`
  (driver `:616`), BEFORE the BUG-4 driver guard (`:809-816`), before CLI override (`:679`),
  before operator construction (`:1135`), and before `SetMixedFluxMode` (`:1245`). The Phase-7
  TPV31 config sets `mixed_flux="adjacent"` in the TOML, so it aborts at config load. The plan
  relaxes only ONE of the TWO series guards. **BUG-22 (HIGH/SCOPE)** — the `[material]` struct
  (`MaterialSpec`, `spatial_friction.hpp:573`), its parse (`spatial_friction.cpp:1626+`), AND
  the BUG-21 guard ALL live in `spatial/code/spatial_friction.{hpp,cpp}` — OUTSIDE the plan's
  declared source-edit scope (`dynamic/` + `drivers/spatial_dyn_driver.cpp` + `tests/` +
  `Makefile`). The plan's hand-wave "the config parse lives with the driver; if the config
  struct is in a separate file, confirm it is within the editable tree" (plan `:612-614`) is
  FACTUALLY WRONG — it does not live with the driver — and silently expands the no-touch
  boundary. **BUG-23 (MEDIUM)** — the Phase-7 metric's overlap window `[t0,t1]` and
  partial-coverage handling are unspecified; a truncated MFEM run (48 h wall / failed restart,
  not reaching `tfinal=15s`) interpolated onto the full reference grid is silently clamped by
  `np.interp` and can either spuriously fail or silently pass over a tiny overlap. **One MODERATE
  gap (BUG-24)**: the plan does not state EXPLICITLY that the `seam_continuous` guard leaves the
  existing TPV31 ADER/`mixed_flux=none` job unaffected (it is, structurally — empty
  `central_flux_face_set_` ⇒ guard never reached — but the byte-exact-contract claim should be
  made explicit). Live source is STILL byte-for-byte unchanged (`git status`: only the markdown
  docs + the moved PDF + `reports/.report-config`); every NEW finding is verified against the
  same live `miniapps/seas/` source. BUG-3 (single-valued central deposit) and BUG-4 (driver
  guard) re-confirmed STILL FIXED.
- **Rev 5** (2026-06-04): Re-review of plan Rev 5, which fixes the four Rev-4 open findings
  (BUG-21 CRITICAL, BUG-22 HIGH, BUG-23 MEDIUM, BUG-24 MODERATE) and — at the user's explicit
  request — adds an "exhaustive guard inventory" to Phase 5 claiming exactly THREE series guards
  (G1/G2/G3) block `matrix+adjacent+RK` with four others correctly KEPT. **This is the FOURTH
  consecutive review; each prior fix round introduced or missed a critical the next caught
  (probe -> BUG-15/16; missed G1 -> BUG-21), so the Rev-5 pass INDEPENDENTLY re-ran the guard
  sweep rather than trusting the inventory.** **Verdict on the four: BUG-21, BUG-22, BUG-23,
  BUG-24 are ALL RESOLVED** — every cited file:line re-verified byte-accurate against the
  (still-unchanged) live source: G1 is exactly `spatial_friction.cpp:1162-1167`; `MaterialSpec`
  is at `spatial_friction.hpp:573`; the `[material]` parse block is `spatial_friction.cpp:1626-1638`;
  `toml_bool` really exists (`spatial_friction.cpp:304`, exact signature
  `toml_bool(const toml::value&, const std::string&, bool)`); `BuildCentralFluxFaceSet_` really
  clears `central_flux_face_set_` and returns on `None` (`wave_operator.inl:1652-1653`), so the
  BUG-24 ADER-byte-exactness claim holds by construction. **The "exactly THREE guards" claim is
  INDEPENDENTLY CONFIRMED** by a fresh grep of the WHOLE relevant tree INCLUDING the places the
  author did not sweep: the RK time-stepper (`rk_time_stepper.{hpp,cpp}`), the ADER substep path,
  the `BimaterialWaveOperator`/base ctors, the integrator config validation, and the
  `fault_iterator="substep"`/checkpoint paths. **NO fourth blocking guard exists**, and the four
  "kept" guards are correctly classified (file:line in the Rev-5 report). Two NEW findings, both
  NON-blocking: **BUG-25 (MEDIUM)** — Test 5.3c's parser-regression corpus is mis-worded and
  partly infeasible-as-written: there is NO `operator==`/dump on `SpatialFrictionConfig` (so
  "byte-identical config" must be realized field-by-field, the established `test_tpv_config_parse`
  pattern), AND **BP5 has no spatial-driver TOML** — it never parses through
  `LoadSpatialFrictionConfig` (it uses the native `seas_bp5_full` driver), so it cannot be in the
  parse corpus; its byte-exactness is protected by leaving its driver untouched, not by this test.
  **BUG-26 (LOW)** — the plan is SILENT on TPV31 being `law="slip_weakening"` (LSW), so the actual
  RK path is `AdvanceRKCoupledLSW_Spatial` (not the rate-state `AdvanceRKCoupled_Spatial` the
  Test-5.1/6.1 framing implies); the central dispatch lives in the law-agnostic `Mult`, so the
  feature works for LSW, but the plan should name which stepper the deliverable exercises. **Neither
  NEW finding is a CRITICAL/HIGH blocker.** Live source is STILL byte-for-byte unchanged
  (`git status`: only markdown docs + the moved PDF + `reports/.report-config`; `0` source files
  changed). **No CRITICAL/HIGH finding remains open after Rev 5.**

---

## Rev-2 Verdict — BUG-1…BUG-9 resolution status

| BUG | Rev-1 sev | Rev-2 status | Verified against live source |
|-----|-----------|--------------|------------------------------|
| BUG-1 | MEDIUM | **RESOLVED** | abort body `bimaterial_wave_operator.inl:567-576` ✓; stale comments `inl:518-520` ✓, `hpp:121` ✓, `wave_operator.inl:1623-1626` ✓ — all four cited locations accurate |
| BUG-2 | HIGH | **RESOLVED** | delete `inl:521-531` ✓; apply factor at per-element `dt_e` (`inl:546`) ✓; MPI_Allreduce(MIN) semantics PRESERVED (proof below) |
| BUG-3 | CRITICAL | **RESOLVED** | side-symmetric single pair `std::unordered_map<int,std::array<DenseMatrix,2>>`; deposit once + `F_h_e2[c]=F_h_e1[c]`; mirrors scalar `wave_operator.inl:1203-1209` ✓ (verified) |
| BUG-4 | CRITICAL | **RESOLVED** | guard at `spatial_dyn_driver.cpp:809-816` ✓; rewrite logically sound + consistent with Phase-4 `ComputeMaxDt` guard |
| BUG-5 | MEDIUM | **RESOLVED** | `sf`-loop shared arm mirroring `inl:375-422` ✓; `pmesh.GetSharedFace(sf)` ✓; `shared_face_neighbour_material_[sf]` keyed by `sf` ✓ |
| BUG-6 | MEDIUM | **NOT RESOLVED → escalated** | abort condition "neighbour ≠ local" can **never fire** because `inl:162` stores local-as-neighbour BY CONSTRUCTION → **BUG-10 (CRITICAL/HIGH)** |
| BUG-7 | MEDIUM | **RESOLVED** | single-valuedness moved to operator-level Test 3.2 on heterogeneous material; Phase 1 explicitly drops the trivial primitive test ✓ |
| BUG-8 | LOW | **RESOLVED** | 4-edit recipe; `makefile:311-317` ✓, `859-893` ✓, `1638` ✓, `5248` ✓, `5252` ✓ — all accurate; primitive vs operator obj-dep lists correct |
| BUG-9 | LOW | **RESOLVED** | scratch `tmp1/tmp2`, `SetSize` output only, `cpp:282-295` ✓; acoustic `MFEM_VERIFY(mu>eps)` added; matches `Central` `cpp:416` construction |

### Phase-7 (NEW) verification summary

| Item | Plan claim | Result | Finding |
|------|-----------|--------|---------|
| B.1 | `tpv31.toml`: matrix / none / order=1 / depth_profile_1d / tfinal=15s / cfl=0.5 / cfl_safety=dg | **CONFIRMED** (lines 240,227,187,97,243,228,231) | — |
| B.2 | template sbatch exist; module/LD block byte-identical | **CONFIRMED** (TPV104 71-89, TPV31 60-78) | — |
| B.3 | CLI flags all exist | **CONFIRMED** (all 9 present, 1 each) | — |
| B.4 | `--verify-dispatch` prints "interior flux:", "mixed flux:", "time integrator:" | **PARTLY FALSE** — verify-dispatch prints only "interior flux :"; the other two come from the unconditional banner | BUG-13 |
| B.5 | `visualize_results.py` accepts `--tol-rms`/`--tol-peak` | **FALSE — flags DO NOT EXIST** | **BUG-11 (CRITICAL)** |
| B.6 | compare vs `scec_seisol/`; README format claim | data IS parseable (`.txt`, 8-col) but **README is self-contradictory** ("dir otherwise empty" vs 30 committed files; `.dat`/header vs `.txt`) | BUG-12 |
| B.6b | TPV31+RK45+adjacent is the FIRST matrix+central job | **CONFIRMED** — all tpv104 mixedflux = scalar; all tpv31 = none | — |
| B.7 | deriv-cache `:1146` "safe on scalar/homogeneous TPV operator"; OMIT for first run | **MISQUOTE** — `:1146` actually says "Works for both the scalar and bimaterial (matrix) operators"; OMIT is over-conservative but harmless | BUG-14 |
| B.8 | Phase 7 = benchmark artifacts only; honors Frontera/RK-CFL/no-local rules | **CONFIRMED** — no source edit; rules honored | — |

---

## Bugs Found — Rev 1

### BUG-1: Abort location cited by the plan is the wrong file and wrong lines

**Status:** **RESOLVED in Rev 2**
**File (plan):** Rev-1 Phase 3 step 1 cited `wave_operator.inl:1623-1626` (a comment).
**Severity:** MEDIUM.
**Original description (audit trail):** The actual functional abort is the
`MFEM_VERIFY(m == MixedFluxMode::None, …)` override at
`bimaterial_wave_operator.inl:567-576`, NOT the comment at
`wave_operator.inl:1623-1626`. Three stale comments (`inl:518-520`, `hpp:121`,
`wave_operator.inl:1623-1626`) also needed updating.

**What was fixed (Rev 2):** Phase 3 step 1 now explicitly says "Replace the ACTUAL
abort — the override body at `bimaterial_wave_operator.inl:567-576`, NOT the comment
at `wave_operator.inl:1623-1626`" and enumerates ALL four stale-comment sites
(`inl:518-520`, `hpp:121`, `wave_operator.inl:1623-1626`). **Re-verified against live
source:** the abort is exactly at `bimaterial_wave_operator.inl:569-574`
(`MFEM_VERIFY(m == MixedFluxMode::None, …)`); the `ComputeMaxDt` stale comment is at
`inl:518-520`; the hpp doc comment "Mixed flux is scalar-only (R-003): abort on any
non-None mode" is at `hpp:121`; the base comment is at `wave_operator.inl:1623-1626`.
All four cites are accurate. The plan correctly keeps `Base::SetMixedFluxMode(m)` first
(preserving the Allgather/Allreduce mode cross-check) and cites R-003 for the
structural-exclusion lift. **Resolved.**

---

### BUG-2: Bi-material `ComputeMaxDt` CFL switch does NOT match the scalar switch the plan extracts (factors + missing `cfl_rk_aware_`)

**Status:** **RESOLVED in Rev 2**
**File (code):** `bimaterial_wave_operator.inl:521-531` vs `wave_operator.inl:5880-5903`.
**Severity:** HIGH.
**Original description (audit trail):** The bimaterial body carries its OWN switch with
ADER factors (`Adjacent→0.9`, `AllContinuous→0.4`) and **no `cfl_rk_aware_` gating** —
divergent from the scalar `cfl_rk_aware_?0.6:0.9 / 0.7:0.4`. The Rev-1 plan called the
refactor "byte-exact" without naming the divergent block to delete.

**What was fixed (Rev 2):** Phase 4 step 2 now explicitly says "DELETE the bi-material
operator's OWN switch at `bimaterial_wave_operator.inl:521-531` … and its hardcoded
`cfl_mixed_flux_factor`. Multiply the per-element walk's `dt_e` (`inl:546`) by
`MixedFluxCflFactor_()` instead." It states the change is NOT byte-exact for the
bimaterial operator and explains why (dead arms → intended 0.6/0.7 under RK). Test 4.1
asserts factor equality for all `(mode, cfl_rk_aware_)` pairs.

**Re-verified against live source AND the placement question raised in the review:**
- `inl:521-531` is exactly the divergent switch (factors 1.0/0.9/0.4, no `cfl_rk_aware_`);
  `inl:546` is `const real_t dt_e = cfl_mixed_flux_factor * cfl * h_e / cp_e;` inside the
  per-element walk `inl:539-548`; the `MPI_Allreduce(…, MPI_MIN, …)` is at `inl:553-555`.
- **Per-element `dt_e` placement does NOT change the MPI_Allreduce(MIN) semantics.** The
  factor `f = MixedFluxCflFactor_()` is a single positive constant, identical on every
  element and every rank (it depends only on `mixed_flux_mode_` + `cfl_rk_aware_`, which
  are rank-symmetric). Therefore `min_e (f · x_e) = f · min_e x_e` (uniform positive
  scaling commutes with min), and `MPI_MIN` over ranks likewise commutes:
  `min_ranks (f · local_min) = f · min_ranks(local_min)`. So multiplying `dt_e` per
  element (the plan's choice) is mathematically identical to applying `f` once after the
  walk/Allreduce. The plan's placement is correct. **Resolved.**

---

### BUG-3 (CRITICAL): Plan's central-matrix storage + per-side dispatch is NOT single-valued under heterogeneity

**Status:** **RESOLVED in Rev 2**
**File (plan):** Rev-1 Phase 2 step 2 / Phase 3 dispatch.
**Severity:** CRITICAL.
**Original description (audit trail):** The Rev-1 plan stored a per-SIDE pair
(`std::array<std::array<DenseMatrix,2>,2>`, `[side][0/1]`) built via the upwind
`ResolveFaceFluxOperands_` swap, so side 1 became `c[1][0]=½A_e2, c[1][1]=½A_e1`. With
both sides consuming the unswapped `(Q_e1,Q_e2)`, this gives
`F_h_e1 = ½(A_e1 Q_e1 + A_e2 Q_e2)` but `F_h_e2 = ½(A_e2 Q_e1 + A_e1 Q_e2)` — equal only
when `A_e1==A_e2`. Non-single-valued / non-conservative under heterogeneity, masked by
the homogeneous gate.

**What was fixed (Rev 2):** The storage is redesigned to a **side-symmetric SINGLE pair
per face**:
```cpp
std::unordered_map<int, std::array<mfem::DenseMatrix,2>> per_face_central_flux_;
// [0]=centralE1 (×Q_e1), [1]=centralE2 (×Q_e2).  NO per-side swap.
```
Phase 3 dispatch deposits `F*` once then copies: `ApplyPerFaceFlux(c[0],c[1],Q_self,
Q_nbr,F_h_e1); for(c) F_h_e2[c]=F_h_e1[c];`. The plan explicitly forbids applying the
upwind `ResolveFaceFluxOperands_` side-swap to the central build (Math §; Phase 2 req 2).

**Scrutiny — is this now genuinely single-valued under heterogeneity? YES.**
1. There is exactly ONE stored pair `(½A_e1, ½A_e2)` per face; no side dimension exists,
   so there is no swapped side-1 pair for the swap to leak into.
2. The runtime passes the UNSWAPPED `(Q_self=Q_e1, Q_nbr=Q_e2)` on both sides — confirmed
   against live `InteriorFaceFlux_` (`bimaterial_wave_operator.inl:471-484`), which calls
   `ApplyPerFaceFlux(…, Q_self, Q_nbr, …)` without swapping operands.
3. `F_h_e1 = ½A_e1·Q_e1 + ½A_e2·Q_e2` is computed once; `F_h_e2 ← F_h_e1` by literal copy
   ⇒ `F_h_e1 == F_h_e2` componentwise by construction, for ANY `A_e1, A_e2`. This is
   exactly the scalar reference structure verified at `wave_operator.inl:1203-1209`
   (`flux_.Central(...)` once into `F_h_e1`, then `for(c) F_h_e2[c]=F_h_e1[c]`).
4. **Shared-face interaction is correct:** `SharedInteriorFaceFlux_` fills a single `F_h`
   (side 0 only) and the shared arm stores the side-0 pair `(½A_e1_local, ½A_e2_nbr)`
   keyed by `mesh_face` — consistent with the upwind shared store and the single-`F_h`
   dispatch (verified `inl:489-497`). The side-symmetric store is a strict superset of the
   information the shared single-`F_h` path needs; no swap path remains.
**Resolved.** (The structural guarantee is enforced by Test 3.2 on heterogeneous
material — see BUG-7.)

---

### BUG-4 (CRITICAL): Plan misses the load-bearing residual driver guard that forbids `interior_flux=matrix` + RK

**Status:** **RESOLVED in Rev 2**
**File (code):** `spatial_dyn_driver.cpp:809-816`.
**Severity:** CRITICAL.
**Original description (audit trail):** The guard
`MFEM_VERIFY(!is_rk || interior_flux==Scalar, …)` aborts `matrix + RK` before
`SetMixedFluxMode` (line 1245). Mixed flux requires RK ⇒ feature cannot run.

**What was fixed (Rev 2):** Phase 5 step 1 now targets `spatial_dyn_driver.cpp:809-816`
directly, with a precise rewrite: ALLOW `matrix + RK` when `mixed_flux != none`; REJECT
`matrix + ADER + mixed_flux != none` with the ADER-instability message; for
`matrix + RK + mixed_flux == none` the plan states the chosen behaviour explicitly
(**default: allow**), and directs updating the stale comment at `806-808`.

**Re-verified against live source:** the guard is exactly at `809-816` (the comment at
`806-808` is the "rk* + matrix aborts (guard kept below)" note). The proposed rewrite is
logically sound:
- `matrix + RK + adjacent` → allowed (reaches `SetMixedFluxMode` at 1245). ✓
- `matrix + ADER + adjacent` → rejected, same condition+message as the Phase-4
  `ComputeMaxDt` integrator guard (Phase 4 step 3). The two guards are consistent
  (both: `mixed_flux_mode_ != None && !cfl_rk_aware_` / `time_integrator == ADER`). ✓
- `matrix + RK + none` → allowed by default; harmless (RK on the bimaterial UPWIND
  operator is fine — only central flux needs RK). ✓
- defaults unchanged (`interior_flux=scalar`, `mixed_flux=none`, `time_integrator=ADER`).
**Resolved.** Test 5.1 (parse/dry-run) is the local proxy; the Frontera dispatch gate
(Phase 7 req 4) is the real-mesh end-to-end gate.

---

### BUG-5: `shared_face_neighbour_material_` is keyed by `sf`, not `mesh_face`

**Status:** **RESOLVED in Rev 2**
**File (code):** `bimaterial_wave_operator.hpp:194`; populated `inl:142-163` (keyed `sf`);
consumed in upwind Pass 2 `inl:381-419` (loops `sf`, maps via `pmesh.GetSharedFace(sf)`).
**Severity:** MEDIUM.
**Original description (audit trail):** Rev-1's `ResolveFaceFluxOperands_(mesh_face, …)`
abstraction could not reach the `sf`-keyed neighbour material; no `mesh_face→sf` map.

**What was fixed (Rev 2):** Phase 2 req 3 now mandates the shared arm "mirror the existing
upwind Pass 2 (`inl:375-422`): iterate shared-face index `sf` (NOT `mesh_face`), skip
fault shared faces (as `inl:393-397` does), map `sf → mesh_face` via
`pmesh.GetSharedFace(sf)`, intersect with `central_flux_face_set_.count(mesh_face)`, look
up `shared_face_neighbour_material_[sf]`, build the side-0 pair, store keyed by
`mesh_face`." It keeps the fail-loud `MFEM_VERIFY` at `inl:405-409`.

**Re-verified against live source:** `shared_face_neighbour_material_` is
`std::unordered_map<int, std::array<real_t,3>>` at `hpp:194`, populated keyed by `sf`
at `inl:162`; upwind Pass 2 is at `inl:381-419` with `GetSharedFace(sf)` at `inl:387`,
fault-skip at `inl:393-397`, missing-entry verify at `inl:405-409`. All cites accurate;
the `sf`-loop design is the correct reuse. **Resolved.**

---

### BUG-6: "fail loud" for the missing-neighbour case never triggers; lateral het silently mispopulates

**Status:** **OPEN — NOT correctly resolved in Rev 2.** The Rev-2 mitigation is
*unimplementable as written* — see the escalation in **BUG-10 (NEW, CRITICAL/HIGH)**.
**File (code):** `bimaterial_wave_operator.inl:120-190` (`ExchangeBiMaterialNeighbours_`,
local-side stub); `inl:162` (stores local-as-neighbour); rank-0 WARNING `inl:165-187`.
**Severity:** MEDIUM (carried) → escalated.
**Original description (audit trail):** `ExchangeBiMaterialNeighbours_` ALWAYS populates
`shared_face_neighbour_material_[sf] = per_elem_lmr_[local_elem]` (the LOCAL element's
material), so the entry is never *missing*; a "fail loud if missing" guard never fires
for genuine lateral heterogeneity.

**Rev-2 attempted mitigation (Phase 2 req 4 / Risk Assessment):** "if the neighbour
material in `shared_face_neighbour_material_[sf]` differs from the local material by more
than a tolerance AND the exchange is the local-side stub, abort … naming the partition
seam." This is the WRONG predicate (see BUG-10): it remains **OPEN**.

**Why still open:** The Rev-2 plan now keys the guard on a *mismatch* (`neighbour ≠
local`) instead of *absence* — addressing the Rev-1 wording — but the mismatch predicate
is structurally dead because `inl:162` makes `neighbour == local` by construction. The
plan moved from "guard never fires because entry never missing" to "guard never fires
because neighbour always equals local." Same dead guard, different reason. Detailed in
BUG-10.

---

### BUG-7: Phase 1 single-valuedness acceptance criterion tests the wrong symmetry

**Status:** **RESOLVED in Rev 2**
**File (plan):** Rev-1 Phase 1 AC bullet 2.
**Severity:** MEDIUM.
**Original description (audit trail):** Rev-1's primitive "swap (self,nbr) and states ⇒
same F*" test is trivial commutativity of addition and would pass with the BUG-3 dispatch
defect present; it never builds `F_h_e2`.

**What was fixed (Rev 2):** Phase 1 AC now carries an explicit "Note (BUG-7)": "the
single-valuedness invariant (`F_h_e1==F_h_e2`) is an *operator-level* property and is
tested in Phase 3 (Test 3.2), NOT here. Do NOT add a primitive 'swap (self,nbr) and
states' test." Phase 3 AC adds Test 3.2 (`Dispatch_CentralFlux_SingleValued_Under
Heterogeneity`) on a `mu_e1 != mu_e2` fixture asserting `F_h_e1 == F_h_e2` ≤1e-12·scale,
labelled "the load-bearing guard." **Resolved.** Phase 1 no longer relies on the trivial
test; the het operator-level guard is correctly placed in Phase 3.

---

### BUG-8: Makefile registration understated (4 edits + object-dep list per test)

**Status:** **RESOLVED in Rev 2**
**File (plan/code):** Phase 6 Makefile; `makefile:311-317, 859-893, 1638, 5248, 5252`.
**Severity:** LOW.
**Original description (audit trail):** Rev-1's one-line "register in the make test group"
omitted the SRC/OBJ variables, the link target, the run alias, and the per-test
object-dependency list (primitive vs operator-level differ).

**What was fixed (Rev 2):** Phase 6 req 2 now spells out the **FOUR coordinated edits per
test** (`TEST_<NAME>_SRC`, `TEST_<NAME>_OBJ`, link target, `test-<name>:` alias) PLUS the
aggregate group inclusion, and gives the exact object-dep lists: primitive
(`test_bimaterial_central_flux`) links `GODUNOV_FLUX_BIMATERIAL_OBJ + GODUNOV_FLUX_OBJ`
(cf. `makefile:5248`); operator-level links the full `WAVE_OPERATOR_OBJ +
GODUNOV_FLUX_BIMATERIAL_OBJ + PRECOMPUTED_FACE_FLUXES_OBJ + GODUNOV_FLUX_OBJ +
PML_LAYER_OBJ + FAULT_FACE_FLUX_OBJ + FRICTION_SOLVER_OBJ` (cf. `makefile:1638, 5252`).

**Re-verified against live source:** `makefile:311-317` is the `TEST_*_SRC`/`TEST_*_OBJ`
pattern; `:5248` is `seas_test_godunov_central_flux: $(TEST_GODUNOV_CENTRAL_FLUX_OBJ)
$(GODUNOV_FLUX_OBJ)` (primitive-style link); `:5252` and `:1638` are operator-level links
with the full obj list; `:859-893` is the aggregate group list. All cites accurate.
*Minor note (LOW):* the reference primitive target at `:5248` links only `GODUNOV_FLUX_OBJ`,
whereas the NEW primitive test additionally needs `GODUNOV_FLUX_BIMATERIAL_OBJ` (it calls
`BimaterialFlux::BuildPerFaceCentralMatricesGlobal`) — the plan correctly adds it; the
`:5248` cite is a *pattern* reference, not a copy-the-deps instruction. **Resolved.**

---

### BUG-9: Central builder `mfem::Mult` non-aliasing / sizing + acoustic guard under-specified

**Status:** **RESOLVED in Rev 2**
**File (plan/code):** Phase 1 step 3; `godunov_flux_bimaterial.cpp:282-295`.
**Severity:** LOW.
**Original description (audit trail):** Rev-1 glossed the scratch-matrix discipline,
the 0.5-scaling placement, and the absence of an acoustic guard on the central path
(which bypasses `BuildGodunovStateFaceLocal`).

**What was fixed (Rev 2):** Phase 1 req 3 now says "using **two distinct scratch
`DenseMatrix` objects** (`tmp1`, `tmp2`) and `SetSize` ONLY the output, exactly as
`cpp:282-295` does," states the 0.5 may be folded into `(AxPlus+AxMinus)` or applied via
`centralSide *= 0.5` (to be documented in-code), and forbids aliasing output with scratch.
Phase 1 Edge Cases (BUG-9) adds the explicit `MFEM_VERIFY(flux_self.mu() > eps &&
flux_nbr.mu() > eps, …)` acoustic guard, noting the central path does NOT route through
`BuildGodunovStateFaceLocal`.

**Re-verified against live source:** `cpp:282-295` is exactly the `tmp1/tmp2 + SetSize
output + mfem::Mult` pattern; the central builder's choice of `GetAxPlus()+GetAxMinus()`
matches `GodunovFlux::Central`'s construction at `godunov_flux.cpp:416`
(`Ax_plus_(i,j)+Ax_minus_(i,j)`), which is the correct basis for ≤1e-11 byte-equivalence
with `Central`. **Resolved.**

---

## Bugs Found — Rev 2

### BUG-10 (CRITICAL/HIGH): The Rev-2 BUG-6 lateral-heterogeneity abort can NEVER fire — its predicate is structurally dead

**Status:** OPEN
**File (plan):** Phase 2 req 4; Phase 2 Edge Cases ("Missing/laterally-varying neighbour
material"); Risk Assessment bullet 4; the "Material-heterogeneity support boundary
(BUG-6)" Constraint.
**File (code):** `bimaterial_wave_operator.inl:162` (the load-bearing line);
`inl:139-190` (`ExchangeBiMaterialNeighbours_` stub); rank-0 WARNING `inl:165-187`.
**Severity:** CRITICAL for the *correctness claim* (the plan asserts a safety guard that
does not exist); HIGH for *TPV31 specifically* (TPV31 is depth-only/seam-continuous, so
the dangerous case does not arise for the deliverable — but the plan presents the guard
as the protection and ships a Test 3.4 that cannot pass as written).

**Description:** Rev-2 Phase 2 req 4 prescribes:
> "if the neighbour material in `shared_face_neighbour_material_[sf]` differs from the
> local material by more than a tolerance AND the exchange is the local-side stub, abort."

But `ExchangeBiMaterialNeighbours_` populates the map at `inl:162` with
```cpp
shared_face_neighbour_material_[sf] = per_elem_lmr_[local_elem];   // LOCAL element's material
```
i.e. **the neighbour entry is the LOCAL element's `(λ,μ,ρ)` BY CONSTRUCTION** — there is no
peer-rank value. Therefore, for every shared face, `shared_face_neighbour_material_[sf]`
is **identically equal** to the local material that the proposed guard compares it against.
The predicate `|neighbour − local| > tol` is **always false**, so the abort **never
fires** — including in exactly the dangerous case (genuinely laterally-varying material
across a seam), where the *true* neighbour material differs but the *stored* (stubbed)
value does not. The guard, as specified, is a no-op: it cannot distinguish "seam-continuous
(safe)" from "laterally varying (unsafe)" because both produce a stored neighbour ==
local. The plan moved the guard from keying on *absence* (Rev-1: never absent) to keying
on *mismatch* (Rev-2: never mismatched) — same dead guard, different reason.

**Consequence for Test 3.4 (Phase 3 AC, BUG-5/BUG-6):** Test 3.4 requires the np=2 test to
"abort on a laterally-varying material (BUG-6 guard)." With the proposed predicate, the
operator will NOT abort (it silently builds a *central* matrix from the wrong, stubbed
`A_nbr`), so Test 3.4's abort assertion **cannot pass** against any implementation that
faithfully follows Phase 2 req 4. The test and the spec are mutually inconsistent.

**Why TPV31 still works (and why this is HIGH not blocking for the deliverable):** TPV31's
`depth_profile_1d` material varies only with depth and is continuous across vertical
partition seams; at such a seam the local centroid material genuinely equals the
neighbour's, so the stub is correct and no abort is needed. The deliverable is safe. The
defect is that the plan's *general* safety guarantee is false and the *test that is
supposed to prove the boundary* cannot pass.

**Expected:** A correct detection cannot compare the stubbed neighbour against the local
material (they are equal by construction). It must instead detect lateral variation
*directly*, e.g.:
- (preferred, fully correct) implement the real cross-rank exchange (`MPI_Allgatherv` of
  the peer's `per_elem_lmr_` over shared faces) and THEN compare true-neighbour vs local;
  OR
- (cheaper, detection-only) inspect the **material field's lateral variation near the
  seam directly** — evaluate `material_->EvalAt` at the seam face's two opposite in-plane
  offsets (or at the local element centroid vs the face centroid projected across) and
  abort if the in-plane gradient exceeds tol; OR
- (most conservative, honest) abort *unconditionally* for `Mode::Coefficient + n_shared>0`
  for a *central* shared face UNLESS the material is declared seam-continuous (e.g. a
  config flag `material.seam_continuous=true` that TPV31 sets), making the supported regime
  explicit rather than silently assumed.

In all three, Test 3.4 must assert the abort using a fixture whose *true* lateral variation
is detectable by whatever mechanism is chosen — not via the dead `neighbour vs local`
compare.

**Trigger:** np>1 + `Mode::Coefficient` material varying laterally across a partition
seam + a fault-adjacent face on that seam + `mixed_flux=adjacent`. (Not TPV31.)
**Suggested fix:** Replace Phase 2 req 4's predicate with one of the three detection
mechanisms above; reconcile Test 3.4's abort assertion with the chosen mechanism;
explicitly state in the Constraint that, until the real exchange exists, the supported
regime is "Constant or seam-continuous (depth-only)" and the *enforcement* is the chosen
detector — not a `neighbour vs local` compare. Cross-reference R-004 (the stub limitation)
when changing the WARNING at `inl:165-187` to an abort.

**Rev-3 update — ROOT CAUSE FIXED, REPLACEMENT NOT IMPLEMENTABLE.** Plan Rev 3 removed
the dead `neighbour vs local` predicate everywhere (Constraint "Material-heterogeneity
support boundary", Phase 2 req 4, Phase 2 Edge Cases, Risk Assessment, Tricky-areas) and
replaced it with a *material-field probe*: "evaluate the material field at `x_loc` (local
centroid) and `x_nbr` (a point reflected a small distance across the shared face along the
face normal); abort if `‖m(x_nbr)−m(x_loc)‖/‖m(x_loc)‖ > tol (~1e-6)`; OR require
`material.seam_continuous=true`." The *diagnosis* is now correct — the dead compare is gone.
**But the probe replacement is itself defective on two counts, both NEW CRITICALs:**
(1) the plan assumes `material_->EvalAt(x)` can be called at an ARBITRARY point, which the
real `MaterialField` API does NOT provide (its only point accessor is
`EvalAt(elem, ElementTransformation& T, IntegrationPoint& ip, …)` —
`heterogeneous_material.hpp:183-209`) — see **BUG-15**; and (2) the probe cannot distinguish
legitimate DEPTH variation from unsupported LATERAL variation, so it FALSE-ABORTS on TPV31's
own `depth_profile_1d` across any seam whose normal has a z-component — see **BUG-16**. The
`material.seam_continuous=true` alternative is implementable and suffers neither defect; if
the plan committed to it EXCLUSIVELY (dropping the probe), the detector would be sound. As
written the plan offers BOTH and says "pick ONE," so the probe arm remains a shippable but
broken option and Test 3.4b is told it may exercise either.
**Status (Rev 3): PARTIALLY RESOLVED — the dead-predicate root cause is FIXED; a faithful
implementation of the probe arm is BLOCKED by BUG-15/BUG-16. Resolve by mandating
`seam_continuous` (or the deferred MPI exchange) and deleting the probe arm.**

---

### BUG-11 (CRITICAL): Phase-7 acceptance criterion depends on `visualize_results.py --tol-rms/--tol-peak` flags that DO NOT EXIST

**Status:** OPEN
**File (plan):** Phase 7 req 5; Phase 7 Acceptance Criteria bullet 5 ("`visualize_results.py`
overlay vs `scec_seisol` within the documented ~5–10% peak band").
**File (code):** `tpv31/visualize_results.py` (argparse block lines 373-437);
`tpv31/benchmark_data/README.md:57-63`.
**Severity:** CRITICAL (the *quantitative* validation gate — the whole point of Phase 7's
"SeisSol reference overlay (quantitative validation)" — cannot be executed; the criterion
is unsatisfiable with the committed tooling).

**Description:** Phase 7 req 5 says to compare on-fault traces against the SeisSol
reference "via `tpv31/visualize_results.py` with its `--tol-rms` / `--tol-peak` flags,"
and the acceptance criterion requires the overlay to be "within the documented ~5–10%
peak band." **Grepping the entire script for `tol-rms`, `tol_rms`, `tol-peak`,
`tol_peak`, `--tol`, `rms`, `peak` returns NOTHING.** The script's complete argparse
surface is: positional `mfem` spec(s), `--eqdyna`, `--seisol`, `--both`, `--no-benchmark`,
`--stations`, `--save`, `--output-dir`, `--closeup-t`, plus the directory/prefix args.
It is a **plotting / overlay tool only** — it produces PNGs of MFEM-vs-reference traces;
it computes **no RMS / peak error metric and has no pass/fail tolerance gate**. The
`benchmark_data/README.md:59-63` likewise claims the tolerances are "documented in the
comparison script's `--tol-rms` / `--tol-peak` CLI flags" — describing flags that were
never implemented (the README is aspirational, copied from a plan §R.5 that postdates the
script). Consequently the Phase-7 quantitative acceptance criterion is unsatisfiable as
written.

**Expected:** Either (a) extend `visualize_results.py` (it is a benchmark artifact / tool,
NOT under the source no-touch tree — see BUG-8/scope) to compute per-station RMS and peak
relative error vs the chosen reference and add `--tol-rms`/`--tol-peak` that exit non-zero
on exceedance; OR (b) reword Phase 7 req 5 + the acceptance criterion to a *qualitative*
visual-overlay check ("manually confirm transmission/reflection structure overlays within
~5–10%") and drop the non-existent flags; OR (c) add a separate small comparison script.
The plan must not cite flags that don't exist as a gating criterion.
**Trigger:** Any attempt to run the documented Phase-7 quantitative validation step.
**Suggested fix:** Add a Phase-7 sub-step "extend `visualize_results.py` with `--tol-rms
<x> --tol-peak <y>` computing per-station error and a non-zero exit on exceedance" (it is
the natural home and is in-scope as a benchmark tool), and make the acceptance criterion
reference the real exit-code gate. Note this is the SAME class of omission as Rev-1 BUG-4
(plan asserts a mechanism that the code does not provide).

**Rev-3 update — PARTIALLY RESOLVED.** Plan Rev 3 took option (a): Phase 7 req 5(a) now
makes implementing `--tol-rms/--tol-peak` part of the deliverable (per-station RMS + peak
relative error, non-zero exit on exceedance), adds `visualize_results.py` to a Phase-7
**Files to Modify** block, adds `tpv31/test_visualize_results_tol.py` (Test NEW-7.1), and
reties the acceptance criterion to the exit-code gate. **The structural omission (citing a
nonexistent flag as a gate) is FIXED.** However the *metric specification itself* is
under-defined in three ways that will block a correct, deterministic implementation: the
"relative error" denominator is unspecified and a naive per-sample form divides by reference
values that cross zero (**BUG-17**); the gated station QUANTITY/COLUMN is never named, and
MFEM's column 8 (`mu_eff`, `visualize_results.py:137,154`) has no reference counterpart
(the SeisSol `mu_eff` channel is all-NaN, `:191`) (**BUG-18**); and the metric requires
resampling the two traces (MFEM RK45 adaptive dt vs SeisSol fixed dt) onto a common time
base — which `visualize_results.py` does NOT do today (no `np.interp` anywhere) and the
plan never mentions (**BUG-19**).
**Status (Rev 3): PARTIALLY RESOLVED — flag-existence gap closed; metric definition
incomplete (BUG-17/18/19).**

---

### BUG-12 (MEDIUM): `benchmark_data/README.md` is self-contradictory about the committed SCEC reference (claims dir empty + `.dat` format; reality is 30 `.txt` files)

**Status:** OPEN
**File (code/docs):** `tpv31/benchmark_data/README.md` (lines 32-43 vs 52-55);
`tpv31/benchmark_data/scec_seisol/*.txt` (30 files), `scec_eqdyna/*.txt` (committed);
`tpv31/visualize_results.py:158-192, 201, 575`.
**Severity:** MEDIUM (the *data + tool DO line up*, so the overlay step itself can read the
files — this is a documentation-vs-reality defect that will mislead the implementer about
what to fetch/parse, not a hard blocker for reading the data).

**Description:** The user flagged a possible inconsistency; confirmed. The README states
two things that contradict the committed branch state:
1. README §32-43: "Each station produces a SCEC-format **`.dat`** file" with header
   `# t  slip-rate-1  slip-rate-2  slip-1  slip-2  traction-1  traction-2` ("8-column").
2. README §52-55: "**this directory is otherwise empty (no `.dat` files committed)**."

But the actual committed files are **`.txt`** (not `.dat`), 30 in `scec_seisol/` and 30 in
`scec_eqdyna/`, each 3101 lines, with a multi-line SeisSol comment header (`# problem=
TPV31`, `# author=Thomas ULRICH`, … `# Column #1 = Time (s)` … `t h-slip h-slip-rate
h-shear-stress v-slip v-slip-rate v-shear-stress n-stress`) and 8 numeric columns
(time, h-slip, h-slip-rate, h-shear-stress, v-slip, v-slip-rate, v-shear-stress,
n-stress). **The script CAN consume these:** `visualize_results.py` builds the reference
filename as `tpv31_{code}_{ref_label}.txt` (line 201) and globs `tpv31_{code}_x2_*_x3_*.txt`
(line 575), and `load_reference_file`/`_parse_numeric_table` parse the 8-column table with
`min_cols=8` (lines 158-192), flipping `n-stress` to compression-positive. So the data and
tool DO line up — but the README's "dir is empty" + "`.dat`/different-header" description
is stale and wrong, and Rev-2 Phase 7 req 5 cites the README as the authority for the
format/tolerance. The implementer following the README would believe no reference exists
and would mis-target a `.dat`/7-named-column layout.

**Expected:** Update `benchmark_data/README.md` to describe the ACTUAL committed state:
60 `.txt` files (`scec_seisol/`, `scec_eqdyna/`), the SeisSol multi-line comment header,
the 8 numeric columns (h=strike, v=dip, n-stress compression-negative), and the glob the
script uses. Remove the "directory is otherwise empty" sentence. (This file is a benchmark
artifact, in-scope to edit.)
**Trigger:** Implementer reads the README to acquire/parse reference data.
**Suggested fix:** Rewrite the README "File layout" + "Acquisition gate" sections to match
reality; have Phase 7 req 5 cite the script's actual parser contract, not the stale README.

**Rev-3 update — RESOLVED.** Phase 7 req 5(b) now mandates rewriting
`tpv31/benchmark_data/README.md` to the committed reality (60 `.txt` files in
`scec_seisol/`+`scec_eqdyna/`, the SeisSol multi-line header, the 8 numeric columns with
`h`=strike / `v`=dip / `n-stress` compression-negative, the `tpv31_{code}_x2_*_x3_*.txt`
glob), removing the "directory is otherwise empty" and `.dat`-format sentences, and citing
the script's parser contract (`visualize_results.py:158-192,201,575`) as the format
authority. README added to Files-to-Modify + Tricky-areas. Re-verified the README is still
stale in the live tree (`.dat` `:34`, 7-col header `:38`, "otherwise empty" `:53`, the
phantom `--tol-rms/--tol-peak` reference `:60-61`) so the fix targets the right text.
**Status (Rev 3): RESOLVED.** (NB: the README `:60-61` also references the nonexistent
tolerance flags — the req 5(b) rewrite must drop those too; folded into BUG-11's metric work.)

---

### BUG-13 (MEDIUM): Phase-7 claims `--verify-dispatch` prints the time-integrator and mixed-flux lines; it prints only the interior-flux line (the other two are unconditional banner output)

**Status:** OPEN
**File (plan):** Phase 7 req 3 / req 4 + Acceptance Criteria ("`--verify-dispatch` prints
the resolved scheme … prints `interior flux: matrix`, `mixed flux: adjacent`, `time
integrator: rk45`").
**File (code):** `spatial_dyn_driver.cpp:2130-2147` (`--verify-dispatch` block);
`:832-866` (unconditional rank-0 startup banner: `time integrator:` line 853, `mixed
flux:` line 859); `:2143-2146` (verify-dispatch interior-flux line).
**Severity:** MEDIUM (the grep gate is still *satisfiable* — the three strings DO all
appear in the `*.out` — but the plan mis-attributes WHERE, the exact string forms differ
from the plan's quoting, and the criterion's premise ("`--verify-dispatch` prints them")
is false for 2 of 3).

**Description:** The Phase-7 dispatch gate (req 4) and req 3's comment assert that
`--verify-dispatch` emits `interior flux:`, `mixed flux:`, and `time integrator:`. In the
live driver:
- `--verify-dispatch` (gated block, `:2130-2147`) prints ONLY
  `"[verify-dispatch] interior flux : "` followed by `"matrix (heterogeneous bimaterial
  Riemann)"` (note: substring is `interior flux :` with a space before the colon, and the
  value is NOT the bare `matrix`). It does NOT print mixed-flux or time-integrator.
- `"time integrator:  "` (line 853) and `"mixed flux:       "` (line 859) are emitted by
  the **unconditional rank-0 startup banner** (`if (rank == 0)` at `:832`), regardless of
  `--verify-dispatch`, with multi-space padding (`mixed flux:       adjacent`).

So all three strings DO land in the `*.out` (banner + verify-dispatch), and a tolerant
grep (`grep -E "interior flux|mixed flux:|time integrator:"`) will find them — this is
exactly what the existing TPV104 mixedflux sbatch does (it greps `*.out` for "time
integrator:", "mixed flux:", "law:" — the BANNER strings — per its own comment at
job line 37). But the plan's statement that `--verify-dispatch` is what prints all three
is wrong, and the exact match strings in the plan (`interior flux: matrix`,
`mixed flux: adjacent`) do not match the emitted text byte-for-byte (`interior flux :
matrix (heterogeneous bimaterial Riemann)`; `mixed flux:       adjacent`). A brittle
exact-string grep in the sbatch pre-flight could fail.

**Expected:** Reword Phase 7 req 3/4 + the acceptance criterion to state that the resolved
scheme is verified by grepping the `*.out` for the BANNER lines (`time integrator:`,
`mixed flux:`) plus the `--verify-dispatch` `interior flux :` line, using tolerant
(substring / regex, whitespace-insensitive) patterns — mirroring the proven TPV104/TPV31
sbatch grep style — rather than asserting `--verify-dispatch` prints all three or quoting
exact-spacing strings.
**Trigger:** Implementing the sbatch pre-flight/post-run grep with the plan's exact strings.
**Suggested fix:** Use the existing TPV31/TPV104 sbatch banner-grep pattern verbatim;
correct the plan's attribution and quoting.

**Rev-3 update — RESOLVED.** Phase 7 req 3 header softened and req 4 rewritten to state the
three scheme lines come from TWO sources — the unconditional rank-0 banner (`time integrator:`
`spatial_dyn_driver.cpp:853`, `mixed flux:` `:859`, with multi-space padding) and the
`--verify-dispatch` block (`interior flux : matrix (heterogeneous bimaterial Riemann)`,
`:2130-2147`) — and to grep the `*.out` with whitespace-tolerant regex for
`time integrator:`, `mixed flux:`, `interior flux`, mirroring the proven TPV104/TPV31
sbatch banner-grep style rather than brittle exact strings. The acceptance criterion was
updated to match. **Status (Rev 3): RESOLVED.**

---

### BUG-14 (LOW): Phase-7 misquotes the `--deriv-cache` driver comment; OMIT is over-conservative (the code does NOT reject `--deriv-cache` on the bimaterial operator)

**Status:** OPEN
**File (plan):** Phase 7 req 3 ("`--deriv-cache` … the TPV104 job notes it is gated 'safe
on this scalar/homogeneous TPV operator', `spatial_dyn_driver.cpp:1146`" → OMIT for first
run unless confirmed safe on `BimaterialWaveOperator`).
**File (code):** `spatial_dyn_driver.cpp:1140-1158`.
**Severity:** LOW (OMITting `--deriv-cache` is a safe, performance-only choice for a first
run, so the *recommendation* is harmless — but the *justification* misquotes the code and
implies a restriction that does not exist).

**Description:** Rev-2 Phase 7 req 3 says the deriv-cache "is gated 'safe on this
scalar/homogeneous TPV operator' (`spatial_dyn_driver.cpp:1146`)" and tells the
implementer to verify "the cache path supports the bimaterial operator before enabling
it." The live comment at `:1140-1148` says the opposite of "scalar-only":
> "Works for both the scalar and bimaterial (matrix) operators (the cache is
> geometry-only). … R-004: SetDerivMode aborts fail-loud if the cache exceeds the per-rank
> budget (… p2 ~18 MB/rank, p3 ~72 MB/rank for TPV31)."
There is **no gate that rejects `--deriv-cache` on `BimaterialWaveOperator`** — `:1149-1158`
unconditionally calls `wave.SetDerivMode(DerivMode::Cached)` when the flag is present
(`wave` is the base reference; `SetDerivMode` is the geometry-only cache, byte-identical
up to round-off per R-002). The TPV31 sizing (~18/72 MB/rank) is even called out for TPV31
specifically. So the code does NOT "outright reject it"; OMITting it loses a ~4-5x ADER
step speedup for no correctness reason (and TPV31 uses RK45 here, not ADER, so the cache's
ADER hot-path benefit may not even apply — another reason the OMIT is moot).

**Expected:** Correct Phase 7 req 3 to quote the actual comment ("works for both operators;
geometry-only; round-off-only per R-002; budget-guarded per R-004"). OMITting for the
*first* run is still a reasonable conservative default (one less variable while calibrating
CFL), but the *reason* should be "isolate the CFL calibration / RK45 doesn't hit the ADER
hot path," NOT "possibly unsafe on the bimaterial operator," which the code contradicts.
**Trigger:** Implementer reads req 3 and believes `--deriv-cache` is gated/unsafe on the
matrix operator.
**Suggested fix:** Replace the misquote with the real comment text; keep OMIT as an
explicit, correctly-justified first-run choice.

**Rev-3 update — RESOLVED.** Phase 7 req 3 now quotes the real driver comment
(`spatial_dyn_driver.cpp:1140-1148`: the deriv-cache "works for both the scalar and
bimaterial (matrix) operators … geometry-only", round-off-only per R-002, budget-guarded per
R-004 with the TPV31 `~18 MB/rank` p2 / `~72 MB/rank` p3 sizing). OMIT is retained for the
FIRST run but correctly re-justified — "isolate the CFL calibration with one fewer variable"
and "RK45 doesn't hit the ADER hot path the cache accelerates" — explicitly NOT "possibly
unsafe." **Status (Rev 3): RESOLVED.**

---

## Bugs Found — Rev 3

### BUG-15 (CRITICAL): The BUG-10 material-field probe is unimplementable against the real `MaterialField` API — there is no arbitrary-point evaluator

**Status:** OPEN
**File (plan):** Phase 2 req 4 (the "material-field probe" arm); the "Material-heterogeneity
support boundary" Constraint; Phase 2 Edge Cases; Risk Assessment bullet 4; Test 3.4b.
**File (code):** `dynamic/heterogeneous_material.hpp:183-209` (`MaterialField::EvalAt`);
`dynamic/bimaterial_wave_operator.hpp:182` (`const MaterialField *material_`);
`dynamic/bimaterial_wave_operator.inl:64-70` (the ONLY material-eval pattern in the operator);
`fem/coefficient.cpp:173-179` (`FunctionCoefficient::Eval` ⇒ `T.Transform(ip, transip)`).
**Severity:** CRITICAL for implementability — the prescribed detector calls an API that does
not exist; a faithful implementer cannot write it. (Not blocking for TPV31 *if* the
`seam_continuous` arm is mandated instead — see Suggested fix.)

**Description:** Phase 2 req 4 says, for `Mode::Coefficient`, "the material is an analytic
function of position (`material_->EvalAt(x)` / the existing per-coordinate evaluator),
evaluable at ANY point on the local rank with no communication. Probe the local element
centroid `x_loc` and a point reflected a small distance across the shared face along the face
normal into the neighbour side, `x_nbr`. Compute `(λ,μ,ρ)` at both." This assumes a
point-evaluable `EvalAt(x)`. **The real `MaterialField` has no such method.** Its only
position accessor is

```cpp
// heterogeneous_material.hpp:183
inline void EvalAt(int /*elem*/, mfem::ElementTransformation& T,
                   const mfem::IntegrationPoint& ip,
                   real_t& lambda_out, real_t& mu_out, real_t& rho_out) const;
```

which in `Mode::Coefficient` calls `lambda_coef->Eval(T, ip)` (`hpp:201-203`).
`FunctionCoefficient::Eval(T, ip)` resolves the physical point via `T.Transform(ip, transip)`
(`fem/coefficient.cpp:173-179`) — i.e. the physical coordinate is whatever the
`ElementTransformation T` maps the reference `IntegrationPoint ip` to. The only material-eval
pattern in the operator (`bimaterial_wave_operator.inl:64-70`) is exactly this: get the
LOCAL element's `T = mesh_.GetElementTransformation(e)`, take the centroid `ip`, call
`material.EvalAt(e, *T, ip, …)`.

Consequences:
1. **`x_loc` is fine** — it is the local element centroid, reachable with the existing
   `(local T, centroid ip)` pattern.
2. **`x_nbr` is NOT reachable.** To evaluate the coefficient at a physical point `x_nbr`
   that lies OUTSIDE the local element (across the seam in the neighbour element), one must
   supply an `(ElementTransformation, IntegrationPoint)` pair whose `Transform` lands at
   `x_nbr`. For a SHARED (cross-rank) face the neighbour element is on another rank and the
   operator holds no transform for it (`FaceElementTransformations` for a shared face gives
   only the local `Elem1` side; `Elem2` is a face-neighbour proxy, not a volume transform you
   can `Transform` an arbitrary `ip` through to an interior point). For a fully-local face one
   could in principle build a throwaway `IsoparametricTransformation`, but the plan gives no
   recipe and the existing code has no such helper.
3. The coordinate-only evaluator the plan alludes to ("the existing per-coordinate
   evaluator") is `DepthProfile1DMaterial::eval_at_xyz` (`heterogeneous_material.cpp:318-334`),
   a `std::function<void(x,y,z,…)>`. **It lives ONLY on the `DepthProfile1DMaterial` WRAPPER,
   not on `MaterialField`.** The operator stores `const MaterialField *material_`
   (`hpp:182`), assigned `&material` in the ctor (`inl:46`) — it has NO handle to the wrapper
   and NO way to recover `eval_at_xyz`. So even the "per-coordinate evaluator" the plan names
   is unreachable from inside `BimaterialWaveOperator`.

In short: the operator CAN re-derive `m(x_loc)` (it already does for every element in the
ctor), but it has no faithful, communication-free way to evaluate `m(x_nbr)` at a point in a
neighbour element across a partition seam. The probe is unimplementable as written.

**Expected:** A detector that uses only APIs the operator actually has. Two sound options:
- **(preferred) Drop the probe; mandate `material.seam_continuous=true`.** Abort any
  `Mode::Coefficient` *central shared* face when the flag is absent. This needs no
  arbitrary-point eval, no cross-rank comm, is deterministic, and TPV31 sets it. It also
  sidesteps BUG-16 entirely. (The plan already lists this as the "equivalent alternative" —
  make it the ONLY mechanism and delete the probe arm.)
- **(deferred, fully correct) the real `MPI_Allgatherv`** of the peer's `per_elem_lmr_` over
  shared faces, then compare true-neighbour vs local — already the plan's deferred follow-up.

If the plan insists on a field probe, it must (a) name a concrete API to evaluate the
coefficient at an arbitrary `x` (e.g. expose `eval_at_xyz` on `MaterialField`, or pass the
`DepthProfile1DMaterial` wrapper to the operator, or build an `IsoparametricTransformation`
for the local element and a documented reference-point lookup), and (b) restrict the probe to
the LOCAL face geometry it can actually reach (see BUG-16, which shows even a reachable probe
gives the wrong answer for depth-varying material).

**Trigger:** Implementer reaches Phase 2 req 4 and tries to write `m(x_nbr)`; there is no
`MaterialField::EvalAt(x)`; the build does not compile / the implementer must invent an API
the plan never specified.
**Suggested fix:** Replace Phase 2 req 4's probe with the `seam_continuous` flag as the sole
detector (delete the probe arm from the Constraint, Phase 2 req 4 + Edge Cases, Risk
Assessment, Tricky-areas, and Test 3.4b). Keep the `MPI_Allgatherv` exchange as the deferred
fully-correct path. This also resolves BUG-16 and BUG-20.
**Test:** Test 3.4b (see Proposed Tests) must be (re)written against the `seam_continuous`
mechanism — a `Mode::Coefficient` central shared face WITHOUT `seam_continuous=true` aborts;
WITH it, builds. A probe-based 3.4b cannot be written against the real API.

**Rev-4 update — RESOLVED.** Plan Rev 4 DELETED the material-field probe everywhere and
mandates the declarative `material.seam_continuous` flag as the SOLE guard (Constraint
`plan:88-99`, Phase 2 req 4 `:272-296`, Edge Cases `:309-313`, Risk `:814-822`, Tricky-areas
`:846-851`). Every remaining "probe" mention is the REJECTED option (grep-confirmed: lines
7,9,89,91,93,273,276,281-282,311,409,414,748,819,849-850 all say "NOT a probe" / "probe is
unimplementable" / "an along-normal probe would wrongly abort"). The guard now reads only the
`seam_continuous_` member — no `MaterialField::EvalAt(x)`, no `eval_at_xyz`, no cross-rank
comm (plan `:289`: "needs NO arbitrary-point material eval"). Test 3.4b is rewritten against
the flag (`:406-410`). **The unimplementable-API defect is gone. RESOLVED.** (NB: the *config
parse* of the new flag, and the field on `MaterialSpec`, land in `spatial/code/` — out of the
plan's declared scope; tracked separately as BUG-22, and the flag's TOML key co-trips an
unrelated parse guard, BUG-21.)

---

### BUG-16 (CRITICAL): The material-field probe cannot separate legitimate DEPTH variation from unsupported LATERAL variation — it FALSE-ABORTS on TPV31's own `depth_profile_1d`

**Status:** OPEN
**File (plan):** Phase 2 req 4 (the probe + its claim "For TPV31's depth-only profile, `x_loc`
and `x_nbr` sit at essentially the same depth across the seam, so `m(x_nbr) ≈ m(x_loc)` and the
abort does NOT fire"); Phase 7 Edge Cases ("the Phase-2 lateral-heterogeneity abort does not
fire for TPV31").
**File (code):** `miniapps/seas/tpv31/configs/tpv31.toml:97-145` (depth_profile_1d:
`depth_axis="z"`, step discontinuities at depth 2400/5000/10000 m, linear interp in
[2400,5000] m); `dynamic/heterogeneous_material.cpp:240-295` (`eval_at_depth` +
`axis_value_to_depth`: depth = max(0,-z)); fault is the vertical `y=0` plane
(`tpv31.toml:52,71-76`).
**Severity:** CRITICAL — the plan's central claim ("the abort does NOT fire for TPV31") is
FALSE for general tet seams, so the *prescribed* detector would block the deliverable; or,
if its tol is loosened to avoid that, it can no longer detect the lateral case it exists for.

**Description:** The detector probes along the FACE NORMAL: `x_nbr = x_loc + δ·n̂`. TPV31's
material is a function of DEPTH only (`m = m(depth)`, depth = max(0,-z)), with a LINEAR ramp
in depth ∈ [2400,5000] m and STEP JUMPS at depth 2400/5000/10000 m (`tpv31.toml:100-145`).
Now consider the partition seams the probe runs on:
- A partition seam is an arbitrary interior face produced by ParMETIS; in a general tet mesh
  its unit normal `n̂` has a nonzero vertical component `n_z` for most faces (only faces lying
  in a horizontal plane have `n_z=0`). The plan's central-set is the *fault-adjacent* faces;
  the fault is the vertical `y=0` plane, but the fault-ADJACENT faces (the other faces of the
  fault-touching tets) are generic tet faces with arbitrary orientation, many with `n_z≠0`.
- For such a face, `x_loc` and `x_nbr = x_loc + δ·n̂` are at DIFFERENT depths
  (`Δdepth = -δ·n_z ≠ 0`). Therefore `m(x_nbr) ≠ m(x_loc)` **by the legitimate depth profile**:
  in the linear-ramp band the relative change is `≈ |∂_depth m|·δ·|n_z| / |m|`, and for a face
  STRADDLING one of the depth discontinuities (2400/5000/10000 m) it is a FULL material JUMP
  (e.g. vs jumps 2250→2550 m/s at 2400 m ⇒ `μ ∝ ρ v_s²` jumps ~28%).
- With the plan's tol `~1e-6` relative, ANY such face (and there are many — every
  fault-adjacent tet face with `n_z≠0` that happens to be a partition seam) trips
  `‖m(x_nbr)−m(x_loc)‖/‖m(x_loc)‖ > tol` and **ABORTS** — even though TPV31 is a fully
  SUPPORTED depth-only material. The plan's assertion that "`x_loc` and `x_nbr` sit at
  essentially the same depth across the seam … the abort does NOT fire" is only true for a
  seam whose normal is horizontal (`n_z=0`); it is FALSE for the general tet seams TPV31's
  ParMETIS partition produces.

So the probe, as specified, is a **false-positive blocker for the deliverable itself**. The
deeper problem: the probe measures the TOTAL material change across the seam, but "seam
continuity" for this stub means the LOCAL-element material equals the TRUE NEIGHBOUR material
*at the seam* — which for a depth-only profile is exactly true (both sides evaluate the same
`m(depth)` at the same physical seam point), regardless of how fast `m` varies with depth. A
correct detector must compare against the depth-profile PREDICTION (or, equivalently, probe
only the IN-PLANE/horizontal offset, holding depth fixed), not the along-normal total change.
The plan's "depth_profile_1d is seam-continuous" statement is accurate for the STUB's needs
(local == true-neighbour at the seam) but is exactly what the along-normal probe FAILS to
verify.

**Expected:** Do not probe along the face normal with a relative tol that conflates depth and
lateral variation. Correct options:
- Adopt the `seam_continuous` flag (BUG-15 preferred fix) — declarative, no false positives.
- If a probe is kept, restrict it to the IN-PLANE direction: pick two points on the seam face
  at the SAME depth but offset horizontally (in the face plane, perpendicular to the depth
  axis) and abort only if THOSE differ — this isolates lateral variation. But this still needs
  an arbitrary-point evaluator (BUG-15) and careful handling when the face plane is itself
  vertical.
**Trigger:** TPV31 (or any depth-varying material) run in parallel with `mixed_flux=adjacent`,
where a fault-adjacent tet face with `n_z≠0` falls on a partition seam — the common case.
**Suggested fix:** Same as BUG-15 — mandate `seam_continuous` and delete the along-normal
probe. Update Phase 2 req 4 and the Phase-7 Edge-Case claim ("the abort does not fire for
TPV31") to reflect the declarative mechanism.
**Test:** Add Test 3.4c (`Dispatch_DepthProfile_NoFalseAbort`): a depth-varying material with
`seam_continuous=true` on an np=2 fixture whose seam normal has a z-component must BUILD (not
abort); the SAME material WITHOUT the flag aborts. (A probe-based detector fails the first
half.)

**Rev-4 update — RESOLVED.** With the along-normal probe deleted (BUG-15 Rev-4 update), the
depth-vs-lateral conflation that false-aborted TPV31 cannot occur: the declarative flag does
not inspect material at all, so legitimate `depth_profile_1d` variation across an `n_z≠0` seam
never trips it. Plan Phase 7 Edge Cases (`:744-749`) now correctly states "There is NO runtime
material probe (it would false-abort on the depth discontinuities at 2400/5000/10000 m)" and
that `seam_continuous=true` satisfies the guard. Test 3.4c added (`:411-414`) as the explicit
no-false-abort reproducer. **RESOLVED.**

---

### BUG-17 (MEDIUM): Phase-7 "relative error" denominator is unspecified — a naive per-sample relative error divides by reference values that cross zero

**Status:** OPEN
**File (plan):** Phase 7 req 5(a) ("per-station RMS and peak **relative** error of the MFEM
trace vs the loaded reference"); Phase 7 acceptance ("within the documented ~5–10% peak band").
**File (code):** `tpv31/visualize_results.py:158-192` (reference reader); the on-fault SCEC
quantities (slip, slip-rate, shear stress) all cross or sit at zero before nucleation arrives.
**Severity:** MEDIUM — the gate is implementable, but "relative error" without a named
denominator is ambiguous and the obvious per-sample form is numerically ill-posed.

**Description:** Phase 7 req 5(a) requires "per-station RMS and peak **relative** error" but
never says relative to WHAT. Candidates: (i) per-sample reference value `|mfem−ref|/|ref|`;
(ii) per-station peak `|mfem−ref|/max|ref|`; (iii) RMS of the reference
`rms(mfem−ref)/rms(ref)`. Option (i) is numerically catastrophic: TPV31 on-fault traces
(slip, slip-rate, the shear-stress perturbation) are ZERO before the rupture front arrives and
slip-rate returns toward zero after — `|ref|→0` makes per-sample relative error blow up
(divide-by-zero / huge spikes at every zero-crossing), so the "peak relative error" would be
dominated by noise near zeros, not by physical disagreement. The "~5–10% peak band" language
only makes sense for a peak-normalized metric (ii).

**Expected:** Name a well-defined, zero-safe metric. Concretely:
`peak_rel = max_t|mfem(t)−ref(t)| / max_t|ref(t)|` and
`rms_rel  = sqrt(mean_t (mfem(t)−ref(t))²) / sqrt(mean_t ref(t)²)` (or
`/ max_t|ref(t)|`), with a small floor on the denominator. These are the standard SCEC-style
peak/RMS-normalized errors and are stable through zero-crossings.
**Trigger:** Implementing req 5(a) literally with a per-sample `|mfem−ref|/|ref|`.
**Suggested fix:** Specify the peak-normalized (and RMS-normalized) metric in req 5(a); state
the denominator floor; cite it in Test NEW-7.1's pass/fail thresholds.

**Rev-4 update — RESOLVED.** Plan Phase 7 req 5(a) now specifies the zero-safe metric
verbatim (`:718-722`): `peak_rel = max_t|mfem(t)-ref(t)| / max(max_t|ref(t)|, floor)` and
`rms_rel = sqrt(mean_t (mfem-ref)^2) / max(sqrt(mean_t ref^2), floor)`, with an explicit
"small denominator `floor`" and an explicit warning that this is "NOT a per-sample
`|mfem-ref|/|ref|`, which blows up at the pre-nucleation zero-crossings." Test NEW-7.4
(`TolGate_PeakNormalized_ZeroCrossingSafe`) added as the zero-crossing reproducer (plan
`:621-625, 726`). **RESOLVED.**

---

### BUG-18 (MEDIUM): Phase-7 tolerance gate does not name the gated station QUANTITY/COLUMN; MFEM's column 8 (`mu_eff`) has no reference counterpart

**Status:** OPEN
**File (plan):** Phase 7 req 5(a) ("for each compared on-fault station, compute the per-station
RMS and peak relative error of the MFEM trace vs the loaded reference").
**File (code):** `tpv31/visualize_results.py:125-155` (MFEM 9-column reader; col 8 = `mu_eff`),
`:158-192` (reference 8-column reader; `mu_eff` set all-NaN at `:191`), `:205-214` (the 8
PANELS keys).
**Severity:** MEDIUM — without a named quantity the gate is ambiguous and could either compare
the wrong channel or attempt to compare a channel with no reference (all-NaN ⇒ NaN error).

**Description:** A "station trace" here is multi-channel: MFEM writes 9 columns
(t, h-slip, h-slip-rate, h-shear, v-slip, v-slip-rate, v-shear, n-stress, mu_eff —
`visualize_results.py:128-137`); the SCEC reference has 8 columns (no `mu_eff`), and the
reader fills `mu_eff` with `np.full(n, np.nan)` (`:191`). The plan says "the MFEM trace vs the
loaded reference" but does not say WHICH of the 7 shared physical channels (h/v slip, slip-rate,
shear stress, plus n-stress) the `--tol-rms/--tol-peak` thresholds apply to. If the gate naively
iterates all 8 PANELS keys it will (a) include `mu_eff`, whose reference is all-NaN ⇒ the error
is NaN ⇒ the `>tol` comparison is `False` (NaN compares false), silently passing; and (b) apply
one scalar tol across channels with wildly different magnitudes (slip ~meters, slip-rate ~m/s,
stress ~MPa) — a single relative-tol per channel is fine, but the plan should say so.

**Expected:** Name the gated quantity/quantities explicitly. The natural primary gate is the
slip-rate channels (`V_strike`, and `V_dip` if non-trivial) and/or peak slip and shear stress —
the standard SCEC on-fault comparison quantities. Exclude `mu_eff` (no reference). State that
the relative error is computed per-channel (peak-normalized per BUG-17) and the station fails if
ANY gated channel exceeds either tol.
**Trigger:** Implementing req 5(a) without choosing a column ⇒ either a NaN-masked false pass or
an apples-to-oranges cross-channel tol.
**Suggested fix:** In req 5(a), name the gated channels (e.g. `V_strike`, `slip_strike`,
`tau_strike`; `V_dip/slip_dip/tau_dip` only if the station has dip motion; `sigma_n`), exclude
`mu_eff`, and require a per-channel peak-normalized error. Reference the
`visualize_results.py:145-155` column map by name.

**Rev-4 update — RESOLVED.** Plan Phase 7 req 5(a) now names the gated channels exactly
(`:714-717`): `V_strike`, `slip_strike`, `tau_strike` (and `_dip` only if the station has dip
motion), plus `sigma_n`, and **excludes `mu_eff` (MFEM col 8)** with the correct rationale
"its reference is all-NaN (reader `:191`), so a NaN error compares false and would *silently
pass*," citing the column map at `visualize_results.py:145-155`. **VERIFIED against live
script:** the named keys are exactly the dict keys in `load_mfem_file` (`:145-155`) and
`load_reference_file` (`:179-192`) and the `PANELS` list (`:205-214`); `mu_eff` is MFEM col 8
(`:137,154`) and the reference fills it `np.full(n, np.nan)` (`:191`). The plan's names are
correct and unambiguous. **RESOLVED.**

---

### BUG-19 (HIGH): Phase-7 RMS/peak comparison requires resampling onto a common time base; MFEM (RK45 adaptive dt) and SeisSol traces are on different samples, and the script does not interpolate

**Status:** OPEN
**File (plan):** Phase 7 req 5(a) (the RMS/peak metric); Testing Strategy / acceptance gate.
**File (code):** `tpv31/visualize_results.py` (NO `np.interp` / resampling anywhere — confirmed
by grep; the plot path at `:235-260` plots each dataset on its OWN `time_s` array);
`:146,180` (MFEM `time_s` vs reference `time_s` are independent column-0 arrays).
**Severity:** HIGH — without a common time base the RMS/peak difference `mfem(t)−ref(t)` is
undefined (the two arrays have different lengths and different sample times); a naive
element-wise subtraction either crashes (length mismatch) or silently compares mismatched
samples, making the gate wrong.

**Description:** The MFEM driver uses adaptive RK45 (`--time-integrator rk45`, no `--dt`), so
its station output is on irregular, run-dependent time samples; the SeisSol reference `.txt`
files are on a fixed SeisSol output cadence (each 3101 lines, regular dt). Computing
`mfem(t)−ref(t)` requires interpolating BOTH onto a common time grid (e.g. interpolate the MFEM
trace onto the reference's `time_s`, restricted to the overlapping `[t0,t1]`). The plotting tool
does NOT do this — it draws two independent `(t, y)` curves — and there is no `np.interp` (or
any resampling) anywhere in `visualize_results.py`. The plan's req 5(a) and Test NEW-7.1 never
mention interpolation, so a literal implementation of "compute RMS/peak error of the MFEM trace
vs the reference" has no defined operation on two differently-sampled arrays.

**Expected:** Specify the resampling step: interpolate the MFEM channel onto the reference
`time_s` (or onto a shared uniform grid) over the overlapping time window via `np.interp`,
THEN compute the per-channel peak-normalized RMS/peak error (BUG-17). Test NEW-7.1's synthetic
fixtures must use DIFFERENT time samples for the MFEM `.dat` vs the reference `.txt` to exercise
the interpolation (otherwise the test passes trivially on aligned samples and the production gate
silently breaks on misaligned ones).
**Trigger:** Running the real Phase-7 gate, where MFEM and SeisSol samples never coincide.
**Suggested fix:** Add an interpolation step to req 5(a) and to Test NEW-7.1's spec; require the
test fixtures to be on offset time grids.

**Rev-4 update — RESOLVED (with a residual gap → BUG-23).** Plan Phase 7 req 5(a) now
mandates the resampling step (`:708-712`): "interpolate the MFEM channel onto the reference
`time_s` over the overlapping `[t0,t1]` window via `np.interp` BEFORE any difference." Test
NEW-7.1 is required to use OFFSET time grids (`:622`, `:725`) so the interpolation is exercised.
**The undefined-difference defect is FIXED.** RESOLVED. *Residual:* the plan names the
`[t0,t1]` window but never DEFINES it (presumably `[max(t_starts), min(t_ends)]`) and does not
handle a truncated MFEM run (partial coverage) — see the NEW finding BUG-23.

---

### BUG-20 (LOW): The probe's reflection geometry and tolerance are under-defined for a general (non-axis-aligned) tet seam face

**Status:** OPEN (moot if BUG-15/BUG-16 are fixed by dropping the probe)
**File (plan):** Phase 2 req 4 ("a point reflected a small distance across the shared face along
the face normal into the neighbour side, `x_nbr`"; "small relative tol, e.g. 1e-6").
**File (code):** `dynamic/bimaterial_wave_operator.inl:330-342` (centroid unit-normal build for
faces); general tet faces are non-axis-aligned.
**Severity:** LOW — a specification gap, secondary to BUG-15/BUG-16.

**Description:** "Reflected a small distance across the face along the face normal" leaves three
things undefined: (1) the distance `δ` — relative to WHAT (element size `h_e`? an absolute
meters value? a fraction of the inscribed diameter)? Too small and round-off dominates; too
large and `x_nbr` overshoots the neighbour element into a third element or out of the domain.
(2) Whether `x_nbr` is GUARANTEED to land inside the neighbour element — for a skewed tet or a
boundary-adjacent seam, `x_loc + δ·n̂` may land in a void or a different element, giving a
material value that is neither the local nor the true neighbour. (3) The `1e-6` relative tol is
far below the round-off floor of a coefficient evaluated through `T.Transform` on a curved/skewed
frame, and (per BUG-16) is far below the legitimate depth-variation signal — so it is
simultaneously too tight for round-off and meaningless for the physics.
**Expected:** If a probe is retained at all, define `δ` as a fraction of `per_elem_h_[e]`
(e.g. `0.25·h_e`), document that `x_nbr` must be verified to lie in the intended neighbour
element (or evaluate AT the seam face centroid from both sides instead of reflecting), and set a
physically-motivated tol (not `1e-6`). **Better:** drop the probe (BUG-15/BUG-16) and this
finding disappears.
**Trigger:** Implementing the probe on a real tet mesh.
**Suggested fix:** Subsumed by the BUG-15 fix (replace the probe with `seam_continuous`).

**Rev-4 update — RESOLVED (dissolved).** The probe is deleted (BUG-15 Rev-4 update), so the
reflection distance `delta`, neighbour-element containment, and `1e-6` tol are no longer part
of any mandated mechanism. The declarative flag has no geometry. **RESOLVED.**

---

## Bugs Found — Rev 4

### BUG-21 (CRITICAL): The plan relaxes only ONE of TWO series guards — a config-PARSE mutual-exclusion guard at `spatial_friction.cpp:1162-1167` aborts `matrix + mixed_flux!="none"` BEFORE the driver guard the plan fixes

**Status:** **RESOLVED in Rev 5**
**File (plan):** Phase 5 (BUG-4 relaxes ONLY `spatial_dyn_driver.cpp:809-816`); Phase 7
config (sets `[numerics].mixed_flux = "adjacent"` in `tpv31_rk_mixedflux.toml`, plan
`:606-608`, `:754-755`); Constraint "Source-edit scope" (`:71-75`).
**File (code):** `spatial/code/spatial_friction.cpp:1162-1167` (the parse-level
mutual-exclusion `MFEM_VERIFY`); `:1152-1157` (the `mixed_flux` value validator);
`drivers/spatial_dyn_driver.cpp:616` (`cfg = LoadSpatialFrictionConfig(config_path)`),
`:679` (CLI `--mixed-flux` override, applied AFTER parse), `:809-816` (the driver guard
BUG-4 relaxes), `:1135` (operator ctor), `:1245` (`SetMixedFluxMode`).
**Severity:** CRITICAL — the Phase-7 deliverable config cannot even be loaded; it aborts at
`LoadSpatialFrictionConfig` before any code the plan touches runs. This is the EXACT same
class of omission as Rev-1 BUG-4 (the plan misses a load-bearing guard), now in the *config
parser* rather than the driver — and it is the THIRD round in a row where the author's fix
for one guard exposed an unaddressed sibling.

**Description:** `LoadSpatialFrictionConfig` validates the `[numerics]` block and aborts:
```cpp
// spatial_friction.cpp:1162-1167  (Phase 6 req 3 mutual-exclusion, R-1203 sibling)
MFEM_VERIFY(cfg.numerics.interior_flux == InteriorFlux::Scalar
            || cfg.numerics.mixed_flux == "none",
            "[numerics] interior_flux=\"matrix\" is incompatible with "
            "mixed_flux=\"" << cfg.numerics.mixed_flux << "\" (mixed-flux "
            "is valid only on the scalar interior-flux path); set "
            "mixed_flux=\"none\" when using matrix.");
```
The driver calls `cfg = LoadSpatialFrictionConfig(config_path)` at `spatial_dyn_driver.cpp:616`
— so this `MFEM_VERIFY` fires during config load. The CLI override
`cfg.numerics.mixed_flux = cli_mixed_flux` is applied LATER at `:679`, and the driver guard
the plan relaxes (BUG-4) is at `:809-816` — both downstream of the abort. The Phase-7 TOML
`tpv31_rk_mixedflux.toml` is specified to set `interior_flux="matrix"` (preserved) AND
`mixed_flux="adjacent"` (the plan's "change TWO keys", `:606-608`; acceptance `:754-755`).
With BOTH keys in the file, `LoadSpatialFrictionConfig` aborts at `:1162` with
`interior_flux="matrix" is incompatible with mixed_flux="adjacent"` — the deliverable never
reaches the time loop, never reaches `SetMixedFluxMode`, never builds a central matrix. The
plan's entire Phase-5 relaxation is necessary but NOT sufficient: there are TWO guards in
series enforcing the same `matrix ⊥ mixed_flux` policy, and the plan names only the second.

(Even the CLI-only escape — keep `mixed_flux="none"` in the TOML, pass `--mixed-flux adjacent`
— does NOT match the plan, whose acceptance criterion `:754` REQUIRES the TOML diff to set
`mixed_flux: none→adjacent`; and that escape would still hit the BUG-4 driver guard, which the
plan does relax. So the config-parse guard is the *uniquely* unaddressed blocker.)

**Expected:** Phase 5 (and Phase 7) must ALSO relax the parse-level guard at
`spatial_friction.cpp:1162-1167` with the same policy as the BUG-4 driver guard: ALLOW
`interior_flux="matrix"` + `mixed_flux != "none"` (the new feature), keep rejecting it on the
ADER path (the integrator guard belongs with BUG-4's `!cfl_rk_aware_`/`time_integrator==ADER`
condition — but note the config parser does not know the resolved integrator until after the
CLI override at `:679`, so the matrix+ADER rejection must stay in the driver, and the parser
guard should simply DROP the `matrix ⊥ mixed_flux` prohibition or downgrade it to a positive
acknowledgement). The plan must name `spatial_friction.cpp:1162` explicitly as a second site
to edit, the same way BUG-4 names `:809-816`.
**Trigger:** Loading `tpv31_rk_mixedflux.toml` (matrix + adjacent) — i.e. the very first
action of the Phase-7 deliverable, local or Frontera.
**Suggested fix:** Add to Phase 5: "Relax the config-parse mutual-exclusion guard at
`spatial/code/spatial_friction.cpp:1162-1167` to allow `interior_flux="matrix"` +
`mixed_flux != "none"` (cite R-1203); keep the integrator (ADER) rejection in the driver
(Phase 5 step 1 / BUG-4), since the parser cannot yet see the resolved `--time-integrator`."
Cross-reference BUG-22 (this edit is in `spatial/code/`, outside the plan's declared scope).
**Test (CRITICAL):** Add **Test 5.3** (`Config_MatrixPlusAdjacent_ParsesWithoutAbort`):
`LoadSpatialFrictionConfig` on a tiny TOML with `interior_flux="matrix"` +
`mixed_flux="adjacent"` (+ a non-Constant material) must return a valid config (NOT abort)
after the `:1162` relaxation; the SAME TOML with `time_integrator="ader"` must still be
rejected — at the driver guard (Test 5.2), since the parser can't see the integrator. This is
the parse-level proxy for the BUG-21 fix and complements Test 5.1 (which is downstream of the
parser).

**What was fixed (Rev 5).** Plan Rev 5 adds an explicit, exhaustive **Guard inventory** to
Phase 5 (`plan:500-518`) that relaxes ALL THREE series guards in firing order — **G1 (BUG-21) =
`spatial/code/spatial_friction.cpp:1162-1167`** (config parse, "fires FIRST inside
`LoadSpatialFrictionConfig`, called at `spatial_dyn_driver.cpp:616` — UPSTREAM of the CLI
override and every later guard"), G2 = `spatial_dyn_driver.cpp:809-816` (driver), G3 =
`bimaterial_wave_operator.inl:569-574` (operator). Phase 5 req 1 (`plan:531-536`) names G1
explicitly, directs relaxing the mutual exclusion + updating the `:1160-1167` comment, cites
R-1203, and correctly keeps the matrix+ADER rejection DOWNSTREAM at G2 + `ComputeMaxDt` (because
"the final integrator is unknown until the CLI override (driver `:679`)"). Test 5.3a is the
parse-level reproducer.
**Re-verified against live source (Rev 5):** G1 is byte-exactly at `spatial_friction.cpp:1162-1167`
(`MFEM_VERIFY(cfg.numerics.interior_flux == InteriorFlux::Scalar || cfg.numerics.mixed_flux ==
"none", ...)`), with the `:1160-1161` "scalar-path-only optimization, R-1203 sibling" comment the
plan targets; the `mixed_flux` value-validator at `:1152-1157` already accepts
`adjacent`/`all_continuous`, so relaxing G1 lets them through; `LoadSpatialFrictionConfig` is
called at `spatial_dyn_driver.cpp:616`, the CLI override at `:679`, G2 at `:811`, ctor at `:1135`,
`SetMixedFluxMode` at `:1245` — all firing-order claims accurate. The "exactly THREE guards" claim
was INDEPENDENTLY re-confirmed (see the Rev-5 report's guard-sweep table — no fourth blocker).
**RESOLVED.**

---

### BUG-22 (HIGH): The `seam_continuous` field + parse + the BUG-21 guard all live in `spatial/code/spatial_friction.{hpp,cpp}` — OUTSIDE the plan's declared source-edit scope; the plan's "config parse lives with the driver" is factually wrong

**Status:** **RESOLVED in Rev 5**
**File (plan):** Constraint "Source-edit scope" (`:71-75`); Phase 7 config wiring note
(`:610-614`: "The driver reads `[material].seam_continuous` … Wiring is the spatial config
parse + one driver line — both in-scope: the config parse lives with the driver; if the
config struct is in a separate file, confirm it is within the editable tree and note it.").
**File (code):** `spatial/code/spatial_friction.hpp:573-584` (`struct MaterialSpec` — where a
`bool seam_continuous` field must be added); `spatial/code/spatial_friction.cpp:1626-1690`
(the `[material]` parse block — where the key must be read); `spatial_friction.cpp:1162-1167`
(the BUG-21 guard to relax); `makefile:132-133` (`SPATIAL_FRICTION_OBJ` is a distinct TU
linked into the driver and into `seas_test_spatial_friction_config`, `:1304-1305`).
**Severity:** HIGH — the plan's declared scope (`dynamic/` + `drivers/spatial_dyn_driver.cpp`
+ `tests/` + `Makefile`) does NOT include `spatial/code/`. Implementing `seam_continuous` and
BUG-21's relaxation REQUIRES editing `spatial/code/spatial_friction.{hpp,cpp}`, silently
expanding the no-touch boundary. The plan's claim that the parse "lives with the driver" is
false and was never resolved.

**Description:** The Rev-4 fix report self-flagged "the spatial-config parse for the new key
must be confirmed in-scope" (fix report `:62`). It is NOT confirmed in-scope; it is OUT of the
plan's *declared* scope. Verified: the `[material]` TOML block is parsed in
`spatial/code/spatial_friction.cpp` (the `if (root.contains("material"))` block at `:1626`,
material `kind` at `:1630-1635`, `depth_axis` at `:1645-1650`), and the struct is
`MaterialSpec` at `spatial_friction.hpp:573`. Neither file is the driver
(`drivers/spatial_dyn_driver.cpp`) nor under `dynamic/`. Adding `bool seam_continuous = false;`
to `MaterialSpec` and a `toml_bool(m, "seam_continuous", false)` parse line both land in
`spatial/code/`. So does the BUG-21 guard relaxation. The plan's scope Constraint (`:71-75`)
must be amended to add `spatial/code/spatial_friction.{hpp,cpp}` (and the wiring note `:610-614`
corrected — it does NOT live with the driver).

**Is this allowed by project policy?** Project memory [C2] no-touch covers
`bp5/bp1/bp2/domain/fault/solver` + `friction/dieterich_ruina.hpp`; `spatial/code/` is NOT on
that list and is already test-covered (`seas_test_spatial_friction_config`), so editing it is
not *forbidden* by project memory — but the PLAN's own narrower self-declared scope excludes
it, and the byte-exact regression contract (TPV102/104/205/BP5 + the existing TPV31 ADER job
all parse through `spatial_friction.cpp`) demands the additive field + the guard relaxation be
proven non-perturbing for every existing config. A new `bool` field defaulting `false` and a
guard that only ADDS an allowed combination is additive and safe IF done carefully, but the
plan must (a) own the scope expansion explicitly, (b) require a parser-regression check
(`seas_test_spatial_friction_config` + `test_tpv_config_parse` stay green), and (c) document
that the default `seam_continuous=false` preserves every existing config's behaviour.

**Expected:** Amend the Constraint scope list to include
`spatial/code/spatial_friction.{hpp,cpp}`; correct the Phase-7 wiring note; add a Phase-5 (or
Phase-2) requirement to add the `MaterialSpec::seam_continuous` field + parse + the BUG-21
guard relaxation, citing the schema doc
`safs/project_7.0_alternative/document/spatial_friction_config_schema.md` (where `[material]`
keys are documented, `:151`) for a doc update.
**Trigger:** Implementer reaches Phase 7 / Phase 2 and discovers the new key has no home in
the declared scope.
**Suggested fix:** As above. Fold the doc update (schema md + `tpv31.toml` comment) in.
**Test:** Covered by Test 5.3 (BUG-21) for the guard + parse; add a one-line assertion to the
config-parse test that `seam_continuous` defaults `false` and round-trips `true`.

**What was fixed (Rev 5).** Plan Rev 5 (a) **amends the source-edit scope Constraint**
(`plan:77-87`) to add `spatial/code/spatial_friction.{hpp,cpp}` explicitly, with the correct
rationale "It is NOT BP5/friction/solver core (project memory [C2] does not forbid
`spatial/code/`), but EVERY TPV*/BP5/TPV31-ADER config parses through it — so any edit MUST be
guarded by a parser-regression test"; (b) **corrects the false "config parse lives with the
driver" claim** — Phase 5 (`plan:521-525`) now says the `MaterialSpec`, its parser, and G1 all
live in that file; (c) adds `bool seam_continuous` to `MaterialSpec` + a `toml_bool(m,
"seam_continuous", false)` parse line (Phase 5 req 3, `plan:543-548`); and (d) requires a
**parser-regression test (Test 5.3c)**. **Re-verified against live source (Rev 5):** `MaterialSpec`
is at `spatial_friction.hpp:573` (fields `kind`/`profile_csv`/`sidecar_path`/`profile_layers`/
`depth_axis` — a clean place for an additive `bool seam_continuous = false;`); the `[material]`
parse block is `spatial_friction.cpp:1626-1638` (`if (root.contains("material"))`, `kind` at
`:1630-1632`, `depth_axis` at `:1645-1650`); the plan's cited `toml_bool` helper EXISTS at
`spatial_friction.cpp:304` with the exact signature `bool toml_bool(const toml::value&, const
std::string&, bool)`; project memory [C2]'s no-touch list (`bp5/bp1/bp2/domain/fault/solver` +
`friction/dieterich_ruina.hpp`) does NOT include `spatial/code/`, so the scope expansion is
project-policy-compatible, and the additive default-`false` field is byte-safe for every existing
config. **RESOLVED.** *(NB: a residual wording problem in Test 5.3c's CORPUS is tracked as the NEW
BUG-25 — it does not reopen BUG-22.)*

---

### BUG-23 (MEDIUM): The Phase-7 metric's overlap window `[t0,t1]` is undefined and partial-coverage (a truncated MFEM run) is unhandled — `np.interp` silently clamps, so the gate can spuriously fail OR silently pass

**Status:** **RESOLVED in Rev 5**
**File (plan):** Phase 7 req 5(a) (`:708-712`: "interpolate the MFEM channel onto the reference
`time_s` over the overlapping `[t0,t1]` window via `np.interp`"); Phase 7 Acceptance
(`:770-772`: the 15 s run "exits 0 within the documented ~5–10% peak band").
**File (code):** `tpv31/visualize_results.py:146,180` (MFEM vs reference independent `time_s`
arrays); `numpy.interp` (clamps out-of-range to endpoint values by default — no extrapolation,
no error). The TPV31 15 s run "spans multiple 48 h segments" via the restart chain (plan
`:682-683`), so a truncated MFEM trace (failed/incomplete restart) is a realistic outcome.
**Severity:** MEDIUM — the metric is now well-defined ON A FULL overlap (BUG-19 resolved), but
the window construction and partial-coverage behaviour are unspecified, which can make the
quantitative gate (the whole point of Phase 7) either falsely red or falsely green on a
truncated run.

**Description:** Rev-4 resolved BUG-19 by mandating `np.interp` onto the reference `time_s`
"over the overlapping `[t0,t1]` window," but never DEFINES `[t0,t1]` and never addresses the
case where the MFEM run is SHORTER than the reference (the common failure mode: 48 h wall, a
restart that didn't resume, or a CFL-driven dt collapse that never reaches `tfinal=15s`). Two
concrete hazards:
1. **Undefined window.** `[t0,t1]` is presumably `[max(mfem_t[0], ref_t[0]), min(mfem_t[-1],
   ref_t[-1])]`, but the plan does not say so. If the implementer instead interpolates the
   MFEM channel onto the FULL reference `time_s`, `np.interp` CLAMPS MFEM values beyond
   `mfem_t[-1]` to the last MFEM sample — comparing a flat-lined MFEM tail against the live
   reference ⇒ a large spurious `peak_rel` ⇒ false FAIL on an otherwise-good (just truncated)
   run.
2. **Silent pass on tiny overlap.** Conversely, if the metric restricts to `[t0,t1]` but the
   overlap is tiny (MFEM died at t=2 s), peak/rms over a 2 s window that may not even include
   the rupture front passes trivially ⇒ false PASS.
There is no minimum-coverage requirement (e.g. "fail if MFEM coverage < X% of the reference
span" or "require `t1 ≥ t_rupture_arrival`").

**Expected:** In req 5(a): (1) define `t0 = max(mfem_t[0], ref_t[0])`, `t1 = min(mfem_t[-1],
ref_t[-1])` and restrict BOTH series to `[t0,t1]` (interpolate the MFEM channel onto the
reference samples that fall in `[t0,t1]`); (2) add a minimum-coverage guard — fail (or
explicitly WARN-and-skip with a non-zero exit) if `t1 - t0` is less than a stated fraction of
the reference span (e.g. < 90%), so a truncated run cannot silently pass; (3) state that the
gate is only meaningful on a completed run, mirroring the Frontera-acceptance criterion.
**Trigger:** A real TPV31 run that does not reach `tfinal=15s` (wall-clock or restart failure),
then `visualize_results.py --tol-rms --tol-peak`.
**Suggested fix:** Add the window definition + coverage guard to req 5(a) and assert it in a
new Test NEW-7.5 (truncated-MFEM fixture ⇒ non-zero exit / explicit insufficient-coverage
failure, NOT a silent pass).

**What was fixed (Rev 5).** Plan Phase 7 req 5(a) now defines the window and coverage handling
fully and self-consistently (`plan:757-769`): (1) explicit window `t0 = max(mfem_t[0], ref_t[0])`,
`t1 = min(mfem_t[-1], ref_t[-1])`, restrict BOTH series to `[t0,t1]`, and **interpolate the MFEM
channel onto the reference samples within `[t0,t1]` via `np.interp` BEFORE** differencing — with
the explicit warning "Do NOT interpolate onto the FULL reference grid: `np.interp` clamps
out-of-range to the endpoint, so a truncated MFEM run would flat-line its tail and spuriously
FAIL"; (2) a **minimum-coverage guard** — "If `(t1 − t0)` is less than a stated fraction of the
reference span (e.g. < 90%, or `t1 < t_rupture_arrival`), fail with an explicit 'insufficient
coverage' message and a non-zero exit"; and (3) Test NEW-7.5 (`TolGate_TruncatedRun_NoSilentPass`)
as the truncated-run reproducer. **Cross-checked for self-consistency across all four sites:** req
5(a) (`:757-769`), the Files-to-Create entry (`:669-674`, "covers far less than the reference span
⇒ non-zero exit with an explicit insufficient-coverage failure, NOT a silent pass"), the Acceptance
Criterion (`:821-825`, "minimum-coverage guard … NEW-7.5 (truncated run ⇒ insufficient-coverage
exit≠0)"), and the Testing-Strategy row (`:850`) — all agree. The peak-normalized, named-channel
metric (BUG-17/18) and the `np.interp` mandate (BUG-19) remain intact. **RESOLVED.**

---

### BUG-24 (MODERATE): The plan does not state EXPLICITLY that the `seam_continuous` guard leaves the existing TPV31 ADER / `mixed_flux=none` job byte-exact (it does, structurally — but the regression claim should be made explicit)

**Status:** **RESOLVED in Rev 5**
**File (plan):** Constraint "Byte-exact regression contract" (`:40-45`, mentions the ADER job
at the *config* level but not the *runtime guard*); Phase 2 req 4 (`:272-296`, describes the
guard but does not assert non-interference with `mixed_flux=none`).
**File (code):** `dynamic/wave_operator.inl:1653-1654` (`BuildCentralFluxFaceSet_` returns
immediately, clearing `central_flux_face_set_`, when `mixed_flux_mode_ == None`);
`wave_operator.inl` `mf_on_ = (m != None)`; Phase-2 `BuildPerFaceCentralFluxMatrices_` iterates
`central_flux_face_set_` (req 3) and the shared arm intersects `central_flux_face_set_.count(...)`
(plan `:264`).
**Severity:** MODERATE — this is a documentation/assurance gap, not a functional defect. The
guard is structurally safe (see below), but given the user's explicit Part-B.2 concern and the
hard byte-exact contract, the plan should say so.

**Description:** The existing parallel TPV31 ADER bimaterial job
(`jobs/tpv31_spatial/tpv31_p2_aderO3_noflux_50m_normal.sbatch`, `tpv31.toml` with
`mixed_flux="none"`, `interior_flux="matrix"`) does NOT set `seam_continuous`. Could the new
Phase-2 guard newly ABORT it? **No — verified structurally:** with `mixed_flux="none"`, the
base `SetMixedFluxMode(None)` sets `mf_on_=false` and `BuildCentralFluxFaceSet_` returns with
an EMPTY `central_flux_face_set_` (`wave_operator.inl:1653-1654`). The Phase-2
`BuildPerFaceCentralFluxMatrices_` iterates that empty set (interior arm) and intersects the
shared arm with `central_flux_face_set_.count(mesh_face)` (always 0) — so the `seam_continuous`
guard is NEVER reached on a `mixed_flux=none` run. The ADER job stays byte-exact, AND the BUG-21
config-parse guard does not fire for it (`mixed_flux="none"` is allowed with `matrix`). So the
guard fires ONLY for `Mode::Coefficient` + central-SHARED face + `mf_on_` true — exactly the new
path. **The concern is correctly handled by construction; it is just not stated.**
**Expected:** Add one sentence to Phase 2 req 4 (and the Constraint) asserting that the
`seam_continuous` guard fires ONLY inside `BuildPerFaceCentralFluxMatrices_` for faces in
`central_flux_face_set_` with `mf_on_` true, so `mixed_flux=none` (the existing TPV31 ADER job)
never triggers it and stays byte-exact.
**Trigger:** N/A (no functional failure) — a clarity/assurance gap only.
**Suggested fix:** As above. Tie it to the existing parity/regression acceptance criteria.

**What was fixed (Rev 5).** Plan Rev 5 makes the no-regression statement explicit in BOTH places
the user's Part-B.2 concern targeted: the **byte-exact regression Constraint** (`plan:46-51`) now
reads "**(BUG-24)** The new `seam_continuous` guard cannot regress that ADER job: with
`mixed_flux="none"`, `BuildCentralFluxFaceSet_` returns an empty set (`wave_operator.inl:1653-1654`)
and `mf_on_=false`, so the Phase-2 `BuildPerFaceCentralFluxMatrices_` (which walks only
`central_flux_face_set_`) is never entered and the guard never fires — and G1 already allows
`matrix + mixed_flux=none`. The guard fires ONLY for `Mode::Coefficient` + a central-set shared
face with `mf_on_` true."; and **Phase 2 req 4** (`plan:304-308`) carries the same `(BUG-24)`
statement. **Re-verified against live source (Rev 5):** `BuildCentralFluxFaceSet_`
(`wave_operator.inl:1650-1653`) does `central_flux_face_set_.clear(); if (mixed_flux_mode_ ==
MixedFluxMode::None) { return; }`, and `mf_on_ = (m != None)` at `:1640` — so on the existing
`mixed_flux="none"` TPV31 ADER job the central set is empty, the Phase-2 build is never entered, the
guard is unreachable, and G1 passes (`matrix + none` is allowed). The plan's cite `:1653-1654` is
one line off from the live `:1652-1653` but substantively exact. **RESOLVED.**

---

## Bugs Found — Rev 5

### BUG-25 (MEDIUM): Test 5.3c's parser-regression CORPUS is mis-worded and partly infeasible-as-written — there is no `SpatialFrictionConfig` equality/dump, and BP5 has no spatial TOML to parse

**Status:** OPEN (Rev 5, NEW)
**File (plan):** Phase 5 Test 5.3c (`plan:567-571`: "the existing TPV102/104/205 + BP5 +
TPV31-ADER TOMLs parse to byte-identical `SpatialFrictionConfig`"); Constraint (`plan:84`: "EVERY
TPV*/BP5/TPV31-ADER config parses through it"); Phase 7 acceptance (`plan:815-818`); Testing
Strategy row (`plan:849`).
**File (code):** `spatial/code/spatial_friction.hpp:573-687` (`MaterialSpec` /
`SpatialFrictionConfig` — NO `operator==`, NO `Dump`/`ToString`/`Serialize`); `tests/unit/
test_tpv_config_parse.cpp:281-283` (the existing parser test loads ONLY
`tpv205/102/104_spatial.toml` via `LoadSpatialFrictionConfig` and asserts FIELD-BY-FIELD with
`TEST_ASSERT`); the worktree has NO `bp5/configs/*.toml` (BP5 runs the native `seas_bp5_full`
driver with CLI/hardcoded params, not `LoadSpatialFrictionConfig`).
**Severity:** MEDIUM — the test is REALIZABLE, but two factual problems in its specification will
mislead the implementer if taken literally. NOT a functional blocker; the underlying guard
(BUG-22 parser-regression) is sound.

**Description:** Test 5.3c is specified as proving the listed configs "parse to byte-identical
`SpatialFrictionConfig` (no field shifted, no new abort)". Two problems:
1. **No equality/dump mechanism.** `SpatialFrictionConfig` (and `MaterialSpec`) have NO
   `operator==`, no `Dump`/`ToString`/`Serialize`. There is no way to compare two parsed configs
   "byte-identically" without either (a) writing a field-by-field comparison (the established
   `test_tpv_config_parse.cpp` pattern — load → `TEST_ASSERT` each field) or (b) adding a
   serializer (out of scope). "Byte-identical `SpatialFrictionConfig`" overstates the available
   mechanism; the realizable form is "parse without abort + the additive `seam_continuous`
   defaults `false` + every previously-asserted field is unchanged."
2. **BP5 has no spatial TOML.** The worktree has no `bp5/configs/*.toml`; BP5 does not parse
   through `LoadSpatialFrictionConfig` at all (it is a separate native driver). So BP5 cannot be in
   a "parse byte-identically" corpus. BP5's byte-exactness against this change is protected by the
   change being purely additive to a parser BP5 never invokes — NOT by Test 5.3c. Including BP5 in
   the parse corpus is a factual error that would send the implementer looking for a non-existent
   file.
   - The TPV31-ADER config (`tpv31/configs/tpv31.toml`) DOES exist and IS parseable (it is already
     loaded by `test_tpv31_canonical_rotation.cpp:45-48`), but `test_tpv_config_parse.cpp` does NOT
     currently include it — Test 5.3c must ADD it (and tpv31_p2/p3) to that test's corpus.

**Expected:** Reword Test 5.3c to: (a) the parse-regression corpus is the configs that actually
parse through `LoadSpatialFrictionConfig` — `tpv205/102/104_spatial.toml` (already in
`test_tpv_config_parse.cpp`) PLUS `tpv31/configs/tpv31.toml` (+ p2/p3) which must be ADDED; DROP
BP5 from the parse corpus (note instead that BP5 is unaffected because it does not use this
parser); (b) realize "byte-identical" as field-by-field `TEST_ASSERT`s (the existing pattern) —
specifically that the additive `seam_continuous` defaults `false` and every field the test already
asserts is unchanged — since no `operator==` exists. Keep the BUG-22 intent (prove the additive
field + G1 relaxation do not perturb existing parses) intact.
**Trigger:** Implementer writes Test 5.3c against the plan's literal corpus — looks for
`bp5/configs/*.toml` (absent) and a `SpatialFrictionConfig::operator==` (absent).
**Suggested fix:** As above. One-line plan edit to the Test 5.3c corpus + the Constraint sentence
("EVERY TPV*/BP5/TPV31-ADER config parses through it" → "every TPV*/TPV31 *spatial* config parses
through it; BP5 uses its own native driver and is unaffected").
**Test:** This finding IS about Test 5.3c — fold the corpus correction into its spec; no new test
ID needed.

---

### BUG-26 (LOW): The plan never states TPV31 is `law="slip_weakening"` (LSW) — so the actual RK path is `AdvanceRKCoupledLSW_Spatial`, not the rate-state stepper the Test-5.1/6.1 framing implies

**Status:** OPEN (Rev 5, NEW)
**File (plan):** Phase 7 (`plan:640-667`, `:705-717` — `--time-integrator rk45` with no statement
of the friction law / which coupled-RK stepper is exercised); Test 5.1
(`Driver_MatrixPlusAdjacentPlusRK_ReachesTimeLoop`); Test 6.1 (homog-equivalence, rate-state-style
framing).
**File (code):** `tpv31/configs/tpv31.toml:58` (`law = "slip_weakening"`), `:77`
(`kind = "tpv31_lsw"`), `:270` (`[friction.slip_weakening]`), `:304`
(`kind = "instantaneous_overstress_circular"`); `drivers/spatial_dyn_driver.cpp:2864-2885` (RK
dispatch: `is_rk && is_lsw` → `AdvanceRKCoupledLSW_Spatial`, else `is_rk` →
`AdvanceRKCoupled_Spatial`); `dynamic/rk_time_stepper.hpp:369-392`
(`AdvanceRKCoupledLSW_Spatial`, hard precondition `GetFaultFrictionLaw()==LSW`).
**Severity:** LOW — the feature WORKS for LSW (the central dispatch is in the law-agnostic `Mult`,
exercised identically by both coupled-RK steppers), and the existing TPV31 ADER job already runs
LSW + matrix. This is a documentation/clarity gap, not a functional defect.

**Description:** TPV31 is a SLIP-WEAKENING (LSW) benchmark (`tpv31.toml:58 law="slip_weakening"`),
not rate-and-state. On the RK path the driver therefore dispatches to `AdvanceRKCoupledLSW_Spatial`
(`spatial_dyn_driver.cpp:2864,2873`), whose only hard precondition is
`GetFaultFrictionLaw()==LSW` (`rk_time_stepper.hpp:389`) — satisfied — with NO flux/matrix/operator
precondition. The plan's Phase-7 CLI uses `--time-integrator rk45` but never states the friction
law or which stepper runs; the check-doc test framing (Test 5.1/6.1) and the fix-report discussion
lean rate-state. An implementer could write the local smoke/Test-5.1 fixture as rate-state and be
surprised the actual deliverable path is the LSW stepper, or could wrongly assume the
`AdvanceRKCoupled_Spatial` RS precondition (`GetFaultFrictionLaw()==RateAndState`,
`rk_time_stepper.hpp:204`) gates TPV31. Both steppers route through `wave.Mult`, where the Phase-3
central branch lives, so the central-flux feature is law-agnostic and correct for TPV31-LSW.
**Expected:** Add one sentence to Phase 7 (and the homog/smoke test notes) stating TPV31 is LSW, so
the RK path is `AdvanceRKCoupledLSW_Spatial`, and that the central dispatch is law-agnostic (lives
in `Mult`) so the feature applies unchanged. Optionally note that `nucleation =
instantaneous_overstress_circular` (parsed at `spatial_friction.cpp:1355`) is applied via
`nuc->ApplyAbsolute` in the LSW RK stepper and that forced-rupture is the only rejected LSW
sub-case (driver `:802-805`) — not relevant here.
**Trigger:** Implementer builds the Phase-7 deliverable / local smoke assuming rate-state.
**Suggested fix:** As above — a one-sentence clarification; no code consequence.
**Test:** The existing Test 5.1 should construct its tiny fixture as LSW (or note both laws reach
the time loop) so it mirrors the real deliverable; no new test ID required.

---

## Proposed Unit Tests

(Rev-1 tests retained; Rev-2 status / new tests flagged.)

### Test Suite 1: `tests/unit/test_bimaterial_central_flux.cpp` (Phase 1 primitive)

**Test 1.1: CentralFlux_HomogeneousLimit_MatchesGodunovCentral** — HIGH.
homogeneous limit equals `GodunovFlux::Central` to ≤1e-11 (Phase 1 AC). *(Plan Test 1.1.)*

**Test 1.2: CentralFlux_Consistency_PhysicalFlux** — MEDIUM.
`Q_self==Q_nbr==Q` ⇒ `F* == A·Q` to ≤1e-11 (Phase 1 AC). *(Plan Test 1.2.)*

**Test 1.3: CentralFlux_HeterogeneousAveraging** — HIGH (anchors BUG-3 at primitive level).
`A_e1 != A_e2`: `ApplyPerFaceFlux(cS,cN,Q,0,F)==½A_e1·Q`, `(cS,cN,0,Q,F)==½A_e2·Q`.
*(Plan Test 1.3.)*

### Test Suite 2: `tests/unit/test_bimaterial_mixed_flux_dispatch.cpp` (Phase 3)

**Test 3.1: Dispatch_CentralVsGodunov_PerFaceSelection** — HIGH. *(Plan Test 3.1.)*

**Test 3.2: Dispatch_CentralFlux_SingleValued_UnderHeterogeneity** — HIGH
(**BUG-3/BUG-7 load-bearing guard**). On `mu_e1 != mu_e2`, `F_h_e1 == F_h_e2` ≤1e-12·scale.
Passes with the Rev-2 side-symmetric store. *(Plan Test 3.2.)*

**Test 3.3: Dispatch_None_ByteExactGodunov** — HIGH. *(Plan Test 3.3.)*

**Test 3.4 (MPI, np=2): Dispatch_SharedCentralFace** — MEDIUM
(BUG-5 OK / **BUG-10 broken**). Builds + dispatches the shared central matrices with the
correct neighbour material (BUG-5 — works). The "aborts on a laterally-varying material"
half **cannot pass** against Phase 2 req 4's dead predicate (BUG-10). **Action:** split
into 3.4a (shared-central build/dispatch, BUG-5) and 3.4b (lateral-het abort) — and 3.4b
must be written against whatever real detector replaces the dead `neighbour vs local`
compare (see BUG-10). *(Plan Test 3.4 — revise.)*

**Test 3.4b (Rev-3 update): MUST be written against the `seam_continuous` flag, NOT the
material-field probe.** Per BUG-15 the probe is unimplementable (no `MaterialField::EvalAt(x)`)
and per BUG-16 it false-aborts on TPV31's own depth profile. A correct 3.4b: a
`Mode::Coefficient` central shared face WITHOUT `material.seam_continuous=true` aborts with
the seam-named message; WITH it, builds. A probe-based 3.4b cannot be authored against the
real API. (NEW in Rev 3.)

**Test 3.4c (NEW in Rev 3): Dispatch_DepthProfile_NoFalseAbort** — HIGH (BUG-16 reproducer).
Target: `BimaterialWaveOperator::BuildPerFaceCentralFluxMatrices_` (np=2,
`tests/parallel/`). A depth-varying `Mode::Coefficient` material with
`seam_continuous=true` on a fixture whose partition seam normal has a z-component must BUILD
(no abort); the SAME material WITHOUT the flag aborts. **Fails against the along-normal
probe** (which aborts on the legitimate depth jump) — proves the probe is a false-positive
blocker and the `seam_continuous` mechanism is correct.

### Test Suite 3: `tests/unit/test_bimaterial_mixed_flux_cfl.cpp` (Phase 4)

**Test 4.1: ComputeMaxDt_FactorMatchesScalar_AllModes** — HIGH (BUG-2 guard). *(Plan 4.1.)*

**Test 4.2: ComputeMaxDt_CentralPlusADER_Aborts** — HIGH. *(Plan 4.2.)*

### Test Suite 4: `tests/unit/test_bimaterial_mixed_flux_homog_equivalence.cpp` (Phase 6 GATE)

**Test 6.1: HomogeneousEquivalence_Adjacent_And_AllContinuous** — HIGH.
necessary, not sufficient (passes with BUG-3 present; 3.2 is the sufficient guard).
*(Plan Test 6.1.)*

### Test Suite 5: Driver (Phase 5)

**Test 5.1 (dry-run): Driver_MatrixPlusAdjacentPlusRK_ReachesTimeLoop** — HIGH (BUG-4).
*(Plan Test 5.1.)* **(Rev-5 note, BUG-26):** the deliverable TPV31 is `law="slip_weakening"`, so
the real RK path is `AdvanceRKCoupledLSW_Spatial` — construct the tiny fixture as LSW (or assert
BOTH laws reach the time loop), not rate-state, so the dry-run mirrors the actual path.

**Test 5.2 (dry-run): Driver_MatrixPlusAdjacentPlusADER_Aborts** — MEDIUM. *(Plan 5.2.)*

**Test 5.3 (NEW in Rev 4; split a/b/c in plan Rev 5): Config_MatrixPlusAdjacent_ParsesWithoutAbort**
— HIGH (**BUG-21/BUG-22 reproducer**) [SPEC FIXED in plan Rev 5]. Target: `LoadSpatialFrictionConfig`
(`spatial/code/spatial_friction.cpp`). Validates: **(a)** a tiny TOML with `interior_flux="matrix"` +
`mixed_flux="adjacent"` + a non-Constant material PARSES (does NOT abort at `:1162`) after the G1
relaxation; **(b)** `[material].seam_continuous` defaults `false` and round-trips `true`;
**(c)** parser-regression — the existing configs that go through `LoadSpatialFrictionConfig` still
parse unchanged. **Fails today — `spatial_friction.cpp:1162-1167` aborts `matrix + mixed_flux!="none"`
at config load, upstream of every site the plan touches.** Runs in `seas_test_spatial_friction_config`
/ `test_tpv_config_parse`.
**(Rev-5 caveat, BUG-25 — MUST correct before implementing 5.3c):** there is NO
`SpatialFrictionConfig::operator==`/dump, so "byte-identical" must be realized field-by-field (the
existing `test_tpv_config_parse.cpp` pattern); and **BP5 has NO spatial TOML** (it never parses
through this function) so it CANNOT be in the parse corpus — the realizable corpus is
`tpv205/102/104_spatial.toml` (already in `test_tpv_config_parse.cpp:281-283`) PLUS
`tpv31/configs/tpv31.toml` (+ p2/p3), which Test 5.3c must ADD. Assert each config parses without
abort and `seam_continuous` defaults `false`.

### Test Suite 6: Phase-7 tooling (NEW in Rev 2)

**Test NEW-7.1: VisualizeResults_TolGate_ExitsNonZeroOnExceedance** *(BUG-11 reproducer)* —
HIGH. Target: `tpv31/visualize_results.py`. Validates: with `--tol-rms`/`--tol-peak`
added, a synthetic MFEM trace that deviates >peak-band from a reference `.txt` makes the
script exit non-zero (and a within-band trace exits 0). **Fails today — the flags do not
exist.** This is the missing quantitative gate the Phase-7 acceptance criterion assumes.
Skeleton:
```python
# build a 2-col synthetic MFEM station .dat + a scec_seisol .txt that differ by 20% peak
# run: visualize_results.py <mfem_dir> --seisol --tol-peak 0.10 --tol-rms 0.10
# assert returncode != 0; then make them differ by 3% and assert returncode == 0
```

**Test NEW-7.1 (Rev-3 update): the synthetic fixtures MUST use a peak-normalized,
interpolated metric, NOT a per-sample relative error.** Per BUG-17 the denominator must be
`max|ref|` (zero-safe), per BUG-18 the gated channel must be named (e.g. `V_strike`), and per
BUG-19 the synthetic MFEM `.dat` and reference `.txt` must be on DIFFERENT time samples so the
test exercises `np.interp` resampling onto a common time base. A fixture on aligned samples
passes trivially and leaves the production gate broken on misaligned (RK45 vs SeisSol) data.

**Test NEW-7.2: BenchmarkReadme_MatchesCommittedData** *(BUG-12 guard)* — LOW.
Validates: the README's described format (extension, columns) matches the committed
`scec_seisol/*.txt` and what `load_reference_file` parses (`min_cols=8`). A drift guard
so the docs and data cannot diverge again.

**Test NEW-7.3: Sbatch_BannerGrep_FindsResolvedScheme** *(BUG-13 guard)* — LOW.
Validates: a `--dry-run` (tiny fixture) `*.out` capture contains regex-matchable
`time integrator:`, `mixed flux:`, and `interior flux :` lines (whitespace-insensitive),
mirroring the proven TPV104/TPV31 sbatch grep. Confirms the gate strings exist where the
sbatch will look. (Driver-level; can be a scripted check.)

**Test NEW-7.4 (NEW in Rev 3): TolGate_PeakNormalized_ZeroCrossingSafe** *(BUG-17/19 guard)*
— MEDIUM. Target: the new `--tol-rms/--tol-peak` metric in `visualize_results.py`. Validates:
a reference channel that passes through zero (pre-nucleation zeros + a slip-rate pulse) does
NOT produce spurious failures; the peak-normalized error of a within-band MFEM trace stays
below tol despite the zero-crossings (a per-sample `|mfem−ref|/|ref|` form would spike to ∞ at
the zeros and falsely fail). Also asserts the metric interpolates correctly when the two
traces are on offset time grids.

**Test NEW-7.5 (NEW in Rev 4): TolGate_TruncatedRun_NoSilentPass** *(BUG-23 guard)* — MEDIUM.
Target: the `--tol-rms/--tol-peak` window/coverage logic in `visualize_results.py`. Validates:
when the synthetic MFEM trace ends well before the reference (a truncated 48 h-wall run), the
gate either (a) restricts to the true overlap `[t0,t1] = [max(starts), min(ends)]` AND fails
the minimum-coverage guard (non-zero exit / explicit insufficient-coverage message), or (b)
does NOT clamp the MFEM tail via `np.interp` and spuriously fail a within-band run. Asserts the
window is the intersection (not the full reference grid) and that a <90%-coverage MFEM trace
cannot silently PASS.

---

## Priority Summary

*(Rev-5 statuses appended. ALL FOUR Rev-4 findings (BUG-21 CRITICAL, BUG-22 HIGH, BUG-23 MEDIUM,
BUG-24 MODERATE) are RESOLVED in plan Rev 5 — every fix re-verified byte-accurate against the
still-unchanged live source, and the "exactly THREE guards" inventory INDEPENDENTLY confirmed by a
fresh sweep of the RK time-stepper / ADER substep / ctor / integrator-config / fault-iterator paths
(no fourth blocker). The six Rev-3 findings (BUG-15/16/17/18/19/20) and the two long-running
CRITICALs (BUG-10, BUG-11) remain RESOLVED. TWO NEW Rev-5 findings, both NON-blocking: BUG-25
MEDIUM (Test 5.3c corpus mis-worded — no config `operator==`, BP5 has no spatial TOML), BUG-26 LOW
(plan silent on TPV31 being LSW). **NO CRITICAL/HIGH finding is OPEN.**)*

| ID | Severity | Status | Description |
|----|----------|--------|-------------|
| BUG-21 | CRITICAL | **RESOLVED in Rev 5** | Plan Rev 5 adds the exhaustive 3-guard inventory (`plan:500-518`) and relaxes G1 `spatial_friction.cpp:1162-1167` (Phase 5 req 1, cite R-1203), keeping matrix+ADER rejection downstream at G2+`ComputeMaxDt`. Re-verified: G1 byte-exact at `:1162-1167`; `LoadSpatialFrictionConfig`@616, CLI@679, G2@811, ctor@1135, SetMixedFluxMode@1245. Three-guard claim INDEPENDENTLY confirmed (no 4th). Test 5.3a |
| BUG-22 | HIGH | **RESOLVED in Rev 5** | Plan Rev 5 amends the scope Constraint (`plan:77-87`) to add `spatial/code/spatial_friction.{hpp,cpp}`, corrects the false "lives with the driver" claim, adds `bool seam_continuous` to `MaterialSpec` + parse, requires parser-regression (Test 5.3c). Re-verified: `MaterialSpec`@hpp:573, parse@cpp:1626-1638, `toml_bool`@cpp:304 (exact sig). [C2] does not forbid `spatial/code/`. (Test-5.3c CORPUS wording → BUG-25) |
| BUG-23 | MEDIUM | **RESOLVED in Rev 5** | Phase 7 req 5(a) (`plan:757-769`) now defines `t0=max(starts)`, `t1=min(ends)`, interpolate-within-`[t0,t1]` via `np.interp` (NOT full grid), + a min-coverage guard (fail/non-zero exit if <~90%). Self-consistent across req 5(a), Files-to-Create (`:669-674`), Acceptance (`:821-825`), Testing-Strategy (`:850`). Test NEW-7.5 |
| BUG-24 | MODERATE | **RESOLVED in Rev 5** | The explicit no-regression statement is now in BOTH the byte-exact Constraint (`plan:46-51`) and Phase 2 req 4 (`plan:304-308`): the guard fires ONLY inside `BuildPerFaceCentralFluxMatrices_` for `mf_on_`+central-set faces; `mixed_flux=none` empties the set (`wave_operator.inl:1652-1653`) ⇒ ADER job byte-exact. Re-verified by construction |
| BUG-25 | MEDIUM | OPEN (Rev 5, NEW) | Test 5.3c corpus mis-worded: NO `SpatialFrictionConfig::operator==`/dump (so "byte-identical" must be field-by-field, the `test_tpv_config_parse` pattern), AND BP5 has no spatial TOML (never parses through `LoadSpatialFrictionConfig`). Realizable corpus = tpv205/102/104 + tpv31 spatial configs; drop BP5. Fix = reword Test 5.3c + the Constraint sentence. NOT a blocker |
| BUG-26 | LOW | OPEN (Rev 5, NEW) | Plan silent on TPV31 being `law="slip_weakening"` (LSW), so the real RK path is `AdvanceRKCoupledLSW_Spatial` (not the RS `AdvanceRKCoupled_Spatial` the Test-5.1/6.1 framing implies). Central dispatch is in the law-agnostic `Mult` ⇒ feature works for LSW. Add a one-sentence clarification; no code consequence |
| BUG-15 | CRITICAL | **RESOLVED in Rev 4** | probe DELETED everywhere; declarative `seam_continuous` flag mandated (needs no arbitrary-point eval, no cross-rank comm). Plan `:88-99,272-296`; Test 3.4b rewritten against the flag |
| BUG-16 | CRITICAL | **RESOLVED in Rev 4** | with the along-normal probe gone, the depth-vs-lateral conflation cannot occur; flag does not inspect material. Plan `:744-749`; Test 3.4c (no-false-abort) added |
| BUG-19 | HIGH | **RESOLVED in Rev 4** (residual → BUG-23) | req 5(a) mandates `np.interp` onto reference `time_s` over `[t0,t1]` BEFORE differencing (`:708-712`); NEW-7.1 uses offset grids. Window/partial-coverage still under-specified → BUG-23 |
| BUG-17 | MEDIUM | **RESOLVED in Rev 4** | req 5(a) specifies zero-safe `peak_rel`/`rms_rel` with denominator floor (`:718-722`); NOT per-sample. Test NEW-7.4 (zero-crossing) added |
| BUG-18 | MEDIUM | **RESOLVED in Rev 4** | req 5(a) names channels `V_strike`/`slip_strike`/`tau_strike`/`sigma_n` (+`_dip` if present), EXCLUDES `mu_eff` (all-NaN ref `:191`). VERIFIED: names == script dict keys `visualize_results.py:145-155,179-192,205-214` |
| BUG-20 | LOW | **RESOLVED in Rev 4** (dissolved) | probe deleted ⇒ reflection `δ`/containment/`1e-6` tol no longer part of any mechanism |
| BUG-10 | CRITICAL/HIGH | **RESOLVED in Rev 4** | dead predicate (Rev 2) → broken probe (Rev 3) → declarative `seam_continuous` flag (Rev 4): implementable, communication-free, immune to depth/lateral conflation. Deferred `MPI_Allgatherv` remains the fully-correct path |
| BUG-11 | CRITICAL | **RESOLVED in Rev 4** | flag gap closed (Rev 3) + metric fully specified (Rev 4: interp + peak-norm + named channels). Implementable end-to-end gate (residual window gap → BUG-23) |
| BUG-12 | MEDIUM | **RESOLVED in Rev 3** | Phase 7 req 5(b) rewrites README to 60 `.txt` / 8-col / SeisSol header; cites script parser contract |
| BUG-13 | MEDIUM | **RESOLVED in Rev 3** | req 4 attributes banner (`:853,:859`) + verify-dispatch (`:2130-2147`) correctly; tolerant whitespace-insensitive grep |
| BUG-14 | LOW | **RESOLVED in Rev 3** | req 3 quotes real `:1140-1148` ("works for BOTH operators; geometry-only; R-002; R-004"); OMIT re-justified |
| BUG-3 | CRITICAL | **FIXED in Rev 2** (re-confirmed Rev 3) | side-symmetric single pair + `F_h_e2←F_h_e1`; single-valued under heterogeneity; dispatch `inl:471-495` unchanged |
| BUG-4 | CRITICAL | **FIXED in Rev 2** (re-confirmed Rev 3) | driver guard `:809-816` relaxed to allow matrix+RK+mixed!=none, reject matrix+ADER+mixed!=none; guard text unchanged in live source |
| BUG-2 | HIGH | **FIXED in Rev 2** | delete `inl:521-531`; route both ops through `MixedFluxCflFactor_()`; per-element `dt_e` placement preserves MPI_MIN |
| BUG-1 | MEDIUM | **FIXED in Rev 2** | abort body `inl:567-576` + 4 stale comments now correctly cited |
| BUG-5 | MEDIUM | **FIXED in Rev 2** | `sf`-loop shared arm mirroring upwind Pass 2 `inl:375-422` |
| BUG-6 | MEDIUM | **RESOLVED in Rev 4** (by supersession) | enforcement is now the declarative `seam_continuous` flag (BUG-10/15/16 Rev-4); the supported regime (Constant / seam-continuous) is explicit and fail-loud. NB the flag PARSE lands in `spatial/code/` (BUG-22) and co-trips BUG-21 |
| BUG-7 | MEDIUM | **FIXED in Rev 2** | single-valuedness moved to operator-level het Test 3.2; trivial primitive test dropped |
| BUG-8 | LOW | **FIXED in Rev 2** | 4-edit Makefile recipe + obj-dep lists; all cites accurate |
| BUG-9 | LOW | **FIXED in Rev 2** | scratch tmp1/tmp2 + SetSize output; acoustic `mu>eps` verify; matches `Central` `cpp:416` |

### Test Coverage Summary

| Suite | File Under Test | Tests | Key Bugs Detected |
|-------|----------------|-------|-------------------|
| 1 | `godunov_flux_bimaterial.*` (primitive) | 3 | BUG-3 (partial, via 1.3) |
| 2 | `bimaterial_wave_operator` dispatch | 6 | **BUG-3 (3.2)**, BUG-5 (3.4a), **BUG-10/15 (3.4b — `seam_continuous`)**, **BUG-16 (3.4c, NEW)** |
| 3 | `bimaterial_wave_operator` CFL | 2 | BUG-2 (4.1) |
| 4 | full-operator equivalence GATE | 1 | (necessary, not sufficient for BUG-3) |
| 5 | driver wiring + config parse | 3 | **BUG-4 (5.1)**, **BUG-21/22 (5.3, NEW Rev 4)** |
| 6 | Phase-7 tooling | 6 | **BUG-11 (7.1)**, BUG-12 (7.2), BUG-13 (7.3), **BUG-17/19 (7.4)**, **BUG-23 (7.5, NEW Rev 4)** |
| **Total** | | **21** | |

*(Suite 2: 3.1, 3.2, 3.3, 3.4a, 3.4b, 3.4c = 6. Suite 5: 5.1, 5.2, 5.3 = 3. Suite 6:
NEW-7.1..7.5 = 5 plus the BUG-18 column-naming assertion folded into 7.1. New tests in Rev 4:
**Test 5.3** (`Config_MatrixPlusAdjacent_ParsesWithoutAbort`, BUG-21/22) and **Test NEW-7.5**
(`TolGate_TruncatedRun_NoSilentPass`, BUG-23). Rev-3-added: Test 3.4c, Test NEW-7.4, the Test
3.4b rewrite against `seam_continuous`.)*

*(Rev-5: NO new test ID added — total stays **21**. The two NEW Rev-5 findings fold into EXISTING
tests: BUG-25 corrects Test 5.3c's CORPUS (drop BP5; add `tpv31.toml`; field-by-field since no
`operator==`); BUG-26 directs Test 5.1's fixture to be built as LSW. Both are spec corrections to
already-counted tests, not new tests.)*

---

## Code Review Report (Rev 1)

*(Unchanged — see the Rev-1 record below for the full original report.)*

**Date:** 2026-06-04
**Scope:** plan `PLAN_mixed_flux_hetero_riemann.md` vs existing `miniapps/seas/` source.
**Reviewer:** chunhui-code-reviewer agent (pre-implementation plan audit).

### Summary (Rev 1)
The plan's overall design — bi-material central flux `F* = ½(A_self Q_self + A_nbr Q_nbr)`
dispatched per-face alongside the bi-material Godunov upwind, abort lifted, CFL de-rated,
gated on a homogeneous-limit equivalence — is sound and matches drdg3d's `get_flux`.
The Rev-1 concrete prescription re-introduced the per-side-Jacobian pitfall (BUG-3) and
missed the load-bearing driver guard (BUG-4); BUG-5/BUG-6 concerned the parallel
shared-face path; BUG-2 (CFL "pure refactor" claim) and BUG-1 (mis-cited abort) rounded
out the correctness/accuracy findings. The six mechanically-verifiable claims (items
2,4,5,6 and the verifiable halves of 8,10) all checked out exactly.

### Code Review Findings (Rev 1, against the plan)

| # | File:line (claim vs actual) | Severity | Description |
|---|------|----------|-------------|
| 1 | plan→`wave_operator.inl:1623-1626` / actual `bimaterial_wave_operator.inl:567-576` | MEDIUM | abort is an override in a different file (BUG-1) |
| 2 | `bimaterial_wave_operator.inl:521-531` vs `wave_operator.inl:5886-5891` | HIGH | divergent CFL switch; no `cfl_rk_aware_` (BUG-2) |
| 3 | plan Phase 2/3 vs `inl:471-484`, `godunov_flux_bimaterial.cpp:319-328` | CRITICAL | non-single-valued central deposit (BUG-3) |
| 4 | `spatial_dyn_driver.cpp:811-816` | CRITICAL | residual matrix+RK guard the plan misses (BUG-4) |
| 5 | `bimaterial_wave_operator.hpp:194`, `inl:142-163,381-419` | MEDIUM | `sf` vs `mesh_face` key mismatch (BUG-5) |
| 6 | `bimaterial_wave_operator.inl:162,165-187,405-409` | MEDIUM | local-side stub + warn, not abort (BUG-6) |
| 7 | plan Phase 1 AC-2 | MEDIUM | single-valuedness test is trivial (BUG-7) |
| 8 | `makefile:311-317,1638-1648,859-893` | LOW | registration understated (BUG-8) |
| 9 | plan Phase 1 step 3 vs `godunov_flux_bimaterial.cpp:282-295` | LOW | scratch/sizing/acoustic under-specified (BUG-9) |

*(Rev-1 "What checked out", completeness, scope, documentation, test status, and
recommended-next-steps sections are preserved verbatim above the Rev-2 report by reference;
they are not repeated here to keep the append-only log readable. Nothing in them is
retracted.)*

---

## Code Review Report (Rev 2)

**Date:** 2026-06-04
**Scope:** revised plan `PLAN_mixed_flux_hetero_riemann.md` (now 7 phases) vs the SAME
`miniapps/seas/` source (byte-for-byte unchanged since Rev 1 — `git status`: only the
markdown docs + the moved PDF), PLUS the new Phase-7 artifacts:
`tpv31/configs/tpv31.toml`, `tpv31/visualize_results.py`, `tpv31/benchmark_data/`,
`jobs/tpv104_spatial/tpv104_p2_rk45_mixedflux_12N_600r_normal.sbatch`,
`jobs/tpv31_spatial/tpv31_p2_aderO3_noflux_50m_normal.sbatch`, and the tpv104/tpv31
config+sbatch corpus.
**Reviewer:** chunhui-code-reviewer agent (continuation pre-implementation plan audit).

### Summary
The Rev-2 plan correctly folds in **8 of the 9** Rev-1 fixes — BUG-1, BUG-2, BUG-3, BUG-4,
BUG-5, BUG-7, BUG-8, BUG-9 are all resolved, with every cited file:line re-verified
accurate against the (unchanged) live source. In particular the two Rev-1 CRITICALs are
genuinely fixed: **BUG-3** is now single-valued under heterogeneity by construction (one
side-symmetric matrix pair per face + `F_h_e2 ← F_h_e1`, mirroring the scalar
`wave_operator.inl:1203-1209` reference — no swap path remains, and the shared single-`F_h`
path consumes a strict subset of the stored pair), and **BUG-4** relaxes the exact
`spatial_dyn_driver.cpp:809-816` guard with logic consistent with the Phase-4 `ComputeMaxDt`
integrator guard. The BUG-2 CFL refactor's per-element `dt_e` placement is proven not to
perturb the `MPI_Allreduce(MIN)` semantics (uniform positive scaling commutes with both
the per-element `min` and the cross-rank `MPI_MIN`).

**However, the Rev-2 edits introduce / leave standing two CRITICAL defects, both in the
NEW material:**
1. **BUG-10 (the BUG-6 "fix" is unimplementable).** The Rev-2 lateral-heterogeneity guard
   aborts "if neighbour material differs from local by more than a tolerance." But
   `ExchangeBiMaterialNeighbours_` stores the LOCAL element's material as the neighbour's
   (`bimaterial_wave_operator.inl:162`), so `neighbour == local` by construction and the
   abort can NEVER fire — least of all in the genuinely-laterally-varying case it is meant
   to catch. The plan replaced a guard that "never fires because the entry is never
   missing" (Rev-1) with one that "never fires because neighbour always equals local"
   (Rev-2): the dangerous case is exactly the one the guard cannot detect. TPV31 itself is
   safe (depth-only / seam-continuous), so the deliverable is unaffected, but the plan's
   general safety claim is false and Test 3.4's abort assertion cannot pass against a
   faithful implementation. A correct detector must inspect the material field's lateral
   variation directly (or implement the real cross-rank exchange) — not compare the stubbed
   neighbour against local.
2. **BUG-11 (Phase-7 quantitative gate cannot run).** Phase 7 req 5 + acceptance criterion
   require `visualize_results.py --tol-rms/--tol-peak` to gate the SeisSol overlay within a
   ~5–10% peak band. Those flags do not exist — the script is a plot-only overlay with no
   error metric or pass/fail exit. The whole quantitative-validation criterion is
   unsatisfiable with the committed tooling (same class of omission as Rev-1 BUG-4: the plan
   asserts a mechanism the code does not provide).

Three further NEW findings are non-blocking but should be corrected: the benchmark README
contradicts the committed reference data (BUG-12, MEDIUM — data IS parseable, docs are
stale), `--verify-dispatch` prints only the interior-flux line while the plan attributes
the time-integrator + mixed-flux lines to it (BUG-13, MEDIUM — they come from the
unconditional banner; the grep gate is still satisfiable with tolerant patterns), and the
`--deriv-cache` driver comment at `:1146` is misquoted as "scalar/homogeneous-only" when it
actually says "works for both the scalar and bimaterial (matrix) operators" (BUG-14, LOW —
the OMIT recommendation is harmless but mis-justified).

### Previously Reported Bugs — Verification (Rev 2)

| BUG | Rev-2 resolution | Verified in code (file:line) | Verdict |
|-----|------------------|------------------------------|---------|
| BUG-1 | replace `inl:567-576` body + 4 stale comments | abort `inl:569-574`; `inl:518-520`; `hpp:121`; `wave_operator.inl:1623-1626` | FIXED ✓ |
| BUG-2 | delete `inl:521-531`; factor at `dt_e` `inl:546` | switch `inl:521-531`; `dt_e` `inl:546`; Allreduce `inl:553-555` (MIN semantics preserved) | FIXED ✓ |
| BUG-3 | side-symmetric single pair + `F_h_e2←F_h_e1` | unswapped operands `inl:471-484`; scalar ref `wave_operator.inl:1203-1209`; apply `godunov_flux_bimaterial.cpp:319-328` | FIXED ✓ |
| BUG-4 | relax `:809-816`; reject matrix+ADER+mixed!=none | guard `spatial_dyn_driver.cpp:809-816`; comment `:806-808`; `ComputeMaxDt` cross-ref Phase 4 | FIXED ✓ |
| BUG-5 | `sf`-loop shared arm | `hpp:194` (sf-keyed); upwind Pass 2 `inl:381-419`; `GetSharedFace` `inl:387`; verify `inl:405-409` | FIXED ✓ |
| BUG-6 | mismatch-predicate abort | `inl:162` stores local-as-nbr ⇒ predicate dead | **NOT FIXED → BUG-10** |
| BUG-7 | het operator-level Test 3.2; drop trivial primitive | Phase 1 AC note; Phase 3 Test 3.2 | FIXED ✓ |
| BUG-8 | 4-edit recipe + obj lists | `makefile:311-317,859-893,1638,5248,5252` | FIXED ✓ |
| BUG-9 | tmp1/tmp2 + SetSize output; acoustic verify | `godunov_flux_bimaterial.cpp:282-295`; `Central` `godunov_flux.cpp:416` | FIXED ✓ |

### Phase-7 factual verification (Rev 2)

| Claim | Evidence | Verdict |
|-------|----------|---------|
| `tpv31.toml` matrix/none/order1/depth_profile_1d/15s/cfl0.5/dg | lines 240,227,187,97,243,228,231 | TRUE |
| `tpv31_rk_mixedflux.toml` diff = only mixed_flux none→adjacent | new file (not yet created); spec is consistent with `tpv31.toml` content | OK (to-create) |
| template sbatch exist; module/LD block byte-identical | TPV104 `:71-89`, TPV31 ADER `:60-78` — identical 7-module list + `export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${LD_LIBRARY_PATH:-}"` | TRUE |
| CLI flags `--mixed-flux/--time-integrator/--verify-dispatch/--print-derived/--paraview-fault-vtu/--checkpoint-every/--restart/--output-dir/--tfinal` all exist | grep: 1 each | TRUE (all exist) |
| `--verify-dispatch` prints interior/mixed/time-integrator | `:2143` interior only; banner `:853,859` for the other two | **PARTIAL (BUG-13)** |
| compare vs `scec_seisol/`; format per README | 30 `.txt` files parseable by `:158-192,575`; README stale/contradictory | data OK / **README wrong (BUG-12)** |
| `visualize_results.py --tol-rms/--tol-peak` | grep: absent | **FALSE (BUG-11)** |
| TPV31+RK45+adjacent is first matrix+central job | tpv104 mixedflux all `interior_flux="scalar"`; tpv31 all `--mixed-flux none` | TRUE |
| deriv-cache `:1146` "safe on scalar/homogeneous" | actual: "works for both scalar and bimaterial (matrix)"; no rejection at `:1149-1158` | **MISQUOTE (BUG-14)** |
| Phase 7 = artifacts only; honors Frontera/no-local/RK-CFL | `tpv31/configs/*.toml` + `jobs/tpv31_spatial/*.sbatch`; no `--dt`; staged-not-submitted | TRUE |

### Completeness assessment (Rev 2)
- **Must change but now under-specified / wrong:** the lateral-het detector (BUG-10 — the
  prescribed predicate is dead); the Phase-7 quantitative gate tooling (BUG-11 — flags do
  not exist). Both are in the NEW material.
- **Phases 1–6 are spec-complete and correctly cited** — all BUG-1..BUG-5,7,8,9 fixes
  verified. The byte-exact regression contract is respected: the ctor-operand refactor is
  parity-guarded (22/22), the CFL extraction is byte-exact for the scalar operator, and the
  side-symmetric central store touches only previously-aborted code.
- **Do the proposed tests catch the highest risks?** BUG-3 → Test 3.2 (het, operator-level)
  ✓. BUG-4 → Test 5.1 + Frontera gate ✓. BUG-2 → Test 4.1 ✓. BUG-10 → Test 3.4b CANNOT
  pass as specified (must be rewritten against a real detector). BUG-11 → no test exists;
  proposed NEW-7.1.

### Scope / project-rule compliance (Rev 2)
- Phases 1–6 stay within `dynamic/` + `drivers/spatial_dyn_driver.cpp` + `tests/unit/`
  (+ `tests/parallel/`) + `Makefile`. Phase 7 adds `tpv31/configs/*.toml` +
  `jobs/tpv31_spatial/*.sbatch` (benchmark artifacts, NOT under the no-touch source tree)
  — **compliant**. The BUG-11/BUG-12 fixes touch `tpv31/visualize_results.py` and
  `tpv31/benchmark_data/README.md`, also benchmark artifacts (not source) — in scope.
  **No BP5/friction/solver-core edit is required or proposed** (project memory [C2]/no-touch
  honored).
- **R-003 (structural exclusion)** is lifted as the feature's purpose, correctly cited; the
  **R-004 stub** limitation (BUG-6/BUG-10) is the one place the plan still over-claims. No
  prior fix is silently reverted.
- **Frontera-approval / no-local-production-run / explicit-RK-CFL** project rules: Phase 7
  stages (does not submit) the job, uses tiny fixtures locally, omits `--dt` (auto-CFL with
  the RK-aware 0.6 Adjacent factor + `cfl=0.5` + `cfl_safety="dg"`). **All honored.**

### Documentation Status (Rev 2)
- Plan-internal: the Rev-2 plan correctly enumerates all four stale source comments to
  update (BUG-1) and adds the central-formula derivation requirement (Phase 1).
- Artifact docs needing correction: `tpv31/benchmark_data/README.md` (BUG-12 — stale
  "directory empty" + `.dat`/header description vs the 30 committed `.txt` files).
- The bi-material central derivation block (why `½A_self·Q_self + ½A_nbr·Q_nbr`, NOT
  `½A_self·(Q_self+Q_nbr)`) is correctly required as a header comment (Phase 1).

### Test Status (Rev 2)
- **Implemented:** 0 of 15 proposed (pre-implementation; existing
  `test_bimaterial_wave_operator_parity` 22/22 + scalar mixed-flux dispatch tests remain
  the regression floor).
- **Test framework:** custom `TEST_ASSERT` harness + `make test` group (no gtest);
  `tpv31/visualize_results.py` tested separately if extended (BUG-11).
- **Priority tests for merge readiness:** Test 3.2 (BUG-3 het single-valuedness),
  Test 5.1 (BUG-4 guard), Test 4.1 (BUG-2 CFL parity), Test 6.1 (homog gate),
  Test 3.4b rewrite (BUG-10), Test NEW-7.1 (BUG-11 tolerance gate).

### Recommended Next Steps (Rev 2)
1. **BUG-10 (CRITICAL):** Replace the dead `neighbour vs local` abort predicate with a real
   lateral-variation detector (inspect `material_->EvalAt` near the seam, OR a
   `material.seam_continuous` declaration, OR implement the real `MPI_Allgatherv` exchange).
   Rewrite Test 3.4b to assert the abort via that detector. Do not ship a guard that cannot
   fire on the case it names.
2. **BUG-11 (CRITICAL):** Either add `--tol-rms`/`--tol-peak` (with RMS/peak error compute
   and non-zero exit) to `tpv31/visualize_results.py`, or reword Phase 7 req 5 + the
   acceptance criterion to a qualitative visual check and drop the non-existent flags.
   Add Test NEW-7.1.
3. **BUG-12 (MEDIUM):** Rewrite `benchmark_data/README.md` to match the 30 committed `.txt`
   files + the script's actual parser contract; have Phase 7 cite the script, not the README.
4. **BUG-13 (MEDIUM):** Reword Phase 7 req 3/4 to grep banner lines (`time integrator:`,
   `mixed flux:`) + the `--verify-dispatch` `interior flux :` line with tolerant patterns,
   mirroring the proven TPV104/TPV31 sbatch grep; fix the attribution and exact-string quoting.
5. **BUG-14 (LOW):** Correct the `--deriv-cache` quote ("works for both operators;
   geometry-only; R-002 round-off; R-004 budget"); keep OMIT as a correctly-justified
   first-run choice.
6. Phases 1–6 are ready to implement as specified; proceed once BUG-10/BUG-11 are addressed
   in the plan, and re-review (Rev 3) is optional for the Phase-7 tooling changes only.

---

## Code Review Report (Rev 3)

**Date:** 2026-06-04
**Scope:** plan `PLAN_mixed_flux_hetero_riemann.md` **Rev 3** vs the SAME `miniapps/seas/`
source (byte-for-byte unchanged since Rev 1/2 — `git status`: only the markdown docs + the
moved PDF `miniapps/seas/document/mixed_flux_dev/PLAN_mixed_flux_hetero_riemann.pdf` are
new/modified; `reports/.report-config` touched). Fix report read:
`PLAN_mixed_flux_hetero_riemann_fix.md` (plan-level edits, no code).
**Reviewer:** chunhui-code-reviewer agent (continuation pre-implementation plan audit; fresh
adversarial pass on the two new specifications, NOT a mere verification of the fixes).

### Summary
Plan Rev 3 set out to fix the five Rev-2 open findings. **Three are cleanly resolved**
(BUG-12 README, BUG-13 verify-dispatch attribution, BUG-14 deriv-cache misquote — all
re-verified accurate against live source). **The two CRITICALs are only partially resolved,
and the two new specifications they introduced each contain a fresh CRITICAL:**

1. **BUG-10 (the lateral-het detector).** Rev 3 correctly *removes* the dead
   `neighbour vs local` predicate (the Rev-2 root cause) and replaces it with a "material-field
   probe at `x_loc`/`x_nbr`." But the probe is **not implementable as written (BUG-15)** —
   the real `MaterialField` exposes only `EvalAt(elem, ElementTransformation& T,
   IntegrationPoint& ip, …)` (`heterogeneous_material.hpp:183-209`), NOT an arbitrary-point
   `EvalAt(x)`; the operator holds only `const MaterialField *material_` (`hpp:182`) with no
   handle to the wrapper's `eval_at_xyz`, and a point in the neighbour element across a
   partition seam has no `ElementTransformation` the local rank can `Transform` through
   (`FunctionCoefficient::Eval` resolves the physical point via `T.Transform(ip, transip)`,
   `fem/coefficient.cpp:173-179`). **Worse, even if it were reachable, the probe would
   FALSE-ABORT on TPV31 itself (BUG-16):** it probes along the face normal, so for any
   fault-adjacent tet seam with a vertical (`n_z≠0`) normal — the common case in a ParMETIS
   tet partition — `x_loc` and `x_nbr` are at different depths and TPV31's depth profile
   (`tpv31.toml:97-145`, step jumps at 2400/5000/10000 m + linear ramp) makes
   `m(x_nbr)≠m(x_loc)` by the *legitimate* depth variation, tripping the `~1e-6` abort. The
   plan's claim "the abort does NOT fire for TPV31" is false for general seams. The
   `seam_continuous=true` alternative the plan also offers is sound and suffers neither defect;
   the fix is to MANDATE it and DELETE the probe arm.

2. **BUG-11 (the Phase-7 tolerance gate).** Rev 3 correctly moves implementing
   `--tol-rms/--tol-peak` into scope (Phase 7 req 5(a) + Files-to-Modify + Test NEW-7.1), so
   the "cites a nonexistent flag" structural defect is fixed. But the *metric itself* is
   under-specified in three ways that block a correct implementation: the "relative error"
   denominator is unspecified and the obvious per-sample form divides by zero-crossing
   reference values (**BUG-17**, MEDIUM); the gated station COLUMN is never named and MFEM's
   col 8 `mu_eff` has no reference (all-NaN ⇒ silent false pass) (**BUG-18**, MEDIUM); and the
   RMS/peak difference requires resampling the MFEM (RK45 adaptive dt) and SeisSol traces onto
   a common time base, which `visualize_results.py` does NOT do (no `np.interp` anywhere) and
   the plan never mentions (**BUG-19**, HIGH).

A LOW finding (**BUG-20**) notes the probe's reflection distance/tolerance/containment are
under-defined; it is moot once the probe is dropped.

**Net:** Phases 1–6 (BUG-1..BUG-9) remain spec-complete and correctly cited; BUG-3 (single-
valued central deposit) and BUG-4 (driver guard) are re-confirmed (the Rev-3 edits were
confined to the Constraint/Phase-2/Phase-3-Test/Phase-7 text and did not touch the BUG-3
dispatch spec at `inl:471-495` or the BUG-4 guard at `:809-816` — both still match live
source). The plan is NOT yet ready to implement the parallel lateral-het guard or the
quantitative Phase-7 gate; the documentation/attribution fixes (BUG-12/13/14) are ready.

### Previously Reported Bugs — Verification (Rev 3)

| BUG | Rev-3 plan change | Verified in code (file:line) | Verdict |
|-----|-------------------|------------------------------|---------|
| BUG-10 | dead predicate removed; replaced by material-field probe OR `seam_continuous` | probe assumes `EvalAt(x)` — absent (`heterogeneous_material.hpp:183`); only `MaterialField*` held (`hpp:182`); depth jumps `tpv31.toml:100-145` | **PARTIAL — root cause fixed, probe broken (→BUG-15/16/20)** |
| BUG-11 | flags moved in-scope; Test NEW-7.1 added; exit-code gate | argparse still lacks flags (`visualize_results.py:373-434`) — to be ADDED; no interp/RMS/peak code present | **PARTIAL — flag gap closed; metric incomplete (→BUG-17/18/19)** |
| BUG-12 | README rewrite to 60 `.txt`/8-col/SeisSol header; cite parser | README still stale (`.dat` `:34`, 7-col `:38`, "empty" `:53`, phantom flags `:60-61`); fix targets correct text | **RESOLVED** |
| BUG-13 | banner (`:853,:859`) + verify-dispatch (`:2130-2147`) attribution; tolerant grep | live `:853` `time integrator:`, `:859` `mixed flux:`, verify-dispatch `interior flux :` (`:2130-2147`) | **RESOLVED** |
| BUG-14 | quote real `:1140-1148` ("works for BOTH operators … geometry-only"); OMIT re-justified | live comment confirms "both scalar and bimaterial (matrix)"; no rejection at `:1149-1158` | **RESOLVED** |
| BUG-3 | (no Rev-3 edit) | dispatch `inl:471-495`, side-symmetric store spec intact; scalar ref `wave_operator.inl:1203-1209` | **STILL FIXED** |
| BUG-4 | (no Rev-3 edit) | guard `spatial_dyn_driver.cpp:809-816` unchanged; relaxation spec consistent w/ Phase-4 | **STILL FIXED** |

### New findings (Rev 3)

| ID | Severity | File:line | One-line |
|----|----------|-----------|----------|
| BUG-15 | CRITICAL | `heterogeneous_material.hpp:183-209`, `bimaterial_wave_operator.hpp:182`, `.inl:46,64-70`, `fem/coefficient.cpp:173-179` | probe needs arbitrary-point `EvalAt(x)`; real API has only `(elem,T,ip)`; no `eval_at_xyz` on `MaterialField`; x_nbr in neighbour element unreachable |
| BUG-16 | CRITICAL | `tpv31/configs/tpv31.toml:97-145`, `heterogeneous_material.cpp:240-295` | along-normal probe conflates depth vs lateral; FALSE-ABORTS on TPV31 depth profile across `n_z≠0` seams |
| BUG-17 | MEDIUM | plan Phase 7 req 5(a); `visualize_results.py:158-192` | "relative error" denominator unspecified; per-sample form divides by zero-crossing refs |
| BUG-18 | MEDIUM | plan Phase 7 req 5(a); `visualize_results.py:125-155,191,205-214` | gated COLUMN unnamed; MFEM col 8 `mu_eff` has all-NaN reference ⇒ silent false pass |
| BUG-19 | HIGH | plan Phase 7 req 5(a); `visualize_results.py` (no `np.interp`) | RMS/peak needs resampling MFEM RK45 vs SeisSol onto common time base; absent in script & plan |
| BUG-20 | LOW | plan Phase 2 req 4; `bimaterial_wave_operator.inl:330-342` | probe δ/containment/`1e-6` tol under-defined for general tet seams (moot if probe dropped) |

### Internal-consistency check (Rev 3)
The Rev-3 fix report claims plan self-consistency was grep-verified. The DEAD-PREDICATE wording
is indeed gone everywhere (Constraint, Phase 2 req 4, Edge Cases, Risk Assessment, Tricky-areas,
Test 3.4b all now describe the probe / `seam_continuous`). **However the SAME detector mechanism
is described INCONSISTENTLY across those sites in a way that masks BUG-15/16:** the Constraint
and Risk Assessment say "evaluable at any point on any rank without communication" (assumes the
nonexistent point API); Phase 2 req 4 offers BOTH the probe AND `seam_continuous` with "pick
ONE"; Test 3.4b says it "MUST exercise the *real* detector (the material-field probe /
`seam_continuous` flag)" — i.e. it permits the broken probe arm. These must be reconciled onto
the single implementable mechanism (`seam_continuous`). The Phase-7 vs Phase-2 agreement on the
gate strings (BUG-13) and the deriv-cache quote (BUG-14) ARE internally consistent. No NEW
cross-section contradiction was introduced beyond the probe-implementability one.

### Scope / project-rule compliance (Rev 3)
- All Rev-3 changes remain plan-document edits; no source written (honors the pre-implementation
  status). The proposed fixes (mandate `seam_continuous`; add interp/peak-norm metric) stay
  within `dynamic/` + `tpv31/` benchmark artifacts — NOT the BP5/friction/solver no-touch tree
  (project memory [C2]). The `seam_continuous` flag is a config/material addition in
  `dynamic/heterogeneous_material.*` + the config parser, both in-scope.
- R-003 (structural exclusion lift), R-004 (stub limitation), R-002 (deriv-cache round-off):
  citations are correct in Rev 3. No prior fix reverted.
- Frontera-approval / no-local-production-run / explicit-RK-CFL: still honored (Phase 7 stages,
  omits `--dt`, tiny local fixtures only).

### Documentation Status (Rev 3)
- `tpv31/benchmark_data/README.md` rewrite (BUG-12) is correctly specified and RESOLVED at the
  plan level; the README is still stale in the tree (to be fixed at implementation).
- The Phase-1 central-derivation header comment requirement is intact.
- NEW doc gap: if `seam_continuous` is adopted, it must be documented in the config schema
  (`safs/.../spatial_friction_config_schema.md`) and set in `tpv31.toml` — fold into the
  BUG-15 fix.

### Test Status (Rev 3)
- **Implemented:** 0 of 19 proposed (pre-implementation; `test_bimaterial_wave_operator_parity`
  22/22 remains the regression floor).
- **Test framework:** custom `TEST_ASSERT` + `make test` (C++); `pytest` for the Python
  `visualize_results.py` tooling tests (NEW-7.1/7.4).
- **Priority tests for merge readiness:** Test 3.2 (BUG-3), Test 5.1 (BUG-4), Test 4.1 (BUG-2),
  Test 6.1 (homog gate), **Test 3.4b rewritten against `seam_continuous` (BUG-10/15)**,
  **Test 3.4c depth-profile no-false-abort (BUG-16)**, **Test NEW-7.1 + NEW-7.4 with
  peak-normalized interpolated metric (BUG-11/17/18/19)**.

### Recommended Next Steps (Rev 3)
1. **BUG-15 + BUG-16 (CRITICAL):** In Phase 2 req 4 (and Constraint / Edge Cases / Risk
   Assessment / Tricky-areas / Test 3.4b), DELETE the material-field probe and mandate
   `material.seam_continuous=true` as the sole detector for `Mode::Coefficient` central shared
   faces; keep the `MPI_Allgatherv` exchange as the deferred fully-correct path. Add the flag to
   the config schema + `tpv31.toml`. Correct the Phase-7 Edge-Case claim about "the abort does
   not fire for TPV31." Add Test 3.4c.
2. **BUG-19 (HIGH) + BUG-17/BUG-18 (MEDIUM):** In Phase 7 req 5(a), specify the metric:
   interpolate the named MFEM channel(s) (e.g. `V_strike`, `slip_strike`, `tau_strike`;
   exclude `mu_eff`) onto the reference `time_s` over the overlap window via `np.interp`, then
   compute peak-normalized (`/max|ref|`, zero-safe) RMS and peak relative error; fail the
   station if any gated channel exceeds either tol. Update Test NEW-7.1 to use offset time grids
   + a zero-crossing reference; add Test NEW-7.4.
3. **BUG-20 (LOW):** Drops out once the probe is removed (step 1).
4. BUG-12/13/14 are RESOLVED at the plan level — ready to implement.
5. Phases 1–6 remain ready to implement as specified. Re-review (Rev 4) recommended only for
   the Phase-2 detector rewrite and the Phase-7 metric spec — the two areas still open.

---

## Code Review Report (Rev 4)

**Date:** 2026-06-04
**Scope:** plan `PLAN_mixed_flux_hetero_riemann.md` **Rev 4** vs the SAME `miniapps/seas/`
source (byte-for-byte unchanged since Rev 1/2/3 — `git status`: only the markdown docs + the
moved PDF `miniapps/seas/document/mixed_flux_dev/PLAN_mixed_flux_hetero_riemann.pdf` +
`reports/.report-config`). Fix report read: `PLAN_mixed_flux_hetero_riemann_fix.md` (Round 2,
plan-level edits only, no code).
**Reviewer:** chunhui-code-reviewer agent (continuation pre-implementation plan audit; FRESH
adversarial pass on the Rev-4 `seam_continuous` mechanism + the BUG-17/18/19 metric, NOT a
rubber-stamp — the author's last TWO critical fixes each contained criticals).

### Summary
Plan Rev 4 set out to fix the six Rev-3 open findings. **All six are genuinely RESOLVED:**
- **BUG-15 / BUG-16 (CRITICAL):** the material-field probe is DELETED everywhere — every
  remaining "probe" mention is the rejected option (grep-confirmed across the plan), and the
  sole guard is now the declarative `material.seam_continuous` flag, which needs no
  arbitrary-point `EvalAt(x)`, no `eval_at_xyz`, no cross-rank comm (plan `:289`), and does not
  inspect material at all ⇒ cannot false-abort on TPV31's depth profile. The BUG-15
  unimplementable-API and BUG-16 depth/lateral-conflation defects are structurally gone.
- **BUG-17 / BUG-18 / BUG-19:** the metric is fully specified — `np.interp` onto reference
  `time_s` over `[t0,t1]` (BUG-19), zero-safe `peak_rel`/`rms_rel` with denominator floor
  (BUG-17), named channels `V_strike`/`slip_strike`/`tau_strike`/`sigma_n` excluding `mu_eff`
  (BUG-18). **The channel names were VERIFIED against the live script** — they are exactly the
  dict keys in `visualize_results.py:145-155` (`load_mfem_file`), `:179-192`
  (`load_reference_file`), and the `PANELS` list `:205-214`; `mu_eff` is MFEM col 8 (`:137,154`)
  with an all-NaN reference (`:191`), so the exclusion is correct. **No channel-name mismatch.**
- **BUG-20 (LOW):** dissolved with the probe.
- **BUG-10 / BUG-11 (the two long-running CRITICALs) are now RESOLVED**, and **BUG-6** is
  resolved by supersession.

**However — and this is the lead finding — the Rev-4 fix introduced / left standing THREE NEW
defects, one a hard CRITICAL functional blocker, exactly the pattern of the prior two rounds:**

1. **BUG-21 (CRITICAL) — the `seam_continuous`/`mixed_flux="adjacent"` TOML co-trips a SECOND,
   config-PARSE-level mutual-exclusion guard the plan never mentions.**
   `spatial/code/spatial_friction.cpp:1162-1167` MFEM_ABORTs `interior_flux="matrix"` +
   `mixed_flux != "none"` inside `LoadSpatialFrictionConfig`, which the driver calls at
   `spatial_dyn_driver.cpp:616` — UPSTREAM of the CLI override (`:679`), the BUG-4 driver guard
   (`:809-816`), the operator ctor (`:1135`), and `SetMixedFluxMode` (`:1245`). The Phase-7
   deliverable TOML sets `mixed_flux="adjacent"` (`plan:606-608,754-755`), so it aborts at
   config load. The plan relaxes only ONE of TWO series guards. This is the THIRD consecutive
   round in which the author's fix exposed an unaddressed sibling guard (Rev-1 BUG-4 driver
   guard; Rev-3 BUG-15/16 probe; now the config-parse guard).
2. **BUG-22 (HIGH/SCOPE) — the new flag's home is outside the plan's declared scope.** The
   `MaterialSpec` struct (`spatial_friction.hpp:573`), its parse (`spatial_friction.cpp:1626+`),
   and the BUG-21 guard ALL live in `spatial/code/spatial_friction.{hpp,cpp}`, which is NOT in
   the plan's scope list (`dynamic/` + `drivers/spatial_dyn_driver.cpp` + `tests/` + `Makefile`,
   plan `:71-75`). The plan's hand-wave "the config parse lives with the driver" (`:613`) is
   factually wrong and silently expands the no-touch boundary. (Project memory [C2] does NOT
   forbid editing `spatial/code/`, so this is a plan-scope-accuracy + regression-assurance
   defect, not a hard policy violation — but it must be owned explicitly and parser-regression
   guarded, since every TPV*/BP5/TPV31-ADER config parses through this file.)
3. **BUG-23 (MEDIUM) — the metric's overlap window `[t0,t1]` is undefined and a truncated MFEM
   run is unhandled.** `np.interp` clamps out-of-range to endpoints, so a run that doesn't reach
   `tfinal=15s` (48 h wall / failed restart) either spuriously fails (full-grid clamp of the
   MFEM tail) or silently passes (tiny overlap). No minimum-coverage guard.

A MODERATE assurance gap (**BUG-24**): the plan does not state explicitly that the
`seam_continuous` guard leaves the existing TPV31 ADER/`mixed_flux=none` job byte-exact. It
DOES, structurally — `BuildCentralFluxFaceSet_` returns with an empty
`central_flux_face_set_` on `mixed_flux_mode_==None` (`wave_operator.inl:1653-1654`), so the
guard (which only walks the central set / central-shared intersection) is never reached on a
`none` run; and BUG-21's config-parse guard allows `matrix + none`. So **Part B.2's regression
concern is correctly handled by construction** — but should be written down.

**Net verdict:** BUG-15/16 and the BUG-11-family (BUG-11/17/18/19) are GENUINELY fixed; the
`seam_continuous` flag is the right mechanism and does NOT introduce a runtime regression on
the existing ADER job (Part B.2 — clean by construction) and has NO channel-name mismatch
(Part C.1 — verified). BUT it DOES introduce a config-parse blocker (BUG-21, Part B.1) and a
scope expansion (BUG-22, Part B.1) that the plan must own before implementation. The plan is
NOT yet clean enough to implement the Phase-5/Phase-7 config path; Phases 1–4 + 6 (the C++
operator/CFL/test machinery) ARE clean and ready.

### Part-by-part answers to the review brief

- **Part A (BUG-15..20 resolved? probe gone?):** YES to all six. The probe is gone as a
  mechanism (only the rejected option). Internally consistent across Constraint / Phase 2 /
  Phase 3 Tests / Phase 7 / Testing Strategy / Risk / Tricky-areas — all describe the SAME
  `seam_continuous` flag and the SAME interp+peak-norm metric. No "pick ONE" / `EvalAt(x)` /
  per-sample-metric residue remains.
- **Part B.1 (config key scope + parse site):** The `[material]` block + `interior_flux` +
  `mixed_flux` are parsed in `spatial/code/spatial_friction.{hpp,cpp}` — OUTSIDE the plan's
  declared scope (**BUG-22**) — and a parse-level guard there aborts the feature config
  (**BUG-21**). Both are NEW and load-bearing.
- **Part B.2 (regression on the existing TPV31 ADER job):** NO regression. The guard fires only
  inside `BuildPerFaceCentralFluxMatrices_` for central-set faces with `mf_on_` true; the ADER
  job (`mixed_flux=none`) has an empty central set (`wave_operator.inl:1653-1654`) and never
  reaches the guard. Plan should STATE this (**BUG-24**, moderate) — but it is correct.
- **Part B.3 (ordering):** Adequately specified — the plan requires `SetSeamContinuous` (or the
  ctor arg) BEFORE `SetMixedFluxMode` (`plan:610-611`); since `BuildPerFaceCentralFluxMatrices_`
  runs inside `SetMixedFluxMode` (driver `:1245`, after the ctor `:1135`), a setter placed in
  `[1135,1245)` is safe and the ctor-arg form is automatically safe. Minor: pin "before
  `SetMixedFluxMode`" → "before the member is read by `BuildPerFaceCentralFluxMatrices_`."
- **Part B.4 (semantics):** `[material].seam_continuous` is the right granularity (it qualifies
  the material's lateral continuity, a material property). It correctly applies only to
  `Mode::Coefficient` central SHARED faces — `Mode::Constant` is trivially seam-continuous and
  the matrix path forbids Constant anyway (driver `:1123`), so the guard never considers it.
  Default `false` (fail-loud) is safe for all existing configs (which omit it) and forces the
  TPV31 config to affirm it. Sound.
- **Part C.1 (channel names match the script):** YES — `V_strike`/`slip_strike`/`tau_strike`/
  `sigma_n` are the exact dict keys (`visualize_results.py:145-155,179-192,205-214`).
- **Part C.2 (`mu_eff` exclusion):** Correct — MFEM col 8 (`:137,154`), reference all-NaN
  (`:191`); the plan's exclusion instruction is unambiguous.
- **Part C.3 (window / partial coverage):** Under-specified → **BUG-23**.
- **Part D (BUG-3/BUG-4 + Phases 1–6 fresh pass):** BUG-3 dispatch (`inl:471-484`, unswapped
  operands) and BUG-4 driver guard (`:809-816`) re-confirmed STILL correct in live source. The
  CFL switches (`inl:521-531` divergent; `wave_operator.inl:5880-5903` `cfl_rk_aware_`-gated),
  the abort body (`inl:567-576`), `ApplyPerFaceFlux` (`cpp:301-329`, `F=M0·Q_self+M1·Q_nbr`),
  and the ctor build order (`ExchangeBiMaterialNeighbours_` `:49`, upwind build `:248-419`) all
  match the plan's citations. No new Phase-1..6 C++ defect found.

### Previously Reported Bugs — Verification (Rev 4)

| BUG | Rev-4 plan change | Verified in code (file:line) | Verdict |
|-----|-------------------|------------------------------|---------|
| BUG-15 | probe DELETED; declarative `seam_continuous` flag is the sole guard | probe absent as mechanism (grep); flag reads only `seam_continuous_` member; no `EvalAt(x)` needed | **RESOLVED** |
| BUG-16 | probe gone ⇒ no depth/lateral conflation; flag does not inspect material | `tpv31.toml:97-145` depth jumps no longer probed; Phase-7 Edge Case corrected (`plan:744-749`) | **RESOLVED** |
| BUG-17 | peak-norm `peak_rel`/`rms_rel` + floor (`plan:718-722`) | n/a (plan spec); `visualize_results.py:158-192` reader unchanged | **RESOLVED** |
| BUG-18 | named channels, exclude `mu_eff` (`plan:714-717`) | names == dict keys `visualize_results.py:145-155,179-192,205-214`; `mu_eff` all-NaN `:191` | **RESOLVED** |
| BUG-19 | `np.interp` onto ref `time_s` over `[t0,t1]` before diff (`plan:708-712`) | `visualize_results.py` still has no interp (to be added); offset-grid test req'd | **RESOLVED** (residual → BUG-23) |
| BUG-20 | dissolved with the probe | n/a | **RESOLVED** |
| BUG-10 | dead/probe → declarative `seam_continuous` flag | `inl:162` stub unchanged (flag bypasses it); `MPI_Allgatherv` deferred | **RESOLVED** |
| BUG-11 | metric fully specified + in-scope tooling | `visualize_results.py:373-434` argparse still lacks flags (to add); metric now complete | **RESOLVED** (residual window → BUG-23) |
| BUG-6 | enforcement = `seam_continuous` flag | superseded; regime explicit + fail-loud | **RESOLVED** (by supersession) |
| BUG-3 | (no Rev-4 edit) | dispatch `inl:471-484` unswapped; side-symmetric store spec intact | **STILL FIXED** |
| BUG-4 | (no Rev-4 edit) | driver guard `spatial_dyn_driver.cpp:809-816` unchanged; relaxation spec intact | **STILL FIXED** (but see BUG-21: a 2nd guard exists) |

### New findings (Rev 4)

| ID | Severity | File:line | One-line |
|----|----------|-----------|----------|
| BUG-21 | CRITICAL | `spatial/code/spatial_friction.cpp:1162-1167`; `drivers/spatial_dyn_driver.cpp:616,679,809-816,1245` | config-parse guard aborts `matrix + mixed_flux!="none"` at `LoadSpatialFrictionConfig`, before every site the plan touches; Phase-7 TOML sets `mixed_flux="adjacent"` ⇒ aborts at load. Plan relaxes only the driver guard, not this one |
| BUG-22 | HIGH | `spatial/code/spatial_friction.hpp:573`, `.cpp:1626+`; `makefile:132-133,1304-1305`; plan `:71-75,610-614` | `seam_continuous` field + parse + the BUG-21 guard live in `spatial/code/` — outside the plan's declared scope; plan's "config parse lives with the driver" is false |
| BUG-23 | MEDIUM | plan Phase 7 req 5(a) `:708-712`; `tpv31/visualize_results.py:146,180` (`np.interp` clamps) | metric overlap window `[t0,t1]` undefined + truncated-MFEM run unhandled ⇒ spurious fail or silent pass; no min-coverage guard |
| BUG-24 | MODERATE | plan `:40-45,272-296`; `dynamic/wave_operator.inl:1653-1654` | plan does not STATE the guard is `mf_on_`/central-set-gated, so the existing TPV31 ADER/`none` job is byte-exact (true by construction, just unwritten) |

### Internal-consistency check (Rev 4)
The probe is gone everywhere as a mechanism (grep over the whole plan: every "probe" token is
in a "NOT a probe" / "unimplementable" / "would wrongly abort" context). The Constraint,
Phase 2 req 4, Phase 2 Edge Cases, Phase 3 Tests 3.4a/b/c, Phase 7 (config + req 5(a) + Edge
Cases + Acceptance), Testing Strategy, Risk Assessment, and Tricky-areas ALL describe the SAME
`seam_continuous` flag and the SAME interp+peak-norm+named-channel metric — no section still
implies a probe or an unnamed/per-sample metric. The Rev-3 "pick ONE" contradiction is gone.
**The only internal-consistency defect is the one BUG-21/BUG-22 expose:** the plan's scope
Constraint and its Phase-7 wiring note disagree with reality about WHERE the config lives, and
the plan is silent on the second (config-parse) guard. No other cross-section contradiction was
introduced.

### Scope / project-rule compliance (Rev 4)
- All Rev-4 changes are plan-document edits; no source written (pre-implementation status
  honored — `git status` confirms only markdown + PDF + `reports/.report-config`).
- **Scope accuracy defect (BUG-22):** the `seam_continuous` implementation requires editing
  `spatial/code/spatial_friction.{hpp,cpp}`, which the plan's declared scope omits. Project
  memory [C2] no-touch (`bp5/bp1/bp2/domain/fault/solver`, `friction/dieterich_ruina.hpp`) does
  NOT list `spatial/code/`, so the edit is not policy-forbidden — but the plan must own it and
  guard the parser regression (`seas_test_spatial_friction_config`, `test_tpv_config_parse`).
- R-003 (structural exclusion lift), R-004 (stub limitation), R-002 (deriv-cache round-off),
  R-1203 (the config-parse mutual-exclusion BUG-21 must cite when relaxing): citations are
  correct where present; R-1203 is the one the plan must ADD (BUG-21). No prior fix reverted.
- Frontera-approval / no-local-production-run / explicit-RK-CFL: still honored (Phase 7 stages,
  omits `--dt`, tiny local fixtures only; the TPV31 production run is Frontera-only and staged).

### Documentation Status (Rev 4)
- The Phase-7 README rewrite (BUG-12) and the central-derivation header comment requirement
  (Phase 1) remain correctly specified.
- NEW doc gap (fold into BUG-22): `material.seam_continuous` must be added to the config schema
  `safs/project_7.0_alternative/document/spatial_friction_config_schema.md` (where `[material]`
  / `mixed_flux` keys are documented, `:151`) and to the `tpv31_rk_mixedflux.toml` comment, with
  the supported-regime semantics (Constant / depth-only seam-continuous) and the deferred
  `MPI_Allgatherv` follow-up.

### Test Status (Rev 4)
- **Implemented:** 0 of 21 proposed (pre-implementation; `test_bimaterial_wave_operator_parity`
  22/22 + `seas_test_spatial_friction_config` + `test_tpv_config_parse` remain the regression
  floor — the latter two now also guard the BUG-22 parser change).
- **Test framework:** custom `TEST_ASSERT` + `make test` (C++); `pytest` for the Python
  `visualize_results.py` tooling tests (NEW-7.1/7.4/7.5).
- **Priority tests for merge readiness:** Test 5.3 (BUG-21/22 config-parse, NEW), Test 3.2
  (BUG-3 het single-valuedness), Test 5.1 (BUG-4 driver guard), Test 4.1 (BUG-2 CFL parity),
  Test 6.1 (homog gate), Test 3.4b (`seam_continuous` abort) + 3.4c (no-false-abort), Test
  NEW-7.1 + NEW-7.4 (interp/peak-norm metric) + NEW-7.5 (truncated-run coverage).

### Recommended Next Steps (Rev 4)
1. **BUG-21 (CRITICAL):** In Phase 5 (and Phase 7), add `spatial/code/spatial_friction.cpp:1162-1167`
   as a SECOND guard to relax — allow `interior_flux="matrix"` + `mixed_flux != "none"` at parse
   (cite R-1203); keep the matrix+ADER rejection in the DRIVER (the parser can't see the resolved
   integrator before the CLI override at `:679`). Add Test 5.3.
2. **BUG-22 (HIGH):** Amend the Constraint scope list to include
   `spatial/code/spatial_friction.{hpp,cpp}`; correct the Phase-7 wiring note (it does NOT live
   with the driver); require the `MaterialSpec::seam_continuous` field (default `false`) + parse,
   a parser-regression check, and the schema-doc + `tpv31.toml` documentation.
3. **BUG-23 (MEDIUM):** In req 5(a), define `[t0,t1]` = intersection of the two time spans, add a
   minimum-coverage guard (fail / explicit insufficient-coverage on a truncated run), and add
   Test NEW-7.5.
4. **BUG-24 (MODERATE):** Add one sentence to Phase 2 req 4 + the Constraint asserting the guard
   is `mf_on_`/central-set-gated, so `mixed_flux=none` (the existing TPV31 ADER job) stays
   byte-exact.
5. Phases 1–4 + 6 (the C++ operator / CFL / unit-test machinery) remain spec-complete and
   correctly cited — ready to implement. Phase 7's Frontera artifacts are ready EXCEPT for the
   BUG-21 config-parse relaxation and the BUG-23 metric window, both of which block a runnable
   deliverable. A check Rev 5 should verify the BUG-21/22 guard+scope edits and the BUG-23
   window spec — the only areas still open.


---

## Code Review Report (Rev 5)

**Date:** 2026-06-04
**Scope:** plan `PLAN_mixed_flux_hetero_riemann.md` **Rev 5** vs the SAME `miniapps/seas/` source
(byte-for-byte unchanged since Rev 1 — `git status`: only the markdown docs + the moved PDF +
`reports/.report-config`; `0` source files changed). This is a pre-implementation plan audit; the
deliverable is this check document, not a `REVIEW.md`.
**Reviewer:** chunhui-code-reviewer agent (FOURTH consecutive continuation review — FINAL
implementability verdict requested).

### Summary
Plan Rev 5 RESOLVES all four Rev-4 findings (BUG-21 CRITICAL, BUG-22 HIGH, BUG-23 MEDIUM, BUG-24
MODERATE), each re-verified byte-accurate against the live source. The headline deliverable of this
round — an exhaustive guard inventory claiming **exactly THREE** series guards (G1/G2/G3) block
`matrix+adjacent+RK` with four others correctly KEPT — **is INDEPENDENTLY CONFIRMED**. Because the
recurring failure mode across four rounds has been a missed sibling guard (probe → BUG-15/16; missed
G1 → BUG-21), this pass did NOT trust the inventory: it re-ran the sweep over the WHOLE relevant tree,
explicitly including the four areas the author's `grep` did not cover (the RK time-stepper, the ADER
substep path, the operator/base ctors, and the integrator-config / `fault_iterator` paths). No fourth
blocking guard surfaced. Two NEW findings emerged, both NON-blocking (BUG-25 MEDIUM, BUG-26 LOW) and
both fold into existing tests. **No CRITICAL or HIGH finding is open after Rev 5.**

### PART B — Independent guard sweep (the priority): the "exactly three guards" claim is CONFIRMED

| Guard | File:line (live) | Verdict | Disposition |
|-------|------------------|---------|-------------|
| **G1** | `spatial/code/spatial_friction.cpp:1162-1167` | CONFIRMED blocking | RELAX (plan Phase 5 req 1) — `MFEM_VERIFY(interior_flux==Scalar \|\| mixed_flux=="none")`, fires in `LoadSpatialFrictionConfig` (driver `:616`), upstream of all |
| **G2** | `drivers/spatial_dyn_driver.cpp:809-816` | CONFIRMED blocking | RELAX (plan Phase 5 req 2 / BUG-4) — `MFEM_VERIFY(!is_rk \|\| interior_flux==Scalar)`; `:806-808` stale comment |
| **G3** | `dynamic/bimaterial_wave_operator.inl:569-574` | CONFIRMED blocking | REPLACE (plan Phase 3 step 1 / BUG-1) — `MFEM_VERIFY(m == MixedFluxMode::None)` |
| keep-1 | `spatial_friction.cpp:1702-1707` | CONFIRMED correct-to-keep | matrix ⇒ non-Constant material; TPV31 is `depth_profile_1d` ⇒ passes |
| keep-2 | `spatial_dyn_driver.cpp:1080` (`:1078-1085`) | CONFIRMED correct-to-keep | scalar ⇒ Constant material; TPV31 is matrix ⇒ N/A |
| keep-3 | `wave_operator.inl:604-610` + `:1613-1616` | CONFIRMED correct-to-keep | R-1203 precomputed-flux ↔ mixed mutual exclusion; `use_precomputed_face_fluxes_` false here |
| keep-4 | central+ADER abort (Phase 4 `ComputeMaxDt` + new G2 ADER arm) | CONFIRMED intended | the integrator constraint, deliberately retained |

**Areas swept beyond the author's `drivers/`+`spatial/code/`+`dynamic/` grep (user Part B 1-4):**
1. **RK path / time-stepper.** `rk_time_stepper.{hpp,cpp}` has NO `matrix ⊥ mixed_flux` abort. Its
   header "SCALAR interior flux only" (`rk_time_stepper.hpp:24-25`) is a STALE SCOPE COMMENT, not an
   enforced guard. `AdvanceRKCoupled_Spatial`/`AdvanceRKCoupledLSW_Spatial` take
   `WaveOperator<MeshType>& wave` and call the VIRTUAL `wave.Mult` (`:256,319`), which dispatches to
   the `BimaterialWaveOperator` override (`Mult` is `virtual ... override`, `wave_operator.hpp:139`;
   the driver holds a base `unique_ptr` and constructs the bimaterial into it, `driver:1099,1135`,
   passing the base ref at `:2883`). The only RK-stepper preconditions are the friction-law verifies
   (`hpp:204` RS / `:389` LSW) — NOT flux/operator. **The matrix operator's RK/Mult path has never run
   in production (all existing matrix jobs are ADER), but there is NO guard or unimplemented branch to
   abort it** — it is a polymorphic `Mult` that will execute. The untested-combination risk is
   correctly staged as the Phase-7 Frontera dispatch gate, not a plan defect.
2. **`BimaterialWaveOperator` ctor + base ctor.** `bimaterial_wave_operator.inl:32-51`: the ONLY ctor
   `MFEM_VERIFY` is the material-mode check (`:37`, Constant/Coefficient). **NO assert couples
   mixed-flux / RK / matrix at construction.** A `SetSeamContinuous(bool)` setter (or a new ctor arg)
   inserts cleanly between the ctor (`driver:1135-1136`) and `SetMixedFluxMode` (`:1245`).
3. **Integrator config validation.** The only `time_integrator`+`interior_flux` guard is G2 (`:811`).
   The other RK-gated `MFEM_VERIFY` (`:2655-2658`) requires a `[friction.rate_state]` block for
   RK+rate_state — NOT a flux guard, and N/A to TPV31 (LSW). The RK Butcher selection is at `:2660`.
4. **Output / checkpoint / fault-iterator.** `fault_iterator="substep"` (tpv31.toml:235) gates ONLY
   the ADER branch (`AdvanceADERWithSubStep_Spatial`, `driver:2888`); on the RK path the substep
   iterator + `nuc_cb` are explicitly NOT used (`driver:2880-2881`). `FaultIteratorSupported`
   (`spatial_friction.hpp:672-675`) checks only `fault_iterator==Substep` — no flux coupling. The
   friction-iterator factory's verifies (`friction_iterator_factory.cpp:36-110`) are rate-state config
   consistency, not flux/operator.

**Conclusion (Part B):** the plan's "exactly three guards" claim is TRUE and the four kept guards are
correctly classified. No Rev-6-blocking fourth guard.

### PART A — BUG-21..24 resolution (verified)

| BUG | Rev-5 resolution | Verified in code (file:line) | Verdict |
|-----|------------------|------------------------------|---------|
| BUG-21 | 3-guard inventory; relax G1 (Phase 5 req 1, cite R-1203) | G1 `spatial_friction.cpp:1162-1167`; load@616, CLI@679, G2@811, ctor@1135, SMFM@1245 | RESOLVED |
| BUG-22 | scope += `spatial/code/`; `MaterialSpec::seam_continuous`+parse; parser-regression | `MaterialSpec` hpp:573; parse cpp:1626-1638; `toml_bool` cpp:304 (exact sig); [C2] OK | RESOLVED |
| BUG-23 | window `[t0,t1]`=intersection + interp-within + min-coverage guard + NEW-7.5 | plan `:757-769`; consistent across `:669-674,821-825,850`; `np.interp` clamp hazard handled | RESOLVED |
| BUG-24 | explicit no-regression in Constraint + Phase 2 req 4 | `BuildCentralFluxFaceSet_` empties on None `wave_operator.inl:1652-1653`; `mf_on_`@1640 | RESOLVED |

### PART C — G1 relaxation side-effects + `seam_continuous` wiring

1. **No silent-no-op landmine on G1 removal.** Grepped every matrix-path read of
   `mixed_flux_mode_`/`mf_on_`: there are exactly THREE — `SetMixedFluxMode` (G3 → replaced in
   Phase 3), `ComputeMaxDt` (the divergent local switch at `bimaterial_wave_operator.inl:518-531`,
   read at `:522` → DELETED in Phase 4/BUG-2), and the NEW Phase-3 dispatch branch. The current
   `InteriorFaceFlux_`/`SharedInteriorFaceFlux_` (`inl:471-497`) do NOT consult
   `central_flux_face_set_` at all (they unconditionally apply Godunov), and BOTH pass the UNSWAPPED
   `(Q_self, Q_nbr)` — so the Phase-3 single-valued copy design is sound and nothing downstream
   silently assumes None once `matrix+adjacent` is allowed through. CONFIRMED CLEAN.
2. **`SetSeamContinuous` ordering achievable.** The ctor (`:1135-1136`) takes no mixed-flux arg and has
   no construction-time coupling; `SetMixedFluxMode` is at `:1245`. A setter call (or a ctor arg) fits
   cleanly BETWEEN them. Adding a `SetSeamContinuous(bool)` setter does NOT conflict with the
   "interfaces that cannot change" list (which freezes `SetMixedFluxMode`/`ComputeMaxDt`/the face-flux
   hooks/`ApplyPerFaceFlux`) — it is a NEW member function, not a change to a frozen signature. ACHIEVABLE.
3. **`toml_bool` exists.** `spatial_friction.cpp:304`, `bool toml_bool(const toml::value&, const
   std::string&, bool)`. The plan's `toml_bool(m, "seam_continuous", false)` is byte-accurate. CONFIRMED.
4. **Test 5.3c feasibility — PARTIALLY accurate (→ BUG-25).** There is NO
   `SpatialFrictionConfig::operator==`/dump, so "parse byte-identically" must be realized field-by-field
   (the `test_tpv_config_parse.cpp:281-283` pattern). AND **BP5 has no `bp5/configs/*.toml`** — it never
   parses through `LoadSpatialFrictionConfig`, so it cannot be in the parse corpus. The realizable corpus
   is tpv205/102/104 (already present) + tpv31.toml (must be ADDED). This is the only inaccuracy in the
   Rev-5 edits → BUG-25 (MEDIUM, NOT a blocker).

### PART D — metric (BUG-23) + fresh full pass
- **Metric self-consistent.** The `[t0,t1]` window, interpolate-within-window, min-coverage guard, and
  peak-normalized named channels (exclude `mu_eff`) agree across req 5(a) (`:757-783`), the test entry
  (`:669-674`), and acceptance (`:821-825`). RESOLVED.
- **Probe gone.** Every "probe" mention in the plan is the explicitly-REJECTED option (`:101,103,105,294`
  — "NOT a probe", "not viable", "would false-abort"). BUG-15/16 remain RESOLVED.
- **BUG-3 / BUG-4 re-confirmed FIXED.** `InteriorFaceFlux_` (`inl:471-484`) passes unswapped operands;
  the side-symmetric store + `F_h_e2←F_h_e1` copy is single-valued by construction. The G2 guard text is
  unchanged in live source; the relaxation is logically consistent with the Phase-4 `ComputeMaxDt` ADER
  rejection.
- **No internal contradiction introduced.** The Phase-5 guard inventory, the Tricky-areas guard list
  (`plan:917-922`), and the Constraints agree; the scope list includes
  `spatial/code/spatial_friction.{hpp,cpp}` in the Constraint (`:77-87`), Phase 5 (`:521-525`), and
  Tricky-areas (`:923-925`).

### Previously Reported Bugs — Verification (Rev 5)

| BUG | Prior status | Rev-5 verification | Verdict |
|-----|--------------|--------------------|---------|
| BUG-1..9 | FIXED in Rev 2 | live source unchanged; all cites still accurate | STILL FIXED |
| BUG-10/15/16/20 | RESOLVED in Rev 4 (probe deleted) | no mandated probe remains; declarative flag only | STILL RESOLVED |
| BUG-11/17/18/19 | RESOLVED in Rev 4 (metric) | interp+peak-norm+named channels intact; window added | STILL RESOLVED |
| BUG-12/13/14 | RESOLVED in Rev 3 | unchanged | STILL RESOLVED |
| BUG-21 | OPEN (Rev 4) | G1 relaxed via 3-guard inventory | RESOLVED |
| BUG-22 | OPEN (Rev 4) | scope amended; field+parse+regression added | RESOLVED |
| BUG-23 | OPEN (Rev 4) | window + coverage guard + NEW-7.5 | RESOLVED |
| BUG-24 | OPEN (Rev 4) | explicit no-regression statement added | RESOLVED |

### Scope / project-rule compliance (Rev 5)
- The Rev-5 scope expansion to `spatial/code/spatial_friction.{hpp,cpp}` is project-policy-compatible:
  memory **[C2]** no-touch covers `bp5/bp1/bp2/domain/fault/solver` + `friction/dieterich_ruina.hpp`,
  NOT `spatial/code/`. The additive `bool seam_continuous=false` field + the purely-permissive G1
  relaxation are byte-safe for every existing config (default off; G1 only ADDS an allowed combination).
- BP5 / solver / friction-core remain untouched — confirmed BP5 does not even route through this parser.
- R-003 (structural exclusion lift), R-1203 (mixed ⊥ precomputed), R-004 (stub limitation) correctly
  cited; no prior fix silently reverted (CLAUDE.md honored).
- Frontera-approval / no-local-production-run / explicit-RK-CFL: Phase 7 stages (does not submit),
  tiny fixtures only, omits `--dt`. All honored.

### Documentation Status (Rev 5)
- Plan-internal documentation is now consistent (guard inventory ↔ Tricky-areas ↔ Constraints).
- Residual doc nits (LOW, fold into implementation): BUG-25 (Test 5.3c corpus wording; the Constraint
  sentence "EVERY TPV*/BP5/TPV31-ADER config parses through it" should drop BP5) and BUG-26 (state TPV31
  is LSW). Neither blocks implementation.

### Test Status (Rev 5)
- **Implemented:** 0 of 21 proposed (pre-implementation; `test_bimaterial_wave_operator_parity` 22/22 +
  `seas_test_spatial_friction_config` + `test_tpv_config_parse` are the regression floor).
- **Priority tests for merge readiness:** Test 5.3 (BUG-21/22, with the BUG-25 corpus correction),
  Test 3.2 (BUG-3 het single-valuedness), Test 5.1 (BUG-4 driver guard, built as LSW per BUG-26),
  Test 4.1 (BUG-2 CFL parity), Test 6.1 (homog gate), Test 3.4b/3.4c (`seam_continuous`), Test
  NEW-7.1/7.4/7.5 (metric).

### FINAL VERDICT — IMPLEMENTABLE
**YES — the plan is implementable end-to-end.** All four Rev-4 findings are resolved; the
"exactly three guards" inventory is independently confirmed complete (no fourth blocker in the RK
stepper, ADER substep, ctors, integrator-config, or fault-iterator paths); the `seam_continuous`
wiring (`toml_bool` exists, `SetSeamContinuous` ordering achievable, `MaterialSpec` field additive
+ byte-safe) is realizable; the G1 relaxation has no silent-no-op side-effect; and the SeisSol-overlay
metric is fully and self-consistently specified. **There are NO open CRITICAL or HIGH findings.**

Residual (acceptable) items, none blocking:
- **BUG-25 (MEDIUM):** correct Test 5.3c's CORPUS during implementation — drop BP5 (no spatial TOML),
  add `tpv31.toml`, realize "byte-identical" field-by-field (no `operator==` exists). This is a
  test-spec wording fix, not a design defect; the BUG-22 parser-regression intent is sound.
- **BUG-26 (LOW):** add a one-sentence note that TPV31 is LSW ⇒ the RK path is
  `AdvanceRKCoupledLSW_Spatial`; the central dispatch is law-agnostic, so no code consequence.
- **Documented residual risk (not a finding):** the matrix-operator RK/Mult path is first-exercised
  by this feature (all prior matrix jobs were ADER). There is no guard to abort it; correctness is a
  RUNTIME (CFL/stability) question the plan correctly stages as the Phase-7 Frontera dispatch gate +
  empirical CFL calibration, NOT a local pass/fail.

### Recommended Next Steps (Rev 5)
1. Proceed to implementation (`/code-implement`). Phases 1-6 are spec-complete and correctly cited;
   Phase 7 artifacts are runnable once the implemented G1 relaxation + `seam_continuous` parse land.
2. When writing Test 5.3c, apply the BUG-25 corpus correction (drop BP5; add `tpv31.toml`;
   field-by-field assertions). When writing Test 5.1 / the local smoke, build the fixture as LSW
   (BUG-26) to mirror TPV31.
3. No further plan revision is required to begin coding. BUG-25/BUG-26 can be folded into the
   implementation commits (or a trivial Rev-6 doc touch-up) — they do not gate the build.


---

# Implementation Review — Phase 1 & 2 (Rev I-1)

**Date:** 2026-06-04
**Reviewer:** chunhui-code-reviewer agent (adversarial implementation pass)
**Scope:** the WRITTEN Phase-1 + Phase-2 source diff (not the plan):
- `miniapps/seas/dynamic/godunov_flux_bimaterial.hpp` / `.cpp` — `BimaterialFlux::BuildPerFaceCentralMatricesGlobal` (new).
- `miniapps/seas/dynamic/bimaterial_wave_operator.hpp` / `.inl` — `ResolveFaceFluxOperands_`, `BuildPerFaceCentralFluxMatrices_`, `per_face_central_flux_`, `seam_continuous_` + `SetSeamContinuous`, `GetPerFaceCentralFlux`, Pass-1 refactor.
- `miniapps/seas/tests/unit/test_bimaterial_central_flux.cpp` (new); `miniapps/seas/Makefile` (new target + `test:` group entry).

This Rev introduces **implementation** finding IDs `IMPL-1..IMPL-8`, distinct from the
plan-defect history `BUG-1..BUG-26` (all of which remain RESOLVED and are left intact above).

## Builds / Tests Re-Run (independent verification — NOT trusting the "all pass" claim)

Built in the worktree against the main checkout's `libmfem.a`
(`MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN ...`), `conda -n mfem-dev`.
Symlinked the gitignored `extern/toml11` from the main checkout first (was empty;
per project memory the parity tests still build without it but the link was restored).

| Binary | Mode | Result |
|--------|------|--------|
| `seas_test_bimaterial_central_flux` | serial | **5/5 PASS** (Test 1.1 = 8.641e-16, 1.2 = 7.133e-16, 1.3 self = 4.359e-16, 1.3 nbr = 2.788e-16, 1.3 contrast = 4.049e-01) |
| `seas_test_bimaterial_wave_operator_parity` | serial | **9/9 PASS** |
| `seas_test_bimaterial_wave_operator_parity` | `mpirun -np 4` | **9/9 PASS** |
| `seas_test_phaseh_wave_operator_constant_parity` | serial | **19/19 PASS** |
| `seas_test_phaseh_wave_operator_constant_parity` | `mpirun -np 4` | **52/52 PASS** (interior 58 + shared 18 faces processed) |

All implementer claims reproduce exactly. Binary mtimes verified newer than the
modified sources, so the runs reflect the current diff.

## Claim-by-claim verification

| Claim | Verdict | Evidence |
|-------|---------|----------|
| Phase-1 central math `½·T·(AxP+AxM)·Tinv` per side, `*=0.5` after, no aliasing | **CONFIRMED** | `godunov_flux_bimaterial.cpp:363-380`; `Asum` is **reassigned** (not `+=`) between sides (`:367` vs `:375`) so no stale data; `Mult(T,Asum,tmp)` then `Mult(tmp,Tinv,out)` gives `T·Asum·Tinv` (correct order); outputs distinct from scratch. |
| Homogeneous reduction exact; the prior `GetAx`-vs-`(AxP+AxM)` ~1e-9 worry was unfounded | **CONFIRMED** | Test 1.2 compares against the **independent** `GetAx`-based global Jacobian and matches at 7.1e-16 — `Ax_plus_+Ax_minus_ == Ax_` is an exact eigen-reconstruction (`godunov_flux.cpp:107-119`: `R·(Λ⁺+Λ⁻)·R⁻¹ = R·Λ·R⁻¹`). |
| Math is `½A_self·Q_self + ½A_nbr·Q_nbr`, NOT `½A_self·(Q_self+Q_nbr)` | **CONFIRMED** | `centralSelf` only ever uses `flux_self`, `centralNbr` only `flux_nbr`; `ApplyPerFaceFlux` (`:404-413`) multiplies each by its own Q. Test 1.3 nbr-term anchor passes (2.8e-16) and would FAIL if `A_nbr` were dropped (min_contrast 0.40). |
| Pass-1 refactor byte-identical | **CONFIRMED** | `ResolveFaceFluxOperands_` (`inl:246-286`) reproduces the exact check order (`!ftr` → `e2<0` → `face_bdr_attr_!=0`), the same centroid CalcOrtho normal, and the same `owned_flux_pool_->At(e1/e2)` refs as the pre-refactor inline code (`git show HEAD:...inl`). Parity 52/52 np4 confirms. |
| `BuildPerFaceCentralFluxMatrices_` dormant (not called from ctor) | **CONFIRMED** | ctor calls only `BuildPerFaceBimaterialFluxMatrices_` (`inl:50`); the subclass `SetMixedFluxMode` aborts on any non-None mode (`inl:795-803`) — see IMPL-8. |
| Shared-arm OOB risk on `At(e2)` for a shared face | **REFUTED (safe)** | For a shared/partition-seam face MFEM leaves `Elem2No < 0` (`wave_operator.inl:2566-2574` comment + the `e2<0 && shared_mesh_face_set_.count` guards at `:2574,:3886`); `ResolveFaceFluxOperands_`'s `if (e2<0) return false` (`inl:262`) skips it → the interior arm `continue`s, no OOB `At(e2)`, no wrong-material build. |
| Double-build (interior arm + shared arm on the same face) | **REFUTED** | Interior arm skips shared faces (returns false on `e2<0`); shared arm handles them by `sf`. Each central face built exactly once. |
| `seam_continuous` guard collective-safe / rank-uniform | **CONFIRMED (acceptable)** | `material_->mode` and `seam_continuous_` are both rank-uniform (same config on every rank), so the predicate is identical on all ranks; `MFEM_VERIFY`→`MPI_Abort` tears down the job (no deadlock). `material_` is non-null post-ctor (`inl:46`). |

---

## Findings

### IMPL-1 — `ResolveFaceFluxOperands_` `side==1` branch is dead and untested
**Severity:** WARNING (dead code / latent)
**File:** `dynamic/bimaterial_wave_operator.inl:283-284` (the `else { flux_self=&fe2; flux_nbr=&fe1; }` branch); declared `hpp:201`.
**Description:** Both call sites pass `side=0` (`inl:407`, `inl:556`). The header/comment
says `side` is "Used by the UPWIND build (per-side swap)", but the upwind build
(`build_face_matrices`, `inl:317-361`) does its own per-side swap via
`compose_side(flux_L.GetAx())` / `compose_side(flux_R.GetAx())` — it never calls
`ResolveFaceFluxOperands_(…, side=1, …)`. So the `side==1` swap path has **zero**
execution and **zero** test coverage; a sign/wiring error there would not be caught.
**Trigger:** A future Phase-3 caller that passes `side=1` would exercise untested code.
**Suggested fix:** Either (a) drop the `side` parameter entirely and hardcode
self=Elem1/nbr=Elem2 (the only used path), or (b) keep it but add a primitive unit
test that asserts `side=1` returns the swapped `(fe2, fe1)` with the SAME `nor`. Prefer (a)
until a real `side=1` caller exists — fewer untested branches.
**Test (if kept):** `Test NEW-I.1 ResolveFaceFluxOperands_SideSwap` — build a 2-element
fixture, call with `side=0` and `side=1`, assert pointers swap and `nor` is identical.

### IMPL-2 — Test 1.3 heterogeneous reference is builder-self-referential (no independent het anchor)
**Severity:** WARNING (test rigor)
**File:** `tests/unit/test_bimaterial_central_flux.cpp:280-296`.
**Description:** Test 1.3 builds the reference `halfA1`/`halfA2` by calling the SAME
`BuildPerFaceCentralMatricesGlobal` in homogeneous form (`(g1,g1)`, `(g2,g2)`). This
correctly isolates the **(self,nbr) WIRING** (the BUG-3 anchor — a dropped `A_nbr` is
caught because `halfA2` comes from `g2`, confirmed by min_contrast 0.40), but it does
**not** independently validate the heterogeneous matrix *values*: a defect common to all
three builder calls (e.g. a wrong rotation order, a `Tinv`/`T` swap) would cancel and the
test would still pass. The homogeneous case IS anchored independently (Test 1.2 vs
`GetAx`), so the suite as a whole is sound — but the het *values* rest on Test 1.1/1.2
transitively, not directly.
**Trigger:** A rotation/transpose bug that affects `centralSelf` and `centralNbr`
identically would pass Test 1.3 (caught only by Test 1.1's `Central` comparison, which is
homogeneous).
**Suggested fix:** Add to Test 1.3 an independent het anchor: for the same `(g1,g2,nor)`,
build `A_g1 = T·g1.GetAx()·Tinv` and `A_g2 = T·g2.GetAx()·Tinv` via the test's own
`BuildGlobalJacobianFromGetAx` (already present, `:148-160`) and assert
`ApplyPerFaceFlux(cS12,cN12,Q,Q,F) == ½(A_g1+A_g2)·Q` to ≤1e-11. This anchors the het
output against the `GetAx` path, not the builder.
**Test:** `Test NEW-I.2 CentralFlux_Heterogeneous_IndependentGetAxAnchor` (HIGH for het
confidence; MEDIUM overall since Test 1.1/1.2 cover the math transitively).

### IMPL-3 — Centroid unit-normal logic is triplicated
**Severity:** SUGGESTION (DRY / maintainability)
**File:** three copies: `compute_centroid_unit_normal` lambda (`inl:374-390`, used by
Godunov Pass 2), `ResolveFaceFluxOperands_` body (`inl:265-279`), and the central
shared arm (`inl:631-641`).
**Description:** The same `GetCenter → SetAllIntPoints → CalcOrtho → /Norml2` sequence is
written three times. The plan (Phase 2 req 3) explicitly offered "extend
`ResolveFaceFluxOperands_` to accept an `sf`" OR "inline the shared lookup"; the impl
inlined, which is acceptable, but the centroid-normal helper should still be shared.
**Trigger:** A future edit to the normal convention (e.g. orientation sign) must be made
in three places; missing one silently desyncs the Godunov and central normals on shared
faces.
**Suggested fix:** Hoist a private `void CentroidUnitNormal_(FaceElementTransformations*, real_t[3]) const`
and call it from all three sites (the lambda, `ResolveFaceFluxOperands_`, and the central
shared arm). No behavior change; parity tests still gate it.

### IMPL-4 — Central shared arm silently `continue`s on an out-of-range mesh-face index where Godunov Pass 2 fails loud
**Severity:** SUGGESTION (consistency / defensive)
**File:** `dynamic/bimaterial_wave_operator.inl:578-579` (central) vs `:431-435` (Godunov).
**Description:** Godunov Pass 2 asserts `MFEM_VERIFY(mesh_face_idx >= 0 && < n_faces, …)`
(fail-loud, "should never happen" bug detector). The central shared arm instead does
`if (mesh_face_idx < 0 || mesh_face_idx >= n_faces) { continue; }` (silent skip). Both
reach the same valid faces in practice, so this is not a correctness defect — but a
genuinely corrupt `sf→mesh_face` map would be silently dropped on the central path and
loudly caught on the Godunov path, an asymmetry that complicates debugging.
**Suggested fix:** Match the Godunov path: use the same `MFEM_VERIFY` (the map is built by
MFEM; an out-of-range index is a real invariant violation, not an expected skip).

### IMPL-5 — Debug disjointness assert omits shared fault faces
**Severity:** SUGGESTION (debug-only assert completeness)
**File:** `dynamic/bimaterial_wave_operator.inl:535-539` + `:547-549`.
**Description:** `fault_face_idx_set` is built only from `fault_interior_faces_`; it does
NOT add `fault_shared_faces_ → GetSharedFace(sf)` (the way `BuildCentralFluxFaceSet_`
itself does at `wave_operator.inl:1703-1708`). So the `MFEM_ASSERT(fault_face_idx_set.count(mesh_face)==0)`
disjointness check at `:547` would NOT catch a shared fault face wrongly present in
`central_flux_face_set_`. The base `BuildCentralFluxFaceSet_` guarantees exclusion of
shared fault faces, so this is only a weakened debug net — but it under-delivers on its
stated intent ("central_flux_face_set_ must exclude fault faces").
**Suggested fix:** In the `IsParallelMesh` branch, also insert
`pmesh.GetSharedFace(fault_shared_faces_[i])` into `fault_face_idx_set` (mirroring
`wave_operator.inl:1699-1710`) before the interior-arm loop.

### IMPL-6 — `std::uniform_real_distribution` makes the "fixed seed → deterministic" test non-portable across stdlib implementations
**Severity:** SUGGESTION (LOW)
**File:** `tests/unit/test_bimaterial_central_flux.cpp:113-115,129,139,255` (with
`std::mt19937 rng(20260604u)`, `:168`).
**Description:** `std::mt19937` is portable, but `std::uniform_real_distribution` is **not**
specified to produce the same sequence across libstdc++ / libc++ / MSVC. The comment
"fixed seed -> deterministic" (`:168`) is true only on a fixed stdlib. Because the test is
a tolerance test (asserts `worst <= 1e-11`, not a golden value), this does not cause false
failures today — but if anyone later adds a golden-value or cross-platform-reproducibility
assertion, it would break.
**Suggested fix:** Either (a) replace `uniform_real_distribution` with an explicit
`(rng() / double(rng.max()))` affine map (fully portable), or (b) add a one-line comment
that determinism holds per-stdlib only and the test asserts a tolerance, not a golden
sequence. (b) is sufficient given the current assertions.

### IMPL-7 — `n_faces` is unused in a serial (`Mesh`) instantiation of the central build
**Severity:** SUGGESTION (LOW — future `-Werror` hazard)
**File:** `dynamic/bimaterial_wave_operator.inl:530` (`const int n_faces = mesh_.GetNumFaces();`).
**Description:** `n_faces` is referenced only inside the `if constexpr (IsParallelMesh<MeshType>::value)`
block (`:579`). For a serial `BimaterialWaveOperator<Mesh>` instantiation (which exists —
`test_phaseh_wave_operator_constant_parity.cpp:137,400`,
`test_wave_operator_spatial_derivative.cpp:295,340`), that block is discarded and `n_faces`
becomes an unused variable. Today no warning fires because the method is dormant (not
ODR-used for serial), but once Phase 3 wires the call, a serial instantiation would emit
`-Wunused-variable` (and break a `-Werror` build).
**Suggested fix:** Move `const int n_faces = mesh_.GetNumFaces();` inside the
`#ifdef MFEM_USE_MPI` / `if constexpr` parallel block, or mark `(void)n_faces;`. (The
Godunov build does not have this issue because it uses `n_faces` in the serial `assign`,
`inl:296`.)

### IMPL-8 — Lifecycle gap: subclass `SetMixedFluxMode` aborts before `central_flux_face_set_` is ever populated, so the dormant central build would today iterate an EMPTY set
**Severity:** WARNING (latent integration gap; informational for Phase 1+2 since dormant)
**File:** `dynamic/bimaterial_wave_operator.inl:795-803` (the override aborts on any
non-None mode) vs the base `WaveOperator::SetMixedFluxMode → BuildCentralFluxFaceSet_`
(`wave_operator.inl:1639-1641`) being the ONLY populator of `central_flux_face_set_`.
**Description:** `BuildPerFaceCentralFluxMatrices_` iterates `central_flux_face_set_`
(`inl:545`) and the shared-arm filters on it (`inl:583`). But on a `BimaterialWaveOperator`
that set is **never** populated, because the subclass override hard-aborts on
`Adjacent`/`AllContinuous` (R-003) before the base `BuildCentralFluxFaceSet_` can run. So if
Phase 3 were to call `BuildPerFaceCentralFluxMatrices_` *without first lifting this abort*,
it would silently build ZERO central matrices (empty map) and the mixed-flux feature would
be a no-op — not a crash, a **silent** no-op. This is consistent with the plan (Phase 3
"lift the mutual-exclusion abort"), so it is not a Phase-1/2 defect; it is flagged so the
Phase-3 implementer does not wire the build before relaxing the guard AND before some path
populates `central_flux_face_set_` on this subclass.
**Trigger:** Phase 3 calls `BuildPerFaceCentralFluxMatrices_` while `SetMixedFluxMode`
still aborts ⇒ empty central map ⇒ dispatch never sees a central face ⇒ feature inert with
no error.
**Suggested fix (for Phase 3):** When lifting the R-003 abort, ensure the relaxed
`SetMixedFluxMode` override calls `Base::SetMixedFluxMode(m)` (to populate
`central_flux_face_set_` via `BuildCentralFluxFaceSet_`) and THEN
`BuildPerFaceCentralFluxMatrices_()`. Add an `MFEM_ASSERT(per_face_central_flux_.size() ==
<expected central-shared+interior count>)` post-build sanity check (the Phase-2 acceptance
criterion "size == number of central-set mesh-faces") so an empty map is loud, not silent.
**Test:** the Phase-3 single-valuedness test (Test 3.2) will exercise this; add an explicit
`per_face_central_flux_` non-empty / size assertion to the Phase-2/3 fixture per the plan's
own acceptance bullet ("`per_face_central_flux_` size == number of central-set mesh-faces").

---

## Proposed Unit Tests (Implementation Rev I-1)

| ID | Name | Target | Validates | Priority |
|----|------|--------|-----------|----------|
| NEW-I.1 | `ResolveFaceFluxOperands_SideSwap` | `BimaterialWaveOperator::ResolveFaceFluxOperands_` | `side=1` returns swapped `(fe2,fe1)` with identical `nor` (covers the IMPL-1 dead branch) — or delete the branch | LOW (or drop branch) |
| NEW-I.2 | `CentralFlux_Heterogeneous_IndependentGetAxAnchor` | `BimaterialFlux::BuildPerFaceCentralMatricesGlobal` | het `F*(Q,Q) == ½(A_g1+A_g2)·Q` vs the `GetAx`-based global Jacobians (independent of the builder; closes the IMPL-2 circularity for het values) | MEDIUM |
| NEW-I.3 | `CentralBuild_SizeAndDisjointness` (Phase 3 fixture) | `BuildPerFaceCentralFluxMatrices_` | after a (relaxed) `SetMixedFluxMode(Adjacent)`, `per_face_central_flux_.size()` == #central faces and ∩ fault == ∅ (closes IMPL-8 silent-empty) | MEDIUM |

Skeleton (NEW-I.2), drop into `test_bimaterial_central_flux.cpp` Test 1.3 loop:
```cpp
// Independent het anchor (does NOT reuse the builder for the reference):
DenseMatrix Ag1, Ag2;
BuildGlobalJacobianFromGetAx(g1, nor, Ag1);
BuildGlobalJacobianFromGetAx(g2, nor, Ag2);
real_t Fhet[NUM_STATE], Fref[NUM_STATE];
BimaterialFlux::ApplyPerFaceFlux(cS12, cN12, Q, Q, Fhet);   // ½A_g1·Q + ½A_g2·Q
for (int i=0;i<NUM_STATE;++i){ real_t s=0; for(int j=0;j<NUM_STATE;++j)
   s += 0.5*(Ag1(i,j)+Ag2(i,j))*Q[j]; Fref[i]=s; }
worst_13_anchor = std::max(worst_13_anchor, WorstRel(Fhet, Fref));
// assert worst_13_anchor <= 1e-11
```

---

## Priority Summary (Implementation Rev I-1)

| ID | Severity | Status | Description |
|----|----------|--------|-------------|
| IMPL-1 | WARNING | OPEN | `ResolveFaceFluxOperands_` `side==1` branch dead/untested |
| IMPL-2 | WARNING | OPEN | Test 1.3 het reference is builder-self-referential (no independent het anchor) |
| IMPL-3 | SUGGESTION | OPEN | Centroid unit-normal logic triplicated |
| IMPL-4 | SUGGESTION | OPEN | Central shared arm silent `continue` vs Godunov fail-loud `MFEM_VERIFY` |
| IMPL-5 | SUGGESTION | OPEN | Debug disjointness assert omits shared fault faces |
| IMPL-6 | SUGGESTION | OPEN | `uniform_real_distribution` non-portable (LOW; tolerance test, not golden) |
| IMPL-7 | SUGGESTION | OPEN | `n_faces` unused in serial instantiation (future `-Werror`) |
| IMPL-8 | WARNING | OPEN | Lifecycle: dormant central build would iterate an empty `central_flux_face_set_` until Phase 3 lifts the R-003 abort |

**No CRITICAL or HIGH (MODERATE) implementation finding.** The two priority adversarial
targets — the shared-arm OOB (#3) and the central math (#1) — are both **clean**: shared
faces are safely rejected by the `e2<0` check (verified against MFEM's documented
`Elem2No==-1`-on-seam convention), and the central composition is mathematically and
numerically correct (independent `GetAx` anchor at 7e-16; het wiring anchor at 2.8e-16).

## Verdict (Phase 1 & 2 implementation)

**PASS WITH FIXES.** The Phase-1 primitive and the Phase-2 precompute/refactor are
correct, byte-exact (52/52 np4 parity), and free of the assumed critical bugs (no OOB, no
double-build, no aliasing, correct `½A_self·Q_self + ½A_nbr·Q_nbr` math, dormant as
claimed). The eight findings are all WARNING/SUGGESTION: dead `side==1` branch (IMPL-1),
a not-fully-independent het test reference (IMPL-2), DRY/consistency/portability nits
(IMPL-3..7), and one latent lifecycle gap the Phase-3 implementer must honor (IMPL-8,
silent-empty-set if the build is wired before the R-003 abort is lifted). None blocks
merging Phase 1+2 in isolation; IMPL-2 and IMPL-8 should be addressed before/with Phase 3
to avoid a silently-inert feature.

---

# Implementation Review — Phase 3 (Rev I-2)

**Date:** 2026-06-04
**Reviewer:** chunhui-code-reviewer agent (adversarial; "assume ≥3 bugs")
**Scope (Phase 3 diff):** `dynamic/bimaterial_wave_operator.inl` (`SetMixedFluxMode`
rewrite, `InteriorFaceFlux_`/`SharedInteriorFaceFlux_` central branch, `ComputeMaxDt`
comment), `dynamic/bimaterial_wave_operator.hpp` (`using Base::mf_on_;`, doc), the base
comment in `dynamic/wave_operator.inl:1623`, NEW `tests/unit/test_bimaterial_mixed_flux_dispatch.cpp`,
`makefile` registration.
**Method:** re-built FRESH (forced `rm` of the four test `.o` + binaries, per the
implementer's stale-`.inl`-tracking note), re-ran all four targets, and built/ran a
throwaway np=2 ParMesh bimaterial-Adjacent probe to exercise the otherwise-untested shared
central path + the IMPL-8 assert.

## Build / test results reproduced (fresh objects)

| Target | Result | Notes |
|--------|--------|-------|
| `seas_test_bimaterial_mixed_flux_dispatch` (serial) | **10 / 10** | central set size 2, `f_c=8`, `f_n=2`; single-valued ≤1e-12, Godunov het ≠ |
| `seas_test_bimaterial_wave_operator_parity` (serial) | **9 / 9** | incl. C-6c ComputeMaxDt parity (homogeneous → None) |
| `seas_test_bimaterial_central_flux` (serial) | **6 / 6** | Phase-1 primitive incl. IMPL-2 independent GetAx anchor 5.2e-16 |
| `seas_test_phaseh_wave_operator_constant_parity` (serial) | **19 / 19** | `mixed_flux=none` byte-exact |
| `seas_test_phaseh_wave_operator_constant_parity` (np=4) | **52 / 52** | `mixed_flux=none` byte-exact, parallel |
| *throwaway* np=2 bimaterial Adjacent ParMesh probe | **SURVIVED** | rank0 set=1/map=1, rank1 set=2/map=2 (incl. 1 SHARED central face); re-entrancy None→0/0→re-Adjacent rebuild OK |

The implementer's "10/10" for the dispatch test is **confirmed independently** on freshly
recompiled objects. The byte-exact `mixed_flux=none` contract (parity 9/9, phaseh 19/19
serial + 52/52 np=4) is **confirmed** — the CENTRAL precompute log does NOT appear on the
None runs (empty set ⇒ no-op), and the central branch is provably unreachable when
`mf_on_==false` (`mf_on_` is `false` for None AND the set is cleared, so the
`&& central_flux_face_set_.count(...)` short-circuit is doubly false).

## Priority adversarial targets — verdicts

1. **Single-valuedness sign (#1, the crux) — CORRECT.** `InteriorFaceFlux_`
   (`bimaterial_wave_operator.inl:711-718`) computes `F_h_e1 = ½A_e1·Q_e1 + ½A_e2·Q_e2`
   via `ApplyPerFaceFlux(c[0],c[1],…)` then copies `F_h_e2[comp]=F_h_e1[comp]` for ALL
   `NUM_STATE` comps and `return`s before the Godunov per-side apply — no fallthrough, no
   double-write. The **same-sign** copy matches the scalar reference
   (`wave_operator.inl:1209`, `F_h_e2[c]=F_h_e1[c]`); the scalar path does NOT negate for
   e2. The DG assembly applies the per-element outward sign downstream
   (`wave_operator.inl:3130-3133`: `rhs[e1] -= w·shape1·F_h_e1`, `rhs[e2] += w·shape2·F_h_e2`),
   so the same-sign deposit is physically correct. Test 3.2 (≤1e-12 under genuine contrast)
   + the "central dispatch == ApplyPerFaceFlux(per_face_central_flux_[f_c]) bit-for-bit"
   value check both pass. **No bug.**
2. **IMPL-8 lifecycle assert (#2) — CORRECT, no false-fire / no false-pass.** The assert
   `per_face_central_flux_.size() == central_flux_face_set_.size()`
   (`bimaterial_wave_operator.inl:813`) is sound because **every map key is provably a
   subset of the set** (interior arm inserts `per_face_central_flux_[mesh_face]` where
   `mesh_face` iterates the set; shared arm inserts only when
   `central_flux_face_set_.count(mesh_face_idx) > 0`). A subset of equal cardinality IS the
   whole set ⇒ exact-equality is NOT too strict and a false-PASS (equal sizes, wrong
   keys) is impossible. Empirically verified at np=2: rank0 1/1, rank1 2/2 including a
   shared central face inserted via the Allgatherv merge (whose LOCAL element is NOT
   fault-adjacent — the remote drove it — and the sf arm still built it). The 2-tet
   degenerate fixture (empty set) no-ops at 0==0. The only theoretical false-fire is a
   central-set SHARED face whose `GetSharedFaceTransformations(sf)` returns null (sf-arm
   `continue` with no map entry); this does not occur for a valid shared face and is not
   reachable on the tested fixtures — see P3-2 for the (low) residual robustness note.
   `<=` would be strictly weaker and is NOT recommended.
3. **Dispatch gate + re-entrancy (#3) — CORRECT.** `mf_on_` is technically redundant with
   a cleared set (base sets `mf_on_=false` AND `central_flux_face_set_.clear()` for None,
   `wave_operator.inl:1643,1655`), but it is a cheap defensive guard that matches the scalar
   reference verbatim. Re-entrancy verified empirically: `Adjacent → None` clears set+map to
   0/0 (base clears the set; `BuildPerFaceCentralFluxMatrices_` clears the map first,
   `inl:282`); `None → Adjacent` rebuilds. No stale central matrices survive a None flip.
4. **`.at()` safety (#4) — CORRECT / consistent.** A central-set face with no map entry
   would `throw std::out_of_range` in the const Mult path — but the IMPL-8 `MFEM_VERIFY`
   (always-on, runs at `SetMixedFluxMode` time, before any Mult) makes that state
   unreachable. The two guards are consistent: the assert is the loud-at-setup failure; the
   `.at()` is a belt-and-suspenders that would only fire if someone mutated the set after
   setup (which the R-1408 lifecycle invariant forbids).
5. **Test rigor (#5) — adequate for the SERIAL central path; SHARED path UNTESTED (P3-1).**
   (a) The 5-hex fixture produces central set `{(1,2),(3,4)}` (size 2) and picks a 2-sided
   non-central, non-fault `f_n`; `f_n` IS heterogeneous (μ(x) gradient) so the Godunov
   discriminator (`max_godunov_diff > 1e-3`) is non-vacuous (passed). (b) The fixed-pair
   `>1e-3` threshold is not in the kernel of (A_e1−A_e2) for THIS fixture (passed with
   margin); it is a fixed-fixture test, not flaky, but see P3-3. (c) Test 3.2 genuinely
   exercises the central branch (the "== ApplyPerFaceFlux(per_face_central_flux_) bit-for-bit"
   check would fail if it routed to Godunov or to the base scalar's poison-(1,1,1)
   `flux_.Central`). (d) `TestableBimat`'s `using …InteriorFaceFlux_` reaches the OVERRIDE
   (proven by the value match). (e) `SetMixedFluxMode(Adjacent)` succeeds with NO
   `SetFaultFlux` — the base needs only `bc_.fault_attr>0` + ctor-populated fault face
   lists; valid.
6. **Byte-exact contract (#6) — CONFIRMED on fresh build.** See the table above.
7. **`ComputeMaxDt` divergent factors (#7) — DEFERRED-but-now-REACHABLE (P3-4, MODERATE).**

## Bugs Found — Phase 3 (Rev I-2)

### P3-1 — Bimaterial SHARED central-face path (Test 3.4a/b/c) has ZERO automated coverage
**Severity:** WARNING (MODERATE — untested production-reachable path)
**File:** test gap; code at `dynamic/bimaterial_wave_operator.inl:574-637` (shared arm of
`BuildPerFaceCentralFluxMatrices_`) + `:740-746` (`SharedInteriorFaceFlux_` central branch).
**Description:** Phase 3's plan lists Tests 3.4a (`Dispatch_SharedCentralFace_Build`), 3.4b
(`Dispatch_SharedCentralFace_RequiresSeamContinuous`), 3.4c (`Dispatch_DepthProfile_NoFalseAbort`)
as MPI np=2 acceptance criteria. NONE are implemented: `tests/unit/test_bimaterial_mixed_flux_dispatch.cpp`
is serial-only (the fixture is a serial `Mesh`), and the existing
`tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp` uses the **scalar** `WaveOperator`
(no `BimaterialWaveOperator`, no `seam_continuous`, no `per_face_central_flux_`). So the
entire bimaterial shared-central machinery — the `sf` loop, the `seam_continuous_` /
`Mode::Constant` abort gate (`:430-442`), the `shared_face_neighbour_material_` lookup, the
side-0 single-valued `SharedInteriorFaceFlux_` deposit — ships with no regression guard.
**Trigger:** any parallel `interior_flux="matrix"` + `mixed_flux="adjacent"` run with a fault
that crosses a partition seam (the production case). A regression in the sf arm (e.g. a
`GetSharedFace` mapping change, an off-by-one in `shared_face_bdr_attr_` sizing, a wrong
`material_->mode` branch) would pass all 87 current Phase-3-adjacent unit checks.
**Mitigation observed:** my throwaway np=2 probe (Constant material, hex-row, partition
{0,1}|{2,3,4}) exercised the path live and it is correct TODAY (rank0 built 1 shared central
face whose local element is not fault-adjacent; sizes matched; re-entrancy clean). So this is
a coverage gap, not a present defect — but it is the single largest risk for Phase 3 because
the most subtle code (cross-rank single-valued central, R-004 stub, seam-continuity gate) is
the part with no test.
**Suggested fix:** Port the planned 3.4a/b/c into `tests/parallel/` as an MPI (np=2)
`BimaterialWaveOperator<ParMesh>` test on a fault-crossing-seam fixture: 3.4a asserts the
shared central matrices build + the side-0 deposit matches `BuildPerFaceCentralMatricesGlobal`;
3.4b asserts a `Mode::Coefficient` shared central face WITHOUT `SetSeamContinuous(true)`
aborts with the seam-named message and WITH it builds; 3.4c asserts a depth-only
(`depth_profile_1d`) material + `seam_continuous=true` on a z-tilted seam builds (no false
abort). **(CRITICAL/MOD test.)**

### P3-4 — `BimaterialWaveOperator::ComputeMaxDt` mixed-flux factor is now REACHABLE but still divergent (no `cfl_rk_aware_` gating, no central+ADER guard)
**Severity:** WARNING (MODERATE — latent stability; reachable since R-003 lifted)
**File:** `dynamic/bimaterial_wave_operator.inl:779-789`.
**Description:** Before Phase 3, the `switch (mixed_flux_mode_)` in the bimaterial
`ComputeMaxDt` was dead (the `SetMixedFluxMode` abort kept `mixed_flux_mode_==None`). Phase 3
lifts that abort, so `mixed_flux_mode_` can now be `Adjacent`/`AllContinuous` here and the
**Adjacent→0.9 / AllContinuous→0.4** arms are LIVE. These diverge from the scalar base
(`wave_operator.inl:5883-5906`), which gates on `cfl_rk_aware_` (Adjacent→`0.6` RK / `0.9`
ADER; AllContinuous→`0.7` RK / `0.4` ADER). Consequences while Phase 4 is pending:
(i) a bimaterial central run under the REQUIRED RK integrator gets `dt·0.9` instead of the
RK target `dt·0.6` — **dt is ~50% too large**, risking the imaginary-axis instability the
RK factor exists to suppress; (ii) there is **no guard** in the bimaterial `ComputeMaxDt`
(nor anywhere on the matrix path) preventing `mixed_flux_mode_ != None && !cfl_rk_aware_`
(central + ADER), which BUILD §5.4/§6.3 documents as a runaway. The plan (Phase 4 BUG-2 +
item 3) explicitly defers the factor unification AND the central+ADER abort to Phase 4, and
the code comment at `:773-778` states the divergence is intentional and "a Frontera
calibration, not a local gate". So this is a documented deferral, NOT a plan violation — but
it is now a *reachable* divergence (it was unreachable dead code before), and a Phase-3
bimaterial central run on Frontera ahead of Phase 4 would be under-de-rated / unguarded.
**Trigger:** `interior_flux="matrix"` + `mixed_flux="adjacent|allcontinuous"` + RK4/RK45
(get 0.9/0.4 not 0.6/0.7) OR + ADER (no abort, certain runaway).
**Suggested fix:** Land Phase 4 BUG-2 (route both operators through the shared
`MixedFluxCflFactor_()`) AND Phase 4 item 3 (the `mixed_flux_mode_!=None && !cfl_rk_aware_`
ADER-instability abort) BEFORE any Frontera bimaterial central run. As an interim Phase-3
safety net, add the central+ADER abort to the bimaterial `ComputeMaxDt` now (one
`MFEM_VERIFY`), since that is the unambiguously-unstable combination. **(MOD test:** a unit
test asserting `ComputeMaxDt` aborts on Adjacent+`!cfl_rk_aware_` once the guard lands.)

### P3-2 — IMPL-8 assert false-fire is theoretically possible if a central-set SHARED face has a null `GetSharedFaceTransformations`
**Severity:** SUGGESTION (LOW — defensive robustness; not reachable on tested fixtures)
**File:** `dynamic/bimaterial_wave_operator.inl:578-580` (`if (!ftr) { continue; }`) vs the
IMPL-8 assert at `:813`.
**Description:** The sf arm skips a shared face on `!ftr` (null transformation) WITHOUT
building a map entry. If such a face were a member of `central_flux_face_set_` (it can be
inserted by the base's Allgatherv key-merge, `wave_operator.inl:1904-1914`, which does NOT
re-verify a non-null local transformation), the interior arm would also skip it
(`ResolveFaceFluxOperands_` returns false for a shared face), leaving a central-set face with
NO map entry ⇒ `per_face_central_flux_.size() < central_flux_face_set_.size()` ⇒ the IMPL-8
`MFEM_VERIFY` aborts a VALID run. `GetSharedFaceTransformations` returning null for a real
shared face is not expected (and did not occur in the np=2 probe), so this is a low-likelihood
defensive note, not a present bug.
**Suggested fix:** Either (a) make the `!ftr` case a loud `MFEM_VERIFY` for a central-set
shared face (consistent with the IMPL-4 fail-loud philosophy: a null transform on a
central-set face IS an invariant violation), or (b) document at `:578` that a null `ftr` on a
central-set shared face is impossible by construction (cite why) so the IMPL-8 assert cannot
false-fire on it.

### P3-3 — Test 3.1's `max_godunov_diff > 1e-3` discriminator is a fixed-fixture margin, not a structural guarantee
**Severity:** SUGGESTION (LOW — test robustness)
**File:** `tests/unit/test_bimaterial_mixed_flux_dispatch.cpp:240-243` (+ `MakeStatePairs`,
`:131-149`).
**Description:** The "Godunov is genuinely per-side" check relies on 4 hand-picked
`(Q_self,Q_nbr)` pairs producing `|F_e1−F_e2| > 1e-3·scale` on a specific μ(x) gradient.
This passed with healthy margin, so it is NOT flaky as written — but it is an empirical
property of the chosen fixture/pairs, not a structural guarantee. If a future edit narrows
the contrast (smaller μ slope) or changes the pairs, the discriminator could silently weaken
toward the kernel of (A_e1−A_e2) and the test would still "pass" while no longer proving the
per-face selection.
**Suggested fix:** Add a comment pinning the contrast magnitude as load-bearing, or assert a
minimum element-to-element μ ratio in the fixture so a future contrast reduction trips an
explicit guard rather than silently de-fanging the discriminator. Optionally cross-check
`max_godunov_diff` against an independent lower bound from the stored
`per_face_bimaterial_flux_[f_n]` matrices.

## Verification of prior IMPL-1..8 fixes (carried into the Phase-3 tree)

| Prior ID | Status in Phase-3 tree | Evidence |
|----------|------------------------|----------|
| IMPL-1 | FIXED | `ResolveFaceFluxOperands_` signature has NO `side` param (`bimaterial_wave_operator.inl:163-181`); every caller uses self=Elem1/nbr=Elem2. |
| IMPL-2 | FIXED | `test_bimaterial_central_flux` Test 1.3 uses an independent GetAx anchor (ran 5.2e-16). |
| IMPL-3 | FIXED | `CentroidUnitNormal_` extracted (`:132-148`); Godunov build + interior operand + central shared build all call it. |
| IMPL-4 | FIXED | central shared arm has the fail-loud `MFEM_VERIFY` on the sf→mesh_face map (`:582-585`). |
| IMPL-5 | FIXED | `fault_face_idx_set` now also inserts `GetSharedFace(fault_shared_faces_[i])` (`:294-304`). |
| IMPL-6 | N/A this round | `uniform_real_distribution` portability nit in `test_bimaterial_central_flux.cpp` (untouched; still LOW). |
| IMPL-7 | FIXED | `const int n_faces` moved inside the `#ifdef MFEM_USE_MPI` / `if constexpr` block (`:336`); serial build clean (no `-Wunused-variable`). |
| IMPL-8 | FIXED (Phase 3) | the abort is lifted; `SetMixedFluxMode` delegates to base → `BuildPerFaceCentralFluxMatrices_` → size-equality `MFEM_VERIFY` (`:838-843`). The chosen size-equality form is *stronger* than the originally-suggested fixed-count form (subset invariant ⇒ no false-pass). No silent-empty: a non-empty set with a missing build aborts loudly. |

None of the IMPL-1..7 fixes broke Phase 3 (all four suites green on fresh objects).

## Proposed Unit Tests (Implementation Rev I-2)

| ID | Name | Target | Validates | Priority |
|----|------|--------|-----------|----------|
| NEW-I.4 (NEW in Rev I-2) | `Dispatch_SharedCentralFace_Build` (MPI np=2) | sf arm of `BuildPerFaceCentralFluxMatrices_` + `SharedInteriorFaceFlux_` central branch | a fault-crossing-seam ParMesh builds the side-0 shared central matrices; the dispatched side-0 `F*` == `BuildPerFaceCentralMatricesGlobal(local,nbr)` apply | **HIGH** (closes P3-1) |
| NEW-I.5 (NEW in Rev I-2) | `Dispatch_SharedCentralFace_RequiresSeamContinuous` (MPI np=2) | the `Mode::Constant \|\| seam_continuous_` gate (`inl:430-442`) | a `Mode::Coefficient` shared central face WITHOUT `SetSeamContinuous(true)` aborts with the seam-named message; WITH it builds | **HIGH** (closes P3-1) |
| NEW-I.6 (NEW in Rev I-2) | `Dispatch_DepthProfile_NoFalseAbort` (MPI np=2) | the seam gate does not false-trip on legitimate depth variation | a depth-only material + `seam_continuous=true` on a z-tilted seam BUILDS (no abort) | MEDIUM (closes P3-1, BUG-16) |
| NEW-I.7 (NEW in Rev I-2) | `BimatComputeMaxDt_CentralAderAborts` | `BimaterialWaveOperator::ComputeMaxDt` | once the Phase-4 / interim guard lands: `mixed_flux_mode_!=None && !cfl_rk_aware_` aborts with the ADER-instability message | MEDIUM (P3-4 regression guard) |

## Priority Summary (Implementation Rev I-2)

| ID | Severity | Status | Description |
|----|----------|--------|-------------|
| P3-1 | WARNING (MODERATE) | OPEN | Bimaterial SHARED central-face path (Test 3.4a/b/c) has zero automated coverage |
| P3-4 | WARNING (MODERATE) | OPEN | `ComputeMaxDt` mixed-flux factor now reachable but ungated (no `cfl_rk_aware_`, no central+ADER abort) — deferred to Phase 4 per plan |
| P3-2 | SUGGESTION (LOW) | OPEN | IMPL-8 assert could false-fire on a central-set shared face with null `GetSharedFaceTransformations` (not reachable on tested fixtures) |
| P3-3 | SUGGESTION (LOW) | OPEN | Test 3.1 Godunov discriminator is a fixed-fixture margin, not a structural guarantee |

**No CRITICAL or HIGH finding.** The four priority crux items — single-valuedness sign
(#1), the IMPL-8 assert (#2), the dispatch gate + re-entrancy (#3), and `.at()` safety (#4) —
are all **clean and verified** (the sign matches the scalar oracle; the assert is robust by
the subset invariant and does not false-fire/false-pass; re-entrancy clears both set and map;
`.at()` is unreachable behind the always-on assert). The two MODERATE findings are a test
coverage gap (P3-1) and a documented-but-now-reachable CFL deferral (P3-4); the two LOW
findings are defensive/test-robustness nits.

### Test Coverage Summary (cumulative)

| Suite | File | Tests | Phase-3 reproduced |
|-------|------|-------|--------------------|
| Dispatch | `test_bimaterial_mixed_flux_dispatch.cpp` | 10 | 10/10 serial |
| Bimat parity | `test_bimaterial_wave_operator_parity.cpp` | 9 | 9/9 serial |
| Central primitive | `test_bimaterial_central_flux.cpp` | 6 | 6/6 serial |
| phaseh parity | `test_phaseh_wave_operator_constant_parity.cpp` | 19 / 52 | 19/19 serial, 52/52 np=4 |
| **Bimaterial shared central (3.4a/b/c)** | — | **0** | **MISSING (P3-1)** |

## Verdict (Phase 3 implementation)

**PASS WITH FIXES.** Phase 3 correctly lifts the R-003 abort and dispatches a SINGLE-VALUED
bi-material central flux on central-set faces alongside the bi-material Godunov upwind. The
hardest correctness questions are all answered correctly: the same-sign `F_h_e2 = F_h_e1`
copy matches the scalar oracle and the downstream `-=/+=` assembly sign; the IMPL-8
size-equality lifecycle assert is robust (no false-fire, no false-pass, proven by the
map⊆set subset invariant and verified live at np=2); re-entrancy (`None → Adjacent`)
correctly clears and rebuilds both the set and the matrix map; `.at()` is unreachable behind
the always-on assert; and `mixed_flux=none` is byte-exact (9/9, 19/19 serial, 52/52 np=4 on
freshly recompiled objects). The assumed "≥3 bugs" did NOT materialize as correctness
defects — the four findings are a MODERATE test-coverage gap on the shared central path
(P3-1, the top action item), a MODERATE documented-but-now-reachable CFL divergence deferred
to Phase 4 (P3-4, add at least the central+ADER abort as an interim guard), and two LOW
defensive/test-robustness nits (P3-2, P3-3). Merge Phase 3 after adding the MPI shared-central
tests (NEW-I.4/5/6) and either landing Phase 4 or adding the interim central+ADER guard.

---

# Implementation Review — Phase 4 (Rev I-4)

**Date:** 2026-06-04
**Reviewer:** chunhui-code-reviewer agent (adversarial; "assume ≥3 bugs")
**Scope (Phase 4 diff):** `dynamic/wave_operator.{hpp,inl}` (extract protected
`MixedFluxCflFactor_()`; scalar `ComputeMaxDt` now has the central+ADER `MFEM_VERIFY` guard
then `return MixedFluxCflFactor_() * cfl * h_min_ / flux_.GetCp()`); `dynamic/bimaterial_wave_operator.inl`
(`ComputeMaxDt` deletes its divergent `0.9/0.4` switch, calls `this->MixedFluxCflFactor_()`, same guard);
`tests/unit/test_wave_operator.cpp` (Test 15 `TestR003MatrixMixedFluxAborts` → `TestMatrixMixedFluxEnabled`);
NEW `tests/unit/test_bimaterial_mixed_flux_cfl.cpp` (Tests 4.1/4.2); `makefile` registration.
**Method:** rebuilt FRESH — forced `rm` of every touched test `.o` + binary AND
`dynamic/wave_operator.o` + `dynamic/bimaterial_wave_operator.o` before each rebuild (per the
implementer's stale-`.inl`-tracking note). Re-ran all suites; did NOT trust the fix report's
counts. Independently swept `tests/ drivers/ jobs/ spatial/ safs/` + every `*.toml` for a
scalar/bimaterial config that reaches `ComputeMaxDt` with `mixed_flux!=none` under ADER.

## Build / test results reproduced (fresh objects)

| Suite | Target | Result reproduced |
|-------|--------|-------------------|
| **CFL (NEW, Tests 4.1/4.2)** | `seas_test_bimaterial_mixed_flux_cfl` (serial) | **17 / 17** |
| Wave operator (Test 15 rewrite + Test 16) | `seas_test_wave_operator` (serial) | **25 / 25** |
| Bimat parity (C-6c ComputeMaxDt) | `seas_test_bimaterial_wave_operator_parity` (serial) | **9 / 9** |
| phaseH constant parity (C-3 bit-identical dt) | `seas_test_phaseh_wave_operator_constant_parity` | **19 / 19** serial, **52 / 52** np=4 |
| Bimat central primitive | `seas_test_bimaterial_central_flux` (serial) | **6 / 6** |
| Bimat mixed-flux dispatch | `seas_test_bimaterial_mixed_flux_dispatch` (serial) | **10 / 10** |
| Bimat shared central (3.4a/c) | `seas_test_bimaterial_mixed_flux_shared` (np=2) | **14 / 14** |
| Scalar adjacent dispatch | `seas_test_mixed_flux_dispatch_adjacent` (serial) | **9 / 9** |

**Reproduced total: 161 passed, 0 failed** (17+25+9+19+52+6+10+14+9). No new compiler warnings.

## Priority adversarial targets — verdicts

1. **Scalar `ComputeMaxDt` byte-exactness (acceptance #1) — CONFIRMED byte-identical.**
   Old: `cfl_mixed_flux_factor = 1.0` (None arm) then `return cfl_mixed_flux_factor * cfl * h_min_ / flux_.GetCp();`.
   New: `MixedFluxCflFactor_()` returns the literal `1.0` for None, then
   `return MixedFluxCflFactor_() * cfl * h_min_ / flux_.GetCp();` (`wave_operator.inl:5904`). Identical
   expression, identical evaluation order (`1.0 * cfl * h_min_ / cp`), no FP reordering. Verified live: phaseH
   C-3 "ComputeMaxDt bit-identical between scalar and hetero ctors" PASS; bimat parity C-6c PASS; both on
   freshly recompiled objects. **No byte-exactness break.**
2. **Scalar + mixed + ADER abort — the riskiest decision. NOT a default-path regression, but a
   real behaviour change with a latent break for direct `*_spatial.toml` runs (see P4-1).**
   - **Existing tests:** every other `ComputeMaxDt` caller in `tests/` (test_wave_bc, test_pml,
     test_wave_operator, parity, parallel, tpv102/tpv205 verification) calls it with the default
     `mixed_flux_mode_==None` and never `SetMixedFluxMode(non-None)` first → guard never fires → no test regresses.
     The scalar mixed-flux dispatch tests (`test_mixed_flux_dispatch_adjacent{,_p2}`, `_all_continuous`,
     `_none`, `test_mixed_flux_adjacent_mpi`, `test_lsw_rk_mixed_flux`, `test_mixed_flux_face_set`) call
     `ComputeMaxDt` **zero** times (grep-verified).
   - **Native TPV102/104/205 drivers (`seas_tpv102_driver` etc.): SAFE by call ORDERING.** `ComputeMaxDt`
     runs at `tpv102_driver.cpp:1266` / `tpv104:1323` / `tpv205:1243` — BEFORE `SetMixedFluxMode` at
     `tpv102:1613` / `tpv104:1672` / `tpv205:1633`. So `mixed_flux_mode_==None` when `ComputeMaxDt` runs, the
     guard does not fire, and the dt was already factor-1.0 both before and after Phase 4 (these drivers never
     applied the mixed-flux de-rating). The gold ADER-O2 + `--mixed-flux adjacent` run
     (`tpv102/gold/results_mixed_flux_adjacent_p1_O2_dev_job7681589`) uses this native driver → unaffected.
   - **Spatial driver production sbatch jobs: SAFE.** `jobs/{tpv102,tpv104,tpv205}_spatial/*_rk45_mixedflux*.sbatch`
     use the `*_rk45_mixedflux.toml` (RK ⇒ `SetCflRkAware(true)` at `spatial_dyn_driver.cpp:1902` ⇒ guard
     passes). `jobs/*_spatial/*aderO*_pureupwind*.sbatch` pass `--mixed-flux none` on the CLI (e.g.
     `tpv104_p2_aderO3_pureupwind_200m_normal.sbatch:232`) overriding the TOML's `adjacent`. SAFS jobs default
     `SAFS_MIXED_FLUX=none` (`spatial_dyn_..._safs.sbatch:165`). None of these trips the guard.
   - **THE GAP (P4-1):** the spatial driver (`WaveOperator<ParMesh>`, scalar branch) sets
     `SetMixedFluxMode` at `:1245` BEFORE `ComputeMaxDt` at `:1903-1906`, and the on-disk configs
     `tpv102/configs/tpv102_spatial{,_p2,_p3}.toml`, `tpv205_spatial{,_p2,_p3}.toml`, and several
     `safs/.../*.toml` ship `mixed_flux="adjacent"` + `interior_flux="scalar"` with the **default
     `time_integrator="ader"`** (`spatial_friction.cpp:1133`). G1 (`spatial_friction.cpp:1162`) only
     forbids `matrix+mixed`; G2 (`spatial_dyn_driver.cpp:813`) only forbids `rk+matrix`; NO guard blocks
     `scalar+adjacent+ADER`. Pre-Phase-4 that path returned `0.9*cfl*h_min/cp` (finite). Post-Phase-4 it
     **aborts**. So `seas_spatial_dyn_driver --config tpv102_spatial.toml` (no override), or any job with
     `SAFS_MIXED_FLUX=adjacent` under ADER, now FAILS where it previously ran. This is a **deliberate
     behaviour change per plan §scope-Integrator + Phase 4 req 3** ("mixed flux REQUIRES RK, exactly as on
     the scalar path"), but it is **incompletely landed**: the shipped default TOMLs were not migrated and
     `spatial_dyn_driver.cpp:1899-1901`'s claim "the ADER branch is byte-unchanged … keeps the pre-Phase-14
     dt to the bit" is now FALSE for scalar+adjacent+ADER (it aborts, not de-rates). Recorded as **P4-1 (HIGH)**.
3. **`MixedFluxCflFactor_` values vs the OLD switch — EXACT match.** New helper
   (`wave_operator.inl:5847-5872`): `None→1.0`, `Adjacent→cfl_rk_aware_?0.6:0.9`,
   `AllContinuous→cfl_rk_aware_?0.7:0.4`, `default→MFEM_ABORT` (R-1600). Identical to the pre-Phase-4
   scalar switch. Test 4.1 asserts all six `(mode×rk)` combos with exact-double compare (`fw==c.expected`),
   reproduced 12/12 PASS. **No drift.**
4. **`return 1.0;` after `MFEM_ABORT` — dead and harmless; compiles clean.** `MFEM_ABORT` expands to
   `mfem_error`, declared `[[noreturn]]` (`mfem/general/error.hpp:62`). The trailing
   `return 1.0;   // unreachable … silences -Wreturn-type` (`wave_operator.inl:5871`) is genuinely
   unreachable; the `[[noreturn]]` attribute already silences `-Wreturn-type`, so the explicit return is
   redundant (not wrong). Fresh `-O2` build emits NO warning. Cannot mask a missing case (every enumerator
   has an explicit arm; the `default` only catches future-added modes and aborts). LOW nit recorded as P4-3.
5. **Guard consistency — IDENTICAL condition on both operators.** Verbatim
   `MFEM_VERIFY(mixed_flux_mode_ == MixedFluxMode::None || cfl_rk_aware_, …)` at
   `wave_operator.inl:5888` and `bimaterial_wave_operator.inl:803`. `cfl_rk_aware_` is reachable in the
   override via `using Base::cfl_rk_aware_;` (`bimaterial_wave_operator.hpp:80`, P3-4). Both place the
   guard BEFORE the factor lookup, so `Adjacent+ADER` aborts identically on both (bimat before its
   per-element walk; scalar before its `h_min_/cp` formula). Messages differ in prose only (same intent).
   **Cannot diverge.**
6. **Test 4.1 rigor — adequate.** Probes BOTH operators' `MixedFluxCflFactor_()` through `using`-re-exported
   protected accessors (`TestableWave` / `TestableBimat`). The exact-double compare is valid (same source
   literals on both sides). `SetMixedFluxMode(Adjacent)` succeeds on the 5-hex fault fixture (BoundaryConfig
   `fault_attr=3`, a tagged fault face at x=3, fault-adjacent interior faces present → R-1205 satisfied;
   observed "interior central faces = 2/3" in the run). Re-entrant cycling None→Adjacent→AllContinuous→…→None
   leaves consistent state (base `BuildCentralFluxFaceSet_` clears + early-returns on None; the override
   clears the map first; IMPL-8 `map.size()==set.size()` holds at 0 for None). Minor: it sets `SetCflRkAware`
   on `b`/`w` per combo but only reads the pure factor (the ComputeMaxDt guard is deliberately not hit), so
   the cfl-aware bit is exercised only as a factor selector — fine for 4.1's intent.
7. **Test 4.2 rigor — genuine.** Forked `RunAbortsInChild` (serial, no MPI ⇒ `MFEM_ABORT`→`abort()`, fork
   safe) reproduces both aborts (bimat + scalar `Adjacent+ADER`). The parent sets the operator state
   (`SetMixedFluxMode(Adjacent); SetCflRkAware(false)`) before the fork; the child inherits it and the guard
   fires → non-zero exit → caught. Control cases are real: `Adjacent+RK` → no abort AND a separately-computed
   `dt_rk` is finite/positive; `None+ADER` → finite/positive dt (byte-exact contract). 5/5 reproduced.
8. **Test 15 update — correct.** `TestMatrixMixedFluxEnabled` on the 8-hex **fault-free** `CreateTestMesh`:
   `SetMixedFluxMode(AllContinuous)` on a `BimaterialWaveOperator` no longer aborts; populates a non-empty
   central set (observed 12 interior central faces — AllContinuous needs no fault) and
   `GetPerFaceCentralFlux().size()==GetCentralFluxFaceSet().size()`. The scalar control still accepts
   AllContinuous. Test 16 (R-006 precomputed-flux abort) is UNCHANGED and PASS. `wave_het` is the concrete
   `BimaterialWaveOperator` type so `GetPerFaceCentralFlux()` resolves (not a sliced base ref).
9. **Other — makefile + doc + FP.** `seas_test_bimaterial_mixed_flux_cfl` link target carries the correct
   operator-level deps (mirrors the parity target) and `test-bimaterial-mixed-flux-cfl` IS in the serial
   `test:` group (`makefile:3848`) — correct, since 4.1/4.2 are serial fork-based death tests, not MPI. The
   `SetCflRkAware` doc-comment (`wave_operator.hpp:212-219`) is updated to state both paths now require RK for
   mixed flux. No integer/FP issues in the factor lookup (all `real_t` literals).

## Bugs Found — Phase 4 (Rev I-4)

### P4-1: scalar `interior_flux=scalar` + `mixed_flux=adjacent` + ADER now ABORTS in the spatial driver — previously a supported, finite-dt path; shipped default TOMLs not migrated

**Status:** OPEN
**Severity:** HIGH (a previously-working, on-disk-default config now hard-aborts; deliberate per plan but incompletely landed)
**Files:**
- guard added: `dynamic/wave_operator.inl:5888` (scalar `ComputeMaxDt`)
- reached via: `drivers/spatial_dyn_driver.cpp:1245` (`SetMixedFluxMode`) then `:1903-1906` (`ComputeMaxDt`, `!is_rk` ⇒ no `SetCflRkAware`)
- stale TOMLs: `tpv102/configs/tpv102_spatial.toml:60-64` (`ader_order=2`, `mixed_flux="adjacent"`, `interior_flux="scalar"`, no `time_integrator`⇒ADER), `tpv102/configs/tpv102_spatial_p2.toml`, `_p3.toml`, `tpv205/configs/tpv205_spatial{,_p2,_p3}.toml`, `safs/project_7.0_alternative/config/spatial_friction_rate_state_safs_projected_stress_depthprofile.toml:81` (+ the other SAFS `mixed_flux="adjacent"` configs)
- stale claim: `drivers/spatial_dyn_driver.cpp:1899-1901` ("the ADER branch is byte-unchanged … keeps the pre-Phase-14 dt to the bit")

**Description:** The new scalar `ComputeMaxDt` guard `MFEM_VERIFY(mixed_flux_mode_==None || cfl_rk_aware_)` aborts on `scalar + Adjacent/AllContinuous + !cfl_rk_aware_`. Pre-Phase-4 that arm returned `0.9*cfl*h_min/cp` (Adjacent) / `0.4*…` (AllContinuous) — a finite, smaller dt that was the *intended* ADER guard band (the `cfl_rk_aware_ ? 0.6 : 0.9` ternary existed precisely to serve the ADER case). In the spatial driver, G1 (`spatial_friction.cpp:1162`) permits `scalar+adjacent`, G2 (`spatial_dyn_driver.cpp:813`) permits `ADER+scalar`, and `SetMixedFluxMode` runs before `ComputeMaxDt`, so this path is reachable and now aborts.

**Expected:** Either (a) the plan's intent — central/mixed flux is non-dissipative ⇒ ADER unstable ⇒ abort — applies and the shipped scalar-Adjacent-ADER TOMLs must be migrated (set `time_integrator="rk45"` or `mixed_flux="none"`) and the `:1899-1901` "byte-unchanged" claim corrected; OR (b) scalar **Adjacent** (≈5–10% central faces, dominantly upwind) was empirically ADER-stable with the 0.9 guard band and should keep returning a finite dt, with the hard abort reserved for **AllContinuous** (≈95% central) and the **bimaterial central** path. The reviewer cannot resolve (a) vs (b) from the source alone (the 0.9/0.4 factors are documented as "interim placeholders / starting calibration; production validation is the §14.5 Frontera run, not a unit test").

**Trigger:** `seas_spatial_dyn_driver --config tpv102/configs/tpv102_spatial.toml` (or any of the listed TOMLs) with no `--mixed-flux`/`--time-integrator` override; or `SAFS_MIXED_FLUX=adjacent` with the default ADER integrator. Aborts at `ComputeMaxDt` (`spatial_dyn_driver.cpp:1905-1906`).

**Suggested fix (pick one, with the maintainer):**
- If intent (a): migrate the shipped scalar-mixed TOMLs to `time_integrator="rk45"` (or `mixed_flux="none"`), correct the `spatial_dyn_driver.cpp:1899-1901` comment, and ADD a parser-or-construction-time diagnostic at the **driver** that names the offending config (so the user gets "set time_integrator=rk45 or mixed_flux=none" at config load, not an opaque `ComputeMaxDt` abort 1900 lines in). This also belongs with the Phase-5 G2 relax (same condition, same message — plan §Phase 4 req 3 / §Phase 5 G2).
- If intent (b): scope the abort to `AllContinuous || matrix-path`, i.e. let the scalar Adjacent ADER arm keep `0.9` (guard only `AllContinuous` and the bimaterial override). This preserves the gold scalar-Adjacent-ADER runs.

**(MODERATE) test:** a driver-level smoke (`Driver_ScalarAdjacentADER_Aborts` or `_RunsWithFiniteDt`, depending on (a)/(b)) on a tiny inline fixture that loads a `scalar+adjacent+ADER` config and asserts the chosen behaviour — closes the gap that NO test exercises the spatial-driver scalar-mixed-ADER reachability (Test 4.2 only probes the operator in isolation, where the guard is unambiguously correct).

### P4-2: stale comment — "Mixed flux is disabled on the matrix path … the switch is kept for parity" was deleted, but the doc on `MixedFluxCflFactor_` and the plan-era comments should name the lifecycle precondition

**Status:** OPEN
**Severity:** LOW (documentation)
**File:** `dynamic/wave_operator.hpp:786-795` (the `MixedFluxCflFactor_` doc-comment)
**Description:** The helper doc says it is "PURE (only the R-1600 fall-through aborts)". True, but it reads `mixed_flux_mode_` and `cfl_rk_aware_` — both **mutable** members set by `SetMixedFluxMode`/`SetCflRkAware`. "Pure" here means "no side effects", not "const-foldable"; a reader could mis-read it as order-independent of the setters. The factor genuinely depends on call order relative to those setters (which is exactly the P4-1 latent issue at the driver level).
**Expected:** A one-line note that the returned factor reflects the CURRENT `mixed_flux_mode_`/`cfl_rk_aware_`, i.e. callers must set both before `ComputeMaxDt`.
**Suggested fix:** Append to the doc-comment: "Reads the current `mixed_flux_mode_` and `cfl_rk_aware_`; set both (via `SetMixedFluxMode`/`SetCflRkAware`) before calling `ComputeMaxDt`."

### P4-3: redundant `return 1.0;` after a `[[noreturn]]` `MFEM_ABORT`; comment overstates its purpose

**Status:** OPEN
**Severity:** LOW (dead code / comment accuracy)
**File:** `dynamic/wave_operator.inl:5871`
**Description:** `MFEM_ABORT`→`mfem_error` is `[[noreturn]]` (`mfem/general/error.hpp:62`), so the `default:` arm never falls through. The trailing `return 1.0;   // unreachable (MFEM_ABORT above) — silences -Wreturn-type.` is dead. The `[[noreturn]]` attribute ALREADY silences `-Wreturn-type` for this path, so the comment's stated reason is inaccurate (the return is harmless belt-and-suspenders, not a necessity). Confirmed: `-O2` build emits no warning with the line, and would emit none without it.
**Expected:** Either drop the line or fix the comment to "defensive; the `default` arm is `[[noreturn]]`".
**Suggested fix:** Trivial — either delete `return 1.0;` or change the comment. Non-blocking.

## Verification of prior findings (carried into the Phase-4 tree)

| ID | Prior status | Phase-4 verification |
|----|--------------|----------------------|
| P3-4 (central+ADER guard interim, bimat) | OPEN (deferred to Phase 4) | **RESOLVED** — the interim bimat guard is now the permanent shared-condition guard on BOTH operators (`bimaterial_wave_operator.inl:803`, `wave_operator.inl:5888`); the divergent `0.9/0.4` switch is deleted; bimat now routes `this->MixedFluxCflFactor_()`. Test 4.2 reproduces the abort on both. (NEW-I.7 `BimatComputeMaxDt_CentralAderAborts` is now IMPLEMENTED as Test 4.2.) |
| P3-1 (shared central coverage) | OPEN | Still satisfied by `test_bimaterial_mixed_flux_shared` (14/14 np=2 reproduced); not a Phase-4 item. |
| P3-2 / P3-3 (LOW nits) | OPEN | Untouched by Phase 4; central 6/6 + dispatch 10/10 reproduced. |
| IMPL-8 lifecycle assert | verified Rev I-2 | Still holds; Test 15 + 4.1 exercise `map.size()==set.size()` across re-entrant mode changes (PASS). |

## Proposed Unit Tests (Implementation Rev I-4)

| ID | Test | Target | Validates | Priority |
|----|------|--------|-----------|----------|
| NEW-I.8 (NEW in Rev I-4) | `Driver_ScalarAdjacentADER_{Aborts\|RunsFiniteDt}` (tiny inline fixture) | `spatial_dyn_driver` config→`SetMixedFluxMode`→`ComputeMaxDt` path | a `scalar+adjacent+ADER` config reaches the intended behaviour (abort with a config-named message, OR finite de-rated dt) — closes the P4-1 reachability gap that Test 4.2 (operator-only) does not cover | **MEDIUM** (closes P4-1) |
| Test 4.1 | `ComputeMaxDt_FactorMatchesScalar_AllModes` | `MixedFluxCflFactor_` (both operators) | `[IMPLEMENTED]` `test_bimaterial_mixed_flux_cfl.cpp` — 6 combos × (scalar==expected, bimat==scalar), 12/12 PASS | HIGH (done) |
| Test 4.2 | `ComputeMaxDt_CentralPlusADER_Aborts` | `ComputeMaxDt` (both operators) | `[IMPLEMENTED]` `test_bimaterial_mixed_flux_cfl.cpp` — bimat+scalar abort on Adjacent+ADER; RK no-abort+finite dt; None+ADER finite dt, 5/5 PASS | HIGH (done) |
| Test 15 | `TestMatrixMixedFluxEnabled` | `BimaterialWaveOperator::SetMixedFluxMode` | `[IMPLEMENTED]` `test_wave_operator.cpp` — R-003 lifted, AllContinuous on fault-free mesh builds 12 central faces, map==set | HIGH (done) |

## Priority Summary (Implementation Rev I-4)

| ID | Severity | Status | Description |
|----|----------|--------|-------------|
| P4-1 | HIGH | OPEN | scalar `interior_flux=scalar`+`mixed_flux=adjacent`+ADER now aborts in `spatial_dyn_driver`; shipped default `*_spatial.toml` not migrated; `spatial_dyn_driver.cpp:1899-1901` "byte-unchanged" claim now false for that combo |
| P4-2 | LOW | OPEN | `MixedFluxCflFactor_` doc says "PURE" but reads mutable members; clarify the set-before-ComputeMaxDt precondition |
| P4-3 | LOW | OPEN | dead `return 1.0;` after `[[noreturn]]` `MFEM_ABORT`; comment overstates `-Wreturn-type` necessity |

**No CRITICAL.** One HIGH (P4-1) — a deliberate-but-incompletely-landed behaviour change that breaks a
previously-working, on-disk-default spatial-driver config (the operator-level guard is correct physics; the
gap is the un-migrated TOMLs + the missing driver-level diagnostic + the now-false "byte-unchanged" comment).
Two LOW doc/dead-code nits. The four byte-exact / single-source-of-truth crux items (acceptance #1 byte
identity; factor-vs-old-switch exactness; guard-condition identity across both operators; `MFEM_ABORT`
noreturn) are all clean and verified live.

### Test Coverage Summary (cumulative, after Phase 4)

| Suite | File | Tests | Phase-4 reproduced |
|-------|------|-------|--------------------|
| **CFL (4.1/4.2)** | `test_bimaterial_mixed_flux_cfl.cpp` | **17** | **17/17 serial** |
| Wave operator (Test 15/16) | `test_wave_operator.cpp` | 25 | 25/25 serial |
| Bimat parity (C-6c) | `test_bimaterial_wave_operator_parity.cpp` | 9 | 9/9 serial |
| phaseH parity (C-3) | `test_phaseh_wave_operator_constant_parity.cpp` | 19 / 52 | 19/19 serial, 52/52 np=4 |
| Bimat central primitive | `test_bimaterial_central_flux.cpp` | 6 | 6/6 serial |
| Bimat dispatch | `test_bimaterial_mixed_flux_dispatch.cpp` | 10 | 10/10 serial |
| Bimat shared central | `test_bimaterial_mixed_flux_shared.cpp` | 14 | 14/14 np=2 |
| Scalar adjacent dispatch | `test_mixed_flux_dispatch_adjacent.cpp` | 9 | 9/9 serial |
| **Driver scalar-mixed-ADER reachability** | — | **0** | **MISSING (P4-1 / NEW-I.8)** |

## Verdict (Phase 4 implementation)

**PASS WITH FIXES.** Phase 4 correctly (i) extracts `MixedFluxCflFactor_()` as the single source of
truth — byte-identical for the scalar operator (acceptance #1 verified live, C-3 bit-identical, C-6c parity),
exact-match to the pre-Phase-4 factors (Test 4.1, 12/12); (ii) deletes the bimaterial operator's divergent
`0.9/0.4` switch so both operators now route the SAME `cfl_rk_aware_`-gated factors; and (iii) installs the
identical central+ADER guard on BOTH `ComputeMaxDt`s (Test 4.2, both abort on Adjacent+ADER, RK/none do not).
Test 15 correctly flips from the R-003 abort to the R-003-lifted enable. All 161 reproduced assertions pass on
freshly recompiled objects (the fix report's claims hold). The assumed "≥3 bugs" did NOT surface as
correctness defects in the refactor itself; the one substantive finding (**P4-1, HIGH**) is a
deliberate-per-plan behaviour change (`scalar+adjacent+ADER` now aborts) that is **incompletely reconciled** —
the shipped default `*_spatial.toml` configs still carry `mixed_flux="adjacent"` under the ADER default, the
spatial driver has no config-named diagnostic for the new failure, and the `spatial_dyn_driver.cpp:1899-1901`
"ADER branch byte-unchanged" comment is now false for that combo. Resolve P4-1 with the maintainer (migrate
the TOMLs + add a driver diagnostic if the abort is intended, or scope the abort to AllContinuous/matrix if
scalar-Adjacent-ADER was a validated guard-band path), add NEW-I.8 (driver reachability test), and clean up
the two LOW nits (P4-2 doc, P4-3 dead return). The operator-level Phase-4 contract is sound; the gap is at the
driver/config boundary, which is Phase 5's surface.

---

# Implementation Review — Phase 5 (Rev I-6)

**Date:** 2026-06-04
**Scope (Phase 5 — driver + config wiring; PLAN lines ~494-586):**
`drivers/spatial_dyn_driver.cpp` (G2 relax, seam wiring, banner, P4 comment touch-ups),
`spatial/code/spatial_friction.{cpp,hpp}` (G1 relax, `seam_continuous` field + parse),
`tests/unit/test_tpv_config_parse.cpp` (Test 5.3). Sanity-checked the P4-fix comment
touch-ups in `dynamic/wave_operator.{hpp,inl}`.
**Reviewer:** chunhui-code-reviewer agent (adversarial, fresh rebuild — implementer pass
counts NOT trusted).

> **Rev numbering:** continues the implementation-review thread (I-1 Phase 1&2, I-2 Phase 3,
> I-4 Phase 4). This is **I-6**. Plan-defect history `BUG-1..BUG-26`, `IMPL-1..8`,
> `P3-1..4`, `P4-1..3` all remain intact above; new IDs here are `P5-1`, `P5-2`, ….

## Builds / Tests Re-Run (independent verification — fresh objects)

All builds in the worktree against the main checkout's `libmfem.a`
(`MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN …`); `.o` + binaries `rm`'d before each rebuild
(unreliable `.inl` dep tracking).

| Target | Result | Notes |
|--------|--------|-------|
| `seas_test_tpv_config_parse` (Phase 5) | **72 / 74 pass, 2 fail** | the 2 fails are `3 VS-border box rules present` (TPV205 :190, TPV104 :238) |
| `seas_test_tpv_config_parse` (worktree HEAD baseline) | **57 / 59 pass, 2 fail** | same 2 fails (`3 VS-border box rules`, TPV205 :166, TPV104 :214) — see byte-exactness ruling |
| `seas_test_bimaterial_mixed_flux_cfl` | 17 / 17 | Phase-4 CFL contract intact |
| `seas_test_wave_operator` | 25 / 25 | |
| `seas_test_bimaterial_wave_operator_parity` | 9 / 9 | |
| `seas_test_bimaterial_mixed_flux_dispatch` | 10 / 10 | |
| `seas_test_phaseh_wave_operator_constant_parity` (serial) | 19 / 19 | |
| `seas_test_phaseh_wave_operator_constant_parity` (`mpirun -np 4`) | 52 / 52 | |
| `seas_test_bimaterial_mixed_flux_shared` (`mpirun -np 2`) | 14 / 14 | incl. Test 3.4c depth-only + `seam_continuous=true` no-false-abort |
| `seas_spatial_dyn_driver` (heavy link) | **compiles + links clean** | no errors, no warnings; not run on a real mesh (per scope) |

**Byte-exactness verification method (acceptance #4/#5).** The G1 region was committed on the
*main* checkout's branch but is an *uncommitted* working-tree modification in this *worktree*
(branch `feat/mixed-flux-hetero-riemann`, HEAD `04674fc`). I established the baseline INSIDE the
worktree by extracting `git show HEAD:…` of the parser + test into the worktree files, rebuilding
(`rm` of `spatial_friction.o` + the test `.o`), running, then restoring the Phase-5 content from
backup. The baseline run (57/59) and the Phase-5 run (72/74) fail the **identical** two pre-existing
assertions (`3 VS-border box rules present`); Phase 5 adds 15 assertions, ALL passing, and changes
ZERO previously-asserted field values. The 2 failures are confirmed PRE-EXISTING and NOT caused by
Phase 5 (the failing `TEST_ASSERT(n_vs == 3, …)` lines are not in the Phase-5 diff). The implementer's
"2 pre-existing fails, unrelated" claim is **verified true**.

## Priority adversarial targets — verdicts

1. **G2 truth table (#1) — CORRECT.** `is_rk` (`:814-815`) is defined BEFORE `matrix_mixed`
   (`:816-818`), which is used at the verify (`:819`). `MFEM_VERIFY(!matrix_mixed || is_rk, …)`
   gives EXACTLY:
   | interior_flux | mixed_flux | integrator | `matrix_mixed` | verdict |
   |---|---|---|---|---|
   | matrix | adjacent/all_continuous | RK | true | **ALLOW** (`!true \|\| true`) |
   | matrix | adjacent/all_continuous | ADER | true | **ABORT** (`!true \|\| false`) |
   | matrix | none | RK | false | ALLOW |
   | matrix | none | ADER | false | ALLOW |
   | scalar | anything | anything | false | ALLOW (G2 no longer constrains scalar) |

   This matches plan req 2 exactly. The OLD guard `!is_rk || interior_flux==Scalar` blocked
   `matrix + RK` for ALL mixed_flux (incl. `matrix+none+RK`); the new guard **intentionally**
   re-admits `matrix+none+RK` (plan req 2: "RK on the bimaterial upwind operator is harmless").
   `scalar+RK` was allowed by BOTH old and new (no regression). `matrix_mixed` is in scope and NOT
   shadowed at the req-4 banner (`:1263`). The CLI `--time-integrator` (`:665-675`) and
   `--mixed-flux` (`:679`) overrides BOTH precede `:814`, so `is_rk`/`matrix_mixed` read the final
   post-override values. G2's message is the actionable variant ("requires an RK integrator. Set
   --time-integrator rk4|rk45, or use mixed_flux=none under ADER"); it is consistent in CONDITION
   with the `ComputeMaxDt` guard (`bimaterial_wave_operator.inl:803`, `wave_operator.inl:5888`) —
   both abort iff `mixed && !is_rk`/`mixed && !cfl_rk_aware_`. (See P5-2 for a cosmetic message
   token mismatch in the operator copy.)

2. **G1 removal (#2) — clean, did NOT over- or under-relax.** The deleted
   `MFEM_VERIFY(interior_flux==Scalar || mixed_flux=="none")` was a standalone statement; its removal
   left no dangling `else`/brace (`spatial_friction.cpp:1159-1170` reads cleanly). The replacement
   comment is accurate. KEPT guards survive and were re-verified: (a) the `[numerics].mixed_flux`
   VALUE validator `none|adjacent|all_continuous` (`:1152-1157`) still rejects typos UPSTREAM of G1;
   (b) the `interior_flux="matrix"` ⇒ non-Constant material guard (`:1707-1712`) intact (TPV31 is
   `depth_profile_1d`, passes). The "companion guard NOTE" (`:1168-1170`) is preserved, not orphaned.

3. **Exhaustiveness / 4th-guard hunt (#3) — NO 4th blocking guard. CONFIRMED 3.** Independent
   re-sweep of `drivers/` + `spatial/code/` + `dynamic/` for `MFEM_VERIFY`/`MFEM_ABORT` gating on
   `interior_flux`/`InteriorFlux`/`mixed_flux`/`MixedFluxMode`/`cfl_rk_aware_`/`is_rk`:
   - Driver: ONLY G2 (`:819`) gates on `(Matrix && mixed!=none)`. No other `Matrix`/`is_rk` abort.
   - `spatial_friction.cpp`: G1 (lifted), the value validator (`:1152`, keep), the non-Constant guard
     (`:1707`, keep). No integrator gate.
   - `bimaterial_wave_operator.inl`: `:803` central+ADER (`ComputeMaxDt`; PASSES on RK because the
     driver calls `SetCflRkAware(true)` at `:1933` BEFORE `ComputeMaxDt` at `:1934`); `:604` the
     null-transform-on-central-shared invariant (IMPL-4, not flux-mode gating); `:635` the
     `seam_continuous_` gate (PASSES because the driver wires it true). `:877` the IMPL-8 size assert
     (holds). The R-003 abort in the override is GONE (delegates to base).
   - `wave_operator.inl`: `:1615` R-1203 precomputed-flux exclusion fires ONLY if
     `use_precomputed_face_fluxes_` — the driver NEVER calls `UsePrecomputedFaceFluxes` (grep: zero
     hits), so it is false on this path (keep, defensive); `:5888` the scalar central+ADER guard
     (same `cfl_rk_aware_` gate). **A `matrix + adjacent + RK + seam_continuous=true` config reaches
     the time loop with no abort.** No CRITICAL "feature can't run" guard remains.

4. **`seam_continuous` wiring (#4) — CORRECT.** `SetSeamContinuous` is a public inline setter on
   `BimaterialWaveOperator` (`bimaterial_wave_operator.hpp:147`, `void SetSeamContinuous(bool v) {
   seam_continuous_ = v; }`). The driver calls it (`:1149-1150`) immediately after the
   `BimaterialWaveOperator` ctor (`:1143`) and BEFORE `SetMixedFluxMode` (`:1259`); `SetMixedFluxMode`
   (`inl:857`) → `BuildPerFaceCentralFluxMatrices_()` (`inl:869`) → reads `seam_continuous_`
   (`inl:636`). Ordering is correct. The `static_cast<BimaterialWaveOperator<ParMesh>&>(*wave_ptr)`
   is SAFE — it is inside the `else // Matrix` branch (`:1118-1151`) where `wave_ptr` was just
   constructed as a `BimaterialWaveOperator`. It runs for `matrix+none` too (harmless: the central
   set is empty, `BuildPerFaceCentralFluxMatrices_` builds nothing, the `:635` guard never fires).
   Default `false` ⇒ byte-exact for every matrix config that omits the key (TPV31 does); regression
   test 5.3c confirms `tpv31.seam_continuous == false`. Confirmed live by `test_bimaterial_mixed_flux_shared`
   Test 3.4c (np=2, depth-only + `seam_continuous=true`, no false abort, 14/14).

5. **Byte-exactness of existing configs (#5) — VERIFIED.** See the verification-method paragraph
   above. The `seam_continuous` parse is purely additive (`toml_bool(m,"seam_continuous",false)`),
   defaults match (`interior_flux=scalar`, `mixed_flux=none`, `time_integrator=ADER`,
   `seam_continuous=false` — req 7/8 confirmed at `spatial_friction.cpp:1123/1096/1133` +
   `.hpp:133/139/592`). The 2 pre-existing `VS-border` failures reproduce IDENTICALLY at the
   worktree HEAD baseline (57/59) — NOT a Phase-5 regression.

6. **Test 5.3 rigor (#6) — adequate; one cleanup nit (P5-3).** 5.3a genuinely exercises G1 (a
   `matrix + adjacent` TOML now `LoadSpatialFrictionConfig`s without abort — it WOULD have aborted
   at `spatial_friction.cpp:1162` pre-Phase-5). 5.3b round-trips `seam_continuous=true`. The string
   mutation is robust: `r1` flips `mixed_flux = "none"` → `"adjacent"`; `r2` is newline-anchored
   (`\ndepth_axis = "z"\n`) which uniquely matches the standalone `[material]` config line 98 and
   NOT the comment at line 25 (`In canonical: depth_axis = "z" with depth…`, not newline-bounded);
   the `r1`/`r2` bool asserts fail loud if a future TOML edit breaks the anchor (no silent no-op).
   The temp file `<path>.matrixadj.TESTGEN.toml` is `std::remove`'d AFTER the parse — see **P5-3**
   for the assertion/abort-leak nit.

7. **Tests 5.1/5.2 deferral (#7) — ACCEPTABLE, but a lighter faithful test was available.** See
   the explicit ruling in the Verdict.

8. **Defaults (#8) + P4-1 driver comment (#8) — accurate.** All four defaults unchanged (above).
   The `SetCflRkAware` doc at `wave_operator.hpp:212-219` was correctly de-staled: the prior
   "Only the scalar path is affected (RK forbids the matrix/bimaterial path)" was FALSE post-Phase-5
   and is replaced by "Both the scalar AND the matrix/bimaterial path use this." The driver
   `:1920-1932` comment correctly scopes the byte-exactness claim to `mixed_flux=none` and documents
   the `mixed_flux!=none` ADER abort. **Accurate.**

## Bugs Found — Phase 5 (Rev I-6)

### P5-1 — Shipped `*_spatial.toml` still abort under their on-disk defaults (P4-1 carried to the driver/config boundary, now Phase-5-owned)

**Status:** OPEN
**File:** `tpv102/configs/tpv102_spatial.toml:61,64`, `tpv104/configs/tpv104_spatial.toml:60,63`,
`tpv205/configs/tpv205_spatial.toml:90,93` (+ driver `spatial_dyn_driver.cpp:1934`)
**Severity:** MEDIUM (was HIGH as P4-1; Phase 5 documents the workaround but does not resolve it)
**Description:** All three shipped spatial TOMLs carry `interior_flux="scalar"` + `mixed_flux="adjacent"`
with NO `time_integrator` key ⇒ default ADER. Phase 4+5 make `scalar+adjacent+ADER` ABORT at
`ComputeMaxDt` (`wave_operator.inl:5888`, reached at driver `:1934`). A bare `seas_spatial_dyn_driver
--config tpv102_spatial.toml` (no `--time-integrator`/`--mixed-flux` CLI override) therefore aborts
on its own default config. Phase 5 did NOT migrate these TOMLs (G2 at `:819` does not catch this — it
gates only the `matrix` path; the scalar+adjacent+ADER abort is the `ComputeMaxDt` guard, which fires
later). This is the same landmine as P4-1; Phase 5 OWNS the config surface and explicitly references it
in the `:1926-1932` comment ("real jobs CLI-override to rk45 or --mixed-flux none") but leaves the
on-disk default broken.
**Expected:** Either (a) migrate the shipped TOMLs to `time_integrator="rk45"` (matching how "real
jobs" run), or (b) add a config-named diagnostic at the driver so the abort points at the TOML +
remediation, or (c) confirm-with-maintainer that the abort is intended and the TOMLs are "templates
that MUST be CLI-overridden" (then document that contract in the TOML header comments).
**Trigger:** `seas_spatial_dyn_driver --config <any of the three>_spatial.toml` with no integrator/flux
CLI flag.
**Suggested fix:** Lowest-risk: add `time_integrator = "rk45"` to the three TOMLs (the gold runs already
use rk45 + adjacent per the tpv205 comment), OR add a `[meta]`-level note + a driver hint string. Confirm
with maintainer whether ADER+adjacent was ever a validated finite-dt path (P4-1's open question).

### P5-2 — `ComputeMaxDt` central+ADER message advertises the non-existent token `"allcontinuous"` (no underscore)

**Status:** OPEN
**File:** `dynamic/bimaterial_wave_operator.inl:~810` and `dynamic/wave_operator.inl:~5895` (the two
`ComputeMaxDt` central+ADER guard messages)
**Severity:** LOW (cosmetic; misleading remediation text)
**Description:** Both `ComputeMaxDt` abort messages suggest `interior_flux="matrix" + mixed_flux=
"adjacent"/"allcontinuous"`, but the parser's value validator (`spatial_friction.cpp:1152`) and
`ParseMixedFlux` (`spatial_dyn_driver.cpp:175`) only accept `"all_continuous"` (WITH underscore). A
user copy-pasting `"allcontinuous"` from the abort message hits the value-validator abort. Not a
functional defect (the abort already fired) but the remediation text is wrong. The Phase-5 G2 driver
message (`:819-824`) is fine — it does not enumerate the modes.
**Expected:** `"all_continuous"` to match the accepted config token.
**Trigger:** User reads the abort message and types the suggested token.
**Suggested fix:** s/`allcontinuous`/`all_continuous`/ in both `ComputeMaxDt` messages. (These are
Phase-4 files, but the typo is in the user-facing remediation path that Phase 5's feature surfaces.)

### P5-3 — Test 5.3a leaves `tpv31.toml.matrixadj.TESTGEN.toml` litter if the parse aborts (no RAII / no cleanup on the failure path)

**Status:** OPEN
**File:** `tests/unit/test_tpv_config_parse.cpp` (`T_TPV31_MatrixAdjacentSeam`, the
`WriteFile(tmp,…); LoadSpatialFrictionConfig(tmp); std::remove(tmp)` sequence)
**Severity:** LOW (test hygiene; only triggers on a regression)
**Description:** `std::remove(tmp.c_str())` runs only AFTER `LoadSpatialFrictionConfig(tmp)` returns.
If a future regression re-introduces an abort in the parser on `matrix+adjacent` (exactly the BUG-21
condition this test guards), `MFEM_VERIFY`/`MFEM_ABORT` calls `mfem::mfem_error` → `abort()`, so the
`std::remove` never runs and a `tpv31.toml.matrixadj.TESTGEN.toml` is left next to the real config in
the source tree. Likewise a failing `TEST_ASSERT` after the write (none currently sit between
`WriteFile` and `remove`, but adding one would leak). The temp name is distinct enough not to clobber
a real file, and BUG-21 is fixed so the abort won't fire today — but the test's own purpose is to
catch a re-introduced abort, and on that exact regression it litters the repo.
**Expected:** Temp artifact removed on every exit path (including abort/throw), or written to a
process-temp dir (`std::tmpnam`/`TMPDIR`) rather than next to the source config.
**Trigger:** A regression re-adds the G1 abort (or any parser abort) on `matrix+adjacent`.
**Suggested fix:** Write the temp copy under `std::getenv("TMPDIR")` (fallback `/tmp`) so litter never
lands in the tracked source tree; or wrap the temp file in a small RAII guard whose dtor `std::remove`s
it (dtor still runs on a C++ exception, though not on `abort()` — the TMPDIR relocation is the robust
fix since MFEM aborts rather than throws).

### P5-4 — `--mixed-flux` CLI override (`:679`) is not re-validated before G2 / `ParseMixedFlux`

**Status:** OPEN — DEFERRED (pre-existing; not introduced by Phase 5)
**File:** `drivers/spatial_dyn_driver.cpp:679` (`cfg.numerics.mixed_flux = cli_mixed_flux;`)
**Severity:** LOW
**Description:** A `--mixed-flux <garbage>` CLI override writes the raw string into
`cfg.numerics.mixed_flux` AFTER the parser's value validator (`spatial_friction.cpp:1152`) has run, so
`matrix_mixed` (`:816`) would treat `matrix + garbage` as `true` and the run would proceed past G2
until `ParseMixedFlux("garbage")` aborts at `:1259` with the (correct) accepted-values message.
Net effect is still a loud abort, just a few lines later than ideal, and with G2's banner possibly not
printed. Pre-existing (the override predates Phase 5); flagged for completeness since Phase 5 made the
`matrix + mixed` path live and `matrix_mixed` now reads this override.
**Suggested fix:** Re-run the `none|adjacent|all_continuous` validator on `cli_mixed_flux` at `:679`
(or route the CLI override through `ParseMixedFlux` immediately to fail at the override site).

## Verification of prior findings (carried into the Phase-5 tree)

| Prior finding | Status after Phase 5 |
|---|---|
| P4-1 (scalar+adjacent+ADER aborts; shipped TOMLs not migrated) | **STILL OPEN — re-scoped as P5-1.** Phase 5 documents the workaround (`:1926-1932`) but does not migrate the on-disk default TOMLs; the abort is the `ComputeMaxDt` scalar guard, untouched by G2. |
| P4-2 (`MixedFluxCflFactor_` "PURE" doc) | OPEN — Phase-5 hpp touch-up (`:782-797`) clarifies it reads CURRENT `mixed_flux_mode_`/`cfl_rk_aware_` set by the setters; the precondition is now documented. Effectively addressed by the doc edit; close on confirmation. |
| P4-3 (dead `return 1.0` after `[[noreturn]]`) | OPEN — untouched by Phase 5 (LOW). |
| P3-1 (shared central coverage) | Satisfied — `test_bimaterial_mixed_flux_shared` 14/14 np=2 reproduced (incl. 3.4c). |
| P3-4 (central+ADER guard) | RESOLVED in Phase 4; re-confirmed here (both `ComputeMaxDt`s gate on `cfl_rk_aware_`; driver sets it true before `ComputeMaxDt` on RK). |
| G1/G2/G3 inventory (plan headline "exactly 3") | **VERIFIED** by independent sweep — no 4th blocking guard. |

## Proposed Unit Tests (Implementation Rev I-6)

| Test ID | Name | Target | Validates | Priority |
|---|---|---|---|---|
| NEW-I.9 (NEW in Rev I-6) | `Driver_G2Predicate_TruthTable` | extract the G2 predicate `(interior_flux==Matrix && mixed_flux!="none") && !is_rk` into a small free `static bool spatial::MatrixMixedRequiresRk(const SpatialFrictionConfig&)` helper and table-test all 5 rows | the G2 driver-main guard (`:819`) without needing a mesh fixture — closes the Test 5.1/5.2 deferral gap with NO `--dry-run`/mesh cost | **MEDIUM** (faithful, mesh-free coverage of the driver guard) |
| NEW-I.10 (NEW in Rev I-6) | `ConfigParse_MixedFluxTypo_Rejected` | `LoadSpatialFrictionConfig` on `mixed_flux="adjcent"` | the value validator (`spatial_friction.cpp:1152`) still rejects typos after the G1 relaxation (guards against an over-relaxation regression) | LOW |
| NEW-I.11 (NEW in Rev I-6) | `Driver_ShippedSpatialToml_ADERDefault_Aborts` (tiny inline fixture or the NEW-I.9 helper) | the P5-1 landmine | a default-integrator load of a `scalar+adjacent` config reaches the intended abort with a config-named message (drives the P5-1 fix) | MEDIUM (closes P5-1) |
| NEW-I.12 (NEW in Rev I-6) | `Test5.3a_TempFile_NoLeakOnAbort` | `T_TPV31_MatrixAdjacentSeam` cleanup | the TESTGEN temp is written under TMPDIR (or removed on all paths) so a parser-abort regression does not litter the source tree | LOW (closes P5-3) |

## Priority Summary (Implementation Rev I-6)

| ID | Severity | Status | Description |
|----|----------|--------|-------------|
| P5-1 | MEDIUM | OPEN | Shipped `tpv102/104/205_spatial.toml` (scalar+adjacent, default ADER) abort on their own on-disk defaults at `ComputeMaxDt`; Phase 5 documents the CLI workaround but does not migrate them (P4-1 carried to the config surface Phase 5 owns) |
| P5-2 | LOW | OPEN | `ComputeMaxDt` central+ADER messages advertise `"allcontinuous"` (no underscore); accepted token is `"all_continuous"` |
| P5-3 | LOW | OPEN | Test 5.3a leaves a `*.matrixadj.TESTGEN.toml` in the source tree if the parser aborts (cleanup runs only post-parse); relocate to TMPDIR |
| P5-4 | LOW | OPEN (DEFERRED, pre-existing) | `--mixed-flux` CLI override not re-validated before G2/`ParseMixedFlux`; garbage value aborts a few lines later than ideal |

**No CRITICAL. No HIGH.** The four Phase-5 requirements (G1 relax, G2 relax, `seam_continuous` field +
wiring, banner) are implemented correctly: the G2 truth table is exact, the G1 removal is clean and
neither over- nor under-relaxes, the independent guard sweep confirms the plan's headline "exactly 3
guards" (no 4th blocker — the feature CAN run), the `seam_continuous` wiring is ordered + type-safe +
byte-exact-by-default, and byte-exactness of every existing config is verified by a fresh
baseline-vs-Phase-5 rebuild (identical 2 pre-existing failures). The four findings are MEDIUM/LOW: the
top one (P5-1) is the P4-1 config landmine surfacing on the driver/config boundary Phase 5 owns; the
rest are cosmetic message/test-hygiene nits.

### Test Coverage Summary (cumulative, after Phase 5)

| Suite | File Under Test | Tests reproduced | Key items |
|-------|-----------------|------------------|-----------|
| config-parse | `spatial_friction.cpp` | 72/74 (2 pre-existing fails) | Test 5.3a/b/c (G1 relax + seam round-trip + TPV31 regression) — all 15 new pass |
| bimaterial CFL | `wave_operator.inl` / `bimaterial_wave_operator.inl` | 17/17 | Phase-4 factor + central+ADER guard |
| wave operator | `wave_operator.inl` | 25/25 | |
| bimaterial parity | `bimaterial_wave_operator.inl` | 9/9 | |
| mixed-flux dispatch | `bimaterial_wave_operator.inl` | 10/10 | |
| phaseh const parity | `wave_operator.inl` | 19/19 + 52/52 (np=4) | |
| bimaterial shared (np=2) | `bimaterial_wave_operator.inl` | 14/14 | Test 3.4c seam_continuous wiring (Phase-5-relevant) |
| spatial driver | `spatial_dyn_driver.cpp` | compiles+links | G1/G2/seam edits build clean |
| **Driver-main G2 guard reachability (5.1/5.2)** | `spatial_dyn_driver.cpp:819` | **0 (deferred)** | covered indirectly by Test 4.2 (operator) + 5.3 (config); see deferral ruling — NEW-I.9 recommended |

## Verdict (Phase 5 implementation)

**PASS WITH FIXES.** Phase 5 lands the driver + config wiring correctly. The adversarial priority
targets all clear: (#1) the G2 truth table is EXACTLY right and `is_rk`/`matrix_mixed` are ordered +
in-scope + not shadowed; (#2) G1 removal is surgically clean (no orphaned brace/else, kept guards
survive, no over/under-relaxation); (#3) the independent re-sweep of `drivers/` + `spatial/code/` +
`dynamic/` confirms the plan's headline claim — EXACTLY three guards (G1/G2/G3) ever blocked
matrix+mixed and all three are relaxed; **NO 4th blocking guard exists**, so a
`matrix+adjacent+RK+seam_continuous=true` config reaches the time loop with no abort; (#4) the
`seam_continuous` setter is public, type-safe (`static_cast` is in the guaranteed-bimaterial branch),
called BEFORE `SetMixedFluxMode`, harmless for `matrix+none`, and byte-exact by default; (#5)
byte-exactness is VERIFIED by a fresh worktree baseline rebuild (57/59) vs Phase-5 (72/74) — the
identical 2 pre-existing `VS-border` failures, zero new failures, zero changed field values. The
assumed "≥3 bugs" did NOT surface as correctness defects in the Phase-5 wiring; the four findings are
one MEDIUM config landmine (P5-1, the P4-1 carry-over) and three LOW nits (P5-2 message typo, P5-3
test litter, P5-4 CLI re-validation).

**Ruling on the Test 5.1/5.2 deferral (#7): ACCEPTABLE, but a lighter faithful test SHOULD have been
added.** The implementer's `--dry-run` rationale is partly inaccurate: `--dry-run` exits at
`spatial_dyn_driver.cpp:2194`, which is AFTER G2 (`:819`), construction (`:1143`), `SetMixedFluxMode`
(`:1259`), AND `ComputeMaxDt` (`:1934`) — so a dry-run WOULD exercise both G2 and the central+ADER
abort. The genuine blocker is that a dry-run needs a real (tiny, correctly-tagged) mesh to construct
the operator (the driver's hardcoded boundary-attribute handling), which is a real fixture cost; on
THAT basis the end-to-end driver test is reasonably deferred to Phase 7. HOWEVER, the G2 guard itself
(`:816-824`) is a PURE boolean over three config fields — it has no mesh dependency. Extracting that
predicate into a free testable helper (NEW-I.9, `spatial::MatrixMixedRequiresRk(cfg)`) and
table-testing the 5 rows would faithfully cover the exact driver-main guard with ZERO mesh/dry-run
cost, and would have closed 5.1/5.2 properly. The current indirect coverage (Test 4.2 operator-level
abort + Test 5.3 config-level G1) leaves the driver-main G2 boolean itself untested in isolation.
Deferral is acceptable for *merge* (the guard is trivially correct by inspection and physics is
covered by 4.2), but NEW-I.9 should be added before relying on the driver guard in CI. The same applies
to the P5-1 landmine (NEW-I.11).


---

# Implementation Review — Phase 6 (Rev I-8)

**Date:** 2026-06-04
**Scope:** Phase 6 code — the homogeneous-equivalence GATE
(`tests/unit/test_bimaterial_mixed_flux_homog_equivalence.cpp`, NEW) + the four-edit
`Makefile` registration of all Phase-1/3/4/6 unit targets.
**Plan:** `PLAN_mixed_flux_hetero_riemann.md` Phase 6 (lines ~589-641).
**Reviewer:** chunhui-code-reviewer agent (adversarial; rebuilt fresh, re-ran, did NOT trust reported counts).

## Build / test results reproduced (fresh objects — `rm`'d `.o` + binary before each build)

| Target | Built | Result reproduced |
|--------|-------|-------------------|
| `seas_test_bimaterial_mixed_flux_homog_equivalence` (the GATE, Test 6.1) | clean | **18 / 18** (Adjacent 9 + AllContinuous 9); central faces = 2, fault QPs = 8 for BOTH modes |
| `seas_test_bimaterial_mixed_flux_dispatch` (Test 3.1-3.3, incl **3.2** het single-valuedness) | clean | **10 / 10** (Test 3.2 PASS — acceptance #3) |
| `seas_test_bimaterial_wave_operator_parity` (C-6) | clean | **9 / 9** |
| `seas_test_bimaterial_central_flux` (Test 1.1-1.3) | clean | **6 / 6** (Test 1.1 homogeneous == scalar central: **1.470e-15**, tol 1e-11) |
| `seas_test_bimaterial_mixed_flux_cfl` (Test 4.1-4.2) | clean | **17 / 17** |
| `seas_test_phaseh_wave_operator_constant_parity` (serial) | clean | **19 / 19** |
| `seas_test_phaseh_wave_operator_constant_parity` (np=4) | clean | **52 / 52** |
| `seas_test_bimaterial_mixed_flux_shared` (Test 3.4, np=2) | clean | **14 / 14** |

All eight suites green. The implementer's "18/18" gate count is HONEST and reproduces.

## Priority adversarial targets — verdicts

**#1 — Does the gate exercise the central flux, or pass vacuously? — VACUOUS-GATE RISK CLEARED.**
The fixture's `GetCentralFluxFaceSet()` contains **exactly the 2 x=1 interior faces** and NOT the y=1
fault faces, verified three ways: (a) the operator log prints `interior central faces = 2, shared
central faces = 0` and `fault QPs = 8` (= 2 fault faces × 4 quad QPs) — fault and central sets are
disjoint and correctly sized; (b) the mesh topology forces it — the 2×2×1 hex block has exactly two
non-fault two-sided interior faces (the x=1 plane separating ix=0/ix=1 in each y-row), and the fault
tag fires only on faces with face-centroid |cy − 1| < 1e-9 (the y=1 plane), whose centroids are y=0.5
/ y=1.5 for x=1 faces ⇒ never mis-tagged; (c) the central matrices are actually applied in `Mult`
(see #2). The implementer's `!empty()` assert is NECESSARY-but-WEAK on its own (a non-empty set could
in principle hold the wrong faces), but it is BACKED by the `size()==2` print + the dispatch-teeth
experiment (#2), so the wrong-face failure mode is foreclosed on this fixture. No bug, but see P6-3
for the assertion hardening.

**#2 — Is the comparison non-trivial / does the gate have TEETH? — YES, PROVEN.**
Scratch experiment (reverted; tree restored to md5 `58d91bef…`): perturbed ONE central matrix factor
in `godunov_flux_bimaterial.cpp:380` (`centralNbr *= 0.5` → `*= 0.5000001`, a ~2e-7 relative
perturbation that lives ONLY on the central-flux path). Rebuilt the bimaterial object + gate;
result: **12/18 (6 FAILED)**. The homogeneous-equivalence relative diff jumped to **3.00e-07** (VX
strike) and **1.00e-07** (VZ anti-plane), driving 64-128 of 288 DOFs over the 1e-9 tolerance. This
proves (i) the central dispatch is genuinely reached inside `Mult`, (ii) the central contribution is
non-negligible relative to the field scale (a 2e-7 matrix perturbation produces a 1-3e-7 field-relative
deviation — five orders above the 1e-9 floor), and (iii) the gate WOULD catch a defect in the
bimaterial central assembly (wrong ½ factor, wrong rotation, missed/extra central face, per-side sign).
The differing VX (3e-7) vs VZ (1e-7) deviations also confirm the three states exercise the central
path to DIFFERENT degrees (VZ anti-plane couples less). The gate is meaningful.

**#3 — Adjacent vs AllContinuous DISTINCT on this fixture? — NO; coverage gap. See P6-1 (MODERATE).**
On the 2×2×1 mesh EVERY element touches the y=1 fault plane, so `E_fault_adj` = all 4 elements ⇒
Adjacent's element-walk visits every interior non-fault face ⇒ **Adjacent's central set ≡
AllContinuous's central set (both = the 2 x=1 faces)**, confirmed by identical `central faces = 2`
for both modes. The "AllContinuous" run is therefore a DUPLICATE of the Adjacent run on this fixture;
it does NOT exercise the all-interior-faces semantics (`AllContinuous`'s distinct
`for(f) is_interior && !is_fault` loop at `wave_operator.inl:1748-1758`) that distinguishes it from
Adjacent's `E_fault_adj` element-walk (`:1764-1806`). Note both paths ARE covered by the upstream
`seas_test_bimaterial_mixed_flux_dispatch` (5-hex fixture, Test 3.1) — so this is a Phase-6 GATE
coverage gap, not a project-wide one. Severity MODERATE (the gate over-claims AllContinuous coverage),
not CRITICAL (the distinguishing logic is tested elsewhere).

**#4 — Fixture correctness. — CORRECT.**
The hex vertex order `{(ix,iy,iz),(ix+1,iy,iz),(ix+1,iy+1,iz),(ix,iy+1,iz), z+1 layer…}` is the
standard MFEM hex (bottom face CCW then top face CCW). No negative-Jacobian / mesh-validity / orient
warnings emitted (grep clean); the central matrices applied without `ApplyPerFaceFlux` SPD/abort,
and the run is finite. `mesh.FinalizeTopology()` (×2), `mesh.Finalize()`, `mesh.SetAttributes()` all
present. The fault tag hits exactly the 2 y=1 faces (fault QPs = 8), x=1 faces left interior (central
= 2). `SetupFault` is byte-identical to the parity test's but is geometry-agnostic via
`IntRules.Get(ftr->GetGeometryType(), 2*order)` — it correctly returns `Geometry::SQUARE` (4 QPs) for
hex faces (vs the parity test's `Geometry::TRIANGLE`); nqp_per_face = 4 reproduced. One pre-existing
inherited pattern (NOT a Phase-6 bug): `SetupFault` reads `nqp_per_face` from `int_faces[0]` only and
assumes all fault faces share it — fine on a uniform hex mesh; flagged as P6-4 (LOW, doc/hardening).

**#5 — 1e-9 tolerance + metric. — APPROPRIATE; comfortably achievable.**
Field-scale `MaxRelDiff` is the SAME metric the C-6 parity test justified, applied here to the
`Mult` RHS (`k_s` vs `k_b`) rather than near-equilibrium fault DOFs — appropriate because the two
operators run the IDENTICAL state and the only intended difference is the central-vs-Godunov face
assembly. The 6-sig-fig material (lam=mu=32e9, rho=2670) keeps the GodunovFluxPool dedup a no-op
(pool == scalar `flux_` bit-for-bit), so the homogeneous bimaterial central `½·T·(AxPlus+AxMinus)·Tinv`
equals the scalar `GodunovFlux::Central` to LU rounding. Verified the floor INSIDE `Mult`: the gate's
observed unperturbed diff is < 1e-9 (n_over = 0) and the standalone Test 1.1 measures **1.470e-15**
this run. (Minor doc note: the gate comment / task cite "8.6e-16"; the reproduced figure is 1.47e-15
— both are sub-1e-9 LU-rounding floors, harmless, but the literal cited constant is stale; P6-4.)

**#6 — Does the NOTE correctly state necessary-not-sufficient (BUG-3 masked by homogeneity)? — YES.**
The gate header (`:13-20`) and the inline comment (`:238-239`) correctly state the gate is NECESSARY
but NOT SUFFICIENT (homogeneous ⇒ A_e1==A_e2, so a BUG-3 per-side swap would also pass here) and name
Test 3.2 as the sufficient guard. Test 3.2 reproduces PASS (single-valuedness under genuine contrast,
F_h_e1 == F_h_e2 to ≤1e-12 rel) — acceptance #3 satisfied.

**#7 — Makefile. — CORRECT.**
The operator-level dep list (`WAVE_OPERATOR_OBJ`, `GODUNOV_FLUX_BIMATERIAL_OBJ`, `GODUNOV_FLUX_OBJ`,
`GODUNOV_FLUX_POOL_OBJ`, `PML_LAYER_OBJ`, `PRECOMPUTED_FACE_FLUXES_OBJ`, `HETEROGENEOUS_MATERIAL_OBJ`,
`FAULT_FACE_FLUX_OBJ`, `FRICTION_SOLVER_OBJ`, `SPATIAL_FRICTION_OBJ`) mirrors the parity/cfl/dispatch
targets exactly. All four coordinated edits present (SRC var, OBJ var, object-compile rule, link
target + `test-…` alias). The gate is in the SERIAL `test:` group (`Makefile:3849`); the np=2 shared
target is correctly EXCLUDED from `test:` (it needs `mpirun -np 2`, `:5354`). Builds clean, links,
runs. (Plan §Phase-6 step 2 prescribes a leaner dep list — `WAVE_OPERATOR + GODUNOV_FLUX_BIMATERIAL +
PRECOMPUTED_FACE_FLUXES + GODUNOV_FLUX + PML_LAYER + FAULT_FACE_FLUX + FRICTION_SOLVER`; the
implementer ALSO links `GODUNOV_FLUX_POOL`, `HETEROGENEOUS_MATERIAL`, `SPATIAL_FRICTION`. These are
real transitive deps of the operator-level objects, not bloat — the link succeeds and `-dead_strip_dylibs`
is set; matching the proven parity-target list over the plan's slightly-stale list is the correct call.)

**#8 — Acceptance #4 deferral (np≤4 het DRIVER smoke, ≥5 RK4 steps). — DEFERRAL IS ACCEPTABLE; ruling below.**
The implementer's reasoning — "a hand-rolled RK4 on Q alone would over-advance the fault state because
the friction solve is inside `Mult`, so faithful coupling needs `AdvanceRKCoupledLSW_Spatial`" — is
SOUND and load-bearing. Verified against `rk_time_stepper.hpp:336-455`: the coupled LSW stepper stages
the per-DOF slip δ into `dof_data[m].slip{1,2}` BEFORE each `wave.Mult` (`:429-441`), so the
slip-stateless `EvaluateLSW` reads the stage-local δ; a naive RK4-on-Q would skip that stage-local
slip write and the friction solve would see stale slip — non-faithful (wrong-physics) dynamics. So a
faithful local smoke genuinely cannot be a hand-rolled RK4. HOWEVER, the same stepper is templated on
`MeshType` and explicitly "serves both the serial-`Mesh` unit tests and the production `ParMesh`
driver" (`:367-368`) — so a FAITHFUL local smoke IS feasible without the full driver by calling
`AdvanceRKCoupledLSW_Spatial` directly on a 2-element het fixture (≥5 steps, assert finite δ/no
dt→0). That is the lighter faithful test the deferral leaves on the table; see P6-2 (proposed
NEW-I.13). **RULING: the np≤4 het DRIVER end-to-end smoke is reasonably DEFERRED to Phase 7** (it
carries real driver-fixture cost — hardcoded boundary attrs, config plumbing — same basis as the
Test 5.1/5.2 deferral ruled acceptable in Rev I-6). Phase 6 acceptance #4 is NOT met by an
implemented test, but the deferral is acceptable FOR MERGE because the central-flux dispatch physics
is covered by the gate (#2) + Test 3.2 + the CFL guard, and the RK-coupling itself is covered by the
existing RK-stepper unit suites. A stepper-level local smoke (NEW-I.13) SHOULD be added before relying
on the het RK path in CI.

**#9 — Regression: parity/byte-exact + Test 3.2. — UNCHANGED.**
Parity C-6 (9/9), phaseh serial (19/19) + np=4 (52/52), central-flux (6/6), dispatch incl 3.2 (10/10),
CFL (17/17), shared np=2 (14/14). Zero new failures, zero changed values. The Phase-6 diff is
ADDITIVE (one new test file + Makefile registration) and touches no operator/flux source — byte-exact
contract intact by construction. Tree restored to exactly the Phase-6 scope after the scratch
experiment (godunov_flux_bimaterial.cpp md5 confirmed `58d91bef8bfe60b00e77c6db423eeab2`; no scratch
files left).

## Bugs Found — Phase 6 (Rev I-8)

### P6-1 — `AllContinuous` run is a DUPLICATE of the `Adjacent` run on the gate fixture (coverage gap)

**Status:** OPEN
**File:** `tests/unit/test_bimaterial_mixed_flux_homog_equivalence.cpp:23-26, 306-307`
**Severity:** MODERATE
**Description:** On the 2×2×1 fixture every element is fault-adjacent (all four hexes touch the y=1
fault plane), so the `Adjacent` element-walk (`wave_operator.inl:1764-1806`) visits every interior
non-fault face — making Adjacent's central set IDENTICAL to AllContinuous's (`:1748-1758`). The
reproduced run confirms `central faces = 2` for BOTH modes. The gate therefore runs the SAME
face-set twice and the "AllContinuous" pass does NOT exercise AllContinuous's distinguishing
all-interior-faces construction (a non-fault interior face NOT adjacent to the fault, which
AllContinuous includes but Adjacent excludes).
**Expected:** The gate should, for AllContinuous, exercise a face set that DIFFERS from Adjacent's —
i.e., contain at least one non-fault interior face not adjacent to the fault.
**Trigger:** A defect specific to AllContinuous's face-set loop (e.g., a wrong `is_interior_face`
predicate that only AllContinuous reaches, or an off-by-one in the `for(f<GetNumFaces())` walk) would
NOT be caught by this gate — it would slip through because Adjacent's element-walk happens to produce
the right answer on this fixture and the gate compares the two modes to the SAME scalar reference.
Note: the AllContinuous-vs-Adjacent DISTINCTION is covered by the 5-hex `*_dispatch` fixture
(Test 3.1), so this is a gate-local over-claim, not an absolute hole.
**Suggested fix:** Use a fixture with a non-fault interior face NOT adjacent to the fault — e.g. a
3×1×1 (or 5-hex, mirroring the dispatch fixture) block with the fault on ONE interior plane, so the
far interior face is in AllContinuous's set but NOT Adjacent's. Then assert
`AllContinuous.central_set.size() > Adjacent.central_set.size()` to PROVE the two modes are distinct
on the fixture, before running the equivalence check for each. Cheaper alternative: keep the current
mesh but add a 1-line `TEST_ASSERT` documenting that on THIS fixture the two sets coincide
(`Adjacent.size() == AllContinuous.size()`), so the duplicate is intentional + visible rather than an
accidental blind spot.
**(MOD) test:** see proposed NEW-I.12.

### P6-2 — Acceptance #4 (het RK local smoke) has ZERO implemented coverage; a faithful stepper-level test is feasible and was not added

**Status:** OPEN — DEFERRED (deferral ruled ACCEPTABLE; see target #8)
**File:** (no file — missing test) plan §Phase 6 req 3 + acceptance #4
**Severity:** MODERATE
**Description:** Phase 6 acceptance #4 ("np≤4 local het smoke advances ≥5 steps without NaN / dt→0")
is documented DEFERRED to Phase 7. The full-driver deferral is justified (real fixture cost). But a
FAITHFUL lighter test exists: `AdvanceRKCoupledLSW_Spatial` (and the RS sibling) is templated on
`MeshType` and explicitly serves serial-`Mesh` unit tests (`rk_time_stepper.hpp:367-368`), so a
2-element het fixture + `BimaterialWaveOperator` + `SetMixedFluxMode(Adjacent)` + `SetCflRkAware(true)`
+ ≥5 invocations of the coupled stepper would faithfully exercise the matrix+adjacent+RK4 path
(including the stage-local slip write that a hand-rolled RK4 cannot replicate) WITHOUT the driver's
config/boundary plumbing.
**Expected:** ≥5 coupled-RK steps on a het 2-material operator complete with finite Q and a
non-collapsing dt.
**Trigger:** A latent defect in the het RK-coupled path (e.g., the central+RK CFL factor, the
per-element het CFL walk under mixed flux, or a NaN from the stub neighbour material on the local
fixture) would not surface until Phase 7 / Frontera.
**Suggested fix:** Add the stepper-level smoke (NEW-I.13) before Phase 7, OR explicitly record in the
Phase-6 acceptance ledger that #4 is intentionally Phase-7-owned (the implementer DID document this,
which is why severity is MODERATE-deferred, not HIGH).
**(MOD) test:** see proposed NEW-I.13.

### P6-3 — `!empty()` gate assertion does not verify the central set is the RIGHT faces (only non-empty + same count)

**Status:** OPEN
**File:** `tests/unit/test_bimaterial_mixed_flux_homog_equivalence.cpp:240-244`
**Severity:** LOW
**Description:** The gate asserts `!GetCentralFluxFaceSet().empty()` (scalar) and
`size()==size()` (scalar vs bimaterial). This catches the degenerate-to-upwind failure (empty set ⇒
re-proving C-6) and a count mismatch, but does NOT assert the count equals the EXPECTED 2, nor that
the member faces are the x=1 interior faces (and not, e.g., a wrongly-classified external face that
slipped past `is_interior_face`). On this fixture the topology forecloses a wrong-face set, and the
teeth experiment (#2) confirms correctness empirically — so this is hardening, not a live bug.
**Expected:** A self-documenting assert pinning the expected central count.
**Trigger:** A future fixture edit (e.g. adding a z-layer) that silently changes the central count
would pass the current `size()==size()` check as long as both operators agree — masking a fixture
regression.
**Suggested fix:** Add `TEST_ASSERT(wave_scalar.GetCentralFluxFaceSet().size() == 2, …)` and,
optionally, verify each central face's centroid is at x=1 (mirror the y-centroid logic the mesh
builder already uses).

### P6-4 — Doc/hardening nits (stale 8.6e-16 constant; SetupFault single-face nqp assumption)

**Status:** OPEN
**File:** `tests/unit/test_bimaterial_mixed_flux_homog_equivalence.cpp` (header comment region; `:119-152`)
**Severity:** LOW
**Description:** (a) The task/plan-era narrative cites Test 1.1 == scalar central at "8.6e-16"; the
reproduced figure is 1.470e-15 this run — both sub-1e-9 LU-rounding floors, but the literal constant
is stale wherever it appears as an exact claim. (b) `SetupFault` derives `nqp_per_face` from
`int_faces[0]` only and assumes all fault faces share it (inherited verbatim from the parity test);
correct on this uniform-hex mesh but would silently under/over-size `dof_data` on a mixed-geometry
fault. Neither is a live Phase-6 bug; both are documentation/robustness.
**Suggested fix:** (a) Soften any exact "8.6e-16" citation to "sub-1e-9 / LU-rounding floor"; (b) add
a one-line comment in `SetupFault` that the single-face nqp read assumes geometry-uniform fault
faces (true for this fixture), or `MFEM_VERIFY` per-face nqp equality if a mixed fault is ever used.

## Verification of prior findings (carried into the Phase-6 tree)

| Prior ID | Status carried | Note |
|----------|----------------|------|
| BUG-3 (per-side single-valuedness) | guarded | Test 3.2 PASS (10/10 dispatch); the gate's homogeneity MASKS BUG-3 by design (necessary-not-sufficient note correct) |
| IMPL-8 (lifecycle: `per_face_central_flux_.size()==central_flux_face_set_.size()`) | holds | The gate's `SetMixedFluxMode` re-entry (Adjacent then AllContinuous in separate RunGate calls, fresh operators) builds size-2 central maps both times; assert at `bimaterial_wave_operator.inl:877` did not fire |
| P3-4 / Phase-4 central+ADER guard | holds | `ComputeMaxDt` guard intact (CFL suite 17/17); gate does not exercise dt (operator-level Mult only) |
| P5-1 shipped-TOML ADER landmine | unchanged | Out of Phase-6 scope; still OPEN |

## Proposed Unit Tests (Implementation Rev I-8)

| ID | Name | Target | Validates | Priority |
|----|------|--------|-----------|----------|
| NEW-I.12 (NEW in Rev I-8) | `HomogEquiv_AllContinuous_DistinctFromAdjacent` | gate fixture / `BuildCentralFluxFaceSet_` AllContinuous arm | a fixture with a non-fault interior face NOT fault-adjacent (3×1×1 or 5-hex) makes `AllContinuous.central_set.size() > Adjacent.central_set.size()`, then re-runs the homogeneous-equivalence check on the LARGER AllContinuous set — closes the P6-1 gate coverage gap so AllContinuous's distinct face-set construction is genuinely exercised by the gate | **MEDIUM** (closes P6-1) |
| NEW-I.13 (NEW in Rev I-8) | `HetRKCoupled_LocalSmoke_5Steps` | `AdvanceRKCoupledLSW_Spatial` on a het `BimaterialWaveOperator` (serial `Mesh`) | a 2-material 2-element fixture + `SetMixedFluxMode(Adjacent)` + `SetCflRkAware(true)` advances ≥5 coupled-RK steps with finite Q + non-collapsing dt — the FAITHFUL local smoke acceptance #4 calls for, feasible at unit scale via the templated stepper (NO driver fixture) | **MEDIUM** (closes the acceptance-#4 gap P6-2 without Phase-7 cost) |
| NEW-I.14 (NEW in Rev I-8) | `HomogEquiv_CentralSet_IsX1Faces` | gate fixture central-set membership | the scalar central set has EXACTLY 2 faces and each is at x=1 (centroid check) — pins the fixture against a silent topology regression (P6-3) | LOW |

## Priority Summary (Implementation Rev I-8)

| ID | Severity | Status | Description |
|----|----------|--------|-------------|
| P6-1 | MODERATE | OPEN | Gate's AllContinuous run is a DUPLICATE of Adjacent (all elements fault-adjacent on 2×2×1) — AllContinuous's distinct face-set construction is NOT exercised by the gate (covered elsewhere by Test 3.1) |
| P6-2 | MODERATE | OPEN — DEFERRED (acceptable) | Acceptance #4 (het RK local smoke) has zero implemented coverage; a faithful stepper-level test (NEW-I.13) is feasible without the driver but was deferred to Phase 7 |
| P6-3 | LOW | OPEN | `!empty()` + `size()==size()` asserts do not pin the expected central count / face identity |
| P6-4 | LOW | OPEN | Stale "8.6e-16" citation (reproduced 1.47e-15); `SetupFault` single-face nqp assumption (inherited) |

**No CRITICAL. No HIGH.** The gate is NOT vacuous (#1 cleared), has TEETH (#2 proven by the
2e-7 → 1-3e-7 perturbation experiment, 6 fails), and the Makefile registration is correct (#7). The
two MODERATE findings are coverage gaps (P6-1 gate-local AllContinuous duplication; P6-2 the
acceptance-#4 deferral), not correctness defects in shipped code.

### Test Coverage Summary (cumulative, after Phase 6)

| Suite | File Under Test | Tests reproduced | Key items |
|-------|-----------------|------------------|-----------|
| **homog-equivalence GATE** | `bimaterial_wave_operator.inl` central path | **18/18** | Test 6.1 (Adjacent + AllContinuous, ≤1e-9); teeth-verified |
| bimaterial central flux | `godunov_flux_bimaterial.cpp` | 6/6 | Test 1.1 homogeneous == scalar central (1.47e-15) |
| mixed-flux dispatch | `bimaterial_wave_operator.inl` | 10/10 | Test 3.2 het single-valuedness (BUG-3 guard) |
| bimaterial CFL | `wave_operator.inl` / `bimaterial_wave_operator.inl` | 17/17 | central+ADER guard |
| bimaterial parity | `bimaterial_wave_operator.inl` | 9/9 | C-6 |
| phaseh const parity | `wave_operator.inl` | 19/19 + 52/52 (np=4) | byte-exact |
| bimaterial shared (np=2) | `bimaterial_wave_operator.inl` | 14/14 | Test 3.4 (NOT in serial group) |
| **Het RK local smoke (accept #4)** | `AdvanceRKCoupledLSW_Spatial` | **0 (deferred)** | NEW-I.13 recommended; deferral acceptable |

## Verdict (Phase 6 implementation)

**PASS WITH FIXES.** The homogeneous-equivalence gate is MEANINGFUL, not vacuous: it exercises the
bi-material central flux on the 2 correct x=1 interior faces (#1), the central contribution is
non-negligible (a 2e-7 central-matrix perturbation drives the field-relative diff to 1-3e-7, failing
6/18 — the gate has genuine TEETH, #2), the 1e-9 field-scale metric is appropriate and achievable
(scalar central == bimaterial central to LU rounding 1.47e-15, #5), the necessary-not-sufficient NOTE
+ the Test 3.2 sufficient-guard cross-reference are correct (#6, acceptance #3 PASS), and the
Makefile registration is complete + correct + serial-grouped (#7). All eight regression/parity suites
are byte-exact unchanged (#9). The assumed "≥3 bugs" did NOT surface as correctness defects in
shipped code; the four findings are two MODERATE coverage gaps (P6-1 the gate-local
Adjacent≡AllContinuous duplication — AllContinuous's distinguishing construction is exercised only by
the upstream Test 3.1 fixture, not by the gate; P6-2 the acceptance-#4 het-RK-smoke deferral) and two
LOW hardening/doc nits (P6-3, P6-4).

**Ruling on the acceptance-#4 deferral (#8): ACCEPTABLE for merge, but a lighter faithful test SHOULD
be added.** The implementer's rationale is sound and verified against `rk_time_stepper.hpp` — a
hand-rolled RK4-on-Q would skip the stage-local slip write (`:429-441`) and feed the friction solve
stale slip, so a faithful smoke genuinely needs the coupled stepper. The full-DRIVER end-to-end smoke
is reasonably Phase-7-owned (real fixture cost, same basis as the Rev-I-6 5.1/5.2 deferral ruling).
HOWEVER, `AdvanceRKCoupledLSW_Spatial` is templated for serial `Mesh` unit tests, so a FAITHFUL
stepper-level local smoke on a 2-element het fixture (NEW-I.13) is feasible at unit scale with no
driver plumbing and would close the acceptance-#4 gap properly. Add NEW-I.13 (and ideally NEW-I.12 to
de-duplicate the gate's AllContinuous arm) before relying on the het RK path in CI.

---

# Final Whole-Feature Review (Rev I-10)

**Date:** 2026-06-04
**Reviewer:** chunhui-code-reviewer agent (adversarial pre-merge gate; "assume ≥3 bugs, hunt CRITICAL")
**Scope:** The ENTIRE Mixed-Flux + Heterogeneous (Bi-material) Riemann feature, Phases 1-7, as a
single cross-cutting audit. NOT a per-phase re-verify — this round re-executes all three passes
(code review, latent-bug hunt, full test re-run) from scratch on the COMPLETE working-tree diff to
catch INTEGRATION bugs across phase boundaries that the per-phase rounds (Rev I-1/I-2/I-4/I-6/I-8)
could not see. Prior history (BUG-1..26, IMPL-1..8, P3-*, P4-*, P5-*, P6-*) left intact above.

**Rev numbering:** continues the implementation-review thread (I-1 Phase 1&2, I-2 Phase 3, I-4
Phase 4, I-6 Phase 5, I-8 Phase 6). This is **I-10** (the whole-feature pre-merge gate). New IDs
here are `F-1`, `F-2`, ….

> **Review history append:** Rev I-10 (2026-06-04): whole-feature cross-cutting audit; 0 CRITICAL,
> 0 HIGH, 0 MODERATE-new; 3 LOW findings (F-1 gate edge, F-2 carried P5-1 reaffirmed, F-3 doc/test
> nits). All prior fixes independently re-verified on a fresh rebuild. **MERGE: PASS WITH FIXES**
> (the 3 fixes are LOW; none blocks merge).

## Builds / Tests Re-Run — fresh objects, counts NOT trusted from prior rounds

All builds in the worktree against the main checkout's `libmfem.a`
(`MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN MFEM_LIB_DIR=$MAIN`). The `.inl`-dependent
objects (`dynamic/wave_operator.o`, `dynamic/bimaterial_wave_operator.o`,
`dynamic/godunov_flux_bimaterial.o`) and every test binary were `rm`'d before rebuilding (unreliable
`.inl` dep tracking). No new compiler warnings on any target.

| Target | Result reproduced | Notes |
|--------|-------------------|-------|
| `seas_test_bimaterial_central_flux` (serial) | **6 / 6** | Test 1.1 homog==scalar central; Test 1.3 het anchors 3.9e-16/6.3e-16/5.2e-16 |
| `seas_test_bimaterial_mixed_flux_dispatch` (serial) | **10 / 10** | Test 3.2 het single-valuedness (BUG-3) PASS; Test 3.1 Godunov per-side contrast 3e-1 |
| `seas_test_bimaterial_mixed_flux_cfl` (serial) | **17 / 17** | Test 4.1 factor parity; Test 4.2 central+ADER abort on both operators; None+ADER no abort |
| `seas_test_bimaterial_mixed_flux_homog_equivalence` (GATE) | **18 / 18** | Test 6.1 Adjacent 9 + AllContinuous 9, ≤1e-9; central faces = 2 |
| `seas_test_wave_operator` (serial) | **25 / 25** | Test 15 `TestMatrixMixedFluxEnabled` (R-003 lifted); Test 16 (R-006) intact |
| `seas_test_bimaterial_wave_operator_parity` (serial) | **9 / 9** | C-6c ComputeMaxDt parity (homogeneous → None, byte-exact) |
| `seas_test_phaseh_wave_operator_constant_parity` (serial) | **19 / 19** | `mixed_flux=none` byte-exact |
| `seas_test_phaseh_wave_operator_constant_parity` (np=4) | **52 / 52** | `mixed_flux=none` byte-exact, parallel |
| `seas_test_bimaterial_mixed_flux_shared` (np=2) | **14 / 14** | Test 3.4a (lateral μ(x)+seam) + 3.4c (depth μ(z)+seam, no false-abort); IMPL-8 lifecycle |
| `seas_test_tpv_config_parse` (serial) | **85 / 87 (2 fail)** | the 2 fails are `3 VS-border box rules present` (:191, :239) — CONFIRMED PRE-EXISTING (see below) |
| `tpv31/test_visualize_results_tol.py` (pytest) | **9 / 9** | NEW-7.1/7.4/7.5 + mu_eff-exclusion + dip-channel + help-flag |
| `seas_spatial_dyn_driver` (heavy link, fresh `.o`) | **compiles + links clean** | no errors/warnings; not run on a real mesh (no-local-mesh rule) |

**Reproduced total: 255 assertions passed, 2 failed** (6+10+17+18+25+9+19+52+14+85+9, minus the 2
documented pre-existing VS-border fails). Every count matches the prior rounds; nothing regressed.

**The 2 config-parse failures are CONFIRMED PRE-EXISTING, not new.** Independently verified:
`grep` of the feature diff (`git diff tests/unit/test_tpv_config_parse.cpp`) shows the failing
assertion text `3 VS-border box rules present` / `n_vs == 3` is **NOT in the diff** — i.e. those two
`TEST_ASSERT`s were neither added nor changed by this feature. They exist verbatim at `git show
HEAD:...test_tpv_config_parse.cpp` lines 166/214; the feature's additive tpv31 assertions shifted
them down to 191/239. They are an unrelated pre-existing TPV205/TPV104 VS-border fixture mismatch.
The feature adds 15+ tpv31 assertions, ALL passing, and changes ZERO previously-asserted field values.

**Working tree unchanged after review:** md5 of all five core source files
(`godunov_flux_bimaterial.cpp 58d91bef…`, `bimaterial_wave_operator.inl 5df460be…`,
`wave_operator.inl e56e3678…`, `spatial_dyn_driver.cpp 5800cf9d…`, `spatial_friction.cpp f9a2c924…`)
match the pre-review baselines exactly; `git diff --stat` is byte-for-byte unchanged (729 ins / 148
del across dynamic/drivers/spatial). No source was perturbed (no perturbation experiment was needed —
the prior rounds' teeth-experiments are accepted; this round verified by inspection + the existing
discriminating tests).

## Priority adversarial cross-cutting targets — verdicts

### #1 — End-to-end single-valuedness (BUG-3), MATRIX ASSEMBLY → SIDE-SYMMETRIC STORE → DISPATCH — CORRECT on BOTH the interior 2-sided path AND the shared rank-seam path.

Traced the full chain on the complete diff:

- **Primitive** (`godunov_flux_bimaterial.cpp:307-381`): `BuildPerFaceCentralMatricesGlobal` builds
  `centralSelf = ½·T·(AxPlus_self+AxMinus_self)·Tinv` and `centralNbr = ½·T·(AxPlus_nbr+AxMinus_nbr)·
  Tinv` — **each side's OWN Jacobian**, ½ applied per-output (`*= 0.5`), NO swap, two distinct scratch
  matrices (`Asum`,`tmp`), output sized via `SetSize` (never aliased). This is the correct
  `½A_self·Q_self + ½A_nbr·Q_nbr` formula, NOT the `½A_self·(Q_self+Q_nbr)` pitfall. Test 1.3 anchors
  it independently (the neighbour term uses A_nbr, ≤6.3e-16). **No defect.**
- **Side-symmetric store** (`bimaterial_wave_operator.inl:512-560` interior arm): one pair per
  `mesh_face`, `per_face_central_flux_[mesh_face] = {½A_e1, ½A_e2}` from `ResolveFaceFluxOperands_`
  (self=Elem1, nbr=Elem2, NO swap parameter — IMPL-1). The upwind per-side swap (`compose_side`) is
  **not** applied to the central build. This is the structural BUG-3 guarantee.
- **Interior dispatch** (`:725-739`): `ApplyPerFaceFlux(c[0],c[1],Q_self,Q_nbr,F_h_e1)` computes
  `F* = ½A_e1·Q_e1 + ½A_e2·Q_e2` ONCE, then `F_h_e2[comp] = F_h_e1[comp]` for all 9 comps and
  `return`s before the Godunov per-side apply (no fallthrough, no double-write). The runtime passes the
  UNSWAPPED `(Q_self=Q_e1, Q_nbr=Q_e2)` to both sides (confirmed `wave_operator.inl:3124`), so `F*` is
  identical on both sides for ANY materials. Test 3.2 enforces `F_h_e1 == F_h_e2` ≤1e-12 on a GENUINE
  μ-contrast fixture (worst rel reproduced ~0); the companion Godunov face gives 3e-1 per-side
  difference, proving the per-face kernel selection is real and the central guard non-vacuous.
- **Shared (rank-seam) dispatch** (`:759-767`): `SharedInteriorFaceFlux_` deposits the SAME single
  `F* = ½A_local·Q_self + ½A_nbr·Q_nbr` for side 0 (local Elem1) — the peer rank deposits its own
  side. The shared store (`:562-640` sf-arm) builds `(½A_local, ½A_nbr)` from the SAME `(nor, flux_local,
  flux_nbr)` ordering as the Godunov shared build. **Test 3.4 (NEW-I.4/5/6) at np=2 verifies the shared
  store routes bit-for-bit through `per_face_central_flux_` (3.4a) and the stored halves `c[0]==c[1]`
  under the R-004 stub (single-valued by construction).** This is the SHARED-path single-valuedness
  evidence the task asks whether Test 3.4 (not just the build) verifies — it does (the bit-for-bit
  store-route check + the c[0]==c[1] entrywise anchor; the cross-rank value anchor was replaced in the
  fix round because the independently re-derived shared normal was run-to-run nondeterministic — see
  fix-doc deviation #1). **No defect on either path.**

### #2 — Conservation / sign at the central face — CORRECT, matches the scalar oracle.

`F_h_e1 == F_h_e2` (same-sign copy, no negation) is deposited by the DG assembly as
`rhs[e1] -= w·shape1·F_h_e1` (`wave_operator.inl:3130`) and `rhs[e2] += w·shape2·F_h_e2` (`:3132`)
— the standard conservative convention where e1's outward normal is `+nor`, e2's is `-nor`, encoded
in the `-=/+=` signs. A single-valued physical `F*` deposited with opposite signs is conservative.
This is **byte-identical to the scalar reference** (`wave_operator.inl:1209`, `F_h_e2[c]=F_h_e1[c]`)
— the scalar central path does NOT negate for e2 either. For the shared path the local-only `-=`
deposit (`:3729`) conserves across the seam because both ranks compute the SAME `F*` from
rank-local opposite `nor` on swapped Q (each rank's `Q_self` is the peer's `Q_nbr`); the face-normal
Jacobian is odd in `nor` (`A(-nor) = -A(nor)`), so the peer's matrices are the negation with self/nbr
swapped — the SAME conservation mechanism the production-validated bi-material UPWIND shared path
already uses (it shares the operand derivation, `ApplyPerFaceFlux`, and `-=` assembly verbatim). The
homogeneous-equivalence GATE (Test 6.1, 18/18) confirms the bimaterial central == the proven-
conservative scalar `GodunovFlux::Central` through the FULL `Mult` assembly to ≤1e-9; under
seam-continuous heterogeneity the materials match across the seam, reducing to the validated case.
**No same-sign-without-matching-assembly defect; no conservation break in the supported regime.**

### #3 — Byte-exactness, WHOLE feature — CONFIRMED on a fresh rebuild.

`mixed_flux=none` on the matrix path AND the entire scalar path are bit-unchanged:
- **Matrix `mixed_flux=none`:** `central_flux_face_set_` empty ⇒ `BuildPerFaceCentralFluxMatrices_`
  no-ops (map empty), the dispatch's `mf_on_ && count(...)>0` short-circuits doubly-false (mf_on_ is
  false AND the set is empty), and the bimaterial `ComputeMaxDt` factor is `MixedFluxCflFactor_()==1.0`
  giving `dt_e = 1.0*cfl*h_e/cp_e` (byte-identical to the deleted `cfl_mixed_flux_factor=1.0` arm).
  Reproduced: parity **9/9** (incl C-6c ComputeMaxDt), phaseh **19/19** serial + **52/52** np=4.
- **Scalar path:** `WaveOperator::ComputeMaxDt` returns `MixedFluxCflFactor_() * cfl * h_min_ /
  flux_.GetCp()` = `1.0 * cfl * h_min_ / cp` for None — IDENTICAL expression + evaluation order to the
  pre-Phase-4 `cfl_mixed_flux_factor * cfl * h_min_ / flux_.GetCp()`. No FP reordering. Verified live
  (phaseh C-3 bit-identical, parity C-6c).
- Reproduced the full required panel: parity 9/9, phaseh 19/19 + 52/52 (np4), wave_operator (Test 15)
  25/25, central 6/6, dispatch 10/10, cfl 17/17, homog-equivalence 18/18, shared 14/14 (np2),
  config-parse 85/87 (2 PRE-EXISTING VS-border, confirmed not new). **No NEW failure; no non-None byte
  change. The byte-exact regression contract holds.**

### #4 — Re-entrancy / lifecycle across phases — CORRECT; IMPL-8 assert robust at np>1 and on the empty-set case.

`SetMixedFluxMode(None)→(Adjacent)→(None)`: the base clears `central_flux_face_set_` + sets
`mf_on_=false` for None (`wave_operator.inl:1643,1655`); the override's
`BuildPerFaceCentralFluxMatrices_` clears `per_face_central_flux_` FIRST (`:506`), so no stale matrices
survive a None flip and a later `.at()` cannot hit a stale key. The IMPL-8 size-equality assert
(`:877`, `per_face_central_flux_.size() == central_flux_face_set_.size()`) holds on EVERY transition
including the degenerate empty-set (0==0) and at np>1 (verified live, 14/14): every map key is a
provable subset of the set (interior arm inserts only keys it iterates from the set; sf-arm inserts
only when `count(mesh_face_idx)>0`), so equal cardinality ⇒ exact set ⇒ no false-PASS. **No false-FIRE
on the one theoretical edge** (a central-set SHARED face with null `GetSharedFaceTransformations`): the
P3-2 fix reorders the sf-arm so `is_central` is computed before the `ftr` lookup and a null `ftr` on a
central face is a loud `MFEM_VERIFY(!is_central, …)` (`:604`), NOT a silent `continue` that would leave
the face unbuilt and trip the size assert. **No false-pass: a shared central face cannot be double-built**
— the interior arm's `ResolveFaceFluxOperands_` returns false for a shared face (Elem2No<0, the same
`e2<0` skip the production-validated upwind Pass-1 build relies on), so only the sf-arm builds it.

### #5 — Exhaustive guard claim + central+ADER guard consistency — CONFIRMED; NO 4th blocking guard; the feature CAN run end-to-end.

Independently re-swept `MFEM_VERIFY`/`MFEM_ABORT` mentioning interior_flux/mixed_flux/matrix/rk/ader/
central across `drivers/spatial_dyn_driver.cpp` + `spatial/code/spatial_friction.cpp` +
`dynamic/{bimaterial_,}wave_operator.inl`. The ONLY guards on the matrix+adjacent+RK path:
- **G1** `spatial_friction.cpp:1162-1167` — RELAXED (the mutual exclusion is gone; only the value
  validator `:1152` remains, which accepts `adjacent`). ✓
- **G2** `spatial_dyn_driver.cpp:823` — `MFEM_VERIFY(!spatial::MatrixMixedFluxUnderAder(cfg), …)`,
  reading the CLI-overridden integrator (override @665-672 BEFORE G2 @823). The predicate
  `interior_flux==Matrix && mixed_flux!="none" && time_integrator==ADER` is LOGICALLY EQUIVALENT to the
  old inline `matrix_mixed && !is_rk` (verified). Blocks ONLY matrix+mixed+ADER; allows matrix+mixed+RK,
  matrix+none+either, scalar+anything. ✓
- **G3** `bimaterial_wave_operator.inl` SetMixedFluxMode — the `MFEM_VERIFY(m==None)` abort is REPLACED
  by `Base::SetMixedFluxMode(m) + BuildPerFaceCentralFluxMatrices_() + IMPL-8 assert`. ✓
- The kept-correct guards (`spatial_friction.cpp:1707` matrix-needs-non-Constant; the base
  `wave_operator.inl:1615` precomputed-flux exclusion, false here) do NOT trip TPV31 (depth_profile_1d,
  non-Constant; use_precomputed_face_fluxes_ false). ✓
- **NO 4th blocking guard.** The driver compiles + links + the operator-construction order is correct;
  the feature reaches the time loop on matrix+adjacent+RK+seam_continuous=true.

**Central+ADER abort consistency:** the operator guard
`MFEM_VERIFY(mixed_flux_mode_ == MixedFluxMode::None || cfl_rk_aware_, …)` is VERBATIM identical
(condition) on both `WaveOperator::ComputeMaxDt` (`:5888`) and `BimaterialWaveOperator::ComputeMaxDt`
(`:803`); messages differ in prose only. The driver G2 uses the equivalent config-level predicate
(matrix path). On the matrix+mixed+RK run, G2 passes (RK), `SetCflRkAware(true)` (`:1933`) runs before
`ComputeMaxDt` (`:1934+`), so the operator guard passes too — G2 is the loud-at-config gate, the
operator guard is belt-and-suspenders. `MixedFluxCflFactor_()` stays a PURE factor (the ONLY abort is
the R-1600 unknown-mode fall-through; the central+ADER abort lives in each `ComputeMaxDt`, not in the
factor), so Test 4.1 probes the `!cfl_rk_aware_` 0.9/0.4 factors without the guard firing (17/17). ✓

### #6 — MPI shared-face correctness (BUG-5/6/10/15/16) — CORRECT at np=2; no observed deadlock.

The neighbour-material local-side stub + the `sf`-keyed shared arm + the declarative `seam_continuous`
gate are exercised by `seas_test_bimaterial_mixed_flux_shared` (np=2, 14/14): a `Mode::Coefficient`
μ(x) shared central face WITH `SetSeamContinuous(true)` builds (3.4a), and a depth-only μ(z)
`Mode::Coefficient` shared central face on a z-split seam builds with NO false-abort (3.4c). The gate
`MFEM_VERIFY(material_->mode == Constant || seam_continuous_, …)` (`:622`) aborts a `Mode::Coefficient`
central shared face when `seam_continuous_` is false (3.4b — construction-verified, not automated:
MPI_Abort terminates the job and is uncatchable without `MFEM_USE_EXCEPTIONS`; the positive 3.4a/3.4c
paths are the exact negative of the abort, so the gate's discrimination is established). **No latent
deadlock:** the only collectives in the central build are the rank-0-gated `MPI_Comm_rank` for the
print (cheap, runs on all ranks) and the byte-accounting; the central build adds NO `MPI_Allgatherv`/
`MPI_Allreduce` inside a rank-conditional branch — the base's `BuildCentralFluxFaceSet_` Allgatherv
runs unconditionally on all ranks (the documented R-001 deadlock-avoidance pattern). np=4 phaseh
(52/52) + np=2 shared (14/14) ran clean (no hang, deterministic).

### #7 — Uninitialized / NaN / resource — CLEAN.

The (1,1,1) sentinel `flux_` on the bimaterial ctor is NEVER read by the central or Godunov dispatch
(central uses precomputed `per_face_central_flux_`; Godunov uses `per_face_bimaterial_flux_`; bimat
`ComputeMaxDt` uses the per-element `per_elem_lmr_` walk, not `flux_.GetCp()`). `per_face_central_flux_`
DenseMatrices are sized via `SetSize` inside `BuildPerFaceCentralMatricesGlobal`. The acoustic guard
`MFEM_VERIFY(flux_self.GetMu() > 1e-12 && flux_nbr.GetMu() > 1e-12, …)` is present at the top of the
central builder (BUG-9). The unit-normal precondition `MFEM_VERIFY(|nor|^2-1 < 1e-10)` guards the
rotation. `std::unordered_map::at(mesh_face)` in the dispatch is unreachable on a missing key (the
IMPL-8 always-on assert at SetMixedFluxMode time forecloses a central-set face with no map entry). No
unguarded `cp/cs` division (the cp_e computation reads validated `per_elem_lmr_`). No NaN propagation
path observed.

### #8 — Phase 7 config/tooling correctness — CORRECT.

- **TOML** (`tpv31/configs/tpv31_rk_mixedflux.toml`): functional diff vs `tpv31.toml` is EXACTLY the
  two keys — `[numerics].mixed_flux "none"→"adjacent"` and `[material].seam_continuous = true` (the rest
  is comments). `interior_flux = "matrix"` preserved; `seam_continuous` correctly under `[material]`
  (line 113), `mixed_flux`/`interior_flux` under `[numerics]`. Parses (config-parse 85/87, the tpv31
  seam_continuous round-trip assertions PASS).
- **sbatch** (`jobs/tpv31_spatial/tpv31_p1_rk45_mixedflux_50m_normal.sbatch`): the `module load` +
  `export LD_LIBRARY_PATH` block is BYTE-IDENTICAL to the working `tpv31_p2_aderO3_noflux_50m_normal.
  sbatch` (`diff` clean). Pre-flight greps match the TOML: `interior_flux="matrix"`,
  `mixed_flux="adjacent"`, `[mesh].order=1` (all present). NO `--dt` in the `ibrun` command (the two
  matches are in comments). Restart chain + per-job-ID OUT present. Whitespace-tolerant banner greps
  (`time integrator:`, `mixed flux:`, `interior flux`) per BUG-13.
- **visualize_results.py `--tol-rms/--tol-peak` gate:** `t0=max(starts)`, `t1=min(ends)`,
  `np.interp` restricted to `ref_t ∈ [t0,t1]` on the SORTED MFEM grid (NOT a full-grid clamp);
  90% min-coverage guard (`covered=False ⇒ FAIL`); `mu_eff` EXCLUDED from gated channels; peak/rms
  normalized with a floor; non-zero exit on any fail / `n_stations_checked==0`. Re-ran the 9 pytest:
  **9/9 PASS** (NEW-7.1 offset grids within/out-of-band end-to-end exit codes, NEW-7.4 zero-crossing
  no-spurious-fail, NEW-7.5 truncated-run coverage fail + end-to-end exit≠0, mu_eff-never-gated,
  dip-channel-gating, help-flags).
- **README** corrected to the real 60 `.txt` / 8-column format (no longer "empty"/".dat").

## Bugs Found — Rev I-10

### F-1 — `visualize_results.py` tolerance gate silently PASSES if it checks ≥1 station but 0 channels

**Status:** OPEN
**Severity:** LOW (robustness; not reachable with the real SeisSol reference)
**File:** `tpv31/visualize_results.py`, `run_tolerance_gate` (the `if n_stations_checked == 0: return 1`
tail). 
**Description:** The gate returns `1` (FAIL) only when `n_stations_checked == 0`. If at least one
station is matched but EVERY gated channel of every station is skipped (`ref_y is None or
_is_all_nan_or_empty(ref_y)` ⇒ printed `SKIP (ref all-NaN)`), then `n_channels_checked == 0`,
`any_fail` stays false, and the function returns `0` (PASS) having validated ZERO channels — a silent
pass. This is the SAME class of defect BUG-18 fixed for `mu_eff` (an all-NaN channel comparing false),
generalized: if a future reference revision or a corrupt write makes ALL named channels (`V_strike`,
`slip_strike`, `tau_strike`, `sigma_n`) all-NaN for the matched stations, the gate greenlights the run.
**Expected:** `n_channels_checked == 0` ⇒ FAIL (exit 1) with an explicit "no gatable channels" message,
symmetric with the `n_stations_checked == 0` guard.
**Trigger:** A reference file (or a partially-written MFEM `.dat`) whose all four gated channels are
all-NaN/empty for every matched station. The committed SeisSol `.txt` reference has valid
V_strike/slip_strike/tau_strike/sigma_n, so this does not fire on the intended TPV31 gate today.
**Suggested fix:** Add after the `n_stations_checked == 0` check: `if n_channels_checked == 0: <print
"GATE FAILED — no gatable channels (all references all-NaN/empty)">; return 1`.
**(MOD test):** a pytest where the reference's four named channels are all-NaN for the only station ⇒
assert the gate exits non-zero (closes the symmetric silent-pass with NEW-7.4's mu_eff precedent).

### F-2 — Shipped `tpv102/104/205_spatial.toml` still abort on their on-disk ADER defaults (P4-1 / P5-1 re-confirmed at the whole-feature gate)

**Status:** OPEN — carried (== P4-1 / P5-1); re-confirmed unresolved at merge
**Severity:** LOW (was HIGH as P4-1, MEDIUM as P5-1; downgraded — it is a deliberate-per-plan behaviour
on configs that "real jobs CLI-override", documented in `spatial_dyn_driver.cpp:1929-1932`, and does NOT
touch the mixed-flux-hetero feature's own deliverable — TPV31 ships a NEW `--time-integrator rk45` job)
**File:** `tpv102/configs/tpv102_spatial.toml`, `tpv104/configs/tpv104_spatial.toml`,
`tpv205/configs/tpv205_spatial.toml` (scalar + `mixed_flux="adjacent"` + default ADER); abort at
`wave_operator.inl:5888` reached via `spatial_dyn_driver.cpp:1934`.
**Description:** Phase 4 made `scalar+adjacent+ADER` abort at the scalar `ComputeMaxDt`. The three
shipped spatial TOMLs still carry `mixed_flux="adjacent"` under the default ADER integrator, so a bare
`seas_spatial_dyn_driver --config tpv102_spatial.toml` (no CLI override) now hard-aborts on its own
default. G2 does NOT catch this (it gates only the matrix path); the scalar `ComputeMaxDt` guard fires
later. This is the P4-1 → P5-1 carry-over, STILL not resolved in the final tree.
**Expected:** Either migrate the three TOMLs to `time_integrator="rk45"` (matching how real jobs run),
OR add a config-named driver diagnostic so the abort points at the TOML + the remediation at config
load (not an opaque `ComputeMaxDt` abort ~1900 lines in), OR document in each TOML header that they are
"templates that MUST be CLI-overridden". Confirm with the maintainer whether scalar-adjacent-ADER was
ever a validated finite-dt path (the open P4-1 question).
**Trigger:** `seas_spatial_dyn_driver --config <tpv102|tpv104|tpv205>_spatial.toml` with no
`--time-integrator`/`--mixed-flux` flag.
**Suggested fix:** Lowest-risk: add `time_integrator = "rk45"` to the three TOMLs (or a `[meta]` note +
a driver hint string). **(MOD test):** NEW-I.11 (`Driver_ShippedSpatialToml_ADERDefault_Aborts`) — a
default-integrator load of a `scalar+adjacent` config reaches the intended abort with a config-named
message.

### F-3 — Whole-feature doc/test-hygiene nits (consolidated LOW)

**Status:** OPEN
**Severity:** LOW (doc/test hygiene; no runtime consequence)
**Items (all carried, re-confirmed present in the final tree):**
- **(P5-2)** Both `ComputeMaxDt` central+ADER messages advertise `"allcontinuous"` (no underscore);
  the parser accepts only `"all_continuous"`. A user copy-pasting the suggested token hits the value
  validator. `dynamic/bimaterial_wave_operator.inl` + `dynamic/wave_operator.inl` ComputeMaxDt
  messages. Fix: `s/allcontinuous/all_continuous/`.
- **(P5-3)** `tests/unit/test_tpv_config_parse.cpp` Test 5.3a writes the TESTGEN temp next to the
  source config; if a future regression re-adds the G1 abort the `std::remove` never runs (MFEM aborts,
  not throws) and litters the source tree. Fix: write under `$TMPDIR`.
- **(P4-3)** Dead `return 1.0;` after the `[[noreturn]]` `MFEM_ABORT` in `MixedFluxCflFactor_`
  (`wave_operator.inl:5874`) — harmless belt-and-suspenders; the comment now correctly says "unreachable …
  Kept for clarity" (the misleading `-Wreturn-type` claim was already corrected).
- **(P6-1)** The homog-equivalence gate's AllContinuous run is a DUPLICATE of Adjacent on the 2×2×1
  fixture (all elements fault-adjacent), so the gate over-claims AllContinuous coverage; AllContinuous's
  distinguishing all-interior-faces construction is exercised only by the upstream Test 3.1 fixture.
  Add NEW-I.12 (a fixture with a non-fault interior face NOT fault-adjacent) before relying on the gate
  for AllContinuous.
**Suggested fix:** batch the four nits; none blocks merge.

## Proposed Unit Tests (Implementation Rev I-10)

| ID | Name | Target | Validates | Priority |
|----|------|--------|-----------|----------|
| NEW-I.15 (NEW in Rev I-10) | `TolGate_AllChannelsNaN_Fails` (pytest) | `visualize_results.run_tolerance_gate` | a reference whose four named channels are all-NaN for the only matched station ⇒ exit non-zero (closes F-1 symmetric silent-pass) | **MEDIUM** (closes F-1) |
| NEW-I.11 (carried from Rev I-6) | `Driver_ShippedSpatialToml_ADERDefault_Aborts` | shipped `*_spatial.toml` default-ADER load | a `scalar+adjacent`+default-ADER config reaches the intended abort with a config-named message (drives the F-2 fix) | MEDIUM (closes F-2) |
| NEW-I.16 (NEW in Rev I-10) | `SharedCentral_CrossRankSingleValued` (MPI np=2, het non-stub) | the shared-seam single-valued F* under a TRUE cross-rank material (once the MPI_Allgatherv exchange lands) | the deferred cross-rank exchange's single-valuedness; not testable today (R-004 stub makes the seam materials equal by construction) — a regression guard for the deferred follow-up | LOW (deferred-feature guard) |
| NEW-I.12 (carried from Rev I-8) | `HomogEquiv_AllContinuous_DistinctFromAdjacent` | gate fixture | exercise an AllContinuous central set LARGER than Adjacent's (closes P6-1 / F-3) | MEDIUM |

## Priority Summary (Implementation Rev I-10)

| ID | Severity | Status | Description |
|----|----------|--------|-------------|
| F-1 | LOW | OPEN | `visualize_results` tolerance gate silently PASSES if ≥1 station but 0 gatable channels (all-NaN ref); add `n_channels_checked==0 ⇒ exit 1` |
| F-2 | LOW (was P4-1 HIGH / P5-1 MEDIUM) | OPEN — carried | shipped `tpv102/104/205_spatial.toml` (scalar+adjacent, default ADER) abort on their on-disk defaults at `ComputeMaxDt`; migrate to rk45 or add a config-named diagnostic |
| F-3 | LOW | OPEN | consolidated doc/test-hygiene nits: P5-2 `allcontinuous` message typo, P5-3 TESTGEN temp litter, P4-3 dead `return 1.0`, P6-1 gate AllContinuous duplication |

**No CRITICAL. No HIGH. No new MODERATE.** The whole-feature cross-cutting audit found NO integration
bug across the phase boundaries: the single-valued central flux is structurally guaranteed and verified
on BOTH the interior and the rank-seam paths (#1), conservation matches the scalar oracle and the
validated upwind shared mechanism (#2), `mixed_flux=none` + the scalar path are byte-exact on a fresh
rebuild (#3), the lifecycle/IMPL-8 assert is robust across re-entrant mode flips at np>1 (#4), there is
NO 4th blocking guard and the three central+ADER guards are consistent (#5), the MPI shared path is
correct with no deadlock (#6), and the sentinel/NaN/resource surface is clean (#7). The three findings
are all LOW (one new gate-robustness edge F-1, one re-confirmed carried config landmine F-2 that does
NOT touch this feature's own TPV31 deliverable, one consolidated doc/test-nit bundle F-3).

### Test Coverage Summary (cumulative, after the whole-feature gate)

| Suite | File Under Test | Tests reproduced | Key items |
|-------|-----------------|------------------|-----------|
| bimaterial central flux | `godunov_flux_bimaterial.cpp` | 6/6 | Test 1.1 homog==scalar central; Test 1.3 het anchors (BUG-3 primitive) |
| mixed-flux dispatch | `bimaterial_wave_operator.inl` | 10/10 | **Test 3.2 het single-valuedness (BUG-3 load-bearing)** |
| bimaterial CFL | `wave_operator.inl` / `bimaterial_wave_operator.inl` | 17/17 | Test 4.1 factor parity; Test 4.2 central+ADER abort both operators |
| homog-equivalence GATE | `bimaterial_wave_operator.inl` central path | 18/18 | Test 6.1 Adjacent+AllContinuous ≤1e-9 (conservation through full Mult) |
| wave operator | `wave_operator.inl` | 25/25 | Test 15 R-003-lifted; Test 16 R-006 intact |
| bimaterial parity | `bimaterial_wave_operator.inl` | 9/9 | C-6c ComputeMaxDt byte-exact (None) |
| phaseh const parity | `wave_operator.inl` | 19/19 + 52/52 (np=4) | `mixed_flux=none` byte-exact |
| bimaterial shared (np=2) | `bimaterial_wave_operator.inl` | 14/14 | Test 3.4a/3.4c (shared single-valued + seam gate + no-false-abort) |
| config-parse | `spatial_friction.cpp` | 85/87 (2 PRE-EXISTING) | G1 relax + seam_continuous round-trip + tpv31 regression |
| visualize tol gate (pytest) | `tpv31/visualize_results.py` | 9/9 | NEW-7.1/7.4/7.5 + mu_eff-exclusion + dip-channel + help |
| spatial driver | `spatial_dyn_driver.cpp` | compiles + links clean | G1/G2/seam wiring; no 4th guard |

## Verdict (Whole-Feature, Rev I-10)

**MERGE: PASS WITH FIXES.** The Mixed-Flux + Heterogeneous (Bi-material) Riemann feature (Phases 1-7)
is CORRECT on every load-bearing cross-cutting invariant: (1) the bi-material central flux is genuinely
single-valued AND conservative under heterogeneity on BOTH the interior 2-sided path and the shared
rank-seam path (the side-symmetric store + same-sign `F_h_e2=F_h_e1` copy + odd-in-`nor` shared
mechanism, verified by Test 3.2 het single-valuedness, the homog-equivalence GATE through full Mult,
and the np=2 shared store-route + c[0]==c[1] anchors); (2) the deposit sign matches the scalar oracle
and the validated upwind shared assembly, so no conservation break; (3) `mixed_flux=none` on the matrix
path and the entire scalar path are byte-exact on a fresh rebuild (parity 9/9, phaseh 19/19+52/52,
central/dispatch/cfl/wave_operator/homog all green; the only 2 config-parse fails are CONFIRMED
pre-existing VS-border, not new); (4) the lifecycle/IMPL-8 assert is robust across re-entrant
None↔Adjacent↔AllContinuous flips at np>1; (5) there is NO 4th blocking guard — G1 relaxed, G2 blocks
only matrix+mixed+ADER (predicate equivalent to the old guard), G3 lifted — and the central+ADER abort
is identical-condition on both `ComputeMaxDt`s + consistent with G2, while `MixedFluxCflFactor_` stays
a pure factor; (6) the MPI shared path builds correctly, the seam_continuous gate discriminates, and no
deadlock exists; (7) the sentinel/NaN/resource surface is clean. The assumed "≥3 CRITICAL bugs" did NOT
materialize — the three findings are all LOW: a gate-robustness edge (F-1, a symmetric silent-pass when
0 channels are gatable), a re-confirmed carried config landmine on the UNRELATED scalar `*_spatial.toml`
defaults (F-2, the P4-1/P5-1 carry-over, which does not touch this feature's TPV31 deliverable), and a
consolidated doc/test-hygiene bundle (F-3). None blocks merge. Land the three LOW fixes (or file them as
fast-follows) and the Frontera TPV31 RK45+adjacent dispatch gate (Phase 7, user-approved) is the final
end-to-end confirmation, per plan.
