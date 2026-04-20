# TPV102 Debug v9.0.0 Plan

> **Last updated 2026-04-20.** Phase 1 executed.  Bug localized to new
> hypothesis **H-V9-J**.  Awaiting user decision on fix direction.
> Full numeric log: `tpv102_debug_v9.0.0_fix.md`.

---

## TODO — remaining actions to close R-1001

### Blocking path (must all pass)

- [ ] **YOU:** Read SeisSol's dynamic-rupture DG assembly at
      `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/` (start
      with `FrictionLaws/` and `Kernels/`) — identify the flux formula
      it uses at fault faces.  SEAS-MFEM's `fault_face_flux.cpp` and
      `godunov_flux.cpp` already cite SeisSol's
      `FrictionSolverCommon.h` and `ElasticSetup.h`, so the port surface
      is narrow. (For Option A; §2.)
- [ ] **YOU:** Decide Option A (reformulate DG flux) vs Option B
      (friction-traction bulk-force term).  See §2 for trade-offs.
- [ ] **YOU:** Authorize the debugger to draft a patch for the chosen
      option.  Per CLAUDE.md, no source edit happens without this.
- [ ] **Debugger:** Commit the failing test `test_fault_face_flux_frame_and_flux`
      **first** (regression gate lock), then the fix as a separate commit.
- [ ] **Debugger:** Apply the chosen patch
      - Option A → edit `wave_operator.inl:807-829` and `:1174-1186`
      - Option B → new bulk-force integrator + call site
- [ ] **Debugger:** Run `make test-fault-face-flux-frame-and-flux`;
      expect `|F_h[VX_global]| ∈ (10, 1e4) m²/s²` on the post-breakaway
      fixture (§3).
- [ ] **Debugger:** Run full `make test` — all other Phase 1 unit tests
      (§3.1, §3.1c, §3.1d) and existing regression tests must still PASS.
- [ ] **YOU:** Approve Frontera FP-5 sbatch submission.
- [ ] **Debugger:** Generate `tpv102_1000m_p1_1.2s_100rank_diag.sbatch`
      from the working-module template (Appendix B.1).
- [ ] **Frontera:** Submit FP-5; expect `max ‖Q‖∞ > 1e-5` (baseline was
      1.66e-8).  Fix confirmed.
- [ ] **Debugger:** Generate `tpv102_200m_p1_2.0s_400rank_postfix.sbatch`
      (Phase 5, Appendix B.4).
- [ ] **Frontera:** Submit Phase 5; expect Pass bar 3 (SCEC TPV102
      receiver time-histories match benchmark dataset visually).  Closes R-1001.

### Documentation

- [ ] **Debugger:** Create `tpv102_debug_v9.0.0_check.md` — R-items list
      updated from v8.0.0 check-doc with H-V9-J finding.
- [ ] **Debugger:** Append each executed step to `tpv102_debug_v9.0.0_fix.md`
      as it happens (fix log, one entry per substantive change).
- [ ] **Debugger:** When R-1001 closes, update this plan to mark Phase 4 & 5
      DONE, move H-V9-J to "Closed" list in Appendix E.1.

### Deferred — only if Option A/B does not close R-1001

These tracks are frozen pre-emptively because H-V9-J explains every
observed symptom.  Re-open only if FP-5 returns `‖Q‖∞ < 1e-5` after
the fix.

- [ ] Phase 2 DIAG instrumentation (C-1…C-4 code, Appendix D).
- [ ] Phase 3 Frontera DIAG at 200 m / 100 rank (Appendix B.3).
- [ ] H-V9-D — MPI-Allreduced face-visit counter test
      `test_fault_face_visit_count.cpp` (Appendix B.2).
- [ ] H-V9-F — IP penalty scaling runs (x0, x1, x10).
- [ ] H-V9-G — mesh refinement to h = 100 m OR order bump to P2.
- [ ] H-V9-H — SeisSol convention audit (blocked on Round 1 Slot B).

---

## Dashboard — read me first

**Where we are.** The R-1001 production bug (V_max pin at 7.69871 m/s,
bulk ‖Q‖∞ stuck at 1.66e-8, exactly-zero dip channels) has a concrete
root cause:

> **H-V9-J — the Godunov flux `flux_.Interior(nor, Q_imp_plus, Q_imp_minus)`
> used at fault faces cancels identically in the tangential sub-block
> when the imposed states come from `FaultFaceFlux::Evaluate` under
> pure strike-slip geometry.**

The cancellation is algebraic, not rounding (§1).  No tangential flux
reaches bulk rhs ⇒ no radiation ⇒ the observed 10¹⁵ bulk amplitude
gap.

**What you decide today.**  Pick Option A or B in §2.  No source code
has been edited.

**What's frozen.**  Phase 0 Frontera baseline, Phase 2 DIAG
instrumentation, Phase 3 Frontera DIAG, and Phase 5 final verification
are all on hold until a fix lands and the §3.1b regression test turns
green.  Full specs preserved in the appendix.

### Phase progress

| Phase | Scope | Status | Detail |
|---|---|---|---|
| **1** | Laptop unit tests (§3.1, §3.1b, §3.1c, §3.1d) | **DONE** — A1-FAIL on §3.1b, H-V9-J identified | Appendix A, `_fix.md` |
| **4** | Apply fix + regression test | **OPEN, blocked on your choice in §2** | §2, §3 |
| 0 | Frontera 1000 m baseline | Frozen | Appendix B.1 |
| 2 | DIAG instrumentation (C-1…C-4) | Frozen — mostly irrelevant under H-V9-J | Appendix B.2, D |
| 3 | Frontera DIAG at 200 m / 100 rank | Frozen | Appendix B.3 |
| 5 | Frontera re-verification (job 7666327 equiv) | Frozen | Appendix B.4 |

### Hypothesis status at a glance

| H-V9- | Status | Evidence | Test |
|---|---|---|---|
| **J** | **RANK 1, ACTIVE** | §3.1b FAIL — F_h[VX_global] = 0 exactly | `test_fault_face_flux_frame_and_flux` |
| A | Falsified | §3.1 12/12 PASS — T_can rotation correct | `test_fault_face_flux_frame` |
| B | Falsified | §3.1 12/12 PASS — Eq. 10 strike/dip decoupled | `test_fault_face_flux_frame` |
| C | Falsified | §3.1d 4/4 PASS — interior & shared branches equivalent | `test_fault_flux_interior_vs_shared_branch_equivalence` |
| I | Falsified | §3.1c 5/5 PASS at ULP — Godunov identity holds | `test_godunov_identity_normal_reversal` |
| D | Open (Phase 2) | Not yet run; likely irrelevant under H-V9-J | — |
| E | Deprioritized | REVIEW R-1003: mass-inverse audit OK | — |
| F | Open (Phase 2) | Not yet run | — |
| G | Open (Phase 3) | Not yet run | — |
| H | Blocked | SeisSol cross-verification unavailable (Round 1 Slot B) | — |

---

## 1. The finding — H-V9-J in plain language

**Setup.** At the TPV102 hypocentre at post-breakaway steady slip
(V = O(m/s), tau ≈ 81 MPa), the DG solver invokes this chain at each
fault-face quadrature point (see `wave_operator.inl:807-829` for the
interior-fault branch, `:1174-1186` for the shared-fault branch):

```
1. Evaluate     (fault-local)       : Q_plus, Q_minus           → Q_imp_plus_local, Q_imp_minus_local
2. T_can        (local → global)    : Q_imp_*_local             → Q_imp_*_g
3. flux.Interior(can_n, Q_imp_plus_g, Q_imp_minus_g)            → F_h  (global DG flux)
4. rhs[Elem1]  -= w · shape · F_h   (DG assembly, interior face)
```

**Unit-test result (§3.1b, `_fix.md §3`):** Steps 1 and 2 give the
expected values (tau_corr ≈ −1.25 MPa; Q_imp_plus_g[VX] = −V/2).
**Step 3 returns F_h = 0 exactly** — not numerical noise but zero to
ULP.

**Why — the 3-line derivation.**  The imposed states are constructed
by Eqs. 11–12 in `fault_face_flux.cpp:117-156`.  For pure strike-slip
with bulk Q = 0 they satisfy three relationships simultaneously:

- `SXY_L = SXY_R = −eta_s·V` (Eqs. 11d / 12d — continuous corrected traction)
- `VY_L − VY_R = −V`          (Eqs. 11a / 12a — velocity jump = slip rate)
- `tau_corr = −eta_s·V`       (Eq. 10)

Plug into the S-wave split-flux formula (derivable from
`godunov_flux.cpp` Ax+/Ax-):

```
F[VY_int] = −(1/(2ρ))·(SXY_L + SXY_R) + (cs/2)·(VY_L − VY_R)
          =  −(1/(2ρ))·(−Zs·V)      + (cs/2)·(−V)
          =  cs·V/2                 − cs·V/2
          =  0                (identity, holds ∀ V, cs, ρ, Zs)
```

Likewise `F[SXY_int] = 0`.  No bulk forcing is injected through the
face.  Per the DG rhs update `rhs[Elem1] −= w · shape · F_h`, nothing
enters the bulk; bulk Q stays at noise floor; V_max never sees the
radiation-damping feedback; V_dip / slip_dip / τ_dip stay at machine
zero on the plot axes.

**Why it matches production exactly.**  Every R-1001 symptom
(§Appendix E.2) flows from F_h = 0:

| Production symptom | H-V9-J explanation |
|---|---|
| V_max pins at 7.69871 m/s, bit-identical next steps | Fault solves friction ignoring bulk feedback; reaches tau_0-only steady V; never revised |
| bulk ‖Q‖∞ = 1.66e-8 (10¹⁵ below physical) | No tangential flux injected; only numerical noise accumulates |
| V_dip, slip_dip, τ_dip exactly 0.000 on DRDG3D plots | Same cancellation in the dip-tangent sub-block; no dip-channel radiation either |
| R-101 verifier clean at 400 ranks | F_h = 0 is symmetric across ranks — no cross-rank inconsistency to detect |
| Nucleation rise phase matches DRDG3D | Rise driven by tau_0 ramp alone; bulk coupling irrelevant until post-peak |

**Falsification of the original hypothesis ladder.**  H-V9-A, -B, -C, -I
are all falsified by Phase 1.  The chain from Evaluate through T_can
through the Godunov identity is individually correct — the cancellation
arises from the algebraic structure of the **combination** of Eqs. 10,
11, 12 with the split-flux formula, not from any local bug.

---

## 2. Decision required — Option A vs Option B

Both paths close H-V9-J.  Both require user approval before any
source edit (per
`/Users/chunhuizhao/projects/seas-mfem/CLAUDE.md` "Proposing Fixes").
§3.1b is the regression gate for whichever option is chosen.

### Option A — reformulate the DG flux at the fault face

**What changes.** `wave_operator.inl:807-829` and `:1174-1186`: replace
```cpp
flux_.Interior(nor, Q_imp_plus_g, Q_imp_minus_g, F_h);   // both sides imposed
```
with a form that includes the bulk self-state explicitly, e.g.
```cpp
flux_.Interior(nor, Q_self, Q_imp_nbr_g, F_h);           // bulk on self, imposed on neighbour
```
or the full DG weak-form residual `∫ shape·(F̂(Q_imp_L, Q_imp_R) − F(Q_self)·n) dA`.
Reference: de la Puente et al. 2009; SeisSol `FrictionSolverCommon.h`,
`DynamicRuptureKernel.h`.  `fault_face_flux.cpp` (Eqs. 7–12) is
untouched.

**Pros.**  Changes are surgical (two sites in `wave_operator.inl`);
no new bulk integrator; aligns with mainstream DR-DG formulations.

**Cons.**  Depends on picking the right bulk-state-vs-imposed-state
routing; subtle sign/convention errors are easy.  Requires reading
SeisSol's exact formula first.

**What you do to unblock:**
1. Read the dynamic-rupture DG assembly in
   `/Users/chunhuizhao/projects/SeisSol/src/DynamicRupture/` (start
   with `FrictionLaws/` for imposed-state construction and `Kernels/`
   for how the bulk rhs consumes the imposed states — SEAS-MFEM gets
   the first right and the second wrong).
2. Identify which flux formula SeisSol uses at fault faces.
3. Authorise the debugger to draft a patch at the two call sites.

### Option B — add friction-traction as a direct bulk-force term

**What changes.**  Leave `flux_.Interior` calls alone.  Add a
supplementary bulk rhs term that applies `tau_corr` as a surface
force on bulk DOFs adjacent to fault faces (Nitsche-style or explicit
surface integration of corrected traction).

**Pros.**  Conceptually clean — friction physics enters as an explicit
surface force, orthogonal to the interior-flux formulation.

**Cons.**  More intrusive: requires a new integrator / assembly path;
larger diff; more risk of regressing unit tests outside §3.1b.

**What you do to unblock:**
1. Authorise the debugger to draft the bulk-force integrator addition.
2. Expect a wider review surface.

### Recommendation

**Option A.**  It is the smaller change, and the Phase 1 evidence points
at exactly two sites in `wave_operator.inl`.  SeisSol has been running
this dynamic-rupture DG formulation for years — its flux form is
battle-tested, and SEAS-MFEM's `Evaluate` / `GodunovFlux` already mirror
SeisSol's `FrictionSolverCommon.h` and `ElasticSetup.h` (explicit
in-source references).  The defect is specifically in how SEAS-MFEM
*consumes* the imposed states, not in how it *computes* them, so the
port from SeisSol should be local.  Option B remains the fallback if
Option A has a convention mismatch we can't reconcile.

---

## 3. Regression gate — `test_fault_face_flux_frame_and_flux`

After any fix, `make test-fault-face-flux-frame-and-flux` must pass
the §3.1b assertion:

```
|F_h[VX_global]| ∈ (10, 1e4) m²/s²   (analytical target ~220 at post-breakaway steady slip)
|F_h[SXY_global]| ≤ 1 Pa·m/s         (stress-flux continuity)
```

Fixture values (post-breakaway, see `_fix.md §1`):
- tau2_0 = 81 MPa, V_operating = 1 m/s (psi derived from friction balance)
- can_n = (0, −1, 0), can_t1 = (0, 0, −1) (dip), can_t2 = (1, 0, 0) (strike)
- Q_plus_local = Q_minus_local = 0 (quiescent bulk)

All other Phase 1 tests (§3.1, §3.1c, §3.1d) must continue to PASS
after the fix — failure of any indicates a regression introduced by
the fix, not a legitimate behaviour change.

After §3.1b turns green, FP-5 (Frontera 1000 m 100-rank) is the next
step: expected `max ‖Q‖∞` > 1e-5 where the baseline was 1e-8.

---

# APPENDICES (reference material — read only as needed)

## Appendix A. Phase 1 — full test details

**Status: all four tests written, compiled, run 2026-04-20.**

### A.1 Executed outcome — Q/I/T (Question / Implement / Test) tracker

| Sub-§ | Implementation file | Build target | Result |
|---|---|---|---|
| §3.1  | `tests/unit/test_fault_face_flux_frame.cpp` | `seas_test_fault_face_flux_frame` | 12/12 PASS |
| §3.1b | `tests/unit/test_fault_face_flux_frame_and_flux.cpp` | `seas_test_fault_face_flux_frame_and_flux` | **7/9 — A1-FAIL on F_h[VX]=0** |
| §3.1c | `tests/unit/test_godunov_identity_normal_reversal.cpp` | `seas_test_godunov_identity_normal_reversal` | 5/5 PASS at ULP |
| §3.1d | `tests/unit/test_fault_flux_interior_vs_shared_branch_equivalence.cpp` | `seas_test_fault_flux_interior_vs_shared_branch_equivalence` | 4/4 PASS |

Makefile: new `seas_test_*` binaries + `test-*` phonies inserted after
the existing `seas_test_fault_face_flux` target.

**Fixture correction applied during execution.**  The original §3.1 /
§3.1b "equilibrium" input (V_ini = 1e-12 m/s) is quasi-static — F_h ≈ 0
by physics, not bug.  Replaced with a post-breakaway operating-point
fixture (V = 1 m/s, tau2_0 = 81 MPa, psi derived analytically from
friction balance).  With this fixture, §3.1 strike-slip decoupling
invariants hold and §3.1b exposes the F_h = 0 structural cancellation.

### A.2 §3.1 — Evaluate in isolation (H-V9-A, H-V9-B)

Drive `FaultFaceFlux::Evaluate` in isolation:

- Equilibrium input Q_plus = Q_minus = 0: dip channel stays exactly 0
  (BP5 pure strike-slip invariant).  Strike channel gets V ≈ 0.27 m/s,
  tau2_corr ≈ −1.25 MPa.
- Perturb Q_plus.VY += δ (1 mm/s): tau1_corr shifts linearly (−4553 Pa
  vs analytical −4624 Pa, within 2 %); tau2 stays at machine noise.
- Perturb Q_plus.VZ += δ: mirror on the tangent-2 channel; tau1 stays at
  noise.

**Smoking gun:** H-V9-A fires if either perturbation leaks across
channels; H-V9-B fires if Q_imp{plus,minus} does not satisfy Eq. 10 to
machine ε.  Neither fired → both falsified.

### A.3 §3.1b — end-to-end chain (FAST-PATH, the one that fails)

Drives the full chain `Evaluate → T_can → flux_.Interior` used by the
wave operator.  Step-wise recorded values (from `_fix.md §3`):

- **Step 1 (Evaluate):** tau2_corr = −1.25e6 Pa (matches Eq. 10);
  Q_imp_plus_local[VZ] = −0.135 m/s; Q_imp_plus_local[SXZ] = −1.25e6 Pa.
- **Step 2 (T_can):** Q_imp_plus_g[VX] = −0.135 m/s (= −V/2);
  Q_imp_plus_g[SXY] = +1.25e6 Pa; all other components zero.
- **Step 3 (`flux.Interior`):** **F_h[VX_global] = 0 exactly**.
  `F_h[SXY_global]` = −2.82e−7 Pa·m/s (noise); all other components 0.

**Analytical cross-check:**  `(eta_s / rho) · V2 = 469 m²/s²` — the
scale the plan used to motivate F_h ≈ 220 m²/s².  The measured value is
10¹¹ below that.  Manual bisection (§1 derivation) isolates Step 3's
`ApplySplitFlux` as the cancellation site.

**Canonical frame used in the test fixture** (BP5 convention for TPV102):
```
can_n  = (0, −1, 0)         // Elem1-outward normal to y=0 fault
can_t1 = (0,  0, −1)        // dip, downward
can_t2 = (1,  0,  0)        // along strike
```

### A.4 §3.1c — Godunov identity under normal reversal (H-V9-I)

Test `F(+n, L, R) + F(−n, R, L) = 0` at five normals (y-axis, x-axis,
z-axis [up-fallback branch], oblique (1,1,1)/√3, tilted (0.6, 0.8, 0))
with random physically-sane Q_A, Q_B.

**Worst relative drift across all 5 cases: 4.59e−16 (single ULP).
H-V9-I falsified.**

### A.5 §3.1d — interior vs shared branch equivalence (H-V9-C)

Computes Elem1's rhs contribution two ways for the same canonical
preamble (canonical frame + Evaluate + T_can):

- Interior-fault branch logic (rank-local nor, routed Q_self_imp /
  Q_nbr_imp, no `accum_sign`).
- Shared-fault branch logic (canonical nor, unrouted Q_imp_plus_g /
  Q_imp_minus_g, explicit `accum_sign = elem1_on_plus ? +1 : −1`).

Cases: elem1_on_plus ∈ {true, false} × bulk ∈ {quiescent, perturbed}.
**All four PASS** (ULP agreement when F_h ≠ 0; noise-floor agreement
when F_h = 0).  **H-V9-C falsified.**

The algebra: matching the two branches on Elem1 reduces to the Godunov
identity (already tested in A.4), with trivial match when
elem1_on_plus = true.

### A.6 §3A — Fast-path FP-1…FP-5 outcomes

| Step | Action | Outcome |
|---|---|---|
| FP-1 | Write & run §3.1b end-to-end test | **DONE.** `F_h[VX_global] = 0` exactly (not noise) |
| FP-2 | Compare to analytical 2.2e2 m²/s² | **Case (c) fires** (exact zero) → manual bisection → **H-V9-J identified** |
| FP-3 | Write & run §3.1c Godunov identity | **PASS** at ULP on 5 normals → H-V9-I eliminated |
| FP-4 | Write & run §3.1d branch equivalence | **PASS** to ULP → H-V9-C eliminated |
| FP-5 | Frontera 1000 m confirmation | **BLOCKED** on H-V9-J fix (user approval) |

---

## Appendix B. Frozen phase plans

All four phases below are on hold.  They will resume only after
the §3.1b regression gate turns green on an H-V9-J fix.  Under
H-V9-J, most of the Phase 2 DIAG bisection is **moot** — if F_h is
already 0 upstream, the downstream C-3 / C-4 checkpoints inherit
zero and tell us nothing new.  Preserved here in case the H-V9-J
fix does not close R-1001 and we need to re-open these tracks.

### B.1 Phase 0 — Frontera reproducer with DIAG at two resolutions (frozen)

**Goal (frozen):** a pair of Frontera sbatch jobs to (a) confirm the
200 m / 400-rank pin still reproduces and (b) check whether the
radiation-failure mode is resolution-independent (1000 m / 100-rank).

**sbatch to generate (when unfrozen):**

| sbatch file | Mesh | Ranks | Nodes | tfinal | Purpose |
|---|---|---|---|---|---|
| `tpv102_200m_p1_1.4s_100rank_diag.sbatch`   | 200 m  | 100 | 2 | 1.4 s | Phase 3 confirmation |
| `tpv102_1000m_p1_1.2s_100rank_diag.sbatch`  | 1000 m | 100 | 2 | 1.2 s | Resolution control |
| `tpv102_200m_p1_2.0s_400rank_postfix.sbatch`| 200 m  | 400 | 8 | 2.0 s | Phase 5 final verification |

The 200 m / 400-rank / 2.0 s reference remains job 7666327; no need to
re-run that exact config unless a fresh baseline is wanted.

**Build flags (all):** `SEAS_DIAG_FAULT_FLUX=1`, `--debug-qnorm`, `--paraview-dt 0.1`.

**Acceptance (when resumed):**
- **A0-PASS** (resolution-independent fail): 1000 m breaks away (V_max > 0.1 m/s by t ≈ 1.0 s) AND bulk ‖Q‖∞ at t = 1.2 s ≤ 1e−6.  Conclusion: bug in fault-to-bulk coupling.
- **A0-PIVOT**: 1000 m breaks away AND bulk ‖Q‖∞ at t = 1.2 s ≥ 1e−3.  Conclusion: bug is 200 m specific (H-V9-G) or topology-specific (H-V9-D).
- **A0-FAIL**: 1000 m never breaks away (V_max < 0.1 m/s at t = 1.2 s).  Rate/state init bug at 1000 m; file a new ticket.

**sbatch skeleton.** Must inherit the working Frontera module list per
`feedback_sbatch_modules.md`; do not hand-prune `module load` lines.
```bash
#!/bin/bash
#SBATCH -J tpv102_phase0_1000m
#SBATCH -N 2
#SBATCH -n 100
#SBATCH -t 0:30:00
#SBATCH -o tpv102_phase0_1000m_%j.out
# ... module load + LD_LIBRARY_PATH copied from working sbatch ...
make seas_tpv102_driver SEAS_EXTRA_CPPFLAGS='-DSEAS_DIAG_FAULT_FLUX=1'
ibrun ./seas_tpv102_driver \
   --mesh tpv102/mesh/tpv102_1000m.msh --order 1 --cfl 0.5 \
   --tfinal 1.2 --output-dir tpv102/results_phase0_1000m_$SLURM_JOB_ID \
   --debug-qnorm
```

### B.2 Phase 2 — DIAG instrumentation (frozen)

**Hypothesis targets (original):** H-V9-D, H-V9-F, and any of -A/-B/-C/-I
surviving Phase 1.  Under H-V9-J these are all irrelevant unless a
fix reveals residual issues.

**Build flags:** `SEAS_DIAG_FAULT_FLUX`, `SEAS_DIAG_GHOST_EXCHANGE`,
new `SEAS_DIAG_BULK_MONITOR`.  C-1…C-4 checkpoint code in Appendix D.

**Watch points (with expected values):**

| ID | Location | Purpose | Expected at t ≈ 1.3 s | §0.5 cross-ref |
|---|---|---|---|---|
| W-0  | Hypo fault QP | Where radiation is born | F_h[VX_global] ≈ +2.2e2 m²/s²; F_h stress ≈ 0 (≤ 1 Pa·m/s) | C-2 |
| W-0a | Hypo QP, Evaluate output | Bulk→fault feedback check | `|tau1_trial| + |tau2_trial|` ∼ 10+ MPa post-peak; exactly-zero dip = §E.2 smoking gun | C-1 |
| W-1  | Fault face, off-hypo by 200 m along strike | First-hop amplitude | ‖Q‖∞ ≈ 1e−2 to 1e−1 (growing) | — |
| W-2  | Bulk element 1-hop normal to fault at hypo depth | Radiation escaping the fault | `|rhs[VX]| ∼ 3.7e5 m⁴/s²` per face QP accumulating; ‖Q‖∞ ≈ O(MPa)/O(m/s) | C-3 |
| W-2b | Post-RK4 bulk update at hypo-adjacent DOFs | Bulk Q genuinely growing | `max|Q[VX]|` monotonic → O(0.1 m/s) by t = 1.2 s | C-4 |
| W-3  | Free-surface element directly above hypocentre | Wave reaching free surface | Arrives at t ≈ 7.5 km / c_s ≈ 2.2 s (LATE; use tfinal ≥ 2.5 s) | — |
| W-4  | Partition-boundary fault face (any rank-pair seam) | H-V9-D double-count check | visit counter = 1 per face per step (MPI-Allreduce) | — |

**Smoking guns:**
- H-V9-A residual: F_h_total at W-0 has non-zero stress beyond 1 Pa·m/s.
- H-V9-C residual: F_h_total from interior vs shared branches differs at a `sign_flipped` face.
- H-V9-D: face-visit counter has min ≠ 1 or max ≠ 1 after Allreduce.
- H-V9-E: ‖Q‖∞ at W-2 grows then shrinks > 1% without boundary absorption.
- H-V9-F: W-3 saturates at ~1e-8 while W-0 is at physical magnitude — radiation arriving but damped.

**Face-visit counter (H-V9-D).**  Per-rank
`std::unordered_map<int64_t, int>` keyed on face-vertex key
(`shared_fault_key.hpp`).  Increment once per `Evaluate`.  At end of
step, MPI-Allreduce; assert max == min == 1.

**Acceptance (when resumed):**
- Routing sub-test: both sides of a `sign_flipped` face sum to zero within machine ε.
- Counter sub-test: exactly 1 visit per face per step globally.
- Radiation sub-test: W-2 and W-3 match the expected magnitudes.  If any watch point is ≥ 10¹⁰ below expected, pivot to the responsible hypothesis.

### B.3 Phase 3 — Frontera DIAG at production scale (frozen)

100 ranks first (smaller / cheaper), mesh 200 m, tfinal = 1.4 s, all
DIAG flags on (throttle to keep stdout < 1 GB).  Scale to 400 ranks
only if the 100-rank log is inconclusive.

**A3-PASS:** 100-rank run reproduces the pin (same location as job
7666327) and DIAG traces follow the §0.5 bisection.
**A3-FAIL:** 100 ranks doesn't show the pin → H-V9-D becomes prime;
re-submit at 400 ranks.

### B.4 Phase 5 — Frontera re-verification (frozen)

Only after Phases 1–4 merge.  Re-run job 7666327 exactly (400 ranks /
200 m / P1 / tfinal = 2.0 s).  Three pass bars:

- **Bar 1 (minimum):** V_max crosses 10 m/s and keeps growing past t = 1.4 s; no pin.
- **Bar 2 (target):** V_max peak within ±20 % of SeisSol reference; bulk ‖Q‖∞ physically reasonable (O(MPa)).
- **Bar 3 (final):** SCEC TPV102 receiver time-histories agree visually with benchmark dataset.

Only Bar 3 closes R-1001.

---

## Appendix C. DRDG3D ground-truth comparison (§0.4 carry-forward)

**Source plot:** `tpv102/plots_results_200m_p1_1.5s_400r_v2_job7665297/tpv102_flt_0_7.5.png`
— MFEM (black) vs DRDG3D (red dashed) at `(along_strike=0, down_dip=7.5 km)`
for `0 ≤ t ≤ 1.5 s` MFEM / `0 ≤ t ≤ 15 s` DRDG3D.

### C.1 Agreement during rise phase (t ≤ 1.0 s)

| Quantity | MFEM | DRDG3D | Agreement |
|---|---|---|---|
| τ_strike            | 75 → 98 MPa ramp      | 75 → 98 MPa ramp      | ≤ 1 % |
| V_strike onset      | takeoff at t ≈ 1.0 s  | takeoff at t ≈ 1.0 s  | match |
| log₁₀(State) drop   | +9 → −2 at t ≈ 1.0–1.2 s | same               | match |
| Normal stress       | 120 MPa flat          | 120 MPa (± 50 kPa)    | match |

### C.2 Deviation post-peak (t ≥ 1.2 s) — bug-driven

| Quantity | MFEM | DRDG3D | Diagnostic value |
|---|---|---|---|
| V_strike peak           | **pins at 7.69871 m/s** | peaks at 4.3 m/s then **decays** | MFEM ~80 % high at peak AND never decays |
| V_strike decay          | **none** (flat plateau) | ~4.3 → 0.5 m/s over 1–2 s | No radiation damping in MFEM |
| τ_strike post-peak      | 65 MPa constant       | 65 MPa with minor oscillations | Magnitude roughly matches |
| **V_dip**               | **exactly 0.000 m/s** | ± 0.003 m/s noise (radiation) | **SMOKING GUN** |
| **slip_dip**            | **exactly 0.000 m**   | ± 1e-4 m oscillations | **SMOKING GUN** |
| **τ_dip**               | **exactly 0.000 MPa** | ± 0.1 MPa oscillations | **SMOKING GUN** |
| Normal stress oscillations | 120.000 MPa constant | 120 ± 0.05 MPa oscillations | MFEM no oscillation |
| log₁₀(State) recovery   | stays at −2 (saturated) | −2 slowly → +0.5 by t = 15 s | MFEM state frozen |

### C.3 Physical interpretation

1. **Rise phase matches** because `τ_strike` ramps under nucleation
   `tau_0(t)` and `V_abs` is driven purely by `tau_0 + tau_trial` with
   `tau_trial ≈ 0` on quiescent bulk.  INDEPENDENT of bulk-coupling
   correctness.
2. **Post-peak deviation comes from radiation damping.**  In a correct
   code, bulk wave radiation produces `Q_plus[VY] ∼ O(m/s)`,
   `Q_plus[SXY] ∼ O(MPa)` near the fault.  Through Eq. 7b–c the
   bulk-induced `tau_trial` reaches O(−20 MPa), cancelling part of
   `tau_0` and reducing `V_abs` (radiation-damping / under-stress).
   MFEM's V_strike never decays ⇒ `tau_trial ≈ 0` throughout ⇒
   bulk Q at hypocentre is NOT receiving physical-scale radiation.
3. **The exactly-zero dip channels** are the most specific diagnostic.
   Any real rupture produces O(1e−3 m/s, 0.1 MPa) noise in the dip
   channels at the hypo QP.  MFEM shows exact zero because bulk Q is
   stuck at `≪ 1e−6` — consistent with measured `‖Q‖∞ = 1.66e−8`.

**Mapping to H-V9-J.**  The F_h = 0 cancellation applies identically
to the dip sub-block as to the strike sub-block, which explains the
V_dip / slip_dip / τ_dip exactly-zero signature.

---

## Appendix D. DIAG instrumentation specification (Phase 2 reference)

**Wall-clock cost:** ~30 min to add the four fprintf blocks + 1
unit-test driver; DIAG log volume stays < 5 MB per rank at tfinal =
1.2 s with the step gate active.

### D.1 Checkpoints C-1…C-4 (bisection of the fault→bulk chain)

| ID | Location | Expected at t ≈ 1.3 s | If wrong → conclusion |
|---|---|---|---|
| C-1 | `fault_face_flux.cpp::Evaluate` after line 82 | `|tau1_trial|+|tau2_trial| ∼ 10+ MPa` post-peak | Bulk→fault feedback broken UPSTREAM of Evaluate |
| C-2 | `wave_operator.inl::ComputeFaceFluxRHS` line ~812 (after `flux_.Interior`) | `|F_h_total[VX_global]| ∼ 2.2e2 m²/s²` | Evaluate correct but T_can+flux.Interior gives wrong F_h → H-V9-A/B/C |
| C-3 | Same file line ~818 (accumulation loop) | `|rhs[hypo-adjacent bulk VX]| ∼ 3.7e5 m⁴/s²` per face QP | F_h correct but accumulation blocked → H-V9-C/D |
| C-4 | `tpv102_driver.cpp::main` post-RK4 Q update | `max|Q[VX_global]|` monotonic → O(0.1 m/s) by t = 1.2 s | Bulk Q genuinely doesn't grow → mass-inverse over-attenuates OR filter zeroes Q |

**Bisection mental model:**
```
Evaluate ←(C-1)← bulk Q feedback ↔(C-2)↔ F_h injection ↔(C-3)↔ rhs accumulation ↔(C-4)↔ bulk Q growth
```

### D.2 MPI deadlock protection (mandatory design rule)

Any collective (MPI_Allreduce, MPI_Barrier, …) in the DIAG blocks
**must be called by EVERY rank, unconditionally**.  A deadlock occurs
if rank A enters Allreduce but rank B does not.

1. **C-1, C-2, C-3 emit printf only — no collectives.**  Gated on
   per-DOF `data.diag_print`; the owning rank prints, others skip.
2. **C-4 uses `MPI_Allreduce` for global max.**  Throttle
   `(step % 100 == 0) || (step % 10 == 0 && t ≥ 1.0 && t < 1.4)`
   depends only on `step` / `t` — identical on every rank — so every
   rank enters Allreduce together.
3. **Do NOT** gate C-4's Allreduce on `diag_hypo_dof >= 0` (true on
   at most one rank).
4. **Do NOT** gate C-4's Allreduce on `rank == 0` (instant deadlock).
5. **No collectives inside `FaultFaceFlux::Evaluate`** — local fault-face
   count varies per rank.

### D.3 Code snippets

**Preamble (header + driver):**
```cpp
// new header or driver .cpp
#ifdef SEAS_DIAG_FAULT_FLUX
extern int g_seas_my_rank;
#endif
```
```cpp
// tpv102_driver.cpp::main, after MPI_Comm_rank
#ifdef SEAS_DIAG_FAULT_FLUX
int g_seas_my_rank = 0;
// later:
g_seas_my_rank = rank;
#endif
```

**C-1** (`fault_face_flux.cpp::Evaluate` after `tau2_total = ...`):
```cpp
#ifdef SEAS_DIAG_FAULT_FLUX
if (data.diag_print) {
   std::fprintf(stderr,
      "[C-1 EVAL] rank=%d  tau1_trial=%+.3e Pa  tau2_trial=%+.3e Pa  "
      "|Q_plus[VY]|=%.3e  |Q_plus[SXY]|=%.3e  |Q_plus[VZ]|=%.3e  "
      "|Q_plus[SXZ]|=%.3e  psi=%.3e\n",
      g_seas_my_rank, tau1_trial, tau2_trial,
      std::abs(Q_plus[VY]),  std::abs(Q_plus[SXY]),
      std::abs(Q_plus[VZ]),  std::abs(Q_plus[SXZ]),
      data.psi);
}
#endif
```

**C-2** (after `flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h_total);`):
```cpp
#ifdef SEAS_DIAG_FAULT_FLUX
if (dof_idx >= 0 &&
    dof_idx < static_cast<int>(fault_dof_data_->size()) &&
    (*fault_dof_data_)[dof_idx].diag_print)
{
   real_t F_v_mag = std::sqrt(F_h_total[VX]*F_h_total[VX]
                              + F_h_total[VY]*F_h_total[VY]
                              + F_h_total[VZ]*F_h_total[VZ]);
   real_t F_s_max = 0.0;
   for (int c = 0; c < 6; c++) F_s_max = std::max(F_s_max, std::abs(F_h_total[c]));
   std::fprintf(stderr,
      "[C-2 FLUX] rank=%d  dof=%d  |F_v|=%.3e m2/s2  max|F_stress|=%.3e Pa·m/s  "
      "F_h[VX]=%+.3e  F_h[SXY]=%+.3e  F_h[VY]=%+.3e  F_h[SXZ]=%+.3e\n",
      g_seas_my_rank, dof_idx, F_v_mag, F_s_max,
      F_h_total[VX], F_h_total[SXY], F_h_total[VY], F_h_total[SXZ]);
}
#endif
```

**C-3** (bracketing the `rhs[...] -= w*shape1(i)*F_h_total[c];` loop):
```cpp
#ifdef SEAS_DIAG_FAULT_FLUX
const bool _c3_diag = (dof_idx >= 0 &&
                       dof_idx < static_cast<int>(fault_dof_data_->size()) &&
                       (*fault_dof_data_)[dof_idx].diag_print);
const int  _c3_probe = VX * ndof_total_ + dof_offset1 + 0;
const real_t _c3_pre = _c3_diag ? rhs[_c3_probe] : 0.0;
#endif

for (int c = 0; c < NUM_STATE; c++)
   for (int i = 0; i < ndof; i++)
      rhs[c * ndof_total_ + dof_offset1 + i] -= w * shape1(i) * F_h_total[c];

#ifdef SEAS_DIAG_FAULT_FLUX
if (_c3_diag) {
   const real_t _c3_post = rhs[_c3_probe];
   std::fprintf(stderr,
      "[C-3 RHS]  rank=%d  dof=%d  rhs[VX,elem1,dof0] pre=%+.3e post=%+.3e  "
      "delta=%+.3e  w=%.3e  shape1(0)=%.3e  F_h_total[VX]=%+.3e\n",
      g_seas_my_rank, dof_idx, _c3_pre, _c3_post,
      _c3_post - _c3_pre, w, shape1(0), F_h_total[VX]);
}
#endif
```

**C-4** (`tpv102_driver.cpp::main` after RK4 update):
```cpp
#ifdef SEAS_DIAG_FAULT_FLUX
{
   const bool c4_fire = (step % 100 == 0) ||
                        (step % 10 == 0 && t >= 1.0 && t < 1.4);
   if (c4_fire) {
      real_t q_vx_local = 0.0, q_sxy_local = 0.0;
      for (int i = 0; i < ndof_total; i++) {
         q_vx_local  = std::max(q_vx_local,
                                std::abs(Q[VX  * ndof_total + i]));
         q_sxy_local = std::max(q_sxy_local,
                                std::abs(Q[SXY * ndof_total + i]));
      }
      real_t q_vx_global  = q_vx_local;
      real_t q_sxy_global = q_sxy_local;
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&q_vx_local,  &q_vx_global,  1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
      MPI_Allreduce(&q_sxy_local, &q_sxy_global, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
      if (rank == 0) {
         std::fprintf(stderr,
            "[C-4 BULK] step=%d  t=%.4f  V_max=%.3e  "
            "max|Q[VX]|=%.3e m/s  max|Q[SXY]|=%.3e Pa\n",
            step, t, V_max_step, q_vx_global, q_sxy_global);
      }
   }
}
#endif
```

### D.4 Expected DIAG output (correct code at t ≈ 1.3 s)

```
[C-1 EVAL] rank=3 tau1_trial=+3.21e+04 Pa  tau2_trial=-1.75e+07 Pa  ...
[C-2 FLUX] rank=3 dof=... |F_v|=2.18e+02 m2/s2  F_h[VX]=-2.18e+02 ...
[C-3 RHS]  rank=3 dof=... pre=+1.23e+05 post=+4.92e+05 delta=+3.69e+05 ...
[C-4 BULK] step=4562 t=1.3020 V_max=4.31e+00 max|Q[VX]|=2.84e+00 max|Q[SXY]|=1.54e+07
```

**On v8.0.0 buggy code:** C-4 reports `~1.66e-8` instead of `~2.84`.
The FIRST of C-1/C-2/C-3 to report a wrong magnitude tells us which
upstream link broke.  Under **H-V9-J**, C-1 is correct (Evaluate
receives and produces correct values) and C-2 shows F_h ≈ 0 — exactly
as we observed in the unit test.

### D.5 Unit-test cover for the MPI deadlock rule

```cpp
// tests/parallel/test_c4_diag_no_deadlock.cpp
TEST_C4_DIAG_CompletesWithoutDeadlock: {
   // 2-rank ParMesh + WaveOperator + driver loop
   // Compile with -DSEAS_DIAG_FAULT_FLUX.
   // Run 300 RK4 steps to exercise all three throttle branches.
   // Wrap mpirun in 30 s timeout.  Hang == deadlock.
}
```

---

## Appendix E. v8.0.0 carry-forward

### E.1 Closed (do NOT revisit)

| ID | Finding | Closure artefact |
|---|---|---|
| R-101     | Shared-fault DOFData verifier pairs cleanly across METIS topologies | v8.0.0 Step 1.10–1.14 |
| R-101-a   | 50-rank job 7666340 — 123 pairs, 0 unpaired, max_rel_diff = 3.33e-16 | job 7666340 |
| R-101-b   | 400-rank job 7666327 — 405 pairs, 0 unpaired, max_rel_diff = 6.66e-16 | job 7666327 |
| BP5 path  | BP5 convention (tangent1=dip, tangent2=strike) propagated throughout dynamic/ | Step 1.2 / R-801 |
| CFL       | Explicit-RK4 picks dt = 2.85e-4 s at CFL = 0.5; no impact on pin | job 7666327 |
| Build     | `SEAS_DIAG_FAULT_FLUX` / `SEAS_DIAG_GHOST_EXCHANGE` build flags in place | Step 1.8 |
| Mass op   | `AssembleElementMassInverse` (`wave_operator.inl:1254-1282`) correct for first-order conservation form | REVIEW R-1003 |

### E.2 Open (v9.0.0 scope; mapped onto H-V9-J)

| ID | Symptom | Source |
|---|---|---|
| **R-1001**   | V_max saturates at 7.69871 m/s from t ≈ 1.30 s onward, bit-identical subsequent steps | job 7666327 |
| R-1001-b   | bulk max ‖Q‖∞ grows to 1.66e-8 — 10¹⁵ below physical | job 7666327 qnorm |
| R-1001-c   | qnorm at hypo rank r3 > other watched ranks; METIS IDs NOT spatial neighbours | job 7666327 qnorm watch |
| R-1001-d   | R-101 verifier clean at 400 ranks while R-1001 reproduces — rules out pairing | job 7666327 stdout |

All four symptoms are consistent with H-V9-J (see §1 mapping table).

### E.3 Ruled out (do NOT re-diagnose)

- **DOFData cross-rank inconsistency** — verifier clean at 400 ranks.
- **CFL / dt over-drive** — pin survives at CFL = 0.5; dt auto-selected.
- **Rank count as root cause** — pin reproduces at 400 ranks; 50-rank only passed because tfinal=0.01 s didn't reach breakaway.
- **Mass-matrix scaling** — audit shows formula correct; 10¹⁵ attenuation not achievable via mass scaling.
- **Friction solver / Evaluate / nucleation ramp** — rise phase matches DRDG3D (Appendix C.1); bug is in post-peak radiation feedback.

---

## Appendix F. Out of scope for v9.0.0

Defer unless triggered by the chosen H-V9-J fix failing:

- Mesh refinement to h = 100 m (H-V9-G).
- Order bump to P2 (H-V9-G).
- SeisSol convention audit (H-V9-H).  Blocked on R-703 Slot B.
- PML layer interaction audit.  Separate ticket if Phase 2 hints at it.
- Quasi-dynamic / dynamic merge path cleanup (policy: duplicate BP5,
  do not edit BP5 — see feedback memory).
- `BuildRotation` / `BuildRotationInverse` stress-block Voigt factor
  for non-axis-aligned canonical frames (per REVIEW Unreviewed Areas).
  Only open if rotated-mesh testing is requested.

---

## Appendix G. Revision history

- **Rev 5 (2026-04-20)** — Plan reorganized per user request for a
  dashboard-first layout.  Same content, restructured: Dashboard +
  §1–§3 (H-V9-J, decision, regression gate) up top; previous Phase 0 /
  Phase 2 / Phase 3 / Phase 5 / DIAG / DRDG3D / carry-forward
  content demoted to labeled appendices.  **No content changes.**
- **Rev 4 (2026-04-20)** — Phase 1 executed.  A1-FAIL on §3.1b,
  F_h[VX_global] = 0 exactly.  Added H-V9-J as new rank-1 hypothesis;
  falsified H-V9-A / -B / -C / -I at unit level.  Test fixture
  corrected from V_ini=1e-12 equilibrium to post-breakaway
  steady-slip operating point (V=1 m/s, tau2_0=81 MPa).  Updated
  hypothesis table, Phase 4 fix options (Option A vs B), decision
  rules.
- **Rev 3 (2026-04-20)** — Applied user feedback (no local reproducer
  runs) + DRDG3D smoking-gun analysis (Appendix C):
  - Phase 0 rewritten as Frontera sbatch generation (no local runs).
  - §3A fast-path FP-5 moved to Frontera-submitted sbatch.
  - Promoted "exactly-zero dip channel" to highest-value diagnostic.
  - C-1…C-4 DIAG checkpoints folded into Phase 2 spec.
- **Rev 2 (2026-04-20)** — Applied REVIEW.md findings R-1001…R-1007:
  - Promoted H-V9-C to rank 1 pre-execution (interior-vs-shared asymmetry).
  - Added H-V9-I (Godunov identity under BuildFrame tangent arbitrariness).
  - Demoted H-V9-E per mass-matrix audit.
  - Added §3.1b end-to-end chain fixture; §3A fast-path sequence.
  - Expected-value column in watch-point table; MPI-Allreduced face-visit counter.
- **Rev 1 (2026-04-20)** — Initial draft.

---

## Appendix H. References

- `tpv102_debug_v9.0.0_fix.md` — Phase 1 execution log, H-V9-J
  statement, Option A vs B detail.
- `REVIEW.md` (2026-04-20) — Rev 2 input; findings R-1001…R-1007.
- `tpv102_debug_v8.0.0_debug_plan.md` — Phase 2–4 specs inherited.
- `tpv102_debug_v8.0.0_debug_fix.md` — Phase 1 execution record.
- `tpv102_debug_v8.0.0_check.md` — R-1001 finding, hypotheses H1–H4.
- Frontera job 7666327 stdout — 400-rank reproducer reference.
- Frontera job 7666340 stdout — 50-rank R-101 closure reference.
- `miniapps/seas/dynamic/wave_operator.inl:697-866` — interior-fault branch.
- `miniapps/seas/dynamic/wave_operator.inl:1155-1194` — shared-fault branch.
- `miniapps/seas/dynamic/fault_face_flux.cpp:60-155` — Evaluate / Eq. 10 block.
- `miniapps/seas/dynamic/fault_basis.hpp:340-452` — `ComputeOrientedFrame` + `sign_flipped`.
- `miniapps/seas/dynamic/godunov_flux.cpp:125-188` — `BuildJacobian` A_x (audited correct).
- `miniapps/seas/dynamic/godunov_flux.cpp:238-300` — `BuildRotation` + `BuildFrame`.
- `miniapps/seas/dynamic/wave_operator.inl:1254-1282` — `AssembleElementMassInverse` (audited correct).
- `miniapps/seas/dynamic/tpv102_setup.hpp:40-100` — BP5 convention init.
- `miniapps/seas/dynamic/shared_fault_key.hpp` — face-vertex key utility.
- `miniapps/seas/drivers/tpv102_driver.cpp:886-901` — qnorm watch ranks (METIS IDs, NOT spatial).
- de la Puente et al. 2009 — "Dynamic rupture modeling on unstructured meshes using a DG method".
