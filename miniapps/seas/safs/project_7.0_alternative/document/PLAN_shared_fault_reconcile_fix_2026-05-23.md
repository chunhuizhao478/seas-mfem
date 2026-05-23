# Implementation Plan: shared-fault cross-rank DOFData consistency (SAFS)

## TL;DR — there are TWO separate problems; this plan is about the second one

| | Problem A — the blow-up | Problem B — the cross-rank inconsistency |
|---|---|---|
| **What** | The fault +/- side was mislabeled on the tilted (curvilinear) fault, so the nucleation pushed with the wrong sign → V_max → 10³ blow-up. | The same fault quad-point is solved **twice** (once per MPI rank) from inputs that differ at the 16th digit; at the rupture onset that hair's-width difference makes the two ranks disagree (slip vs lock). |
| **Scope** | Geometry (which side is "+"). Independent of friction law / dip / strike. | MPI / floating-point. Triggered the same way for any channel; only **runs away** when there is dip slip. |
| **Status** | **FIXED** — commit `800b3281` (see Part A). Validated by the `zerodip` job. | **OPEN** — this plan. |
| **Test** | `spatial_dyn_zerodip_8N_400r_dev_2hr_safs.sbatch` (no dip slip, guard non-fatal). | New local injection test + the Dc2 `XRANK` Frontera run. |

Documents behind this plan: diagnosis `spatial_dynamic_rupture_speckle_blowup_2026-05-22.md`,
design review `spatial_dynamic_rupture_reconcile_review_2026-05-23.md`.

---

## Part A — the blow-up fix (DONE, commit `800b3281`) — documented for the record

> **This part is already landed and is NOT implemented by this plan.** It is
> documented here because the user asked to keep it separate and because it is
> the reason the `zerodip` run is now stable.

### What was wrong (plain language)
On each fault face the code must decide which of the two touching elements is on
the "+" side and which is on the "−" side (this sets the sign of the nucleation
push and the flux). It decided by measuring which element sits "below" a FIXED
reference direction `ref_normal = (0,−1,0)` — i.e. it projected the
element-to-face offset onto that hardcoded direction. That works for a flat
fault whose normal IS (0,±1,0). The meshed SAFS fault is **curvilinear**: where
the fault swings to a roughly N–S strike, its normal becomes nearly
perpendicular to `(0,−1,0)`, so the projection measures mostly the *sideways*
part of the offset instead of the real "above/below." Past a tilt of ≈71.6° the
sideways part wins and the side label **flips** — both ranks of a shared face
then call themselves "−", the "exactly one + side" rule breaks, the nucleation
is applied with the wrong sign, and the rupture blows up.

### The exact computation (old vs new)
The side label is set in `wave_operator.inl` (interior `:459-465`, shared
`:519-525`):
```
elem1_proj = elem1_c · ref_normal       face_proj = face_c · ref_normal
elem1_on_plus = (elem1_proj < face_proj)   ⇔   (elem1_c − face_c) · ref_normal < 0
```
Let `d = elem1_c − face_c` be the vector from the face centroid to the owning
element's centroid. Split it into the part along the TRUE face normal `n̂` and the
in-plane (tangential) part: `d = d_n·n̂ + d_t`, where by definition `d_t · n̂ = 0`.
The OLD test computes `d · ref_normal`:
```
d · ref_normal = d_n (n̂ · ref_normal) + (d_t · ref_normal)
```
- **Flat fault:** `n̂ = ±(0,1,0) = ±ref_normal`, so `n̂·ref_normal = ±1` and the
  tangential term `d_t·ref_normal = 0` (because `d_t ⟂ ref_normal` too). Result:
  `d·ref_normal = ±d_n` — it measures exactly how far the element sits
  above/below the face. **Correct.**
- **Tilted fault (face normal makes angle θ with `ref_normal`):**
  `n̂·ref_normal = cos θ`, and the tangential part now projects with weight `~sin θ`.
  A tetrahedron's centroid is NOT directly over its face centroid, so `d_t ≠ 0`.
  Result: `d·ref_normal ≈ d_n cos θ + (d_t·ref_normal)`. As `θ → 90°` (N–S strike,
  `n̂ ⟂ ref_normal`), `cos θ → 0` and the **tangential** term takes over. When
  `|d_t·ref_normal| > |d_n cos θ|`, the SIGN is set by the sideways offset, not the
  real above/below → **wrong label.**
- On the 2-tet test fixture this margin is `d·ref_normal = kL·(sin θ/12 − cos θ/4)`,
  which crosses zero at `tan θ = 3 → θ = 71.6°`: past that tilt the label flips,
  both ranks of the shared face get the same sign on `(elem1_c−face_c)·ref_normal`
  → both call themselves "−" → "exactly one +" breaks → the fixed-sign nucleation
  vector is added to the wrong side → blow-up.

### What we changed (commit `800b3281`)
Project onto the ACTUAL face normal `n̂` (canonicalized so both ranks agree on its
direction), not the hardcoded `ref_normal`:
```
CalcOrtho(face Jacobian) → n_raw
if NormalNeedsFlipToCanonical(n_raw, ref_normal, dim): n_raw = −n_raw   // → n̂
elem1_on_plus = ((elem1_c − face_c) · n̂ < 0)
```
Because `d_t ⟂ n̂` BY CONSTRUCTION, `d · n̂ = d_n` **exactly** — the tangential
offset can never contaminate the sign, at any tilt. The decision margin is now
`±|d_n| ~ O(element size)`, independent of θ, so it never flips.
- `fault/fault_basis.hpp` — new shared static
  `FaultBasis::NormalNeedsFlipToCanonical(n_raw, ref_normal, dim)`: orient the
  normal LINE so `n̂·ref_normal > 0`; in the degenerate band where `n̂ ⟂ ref_normal`
  (the θ=90° case where even that sign is ill-defined) fall back to "largest-
  magnitude component positive." Both ranks then compute the SAME `n̂` regardless
  of which way `CalcOrtho` happened to point it. `ComputeOrientedFrame` routes its
  sign decision through the same helper.
- `dynamic/wave_operator.inl` — both side blocks compute `n̂` via `CalcOrtho` and
  project the centroid offset onto it instead of `ref_normal`.
- For a flat fault `n̂ = ±ref_normal` so the boolean is identical to before
  (TPV/BP5 byte-exact); the change only differs where `θ > 0`.

### How Part A is validated (per the user's directive)
Job script: `jobs/safs/spatial_dyn_zerodip_8N_400r_dev_2hr_safs.sbatch`.
- `SEAS_ZERO_DIP_PRESTRESS=1` (default) — zeroes the dip pre-stress at the source
  so there is **no dip slip**; the fault is loaded in strike only.
- `SEAS_R101_NONFATAL=1` — the cross-rank guard logs `worst_rel` each check
  instead of aborting, so **Problem B does not force-abort the Problem-A run**
  (this is the "do not force guard the MPI DOFData inconsistency" requirement;
  the mechanism already exists, commit `fc28454`, and the call site is correct:
  `spatial_dyn_driver.cpp:1994` passes `tol, &worst_rel, &worst_field,
  abort_on_fail=!nonfatal`).
- Result (job 7747119, MEASURED): the run advances with **no V_max blow-up**
  (`V_max` peaks ~6.3 m/s then decreases) — **Part A is fixed.** Crucially, the
  R-101 non-fatal log still goes `worst_rel → 1.0` at t≈0.595 s (in the STRIKE
  fields V2/slip2): so Part A removed the *blow-up*, but **issue B (the cross-rank
  speckle) is still present** — direction-blind, just non-catastrophic without dip.

**Acceptance for Part A (met):** `V_max` physical (no 10³, peaks ~6 then decays).
The persistent `worst_rel → 1.0` is **issue B** (confirmed active); the speckled
slip-rate elements coincide with it *in time* (working hypothesis: the speckle is
issue B, to be confirmed/refuted by the Phase-4 post-fix rerun) — either way NOT
a Part-A regression.

---

## Part B — the cross-rank DOFData inconsistency (THIS PLAN)

### B.1 The problem in plain language

A fault quadrature point that lies on the boundary between two MPI ranks is a
**shared** point: rank A owns the element on one side, rank B owns the element
on the other. To advance the fault friction there, each rank needs the stress on
**both** sides. Each rank already has its own side exactly, and gets the other
side as a **ghost copy** from its neighbor. So both ranks have "both sides" and
both run the friction solve for that point. **The same physics is computed
twice, once on each rank.** (This is what R-701 set up: it deleted an older
"one rank computes, tells the other" broadcast, betting that both ranks would
get identical answers anyway.)

That bet is wrong by a hair. The two ranks do **not** feed the friction solve
bit-identical numbers. The reason is pure floating-point bookkeeping, not a
logic bug:

- The shared point is one physical location, but rank A reaches it through its
  element's local coordinates and rank B reaches the *same* point through the
  *neighbor* element's local coordinates. The interpolation weights
  (`shape × DOFs`) are mathematically equal but rounded differently, so the bulk
  stress the two ranks plug in differs at about the **16th significant digit**
  (~1e-14). (Confirmed on Frontera: the two ranks' inputs match to 11 printed
  digits and differ below that.)

Normally a 1e-14 difference in → a 1e-14 difference out: invisible. The trouble
is the **rupture onset is a yes/no switch**. The slip rate is
`V = max(0, (|τ| − strength)/η_s)`: below the strength the point is **locked**
(`V=0`); above it, it **slips**. Right at the threshold `|τ| = strength`, a
hair's-width difference in `|τ|` flips the answer between "locked" and "slips."
The SAFS nucleation ramps the stress up **slowly**, so the point sits **right on
the knife's edge for several sub-steps** — and during that window the 1e-14
difference makes **rank A slip while rank B stays locked** (or vice-versa). From
that instant the two ranks are computing two different earthquakes at that point;
the slip-weakening feedback then drives them apart, and the run is corrupted.

A runtime guard (the R-101 check) compares the two ranks each step and aborts
when they disagree — that is what stops the SAFS run at t≈0.455 s.

**Why benchmarks never caught it:** TPV102/104 use rate-and-state friction,
whose slip rate is a **smooth** function (no yes/no switch) — a 1e-14 input stays
a 1e-14 output forever, so the ranks never visibly disagree. TPV205 uses the
same slip-weakening switch as SAFS, but its rupture front sweeps past a point in
one step (no slow knife's-edge dwell). Only SAFS combines the slip-weakening
switch **and** a slow nucleation that parks a shared point on the threshold.

### B.1.5 Root-cause status — what is CONFIRMED vs what is still INFERRED

Be honest about the chain, because it decides which fix is correct.

**Confirmed (Frontera XRANK trace, job 7747036):**
- The shared-QP friction state is computed **redundantly** on both ranks.
- At the diverging QP the two ranks' inputs match to **11 printed digits** yet the
  friction outputs **split** (one slips, one locks) at the `max(0,·)` kink. Since
  the solve is deterministic, identical-to-11-digits inputs producing different
  outputs **proves there is a sub-11-digit cross-rank input difference**, and that
  the non-smooth onset amplifies it. This is the *mechanism* of the divergence.

**NOT yet measured (the open part of the root cause):**
- The **magnitude** of that input difference ("~1e-14" is an assumption, never
  measured).
- Its **source**, which is one of:
  1. **Interpolation:** the same physical QP is evaluated through two different
     element parametrisations (`shape1∘Loc1` on rank A vs `shape2∘Loc2` on rank
     B's ghost) — mathematically equal, differently rounded. Inherent to DG
     shared faces; unavoidable.
  2. **Ghost-exchange not bit-exact:** rank B's ghost copy of rank A's element
     stress should be bit-identical to rank A's own, but isn't. A real defect —
     potentially fixable AT THE SOURCE.
  3. **Accumulated bulk divergence:** the distributed bulk solve drifts across
     ranks over the run (different reduction/matvec order). Generally unavoidable.

**Why this matters for the fix (decision gate):**
- Source = (1) or (3) → the gap is inherent → the **reconcile** (Phase 2) is the
  correct, general fix.
- Source = (2) → make the ghost copy bit-exact → both ranks get identical inputs
  → identical outputs → **no reconcile needed** (this would *vindicate* R-701's
  no-broadcast design and be cheaper/more local).

**Phase 0 (below) measures the magnitude and localises the source, and chooses
between the two fixes. We do not commit to the reconcile until Phase 0 says the
seed is inherent.**

### B.2 Why the DIP direction is the one flagged — what the equations actually say

I earlier hand-waved "dip couples to normal stress." You asked for equations, and
they are worth being exact about, because **the local fault solve is dip/strike
SYMMETRIC** — there is no built-in per-direction instability, and any claim that
"dip is unstable, strike is stable" is wrong at the solve level.

**The fault-local solve.** In the canonical (fault-aligned) frame the bulk state
splits into THREE independent Riemann channels — normal (local X), dip (Y=t1),
strike (Z=t2) — `ComputeTrialTraction`, `fault_face_flux.cpp:49`:
```
σ_n_trial = η_p (v_n⁻ − v_n⁺ + σ_n⁺/Zp⁺ + σ_n⁻/Zp⁻)        (normal,  Eq.7a)
τ1_trial  = η_s (v_t1⁻ − v_t1⁺ + τ1⁺/Zs⁺ + τ1⁻/Zs⁻)         (dip,     Eq.7b)
τ2_trial  = η_s (v_t2⁻ − v_t2⁺ + τ2⁺/Zs⁺ + τ2⁻/Zs⁻)         (strike,  Eq.7c)
```
Friction (`SolveLSW_TPV205`, `tpv205_friction.hpp:146-165`):
```
|τ|    = sqrt(τ1_total² + τ2_total²)        τ_str = μ_eff·|σ_n_total|
V_abs  = max(0, (|τ| − τ_str)/η_s)
V1 = V_abs·τ1_total/|τ|        V2 = V_abs·τ2_total/|τ|
τ1_corr = τ1_trial − η_s·V1    τ2_corr = τ2_trial − η_s·V2
σ_n_corr = σ_n_trial            ← normal traction NOT changed by slip (:202)
```
Three facts read directly off these:
1. **Dip (1) and strike (2) are identical in form** (same `η_s`, same algebra).
   The solver cannot prefer one over the other — there is NO per-direction
   instability here.
2. **Slip onset is a switch on the magnitude** `|τ|` (the `max(0,·)`). The 1e-14
   cross-rank difference in `|τ|` at the threshold flips `V_abs` between 0 and >0,
   so **`V1` AND `V2` flip across ranks together** (they share `V_abs` and `|τ|`).
   The trigger corrupts the WHOLE slip vector, not one direction.
3. **`σ_n_corr = σ_n_trial`:** slip does NOT change the normal stress *inside the
   solve*. And friction caps only the SHEAR (`|τ_shear| ≤ μ|σ_n|`); there is no
   cap on the normal channel (no "no-opening" limit — this is why the original
   blow-up reached `σ_n = −124 GPa` tension). The normal channel is the
   **unconstrained** one.

**So why is `V1` (dip) the field the guard reports, and why does the oblique run
blow up while zerodip doesn't?** Two parts, kept honest:

- **(reporting) Why `V1` is named.** The guard uses
  `rel_diff = |a−b| / max(|a|,|b|)`. When one rank slips and the other locks, BOTH
  `V1` and `V2` jump to `rel_diff = 1.0`; the verify reports the first field that
  hit the max, and `V1` is checked before `V2`. So "dip" is partly a **reporting
  artifact** — strike desynced too. (Equation-backed: same `rel_diff` formula.)

- **(physics, partly hypothesis) Why the blow-up needs dip.** This is NOT
  derivable from the local solve (which is symmetric, above). It comes from the
  BULK + geometry, and the zerodip job is the experiment built to test it:
  - The local channels are decoupled, but in the BULK the imposed slip velocity
    (`BuildImposedState`, `:254-280`: `v_imp^± = v^± ± (1/Z)(τ_corr − τ^±)`,
    `σ_imp = σ_corr`) radiates and propagates. On a **32°-dipping** fault, dip
    slip has a vertical component that the bulk elastodynamics couple into the
    fault-**normal** stress next step (`σ_n_trial` changes); strike slip (motion
    along the surface) does not. So a dip inconsistency leaks into the
    **unconstrained** normal channel (fact 3) and runs away; a strike
    inconsistency does not feed an un-capped channel and stays bounded.
  - **Empirical result (zerodip job 7747119, MEASURED):** zeroing the dip
    pre-stress separates the two effects cleanly:
    - **Blow-up is gone** — `V_max` stays 0 to t=0.525 s, peaks at **~6.3 m/s**
      (t=0.66 s), then *decreases* to ~3. A healthy rupture, no runaway. So the
      dip load IS what drives the *catastrophic blow-up* (its absence → no
      blow-up), consistent with the dip↔normal-coupling hypothesis (which lives
      in the bulk, not the flux equations; this run is its test).
    - **But the cross-rank inconsistency is STILL THERE** — the R-101 non-fatal
      log shows `worst_rel` at machine precision (≤1.8e-14) through step 1500,
      then jumping to **0.13 (step 1600, t=0.56 s) → 1.0 (step 1700)**, in fields
      **V2 (strike velocity) and slip2 (strike slip)**. So the inconsistency
      diverges in STRIKE here — **it is direction-blind** (corrects my earlier
      "strike works"). It coincides *in time* with the **speckled slip-rate
      elements outside the rupture** seen in ParaView, so the WORKING HYPOTHESIS
      is that the speckle is issue B. Not yet proven *spatially* (`worst_rel` is
      the max over shared QPs; the speckle is scattered) — **the Phase-4 post-fix
      rerun is the confirmation**: reconcile clears the speckle ⇒ it WAS issue B;
      speckle survives at `worst_rel≈0` ⇒ a co-present cause (under-resolution
      `L_nuc/h_min≈3–10 < 10`, or SSO) the reconcile won't fix.
    - **Net:** the cross-rank inconsistency (issue B) is present for any slip
      direction; whether it then BLOWS UP depends on the dip↔normal runaway.
      zerodip keeps it bounded (speckle ~O(1) m/s, no blow-up). Whether the
      OBLIQUE (dip-loaded) run blows up or *also* stays bounded after Part A is
      the **OPEN test** — exactly what `spatial_dyn_resDc2_XRANKdiag_*.sbatch`
      (now non-fatal) checks. The pre-hard-gate Dc2 run (7746601) aborted at the
      R-101 guard ~t=0.455 s, *before* the ~t=1.0 s blow-up window, so we never
      actually observed whether it blows up.

**Bottom line for the fix.** Because the local solve is direction-symmetric and
the trigger corrupts the whole slip vector, the fix must make the slip decision
single-valued for ALL directions — which it does. We do NOT rely on strike being
self-stabilizing or on resolving the bulk dip-vs-normal question; the reconcile
removes the inconsistency before it can feed any channel.

### B.2.5 Scope — do TPV102/104/205 have this too? And the mixed-flux connection

**Q1: yes, all three share the substrate, but only LSW is at risk.** The redundant
per-rank shared-QP solve from ~1e-14-different inputs exists for EVERY config at
np>1. Only the AMPLIFIER differs:
- **TPV102 / TPV104 (rate-and-state):** the solve `V = solve(a·asinh(V·C) =
  τ/σ_n)` (Brent, `friction_solver.cpp:34`) is **smooth & monotonic** — no
  `max(0,·)` switch. A 1e-14 input → ~1e-14 output forever, so the cross-rank
  inconsistency is PRESENT but pinned at the noise floor (`worst_rel ~1e-14 ≪
  1e-10`); never diverges, invisible. **Latent but harmless.**
- **TPV205 (LSW):** SAME `V = max(0,(|τ|−τ_str)/η_s)` switch as SAFS → amplifier
  present → **latent-at-risk.** Not triggered in standard TPV205 because its front
  sweeps each point in ~one step (no slow dwell at the threshold). A TPV205-class
  config with a slow nucleation parking a shared QP on the threshold would diverge.
- **SAFS (LSW + slow `gradual_overstress`):** switch AND slow dwell → **triggered**
  (confirmed: zerodip `worst_rel → 1.0`).

→ This is exactly why the reconcile must be **method-invariant**: the substrate is
shared by all three; gating to LSW would leave rate-state latently inconsistent
and would break the day a non-smooth rate-state cap is added.

**Q1, mixed-flux connection (your memory is half-right).** The mixed-flux fix
(`MIXED_FLUX_PLAN.md`) addressed a DIFFERENT speckle cause — **Spatial Spike
Oscillations**: grid-scale ringing from the ADER predictor under-dissipating
high-frequency modes at **p≥2** (a numerical-DISSIPATION deficit, *consistent*
across ranks). Different axis from this bug:
- SAFS runs at **p=1** (`fe order: 1`), where the SSO deficit is weak, AND an
  SSO/dissipation speckle would be cross-rank *consistent* (`worst_rel` small).
  Since SAFS shows `worst_rel → 1.0`, the speckle here is the **cross-rank**
  issue, not SSO.
- BUT the mixed-flux work flagged the SAME cross-rank PRINCIPLE
  (`MIXED_FLUX_PLAN.md` R5: "shared-face dispatch must apply mixed-flux uniformly
  across rank seams"; suggested an abort guard for `nprocs>1 AND use_substep AND
  mixed_flux≠none`). The reconcile is the general form of that principle, applied
  to the friction DOFData + imposed state instead of the flux mode. So it IS the
  same FAMILY (cross-rank seam consistency); mixed-flux just fixed a different
  symptom (dissipation) on a different field.

### B.3 How the proposed fix works in plain language

Stop computing the shared point twice and hoping the two answers match. Instead,
**pick one rank to be the boss for that point; the boss computes the friction
solve; it tells the other rank the answer; both ranks use the boss's answer.**
Now there is only ONE answer, so the two ranks cannot disagree — the yes/no slip
switch is decided once, by the boss, and copied. The 1e-14 input difference is
irrelevant because the non-boss never uses its own solve.

Concretely (still plain):
1. Both ranks do their friction solve as today (cheap; we don't restructure the
   solver).
2. They **exchange** the results for every shared fault point (a small MPI
   message — only the points on rank boundaries).
3. For each point, the **non-boss overwrites its own answer with the boss's** —
   both the friction state (slip rate, corrected traction, …) **and** the
   "imposed stress" that gets pushed back into the bulk wave solve. (We must copy
   the imposed stress too, not just the stored numbers — otherwise the bulk would
   re-grow a fresh 1e-14 difference next step and we'd be back where we started.)
4. THEN both ranks assemble their side of the flux from that single shared
   answer.

**How the boss is chosen (no negotiation round-trip).** Every shared fault point
appears on exactly two ranks. The reconcile exchange (the R-101 matcher, below)
pairs the two records by their face-vertex key and hands each rank its PEER's
rank number. Both ranks then independently compute the SAME value:
```
boss = min(my_rank, peer_rank)
is_boss = (my_rank == boss)        // exactly one of the two ranks is true
```
The boss keeps its own answer; the non-boss (`my_rank > peer_rank`) overwrites.
"Lower rank wins" is unique, deterministic, and **purely numeric** — so, unlike
the +/- side label that Part A had to fix, it cannot be fooled by the curvilinear
geometry. (Any fixed unique key would do; lower-rank is the cheapest because the
pairing already carries `peer_rank`. We `MFEM_VERIFY` that each shared point pairs
with exactly one peer.)

**What is copied — and why the imposed stress must be too.** The boss sends, per
point: (a) the friction state
`DOFData{V1,V2,τ1_corr,τ2_corr,σ_n_corr,slip1,slip2}` and (b) the canonical-frame
imposed state `Q_imp_plus`, `Q_imp_minus`. The imposed state is what gets pushed
back into the bulk (`BuildImposedState`, `fault_face_flux.cpp:254-280`):
```
Q_imp^±[VY] = v_t1^± ± (1/Zs)(τ1_corr − τ1^±)     Q_imp^±[SXY] = τ1_corr
Q_imp^±[VZ] = v_t2^± ± (1/Zs)(τ2_corr − τ2^±)     Q_imp^±[SXZ] = τ2_corr
Q_imp^±[VX] = v_n^±  ± (1/Zp)(σ_n_corr − σ_n^±)   Q_imp^±[SXX] = σ_n_corr
```
If we copied only the stored `DOFData` but each rank assembled the flux from its
OWN `Q_imp` (built from its own 1e-14-different `τ^±`/`v^±`), the bulk would
re-grow a fresh cross-rank difference next sub-step and we'd be back at the start.
Copying the boss's `Q_imp` makes the flux that BOTH ranks assemble bit-identical.

**Why this is the right kind of fix (general, not a patch):** it does the same
thing for every friction law and every slip direction — it just makes the shared
point single-valued. It does not tune a threshold, special-case the dip channel,
or lean on strike being lucky. (An older "owner-broadcast" doing exactly this
existed and was deleted by R-701; we are restoring it in the path R-701 left
uncovered, now that we've proven R-701's "both ranks agree anyway" assumption is
false.)

**One honest consequence (needs your OK):** making the non-boss adopt the boss's
value changes TPV102/104 by ~1e-14 (they currently differ across ranks by that
much; after the fix they're bit-identical across ranks). So we trade
"bit-identical to the old binary" for "physically identical (≤1e-13) AND now
bit-identical across ranks" — a stronger correctness property, but it does
relax the strict byte-exact rule. Flagged again here.

---

## Constraints
- **Method- and direction-invariant:** the reconcile runs for LSW,
  LSW_ForcedRupture, AND rate-state; no per-law / per-direction / threshold gate.
- **Regression contract (relaxed, see B.3):** TPV102/104/205 + BP5
  physically-exact (`worst_rel ≤ 1e-13` vs pre-fix) AND cross-rank bit-identical
  (R-101 `max_rel_diff == 0`). Overrides CLAUDE.md byte-exact; recorded; needs
  user sign-off.
- **MPI-collective-safe** (R-1600 deadlock class): every rank reaches the
  exchange; ranks with no shared fault faces contribute empty buffers and do not
  hang.
- **Boss = lower MPI rank** (unique, geometry-independent; `MFEM_VERIFY` exactly
  one boss per shared point).
- **Files Requiring Extreme Care:** `wave_operator.inl/.hpp`, `fault_face_flux.cpp`.
- Do not touch the friction solvers; do not add hardcoded constants.

## Phase 0: pin the seed — magnitude + source (DIAGNOSTIC; gates the fix choice)

### Goal
Close the open part of the root cause (B.1.5): measure the cross-rank input
difference at the diverging QP at FULL precision, and localise its source to
interpolation (1) / ghost-exchange bug (2) / bulk divergence (3) — so we know
whether the reconcile (Phase 2) is the right fix or a cheaper source-fix exists.

### Files to Modify
- `dynamic/wave_operator.inl` — extend the existing `SEAS_DIAG_XRANK` trace
  (`ComputeADERSharedFaceFluxRHS`, already committed `a9bd4d2`/`a28a5c5`):
  1. Print `I_self_can` / `I_nbr_can` at **`%.17e`** (full double precision), not
     `%.10e`, so the sub-11-digit difference is visible.
  2. Add a **raw-DOF** comparison: alongside the interpolated `I_*_can`, also dump
     the owning element's stored DOFs that feed the QP (before `shape×DOF`
     interpolation) for both the local element and the ghost copy. This separates
     "DOFs identical but interpolation differs" (source 1) from "ghost DOFs differ
     from the owner's" (source 2/3).

### Detailed Requirements
1. At the target QP, print per rank, per sub-step: `I_self_can[SXX/SXY/SXZ]` and
   `I_nbr_can[SXX/SXY/SXZ]` at `%.17e`.
2. Print the raw face-neighbour DOF buffer entries the QP interpolates from
   (`nbr_data[c][nbr_idx*ndof_per_el_ + i]`) and the local DOFs
   (`Q_data[c*ndof_total_ + dof_offset1 + i]`) at `%.17e`, for one matched
   (local-on-A == ghost-on-B) element, so the two ranks' copies can be diffed.
3. Run `spatial_dyn_resDc2_XRANKdiag_8N_400r_dev_2hr_safs.sbatch` (the OBLIQUE
   Dc2 sim, trace on, auto-extract per rank). It is now wired **non-fatal**
   (`SEAS_R101_NONFATAL=1` + `SEAS_DIAG_BLOWUP=1`) so it does double duty:
   (a) the `%.17e` + `RAWDOF` trace localises the seed (this phase); and
   (b) dropping the hard gate lets the dip-loaded run proceed PAST the
   divergence, so we also see whether it stays bounded like zerodip or blows up
   (the B.2 open test) and it **reproduces the speckle as the pre-fix baseline**
   for the Phase-4 final check.

### Decision (recorded in the diagnosis doc, drives Phase 2 vs a source-fix)
- **Raw ghost DOFs == owner DOFs (bit-exact) but `I_*_can` differ** → source (1)
  interpolation → inherent → proceed to Phase 2 (reconcile).
- **Raw ghost DOFs differ from owner DOFs** → source (2) ghost-exchange bug →
  STOP and write a separate plan to fix the exchange (cheaper; would let R-701's
  no-broadcast design stand). Do NOT build the reconcile yet.
- **DOFs identical at step 0 but drift apart over steps** → source (3) bulk
  divergence → inherent → proceed to Phase 2.

### Acceptance Criteria
- [ ] The `%.17e` trace shows the exact magnitude of the cross-rank input gap at
      the onset QP.
- [ ] The raw-DOF comparison localises the source to (1), (2), or (3).
- [ ] The diagnosis doc records the source and the resulting fix decision.

### Dependencies
- Depends on: nothing (extends the shipped diagnostic). Required by: the choice
  of Phase 2 (reconcile) vs a source-fix.

## Phase 1: reproducible local test (injection-based, np=2, both laws)

### Goal
A deterministic np=2 unit test that reproduces the cross-rank split by injecting
a 1-ULP input difference at the slip-onset threshold and asserts the two ranks'
`DOFData` are bit-identical — RED before Phase 2, GREEN after — for BOTH
rate-state and LSW (proving the fix is method-invariant).

### Files to Create
- `tests/unit/test_shared_fault_reconcile_cross_rank.cpp`.

### Files to Modify
- `Makefile` — 5 entries mirroring the tilted-test wiring (SRC ~323, OBJ ~623,
  link target `seas_test_shared_fault_reconcile_cross_rank` with the tilted
  test's object set, compile rule ~2719, run target
  `test-shared-fault-reconcile-cross-rank` at np=2). Not in `make test`.
- `dynamic/fault_face_flux.hpp` — behind the existing `#ifdef SEAS_TEST_INTERNAL`:
  `static real_t s_seas_test_qplus_xz_perturb_ulp = 0.0;` consumed in
  `ComputeTrialTraction` to add `s · ULP · Q_plus[SXZ]` ONLY when set (compiles
  out in production). (Or, if it reproduces, inject by perturbing one rank's
  bulk Q DOF in the test — preferred; decide at implementation.)

### Detailed Requirements
1. Fixture: reuse `BuildTwoTetFaultMesh()` (planar y=0, 2 tets) + the
   `cy<0→rank0/cy≥0→rank1` partition from
   `test_rupture_multistep_serial_vs_parallel.cpp`. One shared fault face, 3 QPs.
2. Drive `|τ|` across `μ_s·σ_n` during the run (derive load from `TPV205Params` /
   `TPV102Params`); assert at least one step has `0 < V_abs < V_small` on one
   rank (the onset is reached), else FAIL "did not cross slip-onset threshold."
3. On rank 0 only, perturb its `Q_plus[SXZ]` by exactly 1 ULP (`std::nextafter`)
   vs rank 1's matching input — the deterministic stand-in for the 1e-14
   interpolation seed.
4. Run `kNSteps` ADER-2 steps through the threshold.
5. Each step gather both ranks' shared-QP `DOFData` via
   `VerifySharedFaultDOFDataConsistency(tol,&wr,&wf,/*abort=*/false)` (added in
   `a9bd4d2`); assert `wr == 0`.
6. Loop over `{RateState, LSW}` in one executable; the assertion is per-law.

### Acceptance Criteria
- [ ] Builds + runs at np=2.
- [ ] Before Phase 2: RED — LSW `wr→O(1)`; rate-state `wr~1e-14`>0 (shows the
      defect is method-invariant, just sub-threshold for rate-state).
- [ ] After Phase 2: GREEN — `wr == 0` for BOTH laws, every step.
- [ ] Production build (`SEAS_TEST_INTERNAL` off): hook compiles out.

### Dependencies
- Depends on: nothing (built regardless of Phase 0). Required by: Phase 2
  (oracle), Phase 4. Note: if Phase 0 finds a ghost-exchange bug (source 2), the
  oracle's GREEN target becomes the source-fix, not the reconcile — the test
  (cross-rank bit-identity at the threshold) is the right oracle either way.

## Phase 2: the reconcile (the fix) — CONDITIONAL on Phase 0 = inherent seed (source 1 or 3)

> If Phase 0 finds source (2) (ghost-exchange not bit-exact), skip this phase and
> fix the exchange instead (separate plan); the reconcile is for the inherent-seed
> case where both ranks' inputs cannot be made bit-identical.

### Goal
After the friction solve on a shared fault face, both ranks hold bit-identical
`DOFData` and assemble from the boss's imposed stress — for ALL friction laws —
so the slip/lock decision is made once and copied.

### Files to Modify
- `dynamic/wave_operator.inl` `ComputeADERSharedFaceFluxRHS` (~4671-4797): make
  the shared-fault QP handling **two-pass** with a reconcile between them.
- `dynamic/wave_operator.hpp`: declare the helper + assembly-buffer struct.
- Factor the R-101 record gather/pair (`:5470+`) into a shared helper used by
  both the verify and the reconcile (one matcher).

### Detailed Requirements
1. Helper (parallel-only):
   ```cpp
   int WaveOperator<MeshType>::ExchangeAndPairSharedFaultQPs(
        int npay, const std::vector<double>& local_payload,
        const std::function<void(int local_qp, const double* peer_payload,
                                 int peer_rank)>& cb) const;
   ```
   Gathers one record/local-shared-QP (face-vertex-key + qp_idx + rank + payload)
   via the existing `MPI_Allgatherv`, pairs by (face-key, qp_idx), invokes `cb`
   for each paired local QP. The verify uses `npay=8` + a compare callback; the
   reconcile uses `npay = 8 + 2·NUM_STATE` (8 DOFData fields + `I_imp_plus_can` +
   `I_imp_minus_can`) + an overwrite callback.
2. Per-QP assembly buffer:
   ```cpp
   struct SharedFaultQPAssembly {
      int dof_offset1, ndof, dof_idx; bool elem1_on_plus; real_t w;
      real_t can_n[3], can_t1[3], can_t2[3];
      real_t shape1[MAX_NDOF];
      real_t I_imp_plus_can[NUM_STATE], I_imp_minus_can[NUM_STATE];
   };
   ```
3. Pass 1 — compute (no assembly): run the existing friction dispatch (LSW /
   LSW_ForcedRupture / rate-state) → fill `fault_dof_data_[dof_idx]` and
   `I_imp_*`; store the buffer; do NOT touch `rhs`.
4. Reconcile: payload = 8 DOFData fields + `I_imp_plus_can` + `I_imp_minus_can`.
   Call the helper; in the callback, if `peer_rank < my_rank_` (peer is boss),
   overwrite this rank's `fault_dof_data_[dof_idx]` AND the buffer's
   `I_imp_*_can` with the peer's payload. `MFEM_VERIFY` every local shared QP was
   paired.
5. Pass 2 — assemble from the (possibly overwritten) `I_imp_*_can`, rotated to
   global via `T_can` (rebuilt from buffered `can_*`, bit-identical on both
   ranks), exactly as the current inline assembly.
6. No friction-law gate (method-invariant). Rate-state `else` branch included.
7. Collective safety: the helper uses the same `MPI_Allreduce(any_shared)`
   short-circuit + `MPI_Allgatherv` as the verify; empty contribution from ranks
   with no shared fault QPs.

### Edge Cases to Handle
- Unpaired shared QP → `MFEM_VERIFY` fail (same policy as the verify).
- Same-`elem1_on_plus`-on-both-ranks edge (θ≈90°) → boss-by-rank still unique.
- Serial / np=1 → no shared faces → two-pass degenerates to the current path
  (byte-exact).

### Acceptance Criteria
- [ ] Phase-1 test GREEN for both laws (`wr == 0` every step).
- [ ] TPV102/104/205 + BP5: `worst_rel ≤ 1e-13` vs pre-fix AND R-101
      `max_rel_diff == 0`.
- [ ] Local fault/TPV suite green (fault-basis trio, shared-fault role/dof-data,
      interior-flux ×3, godunov identity, tpv102 locked/absorbing/pepper/
      ader-smoke, multistep + tilted serial-vs-parallel, sign-flipped truth
      table).
- [ ] np=2,4,8 incl. a partition with some rank lacking shared fault faces — no
      deadlock.

### Dependencies
- Depends on: Phase 1. Required by: Phase 4.

## Phase 3: correctness guard (runtime, method-invariant)

### Goal
Keep `VerifySharedFaultDOFDataConsistency` as the permanent runtime sentinel so
any future regression aborts loudly; post-fix it must report `max_rel_diff == 0`.

### Files to Modify
- `dynamic/wave_operator.inl`: refactor the verify onto the shared
  `ExchangeAndPairSharedFaultQPs` helper (one matcher); update the stale R-501
  wording in the abort message to point at the Phase-2 reconcile.
- Keep `SEAS_R101_NONFATAL` (commit `fc28454`) as the diagnostic escape hatch;
  keep the existing gated cadence; add `SEAS_VERIFY_XRANK_EVERY=N` (default
  unchanged) for debugging.

### Detailed Requirements
1. The guard stays method-invariant (checks all shared QPs regardless of law).
2. Document that post-fix the expected `max_rel_diff` is `0`; any nonzero is a
   real regression, not roundoff.
3. Negative test (Phase 1, with the reconcile disabled via a test switch): the
   guard trips — proving it still catches a real desync.

### Acceptance Criteria
- [ ] Post-fix Dc2 `XRANK` Frontera run: guard reports `max_rel_diff == 0`
      through and past t≈0.455 s; no abort; run advances past t=1.0 s.
- [ ] Negative test trips the guard.

### Dependencies
- Depends on: Phase 2. Required by: Phase 4.

## Phase 4: full regression + Frontera re-validation

### Goal
End-to-end: the Dc2 mesh advances past t=1.0 s (no abort, no blow-up); fix is
config-agnostic.

### Detailed Requirements
1. Local regression battery (Phase-2 AC) green; document unchanged pre-existing
   failures (adjacent-triangle pepper-bug, r101 missing-precondition, macOS MUMPS
   Bus error).
2. Re-run the `SEAS_DIAG_XRANK` Dc2 sbatch: at the onset QP both ranks now show
   identical V1/V2/τ*_corr through t=0.478 s; R-101 passes; reaches tfinal.
3. Re-run the physical `D_c=1.0` baseline config (not just the inflated Dc2) to
   confirm the fix is not tuned to Dc2.
4. **zerodip as an issue-B oracle (job 7747119 baseline = `worst_rel → 1.0` +
   visible speckle):** re-run zerodip post-fix; `worst_rel` must stay ≤1e-13 (no
   `DIVERGED`) AND the speckled slip-rate elements outside the rupture must
   vanish. This is the cleanest already-set-up confirmation that the reconcile
   fixes issue B in the strike-only (non-catastrophic) case — separate from the
   blow-up.
5. **oblique non-fatal as the dip-channel speckle oracle
   (`spatial_dyn_resDc2_XRANKdiag_*.sbatch`, now non-fatal):** the pre-fix run
   reproduces the speckle WITH dip traction present (and records bounded-vs-blow-
   up — the B.2 open test). Post-fix rerun: `worst_rel ≤ 1e-13` AND the speckle
   clears. If the pre-fix oblique run also stays bounded (no blow-up after Part
   A), this run alone is the full issue-B final check. If the speckle survives
   at `worst_rel≈0`, a co-present cause (under-resolution / SSO) is implicated,
   NOT the reconcile (decide next steps then).

### Acceptance Criteria
- [ ] Local regression green (modulo documented pre-existing).
- [ ] Frontera Dc2 AND `D_c=1.0`: no R-101 abort, no blow-up; `V_max` peaks then
      decreases.
- [ ] **zerodip post-fix:** `worst_rel ≤ 1e-13` for the whole run (was → 1.0 at
      t≈0.595 s) AND no speckle in ParaView (was scattered cyan elements). `V_max`
      unchanged (~6 m/s, still physical).
- [ ] **oblique non-fatal post-fix:** `worst_rel ≤ 1e-13` AND the speckle clears
      in ParaView (pre-fix baseline = this run with the gate non-fatal). The
      bounded-vs-blow-up outcome of the pre-fix oblique run is recorded either way
      (B.2 test); a surviving speckle at `worst_rel≈0` flags a co-present cause.

### Dependencies
- Depends on: Phase 2, Phase 3.

## Testing Strategy

### Q2 — can we reproduce the issue locally as a guard?  Three layers:
1. **Deterministic local guard (primary):** the Phase-1 injection test (np=2). The
   clean 2-tet fixture is bit-identical across ranks, so it CANNOT make the seed
   naturally — we inject the 1-ULP cross-rank difference the real run gets from
   interpolation/bulk drift, drive to the slip-onset threshold, and assert
   cross-rank bit-identity. RED before Phase 2, GREEN after, for BOTH laws. It is
   cheap (np=2, a few steps) and **MUST be added to the parallel CI target**
   (`test-parallel`-class) so a future regression of the reconcile aborts CI —
   this is the standing local guard. (Mark it in the Makefile as a parallel test
   and wire it into the `test-parallel` aggregate, NOT only as a standalone
   target like the tilted test.)
2. **Natural local reproduction (optional, heavier):** a multi-element np≥2 LSW
   fixture with a slow nucleation front crossing a shared QP at the threshold (a
   "mini-SAFS"). It would reproduce the divergence WITHOUT injection, but the seed
   magnitude is mesh-dependent and it is far heavier to stand up, so it is a
   nice-to-have, not the guard. The injection test is the reliable guard.
3. **Production runtime guard:** the Phase-3 R-101 verify (already shipped,
   `SEAS_R101_NONFATAL` for diagnostics) — post-fix `max_rel_diff==0`; a negative
   test (reconcile disabled) trips it. Plus the **zerodip job as a full-physics
   guard** (post-fix `worst_rel ≤ 1e-13` + no speckle).

### Other
- **Local oracle:** Phase-1 injection test (np=2, both laws), RED→GREEN — proves
  the reconcile absorbs the seed for any friction law.
- **Guard:** Phase-3 verify, post-fix `max_rel_diff==0`; negative test trips it.
- **Regression:** TPV102/104/205 + BP5 — capture pre-fix, assert `≤1e-13` AND
  R-101 `max_rel_diff==0`; stash-rebuild-baseline to confirm the ~1e-14 TPV change
  is the reconcile (non-boss→boss), not a physics change.
- **Part-A non-regression:** the zerodip job stays clean.
- **MPI robustness:** np=2/4/8, incl. a no-shared-fault-on-some-ranks partition.

## Risk Assessment
- **Deadlock (R-1600 class):** the new exchange must be reached by every rank;
  detect via np=4/8 smoke hangs.
- **Per-sub-step exchange cost:** correctness-first; if hot on 400r, move the
  helper to point-to-point over the face-nbr topology (same payload). NO
  kink-proximity gate (re-introduces problem-specificity).
- **Two-pass refactor indexing slip** (dof_offset1/shape1/w): the interior path
  is untouched; tilted + multistep + adjacent-triangle tests cover shared-face
  assembly; assert Pass-2 reproduces the pre-refactor RHS when the reconcile is a
  no-op.
- **Byte-exact relaxation:** the ~1e-14 TPV change overrides a CLAUDE.md
  non-negotiable; recorded; needs sign-off.
- **Tricky existing code:** `ComputeADERSharedFaceFluxRHS` (4671-4797, the
  R-1601 inline fallback + assembly), `VerifySharedFaultDOFDataConsistency`
  (5364+, the matcher being factored), the canonical-frame reconstruction
  (4695-4709, must stay bit-identical across ranks for the Pass-2 `T_can`
  rotation).
