# Code review — fault-flux / over-integration / resample / free-surface (dip-slip context)

Date: 2026-06-06
Reviewer: Claude (code-review, adversarial)
Trigger: `/code-review` on the TPV31 residual dip-slip path.
Scope reviewed (read in full or in the cited ranges):
  * `dynamic/fault_face_flux.cpp` — trial traction, dip/strike/normal channels
  * `dynamic/tpv205_friction.hpp` — LSW strength + slip-parallel decomposition
  * `drivers/spatial_dyn_driver.cpp` — over-integration + resample wiring, ADER substep path
  * `dynamic/wave_operator.inl` / `.hpp` — `SetFaultOverint` / `RebuildFaultQuadrature_` / fault residual
  * `dynamic/friction_substep_iterator.cpp` — RS (Δψ) + LSW (magnitude) resample
  * `dynamic/godunov_flux.cpp` — `AbsorbingTotal` / `FreeSurfaceTotal` / `BuildFrame`
  * `tests/unit/test_fault_planar_serial.cpp` — Phase-0 dealiasing harness

Reference for "correct": SCEC TPV31 spec + the project's exact-Riemann fault-flux
derivation (Eq. 7a–c) + the fault-dealiasing plan. Companion diagnosis:
`tpv31_dip_slip_residual_findings_2026-06-06.md`.

Verdict legend: **PASS** (correct), **GAP** (missing capability / limitation),
**FOOTGUN** (correct but easy to misuse), **TEST-GAP** (claimed coverage not real).

---

## R-1 — dip-channel trial traction + LSW decomposition — **PASS (no bug)**

`FaultFaceFlux::ComputeTrialTraction` (`fault_face_flux.cpp:49–85`) builds the
three fault channels by the identical exact-Riemann formula from the across-fault
velocity jump (local frame t1=dip, t2=strike):

```
sigma_n_trial = eta_p*([[v_n]]   + sigma_n^+/Zp^+ + sigma_n^-/Zp^-)   (7a)
tau1_trial    = eta_s*([[v_dip]] + tau1^+/Zs^+    + tau1^-/Zs^-)       (7b)  DIP
tau2_trial    = eta_s*([[v_str]] + tau2^+/Zs^+    + tau2^-/Zs^-)       (7c)
```

`tpv205_friction.hpp:188–205`: `V_abs=max(tau_abs-tau_strength,0)/eta_s`, then
`V1=V_abs*tau1_total/tau_abs`, `V2=V_abs*tau2_total/tau_abs` — slip rate exactly
parallel to the total shear traction; both branches guard `tau_abs>0`.

This is the standard formulation and is **correct**. The dip channel is fully
symmetric with the strike/normal channels — there is **no sign error, no frame
error, no asymmetric treatment** of the dip component. **The dip-slip excess is
NOT a coding bug in this path** (independently corroborated by the data: a
constant frame/sign bug would be depth-uniform and method-independent; the
observed artifact is surface-concentrated, depth-oscillating, and flips sign
between flux schemes — see the diagnosis doc §3).

Minor: the `SEAS_TEST_INTERNAL` cross-rank seed (`fault_face_flux.cpp:82`)
perturbs only `tau2_trial` (strike). There is no analogous dip-channel
(`tau1_trial`) test seed, so cross-rank shared-face consistency of the **dip**
channel has less deterministic test coverage. Low priority (the formula is shared
across channels), but a `tau1` twin would close the gap.

---

## R-2 — over-integration is forbidden on the mixed-flux path — **GAP (material to the cure)**

`WaveOperator::SetFaultOverint` (`wave_operator.inl:634–639`):

```
MFEM_VERIFY(k == 0 ||
            (mixed_flux_mode_ == MixedFluxMode::None && !use_precomputed_face_fluxes_),
            "fault over-integration (k>0) is not compatible with the mixed-flux
             or precomputed-face-flux paths (Phase 1 scope).");
```

So `--fault-overint K` (K>0) **aborts** on any mixed-flux run. Consequence for
this investigation: **the dealiasing cure (if it is one) cannot be applied to the
mixed-flux production path at all** — the very path that already gives the smaller
surface dip slip. The over-integrated fault rule grows `nbf_per_face_`, but the
mixed-flux / precomputed-face caches are built at `2*order_` and would desync.

Recommendation: either (a) extend the over-integrated rule to the precomputed-face
/ mixed-flux caches (rebuild them at `FaultFaceQuadDegree()`), or (b) surface this
limitation in the driver `--help` and the dealiasing plan so nobody plans a
mixed-flux + over-integration A/B that cannot run. (My first-pass §7 in the
diagnosis doc made exactly that mistake; corrected in §9.)

---

## R-3 — LSW `--fault-resample` is magnitude-only (direction preserved) — **FOOTGUN**

`friction_substep_iterator.cpp:284–293` (LSW branch): the resample projects the
accumulated **slip-magnitude path length** `Σ|V|·dt`, then rescales `(slip1,slip2)`
to the resampled magnitude **with direction preserved**. The directional split
comes from the un-resampled friction solve.

Implication: `--fault-resample` **cannot remove a directional dip artifact** for
LSW. A reader who sees "resample dealiases the fault" and reaches for
`--fault-resample` to fix the dip slip will get **no directional correction** —
only over-integration changes the friction-solve direction. This is correct as
designed (it mirrors SeisSol's LSW dealiasing of the slip magnitude), but it is a
real footgun for this specific symptom.

Recommendation: one line in the flag help / schema doc: "`--fault-resample`
dealiases the slip-weakening **magnitude** only; it does not alter slip
**direction**, so it cannot by itself reduce a dip/strike directional bias."

---

## R-4 — Phase-0 harness does NOT test the production dealiasing path — **TEST-GAP**

`test_fault_planar_serial.cpp` is documented (header, line 28) as "the A/B testbed
for Phases 1-3 (over-integration + resample)." Two real limitations make that
claim only partly true:

1. `--fault-resample` is **inert** in the harness (lines 31, 118, 687–688) — it
   prints a NOTE and does nothing. So the **resample half is never exercised.**
2. The harness drives friction by calling `wave.AdvanceADER(Q,dt,2,Q_new)`
   **directly** (line 582) — the self-contained ADER with wave-operator-INTERNAL
   friction. The **production** ADER path
   (`spatial_dyn_driver.cpp:AdvanceADERWithSubStep_Spatial`, line 384) instead
   runs friction through `iterator.Advance(...)` (line 476, where
   `SetFaultResample`/the resample live) and uses `wave.AdvanceADER` only as the
   elastodynamic corrector via a side-channel `I_imp`. **The harness therefore
   exercises a different friction path than production** and cannot validate the
   production iterator + resample at all.

Net: the harness validly tests (a) the elastodynamic + fault-flux **leak
mechanism** (path-independent — this is the solid result) and (b) the
wave-operator-internal over-integration. It does **not** test the production
over-int + resample dealiasing. Recommendation: correct the header to scope the
claim, and — if a true local A/B of the cure is wanted — add a variant that drives
friction through the substep iterator with `--fault-resample` live.

(The 2026-06-06 extension added the `[[v_dip]]`/`[[v_strike]]` probes and two
guards; those are correct and pass. See diagnosis doc §9.)

---

## R-5 — absorbing / free-surface fluxes — **PASS (one known secondary limitation)**

* `GodunovFlux::FreeSurfaceTotal` (`godunov_flux.cpp:~480`): background-aware
  γ-mirror enforcing `(sigma - sigma_bg)·n = 0` (zero **fluctuation** traction),
  so a uniform `Q=Q_bg` reduces to the interior identity — no spurious radiation
  from a pre-stressed free surface. **Correct**, and the right choice for TPV31's
  non-zero surface traction.
* `GodunovFlux::AbsorbingTotal` (`godunov_flux.cpp:459`): background-aware but a
  **basic upwind** absorber (`Interior(nor,Q_self,Q_bg)`), not a characteristic /
  PML absorber — it leaks for oblique P/S. With `use_pml=false` and the 2.6 % box
  margin (prior doc §4), this is a **secondary aggravator** of any pre-existing
  seed, not the dip seed itself (the dip slip onsets at the rupture front, not at
  the ~10.8 s reflection return). Leave as-is for the dip-slip question; revisit
  only if a larger box / PML is on the table for other reasons.
* `BuildFrame` (`godunov_flux.cpp:304`): the non-fault boundary-flux frame
  (`up=(0,0,1)` flip near vertical). Distinct from the fault `FaultBasis` (BP5
  convention); does not touch the fault dip decomposition. Fine.

---

## Summary

| ID | Area | Verdict | Action |
|----|------|---------|--------|
| R-1 | dip trial traction + LSW decomposition | PASS (no bug) | none (optional: `tau1` test seed) |
| R-2 | over-integration vs mixed-flux | GAP | extend to mixed/precomputed, or document |
| R-3 | LSW resample = magnitude-only | FOOTGUN | document "no direction change" |
| R-4 | Phase-0 harness coverage | TEST-GAP | scope the header; add iterator-path variant |
| R-5 | absorbing/free-surface flux | PASS (sec.) | none for dip slip |

**Bottom line:** no correctness bug found in the dip-slip path. The dip slip is a
numerical-scheme artifact (diagnosis doc), not a coding error. The actionable
code findings are about the *cure infrastructure*: over-integration is blocked on
mixed-flux (R-2), resample can't fix a directional bias for LSW (R-3), and the
local harness can't actually test the production dealiasing (R-4).
</content>
