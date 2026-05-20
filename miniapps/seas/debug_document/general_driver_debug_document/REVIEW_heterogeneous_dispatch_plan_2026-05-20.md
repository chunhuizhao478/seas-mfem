# Code Review: PLAN_heterogeneous_volume_bc_fault_dispatch_2026-05-20.md (2026-05-20)

## Review Scope
- Plan reviewed (the "code" under audit is this PLAN):
  `miniapps/seas/debug_document/general_driver_debug_document/PLAN_heterogeneous_volume_bc_fault_dispatch_2026-05-20.md`
- Codebase files checked against the plan's claims:
  `dynamic/wave_operator.inl`, `dynamic/wave_operator.hpp`,
  `dynamic/godunov_flux.hpp`, `dynamic/godunov_flux_pool.hpp`,
  `dynamic/godunov_flux_pool.cpp`, `dynamic/spatial_setup.hpp`,
  `drivers/spatial_dyn_driver.cpp`,
  `tests/unit/test_phaseh_wave_operator_constant_parity.cpp`
- Domain context: `miniapps/seas/CLAUDE.md` (canonical frame,
  fluctuation-Q, byte-parity regression contract), Phase R bimaterial
  notes.
- Findings are factual/structural errors in the plan that would mislead
  the implementing agent. Suggested fixes are edits to the plan markdown.

## Findings

### [R-001] CRITICAL [PLAN §4 Phase 1 / §2 Group A] — ADER CK recursion is NOT a per-element loop; plan misdescribes the structure and understates the work

**Category:** DEVIATION (plan-vs-code structural error)

**Description:**
The plan (§4 Phase 1) states *"All three are per-element loops already,
so threading per-element Jacobians is structurally local."* True only for
**A1** (`ComputeVolumeRHS`, `wave_operator.inl:1588` — a real
`for (int e ...)` loop). FALSE for **A2/A3**:
- `ComputeADERTimeIntegrated` (`:1816`) and `ComputeADERSubStepStates`
  (`:1913`) operate on **global** vectors and apply ONE global Jacobian
  `A_d = flux_.GetReferenceStarMatrix(d)` to the ENTIRE field via the free
  function `ApplyJacobianPerDOF(A_d, dQ_dxd, D_next, ndof_total_, -1.0)`
  (`:1868`, `:1987`).
- `ApplyJacobianPerDOF` (`:1779-1810`) loops
  `for (int i = 0; i < ndof_total; i++)` applying the SAME 9×9 `A` to
  every DOF. No element index is in scope; the CK recursion is not
  element-blocked.

Making A2/A3 heterogeneous requires **restructuring**
`ApplyJacobianPerDOF` (or its callers) to apply each element's star
matrix to that element's DOF block — NOT a "structurally local" `At(e)`
fetch. TPV31 runs ADER (`ader_order = 2` in `tpv31.toml`; `--ader-order 3`
in the p2/O3 sbatch), so the CK predictor is on the critical path.

**Trigger:** Any heterogeneous (`Mode::Coefficient`) ADER run — i.e.
every TPV31 run.

**Actual behavior (plan followed literally):** Implementer looks for a
per-element loop in A2/A3, finds none; the single `ApplyJacobianPerDOF`
keeps applying one Jacobian to all DOFs → the ADER predictor uses one
material everywhere → wrong bulk wave speeds in the 1-D velocity
structure, silently.

**Expected behavior:** Plan specifies a per-element CK application: a
per-element variant of `ApplyJacobianPerDOF` (loop elements; for element
`e`, apply `FluxForElem_(e).GetReferenceStarMatrix(d)` to the
`ndof_per_el_` DOFs of `e`, component stride `ndof_total_`), wired into
both `ComputeADERTimeIntegrated` and `ComputeADERSubStepStates`. A1 stays
the simple in-loop fetch.

**Suggested fix (edit Group A table + Phase 1 in the plan):**
```diff
-| A1 | `wave_operator.inl:1627-1629` (`ComputeVolumeRHS`, RK4) | scalar members `Ax_/Ay_/Az_` ... |
-| A2 | `wave_operator.inl:1867` (`ComputeADERTimeIntegrated`) | `flux_.GetReferenceStarMatrix(d)` |
-| A3 | `wave_operator.inl:1986` (`ComputeADERSubStepStates`) | `flux_.GetReferenceStarMatrix(d)` |
-
-All three are **per-element loops** already (`for e`), so threading
-per-element Jacobians is structurally local.
+| A1 | `wave_operator.inl:1627-1629` (`ComputeVolumeRHS`, RK4) | scalar `Ax_/Ay_/Az_`. **Genuine per-element `for e` loop** — in-loop `FluxForElem_(e).GetReferenceStarMatrix(0/1/2)`. |
+| A2 | `wave_operator.inl:1868` (`ComputeADERTimeIntegrated`) | one global `A_d` applied to ALL DOFs via `ApplyJacobianPerDOF` (`:1779`). **NOT a per-element loop.** |
+| A3 | `wave_operator.inl:1987` (`ComputeADERSubStepStates`) | same global `ApplyJacobianPerDOF(A_d, ...)`. **NOT a per-element loop.** |
+
+A1 is a genuine per-element loop (simple in-loop `At(e)`).  A2/A3 are
+NOT: `ApplyJacobianPerDOF` (`:1779-1810`) applies a single 9×9 `A`
+uniformly across `ndof_total_`.  Heterogeneity requires a per-element
+variant — loop elements, apply `FluxForElem_(e).GetReferenceStarMatrix(d)`
+to element `e`'s `ndof_per_el_` DOFs.  This is the largest piece of
+Phase 1; split it as Phase 1a (A1, RK4 volume) and Phase 1b (A2/A3, ADER
+CK via the new per-element `ApplyJacobianPerDOF`), each gated by the
+Mode::Constant parity test.
```

**Test case (verification):**
```python
def test_R001_ader_ck_is_global_not_per_element():
    body = read("dynamic/wave_operator.inl", 1855, 1872)
    assert "ApplyJacobianPerDOF(A_d" in body
    assert "for (int e" not in body            # no per-element loop here
    apply = read("dynamic/wave_operator.inl", 1779, 1810)
    assert "for (int i = 0; i < ndof_total" in apply  # global over all DOFs
```

---

### [R-002] MODERATE [PLAN §3 invariant #2 / §8 Risk] — pool builds `GodunovFlux` from 6-sig-fig-ROUNDED material; "byte-identical" invariant and Risk-§8 mitigation are false as stated

**Category:** ASSUMPTION (false premise underpinning the regression contract)

**Description:**
Invariant #2 (§3) claims that for `Mode::Constant` the pool builds its
flux "from those constants, so `At(e)` ≡ `flux_` numerically"; Risk §8
says the pool "reproduces `flux_` exactly." But
`GodunovFluxPool::Build` (`godunov_flux_pool.cpp:100-106`) builds from
**rounded** values:
```cpp
const real_t rl = round_sig(lam, dedup_sig_figs);   // 6 sig figs
const real_t rm = round_sig(mu,  dedup_sig_figs);
const real_t rr = round_sig(rho, dedup_sig_figs);
unique_fluxes_.emplace_back(std::make_unique<GodunovFlux>(rl, rm, rr));
```
So `At(e)` uses 6-sig-fig-rounded material; the scalar `flux_` uses exact
material. They match only when the constant is exactly 6-sig-fig
representable. The existing parity test passes only because it uses clean
constants (`k_lambda = k_mu = 32.0e9`, `k_rho = 2670.0` —
`test_phaseh_wave_operator_constant_parity.cpp:88-90`).

Today the heterogeneous ctor's volume term still uses the exact `Ax_`, so
it is byte-identical regardless of sig-figs. Plan Phase 1 (A1) switches it
to `At(e)` (rounded) — a NEW discrepancy vs the scalar ctor for any
>6-sig-fig material. The §6 proposed parity-test extension would FAIL if
the implementer picks a realistic (TPV31 layer-derived) constant.

**Trigger:** Mode::Constant MaterialField ctor with `(λ,μ,ρ)` exceeding 6
sig figs, after Phase 1 routes the volume term through `At(e)`.

**Actual behavior:** `At(e)` ≠ `flux_` at ~1e-6 relative; a bit-equality
parity assertion fails; wasted debugging of a non-bug.

**Expected behavior:** Plan states the truth and picks a resolution
(build pool flux from EXACT values, rounding only the dedup key).

**Suggested fix (edit invariant #2 and Risk §8):**
```diff
-2. **`Mode::Constant` via the MaterialField ctor stays byte-identical.**
-   All elements share one `(λ,μ,ρ)`; the pool dedups to a single
-   `GodunovFlux` built from those constants, so `At(e)` ≡ `flux_`
-   numerically.
+2. **`Mode::Constant` stays byte-identical ONLY IF the pool's cached flux
+   is built from EXACT material.**  `GodunovFluxPool::Build` currently
+   builds from 6-sig-fig-ROUNDED `(λ,μ,ρ)` (`godunov_flux_pool.cpp:103-106`),
+   so `At(e)` differs from exact `flux_` at ~1e-6 for >6-sig-fig
+   constants.  Today this is masked because the volume term uses exact
+   `Ax_`; Phase 1 removes the mask.  **Resolution (Phase 1, before A1):**
+   change `godunov_flux_pool.cpp:106` to build from exact `lam,mu,rho`
+   (keep `make_key` rounded for the dedup KEY only); re-run the
+   interior-path parity test.  Then invariant #2 holds for any constant.
```
```diff
-- **Byte-parity drift on Mode::Constant.** ... Mitigation: the existing
-  interior-path parity test already passes, implying the pool reproduces
-  `flux_` exactly for Mode::Constant; ...
+- **Byte-parity drift on Mode::Constant — CONFIRMED.**
+  `GodunovFluxPool::Build` builds the cached flux from rounded `(λ,μ,ρ)`
+  (`godunov_flux_pool.cpp:103-106`); `At(e)` ≠ exact `flux_` for
+  >6-sig-fig constants.  The existing parity test passes only because it
+  uses clean constants (`32.0e9 / 2670`).  Mandatory mitigation: build
+  the cached flux from EXACT values (round the dedup key only), and have
+  the extended parity test use BOTH a clean and a >6-sig-fig constant so
+  the gate actually catches regressions.
```

**Test case:**
```python
def test_R002_pool_rounds_cached_flux():
    src = read("dynamic/godunov_flux_pool.cpp")
    assert "round_sig(lam" in src and "make_unique<GodunovFlux>(rl, rm, rr)" in src
    # => At(e) uses rounded material; invariant #2 false for >6 sig figs.
```

---

### [R-003] MODERATE [PLAN §2 Group C, row C4] — lines 4921/4924 are in `ComputeADERFaceFluxRHS` (ADER **local**), not `ComputeADERSharedFaceFluxRHS`

**Category:** DEVIATION (incorrect function attribution → a site can be missed)

**Description:**
Group-C row C4 assigns "per-QP variant `:4921-4924`" to
`ComputeADERSharedFaceFluxRHS`. But `ComputeADERFaceFluxRHS` spans
`:4140-5183` and `ComputeADERSharedFaceFluxRHS` starts at `:5183`. Lines
4921/4924 are inside `ComputeADERFaceFluxRHS` (ADER **local**) — a SECOND
ADER-local fault dispatch block alongside `:4603-4606`. The true
ADER-shared fault site is `:5568` only.

**Trigger:** Implementing Phase 3 from the plan's table.

**Actual behavior:** Implementer edits 4603-4606 and 5568, misses
4921-4924 (it is in a function the plan told them not to look in) → one
ADER-local fault path keeps scalar `flux_.Interior` → wrong fault
bulk-side flux on the per-QP-batched path.

**Expected behavior:** 4921/4924 listed under ADER local (C3); note
ADER-local has two fault blocks.

**Suggested fix (edit C3/C4 rows):**
```diff
-| C3 | `wave_operator.inl:4603-4606` (`ComputeADERFaceFluxRHS`) | ADER local fault |
-| C4 | `wave_operator.inl:5568` and per-QP variant `:4921-4924` (`ComputeADERSharedFaceFluxRHS`) | ADER shared fault |
+| C3 | `wave_operator.inl:4603-4606` **and** per-QP-batched variant `:4921-4924` (both in `ComputeADERFaceFluxRHS`) | ADER local fault — TWO blocks |
+| C4 | `wave_operator.inl:5568` (`ComputeADERSharedFaceFluxRHS`) | ADER shared fault |
```

**Test case:**
```python
def test_R003_4921_is_in_ader_local():
    starts = grep_line("dynamic/wave_operator.inl", "::ComputeADERFaceFluxRHS")       # 4140
    shared = grep_line("dynamic/wave_operator.inl", "::ComputeADERSharedFaceFluxRHS") # 5183
    assert starts < 4921 < shared
```

---

### [R-004] LOW [PLAN §3 / §4 Phase 1] — references nonexistent `GetAy()/GetAz()` accessors

**Category:** BUG (API that does not exist)

**Description:**
§3 writes "`FluxForElem_(e).GetAx()/GetAy()/GetAz()`-equivalent for A1".
`GodunovFlux` exposes only `GetAx()` (`godunov_flux.hpp:176`),
`GetAxPlus()`, `GetAxMinus()` — no `GetAy()/GetAz()`. The correct
per-direction accessor is `GetReferenceStarMatrix(0/1/2)`
(`godunov_flux.hpp:216`).

**Suggested fix (edit §3):**
```diff
-`FluxForElem_(e).GetReferenceStarMatrix(d)` (A2/A3) and
-`FluxForElem_(e).GetAx()/GetAy()/GetAz()`-equivalent for A1 (note: A1
-currently reads the cached members `Ax_/Ay_/Az_`; switch to the
-element's `BuildJacobian`/`GetReferenceStarMatrix(d)` ...).
+`FluxForElem_(e).GetReferenceStarMatrix(d)` for d=0,1,2 — the ONLY
+per-direction accessor (`GodunovFlux` has `GetAx()` but NOT
+`GetAy()/GetAz()`).  A1 currently reads cached members `Ax_/Ay_/Az_`
+(== `GetReferenceStarMatrix(0/1/2)` numerically); switch it to the
+element's `GetReferenceStarMatrix(d)`.
```

**Test case:**
```python
def test_R004_no_getay_getaz():
    hdr = read("dynamic/godunov_flux.hpp")
    assert "GetAx()" in hdr
    assert "GetAy()" not in hdr and "GetAz()" not in hdr
    assert "GetReferenceStarMatrix(int dir)" in hdr
```

---

### [R-005] LOW [PLAN §2 inventory completeness] — `ComputeMaxDt` omitted; it is already heterogeneous and the scalar fallback (`:6004`) must NOT be "fixed"

**Category:** QUALITY (incomplete inventory; preempt a wrong "fix")

**Description:**
The plan claims a "full inventory" but omits `ComputeMaxDt`
(`wave_operator.inl:5921`), which has `flux_.GetCp()` at `:6004`. That is
the scalar-ctor fallback; the heterogeneous path
(`:5974` `if (!per_elem_h_.empty())` → `:5985-6001`) already computes the
stable dt per-element from `per_elem_lmr_`/`per_elem_h_`. A blanket
`flux_ → At(e)` sweep might wrongly "fix" `:6004` and break the
scalar-ctor dt.

**Suggested fix (add a note after the Group A/B/C tables in §2):**
```diff
+**Already heterogeneous — do NOT change:** `ComputeMaxDt`
+(`wave_operator.inl:5921`) branches on `!per_elem_h_.empty()` (`:5974`)
+and computes the stable timestep per-element (`:5985-6001`).  The
+`flux_.GetCp()` at `:6004` is the scalar-ctor fallback and MUST stay.
```

**Test case:**
```python
def test_R005_computemaxdt_already_per_element():
    body = read("dynamic/wave_operator.inl", 5921, 6005)
    assert "if (!per_elem_h_.empty())" in body
    assert "per_elem_lmr_[e]" in body and "per_elem_h_[e]" in body
    assert "flux_.GetCp()" in body  # scalar fallback retained
```

---

## Summary
- Critical issues: 1 (R-001)
- Moderate issues: 2 (R-002, R-003)
- Low issues: 2 (R-004, R-005)
- Plan accuracy: **PARTIAL** — the dispatch-site inventory and the
  per-element data model are largely correct, but the ADER-CK structure
  (R-001), the byte-parity premise (R-002), and one site attribution
  (R-003) are wrong and would mislead the implementer.
- Verdict: **PASS WITH FIXES** — apply R-001..R-005 to the plan before
  implementation. R-001 and R-002 are load-bearing: without them the
  implementation would silently use one material in the ADER predictor
  (R-001) and/or fail its own parity gate (R-002).

## Unreviewed Areas
- Bi-material *across-fault* exclusion (plan §7) accepted as correct for
  TPV31 (vertical fault, depth-only material → `Zp_plus == Zp_minus`),
  verified via `spatial_setup.hpp:79-82,125-128`. Not re-derived for
  non-TPV31 geometries (out of scope per plan).
- Performance of the new per-element ADER-CK application (R-001) not
  benchmarked — follow-up, not a correctness finding.
- Absorbing/free-surface background checked and CORRECT (`Q_bg = 0`,
  fluctuation-Q, `spatial_dyn_driver.cpp:1574`) — NOT a finding.
