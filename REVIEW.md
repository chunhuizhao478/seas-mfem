# Code Review (round 2): 2026-05-04 — `safs/smoke_test/PLAN_safs_test.md`

## Review Scope

- **Plan reviewed:** `miniapps/seas/safs/smoke_test/PLAN_safs_test.md` (and
  its PDF render `PLAN_safs_test.pdf`).
- **Reference spec for compliance:** `miniapps/seas/drivers/seas_driver.cpp`
  — the existing BP5 driver. The plan's stated goal (lines 5–7) is "BP5-
  derived ... that exercises the SAFS multi-fault `.msh` end-to-end" with
  the only deltas being the boundary mapping (BoundaryConfig) and the
  per-DOF parameter override (uniform-VW, no-nucleation). Every other
  composition step should mirror BP5.
- **Header surfaces verified:** `domain/{boundary_config,elasticity_operator,
  domain_operator}.hpp`, `fault/{fault_geometry,rate_state_fault}.hpp`,
  `solver/{seas_operator,time_stepper}.hpp`,
  `friction/{dieterich_ruina,state_evolution}.hpp`,
  `config/bp5_params.hpp`, `constitutive/linear_elastic.hpp`.
- **Domain context:** `miniapps/seas/CLAUDE.md`.

## Why R-001..R-003 are deviations even though "the goal is to mirror BP5"

The user's challenge: if the plan mirrors BP5, where do these bugs come
from? Answer: **the plan was rewritten from memory rather than
copy-pasted from `drivers/seas_driver.cpp`, and three of its key code
snippets drifted from the actual BP5 pattern**. The plan IS supposed to
mirror BP5 and IS NOT supposed to have these. Specifically:

| ID | Plan does | BP5 reference does | Drift type |
|---|---|---|---|
| R-001 | `Vector state;` (size 0) | `Vector state(fault_op.StateSize());` (line 364) | Plan dropped the constructor argument. |
| R-002 | `ode_solver.GetNumDtRejects()` (no such API) | counts `!ode_solver.Step(...)` returns into a local int (lines 534, 588, 590) | Plan invented a getter. |
| R-003 | `ListAttrs(bdr_attrs)` (no such helper) | doesn't exist in BP5 (BP5 doesn't pre-validate attrs) | Plan added a SAFS-specific check then forgot a sub-helper. |
| R-008 (new) | `ode_solver.Step(state, t, dt)` (3 args) | `ode_solver.Step(seas_op, state, t, dt)` (4 args, line 590) | Plan dropped the operator argument. |

In other words, the answer to "why we have these deviations if the goal is
to mirror BP5" is: **they were NOT supposed to exist; they are not design
choices, they are transcription errors**. The fix in every case is to
rewrite the snippet to literally match the BP5 reference.

R-003 is slightly different: BP5 doesn't have an analogous check at all.
SAFS legitimately needs one because its tags (1=xm, 2=xp, ..., 100=fault,
10=domain) differ from BP5's hardcoded numbering (1=Natural, 3=Fault,
5=Dirichlet). So adding `AssertSafsAttrs` is reasonable; the bug is that
its body uses `ListAttrs` without defining it.

---

## Findings

### [R-001] [CRITICAL] [PLAN_safs_test.md:501-502] — `Vector state` not pre-sized before `SetInitialCondition`

**Category:** BUG (compile / runtime — abort at startup)

**Description:**
Plan code (lines 501–502):

```cpp
Vector state;
seas_op.SetInitialCondition(state);
```

`SEASQuasiDynamicOperator::SetInitialCondition` (`solver/seas_operator.hpp:217`)
opens with:

```cpp
MFEM_VERIFY(state.Size() == fault_->StateSize(),
            "State vector size mismatch: got " << state.Size()
            << ", expected " << fault_->StateSize());
```

A default-constructed `Vector` has size 0, so this verify fires immediately.

BP5 reference (`drivers/seas_driver.cpp:364`) is the canonical pattern:

```cpp
Vector state(fault_op.StateSize());
seas_op.SetInitialCondition(state);
```

**Trigger:** running the driver as written.

**Actual behavior:** abort at startup with
`"State vector size mismatch: got 0, expected <N>"`.

**Expected behavior:** state pre-sized like BP5.

**Suggested fix:**

```diff
-    // 9. Initial state from operator.
-    Vector state;
-    seas_op.SetInitialCondition(state);
+    // 9. Initial state from operator (mirrors BP5 driver line 364).
+    Vector state(fault_op.StateSize());
+    seas_op.SetInitialCondition(state);
```

**Test case:**
```python
def test_R001_state_vector_pre_sized_before_set_initial_condition():
    pmesh = make_minimal_safs_mesh()
    domain = ElasticityDomainOperator(pmesh, order=1, ...)
    bp5 = BP5Params()
    OverrideToUniformVW(bp5, SafsTestParams())
    fault_geom = FaultGeometry(domain, bp5, mpi_ctx)
    friction = DieterichRuinaFriction(...)
    aging = AgingLawPsi(...)
    fault_op = RateStateFaultOperator(fault_geom, friction, aging, bp5, mpi_ctx)
    seas_op = SEASQuasiDynamicOperator(domain, fault_op, mpi_ctx)
    state = Vector(fault_op.StateSize())   # required pre-size
    seas_op.SetInitialCondition(state)     # must NOT throw
    assert state.Size() == fault_op.StateSize()
```

---

### [R-002] [CRITICAL] [PLAN_safs_test.md:533] — `ode_solver.GetNumDtRejects()` does not exist

**Category:** BUG (compile-time)

**Description:**
Plan diagnostic line 533 reads:

```cpp
n_dt_rejects = ode_solver.GetNumDtRejects()
```

`DormandPrinceRK45` (`solver/time_stepper.hpp`) maintains `int total_rejections_;`
as a *private* member (line 627) and increments it inside `Step` (lines 326,
334, 351, 358, 376, 383, 401, 408, 427) but **exposes no public getter** by
any name (`grep -rn 'GetNumDtRejects\|GetNumRejects\|TotalRejections'
miniapps/seas/solver/time_stepper.hpp` returns nothing).

BP5 reference (`drivers/seas_driver.cpp` lines 534, 587–590):

```cpp
int step_rejections = 0;        // line 534
...
real_t dt;
bool accepted = ode_solver.Step(seas_op, state, t, dt);
if (!accepted) { continue; }    // increment counter here for SAFS
step++;
```

BP5 simply doesn't put rejects into a public getter — it relies on the
boolean return of `Step` (`time_stepper.hpp:226`: "@return true if step
accepted, false if rejected"). The plan invented a non-existent API.

**Trigger:** compiling the driver as written.

**Actual behavior:** compile error `'class DormandPrinceRK45' has no
member named 'GetNumDtRejects'`.

**Expected behavior:** count rejections in the driver via the boolean
return of `Step`.

**Suggested fix:** mirror BP5's pattern with a driver-local counter and an
inner accept loop (so `step` increments only on accepted steps):

```diff
     real_t t = 0.0;
     int step = 0;
+    int n_dt_rejects = 0;
     while (step < n_steps_max) {
         real_t dt = ode_solver.GetDt();
-        ode_solver.Step(state, t, dt);
-        ++step;
+        // Mirror BP5 driver lines 587-590: ode_solver.Step's bool return
+        // is the only rejection signal.  Loop until at least one
+        // accepted step is produced; count intermediate rejects.
+        bool accepted = false;
+        while (!accepted) {
+            accepted = ode_solver.Step(seas_op, state, t, dt);
+            if (!accepted) { ++n_dt_rejects; }
+        }
+        ++step;
```

Then write `n_dt_rejects` (driver-local) to the CSV instead of the
non-existent getter.

**Test case:**
```python
def test_R002_dt_rejects_counted_in_driver_not_solver():
    csv = run_driver(n_steps=5)
    assert all(csv["n_dt_rejects"] >= 0)
    assert "GetNumDtRejects" not in compile_log()
```

---

### [R-003] [CRITICAL] [PLAN_safs_test.md:583] — `ListAttrs()` referenced inside `MFEM_VERIFY` but not defined

**Category:** BUG (compile-time)

**Description:**
`AssertSafsAttrs` (Phase-3 §3.3 lines 568–595) is a SAFS-specific check
that verifies `pmesh.bdr_attributes` contains tags 1..6. Its error
message is:

```cpp
// Plan lines 581-584
MFEM_VERIFY(found,
            "SAFS mesh missing required boundary attribute " << a
            << " (mesh has [" << ListAttrs(bdr_attrs) << "])");
```

`grep -rn ListAttrs miniapps/seas/` returns this single reference. The
helper does not exist.

BP5 reference does not have an analogous check — BP5 mesh tags are
hardcoded (1=Natural, 3=Fault, 5=Dirichlet) and the legacy `BCMode`
constructor enforces them implicitly. So R-003 is an addition to the
plan, not a deviation from BP5; but it is still a compile-stopper.

**Trigger:** compiling the driver.

**Actual behavior:** compile error `'ListAttrs' was not declared in this
scope`.

**Expected behavior:** define the helper or inline the formatting.

**Suggested fix:** add a file-static lambda in the driver and use it:

```diff
+    // File-local helper to format Array<int> for diagnostic messages.
+    auto fmt_attrs = [](const Array<int>& a) {
+        std::ostringstream oss;
+        for (int i = 0; i < a.Size(); ++i) {
+            if (i > 0) { oss << ", "; }
+            oss << a[i];
+        }
+        return oss.str();
+    };
     ...
-        MFEM_VERIFY(found,
-                    "SAFS mesh missing required boundary attribute " << a
-                    << " (mesh has [" << ListAttrs(bdr_attrs) << "])");
+        MFEM_VERIFY(found,
+                    "SAFS mesh missing required boundary attribute " << a
+                    << " (mesh has [" << fmt_attrs(bdr_attrs) << "])");
```

**Test case:**
```python
def test_R003_assert_safs_attrs_compiles_and_reports_attrs():
    out = run_driver_expect_fail(mesh="missing_ztop.msh")
    assert "SAFS mesh missing required boundary attribute 5" in out
    assert "(mesh has [" in out
    body = out.split("(mesh has [")[1].split("])")[0]
    assert all(t.strip().isdigit() for t in body.split(","))
```

---

### [R-008] [CRITICAL] [PLAN_safs_test.md:526] — `ode_solver.Step(state, t, dt)` has wrong arity

**Category:** BUG (compile-time, NEW IN ROUND 2)

**Description:**
Plan time-step loop (line 526):

```cpp
real_t dt = ode_solver.GetDt();
ode_solver.Step(state, t, dt);
++step;
```

`DormandPrinceRK45::Step` signature (`solver/time_stepper.hpp:227`):

```cpp
bool Step(TimeDependentOperator &op, Vector &state, real_t &t, real_t &dt)
```

Step requires **four** arguments (operator + state + t + dt). The plan's
3-argument call `Step(state, t, dt)` will fail to compile.

BP5 reference (`drivers/seas_driver.cpp:590`):

```cpp
bool accepted = ode_solver.Step(seas_op, state, t, dt);
```

— passes `seas_op` as the first argument explicitly. The plan author
likely confused this with the MFEM `ODESolver::Step` style which
captures the operator from a prior `Init(...)` call. `DormandPrinceRK45`
re-takes the operator on every `Step` invocation; the prior `Init` only
allocates stage vectors (`time_stepper.hpp:205-214`).

**Trigger:** compiling the driver.

**Actual behavior:** compile error along the lines of
`no matching function for call to DormandPrinceRK45::Step(Vector&,
real_t&, real_t&)`.

**Expected behavior:** pass `seas_op` as first arg, mirroring BP5.

**Suggested fix:** combined with R-002 above into one corrected loop:

```diff
     while (step < n_steps_max) {
         real_t dt = ode_solver.GetDt();
-        ode_solver.Step(state, t, dt);
-        ++step;
+        bool accepted = false;
+        while (!accepted) {
+            accepted = ode_solver.Step(seas_op, state, t, dt);
+            if (!accepted) { ++n_dt_rejects; }
+        }
+        ++step;
```

(The `seas_op` is the same `PBP5SEASOp` constructed at plan line 497.)

**Test case:**
```python
def test_R008_step_passes_operator_explicitly():
    # Compile must succeed; instrument Step to record that the
    # operator argument is `seas_op`.
    log = run_driver_with_step_trace(n_steps=2)
    assert "Step(seas_op, state, t, dt)" in compile_log()
```

---

### [R-004] [MODERATE] [PLAN_safs_test.md:312-349 vs 358] — `bp5.nucleation_eps = 0.0` promised in narrative but missing from `OverrideToUniformVW` body

**Category:** DEVIATION (plan-internal inconsistency)

**Description:**
Narrative §2.3 line 358: "Set `bp5.nucleation_eps = 0.0` to be safe."
The actual `OverrideToUniformVW` body in §2.2 (lines 312–349) sets
`w_nuc`, `hs`, `ht`, `H`, `l_vw`, `Wf`, `lf`, but does NOT set
`nucleation_eps` (default `1.0e-3` at `config/bp5_params.hpp:140`).

The override is empirically still safe because `hs = 1.0e7` makes the
`x3` condition fail by seven orders of magnitude. But narrative-vs-code
mismatch is exactly the class of issue user-memory
`feedback_complete_sign_sites.md` warns about: silent drift between
documentation and implementation.

**Trigger:** any reader cross-referencing §2.2 with §2.3.

**Actual behavior:** override body does NOT set `nucleation_eps`;
narrative claims it does.

**Expected behavior:** set it. Match.

**Suggested fix:**
```diff
     bp5.smooth_nucleation = false;
+    bp5.nucleation_eps = 0.0;   // narrative §2.3 promises this
```

**Test case:**
```python
def test_R004_nucleation_eps_zeroed_in_override():
    bp5 = BP5Params()                                # default 1e-3
    OverrideToUniformVW(bp5, SafsTestParams())
    assert bp5.nucleation_eps == 0.0
```

---

### [R-005] [MODERATE] [PLAN_safs_test.md:96-97, 511-513] — `dt_init` "clamped to dt_max" claim is false for default values

**Category:** DEVIATION (documentation vs code)

**Description:**
Plan §Constraints (lines 96–97):

> Initial dt = `0.01 * L0 / Vp` = 1.4e6 s ~ 16 days, but **clamped to
> `dt_max`** = 0.1 yr ~ 3.156e6 s.

Phase-3 implementation (lines 511–513):
```cpp
real_t dt_init = (cli_dt_init > 0.0) ? cli_dt_init :
                 std::min(0.01 * params.L0 / params.Vp,
                          0.1 * SafsTestParams::seconds_per_year);
```

Numerically `0.01 * 0.14 / 1e-9 = 1.4e6 s` and `0.1 * 365.25 * 86400 =
3.155e6 s`; since `1.4e6 < 3.155e6` the `min` returns the first arg and
the clamp **never fires**. The narrative is misleading.

**Suggested fix:**
```diff
-- Time stepping: `DormandPrinceRK45` with BP5 tolerances (`atol=1e-7`,
-  `rtol=1e-50`, `dt_min=1e-6`, `dt_max=0.1 yr`). Initial dt = `0.01 * L0 / Vp`
-  = 1.4e6 s ≈ 16 days, but **clamped to `dt_max`** = 0.1 yr ≈ 3.156e6 s.
+- Time stepping: `DormandPrinceRK45` with BP5 tolerances (`atol=1e-7`,
+  `rtol=1e-50`, `dt_min=1e-6`, `dt_max=0.1 yr`). Initial dt =
+  `min(0.01 * L0 / Vp, dt_max)` = 1.4e6 s ≈ 16 days under default
+  parameters; the `dt_max` cap is a safety net for future overrides
+  that increase `L0/Vp` past 36 days.
```

**Test case:**
```python
def test_R005_default_dt_init_is_not_clamped():
    p = SafsTestParams()
    dt_default = min(0.01 * p.L0 / p.Vp, 0.1 * p.seconds_per_year)
    assert abs(dt_default - 1.4e6) < 1.0
    assert dt_default < 0.1 * p.seconds_per_year   # cap not active
```

---

### [R-006] [MODERATE] [PLAN_safs_test.md:511-513] — `dt_init` formula uses `params.Vp` rather than `max(V_init, V_nuc)`

**Category:** ASSUMPTION (silent breakage if invariant changes)

**Description:**
CLAUDE.md mandates `dt_init = 0.01 * L_nuc / V_nuc` ("Too large -> RK45
stage amplification during nucleation, debug v7"). The plan uses
`0.01 * L0 / Vp` (lines 96, 512), which is correct **only under the
implicit invariant `V_nuc == V_init == Vp` and `L_nuc == L0`** that the
override establishes.

If a future tweak (e.g. someone adds a small overstress in
`OverrideToUniformVW` and sets `V_nuc > V_init`) breaks the invariant,
the dt_init formula becomes silently wrong — over-large dt for a now-
fast nucleation seed → RK45 stage amplification.

**Suggested fix:** reference the actually-active velocity scale:

```diff
-    real_t dt_init = (cli_dt_init > 0.0) ? cli_dt_init :
-                     std::min(0.01 * params.L0 / params.Vp,
-                              0.1 * SafsTestParams::seconds_per_year);
+    // dt_init = 0.01 * L_min / V_max_init per CLAUDE.md.  For this plan
+    // V_init == V_nuc == Vp, so the three are equivalent — but write
+    // the safe form so future overrides do not silently break.
+    const real_t V_init_max = std::max(params.V_init, params.Vp);
+    real_t dt_init = (cli_dt_init > 0.0) ? cli_dt_init :
+                     std::min(0.01 * params.L0 / V_init_max,
+                              0.1 * SafsTestParams::seconds_per_year);
```

**Test case:**
```python
def test_R006_dt_init_uses_max_init_velocity():
    p = SafsTestParams()
    p.V_init = 1.0e-7    # someone increases V_init by 100x
    expected = 0.01 * p.L0 / max(p.V_init, p.Vp)
    actual = compute_dt_init(p)
    assert abs(actual - expected) / expected < 1e-12
```

---

### [R-009] [MODERATE] [PLAN_safs_test.md:530] — Plan diagnostic pseudocode references getters that do not exist on `seas_op`

**Category:** DEVIATION (NEW IN ROUND 2)

**Description:**
Plan lines 528–533 describe per-step CSV diagnostics:

```
// Diagnostics on `state`:
//   V_max, V_min, psi_max, psi_min from seas_op
//   traction_max from seas_op.GetTraction()
//   slip_l2 = ||slip components||
//   n_dt_rejects = ode_solver.GetNumDtRejects()
```

Available APIs (`solver/seas_operator.hpp` + `fault/rate_state_fault.hpp`):

| Diagnostic the plan asks for | API actually available |
|---|---|
| `V_max from seas_op` | `seas_op.GetMaxSlipRate()` (exists) ✓ |
| `V_min from seas_op` | NO `GetMinSlipRate()` exists |
| `psi_max from seas_op` | NO `GetMaxPsi()` exists |
| `psi_min from seas_op` | NO `GetMinPsi()` exists |
| `traction_max from seas_op.GetTraction()` | `GetTraction()` returns the full vector; max requires manual reduction |
| `slip_l2` | requires `fault_op.GetSlip(state, slip)` + manual L2 norm |

The CSV columns spec on line 116 (`step, t, dt, V_max, V_min, psi_max,
psi_min, traction_max, slip_l2, n_dt_rejects, success`) commits the
implementer to all six diagnostics, but the plan does not show how to
compute `V_min`, `psi_max`, `psi_min`, `traction_max`, or `slip_l2`.

The implementer would have to:

1. For `V_min`/`V_max` (both): iterate the state vector, extract slip-rate
   components (`state(i*StatePerNode + 0..1)`), compute global max/min
   under MPI reduction.
2. For `psi_min`/`psi_max`: extract `state(i*StatePerNode + PsiIndex)`
   per node, reduce.
3. For `traction_max`: max-reduce over `seas_op.GetTraction()`.
4. For `slip_l2`: call `fault_op.GetSlip(state, slip)` (line 585 of
   `rate_state_fault.hpp`), then compute `slip.Norml2()` and reduce.

This is a meaningful amount of code that the plan describes in one
hand-waving line. Either the plan must show the extraction code, or the
CSV column list must shrink to what the existing API directly exposes
(`V_max`, traction_max via reduction, n_dt_rejects).

**Suggested fix:** flesh out the pseudocode block in §3.2 to show all
six diagnostics explicitly. Example (mirroring BP5's
`bench_out.Write(...)` pattern lightly):

```cpp
// Diagnostics
const real_t V_max = seas_op.GetMaxSlipRate();          // existing API
const Vector& trac = seas_op.GetTraction();             // existing API

real_t local_V_min = std::numeric_limits<real_t>::infinity();
real_t local_psi_min =  std::numeric_limits<real_t>::infinity();
real_t local_psi_max = -std::numeric_limits<real_t>::infinity();
real_t local_trac_max = 0.0;
const int spn = RateStateFaultOperator<ParMesh,2>::StatePerNode;  // = 3
const int n_local = state.Size() / spn;
for (int i = 0; i < n_local; ++i) {
    real_t v0 = state(i*spn + 0), v1 = state(i*spn + 1);
    real_t psi = state(i*spn + 2);
    local_V_min  = std::min(local_V_min, std::hypot(v0, v1));
    local_psi_min = std::min(local_psi_min, psi);
    local_psi_max = std::max(local_psi_max, psi);
}
for (int k = 0; k < trac.Size(); ++k) {
    local_trac_max = std::max(local_trac_max, std::abs(trac(k)));
}
real_t V_min       = mpi.GlobalReduceMin(local_V_min);
real_t psi_min     = mpi.GlobalReduceMin(local_psi_min);
real_t psi_max     = mpi.GlobalReduceMax(local_psi_max);
real_t traction_max = mpi.GlobalReduceMax(local_trac_max);

Vector slip;
fault_op.GetSlip(state, slip);
real_t slip_local2 = slip * slip;            // local sum-of-squares
real_t slip_l2     = std::sqrt(mpi.GlobalSum(slip_local2));
```

**Test case:**
```python
def test_R009_csv_columns_populated_under_mpi():
    csv = run_driver(np=4, n_steps=5)
    for col in ["V_max", "V_min", "psi_max", "psi_min",
                "traction_max", "slip_l2"]:
        assert col in csv.columns
        assert csv[col].notna().all()
        assert (csv[col].abs() < 1e30).all()    # not Inf
```

---

### [R-007] [LOW] [PLAN_safs_test.md:613-614] — Acceptance criterion ambiguous: "5 CSV rows"

**Category:** QUALITY (specification clarity)

**Description:**
> CSV has 5 rows after the run, with monotone non-decreasing `t`.

A typical CSV has a header plus N data rows. With `n_steps=5` the file
should have 6 lines total (1 header + 5 data) or 5 data rows without a
header. The plan does not specify which.

**Suggested fix:**
```diff
-- [ ] CSV has 5 rows after the run, with monotone non-decreasing `t`.
+- [ ] CSV has 6 lines total (1 header line plus 5 data rows) after a
+      `--n-steps 5` run; the 5 data rows have monotone non-decreasing
+      `t`.
```

**Test case:** none — pure spec clarification.

---

## Summary

- Critical issues: **4** (R-001, R-002, R-003, R-008 — all prevent the
  driver from compiling/running)
- Moderate issues: **4** (R-004 plan/code drift; R-005 misleading
  documentation; R-006 hidden invariant; R-009 missing diagnostic
  extraction code)
- Low issues: **1** (R-007 spec clarity)
- Plan compliance: **PARTIAL** — composition strategy is sound and
  correctly identifies the BoundaryConfig route + BP5Params override
  trick, but the actual code snippets diverge from the BP5 reference in
  four compile-stopping ways (R-001, R-002, R-003, R-008) and three
  documentation-vs-code drift sites (R-004, R-005, R-009). The plan was
  most likely re-written from memory rather than mechanically copied
  from `drivers/seas_driver.cpp`.
- **Verdict:** **PASS WITH FIXES**. Must fix R-001, R-002, R-003, R-008
  before the driver can be built; R-004..R-006, R-009 should be cleaned
  up before handoff.

## Direct answer to the user's question

> "but why we have R-001 to R-003 deviations if the goal is to mirror
> BP5 setup but with this SAFS mesh?"

Because the plan author rewrote the BP5 driver pattern from memory and
the resulting C++ snippets dropped or invented details. The deviations
are NOT design choices — they are transcription errors. The fix is, in
each case, to literally match the BP5 reference at
`drivers/seas_driver.cpp:364, 590, 587-590`. The single exception is
R-003 (`AssertSafsAttrs` / `ListAttrs`): BP5 has no analogue, but the
plan introduced this new helper and forgot to define a sub-helper.

Round 2 also found a fourth compile-stopper of the same class (R-008,
3-arg `Step` call) that round 1 missed.

## Unreviewed Areas

- **Phase-4 `check_safs_smoke.py`** — the plan describes acceptance
  metrics but no Python source is included. Cannot review what does not
  exist.
- **Numerical-correctness verification of the override approach.** The
  plan claims "Initial state is identically the steady-state plate-rate
  configuration — by construction nothing should evolve to first
  order." Verifying this claim requires running the driver and
  inspecting the CSV; cannot be checked from the plan alone.
- **Multi-fault FaultBasis behavior.** The plan acknowledges (Risk
  Assessment line 733) that `ref_normal = (0,-1,0)` produces
  inconsistent local frames across the 6 SAFS faults. Documented as
  acceptable shakedown limitation, not a bug.
- **`SetVGuard`** (BP5 driver line 451) — the plan does NOT call
  `ode_solver.SetVGuard(100.0)` even though BP5 reference does. With
  `V_init = Vp = 1e-9` the system should be at steady state and VGuard
  is unlikely to fire, so this is probably intentional. Not flagged as
  a bug. Mention only.
- **`SetStatePerNode`** (plan line 515) — plan calls `SetStatePerNode(3)`,
  BP5 reference does not. The default is 2 (BP2). For BP5 with VGuard
  the correct value is 3. Plan is **more correct than BP5 reference
  here**, not a bug; arguably it exposes a latent BP5 bug, but that's
  out of scope.
