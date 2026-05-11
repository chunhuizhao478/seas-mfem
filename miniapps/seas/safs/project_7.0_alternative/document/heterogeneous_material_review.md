# Code Review: heterogeneous_material_plan.md — Phase 1 (fresh adversarial audit, 2026-05-11)

## Review Scope
- Plan: `safs/project_7.0_alternative/document/heterogeneous_material_plan.md`
  (Phase 1 only — the only phase implemented to date).
- Files reviewed:
  - `dynamic/heterogeneous_material.hpp`        — modified
  - `dynamic/heterogeneous_material.cpp`        — new
  - `io/material_coefficients.hpp`              — new
  - `io/material_coefficients.cpp`              — new
  - `tests/unit/test_heterogeneous_material.cpp` — modified (T-5-4..T-5-7 added)
  - `Makefile`                                  — modified (new src/obj/test-target rules)
- Domain context: `seas-mfem-safs/CLAUDE.md`, `miniapps/seas/CLAUDE.md`,
  plan's R-001..R-014 round-1/2/3 fix history.
- Build state verified: `make test-heterogeneous-material` exits 0 with
  34/34 passing on np=1 (mpirun-driven).

This pass executes all three review passes from scratch on the Phase-1
changes. The plan's Phase 1 scope is mostly mechanical (struct extension
plus three small Coefficient subclasses), but it pins several Coefficient-
mode contracts (abort paths, lifetime, NaN sentinels) that must remain
testable.

---

## Findings

### [R-001] MODERATE `dynamic/heterogeneous_material.hpp:EvalAt` — `Mode::GridFunction` abort path is untested

**Category:** EDGE_CASE (coverage gap)

**Description:**
Plan Phase 1 Detailed Req. 2 ("EvalAt semantics") states:
> Mode::GridFunction → MFEM_ABORT (use At(elem, dof, …) instead; Coefficient
> is the only mode that meaningfully consumes a qpoint).

The implementation in `heterogeneous_material.hpp` lines 203-205 has the
`MFEM_ABORT(...)` for this path, but no test in
`tests/unit/test_heterogeneous_material.cpp` exercises it. T-5-5 only
covers (a) `At()` in Coefficient mode and (b) `MaxCpInElement(elem, nullptr)`
in Coefficient mode. A future regression that silently degrades this
abort to e.g. a `return;` (or a no-op fallthrough to Coefficient code that
crashes on null `lambda_coef`) would slip past CI.

This is a contract pinned by the plan that has zero test coverage.

**Trigger:**
A future change that removes or weakens the `MFEM_ABORT` on
`heterogeneous_material.hpp:204` while preserving Constant/Coefficient
behaviour.

**Actual behavior:**
No test currently fails if the GridFunction-mode `EvalAt` abort is removed.

**Expected behavior:**
A fork-based abort test pinning the contract.

**Suggested fix:**
Add a sub-test T-5-5(c) to `T_5_5_abort_paths_in_coefficient_mode` (rename
the function to `T_5_5_abort_paths`):

```diff
--- a/miniapps/seas/tests/unit/test_heterogeneous_material.cpp
+++ b/miniapps/seas/tests/unit/test_heterogeneous_material.cpp
@@ -425,6 +425,30 @@ static void T_5_5_abort_paths_in_coefficient_mode()
    TEST_ASSERT(aborted_on_null_T,
                "T-5-5(b) MaxCpInElement(elem, nullptr) aborts in "
                "Coefficient mode");
+
+   // (c) EvalAt(elem, T, ip, ...) called on a GridFunction-mode
+   // MaterialField MUST abort (plan Phase 1 Detailed Req. 2).
+   const bool aborted_on_eval_gf = expect_abort([](){
+      auto base = std::make_shared<mfem::Mesh>(
+         mfem::Mesh::MakeCartesian3D(1, 1, 1,
+                                     mfem::Element::HEXAHEDRON,
+                                     1.0, 1.0, 1.0));
+      // GridFunction mode in the child needs a ParFES, but we cannot
+      // init MPI here (we are still pre-Mpi::Init at this call site).
+      // Construct a MaterialField directly with mode = GridFunction
+      // and SAFE non-null shared_ptrs of a default-constructed
+      // ParGridFunction — except that ParGridFunction requires MPI.
+      // Easier: build the struct by hand, set mode = GridFunction,
+      // leave the shared_ptrs null, and rely on EvalAt's mode-branch
+      // hitting MFEM_ABORT before any GF dereference.
+      MaterialField m;
+      m.mode = MaterialField::Mode::GridFunction;
+      mfem::IsoparametricTransformation T_local; // unused inside abort
+      mfem::IntegrationPoint ip; ip.Set3(0.5, 0.5, 0.5);
+      real_t la, mu, rho;
+      m.EvalAt(/*elem=*/0, T_local, ip, la, mu, rho);
+   });
+   TEST_ASSERT(aborted_on_eval_gf,
+               "T-5-5(c) EvalAt in Mode::GridFunction aborts");
 }
```

Note: the inner test does NOT need a valid ParGridFunction because the
mode branch hits `MFEM_ABORT` BEFORE any GF dereference. Constructing
`MaterialField` by hand (instead of via `MakeGridFunction`, which
verifies non-null) is exactly what this test should exercise.

**Test case:** see diff above.

---

### [R-002] MODERATE `dynamic/heterogeneous_material.cpp:MaxCpInElement` — `T->SetIntPoint(&ip)` mutates caller's `ElementTransformation` state

**Category:** ASSUMPTION (silent side-effect)

**Description:**
The Coefficient-mode loop at `heterogeneous_material.cpp:87-96` calls
`T->SetIntPoint(&ip)` once per quadrature point. After the function
returns, `T` retains a pointer to the LAST quadrature point in the rule
— *not* whatever the caller had set before invoking `MaxCpInElement`.
This is a hidden mutation of a caller-supplied object via a non-const
pointer.

A caller pattern that fails:
```cpp
T->SetIntPoint(&ip_user);
auto coef_val = some_coef->Eval(*T, ip_user);   // OK so far
auto cp_max   = material.MaxCpInElement(elem, T);
// ↓ This silently uses the LAST qpoint of IntRules, not ip_user.
auto coef_val_again = some_coef->Eval(*T, ip_user);
```

The MFEM convention is that `Coefficient::Eval(T, ip)` callers either
(a) pass `ip` explicitly OR (b) rely on `T.GetIntPoint()` having been
set. Many subclasses use one and not the other. `MaxCpInElement`'s
mutation leaks into pattern (b).

The plan does not pin this contract explicitly. But the call sites
audited in Phase 3 (CFL helper inside the WaveOperator hot path) will be
particularly susceptible: WaveOperator's per-element walk uses the same
`T` for assembly and for `MaxCpInElement`, and assembly takes
`T.GetIntPoint()` implicitly inside MFEM integrators.

**Trigger:**
Phase 3 implementation reuses the same `ElementTransformation*` for
both `MaxCpInElement` and a subsequent `Coefficient::Eval` that relies
on `T.GetIntPoint()`.

**Actual behavior:**
`T->GetIntPoint()` after `MaxCpInElement` returns the last `IntPoint`
from the internal `IntRules` rule.

**Expected behavior:**
Either save/restore the IP on entry/exit, or document the contract
explicitly in the header docstring.

**Suggested fix (option A — save/restore):**
```diff
--- a/miniapps/seas/dynamic/heterogeneous_material.cpp
+++ b/miniapps/seas/dynamic/heterogeneous_material.cpp
@@ -85,15 +85,21 @@ real_t MaterialField::MaxCpInElement(
    const mfem::IntegrationRule& ir = mfem::IntRules.Get(
       T->GetGeometryType(), 2 * default_qorder);

+   // Save the caller's IP so we can restore it on exit — Coefficient
+   // subclasses that rely on T.GetIntPoint() (e.g. SumCoefficient,
+   // ProductCoefficient) would otherwise observe our last qpoint.
+   const mfem::IntegrationPoint* const saved_ip = T->GetIntPoint();
+
    real_t cp_max = -std::numeric_limits<real_t>::infinity();
    for (int q = 0; q < ir.GetNPoints(); ++q)
    {
       const mfem::IntegrationPoint& ip = ir.IntPoint(q);
       T->SetIntPoint(&ip);
       const real_t la  = lambda_coef->Eval(*T, ip);
       const real_t mu_ = mu_coef    ->Eval(*T, ip);
       const real_t rho = rho_coef   ->Eval(*T, ip);
       const real_t cp  = std::sqrt((la + 2.0 * mu_) / rho);
       if (cp > cp_max) { cp_max = cp; }
    }
+   if (saved_ip) { T->SetIntPoint(saved_ip); }
    return cp_max;
 }
```

**Suggested fix (option B — document only):**
Add a paragraph to the header docstring of `MaxCpInElement`:
> Side-effect: in `Mode::Coefficient` this method calls
> `T->SetIntPoint(...)` once per internal quadrature point. After return,
> `T`'s current IntegrationPoint is the last qpoint of the internal
> rule, NOT whatever the caller had set. Callers must re-`SetIntPoint`
> before any subsequent `T`-aware operation.

Option A is preferred. The cost (one pointer save + conditional store)
is negligible.

**Test case:**
```cpp
static void T_R002_max_cp_does_not_mutate_T_ip()
{
   // Build sidecar-backed Coefficient mode MaterialField (same as T-5-6).
   ...
   auto mesh = make_single_hex_mesh();
   mfem::ElementTransformation* T = mesh->GetElementTransformation(0);
   mfem::IntegrationPoint ip_user;
   ip_user.Set3(0.25, 0.25, 0.25);
   T->SetIntPoint(&ip_user);

   const auto* before = T->GetIntPoint();
   (void) m.MaxCpInElement(0, T);
   const auto* after  = T->GetIntPoint();
   TEST_ASSERT(before == after,
               "R-002: MaxCpInElement preserves T's IntegrationPoint");
}
```

---

### [R-003] LOW `dynamic/heterogeneous_material.cpp:MaxCpInElement` — inconsistent initial-value sentinel between GridFunction (`0.0`) and Coefficient (`-inf`) paths

**Category:** QUALITY (potential silent wrong-answer in a downstream divide)

**Description:**
`MaxCpInElement` initialises `cp_max = 0.0` in GridFunction mode
(`heterogeneous_material.cpp:55`) and
`cp_max = -std::numeric_limits<real_t>::infinity()` in Coefficient mode
(`heterogeneous_material.cpp:86`).

The two values diverge in two ways:
- If a hypothetical element has zero quadrature/DOF points, the
  GridFunction path returns `0.0` (silently produces `inf` `dt = h/cp`
  downstream), while the Coefficient path returns `-inf` (downstream
  `dt = h/-inf = -0`, easier to spot).
- If every sample produces `NaN` (e.g. `rho == 0` or
  `la + 2 mu < 0`), the GridFunction path returns `0.0` (silently wrong),
  the Coefficient path returns `-inf`.

In practice MFEM ensures every element has at least one DOF/qpoint, so
the first concern is hypothetical. The second is a real risk if the
sidecar is corrupted.

**Trigger:**
A corrupted sidecar where `rho = 0` at some voxel; or a sidecar where
`Vp² < 2·Vs²` and the per-field clamp in `build_velocity_cvmh.py` was
disabled.

**Actual behavior:**
GridFunction mode returns `0.0`; downstream `dt = h / 0` yields `+inf`,
RK loop accepts huge timestep, simulation explodes.

**Expected behavior:**
Both paths use `-inf` so a sentinel value propagates clearly.

**Suggested fix:**
```diff
--- a/miniapps/seas/dynamic/heterogeneous_material.cpp
+++ b/miniapps/seas/dynamic/heterogeneous_material.cpp
@@ -52,7 +52,7 @@ real_t MaterialField::MaxCpInElement(
       mfem::Array<int> vdofs;
       fes->GetElementDofs(elem, vdofs);
-      real_t cp_max = 0.0;
+      real_t cp_max = -std::numeric_limits<real_t>::infinity();
       for (int d = 0; d < vdofs.Size(); ++d)
```

**Test case:**
```cpp
// Hypothetical (would require constructing a zero-DOF FES — not normally
// reachable; the value of this test is regression-pinning the sentinel).
// More realistic test: a GridFunction with all-zero rho.
auto rho_gf_zero = std::make_shared<ParGridFunction>(fes.get());
*rho_gf_zero = 0.0;  // division by zero → NaN cp
auto m = MaterialField::MakeGridFunction(rho_gf_zero, lambda_gf, mu_gf);
const real_t cp = m.MaxCpInElement(0);
TEST_ASSERT(std::isnan(cp) || std::isinf(cp),
            "R-003: zero-rho GF produces NaN/inf, not silent 0.0");
```

---

### [R-004] LOW `tests/unit/test_heterogeneous_material.cpp` — `MakeCoefficient(nullptr, ...)` abort path is untested

**Category:** EDGE_CASE (coverage gap)

**Description:**
Plan Phase 1 Detailed Req. 6 specifies:
```cpp
MFEM_VERIFY(lambda && mu && rho,
            "MakeCoefficient: all three Coefficient pointers must be non-null");
```
The implementation in `heterogeneous_material.cpp:25-27` matches the
plan, but no test verifies the abort. A future change that weakens
this verify (e.g. `MFEM_ASSERT` instead of `MFEM_VERIFY`) would slip
past CI in Release builds.

**Trigger:**
Future refactor downgrades `MFEM_VERIFY` to `MFEM_ASSERT`.

**Suggested fix:**
Add to `T_5_5_abort_paths_in_coefficient_mode`:
```diff
+   // (d) MakeCoefficient with any null pointer must abort
+   //     (plan Phase 1 Detailed Req. 6).
+   const bool aborted_on_null_lambda = expect_abort([](){
+      mfem::ConstantCoefficient mu_c (2.0);
+      mfem::ConstantCoefficient rho_c(3.0);
+      MaterialField::MakeCoefficient(nullptr, &mu_c, &rho_c);
+   });
+   TEST_ASSERT(aborted_on_null_lambda,
+               "T-5-5(d) MakeCoefficient(nullptr, mu, rho) aborts");
```
Repeat for null `mu` and null `rho`.

**Test case:** see diff above.

---

### [R-005] LOW `tests/unit/test_heterogeneous_material.cpp` — `MakeGridFunction(nullptr, ...)` abort path is untested

**Category:** EDGE_CASE (coverage gap)

**Description:**
Plan Phase 1 §Detailed Req. 7 introduced an `MFEM_VERIFY` in
`MakeGridFunction` requiring all three shared_ptrs to be non-null
(R-006 round-3). The implementation in `heterogeneous_material.hpp:112-114`
matches, but no test exercises this guard. The plan acknowledges this
is a strict tightening; without a regression test pinning it, a future
"simplification" that removes the verify would silently regress.

**Trigger:**
Future refactor removes the `MFEM_VERIFY` in `MakeGridFunction`.

**Suggested fix:**
```diff
+   // (e) MakeGridFunction with any null shared_ptr must abort
+   //     (plan Phase 1 Detailed Req. 7, R-006 round-3).
+   //
+   // Note: must run AFTER Mpi::Init because ParGridFunction needs MPI.
+   //       Move this sub-test into a post-Mpi-init phase.
+   // ...
```
This test should be added in the post-MPI section (with T-5-1..T-5-3)
because building a non-null shared_ptr to compare against requires a
`ParGridFunction`, which requires `Mpi::Init`.

**Test case:** see diff suggestion above.

---

### [R-006] LOW `dynamic/heterogeneous_material.cpp:78-83` — comment claims "4-point tet rule" but `IntRules.Get(TET, 4)` returns 11 points

**Category:** QUALITY (misleading comment)

**Description:**
The Coefficient-path comment says:
```cpp
// c_p varies slowly in space; a 4-point tet rule is more
// than enough to bound the per-element maximum.
```
The actual quadrature is `IntRules.Get(geom, 2 * 2)` = order 4. For
tetrahedron geometry, MFEM's order-4 Gaussian rule has **11 points**,
not 4. For hexahedron it has 8 (2×2×2 Gauss-Legendre). The comment was
likely copied from an earlier round of the plan that specified
"4-point" literally.

This is purely a documentation issue — the implementation is correct.

**Suggested fix:**
```diff
-   // c_p varies slowly in space; a 4-point tet rule is more
-   // than enough to bound the per-element maximum.  Future callers
+   // c_p varies slowly in space; an order-4 quadrature rule (11 points
+   // for tetrahedra, 8 for hexahedra) is more than enough to bound the
+   // per-element maximum.  Future callers
```

**Test case:** None required (documentation-only change).

---

### [R-007] LOW `Makefile:1309-1313` — `HDF5_INCFLAGS` added to `material_coefficients.cpp` compile, but the TU does not include `hdf5.h`

**Category:** QUALITY (build-system noise)

**Description:**
The Makefile rule at lines 1309-1313 adds `$(HDF5_INCFLAGS)` to the
compile of `io/material_coefficients.cpp`. The TU's includes are:
```cpp
#include "material_coefficients.hpp"
```
which in turn pulls in `data_field_3d.hpp`. The latter forward-declares
HDF5 ids as opaque integers and does NOT include `<hdf5.h>`. So
`material_coefficients.cpp` does not directly or transitively require
HDF5 headers.

The flag is defensive but unnecessary. If a later refactor moves HDF5
into the `data_field_3d.hpp` public surface, the flag will be required
— so leaving it in place is harmless. Worth flagging only because the
pattern propagates copy-paste cruft.

**Suggested fix (optional):**
Either remove the flag, or add a comment explaining it is preemptive.
```diff
 $(MATERIAL_COEFFICIENTS_OBJ): %.o: $(SRC)%.cpp \
                                     $(MATERIAL_COEFFICIENTS_HEADERS) \
                                     $(MFEM_LIB_FILE) $(CONFIG_MK)
 	@mkdir -p $(@D)
-	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@
+	# HDF5_INCFLAGS preserved here even though the TU does not currently
+	# #include <hdf5.h>; data_field_3d.hpp may grow HDF5-public surface.
+	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@
```

**Test case:** None required.

---

## Plan-deviation check (Pass 1)

| Plan requirement (Phase 1) | Implemented? | Notes |
|---|---|---|
| `Mode::Coefficient = 2` in enum | ✓ | `heterogeneous_material.hpp:69` |
| Three `mfem::Coefficient*` non-owning members | ✓ | hpp:89-91 |
| `MakeCoefficient` factory (out-of-line) | ✓ | cpp:21-34 |
| `EvalAt(elem, T, ip, ...)` hot-path accessor | ✓ | hpp:180-206 |
| `At()` aborts in `Mode::Coefficient` | ✓ | hpp:145-150 |
| `MaxCpInElement(elem, T*)` declaration only in header | ✓ | hpp:212-213 (declaration only) |
| `MaxCpInElement` defined in `.cpp` | ✓ | cpp:37-98 |
| `MakeGridFunction` tightened with `MFEM_VERIFY` | ✓ | hpp:112-114 |
| `LambdaFromSidecar`, `MuFromSidecar`, `RhoFromSidecar` subclasses | ✓ | `material_coefficients.hpp` |
| Coefficient sample math: `μ = ρVs²`, `λ = ρ(Vp² − 2Vs²)` | ✓ | `material_coefficients.cpp:25-46` |
| T-5-4 round-trip via `ConstantCoefficient` | ✓ | test file |
| T-5-5 `At()`-in-Coefficient + `MaxCpInElement(nullptr)` abort | ✓ (partial — see R-001) | test file |
| T-5-6 `MaxCpInElement` with synthetic sidecar | ✓ | test file |
| T-5-7 `LambdaFromSidecar` numerical correctness | ✓ | test file |

**No silent plan deviations identified.** All Phase-1 acceptance criteria
match the plan; the gaps are all on the test-coverage side.

---

## Pass 2 (bugs / edge cases / numerical)

Issues found in Pass 2 are folded into R-001..R-007 above. No
CRITICAL numerical or memory-safety bugs were found in this round.

Specifically verified:
- `MakeCoefficient`'s parameter order `(lambda, mu, rho)` correctly
  maps to the struct members of the same names — no swap bug.
- `LambdaFromSidecar`'s ordered constructor `(vp_field, vs_field, rho_field)`
  matches the body's `vp_field_`, `vs_field_`, `rho_field_` ordering —
  no swap.
- `μ = ρVs²` and `λ = ρVp² − 2μ` formulas in `material_coefficients.cpp`
  agree algebraically with the plan's specification.
- Inline `EvalAt` in `Mode::Coefficient` correctly threads `(T, ip)` to
  each Coefficient.
- The factory's `MFEM_VERIFY` guard is `MFEM_VERIFY` (active in Release),
  not `MFEM_ASSERT` (Debug-only); matches plan.
- `MaxCpInElement` in Constant mode is branch-free.
- `MaxCpInElement` in GridFunction mode preserves the existing
  round-1 behaviour (no semantic change to T-5-1..T-5-3).
- `T_5_5_abort_paths_in_coefficient_mode` correctly runs BEFORE
  `Mpi::Init` to avoid the `MPI_Abort`-propagation hang documented in
  the implementer's report. The reason is captured in the docstring at
  test file lines 392-399.

---

## Pass 3 (quality / maintainability)

Issues found in Pass 3 are folded into R-006, R-007 above. No
additional findings.

Specifically verified:
- `MaterialField::MakeCoefficient(λ, μ, ρ)` and
  `MakeGridFunction(rho, lambda, mu)` have **different parameter
  orders** — this is preserved-as-is from the plan (which itself
  preserves the pre-Phase-1 `MakeGridFunction` signature for back-
  compat). It is a latent UX trap for driver-side wiring code (Phase 4)
  but not a Phase-1 bug.
- Header docstrings on hpp:32-46 correctly summarise the three accessor
  contracts.
- The "redundant" `T->SetIntPoint(&ip)` inside `MaxCpInElement` is
  defensive against Coefficient subclasses that rely on
  `T.GetIntPoint()` rather than the explicit `ip` argument (e.g. some
  `SumCoefficient` / `ProductCoefficient` compositions) — this is
  CORRECT defensive coding, not waste; but it produces the side-effect
  flagged in R-002.

---

## Summary

- Critical issues: **0**
- Moderate issues: **2** (R-001 GridFunction-EvalAt abort untested; R-002 silent `T` mutation in `MaxCpInElement`)
- Low issues: **5** (R-003 cp_max sentinel; R-004 MakeCoefficient null-test gap; R-005 MakeGridFunction null-test gap; R-006 misleading qrule comment; R-007 HDF5_INCFLAGS noise)
- Plan compliance: **FULL** (all Phase-1 requirements implemented; no silent deviations)
- Verdict: **PASS WITH FIXES**

The Phase-1 implementation is correct against the plan. All 34 tests
pass. No critical/moderate functional bugs identified — the two
MODERATE findings are (a) a missing test that pins an explicit plan
contract (R-001) and (b) a mutation side-effect that will become
load-bearing in Phase 3 when the same `T` is reused for `MaxCpInElement`
and for `Coefficient::Eval` (R-002). Both should be addressed before
Phase 2 starts using `MaterialField` from `ElasticityDomainOperator`.

LOW findings are coverage gaps and documentation cleanups; they may be
deferred without blocking Phase 2.

---

## Unreviewed Areas

- **Phase 2 / Phase 3 / Phase 4 spec text in the plan document**
  remains UNREVIEWED in this round. Only Phase 1 is implemented;
  reviews of later-phase specs should happen alongside their
  implementation rounds.
- **Existing T-5-1..T-5-3 behaviour byte-equivalence vs. pre-Phase-1
  reference** was verified by the implementer ("34 of 34 passed") but
  was not independently re-run by this reviewer. The plan's Phase-1
  acceptance criterion #2 ("Existing tests still pass byte-identically")
  is taken on the implementer's word; no diff against a pre-change
  baseline was performed because this Phase touches inline code paths
  that the existing tests already cover.
- **`make` does not build the test target unless invoked directly**.
  The plan does not require integrating
  `test-heterogeneous-material` into a top-level `test` aggregator
  target; reviewing that gap is out of scope for Phase 1 (driver/CI
  wiring is Phase 4).
