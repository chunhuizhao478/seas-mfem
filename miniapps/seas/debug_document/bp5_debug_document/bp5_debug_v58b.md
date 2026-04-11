# BP5 Debug v58b: Plan for Safe, Standalone Diagnostics

**Date:** 2026-04-04  
**Status:** Planning document  
**Objective:** Redesign BP5 diagnostics so they are observational, standalone, and incapable of perturbing the live solve path

---

## 1. Problem Summary

The recent BP5/Tandem investigation showed that diagnostic code can easily become part of the bug:

- `dd73ac4` introduced the `[MFEM-TQ]` exact-face diagnostic and the run started stalling
- later changes proved the stall was **not** due to:
  - the one-time print block body
  - the decomposed-to-fast-path transition
  - the later sign-chain cleanup
- the likely trigger is that the diagnostic enabled a different internal traction path inside the live solve

This is the core design failure:

- the diagnostic was **not read-only**
- the diagnostic changed control flow (`need_decomp`)
- the diagnostic reused mutable transform/stateful objects inside the live operator
- the diagnostic logic was spread across:
  - `bp5_verification_full.cpp`
  - `elasticity_operator.hpp`
  - `seas_operator.hpp`
  - `rate_state_fault.hpp`
  - `fault_geometry.hpp`

So the immediate lesson is:

> BP5 diagnostics must be designed as **standalone observers of accepted state**, not as optional branches inside the production solve path.

---

## 2. Design Goals

The new diagnostic system should satisfy all of the following:

1. **No solver-path mutation**
   - enabling a diagnostic must not change:
     - quadrature choice
     - traction path
     - decomposition path
     - time-step acceptance logic
     - MPI communication pattern

2. **No re-entry into the live coupled solve**
   - no extra `domain.Solve(...)`
   - no extra `domain.ComputeTraction(...)`
   - no extra `fault_op.ComputeRHS(...)`
   - no post-step replay of expensive kernels inside the production run

3. **Accepted-state only**
   - diagnostics should operate only on a known accepted state:
     - accepted `step`
     - accepted `t`
     - accepted `dt`
     - accepted `state`
     - accepted `displacement`
     - accepted `traction`

4. **Separation of concerns**
   - collection logic, formatting, and output transport must be separate

5. **One flag, one behavior**
   - a diagnostic flag should enable one clearly defined data product
   - not “also turn on another internal algorithm”

6. **Safe in MPI**
   - no per-face interleaved direct printing in hot loops
   - no hidden collectives that only fire on some ranks

7. **Easy to remove**
   - diagnostics should live in isolated modules or callbacks
   - not mixed into core solver kernels

---

## 3. Root Cause Patterns to Avoid

From the v58 work, the unsafe patterns were:

### 3.1 Diagnostics changing algorithm choice

Example pattern:

- `need_decomp = ... || diag_tnd_tq_`

This is unsafe because the diagnostic is no longer just observing results; it is changing which code path computes them.

**Rule:** diagnostics may inspect outputs of a path, but may not choose the path.

### 3.2 Diagnostics inside hot face loops

Example pattern:

- exact-face matching
- transform mutation via `SetAllIntPoints(...)`
- local face shape evaluation
- per-QP formatted printing

inside `ComputeTractionImpl(...)`

This is unsafe because:

- it couples diagnostics to internal mutable objects
- it increases the chance of stateful transform misuse
- it makes MPI output interleaving likely
- it is hard to reason about side effects

**Rule:** hot kernels should only compute and return data, never print or run ad hoc forensic logic.

### 3.3 Diagnostics implemented as extra live solves

Example pattern:

- recomputing `Solve + ComputeTraction` after accepted PETSc steps

This is unsafe because:

- the replay may not be reentrant
- the replay may perturb solver state or expectations
- the replay can be confused with stage-state logic

**Rule:** no extra live replay in the production executable.

### 3.4 Diagnostics scattered across unrelated layers

The v58 tip/debug hooks were spread across:

- driver
- elasticity domain operator
- friction operator
- fault geometry
- SEAS operator

This makes it difficult to know which flags affect runtime and which are inert.

**Rule:** diagnostic features must have a single owner layer and explicit interfaces.

---

## 4. Proposed Safe Architecture

The proposed design has three layers:

### 4.1 Layer A: Production solver stays pure

Production code computes only:

- accepted state update
- displacement
- traction
- normal traction
- slip rate

No exact-face forensic logic should run here.

The production interfaces should expose only stable read-only outputs:

- accepted `state`
- `slip_owned`
- `slip_local`
- `displacement`
- `traction_owned`
- `traction_local`
- `normal_traction_owned`

### 4.2 Layer B: Snapshot/callback layer

Add a lightweight, read-only snapshot object that is created after an accepted step.

Suggested type:

```cpp
struct BP5AcceptedStepSnapshot
{
   int step;
   real_t t;
   real_t dt;

   const Vector *state_owned;
   const Vector *slip_owned;
   const Vector *slip_local;
   const GridFunction *displacement;
   const Vector *traction_owned;
   const Vector *traction_local;
   const Vector *normal_traction_owned;
};
```

Rules:

- snapshot creation must not allocate large new solver objects unless explicitly requested
- snapshot is read-only
- the callback sees only accepted-state data

### 4.3 Layer C: Standalone diagnostic extractors

All heavy forensic logic moves into explicit helper functions/classes.

Examples:

```cpp
FaceTractionRecord ExtractFaceTractionRecord(
   const ElasticityDomainOperator<MeshType>& domain,
   const BP5AcceptedStepSnapshot& snap,
   int face_idx);

TipMonitorRecord ExtractTipMonitorRecord(
   const FaultOperator& fault,
   const BP5AcceptedStepSnapshot& snap,
   int dof_idx);
```

Properties:

- pure functions or logically pure helpers
- no printing
- no global flags
- no internal control-flow mutation
- output is a data record only

---

## 5. Two Diagnostic Modes

Diagnostics should be split into two modes with different safety levels.

### 5.1 Mode 1: Production-safe monitoring

Purpose:

- low-cost monitoring during real runs

Allowed:

- accepted-step scalar summaries
- selected DOF/station values
- already-computed vectors

Examples:

- `V_max`
- station traction
- station jump residual if already available
- a few selected fault DOF values

Not allowed:

- exact-face per-QP reconstruction
- custom face matching loops over all fault faces
- extra solves

### 5.2 Mode 2: Offline forensic extraction

Purpose:

- exact-face MFEM vs Tandem comparison
- per-QP face pipeline audit

This should run in a **separate executable** or a dedicated postprocess mode:

- input:
  - checkpointed accepted state
  - mesh
  - parameter file
- output:
  - face-level records like `[MFEM-TQ]`

This mode may call heavy face-by-face reconstruction, because it is no longer on the production time-stepping path.

This is the right home for:

- exact-face key matching
- per-QP `u^-`, `u^+`, `slip`, `jump`, `Ty_stress`, `Ty_penalty`, `Ty`
- comparison-target formatting matching Tandem

---

## 6. Recommended Refactor Targets

### 6.1 Remove v58 investigation hooks from live BP5 path

Priority: immediate

Remove or isolate:

- `[MFEM-TQ]` path from `ElasticityDomainOperator::ComputeTractionImpl`
- accepted-step replay logic in `bp5_verification_full.cpp`
- `TIP-MON`, `TIP-FACE-DUMP`, `TIP-STEP1`
- one-off friction tip dumps in `RateStateFault`
- one-off fault geometry tip dumps
- one-off SEAS post-init tip dumps

These should not remain in the production driver/operator path.

### 6.2 Introduce a post-step observer interface

Priority: high

Suggested interface:

```cpp
class BP5StepObserver
{
public:
   virtual ~BP5StepObserver() = default;
   virtual void OnAcceptedStep(const BP5AcceptedStepSnapshot& snap) = 0;
};
```

The driver can own zero or more observers.

Benefits:

- explicit lifecycle
- accepted-state only
- testable
- no hidden interaction with solver internals

### 6.3 Build a standalone face extractor module

Priority: high

Suggested files:

- `miniapps/seas/diagnostics/face_traction_record.hpp`
- `miniapps/seas/diagnostics/face_traction_extractor.hpp`
- `miniapps/seas/diagnostics/tip_monitor.hpp`

These modules should:

- accept explicit inputs
- return records
- never print directly

### 6.4 Centralize output formatting

Priority: medium

Suggested files:

- `miniapps/seas/diagnostics/format_mfem_tq.hpp`
- `miniapps/seas/diagnostics/format_tip_monitor.hpp`

Formatting should be the last step, after data extraction.

This avoids:

- printing from inside operator loops
- MPI interleaving in numerical kernels

### 6.5 Add explicit MPI-safe transport

Priority: medium

For parallel diagnostics:

- gather small structured records to root
- or write rank-local files with deterministic naming

Avoid:

- direct `mfem::out` from many ranks in face loops

---

## 7. Concrete Near-Term Implementation Plan

### Phase 1: Stabilize production path

1. Remove all v58 investigation diagnostics from live runtime code.
2. Confirm normal BP5 jobs run again from the clean baseline.
3. Keep only generic pre-existing diagnostics that are known safe.

Deliverable:

- production BP5 run no longer contaminated by debug instrumentation

### Phase 2: Create accepted-step observer API

1. Add `BP5AcceptedStepSnapshot`
2. Add optional observer registration in the verification driver
3. Call observers only after accepted steps

Deliverable:

- future monitoring no longer requires edits inside the core solve path

### Phase 3: Move tip monitoring into observer

1. Reimplement `TIP-MON` as a post-step observer
2. Restrict it to values already computed by the accepted step
3. No extra solve, no extra traction recomputation

Deliverable:

- production-safe tip monitoring

### Phase 4: Build standalone exact-face forensic tool

1. Add a separate executable or driver mode for exact-face extraction
2. Input:
   - checkpoint/state
   - accepted `t`, `dt`
   - target face selector
3. Output:
   - exact-face per-QP records in stable format

Deliverable:

- safe replacement for `[MFEM-TQ]`

### Phase 5: Regression tests for diagnostic safety

Add tests that verify:

1. enabling a production-safe observer does not change numerical results
2. diagnostic flags do not change chosen algorithm path
3. no extra solve/traction call counts occur when observers are enabled

Suggested test types:

- compare state vectors with observer off/on
- compare traction with observer off/on
- count `Solve()` calls in a mock observer-enabled run

---

## 8. Design Rules Going Forward

These rules should be treated as policy:

1. No diagnostic may alter `need_decomp`, quadrature order, or penalty path.
2. No diagnostic may call `Solve()` or `ComputeTraction()` from inside the live step loop.
3. No diagnostic may print directly from inside element/face hot loops.
4. Diagnostics must be accepted-state based unless explicitly marked offline.
5. Exact-face forensics belong in a standalone tool, not in the production BP5 job.

---

## 9. Immediate Recommendation

Before any more MFEM-vs-Tandem comparison work:

1. strip the v58 live diagnostics out of the production path
2. rerun the clean BP5 baseline
3. create the standalone/offline diagnostic extractor before reintroducing exact-face comparison

That order matters.

If we keep adding new forensic hooks into the live operator path, we will continue to debug interactions between the solver and the diagnostics instead of debugging the BP5 physics mismatch itself.

---

## 10. Bottom Line

The v58 work established an important engineering lesson:

> diagnostics for a tightly coupled nonlinear MPI solve must be treated like tooling, not like optional branches inside the physics kernel.

The correct long-term solution is:

- production-safe accepted-step observers for lightweight monitoring
- standalone offline extractors for heavy exact-face forensic comparison

That split is the safest way to avoid repeating the `dd73ac4` failure mode.
