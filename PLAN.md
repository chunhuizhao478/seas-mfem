# Implementation Plan: Adaptive ParaView Output for BP5 1000m 1800-year Production Run

## Overview

Extend the existing V_max-adaptive ParaView output in `paraview_output.hpp` so the
BP5 1000m Frontera production run produces a single on-fault PVD time series that
is **sparse during the ~225-year interseismic periods** and **dense for the ~100
seconds of each coseismic event**. The current schedule is already V_max-driven,
but (a) the production sbatch does not enable ParaView at all, (b) the three
intervals are hardcoded constants (0.01 s / 1 s / 1 yr) unsuitable for an
1800-year / 400-rank run, and (c) there is no hysteresis, no fault-only mode,
and no way to tune any of it from the command line.

The plan keeps the existing `ParaViewOutput` API, its `Save`/`ShouldWrite`/
`PeekShouldWrite`/`CommitSchedule` protocol, and the `WriteFaultSurfaceVTU`
triangle-PVD output path untouched in spirit — it only generalizes
`OutputInterval(V_max)` into a runtime-configurable schedule, adds hysteresis,
plus a fault-surface-only switch, and wires all of it through CLI flags into
`bp5_verification_full.cpp` and the production sbatch.

## Constraints

### Interface constraints (must not break)
- `ParaViewOutput<MeshType>` public API already consumed by `bp5_verification_full.cpp:1463–1638`, `drivers/tpv102_driver.cpp`, `tests/unit/test_io.cpp`, and `bp2_serial_smoke.cpp`. Do not remove or rename:
  - `Save(int, real_t, real_t) → bool`
  - `ShouldWrite(int, real_t, real_t) → bool`
  - `PeekShouldWrite(int, real_t, real_t) const → bool`
  - `CommitSchedule(real_t) → void` — **kept as a back-compat shim**; a new two-argument overload `CommitSchedule(real_t time, real_t V_max)` is added (see §1 step 3) so fault-only callers can advance `current_regime_` in lockstep.  The one-arg shim delegates to the two-arg form with the last seen V_max (stored in a private `last_v_max_` member updated by `Peek`/`Save`/`ShouldWrite`).
  - `ForceSave(int, real_t) → void`
  - `InitFaultOutputBP5`, `SetFaultParamsBP5`, `UpdateFaultFieldsBP5`, `WriteFaultSurfaceVTU`
  - public members `output_every_n_steps`, `fixed_dt`
- Keep `static real_t OutputInterval(real_t V_max)` callable (test_io.cpp:467–475 uses it) — it becomes a thin wrapper around a default-constructed schedule.
- Default constructor behavior must match today exactly (coseismic 0.01 s, nucleation 1 s, interseismic 1 yr, thresholds 1e-3 / 1e-6), so the BP2 smoke test, combined tests, and bp5 smoke `--paraview-dt` runs remain byte-compatible.

### Dependency constraints
- `paraview_output.hpp` sits under `miniapps/seas/io/`. It already depends only on `mfem.hpp` and `config/bp5_params.hpp` (for `seconds_per_year`). Do not introduce new headers beyond `<algorithm>` / `<optional>` if needed.
- Must compile in both serial and parallel (`#ifdef MFEM_USE_MPI`) branches — class is templated on `MeshType ∈ {Mesh, ParMesh}`.
- Must compile on Frontera with `intel/19.1.1`, `impi/19.0.9`, C++17 (per working sbatch module list, feedback memory).

### Convention constraints
- Follow the `bp5_verification_full.cpp` CLI parsing style: long flags with `--` prefix, accept value via the next argv element (e.g. `--tandem-dt-init 0.01`), collected in a flat `for (int i = 1; i < argc; i++)` loop around lines 733–…
- Never mutate BP5 shared source the way the TPV102 fault-utility rule forbids (memory: feedback_tpv102_bp5_no_shared_edit) — this plan touches `paraview_output.hpp` which is shared. Keep changes **strictly additive and backward-compatible**; no renames, no signature changes to existing methods.
- Never run the 1000m BP5 mesh locally (memory: feedback_no_local_reproducer). All verification runs go through Frontera sbatch; local tests operate on tiny meshes or synthesized (step, time, V_max) sequences only.
- Before any Frontera submission, notify the user and wait for explicit approval (memory: feedback_frontera_approval).

### Numerical / correctness constraints
- The schedule must be **monotonically responsive**: once the decision to write at `(cycle, time, V_max)` is taken, `last_write_time_` advances to `time`. This property already holds via `CommitSchedule` and must be preserved.
- Hysteresis thresholds must satisfy `V_enter_co > V_exit_co >= V_enter_nu > V_exit_nu > 0` to be well-defined.
- No schedule change may silently drop the required initial-condition write at `t=0` (currently done via `paraview_write(0, 0.0, V_init)` at bp5_verification_full.cpp:1638).

---

## Phase 1: Parameterize `OutputInterval` with an `AdaptiveSchedule` struct

### Goal
`ParaViewOutput` carries a runtime-mutable `AdaptiveSchedule` whose three
interval values and two V thresholds can be overridden from the caller, while
the legacy static `OutputInterval(V_max)` still returns today's values.

### Files to Modify
- `miniapps/seas/io/paraview_output.hpp` — add `AdaptiveSchedule` struct, a `adaptive_` member, and route `Save()` / `ShouldWrite()` / `PeekShouldWrite()` through it. Keep all existing public methods and the existing default behavior.

### Files to Create
- (none)

### Detailed Requirements

1. **Introduce a nested struct inside `ParaViewOutput`:**
   ```cpp
   struct AdaptiveSchedule {
      // V thresholds (m/s). V_max >= v_coseismic => coseismic regime.
      real_t v_coseismic   = 1e-3;
      real_t v_nucleation  = 1e-6;
      // Hysteresis margins (factor, multiplicative). Leaving a regime requires
      // V_max to drop below (threshold / hysteresis_factor). 1.0 = no hysteresis.
      real_t hysteresis_factor = 1.0;
      // Output intervals (seconds) in each regime.
      real_t dt_coseismic  = 0.01;
      real_t dt_nucleation = 1.0;
      real_t dt_interseismic = 1.0 * BP5Params::seconds_per_year;
      // Compute the output interval for a given V_max, given the current
      // regime (for hysteresis). regime: 0=inter, 1=nucl, 2=coseismic.
      real_t Interval(real_t V_max, int current_regime) const;
      // Compute the new regime given V_max and the previous regime.
      int NextRegime(real_t V_max, int prev_regime) const;
   };
   ```
   - `Interval` returns `dt_coseismic`/`dt_nucleation`/`dt_interseismic` by regime index.
   - `NextRegime(V, r_prev)`:
     - Coseismic enter when `V >= v_coseismic`; exit when `V < v_coseismic / hysteresis_factor`.
     - Nucleation enter when `V >= v_nucleation`; exit when `V < v_nucleation / hysteresis_factor`.
     - Default `hysteresis_factor = 1.0` reproduces the step function exactly (entry threshold == exit threshold).

2. **Add mutable members to the class:**
   ```cpp
   AdaptiveSchedule adaptive_;
   int current_regime_ = 0;   // 0=interseismic, 1=nucleation, 2=coseismic
   ```
   Add a public accessor `AdaptiveSchedule& GetSchedule() { return adaptive_; }` so the driver can tune thresholds/intervals before the time loop. Initial regime 0 is correct for BP5 (V_init ≈ 1e-9).

3. **Route schedule lookups through the struct, and keep `current_regime_` in lockstep with `last_write_time_` on every committed write — including the fault-only path that uses `PeekShouldWrite` + `CommitSchedule` instead of `Save`.**

   Add a private member `real_t last_v_max_ = 0.0;` (used by the back-compat single-arg `CommitSchedule(time)` shim to remember the most recent V_max seen by `Peek`/`Save`/`ShouldWrite`).

   - `Save(cycle, time, V_max)`:
     - Preserve the `output_every_n_steps > 0` short-circuit as today (paraview_output.hpp:631–639).
     - Preserve the `fixed_dt > 0.0` override (line 641).
     - Adaptive branch: compute `new_regime = adaptive_.NextRegime(V_max, current_regime_)`; update `current_regime_ = new_regime`; then `dt_out = adaptive_.Interval(V_max, new_regime)` instead of calling the static helper.
     - Always update `last_v_max_ = V_max` (used by the one-arg `CommitSchedule` shim).
   - `ShouldWrite(cycle, time, V_max)`: mutates both `last_write_time_` AND `current_regime_`; also updates `last_v_max_`.
   - `PeekShouldWrite(cycle, time, V_max) const`: remains `const`. It must NOT mutate `current_regime_`, `last_write_time_`, or `last_v_max_`. Internally it computes `new_regime` into a **local** variable, uses it to select `dt_out`, and returns the boolean. (Callers that take the write path must subsequently call the new two-arg `CommitSchedule(time, V_max)` to advance both state variables.)
   - `CommitSchedule(real_t time, real_t V_max)` (NEW two-arg overload): advances BOTH state variables in one call:
     ```cpp
     void CommitSchedule(real_t time, real_t V_max) {
        current_regime_ = adaptive_.NextRegime(V_max, current_regime_);
        last_write_time_ = time;
        last_v_max_ = V_max;
     }
     ```
     This is the correct way for the fault-only path to advance the hysteresis state machine when it elects not to call `Save`.  Failing to do so is exactly R-P07: `current_regime_` stays pinned at 0 forever and the Phase 3 hysteresis silently no-ops.
   - `CommitSchedule(real_t time)` (back-compat shim — keeps TPV102 driver and existing call sites unchanged): delegates to the two-arg form using the private `last_v_max_` captured by the most recent `Peek`/`Save`/`ShouldWrite`:
     ```cpp
     void CommitSchedule(real_t time) { CommitSchedule(time, last_v_max_); }
     ```
     This shim is safe for any call site that invokes `Peek`/`Save`/`ShouldWrite` immediately before `CommitSchedule` (which is the documented contract at paraview_output.hpp:690–699) — `last_v_max_` is always current.

4. **Keep the legacy static helper:**
   ```cpp
   static real_t OutputInterval(real_t V_max) {
      AdaptiveSchedule s;               // defaults = legacy behavior
      return s.Interval(V_max, s.NextRegime(V_max, 0));
   }
   ```
   For the 3 test points in `test_io.cpp:467–475` this still returns `1.0*seconds_per_year`, `1.0`, and `0.01` respectively (V=1e-10 → regime 0, V=1e-4 → regime 1, V=1e-2 → regime 2).

5. **Invariants enforced via MFEM_VERIFY in a one-shot `AdaptiveSchedule::Validate()` method:**
   ```cpp
   void Validate() const {
      MFEM_VERIFY(hysteresis_factor >= 1.0,
                  "AdaptiveSchedule: hysteresis_factor must be >= 1.0, got "
                  << hysteresis_factor);
      MFEM_VERIFY(v_coseismic > v_nucleation && v_nucleation > 0.0,
                  "AdaptiveSchedule: require v_coseismic (" << v_coseismic
                  << ") > v_nucleation (" << v_nucleation << ") > 0");
      MFEM_VERIFY(dt_coseismic > 0.0 && dt_nucleation > 0.0 &&
                  dt_interseismic > 0.0,
                  "AdaptiveSchedule: all dt_* must be positive");
   }
   ```
   Call `sched.Validate()` exactly once — in `bp5_verification_full.cpp` immediately after the Phase 2 CLI-override block (see Phase 2 step 3). Do NOT call `Validate()` from inside `NextRegime` or `Interval` — those are hot-path (invoked every step); invariants are properties of the struct, not of the input, so one-shot validation at configuration time is both more efficient and semantically cleaner (a user passing `--paraview-hyst 0.5` sees the abort at CLI parse time, not at the first time-step).

### Interfaces
- New public nested type `ParaViewOutput<MeshType>::AdaptiveSchedule` with a one-shot `Validate() const` method.
- New public accessor `AdaptiveSchedule& GetSchedule()`.
- New public two-arg overload `void CommitSchedule(real_t time, real_t V_max)` — advances `last_write_time_` AND `current_regime_` AND `last_v_max_` atomically. The existing one-arg `CommitSchedule(real_t time)` is kept as a back-compat shim that delegates to the two-arg form using the last captured V_max.
- All other existing signatures unchanged.

### Edge Cases to Handle
- `V_max == NaN` / `Inf`: treat as coseismic (regime 2) so we **do not miss** output during a blow-up. Current code writes too, so behavior is preserved.
- `V_max <= 0` (e.g. pre-initialization): treat as regime 0 (interseismic). Matches today.
- `hysteresis_factor = 1.0` (default): the enter and exit thresholds coincide — `NextRegime` must return exactly the same regime as the thresholds-only lookup `OutputInterval` did before.
- Monotonic `last_write_time_`: if the user calls `PeekShouldWrite(cycle=C, t)` twice at the same `(C, t)`, both must return the same boolean (no state mutation in Peek).
- Writing at `t=0` when `last_write_time_ = -1e30` triggers a write for any schedule, which is the desired behavior for the initial condition.

### Acceptance Criteria
- [ ] `tests/unit/test_io.cpp:TestParaViewOutputInterval` still passes unchanged (no signature break).
- [ ] New unit test `TestParaViewAdaptiveScheduleHysteresis` (see Phase 6) shows that with `hysteresis_factor = 10`, `V` oscillating between 0.5e-6 and 2e-6 does **not** flip regimes on every call.
- [ ] `bp5_verification_full.cpp` still compiles and runs the BP5 smoke test with `--paraview` and no new flags, producing identical PVD frames (byte-compare the .pvd entry list) to the pre-change build.
- [ ] Built binary size change < 1 %.

### Dependencies
- Depends on: nothing.
- Required by: Phase 2, 3, 4, 6.

---

## Phase 2: Add CLI flags to `bp5_verification_full.cpp`

### Goal
The operator can tune thresholds, intervals, hysteresis, and enable a
fault-only mode from the command line, without recompiling.

### Files to Modify
- `miniapps/seas/tests/verification/bp5_verification_full.cpp`

### Detailed Requirements

1. **Add 7 new locals near line 722** (next to `use_paraview`, `paraview_step_interval`, `paraview_dt`):
   ```cpp
   real_t pv_dt_co   = -1.0;   // coseismic interval override (s)
   real_t pv_dt_nu   = -1.0;   // nucleation interval override (s)
   real_t pv_dt_inter = -1.0;  // interseismic interval override (s or yr; see flag)
   real_t pv_v_co    = -1.0;   // coseismic V threshold override (m/s)
   real_t pv_v_nu    = -1.0;   // nucleation V threshold override (m/s)
   real_t pv_hyst    = -1.0;   // hysteresis factor override
   bool   pv_fault_only = false; // if true, skip volume PVD Save()
   ```
   All negative defaults mean "inherit schedule default from Phase 1".

2. **Add flag parsing immediately after the existing `--paraview-dt` block (around line 841):**
   - `--paraview-adaptive` : explicitly enable adaptive schedule (equivalent to `--paraview` with neither `--paraview-every` nor `--paraview-dt`). This flag is a no-op in terms of defaults but makes the sbatch self-documenting.
   - `--paraview-dt-co <seconds>` → `pv_dt_co`
   - `--paraview-dt-nu <seconds>` → `pv_dt_nu`
   - `--paraview-dt-inter-yr <years>` → `pv_dt_inter = value * BP5Params::seconds_per_year`
   - `--paraview-v-co <m/s>` → `pv_v_co`
   - `--paraview-v-nu <m/s>` → `pv_v_nu`
   - `--paraview-hyst <factor>` → `pv_hyst` (must be >= 1)
   - `--paraview-fault-only` → `pv_fault_only = true; use_paraview = true;`
   Follow the existing copy-paste style: short `if (arg == "...") { ... }` chain with manual `++i`.

3. **Apply overrides INSIDE the `if (use_paraview)` guard at bp5_verification_full.cpp:1467–1519.**
   Specific insertion point: AFTER the existing `paraview_step_interval` / `paraview_dt` assignments at `:1512-1519` and BEFORE the log line at `:1521`. Do NOT place this block at `:1520` or later outside the `if` guard — `pv_out` is `nullptr` when `use_paraview` is `false` (declared as `std::unique_ptr` at `:1463`, only constructed inside the guard), so `pv_out->GetSchedule()` outside the guard segfaults.

   ```cpp
   // (inside `if (use_paraview) { ... }`, after :1519 assignments, before :1521 log)
   auto &sched = pv_out->GetSchedule();
   if (pv_dt_co    > 0) sched.dt_coseismic      = pv_dt_co;
   if (pv_dt_nu    > 0) sched.dt_nucleation     = pv_dt_nu;
   if (pv_dt_inter > 0) sched.dt_interseismic   = pv_dt_inter;
   if (pv_v_co     > 0) sched.v_coseismic       = pv_v_co;
   if (pv_v_nu     > 0) sched.v_nucleation      = pv_v_nu;
   if (pv_hyst     > 0) sched.hysteresis_factor = pv_hyst;
   sched.Validate();  // one-shot invariants (R-P04): abort at CLI time on bad inputs.
   ```
   This block must execute **before** the first `paraview_write(0, 0.0, V_init)` call at line 1638 so the initial schedule (and its validated state) is already configured when the t=0 IC write fires.

4. **Update the "ParaView output: ON (adaptive schedule)" log line (1521–1538)** to echo the active thresholds/intervals when adaptive is chosen, e.g.
   `ParaView adaptive: V_co=1e-3 V_nu=1e-6 hyst=10 dt_co=0.5s dt_nu=5s dt_inter=5yr`.

5. **Implement `--paraview-fault-only` inside the `paraview_write` lambda (1580–1633):**

   The fault-only path must (a) populate every `pv_local_*` vector consumed by `WriteFaultSurfaceVTU` — otherwise the VTU holds uninitialized or stale memory (R-P01) — and (b) advance BOTH `last_write_time_` AND `current_regime_` via the new two-arg `CommitSchedule(time, V_max)`, so the Phase 3 hysteresis state machine keeps working (R-P07). What we skip is only `pv_out->UpdateFaultFieldsBP5(...)` (updates volume-PVD GridFunctions) and `pv_out->Save(...)` (writes the volume PVD); the owned→local expansions are **required** for the fault surface VTU and must stay.

   ```cpp
   if (pv_fault_only) {
      if (!pv_out->PeekShouldWrite(step_num, time, V_max)) return;

      // --- Populate pv_local_* (same 5 calls as the volume path at
      //     bp5_verification_full.cpp:1583-1612, MINUS the UpdateFaultFieldsBP5
      //     call which updates volume-PVD GridFunctions we are not writing). ---
      Vector owned_slip;
      fault_op.GetSlip(state, owned_slip);
      domain.ExpandOwnedToLocalFault(owned_slip, pv_local_slip, 2);
      domain.ExpandOwnedToLocalFault(fault_op.GetSlipRate(),
                                     pv_local_slip_rate, 2);
      domain.ExpandOwnedToLocalFault(seas_op.GetTraction(),
                                     pv_local_traction, 2);
      {
         const int spn = 3;  // BP5 ODE layout [slip_dip, slip_strike, psi]
         const int n_owned = fault_op.NumNodes();
         Vector owned_psi(n_owned);
         for (int i = 0; i < n_owned; i++) { owned_psi(i) = state(i * spn + 2); }
         domain.ExpandOwnedToLocalFault(owned_psi, pv_local_state, 1);
      }
      if (seas_op.ElasticSigmaNEnabled() &&
          seas_op.GetNormalTraction().Size() > 0)
      {
         domain.ExpandOwnedToLocalFault(seas_op.GetNormalTraction(),
                                        pv_local_normal_stress, 1);
      }
      else
      {
         pv_local_normal_stress = 0.0;
      }
      // Skip pv_out->UpdateFaultFieldsBP5(...)  (volume-PVD GF update)
      // Skip pv_out->Save(...)                  (volume PVD write)

      // Advance BOTH last_write_time_ AND current_regime_ in lockstep (R-P07).
      pv_out->CommitSchedule(time, V_max);

      // Static friction parameters + coords for the fault surface VTU.
      Vector local_a, local_Dc, local_x2, local_x3;
      domain.ExpandOwnedToLocalFault(fault_geom.GetAValues(),  local_a,  1);
      domain.ExpandOwnedToLocalFault(fault_geom.GetDcValues(), local_Dc, 1);
      domain.ExpandOwnedToLocalFault(fault_geom.GetCoordsX2(), local_x2, 1);
      domain.ExpandOwnedToLocalFault(fault_geom.GetCoordsX3(), local_x3, 1);

      pv_out->WriteFaultSurfaceVTU(output_dir, step_num, time,
                                   mpi.Rank(), mpi.Size(),
                                   pv_local_slip, pv_local_slip_rate,
                                   pv_local_traction, pv_local_state,
                                   pv_local_normal_stress,
                                   local_a, local_Dc, local_x2, local_x3);
      return;
   }
   // else: existing path (Save + WriteFaultSurfaceVTU).
   ```
   Rationale: in fault-only mode we never call `Save()` — so no volume PVD, no displacement grid-function write, no per-rank `.vtu` for the volume mesh. This cuts disk by roughly **10–100× for the 1000m / 400-rank case** (fault triangles ≪ volume cells). What we DO keep is the fault surface VTU (which needs the 5 populated vectors) and the two-arg `CommitSchedule` (which keeps hysteresis working).

### Edge Cases to Handle
- A flag given without a value (e.g. `--paraview-dt-co` at end of argv) — follow the existing style: the `i + 1 < argc` guard prevents overrun, silently ignore.
- Caller passes `--paraview-every 100` AND `--paraview-dt-co 0.5` — document precedence explicitly in the help/log: `output_every_n_steps` > `fixed_dt` > adaptive schedule. (This is already how `Save()` works at paraview_output.hpp:631–647.)
- `--paraview-fault-only` without `--paraview`: implicitly set `use_paraview = true` (same pattern as `--paraview-every` / `--paraview-dt`).
- `--paraview-hyst 0.5` (factor < 1): the `sched.Validate()` call in step 3 aborts at CLI time with `"hysteresis_factor must be >= 1.0"` — validation is one-shot at configuration time (R-P04), NOT inside the hot-path `NextRegime`.

### Acceptance Criteria
- [ ] BP5 smoke test with no new flags: identical PVD frame list (diff of `.pvd` file) versus pre-change build.
- [ ] `seas_bp5_full --inline-mesh --max-steps 20 --paraview-fault-only` finishes without writing any volume `.vtu` file (check: `find ... -name '*.vtu' | grep -v FaultSurface | wc -l` == 0).
- [ ] `--paraview-dt-co 0.5 --paraview-dt-inter-yr 10` run on inline mesh produces the correct regime transitions (verify by parsing timestamps from the .pvd: no two frames within 10 yr interseismic interval and, during the synthesized V_max spike, no two frames within 0.5 s intervals).
- [ ] **R-P04 startup-abort:** `seas_bp5_full --inline-mesh --paraview --paraview-hyst 0.5` aborts with a recognizable `MFEM_VERIFY` message (`"hysteresis_factor must be >= 1.0"`) BEFORE the first time step — validation runs once at the `sched.Validate()` call inserted in step 3, not inside the hot-path `NextRegime`.
- [ ] **R-P02 null-safety:** `seas_bp5_full --inline-mesh --max-steps 5` (no `--paraview` at all) completes without segfault; the CLI override block only touches `pv_out` when the `if (use_paraview)` guard is taken.
- [ ] Usage message (when no args or invalid flag) lists the new flags.

### Dependencies
- Depends on: Phase 1.
- Required by: Phase 5 (sbatch uses the flags).

---

## Phase 3: Add V_max hysteresis to prevent cadence chattering

### Goal
With the default `hysteresis_factor = 1.0`, behavior is identical to today.
When the user sets `--paraview-hyst 10`, the schedule no longer flips regimes
on every infinitesimal V_max oscillation near a threshold.

### Files to Modify
- `miniapps/seas/io/paraview_output.hpp` — implementation of `AdaptiveSchedule::NextRegime` from Phase 1.

### Detailed Requirements

1. Implement `NextRegime` as a state machine.  **Entry comparisons are strict `>` (not `>=`)** to preserve byte-compatibility with the legacy helper at paraview_output.hpp:710–721 which uses `V_max > 1e-3` / `V_max > 1e-6` (R-P08).  Exit comparisons remain strict `<` (an exact-threshold V stays in its current regime, matching legacy boundary behavior).
   ```cpp
   int NextRegime(real_t V, int prev) const {
      const real_t V_co_enter = v_coseismic;
      const real_t V_co_exit  = v_coseismic  / hysteresis_factor;
      const real_t V_nu_enter = v_nucleation;
      const real_t V_nu_exit  = v_nucleation / hysteresis_factor;

      // Treat non-finite V as "coseismic" so we keep dense output during blowup.
      if (!std::isfinite(V)) return 2;

      switch (prev) {
         case 0: // interseismic
            if (V >  V_co_enter) return 2;   // strict > : matches legacy V_max > 1e-3
            if (V >  V_nu_enter) return 1;   // strict > : matches legacy V_max > 1e-6
            return 0;
         case 1: // nucleation
            if (V >  V_co_enter) return 2;   // strict > : matches legacy
            if (V <  V_nu_exit)  return 0;
            return 1;
         case 2: // coseismic
            if (V <  V_co_exit) {
               // Drop-past-threshold semantics (intentional; see R-P03):
               // we treat the coseismic→next transition as if entering a fresh
               // regime.  The target uses the STRICT ENTRY threshold (V_nu_enter,
               // not V_nu_exit), so a V landing in [V_nu_exit, V_nu_enter] goes
               // directly to regime 0 (interseismic) — it does NOT transfer any
               // accumulated nucleation-hysteresis credit from the outbound path.
               // Rationale: after a completed rupture, V is usually already deep
               // in interseismic territory; returning to nucleation only when V
               // is unambiguously above the entry threshold avoids over-sampling
               // the post-rupture ring-down.  The test in Phase 6 locks this
               // semantic in — a V_max trajectory {2e-3, 5e-5, 5e-7} from a
               // coseismic regime yields {2, 2, 0} (NOT {2, 1, 0}).
               if (V >  V_nu_enter) return 1;
               return 0;
            }
            return 2;
         default: return 0;
      }
   }
   ```

2. `Interval(V, regime)` maps regime→interval directly, ignoring V (so the cadence is stable inside a regime):
   ```cpp
   real_t Interval(real_t, int regime) const {
      switch (regime) {
         case 2:  return dt_coseismic;
         case 1:  return dt_nucleation;
         default: return dt_interseismic;
      }
   }
   ```

### Edge Cases
- **Exact-threshold V (R-P08 byte-compat):** `hysteresis_factor = 1` ⇒ entry and exit thresholds coincide.  Entry uses strict `>`, exit uses strict `<`, so V exactly equal to `v_coseismic` stays in whatever regime it is currently in (e.g., an interseismic V that rises to exactly 1e-3 does NOT yet enter coseismic; it must exceed 1e-3).  This matches the legacy `V_max > 1e-3` / `V_max > 1e-6` semantics at paraview_output.hpp:710–721.
- **Coseismic→drop skips nucleation window (R-P03, intentional):** a V_max trajectory that drops from coseismic through the narrow `[V_nu_exit, V_nu_enter]` band within one sample interval transitions directly to regime 0, not through regime 1 — see the comment in `NextRegime` case 2 above.  This is a deliberate design choice to avoid over-sampling post-rupture ring-down.
- Very large `hysteresis_factor` (e.g. 1e6): once coseismic, we stay coseismic essentially forever. This is acceptable — tune-only knob, user's responsibility.

### Acceptance Criteria
- [ ] Unit test: with `hysteresis_factor = 10`, `V_max` sequence `1e-7, 1.1e-6, 1e-7, 1.1e-6, ...` (oscillating around nucleation threshold) yields regime transitions `0 → 1 → 1 → 1 → 1 → ...` (no flipping back to 0 until V falls below `1e-6 / 10 = 1e-7`).
- [ ] Unit test: with `hysteresis_factor = 1` (default), same sequence yields alternating `0 → 1 → 0 → 1 → ...` (legacy behavior).
- [ ] **Boundary test (R-P08):** `ParaViewOutput<Mesh>::OutputInterval(1e-3)` returns `1.0` (nucleation), not `0.01` (coseismic); `OutputInterval(1e-6)` returns `1.0 * BP5Params::seconds_per_year` (interseismic), not `1.0` (nucleation). These are the exact-threshold points at which the legacy static helper uses strict `>` — byte-compat is maintained.
- [ ] **Coseismic-drop test (R-P03):** with `hysteresis_factor = 10`, starting from regime 2, the V_max sequence `{2e-3, 5e-5, 5e-7}` produces regimes `{2, 2, 0}`. The middle sample has `V = 5e-5 > V_co_exit = 1e-4` → false, so it stays in 2. The third sample has `V = 5e-7 < V_co_exit = 1e-4` → drop; then `V = 5e-7 < V_nu_enter = 1e-6` → regime 0. The intermediate nucleation regime is intentionally skipped on the drop path.

### Dependencies
- Depends on: Phase 1.
- Required by: Phase 6 (tests), Phase 5 (sbatch sets `--paraview-hyst 10`).

---

## Phase 4: Document the initial-condition write

### Goal
Ensure `paraview_write(0, 0.0, V_init)` still runs in fault-only mode and with
hysteresis, and is clearly documented so future changes don't silently remove
the t=0 frame.

### Files to Modify
- `miniapps/seas/tests/verification/bp5_verification_full.cpp` around line 1638.

### Detailed Requirements

1. Leave the existing call in place. Verify by inspection that in fault-only
   mode the t=0 call still produces a FaultSurface/fault_surface_c0.vtu (the
   `pv_out->PeekShouldWrite(0, 0.0, V_init)` path returns true because
   `last_write_time_ == -1e30`).
2. No code change beyond a one-line comment stating:
   `// t=0 IC write: always triggers because last_write_time_ starts at -1e30.`

### Acceptance Criteria
- [ ] `seas_bp5_full --inline-mesh --max-steps 1 --paraview-fault-only` produces `FaultSurface/fault_surface_c0.pvtu` and `fault_surface.pvd` with exactly one entry at t=0.

### Dependencies
- Depends on: Phase 2.

---

## Phase 5: Update the production sbatch

### Goal
`bp5_system_update_full_production.sbatch` runs the 1000m 1800-yr job with
adaptive fault-only PVD output tuned for inspection — long interseismic
intervals, dense during EQs, small enough on-disk to finish the 48 h job.

### Files to Modify
- `miniapps/seas/jobs/bp5/bp5_system_update_full_production.sbatch`

### Detailed Requirements

1. Append the following flags to the `ibrun ./seas_bp5_full ...` invocation
   currently at lines 66–80 (after `--verify`):
   ```
         --paraview \
         --paraview-fault-only \
         --paraview-adaptive \
         --paraview-v-co        1e-3 \
         --paraview-v-nu        1e-6 \
         --paraview-hyst        10.0 \
         --paraview-dt-co       0.5 \
         --paraview-dt-nu       10.0 \
         --paraview-dt-inter-yr 5.0
   ```
   Justification of values for a 1800-yr / ~8-event / ~100-s rupture run:
   - `dt-inter = 5 yr` → ~360 interseismic frames (was 1800). Still shows slow-stress build-up.
   - `dt-nu = 10 s` → a few hundred frames across the minutes-to-hours nucleation phase of each event.
   - `dt-co = 0.5 s` → ~200 frames per ~100-s rupture (was ~10,000). Resolves rupture propagation across the 40×80 km fault comfortably.
   - `hyst = 10` → thresholds enter at 1e-3 / 1e-6 but exit at 1e-4 / 1e-7 → no flip-flop near the boundaries.  **This knob is only effective because the fault-only path uses the new two-arg `CommitSchedule(time, V_max)` (Phase 1 §1.3 + Phase 2 §2.5) to advance `current_regime_` in lockstep with `last_write_time_` — without that, R-P07 would pin `current_regime_` at 0 forever and hyst would silently no-op.**

2. Add a short comment block right above the ibrun line explaining the flag choices and expected disk footprint (back-of-envelope: ~8 events × ~200 frames + nucleation + interseismic ≈ ~3000 frames; with ~400 fault faces per rank × 3 vertices × ~12 fields × ~400 ranks ≈ a few GB total — versus hundreds of GB had volume PVD been enabled).

3. Do **not** change any existing flags (mesh, order, solver, petsc-ts, checkpoint interval, etc.) — this phase is additive only.

### Edge Cases
- Some ranks have zero fault faces (common with 400-rank partitioning of a fault that only occupies a strip). They still call `WriteFaultSurfaceVTU` and will emit an empty VTU; rank 0's PVTU still references them. This is OK today and remains OK — see paraview_output.hpp:423–620.
- If the cluster aborts mid-run, the PVD index (`fault_surface.pvd`) is rewritten in full each write (paraview_output.hpp:779–794) so the on-disk file always reflects the last successful write.

### Acceptance Criteria
- [ ] Dry-run locally (tiny inline mesh with `--paraview-fault-only` + the same non-path flags) completes without error.
- [ ] A short Frontera sbatch (not this production one — only after user approval per feedback memory) confirms the flags parse correctly and `FaultSurface/fault_surface.pvd` is produced.

### Dependencies
- Depends on: Phases 1–3.

---

## Phase 6: Tests

### Goal
Lock the new schedule/hysteresis/fault-only behavior behind unit tests so
future edits to `paraview_output.hpp` cannot silently regress the Frontera
run's output cadence.

### Files to Modify
- `miniapps/seas/tests/unit/test_io.cpp`

### Detailed Requirements

1. **Add `TestParaViewAdaptiveScheduleHysteresis`:**
   ```cpp
   void TestParaViewAdaptiveScheduleHysteresis() {
      using Schedule = ParaViewOutput<Mesh>::AdaptiveSchedule;
      Schedule s;
      s.v_coseismic     = 1e-3;
      s.v_nucleation    = 1e-6;
      s.hysteresis_factor = 10.0;
      int r = 0;

      // Enter coseismic
      r = s.NextRegime(2e-3, r); TEST_EQ(r, 2);
      // Stay coseismic down to exit threshold
      r = s.NextRegime(2e-4, r); TEST_EQ(r, 2);
      // Drop below exit → nucleation
      r = s.NextRegime(5e-5, r); TEST_EQ(r, 1);
      // Stay nucleation across thresholds
      r = s.NextRegime(1.1e-6, r); TEST_EQ(r, 1);
      // Drop below nu_exit
      r = s.NextRegime(5e-8, r);  TEST_EQ(r, 0);
   }
   ```
2. **Add `TestParaViewAdaptiveScheduleCustomIntervals`:**
   - Construct `ParaViewOutput<Mesh>` on a tiny inline mesh.
   - Mutate `sched.dt_interseismic = 60.0; sched.dt_coseismic = 0.1;`.
   - Call `Save(0, 0.0, 1e-10)` → must return true.
   - Call `Save(1, 30.0, 1e-10)` → must return false (30 s < 60 s).
   - Call `Save(2, 120.0, 1e-10)` → must return true.
   - Call `Save(3, 120.05, 2e-3)` (coseismic spike) → must return false (0.05 < 0.1).
   - Call `Save(4, 120.15, 2e-3)` → must return true.
3. **Add `TestParaViewLegacyInterval`:** verify `ParaViewOutput<Mesh>::OutputInterval(V)` still returns `1.0*yr / 1.0 / 0.01` for `V = 1e-10 / 1e-4 / 1e-2` — this is the exact same three assertions already in `TestParaViewOutputInterval` and can stay in that function unchanged. **While touching test_io.cpp, also fix the R-P05 inconsistency in the existing `TestParaViewOutputInterval`:** the interseismic assertion at test_io.cpp:468 currently compares against `BP2Params::seconds_per_year`, but the code-under-test uses `BP5Params::seconds_per_year`. Both constants are numerically identical (`365.25*24*3600.0`) today, but the inconsistency is a latent trap if one definition ever drifts. Change the reference to `BP5Params::seconds_per_year` to match the production source of truth.
4. **Add `TestParaViewFaultOnlyAdvancesRegime` (R-P07 coverage):** This test covers the gap Phase 6 step 1 doesn't — it exercises the `PeekShouldWrite` + `CommitSchedule(time, V_max)` pair exactly as the fault-only driver path does, and asserts that `current_regime_` (and therefore `dt_out` from the next `PeekShouldWrite`) advances correctly. Without the R-P07 fix this test fails because the one-arg `CommitSchedule(time)` never advanced `current_regime_`, and `Peek` never mutates.
   ```cpp
   void TestParaViewFaultOnlyAdvancesRegime() {
      // Minimal inline mesh — does NOT run the 1000m BP5 mesh (feedback memory).
      Mesh mesh = Mesh::MakeCartesian3D(2, 2, 1, Element::HEXAHEDRON, 2.0, 2.0, 1.0);
      ParaViewOutput<Mesh> pv("test_pv_fault_only_regime", mesh, 1);
      auto &s = pv.GetSchedule();
      s.v_coseismic = 1e-3; s.v_nucleation = 1e-6;
      s.hysteresis_factor = 10.0;
      s.dt_coseismic = 0.5; s.dt_nucleation = 10.0;
      s.dt_interseismic = 1.0 * BP5Params::seconds_per_year;
      s.Validate();

      // Simulate the fault-only call pattern across 3 consecutive steps.
      // V trajectory {2e-3, 2e-3, 2e-4}: hysteresis_factor=10 ⇒ V_co_exit=1e-4.
      // step 2's V=2e-4 > V_co_exit, so regime MUST stay at 2 (coseismic).

      // Step 0: initial coseismic spike. last_write_time_ is -1e30, so fires.
      TEST_ASSERT(pv.PeekShouldWrite(0, 0.0, 2e-3));
      pv.CommitSchedule(0.0, 2e-3);  // advances last_write_time_ AND current_regime_

      // Step 1: dt_coseismic=0.5s has elapsed. Schedule fires.
      TEST_ASSERT(pv.PeekShouldWrite(1, 0.5, 2e-3));
      pv.CommitSchedule(0.5, 2e-3);

      // Step 2: V drops to 2e-4 (still above V_co_exit=1e-4 thanks to hyst=10).
      // With R-P07 fix (current_regime_=2): dt_out=dt_coseismic=0.5s, elapsed
      //   0.5s from last_write_time_=0.5 to time=1.0 ⇒ fires.
      // Without R-P07 fix (current_regime_ stuck at 0): NextRegime(2e-4, 0)
      //   returns 1 (nucleation); dt_out=dt_nucleation=10s, elapsed 0.5s ⇒
      //   DOES NOT fire. The test assertion below would then fail.
      TEST_ASSERT(pv.PeekShouldWrite(2, 1.0, 2e-4));
   }
   ```
5. **Add `TestParaViewLegacyBoundaryPreserved` (R-P08 coverage):** assert that `OutputInterval(V)` at the legacy thresholds returns the legacy regime's interval (not the adjacent one). Without the R-P08 `>=`→`>` fix this test fails.
   ```cpp
   void TestParaViewLegacyBoundaryPreserved() {
      // V == 1e-3 exactly: legacy uses V_max > 1e-3 (strict), so regime 1 (nu).
      TEST_NEAR(ParaViewOutput<Mesh>::OutputInterval(1e-3),
                1.0, 1e-15, "V=1e-3 (exact): nucleation, not coseismic");
      // V == 1e-6 exactly: legacy uses V_max > 1e-6, so regime 0 (interseismic).
      TEST_NEAR(ParaViewOutput<Mesh>::OutputInterval(1e-6),
                1.0 * BP5Params::seconds_per_year, 1.0,
                "V=1e-6 (exact): interseismic, not nucleation");
   }
   ```
6. **Add `TestParaViewCoseismicDropSkipsNucleation` (R-P03 coverage):** lock in the intentional "coseismic→drop uses strict entry threshold" semantic.
   ```cpp
   void TestParaViewCoseismicDropSkipsNucleation() {
      using Schedule = ParaViewOutput<Mesh>::AdaptiveSchedule;
      Schedule s;
      s.v_coseismic = 1e-3; s.v_nucleation = 1e-6;
      s.hysteresis_factor = 10.0;
      // Derived thresholds: V_co_exit=1e-4, V_nu_enter=1e-6, V_nu_exit=1e-7.
      int r = 2;                                   // start in coseismic
      r = s.NextRegime(2e-4, r); TEST_EQ(r, 2);    // V=2e-4 > V_co_exit=1e-4: stay
      r = s.NextRegime(5e-7, r); TEST_EQ(r, 0);    // drop; V=5e-7 < V_nu_enter=1e-6 → 0
      // Starting fresh in regime 2 with V landing inside [V_nu_exit, V_nu_enter]
      // = [1e-7, 1e-6] — the "intentionally skipped" window — must drop straight
      // to 0, NOT to 1.  This is the R-P03 design choice.
      r = 2;
      r = s.NextRegime(5e-7, r); TEST_EQ(r, 0);
   }
   ```
7. **Add `TestAdaptiveScheduleValidate` (R-P04 coverage):** exercises the valid-input path of `Validate()` (defaults must pass cleanly). Testing the abort-on-invalid side of `MFEM_VERIFY` requires a death-test harness which `test_io.cpp` does not currently carry; adding one is out of scope for this phase. The integration-level check is that passing `--paraview-hyst 0.5` to the driver aborts during startup (BEFORE the first time step), not during the time loop — this is verified manually when the Phase 2 CLI wiring is exercised and documented in the Phase 2 Acceptance Criteria.
   ```cpp
   void TestAdaptiveScheduleValidate() {
      using Schedule = ParaViewOutput<Mesh>::AdaptiveSchedule;
      Schedule s;                         // defaults — must pass.
      s.Validate();                       // no abort, no throw.
      s.hysteresis_factor = 10.0;         // legal user override.
      s.dt_coseismic = 0.5;
      s.dt_nucleation = 10.0;
      s.dt_interseismic = 5.0 * BP5Params::seconds_per_year;
      s.Validate();                       // still passes.
   }
   ```
8. **Add all six tests to `main()` near line 701.** Order them by complexity (Hysteresis → Boundary → CoseismicDrop → CustomIntervals → FaultOnlyAdvancesRegime → Validate).
9. **Do not add any test that runs the 1000m mesh** (feedback memory: feedback_no_local_reproducer).

### Acceptance Criteria
- [ ] `make test-io` (or the equivalent `test_io` binary) passes with all six new tests.
- [ ] `make test` (full unit test suite) still passes.
- [ ] Flipping `>` to `>=` in `NextRegime` entry branches alone causes `TestParaViewLegacyBoundaryPreserved` to fail (confirms R-P08 regression detector is live).
- [ ] Replacing the two-arg `CommitSchedule(time, V_max)` call in the fault-only path with the back-compat one-arg `CommitSchedule(time)` AND setting `last_v_max_ = 0.0` permanently causes `TestParaViewFaultOnlyAdvancesRegime` to fail (confirms R-P07 regression detector is live).

### Dependencies
- Depends on: Phases 1–3.

---

## Testing Strategy

1. **Unit level** (Phase 6): hysteresis state machine, custom interval selection,
   and legacy-behavior preservation — all on synthesized `(cycle, time, V_max)`
   tuples. Runs in < 1 s on a laptop.
2. **Integration level (laptop, tiny mesh only):** `seas_bp5_full --inline-mesh --max-steps 20 --paraview-fault-only --paraview-dt-co 0.1 --paraview-dt-nu 1 --paraview-dt-inter-yr 0.01` — verify:
   - `FaultSurface/fault_surface.pvd` exists and has at least 2 entries.
   - No volume `.vtu` (`find ... -name '*.vtu' -not -path '*FaultSurface*' | wc -l == 0`).
   - `.pvd` entries have strictly increasing timestamps.
3. **Frontera level** (only after user approval per feedback memory): submit the
   updated production sbatch. Success criteria:
   - Output produces `FaultSurface/fault_surface.pvd` with roughly the predicted ~3000 entries.
   - Time series benchmark comparison against golden still passes at `--regression-tolerance 1e-2` — adaptive PVD must not affect numerical output, only disk I/O.
4. **Local oversubscription caveat** (feedback_no_local_oversubscribe): never run >14 ranks locally for debugging; prefer 4-rank smoke for all local checks.

## Risk Assessment

1. **Silent regression of BP5 numerical output**: Adaptive PVD writes are purely
   diagnostic — they pull read-only data from `seas_op` / `fault_op`. No risk to
   the ODE state. The existing `paraview_write` lambda already makes this clear.
   Detection: the 1e-2 regression tolerance vs. golden should flag any mistake.
2. **Disk footprint larger than expected**: If rupture duration is longer than
   assumed (e.g. 200 s instead of 100 s at 0.5 s interval × 400 ranks), fault VTU
   count could grow 2–3× over estimate. Mitigation: compute the worst-case frame
   count from max_steps and log it at startup; abort if it exceeds a user-supplied
   cap (not in this plan, but worth a follow-up flag `--paraview-max-frames N`).
3. **Threshold chattering across regime boundaries without hysteresis**: Default
   `hysteresis_factor = 1.0` preserves today's edge behavior — user who leaves the
   default gets the exact legacy cadence. Setting hysteresis = 10 is explicit.
4. **PVD file rewrite cost on every write** (paraview_output.hpp:779–794): at
   ~3000 entries rewriting is still <1 ms. Not a concern at this scale.
5. **Shared-face rank mismatch in WriteFaultSurfaceVTU**: Not touched by this
   plan. Existing code already handles rank-0-writes-PVTU-only, others write VTU
   (paraview_output.hpp:582–614).
6. **Tricky area — `PeekShouldWrite` vs. `CommitSchedule` pattern**: The
   fault-only path in Phase 2 must call `PeekShouldWrite`, then populate all
   `pv_local_*` vectors (R-P01), then `CommitSchedule(time, V_max)` (two-arg
   overload — R-P07, not the one-arg shim, since the shim relies on
   `last_v_max_` being updated by a recent `Peek`/`Save` and that is already
   the case here, but passing `V_max` explicitly makes the intent obvious and
   is robust against future refactors that might change when `last_v_max_`
   gets updated). Forgetting the `Commit` call means every cycle fires a write
   until the schedule catches up AND `current_regime_` never advances — both
   silent failure modes.
7. **`ForceSave` leaks the volume PVD in fault-only mode**: The existing
   `ForceSave(cycle, time)` at paraview_output.hpp:702–707 unconditionally
   calls `ForceSaveImpl` which writes the volume PVD via `pv_.Save()`.  If
   `--paraview-fault-only` is set but a caller invokes `ForceSave` (e.g., a
   final-state write at simulation end, an event-onset checkpoint for
   golden-file matching), the volume PVD leaks out despite the flag.  This
   plan does NOT gate `ForceSave` on `pv_fault_only` — the existing call
   sites in `bp5_verification_full.cpp` use `ForceWrite` on the
   `bench_out` (BenchmarkOutput), NOT on `pv_out`, so the leak does not
   manifest for the Phase 5 sbatch.  Documented here so a future patch
   adding `pv_out->ForceSave(...)` knows to add a `!pv_fault_only` guard.
7. **Frontera module set** (feedback_sbatch_modules): Phase 5 changes do not
   alter the module stanza; the working `intel/19.1.1 + impi/19.0.9 + hypre +
   mumps + parmetis + petsc/3.15` set at sbatch lines 35–42 stays as-is.
