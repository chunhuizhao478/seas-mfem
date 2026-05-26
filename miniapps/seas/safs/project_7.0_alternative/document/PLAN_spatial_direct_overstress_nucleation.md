# Implementation Plan: `direct_overstress` nucleation kind for `seas_spatial_dyn_driver`

## Overview
Add a second, config-selectable nucleation kind — `direct_overstress` — to the
SAFS spatial dynamic-rupture driver, alongside the existing `gradual_overstress`.
Where `gradual_overstress` ramps the per-DOF target `Δτ · F(r)` into the
time-varying nucleation channel (`tau{1,2}_nuc`) over `[0, T_nuc_s]` via a
SCEC smoothStep, `direct_overstress` adds the **full** `Δτ · F(r)` **statically
to the initial pre-stress** `tau{1,2}_0` at initialization — exactly how SCEC
TPV205 nucleates (the nucleation patch's initial shear simply starts above
yield). The `tau{1,2}_nuc` channels remain zero for `direct_overstress`; the
per-substep accumulator becomes a no-op.

This is the user-confirmed mechanism ("static initial pre-stress, true TPV205"),
wired as a selectable `[nucleation].kind` so the two can be A/B'd from config
alone with no recompile.

## Background / current state (verified by code reading)

- **Consumption path** (`dynamic/tpv205_substep_iterator.cpp:88-90`,
  `dynamic/fault_face_flux.cpp:165-167`): the friction solve sees
  `tau1_total = tau1_0 + tau1_nuc + tau1_trial` (likewise strike, σ_n).
  So a static bump on `tau{1,2}_0` and a held bump on `tau{1,2}_nuc` reach the
  friction law identically; the only difference is whether the bump is present
  at `t = 0` (pre-stress) or stepped/ramped in during dynamics (nuc channel).
- **TPV205** (`dynamic/tpv205_setup.hpp` header doc): nucleation is a *static
  patch pre-stress* written into `tau2_0`; `tau{1,2}_nuc / sigma_n_nuc` are
  zero at init and **never** updated.
- **Pre-stress storage** (`fault/fault_geometry.hpp:459`,
  `fault/fault_geometry_safs_templated.inl:77-126`): `geom.GetTauPre()` returns
  a read-only `Vector` of size `2·N` interleaved `(tau1=dip_i, tau2=strike_i)`.
  The `SEAS_ZERO_DIP_PRESTRESS` env toggle zeroes the dip slot **inside**
  `ComputeSAFSParams`, so anything `GetTauPre()` returns already reflects it.
- **Init consumption** (`dynamic/spatial_setup.hpp:84-85,189-191,240`):
  `InitializeFaultDOFs_Spatial` reads `tau_pre[2·i]→tau1_0`,
  `tau_pre[2·i+1]→tau2_0`, with a hard `tau_pre.Size() == 2·ndof` check.
- **Driver wiring** (`drivers/spatial_dyn_driver.cpp`): resolves
  `nuc_params` at :1162-1167; passes `geom.GetTauPre()` to both
  `InitializeFaultDOFs_Spatial` (:1222) and the `[derived]` gate
  `PrintDerivedAndCheck` (:1302); the per-substep `nuc_cb` (:1749-1757) calls
  `ApplyGradualOverstressIncrement`; PV publishes `nuc_radial / nuc_amplitude`
  from `nuc_params` (:1507-1517); a shared-face ramp-as-step warning fires at
  :1177-1190; a shared-fault consistency tripwire fires through the
  `t ≤ T_nuc_s` window at :1975-1978.
- **The `[derived]` gate** (`dynamic/spatial_print_derived.cpp:184-313`)
  defines the nucleation region from `nuc.gradual_overstress.{center,radius}_*`
  and tests, **inside** the region, the post-nucleation trigger
  `|(tau1_pre + Δτ_dip·F, tau2_pre + Δτ_strike·F)| ≥ μ_s·σ_n_eff` by adding
  `nuc_params.amplitude_{dip,strike}(i)` to the base pre-stress; **outside**, it
  requires `|τ_pre| / (μ_s·σ_n_eff) < 1`. **Crucially this is the same question
  for both kinds** ("after the full overstress, does the patch trigger and does
  the background stay subcritical"), so the gate stays logic-identical if it is
  fed the *base* pre-stress + `nuc_params`; only its region-geometry accessors
  must become kind-aware.

## The math

Per-DOF Gaussian factor in the fault-local `(dip, strike)` frame (unchanged,
`dynamic/spatial_nucleation.cpp:43-74`):

```
dx = (x_i − c) · dip_i ,   ds = (x_i − c) · strike_i
F(r_i) = exp( −[ (dx / radius_dip)² + (ds / radius_strike)² ] )   ∈ [0,1]
```

`direct_overstress` modifies the **initial** pre-stress, in place, per DOF:

```
tau1_0'(i) = tau1_0(i) + F(r_i) · Δτ_dip
tau2_0'(i) = tau2_0(i) + F(r_i) · Δτ_strike
σ_n0       unchanged          tau{1,2}_nuc ≡ 0
```

Gate conditions (already implemented; restated for clarity), evaluated on the
*base* `tau_pre` + amplitudes `a_dip(i)=F·Δτ_dip`, `a_str(i)=F·Δτ_strike`:

```
inside  (r_i ≤ 3·max(radius_dip,radius_strike)):
        trigger  |(tau1_0+a_dip, tau2_0+a_str)| ≥ μ_s·σ_n_eff
outside (r_i  >  3·max(...)):
        subcritical  |(tau1_0, tau2_0)| / (μ_s·σ_n_eff) < 1
```

For `gradual_overstress` `tau1_0` is the un-bumped pre-stress and the bump is
applied later in `tau*_nuc`; for `direct_overstress` the bump is folded into
`tau1_0'` at init. Either way the gate verifies the same final equilibrium, so
**no gate arithmetic changes** — see Phase 3.

## Constraints

- **No recompile to switch kinds**: selection is `[nucleation].kind` in TOML
  only (mirror `gradual_overstress`).
- **Do not perturb the validated `gradual_overstress` path**: its parser branch,
  resolver math, per-substep accumulator, gate arithmetic, and existing unit
  tests (`test_spatial_friction_config.cpp` T-21/22/23, `test_spatial_nucleation`,
  `test_spatial_print_derived`) must remain byte-for-byte equivalent in behavior.
  Refactors that touch shared code must be behavior-preserving and re-verified
  against those tests.
- **TPV*/BP5 byte-exact regression contract is untouched**: this is a
  spatial-driver-only change; do not edit `tpv102/104/205_*` or `bp5_*`
  nucleation paths. `dynamic/spatial_nucleation.{hpp,cpp}` and
  `spatial/code/spatial_friction.{hpp,cpp}` are spatial-only.
- **No hard-coded constants** (project rule): all amplitudes/radii/centers come
  from the TOML spec; the `3·max(radius)` outside factor is the existing
  `PrintDerivedConfig::outside_safety_factor`, not a literal.
- **Frame layout is fixed**: `dof_basis` is `(9, N)` column-major, rows 3-5 =
  dip, rows 6-8 = strike (`spatial_nucleation.hpp:134-142`); `tau_pre` is `2·N`
  interleaved `(dip, strike)`. Re-use these exactly.
- **Sign/`SEAS_ZERO_DIP` interaction** (see Risk): the static bump is applied in
  the driver *after* `ComputeSAFSParams`, hence after `SEAS_ZERO_DIP_PRESTRESS`
  zeroes the dip slot. A non-zero `delta_tau_dip_pa` therefore re-introduces dip
  stress even under the zero-dip diagnostic. Document; default the dip bump to 0.

---

## Phase 1: Config data model + parser

### Goal
A TOML with `[nucleation] kind = "direct_overstress"` + a
`[nucleation.direct_overstress]` sub-block parses into a validated
`SpatialFrictionConfig`; `gradual_overstress` continues to parse unchanged; an
unknown kind still aborts.

### Files to Modify
- `spatial/code/spatial_friction.hpp` — extend the nucleation data model.
- `spatial/code/spatial_friction.cpp` — extend the `[nucleation]` parser branch.
- `dynamic/spatial_nucleation.hpp` — add the `DirectOverstressSpec` struct
  (lives next to `GradualOverstressSpec` since the resolver consumes it).

### Detailed Requirements
1. In `dynamic/spatial_nucleation.hpp`, add (immediately after
   `GradualOverstressSpec`, ~line 116):
   ```cpp
   /// @brief `[nucleation.direct_overstress]` TOML sub-block.  Same spatial
   /// + amplitude fields as GradualOverstressSpec but NO temporal field:
   /// the full Δτ·F(r) is folded statically into tau{1,2}_0 at init.
   /// radius_dip_m, radius_strike_m MUST be > 0 when enabled (parser-checked).
   struct DirectOverstressSpec
   {
      real_t center_x_m          = 0.0;
      real_t center_y_m          = 0.0;
      real_t center_z_m          = 0.0;
      real_t radius_dip_m        = 0.0;   ///< > 0 required when enabled
      real_t radius_strike_m     = 0.0;   ///< > 0 required when enabled
      real_t delta_tau_dip_pa    = 0.0;   ///< may be 0 (recommended for SAFS)
      real_t delta_tau_strike_pa = 0.0;   ///< may be 0
   };
   ```
2. In `spatial/code/spatial_friction.hpp`:
   - Extend the enum (line 219):
     `enum class NucleationKind { GradualOverstress, DirectOverstress };`
   - Update the doc comment at lines 212-218 (it currently says the
     "`Overstress` (TPV205-style one-shot)" kind "was removed") to describe the
     re-added static `direct_overstress`.
   - In `struct NucleationSpec` (lines 225-230) add a member:
     `DirectOverstressSpec direct_overstress;` (alongside the existing
     `GradualOverstressSpec gradual_overstress;`).
   - Add five inline kind-aware geometry accessors used by the gate (Phase 3),
     so no consumer must branch on `kind`:
     ```cpp
     real_t center_x_m()    const { return kind == NucleationKind::DirectOverstress
                                    ? direct_overstress.center_x_m   : gradual_overstress.center_x_m; }
     real_t center_y_m()    const { /* …direct_overstress.center_y_m   : gradual_overstress.center_y_m */ }
     real_t center_z_m()    const { /* …center_z_m */ }
     real_t radius_dip_m()  const { /* …radius_dip_m */ }
     real_t radius_strike_m() const { /* …radius_strike_m */ }
     ```
3. In `spatial/code/spatial_friction.cpp`, replace the single-kind branch at
   lines 756-786 with a kind switch:
   - `kind_s == "gradual_overstress"` → existing behavior verbatim (set
     `kind = GradualOverstress`, require `[nucleation.gradual_overstress]`,
     parse the 8 fields incl. `T_nuc_s`, keep the three `MFEM_VERIFY`s at
     :777-785).
   - `kind_s == "direct_overstress"` → set `kind = DirectOverstress`; require
     `nuc.contains("direct_overstress")` (abort message must name the
     `[nucleation.direct_overstress]` sub-block); parse the 7 fields into
     `cfg.nucleation.direct_overstress`; `MFEM_VERIFY(radius_dip_m > 0)` and
     `MFEM_VERIFY(radius_strike_m > 0)` with precise messages; **abort** if the
     sub-block contains a `T_nuc_s` key (it is meaningless for a static bump —
     fail loudly rather than silently ignore, per the no-silent-fallback rule).
   - any other non-empty `kind_s` → abort listing both valid kinds; empty
     `kind_s` aborts as today ("no implicit default").

### Interfaces
- `DirectOverstressSpec` (new POD, `mfem::seas::spatial`).
- `NucleationKind::DirectOverstress` (new enumerator).
- `NucleationSpec::{center_x_m,center_y_m,center_z_m,radius_dip_m,radius_strike_m}()`
  (new const accessors returning the active spec's geometry).

### Edge Cases to Handle
- `kind="direct_overstress"` but `[nucleation.direct_overstress]` missing → abort.
- `radius_*` ≤ 0 → abort with the field name.
- `T_nuc_s` present under `direct_overstress` → abort (not silently ignored).
- Both `[nucleation.gradual_overstress]` and `[nucleation.direct_overstress]`
  present → parse only the sub-block matching `kind`; the other is ignored
  (document; do not abort — keeps A/B configs editable by flipping `kind`).
- `delta_tau_dip_pa = delta_tau_strike_pa = 0` is allowed (a no-op nucleation;
  the gate will then flag "no DOF triggers" — that is informative, not a parse
  error).

### Acceptance Criteria
- [ ] A `direct_overstress` TOML round-trips into `cfg.nucleation.direct_overstress`
      with all 7 fields correct and `kind == NucleationKind::DirectOverstress`.
- [ ] A `gradual_overstress` TOML parses byte-identically to today (T-21/22/23
      pass unchanged after T-23's typo string is confirmed still-a-typo).
- [ ] Missing `[nucleation.direct_overstress]`, non-positive radius, and a
      stray `T_nuc_s` each abort with a specific message.
- [ ] `make test-spatial-friction-config` (or the bundling target) passes.

### Dependencies
- Depends on: nothing. Required by: Phases 2, 3, 4.

---

## Phase 2: Resolver — per-DOF amplitudes for the static bump

### Goal
`ResolveDirectOverstress(spec, enabled, coords, basis)` returns per-DOF
`amplitude_dip(i)=F(r_i)·Δτ_dip`, `amplitude_strike(i)=F(r_i)·Δτ_strike`,
`radial(i)=F(r_i)` — identical shape/semantics to the gradual resolver output,
with no temporal field.

### Files to Modify
- `dynamic/spatial_nucleation.hpp` / `.cpp` — add the resolver; factor the
  shared per-DOF amplitude loop so the gradual path stays behavior-identical.

### Detailed Requirements
1. Extract the per-DOF loop body shared by both resolvers into a file-local
   static helper in `spatial_nucleation.cpp` (behavior-preserving — the gradual
   path must produce identical numbers; verified by `test_spatial_nucleation`):
   ```cpp
   static GradualOverstressPerDOFParams ResolveOverstressAmplitudes(
      const real_t center[3], real_t radius_dip_m, real_t radius_strike_m,
      real_t dtau_dip_pa, real_t dtau_strike_pa,
      const Vector& dof_coords_3d, const DenseMatrix& dof_basis);
   ```
   It contains exactly the size checks (`dof_coords_3d.Size()==3N`,
   `dof_basis` shape `(9,N)`, radii > 0) and the loop at
   `spatial_nucleation.cpp:104-127`.
2. Re-point `ResolveGradualOverstress` (keep its signature + the `T_nuc_s > 0`
   `MFEM_VERIFY`) to delegate to the helper. Net behavior unchanged.
3. Add:
   ```cpp
   GradualOverstressPerDOFParams ResolveDirectOverstress(
      const DirectOverstressSpec& spec, bool enabled,
      const Vector& dof_coords_3d, const DenseMatrix& dof_basis);
   ```
   Body: `if (!enabled) return {};` then delegate to the helper with the
   direct spec's center/radii/Δτ. **No `T_nuc_s` reference.**
   (The return type is the existing `GradualOverstressPerDOFParams` — the struct
   holds kind-neutral amplitude arrays. Optionally add
   `using OverstressPerDOFParams = GradualOverstressPerDOFParams;` for clarity;
   not required.)

### Interfaces
- `ResolveDirectOverstress(...) -> GradualOverstressPerDOFParams` (new).
- `ResolveOverstressAmplitudes(...)` (new file-local static; not exported).

### Edge Cases to Handle
- `enabled == false` → three zero-sized Vectors (same as gradual).
- A rank with zero local fault DOFs (`N == 0`) → zero-sized Vectors, no loop.
- `Δτ_dip = Δτ_strike = 0` → amplitudes all 0, `radial` still the true Gaussian
  (so PV `nuc_radial` still shows the patch footprint).
- Exponent underflow far from the patch → `GaussianFactorFaceLocal` already
  returns 0 below −700 (re-used unchanged).

### Acceptance Criteria
- [ ] `ResolveDirectOverstress` on a known DOF reproduces
      `F·Δτ` to machine precision (mirror an existing
      `test_spatial_nucleation` gradual-resolver assertion).
- [ ] `ResolveGradualOverstress` outputs are unchanged vs. pre-refactor (the
      existing `test_spatial_nucleation` cases pass without edits).
- [ ] No reference to `T_nuc_s` anywhere on the direct path.

### Dependencies
- Depends on: Phase 1 (`DirectOverstressSpec`). Required by: Phase 4.

---

## Phase 3: `[derived]` gate — kind-aware region geometry

### Goal
`PrintDerivedAndCheck` evaluates the same trigger/subcritical conditions for
`direct_overstress` as for `gradual_overstress`, reading the active kind's
center/radius via the Phase-1 accessors; **no arithmetic change**.

### Files to Modify
- `dynamic/spatial_print_derived.cpp` — replace direct field reads.

### Detailed Requirements
1. Replace, at `spatial_print_derived.cpp:186-192`, the reads of
   `nuc.gradual_overstress.radius_dip_m / radius_strike_m / center_{x,y,z}_m`
   with `nuc.radius_dip_m() / nuc.radius_strike_m() / nuc.center_x_m()` etc.
   (the kind-aware accessors). All downstream math (`r_threshold`, the
   inside/outside split at :200, the trigger overshoot at :246-274 using
   `nuc_params.amplitude_*`) is unchanged.
2. The driver feeds the gate the **base** `geom.GetTauPre()` + `nuc_params`
   for **both** kinds (see Phase 4). This is why no arithmetic changes: the
   gate adds `amplitude_*` to the base pre-stress to form the post-nucleation
   total, which equals the `direct_overstress` `tau{1,2}_0'` it bakes in at
   init, and equals the gradual end-state. Add a one-line comment at the
   top of §3 (~:241) stating this invariant so a future reader does not "fix"
   it into a double-count.
3. No signature change to `PrintDerivedAndCheck` (it already takes
   `const NucleationSpec& nuc` and `const GradualOverstressPerDOFParams&`).

### Edge Cases to Handle
- `nuc.enabled == false` → unchanged (every DOF "outside"; gate asks only the
  subcritical question). `direct_overstress` disabled behaves as gradual
  disabled.
- A `direct_overstress` patch whose `Δτ` is too small to reach yield → the
  inside-trigger fails → gate aborts (abort_on_failure) or warns
  (`SEAS_SKIP_EQUILIBRIUM_GATE=1`), exactly as gradual. This is the desired
  pre-flight protection.

### Acceptance Criteria
- [ ] For a `direct_overstress` config the gate selects the same DOFs
      inside/outside as an equivalent `gradual_overstress` config (same
      center/radii) and returns the same outside-max ratio.
- [ ] `test_spatial_print_derived` gradual cases pass unchanged; a new
      direct-kind case asserts inside-trigger / outside-subcritical PASS and a
      too-small-Δτ case asserts the trigger failure path.

### Dependencies
- Depends on: Phase 1 (accessors). Required by: Phase 4 (driver passes base
  pre-stress + nuc_params to the gate).

---

## Phase 4: Driver wiring in `spatial_dyn_driver.cpp`

### Goal
With `kind="direct_overstress"`, the driver folds `Δτ·F(r)` into `tau{1,2}_0`
at init, leaves `tau{1,2}_nuc ≡ 0` (no-op accumulator), runs the gate on the
base pre-stress, and publishes the patch footprint to PV — all selected purely
by config.

### Files to Modify
- `drivers/spatial_dyn_driver.cpp`.

### Detailed Requirements
1. **Resolve `nuc_params` by kind** (replace :1162-1167):
   ```cpp
   const spatial::GradualOverstressPerDOFParams nuc_params =
      (cfg.nucleation.kind == spatial::NucleationKind::DirectOverstress)
      ? spatial::ResolveDirectOverstress(cfg.nucleation.direct_overstress,
                                         cfg.nucleation.enabled,
                                         dof_coords_3d, dof_basis)
      : spatial::ResolveGradualOverstress(cfg.nucleation.gradual_overstress,
                                          cfg.nucleation.enabled,
                                          dof_coords_3d, dof_basis);
   ```
2. **Build the init pre-stress** between resolver (step 12) and
   `InitializeFaultDOFs_Spatial` (step 14, :1218-1225). The gate (step 15) keeps
   using the **base** `geom.GetTauPre()`:
   ```cpp
   Vector tau_pre_init(geom.GetTauPre());          // copy of base (2·N interleaved)
   if (cfg.nucleation.enabled &&
       cfg.nucleation.kind == spatial::NucleationKind::DirectOverstress &&
       num_fault_total > 0)
   {
      MFEM_VERIFY(nuc_params.amplitude_dip.Size() == num_fault_total &&
                  nuc_params.amplitude_strike.Size() == num_fault_total,
                  "direct_overstress: amplitude size mismatch");
      for (int i = 0; i < num_fault_total; ++i)
      {
         tau_pre_init(2*i)     += nuc_params.amplitude_dip(i);
         tau_pre_init(2*i + 1) += nuc_params.amplitude_strike(i);
      }
   }
   ```
   Pass `tau_pre_init` (not `geom.GetTauPre()`) to `InitializeFaultDOFs_Spatial`
   at :1222. Leave the `PrintDerivedAndCheck` call at :1302 reading
   `geom.GetTauPre()` (base) — Phase 3 invariant.
3. **No-op the accumulator for direct** (`nuc_cb`, :1749-1757): guard the
   `ApplyGradualOverstressIncrement` call with
   `cfg.nucleation.kind == NucleationKind::GradualOverstress`. For
   `direct_overstress` the lambda returns immediately (channels stay 0).
4. **Guard the shared-face ramp-as-step warning** (:1177-1190) with
   `kind == GradualOverstress` — the static bump is not ramped, so the
   1st-order-in-time degradation warning does not apply.
5. **Shared-fault consistency tripwire** (:1975-1978): the `t ≤ T_nuc_s`
   window is gradual-specific. Make it kind-aware:
   - gradual: keep `t <= cfg.nucleation.gradual_overstress.T_nuc_s` window;
   - direct: the static bump is fully present from step 0, so the `step == 0`
     branch already catches a cross-rank strike-sign mismatch in `tau{1,2}_0`.
     Do **not** reference `T_nuc_s` on the direct path.
   Suggested form: keep the `step == 0` unconditional branch; gate the periodic
   `step % 100` branch on `kind == GradualOverstress && t <= …T_nuc_s`.
6. **Startup print** (:740-742): replace the literal `"gradual_overstress
   (enabled)"` with a kind-dependent label
   (`"direct_overstress (enabled, static prestress)"` /
   `"gradual_overstress (enabled, smoothStep ramp)"` / `"disabled"`).
7. **PV publication** (:1507-1517) is unchanged — `nuc_params` is populated for
   both kinds, so `nuc_radial` / `nuc_amplitude` show the patch for direct too.

### Edge Cases to Handle
- `num_fault_total == 0` on a rank → `tau_pre_init` is the size-0 copy; the
  loop is skipped; `InitializeFaultDOFs_Spatial` is already gated on
  `num_fault_total > 0`.
- `SEAS_ZERO_DIP_PRESTRESS=1` + `delta_tau_dip_pa ≠ 0`: the static dip bump is
  re-added after the zero-dip slot — emit a one-time rank-0 WARNING when both
  hold, and document that a strike-only nucleation (`delta_tau_dip_pa = 0`)
  avoids it.
- `--dry-run --print-derived`: the gate must run on base pre-stress and still
  reflect the post-bump trigger (it does, via `nuc_params`). Verify a
  direct-kind dry run reports a triggering patch.

### Acceptance Criteria
- [ ] `seas_spatial_dyn_driver` builds; `kind="direct_overstress"` selects the
      static path, `gradual_overstress` is byte-identical to today (same log
      lines except the startup label).
- [ ] After init with a direct config, `dof_data[i].tau1_0 == base1_i +
      F_i·Δτ_dip` and `tau2_0 == base2_i + F_i·Δτ_strike`, while
      `tau1_nuc == tau2_nuc == 0` and stays 0 through `nsteps`.
- [ ] The `[derived]` gate prints the same inside/outside histogram for
      matched direct vs. gradual configs.
- [ ] A short MPI smoke run (`--tfinal` ~ a few `dt`) shows the patch DOFs
      already above yield at step 0 (rupture onset earlier than the equivalent
      gradual run).

### Dependencies
- Depends on: Phases 1-3. Required by: Phase 5 (docs/tests reference final
  driver behavior).

---

## Phase 5: Docs, example config, and test updates

### Goal
The schema doc and an example config show `direct_overstress`; unit tests cover
parse/resolve/gate/driver behavior; the existing single-kind assertions are
updated.

### Files to Create
- `safs/project_7.0_alternative/config/<example>_direct_overstress.toml` — copy
  of `spatial_friction_slip_weakening_safs_projected_stress.toml` with the
  `[nucleation]` block switched to `kind="direct_overstress"` +
  `[nucleation.direct_overstress]` (strike-only Δτ recommended).

### Files to Modify
- `safs/project_7.0_alternative/document/spatial_friction_config_schema.md` —
  document the new kind, the sub-block fields, the static-prestress semantics,
  the `T_nuc_s`-forbidden rule, and the `SEAS_ZERO_DIP_PRESTRESS` interaction.
- `tests/unit/test_spatial_friction_config.cpp` — keep T-21/22/23; add:
  T-24 `direct_overstress` parses (7 fields, kind enum); T-25 missing
  `[nucleation.direct_overstress]` aborts; T-26 stray `T_nuc_s` under direct
  aborts; T-27 non-positive radius aborts. (T-23's typo string stays a true
  typo, now also distinct from `direct_overstress`.)
- `tests/unit/test_spatial_nucleation.cpp` — add a `ResolveDirectOverstress`
  case (amplitude == F·Δτ) and assert the gradual cases are unchanged.
- `tests/unit/test_spatial_print_derived.cpp` — add a direct-kind PASS case and
  a too-small-Δτ FAIL case.
- The relevant `Makefile` test target list if these are new test fns in
  existing files (no new test binaries needed).

### Detailed Requirements
1. Schema doc: a `### [nucleation.direct_overstress]` section mirroring the
   gradual one, with a table of the 7 fields, the equation
   `tau{1,2}_0' = tau{1,2}_0 + F(r)·Δτ`, and a callout that it is the
   TPV205-style static patch (no time ramp; `tau*_nuc` stay 0).
2. Cross-link the driver runbook
   (`debug_document/spatial_workflow_safs_runbook_2026-05-18.md`) and update the
   CLAUDE.md "Spatial driver nucleation mechanism" section to say the driver
   supports **two** kinds.

### Acceptance Criteria
- [ ] `make test` passes (all spatial unit suites green).
- [ ] The example config passes `--dry-run --print-derived` (gate PASS) on a
      small mesh.
- [ ] Schema doc + CLAUDE.md describe both kinds accurately.

### Dependencies
- Depends on: Phases 1-4.

---

## Testing Strategy
- **Phase 1**: extend `test_spatial_friction_config.cpp` (pure parser; fast).
- **Phase 2**: extend `test_spatial_nucleation.cpp` — assert `ResolveDirect…`
  amplitudes against hand-computed `F·Δτ` at a probe DOF, and a regression guard
  that gradual outputs are unchanged after the shared-helper refactor.
- **Phase 3**: extend `test_spatial_print_derived.cpp` — equivalence of
  inside/outside selection for matched direct vs. gradual; trigger PASS and a
  deliberately-subcritical FAIL.
- **Phase 4**: a driver-level smoke (8 ranks, a few `dt`) asserting
  `tau{1,2}_0` carry the bump and `tau{1,2}_nuc` are zero throughout; compare
  rupture-onset step against the equivalent gradual run (direct should onset at
  step 0 in the patch).
- **Cross-check**: confirm `tpv205_*` / `tpv104_*` / `bp5_*` regression tests are
  untouched (no edits to those files; rerun the byte-exact suites if available).

## Risk Assessment
- **Double-counting at the gate** (highest risk): if someone later passes
  `tau_pre_init` (bumped) *and* `nuc_params` to `PrintDerivedAndCheck`, the
  inside trigger double-adds Δτ. Mitigation: the Phase-3 invariant comment +
  the Phase-4 rule "gate reads `geom.GetTauPre()`, init reads `tau_pre_init`".
  A unit assertion in Phase 3 (matched direct vs gradual gate output) detects a
  regression.
- **`SEAS_ZERO_DIP_PRESTRESS` × `delta_tau_dip_pa`**: silent re-introduction of
  dip stress under the zero-dip diagnostic. Mitigation: default dip Δτ to 0 in
  examples, emit the Phase-4 warning, document.
- **Shared-face strike-sign flip on the curvilinear SAFS fault**: the static
  bump is projected through the per-DOF strike vector, which can disagree across
  ranks on a shared face (the very issue the consistency tripwire guards). For
  direct, the mismatch is baked into `tau{1,2}_0` at init — the `step == 0`
  tripwire branch must remain active for direct (Phase 4 §5). If it fires,
  that is a real frame-consistency bug surfacing, not a new regression.
- **Refactor of `ResolveGradualOverstress`** (Phase 2): extracting the shared
  loop must keep gradual numerically identical. Mitigation: keep the helper a
  literal move of :104-127; rely on existing `test_spatial_nucleation` cases as
  the regression oracle (run before/after).
- **Tricky existing area**: the `tau_pre` interleaving `(2i=dip, 2i+1=strike)`
  vs. the `dof_basis` row layout `(3-5=dip, 6-8=strike)` — both verified above;
  re-check indices in code review (a swapped index silently rotates the bump
  into the wrong component and the gate would still "PASS" on magnitude).
```
