# Implementation Plan: `[FRAME]` canonical-frame diagnostic at the seed shared QPs — 2026-05-23

## Overview

Add a single runtime-env-gated `[FRAME]` stderr trace inside the shared-fault predictor
rotation (`WaveOperator::EvaluateBulkAtFaultQPsCanonical`, shared branch,
`wave_operator.inl:2143–2238`) that, for the seed QPs, reports the canonical fault frame's
**orthonormality**, its **`sign_flipped` / degenerate-band conditioning**, and the
**decomposition of the global bulk velocity jump** into normal vs tangential components.
This is "test #1 / #6" from `DEBUG_speckle_normal_velocity_jump_2026-05-23.md`: it decides
whether the runaway-driving `[[v_n]]` is a frame defect, a `sign_flipped` instability, or a
genuine bulk opening — the three hypotheses that survive the cap-run refutation (job 7747835).

This is a **diagnostic only**. It must be byte-exact when its env switch is unset.

## Background (why this site, why these fields)

- The iterator's `[SLIP]` trace (`tpv205_substep_iterator.cpp:415–448`) reports
  `sn_vjump = eta_p·(Q̃⁻[VX] − Q̃⁺[VX])` — the normal-velocity-jump term, in the **canonical
  frame** (VX = normal). It does **not** expose the frame vectors or the tangential parts of
  the jump.
- Those `Q̃±` are produced by `EvaluateBulkAtFaultQPsCanonical` (shared branch), which:
  - reconstructs `can_n/can_t1/can_t2` from `FaultBasisQPData` (`wave_operator.inl:2196–2205`),
  - rotates the **global** `Q_self`/`Q_nbr` (filled at `:2158–2173`) into the canonical frame
    via `Tinv_can` (`:2207–2222`),
  - routes ± by `elem1_on_plus` (`:2132`, `:2226–2229`) and packs into `Q_plus_flat`/
    `Q_minus_flat` (= the iterator's `Q̃±`).
- `fault_basis.hpp` builds `{normal,tangent1,tangent2}` **per-QP** (`ComputeQPBasisShared`:
  `CalcOrtho` at each `ir.IntPoint(q)`) and **orthonormal by construction** (R-003 explicit
  normalize of both strike and dip). So:
  - orthonormality residuals are expected `~1e-15`; a large residual would be a real bug.
  - the SAFS fault normal is ~⊥ the hardcoded `ref_normal=(0,−1,0)` (memory
    `project_safs_orientation_blowup`), i.e. in the **degenerate band** of
    `NormalNeedsFlipToCanonical` (`fault_basis.hpp:362–379`) where `sign_flipped` is set by a
    largest-component fallback — a plausible instability across adjacent QPs/ranks. Hence we
    print `sign_flipped`, `nl`, and `can_n` so the band can be seen (`n·ref = −can_n[1]`).

### Equations measured
With the same ± routing the iterator uses (so `dv_n` cross-checks `[SLIP]`):
```
vel_plus_g  = elem1_on_plus ? vel(Q_self) : vel(Q_nbr)      // global x/y/z velocity
vel_minus_g = elem1_on_plus ? vel(Q_nbr)  : vel(Q_self)
dvg         = vel_minus_g − vel_plus_g                       // global velocity jump (minus−plus)
dv_n        = dvg · can_n      (≡ Q̃⁻[VX] − Q̃⁺[VX] = sn_vjump/eta_p   — built-in cross-check)
dv_t1       = dvg · can_t1 ,  dv_t2 = dvg · can_t2            // tangential parts of the SAME jump
|dvg|       = ‖dvg‖
orthonormality:  n·t1, n·t2, t1·t2  (→0) ;  |n|−1, |t1|−1, |t2|−1  (→0)
```
Reading: orthonormal frame + `dv_n ≈ |dvg|` ⇒ genuine bulk opening (imposed-state/no-opening
locus). Orthonormal + `dv_n ≪ |dvg|` ⇒ jump is mostly tangential (slip) and the normal leak
is small at *this* sub-step (feedback inflates it later). Large orthonormality residual or
`sign_flipped` flapping across the seed QPs ⇒ frame/`sign_flipped` is the defect.

## Constraints

- **Byte-exact when `SEAS_DIAG_FRAME` is unset.** No change to any computed state
  (`Q_plus_flat`/`Q_minus_flat`, rotations, routing). Only a gated `fprintf` is added. This
  preserves the TPV205/BP5 regression contract (the routine is shared with TPV* via the same
  template).
- **Runtime env-gated** (matches the existing `[SLIP]` / `SEAS_NOOPENING` style —
  `tpv205_substep_iterator.cpp:359–373`), **not** the compile-time `SEAS_DIAG_FAULT_FLUX`
  path. The user toggles via sbatch env; no recompile, no `diag_print` DOF tagging.
- **No struct changes** (no `DOFData`/`FaultBasisQPData` field — R-004 layout invariant).
- **Single insertion point**, shared branch only (the seed QPs are `is_shared=1`).
- **No hardcoded mesh-specific constants in logic**; the target coordinate/radius come from
  env with a documented default equal to the seed (consistent with `SEAS_DIAG_XRANK_QP`).

## Phase 1: the `[FRAME]` trace

### Goal
After this phase, setting `SEAS_DIAG_FRAME=1` prints one `[FRAME]` line per shared-fault QP
within `SEAS_DIAG_FRAME_R` of `SEAS_DIAG_FRAME_XYZ`, every time the shared predictor is
evaluated; the run is byte-identical to baseline when the switch is unset.

### Files to Modify
- `miniapps/seas/dynamic/wave_operator.inl` — add the gated trace in the shared branch of
  `EvaluateBulkAtFaultQPsCanonical`, immediately after `can_n/can_t1/can_t2` are built
  (after `:2205`). Add static env-parse lambdas at the top of the shared per-QP region (or
  just before the `for (int q …)` loop at `:2143`).
- `miniapps/seas/dynamic/wave_operator.cpp` (or `.hpp`) — ensure `<cstdlib>` (getenv/atof),
  `<cstdio>` (fprintf/sscanf), `<cmath>` (sqrt) are included in the TU. `<cstdio>` is already
  used (the `[C-1s INT BASIS]` dump). Add the others if absent.

### Detailed Requirements

1. **Env gates** (parsed once via function-local `static const` lambdas, file-scope-safe):
   - `frame_diag` (bool): `getenv("SEAS_DIAG_FRAME")` truthy (non-empty, not `"0"`).
   - target `(x,y,z,r²)`: parse `getenv("SEAS_DIAG_FRAME_XYZ")` as `"%lf,%lf,%lf"`
     (default `607518, 3706359, −4543`); `getenv("SEAS_DIAG_FRAME_R")` → radius (default
     `300.0` m), stored squared. If `SEAS_DIAG_FRAME_XYZ` fails to parse 3 values, keep the
     default and proceed (no abort).
   - `s_frame_rank` (int): `MPI_Comm_rank(MPI_COMM_WORLD)` under `#ifdef MFEM_USE_MPI` and
     `MPI_Initialized`; else `0`. Mirror `tpv205_substep_iterator.cpp:367–373`.

2. **Per-QP gate**: inside the `q` loop, only when `frame_diag`:
   - compute physical coords `Vector xq(3); ftr->Transform(ip, xq);`
   - if `‖xq − target‖² ≤ r²`, emit the trace; else skip.
   (Guarding the `Transform` call behind `frame_diag` keeps the unset path free of even that
   call.)

3. **Quantities** (all from in-scope locals — `can_n/can_t1/can_t2`, `Q_self`, `Q_nbr`,
   `elem1_on_plus`, `qpd`, `base_dof_idx`, `q`):
   - orthonormality: `n·t1, n·t2, t1·t2`, and `‖n‖−1, ‖t1‖−1, ‖t2‖−1`.
   - `vp = elem1_on_plus ? Q_self : Q_nbr`, `vm = elem1_on_plus ? Q_nbr : Q_self`.
   - `dvg[d] = vm[VX+d] − vp[VX+d]` for `d=0,1,2` (uses `QIndex` VX=6,VY=7,VZ=8 from
     `wave_state.hpp`).
   - `dv_n = Σ can_n[d]·dvg[d]`, `dv_t1 = Σ can_t1[d]·dvg[d]`, `dv_t2 = Σ can_t2[d]·dvg[d]`,
     `dvg_mag = sqrt(Σ dvg[d]²)`.
   - `dof_idx = base_dof_idx + q` (for cross-ref with `[SLIP]` qp index).

4. **Output** — one `fprintf(stderr, …)` + `fflush`, fields in this order:
   ```
   [FRAME] qp=%d c=(%.1f,%.1f,%.1f) rank=%d sign_flipped=%d nl=%.6e
           n.t1=%+.3e n.t2=%+.3e t1.t2=%+.3e d|n|=%+.3e d|t1|=%+.3e d|t2|=%+.3e
           can_n=(%+.6e,%+.6e,%+.6e) |dvg|=%.6e dv_n=%+.6e dv_t1=%+.6e dv_t2=%+.6e
   ```
   (single call, `\n`-terminated; multi-line literal acceptable). `dv_n` must equal the
   matching `[SLIP]` line's `sn_vjump/eta_p` to ~1 ULP (cross-check, see AC).

### Edge Cases to Handle
- **Switch unset** → `frame_diag=false` → no coords computed, no print, byte-exact.
- **Serial build** (`MFEM_USE_MPI` undefined) → `s_frame_rank=0`; compiles and runs.
- **Malformed `SEAS_DIAG_FRAME_XYZ`** → keep default target, do not abort.
- **`SEAS_DIAG_FRAME_R` ≤ 0 or unparsable** → keep default radius.
- **No QP within radius** (target off this rank's partition) → simply no output on that rank;
  the owning rank prints. (Coords are global, so the correct rank self-selects.)
- **`nl` near zero** (degenerate normal) → printed verbatim (it is the diagnostic signal); no
  division by `nl` in the trace.

### Acceptance Criteria
- [ ] With `SEAS_DIAG_FRAME` **unset**, TPV205 regression is byte-exact (parity 14/14) and
      the friction suite (75/75) unchanged — i.e. the build with the patch produces identical
      output to without it when the env is unset.
- [ ] Code compiles under `conda activate mfem-dev` (`make seas_spatial_dyn_driver`).
- [ ] A unit test asserts the trace's arithmetic identity: for a constructed
      `{can_n,can_t1,can_t2}` (orthonormal) and global `Q_self/Q_nbr`,
      `dv_n == can_n·(vel_minus−vel_plus)` and equals `(Q̃⁻−Q̃⁺)[VX]` after `Tinv_can`
      rotation, to `< 1e-12` relative. (Validates the ± routing + projection match the
      iterator's `sn_vjump`.)
- [ ] On Frontera with `SEAS_DIAG_FRAME=1`, `[FRAME]` lines appear at the seed coords on
      rank 104 and `dv_n` matches the `[SLIP]` `sn_vjump/eta_p` at the same `t`/`qp`.

### Dependencies
- Depends on: nothing (self-contained diagnostic).
- Required by: the frame-vs-bulk-opening determination in the debug doc (test #1/#6).

## Testing Strategy
1. **Local compile** (`make seas_spatial_dyn_driver -j`) under `mfem-dev`.
2. **Byte-exact check**: run an existing TPV205 unit/regression target with the patched
   binary, `SEAS_DIAG_FRAME` unset, and confirm parity unchanged (the change is inert).
3. **New unit test** (`miniapps/seas/dynamic/` test dir or the existing frame/rotation test):
   `test_frame_diag_projection_identity` — build an orthonormal frame + two global Q vectors,
   verify `dv_n` (the trace's formula) equals `Tinv_can`-rotated `(Q⁻−Q⁺)[VX]` to 1e-12. This
   guards the ± routing and the projection against the iterator's `sn_vjump`.
4. **Frontera**: add `SEAS_DIAG_FRAME=1` (default target = seed) to the Dc2 sbatch; one short
   run past onset; grep `[FRAME]`, cross-check `dv_n` vs `[SLIP] sn_vjump`.

## Risk Assessment
- **`ftr->Transform` cost** — gated behind `frame_diag`, so zero in production. Per-QP only
  for the ≤~5 seed QPs when enabled.
- **Index/coords mismatch with `[SLIP]`** — mitigated by printing **both** `dof_idx` and
  physical coords; the coords are the source of truth for cross-referencing.
- **`ref_normal` not in scope** for an explicit degenerate-band flag — mitigated by printing
  `can_n` and `nl`; `n·ref = −can_n[1]` for SAFS, readable offline. (Do not hardcode
  `ref_normal` in the trace.)
- **Shared vs interior** — the trace is shared-branch only; if a future seed QP is interior
  (`is_shared=0`) the same block would be added at `:1922`. Out of scope now.
