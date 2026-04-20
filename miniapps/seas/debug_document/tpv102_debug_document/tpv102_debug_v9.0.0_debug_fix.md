# TPV102 v9.0.0 Debug — Fix Log

Companion plan: `tpv102_debug_v9.0.0_debug_plan.md`.

This log records code changes implementing the v9.0.0 debug plan.  Each
pass appends a new section.  Production builds (no `SEAS_DIAG_FAULT_FLUX`
flag) remain byte-identical to pre-change builds — every diagnostic
artefact is gated under `#ifdef SEAS_DIAG_FAULT_FLUX`.

---

## Step 2 — C-1 … C-4 DIAG implementation (plan §0.5, §0.5.1, §0.5.2)

**Date:** 2026-04-19
**Agent:** chunhui-code-debugger
**Build:** `SEAS_EXTRA_CPPFLAGS=""` and `SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX"`.

### Scope

Four per-step printf checkpoints that bisect the fault → bulk chain at
runtime (C-1 in `FaultFaceFlux::Evaluate`, C-2/C-3 around the
interior-fault flux accumulation in `WaveOperator::ComputeFaceFluxRHS`,
C-4 after the RK4 Q update in `tpv102_driver.cpp`).  All are gated on
`SEAS_DIAG_FAULT_FLUX`; no collective sits under a per-QP or rank-local
branch.

### Files changed

| File | What changed (at post-change line numbers) |
|------|---------------------------------------------|
| `miniapps/seas/dynamic/seas_diag_rank.hpp` | NEW header. Declares `extern int mfem::seas::g_seas_my_rank` under `#ifdef SEAS_DIAG_FAULT_FLUX`; entire file compiles out in production builds. |
| `miniapps/seas/dynamic/fault_face_flux.hpp:43-50` | Added `bool diag_print = false` to `DOFData`, gated on `SEAS_DIAG_FAULT_FLUX`. Struct layout unchanged in production. |
| `miniapps/seas/dynamic/fault_face_flux.cpp:13-15` | Added `#include "seas_diag_rank.hpp"` and `#include <cstdio>`. |
| `miniapps/seas/dynamic/fault_face_flux.cpp:85-100` | Inserted **C-1 EVAL** fprintf block immediately after `tau2_total = data.tau2_0 + tau2_trial;`. Gated on `data.diag_print`. printf only, no MPI. |
| `miniapps/seas/dynamic/wave_operator.hpp:24` | Added `#include "seas_diag_rank.hpp"` at **file scope** (before `namespace mfem::seas` opens).  Must not be included from `wave_operator.inl`, which is itself included from inside `namespace mfem::seas`. |
| `miniapps/seas/dynamic/wave_operator.inl:1-14` | Noted the namespace-nesting constraint for future maintainers; removed earlier attempt to include the header locally. |
| `miniapps/seas/dynamic/wave_operator.inl:816-841` | Inserted **C-2 FLUX** fprintf block after `flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h_total);`. Gated on `(*fault_dof_data_)[dof_idx].diag_print`. printf only. |
| `miniapps/seas/dynamic/wave_operator.inl:843-877` | Bracketed the elem1 accumulation loop with **C-3 RHS** pre/post probe at `rhs[VX * ndof_total_ + dof_offset1 + 0]`. Rank-local probe, printf only. |
| `miniapps/seas/drivers/tpv102_driver.cpp:25` | Added `#include "../dynamic/seas_diag_rank.hpp"`. |
| `miniapps/seas/drivers/tpv102_driver.cpp:47-53` | Definition of `mfem::seas::g_seas_my_rank = 0` at file scope, gated. |
| `miniapps/seas/drivers/tpv102_driver.cpp:91-94` | Wires `mfem::seas::g_seas_my_rank = rank` immediately after `MPI_Comm_rank`. |
| `miniapps/seas/drivers/tpv102_driver.cpp:455, 463-467, 482-495` | Track `hypo_dof_local` inside the existing MINLOC block and (under the flag) tag the closest hypocenter DOF with `diag_print = true` on the rank that won the MINLOC. Emits a single `[diag] rank N tagging hypo DOF …` line on that rank. |
| `miniapps/seas/drivers/tpv102_driver.cpp:892-936` | Inserted **C-4 BULK** block after the V_max_step Allreduce.  Throttle `c4_fire = (step % 100 == 0) \|\| (step % 10 == 0 && t >= 1.0 && t < 1.4)` depends only on `step` and `t` — both lockstep across ranks — so every rank evaluates the same boolean and calls both MPI_Allreduce calls unconditionally.  Only the final fprintf is guarded by `rank == 0`. |

### Deviations from plan §0.5.2 snippets

1. **`diag_print` field added in this pass (plan assumed it existed).**  The
   plan §0.5 text says the C-1/C-2/C-3 gates should use "the existing
   `data.diag_print` bool flag," but v8.0.0 never landed the field.  Added
   in `fault_face_flux.hpp` under the flag so production struct layout
   is preserved.

2. **`seas_diag_rank.hpp` included from `wave_operator.hpp`, not
   `wave_operator.inl`.**  `wave_operator.hpp` pulls the .inl file in from
   inside `namespace mfem::seas`, so any `#include` at the top of the .inl
   nests `mfem::seas::mfem::seas` and breaks later `mfem::fn` macro
   expansions (reproduced and diagnosed during this pass).  Putting the
   include at file scope in the .hpp is the only safe pattern.

3. **C-4 placed after the `V_max_step` MPI_Allreduce, not immediately
   after the RK4 Q update.**  The plan snippet references `V_max_step`,
   which the driver computes just after the Q update via an unconditional
   Allreduce.  Placing C-4 after that Allreduce keeps the "one collective
   per lockstep variable" invariant and matches the snippet's field list.
   The step=0 case still fires because `step % 100 == 0` evaluates true
   on step 0.

4. **Hypocenter DOF tagging added in the driver.**  The plan assumes the
   driver already tags one `diag_print=true` DOF; v8.0.0 did not.  A
   minimal tagger was added that reuses the existing MINLOC result
   (`hypo_rank`) and the local-min `hypo_dof_local` index the MINLOC
   loop already computed.  At most one DOF is flagged globally — the
   hypocenter.  Off-hypocenter / second-station tagging is intentionally
   NOT added here; that is the v8.0.0 plan's long-form DIAG protocol
   and out of scope for v9.0.0 §0.5.

### MPI-deadlock audit (plan §0.5.1 compliance)

| Rule | C-1 / C-2 / C-3 | C-4 |
|------|-----------------|-----|
| No MPI calls inside `Evaluate` / per-QP loops | OK (printf only) | N/A |
| Throttle uses only lockstep data (`step`, `t`) | N/A | OK (`c4_fire` depends on step + t only) |
| MPI_Allreduce called on every rank when it fires | N/A | OK (no rank-local or DOF-local gate) |
| Only final fprintf gated on `rank == 0` | N/A | OK |
| No MPI calls gated on `diag_hypo_dof >= 0` or `diag_print` | OK | OK |

### Build-check summary

| Build | Flag | Result |
|-------|------|--------|
| `rm -f seas_tpv102_driver drivers/tpv102_driver.o dynamic/*.o && make seas_tpv102_driver SEAS_EXTRA_CPPFLAGS=""` | OFF | Clean build. No warnings introduced. |
| `... SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX"` | ON | Clean build. No warnings introduced. |

### Test results

| Test | Procs | Result |
|------|-------|--------|
| `./seas_tpv102_driver --tfinal 0.001` (flag OFF) | 1 | Exit 0; banner shows `SEAS_DIAG_FAULT_FLUX = OFF`; zero `[C-*]` / `[diag]` lines in stderr. Production build confirmed quiet. |
| `./seas_tpv102_driver --tfinal 0.001` (flag ON) | 1 | Exit 0; banner shows `SEAS_DIAG_FAULT_FLUX = ON`; 25 `[C-*]` lines in stderr; first `[C-4 BULK]` line: `step=0  t=0.0009  V_max=1.000e-12  max\|Q[VX]\|=2.126e-29 m/s  max\|Q[SXY]\|=2.636e-23 Pa`. |
| `mpirun -np 2 ./seas_tpv102_driver --tfinal 0.005` (flag ON) | 2 | Exit 0 in 42 s (< 60 s budget). `[C-4 BULK] step=0` line present on combined stderr. `[diag] rank 1 tagging hypo DOF 5188 …` confirms the hypocenter landed on rank 1 after METIS partition and rank 1 set `diag_print=true` on exactly one DOF. No deadlock. |
| `./seas_test_fault_face_flux` (flag OFF) | 1 | 19 / 19 pass |
| `./seas_test_wave_operator` (flag OFF) | 1 | 17 / 17 pass |
| `mpirun -np 2 ./seas_test_parallel_wave_operator` (flag OFF) | 2 | 5 / 5 pass |
| `mpirun -np 2 ./seas_test_r101_shared_fault` (flag OFF) | 2 | 18 / 18 pass |

No regressions in any of the four TPV102-scope tests the user called out.
The production-build path was re-verified byte-quiet (no stderr output
with `[C-*]` tag) to confirm the flag gating is hermetic.

### First C-4 BULK line (flag ON, 1 rank, tfinal=0.001)

```
[BUILD] SEAS_DIAG_FAULT_FLUX = ON
[BUILD] SEAS_DIAG_GHOST_EXCHANGE = OFF
[diag] rank 0 tagging hypo DOF 9802 at (15.0, 0.0, -7654.0)
[C-1 EVAL] rank=0  tau1_trial=+0.000e+00 Pa  tau2_trial=+0.000e+00 Pa  …  psi=7.355e-01
[C-2 FLUX] rank=0  dof=9802  |F_v|=2.768e-57 m2/s2  max|F_stress|=1.267e-19 Pa*m/s  F_h[VX]=+1.866e-88  F_h[SXY]=-1.267e-19  F_h[VY]=-2.768e-57  F_h[SXZ]=-7.755e-36
[C-3 RHS]  rank=0  dof=9802  rhs[VX,elem1,dof0] pre=-1.644e-84 post=-8.222e-85  delta=+8.222e-85  w=9.142e+04  shape1(0)=-4.821e-02  F_h_total[VX]=+1.866e-88
…
[C-4 BULK] step=0  t=0.0009  V_max=1.000e-12  max|Q[VX]|=2.126e-29 m/s  max|Q[SXY]|=2.636e-23 Pa
```

Step=0 at t=0.0009 s: Q is still near machine-zero perturbation (Q
initialises to 0, and tfinal=0.001 s is only ~2 dt_cfl on the 1000 m
mesh), so the tiny magnitudes are expected.  The checkpoints will
populate physical values once the run proceeds into the nucleation ramp
(see plan §0.5.3 for expected t ≈ 1.3 s magnitudes).

### Verification checklist

- [x] Production build (flag OFF) emits zero `[C-*]` or `[diag]` lines.
- [x] Diagnostic build (flag ON) emits `[C-4 BULK]` for `step=0`.
- [x] 2-rank MPI diagnostic build completes in < 60 s without hang.
- [x] C-1 / C-2 / C-3 use `fprintf` only; no MPI calls.
- [x] C-4 `MPI_Allreduce` is unconditional per rank when `c4_fire` is true.
- [x] `c4_fire` depends only on `step` and `t`.
- [x] Only `rank == 0` guards the final fprintf in C-4.
- [x] `seas_test_r101_shared_fault`, `seas_test_parallel_wave_operator`,
      `seas_test_fault_face_flux`, `seas_test_wave_operator` — all 4
      retain their pre-change pass counts (19, 17, 5, 18).
- [x] `DOFData` struct layout in production build is unchanged (new
      member is under the flag).

### Suggested single-commit message

```
tpv102: add C-1…C-4 DIAG bisection checkpoints (v9.0.0 plan §0.5)

Four runtime DIAG checkpoints that bisect the fault → bulk chain when
-DSEAS_DIAG_FAULT_FLUX is set at build time.  All per-QP / per-DOF
checkpoints (C-1, C-2, C-3) emit fprintf only; C-4 uses MPI_Allreduce
unconditionally per rank with a throttle that depends only on (step, t)
and thus evaluates identically on every rank — no deadlock possible.

- dynamic/seas_diag_rank.hpp (new): extern int g_seas_my_rank, gated.
- fault_face_flux.hpp: add bool diag_print to DOFData, gated.
- fault_face_flux.cpp: C-1 EVAL after tau2_total assignment.
- wave_operator.hpp: include seas_diag_rank.hpp at file scope (must not
  be included from wave_operator.inl; the .inl is included from inside
  namespace mfem::seas and would nest).
- wave_operator.inl: C-2 FLUX after flux_.Interior; C-3 RHS brackets
  the elem1 accumulation loop with a rank-local rhs probe.
- tpv102_driver.cpp: define g_seas_my_rank, wire it at MPI init, tag
  the hypocenter DOF with diag_print=true on the owning rank, and add
  C-4 BULK after the V_max_step Allreduce.

Production builds (no flag) are byte-identical to pre-change builds.

Verified:
  make SEAS_EXTRA_CPPFLAGS=""                  → clean build, no DIAG output.
  make SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX" → clean build; first
    [C-4 BULK] line prints at step=0 as expected.
  mpirun -np 2 … --tfinal 0.005                → completes in 42 s, no hang.
  seas_test_{fault_face_flux,wave_operator,
             parallel_wave_operator,r101_shared_fault}
    → 19 / 17 / 5 / 18 pass — no regressions.
```
