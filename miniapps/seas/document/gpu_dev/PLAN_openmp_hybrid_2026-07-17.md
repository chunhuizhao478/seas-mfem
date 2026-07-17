# PLAN — OpenMP hybrid (threaded) parallelism for the compute hot path

**Date:** 2026-07-17
**Status:** proposed (not implemented) — larger, higher-risk than the sidecar
shared-memory plan; read that first if the goal is only memory.
**Owner folder:** `document/gpu_dev/` (same `mfem::forall`/`Device`-backend
mechanism as the GPU refactor).
**Companion:** `document/code_optimization_dev/PLAN_sidecar_mpi_shared_memory_2026-07-17.md`
(the cheap memory fix), `document/code_optimization_dev/profiling_analysis_2026-06-24.md`
(the hot-path profile this plan is built on).

---

## 1. Problem / motivation

`seas_spatial_dyn_driver` is **pure-MPI** (no OpenMP anywhere; `config/config.mk:30`
`MFEM_USE_OPENMP = NO`; plain `MPI_Init` at `drivers/spatial_dyn_driver.cpp:750`). To
use all 128 cores of an Expanse node it must run 128 MPI ranks/node, versus SeisSol's
hybrid 8 ranks × 16 threads. Two consequences motivate a hybrid mode:

1. **Memory** — 128 ranks replicate the read-only sidecars 128× (the `--mem` bump in
   the v4_0_0 benchmark). **But this is better solved by the shared-memory plan** — see §3.
2. **Compute/comm scaling** — fewer, fatter ranks mean less MPI surface: fewer halo
   exchanges and fewer ranks in the guard `Allreduce`s that the profile shows dominate
   the post-lever residual (`profiling_analysis_2026-06-24.md`, §5). This is the real,
   memory-independent reason to want a hybrid.

## 2. Honest framing — what OpenMP does and does NOT buy

From `document/code_optimization_dev/profiling_analysis_2026-06-24.md` (443.9 s total):

| Hot region | Share | Threadable? |
|---|---|---|
| `ApplySpatialDerivative` (runs **2×/step**) | **53.4 %** | Yes — element-local, no atomics |
| `ComputeVolumeRHS` | 10.1 % | Yes — element-local |
| Face-flux correctors | ~22.6 % (**mostly MPI wait**) | Partly — see §6 |
| `friction_substep` (Brent) | **0.89 %** | Yes but low value |

- The `--shared-ck-recursion` lever already deletes the redundant second
  `ApplySpatialDerivative` subtree (~221 s) **bit-exact**; the v4_0_0 sbatch enables it.
- **After the levers, the residual is MPI-imbalance-bound** (guard `Allreduce`s + halo
  `Waitall` concentrated on fault-bearing ranks — `profiling_analysis_2026-06-24.md`
  §5). **OpenMP threading does NOT reduce that imbalance.** It attacks the
  compute-bound bulk (`ApplySpatialDerivative` + `ComputeVolumeRHS` ≈ 63 %). So the
  realistic upside is speeding the ~63 % compute-bound fraction and shrinking the MPI
  surface — **not** curing the load imbalance (that is the fault-weighted-partition
  work, `document/code_optimization_dev/action_plan_fault_weighted_partition_*`).

**Therefore:** if the objective is *memory*, do the shared-memory plan (1–2 days, low
risk) — not this. Pursue OpenMP only if the objective is genuine intra-node compute/comm
scaling or execution-model parity with SeisSol.

## 3. Goal / non-goals

**Goal.** Let the driver run `M` MPI ranks × `T` OpenMP threads per node (e.g. 8×16),
threading the element-local compute kernels via MFEM's `omp` `Device` backend, so the
compute-bound bulk scales on-node and the MPI rank count (hence halo/`Allreduce`
surface and sidecar replication) drops.

**Non-goals.** Not a rewrite of the numerics. Not the MPI load-imbalance cure. Not a
GPU port (that is the separate `gpu_full_port_plan`). Not byte-exact with the pure-MPI
OnTheFly path (see §7 — this is the load-bearing caveat).

## 4. Mechanism — reuse the GPU-refactor `forall`, switch the backend to `omp`

MFEM's `mfem::forall` dispatches through the `Device` singleton's active backend; under
the **`omp` backend it becomes an OpenMP-parallel loop automatically**. The
GPU-refactor worktree already wrote the element-local kernels as `mfem::forall`, so the
*same* kernels thread on CPU with `--device omp` once MFEM is built with
`MFEM_USE_OPENMP=YES`. `dynamic/gpu/seas_device.hpp:48,97` already passes any
non-shorthand device string straight to `mfem::Device::Configure`, so `--device omp`
reaches MFEM with no new plumbing.

## 5. Current state (assets that already exist)

**GPU worktree** `.claude/worktrees/spatial-dyn-driver-gpu/miniapps/seas/dynamic/gpu/`
(branch `worktree-spatial-dyn-driver-gpu`, HEAD `4218cfa`) — element-local `mfem::forall`
kernels, validated (per memory: 45/45 cpu+debug, regression 54/54, byte-exact CPU):

- `wave_device_kernels.hpp` — `ApplySpatialDerivativeDevice` (`:142`),
  `ComputeVolumeRHSCachedDevice` (`:194`), `ApplyMassInverseDevice` (`:248`, uses
  **out-of-place** device scratch to avoid an in-place race, `:272,:286`). One
  thread per `(e,c,i)`; DG DOFs are element-unique → **no atomics**.
- `wave_face_kernels.hpp` — `ComputeInteriorFaceFluxCachedDevice` (`:100`), interior
  **non-fault** faces only, probe-extracted 9×9 flux, **`mfem::AtomicAdd`** for the
  two-element scatter (`:172-173`); **ADER-only**, needs `--face-cache`.
- `seas_device.hpp` — `ConfigureSeasDevice` → `mfem::Device::Configure`; driver wiring
  at `drivers/spatial_dyn_driver.cpp:549-551` (`--device`), `:1348-1364`
  (`wave.EnableDeviceKernels(true)`, gated on `--deriv-cache`).

**What stays on the host even in the GPU worktree:** the **fault** face flux
(`dynamic/fault_face_flux.cpp`), the friction Brent solve
(`friction/dieterich_ruina.hpp` via `dynamic/friction_substep_iterator.*`), and
**shared/cross-rank** faces. A hybrid run keeps these serial-per-rank; that is fine at
8 ranks/node (each rank still owns a fault chunk) but means fault-heavy ranks do not
thread their fault work.

**Main checkout** `dynamic/` has **no** `forall`/`Device`/`AtomicAdd` — only
`DerivMode::Cached` hand-loops. So step 1 of any hybrid is to bring the worktree's
`gpu/` kernels (or an OMP-targeted subset) into a buildable branch.

## 6. Phased implementation

**Phase 0 — build + backend prerequisites.**
- Rebuild MFEM and the miniapp with `MFEM_USE_OPENMP = YES` (`config/config.mk:30`).
- Upgrade `MPI_Init` → `MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided)`
  (`drivers/spatial_dyn_driver.cpp:750`) and assert `provided >= FUNNELED`. FUNNELED
  suffices because every MPI call (halo exchange, `Allreduce`, `Allgatherv`) currently
  runs **outside** any parallel region; do **not** thread the comm itself in this plan.
- Land the GPU worktree's `dynamic/gpu/` element-local kernels + `--device`/
  `EnableDeviceKernels` wiring onto the working branch (they are already CPU-validated).

**Phase 1 — element-local, deterministic kernels first (the safe 63 %).**
Enable `--device omp` for the no-atomic kernels only: `ApplySpatialDerivativeDevice`,
`ComputeVolumeRHSCachedDevice`, `ApplyMassInverseDevice`. These have unique `(e,c,i)`
writes → **thread-order-independent → still bit-exact vs the cached serial path**
(same reductions, just parallelized). Gate on the parity tests (§7). This alone threads
the dominant `ApplySpatialDerivative` (53 %) + `ComputeVolumeRHS` (10 %).

**Phase 2 — interior face flux (atomic, tolerance-gated).**
Enable `ComputeInteriorFaceFluxCachedDevice` (needs `--face-cache`). `mfem::AtomicAdd`
threads correctly but its **summation order is non-deterministic** → run-to-run
bitwise variation. This path is **explicitly not byte-exact** (≤1e-11, per
`wave_face_kernels.hpp:22`). Gate on the ≤1e-12/≤1e-11 parity tests, **not** the TOL=0
gates (§7). Optional: a deterministic per-element gather (colour faces by owning
element, sum in fixed order) if reproducibility is required — more work, removes the
non-determinism.

**Phase 3 — BimaterialWaveOperator (the SAFS path) confirmation.**
SAFS uses `interior_flux="matrix"` → **`BimaterialWaveOperator`**
(`dynamic/bimaterial_wave_operator.inl`), a different code path than the homogeneous
operator the GPU kernels were validated against. Confirm the cached/forall kernels are
wired and parity-tested for the bimaterial path
(`test-bimaterial-deriv-cache-parity`, `test-bimaterial-wave-operator-parity`) before
using OMP for a SAFS production run. Per `profiling_analysis_2026-06-24.md:155` the
matrix path needs its own cached-vs-OnTheFly + fingerprint check.

**Phase 4 — (optional, low value) friction Brent.**
`RunSubSteps_` inner loop over fault DOFs (`friction_substep_iterator.hpp:227`) is
per-DOF independent (each touches only `dof_data[i]`, no shared state, no MPI
collective) → embarrassingly parallel. But it is **0.89 %** of runtime; thread it only
if profiling a specific fault-heavy rank shows it matters. Keep Brent (never Newton,
per `CLAUDE.md`).

**Phase 5 — rank/thread sweep + benchmark.**
Sweep `M×T` on 2 nodes (e.g. 8×16, 16×8, 32×4) vs the 128×1 pure-MPI baseline on the
v4_0_0 p3 case; measure steps/s and per-node RSS. Expect the hybrid to (a) cut sidecar
memory as a side effect and (b) win or tie on compute-bound cases, with the crossover
governed by the MPI-imbalance floor (which neither pure-MPI nor hybrid removes).

## 7. Regression contract — the load-bearing caveat

The byte-exact gates **will fail** on any threaded/cached path:

- `make quick-check` (`Makefile:5085`) — 8-rank/20-step BP5 vs golden, **`QUICK_CHECK_TOL=0`**.
- `make regression-check` (`Makefile:5072`) — golden compare, default TOL 0.

These pass **only** on the OnTheFly host path. The threaded kernels are float-reassociated
(≤1e-12 deriv / ≤1e-11 face, per `wave_device_kernels.hpp:12`,
`wave_face_kernels.hpp:22`) and the face AtomicAdd is non-deterministic. So the gate for
this plan is the **tolerance-based parity suite**, not the TOL=0 suite:

- `test-bimaterial-deriv-cache-parity` (`Makefile:5038`, np2, ≤1e-12) — the closest
  existing analog; the primary gate.
- `test-wave-operator-cached-parallel` (`:5024`), `test-wave-operator-spatial-derivative`
  (`:4740`), `test-ader` (`:4776`), `test-ader-tpv102-smoke` (`:4764`),
  `test-bimaterial-wave-operator-parity` (`:4369`), `test-v92-regression-gates`
  (`:4809`), `test-bp5-smoke` (`:5062`).

**The plan must state up front** that hybrid/threaded output is *tolerance-equal, not
byte-equal*, to the pure-MPI baseline, and that Phase 2's atomic path is *not*
run-to-run reproducible unless the deterministic-gather option is taken. Any SAFS
science conclusion drawn from a hybrid run must note which path produced it.

`CLAUDE.md` "Files Requiring Extreme Care" (`friction/dieterich_ruina.hpp`,
`fault/fault_basis.hpp`, `fault/rate_state_fault.hpp`, `solver/*`, …) and its
"do NOT mix structural refactoring with numerical changes in one commit" both apply:
land the `forall` kernels (structural) and the OMP enablement (behavioral) as separate,
individually-gated commits.

## 8. Risks

| Risk | Mitigation |
|---|---|
| Threaded/atomic path not byte-exact → fails TOL=0 gates | Gate on the ≤1e-12 parity suite; document tolerance-equality explicitly (§7). |
| Face AtomicAdd non-deterministic → non-reproducible runs | Phase 2 optional deterministic per-element gather; or restrict OMP to Phase-1 element-local kernels for reproducible runs. |
| SAFS matrix/bimaterial path unvalidated for forall | Phase 3 parity confirmation before any SAFS production use. |
| MPI+threads correctness | `MPI_THREAD_FUNNELED` + keep all MPI outside parallel regions (they already are). |
| Effort creep into the numerically-sensitive fault/friction kernels | Scope OMP to element-local bulk first; fault flux + friction are last and optional. |
| Does not fix the real bottleneck (MPI imbalance) | Set expectations (§2); pair with fault-weighted partition if imbalance dominates. |

## 9. Effort & recommendation

**Scope:** land + wire the GPU-worktree `forall` kernels, MFEM OMP rebuild,
`MPI_Init_thread`, per-phase parity gating, bimaterial confirmation, rank/thread sweep.
**Estimate:** ~1–3 weeks, touching the numerically-sensitive compute path under the
tolerance-parity contract (vs ~1–2 days for the shared-memory plan).

**Recommendation:**
1. **If the goal is memory** → do the **sidecar shared-memory plan**, not this.
2. **If the goal is intra-node compute/comm scaling or SeisSol-hybrid parity** → do this,
   but **only Phases 0–1 initially** (the deterministic, bit-exact-vs-cached element-local
   kernels covering ~63 % of runtime), on the homogeneous path, gated on
   `test-bimaterial-deriv-cache-parity`. Treat Phase 2 (atomic face flux),
   Phase 3 (SAFS bimaterial), and Phase 4 (friction) as separate, individually-justified
   follow-ups. Expect the payoff to be capped by the MPI load-imbalance floor, which is a
   different project (fault-weighted partition).
