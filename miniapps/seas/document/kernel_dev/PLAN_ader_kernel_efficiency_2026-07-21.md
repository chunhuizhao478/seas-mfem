# Implementation Plan: MFEM-native ADER element-kernel efficiency program

**Date:** 2026-07-22 (**rev 3**, post-evidence-round) · **Branch:** `safs-v4_0_0-alt-case1-mfem-speed`
**Scope:** `seas_spatial_dyn_driver` only
**Motivating measurement:** `document/lts_dev/RESULTS_p5_tpv104_200m_2026-07-21.md` (Phase-5 three-way)
**Prior-art verdicts referenced:** `document/code_optimization_dev/action_plan_2026-06-24.md` (the OPT-* IDs)
**Review record:** `document/kernel_dev/REVIEW_plan_ader_kernel_efficiency_2026-07-21.md` (round 1 + round 2)
**Evidence record (rev-3 basis):** `document/kernel_dev/EVIDENCE_microbench_roofline_2026-07-22.md`
— a working microbenchmark (`tests/bench/bench_ader_kernel_variants.cpp`), the Rome roofline
arithmetic, and the end-to-end composition table. **The Phase-0.5 prototype spike demanded by
review R-402 has ALREADY RUN locally**; its results are baked into this rev.
**Companion program:** `document/comm_dev/PLAN_lts_comm_reduction_2026-07-22.md` — the other half
of the runtime claim; see "End-to-end composition" below.

## Summary — read this first

**The problem.** On the TPV104 200 m benchmark at order 3, our solver spends 245.9 µs of core
time per element update while SeisSol spends 14.4 µs — a 17× gap that dwarfs everything else we
have measured. (Both numbers are from pre-rupture windows; see the caveat below.) A live profile
shows where our time goes: about 59 % is the face-flux stage, only 27 % is the time-expansion
predictor, and the fault is ~1 % in that window. Three specific taxes dominate: every basis
evaluation on a tetrahedron secretly performs a dense triangular solve (~15 % of the step, in
the face loop); the face loop re-derives geometry, rotations, and basis tables from scratch
every step and allocates temporaries per quadrature point; and the predictor sweeps the whole
state vector ~40 times per step instead of finishing each element while its data is in cache.

**What the prototype already told us (rev-3 change).** Before committing implementation, we
built and ran a standalone benchmark of the candidate kernel structures against the real MFEM
library. Three results redirect this plan: (1) the current kernel's *arithmetic* structure is
NOT the bottleneck — a faithful replica runs at 5.2 µs/element on a well-fed core (48 % of
peak), ~47× cheaper than the 245.9 µs production figure; (2) MFEM's built-in batched
linear-algebra CPU backend is **2.2× slower** than the current structure and is eliminated as
the mechanism; (3) the surviving case for restructuring is **memory traffic** — the current
structure streams ~8× the DRAM bytes of a fused-tile design, which a bandwidth-starved
128-rank node cannot absorb even though a laptop can. The roofline confirms: with operators
kept resident across the time-expansion levels the target has 2.2–2.8× bandwidth headroom;
without that fusion the target is unreachable. So Phase 2 is re-founded as a
traffic-reduction program with hand-tiled fused loops (optionally small-dgemm), gated by a
Rome-side rerun of the same benchmark.

**The fix.** A four-phase, measure-first program. First, measure-and-activate: turn on the
face cache we already built but never enabled in the benchmark, run the same microbenchmark on
the cluster CPU, decompose a rupture window, and re-baseline. Second, face-stage surgery:
evaluate basis, rotation, and flux tables once at setup and reuse them every step. Third,
rebuild the predictor and volume stages as cache-tiled fused loops over contiguous operator
packs (the traffic fix), per LTS cluster. Fourth, re-run the head-to-head and re-measure.

**End-to-end composition (why this program alone does NOT approach SeisSol).** Kernel success
at the full 6× target still leaves the LTS leg at ~2677 s/sim-s — **11.2× behind SeisSol
as-run** — because the untouched communication wait (2341 s/sim-s) then dominates. Only the
composed state (this program + the companion comm program) reaches **~2.4–3.9× of SeisSol
as-run** (~585–921 s/sim-s), and that composed state is comm-bound. The full scenario table
lives in the evidence record; Phase 4 publishes the measured version of it.

**Expected outcome.** Phase 0 (activation) is bounded by Amdahl at roughly **1.2–1.4×** — real
but modest; it exists to produce the trustworthy baseline ("B0"), not the win. The program
target is **≤ 40 µs·core per element update on the B0 measurement protocol** (≥ 6× vs today's
245.9), which would put MFEM within ~3× of SeisSol's kernels at matched order. Every
downstream numeric gate in this plan is **provisional until B0 exists** and is re-derived in a
mandatory Phase-0 "gate-reset memo" — the current numbers come from an Apple-silicon profile
and are known not to transfer 1:1 to the cluster. We do **not** promise matching SeisSol's
14.4 µs, and if the roofline shows the batched stage is bandwidth-bound on the cluster CPUs,
the honest end target degrades to ~3–3.5× and the plan says so at that gate, not after.

**Main tradeoff / biggest risk.** Most of these changes re-associate floating-point sums, so
they are *round-off-level*, not byte-exact. The project's byte-exact contract is preserved the
way the earlier caching flag did it: every change is opt-in behind a flag, defaults stay
bit-identical, acceptance is a documented tolerance (≤ 1e-12) plus physical metrics. The
performance risk is that data movement (gather/scatter for batching) eats the arithmetic win —
this is exactly why the earlier p1-era GEMM idea was rejected — so every phase carries a
measured stop-and-reassess gate, and the batching design is cache-tiled specifically to avoid
repeating that rejection.

**What this does NOT do.** No fault/friction optimization *yet* — its ~1 % share was measured
pre-rupture and the decision is deferred until Phase 0 measures a rupture window. No
communication work (the exposed-wait problem is a separate program; see "Relationship to the
comm plan"). No default-behavior changes for the TPV*/BP5 standalone drivers. No fp32. No
external code generation. **Version 1 accelerates the TPV104-class path only** (scalar
material, pure upwind); the production-SAFS extension is an explicitly named follow-on phase.

## How to read this plan

- The **Summary** above is the whole idea.
- Each **Phase** opens with a one-sentence goal; **Detailed Requirements** are the
  implementation contract; **Interfaces** give the exact signatures.
- **Phase close-out protocol:** each phase ends with `/code-review` against this plan, findings
  recorded as `document/kernel_dev/REVIEW_<phase>_impl_<date>.md`, `/code-fix` applied, and the
  phase's acceptance boxes ticked only after the review returns no CRITICAL.
- The **Glossary** defines every shorthand used below.

## Glossary

| Label | Plain-language meaning |
| ----- | ---------------------- |
| per-update µs·core | core-seconds per element per timestep × 1e6. Expanse pre-rupture reference: MFEM 245.9, SeisSol 14.4 |
| B0 | the Phase-0 Expanse baseline: per-update µs·core + per-subphase Caliper split, measured with all landed levers ON, on the pinned protocol (mesh, window, flags, binary). Denominator for every later gate |
| gate-reset memo | mandatory Phase-0 exit artifact: re-derives every downstream numeric gate from the measured B0 subphase table |
| CK recursion | the ADER Cauchy–Kovalevskaya predictor: repeated spatial-derivative applications building time-Taylor terms D(k) |
| D(k) | the k-th time-derivative block per element. LTS providers must write raw D(k) to `dk_retain[(slot·order+k)·block+j]` for elements with `retain_slot_of_elem[e] ≥ 0` — the seam bit-copy contract |
| levers | landed opt-in flags: `--deriv-cache` (per-element D_d=M⁻¹K_d matrices), `--shared-ck-recursion` (single CK), `--face-cache` (interior-face geometry/shape cache) |
| Vandermonde-LU tax | MFEM's L2 tet elements evaluate CalcShape/CalcDShape by LU back-substitution against a stored factorization (`fe_l2.cpp`, `Ti` = DenseMatrixInverse) — a dense solve per basis evaluation. Measured ≈ 15 % of the step (face loop; local np=4 profile, levers ON, no --face-cache) |
| byte-exact contract | TPV*/BP5 standalone binaries stay bit-identical; defaults (OnTheFly, two-recursion, no face cache) must not flip |
| round-off lever convention | opt-in flag, legacy path untouched and default, unit parity ≤ 1e-12, end-to-end acceptance on physical metrics — the `--deriv-cache` precedent |
| single-cluster gate | LTS invariant: cluster predictor over all elements byte-identical to the GTS call at fixed kernel mode |
| R-001 / R-003 / R-004 | guardrails: substep path unconditional on every rank (collective); the two factorial-weight accumulations stay distinct; per-element caches carry fail-loud memory budgets with fallback |
| face-cache exclusivity (REVIEW R-001) | `SetUseFaceCache` and `UsePrecomputedFaceFluxes` are mutually exclusive (`wave_operator.inl:871/893`) — distinct from R-1203 below |
| R-1203 | the MixedFluxMode ↔ `UsePrecomputedFaceFluxes` mutual-exclusion guard (`wave_operator.hpp:74`) |
| OPT-* IDs | prior optimization verdicts (adopted/rejected, with reasoning) in `document/code_optimization_dev/action_plan_2026-06-24.md` |
| NATIVE batched backend | `BatchedLinAlg::Get(BatchedLinAlg::NATIVE)` — `mfem::forall` + fixed-order `kernels::AddMult` per matrix slice; bitwise deterministic, runs on the CPU-only builds |

## Technical Overview

The hot path lives in `dynamic/wave_operator.inl` (predictor `ApplySpatialDerivative` +
`ComputeADERSubStepStatesAndIntegral*`; corrector `ComputeADERFaceFluxRHS` →
`ProcessADERFaceToRHS_` per face, `ComputeVolumeRHS`, `ApplyMassInverse`) driven from
`drivers/spatial_dyn_driver.cpp`. Measured split at p3/ADER-O4 with levers ON and
`--face-cache` OFF (local profile, np=4, 1000 m mesh, pre-rupture window): face corrector 59 %
(Vandermonde-LU tax ≈ 15 %, CalcShape/geometry ≈ 6 %, Godunov rotations ≈ 9 % subtree,
allocation churn ≈ 7 %), predictor 27 % (88 % in the cached matvec, which traverses
column-major matrices transposed at stride 20), volume 9 %, mass-inverse 1.4 %, fault 1.3 %
(pre-rupture lower bound). The state vector is component-major (9 strided streams per
element); per-element operator caches are heap-scattered `DenseMatrix` triples. MFEM 4.9.1 in
this checkout provides `DenseTensor`, `BatchedLinAlg` (NATIVE/GPU_BLAS/MAGMA), `kernels::`,
and `mfem::forall`; the GPU worktree holds validated (uncommitted) element-local `forall`
kernels for the same operations. PA/libCEED is a dead end here (simplex DG, no BilinearForm)
per `document/gpu_dev/PA_and_custom_function_gpu_analysis.md`.

## Applicability (which decks this accelerates)

| Path | Phase 0 | Phase 1 | Phase 2 | Phase 3 |
|---|---|---|---|---|
| TPV104/TPV102-class: scalar material, `mixed_flux=none`, upwind | ✔ | ✔ | ✔ | ✔ |
| SAFS production: `BimaterialWaveOperator` (CVM) and/or mixed flux | ✔ (measurement only) | ✖ parse-rejected | ✖ parse-rejected | ✖ |

**The SAFS extension is Phase 5** (explicitly deferred, not forgotten): port the face tables to
the bimaterial per-face matrix layout and the batched kernels to the bimaterial operator cache,
with the bimaterial parity gate (`test_bimaterial_deriv_cache_parity` pattern). Trigger: Phase
1+2 green on TPV104 and the SAFS QD/dynamic program requesting it. Until then, every new flag
is parse-time rejected on the unsupported paths (tested — see Testing Strategy).

## Relationship to the comm plan (rev 3: re-sequenced from the numbers)

The comm program now EXISTS as a plan: `document/comm_dev/PLAN_lts_comm_reduction_2026-07-22.md`
(built from a full exchange inventory of the tick loop). The numbers argue for **co-first-class
sequencing, not kernel-first**: comm removes the larger absolute term (2341 vs the kernel's
~1687 s/sim-s), its Phases 1–2 are **bitwise-identical by construction** (no tolerance
machinery), and its payoff arithmetic depends on message counts already measured — not on any
unmeasured throughput rate. Concretely:

- **Shared Phase 0:** the comm plan's Waitall wire-vs-skew decomposition rides in the SAME
  Expanse job as this plan's B0 (one submission, both baselines).
- **Then in parallel:** comm Phases 1–2 (per-tick merge — touches the LTS steppers + seam
  corrector) and kernel Phase 1 (face tables — touches the face corrector). File overlap is
  limited to `wave_operator.inl` in disjoint regions; both are flag-gated; whichever lands
  second rebases.
- **Revisit trigger:** if the Rome-side spike (Phase 0) shows < 2× realizable kernel headroom,
  kernel Phase 2 is deprioritized below comm Phase 3 (overlap) — the composed end state is
  comm-bound anyway.
- The B0 metric is per-update *compute* µs·core (function-level Caliper region + perf self-time sample, comm excluded) so comm
  changes do not confound kernel gates; if a comm change lands mid-program, B0 is re-measured
  before the next kernel gate is evaluated.

## Constraints

- **No-touch:** `bp5/ bp1/ bp2/ domain/ fault/ solver/ friction/dieterich_ruina.hpp`. Editable:
  `dynamic/`, `drivers/spatial_dyn_driver.cpp`, `jobs/`, `tests/`, **plus `build_expanse.sh`
  (repo root) for the Phase-0 BLAS item only**.
- **Byte-exact contract** (see Glossary); all new kernels opt-in; defaults untouched; TPV*/BP5
  binaries never execute the new paths.
- **Round-off lever convention** for every non-bit-exact change.
- **LTS couplings:** single-cluster gate; D(k) retention to `dk_retain` exactly as the scalar
  retain lambda does it (`wave_operator.inl:1990-2009`), bit-identically, for
  `retain_slot_of_elem[e] ≥ 0` only; retained-CK bitwise rank-independence (no rank-dependent
  reduction order; no divergence between retain pass and main pass); R-001; R-003; the seam
  D(k) epoch guard.
- **Do not hand-optimize the OnTheFly path** — byte-exact reference.
- **Never mix structural refactor and numerics in one commit; `make clean` before verification
  builds** (stale-`.inl`-ABI trap).
- Local validation = unit/parity tests + np≤10 smokes on the 1000 m fixture (**the 200 m
  production mesh cannot run locally with full per-face tables — the memory guard will trip;
  see Phase 1 budget**); Expanse for throughput; **ask before any cluster submission**.
- R-004-style fail-loud memory budget for every new table, with the per-rank estimate printed.

---

## Phase 0: Measure and activate (configuration + measurement only) — IMPLEMENTED

**Status (2026-07-22):** implemented and locally validated; the Expanse job is prepped and
**awaiting user approval to submit** (standing rule). Deliverables committed:
`jobs/kernel_dev/run_p0_baseline_expanse.sbatch` (5 legs), `tests/bench/run_bench.sh`,
`document/kernel_dev/B0_gate_reset_2026-07-22.md` (the exit memo), the `build_expanse.sh` BLAS
fix, and the `--face-cache` default in the speed sbatch.

**In one sentence:** Turn on the already-landed face cache, measure the hot path by *sampling*,
measure a rupture window, settle the BLAS question, and produce the B0 baseline + gate-reset memo.

### Goal
A trustworthy Expanse baseline (B0) with per-subphase attribution, including a rupture-phase
measurement, before any kernel code is written.

### Instrumentation methodology (rev 3.1 change — no numerics file touched)
The rev-2 plan called for `MFEM_PERF_SCOPE` subregions inside the per-QP loops. **Rejected during
implementation:** a Caliper push/pop per face (4.9 M faces × ~250 steps) would add seconds of
overhead and *perturb the very timing it measures* — a self-inflicted measurement artifact. The
sub-function split (LU / CalcShape / rotations) is a *sampling* problem, so Phase 0 uses **`perf
record -g` on rank 0** (the Expanse equivalent of the local macOS `sample` that already produced
the 59 % / LU-15 % / CalcShape-6 % split), plus the existing function-level Caliper region report.
**Consequence: `dynamic/wave_operator.inl` is NOT edited in Phase 0** — the safest possible Phase 0.

### Files created / modified
- `jobs/kernel_dev/run_p0_baseline_expanse.sbatch` (NEW) — the 5-leg B0 job:
  (1) BLAS check (grep the deployed `config.mk`), (2) B0 GTS `--tfinal 1.5` with Caliper + `perf`,
  windows read from the step log at `[0,0.5]` and `[1.0,1.5]`, (3) 64-vs-128 ranks/node layout A/B,
  (4) the kernel bench spike (single + full-occupancy — Phase-2 go/no-go), (5) short LTS leg with
  Caliper `profile.mpi` for the `MPI_Waitall` Max/Avg skew signal (comm-plan Phase 0).
- `tests/bench/run_bench.sh` (NEW) — toolchain-portable bench compile (via `make -n`) + single/
  full-occupancy run.
- `jobs/lts_phase5/tpv104_200m_lts_speed_expanse/run_tpv104_200m_lts_speed_expanse.sbatch` —
  `P5_DERIV_OPT` default now includes `--face-cache`.
- `build_expanse.sh` — **BLAS bug fix (real finding):** the PETSc-openblas wiring only scanned
  `${PETSC_PREFIX}/lib`, so a `lib64` install left `MFEM_USE_LAPACK=NO` (the slow internal LU path,
  a candidate contributor to the 245.9 µs). Now scans `lib` and `lib64` and warns if unresolved.
- `document/kernel_dev/B0_gate_reset_2026-07-22.md` (NEW) — the exit memo template with formulas.

### Detailed Requirements
1. Verify `--face-cache` composes with `--deriv-cache --shared-ck-recursion` (its only
   exclusivity is `UsePrecomputedFaceFluxes` — the face-cache REVIEW R-001 check at
   `wave_operator.inl:871/893`; R-1203 proper is the mixed-flux guard) and is exercised by
   `test-face-geom-cache`.
2. **BLAS (measure-only, no promised gain):** record whether the deployed Expanse binary has
   `MFEM_USE_LAPACK`; if absent, rebuild linking a verified **LP64** BLAS/LAPACK (the
   PETSc-downloaded OpenBLAS; never the ILP64 `cpu/0.17.3b` module OpenBLAS) and A/B a short
   GTS leg. **The relink is a round-off-level build change**: acceptance includes V_max(t) /
   station sanity of the rebuilt binary vs the 52344266 GTS reference under the round-off
   convention. **B0 is measured on the post-rebuild binary** (pin the binary hash in the memo).
   Note: the local profile's LU cost was measured WITH Apple's LAPACK present — BLAS presence
   does not remove the Vandermonde-LU tax, and it becomes moot once Phase 1 removes the solves
   from the hot loop. The local-vs-Expanse per-update gap (58.4 vs 245.9) is expected to be
   largely hardware; do not book it as recoverable.
3. **B0 runs (needs user approval to submit):** GTS leg, levers + `--face-cache`, both windows,
   function-level Caliper region report + `perf record -F 99 -g` self-time sampling on rank 0 (self-time symbol split; intra-function inlined blocks fold into the parent at -O2 — cross-check the local `sample` profile). Deliverables: per-update µs·core (pre-rupture AND rupture windows),
   subphase table, `MPI_Waitall`-free compute shares, HW-FLOP estimate if obtainable
   (`perf`-class counters via a short `LIKWID`/`perf stat` srun wrapper if available on
   Expanse; else the Caliper-derived arithmetic-intensity estimate) — the **roofline check**
   the p1-era GEMM rejection demanded before any batching verdict.
4. **Gate-reset memo** (`document/kernel_dev/B0_gate_reset_<date>.md`): re-derive every
   downstream numeric gate (Phase 1/2/3 targets) from the measured B0 subphase table; record
   the rupture-window fault share and the go/defer decision on fault work.
5. Local: re-run the `sample` profile with `--face-cache` ON (expects the interior-face
   LU/CalcShape block to shrink; fault faces unaffected — they are outside the face cache).
6. **Rome-side kernel spike (review R-402; the local half already ran):** compile and run
   `tests/bench/bench_ader_kernel_variants.cpp` on an Expanse compute node — once single-core,
   once at full 128-rank occupancy (mpirun the bench on every core simultaneously to expose the
   ~2 GB/s/core bandwidth regime). Deliverables: A-vs-B ratio and absolute µs/elem under
   contention (locally: A=5.2 µs, B=4.2, `BatchedLinAlg` C=11.2 — eliminated; the traffic case
   predicts B's edge WIDENS under contention). **Go/no-go for kernel Phase 2:** proceed only if
   contended B beats contended A by ≥ 2×; else Phase 2 is deprioritized below comm Phase 3 per
   the revisit trigger.
7. **Node-layout A/B (review R-409):** one extra GTS leg at 64 ranks/node (same node count,
   half-populated) inside the B0 job; report per-update µs·core at 128 vs 64 ranks/node. This
   separates per-core bandwidth starvation from kernel structure in the 245.9 figure — if 64
   ranks/node recovers a large factor, the traffic case strengthens and an OpenMP-hybrid layout
   becomes a candidate lever (SeisSol runs 8×15 for a reason).
8. The companion comm plan's Phase 0 (Waitall wire-vs-skew decomposition) rides in this same
   job — one submission, both baselines.
9. Expectation setting (honest): activation is Amdahl-bounded at ~**1.2–1.4×** (face-cache
   addresses ≤ the interior-face part of the ~21 % LU+shape share plus some churn). Anything
   more is a bonus, not a plan.

### Interfaces
No new public API and no source edit to `dynamic/` (see the methodology note above). The
measurement contract is the 5-leg sbatch + `run_bench.sh <out_dir>` + the gate-reset memo.

### Rollback
Revert the sbatch default and keep both binaries deployed (pre/post-rebuild paths recorded in
the memo); instrumentation is compiled out without Caliper.

### Acceptance Criteria
- [ ] B0 table (both windows) + gate-reset memo exist; downstream gates re-derived.
- [ ] BLAS status documented; rebuilt binary (if any) passes the round-off physics sanity.
- [ ] Rupture-window fault share measured; fault-work decision recorded.
- [ ] `make test` green; no diffs under no-touch dirs.

### Dependencies
Depends on: nothing (user approval for the Expanse submissions). Required by: all later phases.

---

## Phase 1: Face-corrector surgery (the measured 59 %)

**In one sentence:** Make the per-face, per-QP work table-driven — basis, geometry, rotations,
and flux operators evaluated once at setup and applied every step with zero allocation — in
two sub-flags: a bit-exact table variant and a round-off contraction/fold variant.

### Goal
Remove the Vandermonde-LU tax, per-QP basis/geometry recomputation, per-QP heap allocation, and
per-QP hash lookups from `ProcessADERFaceToRHS_` and the fault corrector's rotation rebuilds.

### Files to Create
- `dynamic/face_kernel_tables.hpp` — table container + builder (owned by `dynamic/`; any BP5
  resemblance is pattern-copied, never shared).
- `tests/unit/test_face_kernel_tables.cpp`.

### Files to Modify
- `dynamic/wave_operator.{hpp,inl}` — flagged face-corrector paths consuming the tables;
  preallocated workspace members replacing in-loop temporaries.
- `drivers/spatial_dyn_driver.cpp` — flags `--face-tables` (bit-exact variant) and
  `--face-tables-fold` (round-off variant; implies the first). Parse-time rules: supersede
  `--face-cache` when both given (warn); mutually exclusive with `UsePrecomputedFaceFluxes`
  (same style as the face-cache REVIEW R-001 check); rejected when `mixed_flux != none`;
  rejected on `BimaterialWaveOperator` (Phase-5 deferral).
- `Makefile` — test targets.

### Detailed Requirements
1. **Table design and memory (corrected per review R-201).** Naive per-face trace tables cost
   `(2·20·12 + 4·81)·8 B ≈ 6.4 KB/face` → **~31.5 GB at np=1 on the 200 m mesh (4.9 M faces)**
   — infeasible locally, ~123 MB/rank at np=256. Therefore the design DE-DUPLICATES:
   - **Trace-basis tables are reference-level, not per-face**: `shape = CalcShape(Loc.Transform(ip))`
     depends only on the (local-face-id, orientation) reference embedding — a fixed catalog of
     ≤ 4 faces × orientations ≈ dozens of small `(nqp × ndof)` matrices TOTAL, built once. The
     LU tax is paid once per catalog entry, not per face.
   - **Per-face storage** is only: unit normal + area weights (constant per planar face),
     optionally the rotation pair `T/T⁻¹` (2·81 doubles) or — fold variant — the pre-rotated
     `A±_rot = T·A±·T⁻¹` pair (2·81 doubles): **≈ 1.3–2.6 KB/face → ~6.4–12.7 GB at np=1,
     25–50 MB/rank at np=256.** R-004 fail-loud budget guard prints the per-rank estimate; the
     guard WILL trip on the production mesh locally — local testing uses the 1000 m fixture
     (documented in Constraints).
2. **Bit-exact variant (`--face-tables`)**: keeps the legacy per-QP scalar loop shape exactly
   (the `ComputeADERFaceFluxRHS_CachedInterior_` precedent) — table lookups replace
   `CalcShape`/`CalcOrtho`/`GetFaceElementTransformations` only; accumulation order into `rhs`
   unchanged. Candidate bit-exact because the replaced values are step-invariant
   recomputations; parity test asserts it (fallback claim: ≤ 1e-12 with the discrepancy
   explained, per the face-cache precedent which documents "≤1e-12, in practice
   bit-identical").
3. **Fold variant (`--face-tables-fold`)**: the `(nqp×ndof)ᵀ`-shaped contraction rewrite with
   `kernels::` primitives AND the `T·A±·T⁻¹` fold — both re-associate sums → round-off lever
   convention (≤ 1e-12 + physical metrics).
4. **Fault-face rotations**: build a **parallel** per-(fb_idx, qp) rotation table owned by
   `face_kernel_tables.hpp`, indexed identically to `FaultBasisData::qp_data` —
   `fault/fault_basis.hpp` is NOT touched. Keyed on the canonical sign-fixed frame
   (`qpd.sign_flipped` semantics), never on rank-local `elem1_on_plus` — the
   OPT-FACEFLUX-FAULT-ROT-REUSE np>1 determinism bug is the cautionary precedent; np=2/np=10
   fault byte gates are mandatory.
5. Loop hygiene (both variants): hoist `fault_face_dof_offset_.find` /
   `LookupInteriorFaultBasisIndex` per face; eliminate `Vector shape1(ndof)` /
   `shape2.SetSize` per-QP constructions via workspace members.
6. LTS: the seam corrector uses the same per-face primitives; seam gates
   (`seas_test_lts_mpi_seam*`) must pass with each flag.

### Interfaces
```cpp
// dynamic/face_kernel_tables.hpp
struct FaceKernelTables
{
   // Reference catalog: trace basis at face QPs per (local_face_id, orientation).
   // cat_index = ref_catalog_index_[f] per mesh face.
   std::vector<DenseMatrix> ref_trace;        // catalog entries: (nqp x ndof)
   Array<int> ref_index_self, ref_index_nbr;  // per interior face -> catalog id
   // Per-face geometry (planar): unit normal, |J_F| weights.
   Vector face_normal;                        // 3 * n_faces
   Vector face_weight;                        // nqp * n_faces (or 1 if constant)
   // Fold variant only:
   Vector Aplus_rot, Aminus_rot;              // 81 * n_faces each
   // Fault path: per-(fb_idx, qp) canonical rotation pair, parallel to
   // FaultBasisData::qp_data (fault/ untouched).
   Vector fault_T_can, fault_Tinv_can;        // 81 * n_fault_qps each
   long long BytesUsed() const;
};

// Builder replicates the legacy per-face rule selection EXACTLY (review R-406):
//   IntRules.Get(geom, face_is_fault ? FaultFaceQuadDegree() : 2*order_)
// (wave_operator.inl:5524-5526) — fault faces use the over-integration branch,
// so fault_T_can/fault_Tinv_can are sized by the FAULT rule's nqp, not 12.
// Guards (review R-407): MFEM_VERIFY(mesh.Conforming()) — nonconforming meshes
// parse-time rejected (none in scope; fail loud); boundary faces (absorbing/
// free-surface) get their own catalog entries or stay on the legacy path in v1
// (decided at implementation; either way stated in the table struct docs).
void BuildFaceKernelTables(const WaveOperator<MeshType>& wave,
                           bool build_fold_tables,
                           FaceKernelTables& out);   // R-004 guard inside

// wave_operator.hpp additions
void SetUseFaceTables(bool enable, bool fold);       // exclusivity checks here
```

### Edge Cases
- **Curved meshes (`mesh.order ≥ 2`)**: normals become QP-dependent → rotation storage goes
  per-QP (budget multiplies by nqp; guard covers it). A `mesh.order=2` tiny-fixture parity test
  is REQUIRED (no silent affine assumption — the OPT-VOL-AFFINE-HOIST lesson).
- Shared (rank-boundary) faces: tables built from face-neighbor data; np>1 parity test.
- Mixed flux / bimaterial: parse-time rejection (Applicability).

### Rollback
Flags off = legacy path byte-identical (legacy code untouched and compiled-in permanently;
`--face-cache` retained as an independent fallback until the Phase-4 disposition decision).

### Acceptance Criteria (numeric values provisional until the B0 gate-reset memo)
- [ ] Local profile with `--face-tables-fold`: LU/CalcShape/GrowSize symbols ~absent from the
      face subtree.
- [ ] Parity: `--face-tables` bit-exact (or ≤ 1e-12 with documented cause); fold ≤ 1e-12;
      fault np=2 byte gate; seam gates green.
- [ ] Expanse vs B0 — **share-based gate, protocol-pinned (review R-405)**: measured on the B0
      protocol + Phase-1 flags ONLY (`--kernel-batch` OFF, no comm changes since the B0 in the
      denominator), pre-rupture window (rupture window reported alongside): the face-stage
      subphase share drops to ≤ 0.5× its B0 share, AND per-update improves ≥ 1.2× vs B0.
      Stop-and-reassess if the share gate passes but per-update moves < 1.1×.
- [ ] Flag-matrix parse-time tests green; `make test` green; no-touch dirs clean.

### Dependencies
Depends on: Phase 0 (B0 + memo). Required by: Phase 3. Parallel-safe with Phase 2 (separate
flags, separate subphases).

---

## Phase 2: Tiled-fused element-local predictor + volume — a TRAFFIC program (rev 3)

**In one sentence:** Replace the whole-vector CK sweeps with cache-tiled fused loops over
contiguous operator packs — cutting DRAM traffic ~8× per element — because the microbenchmark
proved the win is bandwidth, not arithmetic, and eliminated `BatchedLinAlg` NATIVE as the
mechanism (2.2× slower than the current structure).

### Goal
Element-local fusion of the CK recursion over cache-sized element tiles inside cluster
batches: operators read ONCE per element per step (levels fused — the roofline shows this is
load-bearing: fused floor 14–18 µs vs non-fused 45–49 µs at 2 GB/s/core), eliminating the
~40 full-vector passes. Mechanism per the measured spike: **hand-tiled fixed-stride fused
loops (bench variant B), optionally per-element dgemm (variant D) where BLAS is linked** —
NOT `BatchedLinAlg` NATIVE (bench variant C, measured 0.46× — its CPU path is a serial
generic-kernel loop). `DenseTensor` is retained as the contiguous storage container only.
Precondition: the Phase-0 Rome-side contended spike gate (≥ 2× contended B-vs-A), honoring
the OPT-DERIV-GEMM-BATCH revisit clause with measured data.

### Files to Create
- `dynamic/ader_batched_kernels.hpp`, `tests/unit/test_ader_batched_parity.cpp`.

### Files to Modify
- `dynamic/wave_operator.{hpp,inl}` — operator repack at `SetDerivMode(Cached)` into
  per-direction `DenseTensor`s (cluster-order slices when LTS on); flagged batched
  implementations of `ComputeADERSubStepStatesAndIntegral{,Cluster}` and `ComputeVolumeRHS`.
- `drivers/spatial_dyn_driver.cpp` — flag `--kernel-batch` (parse-time: requires
  `--deriv-cache --shared-ck-recursion`; rejected on bimaterial).
- `Makefile`.

### Detailed Requirements
1. **Repack**: `elem_deriv_op_[e][d]` → `DenseTensor deriv_ops_[3]` (ndof, ndof, ne) packed in
   batch order + `elem→slot` map; same for volume ops and `elem_mass_inv_`. Setup-time only;
   R-004 guard (same bytes as today, now contiguous).
2. **Tile structure (review R-203): batch granularity is a cache-blocked tile of E ≈ 32–128
   elements** (panels sized ≤ L2), iterated *within* an LTS cluster (or the whole rank under
   GTS). Per tile and per CK level: gather the state panel once, run the three per-direction
   batched GEMMs, then a fused per-element loop doing the 9×9 Jacobian combination AND both
   factorial accumulations (τᵏ/k! nodal, dtᵏ⁺¹/(k+1)! integral — two accumulators, R-003)
   before the tile's next level — the per-element operator stream is what must not be
   re-streamed per level from DRAM.
3. **D(k) retention (review R-103)**: at each level, for elements with
   `retain_slot_of_elem[e] ≥ 0`, scatter D(k) additionally to
   `dk_retain[(slot·order+k)·block+j]`, **bit-identically to the scalar retain lambda**
   (`wave_operator.inl:1990-2009`). The ping-pong scratch needs no write-back. D(k) stays
   materialized per element (never fused away) — the seam contract.
4. **Determinism/mechanism (rev 3, from the spike)**: the tile kernel is a hand-written
   fixed-stride fused loop (bench variant B — checksum bit-identical to the legacy order by
   construction when the per-element summation order is preserved; else ≤ 1e-12 lever
   convention). Per-element dgemm (variant D) is an optional flag-selected alternative ONLY
   where an LP64 BLAS is verified (Expanse hazard) and NEVER on the seam-contract path unless
   its bitwise reproducibility is proven (BLAS order is implementation-defined). Bitwise
   rank-independence of the retained-D(k) pass is asserted by the existing epoch/parity tests.
   `Xb`/`Yb` are `Vector` panels viewed as (ndof, NUM_STATE, E).
5. Precision: summation order changes (row-dot → column-AXPY) → round-off lever convention;
   parity ≤ 1e-12 vs Cached on unit fixtures across orders 2/3/4, np1/np2; bit-exact checks
   only WITHIN batched mode (retain-vs-main D(k), two-run reproducibility at np=2).
6. Prior-art gate: this phase IS the sanctioned revisit of OPT-DERIV-GEMM-BATCH (rejected at
   p1; revisit clause = compute-bound roofline at higher order). The Phase-0 roofline is a
   PRECONDITION: if it shows the cached matvec bandwidth-bound on Rome at p3, the expectation
   resets per the memo before implementation proceeds.
7. `MFEM_VERIFY` homogeneity check hoists out of the inner loop; `ApplyJacobianPerDOF`'s ~30
   whole-vector sweeps per level disappear into the fused tile body.

### Interfaces
```cpp
// dynamic/ader_batched_kernels.hpp
struct ADERBatchScratch      // per-thread/tile workspace, sized once
{
   Vector Xb, Yb, Dk;        // (ndof * NUM_STATE * E_tile) panels
   void EnsureSize(int ndof, int ns, int e_tile);
};

// One CK level over one tile: batched 3-dir derivative + fused Jacobian +
// factorial accumulations. elem_slots = tile's slots into the DenseTensors;
// retain handled inside per requirement 3.
void BatchedCKLevelTile(const DenseTensor deriv_ops[3],
                        const real_t* A_jac /* 9x9 per dir combination */,
                        const int* elem_ids, const int* elem_slots, int E,
                        int level_k, real_t dt, int order,
                        const int* retain_slot_of_elem, real_t* dk_retain,
                        Vector& Q_state, Vector& Q_per_node, Vector& I_accum,
                        ADERBatchScratch& ws);

// Dispatch hook inside ComputeADERSubStepStatesAndIntegral{,Cluster}:
//   if (kernel_batch_) { for (tile : tiles(elems)) BatchedCKLevelTile(...); }
//   else { legacy path unchanged }
```

### Edge Cases
- Fault-provider and seam-coarse retain passes run the SAME tile kernel (bit-identical inputs
  ⇒ bit-identical D(k)); epoch-guard test extended to batched mode.
- Bimaterial: parse-time rejected (Phase 5).
- Empty clusters / fault-free ranks: zero tiles is a no-op; surrounding collective structure
  untouched (R-001).
- **Ragged tiles (review R-408)**: per-(rank, cluster) populations below `E_min` (default 8,
  flag-tunable — the λ-wiggle's 888-cell finest cluster is ~3.5 elems/rank at np=256) fall
  back to the legacy per-element path for that (rank, cluster) — correctness-neutral because
  every batched operation is per-element independent; the single-cluster byte gate is
  tile-geometry-safe by construction (no cross-element reductions).

### Rollback
Flag off = legacy Cached path untouched.

### Acceptance Criteria (provisional until the memo)
- [ ] Unit parity ≤ 1e-12 (orders 2/3/4, np1/np2); retain-vs-main D(k) bitwise; np=2 two-run
      bitwise reproducibility.
- [ ] LTS gates green with `--kernel-batch`: single-cluster==GTS byte gate AT batched mode,
      seam 58/58 + 111/111, fault np=2/np=10 byte-identity.
- [ ] **Local** microbenchmark: predictor+volume ≥ 2.5× vs the LOCAL levers-on cached path
      (Apple-silicon is compute-bound — this validates the arithmetic).
- [ ] **Expanse**: predictor+volume subphase ≥ 1.5–2× vs B0 (roofline-informed expectation;
      exact gate from the memo); gather+scatter ≤ 25 % of the batched kernel time — else STOP
      and go to Phase 3's layout decision first.
- [ ] Cumulative Expanse per-update ≥ 2.5× vs B0 when combined with Phase 1 (composition of
      the two sub-gates; NOT 3× — review R-101F).

### Dependencies
Depends on: Phase 0 (roofline + memo). Parallel-safe with Phase 1. Required by: Phase 3.

---

## Phase 3: Layout end-game + device convergence

**In one sentence:** With Phase-2's measured gather cost in hand, decide the state-layout
question, port the validated GPU-worktree `forall` kernels forward, and make the batched path
the single CPU/GPU code path.

### Goal
Close the remaining data-movement gap if measurement says it pays, and converge the two kernel
lineages (this plan's batched path; the GPU worktree's forall kernels).

### Files to Create/Modify (provisional pending Phase-2 data)
- `dynamic/wave_operator.{hpp,inl}` — state-panel members if option (b) is chosen.
- `dynamic/ader_batched_kernels.hpp` — device-annotated tile kernels.
- Ported-by-file from the GPU worktree: `dynamic/gpu/wave_device_kernels.hpp`,
  `tests/unit/test_wave_device_kernels.cpp` (re-applied onto current HEAD — the worktree's
  `wave_operator.inl` predates LTS Phases 0–5; NEVER a blind merge).

### Detailed Requirements
1. Layout decision by measurement, in increasing invasiveness: (a) keep component-major +
   per-tile gather (status quo); (b) persistent per-cluster panel mirrors of Q synced at
   sync-points; (c) full relayout — prior rejection OPT-SOA-LAYOUT (60+ sites, halo hazards)
   stands unless evidence is overwhelming; if chosen, lands as its own structural-only commit
   series.
2. GPU worktree port: per-file onto HEAD, byte-gated on `cpu`+`debug` backends per the
   worktree's validated methodology; or an explicit decision to retire it in favor of the
   batched path (recorded in the phase review).
3. Backend policy: NATIVE everywhere the seam bit-copy contract applies; vendor batched GEMM
   only on paths proven outside it.
4. Optional measured A/B: `MFEM_USE_SIMD=YES` rebuild.

### Edge Cases
Mixed clusters (tiles straddling cluster boundaries are forbidden — tiles are cluster-local);
fault-provider elements under any relayout keep `dk_retain` layout unchanged (it is the seam
ABI); shared-face gather under option (b) respects the byNODES halo contract (R-004-adjacent
ghost-unpack hazard from the rejected OPT-SHARED-BATCH-HALO applies).

### Rollback
Option (b)/(c) land as revertible structural-only commit series; GPU port is per-file commits;
flags gate everything.

### Acceptance Criteria (set by the gate-reset memo; the plan pre-commits only the structure)
- [ ] Documented layout decision with Phase-2 gather numbers attached.
- [ ] GPU-worktree kernels merged & green on cpu/debug, or retirement decision recorded.
- [ ] Expanse cumulative per-update: **program target ≤ ~40 µs·core** (≥ 6× vs the 245.9
      pre-activation reference — i.e. ≥ 6/(B0-factor) × vs B0; the memo states the number).
      If the Phase-0 roofline showed bandwidth-bound batching, the recorded fallback target is
      ~3–3.5× vs 245.9 and THIS criterion says so explicitly at memo time.

### Dependencies
Depends on: Phases 1+2. Required by: Phase 4.

---

## Phase 4: Re-benchmark + acceptance + documentation

**In one sentence:** Re-run the Phase-5 TPV104-200m head-to-head with the new kernels, publish
the new three-way table with physics gates, and close out documentation and flag dispositions.

### Goal
Certify the program outcome on the same protocol that motivated it.

### Files to Modify
- `jobs/lts_phase5/tpv104_200m_lts_speed_expanse/` (flags), `document/lts_dev/RESULTS_*.md`
  (extended table), `miniapps/seas/CLAUDE.md` + driver `--help` (new flags + exclusivity
  rules), project memory note.

### Detailed Requirements
1. Same decks as Phase 5. **The certified configuration pins ONE exact flag line** (review
   R-410), decided at Phase-1 close-out (fold in or out, with the measured delta recorded);
   the same line goes to CLAUDE.md as the recommended spatial-driver performance flags. The
   certified LEG (GTS window vs LTS leg) is likewise pinned and stated in the table.
2. Speed readout on BOTH windows (pre-rupture and rupture — the fault share grows
   post-rupture; the certified number states its window). Table: MFEM before/after, SeisSol
   references. **The SeisSol reference is re-measured, not reused (review R-404):** pin the
   existing 52365078 run (binary hash v1.3.1-2135, deck, cadence) in the table AND add ONE
   SeisSol re-run with async output / production cadence (deck-only change) so the comparison
   brackets SeisSol's honest range (~125–239 s/sim-s) instead of freezing its
   favorable-to-MFEM as-run number.
3. Physics gates: V_max(t) vs the validated GTS reference; SCEC station diff legacy-vs-new
   kernels (round-off convention); MFEM-vs-SeisSol qualitative overlay.
4. Dispositions recorded: `--face-cache` (keep as fallback vs deprecate), fault-work decision
   (from the Phase-0 rupture measurement), default-flip question explicitly deferred to a
   dedicated plan (byte-exact contract).
5. Memory note: B0, achieved factor, gate outcomes.

### Edge Cases
Partial adoption (face-tables green but kernel-batch stop-gated): publish the partial result
honestly; the phases are independently valuable.

### Rollback
n/a (measurement + documentation).

### Acceptance Criteria
- [ ] New three-way table published with window provenance; target met or shortfall
      attributed by perf sub-function self-time / function-level Caliper region.
- [ ] All physics gates green; regression suite green; no-touch dirs clean.
- [ ] CLAUDE.md/help/memory updated; dispositions recorded.

### Dependencies
Depends on: Phases 1+2 minimum (3 optional — partial adoption edge case). Expanse submissions
need user approval.

---

## Testing Strategy

- Per phase: unit parity (bit-exact where claimed, ≤ 1e-12 otherwise) + `make test` umbrella +
  LTS gate set (`seas_test_lts_*`, seam np gates, fault np2/np10 byte-identity) + `make clean`
  rebuild before verification.
- **Flag-matrix parse-time tests** (new, Phase 1+2): `--face-tables` vs
  `UsePrecomputedFaceFluxes` exclusivity; `--face-tables` vs `mixed_flux != none`;
  `--kernel-batch` requires `--deriv-cache --shared-ck-recursion`; `--kernel-batch` vs
  bimaterial; `--face-tables` supersedes `--face-cache` (warn path).
- **Curved-mesh fixture** (`mesh.order=2` tiny mesh): face-table parity incl. per-QP rotations.
- Determinism: batched retain-vs-main D(k) bitwise; two-run bitwise reproducibility np=2.
- Performance gates are acceptance items with explicit stop-and-reassess thresholds; all
  numeric values provisional until the B0 gate-reset memo re-derives them.
- Physics: end-to-end TPV104 V_max + station comparisons per the round-off convention, on both
  measurement windows.

## Risk Assessment

| Risk (plain language) | Severity | Mitigation | Related constraint |
|---|---|---|---|
| Round-off changes leak into the byte-exact standalone drivers | High | opt-in flags; defaults untouched; no-touch dirs diff-checked every phase | byte-exact contract |
| Batched kernels break the LTS seam bit-copy | High | NATIVE-pinned backend (`Get(NATIVE)` + setup `MFEM_VERIFY`); fixed op order; `dk_retain` written bit-identically to the scalar lambda; seam + fault np-gates every phase | LTS couplings |
| Gather/scatter eats the GEMM win (p1 rejection repeats) | High | Phase-0 roofline PRECONDITION; cache-tiled batches (≤ L2 panels, fused levels); ≤ 25 % gather gate with stop clause; layout escalation only on evidence | OPT-DERIV-GEMM-BATCH revisit clause |
| Gates promised from Apple-silicon numbers don't transfer to Rome | High | every downstream gate provisional; B0 gate-reset memo is a mandatory Phase-0 exit artifact | Phase 0 |
| Fault cost grows post-rupture and erodes the certified number | Medium | rupture-window measurement in Phase 0 AND Phase 4; fault-work decision deferred, not excluded | review R-301/R-204 |
| Face-table memory (31.5 GB naive at np=1) | Medium | reference-catalog dedup (~1.3–2.6 KB/face residual); R-004 guard printing per-rank estimate; local testing pinned to the 1000 m fixture | corrected budget |
| Fault rotation cache reintroduces the np>1 determinism bug | Medium | canonical sign-fixed keying; parallel table in `dynamic/` (fault/ untouched); np2/np10 byte gates | OPT-FACEFLUX-FAULT-ROT-REUSE |
| Expanse ILP64 BLAS trap | Medium | LP64-only (PETSc OpenBLAS); relink treated as round-off change with physics sanity gate | build_expanse.sh hazard |
| GPU-worktree merge conflicts corrupt LTS code | Medium | port-by-file onto HEAD with byte gates; never blind merge; structural-only commits | CLAUDE.md no-mix rule |
| Comm plan collides in the same files | Medium | declared sequencing (kernel 0–2 first); compute-only B0 metric; B0 re-measure on interleaved landing | Relationship to the comm plan |
| Over-promising vs SeisSol | Low | absolute target ≤ ~40 µs·core with roofline-conditioned fallback (~3–3.5×) stated in advance | Summary |
