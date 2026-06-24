# Profiling Analysis — SEAS ADER spatial-dyn hot path (run 51443560)

- **Date:** 2026-06-24. **Source:** Caliper region report, ALT rate-state upwind-ADER, SAFS mesh, Expanse, 200 ranks (2N × 100r), `tfinal=1s` / 1445 macro-steps.
- **Mesh:** `safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh`, 1,164,469 type-4 (4-node **affine, straight-sided**) tets + 132,152 boundary triangles; no curved (type-11) elements. `triq` is a meshing-scheme tag, not quadratic geometry.
- **Config (verified):** P1 (fe order = 1 ⇒ `ndof_per_el_ = 4`), `ader_order = 2` (CK runs 1 inner level), `NUM_STATE = 9`, `--mixed-flux none`, `interior_flux="matrix"` (⇒ `BimaterialWaveOperator`). `~5,800` elements/rank.
- **Critical run-configuration fact:** the production srun line passes **neither** `--deriv-cache` **nor** `--shared-ck-recursion`. The default `deriv_mode_ = DerivMode::OnTheFly` (`dynamic/wave_operator.hpp:888`) and `use_shared_ck = false` (`drivers/spatial_dyn_driver.cpp:1316`) therefore governed. **This is the unoptimized baseline, not a residual-after-optimization state.**

---

## 1. Executive summary

- **Where the time goes.** The element-local spatial-derivative kernel `ApplySpatialDerivative` is the single dominant cost at **53.4 % of runtime (443.9 s)**. The corrector volume term `ComputeVolumeRHS` adds 10.1 % (83.6 s), and the face-flux correctors add ~22.6 % more (much of it MPI wait). Together the bulk DG operator apply is essentially the entire runtime.

- **The one structural finding.** `ApplySpatialDerivative` runs **TWICE per macro-step over the identical `(Q, dt, order)`** — once inside `ComputeADERSubStepStates` (substep predictor) and again inside `AdvanceADER → ComputeADERTimeIntegrated` (time-integral predictor). The fused single-recursion path (`ComputeADERSubStepStatesAndIntegral`, `dynamic/wave_operator.inl:1510`) that would feed both predictors from one Cauchy-Kovalewski (CK) recursion **is implemented and bit-exact-tested but inert** because the gating flag `--shared-ck-recursion` was not passed (`drivers/spatial_dyn_driver.cpp:459` / `:513`). The redundant pass is ~221 s, ~26.6 % of wall time, removable at zero numerical cost.

- **The load-imbalance story.** The compute regions are well balanced (`ApplySpatialDerivative` Max/Avg = 230.31/222.47 = **1.035**), because METIS balances element count and the kernel is a pure per-local-element loop. **All meaningful imbalance is concentrated in the fault-coupled work**, exposed as wait time at per-step global synchronization points: two 1-int consensus `MPI_Allreduce` guards (42.4 s + 20.2 s = **62.6 s**, min 0.16 s / max 68.2 s) and a substep halo `MPI_Waitall` (14.8 s avg / **42.1 s max**, a 100× spread). These are not reduction cost — they are fast ranks blocking on slow fault-bearing ranks at tiny barriers.

- **The headline opportunity.** Two already-landed, zero-new-code levers were never activated. Activating `--shared-ck-recursion` alone (bit-exact) deletes one ~221 s `ApplySpatialDerivative` subtree; adding `--deriv-cache` (within documented ≤1e-12 tolerance, not bit-exact) replaces the OnTheFly per-QP quadrature in the surviving pass and the volume term with precomputed dense mat-vecs. The genuinely *new* work is the never-implemented "Lever 4" (hoist the two per-step guard `Allreduce`s out of the loop) and fault-aware partitioning.

---

## 2. Measured hotspot table

Caliper region tree, run 51443560. Times in seconds, aggregated over ranks; `%` is of total runtime (avg-rank basis).

| Region (tree) | Min/rank | Max/rank | Avg/rank | % total | Calls/rank |
|---|---:|---:|---:|---:|---:|
| `seas::spatial_dyn::step` | 830.60 | 830.66 | 830.65 | 99.87 | 1445 |
| ` AdvanceADERWithSubStep` | 827.39 | 829.01 | 828.27 | 99.59 | 1445 |
| `  ComputeADERSubStepStates` *(CK #1)* | 236.62 | 257.02 | 246.97 | 29.69 | 1445 |
| `   ApplySpatialDerivative` *(HOT)* | 216.35 | 230.31 | **222.47** | **26.75** | 4335 |
| `  MPI_Waitall (halo)` | 0.42 | 42.10 | 14.76 | 1.77 | — |
| `  friction_substep` | 2.62 | 32.21 | 21.76 | 0.89 | — |
| `  AdvanceADER` *(CK #2 lives here)* | 515.71 | 575.02 | 551.74 | 66.34 | 1445 |
| `   ApplySpatialDerivative` *(SAME HOT kernel, 2nd run)* | 215.45 | 229.56 | **221.43** | **26.62** | 4335 |
| `   ComputeADERVolumeUpdate` | 81.29 | 86.50 | 83.64 | 10.06 | 1445 |
| `    ComputeVolumeRHS` | 81.27 | 86.48 | 83.62 | 10.05 | 1445 |
| `   ComputeADERFaceFluxRHS` | 160.71 | 211.66 | 187.72 | 22.57 | 1445 |
| `    MPI_Allreduce` *(1-int guard)* | 0.16 | 68.16 | 42.43 | 5.10 | 1445 |
| `   ComputeADERSharedFaceFluxRHS` | 22.40 | 47.25 | 37.71 | 4.53 | 1445 |
| `    MPI_Allreduce` *(1-int guard + reconcile)* | 7.09 | 29.37 | 20.24 | 2.43 | 2890 |
| `    MPI_Isend/Irecv` (11.06M total) | — | — | ~0.31 | 0.04 | — |
| `    MPI_Waitall` | 0.47 | 1.15 | 0.84 | 0.10 | — |
| ` MPI_Allreduce` (top-level: ComputeMaxDt etc.) | 0.67 | 2.36 | 1.42 | 0.17 | — |

The `4335 calls/rank = 3 × 1445` for each `ApplySpatialDerivative` subtree = 3 directions × 1445 steps within a single CK recursion; the kernel appears in **two** subtrees because the recursion runs twice.

### Run-totals table (summing the duplicated kernel)

| Aggregate | Time (avg/rank) | % total |
|---|---:|---:|
| `ApplySpatialDerivative` (subtree #1, ComputeADERSubStepStates) | 222.47 s | 26.75 % |
| `ApplySpatialDerivative` (subtree #2, AdvanceADER) | 221.43 s | 26.62 % |
| **`ApplySpatialDerivative` TOTAL (run TWICE/step)** | **443.90 s** | **53.4 %** |
| `ComputeVolumeRHS` | 83.62 s | 10.1 % |
| `ComputeADERFaceFluxRHS` real compute (ex-Allreduce) = 187.72 − 42.43 | 145.29 s | 17.5 % |
| Two 1-int guard `MPI_Allreduce` (42.43 + 20.24) | 62.67 s | 7.5 % |
| Substep halo `MPI_Waitall` | 14.76 s | 1.77 % |
| `friction_substep` | 21.76 s | 0.89 % |

---

## 3. Comparison to the 2026-05-31 baseline (67.6 % → 53.4 %)

The prior baseline placed `ApplySpatialDerivative` at **67.6 %** of runtime; this run measures **53.4 %**. The shift reflects what the landed levers bought *in the codebase* versus what they delivered *in this run*.

**What Levers 1/3 (the `--deriv-cache` geometry caches) and the surrounding restructuring bought:**
- The per-call cost model improved: hoisted scratch allocations (the `KdQ` accumulator at `dynamic/wave_operator.inl:1065` is now function-scope, not per-element), precomputed `elem_mass_inv_[e]`, and the cached-branch infrastructure (`elem_deriv_op_`, `elem_volume_op_`) that *can* collapse the per-QP quadrature into a dense mat-vec. The relative share of `ApplySpatialDerivative` fell from 67.6 % toward 53.4 % as adjacent corrector and MPI costs grew their share and per-call overheads were trimmed.

**What they did NOT do (the un-wired fusion):**
- **Lever 2 (CK fusion) is still inert.** The profile is dispositive: `ApplySpatialDerivative` appears in two subtrees at ~221–222 s each. If `ComputeADERSubStepStatesAndIntegral` were active (one CK recursion feeding both predictors), only one subtree would exist. The driver branch at `drivers/spatial_dyn_driver.cpp:459`/`:513` selects the **separate** `ComputeADERSubStepStates()` + `AdvanceADER(..., nullptr)` path because `use_shared_ck = false`. The fusion is fully implemented (`dynamic/wave_operator.inl:1510`) and proven bit-exact (parity test asserts `max|merged − separate| ≤ 0.0` across orders {2,3,4} and both DerivModes), but the SAFS sbatch never enables it.
- **Lever 1/3 (the geometry cache) is also still inert in this run.** Default `deriv_mode_ = DerivMode::OnTheFly` (`dynamic/wave_operator.hpp:888`); only `--deriv-cache → SetDerivMode(Cached)` (`drivers/spatial_dyn_driver.cpp:1299`–`1302`) flips it, and the sbatch does not pass it. So the 53.4 % `ApplySpatialDerivative` and the 10.1 % `ComputeVolumeRHS` are **OnTheFly per-QP quadrature costs**, not cached mat-vec costs.

**Net:** the baseline-to-now drop is a real measurement of the codebase becoming faster per call and the imbalance/face-flux share rising, **not** evidence that the fusion or the cache ran. Both of the structural levers that would attack the 53.4 % directly remain switched off in production.

---

## 4. Root-cause analysis per region

### 4.1 `ApplySpatialDerivative` — 53.4 %, and why it runs twice

The kernel computes the element-local L2-projected derivative `dQ_c = D_d^e Q_c` per element `e`, per direction `d`, per state component `c` (`dynamic/wave_operator.inl:1004`). Two branches are selected by `deriv_mode_`:

- **Cached branch** (`dynamic/wave_operator.inl:1029`, *verified*): a hand-rolled triple loop `s += D(i,j)*Qc[j]` over a 4×4 `DenseMatrix D = elem_deriv_op_[e][dir]`. No quadrature, no Jacobian. ~252 flop/element/dir. **Not active here.**
- **OnTheFly branch** (`dynamic/wave_operator.inl:1067`, **the branch that ran**): per element it builds `IntRules.Get(geom, 2*order_)` (~4 QPs), heap-allocates `Vector shape(4)` + `DenseMatrix dshape(4,3)`, and per QP calls `Tr->SetIntPoint`, `Tr->Weight()` (3×3 Jacobian determinant), `fe->CalcShape`, and `fe->CalcPhysDShape` (**which inverts the 3×3 Jacobian per QP**). The actual arithmetic is only ~430 flop/element/dir.

**Root cause.** On an **affine P1 tet mesh** the OnTheFly path is pathological: `CalcPhysDShape` re-inverts a Jacobian that is *constant per element* at every quadrature point, `Tr->Weight()` recomputes a constant determinant per QP, and `shape`/`dshape` are heap-allocated per element. Per the findings (`code_refs`: `dynamic/wave_operator.inl:1067,1088,1104,1132`), **>90 % of the cost is fixed per-QP/per-element overhead — Jacobian inverse, MFEM virtual `GetFE`/`GetElementTransformation` dispatch, allocation — not arithmetic.**

**Why twice.** The driver runs the CK Taylor recursion twice over identical `(Q, dt, order)`:
1. `ComputeADERSubStepStates` (`dynamic/wave_operator.inl:1392`, inner derivative loop at `:1465`) builds the substep nodal Taylor states (factorial phase `denom = k+1`, R-003).
2. `AdvanceADER` with `I_precomputed == nullptr` (`dynamic/wave_operator.inl:5721`) falls into `ComputeADERTimeIntegrated` (`:1293`, inner loop at `:1343`) and runs the **same** `D(k+1) = −∑_d A_d ∂_{x_d} D(k)` recursion to build the time integral `I` (factorial phase `denom = k+2`).

Nothing in the numerics forces two passes; `ComputeADERSubStepStatesAndIntegral` (`dynamic/wave_operator.inl:1510`) emits the identical `D(k)` sequence once and accumulates both outputs with their distinct factorial phases. The driver branch that selects fusion (`drivers/spatial_dyn_driver.cpp:459`/`:513`, *verified*) is gated by `use_shared_ck`, default `false`.

Secondary structural note: even with the cache active, the kernel is memory-bound at P1 — the 9 components are looped one 4×4 mat-vec at a time (no BLAS3 reuse of `D`), streaming ~36 Q-doubles + 16 D-doubles for ~252 flop (intensity ~1.2 flop/byte). The adversarial verdicts on `OPT-DERIV-GEMM-BATCH`, `OPT-DERIV-GROUP-GEMM`, `OPT-SOA-LAYOUT`, and `OPT-MIXED-PRECISION` all **reject** further micro-optimization at P1: a flop-neutral repack or a 4×4 dgemm cannot help a bandwidth-bound kernel, and the SoA premise was refuted against the code (`D[e][dir]` is the larger operand, already contiguous; each `Qc` block is L1-resident across its `ndof²` mat-vec).

### 4.2 `ComputeVolumeRHS` — 10.1 %

`ComputeADERVolumeUpdate` (`dynamic/wave_operator.inl:2560`) only zeroes/validates `rhs`, then forwards verbatim to `ComputeVolumeRHS(I, rhs)` (`:2581`), called **exactly once per macro-step** from `AdvanceADER` (`:5736`) — the corrector consumes the time integral `I`, not `Q`. The OnTheFly branch (`dynamic/wave_operator.inl:922`–`991`) again runs the per-QP `CalcShape`/`CalcPhysDShape`/Jacobian-inversion + quadrature sweep, with per-element heap allocs of `shape`/`dshape` (`:944`–`945`). The 9×9 reference-star matrices `Ax/Ay/Az` are *not* recomputed (cached at GodunovFlux construction, `dynamic/godunov_flux.cpp:140`); the 10.1 % is entirely the per-QP transform/quadrature overhead on an affine mesh. The cached branch (`:869`, Lever 3) would replace it with a precomputed `S_d^e` geometry mat-vec plus the per-element star apply.

### 4.3 `ComputeADERFaceFluxRHS` + its per-step guard `Allreduce` — 22.6 %

Called once per macro-step from `AdvanceADER` (`dynamic/wave_operator.inl:5736`), 1445 calls. Two independent root causes:

**(A) The 42.43 s / 5.1 % `MPI_Allreduce` is a per-step barrier guarding a STATIC condition.** *Verified at `dynamic/wave_operator.inl:3848`–`3863`*: an `MPI_Allreduce(MPI_MIN)` of `has_bulk_bg_ ? 1 : 0`, then `MFEM_VERIFY(min_has == 1)` with an R-1505 deadlock-prevention comment. `has_bulk_bg_` is set **exactly once** at setup (`SetAbsorbingBackground`, `drivers/spatial_dyn_driver.cpp:1999`; sole writer at `dynamic/wave_operator.hpp:393`) and never mutated per step. The reduction re-establishes a fact invariant across all 1445 steps. The min/max spread (0.16 s / 68.2 s) proves the reduction itself is trivial — **essentially all 42.43 s is load-imbalance wait** (fast ranks blocking until the slowest fault-flux/friction rank arrives).

**(B) The 145.3 s of real face-flux compute (17.5 %)** is dominated by per-QP recomputation of face-constant quantities on this affine mesh: `CalcOrtho(Face->Jacobian())` and the normal recomputed every QP; per-QP `Vector shape1/shape2` `SetSize`. *Caveat from the verdicts:* the `OPT-FACEFLUX-ROT-CACHE` proposal targeting `GodunovFlux::Interior`'s per-QP rotation rebuilds was **rejected** for this run — because `interior_flux="matrix"` instantiates `BimaterialWaveOperator` (`drivers/spatial_dyn_driver.cpp:1278`), whose non-fault `InteriorFaceFlux_` applies **precomputed per-face matrices** (`dynamic/bimaterial_wave_operator.inl:1044`), so there is no per-QP `BuildFrame`/`BuildRotation` on the non-fault path to hoist. The genuine per-QP rotation rebuilds remain only in the **fault branch** (`dynamic/wave_operator.inl:4109`–`4112`, `T_can`/`Tinv_can` per QP), a heavier but smaller-count path; the fault-rot-reuse proposal was itself **rejected** for a real MPI-determinism hazard (it would key the cache on the rank-local `elem1_on_plus` for both the interior and the shared branch, but the shared branch must negate with the rank-independent `qpd.sign_flipped`).

### 4.4 `ComputeADERSharedFaceFluxRHS` — 4.53 %

The parallel shared-face corrector (`dynamic/wave_operator.inl:4877`), one call/step from `AdvanceADER` (`:5741`). Caliper attributes **2890 `MPI_Allreduce` calls = 2/step** to this region: (i) the mirror `has_bulk_bg_` guard at `dynamic/wave_operator.inl:4920` (*verified*, identical R-1505 pattern) plus (ii) the nested `any_shared` `MPI_Allreduce(MPI_MAX)` inside `ExchangeAndPairSharedFaultQPs` (`:6122`). Both reduce loop-invariant predicates (static topology) — the 20.24 s avg (min 7.09 / max 29.37) is again almost pure imbalance wait. The halo exchange here is the **legacy vdim=1 form**: a loop over `NUM_STATE = 9` components each calling `ExchangeFaceNbrData` (`:4964`–`4975`) — 9 point-to-point rounds/step, even though the sibling substep path was already batched to one `vdim=NUM_STATE` exchange (`ghost_gf_full_state_`, R-1601, `:2242`). The header explicitly flags this as a known, deferred inefficiency (`dynamic/wave_operator.hpp:1134`, R-1504 deep-copy guard). The reconcile (`:6093`) additionally does a **global** `MPI_Allgather` (`:6193`) + `MPI_Allgatherv` (`:6197`) of every rank's shared-fault-QP records.

> Verdict caveats on this region: the shared-face Allreduce-hoist and halo-batch proposals were rated **defer/reject** because (a) the wait relocates to the immediately-following mandatory `Allgather`/`Allgatherv`, recovering only microseconds of genuine 1-int latency; (b) the `MPI_Isend/Irecv` payload is only ~0.31 s, so 9×→1× batching has near-zero direct payoff; and (c) the batched-halo unpack must route through `GetFaceNbrElementVDofs` (R-004) — the exact byNODES layout that previously caused the SAFS normal-velocity runaway. The 2026-05-31 plan itself lists this region under "explicitly NOT worth optimizing."

### 4.5 Load imbalance / partition

Partitioning is plain **unweighted METIS** (`ParMesh pmesh(comm, smesh, part_data)` with `part_data == nullptr` unless `--partition-fault-locality` is passed, `drivers/spatial_dyn_driver.cpp:1007`/`:1030`, *verified flag exists*). METIS balances **element count only** (one unit weight per tet), blind to fault faces and friction cost. `dt` and substep count are globally uniform (CFL uses globally-min `h_min_`, `dynamic/wave_operator.inl:111`; `O = max(1, ader_order)` identical on every rank, `:3002`), so every rank runs the same 1445 steps and the same substeps — **all imbalance is spatial work-per-rank variance, not step-count drift.** See §5.

---

## 5. The load-imbalance section

The compute kernels are balanced; the imbalance is **fault-coupled work surfaced as wait time at per-step collectives.** Evidence:

| Sync point | Min | Max | Avg | Spread | Interpretation |
|---|---:|---:|---:|---:|---|
| `MPI_Waitall (halo)` (substep) | 0.42 | 42.10 | 14.76 | **100×** | imbalance wait on the batched halo exchange |
| `MPI_Allreduce` (face-flux guard) | 0.16 | 68.16 | 42.43 | **426×** | 1-int barrier; ~all of it is wait |
| `MPI_Allreduce` (shared-face guard + reconcile) | 7.09 | 29.37 | 20.24 | 4.1× | 1-int barriers; ~half is wait |

**The 62.6 s of `Allreduce` + 14.8 s of `Waitall` are mostly imbalance wait on tiny reductions, not communication cost.** A 1-int `MPI_MIN`/`MPI_MAX` over 200 ranks is sub-microsecond of actual work. The min-rank values (0.16 s, 7.09 s) are the irreducible collective-latency floor; everything above that is fast ranks blocking until the slowest fault-bearing rank reaches the barrier. The two face guards reduce a **constant boolean** (`has_bulk_bg_`, verified single-writer at setup) — they exist purely as fail-loud R-1505 deadlock prevention and inject a hard global rendezvous twice per step into the hot loop.

**Quantifying the recoverable waste.** Summing the three imbalance waits: `14.76 + 42.43 + 20.24 = 77.43 s` avg. Subtracting the collective-latency floor (the min-rank values `0.42 + 0.16 + 7.09 = 7.67 s`) leaves **~69.7 s of pure imbalance** — ~8.4 % of the 830.65 s wall time that a perfectly balanced partition would recover. The compute regions themselves are the balanced floor: `ApplySpatialDerivative` Max/Avg = 1.035 (3.5 %), because it is a pure per-local-element loop and METIS balances element count.

**Fault-rank concentration hypothesis.** The imbalance is concentrated in fault-bearing ranks: `friction_substep` avg 21.76 / max 32.21 (a ~10.5 s overhang), and the face-flux loop's fault branch (`dynamic/wave_operator.inl:3969`) adds per-QP Riemann + friction Brent solves **only on ranks that own fault faces**. Unweighted METIS cuts the fault plane arbitrarily, so only a subset of ranks carry the extra fault-flux + friction + fault-ghost-packing cost — none of which appears in METIS's unit-per-element weight. This is the performance side of the known SeisSol-parity fault-partition gap (project memory `project_upwind_unstructured_seissol_gap.md`); note that memo concerns the σ_n **numerical** leak (refuted as the accuracy cause), whereas here the *same* missing weighting is the independent, real **load-balance** cause.

> Honest bounding of the partition levers (from the verdicts). `OPT-LB-HOIST-GUARD-ALLREDUCE` was **rejected as a speedup lever**: removing the two guard `Allreduce`s does not delete the per-step rendezvous, because an unconditional per-step `MPI_Allreduce(MAX)` for `V_max` remains at `drivers/spatial_dyn_driver.cpp:3359` plus the `dt` reduction — the wait simply relocates there. So the proposer's "30–50 s reclaimed" is **optimistic**; realistic recovery from the hoist alone is sub-1 s of genuine latency (the hoist remains worthwhile as a one-time-consensus cleanup, just not as a 4–6 % win). `OPT-LB-FAULTLOC-PARTITION` was also **rejected** for this run: the substep `MPI_Waitall` is the wait on the single batched `ExchangeFaceNbrData` over the *entire* shared-face graph (gated on `GetNSharedFaces() > 0`, not on fault faces), so a few-element fault-locality relocation does not shrink per-rank halo volume — and it would shift the σ_n fingerprint on the asymmetric triq mesh, requiring Frontera/Expanse re-validation. `OPT-LB-WEIGHTED-PARTITION` is the *largest* genuine load-balance lever but its claimed 50–70 s was found **~5× optimistic** (the 62.6 s + 14.8 s waits are the *same* stragglers measured at successive barriers, not additive buckets; the recoverable ceiling is the friction overhang, ~5–12 s) and it carries a real bulk-de-balance risk against the already-balanced 53 % `ApplySpatialDerivative`.

---

## 6. What NOT to optimize

- **`friction_substep` (0.89 %, 21.76 s avg).** Do **not** micro-optimize the friction solver itself. Its avg cost is under 1 % of runtime; its *value* in this profile is as an imbalance source (max 32.21 s), which is addressed by partitioning, not by speeding up the Brent solve. The friction substep contains no MPI collective and is closed-form per substep.

- **The OnTheFly per-QP path itself.** Hoisting the affine Jacobian out of `CalcPhysDShape` in OnTheFly (`OPT-VOL-AFFINE-HOIST`) is **subsumed** by `--deriv-cache`, which replaces the entire branch with a precomputed mat-vec. The verdict additionally flagged a correctness hazard: `pmesh.SetCurvature(cfg.mesh.order)` (`drivers/spatial_dyn_driver.cpp:1035`) makes an affine-only hoist silently wrong for config-reachable `order ≥ 2`. And for the TPV*/BP5 byte-exact default path (OnTheFly), re-associating the round-off would break the regression contract. Leave OnTheFly alone; switch the cache on instead.

- **ADER order reduction (`OPT-ADER-ORDER`).** 1.0× — reject. `ader_order = 2` is already the hard floor; all four CK entry points assert `order ∈ {2,3,4}`. Lowering it is not a valid ADER scheme and changes the temporal truncation error and CFL limit, breaking SeisSol parity.

- **Mixed/single precision (`OPT-MIXED-PRECISION`).** Reject. `real_t` is a global MFEM compile-time typedef; fp32 in the CK recursion injects ~1e-7 relative error against the documented ≤1e-12 Lever-1 gate, and the SAFS upwind-ADER fingerprint is documented-sensitive at the 1e-7..1e-12 level (σ_n floor ~1.19e-7 Pa, [[v_n]] seed regime). Far outside any documented tolerance; would change the SCEC-parity fingerprint.

- **GEMM/SoA micro-optimizations of the cached matvec at P1.** Reject/defer. At `ndof = 4` the kernel is memory-bound, the flop count is unchanged by batching, and the mandatory strided gather/scatter (Q is component-major, `c*ndof_total_ + dof_offset`) plus netlib small-DGEMM overhead make these net-neutral or negative. Revisit only at `fe_order ≥ 3` **after** the levers are active and a roofline confirms compute-bound.

- **Shared-face halo batching / reconcile rework / overlap.** Defer. The direct payload (`Isend/Irecv` ~0.31 s, local `Waitall` 0.84 s) is sub-0.1 % of runtime; the dominant cost in that region is the imbalance-wait `Allreduce`, which batching does not touch. The 2026-05-31 plan already designated `ComputeADERSharedFaceFluxRHS` as not worth optimizing.

**The only first-order actions** are therefore the two flag activations: `--shared-ck-recursion` (bit-exact, deletes one ~221 s `ApplySpatialDerivative` subtree) and `--deriv-cache` (within ≤1e-12 tolerance, collapses the surviving OnTheFly quadrature in both `ApplySpatialDerivative` and `ComputeVolumeRHS`). Realistic combined wall-clock from ~830 s toward **~400–480 s (~1.7–2.0×)** — *not* the proposer's optimistic ~350–420 s / ~2.4×, because the residual is increasingly MPI-imbalance-bound (the 62.6 s guard `Allreduce`s and 14.8 s halo `Waitall` are untouched and set the max-rank critical path). One validation caveat from the verdict on `--deriv-cache`: the ≤4.4e-16 / ≤1e-12 equivalence tests cover only the **scalar** `WaveOperator`; the SAFS run uses `interior_flux="matrix"` → `BimaterialWaveOperator` with per-element `A_d^e`, so a bimaterial cached-vs-OnTheFly check and a SAFS fingerprint confirmation should precede adoption of `--deriv-cache`. `--shared-ck-recursion` carries no such caveat (bit-exact, `max|merged − separate| ≤ 0.0`).
