# Review record: PLAN_ader_kernel_efficiency (3-lens adversarial, 2026-07-21)

Reviewed artifact: `PLAN_ader_kernel_efficiency_2026-07-21.md` **rev 1**.
Verdicts: TECHNICAL = PASS WITH FIXES · FEASIBILITY = **FAIL** · COMPLETENESS = PASS WITH FIXES.
All findings below were applied in **rev 2** of the plan (same file). Fix column = what rev 2 does.

## Technical lens

| ID | Sev | Finding (compressed) | Applied fix in rev 2 |
|---|---|---|---|
| R-201 | CRIT | Face-table memory arithmetic 15× too small (6.4 KB/face ⇒ 31.5 GB @np=1, 123 MB/rank @np=256, not 2.1 GB/8 MB) | Corrected numbers; redesigned tables: reference-catalog dedup for trace bases (dozens of matrices total — `CalcShape(Loc.Transform(ip))` depends only on the reference embedding), per-face residual 1.3–2.6 KB; local production-mesh runs declared infeasible; R-004 guard re-derived |
| R-101 | MOD | "Bit-exact" claim contradicted the contraction rewrite (re-associates the QP accumulation) | Split into `--face-tables` (bit-exact: legacy loop shape, lookup-only, CachedInterior_ precedent) and `--face-tables-fold` (contraction + T·A±·T⁻¹ fold, ≤1e-12) |
| R-102 | MOD | "Extend FaultBasisQPData" violates no-touch `fault/` | Parallel per-(fb_idx,qp) rotation table owned by `dynamic/face_kernel_tables.hpp`, indexed like `FaultBasisData::qp_data`; `fault/` untouched |
| R-103 | MOD | D(k) write-back named the wrong buffer (ping-pong scratch) and wrong element set | Corrected: scatter to `dk_retain[(slot·order+k)·block+j]` for `retain_slot_of_elem[e] ≥ 0` only, bit-identical to the scalar retain lambda (`wave_operator.inl:1990-2009`); ping-pong needs no write-back |
| R-301 | MOD | Fault share (1.3 %) and B0 window are pre-nucleation-biased; SeisSol's 14.4 includes rupture DR | Phase 0 + Phase 4 measure a rupture window `[1.0,1.5] s`; fault scope-out downgraded to a deferred decision; Summary annotates window provenance |
| R-104 | LOW | R-1203 mis-attributed (it is the mixed-flux guard; face-cache exclusivity is REVIEW R-001) | Corrected throughout + both defined in the Glossary |
| R-105 | LOW | LU tax 20 % (Summary) vs 15 % (Overview) | Reconciled to 15 % with provenance (local np=4 profile, levers ON, no --face-cache) |

Checked-OK highlights: Vandermonde-LU mechanism confirmed in MFEM source (`fe_l2.cpp` Ti/DenseMatrixInverse
→ `LUFactors::Solve`); BatchedLinAlg NATIVE determinism argument confirmed (fixed-order `kernels::AddMult`
per slice); transposed cached-matvec traversal confirmed (`wave_operator.inl:1244-1249`); face-cache
composability with the deriv levers confirmed; trace-table step-invariance premise confirmed.

## Feasibility lens (the FAIL)

| ID | Sev | Finding | Applied fix in rev 2 |
|---|---|---|---|
| R-101F | CRIT | Phase gates don't compose (P1×P2 sub-gates ≈ 2.7× < the demanded cumulative 3×; Phase 3 had no quantified mechanism for 6×) | Cumulative P2 gate reset to ≥2.5×; program target restated as ABSOLUTE ≤~40 µs·core with a roofline-conditioned fallback (~3–3.5×) declared in advance; Phase-3 criterion states the memo supplies the number |
| R-102F | CRIT | 6× used two denominators (245.9 vs post-activation B0) | Single denominator pinned: ≤40 µs·core absolute on the B0 protocol; every gate restated as memo-derived multiples of B0 |
| R-201F | CRIT | Phase-0 "1.5–2.5×" unjustified: face-cache Amdahl caps ~1.3×; the profiled LU cost was measured WITH Apple LAPACK (BLAS presence ≠ removing the tax); 58.4-vs-245.9 is largely hardware | Phase-0 expectation reset to 1.2–1.4×; BLAS reclassified measure-only (no promised gain, moot after Phase 1); build-handicap implication deleted |
| R-202F | CRIT | Phases 0 and 1 double-sold the same LU/shape savings (B0 includes --face-cache; the old 1.5×-vs-B0 gate would fail by construction) | Phase-1 Expanse gate is now SHARE-based (face subphase ≤0.5× its B0 share, plus ≥1.2× per-update) |
| R-203F | CRIT | Phase-2's 2.5× Expanse gate sat above the Rome bandwidth roofline; whole-cluster batches re-stream operators per level = the exact OPT-DERIV-GEMM-BATCH rejection arithmetic, and its roofline precondition was never scheduled | Batch granularity redesigned to cache-blocked tiles (E≈32–128, ≤L2 panels, levels fused inside the tile); HW-counter roofline added to Phase 0 as a PRECONDITION; Expanse expectation reset to 1.5–2× (memo-derived); 2.5× kept as the LOCAL compute-bound gate only |
| R-103F | MOD | (= tech R-201 memory arithmetic) | Same fix |
| R-204F | MOD | (= tech R-301 rupture-window bias, Phase-4 certification angle) | Certified number states its window; both windows measured in Phase 4 |
| R-301F | LOW | All numeric gates pre-committed from an Apple-silicon profile before B0 exists | All downstream gates marked provisional; gate-reset memo is a mandatory Phase-0 exit artifact |

Checked-OK highlights: baseline numbers real and correctly quoted; ≤40 µs ⇒ within 2.78× of SeisSol
arithmetic consistent; BatchedLinAlg not vaporware in 4.9.1; prior-art citations honest (rejections + revisit
clauses correctly scoped); comm-work scoping correct.

## Completeness lens

| ID | Sev | Finding | Applied fix in rev 2 |
|---|---|---|---|
| R-301C | CRIT | No SAFS applicability statement — as written NO phase applies to production SAFS (bimaterial/mixed-flux parse-rejected) while the Summary implied parity | New **Applicability** section with a per-phase matrix; SAFS extension named as Phase 5 with trigger; Summary states v1 is TPV104-class only |
| R-302C | CRIT | No Interfaces sections; Phase-2 kernels had no signatures | Interfaces added to Phases 0–2 (FaceKernelTables struct + builder incl. the IntegrationRule provenance; `BatchedCKLevelTile` + `ADERBatchScratch` + dispatch hook) |
| R-101C | MOD | (= memory arithmetic) | Same fix |
| R-102C | MOD | Phase-0 "no numerics change" vs LAPACK relink (round-off-level); build_expanse.sh outside stated editable surface; B0-before-or-after-rebuild unpinned | Relink reclassified round-off with physics sanity acceptance vs the 52344266 reference; B0 pinned post-rebuild with binary hash; `build_expanse.sh` added to the editable surface |
| R-201C | MOD | The "separate, already-scoped comm plan" doesn't exist as a document; no sequencing/conflict policy for the shared files | "Relationship to the comm plan" section: doc not yet written (stated), kernel 0–2 land first, compute-only B0 metric isolates, B0 re-measured on interleaved landing |
| R-304C | MOD | Glossary gaps (R-1203, B0, OPT-* pointer) | All added; action_plan doc added to header references |
| R-307C | MOD | No tests for the flag matrix; no curved-mesh fixture | Both added to Testing Strategy + Phase-1 requirements |
| R-305C | MOD | No per-phase review/fix loop | Phase close-out protocol added (REVIEW_<phase>_impl docs; no CRITICAL = close) |
| R-306C | MOD | No rollback stories | Rollback line added to every phase |
| R-303C | MOD | Phase 3 had no file surface/edge cases; Phase 4 had no Goal/Files/Deps | Both restructured (Phase-3 surface marked provisional; Phase-4 partial-adoption edge case) |
| R-103C | LOW | Phase-2 local-vs-Expanse gate ambiguity | Split: local ≥2.5× vs LOCAL cached path; Expanse gates B0-denominated only |
| R-104C | LOW | BatchedLinAlg API glossed (Vector panels, not DenseTensor; global backend singleton) | Call written as `Get(NATIVE).AddMult(..., Xb_vec …)` with setup `MFEM_VERIFY` pin |
| R-308C | LOW | Doc/memory/flag-disposition unassigned | Phase-4 requirements extended (CLAUDE.md, --help, memory note, --face-cache disposition) |
| R-309C | LOW | Summary jargon | De-jargoned; class/flag names moved to the Technical Overview |

## Net verdict after fixes

Rev 2 incorporates every finding. The plan's claims are now either verified against source
(file:line) or explicitly marked provisional pending the Phase-0 measurement, and the
performance model is Amdahl-consistent with declared fallbacks. Status: **READY — Phase 0 is
the next actionable step (needs user approval for the two Expanse baseline submissions).**

---

# Round 2 (2026-07-22): evidence-backed re-review of rev 2 → produced rev 3

Fresh 3-pass review + a 3-agent EVIDENCE round (microbenchmark, roofline, comm pre-scope).
Verdict: **PASS WITH FIXES** — "the plan does NOT yet support the claim that this move brings
the code close to SeisSol runtime" until the two CRITICALs land. All findings applied in rev 3.

| ID | Sev | Finding (compressed) | Disposition in rev 3 |
|---|---|---|---|
| R-401 | CRIT | End-to-end composition unstated: kernel success alone = ~2677 s/sim-s = 11.2× behind SeisSol; only kernel+comm reaches 2.4–3.9× | Composition subsection added to the Summary; scenario table in the EVIDENCE doc; Phase 4 publishes the measured version |
| R-402 | CRIT | ≤40 µs rested on an unmeasured assumption (NATIVE batched GEMM rate on Rome); no prototype gate before implementation | **Spike EXECUTED** (local half): `tests/bench/bench_ader_kernel_variants.cpp` — A 5.2 µs (48 % peak), B 1.25×, **C BatchedLinAlg 0.46× — mechanism eliminated**, D 1.30×; Rome contended rerun added to Phase 0 as the Phase-2 go/no-go (≥2× contended B-vs-A) |
| R-403 | MAJOR | Kernel-first sequencing asserted, not argued; comm plan unwritten | Comm plan WRITTEN (`document/comm_dev/PLAN_lts_comm_reduction_2026-07-22.md`, from the full exchange inventory); sequencing re-argued: shared Phase 0, then parallel comm-merge + face-tables, with a revisit trigger |
| R-404 | MAJOR | SeisSol frozen at its favorable-to-MFEM 239; binary already drifted once | Phase 4 pins the 52365078 reference AND adds one async-IO SeisSol re-run to bracket 125–239 |
| R-405 | MAJOR | Phase-1 share gate not protocol-pinned (parallel Phase-2/comm confound) | Gate pinned: B0 protocol + Phase-1 flags only, pre-rupture window |
| R-406 | MAJOR | REAL SPEC BUG: builder cited `2*order_` rule for all faces; fault faces use `FaultFaceQuadDegree()` (wave_operator.inl:5524-26) | Interfaces comment corrected; fault tables sized by the fault rule's nqp |
| R-407 | MOD | Catalog coverage holes: nonconforming, boundary faces | `MFEM_VERIFY(mesh.Conforming())` + parse-time reject; boundary-face disposition stated |
| R-408 | MOD | No ragged-tile policy (finest cluster ~3.5 elems/rank at np=256) | `E_min` fallback to legacy path per (rank, cluster); correctness-neutral |
| R-409 | MOD | 245.9-vs-14.4 confounds kernels with node layout (128 pure-MPI vs 8×15) | 64-vs-128 ranks/node A/B leg added to the B0 job |
| R-410 | MOD | Certified flag set/leg ambiguous | ONE pinned flag line + pinned leg, decided at Phase-1 close-out |

**Evidence-round numbers now in the record** (EVIDENCE_microbench_roofline_2026-07-22.md):
fused-design bandwidth floor 14–18 µs (level fusion LOAD-BEARING; non-fused floor 45–49 µs
kills the target); 40 µs needs 3.75–4.0 GFLOP/s/core (stretch), fallback 3–3.5× needs 1.8–2.3
(near-certain); comm both-levers → exposed Waitall 2341 → ~150–400 s/sim-s (~1.8–2.0× LTS leg)
with the skew caveat gated by comm-plan Phase 0.
