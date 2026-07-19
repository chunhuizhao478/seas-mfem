# Code Review: LTS Phase 3 primitives (P3-1..P3-3) — 2026-07-19

## Review Scope
- Plan: `document/lts_dev/PLAN_clustered_lts_ader_2026-07-18.md` (A.7, D-3, P-006, Phase 3 Step 0, NORMATIVE §§, Appendix B).
- Files reviewed: `dynamic/friction_substep_iterator.{hpp,cpp}`, `dynamic/friction_iterator.hpp`, `dynamic/spatial_nucleation.{hpp,cpp}`, `dynamic/nucleation_method.hpp`, `dynamic/wave_operator.{hpp,inl}`, `dynamic/bimaterial_wave_operator.{hpp,inl}`, `io/tpv104_checkpoint.hpp`, `tests/unit/test_lts_{friction_range,nucleation_absolute,fault_reorder,predictor}.cpp`.
- Domain context: CLAUDE.md (sign conventions, byte-exact contract); the no-local-reproducer constraint.
- Method: three parallel adversarial reviewers (friction/nucleation range; fault-face reorder; checkpoint canonical order), findings consolidated + de-duplicated below.

**Consensus: no active numerical bug in the reviewed diff.** Verified correct by tracing: memset byte-exactness for `(0,n)`, the range loop-index change, `I_imp` out-of-range untouched, the reorder permutation semantics + placement + no-op path, the D-1 assert, bimaterial ctor forwarding, and the checkpoint write/read being genuine inverses with a layout-independent on-disk order. Findings are integration-boundary, robustness, and test-coverage items.

## Findings

### [R-001] MODERATE [test_lts_friction_range.cpp:main] — B.9 never exercises the SRW global-index path
**Category:** EDGE_CASE (test gap)
**Description:** B.9 builds only a `RateStateAgingIterator`. The one range-body operation that depends on `i` being the GLOBAL fault-QP index (not `iu - qp_begin`) is `StepOneQP_`'s `StatePolicy::UpdatePsi(law_, d, s.V_abs, dt_sub, extra_, i)` — for SRW this is `(*Vw)(i)`. `RateStateAgingPolicy::UpdatePsi` ignores `i`, so B.9 is blind to a local-vs-global mis-index in the range refactor — the highest-risk line.
**Trigger:** if the range lambda passed `iu - qp_begin` instead of global `i`, the production SRW iterator (TPV104/SAFS) would read `V_w[0]` for the first QP of range B instead of `V_w[k]`; B.9 stays green.
**Expected:** an SRW sub-case with per-QP-varying `V_w` so a local index diverges from a standalone reference fed the matching global `V_w`.
**Suggested fix:** add an SRW block mirroring the aging block, `Vw(i)=0.1+0.02*i`, comparing ranged QP `i` to a single-QP `RateStateSlipLawSrwIterator` reference given `V_w = {Vw(i)}`.

### [R-002] MODERATE [nucleation_method.hpp:INucleationMethod] — no range `ApplyAbsolute`; D-3 "route via INucleationMethod" half-built
**Category:** DEVIATION (incomplete contract; latent)
**Description:** D-3 requires per-cluster nucleation to route through the idempotent absolute forms **via `INucleationMethod`**. The free functions got `(qp_begin,qp_end)` range overloads, but `INucleationMethod::ApplyAbsolute` did not — so the range friction `Advance`'s promised "range-apply callback via the absolute range overloads" cannot be built through the strategy object without bypassing it. Latent (no production caller until the Phase-4 driver wiring).
**Suggested fix:** add `virtual void ApplyAbsolute(dof, t, qp_begin, qp_end)` with a full-range-delegating default (`MFEM_VERIFY(full range); ApplyAbsolute(dof,t)`); override in the two gradual kinds (range appliers), `StaticOverstress` (no-op), `InstantaneousOverstressCircular` (range set).

### [R-003] MEDIUM [wave_operator.hpp / io/tpv104_checkpoint.hpp] — per-FACE perm cannot feed the per-QP checkpoint (expansion trap)
**Category:** ASSUMPTION (unit mismatch; latent)
**Description:** `GetFaultFaceCanonicalPerm()` is per-FACE (length `nfi`); the checkpoint's `dof_canonical_perm` must be per-QP (length `nd = nfi*nbf`). No per-QP accessor or expansion helper exists. Phase-4 wiring would trip the size `MFEM_VERIFY` (nbf>1) or silently scatter wrong (nbf==1).
**Suggested fix:** add `WaveOperator::GetFaultQpCanonicalPerm()` returning the QP-expanded perm (`perm[i*nbf+q] = face_perm[i]*nbf + q`), empty when not reordered; document the interior-QP-first + no-shared-QP (D-2) assumption.

### [R-004] MEDIUM [io/tpv104_checkpoint.hpp:ReadTpv104CheckpointV2Impl] — read-path perm guard is not the true mirror of the write guard (silent corruption on a duplicate)
**Category:** BUG (robustness asymmetry)
**Description:** Write validates a full bijection (`inv[c]==-1` catches duplicates + range). Read checks only per-element range (`c>=0 && c<nd`) — a non-bijective perm passes and silently scatters wrong; no abort, no OOB. Header calls the read a "mirror of the write path"; it is not.
**Suggested fix:** add a `seen[]` bijection check on read before the scatter loop.

### [R-005] LOW [wave_operator.inl:reorder block] — D-2 / cluster-size asserts sit inside the `nfi>0` guard
**Category:** BUG (np>1 tripwire defeated)
**Description:** The D-2 assert (`fault_shared_faces_.Size()==0`) and the cluster-id-size assert are inside `if (lts_cluster_id != nullptr && fault_interior_faces_.Size() > 0)`. A rank with `nfi==0` but ≥1 shared fault face skips the block, so the D-2 tripwire never runs — silent pass. Only the sort needs `nfi>0`.
**Suggested fix:** hoist the size + D-2 asserts to `if (lts_cluster_id != nullptr) { …asserts…; if (nfi>0){ …sort… } }`.

### [R-006] LOW [spatial config parse] — `--fault-resample` + `lts` aborts mid-Advance instead of at parse time
**Category:** QUALITY (error surfacing)
**Description:** The range `Advance` guards resample off via a mid-sweep `MFEM_VERIFY`. A user enabling both gets a deep abort rather than a clean setup-time rejection.
**Suggested fix:** add `lts != off ⇒ reject fault-resample` to the existing `lts + {rk, mixed_flux}` config guards.

### [R-007] LOW [test_lts_predictor.cpp:checkpoint_v2_canonical_reorder] — B.8b asserts only 6 of 9 dynamic fields
**Category:** EDGE_CASE (test coverage)
**Description:** Compares psi/slip1/V2/sigma_n_nuc/slip2/tau2_nuc but never slip_rate/V1/tau1_nuc.
**Suggested fix:** compare all 9 dynamic fields in a loop.

### [R-008] LOW [test_lts_friction_range.cpp:224-225] — B.9 uses `1e-15` absolute tol while claiming bit-exact
**Category:** QUALITY (test clarity)
**Suggested fix:** change the two range-vs-reference CHECKs to `== 0.0` (as B.10 does).

### [R-009] LOW [wave_operator.inl:376] — sort tiebreak is the LOCAL mesh face id; plan says "global face id"
**Category:** QUALITY (Phase-4 note)
**Description:** Benign at np=1 (local==global, unique — verified). Add a comment noting the np=1 equivalence and Phase-4 revisit.

### [R-010] LOW [friction_substep_iterator.cpp:LSW range Advance] — LSW diag maxima reset granularity under LTS (Phase-4 note)
**Category:** QUALITY (forward-looking)
**Description:** The range advance accumulates `slip_rate_substep_max`/`sigma_n_substep_min` correctly, but a whole-vector reset per cluster (Phase-4 driver) would re-zero earlier clusters. Document that the Phase-4 driver resets at range granularity.

## Summary
- Critical: 0 · Moderate/Medium: 4 (R-001..R-004) · Low: 6 (R-005..R-010)
- Plan compliance: FULL for the four locally-scoped Phase-3 primitives; driver assembly (reorder activation + fault interleave) documented Frontera-staged.
- Verdict: PASS WITH FIXES — apply R-001..R-005, R-007, R-008 (correctness/coverage), R-006 (UX), R-009/R-010 (doc notes).

## Unreviewed Areas
- Driver fault-half interleave + reorder-activation wiring — not implemented this session (Frontera-staged); physics acceptance is production-mesh-only.
- np>1 reorder behavior (R-005) — not locally constructible at np=1.
