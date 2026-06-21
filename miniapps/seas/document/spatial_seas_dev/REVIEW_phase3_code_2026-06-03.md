# Code Review: Phase 3 code changes — spatial_seas QD driver (2026-06-03)

## Review Scope
- Plan: `document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md` (Phase 3 — fault geometry + spatial per-DOF friction).
- Files reviewed:
  - `domain/elasticity_operator.hpp` + `domain/elasticity_operator_traction.inl` — 3 owned-order getters.
  - `fault/fault_geometry.hpp` — `(domain, seed, mpi, bool)` ctor + template `SetRateStatePerDOF`.
  - `drivers/spatial_seas_driver.cpp` — Phase 3 normative-order wiring (section 4b).
  - `tests/unit/test_spatial_seas_faultgeom_parity.cpp` — new parity test.
  - `Makefile` — parity-test target/obj/vars.
- Domain context: `CLAUDE.md` (no silent fallback), the plan (R-001/R-002/R-003/R-007/R-010), `seas_config_bridge.hpp`/`seas_driver.cpp` (resolver-input precedent), `bp5_params.hpp`, the implementer report, prior reviews (not clobbered).
- Method: three adversarial passes; build + run exercised (parity np=1/np=4; driver np=1/np=4; extreme-care no-regression via stash).

## Build / run evidence
- Parity test: np=1 **PASS** (0/279372), np=4 **PASS** (0/279408, global owned 27936 consistent) — proves getter owned-ordering ≡ geom ordering (R-002) and the `SetRateStatePerDOF` mapping incl. the R-010 V_init decomposition.
- Driver Phase 3 wiring: np=1/np=4 resolve+apply, `geom.HasParams()==false` (R-001).
- Extreme-care no-regression: `seas_driver` + `spatial_dyn_driver` rebuild clean; the `test_bp5_fault_operator` abort is **confirmed pre-existing** (identical with my changes stashed; `setup.inl` untouched) — out of scope.

## Findings

### [R-301] MODERATE [POSSIBLE] [spatial_seas_driver.cpp:4b] — the b/f0/V0 uniformity guard is per-rank LOCAL, not MPI-global

**Category:** BUG (MPI correctness) / ASSUMPTION

**Description:**
The Phase 3 wiring guards against per-DOF `b`/`f0`/`V0` (which the scalar `DieterichRuinaFriction`/`AgingLawPsi` cannot represent) with `assert_uniform(rs.b/f_0/V_0)`. But `assert_uniform` computes `Min/Max` over the **local** owned DOFs only — there is no `MPI_Allreduce`. The resolver DOES allow per-DOF `b` (Phase 11a: a spatial rule can set `r.b`). If `b` is globally non-uniform but happens to be **locally uniform on each rank** (e.g. a region-partitioned `b` whose region boundary aligns with the mesh partition), every rank's local check passes, the guard is bypassed, and Phase 5 then silently bakes the scalar `b_default` — exactly the "silently ignore per-DOF b" outcome the guard exists to prevent. (For `f0`/`V0` the resolver already rejects per-rule overrides — `spatial_friction.cpp:2117` — so those two checks are redundant; only `b` can actually vary.) Harmless for BP5 (globally-uniform `b`), latent for partitioned SAF configs.

**Trigger:** `mpirun -np >1` with a `[friction.rate_state]` spatial rule that makes `b` vary across ranks but uniform within each rank.

**Actual behavior:** local checks pass on all ranks; non-uniform `b` is silently accepted (then ignored in Phase 5).

**Expected behavior:** detect global non-uniformity and abort.

**Suggested fix:** reduce local min/max to global before comparing:
```diff
       auto assert_uniform = [&](const Vector &v, const char *name)
       {
-         if (v.Size() <= 1) { return; }
-         const real_t lo = v.Min(), hi = v.Max();
-         const real_t mean = v.Sum() / v.Size();
+         real_t lo =  std::numeric_limits<real_t>::infinity();
+         real_t hi = -std::numeric_limits<real_t>::infinity();
+         for (int i = 0; i < v.Size(); ++i) { lo = std::min(lo, v(i)); hi = std::max(hi, v(i)); }
+         MPI_Allreduce(MPI_IN_PLACE, &lo, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
+         MPI_Allreduce(MPI_IN_PLACE, &hi, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
+         if (!std::isfinite(lo) || !std::isfinite(hi)) { return; }  // no DOFs anywhere
+         const real_t mean = 0.5 * (lo + hi);
          MFEM_VERIFY(hi - lo <= 1e-12 * (std::abs(mean) + 1.0),
                      "spatial_seas: per-DOF " << name << " is not uniform ...");
       };
```
(All ranks must call the Allreduce — keep `assert_uniform` outside any rank-gated block, which it already is.)

**Test case:**
```cpp
// test_R301_global_b_nonuniform_aborts:
//   np=2; spatial rule sets b on a region that lands entirely on rank 0.
//   EXPECT a collective abort (not a silent pass).
```

---

### [R-302] LOW [POSSIBLE] [fault_geometry.hpp:SetRateStatePerDOF] — V_init dip set to exactly 0, vs BP5's V_zero=1e-20 floor

**Category:** ASSUMPTION (forward-looking, Phase 5)

**Description:**
`SetRateStatePerDOF` with `init_vel_dir=(0,1)` sets `V_init_vec_(2i)=0` (dip). `ComputeBP5Params` sets it to `bp5_params::V_zero = 1e-20` (a deliberate non-zero floor, `bp5_params.hpp:106`). The difference is 1e-20 (≪ the 1e-10 parity tol, so PART A passes), but if the Phase 5 QD fault operator relies on a non-zero per-component velocity floor (e.g. a `log(V)` / `V/V0` term evaluated per component rather than on the magnitude), a hard 0 dip could misbehave. Not a Phase 3 correctness issue; flagged so Phase 5 confirms the rate-state ODE uses the slip-rate magnitude (where 0 vs 1e-20 is immaterial) and not a per-component log.

**Suggested fix:** none for Phase 3. Phase 5: verify the fault-op consumes `|V|`; if a dip floor is needed, set `init_vel_dir` from the magnitude+floor rather than a pure unit strike.

---

### [R-303] LOW [test_spatial_seas_faultgeom_parity.cpp] — GetFaultDOFIntegrationPoints values are unverified (size-only)

**Category:** QUALITY (coverage)

**Description:**
The test calls `GetFaultDOFIntegrationPoints` and asserts only `dof_ips.size() == n_owned`. The actual IP values (the reference coordinates in Elem1's frame) are never checked, so a bug in the IP walk (wrong rule, wrong `GetElement1IntPoint`, mis-restriction) would not be caught until Phase 5 evaluates a Coefficient there. The other two getters are covered indirectly (elem via the resolver run, coords via the BP5 parity).

**Suggested fix:** add a value check — e.g. for each owned DOF, transform its IP through Elem1 and assert the resulting physical coordinate matches `dof_coords_3d(3i..3i+2)` to 1e-10:
```cpp
// for owned dof i: ElementTransformation *T = pmesh.GetElementTransformation(dof_to_elem[i]);
//   T->Transform(dof_ips[i], phys);  assert |phys - dof_coords_3d(3i..)| <= 1e-10
```

**Test case:** the snippet above (folded into the existing test).

---

### [R-304] LOW [fault_geometry.hpp:SetRateStatePerDOF] — second call updates a/dc/eta/V_init but NOT sigma_n (stale-σ_n inconsistency)

**Category:** EDGE_CASE

**Description:**
The `sigma_n_per_dof_` write is guarded by `if (sigma_n_per_dof_.Size() != N)`. On a *second* `SetRateStatePerDOF` call (same geom), `a/dc/eta/V_init` are overwritten but `sigma_n_per_dof_` (already size N) is left at the first call's value — a silent inconsistency. The production path calls it once, so this does not bite today, but the "if not already set" semantics conflate "a stress source set it" with "a prior SetRateStatePerDOF set it".

**Suggested fix:** track provenance explicitly, or document that `SetRateStatePerDOF` is single-shot:
```diff
+      // NOTE: single-shot — call once before ComputeParams*.  The σ_n "if not
+      // already set" guard assumes any pre-existing sigma_n_per_dof_ came from
+      // a stress source, not a prior SetRateStatePerDOF call.
       if (sigma_n_per_dof_.Size() != N)
```

---

### [R-305] LOW [fault_geometry.hpp:SetRateStatePerDOF] — unconstrained template parameter

**Category:** QUALITY

**Description:**
`template <typename RSParams> void SetRateStatePerDOF(const RSParams&, ...)` accepts any type; a wrong argument fails with a deep member-access error at instantiation rather than a clear message. Deliberate (keeps `spatial_friction.hpp` out of the extreme-care header), but a sanity check would improve diagnostics.

**Suggested fix (optional):** a `static_assert` on the required members, e.g.
```cpp
static_assert(std::is_same_v<decltype(rs.a), const mfem::Vector&>
              || std::is_same_v<decltype(rs.a), mfem::Vector>,
              "SetRateStatePerDOF: RSParams must expose mfem::Vector a/Dc/eta/V_init/sigma_n_eff");
```
(or just document the contract — current comment already names the type).

---

## Summary
- Critical issues: **0**
- Moderate issues: 1 (R-301 — local-only uniformity guard, POSSIBLE; bites only under MPI with region-partitioned `b`)
- Low issues: 4 (R-302 V_init dip floor [Phase-5 forward]; R-303 IP values untested; R-304 single-shot σ_n; R-305 unconstrained template)
- Plan compliance: **FULL** for Phase 3 —
  - 3 owned-order getters ✔ (R-002 verified by the np=4 parity); new ctor (compute_bp5_params=false) ✔; `SetRateStatePerDOF` ✔ (a/dc/eta/V_init + R-010 decomposition + σ_n-if-unset + does NOT set `params_computed_`, R-001) ✔; normative-order driver wiring ✔; parity test ✔ (np=1 + np=4).
  - Documented deviations justified: `SetRateStatePerDOF` template (avoids heavy include in the extreme-care header); parity uses direct-bp5 `rs` for the 1e-10 a/dc/eta/V_init parity (the resolver cannot reproduce BP5's 2-D `a(x2,x3)`) while PART B exercises the real `ResolveRateState` — together they cover the plan's stated purpose ("proves the spatial→QD wiring + DOF ordering").
  - Extreme-care invariant held: changes to `elasticity_operator.hpp`/`fault_geometry.hpp` are strictly additive; `seas_driver`/`spatial_dyn_driver` build; BP5 fault-op behavior unchanged (the one abort is pre-existing, verified).
  - R-003 verified: TOTAL σ_n (empty ⇒ `sigma_n_default`) + zero `PorePressureSpec` avoids double-P_p.
- Verdict: **PASS WITH FIXES** — no CRITICAL/MODERATE bug fires on the current build/configs; fix R-301 to make the uniformity guard MPI-correct before any multi-rank SAF config with region-varying `b`. The LOW items are hardening + Phase-5 forward notes.

## Unreviewed Areas
- Phase 5 consumption (RateStateFaultOperator ctor caching the resolver-filled arrays; ComputeParams* overwriting σ_n) — out of scope.
- BP5 physics parity of the assembled traction (Phase 7) — out of scope.
- `seas_test_bp5_fault_operator` (`ValidateFacetBCTables` abort) — pre-existing, verified via stash, explicitly out of scope.
- Full `make all && make test` not run (heavy); non-interference verified by rebuilding both extreme-care consumers + the config test.
