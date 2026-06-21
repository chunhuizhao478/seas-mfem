# Code Review: Phase 2 code changes — spatial_seas QD driver (2026-06-03)

## Review Scope
- Plan: `document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md` (Phase 2 — mesh, boundary, material, domain operator, homogeneous-first).
- Files reviewed (Phase 2 additions):
  - `drivers/spatial_seas_driver.cpp` — section 4 (mesh / `BoundaryConfig` / homogeneous `LinearElastic` / `ElasticityDomainOperator` construction) + restructured dry-run/skeleton exit.
  - `spatial/code/spatial_friction.hpp` — `BoundarySpec.dirichlet_attrs` (new field).
  - `spatial/code/spatial_friction.cpp` — `[boundary].dirichlet_attrs` parse + disjointness guards.
  - `config/safs_qd/bp5_phase2_smoke.toml` — new construction-smoke config.
- Domain context: `CLAUDE.md` (no silent fallback), the plan, the implementer completion report, `domain/elasticity_operator.hpp` + `..._setup.inl`, `seas_driver.cpp` / `seas_config_bridge.hpp` (the reference construction pattern), prior reviews `REVIEW_phase0*`/`REVIEW_phase01*` (not clobbered).
- Method: three adversarial passes; build + run exercised (np=1 and np=4 on the BP5 1000 m mesh; heterogeneous-abort path; full no-regression set).

## Build / run evidence
- QD driver builds clean; link set unchanged (`spatial_seas_driver.o + spatial_friction.o + MFEM_LIBS`) — confirms `ElasticityDomainOperator` is header-only (no new link dep).
- np=1 dry-run on `bp5_tandem_exact.msh`: operator built, getters `27936 / 27936 / 9312 / 3`, **exit 0** (~1.6 s).
- np=4 dry-run: built, rank-0 getters `6378 / 6378 / 2126 / 3`, owned/local split correct (`global_owned=27936` vs `global_local=27945`), **exit 0**.
- Heterogeneous request (`use_sidecar=true` + dataset_root): the Phase 2b abort fires with the documented message, exit 1.
- No-regression: `seas_test_spatial_seas_config` 40/40; `seas_test_spatial_friction_config` 29 + same pre-existing NUM-1 abort; `spatial_dyn_driver` rebuilds clean; QD warning does not leak into it.

## Findings

### [R-201] LOW [POSSIBLE] [spatial_seas_driver.cpp:4.3] — `[boundary]` fallback defaults are BP5-specific and silently apply to any config that omits a key

**Category:** EDGE_CASE / ASSUMPTION

**Description:**
The bc construction falls back to the BP5 convention per-field: `fault_attr→3`, `natural_attrs→{1}`, `dirichlet_attrs→{5}`. These are applied **independently** — so a partially-specified `[boundary]` block (e.g. a SAF-style config that sets `fault_attr=101`, `natural_attrs=[102]`, `absorbing_attrs=[103,104]` but no `dirichlet_attrs`) silently gets `dirichlet_attrs={5}`. If attr 5 is absent from that mesh, `BuildFacetBCTables`/`SetupBoundaryMarkers` aborts with a confusing "Dirichlet attr 5 not in mesh" rather than a clear "QD requires `[boundary].dirichlet_attrs`". Benign for BP5 (the Phase 2 target, which has attr 5), but a foot-gun for the SAF configs Phase 8 will add.

**Trigger:** a `[boundary]` block that sets fault/natural but omits `dirichlet_attrs`, on a mesh without attr 5.

**Actual behavior:** silent BP5 fallback `{5}` → downstream abort with a mesh-attr message.

**Expected behavior:** if `[boundary]` is present but `dirichlet_attrs` is empty, abort with a QD-specific message (Dirichlet plate loading is mandatory); only use the all-BP5 fallback when `[boundary]` is entirely absent.

**Suggested fix:**
```diff
-   bc.dirichlet_attrs = !cfg.boundary.dirichlet_attrs.empty()
-                        ? std::set<int>(cfg.boundary.dirichlet_attrs.begin(),
-                                        cfg.boundary.dirichlet_attrs.end())
-                        : std::set<int>{5};
+   const bool has_boundary_block = (cfg.boundary.fault_attr > 0);
+   if (has_boundary_block)
+   {
+      MFEM_VERIFY(!cfg.boundary.dirichlet_attrs.empty(),
+                  "spatial_seas_driver: [boundary] is present but "
+                  "dirichlet_attrs is empty — the quasi-dynamic elasticity "
+                  "problem requires far-field Dirichlet plate-loading walls.");
+      bc.dirichlet_attrs = std::set<int>(cfg.boundary.dirichlet_attrs.begin(),
+                                         cfg.boundary.dirichlet_attrs.end());
+   }
+   else { bc.dirichlet_attrs = std::set<int>{5}; }  // no [boundary] => BP5 default
```

**Test case:**
```cpp
// test_R201_boundary_without_dirichlet_aborts_clearly:
//   config with [boundary] fault_attr=101 natural_attrs=[102] (no dirichlet_attrs)
//   EXPECT a clear "requires far-field Dirichlet" abort, not a mesh-attr-5 abort.
```

---

### [R-202] LOW [spatial_seas_driver.cpp:4.4] — the `MaterialField` object + `Mode::Constant` assertion are dead/vacuous

**Category:** QUALITY

**Description:**
The heterogeneous gate is `hetero_requested` (config flags), which aborts before the `MaterialField` is built. The subsequently-built `material = MaterialField::MakeConstant(...)` is therefore **always** `Mode::Constant`, so `MFEM_VERIFY(material.mode == ...Constant)` can never fire and `material` is otherwise unused (the operator consumes `le`, the `LinearElastic`). The plan said "build `MaterialField`, for `Mode::Constant` extract λ,μ", but with the flag-based gate the object is redundant. Harmless, but it reads as load-bearing when it isn't.

**Trigger:** always.

**Suggested fix:** drop the redundant object, or keep only as documentation:
```diff
-   const real_t mat_lambda = cfg.material_fallback.lambda;
-   const real_t mat_mu     = cfg.material_fallback.mu;
-   MaterialField material = MaterialField::MakeConstant(
-      mat_lambda, mat_mu, cfg.material_fallback.rho);
-   MFEM_VERIFY(material.mode == MaterialField::Mode::Constant,
-               "spatial_seas_driver: Phase 2 requires MaterialField::Mode::"
-               "Constant.");
-   LinearElastic le(mat_lambda, mat_mu);
+   // hetero_requested (above) already guarantees the constant path.
+   const real_t mat_lambda = cfg.material_fallback.lambda;
+   const real_t mat_mu     = cfg.material_fallback.mu;
+   LinearElastic le(mat_lambda, mat_mu);
```

---

### [R-203] LOW [spatial_seas_driver.cpp:4.1/4.4] — heterogeneous-material abort runs AFTER the (potentially huge) parallel mesh load

**Category:** QUALITY / EDGE_CASE

**Description:**
The `hetero_requested` check (4.4) is after the mesh load + partition (4.1). A heterogeneous config therefore loads and partitions the full SAF/BP5 mesh, then aborts. The check is config-only (no mesh dependency), so it can fail fast before the expensive load.

**Trigger:** a heterogeneous config on a large mesh.

**Suggested fix:** move the `hetero_requested` `MFEM_VERIFY` above the `Mesh smesh(...)` load (it depends only on `cfg`).

---

### [R-204] LOW [spatial_seas_driver.cpp:getter print] — rank-0 print shows rank-0-LOCAL counts, labelled as if absolute

**Category:** QUALITY

**Description:**
`GetNumFaultDOFs()` etc. return **per-rank local** values (np=4 rank-0 printed 6378 while the global owned count is 27936). Printing them only on rank 0 with bare labels (`GetNumFaultDOFs: 6378`) reads as a global total. The plan's acceptance just says "prints" the getters, so this satisfies it literally, but it can mislead. (The operator's own `Owned fault layout` audit already prints the globals.)

**Suggested fix:** label them per-rank, or print the rank too:
```diff
-      std::cout << "[spatial_seas] elasticity domain operator constructed:\n"
-                << "  GetNumFaultDOFs:       " << domain.GetNumFaultDOFs() << "\n"
+      std::cout << "[spatial_seas] elasticity domain operator constructed "
+                << "(rank 0 local counts):\n"
+                << "  GetNumFaultDOFs:       " << domain.GetNumFaultDOFs() << "\n"
```
(Note: the getters are pure member accessors — calling them on rank 0 only is MPI-safe; verified no hang at np=4. Not a correctness issue.)

---

### [R-205] LOW [POSSIBLE] [spatial_seas_driver.cpp:4.5] — `DomainConfig` only sources `check_residual` + `blr_tol`; `face_basis_type` / `penalty_factor` / `match_quad_order` are left at defaults (BP5-parity relevant)

**Category:** ASSUMPTION (forward-looking, Phase 7)

**Description:**
`seas_driver`'s `BuildDomainConfig` sets `face_basis_type`, `penalty_factor`, `blr_tol`, `check_residual` from the solver config. The QD `SolverSpec` (Phase 1) has only `residual_check` + `blr_tol`, so this driver leaves `face_basis_type`, `penalty_factor`, `match_quad_order` at their `DomainConfig` defaults. For Phase 2 (construct-without-abort) this is fine and the defaults happen to match BP5. But if BP5 parity (Phase 7) needs a non-default `penalty_factor` or `face_basis_type`, the QD path cannot express it and parity would silently differ. Not a Phase 2 bug; flagged so Phase 7 either confirms the defaults match or adds the fields to `SolverSpec`.

**Suggested fix:** none required for Phase 2. Track for Phase 7: confirm BP5 uses the `DomainConfig` defaults for `penalty_factor`/`face_basis_type`, or add the missing knobs to `SolverSpec`.

---

### [R-206] LOW [POSSIBLE] [spatial_seas_driver.cpp:4.3] — `absorbing_attrs` forwarded to a quasi-static operator (a wave-only concept)

**Category:** ASSUMPTION

**Description:**
`bc.absorbing_attrs` is populated from `[boundary].absorbing_attrs`. Absorbing BCs model outgoing waves — meaningless in the quasi-static QD elasticity solve. The elasticity setup (`SetupBoundaryMarkers`/`BuildFacetBCTables`) only marks Dirichlet + fault faces, so absorbing (and natural) attrs are effectively unmarked → treated as traction-free (natural). That is harmless for BP5 (empty absorbing list), but a SAF config that lists absorbing walls would have them silently become traction-free rather than the intended far-field condition.

**Trigger:** a QD config with non-empty `[boundary].absorbing_attrs` (none today; Phase 8 SAF).

**Suggested fix:** none for Phase 2. For Phase 8, either reject `absorbing_attrs` in the QD driver (they are dynamic-only) or document that QD maps them to natural.

---

## Observation (NOT a finding against Phase 2 code)
- The operator prints a **"Y=0 Face Classification Audit" warning** ("52168 y=0 faces have NO BC … conflicting with fault slip") on `bp5_tandem_exact.msh`. This is emitted by the shared `ElasticityDomainOperator` setup with exactly the bc the BP5 `seas_driver` verification uses (fault=3 / natural=1 / dirichlet=5) — it is not introduced by Phase 2 and does not abort. It must be reconciled at Phase 7 (BP5 parity): confirm the BP5 golden run produces the same audit and that it is benign for the bounded-fault-patch BP5 geometry.

## Summary
- Critical issues: **0**
- Moderate issues: **0**
- Low issues: 6 (R-201 boundary fallback; R-202 dead MaterialField; R-203 abort ordering; R-204 local-count labelling; R-205 DomainConfig parity fields; R-206 absorbing in QD)
- Plan compliance: **FULL** for Phase 2 —
  - mesh load (3D assert) ✔; `BoundaryConfig` with Dirichlet plate loading ✔; homogeneous `LinearElastic` + heterogeneous abort ✔; `ElasticityDomainOperator` via the line-106 ctor ✔; dry-run prints the four getters ✔; BP5-mesh construction np=1 + np=4 without abort ✔.
  - Documented deviations (added `BoundarySpec.dirichlet_attrs`; `Wf`/`lf` BP5 placeholders; `dg_method=IP`; dry-run now constructs) are all justified and verified dyn-safe (the dynamic driver rebuilds clean and never reads `dirichlet_attrs`).
- Verdict: **PASS** — no CRITICAL/MODERATE issues; construction verified on the BP5 mesh in serial and parallel. The 6 LOW items are robustness/clarity hardening and Phase-7/Phase-8 forward concerns; none block Phase 3.

## Unreviewed Areas
- Physics correctness of the assembled stiffness / BP5 traction (Phase 7 parity) — out of scope; the Y=0 audit warning is flagged for Phase 7.
- The `Wf`/`lf`/`dg_method`/`penalty_factor` config-sourcing gap — functionally inert in Phase 2, relevant to later phases.
- Full `make all && make test` not run (heavy; non-interference verified by rebuilding the dynamic driver + both config tests).
