# Code Review: SAFS QD setup (config + jobs + driver guard + SRW plan) — 2026-06-20

## Review Scope
- Plan: `document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md` (§Phase 8)
  + the new `PLAN_spatial_seas_srw_qd_2026-06-20.md`.
- Files reviewed (this session's changes):
  - `config/safs_qd/safs_qd_500m_rssrw_v1.toml`
  - `drivers/spatial_seas_driver.cpp` (the SRW abort guard, ~line 328)
  - `jobs/safs_qd/safs_qd_500m_rssrw_v1_dev_2hr.sbatch`
  - `jobs/safs_qd/safs_qd_500m_rssrw_v1_production.sbatch`
  - `document/spatial_seas_dev/PLAN_spatial_seas_srw_qd_2026-06-20.md`
- Domain context: `CLAUDE.md` (BP5 loading/sign rules, "no silent fallback", extreme-care files),
  `domain/boundary_config.hpp::MakeBP5DirichletFunc`, `domain/elasticity_operator*.inl`,
  the dynamic mirror `spatial_friction_safs_seisol_v1_0_0_RSSRW.toml`.
- Method: three adversarial passes; verified the loading function body, the constant-tensor
  sign convention, the depth-profile CSVs, and the order>=2 mechanics against source.

## Bonus question — "can the quasi-dynamic code raise polynomial order >= 2?"

**YES.** The QD elasticity + fault stack is order-parametric, not hardcoded to p1:
- FE space: `fec_ = DG_FECollection(order_, 3, GaussLobatto)` (`elasticity_operator_setup.inl:11`);
  `order_` is `[mesh].order` from config. No `order==1` assertion anywhere in the stack.
- DG penalty scales with order: `c_N_1 = order*(order+dim-1)/dim` (`seas_driver.cpp:426`).
- Per-face fault DOF count is order-derived: `nbf_per_face_ = face_quad_->NumBasisFunctions()`
  (`elasticity_operator_setup.inl:660`), so the resolver/per-DOF tables scale automatically.
- **Exercised in production:** the native QD BP5 driver has p2/p4/p6 jobs —
  `jobs/bp5/bp5_v34_1000m_br2_p2_full.sbatch`, `..._v37_2500m_br2_p4_...`,
  `..._v36_4000m_br2_p6_...`, `..._v46c_1000m_ip_p4_...`. `spatial_seas_driver` shares this
  exact `ElasticityDomainOperator`/`RateStateFaultOperator` stack, so it inherits order>=2.
- To run the SAF QD at higher order: set `[mesh].order = 2` (or more) in
  `safs_qd_500m_rssrw_v1.toml`. Caveats: cost grows steeply (~order^3 DOFs in 3D) — the
  `cg_amg` iterative solver matters more at p>=2; CFL/`dt` and the IP penalty change with
  order (memory: MFEM SAFS dt ~4.2x SeisSol is partly p1-vs-p3); BP5 parity (Phase 7) and
  the SAF smoke should be re-confirmed at the chosen order before production.

## Findings

### [R-001] CRITICAL [config + driver] safs_qd_500m_rssrw_v1.toml / MakeBP5DirichletFunc — plate loading DEGENERATES to uniform on the UTM SAF mesh

**Category:** BUG

**Description:**
The QD far-field plate loading is `MakeBP5DirichletFunc` (`domain/boundary_config.hpp:66-83`):
```cpp
real_t y = x(1);
real_t Vh = Vp * t;
if (y > 1000.0)       { Vh *= 0.5; }   // +Vp t/2
else if (y < -1000.0) { Vh *= -0.5; }  // -Vp t/2
// |y| <= 1000: Vh *= 1   (nucleation-zone special case)
```
It keys the antisymmetric +/- sign off the **absolute global** y-coordinate with a +/-1000 m
threshold (BP5's fault sits at y=0). The SAF fault-box mesh is in **UTM 11N** coordinates
(hypocenter (606971, 3707270); the dynamic config nucleates/projects stress at those UTM
coords, so the mesh vertices are UTM). Every mesh y is ~3.7e6 m >> 1000, so the branch
`y > 1000.0` is ALWAYS taken: u_X = +Vp*t/2 UNIFORMLY on every Dirichlet wall. There is NO
antisymmetric loading — the box gets a rigid +x translation, the fault sees ZERO differential
far-field plate motion, no interseismic stress accumulates, and the cycle model never loads.
My config header called this a "planar sgn(Y) approximation the user accepted"; that
understates it — on UTM coordinates it is not an approximation, it produces no loading at all.
This is the Phase-8 open-question-#1 "BLOCKING" item, and it is genuinely blocking.

**Trigger:** run `safs_qd_500m_rssrw_v1.toml` (or any QD config) on a mesh whose y-coordinates
do not straddle 0 (the SAF UTM mesh).

**Actual behavior:** uniform u_X = +Vp*t/2 on all Dirichlet walls → no fault loading.

**Expected behavior:** +Vp*t/2 on the "front" wall and -Vp*t/2 on the "back" wall (differential
plate motion across the fault), as the user requested.

**Suggested fix:** the config alone cannot fix this — it needs a SAF loading function keyed off
a fault-relative (or box-center-relative) coordinate. Add a driver-level loading function, e.g.
center y on the box midpoint y0:
```cpp
// new: MakeSAFDirichletFunc(Vp, y0) with y0 = 0.5*(ymin+ymax) of the mesh (or the fault trace y)
real_t yr = x(1) - y0;
if (yr > 0.0) Vh *= 0.5; else Vh *= -0.5;   // antisymmetric about the box center
```
and wire it in `spatial_seas_driver.cpp` (replace `MakeBP5DirichletFunc(plate_rate)` at line 525
for SAF runs; compute y0 from `pmesh` bounding box). Until then, mark the SAF QD loading as
unimplemented (NOT merely "approximate") in the config + the SRW plan's Phase S4.

**Test case:**
```python
def test_R001_utm_loading_is_degenerate():
    # MakeBP5DirichletFunc semantics: absolute-y threshold.
    def bp5_dirichlet(y, Vp=2e-9, t=1.0):
        Vh = Vp * t
        if y > 1000.0:   Vh *= 0.5
        elif y < -1000.0: Vh *= -0.5
        return Vh
    # UTM coords: two opposite walls of the SAF box
    front = bp5_dirichlet(3.71e6)   # +y wall
    back  = bp5_dirichlet(3.70e6)   # -y wall (still ~3.7e6 in UTM)
    assert front == back            # BUG: both +Vp/2 -> NO antisymmetry (should be opposite)
```

---

### [R-002] MODERATE [jobs] safs_qd_500m_rssrw_v1_production.sbatch — empty-string argv when SAFS_QD_RESTART_ARGS is unset

**Category:** BUG

**Description:**
The production job ends with `"${SAFS_QD_RESTART_ARGS:-}"` as a final argument to the binary.
When the variable is unset (the normal first-run case) this expands to a single EMPTY-STRING
argv element (`""`), not "no argument". The driver's CLI parser then receives an empty token,
which mfem `OptionsParser` (or the driver's own parse) may reject ("unknown option") or
mis-handle as a stray positional, aborting an otherwise-valid production run.

**Trigger:** submit the production job without exporting `SAFS_QD_RESTART_ARGS`.

**Actual behavior:** an empty-string argument is passed to `seas_spatial_seas_driver`.

**Expected behavior:** no extra argument when no restart is requested.

**Suggested fix:** use a bash array and expand it only when non-empty:
```diff
-CKPT_EVERY="${SAFS_QD_CKPT_EVERY:-2000}"
+CKPT_EVERY="${SAFS_QD_CKPT_EVERY:-2000}"
+RESTART_ARGS=()
+if [[ -n "${SAFS_QD_RESTART:-}" ]]; then
+    RESTART_ARGS=(--restart "${SAFS_QD_RESTART}")
+fi
@@
 ibrun "${BINARY}" \
     --config "${CONFIG_TOML}" \
     --output-dir "${OUT_DIR}" \
     --checkpoint-every "${CKPT_EVERY}" \
-    "${SAFS_QD_RESTART_ARGS:-}"
+    "${RESTART_ARGS[@]}"
```

**Test case:**
```bash
# test_R002: SAFS_QD_RESTART unset -> no empty arg reaches the binary
unset SAFS_QD_RESTART
RESTART_ARGS=(); [[ -n "${SAFS_QD_RESTART:-}" ]] && RESTART_ARGS=(--restart "${SAFS_QD_RESTART}")
printf '%s\n' "${#RESTART_ARGS[@]}"   # must print 0 (no stray "")
```

---

### [R-003] MODERATE [POSSIBLE] [config] safs_qd_500m_rssrw_v1.toml — plate_rate_vp=2e-9 may double the intended loading rate

**Category:** ASSUMPTION

**Description:**
The user said "plate loading Vp/2 = 1e-9 on front and back side". I read this as per-side
velocity = +/-1e-9, hence full plate rate Vp = 2e-9, so `plate_rate_vp = 2.0e-9` (the driver
applies u_X = sgn*Vp*t/2 -> +/-1e-9 per side). BUT the standard SEAS/BP5 convention is that
`Vp` IS the plate-rate parameter (= 1e-9 m/s in BP5) and the loading already imposes +/-Vp/2
per side. Under that reading the user meant `plate_rate_vp = 1e-9` and my value is 2x too
fast: Vp=2e-9 m/s ~ 63 mm/yr, whereas the real SAF is ~25-34 mm/yr (Vp~1e-9 ~ 31.5 mm/yr is
the physically standard choice). The interpretation changes the loading rate (and hence
recurrence) by 2x.

**Trigger:** any run with this config.

**Actual behavior:** per-side loading +/-1e-9 m/s (relative plate rate 63 mm/yr).

**Expected behavior:** confirm with the user — likely `plate_rate_vp = 1.0e-9` (per-side
+/-0.5e-9, relative 31.5 mm/yr) if they meant the SCEC `Vp` convention.

**Suggested fix:** confirm intent; if SCEC convention,
```diff
-plate_rate_vp = 2.0e-9
+plate_rate_vp = 1.0e-9   # SCEC Vp convention: per-side +/-0.5e-9 m/s (~31.5 mm/yr)
```
(documented both readings in the config; this finding flags that the default chosen may be 2x.)

**Test case:** N/A (spec-intent, not a code path) — resolve by user confirmation.

---

### [R-004] MODERATE [POSSIBLE] [config/driver] safs_qd_500m_rssrw_v1.toml — constant-material path (no [material] block) UNVALIDATED for the QD driver

**Category:** ASSUMPTION

**Description:**
The config omits a heterogeneous `[material]` block (relying on `[material_constant_fallback]`
+ `[velocity].use_sidecar=false`) to get a constant LinearElastic medium — mirroring the
dynamic v1.0.0 config. But local validation only reached config PARSE: the SRW guard aborts at
driver line ~348, BEFORE Phase-2 material construction, and the only config that constructed
through Phase 2/3 locally was a BP5 config (different material block). So whether the QD driver
correctly builds a constant `LinearElastic(lambda,mu)` from `[material_constant_fallback]` when
NO `[material]` block is present is unverified for this config. If the QD driver requires a
`[material]` block (the skeleton had one) it would abort.

**Trigger:** run the config past the SRW guard (i.e., after the SRW wiring lands, or with an
aging-law variant) on any mesh.

**Suggested fix:** add an aging-law smoke variant of this config (state_evolution="aging_law")
and dry-run it on a small mesh to confirm the constant-material path constructs; OR confirm in
`spatial_seas_driver.cpp` (Phase 2 material section, ~lines 505-586) that a missing `[material]`
falls back to constant. If the driver needs an explicit block, add `[material] kind="constant"`.

**Test case:**
```
// test_R004: dry-run an aging-law copy of safs_qd_500m_rssrw_v1.toml on bp5_1000m.msh;
//   assert it constructs ElasticityDomainOperator with a constant LinearElastic (no abort).
```

---

### [R-005] LOW [config] safs_qd_500m_rssrw_v1.toml — sigma_n_strength_floor omitted vs the dynamic reference

**Category:** ASSUMPTION

**Description:**
The dynamic v1.0.0 mirror sets `[friction] sigma_n_strength_floor_pa = 10.0e6` (an RS-case
stabilizer keeping the strength's effective normal stress off zero). My QD config omits it,
keeping only `[pore_pressure].min_sigma_n_pa = 1.0e6`. If the QD path honors a strength floor,
omitting it could let sigma_n_eff dip near 1 MPa during a cycle and stress the Brent solver.
LOW because it is unclear the QD operator reads `[friction].sigma_n_strength_floor_pa` at all.

**Suggested fix:** confirm whether the QD path uses it; if so, add
`[friction]\nsigma_n_strength_floor_pa = 10.0e6` to match the reference; if not, document that
the QD model relies on `min_sigma_n_pa`.

---

### [R-006] LOW [POSSIBLE] [config] safs_qd_500m_rssrw_v1.toml — bottom (tag 103) as Natural may under-constrain the deep box

**Category:** ASSUMPTION

**Description:**
`natural_attrs = [102, 103]` makes BOTH the z=0 top AND the bottom traction-free, matching
BP5's "Natural on top/bottom, Dirichlet on vertical walls". For a QD SAF box that is correct
for BP5 but a traction-free bottom can leave the deep domain weakly constrained (near-rigid
drift) once the loading function is fixed (R-001). Worth confirming the bottom should be
Natural vs Dirichlet-loaded for the SAF setup.

**Suggested fix:** if drift appears at init/equilibrium, move 103 to `dirichlet_attrs` (loaded)
or pin it; otherwise keep Natural and document the choice.

---

## Summary
- Critical issues: **1** — R-001 (UTM loading degeneracy: the SAF QD fault does NOT load; the
  config's "planar approximation" framing is wrong — it is unimplemented).
- Moderate issues: **3** — R-002 (empty-arg in production job), R-003 (plate-rate 2x intent),
  R-004 (constant-material path unvalidated for this config).
- Low issues: **2** — R-005 (sigma_n floor omitted), R-006 (bottom Natural).
- Plan compliance: **PARTIAL** — the SAFS *friction/stress problem setup* (the v1.0.0 RSSRW
  adoption) is correct and faithful; the *loading* setup is blocked by R-001 (a driver change,
  not a config one); the SRW *engine* is correctly deferred to a plan + abort guard.
- Verdict: **PASS WITH FIXES** — R-001 must be fixed (driver loading function) before any SAF QD
  run can load the fault; R-002 before the first production submit; R-003/R-004 need user/dry-run
  confirmation. The SRW guard + plan + the friction/stress numbers are correct.

## Unreviewed Areas
- The SRW-in-QD wiring itself (not implemented — only planned; `PLAN_spatial_seas_srw_qd_2026-06-20.md`).
- The actual SAF mesh boundary tags (104=all sides vs split) — Frontera-only (mesh gitignored);
  the driver aborts listing present attributes on mismatch, so a wrong tag fails loudly.
- BP5 parity + SAF smoke RUNS — Frontera ([[feedback-no-local-mesh-runs]]).
- Order>=2 for the SAF QD specifically (the BP5 p2/p4/p6 jobs prove the stack; the SAF config
  uses order=1 and a higher-order SAF run is untested).
