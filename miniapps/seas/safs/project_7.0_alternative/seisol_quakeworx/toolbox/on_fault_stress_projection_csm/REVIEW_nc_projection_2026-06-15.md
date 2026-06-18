# Code Review: CSM C1 k=2.39 on-fault stress projection -> ASAGI nc (2026-06-15)

Scope of the question reviewed: **"can we confirm the projection is correct?"**
i.e. the pipeline that produces `safs_seisol_v3_0_0_LSW/safs_stress_csm.nc`
(the field SeisSol actually reads), and the VTU/self-check that validates it.

(IDs use the `P-00x` prefix so they do not collide with the existing
`REVIEW.md`, which audits the separate C2/Luttrell magnitude path and still
has open items R-002..R-004. That file is preserved, not overwritten.)

## Review Scope
- Plan: `PLAN_csm_stress_magnitudes_2026-06-12.md`
- Files reviewed:
  - `csm_stress_to_asagi.py` (the nc builder + self-check) — full
  - `project_csm_stress_to_vtu.py` (the VTU ground truth: `read_csm_csv`,
    `csm_axes_and_shape`, `magnitudes_C1`, `build_tensor_from_axes`,
    `tandem_basis`, `harmonise_normals`, `resolve_tractions`, `main` depth-R) — full
  - Artifacts: `safs_stress_csm.nc`, `*_fault_stress.vtu` (4 variants),
    `*_summary.json`, the PUML mesh fault facets
- Domain context: project `CLAUDE.md`, `miniapps/seas/CLAUDE.md` (sign
  conventions), the plan, the existing `REVIEW.md`.

## What is CONFIRMED correct (positive findings, with evidence)

1. **Orientation** is validated against an independent source (the CSV's own
   V-columns and SHmax): `eigvec_dot_V3_compressional` median 0.99999,
   `eigvec_dot_V1_tensional` median 0.99999, `shmax_recon_vs_model_deg`
   median 0, max 7.1e-14 deg. The eigendecomposition reproduces the CSM axes.
2. **R-convention is self-consistent.** `csm_axes_and_shape` returns
   `R_eig = (w1-w0)/(w2-w0)` = tension-positive phi = (S2-S3)/(S1-S3); feeding
   it to `magnitudes_C1` (`sig3 = Sv_eff/((1-R)k+R)`) and recovering R from the
   built compression-positive tensor gives `R_recovery_max_abs_err = 4.05e-15`.
   The tension->compression sign/sort flip preserves the shape ratio (verified
   algebraically: (sig1-sig2)/(sig1-sig3)_comp == (S2-S3)/(S1-S3)_tension).
3. **Traction projection is standard and sign-correct.** `resolve_tractions`
   computes `t = sigma.n`, `sigma_n_eff = n.sigma.n`, `tau = |shear|`,
   `mu = tau/sigma_n_eff` (sigma_n_eff>0). `sigma_n_eff` and `tau_mag` are
   invariant to the normal-flip and basis-sign, so `harmonise_normals` /
   `tandem_basis` sign choices cannot corrupt the scalar fields.
4. **Coverage + finiteness of the production nc (checked directly):** all 6
   fields are 100% finite (0 NaN over 3,315,519 nodes); the grid
   x[350,630]km y[3680,3850]km z[-17000,0]m fully contains the fault (deepest
   facet centroid -16528 m > nc floor -17000 m), so **0 facets fall back to the
   ConstantMap tensor** and 0 are clamped. k_ratio attr = 2.39.
5. **Self-check ground truth is the right artifact:** the no-suffix VTU is
   `closure=ratio, k=2.39, axes=csm, sampling=linear` (matches the nc), not a
   stale CONSTANT/C2 file (the 4 VTUs are distinct by md5).
6. **The hypocenter / fault-wide numbers downstream code relies on are real
   projection outputs:** hypocenter sigma_n_eff 70.14 MPa, tau 20.77,
   mu_app 0.296; fault-wide mu_app max 0.4494 (this is what pins mu_s=0.47).

**Bottom line up front:** the projection that produced `safs_stress_csm.nc` is
correct and well-validated for the production artifact. The findings below are
(a) gaps in the *verification harness* that make "SELF-CHECK PASS" weaker than
it appears, and (b) one physical-fidelity caveat. None of them are an active
error in the current nc — but P-001/P-002 mean a *future* regeneration could
ship a corrupt nc that still passes, and P-003 is a modelling assumption the
user should consciously accept.

## Findings

### [P-001] [MODERATE] [csm_stress_to_asagi.py:main self-check] — Self-check silently EXCLUDES non-finite / out-of-grid facets, so a corrupt nc can still PASS

**Category:** BUG (verification gap)

**Description:**
The self-check builds the `fin` mask as
`np.isfinite(r["sigma_n_eff"]) & np.isfinite(r["tau_magnitude"]) & (truth>1.0)`
and then computes only the **median** rel err over `fin`. Two consequences:
(1) any facet where the nc produced NaN/Inf is *dropped from the statistic*,
not flagged — so a NaN-corrupted nc passes; (2) facets whose centroid is below
the nc z-floor (or outside x/y) are **clamped** to the edge slice by
`trilinear_sample` (the `locate` helper clips `w` to [0,1] and `i` to
`len-2`), producing a plausible-but-wrong value that is *included* in the
median but, being a minority, does not move it. The gate the user trusts
("SELF-CHECK: PASS") therefore does not guarantee a finite, fully-covering nc.

**Trigger:** Regenerate with a grid whose z-floor is shallower than the fault
(e.g. `--zmin -12000`), or with any degenerate-R grid node (-> NaN, see P-002).

**Actual behavior:** median rel err stays ~1e-3 -> `PASS`, even though a chunk
of the fault is clamped/NaN.

**Expected behavior:** FAIL if any facet is non-finite or out-of-grid.

**Suggested fix:**
```diff
     smp = trilinear_sample(args.out, cent[:, 0], cent[:, 1], cent[:, 2])
+    # explicit coverage check: clamp = facet outside the nc box (silently
+    # snapped to an edge slice by trilinear_sample) -> the nc is NOT the CSM
+    # field there; in SeisSol these fall back to ConstantMap.
+    import numpy as _np
+    with __import__("netCDF4").Dataset(args.out) as _ds:
+        _xg = _np.asarray(_ds["x"][:]); _yg = _np.asarray(_ds["y"][:])
+        _zg = _np.asarray(_ds["z"][:])
+    oob = ((cent[:, 0] < _xg.min()) | (cent[:, 0] > _xg.max())
+           | (cent[:, 1] < _yg.min()) | (cent[:, 1] > _yg.max())
+           | (cent[:, 2] < _zg.min()) | (cent[:, 2] > _zg.max()))
+    if oob.any():
+        sys.exit(f"ERROR: {int(oob.sum())} fault facets are outside the nc "
+                 f"grid box (deepest {cent[:,2].min():.0f} m vs z-floor "
+                 f"{_zg.min():.0f} m); widen --zmin/--xmin/... so the grid "
+                 "covers the whole fault, or they will read ConstantMap.")
     sigma_f = np.zeros((len(tris), 3, 3))
```
and gate the PASS on finiteness + p95 (see P-005), not just the median:
```diff
-    fin = (np.isfinite(r["sigma_n_eff"]) & np.isfinite(r["tau_magnitude"])
-           & (truth["sigma_n_eff_MPa_cell"] > 1.0))
+    n_bad = int(np.sum(~np.isfinite(r["sigma_n_eff"])
+                       | ~np.isfinite(r["tau_magnitude"])))
+    if n_bad:
+        sys.exit(f"ERROR: {n_bad} facets resolve to non-finite traction from "
+                 "the nc (NaN/Inf in the written field) — see P-002.")
+    fin = truth["sigma_n_eff_MPa_cell"] > 1.0
```

**Test case:**
```python
def test_P001_selfcheck_fails_on_uncovered_fault(tmp_path):
    # Build an nc with --zmin -12000 (shallower than the fault's -16.5 km) and
    # assert the converter exits non-zero (coverage), not "SELF-CHECK PASS".
    import subprocess, sys
    r = subprocess.run([sys.executable, "csm_stress_to_asagi.py",
                        "--csm-csv", CSV, "--material-nc", MAT, "--mesh", MESH,
                        "--verify-vtu", VTU, "--out", str(tmp_path/"s.nc"),
                        "--zmin", "-12000"])
    assert r.returncode != 0
```

---

### [P-002] [MODERATE] [csm_stress_to_asagi.py:main grid loop] — No finiteness guard before writing the nc; a degenerate-R grid node writes NaN stress that SeisSol would consume

**Category:** ASSUMPTION / BUG (defensive)

**Description:**
`R_xy = np.clip(R_eig, 0.0, 1.0)` keeps `R_eig = NaN` (np.clip passes NaN
through), and `csm_axes_and_shape` returns NaN where the eigenvalue gap
collapses (`denom <= EPS`). `magnitudes_C1(Sv_eff, NaN, k)` -> NaN ->
`build_tensor_from_axes` -> NaN stress, which `write_stress_asagi` writes with
no check. The current YHSM grid happens to be non-degenerate (verified: 0 NaN),
but a different CSM export, a coarser grid hitting a near-isotropic node, or
the NN-hull fallback landing on a flat tensor would silently ship NaN in the
nc. ASAGI does not reject NaN -> SeisSol reads NaN normal stress -> the run
NaNs out (exactly the kind of silent corruption we just spent effort ruling out
at load time).

**Trigger:** Any grid node with `w1-w0 <= EPS` or `w2-w1 <= EPS` after the
(x,y) tensor interpolation.

**Actual behavior:** NaN written to nc, no error; self-check excludes it (P-001).

**Expected behavior:** abort the build with the offending count.

**Suggested fix:** after the grid loop, before `write_stress_asagi`:
```diff
+    n_nan = sum(int(np.sum(~np.isfinite(comps[f]))) for f in STRESS_FIELDS)
+    if n_nan:
+        bad_xy = int(np.sum(~np.isfinite(R_xy)))
+        sys.exit(f"ERROR: {n_nan} non-finite stress samples "
+                 f"({bad_xy} grid (x,y) nodes have a degenerate CSM tensor "
+                 "-> NaN R/axes). Refusing to write a NaN nc; raise --dx or "
+                 "filter the CSM export.")
     os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
     write_stress_asagi(args.out, gx, gy, gz, comps, attrs)
```

**Test case:**
```python
def test_P002_degenerate_R_aborts(monkeypatch):
    import numpy as np, csm_stress_to_asagi as M
    # force one NaN into the assembled comps and assert the guard fires
    comps = {f: np.zeros((2, 2, 2)) for f in M.STRESS_FIELDS}
    comps["s_xx"][0, 0, 0] = np.nan
    n_nan = sum(int(np.sum(~np.isfinite(comps[f]))) for f in M.STRESS_FIELDS)
    assert n_nan == 1   # the quantity the guard must check (currently unchecked)
```

---

### [P-003] [MODERATE] [POSSIBLE] [project_csm_stress_to_vtu.py:build_tensor_from_axes + magnitudes_C1] — "sig2 = Sv_eff" is the intermediate PRINCIPAL stress, not the vertical stress; sig_zz departs from Sv_eff by up to +/-74 MPa where the CSM intermediate axis plunges

**Category:** ASSUMPTION (physical fidelity)

**Description:**
The docstrings, the nc `attrs`, and `safs_initial_stress.yaml` all describe the
magnitude model as "sig2 = Sv_eff(z) = vertical effective stress". But with
`--axes csm` (the production choice), `magnitudes_C1` assigns `sig2 = Sv_eff` to
the **CSM intermediate eigenvector** (`U_facet[:, :, 1]`), which is vertical
only where the YHSM stress is Andersonian. Where the intermediate axis plunges,
the *vertical* component `sigma_zz` of the built tensor is NOT Sv_eff. The
summary quantifies this directly: `sigma_zz_eff_minus_Sv_eff_MPa` =
median +0.095, **min -53.1, max +73.6 MPa**. So on the tails the lithostatic
anchor the whole closure is built on is violated by tens of MPa, which feeds
straight into the on-fault sigma_n / mu_app / S-ratio (the very quantities used
in the gate and nucleation analyses).

This is a *modelling choice* (Option A vs the `--axes andersonian` Option B that
forces one axis vertical), not a coding error — but it is the single most
important thing to confirm when asking "is the projection correct", and the
"sig2 = Sv_eff" wording oversells it.

**Trigger:** Facets sampling CSM points with a non-trivial V2 (intermediate
axis) plunge — `csm_V2_plunge_deg_cell` far from 90 deg.

**Actual behavior:** sigma_zz != Sv_eff (up to +/-74 MPa) on those facets; the
text claims sig2=Sv_eff is the vertical effective stress.

**Expected behavior:** either (a) explicitly accept Option A and reword the
provenance to "sig2 (intermediate principal) = Sv_eff; vertical stress may
deviate where the intermediate axis plunges", and verify the high-deviation
facets are NOT at the gate / hypocenter / nucleation patch; or (b) use
`--axes andersonian` if "vertical = Sv_eff" is the intended constraint.

**Suggested fix (surface it; choose model separately):** add a warning in the
depth-R `csm` branch right after `sigma_facet` is built:
```diff
         if args.axes == "csm":
             sigma_facet = build_tensor_from_axes(U_facet, sig1, sig2, sig3)
+            dzz = sigma_facet[:, 2, 2] - Sv_eff
+            n_big = int(np.sum(np.abs(dzz) > 0.1 * Sv_eff))
+            if n_big:
+                print(f"warning: {n_big} facets have |sigma_zz - Sv_eff| > "
+                      f"10% of Sv_eff (max {np.nanmax(np.abs(dzz)):.1f} MPa); "
+                      "the CSM intermediate axis is not vertical there, so the "
+                      "'sig2 = vertical Sv_eff' anchor is local only "
+                      "(use --axes andersonian to force vertical).",
+                      file=sys.stderr)
```

**Test case:**
```python
import numpy as np
from project_csm_stress_to_vtu import build_tensor_from_axes
def test_P003_plunging_intermediate_axis_breaks_szz_equals_sv():
    # intermediate axis tilted 30 deg from vertical -> sigma_zz != sig2
    c, s = np.cos(np.radians(30)), np.sin(np.radians(30))
    U = np.array([[[1, 0, 0], [0, c, -s], [0, s, c]]], float)  # (1,3,3)
    sig1 = np.array([85.0]); sig2 = np.array([75.0]); sig3 = np.array([36.0])
    T = build_tensor_from_axes(U, sig1, sig2, sig3)
    assert abs(T[0, 2, 2] - sig2[0]) > 1.0   # vertical stress != Sv_eff anchor
```

---

### [P-004] [LOW] [csm_stress_to_asagi.py:main + project_csm_stress_to_vtu.py:main] — R clipped to the CLOSED interval [0,1]; R in {0,1} collapses principal axes

**Category:** EDGE_CASE

**Description:**
`R_xy = np.clip(R_eig, 0.0, 1.0)` (converter) and
`R_point = np.clip(..., 0.0, 1.0)` (VTU) admit R=0 and R=1. Under C1 these are
finite (no `1/(1-R)` blow-up, unlike the C2 path of prior R-002) but degenerate:
R=1 -> `sig3 = Sv_eff = sig2` (intermediate==least), R=0 -> `sig1 = Sv_eff = sig2`
(greatest==intermediate). The eigenframe of the built tensor is then
ill-defined for the two colliding axes. YHSM R maxes at 0.989 so it is not hit
now, but the col-13/eig fallback can produce exactly 0 or 1.

**Suggested fix (both sites):**
```diff
-    R_xy = np.clip(R_eig, 0.0, 1.0)                      # (nxy,)
+    R_xy = np.clip(R_eig, 1.0e-6, 1.0 - 1.0e-6)          # open: avoid axis collapse
```

**Test case:**
```python
import numpy as np
from project_csm_stress_to_vtu import magnitudes_C1
def test_P004_R_one_collapses_sig2_sig3():
    sig1, sig2, sig3 = magnitudes_C1(np.array([100.0]), np.array([1.0]), 2.39)
    assert abs(sig2[0] - sig3[0]) < 1e-9      # documents the degeneracy clip must avoid
```

---

### [P-005] [LOW] [csm_stress_to_asagi.py:main self-check PASS gate] — PASS gated on the MEDIAN only; a localized wrong region does not fail it

**Category:** BUG (verification gap)

**Description:**
`ok = (np.median(relsn) < tol) and (np.median(reltau) < tol)`. The p95 and max
are printed but not gated, so a coherent patch of the fault (e.g. a few hundred
facets near a CSM-hull edge, or the clamped deep band of P-001) can be 5-10%
wrong while the median stays 1e-3 -> PASS.

**Suggested fix:**
```diff
-    ok = (np.median(relsn) < args.selfcheck_tol
-          and np.median(reltau) < args.selfcheck_tol)
+    ok = (np.percentile(relsn, 95) < args.selfcheck_tol
+          and np.percentile(reltau, 95) < args.selfcheck_tol)
```
(or keep the median gate but add a separate p95 ceiling; tune `--selfcheck-tol`
since p95 is stricter than the median it currently uses.)

**Test case:**
```python
import numpy as np
def test_P005_p95_gate_catches_local_error():
    rel = np.concatenate([np.full(9000, 1e-3), np.full(1000, 0.08)])
    assert np.median(rel) < 0.03           # current gate: PASS (wrong)
    assert np.percentile(rel, 95) > 0.03   # fixed gate: FAIL (correct)
```

---

### [P-006] [LOW] [csm_stress_to_asagi.py:main self-check] — tau rel-err denominator floored at 1.0 MPa makes the self-check lenient in the low-shear gate cores

**Category:** QUALITY

**Description:**
`reltau = |delta| / np.maximum(truth_tau, 1.0)`. Where tau < 1 MPa (the
low-mu_app gate cores, mu_app ~ 0.014 min), the denominator is floored at 1.0,
so a 0.5 MPa absolute discrepancy reads as a 0.5 "rel" error instead of being
measured against the true (sub-MPa) tau. The gate region is exactly where the
S-ratio analysis is most sensitive, so the self-check is weakest where it
matters most. Not wrong, but worth an absolute-tolerance companion.

**Suggested fix:** report and gate an absolute tau error too, e.g.
`abstau = np.abs(r["tau_magnitude"] - truth_tau)` with a small MPa ceiling, in
addition to the floored relative error.

**Test case:** (documentation-level; assert the floor masks a sub-MPa error)
```python
import numpy as np
def test_P006_tau_floor_masks_subMPa_error():
    truth = np.array([0.2]); got = np.array([0.7])   # 0.5 MPa = 250% of true
    reltau = np.abs(got - truth) / np.maximum(truth, 1.0)
    assert reltau[0] < 0.6   # floored: looks like 50%, not 250%
```

---

## Summary
- Critical issues: 0
- Moderate issues: 3 (P-001, P-002, P-003)
- Low issues: 3 (P-004, P-005, P-006)
- Plan compliance: FULL for the C1 k=2.39 production path (orientation
  validated, R recovery 4e-15, ordering/positivity True, nc finite + fully
  covers the fault, self-check vs the correct linear C1 VTU).
- Verdict: **PASS WITH FIXES.** The projection in `safs_stress_csm.nc` is
  correct and validated for the current artifact; the moderate findings are
  (P-001/P-002) verification-harness gaps that could let a *future*
  regeneration ship a clamped/NaN nc that still "passes", and (P-003) a
  modelling assumption ("vertical = Sv_eff") that is local-only and should be
  consciously accepted (check the gate/hypocenter facets' V2 plunge).

## Direct answer: can we confirm the projection is correct?
**Yes, for the production nc — with one caveat to accept.**
- Orientation: independently validated (SHmax recon max 7e-14 deg; eigenvectors
  match the CSV V-columns to 1e-4; R recovery 4e-15).
- Tractions: standard, sign-consistent `sigma_n = n.sigma.n`, `tau = |shear|`,
  `mu = tau/sigma_n`; invariant to the basis/normal sign choices.
- Artifact: 0 NaN, fully covers the fault (0 ConstantMap-fallback facets),
  k=2.39; nc reproduces the C1 linear VTU to median 1.4e-3 (grid-resolution
  error, expected). The hypocenter (70.14 / 20.77 MPa, mu_app 0.296) and
  fault-wide mu_app max (0.4494) the downstream analyses use are genuine
  outputs of this projection.
- Caveat (P-003): the magnitude closure anchors the *intermediate principal*
  stress to Sv_eff, not the vertical stress; sigma_zz deviates from Sv_eff by up
  to +/-74 MPa where the CSM intermediate axis plunges (median ~0). Confirm the
  gate/hypocenter facets have near-vertical intermediate axes before treating
  "sig_n driven by lithostat" as exact there.
- The self-check "PASS" is weaker than it looks (P-001/P-002/P-005): it excludes
  non-finite facets, clamps (not flags) out-of-grid facets, and gates on the
  median only. Those did not bite the current nc, but fix them before the next
  regeneration is trusted blindly.

## Unreviewed Areas
- The C2/Luttrell magnitude path (`magnitudes_C2`, the cap) — covered by the
  existing `REVIEW.md`; not re-audited here (production uses C1).
- `apply_deep_taper` — OFF in production (`taper_zseis_km=None`), not exercised.
- `write_vtu` binary writer — round-trips correctly (the nc self-check reads it
  back); not re-audited.
- The `--axes andersonian` branch — not used by the production artifact (only
  P-003's recommendation references it).
