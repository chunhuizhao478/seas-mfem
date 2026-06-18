# Code Review: depth-R magnitude model in project_csm_stress_to_vtu.py (2026-06-15)

## Review Scope
- Plan: `PLAN_csm_stress_magnitudes_2026-06-12.md` (this folder)
- Files reviewed: `project_csm_stress_to_vtu.py` (depth-R additions; constant
  path verified byte-identical and not re-reviewed)
- Domain context: project `CLAUDE.md`, `miniapps/seas/CLAUDE.md`, the plan
  doc, and the generated artifact `csm_yhsm2013_stress_on_safs_mesh_*`
- User question answered inline: **why mu ≈ 1.0 near the free surface** — see
  R-001 (it is the positivity cap, not a coordinate/sign bug).

## Findings

### [R-001] [CRITICAL] [project_csm_stress_to_vtu.py:magnitudes_C2] — Positivity cap makes the near-surface stress super-critical (mu→1.0)

**STATUS: FIXED (2026-06-15).** `magnitudes_C2` now caps by the
frictional (Byerlee) ceiling `--cap-k-max` (default 5.83) instead of the
`--sig3-floor-frac` floor. Verified on the regenerated VTU: max stress
ratio over all facets = 5.830 (was ~31), mu max 0.92 (was 1.41), 0–500 m
mu max 0.76 (was 1.05); shallowest facet differential tapered from the
raw Luttrell 28.9 MPa to 0.25 MPa; ordering True, R-recovery 1.4e-13,
hypocenter/seismogenic field unchanged. References added to the plan
(§7 refs 23 Sibson 1974, 24 Brace & Kohlstedt 1980, 25 Townend & Zoback
2000; §3/§4 rev 6).

**Category:** BUG (cap design) — this is the answer to "why mu = 1.0 at near free surface"

**Description:**
The closure-C2 positivity cap forces `sig3 = floor_frac * Sv_eff` (default
`floor_frac = 0.05`) on every facet where the Luttrell differential stress
would drive `sig3` negative. That keeps ordering valid (0 < sig3 < sig2) but
pins the *stress ratio* to a fixed, depth-independent value

```
sig1/sig3 = [Sv_eff + R*D_applied] / (floor_frac*Sv_eff)
          = (1/floor_frac) * (1 - floor_frac*R)/(1 - R)
          ≈ 31  for R = 0.37, floor_frac = 0.05
```

A ratio of ~31 is far above the Byerlee frictional ceiling (k = 3.12 at
mu = 0.6, k = 5.83 at mu = 1.0). A stress state that super-critical yields
high apparent friction on the fault. Measured on the produced VTU:

```
depth band   n      mu med   mu max   capped%   sig1/sig3 med
0-500 m     2330    0.511    1.045     84.2%      30.4
500-1000 m  2331    0.437    0.859     58.4%      24.1
1000-2000 m 4503    0.294    1.414     16.7%       2.8
4000-8000 m 20507   0.074    0.358      0.0%       1.2
8000-17000  21595   0.042    0.174      0.0%       1.1
```

So mu ≈ 0.5–1.0 at the free surface is produced *entirely* by the cap: 84%
of the 0–500 m facets are capped to sig1/sig3 ≈ 30. The Luttrell field (D ≈
10–60 MPa, calibrated for ~5 km seismogenic depth) is being imposed on a
near-zero shallow effective overburden (Sv_eff ≈ 0.14 MPa at z = 21 m), which
no real crust could sustain. The physically correct behaviour is that the
sustainable differential stress → 0 as Sv_eff → 0 (frictional strength
vanishes at the free surface), so mu should stay bounded, not climb to 1.0.

**Trigger:**
Any shallow facet where `Sv_eff - (1-R)*D_Luttrell < floor_frac*Sv_eff`
(here, 4165 / 60658 facets, concentrated in the upper ~1 km).

**Actual behavior:**
Capped facets get `sig3 = 0.05*Sv_eff`, fixed ratio sig1/sig3 ≈ 31, and
mu_apparent up to ~1.0 at the free surface. This unphysical near-surface
prestress would seed spurious surface-breaking nucleation in the downstream
dynamic-rupture run.

**Expected behavior:**
Cap the differential stress at the *frictional (Byerlee) limit* instead of a
fixed fraction of Sv_eff, so the stress ratio never exceeds a chosen ceiling
`k_max` and the differential → 0 at the free surface:

```
D_max(z) = (k_max - 1) * Sv_eff(z) / ((1-R)*k_max + R)     [= closure C1 at k_max]
D_applied = min(D_Luttrell, D_max)        => sig1/sig3 <= k_max everywhere
```

This is the plan's own C1 form used as a ceiling; it bounds mu near the
surface and tapers D to zero as Sv_eff → 0.

**Suggested fix:**
Replace the `--sig3-floor-frac` floor with a frictional ratio ceiling
`--cap-k-max` (default 5.83 = Byerlee mu 1.0; use 3.12 for a mu 0.6 ceiling).

In `magnitudes_C2`:
```diff
-def magnitudes_C2(Sv_eff, R, D, floor_frac):
-    """Closure C2 ... cap: where the raw
-    sig3 falls below floor_frac*Sv_eff, reduce D so sig3 = that floor ..."""
-    sig2 = Sv_eff
-    sig3_raw = Sv_eff - (1.0 - R) * D
-    sig3_floor = floor_frac * Sv_eff
-    capped = sig3_raw < sig3_floor
-    one_minus_R = np.maximum(1.0 - R, EPS)
-    D_applied = np.where(capped, (Sv_eff - sig3_floor) / one_minus_R, D)
-    sig3 = np.where(capped, sig3_floor, sig3_raw)
-    sig1 = Sv_eff + R * D_applied
-    return sig1, sig2, sig3, D_applied, capped
+def magnitudes_C2(Sv_eff, R, D, cap_k_max):
+    """Closure C2 (plan Section 4): sig2 = Sv_eff, sig1 = Sv_eff + R*D,
+    sig3 = Sv_eff - (1-R)*D.  Frictional cap: limit D so sig1/sig3 never
+    exceeds the Byerlee ceiling cap_k_max, i.e.
+        D_max = (k-1)*Sv_eff / ((1-R)*k + R)   (= closure C1 at k_max).
+    D_max -> 0 as Sv_eff -> 0, so the near-surface differential vanishes
+    instead of producing a super-critical ratio.  Guarantees
+    0 < sig3 <= sig2 for k_max >= 1, R in [0, 1).
+    Returns (sig1, sig2, sig3, D_applied, capped)."""
+    sig2 = Sv_eff
+    D_max = (cap_k_max - 1.0) * Sv_eff / ((1.0 - R) * cap_k_max + R)
+    capped = D > D_max
+    D_applied = np.where(capped, D_max, D)
+    sig1 = Sv_eff + R * D_applied
+    sig3 = Sv_eff - (1.0 - R) * D_applied
+    return sig1, sig2, sig3, D_applied, capped
```

CLI flag (replace `--sig3-floor-frac`):
```diff
-    ap.add_argument("--sig3-floor-frac", type=float, default=0.05,
-                    dest="sig3_floor_frac",
-                    help="closure diff positivity cap: keep sig3 >= this "
-                    "fraction of Sv_eff (default %(default)s; guarantees "
-                    "0 < sig3 <= sig2 at all depths)")
+    ap.add_argument("--cap-k-max", type=float, default=5.83, dest="cap_k_max",
+                    help="closure diff frictional cap: limit D so the stress "
+                    "ratio sig1/sig3 <= this Byerlee ceiling (default "
+                    "%(default)s = mu 1.0; use 3.12 for mu 0.6). Bounds the "
+                    "near-surface apparent friction and tapers D to 0 at z=0.")
```

Validation guard (replace the `0 < sig3_floor_frac < 1` check):
```diff
-        if not (0.0 < args.sig3_floor_frac < 1.0):
-            sys.exit("ERROR: --sig3-floor-frac must be in (0, 1) "
-                     f"(got {args.sig3_floor_frac})")
+        if args.cap_k_max < 1.0:
+            sys.exit(f"ERROR: --cap-k-max must be >= 1 (got {args.cap_k_max})")
```

Call site + summary param:
```diff
-            sig1, sig2, sig3, D_applied, capped = magnitudes_C2(
-                Sv_eff, R_facet, D_diff, args.sig3_floor_frac)
+            sig1, sig2, sig3, D_applied, capped = magnitudes_C2(
+                Sv_eff, R_facet, D_diff, args.cap_k_max)
```
```diff
-            "sig3_floor_frac": args.sig3_floor_frac,
+            "cap_k_max": args.cap_k_max,
```

**Test case:**
```python
import numpy as np
from project_csm_stress_to_vtu import magnitudes_C2

def test_R001_shallow_cap_is_frictional_not_supercritical():
    # Shallow facet: tiny effective overburden, large Luttrell D
    Sv_eff = np.array([0.14])      # MPa, ~21 m depth
    R      = np.array([0.373])
    D      = np.array([40.0])      # MPa, Luttrell-scale -> must be capped
    cap_k_max = 5.83               # Byerlee mu = 1.0
    sig1, sig2, sig3, D_applied, capped = magnitudes_C2(Sv_eff, R, D, cap_k_max)
    assert capped[0]
    assert sig3[0] > 0.0 and sig3[0] <= sig2[0]
    ratio = sig1[0] / sig3[0]
    # current (floor_frac) code gives ~31; fixed code must respect the ceiling
    assert ratio <= cap_k_max + 1e-6, f"stress ratio {ratio} exceeds k_max"
    # differential must vanish with the overburden (free-surface behaviour)
    Sv_eff0 = np.array([1e-6])
    *_, D_app0, _ = magnitudes_C2(Sv_eff0, R, D, cap_k_max)
    assert D_app0[0] < 1e-3
```

---

### [R-002] [MODERATE] [project_csm_stress_to_vtu.py:main depth-R block] — R is clipped to [0, 1] not [0, 1-eps]; R=1 makes the C2 cap blow up

**Category:** EDGE_CASE

**Description:**
`R_point = np.clip(np.where(np.isfinite(R_csv), R_csv, R_eig), 0.0, 1.0)`
allows R = 1.0 exactly. Both the original cap (`/(1-R)`) and the suggested
C1-ceiling cap (`/((1-R)*k + R)` is fine, but `D_C1` numerator/denominator and
the uncapped `sig3 = Sv_eff - (1-R)*D` degenerate at R→1: sig1 = Sv_eff + D,
sig3 = Sv_eff, so sig1 > sig2 = Sv_eff = sig3 — sig2 and sig3 collide and the
intermediate/least axes become indistinguishable). With the *current* floor
code, R = 1 gives `one_minus_R = EPS` → `D_applied = (0.95*Sv_eff)/1e-9` ≈ 1e9,
i.e. sig1 explodes. YHSM R maxes at 0.989 so it is not hit by this dataset,
but a different CSM export (or the col-13/eig fallback) can produce R = 1.

**Trigger:** A CSM row with R (col 13) = 1.0, or R_eig = 1.0 (w1 == w2).

**Actual behavior:** `1/(1-R)` overflow → absurd sig1 (current code), or
sig2==sig3 degeneracy (either code).

**Expected behavior:** Clip R into a safe open interval so the closures stay
well-conditioned.

**Suggested fix:**
```diff
-        R_point = np.clip(np.where(np.isfinite(R_csv), R_csv, R_eig), 0.0, 1.0)
+        # clip into the OPEN interval: R=1 collapses sig2==sig3 and blows up
+        # the differential closures (1/(1-R)); R=0 collapses sig1==sig2.
+        R_point = np.clip(np.where(np.isfinite(R_csv), R_csv, R_eig),
+                          1.0e-6, 1.0 - 1.0e-6)
```

**Test case:**
```python
import numpy as np
from project_csm_stress_to_vtu import magnitudes_C2

def test_R002_R_near_one_does_not_blow_up():
    Sv_eff = np.array([100.0]); R = np.array([1.0 - 1e-6]); D = np.array([20.0])
    sig1, sig2, sig3, D_applied, capped = magnitudes_C2(Sv_eff, R, D, 5.83)
    assert np.all(np.isfinite([sig1[0], sig2[0], sig3[0], D_applied[0]]))
    assert sig1[0] < 1.0e4   # not an overflow
    assert sig1[0] + 1e-6 >= sig2[0] >= sig3[0] > 0.0
```

---

### [R-003] [LOW] [project_csm_stress_to_vtu.py:build_tensor_andersonian] — NaN SHmax silently produces NaN tensors (no guard)

**Category:** ASSUMPTION

**Description:**
`build_tensor_andersonian(shmax_facet, ...)` takes `np.radians(az_deg)` with no
finiteness check. If a sampled CSM SHmax is NaN, the whole facet tensor is NaN,
which then propagates NaN into tractions/mu for that facet with no warning. The
YHSM csv has 0 NaN SHmax (22131/22131 finite), so the produced artifact is
unaffected — but the `--axes andersonian` path is latently fragile for any CSM
export that has NaN SHmax (the csv schema explicitly allows it: col 11 is the
SHmax *uncertainty* and is empty/NaN for YHSM). The `csm` default path uses
eigenvectors and is not exposed to this.

**Trigger:** `--axes andersonian` with a CSM csv that has any NaN in column 10.

**Actual behavior:** NaN tensor → NaN traction/mu on those facets, no warning.

**Expected behavior:** Count and warn (consistent with the degenerate-basis and
far-nn warnings already in the tool).

**Suggested fix:** after building `sigma_facet` in the andersonian branch:
```diff
         else:
             sigma_facet = build_tensor_andersonian(
                 shmax_facet, sig1, sig2, sig3)
+            n_nan_az = int(np.sum(~np.isfinite(shmax_facet)))
+            if n_nan_az:
+                print(f"warning: {n_nan_az} facets have NaN CSM SHmax; "
+                      "andersonian tensor is NaN there", file=sys.stderr)
```

**Test case:**
```python
import numpy as np
from project_csm_stress_to_vtu import build_tensor_andersonian
def test_R003_andersonian_nan_shmax_is_nan_tensor():
    az = np.array([np.nan]); s = np.array([3.0]); m = np.array([2.0]); l = np.array([1.0])
    T = build_tensor_andersonian(az, s, m, l)
    assert np.isnan(T).all()   # documents the gap the warning must flag
```

---

### [R-004] [LOW] [POSSIBLE] [project_csm_stress_to_vtu.py:muscal_sv_total_profile] — surface guard does not catch above-sea-level (topography) grids

**Category:** ASSUMPTION

**Description:**
The guard `if depth[0] > 1.0e-6: sys.exit(...)` assumes the top of the z grid is
exactly at (or below) sea level. `depth = -zt` with `zt` descending, so for a
grid whose maximum elevation is positive (real topography, e.g. +500 m),
`depth[0] = -500` and the guard does NOT trigger; the lithostatic integral then
starts above sea level and `sv_total_at` clips facet depths to a negative
`depth_grid[0]`. The production `safs_material_cvm.nc` has z_max = 0 exactly, so
this is not triggered now, but the guard is the wrong test for the stated intent
("MUSCAL z grid does not reach the surface").

**Trigger:** A material `.nc` whose z grid includes positive elevations.

**Actual behavior:** Integral starts above sea level; no error.

**Expected behavior:** Either accept topography explicitly (clamp the surface to
z=0) or fail if the grid top is not within tolerance of the surface.

**Suggested fix:**
```diff
-    if depth[0] > 1.0e-6:
-        sys.exit("ERROR: MUSCAL z grid does not reach the surface (z=0)")
+    if abs(depth[0]) > 1.0e-6:
+        sys.exit(f"ERROR: MUSCAL z grid top is at depth {depth[0]:.1f} m, "
+                 "not the surface (z=0); topography handling is not "
+                 "implemented")
```

**Test case:**
```python
# Construct a 2-level profile starting at +500 m elevation and assert the
# loader rejects it (or, post-fix, that integration starts at z=0).
# (Exercises the guard branch; no nc file required if muscal_sv_total_profile
#  is refactored to accept arrays — otherwise verify via a synthetic .nc.)
```

---

## Summary
- Critical issues: 1 (R-001 — the near-surface mu→1.0 the user asked about)
- Moderate issues: 1 (R-002)
- Low issues: 2 (R-003, R-004)
- Plan compliance: FULL — every plan §6 requirement is implemented and all
  eight acceptance criteria pass (R recovery 1.4e-13, ordering/positivity True,
  SHmax 0.069°, Luttrell self-consistency, hypocenter report, orientation
  fields byte-identical to the constant artifact). The cap behaviour (R-001) is
  a consequence of the plan being *underspecified* on the cap ("cap D so
  sig3 >= sig3_min > 0"); the implementer chose a fraction-of-Sv_eff floor that
  is dimensionally valid but physically wrong near the free surface.
- Verdict: **PASS WITH FIXES** — the artifact is correct over the seismogenic
  depth range (where rupture nucleates and propagates), but R-001 must be fixed
  before the shallow (<~1 km) stresses are trusted in a dynamic-rupture run, or
  the shallow band should be excluded/tapered.

## Direct answer: why mu ≈ 1.0 near the free surface
It is the positivity cap (R-001), not a sign/coordinate/units error. The
Luttrell differential stress (~10–60 MPa, a seismogenic-depth quantity) is far
larger than the shallow effective overburden Sv_eff (≈0.14 MPa at 21 m).
`magnitudes_C2` caps by setting sig3 = 0.05·Sv_eff, which pins the stress ratio
at sig1/sig3 ≈ 31 (independent of how shallow the facet is) — vastly
super-critical relative to Byerlee (k = 3.12–5.83). A ratio that high gives
apparent friction up to ~1.0 on the fault. Capping by a frictional ceiling
(D_max = (k_max−1)·Sv_eff/((1−R)k_max+R)) instead makes the differential vanish
as Sv_eff → 0 and bounds mu to the chosen Byerlee value.

## Unreviewed Areas
- The constant-magnitude path: verified byte-identical to the prior artifact
  (VTU `cmp`-identical; summary differs only in the renamed input-csv path), so
  not re-audited for logic.
- The custom VTU binary writer (`write_vtu`): pre-existing code, unchanged by
  this work; round-tripped successfully when re-read for this review.
- `apply_deep_taper`: implemented per plan but OFF by default and not used in
  the deliverable; smoke-tested (`--taper-zseis-km 15` runs, ordering holds) but
  its isotropic-continuation choice (built-tensor mean) is not validated against
  an independent reference.
