# TPV104 σ_n leak — root cause: mesh non-mirror-symmetry

**Date:** 2026-04-24
**Status:** confirmed root cause for the σ_n perturbation channel; secondary candidate for V_strike over-shoot

## Background

The t12 dev production run (`results_t12_dev_job7677022`) and the
dt-halving diagnostic run (`results_dthalf_job7677079`) both showed
σ_n perturbations of 1-4 MPa at every fault station (SeisSol reference
shows 0 perturbation, σ_n flat at 120 MPa).  Mirror-paired stations
across `x2 = 0` showed sign-flipped σ_n peaks of comparable magnitude
(e.g., `x2_-12_x3_3` → +3.7 MPa, `x2_+12_x3_3` → -4.2 MPa).

## Eliminations

| diagnostic | result | rules out |
|---|---|---|
| dt-halving (5.7× finer dt) | σ_n leak grew 0.48 → 0.84 MPa | macro-step ψ cadence |
| D2/D3 unit tests on `ComputeTrialTraction` | 20/20 PASS at 1e-12 | `FaultFaceFlux` arithmetic |
| D1 instrumentation in `wave_operator.inl` | bimodal sign_flipped (52/48), but `can_n` uniform across all 4938 QPs | per-QP rotation pipeline |

## Root cause: mesh asymmetry

The .geo source (`tpv104/mesh/tpv104_1000m.geo`) builds a
**bilaterally-symmetric geometry** (domain `(-60, 60) × (-60, 60) ×
(-60, 0)` km, fault at `y = 0`, all features mirror across `x = 0`),
but Gmsh's unstructured Delaunay tetrahedralisation **does not
preserve the symmetry**.  Verified by `tpv104/scripts/tpv104_mesh_mirror_check.py`
on `tpv104_1000m.msh`:

```
=== Mirror check across x (along-strike) = 0 ===
  Nodes matching mirror partner: 652 / 37146 (1.76%)
  Elements matching mirror tet: 0 / 215143 (0.00%)

=== Mirror check across y (fault-normal) = 0 ===
  Nodes matching mirror partner: 1098 / 37146 (2.96%)
  Elements matching mirror tet: 0 / 215143 (0.00%)

=== Fault-surface (Physical Surface 3) mirror check across x = 0 ===
  Fault triangles: 1646
  Fault triangles with mirror partner across x = 0: 479 / 1646 (29.10%)
  Fault NODES with mirror partner across x = 0: 427 / 878 (48.63%)
```

**Zero of 215,143 tetrahedra** have a mirror-paired tet across either
axis.  Even on the fault surface itself, only 29% of triangles and
49% of nodes have mirror partners.

Implication: a discretized solution on this mesh **cannot** preserve
mirror symmetry exactly.  Any small fault-normal radiation that should
cancel by symmetry on a perfectly-mirror-symmetric mesh leaves an
O(h^{p+1}) residue at p=1 that propagates into σ_n at every fault QP.
For `h ≈ 200 m` and the rupture-front velocity-jump magnitudes we see
(~10 m/s peak), the predicted residue scale is consistent with the
0.5-4 MPa σ_n peaks observed in the t12 run.

## Why dt-halving made this worse

Mesh asymmetry is a **spatial** structural error, not a temporal one.
At finer dt, the wave operator resolves higher-frequency modes more
sharply.  For a mirror-symmetric mesh those modes preserve symmetry;
for an asymmetric mesh they propagate the asymmetric residue with
greater fidelity, so finer dt **amplifies** the σ_n leak instead of
shrinking it.  Observed at the hypocenter: σ_n peak grew from 0.48 to
0.84 MPa under 5.7× finer dt — consistent with this mechanism.

## Why per-QP `sign_flipped` looked like a smoking gun but isn't

`sign_flipped` is bimodal (52% / 48%) across the fault, and it IS
correlated with mesh-element orientation — different on the +y vs -y
side of each fault face depending on Gmsh's enumeration.  But the
wave operator's negation block at `wave_operator.inl:2117-2123`
correctly produces a uniform canonical normal `can_n = (0, -1, 0)`
for every QP, and the routing `elem1_on_plus = !sign_flipped` puts
each element's Q on the right side of the trial traction.  D1 verified
this directly.  The bimodal `sign_flipped` is a CONSEQUENCE of mesh
asymmetry, not a bug in itself.

## Why this matters for the V_strike over-shoot

The 1.7-2.7× peak V_strike over-shoot at the propagating-front
stations is plausibly **also** driven by mesh asymmetry (rupture
fronts on asymmetric tets accelerate non-uniformly across mirror-paired
stations), but the connection is indirect.  Need a controlled
experiment with a symmetric mesh to confirm.  Until then, the V_strike
finding is "consistent with mesh asymmetry but not yet attributed."

## Fix paths

| option | mechanism | cost | symmetry guarantee |
|---|---|---|---|
| **A. Octant-and-reflect** mesh generation | build mesh on `x ≥ 0, y ≥ 0, z ≤ 0` octant, reflect 4× to fill domain | 1-2 d work; needs Python+gmsh | exact |
| **B. Transfinite + structured tet split** in Gmsh | use `Transfinite Surface` / `Transfinite Volume` directives + Kuhn split (`Recombine Mesh = Tetrahedra`) | 0.5-1 d; non-trivial because of fault embedding | high; mesh is symmetric by construction |
| **C. Higher polynomial order** (p = 2 or 3) | reduces O(h^{p+1}) discretization residue | low (CLI flag); but increases compute cost ~5× per order | none — just smaller residue |
| **D. Run-time mesh symmetrisation** at init | iterate mesh nodes; for each, snap to mirror-paired position if a partner exists within tolerance | 1-2 h; fragile (only fixes nodes that already are nearly-mirror-paired) | partial; doesn't fix the tet topology |

Option **A** is the right long-term fix.  Option **C** at p = 2 or 3
is a viable short-term mitigation if compute budget allows.

## Next diagnostic experiment

To **confirm** the mesh-asymmetry → σ_n-leak attribution:

1. Generate a mirror-symmetric mesh via option (A) or (D) at the same
   nominal resolution as `tpv104_1000m.msh`.
2. Run the same dt-halving sbatch (already committed) on the symmetric
   mesh.
3. Compare the σ_n peak magnitude at the hypocenter.

Acceptance: σ_n peak should drop from ~0.5 MPa to ≤ 0.05 MPa (10×
reduction) under a fully-symmetric mesh.

## Tooling

* `tpv104/scripts/tpv104_mesh_mirror_check.py` — runs the mirror
  symmetry check on any `.msh` file.  Default: `tpv104_1000m.msh`.
* `wave_operator.inl` D1 instrumentation block (gated by
  `-DSEAS_DIAG_TPV104_FAULT_BASIS`) — prints per-QP `sign_flipped` +
  basis vectors, useful for verifying that the rotation pipeline
  remains consistent across mesh changes.
* `tests/unit/test_tpv104_sigma_n_invariance.cpp` — D2/D3 regression
  test: `make seas_test_tpv104_sigma_n_invariance`.
