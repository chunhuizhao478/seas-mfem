# CSM Static Equilibration — Implementation Plan

**Scope:** Initialize the bulk stress field for SAFS dynamic / quasi-dynamic
runs from a Community Stress Model (CSM), in a way that leaves bulk Q as a
clean reading of stress oscillations about the background, and avoids
slow drift driven by `∇·σ_bg ≠ 0` in the discrete operator.

**Status:** Plan (not implemented). Target branch: `feature/safs-quasi-dynamic`.

**Why we need this:** TPV102's fluctuation-Q trick (`Q(0) = 0`, background
in `DOFData` pre-stress fields) works only because the TPV102 background is
spatially uniform, so `∇·σ_bg ≡ 0`. For a heterogeneous CSM-derived
background, the divergence is generically nonzero and acts as a persistent
body-force source in the momentum equation. The fix is a one-shot
elastostatic *equilibration* precompute that projects the CSM into the
discrete kernel of the divergence operator under the production BCs.

---

## 0. Inputs and assumptions

- **Mesh.** SAFS tetrahedral mesh (currently from the archived
  `project_7.0_preferred` CGAL-corefine pipeline). DG polynomial degree
  `p` is a parameter; default `p = 2`.
- **Material.** Heterogeneous `ρ(x), λ(x), μ(x)` (CVM-derived) or
  homogeneous half-space depending on production config; treat as a
  `Coefficient` evaluated at quadrature points throughout.
- **CSM source.** SCEC Community Stress Model — most likely the
  Hauksson/Yang 2014 stress inversion or the Hardebeck-style model.
  Format: regular lat/lon/depth grid of six stress components (or
  principal stresses + orientations) in SI units.
- **Gravity.** Carry analytically; do not project lithostatic part into
  the static solve (Phase 2).
- **Production BCs (assumed quasi-dynamic SEAS-style):**
  - Sides (attrs 1–4): Dirichlet plate loading `u_X = sgn(Y)·V_p·t/2`
  - Bottom: Dirichlet plate loading
  - Top: free surface, `σ·n = 0`
  - Fault: rate-and-state friction (locked during equilibration)

If the production BCs change to absorbing far-field (dynamic-rupture
TPV102-style), see Phase 3 §Alt-BC for the modified BC table.

---

## 1. Phase 1 — Sample the CSM into a `Coefficient`

### Deliverable
A `VectorCoefficient` (six components) `sigma_csm(x)` that, when evaluated
at any point `x` in mesh coordinates, returns the CSM stress tensor in
global Cartesian (SI units).

### Steps
1. **Coordinate transform.** Implement `Mesh → CSM` projector. Mesh is in
   UTM Zone 11 meters; CSM is lat/lon/depth. Use a stable projection
   library (PROJ or equivalent) — do NOT roll your own ellipsoidal
   transform.
2. **Trilinear interpolation** at the *quadrature points* used by the
   downstream L2 projection (Phase 2), not at element centroids — for
   `p ≥ 2` centroid sampling loses intra-element variation.
3. **Pre-smooth.** Convolve the CSM grid with a Gaussian kernel of
   width `σ_smooth ≈ 1.5 × h_mesh` (mesh element size) before
   sampling. Document the smoothing — it's the single biggest knob
   for whether the projected field is benign or ringing.
4. **Out-of-CSM regions.** SAFS mesh extends below or beyond the CSM
   support. Decide policy: clamp to nearest in-support value, or
   blend to lithostatic. Default: clamp + log a warning with the
   fraction of quadrature points using the fallback.

### Test (Phase 1)
- Unit test: a synthetic CSM grid with `σ_xx = 10 MPa, σ_yy = -20 MPa,
  σ_xy = 5 MPa` constant. Verify `sigma_csm(x)` returns those values
  at a sample of interior points to machine precision (modulo the
  Gaussian smoothing edge effects).
- Sanity plot: visualize each of the six tensor components on a vertical
  cross-section through the mesh; compare against the SCEC CSM plotting
  utility on the same section.

---

## 2. Phase 2 — Lithostatic / tectonic split

### Deliverable
Two `VectorCoefficient`s:
- `sigma_litho(x)` — analytical lithostatic, depth-only
- `sigma_tect(x)` — `sigma_csm(x) - sigma_litho(x)`

### Why split

- Lithostatic σ_zz reaches ~600 MPa at 20 km depth; tectonic anomalies are
  10–100 MPa. The static solve only needs to chase the tectonic part.
- `∇·σ_litho` is exactly cancelled by gravity (`∂_z σ_litho = -ρg`), so we
  never need to project it into the discrete divergence operator. Carry it
  as a body force from t=0 onward.

### Formula

```
σ_litho(z) = -ρ(z) · g · z · I       (compression-positive convention)
```

For depth-dependent `ρ(z)` from CVM, integrate numerically:
`σ_litho(z) = -g · ∫_0^z ρ(ζ) dζ · I`.

### Test (Phase 2)
- Verify `∇·σ_litho + ρ(z)·g·ê_z = 0` (analytically) on a 1D depth
  profile. This is a hand check, not a unit test.
- Visualize `σ_tect = σ_csm - σ_litho` cross-section: should show
  smooth tectonic anomalies on top of a small residual (mean ~10
  MPa scale).

---

## 3. Phase 3 — Static elastostatic equilibration

### Deliverable
A `GridFunction sigma_bg_eq` (DG stress space, six components) such that
`∇·σ_bg_eq` is in the kernel of the discrete divergence operator under
the production BC set evaluated at t=0.

### Math — reference-stress formulation (after PyLith)

PyLith's governing equation for static elasticity with body forces and a
reference state is

```
∇·(σ_ref + C:ε(u_eq)) + ρ·g = 0
```

with the constitutive law evaluating *total* stress as
`σ_total = σ_ref + C:ε(u_eq)`. If `σ_ref` already balances `ρ·g` under
the chosen BCs (i.e., `∇·σ_ref + ρ·g = 0` weakly), the solver converges
with `u_eq ≡ 0` and `σ_total ≡ σ_ref` — no spurious deformation, no
boundary-layer artifacts.

Map to our notation: identify `σ_ref ↔ σ_litho + σ_tect` (lithostatic +
CSM tectonic) and solve for the residual displacement `u_eq` that
absorbs whatever the CSM does *not* equilibrate analytically. Weak
form, for all `w ∈ V_h` with `w·n = 0` on roller boundaries (only the
normal-component test functions vanish):

```
∫_Ω C(x):∇u_eq : ∇w  =  -∫_Ω (σ_tect):∇w
                       +  ∫_∂Ω_top (-σ_tect·n)·w           (free-surface correction)
                       +  ∫_∂Ω_side+bot (σ_tect·n)·w_tang   (tangential traction balance)
```

where `w_tang` denotes the tangential-component test functions on
roller boundaries (the normal-component contribution is zero because
`w·n = 0` there by construction). The third term is what the full
Dirichlet version was silently throwing away — it is the missing
boundary-layer-killer.

Then define `σ_bg_eq = σ_litho + σ_tect + C:∇u_eq`. If the CSM is
already close to elastic equilibrium with the mesh + BCs, `u_eq` is
small and `σ_bg_eq ≈ σ_csm` to a tight tolerance.

### Boundary conditions — roller pattern (after PyLith)

Adopt PyLith's static-initialization pattern: **roller BCs** (pin only the
normal displacement component) on lateral and bottom faces, plus the
Neumann surface-traction adjustment on top. A full Dirichlet `u_eq = 0`
on sides + bottom would over-constrain the tangential components and
create a shear-traction boundary layer at depth — visible as a drift in
σ_xz, σ_yz near the bottom and sides after equilibration. Roller still
kills all 6 rigid-body modes (3 translations from the 3 mutually
perpendicular normal-pin sets, 3 rotations from non-collinearity of the
faces), so the static system is well-posed without extra anchoring.

| Boundary | Production dynamic BC | Static equilibration BC |
|---|---|---|
| Sides ±x (attrs 1, 2) | Dirichlet `u_X = sgn(Y)·V_p·t/2` | Roller: `u_x = 0`, `u_y, u_z` free |
| Sides ±y (attrs 3, 4) | Dirichlet plate loading | Roller: `u_y = 0`, `u_x, u_z` free |
| Bottom −z (attr 6) | Dirichlet plate loading | Roller: `u_z = 0`, `u_x, u_y` free |
| Top +z (attr 5, free surface) | Neumann `σ·n = 0` | Neumann `C:∇u_eq·n = -σ_tect·n` (cancels CSM surface traction so total = 0) |
| Fault interface | rate-and-state friction | Locked: enforce zero slip jump |

**Rule:** static BCs = production BCs *projected to their normal
components* with the surface-traction adjustment on top, fault locked.
Tangential components are left free on lateral + bottom faces so the
solution can relax to a near-zero displacement everywhere — matching
PyLith's reverse-2d step02 result that `u_eq ≈ 0` after the first
nonlinear residual when `σ_ref` already balances gravity (here
generalized to CSM-derived `σ_tect`).

**Why this is safe even though the dynamic BC is full Dirichlet.** At
`t = 0+` the production plate-loading Dirichlet evaluates to
`u = sgn(Y)·V_p·0 = 0` on all three components, which is *consistent*
with roller (roller pins one component to 0; the dynamic BC pins all
three to 0 — both agree on the normal-component value of 0). The
*operator kernel* changes between static and dynamic, but `σ_bg_eq` is
in the static kernel by construction, and the only thing the dynamic
operator adds at `t = 0+` is the tangential-component pin, which
produces no immediate stress because `u = 0` already satisfies it.
A small `O(V_p · dt_0)` transient appears on the first time step from
the plate loading; that is *physical*, not a discretization artifact.

**Implementation note (MFEM).** Roller = component-wise Dirichlet.
In `mfem::ParGridFunction::ProjectBdrCoefficient` terms, this is a
`VectorArrayCoefficient` with the active component set to zero and
the other two components attached to a `mfem::ConstantCoefficient`
that the solver does not constrain (or, more cleanly, a per-attribute
`Array<int>` of essential DOFs computed from `GetEssentialVDofs` with
the normal-component bit set per boundary attribute). The MFEM
`elasticity_operator.hpp` already supports per-attribute essential BC
arrays via `BoundaryConfig`; the new code path is an extra "roller"
`BCMode` enumerant that builds the per-attribute, per-component
essential-DOF set.

### Fault as locked interface

Same DG operator as the dynamic run, but with the slip jump pinned to
zero (no friction coupling, no `tau_pre`-trial recombination). Concretely:
on fault faces, replace the dynamic friction-flux dispatcher with an
interior-face DG flux that treats the fault like any other element
interface. The result is `[[u_eq]] = 0` to discretization accuracy on
the fault, which is what "locked" means.

### Solver

- AMG-preconditioned CG (MFEM `HypreBoomerAMG` + `CGSolver`).
- Tolerance: `rtol = 1e-10` (this is a precompute, accuracy matters).
- Estimated cost: minutes to tens of minutes on the SAFS mesh; one-shot,
  not in the inner loop.

### Alt-BC: dynamic-rupture / absorbing far-field

If production swaps to TPV102-style absorbing BCs on sides + bottom,
absorbing dampers are a *dynamic-only* concept and have no analog at
`t = 0`. Use roller BCs for the static solve regardless of what the
dynamic operator does at the far-field — same table as above. Roller
already pins all 6 rigid-body modes, so no extra anchoring is needed.

(Earlier draft proposed full Dirichlet on the bottom plus 3-2-1 pinning
as alternatives. Both are superseded by roller, which is strictly less
constraining tangentially and is what PyLith does in production.)

### Test (Phase 3)
- **Discrete divergence check.** After solving, evaluate
  `‖∇·σ_bg_eq‖_L2(Ω) / ‖σ_bg_eq‖_L2(Ω)`. Should be `< 1e-8`. If not,
  the static problem isn't actually solved to tolerance, or the BCs
  are inconsistent with the divergence operator.
- **Free-surface check.** Evaluate `‖σ_bg_eq · n‖_L2(top)` and compare
  to `‖σ_csm · n‖_L2(top)`. Should drop by ≥3 orders of magnitude.
- **No-deep-pollution check (the test that motivates roller).** Sample
  `σ_bg_eq` along a vertical line `x = 0, y = 0, z ∈ [-Z_bot, 0]` from
  surface to bottom. Compute the relative deviation
  `‖σ_bg_eq(z) − σ_csm(z)‖ / ‖σ_csm(z)‖` as a function of `z`. Two
  acceptance criteria:
  - The deviation must not exceed `1e-2` in the deepest 10% of the
    domain (the would-be boundary-layer zone).
  - The deviation must be monotonically *decreasing* with distance
    from the bottom for at least the bottom 30% of the domain.
  Failing either signals that the BCs are still over-constraining
  tangential displacement at depth. Replace bottom Dirichlet with
  roller and re-run.
- **Reference-stress consistency check (after PyLith step02).** With
  a *uniform* synthetic σ_tect that exactly balances ρg under roller
  BCs (e.g., σ_tect = 0, σ_litho = -ρgz·I), confirm that the solver
  converges in one CG iteration with `‖u_eq‖_L2 < 1e-12 · L · ε_machine`.
  This is the analog of PyLith's "first nonlinear residual meets
  convergence criteria" verification.

---

## 4. Phase 4 — Fault interface `tau_pre` extraction

### Deliverable
Per-fault-DOF `tau_pre` populated from `σ_bg_eq` in the fault-local
`(t1, t2, n)` frame:

```
tau_pre[i].sigma_n  =  n_i · σ_bg_eq(x_i) · n_i      (compression-positive)
tau_pre[i].tau_1    =  t1_i · σ_bg_eq(x_i) · n_i
tau_pre[i].tau_2    =  t2_i · σ_bg_eq(x_i) · n_i
```

This is the heterogeneous generalization of BP5's analytic `tau0_vec`
(`config/bp5_params.hpp:333–381`) and TPV102's uniform
`DOFData.sigma_n0 / tau1_0 / tau2_0` (`dynamic/tpv102_setup.hpp:76–78`).

### Frame convention
Use `FaultBasis` (Tandem convention, project-wide): `t1 = dip,
t2 = strike, n = ref_normal`. This is enforced project-wide per CLAUDE.md
"R-801 Option A". Do NOT introduce a competing fault-local frame here.

### Wiring
- Add to `FaultGeometry`: a method `LoadTauPreFromBackground(const GridFunction &sigma_bg_eq)` that evaluates the
  GridFunction at each fault DOF coordinate and projects.
- BP5 driver path keeps its existing analytic `tau0_vec` call.
- SAFS driver calls the new method instead.
- Both populate the same `tau_pre_` member, so downstream
  `RateAndStateFault::ComputeRHS` is unchanged.

### Test (Phase 4)
- Unit test on a planar fault at `y = 0` with synthetic `σ_bg_eq`
  uniform in space: `σ_xy = -75 MPa, σ_yy = -120 MPa` (TPV102 values).
  Verify `tau_pre[i].sigma_n = 120 MPa, tau_pre[i].tau_2 = 75 MPa,
  tau_pre[i].tau_1 = 0` to machine precision at every fault DOF.
- Cross-check against TPV102's `InitializeFaultDOFs` output on the
  same fault and CSM-as-uniform field — should match bit-for-bit.

---

## 5. Phase 5 — End-to-end verification (the test that catches everything)

### The test
1. Load mesh, load CSM, run Phases 1–4 to get `sigma_bg_eq` and `tau_pre`.
2. Run the dynamic / quasi-dynamic operator with:
   - `Q(0) = 0` (fluctuation-Q) or `u(0) = 0` (quasi-dynamic)
   - `tau_pre` from Phase 4
   - Fault locked (friction disabled, slip pinned at 0)
   - No nucleation, no plate loading (`V_p = 0` for this test only)
   - Production BCs otherwise
3. Advance to `t = 100 s` (dynamic) or `t = 10 yr` (quasi-dynamic).
4. Measure `‖Q(t)‖_L2(Ω) / ‖σ_bg_eq‖_L2(Ω)`.

### Pass criterion
`< 1e-10` for dynamic, `< 1e-8` for quasi-dynamic (loose because the
elastostatic solve in the time loop has its own tolerance).

### What failure means
- Uniform drift across the domain → BC mismatch, almost always the
  free surface (Phase 3 Neumann adjustment missing or wrong sign).
- Boundary-layer near sides/bottom → static Dirichlet not matching
  dynamic Dirichlet at t=0.
- Fault-localized drift → "locked" condition not strict (friction
  not disabled, or slip pinned with wrong sign).
- Whole-domain ringing → CSM not smoothed enough relative to mesh,
  high-frequency content beyond what DG can represent.

### Why this test is load-bearing
It validates the *entire* pipeline (sampling + projection + smoothing +
litho-split + equilibration + BCs + frame projection) end-to-end with
a single number. Pass means the dynamics afterward is the same physics
already validated on TPV102 and BP5. Skip this and you'll be chasing
drift in a production run, where it will be entangled with rupture
physics.

---

## 6. Phase 6 — Integration into the SAFS driver

### Deliverable
A configuration flag `safs.stress_init = csm | uniform` in the SAFS TOML
config. When `csm`, the driver:
1. Loads CSM file (`safs.csm_file = ...`).
2. Runs Phases 1–4 at startup.
3. Caches `sigma_bg_eq` (DG `GridFunction`) and `tau_pre` to disk.
4. On restart with the same mesh + CSM file + BC config, reads the
   cached precompute and skips Phases 1–3.

Cache invalidation key: SHA-256 of `(mesh_file, csm_file, p, BC_config,
phase_1_smoothing_params)`.

### Test (Phase 6)
- Integration test: SAFS smoke run with a small CSM region. Verify
  precompute + restart consistency (cached run produces identical
  `sigma_bg_eq` to fresh run).

---

## Open questions / decisions deferred

1. **CSM choice.** Hauksson/Yang vs Hardebeck vs newer SCEC model?
   Pick before Phase 1; affects file format and loader.
2. **Topography.** Does the SAFS mesh include real topography, or
   flat top at z=0? If real, the free-surface Neumann adjustment
   has to be evaluated on a curved surface — works with the proposed
   formulation but worth verifying numerically.
3. **CVM coupling.** Heterogeneous `(λ, μ)` from CVM is already on the
   roadmap; the static equilibration interacts with material
   heterogeneity through `C(x)`. Phase 3 is written to support it
   (use `Coefficient` for `λ, μ`); confirm there's no hidden
   homogeneous assumption in the existing elasticity operator
   (`domain/elasticity_operator.hpp`).
4. **Quasi-dynamic vs dynamic.** SAFS branch is currently
   quasi-dynamic. If we eventually run dynamic rupture from the
   same CSM, Phase 3 stays the same but Phase 6 needs a parallel
   wiring into `tpv102_driver`-style code. Defer until needed.
5. **Equilibration tolerance vs production tolerance.** Production
   elastostatic solve in the quasi-dynamic time loop uses some
   tolerance `rtol_prod`. If Phase 3 uses `rtol = 1e-10` and the
   time loop uses `rtol = 1e-6`, the *first* time step will see
   a `1e-6`-scale residual that wasn't there in the static check.
   Verify this is within the Phase 5 pass criterion. If not, tighten
   `rtol_prod` for the first few steps.

---

## File-level changes (anticipated)

| Path | Change |
|---|---|
| `safs/csm/csm_loader.hpp` | new — CSM file reader + `VectorCoefficient` |
| `safs/csm/csm_smoothing.hpp` | new — Gaussian pre-smoother |
| `safs/csm/lithostatic.hpp` | new — analytical depth profile |
| `safs/equilibration/static_solver.hpp` | new — Phase 3 elastostatic solver |
| `fault/fault_geometry.hpp` | add `LoadTauPreFromBackground()` |
| `solver/seas_operator.hpp` | branch on `stress_init` flag |
| `config/safs_config.hpp` | add `csm_file`, `stress_init`, smoothing params |
| `tests/verification/test_csm_equilibration.cpp` | new — Phase 5 end-to-end test |

No changes to: BP5 driver, TPV102 driver, TPV104/205 paths, friction
solver, time stepper. The CSM pipeline is *additive*; existing
benchmark paths are unaffected.

---

## References — prior art

The static-equilibration design above is modeled directly on PyLith's
established workflow for gravitational body forces with a reference
stress. Two PyLith examples are the load-bearing precedents:

- **Reverse-2D Step 1 (gravity alone, no reference stress).** Applies
  ρg to an undeformed, stress-free domain with roller BCs on the
  lateral and bottom faces and a free top. Produces ~2 km of spurious
  vertical deformation, demonstrating *why* a reference stress is
  needed. → `https://pylith.readthedocs.io/en/stable/user/examples/reverse-2d/step01-gravity.html`

- **Reverse-2D Step 2 (gravity + reference state).** Adds a reference
  stress `σ_ref = -ρgz·I` (depth-dependent lithostatic) and shows that
  "the first nonlinear solver residual evaluation meets the convergence
  criteria" with zero displacement — i.e., `u_eq ≡ 0` when `σ_ref`
  balances `ρg`. The Phase-3 reference-stress-consistency test above
  is the analog of this verification. → `https://pylith.readthedocs.io/en/stable/user/examples/reverse-2d/step02-gravity-refstate.html`

- **Subduction-3D Step 8 (3D gravity + roller).** Confirms the same
  pattern scales to 3D and to heterogeneous materials. → `https://github.com/geodynamics/pylith/blob/main/docs/user/examples/subduction-3d/step08-gravity.md`

Adaptations for CSM-driven SAFS:

1. PyLith's `σ_ref` is a per-material scalar field provided via a
   spatial database; our `σ_tect` is a per-mesh-point tensor field
   sampled from the SCEC CSM. The math is identical — `σ_ref` and
   `σ_tect` enter the constitutive evaluation in the same place
   (`σ_total = σ_ref + C:ε(u)`).
2. PyLith's spatial-database trilinear lookup is replaced by the
   Phase-1 CSM loader, with the same out-of-bbox + smoothing policy.
3. PyLith's roller BCs are component-wise Dirichlet via
   `DirichletTimeDependent` with `constrained_dof = [<normal>]`.
   Our MFEM equivalent is a new `BCMode::Roller` enumerant in
   `domain/boundary_config.hpp` that builds the per-attribute
   essential-DOF set with only the normal-component bit set.

What PyLith does *not* provide and we have to add ourselves:

- A heterogeneous tensorial reference stress with full off-diagonal
  components (PyLith's reverse-2d/step02 uses isotropic `-ρgz·I`).
  This is straightforward — the constitutive contract is unchanged;
  we just need the loader.
- An end-to-end equilibration verification test (Phase 5 above) that
  closes the loop with the dynamic operator. PyLith's example stops
  at the static residual; we need to confirm that `σ_bg_eq` also
  sits in the *dynamic* operator's kernel under the locked-fault,
  zero-loading configuration.
