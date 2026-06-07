# Dip-drift — LOCAL reproduction + hypothesis test (2026-06-06)

Date: 2026-06-06
Author: Claude (code-debug, goal: "reproduce locally + propose fixes")
Worktree: `worktree-fix-dip-drift` @ `b09dbfc`. Build: conda `mfem-dev`, MFEM_* overrides.
Companion: `dip_drift_debug_plan_2026-06-06.md` (this doc CORRECTS its H1-prime ranking).

> All experiments below run LOCALLY on tiny fixtures (no Frontera, no production
> mesh), per project constraints. No no-touch file edited; investigation code is
> in `tests/unit/` and the fix levers are in editable `dynamic/`.

---

## Result 1 (DECISIVE) — H1 (frame) is REFUTED as the dip seed, in compiled code

New compiled probe `tests/unit/test_frame_dip_leak_probe.cpp` (links the REAL
`GodunovFlux::BuildRotationInverse`, `make seas_test_frame_dip_leak_probe`). It
rotates a pure strike-slip global stress (`σ_xy=τ`) through the fault frame under
deliberate frame errors and reads the dip channel `local SXY = n·σ·t1`:

| frame condition | dip channel | note |
|---|---|---|
| ideal | `0.0` | control |
| flip sign of n / t1 / t2 | **`0.0` (all)** | sign flip does NOT leak |
| smooth normal tilt (FaultBasis-rebuilt), ε≤0.1 | **`0.0`** | construction self-corrects |
| in-plane t1↔t2 rotation by ε | `−ε·τ` | leaks — but FaultBasis pins strike to `up`, so this never occurs |
| ideal frame + global `σ_yz = 1e-3·τ` | **`σ_yz` (1:1)** | the ONLY realizable leak path |

**Conclusion:** the rotation cannot manufacture dip-channel traction from a pure
strike-shear state under any sign flip or normal tilt. The dip traction can only
come from a **genuine global `σ_yz`** (vertical/dip shear on the fault plane)
delivered by the bulk wavefield. The clean rotation then faithfully resolves it,
and friction rectifies it into dip slip. **The fault FRAME is exonerated.** This
also means the debug-plan's §6-A "freeze the per-face frame" fix would be a
**no-op** (the probe shows the frame is already clean).

Confirms the earlier dependency-free Python check (`/tmp/frame_dip_leak_probe.py`)
in the actual compiled rotation, and aligns with Agent-A's numerical tilt result.

---

## Result 2 — the [[v_dip]] SEED IS locally reproducible (~46 s with short nsteps)

`seas_test_fault_planar_serial` (structured tet slab, ADER-2; `[[v_dip]]` probe).
The earlier ~30-min cost was purely oversized `nsteps=4000`. The dip SEED peaks at
rupture-front passage (early), so a **short run reproduces it in 46 s**:

`--nx 4 --nz 4 --nsteps 150 --no-symmetry-probe` (serial p1):
```
[main] DIP leak: peak |[[v_dip]]| = 3.50e-04 m/s  (vs strike |[[v_str]]| = 4.03e-02)
```
i.e. a real dip-tangential velocity jump ~0.9% of the strike jump, growing
~linearly with the front. (The σ_n>1e5 acceptance only fails here because 150
steps is too short for the *secular* drift — the seed itself is unambiguous.)

NOTE: this reproduces the SEED only. The production *fingerprint* (anti-symmetric,
linear-in-strike, free-surface-peaked) needs the production mesh + free surface
(absent locally), so spatial structure stays a Frontera question.

## Result 3 — quick A/Bs on the seed (each ~46 s)

| run (150 steps) | peak \|[[v_dip]]\| | vs baseline |
|---|---|---|
| baseline pure upwind, k=0 | 3.50e-4 | — |
| `--fault-overint 2` | 5.64e-4 | **×1.6 WORSE** |
| `SEAS_FORCE_V1_ZERO=1` | 3.38e-4 | ×0.97 (≈unchanged) |

- **Over-integration alone WORSENS the dip seed** (confirms June §9; not the cure).
- **Breaking the dip friction reaction barely changes the seed** ⇒ `[[v_dip]]` is
  a genuine **bulk-wavefield σ_yz**, not a friction-loop artifact. The H2 closed
  loop is a long-secular amplifier, not the short-window seed driver (growth here
  is linear, not exponential).

---

## Result 4 — the ENTIRE pointwise operator stack is symmetry-clean (compiled)

Built `seas_test_offfault_frame_symmetry` (<1 s). Tests flux **equivariance** under
both fault-plane symmetries: `M·Flux(n,Q) == Flux(Mn, M·Q)`.

- **Fault-normal reflection (y→−y):** off-fault `Interior` (upwind) and `Central`
  residual = **0** at all z-tilts.
- **Strike-reversal (x→−x)** — the symmetry the dip drift actually breaks (dip is
  ANTI-symmetric in x): `Interior`, `Central`, **`FreeSurfaceGodunovTotal`, and
  `FreeSurfaceTotal`** residual = **0** at all tilts.

⇒ Every pointwise flux (fault frame, fault embedding [matches SeisSol component-
wise per the bulk-embedding agent], off-fault upwind, central, both free-surface)
is EXACTLY symmetry-equivariant. The agent's "off-fault `BuildFrame` breaks
y-symmetry" lead is **refuted** — the Gram-Schmidt tangent-sign asymmetry cancels
in the full flux. **The symmetry break that produces the dip drift is NOT in any
pointwise operator.**

### Key reframe
The dip drift breaks **strike-reversal (x→−x)**, NOT the fault-normal reflection.
Prior "symmetric mesh" attempts (incl. the user's fully-symmetric/conforming mesh,
and June's harness) were **y-mirror / fault-normal** symmetric — the WRONG symmetry.
A mesh/setup symmetric under x→−x is what would matter.

### What survives (p-independent, operator-symmetric, NOT yet tested)
Since every operator is symmetric and the effect is order-independent, the break
must enter via something non-pointwise AND resolution-independent:
1. **Mesh geometry not symmetric under strike-reversal (x→−x)** — but a pure
   mesh-dissipation asymmetry should shrink with p; the p-independence is a puzzle
   for this unless the effect is O(1).
2. **NUCLEATION source applied with an x-asymmetry** (strong, untested lead): an
   overstress patch seeded with a sub-element or stencil-asymmetric x-profile would
   inject an x-ANTISYMMETRIC, p-INDEPENDENT dip perturbation that no flux/mesh/
   dealiasing change removes. The user has NOT tried initial-equilibrium/nucleation.
3. Prestress / equilibrium setup x-asymmetry.

### Decisive next tests (cheap/local)
- Check the nucleation/overstress for exact x-symmetry: does the seeded `tau2_nuc`
  (and any per-DOF profile) satisfy `f(x)=f(−x)` to machine precision about the
  hypocenter? An x-antisymmetric component is the prime suspect.
- Build an x-reflection-symmetric slab in the harness; if `[[v_dip]]` collapses,
  mesh strike-asymmetry is implicated; if it persists, mesh is excluded too.

## Result 5 — ROOT CAUSE CONFIRMED: the "symmetric mesh" is symmetric in the WRONG axis

The `_symmetric_mesh` gold runs (jobid 7774565/7774567) use
`tpv102/mesh/tpv102_200m_hybrid.msh`. Its `.geo` (Pelties 2012 construction):
the **±y prism strips are exact mirror partners** → symmetric across the FAULT
PLANE (y), which is what tamed the σ_n leak. BUT the fault surface itself is
triangulated with `Mesh.Algorithm=6` (Frontal-Delaunay, UNSTRUCTURED), with **no
strike (x) symmetry constraint**.

Empirical check on the actual 174 MB mesh: of 18 861 fault-plane nodes, **only
1.78% have an exact strike-mirror (−x) partner** → the fault mesh is
**strike(x)-ASYMMETRIC**.

### This resolves everything
- Every pointwise operator is exactly strike-reversal (x→−x) equivariant
  (Result 4). So the assembled operator's only x-asymmetry comes from the **mesh
  geometry** baked into the assembly (element shapes, face normals, QPs).
- An x-asymmetric mesh ⇒ x-asymmetric discrete operator ⇒ x-ANTISYMMETRIC spurious
  `[[v_z]]` ⇒ friction rectifies it ⇒ the observed `dip_slip(x) ≈ θ·x` drift.
- **Immune to flux scheme** (geometry asymmetry is in the assembly, not the flux
  formula — Central is equivariant too) → matches "mixed flux didn't help."
- **Immune to order** at p1/p2 (the seed is the near-singular rupture front, where
  p-convergence is slow) → matches "higher order didn't help."
- **Immune to the y-mirror 'symmetric mesh'** → wrong axis.
- **Immune to over-int/resample** (operators already symmetric).
- Present in BOTH TPV102 and TPV31 → shared mesh geometry, not the friction law.

## THE FIX
Regenerate the **fault-plane triangulation to be strike(x)-symmetric about the
hypocenter (x=0)** — e.g. mesh the x≥0 half and reflect to x<0, or use a structured
x-symmetric fault triangulation — then extrude in y as now (keeping the existing
y-mirror). Pure mesh-generation change (`.geo`), no solver code.

Cheap CONFIRMATION before a full remesh: run the SAME mesh **mirrored in x** and
check the dip drift FLIPS SIGN (it must, if mesh x-asymmetry is the cause). One
Frontera A/B.

## Working diagnosis (revised — supported by local evidence)

- **Frame (H1): REFUTED** (Result 1, compiled). The rotation is clean; a frame
  freeze would be a no-op.
- **Seed: a genuine bulk `σ_yz` / `[[v_dip]]`** radiated into the near-fault
  region by **near-fault upwind-flux dissipation + mesh (Kuhn-dicing) asymmetry**
  (Result 3: friction-reaction-independent; over-int worsens it). In production,
  **free-surface SV/Rayleigh conversion** amplifies it near z=0 and the bilateral
  rupture makes it anti-symmetric and ~linear in strike — the fingerprint. The
  clean V∥τ friction then rectifies it into secular dip slip.
- **Amplifier (H2): the dip closed loop** (`fault_face_flux.cpp:238`) — a
  long-secular multiplier (NOT visible in the short local window; references lack
  a large enough seed to excite it).
- **Same SeisSol-parity gap** already tracked for the σ_n leak; the dip channel is
  the z-tangential survivor that the σ_n fixes (y-mirror mesh, central flux for the
  y-antisymmetric channel) did not address. (memory `project_upwind_unstructured_seissol_gap`:
  gap = over-int + resample + flux-near-fault + fault partition.)

## Proposed fixes (ranked; what the local evidence supports/excludes)

1. **Central/mixed flux near the fault** — PRIMARY. The seed is upwind-dissipation
   driven (over-int worsened it; the reference mixed-flux TPV31 already showed a
   smaller surface dip). Editable `dynamic/`. LOCAL VALIDATION: enable the
   mixed/central path in the slab harness → expect `[[v_dip]]` to drop (cheap, 46 s).
2. **Full dealiasing = over-int + resample TOGETHER** — over-int ALONE is excluded
   locally (×1.6 worse); the resample half is required and is only exercisable via
   the production substep-iterator path. Frontera A/B or a harness extension that
   routes friction through the iterator.
3. **Near-fault mesh symmetry / fault partition** — reduces the mesh half of the
   seed; impractical to make z-symmetric under a free surface.
4. EXCLUDED: frame freeze (no-op, Result 1); over-int alone (worse, Result 3);
   `SEAS_FORCE_V1_ZERO` as a general fix (antiplane-only; doesn't remove the seed).
