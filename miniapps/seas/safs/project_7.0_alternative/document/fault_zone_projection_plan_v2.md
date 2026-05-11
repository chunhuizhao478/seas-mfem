# Fault-Zone Projection-Quality Plan (v2 — implementation)

**Date:** 2026-05-09
**Status:** File-level implementation plan, derived from
`fault_zone_projection_plan_v1.md` + a debug pass that quantified the
streak amplitude on the 500/1000/2000 m meshes.
**Supersedes:** none.  Refines v1 with empirical measurements and
maps each phase to specific files, signatures, and tests.

---

## 0. Confirmed diagnosis (from /code-debug session, 2026-05-09)

The diagnosis stated in v1 §1 is correct in mechanism but not in the
location of the error.  Empirical reproduction on the existing
`data_projected/preview` outputs:

| Mesh   | H1 vertex error vs trilinear sidecar | Per-tet (cell-mean − centroid) L² | Streak amplitude in PNG |
|--------|--------------------------------------|-----------------------------------|-------------------------|
| 500 m  | ≤ 1.7 m/s (rounding noise)           | **258.8 m/s**                     | ~0.5 km/s               |
| 1000 m | ≤ 0.4 m/s                            | **359.0 m/s**                     | ~0.5 km/s               |
| 2000 m | ≤ 0.0 m/s                            | **427.8 m/s**                     | ~0.5 km/s               |

Measured on the z = −1 km ± 100 m slab, fault distance from the cut
STL.  The H1-P1 projection is **bit-correct at the vertices** — every
vertex value equals the trilinear evaluation of the sidecar at that
vertex's coordinates to within rounding noise.

**The streaks are entirely a sub-element interpolation artefact.**
ParaView (and MFEM's downstream consumers) sample the H1-P1 GF at
arbitrary points by linear interpolation between the four vertex
values of the enclosing tet.  Where a tet bridges a basin/basement
gradient, this linear interpolation paints values along the gradient
that exist nowhere in the source — the visible NW-SE streak is the
locus of those tets along the cut fault trace.

This refines the v1 attribution:
- **§1.A (basin/fault coincidence)** — confirmed: the gradient is
  there because the fault cuts through basin edges.
- **§1.B (fault-band smear)** — *the actual mechanism*: linear-in-tet
  interpolation between correctly-sampled vertices.  Cell-mean error
  grows with element size (258 → 428 m/s as mesh coarsens), exactly as
  expected for a representation error proportional to ‖∇Vs‖ × h.
- **§1.C (cross-fault continuity)** — structurally true (H1 cannot
  jump across a face) but **does not contribute extra error at the
  vertex level** in this dataset.  The sidecar is a single
  RectilinearGrid with no actual jump at the fault — there is no
  fault-perpendicular discontinuity for H1 to fail to represent.  C
  becomes a real contributor only if a future sidecar contains
  basin/basement contrasts on opposite fault sides.

**Implication for sequencing:** F-2 (DG-L2(0)) is the dominant fix in
this dataset because the dominant error is in-tet smear.  F-3 is
insurance against future per-side sidecars.  F-4 (size field) is the
correct response if the cell-mean residual at the desired mesh size is
still too large after F-2.

---

## 1. Open questions from v1 §8 — proposed answers

These are best-guess defaults, intended to be confirmed by the user
before implementation begins.

1. **Near-fault accuracy target.** Default: `[0–500 m]` band L² ≤
   200 m/s on Vs (≈ 7 % of median Vs in this dataset).  This is
   roughly 1× the F-2 baseline already achievable on the 500 m mesh
   without further refinement.
2. **Pre-existing fault-perpendicular Vs jumps.** No, not in
   `velocity_safs.h5`.  But the next CVM-H rebuild may add a basin /
   basement contrast across the SAF.  → F-3 is "cheap insurance" and
   should ship in the same batch as F-2.
3. **F-4 (mesh refinement) acceptable?** Defer; pursue only if F-2 +
   F-3 cell-mean residual exceeds the target for any production mesh
   we plan to ship.

---

## 2. Phase F-2 — DG-L²(0) projection (file-level)

**Goal:** add an alternative projection that writes **one constant
value per tet**, evaluated by sampling the trilinear sidecar at the
tet centroid (or by an L² volume average over a small cubature).
Default behaviour of `seas_project_velocity_to_mesh` is unchanged
(stays H1-P1) until the user opts in via a new flag.

### 2.1 `io/field_coefficient.{hpp,cpp}`

Add a sibling to `Project` / `ProjectVelocity` that targets an
L²(0) FES.

```cpp
// field_coefficient.hpp — additions, no removals.
//
// Project a single named field onto an L²(0) ParFiniteElementSpace by
// centroid sampling.  `target_fes` MUST be an L²(0) (= one DOF per
// element) ParFiniteElementSpace, asserted at runtime.
//
// Behaviour vs Project():
//   - bbox containment + post-projection range checks: unchanged.
//   - GF assembly: per element, sample the FieldCoefficient at the
//     element's reference centroid (`Geometries.GetCenter(geom_type)`)
//     and write to that element's single DOF.  No averaging — exact
//     centroid sample, matches the volume integral to first order on
//     a tet.
//   - call_count_ counter: incremented identically.
static Result ProjectDG0(
   const DataField3D& field,
   mfem::ParFiniteElementSpace& target_fes,
   real_t scale = 1.0,
   real_t offset = 0.0);

// Velocity composite, DG-0 variant.
static VelocityFields ProjectVelocityDG0(
   const std::string& sidecar_path,
   mfem::ParFiniteElementSpace& target_fes,
   real_t lambda_min_pa = 0.0,
   real_t lambda_max_pa = 1.0e12,
   real_t mu_min_pa     = 1.0e7,
   real_t mu_max_pa     = 1.0e11);
```

Implementation in `field_coefficient.cpp`:

```cpp
FieldProjector::Result FieldProjector::ProjectDG0(
   const DataField3D& field,
   mfem::ParFiniteElementSpace& target_fes,
   real_t scale, real_t offset)
{
   call_count_.fetch_add(1);

   // Enforce L²(0): one DOF per element.
   const mfem::ParMesh* pmesh = target_fes.GetParMesh();
   MFEM_VERIFY(pmesh != nullptr, "ProjectDG0: target_fes has null ParMesh");
   MFEM_VERIFY(pmesh->GetNE() == 0 ||
               target_fes.GetFE(0)->GetDof() == 1,
               "ProjectDG0: target_fes is not L2(0); got "
               << target_fes.GetFE(0)->GetDof() << " dofs/elem");

   // Standard bbox check (re-using the parallel path).
   real_t mxmin, mxmax, mymin, mymax, mzmin, mzmax;
   ComputeMeshBBoxParallel(*pmesh, mxmin, mxmax, mymin, mymax, mzmin, mzmax);
   if (!field.ContainsBBox(mxmin, mxmax, mymin, mymax, mzmin, mzmax))
   {
      AbortContainmentFailure(field, mxmin, mxmax, mymin, mymax,
                              mzmin, mzmax);
   }

   FieldCoefficient coef(field, scale, offset);
   auto gf = std::make_shared<mfem::ParGridFunction>(&target_fes);
   *gf = 0.0;

   const int ne = pmesh->GetNE();
   mfem::Vector xyz(3);
   mfem::Array<int> vdofs;
   for (int e = 0; e < ne; ++e)
   {
      mfem::ElementTransformation* T = pmesh->GetElementTransformation(e);
      const mfem::Geometry::Type gt = pmesh->GetElementBaseGeometry(e);
      const mfem::IntegrationPoint& ip = mfem::Geometries.GetCenter(gt);
      T->Transform(ip, xyz);
      const real_t v = scale * field.Evaluate(xyz[0], xyz[1], xyz[2])
                     + offset;
      target_fes.GetElementDofs(e, vdofs);
      MFEM_ASSERT(vdofs.Size() == 1, "ProjectDG0: expected 1 dof/elem");
      (*gf)(vdofs[0]) = v;
   }

   real_t lo, hi;
   ComputeGridFunctionMinMaxParallel(*gf, lo, hi);
   const real_t expected_lo = scale * field.MinValue() + offset;
   const real_t expected_hi = scale * field.MaxValue() + offset;
   const real_t lo_band = std::min(expected_lo, expected_hi);
   const real_t hi_band = std::max(expected_lo, expected_hi);
   if (lo < lo_band || hi > hi_band)
   {
      AbortRangeFailure(field.FieldName(), lo, hi, lo_band, hi_band);
   }
   Result r;
   r.gf = gf;
   r.min_value = lo;
   r.max_value = hi;
   return r;
}
```

`ProjectVelocityDG0` is the same skeleton as `ProjectVelocity`,
calling `ProjectDG0` instead of `Project`.  Important detail: in DG-0,
deriving `µ = ρVs²` and `λ = ρVp² − 2µ` per-DOF is **exact** (each tet
holds one centroid sample of each of ρ, Vp, Vs).  No subtle
inconsistency vs sampling-then-deriving in H1 (where the algebra is
applied vertex-wise — also exact, but not invariant under
interpolation).

### 2.2 `dynamic/heterogeneous_material.hpp`

`MaterialField::At(elem, dof, ...)` already works for L²(0):
`fes->GetElementDofs(elem, vdofs)` returns a single index, and `dof =
0` indexes into it.  But the contract is fragile — a future caller
that passes `dof > 0` would assert.  Add an optional fast-path:

```cpp
// New mode tag — strictly equivalent to GridFunction at runtime
// for L²(0) GFs, but documents the invariant and lets the wave op
// short-circuit per-DOF lookups:
enum class Mode : int
{
   Constant         = 0,
   GridFunction     = 1,
   DG_L2_GridFunction = 2,  // NEW; backed by L²(0) GFs (1 dof / elem)
};

static MaterialField MakeDG_L2_GridFunction(
   std::shared_ptr<mfem::ParGridFunction> rho_gf,
   std::shared_ptr<mfem::ParGridFunction> lambda_gf,
   std::shared_ptr<mfem::ParGridFunction> mu_gf);

inline void At(int elem, int dof, real_t& la, real_t& mu, real_t& rho) const
{
   if (mode == Mode::Constant) { /* unchanged */ return; }
   if (mode == Mode::DG_L2_GridFunction)
   {
      // Each GF has 1 DOF per element; vdof = elem in L²(0) ordering.
      la  = (*lambda_gf)(elem);
      mu_ = (*mu_gf)(elem);
      rho = (*rho_gf)(elem);
      return;
   }
   /* fall through to GridFunction path — unchanged */
}
```

`MaxCpInElement(elem)` in `DG_L2_GridFunction` mode returns
`sqrt((la + 2µ) / ρ)` evaluated at the single element value — no
loop.

### 2.3 `drivers/project_velocity_to_mesh.cpp`

Add a `--basis` option, default `h1p1` (back-compat).  When `dg0`,
swap the FEC and call `ProjectVelocityDG0`.

```cpp
// new option:
const char* basis = "h1p1";
args.AddOption(&basis, "-b", "--basis",
               "Material projection basis: h1p1 (default, "
               "smooth Lagrange) | dg0 (per-element constant; "
               "centroid-sampled).");

// after Parse():
const std::string basis_str(basis);
const bool use_dg0 = (basis_str == "dg0");
if (basis_str != "h1p1" && !use_dg0)
{
   if (rank == 0)
   {
      std::cerr << "ERROR: --basis must be 'h1p1' or 'dg0'.\n";
   }
   return 1;
}

std::unique_ptr<mfem::FiniteElementCollection> fec;
if (use_dg0)
{
   fec = std::make_unique<mfem::L2_FECollection>(
            /*order=*/0, pmesh.Dimension(),
            mfem::BasisType::GaussLobatto);
}
else
{
   fec = std::make_unique<mfem::H1_FECollection>(order, pmesh.Dimension());
}
mfem::ParFiniteElementSpace fes(&pmesh, fec.get());

auto vf = use_dg0
   ? mfem::seas::FieldProjector::ProjectVelocityDG0(sidecar_path, fes)
   : mfem::seas::FieldProjector::ProjectVelocity   (sidecar_path, fes);
```

ParaView output is unchanged structurally — the registered fields are
still `Vp / Vs / density / lambda / mu`, just on an L²(0) FES.  The
visual difference is that each tet renders as a single colour,
revealing the mesh tessellation directly (which is the point — see
the regenerated comparison PNG).

### 2.4 Tests

- **T-F2-1 (unit, serial)** in
  `tests/unit/test_field_projector_dg0.cpp`: synthetic constant
  sidecar `Vs ≡ 2500`, mesh = 8-tet unit cube → all 8 cell DOFs equal
  2500 to within 1 ULP.  Mirrors the existing constant-baseline test
  in `test_field_projector.cpp`.
- **T-F2-2 (unit, serial)**: synthetic linear sidecar `Vs(x,y,z) =
  100x + 200y + 300z`, mesh = 8-tet unit cube → each cell DOF equals
  the closed-form centroid value.  Tolerance 1e-12.
- **T-F2-3 (parallel, np=4)**: same sidecar + small unit cube,
  partitioned; verify global L² unchanged.
- **T-F2-4 (regression on real sidecar)**: new
  `tests/integration/test_project_velocity_safs_dg0.cpp` (gated on
  the SAFS sidecar fixture if present): assert that, on a coarse
  20 km × 20 km × 5 km box at 2 km tet, the per-element values are
  within 1.0 m/s of an independent NumPy centroid sample.
- **T-F2-5 (driver smoke)**: `make seas_project_velocity_to_mesh`
  with `--basis dg0` writes a `.pvd` whose `<DataArray Name="Vs">` has
  exactly `pmesh.GetNE()` entries (instead of `GetNV()`).

### 2.5 Acceptance — F-1 metric ported in (replaces v1 §F-1)

`safs/<project>/code_preprocess/projection_quality/fault_zone_metric.py`
(new) computes the F-1 stratified metric over the cut-STL distance
bands.  After F-2 alone, on the **500 m mesh**:

- `[0–500 m]` band cell-mean residual L² (DG-0 GF vs centroid
  trilinear): **target ≤ 50 m/s** (cf. 258.8 m/s for H1-P1
  cell-mean — should drop to numerical-noise floor since each cell
  is *defined* by its centroid sample).
- Bulk band `[> 3000 m]`: unchanged from H1-P1 baseline (also at
  noise floor).

The cell-mean residual drop is the headline number F-2 is judged on.

### 2.6 Cost & risk

- ~120 LOC across `field_coefficient.{hpp,cpp}`,
  `heterogeneous_material.hpp`, `drivers/project_velocity_to_mesh.cpp`.
- ~80 LOC of unit/integration tests.
- No edits to files in the C2-invariant set
  (`bp5/`, `bp1/`, `bp2/`, `domain/`, `friction/dieterich_ruina.hpp`,
  `solver/`, `fault/fault_basis.hpp`).
- Risk: existing wave-op consumers don't yet subscribe to the
  `MaterialField` of the projected GF — that's a separate phase
  (data-projection-feature plan v2 Phase 5).  F-2 lands the
  **projection** path and the **preview** output; consuming the
  projected material is gated by Phase 5 wiring.

---

## 3. Phase F-3 — Per-side fault projection

**Goal:** for every tet that has a face on the cut fault, sample the
sidecar in a half-space that's clearly on that tet's side of the
fault, so a hypothetical basin/basement contrast across the fault is
not averaged together.  Builds on F-2; meaningless on H1.

### 3.1 Source of the fault face list at projection time

The projector currently runs *before* any `WaveOperator` exists
(`project_velocity_to_mesh.cpp` is a standalone driver).  Two routes:

1. **Scan boundary attributes directly from the mesh.** The cut fault
   is tagged with attribute `bc.fault_attr` (typically Physical Surface
   3 in the gmsh file).  Re-implement the same loop already in
   `wave_operator.inl:286–295` here:

   ```cpp
   for (int b = 0; b < pmesh.GetNBE(); ++b) {
      if (pmesh.GetBdrAttribute(b) != fault_attr) continue;
      const int face_idx = pmesh.GetBdrElementFaceIndex(b);
      ...
   }
   ```

   `fault_attr` is a new CLI option, default 3 (Tandem convention; see
   `miniapps/seas/CLAUDE.md` "Tandem mesh tags").

2. **Construct a throwaway WaveOperator just for its fault list.**
   More expensive; not recommended.  Use route (1).

### 3.2 `io/field_coefficient.{hpp,cpp}` — new entry point

```cpp
// Apply F-3 on top of an existing DG-0 projection.
//
// Inputs:
//   - mu_gf, lambda_gf, rho_gf:  L²(0) GFs already populated by
//                                ProjectVelocityDG0.
//   - sidecar_path:              same sidecar (re-loaded for sampling).
//   - fault_attr:                boundary attribute tagging the fault.
//   - eps_normal_m:              half-step into the half-space, metres.
//                                Default = 0.05 × LC_NEAR (e.g. 75 m
//                                for LC_NEAR = 1500 m).  Tunable.
//
// Behaviour:
//   For each (b, face_idx) on the fault:
//     - get face_idx → ElemTrans1, ElemTrans2 from
//       pmesh.GetFaceElementTransformations(face_idx).  In our
//       ParMesh setup interior fault faces have both elem1 and
//       elem2 local; shared fault faces have only elem1 + a
//       `Mesh::FaceInfo` neighbour.  Skip the neighbour side on
//       shared faces (that rank handles its own elem1).
//     - compute outward unit normal n_+ on Elem1 side at the face
//       centroid (use `CalcOrtho` then normalise).
//     - per-side query points:
//         q_+ = centroid(Elem1) + eps × +n
//         q_− = centroid(Elem2) + eps × −n   (where Elem2 is local)
//     - overwrite that element's L²(0) DOF with the trilinear
//       evaluation of the sidecar at q_±.
static void ApplyPerSideFaultProjection(
   const std::string& sidecar_path,
   mfem::ParGridFunction& vp_gf,
   mfem::ParGridFunction& vs_gf,
   mfem::ParGridFunction& rho_gf,
   mfem::ParGridFunction& lambda_gf,
   mfem::ParGridFunction& mu_gf,
   int fault_attr,
   real_t eps_normal_m = 75.0);
```

Implementation outline:

```cpp
mfem::ParMesh* pmesh = vp_gf.ParFESpace()->GetParMesh();
DataField3D vp(sidecar_path, "Vp"), vs(sidecar_path, "Vs"),
            rho(sidecar_path, "density");

mfem::Vector xyz(3), nor(3);
mfem::Array<int> vdofs;
for (int b = 0; b < pmesh->GetNBE(); ++b)
{
   if (pmesh->GetBdrAttribute(b) != fault_attr) continue;
   const int face_idx = pmesh->GetBdrElementFaceIndex(b);
   auto* T = pmesh->GetFaceElementTransformations(face_idx);
   if (!T) continue;

   // Outward face normal at face reference centroid:
   const mfem::Geometry::Type fgeom =
      pmesh->GetFaceBaseGeometry(face_idx);
   const mfem::IntegrationPoint& fip = mfem::Geometries.GetCenter(fgeom);
   T->SetAllIntPoints(&fip);
   mfem::CalcOrtho(T->Jacobian(), nor);
   const real_t nlen = nor.Norml2();
   nor /= nlen;   // outward from Elem1 side

   auto sample_at = [&](const mfem::IntegrationPoint& ipref,
                        mfem::ElementTransformation& Te,
                        const mfem::Vector& step) -> mfem::Vector {
      mfem::Vector x(3);
      Te.Transform(ipref, x);
      x += step;
      mfem::Vector out(3);
      out(0) = vp.Evaluate(x(0), x(1), x(2));
      out(1) = vs.Evaluate(x(0), x(1), x(2));
      out(2) = rho.Evaluate(x(0), x(1), x(2));
      return out;
   };

   auto write = [&](int elem, const mfem::Vector& pvr) {
      vp_gf.ParFESpace()->GetElementDofs(elem, vdofs);
      const int idx = vdofs[0];
      vp_gf(idx)  = pvr(0);
      vs_gf(idx)  = pvr(1);
      rho_gf(idx) = pvr(2);
      const real_t mu_e = pvr(2) * pvr(1) * pvr(1);
      const real_t la_e = pvr(2) * pvr(0) * pvr(0) - 2.0 * mu_e;
      mu_gf(idx)     = mu_e;
      lambda_gf(idx) = la_e;
   };

   const mfem::Geometry::Type gt1 =
      pmesh->GetElementBaseGeometry(T->Elem1No);
   mfem::Vector step1 = nor; step1 *= eps_normal_m;
   write(T->Elem1No, sample_at(mfem::Geometries.GetCenter(gt1),
                               *T->Elem1, step1));

   if (T->Elem2No >= 0)   // local interior fault face
   {
      const mfem::Geometry::Type gt2 =
         pmesh->GetElementBaseGeometry(T->Elem2No);
      mfem::Vector step2 = nor; step2 *= -eps_normal_m;
      write(T->Elem2No, sample_at(mfem::Geometries.GetCenter(gt2),
                                  *T->Elem2, step2));
   }
   // Shared faces: the other rank holds Elem2 — it will run its
   // own loop and will see the symmetric face from its side.
}
```

### 3.3 Tests

- **T-F3-1**: synthetic sidecar with a discontinuous Vs across y=0:
  `Vs(x,y,z) = 1500 if y < 0 else 3500`.  Cube mesh [0,2]³ with a
  `y = 1` internal-Dirichlet face tagged `fault_attr = 7`.  After
  ApplyPerSideFaultProjection: every cell on `y < 1` ≈ 1500, every
  cell on `y > 1` ≈ 3500 (no smear).
- **T-F3-2** (parallel np=2): same fixture, partition cuts the fault
  → check that shared-face elements are still set correctly via the
  rank-local Elem1 side.
- **T-F3-3** (real-data, optional): on the SAFS 1000 m mesh, apply
  F-3 and verify the cell-mean residual on the `[0–500 m]` band
  decreases by ≥ 2× vs F-2 alone (target — depends on whether the
  current sidecar has any fault-perpendicular jump, which v0 does
  not, so this test is initially expected to be a no-op proving F-3
  is bound-preserving on smooth data).

### 3.4 Acceptance

- F-3 is **identity** on smooth sidecars (test T-F3-3 passing
  trivially confirms this).
- F-3 on a synthetic discontinuous sidecar (T-F3-1) eliminates 100 %
  of the cross-fault smear within `eps_normal_m`.
- Bulk error stays at noise floor.

### 3.5 Cost & risk

- ~150 LOC including the test fixtures.
- No edits to C2-invariant files.
- Risk: outward-normal sign convention.  `CalcOrtho` returns
  unnormalised; the *sign* depends on Elem1's orientation.  Robust
  test: T-F3-1 above directly catches a sign flip (the basin/basement
  values would swap).

---

## 4. Phase F-4 — Sidecar-driven background size field

Skipped in v2; pursue only if F-2 + F-3 leave residual error above
target on a production mesh.  v1 §F-4 design unchanged.

---

## 5. Recommended sequence

| Order | Phase | Estimated LOC | Estimated effort |
|-------|-------|--------------:|-----------------:|
| 1 | F-1 metric (Python) | 200 | 0.5 day |
| 2 | F-2 DG-L²(0) projection | 200 | 1 day |
| 3 | (assess) regenerate comparison PNG with F-2; check acceptance | — | 0.5 day |
| 4 | F-3 per-side projection | 250 | 1 day |
| 5 | F-5 convergence study + ship | — | 0.5 day |

**Stop-gate after 3.** If F-2 alone hits the cell-mean residual
target on the 500 m mesh, F-3 is deferred until a per-side sidecar
arrives.  Per the v1 plan §6 recommendation, do **not** jump to F-4
or T6 first.

---

## 6. Non-goals / explicit deferrals

- Wiring the projected `MaterialField` into the `WaveOperator`'s
  ADER/RK4 loops.  That is `data_projection_feature_plan_v2.md`
  Phase 5; it is **logically independent of F-2/F-3** and should be
  scheduled separately.  F-2/F-3 are useful even with the wave-op
  still configured for `Mode::Constant` — they unblock the preview
  workflow and let us track projection quality before we light up
  the heterogeneous run.
- Replacing the H1-P1 path.  H1-P1 stays as the default and remains
  the only path used by tests until F-2 is reviewed and lands.

---

## 7. Files touched (summary)

```
miniapps/seas/io/field_coefficient.hpp                  +35 lines
miniapps/seas/io/field_coefficient.cpp                  +180 lines
miniapps/seas/dynamic/heterogeneous_material.hpp        +30 lines
miniapps/seas/drivers/project_velocity_to_mesh.cpp      +30 lines
miniapps/seas/tests/unit/test_field_projector_dg0.cpp   +200 lines (new)
miniapps/seas/tests/unit/test_field_projector_fault.cpp +180 lines (new)
miniapps/seas/Makefile                                  +6 lines
safs/<project>/code_preprocess/projection_quality/
   fault_zone_metric.py                                 +250 lines (new)
```

C2-invariant set untouched.
