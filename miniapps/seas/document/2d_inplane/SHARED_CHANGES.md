# 2D in-plane branch: changes to shared files

Branch `2d-inplane`, from `fully-dynamic-cycle` at 3926e9c8 (design report rev 3,
`/mnt/project-files/2d-inplane/2d-inplane-design-rev3.pdf`).  Every phase
lists here each shared file it touches, what the change does, what stays
the same without the new inputs, and how that was checked.  The 3D parity
set (`tests/2d_inplane/parity_3d.sh`) must stay bit for bit unchanged after
every phase, apart from the commit stamp in file headers.  The one exception
is a QD solve on more than one rank, which MUMPS does not reproduce bit for
bit from one run to the next even with the same binaries; those runs must
differ from the base by no more than repeated runs of the same binaries
differ from each other.

## Phase A, part 1: one fault orientation system (A1-A3)

New files (branch-owned): `fault/fault_frame_field.hpp`, `cycle/cycle_frames.hpp`,
`tests/2d_inplane/` (field test, parity script, makefile), this file.

| Shared file | Change | Without the new input |
|---|---|---|
| `domain/domain_config.hpp` | `FaultFrameSource` (field or config, forward-declared types) and `DomainConfig::fault_frames` | empty source: the constants |
| `domain/elasticity_operator.hpp` | legacy constructor takes a trailing `FaultFrameSource` (default empty); `FrameField()`; `ComputeSkeletonDirichletSign` takes the creep census normal in field mode; helpers `FaceRefNormal_`, `SlipSignBR2_`; member `up_` | reference normal (0, -1, 0), up (0, 0, 1), BR2 sign `nor_q(1) > 0`, as before |
| `domain/elasticity_operator_setup.inl` | `SetupFrameField_` (collective): builds the field from the operator's own fault and interior Dirichlet face lists, or checks a given field against them; per-face reference normals into `FaultBasis` | not called into: same calls with the same constants |
| `domain/elasticity_operator_traction.inl` | `GetFaultDOFBasis` takes `up_` and the per-face reference normal | the same constants |
| `domain/elasticity_operator_assembly.inl` | BR2 slip sign through `SlipSignBR2_` | `(nor_q(1) > 0) ? 1 : -1` |
| `fault/fault_basis.hpp` | `Compute`, `AppendSharedFaces`, `ComputeQPBasis`, `ComputeQPBasisShared` take an optional per-face reference normal array; `ComputeOrientedFrame` unchanged | `face_ref = nullptr`: the single reference normal |
| `domain/nitsche_wave_operator.hpp` | the elastic copy builds the field from the fault setup's `[fault_geometry]` (or reads one from the `DomainConfig`); the Kelvin-Voigt copy reads the same field; `SetupFault_` builds no census of its own, checks the configuration and `eps_f == sigma_f`; `Orientation()` returns the field's fault census; `FrameField()` | the same census (same faces, order and configuration) and the same frames, computed once inside the field |
| `cycle/qd_phase.hpp` | the QD operator gets `CycleFrameConfig(cfg)` | - |
| `cycle/dyn_phase.hpp` | optional trailing `frames` (the QD regime's field); fault setup from `CycleFrameConfig`; `NormalMismatch <= 1e-12` asserted | null: the Nitsche operator builds the same field |
| `drivers/fd_cycle_driver.cpp` | hands `qd.Domain().FrameField()` to `DynPhase` | - |

Why it is bit for bit on today's meshes: the QD operator's frame arithmetic
(`FaultBasis::ComputeOrientedFrame`, per face, per quadrature point, per
DOF) is unchanged; only its reference normal becomes the face's oriented
normal n from the census and its up vector the deck's.  The arithmetic reads
the reference normal only through the sign of `n_raw . ref`, which is the
same for ref = (0, -1, 0) and ref = n wherever the census orients the fault
along (0, -1, 0) (every BP5 and TPV deck: `ref_normal = [0, -1, 0]`,
`up = [0, 0, 1]`).  The creep sign and the BR2 sign are the same dot product
with the same sign.  The Nitsche census is the call `SetupFault_` made
before, on the same face lists.

How it was checked (3D parity set, cloud container, 10-24): TPV205 and TPV104 on 2 ranks
byte for byte; cycle QD and dynamic regimes on 1 rank and the dynamic regime
on 2 ranks byte for byte apart from the commit stamp; cycle QD on 2 ranks
8.9e-10 from the base (two base runs: 8.6e-10; two branch runs: 3.7e-9);
BP5 quick check on 8 ranks 3.8e-8 from the base (two base reruns: up to
9.3e-8); unit tests unchanged.  Differences are scaled by each column's
largest value.  The fully-dynamic-cycle test targets
(`tests/2d_inplane/cycle_tests.sh`, BLR deck copy): API pin, P2 frame map
on 1, 2, 4, 8 ranks, P2 on 2, P3 and P4 on 1 and 2 pass with both binary
sets; their printed measurements are the same on 1 rank and differ at
rounding level, inside the tests' tolerances, on 2.

Known behaviour changes: (1) the QD regime of the cycle now runs the
orientation census (S3 side test, planarity, pin margin) on its fault and
creep faces and refuses a mesh that fails it, where it ran no check before;
(2) the cycle refuses a fault and a creep surface that the one pin orients
inconsistently across their shared edges (seam rule) or that meet at a
junction edge; (3) the creep-plane jump now takes its side from the deck's
`[fault_geometry].ref_normal`, while the far-field function
`MakeBP5DirichletFunc` stays the y = 0 BP5 function: a deck with
ref_normal (0, 1, 0) would get a creep jump of the opposite sense to its
walls (every deck in config/fd_cycle uses (0, -1, 0); part 2's
`[cycle.loading]` refuses the combination); (4) `cycle-owned-check` (the
fully-dynamic-cycle branch's rule M1) lists the shared files above: expected
on this branch.

## Phase A, part 2: the plate loading (A4), the Okada check, the rotation run

New files (branch-owned): `cycle/loading.hpp`; `tests/2d_inplane/frame_box.hpp`
(the test box and rotations, moved out of the field test),
`test_plate_loading.cpp`, `okada92.hpp` (Okada's DC3D in C++),
`test_okada_dislocation.cpp`, `rotate_case.py`, `rotation_run.sh`,
`rotation_compare.py`.

| Shared file | Change | Without the new input |
|---|---|---|
| `cycle/cycle_config.hpp` | `LoadingConfig` and `CycleConfig::loading`: the `[cycle.loading]` table (`kind = "bp5"` or `"plate"` with `plane_point`, `plane_normal`, `slip_dir`, `band`); `band_sign = sgn(-plane_normal . ref_normal)`; checks (unit vectors, slip direction in the plane, plane normal not in the fault plane); kind `"bp5"` requires `ref_normal = (0, -1, 0)` | no table: kind `"bp5"` |
| `domain/elasticity_operator.hpp` | legacy constructor takes a trailing `dirichlet_func` (default empty) | `MakeBP5DirichletFunc(Vp)`, as before |
| `cycle/qd_phase.hpp` | the QD operator gets `MakeCycleLoading(cfg)` | kind `"bp5"`: `MakeBP5DirichletFunc(plate_rate)`, the function it had |
| `cycle/dyn_phase.hpp` | the dynamic regime's loading is `MakeCycleLoading(cfg)` | the same `MakeBP5DirichletFunc(plate_rate)` |
| `tests/2d_inplane/parity_3d.sh` | `decks OUT` mode (only the decks, for the rotation run) | - |

The band value.  Inside the band the loading value is the jump that the
prescribed-slip (creep) faces impose, and the operator multiplies it by
`ComputeSkeletonDirichletSign` (+1 where the face's Elem1-outward normal
points along the oriented normal n), so the imposed jump is
u(side opposite n) - u(side n points to).  The plate motion's jump across
the plane is u(side of m) - u(other side) = Vp t s, which is that jump when
m points away from n and its negative when m points along n: the band value
is sgn(-m . ref) Vp t s (design rev 3 wrote Vp t s, which is right only for
m . ref < 0, as in BP5).  The value does not change under (m, s) -> (-m, -s).

The ψ0 and initial slip direction items of A3 needed no change: the per-DOF
V_init (dip, strike) is already the deck's (`[cycle.fields]`, H1b
`SetInitialFieldsPerDOF`, which overwrites the constructor's (0, 1)
direction), and the SCEC-mode ψ0 is computed per DOF from the deck's f0, b,
V0 and V_background (uniformity of those is asserted), so neither depends on
the fault's orientation.

How it was checked (cloud container, 10-25).  Plate test (1 and 2 ranks):
kind "plate" with BP5's vectors equals `MakeBP5DirichletFunc` bit for bit at
65156 point-times on the BP5 mesh's Dirichlet boundary; the operator's
solution with the walls and the creep band loaded and the fault slipping by
the same jump is the rigid two-block motion to 4.4e-14 on a flat, a 60°
dipping and a generally rotated plane (O(1) with the band sign reversed);
six refusals.  3D parity (as part 1): TPV205, TPV104, cycle QD and dynamic
on 1 rank and dynamic on 2 ranks byte for byte against the base and part 1
(apart from the commit stamp); cycle QD on 2 ranks 4.8e-9 from the base
(repeat runs of one binary set: 8.6e-10 and 3.7e-9); BP5 quick check on 8
ranks 1.1e-7 (base reruns: up to 9.3e-8); unit tests exit 0, their logs
differing only in the order of ranks' lines and in a 2-rank MUMPS residual
at 1e-12.

Known behaviour changes: (1) a cycle deck without `[cycle.loading]` and with
a `[fault_geometry].ref_normal` other than (0, -1, 0) is now refused (the
BP5 function is the y = 0 plane's); every deck in config/fd_cycle passes.

Okada check (cloud container, 10-25; tests only).  okada92.hpp equals
Okada's Fortran DC3D (REAL*8 build) bit for bit at 1e5 points.
test_okada_dislocation: 39 of 39 on 2 ranks, 29 of 29 in quick mode on 1
rank.  60-degree dip, dip slip, p = 2, n = 4 to 8: surface displacement
3.1 % to 0.33 %, the traction's stress part 32 % to 9.1 %, the full
traction (with the penalty term) 877 % to 88 %: the uniform-slip patch's
edges are singular and pollute the penalty term, so the test bounds the
stress part and asks the full traction to converge.

Rotation run (cloud container, 10-25): the parity set's cycle decks run
unrotated with kind "plate" and BP5's vectors (33 files byte for byte the same as the
parity runs with the default loading, QD and dynamic), rotated by 40
degrees about (1, 2, 3) (mesh nodes and every deck vector), and unrotated
on 2 ranks as the control.  Largest difference, vectors measured together
(rotation_compare.py): QD 10 steps 7.9e-10 rotated / 9.1e-10 control;
dynamic 40 steps 3.6e-11 / 3.6e-12; QD-to-dynamic switch (V_on 3.5e-2 m/s,
then 40 dynamic steps) 1.5e-8 / 2.0e-8, switch at t = 0.421874329 s
unrotated and 0.421874327 s rotated.  The QD regime's adaptive steps carry
roundoff into the step sizes; a change of rank count does the same.  The
checkpoint's volume displacement matches within each tet (3123 of 62874
tets list their DOFs in another order on the rotated mesh: MFEM's
Mesh::MarkTetMeshForRefinement orders a tet's vertices by edge length, and
equal lengths differ in the last bit after the rotation).

## Phase B, part 1: the dimension switch (code and tests)

The dimension now comes from the mesh (design rev 3, sections 5.1 to 5.3,
rules R1 to R4).  On a 2D triangle mesh in the (x, y) plane the mesh loader,
the orientation census, the frame field, the face quadrature, the fault DOF
layout and the frames run with the same code, one dimension down; every
fault-local vector keeps its three slots and the out-of-plane entries are
exact zeros (rule R2).  The 2D stiffness, loads, traction and solves are
phase E: each of those entry points refuses a 2D mesh with a message that
names the phase.  `NumSlipComponents()` stays 2 (decision D2 is open).

New file (branch-owned): `tests/2d_inplane/test_2d_fault_layout.cpp` (and
its make target `test-2d-fault-layout`).

| Shared file | Change | 3D |
|---|---|---|
| `fault/fault_orientation.hpp` | edge mode: a segment face keyed by two vertex keys (slot 2 = -1); census entities are vertices (m = 1 a tip, m >= 3 a junction; new report field `tip_edges`); tau = +1 at the head k1, -1 at the tail k0; nC = t x e_z / \|t\|, length for area; point pin at the segment midpoint; J3 parity written for simplices with dim + 1 vertices (`DropVertex_`); S3 margin with the triangle's expected h_perp / 3; the gathered coordinates padded with z = 0 | the same arithmetic in the same order; new: J3 refuses a non-simplex element next to the fault (before, it read its first 4 vertices) |
| `fault/fault_frame.hpp` | `FaceVertices` with 2 vertices: the F6 fallback uses the segment's one edge | unchanged (3 edges) |
| `fault/face_quadrature.hpp` | the nodal face element is an `H1_SegmentElement` on SEGMENT faces (p + 1 nodes; order 0: the midpoint) | `H1_TriangleElement`, as before |
| `fault/fault_basis.hpp` | `ComputeOrientedFrame`: one rule in both dimensions on the padded normal (nx, ny, 0) and up (z = 0 for a 2-entry up); the old 2D branch (tangent1 = 90-degree rotation of n, no degeneracy guard) is removed; `ProjectTraction`, `EmbedSlip`, `EmbedSlipQP` fill and read both slots in 2D too | the same products: `n_raw(2)` and `up(2)` are read as before |
| `domain/elasticity_operator.hpp` | `dim_` from the mesh (2 or 3, SpaceDimension = Dimension; BR2 refused in 2D); `NumComponents()` and `Dimension()` return it; `ComputeSkeletonDirichletSign` loops over `dim_`; `Require3D_` | `dim_ = 3` |
| `domain/elasticity_operator_setup.inl` | DG space with `dim_` components; fault-face geometry checks, the face quadrature and the IP layout on SEGMENT faces in 2D (barycentric match (1 - x, x) for p >= 2); a 2D mesh with fault or prescribed-slip faces needs the frame field (no constants); the y = 0 audit pads z = 0 | TRIANGLE everywhere, as before |
| `domain/elasticity_operator_traction.inl` | depth and x3 from the last coordinate (2D: -y); DOF coordinates padded with z = 0; `GetFaultDOFBasis` takes 2 or 3; `Require3D_` in `ComputeTraction`, `AssembleSlipOnlyRHS`, `AssembleDirichletOnlyRHS` | `-coords(2)`, as before |
| `domain/elasticity_operator_{assembly,fault_ops,qd,penalty}.inl` | `Require3D_` at the top of `AssembleStiffness`, `BuildFaultFaceOps`, `BuildQDStageOperators`, `BuildQDGreen`, `GetFaultFacePenalties` | no-op |
| `cycle/scaled_mesh.hpp` | `LoadScaledMesh` reads a 2D triangle mesh (in its own space) or a 3D tetrahedral one and scales `dim` coordinates; `BdrAttrCount::triangles` renamed `simplex_faces` (segments in 2D; no user outside this file); `CountBoundary` loops over the space dimension; `FaultDOFsPerFace(order, dim = 3)` | the same reads and products |
| `tests/unit/test_fault_basis.cpp` | the two 2D cases pass 2-slot arrays to `ProjectTraction` / `EmbedSlip` (they passed 1-slot arrays, which the 2-slot calls would overrun), and check that the out-of-plane slot is exactly 0 and adds nothing to the jump; expected values unchanged | - |

The slot rule (derivation 10.2).  With up in the plane (BP3: up = [0, 1],
y up) the frame is s = up x n / \|up x n\| = -+e_z and d = s x n = the
down-dip tangent: slip lives in the dip slot (tangent1).  With up = [0, 0, 1]
(TPV205-2D, TPV105-2D: a horizontal slice) s = e_z x n lies in the plane and
d = -e_z: slip lives in the strike slot (tangent2).  A fault parallel to up
has no frame: the dynamic field uses the F6 fallback, the QD basis refuses
(`ComputeOrientedFrame`'s degeneracy guard), as in 3D.

How it was checked (cloud container, 10-25).  test_2d_fault_layout, 58 of 58
on 1 rank and 55 of 55 on each rank at 2 and 4 (the F6 check runs on rank 0
alone), and seven refusals: the loader on a 2D triangle mesh
(counts, tags, scale bit for bit); a BP3-shaped fault (free surface to
W_f, its creep segment below) at dips 30, 60 and 90 degrees: one patch, 2
tips, 0 junctions, 1 seam vertex, n = (sin, cos, 0) to 4.5e-16, the plus
element on the footwall, the report the same at every rank count; a
Y-shaped fault: 1 junction vertex, 3 tips, 3 patches; the 2D slice of the
phase A box: the same fault census, normals, frames (bit for bit) and plus
side as the 3D box; the F6 fallback where a segment is parallel to an
in-plane up (s along the segment, d = +-e_z); up = [0, 1] and [0, 0, 1]: the out-of-plane slot exactly
0 in the frames, `GetFaultDOFBasis`, `ProjectTraction` and `EmbedSlip`, and
the kernel triad = sigma_f (n, d, s); the fault DOF layout at p = 1, 2, 3:
p + 1 nodes per segment, every segment owned once, the owned -> local
scatter exact on ghost copies, and the (position, n, d, s) table at 2 and 4
ranks equal to the 1-rank one.  Refusals: a quadrilateral mesh, a 2D fault
without the frame field, up off the plane, ref_normal off the plane, a 2D
solve, up parallel to n in the QD basis, BR2 in 2D.

3D parity against phase A (0322c7e7, same container, clean rebuild of both
sets): TPV205 and TPV104 on 2 ranks, the cycle QD and dynamic regimes on 1
rank and the dynamic regime on 2 ranks byte for byte apart from the commit
stamp; cycle QD on 2 ranks 4.0e-9 (repeat runs of one binary set: up to
3.7e-9; part 2 against the base: 4.8e-9); BP5 quick check on 8 ranks 4.4e-8 (repeat runs: up to 9.3e-8); the
unit tests pass, test_fault_basis with its 2 new checks (221 of 221) and
otherwise the same output.  The fully-dynamic-cycle tests pass; their logs
are the same on 1 rank apart from timings and differ at rounding level on 2.
The 2-rank P3 run was killed once by the container's 14.3 GB memory limit;
rerun alone, it passes with both binary sets, and both reach the limit
(page cache included), so it sits at the edge of this container either way.

Known behaviour changes in 3D: (1) J3 refuses an element next to the fault
that is not a tetrahedron (before, it read the element's first 4 vertices);
(2) `BdrAttrCount::triangles` is renamed `simplex_faces`.

## Phase B, part 2: 2D deck vectors and the three 2D meshes

A 2D deck is a 3D deck whose vectors have two entries (design rev 3,
section 5.5).  The readers take [x, y] and pad z = 0, record how many
entries each key had, and a check against the mesh dimension, called by
both drivers right after the mesh is read, refuses a vector of the other
dimension (never truncated or padded silently).  The one 3-entry vector a
2D deck takes is `[fault_geometry].up = [0, 0, +-1]` (out of the plane, a
horizontal slice: the slip in the strike slot); an in-plane up is [ux, uy]
(a vertical section: the slip in the dip slot).  In 2D, ref_normal and up
are required.

New files (branch-owned): `tests/2d_inplane/test_2d_deck.cpp` and its make
target `test-2d-deck`; the meshes below (each `.geo` and the `.msh` it makes,
added with `git add -f` past the folder's `*.msh` ignore rule, as the other
committed meshes); in `tests/2d_inplane/test_2d_fault_layout.cpp` the tests
M1 (the committed meshes) and D2 (the layout on the BP3 mesh).

| Shared file | Change | 3D |
|---|---|---|
| `spatial/code/spatial_friction.hpp` | `FaultGeometrySpec` records `ref_normal_entries`, `up_entries` (0 absent, 2, 3); new inline `CheckFaultGeometryDimension(g, dim)` | the vectors and their defaults unchanged; the entries are read by the check only |
| `spatial/code/spatial_friction.cpp` | `toml_vec3` (used for `[fault_geometry]` ref_normal and up only) takes 2 or 3 entries, pads z = 0, returns the count | 3-entry vectors read, normalized and guarded as before; a 2-entry vector, which used to stop the load ("must be a 3-element array"), now loads and is refused by the check in the driver |
| `cycle/cycle_config.hpp` | `detail::Vec3` takes 2 or 3 entries (z = 0) and records "[section].key" and the count in `CycleConfig::deck_vectors` (`detail::Field` passes the record on); the BP5 `volume_axes` recorded as 3-entry; new `CycleConfig::CheckDimension(dim)` (every recorded vector, then `CheckFaultGeometryDimension`) | every value read as before; a 2-entry vector moves from a load-time refusal to the driver's check |
| `drivers/fd_cycle_driver.cpp` | `cfg.CheckDimension(mesh dim)` after `LoadScaledMesh` | passes for every 3-entry deck |
| `drivers/nitsche_dyn_driver.cpp` | `CheckFaultGeometryDimension(cfg.fault_geometry, mesh dim)` after the mesh read (the driver stays 3-D only until phase C) | passes for every 3-entry deck |

Not changed: `spatial_dyn_driver.cpp` reads the shared deck but not
`[fault_geometry]`; the tests that load committed 3D cycle decks
(`tests/cycle/*`) do not call `CheckDimension`.

### The meshes (Gmsh 4.15.2, `-format msh22`, deterministic: two runs give the same file)

| File | Made by | Tags | Triangles | Fault | sha256 (first 16) | Size |
|---|---|---|---|---|---|---|
| `bp3/mesh/bp3_d30_25m.msh` | `gmsh -2 bp3.geo -format msh22 -setnumber dip 30 -o bp3_d30_25m.msh` | 1 natural 242, 3 fault 1600, 5 Dirichlet 279 (238 creep) | 127843 | 40 km, 25 m | da99995f14a786c9 | 7.0 MB |
| `bp3/mesh/bp3_d60_25m.msh` | the same, dip 60 | 210, 1600, 255 (230 creep) | 212659 | 40 km, 25 m | fd08af92b3b48b78 | 11.7 MB |
| `bp3/mesh/bp3_d90_25m.msh` | the same, dip 90 | 210, 1600, 248 (228 creep) | 244958 | 40 km, 25 m | 41e033c40e1f74ce | 13.8 MB |
| `tpv205/mesh/tpv205_2d_100m.msh` | `gmsh -2 tpv205_2d.geo -format msh22 -o tpv205_2d_100m.msh` | 103 fault 300, 105 absorbing 136 | 40526 | 30 km, 100 m | b1b9467a67d7c53e | 2.1 MB |
| `tpv105/mesh/tpv105_2d_100m.msh` | `gmsh -2 tpv105_2d.geo -format msh22 -o tpv105_2d_100m.msh` | 3 fault 604, 5 absorbing 152 | 56652 | 140 km, 100 m for \|x\| <= 24 km | 184a4d1b0f3ff570 | 3.0 MB |

- `bp3/mesh/bp3.geo` is Tandem's `examples/tandem/2d/bp3.geo` (Tandem
  a26753d6, the file's last commit c0ad90f9) with one change: the default
  resolution Lf = Ls = 0.025 km, the 25 m the BP3 description suggests
  (Tandem's default 0.020).  Geometry, tags and size fields are Tandem's: km,
  y up, the fault from the trace (0, 0) to 40 km down-dip, the creep segment
  below it to the 400 km bottom (tag 5, as the sides), the free surface and
  the bottom natural (tag 1), Tandem's "bottom: free surface" variant.
- `tpv205/mesh/tpv205_2d.geo` is the 2D slice of the folder's 3D meshes
  (SeisSol's TPV5 geometry): the fault y = 0, |x| <= 15 km, embedded (its
  ends are the TPV5 barrier), split at the nucleation and stress patch edges
  (+-1.5, +-6, +-9 km); the box +-60 km as in 3D (45 km beyond the fault
  ends, 60 km from the fault); the 3D size field 0.05 d + (d / 2.5 km)^2 +
  lc_fault capped at 5 km; SeisSol's tags 103 / 105 as the 3D meshes.
- `tpv105/mesh/tpv105_2d.geo`: TPV105 is TPV103 (whole space) with thermal
  pressurization, and its 2D version is the 7.5 km depth line of the 3D
  fault (the seven 2D stations x = 0, +-9, +-15, +-18 km of
  uploadTPV105_2D.pdf).  The description makes the fault velocity
  strengthening outside |x| = 18 km and gives it no end, so the fault
  crosses the whole box (no fault tip inside; its ends lie on the absorbing
  sides).  Split at +-15, +-18 (the box function's transition) and +-24 km
  (the end of the 100 m part, 6 km past the last station); the box x in
  [-70, 70] km, y in [-50, 50] km, more than a P wave travels to a side and
  back from the outermost stations in the 15 s run; tags 3 / 5 as the 3D
  TPV104 meshes.  (The 3D TPV104 meshes end the fault at +-18 km, where the
  2D stations at +-18 km would sit on the tips.)

How it was checked (cloud container, 10-25).  test_2d_deck: 15 of 15 (the
3D first-event deck reads 22 vector keys with 3 entries and passes the 3D
check; its 2D form the same keys with 2, z = 0 everywhere, and passes the 2D
check; up [0, 1] with ref_normal [1, 0], and up [0, 0, -1], pass; the spatial
loader alone normalizes [3, 4] to (0.6, 0.8, 0) and keeps the defaults with
0 entries for absent keys), ten refusals (2-entry ref_normal, up and a
station in 3D; a 3-entry plane_origin, ref_normal, an off-axis 3-entry up,
no up and no ref_normal in 2D; 4-entry vectors in both readers), and both
drivers refuse a 2-entry ref_normal on a 3D mesh right after reading it.
test_2d_fault_layout: 127 of 127 on 1, 2 and 4 ranks, now with M1 on the
five committed meshes (Dimension 2, the triangle count and the tag census of
each file, counted independently by a Python read of the .msh; the fault's
length and bounding box; the census: one patch, 0 junction vertices, 2 tips,
the BP3 creep segment one patch with 1 seam vertex; every oriented normal
the expected one, to 6.4e-13 at dip 60 from the file's 16 digits; the plus
side; the slot of the deck's up: BP3 [0, 1] the dip slot, TPV [0, 0, 1] the
strike slot, the out-of-plane slot exact; the report equal to the 1-rank
one) and D2 (the fault DOF layout on the BP3 dip 60 mesh at p = 1, 2, 3, with
METIS's partition and with one that cuts along the fault so that all 1600
segments are rank-boundary faces: the table equals the 1-rank one bit for
bit).

3D results are unchanged.  A 3D deck goes through the new readers with 3
entries in every vector.  The parity set, rebuilt from clean with review 3's
tree and with this change: TPV205 and TPV104 on 2 ranks, the cycle QD and
dynamic regimes on 1 rank and the dynamic regime on 2 ranks byte for byte;
the cycle QD regime on 2 ranks (1.8e-9) and the BP5 quick check on 8 ranks
(4.4e-8) within MUMPS's run-to-run noise, as at every earlier review; the 15
unit tests, the field, plate loading and Okada tests and the
fully-dynamic-cycle tests (P3 on 2 ranks included) pass, with 1-rank outputs
byte for byte.

## Phase C, part 1: the Nitsche dynamic operator on triangles (C1)

The Nitsche dynamic path now runs on a 2D triangle mesh in the (x, y) plane
(design rev 3, section 5.5).  The volume vectors hold D = 2 components per
node (byNODES, strides D nd); the fault state stays padded: d, V, positions
and tractions are global three-vectors with z = 0 and the frame is 9 entries
per DOF.  The out-of-plane frame slot is an exact zero (Proposition 1), so
`CloseDOF`, the closures, Kelvin-Voigt and Prakash-Clifton are the 3D code.
The stiffness, the slip and Dirichlet loads, the traction, the fault-face
penalties and the fast fault-face operators now accept 2D; the quasi-dynamic
solve (`Solve`, `BuildQDStageOperators`, `BuildQDGreen`) stays phase E and
refuses 2D.  The drivers stay 3D only until C2, which switches their own
stable-constant tables (`drivers/nitsche_dyn_driver.cpp`, `cycle/dyn_phase.hpp`)
to `NitscheStableConstant` and runs LTS in 2D; the LTS edit below is not yet
exercised on a 2D mesh.

New file (branch-owned): `tests/2d_inplane/test_2d_nitsche.cpp` (make target
`test-2d-nitsche`, part of `2d-inplane-all`).

| Shared file | Change | 3D |
|---|---|---|
| `domain/elasticity_operator.hpp` | `Require3D_`'s message names the QD solve only; `SetFirstStepDebugConfig` refuses the first-step dumps in 2D | message text only |
| `domain/elasticity_operator_assembly.inl` | `AssembleStiffness` no longer refuses 2D; the IP integrators (interior and Dirichlet faces), the BlockILU block size and the slip and Dirichlet loads take `dim_` in place of 3 | `dim_ = 3`: the same calls |
| `domain/elasticity_operator_traction.inl` | `ComputeTraction`, `AssembleSlipOnlyRHS`, `AssembleDirichletOnlyRHS` accept 2D (`dim_`); `Solve` gets its own `Require3D_` (it was refused through `AssembleStiffness` before) | no-op checks |
| `domain/elasticity_operator_{fault_ops,penalty,qd,debug}.inl` | `BuildFaultFaceOps` and `GetFaultFacePenalties` accept 2D; the face-traction kernels, the neighbour-exchange row size and the shared-face vdofs use `dim_` | 3 as before |
| `integrator/dg_elasticity_ip_combined_integrator.hpp` | **bug fix**: the trace of grad u in `ComputeTractionAtQuadPoints` and `ComputeTractionAtQuadPointsDecomposed` added `grad(2,2)` also when `dim_ = 2`, reading past the 2 x 2 matrix; now `dim_` terms.  Before the fix the 2D normal traction from pure tangential slip was garbage (N5 aborted with sigma_n_eff = -2.2e9 Pa) | the same three terms |
| `integrator/dg_absorbing_bc_integrator.hpp` | the block dimension is `fe.GetDim()` | 3 |
| `domain/nitsche_mass.hpp` | `dim_` from the mesh, triangles accepted in 2D; component loops and damping blocks D nd x D nd; `Dimension()` | `dim_ = 3` |
| `domain/nitsche_material_sample.hpp` | triangles accepted in 2D | tetrahedra as before |
| `domain/nitsche_wave_operator.hpp` | `dim_`, `Dimension()`; damping-block loops over `dim_`; Kelvin-Voigt face centroids over the face's vertices (2 in 2D) and element centres padded with z = 0; the prestress projection evaluates the deck's stress at (x, y, 0), uses the `dim_`-entry normal and refuses a 2D fault point with sigma_xz or sigma_yz nonzero (no anti-plane unknown); the penalty-spread abort message prints the centroid of 2 or 3 vertices | the same sums in the same order |
| `solver/nitsche_cfl.hpp` | `NitscheElementLength` accepts triangles (h = 4 A / P, twice the inradius, as 6 V / A for tets); new `NitscheStableConstant(dim, p)`: the 3D column (the drivers' `kCeig`) and the measured 2D column | `ComputeMaxDt` unchanged for tets |
| `solver/nitsche_lts_run.hpp` | the LTS block builder and the cluster mass loop use the mesh dimension | 3 |
| `tests/2d_inplane/2d_inplane.mk` | `test-2d-nitsche` | - |

The 2D stable constant.  The 3D column does not carry over to triangles.
N2 measures C = 2 / (omega_max dt_unit) with dt_unit = h_min / (p^2 c_p) on
the 4 x 4 grid of 4-triangle squares, by power iteration and by bisection on
blow-up over 1000 steps; the 8 x 8 grid gives the same C_eig to 0.03%:

| p | C_eig | bisection | 2D column | 3D column | 2D / 3D |
|---|---|---|---|---|---|
| 1 | 0.215484 | 0.215484 | 0.215 | 0.113 | 1.91 |
| 2 | 0.365158 | 0.365158 | 0.365 | 0.212 | 1.72 |
| 3 | 0.464475 | 0.461572 | 0.461 | 0.282 | 1.64 |
| 4 | 0.526223 | 0.526223 | 0.526 | 0.332 | 1.59 |

On the committed unstructured meshes C_eig is higher (N2M): TPV205-2D 100 m
0.244, 0.414, 0.522 and TPV105-2D 100 m 0.256, 0.440, 0.562 at p = 1, 2, 3,
1.13 to 1.22 times the column, so the column is safe for them.

Where the 2D tests differ from the 3D ones they follow.  (1) N3 drift: the
3D slope bound (1e-10 E_0 per step) failed at p = 1 (1.18e-10) although
E_h was conserved to 8e-15: the 2D step is 1.9 times P3's and the dt^2
oscillation of E is larger, so a drift-free E fits that slope.  The test
now asks the drift over the run to stay below 1% of the oscillation bound
(measured 0.10%, 0.11%, 0.66%); E_h is still checked to 1e-10.  (2) N3
temporal order: at the production step the first error ratio was not yet
asymptotic (U 3.08, Udot 5.14 at p = 1); the refinement starts at half the
step (ratios 3.96 to 4.22).  (3) N4 at 45 degrees: the 3D test could only
print the oblique reflection; the 2D box is wide enough that the side
walls' incident-only data cannot reach the fit region, so R_PP and R_PS are
asserted within 0.002 of the analytic values at both resolutions.

How it was checked (cloud container, 10-25 and 10-26).  test_2d_nitsche on
1 rank, 74 of 74 (N1 to N5): N1 plane P and S waves, finest-pair L2 orders 2.10 / 2.05
(p = 1), 2.61 / 2.99 (p = 2; the P value close to its 2.5 floor), 3.98 /
3.89 (p = 3); N3 as above, Newmark loop = `CentralDifferenceSolver` bit for
bit, `ComputeMaxDt` = the serial 4 A / P reference bit for bit; N4 the
damped mass inverse to 8e-16, normal P |D| 7.8e-4 -> 5.2e-5 and S 7.1e-4 ->
7.2e-5 with resolution, at 45 degrees R_PP -0.01381 / -0.01380 against
-0.01385 and R_PS -0.07070 / -0.07070 against -0.07090; N5 at p = 1 (max
slip 7.29 m) and p = 2: the out-of-plane slots exactly 0, |sigma_eff -
120 MPa| <= 8e-7 Pa, the half turn 2.7e-13, the 90-degree rotated run (slip
in the dip slot) 3.7e-13 / 4.7e-12, the fast fault operators against the
CSR path 7.7e-13 / 5.4e-12, Kelvin-Voigt and Prakash-Clifton finite,
Prakash-Clifton with L = 1e-200 bit for bit the run without it.  On 2 and 4
ranks (`--np`, 42 of 42) with a partition that puts the two sides of the
fault on different ranks (every fault face a rank-boundary face): N3 bit
for bit, N5 against the 1-rank run 4.4e-13 / 6.0e-12 (2 ranks) and 5.8e-13
/ 4.0e-12 (4 ranks).

3D parity against review 4's tree (same container, both binary sets rebuilt
from clean; review 4's from a `git archive` export, whose files carry the
commit stamp "unknown"): TPV205 and TPV104 on 2 ranks, the cycle QD and
dynamic regimes on 1 rank byte for byte apart from the commit stamp; the
cycle QD regime on 2 ranks 1.3e-9 (review 4's binaries run twice: 2.7e-9),
the BP5 quick check on 8 ranks 6.0e-8 (twice: 1.1e-7); the cycle dynamic
regime on 2 ranks differs only in the t = 0 quasi-static velocity of the
off-fault stations, 3.7e-9 (twice: 2.2e-9; it needs the container to itself,
2 x 7 GB, and was killed by the memory limit once when run next to other
tests); the 15 unit tests pass with identical outputs.  This phase also
compares the Nitsche verification tests (19 runs of the Makefile's
`test-nitsche-*` targets): identical apart from timings, MPI process names
and, in P5 on 2 ranks, QD-solve diagnostics at MUMPS's run-to-run level.  Two
of their checks fail with both binary sets, identically, and also at the
branch point 3926e9c8: KV (e) at p = 2 (node-scaled L2 order 2.332, floor
2.5) and TPV102 inline T-21(a) (`Mult` against M^{-1}(F - Kx + Gd - Cv) bit
for bit: 6088 and 18422 entries differ); they are not touched here.  The
fully-dynamic-cycle tests pass; their 1-rank logs are identical apart from
the header and timings, their 2-rank logs differ at MUMPS's run-to-run level
(the P3 switch checkpoint 3.4e-11).  `2d-inplane-all` passes.

Known behaviour changes in 3D: none.

## Phase C, part 2: the Nitsche driver on 2D meshes, TPV205-2D (C2)

`seas_nitsche_dyn_driver` now runs a 2D triangle mesh in the (x, y) plane
(plane strain, mode II) with the operator of C1.  A 2D deck is the 3D deck
with 2-entry `[fault_geometry]` vectors (B2); `up = [0, 0, 1]` puts the
in-plane slip in the strike slot, so the station columns, the fault CSVs and
the conventions file are the 3D ones (slip2 = SCEC h-slip).  The step uses
the 2D column of `NitscheStableConstant` (C1) in the global step, the LTS
clusters and the logs; the cycle's dynamic phase reads the same function, so
the drivers no longer keep their own copies of the 3D column.  In 2D the
driver refuses, with a message, what it would otherwise misread: a station
or a nucleation centre off the plane (z != 0), `[[stress_region]]` z bounds
(every fault point has z = 0), and the fault ParaView output (`fault.vtkhdf`
holds triangles; the 2D fault is a polyline, phase F).

New files (branch-owned): `tpv205/configs/tpv205_2d_nitsche.toml` (gate G1),
`tpv205/stations/tpv205_2d_fault_stations.csv` (written by
`scripts/make_fault_stations.py tpv205_2d`), and
`tests/2d_inplane/compare_station_runs.py` (make target `test-2d-driver`,
part of `2d-inplane-all`).

| Shared file | Change | 3D |
|---|---|---|
| `drivers/nitsche_dyn_driver.cpp` | the mesh check takes a 2D mesh in the plane (triangles) or a 3D mesh (tets); `kCeig` removed, `NitscheStableConstant(dim, p)` in its five places; the LTS centroid check sums `SpaceDimension()` coordinates and the DG0 cluster field has the mesh dimension; the station copy average and the fault CSV dump use the edge midpoint and length of a 2D fault face; the four 2D refusals above; "tets" in the log becomes "triangles" in 2D; the no-HDF5 refusal of the fault ParaView output moved after the mesh read, so that a 2D deck gets the 2D message | the same constants (0.113, 0.212, 0.282, 0.332), sums and text; a deck asking for `fault.vtkhdf` from an MFEM without HDF5 is refused after the mesh read instead of before it |
| `cycle/dyn_phase.hpp` | `DynPhase::kCeig` removed; the CFL step and the LTS assignment use `NitscheStableConstant` of the mesh dimension; `StableConstant()` for the log | the same values |
| `drivers/fd_cycle_driver.cpp` | the dynamic-setup log prints `dyn.StableConstant()` | the same text |
| `scripts/make_fault_stations.py` | problem `tpv205_2d`: the seven on-fault stations of the TPV205-2D description, named as the 3D stations at 7.5 km depth, written into `tpv205/stations/` with z = 0 | the other station files unchanged (`--check`) |
| `tpv205/mesh/tpv205_2d.geo`, `tpv105/mesh/tpv105_2d.geo` | comment only: the Nitsche driver has no `--mesh-scale` | the 100 m mesh regenerates byte for byte |
| `tests/2d_inplane/2d_inplane.mk` | `test-2d-driver` | - |

`dynamic/precomputed_face_fluxes.hpp`, which the design's phase C row lists,
belongs to the velocity-stress path; the Nitsche driver does not use it, so
it is not touched.

How it was checked (cloud container, 10-26).  TPV205-2D, the committed deck,
12 s on 4 ranks: p = 1 on 200, 100 and 50 m fault meshes (2.01, 4.11 and
7.68 ms/step) and p = 2 on 200 and 100 m (4.63 and 10.32 ms/step); the
200 m and 50 m meshes from `tpv205_2d.geo` with `-setnumber lc_fault`,
not committed.  At every station the final slip of the three finest runs
(p = 1 50 m, p = 2 200 m and 100 m) agrees within 0.3%, the rupture arrival
(first |V| > 1 mm/s) within 0.04 s except at +12 km (0.13 s, a slow
precursor; at 1 m/s 0.023 s); peak slip rates still change by about 10%
between those runs.  The dip columns (slip1, V1, tau1) are exactly 0 at
every station and step of every run.  `--lts rate2` on the 100 m mesh
(4 clusters, ideal 2.21x): 1.93x faster at p = 1 and 1.94x at p = 2, the
stations within 3.0e-3 (V), 1.5e-3 (tau) and 2.7e-4 (slip) of the global
step, relative to each column's largest value.  1 rank against 4 ranks
(2 s): the same time columns, 8.2e-11 at most, the last printed digit;
`test-2d-driver` asserts it at 0.1 s with 1e-9 (bit for bit in its run),
runs the LTS deck and checks the four refusals.  The PyLith comparison
(gate G1) needs the PyLith file on the Mac and follows the push.

3D parity against review 5's tree (same container, the parity set of C1
rerun with this change): TPV205 and TPV104 on 2 ranks and the cycle QD and
dynamic regimes on 1 rank byte for byte; the cycle dynamic regime on 2 ranks
byte for byte (against C1's rerun of it alone; C1's run in the set had been
killed by the memory limit); the cycle QD regime on 2 ranks 1.9e-9 and the
BP5 quick check on 8 ranks 7.0e-8 (MUMPS run to run; C1: 2.7e-9 and
1.1e-7 between two runs of one binary); the 15 unit tests pass with
identical outputs.  TPV205 1000 m with `--lts rate2` on 2 ranks (1 s):
byte for byte.  The fully-dynamic-cycle tests pass; their 1-rank logs are
identical apart from paths and timings, their 2-rank logs differ at MUMPS's
run-to-run level (the P3 switch checkpoint 8.8e-12).  `2d-inplane-all`
passes, `test-2d-driver` included.

Known behaviour changes in 3D: none in a run; the no-HDF5 refusal comes
after the mesh read.

## Phase D, part 1: the thermal pressurization module (D1)

Thermal pressurization (TPV105, design rev 3 section 5.9, derivation 10.5)
as a module of its own, `fault/thermal_pressurization.hpp` (branch-owned):
per fault point a 1-D grid across the fault, temperature and pore pressure
by finite volumes with an exactly integrated Gaussian source, backward Euler
(default) or Crank-Nicolson, two tridiagonal solves per point per step.  It
is not yet called by any operator or driver: the hook into the fault step,
the deck keys and the station columns are part 2.

New files (branch-owned): `fault/thermal_pressurization.hpp`,
`tests/2d_inplane/test_thermal_pressurization.cpp` (make target
`test-thermal-pressurization`, part of `2d-inplane-all`).

No shared file changes; nothing in the 3D path includes the new header, so
the 3D parity set is not affected (not rerun for this part).

How it was checked (cloud container, 10-26).  test_thermal_pressurization
19 of 19 and its 7 refusals: the grid; zero heating and Lambda = 0 leave T
and p exactly unchanged; the heat content equals the heat put in to 2.1e-13;
the analytic formulas reduce to Rice's sqrt(t) law (h -> 0) and to
p = Lambda theta (alpha_hy -> 0); TPV105 parameters with a constant 30 MW/m^2
source and with the source switched off at 2 s, on a fine grid (dz0 0.1 mm,
growth 1.01, dt 0.1 ms, Crank-Nicolson), T and p at z = 0, h, 2h, 5h within
3.0e-5 of the fault value from 0.01 to 15 s; space error ratios 3.75 and 3.94
per halving; time ratios 2.00 (backward Euler) and 4.00 (Crank-Nicolson); at
the production settings (growth 1.05, z_max 0.4 m, the 2D 100 m steps) T and
p on the fault within 6.1e-4 with dz0 = h/10 and 3.4e-4 with dz0 = h/20;
three points equal three one-point runs bit for bit, and a GetState /
SetState restart continues a run bit for bit.

Known behaviour changes in 3D: none.

## Phase D, part 2: the normal-stress hook and TPV105-2D (D2)

Thermal pressurization (part 1's module) now runs inside a rupture: an
additive hook in the Nitsche operator (hook H-2D1 in
`document/fd_cycle_dev/HOOKS.md`) lets a modifier lower the normal stress
the friction law sees, and the driver builds the module behind it from a
`[friction.thermal_pressurization]` deck table.  At owned fault DOF j the
law sees sigma_eff - (p_j - p_ini), or sigma* - (p_j - p_ini) with
Prakash-Clifton, with p_j from the end of the previous step; after the
closures of every owned DOF the module is advanced once with the step's dt,
tau_j = |X_tan + gamma V| (the shear traction the law balanced, the station
traction) and |V_j^{n+1}| (design rev 3 section 5.9: a one-step lag).  The
bulk keeps its sigma_eff: the modifier changes the strength only.  With
local time stepping the fault cluster's step is the dt.

New files (branch-owned): `fault/fault_normal_stress_modifier.hpp` (the
interface), `fault/thermal_pressurization_modifier.hpp` (part 1's module
behind it), `tests/2d_inplane/test_normal_stress_hook.cpp` (make target
`test-normal-stress-hook`), `tpv105/configs/tpv105_2d_nitsche.toml` (gate
G2) and `tpv105/stations/tpv105_2d_fault_stations.csv` (written by
`scripts/make_fault_stations.py tpv105_2d`); make target `test-2d-driver-tp`.
Both targets are part of `2d-inplane-all`.

| Shared file | Change | Without the new inputs |
|---|---|---|
| `domain/nitsche_wave_operator.hpp` | additions only (57 lines): `SetNormalStressModifier` (checks the modifier's point count against the owned fault DOFs) and `NormalStressModifier()`; in `CloseDOF` the law's normal stress less `NormalStressReduction(j)`, refused with a message when that leaves [sigma_min, inf); in `AdvanceFault_` the step's tau and \|V\| per owned DOF and one `Advance` after the loop | the modifier pointer is null: no new branch is taken, every bit as before |
| `dynamic/fault_station_file.hpp` | `SetExtraColumns(labels, value)`: named columns after the state column, in the header and in every row | no extra columns: the same bytes |
| `drivers/nitsche_dyn_driver.cpp` | the `[friction.thermal_pressurization]` table (alpha_th, alpha_hy, rho_c, Lambda, h, T_ini, p_ini, dz0, growth, z_max, all required; scheme = "backward_euler" (default) or "crank_nicolson"; unknown keys refused); the module on the owned fault DOFs behind the hook and its log line; the station sigma_n column becomes the effective normal stress sigma_eff - (p - p_ini) (TPV105's "effective normal stress"), T and p follow the state column, averaged over the DG copies like the other columns; the conventions file names them; the fault CSV dump gains T and p; the step log and the last line print max T, max p and the smallest sigma_eff - (p - p_ini) | no table: the station copy average is called with no extra values and does the same arithmetic; the same output bytes |
| `scripts/make_fault_stations.py` | problem `tpv105_2d`: the seven 2D stations of TPV105 (x = 0, +-9, +-15, +-18 km on the y = 7.5 km line), named as the 3D stations, written into `tpv105/stations/` with z = 0 | the other station files unchanged (`--check`) |
| `tests/2d_inplane/compare_station_runs.py` | `--common`: compare the columns both runs have, by name, and list the others | without it the same checks and output |
| `tests/2d_inplane/2d_inplane.mk` | `test-normal-stress-hook`, `test-2d-driver-tp` | - |
| `document/fd_cycle_dev/HOOKS.md` | the H-2D1 row (not pinned: the cycle does not call it) | - |

How it was checked (cloud container, 10-26).  test_normal_stress_hook 36 of
36 and its 2 refusals: a modifier that takes off 0 Pa and thermal
pressurization with Lambda = 0 give the run without a hook bit for bit
(p = 1 and 2, Kelvin-Voigt and Prakash-Clifton on); the closure sees
sigma - R_j; one Advance per step with the operator's dt, tau and |V| as
the operator has them, tau = f (sigma - R_j) to 3.7e-16 while slipping;
TPV105 thermal pressurization on the test rupture keeps the heat balance
to 1.7e-15.  TPV105-2D (the committed deck and 100 m mesh, 15 s, 4 ranks):
at p = 2 the arrivals are within 0.008 s of MDSBI and SPECFEM2D (the
benchmark check of 10-25) at six of the seven stations (+18 km: 0.022 s)
and the slip at 12.5 s within 0.5 % from -15 to +9 km (+15 km -2.6 %,
+18 km -6.9 %, -18 km -1.7 %, where the two codes differ by 14 %); p = 1
arrives 0.05 s early, which p = 2 removes.  Without the table the rupture
does not leave the nucleation patch.  Lambda = 0 against no table, 15 s:
the 8 shared station columns bit for bit.  1 rank against 4 ranks (2 s):
8.3e-11, the last printed digit; `test-2d-driver-tp` asserts it at 1 s with
1e-9 and checks the three deck refusals.  Crank-Nicolson against backward
Euler: final slip within 0.19 %, T and p within 0.4 %.  `--lts rate2`
(4 clusters, ideal 1.97x): 1.52x faster, final slip within 0.73 %, the
stations within 1.3e-3 (slip), 3.6e-2 (V), 1.7e-2 (tau), 3.8e-3 (T) and
5.5e-3 (p) of the global step, relative to each column's largest value
(the same size as Crank-Nicolson against backward Euler: this rupture
amplifies small changes at the front).  Thermal pressurization adds 24 %
to the p = 1 step on 4 ranks (6.12 to 7.59 ms; 77 unknowns and two
tridiagonal solves per fault DOF on the ranks that hold the fault).

The bulk sigma_eff at the stations moves by up to 2.5 MPa within about
0.1 s of the rupture front and returns to 120 MPa; it does not shrink from
p = 1 to p = 2 (TPV205-2D's +-1 MPa does).  A planar fault in a uniform
medium has no such change; gate G2's n-stress comparison with MDSBI will
show how much it matters.  Not changed here.

3D parity against review 6's tree 3bd1bc23 (same container, all 31
binaries rebuilt from clean): TPV205 and TPV104 on 2 ranks, the cycle QD and
dynamic regimes on 1 rank and the 15 unit tests byte for byte (apart from
the commit stamp); the cycle QD regime on 2 ranks 9.9e-10 and the BP5 quick
check on 8 ranks 3.7e-8 (MUMPS run to run; C2: 1.9e-9 and 7.0e-8); the
cycle dynamic regime on 2 ranks byte for byte but one off-fault station
file at 2.2e-13, which is run to run: a second run of the same binary
matches 3bd1bc23 byte for byte.  TPV205 1000 m with `--lts rate2` on 2
ranks (1 s): byte for byte.  The fully-dynamic-cycle tests pass; their
1-rank logs are identical apart from timings, the 2-rank P3 switch
checkpoint 6.9e-12 (C2: 8.8e-12).  TPV205-2D at p = 1, 12 s on 4 ranks:
byte for byte with phase C's final run.  `2d-inplane-all` passes, the
three phase D targets included.

Known behaviour changes in 3D: none.

## Phase E, part 1: the quasi-dynamic solve on 2D meshes (E1)

The quasi-dynamic regime's static solve now runs on a 2D triangle mesh.
Phases B and C took every piece it uses to 2D (the 2D spaces, the fault
layout and frames, the stiffness, the slip and Dirichlet loads, the
traction); `ElasticityDomainOperator::Solve` itself still refused a 2D mesh
through `Require3D_`.  That call is removed.  The stage operators
(`BuildQDStageOperators`, `--qd-operators matrix`) and the Green's function
(`--qd-green`) are separate prebuilt paths and still refuse 2D, with a
message that names them; the default regime, the face loops through
`Solve`, is the one BP3-QD needs.

New file (branch-owned): `tests/2d_inplane/test_2d_dislocation.cpp` (make
target `test-2d-dislocation`, part of `2d-inplane-all`): a uniform dip-slip
patch (W = 2 km about 2 km depth, slip 1 m) imposed through the QD
operator's slip embedding on a sheared 2D triangle grid, against Okada
(1992) in the middle plane of a 2000 km long patch (plane strain to 7e-7 m
and 1e-6 of mu |b| / W, P0), for dips 30, 60 and 90 degrees:

- P1, thrust, p = 1 and 2, n = 4 and 8 cells along the patch: the surface
  displacement error e_u and the stress part of the fault traction fall
  with n at every dip and p (dip 60, p = 2: e_u 3.65e-3 to 4.33e-4, stress
  part 0.055 to 0.014, the whole traction 0.54 to 0.023); the strike slot
  of the slip and of the traction is exactly 0.
- P2, normal sense: the tractions are the thrust case's negated, DOF by
  DOF (6e-15 of the largest at most).
- P3, ref_normal flipped: the driving traction reverses, the normal
  traction stays (3e-15).
- P4, the whole problem rotated by 30 degrees in the plane: the fault-local
  tractions are unchanged (1e-14).

This fixes the slip convention BP3 needs: the kernel's slip is the jump of
the side the oriented normal points away from over the other, in the frame
(d, s) with d down-dip for up = [0, 1] (design 10.2, case B).

| Shared file | Change | Without the new inputs |
|---|---|---|
| `domain/elasticity_operator_traction.inl` | `Solve` no longer calls `Require3D_` | on a 3D mesh `Require3D_` passed: the same instructions run |
| `domain/elasticity_operator.hpp` | `Require3D_`'s comment and message: the stage operators and the Green's function are 3D only, the 2D regime runs on the face loops | message text only |

Test changes (branch-owned): `test_2d_fault_layout` refusal 5 now builds the
QD stage operators on a 2D mesh (it called `Solve`, which no longer
refuses); `test_2d_nitsche` N2 accepts a coded C(p) within the bisection's
1 % bracket below the measured value, because the measurement differs by
platform inside that bracket (p = 3: 0.461572 in the cloud container,
0.464472 on the Mac, coded 0.461); the check that C(p) is at most the
measurement is unchanged.

3D parity against the runs of review 8's binaries (b0ce3cb6's code; same
container, all 31 binaries rebuilt from clean): TPV205 and TPV104 on 2
ranks, the cycle QD and dynamic regimes on 1 rank and the 15 unit tests
byte for byte (apart from the commit stamp; the logs differ in the same
timing lines as at D2); the cycle QD regime on 2 ranks 2.5e-9 and the BP5
quick check on 8 ranks 1.1e-7 (MUMPS run to run; C1 measured 1.05e-7
between two runs of one binary); the cycle dynamic regime on 2 ranks byte
for byte but one off-fault station file at 2.2e-13, the run-to-run
difference D2 measured.  TPV205 1000 m with `--lts rate2` on 2 ranks
(1 s): byte for byte.  The fully-dynamic-cycle
tests pass; their 1-rank logs are identical apart from timings, and the
2-rank P3 switch checkpoint differs by 2.9e-11, the size of the test's
own second QD run in the same process (t 4.6e-11; D2: 6.9e-12 and
7.6e-12).  The first try of the 2-rank P3 test was killed by the
container's 14.3 GB memory limit (out of memory, not a failing check);
run again alone, it passes.  TPV205-2D at p = 1, 12 s on 4 ranks:
byte for byte with phase C's final run.  `2d-inplane-all` passes (956 PASSED lines, no FAILED), the
new target included.

Known behaviour changes in 3D: none.
