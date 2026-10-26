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
