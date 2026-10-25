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
