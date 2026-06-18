# Implementation Plan: SAFS LSW — Post-Processing Concept Activities

> Activities 6.1 PGA/PGV/PGD · 6.2 Seismic potency & moment · 6.3 Rupture velocity,
> slip rate & particle velocity · 6.4 Seismogram.
>
> Built on the EXISTING `safs_seisol_v2_2_0_LSW` case (linear slip-weakening, FL=16,
> candidate D, mu_d=0.12, d_c=1.8 m). This plan ADDS outputs + station files + activity
> docs; it does NOT change the physics and does NOT touch the existing `parameters.par`
> (PGV run) or `parameters_kinematics.par` (kinematics run) — both are kept verbatim.

## Overview

The intern activity is reorganized around FOUR post-processing CONCEPTS, each taught by
COMPUTING the quantity two or more independent ways and CROSS-CHECKING the routes against
each other and against SeisSol's own outputs. The goal is conceptual depth — understand
*why* a quantity is what it is and *why* independent routes agree (or disagree) — not just
"load the array and color the mesh".

The four activities share ONE simulation (one QuakeWorx job, one downloaded output set) so
every cross-check is internally consistent: the receiver at station S is the same run as the
surface PGV pixel at S, which is the same run as the fault slip that fed the moment. A single
combined deck `parameters_postprocess.par` emits everything and stays under the 1 GB download
limit.

### Why constant material first (kept from the group-activity design)
Run the CONSTANT medium first: `rho=2670, mu=lambda=3.2e10` => Vs=3462 m/s, Vp=5996 m/s
(Vp/Vs=sqrt(3), Poisson 0.25). Constant mu makes `M0 = mu * P0` EXACT (the harmonic mean of
two equal moduli is mu), and makes Vp/Vs the same everywhere so the seismogram arrival-time
checks (6.4) and supershear test (6.3) have a single known wave speed. THEN switch to CVM
(`safs_material_cvm.yaml`) and observe what heterogeneity changes. The constant case is the
"trust the pipeline" anchor; CVM is the "now do science" case.

### The cross-activity verification web (this is the pedagogy)
The four concepts are not independent silos — they close a loop, and seeing the loop close is
the learning objective:

```
   6.4 SEISMOGRAM  (v(t) at a point) --- peak over t --->  6.1 PGV/PGA/PGD
        |  integrate v -> u (compare to surface u)              ^  (map = field of
        |  differentiate v -> a                                 |   per-station peaks)
        |  first arrival t_P = d/Vp  -> recover Vp -------------+--> verifies MATERIAL
        v
   on-fault PICKPOINT slip-rate(t)
        |  peak = local PSR  -------------------------------->  6.3 PSR
        |  integral over t  = local final slip = ASl ------->  6.2 P0 / 6.3 slip
        v
   6.3 SLIP RATE [[v]] (jump in particle velocity across fault)
        vs particle velocity at an off-fault receiver  ----->  source vs radiated field
   6.3 Vr = 1/|grad RT|  vs  Vs  ---------------------------->  sub/supershear -> 6.1 PGV highs
```

Each activity below explicitly names which other activity its result must agree with.

## Established facts (verified in SeisSol source at `/Users/chunhuizhao/projects/SeisSol`)

### Friction / case (unchanged from the existing LSW deck)
- FL=16 linear slip-weakening; fault model `safs_fault.yaml` (mu_s=0.43, mu_d=0.12, d_c=1.8 m,
  deep barrier z in (-20000,-15000), Tnuc_s Gaussian 10 MPa R=6000 at the hypocenter).
- Hypocenter `(x,y,z) = (606971, 3707270, -4965.62)` m, UTM 11N (x=East, y=North, z=Up).
  Surface epicenter ~ `(606971, 3707270, 0)`.

### Free-surface output (`FreeSurfaceWriterExecutor.cpp:75`)
- Writes `v1,v2,v3,u1,u2,u3` (particle velocity AND displacement) per surface facet.
- `surfacevtkorder = -1` => legacy compact XDMF/HDF5 (geometry once, cell-centered, one growing
  `.h5`); `surfacevtkorder >= 0` is the multi-GB VTKHDF path — FORBIDDEN here.

### Off-fault receivers — `&Output` (`OutputParameters.cpp:140-168`) [NEW for this case]
- `ReceiverOutput = 1` (master on/off; default true).
- `RFileName = '<file>'` — ASCII coordinate file, one receiver per line, 3 whitespace-separated
  columns `x y z` in mesh units (meters), NO header (`ReceiverWriter.cpp:46-62`). Point-in-tet
  search places it; receivers may sit on the free surface or anywhere in the volume.
- `pickdt` — sampling interval [s], default `0.005`. THIS is the high-frequency lever for PGA.
- `ReceiverOutputInterval` — flush interval [s] (default 1e100 = flush at end; set ~10 s).
- `ReceiverComputeRotation`/`ReceiverComputeStrain` — default false; leave off.
- Each receiver file records, per row: `Time` then `s_xx,s_yy,s_zz,s_xy,s_yz,s_xz,v1,v2,v3`
  (`Kernels/Receiver.cpp`, `Datastructures.h:36-37`). **Receivers record VELOCITY + STRESS,
  NOT displacement** — to get displacement at a receiver you INTEGRATE v (contrast the surface
  output, which gives u directly; this contrast is a 6.1 teaching point).

### On-fault pickpoints — `&Pickpoint` (`OutputParameters.cpp:118-138`) [NEW for this case]
- DISTINCT from off-fault receivers; interpolated to FAULT elements (must lie on the fault).
- **REQUIRES `OutputPointType = 5`** in `&DynamicRupture` (`AtPickpointAndElementwise`). With
  `OutputPointType = 4` (Elementwise only) the `&Pickpoint` block is SILENTLY IGNORED and no
  pickpoint files are written (verified: `DRParameters.h:43-47`, `OutputManager.cpp:130-135` — the
  pickpoint writer runs only when `outputPointType == 3` or `== 5`). VERIFIED-BY-RUN 2026-06-17:
  a constant-material run with `OutputPointType = 4` produced fault/surface/energy/receivers but
  ZERO pickpoint files; the deck now uses `= 5`.
- `PPFileName = '<file>'` — same `x y z` ASCII format as receivers.
- `printtimeinterval` — output cadence in NUMBER OF TIMESTEPS (not seconds).
- `OutputMask` — 12-element selector, same groups as `&Elementwise`.
- Writes fault quantities (slip rate group 1, traction group 2, slip group 6, ...) as a time
  series per point — the on-fault analogue of a seismogram (a local source-time function).

### Fault `&Elementwise` (`DataTypes.h` VariableLabels)
- Groups: 1=SRs/SRd 2=Ts/Td/Pn 3=u_n 4=Mud/StV 5=Ts0/Td0/Pn0 6=Sls/Sld 7=Vr 8=ASl 9=PSR
  10=RT 11=DS 12=P_f/Tmp. `refinement=0` => 1 value/facet (smallest; native 500 m).
- `ASl` = integral(|slip rate| dt) (`FrictionSolverCommon.h:603`); `RT` = first t with
  |slip rate|>1e-3 m/s; `Vr = 1/|grad RT|` (degree N-1 polynomial gradient,
  `ReceiverBasedOutput.cpp:463-523`, `Vr=0` where `RT==0`); `PSR` = running max |slip rate|.

### Energy CSV (`EnergyOutput.cpp:711-748`)
- LONG/tidy format. Column layout is build-dependent: the QuakeWorx build (VERIFIED-BY-RUN
  2026-06-17) writes 3 columns `time,variable,measurement` (value is COLUMN 3), NOT the 4-column
  `time,variable,simulation_index,measurement` of multi-sim builds. Worksheet parsing must key on
  the `variable` column, not a fixed column index.
- Ground-truth rows: `variable == potency` (P0 = integral slip dA [m^3],
  `:314-319`) and `variable == seismic_moment` (M0 = potency*mu [N m], `:320`). Also
  `total_frictional_work`, `static_frictional_work`, `elastic_*`, `*_kinetic_energy`, `momentum*`.
- Homogeneous medium => `seismic_moment = 3.2e10 * potency` EXACTLY (constant mu = harmonic mean).
- SeisSol does NOT output PGV/PGA/PGD — those are post-processing products (6.1).

## Constraints
- **Per-job download < 1 GB** (QuakeWorx; no cluster shell — all post-processing is laptop-side
  in ParaView + light Python). The ONE combined run must clear 1 GB on a double build
  (Phase 1 budget proves ~0.77 GB).
- **One simulation for all four activities** — required for the cross-activity web to be exact.
  EndTime = 100 s default (captures the central/southern rupture; far-NW peak may be truncated —
  acceptable: the M0/P0 verification compares ParaView vs CSV AT THE SAME snapshot time, so it is
  exact whether or not the rupture has finished; plateau only matters for "final event size").
- **Mesh fixed**: `safs_mesh.puml.h5` (121,316 fault facets BC 3; 34,096 free-surface triangles
  BC 1). Unchanged.
- **Legacy compact output paths only**: `surfacevtkorder=-1`, fault/surface `refinement=0`.
- **No solver changes.** Decks + station files + docs only.
- **Plain-text ASCII math in all docs** (no LaTeX); PDFs via pandoc+xelatex need ASCII-only
  source (replace `-> => ~= <= x` etc.), per the no-LaTeX project rule.
- **Do not revert** the existing `parameters.par` / `parameters_kinematics.par`; this work is
  additive (new `parameters_postprocess.par`).

---

## Phase 1: Combined output deck + station files (modify the LSW case)

### Goal
A single ready-to-submit deck `parameters_postprocess.par` emits fault `&Elementwise` +
free-surface + energy CSV + off-fault receivers + on-fault pickpoints from one run, provably
under 1 GB; the two station files and the constant-material file exist in the case root.

### Files to Create
- `parameters_postprocess.par` — the combined-output deck (full content below).
- `make_stations.py` — station generator: reads the fault geometry from `safs_mesh.puml.h5` and
  writes `safs_receivers.dat` (off-fault, trace-relative) + `safs_pickpoints.dat` (on-fault,
  facet-snapped). The PRIMARY way to pick points; algorithm in Phase 5 "How to pick station points".
- `safs_receivers.dat` — off-fault / free-surface receiver coordinates (generated by
  `make_stations.py`, or the starter set below as a fallback).
- `safs_pickpoints.dat` — on-fault pickpoint coordinates (generated by `make_stations.py`, or the
  starter set below as a fallback).
- `safs_material.yaml` — `cp archive/safs_material.yaml safs_material.yaml` (constant medium;
  verify `rho 2670.0, mu 3.2e10, lambda 3.2e10`).

### Files to Reuse (unchanged)
- `safs_fault.yaml`, `safs_initial_stress.yaml`, `safs_mesh.puml.h5`,
  `safs_material_cvm.yaml` + `safs_material_cvm.nc`.

### Detailed Requirements

1. Create `parameters_postprocess.par` with EXACTLY this content:

```fortran
! ============================================================================
! SAFS v2_2_0_LSW — COMBINED POST-PROCESSING ACTIVITY deck (one run feeds all
! of 6.1 PGA/PGV/PGD, 6.2 potency/moment, 6.3 Vr/slip-rate/particle-velocity,
! 6.4 seismogram). Physics identical to parameters.par (FL=16 candidate D);
! only the OUTPUT set is expanded and off-fault receivers + on-fault pickpoints
! are added. Submit as ONE QuakeWorx job. Keep parameters.par / _kinematics.par.
!
! MATERIAL: defaults to the CONSTANT medium (M0 = 3.2e10*P0 exactly verifiable,
!   single Vp/Vs for the arrival-time and supershear checks). For the CVM case,
!   change the ONE line below to MaterialFileName = 'safs_material_cvm.yaml'.
! DOWNLOAD-SIZE budget: see POSTPROCESS_ACTIVITY_PLAN.md Phase 1 (~0.77 GB double
!   at EndTime=100 s; all paths legacy/compact, refinement=0). Under 1 GB.
! ============================================================================
&equations
MaterialFileName = 'safs_material.yaml'        ! CONSTANT first. CVM: 'safs_material_cvm.yaml'
Plasticity = 0
numflux = 'godunov'
numfluxnearfault = 'godunov'
/

&IniCondition
/

&DynamicRupture
FL = 16                                         ! linear slip-weakening (candidate D)
ModelFileName = 'safs_fault.yaml'
s_0 = 0.0
t_0 = 1.0                                        ! nucleation ramp [s]
XRef = 0.0
YRef = -1.0
ZRef = 0.0
refPointMethod = 1
OutputPointType = 5                              ! 5 = elementwise AND pickpoint (4 = elementwise ONLY;
                                                 ! the &Pickpoint block is ignored unless this is 3 or 5)
/

! --- Fault surface maps: slip rate (front), Vr, ASl (slip), PSR, RT ---
&Elementwise
printtimeinterval_sec = 2.0                      ! 2 s = good slip-rate front movie for 6.3
! 1=SRs/SRd 7=Vr 8=ASl 9=PSR 10=RT selected (Vr/ASl/PSR/RT are end-state fields):
OutputMask = 1 0 0 0 0 0 1 1 1 1 0 0
refinement_strategy = 2
refinement = 0                                   ! 1 value/facet. KEEP 0 for <1 GB.
/

! --- On-fault PICKPOINTS: local slip-rate / traction / slip time series (6.3, 6.4) ---
&Pickpoint
PPFileName = 'safs_pickpoints.dat'
printtimeinterval = 10                           ! every 10 timesteps (~10-50 ms effective)
! 1=SRs/SRd 2=Ts/Td/Pn 6=Sls/Sld selected (on-fault source-time series):
OutputMask = 1 1 0 0 0 1 0 0 0 0 0 0
/

&SourceType
/

&SpongeLayer
/

&MeshNml
MeshFile = 'safs_mesh.puml.h5'
meshgenerator = 'PUML'
/

&Discretization
CFL = 0.5
ClusteredLTS = 2
FixTimeStep = 0.1
/

&Output
OutputFile = 'output/safs'
Format = 10                                      ! no volume wavefield output (large mesh)
iOutputMask = 0 0 0 0 0 0 1 1 1
iPlasticityMask = 0 0 0 0 0 0 1
TimeInterval = 5.0                               ! moot (Format=10)
refinement = 1                                   ! VOLUME refinement; moot under Format=10

! --- Free surface (6.1 PGV/PGD map; 6.3 surface particle velocity): legacy compact ---
SurfaceOutput = 1
SurfaceOutputRefinement = 0
SurfaceOutputInterval = 0.5                       ! v1/v2/v3 (->PGV) + u1/u2/u3 (->PGD)
surfacevtkorder = -1                             ! MANDATORY -1 (compact). >=0 is multi-GB.

! --- Off-fault RECEIVERS (6.1 PGA via fine pickdt; 6.4 seismograms) ---
ReceiverOutput = 1
RFileName = 'safs_receivers.dat'
pickdt = 0.005                                   ! 0.005 s = 100 Hz Nyquist for PGA (6.1)
ReceiverOutputInterval = 10.0                    ! flush every 10 s

! --- Energy CSV (6.2 ground truth: potency, seismic_moment) ---
EnergyOutput = 1
EnergyTerminalOutput = 1
EnergyOutputInterval = 0.25

Checkpoint = 0
/

&AbortCriteria
EndTime = 100.0                                  ! see Constraints: verification is exact at any
                                                 ! shared snapshot; extend to 150 s (with the
                                                 ! size levers in Phase 1) for a plateaued M0.
/

&Analysis
/

&Debugging
/
```

2. Create `safs_receivers.dat` (off-fault / free-surface receivers; z=0 = free surface).
   STARTER set near the epicenter; off-fault receivers need NOT lie on any feature (the
   point-in-tet search finds the element). The WORKSHEET 6.4 step documents harvesting
   trace-snapped along-strike stations from the run-1 surface output. Content:

```
# SAFS off-fault receivers (free surface, z=0). x y z, meters, UTM 11N. No header needed by
# SeisSol (lines starting with non-numeric are skipped as whitespace? -> keep coords only
# in the shipped file; this banner is for humans, STRIP before submitting if SeisSol errors).
# Epicenter
606971.0   3707270.0   0.0
# Fault-normal-ish cross through the epicenter (map offsets; snap/relabel after run 1):
611971.0   3707270.0   0.0
601971.0   3707270.0   0.0
606971.0   3712270.0   0.0
606971.0   3702270.0   0.0
616971.0   3707270.0   0.0
606971.0   3717270.0   0.0
# Two along-strike stations (NW); VERIFY they are off the fault, relabel by along-strike km:
590000.0   3725000.0   0.0
575000.0   3742000.0   0.0
```

   NOTE for the implementer: confirm whether the QuakeWorx SeisSol receiver parser tolerates
   `#` comment lines (the source skips all-whitespace lines only). If it errors, ship a
   coordinates-only file (no `#` banner) and keep the annotated version as
   `safs_receivers.annotated.txt`.

3. Create `safs_pickpoints.dat` (ON-fault; must lie on the fault surface). STARTER = the
   hypocenter plus two near-vertical neighbors (the SAF is steeply dipping, so same x,y at
   shallower/deeper z is approximately on-fault); WORKSHEET 6.4 documents harvesting exact
   along-strike/down-dip on-fault points from the run-1 fault geometry. Content:

```
606971.0   3707270.0   -4965.62
606971.0   3707270.0   -2000.0
606971.0   3707270.0   -8000.0
```

   NOTE: if SeisSol drops a pickpoint as "not on fault" (search tolerance), harvest replacements
   from the run-1 `*-fault.xdmf` facet centroids (procedure in WORKSHEET 6.4) and resubmit.

4. `cp archive/safs_material.yaml safs_material.yaml`.

5. Create `make_stations.py` (full algorithm in Phase 5 "How to pick station points", Route 1).
   I/O contract: INPUT `safs_mesh.puml.h5` (PUML; select the fault boundary tag = 3) + the
   hardcoded epicenter `(606971, 3707270)` and target distance lists (along-strike +-10/+-20/+40 km;
   fault-normal +-2/+-5/+-10/+-20 km; pickpoint depths 2/5/8/12 km). OUTPUT `safs_receivers.dat`
   (off-fault, z=0) and `safs_pickpoints.dat` (on-fault, snapped to nearest fault-facet centroid).
   Dependencies: `h5py`, `numpy`, `scipy` (`cKDTree` for nearest-facet snapping). Run once before
   submitting; commit the generated `.dat` files. The starter sets in requirements 2-3 are the
   no-scripting fallback (Route 2).

### Download-size budget (EndTime=100 s, double build)

```
Fault &Elementwise (legacy, geometry once, refinement=0, 6 comps SRs,SRd,Vr,ASl,PSR,RT):
   Nf=121,316; per snap 6*121316*8 = 5.82 MB; snaps = 100/2 + 1 = 51  -> 297 MB + 10 MB geom
Free surface (legacy compact, refinement=0, 6 fields v1..u3, geometry once):
   Ns=34,096; per snap 6*34096*8 = 1.64 MB; snaps = 100/0.5 + 1 = 201 -> 329 MB + 1 MB geom
Off-fault receivers (ASCII, 9 chan + time, pickdt 0.005 -> 20,000 rows):
   ~220 B/row * 20,000 = 4.4 MB/receiver * 9 receivers                -> ~40 MB
On-fault pickpoints (ASCII, ~6 chan, every 10 steps ~ 2e4 rows):
   ~13 MB total for 3 points (scales with point count)                -> ~15 MB
Energy CSV (long, 0.25 s * ~13 vars)                                  -> <1 MB
----------------------------------------------------------------------------------------
TOTAL (double) ~= 0.69-0.77 GB   (single build ~= 0.4 GB)             <<  1 GB
```

Levers (document in WORKSHEET): EndTime 100->150 needs `SurfaceOutputInterval 0.5->1.0` AND
`printtimeinterval_sec 2->5` to stay < 1 GB double; tighter quota -> drop receivers/pickpoints
to the few stations a worksheet actually uses; NEVER set `surfacevtkorder>=0` or any
`refinement=1` on a double build.

### Edge Cases to Handle
- Receiver/pickpoint `.dat` comment-line tolerance (requirement 2/3 NOTE) — ship coords-only if
  the parser is strict.
- A pickpoint reported off-fault is dropped silently — verify each appears in the run log; harvest
  replacements from fault geometry if missing.
- ASAGI required for the CVM material line (the QuakeWorx app has it — Kaikoura case uses ASAGI).

### Acceptance Criteria
- [ ] `parameters_postprocess.par` parses; log shows FL=16, constant (or CVM) material, and lists
      ReceiverOutput + Pickpoint enabled with the expected station counts.
- [ ] Output contains `*-fault.xdmf` (SRs,SRd,Vr,ASl,PSR,RT), `*-surface.xdmf` (v1..u3),
      `*-energy.csv` (potency, seismic_moment), receiver files (`*-receiver-*.dat`), pickpoint files.
- [ ] `du -sh output/` for a finished 100 s run is 0.4-0.8 GB.

### Dependencies
- Depends on: nothing. Required by: Phases 2-6.

---

## Phase 2: Activity 6.1 — PGA, PGV, PGD (the ground-motion ladder)

### Goal
Interns produce the PGA/PGV/PGD maps with SeisSol's OFFICIAL ground-motion post-processing script,
then reproduce PGV/PGD in ParaView and CONFRONT the two. The differences between the official
result and the ParaView result are the lesson — they are definition differences, not bugs.

### Physics (plain text)
The three peak ground-motion parameters are the three rungs of one ladder linked by time
derivatives:
```
acceleration a(t) = d v / d t        PGA = max_t |a|
velocity     v(t) = (SeisSol output) PGV = max_t |v|
displacement d(t) = integral v dt    PGD = max_t |d|   (or take u directly)
```
Differentiation amplifies high frequencies (PGA = most sampling-sensitive); integration smooths
(PGD = least sensitive, but accumulates baseline drift).

### PRIMARY tool: SeisSol's official GME script
`/Users/chunhuizhao/projects/SeisSol/postprocessing/science/GroundMotionParametersMaps/`
`ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py`. Run it on the downloaded
`safs-surface.xdmf`; it writes a ParaView-openable map of PGA/PGV/PGD (+ SA(T)/CAV).

**Conda setup + run (laptop, copy-paste). The version pins matter** — the script's
`requirements.txt` pins `numpy==1.22.0`, and the deprecated `gmpe-smtk` uses `np.float`/`np.int`
aliases that numpy >= 1.24 REMOVED, so the env must use Python 3.10 + numpy < 1.24 or the run dies
on import:

```bash
# ===== ONE-TIME setup =====
# 1. Dedicated env. Compiled deps from conda-forge (no MPI/HDF5 compiler hassle);
#    --noMPI means mpi4py is optional, but conda-forge ships it prebuilt so include it.
conda create -n seissol-gme -c conda-forge -y \
    "python=3.10" "numpy=1.22" "scipy=1.10" "h5py" "lxml" "mpi4py"
conda activate seissol-gme

# 2. SeisSol's small pure-python readers/writers (not on conda-forge -> pip)
pip install seissolxdmf==0.0.10 seissolxdmfwriter==0.4.0

# 3. gmpe-smtk (GMRotD50 + response-spectrum routines). Clone it INTO the script's own
#    folder: the script does sys.path.append("<script_dir>/gmpe-smtk") and finds it there,
#    so NO PYTHONPATH export is needed on each run. The deprecated commit needs no openquake.
cd /Users/chunhuizhao/projects/SeisSol/postprocessing/science/GroundMotionParametersMaps
git clone https://github.com/GEMScienceTools/gmpe-smtk
( cd gmpe-smtk && git checkout 4f008173e89f6e4ba4450fb43e95ffdf51a7c2ba^ )

# ===== RUN (every time) =====
conda activate seissol-gme
cd <the folder that holds the downloaded QuakeWorx `output/`>
GME=/Users/chunhuizhao/projects/SeisSol/postprocessing/science/GroundMotionParametersMaps
python3 "$GME/ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py" \
    --noMPI --MP 4 output/safs-surface.xdmf
#   --MP N   : N worker processes; N must be <= your CPU core count (the script asserts it)
#   --lowpass 1.0          : zero-phase 1 Hz low-pass on velocity BEFORE computing (try with/without)
#   --periods 0.3 1.0 3.0  : response-spectrum periods SA(T) [s]
#   --CAV                  : also compute cumulative absolute velocity
# OUTPUT: safs-GME-surface.xdmf (+ .h5), fields PGA, PGV, PGD [, SA..., CAV] -> open in ParaView.
```

If the import still fails with a `numpy has no attribute 'float'` error, the env picked up a newer
numpy — force it: `pip install "numpy<1.24"` inside `seissol-gme`.

### SeisSol's approach, in depth (what the official script computes, and why)

The script turns the raw surface VELOCITY history into the engineering-standard ground-motion
parameters. Verified against the source, step by step:

1. **Horizontal components only.** It reads `v1` (East) and `v2` (North) and DROPS the vertical
   `v3`. Ground-motion engineering characterizes HORIZONTAL shaking because that is what loads and
   damages structures; GMPEs (ground-motion prediction equations) are built on horizontal measures.
2. **(Optional) zero-phase low-pass.** `--lowpass fc` applies a 2-pass, order-2 Butterworth filter
   (`scipy.signal.filtfilt`). "Zero-phase" = run forward then backward so there is NO time shift —
   peaks stay at their true times. Purpose: strip numerical high-frequency content above what the
   mesh can physically carry (~1-2 Hz on a 500 m mesh) before taking peaks. PGA is most affected
   (high-frequency); PGD least.
3. **The acceleration -> velocity -> displacement ladder.** `a = np.gradient(v, dt)` (central
   finite difference of velocity); then it RE-INTEGRATES with `cumulative_trapezoid`:
   `velocity = integral a dt`, `displacement = integral v dt`. One code path yields all three (the
   re-integrated velocity equals the input velocity — self-consistent by construction). Peaks:
   `PGA = max|a|`, `PGV = max|v|`, `PGD = max|d|`.
4. **Orientation-independent reduction (GMRotD50) — the crux of "SeisSol's approach".** A station
   records two horizontals in whatever frame the grid happens to use, so the naive "peak of each
   component" DEPENDS on that arbitrary orientation (a station rotated 45 deg reports different
   numbers). GMRotD50 removes the dependence:
   - rotate the `(v1, v2)` pair through many azimuth angles `theta` in `[0, 90]` deg;
   - at each `theta` take the GEOMETRIC MEAN of the two rotated peaks,
     `sqrt(peak_x(theta) * peak_y(theta))`;
   - take the 50th PERCENTILE (median) of that set over all `theta`.
   The result is ONE orientation-independent number per station — the measure used by the NGA/PEER
   databases and GMPEs, so SeisSol's PGV/PGA compare directly to published attenuation curves.
   `gmrotdpp` is the fast per-period estimate (default); `--ipp` selects `gmrotipp`, the fully
   orientation-independent-across-periods version (~30x slower, more rigorous for SA).
5. **Response spectra SA(T).** For each period `T`, it computes the peak response of a 5%-damped
   single-degree-of-freedom oscillator of natural period `T` driven by the ground acceleration
   (Nigam-Jennings integration). SA(T) over many `T` = the response spectrum, the central design
   tool in earthquake engineering (how a building of period `T` responds). Default ~0.1-5 s.
   ParaView cannot make this — it needs an ODE solve per period per station.
6. **CAV (`--CAV`).** Cumulative absolute velocity `= integral |a| dt`, a cumulative
   shaking-intensity / damage measure.
7. **Units & output.** Internally engineering units (cm/s/s) for the SA step — mind this against
   ParaView's SI. It writes `safs-GME-surface.xdmf` (+ `.h5`) carrying per-surface-cell fields
   `PGA, PGV, PGD, SA(...)[ , CAV]` in SeisSol's own XDMF format — so the RESULT is itself a
   ParaView map; you visualize it like any SeisSol surface field, you just did not compute it
   cell-by-cell yourself.

Net: ParaView's `max_t sqrt(v1^2+v2^2+v3^2)` is the instantaneous 3D speed peak — trivial to make
and perfect for teaching the temporal-max operation, but it is 3D (not horizontal), frame-dependent
(not GMRotD50), unfiltered, and carries no SA/CAV. The official script is the engineering-standard
product; the gap between the two is exactly what the items below enumerate.

### COMPARISON: the same quantities in ParaView (on the SAME `safs-surface.xdmf`)
- PGV: `Calculator vmag=sqrt(v1^2+v2^2+v3^2)` -> `Temporal Statistics (Maximum)`.
- PGD: `Calculator umag=sqrt(u1^2+u2^2+u3^2)` -> `Temporal Statistics (Maximum)` (direct u).
- PGA: `Temporal Derivative` of (v1,v2,v3) -> `Calculator amag` -> `Temporal Statistics (Maximum)`.

### Differences to HIGHLIGHT (official vs ParaView) — the heart of this activity
1. **Components**: official uses the two HORIZONTALS ONLY (v1,v2; vertical dropped — engineering
   convention). ParaView `sqrt(v1^2+v2^2+v3^2)` includes the VERTICAL. => ParaView is generally
   larger; different DEFINITION, not a different answer.
2. **Orientation measure**: official = GMRotD50 (orientation-INDEPENDENT median over rotation
   azimuths) — the GMPE/engineering standard. ParaView (horizontal) = the as-recorded RESULTANT
   peak `sqrt(v1^2+v2^2)`, which is frame-dependent and always >= GMRotD50. Quantify the gap.
3. **PGD route**: official RE-INTEGRATES velocity to displacement (it ignores the output u).
   ParaView uses the directly-output u. Compare: agreement validates `integral v dt = u`; a gap
   exposes integration baseline drift (why direct u is the cleaner PGD).
4. **PGA & sampling**: official differentiates the SURFACE velocity, but at
   `SurfaceOutputInterval = 0.5 s` the Nyquist is 1 Hz, so BOTH the official and the ParaView
   surface PGA are low-frequency-limited and NOT a true PGA. Finer surface sampling costs download
   size; the off-fault RECEIVERS (`pickdt = 0.005 s`) are the high-frequency reference (used in 6.4)
   to show how far below true PGA the 0.5 s map sits.
5. **Filtering**: official offers `--lowpass`; ParaView is raw. Run with and without `--lowpass 1.0`
   and watch PGA move most, PGD least (the ladder again).
6. **Extra products**: official also yields SA(T) response spectra and CAV — ParaView does not.

### Cross-verification (what must agree, and why the rest legitimately differs)
- PGV(official, horizontal GMRotD50) vs PGV(ParaView, horizontal resultant `sqrt(v1^2+v2^2)`): same
  spatial PATTERN; ParaView amplitude >= official by the GMRotD50-vs-resultant factor — state it.
- PGD(official, integrated v) vs PGD(ParaView, direct u, horizontal): agree within integration drift.
- A single station's official PGV == the peak of that station's seismogram horizontal (link to 6.4).

### Deepen-understanding probe
"Your ParaView PGV map is 10-20% higher than the official one everywhere — is your pipeline wrong?"
(No: 3D-vs-horizontal + resultant-vs-GMRotD50; same physics, different definition.) "Which of
PGA/PGV/PGD can you trust from the 0.5 s surface map, and which needs the receivers?" (PGD yes,
PGV ok, PGA no — the sampling ceiling.)

### Acceptance Criteria
- [ ] `safs-GME-surface.xdmf` produced by the official script and opened in ParaView (PGA/PGV/PGD).
- [ ] ParaView PGV and PGD maps produced and compared to the official maps; the numeric gap is
      explained by the component / orientation / route differences (a small table).
- [ ] The PGA sampling ceiling demonstrated (official/ParaView surface 0.5 s vs receiver 0.005 s).

### Dependencies
- Depends on: Phase 1 (surface output + receivers) + the official-script setup. Cross-links: 6.4.

---

## Phase 3: Activity 6.2 — Seismic potency & seismic moment

### Goal
Compute P0 and M0 independent ways, verify against the energy CSV, and decompose the source into
ruptured area x mean slip — establishing what "potency" and "moment" physically are.

### Physics (plain text)
```
potency  P0 = integral over fault of (final slip) dA            [m^3]
moment   M0 = integral over fault of (mu * slip) dA             [N m]
            = mu * P0   ONLY IF mu is constant (constant case!)
mean slip Dbar = P0 / A_rup ;  A_rup = ruptured area
magnitude Mw = (2/3) log10(M0) - 6.07     (M0 in N m)
```

### What these quantities are (read before computing)
- **Potency `P0` [m^3]** = the final slip integrated over the fault, `P0 = integral D(x) dA`. It is
  the total "amount of faulting" — slip summed over area — and is INDEPENDENT of how stiff the rock
  is (slip [m] times area [m^2] => m^3).
- **Ruptured area `A_rup` [m^2]** = the area of fault that actually slipped, i.e. where the final
  slip `ASl > eps` for a small threshold `eps` (excludes numerically creeping cells).
- **Mean slip `Dbar` [m]** = the average slip over the ruptured patch = `P0 / A_rup`. The typical
  offset across the fault.
- **Seismic moment `M0` [N m]** = `integral mu*slip dA = mu*P0` (constant mu) = the textbook
  `M0 = mu * Dbar * A_rup`. `Mw = (2/3) log10(M0) - 6.07`.
- The identity `P0 = Dbar * A_rup` is just "mean = integral / area"; computing both sides
  independently is a consistency check AND it DECOMPOSES the source — a large-area/small-slip event
  and a small-area/large-slip event can share the same `P0` but are very different earthquakes.

### Inputs
- Fault `*-fault.xdmf` `ASl` (final-snapshot accumulated slip).
- Energy CSV `potency`, `seismic_moment` (use the last time for the final scalars).

### Compute & cross-verify (constant material = exact anchor)
1. **P0 two independent ways, plus the area / mean-slip decomposition**:
   - (A) `Integrate Variables` on the final fault snapshot's `ASl`. On CELL data ParaView returns
     the area-weighted sum `sum(ASl_cell * area_cell) = integral ASl dA = P0`, AND a total `Area`.
   - (B) decomposition: `Threshold ASl > eps` (e.g. `eps = 0.05` m) to keep only ruptured cells,
     then `Integrate Variables` on that patch -> its `Area` output = `A_rup`, its integrated `ASl`
     = `P0` (unruptured cells add ~0, so it matches A). Then `Dbar = P0 / A_rup` and verify
     `Dbar * A_rup == P0`. Vary `eps` (0.01, 0.05, 0.1 m) and watch `A_rup` / `Dbar` shift while
     `P0` barely moves — the rupture AREA is threshold-defined (teaching point), the POTENCY is not.
   - (C) ground truth: energy CSV `variable == potency`, last time. A == C within ~1% (quadrature).
2. **M0 and Mw**: `M0 = 3.2e10 * P0 = 3.2e10 * Dbar * A_rup` (constant medium) vs CSV
   `variable == seismic_moment` (last time); report `Mw`. EXACT band (constant mu = harmonic mean).

### Deepen-understanding probe
"P0 came out identical from `Integrate ASl` and from `Dbar*A` — why must it, and what is `Dbar*A`
telling you that the raw integral hides?" (Potency = average slip x area; decouples 'how much' from
'how big'.) "On the CVM run `M0 = 3.2e10*P0` fails — predict whether ParaView-M0 will be above or
below the constant estimate, given CVM mu(z)." (Depends on where slip concentrates vs mu(z).)

### CVM extension (after constant is trusted)
On CVM, `M0 != mu_const*P0`. Options: (a) report P0 (still exact via `Integrate ASl`) and the
SeisSol `seismic_moment`, explain the mu-weighting inequality; (b) ADVANCED — `Resample With
Dataset` mu from `safs_material_cvm.nc` onto the fault, `Calculator mu*ASl`, `Integrate` ->
approximate M0. Do NOT reuse the `3.2e10` factor on CVM (mu varies with depth).

### Acceptance Criteria
- [ ] P0 from `Integrate ASl`, from `Dbar * A_rup`, and from CSV agree within ~1% (constant).
- [ ] `A_rup` and `Dbar` reported with their `eps` sensitivity; `P0` shown insensitive to `eps`.
- [ ] `M0 = 3.2e10 * P0` matches CSV `seismic_moment`; `Mw` reported.

### Dependencies
- Depends on: Phase 1 (fault + energy). Cross-links: 6.3 (slip), 6.4 (local slip = ASl).

---

## Phase 4: Activity 6.3 — Rupture velocity, slip rate, particle velocity

### Goal
Distinguish three quantities that all have velocity units but mean different things: the rupture
FRONT speed (kinematic), the SLIP rate (a jump in particle velocity across the fault), and the
PARTICLE velocity (the radiated ground motion). Compute Vr three ways.

### Physics (plain text)
```
Rupture velocity Vr = 1 / |grad RT|    (speed of the RT=const front; PURELY GEOMETRIC -
                                        the slowness of the rupture-time contour map, NOT a
                                        material wave speed). Classify: Vr<Vs sub-shear,
                                        Vr>Vs supershear (Vs=3462 m/s constant).
Slip rate   = [[v]] = v(+) - v(-)   = the DISCONTINUITY in particle velocity across the fault.
              For symmetric bilateral slip each face moves +-|slip rate|/2.
Particle velocity v(t) = actual ground motion (what a receiver/surface records). On the fault it
              equals the slip-rate jump; off the fault it radiates and decays with distance.
```
"`dRT/dx` physical meaning": `grad RT` is the SLOWNESS VECTOR of the front (units s/m); its
magnitude is 1/Vr and it points along the propagation direction. `dRT/dx` alone = the front's
inverse apparent speed along x (= 1/Vr_app,x = cos(theta)/Vr): how many seconds the front takes to
advance one metre in x. So yes it is physical — a directional slowness component.

### Inputs
- Fault `Vr`, `RT`, `PSR`, `SRs,SRd`; on-fault pickpoint slip-rate; off-fault receiver v;
  surface v.

### Compute & cross-verify
1. **Vr three ways**: (A) SeisSol `Vr` field; (B) ParaView `Threshold RT>0` -> `Cell Data to
   Point Data` -> `Gradient` of RT -> `Calculator 1/mag(gradient)`; (C) MANUAL contour-spacing:
   pick two points on the front, read `RT1,RT2` and the distance `ds` between them,
   `Vr = ds/(RT2-RT1)`. A==B==C in the interior (edge/kink outliers expected); C makes "1/|grad
   RT|" concrete.
2. **Slip rate / PSR**: `|SR| = sqrt(SRs^2+SRd^2)`; `Temporal Statistics (Maximum)` -> compare to
   SeisSol `PSR`. ParaView <= SeisSol because the 2 s snapshot undersamples the pulse (Nyquist
   again — link to 6.1).
3. **Slip rate vs particle velocity**: overlay the on-fault PICKPOINT slip-rate(t) against the
   nearest OFF-fault receiver |v|(t). Same units, different roles: the on-fault trace is the SOURCE
   discontinuity; the off-fault trace is the RADIATED field — note the rupture-arrival time offset
   and the amplitude decay with distance.

### Deepen-understanding probe
"Vr came out 4500 m/s > Vs=3462 in one segment — is the fault breaking faster than its own shear
waves carry information, and how is that possible?" (Supershear; the front is a geometric locus, not
a signal — energy still travels at Vs/Vp; Mach cone forms.) "Slip rate peaked at 6 m/s on the fault
but the receiver 5 km away peaked at 0.4 m/s — which is the earthquake's 'speed' and which is the
ground's, and why are they not equal?" (Source jump vs radiated field via the Green's function.)

### Cross-links
- Supershear segments (Vr>Vs) should sit beneath PGV highs (Mach-cone enhancement) -> overlay the
  6.1 PGV map on the fault Vr (the original correlation question).
- PSR == peak of the pickpoint slip-rate seismogram (6.4); the pickpoint slip-rate integrated over
  time == local final slip `ASl` (6.2 P0 / 6.4).

### Acceptance Criteria
- [ ] Vr from SeisSol, from `1/|grad RT|`, and from manual contour spacing agree (interior median).
- [ ] PSR (ParaView temporal max) vs SeisSol PSR reported with the undersampling explanation.
- [ ] On-fault slip-rate vs off-fault particle-velocity overlaid with arrival-time + decay noted.

### Dependencies
- Depends on: Phase 1 (fault + pickpoints + receivers + surface). Cross-links: 6.1, 6.2, 6.4.

---

## Phase 5: Activity 6.4 — Seismogram

### Goal
Treat the seismogram as the fundamental observable from which the other three concepts derive;
read it, transform it, recover the material from it, and close the verification web.

### Physics (plain text)
A seismogram is the time series of ground motion at a point. Off-fault receiver = `(v, stress)(t)`
(integrate for displacement, differentiate for acceleration). On-fault pickpoint = `(slip rate,
traction, slip)(t)` = a local source-time function. First arrivals encode wave speed:
`t_P = d/Vp`, `t_S = d/Vs` for a receiver at distance `d` from the source.

### Inputs
- Off-fault receiver files `*-receiver-*.dat`; on-fault pickpoint files; surface `u` (for the
  integration cross-check); fault `ASl` (for the pickpoint-slip cross-check).

### Compute & cross-verify (light Python: read ASCII, plot, `numpy`/`scipy`)
1. **Read & plot** one receiver's `v1,v2,v3` vs `Time`. Identify P, S, and surface-wave arrivals.
2. **Transform ladder (closes 6.1)**: `d = cumulative_trapezoid(v, dt)` -> compare `max|d|` to the
   surface `u` map pixel at that station; `a = gradient(v, dt)` -> `PGA`. The station's
   (PGA,PGV,PGD) = the per-pixel values of the 6.1 maps. State the loop closed.
3. **Recover material (verifies the model)**: for a receiver at known map distance `d` from the
   epicenter, read the first-arrival time `t_P`; `Vp_est = d/t_P`; compare to constant Vp=5996 m/s
   (and `Vs` from the S arrival vs 3462 m/s). On CVM expect a path-averaged speed.
4. **On-fault pickpoint (closes 6.2/6.3)**: the pickpoint slip-rate(t) integrated over time = local
   final slip = the `ASl` value at that fault point (links to 6.2 P0 and 6.3); its peak = local
   `PSR` (6.3).

### How to pick station points (answers "how do you pick points")

Two DISTINCT rules, because the two station types have different placement constraints:
- OFF-FAULT receivers (seismograms) may sit ANYWHERE — SeisSol finds the containing tetrahedron by
  point search — so their map coordinates need not be exact.
- ON-FAULT pickpoints MUST lie on the (curved, segmented) fault surface, so you NEVER hand-type
  them on a curved fault; you snap each target to the nearest fault facet.

PLACEMENT RATIONALE (where, and what physics each station captures):
- **Epicentral station**: the surface-trace point nearest the epicenter `(606971, 3707270)` —
  near-source reference; its seismogram peaks ARE the 6.1 map values there.
- **Fault-normal profile**: a line PERPENDICULAR to the local fault trace at the epicenter, at
  +-2, +-5, +-10, +-20 km. Captures distance decay and across-fault asymmetry.
- **Along-strike (directivity) line**: stations following the surface trace in the rupture
  direction (NW), at +10, +20, +40 km, offset ~1-2 km off the trace so they are genuinely off-fault.
  Captures forward directivity (the velocity pulse sharpening ahead of the front) and, where Vr>Vs
  (6.3), the supershear Mach-cone PGV high.
- **On-fault pickpoints**: a down-dip column at the hypocenter (depth 2, 5, 8, 12 km) for depth
  dependence, plus along-strike points at ~7 km depth (s = 0, +10, +20, +40 km) to read the local
  slip-weakening curve (6.3) and Vr from RT differences between them.

ROUTE 1 (preferred, reproducible, PRE-run): `make_stations.py` (Phase 1). It reads the fault
surface from `safs_mesh.puml.h5` (PUML: vertices + connectivity + boundary tags; select the fault
tag = 3) and then:
  (a) extracts the SURFACE TRACE = fault points with z ~ 0, ordered into a polyline by a
      nearest-neighbor walk from the epicenter, with arc length along it;
  (b) places along-strike receivers at the target arc lengths, and the fault-normal profile from
      the local trace tangent (rotate 90 deg in the map plane);
  (c) snaps each on-fault pickpoint target (along-strike s, depth d) to the nearest fault-facet
      CENTROID (`scipy.spatial.cKDTree`) -> guaranteed on-fault (x,y,z);
  (d) writes `safs_receivers.dat` (z=0) and `safs_pickpoints.dat`.

ROUTE 2 (fallback, no scripting, costs one extra run): run once with the Phase-1 starter stations;
open `safs-surface.xdmf` and use ParaView's Spreadsheet View / "Find Data" to read surface-trace
coordinates for the along-strike / fault-normal receivers; open `safs-fault.xdmf` and read facet
centroids near your targets for the on-fault pickpoints; paste into the `.dat` files; re-run.

VERIFY after the run: the log lists each receiver / pickpoint as FOUND. Any pickpoint reported
off-fault (search tolerance) is dropped silently — replace it from a nearer facet centroid.

### Deepen-understanding probe
"From a single receiver at 20 km you measured Vp = 5.9 km/s; the model says 6.0 — what are the two
or three reasons for the 2% gap?" (finite first-break pick, pickdt quantization, mesh dispersion,
path not straight.) "You have the surface PGV map AND this receiver's seismogram — which is primary
and which is derived, and could you rebuild the map from enough seismograms?" (Seismogram primary;
yes — the map is the field of peaks.)

### Acceptance Criteria
- [ ] One receiver seismogram plotted with P/S arrivals identified.
- [ ] Integrate->displacement matches the surface u pixel; differentiate->PGA matches 6.1.
- [ ] Vp recovered from a first-arrival time and compared to the material Vp.
- [ ] Pickpoint slip-rate integral == ASl at that point (links to 6.2/6.3).

### Dependencies
- Depends on: Phase 1 (receivers + pickpoints + surface + fault). Cross-links: 6.1, 6.2, 6.3.

---

## Phase 6: WORKSHEET + ANSWER_KEY + smoke test + budget verification

### Goal
Two participant-facing docs and a facilitator pre-check, all consistent with the LSW case.

### Files to Create
- `WORKSHEET.md` (+ optional `.pdf`) — Activities 6.1-6.4, in order, constant material first then
  CVM, with the official GME-script command (6.1), the exact ParaView filters, and the Python
  snippets for receivers and `make_stations.py`, every step naming the SeisSol ground truth it
  cross-checks and the other activity it links to.
- `ANSWER_KEY.md` (+ optional `.pdf`) — reference numbers + tolerances + every expected
  discrepancy's cause; filled after the official runs.

### Detailed Requirements
1. **Smoke test BEFORE distribution** (the LSW candidate D is known to be touchy — see
   `safs_fault.yaml` History: candidates B/C ARRESTED): run `parameters_postprocess.par` on the
   CONSTANT material and confirm the energy CSV `seismic_moment` GROWS past the nucleation patch
   (rupture propagated, not a confined patch). If it stalls, this is a case-tuning issue, NOT a
   deck bug — report it; do not silently swap parameters (per project rules).
2. Record per material (constant, CVM): final P0, seismic_moment, Mw; P0 (Integrate vs Dbar*A_rup
   vs CSV) agreement %; A_rup + Dbar (with eps); Vr interior median + supershear fraction; PSR
   (ParaView vs SeisSol); PGV/PGD/PGA from the official GME script AND from ParaView, with the
   official-vs-ParaView gap; recovered Vp vs model.
3. Tolerance table: P0 agreement < 1%; M0 = 3.2e10*P0 band (constant); Vr interior median < 10%
   (edges ungraded); PSR(ParaView) <= PSR(SeisSol); PGD(official integrated-v) vs PGD(ParaView
   direct-u) < few %; PGV(official GMRotD50) vs PGV(ParaView resultant) gap REPORTED not graded
   (different definitions); Vp recovery < ~5%.
4. Explain each expected mismatch (quadrature; polynomial vs manual gradient; temporal
   undersampling/Nyquist; mesh frequency ceiling; integration baseline drift; horizontal-GMRotD50
   vs 3D-resultant; CVM mu-weighting).

### Acceptance Criteria
- [ ] Smoke test shows nucleation + propagation on constant material (or the arrest is reported).
- [ ] Every worksheet quantity has a reference value + tolerance + discrepancy cause.
- [ ] `du -sh output/` confirms < 1 GB for the official run(s).

### Dependencies
- Depends on: Phases 1-5 + the official run(s).

---

## Testing Strategy
- **Deck validity**: `parameters_postprocess.par` parses; log shows FL=16, material, receiver and
  pickpoint counts; outputs listed in Phase 1 acceptance appear.
- **Budget**: `du -sh output/` is 0.4-0.8 GB at 100 s.
- **Core self-consistency**: P0(Integrate ASl) == P0(Dbar*A) == P0(CSV) and (constant)
  M0 = 3.2e10*P0 == seismic_moment(CSV), all < 1%.
- **Cross-activity web** (the real test): PGV(map pixel)==PGV(seismogram peak); PGD(surface
  u)==PGD(integral v); PSR==peak(pickpoint slip rate); pickpoint slip-rate integral==ASl;
  recovered Vp ~= material Vp. Each is an independent confirmation the pipeline is correct.
- **Nucleation/propagation**: energy CSV moment grows (smoke test).
- **Reproducibility**: a second intern on the same files gets the same numbers (deterministic).

## Risk Assessment
- **(Highest) LSW candidate D may arrest** (B/C did — `safs_fault.yaml` History). The forced
  nucleation patch always slips, so all four activities have non-zero data regardless, but a
  confined event makes the supershear/PGV story thin. Detect via the smoke test; report arrest
  (do NOT retune silently — project rule).
- **EndTime=100 s may truncate the far-NW rupture** -> M0 not plateaued. The P0/M0 verification is
  still EXACT at any shared snapshot (compare ParaView and CSV at the same t); only the "final
  event size" interpretation needs a longer run (extend to 150 s WITH the size levers).
- **Receiver/pickpoint `.dat` parser strictness** (comments) and **pickpoint off-fault drops** —
  handled in Phase 1 edge cases; verify against the run-1 log before the activity.
- **PGA is mesh-limited even at fine pickdt** (~1-2 Hz on a 500 m mesh) — teach this as a result,
  not a bug; do not chase a "true" PGA the mesh cannot carry.
- **ParaView gradient noise** on the curved/segmented SAF and at fault edges — mask `RT>0`, compare
  interior medians, not per-cell extremes.
- **CVM M0 not exactly verifiable** in ParaView (mu varies, not output) — constant is the exact
  anchor by design; CVM M0 via the optional resample-mu extension only.
- **Stale comment in `safs_material.yaml`** ("Vp ~= 5.39 km/s") — the correct constant Vp is
  5996 m/s (used by the 6.4 arrival-time check). Optionally fix the comment; do not let interns
  cross-check against the wrong number.

## Deliverables checklist (in safs_seisol_v2_2_0_LSW/)
- [ ] `parameters_postprocess.par`        (Phase 1)
- [ ] `make_stations.py`                   (Phase 1, station generator)
- [ ] `safs_receivers.dat`                 (Phase 1, from make_stations.py)
- [ ] `safs_pickpoints.dat`                (Phase 1, from make_stations.py)
- [ ] `safs_material.yaml`                 (Phase 1, copied from archive/)
- [ ] `WORKSHEET.md` (+ `.pdf`)            (Phases 2-5)
- [ ] `ANSWER_KEY.md` (+ `.pdf`)           (Phase 6, after official run + smoke test)
- [ ] `POSTPROCESS_ACTIVITY_PLAN.md` (this file)
- [ ] (unchanged) `parameters.par`, `parameters_kinematics.par`, `safs_fault.yaml`,
      `safs_initial_stress.yaml`, `safs_material_cvm.yaml`, `safs_material_cvm.nc`, `safs_mesh.puml.h5`
```
