# SeisSol on Frontera — TPV102 reference run guide

**Goal:** build SeisSol on TACC Frontera and run the SCEC TPV102 benchmark
(rate-and-state aging law, 3-D vertical strike-slip fault) so we can
compare the output step-by-step against the seas-mfem TPV102 driver and
debug divergences.

**Status of research:** online exploration only.  All commands below are
taken from the official SeisSol docs or the SeisSol training repository.
Treat as a working checklist — anything marked `[verify on-machine]` is
drawn from a secondary source and should be confirmed the first time the
step is run.

---

## 0. Why SeisSol for the TPV102 reference

- SeisSol is the SCEC TPV101/102 aging-law reference code (LMU Munich /
  TUM, used extensively for community verification).  A SeisSol run
  that passes the SCEC TPV102 output format is the exact thing the
  Code Verification Project compares seas-mfem against.
- **TPV101 vs TPV102 per the SCEC spec:** both benchmarks are 3-D, both
  use the rate-and-state aging law on a vertical strike-slip fault,
  and they share identical friction parameters, fault dimensions,
  material, and nucleation pulse.  **The only difference is the
  medium: TPV101 = whole space (symmetric about z=0, no free surface);
  TPV102 = half space (free surface at z=0).**  That means SeisSol's
  `Examples/tpv101/` is the direct basis — we only have to swap the
  mesh so the top of the domain is a free surface and extend the
  symmetric box to a half-space box.  `parameters.par` and
  `tpv101_fault.yaml` carry over essentially verbatim (`FL = 3`).
- A single-rank SeisSol run on a 1000 m TPV102 mesh (~30k tets) takes
  under 1 hour on a Cascade Lake node and produces per-station traces
  we can diff against seas-mfem's `station_*.dat` output directly.

## 1. Three build options on Frontera

Pick one based on how much customisation you need.  **Option A is the
recommended starting point** — running a TPV102 mesh end-to-end inside
the container validates the whole pipeline before you touch a compiler.

### Option A: Apptainer (Singularity) container — fastest

The official SeisSol training image has PUMGen + Gmsh + SeisSol + rconv
+ easi + ASAGI pre-baked for Cascade Lake (Frontera's CLX nodes).

```bash
# On a Frontera login node
module load tacc-apptainer
apptainer pull -F docker://seissol/training:latest
ln -sf $(realpath latest.sif) ~/seissol.sif
```

To get an interactive compute node (required — containers cannot launch
MPI from the login node):

```bash
idev -m 30 -N 1 --tasks-per-node 2 -p development
```

Smoke-test with the bundled TPV13 scenario:

```bash
cd seissol-training/tpv13        # comes with the image
mpirun apptainer run ~/seissol.sif gmsh -3 tpv13_training.geo
mpirun apptainer run ~/seissol.sif pumgen -s msh2 tpv13_training.msh
OMP_NUM_THREADS=26 mpirun -n 2 apptainer run ~/seissol.sif seissol parameters.par
```

If that run completes and produces `output/tpv13-fault.xdmf`, the
container is good.

### Option B: Native build against the preinstalled `seissol-env` module

Frontera has a Spack-managed `seissol-env` module with all binary deps
(HDF5 parallel, NetCDF parallel, Eigen, easi+ASAGI, libxsmm, PSpaMM,
ImpalaJIT, Lua).  This is the recommended path once you need to patch
SeisSol source or run non-default configs.

```bash
# Load env (Frontera documentation path; confirm with `ls` once)
module unload xalt
module switch python3 python3/3.9.2
module use /work2/09160/ulrich/frontera/spack/share/spack/modules/linux-centos7-cascadelake
module load seissol-env
export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort
```

Build SeisSol:

```bash
cd $WORK
git clone --recursive https://github.com/SeisSol/SeisSol.git
cd SeisSol
mkdir build && cd build
cmake ..  \
    -DHOST_ARCH=skx            \
    -DORDER=4                  \
    -DNUMBER_OF_MECHANISMS=0   \
    -DEQUATIONS=elastic        \
    -DPRECISION=double         \
    -DCMAKE_BUILD_TYPE=Release \
    -DDR_QUAD_RULE=dunavant
make -j 16
# Binary: build/SeisSol_Release_dskx_4_elastic   (name depends on flags)
```

`HOST_ARCH=skx` (Skylake) is the correct value for Frontera's CLX nodes —
SeisSol's cookbook uses this value, and Cascade Lake is backward-
compatible with the Skylake code generator kernels.

### Option C: From-scratch Spack install (reference — ~half-day)

Only needed if `seissol-env` is unavailable or you want specific versions.

```bash
module load python3/3.9.2 gcc/9 intel/19.1.1 impi/19.0.9 cmake
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
source spack/share/spack/setup-env.sh
git clone https://github.com/SeisSol/seissol-spack-aid.git
spack repo add seissol-spack-aid/spack
spack install seissol-env %intel@19 target=cascadelake ^[virtuals=mpi] intel-oneapi-mpi
spack load seissol-env
```

Then proceed with the CMake from Option B.  `[verify on-machine]` —
`seissol-spack-aid` is the upstream recipe; the exact Spack `install`
syntax above is from the SeisSol docs and may need a version tag.

## 2. Verify the container / binary works end-to-end

```bash
# Inside idev session, after Option A or B
which seissol        # should point at the sif or the build dir
seissol --help       # prints parameter layout
```

If `--help` lists `FL` / `DynamicRupture` / `Elementwise` namelists,
you're good.

## 3. TPV102 input files

SeisSol ships TPV101 (whole space); TPV102 (half space) is identical
**except the mesh**.  The physics config — `parameters.par`,
`tpv101_fault.yaml`, `tpv101_material.yaml`, station files — is reused
as-is.

- **`SeisSol/Examples/tpv101`** — 3-D aging-law rate-state reference
  (whole space).  Provides `parameters.par`, `tpv101_fault.yaml`,
  `tpv101_material.yaml`, `tpv101_receivers.dat`, and
  `tpv101_faultreceivers.dat`.
- **SCEC TPV101/102 official spec** at
  `https://strike.scec.org/cvws/tpv101_102docs.html` (links to PDFs for
  problem description, station list, and file-format instructions).

Workflow:

```bash
mkdir -p $WORK/tpv102_seissol && cd $_
git clone https://github.com/SeisSol/Examples.git
cp -r Examples/tpv101/* .
# Rename for clarity (optional, pure bookkeeping):
mv tpv101_fault.yaml         tpv102_fault.yaml
mv tpv101_material.yaml      tpv102_material.yaml
mv tpv101_receivers.dat      tpv102_receivers.dat
mv tpv101_faultreceivers.dat tpv102_faultreceivers.dat
sed -i 's/tpv101_/tpv102_/g' parameters.par
# The only physics-relevant edit is the mesh — see §3.4.
```

### 3.1 TPV102 physical parameters (from SCEC spec)

Keep these as the authoritative source; `parameters.par` and
`tpv102_fault.yaml` must match.  **All of these are identical to TPV101
— the only TPV101→TPV102 change is the domain topology (whole → half
space) and therefore the mesh.**

| Quantity                       | Value                                             |
|---|---|
| Fault plane                    | vertical strike-slip, y=0                         |
| Fault dimensions               | 30 km along-strike × 15 km down-dip               |
| Domain (TPV101 — whole space)  | symmetric about z=0, 60 × 40 × 30 km, absorbing   |
| Domain (TPV102 — half space)   | z ∈ [-30, 0] km, top is a **free surface**        |
| Element size (recommended)     | 100 m on fault, grading outward                   |
| Material (homogeneous)         | ρ = 2670 kg/m³, Vp = 6000 m/s, Vs = 3464 m/s     |
| Effective normal stress σₙ₀    | 120 MPa (compression)                             |
| Initial shear τ₀ (along-strike)| 75 MPa                                            |
| Rate-state: f₀                 | 0.6                                               |
| Rate-state: a                  | 0.008 (VW patch), 0.016 (VS outside)              |
| Rate-state: b                  | 0.012                                             |
| Rate-state: L (Dc)             | 0.02 m                                            |
| Rate-state: V₀ (ref slip rate) | 1e-6 m/s                                          |
| Initial slip rate V_ini        | 1e-12 m/s along-strike                            |
| Nucleation zone                | sphere radius 3 km at (x=0, z=-7.5 km)            |
| Nucleation perturbation        | Δτ_max = 25 MPa, t_ramp = 1 s (Gaussian in time)  |
| End time                       | 12–15 s                                           |

Both TPV101 and TPV102 use the **aging law** in SeisSol: `FL = 3` in
`parameters.par`.  (TPV103 and TPV104 use the slip-law variants, `FL = 4`
and `FL = 103` respectively — not relevant here.)

### 3.2 Minimal `parameters.par` edits (derived from SeisSol tpv101)

Starting from the SeisSol `tpv101/parameters.par`, keep:

```
&Equations
FrictionLaw = 3                 ! aging law (TPV101 + TPV102)
/

&DynamicRupture
OutputPointType = 4             ! on-fault receivers + elementwise
FL = 3
/

&Elementwise
printtimeinterval_sec = 0.5
OutputMask = 1 1 1 1 1 1 1 1 1 1 1 1
refinement_strategy = 2
/

&Pickpoint
PPFileName = 'tpv102_faultreceivers.dat'
/

&Output
FileName = 'output/tpv102'
EndTime = 12.0
printIntervalCriterion = 2
printtimeinterval_sec = 0.5
OutputRegionBounds = -30000. 30000. -15000. 15000. -30000. 0.
ReceiverOutputInterval = 0.005
RFileName = 'tpv102_receivers.dat'
/

&Mesh
MeshFile = 'mesh/tpv102_100m'   ! PUML format, produced by pumgen
/
```

### 3.3 `tpv102_fault.yaml` (adapted from tpv101_fault.yaml)

Essential adjustments for TPV102:

```yaml
!LuaMap
returns: [s_xx, s_yy, s_zz, s_xy, s_yz, s_xz,
          rs_a, rs_sl0, rs_srW, rs_f0, rs_b,
          cohesion, forced_rupture_time,
          nuc_s_xx, nuc_s_yy, nuc_s_zz,
          nuc_s_xy, nuc_s_yz, nuc_s_xz]
function: |
  function f (x)
    local along = x.x
    local down  = math.abs(x.z)
    local r     = math.sqrt(along^2 + (down - 7500)^2)
    local nuc_xy = 0.0
    if r < 3000.0 then
       nuc_xy = 25.0e6 * math.exp(r*r / (r*r - 9.0e6))
    end
    -- VW/VS transition: box of 30 × 15 km centred at (0, -7.5 km)
    local inside = 1.0
    if math.abs(along) > 15000 then inside = 0.0 end
    if down > 15000 then inside = 0.0 end
    local rs_a = inside == 1.0 and 0.008 or 0.016
    return {
      s_xx = 0.0, s_yy = -120.0e6, s_zz = 0.0,
      s_xy = 75.0e6, s_yz = 0.0, s_xz = 0.0,
      rs_a = rs_a, rs_sl0 = 0.02,
      rs_srW = 1.0e-12, rs_f0 = 0.6, rs_b = 0.012,
      cohesion = 0.0, forced_rupture_time = 0.0,
      nuc_s_xy = nuc_xy,
      nuc_s_xx = 0.0, nuc_s_yy = 0.0, nuc_s_zz = 0.0,
      nuc_s_yz = 0.0, nuc_s_xz = 0.0
    }
  end
```

### 3.4 Mesh generation — the TPV101→TPV102 edit

**This is the only content-bearing change from TPV101.**  TPV101's mesh
is a box symmetric about z=0 (so the top "surface" is just an internal
plane with continuation, i.e., a numerical absorbing or reflective
boundary); TPV102's mesh is a half-space box whose top face at z=0 is
tagged as a free-surface boundary.

In Gmsh tags (SeisSol convention), boundary surfaces are numbered:

| tag | BC type          | TPV101 usage            | TPV102 usage           |
|---|---|---|---|
| 1   | free surface     | (not used)              | **top face z=0**       |
| 3   | dynamic rupture  | fault plane y=0         | fault plane y=0        |
| 5   | absorbing        | all 6 outer faces       | 5 outer faces (not top)|

**3.4a. Gmsh + PUMGen (simplest, what the training uses):**

```bash
# tpv102_100m.geo:
#   - box z ∈ [-30000, 0] m (half-space), x ∈ [-30000, 30000], y ∈ [-20000, 20000]
#   - planar fault stamped at y=0 from (-15000, 0, -15000) to (15000, 0, 0)
#   - Physical Surface 1 = top face at z=0           (free surface)
#   - Physical Surface 3 = fault plane y=0 within    (dynamic rupture)
#   - Physical Surface 5 = remaining outer faces     (absorbing)
#   - fault face mesh size 100 m, growing to 2000 m at outer boundary
gmsh -3 tpv102_100m.geo
apptainer run ~/seissol.sif pumgen -s msh2 tpv102_100m.msh tpv102_100m
# Output: tpv102_100m.xdmf + tpv102_100m.h5
```

For TPV101, the only differences in the `.geo` would be: the z-range is
symmetric `[-15000, 15000]`, and **no** Physical Surface 1 is created —
all six outer faces get tag 5 (absorbing).  Keeping both geo files side
by side makes the whole-space vs half-space diff obvious in review.

**3.4b. SimModeler** (if you have a licence) — not covered here.

### 3.5 Station file (`tpv102_faultreceivers.dat`)

SCEC TPV102 requires 9 on-fault stations (3 along-strike × 3 down-dip).
Use the same (along, depth) grid we already have in
`miniapps/seas/dynamic/tpv102_setup.hpp::DefaultStations` to keep the
diff trivial.  The file format is one line per station, `x y z`
physical coordinates (z negative for depth).

## 4. Running TPV102

Minimal sbatch for Frontera dev queue (~2 hr):

```bash
#!/bin/bash
#SBATCH -J seissol_tpv102
#SBATCH -o seissol_tpv102_%j.out
#SBATCH -e seissol_tpv102_%j.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -A EAR20006       # same allocation seas-mfem uses

module load tacc-apptainer
export I_MPI_SHM_HEAP_VSIZE=32768
export OMP_NUM_THREADS=26
export OMP_PLACES="cores(27)"
export OMP_PROC_BIND="close"

cd $WORK/tpv102_seissol
mpirun -n 2 apptainer run ~/seissol.sif seissol parameters.par
```

**Expected runtime:** 1–2 hr for 100 m, 12-second window, on 2 ranks × 26
OMP threads per node (the canonical Frontera configuration for SeisSol).

**Output:**

- `output/tpv102-fault.xdmf` (+ `.h5`) — full fault-surface time series.
- `output/tpv102-surface.xdmf` — free-surface velocities for comparison
  with seas-mfem's surface-station writer.
- `output/tpv102-receiver-NNNNN-NNNNN.dat` — per-station ASCII time series.
  Direct diff target for `seas-mfem/tpv102/results_*/tpv102_station_*.dat`.

## 5. Comparing SeisSol vs seas-mfem bit by bit

Proposed staircase — each level diff'd before moving up.

| Level | Quantity                                    | SeisSol source                        | seas-mfem source                                |
|---|---|---|---|
| 0     | Mesh nodes on fault                         | PUML h5 → extract fault verts         | Gmsh .msh → extract fault verts                |
| 1     | Initial DOFData (τ, σₙ, ψ, V)               | `output/tpv102-fault-receiver-0.dat` row 0 | `station_flt_0_3.dat` row 0                |
| 2     | t=0.5 s fault state (post-nucleation)       | fault receiver series at t=0.5        | station series at t=0.5                         |
| 3     | Peak slip rate, peak traction drop          | fault receiver max over 0..12 s       | station max over 0..12 s                        |
| 4     | Free-surface `v_x` at (0, 0, 0) over time   | off-fault receiver                    | `surface_station_0_0.dat`                       |
| 5     | Integrated moment rate M₀'(t)               | derived from fault output             | derived from DOFData V integral                 |

Run levels 0–2 BEFORE submitting a 400-rank Frontera job — they catch
most configuration errors cheaply.

## 6. Recommended first-run scripts to commit to the repo

Place under `miniapps/seas/jobs/seissol_reference/`:

- `seissol_tpv102_1000m_container.sbatch` — Option A container, 1000 m
  mesh (coarse smoke, ~30 min).
- `seissol_tpv102_200m_native.sbatch` — Option B native build, 200 m
  mesh (2 hr dev, matches our own 200m production grid).
- `tpv102_fault.yaml`, `tpv102_material.yaml`, `tpv102_receivers.dat`,
  `tpv102_faultreceivers.dat` — config files generated per §3.
- `compare_seissol_vs_seasmfem.py` — Python diff over station output;
  wraps numpy/matplotlib.  A good fit for the existing
  `miniapps/seas/jobs/tpv102/test_result_regex.sh` pattern (shell +
  one diagnostic per test).

## 7. Troubleshooting quick-reference

| Symptom                                                | Likely cause / fix                                                         |
|---|---|
| `module load seissol-env` fails                        | Path may have moved; `find /work2 -name seissol-env.lua 2>/dev/null`       |
| `seissol` exits with `Unknown keyword FrictionLaw`     | Version mismatch — the older flag was `FL`; try both                       |
| Per-rank output sits empty                             | Forgot `OutputRegionBounds` or bounds clip the fault                       |
| Nucleation never triggers                              | `nuc_s_xy` magnitude too low or `r > 3 km` for all QPs (check mesh origin) |
| Container hangs at startup                             | Must be inside `idev` / sbatch — container MPI needs a compute node        |
| PUMGen segfault                                        | Input `.msh` must be msh v2 format: `gmsh -3 -format msh2` or `-s msh2`    |

## 8. Open items / things to verify on machine

- [ ] `seissol-env` exact path still current — test via `module avail seissol`.
- [ ] SeisSol `FL` value accepted for TPV102 (aging) — confirm with a
      working tpv101 run first; TPV102 uses the same FL.
- [ ] Container image tag (`seissol/training:latest`) is not stale; may
      need `seissol/seissol:v1.3.0` or similar for a published release.
- [ ] TPV101's `mesh_file` in the bundled `parameters.par` references a
      mesh we have not built; TPV102 needs its own PUML mesh from §3.4.

---

**Next concrete action when you're on Frontera:**

1. `idev` → `module load tacc-apptainer` → `apptainer pull` (§1 Option A).
2. Run the bundled `tpv13` scenario (§1) to prove the container works.
3. Clone `SeisSol/Examples`, copy `tpv101/` → `tpv102/`, edit per §3.
4. Produce a 1000 m mesh via Gmsh+PUMGen (§3.4a), launch a short
   `tfinal=1.5s` run (§4).
5. Fetch the fault receiver at `(0, 0, -3000)` and line it up with
   our `station_flt_0_3.dat` row-for-row (§5 level 1–2).

## 9. Source-code cross-reference (local clone for reading)

The container packages SeisSol for **execution** on Frontera; it does
NOT prevent you from cloning the source locally for reading and diff.
Recommended workflow: local clone on your laptop for grep / editor
side-by-side with `seas-mfem`, container on Frontera for runs.  The
only requirement is to **pin the local clone to the same git tag** the
container was built from — otherwise numerical details may differ
between what you're reading locally and what's actually executing on
Frontera.

### 9.1 What to clone locally

```bash
mkdir -p ~/seissol-reference && cd ~/seissol-reference
git clone --recursive https://github.com/SeisSol/SeisSol.git
git clone            https://github.com/SeisSol/Examples.git
git clone            https://github.com/SeisSol/easi.git
git clone            https://github.com/SeisSol/PUMGen.git
git clone            https://github.com/SeisSol/Training.git
```

### 9.2 Pin the local clone to the container's version

```bash
# On Frontera, after pulling the container:
apptainer inspect ~/seissol.sif | grep -i version
# Prints e.g.  org.opencontainers.image.version=1.3.0

# On your laptop:
cd ~/seissol-reference/SeisSol
git checkout v1.3.0
git submodule update --init --recursive
```

If the training image lacks a version label, match the image build date
to a SeisSol tag via <https://github.com/SeisSol/SeisSol/releases>.

### 9.3 Highest-value files for direct comparison with seas-mfem

| SeisSol file                                                                                | seas-mfem counterpart                                                 | Why                                                                                  |
|---|---|---|
| `src/Physics/Evaluate_friction_law.f90` (older) or `src/DynamicRupture/FrictionLaws/RateAndState*.cpp` (newer) | `friction/dieterich_ruina.hpp`, `dynamic/friction_solver.cpp`         | Aging-law Brent/Newton solve — the numerical heart of the comparison                 |
| `src/DynamicRupture/Initializer/BaseDRInitializer.cpp`                                      | `dynamic/tpv102_setup.hpp::InitializeFaultDOFs`                       | How initial ψ is computed from stress equilibrium                                    |
| `src/DynamicRupture/Output/OutputManager.cpp`                                               | `dynamic/tpv102_setup.hpp::TPV102StationWriter`                       | On-fault receiver sampling convention (critical for diffs)                           |
| `src/Numerical_aux/ODEInt.*`                                                                | `drivers/tpv102_driver.cpp` (RK4 stage loop)                          | State-variable time integration (aging- vs slip-law ODE)                             |
| `src/Initializer/BoundaryConditions/FreeSurfaceIntegrator.*`                                | `dynamic/wave_operator.inl` `FaceBC::FreeSurface` branch              | Free-surface BC implementation (the TPV101→TPV102 diff!)                             |
| `src/Solver/FrictionIntegrator.cpp`                                                         | `dynamic/wave_operator.inl` fault-flux branch                          | How friction traction feeds back into the wave equation                              |
| `src/Kernels/Time.cpp` / `src/Model/PlasticityEqs.*`                                        | `dynamic/godunov_flux.cpp`                                            | Godunov upwind numerical flux — our `Interior/Absorbing/FreeSurface` mirrors this    |
| `Examples/tpv101/tpv101_fault.yaml`                                                         | `config/tpv102_params.hpp` + `dynamic/tpv102_setup.hpp`               | Parameter authority — `rs_a, rs_b, rs_sl0, rs_f0, rs_srW` numerical values          |
| `Examples/tpv101/parameters.par`                                                            | `drivers/tpv102_driver.cpp` CLI arg parsing                           | Output cadence, CFL, mesh filename, friction-law selector                            |

### 9.4 Mounting your local source into the container (for patch testing)

If during diff-debugging you find a SeisSol bug or want to try a
modification without rebuilding the whole image, bind-mount your local
source into the container at runtime:

```bash
apptainer run -B /your/local/SeisSol:/src ~/seissol.sif bash -c \
  "cd /src && mkdir build && cd build && cmake .. -DHOST_ARCH=skx \
   -DORDER=4 -DCMAKE_BUILD_TYPE=Release && make -j"
```

This lets you test a patch against the container's pinned dependencies
before upstreaming.  Useful once the baseline TPV101 reproduction is
working and you want to explore specific numerical paths (e.g., the
aging-law Brent solve) in isolation.

## Sources

- [SeisSol Frontera docs](https://seissol.readthedocs.io/en/latest/frontera.html)
- [SeisSol Training — frontera.md](https://github.com/SeisSol/Training/blob/main/frontera.md)
- [SeisSol Training — setup_modules_Frontera_vnc.sh](https://github.com/SeisSol/Training/blob/main/setup_modules_Frontera_vnc.sh)
- [SeisSol Build overview](https://seissol.readthedocs.io/en/latest/build-overview.html)
- [SeisSol Build dependencies](https://seissol.readthedocs.io/en/latest/build-dependencies.html)
- [SeisSol Dynamic rupture](https://seissol.readthedocs.io/en/latest/dynamic-rupture.html)
- [SeisSol A first example](https://seissol.readthedocs.io/en/latest/a-first-example.html)
- [SeisSol Examples repo — TPV101](https://github.com/SeisSol/Examples/tree/master/tpv101)
- [SeisSol Examples repo — TPV104 (rate-state reference)](https://github.com/SeisSol/Examples/tree/master/tpv104)
- [PUMGen mesh conversion](https://seissol.readthedocs.io/en/latest/meshing-with-pumgen.html)
- [SCEC TPV101/102 problem description](https://strike.scec.org/cvws/tpv101_102docs.html)
- [Frontera TACC user guide](https://docs.tacc.utexas.edu/hpc/frontera/)
