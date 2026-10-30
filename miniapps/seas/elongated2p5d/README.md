# Level 2: 2.5D map view against a 3-D elongated rupture

The 2.5D note's Level 2 test (Chunhui Zhao, 2026-10-30; note rev8,
`/mnt/project-files/2d-inplane/nonplanar-literature/inplane-2p5d-derivation-rev8.pdf`,
Section 7.2).

- **Question:** does the map-view model with the 2.5D finite-width term,
  `beta = 4 c^2 mu / W^2`, give the depth-averaged slip and rupture arrival of
  a 3-D rupture of width W?
- **Rupture:** a buried vertical strike-slip fault, L = 120 km by W = 15 km
  (L / W = 8). TPV205 material and background stress (sigma_yy = -120 MPa,
  sigma_xy = 70 MPa), linear slip weakening (0.677 / 0.525, Dc 0.40 m), and
  only the 3 km nucleation patch at 81.6 MPa. No free surface, so the buried
  values of c apply.
- **Pass (the note):** for every station with |x| between 15 and 57.5 km
  (where the rupture is longer than 2 W):
  - map-view slip within 10% of the 3-D slip averaged over the width;
  - arrival within 5% of the 3-D arrival.
- **Also reported:** the implied c = c_run x D_2D / D_3D. Weng & Ampuero
  (2019) fit their 3-D runs with gamma = 0.96 / pi.

## Runs

| Deck | Model | c | beta (Pa/m^2) | Where |
|---|---|---|---|---|
| `configs/elongated_3d_buried.toml` | 3-D reference, 400 m mesh, p = 2 | (none) | (none) | Vista, 16 nodes |
| `configs/elongated_2d_2p5d_wa2020.toml` | map view + 2.5D | (1 - nu) pi / 2 = 1.178 | 790.5 | Mac |
| `configs/elongated_2d_2p5d_lehner.toml` | map view + 2.5D | 4 / pi = 1.273 | 923.3 | Mac |
| `configs/elongated_2d_2p5d_liang.toml` | map view + 2.5D | pi / 2 = 1.571 (for information) | 1405 | Mac |
| `configs/elongated_2d_plane_strain.toml` | map view, no term (the control) | (none) | (none) | Mac |

- **Map-view decks:** all four are the same but for `[finite_width]`, the tag
  and the output folder. They use p = 2 on the committed 200 m mesh and run
  to 40 s.
- **Shared settings:** all five decks have Kelvin-Voigt kv_eta 0.1 and
  absorbing walls on every side.
- **Stations:**
  - map view: every 2.5 km along strike (`stations/elongated_2d_stations.csv`);
  - 3-D: at the same x, at 15 depths from z = -7 to 7 km, 1 km apart
    (`stations/elongated_3d_stations.csv`).

## Files

- `mesh/elongated_2d.geo`, `mesh/elongated_2d_200m.msh`: the map-view mesh.
  - It is committed; sha256 `5d74930bb73a9145302190657c36213da9342e8bf814b7a9f9fd8c472991afbe`.
  - It has 62,068 triangles and a 210 x 120 km box.
- `mesh/elongated_3d.geo`: the 3-D mesh, made on Vista by the build job.
  - The mesh file (128 MB) is not committed.
  - Its sha256 goes to `MESH_SHA256_VISTA`, which the run job checks.
  - Cloud gmsh 4.15.2 gave 2,409,373 tets and 25,660 fault faces at 400 m.
- `jobs/vista_l2_build.sbatch`: Vista part 1 (1 node, 2 h or less, 0.66 SU or less, EAR20006). It:
  - builds PETSc, MFEM and `seas_nitsche_dyn_driver` at the checkout's commit
    and writes `BUILD_COMMIT_VISTA_NITSCHE`;
  - runs `make test-2d5-dyn-driver`;
  - meshes 1000 m and 400 m with TACC gmsh;
  - runs a 0.5 s smoke on the 1000 m mesh;
  - runs a load check on one node with the 400 m mesh (node memory, dt).
- `jobs/elongated_3d_vista.sbatch`: Vista part 2, the 3-D run (16 nodes x 144 ranks).
  - Output goes to `<checkout>-runs/l2_3d_<job>`: stations every 4th step,
    fault CSV dumps every 4000 steps, `RUN_INFO.txt` and the memory log.
  - The 3 h limit costs 15.8 SU at most.
  - Estimate: about 1.7 h, or 9 SU. This is scaled from the TPV205 200 m p = 2
    gold and not measured; the build job's load check gives dt.
- `scripts/run_2d_mac.sh`: the four map-view runs on the Mac, one after the other.
  - It holds the Mac run lock and uses OMP_NUM_THREADS=1.
  - About 48 min each on 4 ranks, about 3.2 h in all (from the cloud's 0.2 s
    smoke: 340 steps in 14.5 s).
- `scripts/compare_2p5d_3d.py`: the comparison. It writes `summary.txt`,
  `summary.json` and `compare_2p5d_3d.pdf`.

## Workflow

Every step runs a pushed commit from a clean checkout (`CLAUDE.md`). Each
Vista job needs Chunhui Zhao's typed go.

1. **Vista build.** In the `$WORK/seas-mfem-2d-inplane` checkout at the
   commit, from `miniapps/seas/elongated2p5d/jobs/`:
   `sbatch vista_l2_build.sbatch`.
   - Read its `.out`: test results, the mesh sizes, the smoke, the node peak
     memory and the 400 m dt.
2. **Vista 3-D run.** From the same folder:
   `sbatch --dependency=afterok:<build job> --kill-on-invalid-dep=yes elongated_3d_vista.sbatch`.
   - Or submit after the build has passed.
   - It refuses if code changed since the build or if the mesh sha256 does
     not match.
3. **Mac map-view runs.** From `miniapps/seas` of the Mac checkout at the
   commit (`conda activate mfem-dev`):
   ```bash
   make -f tests/2d_inplane/2d_inplane.mk seas_nitsche_dyn_driver
   git rev-parse HEAD > BUILD_COMMIT_MAC_NITSCHE
   elongated2p5d/scripts/run_2d_mac.sh /Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/level2 4
   ```
4. **Copy the 3-D output to the Mac** (needs Chunhui Zhao's words). The
   comparison needs only these files from `l2_3d_<job>`:
   - `*_station_*.dat` without the `_copies` files, about 2.3 GB, and
     `*_station_conventions.txt`;
   - `*_rupture_time_r*.csv`;
   - `*_fault_r*.csv`, about 1.8 GB (the script reads only the last dump);
   - `RUN_INFO.txt`.
   The sizes are inferred from the cloud smoke's bytes per row. The run
   writes about 7 GB in all on Vista.
5. **Compare** (from `miniapps/seas`):
   ```bash
   R=/Users/chunhuizhao/projects/seas-mfem-2d-inplane-results/level2
   python3 elongated2p5d/scripts/compare_2p5d_3d.py --run3d $R/l2_3d_<job> \
       --run2d wa2020=$R/elongated_2d_wa2020_<commit>:1.17806 \
       --run2d lehner=$R/elongated_2d_lehner_<commit>:1.27324 \
       --run2d liang=$R/elongated_2d_liang_<commit>:1.57080 \
       --run2d "plane strain=$R/elongated_2d_plane_strain_<commit>" \
       --out $R/compare_<commit>
   ```
   - The 2.5D runs PASS when slip and arrival are within tolerance at every
     station in the window.
   - The plane-strain control is listed but not judged. It shows how far the
     map view is off without the term.

## Checks done before the runs (cloud, 2026-10-30)

- **`make test-2d5-dyn-driver`** (`tests/2d_inplane/nitsche_finite_width_decks.sh`).
  - Refusals:
    - zero W;
    - calibration and c together;
    - an unknown calibration or key;
    - LTS;
    - the energy balance;
    - a 3-D mesh.
  - Accepted (each runs 0.02 s with the expected c and beta): wa2020, lehner,
    and c as a number.
  - The term changes the stations.
- **Without `[finite_width]` the driver is unchanged bit for bit:**
  - 3-D TPV205: 66 files (`tests/2d_inplane/parity_3d.sh`);
  - TPV205-2D, 0.5 s: 20 files.
- **`compare_2p5d_3d.py` on synthetic input.** The 3-D input was the
  map-view smoke station copied to every depth.
  - Station depth average: 14.5 / 15 of the 2-D slip, as expected.
  - Dump average: 0.9833, as expected.
  - Implied c: 1.219, as expected.
  - The PDF has 6 pages.
