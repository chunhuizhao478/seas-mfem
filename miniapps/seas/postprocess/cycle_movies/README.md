# Standard videos and magnitude catalog of a fully dynamic cycle run

Chunhui Zhao, 2026-11-16: "let's develop this, later this could be a
standardlized video generation process, and I want the branch to save this as
well".

For runs of `seas_fd_cycle_driver` on a 2-D map-view mesh (triangles, fault
faces = line cells).  Everything here only reads the run folder.

## Workflow

1. Once per cluster: a Python environment with numpy, h5py, matplotlib,
   imageio-ffmpeg; for the line map, `fetch_natural_earth.sh FOLDER` on a
   login node (38.5 MB, public domain).
2. Push the commit, pull it on the cluster (CLAUDE.md, GitHub first).
3. Submit the job (Lonestar6, 1 development node, at most 1 SU):

   ```bash
   cd miniapps/seas/jobs/postprocess
   MOVIE_RUN=$SCRATCH/<runs>/<run folder> MOVIE_VENV=$SCRATCH/venv-post-ls6 \
   MOVIE_MAP_DATA=$SCRATCH/natural_earth \
   MAG_MU=32038120320 MAG_WIDTH_M=15000 \
      sbatch -A <allocation> ls6_cycle_movies.sbatch
   ```

4. Results in `$SCRATCH/cycle_movies/<run folder>_<job id>/`:

   | file | content |
   |---|---|
   | `particle_velocity.mp4` | particle speed, linear from 0, fault traces on top |
   | `fault_slip.mp4` | slip rate (linear from 0) beside cumulative slip, on the fault traces |
   | `particle_velocity_map.mp4`, `fault_slip_map.mp4` | the same on the line map (with `MOVIE_MAP_DATA`) |
   | `magnitude_catalog.csv`, `.png` | Mw of every catalog row (with `MAG_MU`, `MAG_WIDTH_M`) |
   | `RENDER_INFO.txt`, `DONE` | commit, run folder, options; mark of exit 0 |

By hand (a look at one frame before a long render, few workers on a login node):

```bash
python3 cycle_movies.py --run RUN --out OUT --frame-png 100
python3 cycle_movies.py --run RUN --out OUT --workers 4 --last-phase 3
```

## What a video shows

- The dynamic phase of every finished event, back to back.
- One frame per 4 s of simulated time (`--dt`), 24 frames per second (`--fps`).
- Label: phase number, model year of the switch, seconds since the switch.
- Fixed linear colour scales: `--v-max` 2 m/s, `--rate-max` 10 m/s,
  `--slip-max` = largest slip at the last frame.
- The last phase of the file is taken only when the run folder has `DONE`
  (or with `--include-last`).

## Why it is fast

- Frames are cut into segments of 120 (`--segment-frames`).
- `--workers` processes render segments at the same time, each into its own
  H.264 file (1-thread ffmpeg).
- The segment files are joined without re-encoding; the frame count of the
  result is checked.
- Per frame only the layers that change are drawn (`Blit`); `--full-redraw`
  gives the same frames the slow way.
- The pixel-to-triangle map of the volume video is computed once
  (`raster.locate`) and shared by all workers.
- The result does not depend on the number of workers (same bytes; test 4).

## Files

| file | job | in | out |
|---|---|---|---|
| [`cycle_movies.py`](cycle_movies.py) | the videos | run folder, options | MP4 files |
| [`magnitude_catalog.py`](magnitude_catalog.py) | Mw per catalog row | run folder, mu, W | CSV, PNG |
| [`cycle_run.py`](cycle_run.py) | readers: snapshot series, switch log, catalog, frames of each phase | run folder | arrays |
| [`raster.py`](raster.py) | which triangle holds each pixel centre, and its weights | points, triangles, grid | triangle id, 3 weights per pixel |
| [`basemap.py`](basemap.py) | line map in model x, y | region file, Natural Earth folder | drawn layers |
| [`regions/socal_ucerf2d.json`](regions/socal_ucerf2d.json) | where the UCERF 2-D model sits on the Earth, window, cities | | |
| [`fetch_natural_earth.sh`](fetch_natural_earth.sh) | downloads the 4 map files once | folder | 4 GeoJSON files |
| [`../../jobs/postprocess/ls6_cycle_movies.sbatch`](../../jobs/postprocess/ls6_cycle_movies.sbatch) | the Lonestar6 job | environment variables | results folder |
| [`../../tests/postprocess/test_cycle_movies.py`](../../tests/postprocess/test_cycle_movies.py) | checks on a small made-up run | | pass / fail |

## Pseudocode

`cycle_movies.py`

```
frames = dynamic snapshots of every finished phase, every n-th (n = dt / file cadence)
cut frames into segments
in parallel, for each segment:
    for each frame: read fields -> update the changing layers -> draw -> pipe RGB to ffmpeg
join the segment files (no re-encoding); count the frames; optional parts of --max-mb
```

Volume frame: `v(pixel) = sum of 3 weights x v(triangle node)`, `|v| = hypot(vx, vy)`.
Fault frame: per cell, mean of its 2 points of `hypot(strike, dip)` for slip rate and slip.

`magnitude_catalog.py`

```
for each catalog row (t_start, t_end):
    slip = |slip(first snapshot at or after t_end) - slip(start snapshot)|
    M0   = sum over fault cells of  mu x W x cell length x mean slip of its 2 points
    Mw   = (2/3) (log10 M0 - 9.1)
```

- W comes from outside (a 2-D run has no rupture width; for a 2.5-D run the
  deck's `[cycle.finite_width] W_m`).
- For an event that opens a phase the start snapshot is the one at the
  switch; slip before the switch is not counted (CSV column `start_gap_s`).

## Line map

- Model x, y (km) = azimuthal equidistant projection about `lat0`, `lon0`
  on a sphere of `radius_km` (region file).
- Layers: coastline, state lines, country border, lakes, graticule, cities.

## Limits

- 2-D map-view runs only (triangles in the volume file, 2-point cells in the
  fault file); anything else stops with a message.
- One shear modulus for the whole fault in the magnitudes.
- Colour scales are linear; there is no log option.

## Test

```bash
python3 miniapps/seas/tests/postprocess/test_cycle_movies.py
```
