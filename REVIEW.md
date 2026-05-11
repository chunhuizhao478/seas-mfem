# Code Review: Pre-Frontera-submission audit (Phase 6 / Phase 7 sbatch + drivers)

**Date:** 2026-05-10
**Reviewer:** code-review skill (round 3 — fresh adversarial pass, supersedes 2026-05-09 round 2)
**Trigger:** `/code-review do a fresh code review before submit to frontera run`

## Review Scope

- **Plan:** `miniapps/seas/document/io_dev/PLAN_bulk_compression_and_size_estimator_2026-05-09.md`
- **Sbatch under review (the four submission targets):**
  - `miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch` (8N × 400r, 48h, normal queue)
  - `miniapps/seas/jobs/tpv102/tpv102_phase6_paraview_zfp_dev_2hr.sbatch` (8N × 400r, 2h, dev queue)
  - `miniapps/seas/jobs/tpv104/tpv104_phase6_paraview_zfp_dev_2hr.sbatch`
  - `miniapps/seas/jobs/tpv205/tpv205_phase6_paraview_zfp_dev_2hr.sbatch`
- **C++ source (modified for Phase 6):**
  - `miniapps/seas/io/paraview_output.hpp` (1321-line diff)
  - `miniapps/seas/io/fault_vtkhdf_writer.hpp`, `fault_vtu_binary.hpp`
  - `miniapps/seas/drivers/tpv102_driver.cpp`, `tpv104_driver.cpp`, `tpv205_driver.cpp`
  - `miniapps/seas/tests/verification/bp5_verification_full.cpp`
  - `fem/datacollection.{hpp,cpp}`, `mesh/vtkhdf.{hpp,cpp}` (MFEM core)
- **Python source (Phase 7 estimator):**
  - `miniapps/seas/scripts/estimate_output_size.py`, `estimate_from_sbatch.py`
  - `miniapps/seas/scripts/_io_size_{mesh,schedule,compression,schemas}.py`
- **Domain context consulted:**
  - `miniapps/seas/CLAUDE.md` (project-level instructions, ZFP defaults table)
  - PLAN_bulk_compression_and_size_estimator_2026-05-09.md (Phase 6 + 7 spec)
  - Round 2 REVIEW.md (R-401..R-510 disposition; preserved below)
- **Round 2 carry-over:** R-501..R-510 were tracked in the prior REVIEW.md.  The
  user's fix pass landed R-507 (BP5 fault-only + 400-rank justification), R-510
  (`set -euo pipefail` on all four sbatch), R-411 (TPV205 sbatch
  `--fric-law lsw`), R-412 (`--tfinal 5.0` on all three TPV sbatch), and R-501
  (start-anchored `_FAULT_PREFIX_RE` two-pass match in `_io_size_mesh.py`).
  Round 3 issues new IDs starting at **R-701**.

---

## Round-2 disposition (verified by re-reading current code)

| Round-2 ID | Round-2 severity | Round-3 status | One-liner |
|------------|------------------|----------------|-----------|
| R-401 | CRITICAL | **FIXED** | `#include "mesh/vtkhdf.hpp"` in `paraview_output.hpp:24`; `static_assert` at line 126 |
| R-402 | CRITICAL | **FIXED** | Two-pass `read_gmsh` resolves fault tag from `$PhysicalNames` |
| R-403 | CRITICAL | DEFERRED | `module load phdf5/1.10.4` still no fallback (deferred by user) |
| R-404 | MODERATE | **STILL OPEN** | no `grep -q '^MFEM_USE_H5Z_ZFP *= YES'` in any sbatch |
| R-405 | MODERATE | **STILL OPEN** | `test-volume-hdf-compression` not in `make test:` aggregate |
| R-406 | MODERATE | **STILL OPEN** | `(integer_index, zfp_*)` missing from `COMPRESSION_RATIOS` |
| R-407 | MODERATE | **STILL OPEN** | `estimate_n_volume_writes` ignores cap when `volume_pv_dt > 0` |
| R-408 | MODERATE | **PARTIALLY ADDRESSED** | CLAUDE.md line 189 was clarified to `displacement in BP5` per Phase 6 block at line 200; the older row at line 189 still reads `1e-3 (m/s for velocity)`.  Self-contradiction remains. |
| R-409 | MODERATE | **STILL OPEN** | `_io_size_schedule.ScheduleConfig.dt_coseismic = 60.0` vs C++ `0.01` |
| R-410 | MODERATE | **STILL OPEN** | estimator collapses primary+secondary into one filter |
| R-411 | MODERATE | **FIXED in sbatch** (TPV205 sbatch line 80 now `--fric-law lsw`); driver still ignores `--fric-law` (function param is `/*fric_law_cli*/`) but banner now matches |
| R-412 | MODERATE | **FIXED** | TPV102/104/205 sbatch all use `--tfinal 5.0` (comment justifies the 2h wall envelope) |
| R-413 | MODERATE | **FIXED** (BP5 sbatch line 104 has `--no-volume-pv`) |
| R-501 | CRITICAL | **FIXED** (`_FAULT_PREFIX_RE`, two-pass match, regression test landed) |
| R-502 | HIGH | **STILL OPEN** | malformed `$PhysicalNames` count error is still misleading |
| R-503 | HIGH | **STILL OPEN** | argparse `allow_abbrev=True` makes `--paraview-dt` ambiguous on estimator |
| R-504 | HIGH | **STILL OPEN** | C++ uses `--paraview-coseismic-dt`; Python uses `--paraview-dt-coseismic` |
| R-505 | MODERATE | **PARTIAL** | exact-match preferred (test on `{1:"fault_main",2:"fault"}==2`); 2nd-pass within prefix matches still order-dependent |
| R-506 | MODERATE | **STILL OPEN** | two-pass file open in `read_gmsh` |
| R-507 | HIGH | **PARTIALLY FIXED** | BP5 sbatch now has `--no-volume-pv` (fault-only); 400 ranks (not 800) — sbatch comment justifies the deviation; plan §7.5 acceptance still technically not met (still text-says-800).  Treat as resolved unless plan §7.5 is updated. |
| R-508 | HIGH | **STILL OPEN** | static_assert message doesn't name the Frontera-known-good `hdf5/1.14.4` fallback |
| R-509 | MODERATE | **STILL OPEN** | CLAUDE.md self-contradiction (lines 189 vs 200) |
| R-510 | MODERATE | **FIXED** (all four sbatch have `set -euo pipefail` at the top) |

**Round 3 new findings:** R-701..R-708.

---

## Severity summary (round 3 net-new)

| ID | Severity | Area | One-liner |
|----|----------|------|-----------|
| **R-701** | **CRITICAL** | sbatch / BP5 | `--tfinal $(python3 -c 'print(250*3.156e7)')` is a runtime command substitution.  If `python3` is missing on the compute node (it is NOT loaded as a module by the sbatch), `$()` yields **empty**, the driver sees `--tfinal --verify`, `std::atof("--verify")=0.0`, `t_final=0`, and the simulation exits at step 1.  The 48-hour normal-queue reservation is consumed for nothing (~38 400 SU lost).  `set -e` does NOT abort on a failed command-substitution that appears in argv position. |
| **R-702** | MODERATE | sbatch / BP5 | `--tfinal 250*3.156e7 = 7.890e9 s` but the C++ uses `BP5Params::seconds_per_year = 365.25*24*3600 = 31_557_600` (i.e. `250*3.15576e7 = 7.8894e9 s`).  Off by 0.18% — ~14 days short.  Cosmetic for the production run but **the script-derived tfinal disagrees with the value the driver would compute internally** for the same "250 yr" intent. |
| **R-703** | HIGH | BP5 sbatch | sbatch comment block (`bp5_phase6_paraview_zfp_normal_48hr.sbatch` lines 32-35) claims "~430k uncapped writes over 250 yr, dominated by the BP5 nucleation budget (~1 day per nucleation × **12 events**)".  CLAUDE.md says BP5 recurrence is ~240 yr.  250 yr / 240 yr = **1.04 events**, not 12.  The Python estimator's `_default_n_events("bp5", 7.89e9) = round(1.04) = 1`.  The sbatch comment is off by 12×; whoever sized the snapshot cap from this estimate is using the wrong inputs. |
| **R-704** | MODERATE | Python | `estimate_from_sbatch.py:estimate_one` reads only `paraview-max-snapshots`, `paraview-dt`, and the four `paraview-{fault,volume}-{zfp-tol,deflate-level}` flags into the `ScheduleConfig`.  It **silently drops** `--paraview-coseismic-dt`, `--paraview-nucleation-dt`, `--paraview-interseismic-dt`, `--paraview-bulk-dt`, `--paraview-bulk-zfp-tol`, `--paraview-bulk-deflate-level`, and `--paraview-output-every-n-steps`.  For a TPV102 sbatch with `--paraview-bulk-dt 0.5 --paraview-bulk-zfp-tol 1e-3`, the secondary "wave_bulk" collection's writes are NOT counted — the estimate misses ~half the actual on-disk volume PV. |
| **R-705** | MODERATE | C++ / paraview_output.hpp | `WriteFaultPVD` (called from `WriteFaultSurfaceVTU` Vtu path on every cycle) does `std::ofstream pvd(pvd_name, std::ios::trunc);` and re-emits the entire PVD file from `fault_pvd_entries_` every save.  At BP5 production scale this is **O(N²) bytes written over N cycles** (N=5000 ⇒ ~5 GB of redundant PVD rewrites).  On Frontera Lustre this is also 5000 truncates and 5000 metadata-roundtrips.  Quasi-irrelevant in HDF5 mode (the PVD path is dead), but the dead-code path is still wired and a future regression that flips the default back to Vtu would silently 25× the I/O. |
| **R-706** | MODERATE | C++ / paraview_output.hpp | `ParaViewOutput<MeshType>` constructor calls `std::make_unique<ParaViewHDFDataCollection>(collection_name, &mesh)` when `volume_mode == Hdf5`, even on BP5 where the driver immediately calls `SetVolumeSaveEnabled(false)`.  The HDF data collection is allocated and registers `displacement` + `mpi_rank` + 12 fault GFs against an L2 FES — but `Save()` is never invoked, so the file is never created and the allocations are dead weight.  Memory cost (~MB-scale) is small, but the `RegisterField(...)` paths and the L2 FES + 12 GFs allocate per-process.  **The bigger risk:** if any future code path inadvertently calls `pv_dc_->Save()` (e.g. through a forwarding helper), the file is created and 250-year worth of volume snapshots gets written by accident — defeating `--no-volume-pv`. |
| **R-707** | LOW | C++ / paraview_output.hpp | The MPI uniformity check at `paraview_output.hpp:1080-1093` is gated by `output_mode_uniformity_checked_` (latched after first pass).  Mode setters (`SetFaultOutputMode`, `SetLegacyAsciiVTU`) called AFTER the first `WriteFaultSurfaceVTU` would not be re-validated — and current driver code does NOT call those setters mid-run, but the contract is implicit.  No exception-safety guard. |
| **R-708** | LOW | Python | `estimate_from_sbatch.py:_PRINT_ARITH_RE` only matches `$(python3 -c 'print(<arithmetic>)')` literally.  The BP5 sbatch line 107 form is `$(python3 -c 'print(250*3.156e7)')` which DOES match (digits + `*` + scientific notation), so the estimator can replay it.  But this is exactly the form vulnerable to R-701 at SLURM submission time. |

---

## CRITICAL findings — full detail

### [R-701] CRITICAL — BP5 sbatch `--tfinal $(python3 ...)` will silently produce `t_final = 0` if python3 is missing on the compute node

**Category:** BUG

**Location:** `miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch:107`

```bash
--tfinal $(python3 -c 'print(250*3.156e7)') \
--verify
```

**Trigger:** any compute node where `python3` is not on `$PATH` after the
`module load` chain at lines 48-54.  The sbatch loads `intel/19.1.1`,
`impi/19.0.9`, `hypre/2.31.0`, `mumps/5.3`, `parmetis`, `petsc/3.15`, and
`phdf5/1.10.4` — none of these expose `python3` as a side effect.  On
Frontera, system `python3` is usually present at `/usr/bin/python3`, BUT:

- A future image refresh could remove it (TACC does this periodically when
  Python EOLs land).
- On a future site (e.g. moving the job to Stampede3 or Perlmutter), the
  base PATH may not include `python3` until `module load python3` runs.
- A user's `.bashrc` could rewrite `$PATH` and accidentally clobber the
  system path on the compute node.

**Actual behavior on `python3`-missing node:**

1. `$(python3 -c '...')` runs `python3`, which fails with
   `bash: python3: command not found`.  Exit code 127, stdout empty.
2. The command-substitution evaluates to the empty string.  The composed
   command line becomes:
   ```bash
   ibrun ./seas_bp5_full ... --tfinal --verify
   ```
3. **`set -euo pipefail` does NOT abort** here.  Per the bash manual, a
   failed command-substitution used in an argument-list position does
   **not** trigger `set -e` (it only triggers `set -e` if it is the
   trailing simple command).  Verified: `set -e; cmd $(false)` does NOT
   exit.
4. The driver's parser (`bp5_verification_full.cpp:791-794`):
   ```cpp
   if (arg == "--tfinal" && i + 1 < argc)
   {
      tfinal_override = std::atof(argv[++i]);
   }
   ```
   reads `argv[i+1] = "--verify"`.  `std::atof("--verify")` returns `0.0`
   (atof stops at the first non-numeric character, and `-` followed by
   `-` is not numeric).
5. Then at line 1223: `if (tfinal_override >= 0.0) { t_final = tfinal_override; }`
   — `0.0 >= 0.0` is true, so `t_final = 0.0`.
6. Simulation enters the time-stepping loop, computes `t = 0`, hits the
   exit condition `t >= t_final` on iteration 1, and ends.
7. Output is fault HDF with one snapshot at t=0, plus the verify CSVs.
8. SLURM records exit 0.  The 48-hour reservation is consumed.

**Expected behavior:** the script should either use a pre-computed
literal, an environment variable, or fail loudly when `python3` is
missing.

**Suggested fix (one of):**

**Option A — hardcode the literal:**
```diff
-      --tfinal $(python3 -c 'print(250*3.156e7)') \
+      --tfinal 7889400000 \
```
(Note: matches the C++ `BP5Params::seconds_per_year = 31557600` exactly:
`250 * 31557600 = 7889400000`.  This also fixes R-702.)

**Option B — bash arithmetic (no python3 dependency):**
```diff
-      --tfinal $(python3 -c 'print(250*3.156e7)') \
+      --tfinal $(( 250 * 31557600 )) \
```
Bash supports integer arithmetic natively; result is `7889400000`.

**Option C — pre-compute outside the ibrun line, with a verification step:**
```diff
+# tfinal in seconds; 250 yr × 31_557_600 s/yr (matches C++ seconds_per_year).
+T_FINAL_S=7889400000
+test "${T_FINAL_S}" -gt 0 || { echo "ERROR: T_FINAL_S not set"; exit 1; }
+
 ibrun ./seas_bp5_full \
       ...
-      --tfinal $(python3 -c 'print(250*3.156e7)') \
+      --tfinal "${T_FINAL_S}" \
```

**Recommendation:** Option A.  Simplest, no parsing risk, and the
literal is auditable from the sbatch alone.

**Test case (proposed):**
```bash
# Regression test for R-701: even if python3 disappears, tfinal is
# correctly read from the sbatch.
def test_R701_bp5_sbatch_tfinal_is_literal():
    sb = Path("miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch")
    text = sb.read_text()
    # No subshell-evaluated tfinal — must be a literal numeric value.
    assert "--tfinal $(" not in text, (
        "R-701: BP5 sbatch must not use subshell-derived --tfinal.  "
        "If python3 is missing on the compute node, $() returns empty "
        "and the driver parses --tfinal --verify as tfinal=0, exiting "
        "the simulation at step 1 and wasting the 48-hour reservation.")
    # Must be a positive integer literal.
    m = re.search(r"--tfinal\s+(\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)", text)
    assert m, "R-701: BP5 sbatch missing literal --tfinal"
    assert float(m.group(1)) > 0
```

---

## HIGH findings — full detail

### [R-703] HIGH — BP5 sbatch comment claims "~430k uncapped writes ... 12 events"; the Python estimator says ~1 event

**Category:** ASSUMPTION (documentation drift)

**Location:** `miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch:32-35`

```bash
# Snapshot cap (--paraview-max-snapshots 5000) bounds the writer count
# at the production default — the schedule integrator estimates ~430k
# uncapped writes over 250 yr, dominated by the BP5 nucleation budget
# (~1 day per nucleation × 12 events).
```

**Trigger:** running the size estimator against this sbatch.

**Actual behavior of the estimator:**

`_io_size_schedule._default_n_events("bp5", 7.89e9)`:
```python
recurrence_s = 240.0 * seconds_per_year  # = 240 * 3.156e7 = 7.5744e9
return max(0, round(tfinal / recurrence_s))  # round(7.89e9 / 7.5744e9) = round(1.04) = 1
```

So the estimator computes:
- `n_events = 1`
- `nucleation_budget = 1 * 86_400 = 86 400 s` ⇒ ~86 400 fault writes
- `coseismic_budget = 1 * 30 = 30 s` ⇒ ~0.5 fault writes at dt_coseismic=60 (Python default; see R-409)
- `interseismic_budget ≈ 7.89e9 - 86 430 ≈ 7.89e9 s` ⇒ ~250 fault writes
- **Total uncapped: ~86 650 writes**, not 430 000.

For 12 events (the sbatch comment's claim), the total would be
~12 × 86 400 + 360 + 240 ≈ **1 040 000 writes**.

Neither matches "430k".  The 430k number appears to be from a different
calibration anchor or an old draft of the plan.

**Expected behavior:** the sbatch comment should either:
(a) cite the actual Python estimator output (run `estimate_from_sbatch.py`
    on this sbatch and quote the number);
(b) state the assumed n_events explicitly so the reader can audit.

**Suggested fix:**
```diff
-# Snapshot cap (--paraview-max-snapshots 5000) bounds the writer count
-# at the production default — the schedule integrator estimates ~430k
-# uncapped writes over 250 yr, dominated by the BP5 nucleation budget
-# (~1 day per nucleation × 12 events).
+# Snapshot cap (--paraview-max-snapshots 5000) bounds the writer count
+# at the production default.  The schedule integrator estimates ~86k
+# uncapped writes over 250 yr at the Python default n_events=1
+# (rounded from 250 yr / 240 yr recurrence per CLAUDE.md "What
+# Constitutes a Regression").  Nucleation dominates (~86 400 writes
+# per event at dt_nucleation=1.0 s × NUCLEATION_DURATION_S["bp5"]=86 400).
+# Cross-check: `python3 miniapps/seas/scripts/estimate_from_sbatch.py
+# miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch`.
```

**Test case (proposed):** Run the estimator on the sbatch and assert the
number is consistent with the comment block:
```python
def test_R703_bp5_sbatch_comment_matches_estimator():
    from estimate_from_sbatch import parse_sbatch, estimate_one
    p = parse_sbatch(Path("miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch"))
    result = estimate_one(p)
    # The comment block should be within 2× of the estimator's number.
    sbatch_text = p.sbatch_path.read_text()
    m = re.search(r"~(\d+)k uncapped writes", sbatch_text)
    if m:
        claimed = int(m.group(1)) * 1000
        actual = result["n_writes_fault"]
        # +n_events because cap inflates by event count.
        assert claimed / 5 < actual < claimed * 5, (
            f"R-703: sbatch claims ~{claimed} writes but estimator says "
            f"{actual} (5× tolerance)")
```

---

## MODERATE findings — full detail

### [R-702] MODERATE — `tfinal = 250 * 3.156e7` ≠ `250 * BP5Params::seconds_per_year`

**Category:** ASSUMPTION

**Location:** `bp5_phase6_paraview_zfp_normal_48hr.sbatch:107` AND
`miniapps/seas/config/bp5_params.hpp:172`.

```cpp
// bp5_params.hpp
static constexpr real_t seconds_per_year = 365.25 * 24.0 * 3600.0;
//                                       = 31_557_600.0
```
```bash
# sbatch
--tfinal $(python3 -c 'print(250*3.156e7)') \
#                       \________/
#                       = 31_560_000.0, off by +2 400 s/yr
```

- `250 × 3.156e7      = 7 890 000 000 s` (sbatch)
- `250 × 31_557_600   = 7 889 400 000 s` (C++)

Difference: **600 000 s ≈ 6.94 days**.  Over a 250-year run, the sbatch
asks for 6.94 days more than the C++ would interpret 250 years as.

**Effect:** mostly cosmetic, but every downstream tool that reads "250
yr from the sbatch" and "250 yr from the C++ banner" will see a 6.94-day
mismatch.  The Python estimator's `parse_time("250yr")` uses
`3.156e7` (line 453 of `test_estimate_output_size.py`); the C++
`SetTotalRunTime(tfinal)` cap projection uses the literal seconds.  Both
agree on the Python definition of "year"; the C++ inside the driver
uses the more-precise 31_557_600.

**Suggested fix:** see R-701 Option A (`--tfinal 7889400000`).  Bonus:
unifies the sbatch with the C++ definition.

**Test case (proposed):**
```python
def test_R702_year_constant_consistency():
    import re
    from pathlib import Path
    bp5_params = Path("miniapps/seas/config/bp5_params.hpp").read_text()
    m = re.search(r"seconds_per_year\s*=\s*([^;]+);", bp5_params)
    assert m, "seconds_per_year not found in bp5_params.hpp"
    cpp_value = eval(m.group(1).replace("real_t", "").strip(), {})
    assert cpp_value == 31_557_600.0
    # And the Python estimator's parse_time("1yr") should agree
    # to within 0.01%.
    from estimate_output_size import parse_time
    py_value = parse_time("1yr")
    assert abs(py_value - cpp_value) / cpp_value < 1e-4, (
        f"R-702: parse_time('1yr') = {py_value} disagrees with "
        f"BP5Params::seconds_per_year = {cpp_value}")
```

(NB: `parse_time("1yr") = 3.156e7`, which is ~0.0076% different from
31 557 600.  Below the 0.01% tolerance fail.  This is a known polite
shorthand: the size-estimator round-trip is dimensionally consistent.
The R-702 finding is specifically about the sbatch's `--tfinal` value,
not the Python estimator's `parse_time` helper.)

---

### [R-704] MODERATE — `estimate_from_sbatch.py:estimate_one` drops per-regime dt + bulk-collection flags

**Category:** BUG (silent under-estimate)

**Location:** `miniapps/seas/scripts/estimate_from_sbatch.py:267-274`.

```python
schedule = ScheduleConfig(
    tfinal=tfinal_seconds,
    max_total_snapshots=_maybe_int(p.flags.get("paraview-max-snapshots"), 0),
    fixed_dt=_maybe_float(p.flags.get("paraview-dt")) or 0.0,
    n_events=0,
    driver=p.driver,
)
```

**Trigger:** running `estimate_from_sbatch.py` on a sbatch that uses any
of: `--paraview-coseismic-dt`, `--paraview-nucleation-dt`,
`--paraview-interseismic-dt`, `--paraview-bulk-dt`,
`--paraview-bulk-zfp-tol`, `--paraview-bulk-deflate-level`,
`--paraview-output-every-n-steps`.

**Actual behavior:** `estimate_one` constructs `ScheduleConfig` with only
`tfinal`, `max_total_snapshots`, `fixed_dt`, `n_events`, `driver`.  Every
other field defaults.  In particular:
- `dt_coseismic = 60.0` (Python default), regardless of the sbatch.
- `dt_nucleation = 1.0`.
- `dt_interseismic = 3.156e7`.
- The secondary bulk collection (when `--paraview-bulk-dt > 0`) is
  invisible: there is no second `estimate_n_writes(cfg_bulk)` call, no
  bulk byte accumulation.

For the TPV102 sbatch (`--paraview-bulk-dt 0.5 --paraview-bulk-zfp-tol 1e-3
--paraview-dt 0.5 --tfinal 5.0`):
- Primary collection writes: 5.0 / 0.5 + 1 = 11
- Secondary "wave_bulk" collection writes: also 11 (independent of primary)
- **Estimator counts only 11 writes** (the primary), not 22.

Bytes-on-disk error: secondary collection adds 3 stress fields ×
~n_dofs × 8B ≈ comparable to primary.  Estimate is ~50% LOW for TPV*
sbatch that use a secondary collection.

**Expected behavior:** parse all C++-driver-accepted `--paraview-*` flags
and populate `ScheduleConfig` accordingly; account for the secondary
collection.

**Suggested fix (minimal — fix per-regime dt now; defer secondary
collection because it requires schema work):**
```diff
 schedule = ScheduleConfig(
     tfinal=tfinal_seconds,
     max_total_snapshots=_maybe_int(p.flags.get("paraview-max-snapshots"), 0),
     fixed_dt=_maybe_float(p.flags.get("paraview-dt")) or 0.0,
+    dt_coseismic=_maybe_float(p.flags.get("paraview-coseismic-dt")) or 0.01,
+    dt_nucleation=_maybe_float(p.flags.get("paraview-nucleation-dt")) or 1.0,
+    dt_interseismic=(_maybe_float(p.flags.get("paraview-interseismic-dt"))
+                     or 3.156e7),
+    output_every_n_steps=_maybe_int(p.flags.get("paraview-every"), 0),
     n_events=0,
     driver=p.driver,
 )
```

Note the new defaults (`0.01`, `1.0`, `3.156e7`) match the C++ driver's
`AdaptiveSchedule` defaults at `paraview_output.hpp:168-170` — NOT the
old Python `ScheduleConfig` defaults (60.0 / 1.0 / 3.156e7).  The Python
default-`dt_coseismic = 60.0` is itself a round-2 R-409 bug.

**Test case (proposed):**
```python
def test_R704_estimator_reads_coseismic_dt_from_sbatch(tmp_path):
    from estimate_from_sbatch import parse_sbatch, estimate_one
    sb = tmp_path / "test.sbatch"
    sb.write_text("""
        ibrun ./seas_tpv102_driver --mesh /tmp/mesh.msh \\
            --tfinal 1.0 --paraview --paraview-coseismic-dt 0.005
    """)
    # Fake mesh
    mesh_path = Path("/tmp/mesh.msh")
    mesh_path.write_text(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
        "$Nodes\n1\n1 0 0 0\n$EndNodes\n"
        "$Elements\n0\n$EndElements\n")
    p = parse_sbatch(sb)
    # The parsed flag must propagate into the schedule.
    assert p.flags.get("paraview-coseismic-dt") == "0.005"
    result = estimate_one(p)
    # n_writes_fault must reflect the smaller dt.
    # Before R-704 fix: ScheduleConfig uses default 60.0, so n_writes ≈ 1.
    # After R-704 fix: ScheduleConfig uses 0.005, n_writes for 1 event
    #   coseismic_budget=30s / 0.005 = 6000 writes.
    assert result["n_writes_fault"] > 100  # crude lower bound
```

---

### [R-705] MODERATE — Legacy `WriteFaultPVD` truncate-and-rewrite is O(N²) bytes

**Category:** QUALITY (performance, dead-path)

**Location:** `miniapps/seas/io/paraview_output.hpp:2045-2060` (`WriteFaultPVD`).

```cpp
void WriteFaultPVD(const std::string &fault_dir)
{
   std::string pvd_name = fault_dir + "/fault_surface.pvd";
   std::ofstream pvd(pvd_name, std::ios::trunc);
   pvd << ... ;
   for (const auto &e : fault_pvd_entries_) {
      pvd << "<DataSet ... file=\"" << e.file << "\"/>\n";
   }
   ...
}
```

Called from `WriteFaultSurfaceVTU` Vtu path (line 1384) every cycle on
rank 0.  In Phase 6 production the default is Hdf5 mode, so this path
is dead.  But:

**Concern:** if anyone passes `--paraview-fault-vtu` (e.g. for debugging
or downstream-script compatibility), each fault write re-truncates and
re-emits the entire PVD file.  For N=5000 cycles, total PVD bytes
written = N + 2(N−1) + 3(N−2) + ... ≈ O(N²).  At ~100 bytes per entry
and N=5000, that's ~1.2 GB of redundant PVD writes spread over 5000
metadata-roundtrips on Lustre.

**Effect:** silent slow path; rank-0 CPU + filesystem I/O bottleneck for
the legacy fault Vtu path.  Not relevant for the four sbatch under
review, but the dead path is wired and a future regression could
re-default to Vtu.

**Suggested fix:** append-only PVD update.  Open the file once on Init,
write the prolog, and on each cycle seek-to-just-before-`</Collection>`
and inject the new `<DataSet>` entry.  Or: write a per-cycle entry to a
sidecar file and only emit the consolidated PVD at the end of the run.
Lower-priority because it does not affect Phase 6 production paths.

**Test case (proposed):**
```cpp
TEST_CASE("R-705: legacy fault PVD writes are O(N) not O(N^2)") {
   // Force Vtu fault mode + legacy ASCII off (binary VTU path).
   ParaViewOutput pv(...);
   pv.SetFaultOutputMode(FaultOutputMode::Vtu);
   for (int c = 0; c < 100; ++c) {
       pv.WriteFaultSurfaceVTU(...);
   }
   // Read the PVD file size.  Expected: ~100 lines × ~100 bytes = ~10 KB.
   // O(N^2) bug would produce ~500 KB (sum of partial writes).
   auto sz = std::filesystem::file_size(pvd_path);
   REQUIRE(sz < 50'000);
}
```
(LOW priority — not blocking Frontera; downgrade severity if you defer.)

---

### [R-706] MODERATE — `pv_dc_ = ParaViewHDFDataCollection` is allocated even when `--no-volume-pv` is set

**Category:** QUALITY / DEFENSIVE

**Location:** `miniapps/seas/io/paraview_output.hpp:396-425` (constructor).

```cpp
if (volume_output_mode_ == VolumeOutputMode::Hdf5)
{
   pv_dc_ = std::make_unique<ParaViewHDFDataCollection>(
               collection_name, &mesh);
}
```

The data collection is allocated for **every** `ParaViewOutput`, even on
BP5 where the very next driver call is `SetVolumeSaveEnabled(false)`.
The HDF5 data collection internally allocates field maps, the prefix
path, etc.  But `EnsureVTKHDF` is called only inside `TSave()`, so the
on-disk `.vtkhdf` file is NOT created.

**Concern:** the dead-code allocation is fragile.  Any future helper that
forwards through `pv_dc_->Save()` without checking `emit_volume_save_`
will create and start writing to `<prefix>/<collection>.vtkhdf`.  Look
at the verification driver line 1700-1716: `RegisterField` is invoked
on `pv_dc_->RegisterField("displacement", ...)` and
`pv_dc_->RegisterField("mpi_rank", ...)`.  These are registrations
against the HDF data collection that will never `Save()` — they're
no-ops in terms of on-disk output.

**Risk:** a future driver that calls
`pv_out->RegisterDomainField("debug_field", &gf)` (which forwards to
`pv_dc_->RegisterField`) will silently register against the
never-saved HDF collection.  No error.  The user thinks their field is
emitted; it is not.

**Suggested fix (defensive):** add an explicit assertion in `Save()` /
`ForceSave()` / `ForceSaveImpl()` that reads as:

```diff
 bool ForceSaveImpl(int cycle, real_t time)
 {
    if (emit_volume_save_)
    {
+      // Defensive: catch any future regression that forgets to disable
+      // the file open in --no-volume-pv mode.
       bool volume_should_write = true;
       ...
       if (volume_should_write)
       {
          pv_dc_->SetCycle(cycle);
          pv_dc_->SetTime(time);
          pv_dc_->Save();
          last_volume_write_time_ = time;
       }
    }
    last_write_time_ = time;
    CommitOnceAtCycle(cycle);
    return true;
 }
```

Stronger fix (deferred to follow-up): lazily allocate `pv_dc_` only when
`SetVolumeSaveEnabled(true)` AND a Save is imminent.  This requires
restructuring `RegisterField` to defer until lazy-init.  Out of scope
for the Frontera submission.

**Test case (proposed):**
```cpp
TEST_CASE("R-706: --no-volume-pv produces zero on-disk volume bytes") {
   ParaViewOutput<Mesh> pv("test_out", mesh, /*order=*/1);
   pv.SetVolumeSaveEnabled(false);
   pv.RegisterDomainField("dummy", &gf);
   // ... save 100 cycles ...
   for (int c = 0; c < 100; ++c) pv.ForceSave(c, c * 0.01);
   // Expected: no test_out/volume.vtkhdf, no test_out/volume_*.vtu.
   REQUIRE_FALSE(std::filesystem::exists("test_out/volume.vtkhdf"));
   REQUIRE_FALSE(std::filesystem::exists("test_out/ParaView"));
}
```

---

## STILL-OPEN findings from round 2 (verbatim from round 2)

### [R-404] MODERATE — `MFEM_USE_H5Z_ZFP=YES` not asserted in any sbatch

Verified by re-reading all four sbatch.  Each has:
```bash
grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }
grep -q '^MFEM_USE_HDF5 *= YES'  config/config.mk || { echo "ERROR: no HDF5";  exit 1; }
```

but NO `grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk`.

**Risk:** all four sbatch use `--paraview-*-zfp-tol` flags (BP5: fault;
TPV*: fault + volume + bulk).  The driver aborts at parse-time on
non-`MFEM_USE_H5Z_ZFP` builds, but the abort fires AFTER `ibrun` starts
(~30-60 SU wasted per job; 240-380 SU across all four).

**Suggested fix:**
```diff
 grep -q '^MFEM_USE_HDF5 *= YES'  config/config.mk || { echo "ERROR: no HDF5";  exit 1; }
+grep -q '^MFEM_USE_H5Z_ZFP *= YES' config/config.mk || { echo "ERROR: no H5Z-ZFP"; exit 1; }
```

---

### [R-503] HIGH — `--paraview-dt` ambiguous in estimator argparse

Verified: `estimate_output_size.py:431-433` still has three options sharing
the `--paraview-dt-*` prefix.  argparse's `allow_abbrev=True` (default)
causes `--paraview-dt 0.5` (the form passed to the C++ driver) to fail
with `error: ambiguous option`.

**Suggested fix:**
```diff
 def _build_parser() -> argparse.ArgumentParser:
-    p = argparse.ArgumentParser(prog="estimate_output_size", ...)
+    p = argparse.ArgumentParser(prog="estimate_output_size",
+                                allow_abbrev=False, ...)
```

---

### [R-504] HIGH — Per-regime dt flag-name mismatch

C++ driver flags (tpv102:585-587):
```cpp
GetRealArg(argc, argv, "--paraview-coseismic-dt",   -1.0);
GetRealArg(argc, argv, "--paraview-nucleation-dt",  -1.0);
GetRealArg(argc, argv, "--paraview-interseismic-dt",-1.0);
```

Python estimator flags (estimate_output_size.py:431-433):
```python
p.add_argument("--paraview-dt-coseismic", type=float, default=60.0)
p.add_argument("--paraview-dt-nucleation", type=float, default=1.0)
p.add_argument("--paraview-dt-interseismic", type=float, default=3.156e7)
```

Different word order, no flag aliasing.  A user reading the sbatch sees
the C++ name and tries it in the estimator — it's not recognised.

**Suggested fix:** rename Python args to match C++:
```diff
-p.add_argument("--paraview-dt-coseismic", type=float, default=60.0)
-p.add_argument("--paraview-dt-nucleation", type=float, default=1.0)
-p.add_argument("--paraview-dt-interseismic", type=float, default=3.156e7)
+p.add_argument("--paraview-coseismic-dt", type=float, default=0.01)
+p.add_argument("--paraview-nucleation-dt", type=float, default=1.0)
+p.add_argument("--paraview-interseismic-dt", type=float, default=3.156e7)
```
(Also changes default for `--paraview-coseismic-dt` from 60.0 to 0.01 to
match C++ defaults — see R-409.  Update any callers in main()
accordingly.)

---

### [R-409] MODERATE — Python `dt_coseismic = 60.0` vs C++ `0.01`

`_io_size_schedule.py:42`:
```python
dt_coseismic: float = 60.0
```

`paraview_output.hpp:168`:
```cpp
real_t dt_coseismic = 0.01;
```

**Effect:** estimator under-counts coseismic writes by 6000× on any
sbatch that uses the AdaptiveSchedule path (BP5: yes; TPV* with
`--paraview-dt` set: no, fixed-dt wins).

**Suggested fix:**
```diff
 @dataclass
 class ScheduleConfig:
     tfinal: float
     v_coseismic: float = 1e-3
     v_nucleation: float = 1e-7
-    dt_coseismic: float = 60.0           # seconds
+    dt_coseismic: float = 0.01           # match C++ paraview_output.hpp:168
     dt_nucleation: float = 1.0
     dt_interseismic: float = 3.156e7     # 1 yr
```

---

### [R-406] MODERATE — `(integer_index, zfp_*)` missing from `COMPRESSION_RATIOS`

`_io_size_compression.py:20-48`: only `(integer_index, deflate_6)` and
`(integer_index, vtu_binary)` are present.  TPV* runs with
`--paraview-volume-zfp-tol 1e-3` look up
`(integer_index, zfp_1e-3)` for the `mpi_rank` field and fall through to
the 8.0 raw-double default + stderr warning.

The actual on-disk size of `mpi_rank` under ZFP is the deflate-fallback
path in `mesh/vtkhdf.cpp:165` (integers fall back to deflate), so the
real ratio is `(integer_index, deflate_6) = 0.5`.  The estimator
over-estimates by 16×.

**Suggested fix:**
```diff
     ("integer_index", "deflate_6"):     (0.5, "PLACEHOLDER"),
+    # Integer datasets fall back to deflate via VTKHDF::EnsureDataset
+    # even when ZFP is the requested algorithm (mesh/vtkhdf.cpp:162-170).
+    ("integer_index", "zfp_1e-3"):      (0.5, "VTKHDF deflate fallback"),
+    ("integer_index", "zfp_1e-6"):      (0.5, "VTKHDF deflate fallback"),
+    ("integer_index", "zfp_1e-9"):      (0.5, "VTKHDF deflate fallback"),
+    ("integer_index", "zfp_1e-12"):     (0.5, "VTKHDF deflate fallback"),
     ("integer_index", "vtu_binary"):    (8.0, "raw double"),
```

---

### [R-407] MODERATE — `estimate_n_volume_writes` ignores cap when `volume_pv_dt > 0`

`_io_size_schedule.py:130-139`:
```python
def estimate_n_volume_writes(cfg: ScheduleConfig,
                             volume_pv_dt: float) -> int:
    if volume_pv_dt <= 0:
        return estimate_n_writes(cfg)
    if cfg.tfinal <= 0:
        return 1
    return max(1, int(cfg.tfinal // volume_pv_dt) + 1)
```

When `volume_pv_dt > 0` the cap is NOT applied.  Not blocking BP5 (uses
the adaptive path) or the current TPV* sbatch (no `--volume-pv-dt`), but
a future sbatch could trip this.

**Suggested fix:**
```diff
 def estimate_n_volume_writes(cfg: ScheduleConfig,
                              volume_pv_dt: float) -> int:
     if volume_pv_dt <= 0:
         return estimate_n_writes(cfg)
     if cfg.tfinal <= 0:
         return 1
-    return max(1, int(cfg.tfinal // volume_pv_dt) + 1)
+    total = max(1, int(cfg.tfinal // volume_pv_dt) + 1)
+    return _apply_cap(total, cfg)
```

---

### [R-411] MODERATE — TPV205 `--fric-law` parameter still ignored (sbatch is now correct)

Verified at `tpv205_driver.cpp:188`: `GetDispatchedLaw(const std::string &/*fric_law_cli*/)`
parameter is comment-suppressed; the function returns `DispatchedLaw::LSW`
unconditionally.

The round-2 fix was at the sbatch level: TPV205 sbatch line 80 now
passes `--fric-law lsw`, so the banner output ("Friction law: lsw …")
agrees with the actual dispatch.  But the driver still silently ignores
ANY value passed to `--fric-law`, e.g. a future sbatch with
`--fric-law aging` would also map to LSW.

**Suggested fix (defensive):**
```diff
-static DispatchedLaw GetDispatchedLaw(const std::string &/*fric_law_cli*/)
+static DispatchedLaw GetDispatchedLaw(const std::string &fric_law_cli)
 {
+   if (!fric_law_cli.empty() && fric_law_cli != "lsw")
+   {
+      MFEM_ABORT("TPV205 only supports --fric-law lsw (slip-weakening "
+                 "via SCEC TPV5 §7-11 closed-form).  Got: '"
+                 << fric_law_cli << "'.  Use --fric-law lsw or omit "
+                 "the flag.");
+   }
    return DispatchedLaw::LSW;
 }
```

(Comment in the existing code says the param is "accepted for symmetry
with TPV102/TPV104 launch scripts but intentionally ignored."  That's
defensible IF the banner accurately reports the dispatch.  Either keep
silent + clearer banner, or reject typos.  Either is an
improvement over silently ignoring.)

---

### [R-408 / R-509] MODERATE — CLAUDE.md self-contradiction (BP5 volume tol unit)

`miniapps/seas/CLAUDE.md`:
- **Line 189** (in the table): `| seas (BP5) | 1e-3 (m/s for velocity) ...`
- **Line 200** (Phase 6 narrative): `--paraview-volume-zfp-tol ... (e.g.,
  velocity in TPV* / **displacement in BP5**)`

The BP5 driver actually registers `displacement` (m), not velocity
(m/s).  Verified at `bp5_verification_full.cpp:1704`: `pv_out->RegisterDomainField("displacement", ...)`.

**Suggested fix:**
```diff
-| seas (BP5)         | 1e-3 (m/s for velocity)      | n/a (no secondary collection) | 1e-12 (slip-rate floor)    |
+| seas (BP5)         | 1e-3 (m, displacement)       | n/a (no secondary collection) | 1e-12 (slip-rate floor)    |
```

---

## LOW findings — abbreviated

### [R-707] LOW — `output_mode_uniformity_checked_` latch is permanent

After first write, mode-setter mutations are not re-validated.  Current
drivers don't mutate mid-run, but contract is implicit.  No fix
recommended unless mid-run mode changes become a feature.

### [R-708] LOW — `_PRINT_ARITH_RE` only matches `python3 -c print(...)`

Confirms the round-2 limitation: the BP5 sbatch's
`$(python3 -c 'print(250*3.156e7)')` form DOES match the regex, so the
estimator can replay it offline.  No fix needed unless other subshell
forms appear in sbatch.

### [R-414..R-417] LOW (round 2, unchanged)

Re-confirmed by reading current code: no ParMesh test for 3-arg
back-compat ctor; `_filter_name_from_args` docstring incomplete;
hard-coded `/Users/...` path in scripts; `--verify-dispatch` log noise
on TPV* drivers.  Defer.

---

## Verification commands (round 3)

```bash
# R-701: confirm the BP5 sbatch uses a subshell tfinal.
grep -E '^\s+--tfinal\s+\$\(' miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch
# Currently: 1 match.  After fix: 0 matches.

# R-701: confirm command-substitution-empty + atof behaviour.
python3 -c 'print(repr(0.0))'  # the value tfinal_override gets when atof("--verify") is called
# Output: '0.0' — confirming the simulation would exit at step 1.

# R-702: cross-check seconds_per_year.
echo "C++ seconds_per_year = $(python3 -c 'print(365.25*24*3600)')"
echo "sbatch value = $(python3 -c 'print(250*3.156e7)')"
echo "C++ value for 250 yr = $(python3 -c 'print(250*365.25*24*3600)')"

# R-703: confirm the comment claims 12 events.
grep -E '12 events|430k uncapped' miniapps/seas/jobs/bp5/bp5_phase6_paraview_zfp_normal_48hr.sbatch
# Expected after fix: empty (or rewritten).

# R-704: confirm estimate_from_sbatch drops per-regime dt.
python3 -c "
import sys; sys.path.insert(0,'miniapps/seas/scripts')
from estimate_from_sbatch import parse_sbatch
from _io_size_schedule import ScheduleConfig
from pathlib import Path
p = parse_sbatch(Path('miniapps/seas/jobs/tpv102/tpv102_phase6_paraview_zfp_dev_2hr.sbatch'))
print('flags:', sorted(p.flags.keys()))
print('paraview-bulk-dt:', p.flags.get('paraview-bulk-dt'))
"
# Expected after fix: the flag is read into ScheduleConfig.

# R-404: confirm no H5Z_ZFP grep in sbatch.
grep -l 'MFEM_USE_H5Z_ZFP' miniapps/seas/jobs/*/*phase6*.sbatch
# Currently: empty.  After fix: 4 lines.

# R-503: confirm --paraview-dt is ambiguous.
python3 miniapps/seas/scripts/estimate_output_size.py \
    --driver tpv102 --inline-mesh --tfinal 1.0 --paraview \
    --paraview-dt 0.5 2>&1 | tail -1
# Expected: error: ambiguous option: --paraview-dt ...

# R-504: confirm flag-name divergence.
diff <(grep -oE 'paraview-(co|nu|inter)seismic-dt' \
       miniapps/seas/drivers/tpv102_driver.cpp | sort -u) \
     <(grep -oE 'paraview-dt-(co|nu|inter)seismic' \
       miniapps/seas/scripts/estimate_output_size.py | sort -u)
# Currently: 3 lines on each side, zero overlap.
```

---

## Recommended PRE-SUBMIT actions (priority order)

1. **R-701 (CRITICAL)**: change the BP5 sbatch `--tfinal` to a hard
   literal (`7889400000`).  This is a 1-line edit that also fixes R-702.
   **Do this before submission.**
2. **R-703 (HIGH)**: update the BP5 sbatch comment block — either
   regenerate from the estimator output or pick a defensible n_events
   value.  Either way the comment must not advertise "12 events" if the
   plan is to size the cap from that assumption.
3. **R-404 (round 2 MODERATE)**: add `grep -q '^MFEM_USE_H5Z_ZFP *= YES'`
   to all four sbatch.  Saves ~30 SU on a bad build.
4. **R-503, R-504, R-409 (round 2 HIGH/MOD)**: fix the Python estimator
   flag naming + defaults.  Required if Phase 7.6 calibration is going
   to happen this week.  Not strictly blocking the Frontera submission.
5. **R-704 (round 3 MOD)**: wire per-regime dt + bulk flags into
   `estimate_from_sbatch.py`.  Required for accurate post-run
   calibration against `REFERENCE_RUNS.md`.
6. **R-706 (round 3 MOD; defensive)**: add a defensive note or
   assertion that the volume HDF dc is never `Save()`-ed when
   `emit_volume_save_ == false`.  No file system effect for the
   current sbatch, but cheap insurance.
7. **R-408/R-509 (round 2 MOD)**: fix CLAUDE.md line 189 to read
   `1e-3 (m, displacement)`.
8. **R-411 (round 2 MOD)**: add a soft-reject on `--fric-law` other than
   `lsw` in TPV205.  Cosmetic for the current sbatch (already
   `--fric-law lsw`).
9. **R-406, R-407 (round 2 MOD)**: complete `COMPRESSION_RATIOS` and
   apply cap in `estimate_n_volume_writes` — both Python-only,
   non-blocking for Frontera but blocking for accurate calibration.
10. **R-705, R-707, R-708 (round 3 LOW)**: defer post-Frontera.

---

## Final assessment

**Recommend HOLD on the BP5 sbatch submission until R-701 is fixed.**
The python3-subshell pattern is a single-line change but the failure
mode is silent and consumes the entire 48-hour reservation.  All other
findings are non-blocking for the four sbatch submission, though
R-703's documentation drift would mislead any future operator reading
the sbatch.

The three TPV* dev-queue sbatch (2-hour wall) can be submitted
**immediately** (R-510, R-411, R-412 are all addressed; the Python
estimator findings R-503/R-504/R-704 only affect the calibration step,
not the runs themselves).  The 2-hour TPV* runs are also the right
proving ground for R-403 (deferred phdf5 module-load risk) — if the
Frontera build dies at the `static_assert` in `paraview_output.hpp:126`
or at `ProbeH5ZZfpPluginOrAbort`, that surfaces fast and cheap.

The C++ side is **correctly implemented**: Phase 6.1 (R-401) static_assert
is template-parameter-dependent and only fires for ParMesh on builds
without parallel HDF5 (verified); Phase 6.2 (R-301) volume-side
compression routes through `SetVolumeHDFCompression` for both the
primary `pv_out` and the secondary `pv_bulk_out`; Phase 6.4 wires the
`--paraview-bulk-*` flags into `pv_bulk_out` (TPV*) and emits the
"no effect" warning on BP5.  Phase 3 snapshot cap is wired into the
schedule with the `SetTotalRunTime(tfinal)` call from both BP5 and TPV*
drivers.

The Python side has the **correct R-402 architecture** (auto-detect from
`$PhysicalNames`) and R-501 is properly addressed by the two-pass
exact-then-prefix match.  The remaining findings (R-503, R-504, R-409,
R-704, R-406, R-407) are all gaps between the Python estimator and the
C++ driver CLI surface — they affect calibration accuracy but not the
Frontera runs.

---

## Summary

- Critical issues: **1** (R-701)
- High issues: **2** (R-703 [new] + R-503, R-504 from round 2 carry)
- Moderate issues: **9** (R-702, R-704, R-706, R-705 [new] + R-404, R-406, R-407, R-409, R-411 carry; R-408/R-509 carry counted once)
- Low issues: **3** (R-707, R-708 [new] + R-414..R-417 carry)
- Plan compliance: **PARTIAL** (R-507 partially-resolved: BP5 sbatch fault-only + 400 ranks; plan §7.5 still text-says-800)
- Verdict: **PASS WITH FIXES — must fix R-701 before BP5 submission;
  the three TPV* dev-queue sbatch can submit now.**

## Unreviewed Areas

- `fem/datacollection.cpp` / `mesh/vtkhdf.cpp`: I read the relevant
  Phase 6.2 dispatch (TSave + EnsureDataset) and verified the ZFP
  cd_values encoding, but did NOT audit the full MFEM patch for
  per-rank chunking edge cases (e.g. rank with zero local DOFs writing
  to a parallel HDF dataset).  Plan §Risk #2 calls out this risk; it's
  covered by `test_paraview_schedule_cap.cpp` and the `np=4` zero-DOF
  unit test (not re-run in this review).
- `miniapps/seas/io/fault_vtu_binary.hpp` (735 lines): the binary-VTU
  back end is not exercised by the four Phase 6 sbatch (all four use
  the HDF5 default).  Skipped.
- `miniapps/seas/scripts/pack_fault_output.py` (630 lines) and its
  test (197 lines): a post-run compaction tool not exercised by the
  four sbatch.  Skipped.
- Unit tests under `miniapps/seas/tests/unit/` (test_paraview_schedule_cap,
  test_volume_hdf_compression, test_vtkhdf_zfp, etc.): the existence of
  each test was confirmed by directory listing; their pass/fail status
  was NOT re-verified in this review (round 2 reported 31 passed, 3
  skipped — unchanged).
