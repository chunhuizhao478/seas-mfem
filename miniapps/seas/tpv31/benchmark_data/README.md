# TPV31 SCEC reference traces

This directory holds the SCEC TPV31 on-fault station traces used as the
reference for `tpv31/visualize_results.py` (the MFEM-vs-SCEC overlay
tool; plan §R.5 step 1, `PLAN_phase_R_exact_bimaterial_riemann_rev3.md`).
Two community submissions are committed: `scec_eqdyna/` (EQdyna, Benchun
Duan) and `scec_seisol/` (SeisSol ADER-DG, Thomas Ulrich), 30 on-fault
stations each.

## Source

The committed traces are the on-fault station submissions from the
SCEC Code Verification Working Group submission area for problem TPV31:

  https://strike.scec.org/cvws/cgi-bin/cvws.cgi?problem=31

The TPV31 problem definition (see `tpv31/benchmark_document/TPV31_32_Description_v03.pdf`,
spec p. 12) requires submissions to provide:

  * **30 on-fault stations** (≥ 4 of which are required for the
    automated comparison)
  * Off-fault stations (depth and surface)
  * Rupture-time contour file

The full 30-station on-fault set from each of two community-verified
submissions is committed here — SeisSol (`scec_seisol/`) and EQdyna
(`scec_eqdyna/`).  Stations span the 3 along-strike × 10 down-dip grid
(strike {0, 6, 12} km × depth {0, 0.2, 0.5, 1.0, 2.4, 3.0, 5.0, 7.5,
10.0, 12.0} km); e.g. the SeisSol trace 0 m along strike / 7500 m
down-dip is `scec_seisol/tpv31_seisol_x2_0_x3_7.5.txt`.

## File layout

The committed reference set is **60 `.txt` files** — 30 on-fault
stations per code:

```
scec_seisol/   30 files   (SeisSol ADER-DG, Thomas Ulrich)
scec_eqdyna/   30 files   (EQdyna, Benchun Duan)
```

Each filename follows the glob that `tpv31/visualize_results.py` uses
(`visualize_results.py:201,575`, helper `reference_filename`):

```
tpv31_{code}_x2_<strike_km>_x3_<depth_km>.txt
```

where `{code}` is `seisol` or `eqdyna`, and the `x2_*`/`x3_*` suffixes
are the along-strike and down-dip station coordinates in km with
trailing zeros stripped (e.g. `tpv31_seisol_x2_0_x3_7.5.txt`,
`tpv31_eqdyna_x2_12_x3_10.txt`).  The 3 along-strike × 10 down-dip grid
is built in `_build_stations()` (`visualize_results.py:71`).

### File format (the parser is the authority)

The format authority is the script's reader
`load_reference_file` (`visualize_results.py:158`), backed by
`_parse_numeric_table` (`:100`) — NOT this README.  In summary:

* A **multi-line `#`-prefixed comment header** (possibly indented, as in
  the EQdyna files) carries provenance and column documentation
  (`# code=...`, `# Column #1 = Time (s)`, etc.).  `_parse_numeric_table`
  skips every line that is blank or starts with `#`.
* A **field-name line** `t h-slip h-slip-rate h-shear-stress v-slip
  v-slip-rate v-shear-stress n-stress` follows the header.  Because it
  begins with the letter `t`, the `float(...)` parse raises `ValueError`
  and the line is skipped automatically (`:116-119`).
* **8 numeric columns** per data row (`min_cols=8`):

  | Col | Quantity                      | Units | Notes |
  |-----|-------------------------------|-------|-------|
  | 1   | time                          | s     | fixed cadence (SeisSol ≈ 4.875e-3 s; EQdyna 5e-3 s) |
  | 2   | horizontal (h) slip           | m     | h = strike |
  | 3   | horizontal (h) slip rate      | m/s   | h = strike |
  | 4   | horizontal (h) shear stress   | MPa   | h = strike |
  | 5   | vertical (v) slip             | m     | v = down-dip |
  | 6   | vertical (v) slip rate        | m/s   | v = down-dip |
  | 7   | vertical (v) shear stress     | MPa   | v = down-dip |
  | 8   | normal stress                 | MPa   | **compression-NEGATIVE in the file** (SCEC extension-positive) |

* Stresses are already in **MPa** (no Pa→MPa conversion, unlike the
  MFEM `.dat`).  The reader **negates column 8** so σ_n overlays the
  MFEM compression-positive convention (`:188`).
* The reference does **not** provide μ_eff; the reader fills that channel
  with NaN (`:191`).  The quantitative gate (`--tol-rms` / `--tol-peak`)
  therefore EXCLUDES `mu_eff` (a NaN comparison would silently pass).

The MFEM side is a 9-column `.dat` (`<prefix>_station_<name>.dat`,
stresses in **Pa**, σ_n compression-positive, plus a μ_eff column);
see `load_mfem_file` (`visualize_results.py:125`).

## Acquisition status

The reference set **is committed** to this repository: 60 `.txt` files
(`scec_seisol/` 30 + `scec_eqdyna/` 30), one per on-fault station per
code, covering the full 3 × 10 grid.  No download step is required to
run the overlay or the quantitative gate.

If a future re-acquisition is ever needed and fails (login wall,
network unavailable, or the SCEC server offline), the plan §R.5 Edge
Case says **STOP and ask the user for an alternate data source** rather
than proceed without reference data.

## Acceptance tolerance (the quantitative gate)

`tpv31/visualize_results.py` exposes a pass/fail gate via the
`--tol-peak <x>` and `--tol-rms <y>` flags.  Supplying EITHER flag
switches the script from plotting to the gate (it then runs headless —
no matplotlib / display needed) and `sys.exit`s non-zero on any
exceedance.  The gate compares each MFEM station trace (`--mfem DIR`,
prefix auto-detected) against the SeisSol reference (code `seisol`):

* **Overlap window + interpolation.**  MFEM (adaptive RK45) and the
  SeisSol `.txt` (fixed cadence) are on different time grids.  For each
  channel the gate restricts to `[t0, t1] = [max(mfem_t[0], ref_t[0]),
  min(mfem_t[-1], ref_t[-1])]` and interpolates the MFEM channel onto
  the **reference samples inside `[t0, t1]`** (`np.interp`) before
  differencing — never onto the full reference grid (which would
  flat-line a truncated run's tail).
* **Minimum-coverage guard.**  If the overlap spans less than
  `MIN_COVERAGE_FRAC` (= 0.90) of the reference span, the station FAILS
  with an explicit "insufficient coverage" message — a wall-truncated
  run cannot masquerade as agreement.
* **Peak-normalized, zero-safe metric** (per gated channel):
  `peak_rel = max|mfem-ref| / max(max|ref|, floor)` and
  `rms_rel  = sqrt(mean (mfem-ref)^2) / max(sqrt(mean ref^2), floor)`,
  with a small denominator `floor` so pre-nucleation zero-crossings do
  not blow the relative error up.
* **Gated channels (by name):** `V_strike`, `slip_strike`,
  `tau_strike`, `sigma_n`; plus `V_dip`, `slip_dip`, `tau_dip` only for
  stations with non-trivial dip motion.  `mu_eff` is always EXCLUDED
  (its reference is all-NaN).

The pass-band is whatever the user supplies on the flags; ~5–10% peak
relative error is appropriate because SeisSol / EQdyna are different
codes from this implementation.  Gate-metric tests live in
`tpv31/test_visualize_results_tol.py` (NEW-7.1 / 7.4 / 7.5).
