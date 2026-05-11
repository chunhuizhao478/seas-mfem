# Reference Runs for Output-Size Estimator (Phase 7.6)

This file is the ground-truth table that
`test_estimate_output_size.py::test_reference_run_*` reads to compute
`(estimate − measured) / measured`.  Each row's "measured" column is
populated from a real Frontera run; until then the corresponding test
case is `@pytest.mark.skip`.

The acceptance gate (per plan §Phase 7.6 step 4) is:

- Estimate within ±30% of measured: PASS.
- Estimate underestimates by more than 50%: FAIL hard (the user
  would provision insufficient `$SCRATCH`).
- Overestimate up to +100%: tolerated.

## Reference run table

| Run name   | Driver  | Mesh                       | tfinal | Filter (fault) | Filter (volume) | Cap  | Ranks | Measured  | Source / sbatch ID  |
|------------|---------|----------------------------|--------|----------------|-----------------|------|-------|-----------|---------------------|
| ref1_400r  | bp5     | bp5_tandem_exact (1000 m)  | 250yr  | zfp_1e-12      | --no-volume-pv  | 5000 | 400   | (pending) | bp5_phase6_paraview_zfp_normal_48hr.sbatch |
| ref1_800r  | bp5     | bp5_tandem_exact (1000 m)  | 250yr  | deflate-6      | --no-volume-pv  | 5000 | 800   | (pending) | (future row — needs 16-node sbatch) |
| ref2       | tpv102  | tpv102_200m                | 5s     | zfp_1e-12      | zfp_1e-3        | 200  | 400   | (pending) | tpv102_phase6_paraview_zfp_dev_2hr.sbatch |
| ref3       | tpv205  | tpv205_200m                | 5s     | zfp_1e-12      | zfp_1e-3        | 200  | 400   | (pending) | tpv205_phase6_paraview_zfp_dev_2hr.sbatch |

Once a reference run completes on Frontera:

1. Run `du -sb $OUTPUT_DIR` and record the byte count.
2. Replace `(pending)` in the **Measured** column with the formatted
   value (e.g., `12.0 GB`, `3.2 GB`, `0.4 GB`) — must be in a format
   `parse_bytes()` understands.
3. Replace `(pending)` in the **Source** column with the sbatch job
   id and the absolute path on Frontera.
4. Remove the `@pytest.mark.skip(...)` decorator on the matching
   `test_reference_run_*_within_30pct` test.

The test parser (in `test_estimate_output_size.py::_reference_measured_bytes`)
matches the row by `run_name` and reads the **Measured** column.
Keep the table format intact (Markdown pipe-delimited, no extra
columns inserted to the LEFT of "Measured").

## Calibration of compression ratios

Phase 7.4's `COMPRESSION_RATIOS` table contains entries marked
`PLACEHOLDER`.  After each reference run, run:

```
python3 scripts/calibrate_compression.py \
    --hdf $OUTPUT_DIR/fault_surface.vtkhdf \
    --field slip_dip --filter zfp_1e-12
```

(That helper is not yet written; for now the calibration is manual:
run `h5dump -p` on a representative dataset, divide the on-disk
storage size by the raw `n_dofs × 8` byte count, replace the
placeholder entry, and re-run the test suite.)

When every entry in `COMPRESSION_RATIOS` is sourced (i.e., not
`PLACEHOLDER`), Phase 7.4's calibration acceptance gate flips to
PASS and the runtime warning stops firing.
