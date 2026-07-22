# Phase-0 implementation review + fix record (2026-07-22)

Two-lens adversarial review (ops-correctness + methodology/consistency) of the Phase-0
measurement harness. Ops verdict: **"will NOT run cleanly as written"** (1 CRITICAL). Method
verdict: **sound arithmetic, 1 CRITICAL formula defect + 4 doc issues.** All 11 findings applied.

## Ops lens

| ID | Sev | Finding | Fix applied |
|---|---|---|---|
| P0-1 | **CRIT** | `perf_wrap.sh` used `exec perf record` with a fallback only for perf *absence*; if perf exists but `perf_event_paranoid≥2` (Expanse norm) it exits before exec'ing the driver → rank 0 never enters MPI_Init → all 256 ranks hang for the full 8 h walltime | Runtime probe `perf record -o /dev/null -- true` before wrapping; always fall through to `exec "$@"`; `--kill-on-bad-exit=1` on the leg-2 srun as backstop |
| P0-2 | MOD | Leg-3 `srun -n 128 --ntasks-per-node=64` with default block distribution packs the 64 ranks onto one socket → still bandwidth-contended → the "2× per-core bandwidth" premise false → memory-vs-kernel attribution corrupted | `--distribution=cyclic:cyclic --cpu-bind=cores` spreads ranks across both sockets/8 NUMA domains |
| P0-3 | MOD | Uncapped `perf record --call-graph dwarf` over a ~3 h run → tens–hundreds of GB perf.data + I/O perturbs the timing | `-F 99 --call-graph fp` (bounded; self-time split needs no frame pointers); noted OUT_ROOT is on lustre |
| P0-4 | LOW | bench compile relies on `MFEM_LINK_FLAGS` carrying `-I` paths (version-dependent) | explicit `-I<mfem-root>` appended in `run_bench.sh` |
| P0-5 | LOW | window step-count one-ended + %100 print granularity | grep BOTH endpoints per window; `Δsteps` differenced; ±100 noted |
| P0-6 | LOW | `profile.mpi` + perf both on rank 0 bias it into the MPI straggler | `profile.mpi` dropped from leg-2 Caliper (kept in leg-5, where comm is the point) |

**Ops checked-OK:** window regex matches the driver's `step N/M  t = …` print; `CALI_CONFIG`
exported (comma-safe) not via `--export`; SLURM permits leg-3 under-subscription; the non-MPI
bench runs fine under `srun --mpi=pmi2` (independent copies, PMI unused); `run_bench` link-as-compile
works on standard config.mk with `-O3` overriding `-O2`; `SLURM_PROCID==0` selects exactly rank 0.

## Methodology lens

| ID | Sev | Finding | Fix applied |
|---|---|---|---|
| F1 | **CRIT** | Gate-reset formula hardcoded `cores(256)`, but Leg 3 runs 128 ranks → applying it verbatim doubles Leg-3's per-update and self-corrupts the memory-starvation conclusion driving the Phase-2 case | Formula parameterized `ranks_in_window` (256 Leg2/5, 128 Leg3); Leg-3 row annotated; anchor check `6465×256/(2.46e6×2732)=245.8` added |
| F2 | MOD | Rupture window [1.0,1.5] *ends at* the V_max peak (~1.48 s) → most of it is the low-V pre-peak ramp → understates the fault share the defer decision hinges on | Window → **[1.4,1.8]** (brackets+trails the peak); sbatch `--tfinal 1.5→1.8` (3.2 h, fits 8 h); memo asks for the peak-step share too |
| F3 | MOD | Two dangling references to the removed per-QP subregions ("Caliper subregions on", "attributed by subregion") contradicted the rev-3.1 perf methodology | Both reworded to perf self-time / function-level Caliper (plan lines 159/241/630) |
| F4 | MOD | perf methodology didn't note the -O2 inline-attribution limit; some memo rows perf can't populate | Noted: out-of-line hotspots (LU, CalcShape, GodunovFlux::Interior, ComputeVolumeRHS, ApplyMassInverse) attribute by self-time; inlined blocks fold to parent → cross-check the local `sample` |
| F5 | LOW | Phase-2 go/no-go `B≥2×` applied to *partial*-fusion bench B (a lower bound on full-fusion) → could kill a viable phase | Bar restated: `R<1.6` clean no-go; `[1.6,2.0)` requires a full-fusion bench variant before deciding; B labeled a lower bound |

**Method checked-OK:** 245.9 formula reproduces (245.8); stepping-wall consistency; cfl 0.5 from
the deck; t=1.8 reachable in 8 h; `build_expanse.sh` lib64 fix correctly guarded, idempotent, no
MUMPS/non-PETSc/local regression; all 5 deliverables exist and match the plan's status.

## Net

Harness re-verified after fixes: `bash -n` clean; `run_bench.sh` compiles+runs locally
(A 5.26 / B 1.23× / C 0.45× / D 1.31×, bit-identical checksums); `--face-cache` composes.
**Phase-0 is implementation-complete and locally validated; the Expanse job awaits user approval.**
