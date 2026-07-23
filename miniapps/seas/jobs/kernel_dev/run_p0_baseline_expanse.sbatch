#!/bin/bash
# =============================================================================
# KERNEL PLAN — PHASE 0 baseline (B0) measurement job.  SDSC Expanse.
#   document/kernel_dev/PLAN_ader_kernel_efficiency_2026-07-21.md (rev 3), Phase 0.
#
# Produces the B0 gate-reset inputs — the numbers every later phase gate is
# derived from — in ONE job, on TPV104-200m, 2 nodes.  NOTHING here changes
# numerics: it is measurement + a standalone kernel benchmark.
#
# LEGS (toggle with P0_LEG_*):
#   1  BLAS check            grep the deployed config.mk for MFEM_USE_LAPACK.
#   2  B0 GTS, two windows   lts=off, levers+--face-cache, EndTime 1.8 s; the
#                            step log yields per-update us-core over [0,0.5]
#                            (pre-rupture) AND [1.4,1.8] (rupture — brackets AND
#                            trails the V_max peak ~1.48s, so the fault share is
#                            not diluted by the pre-peak low-V ramp).  Caliper
#                            region report (function level) + `perf record` on
#                            rank 0 (SUB-function attribution: LU / CalcShape /
#                            rotations — the perf-sampling methodology REPLACES
#                            the plan's per-QP Caliper subregions, which would
#                            have perturbed the very timing they measured).
#   3  layout A/B            one short GTS leg at 64 ranks/node (half-populated)
#                            -> per-update us-core at 128 vs 64 ranks/node,
#                            separating per-core memory starvation from kernel
#                            structure in the 245.9 us figure.
#   4  kernel bench spike    compile tests/bench/bench_ader_kernel_variants.cpp
#                            on Rome; run single-core AND full-occupancy (one
#                            copy per core) -> contended A-vs-B ratio = the
#                            Phase-2 go/no-go (>= 2x contended B-vs-A).
#   5  comm skew capture     one short LTS leg (lts=rate2, fault interleave) with
#                            Caliper profile.mpi -> MPI_Waitall Max/Avg (the
#                            wire-vs-skew first-order signal for the comm plan).
#
# APPROVAL: submitting this needs the user's go-ahead (standing rule).  The
# heavy pole is leg 2 (GTS to 1.8 s ~ 3.2 h stepping at cfl=0.5; fits the 8 h
# walltime with the short legs); the rest are minutes.  OUT_ROOT is under
# SEAS_MFEM_ROOT which lives on /expanse/lustre — perf.data + snapshots land there.
# =============================================================================
#SBATCH --job-name=kernel-p0-b0
#SBATCH --account=lbl107
#SBATCH --partition=compute
#SBATCH --constraint="lustre"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=140000M
#SBATCH -t 08:00:00
#SBATCH -o kernel-p0-b0.%j.out
#SBATCH -e kernel-p0-b0.%j.out
#SBATCH --export=ALL

set -o pipefail
export LC_ALL=C LANG=C

module purge
module load slurm 2>/dev/null || true
module load cpu/0.17.3b
module load gcc/10.2.0
for _mpi in openmpi/4.1.5 openmpi/4.1.3 openmpi/4.1.1 openmpi; do
    module load "${_mpi}" 2>/dev/null && module is-loaded "${_mpi}" 2>/dev/null && { echo "OpenMPI module: ${_mpi}"; break; }
done

# --- resolve the seas-mfem root ---------------------------------------------
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SUBMIT_DIR}}"
_walk="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5 6 7; do
    [[ -f "${_walk}/miniapps/seas/Makefile" ]] && { SEAS_MFEM_ROOT="$(cd "${_walk}" && pwd)"; break; }
    _walk="${_walk}/.."
done
[[ -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]] || { echo "ERROR: cannot locate seas-mfem root from ${SUBMIT_DIR}"; exit 1; }

export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_USE_FILE_LOCKING=FALSE

cd "${SEAS_MFEM_ROOT}/miniapps/seas"
DRIVER=./seas_spatial_dyn_driver
[[ -x "${DRIVER}" ]] || { echo "ERROR: build ${DRIVER} first (bash ${SEAS_MFEM_ROOT}/build_expanse.sh)"; exit 1; }

DECK=jobs/lts_phase5/tpv104_200m_lts_speed_expanse
CFG_OFF="${DECK}/tpv104_200m_lts_off.toml"
CFG_RATE2="${DECK}/tpv104_200m_lts_rate2.toml"
MESH_LN="${DECK}/mesh_tpv104_200m.msh"
[[ -e "${MESH_LN}" ]] || { echo "ERROR: 200m mesh not staged at ${MESH_LN} (gitignored; scp it)"; exit 1; }

OUT_ROOT="${SEAS_MFEM_ROOT}/miniapps/seas/${DECK}/runs/p0_${SLURM_JOB_ID:-local}"
mkdir -p "${OUT_ROOT}"
COMMON="--mixed-flux none --time-integrator ader --ader-order 4 --cfl-dg-safety 1.0 --deriv-cache --shared-ck-recursion --face-cache --print-derived"
P0_LEG_1="${P0_LEG_1:-1}"; P0_LEG_2="${P0_LEG_2:-1}"; P0_LEG_3="${P0_LEG_3:-1}"; P0_LEG_4="${P0_LEG_4:-1}"; P0_LEG_5="${P0_LEG_5:-1}"

echo "############ KERNEL PHASE-0 BASELINE | job ${SLURM_JOB_ID}"
echo "=== git $(git rev-parse --short HEAD 2>/dev/null) | nodes=${SLURM_NNODES} ranks=${SLURM_NTASKS} | out=${OUT_ROOT}"

# ---- LEG 1: BLAS check ------------------------------------------------------
if [[ "${P0_LEG_1}" == "1" ]]; then
  echo ""; echo "===== LEG 1: MFEM_USE_LAPACK in the deployed build ====="
  grep -E "^MFEM_USE_LAPACK" "${SEAS_MFEM_ROOT}/config/config.mk" || echo "  (MFEM_USE_LAPACK not set -> defaults NO; the slow internal LU path)"
  ldd "${DRIVER}" 2>/dev/null | grep -iE "openblas|lapack|libblas" \
    || echo "  (nothing in ldd — expected when openblas is STATICALLY linked; check nm below)"
  nm "${DRIVER}" 2>/dev/null | grep -m2 -E " T (dgemm_|dgetrs_)" \
    || echo "  (no dgemm_ symbol — if config.mk also says NO, the LU tax runs MFEM internal LSolve/USolve)"
  echo "  Interpretation: config.mk YES + a defined dgemm_ (T) = optimized LAPACK LU path."
fi

# ---- LEG 2: B0 GTS, two windows, perf sub-function attribution ---------------
if [[ "${P0_LEG_2}" == "1" ]]; then
  _T0=$SECONDS
  echo ""; echo "===== LEG 2: B0 GTS (EndTime 1.8s; windows [0,0.5] pre-rupture, [1.4,1.8] rupture) ====="
  OUT_B0="${OUT_ROOT}/b0_gts"; mkdir -p "${OUT_B0}"; LOG_B0="${OUT_ROOT}/b0_gts.log"
  # profile.mpi DROPPED here (review P0-6): perf covers sub-function attribution;
  # profile.mpi + perf both on rank 0 would bias rank 0 into the MPI straggler.
  export CALI_CONFIG="runtime-report(calc.inclusive,region.count,max_column_width=60,output=${OUT_ROOT}/b0_gts.cali-region-report.txt)"
  # perf on rank 0 ONLY, and ONLY if perf can actually open a perf_event (a
  # runtime probe — 'perf record -- true' — catches paranoid>=2 / missing
  # kernel-tools).  ANY failure falls through to 'exec "$@"' so the MPI job never
  # hangs (review P0-1 CRITICAL: bare 'exec perf' would kill rank 0 before
  # MPI_Init and hang all 256 ranks to walltime).  -F 99 --call-graph fp is
  # cheap + bounded (P0-3): self-time symbol split works without frame pointers;
  # the call tree is best-effort.
  cat > "${OUT_ROOT}/perf_wrap.sh" <<'WRAP'
#!/bin/bash
if [[ "${SLURM_PROCID:-1}" == "0" ]] && command -v perf >/dev/null 2>&1 \
   && perf record -F 99 -o /dev/null -- true >/dev/null 2>&1; then
    exec perf record -F 99 -g --call-graph fp -o "${P0_PERF_OUT}" -- "$@"
fi
exec "$@"
WRAP
  chmod +x "${OUT_ROOT}/perf_wrap.sh"
  export P0_PERF_OUT="${OUT_ROOT}/b0_gts.perf.data"
  # --kill-on-bad-exit=1: if rank 0 dies anyway, tear the step down (do not hang).
  srun --mpi=pmi2 -n "${SLURM_NTASKS}" --kill-on-bad-exit=1 --unbuffered stdbuf -oL -eL \
       "${OUT_ROOT}/perf_wrap.sh" "${DRIVER}" --config "${CFG_OFF}" ${COMMON} \
       --tfinal 1.8 --paraview --paraview-fault-vtu --paraview-max-snapshots 40 \
       --output-dir "${OUT_B0}" > "${LOG_B0}" 2>&1
  echo "B0 rc=$?  (steps + t + V_max in ${LOG_B0})"
  echo "--- window step boundaries (per-update us-core = wall*ranks(256)*1e6/(elems*Δsteps)) ---"
  echo "  pre-rupture [0,0.5]  start:"; grep -E "^step " "${LOG_B0}" | head -1
  echo "                        end:";  grep -E "^step .* t = 0\.[45]" "${LOG_B0}" | tail -1
  echo "  rupture [1.4,1.8] (brackets+trails the V_max peak ~1.48s) start:"; grep -E "^step .* t = 1\.[34]" "${LOG_B0}" | head -1
  echo "                                                             end:"; grep -E "^step .* t = 1\.[78]" "${LOG_B0}" | tail -1
  echo "--- Caliper region head ---"; head -20 "${OUT_ROOT}/b0_gts.cali-region-report.txt" 2>/dev/null
  if [[ -f "${P0_PERF_OUT}" ]]; then
      echo "--- perf sub-function self-time (top 25; the LU/CalcShape/rotation split) ---"
      perf report -i "${P0_PERF_OUT}" --stdio --sort=symbol -g none 2>/dev/null | grep -E "^\s+[0-9]" | head -25 \
        || echo "  (perf report unavailable on this node; copy ${P0_PERF_OUT} to a perf-enabled host)"
  else
      echo "  NOTE: perf record produced no data (perf absent or paranoid>1); Caliper function-level only."
  fi
  echo "LEG 2 WALL = $((SECONDS-_T0)) s   <-- per-update = WALL*256*1e6/(2464689*Δsteps)"
  unset CALI_CONFIG
fi

# ---- LEG 3: node-layout A/B (128 vs 64 ranks/node) --------------------------
if [[ "${P0_LEG_3}" == "1" ]]; then
  _T0=$SECONDS
  echo ""; echo "===== LEG 3: layout A/B — 64 ranks/node (half-populated), short GTS ====="
  OUT_L="${OUT_ROOT}/layout64"; mkdir -p "${OUT_L}"; LOG_L="${OUT_ROOT}/layout64.log"
  # 2 nodes x 64 = 128 ranks.  --distribution=cyclic:cyclic + --cpu-bind=cores
  # SPREADS the 64 ranks across BOTH sockets / all 8 NUMA domains (review P0-2:
  # default block packing would cram them onto one socket, so they'd still
  # contend on one memory controller and the "~2x bandwidth" premise would be
  # false, corrupting the memory-vs-kernel attribution this leg exists to make).
  # NB: its per-update uses ranks=128, NOT 256 (memo formula parameter).
  srun --mpi=pmi2 -n 128 --ntasks-per-node=64 --distribution=cyclic:cyclic --cpu-bind=cores \
       --unbuffered stdbuf -oL -eL \
       "${DRIVER}" --config "${CFG_OFF}" ${COMMON} --tfinal 0.5 \
       --output-dir "${OUT_L}" > "${LOG_L}" 2>&1
  echo "layout64 rc=$?  (per-update us-core = wall*ranks(128)*1e6/(elems*Δsteps);"
  echo "  compare to LEG 2's [0,0.5] window: a large speedup at 64/node => the"
  echo "  245.9 us is per-core memory starvation, not kernel structure)"
  grep -E "^step " "${LOG_L}" | head -1; grep -E "^step .* t = 0\.[45]" "${LOG_L}" | tail -1
  echo "LEG 3 WALL = $((SECONDS-_T0)) s   <-- per-update = WALL*128*1e6/(2464689*Δsteps)  [ranks=128!]"
fi

# ---- LEG 4: kernel bench spike (contended A-vs-B = Phase-2 go/no-go) ---------
if [[ "${P0_LEG_4}" == "1" ]]; then
  echo ""; echo "===== LEG 4: kernel microbenchmark on Rome (single-core + full-occupancy) ====="
  bash "${SEAS_MFEM_ROOT}/miniapps/seas/tests/bench/run_bench.sh" "${OUT_ROOT}" > "${OUT_ROOT}/bench.log" 2>&1
  BENCH_RC=$?
  tail -60 "${OUT_ROOT}/bench.log"
  [[ ${BENCH_RC} -eq 0 ]] || echo "  (bench rc=${BENCH_RC} — see ${OUT_ROOT}/bench.log; non-fatal to the rest of the job)"
fi

# ---- LEG 5: comm Waitall skew capture (LTS leg, short) ----------------------
if [[ "${P0_LEG_5}" == "1" ]]; then
  _T0=$SECONDS
  echo ""; echo "===== LEG 5: comm skew — short LTS leg, Caliper profile.mpi (Waitall Max/Avg) ====="
  OUT_C="${OUT_ROOT}/comm_lts"; mkdir -p "${OUT_C}"; LOG_C="${OUT_ROOT}/comm_lts.log"
  export CALI_CONFIG="runtime-report(profile.mpi,calc.inclusive,region.count,max_column_width=60,output=${OUT_ROOT}/comm_lts.cali-region-report.txt)"
  SEAS_LTS_FAULT_INTERLEAVE=1 srun --mpi=pmi2 -n "${SLURM_NTASKS}" --unbuffered stdbuf -oL -eL \
       "${DRIVER}" --config "${CFG_RATE2}" ${COMMON} --tfinal 0.5 \
       --output-dir "${OUT_C}" > "${LOG_C}" 2>&1
  echo "comm_lts rc=$?"
  echo "--- MPI_Waitall Max vs Avg (Max >> Avg = late-sender SKEW, not wire time) ---"
  grep -iE "MPI_Waitall|MPI_Isend|MPI_Irecv" "${OUT_ROOT}/comm_lts.cali-region-report.txt" 2>/dev/null
  echo "  If Max/Avg >> 1 and skew dominates, the comm plan's payoff degrades and the"
  echo "  fault-weighted partition (METIS ncon=2 v3) becomes co-requisite."
  echo "LEG 5 WALL = $((SECONDS-_T0)) s"
  unset CALI_CONFIG
fi

echo ""; echo "############ PHASE-0 done.  Fill the gate-reset memo from these legs:"
echo "  document/kernel_dev/B0_gate_reset_2026-07-22.md"
exit 0
