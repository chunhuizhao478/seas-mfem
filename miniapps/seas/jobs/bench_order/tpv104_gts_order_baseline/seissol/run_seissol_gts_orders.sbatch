#!/bin/bash
# =============================================================================
# TPV104 200 m  ---  GTS ORDER SPEED BASELINE, SeisSol SIDE (o2, o3, o4)
#   *** EXPANSE (SDSC).  Partner job: ../run_mfem_gts_orders.sbatch ***
#
# WHAT THIS IS ----------------------------------------------------------------
#   Three SeisSol legs, o2 / o3 / o4, on the SAME 200 m mesh as the MFEM legs
#   (tpv104.puml.h5 = pumgen's conversion of a byte-for-byte copy of
#   miniapps/seas/tpv104/mesh/tpv104_200m.msh; 2,464,689 tets, 37,904 DR faces),
#   GLOBAL time stepping (ClusteredLTS = 1), CFL 0.5, EndTime 2.0 s, ALL OUTPUT
#   OFF.  R2 ACCURACY-MATCHED PAIRING: SeisSol o = MFEM p + 1, so
#       o2 <-> MFEM p1 ,  o3 <-> MFEM p2 ,  o4 <-> MFEM p3
#   (both codes then advance 4 / 10 / 20 modes per element x 9 fields, with the
#   same dt and the same step count).  This is a COST baseline at NOMINALLY
#   matched order — NO accuracy claim is made.
#
# WHY A SEPARATE JOB FROM THE MFEM LEGS ---------------------------------------
#   The rank shape is fixed at allocation time and the two codes' validated
#   shapes are incompatible: MFEM = 2 x 128 ranks x 1 core (pure MPI); SeisSol =
#   2 x 8 ranks x 16 cores (15 OMP threads + 1 reserved core).  This job's
#   prologue DERIVES OMP_NUM_THREADS and SRUN_CPUS_PER_TASK from
#   $SLURM_CPUS_PER_TASK (SLURM 23.02 does not propagate -c to srun), so it
#   cannot run under MFEM's -c 1 header.  R5 forbids changing either shape.
#   Also: different PMI (pmix here, pmi2 there), different filesystem root (this
#   deck + the 103 MiB .puml.h5 live on /expanse/lustre, not in the repo),
#   different account, and a ~1 h vs ~5 h walltime.
#
# THREE DECKS WITH BYTE-IDENTICAL BODIES, THREE BINARIES.  SeisSol's ORDER IS
#   COMPILE-TIME — no key in any .par sets it — so the physics config is the same
#   for all three legs and ONLY the executable differs:
#     parameters_gts_o2.par  ->  $SEISSOL_ROOT/SeisSol/build-rome-o2-elastic/seissol-elastic-o2-f64
#     parameters_gts_o3.par  ->  $SEISSOL_ROOT/SeisSol/build-rome-o3-elastic/seissol-elastic-o3-f64
#     parameters_gts_o4.par  ->  $SEISSOL_ROOT/SeisSol/build-rome-o4-elastic/seissol-elastic-o4-f64  (pre-existing)
#   The three .par files differ ONLY in their leading '!' comment header (each
#   states its own binary, its own pre-registered dt and step count).  That is
#   deliberate PROVENANCE: a leg's log, deck and binary can be tied together
#   after the fact.  The identical-body contract is machine-checkable and this
#   job checks it before running anything:
#     diff <(grep -v '^!' parameters_gts_o2.par) <(grep -v '^!' parameters_gts_o4.par)
#   Every leg's deck is independently run through the full 13-assertion contract
#   pre-flight below, so a drifted body cannot slip through on one order only.
#
#   o2/o3 must be built first with the CO-LOCATED build job:
#     jobs/bench_order/tpv104_gts_order_baseline/seissol/build_seissol_o2_o3_expanse.sbatch
#   (that copy — and ONLY that copy — carries the CMakeCache ORDER-only fairness
#   diff, the ORDER:STRING assertion, the stale-binary mtime check and the
#   truncated-build-dir cleanup.  Do not use any other script of that name.)
#   This job REFUSES (skips, and exits non-zero) any leg whose binary is missing,
#   and asserts all three binaries carry the SAME source commit — otherwise the
#   study measures source revisions, not orders.
#
# DECK CONTRACT (the .par must already satisfy this; asserted below) ----------
#   ClusteredLTS = 1     GTS.  (The old speed deck had 2 = rate-2 LTS.)
#   CFL          = 0.5   matched to the MFEM decks' [numerics].cfl
#   EndTime      = 2.0 , t_0 = 1.0   matched window + nucleation ramp
#   OutputPointType = 0  kills BOTH &Elementwise and &Pickpoint fault output
#   SurfaceOutput   = 0  kills free-surface output AND the PGV fold, which runs
#                        INSIDE the compute stopwatch at EVERY local timestep
#   ReceiverOutput  = 0  off-fault receivers (the key DEFAULTS TO TRUE)
#   EnergyOutput = 0 , EnergyTerminalOutput = 0
#                        NOT merely an IO knob: energyoutput=true switches on
#                        extra per-GP work inside the friction kernel
#                        (DRParameters.cpp:154 -> BaseFrictionLaw.h:136,173).
#                        The resulting compute figure is therefore NOT comparable
#                        to the historical 111 s/sim-s, which had energy ON.
#   Checkpoint   = 0     explicit (the key defaults to TRUE in v1.3.x)
#   Format       = 10    no wavefield; iOutputMask must stay (readOrFail).
#   With all periodic output off there is EXACTLY ONE simulation epoch, so no
#   sync point clips dt and steps = ceil(EndTime/dt) exactly (the old 0.05 s
#   cadence added ~40 truncated steps, +0.7%, purely from IO).
#
# NOTE ON THE BINARY VERSION: the qstore build is v1.3.1-2135-gdc6db6513, NOT
#   v1.1.3.  The v1.1.3-era `numfluxnearfault` FATAL pre-flight from the Phase-5
#   sbatch is therefore RETIRED here; `surfacevtkorder` is honoured on this
#   binary (which is exactly why SurfaceOutput must be 0 for a clean wall).
#
# DEPLOY (this folder is NOT run from the repo tree on Expanse) ----------------
#   This folder is SELF-CONTAINED except for the 103 MiB mesh: the three .par
#   files AND tpv104_{fault,material}.yaml + tpv104_{receivers,faultreceivers}.dat
#   all live here, so ONE rsync stages a complete run directory.  (The two .dat
#   files are inert under ReceiverOutput = 0 / OutputPointType = 0 but the
#   pre-flight still requires them, and the fault yaml carries the FL=103
#   spatial rs_a / rs_srW that the physics actually needs.)
#   Use a FRESH directory — never reuse an existing SeisSol run dir.
#
#   rsync -avP <repo>/miniapps/seas/jobs/bench_order/tpv104_gts_order_baseline/seissol/ \
#         czhao1@login.expanse.sdsc.edu:/expanse/lustre/projects/lbl107/czhao1/safs/tpv104_gts_order/
#   rsync -avP ~/Downloads/tpv104/tpv104.puml.h5 \
#         czhao1@login.expanse.sdsc.edu:/expanse/lustre/projects/lbl107/czhao1/safs/tpv104_gts_order/
#   cd /expanse/lustre/projects/lbl107/czhao1/safs/tpv104_gts_order
#   sbatch run_seissol_gts_orders.sbatch
# WATCH:  squeue -u $USER ; tail -f tpv104_gts_order_ss.<jobid>.out
#
# Each leg runs in its OWN subdirectory leg_o<N>.<jobid>/ with the inputs
# symlinked in, so no two legs can collide on outputs/ and no stop-file from an
# earlier run can silently short-circuit a leg (the Phase-5 sbatch's
# outputs/.SIMULATION_DONE marker is NOT used here — a repeat run for run-to-run
# spread must actually re-run).
#
# ENV OVERRIDES ---------------------------------------------------------------
#   BENCH_SS_LEGS  (default "2 3 4")  orders to run, cheapest first.
#   BENCH_SS_PAR   (default: per leg, parameters_gts_o<N>.par)  Setting this
#                  forces ONE deck onto ALL legs — only for a deliberate control.
#   BENCH_MD5      (default 1)  0 => skip the mesh md5 pin.
#   SEISSOL_ROOT   (default the qstore build root)
#
# COST: ~25 min (o4) + ~10 min (o3) + ~4 min (o2) of stepping plus PUML read and
#   partition per leg; 3 h is requested as headroom.  ~150 SU total.
# =============================================================================
#SBATCH --job-name=tpv104_gts_order_ss
#SBATCH --account=usc143              # the account the SeisSol qstore build ran under
                                      # (~291k SU free 2026-07-24); override: sbatch -A lbl107
#SBATCH --partition=compute
#SBATCH --constraint="lustre"         # REQUIRED: this run lives on /expanse/lustre
#SBATCH --nodes=2                     # R5: SAME node count as the MFEM job (2 nodes)
#SBATCH --ntasks-per-node=8           # 8 ranks x 16 cores = 128 cores/node = THE WHOLE NODE.
                                      # EPYC 7742 = 8 NUMA x 16 cores -> 1 rank per NUMA domain.
#SBATCH --cpus-per-task=16            # SeisSol's VALIDATED shape.  DO NOT change (R5).
#SBATCH --mem=140000M
#SBATCH -t 03:00:00
#SBATCH -o tpv104_gts_order_ss.%j.out
#SBATCH -e tpv104_gts_order_ss.%j.out
#SBATCH --export=ALL

set -o pipefail
export LC_ALL=C
export LANG=C

# never dump core: SeisSol's watchdog abort() writes ~2.5 GB/rank otherwise.
ulimit -c 0

SEISSOL_ROOT=${SEISSOL_ROOT:-/expanse/projects/qstore/usc143/qwxdev/apps/expanse/rocky8.8/seisol_build}
RUN_DIR=${SLURM_SUBMIT_DIR:-$PWD}
cd "$RUN_DIR" || exit 1

# Deck selection is PER LEG: order is compile-time, so each order gets its own
# .par whose '!' header names its binary, its pre-registered dt and its step
# count.  BENCH_SS_PAR forces one deck onto every leg (control runs only).
par_for () { if [ -n "$BENCH_SS_PAR" ]; then echo "$BENCH_SS_PAR"; else echo "parameters_gts_o${1}.par"; fi; }
LEGS=${BENCH_SS_LEGS:-2 3 4}
NELEM=2464689                                  # md5-pinned; cross-checked against the log
NDRFACE=37904
MESH=tpv104.puml.h5
MESH_MD5_PIN=5911b0b4ec29a8f61968d5b32ed0a11c
CORES_ALLOC=$(( ${SLURM_NTASKS:-16} * ${SLURM_CPUS_PER_TASK:-16} ))          # 256
CORES_ACTIVE=$(( ${SLURM_NTASKS:-16} * (${SLURM_CPUS_PER_TASK:-16} - 1) ))   # 240 (1 idle core/rank)

echo "============================================================================="
echo "TPV104 200 m GTS ORDER BASELINE — SeisSol side | job ${SLURM_JOB_ID:-local}"
echo "============================================================================="
echo "=== nodes=$SLURM_NNODES ranks=$SLURM_NTASKS cpus/rank=$SLURM_CPUS_PER_TASK"
echo "=== cores: ${CORES_ALLOC} allocated / ${CORES_ACTIVE} active (1 core per rank reserved and"
echo "===        left IDLE: OMP_NUM_THREADS=cpus-1 with SEISSOL_COMMTHREAD=0).  Publish the"
echo "===        per-element-update metric with BOTH denominators — MFEM uses all 256."
echo "=== dir : $RUN_DIR"
echo "=== legs: o{$LEGS}   decks: $(for o in $LEGS; do printf '%s ' "$(par_for "$o")"; done)"

# --- deck + input pre-flights -------------------------------------------------
for ORD in $LEGS; do
    P="$(par_for "$ORD")"
    [ -f "$P" ] || { echo "FATAL: deck $P not found in $RUN_DIR"; ls -1; exit 2; }
done
if [ ! -f "$MESH" ]; then
    echo "FATAL: $MESH missing (103 MiB, over GitHub's limit — it travels by rsync):"
    echo "  rsync -avP ~/Downloads/tpv104/tpv104.puml.h5 <expanse>:$RUN_DIR/"
    exit 2
fi
if [ "${BENCH_MD5:-1}" = "1" ]; then
    MD5="$(md5sum "$MESH" 2>/dev/null | cut -d' ' -f1)"
    if [ -n "$MD5" ] && [ "$MD5" != "$MESH_MD5_PIN" ]; then
        echo "FATAL: $MESH md5 $MD5 != pinned $MESH_MD5_PIN"
        echo "       R1 (same mesh both codes) VOID — this must be the conversion of the"
        echo "       MFEM tpv104_200m.msh the partner job runs."
        exit 2
    fi
    [ -z "$MD5" ] && echo "WARN: md5sum unavailable — mesh integrity pin SKIPPED"
    echo "=== mesh md5: ${MD5:-<skipped>}  (pin $MESH_MD5_PIN; $NELEM tets, $NDRFACE DR faces)"
fi
for f in tpv104_fault.yaml tpv104_material.yaml tpv104_receivers.dat tpv104_faultreceivers.dat; do
    [ -f "$f" ] || { echo "FATAL: missing $f"; exit 2; }
done

# Deck contract, applied to EVERY leg's deck (not just one).  Strip '!' comments
# before grepping — the decks' own headers MENTION these keys and an unstripped
# grep would false-positive on every one of them.
par_assert () {   # $1 = regex, $2 = human description of the contract it enforces
    grep -qiE "$1" "$PARSTRIP" || { echo "FATAL[deck $PARNAME]: $2   (regex: $1)"; PAR_RC=2; }
}
check_par () {    # $1 = .par path -> 0 if it satisfies the whole baseline contract
    PARNAME="$1"
    PARSTRIP="$(mktemp)"; sed 's/!.*//' "$PARNAME" > "$PARSTRIP"
    PAR_RC=0
    par_assert '^[[:space:]]*ClusteredLTS[[:space:]]*=[[:space:]]*1([^0-9]|$)'  "R1: ClusteredLTS must be 1 (GTS). 2 = rate-2 LTS and voids the GTS-vs-GTS baseline."
    par_assert '^[[:space:]]*CFL[[:space:]]*=[[:space:]]*0\.5([^0-9]|$)'        "R1: CFL must be 0.5 (matched to the MFEM decks)."
    par_assert '^[[:space:]]*EndTime[[:space:]]*=[[:space:]]*2\.0([^0-9]|$)'    "R1: EndTime must be 2.0 s."
    par_assert '^[[:space:]]*t_0[[:space:]]*=[[:space:]]*1\.0([^0-9]|$)'        "R1: t_0 must be 1.0 s (parity with the MFEM T_nuc)."
    par_assert '^[[:space:]]*OutputPointType[[:space:]]*=[[:space:]]*0([^0-9]|$)' "R4: OutputPointType must be 0 (kills &Elementwise AND &Pickpoint; both write inside the compute stopwatch)."
    par_assert '^[[:space:]]*SurfaceOutput[[:space:]]*=[[:space:]]*0([^0-9]|$)' "R4: SurfaceOutput must be 0 — it also enables the per-timestep PGV fold inside the compute region."
    par_assert '^[[:space:]]*ReceiverOutput[[:space:]]*=[[:space:]]*0([^0-9]|$)' "R4: ReceiverOutput must be EXPLICITLY 0 — the key defaults to TRUE."
    par_assert '^[[:space:]]*EnergyOutput[[:space:]]*=[[:space:]]*0([^0-9]|$)'  "R4: EnergyOutput must be 0 — it changes the DR friction KERNEL, not just IO."
    par_assert '^[[:space:]]*EnergyTerminalOutput[[:space:]]*=[[:space:]]*0([^0-9]|$)' "R4: EnergyTerminalOutput must be 0."
    par_assert '^[[:space:]]*Checkpoint[[:space:]]*=[[:space:]]*0([^0-9]|$)'    "R4: Checkpoint must be EXPLICITLY 0 (the key defaults to TRUE in v1.3.x)."
    par_assert '^[[:space:]]*Format[[:space:]]*=[[:space:]]*10([^0-9]|$)'       "R4: Format must be 10 (no wavefield output)."
    par_assert '^[[:space:]]*iOutputMask[[:space:]]*='                          "iOutputMask must be present — it is read with readOrFail even when output is disabled."
    par_assert "^[[:space:]]*MeshFile[[:space:]]*=[[:space:]]*'?${MESH}'?"      "R1: MeshFile must be $MESH."
    # Keys that MUST NOT be present (they would re-introduce sync points / LTS).
    for bad in 'LtsWiggleFactorMin' 'LtsAutoMergeClusters' 'LtsMaxNumberOfClusters'; do
        grep -qiE "^[[:space:]]*${bad}[[:space:]]*=" "$PARSTRIP" && {
            echo "FATAL[deck $PARNAME]: $bad present — under ClusteredLTS=1 it is inert at best"
            echo "             and can trip the isWiggleFactorValid logError.  Remove it."; PAR_RC=2; }
    done
    rm -f "$PARSTRIP"
    return $PAR_RC
}

DECK_RC=0
REF_BODY=""
for ORD in $LEGS; do
    P="$(par_for "$ORD")"
    check_par "$P" || DECK_RC=2
    # IDENTICAL-BODY CONTRACT (R4): every leg's deck must be byte-identical once
    # the '!' header is stripped, so that ORDER is the only variable on the
    # SeisSol side.  Compared against the first requested leg's body.
    B="$(grep -v '^!' "$P" | md5sum | cut -d' ' -f1)"
    if [ -z "$REF_BODY" ]; then REF_BODY="$B"; REF_NAME="$P"
    elif [ "$B" != "$REF_BODY" ]; then
        echo "FATAL[deck]: $P body differs from $REF_NAME once '!' comments are stripped."
        echo "             The three decks must differ ONLY in their header — otherwise the"
        echo "             order legs do not run the same physics/IO configuration."
        diff <(grep -v '^!' "$REF_NAME") <(grep -v '^!' "$P") | head -20 | sed 's/^/             /'
        DECK_RC=2
    fi
done
[ $DECK_RC -eq 0 ] || { echo "=== deck pre-flights FAILED — refusing to run."; exit $DECK_RC; }
echo "=== deck pre-flights OK for all legs (GTS, CFL 0.5, EndTime 2.0, all output OFF,"
echo "===   checkpoint off) and all requested decks share one byte-identical body (md5 ${REF_BODY})"

# --- binaries: existence + SOURCE-COMMIT parity (fairness gate) ---------------
exe_for () { echo "$SEISSOL_ROOT/SeisSol/build-rome-o${1}-elastic/seissol-elastic-o${1}-f64"; }
commit_of () { command -v strings >/dev/null 2>&1 && strings "$1" 2>/dev/null | grep -m1 -E '^[0-9a-f]{40}$'; }
version_of () { command -v strings >/dev/null 2>&1 && strings "$1" 2>/dev/null | grep -m1 -E '^v1\.[0-9]+\.[0-9]+-'; }

echo "=== binary inventory"
RUN_LEGS=""
MISSING=0
REF_COMMIT=""
for ORD in $LEGS; do
    B="$(exe_for "$ORD")"
    if [ -x "$B" ]; then
        C="$(commit_of "$B")"; V="$(version_of "$B")"
        printf "  OK   o%s  %s\n       version %s  commit %s\n" "$ORD" "$B" "${V:-?}" "${C:-?}"
        [ -z "$REF_COMMIT" ] && REF_COMMIT="$C"
        if [ -n "$C" ] && [ -n "$REF_COMMIT" ] && [ "$C" != "$REF_COMMIT" ]; then
            echo "  FATAL: o${ORD} was built from commit $C but another leg's binary is $REF_COMMIT."
            echo "         The study would be comparing SOURCE REVISIONS, not orders.  Rebuild all"
            echo "         orders from one commit (build_seissol_o2_o3_expanse.sbatch pins"
            echo "         dc6db65133a999eb38c578c14cc443df1000dfd0)."
            exit 3
        fi
        RUN_LEGS="$RUN_LEGS $ORD"
    else
        printf "  MISS o%s  %s\n" "$ORD" "$B"
        echo "       REFUSING this leg.  Build it first with the CO-LOCATED build job"
        echo "       (it is the only copy carrying the CMakeCache ORDER-only fairness diff,"
        echo "        the ORDER:STRING assertion and the stale-binary check):"
        echo "         scp <repo>/miniapps/seas/jobs/bench_order/tpv104_gts_order_baseline/seissol/build_seissol_o2_o3_expanse.sbatch \\"
        echo "             czhao1@login.expanse.sdsc.edu:$SEISSOL_ROOT/"
        echo "         ssh czhao1@login.expanse.sdsc.edu 'cd $SEISSOL_ROOT && sbatch build_seissol_o2_o3_expanse.sbatch'"
        MISSING=1
    fi
done
[ -n "$RUN_LEGS" ] || { echo "FATAL: no requested leg has a binary — nothing to run."; exit 2; }

# --- toolchain (MUST match the build) ----------------------------------------
module purge
module load slurm cpu/0.17.3b gcc/10.2.0
module load openmpi/4.1.3 || module load openmpi/4.1.1
module list 2>&1

# --- runtime env (verbatim from the proven SAFS / Phase-5 SeisSol decks) -----
# HDF5 >=1.10 POSIX file locking FAILS on Lustre (errno 524) and kills PUML IO.
export HDF5_USE_FILE_LOCKING=FALSE
# 1 core per rank is reserved for SeisSol's communication thread -> OMP gets N-1.
# R5: this is the VALIDATED shape.  Do NOT "fix" the idle core by setting 16 —
# the 240-vs-256 spread is finding M4 and is reported, not removed.
export OMP_NUM_THREADS=$((SLURM_CPUS_PER_TASK - 1))
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export ASYNC_MODE=SYNC
export ASYNC_BUFFER_ALIGNMENT=8388608
export XDMFWRITER_ALIGNMENT=8388608
export XDMFWRITER_BLOCK_SIZE=8388608
export SEISSOL_ASAGI_MPI_MODE=OFF
# comm thread off: confirmed fix for the 64-rank LTS ghost timeout.
export SEISSOL_COMMTHREAD=0

# --- UCX robustness (ride out fabric hiccups under MPI_THREAD_MULTIPLE) ------
export UCX_RC_MLX5_TIMEOUT=35000000.00us
export UCX_RC_MLX5_RNR_TIMEOUT=35000000.00us
export UCX_RC_MLX5_RETRY_COUNT=180
export UCX_RC_MLX5_RNR_RETRY_COUNT=180
export UCX_DC_MLX5_TIMEOUT=35000000.00us
export UCX_DC_MLX5_RNR_TIMEOUT=35000000.00us
export UCX_DC_MLX5_RETRY_COUNT=180
export UCX_DC_MLX5_RNR_RETRY_COUNT=180
export UCX_UD_MLX5_TIMEOUT=35000000.00us
export UCX_UD_MLX5_RETRY_COUNT=180

# *** SLURM 23.02 DOES NOT propagate --cpus-per-task from sbatch to srun. ***
# Without BOTH the env var and the explicit -c below each rank gets ONE cpu and
# the run crawls at ~1/15 speed.
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
echo "=== OMP_NUM_THREADS=$OMP_NUM_THREADS  SRUN_CPUS_PER_TASK=$SRUN_CPUS_PER_TASK"
numactl -H 2>/dev/null | head -3

# --- pre-registered per-order predictions (R1/R2 gates) ----------------------
# dt = CFL * 2*insphere / (vp * (2*ConvergenceOrder - 1)), nothing clamps it for
# elastic.  Scaled exactly by 7/(2o-1) from the verified o4 value.
dt_pred     () { case "$1" in 2) echo 8.542987e-04;; 3) echo 5.125792e-04;; 4) echo 3.661280e-04;; esac; }
steps_pred  () { case "$1" in 2) echo 2342;;         3) echo 3902;;         4) echo 5463;;        esac; }
modes_of    () { case "$1" in 2) echo 4;;            3) echo 10;;           4) echo 20;;          esac; }
mfem_partner() { case "$1" in 2) echo p1;;           3) echo p2;;           4) echo p3;;          esac; }

SCORE="$RUN_DIR/scoreboard_ss.${SLURM_JOB_ID:-local}.tsv"
: > "$SCORE"

# -----------------------------------------------------------------------------
# extract_leg <order> <log>   — every R6 quantity, plus the fairness gates.
# -----------------------------------------------------------------------------
extract_leg () {
    ORD="$1"; LOG="$2"
    echo ""
    echo "--- [bench] LEG o${ORD} readout (MFEM partner $(mfem_partner "$ORD")) ---------------"
    [ -f "$LOG" ] || { echo "  (no log)"; return; }

    # (1) GTS gates — three independent log signatures.
    NGTS=$(grep -c 'GTS has been selected' "$LOG")
    echo "  'GTS has been selected.' x${NGTS}   $( [ "$NGTS" -ge 1 ] && echo OK || echo 'FAIL -> ClusteredLTS != 1, leg VOID' )"
    grep -m1 'Theoretical speedup to GTS' "$LOG" | sed 's/^/  /'
    NHIST=$(awk '/Cluster histogram/{f=1;next} f&&/[0-9]+ : [0-9]+ , [0-9]+/{n++;next} f&&n>0{exit} END{print n+0}' "$LOG")
    echo "  cluster histogram rows: ${NHIST}  $( [ "$NHIST" = "1" ] && echo 'OK (single cluster = GTS)' || echo 'FAIL -> more than one cluster, leg VOID' )"
    grep -A2 'Cluster histogram' "$LOG" | sed 's/^/    /' | head -4

    # (2) exactly ONE epoch (no periodic output => no sync-point dt clipping).
    NEP=$(grep -c 'Start simulation epoch' "$LOG")
    echo "  simulation epochs: ${NEP}  $( [ "$NEP" = "1" ] && echo 'OK (IO fully off)' || echo 'FAIL[R4] -> a periodic output survived; steps are clipped' )"
    grep -m1 'End simulation epoch' "$LOG" | sed 's/^/  /'
    grep -qiE 'Peak ground motion \(PGV\) enabled' "$LOG" && \
        echo "  ERROR[R4]: PGV sampling is ON — it folds over every free-surface face at EVERY local timestep INSIDE the compute region."

    # (3) elements (the histogram row prints the raw int; 'Cell count:' uses ' separators).
    NEL=$(awk '/Cluster histogram/{f=1;next} f&&match($0,/[0-9]+ : [0-9]+ , [0-9]+/){s=substr($0,RSTART,RLENGTH); split(s,a," "); print a[3]; exit}' "$LOG")
    grep -m1 'Cell count:' "$LOG" | sed 's/^/  /'
    if [ -n "$NEL" ] && [ "$NEL" != "$NELEM" ]; then
        echo "  ERROR[R1]: element count $NEL != pinned $NELEM — NOT the same mesh as the MFEM legs."
    else
        echo "  elements: ${NEL:-$NELEM} (pin $NELEM) OK"
    fi

    # (4) dt — 'Minimum timestep: 366.1280 µs'.  µ is UTF-8 (0xC2 0xB5); map it to
    #     'u' so the unit is parseable under LC_ALL=C.
    DTLINE="$(grep -m1 'Minimum timestep' "$LOG" | sed -e 's/\xc2\xb5/u/g' -e 's/µ/u/g')"
    echo "  ${DTLINE:-Minimum timestep <absent>}"
    DT=$(printf '%s\n' "$DTLINE" | awk '{
        for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?$/) { v=$i+0; u=$(i+1);
            m = (u ~ /^us/) ? 1e-6 : (u ~ /^ms/) ? 1e-3 : (u ~ /^ns/) ? 1e-9 : (u ~ /^s/) ? 1 : 0;
            if (m > 0) { printf "%.6e", v*m; exit } } }')
    DTP="$(dt_pred "$ORD")"
    if [ -n "$DT" ]; then
        awk -v a="$DT" -v b="$DTP" -v o="$ORD" 'BEGIN{
            d=(a-b)/b; if(d<0) d=-d;
            if (d < 5e-3) printf "  dt GATE OK   : %g s vs pre-registered %g s (%.3f%%)\n", a, b, 100*d;
            else          printf "  dt GATE FAIL : %g s vs pre-registered %g s (%.3f%%) -> o%d leg is VOID (wrong binary? drifted mesh/CFL?)\n", a, b, 100*d, o;
        }'
    fi
    grep -m1 -iE 'wiggle' "$LOG" | sed 's/^/  WARN: wiggle line present (must not appear under GTS): /' || true

    # (5) THE wall: the TERMINAL accumulated compute stopwatch.  The PER-EPOCH
    #     'Time spent this epoch (compute)' line uses split() and double-counts
    #     that epoch's IO — never scrape it.
    for k in 'Simulation time (total)' 'Simulation time (compute)' 'Simulation time (blocking IO)' 'Time spent for initial IO'; do
        grep -m1 -F "$k" "$LOG" | sed 's/^/  /'
    done
    # Stopwatch::print (Monitoring/Stopwatch.cpp:69-76) emits
    #   <mean> (per rank: (<mean> ± <std>) s ; range: [<min>, <max>])
    # so the LEADING value is the MEAN over ranks.  MFEM's Caliper wall is a MAX
    # over ranks.  Quoting SeisSol-mean against MFEM-max would credit SeisSol with
    # its whole load imbalance (~6.4% on this problem) — the same magnitude as the
    # 240-vs-256 core caveat.  Both statistics are therefore extracted and the
    # scoreboard publishes both; use max-vs-max as the headline (a bulk-
    # synchronous step ends when the LAST rank arrives).
    TCOMP=$(grep -m1 -F 'Simulation time (compute)' "$LOG" | sed 's/.*per rank:[^0-9-]*//' \
            | grep -oE '^[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?' | head -1)
    TCOMPMAX=$(grep -m1 -F 'Simulation time (compute)' "$LOG" \
            | sed -E 's/.*range:[[:space:]]*\[[^,]*,[[:space:]]*//' \
            | grep -oE '^[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?' | head -1)
    echo "  T_compute statistic: mean=${TCOMP:-?} s   max-over-ranks=${TCOMPMAX:-?} s"
    echo "    (MFEM's Caliper wall is a MAX; publish max-vs-max AND mean-vs-mean, never mixed.)"
    TTOT=$(grep -m1 -F 'Simulation time (total)' "$LOG" | sed 's/.*per rank:[^0-9-]*//' \
            | grep -oE '^[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?' | head -1)
    TIO=$(grep -m1 -F 'Simulation time (blocking IO)' "$LOG" | sed 's/.*per rank:[^0-9-]*//' \
            | grep -oE '^[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?' | head -1)
    echo "  parsed: T_compute=${TCOMP:-?} s  T_total=${TTOT:-?} s  T_blockingIO=${TIO:-?} s"
    if [ -n "$TCOMP" ] && [ -n "$TTOT" ]; then
        awk -v t="$TTOT" -v c="${TCOMP:-0}" -v i="${TIO:-0}" 'BEGIN{
            r = t - c - i; if (r<0) r=-r;
            if (t>0 && r/t > 0.05) printf "  WARN: total - compute - blockingIO = %.1f s (%.1f%% of total) — check the parse.\n", r, 100*r/t;
        }'
    fi

    # (6) steps: SeisSol never prints a total.  DERIVE FROM THE MEASURED dt, not
    #     from the pre-registered prediction — if the dt gate above failed (wrong
    #     binary, drifted mesh/CFL, a surviving output module clipping steps at
    #     sync points), using the prediction would still produce a plausible-
    #     looking us-core/update number against a fictional step count.
    if [ -n "$DT" ]; then
        STP=$(awk -v dt="$DT" 'BEGIN{ printf "%d", int(2.0/dt) + ((2.0/dt > int(2.0/dt)) ? 1 : 0) }')
        echo "  steps = ceil(EndTime/dt_measured) = ${STP}   (pre-registered: $(steps_pred "$ORD"))"
    else
        STP=""
        echo "  steps: CANNOT DERIVE — no 'Minimum timestep' line was parsed.  R6 is unavailable"
        echo "         for this leg; refusing to substitute the pre-registered step count."
    fi
    LASTPROG="$(grep -E 'updates since sync' "$LOG" | tail -1)"
    echo "  last progress line: ${LASTPROG:-<none>}"
    NLAST=$(printf '%s\n' "$LASTPROG" | grep -oE 'sync: *[0-9]+' | grep -oE '[0-9]+' | head -1)
    if [ -n "$NLAST" ] && [ -n "$STP" ]; then
        # The progress line fires every 100 steps under GTS, so the LAST one must
        # sit in [derived-100, derived] — strictly BELOW the total, never above.
        # A symmetric |diff| < 100 test would accept an over-run; at o2 that
        # window is 4.3% of the whole step count, i.e. 4.3% on the metric.
        awk -v a="$NLAST" -v b="$STP" 'BEGIN{ d=b-a;
            if (d >= 0 && d < 100) printf "  step audit OK: last printed %d, derived total %d (gap %d, expected 0..99)\n", a, b, d;
            else                   printf "  step audit FAIL: last printed %d vs derived %d (gap %d outside 0..99) — truncated run, or a sync point clipped steps\n", a, b, d; }'
    fi

    # (7) free cross-checks: per-element kernel regression, imbalance, FLOP.
    echo "  --- LoopStatistics / FLOP cross-checks (independent of the wall above) ---"
    grep -m1 'Load imbalance' "$LOG" | sed 's/^/    /'
    grep -m1 'Total time spent in compute kernels' "$LOG" | sed 's/^/    /'
    grep -E '\( per element \)|\(total time\)' "$LOG" | sed 's/^/    /' | head -12
    grep -E 'Total calculated (HW|NZ)-FLOP|NZ part of HW-FLOP' "$LOG" | sed 's/^/    /'
    echo "    NOTE: '( per element )' is a per-RANK-TEAM (15 threads) kernel-only slope — use it"
    echo "          as an internal consistency check, NEVER compared directly to an MFEM number."
    echo "          Recompute FLOP/s yourself against T_compute; SeisSol's printed HW-FLOP/s uses"
    echo "          compute+blockingIO as its divisor."

    # (8) the R6 metric — both core denominators AND both rank statistics.
    if [ -n "$TCOMP" ] && [ -n "$STP" ]; then
        awk -v t="$TCOMP" -v tx="${TCOMPMAX:-0}" -v st="$STP" -v ne="$NELEM" \
            -v ca="$CORES_ALLOC" -v cc="$CORES_ACTIVE" -v md="$(modes_of "$ORD")" 'BEGIN{
            upd = st*ne;
            printf "  R6 (MEAN over ranks): T %.1f s | %.1f s/sim-s | %.3f us-core/upd @%d | %.3f @%d | /mode %.4f\n",
                   t, t/2.0, t*ca*1e6/upd, ca, t*cc*1e6/upd, cc, t*ca*1e6/upd/md;
            if (tx > 0)
              printf "  R6 (MAX  over ranks): T %.1f s | %.1f s/sim-s | %.3f us-core/upd @%d | %.3f @%d | /mode %.4f   <- pair with MFEM Caliper Max\n",
                   tx, tx/2.0, tx*ca*1e6/upd, ca, tx*cc*1e6/upd, cc, tx*ca*1e6/upd/md;
        }'
    else
        echo "  R6: NOT COMPUTED — missing T_compute and/or a measured step count.  This leg"
        echo "      contributes no baseline point; do not substitute predicted values."
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$ORD" "$(mfem_partner "$ORD")" "${DT:-NA}" "${STP:-NA}" "${TCOMP:-NA}" "${TTOT:-NA}" "$(modes_of "$ORD")" "${TCOMPMAX:-NA}" >> "$SCORE"

    grep -m1 -iE 'nan|error|abort' "$LOG" | sed 's/^/  WARNING: /' || true
}

# -----------------------------------------------------------------------------
# MAIN — one leg per binary, each in its own directory so the (few) files SeisSol
# still opens cannot collide between legs.  The .par and all inputs are SHARED by
# symlink: the three legs are byte-identically configured by construction.
# -----------------------------------------------------------------------------
RC_TOTAL=0
[ "$MISSING" = "1" ] && RC_TOTAL=5     # a refused leg must not look like success

for ORD in $RUN_LEGS; do
    EXE="$(exe_for "$ORD")"
    PAR="$(par_for "$ORD")"
    LEGDIR="$RUN_DIR/leg_o${ORD}.${SLURM_JOB_ID:-local}"
    LOG="$RUN_DIR/run_o${ORD}.${SLURM_JOB_ID:-local}.log"
    mkdir -p "$LEGDIR/outputs"          # SeisSol will NOT create the output dir
    for f in "$PAR" "$MESH" tpv104_fault.yaml tpv104_material.yaml tpv104_receivers.dat tpv104_faultreceivers.dat; do
        ln -sf "$RUN_DIR/$f" "$LEGDIR/$f"
    done
    # Deck <-> binary provenance: the deck's own header names the order it is for.
    # A mismatch means the wrong .par was forced onto this leg (BENCH_SS_PAR), and
    # since the bodies are identical it would otherwise be undetectable except by
    # the pre-registered dt.  Warn loudly; the dt gate is the hard check.
    case "$PAR" in
        *_o${ORD}.par) : ;;
        *) echo "WARN[provenance]: leg o${ORD} is running deck '$PAR', whose name does not"
           echo "     match the order.  The BODY is contract-checked and identical across decks,"
           echo "     so the physics is right, but this leg's log header will name another order."
           echo "     The authoritative order check is the 'Minimum timestep' gate below." ;;
    esac
    # SeisSol takes the deck path as argv[1]; run with cwd = the leg dir so every
    # relative path in the .par (yaml, mesh, outputs/) resolves inside it.
    echo ""
    echo "============================================================================="
    echo "[bench] LEG o${ORD}  (R2 partner: MFEM $(mfem_partner "$ORD"))"
    echo "        exe : $EXE"
    echo "        deck: $PAR"
    echo "        dir : $LEGDIR"
    echo "        pred: dt $(dt_pred "$ORD") s   steps $(steps_pred "$ORD")   modes/elem $(modes_of "$ORD")"
    echo "============================================================================="
    echo "=========================== o${ORD} START $(date) ==========================="
    START=$SECONDS
    ( cd "$LEGDIR" && \
      srun --mpi=pmix -n "$SLURM_NTASKS" -c "$SLURM_CPUS_PER_TASK" --cpu-bind=cores --unbuffered \
           stdbuf -oL -eL "$EXE" "$PAR" ) > "$LOG" 2>&1
    rc=$?
    ELAPSED=$((SECONDS - START))
    echo "=========================== o${ORD} END   $(date)  exit=$rc  (${ELAPSED}s job wall) ==="
    echo "NOTE: ${ELAPSED}s is the JOB wall (PUML read + partition + init included).  It is a"
    echo "      FOOTNOTE, not the baseline number — R6 uses 'Simulation time (compute)'."
    [ $rc -eq 0 ] || RC_TOTAL=3
    extract_leg "$ORD" "$LOG"
    # R4 audit: with every output module off, outputs/ must stay essentially empty.
    NOUT=$(ls -1 "$LEGDIR/outputs" 2>/dev/null | wc -l | tr -d ' ')
    echo "  outputs/ file count: ${NOUT}  $( [ "$NOUT" = "0" ] && echo 'OK (fully IO-free)' || echo '<- inspect; something still wrote' )"
    ls -lh "$LEGDIR/outputs" 2>/dev/null | head -8
done

# -----------------------------------------------------------------------------
# SCOREBOARD (R3 + R6)
# -----------------------------------------------------------------------------
echo ""
echo "============================================================================="
echo "[bench] SeisSol GTS ORDER SCOREBOARD — TPV104 200 m, ${NELEM} tets"
echo "  T_loop = terminal 'Simulation time (compute)' (NOT the job wall, NOT the per-epoch"
echo "  split() line, which double-counts IO).  BOTH rank statistics are shown: SeisSol's"
echo "  line leads with the MEAN, while MFEM's Caliper wall is a MAX — pair like with like."
echo "  steps are derived from the MEASURED dt, not from the pre-registered value."
echo "============================================================================="
if [ -s "$SCORE" ]; then
    awk -v NE="$NELEM" -v CA="$CORES_ALLOC" -v CC="$CORES_ACTIVE" -F'\t' '
    BEGIN { printf "%-4s %-8s %-13s %-7s %-11s %-11s %-12s %-14s %-14s %-14s\n",
                   "o", "partner", "dt [s]", "steps", "Tmean [s]", "Tmax [s]", "s / sim-s",
                   "us-core@'"$CORES_ALLOC"'", "us-core@'"$CORES_ACTIVE"'", "us-core(max)@'"$CORES_ALLOC"'"; }
    {
        o=$1; part=$2; dt=$3; st=$4; tc=$5; modes=$7; tx=$8;
        upd=st*NE;
        ua=(upd>0 && tc+0>0)? tc*CA*1e6/upd : 0;
        uc=(upd>0 && tc+0>0)? tc*CC*1e6/upd : 0;
        ux=(upd>0 && tx+0>0)? tx*CA*1e6/upd : 0;
        printf "%-4s %-8s %-13s %-7s %-11.1f %-11.1f %-12.1f %-14.3f %-14.3f %-14.3f\n",
               o, part, dt, st, tc, tx+0, tc/2.0, ua, uc, ux;
        w[o]=tc+0; s[o]=st+0;
    }
    END {
        print "";
        print "  Within-code order ratios: the (2o-1) step penalty is exact and known; only the";
        print "  per-step factor carries information.";
        if (w[2] > 0) for (k=2;k<=4;k++) if (w[k] > 0) {
            printf "    o%d: total %.2fx o2", k, w[k]/w[2];
            if (s[k]>0 && s[2]>0) printf " ; per-step %.2fx o2 (steps %.2fx)", (w[k]/s[k])/(w[2]/s[2]), s[k]/s[2];
            printf "\n";
        }
    }' "$SCORE"
else
    echo "  (no completed legs)"
fi
echo ""
echo "Raw scoreboard TSV (o, mfem_partner, dt, steps, T_compute_mean, T_total, modes, T_compute_max): $SCORE"
echo "CROSS-CODE: pair each row with the MFEM p(o-1) row from ../run_mfem_gts_orders.sbatch."
echo "  MATCH THE STATISTIC: MFEM Caliper Max <-> the 'us-core(max)' column here; MFEM Caliper"
echo "  Avg <-> the mean columns.  Mixing them credits SeisSol with its ~6% load imbalance."
echo "  Publish BOTH core denominators (256 allocated for both codes = billing-fair headline;"
echo "  240 active here vs 256 there = the ~6.7% idle-core self-handicap, finding M4)."
echo "  Footnote the o2/p1 pair as the least reliable: at 4 modes SeisSol pads vector lanes,"
echo "  so its HW-FLOP >> NZ-FLOP there (read the NZ-vs-HW ratio printed per leg)."
echo "  This compute figure is NOT comparable to the historical 111 s/sim-s (that run had"
echo "  EnergyOutput=1, PGV sampling, per-step pickpoints and 0.005 s receivers all inside"
echo "  the compute stopwatch)."
echo ""
if [ $RC_TOTAL -ne 0 ]; then
    echo "[bench] FINISHED WITH FAILURES / REFUSED LEGS (rc=$RC_TOTAL).  Logs: $RUN_DIR/run_o*.log"
    exit $RC_TOTAL
fi
echo "[bench] done — all requested SeisSol legs completed."
exit 0
