#!/bin/bash
#SBATCH --job-name=tpv104_gts_order_mfem
#SBATCH --account=usc143             # SAME account as the SeisSol partner job — running the two
                                     # halves of one baseline on different projects is an avoidable
                                     # asymmetry in queue priority, and lbl107 (the old default) had
                                     # only ~4.8k SU left on 2026-07-24 against ~1.2-1.8k SU for
                                     # these three legs alone, expiring 2026-08-27.  usc143 had
                                     # ~291k SU free; ddp408 ~719k.  Override: sbatch -A ddp408 ...
#SBATCH --partition=compute
#SBATCH --constraint="lustre"        # match the SeisSol side (both runs live on /expanse/lustre)
#SBATCH --nodes=2                    # R5: MFEM's VALIDATED shape — 2 nodes x 128 pure-MPI ranks
#SBATCH --ntasks-per-node=128        #     = 256 ranks x 1 core.  DO NOT change: it is the shape
#SBATCH --cpus-per-task=1            #     every prior MFEM TPV104 number was measured at.
#SBATCH --mem=140000M                # TPV104 has NO sidecars
#SBATCH -t 08:00:00                  # p1 ~0.2h + p2 ~0.9h + p3 ~3.6h + 3x setup; see COST below
#SBATCH -o tpv104_gts_order_mfem.%j.out
#SBATCH -e tpv104_gts_order_mfem.%j.out
#SBATCH --export=ALL

# =============================================================================
# TPV104 200 m  ---  GTS ORDER SPEED BASELINE, MFEM SIDE (p1, p2, p3)
#   *** EXPANSE (SDSC).  Partner job: seissol/run_seissol_gts_orders.sbatch ***
#
# WHAT THIS IS ----------------------------------------------------------------
#   A COST baseline at NOMINALLY-MATCHED ORDER.  Three MFEM legs (p1/p2/p3, each
#   with ADER order p+1) on the SCEC TPV104 200 m production mesh (2,464,689
#   tets), GLOBAL time stepping, run to t = 2.0 s, with ALL output OFF.  The
#   SeisSol partner job runs o2/o3/o4 on the SAME mesh, same CFL, same window.
#   NO ACCURACY CLAIM IS MADE — there is no error functional in this experiment.
#
# WHY TWO SBATCH FILES (not one job for all six legs) --------------------------
#   1. RANK SHAPE IS AN ALLOCATION-TIME PROPERTY.  MFEM's validated shape is
#      --ntasks-per-node=128 --cpus-per-task=1 (256 pure-MPI ranks); SeisSol's is
#      --ntasks-per-node=8 --cpus-per-task=16 (16 ranks x 15 OMP + 1 reserved
#      core).  A single #SBATCH header can declare only one, and the SeisSol
#      prologue DERIVES OMP_NUM_THREADS and SRUN_CPUS_PER_TASK from
#      $SLURM_CPUS_PER_TASK (SLURM 23.02 does not propagate -c to srun) — under a
#      -c 1 allocation that silently gives OMP_NUM_THREADS=0.  Fairness rule R5
#      forbids "fixing" that by changing either shape.
#   2. Different PMI: MFEM is proven under `srun --mpi=pmi2`, SeisSol under
#      `srun --mpi=pmix` (MPI_THREAD_MULTIPLE + comm-thread machinery).
#   3. Different filesystem roots: MFEM runs out of the seas-mfem repo tree;
#      the SeisSol deck + its 103 MiB tpv104.puml.h5 are rsync'd to
#      /expanse/lustre/projects/lbl107/czhao1/safs/... and are NOT in the repo.
#   4. Different accounts and wildly different walltimes (~5 h vs ~1 h).  One
#      combined job means a walltime kill or a crash in either code destroys the
#      other code's measurements — and a killed job produces NO Caliper report.
#
# FAIRNESS CONTRACT ENFORCED HERE (R1-R6) -------------------------------------
#   R1  Same mesh (md5-pinned), EndTime 2.0 s, CFL 0.5, GTS ([numerics].lts="off"
#       asserted per deck;  SeisSol side asserts ClusteredLTS = 1).
#   R2  ACCURACY-MATCHED PAIRING  p_k <-> o_(k+1):  p1<->o2, p2<->o3, p3<->o4.
#       SeisSol's convention is o = p+1 (o modes = (o)(o+1)(o+2)/6 = MFEM's p
#       modes), so at a matched pair BOTH codes advance 4 / 10 / 20 modes per
#       element x 9 fields, with the SAME dt and the SAME step count.  Nothing in
#       the MFEM code ties [mesh].order to [numerics].ader_order, so this sbatch
#       ASSERTS ader_order == order + 1 per leg (exit 2 otherwise).
#   R3  dt is NOT held equal across orders — each order gets its own CFL dt
#       (dt ~ 1/(2p+1)).  Both s/sim-s AND per-element-update us-core are printed.
#   R4  IO OFF: no --paraview, no --paraview-fault-vtu, no --checkpoint-every.
#       The deck must set [output].paraview_fault = "off", paraview_free_surface
#       = "off" (that one is NOT under the paraview_enabled master gate and
#       DEFAULTS TO "vtu" @0.05 s — spatial_dyn_driver.cpp:3865-3871 vs the gate
#       at :3628-3632), [problem].tag = "" (a "tpv104" tag wires the per-step
#       SCEC station writer at :4467 -> :5126 INSIDE the timed scope, flushing 9
#       files per step on <=9 of 256 ranks, with a row count that GROWS with
#       order — while the SeisSol legs write nothing at all in their compute
#       region), and checkpoint_every_steps > nsteps.  Enforced by assert_deck;
#       BENCH_STRICT_IO=0 downgrades the tag check to a warning (NOT recommended
#       — the default is 1 and a tag-on leg is not comparable to the SeisSol row).
#       The output dir is scanned afterwards for .vtu/.pvd/.vtkhdf (must be ZERO)
#       and for cp_checkpoint_r*.txt (must be EXACTLY nranks: the driver's
#       post-loop final checkpoint at :5265-5274 ALWAYS fires and cannot be
#       disabled, but it is outside the Caliper scope; an EXCESS means an in-loop
#       write did land inside the timed region).
#   R5  Rank shape = MFEM's validated 2 x 128 pure MPI (see #SBATCH above).
#   R6  Per leg this job records: Caliper wall for seas::spatial_dyn::step,
#       steps, dt, elements, and derives us-core per element-update.
#
# THE ONE THAT SILENTLY VOIDS EVERYTHING --------------------------------------
#   [numerics].cfl_dg_safety defaults to 3.0 (spatial/code/spatial_friction.hpp).
#   At 3.0 the run takes 3x the steps and MFEM looks 3x slower than SeisSol while
#   every log line still looks sane.  This job (a) requires `cfl_dg_safety = 1.0`
#   IN the TOML, (b) passes --cfl-dg-safety 1.0 on the CLI anyway, and (c)
#   hard-checks the printed `[time] dt` against the pre-registered value.
#
# PRE-REGISTERED PREDICTIONS (fail the leg if the log disagrees) ---------------
#   dt(p1) = 8.542987e-04 s  nsteps = 2342   element-updates = 5.7723e9
#   dt(p2) = 5.125792e-04 s  nsteps = 3902   element-updates = 9.6172e9
#   dt(p3) = 3.661280e-04 s  nsteps = 5463   element-updates = 1.34646e10
#   (from the verified p3 dt at CFL 0.5, scaled exactly by 7/(2p+1);
#    identical to SeisSol's GlobalTimestep formula at cfl_dg_safety = 1.0)
#
# THE BINARY MUST BE BUILT WITH CALIPER — MFEM_PERF_SCOPE compiles to NOTHING
# without it (mfem general/annotation.hpp:31) and the driver has NO other timing
# instrumentation, so R6 would be unsatisfiable and the SLURM wall (which
# includes a 126 MB .msh read + serial partition) is NOT a substitute:
#   USE_CALIPER=YES bash <root>/build_expanse.sh
# This job EXITS if libcaliper is not linked (the Phase-5 template only warned).
#
# COST / SCHEDULING -----------------------------------------------------------
#   Projected from the measured p3 GTS anchor (2.366 s/step, 6465 s/sim-s, job
#   52344266) and the p1:p2:p3 = 1 : 2.74 : 7.45 per-step ratio:
#     p1 ~744 s, p2 ~3395 s, p3 ~12925 s of stepping (~4.8 h total) + ~3 x setup.
#   Legs run CHEAPEST FIRST so a p3 overrun cannot destroy p1/p2.  To split into
#   three jobs (recommended if the queue is tight):
#     sbatch --export=ALL,BENCH_LEGS=1 -t 1:00:00 run_mfem_gts_orders.sbatch
#     sbatch --export=ALL,BENCH_LEGS=2 -t 2:00:00 run_mfem_gts_orders.sbatch
#     sbatch --export=ALL,BENCH_LEGS=3 -t 6:00:00 run_mfem_gts_orders.sbatch
#
# ENV OVERRIDES ---------------------------------------------------------------
#   BENCH_LEGS       (default "1 2 3")  which orders to run, in this order.
#   BENCH_CFG_P1/2/3 (default tpv104_gts_p<N>.toml in this folder)
#   BENCH_TFINAL     (default 2.0)  SECONDS, BARE NUMBER — the CLI parser is
#                    std::stod (drivers/spatial_dyn_driver.cpp:157-162), which
#                    SILENTLY DISCARDS a unit suffix: "500ms" would become 500 s.
#                    The decks carry the unit-aware tfinal = "2.0s"; this flag
#                    only re-states it.  Do not change for a baseline.
#   BENCH_LEVERS     (default "--deriv-cache --shared-ck-recursion --face-cache")
#                    ONE lever string for ALL legs (the lever choice is itself a
#                    first-order result: p3 measured 6465 s/sim-s levers-partial
#                    vs 4930 with --face-cache).  Set "" for a levers-off control.
#   BENCH_STRICT_IO  (default 1)  a non-empty [problem].tag is FATAL.  Set 0 only
#                    for a deliberately tag-on control leg, whose number is then
#                    NOT comparable to the SeisSol row.
#   BENCH_MD5        (default 1)  0 => skip the 126 MB mesh md5 pin.
#   BENCH_OUT_BASE   (default <folder>/runs)   job-ID-keyed output root.  POINT
#                    THIS AT LUSTRE SCRATCH: the unavoidable post-loop final
#                    checkpoint writes 256 ASCII files per leg, ~0.7 / 1.8 /
#                    3.6 GB at p1 / p2 / p3 (~6 GB for all three).  e.g.
#                    BENCH_OUT_BASE=/expanse/lustre/scratch/$USER/temp_project/tpv104_gts
# =============================================================================

set -o pipefail
export LC_ALL=C
export LANG=C

# --- Expanse modules (MUST match build_expanse.sh; verbatim from the proven
#     jobs/comm_dev/a0_wiggle_tpv104_200m_expanse prologue) ---------------------
module purge
module load slurm 2>/dev/null || true
module load cpu/0.17.3b
module load gcc/10.2.0
for _mpi in openmpi/4.1.5 openmpi/4.1.3 openmpi/4.1.1 openmpi; do
    if module load "${_mpi}" 2>/dev/null && module is-loaded "${_mpi}" 2>/dev/null; then
        echo "OpenMPI module: ${_mpi}"; break
    fi
done

# --- resolve the seas-mfem root (walk UP to miniapps/seas/Makefile) ----------
# This sbatch lives at <root>/miniapps/seas/jobs/bench_order/<folder>/.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SEAS_MFEM_ROOT="${SEAS_MFEM_ROOT:-${SUBMIT_DIR}}"
_walk="${SEAS_MFEM_ROOT}"
for _ in 0 1 2 3 4 5 6 7; do
    if [[ -f "${_walk}/miniapps/seas/Makefile" ]]; then
        SEAS_MFEM_ROOT="$(cd "${_walk}" && pwd)"; break
    fi
    _walk="${_walk}/.."
done
if [[ ! -f "${SEAS_MFEM_ROOT}/miniapps/seas/Makefile" ]]; then
    echo "ERROR: cannot locate seas-mfem root from ${SUBMIT_DIR}."
    echo "  Set SEAS_MFEM_ROOT via sbatch --export=ALL,SEAS_MFEM_ROOT=/path/to/seas-mfem ..."
    exit 1
fi

# --- runtime shared libs (caliper FIRST — required here; then hdf5/zfp/petsc) -
export LD_LIBRARY_PATH="${SEAS_MFEM_ROOT}/extern/caliper/install/lib64:${SEAS_MFEM_ROOT}/extern/caliper/install/lib:${SEAS_MFEM_ROOT}/extern/hdf5/install/lib:${SEAS_MFEM_ROOT}/extern/zfp/install/lib64:${SEAS_MFEM_ROOT}/extern/zfp/install/lib:${SEAS_MFEM_ROOT}/extern/petsc/install/lib:${LD_LIBRARY_PATH}"
export HDF5_PLUGIN_PATH="${SEAS_MFEM_ROOT}/extern/h5z-zfp/install/plugin"
# HDF5 POSIX file locking FAILS on Lustre (errno 524).
export HDF5_USE_FILE_LOCKING=FALSE

cd "${SEAS_MFEM_ROOT}"
grep -q '^MFEM_USE_PETSC *= YES'  config/config.mk || { echo "ERROR: no PETSc in config.mk"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES'  config/config.mk || { echo "ERROR: no MUMPS in config.mk"; exit 1; }

cd miniapps/seas
DRIVER=./seas_spatial_dyn_driver
if [[ ! -x "${DRIVER}" ]]; then
    echo "ERROR: ${DRIVER} not found.  Build it FIRST on Expanse (SDSC forbids login-node compiles):"
    echo "   USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 1
fi

# --- Caliper is MANDATORY (R6: it is the only wall-clock instrument) ----------
if ! ldd "${DRIVER}" 2>/dev/null | grep -qi 'libcaliper'; then
    echo "ERROR: the binary lacks libcaliper.  MFEM_PERF_SCOPE compiles to a NO-OP"
    echo "       without it, so there is NO time-loop wall at all and R6 cannot be"
    echo "       satisfied (the SLURM wall includes mesh read + partition and biases"
    echo "       the cheap p1 leg far more than p3).  Rebuild:"
    echo "         USE_CALIPER=YES bash ${SEAS_MFEM_ROOT}/build_expanse.sh"
    exit 4
fi
echo "Caliper:     ON  $(ldd "${DRIVER}" 2>/dev/null | grep -m1 -i libcaliper | awk '{print $1" -> "$3}')"

# --- inputs (all paths resolved from miniapps/seas) --------------------------
JOB_REL="jobs/bench_order/tpv104_gts_order_baseline"
MESH_LN="${JOB_REL}/mesh_tpv104_200m.msh"     # the decks reference this via [mesh].path
MESH_MD5_PIN="e88ae223205530548b1d3438e6eb17c1"   # tpv104_200m.msh (the .puml.h5 partner is
                                                  # 5911b0b4ec29a8f61968d5b32ed0a11c — see the
                                                  # SeisSol deck; both are the SAME 2,464,689 tets)
NELEM=2464689           # md5-pinned element count (R6 denominator; the GTS path does
                        # not print a global element count, so it is pinned, not scraped)
NDRFACE=37904

# Re-create the in-folder mesh symlink if a fresh checkout dropped it (*.msh is gitignored).
if [[ ! -e "${MESH_LN}" ]]; then
    ln -sf ../../../tpv104/mesh/tpv104_200m.msh "${MESH_LN}" 2>/dev/null || true
fi
if [[ ! -e "${MESH_LN}" ]]; then
    echo "ERROR: mesh missing / dangling symlink: $(pwd)/${MESH_LN} -> $(readlink "${MESH_LN}" 2>/dev/null)"
    echo "       The 200 m TPV104 mesh (tpv104/mesh/tpv104_200m.msh, ~126 MB) is *.msh and"
    echo "       GITIGNORED — it does NOT arrive via git clone/pull.  scp/rsync it into"
    echo "       tpv104/mesh/ on Expanse before submitting."
    exit 1
fi
if [[ "${BENCH_MD5:-1}" == "1" ]]; then
    _md5="$(md5sum "${MESH_LN}" 2>/dev/null | cut -d' ' -f1)"
    if [[ -n "${_md5}" && "${_md5}" != "${MESH_MD5_PIN}" ]]; then
        echo "ERROR: mesh md5 ${_md5} != pinned ${MESH_MD5_PIN}."
        echo "       R1 (same mesh both codes) is VOID — the SeisSol partner runs the"
        echo "       pumgen conversion of exactly this file.  Refusing to run."
        exit 1
    fi
    [[ -z "${_md5}" ]] && echo "WARN: md5sum unavailable — mesh integrity pin SKIPPED"
    echo "Mesh md5:    ${_md5:-<skipped>}  (pin ${MESH_MD5_PIN}; ${NELEM} tets, ${NDRFACE} DR faces)"
fi

# --- knobs -------------------------------------------------------------------
LEGS="${BENCH_LEGS:-1 2 3}"          # cheapest first
TFINAL="${BENCH_TFINAL:-2.0}"        # BARE SECONDS — --tfinal is std::stod'd and would
                                     # silently drop a unit suffix (spatial_dyn_driver.cpp:157)
LEVERS="${BENCH_LEVERS---deriv-cache --shared-ck-recursion --face-cache}"
STRICT_IO="${BENCH_STRICT_IO:-1}"    # default FATAL: an MFEM leg with the per-step station
                                     # writer on is not IO-comparable to the SeisSol legs.

# R-101 shared-fault consistency tripwire: an MPI collective that fires at step 0
# and then every 100th step while t <= T_nuc (= 1.0 s), i.e. 12 / 20 / 28 times at
# p1 / p2 / p3 — an ORDER-DEPENDENT global sync inside the timed region with no
# SeisSol counterpart.  TPV104's planar fault has |n.n_ref| = 1 everywhere so the
# check is a provable no-op here (drivers/spatial_dyn_driver.cpp:5147-5150).
# Exported UNIFORMLY for all three legs and recorded in the banner.
export SEAS_R101_SKIP=1
CORES="${SLURM_NTASKS:-256}"         # ALLOCATED = ACTIVE for MFEM (pure MPI, no idle core).
                                     # SeisSol's partner leg allocates 256 but works 240
                                     # (OMP_NUM_THREADS=15 + SEISSOL_COMMTHREAD=0) — the
                                     # cross-code table must publish BOTH denominators.

JOB_ABS="${SEAS_MFEM_ROOT}/miniapps/seas/${JOB_REL}"
OUT_BASE="${BENCH_OUT_BASE:-${JOB_ABS}/runs}"
OUT_ROOT="${OUT_BASE}/${SLURM_JOB_ID:-local}"
mkdir -p "${OUT_ROOT}"
SCORE="${OUT_ROOT}/scoreboard.tsv"
: > "${SCORE}"

echo "============================================================================="
echo "TPV104 200 m GTS ORDER BASELINE — MFEM side | job ${SLURM_JOB_ID:-local}"
echo "============================================================================="
echo "Git commit:  $(git rev-parse HEAD 2>/dev/null)"
echo "Cluster:     Expanse (SDSC)   Launcher: srun --mpi=pmi2"
echo "Shape (R5):  ${SLURM_NNODES:-2}N x ${SLURM_NTASKS_PER_NODE:-128}r = ${CORES} ranks, 1 core/rank, pure MPI"
echo "             -> ${NELEM} tets / ${CORES} ranks = $((NELEM / CORES)) elements per rank"
echo "Legs:        p in {${LEGS}}   (ADER order = p+1; R2 accuracy-matched to SeisSol o=p+1)"
echo "Window (R1): tfinal=${TFINAL}   CFL=0.5 (deck)   cfl_dg_safety=1.0 (deck + CLI)   GTS (lts=off)"
echo "Levers:      ${LEVERS:-<none>}   (IDENTICAL on all legs; the lever set is part of the result)"
echo "IO (R4):     ParaView OFF (fault + free surface), no in-loop checkpoint, [problem].tag"
echo "             must be empty (STRICT_IO=${STRICT_IO}).  The post-loop FINAL checkpoint is"
echo "             unavoidable but sits OUTSIDE the timed scope; it is counted, not forbidden."
echo "Metric (R6): Caliper region seas::spatial_dyn::step.  BOTH Max and Avg across ranks are"
echo "             recorded: MFEM's Max is the honest bulk-synchronous wall, but SeisSol's"
echo "             'Simulation time (compute):' LEADS WITH THE MEAN — so publish the cross-code"
echo "             ratio max-vs-max (headline) AND mean-vs-mean, never mixed."
echo "Tripwires:   SEAS_R101_SKIP=1 exported on ALL legs (removes an order-dependent MPI"
echo "             collective that fires 12/20/28x at p1/p2/p3).  The per-step Q.Norml2() +"
echo "             MPI_Allreduce NaN guard (spatial_dyn_driver.cpp:5190-5199) CANNOT be"
echo "             disabled and remains a disclosed MFEM-only per-step global sync."
echo "Output root: ${OUT_ROOT}"
echo "SCOPE:       COST baseline at NOMINALLY-matched order.  GTS only, one mesh, one"
echo "             node count, three order points.  NO accuracy claim, no scaling study."
echo "CAVEAT:      pmesh.SetCurvature([mesh].order) (drivers/spatial_dyn_driver.cpp:1733) gives"
echo "             p3 a 20-node H1 nodal transformation on a STRAIGHT-SIDED tet mesh where p1"
echo "             gets 4 — an order-scaling geometry tax SeisSol (always affine) does not pay."
echo "             --face-cache removes it on non-fault interior faces only."

# -----------------------------------------------------------------------------
# Per-order constants + deck assertions.
# -----------------------------------------------------------------------------
leg_dt_pred () { case "$1" in 1) echo 8.542987e-04;; 2) echo 5.125792e-04;; 3) echo 3.661280e-04;; esac; }
leg_nsteps_pred () { case "$1" in 1) echo 2342;; 2) echo 3902;; 3) echo 5463;; esac; }
leg_modes () { case "$1" in 1) echo 4;; 2) echo 10;; 3) echo 20;; esac; }   # (p+1)(p+2)(p+3)/6

# assert_deck <order> <ader> <cfg>  — every gate that, if missed, voids a rule.
assert_deck () {
    local ord="$1" ader="$2" cfg="$3" rc=0
    [[ -f "${cfg}" ]] || {
        echo "ERROR: deck not found: $(pwd)/${cfg}"
        echo "       Decks expected in ${JOB_REL}/ (override with BENCH_CFG_P${ord}=<path>).  Folder holds:"
        ls -1 "${JOB_REL}" 2>/dev/null | sed 's/^/         /'
        return 2; }

    local v
    v="$(grep -E '^[[:space:]]*order[[:space:]]*=' "${cfg}" | head -1 | sed -E 's/[^0-9]*([0-9]+).*/\1/')"
    [[ "${v}" == "${ord}" ]] || { echo "ERROR[R2]: ${cfg} sets [mesh].order=${v:-<absent>}, expected ${ord}"; rc=2; }
    v="$(grep -E '^[[:space:]]*ader_order[[:space:]]*=' "${cfg}" | head -1 | sed -E 's/[^0-9]*([0-9]+).*/\1/')"
    [[ "${v}" == "${ader}" ]] || { echo "ERROR[R2]: ${cfg} sets ader_order=${v:-<absent>}, expected ${ader} (= order+1;"
                                   echo "           NOTHING in the driver ties the two — a desync silently voids the pairing)"; rc=2; }
    grep -qE '^[[:space:]]*cfl_dg_safety[[:space:]]*=[[:space:]]*1\.0' "${cfg}" \
        || { echo "ERROR[R1]: ${cfg} must set cfl_dg_safety = 1.0 — the struct default is 3.0"
             echo "           (spatial/code/spatial_friction.hpp:150) which triples the step count"
             echo "           and inflates MFEM 3x against SeisSol.  Refusing to run."; rc=2; }
    grep -qE '^[[:space:]]*cfl[[:space:]]*=[[:space:]]*0\.5([^0-9]|$)' "${cfg}" \
        || { echo "ERROR[R1]: ${cfg} must set [numerics].cfl = 0.5 (matched to SeisSol CFL=0.5)"; rc=2; }
    grep -qE '^[[:space:]]*lts[[:space:]]*=[[:space:]]*"off"' "${cfg}" \
        || { echo "ERROR[R1]: ${cfg} must set [numerics].lts = \"off\" (this baseline is GTS-vs-GTS)"; rc=2; }
    grep -qE '^[[:space:]]*mixed_flux[[:space:]]*=[[:space:]]*"none"' "${cfg}" \
        || { echo "ERROR: ${cfg} must set mixed_flux = \"none\" (pure upwind; \"adjacent\" ABORTS"
             echo "       under ADER at ComputeMaxDt, dynamic/wave_operator.inl:7717)"; rc=2; }
    grep -qE '^[[:space:]]*time_integrator[[:space:]]*=[[:space:]]*"ader"' "${cfg}" \
        || { echo "ERROR[R1]: ${cfg} must set [numerics].time_integrator = \"ader\" EXPLICITLY."
             echo "           The RK path uses RkCflFactor = 3/(2p+1) instead of 1/(2p+1)"
             echo "           (spatial/code/spatial_friction.hpp:767 vs :744), i.e. dt x3 and 3x"
             echo "           FEWER steps than SeisSol — while 'cfl_dg_safety = 1' still prints."; rc=2; }
    grep -qE '^[[:space:]]*paraview_fault[[:space:]]*=[[:space:]]*"off"' "${cfg}" \
        || { echo "ERROR[R4]: ${cfg} must set [output].paraview_fault = \"off\" for an IO-free wall"; rc=2; }
    # The free-surface slice is the one output collection that is NOT under the
    # paraview_enabled master gate and whose DEFAULT IS ON ("vtu" @ 0.05 s):
    # spatial_dyn_driver.cpp:3865-3871 builds fs_out independently of the gate at
    # :3628-3632, and paraview_write (:5121, inside the timed scope) saves it at
    # :4325 above the `if (!fault_wants) return;` guard.  The driver still prints
    # "ParaView output: OFF" while it writes (:3960 vs the "  FreeSurf:" line at
    # :3966), so no log gate can catch this — only the deck gate can.
    grep -qE '^[[:space:]]*paraview_free_surface[[:space:]]*=[[:space:]]*"off"' "${cfg}" \
        || { echo "ERROR[R4]: ${cfg} must set [output].paraview_free_surface = \"off\"."
             echo "           It is NOT under the paraview_enabled master gate and DEFAULTS to"
             echo "           \"vtu\" @ 0.05 s, writing ~40 free-surface VTUs per leg from all"
             echo "           ${CORES} ranks INSIDE seas::spatial_dyn::step, at a per-write cost"
             echo "           that grows with [mesh].order."; rc=2; }
    grep -qE '^[[:space:]]*tfinal[[:space:]]*=[[:space:]]*"2\.0s"' "${cfg}" \
        || { echo "ERROR[R1]: ${cfg} must set [time].tfinal = \"2.0s\" (matched to SeisSol EndTime)"; rc=2; }

    # checkpoint: the in-loop write sits INSIDE the timed region and, at the
    # Phase-5 template's 4000, fires ONLY on the p3 leg (5463 steps) — a
    # full-state parallel write charged to the highest-order leg alone.
    v="$(grep -E '^[[:space:]]*checkpoint_every_steps[[:space:]]*=' "${cfg}" | head -1 | sed -E 's/[^0-9]*([0-9]+).*/\1/')"
    if [[ -n "${v}" ]] && [[ "${v}" -le "$(leg_nsteps_pred "${ord}")" ]]; then
        echo "ERROR[R4]: ${cfg} sets checkpoint_every_steps=${v} <= nsteps=$(leg_nsteps_pred "${ord}") — an"
        echo "           in-loop checkpoint would land inside the timed region.  Leave it unset"
        echo "           (default 10000 exceeds every leg); the parser forbids 0."; rc=2
    fi
    # dt_max must not clip the CFL dt (that would break dt-equality with SeisSol).
    v="$(grep -E '^[[:space:]]*dt_max[[:space:]]*=' "${cfg}" | head -1)"
    [[ -n "${v}" ]] && echo "  note: deck has ${v} (must exceed dt=$(leg_dt_pred "${ord}") s; the run log is checked for clipping)"

    # R4 SYMMETRY: [problem].tag must be EMPTY.  A recognised tag (tpv6 / tpv31 /
    # tpv102 / tpv104 / tpv205) wires a station writer at spatial_dyn_driver.cpp
    # :4466-4470 whose WriteStep + per-station fflush (dynamic/tpv104_stations.hpp
    # :229) runs EVERY step INSIDE the MFEM_PERF_SCOPE opened at :4874, on <=9 of
    # ${CORES} ranks, in a bulk-synchronous loop — with a row count that grows
    # with the step count (2342 / 3902 / 5463 at p1/p2/p3).  The SeisSol partner
    # legs run OutputPointType = 0 + ReceiverOutput = 0: ZERO writes inside their
    # compute stopwatch.  Keeping the tag is therefore both a cross-code and an
    # order asymmetry, and MFEM has no compute/IO split to subtract it with.
    local tagv
    tagv="$(grep -E '^[[:space:]]*tag[[:space:]]*=' "${cfg}" | head -1 | sed -E 's/.*=[[:space:]]*"([^"]*)".*/\1/')"
    if [[ -n "${tagv}" ]]; then
        if [[ "${STRICT_IO}" == "1" ]]; then
            echo "ERROR[R4]: ${cfg} sets [problem].tag = \"${tagv}\" -> a per-step station writer"
            echo "           runs INSIDE the timed region on <=9 of ${CORES} ranks (flushing 9 files"
            echo "           per step) and writes MORE rows at higher order, while the SeisSol legs"
            echo "           write NOTHING inside their compute stopwatch.  Set tag = \"\"."
            echo "           (BENCH_STRICT_IO=0 downgrades this to a warning for a deliberate"
            echo "            tag-on control leg — whose number is then NOT cross-code comparable.)"; rc=2
        else
            echo "  WARN[R4]: [problem].tag = \"${tagv}\" — per-step SCEC station writes are INSIDE"
            echo "            the timed region.  This leg is a CONTROL, not a baseline point; quote"
            echo "            its delta against the tag-off leg as the station-write tax."
        fi
    else
        echo "  R4 OK: [problem].tag is empty -> no station writer wired; stations_write stays a"
        echo "         null std::function and spatial_dyn_driver.cpp:5126 is a branch-not-taken."
        echo "         Physics evidence for this leg = the V_max progress prints (:5223) and the"
        echo "         final V_max_global (:5288), both stdout-only."
    fi
    return ${rc}
}

# -----------------------------------------------------------------------------
# run_leg <order> <ader> <cfg> <out> <log> <cali_stem>
# -----------------------------------------------------------------------------
run_leg () {
    local ord="$1" ader="$2" cfg="$3" out="$4" log="$5" stem="$6"
    # CALI_CONFIG contains commas/parens -> export it (propagated by --export=ALL)
    # rather than passing it on srun --export (comma-separated, would corrupt it).
    #
    # profile.mpi is DELIBERATELY OMITTED.  It wraps every MPI call inside the
    # timed region with an annotation, and MFEM is decomposed 256 ways against
    # SeisSol's 16 — so the wrapper tax lands disproportionately on MFEM while
    # SeisSol's compute stopwatch (Simulator.cpp:106-108) is a bare start/pause
    # pair with no instrumentation at all.  That is a one-sided overhead in
    # exactly the quantity this baseline measures.  region.count is kept: it is
    # the independent step-count read-back that cross-checks [time] nsteps.
    export CALI_CONFIG="runtime-report(calc.inclusive,region.count,max_column_width=60,output=${stem}.cali-region-report.txt),spot(output=${stem}.cali)"
    echo ""
    echo "[bench] === LEG p${ord} / ADER ${ader}  cfg=${cfg}  tfinal=${TFINAL}  out=${out} ==="
    date
    srun --mpi=pmi2 -n "${CORES}" --export=ALL --unbuffered stdbuf -oL -eL "${DRIVER}" \
        --config "${cfg}" \
        --mixed-flux none \
        --time-integrator ader \
        --ader-order "${ader}" \
        --cfl-dg-safety 1.0 \
        ${LEVERS} \
        --tfinal "${TFINAL}" \
        --print-derived \
        --output-dir "${out}" \
        > "${log}" 2>&1
    local rc=$?
    date
    echo "[bench] leg p${ord} rc=${rc}"
    return ${rc}
}

# --- numeric helpers ----------------------------------------------------------
# First float on a line (handles 3.66128e-04 / 0.000366128 / 5463).
first_num () { sed -E 's/[^0-9eE+.-]+/ /g' <<< "$1" | tr ' ' '\n' \
               | grep -E '^[0-9]+(\.[0-9]*)?([eE][-+]?[0-9]+)?$' | head -1; }
# Caliper runtime-report column layout.  With this CALI_CONFIG the header reads
#   Path   Min time/rank   Max time/rank   Avg time/rank   Time %   <count>
# (verified against document/code_optimization_dev/profiling_analysis_2026-06-24
# .md:26-28, where the step row is 830.60 | 830.66 | 830.65 | 99.87 | 1445).
# We do NOT trust that positionally: cali_ncols counts the numeric tokens on the
# row and cali_header echoes the header, and the caller REFUSES to publish if the
# two disagree with the assumed layout.  Note that if calc.inclusive ever emits a
# separate "(inc)" column group, token #2 would be the EXCLUSIVE max (the loop
# wall minus its annotated children) — a silently understated MFEM number.  That
# is why the layout is asserted rather than assumed.
cali_ncols () {      # $1 = report -> number of numeric tokens on the step row
    local rpt="$1"; [[ -f "${rpt}" ]] || { echo 0; return; }
    awk '/seas::spatial_dyn::step/ { n=0; s=$0;
        while (match(s, /[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?/)) {
            n++; s=substr(s, RSTART+RLENGTH); }
        if (n > best) best=n } END { print best+0 }' "${rpt}"
}
cali_has_inc () {    # 1 if the header advertises a separate inclusive column group
    local rpt="$1"; [[ -f "${rpt}" ]] || { echo 0; return; }
    grep -qiE 'time/rank[[:space:]]*\(inc\)|\(inc\)' "${rpt}" && echo 1 || echo 0
}
cali_step_col () {   # $1 = report, $2 = column index (1=Min 2=Max 3=Avg 4=Time% 5=Count)
    local rpt="$1" col="${2:-2}"; [[ -f "${rpt}" ]] || return 0
    awk -v col="${col}" '/seas::spatial_dyn::step/ {
        n=0; s=$0;
        while (match(s, /[0-9]+\.?[0-9]*([eE][-+]?[0-9]+)?/)) {
            n++; tok[n]=substr(s, RSTART, RLENGTH)+0; s=substr(s, RSTART+RLENGTH);
        }
        if (n >= 3 && tok[2] > bestmax) { bestmax=tok[2]; for(i=1;i<=n && i<=6;i++) best[i]=tok[i]; nb=n }
    } END { if (bestmax > 0 && col <= nb) printf "%.6f\n", best[col] }' "${rpt}"
}

# -----------------------------------------------------------------------------
# extract_leg <order> <ader> <log> <report> <out>  -> appends one scoreboard row
# -----------------------------------------------------------------------------
extract_leg () {
    local ord="$1" ader="$2" log="$3" rpt="$4" out="$5"
    echo ""
    echo "--- [bench] LEG p${ord} readout -------------------------------------------"
    [[ -f "${log}" ]] || { echo "  (no log)"; return; }

    # Banner: the driver's own view of the order pairing (R2 audit trail).
    grep -m1 -E '^fe order:'    "${log}" | sed 's/^/  /'
    grep -m1 -E '^ader order:'  "${log}" | sed 's/^/  /'
    grep -m1 -E '^cfl:'         "${log}" | sed 's/^/  /'

    # R1 gate: cfl_dg_safety must have taken effect.
    if ! grep -q 'cfl_dg_safety = 1 ' "${log}"; then
        echo "  ERROR[R1]: the run did NOT use cfl_dg_safety=1.0 -> timing inflated ~3x:"
        grep -m1 'cfl_dg_safety' "${log}" | sed 's/^/    /'
    fi
    # dt_max clipping would break dt-equality with SeisSol.
    grep -m1 'tightened by user dt_max' "${log}" | sed 's/^/  ERROR[R1]: /' || true

    local dt_line nst_line dt nst dtp nstp
    dt_line="$(grep -m1 -E '^\[time\] dt ' "${log}")"
    nst_line="$(grep -m1 -E '^\[time\] nsteps' "${log}")"
    echo "  ${dt_line:-[time] dt <absent>}"
    echo "  ${nst_line:-[time] nsteps <absent>}"
    grep -m1 -E '^\[time\]   SeisSol-equivalent dt_cfl' "${log}" | sed 's/^/  /' || true
    dt="$(first_num "${dt_line#*=}")"
    nst="$(first_num "${nst_line#*=}")"
    dtp="$(leg_dt_pred "${ord}")"; nstp="$(leg_nsteps_pred "${ord}")"
    if [[ -n "${dt}" ]]; then
        awk -v a="${dt}" -v b="${dtp}" -v p="${ord}" 'BEGIN{
            d = (a-b)/b; if (d<0) d=-d;
            if (d < 5e-3) printf "  dt GATE OK   : %g s vs pre-registered %g s (%.3f%%)\n", a, b, 100*d;
            else          printf "  dt GATE FAIL : %g s vs pre-registered %g s (%.3f%%) -> p%d leg is VOID\n", a, b, 100*d, p;
        }'
    fi
    [[ -n "${nst}" && "${nst}" != "${nstp}" ]] && echo "  WARN: nsteps ${nst} != pre-registered ${nstp}"

    # R6 wall from Caliper.
    echo "  Caliper report: ${rpt}"
    if [[ -f "${rpt}" ]]; then
        grep -m1 -iE 'time/rank|Path' "${rpt}" | sed 's/^/    [header] /'
        grep -m1 -E 'seas::spatial_dyn::step' "${rpt}" | sed 's/^/    [row]    /'
    else
        echo "    MISSING — no Caliper report was written (walltime kill? report is written at clean exit)"
    fi
    local wmin wmax wavg wcnt ncols hasinc
    wmin="$(cali_step_col "${rpt}" 1)"; wmax="$(cali_step_col "${rpt}" 2)"
    wavg="$(cali_step_col "${rpt}" 3)"; wcnt="$(cali_step_col "${rpt}" 5)"
    ncols="$(cali_ncols "${rpt}")"; hasinc="$(cali_has_inc "${rpt}")"
    echo "    parsed: Min=${wmin:-?} Max=${wmax:-?} Avg=${wavg:-?} s/rank   region.count=${wcnt:-?} (should equal nsteps)"
    # Layout assertion: 5 numeric tokens (Min/Max/Avg/Time%/count) and NO separate
    # inclusive column group.  Anything else means the positional map is wrong and
    # the "Max" above may actually be an EXCLUSIVE time (loop wall minus annotated
    # children) — a silently understated MFEM wall.  Refuse to publish, don't guess.
    if [[ "${ncols}" != "5" || "${hasinc}" != "0" ]]; then
        echo "    ERROR[R6]: unexpected Caliper column layout (numeric tokens=${ncols}, expected 5;"
        echo "               separate '(inc)' group present=${hasinc}, expected 0).  The positional"
        echo "               Min/Max/Avg/Time%/count map does NOT hold, so T_loop may be an"
        echo "               EXCLUSIVE time.  Re-read the raw row above and fix cali_step_col"
        echo "               BEFORE publishing any number from this leg."
    else
        echo "    layout OK: 5 numeric tokens, no separate inclusive group -> Min/Max/Avg/Time%/count."
    fi
    echo "    NOTE: 'Max across ranks' is MFEM's honest bulk-synchronous wall.  Avg is recorded"
    echo "          too because SeisSol's 'Simulation time (compute):' leads with the MEAN —"
    echo "          the cross-code ratio must be published max-vs-max AND mean-vs-mean."

    # steps actually executed: prefer Caliper's region.count, else [time] nsteps.
    local steps="${nst}"
    if [[ -n "${wcnt}" ]] && awk -v c="${wcnt}" -v n="${nst:-0}" 'BEGIN{exit !(c>0 && (c-n<3) && (n-c<3))}'; then
        steps="${wcnt%.*}"
    fi

    # Physics sanity (R-soft): p1 at 200 m may simply not rupture like p3.
    local vmax
    vmax="$(awk '/V_max = /{ for(i=1;i<=NF;i++) if($i=="V_max"){v=$(i+2)+0; if(v>m) m=v} } END{ if(m>0) printf "%.4g", m }' "${log}")"
    echo "  V_max(peak over run) = ${vmax:-?} m/s   (p3 reference: ~17.3 m/s at t~1.48 s)"
    grep -m1 -E '^step .*/' "${log}" | sed 's/^/  first progress: /' || true
    tail -3 "${log}" | sed 's/^/  tail: /'
    grep -m1 -iE 'nan|abort|diverg' "${log}" | sed 's/^/  WARNING: /' || true

    # R4 verification, in TWO parts — ParaView artifacts and checkpoints have
    # completely different expectations and must not be lumped together.
    #
    # (a) ParaView: ANY .vtu/.pvd/.vtkhdf means a writer ran inside the timed
    #     region (paraview_write is called at spatial_dyn_driver.cpp:5121, inside
    #     the scope).  Expected count: ZERO.
    local npv
    npv="$(find "${out}" \( -name '*.vtu' -o -name '*.pvd' -o -name '*.vtkhdf' \) 2>/dev/null | wc -l | tr -d ' ')"
    if [[ "${npv}" != "0" ]]; then
        echo "  ERROR[R4]: ${npv} ParaView artifact(s) under ${out} — a writer ran INSIDE the timed region."
        find "${out}" \( -name '*.vtu' -o -name '*.pvd' -o -name '*.vtkhdf' \) 2>/dev/null | head -5 | sed 's/^/    /'
    else
        echo "  R4 OK: no .vtu/.pvd/.vtkhdf under ${out}"
    fi
    # (b) Checkpoints: the driver ALWAYS writes ONE post-loop final checkpoint
    #     (spatial_dyn_driver.cpp:5265-5274; the parser forbids
    #     checkpoint_every_steps = 0, spatial_friction.cpp:1408), producing
    #     exactly ${CORES} files named cp_checkpoint_r<rank>.txt.  That write is
    #     OUTSIDE the Caliper scope and is EXPECTED — flagging it as an R4
    #     violation on every leg would train the reader to ignore this check.
    #     Only an EXCESS is the in-loop-write signature.
    local ncp
    ncp="$(find "${out}" -name 'cp_checkpoint_r*.txt' 2>/dev/null | wc -l | tr -d ' ')"
    if [[ "${ncp}" == "${CORES}" ]]; then
        echo "  R4 OK: ${ncp} cp_checkpoint_r*.txt = exactly one post-loop final checkpoint set"
        echo "         (expected and unavoidable; written AFTER the timed scope closed)."
    elif [[ "${ncp}" == "0" ]]; then
        echo "  WARN: no final checkpoint under ${out} — the run may not have reached the"
        echo "        post-loop section (killed?).  Cross-check '[spatial_dyn] done.' in the log."
    else
        echo "  ERROR[R4]: ${ncp} cp_checkpoint_r*.txt but ${CORES} ranks — an IN-LOOP checkpoint"
        echo "             fired inside the timed region.  Check [output].checkpoint_every_steps"
        echo "             and that no --checkpoint-every reached the CLI."
    fi
    # Free the disk: at p3 this set is ~3.6 GB of ASCII and nothing downstream reads it.
    if [[ "${BENCH_KEEP_CHECKPOINT:-0}" != "1" && "${ncp}" != "0" ]]; then
        find "${out}" -name 'cp_checkpoint_r*.txt' -delete 2>/dev/null
        echo "         (deleted after counting; BENCH_KEEP_CHECKPOINT=1 to retain)"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${ord}" "${ader}" "${dt:-NA}" "${steps:-NA}" "${wmax:-NA}" "${wavg:-NA}" "$(leg_modes "${ord}")" "${vmax:-NA}" \
        >> "${SCORE}"
}

# -----------------------------------------------------------------------------
# MAIN — legs in the given order (cheapest first).
# -----------------------------------------------------------------------------
RC_TOTAL=0
for ORD in ${LEGS}; do
    case "${ORD}" in 1|2|3) ;; *) echo "ERROR: BENCH_LEGS contains '${ORD}'; only 1 2 3 are defined."; exit 1;; esac
    ADER=$(( ORD + 1 ))
    eval "CFG=\${BENCH_CFG_P${ORD}:-${JOB_REL}/tpv104_gts_p${ORD}.toml}"
    echo ""
    echo "============================================================================="
    echo "[bench] LEG p${ORD}  (ADER ${ADER}; R2 partner = SeisSol o${ADER})"
    echo "        deck   : ${CFG}"
    echo "        dt pred: $(leg_dt_pred "${ORD}") s   nsteps pred: $(leg_nsteps_pred "${ORD}")   modes/elem: $(leg_modes "${ORD}")"
    echo "============================================================================="
    if ! assert_deck "${ORD}" "${ADER}" "${CFG}"; then
        echo "[bench] LEG p${ORD} REFUSED — deck gate failed (see above).  Continuing to the next leg."
        RC_TOTAL=2
        continue
    fi
    OUT="${OUT_ROOT}/p${ORD}"
    LOG="${OUT_ROOT}/run_p${ORD}.log"
    mkdir -p "${OUT}"
    run_leg "${ORD}" "${ADER}" "${CFG}" "${OUT}" "${LOG}" "${LOG%.log}"
    rc=$?
    [[ ${rc} -ne 0 ]] && RC_TOTAL=3
    extract_leg "${ORD}" "${ADER}" "${LOG}" "${LOG%.log}.cali-region-report.txt" "${OUT}"
done

# -----------------------------------------------------------------------------
# SCOREBOARD (R3 + R6).  Both throughput and per-element-update cost, because dt
# is NOT held equal across orders and s/sim-s alone is illegible.
# -----------------------------------------------------------------------------
echo ""
echo "============================================================================="
echo "[bench] MFEM GTS ORDER SCOREBOARD — TPV104 200 m, ${NELEM} tets, ${CORES} cores"
echo "  T_loop = Caliper seas::spatial_dyn::step (calc.inclusive requested; the per-leg"
echo "           readout ASSERTS the 5-column Min/Max/Avg/Time%/count layout before this"
echo "           table is believed).  Max across ranks is the headline; Avg is printed too"
echo "           because SeisSol reports its compute wall as a MEAN over ranks."
echo "  us-core/elem-update = 1e6 * T_loop * cores / (steps * elements)."
echo "  /mode divides that by the modes per element (4 / 10 / 20) so the reader can"
echo "  separate dof growth from the extra ADER sub-step (2/3/4) and fault QPs (3/6/12)."
echo "============================================================================="
if [[ -s "${SCORE}" ]]; then
    awk -v NE="${NELEM}" -v C="${CORES}" -F'\t' '
    BEGIN {
        printf "%-4s %-5s %-13s %-7s %-10s %-10s %-11s %-13s %-13s %-10s %-9s\n",
               "p", "ader", "dt [s]", "steps", "Tmax [s]", "Tavg [s]", "s / sim-s",
               "us-core/upd", "us-core(avg)", "/mode", "Vmax";
    }
    {
        p=$1; ader=$2; dt=$3; st=$4; wmax=$5; wavg=$6; modes=$7; vmax=$8;
        sps = (wmax+0 > 0) ? wmax/2.0 : 0;            # 2.0 sim-s window
        upd = st * NE;
        us   = (upd > 0) ? wmax * C * 1e6 / upd : 0;
        usav = (upd > 0 && wavg+0 > 0) ? wavg * C * 1e6 / upd : 0;
        printf "%-4s %-5s %-13s %-7s %-10.1f %-10.1f %-11.1f %-13.3f %-13.3f %-10.4f %-9s\n",
               p, ader, dt, st, wmax, wavg, sps, us, usav, (modes>0? us/modes : 0), vmax;
        w[p]=wmax; s[p]=st;
    }
    END {
        print "";
        print "  Within-code order ratios (R3): total cost splits into the EXACTLY-KNOWN";
        print "  (2p+1) step penalty and the measured per-step work.  The (2p+1)*modes line";
        print "  (12 : 50 : 140, i.e. 4.17x and 11.67x vs p1) is a CRUDE LINEAR-IN-DOF MODEL,";
        print "  NOT a prediction — ADER volume kernels are closer to O(modes^2), so a";
        print "  measured p3/p1 above 11.67x is expected, not a defect.";
        if (w[1] > 0) {
            for (k = 1; k <= 3; k++) if (w[k] > 0) {
                printf "    p%d: total %.2fx p1", k, w[k]/w[1];
                if (s[k] > 0 && s[1] > 0) printf " ; per-step %.2fx p1 (steps %.2fx)", (w[k]/s[k])/(w[1]/s[1]), s[k]/s[1];
                printf "\n";
            }
        }
    }' "${SCORE}"
else
    echo "  (no completed legs)"
fi
echo ""
echo "Raw scoreboard TSV (p, ader, dt, steps, T_max, T_avg, modes, Vmax): ${SCORE}"
echo "STATISTIC (must be stated, not assumed): the headline column above is Max-across-ranks."
echo "  SeisSol's 'Simulation time (compute):' LEADS WITH THE MEAN over its 16 ranks, and its"
echo "  documented load imbalance on this problem is ~6.4%.  Quoting MFEM-max against"
echo "  SeisSol-mean would credit SeisSol with that imbalance — the same magnitude as the"
echo "  240-vs-256 core caveat.  Publish BOTH pairs: max-vs-max (headline; a bulk-synchronous"
echo "  step ends when the LAST rank arrives) and mean-vs-mean beside it."
echo "Cross-code step: pair each row with the SeisSol o(p+1) row from the partner job"
echo "  (seissol/run_seissol_gts_orders.sbatch).  Publish the per-element-update cost"
echo "  with BOTH core denominators: 256 allocated for both codes (billing-fair headline)"
echo "  and 256 MFEM / 240 SeisSol active (SeisSol reserves 1 core/rank and leaves it idle"
echo "  — OMP_NUM_THREADS=15 + SEISSOL_COMMTHREAD=0; that ~6.7% self-handicap is finding M4,"
echo "  not something to 'fix')."
echo "Also disclose: MFEM is decomposed 256 ways vs SeisSol 16 -> ~2.5x more partition"
echo "  interface for MFEM (all MPI), worst at p1.  Per-thread work is matched to 6.7%"
echo "  ($((NELEM / CORES)) elem/rank MFEM vs 10270 elem/thread SeisSol)."

echo ""
if [[ ${RC_TOTAL} -ne 0 ]]; then
    echo "[bench] FINISHED WITH FAILURES (rc=${RC_TOTAL}).  Inspect logs under ${OUT_ROOT}."
    exit ${RC_TOTAL}
fi
echo "[bench] done — all requested MFEM legs completed."
exit 0
