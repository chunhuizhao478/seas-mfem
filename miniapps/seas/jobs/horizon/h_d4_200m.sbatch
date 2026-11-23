#!/bin/bash
#SBATCH -J h_d4_200m
#SBATCH -o h_d4_200m_%j.out
#SBATCH -e h_d4_200m_%j.err
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 06:00:00
#SBATCH -A EAR26024

# GPU dynamic-rupture plan D4 (plan rev2 §4.5) on the 200 m meshes, run here
# with the D2 device step (the fault on the host, one rank per GPU) for the
# GPU-against-CPU speed of a full run: PROB on the 200 m gold mesh, order P,
# 0 to 12 s, with the committed deck, --device cuda --matrix-free 1, on the
# 4 GB200 of one Horizon node (or -N NODES -n NP, one rank per GPU).  Then
# the stations are compared
#   - with the Vista CPU gold of the same case (81759da, global stepping,
#     576 to 2304 ranks: tpv205 1076378 / 1075916, tpv104 1075919 / 1075917)
#     by scripts/nitsche_gate_check.py: the relative L2 error of every station
#     field, reported (how far round-off grows through the rupture);
#   - with the SCEC references by scripts/compare_nitsche_fault_stations.py,
#     as in that gold's README, the gold alongside (needs numpy in python3;
#     skipped, and said, when it has none).
#
#   sbatch jobs/horizon/h_d4_200m.sbatch PROB P   (PROB tpv205 or tpv104, P 1
#      or 2; from miniapps/seas of the Horizon clone after h_d2_build passed
#      at the same commit; or H_CHAIN="build d4" bash jobs/horizon/h_d2_submit.sh)
#
# DEVICE_FAULT=1 adds --device-fault 1 (GPU plan D3: the fault on the GPU
# too, the first full-GPU timing); the job then also fails when the summary
# does not say the fault ran on the device, and the output folder is
# d4_200m_<PROB>_p<P>_df1_<jobid>.
# The mesh is the clone's, or else MESH_SRC's (h_dyn_env.sh) read in place,
# and must have the gold's sha256 (D4_MESH_SHA256).  CPU_GOLD=<folder of a CPU run's
# station files and driver.log> compares with that run instead of the gold.
# The job fails when the run fails, when the driver does not say the device
# step ran, or when the station files do not line up with the gold's
# (station set, rows, times); the errors and the SCEC tests are reported,
# not gated.
# Cost (inferred from the 1000 m gates scaled by tets, fault faces and
# steps): TPV205 minutes at p = 1 and under an hour at p = 2; TPV104, whose
# rate-and-state fault runs on the 4 host cores, up to about 4 h at p = 2.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PROB=${1:?usage: sbatch jobs/horizon/h_d4_200m.sbatch PROB P}
P=${2:?usage: sbatch jobs/horizon/h_d4_200m.sbatch PROB P}
case "$PROB" in
   tpv205) MESH=tpv205/mesh/tpv2053d_200m.msh
           SCEC=(--stations tpv205/stations/tpv205_fault_stations.csv
                 --ref drdg3d_100m=tpv205/benchmark_data/DRDG3D_100m_O4
                 --ref drdg3d_200m=tpv205/benchmark_data/DRDG3D_200m_O4 --peak-min-stations 14) ;;
   tpv104) MESH=tpv104/mesh/tpv104_200m.msh
           SCEC=(--stations tpv104/stations/tpv104_fault_stations.csv
                 --ref drdg3d=tpv104/benchmark_data/DRDG3D --ref seissol=tpv104/benchmark_data/seisol) ;;
   *) echo "ERROR: PROB is tpv205 or tpv104, not '$PROB'"; exit 2 ;;
esac
case "$P" in 1|2) ;; *) echo "ERROR: P is 1 or 2, not '$P'"; exit 2 ;; esac
DECK=$PROB/configs/${PROB}_nitsche.toml
MF_LAYOUT=${MF_LAYOUT:-auto}
DEVICE_FAULT=${DEVICE_FAULT:-0}
case "$DEVICE_FAULT" in 0|1) ;; *) echo "ERROR: DEVICE_FAULT is 0 or 1, not '$DEVICE_FAULT'"; exit 2 ;; esac
DF_TAG=""; [ "$DEVICE_FAULT" = 0 ] || DF_TAG="_df1"
ROOT=$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)              # submitted from miniapps/seas
# shellcheck source=jobs/horizon/h_dyn_env.sh
source "$SLURM_SUBMIT_DIR/jobs/horizon/h_dyn_env.sh" || exit 1
# shellcheck source=jobs/vista/d1_env.sh
source "$SLURM_SUBMIT_DIR/jobs/vista/d1_env.sh"
NP=${SLURM_NTASKS:?run with sbatch -n NP}
NGPU=$(nvidia-smi -L | grep -c '^GPU ')
RPN=$(( (NP + ${SLURM_NNODES:-1} - 1) / ${SLURM_NNODES:-1} ))
[ "$RPN" -le "$NGPU" ] || { echo "ERROR: one rank per GPU: $NP tasks on ${SLURM_NNODES:-1} node(s) of $NGPU GPU(s)"; exit 2; }
T0=$(date +%s)

d1_guard "$ROOT"
[ "$(git -C "$ROOT" rev-parse HEAD)" = "$(git -C "$ROOT" rev-parse '@{u}' 2>/dev/null)" ] \
   || { echo "ERROR: $ROOT HEAD is not its upstream; pull first"; exit 2; }
grep -q '^MFEM_USE_CUDA *= YES' "$ROOT/config/config.mk" \
   || { echo "ERROR: $ROOT is not a CUDA build (h_d2_build.sbatch)"; exit 2; }
d1_check_build "$ROOT" BUILD_COMMIT_CUDA
SHA=$(git -C "$ROOT" rev-parse HEAD)
cd "$ROOT/miniapps/seas" || exit 2
MDIR=.; [ -f "$MESH" ] || MDIR=$MESH_SRC        # read in place: parallel jobs copy nothing
d4_check_mesh "$MDIR" "$MESH" || exit 2
if [ -n "${CPU_GOLD:-}" ]; then
   GOLD=$CPU_GOLD
else
   mapfile -t _golds < <(ls -d "$PROB/gold/nitsche-kv0p1_p${P}_cd_pc_200m_jobid"*_git81759da 2>/dev/null)
   [ "${#_golds[@]}" = 1 ] || { echo "ERROR: want one $PROB p$P 200 m CPU gold, found ${#_golds[@]}"; exit 2; }
   GOLD=${_golds[0]}/results
fi
[ -f "$GOLD/driver.log" ] || { echo "ERROR: no CPU run in $GOLD (driver.log)"; exit 2; }
OUT=$SCRATCH/gpu_dyn/d4_200m_${PROB}_p${P}${DF_TAG}_$SLURM_JOB_ID
mkdir -p "$OUT/run"
d1_gh_modules
{ echo "job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ), $NP rank(s) on ${SLURM_NNODES:-1} node(s), commit $SHA ($ROOT)"
  echo "case $PROB p$P, deck $DECK, mesh $MDIR/$MESH $(cd "$MDIR" && sha256sum "$MESH" | cut -d' ' -f1), --mf-layout $MF_LAYOUT, --device-fault $DEVICE_FAULT"
  echo "CPU gold $GOLD"
  echo "modules: $(module -t list 2>&1 | tr '\n' ' ')"
  nvidia-smi --query-gpu=name,compute_mode,driver_version --format=csv,noheader; } | tee "$OUT/RUN_INFO.txt"
cp "$DECK" "$OUT/"

O=$OUT/run
t=$(date +%s)
ibrun -n "$NP" ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MDIR/$MESH" --order "$P" --tfinal 12.0 \
      --output-dir "$O" --log-every 1000 --device cuda --matrix-free 1 --mf-layout "$MF_LAYOUT" \
      --device-fault "$DEVICE_FAULT" > "$O/driver.log" 2>&1
rc=$?
echo $rc > "$O/exit_code"
echo "RUN $PROB 200 m p$P on $NP GPU(s): exit $rc in $(( $(date +%s) - t )) s"
grep -a -E "=== done|time in ode.Step|fault work:" "$O/driver.log" | sed 's/^ */   /'
grep -a -E "=== done|time in ode.Step|fault work:" "$GOLD/driver.log" | sed 's/^ */   CPU gold: /'
fail=0
if [ "$rc" != 0 ] || ! grep -a -q "^=== done: " "$O/driver.log"; then
   echo "RUN: FAIL (exit $rc; last lines of $O/driver.log:)"; tail -20 "$O/driver.log"; fail=1
elif ! grep -a -q "device step on" "$O/driver.log"; then
   echo "RUN: FAIL (the summary does not say the device step ran)"; fail=1
elif [ "$DEVICE_FAULT" = 1 ] && ! grep -a -q "fault on the device: " "$O/driver.log"; then
   echo "RUN: FAIL (the summary does not say the fault ran on the device)"; fail=1
fi

if [ "$fail" = 0 ]; then
   python3 scripts/nitsche_gate_check.py --run-dir "$O" --run-prefix "$PROB" \
           --golden-dir "$GOLD" --golden-prefix "$PROB" --json "$OUT/vs_cpu_gold.json" > "$OUT/vs_cpu_gold.txt" 2>&1
   crc=$?
   sed 's/^/   /' "$OUT/vs_cpu_gold.txt"
   python3 - "$OUT/vs_cpu_gold.json" <<'EOF' | sed 's/^/   /'
import json, sys
st = json.load(open(sys.argv[1]))["errors"]["stations"]
fields = next(iter(st.values())).keys() if st else []
for f in fields:
    s, e = max(((k, v[f]) for k, v in st.items()), key=lambda x: x[1])
    print(f"largest relative L2 error of {f}: {e:.3e} ({s})")
EOF
   if [ "$crc" = 0 ]; then echo "STATIONS vs CPU gold: compared (report only)"
   else echo "STATIONS vs CPU gold: FAIL (the files do not line up, $OUT/vs_cpu_gold.txt)"; fail=1; fi

   if python3 -c 'import numpy' 2>/dev/null; then
      python3 -I scripts/compare_nitsche_fault_stations.py "${SCEC[@]}" --prefix "$PROB" \
              --run "gpu_p${P}_200m=$O" --gold "vista_p${P}_200m=$GOLD:nitsche" \
              --json "$OUT/compare_summary.json" > "$OUT/compare.txt" 2> "$OUT/compare.err"
      echo "SCEC comparison (exit $?): $OUT/compare.txt"
      sed -n '/^# median noise RMS/,/^# wrote/p' "$OUT/compare.txt" | grep -v '^# wrote' | sed 's/^/   /'
   else
      echo "SCEC comparison: NOT RUN (python3 has no numpy here); run compare_nitsche_fault_stations.py on $O elsewhere"
   fi
fi

echo "D4 200 m $PROB p$P, $NP GPU(s), --device-fault $DEVICE_FAULT: $([ $fail = 0 ] && echo DONE || echo FAIL); commit $SHA; out $OUT; total $(( $(date +%s) - T0 )) s"
exit $fail
