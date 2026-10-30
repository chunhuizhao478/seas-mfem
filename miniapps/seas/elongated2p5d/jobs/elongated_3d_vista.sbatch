#!/bin/bash
#SBATCH -J l2_3d
#SBATCH -o l2_3d_%j.out
#SBATCH -e l2_3d_%j.err
#SBATCH -p gg
#SBATCH -N 16
#SBATCH -n 2304
#SBATCH -t 03:00:00
#SBATCH -A EAR20006
# The 2.5D note's Level 2 test on Vista, part 2 (Chunhui Zhao, 2026-10-30;
# /mnt/project-files/2d-inplane/nonplanar-literature/
# inplane-2p5d-derivation-rev8.pdf, Section 7.2): the 3-D reference
# elongated2p5d/configs/elongated_3d_buried.toml, a buried strike-slip
# rupture L = 120 km by W = 15 km, on the 400 m mesh at p = 2 to 40 s.
# Submit from miniapps/seas/elongated2p5d/jobs/ after vista_l2_build.sbatch
# passed at this code (its load check gives the node memory and dt), e.g.
#   sbatch --dependency=afterok:<build job> --kill-on-invalid-dep=yes elongated_3d_vista.sbatch
# Output: <checkout>-runs/l2_3d_<jobid>: the driver's station files (390
# stations, every 4th step), rupture times, and the fault CSV dumps every
# 4000 steps (about 1 s) for the slip over the whole fault; RUN_INFO.txt
# (commit, mesh sha256, deck); memlog (Committed_AS per node every minute).
# The map-view runs (elongated_2d_*.toml) are on the Mac; the comparison is
# elongated2p5d/scripts/compare_2p5d_3d.py (README.md).
# Estimate (inferred, not measured): the TPV205 200 m p = 2 gold (1.89 M
# tets, 88,725 steps, 16 nodes, 2670 s; job 1075916) scaled by tets (2.41 M)
# and steps (dt about 2.6e-4 s from the cloud's 1000 m dt and the meshes'
# min 6V/A, about 155,000 steps): about 1.7 h on 16 nodes, 9 SU.
# Requested: 16 nodes x 3 h x 0.33 SU = 15.8 SU at most (EAR20006, gg).
# 144 ranks per node: the build job's 400 m load check runs the whole mesh
# on one node.
# Submit only on Chunhui Zhao's typed go for this job.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NVPL_BLAS_NUM_THREADS=1 NVPL_LAPACK_NUM_THREADS=1
MODULES="gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5"
DECK=elongated2p5d/configs/elongated_3d_buried.toml
MESH=elongated2p5d/mesh/elongated_3d_400m.msh
ORDER=2
TFINAL=40.0
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd "$ROOT" || exit 2

# ---- GitHub-first guard (CLAUDE.md), and the binary matches this code ----
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
HEAD=$(git rev-parse HEAD)
BUILT=$(cat miniapps/seas/BUILD_COMMIT_VISTA_NITSCHE 2>/dev/null)
git cat-file -e "$BUILT^{commit}" 2>/dev/null \
   || { echo "ERROR: seas_nitsche_dyn_driver built at '$BUILT', not a commit here; run vista_l2_build.sbatch"; exit 2; }
CODE=$(git diff --name-only "$BUILT" HEAD | grep -v -E '^miniapps/seas/(elongated2p5d/(configs|jobs|scripts)|config/fd_cycle|jobs/2d_inplane)/')
[ -z "$CODE" ] || { echo "ERROR: code changed since the build at $BUILT:"; echo "$CODE"; exit 2; }
echo "commit $HEAD, seas_nitsche_dyn_driver built at $BUILT"

cd miniapps/seas
[ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
[ -f "$MESH" ] || { echo "ERROR: $MESH missing; vista_l2_build.sbatch makes it"; exit 2; }
SHA=$(sha256sum "$MESH" | cut -d' ' -f1)
grep -q "^$SHA  $MESH\$" elongated2p5d/MESH_SHA256_VISTA 2>/dev/null \
   || { echo "ERROR: $MESH sha256 $SHA not in elongated2p5d/MESH_SHA256_VISTA"; exit 2; }

OUT=$ROOT-runs/l2_3d_$SLURM_JOB_ID
[ ! -e "$OUT" ] || { echo "ERROR: $OUT exists"; exit 2; }
mkdir -p "$OUT" || exit 2
cp "$DECK" "$OUT/"
{ echo "commit $HEAD"; echo "built $BUILT"; echo "mesh $MESH $SHA"; echo "deck $DECK"
  echo "order $ORDER tfinal $TFINAL"
  echo "job $SLURM_JOB_ID nodes $SLURM_JOB_NUM_NODES ranks $SLURM_NTASKS"; } > "$OUT/RUN_INFO.txt"
cat "$OUT/RUN_INFO.txt"

module purge
for m in $MODULES; do module load $m || exit 2; done
srun --overlap -N "$SLURM_NNODES" -n "$SLURM_NNODES" --ntasks-per-node=1 \
   "$ROOT/miniapps/seas/jobs/fd_cycle/vista_memlog.sh" 60 > "$OUT/memlog_$SLURM_JOB_ID.txt" 2>&1 &
MEMLOG=$!
trap 'kill $MEMLOG 2>/dev/null' EXIT
s=$(date +%s)
ibrun ./seas_nitsche_dyn_driver --config "$DECK" --mesh "$MESH" --order $ORDER --tfinal $TFINAL \
   --output-dir "$OUT" --station-every 4 --dump-every 4000 --log-every 1000 \
   > "$OUT/driver.log" 2>&1
rc=$?
echo "seas_nitsche_dyn_driver exit $rc after $(( $(date +%s) - s )) s"
grep -m1 -E "^  dt = " "$OUT/driver.log"
grep "=== done" "$OUT/driver.log" || tail -20 "$OUT/driver.log"
awk '{ if ($6 > p) { p = $6; lim = $4 } } END { printf "node peak committed %.1f GiB of %.1f GiB\n", p, lim }' \
   "$OUT/memlog_$SLURM_JOB_ID.txt"
du -sh "$OUT"
exit $rc
