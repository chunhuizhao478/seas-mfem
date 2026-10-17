#!/bin/bash
#SBATCH -J nitsche_mirror_meshes
#SBATCH -o nitsche_mirror_meshes_%j.out
#SBATCH -e nitsche_mirror_meshes_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:30:00
#SBATCH -A EAR26024
# Mirror-symmetric 200 m and 100 m meshes for TPV104 (22 km fault) and TPV205 on Vista (Chunhui Zhao,
# 2026-10-17): scripts/make_mirror_meshes.sh with TACC gmsh, sha256 to miniapps/seas/MESH_SHA256_MIRROR_VISTA,
# which run_nitsche_vista.sh accepts like MESH_SHA256_VISTA.  Submit from miniapps/seas/jobs/vista_nitsche/.
# Cloud (gmsh 4.15.2, one core): 200 m 3,453,310 tets (TPV104) and 1,976,500 (TPV205) in about 2 and
# 1 min; 100 m 11,473,766 and 6,036,560 tets in 10.5 and 4.5 min; meshes 1.9 GB in all.
# Cost: 1 node x <= 1.5 h x 0.33 SU = <= 0.5 SU.
export LC_ALL=C LANG=C
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd $ROOT
for i in 1 2 3 4; do                                   # retry as run_nitsche_vista.sh (ref lock)
   git fetch -q origin && break
   [ $i = 4 ] && { echo "ERROR: git fetch failed"; exit 2; }
   sleep $(( 5 + RANDOM % 20 ))
done
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git rev-parse HEAD)"
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
# same modules as build_verify_vista.sbatch (gmsh is the one used here)
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6 gmsh/4.15.0
cd miniapps/seas
T0=$(date +%s)
SHA_FILE=MESH_SHA256_MIRROR_VISTA bash scripts/make_mirror_meshes.sh 200 100
rc=$?
echo "make_mirror_meshes.sh exit $rc after $(( $(date +%s) - T0 )) s"
exit $rc
