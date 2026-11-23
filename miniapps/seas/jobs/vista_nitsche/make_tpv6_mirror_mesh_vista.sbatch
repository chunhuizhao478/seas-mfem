#!/bin/bash
#SBATCH -J tpv6_mirror_mesh
#SBATCH -o tpv6_mirror_mesh_%j.out
#SBATCH -e tpv6_mirror_mesh_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:30:00
#SBATCH -A EAR26024
# Mirror-symmetric 100 m TPV6/TPV7 mesh on Vista (Chunhui Zhao, 2026-11-22): scripts/make_mirror_meshes.sh
# case tpv6 (tpv6/mesh/tpv6_nitsche_half.geo, built-in kernel) with TACC gmsh, sha256 to
# miniapps/seas/MESH_SHA256_TPV6_MIRROR_VISTA, which run_nitsche_vista.sh accepts like MESH_SHA256_VISTA.
# Its own sha file, so it does not overwrite the TPV104 one of make_mirror_meshes_vista.sbatch.
# Why: TPV7 100 m p = 3 (1159765) aborted at t = 5.99 s, sigma_n_eff -55 kPa at the fault end at the free
# surface, on the gmsh mesh tpv6_nitsche_100m.msh (two halves meshed separately).
# Size (inferred from tpv6_nitsche_100m.msh, 6,199,161 tets): about 3.1 M tets per half, 6.2 M mirrored.
# Submit from miniapps/seas/jobs/vista_nitsche/.  Cost: 1 node x <= 1.5 h x 0.33 SU = <= 0.5 SU.
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
module load openmpi/5.0.11 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6 gmsh/4.15.0
cd miniapps/seas
T0=$(date +%s)
MIRROR_CASES=tpv6 SHA_FILE=MESH_SHA256_TPV6_MIRROR_VISTA bash scripts/make_mirror_meshes.sh 100
rc=$?
echo "make_mirror_meshes.sh exit $rc after $(( $(date +%s) - T0 )) s"
exit $rc
