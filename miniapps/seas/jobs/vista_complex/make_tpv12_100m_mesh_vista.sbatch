#!/bin/bash
#SBATCH -J tpv12_100m_mesh
#SBATCH -o tpv12_100m_mesh_%j.out
#SBATCH -e tpv12_100m_mesh_%j.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -A EAR26024
# The 100 m TPV12 mesh (plan Q2: p = 1 at 100 m for TPV12, whose cohesive zone is 210 m at its base) on
# Vista (Chunhui Zhao, 2026-10-21): gmsh on tpv12/mesh/tpv12_100m.geo, msh 2.2, sha256 to
# miniapps/seas/MESH_SHA256_TPV12_100M_VISTA, which run_complex_vista.sh checks (a file of its own, so a
# later build-verify, which rewrites MESH_SHA256_COMPLEX_VISTA, keeps it).  Submit from
# miniapps/seas/jobs/vista_complex/.
# The cloud could not hold this mesh.  Cloud gmsh times for the 200 m TPV12 mesh (1.88 M tets): 3 min; the
# TPV205 100 m mesh has 2.9 times the tets of its 200 m one, so about 5.4 M tets and 10-20 min here (inferred).
# Cost: 1 node x <= 2 h x 0.33 SU = <= 0.66 SU.
export LC_ALL=C LANG=C
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd $ROOT
for i in 1 2 3 4; do                                   # retry as run_complex_vista.sh (ref lock)
   git fetch -q origin && break
   [ $i = 4 ] && { echo "ERROR: git fetch failed"; exit 2; }
   sleep $(( 5 + RANDOM % 20 ))
done
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git rev-parse HEAD)"
UNTRACKED=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.toml' '*.sbatch' '*.geo' '*.csv')
[ -z "$UNTRACKED" ] || { echo "ERROR: untracked sources in $ROOT:"; echo "$UNTRACKED"; exit 2; }
module purge; module load TACC
module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6 gmsh/4.15.0
cd miniapps/seas
G=tpv12/mesh/tpv12_100m
T0=$(date +%s)
gmsh -3 $G.geo -format msh22 -o $G.msh > tpv12_100m_gmsh_$SLURM_JOB_ID.log 2>&1
rc=$?
ne=$(grep -a -m1 -A1 '^\$Elements' $G.msh 2>/dev/null | tail -1)
echo "gmsh exit $rc after $(( $(date +%s) - T0 )) s; ${ne:-0} elements (all types)"
[ $rc -eq 0 ] && [ "${ne:-0}" -gt 0 ] || { echo "ERROR: $G.msh not made"; grep -a "Error" tpv12_100m_gmsh_$SLURM_JOB_ID.log | head -5; exit 3; }
sha256sum $G.msh | tee MESH_SHA256_TPV12_100M_VISTA
