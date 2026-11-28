#!/bin/bash
# pff cmame 3D hydromech em_5J (cmame_revision/configs/3d_hydromech_em_5J.toml), first pulses, on Vista Grace nodes.
# MOOSE input: farms@cdms pulsepower/cmame_revision/3d_hydromech/em_5J/elasticity_mesh2x.i (no MOOSE output of it
# is reachable, 2026-11-27: this run is the pff timing and a physical check, not yet a comparison).
# Steps: commit guard -> mesh (gmsh module, serial on the first node, once; kept in $SCRATCH/pff/mesh)
#        -> pff_partition (serial) -> pff_driver -pm.
# TF = end time [s] (default 5e-5 = 5 pulses of 10 us; the case file has 50 pulses).
# Submit from Vista: mkdir -p $SCRATCH/pff/logs && sbatch $WORK/pff_build/seas-mfem/miniapps/pff/jobs/vista/pff_cmame3d_hydromech_em_5J.sbatch
# Build: $WORK/pff_build (BUILD_NOTES.txt). Results: $SCRATCH/pff/<run>_<jobid> only.
#SBATCH -J pff_cmame3d_hm
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 06:00:00
#SBATCH -A EAR20006
#SBATCH -o /scratch/10024/zhaochun/pff/logs/%x_%j.out
#SBATCH -e /scratch/10024/zhaochun/pff/logs/%x_%j.err
set -eo pipefail   # -u after the module commands (Lmod reads unset variables)

SRC=$WORK/pff_build/seas-mfem
PFF=$SRC/miniapps/pff
NAME=3d_hydromech_em_5J
GEO=$PFF/cmame_revision/3d_mesh/cylinder_sample.geo
MESH=$SCRATCH/pff/mesh/cmame_cylinder_sample_gmsh4.15.0_msh22.msh
TF=${TF:-5e-5}
RUN=$SCRATCH/pff/cmame_${NAME}_${SLURM_JOB_ID}

module reset
module load gcc/14.2.0 openmpi/5.0.11 hypre/2.31.0 metis/5.1.0.3 gmsh/4.15.0
module list 2>&1
# OpenMPI 5.0.11 coll/han deadlocked in MPI_Allreduce in the MOOSE 3D chain (Vista 1164144);
# pff reduces every step (PCG), so keep HAN and hcoll off for this job only.
export OMPI_MCA_coll=^hcoll,han
set -u

# ---- commit guard: run only a pushed, clean checkout --------------------------
head=$(git -C "$SRC" rev-parse HEAD)
[ -z "$(git -C "$SRC" status --porcelain --untracked-files=no)" ] || { echo "commit guard: tracked changes in $SRC"; git -C "$SRC" status --short; exit 1; }
git -C "$SRC" branch -r --contains "$head" | grep -q "origin/cdms" || { echo "commit guard: $head not on origin/cdms"; exit 1; }
echo "commit $head (clean, on origin/cdms)"

mkdir -p "$RUN" "$(dirname "$MESH")"
cd "$RUN"
echo "$head" > COMMIT
cp "$0" job.sbatch 2>/dev/null || cp "$PFF/jobs/vista/pff_cmame3d_hydromech_em_5J.sbatch" job.sbatch

# ---- mesh: the farms .geo unchanged (lc = 0.2 mm), Gmsh format 2.2 for MFEM ----
if [ ! -s "$MESH" ]; then
   gmsh --version
   gmsh "$GEO" -3 -format msh22 -o "$MESH.tmp_$SLURM_JOB_ID.msh" 2>&1 | tail -20
   mv "$MESH.tmp_$SLURM_JOB_ID.msh" "$MESH"
fi
md5sum "$MESH" | tee MESH_MD5
ls -l "$MESH"

# Run config = the committed case file with [mesh] file pointed at the mesh (only change).
sed "s#^file = \"../3d_mesh/cylinder_sample.msh\"#file = \"$MESH\"#" \
    "$PFF/cmame_revision/configs/$NAME.toml" > "$NAME.toml"
diff "$PFF/cmame_revision/configs/$NAME.toml" "$NAME.toml" || true
grep -q "^file = \"$MESH\"" "$NAME.toml"

# ---- partition once (serial), then run ----------------------------------------
"$PFF/pff_partition" -m "$MESH" -np "$SLURM_NTASKS" -o "$RUN/parts_np$SLURM_NTASKS"
ibrun "$PFF/pff_driver" -c "$RUN/$NAME.toml" -o "$RUN/out" -pm "$RUN/parts_np$SLURM_NTASKS" -tf "$TF" -pv 2400
echo "finished $(date)"
