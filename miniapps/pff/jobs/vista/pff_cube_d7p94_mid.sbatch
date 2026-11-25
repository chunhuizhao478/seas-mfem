#!/bin/bash
# pff 3D SPE cube cube_d7p94_mid_eta0p45_absorb, 20 pulses (1.2 ms), on Vista Grace nodes.
# pff counterpart of the MOOSE Vista run $SCRATCH/spe3d_vista/cube_d7p94_mid_eta0p45_absorb
# (farms@cdms 46ad2eb code_3d/elasticity_absorbing.i), on the same mesh: MOOSE's
# spe_cube_d7p94_h0p3.msh (md5 9109e05c) converted to msh 2.2 (md5 482e9a25).
# Steps: commit guard -> pff_partition (serial on the first node, ~6 GB, ~3 min) -> pff_driver -pm.
# Submit from Vista: mkdir -p $SCRATCH/pff/logs && sbatch $WORK/pff_build/seas-mfem/miniapps/pff/jobs/vista/pff_cube_d7p94_mid.sbatch
# Build: $WORK/pff_build (BUILD_NOTES.txt). Results: $SCRATCH/pff/<run>_<jobid> only.
#SBATCH -J pff_cube3d
#SBATCH -p gg
#SBATCH -N 4
#SBATCH -n 576
#SBATCH -t 12:00:00
#SBATCH -A EAR20006
#SBATCH -o /scratch/10024/zhaochun/pff/logs/%x_%j.out
#SBATCH -e /scratch/10024/zhaochun/pff/logs/%x_%j.err
set -eo pipefail   # -u after the module commands (Lmod reads unset variables)

SRC=$WORK/pff_build/seas-mfem
PFF=$SRC/miniapps/pff
MESH=$SCRATCH/pff/mesh/spe_cube_d7p94_h0p3_moose_msh22.msh
MESH_MD5=482e9a2520cf39c02bf24fdb35b8ae92
NAME=cube_d7p94_mid_eta0p45_absorb
RUN=$SCRATCH/pff/${NAME}_${SLURM_JOB_ID}

module reset
module load gcc/14.2.0 openmpi/5.0.11 hypre/2.31.0 metis/5.1.0.3
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

mkdir -p "$RUN"
cd "$RUN"
echo "$head" > COMMIT
cp "$0" job.sbatch 2>/dev/null || cp "$PFF/jobs/vista/pff_cube_d7p94_mid.sbatch" job.sbatch
echo "$MESH_MD5  $MESH" | md5sum -c -

# Run config = the committed case file with [mesh] file pointed at the MOOSE mesh (only change).
sed "s#^file = \"../mesh/spe_cube_d7p94_h0p3.msh\"#file = \"$MESH\"#" \
    "$PFF/spe_hftc/configs/$NAME.toml" > "$NAME.toml"
diff "$PFF/spe_hftc/configs/$NAME.toml" "$NAME.toml" || true
grep -q "^file = \"$MESH\"" "$NAME.toml"

# ---- partition once (serial), then run ----------------------------------------
"$PFF/pff_partition" -m "$MESH" -np "$SLURM_NTASKS" -o "$RUN/parts_np$SLURM_NTASKS"
ibrun "$PFF/pff_driver" -c "$RUN/$NAME.toml" -o "$RUN/out" -pm "$RUN/parts_np$SLURM_NTASKS" -pv 3300
echo "finished $(date)"
