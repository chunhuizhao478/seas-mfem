#!/bin/bash
#SBATCH --job-name=node_fs_check
#SBATCH --clusters=condo
#SBATCH --partition=scec
#SBATCH --account=scec_608
#SBATCH --nodelist=e19-[06,08,10,12,14,16,18,20,22,24-33,39]
#SBATCH --nodes=20 --ntasks-per-node=14 --cpus-per-task=1
#SBATCH --mem=56G
#SBATCH --time=00:30:00
#SBATCH --output=node_fs_check_%j.out
#SBATCH --error=node_fs_check_%j.err
# File-read check for the condo Bus errors (Chunhui Zhao, 2026-10-16).  The debug run
# tpv104_200m_p1 (6493451) died in setup: 13 of 266 ranks could not open the 418 B station file
# that the other ranks read, then SIGBUS inside MPI_Abort on e19-20 (ld-linux paging in code).
# Earlier: 6490983 / 6491334 died on e19-20, 6491016 on e19-29 after METIS errors, TPV205 6490984
# with SIGBUS on e19-39 at t = 4.15 s.  Question: do reads of /scratch1 (and of the module
# libraries) fail on some of these nodes?  Nodes: the 19 of 6493451 plus e19-39, same layout
# (14 ranks per node, ppr:14:node).  Nothing is written outside $OUT.
#   1. on the batch node: reference sha256 of the driver, every shared library ldd resolves, the
#      dynamic loader, the TPV104 deck, station file and 200 m mesh; file system of each
#   2. every rank: sha256sum -c of that list, then 200 open + read + compare of the station file
#   3. every rank: start the driver with no arguments (MPI + hypre init, then exit 1 with the
#      usage line), which pages in the binary and all its libraries through the loader
#   4. summary per node
# Cost: 20 nodes x <= 30 min on the SCEC condo (scec_608); expected a few minutes.
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ROOT=$(cd "$SLURM_SUBMIT_DIR/../../../.." && pwd)          # repository root
cd $ROOT
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] \
   || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] \
   || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "commit $(git rev-parse HEAD)"
UNTRACKED=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.toml' '*.sbatch')
[ -z "$UNTRACKED" ] || { echo "ERROR: untracked sources in $ROOT:"; echo "$UNTRACKED"; exit 2; }
HEAD=$(git rev-parse HEAD)
cd miniapps/seas
[ -x ./seas_nitsche_dyn_driver ] || { echo "ERROR: seas_nitsche_dyn_driver not built"; exit 2; }
[ "$(cat BUILD_COMMIT_DISCOVERY 2>/dev/null)" = "$HEAD" ] \
   || { echo "ERROR: driver built at '$(cat BUILD_COMMIT_DISCOVERY 2>/dev/null)', checkout is $HEAD"; exit 2; }
module purge
module load gcc/13.3.0 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0 hdf5/1.14.3
echo "host $(hostname)  nodes $SLURM_JOB_NUM_NODES ($SLURM_JOB_NODELIST)  ranks $SLURM_NTASKS"
OUT=/scratch1/$USER/nitsche_benchmarks/node_fs_check_$SLURM_JOB_ID
mkdir -p $OUT/ranks
STATION=tpv104/stations/tpv104_fault_stations.csv
cp $STATION $OUT/station.ref

echo "=== 1. reference files ==="
{ echo $PWD/seas_nitsche_dyn_driver
  ldd ./seas_nitsche_dyn_driver | awk '$2 == "=>" && $3 ~ /^\// {print $3}'
  ldd ./seas_nitsche_dyn_driver | awk '$1 ~ /^\/.*ld-linux/ {print $1}'
  echo $PWD/tpv104/configs/tpv104_nitsche.toml
  echo $PWD/$STATION
  echo $PWD/tpv104/mesh/tpv104_200m.msh; } | xargs readlink -f | sort -u > $OUT/files.txt
xargs sha256sum < $OUT/files.txt > $OUT/reference.sha256 || { echo "ERROR: reference sha256 failed"; exit 2; }
echo "$(wc -l < $OUT/files.txt) files"
echo "file systems:"; xargs df -PT < $OUT/files.txt | awk 'NR > 1 {print "  " $2, $7}' | sort | uniq -c

cat > $OUT/rank_check.sh <<'EOF'
#!/bin/bash
r=${OMPI_COMM_WORLD_RANK:-$PMIX_RANK}; h=$(hostname -s)
sha256sum -c --quiet $OUT/reference.sha256 > $OUT/ranks/$r.sha 2>&1; sha=$?
bad=0
for i in $(seq 200); do
   cmp -s $STATION $OUT/station.ref 2>> $OUT/ranks/$r.err || bad=$((bad + 1))
done
echo "$h rank $r sha_rc $sha station_bad $bad"
EOF
export OUT STATION

echo "=== 2. every rank: sha256 of the list + 200 station-file reads ==="
s=$(date +%s)
mpirun -np $SLURM_NTASKS --map-by ppr:$SLURM_NTASKS_PER_NODE:node --bind-to core \
       -x OUT -x STATION bash $OUT/rank_check.sh > $OUT/ranks.txt 2> $OUT/ranks_mpirun.err
echo "mpirun exit $? in $(( $(date +%s) - s )) s; $(wc -l < $OUT/ranks.txt) of $SLURM_NTASKS ranks reported"

echo "=== 3. every rank: start the driver with no arguments (expect exit 1, usage once) ==="
s=$(date +%s)
mpirun -np $SLURM_NTASKS --map-by ppr:$SLURM_NTASKS_PER_NODE:node --bind-to core \
       ./seas_nitsche_dyn_driver > $OUT/driver_noargs.log 2>&1
echo "mpirun exit $? in $(( $(date +%s) - s )) s; usage lines $(grep -c '^usage:' $OUT/driver_noargs.log)"
grep -iE 'signal|bus error|segmentation|abort' $OUT/driver_noargs.log | head -20

echo "=== 4. summary per node (ranks with a sha mismatch / unreadable file, or a failed station read) ==="
awk '{n[$1]++; if ($5 != 0) s[$1]++; if ($7 != 0) b[$1]++}
     END {for (h in n) printf "%s ranks %d sha_bad %d station_bad %d\n", h, n[h], s[h] + 0, b[h] + 0}' \
    $OUT/ranks.txt | sort
echo "first sha errors:"; cat $OUT/ranks/*.sha 2>/dev/null | head -20
echo "first read errors:"; cat $OUT/ranks/*.err 2>/dev/null | head -20
echo "mpirun stderr:"; head -20 $OUT/ranks_mpirun.err
