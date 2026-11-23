# shellcheck shell=bash
# d1_env.sh: shared pieces of the GPU dynamic-rupture plan's D1 jobs on Vista
# (d1_gg_build, d1_gg_check, d1_cuda_build, d1_gh_check).  Sourced, not a job.
#
# The four D1 cases are TPV205 and TPV104 on the 1000 m meshes at p = 1 and
# p = 2 with the committed decks (D0a's cases), run to the decks' tfinal of
# 12 s (Vista gates run the full 12 s, Chunhui Zhao 2026-10-19).

# shellcheck disable=SC2034   # read by the jobs that source this file
D1_CASES=("tpv205 tpv205/configs/tpv205_nitsche.toml tpv205/mesh/tpv2053d_1000m.msh"
          "tpv104 tpv104/configs/tpv104_nitsche.toml tpv104/mesh/tpv104_1000m.msh")
D1_ORDERS=(1 2)
# The meshes D0a's job 1080349 generated with TACC gmsh/4.15.0; every D1 job
# checks it reads the same files.
D1_MESH_SHA256="c6eeaaa5350d06e6801f6563e9c5b47b7486dd5a647ee7f55237b0b664c04c3e  tpv205/mesh/tpv2053d_1000m.msh
f314e3af62591f8f81a38cf226305adb1deb020542b5b60f18660685ab490c08  tpv104/mesh/tpv104_1000m.msh"

# Lmod in a batch shell
# shellcheck disable=SC1091   # Vista's Lmod init files
if ! type module >/dev/null 2>&1; then
   [ -f /etc/profile.d/z00_lmod.sh ] && source /etc/profile.d/z00_lmod.sh
   [ -f /usr/share/lmod/lmod/init/bash ] && source /usr/share/lmod/lmod/init/bash
fi

# The repo CLAUDE.md guard: fetched, no tracked changes, HEAD on GitHub, no
# untracked sources; prints the commit.
d1_guard () {
   local d=$1 i U
   cd "$d" || exit 2
   for i in 1 2 3 4; do
      git fetch -q origin && break
      [ "$i" = 4 ] && { echo "ERROR: git fetch failed in $d"; exit 2; }
      sleep $(( 5 + RANDOM % 20 ))
   done
   [ -z "$(git status --porcelain --untracked-files=no)" ] \
      || { echo "ERROR: uncommitted changes in $d"; git status --short; exit 2; }
   [ -n "$(git branch -r --contains HEAD)" ] \
      || { echo "ERROR: HEAD $(git rev-parse --short HEAD) in $d is not on GitHub"; exit 2; }
   U=$(git ls-files --others --exclude-standard -- '*.cpp' '*.hpp' '*.inl' Makefile '*.toml' '*.sbatch')
   [ -z "$U" ] || { echo "ERROR: untracked sources in $d:"; echo "$U"; exit 2; }
   echo "commit $(git rev-parse HEAD) in $d"
}

# The binary in clone $1 was built at its HEAD: BUILD_FILE ($2) holds HEAD.
d1_check_build () {
   local d=$1 f=$2
   [ -x "$d/miniapps/seas/seas_nitsche_dyn_driver" ] \
      || { echo "ERROR: no seas_nitsche_dyn_driver in $d; run the build job first"; exit 2; }
   [ "$(cat "$d/miniapps/seas/$f" 2>/dev/null)" = "$(git -C "$d" rev-parse HEAD)" ] \
      || { echo "ERROR: $d/miniapps/seas/$f does not hold HEAD $(git -C "$d" rev-parse HEAD); rebuild"; exit 2; }
   echo "binary in $d built at $(cat "$d/miniapps/seas/$f")"
}

# The meshes in miniapps/seas ($1) are D0a's; prints their sha256.
d1_check_meshes () {
   local s
   s=$(cd "$1" && sha256sum tpv205/mesh/tpv2053d_1000m.msh tpv104/mesh/tpv104_1000m.msh)
   echo "$s"
   [ "$s" = "$D1_MESH_SHA256" ] || { echo "MESHES: FAIL (not D0a's sha256)"; return 1; }
   echo "MESHES: same as D0a 1080349"
}

# The modules of build_vista.sh (gg runs) plus gmsh.
d1_gg_modules () {
   module purge
   module load TACC
   module swap nvidia gcc/14.2.0 2>/dev/null || module load gcc/14.2.0
   module load cuda/12.6 openmpi/5.0.5 hypre/2.31.0 metis/5.1.0.3 phdf5/1.14.6 gmsh/4.15.0
}

# The modules of build_horizon.sh MODE=cuda (gh runs), or G0_MODULES when set
# (jobs/horizon/h_env.sh: Horizon's GB200 nodes).
d1_gh_modules () {
   module purge
   # shellcheck disable=SC2086   # a list of modules
   module load ${G0_MODULES:-gcc/13.2.0 ${CUDA_MODULE:-cuda/13.0} openmpi/5.0.11 nvpl/26.9 cmake/3.31.5}
}

# Bitwise check of run dir $1 against run dir $2: same exit code (0), same
# output file set, and every station file, station-copy file, station snap
# list, rupture-time and initial-state file byte-identical (fault.vtkhdf,
# logs and gate files excluded: HDF5 stores object times, logs hold timings).
# Prints one line ("yes (n files)" or "NO" with what differs); returns 1 on a
# difference.
d1_bitwise () {
   local N=$1 R=$2 ok=1 n=0 f why=""
   [ "$(cat "$R/exit_code" 2>/dev/null)" = 0 ] || { why="$why; ref exit $(cat "$R/exit_code" 2>/dev/null)"; ok=0; }
   cmp -s "$N/exit_code" "$R/exit_code" || { why="$why; exit codes differ"; ok=0; }
   find "$N" -maxdepth 1 -type f -printf '%f\n' | grep -v -E '\.vtkhdf$|\.log$|_gate\.json$|^exit_code$' | sort > "$N.files"
   find "$R" -maxdepth 1 -type f -printf '%f\n' | grep -v -E '\.vtkhdf$|\.log$|_gate\.json$|^exit_code$' | sort > "$R.files"
   cmp -s "$N.files" "$R.files" || { why="$why; file sets differ"; ok=0; }
   while read -r f; do
      [ -f "$R/$f" ] || continue
      n=$((n + 1))
      cmp -s "$N/$f" "$R/$f" || { why="$why; differs: $f"; ok=0; }
   done < "$N.files"
   [ "$n" -gt 0 ] || { why="$why; no files"; ok=0; }
   if [ "$ok" = 1 ]; then echo "yes ($n files byte-identical)"; else echo "NO ($n files compared$why)"; return 1; fi
}
