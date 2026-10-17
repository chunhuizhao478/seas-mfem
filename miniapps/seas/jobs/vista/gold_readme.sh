#!/bin/bash
# gold_readme.sh — README.md for a results folder on Vista: the GitHub commit
# the results came from, the Slurm job, the job script, the flags, the
# meshes with their sha256, the binary and its build, and the commands that
# reproduce them (project rule, Chunhui Zhao, Oct 17: every gold result
# names the commit it was made with).
#
# Two uses:
#   - sourced by g0_env.sh: g0_provenance DIR calls gold_readme_job DIR, so
#     every job that records provenance also writes DIR/README.md;
#   - run on a login node for folders made before this file existed:
#       bash jobs/vista/gold_readme.sh DIR [DIR ...]
#     Each DIR needs the provenance.txt its job wrote.  The commit and the
#     job come from it, the job's details from sacct, and the meshes and
#     flags from that commit (git show), so the README describes what ran,
#     not today's checkout.  An existing README.md is left alone.

GOLD_README_REPO="https://github.com/chunhuizhao478/seas-mfem"

# Meshes the G-phase jobs use (relative to miniapps/seas), as in g0_env.sh.
GOLD_README_MESHES=(
   bp5/mesh/reference/bp5_tandem_exact.msh
   bp5/mesh/reference/bp5_tandem_500m.msh
   bp5/mesh/reference/bp5_tandem_250m.msh
)

# Settings a job may take from sbatch --export.
GOLD_README_KNOBS=(ORDER RES NP TFINAL_YR ALLOW_P2_MULTIRANK CUDA_MODULE
                   G1_GOLDEN_ROOT)

# _gold_readme_emit: print the README from the GR_* variables below.
#   GR_TITLE GR_HASH GR_BRANCH GR_DIRTY GR_JOB GR_JOBNAME GR_HOST GR_START
#   GR_PARTITION GR_ACCOUNT GR_NODES GR_TASKS GR_SCRIPT GR_SUBMIT GR_KNOBS
#   GR_MODE GR_BINARY GR_FLAGS GR_MESHES GR_BUILD GR_NOTE
_gold_readme_emit() {
   local fence='```'
   echo "# ${GR_TITLE}"
   echo
   echo "Produced from commit \`${GR_HASH}\` (branch \`${GR_BRANCH}\`) of ${GOLD_README_REPO}:"
   echo "${GOLD_README_REPO}/commit/${GR_HASH}"
   echo
   if [ "${GR_DIRTY}" = "0" ]; then
      echo "The checkout had no uncommitted changes in tracked files, so these results come from that commit exactly."
   else
      echo "WARNING: the checkout had ${GR_DIRTY} tracked file(s) with uncommitted changes; these results do not come from the commit alone."
   fi
   if [ -n "${GR_NOTE}" ]; then
      echo
      echo "${GR_NOTE}"
   fi
   echo
   echo "## Job"
   echo
   echo "- Cluster: Vista (${GR_HOST}), Slurm job ${GR_JOB} (${GR_JOBNAME}), started ${GR_START}"
   echo "- Partition ${GR_PARTITION}, account ${GR_ACCOUNT}, ${GR_NODES} node(s), ${GR_TASKS} task(s)"
   echo "- Job script: \`${GR_SCRIPT}\` at that commit"
   echo "- Submit line: \`${GR_SUBMIT}\`"
   echo "- Settings from the environment: ${GR_KNOBS:-none}"
   echo
   echo "## Inputs"
   echo
   echo "- Flags shared by the G-phase runs (\`G0_V63_FLAGS\` in \`miniapps/seas/jobs/vista/g0_env.sh\`; the job script adds the mesh, order and per-run flags):"
   echo
   echo "${fence}"
   echo "${GR_FLAGS}"
   echo "${fence}"
   echo
   echo "- Meshes at that commit (paths relative to \`miniapps/seas\`), sha256:"
   echo
   echo "${fence}"
   echo "${GR_MESHES}"
   echo "${fence}"
   echo
   echo "## Build"
   echo
   echo "- Binary: ${GR_BINARY}"
   echo "- Built with \`MODE=${GR_MODE} bash build_horizon.sh\` (PETSc and its packages, MFEM, then seas). What the build used:"
   echo
   echo "${fence}"
   echo "${GR_BUILD}"
   echo "${fence}"
   echo
   echo "## Reproduce"
   echo
   echo "${fence}"
   echo "git clone ${GOLD_README_REPO}.git && cd seas-mfem"
   echo "git checkout ${GR_HASH}"
   if [ "${GR_MODE}" = "cuda" ]; then
      echo "MODE=cuda GPU_ARCH=90 STAGE=fetch bash build_horizon.sh   # login node"
      echo "cd miniapps/seas && sbatch jobs/vista/g1_cuda_build.sbatch  # deps, MFEM, seas"
   else
      echo "STAGE=fetch bash build_horizon.sh                          # login node"
      echo "cd miniapps/seas && sbatch jobs/vista/g0_build.sbatch       # deps, MFEM, seas"
   fi
   echo "${GR_SUBMIT}   # from miniapps/seas, after the build"
   echo "${fence}"
}

# _gold_readme_meshes ROOT [HASH]: "sha256  path" for each mesh, from the
# checkout ROOT, or from commit HASH of it when given.
_gold_readme_meshes() {
   local root="$1" hash="${2:-}" m sum
   for m in "${GOLD_README_MESHES[@]}"; do
      if [ -n "${hash}" ]; then
         if git -C "${root}" cat-file -e "${hash}:miniapps/seas/${m}" 2>/dev/null; then
            sum="$(git -C "${root}" show "${hash}:miniapps/seas/${m}" | sha256sum | cut -d' ' -f1)"
         else
            sum="(not in this commit)"
         fi
      elif [ -f "${root}/miniapps/seas/${m}" ]; then
         sum="$(sha256sum "${root}/miniapps/seas/${m}" | cut -d' ' -f1)"
      else
         sum="(missing)"
      fi
      echo "${sum}  ${m}"
   done
}

# _gold_readme_flags FILE_TEXT: the G0_V63_FLAGS array body of g0_env.sh.
_gold_readme_flags() {
   sed -n '/^G0_V63_FLAGS=(/,/^)/p' | sed '1d;$d' | sed 's/^ *//'
}

# _gold_readme_mode ROOT: cuda or cpu, from MFEM's config.
_gold_readme_mode() {
   if grep -q '^MFEM_USE_CUDA *= YES' "$1/config/config.mk" 2>/dev/null; then
      echo cuda
   else
      echo cpu
   fi
}

# gold_readme_job DIR: README.md for DIR from inside a running job (called
# by g0_provenance).  Needs SEAS_MFEM_ROOT and SEAS_DIR from g0_env.sh.
gold_readme_job() {
   local dir="$1" k
   local root="${SEAS_MFEM_ROOT}"
   GR_TITLE="$(basename "$(cd "${dir}" && pwd)")"
   GR_HASH="$(git -C "${root}" rev-parse HEAD)"
   GR_BRANCH="$(git -C "${root}" rev-parse --abbrev-ref HEAD)"
   GR_DIRTY="$(git -C "${root}" status --porcelain --untracked-files=no | wc -l)"
   GR_JOB="${SLURM_JOB_ID:-none}"
   GR_JOBNAME="${SLURM_JOB_NAME:-}"
   GR_HOST="$(hostname)"
   GR_START="$(date -u +%FT%TZ)"
   GR_PARTITION="${SLURM_JOB_PARTITION:-?}"
   GR_ACCOUNT="${SLURM_JOB_ACCOUNT:-?}"
   GR_NODES="${SLURM_JOB_NUM_NODES:-?}"
   GR_TASKS="${SLURM_NTASKS:-?}"
   GR_SCRIPT="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/^ *Command=//p')"
   GR_SCRIPT="${GR_SCRIPT#"${root}"/}"
   GR_SUBMIT="$(sacct -j "${SLURM_JOB_ID:-0}" -X -n -P -o SubmitLine 2>/dev/null | head -n 1)"
   [ -n "${GR_SUBMIT}" ] || GR_SUBMIT="sbatch ${GR_SCRIPT#miniapps/seas/}"
   GR_KNOBS=""
   for k in "${GOLD_README_KNOBS[@]}"; do
      [ -n "${!k:-}" ] && GR_KNOBS="${GR_KNOBS}${GR_KNOBS:+, }${k}=${!k}"
   done
   GR_MODE="$(_gold_readme_mode "${root}")"
   GR_BINARY="\`miniapps/seas/seas_bp5_full\`, sha256 $(sha256sum "${SEAS_DIR}/seas_bp5_full" 2>/dev/null | cut -d' ' -f1)"
   GR_FLAGS="$(printf '%s\n' "${G0_V63_FLAGS[@]}" | paste -sd' ')"
   GR_MESHES="$(_gold_readme_meshes "${root}")"
   GR_BUILD="$(cat "${G0_DEPS_PROVENANCE:-${WORK:-$HOME}/seas-gpu-deps/build_provenance_${GR_MODE}.txt}" 2>/dev/null || echo "(no build record found)")"
   GR_NOTE=""
   _gold_readme_emit > "${dir}/README.md"
}

# gold_readme_backfill DIR: README.md for a folder a job wrote before this
# file existed, from its provenance.txt.
gold_readme_backfill() {
   local dir="$1" root prov line
   prov="${dir}/provenance.txt"
   if [ ! -f "${prov}" ]; then
      echo "SKIP ${dir}: no provenance.txt" >&2
      return 1
   fi
   if [ -e "${dir}/README.md" ]; then
      echo "SKIP ${dir}: README.md exists" >&2
      return 0
   fi
   root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
   # provenance.txt lines (g0_provenance):
   #   job: ID NAME on HOST, DATE
   #   nodes: N, tasks: T, partition: P
   #   source: HASH (BRANCH), dirty files: D
   #   binary: DATE SIZE sha256 PREFIX
   line="$(sed -n 's/^job: //p' "${prov}")"
   GR_JOB="$(echo "${line}" | awk '{print $1}')"
   GR_JOBNAME="$(echo "${line}" | awk '{print $2}')"
   GR_HOST="$(echo "${line}" | sed -n 's/.* on \([^,]*\),.*/\1/p')"
   GR_START="$(echo "${line}" | sed -n 's/.*, //p')"
   line="$(sed -n 's/^nodes: //p' "${prov}")"
   GR_NODES="$(echo "${line}" | sed -n 's/^\([^,]*\),.*/\1/p')"
   GR_TASKS="$(echo "${line}" | sed -n 's/.*tasks: \([^,]*\),.*/\1/p')"
   GR_PARTITION="$(echo "${line}" | sed -n 's/.*partition: //p')"
   line="$(sed -n 's/^source: //p' "${prov}")"
   GR_HASH="$(echo "${line}" | awk '{print $1}')"
   GR_BRANCH="$(echo "${line}" | sed -n 's/^[^ ]* (\([^)]*\)).*/\1/p')"
   GR_DIRTY="$(echo "${line}" | sed -n 's/.*dirty files: *//p')"
   if ! git -C "${root}" cat-file -e "${GR_HASH}^{commit}" 2>/dev/null; then
      echo "SKIP ${dir}: commit ${GR_HASH} is not in ${root}; git fetch first" >&2
      return 1
   fi
   GR_TITLE="$(basename "$(cd "${dir}" && pwd)")"
   GR_ACCOUNT="$(sacct -j "${GR_JOB}" -X -n -P -o Account 2>/dev/null | head -n 1)"
   GR_SUBMIT="$(sacct -j "${GR_JOB}" -X -n -P -o SubmitLine 2>/dev/null | head -n 1)"
   GR_SCRIPT="$(echo "${GR_SUBMIT}" | grep -o 'jobs/[^ ]*\.sbatch' | head -n 1)"
   [ -n "${GR_SCRIPT}" ] && GR_SCRIPT="miniapps/seas/${GR_SCRIPT}" || GR_SCRIPT="(see the submit line)"
   [ -n "${GR_ACCOUNT}" ] || GR_ACCOUNT="?"
   [ -n "${GR_SUBMIT}" ] || GR_SUBMIT="(sacct has no submit line for job ${GR_JOB})"
   GR_KNOBS="(not recorded by this job; see the submit line)"
   GR_MODE="cpu"
   grep -q 'MFEM_USE_CUDA *= YES' "${prov}" && GR_MODE="cuda"
   GR_BINARY="\`miniapps/seas/seas_bp5_full\`, $(sed -n 's/^binary: [^ ]* [^ ]* //p' "${prov}") (first 16 hex digits, from provenance.txt)"
   GR_FLAGS="$(git -C "${root}" show "${GR_HASH}:miniapps/seas/jobs/vista/g0_env.sh" 2>/dev/null | _gold_readme_flags | paste -sd' ')"
   GR_MESHES="$(_gold_readme_meshes "${root}" "${GR_HASH}")"
   GR_BUILD="$(sed -n '/^--- build provenance ---$/,$p' "${prov}" | sed '1d')"
   [ -n "${GR_BUILD}" ] || GR_BUILD="(not in provenance.txt)"
   GR_NOTE="This README was written after the run, on $(date -u +%F), from the job's provenance.txt, sacct and the commit; the run itself recorded only provenance.txt."
   _gold_readme_emit > "${dir}/README.md"
   echo "WROTE ${dir}/README.md (commit ${GR_HASH}, job ${GR_JOB})"
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
   if [ $# -eq 0 ]; then
      echo "usage: bash jobs/vista/gold_readme.sh DIR [DIR ...]" >&2
      exit 2
   fi
   rc=0
   for d in "$@"; do
      gold_readme_backfill "${d}" || rc=1
   done
   exit "${rc}"
fi
