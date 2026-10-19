#!/bin/bash
# GPU plan G0 on Vista: print (default) or submit every G0 job.
#
#   cd miniapps/seas
#   bash jobs/vista/g0_submit.sh                    # dry run: print the sbatch lines
#   SUBMIT=1 bash jobs/vista/g0_submit.sh           # submit them
#   SUBMIT=1 STEPS="goldens timing" BUILD_JOB=none bash jobs/vista/g0_submit.sh
#   SUBMIT=1 STEPS="counts factor" MESHES="500m 250m" BUILD_JOB=none \
#      bash jobs/vista/g0_submit.sh
#
# STEPS (default "build goldens timing counts factor accst"):
#   build    jobs/vista/g0_build.sbatch (deps, MFEM, driver, unit test)
#   goldens  50-step goldens and floors: p=1 on 1, 2, 4, 8 ranks; p=2 on 1
#   timing   1,000 steps on 1 node / 144 ranks at p=1 and p=2
#   counts   --exit-after-assembly on prod, 500 m, 250 m at p=1 and p=2
#   factor   factor size: prod at p=1 and p=2 on the timing layout, then
#            500m_p1 -> 250m_p1 -> 500m_p2 -> 250m_p2, each after the last
#            succeeded (a failure cancels the rest)
#   accst    accepted-state check, 50 fresh steps, p=1 on 8 ranks
# MESHES (default "prod 500m 250m") limits the counts and factor jobs to
# those meshes.
# Every job after the build waits for it (afterok) unless BUILD_JOB=none
# (build already done) or BUILD_JOB=<jobid>.
#
# Node counts and wall times for the factor-size and p=2 jobs are estimates
# (scaled from the production mesh's ~18 GB MUMPS workspace on one rank, not
# measured on Vista); the jobs report MUMPS's own estimate if they stop.

set -euo pipefail
SUBMIT="${SUBMIT:-0}"
STEPS="${STEPS:-build goldens timing counts factor accst}"
BUILD_JOB="${BUILD_JOB:-}"
MESHES="${MESHES:-prod 500m 250m}"

if [ ! -f jobs/vista/g0_env.sh ]; then
   echo "ERROR: run from miniapps/seas" >&2
   exit 1
fi

n=0
JOB=""
# sub NAME [sbatch args...] : prints the sbatch line, or submits it; the job
# id (a placeholder in a dry run) is left in JOB.
sub() {
   local name="$1"; shift
   n=$((n + 1))
   if [ "${SUBMIT}" = 1 ]; then
      # Vista's sbatch prints a banner on stdout even with --parsable; the
      # job id ("id" or "id;cluster") is the last line.
      local out
      out="$(sbatch --parsable "$@")"
      JOB="$(printf '%s\n' "${out}" | tail -n 1)"
      JOB="${JOB%%;*}"
      if ! [[ "${JOB}" =~ ^[0-9]+$ ]]; then
         echo "ERROR: ${name}: no job id in sbatch output:" >&2
         printf '%s\n' "${out}" >&2
         exit 1
      fi
      echo "${name}: job ${JOB}"
   else
      JOB="JOB${n}"
      printf '%s (%s): sbatch' "${name}" "${JOB}"
      printf ' %q' "$@"
      printf '\n'
   fi
}

has() { [[ " ${STEPS} " == *" $1 "* ]]; }
has_mesh() { [[ " ${MESHES} " == *" $1 "* ]]; }

if has build; then
   sub build jobs/vista/g0_build.sbatch
   BUILD_JOB="${JOB}"
fi
DEP=()
if [ -n "${BUILD_JOB}" ] && [ "${BUILD_JOB}" != none ]; then
   DEP=(--dependency="afterok:${BUILD_JOB}" --kill-on-invalid-dep=yes)
fi

if has goldens; then
   for r in 1 2 4 8; do
      sub "golden p1 r${r}" "${DEP[@]}" -n "${r}" -t 04:00:00 \
         --export=ALL,ORDER=1 jobs/vista/g0_goldens.sbatch
   done
   sub "golden p2 r1" "${DEP[@]}" -n 1 -t 24:00:00 \
      --export=ALL,ORDER=2 jobs/vista/g0_goldens.sbatch
fi

if has timing; then
   for p in 1 2; do
      sub "timing p${p}" "${DEP[@]}" --export=ALL,ORDER="${p}" \
         jobs/vista/g0_timing.sbatch
   done
fi

if has counts; then
   for m in prod 500m 250m; do
      has_mesh "${m}" || continue
      for p in 1 2; do
         if [ "${m}" = 250m ] && [ "${p}" = 2 ]; then
            sub "counts ${m} p${p}" "${DEP[@]}" -N 4 -n 8 -t 06:00:00 \
               --export=ALL,MESH="${m}",ORDER="${p}",RANKS="4 8" \
               jobs/vista/g0_matrix_counts.sbatch
         else
            sub "counts ${m} p${p}" "${DEP[@]}" \
               --export=ALL,MESH="${m}",ORDER="${p}" \
               jobs/vista/g0_matrix_counts.sbatch
         fi
      done
   done
fi

if has factor; then
   for c in prod_p1 prod_p2; do
      has_mesh prod || continue
      sub "factor ${c}" "${DEP[@]}" -N 1 -n 144 -t 02:00:00 \
         --export=ALL,CASE="${c}" jobs/vista/g0_factor_size.sbatch
   done
   # case:nodes:ranks:walltime, smallest first.
   prev="${BUILD_JOB}"
   for spec in 500m_p1:2:72:02:00:00 250m_p1:8:288:04:00:00 \
               500m_p2:8:288:04:00:00 250m_p2:32:1152:06:00:00; do
      IFS=: read -r c nodes ranks hh mm ss <<< "${spec}"
      has_mesh "${c%_p*}" || continue
      chain=()
      if [ -n "${prev}" ] && [ "${prev}" != none ]; then
         chain=(--dependency="afterok:${prev}" --kill-on-invalid-dep=yes)
      fi
      sub "factor ${c}" "${chain[@]}" -N "${nodes}" -n "${ranks}" \
          -t "${hh}:${mm}:${ss}" --export=ALL,CASE="${c}" \
          jobs/vista/g0_factor_size.sbatch
      prev="${JOB}"
   done
fi

if has accst; then
   sub "accepted-state p1 r8" "${DEP[@]}" -n 8 --export=ALL,ORDER=1 \
      jobs/vista/g0_accepted_state.sbatch
fi

if [ "${SUBMIT}" != 1 ]; then
   echo "dry run: ${n} jobs printed, nothing submitted (SUBMIT=1 submits)."
fi
