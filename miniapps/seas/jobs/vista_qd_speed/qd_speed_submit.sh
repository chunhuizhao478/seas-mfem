#!/bin/bash
# QD speed items 1 and 2 on Vista: print (default) or submit the jobs.
#
#   cd $WORK/seas-mfem-qd-speed/miniapps/seas
#   bash jobs/vista_qd_speed/qd_speed_submit.sh            # dry run
#   SUBMIT=1 bash jobs/vista_qd_speed/qd_speed_submit.sh   # submit
#   SUBMIT=1 STEPS="timing" BUILD_JOB=none bash jobs/vista_qd_speed/qd_speed_submit.sh
#
# The clone is its own: $WORK/seas-mfem-qd-speed on branch
# quasi-dynamic-profiling-speedup.  The GPU thread's clones
# ($WORK/seas-mfem-gpu, $WORK/seas-mfem-gpu-g2) are only read, for G0's
# goldens (jobs/vista/cpu_clone_build.sbatch copies them once); PETSc and
# MUMPS are G0's in $WORK/seas-gpu-deps, not rebuilt.
#
# STEPS (default "build gate timing"):
#   build   jobs/vista/cpu_clone_build.sbatch: MFEM and seas_bp5_full in
#           this clone (1 node, up to 1.5 h)
#   gate    jobs/vista_qd_speed/qd_default_gate.sbatch: the new default
#           against G0's goldens, and which path it picks (1 node, 8 ranks)
#   timing  jobs/vista_qd_speed/qd_speed_timing.sbatch at every layout
#           below, 1,000 steps per run at p=1
#   green   the same script with MODES="green default" (QD speed item 3,
#           --qd-green) on 1, 2 and 4 full nodes: g_n1_r144 builds G
#           (37.5 GB at p=1) and writes the cache file (GREEN_CACHE in the
#           script); g_n2_r288 and g_n4_r576 wait for it (afterok) and
#           read it.  Each also runs the default solve at its layout.
#           Not in the default STEPS; GREEN_LAYOUTS limits it
#   events  the accuracy gate through earthquakes: the same script with
#           RUN_STEPS=62000 (gold run 1074099's events 1 and 2 end at its
#           accepted steps 10,538 and 57,742) and EVENT_REF set to gold
#           1074099 (bp5/benchmark_data/golden_vista_144r_p1_1000m), so
#           each run is compared with it by qd_event_compare.py.  One
#           mode per job, so the four run side by side:
#             e_n1_r144_loops    the old path: its difference to gold is
#                                the run-to-run floor (gold is loops, 144)
#             e_n1_r144_default  the new default (item 1)
#             e_n1_r144_green    dense G (item 3), same ranks as gold
#             e_n4_r576_green    dense G on 4 nodes (dip-0 stations
#                                depend on the rank count)
#           The green runs read GREEN_CACHE and must not build it (two
#           builds of one file at once): the step stops unless the file
#           exists.  Not in the default STEPS; EVENT_LAYOUTS limits it
#   compress  jobs/vista_qd_speed/qd_green_compress.sbatch: how far the
#           p = 1 cache compresses at round-off-level product error (QD
#           speed item 4, measurement only; 1 node, up to 3 h).  Reads
#           GREEN_CACHE; stops if it is missing.  Not in the default STEPS
#   hgreen  QD speed item 4, design rev2 T2: the timing script with
#           MODES="hgreen green" (compressed G, then dense G at the same
#           layout) on 1, 2 and 4 full nodes, 1,000 steps each.
#           h_n1_r144 compresses GREEN_CACHE at HEPS (default 1e-10) and
#           writes the H file (HCACHE in the timing script,
#           <GREEN_CACHE dir>/bp5_prod_p1_h<HEPS>.bin); h_n2_r288 and
#           h_n4_r576 wait for it (afterok) and read it.  Stops unless
#           GREEN_CACHE exists (the dense G is read, never built here).
#           Not in the default STEPS; HGREEN_LAYOUTS limits it
#   hevents design rev2 T3: the events step with MODES=hgreen, on 1 node
#           (gold's 144 ranks) and on 4 nodes, against gold 1074099.  Reads
#           the H file: waits for h_n1_r144 when hgreen is in the same
#           STEPS, else stops unless the file exists.  Not in the default
#           STEPS; HEVENT_LAYOUTS limits it
# Layouts (nodes x ranks, CPUs per rank):
#   n1_r144  1 x 144, 1   loops and the new default, twice each
#   n1_r72   1 x 72,  2   default
#   n1_r36   1 x 36,  4   default
#   n2_r288  2 x 144, 1   default
#   n2_r144  2 x 72,  2   default
#   n2_r72   2 x 36,  4   default
# LAYOUTS limits the timing jobs to those names.  Every job after the build
# waits for it (afterok) unless BUILD_JOB=none (build done) or
# BUILD_JOB=<jobid>.  Wall times are estimates from G0's 0.366 s/step at
# 144 ranks; fewer ranks per node are expected to be slower per step.

set -euo pipefail
SUBMIT="${SUBMIT:-0}"
STEPS="${STEPS:-build gate timing}"
BUILD_JOB="${BUILD_JOB:-}"
LAYOUTS="${LAYOUTS:-n1_r144 n1_r72 n1_r36 n2_r288 n2_r144 n2_r72}"
GREEN_LAYOUTS="${GREEN_LAYOUTS:-g_n1_r144 g_n2_r288 g_n4_r576}"
EVENT_LAYOUTS="${EVENT_LAYOUTS:-e_n1_r144_loops e_n1_r144_default e_n1_r144_green e_n4_r576_green}"
HGREEN_LAYOUTS="${HGREEN_LAYOUTS:-h_n1_r144 h_n2_r288 h_n4_r576}"
HEVENT_LAYOUTS="${HEVENT_LAYOUTS:-e_n1_r144_hgreen e_n4_r576_hgreen}"

if [ ! -f jobs/vista_qd_speed/qd_speed_timing.sbatch ]; then
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
has_layout() { [[ " ${LAYOUTS} " == *" $1 "* ]]; }

if has build; then
   sub build jobs/vista/cpu_clone_build.sbatch
   BUILD_JOB="${JOB}"
fi
DEP=()
if [ -n "${BUILD_JOB}" ] && [ "${BUILD_JOB}" != none ]; then
   DEP=(--dependency="afterok:${BUILD_JOB}" --kill-on-invalid-dep=yes)
fi

if has gate; then
   sub "default gate" "${DEP[@]}" jobs/vista_qd_speed/qd_default_gate.sbatch
fi

if has timing; then
   # name:nodes:ranks:cpus-per-task:walltime
   for spec in n1_r144:1:144:1:01:30:00 n1_r72:1:72:2:01:30:00 \
               n1_r36:1:36:4:02:00:00 n2_r288:2:288:1:01:00:00 \
               n2_r144:2:144:2:01:30:00 n2_r72:2:72:4:01:30:00; do
      IFS=: read -r name nodes ranks cpt hh mm ss <<< "${spec}"
      has_layout "${name}" || continue
      sub "timing ${name}" "${DEP[@]}" -N "${nodes}" -n "${ranks}" \
          --ntasks-per-node="$((ranks / nodes))" -c "${cpt}" \
          -t "${hh}:${mm}:${ss}" jobs/vista_qd_speed/qd_speed_timing.sbatch
   done
fi

if has green; then
   # name:nodes:ranks:walltime; the first one builds G, the others read it.
   GREEN_DEP=("${DEP[@]}")
   for spec in g_n1_r144:1:144:02:00:00 g_n2_r288:2:288:01:00:00 \
               g_n4_r576:4:576:01:00:00; do
      IFS=: read -r name nodes ranks hh mm ss <<< "${spec}"
      [[ " ${GREEN_LAYOUTS} " == *" ${name} "* ]] || continue
      sub "green ${name}" "${GREEN_DEP[@]}" -J qd_green_timing \
          --export="ALL,MODES=green default" -N "${nodes}" -n "${ranks}" \
          --ntasks-per-node="$((ranks / nodes))" -c 1 \
          -t "${hh}:${mm}:${ss}" jobs/vista_qd_speed/qd_speed_timing.sbatch
      if [ "${name}" = g_n1_r144 ]; then
         GREEN_DEP=(--dependency="afterok:${JOB}" --kill-on-invalid-dep=yes)
      fi
   done
fi

if has events; then
   GREEN_CACHE_P1="${GREEN_CACHE:-${SCRATCH:-${WORK}}/seas-qd-green/bp5_prod_p1.bin}"
   if [[ " ${EVENT_LAYOUTS} " == *"_green "* ]] && [ ! -f "${GREEN_CACHE_P1}" ]; then
      echo "ERROR: events: no Green's cache ${GREEN_CACHE_P1}; the green event" \
           "runs must read it, not build it." >&2
      exit 1
   fi
   EVENT_REF_SPEC="bp5/benchmark_data/golden_vista_144r_p1_1000m bp5_cycle_p1"
   # name:nodes:ranks:mode:walltime; estimates at 62,000 steps from the
   # 1,000-step timings (0.38, 0.26, 0.26 and 0.076 s/step) plus G's read.
   for spec in e_n1_r144_loops:1:144:loops:08:30:00 \
               e_n1_r144_default:1:144:default:06:00:00 \
               e_n1_r144_green:1:144:green:06:00:00 \
               e_n4_r576_green:4:576:green:02:30:00; do
      IFS=: read -r name nodes ranks mode hh mm ss <<< "${spec}"
      [[ " ${EVENT_LAYOUTS} " == *" ${name} "* ]] || continue
      sub "events ${name}" "${DEP[@]}" -J qd_event_gate \
          --export="ALL,MODES=${mode},RUN_STEPS=62000,EVENT_REF=${EVENT_REF_SPEC},GREEN_CACHE=${GREEN_CACHE_P1}" \
          -N "${nodes}" -n "${ranks}" --ntasks-per-node="$((ranks / nodes))" -c 1 \
          -t "${hh}:${mm}:${ss}" jobs/vista_qd_speed/qd_speed_timing.sbatch
   done
fi

if has compress; then
   GREEN_CACHE_P1="${GREEN_CACHE:-${SCRATCH:-${WORK}}/seas-qd-green/bp5_prod_p1.bin}"
   if [ ! -f "${GREEN_CACHE_P1}" ]; then
      echo "ERROR: compress: no Green's cache ${GREEN_CACHE_P1}." >&2
      exit 1
   fi
   sub "compress" --export="ALL,GREEN_CACHE=${GREEN_CACHE_P1}" \
       jobs/vista_qd_speed/qd_green_compress.sbatch
fi

# Item 4 (compressed G).  The H file's name follows the timing script's
# default, so both read the same file.
GREEN_CACHE_P1="${GREEN_CACHE:-${SCRATCH:-${WORK}}/seas-qd-green/bp5_prod_p1.bin}"
HEPS="${HEPS:-1e-10}"
HCACHE_P1="${HCACHE:-$(dirname "${GREEN_CACHE_P1}")/bp5_prod_p1_h${HEPS}.bin}"
H_WRITE_JOB=""
if has hgreen; then
   if [ ! -f "${GREEN_CACHE_P1}" ]; then
      echo "ERROR: hgreen: no Green's cache ${GREEN_CACHE_P1} to compress." >&2
      exit 1
   fi
   # name:nodes:ranks:walltime; h_n1_r144 also compresses (estimate: the
   # Python measurement took 2051 s on 36 processes) and writes the file.
   H_DEP=("${DEP[@]}")
   for spec in h_n1_r144:1:144:04:00:00 h_n2_r288:2:288:01:00:00 \
               h_n4_r576:4:576:01:00:00; do
      IFS=: read -r name nodes ranks hh mm ss <<< "${spec}"
      [[ " ${HGREEN_LAYOUTS} " == *" ${name} "* ]] || continue
      sub "hgreen ${name}" "${H_DEP[@]}" -J qd_hgreen_timing \
          --export="ALL,MODES=hgreen green,GREEN_CACHE=${GREEN_CACHE_P1},HEPS=${HEPS},HCACHE=${HCACHE_P1}" \
          -N "${nodes}" -n "${ranks}" --ntasks-per-node="$((ranks / nodes))" -c 1 \
          -t "${hh}:${mm}:${ss}" jobs/vista_qd_speed/qd_speed_timing.sbatch
      if [ "${name}" = h_n1_r144 ]; then
         H_WRITE_JOB="${JOB}"
         H_DEP=(--dependency="afterok:${JOB}" --kill-on-invalid-dep=yes)
      fi
   done
fi

if has hevents; then
   if [ -n "${H_WRITE_JOB}" ]; then
      HE_DEP=(--dependency="afterok:${H_WRITE_JOB}" --kill-on-invalid-dep=yes)
   elif [ -f "${HCACHE_P1}" ]; then
      HE_DEP=("${DEP[@]}")
   else
      echo "ERROR: hevents: no H file ${HCACHE_P1}, and no hgreen h_n1_r144" \
           "in this submission to write it." >&2
      exit 1
   fi
   EVENT_REF_SPEC="bp5/benchmark_data/golden_vista_144r_p1_1000m bp5_cycle_p1"
   # name:nodes:ranks:walltime; estimates at 62,000 steps from dense G's
   # event gate (0.293 s/step on 1 node, 0.076 on 4) plus the H read.
   for spec in e_n1_r144_hgreen:1:144:06:30:00 e_n4_r576_hgreen:4:576:02:00:00; do
      IFS=: read -r name nodes ranks hh mm ss <<< "${spec}"
      [[ " ${HEVENT_LAYOUTS} " == *" ${name} "* ]] || continue
      sub "hevents ${name}" "${HE_DEP[@]}" -J qd_hevent_gate \
          --export="ALL,MODES=hgreen,RUN_STEPS=62000,EVENT_REF=${EVENT_REF_SPEC},GREEN_CACHE=${GREEN_CACHE_P1},HEPS=${HEPS},HCACHE=${HCACHE_P1}" \
          -N "${nodes}" -n "${ranks}" --ntasks-per-node="$((ranks / nodes))" -c 1 \
          -t "${hh}:${mm}:${ss}" jobs/vista_qd_speed/qd_speed_timing.sbatch
   done
fi
