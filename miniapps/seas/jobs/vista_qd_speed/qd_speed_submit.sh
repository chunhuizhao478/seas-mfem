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
#   integrator  QD speed-up plan v2, item A: jobs/vista_qd_speed/
#           qd_integrator.sbatch once per run below (1 node, 144 ranks,
#           compressed G read from the H file, t = 0 to 7.2e9 s), then
#           qd_integrator_compare.sbatch after all of them (afterany).
#           Reads the H file and stops unless it exists.  Not in the
#           default STEPS; INTEGRATOR_RUNS limits it (ref always runs).
#             ref  5dp, atol 1e-10        the reference
#             a0   5dp, atol 1e-7         today's setting
#             a1   5bs, atol 1e-7
#             a2   6vr, atol 1e-7
#             a3   5dp, atol 1e-8         with a0, today's cost-accuracy line
#             a4   5dp + DSP PI42, 1e-7
#   layout  QD speed item C (qd-c-design-rev1 section 5.2): the timing
#           script with MODES="hcount hbytes" and LOG_VIEW=1, so each job
#           runs the count layout and then the bytes layout on the same
#           nodes, compares their output files byte for byte and writes
#           PETSc's log (item D1).  The H files are only read, and the step
#           stops unless they exist.  Not in the default STEPS;
#           C_LAYOUTS limits it:
#             c_n1_r144     p = 1, 1,000 steps, 1 node   (count run also
#             c_n2_r288     p = 1, 1,000 steps, 2 nodes   compared with
#             c_n4_r576     p = 1, 1,000 steps, 4 nodes   1084131-33)
#             c_p2_n4_r576  p = 2, 200 steps, 4 nodes    (with 1086777)
#           EXCLUDE (default i617-131) keeps those jobs off nodes that
#           failed at startup.
#   gather  QD speed item B (qd-b-design-rev1 section 7): the timing script
#           with MODES="hrank hnode hnode8 hnode16", so each job runs the
#           bytes layout with the slip all-gather to every rank, then with
#           one slip copy per node, then that with V^T piece divisors 8 and
#           16, on the same nodes; every node run's output files must be
#           byte-identical to the rank run's.  The H files are only read,
#           and the step stops unless they exist.  Not in the default
#           STEPS; B_LAYOUTS limits it:
#             b_n1_r144     p = 1, 1,000 steps, 1 node
#             b_n2_r288     p = 1, 1,000 steps, 2 nodes
#             b_n4_r576     p = 1, 1,000 steps, 4 nodes
#             b_p2_n4_r576  p = 2, 200 steps, 4 nodes
#           EXCLUDE applies as for layout.
# Layouts (nodes x ranks, CPUs per rank):
#   n1_r144  1 x 144, 1   loops and the new default, twice each
#   n1_r72   1 x 72,  2   default
#   n1_r36   1 x 36,  4   default
#   n2_r288  2 x 144, 1   default
#   n2_r144  2 x 72,  2   default
#   n2_r72   2 x 36,  4   default
# LAYOUTS limits the timing jobs to those names.  Every job after the build
# waits for it (afterok) unless BUILD_JOB=none (build done) or
# BUILD_JOB=<jobid>.  AFTER="id:id..." holds the build job, and with it
# every job that waits for the build, until those jobs have started (Slurm
# after:), so a new study does not take queue priority from jobs already
# waiting.  Wall times are estimates from G0's 0.366 s/step at
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
INTEGRATOR_RUNS="${INTEGRATOR_RUNS:-a0 a1 a2 a3 a4}"
C_LAYOUTS="${C_LAYOUTS:-c_n1_r144 c_n2_r288 c_n4_r576 c_p2_n4_r576}"
B_LAYOUTS="${B_LAYOUTS:-b_n1_r144 b_n2_r288 b_n4_r576 b_p2_n4_r576}"
EXCLUDE="${EXCLUDE:-i617-131}"

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

AFTER="${AFTER:-}"
if has build; then
   if [ -n "${AFTER}" ]; then
      sub build --dependency="after:${AFTER}" jobs/vista/cpu_clone_build.sbatch
   else
      sub build jobs/vista/cpu_clone_build.sbatch
   fi
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

if has integrator; then
   if [ ! -f "${HCACHE_P1}" ]; then
      echo "ERROR: integrator: no H file ${HCACHE_P1}; the runs only read it." >&2
      exit 1
   fi
   # tag:options file:walltime.  Wall times from compressed G's 0.116
   # s/step at 6.2 stages per step (job 1084131) and gold 1074099's 63,557
   # steps to 7.2e9 s; ref takes about 4x the steps of a0 (atol^(1/5)).
   RUNS_DONE=""
   for spec in ref:tests/verification/petsc_ts_rk5dp_atol1e-10.cfg:12:00:00 \
               a0:tests/verification/petsc_ts_rk45_tandem.cfg:04:00:00 \
               a1:tests/verification/petsc_ts_rk5bs.cfg:05:00:00 \
               a2:tests/verification/petsc_ts_rk6vr.cfg:05:00:00 \
               a3:tests/verification/petsc_ts_rk5dp_atol1e-8.cfg:06:00:00 \
               a4:tests/verification/petsc_ts_rk5dp_dsp_pi42.cfg:04:00:00; do
      IFS=: read -r tag cfg hh mm ss <<< "${spec}"
      [ "${tag}" = ref ] || [[ " ${INTEGRATOR_RUNS} " == *" ${tag} "* ]] || continue
      sub "integrator ${tag}" "${DEP[@]}" -J "qd_int_${tag}" \
          --export="ALL,TSCFG=${cfg},TAG=${tag},HEPS=${HEPS},HCACHE=${HCACHE_P1}" \
          -N 1 -n 144 --ntasks-per-node=144 -c 1 \
          -t "${hh}:${mm}:${ss}" jobs/vista_qd_speed/qd_integrator.sbatch
      RUNS_DONE="${RUNS_DONE} ${tag}:${JOB}"
   done
   ids="$(for x in ${RUNS_DONE}; do printf ':%s' "${x##*:}"; done)"
   sub "integrator compare" --dependency="afterany${ids}" \
       --export="ALL,RUNS=${RUNS_DONE# },LINE=a0 a3" \
       jobs/vista_qd_speed/qd_integrator_compare.sbatch
fi

if has layout; then
   HCACHE_P2="$(dirname "${GREEN_CACHE_P1}")/bp5_prod_p2_h${HEPS}.bin"
   # name:order:steps:nodes:ranks:reference job:walltime.  Each job reads
   # the H file twice (p = 1 about 10 s, p = 2 about 19 s) and runs two
   # timing runs (p = 1 0.08-0.12 s/step, p = 2 0.22 s/step).
   C_SPECS="c_n1_r144:1:1000:1:144:1084131:00:45:00
            c_n2_r288:1:1000:2:288:1084132:00:30:00
            c_n4_r576:1:1000:4:576:1084133:00:30:00
            c_p2_n4_r576:2:200:4:576:1086777:00:45:00"
   # Every H file first, so a missing one submits none of these jobs
   for spec in ${C_SPECS}; do
      IFS=: read -r name order _ <<< "${spec}"
      [[ " ${C_LAYOUTS} " == *" ${name} "* ]] || continue
      hfile="${HCACHE_P1}"
      [ "${order}" = 1 ] || hfile="${HCACHE_P2}"
      if [ ! -f "${hfile}" ]; then
         echo "ERROR: layout ${name}: no H file ${hfile}; these jobs only" \
              "read it." >&2
         exit 1
      fi
   done
   for spec in ${C_SPECS}; do
      IFS=: read -r name order steps nodes ranks ref hh mm ss <<< "${spec}"
      [[ " ${C_LAYOUTS} " == *" ${name} "* ]] || continue
      hfile="${HCACHE_P1}"
      [ "${order}" = 1 ] || hfile="${HCACHE_P2}"
      sub "layout ${name}" "${DEP[@]}" -J qd_layout --exclude="${EXCLUDE}" \
          --export="ALL,MODES=hcount hbytes,ORDER=${order},RUN_STEPS=${steps},LOG_VIEW=1,LAYOUT_REF_JOB=${ref},HEPS=${HEPS},HCACHE=${hfile}" \
          -N "${nodes}" -n "${ranks}" --ntasks-per-node="$((ranks / nodes))" -c 1 \
          -t "${hh}:${mm}:${ss}" jobs/vista_qd_speed/qd_speed_timing.sbatch
   done
fi

if has gather; then
   HCACHE_P2="$(dirname "${GREEN_CACHE_P1}")/bp5_prod_p2_h${HEPS}.bin"
   # name:order:steps:nodes:ranks:walltime.  Each job reads the H file four
   # times (p = 1 about 10 s, p = 2 about 19 s) and runs four timing runs
   # (item C's bytes layout: p = 1 0.04-0.11 s/step, p = 2 0.14 s/step).
   B_SPECS="b_n1_r144:1:1000:1:144:00:45:00
            b_n2_r288:1:1000:2:288:00:30:00
            b_n4_r576:1:1000:4:576:00:30:00
            b_p2_n4_r576:2:200:4:576:00:45:00"
   # Every H file first, so a missing one submits none of these jobs
   for spec in ${B_SPECS}; do
      IFS=: read -r name order _ <<< "${spec}"
      [[ " ${B_LAYOUTS} " == *" ${name} "* ]] || continue
      hfile="${HCACHE_P1}"
      [ "${order}" = 1 ] || hfile="${HCACHE_P2}"
      if [ ! -f "${hfile}" ]; then
         echo "ERROR: gather ${name}: no H file ${hfile}; these jobs only" \
              "read it." >&2
         exit 1
      fi
   done
   for spec in ${B_SPECS}; do
      IFS=: read -r name order steps nodes ranks hh mm ss <<< "${spec}"
      [[ " ${B_LAYOUTS} " == *" ${name} "* ]] || continue
      hfile="${HCACHE_P1}"
      [ "${order}" = 1 ] || hfile="${HCACHE_P2}"
      sub "gather ${name}" "${DEP[@]}" -J qd_gather --exclude="${EXCLUDE}" \
          --export="ALL,MODES=hrank hnode hnode8 hnode16,ORDER=${order},RUN_STEPS=${steps},HEPS=${HEPS},HCACHE=${hfile}" \
          -N "${nodes}" -n "${ranks}" --ntasks-per-node="$((ranks / nodes))" -c 1 \
          -t "${hh}:${mm}:${ss}" jobs/vista_qd_speed/qd_speed_timing.sbatch
   done
fi
