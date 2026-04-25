#!/bin/bash
#SBATCH -J tpv104_normaltrace
#SBATCH -o tpv104_normaltrace_%j.out
#SBATCH -e tpv104_normaltrace_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 σ_n trial-traction trace — root-cause diagnostic for the
# fault-normal stress perturbation observed at the hypocenter.
# =====================================================
#
# Background
# ----------
# The f_w = 0.1 → 0.2 fix (commit 4d05df2) collapsed the V_strike
# over-shoot at the hypocenter to within 0.35% of the SeisSol reference.
# But three secondary symptoms remain at the hypocenter:
#   (a) σ_n perturbation ±0.5 MPa during rupture transit (SeisSol shows
#       flat 120 MPa).
#   (b) V_dip leakage ±0.03 m/s (SeisSol ~0).
#   (c) slip_dip accumulates 3.5 mm by t=2s (SeisSol ~0).
# All three are signatures of strike→{normal, dip} channel coupling that
# should not exist for a planar strike-slip rupture.
#
# This run dumps the σ_n trial-traction decomposition at every macro-
# step at the hypocenter QP using the new C-1n NORMAL probe in
# `dynamic/fault_face_flux.cpp:271-303`.  σ_n_trial =
#   etaP * [(Q_minus[VX] − Q_plus[VX])               ← v_jump_term
#         + (Q_plus[SXX]/Zp + Q_minus[SXX]/Zp_neig)] ← stress_term
# Tells us which channel of bulk Q is being populated with non-zero
# fault-normal content (velocity jump or stress average).
#
# Hypotheses & diagnostic outcome:
#   v_jump_term dominant   → bulk wave operator producing spurious
#                            fault-normal velocity at fault QPs (most
#                            likely from mesh non-mirror-symmetry; the
#                            mesh-mirror script confirmed 0/215143 tets
#                            have a mirror partner).
#   stress_term dominant   → bulk Q is acquiring non-zero σ_yy
#                            fluctuation near the fault, either from
#                            the imposed-state back-flux (Riemann path)
#                            or from numerical dispersion.
#   Both ~equal magnitude  → both mechanisms active simultaneously;
#                            need a third diagnostic to disambiguate
#                            (e.g. reset Q_plus/Q_minus normal-channel
#                            entries before friction call, see if the
#                            σ_n leak persists).
#
# Run parameters
# --------------
# Same setup as `tpv104_fw_fix_validation.sbatch` (which produced the
# clean V_strike PASS we want to preserve).  Only changes:
#   - Build with `MFEM_CPPFLAGS=-DSEAS_DIAG_FAULT_FLUX` (enables the
#     C-1n probe; otherwise the printf is #ifdef'd out and run looks
#     like the production path).
#   - Stderr captured into the result dir as `c1n_probe.log` so we can
#     post-process the trace without interleaving from regular driver
#     stdout.
#   - --output-prefix tpv104_normaltrace (avoids clobbering fwfix run).
#   - --tfinal 2.0 (covers full hypocenter rupture transit; matches
#     fwfix run's window for direct comparison of the same QP).
#
# Cost: same as fwfix run (~30 min wall on dev queue at 400 ranks).
#
# C-1n probe MPI safety
# ---------------------
# `data.diag_print` is set on exactly ONE DOF on the rank that owns the
# closest local QP to the hypocenter (via MPI_Allreduce MINLOC in the
# driver's SEAS_DIAG_FAULT_FLUX init block).  Other ranks have zero
# tagged DOFs.  The probe block runs only when data.diag_print is true,
# so:
#   - exactly one rank emits C-1n NORMAL lines,
#   - probe is single-process stderr printf, no MPI calls,
#   - no collective communication, no deadlock possible at any rank
#     count.
# Probe overhead ~9 MB stderr output total over t=2s — negligible.
#
# Pre-requisites
# --------------
#   1. tpv104/mesh/tpv104_200m.msh exists (built locally + scp'd, or
#      via tpv104_mesh_build.sbatch).
#   2. tpv104_init_stations.sbatch PASS — confirms ψ_ini gate.
#   3. config/tpv104_params.hpp has f_w = 0.2 (post-fix; verified at
#      startup).

export LC_ALL=C
export LANG=C

SBATCH_LOG_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

module load intel/19.1.1
module load impi/19.0.9
module load hypre/2.31.0
module load mumps/5.3
module load parmetis
module load petsc/3.15
module load fftw3/3.3.8

export LD_LIBRARY_PATH="${TACC_HYPRE_LIB}:${TACC_PARMETIS_LIB}:${TACC_MUMPS_LIB}:${TACC_PETSC_LIB}:${TACC_FFTW3_LIB}:${LD_LIBRARY_PATH}"

cd /scratch2/10024/zhaochun/seas-project/seas-mfem

JOBS=8

grep -q '^MFEM_USE_PETSC *= YES' config/config.mk || { echo "ERROR: no PETSc"; exit 1; }
grep -q '^MFEM_USE_MUMPS *= YES' config/config.mk || { echo "ERROR: no MUMPS"; exit 1; }

cd miniapps/seas

# Verify f_w fix landed (post-4d05df2 source).
F_W_LINE=$(grep -E '^\s*static\s+constexpr\s+real_t\s+f_w' config/tpv104_params.hpp | head -1)
if ! echo "${F_W_LINE}" | grep -q '0\.2'; then
   echo "ERROR: config/tpv104_params.hpp f_w is NOT 0.2 — pull latest source."
   exit 1
fi
echo "f_w confirmed: ${F_W_LINE}"

PREV_INIT_RESULT=$(ls -dt tpv104/results_init_stations_job* 2>/dev/null | head -1)
if [ -z "${PREV_INIT_RESULT}" ] \
   || ! grep -q "^STATUS: PASS" "${PREV_INIT_RESULT}/RESULT.txt" 2>/dev/null; then
   echo "ERROR: no passing tpv104_init_stations.sbatch result found."
   exit 1
fi
echo "Found init PASS at: ${PREV_INIT_RESULT}"

RESULT_DIR="tpv104/results_normaltrace_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

# Force rebuild with SEAS_DIAG_FAULT_FLUX defined.  The C-1n NORMAL
# probe is gated by this preprocessor flag in fault_face_flux.cpp, AND
# the hypocenter-tagging block (data.diag_print = true) is gated by it
# in tpv104_driver.cpp — both .o files must be rebuilt with the flag.
rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o
make SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX" seas_tpv104_driver -j"${JOBS}" \
   || { echo "BUILD FAILED"; exit 1; }

if [ ! -s tpv104/mesh/tpv104_200m.msh ]; then
   echo "ERROR: tpv104/mesh/tpv104_200m.msh missing."
   exit 1
fi

echo "========================================"
echo "Git commit:   $(git rev-parse HEAD)"
echo "Git branch:   $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:   ${RESULT_DIR}"
echo "Target:       tfinal = 2.0 s on tpv104_200m.msh at 400 ranks"
echo "Probe:        C-1n NORMAL  (σ_n_trial decomposition at hypocenter)"
echo "========================================"

ibrun ./seas_tpv104_driver \
      --mesh tpv104/mesh/tpv104_200m.msh \
      --mesh-scale 1 \
      --order 1 \
      --bc-mode absorbing \
      --cfl 0.5 \
      --tfinal 2.0 \
      --ader-order 2 \
      --friction-solver newton-stable \
      --fric-law slip-srw \
      --fault-iterator substep \
      --output-dir "${RESULT_DIR}" \
      --output-prefix tpv104_normaltrace \
      --output-dt 0.05

# Capture the C-1n probe output from the slurm .err into a dedicated
# log inside the result dir.  Slurm appends all ranks' stderr to a
# single .err file; the C-1n probe only fires from one rank (the
# diag-tagged rank), so all C-1n NORMAL lines come from that rank.
JOB_ERR="${SBATCH_LOG_DIR}/tpv104_normaltrace_${SLURM_JOB_ID}.err"
PROBE_LOG="${RESULT_DIR}/c1n_probe.log"
if [ -f "${JOB_ERR}" ]; then
   grep -E "C-1 EVAL|C-1n NORMAL|\[diag\]" "${JOB_ERR}" > "${PROBE_LOG}"
   echo "Captured $(wc -l < "${PROBE_LOG}") probe lines into ${PROBE_LOG}"
else
   echo "WARN: ${JOB_ERR} not found at post-run check"
fi

# Post-run analysis: classify the σ_n_trial decomposition at the
# moment of peak σ_n perturbation.  Awk over the C-1n NORMAL lines,
# track max |sigma_n_trial| and the corresponding v_jump_term /
# stress_term values.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_normaltrace_station_x2_0_x3_7.5.dat"

PEAK_V_HYPO="N/A"
if [ -s "${HYPO_FILE}" ]; then
   PEAK_V_HYPO=$(awk '!/^#/ && NF>=3 && $1+0==$1 {v=$3; if(v<0)v=-v; if(v>m)m=v} END{print m+0}' "${HYPO_FILE}")
fi

# Probe-line analysis.  Format:
#   [C-1n NORMAL] rank=R  sigma_n_trial=±X Pa  v_jump_term=±Y Pa  stress_term=±Z Pa  | Q_plus[VX]=±A  Q_minus[VX]=±B  Q_plus[SXX]=±C  Q_minus[SXX]=±D
# Track |sigma_n_trial| max and the corresponding v_jump / stress fractions.
PROBE_REPORT="(no probe data)"
if [ -s "${PROBE_LOG}" ]; then
   PROBE_REPORT=$(awk '
      /\[C-1n NORMAL\]/ {
         # Extract numeric values via gsub.
         line = $0;
         match(line, /sigma_n_trial=([+-][0-9.e+-]+)/, a); sn = a[1]+0;
         match(line, /v_jump_term=([+-][0-9.e+-]+)/, a);  vj = a[1]+0;
         match(line, /stress_term=([+-][0-9.e+-]+)/, a);  st = a[1]+0;
         abs_sn = sn<0 ? -sn : sn;
         if (abs_sn > maxabs) {
            maxabs = abs_sn;
            sn_at_peak = sn;
            vj_at_peak = vj;
            st_at_peak = st;
         }
         n++;
      }
      END {
         if (n == 0) { print "no C-1n lines parsed"; exit }
         abs_vj = vj_at_peak<0 ? -vj_at_peak : vj_at_peak;
         abs_st = st_at_peak<0 ? -st_at_peak : st_at_peak;
         frac_vj = abs_vj / maxabs;
         frac_st = abs_st / maxabs;
         dominant = (abs_vj > abs_st) ? "v_jump_term (velocity-jump)" : "stress_term (stress-average)";
         printf "  probe lines           : %d\n", n
         printf "  peak |sigma_n_trial|  : %.4e Pa\n", maxabs
         printf "  at peak: v_jump_term  : %+.4e Pa  (%.1f%% of |sigma_n_trial|)\n",
                vj_at_peak, 100*frac_vj
         printf "  at peak: stress_term  : %+.4e Pa  (%.1f%% of |sigma_n_trial|)\n",
                st_at_peak, 100*frac_st
         printf "  dominant channel      : %s\n", dominant
      }
   ' "${PROBE_LOG}")
fi

{
   echo "TPV104 σ_n trial-traction trace run"
   echo "Mesh:                  tpv104_200m"
   echo "tfinal:                2.0 s"
   echo "f_w (compiled):        $(echo "${F_W_LINE}" | grep -oE '[01]\.[0-9]+' | head -1)"
   echo ""
   echo "peak |V_strike| at hypo : ${PEAK_V_HYPO} m/s"
   echo "  (post-fix expected ~7.4 m/s, matching fwfix run)"
   echo ""
   echo "C-1n NORMAL probe trace at hypocenter:"
   echo "${PROBE_REPORT}"
   echo ""
   echo "Probe log: ${PROBE_LOG}"
   echo "  Inspect with:  grep 'C-1n NORMAL' ${PROBE_LOG} | head -50"
   echo ""
   echo "Diagnostic interpretation (per sbatch comment):"
   echo "  v_jump_term dominant  → mesh non-mirror-symmetry / bulk-wave-operator"
   echo "                          spurious fault-normal velocity at fault QPs."
   echo "  stress_term dominant  → bulk Q acquiring σ_yy fluctuation, likely"
   echo "                          from imposed-state back-flux Riemann path."
   echo "  both comparable       → both mechanisms active; needs further isolation."
   echo ""
   if [ -s "${PROBE_LOG}" ] && [ "${PEAK_V_HYPO}" != "N/A" ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
