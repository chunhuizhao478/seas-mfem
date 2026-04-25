#!/bin/bash
#SBATCH -J tpv104_bulk_asym
#SBATCH -o tpv104_bulk_asym_%j.out
#SBATCH -e tpv104_bulk_asym_%j.err
#SBATCH -p development
#SBATCH -N 8
#SBATCH -n 400
#SBATCH -t 02:00:00
#SBATCH -A EAR20006

# TPV104 bulk-side asymmetry probe — t = 2 s, 200 m, 400 ranks
# ============================================================
#
# Purpose
# -------
# Identify whether the σ_n stress-channel asymmetry is born in the
# bulk wave operator (interior non-fault face Godunov flux on the
# non-mirror mesh) or self-perpetuated by the fault-Riemann back-
# rotation injecting asymmetric SYY/SZZ into the bulk.
#
# Background (from earlier analysis of normaltrace job 7677547)
# -------------------------------------------------------------
# Reading the existing C-1s ADER-INT trace at the hypocenter QP shows:
#   - GLOB-IN[SXX]+ vs - already differs by 5.97 MPa at peak step 3734
#   - LOC-IN[SXX]   asymmetry  ≈ 0.35 MPa  (rotation attenuates 17×)
#   - LOC-IMP[SXX]  asymmetry  = 0  (BuildImposedState symmetric — clean)
#   - GLOB-OUT[SXX] asymmetry  = 5.97 MPa  (re-acquires GLOB-IN level
#                                          via untouched SYY/SZZ channels)
# Conclusion: stages 1, 2, 3 of the fault-Riemann pipeline are correct;
# the asymmetry is in the bulk Q before fault code touches it.  This
# sbatch instruments the bulk-side path to localize WHICH non-fault
# face delivers the asymmetric flux.
#
# New probes (gated by SEAS_DIAG_FAULT_FLUX, all silent unset)
# ------------------------------------------------------------
#   C-2A BULK-DOF        : per-Mult-call dump of bulk Q and the pre-
#                          mass-inverse RHS k at the diag DOFs of the
#                          two adjacent tets (E_plus / E_minus of the
#                          hypocenter face).  Reads raw bulk DOFs (no
#                          face-QP interpolation); single rank only.
#   C-2B NONFAULT-FACE   : on each of the 6 non-fault interior faces of
#                          the two diag tets, dump (face_id, e1, e2,
#                          nor, w, Q_self, Q_nbr, F_h) on the FIRST QP
#                          per face per Mult call.  Identifies which
#                          face's flux is asymmetric and by how much.
#   C-2C BULK-DELTA      : per-macro-step delta Q^{n+1} - Q^n at the
#                          two diag DOFs.  One line per macro-step;
#                          enables direct correlation with C-1n peak.
#
# Together with the existing C-1s ADER-INT trace, the four probes give
# a mass balance:
#
#   Δ_E_plus_SXX_dof  =  Σ_{6 non-fault faces} (E_plus side flux contrib)
#                       + (fault-face GLOB-OUT+ contribution, from C-1s)
#   Δ_E_minus_SXX_dof =  same on E_minus side
#
#   Δ(E_plus - E_minus)_SXX_dof = a + b
#       a = Σ_nonfault (E+ flux - E- flux)  ← bulk source
#       b = fault GLOB-OUT (+ vs -)        ← fault back-injection
#
# Expected outcomes
# -----------------
#   a >> b  : bulk wave operator on non-mirror mesh dominates →
#             σ_n freeze does not address root cause; need mirror
#             mesh or a different decoupling test.
#   b >> a  : fault back-injection self-perpetuates → σ_n freeze
#             intercepts the wrong channel, but a "freeze SYY/SZZ"
#             override would address it.
#   a ~ b   : both contribute; per-face C-2B identifies which non-
#             fault tet pair drives `a`.
#
# Mesh / dispatch parity with the normaltrace baseline
# ----------------------------------------------------
#   - tpv104_200m.msh (same mesh as normaltrace_dev_7677547 →
#     mesh-asymmetry contribution is HELD CONSTANT, only the probe
#     surface changes).
#   - 8 nodes × 400 ranks (same).
#   - --ader-order 2, --order 1, --cfl 0.5 (same).
#   - --tfinal 2.0 (same window, direct overlay vs normaltrace).
#   - --friction-solver newton-stable, --fric-law slip-srw,
#     --fault-iterator substep (banner-only; production path is
#     Brent + one-shot ADER + macro-step ψ per R7-001(b)).
#
# Run economics
# -------------
#   Stderr volume estimate:
#     C-2A : 4 lines × ~7000 macro-steps × ~250 chars =  7 MB
#     C-2B : 6 faces × 4 lines × 7000 steps × ~300 chars = 50 MB
#     C-2C : 1 line × 7000 steps × 300 chars = 2 MB
#     C-1n / C-1s ADER-INT (existing) ≈ 18 MB
#     ----------------------------------------- TOTAL ≈ 77 MB
#   Within Frontera limits.  Wall time matches normaltrace (~30 min
#   at 400 ranks with the diag printf overhead).
#
# Pre-requisites
# --------------
#   1. Source on Frontera carries the C-2 probe blocks.  Verified at
#      startup by grep'ing for "C-2A BULK-DOF" in wave_operator.inl.
#   2. tpv104/mesh/tpv104_200m.msh exists.
#   3. tpv104_init_stations.sbatch PASS.
#   4. f_w = 0.2 (post-fix).

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

# Source verification — make sure the C-2 probe blocks are present.
if ! grep -q 'C-2A BULK-DOF' dynamic/wave_operator.inl; then
   echo "ERROR: dynamic/wave_operator.inl does not carry the C-2A probe."
   echo "       git pull on Frontera before submitting."
   exit 1
fi
if ! grep -q 'C-2B NONFAULT-FACE' dynamic/wave_operator.inl; then
   echo "ERROR: dynamic/wave_operator.inl does not carry the C-2B probe."
   exit 1
fi
if ! grep -q 'C-2C BULK-DELTA' drivers/tpv104_driver.cpp; then
   echo "ERROR: drivers/tpv104_driver.cpp does not carry the C-2C probe."
   exit 1
fi
echo "C-2 bulk probes confirmed in source."

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

RESULT_DIR="tpv104/results_bulk_asym_job${SLURM_JOB_ID}"
mkdir -p "${RESULT_DIR}"

# Force rebuild with SEAS_DIAG_FAULT_FLUX defined.  All four probe
# blocks (C-1 EVAL, C-1n NORMAL, C-1s ADER-INT, C-2A/B/C) live behind
# this one flag.
rm -f seas_tpv104_driver drivers/tpv104_driver.o \
      dynamic/tpv104_substep_iterator.o dynamic/friction_solver.o \
      dynamic/fault_face_flux.o dynamic/wave_operator.o
make SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX" seas_tpv104_driver -j"${JOBS}" \
   || { echo "BUILD FAILED"; exit 1; }

if [ ! -s tpv104/mesh/tpv104_200m.msh ]; then
   echo "ERROR: tpv104/mesh/tpv104_200m.msh missing."
   exit 1
fi

# IMPORTANT: SEAS_TPV104_FREEZE_SIGMA_N is UNSET — we want the natural
# asymmetry, not a frozen one.
unset SEAS_TPV104_FREEZE_SIGMA_N

echo "========================================"
echo "Git commit:   $(git rev-parse HEAD)"
echo "Git branch:   $(git rev-parse --abbrev-ref HEAD)"
echo "Output dir:   ${RESULT_DIR}"
echo "Target:       tfinal = 2.0 s on tpv104_200m.msh at 400 ranks"
echo "Probes:       C-1 EVAL + C-1n NORMAL + C-1s ADER-INT + C-2A/B/C"
echo "freeze:       UNSET (natural asymmetry capture)"
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
      --output-prefix tpv104_bulk_asym \
      --output-dt 0.05

# Capture the C-2 probe output from the slurm .err into dedicated
# per-stream logs inside the result dir.  Slurm appends all ranks'
# stderr to a single .err; the C-2 probes only fire from one rank.
JOB_ERR="${SBATCH_LOG_DIR}/tpv104_bulk_asym_${SLURM_JOB_ID}.err"
PROBE_LOG_C1n="${RESULT_DIR}/c1n_probe.log"
PROBE_LOG_C1s="${RESULT_DIR}/c1s_probe.log"
PROBE_LOG_C2A="${RESULT_DIR}/c2a_bulk_dof.log"
PROBE_LOG_C2B="${RESULT_DIR}/c2b_nonfault_face.log"
PROBE_LOG_C2C="${RESULT_DIR}/c2c_bulk_delta.log"

if [ -f "${JOB_ERR}" ]; then
   grep -E "C-1n NORMAL"            "${JOB_ERR}" > "${PROBE_LOG_C1n}"
   grep -E "C-1s ADER-INT|diag-c2"  "${JOB_ERR}" > "${PROBE_LOG_C1s}"
   grep -E "C-2A BULK-DOF"          "${JOB_ERR}" > "${PROBE_LOG_C2A}"
   grep -E "C-2B NONFAULT-FACE"     "${JOB_ERR}" > "${PROBE_LOG_C2B}"
   grep -E "C-2C BULK-DELTA"        "${JOB_ERR}" > "${PROBE_LOG_C2C}"
   echo "Captured probe logs into ${RESULT_DIR}"
   echo "  c1n: $(wc -l < ${PROBE_LOG_C1n}) lines"
   echo "  c1s: $(wc -l < ${PROBE_LOG_C1s}) lines"
   echo "  c2a: $(wc -l < ${PROBE_LOG_C2A}) lines"
   echo "  c2b: $(wc -l < ${PROBE_LOG_C2B}) lines"
   echo "  c2c: $(wc -l < ${PROBE_LOG_C2C}) lines"
else
   echo "WARN: ${JOB_ERR} not found at post-run check"
fi

# Mass-balance analysis: at peak step (max |sigma_n_trial| from C-1n),
# compute the per-component asymmetry contributions.
RESULT_FILE="${RESULT_DIR}/RESULT.txt"
HYPO_FILE="${RESULT_DIR}/tpv104_bulk_asym_station_x2_0_x3_7.5.dat"

PEAK_V_HYPO="N/A"
if [ -s "${HYPO_FILE}" ]; then
   PEAK_V_HYPO=$(awk '!/^#/ && NF>=3 && $1+0==$1 {v=$3; if(v<0)v=-v; if(v>m)m=v} END{print m+0}' "${HYPO_FILE}")
fi

# C-2C bulk-delta scan: max |Q+_SXX - Q-_SXX| at face DOFs over the run.
C2C_REPORT="(no C-2C data)"
if [ -s "${PROBE_LOG_C2C}" ]; then
   C2C_REPORT=$(awk '
      /\[C-2C BULK-DELTA\]/ {
         for (i=1;i<=NF;i++) {
            if ($i ~ /^diff_SXX=/)   { gsub("diff_SXX=","",$i);   d=$i+0; ad=(d<0?-d:d); if(ad>mx_sxx){mx_sxx=ad;t_mx_sxx=$3} }
            if ($i ~ /^diff_SYY=/)   { gsub("diff_SYY=","",$i);   d=$i+0; ad=(d<0?-d:d); if(ad>mx_syy){mx_syy=ad;t_mx_syy=$3} }
            if ($i ~ /^diff_SZZ=/)   { gsub("diff_SZZ=","",$i);   d=$i+0; ad=(d<0?-d:d); if(ad>mx_szz){mx_szz=ad;t_mx_szz=$3} }
         }
         n++
      }
      END {
         if (n==0) { print "no C-2C lines parsed"; exit }
         printf "  macro-steps             : %d\n", n;
         printf "  max |Q+_SXX - Q-_SXX|   : %.4e Pa  at t %s\n", mx_sxx, t_mx_sxx;
         printf "  max |Q+_SYY - Q-_SYY|   : %.4e Pa  at t %s\n", mx_syy, t_mx_syy;
         printf "  max |Q+_SZZ - Q-_SZZ|   : %.4e Pa  at t %s\n", mx_szz, t_mx_szz;
      }
   ' "${PROBE_LOG_C2C}")
fi

# C-2B per-face flux summary: which face delivers the largest |F_h|
# asymmetry over the run?  Group by face_id, max abs(F_h[SXX]).
C2B_REPORT="(no C-2B data)"
if [ -s "${PROBE_LOG_C2B}" ]; then
   C2B_REPORT=$(awk '
      /\[C-2B NONFAULT-FACE F_h\]/ {
         face=""
         for (i=1;i<=NF;i++) if ($i ~ /^face=/) { gsub("face=","",$i); face=$i }
         for (i=1;i<=NF;i++) if ($i ~ /^SXX=/) { gsub("SXX=","",$i); v=$i+0; av=(v<0?-v:v); if(av>face_max[face]){face_max[face]=av} }
         face_count[face]++
      }
      END {
         printf "  unique non-fault faces seen : %d\n", length(face_max);
         for (f in face_max) {
            printf "    face=%s  max|F_h_SXX| = %.4e Pa  (%d samples)\n",
                   f, face_max[f], face_count[f];
         }
      }
   ' "${PROBE_LOG_C2B}")
fi

{
   echo "TPV104 bulk-side asymmetry probe run"
   echo "Mesh:                  tpv104_200m"
   echo "tfinal:                2.0 s"
   echo "freeze:                UNSET"
   echo "f_w (compiled):        $(echo "${F_W_LINE}" | grep -oE '[01]\.[0-9]+' | head -1)"
   echo ""
   echo "peak |V_strike| at hypo  : ${PEAK_V_HYPO} m/s"
   echo "  (post-fix expected ~7.4 m/s, matching fwfix run)"
   echo ""
   echo "C-2C bulk-delta summary:"
   echo "${C2C_REPORT}"
   echo ""
   echo "C-2B per-face summary:"
   echo "${C2B_REPORT}"
   echo ""
   echo "Probe logs:"
   echo "  ${PROBE_LOG_C1n}"
   echo "  ${PROBE_LOG_C1s}"
   echo "  ${PROBE_LOG_C2A}"
   echo "  ${PROBE_LOG_C2B}"
   echo "  ${PROBE_LOG_C2C}"
   echo ""
   echo "Recommended post-run analysis:"
   echo "  1. Pair C-2A and C-1s ADER-INT GLOB-IN per macro-step:"
   echo "     C-2A reads bulk Q at the face DOF; C-1s GLOB-IN reads"
   echo "     bulk Q at the face QP (shape-interpolated).  The two"
   echo "     should agree on the asymmetry signal."
   echo "  2. Sum C-2B F_h × w contributions to elem_plus across the"
   echo "     6 non-fault faces; compare to C-2A k+ delta per Mult."
   echo "     Discrepancy = Mult contributions outside this face loop"
   echo "     (volume integral, fault face, BCs)."
   echo "  3. Compute a vs b mass-balance:"
   echo "       a = Σ_nonfault (E+ flux - E- flux)       [bulk-born]"
   echo "       b = fault_GLOB-OUT (+ vs -)              [fault-born]"
   echo "     a/b ratio at peak → the diagnostic verdict."
   echo ""
   if [ -s "${PROBE_LOG_C2A}" ] && [ -s "${PROBE_LOG_C2B}" ] && \
      [ -s "${PROBE_LOG_C2C}" ] && [ "${PEAK_V_HYPO}" != "N/A" ]; then
      echo "STATUS: PASS"
   else
      echo "STATUS: FAIL"
   fi
} > "${RESULT_FILE}"
cat "${RESULT_FILE}"
