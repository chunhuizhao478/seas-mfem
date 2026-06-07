// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// godunov_flux_bimaterial.cpp — Phase R.1 of
// PLAN_phase_R_exact_bimaterial_riemann_rev3.md.
//
// Implementation note (deviation from godunov_flux_bimaterial.hpp's
// initial header comment): after empirical testing against the R.1.T-2
// bi-material P-wave transmission test, switched from MFEM's
// GodunovFlux R sign convention to SeisSol's matR sign convention.
//
// Why: MFEM's R for the LOCAL +cp eigenvector has SXX = -lp (negative);
// SeisSol's matR has SXX = +lp (positive).  They differ by a column
// sign flip.  This sign matters in the bi-material case because the
// matrix is built from MIXED material parameters (cols 0..2 use LOCAL,
// cols 6..8 use NEIGHBOUR) — a sign flip on individual columns
// changes the SUBSPACE structure and therefore the projector
// matR · chi · matR^{-1}.
//
// With SeisSol's sign convention AND the formula
//   qGodLocal    = I - godunov     (multiplier of Q_L)
//   qGodNeighbor =     godunov     (multiplier of Q_R)
//   Q* = qGodLocal · Q_L + qGodNeighbor · Q_R
// the bi-material formula:
//   (a) reduces to F = A^+ Q_L + A^- Q_R in the homogeneous limit
//       (byte-exact with `GodunovFlux::Interior` modulo LU rounding),
//   (b) reproduces the analytic interface state in the bi-material
//       case for a right-going P-wave incident on the LEFT side:
//         Q*[VX]  = (2 Z_L      / (Z_L + Z_R)) · v_n^L   (velocity)
//         Q*[SXX] = -(2 Z_L Z_R / (Z_L + Z_R)) · v_n^L   (stress)
//       The velocity factor 2 Z_L / (Z_L + Z_R) is "Tv" (velocity
//       transmission); the stress factor 2 Z_R / (Z_L + Z_R) is
//       "Tp" (stress / pressure transmission).  The L and R
//       impedances swap places between the two channels — both
//       formulas are needed, depending on whether the test reads
//       Q*[VX] or Q*[SXX].
//
// SUBTLETY about WHY the SeisSol sign gives standard upwinding:
// using SeisSol's matR (POSITIVE stress entries), col 0 — labelled
// "LOCAL +cp outgoing" — is actually an eigenvector of MFEM's A_x
// with eigenvalue −cp_L (not +cp_L) under MFEM's stress-flux sign
// convention `dQ/dt + A · ∂Q/∂x = 0`.  The chi-projection then selects
// the MFEM-NEGATIVE-eigenvalue subspace, so godunov = matR · chi ·
// matR^{-1} satisfies A · godunov = A^- (NOT A^+).  Hence
// A · (I - godunov) = A - A^- = A^+, and
//   F = A · Q* = A · ((I-godunov) Q_L + godunov Q_R) = A^+ Q_L + A^- Q_R
// — standard upwind.  The "labelling vs. physics" inversion is purely
// a notational artefact of SeisSol's eigenvector sign choice; the
// resulting math gives the correct upwind flux + correct bi-material
// physics.
//
// References (re-verified at implementation time):
//   /Users/chunhuizhao/projects/SeisSol/src/Equations/elastic/Model/
//   ElasticSetup.h:79-166 — `getTransposedGodunovState`.

#include "godunov_flux_bimaterial.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace mfem
{
namespace seas
{

// ---------------------------------------------------------------------------
// BuildGodunovStateFaceLocal — see file header for the sign-convention
// derivation.  Implementation mirrors SeisSol's getTransposedGodunovState
// non-acoustic branch (ElasticSetup.h:102-139), translated to MFEM's
// column-vector convention.
// ---------------------------------------------------------------------------
/* static */ void BimaterialFlux::BuildGodunovStateFaceLocal(
   real_t lam_L, real_t mu_L, real_t rho_L,
   real_t lam_R, real_t mu_R, real_t rho_R,
   DenseMatrix& qGodLocal,
   DenseMatrix& qGodNeighbor)
{
   // Acoustic-input guard.  SeisSol handles mu=0 (acoustic) with a
   // separate branch in getTransposedGodunovState
   // (ElasticSetup.h:92-99 + :119-127); SAFS does not have acoustic
   // regions so we abort.
   {
      // mu_eps = 1e-12 Pa is "essentially zero" for solid earth
      // materials (typical mu ~ 1e9..1e11 Pa, so 21+ orders of
      // margin).  SAFS uses SI units throughout (see
      // miniapps/seas/CLAUDE.md "All physical units SI").
      const real_t mu_eps = 1e-12;
      MFEM_VERIFY(mu_L > mu_eps && mu_R > mu_eps,
                  "BimaterialFlux::BuildGodunovStateFaceLocal: acoustic "
                  "input (mu_L=" << mu_L << ", mu_R=" << mu_R
                  << ") not supported; SAFS does not have acoustic "
                  "regions.  Add the acoustic branch from "
                  "SeisSol/src/Equations/elastic/Model/ElasticSetup.h:92-139 "
                  "to enable.");
   }

   MFEM_VERIFY(rho_L > 0.0,
               "BimaterialFlux::BuildGodunovStateFaceLocal: rho_L must "
               "be positive, got " << rho_L);
   MFEM_VERIFY(rho_R > 0.0,
               "BimaterialFlux::BuildGodunovStateFaceLocal: rho_R must "
               "be positive, got " << rho_R);
   MFEM_VERIFY(lam_L + 2.0 * mu_L > 0.0,
               "BimaterialFlux::BuildGodunovStateFaceLocal: "
               "lam_L+2*mu_L must be positive, got "
               << lam_L + 2.0 * mu_L);
   MFEM_VERIFY(lam_R + 2.0 * mu_R > 0.0,
               "BimaterialFlux::BuildGodunovStateFaceLocal: "
               "lam_R+2*mu_R must be positive, got "
               << lam_R + 2.0 * mu_R);

   const real_t cp_L = std::sqrt((lam_L + 2.0 * mu_L) / rho_L);
   const real_t cs_L = std::sqrt(mu_L / rho_L);
   const real_t cp_R = std::sqrt((lam_R + 2.0 * mu_R) / rho_R);
   const real_t cs_R = std::sqrt(mu_R / rho_R);
   const real_t lp_L = lam_L + 2.0 * mu_L;
   const real_t lp_R = lam_R + 2.0 * mu_R;

   // ----------------------------------------------------------------
   // Build matR (9×9) in SeisSol's exact sign convention
   // (ElasticSetup.h:102-139, non-acoustic elastic branch).
   //
   // Index mapping is identical between SeisSol and MFEM's wave_state
   // ordering: (σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz, v_x, v_y, v_z)
   // = (SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ).
   //
   // SeisSol uses POSITIVE stress entries for the "LOCAL outgoing"
   // cols 0..2 — these are eigenvectors of MFEM's A_x with eigenvalue
   // −cp_L / −cs_L (not +cp_L) under MFEM's convention.  The chi-
   // projection-based formula compensates so that the final flux is
   // the standard upwind A^+ Q_L + A^- Q_R.  See file header.
   // ----------------------------------------------------------------
   DenseMatrix matR(NUM_STATE, NUM_STATE);
   matR = 0.0;

   // Cols 0..2: LOCAL material columns (SeisSol ElasticSetup.h:102-111).

   // Col 0: P-wave at LOCAL material.
   matR(SXX, 0) = lam_L + 2.0 * mu_L;       //  σ_xx
   matR(SYY, 0) = lam_L;                    //  σ_yy
   matR(SZZ, 0) = lam_L;                    //  σ_zz
   matR(VX,  0) = std::sqrt((lam_L + 2.0 * mu_L) / rho_L);  // v_x = +cp_L

   // Col 1: S-wave at LOCAL material, σ_xy / v_y polarisation.
   matR(SXY, 1) = mu_L;                     //  σ_xy
   matR(VY,  1) = std::sqrt(mu_L / rho_L);  //  v_y = +cs_L

   // Col 2: S-wave at LOCAL material, σ_xz / v_z polarisation.
   matR(SXZ, 2) = mu_L;                     //  σ_xz
   matR(VZ,  2) = std::sqrt(mu_L / rho_L);  //  v_z = +cs_L

   // Cols 3..5: zero-eigenvalue modes (passive); SeisSol uses the
   // (lam+2mu) scaling for matR conditioning (ElasticSetup.h:115-117).
   // We follow SeisSol's permutation:
   //   col 3 puts (lam_L + 2*mu_L) into row 4 (σ_yz),
   //   col 4 puts (lam_L + 2*mu_L) into row 1 (σ_yy),
   //   col 5 puts (lam_L + 2*mu_L) into row 2 (σ_zz).
   matR(SYZ, 3) = lp_L;                    // σ_yz passive mode
   matR(SYY, 4) = lp_L;                    // σ_yy passive mode
   matR(SZZ, 5) = lp_L;                    // σ_zz passive mode

   // Cols 6..8: NEIGHBOUR material columns (SeisSol ElasticSetup.h:128-139).
   // These have NEGATIVE wave-speed entries on the velocity rows
   // (i.e., −cs_R, −cs_R, −cp_R) — the "incoming from neighbour"
   // characteristic family.

   // Col 6: S-wave at NEIGHBOUR material, σ_xz / v_z polarisation.
   matR(SXZ, 6) =  mu_R;
   matR(VZ,  6) = -std::sqrt(mu_R / rho_R);

   // Col 7: S-wave at NEIGHBOUR material, σ_xy / v_y polarisation.
   matR(SXY, 7) =  mu_R;
   matR(VY,  7) = -std::sqrt(mu_R / rho_R);

   // Col 8: P-wave at NEIGHBOUR material.
   matR(SXX, 8) =  lam_R + 2.0 * mu_R;
   matR(SYY, 8) =  lam_R;
   matR(SZZ, 8) =  lam_R;
   matR(VX,  8) = -std::sqrt((lam_R + 2.0 * mu_R) / rho_R);

   // ----------------------------------------------------------------
   // Invert matR via mfem::DenseMatrixInverse.
   // ----------------------------------------------------------------
   DenseMatrix matR_copy(matR);
   DenseMatrixInverse matR_inv_solver(matR_copy);
   DenseMatrix matR_inv(NUM_STATE, NUM_STATE);
   matR_inv_solver.GetInverseMatrix(matR_inv);

   // ----------------------------------------------------------------
   // chi · matR^{-1}: zero rows 3..8 of matR_inv (chi selects rows
   // 0..2 via left-multiplication by the diagonal selector
   // diag(1,1,1,0,0,0,0,0,0)).
   // ----------------------------------------------------------------
   DenseMatrix chi_mRinv(matR_inv);
   for (int r = 3; r < NUM_STATE; ++r)
   {
      for (int c = 0; c < NUM_STATE; ++c)
      {
         chi_mRinv(r, c) = 0.0;
      }
   }

   // ----------------------------------------------------------------
   // godunov = matR · (chi · matR^{-1}).  In MFEM's sign convention
   // this projects onto the NEGATIVE-eigenvalue subspace of A
   // (because SeisSol's "LOCAL outgoing" eigenvectors are actually
   // MFEM's "LOCAL incoming"; see file header).
   //
   // qGodNeighbor = godunov          (multiplier of Q_R)
   // qGodLocal    = I - godunov      (multiplier of Q_L)
   //
   // Final flux: F^self = A · Q* = A · (qGodLocal Q_L + qGodNeighbor Q_R)
   //                            = A · ((I - godunov) Q_L + godunov Q_R)
   //                            = (A - A^-) Q_L + A^- Q_R
   //                            = A^+ Q_L + A^- Q_R    (standard upwind).
   // ----------------------------------------------------------------
   DenseMatrix godunov(NUM_STATE, NUM_STATE);
   mfem::Mult(matR, chi_mRinv, godunov);

   qGodNeighbor.SetSize(NUM_STATE, NUM_STATE);
   qGodLocal.SetSize(NUM_STATE, NUM_STATE);
   for (int i = 0; i < NUM_STATE; ++i)
   {
      for (int j = 0; j < NUM_STATE; ++j)
      {
         qGodNeighbor(i, j) = godunov(i, j);
         qGodLocal(i, j) = (i == j ? 1.0 : 0.0) - godunov(i, j);
      }
   }
}

// ---------------------------------------------------------------------------
// BuildPerFaceFluxMatricesGlobal — plan §R.1 step 3.
//
//   fluxLocal    = T · A_self_FL · qGodLocal    · T^{-1}
//   fluxNeighbor = T · A_self_FL · qGodNeighbor · T^{-1}
//
// At runtime, F_self_global = fluxLocal · Q_self_global
//                           + fluxNeighbor · Q_nbr_global.
// ---------------------------------------------------------------------------
/* static */ void BimaterialFlux::BuildPerFaceFluxMatricesGlobal(
   const real_t* nor,
   const GodunovFlux& flux_self,
   const GodunovFlux& flux_nbr,
   DenseMatrix& fluxLocal,
   DenseMatrix& fluxNeighbor)
{
   // Precondition: `nor` is unit.  GodunovFlux::BuildFrame and the
   // subsequent T/Tinv assume |nor|=1; an un-normalised normal
   // silently produces non-orthonormal rotation matrices (T·Tinv ≠ I)
   // and corrupts the flux composition without any visible error.
   // The R.2 face-loop caller MUST normalise CalcOrtho's output
   // before forwarding (CalcOrtho returns the un-normalised normal
   // with magnitude |J_F|).  Hard check here for defence in depth.
   const real_t n2 = nor[0]*nor[0] + nor[1]*nor[1] + nor[2]*nor[2];
   MFEM_VERIFY(std::abs(n2 - 1.0) < 1e-10,
               "BimaterialFlux::BuildPerFaceFluxMatricesGlobal: `nor` "
               "must be a unit vector (got |nor|^2 = " << n2 << ").");

   // 1. Build orthonormal frame and rotation matrices using
   //    GodunovFlux's existing helpers.
   real_t t1[3], t2[3];
   GodunovFlux::BuildFrame(nor, t1, t2);

   DenseMatrix T(NUM_STATE, NUM_STATE);
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   GodunovFlux::BuildRotation(nor, t1, t2, T);
   GodunovFlux::BuildRotationInverse(nor, t1, t2, Tinv);

   // 2. Build qGodLocal / qGodNeighbor in the face-local frame.
   DenseMatrix qGodL_FL, qGodN_FL;
   BuildGodunovStateFaceLocal(flux_self.GetLambda(), flux_self.GetMu(),
                              flux_self.GetRho(),
                              flux_nbr.GetLambda(), flux_nbr.GetMu(),
                              flux_nbr.GetRho(),
                              qGodL_FL, qGodN_FL);

   // 3. A_self_FL is the LOCAL-side x-direction Jacobian (face-local
   //    frame's "x" = face normal).
   const DenseMatrix& A_self_FL = flux_self.GetAx();

   // 4. Compose the per-side flux matrices in the global frame.
   DenseMatrix tmp1(NUM_STATE, NUM_STATE);
   DenseMatrix tmp2(NUM_STATE, NUM_STATE);

   // fluxLocal = T · A · qGodLocal · T^{-1}
   mfem::Mult(A_self_FL, qGodL_FL, tmp1);
   mfem::Mult(T, tmp1, tmp2);
   fluxLocal.SetSize(NUM_STATE, NUM_STATE);
   mfem::Mult(tmp2, Tinv, fluxLocal);

   // fluxNeighbor = T · A · qGodNeighbor · T^{-1}
   mfem::Mult(A_self_FL, qGodN_FL, tmp1);
   mfem::Mult(T, tmp1, tmp2);
   fluxNeighbor.SetSize(NUM_STATE, NUM_STATE);
   mfem::Mult(tmp2, Tinv, fluxNeighbor);
}

// ---------------------------------------------------------------------------
// BuildPerFaceCentralMatricesGlobal — PLAN_mixed_flux_hetero_riemann.md Phase 1.
//
// Bi-material CENTRAL flux for fault-adjacent faces.  Produces the two
// per-face matrices ½·A_self and ½·A_nbr (global frame) so that
//   F* = ½·A_self·Q_self + ½·A_nbr·Q_nbr
// is single-valued (deposited identically to both sides — see the header
// doc + PLAN BUG-3).  Mirrors BuildPerFaceFluxMatricesGlobal's rotation +
// scratch-matrix machinery, but composes ½·T·(AxPlus+AxMinus)·T^{-1} per
// side (the FULL face-normal Jacobian, exactly as GodunovFlux::Central)
// instead of the upwind T·A·qGod·T^{-1}.  There is NO bi-material
// characteristic projection here: the central flux is the average of the
// two physical fluxes and carries no shared Riemann state.
// ---------------------------------------------------------------------------
/* static */ void BimaterialFlux::BuildPerFaceCentralMatricesGlobal(
   const real_t* nor,
   const GodunovFlux& flux_self,
   const GodunovFlux& flux_nbr,
   DenseMatrix& centralSelf,
   DenseMatrix& centralNbr)
{
   // Precondition: `nor` is unit (same contract + rationale as
   // BuildPerFaceFluxMatricesGlobal — an un-normalised normal silently
   // produces non-orthonormal T/Tinv and corrupts the flux composition).
   const real_t n2 = nor[0]*nor[0] + nor[1]*nor[1] + nor[2]*nor[2];
   MFEM_VERIFY(std::abs(n2 - 1.0) < 1e-10,
               "BimaterialFlux::BuildPerFaceCentralMatricesGlobal: `nor` "
               "must be a unit vector (got |nor|^2 = " << n2 << ").");

   // Acoustic-input guard (PLAN Phase 1 Edge Case / BUG-9).  Unlike the
   // upwind builder, the central path does NOT route through
   // BuildGodunovStateFaceLocal, so it inherits no acoustic check — add one
   // explicitly with the SAME mu_eps = 1e-12 contract used there (SAFS has
   // no acoustic regions; the acoustic branch is a documented follow-up).
   const real_t mu_eps = 1e-12;
   MFEM_VERIFY(flux_self.GetMu() > mu_eps && flux_nbr.GetMu() > mu_eps,
               "BimaterialFlux::BuildPerFaceCentralMatricesGlobal: acoustic "
               "input (mu_self=" << flux_self.GetMu() << ", mu_nbr="
               << flux_nbr.GetMu() << ") not supported; SAFS does not have "
               "acoustic regions.  Add the acoustic branch from "
               "SeisSol/src/Equations/elastic/Model/ElasticSetup.h:92-139 "
               "to enable.");

   // 1. Orthonormal frame + rotation matrices — the SAME helpers (and hence
   //    the SAME frame convention) used by BuildPerFaceFluxMatricesGlobal
   //    and GodunovFlux::Central, so the central and Godunov per-face
   //    contributions stay mutually consistent.
   real_t t1[3], t2[3];
   GodunovFlux::BuildFrame(nor, t1, t2);

   DenseMatrix T(NUM_STATE, NUM_STATE);
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   GodunovFlux::BuildRotation(nor, t1, t2, T);
   GodunovFlux::BuildRotationInverse(nor, t1, t2, Tinv);

   // 2. Per side, build the GLOBAL-frame ½·A_side.
   //    A_side_facelocal = AxPlus_side + AxMinus_side  is the FULL
   //    face-normal Jacobian — IDENTICAL to GodunovFlux::Central's
   //    (Ax_plus_ + Ax_minus_) construction (godunov_flux.cpp), NOT GetAx()
   //    (built by a separate path; would only agree to the ~1e-9
   //    sum-of-splits floor).  centralSide = ½ · T · A_side_facelocal · T^{-1}.
   //    Two distinct scratch matrices (Asum, tmp); SetSize ONLY the output
   //    (mfem::Mult requires a separate, correctly-sized output — never
   //    alias the output with a scratch).  The ½ is applied to the output
   //    via `*= 0.5` AFTER the rotations.
   DenseMatrix Asum(NUM_STATE, NUM_STATE);
   DenseMatrix tmp(NUM_STATE, NUM_STATE);

   // --- self side: centralSelf = ½ · T · (AxPlus_self + AxMinus_self) · T^{-1}
   Asum  = flux_self.GetAxPlus();
   Asum += flux_self.GetAxMinus();
   mfem::Mult(T, Asum, tmp);
   centralSelf.SetSize(NUM_STATE, NUM_STATE);
   mfem::Mult(tmp, Tinv, centralSelf);
   centralSelf *= 0.5;

   // --- neighbour side: centralNbr = ½ · T · (AxPlus_nbr + AxMinus_nbr) · T^{-1}
   Asum  = flux_nbr.GetAxPlus();
   Asum += flux_nbr.GetAxMinus();
   mfem::Mult(T, Asum, tmp);
   centralNbr.SetSize(NUM_STATE, NUM_STATE);
   mfem::Mult(tmp, Tinv, centralNbr);
   centralNbr *= 0.5;
}

// ---------------------------------------------------------------------------
// ApplyPerFaceFlux — plan §R.1 step 4.
// ---------------------------------------------------------------------------
/* static */ void BimaterialFlux::ApplyPerFaceFlux(
   const DenseMatrix& fluxLocal,
   const DenseMatrix& fluxNeighbor,
   const real_t* Q_self,
   const real_t* Q_nbr,
   real_t* F_h_self)
{
   // Hard runtime check (not MFEM_ASSERT): R.2's per-(face, side)
   // storage leaves shared-face side=1 default-constructed; a
   // mis-indexed lookup would silently corrupt output in release
   // builds.  Per-QP cost is one int comparison vs the 162 multiplies
   // in this function — negligible.
   MFEM_VERIFY(fluxLocal.Height() == NUM_STATE
               && fluxLocal.Width() == NUM_STATE
               && fluxNeighbor.Height() == NUM_STATE
               && fluxNeighbor.Width() == NUM_STATE,
               "BimaterialFlux::ApplyPerFaceFlux: input matrices must "
               "both be NUM_STATE × NUM_STATE.");
   for (int i = 0; i < NUM_STATE; ++i)
   {
      real_t s = 0.0;
      for (int j = 0; j < NUM_STATE; ++j)
      {
         s += fluxLocal(i, j) * Q_self[j]
              + fluxNeighbor(i, j) * Q_nbr[j];
      }
      F_h_self[i] = s;
   }
}

// ---------------------------------------------------------------------------
// IsStrongContrast — central-flux corridor guard predicate (Part A, R-007).
// Returns true iff tol >= 0 and the relative impedance jump exceeds tol,
// measured as max(|Zp jump|/maxZp, |Zs jump|/maxZs).  tol < 0 disables; a
// non-positive max impedance (acoustic / unset) returns false.
// ---------------------------------------------------------------------------
/* static */ real_t BimaterialFlux::ContrastValue(const GodunovFlux& a,
                                                  const GodunovFlux& b)
{
   const real_t zp_a = a.GetZp(), zp_b = b.GetZp();
   const real_t zs_a = a.GetZs(), zs_b = b.GetZs();
   const real_t max_zp = std::max(zp_a, zp_b);
   const real_t max_zs = std::max(zs_a, zs_b);
   if (max_zp <= 0.0 || max_zs <= 0.0) { return -1.0; }   // acoustic / unset
   const real_t c_zp = std::abs(zp_a - zp_b) / max_zp;
   const real_t c_zs = std::abs(zs_a - zs_b) / max_zs;
   return std::max(c_zp, c_zs);
}

/* static */ bool BimaterialFlux::IsStrongContrast(const GodunovFlux& a,
                                                   const GodunovFlux& b,
                                                   real_t tol)
{
   if (tol < 0.0) { return false; }
   return ContrastValue(a, b) > tol;   // ContrastValue == -1 (acoustic) => false
}

}  // namespace seas
}  // namespace mfem
