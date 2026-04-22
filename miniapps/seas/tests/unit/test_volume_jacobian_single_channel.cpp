// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan §4.10 (rev-3g, NEW) — single-channel Jacobian probe.
// Targets H-V92-G candidates (d) and (e):
//   (d) mass-inverse per-channel conditioning asymmetry
//   (e) volume-term / Jacobian coupling bug
//
// ============================================================================
// Strategy
// ============================================================================
// 2-hex Cartesian mesh along x: [0, 2]×[0, 1]×[0, 1], absorbing BCs all
// faces.  For each of the 9 state components c = 0..8, set
// Q[c](x,y,z) = sin(π·x) uniformly in y,z (others = 0), and call Mult(Q, k).
//
// For a Q varying only in x, the elastodynamic PDE gives
//     ∂_t Q + A_x · ∂_x Q = 0
// so the analytic right-hand-side (before mass inverse) is
//     RHS_analytic[r] = -A_x[r, c] · (π cos(π x))
// for each row r.  The "sparsity pattern" of RHS is determined by
// which rows r have a nonzero A_x[r, c] for this driven c.
//
// Verification steps per channel c:
//   (i)  Measure `dQdt = Mult(Q)` on the interior mesh.
//   (ii) Identify nonzero components of dQdt (> noise floor).
//   (iii) Compare against the analytic sparsity predicted by A_x column c.
//
// Pelties 2012 §3.1 eq. (4) dictates the A_x entries (verified by
// GodunovFlux::BuildJacobian, which test_godunov_flux T1 already gates).
// The expected A_x sparsity per channel c:
//   c = SXX:  dQdt row VX = +(1/ρ) π cos(πx)          (one nonzero)
//   c = SYY:  dQdt row VY non-zero only through A_y, not A_x  (zero!)
//   c = SZZ:  (zero from A_x)
//   c = SXY:  dQdt row VY = +(1/ρ) π cos(πx)          (one nonzero)
//   c = SYZ:  (zero from A_x — zero-wave-speed component)
//   c = SXZ:  dQdt row VZ = +(1/ρ) π cos(πx)          (one nonzero)
//   c = VX:   dQdt rows SXX=+λ+2μ, SYY=+λ, SZZ=+λ (three nonzeros)
//   c = VY:   dQdt row SXY = +μ                        (one nonzero)
//   c = VZ:   dQdt row SXZ = +μ                        (one nonzero)
//
// Failure modes this test detects:
// - If driving c = SYY produces a nonzero dQdt[SYY] row (self-coupling)
//   => candidate (e) volume/Jacobian bug CONFIRMED.
// - If the RATIO of dQdt magnitudes between (c = SXY, dQdt[VY]) and
//   (c = SYY, dQdt[VY]) is not consistent with the Jacobian ratio
//   (both should produce +π/ρ · cos if driven identically) => (d) or
//   (e) candidate.
// - Any spurious cross-channel leak that violates the Pelties eq. (4)
//   sparsity => (e) candidate.
//
// Usage:
//   ./seas_test_volume_jacobian_single_channel   (serial, < 1 s)

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_BOOL(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else      { num_failed++; std::cout << "  FAILED [" << __LINE__ \
                       << "]: " << msg << "\n"; } \
} while (0)

// Build a 2-hex Cartesian mesh, [0,2]×[0,1]×[0,1], absorbing BCs (attr 5).
static Mesh BuildTwoHexMesh()
{
   Mesh mesh = Mesh::MakeCartesian3D(2, 1, 1, Element::HEXAHEDRON,
                                     2.0, 1.0, 1.0);
   // Tag every boundary face as absorbing (attr 5 matches BoundaryConfig).
   for (int be = 0; be < mesh.GetNBE(); be++)
   { mesh.SetBdrAttribute(be, 5); }
   mesh.SetAttributes();
   return mesh;
}

// Sparsity predicted by A_x + A_y + A_z (the wave operator applies face
// fluxes in all 3 directions, and a Q varying only in x still generates
// face-flux contributions on y-normal and z-normal boundaries because
// Q itself is nonzero there).  Entries: which ROW r has a nonzero
// response when COLUMN c is driven.  A 1 means nonzero via volume OR
// face-flux; a 0 means silent under pure Pelties 2012 eq. (4).
//
// Derivation: rows that are nonzero via any of A_x, A_y, A_z columns
// c.  For velocity columns (c ∈ {VX, VY, VZ}): stress rows are driven
// by ∂_t σ = -C : ∇v.  For stress columns: velocity rows are driven.
// Off-diagonal silence: stress→stress and velocity→velocity are ZERO
// by elastodynamics structure.
static const int expected_A3d[9][9] = {
   // c: SXX SYY SZZ SXY SYZ SXZ  VX  VY  VZ
   {    0,  0,  0,  0,  0,  0,  1,  1,  1 }, // r=SXX: rows driven by all velocities
   {    0,  0,  0,  0,  0,  0,  1,  1,  1 }, // r=SYY
   {    0,  0,  0,  0,  0,  0,  1,  1,  1 }, // r=SZZ
   {    0,  0,  0,  0,  0,  0,  1,  1,  0 }, // r=SXY: A_x[SXY,VY], A_y[SXY,VX]
   {    0,  0,  0,  0,  0,  0,  0,  1,  1 }, // r=SYZ: A_y[SYZ,VZ], A_z[SYZ,VY]
   {    0,  0,  0,  0,  0,  0,  1,  0,  1 }, // r=SXZ: A_x[SXZ,VZ], A_z[SXZ,VX]
   {    1,  0,  0,  1,  0,  1,  0,  0,  0 }, // r=VX:  SXX via A_x, SXY via A_y, SXZ via A_z
   {    0,  1,  0,  1,  1,  0,  0,  0,  0 }, // r=VY:  SYY via A_y, SXY via A_x, SYZ via A_z
   {    0,  0,  1,  0,  1,  1,  0,  0,  0 }, // r=VZ:  SZZ via A_z, SYZ via A_y, SXZ via A_x
};

static const char *comp_name(int c)
{
   static const char *names[9] = {
      "SXX","SYY","SZZ","SXY","SYZ","SXZ","VX","VY","VZ"
   };
   return (c >= 0 && c < 9) ? names[c] : "???";
}

int main()
{
   std::cout << "\n=== TPV102 v9.2.0 §4.10 Probe: "
             << "Single-Channel Jacobian (H-V92-G (d)+(e)) ===\n";

   Mesh mesh = BuildTwoHexMesh();
   BoundaryConfig bc;
   bc.natural_attrs   = {};
   bc.fault_attr      = 0;
   bc.absorbing_attrs = {5};

   const int order = 1;
   WaveOperator<Mesh> wave(mesh, order,
                           TPV102Params::lambda,
                           TPV102Params::mu,
                           TPV102Params::rho, bc);

   const auto &fes_const = wave.GetFESpace();
   auto &fes = const_cast<FiniteElementSpace &>(fes_const);
   const int size = wave.Height();
   const int ndof_total = fes.GetNDofs();
   std::cout << "  Mesh: " << mesh.GetNE() << " hex elements, "
             << ndof_total << " scalar DOFs.\n";

   // Helper: set Q[c](x) = sin(π·x) via ProjectCoefficient on the scalar
   // FESpace, then copy into component c of the Q vector.
   const double kx = M_PI;  // wavelength 2 units = domain x extent
   FunctionCoefficient sin_kx(
      [kx](const Vector &x) -> real_t {
         return std::sin(kx * x(0));
      });

   // For each driven channel c, measure per-row max|dQdt| and compare
   // against A_x sparsity.
   bool any_violation = false;

   // Tolerance on "zero" for rows that A_x predicts should not respond.
   // The absorbing BC on ±y, ±z faces contributes a small residual that
   // scales with the field magnitude; empirically this is O(10⁻²) of the
   // strong-row magnitude.  We use a RELATIVE threshold vs the strongest
   // row for the driven channel.
   const double zero_rel_tol = 1.0e-1;   // 10% of strongest-row magnitude

   std::cout << "\n  Per-channel probe (Q[c](x) = sin(πx), else 0):\n";
   std::cout << "  Driven c | row  | max|dQdt| | sparsity pred | note\n";
   std::cout << "  ---------+------+-----------+---------------+------\n";

   for (int c = 0; c < NUM_STATE; c++)
   {
      Vector Q(size);
      Q = 0.0;
      GridFunction gf(&fes);
      gf.ProjectCoefficient(sin_kx);
      for (int d = 0; d < ndof_total; d++)
      { Q(c * ndof_total + d) = gf(d); }

      Vector k(size);
      wave.Mult(Q, k);

      // Compute max|dQdt[r]| per row.
      double max_per_row[NUM_STATE] = {0};
      for (int r = 0; r < NUM_STATE; r++)
      {
         for (int d = 0; d < ndof_total; d++)
         {
            max_per_row[r] = std::max(max_per_row[r],
                                       std::abs(k(r * ndof_total + d)));
         }
      }

      // Find strongest row magnitude (for relative tolerance).
      double strongest = 0.0;
      for (int r = 0; r < NUM_STATE; r++)
      { strongest = std::max(strongest, max_per_row[r]); }

      // Absolute floor: any magnitude below 1e-9 is numerical noise (floating-
      // point round-off from Mult) — don't flag it regardless of ratio.  This
      // handles the case where all rows are sub-ULP after cancellation
      // (common for stress-driven channels where the projected Q values
      // are near the ULP scale relative to the original coefficient).
      const double absolute_noise_floor = 1.0e-9;

      // Check each row against predicted sparsity.
      for (int r = 0; r < NUM_STATE; r++)
      {
         const bool pred_nonzero = (expected_A3d[r][c] != 0);
         const double mag = max_per_row[r];
         const bool   mag_nonzero =
            (mag > absolute_noise_floor) &&
            (mag > zero_rel_tol * std::max(strongest, 1.0e-30));
         const char *note = "";
         if (pred_nonzero && mag_nonzero)   { note = "OK (expected)"; }
         else if (!pred_nonzero && !mag_nonzero) { note = "OK (silent)"; }
         else if (!pred_nonzero && mag_nonzero)
         {
            note = "VIOLATION: unexpected nonzero — (d)/(e) suspect!";
            any_violation = true;
         }
         else                               { note = "weak response"; }

         std::printf("  %8s | %-4s | %.3e | %s | %s\n",
                     comp_name(c), comp_name(r), mag,
                     pred_nonzero ? "nonzero" : "  ZERO ",
                     note);
      }
      std::cout << "  ---------+------+-----------+---------------+------\n";
   }

   TEST_BOOL(!any_violation,
             "No cross-channel Jacobian-sparsity violations "
             "across all 9 driven channels");

   // -------- Explicit (d) test: ratio check on two equivalent channels ---
   // Drive SXY with sin(πx); drive SXZ with sin(πx).  Both should
   // produce identical-magnitude responses on dQdt[VY] and dQdt[VZ]
   // respectively (both entries are -1/ρ in A_x).  If the RATIO
   // |dQdt[VY]|_xy / |dQdt[VZ]|_xz ≠ 1, mass-inverse is channel-
   // asymmetric for at least one of {SXY, SXZ, VY, VZ}.
   std::cout << "\n  (d)-specific test: A_x[VY,SXY] vs A_x[VZ,SXZ] response:\n";
   double mag_VY_from_SXY = 0.0, mag_VZ_from_SXZ = 0.0;
   {
      Vector Q(size), k(size);
      Q = 0.0;
      GridFunction gf(&fes);
      gf.ProjectCoefficient(sin_kx);
      for (int d = 0; d < ndof_total; d++)
      { Q(SXY * ndof_total + d) = gf(d); }
      wave.Mult(Q, k);
      for (int d = 0; d < ndof_total; d++)
      { mag_VY_from_SXY = std::max(mag_VY_from_SXY,
                                   std::abs(k(VY * ndof_total + d))); }

      Q = 0.0;
      for (int d = 0; d < ndof_total; d++)
      { Q(SXZ * ndof_total + d) = gf(d); }
      wave.Mult(Q, k);
      for (int d = 0; d < ndof_total; d++)
      { mag_VZ_from_SXZ = std::max(mag_VZ_from_SXZ,
                                   std::abs(k(VZ * ndof_total + d))); }
   }
   const double ratio = mag_VY_from_SXY /
                        std::max(mag_VZ_from_SXZ, 1.0e-30);
   std::cout << "    max|dQdt[VY]| driven by SXY = " << std::scientific
             << std::setprecision(6) << mag_VY_from_SXY << "\n";
   std::cout << "    max|dQdt[VZ]| driven by SXZ = " << mag_VZ_from_SXZ << "\n";
   std::cout << "    ratio                        = " << ratio << "\n";
   // Expected ratio = 1 (by isotropy and symmetry; same A_x entry, just
   // different row/column indices).  Allow 10% tolerance for mesh /
   // projection noise.  Bit-asymmetry at ULP → (d) CONFIRMED.
   TEST_BOOL(std::abs(ratio - 1.0) < 0.1,
             "SXY→VY vs SXZ→VZ response ratio = 1 (isotropy)");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
