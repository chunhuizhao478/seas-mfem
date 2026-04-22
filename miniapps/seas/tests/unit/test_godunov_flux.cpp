// Phase 1 unit tests for GodunovFlux (Tests 1-8 from plan Section 3.1.5).

#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) < t_) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", expected " << e_ << ", diff " << std::abs(v_-e_) << ")\n"; } \
} while(0)

// ===== Test 1: Jacobian A,B,C symmetry and nonzero entries =====
void TestJacobianSymmetryABC()
{
   std::cout << "Test 1: TestJacobianSymmetryABC\n";
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   const DenseMatrix &Ax = flux.GetAx();

   // Check known entries: A[SXX][VX] = -(lambda+2*mu)
   TEST_NEAR(Ax(SXX, VX), -(lambda + 2.0*mu), 1e-6,
             "A[SXX][VX] = -(lambda+2mu)");

   // A[VX][SXX] = -1/rho
   TEST_NEAR(Ax(VX, SXX), -1.0/rho, 1e-20,
             "A[VX][SXX] = -1/rho");

   // A[SXY][VY] = -mu
   TEST_NEAR(Ax(SXY, VY), -mu, 1e-6, "A[SXY][VY] = -mu");

   // A[VY][SXY] = -1/rho
   TEST_NEAR(Ax(VY, SXY), -1.0/rho, 1e-20, "A[VY][SXY] = -1/rho");

   // B (y-direction): B[SYY][VY] = -(lambda+2*mu), B[VY][SYY] = -1/rho
   DenseMatrix By(NUM_STATE, NUM_STATE);
   flux.BuildJacobian(1, By);
   TEST_NEAR(By(SYY, VY), -(lambda + 2.0*mu), 1e-6, "B[SYY][VY] = -(lambda+2mu)");
   TEST_NEAR(By(VY, SYY), -1.0/rho, 1e-20,          "B[VY][SYY] = -1/rho");
   TEST_NEAR(By(SXY, VX), -mu, 1e-6,                 "B[SXY][VX] = -mu");
   TEST_NEAR(By(VX, SXY), -1.0/rho, 1e-20,           "B[VX][SXY] = -1/rho");

   // C (z-direction): C[SZZ][VZ] = -(lambda+2*mu), C[VZ][SZZ] = -1/rho
   DenseMatrix Cz(NUM_STATE, NUM_STATE);
   flux.BuildJacobian(2, Cz);
   TEST_NEAR(Cz(SZZ, VZ), -(lambda + 2.0*mu), 1e-6, "C[SZZ][VZ] = -(lambda+2mu)");
   TEST_NEAR(Cz(VZ, SZZ), -1.0/rho, 1e-20,          "C[VZ][SZZ] = -1/rho");
   TEST_NEAR(Cz(SXZ, VX), -mu, 1e-6,                 "C[SXZ][VX] = -mu");
   TEST_NEAR(Cz(VX, SXZ), -1.0/rho, 1e-20,           "C[VX][SXZ] = -1/rho");
}

// ===== Test 2: Eigenvalues via split flux =====
void TestEigenvalues()
{
   std::cout << "Test 2: TestEigenvalues\n";
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   real_t cp = flux.GetCp();
   real_t cs = flux.GetCs();

   // Verify: A = A^+ + A^-
   const DenseMatrix &Ax = flux.GetAx();
   const DenseMatrix &Ap = flux.GetAxPlus();
   const DenseMatrix &Am = flux.GetAxMinus();

   real_t max_err = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
      for (int j = 0; j < NUM_STATE; j++)
         max_err = std::max(max_err, std::abs(Ax(i,j) - Ap(i,j) - Am(i,j)));

   TEST_ASSERT(max_err < 1e-10,
               "A = A^+ + A^- (max error " + std::to_string(max_err) + ")");

   // Check eigenvalues: apply A^+ and A^- to known eigenvectors.
   // A P-wave eigenvector for +cp: e = [-Zp, -lambda/cp, -lambda/cp, 0, 0, 0, cp, 0, 0]
   // should give: A * e = cp * e
   real_t Zp = flux.GetZp();
   real_t lp = lambda + 2.0*mu;
   real_t e[NUM_STATE] = {-lp, -lambda, -lambda, 0, 0, 0, cp, 0, 0};
   real_t Ae[NUM_STATE];
   Ax.Mult(e, Ae);

   real_t eigval_err = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
      eigval_err = std::max(eigval_err, std::abs(Ae[i] - cp * e[i]));

   TEST_ASSERT(eigval_err < 1e-6,
               "P-wave eigenvalue cp verified (error " + std::to_string(eigval_err) + ")");

   // S-wave eigenvector for +cs (y-polarized): [0,0,0,-mu,0,0,0,cs,0]
   // Values are O(mu*cs) ≈ 1e14, so use relative tolerance.
   real_t e_s[NUM_STATE] = {0, 0, 0, -mu, 0, 0, 0, cs, 0};
   real_t Ae_s[NUM_STATE];
   Ax.Mult(e_s, Ae_s);
   real_t s_err = 0.0, s_scale = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
   {
      s_err = std::max(s_err, std::abs(Ae_s[i] - cs * e_s[i]));
      s_scale = std::max(s_scale, std::abs(cs * e_s[i]));
   }
   TEST_ASSERT(s_err / s_scale < 1e-10,
               "S-wave eigenvalue cs verified (rel error " + std::to_string(s_err/s_scale) + ")");

   // Leftgoing P-wave eigenvector for -cp: [+lp, +lambda, +lambda, 0,0,0, cp, 0, 0]
   real_t e_left[NUM_STATE] = {lp, lambda, lambda, 0, 0, 0, cp, 0, 0};
   real_t Ae_left[NUM_STATE];
   Ax.Mult(e_left, Ae_left);
   real_t left_err = 0.0, left_scale = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
   {
      left_err = std::max(left_err, std::abs(Ae_left[i] - (-cp) * e_left[i]));
      left_scale = std::max(left_scale, std::abs(cp * e_left[i]));
   }
   TEST_ASSERT(left_err / left_scale < 1e-10,
               "Leftgoing P-wave eigenvalue -cp verified (rel error " + std::to_string(left_err/left_scale) + ")");
}

// ===== Test 3: Rotation T * T^{-1} = Identity =====
void TestRotationOrthogonality()
{
   std::cout << "Test 3: TestRotationOrthogonality\n";

   // Test with several normals
   real_t normals[][3] = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
      {1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)}
   };

   for (int n = 0; n < 4; n++)
   {
      real_t t1[3], t2[3];
      GodunovFlux::BuildFrame(normals[n], t1, t2);

      DenseMatrix T(NUM_STATE, NUM_STATE);
      DenseMatrix Tinv(NUM_STATE, NUM_STATE);
      GodunovFlux::BuildRotation(normals[n], t1, t2, T);
      GodunovFlux::BuildRotationInverse(normals[n], t1, t2, Tinv);

      // Product T * T^{-1} should be identity
      DenseMatrix prod(NUM_STATE, NUM_STATE);
      mfem::Mult(T, Tinv, prod);

      real_t max_err = 0.0;
      for (int i = 0; i < NUM_STATE; i++)
         for (int j = 0; j < NUM_STATE; j++)
         {
            real_t expected = (i == j) ? 1.0 : 0.0;
            max_err = std::max(max_err, std::abs(prod(i,j) - expected));
         }

      TEST_ASSERT(max_err < 1e-14,
                  "T * T^{-1} = I for normal " + std::to_string(n) +
                  " (error " + std::to_string(max_err) + ")");
   }
}

// ===== Test 4: Interior flux conservation (uniform Q => F = A_n Q) =====
void TestInteriorFluxConstantQ()
{
   std::cout << "Test 4: TestInteriorFluxConstantQ\n";
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   // For uniform Q: Q_L = Q_R => F = (A^+ + A^-) Q = A_n Q
   // This means no dissipation for uniform flow.
   srand(42);
   real_t Q_uniform[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
      Q_uniform[c] = (double)rand() / RAND_MAX * 1e6;

   real_t nor[3] = {0.0, 1.0, 0.0};
   real_t F_h[NUM_STATE];
   flux.Interior(nor, Q_uniform, Q_uniform, F_h);

   // Compare with A_n * Q (A_n = n_y * B)
   DenseMatrix Ay(NUM_STATE, NUM_STATE);
   flux.BuildJacobian(1, Ay);
   real_t F_exact[NUM_STATE];
   Ay.Mult(Q_uniform, F_exact);

   real_t max_err = 0.0;
   for (int c = 0; c < NUM_STATE; c++)
      max_err = std::max(max_err, std::abs(F_h[c] - F_exact[c]));

   TEST_ASSERT(max_err < 1e-6,
               "Uniform Q: F_godunov = A_n Q (error " + std::to_string(max_err) + ")");
}

// ===== Test 5: Interior flux conservation: F(Q_R, Q_L; -n) = -F(Q_L, Q_R; n) =====
void TestInteriorFluxConservation()
{
   std::cout << "Test 5: TestInteriorFluxConservation\n";
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   srand(123);
   real_t Q_L[NUM_STATE], Q_R[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      Q_L[c] = (double)rand() / RAND_MAX * 1e6;
      Q_R[c] = (double)rand() / RAND_MAX * 1e6;
   }

   real_t nor[3] = {1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
   real_t nor_neg[3] = {-nor[0], -nor[1], -nor[2]};

   real_t F1[NUM_STATE], F2[NUM_STATE];
   flux.Interior(nor, Q_L, Q_R, F1);
   flux.Interior(nor_neg, Q_R, Q_L, F2);

   // Should sum to zero (flux conservation)
   real_t max_err = 0.0;
   for (int c = 0; c < NUM_STATE; c++)
      max_err = std::max(max_err, std::abs(F1[c] + F2[c]));

   TEST_ASSERT(max_err < 1e-6,
               "F(QL,QR;n) + F(QR,QL;-n) = 0 (error " + std::to_string(max_err) + ")");
}

// ===== Test 6: Absorbing BC: outgoing P-wave => full flux, incoming = 0 =====
void TestAbsorbingOutgoingPWave()
{
   std::cout << "Test 6: TestAbsorbingOutgoingPWave\n";
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   real_t cp = flux.GetCp();
   real_t lp = lambda + 2.0*mu;

   // Outgoing P-wave in +x direction at x=+Lx boundary (normal = +x)
   real_t Q_pwave[NUM_STATE] = {-lp, -lambda, -lambda, 0, 0, 0, cp, 0, 0};

   real_t nor[3] = {1.0, 0.0, 0.0};
   real_t F_abs[NUM_STATE];
   flux.Absorbing(nor, Q_pwave, F_abs);

   // For a purely outgoing wave, F_abs = A_n Q (full flux, no reflected part)
   DenseMatrix Ax(NUM_STATE, NUM_STATE);
   flux.BuildJacobian(0, Ax);
   real_t F_exact[NUM_STATE];
   Ax.Mult(Q_pwave, F_exact);

   real_t max_err = 0.0;
   for (int c = 0; c < NUM_STATE; c++)
      max_err = std::max(max_err, std::abs(F_abs[c] - F_exact[c]));

   TEST_ASSERT(max_err < 1e-6,
               "Absorbing: outgoing P-wave flux = AQ (error " + std::to_string(max_err) + ")");
}

// ===== Test 6b: Absorbing BC: incoming P-wave produces zero flux =====
void TestAbsorbingIncomingPWave()
{
   std::cout << "Test 6b: TestAbsorbingIncomingPWave\n";
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   real_t cp = flux.GetCp();
   real_t lp = lambda + 2.0*mu;

   // Leftgoing P-wave eigenvector (eigenvalue -cp of A_x):
   // [+lp, +lambda, +lambda, 0, 0, 0, cp, 0, 0]
   // At x=+Lx boundary (normal = +x), this is an INCOMING wave.
   real_t e_left[NUM_STATE] = {lp, lambda, lambda, 0, 0, 0, cp, 0, 0};

   real_t nor[3] = {1.0, 0.0, 0.0};
   real_t F_abs[NUM_STATE];
   flux.Absorbing(nor, e_left, F_abs);

   // A^+ applied to the -cp eigenvector should give zero
   // (the entire wave is in the negative eigenspace of A_x).
   // Values are O(lp*cp) ≈ 1e14, so use relative tolerance.
   real_t F_norm = 0.0;
   for (int c = 0; c < NUM_STATE; c++)
      F_norm += F_abs[c] * F_abs[c];
   real_t Q_norm = 0.0;
   for (int c = 0; c < NUM_STATE; c++)
      Q_norm += e_left[c] * e_left[c];

   TEST_ASSERT(std::sqrt(F_norm) / std::sqrt(Q_norm) < 1e-10,
               "Absorbing: incoming P-wave gives zero flux (rel " +
               std::to_string(std::sqrt(F_norm) / std::sqrt(Q_norm)) + ")");
}

// ===== Test 7: Free surface: Godunov state has sigma . n = 0 =====
void TestFreeSurfaceZeroTraction()
{
   std::cout << "Test 7: TestFreeSurfaceZeroTraction\n";
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   srand(777);
   real_t Q_rand[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
      Q_rand[c] = (double)rand() / RAND_MAX * 1e6;

   // Free surface with normal in z-direction
   real_t nor[3] = {0.0, 0.0, 1.0};
   real_t t1[3], t2[3];
   GodunovFlux::BuildFrame(nor, t1, t2);

   // Rotate Q to face-local frame
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   GodunovFlux::BuildRotationInverse(nor, t1, t2, Tinv);
   real_t Q_rot[NUM_STATE];
   Tinv.Mult(Q_rand, Q_rot);

   // Build ghost state using gamma
   static const real_t gamma[NUM_STATE] = {-1, 1, 1, -1, 1, -1, 1, 1, 1};
   real_t Q_ghost_rot[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
      Q_ghost_rot[c] = gamma[c] * Q_rot[c];

   // Verify Godunov state has zero normal traction.
   // sigma_nn* = eta_p * (sigma_nn/Zp + sigma_nn_ghost/Zp + v_n_ghost - v_n)
   real_t Zp = flux.GetZp(), Zs = flux.GetZs();
   real_t eta_p = Zp / 2.0;  // homogeneous
   real_t eta_s = Zs / 2.0;

   real_t sigma_nn_star = eta_p * (Q_rot[SXX]/Zp + Q_ghost_rot[SXX]/Zp
                                   + Q_ghost_rot[VX] - Q_rot[VX]);
   TEST_NEAR(sigma_nn_star, 0.0, 1e-6, "Free surface: sigma_nn* = 0");

   real_t tau_1_star = eta_s * (Q_rot[SXY]/Zs + Q_ghost_rot[SXY]/Zs
                                + Q_ghost_rot[VY] - Q_rot[VY]);
   TEST_NEAR(tau_1_star, 0.0, 1e-6, "Free surface: tau_nt1* = 0");

   real_t tau_2_star = eta_s * (Q_rot[SXZ]/Zs + Q_ghost_rot[SXZ]/Zs
                                + Q_ghost_rot[VZ] - Q_rot[VZ]);
   TEST_NEAR(tau_2_star, 0.0, 1e-6, "Free surface: tau_nt2* = 0");

   // End-to-end: call FreeSurface() and verify it matches Interior() with ghost
   DenseMatrix T(NUM_STATE, NUM_STATE);
   GodunovFlux::BuildRotation(nor, t1, t2, T);
   real_t Q_ghost_global[NUM_STATE];
   T.Mult(Q_ghost_rot, Q_ghost_global);

   real_t F_free[NUM_STATE], F_interior_ghost[NUM_STATE];
   flux.FreeSurface(nor, Q_rand, F_free);
   flux.Interior(nor, Q_rand, Q_ghost_global, F_interior_ghost);

   real_t max_diff = 0.0;
   for (int c = 0; c < NUM_STATE; c++)
      max_diff = std::max(max_diff, std::abs(F_free[c] - F_interior_ghost[c]));
   TEST_ASSERT(max_diff < 1e-6,
               "FreeSurface() matches Interior(Q, ghost) end-to-end (diff " +
               std::to_string(max_diff) + ")");
}

// ===== Test 8: Split flux A^+ + A^- = A =====
void TestSplitFluxReconstructsA()
{
   std::cout << "Test 8: TestSplitFluxReconstructsA\n";
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   const DenseMatrix &Ax = flux.GetAx();
   const DenseMatrix &Ap = flux.GetAxPlus();
   const DenseMatrix &Am = flux.GetAxMinus();

   real_t max_err = 0.0;
   for (int i = 0; i < NUM_STATE; i++)
      for (int j = 0; j < NUM_STATE; j++)
         max_err = std::max(max_err, std::abs(Ap(i,j) + Am(i,j) - Ax(i,j)));

   TEST_ASSERT(max_err < 1e-6,
               "A^+ + A^- = A (max error " + std::to_string(max_err) + ")");
}

// ADER I-05 Phase 2: GetReferenceStarMatrix(dir) must equal BuildJacobian(dir).
void TestReferenceStarMatricesMatchJacobians()
{
   std::cout << "Test 9: TestReferenceStarMatricesMatchJacobians\n";
   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   GodunovFlux flux(lambda, mu, rho);

   for (int d = 0; d < 3; d++)
   {
      DenseMatrix A_ref(NUM_STATE, NUM_STATE);
      flux.BuildJacobian(d, A_ref);
      const DenseMatrix &A_star = flux.GetReferenceStarMatrix(d);
      real_t max_err = 0.0;
      for (int i = 0; i < NUM_STATE; i++)
      {
         for (int j = 0; j < NUM_STATE; j++)
         {
            max_err = std::max(max_err, std::abs(A_ref(i,j) - A_star(i,j)));
         }
      }
      TEST_ASSERT(max_err == 0.0,
                  "GetReferenceStarMatrix(" + std::to_string(d) +
                  ") == BuildJacobian(" + std::to_string(d) + ") (max error "
                  + std::to_string(max_err) + ")");
   }
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "GodunovFlux Unit Tests (Phase 1)\n";
   std::cout << "========================================\n\n";

   TestJacobianSymmetryABC();
   TestEigenvalues();
   TestRotationOrthogonality();
   TestInteriorFluxConstantQ();
   TestInteriorFluxConservation();
   TestAbsorbingOutgoingPWave();
   TestAbsorbingIncomingPWave();
   TestFreeSurfaceZeroTraction();
   TestSplitFluxReconstructsA();
   TestReferenceStarMatricesMatchJacobians();

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";

   return (num_failed > 0) ? 1 : 0;
}
