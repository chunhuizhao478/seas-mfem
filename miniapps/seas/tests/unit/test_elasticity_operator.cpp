// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Test suite for ElasticityDomainOperator

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../config/bp5_params.hpp"
#include "../../fault/fault_basis.hpp"
#include "../../fault/fault_geometry.hpp"
#include "test_macros.hpp"

#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

// Helper: Create a simple 3D hex mesh with fault at x1=0
// Domain: [-Lx, Lx] x [-Ly, Ly] x [0, Lz]
// Boundary attributes:
//   1 = x1=-Lx (normal-)
//   2 = x1=+Lx (normal+)
//   3 = x2=+Ly (strike+)
//   4 = x2=-Ly (strike-)
//   5 = x3=0 (free surface)
//   6 = x3=Lz (bottom)
Mesh CreateTestMesh3D(int nx, int ny, int nz,
                       real_t Lx, real_t Ly, real_t Lz)
{
   // Create mesh centered at origin: [-Lx, Lx] x [-Ly, Ly] x [0, Lz]
   Mesh mesh = Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz,
                                      Element::HEXAHEDRON,
                                      2.0 * Lx, 2.0 * Ly, Lz);

   // Shift so x1 ranges [-Lx, Lx] and x2 ranges [-Ly, Ly]
   Vector shift(3);
   shift(0) = -Lx;
   shift(1) = -Ly;
   shift(2) = 0.0;

   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *v = mesh.GetVertex(i);
      for (int d = 0; d < 3; d++)
      {
         v[d] += shift(d);
      }
   }

   // Set boundary attributes based on face location
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      ElementTransformation *T = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(T->GetGeometryType());
      T->SetIntPoint(&ip);
      Vector center(3);
      T->Transform(ip, center);

      real_t tol = 1e-6;
      if (std::abs(center(0) - (-Lx)) < tol)      { mesh.SetBdrAttribute(be, 1); }
      else if (std::abs(center(0) - Lx) < tol)     { mesh.SetBdrAttribute(be, 2); }
      else if (std::abs(center(1) - Ly) < tol)     { mesh.SetBdrAttribute(be, 3); }
      else if (std::abs(center(1) - (-Ly)) < tol)  { mesh.SetBdrAttribute(be, 4); }
      else if (std::abs(center(2) - 0.0) < tol)    { mesh.SetBdrAttribute(be, 5); }
      else if (std::abs(center(2) - Lz) < tol)     { mesh.SetBdrAttribute(be, 6); }
   }

   // Update boundary attribute list
   mesh.SetAttributes();

   return mesh;
}

// =============================================================================
// Test 1: Constructor and basic properties
// =============================================================================
void TestConstruction()
{
   std::cout << "\n--- Test: ElasticityDomainOperator Construction ---\n";

   real_t Lx = 10.0, Ly = 10.0, Lz = 10.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   BP5Params params;
   real_t lambda = params.lambda();
   real_t mu = params.mu();
   real_t Vp = params.Vp;
   real_t Wf = Lz;  // All within fault depth
   real_t lf = 2.0 * Ly;

   ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, Vp, Wf, lf,
                                      DGMethod::BR2);

   TEST_ASSERT(op.NumComponents() == 3, "NumComponents() == 3");
   TEST_ASSERT(op.Dimension() == 3, "Dimension() == 3");
   TEST_ASSERT(op.NumSlipComponents() == 2, "NumSlipComponents() == 2");
   TEST_NEAR(op.GetShearModulus(), mu, 1e-6, "Shear modulus matches");
   TEST_ASSERT(op.GetFaultBasis() != nullptr, "FaultBasis is available");
}

// =============================================================================
// Test 2: Fault detection at x1=0
// =============================================================================
void TestFaultDetection()
{
   std::cout << "\n--- Test: Fault Detection at x1=0 ---\n";

   // 4x2x2 mesh: [-2, 2] x [-1, 1] x [0, 1]
   // With fault at x1=0, we expect 2*2 = 4 interior fault faces
   real_t Lx = 2.0, Ly = 1.0, Lz = 1.0;
   Mesh mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);

   ElasticityDomainOperator<Mesh> op(mesh, 1, 1.0, 1.0, 1e-9, Lz, 2.0 * Ly,
                                      DGMethod::BR2);

   int nf = op.GetNumFaultDOFs();
   TEST_ASSERT(nf > 0, "Found fault faces at x1=0");

   // Check fault depths
   Vector depths;
   op.GetFaultDepths(depths);
   TEST_ASSERT(depths.Size() == nf, "Depths vector has correct size");

   // Check 2D coordinates
   Vector x2, x3;
   op.GetFaultCoords2D(x2, x3);
   TEST_ASSERT(x2.Size() == nf, "x2 coordinates have correct size");
   TEST_ASSERT(x3.Size() == nf, "x3 coordinates have correct size");

   // All fault faces should be at x1=0
   std::cout << "  Found " << nf << " fault faces\n";
}

// =============================================================================
// Test 3: FaultBasis properties
// =============================================================================
void TestFaultBasis()
{
   std::cout << "\n--- Test: FaultBasis for ElasticityDomainOperator ---\n";

   real_t Lx = 2.0, Ly = 1.0, Lz = 1.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   ElasticityDomainOperator<Mesh> op(mesh, 1, 1.0, 1.0, 1e-9, Lz, 2.0 * Ly,
                                      DGMethod::BR2);

   const FaultBasis *fb = op.GetFaultBasis();
   TEST_ASSERT(fb != nullptr, "FaultBasis is non-null");

   if (fb && fb->NumFaces() > 0)
   {
      TEST_ASSERT(fb->Dimension() == 3, "FaultBasis dimension is 3");
      TEST_ASSERT(fb->NumTangentComponents() == 2, "2 tangential components");

      // Check first face basis
      const auto &b = fb->GetBasis(0);

      // Normal should be approximately (1, 0, 0) (fault at x1=0)
      real_t n_len = std::sqrt(b.normal[0] * b.normal[0] +
                                b.normal[1] * b.normal[1] +
                                b.normal[2] * b.normal[2]);
      TEST_NEAR(n_len, 1.0, 1e-12, "Normal is unit vector");
      TEST_NEAR(std::abs(b.normal[0]), 1.0, 1e-6,
                "Normal is approximately along x1");

      // tangent1 (dip) and tangent2 (strike) should be orthonormal
      real_t t1_len = std::sqrt(b.tangent1[0] * b.tangent1[0] +
                                 b.tangent1[1] * b.tangent1[1] +
                                 b.tangent1[2] * b.tangent1[2]);
      real_t t2_len = std::sqrt(b.tangent2[0] * b.tangent2[0] +
                                 b.tangent2[1] * b.tangent2[1] +
                                 b.tangent2[2] * b.tangent2[2]);
      TEST_NEAR(t1_len, 1.0, 1e-12, "Tangent1 is unit vector");
      TEST_NEAR(t2_len, 1.0, 1e-12, "Tangent2 is unit vector");

      // Orthogonality: n·t1 = 0, n·t2 = 0, t1·t2 = 0
      real_t n_dot_t1 = b.normal[0]*b.tangent1[0] + b.normal[1]*b.tangent1[1] +
                         b.normal[2]*b.tangent1[2];
      real_t n_dot_t2 = b.normal[0]*b.tangent2[0] + b.normal[1]*b.tangent2[1] +
                         b.normal[2]*b.tangent2[2];
      real_t t1_dot_t2 = b.tangent1[0]*b.tangent2[0] +
                          b.tangent1[1]*b.tangent2[1] +
                          b.tangent1[2]*b.tangent2[2];
      TEST_NEAR(n_dot_t1, 0.0, 1e-12, "Normal perpendicular to tangent1");
      TEST_NEAR(n_dot_t2, 0.0, 1e-12, "Normal perpendicular to tangent2");
      TEST_NEAR(t1_dot_t2, 0.0, 1e-12, "Tangent1 perpendicular to tangent2");
   }
}

// =============================================================================
// Test 4: Zero slip → zero displacement (no loading)
// =============================================================================
void TestZeroSlipEquilibrium()
{
   std::cout << "\n--- Test: Zero Slip + Zero Loading → Zero Displacement ---\n";

   real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   // Use IP method for simpler testing
   ElasticityDomainOperator<Mesh> op(mesh, 1, 1.0, 1.0, 0.0, Lz, 2.0 * Ly,
                                      DGMethod::IP);

   int nf = op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   Vector slip_bc(2 * nf);
   slip_bc = 0.0;

   GridFunction u(&op.GetFESpace());
   u = 0.0;

   op.Solve(0.0, slip_bc, u);

   real_t u_norm = u.Norml2();
   TEST_ASSERT(u_norm < 1e-10, "Zero slip + zero loading gives zero displacement");
}

// =============================================================================
// Test 5: Traction extraction
// =============================================================================
void TestTractionExtraction()
{
   std::cout << "\n--- Test: Traction Extraction ---\n";

   real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   ElasticityDomainOperator<Mesh> op(mesh, 1, 1.0, 1.0, 0.0, Lz, 2.0 * Ly,
                                      DGMethod::IP);

   int nf = op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   Vector slip_bc(2 * nf);
   slip_bc = 0.0;

   GridFunction u(&op.GetFESpace());
   u = 0.0;

   op.Solve(0.0, slip_bc, u);

   Vector traction;
   op.ComputeTraction(u, slip_bc, traction);

   TEST_ASSERT(traction.Size() == 2 * nf,
               "Traction vector has correct size (2 * num_fault_dofs)");

   // With zero slip and zero loading, traction should be zero
   real_t trac_norm = traction.Norml2();
   TEST_ASSERT(trac_norm < 1e-8,
               "Zero slip gives near-zero traction");
}

// =============================================================================
// Test 6: Stiffness matrix assembly (both methods)
// =============================================================================
void TestStiffnessAssembly()
{
   std::cout << "\n--- Test: Stiffness Assembly (IP and BR2) ---\n";

   real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   // IP method
   {
      ElasticityDomainOperator<Mesh> op_ip(mesh, 1, 1.0, 1.0, 0.0, Lz, 2.0 * Ly,
                                            DGMethod::IP);
      int nf = op_ip.GetNumFaultDOFs();
      if (nf > 0)
      {
         Vector slip(2 * nf);
         slip = 0.0;
         GridFunction u(&op_ip.GetFESpace());
         u = 0.0;
         op_ip.Solve(0.0, slip, u);  // Triggers assembly
         TEST_ASSERT(true, "IP stiffness assembly succeeds");
      }
   }

   // BR2 method
   {
      ElasticityDomainOperator<Mesh> op_br2(mesh, 1, 1.0, 1.0, 0.0, Lz, 2.0 * Ly,
                                              DGMethod::BR2);
      int nf = op_br2.GetNumFaultDOFs();
      if (nf > 0)
      {
         Vector slip(2 * nf);
         slip = 0.0;
         GridFunction u(&op_br2.GetFESpace());
         u = 0.0;
         op_br2.Solve(0.0, slip, u);  // Triggers assembly
         TEST_ASSERT(true, "BR2 stiffness assembly succeeds");
      }
   }
}

// =============================================================================
// Test 7: FaultGeometry 3D construction with BP5Params
// =============================================================================
void TestFaultGeometry3D()
{
   std::cout << "\n--- Test: FaultGeometry 3D with BP5Params ---\n";

   real_t Lx = 50e3, Ly = 60e3, Lz = 40e3;
   // Coarse mesh
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   BP5Params params;

   ElasticityDomainOperator<Mesh> op(mesh, 1, params.lambda(), params.mu(),
                                      params.Vp, params.Wf, params.lf,
                                      DGMethod::IP);

   int nf = op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found on this mesh)\n";
      return;
   }

   FaultGeometry<Mesh> geom(op, params);

   TEST_ASSERT(geom.NumFaultDOFs() == nf, "FaultGeometry DOF count matches");
   TEST_ASSERT(geom.IsBP5(), "FaultGeometry identifies as BP5");

   const Vector &a_vals = geom.GetAValues();
   const Vector &eta_vals = geom.GetEtaValues();
   const Vector &dc_vals = geom.GetDcValues();
   const Vector &tau_pre = geom.GetTauPre();

   TEST_ASSERT(a_vals.Size() == nf, "a values have correct size");
   TEST_ASSERT(eta_vals.Size() == nf, "eta values have correct size");
   TEST_ASSERT(dc_vals.Size() == nf, "Dc values have correct size");
   TEST_ASSERT(tau_pre.Size() == 2 * nf, "tau_pre has correct size (2*N)");

   // eta should be constant
   real_t eta_expected = params.eta();
   for (int i = 0; i < nf; i++)
   {
      TEST_NEAR(eta_vals(i), eta_expected, 1.0,
                "eta matches BP5 eta()");
   }
}

// =============================================================================
// Test 7b: FaultGeometry 3D value verification (§6 test gaps)
// =============================================================================
void TestFaultGeometry3DValues()
{
   std::cout << "\n--- Test: FaultGeometry 3D Value Verification ---\n";

   real_t Lx = 50e3, Ly = 60e3, Lz = 40e3;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   BP5Params params;

   ElasticityDomainOperator<Mesh> op(mesh, 1, params.lambda(), params.mu(),
                                      params.Vp, params.Wf, params.lf,
                                      DGMethod::IP);

   int nf = op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found on this mesh)\n";
      return;
   }

   FaultGeometry<Mesh> geom(op, params);

   // GetVInit size and non-zero
   const Vector &v_init = geom.GetVInit();
   TEST_ASSERT(v_init.Size() == 2 * nf, "V_init has correct size (2*N)");
   TEST_ASSERT(v_init.Norml2() > 0.0, "V_init is non-zero");

   // Coordinates within bounds
   const Vector &x2 = geom.GetCoordsX2();
   const Vector &x3 = geom.GetCoordsX3();
   TEST_ASSERT(x2.Size() == nf, "x2 coords correct size");
   TEST_ASSERT(x3.Size() == nf, "x3 coords correct size");
   for (int i = 0; i < nf; i++)
   {
      TEST_ASSERT(std::abs(x2(i)) <= Ly + 1.0, "x2 within mesh bounds");
      TEST_ASSERT(x3(i) >= -1.0 && x3(i) <= Lz + 1.0, "x3 within mesh bounds");
   }

   // a values within [a0, amax]
   const Vector &a = geom.GetAValues();
   for (int i = 0; i < nf; i++)
   {
      TEST_ASSERT(a(i) >= params.a0 - 1e-15 && a(i) <= params.amax + 1e-15,
                  "a in [a0, amax]");
   }

   // dc values are L0 or L_nuc
   const Vector &dc = geom.GetDcValues();
   for (int i = 0; i < nf; i++)
   {
      TEST_ASSERT(std::abs(dc(i) - params.L0) < 1e-15 ||
                  std::abs(dc(i) - params.L_nuc) < 1e-15,
                  "dc is L0 or L_nuc");
   }

   // tau_pre physically reasonable (non-zero, < 100 MPa per-DOF average)
   const Vector &tau = geom.GetTauPre();
   real_t tau_norm = tau.Norml2();
   TEST_ASSERT(tau_norm > 0.0, "tau_pre is non-zero");
   TEST_ASSERT(tau_norm / std::sqrt(2.0 * nf) < 1e8,
               "tau_pre per-DOF < 100 MPa average");

   // GetBP5Params accessor
   const BP5Params &p = geom.GetBP5Params();
   TEST_NEAR(p.b, params.b, 1e-15, "GetBP5Params returns correct params");

   // Print() for BP5 path (exercises §2.1 fix)
   std::ostringstream oss;
   geom.Print(oss);
   std::string output = oss.str();
   TEST_ASSERT(output.find("VW DOFs") != std::string::npos,
               "Print() outputs VW DOFs for BP5");
   std::cout << "  Print output:\n" << output;
}

// =============================================================================
// Test 8: Non-zero slip produces non-zero displacement (IP method)
// =============================================================================
void TestNonZeroSlipIP()
{
   std::cout << "\n--- Test: Non-Zero Slip Displacement (IP) ---\n";

   real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   real_t lambda = 1.0, mu = 1.0;
   ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, 0.0, Lz, 2.0 * Ly,
                                      DGMethod::IP);

   int nf = op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   // Apply uniform dip-direction slip of 1.0 on all fault faces
   Vector slip_bc(2 * nf);
   slip_bc = 0.0;
   for (int i = 0; i < nf; i++)
   {
      slip_bc(2 * i) = 1.0;      // dip slip = 1.0
      slip_bc(2 * i + 1) = 0.0;  // strike slip = 0.0
   }

   GridFunction u(&op.GetFESpace());
   u = 0.0;
   op.Solve(0.0, slip_bc, u);

   real_t u_norm = u.Norml2();
   TEST_ASSERT(u_norm > 1e-6,
               "IP: non-zero slip produces non-zero displacement");

   // Extract traction
   Vector traction;
   op.ComputeTraction(u, slip_bc, traction);
   TEST_ASSERT(traction.Size() == 2 * nf,
               "IP: traction vector has correct size");

   // With slip applied, traction should be non-zero
   real_t trac_norm = traction.Norml2();
   TEST_ASSERT(trac_norm > 1e-6,
               "IP: non-zero slip produces non-zero traction");
   std::cout << "  IP: |u| = " << u_norm << ", |traction| = " << trac_norm << "\n";
}

// =============================================================================
// Test 9: Non-zero slip produces non-zero displacement (BR2 method)
// =============================================================================
void TestNonZeroSlipBR2()
{
   std::cout << "\n--- Test: Non-Zero Slip Displacement (BR2) ---\n";

   real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   real_t lambda = 1.0, mu = 1.0;
   ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, 0.0, Lz, 2.0 * Ly,
                                      DGMethod::BR2);

   int nf = op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   // Apply uniform dip-direction slip of 1.0 on all fault faces
   Vector slip_bc(2 * nf);
   slip_bc = 0.0;
   for (int i = 0; i < nf; i++)
   {
      slip_bc(2 * i) = 1.0;      // dip slip = 1.0
      slip_bc(2 * i + 1) = 0.0;  // strike slip = 0.0
   }

   GridFunction u(&op.GetFESpace());
   u = 0.0;
   op.Solve(0.0, slip_bc, u);

   real_t u_norm = u.Norml2();
   TEST_ASSERT(u_norm > 1e-6,
               "BR2: non-zero slip produces non-zero displacement");

   // Extract traction
   Vector traction;
   op.ComputeTraction(u, slip_bc, traction);
   TEST_ASSERT(traction.Size() == 2 * nf,
               "BR2: traction vector has correct size");

   real_t trac_norm = traction.Norml2();
   TEST_ASSERT(trac_norm > 1e-6,
               "BR2: non-zero slip produces non-zero traction");
   std::cout << "  BR2: |u| = " << u_norm << ", |traction| = " << trac_norm << "\n";
}

// =============================================================================
// Test 10: Dirichlet loading produces non-zero displacement
// =============================================================================
void TestDirichletLoading()
{
   std::cout << "\n--- Test: Dirichlet Loading (time > 0) ---\n";

   real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   real_t Vp = 1.0;  // Large Vp so loading is significant
   real_t lambda = 1.0, mu = 1.0;

   // Test with IP
   {
      ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, Vp, Lz, 2.0 * Ly,
                                         DGMethod::IP);

      int nf = op.GetNumFaultDOFs();
      if (nf == 0) { std::cout << "  (Skipped IP: no fault faces)\n"; }
      else
      {
         Vector slip_bc(2 * nf);
         slip_bc = 0.0;  // Zero slip, loading from boundary only

         GridFunction u(&op.GetFESpace());
         u = 0.0;

         real_t time = 1.0;  // Non-zero time
         op.Solve(time, slip_bc, u);

         real_t u_norm = u.Norml2();
         TEST_ASSERT(u_norm > 1e-6,
                     "IP: Dirichlet loading at t>0 produces non-zero displacement");
         std::cout << "  IP with Vp*t/2=" << Vp * time / 2.0
                   << ": |u| = " << u_norm << "\n";
      }
   }

   // Test with BR2
   {
      ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, Vp, Lz, 2.0 * Ly,
                                         DGMethod::BR2);

      int nf = op.GetNumFaultDOFs();
      if (nf == 0) { std::cout << "  (Skipped BR2: no fault faces)\n"; }
      else
      {
         Vector slip_bc(2 * nf);
         slip_bc = 0.0;

         GridFunction u(&op.GetFESpace());
         u = 0.0;

         real_t time = 1.0;
         op.Solve(time, slip_bc, u);

         real_t u_norm = u.Norml2();
         TEST_ASSERT(u_norm > 1e-6,
                     "BR2: Dirichlet loading at t>0 produces non-zero displacement");
         std::cout << "  BR2 with Vp*t/2=" << Vp * time / 2.0
                   << ": |u| = " << u_norm << "\n";
      }
   }
}

// =============================================================================
// Test 10b: Dirichlet loading at fault-normal walls produces shear traction
// =============================================================================
void TestDirichletLoadingShearTraction()
{
   std::cout << "\n--- Test: Dirichlet Loading Produces Shear Traction ---\n";

   // Use a mesh large enough to have fault faces
   real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   real_t Vp = 1.0;
   real_t lambda = 1.0, mu = 1.0;

   for (int method = 0; method < 2; method++)
   {
      DGMethod dg = (method == 0) ? DGMethod::IP : DGMethod::BR2;
      std::string label = (method == 0) ? "IP" : "BR2";

      ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, Vp, Lz, 2.0 * Ly,
                                         dg);

      int nf = op.GetNumFaultDOFs();
      if (nf == 0)
      {
         std::cout << "  (Skipped " << label << ": no fault faces)\n";
         continue;
      }

      // Zero slip, loading from boundary only at t=1
      Vector slip_bc(2 * nf);
      slip_bc = 0.0;

      GridFunction u(&op.GetFESpace());
      u = 0.0;

      real_t time = 1.0;
      op.Solve(time, slip_bc, u);

      // Extract traction on the fault
      Vector traction;
      op.ComputeTraction(u, slip_bc, traction);

      // Traction layout (Tandem convention): [dip_0, strike_0, dip_1, strike_1, ...]
      //   traction(2*i)   = tau_local[0] = dip component (tangent1)
      //   traction(2*i+1) = tau_local[1] = strike component (tangent2)
      //
      // Dirichlet loading u_y = ±Vp*t/2 at ±x walls creates σ_xy.
      // With up=(0,0,1), n=(1,0,0):
      //   strike = up × n = (0,1,0) → tangent2
      //   dip = -(strike × n) = (0,0,1) → tangent1 (negated for Tandem/SCEC convention)
      // σ_xy projects onto strike (tangent2), so traction(2*i+1) should dominate.

      // Check traction is non-zero
      real_t trac_norm = traction.Norml2();
      TEST_ASSERT(trac_norm > 1e-6,
                  (label + ": Dirichlet loading produces non-zero traction").c_str());

      // Separate dip (index 2*i) and strike (index 2*i+1) components
      real_t strike_sum = 0.0, dip_sum = 0.0;
      for (int i = 0; i < nf; i++)
      {
         dip_sum += traction(2 * i) * traction(2 * i);
         strike_sum += traction(2 * i + 1) * traction(2 * i + 1);
      }
      strike_sum = std::sqrt(strike_sum);
      dip_sum = std::sqrt(dip_sum);

      std::cout << "  " << label << ": |trac_strike| = " << strike_sum
                << ", |trac_dip| = " << dip_sum << "\n";

      TEST_ASSERT(strike_sum > 1e-6,
                  (label + ": Strike traction is non-zero from shear loading").c_str());

      // Strike component should dominate over dip for pure shear loading
      if (dip_sum > 1e-10)
      {
         TEST_ASSERT(strike_sum > dip_sum,
                     (label + ": Strike traction dominates over dip").c_str());
      }
      else
      {
         TEST_ASSERT(true,
                     (label + ": Dip traction is negligible (pure shear)").c_str());
      }

      // Sign check: with up=(0,0,1), strike=(0,1,0) [unchanged by dip negation],
      // the strike traction component should be positive.
      real_t avg_strike = 0.0;
      for (int i = 0; i < nf; i++)
      {
         avg_strike += traction(2 * i + 1);
      }
      avg_strike /= nf;
      std::cout << "  " << label << ": avg strike traction = " << avg_strike << "\n";

      TEST_ASSERT(avg_strike > 0.0,
                  (label + ": Strike traction sign is positive (right-lateral)").c_str());
   }
}

// =============================================================================
// Test 11: BR2 vs IP comparison — both methods should produce similar results
// =============================================================================
void TestBR2vsIP()
{
   std::cout << "\n--- Test: BR2 vs IP Comparison ---\n";

   real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   real_t lambda = 1.0, mu = 1.0;

   // IP operator
   ElasticityDomainOperator<Mesh> op_ip(mesh, 1, lambda, mu, 0.0, Lz, 2.0 * Ly,
                                         DGMethod::IP);

   // BR2 operator
   ElasticityDomainOperator<Mesh> op_br2(mesh, 1, lambda, mu, 0.0, Lz, 2.0 * Ly,
                                          DGMethod::BR2);

   int nf_ip = op_ip.GetNumFaultDOFs();
   int nf_br2 = op_br2.GetNumFaultDOFs();
   TEST_ASSERT(nf_ip == nf_br2,
               "IP and BR2 detect same number of fault faces");

   if (nf_ip == 0 || nf_br2 == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   // Apply the same uniform slip
   Vector slip_ip(2 * nf_ip), slip_br2(2 * nf_br2);
   slip_ip = 0.0;
   slip_br2 = 0.0;
   for (int i = 0; i < nf_ip; i++)
   {
      slip_ip(2 * i) = 1.0;
      slip_br2(2 * i) = 1.0;
   }

   GridFunction u_ip(&op_ip.GetFESpace());
   GridFunction u_br2(&op_br2.GetFESpace());
   u_ip = 0.0;
   u_br2 = 0.0;

   op_ip.Solve(0.0, slip_ip, u_ip);
   op_br2.Solve(0.0, slip_br2, u_br2);

   real_t norm_ip = u_ip.Norml2();
   real_t norm_br2 = u_br2.Norml2();

   // Both should be non-zero
   TEST_ASSERT(norm_ip > 1e-6, "IP displacement is non-zero");
   TEST_ASSERT(norm_br2 > 1e-6, "BR2 displacement is non-zero");

   // The two solutions should be in the same ballpark (same order of magnitude).
   // On a coarse mesh they won't match exactly since the DG methods differ,
   // but the ratio should be reasonable (within a factor of 10).
   real_t ratio = norm_br2 / norm_ip;
   std::cout << "  |u_ip| = " << norm_ip
             << ", |u_br2| = " << norm_br2
             << ", ratio = " << ratio << "\n";
   TEST_ASSERT(ratio > 0.1 && ratio < 10.0,
               "BR2 and IP displacement norms are within factor of 10");

   // Compare tractions
   Vector trac_ip, trac_br2;
   op_ip.ComputeTraction(u_ip, slip_ip, trac_ip);
   op_br2.ComputeTraction(u_br2, slip_br2, trac_br2);

   real_t trac_norm_ip = trac_ip.Norml2();
   real_t trac_norm_br2 = trac_br2.Norml2();
   std::cout << "  |trac_ip| = " << trac_norm_ip
             << ", |trac_br2| = " << trac_norm_br2 << "\n";
}

// =============================================================================
// Test 12: BR2 default method verification
// =============================================================================
void TestBR2Default()
{
   std::cout << "\n--- Test: BR2 is Default Method ---\n";

   real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   // Construct without specifying method — should default to BR2
   ElasticityDomainOperator<Mesh> op(mesh, 1, 1.0, 1.0, 0.0, Lz, 2.0 * Ly);

   // Verify it works (solves correctly with default)
   int nf = op.GetNumFaultDOFs();
   if (nf == 0) { std::cout << "  (Skipped: no fault faces)\n"; return; }

   Vector slip(2 * nf);
   slip = 0.0;
   for (int i = 0; i < nf; i++) { slip(2 * i) = 0.5; }

   GridFunction u(&op.GetFESpace());
   u = 0.0;
   op.Solve(0.0, slip, u);

   real_t u_norm = u.Norml2();
   TEST_ASSERT(u_norm > 1e-6,
               "Default (BR2) method produces non-zero displacement");
   std::cout << "  Default method: |u| = " << u_norm << "\n";
}

// =============================================================================
// Test 15: Traction with DG penalty correction
// =============================================================================
void TestTractionWithPenaltyCorrection()
{
   std::cout << "\n--- Test: Traction With DG Penalty Correction ---\n";

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);
   real_t lambda = 1.0, mu = 1.0;

   real_t ip_avg_strike = 0.0, br2_avg_strike = 0.0;

   for (int method = 0; method < 2; method++)
   {
      DGMethod dg = (method == 0) ? DGMethod::IP : DGMethod::BR2;
      std::string label = (method == 0) ? "IP" : "BR2";

      ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, 0.0, Lz, 2*Ly, dg);
      int nf = op.GetNumFaultDOFs();
      if (nf == 0)
      {
         std::cout << "  (Skipped " << label << ": no fault faces)\n";
         continue;
      }

      // --- Sub-test A: Zero slip → zero traction ---
      {
         Vector slip0(2 * nf);
         slip0 = 0.0;
         GridFunction u0(&op.GetFESpace());
         u0 = 0.0;
         op.Solve(0.0, slip0, u0);
         Vector trac0;
         op.ComputeTraction(u0, slip0, trac0);
         real_t trac0_norm = trac0.Norml2();
         TEST_ASSERT(trac0_norm < 1e-8,
                     (label + ": Zero slip gives zero traction").c_str());
         std::cout << "    " << label << " zero-slip traction norm: "
                   << trac0_norm << "\n";
      }

      // --- Sub-test B: Uniform strike slip → stress drop ---
      {
         Vector slip1(2 * nf);
         slip1 = 0.0;
         for (int i = 0; i < nf; i++) { slip1(2*i+1) = 1.0; }  // strike slip

         GridFunction u1(&op.GetFESpace());
         u1 = 0.0;
         op.Solve(0.0, slip1, u1);
         Vector trac1;
         op.ComputeTraction(u1, slip1, trac1);

         real_t avg_strike = 0.0;
         for (int i = 0; i < nf; i++) { avg_strike += trac1(2*i+1); }
         avg_strike /= nf;

         // Positive slip should produce negative traction (stress drop)
         TEST_ASSERT(avg_strike < 0,
                     (label + ": Strike slip causes negative traction (stress drop)").c_str());
         std::cout << "    " << label << " avg strike traction: "
                   << avg_strike << "\n";

         // Traction should be bounded
         real_t trac1_max = trac1.Normlinf();
         TEST_ASSERT(trac1_max < 100.0,
                     (label + ": Traction is bounded").c_str());
         std::cout << "    " << label << " max traction: "
                   << trac1_max << "\n";

         if (method == 0) { ip_avg_strike = avg_strike; }
         else { br2_avg_strike = avg_strike; }
      }
   }

   // --- Sub-test C: IP vs BR2 consistency ---
   if (std::abs(ip_avg_strike) > 1e-12 && std::abs(br2_avg_strike) > 1e-12)
   {
      // Both should have same sign (negative)
      TEST_ASSERT(ip_avg_strike * br2_avg_strike > 0,
                  "IP and BR2 strike traction have same sign");
      // On coarse meshes with order 1, IP has much larger penalty than BR2,
      // so allow a wide ratio. The key check is same sign (stress drop).
      real_t ratio = std::abs(ip_avg_strike / br2_avg_strike);
      TEST_ASSERT(ratio > 0.001 && ratio < 1000.0,
                  "IP and BR2 strike traction within 3 orders of magnitude");
      std::cout << "    IP/BR2 ratio: " << ratio << "\n";
   }
}

// =============================================================================
// Test 16: BR2 traction correction consistency — correct scaling under mesh refinement
//
// The BR2 traction correction in ComputeTraction uses unnormalized normals
// (from CalcOrtho) with the centroid quadrature weight for face_int (matching
// the bilinear form integrator), and unit normals for TestNormal (point eval).
// The correction scales as μ*g/h, so it should not blow up or vanish with
// mesh size changes.
// =============================================================================
void TestBR2TractionCorrectionConsistency()
{
   std::cout << "\n--- Test: BR2 Traction Correction Consistency (Normal Invariance) ---\n";

   real_t lambda = 1.0, mu = 1.0;

   // Mesh 1: small domain
   real_t Lx1 = 2.0, Ly1 = 2.0, Lz1 = 2.0;
   Mesh mesh1 = CreateTestMesh3D(1, 1, 1, Lx1, Ly1, Lz1);
   ElasticityDomainOperator<Mesh> op1(mesh1, 1, lambda, mu, 0.0, Lz1, 2.0 * Ly1,
                                       DGMethod::BR2);

   // Mesh 2: 2x larger domain (face areas 4x larger)
   real_t Lx2 = 4.0, Ly2 = 4.0, Lz2 = 4.0;
   Mesh mesh2 = CreateTestMesh3D(1, 1, 1, Lx2, Ly2, Lz2);
   ElasticityDomainOperator<Mesh> op2(mesh2, 1, lambda, mu, 0.0, Lz2, 2.0 * Ly2,
                                       DGMethod::BR2);

   int nf1 = op1.GetNumFaultDOFs();
   int nf2 = op2.GetNumFaultDOFs();
   if (nf1 == 0 || nf2 == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   // Apply same unit dip slip on both
   Vector slip1(2 * nf1), slip2(2 * nf2);
   slip1 = 0.0; slip2 = 0.0;
   for (int i = 0; i < nf1; i++) { slip1(2 * i) = 1.0; }
   for (int i = 0; i < nf2; i++) { slip2(2 * i) = 1.0; }

   GridFunction u1(&op1.GetFESpace()), u2(&op2.GetFESpace());
   u1 = 0.0; u2 = 0.0;
   op1.Solve(0.0, slip1, u1);
   op2.Solve(0.0, slip2, u2);

   Vector trac1, trac2;
   op1.ComputeTraction(u1, slip1, trac1);
   op2.ComputeTraction(u2, slip2, trac2);

   // Compare average traction magnitudes per fault DOF
   real_t avg_trac1 = trac1.Norml2() / std::sqrt(2.0 * nf1);
   real_t avg_trac2 = trac2.Norml2() / std::sqrt(2.0 * nf2);

   std::cout << "  Mesh1 avg |trac|/DOF = " << avg_trac1
             << ", Mesh2 avg |trac|/DOF = " << avg_trac2 << "\n";

   // The traction per DOF should scale correctly with element size.
   // The BR2 correction scales as μ*g/h, so for 2x larger mesh the
   // correction is ~0.5x. Combined with the average gradient term,
   // the ratio should stay bounded and reasonable.
   if (avg_trac1 > 1e-10 && avg_trac2 > 1e-10)
   {
      real_t ratio = avg_trac2 / avg_trac1;
      std::cout << "  Traction ratio (mesh2/mesh1) = " << ratio << "\n";
      // Allow wide tolerance since the solutions differ on different meshes
      TEST_ASSERT(ratio < 3.5,
                  "BR2 traction correction scaling is reasonable");
   }
}

// =============================================================================
// Test 17: Patch test — constant strain field gives exact stress with BR2
//
// Apply a displacement field u(x) = ε·x corresponding to constant strain.
// The traction from ComputeTraction should match the analytical stress.
// =============================================================================
void TestBR2PatchTestTraction()
{
   std::cout << "\n--- Test: BR2 Patch Test (Constant Strain Traction) ---\n";

   real_t Lx = 2.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(1, 1, 1, Lx, Ly, Lz);

   real_t lambda = 1.0, mu = 1.0;
   ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, 0.0, Lz, 2.0 * Ly,
                                      DGMethod::BR2);

   int nf = op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   // Set displacement to u = (ε_11 * x1, 0, 0) with ε_11 = 1
   // This gives strain ε = diag(1,0,0), stress σ = diag(λ+2μ, λ, λ)
   // On fault (normal = ±x1): T = σ · n = ((λ+2μ)·n1, 0, 0)
   FiniteElementSpace &fes = op.GetFESpace();
   GridFunction u(&fes);
   u = 0.0;

   // Get scalar FES info for byNODES ordering
   int ndof_scalar = fes.GetVSize() / 3;

   // Use DG_FECollection to get element DOFs
   DG_FECollection fec(1, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);

   for (int e = 0; e < mesh.GetNE(); e++)
   {
      const FiniteElement *fe = scalar_fes.GetFE(e);
      ElementTransformation *T = mesh.GetElementTransformation(e);
      Array<int> sdofs;
      scalar_fes.GetElementDofs(e, sdofs);
      for (int k = 0; k < fe->GetDof(); k++)
      {
         const IntegrationPoint &ip = fe->GetNodes().IntPoint(k);
         T->SetIntPoint(&ip);
         Vector x(3);
         T->Transform(ip, x);
         int dof_idx = sdofs[k];
         u(dof_idx) = x(0);                    // u_x = x1
         u(ndof_scalar + dof_idx) = 0.0;       // u_y = 0
         u(2 * ndof_scalar + dof_idx) = 0.0;   // u_z = 0
      }
   }

   // Zero slip — the displacement field represents continuous deformation
   Vector slip(2 * nf);
   slip = 0.0;

   Vector traction;
   op.ComputeTraction(u, slip, traction);

   // For a continuous field with zero slip, the DG penalty correction
   // should be zero ([[u]] = δ = 0), so traction comes from the
   // average gradient term only.
   // Expected: σ·n̂ where σ = diag(λ+2μ, λ, λ) and n̂ ≈ (±1,0,0)
   // → T = (±(λ+2μ), 0, 0)
   // Projected to fault frame: dip ≈ T·dip_dir, strike ≈ T·strike_dir
   // With n=(1,0,0), dip=(0,0,-1), strike=(0,1,0):
   //   T_dip = 0, T_strike = 0 (since T is along normal)
   // So both traction components should be approximately zero.

   real_t trac_norm = traction.Norml2();
   std::cout << "  Patch test traction norm: " << trac_norm << "\n";

   // For DG order 1 on a uniform mesh, the patch test should give
   // near-zero tangential traction (normal stress only)
   for (int i = 0; i < nf; i++)
   {
      std::cout << "    DOF " << i << ": dip=" << traction(2*i)
                << " strike=" << traction(2*i+1) << "\n";
   }

   // The tangential components should be small relative to (λ+2μ)
   real_t scale = lambda + 2.0 * mu;
   real_t rel_trac = trac_norm / (scale * std::sqrt(2.0 * nf));
   TEST_ASSERT(rel_trac < 0.1,
               "BR2 patch test: tangential traction small for uniaxial strain");
}

// =============================================================================
// Test 18: Slip sign convention — positive strike slip → negative traction
//
// Verifies that the BR2 traction has the correct sign: positive slip causes
// a stress drop (negative traction in the slip direction).
// =============================================================================
void TestBR2SlipSignConvention()
{
   std::cout << "\n--- Test: BR2 Slip Sign Convention ---\n";

   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;
   Mesh mesh = CreateTestMesh3D(2, 1, 1, Lx, Ly, Lz);

   real_t lambda = 1.0, mu = 1.0;
   ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, 0.0, Lz, 2.0 * Ly,
                                      DGMethod::BR2);

   int nf = op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   // Apply uniform positive strike slip
   Vector slip(2 * nf);
   slip = 0.0;
   for (int i = 0; i < nf; i++)
   {
      slip(2 * i + 1) = 1.0;  // positive strike slip
   }

   GridFunction u(&op.GetFESpace());
   u = 0.0;
   op.Solve(0.0, slip, u);

   Vector traction;
   op.ComputeTraction(u, slip, traction);

   // Average strike traction should be negative (stress drop)
   real_t avg_strike = 0.0;
   for (int i = 0; i < nf; i++)
   {
      avg_strike += traction(2 * i + 1);
   }
   avg_strike /= nf;

   std::cout << "  BR2 avg strike traction = " << avg_strike << "\n";
   TEST_ASSERT(avg_strike < 0.0,
               "BR2: positive strike slip → negative traction (stress drop)");

   // Also check dip direction: with pure strike slip, dip traction should be small
   real_t avg_dip = 0.0;
   for (int i = 0; i < nf; i++)
   {
      avg_dip += std::abs(traction(2 * i));
   }
   avg_dip /= nf;

   std::cout << "  BR2 avg |dip traction| = " << avg_dip << "\n";
   TEST_ASSERT(avg_dip < std::abs(avg_strike),
               "BR2: dip traction smaller than strike for pure strike slip");

   // Traction magnitude should be bounded
   real_t trac_max = traction.Normlinf();
   TEST_ASSERT(trac_max < 100.0, "BR2: traction is bounded");
}

// =============================================================================
// Main
// =============================================================================
// =============================================================================
// Test: Tag-based fault detection on Gmsh mesh with Physical Surface 100
// =============================================================================
void TestTagBasedFaultDetection()
{
   std::cout << "\n--- Test: Tag-Based Fault Detection ---\n";

   // Load the actual BP5 mesh which has Physical Surface 100 (fault)
   const std::string mesh_file = "bp5/mesh/bp5_1000m.msh";
   std::ifstream f(mesh_file);
   if (!f.good())
   {
      std::cout << "  (Skipped: " << mesh_file << " not found)\n";
      return;
   }
   f.close();

   Mesh mesh(mesh_file.c_str(), 1, 1);

   // Scale to meters (mesh is in km)
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *v = mesh.GetVertex(i);
      v[0] *= 1000.0;
      v[1] *= 1000.0;
      v[2] *= 1000.0;
   }
   mesh.SetAttributes();

   // Check that the mesh has attribute 100
   bool has_100 = false;
   for (int i = 0; i < mesh.bdr_attributes.Size(); i++)
   {
      if (mesh.bdr_attributes[i] == 100) { has_100 = true; break; }
   }
   TEST_ASSERT(has_100, "Gmsh mesh has boundary attribute 100 (fault)");

   // Count boundary elements with attr 100
   int bdr_100_count = 0;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      if (mesh.GetBdrAttribute(be) == 100) { bdr_100_count++; }
   }
   TEST_ASSERT(bdr_100_count > 0, "Mesh has fault boundary elements (attr 100)");
   std::cout << "  Fault boundary elements: " << bdr_100_count << "\n";

   // Build the operator (this calls SetupFaultInfo → BuildFaultTaggedFaces)
   real_t lambda = 32.04e9, mu = 32.04e9;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;
   ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, Vp, Wf, lf,
                                       DGMethod::IP);

   int nf = op.GetNumFaultDOFs();
   TEST_ASSERT(nf > 0, "Tag-based detection found fault DOFs");
   std::cout << "  Fault DOFs (tag-based): " << nf << "\n";

   // Verify no fault faces at z=0 or z=Wf boundaries
   Vector depths;
   op.GetFaultDepths(depths);
   real_t z_min = depths.Min();
   real_t z_max = depths.Max();
   std::cout << "  Fault depth range: [" << z_min << ", " << z_max << "] m\n";

   // z_min should be > 0 (no faces at surface z=0)
   // z_max should be < Wf (no faces at fault base z=40km)
   TEST_ASSERT(z_min > 100.0,
               "No fault faces at z=0 surface (z_min > 100m)");
   TEST_ASSERT(z_max < Wf - 100.0,
               "No fault faces at z=Wf base (z_max < Wf-100m)");

   // Tag-based detection should find fewer or equal DOFs compared to
   // total tagged boundary elements (some tagged faces may not correspond
   // to interior faces, e.g. at mesh boundaries)
   TEST_ASSERT(nf > 0 && nf <= bdr_100_count,
               "Tag detection DOFs <= total tagged boundary elements");
}

// =============================================================================
// Test: Tag-based detection excludes fault-boundary intersection faces
// =============================================================================
void TestTagExcludesBoundaryFaces()
{
   std::cout << "\n--- Test: Tag Detection Excludes Boundary Faces ---\n";

   const std::string mesh_file = "bp5/mesh/bp5_1000m.msh";
   std::ifstream f(mesh_file);
   if (!f.good())
   {
      std::cout << "  (Skipped: " << mesh_file << " not found)\n";
      return;
   }
   f.close();

   Mesh mesh(mesh_file.c_str(), 1, 1);
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *v = mesh.GetVertex(i);
      v[0] *= 1000.0;
      v[1] *= 1000.0;
      v[2] *= 1000.0;
   }
   mesh.SetAttributes();

   real_t lambda = 32.04e9, mu = 32.04e9;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;
   ElasticityDomainOperator<Mesh> op(mesh, 1, lambda, mu, Vp, Wf, lf,
                                       DGMethod::IP);

   // Get fault face coordinates
   Vector coords_x2, coords_x3;
   op.GetFaultCoords2D(coords_x2, coords_x3);

   int nf = op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault DOFs found)\n";
      return;
   }

   // Check that no fault face is at the boundary edges
   real_t y_min = coords_x2.Min();
   real_t y_max = coords_x2.Max();
   real_t z_min = coords_x3.Min();
   real_t z_max = coords_x3.Max();

   std::cout << "  y range: [" << y_min << ", " << y_max << "] m\n";
   std::cout << "  z range: [" << z_min << ", " << z_max << "] m\n";

   // The fault extends y in [-50km, +50km] and z in [0, 40km]
   // Tagged faces should NOT be at the exact boundary
   TEST_ASSERT(y_min > -lf/2.0 + 100.0,
               "No fault faces at y=-lf/2 edge");
   TEST_ASSERT(y_max < lf/2.0 - 100.0,
               "No fault faces at y=+lf/2 edge");
   TEST_ASSERT(z_min > 100.0,
               "No fault faces at z=0 edge");
   TEST_ASSERT(z_max < Wf - 100.0,
               "No fault faces at z=Wf edge");
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "  Elasticity Domain Operator Tests\n";
   std::cout << "========================================\n";

   TestConstruction();
   TestFaultDetection();
   TestFaultBasis();
   TestZeroSlipEquilibrium();
   TestTractionExtraction();
   TestStiffnessAssembly();
   TestFaultGeometry3D();
   TestFaultGeometry3DValues();
   TestNonZeroSlipIP();
   TestNonZeroSlipBR2();
   TestDirichletLoading();
   TestDirichletLoadingShearTraction();
   TestBR2vsIP();
   TestBR2Default();
   TestTractionWithPenaltyCorrection();
   TestBR2TractionCorrectionConsistency();
   TestBR2PatchTestTraction();
   TestBR2SlipSignConvention();

   // Tag-based fault detection tests
   TestTagBasedFaultDetection();
   TestTagExcludesBoundaryFaces();

   TEST_PRINT_RESULTS();

   return num_failed;
}
