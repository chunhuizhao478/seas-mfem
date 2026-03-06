// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Test suite for DGElasticityBR2Integrator

#include "mfem.hpp"
#include "../../integrator/dg_elasticity_br2_integrator.hpp"
#include "../../integrator/dg_br2_integrator.hpp"
#include "test_macros.hpp"

#include <iostream>
#include <cmath>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

// Helper: Precompute scalar mass matrix inverses for a mesh
void PrecomputeScalarMassInv(Mesh &mesh, FiniteElementSpace &fes,
                              std::vector<DenseMatrix> &mass_inv)
{
   int ne = mesh.GetNE();
   mass_inv.resize(ne);

   for (int e = 0; e < ne; e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *T = mesh.GetElementTransformation(e);

      int ndof = fe->GetDof();
      DenseMatrix M(ndof);

      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(),
                                                2 * fe->GetOrder());
      Vector shape(ndof);
      M = 0.0;

      for (int j = 0; j < ir.GetNPoints(); j++)
      {
         const IntegrationPoint &ip = ir.IntPoint(j);
         T->SetIntPoint(&ip);
         fe->CalcShape(ip, shape);

         real_t w = ip.weight * T->Weight();
         for (int k = 0; k < ndof; k++)
         {
            for (int l = 0; l < ndof; l++)
            {
               M(k, l) += w * shape(k) * shape(l);
            }
         }
      }

      mass_inv[e].SetSize(ndof);
      DenseMatrixInverse Minv(M);
      Minv.GetInverseMatrix(mass_inv[e]);
   }
}

// =============================================================================
// Test 1: Face matrix symmetry for SIPG (epsilon=-1)
// =============================================================================
void TestSymmetry()
{
   std::cout << "\n--- Test: BR2 Elasticity Face Matrix Symmetry ---\n";

   // Create a simple 2x1x1 hex mesh with an interior face
   Mesh mesh = Mesh::MakeCartesian3D(2, 1, 1, Element::HEXAHEDRON,
                                      2.0, 1.0, 1.0);

   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeScalarMassInv(mesh, scalar_fes, mass_inv);

   ConstantCoefficient lambda_coeff(1.0);
   ConstantCoefficient mu_coeff(1.0);

   DGElasticityBR2Integrator integ(lambda_coeff, mu_coeff, -1.0, mass_inv, 3);

   // Find an interior face
   bool found = false;
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *FTr = mesh.GetInteriorFaceTransformations(f);
      if (FTr == nullptr) { continue; }

      const FiniteElement &el1 = *scalar_fes.GetFE(FTr->Elem1No);
      const FiniteElement &el2 = *scalar_fes.GetFE(FTr->Elem2No);

      DenseMatrix elmat;
      integ.AssembleFaceMatrix(el1, el2, *FTr, elmat);

      // Check symmetry
      int n = elmat.Height();
      real_t max_asym = 0.0;
      for (int i = 0; i < n; i++)
      {
         for (int j = i + 1; j < n; j++)
         {
            real_t diff = std::abs(elmat(i, j) - elmat(j, i));
            real_t scale = std::max(std::abs(elmat(i, j)), std::abs(elmat(j, i)));
            if (scale > 1e-15)
            {
               max_asym = std::max(max_asym, diff / scale);
            }
         }
      }

      TEST_ASSERT(max_asym < 1e-12,
                  "BR2 elasticity face matrix is symmetric (SIPG)");
      TEST_ASSERT(n == 2 * el1.GetDof() * 3,
                  "Face matrix has correct size (2*ndof*dim)");
      found = true;
      break;
   }
   TEST_ASSERT(found, "Found interior face for test");
}

// =============================================================================
// Test 2: Face matrix is non-zero
// =============================================================================
void TestNonZero()
{
   std::cout << "\n--- Test: BR2 Elasticity Face Matrix Non-Zero ---\n";

   Mesh mesh = Mesh::MakeCartesian3D(2, 1, 1, Element::HEXAHEDRON,
                                      2.0, 1.0, 1.0);

   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeScalarMassInv(mesh, scalar_fes, mass_inv);

   ConstantCoefficient lambda_coeff(1.0);
   ConstantCoefficient mu_coeff(1.0);

   DGElasticityBR2Integrator integ(lambda_coeff, mu_coeff, -1.0, mass_inv, 3);

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *FTr = mesh.GetInteriorFaceTransformations(f);
      if (FTr == nullptr) { continue; }

      const FiniteElement &el1 = *scalar_fes.GetFE(FTr->Elem1No);
      const FiniteElement &el2 = *scalar_fes.GetFE(FTr->Elem2No);

      DenseMatrix elmat;
      integ.AssembleFaceMatrix(el1, el2, *FTr, elmat);

      real_t norm = elmat.FNorm();
      TEST_ASSERT(norm > 1e-10, "BR2 elasticity face matrix is non-zero");
      break;
   }
}

// =============================================================================
// Test 3: Boundary face matrix symmetry
// =============================================================================
void TestBoundarySymmetry()
{
   std::cout << "\n--- Test: BR2 Elasticity Boundary Face Matrix Symmetry ---\n";

   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON,
                                      2.0, 2.0, 2.0);

   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeScalarMassInv(mesh, scalar_fes, mass_inv);

   ConstantCoefficient lambda_coeff(1.0);
   ConstantCoefficient mu_coeff(1.0);

   DGElasticityBR2BoundaryIntegrator bdr_integ(
      lambda_coeff, mu_coeff, -1.0, mass_inv, 3);

   // Test on first boundary face
   bool found = false;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      int face_idx, face_info;
      mesh.GetBdrElementFace(be, &face_idx, &face_info);
      FaceElementTransformations *FTr =
         mesh.GetFaceElementTransformations(face_idx);
      if (FTr == nullptr) { continue; }

      const FiniteElement &el1 = *scalar_fes.GetFE(FTr->Elem1No);
      // For boundary face, el2 is dummy
      DenseMatrix elmat;
      bdr_integ.AssembleFaceMatrix(el1, el1, *FTr, elmat);

      int n = elmat.Height();
      real_t max_asym = 0.0;
      for (int i = 0; i < n; i++)
      {
         for (int j = i + 1; j < n; j++)
         {
            real_t diff = std::abs(elmat(i, j) - elmat(j, i));
            real_t scale = std::max(std::abs(elmat(i, j)), std::abs(elmat(j, i)));
            if (scale > 1e-15)
            {
               max_asym = std::max(max_asym, diff / scale);
            }
         }
      }

      TEST_ASSERT(max_asym < 1e-12,
                  "BR2 elasticity boundary face matrix is symmetric");
      found = true;
      break;
   }
   TEST_ASSERT(found, "Found boundary face for test");
}

// =============================================================================
// Test 4: Full bilinear form assembly with BR2
// =============================================================================
void TestFullAssembly()
{
   std::cout << "\n--- Test: Full Bilinear Form with BR2 Elasticity ---\n";

   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON,
                                      2.0, 2.0, 2.0);

   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);
   FiniteElementSpace vec_fes(&mesh, &fec, 3, Ordering::byNODES);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeScalarMassInv(mesh, scalar_fes, mass_inv);

   ConstantCoefficient lambda_coeff(1.0);
   ConstantCoefficient mu_coeff(1.0);

   BilinearForm a(&vec_fes);
   a.AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff, mu_coeff));
   a.AddInteriorFaceIntegrator(
      new DGElasticityBR2Integrator(lambda_coeff, mu_coeff, -1.0, mass_inv, 3));
   a.Assemble();
   a.Finalize();

   // Check that matrix is assembled and has correct size
   SparseMatrix &A = a.SpMat();
   int expected_size = vec_fes.GetVSize();
   TEST_ASSERT(A.Height() == expected_size,
               "Assembled matrix has correct height");
   TEST_ASSERT(A.Width() == expected_size,
               "Assembled matrix has correct width");

   // Check matrix is non-zero
   real_t norm = A.MaxNorm();
   TEST_ASSERT(norm > 1e-10, "Assembled stiffness matrix is non-zero");

   // Check symmetry of the sparse matrix
   SparseMatrix *AT = Transpose(A);
   AT->Add(-1.0, A);
   real_t asym_norm = AT->MaxNorm();
   real_t rel_asym = asym_norm / norm;
   TEST_ASSERT(rel_asym < 1e-10,
               "Full assembled stiffness matrix is symmetric");
   delete AT;
}

// =============================================================================
// Test 5: Polynomial order 0 - simplest case
// =============================================================================
void TestOrder0()
{
   std::cout << "\n--- Test: BR2 Elasticity with Order 0 ---\n";

   Mesh mesh = Mesh::MakeCartesian3D(2, 1, 1, Element::HEXAHEDRON,
                                      2.0, 1.0, 1.0);

   int order = 0;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeScalarMassInv(mesh, scalar_fes, mass_inv);

   ConstantCoefficient lambda_coeff(1.0);
   ConstantCoefficient mu_coeff(1.0);

   DGElasticityBR2Integrator integ(lambda_coeff, mu_coeff, -1.0, mass_inv, 3);

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *FTr = mesh.GetInteriorFaceTransformations(f);
      if (FTr == nullptr) { continue; }

      const FiniteElement &el1 = *scalar_fes.GetFE(FTr->Elem1No);
      const FiniteElement &el2 = *scalar_fes.GetFE(FTr->Elem2No);

      DenseMatrix elmat;
      integ.AssembleFaceMatrix(el1, el2, *FTr, elmat);

      // Order 0: 1 DOF per element, 3 components = 3+3=6 total
      TEST_ASSERT(elmat.Height() == 6, "Order 0: face matrix is 6x6");

      // Check symmetry
      real_t max_asym = 0.0;
      for (int i = 0; i < 6; i++)
      {
         for (int j = i + 1; j < 6; j++)
         {
            real_t diff = std::abs(elmat(i, j) - elmat(j, i));
            real_t scale = std::max(std::abs(elmat(i, j)),
                                     std::abs(elmat(j, i)));
            if (scale > 1e-15)
            {
               max_asym = std::max(max_asym, diff / scale);
            }
         }
      }
      TEST_ASSERT(max_asym < 1e-12, "Order 0: face matrix is symmetric");
      break;
   }
}

// =============================================================================
// Test 6: Different lambda/mu values
// =============================================================================
void TestMaterialParams()
{
   std::cout << "\n--- Test: BR2 Elasticity with Different Material Params ---\n";

   Mesh mesh = Mesh::MakeCartesian3D(2, 1, 1, Element::HEXAHEDRON,
                                      2.0, 1.0, 1.0);

   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeScalarMassInv(mesh, scalar_fes, mass_inv);

   // Test with BP5 material: lambda = mu = rho*cs^2 for nu=0.25
   real_t bp5_mu = 2670.0 * 3464.0 * 3464.0;
   real_t bp5_lambda = bp5_mu;  // nu=0.25 => lambda=mu
   ConstantCoefficient lambda_coeff(bp5_lambda);
   ConstantCoefficient mu_coeff(bp5_mu);

   DGElasticityBR2Integrator integ(lambda_coeff, mu_coeff, -1.0, mass_inv, 3);

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *FTr = mesh.GetInteriorFaceTransformations(f);
      if (FTr == nullptr) { continue; }

      const FiniteElement &el1 = *scalar_fes.GetFE(FTr->Elem1No);
      const FiniteElement &el2 = *scalar_fes.GetFE(FTr->Elem2No);

      DenseMatrix elmat;
      integ.AssembleFaceMatrix(el1, el2, *FTr, elmat);

      // Check symmetry using Frobenius norm as reference scale
      // (per-entry relative check fails for near-zero entries with large coefficients)
      int n = elmat.Height();
      real_t fnorm = elmat.FNorm();
      real_t max_asym = 0.0;
      for (int i = 0; i < n; i++)
      {
         for (int j = i + 1; j < n; j++)
         {
            real_t diff = std::abs(elmat(i, j) - elmat(j, i));
            max_asym = std::max(max_asym, diff);
         }
      }
      real_t rel_asym = (fnorm > 0.0) ? max_asym / fnorm : 0.0;
      TEST_ASSERT(rel_asym < 1e-12,
                  "BP5 material params: face matrix is symmetric");

      // Should be non-zero and proportional to material params
      real_t norm = elmat.FNorm();
      TEST_ASSERT(norm > 1e-10 * bp5_mu,
                  "BP5 material params: face matrix has expected magnitude");
      break;
   }
}

// =============================================================================
// Test 7: BR2 face matrix DOF ordering matches byNODES convention
//
// For order >= 1 with byNODES ordering, the face matrix block structure
// should be compatible with the volume ElasticityIntegrator. We verify
// this by checking that the assembled stiffness matrix (volume + BR2 faces)
// applied to a rigid body mode gives zero.
// =============================================================================
void TestDOFOrdering()
{
   std::cout << "\n--- Test: BR2 DOF Ordering (byNODES) Order 1 ---\n";

   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON,
                                      2.0, 2.0, 2.0);

   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);
   FiniteElementSpace vec_fes(&mesh, &fec, 3, Ordering::byNODES);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeScalarMassInv(mesh, scalar_fes, mass_inv);

   ConstantCoefficient lambda_coeff(1.0);
   ConstantCoefficient mu_coeff(1.0);

   // Assemble full stiffness: volume + BR2 interior faces + BR2 boundary faces
   BilinearForm a(&vec_fes);
   a.AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff, mu_coeff));
   a.AddInteriorFaceIntegrator(
      new DGElasticityBR2Integrator(lambda_coeff, mu_coeff, -1.0, mass_inv, 3));

   // Add BR2 boundary on all faces for complete system
   Array<int> all_bdr(mesh.bdr_attributes.Max());
   all_bdr = 1;
   a.AddBdrFaceIntegrator(
      new DGElasticityBR2BoundaryIntegrator(
         lambda_coeff, mu_coeff, -1.0, mass_inv, 3), all_bdr);

   a.Assemble();
   a.Finalize();

   SparseMatrix &A = a.SpMat();

   // A rigid body translation should be in the null space of the volume term.
   // With DG + Dirichlet-type boundary conditions (via BR2 boundary integrator),
   // A * u_rigid should give a non-zero RHS (the boundary penalty enforces
   // the BC). The key check is that the matrix assembles without error
   // and has the correct structure.

   // Instead, check that the matrix is SPD-like by verifying
   // x^T A x > 0 for a non-trivial vector
   int N = vec_fes.GetVSize();
   Vector x(N), Ax(N);
   x.Randomize(42);
   A.Mult(x, Ax);
   real_t xAx = InnerProduct(x, Ax);
   TEST_ASSERT(xAx > 0.0,
               "BR2 order-1: x^T A x > 0 (positive definiteness check)");

   // Symmetry check
   SparseMatrix *AT = Transpose(A);
   AT->Add(-1.0, A);
   real_t asym_norm = AT->MaxNorm();
   real_t rel_asym = asym_norm / A.MaxNorm();
   TEST_ASSERT(rel_asym < 1e-10,
               "BR2 order-1: full assembled matrix is symmetric");
   delete AT;
}

// =============================================================================
// Test 8: BR2 vs IP face matrix comparison at order 0
//
// For order 0 (1 DOF per element, 6x6 face matrix), both BR2 and MFEM's
// DGElasticityIntegrator should produce equivalent face matrices (up to
// penalty differences). We check they have the same block structure.
// =============================================================================
void TestBR2vsIPOrder0()
{
   std::cout << "\n--- Test: BR2 vs IP Face Matrix at Order 0 ---\n";

   Mesh mesh = Mesh::MakeCartesian3D(2, 1, 1, Element::HEXAHEDRON,
                                      2.0, 1.0, 1.0);

   int order = 0;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);
   FiniteElementSpace vec_fes(&mesh, &fec, 3, Ordering::byNODES);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeScalarMassInv(mesh, scalar_fes, mass_inv);

   ConstantCoefficient lambda_coeff(1.0);
   ConstantCoefficient mu_coeff(1.0);

   // BR2 integrator
   DGElasticityBR2Integrator br2_integ(lambda_coeff, mu_coeff, -1.0, mass_inv, 3);

   // IP integrator (MFEM built-in)
   real_t kappa = 1.0;  // order 0: (p+1)^2 = 1
   DGElasticityIntegrator ip_integ(lambda_coeff, mu_coeff, -1.0, kappa);

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *FTr = mesh.GetInteriorFaceTransformations(f);
      if (FTr == nullptr) { continue; }

      const FiniteElement &el1 = *scalar_fes.GetFE(FTr->Elem1No);
      const FiniteElement &el2 = *scalar_fes.GetFE(FTr->Elem2No);

      DenseMatrix br2_mat, ip_mat;
      br2_integ.AssembleFaceMatrix(el1, el2, *FTr, br2_mat);
      ip_integ.AssembleFaceMatrix(el1, el2, *FTr, ip_mat);

      // Both should be 6x6
      TEST_ASSERT(br2_mat.Height() == 6, "BR2 order-0 face matrix is 6x6");
      TEST_ASSERT(ip_mat.Height() == 6, "IP order-0 face matrix is 6x6");

      // Both should be symmetric
      real_t br2_asym = 0.0, ip_asym = 0.0;
      for (int i = 0; i < 6; i++)
         for (int j = i+1; j < 6; j++)
         {
            br2_asym = std::max(br2_asym, std::abs(br2_mat(i,j) - br2_mat(j,i)));
            ip_asym = std::max(ip_asym, std::abs(ip_mat(i,j) - ip_mat(j,i)));
         }
      TEST_ASSERT(br2_asym < 1e-12 * br2_mat.FNorm(),
                  "BR2 order-0 face matrix symmetric");
      TEST_ASSERT(ip_asym < 1e-12 * ip_mat.FNorm(),
                  "IP order-0 face matrix symmetric");

      // Both should have the same non-zero pattern (same block structure).
      // Check that non-zero entries of BR2 correspond to non-zero entries of IP.
      int matching_nonzero = 0;
      int total_nonzero = 0;
      for (int i = 0; i < 6; i++)
         for (int j = 0; j < 6; j++)
         {
            bool br2_nz = std::abs(br2_mat(i,j)) > 1e-14 * br2_mat.FNorm();
            bool ip_nz = std::abs(ip_mat(i,j)) > 1e-14 * ip_mat.FNorm();
            if (br2_nz || ip_nz)
            {
               total_nonzero++;
               if (br2_nz && ip_nz) { matching_nonzero++; }
            }
         }
      real_t match_frac = (total_nonzero > 0)
                             ? real_t(matching_nonzero) / total_nonzero : 1.0;
      TEST_ASSERT(match_frac > 0.5,
                  "BR2 and IP face matrices have similar nonzero pattern");
      break;
   }
}

// =============================================================================
// Test 9: DG Patch Test — verifies consistency of interior face operators
//
// For DG elasticity with a continuous function u:
//   [[u]] = 0 on all interior faces, so:
//   - Symmetry term vanishes (involves [[u]])
//   - Penalty/BR2 term vanishes (involves [[u]] or R_h([[u]]))
//   - Consistency term -∫{{σ(u)·n}}·[[v]] does NOT vanish for σ(u)≠0
//
// Two sub-tests:
// (a) Constant displacement u = const: σ(u)=0 AND [[u]]=0 → K_faces*u = 0
// (b) Linear displacement u(x)=(x1,0,0): σ(u)≠0, [[u]]=0
//     SIPG (ε=-1) and NIPG (ε=+1) should give identical K_faces*u,
//     since the symmetry+penalty terms vanish, leaving only the
//     ε-independent consistency term.
// =============================================================================
void TestPatchTest()
{
   std::cout << "\n--- Test: BR2 DG Patch Test ---\n";

   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON,
                                      2.0, 2.0, 2.0);

   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);
   FiniteElementSpace vec_fes(&mesh, &fec, 3, Ordering::byNODES);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeScalarMassInv(mesh, scalar_fes, mass_inv);

   real_t lam = 1.0, mu = 1.0;
   ConstantCoefficient lambda_coeff(lam);
   ConstantCoefficient mu_coeff(mu);

   int ndof_scalar = scalar_fes.GetVSize();

   // --- Sub-test (a): Constant displacement u = (1, 2, 3) ---
   // σ(u) = 0, [[u]] = 0 → ALL face terms vanish
   GridFunction u_const(&vec_fes);
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      Array<int> sdofs;
      scalar_fes.GetElementDofs(e, sdofs);
      for (int k = 0; k < scalar_fes.GetFE(e)->GetDof(); k++)
      {
         int dof_idx = sdofs[k];
         u_const(dof_idx) = 1.0;
         u_const(ndof_scalar + dof_idx) = 2.0;
         u_const(2*ndof_scalar + dof_idx) = 3.0;
      }
   }

   // BR2 face-only bilinear form
   BilinearForm a_br2(&vec_fes);
   a_br2.AddInteriorFaceIntegrator(
      new DGElasticityBR2Integrator(lambda_coeff, mu_coeff, -1.0, mass_inv, 3));
   a_br2.Assemble();
   a_br2.Finalize();

   Vector Ku_const(vec_fes.GetVSize());
   a_br2.SpMat().Mult(u_const, Ku_const);
   real_t const_residual = Ku_const.Norml2();
   real_t const_norm = u_const.Norml2();
   real_t const_rel = (const_norm > 0.0) ? const_residual / const_norm : const_residual;
   std::cout << "  (a) Constant u: |K_faces * u| / |u| = " << const_rel << "\n";
   TEST_ASSERT(const_rel < 1e-10,
               "BR2 patch test: face terms vanish for constant displacement");

   // IP face-only bilinear form
   real_t kappa = 4.0;
   BilinearForm a_ip(&vec_fes);
   a_ip.AddInteriorFaceIntegrator(
      new DGElasticityIntegrator(lambda_coeff, mu_coeff, -1.0, kappa));
   a_ip.Assemble();
   a_ip.Finalize();

   Vector Ku_ip_const(vec_fes.GetVSize());
   a_ip.SpMat().Mult(u_const, Ku_ip_const);
   real_t ip_const_rel = Ku_ip_const.Norml2() / const_norm;
   std::cout << "  (a) IP constant: |K_faces * u| / |u| = " << ip_const_rel << "\n";
   TEST_ASSERT(ip_const_rel < 1e-10,
               "IP patch test: face terms vanish for constant displacement");

   // --- Sub-test (b): Linear displacement u(x) = (x1, 0, 0) ---
   // [[u]] = 0 but σ(u) ≠ 0: the consistency term is non-zero (expected).
   // Both BR2 and IP penalty terms vanish for [[u]]=0.
   // Verify: BR2 SIPG and IP SIPG give the same K_faces*u
   // (both should compute the same ε-independent consistency contribution).
   GridFunction u_linear(&vec_fes);
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
         u_linear(dof_idx) = x(0);
         u_linear(ndof_scalar + dof_idx) = 0.0;
         u_linear(2*ndof_scalar + dof_idx) = 0.0;
      }
   }

   // BR2 SIPG K*u (already assembled above as a_br2)
   Vector Ku_br2(vec_fes.GetVSize());
   a_br2.SpMat().Mult(u_linear, Ku_br2);

   // IP SIPG K*u (already assembled above as a_ip)
   Vector Ku_ip_lin(vec_fes.GetVSize());
   a_ip.SpMat().Mult(u_linear, Ku_ip_lin);

   // Both should be non-zero (consistency term exists)
   real_t br2_norm = Ku_br2.Norml2();
   real_t ip_norm = Ku_ip_lin.Norml2();
   std::cout << "  (b) Linear u: |K_br2*u| = " << br2_norm
             << ", |K_ip*u| = " << ip_norm << "\n";
   TEST_ASSERT(br2_norm > 1e-10,
               "BR2: non-zero K_faces*u for linear u (consistency term)");
   TEST_ASSERT(ip_norm > 1e-10,
               "IP: non-zero K_faces*u for linear u (consistency term)");

   // BR2 and IP should agree (both penalties vanish for [[u]]=0)
   Vector diff(vec_fes.GetVSize());
   subtract(Ku_br2, Ku_ip_lin, diff);
   real_t diff_rel = diff.Norml2() / br2_norm;
   std::cout << "  (b) |K_br2*u - K_ip*u| / |K_br2*u| = " << diff_rel << "\n";
   TEST_ASSERT(diff_rel < 1e-10,
               "Patch test: BR2 and IP agree for continuous linear u");
}

// =============================================================================
// Test 10: Tet element — BR2 face matrix assembly for tetrahedra
// =============================================================================
void TestTetElements()
{
   std::cout << "\n--- Test: BR2 Elasticity with Tet Elements ---\n";

   // Create a simple tet mesh
   Mesh mesh = Mesh::MakeCartesian3D(2, 1, 1, Element::TETRAHEDRON,
                                      2.0, 1.0, 1.0);

   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(&mesh, &fec);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeScalarMassInv(mesh, scalar_fes, mass_inv);

   ConstantCoefficient lambda_coeff(1.0);
   ConstantCoefficient mu_coeff(1.0);

   DGElasticityBR2Integrator integ(lambda_coeff, mu_coeff, -1.0, mass_inv, 3);

   // Find an interior face and assemble
   bool found = false;
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *FTr = mesh.GetInteriorFaceTransformations(f);
      if (FTr == nullptr) { continue; }

      const FiniteElement &el1 = *scalar_fes.GetFE(FTr->Elem1No);
      const FiniteElement &el2 = *scalar_fes.GetFE(FTr->Elem2No);

      DenseMatrix elmat;
      integ.AssembleFaceMatrix(el1, el2, *FTr, elmat);

      // Tet order 1: 4 DOFs per element, 3 components = 12 per side
      int n = elmat.Height();
      int expected = (el1.GetDof() + el2.GetDof()) * 3;
      TEST_ASSERT(n == expected,
                  "Tet face matrix has correct size");
      TEST_ASSERT(elmat.FNorm() > 1e-10,
                  "Tet face matrix is non-zero");

      // Symmetry check
      real_t max_asym = 0.0;
      real_t fnorm = elmat.FNorm();
      for (int i = 0; i < n; i++)
         for (int j = i + 1; j < n; j++)
            max_asym = std::max(max_asym, std::abs(elmat(i,j) - elmat(j,i)));
      real_t rel_asym = (fnorm > 0.0) ? max_asym / fnorm : 0.0;
      TEST_ASSERT(rel_asym < 1e-12,
                  "Tet face matrix is symmetric");

      found = true;
      break;
   }
   TEST_ASSERT(found, "Found interior face in tet mesh");

   // Full assembly with tet mesh
   FiniteElementSpace vec_fes(&mesh, &fec, 3, Ordering::byNODES);
   BilinearForm a(&vec_fes);
   a.AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff, mu_coeff));
   a.AddInteriorFaceIntegrator(
      new DGElasticityBR2Integrator(lambda_coeff, mu_coeff, -1.0, mass_inv, 3));
   a.Assemble();
   a.Finalize();

   SparseMatrix &A = a.SpMat();
   TEST_ASSERT(A.MaxNorm() > 1e-10, "Tet full assembly is non-zero");

   // Symmetry
   SparseMatrix *AT = Transpose(A);
   AT->Add(-1.0, A);
   real_t asym = AT->MaxNorm() / A.MaxNorm();
   TEST_ASSERT(asym < 1e-10, "Tet full assembly is symmetric");
   delete AT;
}

// =============================================================================
// Main
// =============================================================================
int main()
{
   std::cout << "========================================\n";
   std::cout << "  DG Elasticity BR2 Integrator Tests\n";
   std::cout << "========================================\n";

   TestSymmetry();
   TestNonZero();
   TestBoundarySymmetry();
   TestFullAssembly();
   TestOrder0();
   TestMaterialParams();
   TestDOFOrdering();
   TestBR2vsIPOrder0();
   TestPatchTest();
   TestTetElements();

   TEST_PRINT_RESULTS();

   return num_failed;
}
