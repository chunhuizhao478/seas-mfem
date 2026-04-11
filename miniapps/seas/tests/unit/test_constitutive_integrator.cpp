// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// KEY equivalence test (plan item 3c-ii):
// Verify that DG integrators constructed with old (Coefficient&) and new
// (ConstitutiveModel&) paths produce identical element matrices.
//
// Run: ./seas_test_constitutive_integrator

#include "mfem.hpp"
#include "../../constitutive/linear_elastic.hpp"
#include "../../integrator/dg_elasticity_br2_integrator.hpp"
#include "../../integrator/dg_elasticity_ip_combined_integrator.hpp"
#include "test_macros.hpp"

#include <cmath>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

/// Create a small 2x1x1 hex mesh for face matrix assembly tests.
/// Returns a mesh with at least one interior face.
std::unique_ptr<Mesh> MakeSmallMesh()
{
   auto mesh = std::make_unique<Mesh>(
      Mesh::MakeCartesian3D(2, 1, 1, Element::HEXAHEDRON, 2.0, 1.0, 1.0));
   return mesh;
}

/// Find the first interior face in the mesh.
int FindInteriorFace(Mesh &mesh)
{
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      auto *tr = mesh.GetInteriorFaceTransformations(f);
      if (tr) { return f; }
   }
   MFEM_ABORT("No interior face found");
   return -1;
}

/// Compare two dense matrices entry-by-entry. Returns max absolute difference.
real_t MaxMatrixDiff(const DenseMatrix &A, const DenseMatrix &B)
{
   MFEM_VERIFY(A.Height() == B.Height() && A.Width() == B.Width(),
               "Matrix dimension mismatch");
   real_t max_diff = 0.0;
   for (int i = 0; i < A.Height(); i++)
   {
      for (int j = 0; j < A.Width(); j++)
      {
         real_t diff = std::abs(A(i, j) - B(i, j));
         if (diff > max_diff) { max_diff = diff; }
      }
   }
   return max_diff;
}

/// Precompute scalar element mass matrix inverses for BR2.
void PrecomputeMassInv(Mesh &mesh, FiniteElementSpace &fes,
                       std::vector<DenseMatrix> &mass_inv)
{
   mass_inv.resize(mesh.GetNE());
   MassIntegrator mi;
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      const FiniteElement &fe = *fes.GetFE(e);
      ElementTransformation &eltrans = *fes.GetElementTransformation(e);
      DenseMatrix M(fe.GetDof());
      mi.AssembleElementMatrix(fe, eltrans, M);
      mass_inv[e].SetSize(fe.GetDof());
      DenseMatrixInverse Minv(M);
      Minv.GetInverseMatrix(mass_inv[e]);
   }
}

void TestBR2Equivalence()
{
   std::cout << "\n=== Test: BR2 integrator old vs new path ===\n";

   auto mesh = MakeSmallMesh();
   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(mesh.get(), &fec);

   // Precompute mass inverses
   std::vector<DenseMatrix> mass_inv;
   PrecomputeMassInv(*mesh, scalar_fes, mass_inv);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t epsilon = -1.0;

   // Old path: Coefficient references
   ConstantCoefficient lam_coeff(lambda), mu_coeff(mu);
   DGElasticityBR2Integrator old_integ(lam_coeff, mu_coeff, epsilon, mass_inv, 3);

   // New path: ConstitutiveModel
   LinearElastic model(lambda, mu);
   DGElasticityBR2Integrator new_integ(model, epsilon, mass_inv, 3);

   // Find interior face and assemble
   int face_idx = FindInteriorFace(*mesh);
   auto *tr = mesh->GetInteriorFaceTransformations(face_idx);

   const FiniteElement &fe1 = *scalar_fes.GetFE(tr->Elem1No);
   const FiniteElement &fe2 = *scalar_fes.GetFE(tr->Elem2No);

   DenseMatrix elmat_old, elmat_new;
   old_integ.AssembleFaceMatrix(fe1, fe2, *tr, elmat_old);
   new_integ.AssembleFaceMatrix(fe1, fe2, *tr, elmat_new);

   TEST_ASSERT(elmat_old.Height() == elmat_new.Height(),
               "BR2 matrix dimensions match (height)");
   TEST_ASSERT(elmat_old.Width() == elmat_new.Width(),
               "BR2 matrix dimensions match (width)");

   real_t max_diff = MaxMatrixDiff(elmat_old, elmat_new);
   std::cout << "  BR2 max entry diff: " << max_diff << "\n";
   TEST_ASSERT(max_diff == 0.0, "BR2 old vs new: bit-for-bit identical");
}

void TestBR2BoundaryEquivalence()
{
   std::cout << "\n=== Test: BR2 boundary integrator old vs new path ===\n";

   auto mesh = MakeSmallMesh();
   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(mesh.get(), &fec);

   std::vector<DenseMatrix> mass_inv;
   PrecomputeMassInv(*mesh, scalar_fes, mass_inv);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t epsilon = -1.0;

   ConstantCoefficient lam_coeff(lambda), mu_coeff(mu);
   DGElasticityBR2BoundaryIntegrator old_integ(
      lam_coeff, mu_coeff, epsilon, mass_inv, 3);

   LinearElastic model(lambda, mu);
   DGElasticityBR2BoundaryIntegrator new_integ(
      model, epsilon, mass_inv, 3);

   // Use a boundary face
   const FiniteElement &fe1 = *scalar_fes.GetFE(0);
   // For boundary faces, el2 is a dummy — get the first boundary face
   for (int bf = 0; bf < mesh->GetNBE(); bf++)
   {
      int f = mesh->GetBdrElementFaceIndex(bf);
      auto *tr = mesh->GetFaceElementTransformations(f);
      if (tr && tr->Elem2No < 0)
      {
         const FiniteElement &bfe = *scalar_fes.GetFE(tr->Elem1No);
         // el2 = dummy for boundary
         DenseMatrix elmat_old, elmat_new;
         old_integ.AssembleFaceMatrix(bfe, bfe, *tr, elmat_old);
         new_integ.AssembleFaceMatrix(bfe, bfe, *tr, elmat_new);

         real_t max_diff = MaxMatrixDiff(elmat_old, elmat_new);
         std::cout << "  BR2 bdr max entry diff: " << max_diff << "\n";
         TEST_ASSERT(max_diff == 0.0,
                     "BR2 boundary old vs new: bit-for-bit identical");
         return;
      }
   }
   TEST_ASSERT(false, "BR2 boundary: found boundary face");
}

void TestIPCombinedEquivalence()
{
   std::cout << "\n=== Test: IP combined integrator old vs new path ===\n";

   auto mesh = MakeSmallMesh();
   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(mesh.get(), &fec);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t epsilon = -1.0;
   real_t penalty_factor = 1.0;

   // Old path
   ConstantCoefficient lam_coeff(lambda), mu_coeff(mu);
   DGElasticityIPCombinedIntegrator old_integ(
      lam_coeff, mu_coeff, 3, epsilon, penalty_factor);

   // New path
   LinearElastic model(lambda, mu);
   DGElasticityIPCombinedIntegrator new_integ(
      model, 3, epsilon, penalty_factor);

   int face_idx = FindInteriorFace(*mesh);
   auto *tr = mesh->GetInteriorFaceTransformations(face_idx);

   const FiniteElement &fe1 = *scalar_fes.GetFE(tr->Elem1No);
   const FiniteElement &fe2 = *scalar_fes.GetFE(tr->Elem2No);

   DenseMatrix elmat_old, elmat_new;
   old_integ.AssembleFaceMatrix(fe1, fe2, *tr, elmat_old);
   new_integ.AssembleFaceMatrix(fe1, fe2, *tr, elmat_new);

   TEST_ASSERT(elmat_old.Height() == elmat_new.Height(),
               "IP matrix dimensions match (height)");
   TEST_ASSERT(elmat_old.Width() == elmat_new.Width(),
               "IP matrix dimensions match (width)");

   real_t max_diff = MaxMatrixDiff(elmat_old, elmat_new);
   std::cout << "  IP combined max entry diff: " << max_diff << "\n";
   TEST_ASSERT(max_diff == 0.0, "IP combined old vs new: bit-for-bit identical");
}

void TestIPCombinedBoundaryEquivalence()
{
   std::cout << "\n=== Test: IP combined boundary old vs new path ===\n";

   auto mesh = MakeSmallMesh();
   int order = 1;
   DG_FECollection fec(order, 3, BasisType::GaussLobatto);
   FiniteElementSpace scalar_fes(mesh.get(), &fec);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t epsilon = -1.0;

   ConstantCoefficient lam_coeff(lambda), mu_coeff(mu);
   DGElasticityIPCombinedIntegrator old_integ(lam_coeff, mu_coeff, 3, epsilon);
   LinearElastic model(lambda, mu);
   DGElasticityIPCombinedIntegrator new_integ(model, 3, epsilon);

   for (int bf = 0; bf < mesh->GetNBE(); bf++)
   {
      int f = mesh->GetBdrElementFaceIndex(bf);
      auto *tr = mesh->GetFaceElementTransformations(f);
      if (tr && tr->Elem2No < 0)
      {
         const FiniteElement &bfe = *scalar_fes.GetFE(tr->Elem1No);
         DenseMatrix elmat_old, elmat_new;
         old_integ.AssembleFaceMatrix(bfe, bfe, *tr, elmat_old);
         new_integ.AssembleFaceMatrix(bfe, bfe, *tr, elmat_new);

         real_t max_diff = MaxMatrixDiff(elmat_old, elmat_new);
         std::cout << "  IP bdr max entry diff: " << max_diff << "\n";
         TEST_ASSERT(max_diff == 0.0,
                     "IP combined boundary old vs new: bit-for-bit identical");
         return;
      }
   }
   TEST_ASSERT(false, "IP bdr: found boundary face");
}

int main(int argc, char *argv[])
{
   TestBR2Equivalence();
   TestBR2BoundaryEquivalence();
   TestIPCombinedEquivalence();
   TestIPCombinedBoundaryEquivalence();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
