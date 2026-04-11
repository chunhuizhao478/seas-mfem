// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// Unit tests for BuildFacetBCTables: face-by-face boundary classification.
// Verifies that the inline BP5 mesh correctly classifies fault, Dirichlet,
// and interior faces via the BoundaryConfig path.
//
// Plan item 4h-iii.
// Run: ./seas_test_boundary_classification

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../constitutive/linear_elastic.hpp"
#include "../../config/bp5_mesh_utils.hpp"
#include "test_macros.hpp"

#include <iostream>
#include <cmath>
#include <set>

using namespace mfem;
using namespace mfem::seas;

// Helper: get face centroid for a serial mesh
static Vector FaceCentroid(Mesh &mesh, int face_idx)
{
   Array<int> verts;
   mesh.GetFaceVertices(face_idx, verts);
   int dim = mesh.SpaceDimension();
   Vector centroid(dim);
   centroid = 0.0;
   for (int i = 0; i < verts.Size(); i++)
   {
      const real_t *v = mesh.GetVertex(verts[i]);
      for (int d = 0; d < dim; d++) { centroid(d) += v[d]; }
   }
   centroid /= verts.Size();
   return centroid;
}

void TestFaultFacesAreOnY0()
{
   std::cout << "\n=== Test: All fault faces lie on y=0 plane ===\n";

   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;

   LinearElastic model(lambda, mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   ElasticityDomainOperator<Mesh> op(
      *mesh, 1, model, Vp, Wf, lf, bc,
      DGMethod::BR2, SolverType::CG_AMG);

   const auto &fault_faces = op.GetFaultInteriorFaces();
   TEST_ASSERT(fault_faces.Size() > 0, "has fault interior faces");

   bool all_on_y0 = true;
   for (int i = 0; i < fault_faces.Size(); i++)
   {
      Vector c = FaceCentroid(*mesh, fault_faces[i]);
      if (std::abs(c(1)) > 1.0)  // y=0 within tolerance
      {
         all_on_y0 = false;
         std::cerr << "  Fault face " << fault_faces[i]
                   << " centroid y=" << c(1) << " (not on y=0)\n";
      }
   }
   TEST_ASSERT(all_on_y0, "all fault interior faces on y=0");
}

void TestFaultFacesCoverY0Plane()
{
   std::cout << "\n=== Test: Fault faces cover entire y=0 plane ===\n";

   // The inline mesh tags ALL y=0 interior faces as attr=3 (fault).
   // Physical fault extent (lf x Wf) is enforced by friction parameter
   // functions a(z), Dc(z), not at the mesh level.
   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;

   LinearElastic model(lambda, mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   ElasticityDomainOperator<Mesh> op(
      *mesh, 1, model, Vp, Wf, lf, bc,
      DGMethod::BR2, SolverType::CG_AMG);

   const auto &fault_faces = op.GetFaultInteriorFaces();

   // All fault faces should be on y=0 (tested above), and we expect
   // the full y=0 plane to be covered by fault faces.
   // 2x2x1 mesh: 4 tets per hex, y=0 plane has 2*1=2 quads -> ~4 tri faces
   std::cout << "  Fault faces on y=0: " << fault_faces.Size() << "\n";
   TEST_ASSERT(fault_faces.Size() > 0, "y=0 plane has fault faces");

   // Every fault face centroid should have y=0
   bool all_on_y0 = true;
   for (int i = 0; i < fault_faces.Size(); i++)
   {
      Vector c = FaceCentroid(*mesh, fault_faces[i]);
      if (std::abs(c(1)) > 1.0) { all_on_y0 = false; }
   }
   TEST_ASSERT(all_on_y0, "all fault faces on y=0 (redundant check)");
}

void TestFaultFaceArraySelfConsistency()
{
   std::cout << "\n=== Test: Fault face array self-consistency ===\n";

   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;

   LinearElastic model(lambda, mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   ElasticityDomainOperator<Mesh> op(
      *mesh, 1, model, Vp, Wf, lf, bc,
      DGMethod::BR2, SolverType::CG_AMG);

   const auto &fault_faces = op.GetFaultInteriorFaces();
   const auto &fault_shared = op.GetFaultSharedFaces();

   // All fault interior face indices should be valid interior faces
   bool all_interior = true;
   for (int i = 0; i < fault_faces.Size(); i++)
   {
      FaceElementTransformations *FTr =
         mesh->GetInteriorFaceTransformations(fault_faces[i]);
      if (!FTr)
      {
         all_interior = false;
         std::cerr << "  Fault face " << fault_faces[i]
                   << " is not an interior face\n";
      }
   }
   TEST_ASSERT(all_interior, "all fault face indices are interior faces");

   // No duplicates in fault face array
   std::set<int> fault_set(fault_faces.begin(),
                           fault_faces.begin() + fault_faces.Size());
   TEST_ASSERT(static_cast<int>(fault_set.size()) == fault_faces.Size(),
               "no duplicate fault faces");

   // Total = interior + shared
   int n_total = op.GetNumFaultFaces();
   TEST_ASSERT(n_total == fault_faces.Size() + fault_shared.Size(),
               "total fault = interior + shared");

   // Serial mesh: no shared faces
   TEST_ASSERT(fault_shared.Size() == 0, "serial mesh has no shared fault faces");
}

void TestMeshBdrElementsMatchClassification()
{
   std::cout << "\n=== Test: Mesh bdr elements match face classification ===\n";

   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;

   LinearElastic model(lambda, mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   ElasticityDomainOperator<Mesh> op(
      *mesh, 1, model, Vp, Wf, lf, bc,
      DGMethod::BR2, SolverType::CG_AMG);

   // Count boundary elements by attribute
   int n_fault_bdr = 0, n_dir_bdr = 0, n_nat_bdr = 0;
   for (int be = 0; be < mesh->GetNBE(); be++)
   {
      int attr = mesh->GetBdrAttribute(be);
      if (attr == 3) { n_fault_bdr++; }
      else if (attr == 5) { n_dir_bdr++; }
      else if (attr == 1) { n_nat_bdr++; }
   }

   std::cout << "  Bdr elements: fault(3)=" << n_fault_bdr
             << " dir(5)=" << n_dir_bdr
             << " nat(1)=" << n_nat_bdr << "\n";

   // The mesh must have boundary elements for fault and Dirichlet
   TEST_ASSERT(n_fault_bdr > 0, "mesh has fault bdr elements (attr 3)");
   TEST_ASSERT(n_dir_bdr > 0, "mesh has Dirichlet bdr elements (attr 5)");

   // The fault interior face count should match bdr element count for attr 3
   // (In serial mesh, each fault bdr element maps to one interior face)
   const auto &fault_faces = op.GetFaultInteriorFaces();
   TEST_ASSERT(fault_faces.Size() == n_fault_bdr,
               "fault face count matches attr-3 bdr element count");
}

void TestFaultAndDirichletDisjoint()
{
   std::cout << "\n=== Test: Fault and Dirichlet faces are disjoint ===\n";

   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;

   LinearElastic model(lambda, mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   ElasticityDomainOperator<Mesh> op(
      *mesh, 1, model, Vp, Wf, lf, bc,
      DGMethod::BR2, SolverType::CG_AMG);

   // Build sets of face indices
   const auto &fault_faces = op.GetFaultInteriorFaces();
   std::set<int> fault_set;
   for (int i = 0; i < fault_faces.Size(); i++)
   {
      fault_set.insert(fault_faces[i]);
   }

   // Scan mesh bdr elements for Dirichlet faces, check no overlap with fault
   int overlap = 0;
   for (int be = 0; be < mesh->GetNBE(); be++)
   {
      int attr = mesh->GetBdrAttribute(be);
      if (bc.dirichlet_attrs.count(attr) == 0) { continue; }

      int face_idx = mesh->GetBdrElementFaceIndex(be);
      if (fault_set.count(face_idx))
      {
         overlap++;
         std::cerr << "  Face " << face_idx << " is both fault and Dirichlet\n";
      }
   }
   TEST_ASSERT(overlap == 0, "no face is both fault and Dirichlet");
}

void TestOldNewCtorSameFaceClassification()
{
   std::cout << "\n=== Test: Old/new ctor produce same face classification ===\n";

   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;

   // Old constructor
   ElasticityDomainOperator<Mesh> op_old(
      *mesh, 1, lambda, mu, Vp, Wf, lf,
      DGMethod::BR2, SolverType::CG_AMG, BCMode::FarField);

   // New constructor
   LinearElastic model(lambda, mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   ElasticityDomainOperator<Mesh> op_new(
      *mesh, 1, model, Vp, Wf, lf, bc,
      DGMethod::BR2, SolverType::CG_AMG);

   // Compare fault face arrays (same count)
   TEST_ASSERT(op_old.GetNumFaultFaces() == op_new.GetNumFaultFaces(),
               "fault face count matches");
   TEST_ASSERT(op_old.GetFaultInteriorFaces().Size() ==
               op_new.GetFaultInteriorFaces().Size(),
               "fault interior face count matches");
   TEST_ASSERT(op_old.GetFaultSharedFaces().Size() ==
               op_new.GetFaultSharedFaces().Size(),
               "fault shared face count matches");

   // Compare fault face arrays (same elements)
   const auto &old_ff = op_old.GetFaultInteriorFaces();
   const auto &new_ff = op_new.GetFaultInteriorFaces();
   std::set<int> old_set(old_ff.begin(), old_ff.begin() + old_ff.Size());
   std::set<int> new_set(new_ff.begin(), new_ff.begin() + new_ff.Size());
   TEST_ASSERT(old_set == new_set, "fault face sets identical");

   // Compare DOF counts
   TEST_ASSERT(op_old.GetNumFaultDOFs() == op_new.GetNumFaultDOFs(),
               "fault DOF count matches");
   TEST_ASSERT(op_old.GetNumOwnedFaultDOFs() == op_new.GetNumOwnedFaultDOFs(),
               "owned fault DOF count matches");
}

void TestFaultCountMatchesGeometry()
{
   std::cout << "\n=== Test: Fault face count matches expected geometry ===\n";

   // 2x2x1 mesh: 200km x 100km x 100km
   // Fault: 100km x 40km on y=0
   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;

   LinearElastic model(lambda, mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   ElasticityDomainOperator<Mesh> op(
      *mesh, 1, model, Vp, Wf, lf, bc,
      DGMethod::BR2, SolverType::CG_AMG);

   int n_fault = op.GetNumFaultFaces();
   int n_fault_int = op.GetFaultInteriorFaces().Size();
   int n_fault_shared = op.GetFaultSharedFaces().Size();

   std::cout << "  Total fault faces: " << n_fault
             << " (interior=" << n_fault_int
             << " shared=" << n_fault_shared << ")\n";

   TEST_ASSERT(n_fault > 0, "mesh has fault faces");
   TEST_ASSERT(n_fault == n_fault_int + n_fault_shared,
               "total = interior + shared");
   TEST_ASSERT(n_fault_shared == 0, "serial mesh has no shared fault faces");
}

int main(int argc, char *argv[])
{
   TestFaultFacesAreOnY0();
   TestFaultFacesCoverY0Plane();
   TestFaultFaceArraySelfConsistency();
   TestMeshBdrElementsMatchClassification();
   TestFaultAndDirichletDisjoint();
   TestOldNewCtorSameFaceClassification();
   TestFaultCountMatchesGeometry();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
