// Isolated equivalence test for the non-fault interior face geometry cache
// (Phase 1 of action_plan_faceflux_cache_2026-06-24.md).
//
//   1. On a fault-free tet box: the cached {nor (signed), w, shape1, shape2}
//      reproduce a fresh on-the-fly GetFaceElementTransformations + CalcOrtho +
//      CalcShape to <= 1e-12, and the cached face set equals the runtime
//      predicate's interior-non-fault set (boundary faces excluded).
//   2. On a 2-tet FAULT mesh: the fault face is ABSENT from the cache
//      (fault-exclusion — REVIEW R-002).
//
// No runtime code path is touched; this exercises only the builder.

#include "mfem.hpp"
#include "../../dynamic/face_geom_cache.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <unordered_map>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_failed = 0;
#define CHECK(cond, msg) do {                                               \
   if (cond) { std::cout << "  PASSED: " << msg << "\n"; }                  \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: "      \
                                  << msg << "\n"; } } while (0)

// face -> boundary attr, identical to WaveOperator (wave_operator.inl:119-124).
static std::vector<int> BuildFaceBdrAttr(Mesh &mesh)
{
   std::vector<int> fba(mesh.GetNumFaces(), 0);
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      const int fi = mesh.GetBdrElementFaceIndex(b);
      fba[fi] = mesh.GetBdrAttribute(b);
   }
   return fba;
}

// 2-tet mesh sharing one interior face at y=0 tagged FAULT(3); external FREE(1).
// Copied verbatim from test_bimaterial_wave_operator_parity::BuildTwoTetFaultMesh.
static Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
      {0.0, 1.0, 0.0}, {0.0, -1.0, 0.0},
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);
   mesh.AddTet(0, 1, 2, 3, 1);
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-10)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);   // FAULT
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);          // FREE
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

int main()
{
   std::cout << "=== FaceGeomCache builder equivalence + exclusion ===\n";

   // ---- Test 1: fault-free box, geometry equivalence + set parity ----------
   for (int order : {1, 2})
   {
      Mesh mesh = Mesh::MakeCartesian3D(4, 4, 4, Element::TETRAHEDRON,
                                        1.0, 1.0, 1.0);
      for (int b = 0; b < mesh.GetNBE(); b++) { mesh.SetBdrAttribute(b, 5); }
      mesh.SetAttributes();

      L2_FECollection fec(order, mesh.Dimension());
      FiniteElementSpace fes(&mesh, &fec);            // vdim=1 (geometry-only)

      std::vector<int> fba = BuildFaceBdrAttr(mesh);
      std::set<int> no_shared;
      const int quad_order = 2 * order;

      std::unordered_map<int, FaceGeomEntry> cache;
      BuildNonFaultInteriorFaceGeomCache(mesh, fes, fba, /*fault_attr=*/0,
                                         no_shared, quad_order, cache);

      // Independent recomputation of the runtime predicate + geometry.
      std::set<int> expected_faces;
      double max_geom_err = 0.0;
      int    n_boundary_in_cache = 0;
      for (int f = 0; f < mesh.GetNumFaces(); f++)
      {
         FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
         if (!ftr) { continue; }
         const int e1 = ftr->Elem1No, e2 = ftr->Elem2No;
         if (e2 < 0) { continue; }              // boundary -> not cached
         expected_faces.insert(f);

         auto it = cache.find(f);
         if (it == cache.end()) { continue; }   // set mismatch caught below
         const FaceGeomEntry &fc = it->second;
         if (e2 < 0) { n_boundary_in_cache++; }

         const FiniteElement *fe1 = fes.GetFE(e1);
         const FiniteElement *fe2 = fes.GetFE(e2);
         const IntegrationRule &ir =
            IntRules.Get(ftr->GetGeometryType(), quad_order);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            Vector nor_vec(3);
            CalcOrtho(ftr->Face->Jacobian(), nor_vec);
            const double nor_len = nor_vec.Norml2();
            if (nor_len > 0) { nor_vec /= nor_len; }
            const double w = ip.weight * nor_len;

            // signed, component-wise normal comparison (REVIEW R-004).
            for (int d = 0; d < 3; d++)
            {
               max_geom_err = std::max(max_geom_err,
                                       std::abs(fc.nor[q][d] - nor_vec(d)));
            }
            max_geom_err = std::max(max_geom_err, std::abs(fc.w[q] - w));

            Vector s1(fe1->GetDof()); IntegrationPoint ip1;
            ftr->Loc1.Transform(ip, ip1); fe1->CalcShape(ip1, s1);
            Vector s2(fe2->GetDof()); IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2); fe2->CalcShape(ip2, s2);
            for (int i = 0; i < fe1->GetDof(); i++)
            {
               max_geom_err = std::max(max_geom_err,
                                       std::abs(fc.shape1[q](i) - s1(i)));
               max_geom_err = std::max(max_geom_err,
                                       std::abs(fc.shape2[q](i) - s2(i)));
            }
         }
      }

      std::set<int> cached_faces;
      for (const auto &kv : cache) { cached_faces.insert(kv.first); }

      CHECK(max_geom_err <= 1e-12,
            "order=" + std::to_string(order) +
            ": cached geometry == on-the-fly (signed), max err " +
            std::to_string(max_geom_err));
      CHECK(cached_faces == expected_faces,
            "order=" + std::to_string(order) +
            ": cached set == interior-non-fault predicate set (" +
            std::to_string(cached_faces.size()) + " faces)");
      CHECK(n_boundary_in_cache == 0,
            "order=" + std::to_string(order) + ": no boundary face cached");
   }

   // ---- Test 2: fault mesh, fault-exclusion (REVIEW R-002) -----------------
   {
      Mesh mesh = BuildTwoTetFaultMesh();
      L2_FECollection fec(1, mesh.Dimension());
      FiniteElementSpace fes(&mesh, &fec);
      std::vector<int> fba = BuildFaceBdrAttr(mesh);
      std::set<int> no_shared;

      // locate the fault face (the interior face with bdr_attr == 3).
      int fault_face = -1;
      for (int f = 0; f < mesh.GetNumFaces(); f++)
      {
         FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
         if (ftr && ftr->Elem2No >= 0 && fba[f] == 3) { fault_face = f; }
      }
      CHECK(fault_face >= 0, "fault face located (interior, bdr_attr==3)");

      std::unordered_map<int, FaceGeomEntry> cache;
      BuildNonFaultInteriorFaceGeomCache(mesh, fes, fba, /*fault_attr=*/3,
                                         no_shared, /*quad_order=*/2, cache);

      CHECK(cache.find(fault_face) == cache.end(),
            "fault face absent from cache (fault-exclusion)");
   }

   std::cout << (num_failed == 0 ? "=== PASSED ===\n" : "=== FAILED ===\n");
   return num_failed == 0 ? 0 : 1;
}
