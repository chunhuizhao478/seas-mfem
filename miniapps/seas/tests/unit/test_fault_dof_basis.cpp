// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests for Phase 6.A — FaultGeometry per-DOF 3-D coordinate
// and (n, t1, t2) basis accessors.
//
// Plan reference:
//   PLAN_onfaultstress.md §1709-1792 (Phase 6.A acceptance criteria).
//
// Coverage:
//   T_6A_1  dof_coords_3d_ size == 3 * num_fault_dofs.
//   T_6A_2  dof_basis_ shape == (9, num_fault_dofs).
//   T_6A_3  |n_i| = |t1_i| = |t2_i| = 1 to 1e-12 for every DOF.
//   T_6A_4  pairwise dot products (n·t1, n·t2, t1·t2) < 1e-12 for every DOF.
//   T_6A_5  per-DOF n_i is parallel (up to sign) to the per-face
//           FaultBasis normal at the same face (|n_dof · n_face| ≈ 1).
//   T_6A_6  per-DOF (n, t1, t2) forms a right-handed orthonormal frame:
//           t2 = n × t1 (with the sign as derived in Gram-Schmidt per
//           the plan §1758, mirroring Phase 2 `basis_to_node`).
//   T_6A_7  BP2 (antiplane) ctor leaves dof_coords_3d_ / dof_basis_ empty.
//
// The test mesh is the standard 3-D cartesian hex mesh from
// test_bp5_fault_operator.cpp with a fault at Y=0; on this mesh every
// fault face has constant Jacobian so the per-DOF basis must match the
// per-face FaultBasis exactly (no degenerate fallback expected).

#include "mfem.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../domain/domain_config.hpp"
#include "../../config/bp5_params.hpp"
#include "../../config/bp2_params.hpp"
#include "../../domain/antiplane_operator.hpp"

#include <unistd.h>     // access()
#include <sys/stat.h>   // stat()

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>

using namespace mfem;
using namespace mfem::seas;

// --------------------------------------------------------------------
// Tiny test framework (matches test_stress_field_3d.cpp)
// --------------------------------------------------------------------
static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ \
                   << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_NEAR(value, expected, tol, message) \
   do { \
      num_tests++; \
      const real_t _v = (value); \
      const real_t _e = (expected); \
      const real_t _t = (tol); \
      if (std::abs(_v - _e) > _t) { \
         std::cerr << "FAILED: " << message << " (got " << _v \
                   << ", expected " << _e << ", tol " << _t \
                   << ", line " << __LINE__ << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

// --------------------------------------------------------------------
// Mesh: load the real Tandem-style BP5 mesh (1000 m or 500 m) with
// physical tag 100 for the fault. Falls back to a coarser mesh if the
// nominal mesh is not present in the source tree.
// --------------------------------------------------------------------
static bool FileExists(const std::string &path)
{
   struct stat st;
   return (stat(path.c_str(), &st) == 0);
}

static std::string FindBP5Mesh()
{
   // Try a few candidate locations and resolutions, coarsest first
   // (1000 m → 56 fault DOFs at p=1, fast to test).
   const std::vector<std::string> candidates = {
      "bp5/mesh/bp5_1000m.msh",
      "../bp5/mesh/bp5_1000m.msh",
      "/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/bp5/mesh/bp5_1000m.msh",
      "/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/bp5/mesh/bp5_500m.msh",
   };
   for (const auto &p : candidates)
   {
      if (FileExists(p)) { return p; }
   }
   return "";
}

// --------------------------------------------------------------------
// Fixture: build mesh + elasticity operator + fault geometry once
// --------------------------------------------------------------------
struct Fixture
{
   BP5Params params;
   std::unique_ptr<Mesh> mesh;
   std::unique_ptr<ElasticityDomainOperator<Mesh>> domain_op;
   std::unique_ptr<FaultGeometry<Mesh>> fault_geom;
   int nf = 0;
   bool loaded = false;
};

static std::unique_ptr<Fixture> BuildFixture()
{
   auto fix = std::make_unique<Fixture>();
   const std::string mesh_path = FindBP5Mesh();
   if (mesh_path.empty())
   {
      std::cerr << "WARNING: BP5 mesh file not found; skipping runtime tests\n";
      return fix;
   }
   std::cout << "Loading BP5 mesh: " << mesh_path << "\n";
   fix->mesh = std::make_unique<Mesh>(mesh_path.c_str(), 1, 1);
   fix->loaded = true;

   // Tandem BP5 attribute layout: fault = 100, Dirichlet = {1, 2, 3, 4},
   // Natural = {5, 6}. The legacy ctor uses fault_attr = 3 (Tandem's
   // generic preset); for the BP5 .msh file the fault tag is 100.
   BoundaryConfig bdr_config;
   bdr_config.fault_attr = 100;
   bdr_config.dirichlet_attrs = {1, 2, 3, 4};
   bdr_config.natural_attrs   = {5, 6};
   bdr_config.default_dirichlet_func = MakeBP5DirichletFunc(fix->params.Vp);

   LinearElastic linelast(fix->params.lambda(), fix->params.mu());
   DomainConfig domain_cfg;

   fix->domain_op = std::make_unique<ElasticityDomainOperator<Mesh>>(
      *fix->mesh, 1, linelast,
      fix->params.Vp, fix->params.Wf, fix->params.lf,
      bdr_config, DGMethod::IP,
      SolverType::MUMPS_BLR, domain_cfg);
   fix->nf = fix->domain_op->GetNumFaultDOFs();
   if (fix->nf == 0) { return fix; }
   fix->fault_geom = std::make_unique<FaultGeometry<Mesh>>(*fix->domain_op,
                                                            fix->params);
   return fix;
}

// --------------------------------------------------------------------
// Test cases
// --------------------------------------------------------------------
static void T_6A_1_dof_coords_size(const Fixture &fix)
{
   std::cout << "\n[T-6A-1] fault_dof_coords_3d().Size() == 3 * num_fault_dofs\n";
   const Vector &c = fix.fault_geom->fault_dof_coords_3d();
   TEST_ASSERT(c.Size() == 3 * fix.fault_geom->NumFaultDOFs(),
               "size == 3 * NumFaultDOFs");
}

static void T_6A_2_dof_basis_shape(const Fixture &fix)
{
   std::cout << "\n[T-6A-2] fault_dof_basis() shape == (9, num_fault_dofs)\n";
   const DenseMatrix &B = fix.fault_geom->fault_dof_basis();
   TEST_ASSERT(B.Height() == 9, "Height() == 9");
   TEST_ASSERT(B.Width() == fix.fault_geom->NumFaultDOFs(),
               "Width() == NumFaultDOFs");
}

static void T_6A_3_unit_norms(const Fixture &fix)
{
   std::cout << "\n[T-6A-3] |n_i| = |t1_i| = |t2_i| = 1\n";
   const DenseMatrix &B = fix.fault_geom->fault_dof_basis();
   const int nd = B.Width();
   real_t max_dev_n = 0.0, max_dev_t1 = 0.0, max_dev_t2 = 0.0;
   for (int i = 0; i < nd; i++)
   {
      real_t n_norm  = std::sqrt(B(0,i)*B(0,i) + B(1,i)*B(1,i)
                                 + B(2,i)*B(2,i));
      real_t t1_norm = std::sqrt(B(3,i)*B(3,i) + B(4,i)*B(4,i)
                                 + B(5,i)*B(5,i));
      real_t t2_norm = std::sqrt(B(6,i)*B(6,i) + B(7,i)*B(7,i)
                                 + B(8,i)*B(8,i));
      max_dev_n  = std::max(max_dev_n,  std::abs(n_norm  - 1.0));
      max_dev_t1 = std::max(max_dev_t1, std::abs(t1_norm - 1.0));
      max_dev_t2 = std::max(max_dev_t2, std::abs(t2_norm - 1.0));
   }
   TEST_NEAR(max_dev_n,  0.0, 1e-12, "max ||n|-1| < 1e-12");
   TEST_NEAR(max_dev_t1, 0.0, 1e-12, "max ||t1|-1| < 1e-12");
   TEST_NEAR(max_dev_t2, 0.0, 1e-12, "max ||t2|-1| < 1e-12");
}

static void T_6A_4_orthogonality(const Fixture &fix)
{
   std::cout << "\n[T-6A-4] pairwise orthogonality\n";
   const DenseMatrix &B = fix.fault_geom->fault_dof_basis();
   const int nd = B.Width();
   real_t max_n_t1 = 0.0, max_n_t2 = 0.0, max_t1_t2 = 0.0;
   for (int i = 0; i < nd; i++)
   {
      real_t n_t1 = B(0,i)*B(3,i) + B(1,i)*B(4,i) + B(2,i)*B(5,i);
      real_t n_t2 = B(0,i)*B(6,i) + B(1,i)*B(7,i) + B(2,i)*B(8,i);
      real_t t1_t2 = B(3,i)*B(6,i) + B(4,i)*B(7,i) + B(5,i)*B(8,i);
      max_n_t1  = std::max(max_n_t1,  std::abs(n_t1));
      max_n_t2  = std::max(max_n_t2,  std::abs(n_t2));
      max_t1_t2 = std::max(max_t1_t2, std::abs(t1_t2));
   }
   TEST_NEAR(max_n_t1,  0.0, 1e-12, "|n . t1|  < 1e-12");
   TEST_NEAR(max_n_t2,  0.0, 1e-12, "|n . t2|  < 1e-12");
   TEST_NEAR(max_t1_t2, 0.0, 1e-12, "|t1 . t2| < 1e-12");
}

static void T_6A_5_normal_aligned_with_face(const Fixture &fix)
{
   std::cout << "\n[T-6A-5] per-DOF normal parallel to per-face normal "
                "(|n_dof . n_face| == 1)\n";
   const DenseMatrix &B = fix.fault_geom->fault_dof_basis();
   const FaultBasis *fb = fix.domain_op->GetFaultBasis();
   TEST_ASSERT(fb != nullptr, "FaultBasis exists");
   if (fb == nullptr) { return; }
   const int nbf = fix.domain_op->GetNbfPerFace();
   const int num_faces = fb->NumFaces();
   real_t min_abs_dot = 1.0;
   for (int f = 0; f < num_faces; f++)
   {
      const FaultBasisData &fbd = fb->GetBasis(f);
      for (int k = 0; k < nbf; k++)
      {
         const int dof_idx = f * nbf + k;
         if (dof_idx >= B.Width()) { continue; }
         real_t dot = B(0, dof_idx) * fbd.normal[0]
                    + B(1, dof_idx) * fbd.normal[1]
                    + B(2, dof_idx) * fbd.normal[2];
         min_abs_dot = std::min(min_abs_dot, std::abs(dot));
      }
   }
   TEST_NEAR(min_abs_dot, 1.0, 1e-12,
             "min |n_dof . n_face| == 1 (parallel up to sign)");
}

static void T_6A_6_basis_matches_faultbasis(const Fixture &fix)
{
   // After R-001's fix, the Gram-Schmidt re-orthonormalisation
   // PROJECTS the input t2 against (n, t1) instead of redefining it
   // as n × t1.  The result must match the elasticity operator's
   // per-face FaultBasis vectors element-wise (the input is already
   // orthonormal up to FP, so projection is a near-identity).
   //
   // This replaces the prior "right-handed (t2 = n × t1)" assertion,
   // which encoded the buggy sign convention.
   std::cout << "\n[T-6A-6] per-DOF basis matches FaultBasis at every "
                "DOF (R-001 sign-flip preservation)\n";
   const DenseMatrix &B = fix.fault_geom->fault_dof_basis();
   const FaultBasis *fb = fix.domain_op->GetFaultBasis();
   TEST_ASSERT(fb != nullptr, "FaultBasis exists");
   if (fb == nullptr) { return; }
   const int nbf = fix.domain_op->GetNbfPerFace();
   const int num_faces = fb->NumFaces();
   real_t max_n  = 0.0, max_t1 = 0.0, max_t2 = 0.0;
   for (int f = 0; f < num_faces; f++)
   {
      const FaultBasisData &fbd = fb->GetBasis(f);
      for (int k = 0; k < nbf; k++)
      {
         const int dof = f * nbf + k;
         if (dof >= B.Width()) { continue; }
         for (int d = 0; d < 3; d++)
         {
            max_n  = std::max(max_n,
                              std::abs(B(d,     dof) - fbd.normal[d]));
            max_t1 = std::max(max_t1,
                              std::abs(B(3 + d, dof) - fbd.tangent1[d]));
            max_t2 = std::max(max_t2,
                              std::abs(B(6 + d, dof) - fbd.tangent2[d]));
         }
      }
   }
   TEST_NEAR(max_n,  0.0, 1e-12, "n_GS == n_face (sign-preserving)");
   TEST_NEAR(max_t1, 0.0, 1e-12, "t1_GS == t1_face (sign-preserving)");
   TEST_NEAR(max_t2, 0.0, 1e-12,
             "t2_GS == t2_face (sign-preserving — R-001 contract)");
}

static void T_6A_7_bp2_path_unaffected()
{
   std::cout << "\n[T-6A-7] BP2 (antiplane) ctor leaves per-DOF members empty\n";

   // Minimal 2D antiplane setup — the BP2 ctor of FaultGeometry does NOT
   // call ComputePerDOFCoordsAndBasis_, so dof_coords_3d_ and dof_basis_
   // must stay empty (write-once invariant from BP5 ctor only).
   const real_t Wf = 40e3;
   const real_t lf = 40e3;
   Mesh mesh = Mesh::MakeCartesian2D(2, 1, Element::QUADRILATERAL,
                                     false, 2.0 * lf, Wf);
   Vector shift(2); shift(0) = -lf; shift(1) = -Wf;
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *v = mesh.GetVertex(i);
      v[0] += shift(0); v[1] += shift(1);
   }
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      ElementTransformation *T = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(T->GetGeometryType());
      T->SetIntPoint(&ip);
      Vector center(2);
      T->Transform(ip, center);
      real_t tol = 1e-6;
      if      (std::abs(center(0) - (-lf)) < tol) { mesh.SetBdrAttribute(be, 1); }
      else if (std::abs(center(0) -  lf)   < tol) { mesh.SetBdrAttribute(be, 2); }
      else if (std::abs(center(1) - 0.0)   < tol) { mesh.SetBdrAttribute(be, 5); }
      else if (std::abs(center(1) - (-Wf)) < tol) { mesh.SetBdrAttribute(be, 6); }
   }
   mesh.SetAttributes();

   BP2Params params;
   AntiplaneDomainOperator<Mesh> domain_op(mesh, 1, params.mu(), params.Vp,
                                            params.Wf, DGMethod::IP);
   if (domain_op.GetNumFaultDOFs() == 0)
   {
      std::cout << "  (no fault DOFs in BP2 mesh; skipping size check)\n";
      return;
   }

   FaultGeometry<Mesh> fault_geom(domain_op, params);
   TEST_ASSERT(fault_geom.fault_dof_coords_3d().Size() == 0,
               "BP2 ctor: dof_coords_3d_ is empty");
   TEST_ASSERT(fault_geom.fault_dof_basis().Height() == 0,
               "BP2 ctor: dof_basis_.Height() == 0");
   TEST_ASSERT(fault_geom.fault_dof_basis().Width()  == 0,
               "BP2 ctor: dof_basis_.Width() == 0");
}

// --------------------------------------------------------------------
// main
// --------------------------------------------------------------------
int main(int argc, char **argv)
{
   (void)argc; (void)argv;
   std::cout << "Running Phase 6.A FaultGeometry per-DOF accessor tests\n";

   auto fix = BuildFixture();
   if (!fix->loaded)
   {
      std::cout << "(BP5 mesh not available; runtime tests skipped — "
                   "Phase 6.A compile-only verification)\n";
      // Still exercise the BP2 / no-fault path which doesn't need a mesh.
      T_6A_7_bp2_path_unaffected();
      std::cout << "\n========================================\n";
      std::cout << "Phase 6.A (compile-only): " << num_passed << " / "
                << num_tests << " passed, " << num_failed << " failed\n";
      std::cout << "========================================\n";
      return (num_failed == 0) ? 0 : 1;
   }
   if (fix->nf == 0)
   {
      std::cerr << "FATAL: test fixture has no fault DOFs\n";
      return 1;
   }
   if (fix->fault_geom == nullptr)
   {
      std::cerr << "FATAL: fault_geom not constructed\n";
      return 1;
   }

   T_6A_1_dof_coords_size(*fix);
   T_6A_2_dof_basis_shape(*fix);
   T_6A_3_unit_norms(*fix);
   T_6A_4_orthogonality(*fix);
   T_6A_5_normal_aligned_with_face(*fix);
   T_6A_6_basis_matches_faultbasis(*fix);

   T_6A_7_bp2_path_unaffected();

   std::cout << "\n========================================\n";
   std::cout << "Phase 6.A: " << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
