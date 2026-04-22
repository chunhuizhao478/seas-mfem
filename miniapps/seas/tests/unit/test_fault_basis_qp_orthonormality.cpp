// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan §4.3 (regression gate): per-QP fault-basis
// orthonormality on a TPV102 fault mesh.
//
// For every quadrature point on every fault face, verifies:
//   1. |n|=1, |t1|=1, |t2|=1          to <= 16 ULP
//   2. n · t1 = 0, n · t2 = 0,
//      t1 · t2 = 0                    to <= 16 ULP
//
// Default mesh: `tpv102/mesh/tpv102_1000m.msh` (10 MB).  The plan's
// rev-2 §4.3 text says "production 200 m mesh", but feedback memory
// "TPV102 no local reproducer runs on production mesh" argues against
// loading the 126 MB file for a regression gate.  A failure on 1000 m
// would reflect a FaultBasis bug that is mesh-independent.  Pass
// `--mesh <path>` to override if the 200 m mesh is required.
//
// Rationale (plan §3 H-V92-F, §A): REVIEW R-004 identified the 16-ULP
// budget (9-term FMA round-off).  Axis-aligned BP5 faces are
// bit-exact; Gmsh-generated tet meshes may show up to 5 ULP drift on
// the |dip| norm (pre-v9.1.0 R-003) — this test now guards against
// regression of the R-003 dip-normalization patch.

#include "mfem.hpp"
#include "../../fault/fault_basis.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static constexpr double kEps       = 2.220446049250313e-16;
static constexpr int    kUlpBudget = 16;
static constexpr double kTolUlp    = kUlpBudget * kEps;

static inline double dot3(const real_t a[3], const real_t b[3])
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline double norm3(const real_t v[3])
{
   return std::sqrt(dot3(v, v));
}

// Check one triple of basis vectors; return worst deviation from the
// orthonormality ideal.  Reports both |v|=1 drift and off-diagonal dot.
struct FrameStats
{
   double worst_norm_dev = 0.0;   // max(||·| − 1|)
   double worst_dot_dev  = 0.0;   // max(|n·t1|, |n·t2|, |t1·t2|)
};

static FrameStats CheckFrame(const real_t n[3], const real_t t1[3],
                              const real_t t2[3])
{
   FrameStats s;
   s.worst_norm_dev = std::max({std::abs(norm3(n)  - 1.0),
                                std::abs(norm3(t1) - 1.0),
                                std::abs(norm3(t2) - 1.0)});
   s.worst_dot_dev  = std::max({std::abs(dot3(n, t1)),
                                std::abs(dot3(n, t2)),
                                std::abs(dot3(t1, t2))});
   return s;
}

static std::string GetArg(int argc, char **argv, const std::string &flag,
                          const std::string &default_value)
{
   for (int i = 1; i < argc - 1; i++)
   {
      if (flag == argv[i]) { return argv[i+1]; }
   }
   return default_value;
}

static int GetIntArg(int argc, char **argv, const std::string &flag, int def)
{
   for (int i = 1; i < argc - 1; i++)
   {
      if (flag == argv[i]) { return std::atoi(argv[i+1]); }
   }
   return def;
}

int main(int argc, char **argv)
{
   std::cout << "\n=== TPV102 v9.2.0 §4.3 Regression: "
             << "Per-QP Fault-Basis Orthonormality ===\n";

   const std::string mesh_file = GetArg(argc, argv, "--mesh",
      "tpv102/mesh/tpv102_1000m.msh");
   const int fault_attr = GetIntArg(argc, argv, "--fault-attr", 3);
   const int order      = GetIntArg(argc, argv, "--order", 1);

   std::cout << "  mesh        : " << mesh_file << "\n";
   std::cout << "  fault_attr  : " << fault_attr << "\n";
   std::cout << "  order       : " << order
             << "  (face ir = 2*order = " << 2*order << ")\n";

   // -------- Load mesh ------------------------------------------------------
   Mesh mesh(mesh_file.c_str(), 1, 1);
   MFEM_VERIFY(mesh.Dimension() == 3, "Test needs a 3D mesh");

   std::cout << "  mesh dim    : " << mesh.Dimension() << "\n";
   std::cout << "  elements    : " << mesh.GetNE() << "\n";
   std::cout << "  faces       : " << mesh.GetNumFaces() << "\n";
   std::cout << "  bdr elems   : " << mesh.GetNBE() << "\n";

   // -------- Identify fault interior faces (same pattern as wave_op) -------
   Array<int> fault_faces;
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      if (mesh.GetBdrAttribute(b) != fault_attr) { continue; }
      int face_idx = mesh.GetBdrElementFaceIndex(b);
      if (mesh.GetInteriorFaceTransformations(face_idx) == nullptr) { continue; }
      fault_faces.Append(face_idx);
   }
   std::cout << "  fault faces : " << fault_faces.Size() << "\n";

   if (fault_faces.Size() == 0)
   {
      std::cout << "  SKIPPED (no fault faces found with attribute "
                << fault_attr << ").  Ensure mesh has BooleanFragments "
                << "+ fault boundary tag, or override with --fault-attr.\n";
      return 77;  // "skipped" sentinel (automake convention)
   }

   // -------- Compute centroid + QP bases -----------------------------------
   Vector ref_normal(3); ref_normal = 0.0; ref_normal(1) = -1.0;
   Vector up(3);         up = 0.0;         up(2) = 1.0;

   FaultBasis fb;
   fb.Compute(mesh, fault_faces, ref_normal, up);

   // Determine face geometry from the first fault face.  ComputeQPBasis
   // needs an IntegrationRule matching the face geometry.
   Geometry::Type face_geom = Geometry::TRIANGLE;
   {
      auto *ftr0 = mesh.GetInteriorFaceTransformations(fault_faces[0]);
      if (ftr0) { face_geom = ftr0->GetGeometryType(); }
   }
   const IntegrationRule &face_ir = IntRules.Get(face_geom, 2*order);
   fb.ComputeQPBasis(mesh, fault_faces, ref_normal, up, face_ir);

   const int nq = face_ir.GetNPoints();
   std::cout << "  QPs/face    : " << nq << "\n";

   // -------- Walk every (face, QP) pair and aggregate stats ---------------
   num_tests++;
   double worst_norm = 0.0;
   double worst_dot  = 0.0;
   int    worst_face = -1, worst_qp = -1;
   std::string worst_kind;
   long long total_qp = 0;

   for (int fi = 0; fi < fb.NumFaces(); fi++)
   {
      const FaultBasisData &bd = fb.GetBasis(fi);
      MFEM_VERIFY((int)bd.qp_data.size() == nq,
                  "FaultBasis qp_data size mismatch on face " << fi);

      for (int q = 0; q < nq; q++)
      {
         const FaultBasisQPData &qpd = bd.qp_data[q];
         FrameStats s = CheckFrame(qpd.normal, qpd.tangent1, qpd.tangent2);
         total_qp++;
         if (s.worst_norm_dev > worst_norm)
         {
            worst_norm = s.worst_norm_dev;
            worst_face = fi; worst_qp = q; worst_kind = "|v|=1";
         }
         if (s.worst_dot_dev  > worst_dot)
         {
            worst_dot  = s.worst_dot_dev;
            if (s.worst_dot_dev > worst_norm)
            {
               worst_face = fi; worst_qp = q; worst_kind = "orth";
            }
         }
      }
   }

   std::cout << "  Total QPs checked : " << total_qp << "\n";
   std::cout << "  worst_norm_dev    : " << std::scientific
             << std::setprecision(3) << worst_norm
             << " (" << worst_norm / kEps << " ulp) at face "
             << worst_face << " qp " << worst_qp << "\n";
   std::cout << "  worst_dot_dev     : " << std::scientific
             << std::setprecision(3) << worst_dot
             << " (" << worst_dot / kEps << " ulp)\n";

   const double worst_overall = std::max(worst_norm, worst_dot);
   if (worst_overall <= kTolUlp)
   {
      num_passed++;
      std::cout << "  PASSED: per-QP frame orthonormality <= 16 ULP "
                << "on all " << total_qp << " QPs of "
                << fb.NumFaces() << " fault faces\n";
   }
   else
   {
      num_failed++;
      std::cout << "  FAILED: per-QP frame orthonormality exceeded 16 ULP "
                << "(worst " << worst_overall / kEps << " ulp, kind "
                << worst_kind << ")\n";
   }

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";

   return (num_failed == 0) ? 0 : 1;
}
