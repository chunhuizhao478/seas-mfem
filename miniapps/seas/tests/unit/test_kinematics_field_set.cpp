// PLAN_split_bulk_solutions_2026-05-12 — acceptance test for the field
// set registered with the primary (kinematics) volume PV collection.
//
// Guards two REVIEW.md findings from the post-split-bulk audit:
//
//   R-001 (MODERATE) — `InitFaultOutputBP5` used to register 12 L2-p0
//       fault projection fields with `pv_dc_` unconditionally.  After
//       the fix, the registrations are gated behind
//       `SetRegisterFaultProjectionsInVolumePV(true)`.  Default OFF
//       means only the driver-side explicit `RegisterDomainField` calls
//       (e.g. `velocity` + `mpi_rank` for TPV*, `displacement` +
//       `mpi_rank` for BP5) land in kinematics.vtkhdf.
//
// What this test verifies:
//   Scenario A (default toggle OFF):
//     - Register two probe fields ("velocity", "mpi_rank") with a
//       freshly-constructed `seas::ParaViewOutput<Mesh>`.
//     - Call `InitFaultOutputBP5(...)`.
//     - Assert the underlying DataCollection has the two probe fields
//       AND NO L2-p0 fault projection (the 12 names from
//       paraview_output.hpp:784-801).
//
//   Scenario B (toggle ON):
//     - Same setup, but call `SetRegisterFaultProjectionsInVolumePV(true)`
//       BEFORE `InitFaultOutputBP5`.
//     - Assert ALL 12 L2-p0 projection names are present.
//
// The test runs on a synthetic 2x1x1 hex serial mesh with a single
// fault face — the smallest possible setup that satisfies the
// InitFaultOutputBP5 contract (one interior fault face, valid
// FaceElementTransformations on both sides).

#include "mfem.hpp"
#include "../../io/paraview_output.hpp"

#include <array>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_TRUE(cond, msg) do {                          \
   num_tests++;                                            \
   if (cond) {                                             \
      num_passed++;                                        \
      std::cout << "  PASSED: " << msg << "\n";            \
   } else {                                                \
      num_failed++;                                        \
      std::cout << "  FAILED [" << __LINE__ << "]: "       \
                << msg << "\n";                            \
   }                                                       \
} while (0)

// -------------------------------------------------------------------
// Build a minimal serial hex mesh with one interior face that we'll
// declare as "fault" for the purpose of InitFaultOutputBP5.  Layout:
// two unit hexes sharing the face x=1 (face index 2 by Mesh::AddHex
// convention; we discover it programmatically).
// -------------------------------------------------------------------
static Mesh BuildTwoHexMesh()
{
   const int dim = 3;
   const int nverts = 12;
   const int nhex   = 2;
   Mesh mesh(dim, nverts, nhex, /*NBdrElem=*/0, /*spaceDim=*/3);
   // 12 vertices of two stacked unit cubes along x.
   const double verts[12][3] = {
      {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
      {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1},
      {2,0,0}, {2,1,0}, {2,0,1}, {2,1,1},
   };
   for (int i = 0; i < nverts; ++i)
   {
      mesh.AddVertex(verts[i][0], verts[i][1], verts[i][2]);
   }
   const int hex_a[8] = {0, 1, 2, 3, 4, 5, 6, 7};
   const int hex_b[8] = {1, 8, 9, 2, 5, 10, 11, 6};
   mesh.AddHex(hex_a, /*attr=*/1);
   mesh.AddHex(hex_b, /*attr=*/1);
   mesh.FinalizeHexMesh(/*generate_edges=*/0, /*refine=*/0,
                        /*fix_orientation=*/true);
   return mesh;
}

// Locate the (one) interior face shared by both hexes.
static int FindInteriorFace(Mesh &mesh)
{
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      FaceElementTransformations *FTr =
         mesh.GetInteriorFaceTransformations(f);
      if (FTr) { return f; }
   }
   return -1;
}

// Names of the 12 L2-p0 fault projection fields that
// `InitFaultOutputBP5` registers when the opt-in toggle is set.  Must
// stay in sync with paraview_output.hpp:784-801.
static const std::vector<std::string> kFaultProjectionNames = {
   "slip_dip", "slip_strike",
   "slip_rate_dip", "slip_rate_strike",
   "traction_dip", "traction_strike",
   "state_variable", "normal_stress",
   "param_a", "param_Dc", "fault_x2", "fault_x3",
};

static void RunScenarioToggleOff(Mesh &mesh, int fault_face)
{
   std::cout << "\n=== Scenario A: toggle DEFAULT OFF ===\n";
   ParaViewOutput<Mesh> pv("test_kinematics_field_set_out_A", mesh,
                           /*order=*/1, /*collection_name=*/"kinematics");

   // Driver-style explicit registrations (mimics what tpv*_driver.cpp
   // does at line ~1933).
   L2_FECollection vel_fec(1, mesh.Dimension(), BasisType::GaussLobatto);
   FiniteElementSpace vel_fes(&mesh, &vel_fec,
                              mesh.SpaceDimension(),
                              Ordering::byNODES);
   GridFunction velocity(&vel_fes);
   velocity = 0.0;
   pv.RegisterDomainField("velocity", &velocity);

   L2_FECollection rank_fec(0, mesh.Dimension());
   FiniteElementSpace rank_fes(&mesh, &rank_fec);
   GridFunction mpi_rank(&rank_fes);
   mpi_rank = 0.0;
   pv.RegisterDomainField("mpi_rank", &mpi_rank);

   // Default toggle is OFF — no opt-in call.
   TEST_TRUE(pv.GetRegisterFaultProjectionsInVolumePV() == false,
             "default `register_fault_projections_in_volume_pv_` is false");

   // Construct dummy fault-face arrays.  One interior fault face, no
   // shared faces (serial mesh).
   Array<int> fault_int_faces(1);
   fault_int_faces[0] = fault_face;
   Array<int> fault_shr_faces;
   pv.InitFaultOutputBP5(fault_int_faces, fault_shr_faces,
                         /*nbf_per_face=*/4);

   const auto *dc = pv.GetVolumeDataCollection();
   TEST_TRUE(dc != nullptr,
             "GetVolumeDataCollection() returns non-null");
   if (!dc) { return; }

   // Required fields:
   TEST_TRUE(dc->HasField("velocity"),
             "kinematics.vtkhdf contains `velocity`");
   TEST_TRUE(dc->HasField("mpi_rank"),
             "kinematics.vtkhdf contains `mpi_rank`");

   // FORBIDDEN fields (the 12 L2-p0 projections — duplicates of
   // fault.vtkhdf, scheduled to NOT appear in kinematics.vtkhdf by
   // default per PLAN_split_bulk_solutions_2026-05-12).
   for (const auto &name : kFaultProjectionNames)
   {
      const bool present = dc->HasField(name);
      TEST_TRUE(!present,
                "kinematics.vtkhdf does NOT contain `" + name +
                "` when toggle is OFF (default)");
   }
}

static void RunScenarioToggleOn(Mesh &mesh, int fault_face)
{
   std::cout << "\n=== Scenario B: toggle EXPLICITLY ON ===\n";
   ParaViewOutput<Mesh> pv("test_kinematics_field_set_out_B", mesh,
                           /*order=*/1, /*collection_name=*/"kinematics");

   L2_FECollection vel_fec(1, mesh.Dimension(), BasisType::GaussLobatto);
   FiniteElementSpace vel_fes(&mesh, &vel_fec,
                              mesh.SpaceDimension(),
                              Ordering::byNODES);
   GridFunction velocity(&vel_fes);
   velocity = 0.0;
   pv.RegisterDomainField("velocity", &velocity);

   // Opt in BEFORE InitFaultOutputBP5 (registrations happen inside
   // that call; setting the toggle after would be too late).
   pv.SetRegisterFaultProjectionsInVolumePV(true);
   TEST_TRUE(pv.GetRegisterFaultProjectionsInVolumePV() == true,
             "toggle is now ON after the setter");

   Array<int> fault_int_faces(1);
   fault_int_faces[0] = fault_face;
   Array<int> fault_shr_faces;
   pv.InitFaultOutputBP5(fault_int_faces, fault_shr_faces,
                         /*nbf_per_face=*/4);

   const auto *dc = pv.GetVolumeDataCollection();
   TEST_TRUE(dc != nullptr,
             "GetVolumeDataCollection() returns non-null");
   if (!dc) { return; }

   TEST_TRUE(dc->HasField("velocity"),
             "kinematics.vtkhdf contains `velocity` (opt-in path)");

   // Now all 12 must be present.
   for (const auto &name : kFaultProjectionNames)
   {
      TEST_TRUE(dc->HasField(name),
                "kinematics.vtkhdf contains `" + name +
                "` when toggle is ON");
   }
}

int main(int /*argc*/, char ** /*argv*/)
{
   std::cout << "test_kinematics_field_set: PLAN_split_bulk_solutions_"
                "2026-05-12 acceptance test for R-001\n";
   Mesh mesh = BuildTwoHexMesh();
   const int fault_face = FindInteriorFace(mesh);
   if (fault_face < 0)
   {
      std::cerr << "FATAL: synthetic two-hex mesh has no interior face; "
                   "test cannot proceed.\n";
      return 2;
   }

   RunScenarioToggleOff(mesh, fault_face);
   RunScenarioToggleOn(mesh, fault_face);

   std::cout << "\n=== Summary: " << num_passed << " / " << num_tests
             << " passed; " << num_failed << " failed ===\n";
   return num_failed == 0 ? 0 : 1;
}
