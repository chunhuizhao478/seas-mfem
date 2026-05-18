// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_paraview_set_fault_params_spatial.cpp — R-005 round-3 regression
// for the new `seas::ParaViewOutput::SetFaultParamsSpatial(...)` method
// added in Parity Phase 6.
//
// Asserts that calling `SetFaultParamsSpatial` with the 8 SAFS-semantic
// named arrays registers each one under its EXACT name (`lsw_mu_s`,
// `lsw_mu_d`, `lsw_d_c`, `nuc_amplitude`, `nuc_radial_factor`,
// `sigma_n_init`, `tau1_init`, `tau2_init`) on the underlying
// ParaView DataCollection — and does NOT mislabel them as TPV104's
// `param_a` / `param_Dc` / `fault_x2` / `fault_x3` (the round-2 R-002
// bug regression).
//
// Uses a 1-element ParMesh + empty fault-face lists so the L2-p0 face-
// average loop is a no-op; the test focuses on the registration
// pathway, not the averaging math.

#include "mfem.hpp"
#include "../../io/paraview_output.hpp"

#include <iostream>
#include <string>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

// =====================================================================
// R-005 / Parity Phase 6: SetFaultParamsSpatial registers 8 SAFS fields
// under SAFS-semantic names, NOT TPV104's a/Dc/x2/x3.
// =====================================================================
static void R005_SetFaultParamsSpatial_registers_8_safs_fields()
{
   std::cout << "\n[R-005] SetFaultParamsSpatial registers 8 SAFS fields\n";

   // 1-element hex ParMesh.  No fault DOFs (empty face lists below).
   Mesh smesh = Mesh::MakeCartesian3D(1, 1, 1, Element::HEXAHEDRON);
   ParMesh pmesh(MPI_COMM_WORLD, smesh);
   smesh.Clear();

   // Build the ParaViewOutput in VTU mode so we don't need HDF5 at link
   // time.  collection_name = "volume" matches the spatial driver.
   ParaViewOutput<ParMesh> pv("/tmp", pmesh, /*order=*/1,
                              "volume",
                              ParaViewOutput<ParMesh>::VolumeOutputMode::Vtu);

   // Empty fault face lists — InitFaultOutputBP5 still allocates the
   // L2-p0 FES and sets has_fault_output_ = true, so SetFaultParamsSpatial
   // can proceed.
   Array<int> empty_int, empty_shr;
   pv.InitFaultOutputBP5(empty_int, empty_shr, /*nbf_per_face=*/1);

   // SetFaultParamsSpatial expects per-DOF arrays of size
   // (n_int + n_shr) * nbf = 0.  Empty Vectors are accepted by the
   // "n_dofs_local == 0 && data->Size() == 0" branch of the validator.
   Vector empty(0);
   pv.SetFaultParamsSpatial({
      {"lsw_mu_s",          &empty},
      {"lsw_mu_d",          &empty},
      {"lsw_d_c",           &empty},
      {"nuc_amplitude",     &empty},
      {"nuc_radial_factor", &empty},
      {"sigma_n_init",      &empty},
      {"tau1_init",         &empty},
      {"tau2_init",         &empty},
   });

   // Inspect registered field names on the underlying collection.
   const auto* dc = pv.GetVolumeDataCollection();
   TEST_ASSERT(dc != nullptr, "pv_dc_ non-null after InitFaultOutputBP5");
   if (dc == nullptr) { return; }

   // All 8 SAFS-semantic names must be present.
   TEST_ASSERT(dc->HasField("lsw_mu_s"),
               "lsw_mu_s registered under its own name");
   TEST_ASSERT(dc->HasField("lsw_mu_d"),
               "lsw_mu_d registered under its own name");
   TEST_ASSERT(dc->HasField("lsw_d_c"),
               "lsw_d_c registered under its own name");
   TEST_ASSERT(dc->HasField("nuc_amplitude"),
               "nuc_amplitude registered under its own name");
   TEST_ASSERT(dc->HasField("nuc_radial_factor"),
               "nuc_radial_factor registered under its own name");
   TEST_ASSERT(dc->HasField("sigma_n_init"),
               "sigma_n_init registered under its own name");
   TEST_ASSERT(dc->HasField("tau1_init"),
               "tau1_init registered under its own name");
   TEST_ASSERT(dc->HasField("tau2_init"),
               "tau2_init registered under its own name");

   // Critically: the SAFS arrays must NOT be mislabelled as TPV104's
   // BP5 names.  SetFaultParamsBP5 is NOT called by the spatial driver
   // (R-002 round-2 fix), so the BP5-named fields are only present if
   // register_fault_projections_in_volume_pv_ was set (default false).
   // Verify they are absent in this default-OFF configuration.
   TEST_ASSERT(!dc->HasField("param_a"),
               "lsw_mu_s NOT mislabelled as 'param_a' (R-002 regression)");
   TEST_ASSERT(!dc->HasField("param_Dc"),
               "lsw_d_c NOT mislabelled as 'param_Dc' (R-002 regression)");
   TEST_ASSERT(!dc->HasField("fault_x2"),
               "tau1_init NOT mislabelled as 'fault_x2' (R-002 regression)");
   TEST_ASSERT(!dc->HasField("fault_x3"),
               "tau2_init NOT mislabelled as 'fault_x3' (R-002 regression)");
}

// =====================================================================
// Idempotence: a SECOND call to SetFaultParamsSpatial with the same
// names should NOT double-register (each name maps to ONE GF).
// =====================================================================
static void R005_SetFaultParamsSpatial_idempotent_on_repeat_call()
{
   std::cout << "\n[R-005-IDEMPOTENT] repeat call does not double-register\n";
   Mesh smesh = Mesh::MakeCartesian3D(1, 1, 1, Element::HEXAHEDRON);
   ParMesh pmesh(MPI_COMM_WORLD, smesh);
   smesh.Clear();
   ParaViewOutput<ParMesh> pv("/tmp", pmesh, /*order=*/1, "volume",
                              ParaViewOutput<ParMesh>::VolumeOutputMode::Vtu);
   Array<int> empty_int, empty_shr;
   pv.InitFaultOutputBP5(empty_int, empty_shr, /*nbf_per_face=*/1);
   Vector empty(0);
   const std::vector<std::pair<std::string, const Vector*>> named = {
      {"lsw_mu_s",          &empty},
      {"lsw_mu_d",          &empty},
      {"lsw_d_c",           &empty},
      {"nuc_amplitude",     &empty},
      {"nuc_radial_factor", &empty},
      {"sigma_n_init",      &empty},
      {"tau1_init",         &empty},
      {"tau2_init",         &empty},
   };
   pv.SetFaultParamsSpatial(named);
   const auto* dc = pv.GetVolumeDataCollection();
   const int n_after_first = static_cast<int>(dc->GetFieldMap().size());
   pv.SetFaultParamsSpatial(named);
   const int n_after_second = static_cast<int>(dc->GetFieldMap().size());
   TEST_ASSERT(n_after_first == n_after_second,
               "field count unchanged across repeat SetFaultParamsSpatial");
   TEST_ASSERT(dc->HasField("lsw_mu_s"),
               "lsw_mu_s still registered after repeat call");
   TEST_ASSERT(dc->HasField("tau2_init"),
               "tau2_init still registered after repeat call");
}

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
#else
   (void)argc; (void)argv;
#endif

   std::cout << "Running test_paraview_set_fault_params_spatial\n";
   R005_SetFaultParamsSpatial_registers_8_safs_fields();
   R005_SetFaultParamsSpatial_idempotent_on_repeat_call();

   std::cout << "\n========================================\n";
   std::cout << "test_paraview_set_fault_params_spatial: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";

#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return num_failed == 0 ? 0 : 1;
}
