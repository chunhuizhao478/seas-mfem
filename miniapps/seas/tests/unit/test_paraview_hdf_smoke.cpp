// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 2a smoke test (PLAN_paraview_compaction_2026-04-28.md §Phase 2a.6).
//
// Verifies MFEM's ParaViewHDFDataCollection can write a minimal one-element
// mesh + one L2-p0 GridFunction to a single .vtkhdf file and that the
// resulting file's first 8 bytes match the HDF5 superblock magic.  No
// SEAS-specific code paths are exercised — this test isolates the build
// configuration change (MFEM_USE_HDF5=YES + parallel HDF5 link) from the
// fault-writer integration in Phase 2b.

#include "mfem.hpp"

#include <fstream>
#include <iostream>

#ifdef MFEM_USE_HDF5
#include <hdf5.h>     // H5_VERS_*
#endif

int main(int argc, char *argv[])
{
   mfem::Mpi::Init(argc, argv);
   const int rank = mfem::Mpi::WorldRank();
#ifndef MFEM_USE_HDF5
   // R-010 (REVIEW.md 2026-04-28): on a build without MFEM_USE_HDF5,
   // there is nothing to test.  Exit 0 (SKIP) so the Makefile target
   // does not flag the run as failed.  Plan §Constraints commits to
   // "Phase 2 code is #ifdef MFEM_USE_HDF5-guarded so a non-HDF5
   // build still compiles and runs".
   if (rank == 0)
   {
      std::cout << "SKIP: MFEM_USE_HDF5 is OFF in this build — Phase 2a "
                   "build prerequisite has not been applied.\n";
   }
   mfem::Mpi::Finalize();
   return 0;
#else
   if (rank == 0)
   {
      std::cout << "=== test_paraview_hdf_smoke (Phase 2a gate) ===\n";
      std::cout << "HDF5 version compiled against: "
                << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "."
                << H5_VERS_RELEASE << "\n";
   }

   // Minimal mesh: one tetrahedron from MakeCartesian3D, then promoted
   // to ParMesh.  Single L2-p0 grid function on the volume.
   mfem::Mesh smesh = mfem::Mesh::MakeCartesian3D(
      /*nx=*/1, /*ny=*/1, /*nz=*/1, mfem::Element::TETRAHEDRON);
   mfem::ParMesh mesh(MPI_COMM_WORLD, smesh);
   mfem::L2_FECollection fec(0, mesh.Dimension());
   mfem::ParFiniteElementSpace fes(&mesh, &fec);
   mfem::ParGridFunction gf(&fes);
   gf = 3.14;

   const std::string out_dir = "/tmp/test_paraview_hdf_smoke_out";
   const std::string vtkhdf_path = out_dir + "/smoke.vtkhdf";

   {
      mfem::ParaViewHDFDataCollection dc("smoke", &mesh);
      dc.SetPrefixPath(out_dir);
      dc.SetDataFormat(mfem::VTKFormat::BINARY);
      dc.SetCompression(true);
      dc.SetCompressionLevel(3);
      dc.RegisterField("scalar", &gf);
      dc.SetCycle(0);
      dc.SetTime(0.0);
      dc.Save();
   }

   int local_failed = 0;
   if (rank == 0)
   {
      std::ifstream f(vtkhdf_path, std::ios::binary);
      if (!f.is_open())
      {
         std::cerr << "FAIL: " << vtkhdf_path << " not created\n";
         local_failed = 1;
      }
      else
      {
         f.close();
         // R-009 (REVIEW.md 2026-04-28): use HDF5's own predicate
         // instead of an 8-byte signature read.  H5Fis_hdf5 handles
         // user-block padding (up to 512 bytes prepended before the
         // superblock) which a raw byte-position-0 read would miss.
         const htri_t is_h5 = H5Fis_hdf5(vtkhdf_path.c_str());
         if (is_h5 <= 0)
         {
            std::cerr << "FAIL: H5Fis_hdf5(" << vtkhdf_path
                      << ") returned " << is_h5 << "\n";
            local_failed = 1;
         }
         else
         {
            std::cout << "PASSED: " << vtkhdf_path
                      << " written with valid HDF5 superblock "
                         "(H5Fis_hdf5 == 1)\n";
         }
      }
   }

   int total_failed = 0;
   MPI_Allreduce(&local_failed, &total_failed, 1, MPI_INT, MPI_MAX,
                 MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::cout << "=== Result: "
                << (total_failed == 0 ? "PASS" : "FAIL")
                << " ===\n";
   }
   mfem::Mpi::Finalize();
   return total_failed == 0 ? 0 : 1;
#endif
}
