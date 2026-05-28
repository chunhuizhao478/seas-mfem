// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Round-trip test for the pre-split mesh path used by seas_partition_mesh and
// seas_spatial_dyn_driver --mesh-mode presplit.
//
// It exercises the exact MFEM round-trip the feature relies on:
//     ParMesh(comm, serial_mesh)  ->  ParPrint(file)  ->  ParMesh(comm, file)
// and asserts the reloaded ParMesh has the SAME global element count, global
// boundary-element count, dimension, and attribute sets (including a
// fault-like boundary attribute 101) as the partitioned-in-memory mesh.  This
// guards the core correctness claim: pre-splitting changes only WHERE the
// partition is built, not the resulting mesh.
//
// Run with: mpirun -np 2 ./seas_test_presplit_mesh
//       or: mpirun -np 4 ./seas_test_presplit_mesh

#include "mfem.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

using namespace mfem;

namespace
{
// Sum a local count across ranks (returns the global total on all ranks).
long long GlobalSum(long long local, MPI_Comm comm)
{
   long long g = local;
   MPI_Allreduce(MPI_IN_PLACE, &g, 1, MPI_LONG_LONG, MPI_SUM, comm);
   return g;
}

// Max of a local integer across ranks.
int GlobalMax(int local, MPI_Comm comm)
{
   int g = local;
   MPI_Allreduce(MPI_IN_PLACE, &g, 1, MPI_INT, MPI_MAX, comm);
   return g;
}
} // namespace

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);

   int fails = 0;
   auto check = [&](const char *name, bool ok)
   {
      if (!ok)
      {
         fails++;
         if (rank == 0) { std::cerr << "FAIL: " << name << "\n"; }
      }
   };

   // --- serial reference mesh, with a fault-like boundary attribute (101) so
   //     the test pins attribute preservation, not just counts. ---
   Mesh smesh = Mesh::MakeCartesian3D(6, 6, 6, Element::TETRAHEDRON);
   smesh.GetBdrElement(0)->SetAttribute(101);
   smesh.SetAttributes();                 // rebuild element/bdr attribute tables
   const int dim_ref = smesh.Dimension();

   // --- reference invariants from the partitioned-in-memory ParMesh ---
   ParMesh pmesh(comm, smesh);
   smesh.Clear();
   const long long gne_ref  = pmesh.GetGlobalNE();
   const long long gnbe_ref = GlobalSum(pmesh.GetNBE(), comm);
   const int battr_ref =
      GlobalMax(pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0,
                comm);
   const int eattr_ref =
      GlobalMax(pmesh.attributes.Size() ? pmesh.attributes.Max() : 0, comm);

   // --- write this rank's partition exactly as seas_partition_mesh does ---
   const std::string fname = MakeParFilename("test_presplit_tmp.", rank);
   {
      std::ofstream ofs(fname.c_str());
      check("open ParPrint output file", ofs.good());
      ofs.precision(16);
      pmesh.ParPrint(ofs);
   }
   MPI_Barrier(comm);                      // all files complete before reading

   // --- read back via the istream ctor exactly as --mesh-mode presplit does ---
   long long gne_rt = -1, gnbe_rt = -1;
   int battr_rt = -1, eattr_rt = -1, dim_rt = -1;
   {
      std::ifstream ifs(fname.c_str());
      check("open ParPrint input file", ifs.good());
      ParMesh pmesh2(comm, ifs);
      dim_rt   = pmesh2.Dimension();
      gne_rt   = pmesh2.GetGlobalNE();
      gnbe_rt  = GlobalSum(pmesh2.GetNBE(), comm);
      battr_rt = GlobalMax(
                    pmesh2.bdr_attributes.Size() ? pmesh2.bdr_attributes.Max() : 0,
                    comm);
      eattr_rt = GlobalMax(
                    pmesh2.attributes.Size() ? pmesh2.attributes.Max() : 0, comm);
   }

   check("dimension preserved",                 dim_rt   == dim_ref);
   check("global element count preserved",      gne_rt   == gne_ref);
   check("global boundary count preserved",     gnbe_rt  == gnbe_ref);
   check("fault-like bdr attr 101 preserved",
         battr_rt == battr_ref && battr_ref == 101);
   check("max element attribute preserved",     eattr_rt == eattr_ref);

   std::remove(fname.c_str());             // each rank cleans its own temp file

   fails = static_cast<int>(GlobalSum(fails, comm));
   if (rank == 0)
   {
      if (fails == 0)
      {
         std::cout << "PASS: presplit mesh round-trip (np=" << nprocs
                   << "): GlobalNE=" << gne_ref << " GlobalNBE=" << gnbe_ref
                   << " fault_attr=" << battr_ref << "\n";
      }
      else
      {
         std::cerr << fails << " presplit round-trip check(s) FAILED (np="
                   << nprocs << ").\n";
      }
   }
   MPI_Finalize();
   return fails == 0 ? 0 : 1;
}
