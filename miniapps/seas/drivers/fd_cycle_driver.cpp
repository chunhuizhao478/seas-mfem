// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fd_cycle_driver.cpp — fully dynamic BP5 cycle driver, seas_fd_cycle_driver
// (plan rev 3, item j).  Quasi-dynamic between earthquakes, Nitsche
// elastodynamics during them, switched on the slip rate.
//
// Phase P0 stub.  It only loads the mesh and reports what the later phases
// rely on:
//   - the commit the binary was built from (cycle/build_info.hpp, rule M10);
//   - per boundary attribute, the number of faces, how many are exterior (one
//     element) and interior (two elements), and their bounding box after
//     scaling;
//   - the fault faces (--fault-attr) and the fault DOFs of the IP layout,
//     (p+1)(p+2)/2 nodal points per triangular face, and after the partition
//     how many fault faces each rank holds and how many lie on a face shared
//     between ranks.
// The counts are mesh-level; the operators' own counts are printed from P1
// on.  No solver header is included here (rule M4): QDPhase arrives in P1,
// DynPhase in P2.
//
// Usage (from miniapps/seas, built with make -f cycle/cycle.mk):
//   mpirun -np N ./seas_fd_cycle_driver \
//      --mesh bp5/mesh/reference/bp5_tandem_exact.msh --mesh-scale 1000 \
//      --order 1 [--fault-attr 3]

#include "mfem.hpp"
#include "../cycle/build_info.hpp"
#include "../cycle/scaled_mesh.hpp"

#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace mfem;
using namespace mfem::seas;

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   const int rank = Mpi::WorldRank();
   const int nranks = Mpi::WorldSize();

   const char *mesh_file = "bp5/mesh/reference/bp5_tandem_exact.msh";
   real_t mesh_scale = 1000.0;
   int order = 1;
   int fault_attr = 3;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file (km for the BP5 mesh).");
   args.AddOption(&mesh_scale, "-s", "--mesh-scale",
                  "Vertex scale applied before the partition (1000: km -> m).");
   args.AddOption(&order, "-o", "--order", "Polynomial order p of the DG space.");
   args.AddOption(&fault_attr, "-fa", "--fault-attr", "Boundary attribute of the fault.");
   args.Parse();
   if (!args.Good())
   {
      if (rank == 0) { args.PrintUsage(std::cout); }
      return 1;
   }
   MFEM_VERIFY(order >= 1, "--order must be >= 1");

   if (rank == 0)
   {
      std::cout << "seas_fd_cycle_driver (P0 stub)\n"
                << "  commit: " << cycle::BuildCommit()
                << (cycle::BuildIsClean() ? "" : "   [not a clean build of a commit]")
                << "\n"
                << "  ranks: " << nranks << "\n"
                << "  mesh: " << mesh_file << " x " << mesh_scale << "\n"
                << "  order: " << order << "\n";
   }

   std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(mesh_file, mesh_scale);
   const std::map<int, cycle::BdrAttrCount> bdr = cycle::CountBoundary(*smesh);
   const long long ne_global = smesh->GetNE();

   long long fault_faces = 0;
   if (rank == 0)
   {
      std::cout << "  elements: " << ne_global << ", vertices: " << smesh->GetNV()
                << ", boundary faces: " << smesh->GetNBE() << "\n"
                << "  boundary attributes (faces = exterior + interior; bounding box, m):\n";
      for (const auto &kv : bdr)
      {
         const cycle::BdrAttrCount &c = kv.second;
         std::cout << "    attr " << std::setw(2) << c.attr << ": " << std::setw(6) << c.faces
                   << " = " << std::setw(6) << c.exterior << " + " << std::setw(6) << c.interior
                   << std::setprecision(6)
                   << "   x [" << c.lo[0] << ", " << c.hi[0] << "]"
                   << " y [" << c.lo[1] << ", " << c.hi[1] << "]"
                   << " z [" << c.lo[2] << ", " << c.hi[2] << "]"
                   << (c.attr == fault_attr ? "   <- fault" : "") << "\n";
      }
   }
   const auto it = bdr.find(fault_attr);
   MFEM_VERIFY(it != bdr.end(), "no boundary faces carry --fault-attr " << fault_attr);
   const cycle::BdrAttrCount &fc = it->second;
   MFEM_VERIFY(fc.exterior == 0, fc.exterior << " fault faces (attribute " << fault_attr
               << ") are on the outer boundary; a fault face must be interior");
   MFEM_VERIFY(fc.triangles == fc.faces, "the IP fault layout needs triangular fault "
               "faces; " << fc.faces - fc.triangles << " are not triangles");
   fault_faces = fc.faces;
   const int nbf = cycle::FaultDOFsPerFace(order);

   ParMesh pmesh(MPI_COMM_WORLD, *smesh);
   smesh.reset();

   const cycle::LocalFaultFaces loc = cycle::CountLocalFaultFaces(pmesh, fault_attr);
   long long sums[3] = {loc.faces, loc.shared, pmesh.GetNE()};
   MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   long long fmin = loc.faces, fmax = loc.faces;
   MPI_Allreduce(MPI_IN_PLACE, &fmin, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &fmax, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
   MFEM_VERIFY(sums[0] == fault_faces && sums[2] == ne_global,
               "partition lost faces or elements: fault faces " << sums[0] << " vs "
               << fault_faces << ", elements " << sums[2] << " vs " << ne_global);

   if (rank == 0)
   {
      std::cout << "  fault: " << fault_faces << " faces x " << nbf
                << " DOFs per face (p = " << order << ") = " << fault_faces * nbf
                << " fault DOFs\n"
                << "  partition: " << fmin << " .. " << fmax
                << " fault faces per rank; " << sums[1]
                << " fault faces on faces shared between ranks\n"
                << "P0 stub done.\n";
   }
   return 0;
}
