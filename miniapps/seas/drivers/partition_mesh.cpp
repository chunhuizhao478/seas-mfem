// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// seas_partition_mesh — pre-partition a (large) serial gmsh .msh into MFEM
// per-rank parallel mesh files, so production runs can read their own
// partition directly via ParMesh(MPI_Comm, std::istream&) WITHOUT every rank
// building the full serial mesh.
//
// Why this exists.  The serial mesh-read path
//   Mesh smesh(file);  ParMesh pmesh(comm, smesh);   (spatial_dyn_driver.cpp)
// requires the COMPLETE serial mesh on EVERY rank, partitioned in-memory with
// METIS.  Per-node startup memory therefore scales as
// (ranks_per_node x total_elements).  On the 250 m triqsubdiv mesh (16 M
// elements) at ~38-50 ranks/node this exhausts node RAM and aborts with
// std::bad_alloc before time-stepping.  Pre-splitting moves that one-time
// replicated build into a single short job that can be spread over enough
// nodes to fit, after which production reads ~1/np of the mesh per rank.
//
// This writes, for an MPI run of size np:
//   <out>.000000 ... <out>.<np-1>   one ParPrint() file per rank
//   <out>.np                        manifest: the integer np (read-side guard)
//
// It MUST be run at the SAME rank count that production will use (the per-rank
// files and the baked partition are valid ONLY at that np), but spread over
// ENOUGH nodes that the one-time replicated build fits in node memory (use few
// ranks/node).  Production then runs:
//   ibrun -n <np> seas_spatial_dyn_driver --mesh-mode presplit --mesh <out> ...
//
// File format note.  We use ParMesh::ParPrint (NOT Save/Print, which write the
// mfem v1.0 *visualization* format and are NOT round-trippable into a ParMesh).
// ParPrint terminates the serial section with 'mfem_serial_mesh_end' and
// appends the parallel topology; ParMesh(comm, istream) reads exactly that.
// This mirrors MFEM's own checkpoint idiom (examples/ex6p.cpp) and the parallel
// DataCollection round-trip (fem/datacollection.cpp: ParPrint -> ParMesh(comm,
// file)).
//
// Usage (conda/Frontera: built by the seas Makefile):
//   ibrun -n 400 seas_partition_mesh --mesh BIG.msh --out path/to/prefix
//
// The optional finite-element order is NOT applied here: SetCurvature is left
// to the read side (spatial_dyn_driver) so the serial and presplit paths build
// an identical ParMesh and then curve it identically.

#include "mfem.hpp"

#include <fstream>
#include <iostream>
#include <string>

using namespace mfem;

int main(int argc, char *argv[])
{
   // Raw MPI_Init to match the other seas drivers (mfem::Mpi::Init +
   // Hypre::Init segfaults on this MFEM build path; see
   // drivers/project_velocity_to_mesh.cpp).
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);

   const char *mesh_path  = "";
   const char *out_prefix = "";
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_path, "-m", "--mesh",
                  "Serial input mesh (gmsh .msh, required).");
   args.AddOption(&out_prefix, "-o", "--out",
                  "Output per-rank prefix (required).  Writes <out>.NNNNNN "
                  "(one ParPrint file per rank) plus <out>.np (manifest).");
   args.Parse();
   if (!args.Good() || std::string(mesh_path).empty()
       || std::string(out_prefix).empty())
   {
      if (rank == 0) { args.PrintUsage(std::cout); }
      MPI_Finalize();
      return 1;
   }

   if (rank == 0)
   {
      std::cout << "seas_partition_mesh: reading serial mesh\n  " << mesh_path
                << "\n  and partitioning into " << nprocs << " parts -> '"
                << out_prefix << ".NNNNNN'\n";
   }

   // --- Build the ParMesh the replicated way.  This is the memory-heavy step
   //     we are intentionally confining to one short, wide (few-ranks/node)
   //     job.  Identical to spatial_dyn_driver.cpp's serial path. ---
   Mesh smesh(mesh_path, 1, 1);
   MFEM_VERIFY(smesh.Dimension() == 3,
               "seas_partition_mesh: only 3D meshes supported; got dim="
               << smesh.Dimension());
   ParMesh pmesh(comm, smesh);
   smesh.Clear();

   // --- Write THIS rank's partition in the round-trippable parallel format. ---
   const std::string fname = MakeParFilename(std::string(out_prefix) + ".",
                                             rank);
   std::ofstream ofs(fname);
   MFEM_VERIFY(ofs.good(),
               "seas_partition_mesh: cannot open output file '" << fname
               << "' (rank " << rank << ").  Does the output directory exist "
               "and is it writable?");
   ofs.precision(16);
   pmesh.ParPrint(ofs);
   ofs.flush();
   MFEM_VERIFY(ofs.good(),
               "seas_partition_mesh: error while writing '" << fname << "'.");
   ofs.close();

   // --- rank 0 writes the manifest AFTER every rank has finished its file, so
   //     the manifest's presence implies a complete partition set. ---
   MPI_Barrier(comm);
   if (rank == 0)
   {
      const std::string manifest = std::string(out_prefix) + ".np";
      std::ofstream mfs(manifest);
      MFEM_VERIFY(mfs.good(),
                  "seas_partition_mesh: cannot open manifest '" << manifest
                  << "'.");
      mfs << nprocs << "\n";
      mfs.flush();
      MFEM_VERIFY(mfs.good(),
                  "seas_partition_mesh: error while writing manifest '"
                  << manifest << "'.");
      mfs.close();
      std::cout << "seas_partition_mesh: wrote " << nprocs
                << " partition files + manifest '" << manifest << "'.\n"
                << "  Run production at EXACTLY -n " << nprocs << ":\n"
                << "    seas_spatial_dyn_driver --mesh-mode presplit --mesh "
                << out_prefix << " ...\n";
   }

   MPI_Finalize();
   return 0;
}
