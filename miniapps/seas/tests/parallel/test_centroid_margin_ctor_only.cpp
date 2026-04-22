// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan §17 rev-3d (post-REVIEW R-V92-C04) —
// ctor-only centroid-margin probe.
//
// ============================================================================
// Purpose
// ============================================================================
// Loads a mesh file, partitions to ParMesh, instantiates WaveOperator<ParMesh>
// with SEAS_DIAG_CENTROID_MARGIN compiled in, and EXITS immediately after
// ctor.  No time-stepping, no Mult calls, no simulation.  The WaveOperator
// ctor emits one stderr line per shared fault face with the FP margin of
// the `elem1_on_plus` comparison.
//
// After the run, post-process the stderr with `grep` / `awk` to extract
// the distribution of margin values.  A histogram that extends to O(ε_FP)
// indicates FP-fragility at those faces — H-V92-P mechanism CONFIRMED at
// the production-mesh geometry.
//
// Per feedback_no_local_reproducer.md: "local unit tests OK only on tiny
// fixtures, never the production mesh".  This probe does NOT run the
// driver's time-stepping (no `Mult` call), so it is categorically
// ctor-only + I/O.  Safe to run on the 1000 m mesh locally to measure
// the FP margin distribution before requesting Frontera for 200 m.
//
// Usage:
//   mpirun -np 4 ./seas_test_centroid_margin_ctor_only \
//       --mesh tpv102/mesh/tpv102_1000m.msh    2> margins_1000m.log
//   grep '\[CENTROID-MARGIN\]' margins_1000m.log | awk '{...}'
//
// Build with the flag:
//   make seas_test_centroid_margin_ctor_only \
//       SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_CENTROID_MARGIN"

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <iostream>
#include <string>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static std::string GetArg(int argc, char **argv, const std::string &flag,
                          const std::string &def)
{
   for (int i = 1; i < argc - 1; i++)
   { if (flag == argv[i]) { return argv[i+1]; } }
   return def;
}
static int GetIntArg(int argc, char **argv, const std::string &flag, int def)
{
   for (int i = 1; i < argc - 1; i++)
   { if (flag == argv[i]) { return std::atoi(argv[i+1]); } }
   return def;
}

int main(int argc, char **argv)
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   const std::string mesh_file = GetArg(argc, argv, "--mesh",
      "tpv102/mesh/tpv102_1000m.msh");
   const int order      = GetIntArg(argc, argv, "--order", 1);
   const int fault_attr = GetIntArg(argc, argv, "--fault-attr", 3);

   if (rank == 0)
   {
      std::cout << "=== TPV102 v9.2.0 §17 CENTROID_MARGIN (ctor-only) ===\n";
      std::cout << "  mesh        : " << mesh_file << "\n";
      std::cout << "  ranks       : " << nprocs << "\n";
      std::cout << "  order       : " << order << "\n";
      std::cout << "  fault_attr  : " << fault_attr << "\n";
   }

   Mesh serial_mesh(mesh_file.c_str(), 1, 1);
   MFEM_VERIFY(serial_mesh.Dimension() == 3, "Need a 3D mesh");

   ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = fault_attr;
   bc.absorbing_attrs = {5};

   // Instantiate WaveOperator — emits [CENTROID-MARGIN] lines on stderr.
   WaveOperator<ParMesh> wave(pmesh, order,
                              TPV102Params::lambda,
                              TPV102Params::mu,
                              TPV102Params::rho, bc);

   const int n_shr = wave.GetFaultSharedFaces().Size();
   int tot_shr = 0;
   MPI_Allreduce(&n_shr, &tot_shr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::cout << "  total shared fault faces (summed over ranks): "
                << tot_shr << "\n";
      std::cout << "  see stderr for [CENTROID-MARGIN] lines (requires "
                << "-DSEAS_DIAG_CENTROID_MARGIN at build time)\n";
      std::cout << "  post-process: grep '\\[CENTROID-MARGIN\\]' <stderr> | "
                << "awk '{print $4}' | sort -g | head -5   "
                << "(smallest margins)\n";
   }

   MPI_Finalize();
   return 0;
}
