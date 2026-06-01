// Parallel (ParMesh) DerivMode Cached-vs-OnTheFly equivalence test.
// REVIEW R-002 (Lever 1/3): the cached ApplySpatialDerivative (Lever 1) and
// ComputeVolumeRHS (Lever 3) kernels are element-LOCAL (no MPI), so on a
// partitioned ParMesh each rank runs them on its local elements.  Assert
// cached == OnTheFly per rank to <= 1e-12 (machine-eps, NOT bit-exact).
//
// Standalone (does NOT use Mult / the fault flux), so it needs no absorbing
// background or fault setup — it exercises exactly the two cached kernels in a
// partitioned, ghost-face topology.
//
// Run with: mpirun -np 2 ./seas_test_wave_operator_cached_parallel

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../domain/boundary_config.hpp"

#include <iostream>
#include <cmath>
#include <algorithm>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);

   if (rank == 0)
   {
      std::cout << "=== Parallel DerivMode Cached==OnTheFly (np=" << nprocs
                << ") ===\n";
   }

   Mesh serial_mesh = Mesh::MakeCartesian3D(4, 4, 4, Element::TETRAHEDRON,
                                            1.0, 1.0, 1.0);
   for (int b = 0; b < serial_mesh.GetNBE(); b++)
   {
      serial_mesh.SetBdrAttribute(b, 5);
   }
   serial_mesh.SetAttributes();

   ParMesh pmesh(comm, serial_mesh);
   BoundaryConfig bc;
   bc.absorbing_attrs = {5};
   bc.fault_attr = 0;   // no fault (the cached kernels do not need one)

   const int order = 2;
   WaveOperator<ParMesh> wave(pmesh, order, 32.04e9, 32.04e9, 2670.0, bc);

   const int N = wave.Height();
   Vector Q(N);
   // Deterministic, rank-shifted fill so ranks carry distinct local state.
   unsigned s = 12345u + 7919u * static_cast<unsigned>(rank);
   for (int i = 0; i < N; i++)
   {
      s = 1664525u * s + 1013904223u;
      Q[i] = (static_cast<double>(s) / 4294967296.0) - 0.5;
   }

   auto rel = [](const Vector &a, const Vector &b)
   {
      double mr = 0.0, md = 0.0;
      for (int i = 0; i < a.Size(); i++)
      {
         mr = std::max(mr, std::abs(a[i]));
         md = std::max(md, std::abs(a[i] - b[i]));
      }
      return md / (mr + 1e-300);
   };

   double local_max = 0.0;

   // Lever 1: ApplySpatialDerivative, all three directions.
   for (int d = 0; d < 3; d++)
   {
      wave.SetDerivMode(DerivMode::OnTheFly);
      Vector r1;
      wave.ApplySpatialDerivative(d, Q, r1);
      wave.SetDerivMode(DerivMode::Cached);
      Vector r2;
      wave.ApplySpatialDerivative(d, Q, r2);
      local_max = std::max(local_max, rel(r1, r2));
   }

   // Lever 3: ComputeVolumeRHS.
   wave.SetDerivMode(DerivMode::OnTheFly);
   Vector v1(N);
   v1 = 0.0;
   wave.ComputeVolumeRHS_ForTest(Q, v1);
   wave.SetDerivMode(DerivMode::Cached);
   Vector v2(N);
   v2 = 0.0;
   wave.ComputeVolumeRHS_ForTest(Q, v2);
   wave.SetDerivMode(DerivMode::OnTheFly);
   local_max = std::max(local_max, rel(v1, v2));

   double global_max = 0.0;
   MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm);

   const int failed = (global_max <= 1e-12) ? 0 : 1;
   if (rank == 0)
   {
      std::cout << (failed ? "  FAILED" : "  PASSED")
                << ": ApplySpatialDerivative + ComputeVolumeRHS cached==onthefly "
                << "(max rel " << global_max << ", tol 1e-12)\n";
   }

   MPI_Finalize();
   return failed;
#else
   (void)argc; (void)argv;
   std::cout << "Parallel test requires MFEM_USE_MPI.\n";
   return 0;
#endif
}
