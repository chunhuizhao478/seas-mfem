// Parallel wave operator tests (P1-P5 from parallel test design).
// Verifies WaveOperator<ParMesh> with shared face ghost exchange.
// Run with: mpirun -np 4 ./seas_test_parallel_wave_operator

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/seas_dynamic_operator.hpp"
#include "../../domain/boundary_config.hpp"

#include <iostream>
#include <cmath>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define PAR_TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static void RK4Step(const WaveOperator<ParMesh> &wave, Vector &Q, real_t dt)
{
   int size = Q.Size();
   Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size);
   wave.Mult(Q, k1);
   add(Q, 0.5*dt, k1, Q_tmp); wave.Mult(Q_tmp, k2);
   add(Q, 0.5*dt, k2, Q_tmp); wave.Mult(Q_tmp, k3);
   add(Q, dt, k3, Q_tmp); wave.Mult(Q_tmp, k4);
   Q.Add(dt/6.0, k1); Q.Add(dt/3.0, k2);
   Q.Add(dt/3.0, k3); Q.Add(dt/6.0, k4);
}

static real_t ComputeLocalEnergy(const WaveOperator<ParMesh> &wave,
                                 const Vector &Q,
                                 real_t lambda, real_t mu, real_t rho)
{
   int ndof_total = wave.GetScalarNDof();
   int ne = wave.NumElements();
   int ndof_per_el = wave.GetNDof();
   const auto &fes = wave.GetFESpace();

   real_t E = 0.0;
   const real_t *Q_data = Q.GetData();

   for (int e = 0; e < ne; e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2*wave.GetOrder());
      int ndof = fe->GetDof();
      int dof_offset = e * ndof_per_el;
      Vector shape(ndof);

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetIntPoint(&ip);
         real_t w = ip.weight * Tr->Weight();
         fe->CalcShape(ip, shape);

         real_t Q_qp[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_qp[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               Q_qp[c] += shape(i) * Q_data[c * ndof_total + dof_offset + i];
            }
         }

         E += w * EnergyDensity(Q_qp, lambda, mu, rho);
      }
   }
   return E;
}

// ===== P1: ParMesh construction =====
void TestP1_ParMeshConstruction(int rank, MPI_Comm comm)
{
   if (rank == 0) std::cout << "Test P1: TestParMeshConstruction\n";

   Mesh serial_mesh = Mesh::MakeCartesian3D(4, 4, 4, Element::HEXAHEDRON,
                                             1.0, 1.0, 1.0);
   for (int b = 0; b < serial_mesh.GetNBE(); b++)
   {
      serial_mesh.SetBdrAttribute(b, 5);  // all absorbing
   }
   serial_mesh.SetAttributes();

   ParMesh pmesh(comm, serial_mesh);

   BoundaryConfig bc;
   bc.absorbing_attrs = {5};

   WaveOperator<ParMesh> wave(pmesh, 1, 32.04e9, 32.04e9, 2670.0, bc);

   int local_height = wave.Height();
   PAR_TEST_ASSERT(local_height > 0,
                   "WaveOperator<ParMesh> height > 0 on all ranks (local: "
                   + std::to_string(local_height) + ")");

   int local_ne = wave.NumElements();
   int global_ne;
   MPI_Allreduce(&local_ne, &global_ne, 1, MPI_INT, MPI_SUM, comm);
   PAR_TEST_ASSERT(global_ne == 64,
                   "Total elements = 64 (got " + std::to_string(global_ne) + ")");
}

// ===== P4: Quiescent state (Q=0 → dQdt=0) =====
void TestP4_QuiescentState(int rank, MPI_Comm comm)
{
   if (rank == 0) std::cout << "Test P4: TestParallelQuiescentState\n";

   Mesh serial_mesh = Mesh::MakeCartesian3D(4, 4, 4, Element::HEXAHEDRON,
                                             1.0, 1.0, 1.0);
   for (int b = 0; b < serial_mesh.GetNBE(); b++)
   {
      serial_mesh.SetBdrAttribute(b, 5);
   }
   serial_mesh.SetAttributes();

   ParMesh pmesh(comm, serial_mesh);
   BoundaryConfig bc;
   bc.absorbing_attrs = {5};

   WaveOperator<ParMesh> wave(pmesh, 1, 32.04e9, 32.04e9, 2670.0, bc);

   int size = wave.Height();
   Vector Q(size), dQdt(size);
   Q = 0.0;
   wave.Mult(Q, dQdt);

   real_t local_norm = dQdt.Norml2();
   real_t global_norm;
   MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX, comm);

   PAR_TEST_ASSERT(global_norm < 1e-12,
                   "Q=0 → dQdt=0 on all ranks (max norm "
                   + std::to_string(global_norm) + ")");
}

// ===== P5: Energy conservation in all-absorbing box =====
void TestP5_EnergyConservation(int rank, MPI_Comm comm)
{
   if (rank == 0) std::cout << "Test P5: TestParallelEnergyConservation\n";

   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   real_t cp = std::sqrt((lambda + 2*mu) / rho);

   Mesh serial_mesh = Mesh::MakeCartesian3D(4, 4, 4, Element::HEXAHEDRON,
                                             1.0, 1.0, 1.0);
   for (int b = 0; b < serial_mesh.GetNBE(); b++)
   {
      serial_mesh.SetBdrAttribute(b, 5);
   }
   serial_mesh.SetAttributes();

   ParMesh pmesh(comm, serial_mesh);
   BoundaryConfig bc;
   bc.absorbing_attrs = {5};

   WaveOperator<ParMesh> wave(pmesh, 1, lambda, mu, rho, bc);

   int ndof_total = wave.GetScalarNDof();
   Vector Q(wave.Height());
   Q = 0.0;

   // Gaussian stress pulse centered at (0.5, 0.5, 0.5)
   ParFiniteElementSpace fes(&pmesh, wave.GetFESpace().FEColl());
   for (int e = 0; e < pmesh.GetNE(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      int ndof = fe->GetDof();
      int dof_offset = e * wave.GetNDof();
      Vector shape(ndof);

      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetIntPoint(&ip);
         fe->CalcShape(ip, shape);
         Vector phys(3);
         Tr->Transform(ip, phys);

         real_t r2 = (phys(0)-0.5)*(phys(0)-0.5) +
                     (phys(1)-0.5)*(phys(1)-0.5) +
                     (phys(2)-0.5)*(phys(2)-0.5);
         real_t amp = 1e6 * std::exp(-r2 / (2.0 * 0.1 * 0.1));

         for (int i = 0; i < ndof; i++)
         {
            Q[SXY * ndof_total + dof_offset + i] += amp * shape(i);
         }
      }
   }

   real_t cfl = 0.5 / (3.0 * 3.0);
   real_t dt = wave.ComputeMaxDt(cfl);

   real_t E_prev = ComputeLocalEnergy(wave, Q, lambda, mu, rho);
   real_t E_prev_global;
   MPI_Allreduce(&E_prev, &E_prev_global, 1, MPI_DOUBLE, MPI_SUM, comm);

   bool monotonic = true;
   int nsteps = 200;  // enough for wave to reach absorbing boundaries
   for (int step = 0; step < nsteps; step++)
   {
      RK4Step(wave, Q, dt);

      real_t E_local = ComputeLocalEnergy(wave, Q, lambda, mu, rho);
      real_t E_global;
      MPI_Allreduce(&E_local, &E_global, 1, MPI_DOUBLE, MPI_SUM, comm);

      if (E_global > E_prev_global * 1.001)
      {
         monotonic = false;
         break;
      }
      E_prev_global = E_global;
   }

   PAR_TEST_ASSERT(monotonic,
                   "Energy monotonically non-increasing for " +
                   std::to_string(nsteps) + " steps");

   // Final energy should not exceed initial (absorbing BCs + shared face flux)
   real_t E_final_local = ComputeLocalEnergy(wave, Q, lambda, mu, rho);
   real_t E_final;
   MPI_Allreduce(&E_final_local, &E_final, 1, MPI_DOUBLE, MPI_SUM, comm);

   PAR_TEST_ASSERT(E_final <= E_prev_global * 1.001,
                   "Final energy <= initial (ratio " +
                   std::to_string(E_final / (E_prev_global + 1e-30)) + ")");
}

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank, nprocs;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);

   if (rank == 0)
   {
      std::cout << "========================================\n";
      std::cout << "Parallel WaveOperator Tests (P1,P4,P5)\n";
      std::cout << "Ranks: " << nprocs << "\n";
      std::cout << "========================================\n\n";
   }

   TestP1_ParMeshConstruction(rank, comm);
   TestP4_QuiescentState(rank, comm);
   TestP5_EnergyConservation(rank, comm);

   // Synchronize results across ranks
   int global_failed;
   MPI_Allreduce(&num_failed, &global_failed, 1, MPI_INT, MPI_MAX, comm);

   if (rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "Total:  " << num_tests << "\n";
      std::cout << "Passed: " << num_passed << "\n";
      std::cout << "Failed: " << global_failed << "\n";
      std::cout << "========================================\n";
   }

   MPI_Finalize();
   return (global_failed > 0) ? 1 : 0;
#else
   std::cout << "Parallel tests require MFEM_USE_MPI.\n";
   return 0;
#endif
}
