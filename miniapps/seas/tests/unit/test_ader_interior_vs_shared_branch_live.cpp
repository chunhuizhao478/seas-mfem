// ADER live bulk-fault coupling cross-check on the same physical face.
//
// This is the ADER analogue of test_interior_vs_shared_branch_live.cpp:
// rank 0 runs a serial 2-tet fixture where the fault uses the INTERIOR
// branch, then all ranks run the same geometry as a 2-rank ParMesh where
// the fault uses the SHARED branch. We compare the one-step ADER update
// slope (Q_new - Q) / dt on the minus-side element.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

namespace
{

constexpr real_t kFixtureSize = 1000.0;
constexpr real_t kVelAmp      = 2.0e-5;
constexpr real_t kStressAmp   = 5.0e4;
constexpr real_t kDt          = 1.0e-4;

Mesh BuildTwoTetFaultMesh(real_t L = kFixtureSize)
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {L, 0.0, 0.0}, {0.0, 0.0, L},
      {0.0,  L, 0.0}, {0.0, -L, 0.0},
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);
   mesh.AddTet(0, 1, 2, 3, 1);
   mesh.FinalizeTopology();

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0.0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-8)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }

   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

template <typename MeshT>
int SetupFault(WaveOperator<MeshT> &wave, MeshT &mesh, int order,
               std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
               std::vector<Vector> &fault_coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();

   int nqp_per_face = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault face FTR null");
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2 * order).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      if (nqp_per_face == 0 && shr_faces.Size() > 0)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[0]);
         MFEM_VERIFY(ftr, "shared fault face FTR null");
         nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2 * order).GetNPoints();
      }
   }
#endif

   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      for (int i = 0; i < shr_faces.Size(); i++)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[i]);
         const IntegrationRule &ir =
            IntRules.Get(ftr->GetGeometryType(), 2 * order);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            Vector phys(3);
            ftr->Face->Transform(ip, phys);
            fault_coords.push_back(phys);
         }
      }
   }
#endif

   const int nfault = (int_faces.Size() + shr_faces.Size()) * nqp_per_face;
   InitializeFaultDOFs(dof_data, nfault, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, nfault);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);
   return nfault;
}

template <typename MeshT>
void AddSideConstant(const MeshT &mesh, const FiniteElementSpace &fes,
                     Vector &Q, int comp, real_t val_plus, real_t val_minus)
{
   const int ndof_total = fes.GetNDofs();
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      Array<int> ev;
      mesh.GetElementVertices(e, ev);
      real_t cy = 0.0;
      for (int v = 0; v < ev.Size(); v++) { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();

      const real_t val = (cy > 0.0) ? val_plus : val_minus;
      Array<int> edofs;
      fes.GetElementDofs(e, edofs);
      for (int j = 0; j < edofs.Size(); j++)
      {
         Q(comp * ndof_total + edofs[j]) += val;
      }
   }
}

template <typename MeshT>
void BuildExcitedState(const MeshT &mesh, const WaveOperator<MeshT> &wave, Vector &Q)
{
   const FiniteElementSpace &fes = wave.GetFESpace();
   InitializeStateTotal(Q, fes.GetNDofs(),
                        TPV102Params::sigma_n, TPV102Params::tau_ini);
   AddSideConstant(mesh, fes, Q, VX,  +kVelAmp,        -0.5 * kVelAmp);
   AddSideConstant(mesh, fes, Q, VY,  -0.7 * kVelAmp,  +0.4 * kVelAmp);
   AddSideConstant(mesh, fes, Q, VZ,  +0.3 * kVelAmp,  -0.2 * kVelAmp);
   AddSideConstant(mesh, fes, Q, SXY, +kStressAmp,     -0.6 * kStressAmp);
   AddSideConstant(mesh, fes, Q, SXZ, -0.5 * kStressAmp, +0.4 * kStressAmp);
}

template <typename MeshT>
void ExtractMinusSideSlope(MeshT &mesh, WaveOperator<MeshT> &wave, Vector &slope)
{
   const auto &fes = wave.GetFESpace();
   Vector Q;
   BuildExcitedState(mesh, wave, Q);
   Vector Q_new(Q.Size());
   wave.AdvanceADER(Q, kDt, /*order=*/2, Q_new);

   const int ndof_total = fes.GetNDofs();
   int e_minus = -1;
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      Array<int> ev;
      mesh.GetElementVertices(e, ev);
      real_t cy = 0.0;
      for (int v = 0; v < ev.Size(); v++) { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      if (cy < 0.0) { e_minus = e; break; }
   }
   slope.SetSize(NUM_STATE);
   slope = 0.0;
   if (e_minus < 0)
   {
      return;
   }

   Array<int> edofs;
   fes.GetElementDofs(e_minus, edofs);
   const int dof0 = edofs[0];
   for (int c = 0; c < NUM_STATE; c++)
   {
      slope(c) = (Q_new(c * ndof_total + dof0) - Q(c * ndof_total + dof0)) / kDt;
   }
}

} // namespace

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (nprocs != 2)
   {
      if (rank == 0)
      {
         std::cout << "SKIPPED: test_ader_interior_vs_shared_branch_live requires np=2\n";
      }
      MPI_Finalize();
      return 77;
   }

   const int order = 1;
   real_t bulk_bg[NUM_STATE] = {0.0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;

   Vector slope_serial(NUM_STATE);
   if (rank == 0)
   {
      Mesh mesh = BuildTwoTetFaultMesh();
      BoundaryConfig bc;
      bc.natural_attrs = {1};
      bc.fault_attr = 3;
      bc.absorbing_attrs = {};
      WaveOperator<Mesh> wave(mesh, order,
                              TPV102Params::lambda,
                              TPV102Params::mu,
                              TPV102Params::rho, bc);
      wave.SetAbsorbingBackground(bulk_bg);

      FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
      std::vector<DOFData> dof_data;
      std::vector<Vector> fault_coords;
      SetupFault(wave, mesh, order, dof_data, ff, fault_coords);
      ExtractMinusSideSlope(mesh, wave, slope_serial);
   }

   MPI_Bcast(slope_serial.GetData(), NUM_STATE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   Mesh serial_mesh = BuildTwoTetFaultMesh();
   std::vector<int> partitioning(serial_mesh.GetNE());
   for (int e = 0; e < serial_mesh.GetNE(); e++)
   {
      Array<int> ev;
      serial_mesh.GetElementVertices(e, ev);
      real_t cy = 0.0;
      for (int v = 0; v < ev.Size(); v++) { cy += serial_mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      partitioning[e] = (cy < 0.0) ? 0 : 1;
   }

   ParMesh pmesh(MPI_COMM_WORLD, serial_mesh, partitioning.data());
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};
   WaveOperator<ParMesh> wave(pmesh, order,
                              TPV102Params::lambda,
                              TPV102Params::mu,
                              TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bulk_bg);

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   SetupFault(wave, pmesh, order, dof_data, ff, fault_coords);

   Vector slope_parallel(NUM_STATE);
   ExtractMinusSideSlope(pmesh, wave, slope_parallel);

   int exit_code = 0;
   if (rank == 0)
   {
      std::cout << "\n=== ADER live interior-vs-shared branch cross-check ===\n";
      int n_fail = 0;
      real_t worst_rel = 0.0;
      for (int c = 0; c < NUM_STATE; c++)
      {
         const real_t scale =
            std::max({std::abs(slope_serial(c)),
                      std::abs(slope_parallel(c)),
                      real_t(1.0)});
         const real_t diff = std::abs(slope_serial(c) - slope_parallel(c));
         const real_t rel = diff / scale;
         worst_rel = std::max(worst_rel, rel);
         std::cout << "  " << c
                   << "  serial=" << std::scientific << std::setprecision(6)
                   << slope_serial(c)
                   << "  parallel=" << slope_parallel(c)
                   << "  rel=" << rel << "\n";
         if (rel > 1.0e-10) { n_fail++; }
      }
      std::cout << "  worst relative mismatch = " << worst_rel << "\n";
      if (n_fail == 0)
      {
         std::cout << "  PASSED: ADER interior branch matches shared branch\n";
      }
      else
      {
         std::cout << "  FAILED: ADER interior/shared branch divergence\n";
         exit_code = 1;
      }
      std::cout << "========================================\n";
   }

   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
   return exit_code;
}
