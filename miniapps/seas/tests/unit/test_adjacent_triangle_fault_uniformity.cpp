// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug: ADJACENT-TRIANGLE UNIFORMITY test.
//
// ============================================================================
// Motivation
// ============================================================================
// The 400-rank Frontera run 7672897 shows per-cell pepper on the fault
// surface inside the ruptured area.  WriteFaultSurfaceVTU emits one
// scalar per triangle (per-cell QP average), so adjacent triangles
// genuinely disagree in physics.
//
// Existing pepper-bug tests and the new
// `test_rupture_multistep_serial_vs_parallel` all use a 2-tet fixture
// with exactly ONE fault triangle — so triangle-to-triangle
// inconsistency cannot be probed there.
//
// This test fills that gap:
//   1. Build an 8-fault-triangle fixture (2×2×2 Cartesian hex → 48 tets,
//      fault plane at y = L/2 splits 4 hexes each giving 2 triangles).
//   2. Set DOFData.tau2_nuc = nuc_dtau UNIFORMLY on every fault QP.
//   3. Run N=20 ADER-2 steps at a CFL-friendly dt.
//   4. Under uniform material + uniform forcing + planar fault + symmetric
//      initial Q, EVERY fault QP must see physically IDENTICAL state.
//      Assert max across-QP spread of slip_rate, tau1_corr, tau2_corr,
//      sigma_n_corr stays < 1e-10 relative on every step.
//
// Interpretation:
//   - PASS: all 24 fault QPs (8 triangles × 3 QPs) stay bit-identical →
//     per-triangle dispatch is consistent on a 2x2-plane-of-fault fixture.
//   - FAIL at triangle i, component C, step K: triangles i and j disagree.
//     The test dumps which triangles diverge and by how much → directly
//     points to the per-triangle code path responsible.
//
// Runs in serial (np=1) first; np=4 parallel version follows the same
// logic but partitions the fault across ranks.
//
// Usage:
//   serial:   ./seas_test_adjacent_triangle_fault_uniformity
//   parallel: mpirun -np 4 ./seas_test_adjacent_triangle_fault_uniformity

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

namespace {

constexpr real_t kL = 1000.0;       // 1 km cube; dt safe for p=1 ADER-2
constexpr real_t kDt = 5.0e-5;
constexpr int    kNSteps = 20;
constexpr int    kOrder = 1;

// 2x2x2 Cartesian hex, then convert to tets.  MFEM's MakeCartesian3D
// splits each hex into 6 tets.  Fault plane y = L/2 cuts 2x2 = 4 hexes,
// giving 4 hex-faces = 8 fault triangles (each hex face = 2 tet triangles).
Mesh BuildCartesianFaultMesh()
{
   // Cartesian 2 × 2 × 2 hex mesh spanning [0, L]^3; y = L/2 cuts the
   // middle.  Use Hexahedron elements and then FinalizeTet-style split.
   // MFEM's Mesh::MakeCartesian3D with Element::TETRAHEDRON already
   // produces a tet mesh with conforming fault-plane triangles.
   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::TETRAHEDRON,
                                     kL, kL, kL, /*sfc_ordering=*/false);
   mesh.FinalizeTopology();
   mesh.Finalize();

   // Set fault attribute 3 on triangles at y = L/2 (interior faces only);
   // free-surface attribute 1 on all exterior triangles.
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy - 0.5 * kL) < 1e-8)
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

   int nqp = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault FTR null");
      nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      if (nqp == 0 && shr_faces.Size() > 0)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[0]);
         MFEM_VERIFY(ftr, "shared fault FTR null");
         nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
      }
   }
#endif

   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      for (int i = 0; i < shr_faces.Size(); i++)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[i]);
         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                   2*order);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            Vector phys(3); ftr->Face->Transform(ip, phys);
            fault_coords.push_back(phys);
         }
      }
   }
#endif

   const int n_fault =
      (int_faces.Size() + shr_faces.Size()) * nqp;
   if (n_fault > 0)
   {
      InitializeFaultDOFs(dof_data, n_fault, fault_coords);
      ZeroDOFDataPreStressTotal(dof_data, n_fault);
      // Persistent nucleation uniform across every fault QP.
      for (int i = 0; i < n_fault; i++)
      {
         dof_data[i].tau2_nuc = TPV102Params::nuc_dtau;
#ifdef SEAS_DIAG_FAULT_FLUX
         dof_data[i].diag_print = true;
#endif
      }
   }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);
   return n_fault;
}

struct UniformityStats
{
   int    n_qp_local = 0;
   real_t min_slip = std::numeric_limits<real_t>::max();
   real_t max_slip = std::numeric_limits<real_t>::lowest();
   real_t min_tau1 = std::numeric_limits<real_t>::max();
   real_t max_tau1 = std::numeric_limits<real_t>::lowest();
   real_t min_tau2 = std::numeric_limits<real_t>::max();
   real_t max_tau2 = std::numeric_limits<real_t>::lowest();
   real_t min_sn   = std::numeric_limits<real_t>::max();
   real_t max_sn   = std::numeric_limits<real_t>::lowest();
};

UniformityStats CollectStats(const std::vector<DOFData> &dof_data)
{
   UniformityStats s;
   s.n_qp_local = static_cast<int>(dof_data.size());
   for (const DOFData &d : dof_data)
   {
      s.min_slip = std::min(s.min_slip, d.slip_rate);
      s.max_slip = std::max(s.max_slip, d.slip_rate);
      s.min_tau1 = std::min(s.min_tau1, d.tau1_corr);
      s.max_tau1 = std::max(s.max_tau1, d.tau1_corr);
      s.min_tau2 = std::min(s.min_tau2, d.tau2_corr);
      s.max_tau2 = std::max(s.max_tau2, d.tau2_corr);
      s.min_sn   = std::min(s.min_sn,   d.sigma_n_corr);
      s.max_sn   = std::max(s.max_sn,   d.sigma_n_corr);
   }
   return s;
}

#ifdef MFEM_USE_MPI
void GlobalReduce(UniformityStats &s, MPI_Comm comm)
{
   real_t local_min[4] = {s.min_slip, s.min_tau1, s.min_tau2, s.min_sn};
   real_t local_max[4] = {s.max_slip, s.max_tau1, s.max_tau2, s.max_sn};
   real_t global_min[4], global_max[4];
   if (s.n_qp_local == 0)
   {
      // Set sentinels that never win min/max reductions.
      local_min[0] = local_min[1] = local_min[2] = local_min[3] =
         std::numeric_limits<real_t>::max();
      local_max[0] = local_max[1] = local_max[2] = local_max[3] =
         std::numeric_limits<real_t>::lowest();
   }
   MPI_Allreduce(local_min, global_min, 4, MPI_DOUBLE, MPI_MIN, comm);
   MPI_Allreduce(local_max, global_max, 4, MPI_DOUBLE, MPI_MAX, comm);
   int total_qp = 0;
   MPI_Allreduce(&s.n_qp_local, &total_qp, 1, MPI_INT, MPI_SUM, comm);
   s.n_qp_local = total_qp;
   s.min_slip = global_min[0];  s.max_slip = global_max[0];
   s.min_tau1 = global_min[1];  s.max_tau1 = global_max[1];
   s.min_tau2 = global_min[2];  s.max_tau2 = global_max[2];
   s.min_sn   = global_min[3];  s.max_sn   = global_max[3];
}
#endif

} // anonymous

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#else
   const int rank = 0, nprocs = 1;
#endif
   const bool is_parallel = (nprocs > 1);

   if (rank == 0)
   {
      std::cout << "\n=== TPV102 pepper-bug: ADJACENT-TRIANGLE "
                << "UNIFORMITY test ===\n"
                << "  MPI ranks: " << nprocs
                << "  (is_parallel=" << (is_parallel?"yes":"no") << ")\n"
                << "  fixture: 2x2x2 Cartesian hex → tets, fault at y=L/2\n"
                << "  L = " << kL << " m, dt = " << std::scientific
                << std::setprecision(3) << kDt << " s, N = " << kNSteps
                << " ADER-2 steps\n"
                << "  UNIFORM tau2_nuc = " << TPV102Params::nuc_dtau
                << " Pa on every fault QP\n"
                << "  expected: every fault QP → identical state (uniform "
                << "material + forcing + planar fault)\n";
   }

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;

   Mesh serial_mesh = BuildCartesianFaultMesh();
   if (rank == 0)
   {
      std::cout << "  serial mesh: NE = " << serial_mesh.GetNE()
                << ", NumFaces = " << serial_mesh.GetNumFaces() << "\n";
   }

   // Choose serial or parallel construction.
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   int n_fault_local = 0;
   int ndof_total = 0;
   int size = 0;
   Vector Q, Q_new;

#ifdef MFEM_USE_MPI
   if (is_parallel)
   {
      ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);
      if (rank == 0)
      {
         std::cout << "  parallel mesh: pmesh NE (rank 0) = "
                   << pmesh.GetNE() << "\n";
      }
      WaveOperator<ParMesh> wave(pmesh, kOrder, TPV102Params::lambda,
                                  TPV102Params::mu, TPV102Params::rho, bc);
      wave.SetAbsorbingBackground(bulk_bg);
      n_fault_local = SetupFault(wave, pmesh, kOrder, dof_data, ff,
                                  fault_coords);
      int total_n_fault = 0;
      MPI_Allreduce(&n_fault_local, &total_n_fault, 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
      if (rank == 0)
      {
         std::cout << "  total fault QPs (global): " << total_n_fault
                   << "  (each rank has " << n_fault_local << " on rank 0)\n";
      }

      ndof_total = wave.GetFESpace().GetNDofs();
      size = wave.Height();
      Q.SetSize(size); Q_new.SetSize(size);
      InitializeStateTotal(Q, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);

      // Run-and-check loop.
      real_t worst_slip_spread = 0.0, worst_tau1_spread = 0.0;
      real_t worst_tau2_spread = 0.0, worst_sn_spread = 0.0;
      int    worst_step = -1;
      for (int step = 0; step < kNSteps; step++)
      {
         wave.AdvanceADER(Q, kDt, /*ader_order=*/2, Q_new);
         Q.Swap(Q_new);
         UniformityStats s = CollectStats(dof_data);
         GlobalReduce(s, MPI_COMM_WORLD);
         const real_t slip_mid = 0.5 * (s.max_slip + s.min_slip);
         const real_t tau1_mid = 0.5 * (s.max_tau1 + s.min_tau1);
         const real_t tau2_mid = 0.5 * (s.max_tau2 + s.min_tau2);
         const real_t sn_mid   = 0.5 * (s.max_sn   + s.min_sn  );
         const real_t slip_spread = (s.max_slip - s.min_slip) /
            std::max(std::abs(slip_mid), real_t(1.0));
         const real_t tau1_spread = (s.max_tau1 - s.min_tau1) /
            std::max(std::abs(tau1_mid), real_t(1.0));
         const real_t tau2_spread = (s.max_tau2 - s.min_tau2) /
            std::max(std::abs(tau2_mid), real_t(1.0));
         const real_t sn_spread   = (s.max_sn   - s.min_sn) /
            std::max(std::abs(sn_mid),   real_t(1.0));
         if (slip_spread > worst_slip_spread)
         { worst_slip_spread = slip_spread; worst_step = step; }
         worst_tau1_spread = std::max(worst_tau1_spread, tau1_spread);
         worst_tau2_spread = std::max(worst_tau2_spread, tau2_spread);
         worst_sn_spread   = std::max(worst_sn_spread,   sn_spread);
         if (rank == 0)
         {
            std::cout << "  step " << std::setw(3) << step
                      << "  slip_rate [" << std::scientific
                      << std::setprecision(3) << s.min_slip << " ... "
                      << s.max_slip << "]  spread = " << slip_spread
                      << "  |  tau2_corr [" << s.min_tau2 << " ... "
                      << s.max_tau2 << "]  spread = " << tau2_spread
                      << "\n";
         }
      }
      if (rank == 0)
      {
         std::cout << "\n  worst spreads over " << kNSteps << " steps:\n"
                   << "    slip_rate : " << std::scientific
                   << std::setprecision(6) << worst_slip_spread
                   << " (step " << worst_step << ")\n"
                   << "    tau1_corr : " << worst_tau1_spread << "\n"
                   << "    tau2_corr : " << worst_tau2_spread << "\n"
                   << "    sigma_n_corr: " << worst_sn_spread << "\n";
         TEST_LE(worst_slip_spread, 1.0e-10,
                 "parallel: slip_rate uniform across all fault QPs");
         TEST_LE(worst_tau2_spread, 1.0e-10,
                 "parallel: tau2_corr uniform across all fault QPs");
         TEST_LE(worst_tau1_spread, 1.0e-10,
                 "parallel: tau1_corr uniform across all fault QPs");
         TEST_LE(worst_sn_spread, 1.0e-10,
                 "parallel: sigma_n_corr uniform across all fault QPs");
      }
   }
   else
#endif
   {
      WaveOperator<Mesh> wave(serial_mesh, kOrder, TPV102Params::lambda,
                               TPV102Params::mu, TPV102Params::rho, bc);
      wave.SetAbsorbingBackground(bulk_bg);
      n_fault_local = SetupFault(wave, serial_mesh, kOrder, dof_data, ff,
                                  fault_coords);
      if (rank == 0)
      {
         std::cout << "  serial fault QPs: " << n_fault_local << "\n";

         // DIAGNOSTIC: dump per-QP face geometry + which element DOFs
         // shape1/shape2 sample.  If two triangles share the same physical
         // Q_bulk but sample DIFFERENT element DOFs with different shape
         // weights, that's the per-triangle asymmetry source.
         const Array<int> &int_faces = wave.GetFaultInteriorFaces();
         const auto &fes = wave.GetFESpace();
         std::cout << "\n  per-QP geometry dump:\n";
         for (int i = 0; i < int_faces.Size(); i++)
         {
            const int face_idx = int_faces[i];
            auto *ftr = serial_mesh.GetInteriorFaceTransformations(face_idx);
            const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                      2*kOrder);
            const FiniteElement *fe1 = fes.GetFE(ftr->Elem1No);
            const FiniteElement *fe2 = fes.GetFE(ftr->Elem2No);
            const int nd1 = fe1->GetDof();
            const int nd2 = fe2->GetDof();
            for (int q = 0; q < ir.GetNPoints(); q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               ftr->SetAllIntPoints(&ip);
               Vector phys(3); ftr->Face->Transform(ip, phys);
               Vector n_raw(3); CalcOrtho(ftr->Face->Jacobian(), n_raw);
               const real_t nl = n_raw.Norml2();

               // Reference-coord of ip on elem1 and elem2.
               IntegrationPoint ip1, ip2;
               ftr->Loc1.Transform(ip, ip1);
               ftr->Loc2.Transform(ip, ip2);
               Vector s1(nd1), s2(nd2);
               fe1->CalcShape(ip1, s1);
               fe2->CalcShape(ip2, s2);
               const int qp_idx = i * ir.GetNPoints() + q;
               std::cout << "    qp=" << std::setw(3) << qp_idx
                         << "  face=" << std::setw(3) << face_idx
                         << "  q=" << q
                         << "  nl=" << std::scientific
                         << std::setprecision(5) << nl
                         << "  e1=" << std::noshowpos << ftr->Elem1No
                         << "  e2=" << ftr->Elem2No
                         << "\n      ip1_ref=(" << std::fixed
                         << std::setprecision(4)
                         << ip1.x << "," << ip1.y << "," << ip1.z
                         << ")  ip2_ref=(" << ip2.x << "," << ip2.y
                         << "," << ip2.z << ")\n"
                         << "      shape1=[";
               for (int k = 0; k < nd1; k++)
               { std::cout << std::setprecision(4) << s1(k)
                           << (k < nd1-1 ? "," : ""); }
               std::cout << "]  shape2=[";
               for (int k = 0; k < nd2; k++)
               { std::cout << std::setprecision(4) << s2(k)
                           << (k < nd2-1 ? "," : ""); }
               std::cout << "]\n";
            }
         }
      }
      ndof_total = wave.GetFESpace().GetNDofs();
      size = wave.Height();
      Q.SetSize(size); Q_new.SetSize(size);
      InitializeStateTotal(Q, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);

      real_t worst_slip_spread = 0.0, worst_tau1_spread = 0.0;
      real_t worst_tau2_spread = 0.0, worst_sn_spread = 0.0;
      int    worst_step = -1;
      auto dump_per_qp = [&](int step) {
         std::cout << "\n  PER-QP STATE at step " << step << ":\n";
         for (int qp = 0; qp < n_fault_local; qp++)
         {
            const DOFData &d = dof_data[qp];
            std::cout << "    qp=" << std::setw(3) << qp
                      << "  V1=" << std::scientific << std::setprecision(6)
                      << d.V1 << "  V2=" << d.V2
                      << "  |V|=" << d.slip_rate
                      << "  tau1_c=" << std::showpos
                      << std::setprecision(3) << d.tau1_corr
                      << "  tau2_c=" << d.tau2_corr
                      << "  sn_c=" << std::noshowpos << d.sigma_n_corr
                      << "\n";
         }
      };

      // Dump Q at EVERY DOF of every fault-adjacent tet.  If bulk Q
      // within a single tet becomes non-uniform, we've localized the
      // per-DOF asymmetry.  Also dump Q[SXY] and Q[SXZ] (which map to
      // dip/strike in global Cartesian — fault-local rotation aside,
      // non-uniformity in either = per-DOF flux deposition asymmetry).
      const auto &fes = wave.GetFESpace();
      auto dump_bulk_at_fault_tets = [&](int step) {
         std::cout << "\n  BULK Q at fault-adjacent tets at step "
                   << step << ":\n";
         // Track tets that touch the fault plane (y=L/2).
         std::vector<int> fault_tets;
         const Array<int> &int_faces = wave.GetFaultInteriorFaces();
         for (int i = 0; i < int_faces.Size(); i++)
         {
            auto *ftr = serial_mesh.GetInteriorFaceTransformations(
               int_faces[i]);
            fault_tets.push_back(ftr->Elem1No);
            fault_tets.push_back(ftr->Elem2No);
         }
         std::sort(fault_tets.begin(), fault_tets.end());
         fault_tets.erase(std::unique(fault_tets.begin(),
                                       fault_tets.end()),
                          fault_tets.end());
         for (int e : fault_tets)
         {
            Array<int> edofs; fes.GetElementDofs(e, edofs);
            // Compute centroid.
            Array<int> ev; serial_mesh.GetElementVertices(e, ev);
            real_t cx = 0, cy = 0, cz = 0;
            for (int v = 0; v < ev.Size(); v++)
            {
               cx += serial_mesh.GetVertex(ev[v])[0];
               cy += serial_mesh.GetVertex(ev[v])[1];
               cz += serial_mesh.GetVertex(ev[v])[2];
            }
            cx /= ev.Size(); cy /= ev.Size(); cz /= ev.Size();
            std::cout << "    tet=" << std::setw(3) << e
                      << "  cy=" << std::fixed << std::setprecision(1)
                      << cy << "  SXY=[";
            for (int k = 0; k < edofs.Size(); k++)
            {
               std::cout << std::scientific << std::setprecision(4)
                         << std::showpos
                         << Q(SXY * ndof_total + edofs[k])
                         << (k < edofs.Size() - 1 ? "," : "");
            }
            std::cout << "]\n             SXZ=[";
            for (int k = 0; k < edofs.Size(); k++)
            {
               std::cout << std::scientific << std::setprecision(4)
                         << std::showpos
                         << Q(SXZ * ndof_total + edofs[k])
                         << (k < edofs.Size() - 1 ? "," : "");
            }
            std::cout << std::noshowpos << "]\n";
         }
      };
      for (int step = 0; step < kNSteps; step++)
      {
         wave.AdvanceADER(Q, kDt, /*ader_order=*/2, Q_new);
         Q.Swap(Q_new);
         UniformityStats s = CollectStats(dof_data);
         const real_t slip_mid = 0.5 * (s.max_slip + s.min_slip);
         const real_t tau1_mid = 0.5 * (s.max_tau1 + s.min_tau1);
         const real_t tau2_mid = 0.5 * (s.max_tau2 + s.min_tau2);
         const real_t sn_mid   = 0.5 * (s.max_sn   + s.min_sn  );
         const real_t slip_spread = (s.max_slip - s.min_slip) /
            std::max(std::abs(slip_mid), real_t(1.0));
         const real_t tau1_spread = (s.max_tau1 - s.min_tau1) /
            std::max(std::abs(tau1_mid), real_t(1.0));
         const real_t tau2_spread = (s.max_tau2 - s.min_tau2) /
            std::max(std::abs(tau2_mid), real_t(1.0));
         const real_t sn_spread   = (s.max_sn   - s.min_sn  ) /
            std::max(std::abs(sn_mid),   real_t(1.0));
         if (slip_spread > worst_slip_spread)
         { worst_slip_spread = slip_spread; worst_step = step; }
         worst_tau1_spread = std::max(worst_tau1_spread, tau1_spread);
         worst_tau2_spread = std::max(worst_tau2_spread, tau2_spread);
         worst_sn_spread   = std::max(worst_sn_spread,   sn_spread);
         std::cout << "  step " << std::setw(3) << step
                   << "  slip_rate [" << std::scientific
                   << std::setprecision(3) << s.min_slip << " ... "
                   << s.max_slip << "]  spread = " << slip_spread
                   << "  |  tau2_corr [" << s.min_tau2 << " ... "
                   << s.max_tau2 << "]  spread = " << tau2_spread
                   << "  |  tau1_corr [" << s.min_tau1 << " ... "
                   << s.max_tau1 << "]\n";
         if (step == 0 || step == 1 || step == 2 || step == 5 || step == 19)
         {
            dump_per_qp(step);
            dump_bulk_at_fault_tets(step);
         }
      }
      std::cout << "\n  worst spreads over " << kNSteps << " steps:\n"
                << "    slip_rate : " << std::scientific
                << std::setprecision(6) << worst_slip_spread
                << " (step " << worst_step << ")\n"
                << "    tau1_corr : " << worst_tau1_spread << "\n"
                << "    tau2_corr : " << worst_tau2_spread << "\n"
                << "    sigma_n_corr: " << worst_sn_spread << "\n";
      TEST_LE(worst_slip_spread, 1.0e-10,
              "serial: slip_rate uniform across all fault QPs");
      TEST_LE(worst_tau2_spread, 1.0e-10,
              "serial: tau2_corr uniform across all fault QPs");
      TEST_LE(worst_tau1_spread, 1.0e-10,
              "serial: tau1_corr uniform across all fault QPs");
      TEST_LE(worst_sn_spread, 1.0e-10,
              "serial: sigma_n_corr uniform across all fault QPs");
   }

   if (rank == 0)
   {
      std::cout << "\n========================================\n"
                << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests << " tests\n"
                << "========================================\n";
   }
   int exit_code = (num_failed == 0) ? 0 : 1;
#ifdef MFEM_USE_MPI
   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
#endif
   return exit_code;
}
