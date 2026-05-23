// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// SAFS curvilinear-fault blow-up: H1-vs-H3 DISCRIMINATOR.
//
// ============================================================================
// Motivation
// ============================================================================
// The SAFS dynamic-rupture run blows up on shared (cross-rank) fault faces,
// invariant to D_c (Dc=10 over-resolved still blows up), cfl, mu_s,
// mixed-flux, and dip-prestress.  Two hypotheses remain:
//
//   H1 (MPI / cross-rank inconsistency): the per-QP fault frame
//       (FaultBasis::ComputeOrientedFrame) canonicalizes orientation via
//       sign(n_raw . ref_normal) with ref_normal = (0,-1,0) HARDCODED
//       (wave_operator.inl:349), and the +/- side via (elem1_c - face_c) .
//       ref_normal (wave_operator.inl:459-465).  BOTH degrade as the local
//       fault normal approaches PERPENDICULAR to ref_normal — the margin
//       collapses to O(|n . ref_normal|) -> FP-noise-determined -> the two
//       ranks owning a shared face can pick OPPOSITE strike directions, and
//       the fixed-sign nucleation increment (+F.dtau into tau2_nuc) then
//       forces the two sides oppositely.  PLANAR TPV faults have
//       |n . ref_normal| = 1 and never hit this; a CURVILINEAR fault does.
//
//   H3 (intrinsic, NOT MPI): the undamped fault-normal/opening channel is
//       unstable regardless of partitioning; it would blow up in serial too.
//
// CLEAN SEPARATOR:  H1 needs shared faces (parallel only).  H3 does not.
// So: run the SAME tilted-fault problem SERIAL and np=2, sweeping the fault
// tilt angle theta through the degenerate regime, and compare.
//
//   - serial == parallel to 1e-10 at EVERY theta (incl. theta -> 90 deg,
//     n perp ref_normal)  =>  MPI is consistent at all fault angles =>
//     H1 ELIMINATED, including for varying fault trace angle.  The blow-up
//     is then intrinsic (H3).
//   - serial != parallel at some theta  =>  H1 CONFIRMED at that angle; the
//     failing component (SXX=normal vs SXY/SXZ=tangential vs VX/VY/VZ) and
//     the angle localize the orientation/+side bug.
//
// The background stress tensor is ROTATED with the fault (R_z(theta)) so the
// fault sees the SAME physical (sigma_n, tau) at every angle: the ONLY thing
// that changes across the sweep is the conditioning of the ref_normal-based
// canonicalization.  This isolates the mechanism.
//
// ============================================================================
// Usage:  mpirun -np 2 ./seas_test_rupture_tilted_fault_serial_vs_parallel
//         Requires MFEM built with MPI.

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

constexpr real_t kL = 1000.0;
constexpr real_t kDt = 1.0e-4;
constexpr int    kNSteps = 30;
constexpr int    kOrder = 1;

const char *StateName(int c)
{
   static const char *kNames[NUM_STATE] =
   {"SXX", "SYY", "SZZ", "SXY", "SYZ", "SXZ", "VX", "VY", "VZ"};
   return (c >= 0 && c < NUM_STATE) ? kNames[c] : "?";
}

// 2-tet fault fixture, fault plane rotated about z by theta.  At theta=0 the
// fault normal is (0,1,0) (|n . ref_normal| = 1, planar TPV regime); at
// theta=90 deg the fault normal is (-1,0,0) (n PERP ref_normal=(0,-1,0), the
// degenerate regime).  Rotation about z keeps n in the xy-plane so it is
// never parallel to up=(0,0,1) (no up-degeneracy is introduced — only the
// ref_normal-perpendicularity is swept).
Mesh BuildTiltedTwoTetFaultMesh(real_t theta_deg)
{
   const real_t th = theta_deg * M_PI / 180.0;
   const real_t c = std::cos(th), s = std::sin(th);
   const real_t base[5][3] = {
      {0.0, 0.0, 0.0}, {kL, 0.0, 0.0}, {0.0, 0.0, kL},
      {0.0,  kL, 0.0}, {0.0, -kL, 0.0},
   };
   Mesh mesh(3, 5, 2, 0);
   for (int v = 0; v < 5; v++)
   {
      const real_t x = base[v][0], y = base[v][1], z = base[v][2];
      real_t rv[3] = { c*x - s*y, s*x + c*y, z };   // R_z(theta)
      mesh.AddVertex(rv);
   }
   mesh.AddTet(0, 1, 2, 4, 1);   // tet 0 -> rank 0
   mesh.AddTet(0, 1, 2, 3, 1);   // tet 1 -> rank 1
   mesh.FinalizeTopology();
   // Tag the single interior face (the shared fault) as attr 3; every
   // boundary face as natural attr 1.  Rotation-robust (no cy test).
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      const int attr = (ftr && ftr->Elem2No >= 0) ? 3 : 1;
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], attr);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Background stress tensor rotated with the fault: sigma' = R_z sigma R_z^T,
// sigma_0 = {SYY=sigma_n, SXY=-tau_ini}.  Keeps the fault loading physically
// identical at every theta.
void RotatedBackground(real_t theta_deg, real_t bg[NUM_STATE])
{
   const real_t th = theta_deg * M_PI / 180.0;
   const real_t c = std::cos(th), s = std::sin(th);
   const real_t b = TPV102Params::sigma_n;   // sigma_yy0
   const real_t h = -TPV102Params::tau_ini;   // sigma_xy0
   for (int i = 0; i < NUM_STATE; i++) { bg[i] = 0.0; }
   bg[SXX] = b*s*s - 2.0*h*c*s;
   bg[SYY] = b*c*c + 2.0*h*c*s;
   bg[SXY] = (0.0 - b)*c*s + h*(c*c - s*s);
}

void FillStateTotal(Vector &Q, int ndof, const real_t bg[NUM_STATE])
{
   Q = 0.0;
   for (int c = 0; c < NUM_STATE; c++)
   {
      for (int d = 0; d < ndof; d++) { Q(c*ndof + d) = bg[c]; }
   }
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
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
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

   const int n_fault = (int_faces.Size() + shr_faces.Size()) * nqp;
   if (n_fault > 0)
   {
      InitializeFaultDOFs(dof_data, n_fault, fault_coords);
      ZeroDOFDataPreStressTotal(dof_data, n_fault);
      // Uniform full nucleation in the local STRIKE channel at every QP —
      // identical scalar on serial and parallel.  The ONLY way serial and
      // parallel can disagree is if they assign a different physical strike
      // direction to the shared QP (the H1 orientation flip).
      for (int i = 0; i < n_fault; i++)
      {
         dof_data[i].tau2_nuc = TPV102Params::nuc_dtau;
      }
   }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);
   return n_fault;
}

struct RunResult
{
   std::vector<std::vector<real_t>> Q_per_step;   // [step][comp] at elem0/dof0
   std::vector<real_t> slip_per_step;             // |slip_rate| at fault QP 0
   real_t max_V = 0.0;
};

// Serial run (rank 0 only).  Records bulk Q at local element 0 / dof 0.
RunResult RunSerial(real_t theta_deg, const BoundaryConfig &bc,
                    const real_t bg[NUM_STATE])
{
   RunResult r;
   r.Q_per_step.assign(kNSteps, std::vector<real_t>(NUM_STATE, 0.0));
   r.slip_per_step.assign(kNSteps, 0.0);

   Mesh mesh = BuildTiltedTwoTetFaultMesh(theta_deg);
   WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                           TPV102Params::mu, TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bg);

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   const int n_fault = SetupFault(wave, mesh, kOrder, dof_data, ff,
                                  fault_coords);
   MFEM_VERIFY(n_fault > 0, "serial: no fault QPs");

   const auto &fes = wave.GetFESpace();
   const int ndof = fes.GetNDofs();
   Array<int> edofs; fes.GetElementDofs(0, edofs);
   const int dof0 = edofs[0];

   Vector Q(wave.Height());
   FillStateTotal(Q, ndof, bg);
   Vector Q_new(wave.Height());
   for (int step = 0; step < kNSteps; step++)
   {
      wave.AdvanceADER(Q, kDt, /*ader_order=*/2, Q_new);
      Q.Swap(Q_new);
      for (int c = 0; c < NUM_STATE; c++)
      { r.Q_per_step[step][c] = Q(c*ndof + dof0); }
      r.slip_per_step[step] = dof_data[0].slip_rate;
      r.max_V = std::max(r.max_V, std::abs(dof_data[0].slip_rate));
   }
   return r;
}

} // anonymous

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (rank == 0)
   {
      std::cout << "\n=== SAFS curvilinear-fault H1-vs-H3 discriminator ===\n"
                << "  tilted-fault serial-vs-parallel sweep, np = " << nprocs
                << ", N = " << kNSteps << " ADER-2 steps\n"
                << "  ref_normal = (0,-1,0) HARDCODED; theta=90 deg => fault "
                << "normal PERP ref_normal (degenerate)\n";
   }
   if (nprocs != 2)
   {
      if (rank == 0) { std::cout << "  SKIPPED: requires np == 2\n"; }
      MPI_Finalize();
      return 77;
   }

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};

   const real_t thetas[] = {0.0, 30.0, 45.0, 60.0, 75.0, 85.0, 89.0, 90.0};
   const int n_theta = static_cast<int>(sizeof(thetas)/sizeof(thetas[0]));

   if (rank == 0)
   {
      std::cout << "\n" << std::setw(8) << "theta"
                << std::setw(16) << "|n.ref_n|"
                << std::setw(16) << "worst_rel"
                << std::setw(7) << "comp"
                << std::setw(15) << "max|V|_ser"
                << std::setw(15) << "max|V|_par" << "\n";
   }

   real_t worst_rel_global = 0.0;
   real_t worst_theta = -1.0;

   for (int ti = 0; ti < n_theta; ti++)
   {
      const real_t theta = thetas[ti];
      real_t bg[NUM_STATE];
      RotatedBackground(theta, bg);

      // |n . ref_normal| = |cos(theta)| for this rotation.
      const real_t ndot = std::abs(std::cos(theta * M_PI / 180.0));

      // ---- serial (rank 0) ----
      RunResult ser;
      if (rank == 0) { ser = RunSerial(theta, bc, bg); }

      // broadcast serial Q + slip + max_V to all ranks.
      std::vector<real_t> ser_flat(kNSteps * NUM_STATE, 0.0);
      std::vector<real_t> ser_slip(kNSteps, 0.0);
      real_t ser_maxV = 0.0;
      if (rank == 0)
      {
         for (int s = 0; s < kNSteps; s++)
         {
            for (int c = 0; c < NUM_STATE; c++)
            { ser_flat[s*NUM_STATE + c] = ser.Q_per_step[s][c]; }
            ser_slip[s] = ser.slip_per_step[s];
         }
         ser_maxV = ser.max_V;
      }
      MPI_Bcast(ser_flat.data(), kNSteps*NUM_STATE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(ser_slip.data(), kNSteps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&ser_maxV, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      // ---- parallel (np=2, explicit partition tet0->rank0, tet1->rank1) ----
      Mesh serial_mesh = BuildTiltedTwoTetFaultMesh(theta);
      std::vector<int> part(serial_mesh.GetNE(), 0);
      for (int e = 0; e < serial_mesh.GetNE(); e++) { part[e] = e; }  // e->rank e
      ParMesh pmesh(MPI_COMM_WORLD, serial_mesh, part.data());

      WaveOperator<ParMesh> wave_p(pmesh, kOrder, TPV102Params::lambda,
                                   TPV102Params::mu, TPV102Params::rho, bc);
      wave_p.SetAbsorbingBackground(bg);
      FaultFaceFlux ff_p(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
      std::vector<DOFData> dof_p;
      std::vector<Vector> fc_p;
      const int n_fault_p = SetupFault(wave_p, pmesh, kOrder, dof_p, ff_p, fc_p);

      const auto &fes_p = wave_p.GetFESpace();
      const int ndof_p = fes_p.GetNDofs();
      // rank 0 owns global tet 0 as its local element 0.
      int dof0_p = -1;
      if (rank == 0 && pmesh.GetNE() > 0)
      {
         Array<int> ed; fes_p.GetElementDofs(0, ed);
         dof0_p = ed.Size() > 0 ? ed[0] : -1;
      }

      Vector Qp(wave_p.Height());
      FillStateTotal(Qp, ndof_p, bg);
      Vector Qp_new(wave_p.Height());

      std::vector<std::vector<real_t>> par(kNSteps,
                                           std::vector<real_t>(NUM_STATE, 0.0));
      std::vector<real_t> par_slip(kNSteps, 0.0);
      real_t par_maxV = 0.0;
      for (int step = 0; step < kNSteps; step++)
      {
         wave_p.AdvanceADER(Qp, kDt, /*ader_order=*/2, Qp_new);
         Qp.Swap(Qp_new);
         if (rank == 0 && dof0_p >= 0)
         {
            for (int c = 0; c < NUM_STATE; c++)
            { par[step][c] = Qp(c*ndof_p + dof0_p); }
            if (n_fault_p > 0)
            {
               par_slip[step] = dof_p[0].slip_rate;
               par_maxV = std::max(par_maxV, std::abs(dof_p[0].slip_rate));
            }
         }
      }

      // ---- compare on rank 0 ----
      if (rank == 0)
      {
         real_t worst_rel = 0.0; int worst_comp = -1;
         for (int step = 0; step < kNSteps; step++)
         {
            for (int c = 0; c < NUM_STATE; c++)
            {
               const real_t a = ser_flat[step*NUM_STATE + c];
               const real_t b = par[step][c];
               const real_t denom = std::max({std::abs(a), std::abs(b),
                                               real_t(1.0)});
               const real_t rel = std::abs(a - b) / denom;
               if (rel > worst_rel) { worst_rel = rel; worst_comp = c; }
            }
            const real_t sd =
               std::abs(ser_slip[step] - par_slip[step])
               / std::max({std::abs(ser_slip[step]), std::abs(par_slip[step]),
                           real_t(1.0)});
            if (sd > worst_rel) { worst_rel = sd; worst_comp = -2; }  // -2: slip
         }
         if (worst_rel > worst_rel_global)
         { worst_rel_global = worst_rel; worst_theta = theta; }

         std::cout << std::setw(8) << std::fixed << std::setprecision(1) << theta
                   << std::setw(16) << std::scientific << std::setprecision(3) << ndot
                   << std::setw(16) << worst_rel
                   << std::setw(7) << (worst_comp == -2 ? "slip"
                                       : StateName(worst_comp))
                   << std::setw(15) << ser_maxV
                   << std::setw(15) << par_maxV << "\n";
      }
   }

   if (rank == 0)
   {
      std::cout << "\n  worst serial-vs-parallel rel diff = "
                << std::scientific << std::setprecision(6) << worst_rel_global
                << "  at theta = " << std::fixed << std::setprecision(1)
                << worst_theta << " deg\n\n";
      TEST_LE(worst_rel_global, 1.0e-10,
              "tilted-fault serial==parallel at ALL fault trace angles "
              "(PASS => H1 eliminated, blow-up is intrinsic/H3; "
              "FAIL => H1 confirmed at the reported theta)");
      std::cout << "\n========================================\n"
                << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests << " tests\n"
                << "========================================\n";
   }

   int exit_code = (num_failed == 0) ? 0 : 1;
   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
   return exit_code;
}
