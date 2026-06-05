// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 1 of the fault-dealiasing plan
// (document/fault_dealiasing_dev/seisol_overintegration_resample_speckle_2026-06-02.md):
// UNIT TEST for fault-flux OVER-INTEGRATION (WaveOperator::SetFaultOverint).
//
// Runs on a TINY 2-tet fault fixture (one fault triangle) — never a production
// mesh.  Verifies the plumbing and the two contracts that matter:
//
//   [T1] default (k=0): GetNbfPerFace() == IntRules.Get(TRIANGLE, 2*order)
//        count, and FaultFaceQuadDegree() == 2*order.   (byte-exact-when-off)
//   [T2] SetFaultOverint(k): FaultFaceQuadDegree() == 2*(order+k) and
//        GetNbfPerFace() grows to the 2*(order+k) rule count.  (knob is LIVE)
//   [T3] the per-QP FaultBasis is resized to the grown count
//        (GetBasis(0).qp_data.size() == GetNbfPerFace()).
//   [T4] SetFaultOverint(0) is idempotent: nbf returns to the baseline.
//   [T5] BYTE-EXACT-WHEN-OFF: one AdvanceADER step on an operator that called
//        SetFaultOverint(0) is bit-identical to one that never called it.
//        (The rebuild-at-k=0 path must reproduce the ctor's fault quadrature
//        exactly — degree 2*order, same FaultBasis ordering.)
//   [T6] knob is LIVE on the flux: one AdvanceADER step with k=1 DIFFERS from
//        the k=0 step (over-integration actually changes the radiated flux).
//
// Usage:  make seas_test_fault_overint ; ./seas_test_fault_overint

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
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define CHECK(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace {

constexpr real_t kL = 1000.0;
constexpr real_t kDt = 1.0e-4;
constexpr int    kOrder = 1;

// Two tets sharing one fault triangle at y = 0 (attr 3); exterior attr 1
// (natural).  Lifted from test_rupture_multistep_serial_vs_parallel.cpp.
Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {kL, 0.0, 0.0}, {0.0, 0.0, kL},
      {0.0,  kL, 0.0},   // V3 on +y side
      {0.0, -kL, 0.0},   // V4 on -y side
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);
   mesh.AddTet(0, 1, 2, 3, 1);
   mesh.FinalizeTopology();
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
         if (std::abs(cy) < 1e-8) { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Build fault DOFData sized to wave.GetNbfPerFace() (which reflects the current
// over-integration factor), uniform TPV102 nucleation drive at every QP.
int SetupFault(WaveOperator<Mesh> &wave, Mesh &mesh,
               std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
               std::vector<Vector> &fault_coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const int nqp = wave.GetNbfPerFace();
   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), wave.FaultFaceQuadDegree());
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   const int n_fault = int_faces.Size() * nqp;
   if (n_fault > 0)
   {
      InitializeFaultDOFs(dof_data, n_fault, fault_coords);
      ZeroDOFDataPreStressTotal(dof_data, n_fault);
      for (int i = 0; i < n_fault; i++) { dof_data[i].tau2_nuc = TPV102Params::nuc_dtau; }
   }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);
   return n_fault;
}

// One ADER-2 step; returns the post-step bulk Q (for byte-exact comparison).
Vector OneStep(int overint_k)
{
   Mesh mesh = BuildTwoTetFaultMesh();
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;

   WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                           TPV102Params::mu, TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bulk_bg);
   if (overint_k > 0) { wave.SetFaultOverint(overint_k); }

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   SetupFault(wave, mesh, dof_data, ff, fault_coords);

   const int ndof_total = wave.GetFESpace().GetNDofs();
   Vector Q(wave.Height());
   InitializeStateTotal(Q, ndof_total, TPV102Params::sigma_n, TPV102Params::tau_ini);
   Vector Q_new(wave.Height());
   wave.AdvanceADER(Q, kDt, /*ader_order=*/2, Q_new);
   return Q_new;
}

real_t MaxAbsDiff(const Vector &a, const Vector &b)
{
   real_t d = 0.0;
   for (int i = 0; i < a.Size(); i++) { d = std::max(d, std::abs(a(i) - b(i))); }
   return d;
}

} // anonymous namespace

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
#else
   (void)argc; (void)argv;
#endif

   std::cout << "\n=== Phase 1: fault-flux over-integration unit test ===\n";

   const int tri_deg_base = 2 * kOrder;        // baseline fault rule degree
   const int tri_deg_k1   = 2 * (kOrder + 1);  // k=1 over-integrated degree
   const int nbf_base = IntRules.Get(Geometry::TRIANGLE, tri_deg_base).GetNPoints();
   const int nbf_k1   = IntRules.Get(Geometry::TRIANGLE, tri_deg_k1).GetNPoints();

   Mesh mesh = BuildTwoTetFaultMesh();
   BoundaryConfig bc; bc.natural_attrs = {1}; bc.fault_attr = 3;
   WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                           TPV102Params::mu, TPV102Params::rho, bc);

   // [T1] default (k=0) is byte-exact: degree 2*order, baseline QP count.
   CHECK(wave.GetFaultOverint() == 0, "default over-integration factor is 0");
   CHECK(wave.FaultFaceQuadDegree() == tri_deg_base,
         "default FaultFaceQuadDegree() == 2*order");
   CHECK(wave.GetNbfPerFace() == nbf_base,
         "default nbf_per_face == IntRules(TRIANGLE, 2*order) point count");

   // [T2] SetFaultOverint(1): degree and QP count grow.
   wave.SetFaultOverint(1);
   CHECK(wave.GetFaultOverint() == 1, "SetFaultOverint(1) records k=1");
   CHECK(wave.FaultFaceQuadDegree() == tri_deg_k1,
         "FaultFaceQuadDegree() == 2*(order+1) after SetFaultOverint(1)");
   CHECK(wave.GetNbfPerFace() == nbf_k1 && nbf_k1 > nbf_base,
         "nbf_per_face grows to the 2*(order+1) rule count (over-integrated)");

   // [T3] the per-QP FaultBasis tracks the grown count.
   const FaultBasis *fb = wave.GetFaultBasis();
   CHECK(fb != nullptr && fb->NumFaces() > 0 &&
         static_cast<int>(fb->GetBasis(0).qp_data.size()) == nbf_k1,
         "FaultBasis qp_data resized to the over-integrated QP count");

   // [T4] SetFaultOverint(0) restores the baseline (idempotent off).
   wave.SetFaultOverint(0);
   CHECK(wave.FaultFaceQuadDegree() == tri_deg_base &&
         wave.GetNbfPerFace() == nbf_base,
         "SetFaultOverint(0) restores the baseline degree + QP count");

   // [T5] BYTE-EXACT-WHEN-OFF: rebuild-at-k=0 == never-rebuilt, to the last bit.
   {
      Vector q_ctor = OneStep(/*overint_k=*/0);   // never calls SetFaultOverint
      // Operator that explicitly rebuilds at k=0 then steps:
      Mesh m2 = BuildTwoTetFaultMesh();
      BoundaryConfig bc2; bc2.natural_attrs = {1}; bc2.fault_attr = 3;
      real_t bg[NUM_STATE] = {0};
      bg[SYY] = TPV102Params::sigma_n; bg[SXY] = -TPV102Params::tau_ini;
      WaveOperator<Mesh> w2(m2, kOrder, TPV102Params::lambda, TPV102Params::mu,
                            TPV102Params::rho, bc2);
      w2.SetAbsorbingBackground(bg);
      w2.SetFaultOverint(0);   // explicit rebuild at the baseline degree
      FaultFaceFlux ff2(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
      std::vector<DOFData> dd2; std::vector<Vector> fc2;
      SetupFault(w2, m2, dd2, ff2, fc2);
      const int ndof2 = w2.GetFESpace().GetNDofs();
      Vector Q2(w2.Height());
      InitializeStateTotal(Q2, ndof2, TPV102Params::sigma_n, TPV102Params::tau_ini);
      Vector Q2n(w2.Height());
      w2.AdvanceADER(Q2, kDt, 2, Q2n);
      const real_t d = MaxAbsDiff(q_ctor, Q2n);
      CHECK(d == 0.0,
            "AdvanceADER with SetFaultOverint(0) is BIT-IDENTICAL to never "
            "calling it (max|dQ| == 0)");
   }

   // [T6] knob is LIVE on the flux: k=1 step differs from the k=0 step.
   {
      Vector q0 = OneStep(0);
      Vector q1 = OneStep(1);
      CHECK(MaxAbsDiff(q0, q1) > 0.0,
            "AdvanceADER with --fault-overint 1 changes the radiated flux "
            "(over-integration is actually applied)");
   }

   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n"
             << "========================================\n";

   int rc = (num_failed == 0) ? 0 : 1;
#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return rc;
}
