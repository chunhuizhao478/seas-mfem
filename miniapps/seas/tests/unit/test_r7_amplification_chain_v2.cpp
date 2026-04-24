// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Round-7 review (2026-04-23): strengthened amplification-chain + per-QP
// stage-resolved diagnostic.
//
// Round-6 §H's T1 probe used centroid-bucket orbit pairing and covered
// only 2 of 8 fault-triangle pairs (R-001).  It also skipped T2
// (trial-traction stage between Q reconstruction and friction output,
// R-002).  This test fixes both:
//
//   T1 (fixed coverage): pair ALL 8 fault-triangle Elem1/Elem2 pairs.
//       Measure orbit spread across all 8 upper-side fault-adjacent
//       tets (and separately lower-side).  Q_self at each QP is
//       computed via shape·Q.  Symmetry-aware: expect Q_self_plus
//       and Q_self_minus to be related by y-reflection symmetry
//       (SXY, SYZ, VY flip sign; others preserved).
//
//   T2 (NEW): trial-traction spread.  Call FaultFaceFlux::ComputeTrialTraction
//       with each fault triangle's (Q_self_plus, Q_self_minus)
//       in the canonical frame.  Report orbit spread of
//       (sigma_n_trial, tau1_trial, tau2_trial) across the 8
//       fault triangles.
//
//   T3: DOFData tau1_corr / tau2_corr / sigma_n_corr spread (global,
//       matches pepper guard's v9.4.0 §11 metric).
//
//   T3b (sub-stage): DOFData V1, V2 spread — reveals whether the
//       amplification is upstream of the V-decomposition (V1 clean,
//       tau1_corr dirty → writeback) or at the V-decomposition itself.
//
// Decision tree (reviewer round-7):
//   T1 ≈ 0 with full coverage → DG layer genuinely clean.
//   T1 > 0 → DG layer NOT clean; reopen at DG.
//   T2 ≈ 0 and T3 > 0 → asymmetry in Brent / V-decomposition / writeback
//                        (stages 4-6 of FaultFaceFlux::Evaluate).
//   T2 > 0 and T3 ≈ T2 → asymmetry in ComputeTrialTraction or rotation
//                         (stages 1-3).

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/friction_solver.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/d4_tet_mesh.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../fault/fault_basis.hpp"
#include <cstdlib>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

namespace
{

constexpr real_t kL      = 1000.0;
constexpr real_t kDt     = 5.0e-5;
constexpr int    kNSteps = 20;
constexpr int    kOrder  = 1;
constexpr int    kAderOrder = 2;

const char *kCompName[NUM_STATE] = {"SXX","SYY","SZZ","SXY","SYZ","SXZ",
                                    "VX","VY","VZ"};

Mesh BuildKuhnFaultMesh()
{
   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::TETRAHEDRON,
                                     kL, kL, kL, /*sfc_ordering=*/false);
   mesh.FinalizeTopology();
   mesh.Finalize();
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
         { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Round-8 Step 0 (R-005): honor SEAS_TEST_FIXTURE env var for
// cross-fixture verification.
Mesh BuildFaultMeshByEnv()
{
   const char *fx = std::getenv("SEAS_TEST_FIXTURE");
   if (fx && std::string(fx) == "d4")
   {
      std::cout << "  [SEAS_TEST_FIXTURE=d4] using BuildD4Mesh(true).\n";
      return BuildD4Mesh(/*add_fault=*/true, kL);
   }
   std::cout << "  [fixture=kuhn] (default)\n";
   return BuildKuhnFaultMesh();
}

// Per-triangle per-QP raw and rotated Q_self.  "Plus" / "Minus" refer
// to the (Elem1, Elem2) convention.
struct TriangleProbe
{
   int face_idx;
   int elem_plus;
   int elem_minus;
   std::vector<real_t> Q_self_plus_rot;    // [NUM_STATE * nqp]
   std::vector<real_t> Q_self_minus_rot;   // [NUM_STATE * nqp]
   // Trial tractions (sigma_n_trial, tau1_trial, tau2_trial) per QP
   std::vector<real_t> sigma_n_trial_qp;
   std::vector<real_t> tau1_trial_qp;
   std::vector<real_t> tau2_trial_qp;
};

// Identify fault triangles as (face_idx, Elem1, Elem2) triples.
std::vector<std::tuple<int,int,int>> IdentifyFaultTriangles(Mesh &mesh,
                                                             int fault_attr = 3)
{
   std::vector<std::tuple<int,int,int>> out;
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      if (mesh.GetBdrAttribute(b) != fault_attr) { continue; }
      Array<int> bv; mesh.GetBdrElementVertices(b, bv);
      for (int f = 0; f < mesh.GetNumFaces(); f++)
      {
         Array<int> fv; mesh.GetFaceVertices(f, fv);
         if (fv.Size() != bv.Size()) { continue; }
         std::vector<int> a(bv.begin(), bv.end()), c(fv.begin(), fv.end());
         std::sort(a.begin(), a.end()); std::sort(c.begin(), c.end());
         if (a != c) { continue; }
         auto *ftr = mesh.GetFaceElementTransformations(f);
         if (!ftr || ftr->Elem1No < 0 || ftr->Elem2No < 0) { continue; }
         // Canonical: elem_plus = side with higher cy (upper-y).
         Array<int> ev1; mesh.GetElementVertices(ftr->Elem1No, ev1);
         real_t cy1 = 0;
         for (int v = 0; v < ev1.Size(); v++)
         { cy1 += mesh.GetVertex(ev1[v])[1]; }
         cy1 /= ev1.Size();
         const int e_plus  = (cy1 > 0.5 * kL) ? ftr->Elem1No : ftr->Elem2No;
         const int e_minus = (cy1 > 0.5 * kL) ? ftr->Elem2No : ftr->Elem1No;
         out.emplace_back(f, e_plus, e_minus);
         break;
      }
   }
   return out;
}

// Max - min across a vector of real values.
real_t SpreadMaxMin(const std::vector<real_t> &v)
{
   if (v.empty()) { return 0.0; }
   real_t mn = v[0], mx = v[0];
   for (real_t x : v) { mn = std::min(mn, x); mx = std::max(mx, x); }
   return mx - mn;
}

// Compute Q_self at a fault QP on one side (element e), rotating into
// the canonical frame via the FaultBasis (same convention as
// wave_operator.inl's fault dispatch).
void ComputeQSelfCan(const Vector &Q, const FiniteElementSpace &fes,
                      FaceElementTransformations *ftr, int e,
                      const IntegrationPoint &ip_face,
                      const FaultBasisQPData &qpd,
                      real_t Q_self_can[NUM_STATE])
{
   const int ndof_total = fes.GetNDofs();
   const FiniteElement *fe = fes.GetFE(e);
   const int ndof = fe->GetDof();
   const int dof_offset = e * ndof;
   IntegrationPoint ip_elem;
   if (ftr->Elem1No == e) { ftr->Loc1.Transform(ip_face, ip_elem); }
   else                    { ftr->Loc2.Transform(ip_face, ip_elem); }
   Vector shape(ndof);
   fe->CalcShape(ip_elem, shape);

   real_t Q_self[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      real_t s = 0.0;
      for (int i = 0; i < ndof; i++)
      {
         s += shape(i) * Q(c * ndof_total + dof_offset + i);
      }
      Q_self[c] = s;
   }

   // Canonical frame from FaultBasis data (matches wave_operator.inl).
   real_t can_n[3], can_t1[3], can_t2[3];
   for (int d = 0; d < 3; d++)
   {
      can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]    : qpd.normal[d];
      can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d]  : qpd.tangent1[d];
      can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d]  : qpd.tangent2[d];
   }
   DenseMatrix Tinv_can(NUM_STATE);
   GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv_can);
   for (int c = 0; c < NUM_STATE; c++)
   {
      real_t v = 0.0;
      for (int k = 0; k < NUM_STATE; k++)
      { v += Tinv_can(c, k) * Q_self[k]; }
      Q_self_can[c] = v;
   }
}

struct ChainRecord
{
   int step;
   // T1: max component-wise spread of Q_self_plus_rot and Q_self_minus_rot
   //     across the 8 fault triangles, at each QP.
   real_t T1_plus_spread_max;   // worst over (comp, qp)
   real_t T1_minus_spread_max;
   int    T1_plus_worst_c, T1_plus_worst_qp;
   int    T1_minus_worst_c, T1_minus_worst_qp;

   // T2: spread of (sigma_n_trial, tau1_trial, tau2_trial) across triangles, per QP.
   real_t T2_sigma_n_spread;
   real_t T2_tau1_spread;
   real_t T2_tau2_spread;

   // T3: DOFData tau1_corr / tau2_corr / sigma_n_corr spread (global).
   real_t T3_tau1_corr_spread;
   real_t T3_tau2_corr_spread;
   real_t T3_sigma_n_corr_spread;

   // T3b: DOFData V1 / V2 / slip_rate spread (sub-stage).
   real_t T3b_V1_spread;
   real_t T3b_V2_spread;
   real_t T3b_slip_rate_spread;
};

real_t GlobalSpread(const std::vector<DOFData> &dof,
                     std::function<real_t(const DOFData&)> accessor)
{
   if (dof.empty()) { return 0.0; }
   real_t mn = accessor(dof[0]), mx = mn;
   for (const DOFData &d : dof)
   { const real_t v = accessor(d); mn = std::min(mn,v); mx = std::max(mx,v); }
   return mx - mn;
}

} // anonymous

int main()
{
   std::cout << "\n=== Round-7: strengthened amplification chain ===\n";
   std::cout << std::scientific << std::setprecision(6);

   Mesh mesh = BuildFaultMeshByEnv();
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;
   real_t bulk_bg[NUM_STATE] = {0};

   WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                            TPV102Params::mu, TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bulk_bg);
   const auto &fes = wave.GetFESpace();

   // Fault DOFData setup (same as pepper guard).
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   int nqp = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      nqp = IntRules.Get(ftr->GetGeometryType(), 2*kOrder).GetNPoints();
   }
   std::vector<Vector> fault_coords;
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*kOrder);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   const int n_fault = int_faces.Size() * nqp;
   std::vector<DOFData> dof_data;
   InitializeFaultDOFs(dof_data, n_fault, fault_coords);
   for (int i = 0; i < n_fault; i++)
   { dof_data[i].tau2_nuc = TPV102Params::nuc_dtau; }

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);

   // Shared-face triangle pairing (R-001 fix).
   std::vector<std::tuple<int,int,int>> triangles = IdentifyFaultTriangles(mesh);
   std::cout << "  fault triangles (shared-face pairs): " << triangles.size()
             << "  (expected 8 on Kuhn 2x2x2)\n";
   std::cout << "  QPs per face: " << nqp
             << "   total fault QPs: " << n_fault << "\n\n";

   // Build FaultBasis-like lookup.  The wave operator already built
   // `fault_basis_`; re-use via GetFaultBasis().
   const FaultBasis *fb = wave.GetFaultBasis();
   if (!fb)
   {
      std::cout << "FAIL: FaultBasis not built by WaveOperator\n";
      return 1;
   }

   Vector Q(wave.Height()), Q_new(wave.Height());
   Q = 0.0;

   std::vector<ChainRecord> records;
   std::cout << "  step | T1+ max  | T1- max  | T2 sigma_n | T2 tau1  | T2 tau2  |"
             << " T3 tau1  | T3 tau2  | T3 sigma_n| T3b V1  | T3b V2  | T3b slip_rate\n";
   std::cout << "  -----|----------|----------|------------|----------|----------|"
             << "----------|----------|-----------|---------|---------|---------------\n";

   for (int step = 0; step < kNSteps; step++)
   {
      wave.AdvanceADER(Q, kDt, kAderOrder, Q_new);
      Q.Swap(Q_new);

      // Build per-triangle probes for this step.
      std::vector<TriangleProbe> probes;
      probes.reserve(triangles.size());

      for (std::size_t ti = 0; ti < triangles.size(); ti++)
      {
         int f, ep, em;
         std::tie(f, ep, em) = triangles[ti];
         TriangleProbe p;
         p.face_idx = f; p.elem_plus = ep; p.elem_minus = em;

         // Find the FaultBasis index for this face (matches
         // wave_operator.inl's LookupInteriorFaultBasisIndex logic).
         int fb_idx = -1;
         for (int i = 0; i < int_faces.Size(); i++)
         {
            if (int_faces[i] == f) { fb_idx = i; break; }
         }
         if (fb_idx < 0) { continue; }

         auto *ftr = mesh.GetInteriorFaceTransformations(f);
         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*kOrder);

         p.Q_self_plus_rot.assign(NUM_STATE * nqp, 0.0);
         p.Q_self_minus_rot.assign(NUM_STATE * nqp, 0.0);
         p.sigma_n_trial_qp.assign(nqp, 0.0);
         p.tau1_trial_qp.assign(nqp, 0.0);
         p.tau2_trial_qp.assign(nqp, 0.0);

         const FaultBasisData &bd = fb->GetBasis(fb_idx);

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);

            if (q >= static_cast<int>(bd.qp_data.size())) { continue; }
            const FaultBasisQPData &qpd = bd.qp_data[q];

            // Plus side (upper y by our convention).  In the bulk
            // Q path, wave_operator.inl's interior-fault branch
            // identifies which side is plus based on
            // `qpd.sign_flipped` (elem1_on_plus = !qpd.sign_flipped).
            // We use the centroid-based plus/minus assignment from
            // IdentifyFaultTriangles; this may differ from MFEM's
            // Elem1/Elem2 assignment but is consistent across the
            // per-triangle sweep and captures the orbit-spread
            // signal regardless of which label is "+" vs "-".
            real_t Q_self_plus_can[NUM_STATE], Q_self_minus_can[NUM_STATE];
            ComputeQSelfCan(Q, fes, ftr, ep, ip, qpd, Q_self_plus_can);
            ComputeQSelfCan(Q, fes, ftr, em, ip, qpd, Q_self_minus_can);

            for (int c = 0; c < NUM_STATE; c++)
            {
               p.Q_self_plus_rot [c * nqp + q] = Q_self_plus_can[c];
               p.Q_self_minus_rot[c * nqp + q] = Q_self_minus_can[c];
            }

            // T2: trial traction via public ComputeTrialTraction.
            // Need a DOFData with eta_p/eta_s filled in.  The test's
            // dof_data vector has these (set by InitializeFaultDOFs).
            // Find any dof_data entry that matches this (face, QP) —
            // the face index ordering matches int_faces[0..N-1]
            // concatenation, each with nqp entries.  Entry index = fb_idx * nqp + q.
            const int dof_idx = fb_idx * nqp + q;
            if (dof_idx < static_cast<int>(dof_data.size()))
            {
               real_t sigma_n_trial, tau1_trial, tau2_trial;
               FaultFaceFlux::ComputeTrialTraction(dof_data[dof_idx],
                                                    Q_self_plus_can,
                                                    Q_self_minus_can,
                                                    sigma_n_trial,
                                                    tau1_trial, tau2_trial);
               p.sigma_n_trial_qp[q] = sigma_n_trial;
               p.tau1_trial_qp[q]    = tau1_trial;
               p.tau2_trial_qp[q]    = tau2_trial;
            }
         }
         probes.push_back(std::move(p));
      }

      // T1: max over (c, q) of spread across triangles of Q_self_{plus,minus}_rot[c,q].
      ChainRecord r{step, 0,0,-1,-1, -1,-1, 0,0,0, 0,0,0, 0,0,0};
      for (int c = 0; c < NUM_STATE; c++)
      {
         for (int q = 0; q < nqp; q++)
         {
            std::vector<real_t> vp, vm;
            for (const auto &p : probes)
            {
               if (!p.Q_self_plus_rot.empty())
               { vp.push_back(p.Q_self_plus_rot[c * nqp + q]); }
               if (!p.Q_self_minus_rot.empty())
               { vm.push_back(p.Q_self_minus_rot[c * nqp + q]); }
            }
            const real_t sp = SpreadMaxMin(vp);
            const real_t sm = SpreadMaxMin(vm);
            if (sp > r.T1_plus_spread_max)
            { r.T1_plus_spread_max = sp; r.T1_plus_worst_c = c; r.T1_plus_worst_qp = q; }
            if (sm > r.T1_minus_spread_max)
            { r.T1_minus_spread_max = sm; r.T1_minus_worst_c = c; r.T1_minus_worst_qp = q; }
         }
      }

      // T2: max across QPs of per-QP spread of trial tractions across triangles.
      for (int q = 0; q < nqp; q++)
      {
         std::vector<real_t> v_sn, v_t1, v_t2;
         for (const auto &p : probes)
         {
            if (!p.sigma_n_trial_qp.empty())
            {
               v_sn.push_back(p.sigma_n_trial_qp[q]);
               v_t1.push_back(p.tau1_trial_qp[q]);
               v_t2.push_back(p.tau2_trial_qp[q]);
            }
         }
         r.T2_sigma_n_spread = std::max(r.T2_sigma_n_spread, SpreadMaxMin(v_sn));
         r.T2_tau1_spread    = std::max(r.T2_tau1_spread,    SpreadMaxMin(v_t1));
         r.T2_tau2_spread    = std::max(r.T2_tau2_spread,    SpreadMaxMin(v_t2));
      }

      // T3: DOFData spread (global, matches v9.4.0 §11).
      r.T3_tau1_corr_spread =
         GlobalSpread(dof_data, [](const DOFData &d){return d.tau1_corr;});
      r.T3_tau2_corr_spread =
         GlobalSpread(dof_data, [](const DOFData &d){return d.tau2_corr;});
      r.T3_sigma_n_corr_spread =
         GlobalSpread(dof_data, [](const DOFData &d){return d.sigma_n_corr;});
      // T3b: V-decomposition sub-stage.
      r.T3b_V1_spread =
         GlobalSpread(dof_data, [](const DOFData &d){return d.V1;});
      r.T3b_V2_spread =
         GlobalSpread(dof_data, [](const DOFData &d){return d.V2;});
      r.T3b_slip_rate_spread =
         GlobalSpread(dof_data, [](const DOFData &d){return d.slip_rate;});

      std::cout << "  " << std::setw(4) << step
                << " | " << std::setw(8) << r.T1_plus_spread_max
                << " | " << std::setw(8) << r.T1_minus_spread_max
                << " | " << std::setw(10) << r.T2_sigma_n_spread
                << " | " << std::setw(8) << r.T2_tau1_spread
                << " | " << std::setw(8) << r.T2_tau2_spread
                << " | " << std::setw(8) << r.T3_tau1_corr_spread
                << " | " << std::setw(8) << r.T3_tau2_corr_spread
                << " | " << std::setw(9) << r.T3_sigma_n_corr_spread
                << " | " << std::setw(7) << r.T3b_V1_spread
                << " | " << std::setw(7) << r.T3b_V2_spread
                << " | " << std::setw(13) << r.T3b_slip_rate_spread << "\n";
      records.push_back(r);

      // R-001/R-003 fix (round-8 Step 0): per-component T1 dump at
      // MULTIPLE steps (1, 5, 10, 15, 19) — resolves whether SXY/VY
      // stay ULP-clean or grow into Pa-scale, and whether the
      // "y-antisymmetric clean" narrative holds across timesteps.
      const bool dump_this_step = (step == 1 || step == 5 || step == 10
                                    || step == 15 || step == 19);
      if (dump_this_step)
      {
         std::cout << "\n  [R8 step-" << step << " per-component T1 breakdown]\n";
         std::cout << "    per-side, per-component Q_self canonical-frame spread\n";
         std::cout << "    comp  |  plus-side spread  |  minus-side spread\n";
         std::cout << "    ------|--------------------|-------------------\n";
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t sp_max = 0, sm_max = 0;
            for (int q = 0; q < nqp; q++)
            {
               std::vector<real_t> vp, vm;
               for (const auto &p : probes)
               {
                  if (!p.Q_self_plus_rot.empty())
                  { vp.push_back(p.Q_self_plus_rot[c * nqp + q]); }
                  if (!p.Q_self_minus_rot.empty())
                  { vm.push_back(p.Q_self_minus_rot[c * nqp + q]); }
               }
               sp_max = std::max(sp_max, SpreadMaxMin(vp));
               sm_max = std::max(sm_max, SpreadMaxMin(vm));
            }
            std::cout << "    " << std::setw(5) << kCompName[c]
                      << " |  " << std::scientific << std::setprecision(3)
                      << std::setw(17) << sp_max
                      << " |  " << std::setw(17) << sm_max << "\n";
         }
         std::cout << std::setprecision(6) << "\n";
      }
   }

   // Worst-case summary.
   std::cout << "\n=== Worst spreads across " << kNSteps << " steps ===\n";
   real_t w_T1p = 0, w_T1m = 0, w_T2_sn = 0, w_T2_t1 = 0, w_T2_t2 = 0;
   real_t w_T3_t1 = 0, w_T3_t2 = 0, w_T3_sn = 0;
   real_t w_T3b_V1 = 0, w_T3b_V2 = 0, w_T3b_sr = 0;
   int w_T1p_c = -1, w_T1p_q = -1, w_T1m_c = -1, w_T1m_q = -1;
   int w_T1p_step = -1, w_T1m_step = -1, w_T2_step = -1, w_T3_step = -1;
   for (const auto &r : records)
   {
      if (r.T1_plus_spread_max > w_T1p)
      { w_T1p = r.T1_plus_spread_max; w_T1p_c = r.T1_plus_worst_c;
        w_T1p_q = r.T1_plus_worst_qp; w_T1p_step = r.step; }
      if (r.T1_minus_spread_max > w_T1m)
      { w_T1m = r.T1_minus_spread_max; w_T1m_c = r.T1_minus_worst_c;
        w_T1m_q = r.T1_minus_worst_qp; w_T1m_step = r.step; }
      if (r.T2_sigma_n_spread > w_T2_sn || r.T2_tau1_spread > w_T2_t1
          || r.T2_tau2_spread > w_T2_t2)
      { w_T2_step = r.step; }
      w_T2_sn = std::max(w_T2_sn, r.T2_sigma_n_spread);
      w_T2_t1 = std::max(w_T2_t1, r.T2_tau1_spread);
      w_T2_t2 = std::max(w_T2_t2, r.T2_tau2_spread);
      if (r.T3_tau1_corr_spread > w_T3_t1) { w_T3_step = r.step; }
      w_T3_t1 = std::max(w_T3_t1, r.T3_tau1_corr_spread);
      w_T3_t2 = std::max(w_T3_t2, r.T3_tau2_corr_spread);
      w_T3_sn = std::max(w_T3_sn, r.T3_sigma_n_corr_spread);
      w_T3b_V1 = std::max(w_T3b_V1, r.T3b_V1_spread);
      w_T3b_V2 = std::max(w_T3b_V2, r.T3b_V2_spread);
      w_T3b_sr = std::max(w_T3b_sr, r.T3b_slip_rate_spread);
   }

   // R-002 fix: growth-rate characterization.
   {
      std::cout << "\n=== Growth-rate characterization (R-002 fix) ===\n";
      auto analyze_growth = [](const char *label, const std::vector<real_t> &xs)
      {
         std::cout << "  " << label << ":\n";
         if (xs.size() < 2) { std::cout << "    (need >=2 steps)\n"; return; }
         std::cout << "    per-step ratios x[n]/x[n-1]:";
         for (std::size_t n = 1; n < xs.size(); n++)
         {
            real_t r = (xs[n-1] > 0) ? xs[n]/xs[n-1] : 0.0;
            std::cout << " " << std::fixed << std::setprecision(3) << r;
         }
         std::cout << std::scientific << std::setprecision(6) << "\n";
         // Power-law fit: x ~ a * t^b (log-log linear regression).
         // Skip t=0 if x[0] == 0.
         std::vector<real_t> lt, lx;
         for (std::size_t n = 0; n < xs.size(); n++)
         {
            if (xs[n] > 0)
            { lt.push_back(std::log(real_t(n+1))); lx.push_back(std::log(xs[n])); }
         }
         if (lt.size() >= 2)
         {
            const int N = (int)lt.size();
            real_t sx=0, sy=0, sxx=0, sxy=0;
            for (int i = 0; i < N; i++)
            { sx += lt[i]; sy += lx[i]; sxx += lt[i]*lt[i]; sxy += lt[i]*lx[i]; }
            const real_t b = (N*sxy - sx*sy) / (N*sxx - sx*sx);
            const real_t a = std::exp((sy - b*sx)/N);
            std::cout << "    power-law fit: x ~ " << a << " * t^"
                      << std::fixed << std::setprecision(3) << b
                      << std::scientific << std::setprecision(6) << "\n";
            const real_t linear_rel = (xs[0] * xs.size() > 0)
               ? xs.back() / (xs[0] * xs.size()) : 0.0;
            const real_t power_pred = a * std::pow(real_t(xs.size()), b);
            std::cout << "    observed x[N]: " << xs.back()
                      << "   linear pred: " << xs[0] * xs.size()
                      << "   (linear/obs ratio: " << std::fixed
                      << std::setprecision(2) << linear_rel
                      << std::scientific << std::setprecision(6) << ")\n";
            std::cout << "    power-law pred at t=" << xs.size()
                      << ": " << power_pred << "\n";
         }
      };
      // Gather per-step series.
      std::vector<real_t> s_T1p, s_T1m, s_T2sn, s_T2t1, s_T2t2,
                            s_T3t1, s_T3t2, s_T3sn;
      for (const auto &r : records)
      {
         s_T1p.push_back(r.T1_plus_spread_max);
         s_T1m.push_back(r.T1_minus_spread_max);
         s_T2sn.push_back(r.T2_sigma_n_spread);
         s_T2t1.push_back(r.T2_tau1_spread);
         s_T2t2.push_back(r.T2_tau2_spread);
         s_T3t1.push_back(r.T3_tau1_corr_spread);
         s_T3t2.push_back(r.T3_tau2_corr_spread);
         s_T3sn.push_back(r.T3_sigma_n_corr_spread);
      }
      analyze_growth("T1 plus-side max",  s_T1p);
      analyze_growth("T1 minus-side max", s_T1m);
      analyze_growth("T2 sigma_n_trial",  s_T2sn);
      analyze_growth("T2 tau1_trial",     s_T2t1);
      analyze_growth("T2 tau2_trial",     s_T2t2);
      analyze_growth("T3 tau1_corr",      s_T3t1);
      analyze_growth("T3 tau2_corr",      s_T3t2);
      analyze_growth("T3 sigma_n_corr",   s_T3sn);
   }

   auto cn = [&](int c){ return (c>=0 && c<NUM_STATE) ? kCompName[c] : "?"; };
   std::cout << "  T1 plus-side max spread   : " << w_T1p
             << "  (comp=" << cn(w_T1p_c) << " qp=" << w_T1p_q
             << " step=" << w_T1p_step << ")\n";
   std::cout << "  T1 minus-side max spread  : " << w_T1m
             << "  (comp=" << cn(w_T1m_c) << " qp=" << w_T1m_q
             << " step=" << w_T1m_step << ")\n";
   std::cout << "  T2 sigma_n_trial spread   : " << w_T2_sn
             << "  (step=" << w_T2_step << ")\n";
   std::cout << "  T2 tau1_trial spread      : " << w_T2_t1 << "\n";
   std::cout << "  T2 tau2_trial spread      : " << w_T2_t2 << "\n";
   std::cout << "  T3 tau1_corr spread       : " << w_T3_t1
             << "  (step=" << w_T3_step << ")\n";
   std::cout << "  T3 tau2_corr spread       : " << w_T3_t2 << "\n";
   std::cout << "  T3 sigma_n_corr spread    : " << w_T3_sn << "\n";
   std::cout << "  T3b V1 spread             : " << w_T3b_V1 << "\n";
   std::cout << "  T3b V2 spread             : " << w_T3b_V2 << "\n";
   std::cout << "  T3b slip_rate spread      : " << w_T3b_sr << "\n";

   std::cout << "\n=== Decision-tree verdict ===\n";
   const real_t kTolT1 = 1e-6;  // relative to 1e+6 Pa stress scale, 1e-12 rel
   const bool T1_clean = (w_T1p < kTolT1) && (w_T1m < kTolT1);
   const bool T2_clean = (w_T2_sn < kTolT1) && (w_T2_t1 < kTolT1)
                          && (w_T2_t2 < kTolT1);
   const bool T3_dirty = (w_T3_t1 > 1e-10) || (w_T3_t2 > 1e-10)
                          || (w_T3_sn > 1e-10);
   std::cout << "  T1 clean? " << (T1_clean ? "yes" : "NO") << "\n";
   std::cout << "  T2 clean? " << (T2_clean ? "yes" : "NO") << "\n";
   std::cout << "  T3 dirty? " << (T3_dirty ? "yes" : "no") << "\n";

   if (T1_clean && T2_clean && T3_dirty)
   {
      std::cout << "\n  → Branch: T1 ≈ 0 with full coverage AND T2 ≈ 0 AND T3 > 0.\n"
                << "    Asymmetry is in stages 4-6 of FaultFaceFlux::Evaluate:\n"
                << "      stage 4: Brent friction solve (fault_face_flux.cpp:152-157)\n"
                << "      stage 5: V-decomposition (fault_face_flux.cpp:163-177)\n"
                << "      stage 6: tau1_corr writeback (fault_face_flux.cpp:222-224)\n"
                << "    Further disambiguation via T3b (V1/V2 vs tau1_corr):\n";
      if (w_T3b_V1 < kTolT1 && w_T3_t1 > kTolT1)
      {
         std::cout << "      T3b V1 clean, T3 tau1_corr dirty → writeback line 222-223\n"
                   << "      (tau1_corr = tau1_0 + tau1_nuc + tau1_corr).\n";
      }
      else if (w_T3b_V1 > kTolT1 && w_T3b_V2 > kTolT1)
      {
         std::cout << "      T3b V1 and V2 both dirty → V-decomposition\n"
                   << "      (lines 171-172: V1 = V_abs * tau1_total / (strength + eta_s*V_abs)).\n";
      }
   }
   else if (T1_clean && !T2_clean)
   {
      std::cout << "\n  → Branch: T1 ≈ 0 AND T2 > 0.\n"
                << "    Asymmetry is in stages 1-3:\n"
                << "      stage 1: ComputeTrialTraction (fault_face_flux.cpp:39-68)\n"
                << "      stage 2: Total-traction sum (fault_face_flux.cpp:126-128)\n"
                << "      stage 3: Theta = sqrt(tau1² + tau2²) (fault_face_flux.cpp:149).\n";
   }
   else if (!T1_clean)
   {
      std::cout << "\n  → Branch: T1 > 0 — DG layer NOT clean.\n"
                << "    Reopen at DG layer.  Target: shape·Q sampling at fault QPs,\n"
                << "    or the per-face CalcOrtho normal (wave_operator.inl:1186-1204).\n";
   }

   return 0;
}
