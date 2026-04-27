// Round-11 R-1102 — end-to-end Adjacent-mode dispatch correctness.
//
// All other mixed-flux tests stop at set construction or the flux
// primitive in isolation; none calls `wave.AdvanceADER` with
// `MixedFluxMode != None`.  This test exercises the full dispatch
// pipeline on a non-uniform initial state across a fixture where
// `|central_flux_face_set_| > 0` and verifies that:
//
//   Gate 1: Q_new_adj differs from Q_new_none by a non-trivial amount
//           (the central dispatch is actually firing).  Quantitatively,
//           max |Q_new_adj − Q_new_none| > scale_dependent_floor.
//
//   Gate 2: the difference is bounded above by the analytic estimate
//           |0.5·|A_n|·jump| × |central_flux_face_set_| × dt × M^{-1}
//           order of magnitude — basically a sanity ceiling so the
//           test catches blowups that masquerade as "dispatch fired".
//
//   Gate 3: switching back to None on the SAME wave operator and
//           re-stepping produces a Q_new bit-identical to a fresh
//           never-touched-Adjacent run (idempotency / clean-up
//           regression — guards against latent state from the
//           Adjacent build that survives the None reset).
//
// Fixture: BuildSmallFaultTetMesh (48 tets — 8 hexes × 6-tet split;
// fault at y=L/2; |central_flux_face_set_| in Adjacent mode is 32 —
// the same fixture already verified in test_mixed_flux_face_set).  Initial state has
// per-DOF stress that varies linearly with (x, y, z), producing
// non-trivial jumps across every interior face.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
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

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_GE(v, lo, msg) do { \
   num_tests++; \
   double vv = (v), ll = (lo); \
   if (vv >= ll) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", floor " << ll << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected >= " << ll << ")\n"; } \
} while (0)

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

namespace
{

// Same as BuildSmallFaultTetMesh from test_mixed_flux_face_set.cpp:
// 48-tet 3x3x3 vertex grid with interior fault plane at y = L/2.
Mesh BuildSmallFaultTetMesh(real_t L = 1000.0)
{
   const int Nv = 3 * 3 * 3;
   Mesh mesh(3, Nv, 0, 0);
   auto vid = [&](int i, int j, int k) { return i + 3*j + 9*k; };
   for (int k = 0; k < 3; k++)
      for (int j = 0; j < 3; j++)
         for (int i = 0; i < 3; i++)
         {
            real_t v[3] = {i*0.5*L, j*0.5*L, k*0.5*L};
            mesh.AddVertex(v);
         }
   auto vfunc = [&](int ci, int cj, int ck, int lv) {
      return vid(ci + (lv & 1), cj + ((lv >> 1) & 1), ck + ((lv >> 2) & 1));
   };
   const int pat[6][4] = {
      {0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7},
      {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}
   };
   for (int ck = 0; ck < 2; ck++)
      for (int cj = 0; cj < 2; cj++)
         for (int ci = 0; ci < 2; ci++)
            for (int t = 0; t < 6; t++)
            {
               mesh.AddTet(vfunc(ci, cj, ck, pat[t][0]),
                           vfunc(ci, cj, ck, pat[t][1]),
                           vfunc(ci, cj, ck, pat[t][2]),
                           vfunc(ci, cj, ck, pat[t][3]),
                           1);
            }
   mesh.FinalizeTopology();
   const real_t eps = 1e-6;
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      real_t cy = 0.0;
      for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
      cy /= fv.Size();
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cy - 0.5 * L) < eps)
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

real_t MaxAbsDiff(const Vector &a, const Vector &b)
{
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      const real_t d = std::abs(a(i) - b(i));
      if (d > m) { m = d; }
   }
   return m;
}

real_t MaxAbs(const Vector &a)
{
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      const real_t d = std::abs(a(i));
      if (d > m) { m = d; }
   }
   return m;
}

// Set up a non-uniform initial state that is DISCONTINUOUS across
// element interfaces.  In DG with order=1, vertex DOFs at a shared vertex
// have the same coordinate → a coord-based IC is L2-continuous and
// produces zero jump.  We instead set Q PER ELEMENT (constant within each
// element, varying across elements via a deterministic per-element
// formula), guaranteeing jumps at every interior face.
void InitializeNonUniformState(Vector &Q, const FiniteElementSpace &fes,
                               int ndof_total, real_t /*L*/)
{
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
   real_t *Q_data = Q.GetData();
   for (int e = 0; e < fes.GetNE(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      const int ndof_e = fe->GetDof();
      Array<int> dof_idx;
      fes.GetElementDofs(e, dof_idx);

      // Per-element pseudo-random scaling (deterministic).
      const real_t s = 1.0 + 0.1 * static_cast<real_t>((e * 17 + 5) % 23);
      const real_t t = 1.0 + 0.05 * static_cast<real_t>((e * 31 + 11) % 19);

      for (int i = 0; i < ndof_e; i++)
      {
         const int dof = dof_idx[i];
         Q_data[SXX * ndof_total + dof] = -TPV102Params::sigma_n * s;
         Q_data[SYY * ndof_total + dof] =  TPV102Params::sigma_n * t;
         Q_data[SXY * ndof_total + dof] = -TPV102Params::tau_ini * s;
         Q_data[SZZ * ndof_total + dof] =  1.0e6 * t;
         Q_data[SYZ * ndof_total + dof] =  1.0e6 * s;
         Q_data[SXZ * ndof_total + dof] =  1.0e6 * t;
         Q_data[VX  * ndof_total + dof] =  0.1 * s;
         Q_data[VY  * ndof_total + dof] =  0.1 * t;
         Q_data[VZ  * ndof_total + dof] =  0.1 * (s + t);
      }
   }
}

void SetupWaveForRun(WaveOperator<Mesh> &wave,
                     std::vector<DOFData> &dof_data,
                     const Mesh &mesh, int order)
{
   const int n_local_qps = wave.GetNumLocalFaultQPs();
   const int nbf = wave.GetNbfPerFace();

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   std::vector<Vector> fault_coords;
   fault_coords.reserve(n_local_qps);
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr =
         const_cast<Mesh &>(mesh).GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(
         ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   InitializeFaultDOFs(dof_data, n_local_qps, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, n_local_qps);
   wave.SetFaultDOFData(&dof_data, nbf);

   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);
}

void RunOneStepGivenIC(WaveOperator<Mesh> &wave,
                       std::vector<DOFData> &dof_data,
                       const Mesh &mesh, int order,
                       const Vector &Q_init,
                       Vector &Q_new_out)
{
   SetupWaveForRun(wave, dof_data, mesh, order);
   const real_t dt = 1e-5;
   wave.AdvanceADER(Q_init, dt, /*ader_order=*/2, Q_new_out);
}

// R-1200: assemble the analytic prediction of dQdt_adj − dQdt_none
// element-locally.  For each face f in central_flux_face_set_, for
// each QP q on f, compute F_int(Q_self, Q_nbr) − F_central(Q_self,
// Q_nbr).  Accumulate +w·shape1·δF into elem1's residual block and
// −w·shape2·δF into elem2's residual block.  Then apply the
// element-local mass inverse via WaveOperator::GetElementMassInverse.
//
// The result must equal Mult(Q_init, mode=Adjacent) − Mult(Q_init, mode=None)
// to FP precision (1e-12 relative) since the only mode-dependent term
// in Mult is the per-face F_h selection at central faces.  Volume RHS
// and mass inverse are both mode-independent.
void ComputeAnalyticAdjacentDelta(const WaveOperator<Mesh> &wave,
                                  const Mesh &mesh, int order,
                                  const Vector &Q_init,
                                  Vector &delta_out)
{
   const auto &fes = wave.GetFESpace();
   const int ndof_total  = fes.GetNDofs();
   const int ndof_per_el = wave.GetNDof();
   delta_out.SetSize(NUM_STATE * ndof_total);
   delta_out = 0.0;

   const real_t *Q_data = Q_init.GetData();
   const auto &flux = wave.GetFlux();
   const auto &central_set = wave.GetCentralFluxFaceSet();

   for (int f : central_set)
   {
      auto *ftr =
         const_cast<Mesh &>(mesh).GetInteriorFaceTransformations(f);
      if (!ftr || ftr->Elem2No < 0) { continue; }
      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      const FiniteElement *fe1 = fes.GetFE(e1);
      const FiniteElement *fe2 = fes.GetFE(e2);
      const int ndof_e1 = fe1->GetDof();
      const int ndof_e2 = fe2->GetDof();
      const int dof_offset1 = e1 * ndof_per_el;
      const int dof_offset2 = e2 * ndof_per_el;

      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0) { nor_vec /= nor_len; }
         const real_t w = ip.weight * nor_len;
         const real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         Vector shape1(ndof_e1), shape2(ndof_e2);
         fe1->CalcShape(ftr->GetElement1IntPoint(), shape1);
         fe2->CalcShape(ftr->GetElement2IntPoint(), shape2);

         real_t Q_self[NUM_STATE], Q_nbr[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_self[c] = 0.0;
            Q_nbr[c]  = 0.0;
            for (int i = 0; i < ndof_e1; i++)
            {
               Q_self[c] += shape1(i) *
                  Q_data[c * ndof_total + dof_offset1 + i];
            }
            for (int i = 0; i < ndof_e2; i++)
            {
               Q_nbr[c]  += shape2(i) *
                  Q_data[c * ndof_total + dof_offset2 + i];
            }
         }

         real_t F_int[NUM_STATE], F_ce[NUM_STATE];
         flux.Interior(nor, Q_self, Q_nbr, F_int);
         flux.Central (nor, Q_self, Q_nbr, F_ce);

         for (int c = 0; c < NUM_STATE; c++)
         {
            const real_t dF = F_int[c] - F_ce[c];
            for (int i = 0; i < ndof_e1; i++)
            {
               delta_out(c * ndof_total + dof_offset1 + i) +=
                  w * shape1(i) * dF;
            }
            for (int i = 0; i < ndof_e2; i++)
            {
               delta_out(c * ndof_total + dof_offset2 + i) -=
                  w * shape2(i) * dF;
            }
         }
      }
   }

   // Apply element-local mass inverse: delta_out ← M_e^{-1} · delta_out per e.
   const int ne = mesh.GetNE();
   for (int e = 0; e < ne; e++)
   {
      const DenseMatrix &Minv = wave.GetElementMassInverse(e);
      const int dof_off = e * ndof_per_el;
      Vector tmp_in(ndof_per_el), tmp_out(ndof_per_el);
      for (int c = 0; c < NUM_STATE; c++)
      {
         for (int i = 0; i < ndof_per_el; i++)
         {
            tmp_in(i) = delta_out(c * ndof_total + dof_off + i);
         }
         Minv.Mult(tmp_in, tmp_out);
         for (int i = 0; i < ndof_per_el; i++)
         {
            delta_out(c * ndof_total + dof_off + i) = tmp_out(i);
         }
      }
   }
}

} // anonymous

int main()
{
   std::cout << "\n=== Round-11 R-1102: Adjacent-mode dispatch correctness ===\n";

   const int order = 1;
   const real_t L = 1000.0;

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;
   bc.absorbing_attrs = {};

   // ----------------------------------------------------------------
   // Path A: mode None.  Q_new_none = baseline.
   // ----------------------------------------------------------------
   Mesh mesh_a = BuildSmallFaultTetMesh(L);
   WaveOperator<Mesh> wave_a(mesh_a, order,
                             TPV102Params::lambda, TPV102Params::mu,
                             TPV102Params::rho, bc);
   FaultFaceFlux ff_a(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave_a.SetFaultFlux(&ff_a);
   std::vector<DOFData> dof_a;

   const auto &fes_a = wave_a.GetFESpace();
   const int ndof_total = fes_a.GetNDofs();

   Vector Q_init(NUM_STATE * ndof_total);
   InitializeNonUniformState(Q_init, fes_a, ndof_total, L);

   Vector Q_new_none(NUM_STATE * ndof_total);
   RunOneStepGivenIC(wave_a, dof_a, mesh_a, order, Q_init, Q_new_none);

   const real_t scale_q_new = MaxAbs(Q_new_none);
   std::cout << "  max |Q_new_none| = " << std::scientific
             << std::setprecision(3) << scale_q_new << "\n";

   // ----------------------------------------------------------------
   // Path B: mode Adjacent on a fresh wave.  Q_new_adj.
   // ----------------------------------------------------------------
   Mesh mesh_b = BuildSmallFaultTetMesh(L);
   WaveOperator<Mesh> wave_b(mesh_b, order,
                             TPV102Params::lambda, TPV102Params::mu,
                             TPV102Params::rho, bc);
   FaultFaceFlux ff_b(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave_b.SetFaultFlux(&ff_b);
   std::vector<DOFData> dof_b;

   wave_b.SetMixedFluxMode(MixedFluxMode::Adjacent);
   const std::size_t central_set_size = wave_b.GetCentralFluxFaceSet().size();
   std::cout << "  |central_flux_face_set_| (Adjacent) = "
             << central_set_size << "\n";

   Vector Q_new_adj(NUM_STATE * ndof_total);
   RunOneStepGivenIC(wave_b, dof_b, mesh_b, order, Q_init, Q_new_adj);

   const real_t diff_adj_vs_none = MaxAbsDiff(Q_new_adj, Q_new_none);
   std::cout << "  max |Q_new_adj − Q_new_none| = " << std::scientific
             << std::setprecision(3) << diff_adj_vs_none << "\n";

   // Gate 1: dispatch fired (non-trivial diff).  Floor: 1e-9 × scale.
   // Pure FP noise from re-running AdvanceADER would be ~1e-15 × scale.
   // A genuine central-flux delta on this fixture should be well above
   // that — we expect at least ~|A_n|·jump·dt scale, i.e., 1e3 · 1e7 ·
   // 1e-5 = 1e5 absolute, which is ~1e-3 × scale_q_new (since scale ~1e8).
   // Use 1e-9 × scale_q_new as a conservative lower bound that still
   // distinguishes real dispatch from FP noise.
   TEST_GE(diff_adj_vs_none / scale_q_new, 1e-9,
           "Gate 1: Q_new_adj differs from Q_new_none — Adjacent dispatch fires");
   TEST_GE(static_cast<double>(central_set_size), 1.0,
           "Gate 1b: central_flux_face_set_ is non-empty");

   // Gate 2: bounded above by reasonable analytic ceiling.
   // |delta| <= |A_n| · |jump| · dt · |M^{-1}|.  |A_n| ~ cp = 6e3,
   // |jump| ~ 1e7 (stress), dt = 1e-5, |M^{-1}| ~ 1/h^3 / vol ~ O(1)
   // for a unit element.  So |delta| < 1e6 absolute.  Use 1.0 × scale_q_new
   // (i.e., delta < ~10× of Q_new itself) as a sanity ceiling — catches
   // unbounded blowups but allows substantial difference.
   TEST_LE(diff_adj_vs_none / scale_q_new, 1.0,
           "Gate 2: Q_new_adj − Q_new_none is bounded (no blowup)");

   // -----------------------------------------------------------------
   // R-1200: Gate 2-analytic — verify the Phase 4 plan acceptance
   // criterion (algebraic identity).  Through the RK1 path Mult, the
   // ONLY mode-dependent term is the central-face F_h selection, so:
   //   dQdt_adj − dQdt_none
   //     = M^{-1} · sum_{f in central_set, q on f}
   //         { +w · shape1 · (F_int − F_central)   for elem1
   //           −w · shape2 · (F_int − F_central)   for elem2 }
   // The coarse Gate 1 (1e-9 floor) and Gate 2 (1.0 ceiling) cannot
   // detect a sign-flipped or wrong-magnitude Central; Gate 2-analytic
   // closes that hole at 1e-12 relative.
   // -----------------------------------------------------------------
   Mesh mesh_d_none = BuildSmallFaultTetMesh(L);
   WaveOperator<Mesh> wave_d_none(mesh_d_none, order,
                                   TPV102Params::lambda, TPV102Params::mu,
                                   TPV102Params::rho, bc);
   FaultFaceFlux ff_d_none(TPV102Params::rho, TPV102Params::cp,
                           TPV102Params::cs);
   wave_d_none.SetFaultFlux(&ff_d_none);
   std::vector<DOFData> dof_d_none;
   SetupWaveForRun(wave_d_none, dof_d_none, mesh_d_none, order);
   // mode None (default).
   Vector dQdt_none(NUM_STATE * ndof_total);
   wave_d_none.Mult(Q_init, dQdt_none);

   Mesh mesh_d_adj = BuildSmallFaultTetMesh(L);
   WaveOperator<Mesh> wave_d_adj(mesh_d_adj, order,
                                  TPV102Params::lambda, TPV102Params::mu,
                                  TPV102Params::rho, bc);
   FaultFaceFlux ff_d_adj(TPV102Params::rho, TPV102Params::cp,
                          TPV102Params::cs);
   wave_d_adj.SetFaultFlux(&ff_d_adj);
   std::vector<DOFData> dof_d_adj;
   SetupWaveForRun(wave_d_adj, dof_d_adj, mesh_d_adj, order);
   wave_d_adj.SetMixedFluxMode(MixedFluxMode::Adjacent);
   Vector dQdt_adj(NUM_STATE * ndof_total);
   wave_d_adj.Mult(Q_init, dQdt_adj);

   Vector observed_delta(NUM_STATE * ndof_total);
   for (int i = 0; i < observed_delta.Size(); i++)
   {
      observed_delta(i) = dQdt_adj(i) - dQdt_none(i);
   }
   const real_t observed_scale = MaxAbs(observed_delta);
   std::cout << "  max |dQdt_adj − dQdt_none| (Mult/RK1) = "
             << std::scientific << std::setprecision(3)
             << observed_scale << "\n";

   Vector expected_delta(NUM_STATE * ndof_total);
   ComputeAnalyticAdjacentDelta(wave_d_adj, mesh_d_adj, order, Q_init,
                                expected_delta);

   const real_t analytic_residual = MaxAbsDiff(observed_delta, expected_delta);
   const real_t analytic_rel = analytic_residual /
                               std::max<real_t>(observed_scale,
                                                static_cast<real_t>(1.0));
   std::cout << "  max |observed − expected| / |observed|_max = "
             << std::scientific << std::setprecision(3) << analytic_rel
             << "\n";
   TEST_LE(analytic_rel, 1.0e-12,
           "Gate 2-analytic (R-1200): dQdt_adj − dQdt_none equals the "
           "per-face Σ M^{-1} w shape (F_int − F_central) sum to 1e-12 "
           "relative (Phase 4 algebraic identity; catches sign-flip / "
           "wrong-magnitude Central bugs the coarse gates miss)");
   TEST_GE(observed_scale, 1.0,
           "Gate 2-analytic-floor: dispatch actually fires at the Mult "
           "(RK1) level too — observed delta is non-trivial");

   // ----------------------------------------------------------------
   // R-1501: Gate 2-analytic-ADER — verify the same algebraic identity
   // through the AdvanceADER (production) code path.  AdvanceADER's
   // flux assembly at wave_operator.inl:3993 (local) and :4478 (shared)
   // is structurally distinct from Mult's at L2510 / L3076; a
   // sign-flipped or factor-of-(1+ε) bug specific to those sites slips
   // past the Mult-only Gate 2-analytic above.
   //
   // Linearity argument: AdvanceADER computes
   //   I = ∫_0^dt Q(t+τ) dτ      (mode-independent — predictor is
   //                              volume-only via CK recursion)
   //   rhs = vol(I) − ∫ φ · F_h(I) dA   (per-face)
   //   Q_new = Q + M⁻¹ · rhs
   // The ONLY mode-dependent term is F_h at central faces.  So
   //   Q_new_adj − Q_new_none =
   //     M⁻¹ · sum_{f ∈ central} sum_q ±w · shape · (F_int(I) − F_central(I))
   // The same `ComputeAnalyticAdjacentDelta` helper computes this when
   // fed `I` instead of `Q_init` (the helper is mode-/state-neutral).
   // ----------------------------------------------------------------
   {
      // R-1604: construct FRESH WaveOperator instances rather than
      // reusing wave_d_none / wave_d_adj which already had Mult called.
      // Mult's mutable scratch buffers are designed to be reset per
      // call, but a future refactor that adds non-idempotent state
      // (e.g., a per-call counter or accumulator) would silently
      // corrupt the AdvanceADER reference.  Test isolation is cheap
      // insurance.
      Mesh mesh_e_none = BuildSmallFaultTetMesh(L);
      WaveOperator<Mesh> wave_e_none(mesh_e_none, order,
                                      TPV102Params::lambda, TPV102Params::mu,
                                      TPV102Params::rho, bc);
      FaultFaceFlux ff_e_none(TPV102Params::rho, TPV102Params::cp,
                               TPV102Params::cs);
      wave_e_none.SetFaultFlux(&ff_e_none);
      std::vector<DOFData> dof_e_none;
      SetupWaveForRun(wave_e_none, dof_e_none, mesh_e_none, order);

      Mesh mesh_e_adj = BuildSmallFaultTetMesh(L);
      WaveOperator<Mesh> wave_e_adj(mesh_e_adj, order,
                                     TPV102Params::lambda, TPV102Params::mu,
                                     TPV102Params::rho, bc);
      FaultFaceFlux ff_e_adj(TPV102Params::rho, TPV102Params::cp,
                              TPV102Params::cs);
      wave_e_adj.SetFaultFlux(&ff_e_adj);
      std::vector<DOFData> dof_e_adj;
      SetupWaveForRun(wave_e_adj, dof_e_adj, mesh_e_adj, order);
      wave_e_adj.SetMixedFluxMode(MixedFluxMode::Adjacent);

      const real_t dt_ader = 1.0e-5;
      Vector Q_new_none_ader(NUM_STATE * ndof_total);
      Vector Q_new_adj_ader(NUM_STATE * ndof_total);
      wave_e_none.AdvanceADER(Q_init, dt_ader, /*ader_order=*/2, Q_new_none_ader);
      wave_e_adj .AdvanceADER(Q_init, dt_ader, /*ader_order=*/2, Q_new_adj_ader);

      Vector observed_delta_ader(NUM_STATE * ndof_total);
      for (int i = 0; i < observed_delta_ader.Size(); i++)
      {
         observed_delta_ader(i) = Q_new_adj_ader(i) - Q_new_none_ader(i);
      }
      const real_t scale_ader = MaxAbs(observed_delta_ader);

      Vector I(NUM_STATE * ndof_total);
      wave_e_adj.ComputeADERTimeIntegrated(Q_init, dt_ader, 2, I);
      Vector expected_delta_ader(NUM_STATE * ndof_total);
      ComputeAnalyticAdjacentDelta(wave_e_adj, mesh_e_adj, order, I,
                                   expected_delta_ader);

      const real_t residual_ader =
         MaxAbsDiff(observed_delta_ader, expected_delta_ader);
      const real_t rel_ader = residual_ader /
         std::max<real_t>(scale_ader, static_cast<real_t>(1.0));
      std::cout << "  max |Q_new_adj − Q_new_none|       (AdvanceADER) = "
                << std::scientific << std::setprecision(3) << scale_ader
                << "\n  max |observed − expected| / |observed|_max  (ADER)= "
                << std::scientific << std::setprecision(3) << rel_ader
                << "\n";
      // R-1603: tightened from 1e-10 → 1e-12, matching the Mult/RK1
      // gate.  The actual residual on this fixture is ~3.98e-13 (well
      // inside 1e-12); there is no known FP source that would push the
      // residual above 1e-12 — the helper feeds I (the same time-
      // integrated state used by ComputeADERFaceFluxRHS), so the
      // analytic computation matches AdvanceADER's flux assembly
      // step-for-step.  A 100× looser ADER tolerance was masking
      // small relative bugs (1e-11 ULP-or-two-scale defects) that the
      // Mult gate would catch.
      TEST_LE(rel_ader, 1.0e-12,
              "Gate 2-analytic-ADER (R-1501 + R-1603): Q_new_adj − Q_new_none "
              "equals Σ_f M⁻¹·w·shape·(F_int(I)−F_central(I)) at 1e-12 "
              "relative through the AdvanceADER (production) path; "
              "covers wave_operator.inl:3993/4478 dispatch sites");
      TEST_GE(scale_ader, 1.0,
              "Gate 2-analytic-ADER-floor: AdvanceADER central dispatch "
              "actually fires — observed delta is non-trivial");
   }

   // ----------------------------------------------------------------
   // Gate 3: idempotency under Adjacent → None transition on the
   // SAME wave operator (R-1103-companion).  Latent state from the
   // Adjacent build (e.g., a flag, a stale set) must NOT survive the
   // mode reset.
   // ----------------------------------------------------------------
   Mesh mesh_c = BuildSmallFaultTetMesh(L);
   WaveOperator<Mesh> wave_c(mesh_c, order,
                             TPV102Params::lambda, TPV102Params::mu,
                             TPV102Params::rho, bc);
   FaultFaceFlux ff_c(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave_c.SetFaultFlux(&ff_c);
   std::vector<DOFData> dof_c;

   wave_c.SetMixedFluxMode(MixedFluxMode::Adjacent);
   wave_c.SetMixedFluxMode(MixedFluxMode::None);
   TEST_LE(static_cast<double>(wave_c.GetCentralFluxFaceSet().size()), 0.0,
           "Gate 3a: SetMixedFluxMode(None) clears the set after Adjacent");

   Vector Q_new_c(NUM_STATE * ndof_total);
   RunOneStepGivenIC(wave_c, dof_c, mesh_c, order, Q_init, Q_new_c);

   const real_t diff_round_trip = MaxAbsDiff(Q_new_c, Q_new_none);
   std::cout << "  max |Q_new (Adjacent->None) − Q_new_none| = "
             << std::scientific << std::setprecision(3) << diff_round_trip << "\n";
   TEST_LE(diff_round_trip, 0.0,
           "Gate 3b: Q_new bit-identical after Adjacent->None transition");

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
