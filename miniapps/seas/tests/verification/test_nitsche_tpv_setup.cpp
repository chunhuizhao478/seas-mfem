// Nitsche dynamic rupture: the deck material, the bimaterial penalty and time
// step, the mu-scaled prestress, the instantaneous nucleation and the per-side
// stations (TPV6/TPV7/TPV31 plan rev 1, phase A).
//
//   (a) material: the TPV6 and TPV31 decks' materials, sampled per element by
//       NitscheMaterialSample, equal the deck values at each centroid; the
//       TPV31 jumps are found at 2400, 5000 and 10000 m; elements across a
//       jump are counted (none when it is a mesh surface, some when not);
//       the fault penalty sigma_F of every TPV6 fault face takes each side's
//       moduli, and the time step is the minimum over elements of
//       cfl h_e / (p^2 c_p,e) with each element's own c_p (on a copy of the
//       box with the far side squashed, so the slow side sets the step).
//   (b) prestress: the TPV31 deck's tensor times mu(x)/mu_0 through the
//       operator's projection gives along-strike shear 30 MPa mu/mu_0 and
//       sigma_n 60 MPa mu/mu_0 at every owned fault DOF of the constant
//       layers (each fault face lies in one layer of the box below).
//   (c) nucleation: 4.95 MPa mu/mu_0 x (1, 1, 0.5, 0, 0) at r = 0, 1400,
//       1700, 2000, 2100 m from the hypocentre, along strike and along dip;
//       no normal or dip traction.  The planar-fault check the driver applies
//       to both nucleation kinds: the largest normal deviation is round-off
//       on the box and above the 1e-6 limit on a bent copy (fault Y = a X^2).
//   (d) cohesion: the TPV31 deck's rule, 425 Pa/m x (2400 m - depth) above
//       2400 m (1.02, 0.51 and 0 MPa at 0, 1.2 and 2.4 km), reaches the
//       closure's parameters at every owned fault DOF.
//   (e) per-side writer: displacement and velocity fields linear on each
//       side with a jump at Y = 0 come out to round-off on each side, far =
//       Y > 0 for the reference normal (0, -1, 0), with h = +X, v = -Z,
//       n = +Y; the files hold the same numbers.
//
// The box: X, Y in [-30, 30] km, Z in [-30, 0] km, 24 tets per hex, the
// fault (attribute 3) on the whole plane Y = 0, its hex layers moved to the
// depths 0, 2.4, 5, 10, 20, 30 km so the TPV31 jumps are mesh planes.  On 2
// ranks the partition is the fault plane, so every fault face is shared and
// each side of a station lives on its own rank.  Run from miniapps/seas (the
// test reads the TPV6 and TPV31 decks).
//
// Scope: Chunhui Zhao, 2026-10-19.

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../dynamic/fault_side_stations.hpp"
#include "../../dynamic/nitsche_deck_fields.hpp"
#include "../../friction/slip_weakening_closure.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../spatial/code/spatial_friction.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kL = 30.0e3;
static const real_t kMu0 = 32.03812032e9;
static const char *kDeckTPV6 = "tpv6/configs/tpv6_nitsche.toml";
static const char *kDeckTPV31 = "tpv31/configs/tpv31_nitsche.toml";

static real_t RMax(real_t v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return v; }
static long long LSum(long long v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD); return v; }

// The layered fault box (header).
static Mesh *LayeredFaultBox(int nx, int ny)
{
   const real_t depth[6] = {0.0, 2400.0, 5000.0, 10000.0, 20000.0, 30000.0};
   const int nz = 5;
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(nx, ny, nz, 2 * kL, 2 * kL, kL));
   const real_t dz = kL / nz;
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL;
      c[1] -= kL;
      // grid depth g = kL - z in [0, kL], piecewise linear onto the layer depths
      const real_t g = kL - c[2];
      int k = std::min(nz - 1, static_cast<int>(std::floor(g / dz + 1e-9)));
      const real_t s = (g - k * dz) / dz;
      c[2] = -((1.0 - s) * depth[k] + s * depth[k + 1]);
   }
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      real_t zc = 0.0;
      for (int i = 0; i < verts.Size(); i++) { zc += mesh->GetVertex(verts[i])[2]; }
      zc /= verts.Size();
      mesh->SetBdrAttribute(b, (std::abs(zc) < 1.0) ? 1 : 5);
   }
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool on_plane = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         if (std::abs(mesh->GetVertex(fv[i])[1]) > 1e-6 * kL) { on_plane = false; }
      }
      if (on_plane) { tris.push_back({fv[0], fv[1], fv[2]}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   return mesh;
}

static BoundaryConfig BC()
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.absorbing_attrs = {5};
   return bc;
}

static void Centroid(const Mesh &m, int e, real_t c[3])
{
   Array<int> ev;
   m.GetElementVertices(e, ev);
   c[0] = c[1] = c[2] = 0.0;
   for (int a = 0; a < ev.Size(); a++)
      for (int d = 0; d < 3; d++) { c[d] += m.GetVertex(ev[a])[d] / ev.Size(); }
}

static real_t RelDiff(real_t a, real_t b) { return std::abs(a - b) / std::max(std::abs(b), real_t(1e-300)); }

/// Fault setup with the deck's background tensor (scaled by mu when asked)
/// and its slip-weakening parameters.
static NitscheFaultSetup DeckSetup(const spatial::SpatialFrictionConfig &cfg,
                                   const NitscheDeckMaterial &dm, bool scale_by_mu,
                                   const IFaultClosure &cl, spatial::SpatialFrictionResolver &res)
{
   NitscheFaultSetup fs;
   fs.closure = &cl;
   const auto &S = cfg.stress;
   const real_t s0[9] = {S.sigma_xx_pa, S.sigma_xy_pa, S.sigma_xz_pa, S.sigma_xy_pa, S.sigma_yy_pa,
                         S.sigma_yz_pa, S.sigma_xz_pa, S.sigma_yz_pa, S.sigma_zz_pa};
   std::array<real_t, 9> sv;
   for (int i = 0; i < 9; i++) { sv[i] = s0[i]; }
   fs.prestress = [sv](const real_t[3], real_t s[9]) { for (int i = 0; i < 9; i++) { s[i] = sv[i]; } };
   if (scale_by_mu) { fs.prestress = ScaleByShearModulus(fs.prestress, dm.mu_at, kMu0); }
   const spatial::SlipWeakeningBlock sw = *cfg.slip_weakening;
   fs.parameters = [sw, &res](const real_t x[3], FaultPointState &st) {
      Vector xc(3);
      for (int i = 0; i < 3; i++) { xc(i) = x[i]; }
      Array<int> at(1);
      at[0] = 3;
      spatial::SlipWeakeningPerDOFParams r = res.ResolveSlipWeakening(sw, xc, at);
      st.mu_s = r.mu_s(0);
      st.mu_d = r.mu_d(0);
      st.Dc = r.d_c(0);
      st.cohesion = r.cohesion(0);
   };
   return fs;
}

// ---------------------------------------------------------------- (a)
static void TestMaterial(ParMesh &pm)
{
   if (g_rank == 0) { std::cout << "\n=== (a) deck material per element, jumps ===\n"; }
   {
      const spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(kDeckTPV6);
      NitscheDeckMaterial dm = BuildNitscheDeckMaterial(cfg);
      NitscheMaterialSample sample(pm, dm.field);
      const auto &hs = cfg.material.halfspace;
      real_t err = 0.0;
      long long nfar = 0, nnear = 0;
      for (int e = 0; e < pm.GetNE(); e++)
      {
         real_t c[3];
         Centroid(pm, e, c);
         const bool near = c[1] < 0.0;   // n = (0, -1, 0): near = sign(-y) >= 0
         const real_t vp = near ? hs.vp_near : hs.vp_far, vs = near ? hs.vs_near : hs.vs_far,
                      rho = near ? hs.rho_near : hs.rho_far;
         const real_t mu = rho * vs * vs, lam = rho * (vp * vp - 2.0 * vs * vs);
         err = std::max({err, RelDiff(sample.MuGF()(e), mu), RelDiff(sample.LambdaGF()(e), lam),
                         RelDiff(sample.RhoGF()(e), rho)});
         (near ? nnear : nfar)++;
      }
      err = RMax(err);
      nfar = LSum(nfar);
      nnear = LSum(nnear);
      TEST_ASSERT(err < 1e-14 && nfar > 0 && nnear > 0,
                  "TPV6 deck: element samples = near side (Y < 0) / far side (Y > 0) values, max rel "
                  "diff " << err);
      TEST_ASSERT(!dm.mu_at && dm.jumps.has_plane && dm.jumps.depths.empty(),
                  "TPV6 deck: no point shear modulus on the fault; one jump plane");
      std::vector<long long> c0 = CountElementsAcrossJumps(pm, dm.jumps, 1e-6);
      MaterialJumps off = dm.jumps;
      off.x0[1] = 1000.0;
      std::vector<long long> c1 = CountElementsAcrossJumps(pm, off, 1e-6);
      const long long a0 = LSum(c0[0]), a1 = LSum(c1[0]);
      TEST_ASSERT(a0 == 0 && a1 > 0, "elements across the plane: " << a0
                  << " at Y = 0 (a mesh plane), " << a1 << " at Y = 1000 m");
   }
   {
      const spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(kDeckTPV31);
      NitscheDeckMaterial dm = BuildNitscheDeckMaterial(cfg);
      NitscheMaterialSample sample(pm, dm.field);
      real_t err = 0.0;
      for (int e = 0; e < pm.GetNE(); e++)
      {
         real_t c[3], lam, mu, rho;
         Centroid(pm, e, c);
         dm.depth->eval_at_xyz(c[0], c[1], c[2], lam, mu, rho);
         err = std::max({err, RelDiff(sample.MuGF()(e), mu), RelDiff(sample.LambdaGF()(e), lam),
                         RelDiff(sample.RhoGF()(e), rho)});
      }
      err = RMax(err);
      TEST_ASSERT(err < 1e-12, "TPV31 deck: element samples = profile at the centroid, max rel diff "
                  << err);
      const std::vector<real_t> want = {2400.0, 5000.0, 10000.0};
      TEST_ASSERT(dm.jumps.depths == want, "TPV31 deck: jumps at 2400, 5000, 10000 m (found "
                  << dm.jumps.depths.size() << ")");
      // spec table values just above / below each jump
      auto mu_at_depth = [&](real_t d) { real_t x[3] = {0.0, 0.0, -d}; return dm.mu_at(x); };
      const real_t m1 = mu_at_depth(1000.0), m2 = mu_at_depth(2400.0 + 1e-3),
                   m3 = mu_at_depth(5000.0 - 1e-3), m4 = mu_at_depth(7000.0), m5 = mu_at_depth(12000.0);
      TEST_ASSERT(RelDiff(m1, 2580.0 * 2250.0 * 2250.0) < 1e-12
                  && RelDiff(m2, 2600.0 * 2550.0 * 2550.0) < 1e-6
                  && RelDiff(m3, 2620.0 * 3050.0 * 3050.0) < 1e-6
                  && RelDiff(m4, 2720.0 * 3450.0 * 3450.0) < 1e-12
                  && RelDiff(m5, 3000.0 * 3800.0 * 3800.0) < 1e-12,
                  "TPV31 deck: mu = rho vs^2 of the description's table at 1, 2.4+, 5-, 7, 12 km");
      std::vector<long long> c0 = CountElementsAcrossJumps(pm, dm.jumps, 1e-6);
      MaterialJumps off = dm.jumps;
      off.depths = {3000.0};
      std::vector<long long> c1 = CountElementsAcrossJumps(pm, off, 1e-6);
      const long long a0 = LSum(c0[1]) + LSum(c0[2]) + LSum(c0[3]), a1 = LSum(c1[1]);
      TEST_ASSERT(a0 == 0 && a1 > 0, "elements across the jumps: none at 2400/5000/10000 m "
                  "(mesh planes), " << a1 << " across 3000 m");
   }
}

// ------------------------------------------------- (a) penalty and time step
static void TestPenaltyAndDt(ParMesh &pm, int p)
{
   if (g_rank == 0) { std::cout << "\n=== (a) bimaterial penalty and time step (p = " << p << ") ===\n"; }
   const spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(kDeckTPV6);
   NitscheDeckMaterial dm = BuildNitscheDeckMaterial(cfg);
   SlipWeakeningClosure cl;
   spatial::SpatialFrictionResolver res;
   NitscheFaultSetup fs = DeckSetup(cfg, dm, false, cl, res);
   NitscheWaveOperator op(pm, p, dm.field, BC(), {}, &fs);
   ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
   Vector sigF, spread;
   E.GetFaultFacePenalties(sigF, spread);
   const auto &hs = cfg.material.halfspace;
   auto moduli = [&](bool near, real_t &lam, real_t &mu) {
      const real_t vp = near ? hs.vp_near : hs.vp_far, vs = near ? hs.vs_near : hs.vs_far,
                   rho = near ? hs.rho_near : hs.rho_far;
      mu = rho * vs * vs;
      lam = rho * (vp * vp - 2.0 * vs * vs);
   };
   ConstantCoefficient one(1.0);
   DGElasticityIPCombinedIntegrator integ(one, one, 3, -1.0, 1.0);
   const ParFiniteElementSpace &fes = E.GetFESpace();
   const int n_int = E.GetFaultInteriorFaces().Size();
   const int nloc = n_int + E.GetFaultSharedFaces().Size();
   real_t err = 0.0, min_vs_near = 1e300;
   for (int fi = 0; fi < nloc; fi++)
   {
      FaceElementTransformations *Tr =
         (fi < n_int) ? pm.GetInteriorFaceTransformations(E.GetFaultInteriorFaces()[fi])
                      : pm.GetSharedFaceTransformations(E.GetFaultSharedFaces()[fi - n_int]);
      const IntegrationPoint &ipc = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetAllIntPoints(&ipc);
      Vector nor(3), c1(3), c2(3);
      CalcOrtho(Tr->Jacobian(), nor);
      const IntegrationPoint &ipe = Geometries.GetCenter(Geometry::TETRAHEDRON);
      Tr->Elem1->Transform(ipe, c1);
      Tr->Elem2->Transform(ipe, c2);
      real_t l1, m1, l2, m2, ln, mn;
      moduli(c1(1) < 0.0, l1, m1);
      moduli(c2(1) < 0.0, l2, m2);
      moduli(true, ln, mn);
      Tr->SetAllIntPoints(&ipc);
      const real_t w1 = Tr->Elem1->Weight(), w2 = Tr->Elem2->Weight();
      const FiniteElement &fe1 = *fes.GetFE(Tr->Elem1No);
      const FiniteElement &fe2 = (fi < n_int) ? *fes.GetFE(Tr->Elem2No) : fe1;
      const real_t want = integ.GetPenalty(fe1, fe2, w1, w2, l1, m1, l2, m2, nor.Norml2());
      const real_t near_only = integ.GetPenalty(fe1, fe2, w1, w2, ln, mn, ln, mn, nor.Norml2());
      err = std::max(err, RelDiff(sigF(fi), want));
      min_vs_near = std::min(min_vs_near, std::abs(sigF(fi) / near_only - 1.0));
   }
   err = RMax(err);
   min_vs_near = -RMax(-min_vs_near);
   TEST_ASSERT(err < 1e-12 && min_vs_near > 0.1, "TPV6 fault faces: sigma_F = the IP penalty with "
               "each side's moduli (max rel diff " << err << "), at least " << 100 * min_vs_near
               << " % off the near-side-only value");
   // Time step on a copy whose far side (Y > 0) is squashed to 1/20 of its
   // width, so the slow far side sets dt and a single c_p would show.
   ParMesh sq(pm);
   sq.Transform([](const Vector &x, Vector &y) { y = x; if (x(1) > 0.0) { y(1) *= 0.05; } });
   NitscheMaterialSample sample(sq, dm.field);
   const real_t cfl = 0.5;
   const real_t dt = ComputeMaxDt(sample, sq, p, cfl);
   real_t dt_side[2] = {1e300, 1e300};   // far, near
   for (int e = 0; e < sq.GetNE(); e++)
   {
      real_t c[3], lam, mu;
      Centroid(sq, e, c);
      const bool near = c[1] < 0.0;
      moduli(near, lam, mu);
      const real_t rho = near ? hs.rho_near : hs.rho_far;
      dt_side[near] = std::min(dt_side[near], cfl * NitscheElementLength(sq, e) /
                               (p * p * std::sqrt((lam + 2.0 * mu) / rho)));
   }
   dt_side[0] = -RMax(-dt_side[0]);
   dt_side[1] = -RMax(-dt_side[1]);
   const real_t dt_hand = std::min(dt_side[0], dt_side[1]);
   TEST_ASSERT(RelDiff(dt, dt_hand) < 1e-14 && dt_side[0] < 0.9 * dt_side[1],
               "time step " << dt << " s = min over elements of cfl h_e / (p^2 c_p,e), each "
               "element's c_p (hand " << dt_hand << " s; far side " << dt_side[0]
               << " s sets it, near side " << dt_side[1] << " s)");
}

// ---------------------------------------------------------------- (b), (d)
static void TestPrestressAndCohesion(ParMesh &pm, int p)
{
   if (g_rank == 0) { std::cout << "\n=== (b) mu-scaled prestress, (d) cohesion (p = " << p << ") ===\n"; }
   const spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(kDeckTPV31);
   NitscheDeckMaterial dm = BuildNitscheDeckMaterial(cfg);
   SlipWeakeningClosure cl;
   spatial::SpatialFrictionResolver res;
   NitscheFaultSetup fs = DeckSetup(cfg, dm, true, cl, res);
   NitscheWaveOperator op(pm, p, dm.field, BC(), {}, &fs);
   const int no = op.NumOwnedFaultDOFs();
   const Vector &S0 = op.PrestressTraction(), &X = op.FaultPositions(), &F = op.FaultFrames();
   real_t err_c = 0.0, err_g = 0.0, err_coh = 0.0;
   long long n_c = 0, n_g = 0, n_if = 0, n_if_bad = 0;
   long long have[3] = {0, 0, 0};
   const real_t tops[3] = {2400.0, 5000.0, 10000.0};
   for (int j = 0; j < no; j++)
   {
      const real_t *x = &X(3 * j), *n = &F(9 * j), *s = &F(9 * j + 6);
      const real_t *T = &S0(3 * j);
      const real_t tau2 = -(T[0] * s[0] + T[1] * s[1] + T[2] * s[2]);
      const real_t sig = -(T[0] * n[0] + T[1] * n[1] + T[2] * n[2]);
      const real_t depth = -x[2];
      bool on_if = false;
      for (real_t d : tops) { on_if = on_if || std::abs(depth - d) < 1.0; }
      const real_t f = dm.mu_at(x) / kMu0;
      if (on_if)
      {
         // the DOF belongs to a face above or below: either side's value
         // (round-off next to two constant layers; next to the gradient
         // layer the face projection of a varying mu, as below)
         n_if++;
         real_t xa[3] = {x[0], x[1], x[2] + 1.0}, xb[3] = {x[0], x[1], x[2] - 1.0};
         const real_t fa = dm.mu_at(xa) / kMu0, fb = dm.mu_at(xb) / kMu0;
         const real_t tol = (std::abs(depth - 10000.0) < 1.0) ? 1e-11 : 0.05;
         const bool ok = std::min(RelDiff(tau2, 30e6 * fa), RelDiff(tau2, 30e6 * fb)) < tol
                         && std::min(RelDiff(sig, 60e6 * fa), RelDiff(sig, 60e6 * fb)) < tol;
         if (!ok) { n_if_bad++; }
      }
      else if (depth > 2400.0 && depth < 5000.0)
      {
         // the gradient layer: mu varies on the face, the projection averages it
         n_g++;
         err_g = std::max({err_g, RelDiff(tau2, 30e6 * f), RelDiff(sig, 60e6 * f)});
      }
      else
      {
         n_c++;
         err_c = std::max({err_c, RelDiff(tau2, 30e6 * f), RelDiff(sig, 60e6 * f)});
      }
      const real_t coh = (depth <= 2400.0) ? 425.0 * (2400.0 - depth) : 0.0;
      err_coh = std::max(err_coh, std::abs(op.FaultParams(j).cohesion - coh));
      if (std::abs(depth) < 1.0) { have[0]++; }
      if (std::abs(depth - 1200.0) < 1.0) { have[1]++; }
      if (std::abs(depth - 2400.0) < 1.0) { have[2]++; }
   }
   err_c = RMax(err_c);
   err_g = RMax(err_g);
   err_coh = RMax(err_coh);
   n_c = LSum(n_c);
   n_g = LSum(n_g);
   n_if = LSum(n_if);
   n_if_bad = LSum(n_if_bad);
   for (int k = 0; k < 3; k++) { have[k] = LSum(have[k]); }
   TEST_ASSERT(n_c > 0 && err_c < 1e-11, "constant layers: tau2 = 30 MPa mu/mu_0, sigma_n = "
               "60 MPa mu/mu_0 at " << n_c << " DOFs, max rel diff " << err_c);
   TEST_ASSERT(n_g > 0 && err_g < 0.05, "gradient layer (2.4-5 km): within 5 % of the point "
               "value at " << n_g << " DOFs (face projection of a varying mu), max rel diff "
               << err_g);
   TEST_ASSERT(n_if > 0 && n_if_bad == 0, "layer interfaces: each DOF carries its face's layer "
               "value (" << n_if << " DOFs, " << n_if_bad << " off)");
   TEST_ASSERT(err_coh < 1e-6 && have[0] > 0 && have[1] > 0 && have[2] > 0,
               "cohesion 425 Pa/m x (2400 m - depth) above 2400 m, 0 below, at every DOF "
               "(DOFs at 0, 1.2, 2.4 km: " << have[0] << ", " << have[1] << ", " << have[2]
               << "), max abs diff " << err_coh << " Pa");
   Vector xc(3);
   xc(0) = 0.0; xc(1) = 0.0; xc(2) = -3000.0;
   Array<int> at(1);
   at[0] = 3;
   const real_t c3 = res.ResolveSlipWeakening(*cfg.slip_weakening, xc, at).cohesion(0);
   xc(2) = 0.0;
   const real_t c0 = res.ResolveSlipWeakening(*cfg.slip_weakening, xc, at).cohesion(0);
   TEST_ASSERT(c3 == 0.0 && std::abs(c0 - 1.02e6) < 1e-6, "cohesion 1.02 MPa at the surface, 0 at "
               "3 km (" << c0 << ", " << c3 << ")");
}

// ---------------------------------------------------------------- (c)
static void TestNucleation()
{
   if (g_rank == 0) { std::cout << "\n=== (c) instantaneous circular nucleation ===\n"; }
   const spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(kDeckTPV31);
   NitscheDeckMaterial dm = BuildNitscheDeckMaterial(cfg);
   TEST_ASSERT(cfg.nucleation.enabled && cfg.nucleation.kind
               == spatial::NucleationKind::InstantaneousOverstressCircular,
               "TPV31 deck: instantaneous_overstress_circular");
   const auto &spec = cfg.nucleation.instantaneous_circular;
   const real_t n[3] = {0.0, -1.0, 0.0}, s[3] = {1.0, 0.0, 0.0}, d[3] = {0.0, 0.0, -1.0};
   StressFn nuc = InstantaneousCircularOverstress(spec, s, n, d, dm.mu_at);
   const real_t r[5] = {0.0, 1400.0, 1700.0, 2000.0, 2100.0}, fac[5] = {1.0, 1.0, 0.5, 0.0, 0.0};
   real_t err = 0.0, other = 0.0;
   for (int dir = 0; dir < 2; dir++)
   {
      for (int k = 0; k < 5; k++)
      {
         // along strike, then down dip from the hypocentre (0, 0, -7500)
         const real_t x[3] = {dir == 0 ? r[k] : 0.0, 0.0, dir == 0 ? -7500.0 : -7500.0 - r[k]};
         real_t S[9];
         nuc(x, S);
         real_t T[3];
         for (int i = 0; i < 3; i++) { T[i] = S[3 * i] * n[0] + S[3 * i + 1] * n[1] + S[3 * i + 2] * n[2]; }
         const real_t tau2 = -(T[0] * s[0] + T[1] * s[1] + T[2] * s[2]);
         const real_t tau1 = -(T[0] * d[0] + T[1] * d[1] + T[2] * d[2]);
         const real_t sn = -(T[0] * n[0] + T[1] * n[1] + T[2] * n[2]);
         const real_t want = 4.95e6 * fac[k] * dm.mu_at(x) / kMu0;
         err = std::max(err, std::abs(tau2 - want) / 4.95e6);
         other = std::max({other, std::abs(tau1), std::abs(sn)});
      }
   }
   TEST_ASSERT(err < 1e-14 && other == 0.0, "tau2 = 4.95 MPa mu/mu_0 x (1, 1, 0.5, 0, 0) at r = 0, "
               "1400, 1700, 2000, 2100 m (strike and dip), max diff " << err * 4.95 << " MPa; "
               "dip and normal components " << other);
}

/// (c) planar-fault check (driver, both nucleation kinds): sin of the largest
/// angle between a fault face normal and the reference normal.
static void TestPlanarCheck(ParMesh &pm)
{
   if (g_rank == 0) { std::cout << "\n=== (c) planar-fault check for the nucleation ===\n"; }
   const spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(kDeckTPV6);
   NitscheDeckMaterial dm = BuildNitscheDeckMaterial(cfg);
   SlipWeakeningClosure cl;
   spatial::SpatialFrictionResolver res;
   const real_t n_ref[3] = {0.0, -1.0, 0.0};
   NitscheFaultSetup fs = DeckSetup(cfg, dm, false, cl, res);
   NitscheWaveOperator op(pm, 1, dm.field, BC(), {}, &fs);
   const real_t flat = RMax(MaxFaultNormalDeviation(op.Orientation(), n_ref));
   // bend the fault plane into Y = 1e-6 X^2 (slope up to 0.06 at |X| = 30 km)
   ParMesh bent(pm);
   bent.Transform([](const Vector &x, Vector &y) { y = x; y(1) += 1e-6 * x(0) * x(0); });
   NitscheFaultSetup fsb = DeckSetup(cfg, dm, false, cl, res);
   NitscheWaveOperator opb(bent, 1, dm.field, BC(), {}, &fsb);
   const real_t curved = RMax(MaxFaultNormalDeviation(opb.Orientation(), n_ref));
   TEST_ASSERT(flat < 1e-12 && curved > 1e-6, "largest |n x n_ref| over the fault faces: " << flat
               << " on the planar fault (accepted), " << curved << " on the bent copy (refused, "
               "limit 1e-6)");
}

// ---------------------------------------------------------------- (e)
/// u = A_side + B_side x, the side from the element centroid (far: Y > 0).
class SidedLinear : public VectorCoefficient
{
public:
   SidedLinear(const Mesh &m, const real_t Af[3], const real_t Bf[9], const real_t An[3],
               const real_t Bn[9])
      : VectorCoefficient(3)
   {
      far_.resize(m.GetNE());
      for (int e = 0; e < m.GetNE(); e++)
      {
         real_t c[3];
         Centroid(m, e, c);
         far_[e] = c[1] > 0.0;
      }
      for (int i = 0; i < 3; i++) { Af_[i] = Af[i]; An_[i] = An[i]; }
      for (int i = 0; i < 9; i++) { Bf_[i] = Bf[i]; Bn_[i] = Bn[i]; }
   }
   static void At(const real_t A[3], const real_t B[9], const real_t x[3], real_t u[3])
   {
      for (int i = 0; i < 3; i++) { u[i] = A[i] + B[3 * i] * x[0] + B[3 * i + 1] * x[1] + B[3 * i + 2] * x[2]; }
   }
   void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
   {
      Vector x(3);
      T.Transform(ip, x);
      V.SetSize(3);
      const bool f = far_[T.ElementNo];
      At(f ? Af_ : An_, f ? Bf_ : Bn_, x.GetData(), V.GetData());
   }
   real_t Af_[3], Bf_[9], An_[3], Bn_[9];
private:
   std::vector<char> far_;
};

static void TestPerSide(ParMesh &pm, int p)
{
   if (g_rank == 0) { std::cout << "\n=== (e) per-side station writer (p = " << p << ") ===\n"; }
   const spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(kDeckTPV6);
   NitscheDeckMaterial dm = BuildNitscheDeckMaterial(cfg);
   SlipWeakeningClosure cl;
   spatial::SpatialFrictionResolver res;
   NitscheFaultSetup fs = DeckSetup(cfg, dm, false, cl, res);
   NitscheWaveOperator op(pm, p, dm.field, BC(), {}, &fs);
   ParFiniteElementSpace &fes = op.Elasticity().GetFESpace();
   const real_t Af[3] = {0.7, -0.2, 0.31}, An[3] = {-0.4, 0.15, -0.05};
   const real_t Bf[9] = {1e-5, 2e-6, -3e-6, 4e-6, -1e-5, 5e-6, -2e-6, 3e-6, 7e-6};
   const real_t Bn[9] = {-3e-6, 1e-6, 2e-6, 6e-6, 2e-6, -4e-6, 1e-6, -5e-6, 3e-6};
   const real_t Vf[3] = {2.0, 0.3, -0.1}, Vn[3] = {-1.0, -0.2, 0.4};
   const real_t Wf[9] = {3e-5, 0, 0, 0, 1e-5, 0, 0, 0, -2e-5}, Wn[9] = {0, 1e-5, 0, 0, 0, 2e-5, 1e-5, 0, 0};
   SidedLinear cu(pm, Af, Bf, An, Bn), cv(pm, Vf, Wf, Vn, Wn);
   ParGridFunction ug(&fes), vg(&fes);
   ug.ProjectCoefficient(cu);
   vg.ProjectCoefficient(cv);
   Vector U(ug.GetData(), fes.GetVSize()), Ud(vg.GetData(), fes.GetVSize());
   std::vector<FaultStation> st = {{"surface_x0", {0.0, 0.0, 0.0}},
                                   {"x7.5_z5", {7500.0, 0.0, -5000.0}},
                                   {"xm13_z2.4", {-13000.0, 0.0, -2400.0}},
                                   {"x3_z12", {3000.0, 0.0, -12000.0}}};
   const std::vector<FaultStationSnap> snap = SnapFaultStations(st, op.FaultPositions(), MPI_COMM_WORLD);
   const std::string dir = "test_nitsche_tpv_setup_out";
   if (g_rank == 0) { mkdir(dir.c_str(), 0755); }
   MPI_Barrier(MPI_COMM_WORLD);
   const real_t up[3] = {0.0, 0.0, 1.0};
   FaultSideStationWriter w(op, st, snap, up, dir, "tps");
   std::vector<DOFData> dd(op.NumOwnedFaultDOFs());
   w.WriteStep(0.0, U, Ud, dd);
   w.Close();
   bool copies_ok = true;
   for (size_t s = 0; s < st.size(); s++)
   {
      copies_ok = copies_ok && w.NumCopies(s, 0) > 0 && w.NumCopies(s, 0) == w.NumCopies(s, 1);
   }
   TEST_ASSERT(copies_ok, "every station has the same number (> 0) of far and near copies");
   if (g_rank != 0) { return; }
   real_t err = 0.0, ferr = 0.0;
   for (size_t s = 0; s < st.size(); s++)
   {
      for (int side = 0; side < 2; side++)
      {
         real_t u[3], v[3];
         const bool far = (side == FaultSideStationWriter::kFar);
         SidedLinear::At(far ? Af : An, far ? Bf : Bn, snap[s].node, u);
         SidedLinear::At(far ? Vf : Vn, far ? Wf : Wn, snap[s].node, v);
         // h = +X, v = -Z (down dip), n = +Y (near side to far side)
         const real_t want[6] = {u[0], -u[2], u[1], v[0], -v[2], v[1]};
         const auto &got = w.LastValues(s, side);
         for (int c = 0; c < 6; c++) { err = std::max(err, std::abs(got[c] - want[c])); }
         std::ifstream f(dir + "/tps_station_" + st[s].name + (far ? "_farside.dat" : "_nearside.dat"));
         std::string line, last;
         while (std::getline(f, line)) { if (!line.empty() && line[0] != '#') { last = line; } }
         std::istringstream is(last);
         real_t col[10];
         for (int c = 0; c < 10; c++) { is >> col[c]; }
         const real_t fw[6] = {col[1], col[4], col[7], col[2], col[5], col[8]};
         for (int c = 0; c < 6; c++) { ferr = std::max(ferr, std::abs(fw[c] - want[c])); }
      }
   }
   std::cout << "  (stations snapped to";
   for (size_t s = 0; s < st.size(); s++)
   {
      std::cout << " (" << snap[s].node[0] << ", " << snap[s].node[1] << ", " << snap[s].node[2] << ")";
   }
   std::cout << ")\n";
   TEST_ASSERT(err < 1e-12, "far (Y > 0) and near (Y < 0) displacement and velocity, h = +X, "
               "v = -Z, n = +Y: max abs diff " << err);
   TEST_ASSERT(ferr < 1e-9, "the per-side files hold the same values (10 significant digits): max "
               "abs diff " << ferr);
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   if (g_rank == 0)
   {
      std::cout << "=== Nitsche TPV6/TPV7/TPV31 setup on " << Mpi::WorldSize() << " rank(s) ===\n";
   }
   std::unique_ptr<Mesh> m(LayeredFaultBox(4, 4));
   // on 2 ranks the partition is the fault plane itself (Y < 0 on rank 0,
   // Y > 0 on rank 1), so every fault face is a shared face; otherwise METIS
   std::vector<int> part;
   if (Mpi::WorldSize() == 2)
   {
      part.resize(m->GetNE());
      for (int e = 0; e < m->GetNE(); e++)
      {
         real_t c[3];
         Centroid(*m, e, c);
         part[e] = (c[1] > 0.0) ? 1 : 0;
      }
   }
   ParMesh pm(MPI_COMM_WORLD, *m, part.empty() ? nullptr : part.data());
   TestMaterial(pm);
   TestPenaltyAndDt(pm, 1);
   TestPenaltyAndDt(pm, 2);
   TestPrestressAndCohesion(pm, 1);
   TestPrestressAndCohesion(pm, 2);
   TestNucleation();
   TestPlanarCheck(pm);
   TestPerSide(pm, 1);
   TestPerSide(pm, 2);
   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests << " passed, " << num_failed
                << " failed ===\n";
   }
   return (num_failed == 0) ? 0 : 1;
}
