// Part B / Phase B1: per-side fault impedance assignment via the eps-offset rule
// (BimaterialWaveOperator::AssignFaultSidePerMaterialImpedances).
//
// Fixture: a row of 6 unit hexes along x; the interior face at x=3 is the FAULT
// (between elem 2 [centroid 2.5] and elem 3 [centroid 3.5]).  Two materials:
//   - ACROSS-FAULT (halfspace): mu/lambda/rho step at x=3 -> a bi-material fault.
//     After assignment each fault DOF must have {Zp_plus,Zp_minus} == {Zp_far,Zp_near}
//     (a real contrast), eta the harmonic mean.
//   - SYMMETRIC (constant): -> Zp_plus == Zp_minus byte-exact (no spurious contrast;
//     the eps-offset along the fault normal does not change the constant value).
//
// This verifies the eps-offset disambiguates the across-fault step (evaluating
// exactly at the face is ambiguous) AND keeps a fault-symmetric material equal.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
#define TEST_ASSERT(cond, msg) do {                                          \
   num_tests++;                                                              \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; }     \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: "       \
                                  << msg << "\n"; }                          \
} while (0)

namespace
{
// Row of `nx` unit hexes along x; interior face at x==fault_x tagged FAULT (attr 3).
Mesh BuildHexRowFaultMesh(int nx, int fault_x)
{
   Mesh mesh(3, (nx + 1) * 4, nx, 0);
   auto vid = [](int ix, int iy, int iz) { return ix * 4 + iy * 2 + iz; };
   for (int ix = 0; ix <= nx; ++ix)
      for (int iy = 0; iy < 2; ++iy)
         for (int iz = 0; iz < 2; ++iz)
         { real_t v[3] = {(real_t)ix, (real_t)iy, (real_t)iz}; mesh.AddVertex(v); }
   for (int e = 0; e < nx; ++e)
   {
      const int v[8] = { vid(e,0,0), vid(e+1,0,0), vid(e+1,1,0), vid(e,1,0),
                         vid(e,0,1), vid(e+1,0,1), vid(e+1,1,1), vid(e,1,1) };
      mesh.AddHex(v, 1);
   }
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0;
      for (int i = 0; i < fv.Size(); ++i) { cx += mesh.GetVertex(fv[i])[0]; }
      cx /= fv.Size();
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cx - fault_x) < 1e-9)
         { mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 3); }
         continue;
      }
      mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Lame/density from (vp, vs, rho).
struct Mat { real_t vp, vs, rho; };
real_t Zp(const Mat &m) { return m.rho * m.vp; }
real_t Zs(const Mat &m) { return m.rho * m.vs; }
void lame(const Mat &m, real_t &lam, real_t &mu, real_t &rho)
{ mu = m.rho * m.vs * m.vs; lam = m.rho * m.vp * m.vp - 2.0 * mu; rho = m.rho; }

const Mat kFar  = {3750.0, 2165.0, 2225.0};   // x < 3
const Mat kNear = {6000.0, 3464.0, 2670.0};   // x > 3
const real_t kFaultX = 3.0;

bool near_rel(real_t a, real_t b, real_t rel)
{ return std::abs(a - b) <= rel * std::max(1.0, std::max(std::abs(a), std::abs(b))); }

// Size + seed the fault DOFData like the driver (single-material default), then
// run the per-side assignment.
int RunAssign(Mesh &mesh, BimaterialWaveOperator<Mesh> &op,
              std::vector<DOFData> &dof, const Mat &seed)
{
   const int qdeg = op.FaultFaceQuadDegree();
   const int nfaces = op.GetFaultInteriorFaces().Size();
   int nbf = 0;
   if (nfaces > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(op.GetFaultInteriorFaces()[0]);
      nbf = IntRules.Get(ftr->GetGeometryType(), qdeg).GetNPoints();
   }
   const int ntot = nfaces * nbf;
   dof.assign(ntot, DOFData());
   real_t lamS, muS, rhoS; lame(seed, lamS, muS, rhoS);
   const real_t ZpS = std::sqrt((lamS + 2.0 * muS) * rhoS);
   const real_t ZsS = std::sqrt(muS * rhoS);
   for (auto &d : dof)
   {
      d.Zp_plus = d.Zp_minus = ZpS;
      d.Zs_plus = d.Zs_minus = ZsS;
      d.eta_p = 0.5 * ZpS; d.eta_s = 0.5 * ZsS;
   }
   op.SetFaultDOFData(&dof, nbf);
   op.AssignFaultSidePerMaterialImpedances(dof);
   return ntot;
}
} // anonymous namespace

int main()
{
   std::cout << "\n=== Part B / B1: per-side fault impedance (eps-offset) ===\n";

   Mesh mesh = BuildHexRowFaultMesh(/*nx=*/6, /*fault_x=*/3);
   BoundaryConfig bc; bc.natural_attrs = {1}; bc.fault_attr = 3; bc.absorbing_attrs = {};
   const int order = 1;

   const real_t ZpFar = Zp(kFar), ZpNear = Zp(kNear);
   const real_t ZsFar = Zs(kFar), ZsNear = Zs(kNear);

   // ---- ACROSS-FAULT (halfspace) ----
   std::cout << "\n-- across-fault halfspace: bi-material contrast --\n";
   {
      auto lam_fn = [](const Vector &x){ real_t l,m,r; lame(x(0)<kFaultX?kFar:kNear,l,m,r); return l; };
      auto mu_fn  = [](const Vector &x){ real_t l,m,r; lame(x(0)<kFaultX?kFar:kNear,l,m,r); return m; };
      auto rho_fn = [](const Vector &x){ return x(0)<kFaultX ? kFar.rho : kNear.rho; };
      FunctionCoefficient lam_c(lam_fn), mu_c(mu_fn), rho_c(rho_fn);
      // The operator stores a POINTER to the MaterialField (material_ = &material),
      // and AssignFault runs after the ctor, so keep it alive in a named local.
      MaterialField mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
      BimaterialWaveOperator<Mesh> op(mesh, order, mat, bc);

      std::vector<DOFData> dof;
      const int ntot = RunAssign(mesh, op, dof, kNear);   // seed = single (near) material
      TEST_ASSERT(ntot > 0, "fault DOFs exist (ntot > 0)");

      bool all_contrast = true, set_ok = true, eta_ok = true, finite = true;
      for (const auto &d : dof)
      {
         const real_t zpmin = std::min(d.Zp_plus, d.Zp_minus);
         const real_t zpmax = std::max(d.Zp_plus, d.Zp_minus);
         const real_t zsmin = std::min(d.Zs_plus, d.Zs_minus);
         const real_t zsmax = std::max(d.Zs_plus, d.Zs_minus);
         all_contrast = all_contrast && (zpmax - zpmin > 1.0);
         set_ok = set_ok && near_rel(zpmin, ZpFar, 1e-12) && near_rel(zpmax, ZpNear, 1e-12)
                         && near_rel(zsmin, ZsFar, 1e-12) && near_rel(zsmax, ZsNear, 1e-12);
         const real_t etap = d.Zp_plus * d.Zp_minus / (d.Zp_plus + d.Zp_minus);
         const real_t etas = d.Zs_plus * d.Zs_minus / (d.Zs_plus + d.Zs_minus);
         eta_ok = eta_ok && near_rel(d.eta_p, etap, 1e-12) && near_rel(d.eta_s, etas, 1e-12);
         finite = finite && std::isfinite(d.Zp_plus) && std::isfinite(d.Zp_minus);
      }
      TEST_ASSERT(finite, "all per-side impedances finite");
      TEST_ASSERT(all_contrast, "Zp_plus != Zp_minus at every fault DOF (real contrast)");
      TEST_ASSERT(set_ok,
                  "{Zp_plus,Zp_minus} == {Zp_far,Zp_near} and {Zs...} (correct per-side "
                  "material; eps-offset disambiguated the step)");
      TEST_ASSERT(eta_ok, "eta_p/eta_s == harmonic mean of the per-side impedances");
   }

   // ---- SYMMETRIC (constant) ----
   std::cout << "\n-- fault-symmetric (constant) material: byte-exact equal --\n";
   {
      real_t l,m,r; lame(kNear, l, m, r);
      ConstantCoefficient lam_c(l), mu_c(m), rho_c(r);
      MaterialField mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
      BimaterialWaveOperator<Mesh> op(mesh, order, mat, bc);
      std::vector<DOFData> dof;
      const int ntot = RunAssign(mesh, op, dof, kFar);   // seed with the OTHER material...
      bool equal = (ntot > 0);
      for (const auto &d : dof)
      {
         equal = equal && (d.Zp_plus == d.Zp_minus) && (d.Zs_plus == d.Zs_minus);
         // ...and the assignment overwrote the seed with the constant material:
         equal = equal && near_rel(d.Zp_plus, ZpNear, 1e-12);
      }
      TEST_ASSERT(equal,
                  "constant material => Zp_plus == Zp_minus (byte-exact, no spurious "
                  "contrast) and overwritten to the constant value");
   }

   // ---- B3: MakeHalfspaceAcrossFaultMaterial factory ----
   std::cout << "\n-- B3 factory: MakeHalfspaceAcrossFaultMaterial --\n";
   {
      // n = (0,-1,0): near side = sign((x-x0).n) = -y >= 0 => y <= 0.
      const real_t x0[3] = {0.0, 0.0, 0.0};
      const real_t nrm[3] = {0.0, -1.0, 0.0};
      auto w = MakeHalfspaceAcrossFaultMaterial(
         kNear.vp, kNear.vs, kNear.rho, kFar.vp, kFar.vs, kFar.rho, x0, nrm);
      TEST_ASSERT(w && w->field.mode == MaterialField::Mode::Coefficient,
                  "factory returns a Mode::Coefficient MaterialField");
      real_t lamN, muN, rhoN, lamF, muF, rhoF;
      lame(kNear, lamN, muN, rhoN);
      lame(kFar,  lamF, muF, rhoF);
      real_t l, m, r;
      w->eval_at_xyz(0.0, -1.0, 0.0, l, m, r);   // y<0 -> NEAR
      const bool near_ok = near_rel(l, lamN, 1e-12) && near_rel(m, muN, 1e-12)
                        && near_rel(r, rhoN, 1e-12);
      w->eval_at_xyz(0.0, +1.0, 0.0, l, m, r);   // y>0 -> FAR
      const bool far_ok = near_rel(l, lamF, 1e-12) && near_rel(m, muF, 1e-12)
                       && near_rel(r, rhoF, 1e-12);
      TEST_ASSERT(near_ok, "near side (sign((x-x0).n)>=0) -> near (vp,vs,rho)");
      TEST_ASSERT(far_ok,  "far side  (sign((x-x0).n)<0)  -> far  (vp,vs,rho)");
   }

   // ---- B2 wiring (R-001): operator SetFaultFlux affirms (or not) per-side-A ----
   std::cout << "\n-- B2 wiring: SetFaultFlux affirms per-side-A on the matrix op only --\n";
   {
      FunctionCoefficient lc([](const Vector&){ return 3.0e10; });
      FunctionCoefficient mc([](const Vector&){ return 3.0e10; });
      ConstantCoefficient  rc(2670.0);
      MaterialField mat = MaterialField::MakeCoefficient(&lc, &mc, &rc);  // outlives opb
      BimaterialWaveOperator<Mesh> opb(mesh, order, mat, bc);
      FaultFaceFlux ffb(2670.0, 6000.0, 3464.0);
      opb.SetFaultFlux(&ffb);
      TEST_ASSERT(ffb.GetPerSideFluxApplied(),
                  "BimaterialWaveOperator::SetFaultFlux affirms per-side-A "
                  "(virtual override sets the flag => matrix bi-material fault runs)");

      WaveOperator<Mesh> ops(mesh, order, 3.0e10, 3.0e10, 2670.0, bc);
      FaultFaceFlux ffs(2670.0, 6000.0, 3464.0);
      ops.SetFaultFlux(&ffs);
      TEST_ASSERT(!ffs.GetPerSideFluxApplied(),
                  "scalar WaveOperator::SetFaultFlux leaves the flag false "
                  "(scalar-path bi-material guard stays armed)");
   }

   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, " << num_failed
             << " failed out of " << num_tests << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
