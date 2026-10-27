// Nitsche dynamic rupture: Prakash-Clifton regularization of the normal
// stress the friction law sees (NitscheWaveOperator::SetPrakashClifton).
//
//   PC-b  exact update: one owned fault DOF, a constant slip rate V and a
//         step in sigma_eff from sigma*^0 to sigma_1; after k box-9 closes
//         sigma* = sigma_1 + (sigma*^0 - sigma_1) exp(-(V + v_star) k dt / L)
//         to round-off.
//   PC-c  L -> tiny (exp underflows to 0): (d, V, psi, sigma_eff, X) after
//         N steps of TPV102 with nucleation are bitwise those of PC off, and
//         sigma* equals sigma_eff bitwise.
//   PC-d  equilibrium without nucleation: |V| stays V_ini and sigma* stays
//         the initial 120 MPa (v_star = 1 m/s, L = 0.01 m).
//   PC-s  sanity on a JITTERED box (interior vertices off Y = 0 moved by up
//         to 15% of the cell size, so the Y -> -Y mirror symmetry is broken
//         and slip perturbs sigma_eff): with v_star = 1 m/s, L = 0.01 m and
//         nucleation the run differs from PC off and stays finite, and sigma*
//         lies within the range sigma_eff took (a convex average).
//   PC-a  abort: Prakash-Clifton on but sigma* never seeded (--abort-unseeded,
//         must exit non-zero with its message).
//
// Fixture: the TPV102 inline box of test_nitsche_tpv102_inline.cpp
// (24 tets per hex, fault on Y = 0), p = 1 and 2.
//
// Usage: mpirun -np 1|2|4 ./seas_test_nitsche_pc [--abort-unseeded]
//
// Scope (Chunhui Zhao, 2026-10-14).

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../friction/rate_state_closure.hpp"
#include "../../config/tpv102_params.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
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
static const real_t kCuse[3] = {0.0, 0.057, 0.11};   // 0.5 C_eig (P3)
static const real_t kVStar = 1.0, kLength = 0.01;     // SeisSol TPV6/7

static int ISum(int v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }

// ---------------------------------------------------------------- fixture
// (as test_nitsche_tpv102_inline.cpp: TPV102Box, TPV102BC, TPV102Material,
// TPV102Friction, TPV102Setup, TPV102InitialState)

static Mesh *TPV102Box(int nx, int ny, int nz, bool jitter = false)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2 * kL, 2 * kL, kL));
   const real_t hx = 2 * kL / nx, hy = 2 * kL / ny, hz = kL / nz;
   unsigned seed = 12345u;
   auto rnd = [&seed]() {   // deterministic in [-1, 1], same on every rank
      seed = seed * 1664525u + 1013904223u;
      return 2.0 * (seed >> 8) / 16777216.0 - 1.0;
   };
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
      const real_t r0 = rnd(), r1 = rnd(), r2 = rnd();
      // jitter only vertices strictly inside the box and off the fault plane
      const real_t e = 1e-6 * kL;
      if (jitter && std::abs(c[1]) > e && std::abs(std::abs(c[0]) - kL) > e
          && std::abs(std::abs(c[1]) - kL) > e && std::abs(c[2]) > e
          && std::abs(c[2] + kL) > e)
      {
         c[0] += 0.15 * hx * r0;
         c[1] += 0.15 * hy * r1;
         c[2] += 0.15 * hz * r2;
      }
   }
   const real_t tol = kL / nz * 0.1;
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      real_t zc = 0.0;
      for (int i = 0; i < verts.Size(); i++) { zc += mesh->GetVertex(verts[i])[2]; }
      zc /= verts.Size();
      mesh->SetBdrAttribute(b, (std::abs(zc) < tol) ? 1 : 5);
   }
   const real_t ytol = 1e-6 * kL;
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool on_plane = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         if (std::abs(mesh->GetVertex(fv[i])[1]) > ytol) { on_plane = false; }
      }
      if (on_plane) { tris.push_back({fv[0], fv[1], fv[2]}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   return mesh;
}

static BoundaryConfig TPV102BC()
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.absorbing_attrs = {5};
   return bc;
}

static MaterialField TPV102Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu,
                                      TPV102Params::rho);
}

struct TPV102Friction
{
   DieterichRuinaFriction fr{{TPV102Params::V0, TPV102Params::f0, TPV102Params::b,
                              TPV102Params::Dc}};
   RateStateClosure cl{fr, RateStateClosure::StateLaw::Aging};
};

static NitscheFaultSetup TPV102Setup(const IFaultClosure &cl, bool nucleation)
{
   NitscheFaultSetup fs;
   fs.closure = &cl;
   fs.prestress = [](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = TPV102Params::tau_ini;
      s[4] = -TPV102Params::sigma_n;
   };
   if (nucleation)
   {
      fs.nucleation_spatial = [](const real_t x[3], real_t s[9]) {
         for (int i = 0; i < 9; i++) { s[i] = 0.0; }
         const real_t dx = x[0] - TPV102Params::hypo_along_strike;
         const real_t dz = -x[2] - TPV102Params::hypo_down_dip;
         const real_t v = TPV102Params::nuc_dtau * NucleationSpatial(std::sqrt(dx * dx + dz * dz));
         s[1] = s[3] = v;
      };
      fs.nucleation_temporal = [](real_t t) { return NucleationTemporal(t); };
   }
   fs.parameters = [](const real_t x[3], FaultPointState &st) {
      st.a = ComputeA(x[0], -x[2]);
      st.Dc = TPV102Params::Dc;
   };
   fs.sigma_min = 0.0;
   return fs;
}

static void TPV102InitialState(NitscheWaveOperator &op)
{
   const int no = op.NumOwnedFaultDOFs();
   Vector d(3 * no), V(3 * no), psi(no);
   d = 0.0;
   for (int j = 0; j < no; j++)
   {
      const real_t *s = &op.FaultFrames()(9 * j + 6);
      for (int c = 0; c < 3; c++) { V(3 * j + c) = TPV102Params::V_ini * s[c]; }
      psi(j) = ComputeInitialPsi(op.FaultA()(j));
   }
   op.SetFaultState(d, V, psi);
}

/// sigma*^0 = the initial sigma_n per owned DOF (as the driver seeds it).
static void SeedPC(NitscheWaveOperator &op, real_t v_star, real_t length)
{
   const int no = op.NumOwnedFaultDOFs();
   Vector S, sig0(no);
   op.BackgroundTraction(0.0, S);
   for (int j = 0; j < no; j++)
   {
      const real_t *f = &op.FaultFrames()(9 * j);
      real_t X0, X1;
      fault_frame::ToLocal(&S(3 * j), 0.0, f, f + 3, f + 6, X0, X1, sig0(j));
   }
   op.SetPrakashClifton(v_star, length);
   op.SetPrakashCliftonState(sig0);
}

struct RunState
{
   Vector d, V, psi, sig, X, pc;
   real_t min_sig = 1e300, max_sig = -1e300;
   int nsteps = 0;
};

/// Central-difference steps from t = 0 to tfinal; PC off when length <= 0.
static RunState Run(ParMesh &pm, int p, bool nucleation, real_t v_star, real_t length,
                    real_t tfinal)
{
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl, nucleation);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   op.SetTimeStep(dt);
   TPV102InitialState(op);
   if (length > 0.0) { SeedPC(op, v_star, length); }
   const int nsteps = static_cast<int>(std::ceil(tfinal / dt - 1e-9));
   const int n = op.Height();
   Vector U(n), Ud(n);
   U = 0.0;
   Ud = 0.0;
   real_t t = 0.0;
   op.SetTime(t);
   CentralDifferenceSolver ode;
   ode.Init(op);
   RunState r;
   for (int s = 0; s < nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
      r.min_sig = std::min(r.min_sig, op.LastStep().min_sigma_eff);
      r.max_sig = std::max(r.max_sig, op.LastStep().max_sigma_eff);
   }
   r.nsteps = nsteps;
   r.d = op.Slip(); r.V = op.SlipRate(); r.psi = op.State();
   r.sig = op.LastSigmaEff(); r.X = op.LastShear();
   if (length > 0.0) { r.pc = op.PrakashCliftonSigma(); }
   return r;
}

static bool Same(const Vector &a, const Vector &b)
{
   return a.Size() == b.Size()
          && std::memcmp(a.GetData(), b.GetData(), sizeof(real_t) * a.Size()) == 0;
}

// ------------------------------------------------------------------ PC-b
static void TestExactUpdate(ParMesh &pm, int p)
{
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl, false);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   op.SetTimeStep(dt);
   TPV102InitialState(op);
   SeedPC(op, kVStar, kLength);
   real_t worst = 0.0;
   int checked = 0;
   if (op.NumOwnedFaultDOFs() > 0)
   {
      const int j = 0;
      Vector S;
      op.BackgroundTraction(0.0, S);
      const real_t *n = &op.FaultFrames()(9 * j);
      const real_t sig0 = op.PrakashCliftonSigma()(j), sig1 = sig0 + 10.0e6;
      // sigma_eff = -(TS . n): shift TS along n by (sig0 - sig1)
      real_t TS[3];
      for (int c = 0; c < 3; c++) { TS[c] = S(3 * j + c) + (sig0 - sig1) * n[c]; }
      const real_t V = 0.5;   // m/s, held constant
      const real_t x = -(V + kVStar) * dt / kLength;
      const real_t E = std::exp(x), omE = -std::expm1(x);
      real_t sig_pc = sig0;
      for (int k = 1; k <= 200; k++)
      {
         real_t Vn[3], sig, X[2];
         op.CloseDOF(j, TS, op.State()(j), Vn, sig, X, &sig_pc, E, omE);
         const real_t exact = sig1 + (sig0 - sig1) * std::exp(-(V + kVStar) * k * dt / kLength);
         worst = std::max(worst, std::abs(sig_pc - exact) / sig0);
         checked++;
      }
   }
   worst = RMax(worst, pm.GetComm());
   checked = ISum(checked, pm.GetComm());
   if (g_rank == 0)
   {
      std::cout << "  PC-b p = " << p << ": " << checked << " closes, max |sigma* - exact| / sigma*^0 = "
                << worst << "\n";
   }
   TEST_ASSERT(checked > 0 && worst < 1e-13,
               "PC-b p=" << p << " sigma* follows the exact exponential relaxation to round-off");
}

// ------------------------------------------------------------------ PC-c
static void TestTinyLength(ParMesh &pm, int p, real_t T)
{
   RunState off = Run(pm, p, true, 0.0, 0.0, T);
   RunState on = Run(pm, p, true, kVStar, 1e-200, T);
   const int nsteps = on.nsteps;
   int bad = 0;
   if (!(Same(off.d, on.d) && Same(off.V, on.V) && Same(off.psi, on.psi)
         && Same(off.sig, on.sig) && Same(off.X, on.X))) { bad = 1; }
   if (!Same(on.pc, on.sig)) { bad += 2; }
   bad = ISum(bad, pm.GetComm());
   real_t vmax = RMax(off.V.Normlinf(), pm.GetComm());
   if (g_rank == 0)
   {
      std::cout << "  PC-c p = " << p << ": " << nsteps << " steps, max|V| = " << vmax
                << " m/s; mismatch flags (sum over ranks) " << bad << "\n";
   }
   TEST_ASSERT(bad == 0, "PC-c p=" << p << " L = 1e-200: state bitwise equal to PC off, sigma* == sigma_eff");
}

// ------------------------------------------------------------------ PC-d
static void TestEquilibrium(ParMesh &pm, int p, real_t T)
{
   RunState on = Run(pm, p, false, kVStar, kLength, T);
   const int nsteps = on.nsteps;
   real_t dV = 0.0, dS = 0.0;
   const int no = on.psi.Size();
   for (int j = 0; j < no; j++)
   {
      const real_t *V = &on.V(3 * j);
      const real_t Va = std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
      dV = std::max(dV, std::abs(Va - TPV102Params::V_ini) / TPV102Params::V_ini);
      dS = std::max(dS, std::abs(on.pc(j) - TPV102Params::sigma_n) / TPV102Params::sigma_n);
   }
   dV = RMax(dV, pm.GetComm());
   dS = RMax(dS, pm.GetComm());
   if (g_rank == 0)
   {
      std::cout << "  PC-d p = " << p << ": " << nsteps << " steps, max ||V| - V_ini| / V_ini = " << dV
                << ", max |sigma* - 120 MPa| / 120 MPa = " << dS << "\n";
   }
   TEST_ASSERT(dV < 1e-6 && dS < 1e-9,
               "PC-d p=" << p << " equilibrium: |V| stays V_ini, sigma* stays 120 MPa");
}

// ------------------------------------------------------------------ PC-s
static void TestActive(ParMesh &pm, int p, real_t T)
{
   // spread of sigma* and of the run's last sigma_eff (max - min over DOFs)
   RunState off = Run(pm, p, true, 0.0, 0.0, T);
   RunState on = Run(pm, p, true, kVStar, kLength, T);
   const int nsteps = on.nsteps;
   real_t dV = 0.0, lo = 1e300, hi = -1e300;
   for (int i = 0; i < on.V.Size(); i++) { dV = std::max(dV, std::abs(on.V(i) - off.V(i))); }
   for (int j = 0; j < on.pc.Size(); j++) { lo = std::min(lo, on.pc(j)); hi = std::max(hi, on.pc(j)); }
   dV = RMax(dV, pm.GetComm());
   const real_t vmax = RMax(on.V.Normlinf(), pm.GetComm());
   hi = RMax(hi, pm.GetComm());
   lo = -RMax(-lo, pm.GetComm());
   const real_t smin = -RMax(-on.min_sig, pm.GetComm()), smax = RMax(on.max_sig, pm.GetComm());
   real_t elo = 1e300, ehi = -1e300;
   for (int j = 0; j < on.sig.Size(); j++) { elo = std::min(elo, on.sig(j)); ehi = std::max(ehi, on.sig(j)); }
   ehi = RMax(ehi, pm.GetComm());
   elo = -RMax(-elo, pm.GetComm());
   if (g_rank == 0)
   {
      std::cout << "  PC-s p = " << p << ": " << nsteps << " steps, max|V| = " << vmax
                << " m/s, max|V_on - V_off| = " << dV << " m/s; sigma* in [" << lo << ", " << hi
                << "], sigma_eff over the run in [" << smin << ", " << smax << "], last sigma_eff in ["
                << elo << ", " << ehi << "]\n";
   }
   TEST_ASSERT(std::isfinite(vmax) && dV > 0.0,
               "PC-s p=" << p << " v_star = 1 m/s, L = 0.01 m changes the solution and stays finite");
   TEST_ASSERT(lo >= smin && hi <= smax,
               "PC-s p=" << p << " sigma* stays within the range sigma_eff took (a convex average)");
   // (sigma* is a lagged average, so its spread at one instant need not be
   // smaller than sigma_eff's; the smoothing is judged on the TPV102 runs.)
   // p = 1 does not rupture on this coarse box by T, sigma_eff stays flat.
   TEST_ASSERT(p == 1 || ehi - elo > 0.0,
               "PC-s p=" << p << " sigma_eff varies on the jittered box (mirror symmetry broken)");
}

// ------------------------------------------------------------------ PC-a
static int AbortUnseeded(ParMesh &pm)
{
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl, false);
   NitscheWaveOperator op(pm, 1, TPV102Material(), TPV102BC(), {}, &fs);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, 1, kCuse[1]);
   op.SetTimeStep(dt);
   TPV102InitialState(op);
   op.SetPrakashClifton(kVStar, kLength);   // no SetPrakashCliftonState
   const int n = op.Height();
   Vector U(n), Ud(n);
   U = 0.0;
   Ud = 0.0;
   real_t t = 0.0, h = dt;
   op.SetTime(t);
   CentralDifferenceSolver ode;
   ode.Init(op);
   ode.Step(U, Ud, t, h);   // must abort
   return 0;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   std::unique_ptr<Mesh> m(TPV102Box(4, 4, 2));
   ParMesh pm(MPI_COMM_WORLD, *m);
   std::unique_ptr<Mesh> mj(TPV102Box(4, 4, 2, true));
   ParMesh pmj(MPI_COMM_WORLD, *mj);
   if (argc > 1 && std::string(argv[1]) == "--abort-unseeded") { return AbortUnseeded(pm); }
   if (g_rank == 0)
   {
      std::cout << "=== Prakash-Clifton (v_star = " << kVStar << " m/s, L = " << kLength
                << " m) on " << Mpi::WorldSize() << " rank(s) ===\n";
   }
   for (int p = 1; p <= 2; p++)
   {
      TestExactUpdate(pm, p);
      TestTinyLength(pm, p, 2.0);    // nucleation ramps over 1 s
      TestEquilibrium(pm, p, 2.0);
      TestActive(pmj, p, 3.0);
   }
   if (g_rank == 0)
   {
      std::cout << "=== Summary: " << num_passed << "/" << num_tests << " passed, "
                << num_failed << " failed ===\n";
   }
   return num_failed == 0 ? 0 : 1;
}
