// Nitsche dynamic rupture: linear slip weakening on the inline box
// (benchmark plan rev 1.3, phase B3; the slip-weakening counterpart of
// test_nitsche_tpv102_inline.cpp).
//
// Mesh: the mirror-symmetric 60 x 60 x 30 km box of the P6 fixture
// (24 tets per hex, fault on the whole plane Y = 0 (attr 3), free surface
// Z = 0 (attr 1), absorbing walls (attr 5)).  Material as TPV102 (= TPV205).
// Friction: SlipWeakeningClosure with the TPV205 values mu_s = 0.677,
// mu_d = 0.525, Dc = 0.40 m, no cohesion; prestress sigma_XY = 70 MPa,
// sigma_YY = -120 MPa.  Nucleation smooth in time, so the dt studies
// measure the scheme and not the start: dtau F(r) G(t) with the TPV102
// shapes F(r) = exp(r^2/(r^2 - R^2)), G(t) of rise time 1 s, at the TPV102
// hypocentre (0, -7.5 km), R = 9 km, dtau = 15 MPa (peak 85 MPa against a
// static strength of 81.24 MPa).  Start at rest: d = 0, V = 0, slip path 0.
//
// Modes:
//   (default)  at p = 1 and 2 on the 4x4x2 box to 3 s:
//              L-05  mirror symmetry: |sigma_eff - 120 MPa| < 1e-8 x 120 MPa
//                    at every DOF and step (T-13(a) analogue); the first
//                    DOFs to slip (|V| > 0) lie inside the nucleation patch
//                    and start inside the nucleation time.
//              L-08  rank invariance: fault data (d, V, slip path, sigma_eff)
//                    at the end agree with a one-rank run of the same
//                    problem (done on rank 0 on MPI_COMM_SELF) to 1e-11
//                    relative; DOFs matched by position and face centroid.
//              L-09  Prakash-Clifton (v* = 1 m/s, L = 0.01 m) and
//                    Kelvin-Voigt (eta_bar = 0.1, whole domain) switch on
//                    and run finite; with L = 1e-200 Prakash-Clifton is bit
//                    for bit the run without it, and sigma* == sigma_eff.
//   --l06 nx ny nz p tfinal
//              L-06  energy, Kelvin-Voigt off: nodal E4 (T_h + S).V <= 1e-12
//                    x scale at every DOF and step; the accumulated balance
//                    residual falls by a factor in [3.5, 4.5] per dt halving
//                    (the T-15 criterion).  If not, the test fails and
//                    reports the ratios; the bound is not loosened.
//   --l07 nx ny nz p tfinal [base] [kv_eta]
//              L-07  coupled dt refinement at dt0/base x {1, 2, 4, 8}: self-
//                    convergence ratios |u_k - u_2k| / |u_2k - u_4k| of slip,
//                    slip rate, endpoint shear and station velocity (max over
//                    DOFs and the dt0 output times).  Order 1 is expected
//                    (plan A.1 Remark); the observed order is printed.
//                    kv_eta (default 0.1, the benchmark decks' Kelvin-Voigt
//                    eta_bar) > 0: the gate, every slip and slip-rate ratio
//                    >= 1.6.  kv_eta = 0: undamped diagnostic, every slip,
//                    shear and station-velocity ratio >= 1.6; the slip rate
//                    is reported, not asserted: without damping the LSW step
//                    carries a 2-step zigzag in V (V^{n+1} ~ -V^n + smooth,
//                    since X* has slope -gamma in V^n and the LSW strength
//                    does not depend on V) that does not shrink with dt.
//                    Decided by Chunhui Zhao, 2026-10-16: gate L-07 damped.
//
// Usage: mpirun -np 1|2|4 ./seas_test_nitsche_lsw_inline [mode ...]
//
// Scope (Chunhui Zhao, 2026-10-16).

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../domain/nitsche_energy.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../friction/slip_weakening_closure.hpp"
#include "../../config/tpv102_params.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
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
static const real_t kCuse[5] = {0.0, 0.057, 0.11, 0.14, 0.166};   // 0.5 C_eig (P3)

// LSW problem constants (TPV205 friction; the nucleation is this test's own)
static const real_t kMuS = 0.677, kMuD = 0.525, kDc = 0.40;
static const real_t kTau0 = 70.0e6, kSigma0 = 120.0e6;
static const real_t kNucR = 9.0e3, kNucDtau = 15.0e6;
static const real_t kVStar = 1.0, kLength = 0.01;     // Prakash-Clifton, SeisSol TPV6/7

static int ISum(int v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static real_t RMin(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MIN, c); return v; }

// ---------------------------------------------------------------- fixture
// (as test_nitsche_pc.cpp: TPV102Box without jitter, TPV102BC, TPV102Material)

static Mesh *Box(int nx, int ny, int nz)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2 * kL, 2 * kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
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

static BoundaryConfig BC()
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.absorbing_attrs = {5};
   return bc;
}

static MaterialField Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu,
                                      TPV102Params::rho);
}

static void Hypo(real_t h[3])
{
   h[0] = TPV102Params::hypo_along_strike;
   h[1] = 0.0;
   h[2] = -TPV102Params::hypo_down_dip;
}

static real_t NucF(real_t r)
{
   if (r >= kNucR) { return 0.0; }
   const real_t r2 = r * r, R2 = kNucR * kNucR;
   return std::exp(r2 / (r2 - R2));
}

static NitscheFaultSetup LSWSetup(const IFaultClosure &cl)
{
   NitscheFaultSetup fs;
   fs.closure = &cl;
   fs.prestress = [](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = kTau0;
      s[4] = -kSigma0;
   };
   fs.nucleation_spatial = [](const real_t x[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      real_t h[3];
      Hypo(h);
      const real_t v = kNucDtau * NucF(std::hypot(x[0] - h[0], x[2] - h[2]));
      s[1] = s[3] = v;
   };
   fs.nucleation_temporal = [](real_t t) { return NucleationTemporal(t); };
   fs.parameters = [](const real_t[3], FaultPointState &st) {
      st.mu_s = kMuS;
      st.mu_d = kMuD;
      st.Dc = kDc;
      st.cohesion = 0.0;
   };
   fs.sigma_min = 0.0;
   return fs;
}

/// At rest: d = 0, V = 0, slip path from SeedEquilibrium (0).
static void InitialState(NitscheWaveOperator &op)
{
   const int no = op.NumOwnedFaultDOFs();
   Vector d(3 * no), V(3 * no), psi(no);
   d = 0.0;
   V = 0.0;
   for (int j = 0; j < no; j++)
   {
      FaultPointState st = op.FaultParams(j);
      op.Closure().SeedEquilibrium(kTau0, 0.0, kSigma0, st);
      psi(j) = st.psi;
   }
   op.SetFaultState(d, V, psi);
}

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

/// Face centroid per owned fault DOF (3 no): with the position it keys a DG
/// fault DOF across rank counts (as test_nitsche_tpv102_inline.cpp).
static void OwnedFaceCentroids(NitscheWaveOperator &op, Vector &c_own)
{
   ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
   const int nd = E.GetNumFaultDOFs(), nbf = E.GetNbfPerFace();
   Vector c_loc(3 * nd);
   for (int k = 0; k < nd; k++)
   {
      long long key[3];
      real_t xv[3][3];
      op.Orientation().FaceVertexData(k / nbf, key, xv);
      for (int d = 0; d < 3; d++) { c_loc(3 * k + d) = (xv[0][d] + xv[1][d] + xv[2][d]) / 3.0; }
   }
   E.RestrictToOwnedFault(c_loc, c_own, 3);
}

// ------------------------------------------------------------ one run

struct RunConfig
{
   int nx = 4, ny = 4, nz = 2, p = 1;
   real_t tfinal = 3.0;
   real_t dt_scale = 1.0;
   real_t pc_length = 0.0;    ///< > 0: Prakash-Clifton on with this L
   real_t kv_eta = 0.0;       ///< > 0: Kelvin-Voigt on, whole domain
   bool self = false;         ///< run on MPI_COMM_SELF (rank 0 only)
};

/// One record per owned fault DOF: key (position, face centroid; mm) and
/// values (d 3, V 3, slip path, sigma_eff).
struct DOFRow
{
   std::array<long long, 6> key;
   std::array<real_t, 8> val;
};

struct RunResult
{
   bool finite = true;
   long steps = 0;
   real_t dt = 0.0;
   int fault_dofs = 0;
   real_t sig_min = 1e300, sig_max = -1e300, Vmax = 0.0;
   real_t hypo_dist = 0.0;
   real_t first_t = -1.0, first_dist = -1.0;   ///< first |V| > 0; farthest such DOF
   int slipped = 0, stuck = 0;           ///< |d| > 1 mm; V == 0 exactly at the end
   real_t x_min = 0.0, x_max = 0.0, dmax = 0.0;
   Vector d, V, psi, sig, X, pc;         ///< owned, end of run
   std::vector<DOFRow> rows;             ///< gathered on rank 0 of the run's comm
};

static RunResult Run(const RunConfig &rc)
{
   RunResult R;
   MPI_Comm comm = rc.self ? MPI_COMM_SELF : MPI_COMM_WORLD;
   std::unique_ptr<Mesh> m(Box(rc.nx, rc.ny, rc.nz));
   ParMesh pm(comm, *m);
   m.reset();
   const int crank = [&] { int r; MPI_Comm_rank(comm, &r); return r; }();
   SlipWeakeningClosure cl;
   NitscheFaultSetup fs = LSWSetup(cl);
   NitscheWaveOperator op(pm, rc.p, Material(), BC(), {}, &fs);
   const real_t dt = rc.dt_scale * ComputeMaxDt(op.Sample(), pm, rc.p, kCuse[rc.p]);
   op.SetTimeStep(dt);
   if (rc.kv_eta > 0.0)
   {
      op.SetKelvinVoigt(rc.kv_eta, NitscheWaveOperator::KVScale::Node, -1.0);
   }
   InitialState(op);
   if (rc.pc_length > 0.0) { SeedPC(op, kVStar, rc.pc_length); }
   R.dt = dt;
   const int no = op.NumOwnedFaultDOFs();
   R.fault_dofs = ISum(no, comm);

   real_t h[3];
   Hypo(h);
   const Vector &Xp = op.FaultPositions();
   std::vector<real_t> hdist(no);
   real_t hmin = 1e300;
   for (int j = 0; j < no; j++)
   {
      hdist[j] = std::hypot(Xp(3 * j) - h[0], Xp(3 * j + 1) - h[1], Xp(3 * j + 2) - h[2]);
      hmin = std::min(hmin, hdist[j]);
   }
   R.hypo_dist = RMin(hmin, comm);
   std::vector<long> first_step(no, -1);

   const int n = op.Height();
   Vector U(n), Ud(n);
   U = 0.0;
   Ud = 0.0;
   real_t t = 0.0;
   op.SetTime(t);
   CentralDifferenceSolver ode;
   ode.Init(op);
   const long nsteps = std::lround(std::ceil(rc.tfinal / dt - 1e-9));
   for (long s = 0; s < nsteps; s++)
   {
      real_t h_s = dt;
      ode.Step(U, Ud, t, h_s);
      const NitscheStepRecord &rec = op.LastStep();
      R.sig_min = std::min(R.sig_min, rec.min_sigma_eff);
      R.sig_max = std::max(R.sig_max, rec.max_sigma_eff);
      R.Vmax = std::max(R.Vmax, rec.max_V);
      for (int j = 0; j < no; j++)
      {
         const real_t *V = &op.SlipRate()(3 * j);
         if (first_step[j] < 0 && (V[0] != 0.0 || V[1] != 0.0 || V[2] != 0.0)) { first_step[j] = s + 1; }
      }
      if (!std::isfinite(rec.max_V) || !std::isfinite(rec.min_sigma_eff))
      {
         R.finite = false;
         break;
      }
   }
   R.steps = op.LastStep().step;
   R.finite = (ISum(R.finite ? 0 : 1, comm) == 0) && std::isfinite(RMax(U.Normlinf(), comm));
   {
      // first step with |V| > 0 over all DOFs, and the farthest DOF that
      // starts at that step (ties included: partition-invariant)
      long smin = std::numeric_limits<long>::max();
      for (int j = 0; j < no; j++) { if (first_step[j] > 0) { smin = std::min(smin, first_step[j]); } }
      MPI_Allreduce(MPI_IN_PLACE, &smin, 1, MPI_LONG, MPI_MIN, comm);
      real_t far = -1.0;
      for (int j = 0; j < no; j++) { if (first_step[j] == smin) { far = std::max(far, hdist[j]); } }
      far = RMax(far, comm);
      if (smin != std::numeric_limits<long>::max())
      {
         R.first_t = smin * dt;
         R.first_dist = far;
      }
   }
   R.d = op.Slip(); R.V = op.SlipRate(); R.psi = op.State();
   R.sig = op.LastSigmaEff(); R.X = op.LastShear();
   if (rc.pc_length > 0.0) { R.pc = op.PrakashCliftonSigma(); }
   real_t xmin = 1e300, xmax = -1e300, dmax = 0.0;
   int slipped = 0, stuck = 0;
   for (int j = 0; j < no; j++)
   {
      const real_t *d = &R.d(3 * j), *V = &R.V(3 * j);
      const real_t da = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      dmax = std::max(dmax, da);
      if (da > 1e-3)
      {
         slipped++;
         xmin = std::min(xmin, Xp(3 * j));
         xmax = std::max(xmax, Xp(3 * j));
      }
      if (V[0] == 0.0 && V[1] == 0.0 && V[2] == 0.0) { stuck++; }
   }
   R.slipped = ISum(slipped, comm);
   R.stuck = ISum(stuck, comm);
   R.x_min = RMin(xmin, comm);
   R.x_max = RMax(xmax, comm);
   R.dmax = RMax(dmax, comm);

   // gather keyed rows on the comm's rank 0
   Vector cf;
   OwnedFaceCentroids(op, cf);
   std::vector<double> buf(14 * no);
   for (int j = 0; j < no; j++)
   {
      double *b = &buf[14 * j];
      for (int c = 0; c < 3; c++) { b[c] = Xp(3 * j + c); b[3 + c] = cf(3 * j + c); }
      for (int c = 0; c < 3; c++) { b[6 + c] = R.d(3 * j + c); b[9 + c] = R.V(3 * j + c); }
      b[12] = R.psi(j);
      b[13] = R.sig(j);
   }
   int nloc = 14 * no, csize;
   MPI_Comm_size(comm, &csize);
   std::vector<int> counts(csize), displs(csize, 0);
   MPI_Gather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
   std::vector<double> all;
   if (crank == 0)
   {
      for (int r = 1; r < csize; r++) { displs[r] = displs[r - 1] + counts[r - 1]; }
      all.resize(displs[csize - 1] + counts[csize - 1]);
   }
   MPI_Gatherv(buf.data(), nloc, MPI_DOUBLE, all.data(), counts.data(), displs.data(),
               MPI_DOUBLE, 0, comm);
   if (crank == 0)
   {
      const size_t nr = all.size() / 14;
      R.rows.resize(nr);
      for (size_t k = 0; k < nr; k++)
      {
         const double *b = &all[14 * k];
         for (int c = 0; c < 6; c++) { R.rows[k].key[c] = std::llround(b[c] * 1e3); }
         for (int c = 0; c < 8; c++) { R.rows[k].val[c] = b[6 + c]; }
      }
      std::sort(R.rows.begin(), R.rows.end(),
                [](const DOFRow &a, const DOFRow &b) { return a.key < b.key; });
   }
   return R;
}

static bool Same(const Vector &a, const Vector &b)
{
   return a.Size() == b.Size()
          && std::memcmp(a.GetData(), b.GetData(), sizeof(real_t) * a.Size()) == 0;
}

static void Report(const char *label, const RunConfig &rc, const RunResult &R)
{
   if (g_rank != 0) { return; }
   std::cout << std::setprecision(6) << "  " << label << ": p = " << rc.p << ", dt = " << R.dt
             << " s, " << R.steps << " steps, " << R.fault_dofs << " fault DOFs (nearest "
             << R.hypo_dist << " m from the hypocentre); first slip at t = " << R.first_t
             << " s, up to " << R.first_dist << " m from the hypocentre; max |V| = " << R.Vmax << " m/s; " << R.slipped << " DOFs slipped > 1 mm (max "
             << R.dmax << " m, x in [" << R.x_min << ", " << R.x_max << "] m), " << R.stuck
             << " stuck (V == 0) at the end; sigma_eff in [" << R.sig_min << ", " << R.sig_max
             << "] Pa\n";
}

// ------------------------------------------------------- L-05, L-08, L-09
static void TestDefault(int p)
{
   if (g_rank == 0)
   {
      std::cout << "\n=== L-05, L-08, L-09: 4x4x2 box, p = " << p << ", 3 s, "
                << Mpi::WorldSize() << " rank(s) ===\n";
   }
   const std::string tag = "p=" + std::to_string(p);
   RunConfig rc;
   rc.p = p;
   const RunResult R = Run(rc);
   Report("plain", rc, R);
   const real_t dev = std::max(std::abs(R.sig_min - kSigma0), std::abs(R.sig_max - kSigma0));
   TEST_ASSERT(R.finite, "L-05 " << tag << " finite for the whole run");
   TEST_ASSERT(dev < 1e-8 * kSigma0,
               "L-05 " << tag << " |sigma_eff - 120 MPa| < 1e-8 x 120 MPa at every DOF and step");
   TEST_ASSERT(R.first_t > 0.0 && R.first_t <= TPV102Params::nuc_T && R.first_dist >= 0.0
               && R.first_dist < kNucR && R.slipped > 0,
               "L-05 " << tag << " slip starts inside the nucleation patch within the nucleation time");

   // L-08: against a one-rank run on rank 0
   {
      int bad_keys = 0;
      real_t rel[4] = {0.0, 0.0, 0.0, 0.0};   // d, V, slip path, sigma_eff
      if (g_rank == 0)
      {
         RunResult S;
         if (Mpi::WorldSize() > 1)
         {
            RunConfig rs = rc;
            rs.self = true;
            S = Run(rs);
         }
         const std::vector<DOFRow> &ref = (Mpi::WorldSize() > 1) ? S.rows : R.rows;
         if (ref.size() != R.rows.size()) { bad_keys = 1; }
         else
         {
            real_t scale[4] = {0.0, 0.0, 0.0, 0.0}, diff[4] = {0.0, 0.0, 0.0, 0.0};
            const int f0[4] = {0, 3, 6, 7}, f1[4] = {3, 6, 7, 8};
            for (size_t k = 0; k < ref.size(); k++)
            {
               if (ref[k].key != R.rows[k].key) { bad_keys++; continue; }
               for (int g = 0; g < 4; g++)
                  for (int c = f0[g]; c < f1[g]; c++)
                  {
                     scale[g] = std::max(scale[g], std::abs(ref[k].val[c]));
                     diff[g] = std::max(diff[g], std::abs(ref[k].val[c] - R.rows[k].val[c]));
                  }
            }
            for (int g = 0; g < 4; g++) { rel[g] = diff[g] / std::max(scale[g], 1e-300); }
         }
         std::cout << std::setprecision(4) << "  L-08 against one rank: " << ref.size()
                   << " DOFs, key mismatches " << bad_keys << "; relative max differences: d "
                   << rel[0] << ", V " << rel[1] << ", slip path " << rel[2] << ", sigma_eff "
                   << rel[3] << "\n";
      }
      MPI_Bcast(&bad_keys, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(rel, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      TEST_ASSERT(bad_keys == 0 && rel[0] <= 1e-11 && rel[1] <= 1e-11 && rel[2] <= 1e-11
                  && rel[3] <= 1e-11,
                  "L-08 " << tag << " fault data agree with one rank to 1e-11 relative");
   }

   // L-09: Prakash-Clifton and Kelvin-Voigt
   {
      RunConfig rp = rc;
      rp.pc_length = kLength;
      const RunResult P = Run(rp);
      Report("Prakash-Clifton", rp, P);
      real_t lo = 1e300, hi = -1e300;
      for (int j = 0; j < P.pc.Size(); j++) { lo = std::min(lo, P.pc(j)); hi = std::max(hi, P.pc(j)); }
      lo = RMin(lo, MPI_COMM_WORLD);
      hi = RMax(hi, MPI_COMM_WORLD);
      TEST_ASSERT(P.finite && P.slipped > 0 && lo >= P.sig_min - 1e-6 * kSigma0
                  && hi <= P.sig_max + 1e-6 * kSigma0,
                  "L-09 " << tag << " Prakash-Clifton on: finite, slips, sigma* inside the sigma_eff range");

      RunConfig rt = rc;
      rt.pc_length = 1e-200;
      const RunResult T = Run(rt);
      int bad = 0;
      if (!(Same(R.d, T.d) && Same(R.V, T.V) && Same(R.psi, T.psi) && Same(R.sig, T.sig)
            && Same(R.X, T.X))) { bad = 1; }
      if (!Same(T.pc, T.sig)) { bad += 2; }
      bad = ISum(bad, MPI_COMM_WORLD);
      TEST_ASSERT(bad == 0, "L-09 " << tag << " Prakash-Clifton L = 1e-200: bit for bit the run "
                  "without it, sigma* == sigma_eff");

      RunConfig rk = rc;
      rk.kv_eta = 0.1;
      const RunResult K = Run(rk);
      Report("Kelvin-Voigt 0.1", rk, K);
      TEST_ASSERT(K.finite && K.slipped > 0, "L-09 " << tag << " Kelvin-Voigt on: finite, slips");
   }
}

// ------------------------------------------------------------------ L-06
struct EnergyRun
{
   real_t dt = 0.0, sum_rho = 0.0, sum_abs = 0.0, Emax = 0.0;
   real_t node_max = 0.0, node_scale = 0.0;
   real_t W_pre = 0.0, W_nuc = 0.0, W_fric = 0.0, W_abs = 0.0;
   NitscheEnergyRecord last;
   long steps = 0;
};

static EnergyRun RunEnergy(int nx, int ny, int nz, int p, real_t tfinal, real_t dt_scale)
{
   EnergyRun R;
   std::unique_ptr<Mesh> m(Box(nx, ny, nz));
   ParMesh pm(MPI_COMM_WORLD, *m);
   m.reset();
   SlipWeakeningClosure cl;
   NitscheFaultSetup fs = LSWSetup(cl);
   NitscheWaveOperator op(pm, p, Material(), BC(), {}, &fs);
   const real_t dt = dt_scale * ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   R.dt = dt;
   op.SetTimeStep(dt);
   InitialState(op);
   NitscheEnergyMonitor mon(op, false);
   const int n = op.Height();
   Vector U(n), Ud(n), A0(n);
   U = 0.0; Ud = 0.0;
   op.SetTime(0.0);
   op.Mult(U, Ud, A0);
   NitscheEnergyRecord r = mon.Evaluate(U, Ud, A0, 0.0);
   mon.Start(r);
   CentralDifferenceSolver ode;
   ode.Init(op);
   real_t t = 0.0;
   const long nsteps = std::lround(std::ceil(tfinal / dt - 1e-9));
   NitscheEnergyRecord prev = r;
   for (long s = 0; s < nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
      r = mon.Evaluate(U, Ud, op.LastAcceleration(), t);
      mon.Step(r);
      R.W_pre += 0.5 * dt * (r.P_pre + prev.P_pre);
      R.W_nuc += 0.5 * dt * (r.P_nuc + prev.P_nuc);
      R.W_fric += 0.5 * dt * (r.P_fric + prev.P_fric);
      R.W_abs += 0.5 * dt * (r.P_abs + prev.P_abs);
      R.Emax = std::max(R.Emax, std::abs(r.E_h));
      prev = r;
   }
   R.sum_rho = mon.SumRho();
   R.sum_abs = mon.SumAbsRho();
   R.node_max = mon.NodalFrictionMax();
   R.node_scale = mon.NodalFrictionScale();
   R.last = r;
   R.steps = nsteps;
   return R;
}

static void TestL06(int nx, int ny, int nz, int p, real_t tfinal)
{
   if (g_rank == 0)
   {
      std::cout << "\n=== L-06 energy balance (Kelvin-Voigt off): " << nx << "x" << ny << "x"
                << nz << " box, p = " << p << ", " << tfinal << " s ===\n";
   }
   EnergyRun a = RunEnergy(nx, ny, nz, p, tfinal, 1.0);
   EnergyRun b = RunEnergy(nx, ny, nz, p, tfinal, 0.5);
   EnergyRun c = RunEnergy(nx, ny, nz, p, tfinal, 0.25);
   const real_t r1 = std::abs(a.sum_rho / b.sum_rho), r2 = std::abs(b.sum_rho / c.sum_rho);
   const real_t scale = std::max({std::abs(a.W_pre), std::abs(a.W_fric), a.Emax});
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6);
      for (const EnergyRun *e : {&a, &b, &c})
      {
         std::cout << "  dt = " << e->dt << " (" << e->steps << " steps): E_h(T) = " << e->last.E_h
                   << " J; work: prestress " << e->W_pre << ", nucleation " << e->W_nuc
                   << ", friction " << e->W_fric << ", dashpot " << e->W_abs << "; sum rho = "
                   << e->sum_rho << " (sum |rho| " << e->sum_abs << "); nodal max (T+S).V = "
                   << e->node_max << " (scale " << e->node_scale << ")\n";
      }
      std::cout << "  accumulated residual ratios: dt0/dt1 = " << r1 << ", dt1/dt2 = " << r2
                << "; |sum rho| / work scale = " << std::abs(a.sum_rho) / scale << "\n";
   }
   TEST_ASSERT(r1 >= 3.5 && r1 <= 4.5 && r2 >= 3.5 && r2 <= 4.5,
               "L-06 accumulated residual falls by [3.5, 4.5] per dt halving (T-15 criterion)");
   const real_t nmax = std::max({a.node_max, b.node_max, c.node_max});
   const real_t nsc = std::max({a.node_scale, b.node_scale, c.node_scale});
   TEST_ASSERT(nmax <= 1e-12 * nsc, "L-06 E4: (T_h + S).V <= 0 at every fault DOF and step");
}

// ------------------------------------------------------------------ L-07
struct Series
{
   std::vector<Vector> d, V, tau, st;
   std::vector<Vector> Vm;   ///< (V^{n-1} + V^n)/2 at the output step n
   std::vector<Vector> zz;   ///< |V^{n+1} - 2 V^n + V^{n-1}| / 4 at the output step (per DOF, max comp)
   Vector X;   ///< owned fault DOF positions (3 no)
};

static const real_t kStations[3][3] = {
   {0.0, 1000.0, -7500.0}, {6000.0, -1000.0, -7500.0}, {0.0, 1000.0, -2000.0}};

static Series RunSeries(ParMesh &pm, int p, real_t dt0, int refine, real_t tfinal,
                        real_t kv_eta = 0.0)
{
   Series R;
   SlipWeakeningClosure cl;
   NitscheFaultSetup fs = LSWSetup(cl);
   NitscheWaveOperator op(pm, p, Material(), BC(), {}, &fs);
   const real_t dt = dt0 / refine;
   op.SetTimeStep(dt);
   if (kv_eta > 0.0) { op.SetKelvinVoigt(kv_eta, NitscheWaveOperator::KVScale::Node, -1.0); }
   InitialState(op);
   const int no = op.NumOwnedFaultDOFs(), n = op.Height();
   R.X = op.FaultPositions();
   DenseMatrix pts(3, 3);
   for (int s = 0; s < 3; s++) for (int c = 0; c < 3; c++) { pts(c, s) = kStations[s][c]; }
   Array<int> el;
   Array<IntegrationPoint> ips;
   pm.FindPoints(pts, el, ips, false);
   ParGridFunction vg(&op.Elasticity().GetFESpace());
   Vector U(n), Ud(n);
   U = 0.0; Ud = 0.0;
   op.SetTime(0.0);
   CentralDifferenceSolver ode;
   ode.Init(op);
   real_t t = 0.0;
   const long nsteps = std::lround(std::ceil(tfinal / dt0 - 1e-9)) * refine;
   Vector Vprev = op.SlipRate(), Vprev2 = op.SlipRate();
   for (long s = 0; s < nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
      {
         // second difference of the last three steps (zigzag indicator)
         const Vector &Vc = op.SlipRate();
         if ((s + 1) % refine == 0)
         {
            Vector z(no);
            for (int j = 0; j < no; j++)
            {
               real_t e = 0.0;
               for (int c = 0; c < 3; c++)
               {
                  e = std::max(e, std::abs(Vc(3 * j + c) - 2 * Vprev(3 * j + c) + Vprev2(3 * j + c)) / 4);
               }
               z(j) = e;
            }
            R.zz.push_back(z);
            Vector vm(Vc);
            vm += Vprev;
            vm *= 0.5;
            R.Vm.push_back(vm);
         }
         Vprev2 = Vprev;
         Vprev = Vc;
      }
      if ((s + 1) % refine != 0) { continue; }
      R.d.push_back(op.Slip());
      R.V.push_back(op.SlipRate());
      Vector tau(2 * no);
      for (int j = 0; j < no; j++)
      {
         const real_t *dv = &op.FaultFrames()(9 * j + 3), *sv = &op.FaultFrames()(9 * j + 6);
         const real_t *V = &op.SlipRate()(3 * j);
         const real_t g = op.ClosureGamma(j);
         tau(2 * j) = op.LastShear()(2 * j) + g * fault_frame::Dot(V, dv);
         tau(2 * j + 1) = op.LastShear()(2 * j + 1) + g * fault_frame::Dot(V, sv);
      }
      R.tau.push_back(tau);
      vg = Ud;
      Vector out(9);
      out = 0.0;
      Vector val(3);
      for (int k = 0; k < 3; k++)
      {
         if (el[k] < 0) { continue; }
         vg.GetVectorValue(el[k], ips[k], val);
         for (int c = 0; c < 3; c++) { out(3 * k + c) = val(c); }
      }
      MPI_Allreduce(MPI_IN_PLACE, out.GetData(), 9, MPI_DOUBLE, MPI_SUM, pm.GetComm());
      R.st.push_back(out);
   }
   return R;
}

static real_t SeriesDiff(const std::vector<Vector> &a, const std::vector<Vector> &b, MPI_Comm c)
{
   real_t e = 0.0;
   const size_t n = std::min(a.size(), b.size());
   for (size_t k = 0; k < n; k++)
      for (int i = 0; i < a[k].Size(); i++) { e = std::max(e, std::abs(a[k](i) - b[k](i))); }
   return RMax(e, c);
}

static void TestL07(int nx, int ny, int nz, int p, real_t tfinal, int base, real_t kv_eta = 0.0)
{
   std::unique_ptr<Mesh> m(Box(nx, ny, nz));
   ParMesh pm(MPI_COMM_WORLD, *m);
   MPI_Comm comm = pm.GetComm();
   real_t dt0;
   {
      NitscheMaterialSample smp(pm, Material());
      dt0 = ComputeMaxDt(smp, pm, p, kCuse[p]);
   }
   if (g_rank == 0)
   {
      std::cout << "\n=== L-07 coupled dt refinement: " << nx << "x" << ny << "x" << nz
                << " box, p = " << p << ", to t = " << tfinal << " s, dt0 = " << dt0
                << "; runs at dt0/" << base << ", /" << 2 * base << ", /" << 4 * base << ", /"
                << 8 * base << (kv_eta > 0.0 ? "; Kelvin-Voigt eta_bar " + std::to_string(kv_eta)
                                             : std::string()) << " ===\n";
   }
   std::vector<Series> runs;
   for (int k : {1, 2, 4, 8}) { runs.push_back(RunSeries(pm, p, dt0, k * base, tfinal, kv_eta)); }
   const char *names[4] = {"slip", "slip rate", "shear (endpoint)", "station velocity"};
   real_t D[4][3];
   for (int i = 0; i < 3; i++)
   {
      D[0][i] = SeriesDiff(runs[i].d, runs[i + 1].d, comm);
      D[1][i] = SeriesDiff(runs[i].V, runs[i + 1].V, comm);
      D[2][i] = SeriesDiff(runs[i].tau, runs[i + 1].tau, comm);
      D[3][i] = SeriesDiff(runs[i].st, runs[i + 1].st, comm);
   }
   // Damped (the gate): slip and slip rate.  Undamped (diagnostic): slip,
   // shear and station velocity; the slip rate is reported only.
   const bool damped = kv_eta > 0.0;
   const bool asserted[4] = {true, damped, !damped, !damped};
   bool ok = true;
   for (int f = 0; f < 4; f++)
   {
      const real_t q1 = D[f][0] / D[f][1], q2 = D[f][1] / D[f][2];
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "  " << std::setw(18) << names[f]
                   << ": |u_k - u_2k| = " << D[f][0] << ", " << D[f][1] << ", " << D[f][2]
                   << "  ratios " << q1 << ", " << q2 << "  (observed order "
                   << std::log2(q1) << ", " << std::log2(q2) << ")"
                   << (asserted[f] ? "" : "  [reported only]") << "\n";
      }
      if (asserted[f]) { ok = ok && q1 >= 1.6 && q2 >= 1.6; }
   }
   {
      real_t Dm[3], Z[4];
      for (int i = 0; i < 3; i++) { Dm[i] = SeriesDiff(runs[i].Vm, runs[i + 1].Vm, comm); }
      for (int i = 0; i < 4; i++)
      {
         real_t z = 0.0;
         for (const Vector &v : runs[i].zz) { z = std::max(z, v.Normlinf()); }
         Z[i] = RMax(z, comm);
      }
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "    two-step mean slip rate: |u_k - u_2k| = " << Dm[0]
                   << ", " << Dm[1] << ", " << Dm[2] << "  ratios " << Dm[0] / Dm[1] << ", "
                   << Dm[1] / Dm[2] << "\n    max |V^{n+1} - 2V^n + V^{n-1}| / 4 at the output "
                   << "steps, per run: " << Z[0] << ", " << Z[1] << ", " << Z[2] << ", " << Z[3]
                   << " m/s\n";
      }
   }
   // Diagnostics for the slip rate (information only): where the max
   // difference sits, the max over samples where both runs slip (V != 0),
   // and the time-L1 norm max_j sum_k dt0 |V_a - V_b|.
   for (int i = 0; i < 3; i++)
   {
      const Series &a = runs[i], &b = runs[i + 1];
      struct { double v; int r; } loc = {-1.0, g_rank}, glob;
      real_t where[6] = {0, 0, 0, 0, 0, 0};   // t, x, z, |V_a|, |V_b|, -
      real_t both = 0.0, l1 = 0.0;
      const int no = a.X.Size() / 3;
      std::vector<real_t> acc(no, 0.0);
      for (size_t k = 0; k < std::min(a.V.size(), b.V.size()); k++)
         for (int j = 0; j < no; j++)
         {
            const real_t *va = &a.V[k](3 * j), *vb = &b.V[k](3 * j);
            const real_t na = std::sqrt(va[0] * va[0] + va[1] * va[1] + va[2] * va[2]);
            const real_t nb = std::sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);
            real_t e = 0.0;
            for (int c = 0; c < 3; c++) { e = std::max(e, std::abs(va[c] - vb[c])); }
            if (e > loc.v)
            {
               loc.v = e;
               where[0] = (k + 1) * dt0; where[1] = a.X(3 * j); where[2] = a.X(3 * j + 2);
               where[3] = na; where[4] = nb;
            }
            if (na > 0.0 && nb > 0.0) { both = std::max(both, e); }
            acc[j] += dt0 * e;
         }
      for (int j = 0; j < no; j++) { l1 = std::max(l1, acc[j]); }
      MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
      MPI_Bcast(where, 6, MPI_DOUBLE, glob.r, comm);
      both = RMax(both, comm);
      l1 = RMax(l1, comm);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "    slip rate, pair " << i << ": max diff "
                   << glob.v << " at t = " << where[0] << " s, (x, z) = (" << where[1] << ", "
                   << where[2] << ") m, |V| " << where[3] << " vs " << where[4]
                   << "; max diff where both slip " << both << "; time-L1 " << l1 << " m\n";
      }
   }
   if (damped)
   {
      TEST_ASSERT(ok, "L-07 (Kelvin-Voigt on) slip and slip rate converge in dt "
                      "(every ratio >= 1.6; order 1 expected)");
   }
   else
   {
      TEST_ASSERT(ok, "L-07 undamped diagnostic: slip, shear and station velocity converge in "
                      "dt (every ratio >= 1.6); slip rate reported only");
   }
}

// ------------------------------------------------------------------ main

static int Summary()
{
   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests << " passed, "
                << num_failed << " failed ===\n";
   }
   return num_failed == 0 ? 0 : 1;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   const std::string mode = (argc > 1) ? argv[1] : "";
   if (g_rank == 0)
   {
      std::cout << "=== Nitsche linear slip weakening, inline box, on " << Mpi::WorldSize()
                << " rank(s) " << mode << " ===\n";
   }
   if (mode == "--l06" && argc > 6)
   {
      TestL06(std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]),
              std::atof(argv[6]));
      return Summary();
   }
   if (mode == "--l07" && argc > 6)
   {
      TestL07(std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]),
              std::atof(argv[6]), (argc > 7) ? std::atoi(argv[7]) : 1,
              (argc > 8) ? std::atof(argv[8]) : 0.1);
      return Summary();
   }
   for (int p = 1; p <= 2; p++) { TestDefault(p); }
   return Summary();
}
