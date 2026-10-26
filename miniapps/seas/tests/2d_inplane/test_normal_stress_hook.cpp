// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_normal_stress_hook.cpp — 2D in-plane design rev 3,
// phase D (part 2, section 5.9): the Nitsche operator's normal-stress hook
// (NitscheWaveOperator::SetNormalStressModifier) and thermal pressurization
// behind it (fault/thermal_pressurization_modifier.hpp).  The problem is the
// in-plane slip-weakening rupture of test_2d_nitsche.cpp N5: the 24 x 24 grid
// of [-30, 30]^2 km (each cell split into 4 triangles), the fault on y = 0,
// |x| < 15 km, absorbing walls, TPV205 friction, tau_0 = 70 MPa, sigma_0 =
// 120 MPa, the nucleation 15 MPa F(|x|) G(t) over 9 km, up = [0, 0, 1].
//
//   H1  off is off: a modifier that takes off 0 Pa, and thermal
//       pressurization with Lambda = 0, give the run without a modifier bit
//       for bit (slip, slip rate, state and sigma_eff at every owned DOF and
//       step), with Kelvin-Voigt and Prakash-Clifton on.  p = 1, 2.
//   H2  what the law sees: at every owned DOF, CloseDOF with a modifier that
//       takes off R_j gives the slip rate of CloseDOF without it at the
//       trial traction moved by R_j along the normal (sigma_eff - R_j) to
//       1e-12 of the largest slip rate, and returns the same sigma_eff as
//       without the modifier, bit for bit.
//   H3  what the modifier gets: one Advance per step with the operator's dt;
//       V[j] = |V^{n+1}| and tau[j] = |X + gamma V| (the station traction)
//       bit for bit, and at every slipping DOF tau[j] = mu(|d^n|) (sigma_eff
//       - R_j), the strength the law used, to 1e-12.
//   H4  thermal pressurization on the rupture (TPV105's parameters, Lambda
//       0.1 MPa/K): the run is finite; T and p rise only where the fault
//       slipped, and the effective normal stress the law sees falls there;
//       at every DOF the heat content equals the heat put in (sum of tau V
//       dt times the grid's source fraction) to 1e-12; the slip is larger
//       than without thermal pressurization.
//
// Refusals (--refuse n, each must abort with its message):
//   1  a modifier that covers another number of points than the owned DOFs
//   2  a modifier that takes off more than sigma_eff
//
// Usage: mpirun -np 1 ./seas_test_normal_stress_hook [--refuse n]
//
// Scope (Chunhui Zhao, 2026-10-26).

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../friction/slip_weakening_closure.hpp"
#include "../../fault/thermal_pressurization_modifier.hpp"
#include "../../config/tpv102_params.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kL = 30.0e3;
static const real_t kMuS = 0.677, kMuD = 0.525, kDc = 0.40;
static const real_t kTau0 = 70.0e6, kSigma0 = 120.0e6;
static const real_t kNucR = 9.0e3, kNucDtau = 15.0e6;
static const real_t kLf = 15.0e3;
static const int kNf = 24;

/// [-L, L]^2 in kNf x kNf cells split into 4 triangles through the centre;
/// walls attr 5, the grid edges on y = 0 with |x| < kLf attr 3 (the mesh of
/// test_2d_nitsche.cpp N5).
static Mesh *FaultBox()
{
   const int nx = kNf, ny = kNf;
   const real_t h = 2 * kL / kNf, x0 = -kL, y0 = -kL;
   const int ngv = (nx + 1) * (ny + 1);
   Mesh *m = new Mesh(2, ngv + nx * ny, 4 * nx * ny, 2 * (nx + ny), 2);
   auto gv = [&](int i, int j) { return i + (nx + 1) * j; };
   auto cv = [&](int i, int j) { return ngv + i + nx * j; };
   for (int j = 0; j <= ny; j++)
      for (int i = 0; i <= nx; i++) { m->AddVertex(x0 + i * h, y0 + j * h); }
   for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) { m->AddVertex(x0 + (i + 0.5) * h, y0 + (j + 0.5) * h); }
   for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
      {
         const int c = cv(i, j);
         m->AddTriangle(gv(i, j), gv(i + 1, j), c, 1);
         m->AddTriangle(gv(i + 1, j), gv(i + 1, j + 1), c, 1);
         m->AddTriangle(gv(i + 1, j + 1), gv(i, j + 1), c, 1);
         m->AddTriangle(gv(i, j + 1), gv(i, j), c, 1);
      }
   for (int i = 0; i < nx; i++)
   {
      m->AddBdrSegment(gv(i, 0), gv(i + 1, 0), 5);
      m->AddBdrSegment(gv(i + 1, ny), gv(i, ny), 5);
   }
   for (int j = 0; j < ny; j++)
   {
      m->AddBdrSegment(gv(nx, j), gv(nx, j + 1), 5);
      m->AddBdrSegment(gv(0, j + 1), gv(0, j), 5);
   }
   const int jf = ny / 2;   // y = 0
   for (int i = 0; i < nx; i++)
   {
      if (std::abs(x0 + (i + 0.5) * h) < kLf) { m->AddBdrSegment(gv(i, jf), gv(i + 1, jf), 3); }
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

static real_t NucF(real_t r)
{
   if (r >= kNucR) { return 0.0; }
   const real_t r2 = r * r, R2 = kNucR * kNucR;
   return std::exp(r2 / (r2 - R2));
}

static NitscheFaultSetup FaultSetup(const IFaultClosure &cl)
{
   NitscheFaultSetup fs;
   fs.closure = &cl;
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   for (int d = 0; d < 3; d++)
   {
      fs.orientation.ref_normal[d] = ref[d];
      fs.orientation.up[d] = up[d];
      fs.up[d] = up[d];
   }
   fs.orientation.fault_attrs = {3};
   fs.prestress = [](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = kTau0;
      s[4] = -kSigma0;
   };
   fs.nucleation_spatial = [](const real_t x[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = kNucDtau * NucF(std::abs(x[0]));
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

/// A modifier that takes off red[j] (or forwards to `inner`) and records
/// what Advance receives.
class SpyModifier : public IFaultNormalStressModifier
{
public:
   SpyModifier(int n, IFaultNormalStressModifier *inner = nullptr)
      : n_(n), inner_(inner), red(n, 0.0), heat_in(n, 0.0) {}
   int NumPoints() const override { return n_; }
   real_t NormalStressReduction(int j) const override
   {
      return inner_ ? inner_->NormalStressReduction(j) : red[j];
   }
   void Advance(real_t dt, const real_t *t, const real_t *v) override
   {
      calls++;
      last_dt = dt;
      tau.assign(t, t + n_);
      V.assign(v, v + n_);
      for (int j = 0; j < n_; j++) { heat_in[j] += t[j] * v[j] * dt; }
      if (inner_) { inner_->Advance(dt, t, v); }
   }
   int n_;
   IFaultNormalStressModifier *inner_;
   std::vector<real_t> red, heat_in, tau, V;
   long calls = 0;
   real_t last_dt = 0.0;
};

static TPParameters TPV105Params(real_t Lambda)
{
   TPParameters P;
   P.alpha_th = 1.0e-6;
   P.alpha_hy = 1.0e-4;
   P.rho_c = 2.7e6;
   P.Lambda = Lambda;
   P.h = 0.01;
   P.T_ini = 483.15;
   P.p_ini = 80.0e6;
   return P;
}

enum class Mode { None, Zero, Reduce, TP };

struct RunConfig
{
   int p = 1;
   real_t tfinal = 2.0;
   bool kvpc = false;        ///< Kelvin-Voigt (eta_bar 0.1) and Prakash-Clifton (v* 1 m/s, L 0.01 m)
   Mode mode = Mode::None;
   real_t Lambda = 0.0;      ///< Mode::TP
   bool checks = false;      ///< H2 and H3 (Mode::Reduce)
   int refuse = 0;
};

struct RunResult
{
   bool finite = true;
   long steps = 0;
   std::vector<real_t> trace;   ///< per step and owned DOF: d (3), V (3), psi, sigma_eff
   // H2, H3
   long adv_calls = 0;
   bool adv_dt_ok = true;
   long adv_bad = 0;            ///< tau or V not the operator's, bit for bit
   real_t str_err = 0.0;        ///< max |tau - mu (sigma_eff - R)| / tau at slipping DOFs
   long slipping = 0;
   real_t h2_err = -1.0;        ///< max |V_hook - V_moved| / scale
   long h2_sig_bad = 0;
   // H4
   real_t dmax = 0.0, Tmax = 0.0, pmax = 0.0, min_law = 1e300, heat_err = 0.0;
   long tp_quiet_bad = 0;       ///< T or p changed where the fault never slipped
   long law_bad = 0;            ///< sigma the law sees not below sigma_eff where p rose
};

static RunResult Run(const RunConfig &rc)
{
   RunResult R;
   std::unique_ptr<Mesh> m(FaultBox());
   ParMesh pm(MPI_COMM_WORLD, *m);
   m.reset();
   SlipWeakeningClosure cl;
   NitscheFaultSetup fs = FaultSetup(cl);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.absorbing_attrs = {5};
   const MaterialField mat = MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu,
                                                         TPV102Params::rho);
   NitscheWaveOperator op(pm, rc.p, mat, bc, {}, &fs);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, rc.p, 0.5 * NitscheStableConstant(2, rc.p));
   op.SetTimeStep(dt);
   if (rc.kvpc) { op.SetKelvinVoigt(0.1, NitscheWaveOperator::KVScale::Node, -1.0); }
   const int no = op.NumOwnedFaultDOFs();
   Vector sig0(no);
   {
      Vector d(3 * no), V(3 * no), psi(no), S;
      d = 0.0;
      V = 0.0;
      op.BackgroundTraction(0.0, S);
      for (int j = 0; j < no; j++)
      {
         FaultPointState st = op.FaultParams(j);
         op.Closure().SeedEquilibrium(kTau0, 0.0, kSigma0, st);
         psi(j) = st.psi;
         const real_t *f = &op.FaultFrames()(9 * j);
         real_t X0, X1;
         fault_frame::ToLocal(&S(3 * j), 0.0, f, f + 3, f + 6, X0, X1, sig0(j));
      }
      op.SetFaultState(d, V, psi);
   }
   if (rc.kvpc)
   {
      op.SetPrakashClifton(1.0, 0.01);
      op.SetPrakashCliftonState(sig0);
   }
   // the modifier
   std::unique_ptr<ThermalPressurization> tp;
   std::unique_ptr<ThermalPressurizationModifier> tpm;
   std::unique_ptr<SpyModifier> spy;
   if (rc.mode == Mode::TP)
   {
      TPGridSpec g;
      g.dz0 = 0.5e-3;
      g.growth = 1.05;
      g.z_max = 0.4;
      tp = std::make_unique<ThermalPressurization>(no, g);
      tp->SetParameters(TPV105Params(rc.Lambda));
      tpm = std::make_unique<ThermalPressurizationModifier>(*tp);
      spy = std::make_unique<SpyModifier>(no, tpm.get());
   }
   else if (rc.mode != Mode::None)
   {
      spy = std::make_unique<SpyModifier>(no + (rc.refuse == 1 ? 1 : 0));
      if (rc.mode == Mode::Reduce)
      {
         // 4 MPa plus a per-DOF part, so that every DOF has its own value
         for (int j = 0; j < no; j++) { spy->red[j] = 4.0e6 + 1.0e3 * j; }
      }
      if (rc.refuse == 2) { for (int j = 0; j < no; j++) { spy->red[j] = 2.0 * kSigma0; } }
   }
   if (spy) { op.SetNormalStressModifier(spy.get()); }

   const int n = op.Height();
   Vector U(n), Ud(n), dprev;
   U = 0.0;
   Ud = 0.0;
   real_t t = 0.0;
   op.SetTime(t);
   CentralDifferenceSolver ode;
   ode.Init(op);
   const long nsteps = std::lround(std::ceil(rc.tfinal / dt - 1e-9));
   std::vector<char> slipped(no, 0);
   for (long s = 0; s < nsteps; s++)
   {
      dprev = op.Slip();
      const long calls0 = spy ? spy->calls : 0;
      real_t h_s = dt;
      ode.Step(U, Ud, t, h_s);
      const NitscheStepRecord &rec = op.LastStep();
      if (!std::isfinite(rec.max_V) || !std::isfinite(rec.min_sigma_eff)) { R.finite = false; break; }
      const Vector &d = op.Slip(), &V = op.SlipRate(), &psi = op.State();
      const Vector &X = op.LastShear(), &sg = op.LastSigmaEff();
      for (int j = 0; j < no; j++)
      {
         for (int c = 0; c < 3; c++) { R.trace.push_back(d(3 * j + c)); }
         for (int c = 0; c < 3; c++) { R.trace.push_back(V(3 * j + c)); }
         R.trace.push_back(psi(j));
         R.trace.push_back(sg(j));
         const real_t *v = &V(3 * j);
         if (v[0] != 0.0 || v[1] != 0.0 || v[2] != 0.0) { slipped[j] = 1; }
      }
      if (spy && rc.checks)
      {
         // H3: once per step, the operator's dt, its |V| and station traction
         R.adv_calls += spy->calls - calls0;
         R.adv_dt_ok = R.adv_dt_ok && spy->last_dt == dt;
         for (int j = 0; j < no; j++)
         {
            const real_t *v = &V(3 * j), *dv = &op.FaultFrames()(9 * j + 3),
                          *sv = &op.FaultFrames()(9 * j + 6);
            const real_t Vabs = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            const real_t g = op.ClosureGamma(j);
            const real_t tau = std::hypot(X(2 * j) + g * fault_frame::Dot(v, dv),
                                          X(2 * j + 1) + g * fault_frame::Dot(v, sv));
            if (spy->V[j] != Vabs || spy->tau[j] != tau) { R.adv_bad++; }
            if (Vabs > 0.0)
            {
               const real_t *dn = &dprev(3 * j);
               const real_t D = std::sqrt(dn[0] * dn[0] + dn[1] * dn[1] + dn[2] * dn[2]);
               const real_t mu = kMuS - (kMuS - kMuD) * std::min(D / kDc, real_t(1.0));
               const real_t sig_law = (rc.kvpc ? op.PrakashCliftonSigma()(j) : sg(j)) - spy->red[j];
               R.str_err = std::max(R.str_err, std::abs(tau - mu * sig_law) / tau);
               R.slipping++;
            }
         }
      }
   }
   R.steps = op.LastStep().step;
   for (real_t q : R.trace) { if (!std::isfinite(q)) { R.finite = false; } }

   if (rc.checks && spy)
   {
      // H2: CloseDOF with the modifier against CloseDOF without it at the
      // trial traction moved by R_j along the normal; the state of the last step
      Vector S;
      op.BackgroundTraction(t, S);
      const Vector &T = op.LastTrialTraction();
      std::vector<real_t> Va(3 * no), Vb(3 * no);
      real_t scale = 1e-3;
      for (int j = 0; j < no; j++)
      {
         const real_t *nv = &op.FaultFrames()(9 * j);
         real_t TS[3], TSm[3];
         for (int c = 0; c < 3; c++)
         {
            TS[c] = T(3 * j + c) + S(3 * j + c);
            TSm[c] = TS[c] + spy->red[j] * nv[c];
         }
         const real_t *dn = &op.Slip()(3 * j);
         const real_t D = std::sqrt(dn[0] * dn[0] + dn[1] * dn[1] + dn[2] * dn[2]);
         real_t sa, sb, sc, X[2];
         op.SetNormalStressModifier(spy.get());
         op.CloseDOF(j, TS, op.State()(j), &Va[3 * j], sa, X, nullptr, 0.0, 1.0, D);
         op.SetNormalStressModifier(nullptr);
         op.CloseDOF(j, TSm, op.State()(j), &Vb[3 * j], sb, X, nullptr, 0.0, 1.0, D);
         real_t Vc[3];
         op.CloseDOF(j, TS, op.State()(j), Vc, sc, X, nullptr, 0.0, 1.0, D);
         if (sa != sc) { R.h2_sig_bad++; }
         scale = std::max(scale, std::hypot(Vb[3 * j], Vb[3 * j + 1], Vb[3 * j + 2]));
      }
      real_t e = 0.0;
      for (int i = 0; i < 3 * no; i++) { e = std::max(e, std::abs(Va[i] - Vb[i])); }
      R.h2_err = e / scale;
      op.SetNormalStressModifier(spy.get());
   }
   if (tp)
   {
      // H4
      const Vector &sg = op.LastSigmaEff();
      for (int j = 0; j < no; j++)
      {
         const real_t *d = &op.Slip()(3 * j);
         R.dmax = std::max(R.dmax, std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
         R.Tmax = std::max(R.Tmax, tp->Temperature(j));
         R.pmax = std::max(R.pmax, tp->PorePressure(j));
         const real_t law = sg(j) - tpm->NormalStressReduction(j);
         R.min_law = std::min(R.min_law, law);
         if (!slipped[j] && (tp->TemperatureChange(j) != 0.0 || tp->PorePressureChange(j) != 0.0))
         {
            R.tp_quiet_bad++;
         }
         if (tp->PorePressureChange(j) > 0.0 && !(law < sg(j))) { R.law_bad++; }
         const real_t in = spy->heat_in[j] * tp->SourceFraction(j);
         if (in > 0.0) { R.heat_err = std::max(R.heat_err, std::abs(tp->HeatContent(j) - in) / in); }
      }
   }
   else
   {
      for (int j = 0; j < no; j++)
      {
         const real_t *d = &op.Slip()(3 * j);
         R.dmax = std::max(R.dmax, std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
      }
   }
   return R;
}

static real_t TraceDiff(const RunResult &a, const RunResult &b)
{
   if (a.trace.size() != b.trace.size()) { return -1.0; }
   real_t e = 0.0;
   for (size_t i = 0; i < a.trace.size(); i++) { e = std::max(e, std::abs(a.trace[i] - b.trace[i])); }
   return e;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   int refuse = 0;
   for (int i = 1; i + 1 < argc; i++) { if (!std::strcmp(argv[i], "--refuse")) { refuse = std::atoi(argv[i + 1]); } }
   MFEM_VERIFY(Mpi::WorldSize() == 1, "test_normal_stress_hook runs on 1 rank");
   std::cout << std::setprecision(6);
   if (refuse)
   {
      RunConfig rc;
      rc.tfinal = 0.5;
      rc.mode = Mode::Zero;
      rc.refuse = refuse;
      Run(rc);
      std::cout << "refusal " << refuse << " did not abort\n";
      return 0;
   }
   std::cout << "test_normal_stress_hook (2D in-plane design rev 3, phase D part 2)\n";

   for (int p = 1; p <= 2; p++)
   {
      std::cout << "\nH1: off is off, p = " << p << " (Kelvin-Voigt and Prakash-Clifton on)\n";
      RunConfig rc;
      rc.p = p;
      rc.kvpc = true;
      const RunResult base = Run(rc);
      rc.mode = Mode::Zero;
      const RunResult zero = Run(rc);
      rc.mode = Mode::TP;
      rc.Lambda = 0.0;
      const RunResult lam0 = Run(rc);
      std::cout << "  " << base.steps << " steps, max slip " << base.dmax << " m\n";
      TEST_ASSERT(base.finite && base.steps > 0 && base.dmax > 1.0,
                  "H1 p = " << p << ": the rupture runs and slips more than 1 m");
      TEST_ASSERT(TraceDiff(base, zero) == 0.0, "H1 p = " << p << ": a modifier that takes "
                  "off 0 Pa gives the run without one bit for bit (" << base.trace.size()
                  << " values)");
      TEST_ASSERT(TraceDiff(base, lam0) == 0.0, "H1 p = " << p << ": thermal pressurization "
                  "with Lambda = 0 gives the run without it bit for bit (T max " << lam0.Tmax
                  << " K, p max " << lam0.pmax << " Pa)");
      TEST_ASSERT(lam0.Tmax > 483.15 && lam0.pmax == 80.0e6, "H1 p = " << p << ": with Lambda "
                  "= 0 the fault heats and p stays exactly p_ini");
   }

   for (int p = 1; p <= 2; p++)
   {
      for (int kvpc = 0; kvpc <= 1; kvpc++)
      {
         std::cout << "\nH2, H3: a modifier that takes off 4 MPa + 1 kPa j, p = " << p
                   << (kvpc ? ", Kelvin-Voigt and Prakash-Clifton on" : "") << "\n";
         RunConfig rc;
         rc.p = p;
         rc.kvpc = kvpc != 0;
         rc.mode = Mode::Reduce;
         rc.checks = true;
         const RunResult r = Run(rc);
         std::cout << "  " << r.steps << " steps, max slip " << r.dmax << " m; CloseDOF "
                   << "difference " << r.h2_err << " of the largest slip rate; strength error "
                   << r.str_err << " over " << r.slipping << " slipping DOF-steps\n";
         const std::string tag = "p = " + std::to_string(p) + (kvpc ? ", KV + PC" : "");
         TEST_ASSERT(r.finite && r.dmax > 1.0, "H2 " << tag << ": the rupture runs");
         TEST_ASSERT(r.h2_err >= 0.0 && r.h2_err < 1e-12 && r.h2_sig_bad == 0,
                     "H2 " << tag << ": the law sees sigma_eff - R_j (slip rate to "
                     << r.h2_err << "); sigma_eff on return unchanged");
         TEST_ASSERT(r.adv_calls == r.steps && r.adv_dt_ok, "H3 " << tag << ": one Advance per "
                     "step (" << r.adv_calls << " for " << r.steps << "), with the operator's dt");
         TEST_ASSERT(r.adv_bad == 0, "H3 " << tag << ": tau = |X + gamma V| and V = |V^{n+1}| "
                     "bit for bit");
         TEST_ASSERT(r.slipping > 0 && r.str_err < 1e-12, "H3 " << tag << ": tau = mu (sigma - "
                     "R_j) at slipping DOFs (" << r.str_err << ")");
      }
   }

   for (int p = 1; p <= 2; p++)
   {
      std::cout << "\nH4: thermal pressurization on the rupture, TPV105 parameters, p = " << p << "\n";
      RunConfig rc;
      rc.p = p;
      rc.kvpc = true;
      const RunResult base = Run(rc);
      rc.mode = Mode::TP;
      rc.Lambda = 0.1e6;
      const RunResult r = Run(rc);
      std::cout << "  max slip " << r.dmax << " m (without thermal pressurization " << base.dmax
                << " m); max T " << r.Tmax << " K, max p " << r.pmax << " Pa, min sigma the law "
                "sees " << r.min_law << " Pa; heat balance " << r.heat_err << "\n";
      TEST_ASSERT(r.finite, "H4 p = " << p << ": finite");
      TEST_ASSERT(r.Tmax > 483.15 && r.pmax > 80.0e6 && r.tp_quiet_bad == 0 && r.law_bad == 0,
                  "H4 p = " << p << ": T and p rise only where the fault slipped, and the law "
                  "sees less than sigma_eff there");
      TEST_ASSERT(r.heat_err < 1e-12, "H4 p = " << p << ": heat content = heat put in at every "
                  "DOF (" << r.heat_err << ")");
      TEST_ASSERT(r.dmax > base.dmax, "H4 p = " << p << ": more slip with thermal "
                  "pressurization");
   }

   std::cout << "\ntest_normal_stress_hook: " << (num_tests - num_failed) << " / " << num_tests
             << " passed\n";
   return num_failed ? 1 : 0;
}
