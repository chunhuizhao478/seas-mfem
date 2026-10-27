// Nitsche dynamic rupture: hook H-B1 (per-point parameters and |d^n| reach
// the closure; benchmark plan rev 1.3, section 4.2, gate L-10).
//
//   L-10a  A test-only spy closure (the TPV102 aging closure underneath)
//          records st.slip_abs at every box-9 SolveSlipRate call.  In every
//          step it is called exactly once per owned fault DOF, in DOF order,
//          and st.slip_abs equals |d^n| of that DOF bit for bit (d^n read
//          from Slip() before the step).  Run with Prakash-Clifton off and
//          on, since AdvanceFault_ has one CloseDOF call for each.
//   L-10b  The spy sees the stored per-point parameters: st.a and st.Dc are
//          FaultA()(j), FaultDc()(j) bit for bit, and the fields the TPV102
//          setup leaves unset (mu_s, mu_d, cohesion, V_w) are zero.  Box-4
//          and box-10 StateRate calls see slip_abs = 0.
//   L-10c  Abort: a closure whose CheckParameters rejects every point stops
//          the operator setup (--abort-reject, must exit non-zero with its
//          message).
//
// Fixture: the TPV102 inline box of test_nitsche_tpv102_inline.cpp
// (24 tets per hex, fault on Y = 0), p = 1 and 2, with nucleation, 3 s.
//
// Usage: mpirun -np 1|2|4 ./seas_test_nitsche_hook_spy [--abort-reject]
//
// Scope (Chunhui Zhao, 2026-10-16).

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
// (as test_nitsche_pc.cpp: TPV102Box without jitter, TPV102BC,
// TPV102Material, TPV102Setup, TPV102InitialState, SeedPC)

static Mesh *TPV102Box(int nx, int ny, int nz)
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

static NitscheFaultSetup TPV102Setup(const IFaultClosure &cl)
{
   NitscheFaultSetup fs;
   fs.closure = &cl;
   fs.prestress = [](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = TPV102Params::tau_ini;
      s[4] = -TPV102Params::sigma_n;
   };
   fs.nucleation_spatial = [](const real_t x[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      const real_t dx = x[0] - TPV102Params::hypo_along_strike;
      const real_t dz = -x[2] - TPV102Params::hypo_down_dip;
      const real_t v = TPV102Params::nuc_dtau * NucleationSpatial(std::sqrt(dx * dx + dz * dz));
      s[1] = s[3] = v;
   };
   fs.nucleation_temporal = [](real_t t) { return NucleationTemporal(t); };
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

// ------------------------------------------------------------ spy closure

/// Forwards to the TPV102 aging closure and records what the operator hands
/// it.  Test-only.
class SpyClosure : public IFaultClosure
{
public:
   explicit SpyClosure(const IFaultClosure &inner) : inner_(inner) {}

   real_t SolveSlipRate(const real_t X_tan[2], real_t gamma, real_t sigma_n,
                        const FaultPointState &s, real_t V_tan[2]) const override
   {
      solve_.push_back(s);
      return inner_.SolveSlipRate(X_tan, gamma, sigma_n, s, V_tan);
   }
   real_t StateRate(real_t V, const FaultPointState &s) const override
   {
      if (s.slip_abs != 0.0) { ++rate_nonzero_slip_; }
      ++rate_calls_;
      return inner_.StateRate(V, s);
   }
   void SeedEquilibrium(real_t tau0, real_t V_ini, real_t sigma_n,
                        FaultPointState &s) const override
   {
      inner_.SeedEquilibrium(tau0, V_ini, sigma_n, s);
   }

   void Clear() const { solve_.clear(); rate_calls_ = 0; rate_nonzero_slip_ = 0; }
   const std::vector<FaultPointState> &Solves() const { return solve_; }
   long RateCalls() const { return rate_calls_; }
   long RateNonzeroSlip() const { return rate_nonzero_slip_; }

private:
   const IFaultClosure &inner_;
   mutable std::vector<FaultPointState> solve_;
   mutable long rate_calls_ = 0, rate_nonzero_slip_ = 0;
};

/// A closure that rejects every point (L-10c).
class RejectingClosure : public RateStateClosure
{
public:
   using RateStateClosure::RateStateClosure;
   bool CheckParameters(const FaultPointState &) const override { return false; }
};

static bool SameBits(real_t a, real_t b) { return std::memcmp(&a, &b, sizeof(real_t)) == 0; }

// ------------------------------------------------------------- L-10a, b
static void TestSpy(ParMesh &pm, int p, bool pc, real_t T)
{
   DieterichRuinaFriction fr{{TPV102Params::V0, TPV102Params::f0, TPV102Params::b,
                              TPV102Params::Dc}};
   RateStateClosure aging(fr, RateStateClosure::StateLaw::Aging);
   SpyClosure spy(aging);
   NitscheFaultSetup fs = TPV102Setup(spy);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   op.SetTimeStep(dt);
   TPV102InitialState(op);
   if (pc) { SeedPC(op, kVStar, kLength); }
   const int no = op.NumOwnedFaultDOFs();
   const int nsteps = static_cast<int>(std::ceil(T / dt - 1e-9));
   const int n = op.Height();
   Vector U(n), Ud(n), dn;
   U = 0.0;
   Ud = 0.0;
   real_t t = 0.0;
   op.SetTime(t);
   CentralDifferenceSolver ode;
   ode.Init(op);

   long count_bad = 0, slip_bad = 0, par_bad = 0, checked = 0, nonzero = 0;
   long rate_calls = 0, rate_nonzero = 0;
   real_t max_slip = 0.0;
   for (int s = 0; s < nsteps; s++)
   {
      dn = op.Slip();   // d^n
      spy.Clear();
      real_t h = dt;
      ode.Step(U, Ud, t, h);
      const std::vector<FaultPointState> &rec = spy.Solves();
      if (static_cast<int>(rec.size()) != no) { ++count_bad; continue; }
      for (int j = 0; j < no; j++)
      {
         const real_t *d = &dn(3 * j);
         const real_t ref = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
         if (!SameBits(rec[j].slip_abs, ref)) { ++slip_bad; }
         if (ref > 0.0) { ++nonzero; }
         max_slip = std::max(max_slip, ref);
         if (!SameBits(rec[j].a, op.FaultA()(j)) || !SameBits(rec[j].Dc, op.FaultDc()(j))
             || rec[j].mu_s != 0.0 || rec[j].mu_d != 0.0 || rec[j].cohesion != 0.0
             || rec[j].V_w != 0.0)
         {
            ++par_bad;
         }
         ++checked;
      }
      rate_calls += spy.RateCalls();
      rate_nonzero += spy.RateNonzeroSlip();
   }
   long red[7] = {count_bad, slip_bad, par_bad, checked, nonzero, rate_calls, rate_nonzero};
   MPI_Allreduce(MPI_IN_PLACE, red, 7, MPI_LONG, MPI_SUM, pm.GetComm());
   max_slip = RMax(max_slip, pm.GetComm());
   const int no_glob = ISum(no, pm.GetComm());
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  p = " << p << ", PC " << (pc ? "on" : "off")
                << ": " << nsteps << " steps, " << no_glob << " owned fault DOFs; "
                << red[3] << " box-9 calls checked, " << red[4] << " with |d^n| > 0, max |d^n| = "
                << max_slip << " m; call-count mismatches " << red[0] << ", slip_abs mismatches "
                << red[1] << ", parameter mismatches " << red[2] << "; StateRate calls "
                << red[5] << " (" << red[6] << " with slip_abs != 0)\n";
   }
   const std::string tag = "p=" + std::to_string(p) + (pc ? " PC on" : " PC off");
   TEST_ASSERT(red[0] == 0 && red[3] == static_cast<long>(nsteps) * no_glob,
               "L-10a " << tag << ": one box-9 closure call per owned DOF per step");
   TEST_ASSERT(red[1] == 0 && red[4] > 0,
               "L-10a " << tag << ": st.slip_abs == |d^n| bit for bit (nonzero slip seen)");
   TEST_ASSERT(red[2] == 0,
               "L-10b " << tag << ": st carries FaultA, FaultDc bitwise; unset fields zero");
   TEST_ASSERT(red[5] == 2L * nsteps * no_glob && red[6] == 0,
               "L-10b " << tag << ": box-4/10 StateRate calls see slip_abs = 0");
}

// ------------------------------------------------------------------ L-10c
static int AbortReject(ParMesh &pm)
{
   DieterichRuinaFriction fr{{TPV102Params::V0, TPV102Params::f0, TPV102Params::b,
                              TPV102Params::Dc}};
   RejectingClosure cl(fr, RateStateClosure::StateLaw::Aging);
   NitscheFaultSetup fs = TPV102Setup(cl);
   NitscheWaveOperator op(pm, 1, TPV102Material(), TPV102BC(), {}, &fs);   // must abort
   return 0;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   std::unique_ptr<Mesh> m(TPV102Box(4, 4, 2));
   ParMesh pm(MPI_COMM_WORLD, *m);
   if (argc > 1 && std::string(argv[1]) == "--abort-reject") { return AbortReject(pm); }
   if (g_rank == 0)
   {
      std::cout << "=== Hook H-B1 spy (L-10) on " << Mpi::WorldSize() << " rank(s) ===\n";
   }
   for (int p = 1; p <= 2; p++)
   {
      TestSpy(pm, p, false, 3.0);
      TestSpy(pm, p, true, 3.0);
   }
   if (g_rank == 0)
   {
      std::cout << "=== Summary: " << num_passed << "/" << num_tests << " passed, "
                << num_failed << " failed ===\n";
   }
   return num_failed == 0 ? 0 : 1;
}
