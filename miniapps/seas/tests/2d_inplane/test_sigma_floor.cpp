// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_sigma_floor.cpp — the opt-in floor of the normal
// stress the friction law sees ([cycle].sigma_n_floor_pa; Chunhui Zhao,
// 2026-10-30: the UCERF 2-D map-view run, sigma_n' < 0 at trace bends, "Opt-in
// sigma_n' floor in the friction law", 1 MPa in both regimes).
//
//   K1  the QD kernels (friction/device_rate_state.hpp, on the host):
//       RSFaultRHS and RSFaultSlipRate with sigma_floor = F give, bit for bit,
//       what they give without the floor at sigma_n_eff = max(sigma_n_eff, F);
//       nodes at or above F are unchanged bit for bit; kRSSigmaNFloored is set
//       and kRSSigmaNNonPositive is not when a node is below F; without the
//       floor the same inputs set kRSSigmaNNonPositive and not
//       kRSSigmaNFloored; a NaN normal stress is still kRSBadInput with the
//       floor on.
//   D1  the dynamic operator (NitscheWaveOperator::CloseDOF, the in-plane
//       slip-weakening fault of test_normal_stress_hook.cpp, p = 1, 2): at
//       trial tractions moved along the normal to sigma_t from 2F down to
//       -2F, the operator with sigma_floor = F gives the slip rate of the
//       operator without it at the traction moved to F (1e-12 of the largest
//       slip rate) where sigma_t < F, and the same slip rate and sigma_eff bit
//       for bit where sigma_t >= F; sigma_eff on return is the true sigma_t;
//       FloorHits / FlooredDOFs count the floored DOF-calls and DOFs, and
//       ResetFloorCount clears them.
//
// Refusals (--refuse n, each must abort with its message):
//   1  NitscheFaultSetup sigma_floor < 0
//   2  sigma_floor > 0 with sigma_min > 0
//   3  CloseDOF with a non-finite sigma_n_eff and the floor on
//   4  without the floor, CloseDOF below sigma_min still aborts (unchanged)
// The deck keys are checked through seas_fd_cycle_driver (test-sigma-floor in
// 2d_inplane.mk).
//
// Usage: mpirun -np 1 ./seas_test_sigma_floor [--refuse n]

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../friction/slip_weakening_closure.hpp"
#include "../../friction/device_rate_state.hpp"
#include "../../config/tpv102_params.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
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

static const real_t kF = 1.0e6;   // the floor [Pa]

// ---------------------------------------------------------------- K1

/// Per-node inputs of the BP5 fault kernels (state [slip, slip, psi]).
struct KernelInputs
{
   Vector T, X, TN, TP, SN, A, ETA, DC;
   void Resize(int n)
   {
      T.SetSize(2 * n); X.SetSize(3 * n); TN.SetSize(n); TP.SetSize(2 * n);
      SN.SetSize(n); A.SetSize(n); ETA.SetSize(n); DC.SetSize(n);
   }
   RSFaultInputs In(real_t floor, bool normal) const
   {
      RSFaultInputs in;
      in.traction = &T;
      in.state = &X;
      in.normal = normal ? &TN : nullptr;
      in.tau_pre = &TP;
      in.sigma_n_dof = &SN;
      in.sigma_n = 0.0;
      in.a = &A;
      in.eta = &ETA;
      in.sigma_floor = floor;
      return in;
   }
};

static bool SameBits(const Vector &a, const Vector &b, const std::vector<int> &idx, int stride)
{
   for (int i : idx)
   {
      for (int c = 0; c < stride; c++)
      {
         const real_t x = a(stride * i + c), y = b(stride * i + c);
         if (!(x == y || (std::isnan(x) && std::isnan(y)))) { return false; }
      }
   }
   return true;
}

static void TestKernels()
{
   std::cout << "\nK1: the QD fault kernels with sigma_floor = " << kF << " Pa\n";
   // sigma_n_eff = SN + TN: above, at, below the floor, zero, in tension
   const std::vector<real_t> sig = {80.0e6, kF, 0.5e6, 1.0e-3, 0.0, -0.08e6, -5.0e6, 45.4e6};
   const int n = static_cast<int>(sig.size());
   KernelInputs k;
   k.Resize(n);
   for (int i = 0; i < n; i++)
   {
      k.TN(i) = -20.0e6 + 3.0e6 * i;       // elastic normal traction
      k.SN(i) = sig[i] - k.TN(i);           // pre-stress, so SN + TN ~ sig
      k.TP(2 * i) = 0.1e6 * i;
      k.TP(2 * i + 1) = -(sig[i] > 0.0 ? 0.6 * sig[i] : 2.0e6) - 0.2e6 * i;
      k.T(2 * i) = 1.0e3 * i;
      k.T(2 * i + 1) = -3.0e3;
      k.X(3 * i) = 0.1 * i;
      k.X(3 * i + 1) = 0.2;
      k.X(3 * i + 2) = 0.55 + 0.02 * i;     // psi
      k.A(i) = (i % 2) ? 0.04 : 0.004;
      k.ETA(i) = 4.62e6;
      k.DC(i) = 0.47;
   }
   std::vector<int> above, below;
   for (int i = 0; i < n; i++)
   {
      ((k.SN(i) + k.TN(i)) < kF ? below : above).push_back(i);
   }
   // the same nodes without the floor, at max(SN + TN, F), with no normal
   KernelInputs m = k;
   for (int i = 0; i < n; i++) { m.SN(i) = std::max(k.SN(i) + k.TN(i), kF); }

   RSDeviceParams prm;
   prm.V0 = 1.0e-6;
   prm.aging_b = 0.03;
   prm.aging_V0 = 1.0e-6;
   prm.aging_f0 = 0.6;
   Array<RSReduceValue> ws;
   Vector Ra(3 * n), Va(2 * n), Rb(3 * n), Vb(2 * n), Rc(3 * n), Vc(2 * n);
   const RSReduceValue a = RSFaultRHS(n, false, prm, k.In(kF, true), k.DC, Ra, Va, ws);
   const RSReduceValue b = RSFaultRHS(n, false, prm, m.In(0.0, false), m.DC, Rb, Vb, ws);
   const RSReduceValue c = RSFaultRHS(n, false, prm, k.In(0.0, true), k.DC, Rc, Vc, ws);
   std::vector<int> all(n);
   for (int i = 0; i < n; i++) { all[i] = i; }
   std::cout << "  " << below.size() << " of " << n << " nodes below the floor; status with the "
             << "floor " << a.second << ", at max(sigma, F) without it " << b.second
             << ", the same inputs without it " << c.second << "\n";
   TEST_ASSERT(!below.empty() && !above.empty(), "K1: nodes on both sides of the floor");
   TEST_ASSERT(SameBits(Ra, Rb, all, 3) && SameBits(Va, Vb, all, 2) && a.first == b.first,
               "K1 RSFaultRHS: the floor = the law at max(sigma_n_eff, F), bit for bit "
               "(rate, slip rate, max V)");
   TEST_ASSERT(SameBits(Ra, Rc, above, 3) && SameBits(Va, Vc, above, 2),
               "K1 RSFaultRHS: nodes at or above the floor unchanged, bit for bit");
   TEST_ASSERT((a.second & kRSSigmaNFloored) && !(a.second & kRSSigmaNNonPositive) &&
               b.second == 0, "K1 RSFaultRHS: kRSSigmaNFloored set, kRSSigmaNNonPositive not");
   TEST_ASSERT((c.second & kRSSigmaNNonPositive) && !(c.second & kRSSigmaNFloored),
               "K1 RSFaultRHS without the floor: kRSSigmaNNonPositive as before, no "
               "kRSSigmaNFloored");

   Vector Sa(2 * n), Sb(2 * n), Sc(2 * n);
   const int sa = RSFaultSlipRate(n, false, prm.V0, k.In(kF, true), Sa, ws);
   const int sb = RSFaultSlipRate(n, false, prm.V0, m.In(0.0, false), Sb, ws);
   const int sc = RSFaultSlipRate(n, false, prm.V0, k.In(0.0, true), Sc, ws);
   TEST_ASSERT(SameBits(Sa, Sb, all, 2) && SameBits(Sa, Va, all, 2),
               "K1 RSFaultSlipRate: the floor = the law at max(sigma_n_eff, F), bit for bit, "
               "and = RSFaultRHS's slip rate");
   TEST_ASSERT(SameBits(Sa, Sc, above, 2) && (sa & kRSSigmaNFloored) && sb == 0 &&
               !(sc & kRSSigmaNFloored), "K1 RSFaultSlipRate: nodes at or above unchanged; "
               "kRSSigmaNFloored only with the floor");

   KernelInputs q = k;
   q.TN(2) = std::numeric_limits<real_t>::quiet_NaN();
   Vector Rq(3 * n), Vq(2 * n);
   const RSReduceValue r = RSFaultRHS(n, false, prm, q.In(kF, true), q.DC, Rq, Vq, ws);
   TEST_ASSERT((r.second & kRSBadInput) && std::isnan(Rq(3 * 2)),
               "K1: a NaN normal stress is still kRSBadInput with the floor on");
}

// ---------------------------------------------------------------- D1

static const real_t kL = 30.0e3;
static const real_t kMuS = 0.677, kMuD = 0.525, kDc = 0.40;
static const real_t kTau0 = 70.0e6, kSigma0 = 120.0e6;
static const real_t kLf = 15.0e3;
static const int kNf = 24;

/// [-L, L]^2 in kNf x kNf cells split into 4 triangles through the centre;
/// walls attr 5, the grid edges on y = 0 with |x| < kLf attr 3 (the mesh of
/// test_normal_stress_hook.cpp).
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

static NitscheFaultSetup FaultSetup(const IFaultClosure &cl, real_t floor, real_t sigma_min)
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
   fs.parameters = [](const real_t[3], FaultPointState &st) {
      st.mu_s = kMuS;
      st.mu_d = kMuD;
      st.Dc = kDc;
      st.cohesion = 0.0;
   };
   fs.sigma_min = sigma_min;
   fs.sigma_floor = floor;
   return fs;
}

struct Ops
{
   std::unique_ptr<ParMesh> pm;
   SlipWeakeningClosure cl;
   std::unique_ptr<NitscheWaveOperator> off, on;
   BoundaryConfig bc;
   std::unique_ptr<MaterialField> mat;
   Ops(int p, real_t floor, real_t sigma_min = 0.0)
   {
      std::unique_ptr<Mesh> m(FaultBox());
      pm = std::make_unique<ParMesh>(MPI_COMM_WORLD, *m);
      bc.fault_attr = 3;
      bc.absorbing_attrs = {5};
      mat = std::make_unique<MaterialField>(MaterialField::MakeConstant(
         TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho));
      NitscheFaultSetup f0 = FaultSetup(cl, 0.0, 0.0), f1 = FaultSetup(cl, floor, sigma_min);
      off = std::make_unique<NitscheWaveOperator>(*pm, p, *mat, bc, DomainConfig{}, &f0);
      on = std::make_unique<NitscheWaveOperator>(*pm, p, *mat, bc, DomainConfig{}, &f1);
      for (NitscheWaveOperator *op : {off.get(), on.get()})
      {
         op->SetTimeStep(ComputeMaxDt(op->Sample(), *pm, p, 0.5 * NitscheStableConstant(2, p)));
      }
   }
};

/// The background traction of DOF j moved along the normal to sigma_n_eff = s.
static void MovedTraction(const NitscheWaveOperator &op, const Vector &S, int j, real_t s,
                          real_t TS[3])
{
   const real_t *nv = &op.FaultFrames()(9 * j);
   real_t sig = 0.0;
   for (int c = 0; c < 3; c++) { sig -= S(3 * j + c) * nv[c]; }
   for (int c = 0; c < 3; c++) { TS[c] = S(3 * j + c) + (sig - s) * nv[c]; }
}

static void TestDynamic(int p)
{
   std::cout << "\nD1: CloseDOF with sigma_floor = " << kF << " Pa, p = " << p << "\n";
   Ops o(p, kF);
   NitscheWaveOperator &on = *o.on, &off = *o.off;
   const int no = on.NumOwnedFaultDOFs();
   Vector S;
   on.BackgroundTraction(0.0, S);
   long below = 0, above_bad = 0, sig_bad = 0;
   real_t err = 0.0, scale = 1e-3, worst_sig = 0.0;
   for (int j = 0; j < no; j++)
   {
      const real_t u = no > 1 ? real_t(j) / (no - 1) : 0.0;
      const real_t st = kF * (2.0 - 4.0 * u);   // 2F .. -2F
      real_t TS[3], Va[3], Vb[3], X[2], sa, sb;
      MovedTraction(on, S, j, st, TS);
      on.CloseDOF(j, TS, 0.0, Va, sa, X);
      worst_sig = std::max(worst_sig, std::abs(sa - st));
      if (sa >= kF)
      {
         off.CloseDOF(j, TS, 0.0, Vb, sb, X);
         for (int c = 0; c < 3; c++) { if (Va[c] != Vb[c]) { above_bad++; } }
         if (sa != sb) { sig_bad++; }
      }
      else
      {
         below++;
         real_t TSf[3];
         MovedTraction(off, S, j, kF, TSf);
         off.CloseDOF(j, TSf, 0.0, Vb, sb, X);
         for (int c = 0; c < 3; c++) { err = std::max(err, std::abs(Va[c] - Vb[c])); }
         scale = std::max(scale, std::hypot(Vb[0], Vb[1], Vb[2]));
         if (!(sa < kF)) { sig_bad++; }
      }
   }
   const long hits = on.FloorHits();
   const int dofs = on.FlooredDOFs();
   std::cout << "  " << no << " owned fault DOFs, " << below << " below the floor; slip-rate "
             << "difference " << err / scale << " of the largest; |sigma_eff - sigma_t| <= "
             << worst_sig << " Pa; FloorHits " << hits << ", FlooredDOFs " << dofs << "\n";
   TEST_ASSERT(below > 0 && below < no, "D1 p = " << p << ": DOFs on both sides of the floor");
   TEST_ASSERT(err / scale < 1e-12, "D1 p = " << p << ": below the floor the law sees F (slip "
               "rate of the operator without the floor at F, to " << err / scale << ")");
   TEST_ASSERT(above_bad == 0, "D1 p = " << p << ": at or above the floor, the slip rate of the "
               "operator without the floor, bit for bit");
   TEST_ASSERT(sig_bad == 0 && worst_sig < 1e-3, "D1 p = " << p << ": sigma_eff on return is "
               "the true one (below F where floored; the same bits where not)");
   TEST_ASSERT(hits == below && dofs == below && off.FloorHits() == 0,
               "D1 p = " << p << ": FloorHits = FlooredDOFs = the floored DOFs; 0 without "
               "the floor");
   {
      real_t TS[3], V[3], X[2], s;
      MovedTraction(on, S, no - 1, -kF, TS);
      on.CloseDOF(no - 1, TS, 0.0, V, s, X);
   }
   TEST_ASSERT(on.FloorHits() == below + 1 && on.FlooredDOFs() == dofs,
               "D1 p = " << p << ": a second call at a floored DOF adds a hit, not a DOF");
   on.ResetFloorCount();
   TEST_ASSERT(on.FloorHits() == 0 && on.FlooredDOFs() == 0 && on.SigmaFloor() == kF,
               "D1 p = " << p << ": ResetFloorCount clears both counts");
}

static int Refuse(int n)
{
   std::cout << "--refuse " << n << "\n";
   switch (n)
   {
      case 1: { Ops o(1, -1.0); break; }
      case 2: { Ops o(1, kF, 1.0e5); break; }
      case 3:
      {
         Ops o(1, kF);
         real_t TS[3] = {0.0, std::numeric_limits<real_t>::quiet_NaN(), 0.0}, V[3], X[2], s;
         o.on->CloseDOF(0, TS, 0.0, V, s, X);
         break;
      }
      case 4:
      {
         Ops o(1, kF);
         Vector S;
         o.off->BackgroundTraction(0.0, S);
         real_t TS[3], V[3], X[2], s;
         MovedTraction(*o.off, S, 0, -kF, TS);
         o.off->CloseDOF(0, TS, 0.0, V, s, X);
         break;
      }
      default: MFEM_ABORT("--refuse N: N = 1..4");
   }
   std::cout << "--refuse " << n << ": accepted\n";
   return 0;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   int refuse = 0;
   for (int i = 1; i < argc; i++)
   {
      const std::string a = argv[i];
      if (a == "--refuse" && i + 1 < argc) { refuse = std::atoi(argv[++i]); }
   }
   if (refuse) { return Refuse(refuse); }

   TestKernels();
   for (int p = 1; p <= 2; p++) { TestDynamic(p); }

   std::cout << "\ntest_sigma_floor: " << (num_tests - num_failed) << " / " << num_tests
             << " passed\n";
   return num_failed ? 1 : 0;
}
