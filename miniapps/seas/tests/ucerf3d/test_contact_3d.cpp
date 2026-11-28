// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/ucerf3d/test_contact_3d.cpp — UCERF 3-D plan rev5, P1.3 (Chunhui
// Zhao, 2026-11-25: "start P1"; 2026-11-26 "proceed next"): unilateral
// contact (fault opening version B, [cycle.opening] mode = "contact") in the
// 3-D Nitsche dynamic phase, three open unknowns per fault node (dip slip,
// strike slip, gap along n).
//
//   U-C1  a block pulled apart opens with zero traction: a box cut through
//         by a planar fault and turned by a general rotation; the two faces
//         parallel to the fault are pulled at u = +-v0 t n (Dirichlet; the
//         operator needs a Dirichlet or absorbing face, and a dashpot would
//         brake the blocks), the rest are free; the two halves start at
//         +-v0 n.  Each half must move as
//         a rigid block: after N steps every owned fault DOF is open with
//         sigma_eff = 0, V = 2 v0 n, the gap d . n = 2 v0 t, and the volume
//         velocity is still +-v0 n (a traction on the fault would slow it).
//   U-C2  a closed fault under compression is bit-identical to the
//         no-contact path: the slip-weakening rupture of
//         test_nitsche_lsw_inline (TPV205 friction, 120 MPa compression,
//         nucleation patch) run with SetContactOpening(true) and without;
//         no DOF ever opens; every field equal by value (-0 == +0) until
//         slip starts, then within 1e-10 relative after N steps.  Not bit
//         for bit: the contact path measures the slip as hypot(d . dip,
//         d . strike) (d carries the gap; 363db260), the plain path
//         sqrt(d . d), so the last bits part once a DOF slips (measured
//         1e-12 after 300 steps).
//   U-C3  energy balance with opening on T-bend (ucerf3d/mesh, the 20 deg
//         bend at x = y = 0): sigma_yy = -5 MPa clamps both planes; the shear
//         sigma_xy = -10 MPa comes on smoothly over 0.5 s.  It makes the
//         west plane y = 0 slip and pulls the east plane (normal (-sin 20,
//         cos 20, 0)) apart, so the east plane opens.
//         a  the monitor's fault power int_F T_h . W equals dE_el/dd (a
//            central difference) along n, dip and strike (1e-6).  This is
//            the check that found the reversed normal part (nitsche_energy
//            .hpp, contact only).
//         b  NitscheEnergyMonitor's stepwise residual rho = E_h^{n+1} -
//            E_h^n - dt Pbar summed over 1 s (the T-15 criterion of
//            test_nitsche_lsw_inline L-06) falls by about 4 per halving of
//            dt: the discrete balance holds to O(dt^2) with DOFs open.
//
//   seas_test_contact_3d [--only C1|C2|C3]   (1 or 2 ranks; CT3_DUMP=DIR
//                                            writes U-C3's fault fields)
//   make -f tests/ucerf3d/ucerf3d.mk test-contact-3d
//
// Scope (Chunhui Zhao, 2026-11-25): tetrahedral meshes, the host step.
//
// P1.9, U-G3 (plan rev5: "U-C1 and U-C3 on the device against the host"):
// the device step's contact branch.  One device per process, so the host
// run writes its end state and the device run compares:
//
//   seas_test_contact_3d --device cpu   --write PREFIX
//   seas_test_contact_3d --device debug --compare PREFIX
//   make -f tests/ucerf3d/ucerf3d.mk test-device-contact-3d
//
//   U-G3a  U-C1's turned block pulled apart, p = 2, 200 steps.  The device
//          step has no Dirichlet faces yet (P1.11), so the two faces
//          parallel to the fault are absorbing and the halves start at
//          +-v0 n; every fault DOF stays open with sigma_eff = 0.
//   U-G3b  U-C3's T-bend run (p = 1, 1 s): the east plane opens, the west
//          plane slips.
//   U-G3c  U-G3b with the relaxed contact condition (r = +2 MPa west,
//          -1 MPa east), a 1 MPa floor and the floor taper.  With the
//          relaxed condition the floor has nothing left to do (0 hits).
//   U-G3d  U-G3b with the same r, floor and taper but the plain contact
//          condition: here the floor arm and the taper act.
//   Checks: the device step ran; the case opens; slip, slip rate, state,
//   sigma_eff, U and Udot within --tol (1e-12) of the host, relative to the
//   field's largest value, the slip rate within --tol-v (1e-10; G3Args);
//   the open flags and the floor counts equal.

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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kL = 30.0e3;
static const real_t kCuse[3] = {0.0, 0.057, 0.11};   // 0.5 C_eig (as test_nitsche_lsw_inline)
static const real_t kMuS = 0.677, kMuD = 0.525, kDc = 0.40;
static const real_t kTau0 = 70.0e6, kSigma0 = 120.0e6;
static const real_t kNucR = 9.0e3, kNucDtau = 15.0e6;

static real_t RMax(real_t v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return v; }
static int ISum(int v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); return v; }

static MaterialField Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
}

/// [-L, L]^2 x [-L, 0] in nx x ny x nz hexes of 24 tets; top z = 0 attr 1,
/// other outer faces attr 5 (the faces y = +-L attr 4 with ends4), the plane
/// y = 0 attr 3 (the box of test_nitsche_lsw_inline).  R (row-major 3 x 3,
/// or null) turns the vertices afterwards.
static Mesh *Box(int nx, int ny, int nz, const real_t *R = nullptr, bool ends4 = false)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(nx, ny, nz, 2 * kL, 2 * kL, kL));
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
      real_t zc = 0.0, yc = 0.0;
      for (int i = 0; i < verts.Size(); i++)
      {
         zc += mesh->GetVertex(verts[i])[2];
         yc += mesh->GetVertex(verts[i])[1];
      }
      zc /= verts.Size();
      yc /= verts.Size();
      mesh->SetBdrAttribute(b, (std::abs(zc) < tol) ? 1 : (ends4 && std::abs(std::abs(yc) - kL) < tol) ? 4 : 5);
   }
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool on = true;
      for (int i = 0; i < fv.Size(); i++) { on = on && std::abs(mesh->GetVertex(fv[i])[1]) < 1e-6 * kL; }
      if (on) { tris.push_back({fv[0], fv[1], fv[2]}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   if (R)
   {
      for (int v = 0; v < mesh->GetNV(); v++)
      {
         real_t *c = mesh->GetVertex(v), x[3] = {c[0], c[1], c[2]};
         for (int i = 0; i < 3; i++) { c[i] = R[3 * i] * x[0] + R[3 * i + 1] * x[1] + R[3 * i + 2] * x[2]; }
      }
   }
   return mesh;
}

static void Apply(const real_t *R, const real_t a[3], real_t b[3])
{
   for (int i = 0; i < 3; i++) { b[i] = R ? R[3 * i] * a[0] + R[3 * i + 1] * a[1] + R[3 * i + 2] * a[2] : a[i]; }
}

/// R sigma R^T (row-major 3 x 3).
static void Turn(const real_t *R, const real_t s[9], real_t o[9])
{
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         real_t v = 0.0;
         for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++) { v += R[3 * i + k] * s[3 * k + l] * R[3 * j + l]; }
         o[3 * i + j] = v;
      }
}

static void Params(FaultPointState &st)
{
   st.mu_s = kMuS;
   st.mu_d = kMuD;
   st.Dc = kDc;
   st.cohesion = 0.0;
}

static void AtRest(NitscheWaveOperator &op, real_t tau, real_t sigma)
{
   const int no = op.NumOwnedFaultDOFs();
   Vector d(3 * no), V(3 * no), psi(no);
   d = 0.0;
   V = 0.0;
   for (int j = 0; j < no; j++)
   {
      FaultPointState st = op.FaultParams(j);
      op.Closure().SeedEquilibrium(tau, 0.0, sigma, st);
      psi(j) = st.psi;
   }
   op.SetFaultState(d, V, psi);
}

// ------------------------------------------------------------------ U-G3
struct G3Args
{
   std::string device = "cpu", write, compare;
   real_t tol = 1e-12;
   /// The slip rate's limit.  With points open the slip rate of the device
   /// and the host part more than the other fields and the difference grows
   /// with the steps (measured, 1 rank: U-G3a 0 at 1 step, 1e-12 m/s at 50,
   /// 6.8e-12 at 200, 1.8e-11 at 400 of V = 2 m/s; U-G3b 4e-12 at 200, 2e-11
   /// at 400, open and closed points alike).  The open flags stay equal.
   /// The cause is not isolated; a guess is V = -X / gamma at an open point,
   /// X being the difference of two large penalty terms.
   real_t tol_v = 1e-10;
   long steps = 0;           ///< --steps N: instead of the case's own number (a measurement)
   bool host_fast = false;   ///< --host-fast 1: a second host run with the per-face fault kernels
   bool On() const { return !write.empty() || !compare.empty(); }
};

/// A reduction that does not change (the cycle's frozen relaxation r).
class ConstantReduction : public IFaultNormalStressModifier
{
public:
   explicit ConstantReduction(const Vector &r) : r_(r) {}
   int NumPoints() const override { return r_.Size(); }
   real_t NormalStressReduction(int j) const override { return r_(j); }
   void Advance(real_t, const real_t *, const real_t *) override {}
   bool Constant() const override { return true; }
private:
   Vector r_;
};

static void G3DeviceOn(NitscheWaveOperator &op, const G3Args &a)
{
   if (a.device == "cpu")
   {
      // a host run in another order of sums (the round-off reference)
      if (a.host_fast) { op.SetFastFaultKernels(true); }
      return;
   }
   op.SetMatrixFree(true);
   op.SetFastFaultKernels(true);
   op.SetDeviceFault(true);
}

/// The end state of this rank: written (host run) or compared (device run).
static void G3Finish(const std::string &name, NitscheWaveOperator &op, const Vector &U,
                     const Vector &Ud, long steps, const G3Args &a, bool all_open)
{
   const int no = op.NumOwnedFaultDOFs(), n = U.Size();
   U.HostRead();
   Ud.HostRead();
   std::vector<double> buf;
   auto put = [&buf](const Vector &v) { for (int i = 0; i < v.Size(); i++) { buf.push_back(v(i)); } };
   const Vector &d = op.Slip(), &V = op.SlipRate(), &psi = op.State(), &sig = op.LastSigmaEff();
   put(d); put(V); put(psi); put(sig); put(U); put(Ud);
   const std::vector<char> &open = op.OpenFlags();
   for (int j = 0; j < no; j++) { buf.push_back(open[j] ? 1.0 : 0.0); }
   buf.push_back(static_cast<double>(op.FloorHits()));
   buf.push_back(static_cast<double>(op.FlooredDOFs()));
   real_t gmax = 0.0, smax = 0.0;
   int nopen = 0;
   for (int j = 0; j < no; j++)
   {
      const real_t *nv = &op.FaultFrames()(9 * j);
      gmax = std::max(gmax, std::abs(d(3 * j) * nv[0] + d(3 * j + 1) * nv[1] + d(3 * j + 2) * nv[2]));
      smax = std::max(smax, std::abs(sig(j)));
      nopen += open[j] ? 1 : 0;
   }
   gmax = RMax(gmax);
   smax = RMax(smax);
   const int hits_g = ISum(static_cast<int>(op.FloorHits()));
   const int nopen_g = ISum(nopen), ntot = ISum(no);
   const std::string file = (a.write.empty() ? a.compare : a.write) + "." + name + "."
                            + std::to_string(g_rank);
   const bool dev = a.device != "cpu";
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  " << name << " on device \"" << a.device << "\": "
                << steps << " steps, open " << nopen_g << " / " << ntot << " fault DOFs, max gap "
                << gmax << " m, floor " << hits_g << " DOF-steps\n";
   }
   TEST_ASSERT(op.MatrixFreeActive() == dev && op.DeviceStepActive() == dev
               && op.DeviceFaultActive() == dev,
               name << " G1: the device step " << (dev ? "ran" : "is off") << ": matrix-free "
               << op.MatrixFreeActive() << ", device step " << op.DeviceStepActive()
               << ", device fault " << op.DeviceFaultActive());
   // --steps N (a shorter run, e.g. the one-step CUDA check V-G1): the bend
   // has not opened yet after a few steps, so the count is reported only
   const bool short_run = (a.steps > 0 && !all_open);
   TEST_ASSERT(short_run
               || (nopen_g > 0 && gmax > 0.0 && (!all_open || (nopen_g == ntot && smax == 0.0))),
               name << " G2: the fault opens (" << nopen_g << " of " << ntot << " DOFs, gap "
               << gmax << " m" << (all_open ? "; every DOF open with sigma_eff = 0" : "")
               << (short_run ? "; a short run (--steps), not required" : "") << ")");
   if (!a.write.empty())
   {
      std::ofstream o(file, std::ios::binary);
      const long m = static_cast<long>(buf.size());
      o.write(reinterpret_cast<const char *>(&m), sizeof(m));
      o.write(reinterpret_cast<const char *>(buf.data()), sizeof(double) * buf.size());
      return;
   }
   std::ifstream in(file, std::ios::binary);
   long m = -1;
   in.read(reinterpret_cast<char *>(&m), sizeof(m));
   std::vector<double> ref(m > 0 ? m : 0);
   if (m > 0) { in.read(reinterpret_cast<char *>(ref.data()), sizeof(double) * ref.size()); }
   int ok = (in && m == static_cast<long>(buf.size())) ? 1 : 0;
   MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   TEST_ASSERT(ok == 1, name << " G0: the host file " << file << " has this run's layout");
   if (ok != 1) { return; }
   // d, V, psi, sigma_eff, U, Udot: the largest difference over the largest value
   const int len[6] = {3 * no, 3 * no, no, no, n, n};
   real_t rel[6];
   int off = 0;
   for (int f = 0; f < 6; f++)
   {
      real_t dmax = 0.0, vmax = 0.0;
      for (int i = 0; i < len[f]; i++)
      {
         dmax = std::max(dmax, std::abs(buf[off + i] - ref[off + i]));
         vmax = std::max(vmax, std::abs(ref[off + i]));
      }
      off += len[f];
      dmax = RMax(dmax);
      vmax = RMax(vmax);
      rel[f] = (vmax > 0.0) ? dmax / vmax : dmax;
   }
   // the slip rate's difference over the DOFs open and closed at the end
   real_t dV[2] = {0.0, 0.0};
   {
      const int oo = 3 * no + 3 * no + no + no + n + n;
      for (int j = 0; j < no; j++)
         for (int c = 0; c < 3; c++)
         {
            const int i = 3 * no + 3 * j + c, w = (ref[oo + j] != 0.0) ? 0 : 1;
            dV[w] = std::max(dV[w], std::abs(buf[i] - ref[i]));
         }
      dV[0] = RMax(dV[0]);
      dV[1] = RMax(dV[1]);
      if (g_rank == 0)
      {
         std::cout << "  " << name << ": largest slip-rate difference " << dV[0]
                   << " m/s over the open DOFs, " << dV[1] << " m/s over the closed DOFs\n";
      }
   }
   int flips = 0, host_open = 0;
   for (int j = 0; j < no; j++)
   {
      flips += (buf[off + j] != ref[off + j]) ? 1 : 0;
      host_open += (ref[off + j] != 0.0) ? 1 : 0;
   }
   off += no;
   flips = ISum(flips);
   host_open = ISum(host_open);
   const long hits = static_cast<long>(buf[off]), hits_h = static_cast<long>(ref[off]);
   const int fl = static_cast<int>(buf[off + 1]), fl_h = static_cast<int>(ref[off + 1]);
   int cnt = (hits == hits_h && fl == fl_h) ? 1 : 0;
   MPI_Allreduce(MPI_IN_PLACE, &cnt, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   TEST_ASSERT(rel[0] <= a.tol && rel[1] <= a.tol_v && rel[2] <= a.tol && rel[3] <= a.tol,
               name << " G3: device against host on the fault: slip " << rel[0] << ", state "
               << rel[2] << ", sigma_eff " << rel[3] << " (relative, limit " << a.tol
               << "); slip rate " << rel[1] << " (limit " << a.tol_v << ")");
   TEST_ASSERT(rel[4] <= a.tol && rel[5] <= a.tol,
               name << " G4: the volume fields: U " << rel[4] << ", Udot " << rel[5]
               << " (relative, limit " << a.tol << ")");
   TEST_ASSERT(flips == 0 && cnt == 1,
               name << " G5: the open flags are the host's (" << flips << " differ; host open "
               << host_open << ", device " << nopen_g << "), and so are the floor counts on "
               "every rank (rank 0: " << hits << " DOF-steps on " << fl << " DOFs, host "
               << hits_h << " on " << fl_h << ")");
}

static void TestG3Block(const G3Args &a);
static void TestG3Bend(const G3Args &a, int variant);

// ------------------------------------------------------------------ U-C1
static void TestC1(int p)
{
   // a general rotation: 30 deg about z, then 40 deg about x
   const real_t a = M_PI / 6, b = 2 * M_PI / 9;
   const real_t Rz[9] = {std::cos(a), -std::sin(a), 0, std::sin(a), std::cos(a), 0, 0, 0, 1};
   const real_t Rx[9] = {1, 0, 0, 0, std::cos(b), -std::sin(b), 0, std::sin(b), std::cos(b)};
   real_t R[9];
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         R[3 * i + j] = 0.0;
         for (int k = 0; k < 3; k++) { R[3 * i + j] += Rx[3 * i + k] * Rz[3 * k + j]; }
      }
   const real_t e_y[3] = {0, 1, 0}, e_z[3] = {0, 0, 1};
   real_t ny[3], up[3];
   Apply(R, e_y, ny);
   Apply(R, e_z, up);
   std::unique_ptr<Mesh> m(Box(4, 4, 2, R, true));
   ParMesh pm(MPI_COMM_WORLD, *m);
   m.reset();
   // the two faces parallel to the fault pull at +-v0 n_y (u = +-v0 t n_y),
   // the rest free: the exact motion is two rigid blocks
   const real_t v0 = 1.0;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1, 5};
   bc.dirichlet_attrs = {4};
   const std::array<real_t, 3> nyv = {ny[0], ny[1], ny[2]};
   bc.default_dirichlet_func = [nyv, v0](const Vector &x, real_t t, Vector &u) {
      u.SetSize(3);
      const real_t s = (x(0) * nyv[0] + x(1) * nyv[1] + x(2) * nyv[2] > 0.0) ? 1.0 : -1.0;
      for (int k = 0; k < 3; k++) { u(k) = s * v0 * t * nyv[k]; }
   };
   SlipWeakeningClosure cl;
   NitscheFaultSetup fs;
   fs.closure = &cl;
   for (int d = 0; d < 3; d++)
   {
      fs.orientation.ref_normal[d] = -ny[d];
      fs.orientation.up[d] = up[d];
      fs.up[d] = up[d];
   }
   fs.orientation.fault_attrs = {3};
   fs.prestress = [](const real_t[3], real_t s[9]) { for (int i = 0; i < 9; i++) { s[i] = 0.0; } };
   fs.parameters = [](const real_t[3], FaultPointState &st) { Params(st); };
   NitscheWaveOperator op(pm, p, Material(), bc, {}, &fs);
   op.SetContactOpening(true);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   op.SetTimeStep(dt);
   AtRest(op, 0.0, 1.0e6);
   // the halves at +-v0 n_y (n_y the turned +y, from the element centroid)
   const ParFiniteElementSpace &fes = op.Elasticity().GetFESpace();
   const int n = op.Height();
   Vector U(n), Ud(n), side(n);
   U = 0.0;
   Ud = 0.0;
   for (int e = 0; e < pm.GetNE(); e++)
   {
      Vector c(3);
      pm.GetElementCenter(e, c);
      const real_t s = (c(0) * ny[0] + c(1) * ny[1] + c(2) * ny[2] > 0.0) ? 1.0 : -1.0;
      Array<int> vd;
      fes.GetElementVDofs(e, vd);
      const int nd = vd.Size() / 3;
      for (int i = 0; i < nd; i++)
         for (int k = 0; k < 3; k++)
         {
            const int q = fes.GetOrdering() == Ordering::byNODES ? k * nd + i : 3 * i + k;
            Ud(vd[q]) = s * v0 * ny[k];
            side(vd[q]) = s;
         }
   }
   // the consistent fault rate V^0 = 2 v0 n (the jump rate of the start).
   // V^0 = 0 instead gives the step's 2-step zigzag in V (V . n = 4 v0, 0,
   // 4 v0, ...; test_nitsche_lsw_inline L-07) with the gap still exact.
   {
      const int no = op.NumOwnedFaultDOFs();
      Vector d0(3 * no), V0(3 * no);
      d0 = 0.0;
      for (int j = 0; j < no; j++)
      {
         const real_t *nv = &op.FaultFrames()(9 * j);
         const real_t sg = (nv[0] * ny[0] + nv[1] * ny[1] + nv[2] * ny[2]) > 0 ? 1.0 : -1.0;
         for (int c = 0; c < 3; c++) { V0(3 * j + c) = sg * 2 * v0 * ny[c]; }
      }
      Vector psi(op.State());
      op.SetFaultState(d0, V0, psi);
   }
   op.SetTime(0.0);
   CentralDifferenceSolver ode;
   ode.Init(op);
   real_t t = 0.0;
   const int nsteps = 200;
   int ever_closed = 0;
   for (int s = 0; s < nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
      if (s > 0) { ever_closed += op.NumOwnedFaultDOFs() - static_cast<int>(std::count(op.OpenFlags().begin(), op.OpenFlags().end(), char(1))); }
   }
   ever_closed = ISum(ever_closed);
   // volume: still +-v0 n_y
   real_t ev = 0.0;
   for (int e = 0; e < pm.GetNE(); e++)
   {
      Array<int> vd;
      fes.GetElementVDofs(e, vd);
      const int nd = vd.Size() / 3;
      for (int i = 0; i < nd; i++)
         for (int k = 0; k < 3; k++)
         {
            const int q = fes.GetOrdering() == Ordering::byNODES ? k * nd + i : 3 * i + k;
            ev = std::max(ev, std::abs(Ud(vd[q]) - side(vd[q]) * v0 * ny[k]));
         }
   }
   // fault: open, sigma_eff 0, V = 2 v0 n, gap 2 v0 t (n the oriented frame normal)
   const int no = op.NumOwnedFaultDOFs();
   const Vector &d = op.Slip(), &V = op.SlipRate(), &sig = op.LastSigmaEff();
   real_t eV = 0.0, eg = 0.0, es = 0.0, et = 0.0;
   int nopen = 0;
   for (int j = 0; j < no; j++)
   {
      const real_t *nv = &op.FaultFrames()(9 * j);
      const real_t sgn = (nv[0] * ny[0] + nv[1] * ny[1] + nv[2] * ny[2]) > 0 ? 1.0 : -1.0;
      real_t g = 0.0, vt = 0.0;
      for (int c = 0; c < 3; c++)
      {
         g += d(3 * j + c) * nv[c];
         eV = std::max(eV, std::abs(V(3 * j + c) - sgn * 2 * v0 * ny[c]));
      }
      for (int c = 0; c < 3; c++) { const real_t r = d(3 * j + c) - g * nv[c]; vt += r * r; }
      eg = std::max(eg, std::abs(std::abs(g) - 2 * v0 * t));
      et = std::max(et, std::sqrt(vt));
      es = std::max(es, std::abs(sig(j)));
      nopen += op.OpenFlags()[j] ? 1 : 0;
   }
   ev = RMax(ev); eV = RMax(eV); eg = RMax(eg); es = RMax(es); et = RMax(et);
   nopen = ISum(nopen);
   const int ntot = ISum(no);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(4) << "  p = " << p << ", " << nsteps << " steps of " << dt
                << " s: open " << nopen << " / " << ntot << " fault DOFs (closed at any step after "
                << "the first: " << ever_closed << "); max |Udot - (+-v0 n)| = " << ev
                << " m/s, max |V - 2 v0 n| = " << eV << " m/s, max ||g| - 2 v0 t| = " << eg
                << " m (2 v0 t = " << 2 * v0 * t << " m), tangential slip " << et
                << " m, max |sigma_eff| = " << es << " Pa\n";
   }
   TEST_ASSERT(nopen == ntot && ever_closed == 0 && es == 0.0,
               "U-C1 p = " << p << ": every fault DOF open at every step, sigma_eff = 0");
   TEST_ASSERT(ev <= 1e-9 * v0 && eV <= 1e-9 * v0,
               "U-C1 p = " << p << ": both halves still move rigidly at +-v0 n (" << ev
               << " m/s) and V = 2 v0 n (" << eV << " m/s): zero traction on the open fault");
   TEST_ASSERT(eg <= 1e-9 * 2 * v0 * t && et <= 1e-9 * 2 * v0 * t,
               "U-C1 p = " << p << ": the gap is 2 v0 t (" << eg << " m) with no tangential slip ("
               << et << " m)");
}

// ------------------------------------------------------------------ U-C2
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
      const real_t v = kNucDtau * NucF(std::hypot(x[0] - TPV102Params::hypo_along_strike,
                                                  x[2] + TPV102Params::hypo_down_dip));
      s[1] = s[3] = v;
   };
   fs.nucleation_temporal = [](real_t t) { return NucleationTemporal(t); };
   fs.parameters = [](const real_t[3], FaultPointState &st) { Params(st); };
   return fs;
}

static int BitDiff(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return 1 << 30; }
   int n = 0;
   for (int i = 0; i < a.Size(); i++) { n += !(a(i) == b(i)); }   // by value: -0 == +0
   return n;
}

static void TestC2(int p)
{
   std::unique_ptr<Mesh> m(Box(6, 6, 3));
   ParMesh pm(MPI_COMM_WORLD, *m);
   m.reset();
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.absorbing_attrs = {5};
   SlipWeakeningClosure cl;
   NitscheFaultSetup fa = LSWSetup(cl), fb = LSWSetup(cl);
   NitscheWaveOperator ct(pm, p, Material(), bc, {}, &fa);
   NitscheWaveOperator pl(pm, p, Material(), bc, {}, &fb);
   ct.SetContactOpening(true);
   const real_t dt = ComputeMaxDt(ct.Sample(), pm, p, kCuse[p]);
   const int n = ct.Height(), nsteps = 300;
   Vector Ua(n), Va(n), Ub(n), Vb(n);
   Ua = 0.0; Va = 0.0; Ub = 0.0; Vb = 0.0;
   real_t ta = 0.0, tb = 0.0;
   int ever_open = 0, first = -1;
   real_t vmax = 0.0, smin = 1e300, vfirst = 0.0, v_at_first = 0.0;
   for (NitscheWaveOperator *op : {&ct, &pl})
   {
      op->SetTimeStep(dt);
      AtRest(*op, kTau0, kSigma0);
      op->SetTime(0.0);
   }
   CentralDifferenceSolver oa, ob;
   oa.Init(ct);
   ob.Init(pl);
   for (int s = 0; s < nsteps; s++)
   {
      real_t h = dt;
      oa.Step(Ua, Va, ta, h);
      h = dt;
      ob.Step(Ub, Vb, tb, h);
      // the first step where any field differs by value (-0 == +0), and
      // the largest slip rate at or before it
      vfirst = std::max(vfirst, [&]() { real_t m = 0.0; for (int i = 0; i < pl.SlipRate().Size(); i++) { m = std::max(m, std::abs(pl.SlipRate()(i))); } return RMax(m); }());
      if (first < 0)
      {
         int dd = BitDiff(Ua, Ub) + BitDiff(Va, Vb) + BitDiff(ct.Slip(), pl.Slip())
                  + BitDiff(ct.SlipRate(), pl.SlipRate()) + BitDiff(ct.State(), pl.State());
         dd = ISum(dd);
         if (dd > 0) { first = s + 1; v_at_first = vfirst; }
      }
      ever_open += static_cast<int>(std::count(ct.OpenFlags().begin(), ct.OpenFlags().end(), char(1)));
      for (int j = 0; j < ct.LastSigmaEff().Size(); j++) { smin = std::min(smin, ct.LastSigmaEff()(j)); }
   }
   for (int i = 0; i < ct.SlipRate().Size(); i++) { vmax = std::max(vmax, std::abs(ct.SlipRate()(i))); }
   // relative differences: max |a - b| / max |b| per field
   auto rel = [](const Vector &a, const Vector &b) {
      real_t num = 0.0, den = 0.0;
      for (int i = 0; i < a.Size(); i++) { num = std::max(num, std::abs(a(i) - b(i))); den = std::max(den, std::abs(b(i))); }
      num = RMax(num); den = RMax(den);
      return den > 0.0 ? num / den : num;
   };
   const real_t rU = rel(Ua, Ub), rUd = rel(Va, Vb), rd = rel(ct.Slip(), pl.Slip()),
                rV = rel(ct.SlipRate(), pl.SlipRate()), rpsi = rel(ct.State(), pl.State()),
                rs = rel(ct.LastSigmaEff(), pl.LastSigmaEff());
   const real_t rmax = std::max({rU, rUd, rd, rV, rpsi, rs});
   if (g_rank == 0)
   {
      std::cout << std::setprecision(3) << "  p = " << p << " relative differences: U " << rU << ", Udot " << rUd
                << ", d " << rd << ", V " << rV << ", psi " << rpsi << ", sigma_eff " << rs << "\n";
   }
   int diff[6] = {BitDiff(Ua, Ub), BitDiff(Va, Vb), BitDiff(ct.Slip(), pl.Slip()),
                  BitDiff(ct.SlipRate(), pl.SlipRate()), BitDiff(ct.State(), pl.State()),
                  BitDiff(ct.LastSigmaEff(), pl.LastSigmaEff())};
   MPI_Allreduce(MPI_IN_PLACE, diff, 6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   ever_open = ISum(ever_open);
   vmax = RMax(vmax);
   smin = -RMax(-smin);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  p = " << p << ", " << nsteps << " steps to t = " << ta
                << " s: max slip rate " << vmax << " m/s, min sigma_eff " << smin / 1e6
                << " MPa, DOF-steps open " << ever_open << "; differing entries U " << diff[0]
                << ", Udot " << diff[1] << ", d " << diff[2] << ", V " << diff[3] << ", psi "
                << diff[4] << ", sigma_eff " << diff[5] << "; first difference at step " << first << "\n";
   }
   TEST_ASSERT(vmax > 0.1 && ever_open == 0,
               "U-C2 p = " << p << ": the rupture runs (max V " << vmax << " m/s) and no DOF opens");
   TEST_ASSERT(first < 0 || v_at_first > 0.0,
               "U-C2 p = " << p << ": contact on = contact off by value until slip starts (first "
               "difference at step " << first << ", max slip rate then " << v_at_first << " m/s)");
   TEST_ASSERT(rmax <= 1e-10,
               "U-C2 p = " << p << ": after " << nsteps << " steps every field within 1e-10 "
               "relative (max " << rmax << "; the contact path measures the slip as hypot(d.dip, "
               "d.strike), the plain path sqrt(d.d), so the last bits differ once slip starts)");
}

// ------------------------------------------------------------------ U-C3
struct C3Run
{
   real_t dt = 0.0, sum_rho = 0.0, sum_abs = 0.0, Emax = 0.0, W_pre = 0.0, W_fric = 0.0;
   real_t gap_max = 0.0, rho_first = 0.0;
   real_t pair_dE[3] = {0, 0, 0}, pair_P[3] = {0, 0, 0};   ///< n, dip, strike (the dt_scale 1 run)
   int open_max = 0, slip_dofs = 0;
   long steps = 0;
};

static const real_t kS = 5.0e6, kTxy = -10.0e6, kRamp = 0.5;
static const bool kC3Contact = true;

static C3Run RunC3(ParMesh &pm, int p, real_t tfinal, real_t dt_scale)
{
   C3Run R;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.absorbing_attrs = {2, 5};
   SlipWeakeningClosure cl;
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
   // the clamp sigma_yy = -kS from the start; the shear sigma_xy = kTxy
   // comes on smoothly over kRamp (the SCEC smooth step), so slip and
   // opening begin from rest
   fs.prestress = [](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[4] = -kS;
      if (kRamp <= 0.0) { s[1] = s[3] = kTxy; }
   };
   if (kRamp > 0.0)
   {
      fs.nucleation_spatial = [](const real_t[3], real_t s[9]) {
         for (int i = 0; i < 9; i++) { s[i] = 0.0; }
         s[1] = s[3] = kTxy;
      };
      fs.nucleation_temporal = [](real_t t) {
         if (t <= 0.0) { return real_t(0); }
         if (t >= kRamp) { return real_t(1); }
         return std::exp((t - kRamp) * (t - kRamp) / (t * (t - 2.0 * kRamp)));
      };
   }
   fs.parameters = [](const real_t[3], FaultPointState &st) { Params(st); };
   NitscheWaveOperator op(pm, p, Material(), bc, {}, &fs);
   if (kC3Contact) { op.SetContactOpening(true); }
   const real_t dt = dt_scale * ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   R.dt = dt;
   op.SetTimeStep(dt);
   AtRest(op, 0.0, kS);
   NitscheEnergyMonitor mon(op, false);
   const int n = op.Height();
   Vector U(n), Ud(n), A0(n);
   U = 0.0;
   Ud = 0.0;
   op.SetTime(0.0);
   op.Mult(U, Ud, A0);
   NitscheEnergyRecord r = mon.Evaluate(U, Ud, A0, 0.0), prev = r;
   mon.Start(r);
   CentralDifferenceSolver ode;
   ode.Init(op);
   real_t t = 0.0;
   R.steps = std::lround(std::ceil(tfinal / dt - 1e-9));
   const int no = op.NumOwnedFaultDOFs();
   for (long s = 0; s < R.steps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
      r = mon.Evaluate(U, Ud, op.LastAcceleration(), t);
      const real_t rho = mon.Step(r);
      if (s == 0) { R.rho_first = rho; }
      R.W_pre += 0.5 * dt * (r.P_pre + prev.P_pre + r.P_nuc + prev.P_nuc);
      R.W_fric += 0.5 * dt * (r.P_fric + prev.P_fric);
      R.Emax = std::max(R.Emax, std::abs(r.E_h));
      if (kC3Contact)
      {
         R.open_max = std::max(R.open_max, ISum(static_cast<int>(
            std::count(op.OpenFlags().begin(), op.OpenFlags().end(), char(1)))));
      }
      prev = r;
   }
   if (dt_scale == 1.0)
   {
      // dE_el/dd along a unit field W against the monitor's int_F T_h . W
      const Vector d0(op.Slip()), V0(op.SlipRate()), psi0(op.State());
      for (int dir = 0; dir < 3; dir++)   // 0 = n, 1 = dip, 2 = strike
      {
         Vector W(3 * no);
         for (int j = 0; j < no; j++)
            for (int c = 0; c < 3; c++) { W(3 * j + c) = op.FaultFrames()(9 * j + 3 * dir + c); }
         const real_t eps = 1e-3;
         Vector dp(d0), dm(d0);
         dp.Add(eps, W);
         dm.Add(-eps, W);
         op.SetFaultState(dp, W, psi0);
         const real_t Ep = mon.Evaluate(U, Ud, op.LastAcceleration(), t).E_el;
         op.SetFaultState(dm, W, psi0);
         const real_t Em = mon.Evaluate(U, Ud, op.LastAcceleration(), t).E_el;
         op.SetFaultState(d0, W, psi0);
         const NitscheEnergyRecord q = mon.Evaluate(U, Ud, op.LastAcceleration(), t);
         R.pair_dE[dir] = (Ep - Em) / (2 * eps);
         R.pair_P[dir] = q.P_fric + q.P_pre + q.P_nuc;
      }
      op.SetFaultState(d0, V0, psi0);
   }
   real_t gmax = 0.0;
   int nslip = 0;
   for (int j = 0; j < no; j++)
   {
      const real_t *nv = &op.FaultFrames()(9 * j);
      real_t g = 0.0, tt = 0.0;
      for (int c = 0; c < 3; c++) { g += op.Slip()(3 * j + c) * nv[c]; }
      for (int c = 0; c < 3; c++) { const real_t q = op.Slip()(3 * j + c) - g * nv[c]; tt += q * q; }
      gmax = std::max(gmax, std::abs(g));
      nslip += (!(kC3Contact && op.OpenFlags()[j]) && std::sqrt(tt) > 1e-3) ? 1 : 0;
   }
   if (const char *dir = std::getenv("CT3_DUMP"))   // x y z gap |tangential slip| open, for the review figure
   {
      if (dt_scale == 1.0)
      {
         std::ofstream o(std::string(dir) + "/c3_fault." + std::to_string(g_rank) + ".txt");
         o.precision(10);
         for (int j = 0; j < no; j++)
         {
            const real_t *nv = &op.FaultFrames()(9 * j), *x = &op.FaultPositions()(3 * j);
            real_t g = 0.0, tt = 0.0;
            for (int c = 0; c < 3; c++) { g += op.Slip()(3 * j + c) * nv[c]; }
            for (int c = 0; c < 3; c++) { const real_t q = op.Slip()(3 * j + c) - g * nv[c]; tt += q * q; }
            o << x[0] << " " << x[1] << " " << x[2] << " " << g << " " << std::sqrt(tt) << " "
              << int(op.OpenFlags()[j]) << "\n";
         }
      }
   }
   R.gap_max = RMax(gmax);
   R.slip_dofs = ISum(nslip);
   R.sum_rho = mon.SumRho();
   R.sum_abs = mon.SumAbsRho();
   return R;
}

static void TestC3(int p, real_t tfinal)
{
   Mesh sm("ucerf3d/mesh/T-bend.msh", 1, 1);
   ParMesh pm(MPI_COMM_WORLD, sm);
   sm.Clear();
   C3Run a = RunC3(pm, p, tfinal, 1.0), b = RunC3(pm, p, tfinal, 0.5), c = RunC3(pm, p, tfinal, 0.25);
   const real_t r1 = std::abs(a.sum_rho / b.sum_rho), r2 = std::abs(b.sum_rho / c.sum_rho);
   const real_t scale = std::max({std::abs(a.W_pre), std::abs(a.W_fric), a.Emax});
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6);
      for (const C3Run *e : {&a, &b, &c})
      {
         std::cout << "  dt = " << e->dt << " (" << e->steps << " steps): open DOFs up to "
                   << e->open_max << ", max gap " << e->gap_max << " m, closed DOFs slipping > 1 mm "
                   << e->slip_dofs << "; work: prestress " << e->W_pre << " J, friction "
                   << e->W_fric << " J; max |E_h| " << e->Emax << " J; sum rho = " << e->sum_rho
                   << " (sum |rho| " << e->sum_abs << "), of it step 1: " << e->rho_first << "\n";
      }
      std::cout << "  residual ratios: dt0/dt1 = " << r1 << ", dt1/dt2 = " << r2
                << "; |sum rho| / work scale = " << std::abs(a.sum_rho) / scale << "\n";
   }
   real_t pair_err = 0.0, pair_scale = 0.0;
   for (int k = 0; k < 3; k++)
   {
      pair_err = std::max(pair_err, std::abs(a.pair_dE[k] - a.pair_P[k]));
      pair_scale = std::max(pair_scale, std::abs(a.pair_dE[k]));
   }
   if (g_rank == 0)
   {
      std::cout << "  pairing at t = " << tfinal << " s, dE_el/dd against the monitor's int_F T_h . W: n "
                << a.pair_dE[0] << " / " << a.pair_P[0] << ", dip " << a.pair_dE[1] << " / "
                << a.pair_P[1] << ", strike " << a.pair_dE[2] << " / " << a.pair_P[2] << "\n";
   }
   TEST_ASSERT(pair_err <= 1e-6 * pair_scale && std::abs(a.pair_dE[0]) > 0.01 * pair_scale,
               "U-C3 p = " << p << ": the monitor's fault power is dE_el/dd along n, dip and "
               "strike (max difference " << pair_err / pair_scale << " of the largest)");
   TEST_ASSERT(a.open_max > 0 && c.open_max > 0 && a.gap_max > 0.0 && a.slip_dofs > 0,
               "U-C3 p = " << p << ": the east plane opens (" << a.open_max << " DOFs, gap "
               << a.gap_max << " m) and the west plane slips (" << a.slip_dofs << " DOFs)");
   TEST_ASSERT(r1 >= 3.5 && r1 <= 4.5 && r2 >= 3.5 && r2 <= 4.5,
               "U-C3 p = " << p << ": the energy residual falls by [3.5, 4.5] per dt halving ("
               << r1 << ", " << r2 << ")");
}

// ------------------------------------------------------------------ U-G3
static void TestG3Block(const G3Args &a)
{
   const int p = 2;
   // U-C1's box and rotation
   const real_t al = M_PI / 6, be = 2 * M_PI / 9;
   const real_t Rz[9] = {std::cos(al), -std::sin(al), 0, std::sin(al), std::cos(al), 0, 0, 0, 1};
   const real_t Rx[9] = {1, 0, 0, 0, std::cos(be), -std::sin(be), 0, std::sin(be), std::cos(be)};
   real_t R[9];
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         R[3 * i + j] = 0.0;
         for (int k = 0; k < 3; k++) { R[3 * i + j] += Rx[3 * i + k] * Rz[3 * k + j]; }
      }
   const real_t e_y[3] = {0, 1, 0}, e_z[3] = {0, 0, 1};
   real_t ny[3], up[3];
   Apply(R, e_y, ny);
   Apply(R, e_z, up);
   std::unique_ptr<Mesh> m(Box(4, 4, 2, R, true));
   ParMesh pm(MPI_COMM_WORLD, *m);
   m.reset();
   const real_t v0 = 1.0;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1, 5};
   bc.absorbing_attrs = {4};   // U-C1 pulls here (Dirichlet: not on the device before P1.11)
   SlipWeakeningClosure cl;
   NitscheFaultSetup fs;
   fs.closure = &cl;
   for (int d = 0; d < 3; d++)
   {
      fs.orientation.ref_normal[d] = -ny[d];
      fs.orientation.up[d] = up[d];
      fs.up[d] = up[d];
   }
   fs.orientation.fault_attrs = {3};
   fs.prestress = [](const real_t[3], real_t s[9]) { for (int i = 0; i < 9; i++) { s[i] = 0.0; } };
   fs.parameters = [](const real_t[3], FaultPointState &st) { Params(st); };
   NitscheWaveOperator op(pm, p, Material(), bc, {}, &fs);
   op.SetContactOpening(true);
   G3DeviceOn(op, a);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   op.SetTimeStep(dt);
   AtRest(op, 0.0, 1.0e6);
   const ParFiniteElementSpace &fes = op.Elasticity().GetFESpace();
   const int n = op.Height();
   Vector U(n), Ud(n);
   U.UseDevice(true);
   Ud.UseDevice(true);
   U = 0.0;
   Ud = 0.0;
   Ud.HostReadWrite();
   for (int e = 0; e < pm.GetNE(); e++)
   {
      Vector c(3);
      pm.GetElementCenter(e, c);
      const real_t s = (c(0) * ny[0] + c(1) * ny[1] + c(2) * ny[2] > 0.0) ? 1.0 : -1.0;
      Array<int> vd;
      fes.GetElementVDofs(e, vd);
      const int nd = vd.Size() / 3;
      for (int i = 0; i < nd; i++)
         for (int k = 0; k < 3; k++)
         {
            const int q = fes.GetOrdering() == Ordering::byNODES ? k * nd + i : 3 * i + k;
            Ud(vd[q]) = s * v0 * ny[k];
         }
   }
   {
      const int no = op.NumOwnedFaultDOFs();
      Vector d0(3 * no), V0(3 * no);
      d0 = 0.0;
      for (int j = 0; j < no; j++)
      {
         const real_t *nv = &op.FaultFrames()(9 * j);
         const real_t sg = (nv[0] * ny[0] + nv[1] * ny[1] + nv[2] * ny[2]) > 0 ? 1.0 : -1.0;
         for (int c = 0; c < 3; c++) { V0(3 * j + c) = sg * 2 * v0 * ny[c]; }
      }
      Vector psi(op.State());
      op.SetFaultState(d0, V0, psi);
   }
   op.SetTime(0.0);
   CentralDifferenceSolver ode;
   ode.Init(op);
   real_t t = 0.0;
   const int nsteps = a.steps > 0 ? static_cast<int>(a.steps) : 200;
   for (int s = 0; s < nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
   }
   if (g_rank == 0) { std::cout << "  U-G3a: dt = " << dt << " s, t = " << t << " s\n"; }
   G3Finish("U-G3a", op, U, Ud, nsteps, a, true);
}

/// variant 0: U-G3b; 1: U-G3c (relaxed contact); 2: U-G3d (the floor and
/// its taper with a reduction the contact condition does not see).
static void TestG3Bend(const G3Args &a, int variant)
{
   const bool relaxed = variant == 1, with_r = variant >= 1;
   const int p = 1;
   Mesh sm("ucerf3d/mesh/T-bend.msh", 1, 1);
   ParMesh pm(MPI_COMM_WORLD, sm);
   sm.Clear();
   // U-C3's setup (RunC3)
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.absorbing_attrs = {2, 5};
   SlipWeakeningClosure cl;
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
      s[4] = -kS;
   };
   fs.nucleation_spatial = [](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = kTxy;
   };
   fs.nucleation_temporal = [](real_t t) {
      if (t <= 0.0) { return real_t(0); }
      if (t >= kRamp) { return real_t(1); }
      return std::exp((t - kRamp) * (t - kRamp) / (t * (t - 2.0 * kRamp)));
   };
   fs.parameters = [](const real_t[3], FaultPointState &st) { Params(st); };
   if (with_r) { fs.sigma_floor = 1.0e6; }
   NitscheWaveOperator op(pm, p, Material(), bc, {}, &fs);
   op.SetContactOpening(true);
   std::unique_ptr<ConstantReduction> red;
   if (with_r)
   {
      // r = +2 MPa on the west plane (x < 0), -1 MPa on the east plane
      const int no = op.NumOwnedFaultDOFs();
      Vector r(no);
      for (int j = 0; j < no; j++) { r(j) = (op.FaultPositions()(3 * j) < 0.0) ? 2.0e6 : -1.0e6; }
      red = std::make_unique<ConstantReduction>(r);
      op.SetNormalStressModifier(red.get());
      if (relaxed) { op.SetRelaxedContact(true); }
      op.SetFloorTaper(true);
   }
   G3DeviceOn(op, a);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   op.SetTimeStep(dt);
   AtRest(op, 0.0, kS);
   const int n = op.Height();
   Vector U(n), Ud(n);
   U.UseDevice(true);
   Ud.UseDevice(true);
   U = 0.0;
   Ud = 0.0;
   op.SetTime(0.0);
   CentralDifferenceSolver ode;
   ode.Init(op);
   real_t t = 0.0;
   const long steps = a.steps > 0 ? a.steps : std::lround(std::ceil(1.0 / dt - 1e-9));
   for (long s = 0; s < steps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
   }
   G3Finish(variant == 0 ? "U-G3b" : relaxed ? "U-G3c" : "U-G3d", op, U, Ud, steps, a, false);
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   std::string only = (argc == 3 && std::strcmp(argv[1], "--only") == 0) ? argv[2] : "";
   G3Args g3;
   for (int i = 1; i + 1 < argc; i += 2)
   {
      const std::string k = argv[i], v = argv[i + 1];
      if (k == "--device") { g3.device = v; }
      else if (k == "--write") { g3.write = v; }
      else if (k == "--compare") { g3.compare = v; }
      else if (k == "--tol") { g3.tol = std::atof(v.c_str()); }
      else if (k == "--tol-v") { g3.tol_v = std::atof(v.c_str()); }
      else if (k == "--case") { only = v; }
      else if (k == "--steps") { g3.steps = std::atol(v.c_str()); }
      else if (k == "--host-fast") { g3.host_fast = (v == "1"); }
   }
   if (g3.On())
   {
      // U-G3 (P1.9): one device per process
      MFEM_VERIFY(g3.write.empty() != g3.compare.empty()
                  && ((g3.device == "cpu") == g3.compare.empty() || g3.host_fast),
                  "test_contact_3d: --device cpu --write PREFIX, or --device debug|cuda "
                  "--compare PREFIX");
      Device device(g3.device);
#if defined(HYPRE_USING_GPU) && (MFEM_HYPRE_VERSION >= 23100)
      HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);   // as seas_nitsche_dyn_driver
      HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
#endif
      if (g_rank == 0)
      {
         std::cout << "test_contact_3d (UCERF 3-D plan rev5, P1.9, U-G3) on " << Mpi::WorldSize()
                   << " rank(s), ";
         device.Print();
      }
      if (only.empty() || only == "G3a") { TestG3Block(g3); }
      if (only.empty() || only == "G3b") { TestG3Bend(g3, 0); }
      if (only.empty() || only == "G3c") { TestG3Bend(g3, 1); }
      if (only.empty() || only == "G3d") { TestG3Bend(g3, 2); }
      if (g_rank == 0)
      {
         std::cout << "\n" << (num_tests - num_failed) << "/" << num_tests << " tests passed\n";
      }
      return num_failed == 0 ? 0 : 1;
   }
   if (g_rank == 0)
   {
      std::cout << "test_contact_3d (UCERF 3-D plan rev5, P1.3) on " << Mpi::WorldSize()
                << " rank(s)\n";
   }
   if (only.empty() || only == "C1")
   {
      if (g_rank == 0) { std::cout << "\n=== U-C1 a block pulled apart (turned box, free faces) ===\n"; }
      for (int p : {1, 2}) { TestC1(p); }
   }
   if (only.empty() || only == "C2")
   {
      if (g_rank == 0) { std::cout << "\n=== U-C2 closed under compression: contact on = off ===\n"; }
      for (int p : {1, 2}) { TestC2(p); }
   }
   if (only.empty() || only == "C3")
   {
      if (g_rank == 0) { std::cout << "\n=== U-C3 energy balance with opening, T-bend ===\n"; }
      TestC3(1, 1.0);
   }
   if (g_rank == 0)
   {
      std::cout << "\n" << (num_tests - num_failed) << "/" << num_tests << " tests passed\n";
   }
   return num_failed == 0 ? 0 : 1;
}
