// Nitsche dynamic rupture: the standalone options of the complex-geometries
// plan (rev 2, Section 8, gates G-unit and G-inactive; Chunhui Zhao
// 2026-10-20).  Pointwise and serial-mesh checks, no time stepping:
//
//   U1  O1 [prestress_depth_table]: the TPV12 effective tensor (plan 10.1)
//       on a 60 deg fault 7.5 km down-dip gives sigma_n' = 8533.25 z_d and
//       a down-dip shear of 4692.03 z_d (55.43 and 30.48 MPa, the DRDG3D
//       t = 0 values); the value just above a jump, at it and below it;
//       linear interpolation between nodes; the table checks refuse
//       non-increasing depths, a jump that is not a node and a missing
//       _below list; CountFaultFacesAcrossDepths finds faces across a depth
//       that is not a mesh line and none across one that is.
//   I1  O1 inactive: a constant table returns the constant tensor bit for bit.
//   U2  O2 SetTensionFloor: at a DOF in tension (sigma_eff = -1 MPa) the law
//       sees 0, so the strength is the cohesion: stick below it, slip at
//       (|X| - c) / gamma above it; the returned sigma_eff and sigma* keep
//       the true value; the floored count grows.  I2: in compression one
//       fault advance with the floor on equals the advance without it, bit
//       for bit.  (The normal-jump scaling sigma / sigma_F of plan U2 needs
//       a time-stepped run and is left to the smoke runs.)
//
// Usage: mpirun -np 1 ./seas_test_nitsche_options

#include "mfem.hpp"
#include "../../dynamic/nitsche_opt_depth_prestress.hpp"
#include "../../dynamic/nitsche_opt_tension.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../friction/slip_weakening_closure.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using namespace mfem::seas::nitsche_opt;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static bool RelClose(real_t a, real_t b, real_t rtol)
{
   return std::abs(a - b) <= rtol * std::max(std::abs(a), std::abs(b));
}

// ---------------------------------------------------------------- U1 / I1

/// TPV12 effective tensor (plan 10.1): z_d < 13800 sin 60 = 11951.15 m:
/// sigma_zz = -16660 z_d, sigma_yy = -0.3496 x 16660 z_d = -5824.336 z_d,
/// sigma_xx = -(26460 + 15624.336) / 2 + 9800 = -11242.168 z_d; below, all
/// three -16660 z_d; no shear components.  The jump depth is the double the
/// mesh line has (tpv12/configs/tpv12_nitsche.toml).
static const real_t kJump = 13800.0 * std::sin(M_PI / 3.0);
static const real_t kSyy = 5824.336, kSxx = 11242.168;

static DepthPrestressTable Tpv12Table()
{
   const real_t dj = kJump, d1 = 20000.0;
   std::ostringstream deck;
   deck << std::setprecision(17)
        << "[prestress_depth_table]\n"
        << "stress_measure = \"effective\"\n"
        << "depth_m        = [0.0, " << dj << ", " << d1 << "]\n"
        << "jump_depths_m  = [" << dj << "]\n"
        << "sigma_xx_pa    = [0.0, " << -kSxx * dj << ", " << -16660.0 * d1 << "]\n"
        << "sigma_yy_pa    = [0.0, " << -kSyy * dj << ", " << -16660.0 * d1 << "]\n"
        << "sigma_zz_pa    = [0.0, " << -16660.0 * dj << ", " << -16660.0 * d1 << "]\n"
        << "sigma_xy_pa    = [0.0, 0.0, 0.0]\n"
        << "sigma_yz_pa    = [0.0, 0.0, 0.0]\n"
        << "sigma_xz_pa    = [0.0, 0.0, 0.0]\n"
        << "sigma_xx_pa_below = [" << -16660.0 * dj << "]\n"
        << "sigma_yy_pa_below = [" << -16660.0 * dj << "]\n"
        << "sigma_zz_pa_below = [" << -16660.0 * dj << "]\n"
        << "sigma_xy_pa_below = [0.0]\n"
        << "sigma_yz_pa_below = [0.0]\n"
        << "sigma_xz_pa_below = [0.0]\n";
   std::istringstream is(deck.str());
   const toml::value root = toml::parse(is, "<tpv12>");
   const real_t up[3] = {0.0, 0.0, 1.0};
   return ParseDepthPrestress(root.at("prestress_depth_table"), up);
}

static void NormalShear(const real_t s[9], const real_t n[3], const real_t d[3],
                        real_t &sig_n, real_t &tau_d)
{
   real_t T[3] = {0.0, 0.0, 0.0};
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) { T[i] += s[3 * i + j] * n[j]; }
   sig_n = -(T[0] * n[0] + T[1] * n[1] + T[2] * n[2]);
   tau_d = T[0] * d[0] + T[1] * d[1] + T[2] * d[2];
}

static void TestU1()
{
   std::cout << "U1: O1 depth table, TPV12 effective tensor\n";
   const DepthPrestressTable tab = Tpv12Table();
   TEST_ASSERT(tab.measure == DepthPrestressTable::Measure::Effective, "stress_measure parsed");
   TEST_ASSERT(tab.Check().empty(), "TPV12 table passes Check()");

   // 60 deg fault dipping toward +Y, trace on the X axis; 7.5 km down-dip
   const real_t dip = M_PI / 3.0, w = 7500.0;
   const real_t x[3] = {0.0, w * std::cos(dip), -w * std::sin(dip)};
   const real_t n[3] = {0.0, -std::sin(dip), -std::cos(dip)};
   const real_t dd[3] = {0.0, std::cos(dip), -std::sin(dip)};
   real_t s[9], sig, tau;
   tab.Eval(x, s);
   NormalShear(s, n, dd, sig, tau);
   const real_t zd = -x[2];
   const real_t sig_ref = (0.75 * kSyy + 0.25 * 16660.0) * zd;
   const real_t tau_ref = (16660.0 - kSyy) * std::sin(dip) * std::cos(dip) * zd;
   std::cout << std::setprecision(10) << "    sigma_n' = " << sig / 1e6 << " MPa, |tau_dip| = "
             << std::abs(tau) / 1e6 << " MPa at z_d = " << zd << " m\n";
   TEST_ASSERT(RelClose(sig, sig_ref, 1e-12), "sigma_n' = 8533.252 z_d to 1e-12");
   TEST_ASSERT(RelClose(std::abs(tau), tau_ref, 1e-12), "|tau_dip| = 4691.98 z_d (0.549847 sigma_n') to 1e-12");
   TEST_ASSERT(std::abs(sig / 1e6 - 55.43) < 0.01 && std::abs(std::abs(tau) / 1e6 - 30.48) < 0.01,
               "55.43 / 30.48 MPa, the DRDG3D t = 0 values at 7.5 km down-dip");

   // jump: above (depth < jump), at it and below (depth >= jump) take their branch
   const real_t ya[3] = {0.0, 0.0, -(kJump - 1.0)}, yj[3] = {0.0, 0.0, -kJump},
                yb[3] = {0.0, 0.0, -(kJump + 1.0)};
   real_t sa[9], sj[9], sb[9];
   tab.Eval(ya, sa);
   tab.Eval(yj, sj);
   tab.Eval(yb, sb);
   TEST_ASSERT(RelClose(sa[0], -kSxx * (kJump - 1.0), 1e-9)
               && RelClose(sa[4], -kSyy * (kJump - 1.0), 1e-9),
               "1 m above the jump: sigma_xx, sigma_yy on the upper branch");
   const real_t vb = -16660.0 * kJump;
   TEST_ASSERT(sj[0] == vb && sj[4] == vb && sj[8] == vb, "at the jump depth: the value below, exactly");
   TEST_ASSERT(RelClose(sb[0], -16660.0 * (kJump + 1.0), 1e-12) && sb[0] == sb[4] && sb[4] == sb[8],
               "1 m below the jump: isotropic -16660 z_d");
   TEST_ASSERT(sa[1] == 0.0 && sa[2] == 0.0 && sa[5] == 0.0 && sb[1] == 0.0,
               "no shear components in the global frame");
   const real_t ym[3] = {0.0, 0.0, -5000.0};
   tab.Eval(ym, s);
   TEST_ASSERT(RelClose(s[8], -16660.0 * 5000.0, 1e-12), "linear between nodes (z_d = 5 km)");
   bool sym = true;
   for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) { sym = sym && s[3 * i + j] == s[3 * j + i]; }
   TEST_ASSERT(sym, "symmetric tensor");

   // Check() refusals
   DepthPrestressTable bad = tab;
   bad.depth = {0.0, 20000.0, kJump};
   TEST_ASSERT(!bad.Check().empty(), "refuses non-increasing depth_m: " + bad.Check());
   bad = tab;
   bad.jump = {12000.0};
   TEST_ASSERT(!bad.Check().empty(), "refuses a jump that is not a node: " + bad.Check());
   bad = tab;
   bad.below[2].clear();
   TEST_ASSERT(!bad.Check().empty(), "refuses a missing _below value: " + bad.Check());
   bad = tab;
   bad.jump = {0.0};
   TEST_ASSERT(!bad.Check().empty(), "refuses a jump at the first node: " + bad.Check());
}

static void TestI1()
{
   std::cout << "I1: O1 inactive (constant table)\n";
   DepthPrestressTable tab;
   tab.measure = DepthPrestressTable::Measure::Total;
   tab.depth = {0.0, 7300.0, 100000.0};
   const real_t c[6] = {-1.1e6, -120.0e6, -3.3e7, 70.0e6, 0.25e6, -0.75e6};
   for (int k = 0; k < 6; k++) { tab.node[k] = {c[k], c[k], c[k]}; }
   TEST_ASSERT(tab.Check().empty(), "constant table passes Check()");
   const real_t full[9] = {c[0], c[3], c[5], c[3], c[1], c[4], c[5], c[4], c[2]};
   bool same = true;
   for (real_t z : {0.0, -0.3, -1.0, -999.9, -7300.0, -12345.678, -99999.0, -100000.0})
   {
      const real_t x[3] = {3.0, -7.0, z};
      real_t s[9];
      tab.Eval(x, s);
      for (int i = 0; i < 9; i++) { same = same && (std::memcmp(&s[i], &full[i], sizeof(real_t)) == 0); }
   }
   TEST_ASSERT(same, "constant table == constant tensor, bit for bit, at 8 depths");
}

// ---------------------------------------------------------------- meshes

static const real_t kL = 3000.0;   // box [-L, L] x [-L, L] x [-L, 0]

/// 24 tets per hex box; free surface Z = 0 (attr 1), walls (attr 5), fault
/// faces on Y = 0 with |X| <= a and Z >= -b (attr 3; a, b on grid lines).
static Mesh *Box(int nx, int ny, int nz, real_t a, real_t b)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(nx, ny, nz, 2 * kL, 2 * kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
   }
   const real_t tol = 1e-6 * kL;
   for (int be = 0; be < mesh->GetNBE(); be++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(be, verts);
      bool top = true;
      for (int i = 0; i < verts.Size(); i++) { top = top && std::abs(mesh->GetVertex(verts[i])[2]) < tol; }
      mesh->SetBdrAttribute(be, top ? 1 : 5);
   }
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool on = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         const real_t *x = mesh->GetVertex(fv[i]);
         on = on && std::abs(x[1]) < tol && std::abs(x[0]) <= a + tol && x[2] >= -b - tol;
      }
      if (on) { tris.push_back({fv[0], fv[1], fv[2]}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   return mesh;
}

static void TestU1Mesh()
{
   std::cout << "U1: O1 faces across a jump depth (serial mesh)\n";
   // 6 x 6 x 3 hexes of 1 km; fault |X| <= 2 km, Z >= -3 km
   std::unique_ptr<Mesh> m(Box(6, 6, 3, 2000.0, 3000.0));
   const real_t up[3] = {0.0, 0.0, 1.0};
   std::string first;
   const long long n_line = CountFaultFacesAcrossDepths(*m, 3, {1000.0, 2000.0}, up, 1e-6, &first);
   const long long n_mid = CountFaultFacesAcrossDepths(*m, 3, {1500.0}, up, 1e-6, &first);
   std::cout << "    across 1 and 2 km (mesh lines): " << n_line << "; across 1.5 km: " << n_mid
             << " (" << first << ")\n";
   TEST_ASSERT(n_line == 0, "no fault face crosses a depth that is a mesh line");
   // each 1 km hex face of the 1-2 km layer is 4 triangles around its centre
   // (1.5 km deep); the two with a top and a bottom corner cross 1.5 km, the
   // top and bottom ones only touch it: 2 x 4 hex faces across |X| <= 2 km
   TEST_ASSERT(n_mid == 8, "8 faces cross a mid-layer depth");
}

// ---------------------------------------------------------------- operator

// TPV205 values (prestress 70 / 120 MPa; mu_s 0.677, mu_d 0.525, Dc 0.4 m)
static const real_t kTau0 = 70.0e6, kSig0 = 120.0e6, kMuS = 0.677, kMuD = 0.525, kDc = 0.40;

/// Operator on the box with the fault on the whole plane Y = 0, linear slip
/// weakening, uniform prestress (sigma_xy, sigma_yy), at rest, dt = 1 ms.
struct OpFixture
{
   std::unique_ptr<ParMesh> pm;
   SlipWeakeningClosure cl;
   std::unique_ptr<NitscheWaveOperator> op;
   Vector u;
   OpFixture(real_t sxy, real_t syy, real_t cohesion)
   {
      std::unique_ptr<Mesh> m(Box(4, 4, 2, kL, kL));
      pm.reset(new ParMesh(MPI_COMM_WORLD, *m));
      NitscheFaultSetup fs;
      fs.closure = &cl;
      fs.prestress = [sxy, syy](const real_t[3], real_t s[9]) {
         for (int i = 0; i < 9; i++) { s[i] = 0.0; }
         s[1] = s[3] = sxy;
         s[4] = syy;
      };
      fs.parameters = [cohesion](const real_t[3], FaultPointState &st) {
         st.mu_s = kMuS; st.mu_d = kMuD; st.Dc = kDc; st.cohesion = cohesion;
      };
      fs.sigma_min = 0.0;
      BoundaryConfig bc;
      bc.fault_attr = 3;
      bc.natural_attrs = {1};
      bc.absorbing_attrs = {5};
      op.reset(new NitscheWaveOperator(*pm, 1, MaterialField::MakeConstant(32043759360.0,
                                       32038120320.0, 2670.0), bc, {}, &fs));
      op->SetTimeStep(1e-3);
      u.SetSize(op->Height());
      u = 0.0;
   }
   int No() const { return op->NumOwnedFaultDOFs(); }
   /// d = 0, V = 0, psi as given.
   void State(real_t psi0)
   {
      const int no = No();
      Vector d(3 * no), V(3 * no), psi(no);
      d = 0.0; V = 0.0; psi = psi0;
      op->SetFaultState(d, V, psi);
   }
   real_t MaxV() const
   {
      real_t m = 0.0;
      for (int i = 0; i < op->SlipRate().Size(); i++) { m = std::max(m, std::abs(op->SlipRate()(i))); }
      MPI_Allreduce(MPI_IN_PLACE, &m, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, MPI_COMM_WORLD);
      return m;
   }
   real_t MinV() const   ///< min over DOFs of |V|
   {
      real_t m = 1e300;
      for (int j = 0; j < No(); j++)
      {
         const real_t *V = &op->SlipRate()(3 * j);
         m = std::min(m, std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]));
      }
      MPI_Allreduce(MPI_IN_PLACE, &m, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, MPI_COMM_WORLD);
      return m;
   }
};

static bool SameVec(const Vector &a, const Vector &b)
{
   int same = (a.Size() == b.Size()) && (a.Size() == 0 ||
              std::memcmp(a.GetData(), b.GetData(), sizeof(real_t) * a.Size()) == 0);
   MPI_Allreduce(MPI_IN_PLACE, &same, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
   return same != 0;
}

static void TestU2()
{
   std::cout << "U2: O2 tension floor\n";
   const real_t c0 = 1.0e6;
   OpFixture F(kTau0, -kSig0, c0);
   NitscheWaveOperator &op = *F.op;
   F.State(0.0);
   op.SetTensionFloor(true);
   TEST_ASSERT(op.HasTensionFloor() && op.TensionFloorCount() == 0, "SetTensionFloor(true)");
   bool stick_ok = true, slip_ok = true, sig_ok = true, pc_ok = true;
   for (int j = 0; j < F.No(); j++)
   {
      const real_t *n = &op.FaultFrames()(9 * j), *dd = n + 3, *ss = n + 6;
      for (real_t Xs : {0.5e6, 2.0e6})
      {
         // T = -sigma n + X_s s with sigma = -1 MPa (tension)
         real_t TS[3], V[3], X[2], sig, sigpc = 3.0e6;
         for (int c = 0; c < 3; c++) { TS[c] = 1.0e6 * n[c] + Xs * ss[c]; }
         op.CloseDOF(j, TS, 0.0, V, sig, X, &sigpc, 0.5, 0.5, 0.0);
         const real_t Vabs = std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
         sig_ok = sig_ok && RelClose(sig, -1.0e6, 1e-12);
         pc_ok = pc_ok && RelClose(sigpc, 0.5 * 3.0e6 + 0.5 * (-1.0e6), 1e-12);
         // with sigma* = 1 MPa > 0 the law sees 1 MPa: strength c0 + mu_s 1 MPa
         const real_t tau_s = kMuS * sigpc + c0;
         if (Xs < tau_s) { stick_ok = stick_ok && Vabs == 0.0; }
         else { slip_ok = slip_ok && RelClose(Vabs, (std::abs(X[1]) - tau_s) / op.ClosureGamma(j), 1e-12); }
         (void)dd;
      }
      // sigma* also negative: the law sees 0, strength = c0
      real_t TS[3], V[3], X[2], sig, sigpc = -2.0e6;
      for (int c = 0; c < 3; c++) { TS[c] = 1.0e6 * n[c] + 2.0e6 * ss[c]; }
      op.CloseDOF(j, TS, 0.0, V, sig, X, &sigpc, 0.5, 0.5, 0.0);
      const real_t Vabs = std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
      slip_ok = slip_ok && RelClose(Vabs, (2.0e6 - c0) / op.ClosureGamma(j), 1e-12);
      pc_ok = pc_ok && sigpc == -1.5e6;
      // no Prakash-Clifton, sigma = -1 MPa: strength c0, stick at 0.5 MPa
      for (int c = 0; c < 3; c++) { TS[c] = 1.0e6 * n[c] + 0.5e6 * ss[c]; }
      op.CloseDOF(j, TS, 0.0, V, sig, X, nullptr, 0.0, 1.0, 0.0);
      stick_ok = stick_ok && V[0] == 0.0 && V[1] == 0.0 && V[2] == 0.0;
   }
   TEST_ASSERT(sig_ok, "the returned sigma_eff keeps the true -1 MPa");
   TEST_ASSERT(pc_ok, "sigma* follows the true sigma_eff (may go negative)");
   TEST_ASSERT(stick_ok, "stick below the strength (cohesion + mu_s max(sigma*, 0))");
   TEST_ASSERT(slip_ok, "slip rate (|X| - c0) / gamma when the law sees 0");
   long long nfl = op.TensionFloorCount();
   MPI_Allreduce(MPI_IN_PLACE, &nfl, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   long long no = F.No();
   MPI_Allreduce(MPI_IN_PLACE, &no, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   TEST_ASSERT(nfl == 2 * no, "floored count: the two calls per DOF whose law sigma was < 0");

   // I2: compression, floor on vs off, one advance, bit for bit
   OpFixture A(kTau0, -kSig0, 0.0), B(kTau0, -kSig0, 0.0);
   A.State(0.0); B.State(0.0);
   B.op->SetTensionFloor(true);
   A.op->LtsAdvanceFault(A.u, 0.0, nullptr);
   B.op->LtsAdvanceFault(B.u, 0.0, nullptr);
   TEST_ASSERT(SameVec(A.op->SlipRate(), B.op->SlipRate()) && SameVec(A.op->Slip(), B.op->Slip())
               && SameVec(A.op->State(), B.op->State())
               && SameVec(A.op->LastSigmaEff(), B.op->LastSigmaEff()),
               "I2: in compression the floor changes no bit");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   std::cout << "=== test_nitsche_options (complex-geometries plan rev 2) ===\n";
   TestU1();
   TestI1();
   TestU1Mesh();
   TestU2();
   std::cout << "\n" << num_passed << " / " << num_tests << " passed, " << num_failed
             << " failed\n";
   return num_failed == 0 ? 0 : 1;
}
