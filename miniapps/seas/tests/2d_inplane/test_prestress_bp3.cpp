// tests/2d_inplane/test_prestress_bp3.cpp — 2D in-plane design rev 3,
// phase E (decision D5, derivation 10.4): BP3's prestress,
// [cycle.prestress].a = "amax" (cycle/prestress.hpp), and the initial state
// of the BP3-QD decks on a 2D mesh.
//
// BP3 Table 1: sigma_n 50 MPa, a0 0.010, a_max 0.025, b 0.015, L 0.008 m,
// V0 1e-6 m/s, f0 0.6, V_init = V_p = 1e-9 m/s, eta = mu / (2 cs).  With
// X = V_init / (2 V0) exp(psi_bg / a_max) ~ 8e8, asinh X = ln 2X to 1e-17,
// so eqs 25-26 give in closed form
//   tau0  = sigma_n [f0 + (a_max - b) ln(V_init / V0)] + eta V_init  (26.546 MPa)
//   psi0  = f0 + (a_max - b - a) ln(V_init / V0)
// and psi0 = f0, theta0 = L / V0 = 8000 s where a = a0 = a_max - b, psi0 =
// psi_bg, theta0 = L / V_init = 8e6 s where a = a_max (design 3.5).
//
//   B1  no mesh, a = "amax", both senses (V_init = +-1e-9 m/s in the dip
//       slot, 0 in the strike slot), a = 0.010, 0.0175, 0.025: tau0 the
//       same at every DOF bit for bit and equal to the closed form (1e-12
//       rel.); anti-parallel to V_init (dip slot -sgn(V) |tau0|, strike slot
//       0); psi0 = the closed form (1e-12); theta0 = 8000 s and 8e6 s (1e-9
//       rel.); the law at (V_init, psi0, the DOF's a) plus eta V_init gives
//       tau0 back (1e-12 rel.)
//   B2  no mesh, a = "local" (the default): bit for bit the expression the
//       function had before the switch (written out below), and Table 1's
//       local column: tau0 = sigma_n [f0 + (a - b) ln(V_init / V0)] + eta
//       V_init (31.727 MPa where a = 0.010) and psi0 = psi_bg everywhere
//   B3  each BP3 deck (config/fd_cycle/bp3_qd_*.toml) on its 1 km mesh
//       (bp3/mesh/bp3_d*_1km.msh, the deck's own geometry at 1 km), p = 1
//       and 2: CheckDimension(2) passes; IcLabel "gold+a_max"; GlobalMaxA =
//       0.025 exactly; tau0 the same at every DOF on every rank; after
//       QDPhase::Setup (one domain solve at t = 0) the kernel's psi equals
//       BuildPrestress's psi0 (1e-12) and its slip rate is V_init at every
//       owned DOF: the dip slot to 1e-9 relative, the strike slot exactly 0
//   B4  the BP3 outputs of the same set-ups at t = 0, written to a scratch
//       folder: the 12 on-fault stations each on a fault segment of its
//       owner rank (LocateOnFaces with segment faces, distance 0); their
//       rows (FdCycleOutput::Sample, 9 columns with station_sigma_n): slip
//       0, log10 |V_dip| = -9 (1e-9), the strike slot 0, tau_dip = sgn(V)
//       sigma_n [f0 + (a_max - b) ln(V_p / V0)] (1e-9 rel.: tau0 less eta
//       V), sigma_n = 50 MPa (1e-9 MPa); the fault length 40 km as
//       BenchmarkOutput's fault area and domain A (1e-9 rel.); the 81
//       profile points on segments; the profile rows: x_d = 0, 500, ..
//       40000 m, slip 0, stress_3 the stations' tau_dip, sigma_n 50 MPa; the
//       8 surface stations at t = 0 (no slip yet, V = V_p on the whole
//       fault: two rigid blocks): displacement 0 and velocity +-V_p / 2 s
//       (s: [cycle.loading].slip_dir, + on the side x > 0, into which n
//       points, - on the other; 0+ and 0- by their `toward`) along e1 =
//       (1, 0), e2 = (0, -1), to 1e-8 of V_p / 2
//   B5  no mesh: LocateOnFaces with segment faces (face_dim = 1) on six
//       collinear segments of 1 km at p = 1, 2, 3, every other one listed
//       from its far end, in two orders: a degree-p polynomial of x2 at
//       the DOFs (GaussLobatto nodes) is interpolated to 1e-12 at points
//       inside, at the ends of the line and at the shared vertices; at a
//       shared vertex the segment with the smaller centre x2 wins in both
//       orders; a point 1 m off the line gets the nearest DOF (its
//       distance), not a segment
//
// Modes: (default) B1-B5;  --refuse N (N = 1..7) must abort with:
//   1 a = "amax" with ic = "scec_fd"        "needs ic = \"gold\""
//   2 a = "maximum"                         "must be \"local\" or \"amax\""
//   3 BuildPrestress, amax, a_max = 0       "needs a_max > 0"
//   4 BuildPrestress, amax, no eta |V|      "needs eta_V_init = true"
//   5 2-D deck, a field with axis2, half2   "takes no axis2 or half2 on a 2-D mesh"
//   6 3-D deck, a field without them        "needs axis2 and half2 on a 3-D mesh"
//   7 a field with axis2 but no half2       "half2 is required"
//   8 B5's line at p = 3 with the interior   "from its node (DOF order != basis order)"
//     DOFs of a segment swapped
//
//   make -f tests/2d_inplane/2d_inplane.mk test-prestress-bp3
//
// Scope (Chunhui Zhao, 2026-10-26).

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/fd_cycle_output.hpp"
#include "../../cycle/bench_output.hpp"

#include <unistd.h>

#include <cmath>
#include <filesystem>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
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

// BP3 Table 1
static const real_t kRho = 2670.0, kMu = 32038120320.0;
static const real_t kSigma = 50.0e6, kA0 = 0.010, kAmax = 0.025, kB = 0.015, kL = 0.008;
static const real_t kV0 = 1.0e-6, kF0 = 0.6, kVp = 1.0e-9;

static real_t Eta() { return 0.5 * std::sqrt(kMu * kRho); }   // mu / (2 cs)

/// FaultParams of N DOFs with the given a, Table 1 otherwise, V_init = (v, 0).
static cycle::FaultParams Params(const std::vector<real_t> &a, real_t v)
{
   const int N = static_cast<int>(a.size());
   cycle::FaultParams par;
   spatial::RateStatePerDOFParams &rs = par.rs;
   for (Vector *x : {&rs.a, &rs.b, &rs.Dc, &rs.V_init, &rs.f_0, &rs.V_0, &rs.eta,
                     &rs.sigma_n_eff})
   {
      x->SetSize(N);
   }
   par.V_init.SetSize(2 * N);
   for (int i = 0; i < N; i++)
   {
      rs.a(i) = a[i];
      rs.b(i) = kB;
      rs.Dc(i) = kL;
      rs.V_init(i) = std::abs(v);
      rs.f_0(i) = kF0;
      rs.V_0(i) = kV0;
      rs.eta(i) = Eta();
      rs.sigma_n_eff(i) = kSigma;
      par.V_init(2 * i) = v;
      par.V_init(2 * i + 1) = 0.0;
   }
   return par;
}

static cycle::PrestressConfig Amax()
{
   cycle::PrestressConfig c;
   c.include_eta = true;
   c.V_background = kVp;
   c.a = "amax";
   c.a_max = kAmax;
   return c;
}

static real_t Rel(real_t x, real_t ref) { return std::abs(x - ref) / std::abs(ref); }

// ------------------------------------------------------------------- B1
static void TestAmax()
{
   if (g_rank == 0) { std::cout << "\n=== B1: a = \"amax\", BP3 Table 1, no mesh ===\n"; }
   const std::vector<real_t> a = {kA0, kA0, 0.0175, kAmax, kAmax};
   const real_t ln = std::log(kVp / kV0);
   const real_t tau_cf = kSigma * (kF0 + (kAmax - kB) * ln) + Eta() * kVp;
   for (const real_t v : {kVp, -kVp})
   {
      const std::string sense = (v > 0.0) ? "V_init +1e-9" : "V_init -1e-9";
      const cycle::FaultParams par = Params(a, v);
      const cycle::Prestress pre = cycle::BuildPrestress(par, Amax());
      const int N = static_cast<int>(a.size());
      bool same = true, anti = true;
      real_t e_tau = 0.0, e_psi = 0.0, e_back = 0.0;
      const real_t tau = -pre.tau0(0) * ((v > 0.0) ? 1.0 : -1.0);
      for (int i = 0; i < N; i++)
      {
         same = same && pre.tau0(2 * i) == pre.tau0(0) && pre.tau0(2 * i + 1) == pre.tau0(1);
         anti = anti && pre.tau0(2 * i) * v < 0.0 && pre.tau0(2 * i + 1) == 0.0;
         e_tau = std::max(e_tau, Rel(std::abs(pre.tau0(2 * i)), tau_cf));
         const real_t psi_cf = kF0 + (kAmax - kB - a[i]) * ln;
         e_psi = std::max(e_psi, std::abs(pre.psi0(i) - psi_cf));
         // the law at (V_init, psi0, the DOF's a) + eta V_init
         const real_t back = kSigma * a[i] * std::asinh(kVp / (2.0 * kV0) *
                                                         std::exp(pre.psi0(i) / a[i]))
                             + Eta() * kVp;
         e_back = std::max(e_back, Rel(back, tau));
      }
      TEST_ASSERT(same, "B1 " << sense << ": tau0 the same at every DOF, bit for bit ("
                  << tau * 1e-6 << " MPa)");
      TEST_ASSERT(e_tau <= 1e-12, "B1 " << sense << ": |tau0| = sigma_n [f0 + (a_max - b) "
                  "ln(V_init / V0)] + eta V_init = " << tau_cf * 1e-6 << " MPa (rel. " << e_tau
                  << ")");
      TEST_ASSERT(anti, "B1 " << sense << ": tau0 anti-parallel to V_init (dip slot "
                  "-sgn(V) |tau0|, strike slot 0)");
      TEST_ASSERT(e_psi <= 1e-12, "B1 " << sense << ": psi0 = f0 + (a_max - b - a) "
                  "ln(V_init / V0) (" << e_psi << ")");
      const real_t th_vw = kL / kV0 * std::exp((pre.psi0(0) - kF0) / kB);
      const real_t th_vs = kL / kV0 * std::exp((pre.psi0(4) - kF0) / kB);
      TEST_ASSERT(Rel(th_vw, 8000.0) <= 1e-9 && Rel(th_vs, 8.0e6) <= 1e-9,
                  "B1 " << sense << ": theta0 = L / V0 = 8000 s where a = a0 ("
                  << th_vw << " s), L / V_init = 8e6 s where a = a_max (" << th_vs << " s)");
      TEST_ASSERT(e_back <= 1e-12, "B1 " << sense << ": the law at (V_init, psi0, local a) "
                  "+ eta V_init gives tau0 back (rel. " << e_back << ")");
   }
}

// ------------------------------------------------------------------- B2
static void TestLocal()
{
   if (g_rank == 0) { std::cout << "\n=== B2: a = \"local\" (default), no mesh ===\n"; }
   const std::vector<real_t> a = {kA0, 0.0175, kAmax};
   const cycle::FaultParams par = Params(a, kVp);
   cycle::PrestressConfig c;
   c.include_eta = true;
   c.V_background = kVp;
   const cycle::Prestress pre = cycle::BuildPrestress(par, c);
   // The function before the switch, statement for statement.
   const spatial::RateStatePerDOFParams &rs = par.rs;
   bool bits = true;
   for (int i = 0; i < 3; i++)
   {
      const real_t Vi[2] = {par.V_init(2 * i), par.V_init(2 * i + 1)};
      const real_t Vi_abs = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);
      const real_t ai = rs.a(i);
      const real_t sigma_n = rs.sigma_n_eff(i);
      const real_t V0 = rs.V_0(i);
      const real_t psi0 = rs.f_0(i) + rs.b(i) * std::log(V0 / c.V_background);
      const real_t e = std::exp(psi0 / ai);
      const real_t tau =
         c.include_eta
         ? sigma_n * ai * std::asinh((Vi_abs / (2.0 * V0)) * e) + rs.eta(i) * Vi_abs
         : sigma_n * ai * std::asinh((Vi_abs / (2.0 * V0)) * e);
      bits = bits && pre.tau0(2 * i) == -tau * Vi[0] / Vi_abs
             && pre.tau0(2 * i + 1) == -tau * Vi[1] / Vi_abs
             && pre.sigma0(i) == sigma_n && pre.psi0(i) == psi0;
   }
   TEST_ASSERT(bits, "B2 tau0, sigma0, psi0 bit for bit the pre-switch expression");
   const real_t ln = std::log(kVp / kV0);
   const real_t tau_vw = kSigma * (kF0 + (kA0 - kB) * ln) + Eta() * kVp;
   const real_t psi_bg = kF0 + kB * std::log(kV0 / kVp);
   TEST_ASSERT(Rel(-pre.tau0(0), tau_vw) <= 1e-12, "B2 tau0 where a = 0.010: " << -pre.tau0(0) * 1e-6
               << " MPa = Table 1's local column " << tau_vw * 1e-6 << " MPa (5.18 MPa above amax's)");
   TEST_ASSERT(pre.psi0(0) == psi_bg && pre.psi0(1) == psi_bg && pre.psi0(2) == psi_bg,
               "B2 psi0 = psi_bg = f0 + b ln(V0 / V_init) = " << psi_bg << " at every DOF");
}

// ------------------------------------------------------------------- B3
static const char *kDecks[] = {"config/fd_cycle/bp3_qd_d30_thrust.toml",
                               "config/fd_cycle/bp3_qd_d30_normal.toml",
                               "config/fd_cycle/bp3_qd_d60_thrust.toml",
                               "config/fd_cycle/bp3_qd_d60_normal.toml",
                               "config/fd_cycle/bp3_qd_d90.toml"
                              };

/// The 1 km mesh of a deck's 25 m mesh: bp3_dNN_25m.msh -> bp3_dNN_1km.msh.
static std::string CoarseMesh(const std::string &path)
{
   const std::string fine = "_25m.msh";
   MFEM_VERIFY(path.size() > fine.size() && path.compare(path.size() - fine.size(),
                                                          fine.size(), fine) == 0,
               "test_prestress_bp3: the deck's mesh " << path << " is not a _25m.msh");
   return path.substr(0, path.size() - fine.size()) + "_1km.msh";
}

static void TestOutputs(cycle::CycleConfig cfg, cycle::QDPhase &qd, ParMesh &pmesh,
                        MPIContext &mpi, const std::string &tag);

static void TestDeck(const std::string &deck, int p, MPIContext &mpi)
{
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   cfg.shared.mesh.path = CoarseMesh(cfg.shared.mesh.path);
   cfg.shared.mesh.order = p;
   std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
   cfg.CheckDimension(smesh->Dimension());
   ParMesh pmesh(mpi.GetComm(), *smesh);
   smesh.reset();
   std::ostringstream tag;
   tag << "B3 " << deck.substr(deck.rfind('/') + 1) << " p = " << p << ": ";
   TEST_ASSERT(pmesh.Dimension() == 2 && cfg.IcLabel() == "gold+a_max",
               tag.str() << "2-D mesh, CheckDimension passed, IC set " << cfg.IcLabel());

   cycle::QDPhase qd(pmesh, mpi, cfg);
   const spatial::SpatialFrictionConfig &s = cfg.shared;
   const MaterialField material = MaterialField::MakeConstant(
      s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
   const spatial::PorePressureSpec no_pp;
   const cycle::FaultParams par = cycle::BuildFaultParams(
      qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
   cycle::PrestressConfig pc = cfg.prestress;
   pc.a_max = cycle::GlobalMaxA(par, mpi.GetComm());
   TEST_ASSERT(pc.a_max == kAmax, tag.str() << "GlobalMaxA = " << pc.a_max << " (0.025)");
   const cycle::Prestress pre = cycle::BuildPrestress(par, pc);
   const int N = par.rs.a.Size();

   // tau0 uniform over every rank: the min and max of each slot agree.
   real_t lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
   for (int i = 0; i < N; i++)
   {
      for (int k = 0; k < 2; k++)
      {
         lo[k] = std::min(lo[k], pre.tau0(2 * i + k));
         hi[k] = std::max(hi[k], pre.tau0(2 * i + k));
      }
   }
   for (int k = 0; k < 2; k++) { lo[k] = mpi.GlobalMin(lo[k]); hi[k] = mpi.GlobalMax(hi[k]); }
   TEST_ASSERT(lo[0] == hi[0] && lo[1] == hi[1], tag.str() << "tau0 = (" << lo[0] * 1e-6
               << ", " << lo[1] * 1e-6 << ") MPa at every DOF");

   qd.Setup(par, pre);
   Vector psi;
   qd.Psi(psi);
   const Vector &V = qd.SlipRate();
   real_t e_psi = 0.0, e_v = 0.0, v_s = 0.0;
   int n_vw = 0, n_vs = 0;
   for (int i = 0; i < N; i++)
   {
      e_psi = std::max(e_psi, std::abs(psi(i) - pre.psi0(i)));
      e_v = std::max(e_v, std::abs(V(2 * i) - par.V_init(2 * i)) / std::abs(par.V_init(2 * i)));
      v_s = std::max(v_s, std::abs(V(2 * i + 1)));
      if (par.rs.a(i) == kA0) { n_vw++; }
      if (par.rs.a(i) == kAmax) { n_vs++; }
   }
   e_psi = mpi.GlobalMax(e_psi);
   e_v = mpi.GlobalMax(e_v);
   v_s = mpi.GlobalMax(v_s);
   int counts[2] = {n_vw, n_vs};
   MPI_Allreduce(MPI_IN_PLACE, counts, 2, MPI_INT, MPI_SUM, mpi.GetComm());
   TEST_ASSERT(counts[0] > 0 && counts[1] > 0, tag.str() << counts[0] << " DOFs with a = a0, "
               << counts[1] << " with a = a_max, " << qd.NumGlobalFaultDOFs() << " in all");
   TEST_ASSERT(e_psi <= 1e-12, tag.str() << "kernel psi = BuildPrestress psi0 (" << e_psi << ")");
   TEST_ASSERT(e_v <= 1e-9 && v_s == 0.0, tag.str() << "kernel slip rate at t = 0: dip slot = "
               "V_init (rel. " << e_v << "), strike slot " << v_s);
   TestOutputs(cfg, qd, pmesh, mpi, tag.str());
}

// ------------------------------------------------------------------- B4
/// The last row of a RowFile (fd_cycle_output.hpp).
static std::vector<real_t> LastRow(const std::string &path, std::vector<real_t> *first = nullptr)
{
   std::ifstream f(path);
   MFEM_VERIFY(f, "test_prestress_bp3: cannot read " << path);
   std::string line;
   std::vector<real_t> last;
   bool have_first = false;
   while (std::getline(f, line))
   {
      if (line.empty() || line[0] == '#') { continue; }
      std::istringstream is(line);
      std::vector<real_t> r;
      real_t v;
      while (is >> v) { r.push_back(v); }
      if (!have_first && first) { *first = r; }
      have_first = true;
      last = r;
   }
   return last;
}

static void TestOutputs(cycle::CycleConfig cfg, cycle::QDPhase &qd, ParMesh &pmesh,
                        MPIContext &mpi, const std::string &tag3)
{
   const std::string tag = "B4" + tag3.substr(2);   // "B3 <deck> p = <p>: "
   const MPI_Comm comm = mpi.GetComm();
   const char *tmp = std::getenv("TMPDIR");
   long long id = (g_rank == 0) ? static_cast<long long>(getpid()) : 0;
   MPI_Bcast(&id, 1, MPI_LONG_LONG, 0, comm);
   const std::string dir = std::string(tmp ? tmp : "/tmp") + "/seas_test_prestress_bp3_out_"
                           + std::to_string(id);
   if (g_rank == 0) { std::filesystem::create_directories(dir); }
   MPI_Barrier(comm);
   cfg.output.dir = dir;
   const cycle::OutputConfig &oc = cfg.output;
   const cycle::BenchmarkConfig &bc = cfg.bench;
   const int ns = static_cast<int>(oc.stations.size());
   TEST_ASSERT(ns == 12 && oc.station_sigma_n && bc.enabled && bc.volume_stations.size() == 8 &&
               bc.volume_dim == 2 && bc.profiles.size() == 1 && bc.profiles[0].points == 81 &&
               bc.profile_sigma_n, tag << "deck outputs: " << ns << " fault stations, "
               << bc.volume_stations.size() << " surface stations, " << bc.profiles.size()
               << " profile");

   // The stations: each on a segment of its owner.
   const std::vector<cycle::StationLocation> locs = qd.LocateStations(oc);
   std::vector<real_t> gmin;
   const std::vector<int> owner = cycle::StationOwners(locs, comm, gmin);
   int exact = 0;
   real_t dmax = 0.0;
   for (int s = 0; s < ns; s++)
   {
      if (owner[s] == g_rank && locs[s].exact) { exact++; }
      dmax = std::max(dmax, gmin[s]);
   }
   MPI_Allreduce(MPI_IN_PLACE, &exact, 1, MPI_INT, MPI_SUM, comm);
   TEST_ASSERT(exact == ns && dmax == 0.0, tag << "stations on fault segments: " << exact
               << " of " << ns << ", largest distance " << dmax << " m");

   // Their rows at t = 0.
   // the deck's sign (a rank may own no fault DOF)
   const real_t sgn = (cfg.fields.V_init_dip > 0.0) ? 1.0 : -1.0;
   const real_t tau_c = sgn * kSigma * (kF0 + (kAmax - kB) * std::log(kVp / kV0)) / 1e6;
   cycle::FdCycleOutput out(oc, comm, locs, "test");
   auto fill = [&qd](cycle::StationFields &f) { qd.FillStationFields(f); };
   const std::vector<real_t> rows = out.Sample(0.0, fill);
   const int nc = out.NumStationColumns();
   real_t e_slip = 0.0, e_lv = 0.0, e_s = 0.0, e_tau = 0.0, e_sn = 0.0;
   if (g_rank == 0)
   {
      for (int s = 0; s < ns; s++)
      {
         const real_t *r = &rows[nc * s];
         e_slip = std::max({e_slip, std::abs(r[1]), std::abs(r[2])});
         e_lv = std::max(e_lv, std::abs(r[4] - std::log10(kVp)));
         e_s = std::max(e_s, std::abs(r[5]));
         e_tau = std::max(e_tau, Rel(r[6], tau_c));
         e_sn = std::max(e_sn, std::abs(r[8] - kSigma / 1e6));
      }
   }
   TEST_ASSERT(g_rank != 0 || (nc == 9 && e_slip == 0.0 && e_lv <= 1e-9 && e_s == 0.0 &&
                               e_tau <= 1e-9 && e_sn <= 1e-9),
               tag << "station rows at t = 0 (" << nc << " columns): slip " << e_slip
               << ", log10 V_dip + 9: " << e_lv << ", tau_strike " << e_s << ", tau_dip vs "
               << tau_c << " MPa rel. " << e_tau << ", sigma_n - 50 MPa: " << e_sn);

   // The fault length, the profiles and the surface stations.
   const cycle::BenchFaultGeometry bg = qd.BenchGeometry(cfg);
   cycle::BenchmarkOutput bench(cfg, comm, pmesh, qd.Domain().GetFESpace(), bg, "test");
   TEST_ASSERT(Rel(bench.FaultArea(), 40.0e3) <= 1e-9 && Rel(bench.AreaA(), 40.0e3) <= 1e-9,
               tag << "fault length " << bench.FaultArea() << " m, domain A "
               << bench.AreaA() << " m (40 km)");
   std::vector<real_t> pmin;
   const std::vector<int> powner = cycle::StationOwners(bg.profile_locs, comm, pmin);
   int pexact = 0;
   real_t pdmax = 0.0;
   for (size_t i = 0; i < bg.profile_locs.size(); i++)
   {
      if (powner[i] == g_rank && bg.profile_locs[i].exact) { pexact++; }
      pdmax = std::max(pdmax, pmin[i]);
   }
   MPI_Allreduce(MPI_IN_PLACE, &pexact, 1, MPI_INT, MPI_SUM, comm);
   TEST_ASSERT(pexact == 81 && pdmax == 0.0, tag << "profile points on fault segments: "
               << pexact << " of 81, largest distance " << pdmax << " m");
   Vector qv;
   qd.Velocity(qv);
   auto vol = [&qd, &qv](const Vector *&u, const Vector *&v) { u = &qd.Displacement(); v = &qv; };
   bench.Rows(0.0, fill, vol);
   bench.Profiles(0.0, kVp, fill, true);
   bench.Flush();
   MPI_Barrier(comm);
   if (g_rank == 0)
   {
      const std::string pre = dir + "/" + oc.prefix;
      std::vector<real_t> xd;
      const std::vector<real_t> sl = LastRow(pre + "_profile_dip_slip_3.txt", &xd);
      const std::vector<real_t> st = LastRow(pre + "_profile_dip_stress_3.txt");
      const std::vector<real_t> sn = LastRow(pre + "_profile_dip_sigma_n.txt");
      real_t e_x = 0.0, p_slip = 0.0, p_tau = 0.0, p_sn = 0.0;
      bool sizes = xd.size() == 83 && sl.size() == 83 && st.size() == 83 && sn.size() == 83;
      for (int i = 0; sizes && i < 81; i++)
      {
         e_x = std::max(e_x, std::abs(xd[2 + i] - 500.0 * i));
         p_slip = std::max(p_slip, std::abs(sl[2 + i]));
         p_tau = std::max(p_tau, Rel(st[2 + i], tau_c));
         p_sn = std::max(p_sn, std::abs(sn[2 + i] - kSigma / 1e6));
      }
      TEST_ASSERT(sizes && e_x <= 1e-9 && p_slip == 0.0 && p_tau <= 1e-9 && p_sn <= 1e-9,
                  tag << "profile rows: x_d - 500 i: " << e_x << " m, slip " << p_slip
                  << ", stress_3 rel. " << p_tau << ", sigma_n - 50 MPa: " << p_sn);
      const std::array<real_t, 3> &sd = cfg.loading.slip_dir;
      real_t e_u = 0.0, e_v = 0.0;
      for (const cycle::StationSpec &v : bc.volume_stations)
      {
         const std::vector<real_t> r = LastRow(pre + "_" + v.name + ".txt");
         const real_t side = (v.x[0] > 0.0 || v.toward[0] > 0.0) ? 1.0 : -1.0;
         const real_t vx = side * 0.5 * kVp * sd[0], vy = side * 0.5 * kVp * sd[1];
         if (r.size() != 5) { e_u = 1e300; continue; }
         e_u = std::max({e_u, std::abs(r[0]), std::abs(r[1]), std::abs(r[2])});
         e_v = std::max({e_v, std::abs(r[3] - vx), std::abs(r[4] + vy)});
      }
      TEST_ASSERT(e_u == 0.0 && e_v <= 1e-8 * 0.5 * kVp, tag << "surface stations at t = 0: "
                  "displacement " << e_u << " m, velocity - (+-V_p / 2 s) " << e_v << " m/s");
      std::filesystem::remove_all(dir);
   }
   else
   {
      // the other ranks' TEST_ASSERT count must match the root's
      num_tests += 2;
   }
}

// ------------------------------------------------------------------- B5
/// Six segments of 1 km on the line x3 = 0 from x2 = 0, nbf = p + 1 DOFs
/// each (the two ends, then the GaussLobatto interior nodes from the first
/// end); odd segments start at their far end; `order` lists the segments.
/// swap: the first segment's two interior DOFs exchanged (p = 3).
static void SegmentLine(int p, const std::vector<int> &order, Vector &x2, Vector &x3,
                        bool swap = false)
{
   const int nbf = p + 1;
   H1_SegmentElement fe(p, BasisType::GaussLobatto);
   const IntegrationRule &nodes = fe.GetNodes();
   x2.SetSize(nbf * static_cast<int>(order.size()));
   x3.SetSize(x2.Size());
   x3 = 0.0;
   for (size_t f = 0; f < order.size(); f++)
   {
      const int k = order[f];
      const real_t a = 1000.0 * ((k % 2) ? k + 1 : k), b = 1000.0 * ((k % 2) ? k : k + 1);
      for (int j = 0; j < nbf; j++)
      {
         x2(nbf * f + j) = a + nodes.IntPoint(j).x * (b - a);
      }
      if (swap && f == 0 && nbf == 4) { std::swap(x2(2), x2(3)); }
   }
}

static real_t Poly(real_t x, int p)
{
   const real_t s = x / 1000.0;
   real_t v = 1.0, m = 1.0;
   for (int k = 1; k <= p; k++) { m *= s; v += (k + 1) * m * (k % 2 ? 1.0 : -0.5); }
   return v;
}

static void TestSegments()
{
   const std::vector<std::array<real_t, 2>> pts = {{{0.0, 0.0}}, {{1234.5, 0.0}},
      {{2000.0, 0.0}}, {{3000.0, 0.0}}, {{4321.0, 0.0}}, {{5999.9, 0.0}}, {{6000.0, 0.0}}};
   const std::vector<std::array<real_t, 2>> off = {{{2500.0, 1.0}}};
   for (int p = 1; p <= 3; p++)
   {
      std::array<real_t, 3> key2[2];
      real_t err = 0.0, worst_val = 0.0;
      bool all_exact = true, off_ok = true;
      for (int o = 0; o < 2; o++)
      {
         const std::vector<int> order = (o == 0) ? std::vector<int> {0, 1, 2, 3, 4, 5}
                                        : std::vector<int> {4, 2, 5, 0, 3, 1};
         Vector x2, x3;
         SegmentLine(p, order, x2, x3);
         Vector val(x2.Size());
         for (int i = 0; i < x2.Size(); i++) { val(i) = Poly(x2(i), p); }
         const std::vector<cycle::StationLocation> L =
            cycle::LocateOnFaces(x2, x3, p + 1, pts, 1);
         for (size_t i = 0; i < pts.size(); i++)
         {
            all_exact = all_exact && L[i].exact && L[i].distance == 0.0;
            const real_t ref = Poly(pts[i][0], p);
            err = std::max(err, std::abs(cycle::FdCycleOutput::Eval(val, L[i]) - ref));
            worst_val = std::max(worst_val, std::abs(ref));
         }
         key2[o] = L[2].key;
         const std::vector<cycle::StationLocation> Lo =
            cycle::LocateOnFaces(x2, x3, p + 1, off, 1);
         real_t dmin = 1e300;
         for (int i = 0; i < x2.Size(); i++)
         {
            dmin = std::min(dmin, std::hypot(x2(i) - off[0][0], x3(i) - off[0][1]));
         }
         off_ok = off_ok && !Lo[0].exact && Lo[0].distance == dmin && Lo[0].nearest >= 0 &&
                  std::hypot(x2(Lo[0].nearest) - off[0][0], x3(Lo[0].nearest) - off[0][1]) == dmin;
      }
      TEST_ASSERT(all_exact && err <= 1e-12 * worst_val, "B5 p = " << p << ": "
                  << pts.size() << " points on segments in both orders, degree-" << p
                  << " polynomial to " << err << " (of " << worst_val << ")");
      TEST_ASSERT(key2[0] == key2[1] && key2[0][1] == 1500.0, "B5 p = " << p << ": the "
                  "vertex x2 = 2000 m takes the segment centred at " << key2[0][1]
                  << " m in both orders (" << key2[1][1] << ")");
      TEST_ASSERT(off_ok, "B5 p = " << p << ": a point 1 m off the line gets the nearest DOF");
   }
}

// --------------------------------------------------------------- refusals
static std::string ReadFile(const std::string &p)
{
   std::ifstream f(p);
   MFEM_VERIFY(f, "cannot read " << p);
   std::stringstream ss;
   ss << f.rdbuf();
   return ss.str();
}

static std::string WriteDeck(const std::string &tag, const std::string &text)
{
   const char *tmp = std::getenv("TMPDIR");
   const std::string path = std::string(tmp ? tmp : "/tmp") + "/seas_test_prestress_bp3_" + tag
                            + "_" + std::to_string(getpid()) + "_" + std::to_string(g_rank) + ".toml";
   std::ofstream f(path);
   f << text;
   MFEM_VERIFY(f, "cannot write " << path);
   return path;
}

/// text with the one line `from` replaced by `to`.
static std::string Swap(std::string text, const std::string &from, const std::string &to)
{
   const size_t at = text.find(from);
   MFEM_VERIFY(at != std::string::npos && text.find(from, at + 1) == std::string::npos,
               "--refuse: the deck line \"" << from << "\" moved");
   return text.replace(at, from.size(), to);
}

static int Refuse(int n)
{
   const std::string bp3 = ReadFile("config/fd_cycle/bp3_qd_d60_thrust.toml");
   const std::string bp5 = ReadFile("config/fd_cycle/bp5_fd_first_event.toml");
   const std::string a_line = "a            = \"amax\"            # BP3: a_max = 0.025 in tau0 (design D5)";
   std::string text;
   switch (n)
   {
      case 1: text = Swap(bp3, "ic                   = \"gold\"", "ic                   = \"scec_fd\""); break;
      case 2: text = Swap(bp3, a_line, "a            = \"maximum\""); break;
      case 3:
      case 4:
      {
         cycle::PrestressConfig c = Amax();
         if (n == 3) { c.a_max = 0.0; }
         else { c.include_eta = false; }
         if (g_rank == 0) { std::cout << "--refuse " << n << ": BuildPrestress\n"; }
         cycle::BuildPrestress(Params({kA0, kAmax}, kVp), c);   // must abort
         if (g_rank == 0) { std::cout << "--refuse " << n << ": accepted\n"; }
         return 0;
      }
      case 5: text = Swap(bp3, "half1   = 15.0e3", "half1   = 15.0e3\naxis2   = [0.8660254037844386, 0.5]\n"
                          "half2   = 1.0e6"); break;
      case 6: text = Swap(bp5, "axis2   = [0.0, 0.0, 1.0]\nhalf1   = 30.0e3\nhalf2   = 6.0e3",
                          "half1   = 30.0e3"); break;
      case 7: text = Swap(bp5, "half2   = 6.0e3\ntaper   = 2.0e3", "taper   = 2.0e3"); break;
      case 8:
      {
         Vector x2, x3;
         SegmentLine(3, {0, 1, 2, 3, 4, 5}, x2, x3, true);
         if (g_rank == 0) { std::cout << "--refuse 8: LocateOnFaces\n"; }
         cycle::LocateOnFaces(x2, x3, 4, {{{500.0, 0.0}}}, 1);   // must abort
         if (g_rank == 0) { std::cout << "--refuse 8: accepted\n"; }
         return 0;
      }
      default: MFEM_ABORT("--refuse N: N = 1..8");
   }
   const std::string p = WriteDeck("refuse" + std::to_string(n), text);
   if (g_rank == 0) { std::cout << "--refuse " << n << ": loading " << p << "\n"; }
   cycle::CycleConfig c = cycle::CycleConfig::Load(p);   // 1, 2, 7 abort here
   std::remove(p.c_str());
   c.CheckDimension(n == 5 ? 2 : 3);                     // 5, 6 abort here
   if (g_rank == 0) { std::cout << "--refuse " << n << ": the deck was accepted\n"; }
   return 0;
}

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   g_rank = mpi.Rank();
   int refuse = 0;
   for (int i = 1; i < argc; i++)
   {
      if (std::strcmp(argv[i], "--refuse") == 0 && i + 1 < argc) { refuse = std::atoi(argv[++i]); }
   }
   if (refuse) { return Refuse(refuse); }
   if (g_rank == 0)
   {
      std::cout << "test_prestress_bp3 (2D in-plane design rev 3, phase E, decision D5) on "
                << mpi.Size() << " rank(s)\n";
   }
   TestAmax();
   TestLocal();
   if (g_rank == 0) { std::cout << "\n=== B5: segment faces (LocateOnFaces, face_dim = 1) ===\n"; }
   TestSegments();
   if (g_rank == 0) { std::cout << "\n=== B3: the BP3 decks on the 1 km meshes ===\n"; }
   MFEMInitializePetsc(&argc, &argv, nullptr, nullptr);
   for (const char *d : kDecks)
   {
      for (int p : {1, 2}) { TestDeck(d, p, mpi); }
   }
   MFEMFinalizePetsc();
   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_tests - num_failed << "/" << num_tests
                << " passed ===\n";
   }
   return num_failed ? 1 : 0;
}
