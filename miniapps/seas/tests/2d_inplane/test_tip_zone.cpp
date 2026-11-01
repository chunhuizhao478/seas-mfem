// tests/2d_inplane/test_tip_zone.cpp — the velocity-strengthening tip-zone
// field, [cycle.fields.a] kind = "tip_zone" (cycle/fault_chains.hpp; Chunhui
// Zhao, 2026-10-29: "let's write this in the document, and prepare to add
// this function"; UCERF 2-D report rev 15, section 3.1).
//
// Meshes built here: the square |x|, |y| <= 20 km in 40 x 40 cells of 1 km
// (diagonal (i, j)-(i+1, j+1)), the sides attr 5, fault faces attr 3 along
// cell edges, the vertices turned by R about z after tagging.  Layouts (in
// the unrotated frame, km):
//   three     A: y = 0, -10 <= x <= 10 (L = 20, two free tips);
//             B: y = 6, -9 <= x <= 0, then the 45 degree diagonal to (4, 10)
//                (L = 9 + 4 sqrt 2, a 45 degree bend);
//             C: y = -1, 1 <= x <= 9 (L = 8): its tip (1, -1) lies 1.41 km
//                from A's middle (0, 0), as the San Jacinto tip lies 4.6 km
//                from the San Andreas;
//   boundary  A, and D: x = -15 from the bottom side (y = -20) to y = -10
//             (one free tip; the end on the outer boundary is not a tip);
//   junction  A and a branch x = 0, 0 <= y <= 5;  loop  A and the closed
//             square (-15, 5)-(-10, 5)-(-10, 10)-(-15, 10);  through  x =
//             -15 across the box (both ends on the outer boundary).
//
//   G1  FaultChains on three, boundary at R = 0 and 30 degrees: chain count,
//       lengths (1e-9 rel.), free tips; the first end is the lexicographically
//       smaller one; Locate at every vertex and edge midpoint of the faults:
//       the chain and d = the distance along it to the nearer free tip
//       (1e-9 m), D measured from its free tip only
//   D1  the deck path (CycleConfig::Load, QDPhase, BuildFaultParams, as
//       seas_fd_cycle_driver) on three at p = 1 and 2, R = 0 and 30 degrees,
//       width 4 km, taper 1.5 km, inside 0.004, tip 0.04: at every owned
//       fault DOF a = the tip-zone value of the analytic d (written out here,
//       1e-12 abs.); with max_fraction = 0.25 W_eff = min(4 km, L / 4)
//   D2  A's DOF at its middle has a = inside although C's tip is 1.41 km
//       away in a straight line (the distance is along the own chain)
//   D3  every regime (tip, ramp, inside) is hit on A and B; the report names
//       3 chains and 6 free tips
//   D4  [cycle.fields.a] without kind and with kind = "rect_taper" give the
//       same a, bit for bit, = RectTaperField at the DOFs (the default path)
//   D5  on a straight chain the tip zone is the rectangle taper it replaces
//       (Chunhui Zhao, 2026-11-01: the velocity-strengthening patch near the
//       tips): on A (L = 20 km) width 5 km, taper 2 km gives, at every owned
//       DOF of A, the a of the rect taper half1 = L / 2 - width = 5 km,
//       taper 2 km, inside 0.004, outside 0.04 (1e-15 abs.)
//   On 2 ranks the cut x = 2.3 km (unrotated) crosses A, B and C.
//
// Modes: (default) G1, D1-D5;  --refuse N (N = 1..10) must abort with:
//   1 a junction                          "three or more fault faces"
//   2 a closed chain                      "form closed chains"
//   3 a chain with both ends on the outer boundary  "has no free tip"
//   4 kind = "tip_zone" on Dc             "\"tip_zone\" is for [cycle.fields.a] only"
//   5 a 3-D deck and mesh                 "needs a 2-D mesh (fault lines with tips)"
//   6 max_fraction 0.1 (C: 0.8 km < taper)  "is shorter than the taper"
//   7 taper > width                       "taper must be in [0, width]"
//   8 an unknown key (half1)              "has an unknown key"
//   9 a point off the fault               "lies on no fault segment"
//  10 max_fraction = 0                    "max_fraction must be > 0"
// --dump-mesh PATH --order P --out PREFIX [--width W --taper T --inside A
//   --tip A --max-fraction F]: the chains of a mesh file (fault attr 3) and,
//   per local fault face, its p + 1 Gauss-Lobatto nodes (the QD operator's
//   fault DOFs) with chain, s, d and the tip-zone value, to PREFIX.<rank>
//   (tests/2d_inplane/tip_zone_check.py compares them with its own chains).
// --dump-deck DECK --out PREFIX: the deck path (QDPhase, BuildFaultParams)
//   on the deck's mesh: x, y, a at every owned fault DOF to PREFIX.<rank>.
//
//   make -f tests/2d_inplane/2d_inplane.mk test-tip-zone
//
// Scope (Chunhui Zhao, 2026-10-29).

#include "mfem.hpp"
#include "frame_box.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_chains.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/scaled_mesh.hpp"

#include <unistd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using frame_box::Rot;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kRho = 2670.0, kMu = 32038120320.0;
static const real_t kL2 = 20.0e3;        // box half-width [m]
static const int kN = 40;                // cells per side (1 km)
static const real_t kW = 4.0e3, kT = 1.5e3, kIn = 0.004, kTip = 0.04;
static const real_t kLB = 9.0e3 + 4.0e3 * std::sqrt(2.0);   // chain B

// ------------------------------------------------------------------- meshes
enum class Layout { Three, Boundary, Junction, Loop, Through };

/// Grid point (i, j) of the unrotated box [m].
static void GridXY(int i, int j, real_t &x, real_t &y)
{
   x = -kL2 + 2.0 * kL2 * i / kN;
   y = -kL2 + 2.0 * kL2 * j / kN;
}

/// The fault paths of a layout as grid points (km offsets: i = 20 + x, j = 20 + y).
static std::vector<std::vector<std::array<int, 2>>> Paths(Layout l)
{
   auto line = [](int i0, int j0, int i1, int j1)
   {
      std::vector<std::array<int, 2>> p;
      const int n = std::max(std::abs(i1 - i0), std::abs(j1 - j0));
      for (int k = 0; k <= n; k++)
      {
         p.push_back({{i0 + (i1 - i0) * k / n, j0 + (j1 - j0) * k / n}});
      }
      return p;
   };
   auto join = [](std::vector<std::array<int, 2>> a, const std::vector<std::array<int, 2>> &b)
   { a.insert(a.end(), b.begin() + 1, b.end()); return a; };
   const auto A = line(10, 20, 30, 20);
   switch (l)
   {
      case Layout::Three:
         return {A, join(line(11, 26, 20, 26), line(20, 26, 24, 30)), line(21, 19, 29, 19)};
      case Layout::Boundary: return {A, line(5, 0, 5, 10)};
      case Layout::Junction: return {A, line(20, 20, 20, 25)};
      case Layout::Loop:
         return {A, join(join(join(line(5, 25, 10, 25), line(10, 25, 10, 30)),
                              line(10, 30, 5, 30)), line(5, 30, 5, 25))};
      case Layout::Through: return {line(5, 0, 5, 40)};
   }
   return {};
}

/// The box with the layout's faults; part = the 2-rank cut x < 2.3 km
/// (unrotated); vertices turned by R.
static Mesh *Box(Layout l, const Rot &r, std::vector<int> &part)
{
   const int n = kN;
   Mesh *m = new Mesh(2, (n + 1) * (n + 1), 2 * n * n, 0, 2);
   for (int j = 0; j <= n; j++)
      for (int i = 0; i <= n; i++)
      {
         real_t x, y;
         GridXY(i, j, x, y);
         m->AddVertex(x, y);
      }
   auto v = [&](int i, int j) { return i + (n + 1) * j; };
   for (int j = 0; j < n; j++)
      for (int i = 0; i < n; i++)
      {
         m->AddTriangle(v(i, j), v(i + 1, j), v(i + 1, j + 1), 1);
         m->AddTriangle(v(i, j), v(i + 1, j + 1), v(i, j + 1), 1);
      }
   for (int k = 0; k < n; k++)
   {
      m->AddBdrSegment(v(k, 0), v(k + 1, 0), 5);
      m->AddBdrSegment(v(k, n), v(k + 1, n), 5);
      m->AddBdrSegment(v(0, k), v(0, k + 1), 5);
      m->AddBdrSegment(v(n, k), v(n, k + 1), 5);
   }
   for (const auto &p : Paths(l))
      for (size_t k = 0; k + 1 < p.size(); k++)
      {
         m->AddBdrSegment(v(p[k][0], p[k][1]), v(p[k + 1][0], p[k + 1][1]), 3);
      }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   part.resize(m->GetNE());
   for (int e = 0; e < m->GetNE(); e++)
   {
      Vector c;
      m->GetElementCenter(e, c);
      part[e] = (c(0) < 2.3e3) ? 0 : 1;
   }
   for (int i = 0; i < m->GetNV(); i++)
   {
      real_t *c = m->GetVertex(i);
      const real_t x[3] = {c[0], c[1], 0.0};
      real_t y[3];
      r.Apply(x, y);
      c[0] = y[0];
      c[1] = y[1];
   }
   return m;
}

static std::unique_ptr<ParMesh> ParBox(Layout l, const Rot &r, MPIContext &mpi)
{
   std::vector<int> part;
   std::unique_ptr<Mesh> m(Box(l, r, part));
   return std::make_unique<ParMesh>(mpi.GetComm(), *m, mpi.Size() == 2 ? part.data() : nullptr);
}

/// R^T x: the unrotated frame.
static void Unrotate(const Rot &r, const real_t x[3], real_t u[3])
{
   for (int i = 0; i < 3; i++) { u[i] = r.R[0][i] * x[0] + r.R[1][i] * x[1] + r.R[2][i] * x[2]; }
}

/// The analytic chain of an unrotated fault point of the three layout:
/// returns false off the faults; else L and d (distance to the nearer tip).
static bool ThreeChain(real_t u, real_t v, real_t &L, real_t &d, char &name)
{
   const real_t tol = 1e-6;
   if (std::abs(v) < tol && std::abs(u) <= 10.0e3 + tol)
   { L = 20.0e3; d = 10.0e3 - std::abs(u); name = 'A'; return true; }
   if (std::abs(v + 1.0e3) < tol && u >= 1.0e3 - tol && u <= 9.0e3 + tol)
   { L = 8.0e3; d = 4.0e3 - std::abs(u - 5.0e3); name = 'C'; return true; }
   real_t s = -1.0;
   if (std::abs(v - 6.0e3) < tol && u >= -9.0e3 - tol && u <= tol) { s = u + 9.0e3; }
   else if (std::abs(v - 6.0e3 - u) < tol && u >= -tol && u <= 4.0e3 + tol)
   { s = 9.0e3 + std::sqrt(2.0) * u; }
   if (s < 0.0) { return false; }
   L = kLB;
   d = std::min(s, L - s);
   name = 'B';
   return true;
}

/// The tip-zone value written out (independent of TipZoneSpec::Value).
static real_t Expected(real_t d, real_t L, real_t max_fraction)
{
   const real_t W = (max_fraction > 0.0) ? std::min(kW, max_fraction * L) : kW;
   if (d <= W - kT) { return kTip; }
   if (d >= W) { return kIn; }
   return kTip + (kIn - kTip) * (d - (W - kT)) / kT;
}

// --------------------------------------------------------------------- deck
static std::string Num(real_t v)
{
   std::ostringstream o;
   o << std::setprecision(17) << v;
   std::string s = o.str();
   if (s.find_first_of(".eE") == std::string::npos) { s += ".0"; }
   return s;
}

static std::string DeckVec(const Rot &r, real_t x, real_t y, real_t z, int dim)
{
   const real_t a[3] = {x, y, z};
   real_t b[3];
   r.Apply(a, b);
   std::string s = "[" + Num(b[0]) + ", " + Num(b[1]);
   if (dim == 3) { s += ", " + Num(b[2]); }
   return s + "]";
}

struct Spec
{
   int dim = 2;
   int p = 1;
   Rot R;
   std::string a_table;   ///< the [cycle.fields.a] body
   std::string dc_extra;  ///< lines appended to [cycle.fields.Dc]
};

static std::string TipZoneTable(real_t max_fraction = 0.0, const std::string &extra = "")
{
   std::ostringstream o;
   o << "kind = \"tip_zone\"\nwidth = " << Num(kW) << "\ntaper = " << Num(kT)
     << "\ninside = " << Num(kIn) << "\ntip = " << Num(kTip) << "\n";
   if (max_fraction > 0.0) { o << "max_fraction = " << Num(max_fraction) << "\n"; }
   return o.str() + extra;
}

static std::string RectTable(const Spec &sp, real_t half, real_t taper, real_t in, real_t out,
                             bool kind)
{
   std::ostringstream o;
   if (kind) { o << "kind = \"rect_taper\"\n"; }
   o << "center = " << DeckVec(sp.R, 0.0, 0.0, 0.0, sp.dim) << "\n"
     << "axis1 = " << DeckVec(sp.R, 1.0, 0.0, 0.0, sp.dim) << "\n"
     << "half1 = " << Num(half) << "\n";
   if (sp.dim == 3)
   {
      o << "axis2 = " << DeckVec(sp.R, 0.0, 0.0, 1.0, 3) << "\nhalf2 = " << Num(half) << "\n";
   }
   if (taper > 0.0) { o << "taper = " << Num(taper) << "\n"; }
   o << "inside = " << Num(in) << "\noutside = " << Num(out) << "\n";
   return o.str();
}

static std::string DeckText(const Spec &sp)
{
   const int d = sp.dim;
   std::ostringstream o;
   o << "[meta]\nschema_version = 1\nlaw = \"rate_state\"\ndescription = \"test_tip_zone\"\n\n"
     << "[material_constant_fallback]\nlambda = " << Num(kMu) << "\nmu = " << Num(kMu)
     << "\nrho = " << Num(kRho) << "\n\n"
     << "[pore_pressure]\nP_p_pa = 0.0\nP_p_grad_pa_per_m = 0.0\nmin_sigma_n_pa = 0.0\n\n"
     << "[mesh]\npath = \"built_in_the_test.msh\"\norder = " << sp.p << "\n\n"
     << "[velocity]\nuse_sidecar = false\n\n"
     << "[stress]\nkind = \"fault_local_prestress\"\ntau_strike_pa = 0.0\n"
     << "tau_dip_pa = 0.0\nsigma_n_pa = 25.0e6\n\n"
     << "[numerics]\ncfl = 0.5\nkv_eta = 0.1\nkv_scale = \"node\"\nkv_length = \"inscribed\"\n\n"
     << "[time]\ntfinal = 1.0e9\ndt_initial = 0.01\n\n"
     << "[output]\noutput_dir = \".\"\n\n"
     << "[boundary]\nfault_attr = 3\nnatural_attrs = [1]\n\n"
     << "[fault_geometry]\nref_normal = " << DeckVec(sp.R, 0.0, -1.0, 0.0, d) << "\n"
     << "up = " << (d == 2 ? std::string("[0.0, 0.0, 1.0]") : DeckVec(sp.R, 0.0, 0.0, 1.0, 3))
     << "\nkind = \"vertical_strike_slip\"\norientation = \"propagated\"\npin_by = \"reference\"\n"
     << "tangent = \"pointwise\"\njunction_policy = \"refuse\"\npin_min_margin = 0.1\n\n"
     << "[friction.rate_state]\na_default = 0.004\nb_default = 0.03\nDc_default = 0.14\n"
     << "V_0_default = 1.0e-6\nf_0_default = 0.6\nV_init_default = 1.0e-9\n"
     << "sigma_n_default = 25.0e6\nstate_evolution = \"aging_law\"\neta = \"auto\"\n\n"
     << "[cycle]\nic = \"scec_fd\"\nswitching = true\n"
     << "qd_radiation_damping = true\nmesh_scale = 1.0\nplate_rate = 1.0e-9\n"
     << "penalty_factor = 1.0\n\n"
     << "[cycle.qd]\nsolver = \"mumps\"\n"
     << "petsc_options = \"tests/verification/petsc_ts_rk45_tandem.cfg\"\n"
     << "ts_atol = 1.0e-7\nts_rtol = 1.0e-50\n\n"
     << "[cycle.prestress]\nV_background = 1.0e-9\n\n"
     << "[cycle.loading]\nkind = \"plate\"\nplane_point = " << DeckVec(sp.R, 0.0, 0.0, 0.0, d) << "\n"
     << "plane_normal = " << DeckVec(sp.R, 0.0, 1.0, 0.0, d) << "\n"
     << "slip_dir = " << DeckVec(sp.R, 1.0, 0.0, 0.0, d) << "\n"
     << "band = " << (d == 2 ? "0.0" : "1.0") << "\n\n"
     << "[cycle.fields]\nV_init_dip = 0.0\n\n"
     << "[cycle.fields.a]\n" << sp.a_table << "\n"
     << "[cycle.fields.Dc]\n" << RectTable(sp, 0.0, 0.0, 0.14, 0.14, false) << sp.dc_extra << "\n"
     << "[cycle.fields.V_init]\n" << RectTable(sp, 0.0, 0.0, 1.0e-9, 1.0e-9, false) << "\n"
     << "[cycle.output]\nprefix = \"t\"\nplane_origin = " << DeckVec(sp.R, 0.0, 0.0, 0.0, d) << "\n"
     << "plane_axis1 = " << DeckVec(sp.R, 1.0, 0.0, 0.0, d) << "\n"
     << "plane_axis2 = " << DeckVec(sp.R, 0.0, d == 2 ? 1.0 : 0.0, d == 2 ? 0.0 : 1.0, d) << "\n";
   return o.str();
}

static cycle::CycleConfig LoadText(const std::string &text)
{
   const char *tmp = std::getenv("TMPDIR");
   const std::string path = std::string(tmp ? tmp : "/tmp") + "/seas_test_tip_zone_"
                            + std::to_string(getpid()) + "_" + std::to_string(g_rank) + ".toml";
   {
      std::ofstream f(path);
      f << text;
      MFEM_VERIFY(f, "cannot write " << path);
   }
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(path);
   std::remove(path.c_str());
   return cfg;
}

/// QDPhase + BuildFaultParams on the deck (the driver's path).
static cycle::FaultParams Build(const cycle::CycleConfig &cfg, ParMesh &pmesh, MPIContext &mpi,
                                std::unique_ptr<cycle::QDPhase> &qd)
{
   qd = std::make_unique<cycle::QDPhase>(pmesh, mpi, cfg);
   const spatial::SpatialFrictionConfig &s = cfg.shared;
   const MaterialField material = MaterialField::MakeConstant(
      s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
   const spatial::PorePressureSpec no_pp;
   return cycle::BuildFaultParams(qd->Table(), *s.rate_state, cfg.fields, material, pmesh,
                                  no_pp, cfg.fault_attr);
}

// ---------------------------------------------------------------------- G1
static void TestChains(Layout l, const std::string &lname, real_t deg, MPIContext &mpi)
{
   const Rot r = Rot::Axis(0, 0, 1, deg);
   std::unique_ptr<ParMesh> pm = ParBox(l, r, mpi);
   cycle::FaultChains ch(*pm, 3);
   const std::string tag = "G1 " + lname + " " + Num(deg) + " deg: ";
   const auto paths = Paths(l);
   TEST_ASSERT(ch.NumChains() == static_cast<int>(paths.size()),
               tag << ch.NumChains() << " chains");
   int bad_len = 0, bad_end = 0, bad_d = 0, npts = 0;
   real_t worst = 0.0;
   for (const auto &p : paths)
   {
      // this path's geometry, unrotated: arc lengths along it
      std::vector<real_t> s(p.size(), 0.0);
      std::vector<std::array<real_t, 2>> xy(p.size());
      for (size_t k = 0; k < p.size(); k++)
      {
         GridXY(p[k][0], p[k][1], xy[k][0], xy[k][1]);
         if (k > 0) { s[k] = s[k - 1] + std::hypot(xy[k][0] - xy[k - 1][0], xy[k][1] - xy[k - 1][1]); }
      }
      const real_t L = s.back();
      // ends on the outer boundary (|x| or |y| = 20 km) are not tips
      auto onb = [](const std::array<real_t, 2> &q)
      { return std::abs(std::abs(q[0]) - kL2) < 1e-6 || std::abs(std::abs(q[1]) - kL2) < 1e-6; };
      const bool free0 = !onb(xy.front()), free1 = !onb(xy.back());
      for (size_t k = 0; k + 1 < 2 * p.size(); k++)   // vertices and edge midpoints
      {
         const size_t a = k / 2, b = std::min(a + 1, p.size() - 1);
         const real_t w = (k % 2 == 0) ? 0.0 : 0.5;
         const real_t u[3] = {(1 - w) * xy[a][0] + w * xy[b][0], (1 - w) * xy[a][1] + w * xy[b][1], 0.0};
         const real_t sa = (1 - w) * s[a] + w * s[b];
         const real_t dw = std::min(free0 ? sa : 1e300, free1 ? L - sa : 1e300);
         real_t x[3];
         r.Apply(u, x);
         x[2] = 0.0;
         int c = -1;
         real_t sc = 0.0, dc = 0.0;
         ch.Locate(x, c, sc, dc);
         npts++;
         if (std::abs(ch.Length(c) - L) > 1e-9 * L) { bad_len++; }
         // the chain's first end is the lexicographically smaller of the two;
         // its free flags are the path's, matched end to end
         const auto &e0 = ch.FirstEnd(c), &e1 = ch.LastEnd(c);
         real_t f0[3];
         const real_t uf[3] = {xy.front()[0], xy.front()[1], 0.0};
         r.Apply(uf, f0);
         const bool front_first = std::hypot(f0[0] - e0[0], f0[1] - e0[1]) < 1e-6;
         const bool flags_ok = front_first ? (ch.FreeFirst(c) == free0 && ch.FreeLast(c) == free1)
                                           : (ch.FreeFirst(c) == free1 && ch.FreeLast(c) == free0);
         if (!(e0 < e1) || !flags_ok) { bad_end++; }
         worst = std::max(worst, std::abs(dc - dw));
         if (std::abs(dc - dw) > 1e-9) { bad_d++; }
      }
   }
   TEST_ASSERT(bad_len == 0 && bad_end == 0, tag << npts << " points: chain length and ends");
   TEST_ASSERT(bad_d == 0, tag << "d along the own chain to the nearer free tip, max error "
               << worst << " m");
   if (l == Layout::Boundary)
   {
      int nfree = 0;
      for (int c = 0; c < ch.NumChains(); c++) { nfree += ch.FreeFirst(c) + ch.FreeLast(c); }
      TEST_ASSERT(nfree == 3, tag << "3 free tips (D's end on the outer boundary is not one)");
   }
}

// ------------------------------------------------------------------- D1-D4
static void TestDeck(int p, real_t deg, real_t max_fraction, MPIContext &mpi)
{
   Spec sp;
   sp.p = p;
   sp.R = Rot::Axis(0, 0, 1, deg);
   sp.a_table = TipZoneTable(max_fraction);
   const std::string tag = "D1 p = " + std::to_string(p) + ", " + Num(deg) + " deg"
                           + (max_fraction > 0.0 ? ", max_fraction " + Num(max_fraction) : "")
                           + ": ";
   cycle::CycleConfig cfg = LoadText(DeckText(sp));
   TEST_ASSERT(cfg.fields.a_tip_zone.set && !cfg.fields.a.set, tag << "deck loads kind = tip_zone");
   std::unique_ptr<ParMesh> pm = ParBox(Layout::Three, sp.R, mpi);
   cfg.CheckDimension(2);
   std::unique_ptr<cycle::QDPhase> qd;
   const cycle::FaultParams par = Build(cfg, *pm, mpi, qd);
   const cycle::FaultDOFTable &t = qd->Table();
   int off = 0, bad = 0, n = 0;
   int regime[3][3] = {{0}};   // chain A, B, C x tip, ramp, inside
   real_t worst = 0.0;
   int mid_ok = 0, mid_seen = 0;
   for (int i = 0; i < t.Size(); i++)
   {
      real_t u[3];
      Unrotate(sp.R, t.coords.GetData() + 3 * i, u);
      real_t L = 0.0, d = 0.0;
      char name = 0;
      if (!ThreeChain(u[0], u[1], L, d, name)) { off++; continue; }
      n++;
      const real_t want = Expected(d, L, max_fraction);
      worst = std::max(worst, std::abs(par.rs.a(i) - want));
      if (std::abs(par.rs.a(i) - want) > 1e-12) { bad++; }
      const int k = (name == 'A') ? 0 : (name == 'B') ? 1 : 2;
      regime[k][(want == kTip) ? 0 : (want == kIn) ? 2 : 1]++;
      if (name == 'A' && std::abs(u[0]) < 1e-6)
      {
         mid_seen++;
         if (par.rs.a(i) == kIn) { mid_ok++; }
      }
   }
   const int gn = mpi.GlobalSumInt(n), goff = mpi.GlobalSumInt(off), gbad = mpi.GlobalSumInt(bad);
   real_t gworst = 0.0;
   MPI_Allreduce(&worst, &gworst, 1, MPI_DOUBLE, MPI_MAX, mpi.GetComm());
   TEST_ASSERT(goff == 0 && gbad == 0 && gn > 0, tag << gn << " owned fault DOFs on A, B, C; a = "
               "the tip-zone value of the analytic d, max |error| " << gworst);
   if (max_fraction == 0.0)
   {
      const int gmid = mpi.GlobalSumInt(mid_seen), gmid_ok = mpi.GlobalSumInt(mid_ok);
      TEST_ASSERT(gmid > 0 && gmid_ok == gmid, "D2 p = " << p << ", " << Num(deg) << " deg: A's "
                  << gmid << " DOF(s) at its middle have a = inside; C's tip is "
                  << std::sqrt(2.0) << " km away in a straight line, < W - taper");
      int g[3][3];
      MPI_Allreduce(&regime[0][0], &g[0][0], 9, MPI_INT, MPI_SUM, mpi.GetComm());
      bool all = true;
      for (int k = 0; k < 2; k++)
         for (int m = 0; m < 3; m++) { all = all && g[k][m] > 0; }
      TEST_ASSERT(all, "D3 p = " << p << ", " << Num(deg) << " deg: tip / ramp / inside DOFs "
                  "on A " << g[0][0] << " / " << g[0][1] << " / " << g[0][2] << ", on B "
                  << g[1][0] << " / " << g[1][1] << " / " << g[1][2] << ", on C " << g[2][0]
                  << " / " << g[2][1] << " / " << g[2][2]);
      TEST_ASSERT(par.a_report.find("3 chains") != std::string::npos &&
                  par.a_report.find("6 free tips") != std::string::npos,
                  "D3 p = " << p << ", " << Num(deg) << " deg: the report: "
                  << par.a_report.substr(0, par.a_report.find('\n')));
   }
}

static void TestDefaultPath(int p, MPIContext &mpi)
{
   Spec sp;
   sp.p = p;
   sp.R = Rot::Axis(0, 0, 1, 30.0);
   std::vector<real_t> a[2];
   cycle::RectTaperSpec rect;
   for (int k = 0; k < 2; k++)
   {
      sp.a_table = RectTable(sp, 5.0e3, 2.0e3, kIn, kTip, k == 1);
      cycle::CycleConfig cfg = LoadText(DeckText(sp));
      rect = cfg.fields.a;
      std::unique_ptr<ParMesh> pm = ParBox(Layout::Three, sp.R, mpi);
      std::unique_ptr<cycle::QDPhase> qd;
      const cycle::FaultParams par = Build(cfg, *pm, mpi, qd);
      int bad = 0;
      for (int i = 0; i < qd->Table().Size(); i++)
      {
         a[k].push_back(par.rs.a(i));
         if (par.rs.a(i) != cycle::RectTaperField(rect, qd->Table().coords.GetData() + 3 * i)) { bad++; }
      }
      TEST_ASSERT(!cfg.fields.a_tip_zone.set && cfg.fields.a.set && mpi.GlobalSumInt(bad) == 0
                  && par.a_report.empty(),
                  "D4 p = " << p << ": " << (k ? "kind = \"rect_taper\"" : "no kind")
                  << ": a = RectTaperField at every DOF, bit for bit");
   }
   int same = (a[0] == a[1]) ? 1 : 0, all = 0;
   MPI_Allreduce(&same, &all, 1, MPI_INT, MPI_MIN, mpi.GetComm());
   TEST_ASSERT(all == 1, "D4 p = " << p << ": no kind and kind = \"rect_taper\" give the same a");
}

// ---------------------------------------------------------------------- D5
static void TestRectEquivalent(int p, MPIContext &mpi)
{
   Spec sp;
   sp.p = p;
   sp.R = Rot::Axis(0, 0, 1, 30.0);
   const real_t W = 5.0e3, T = 2.0e3, LA = 20.0e3;
   std::vector<real_t> a[2];
   std::vector<char> on_a;
   for (int k = 0; k < 2; k++)
   {
      if (k == 0)
      {
         std::ostringstream o;
         o << "kind = \"tip_zone\"\nwidth = " << Num(W) << "\ntaper = " << Num(T)
           << "\ninside = " << Num(kIn) << "\ntip = " << Num(kTip) << "\n";
         sp.a_table = o.str();
      }
      else { sp.a_table = RectTable(sp, 0.5 * LA - W, T, kIn, kTip, true); }
      cycle::CycleConfig cfg = LoadText(DeckText(sp));
      std::unique_ptr<ParMesh> pm = ParBox(Layout::Three, sp.R, mpi);
      std::unique_ptr<cycle::QDPhase> qd;
      const cycle::FaultParams par = Build(cfg, *pm, mpi, qd);
      for (int i = 0; i < qd->Table().Size(); i++)
      {
         a[k].push_back(par.rs.a(i));
         if (k == 0)
         {
            real_t u[3];
            Unrotate(sp.R, qd->Table().coords.GetData() + 3 * i, u);
            real_t L = 0.0, d = 0.0;
            char name = 0;
            on_a.push_back(ThreeChain(u[0], u[1], L, d, name) && name == 'A');
         }
      }
   }
   int n = 0, bad = 0, nt = 0, nr = 0, ni = 0;
   real_t worst = 0.0;
   const bool same_size = a[0].size() == a[1].size();
   for (size_t i = 0; same_size && i < a[0].size(); i++)
   {
      if (!on_a[i]) { continue; }
      n++;
      const real_t e = std::abs(a[0][i] - a[1][i]);
      worst = std::max(worst, e);
      if (e > 1e-15) { bad++; }
      if (a[1][i] == kTip) { nt++; } else if (a[1][i] == kIn) { ni++; } else { nr++; }
   }
   real_t gworst = 0.0;
   MPI_Allreduce(&worst, &gworst, 1, MPI_DOUBLE, MPI_MAX, mpi.GetComm());
   const int gn = mpi.GlobalSumInt(n), gbad = mpi.GlobalSumInt(bad) + mpi.GlobalSumInt(!same_size);
   const int gt = mpi.GlobalSumInt(nt), gr = mpi.GlobalSumInt(nr), gi = mpi.GlobalSumInt(ni);
   TEST_ASSERT(gn > 0 && gbad == 0 && gt > 0 && gr > 0 && gi > 0,
               "D5 p = " << p << ": on A's " << gn << " owned DOFs (tip / ramp / inside " << gt
               << " / " << gr << " / " << gi << ") tip zone width 5 km, taper 2 km = the rect "
               "taper half1 5 km, taper 2 km, max |difference| " << gworst);
}

// ----------------------------------------------------------------- refusals
static int Refuse(int n, MPIContext &mpi)
{
   Spec sp;
   sp.a_table = TipZoneTable();
   Layout l = Layout::Three;
   switch (n)
   {
      case 1: case 2: case 3:
      {
         l = (n == 1) ? Layout::Junction : (n == 2) ? Layout::Loop : Layout::Through;
         std::unique_ptr<ParMesh> pm = ParBox(l, Rot(), mpi);
         cycle::FaultChains ch(*pm, 3);   // aborts here
         return 0;
      }
      case 4: sp.dc_extra = "kind = \"tip_zone\"\n"; break;
      case 5: sp.dim = 3; break;
      case 6: sp.a_table = TipZoneTable(0.1); break;
      case 7: sp.a_table = "kind = \"tip_zone\"\nwidth = 1000.0\ntaper = 1500.0\ninside = 0.004\n"
                           "tip = 0.04\n"; break;
      case 8: sp.a_table = TipZoneTable(0.0, "half1 = 1.0\n"); break;
      case 9:
      {
         std::unique_ptr<ParMesh> pm = ParBox(l, Rot(), mpi);
         cycle::FaultChains ch(*pm, 3);
         const real_t x[3] = {0.0, 3.0e3, 0.0};
         int c;
         real_t s, d;
         ch.Locate(x, c, s, d);   // aborts here
         return 0;
      }
      case 10: sp.a_table = TipZoneTable(0.0, "max_fraction = 0.0\n"); break;
      default: MFEM_ABORT("--refuse N: N = 1..10");
   }
   if (g_rank == 0) { std::cout << "--refuse " << n << ": loading the deck\n"; }
   cycle::CycleConfig cfg = LoadText(DeckText(sp));   // 4, 7, 8, 10 abort here
   cfg.CheckDimension(sp.dim);                         // 5 aborts here
   std::unique_ptr<ParMesh> pm = ParBox(l, Rot(), mpi);
   std::unique_ptr<cycle::QDPhase> qd;
   Build(cfg, *pm, mpi, qd);                           // 6 aborts here
   if (g_rank == 0) { std::cout << "--refuse " << n << ": accepted\n"; }
   return 0;
}

// -------------------------------------------------------------------- dumps
static std::string RankFile(const std::string &prefix)
{
   return prefix + "." + std::to_string(g_rank);
}

/// The chains of a mesh file and the tip-zone value at every local fault
/// face's p + 1 Gauss-Lobatto nodes.
static int DumpMesh(const std::string &path, int p, const std::string &prefix,
                    const cycle::TipZoneSpec &tz, MPIContext &mpi)
{
   std::unique_ptr<Mesh> sm = cycle::LoadScaledMesh(path, 1.0);
   ParMesh pm(mpi.GetComm(), *sm);
   sm.reset();
   cycle::FaultChains ch(pm, 3);
   cycle::CheckTipZone(tz, ch, "a");
   if (g_rank == 0) { std::cout << "chains: " << ch.Report(&tz); }
   IntegrationRule nodes;
   QuadratureFunctions1D::GaussLobatto(p + 1, &nodes);
   std::ofstream f(RankFile(prefix));
   f << std::setprecision(17);
   int rows = 0;
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      if (pm.GetBdrAttribute(be) != 3) { continue; }
      ElementTransformation *T = pm.GetBdrElementTransformation(be);
      for (int q = 0; q < nodes.GetNPoints(); q++)
      {
         Vector X;
         T->Transform(nodes.IntPoint(q), X);
         const real_t x[3] = {X(0), X(1), 0.0};
         int c = -1;
         real_t s = 0.0, d = 0.0;
         ch.Locate(x, c, s, d);
         f << x[0] << " " << x[1] << " " << c << " " << s << " " << d << " "
           << tz.Value(d, ch.Length(c)) << "\n";
         rows++;
      }
   }
   const int g = mpi.GlobalSumInt(rows);
   if (g_rank == 0) { std::cout << "dump-mesh: " << g << " rows (p = " << p << ", " << mpi.Size()
                                << " rank(s)) to " << prefix << ".<rank>\n"; }
   return 0;
}

/// The deck path on the deck's own mesh: x, y, a at every owned fault DOF.
static int DumpDeck(const std::string &deck, const std::string &prefix, MPIContext &mpi)
{
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   std::unique_ptr<Mesh> sm = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
   cfg.CheckDimension(sm->Dimension());
   ParMesh pm(mpi.GetComm(), *sm);
   sm.reset();
   std::unique_ptr<cycle::QDPhase> qd;
   const cycle::FaultParams par = Build(cfg, pm, mpi, qd);
   if (g_rank == 0) { std::cout << "chains: " << par.a_report; }
   std::ofstream f(RankFile(prefix));
   f << std::setprecision(17);
   const cycle::FaultDOFTable &t = qd->Table();
   for (int i = 0; i < t.Size(); i++)
   {
      f << t.coords(3 * i) << " " << t.coords(3 * i + 1) << " " << par.rs.a(i) << "\n";
   }
   const int g = mpi.GlobalSumInt(t.Size());
   if (g_rank == 0) { std::cout << "dump-deck: " << g << " owned fault DOFs (p = "
                                << cfg.shared.mesh.order << ", " << mpi.Size()
                                << " rank(s)) to " << prefix << ".<rank>\n"; }
   return 0;
}

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   g_rank = mpi.Rank();
   int refuse = 0, order = 1;
   std::string mesh, deck, out;
   cycle::TipZoneSpec tz;
   tz.width = 20.0e3;
   tz.taper = 2.0e3;
   tz.inside = 0.004;
   tz.tip = 0.04;
   tz.set = true;
   for (int i = 1; i < argc; i++)
   {
      const std::string a = argv[i];
      const bool more = i + 1 < argc;
      if (a == "--refuse" && more) { refuse = std::atoi(argv[++i]); }
      else if (a == "--dump-mesh" && more) { mesh = argv[++i]; }
      else if (a == "--dump-deck" && more) { deck = argv[++i]; }
      else if (a == "--order" && more) { order = std::atoi(argv[++i]); }
      else if (a == "--out" && more) { out = argv[++i]; }
      else if (a == "--width" && more) { tz.width = std::atof(argv[++i]); }
      else if (a == "--taper" && more) { tz.taper = std::atof(argv[++i]); }
      else if (a == "--inside" && more) { tz.inside = std::atof(argv[++i]); }
      else if (a == "--tip" && more) { tz.tip = std::atof(argv[++i]); }
      else if (a == "--max-fraction" && more) { tz.max_fraction = std::atof(argv[++i]); }
      else { MFEM_ABORT("unknown argument " << a); }
   }
   MFEMInitializePetsc(&argc, &argv, nullptr, nullptr);
   int rc = 0;
   if (refuse) { rc = Refuse(refuse, mpi); }
   else if (!mesh.empty())
   {
      MFEM_VERIFY(!out.empty(), "--dump-mesh needs --out PREFIX");
      tz.Validate("a");
      rc = DumpMesh(mesh, order, out, tz, mpi);
   }
   else if (!deck.empty())
   {
      MFEM_VERIFY(!out.empty(), "--dump-deck needs --out PREFIX");
      rc = DumpDeck(deck, out, mpi);
   }
   else
   {
      if (g_rank == 0)
      {
         std::cout << "test_tip_zone ([cycle.fields.a] kind = \"tip_zone\") on " << mpi.Size()
                   << " rank(s)\n";
      }
      for (real_t deg : {0.0, 30.0})
      {
         TestChains(Layout::Three, "three", deg, mpi);
         TestChains(Layout::Boundary, "boundary", deg, mpi);
      }
      for (int p : {1, 2})
      {
         for (real_t deg : {0.0, 30.0}) { TestDeck(p, deg, 0.0, mpi); }
         TestDeck(p, 30.0, 0.25, mpi);
         TestDefaultPath(p, mpi);
         TestRectEquivalent(p, mpi);
      }
      if (g_rank == 0)
      {
         std::cout << "\n=== Summary: " << num_tests - num_failed << "/" << num_tests
                   << " passed ===\n";
      }
      rc = (num_failed == 0) ? 0 : 1;
   }
   MFEMFinalizePetsc();
   return rc;
}
