// tests/2d_inplane/test_plate_loading.cpp — 2D in-plane design rev 3,
// phase A (A4): the fully dynamic cycle's Dirichlet loading [cycle.loading]
// (cycle/loading.hpp), one function for the QD and the dynamic regime.
//
//   L1  kind "plate" with BP5's vectors (p = 0, m = (0, 1, 0), s = (1, 0, 0),
//       b = 1 km; ref (0, -1, 0), band_sign +1) equals MakeBP5DirichletFunc
//       bit for bit at every vertex and quadrature point of the Dirichlet
//       boundary elements (attr 5: walls and creep plane) of the BP5 mesh
//       bp5/mesh/reference/bp5_tandem_exact.msh (x 1000, as the cycle decks
//       scale it), and at the band edges y = +-1000 m and their floating-point
//       neighbours, at t = 0, 1 s, 1 yr and 1800 yr
//   L2  (m, s) -> (-m, -s), the same plate motion, gives the same values bit
//       for bit
//   L3  rotation: with p, m, s and ref rotated by R (40 degrees about
//       (1, 2, 3)), u_R(R x) = R u(x) to 1e-14 relative (points more than
//       1e-6 b from the band edges)
//   L4  u(x, 0) = 0 and u(x, 2 t) = 2 u(x, t) exactly (the premise that
//       QDPhase::CheckLinearLoading_ checks on the assembled right-hand side)
//   L5  the band value agrees with the operator's prescribed-slip sign: on
//       the test box (frame_box.hpp: fault patch and coplanar creep band),
//       flat, dipping 60 degrees and generally rotated, with the walls and
//       the creep band loaded by the plate function and the fault slipping
//       by the same jump (u+ - u- = band_sign Vp t s in the field's frame),
//       the QD operator's solution is the rigid two-block motion
//       u = sgn((x - p) . m) (Vp t / 2) s to 1e-9 of Vp t / 2 at every
//       element centre.  With band_sign reversed the error is O(1) (> 0.1):
//       the factor sgn(-m . ref) is needed
//   L6  CycleConfig::Load of the first-event deck with a [cycle.loading]
//       table: BP5's vectors give band_sign +1; plane_normal left out takes
//       ref_normal (0, -1, 0) and gives band_sign -1, with the same plate
//       motion when slip_dir is reversed (bit for bit with BP5's values)
//   L7  a 2-D point (phase E): with BP3's vectors (dip 60, thrust: p = 0,
//       m = (sin 60, cos 60), s = (-cos 60, sin 60), b = 1 km, z = 0), the
//       value at (x, y) has 2 entries, equal bit for bit to the first two
//       of the value at (x, y, 0), at points beyond, on and inside the band,
//       at t = 0, 1 s and 1500 yr
//
// Modes: (default) L1-L7;  --refuse N (N = 1..6) writes a deck that
// CycleConfig::Load must refuse and loads it (the makefile checks the exit
// status and the message):
//   1 kind "bp5" with another key        "takes no other key"
//   2 slip_dir not a unit vector         "must be unit vectors"
//   3 slip_dir not in the plane          "slip_dir must lie in the plane"
//   4 plane_normal in the fault plane    "lies in the fault plane"
//   5 default "bp5", ref_normal (0,1,0)  "needs [fault_geometry].ref_normal"
//   6 kind "slab"                        "must be \"bp5\" or \"plate\""
//
//   make -f tests/2d_inplane/2d_inplane.mk test-plate-loading
//
// Scope (Chunhui Zhao, 2026-10-25).

#include "mfem.hpp"
#include "../../cycle/loading.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_frame_field.hpp"
#include "frame_box.hpp"

#include <unistd.h>

#include <array>
#include <cmath>
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
using frame_box::FrameBox;
using frame_box::kL;
using frame_box::Rot;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kLambda = 32.04e9, kMu = 32.04e9;
static const real_t kYear = 3.15576e7;
static const char *kDeck = "config/fd_cycle/bp5_fd_first_event.toml";
static const char *kBP5Mesh = "bp5/mesh/reference/bp5_tandem_exact.msh";

static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }

static cycle::LoadingConfig BP5Plate()
{
   cycle::LoadingConfig L;
   L.kind = "plate";
   L.plane_point = {{0.0, 0.0, 0.0}};
   L.plane_normal = {{0.0, 1.0, 0.0}};
   L.slip_dir = {{1.0, 0.0, 0.0}};
   L.band = 1000.0;
   L.band_sign = 1.0;   // sgn(-m . ref), ref = (0, -1, 0)
   return L;
}

static bool SameBits(const Vector &a, const Vector &b)
{
   return a.Size() == b.Size() && std::memcmp(a.GetData(), b.GetData(), sizeof(real_t) * a.Size()) == 0;
}

/// Every vertex and quadrature point of the attr-5 boundary elements of the
/// BP5 mesh (x 1000), and the band-edge points.
static std::vector<std::array<real_t, 3>> BP5Points()
{
   std::vector<std::array<real_t, 3>> pts;
   Mesh mesh(kBP5Mesh, 1, 1);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      real_t *c = mesh.GetVertex(v);
      for (int d = 0; d < 3; d++) { c[d] *= 1000.0; }
   }
   int nbe5 = 0;
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      if (mesh.GetBdrAttribute(b) != 5) { continue; }
      nbe5++;
      Array<int> verts;
      mesh.GetBdrElementVertices(b, verts);
      for (int i = 0; i < verts.Size(); i++)
      {
         const real_t *x = mesh.GetVertex(verts[i]);
         pts.push_back({{x[0], x[1], x[2]}});
      }
      ElementTransformation *T = mesh.GetBdrElementTransformation(b);
      const IntegrationRule &ir = IntRules.Get(mesh.GetBdrElementGeometry(b), 5);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         Vector x(3);
         T->Transform(ir.IntPoint(q), x);
         pts.push_back({{x(0), x(1), x(2)}});
      }
   }
   const real_t edges[] = {1000.0, -1000.0};
   for (real_t e : edges)
   {
      for (real_t y : {e, std::nextafter(e, 1e300), std::nextafter(e, -1e300)})
      {
         pts.push_back({{-12.5e3, y, -7.0e3}});
      }
   }
   pts.push_back({{3.0e3, 0.0, -2.0e3}});
   pts.push_back({{3.0e3, -0.0, -2.0e3}});
   pts.push_back({{3.0e3, 1e-300, -2.0e3}});
   if (g_rank == 0)
   {
      std::cout << "  BP5 mesh: " << nbe5 << " Dirichlet (attr 5) boundary elements, "
                << pts.size() << " points\n";
   }
   return pts;
}

static void TestFunction()
{
   if (g_rank == 0) { std::cout << "\n=== L1-L4: the plate function ===\n"; }
   const real_t Vp = 1.0e-9;
   const std::vector<real_t> times = {0.0, 1.0, kYear, 1800.0 * kYear};
   const std::vector<std::array<real_t, 3>> pts = BP5Points();
   const cycle::LoadingConfig L = BP5Plate();
   const DirichletFunc bp5 = MakeBP5DirichletFunc(Vp);
   const DirichletFunc plate = cycle::MakePlateLoading(L, Vp);
   cycle::LoadingConfig Lneg = L;
   for (int d = 0; d < 3; d++) { Lneg.plane_normal[d] = -L.plane_normal[d]; Lneg.slip_dir[d] = -L.slip_dir[d]; }
   Lneg.band_sign = -1.0;   // sgn(-(-m) . ref)
   const DirichletFunc plate_neg = cycle::MakePlateLoading(Lneg, Vp);

   long n1 = 0, bad1 = 0, bad2 = 0, bad4 = 0;
   Vector x(3), a, b, c, a2;
   for (const auto &p : pts)
   {
      x(0) = p[0]; x(1) = p[1]; x(2) = p[2];
      for (real_t t : times)
      {
         bp5(x, t, a);
         plate(x, t, b);
         plate_neg(x, t, c);
         n1++;
         if (!SameBits(a, b)) { bad1++; }
         if (!SameBits(b, c)) { bad2++; }
         plate(x, 2.0 * t, a2);
         for (int d = 0; d < 3; d++)
         {
            if (a2(d) != 2.0 * b(d) || (t == 0.0 && b(d) != 0.0)) { bad4++; }
         }
      }
   }
   TEST_ASSERT(bad1 == 0, "L1 \"plate\" with BP5's vectors = MakeBP5DirichletFunc bit for bit ("
               << n1 << " point-times, " << bad1 << " differ)");
   TEST_ASSERT(bad2 == 0, "L2 (m, s) -> (-m, -s): the same values bit for bit (" << bad2 << " differ)");
   TEST_ASSERT(bad4 == 0, "L4 u(x, 0) = 0 and u(x, 2t) = 2 u(x, t) exactly (" << bad4 << " differ)");

   // L3: rotated vectors
   const Rot R = Rot::Axis(1, 2, 3, 40.0);
   cycle::LoadingConfig LR = L;
   {
      real_t m0[3] = {0, 1, 0}, s0[3] = {1, 0, 0}, mR[3], sR[3];
      R.Apply(m0, mR);
      R.Apply(s0, sR);
      for (int d = 0; d < 3; d++) { LR.plane_normal[d] = mR[d]; LR.slip_dir[d] = sR[d]; }
      LR.band_sign = 1.0;   // m . ref = R(0,1,0) . R(0,-1,0) = -1
   }
   const DirichletFunc plateR = cycle::MakePlateLoading(LR, Vp);
   real_t err3 = 0.0;
   long n3 = 0;
   for (const auto &p : pts)
   {
      if (std::abs(std::abs(p[1]) - L.band) < 1e-6 * L.band) { continue; }
      real_t xr[3];
      R.Apply(p.data(), xr);
      x(0) = p[0]; x(1) = p[1]; x(2) = p[2];
      Vector xR(xr, 3);
      const real_t t = 1800.0 * kYear;
      plate(x, t, a);
      plateR(xR, t, b);
      real_t Ra[3];
      R.Apply(a.GetData(), Ra);
      for (int d = 0; d < 3; d++) { err3 = std::max(err3, std::abs(Ra[d] - b(d)) / (Vp * t)); }
      n3++;
   }
   TEST_ASSERT(err3 <= 1e-14, "L3 rotated vectors: u_R(R x) = R u(x) (" << err3 << " relative, "
               << n3 << " points)");
}

// ---------------------------------------------------------------- L5
using Op = ElasticityDomainOperator<ParMesh>;

static real_t RigidError(MPI_Comm comm, const Rot &rot, const real_t up[3], real_t band_sign_factor)
{
   std::unique_ptr<Mesh> m(FrameBox(0, rot));
   ParMesh pm(comm, *m);
   const real_t Vp = 2.0, t = 1.0;   // |u| = Vp t / 2 = 1 m
   real_t m0[3] = {0, 1, 0}, s0[3] = {1, 0.0, 0}, r0[3] = {0, -1, 0}, mv[3], sv[3], ref[3];
   rot.Apply(m0, mv);
   rot.Apply(s0, sv);
   rot.Apply(r0, ref);
   cycle::LoadingConfig L;
   L.kind = "plate";
   for (int d = 0; d < 3; d++) { L.plane_normal[d] = mv[d]; L.slip_dir[d] = sv[d]; }
   L.band = 1e-3;   // the plane only: the walls are beyond it
   const real_t mr = mv[0] * ref[0] + mv[1] * ref[1] + mv[2] * ref[2];
   L.band_sign = ((mr < 0.0) ? 1.0 : -1.0) * band_sign_factor;

   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5, 6};
   bc.default_dirichlet_func = cycle::MakePlateLoading(L, Vp);
   auto fc = std::make_shared<FaultFrameConfig>();
   for (int d = 0; d < 3; d++)
   {
      fc->orientation.ref_normal[d] = ref[d];
      fc->orientation.up[d] = up[d];
      fc->up[d] = up[d];
   }
   fc->orientation.fault_attrs = {3};
   DomainConfig dc;
   dc.fault_frames.config = fc;
   LinearElastic model(kLambda, kMu);
   Op op(pm, 1, model, 0.0, kL, 2.0 * kL, bc, DGMethod::IP, SolverType::MUMPS, dc);
   auto F = op.FrameField();
   MFEM_VERIFY(F, "L5: no frame field");

   // fault slip u+ - u- = J (the plate jump, always with the correct sign),
   // in the field's (d, s): S_dip = J . d, S_strike = J . s
   const real_t J[3] = {(mr < 0.0 ? 1.0 : -1.0) * Vp * t * sv[0],
                        (mr < 0.0 ? 1.0 : -1.0) * Vp * t * sv[1],
                        (mr < 0.0 ? 1.0 : -1.0) * Vp * t * sv[2]};
   const int nd = op.GetNumFaultDOFs();
   const int nf = F->NumFaultFaces();
   MFEM_VERIFY(nf == 0 || nd % nf == 0, "L5: DOFs per face");
   const int nbf = nf ? nd / nf : 1;
   Vector slip(2 * nd);
   for (int k = 0; k < nd; k++)
   {
      const fault_frame::Frame &fr = F->Frame(k / nbf);
      slip(2 * k) = J[0] * fr.d[0] + J[1] * fr.d[1] + J[2] * fr.d[2];
      slip(2 * k + 1) = J[0] * fr.s[0] + J[1] * fr.s[1] + J[2] * fr.s[2];
   }
   ParGridFunction u(&op.GetFESpace());
   op.Solve(t, slip, u);

   real_t err = 0.0;
   Vector val(3), xc(3);
   for (int e = 0; e < pm.GetNE(); e++)
   {
      ElementTransformation *T = pm.GetElementTransformation(e);
      const IntegrationPoint &ic = Geometries.GetCenter(pm.GetElementBaseGeometry(e));
      T->SetIntPoint(&ic);
      T->Transform(ic, xc);
      u.GetVectorValue(*T, ic, val);
      const real_t dd = xc(0) * mv[0] + xc(1) * mv[1] + xc(2) * mv[2];
      const real_t side = (dd > 0.0) ? 1.0 : -1.0;
      for (int d = 0; d < 3; d++) { err = std::max(err, std::abs(val(d) - side * 0.5 * Vp * t * sv[d])); }
   }
   return RMax(err, comm) / (0.5 * Vp * t);
}

static void TestOperatorSign(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== L5: band value against the operator's prescribed-slip sign ===\n"; }
   struct Case { const char *name; Rot R; bool rot_up; };
   const Case cases[] = {{"flat", Rot(), false},
                         {"dip60", Rot::Axis(1, 0, 0, 30.0), false},
                         {"general", Rot::Axis(1, 2, 3, 40.0), true}};
   for (const Case &c : cases)
   {
      const real_t up0[3] = {0, 0, 1};
      real_t up[3] = {0, 0, 1};
      if (c.rot_up) { c.R.Apply(up0, up); }
      const real_t e_ok = RigidError(comm, c.R, up, 1.0);
      const real_t e_bad = RigidError(comm, c.R, up, -1.0);
      TEST_ASSERT(e_ok <= 1e-9, "L5 " << c.name << ": rigid two-block motion (" << e_ok
                  << " of Vp t / 2)");
      TEST_ASSERT(e_bad > 0.1, "L5 " << c.name << ": with band_sign reversed the solution is not "
                  "rigid (" << e_bad << " of Vp t / 2)");
   }
}

// ---------------------------------------------------------------- L6, refusals
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
   const std::string path = std::string(tmp ? tmp : "/tmp") + "/seas_test_plate_loading_" + tag
                            + "_" + std::to_string(getpid()) + "_" + std::to_string(g_rank) + ".toml";
   std::ofstream f(path);
   f << text;
   MFEM_VERIFY(f, "cannot write " << path);
   return path;
}

static std::string WithLoading(const std::string &table)
{
   return ReadFile(kDeck) + "\n[cycle.loading]\n" + table;
}

static void TestDeck()
{
   if (g_rank == 0) { std::cout << "\n=== L6: [cycle.loading] in the deck ===\n"; }
   const std::string p1 = WriteDeck("bp5vec", WithLoading(
      "kind = \"plate\"\nplane_point = [0.0, 0.0, 0.0]\nplane_normal = [0.0, 1.0, 0.0]\n"
      "slip_dir = [1.0, 0.0, 0.0]\nband = 1000.0\n"));
   const cycle::CycleConfig c1 = cycle::CycleConfig::Load(p1);
   std::remove(p1.c_str());
   TEST_ASSERT(c1.loading.Plate() && c1.loading.band_sign == 1.0,
               "L6 BP5's vectors: kind plate, band_sign +1 (" << c1.loading.band_sign << ")");

   const std::string p2 = WriteDeck("defnormal", WithLoading(
      "kind = \"plate\"\nplane_point = [0.0, 0.0, 0.0]\nslip_dir = [-1.0, 0.0, 0.0]\nband = 1000.0\n"));
   const cycle::CycleConfig c2 = cycle::CycleConfig::Load(p2);
   std::remove(p2.c_str());
   const bool m_ref = c2.loading.plane_normal[0] == 0.0 && c2.loading.plane_normal[1] == -1.0
                      && c2.loading.plane_normal[2] == 0.0;
   TEST_ASSERT(m_ref && c2.loading.band_sign == -1.0,
               "L6 plane_normal left out: m = ref_normal (0, -1, 0), band_sign -1 ("
               << c2.loading.band_sign << ")");
   // the same plate motion as BP5's: m and s both reversed
   const DirichletFunc f2 = cycle::MakeCycleLoading(c2), f0 = MakeBP5DirichletFunc(c2.plate_rate);
   long bad = 0, n = 0;
   Vector x(3), a, b;
   for (real_t y : {-5000.0, -1000.0, -10.0, 0.0, 10.0, 1000.0, 5000.0})
   {
      x(0) = 1234.0; x(1) = y; x(2) = -4321.0;
      for (real_t t : {1.0, 1800.0 * kYear})
      {
         f2(x, t, a);
         f0(x, t, b);
         n++;
         if (!SameBits(a, b)) { bad++; }
      }
   }
   TEST_ASSERT(bad == 0, "L6 that deck loads as the BP5 function bit for bit (" << n
               << " point-times, " << bad << " differ)");
   const cycle::CycleConfig c0 = cycle::CycleConfig::Load(kDeck);
   TEST_ASSERT(!c0.loading.Plate(), "L6 no [cycle.loading]: kind \"bp5\"");
}

static void TestPoint2D()
{
   if (g_rank == 0) { std::cout << "\n=== L7: the plate function on a 2-D point ===\n"; }
   cycle::LoadingConfig L;
   L.kind = "plate";
   L.plane_point = {{0.0, 0.0, 0.0}};
   L.plane_normal = {{0.8660254037844386, 0.5, 0.0}};
   L.slip_dir = {{-0.5, 0.8660254037844386, 0.0}};
   L.band = 1000.0;
   L.band_sign = -1.0;   // sgn(-m . ref), ref = (1, 0, 0)
   const real_t Vp = 1.0e-9;
   const DirichletFunc f = cycle::MakePlateLoading(L, Vp);
   long n = 0, bad = 0, sizes = 0;
   Vector x2(2), x3(3), u2, u3;
   // along the plane normal: beyond, on and inside the band, either side
   for (real_t d : {-5.0e5, -1000.5, -1000.0, -10.0, 0.0, 10.0, 1000.0, 1000.5, 5.0e5})
   {
      for (real_t along : {-3.0e4, 0.0, 4.0e4})
      {
         x2(0) = d * L.plane_normal[0] - along * L.slip_dir[0];
         x2(1) = d * L.plane_normal[1] - along * L.slip_dir[1];
         x3(0) = x2(0); x3(1) = x2(1); x3(2) = 0.0;
         for (real_t t : {0.0, 1.0, 1500.0 * kYear})
         {
            f(x2, t, u2);
            f(x3, t, u3);
            n++;
            if (u2.Size() != 2 || u3.Size() != 3) { sizes++; continue; }
            if (std::memcmp(u2.GetData(), u3.GetData(), 2 * sizeof(real_t)) != 0 || u3(2) != 0.0)
            {
               bad++;
            }
         }
      }
   }
   TEST_ASSERT(sizes == 0 && bad == 0, "L7 (x, y): 2 entries, bit for bit the first two of "
               "(x, y, 0) (" << n << " point-times, " << sizes << " wrong sizes, " << bad
               << " differ)");
}

static int Refuse(int n)
{
   std::string text;
   switch (n)
   {
      case 1: text = WithLoading("kind = \"bp5\"\nband = 1000.0\n"); break;
      case 2: text = WithLoading("kind = \"plate\"\nplane_point = [0.0, 0.0, 0.0]\n"
                                 "slip_dir = [2.0, 0.0, 0.0]\nband = 1000.0\n"); break;
      case 3: text = WithLoading("kind = \"plate\"\nplane_point = [0.0, 0.0, 0.0]\n"
                                 "slip_dir = [0.0, 1.0, 0.0]\nband = 1000.0\n"); break;
      case 4: text = WithLoading("kind = \"plate\"\nplane_point = [0.0, 0.0, 0.0]\n"
                                 "plane_normal = [1.0, 0.0, 0.0]\nslip_dir = [0.0, 1.0, 0.0]\n"
                                 "band = 1000.0\n"); break;
      case 5:
      {
         text = ReadFile(kDeck);
         const std::string a = "ref_normal      = [0.0, -1.0, 0.0]";
         const size_t at = text.find(a);
         MFEM_VERIFY(at != std::string::npos && text.find(a, at + 1) == std::string::npos,
                     "--refuse 5: the deck's ref_normal line moved");
         text.replace(at, a.size(), "ref_normal      = [0.0, 1.0, 0.0]");
         break;
      }
      case 6: text = WithLoading("kind = \"slab\"\n"); break;
      default: MFEM_ABORT("--refuse N: N = 1..6");
   }
   const std::string p = WriteDeck("refuse" + std::to_string(n), text);
   if (g_rank == 0) { std::cout << "--refuse " << n << ": loading " << p << "\n"; }
   cycle::CycleConfig c = cycle::CycleConfig::Load(p);   // must abort
   std::remove(p.c_str());
   if (g_rank == 0) { std::cout << "--refuse " << n << ": the deck was accepted\n"; }
   return 0;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   g_rank = Mpi::WorldRank();
   int refuse = 0;
   for (int i = 1; i < argc; i++)
   {
      if (std::strcmp(argv[i], "--refuse") == 0 && i + 1 < argc) { refuse = std::atoi(argv[++i]); }
   }
   if (refuse) { return Refuse(refuse); }
   if (g_rank == 0)
   {
      std::cout << "test_plate_loading (2D in-plane design rev 3, phase A, A4) on "
                << Mpi::WorldSize() << " rank(s)\n";
   }
   TestFunction();
   TestOperatorSign(comm);
   TestDeck();
   TestPoint2D();
   if (g_rank == 0)
   {
      std::cout << "\ntest_plate_loading: " << (num_tests - num_failed) << " / " << num_tests
                << " passed" << (num_failed ? " -- FAILED" : "") << "\n";
   }
   return num_failed ? 1 : 0;
}
