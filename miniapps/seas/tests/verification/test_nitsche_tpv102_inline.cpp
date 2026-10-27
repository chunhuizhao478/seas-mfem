// Nitsche dynamic rupture, phase P6 (build sequence rev 5.2): first rupture,
// TPV102 on the inline mesh.  Plan item N-10, the fixture: the Nitsche
// counterpart of tests/verification/test_tpv102_local.cpp.
//
// Mesh: the 60 x 60 x 30 km box of test_tpv102_local.cpp (CreateTestMesh),
// split into 24 tetrahedra per hexahedron, with the fault on Y = 0 tagged
// as in P0 (attr 3), the free surface Z = 0 (attr 1) and absorbing walls
// (attr 5).  Material, prestress and nucleation are set in code from the
// constants the deck would carry (config/tpv102_params.hpp verbatim; plan
// N-07: global tensors, sigma_XY = +75 MPa, sigma_YY = -120 MPa).
//
// Modes:
//   (default)          M1: C-04, T-04 at p = 1, 2, 3 (uniform and x-graded
//                      box), T-21(a), shared fault-face count.
//   --t04-curved       T-04 negative control on a curved order-2 tet mesh:
//                      GetFaultFacePenalties spread and the would-be
//                      off-diagonal Jacobian norm (report, no abort).
//   --t04-curved-n01   the same mesh given to N-01: must abort (spread).
//   --abort-fac0 | --abort-fac1 | --abort-twice | --abort-linacc
//                      T-21(b), T-21(d): each must abort.
//   --m2 nx ny nz p [tfinal] [prefix | -] [dt_scale]
//                      M2: L1, L2, T-13(a); optional per-rank CSV of the
//                      fault state every step (full-run rank comparison).
//   --onestep prefix [p nx ny nz]
//                      one step from a smooth nontrivial state at t0 = 0.8 s,
//                      per-rank CSV (one-step rank comparison).
//
// Scope (Chunhui Zhao, 2026-10-14).

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../domain/nitsche_energy.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../friction/rate_state_closure.hpp"
#include "../../config/tpv102_params.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <set>
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

static const real_t kL = 30.0e3;   // half width; the box is [-L,L]x[-L,L]x[-L,0]
static const real_t kCuseM1[5] = {0.0, 0.057, 0.11, 0.14, 0.166};   // 0.5 C_eig (P3)

static int ISum(int v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static real_t RMin(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MIN, c); return v; }

// ---------------------------------------------------------------- fixture

/// The TPV102 box, 24 tets per hex; the fault (attr 3) on Y = 0 reaching the
/// free surface.  graded: y -> y (1 + 0.3 cos(pi x / 2L)) and
/// x -> x + 0.15 L sin(pi x / L): Y = 0 stays a plane, the Y -> -Y mirror
/// symmetry is kept, and the fault-normal element height varies along
/// strike, so sigma_F differs from face to face.
static Mesh *TPV102Box(int nx, int ny, int nz, bool graded = false)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2 * kL, 2 * kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
      if (graded)
      {
         // the fault-normal element height varies along strike (sigma_F
         // ~ nl / detJ does too); odd in Y, so the mirror symmetry stays
         c[1] *= 1.0 + 0.3 * std::cos(M_PI * c[0] / (2.0 * kL));
         c[0] += 0.15 * kL * std::sin(M_PI * c[0] / kL);
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

/// Partition `m` over MPI_COMM_WORLD.  Decision 1d (Chunhui Zhao,
/// 2026-10-14): the cross-rank fault path must run, so if METIS's default
/// partition leaves no fault face with its two elements on different ranks,
/// a forced partition is used instead: elements with centroid Y < 0 go to
/// ranks [0, np/2), Y > 0 to [np/2, np), each half split into np/2 slabs
/// along X.  Every fault face is then shared.  The choice is printed.
static std::string g_partition = "metis";
/// P7 regression (Chunhui Zhao, 2026-10-14): when set, the LAST rank gets
/// every element that has no face on the domain boundary (so it owns no
/// absorbing or free-surface face) and the others split the rest in X
/// slabs.  NitscheWaveOperator::SetTimeStep used to build the collective
/// damped mass only on ranks with absorbing faces and deadlocked here.
static bool g_interior_rank = false;
static std::unique_ptr<ParMesh> MakeParMesh(Mesh &m)
{
   const int np = Mpi::WorldSize();
   if (g_interior_rank && np >= 2)
   {
      std::vector<bool> on_bdr(m.GetNE(), false);
      for (int b = 0; b < m.GetNBE(); b++)
      {
         int e1, e2;
         m.GetFaceElements(m.GetBdrElementFaceIndex(b), &e1, &e2);
         if (e2 < 0 && m.GetBdrAttribute(b) != 3) { on_bdr[e1] = true; }
      }
      int *part = new int[m.GetNE()];
      int ninner = 0;
      for (int e = 0; e < m.GetNE(); e++)
      {
         if (!on_bdr[e]) { part[e] = np - 1; ninner++; continue; }
         Vector c(3);
         m.GetElementCenter(e, c);
         int slab = static_cast<int>(std::floor((c(0) + kL) / (2.0 * kL) * (np - 1)));
         part[e] = std::min(std::max(slab, 0), np - 2);
      }
      MFEM_VERIFY(ninner > 0, "interior-rank partition: the mesh has no interior element");
      g_partition = "interior-only last rank (" + std::to_string(ninner) + " elements)";
      auto pm = std::make_unique<ParMesh>(MPI_COMM_WORLD, m, part);
      delete [] part;
      return pm;
   }
   if (np == 1)
   {
      g_partition = "serial";
      return std::make_unique<ParMesh>(MPI_COMM_WORLD, m);
   }
   int *part = m.GeneratePartitioning(np);
   int cut = 0;
   for (int b = 0; b < m.GetNBE(); b++)
   {
      if (m.GetBdrAttribute(b) != 3) { continue; }
      int e1, e2;
      m.GetFaceElements(m.GetBdrElementFaceIndex(b), &e1, &e2);
      if (e2 >= 0 && part[e1] != part[e2]) { cut++; }
   }
   if (cut == 0)
   {
      g_partition = "forced Y-split";
      const int half = std::max(1, np / 2);
      for (int e = 0; e < m.GetNE(); e++)
      {
         Vector c(3);
         m.GetElementCenter(e, c);
         int slab = static_cast<int>(std::floor((c(0) + kL) / (2.0 * kL) * half));
         slab = std::min(std::max(slab, 0), half - 1);
         part[e] = ((c(1) < 0.0) ? 0 : half) + slab;
         if (part[e] >= np) { part[e] = np - 1; }
      }
   }
   else
   {
      g_partition = "metis (" + std::to_string(cut) + " fault faces cut)";
   }
   auto pm = std::make_unique<ParMesh>(MPI_COMM_WORLD, m, part);
   delete [] part;
   return pm;
}

// ------------------------------------------------------ graded box

/// Node coordinates of one axis on [a, b]: uniform cells of size h on the
/// fine interval [f0, f1] (f0 <= f1, both nodes), cells growing by the
/// factor g per cell outside it, each side then rescaled uniformly so that
/// it ends exactly on a (or b): the ratio between neighbours inside a side
/// stays g, and at the fine/coarse junction it is at most g.
static std::vector<real_t> AxisNodes(real_t a, real_t b, real_t f0, real_t f1,
                                     real_t h, real_t g)
{
   std::vector<real_t> fine;
   const int nf = static_cast<int>(std::lround((f1 - f0) / h));
   for (int i = 0; i <= nf; i++) { fine.push_back(f0 + (f1 - f0) * i / std::max(nf, 1)); }
   if (nf == 0) { fine.resize(1); }
   auto side = [h, g](real_t L) {
      std::vector<real_t> sz;
      real_t sum = 0.0, c = h;
      while (L - sum > 1e-9 * L)
      {
         c *= g;
         sz.push_back(c);
         sum += c;
      }
      for (real_t &v : sz) { v *= L / sum; }
      return sz;
   };
   std::vector<real_t> left = (f0 > a) ? side(f0 - a) : std::vector<real_t>();
   std::vector<real_t> right = (b > f1) ? side(b - f1) : std::vector<real_t>();
   std::vector<real_t> nodes, lnodes;
   real_t x = f0;
   for (real_t c : left) { x -= c; lnodes.push_back(x); }
   if (!lnodes.empty()) { lnodes.back() = a; }
   for (auto it = lnodes.rbegin(); it != lnodes.rend(); ++it) { nodes.push_back(*it); }
   for (real_t v : fine) { nodes.push_back(v); }
   x = f1;
   for (real_t c : right) { x += c; nodes.push_back(x); }
   if (!right.empty()) { nodes.back() = b; }
   return nodes;
}

/// The graded TPV102 box (decision, Chunhui Zhao, 2026-10-14: "Add the
/// graded box"): the same 60 x 60 x 30 km domain, MakeCartesian3D split into
/// 24 tets per hex on a logical grid, then a monotone piecewise-linear map of
/// the vertices onto per-axis node arrays.  Fine hexes of size h on
/// |x| <= x_fine and z in [z_bot, z_top] and next to Y = 0; cells grow by g
/// outside.  The 24-tet face and body centres sit at half-indices and map to
/// the midpoints of their cell, so every hex stays a rectangular box and
/// every tet stays affine.  The Y nodes are built for y > 0 and negated, so
/// the mesh is bitwise symmetric under Y -> -Y with Y = 0 exactly on the
/// fault nodes.  With 24 tets per hex the fault nodes are h/2 apart at p = 1
/// (hex corners and face centres).
struct GradedSpec
{
   real_t h = 1000.0, g = 1.2, x_fine = 6000.0, z_top = 0.0, z_bot = -12000.0;
};

static real_t MaxCell(const std::vector<real_t> &N)
{
   real_t m = 0.0;
   for (size_t i = 1; i < N.size(); i++) { m = std::max(m, N[i] - N[i - 1]); }
   return m;
}
static real_t MaxRatio(const std::vector<real_t> &N)
{
   real_t gm = 1.0;
   for (size_t i = 1; i + 1 < N.size(); i++)
   {
      const real_t a = N[i] - N[i - 1], b = N[i + 1] - N[i];
      gm = std::max(gm, std::max(a / b, b / a));
   }
   return gm;
}

static Mesh *TPV102GradedBox(const GradedSpec &gs)
{
   const std::vector<real_t> X = AxisNodes(-kL, kL, -gs.x_fine, gs.x_fine, gs.h, gs.g);
   const std::vector<real_t> Z = AxisNodes(-kL, 0.0, gs.z_bot, gs.z_top, gs.h, gs.g);
   // y >= 0: one fine cell at the fault, then growth (the first cell is h)
   std::vector<real_t> Yp = AxisNodes(0.0, kL, 0.0, gs.h, gs.h, gs.g);
   std::vector<real_t> Y;
   for (auto it = Yp.rbegin(); it != Yp.rend(); ++it) { Y.push_back(-*it); }
   for (size_t i = 1; i < Yp.size(); i++) { Y.push_back(Yp[i]); }
   Y[Yp.size() - 1] = 0.0;
   const int nx = static_cast<int>(X.size()) - 1, ny = static_cast<int>(Y.size()) - 1,
             nz = static_cast<int>(Z.size()) - 1;
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(nx, ny, nz, nx, ny, nz));
   auto map = [](const std::vector<real_t> &N, real_t s) {
      const long k = std::lround(2.0 * s);     // half-index
      return (k % 2 == 0) ? N[k / 2] : 0.5 * (N[(k - 1) / 2] + N[(k + 1) / 2]);
   };
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] = map(X, c[0]);
      c[1] = map(Y, c[1]);
      c[2] = map(Z, c[2]);
   }
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      bool top = true;
      for (int i = 0; i < verts.Size(); i++) { if (mesh->GetVertex(verts[i])[2] != 0.0) { top = false; } }
      mesh->SetBdrAttribute(b, top ? 1 : 5);
   }
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool on_plane = true;
      for (int i = 0; i < fv.Size(); i++) { if (mesh->GetVertex(fv[i])[1] != 0.0) { on_plane = false; } }
      if (on_plane) { tris.push_back({fv[0], fv[1], fv[2]}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   if (g_rank == 0)
   {
      std::cout << "  graded box: " << nx << " x " << ny << " x " << nz << " = "
                << nx * ny * nz << " hexes; fine hex " << gs.h << " m on |x| <= "
                << gs.x_fine << ", z in [" << gs.z_bot << ", " << gs.z_top
                << "] and next to Y = 0; max neighbour ratio x/y/z = " << MaxRatio(X)
                << "/" << MaxRatio(Y) << "/" << MaxRatio(Z) << "; largest cell x/y/z = "
                << MaxCell(X) << "/" << MaxCell(Y) << "/" << MaxCell(Z) << " m\n";
   }
   return mesh;
}

/// Tet quality: min/max inscribed-diameter length (the CFL length) and the
/// worst aspect ratio (longest edge / inscribed diameter).
static void MeshQuality(ParMesh &pm, real_t &hmin, real_t &hmax, real_t &aspect)
{
   hmin = 1e300; hmax = 0.0; aspect = 0.0;
   for (int e = 0; e < pm.GetNE(); e++)
   {
      const real_t he = NitscheElementLength(pm, e);
      Array<int> v;
      pm.GetElementVertices(e, v);
      real_t lmax = 0.0;
      for (int a = 0; a < 4; a++)
         for (int b = a + 1; b < 4; b++)
         {
            const real_t *pa = pm.GetVertex(v[a]), *pb = pm.GetVertex(v[b]);
            lmax = std::max(lmax, std::hypot(pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]));
         }
      hmin = std::min(hmin, he);
      hmax = std::max(hmax, he);
      aspect = std::max(aspect, lmax / he);
   }
   hmin = RMin(hmin, pm.GetComm());
   hmax = RMax(hmax, pm.GetComm());
   aspect = RMax(aspect, pm.GetComm());
}

static BoundaryConfig TPV102BC()
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};      // free surface
   bc.absorbing_attrs = {5};    // dashpot (P4)
   return bc;
}

static MaterialField TPV102Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu,
                                      TPV102Params::rho);
}

/// Friction law and closure of the deck (aging law, SCEC TPV102).
struct TPV102Friction
{
   DieterichRuinaFriction fr{{TPV102Params::V0, TPV102Params::f0, TPV102Params::b,
                              TPV102Params::Dc}};
   RateStateClosure cl{fr, RateStateClosure::StateLaw::Aging};
};

/// N-07 inputs as global tensors (code frame, tension positive).  Along
/// strike = X, down dip = -Z.
static NitscheFaultSetup TPV102Setup(const IFaultClosure &cl)
{
   NitscheFaultSetup fs;
   fs.closure = &cl;
   fs.prestress = [](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = TPV102Params::tau_ini;     // sigma_XY = sigma_YX = +75 MPa
      s[4] = -TPV102Params::sigma_n;           // sigma_YY = -120 MPa
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
   fs.sigma_min = 0.0;   // [pore_pressure] min_sigma_n_pa default (N-01)
   return fs;
}

/// State at t = 0 (N-10): d = 0, V = V_ini shat, psi = ComputeInitialPsi(a).
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

// ------------------------------------------------------------ C-04, T-04

/// Full finite-difference Jacobian of the projected traction (kernel
/// components dip, strike and sigma_n) with respect to every tangential slip
/// coefficient, frozen random u.  Returns max |J - sigma_F_bar delta| /
/// sigma_F_bar over all rows and columns on every rank; `offdiag` gets the
/// largest off-diagonal entry relative to sigma_F_bar.
static real_t FullJacobian(ElasticityDomainOperator<ParMesh> &E, const Vector &sigF,
                           real_t &offdiag, int &ncols_global, int &order_bad)
{
   MPI_Comm comm = E.GetMesh().GetComm();
   const int nd = E.GetNumFaultDOFs(), nbf = E.GetNbfPerFace();
   const int nranks = Mpi::WorldSize();
   ParGridFunction u(&E.GetFESpace());
   u.Randomize(1234 + g_rank);
   u -= 0.5;
   Vector slip(2 * nd);
   slip.Randomize(99 + g_rank);
   slip -= 0.5;
   const real_t h = 1.0;   // the kernel is affine in the slip
   Vector Tp, Tm, snp, snm;
   real_t worst = 0.0;
   offdiag = 0.0;
   order_bad = 0;
   ncols_global = 0;
   for (int r = 0; r < nranks; r++)
   {
      int ncols = (g_rank == r) ? 2 * nd : 0;
      MPI_Bcast(&ncols, 1, MPI_INT, r, comm);
      ncols_global += ncols;
      for (int col = 0; col < ncols; col++)
      {
         const bool mine = (g_rank == r);
         if (mine) { slip(col) += h; }
         E.ComputeTraction(u, slip, Tp, &snp);
         if (mine) { slip(col) -= 2 * h; }
         E.ComputeTraction(u, slip, Tm, &snm);
         if (mine) { slip(col) += h; }
         if (!mine)
         {
            for (int i = 0; i < Tp.Size(); i++) { if (Tp(i) != Tm(i)) { worst = 1e300; } }
            for (int i = 0; i < snp.Size(); i++) { if (snp(i) != snm(i)) { worst = 1e300; } }
            continue;
         }
         const int kc = col / 2, cc = col % 2, fk = kc / nbf;
         const real_t sb = sigF(fk);
         for (int i = 0; i < nd; i++)
         {
            for (int c = 0; c < 3; c++)
            {
               const real_t Jv = (c < 2) ? (Tp(2 * i + c) - Tm(2 * i + c)) / (2 * h)
                                         : (snp(i) - snm(i)) / (2 * h);
               const real_t expect = (i == kc && c == cc) ? sb : 0.0;
               const real_t e = std::abs(Jv - expect) / sb;
               worst = std::max(worst, e);
               if (expect == 0.0) { offdiag = std::max(offdiag, std::abs(Jv) / sb); }
               // order check: the diagonal must be THIS face's sigma_F and
               // not another face's
               if (i == kc && c == cc)
               {
                  for (int g = 0; g < sigF.Size(); g++)
                  {
                     if (g != fk && std::abs(sigF(g) - sb) > 1e-12 * sb &&
                         std::abs(Jv - sigF(g)) < std::abs(Jv - sb)) { order_bad++; }
                  }
               }
            }
         }
      }
   }
   offdiag = RMax(offdiag, comm);
   order_bad = ISum(order_bad, comm);
   return RMax(worst, comm);
}

static void TestC04T04(int nx, int ny, int nz, bool graded, int p, real_t &worst_out)
{
   std::unique_ptr<Mesh> m(TPV102Box(nx, ny, nz, graded));
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   MPI_Comm comm = pm.GetComm();
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
   ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
   const Vector &sigF = op.FacePenalty();
   const Vector &spread = op.FacePenaltySpread();
   real_t smax = 0.0, smin = 1e300, smx = 0.0;
   std::set<long long> distinct;
   for (int f = 0; f < sigF.Size(); f++)
   {
      smax = std::max(smax, spread(f));
      smin = std::min(smin, sigF(f));
      smx = std::max(smx, sigF(f));
      distinct.insert(std::llround(sigF(f) * 1e-3));   // kPa/m bins
   }
   smax = RMax(smax, comm);
   smin = RMin(smin, comm);
   smx = RMax(smx, comm);
   int ndist = ISum(static_cast<int>(distinct.size()), comm);
   // C-04 vs an independent evaluation of the Tandem formula at the face
   // centroid: interior penalty = (p0 + p1)/4, p_i = (D+1) c_N (D nl/detJ_i)
   // (c1^2/c0) with c_N = p (p + D - 1)/D.
   real_t worst_formula = 0.0;
   {
      const real_t lam = TPV102Params::lambda, mu = TPV102Params::mu, D = 3.0;
      const real_t cN = p * (p + D - 1.0) / D;
      const real_t ratio = (D * lam + 2 * mu) * (D * lam + 2 * mu) / (2 * mu);
      const Array<int> &fint = E.GetFaultInteriorFaces();
      const Array<int> &fsh = E.GetFaultSharedFaces();
      for (int fi = 0; fi < sigF.Size(); fi++)
      {
         FaceElementTransformations *Tr = (fi < fint.Size())
            ? pm.GetInteriorFaceTransformations(fint[fi])
            : pm.GetSharedFaceTransformations(fsh[fi - fint.Size()]);
         const IntegrationPoint &ip = Geometries.GetCenter(Tr->GetGeometryType());
         Tr->SetAllIntPoints(&ip);
         Vector nor(3);
         CalcOrtho(Tr->Jacobian(), nor);
         const real_t nl = nor.Norml2();
         const real_t p0 = (D + 1) * cN * (D * nl / Tr->Elem1->Weight()) * ratio;
         const real_t p1 = (D + 1) * cN * (D * nl / Tr->Elem2->Weight()) * ratio;
         const real_t ref = (p0 + p1) / 4.0;
         worst_formula = std::max(worst_formula, std::abs(sigF(fi) - ref) / ref);
      }
      worst_formula = RMax(worst_formula, comm);
   }
   real_t off = 0.0;
   int ncols = 0, order_bad = 0;
   const real_t worst = FullJacobian(E, sigF, off, ncols, order_bad);
   worst_out = worst;
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  " << (graded ? "x-graded" : "uniform")
                << " " << nx << "x" << ny << "x" << nz << " box, p = " << p
                << ": fault faces " << op.Orientation().GetReport().faces
                << ", shared fault faces " << op.NumSharedFaultFacesGlobal()
                << " (partition " << g_partition << ")"
                << ", sigma_F in [" << smin << ", " << smx << "] Pa/m ("
                << ndist << " distinct values summed over ranks), max spread "
                << smax << ", vs Tandem formula " << worst_formula
                << "\n     T-04: " << ncols << " slip columns, max |J - sigma_F I| / sigma_F = "
                << worst << ", max off-diagonal / sigma_F = " << off
                << ", diagonal closer to another face's sigma_F: " << order_bad << "\n";
   }
   TEST_ASSERT(smax <= NitscheWaveOperator::kPenaltySpreadTol,
               "C-04 spread <= 1e-10 on every fault face");
   TEST_ASSERT(worst_formula < 1e-13, "C-04 sigma_F == Tandem penalty formula (1e-13)");
   TEST_ASSERT(worst < 1e-6, "T-04 full Jacobian == sigma_F_bar I (1e-6), off-diagonals included");
   TEST_ASSERT(order_bad == 0, "T-04 face order of sigma_F == ComputeTraction's");
}

// ------------------------------------- fault-DOF scatter consistency

/// Open question 5: ExpandOwnedToLocalFault must put on every local fault
/// DOF (ghost copies of shared faces included) the value owned at the SAME
/// physical node.  Field f(x) = x + 2 z + 1e-3 x z on owned DOFs, expanded,
/// compared with f at the local DOF positions of GetFaultDOFCoords3D.
static void TestScatter(int p)
{
   std::unique_ptr<Mesh> m(TPV102Box(8, 8, 4));
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
   ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
   auto f = [](const real_t *x) { return x[0] + 2.0 * x[2] + 1e-3 * x[0] * x[2]; };
   const int no = op.NumOwnedFaultDOFs(), nd = E.GetNumFaultDOFs();
   Vector fo(no), floc, xloc;
   for (int j = 0; j < no; j++) { fo(j) = f(&op.FaultPositions()(3 * j)); }
   E.ExpandOwnedToLocalFault(fo, floc, 1);
   E.GetFaultDOFCoords3D(xloc);
   const int n_int = E.GetFaultInteriorFaces().Size(), nbf = E.GetNbfPerFace();
   int bad_int = 0, bad_sh = 0;
   real_t worst = 0.0;
   for (int k = 0; k < nd; k++)
   {
      const real_t e = std::abs(floc(k) - f(&xloc(3 * k)));
      worst = std::max(worst, e);
      if (e > 1e-6) { ((k / nbf < n_int) ? bad_int : bad_sh)++; }
   }
   MPI_Comm comm = pm.GetComm();
   bad_int = ISum(bad_int, comm); bad_sh = ISum(bad_sh, comm);
   worst = RMax(worst, comm);
   const int nsh = op.NumSharedFaultFacesGlobal();
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  p = " << p << ": shared fault faces " << nsh
                << " (partition " << g_partition << "); local DOFs whose expanded value is not "
                   "f(own position): interior faces " << bad_int << ", shared faces " << bad_sh
                << "; max |error| = " << worst << " m\n";
   }
   TEST_ASSERT(bad_int + bad_sh == 0, "fault-DOF scatter: every local copy holds its own node's value");
}

// ------------------------------------------ N-07 SCEC prestress check

/// T_0 = sigma_0 nhat = (-75, 120, 0) MPa at every owned DOF; in the frame
/// X_tan = -75 MPa shat (-T_0 . shat = 75 MPa), sigma_n_eff = 120 MPa, and
/// the closure at the seeded psi^0 gives V parallel to +shat.
static void TestPrestress(int p)
{
   std::unique_ptr<Mesh> m(TPV102Box(4, 4, 2));
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
   op.SetTimeStep(ComputeMaxDt(op.Sample(), pm, p, kCuseM1[p]));
   TPV102InitialState(op);
   Vector S;
   op.BackgroundTraction(0.0, S);
   const real_t T0[3] = {-TPV102Params::tau_ini, TPV102Params::sigma_n, 0.0};
   real_t worst = 0.0, worstX = 0.0, worst_sig = 0.0;
   int bad_dir = 0;
   for (int j = 0; j < op.NumOwnedFaultDOFs(); j++)
   {
      for (int c = 0; c < 3; c++) { worst = std::max(worst, std::abs(S(3 * j + c) - T0[c]) / 120e6); }
      real_t V[3], sig, X[2];
      op.CloseDOF(j, &S(3 * j), op.State()(j), V, sig, X);
      worstX = std::max({worstX, std::abs(X[0]) / 75e6, std::abs(X[1] + 75e6) / 75e6});
      worst_sig = std::max(worst_sig, std::abs(sig - 120e6) / 120e6);
      const real_t *sv = &op.FaultFrames()(9 * j + 6);
      if (!(fault_frame::Dot(V, sv) > 0.0)) { bad_dir++; }
   }
   MPI_Comm comm = pm.GetComm();
   worst = RMax(worst, comm); worstX = RMax(worstX, comm); worst_sig = RMax(worst_sig, comm);
   bad_dir = ISum(bad_dir, comm);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  p = " << p << ": max |S - (-75, 120, 0) MPa| / 120 MPa = "
                << worst << "; X_tan vs (0, -75 MPa): " << worstX << "; sigma_n_eff vs 120 MPa: "
                << worst_sig << "; DOFs with V.shat <= 0: " << bad_dir << "\n";
   }
   TEST_ASSERT(worst < 1e-14, "N-07 projected prestress == (-75, 120, 0) MPa at every DOF");
   TEST_ASSERT(worstX < 1e-14 && worst_sig < 1e-14,
               "N-07 X_tan = -75 MPa shat, sigma_n_eff = 120 MPa (SCEC values)");
   TEST_ASSERT(bad_dir == 0, "C-07 V^0 parallel to +shat (V_x > 0) at every DOF");
}

// -------------------------------------------------------------- T-21(a)

static void TestT21a(int p)
{
   std::unique_ptr<Mesh> m(TPV102Box(4, 4, 2));
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuseM1[p]);
   op.SetTimeStep(dt);
   TPV102InitialState(op);
   // a nonzero stored slip so that G d is exercised
   Vector d = op.Slip(), V = op.SlipRate(), psi = op.State();
   for (int j = 0; j < op.NumOwnedFaultDOFs(); j++)
   {
      const real_t *s = &op.FaultFrames()(9 * j + 6);
      for (int c = 0; c < 3; c++) { d(3 * j + c) = 0.01 * std::sin(0.37 * j) * s[c]; }
   }
   op.SetFaultState(d, V, psi);
   const int n = op.Height();
   Vector x(n), v(n), y(n);
   x.Randomize(5 + g_rank); x -= 0.5; x *= 1e-3;
   v.Randomize(6 + g_rank); v -= 0.5;
   op.SetTime(0.3);
   op.Mult(x, v, y);
   const bool unchanged = std::memcmp(op.Slip().GetData(), d.GetData(), sizeof(real_t) * d.Size()) == 0
                          && std::memcmp(op.SlipRate().GetData(), V.GetData(), sizeof(real_t) * V.Size()) == 0
                          && std::memcmp(op.State().GetData(), psi.GetData(), sizeof(real_t) * psi.Size()) == 0;
   // independent M^{-1}(F - K x + G d - C v) (here F = 0: no Dirichlet faces)
   Vector Kx, Gd, Cv, r(n), a;
   op.Elasticity().ApplyStiffness(x, Kx);
   op.SlipLoad(d, Gd);
   op.ApplyDamping(v, Cv);
   Vector Fd;
   op.Elasticity().AssembleDirichletOnlyRHS(Fd, 0.3);
   r = Fd; r -= Kx; r += Gd; r -= Cv;
   op.Mass().Apply(r, a);
   int diff = 0;
   for (int i = 0; i < n; i++) { if (a(i) != y(i)) { diff++; } }
   diff = ISum(diff, pm.GetComm());
   const int nbad = ISum(unchanged ? 0 : 1, pm.GetComm());
   real_t gdn = RMax(Gd.Normlinf(), pm.GetComm());
   if (g_rank == 0)
   {
      std::cout << "  T-21(a) p = " << p << ": |G d|_inf = " << gdn
                << ", entries of Mult differing from M^{-1}(F - Kx + Gd - Cv): " << diff << "\n";
   }
   TEST_ASSERT(nbad == 0, "T-21(a) Mult leaves (delta, V, psi) bitwise unchanged");
   TEST_ASSERT(diff == 0 && gdn > 0.0, "T-21(a) Mult == M^{-1}(F - K x + G d - C v) at the stored state, bit for bit");
}

// ------------------------------------------------------- abort modes

static int AbortMode(const std::string &mode)
{
   std::unique_ptr<Mesh> m(TPV102Box(4, 4, 2));
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl);
   NitscheWaveOperator op(pm, 1, TPV102Material(), TPV102BC(), {}, &fs);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, 1, 0.5 * 0.113);
   op.SetTimeStep(dt);
   TPV102InitialState(op);
   const int n = op.Height();
   Vector U(n), Ud(n), k(n);
   U = 0.0; Ud = 0.0;
   if (g_rank == 0) { std::cout << "  abort mode " << mode << ": expecting an abort\n" << std::flush; }
   if (mode == "--abort-fac0")
   {
      op.SetTime(dt);
      op.ImplicitSolve(dt * dt / 6.0, 0.5 * dt, U, Ud, k);
   }
   else if (mode == "--abort-fac1")
   {
      op.SetTime(dt);
      op.ImplicitSolve(0.0, 0.75 * dt, U, Ud, k);
   }
   else if (mode == "--abort-twice")
   {
      op.SetTime(dt);
      op.ImplicitSolve(0.0, 0.5 * dt, U, Ud, k);
      if (g_rank == 0) { std::cout << "  first call at t = dt returned\n" << std::flush; }
      op.ImplicitSolve(0.0, 0.5 * dt, U, Ud, k);
   }
   else if (mode == "--abort-linacc")
   {
      LinearAccelerationSolver ode;
      ode.Init(op);
      real_t t = 0.0, h = dt;
      ode.Step(U, Ud, t, h);
   }
   if (g_rank == 0) { std::cout << "  NO ABORT in " << mode << "\n"; }
   return 0;   // exit 0 = the expected abort did not happen (test failure)
}

// ----------------------------------------------- T-04 negative control

/// Curved order-2 tets: y -> y (1 + 0.2 sin(pi z / 7.5 km)) on the nodes.
/// The factor is 1 at every vertex (the 24-tet vertices sit on multiples of
/// 7.5 km in z) and the map keeps Y = 0 fixed, so the fault stays planar
/// with its vertices on Y = 0, while the fault-normal element height --
/// and with it detJ, hence sigma_F ~ nl / detJ -- varies over each face.
/// (An in-plane map would not do: it scales nl and detJ alike.)
static Mesh *CurvedBox()
{
   Mesh *m = TPV102Box(4, 4, 2);
   m->SetCurvature(2);
   GridFunction *nodes = m->GetNodes();
   const FiniteElementSpace *nfes = nodes->FESpace();
   const int ndofs = nfes->GetNDofs();
   for (int i = 0; i < ndofs; i++)
   {
      const int iy = nfes->DofToVDof(i, 1), iz = nfes->DofToVDof(i, 2);
      (*nodes)(iy) *= 1.0 + 0.2 * std::sin(M_PI * (*nodes)(iz) / 7.5e3);
   }
   return m;
}

static int CurvedReport()
{
   std::unique_ptr<Mesh> m(CurvedBox());
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   MPI_Comm comm = pm.GetComm();
   // the operator N-01 would build, without N-01
   MaterialField mat = TPV102Material();
   NitscheMaterialSample sample(pm, mat);
   ElasticityDomainOperator<ParMesh> E(pm, 1, sample.Lambda(), sample.Mu(),
                                       TPV102Params::lambda, TPV102Params::mu,
                                       0.0, 0.0, 0.0, TPV102BC(), DGMethod::IP,
                                       SolverType::CG_AMG, DomainConfig{});
   Vector sigF, spread;
   E.GetFaultFacePenalties(sigF, spread);
   real_t smax = 0.0;
   for (int f = 0; f < spread.Size(); f++) { smax = std::max(smax, spread(f)); }
   smax = RMax(smax, comm);
   real_t off = 0.0;
   int ncols = 0, order_bad = 0;
   const real_t worst = FullJacobian(E, sigF, off, ncols, order_bad);
   // N-08 on this mesh (vertices planar): does it refuse?
   FaultOrientation fo(pm, E.GetFaultInteriorFaces(), E.GetFaultSharedFaces());
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6)
                << "  curved order-2 tet box (p = 1): N-08 accepted the fault ("
                << fo.GetReport().faces << " faces, max d_F " << fo.GetReport().max_face_departure
                << " m); C-04 max spread = " << smax << " (tolerance "
                << NitscheWaveOperator::kPenaltySpreadTol << ")\n     would-be T-04: max |J - sigma_F_bar I| / sigma_F_bar = "
                << worst << ", max off-diagonal / sigma_F_bar = " << off << "\n";
   }
   TEST_ASSERT(smax > NitscheWaveOperator::kPenaltySpreadTol,
               "T-04 negative control: curved mesh spread above tolerance");
   TEST_ASSERT(worst > 1e-6, "T-04 negative control: Jacobian is not sigma_F I on the curved mesh");
   return (num_failed == 0) ? 0 : 1;
}

static int CurvedN01()
{
   std::unique_ptr<Mesh> m(CurvedBox());
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl);
   if (g_rank == 0) { std::cout << "  curved mesh given to N-01: expecting an abort\n" << std::flush; }
   NitscheWaveOperator op(pm, 1, TPV102Material(), TPV102BC(), {}, &fs);
   if (g_rank == 0) { std::cout << "  NO ABORT: N-01 accepted the curved mesh\n"; }
   return 0;   // exit 0 = the expected abort did not happen (test failure)
}

/// Centroid of the fault face each owned fault DOF belongs to (3 per DOF).
/// Partition-invariant; with the DOF position it keys a DG fault DOF
/// uniquely across rank counts (the same point sits on several faces).
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

/// One CSV row per owned fault DOF: step, t, x, y, z, face centroid (3),
/// d (3), V (3), psi, X_dip, X_strike, sigma_eff.
static void DumpFault(std::ostream &os, NitscheWaveOperator &op, const Vector &cf,
                      long step, real_t t)
{
   const Vector &X = op.FaultPositions(), &d = op.Slip(), &V = op.SlipRate();
   const Vector &psi = op.State(), &Xs = op.LastShear(), &sg = op.LastSigmaEff();
   for (int j = 0; j < op.NumOwnedFaultDOFs(); j++)
   {
      os << step << "," << t << "," << X(3 * j) << "," << X(3 * j + 1) << "," << X(3 * j + 2)
         << "," << cf(3 * j) << "," << cf(3 * j + 1) << "," << cf(3 * j + 2)
         << "," << d(3 * j) << "," << d(3 * j + 1) << "," << d(3 * j + 2)
         << "," << V(3 * j) << "," << V(3 * j + 1) << "," << V(3 * j + 2)
         << "," << psi(j) << "," << Xs(2 * j) << "," << Xs(2 * j + 1) << "," << sg(j) << "\n";
   }
}

// ------------------------------------------------- rupture run (M2, M3)

static const real_t kCuse[5] = {0.0, 0.057, 0.11, 0.14, 0.166};   // 0.5 C_eig (P3)

struct RunConfig
{
   int nx = 4, ny = 4, nz = 2, p = 2;
   bool graded = false;
   bool graded_box = false;        ///< use TPV102GradedBox(gspec)
   GradedSpec gspec;
   real_t tfinal = 2.0;
   real_t dt_scale = 1.0;          ///< dt = dt_scale * (0.5 C_eig CFL step)
   std::string dump;               ///< per-rank CSV prefix ("" = none)
   int dump_every = 1;
};

struct RunResult
{
   bool finite = true;
   long steps = 0;
   real_t dt = 0.0, t_end = 0.0;
   long long hexes = 0, tets = 0;
   int fault_faces = 0, shared_faces = 0, fault_dofs = 0;
   real_t Vmax = 0.0;                    ///< max over run and DOFs
   real_t hypo_dist = 0.0;               ///< distance of the hypocentre DOF
   real_t hypo_Vpeak = 0.0, hypo_tinit = -1.0;
   real_t first_t = -1.0, first_dist = -1.0;   ///< first DOF with |V| >= 1e-3
   real_t x_min = 0.0, x_max = 0.0;      ///< along-strike extent of ruptured DOFs
   real_t sig_min = 1e300, sig_max = -1e300;
   real_t max_ratio = 0.0;
   real_t hmin = 0.0, hmax = 0.0, aspect = 0.0;
   long long vdofs = 0;
};

static void HypoXYZ(real_t h[3])
{
   h[0] = TPV102Params::hypo_along_strike;
   h[1] = 0.0;
   h[2] = -TPV102Params::hypo_down_dip;
}

static RunResult RunRupture(const RunConfig &rc)
{
   RunResult R;
   std::unique_ptr<Mesh> m(rc.graded_box ? TPV102GradedBox(rc.gspec)
                                         : TPV102Box(rc.nx, rc.ny, rc.nz, rc.graded));
   R.tets = m->GetNE();
   R.hexes = R.tets / 24;
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   m.reset();
   MeshQuality(pm, R.hmin, R.hmax, R.aspect);
   MPI_Comm comm = pm.GetComm();
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl);
   NitscheWaveOperator op(pm, rc.p, TPV102Material(), TPV102BC(), {}, &fs);
   const real_t dt = rc.dt_scale * ComputeMaxDt(op.Sample(), pm, rc.p, kCuse[rc.p]);
   op.SetTimeStep(dt);
   TPV102InitialState(op);
   R.dt = dt;
   R.fault_faces = op.Orientation().GetReport().faces;
   R.shared_faces = op.NumSharedFaultFacesGlobal();
   const int no = op.NumOwnedFaultDOFs();
   R.fault_dofs = ISum(no, comm);
   {
      long long nloc = op.Height();
      MPI_Allreduce(MPI_IN_PLACE, &nloc, 1, MPI_LONG_LONG, MPI_SUM, comm);
      R.vdofs = nloc;
   }

   // hypocentre DOF: nearest owned DOF to (0, 0, -7.5 km), global MINLOC
   real_t h[3];
   HypoXYZ(h);
   struct { double d; int r; } loc = {1e300, g_rank}, glob;
   int jh = -1;
   const Vector &X = op.FaultPositions();
   for (int j = 0; j < no; j++)
   {
      const real_t dd = std::hypot(X(3 * j) - h[0], X(3 * j + 1) - h[1], X(3 * j + 2) - h[2]);
      if (dd < loc.d) { loc.d = dd; jh = j; }
   }
   MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
   R.hypo_dist = glob.d;
   const bool hypo_mine = (glob.r == g_rank);

   const int n = op.Height();
   Vector U(n), Ud(n);
   U = 0.0;
   Ud = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(op);
   real_t t = 0.0;
   const long nsteps = std::lround(std::ceil(rc.tfinal / dt - 1e-9));
   std::vector<real_t> trupt(no, -1.0);
   Vector cf;
   OwnedFaceCentroids(op, cf);
   std::unique_ptr<std::ofstream> os;
   if (!rc.dump.empty())
   {
      os = std::make_unique<std::ofstream>(rc.dump + "_r" + std::to_string(g_rank) + ".csv");
      *os << std::setprecision(17);
   }
   for (long s = 0; s < nsteps; s++)
   {
      real_t dt_s = dt;
      ode.Step(U, Ud, t, dt_s);
      const NitscheStepRecord &rec = op.LastStep();
      R.sig_min = std::min(R.sig_min, rec.min_sigma_eff);
      R.sig_max = std::max(R.sig_max, rec.max_sigma_eff);
      R.Vmax = std::max(R.Vmax, rec.max_V);
      R.max_ratio = std::max(R.max_ratio, rec.max_V_dt_over_Dc);
      const Vector &V = op.SlipRate();
      for (int j = 0; j < no; j++)
      {
         const real_t va = std::sqrt(V(3 * j) * V(3 * j) + V(3 * j + 1) * V(3 * j + 1)
                                     + V(3 * j + 2) * V(3 * j + 2));
         if (trupt[j] < 0.0 && va >= 1e-3) { trupt[j] = t; }
         if (hypo_mine && j == jh) { R.hypo_Vpeak = std::max(R.hypo_Vpeak, va); }
      }
      if (os && (s % rc.dump_every == 0 || s == nsteps - 1))
      {
         DumpFault(*os, op, cf, s + 1, t);
      }
      if (!std::isfinite(rec.max_V) || !std::isfinite(rec.min_sigma_eff))
      {
         R.finite = false;
         break;
      }
      if (g_rank == 0 && (s % 50 == 0 || s == nsteps - 1))
      {
         std::cout << std::setprecision(6) << "    step " << rec.step << " t = " << rec.t
                   << "  max|V| = " << rec.max_V << "  sigma_eff in [" << rec.min_sigma_eff
                   << ", " << rec.max_sigma_eff << "]  max|V|dt/Dc = "
                   << rec.max_V_dt_over_Dc << "\n";
      }
   }
   R.steps = op.LastStep().step;
   R.t_end = t;
   R.finite = (ISum(R.finite ? 0 : 1, comm) == 0) && std::isfinite(RMax(U.Normlinf(), comm));
   // hypocentre history, first rupture, along-strike extent
   real_t hv[2] = {R.hypo_Vpeak, -1.0};
   if (hypo_mine) { hv[1] = trupt[jh]; }
   MPI_Bcast(hv, 2, MPI_DOUBLE, glob.r, comm);
   R.hypo_Vpeak = hv[0];
   R.hypo_tinit = hv[1];
   struct { double t; int r; } fl = {1e300, g_rank}, fg;
   int jf = -1;
   real_t xmin = 1e300, xmax = -1e300;
   for (int j = 0; j < no; j++)
   {
      if (trupt[j] < 0.0) { continue; }
      if (trupt[j] < fl.t) { fl.t = trupt[j]; jf = j; }
      xmin = std::min(xmin, X(3 * j));
      xmax = std::max(xmax, X(3 * j));
   }
   MPI_Allreduce(&fl, &fg, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
   real_t fd = -1.0;
   if (fg.r == g_rank && jf >= 0)
   {
      fd = std::hypot(X(3 * jf) - h[0], X(3 * jf + 1) - h[1], X(3 * jf + 2) - h[2]);
   }
   MPI_Bcast(&fd, 1, MPI_DOUBLE, fg.r, comm);
   R.first_t = (fg.t < 1e299) ? fg.t : -1.0;
   R.first_dist = fd;
   R.x_min = RMin(xmin, comm);
   R.x_max = RMax(xmax, comm);
   return R;
}

/// M2: L1, L2 and T-13(a) on one mesh at one order.
static void TestM2(const RunConfig &rc)
{
   if (g_rank == 0)
   {
      if (rc.graded_box)
      {
         std::cout << "\n=== M2 rupture: graded box (h_fine " << rc.gspec.h << " m), p = "
                   << rc.p << ", " << rc.tfinal << " s ===\n";
      }
      else
      {
         std::cout << "\n=== M2 rupture: " << rc.nx << "x" << rc.ny << "x" << rc.nz
                   << " box, p = " << rc.p << ", " << rc.tfinal << " s ===\n";
      }
   }
   const RunResult R = RunRupture(rc);
   const real_t sn = TPV102Params::sigma_n;
   const real_t dev = std::max(std::abs(R.sig_min - sn), std::abs(R.sig_max - sn));
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6)
                << "  mesh: " << R.hexes << " hexes, " << R.tets << " tets (CFL length "
                << R.hmin << " - " << R.hmax << " m, worst aspect " << R.aspect << "), "
                << R.vdofs << " volume DOFs, "
 << R.fault_faces << " fault faces, " << R.fault_dofs << " fault DOFs; "
                << "shared fault faces " << R.shared_faces << " (partition "
                << g_partition << "); dt = " << R.dt
                << " s, " << R.steps << " steps to t = " << R.t_end << "\n"
                << "  hypocentre DOF at distance " << R.hypo_dist << " m: |V| >= 1e-3 at t = "
                << R.hypo_tinit << " s, peak |V| = " << R.hypo_Vpeak << " m/s\n"
                << "  first DOF with |V| >= 1e-3: t = " << R.first_t << " s, "
                << R.first_dist << " m from the hypocentre\n"
                << "  ruptured along strike x in [" << R.x_min << ", " << R.x_max
                << "] m; max |V| over run = " << R.Vmax << " m/s; max |V|dt/Dc = "
                << R.max_ratio << "\n"
                << "  sigma_eff in [" << R.sig_min << ", " << R.sig_max
                << "] Pa; max |sigma_eff - 120 MPa| / 120 MPa = " << dev / sn << "\n";
   }
   const real_t R_nuc = TPV102Params::nuc_radius;
   TEST_ASSERT(R.finite, "L2 finite and stable for the whole run");
   TEST_ASSERT(R.hypo_tinit >= 0.0 && R.hypo_tinit <= TPV102Params::nuc_T,
               "L1 rupture initiates at the hypocentre DOF within the nucleation time (1 s)");
   TEST_ASSERT(R.first_dist >= 0.0 && R.first_dist <= R_nuc + R.hypo_dist,
               "L1 the first DOF to rupture lies in the nucleation patch");
   TEST_ASSERT(R.hypo_Vpeak > 0.1, "L1 peak slip rate at the hypocentre > 0.1 m/s");
   TEST_ASSERT(R.x_min < -R_nuc && R.x_max > R_nuc,
               "L2 bilateral: ruptured DOFs beyond the nucleation radius on both sides");
   TEST_ASSERT(dev < 1e-8 * sn, "T-13(a) |sigma_n_eff - 120 MPa| < 1e-8 x 120 MPa at all steps and DOFs");
}

/// Smooth, partition-independent nontrivial state for the one-step rank
/// comparison (t0 = 0.8 s, nucleation active).
static void SmoothU(const Vector &x, Vector &u)
{
   u.SetSize(3);
   u(0) = 1e-2 * std::sin(x(0) / 7e3) * std::cos(x(2) / 5e3) * std::exp(-x(1) * x(1) / 4e8);
   u(1) = 4e-3 * std::cos(x(0) / 9e3) * std::sin(x(1) / 6e3);
   u(2) = -6e-3 * std::sin(x(2) / 8e3) * std::cos(x(1) / 11e3);
}
static void SmoothUd(const Vector &x, Vector &u)
{
   u.SetSize(3);
   u(0) = 0.2 * std::cos(x(0) / 6e3) * std::exp(-x(1) * x(1) / 1e8);
   u(1) = 0.05 * std::sin(x(2) / 7e3);
   u(2) = -0.1 * std::sin(x(0) / 5e3) * std::cos(x(1) / 9e3);
}

static void OneStep(const std::string &prefix, int p, int nx, int ny, int nz)
{
   std::unique_ptr<Mesh> m(TPV102Box(nx, ny, nz));
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   m.reset();
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   op.SetTimeStep(dt);
   const int no = op.NumOwnedFaultDOFs();
   const Vector &X = op.FaultPositions();
   real_t h[3];
   HypoXYZ(h);
   Vector d(3 * no), V(3 * no), psi(no);
   for (int j = 0; j < no; j++)
   {
      const real_t r2 = std::pow(X(3 * j) - h[0], 2) + std::pow(X(3 * j + 2) - h[2], 2);
      const real_t *s = &op.FaultFrames()(9 * j + 6), *dd = &op.FaultFrames()(9 * j + 3);
      for (int c = 0; c < 3; c++)
      {
         // small seeded slip: with no matching displacement jump the penalty
         // turns d into a stress drop sigma_F d (~1 MPa per mm at p = 2)
         d(3 * j + c) = 5e-4 * std::exp(-r2 / 2.5e7) * s[c] + 1e-4 * std::exp(-r2 / 1e8) * dd[c];
         V(3 * j + c) = (TPV102Params::V_ini + 0.3 * std::exp(-r2 / 1.6e7)) * s[c];
      }
      // state lowered by up to 0.2 near the hypocentre: V^{n+1} ~ 0.1 m/s
      // there, so the closure runs in the sliding regime on every rank count
      psi(j) = ComputeInitialPsi(op.FaultA()(j)) - 0.2 * std::exp(-r2 / 3e7);
   }
   op.SetFaultState(d, V, psi);
   ParGridFunction gu(&op.Elasticity().GetFESpace()), gv(&op.Elasticity().GetFESpace());
   VectorFunctionCoefficient cu(3, SmoothU), cv(3, SmoothUd);
   gu.ProjectCoefficient(cu);
   gv.ProjectCoefficient(cv);
   Vector U(gu), Ud(gv);
   CentralDifferenceSolver ode;
   ode.Init(op);
   real_t t = 0.8, dts = dt;
   ode.Step(U, Ud, t, dts);
   std::ofstream os(prefix + "_r" + std::to_string(g_rank) + ".csv");
   os << std::setprecision(17);
   Vector cf;
   OwnedFaceCentroids(op, cf);
   DumpFault(os, op, cf, 1, t);
   if (g_rank == 0)
   {
      const NitscheStepRecord &rec = op.LastStep();
      std::cout << std::setprecision(6) << "  one step (p = " << p << ", dt = " << dt
                << ") from the smooth state at t0 = 0.8: shared fault faces "
                << op.NumSharedFaultFacesGlobal() << " (partition " << g_partition
                << "), max|V| = " << rec.max_V
                << ", sigma_eff in [" << rec.min_sigma_eff << ", " << rec.max_sigma_eff
                << "]; wrote " << prefix << "_r*.csv\n";
   }
}

// ------------------------------------------------------------- T-15

struct EnergyRun
{
   real_t dt = 0.0, sum_rho = 0.0, sum_abs = 0.0, E0 = 0.0, Emax = 0.0;
   real_t fric_min = 0.0, fric_max = 0.0, node_max = 0.0, node_scale = 0.0;
   real_t W_pre = 0.0, W_nuc = 0.0, W_fric = 0.0, W_abs = 0.0;   // time integrals
   NitscheEnergyRecord last;
   long steps = 0;
};

static EnergyRun RunEnergy(int nx, int ny, int nz, int p, real_t tfinal,
                           real_t dt_scale, bool drop_cons, const std::string &csv = "")
{
   EnergyRun R;
   std::unique_ptr<Mesh> m(TPV102Box(nx, ny, nz));
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   m.reset();
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
   const real_t dt = dt_scale * ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   R.dt = dt;
   op.SetTimeStep(dt);
   TPV102InitialState(op);
   NitscheEnergyMonitor mon(op, drop_cons);
   const int n = op.Height();
   Vector U(n), Ud(n), A0(n);
   U = 0.0; Ud = 0.0;
   op.SetTime(0.0);
   op.Mult(U, Ud, A0);
   NitscheEnergyRecord r = mon.Evaluate(U, Ud, A0, 0.0);
   mon.Start(r);
   R.E0 = r.E_h;
   std::unique_ptr<std::ofstream> os;
   if (!csv.empty() && g_rank == 0)
   {
      os = std::make_unique<std::ofstream>(csv);
      *os << std::setprecision(17)
          << "t,E_h,E_kin,E_el,E_sipg,E_cons,E_pen,E_mod,P_ext,P_abs,P_pre,P_nuc,P_fric,rho\n";
   }
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
      const real_t rho = mon.Step(r);
      R.W_pre += 0.5 * dt * (r.P_pre + prev.P_pre);
      R.W_nuc += 0.5 * dt * (r.P_nuc + prev.P_nuc);
      R.W_fric += 0.5 * dt * (r.P_fric + prev.P_fric);
      R.W_abs += 0.5 * dt * (r.P_abs + prev.P_abs);
      R.Emax = std::max(R.Emax, std::abs(r.E_h));
      prev = r;
      if (os)
      {
         *os << r.t << "," << r.E_h << "," << r.E_kin << "," << r.E_el << "," << r.E_sipg << ","
             << r.E_cons << "," << r.E_pen << "," << r.E_mod << "," << r.P_ext << ","
             << r.P_abs << "," << r.P_pre << "," << r.P_nuc << "," << r.P_fric << "," << rho << "\n";
      }
   }
   R.sum_rho = mon.SumRho();
   R.sum_abs = mon.SumAbsRho();
   R.fric_min = mon.MinFrictionPower();
   R.fric_max = mon.MaxFrictionPower();
   R.node_max = mon.NodalFrictionMax();
   R.node_scale = mon.NodalFrictionScale();
   R.last = r;
   R.steps = nsteps;
   return R;
}

static void TestT15(int nx, int ny, int nz, int p, real_t tfinal, const std::string &csv)
{
   if (g_rank == 0)
   {
      std::cout << "\n=== T-15 energy balance: " << nx << "x" << ny << "x" << nz << " box, p = "
                << p << ", " << tfinal << " s ===\n";
   }
   EnergyRun a = RunEnergy(nx, ny, nz, p, tfinal, 1.0, false, csv.empty() ? "" : csv + "_dt0.csv");
   EnergyRun b = RunEnergy(nx, ny, nz, p, tfinal, 0.5, false, csv.empty() ? "" : csv + "_dt1.csv");
   EnergyRun c = RunEnergy(nx, ny, nz, p, tfinal, 0.25, false, csv.empty() ? "" : csv + "_dt2.csv");
   EnergyRun nc = RunEnergy(nx, ny, nz, p, tfinal, 1.0, true);
   EnergyRun nc2 = RunEnergy(nx, ny, nz, p, tfinal, 0.5, true);
   const real_t r1 = std::abs(a.sum_rho / b.sum_rho), r2 = std::abs(b.sum_rho / c.sum_rho);
   const real_t rnc = std::abs(nc.sum_rho / nc2.sum_rho);
   const real_t scale = std::max({std::abs(a.W_pre), std::abs(a.W_fric), a.Emax});
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6);
      for (const EnergyRun *e : {&a, &b, &c})
      {
         std::cout << "  dt = " << e->dt << " (" << e->steps << " steps): E_h(T) = " << e->last.E_h
                   << " J [kin " << e->last.E_kin << ", SIPG " << e->last.E_sipg << ", fault cons "
                   << e->last.E_cons << ", fault pen " << e->last.E_pen << ", -dt^2/8 AMA "
                   << -e->last.E_mod << "]\n     work: prestress " << e->W_pre << ", nucleation "
                   << e->W_nuc << ", friction " << e->W_fric << ", dashpot " << e->W_abs
                   << "; sum rho = " << e->sum_rho << " (sum |rho| " << e->sum_abs
                   << "); P_fric (integral) in [" << e->fric_min << ", " << e->fric_max
                   << "]; nodal max (T+S).V = " << e->node_max << " (scale "
                   << e->node_scale << ")\n";
      }
      std::cout << "  accumulated residual ratios: dt0/dt1 = " << r1 << ", dt1/dt2 = " << r2
                << "; |sum rho| / work scale = " << std::abs(a.sum_rho) / scale << "\n"
                << "  negative control (fault consistency dropped from E_el): sum rho = "
                << nc.sum_rho << " at dt0, " << nc2.sum_rho << " at dt0/2 (ratio " << rnc
                << "); " << std::abs(nc.sum_rho / a.sum_rho) << " x the true residual\n";
   }
   TEST_ASSERT(r1 >= 3.5 && r1 <= 4.5 && r2 >= 3.5 && r2 <= 4.5,
               "T-15 accumulated residual falls by [3.5, 4.5] per dt halving");
   const real_t nmax = std::max({a.node_max, b.node_max, c.node_max});
   const real_t nsc = std::max({a.node_scale, b.node_scale, c.node_scale});
   TEST_ASSERT(nmax <= 1e-12 * nsc, "T-15 E4: (T_h + S).V <= 0 at every fault DOF and step");
   TEST_ASSERT(rnc < 1.5 && std::abs(nc.sum_rho) > 3.0 * std::abs(a.sum_rho),
               "T-15 negative control: without the consistency term the residual does not converge");
}

// --------------------------------------- hand loop (T-21(c), T-15 NIPG)

/// Boxes 1-14 written out by hand, through the operator's public pieces
/// (StateRate, TrialTraction, BackgroundTraction, CloseDOF, SlipLoad,
/// ApplyDamping, the mass inverses), in MFEM NewmarkSolver's arithmetic
/// order (linalg/ode.cpp:1271-1304, beta = 0, gamma = 1/2).  With K the
/// operator's own stiffness this is T-21(c)'s hand-written loop; with
/// another K it is the T-15 NIPG negative control.
struct HandLoop
{
   NitscheWaveOperator &op;
   std::function<void(const Vector &, Vector &)> K;
   Vector U, Ud, A;
   real_t t = 0.0;
   /// 0: Heun (the operator's boxes); 1: rev 2-6 ordering (closure with
   /// psi^n, then psi^{n+1} = psi^n + dt r(|V^{n+1}|, psi^n)) -- T-20's
   /// negative control; 2: Heun plus a second corrector pass (re-close with
   /// psi^{n+1}, repeat box 10) -- T-20's knob.
   int scheme = 0;
   Vector Xlast, siglast;    ///< X* in the frame (2 no), sigma_eff (no)
   real_t max_ratio = 0.0;   ///< max |V| dt / Dc over the steps

   HandLoop(NitscheWaveOperator &o, std::function<void(const Vector &, Vector &)> k)
      : op(o), K(std::move(k)) {}

   /// r = F(t) - K x + G d - C v, in the operator's order
   void Residual(const Vector &x, const Vector *v, Vector &r) const
   {
      Vector Kx, Gd, Cv;
      K(x, Kx);
      op.Elasticity().AssembleDirichletOnlyRHS(r, op.GetTime());
      r -= Kx;
      op.SlipLoad(op.Slip(), Gd);
      r += Gd;
      if (op.HasDamping() && v)
      {
         op.ApplyDamping(*v, Cv);
         r -= Cv;
      }
   }

   /// A^0 as NewmarkSolver's first pass: Mult = M^{-1}(F - K U + G d - C Udot).
   void Init(const Vector &U0, const Vector &Ud0, real_t t0)
   {
      U = U0;
      Ud = Ud0;
      t = t0;
      op.SetTime(t);
      Vector r;
      Residual(U, &Ud, r);
      A.SetSize(r.Size());
      op.Mass().Apply(r, A);
   }

   void Step()
   {
      const real_t dt = op.GetTimeStep();
      const real_t fac0 = 0.5, fac2 = 0.5, fac3 = 0.0, fac4 = 0.5;
      op.SetTime(t + dt);
      U.Add(dt, Ud);                          // 1
      U.Add(fac0 * dt * dt, A);
      Ud.Add(fac2 * dt, A);                   // 2
      // 3-11
      const int no = op.NumOwnedFaultDOFs();
      Vector d = op.Slip(), V = op.SlipRate(), psi = op.State();
      Vector dstar(3 * no), psistar(no), r0(no), T, S;
      for (int j = 0; j < no; j++)
      {
         for (int c = 0; c < 3; c++) { dstar(3 * j + c) = d(3 * j + c) + 0.5 * dt * V(3 * j + c); }
         const real_t Va = std::sqrt(V(3 * j) * V(3 * j) + V(3 * j + 1) * V(3 * j + 1)
                                     + V(3 * j + 2) * V(3 * j + 2));
         r0(j) = op.StateRate(j, Va, psi(j));
         psistar(j) = psi(j) + dt * r0(j);
      }
      op.TrialTraction(U, dstar, T);
      op.BackgroundTraction(op.GetTime(), S);
      Xlast.SetSize(2 * no);
      siglast.SetSize(no);
      for (int j = 0; j < no; j++)
      {
         real_t TS[3], Vn[3], X[2], sig;
         for (int c = 0; c < 3; c++) { TS[c] = T(3 * j + c) + S(3 * j + c); }
         if (scheme == 1)
         {
            op.CloseDOF(j, TS, psi(j), Vn, sig, X);
            const real_t Va = std::sqrt(Vn[0] * Vn[0] + Vn[1] * Vn[1] + Vn[2] * Vn[2]);
            psi(j) = psi(j) + dt * op.StateRate(j, Va, psi(j));
            max_ratio = std::max(max_ratio, Va * dt / op.FaultDc()(j));
         }
         else
         {
            op.CloseDOF(j, TS, psistar(j), Vn, sig, X);
            real_t Va = std::sqrt(Vn[0] * Vn[0] + Vn[1] * Vn[1] + Vn[2] * Vn[2]);
            real_t psi1 = psi(j) + 0.5 * dt * (r0(j) + op.StateRate(j, Va, psistar(j)));
            if (scheme == 2)
            {
               op.CloseDOF(j, TS, psi1, Vn, sig, X);
               Va = std::sqrt(Vn[0] * Vn[0] + Vn[1] * Vn[1] + Vn[2] * Vn[2]);
               psi1 = psi(j) + 0.5 * dt * (r0(j) + op.StateRate(j, Va, psi1));
            }
            psi(j) = psi1;
            max_ratio = std::max(max_ratio, Va * dt / op.FaultDc()(j));
         }
         Xlast(2 * j) = X[0];
         Xlast(2 * j + 1) = X[1];
         siglast(j) = sig;
         for (int c = 0; c < 3; c++)
         {
            V(3 * j + c) = Vn[c];
            d(3 * j + c) = dstar(3 * j + c) + 0.5 * dt * Vn[c];
         }
      }
      op.SetFaultState(d, V, psi);
      // 12-13
      Vector r, k(U.Size());
      if (op.HasDamping())
      {
         Residual(U, &Ud, r);
         op.DampedMass()->Apply(r, k);
      }
      else
      {
         Residual(U, nullptr, r);
         op.Mass().Apply(r, k);
      }
      A = k;
      U.Add(fac3 * dt * dt, k);
      Ud.Add(fac4 * dt, k);                   // 14
      t += dt;
   }
};

/// The operator's K, or the same assembly with the face integrator's
/// epsilon set to `eps` (T-15 NIPG control, decision Chunhui Zhao,
/// 2026-10-14: built inside the test, no existing file changes): the same
/// mesh, FE space, per-element (lambda, mu) coefficients of the material
/// sample, penalty factor 1, KeepNbrBlock, Assemble(0), as
/// ElasticityDomainOperator::AssembleStiffness (no Dirichlet faces here).
/// Note: DGElasticityIPCombinedIntegrator symmetrizes every face matrix
/// (dg_elasticity_ip_combined_integrator.hpp:183-195) whatever epsilon is,
/// so eps = +1 gives the symmetric part of the NIPG face form: the
/// consistency and symmetry terms cancel and only the penalty is left.
struct AssembledK
{
   std::unique_ptr<ParBilinearForm> form;
   std::unique_ptr<HypreParMatrix> K;
   AssembledK(NitscheWaveOperator &op, real_t eps)
   {
      ParMesh &pm = op.GetMesh();
      ParFiniteElementSpace &fes = op.Elasticity().GetFESpace();
      pm.ExchangeFaceNbrData();
      fes.ExchangeFaceNbrData();
      auto &smp = const_cast<NitscheMaterialSample &>(op.Sample());
      form = std::make_unique<ParBilinearForm>(&fes);
      form->KeepNbrBlock(true);
      form->AddDomainIntegrator(new ElasticityIntegrator(smp.Lambda(), smp.Mu()));
      form->AddInteriorFaceIntegrator(
         new DGElasticityIPCombinedIntegrator(smp.Lambda(), smp.Mu(), 3, eps, 1.0));
      form->Assemble(0);
      form->Finalize();
      K.reset(form->ParallelAssemble());
   }
   void Mult(const Vector &x, Vector &y) const { y.SetSize(x.Size()); K->Mult(x, y); }
};

/// T-21(c): CentralDifferenceSolver::Step vs the hand loop, bit for bit,
/// over n steps of the N-10 fixture.
static void TestT21c(int nx, int ny, int nz, int p, int nsteps)
{
   std::unique_ptr<Mesh> m(TPV102Box(nx, ny, nz));
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   TPV102Friction F1, F2;
   NitscheFaultSetup fs1 = TPV102Setup(F1.cl), fs2 = TPV102Setup(F2.cl);
   NitscheWaveOperator a(pm, p, TPV102Material(), TPV102BC(), {}, &fs1);
   NitscheWaveOperator b(pm, p, TPV102Material(), TPV102BC(), {}, &fs2);
   const real_t dt = ComputeMaxDt(a.Sample(), pm, p, kCuse[p]);
   a.SetTimeStep(dt);
   b.SetTimeStep(dt);
   TPV102InitialState(a);
   TPV102InitialState(b);
   const int n = a.Height();
   Vector U(n), Ud(n);
   U = 0.0; Ud = 0.0;
   CentralDifferenceSolver ode;
   ode.Init(a);
   real_t t = 0.0;
   ElasticityDomainOperator<ParMesh> &Eb = b.Elasticity();
   HandLoop hl(b, [&Eb](const Vector &x, Vector &y) { Eb.ApplyStiffness(x, y); });
   hl.Init(U, Ud, 0.0);
   int first_bad = -1;
   for (int s = 0; s < nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
      hl.Step();
      auto same = [](const Vector &x, const Vector &y) {
         return x.Size() == y.Size() && std::memcmp(x.GetData(), y.GetData(), sizeof(real_t) * x.Size()) == 0;
      };
      int bad = !(same(a.Slip(), b.Slip()) && same(a.SlipRate(), b.SlipRate())
                  && same(a.State(), b.State()) && same(U, hl.U) && same(Ud, hl.Ud)
                  && a.GetTime() == t && hl.t == t);
      bad = ISum(bad, pm.GetComm());
      if (bad && first_bad < 0) { first_bad = s + 1; break; }
   }
   if (g_rank == 0)
   {
      std::cout << "  T-21(c) " << nx << "x" << ny << "x" << nz << " box, p = " << p << ", "
                << nsteps << " steps (t = " << t << " s, max|V| = " << a.LastStep().max_V
                << "): first step where solver and hand loop differ: "
                << (first_bad < 0 ? std::string("none") : std::to_string(first_bad))
                << "; operator time " << a.GetTime() << ", stamp of last advance "
                << a.LastStep().t << ", t " << t << "\n";
   }
   TEST_ASSERT(first_bad < 0, "T-21(c) fault state, U, Udot after CentralDifferenceSolver::Step == hand loop, bit for bit");
   TEST_ASSERT(a.GetTime() == t && a.LastStep().t == t,
               "T-21(c) operator time == stamp of the last fault advance == t");
}

/// T-15 NIPG control: the SIPG assembly reproduces the operator's K; the
/// dynamics with the eps = +1 K break the SIPG energy balance.
static void TestT15Nipg(int nx, int ny, int nz, int p, real_t tfinal)
{
   if (g_rank == 0)
   {
      std::cout << "\n=== T-15 NIPG negative control (eps = +1 K in the step, G keeps eps = -1) ===\n";
   }
   std::unique_ptr<Mesh> m(TPV102Box(nx, ny, nz));
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   std::vector<real_t> sums;
   for (int halve = 0; halve < 2; halve++)
   {
      TPV102Friction F;
      NitscheFaultSetup fs = TPV102Setup(F.cl);
      NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
      const real_t dt = (halve ? 0.5 : 1.0) * ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
      op.SetTimeStep(dt);
      TPV102InitialState(op);
      if (halve == 0)
      {
         AssembledK Ks(op, -1.0);
         const int n = op.Height();
         real_t worst = 0.0, scale = 0.0;
         for (int r = 0; r < 3; r++)
         {
            Vector x(n), y1, y2;
            x.Randomize(100 + r + 7 * g_rank);
            x -= 0.5;
            op.Elasticity().ApplyStiffness(x, y1);
            Ks.Mult(x, y2);
            y2 -= y1;
            worst = std::max(worst, RMax(y2.Normlinf(), pm.GetComm()));
            scale = std::max(scale, RMax(y1.Normlinf(), pm.GetComm()));
         }
         if (g_rank == 0)
         {
            std::cout << "  SIPG assembly in the test vs the operator's K (3 random vectors): max |dK x| / |K x| = "
                      << worst / scale << "\n";
         }
         TEST_ASSERT(worst <= 1e-12 * scale, "T-15 NIPG control: the test's eps = -1 assembly == operator K");
      }
      AssembledK Kn(op, +1.0);
      NitscheEnergyMonitor mon(op);
      HandLoop hl(op, [&Kn](const Vector &x, Vector &y) { Kn.Mult(x, y); });
      const int n = op.Height();
      Vector U0(n), Ud0(n);
      U0 = 0.0; Ud0 = 0.0;
      hl.Init(U0, Ud0, 0.0);
      mon.Start(mon.Evaluate(hl.U, hl.Ud, hl.A, 0.0));
      const long nsteps = std::lround(std::ceil(tfinal / dt - 1e-9));
      for (long s = 0; s < nsteps; s++)
      {
         hl.Step();
         mon.Step(mon.Evaluate(hl.U, hl.Ud, hl.A, hl.t));
      }
      sums.push_back(mon.SumRho());
      if (g_rank == 0)
      {
         std::cout << "  dt = " << dt << " (" << nsteps << " steps): sum rho = " << mon.SumRho()
                   << " (sum |rho| " << mon.SumAbsRho() << ")\n";
      }
   }
   const real_t ratio = std::abs(sums[0] / sums[1]);
   if (g_rank == 0) { std::cout << "  NIPG control residual ratio dt0/dt1 = " << ratio << "\n"; }
   TEST_ASSERT(ratio < 1.5, "T-15 NIPG control: the residual does not fall with dt (balance broken)");
}

// ------------------------------------------------------------- T-20

/// One run of the coupled dt-refinement study: fault slip, state,
/// endpoint shear traction (X* + gamma V in the frame: the traction at
/// (U^{n+1}, d^{n+1}); the trial X* alone differs from it by an O(dt) term),
/// sigma_eff, and the velocity at a few volume stations, sampled every
/// `every` steps (the output times are the dt0 grid).
struct T20Series
{
   std::vector<Vector> d, psi, tau, sig, st;
   real_t max_ratio = 0.0;
};

enum class T20Scheme { Solver, Heun, OldOrder, TwoCorrector };

/// Station positions (code frame): off the fault, Y = +-1 km.
static const real_t kT20Stations[3][3] = {
   {0.0, 1000.0, -7500.0}, {6000.0, -1000.0, -7500.0}, {0.0, 1000.0, -2000.0}};

static void StationVelocity(ParGridFunction &v, const Array<int> &el,
                            const Array<IntegrationPoint> &ip, Vector &out)
{
   const int ns = el.Size();
   out.SetSize(3 * ns);
   out = 0.0;
   Vector val(3);
   for (int s = 0; s < ns; s++)
   {
      if (el[s] < 0) { continue; }
      v.GetVectorValue(el[s], ip[s], val);
      for (int c = 0; c < 3; c++) { out(3 * s + c) = val(c); }
   }
   MPI_Allreduce(MPI_IN_PLACE, out.GetData(), 3 * ns, MPI_DOUBLE, MPI_SUM,
                 v.ParFESpace()->GetComm());
}

static T20Series RunT20(ParMesh &pm, int p, real_t dt0, int refine, real_t tfinal,
                        T20Scheme scheme)
{
   T20Series R;
   TPV102Friction F;
   NitscheFaultSetup fs = TPV102Setup(F.cl);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(), {}, &fs);
   const real_t dt = dt0 / refine;
   op.SetTimeStep(dt);
   TPV102InitialState(op);
   const int no = op.NumOwnedFaultDOFs(), n = op.Height();
   // stations: each found by exactly one rank
   DenseMatrix pts(3, 3);
   for (int s = 0; s < 3; s++) for (int c = 0; c < 3; c++) { pts(c, s) = kT20Stations[s][c]; }
   Array<int> el;
   Array<IntegrationPoint> ips;
   pm.FindPoints(pts, el, ips, false);
   ParGridFunction vg(&op.Elasticity().GetFESpace());
   Vector U(n), Ud(n);
   U = 0.0; Ud = 0.0;
   const long nsteps = std::lround(std::ceil(tfinal / dt0 - 1e-9)) * refine;
   auto sample = [&](NitscheWaveOperator &o, const Vector &Udot, const Vector &Xs,
                     const Vector &sg) {
      R.d.push_back(o.Slip());
      R.psi.push_back(o.State());
      Vector tau(2 * no);
      for (int j = 0; j < no; j++)
      {
         const real_t *dv = &o.FaultFrames()(9 * j + 3), *sv = &o.FaultFrames()(9 * j + 6);
         const real_t *V = &o.SlipRate()(3 * j);
         const real_t g = 0.5 * o.FaultSigmaF()(j) * dt;
         tau(2 * j) = Xs(2 * j) + g * fault_frame::Dot(V, dv);
         tau(2 * j + 1) = Xs(2 * j + 1) + g * fault_frame::Dot(V, sv);
      }
      R.tau.push_back(tau);
      R.sig.push_back(sg);
      vg = Udot;
      Vector sv;
      StationVelocity(vg, el, ips, sv);
      R.st.push_back(sv);
   };
   if (scheme == T20Scheme::Solver)
   {
      CentralDifferenceSolver ode;
      ode.Init(op);
      real_t t = 0.0;
      for (long s = 0; s < nsteps; s++)
      {
         real_t h = dt;
         ode.Step(U, Ud, t, h);
         R.max_ratio = std::max(R.max_ratio, op.LastStep().max_V_dt_over_Dc);
         if ((s + 1) % refine == 0) { sample(op, Ud, op.LastShear(), op.LastSigmaEff()); }
      }
   }
   else
   {
      ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
      HandLoop hl(op, [&E](const Vector &x, Vector &y) { E.ApplyStiffness(x, y); });
      hl.Init(U, Ud, 0.0);
      hl.scheme = (scheme == T20Scheme::OldOrder) ? 1 : (scheme == T20Scheme::TwoCorrector) ? 2 : 0;
      for (long s = 0; s < nsteps; s++)
      {
         hl.Step();
         if ((s + 1) % refine == 0) { sample(op, hl.Ud, hl.Xlast, hl.siglast); }
      }
      R.max_ratio = std::max(R.max_ratio, RMax(hl.max_ratio, pm.GetComm()));
   }
   return R;
}

/// max over output times and DOFs of |a - b| (global)
static real_t SeriesErr(const std::vector<Vector> &a, const std::vector<Vector> &b, MPI_Comm c)
{
   real_t e = 0.0;
   const size_t n = std::min(a.size(), b.size());
   for (size_t k = 0; k < n; k++)
      for (int i = 0; i < a[k].Size(); i++) { e = std::max(e, std::abs(a[k](i) - b[k](i))); }
   return RMax(e, c);
}

static void TestT20(int nx, int ny, int nz, int p, real_t tfinal, int base = 1)
{
   std::unique_ptr<Mesh> m(TPV102Box(nx, ny, nz));
   std::unique_ptr<ParMesh> pmp = MakeParMesh(*m);
   ParMesh &pm = *pmp;
   MPI_Comm comm = pm.GetComm();
   real_t dt0;
   {
      NitscheMaterialSample smp(pm, TPV102Material());
      dt0 = ComputeMaxDt(smp, pm, p, kCuse[p]);
   }
   if (g_rank == 0)
   {
      std::cout << "\n=== T-20 coupled dt refinement: " << nx << "x" << ny << "x" << nz
                << " box, p = " << p << ", to t = " << tfinal << " s, dt0 = " << dt0
                << "; runs at dt0/" << base << ", /" << 2 * base << ", /" << 4 * base
                << " against dt0/" << 16 * base << " ===\n";
   }
   const char *names[5] = {"slip", "psi", "shear (endpoint)", "sigma_eff", "station velocity"};
   auto study = [&](T20Scheme sch, const char *label, std::vector<std::array<real_t, 5>> &errs,
                    real_t &ratio_max) {
      T20Series ref = RunT20(pm, p, dt0, 16 * base, tfinal, sch);
      ratio_max = ref.max_ratio;
      errs.clear();
      for (int k : {1, 2, 4})
      {
         T20Series r = RunT20(pm, p, dt0 * 1.0, k * base, tfinal, sch);
         ratio_max = std::max(ratio_max, r.max_ratio);
         errs.push_back({SeriesErr(r.d, ref.d, comm), SeriesErr(r.psi, ref.psi, comm),
                         SeriesErr(r.tau, ref.tau, comm), SeriesErr(r.sig, ref.sig, comm),
                         SeriesErr(r.st, ref.st, comm)});
      }
      if (g_rank == 0)
      {
         std::cout << "  " << label << " (max |V| dt / Dc over the runs = " << ratio_max << "):\n";
         for (int f = 0; f < 5; f++)
         {
            std::cout << std::setprecision(4) << "    " << std::setw(18) << names[f] << ": err "
                      << errs[0][f] << ", " << errs[1][f] << ", " << errs[2][f] << "  ratios "
                      << errs[0][f] / errs[1][f] << ", " << errs[1][f] / errs[2][f] << "\n";
         }
      }
   };
   std::vector<std::array<real_t, 5>> eS, eH, eO, eC;
   real_t rS, rH, rO, rC;
   study(T20Scheme::Solver, "operator (CentralDifferenceSolver, Heun)", eS, rS);
   study(T20Scheme::OldOrder, "negative control: rev 2-6 ordering (closure with psi^n)", eO, rO);
   study(T20Scheme::TwoCorrector, "knob: second corrector pass", eC, rC);
   auto in = [](real_t r) { return r >= 3.5 && r <= 4.5; };
   // a field whose errors are at round-off has no convergence rate: sigma_eff
   // is exact on this mirror-symmetric mesh (T-13(a)); skipped and reported
   const real_t floor[5] = {0.0, 0.0, 0.0, 1e-12 * TPV102Params::sigma_n, 0.0};
   bool ok = true;
   for (int f = 0; f < 5; f++)
   {
      if (eS[0][f] <= floor[f])
      {
         if (g_rank == 0) { std::cout << "  (" << names[f] << " errors at round-off: no rate)\n"; }
         continue;
      }
      ok = ok && in(eS[0][f] / eS[1][f]) && in(eS[1][f] / eS[2][f]);
   }
   TEST_ASSERT(ok, "T-20 every error ratio in [3.5, 4.5] (second order in the coupled state)");
   const real_t od = eO[1][0] / eO[2][0], op_ = eO[1][1] / eO[2][1];
   TEST_ASSERT(od < 3.0 && op_ < 3.0, "T-20 negative control: rev 2-6 ordering is first order in slip and psi");
   bool okc = true;
   for (int f = 0; f < 5; f++)
   {
      if (eC[0][f] <= floor[f]) { continue; }
      okc = okc && in(eC[0][f] / eC[1][f]) && in(eC[1][f] / eC[2][f]);
   }
   TEST_ASSERT(okc, "T-20 knob: second corrector keeps second order");
   TEST_ASSERT(rS <= 2.0, "T-20 predictor stability: max |V| dt / Dc <= 2");
}

// ------------------------------------------------------------------ main

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   const std::string mode = (argc > 1) ? argv[1] : "";
   if (g_rank == 0)
   {
      std::cout << "=== Nitsche P6 (TPV102 inline fixture) on " << Mpi::WorldSize()
                << " rank(s) " << mode << " ===\n";
   }
   if (mode.rfind("--abort-", 0) == 0) { return AbortMode(mode); }
   if (mode == "--t04-curved") { const int r = CurvedReport(); return r; }
   if (mode == "--t04-curved-n01") { return CurvedN01(); }
   if (mode == "--graded")
   {
      // --graded h_fine p tfinal [dump_prefix | -] [dump_every] [growth] [x_fine]
      RunConfig rc;
      rc.graded_box = true;
      rc.gspec.h = std::atof(argv[2]);
      rc.p = std::atoi(argv[3]);
      rc.tfinal = std::atof(argv[4]);
      if (argc > 5 && std::string(argv[5]) != "-") { rc.dump = argv[5]; }
      if (argc > 6) { rc.dump_every = std::atoi(argv[6]); }
      if (argc > 7) { rc.gspec.g = std::atof(argv[7]); }
      if (argc > 8) { rc.gspec.x_fine = std::atof(argv[8]); }
      TestM2(rc);
      if (g_rank == 0)
      {
         std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                   << " passed, " << num_failed << " failed ===\n";
      }
      return (num_failed == 0) ? 0 : 1;
   }
   if (mode == "--t15")
   {
      // --t15 nx ny nz p tfinal [csv_prefix]
      TestT15(std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]),
              std::atof(argv[6]), (argc > 7) ? argv[7] : "");
      if (g_rank == 0)
      {
         std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                   << " passed, " << num_failed << " failed ===\n";
      }
      return (num_failed == 0) ? 0 : 1;
   }
   if (mode == "--t21c")
   {
      // --t21c nx ny nz p nsteps
      TestT21c(std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]),
               std::atoi(argv[6]));
      return (num_failed == 0) ? 0 : 1;
   }
   if (mode == "--t15-nipg")
   {
      // --t15-nipg nx ny nz p tfinal
      TestT15Nipg(std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]),
                  std::atof(argv[6]));
      return (num_failed == 0) ? 0 : 1;
   }
   if (mode == "--interior-rank")
   {
      // P7 regression: a rank with no absorbing face must not deadlock
      // SetTimeStep (collective damped mass); 20 steps of the rupture run.
      g_interior_rank = true;
      RunConfig rc;
      rc.nx = 8; rc.ny = 8; rc.nz = 4; rc.p = 1;
      const real_t dt1 = 0.0;
      (void)dt1;
      rc.tfinal = 0.2;
      const RunResult R = RunRupture(rc);
      if (g_rank == 0)
      {
         std::cout << "  interior-rank partition: " << g_partition << "; " << R.steps
                   << " steps completed, finite " << R.finite << "\n";
      }
      TEST_ASSERT(R.finite && R.steps > 0,
                  "P7 regression: a rank with no absorbing face runs (no deadlock in SetTimeStep)");
      return (num_failed == 0) ? 0 : 1;
   }
   if (mode == "--t20")
   {
      // --t20 nx ny nz p tfinal [base]  (runs at dt0/base, /2base, /4base vs /16base)
      TestT20(std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]), std::atoi(argv[5]),
              std::atof(argv[6]), (argc > 7) ? std::atoi(argv[7]) : 1);
      if (g_rank == 0)
      {
         std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                   << " passed, " << num_failed << " failed ===\n";
      }
      return (num_failed == 0) ? 0 : 1;
   }
   if (mode == "--scatter")
   {
      for (int p = 1; p <= 3; p++) { TestScatter(p); }
      return (num_failed == 0) ? 0 : 1;
   }
   if (mode == "--onestep" && argc > 2)
   {
      const int p = (argc > 3) ? std::atoi(argv[3]) : 2;
      const int nx = (argc > 6) ? std::atoi(argv[4]) : 4;
      const int ny = (argc > 6) ? std::atoi(argv[5]) : 4;
      const int nz = (argc > 6) ? std::atoi(argv[6]) : 2;
      OneStep(argv[2], p, nx, ny, nz);
      return 0;
   }
   if (mode == "--m2")
   {
      // --m2 nx ny nz p [tfinal] [dump_prefix | -] [dt_scale]
      RunConfig rc;
      if (argc > 5)
      {
         rc.nx = std::atoi(argv[2]); rc.ny = std::atoi(argv[3]);
         rc.nz = std::atoi(argv[4]); rc.p = std::atoi(argv[5]);
      }
      if (argc > 6) { rc.tfinal = std::atof(argv[6]); }
      if (argc > 7 && std::string(argv[7]) != "-") { rc.dump = argv[7]; }
      if (argc > 8) { rc.dt_scale = std::atof(argv[8]); }
      TestM2(rc);
      if (g_rank == 0)
      {
         std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                   << " passed, " << num_failed << " failed ===\n";
      }
      return (num_failed == 0) ? 0 : 1;
   }

   if (g_rank == 0) { std::cout << "\n=== C-04 + T-04 (full traction Jacobian) ===\n"; }
   for (int p = 1; p <= 3; p++)
   {
      real_t w;
      TestC04T04(4, 4, 2, false, p, w);
      TestC04T04(4, 4, 2, true, p, w);
   }
   if (g_rank == 0) { std::cout << "\n=== fault-DOF scatter (open question 5) ===\n"; }
   for (int p = 1; p <= 3; p++) { TestScatter(p); }
   if (g_rank == 0) { std::cout << "\n=== N-07 SCEC prestress check ===\n"; }
   TestPrestress(1);
   TestPrestress(3);
   if (g_rank == 0) { std::cout << "\n=== T-21(a) ===\n"; }
   TestT21a(1);
   TestT21a(2);

   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                << " passed, " << num_failed << " failed ===\n";
   }
   return (num_failed == 0) ? 0 : 1;
}
