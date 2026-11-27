// tests/ucerf3d/test_pml_lts_3d.cpp - UCERF 3-D plan rev5, P1.4: the CFS-PML
// band and local time stepping on a 3-D mesh (tests U-P1 and U-P2).
//
// U-P1  a plane wave into the band.  Material with lambda = 0, so the P wave
//       u = (f(x - c t), 0, 0), c = sqrt(2 mu / rho), has sigma_xx only: every
//       face parallel to x is traction free, in the band as well, and the
//       natural condition there is exact.  The band sits on the two walls
//       x = +-X (PmlBand with the walls given), which are held (u = 0); the
//       other faces are natural.  The pulse f is a Gaussian of width w at
//       x = 0 moving to +x.  Once it has left the inner box, whatever is in
//       the inner box came back from the band:
//          R = max over t in the return window of max |u_x| (inner box) / A.
//       The same run without the band (the held wall alone) gives R near 1.
// U-P2  local time stepping with the band (the 3-D version of the 2-D
//       make test-lts-pml-band, tests/2d_inplane/lts_pml_band_check.py).  The
//       deck tests/ucerf3d/configs/T-pml_lts_band.toml: the cycle's dynamic
//       phase (DynPhase, with PlanLts / MapLts as the driver) on T-pml, the
//       band on the four walls and the bottom, an event started at t = 0 by
//       a V_init patch.  Three runs to the same time T:
//          glob   global stepping
//          same   [cycle.lts] rate2 with the test key same_step (the LTS
//                 run's clusters, every cluster at dt0)
//          rate2  [cycle.lts] rate2
//       Checks, with the 2-D test's tolerances:
//        1. the band sits in one cluster and the three runs have the same dt0;
//        2. the waves are inside the band at T (|v| there > 1e-4 m/s);
//        3. same against glob: slip within 1e-9 m and slip rate within 1e-9
//           m/s at every fault point, velocity within 1e-10 m/s at every
//           element centroid;
//        4. rate2 against glob: slip within 0.1 % of the largest slip, the
//           phase's peak slip rate within 0.1 %, the band's velocity within
//           5 % of its largest.
//
// Usage (from miniapps/seas):
//   mpirun -np N ./seas_test_pml_lts_3d [--only P1,P2] [--mesh ucerf3d/mesh/T-pml.msh]
//          [--deck tests/ucerf3d/configs/T-pml_lts_band.toml] [--sync N]

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"
#include "../../common/mpi_context.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int g_rank = 0, num_tests = 0, num_failed = 0;
std::set<std::string> g_only;
bool Want(const char *p) { return g_only.empty() || g_only.count(p); }
void Check(bool ok, const std::string &what)
{
   num_tests++;
   if (!ok) { num_failed++; }
   if (g_rank == 0) { std::cout << (ok ? "  PASSED: " : "  FAILED: ") << what << std::endl; }
}
real_t GMax(real_t v, MPI_Comm c)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, c);
   return v;
}

// ------------------------------------------------------------------ U-P1
struct WaveResult
{
   real_t R = 0.0, t_R = 0.0, dt = 0.0, c = 0.0, A0 = 0.0, kappa_max = 1.0;
   long steps = 0;
   bool finite = true;
   std::vector<std::array<real_t, 2>> series;   // (t, max |u_x| inner / A)
};

/// The x walls (|n_x| > 0.9) of `wall_attr` keep attribute 5; its other
/// faces get attribute 7.
void SplitWalls(Mesh &m, int wall_attr)
{
   Vector nor(3);
   for (int b = 0; b < m.GetNBE(); b++)
   {
      if (m.GetBdrAttribute(b) != wall_attr) { continue; }
      ElementTransformation *T = m.GetBdrElementTransformation(b);
      T->SetIntPoint(&Geometries.GetCenter(T->GetGeometryType()));
      CalcOrtho(T->Jacobian(), nor);
      m.SetBdrAttribute(b, std::abs(nor(0)) > 0.9 * nor.Norml2() ? 5 : 7);
   }
   m.SetAttributes();
}

/// U-G4 (P1.10): the wave run on the mfem device against the host.  One
/// device per process: the host run (--device cpu) writes U and U' of every
/// rank at three times, the device run (--device debug | cuda) compares.
struct G4Opt
{
   std::string device = "cpu", write, compare, name;
   // U' limit: measured on 1 rank (column, p = 2), device against host:
   // 1.6e-12 after 14,394 steps with the band, 5.8e-13 after 2,733 steps of
   // the same run without it (the plain device step, P1.7); U 5e-14 in both
   real_t tol = 1e-12, tol_v = 1e-11;
   bool band = true;      ///< false: the same run without the band (round-off reference)
   std::string only;      ///< "col" or "mesh": one case
   real_t snap_dv[3] = {0.0, 0.0, 0.0};
   // results (this rank's part reduced over the ranks)
   real_t du = 0.0, dv = 0.0, in_band = 0.0;
   bool on_device = false, sizes_ok = true;
   int snaps = 0;
};

/// X: the walls x = +-X; L: band thickness; Xin: half width of the inner
/// box in x (the measurement region is |x| <= Xin); w: pulse width.
/// g4: the U-G4 form: the walls behind the band are natural (the device
/// step has no Dirichlet faces before P1.11) and the run goes to the mfem
/// device when g4->device is not "cpu".
WaveResult RunWave(Mesh &sm, int p, bool with_band, real_t X, real_t L, real_t w,
                   real_t f_alpha, real_t R0, real_t cfl, bool ricker = false,
                   G4Opt *g4 = nullptr)
{
   const real_t mu = 3.2e10, rho = 2670.0;
   const real_t c = std::sqrt(2.0 * mu / rho);
   ParMesh pm(MPI_COMM_WORLD, sm);
   BoundaryConfig bc;
   bc.fault_attr = 99;   // none: the mesh's fault and creep faces are welded
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1, 2, 7};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
   if (g4)
   {
      bc.dirichlet_attrs.clear();
      bc.natural_attrs = {1, 2, 5, 7};
   }
   NitscheWaveOperator op(pm, p, MaterialField::MakeConstant(0.0, mu, rho), bc);
   WaveResult r;
   r.c = c;
   if (with_band)
   {
      PmlBandParams prm;
      prm.thickness = L;
      prm.power = 2;
      prm.R0 = R0;
      prm.c_ref = c;
      prm.f_alpha = f_alpha;
      PmlWall wp, wm;
      wp.n = {1.0, 0.0, 0.0};
      wp.h = X;
      wp.attr = 5;
      wm.n = {-1.0, 0.0, 0.0};
      wm.h = X;
      wm.attr = 5;
      op.SetPml(std::make_shared<PmlBand>(3, std::vector<PmlWall>{wp, wm}, prm));
      r.kappa_max = op.Pml()->KappaMax();
   }
   if (g4 && g4->device != "cpu")
   {
      op.SetMatrixFree(true);
      op.SetFastFaultKernels(true);
      op.SetDeviceFault(true);
   }
   const real_t ceig = NitscheStableConstant(3, p);
   r.dt = ComputeMaxDt(op.Sample(), pm, p, cfl * ceig);
   if (with_band)
   {
      // the cycle's [cycle.pml] dt_rule = "element" (DynPhase): each element's
      // own step, the band's divided by sqrt(kappa_e)
      Vector ke;
      op.Pml()->ElementKappa(ke);
      real_t dmin = std::numeric_limits<real_t>::infinity();
      for (int e = 0; e < pm.GetNE(); e++)
      {
         dmin = std::min(dmin, cfl * ceig / (real_t(p) * real_t(p)) * NitscheElementLength(pm, e)
                               / (c * std::sqrt(ke(e))));
      }
      MPI_Allreduce(MPI_IN_PLACE, &dmin, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, pm.GetComm());
      r.dt = std::min(r.dt, dmin);
   }
   op.SetTimeStep(r.dt);
   op.SetTime(0.0);
   auto &fes = op.Elasticity().GetFESpace();
   ParGridFunction gu(&fes), gv(&fes);
   const real_t A = 1e-3;
   // f(q) = A exp(-q^2) (Gaussian) or A (1 - 2 q^2) exp(-q^2) (Ricker: no
   // content at zero frequency), q = x / w; u = f(x - c t), u' = -c f'
   VectorFunctionCoefficient cu(3, [=](const Vector &x, Vector &u)
   {
      const real_t q = x(0) / w, g = std::exp(-q * q);
      u = 0.0;
      u(0) = A * (ricker ? (1.0 - 2.0 * q * q) * g : g);
   });
   VectorFunctionCoefficient cv(3, [=](const Vector &x, Vector &u)
   {
      const real_t q = x(0) / w, g = std::exp(-q * q);
      u = 0.0;
      u(0) = -c * A / w * (ricker ? (4.0 * q * q * q - 6.0 * q) * g : -2.0 * q * g);
   });
   gu.ProjectCoefficient(cu);
   gv.ProjectCoefficient(cv);
   Vector U(gu), V(gv), us(U.Size());
   if (g4)
   {
      U.UseDevice(true);
      V.UseDevice(true);
   }
   us = 0.0;
   if (with_band) { op.PmlEnter(us); }
   // inner DOFs: x component at nodes with |x| <= Xin
   const real_t Xin = X - L;
   std::vector<int> inner;
   {
      ParGridFunction gx(&fes);
      VectorFunctionCoefficient cx(3, [](const Vector &x, Vector &u) { u = x; });
      gx.ProjectCoefficient(cx);
      const int nd = fes.GetNDofs();
      for (int i = 0; i < nd; i++)
      {
         const int ix = fes.DofToVDof(i, 0);
         if (std::abs(gx(ix)) <= Xin) { inner.push_back(ix); }
      }
   }
   auto inner_max = [&](const Vector &u)
   {
      real_t m = 0.0;
      u.HostRead();
      for (int i : inner) { m = std::max(m, std::abs(u(i))); }
      return GMax(m, pm.GetComm());
   };
   // U-G4: U and U' at three times, written or compared per rank
   std::FILE *g4f = nullptr;
   long g4_at[3] = {-1, -1, -1};
   const real_t g4_v0 = c * A / w;
   auto g4_snap = [&]()
   {
      const int n = U.Size();
      const real_t *u = U.HostRead(), *v = V.HostRead();
      if (!g4->write.empty())
      {
         std::fwrite(&n, sizeof(int), 1, g4f);
         std::fwrite(u, sizeof(real_t), n, g4f);
         std::fwrite(v, sizeof(real_t), n, g4f);
         real_t m = 0.0;
         for (int i = 0; i < n; i++) { if (!std::binary_search(inner.begin(), inner.end(), i)) { m = std::max(m, std::abs(u[i])); } }
         g4->in_band = std::max(g4->in_band, m / A);
      }
      else
      {
         int n0 = -1;
         std::vector<real_t> hu(n), hv(n);
         if (std::fread(&n0, sizeof(int), 1, g4f) != 1 || n0 != n
             || std::fread(hu.data(), sizeof(real_t), n, g4f) != static_cast<size_t>(n)
             || std::fread(hv.data(), sizeof(real_t), n, g4f) != static_cast<size_t>(n))
         { g4->sizes_ok = false; return; }
         real_t m = 0.0;
         for (int i = 0; i < n; i++)
         {
            g4->du = std::max(g4->du, std::abs(u[i] - hu[i]) / A);
            g4->dv = std::max(g4->dv, std::abs(v[i] - hv[i]) / g4_v0);
            if (g4->snaps < 3)
            { g4->snap_dv[g4->snaps] = std::max(g4->snap_dv[g4->snaps], std::abs(v[i] - hv[i]) / g4_v0); }
            if (!std::binary_search(inner.begin(), inner.end(), i)) { m = std::max(m, std::abs(u[i])); }
         }
         g4->in_band = std::max(g4->in_band, m / A);
      }
      g4->snaps++;
   };
   r.A0 = A;
   CentralDifferenceSolver ode;
   ode.Init(op);
   const real_t t_out = (Xin + 3.5 * w) / c;             // the pulse has left the inner box
   const real_t t_end = (2.0 * X + 3.5 * w) / c;         // the echo's tail has passed x = 0
   r.steps = std::lround(std::ceil(t_end / r.dt));
   const long every = std::max(1L, r.steps / 60);
   if (g4)
   {
      // the pulse's centre a quarter and a half of the way into the band,
      // and at the wall, where the run ends
      g4_at[0] = std::lround((Xin + 0.25 * L) / c / r.dt);
      g4_at[1] = std::lround((Xin + 0.5 * L) / c / r.dt);
      g4_at[2] = std::lround(X / c / r.dt);
      r.steps = g4_at[2];
      if (g_rank == 0) { std::printf("  U-G4: %ld steps, samples at steps %ld, %ld, %ld\n", r.steps, g4_at[0], g4_at[1], g4_at[2]); }
      const std::string fn = (g4->write.empty() ? g4->compare : g4->write) + "." + g4->name
                             + "." + std::to_string(Mpi::WorldRank());
      g4f = std::fopen(fn.c_str(), g4->write.empty() ? "rb" : "wb");
      MFEM_VERIFY(g4f, "U-G4: cannot open " << fn);
   }
   if (g_rank == 0)
   {
      std::printf("  %s: %d DOFs here, dt = %.5f s, %ld steps to t = %.1f s (kappa_max %.3f)\n",
                  with_band ? "band" : "held wall", U.Size(), r.dt, r.steps, t_end, r.kappa_max);
      std::fflush(stdout);
   }
   real_t t = 0.0;
   for (long s = 1; s <= r.steps; s++)
   {
      real_t h = r.dt;
      ode.Step(U, V, t, h);
      if (g4 && (s == g4_at[0] || s == g4_at[1] || s == g4_at[2])) { g4_snap(); }
      if (s % every == 0 || s == r.steps)
      {
         const real_t m = inner_max(U) / r.A0;
         r.series.push_back({t, m});
         if (g_rank == 0 && std::getenv("PL3_TRACE"))
         { std::printf("    step %ld t %.2f inner %.3e\n", s, t, m); std::fflush(stdout); }
         if (!std::isfinite(m)) { r.finite = false; break; }
         if (t >= t_out && m > r.R) { r.R = m; r.t_R = t; }
      }
   }
   if (g4)
   {
      std::fclose(g4f);
      g4->on_device = op.PmlOnDevice();
      g4->du = GMax(g4->du, pm.GetComm());
      g4->dv = GMax(g4->dv, pm.GetComm());
      g4->in_band = GMax(g4->in_band, pm.GetComm());
      for (real_t &d : g4->snap_dv) { d = GMax(d, pm.GetComm()); }
      int ok = (g4->sizes_ok && g4->snaps == 3) ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, pm.GetComm());
      g4->sizes_ok = (ok == 1);
   }
   return r;
}

void PrintWave(const char *name, const WaveResult &r)
{
   if (g_rank != 0) { return; }
   std::printf("  %s: c = %.1f m/s, dt = %.5f s, %ld steps, kappa_max %.3f; R = %.4e at t = "
               "%.2f s\n    max |u_x| (inner) / A against t:", name, r.c, r.dt, r.steps,
               r.kappa_max, r.R, r.t_R);
   for (size_t k = 0; k < r.series.size(); k += 3)
   { std::printf(" %.1f:%.2e", r.series[k][0], r.series[k][1]); }
   std::printf("\n");
}

void TestP1(const char *mesh_file, int p, real_t w)
{
   const real_t X = 90e3, L = 30e3, fa = 0.01, R0 = 1e-3;
   if (g_rank == 0)
   {
      std::cout << "\n=== U-P1 plane P wave into the band, " << mesh_file << ", p = " << p
                << ", pulse width " << w / 1e3 << " km ===" << std::endl;
   }
   Mesh sm(mesh_file, 1, 1);
   SplitWalls(sm, 5);
   WaveResult band = RunWave(sm, p, true, X, L, w, fa, R0, 0.5);
   PrintWave("band", band);
   WaveResult bare = RunWave(sm, p, false, X, L, w, fa, R0, 0.5);
   PrintWave("held wall, no band", bare);
   std::ostringstream os;
   os << "U-P1 p = " << p << ": the band returns R = " << band.R << " of the pulse (design R0 = "
      << R0 << "); the held wall alone returns " << bare.R;
   Check(band.finite && bare.finite && band.R <= R0 && bare.R > 0.5, os.str());
}

/// The same wave on a uniform column [-X, X] x [0, B]^2 of tets, h across.
WaveResult ColumnWave(int p, real_t h, real_t w, bool with_band, real_t fa, real_t R0, bool ricker)
{
   const real_t X = 90e3, L = 30e3, B = h;
   const int nx = static_cast<int>(std::lround(2 * X / h));
   Mesh sm = Mesh::MakeCartesian3D(nx, 1, 1, Element::TETRAHEDRON, 2 * X, B, B);
   for (int v = 0; v < sm.GetNV(); v++) { sm.GetVertex(v)[0] -= X; }
   for (int b = 0; b < sm.GetNBE(); b++) { sm.SetBdrAttribute(b, 5); }
   sm.SetAttributes();
   SplitWalls(sm, 5);
   return RunWave(sm, p, with_band, X, L, w, fa, R0, 0.5, ricker);
}

/// The design echo of the pulse (PML plan rev2, (2.13)-(2.15)): a plane wave
/// of angular frequency omega = c k that crosses the band, turns at the held
/// wall and crosses back returns with the factor
///    Rc(omega) = - exp(Lam (i omega alpha - omega^2) / (alpha^2 + omega^2)),
/// Lam = ln(1 / R0) (the profile's integral: 2 int d / c = Lam), so |Rc| =
/// R0^(omega^2 / (omega^2 + alpha^2)): R0 well above f_alpha, 1 at zero
/// frequency.  The echo is u(xi) = 2 Re int_0^inf F(k) Rc(c k) e^(i k xi) dk,
/// xi = 2 X - x - c t; returned: max |u| / A over the xi the test's window
/// and inner box see.
real_t DesignEcho(real_t c, real_t w, real_t X, real_t Xin, real_t R0, real_t fa, bool ricker)
{
   const real_t al = 2.0 * M_PI * fa, Lam = std::log(1.0 / R0);
   const int nk = 40000, nx = 3000;
   const real_t kmax = 12.0 / w, dk = kmax / nk;
   const real_t xa = -Xin - 3.5 * w, xb = 2.0 * X - 3.5 * w;
   std::vector<std::complex<real_t>> G(nk);
   for (int i = 0; i < nk; i++)
   {
      const real_t k = (i + 0.5) * dk, om = c * k, q = 0.25 * k * k * w * w;
      const real_t F = w / (2.0 * std::sqrt(M_PI)) * (ricker ? 2.0 * q : 1.0) * std::exp(-q);
      const real_t den = al * al + om * om;
      G[i] = -F * std::exp(std::complex<real_t>(-Lam * om * om / den, Lam * om * al / den));
   }
   real_t best = 0.0;
   for (int j = 0; j <= nx; j++)
   {
      const real_t xi = xa + (xb - xa) * j / nx;
      real_t u = 0.0;
      for (int i = 0; i < nk; i++)
      {
         const real_t k = (i + 0.5) * dk;
         u += G[i].real() * std::cos(k * xi) - G[i].imag() * std::sin(k * xi);
      }
      best = std::max(best, std::abs(2.0 * u * dk));
   }
   return best;
}

void TestP1C(int p, real_t h, real_t w, real_t fa)
{
   const real_t R0 = 1e-3;
   for (int ricker = 1; ricker >= 0; ricker--)
   {
      const char *name = ricker ? "Ricker" : "Gaussian";
      if (g_rank == 0)
      {
         std::cout << "\n=== U-P1 plane P wave into the band: column of tets, h = " << h / 1e3
                   << " km (" << 30e3 / h << " across the band), p = " << p << ", " << name
                   << " pulse, width " << w / 1e3 << " km, f_alpha " << fa << " Hz ===" << std::endl;
      }
      WaveResult band = ColumnWave(p, h, w, true, fa, R0, ricker);
      PrintWave("band", band);
      WaveResult bare = ColumnWave(p, h, w, false, fa, R0, ricker);
      PrintWave("held wall, no band", bare);
      const real_t Rd = DesignEcho(band.c, w, 90e3, 60e3, R0, fa, ricker);
      std::ostringstream os;
      if (ricker)
      {
         os << "U-P1 " << name << " pulse (no content below f_alpha): the band returns R = "
            << band.R << " <= the design R0 = " << R0 << " (design echo of this pulse " << Rd
            << "); the held wall alone returns " << bare.R;
         Check(band.finite && bare.finite && band.R <= R0 && bare.R > 0.5, os.str());
      }
      os.str("");
      // Gaussian: 0.1 % at h = 1.5 km, 0.3 % at 3 km.  Ricker (shorter
      // waves): 12 % at 1.5 km, 98 % at 3 km (the mesh's own echo, about 1e-3
      // at 10 elements across the band)
      const real_t tol = ricker ? 0.20 : 0.05;
      os << "U-P1 " << name << " pulse: R = " << band.R << " is the CFS design echo " << Rd
         << " (R0^(w^2 / (w^2 + alpha^2)) per frequency) within " << 100 * tol << " % ("
         << 100 * std::abs(band.R - Rd) / Rd << " %)";
      Check(band.finite && std::abs(band.R - Rd) <= tol * Rd, os.str());
   }
}

// ------------------------------------------------------------------ U-P2
struct PhaseResult
{
   // rank 0, after the gather: rows of 9 (x, slip, V) and of 7 (centroid, v, in band)
   std::vector<double> fault, vol;
   real_t dt0 = 0.0, T = 0.0, peak_V = 0.0, sync = 0.0, dt = 0.0, dt_cfl = 0.0;
   int pml_cluster = -1, num_clusters = 1;
   long fault_steps = 0;
   double seconds = 0.0;
   std::string report;
};

void GatherRows(MPI_Comm comm, const std::vector<double> &loc, std::vector<double> &all)
{
   int np = 1, rank = 0;
   MPI_Comm_size(comm, &np);
   MPI_Comm_rank(comm, &rank);
   int n = static_cast<int>(loc.size());
   std::vector<int> cnt(np), dsp(np, 0);
   MPI_Gather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, comm);
   if (rank == 0)
   {
      for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
      all.assign(dsp[np - 1] + cnt[np - 1], 0.0);
   }
   MPI_Gatherv(loc.data(), n, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(), MPI_DOUBLE, 0, comm);
}

/// mode: "glob", "same" or "rate2".  n_sync: the run's length in sync
/// intervals of the rate2 run (T = n_sync * sync); sync_rate2 <= 0: this is
/// the rate2 run and it sets T.
PhaseResult RunPhase(MPIContext &mpi, cycle::CycleConfig cfg, const std::string &mode,
                     real_t T_want, real_t dt_glob = 0.0)
{
   cfg.lts.mode = (mode == "glob") ? "off" : "rate2";
   cfg.lts.same_step = (mode == "same");
   cfg.lts.static_split = false;
   PhaseResult r;
   std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
   cycle::LtsLayout lts = cycle::DynPhase::PlanLts(*smesh, cfg, mpi.GetComm());
   ParMesh pmesh(mpi.GetComm(), *smesh, lts.part.empty() ? nullptr : lts.part.data());
   cycle::DynPhase::MapLts(lts, *smesh, pmesh);
   r.report = lts.report;
   r.pml_cluster = lts.on ? lts.pml_cluster : -1;
   r.num_clusters = lts.on ? lts.num_clusters : 1;
   smesh.reset();
   // the QD regime's operator, built but not set up: the fault DOF table, the
   // frames and the band, as the driver hands them to the dynamic phase
   cycle::QDPhase qd(pmesh, mpi, cfg);
   const spatial::SpatialFrictionConfig &s = cfg.shared;
   const MaterialField material = MaterialField::MakeConstant(
      s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
   const spatial::PorePressureSpec no_pp;
   const cycle::FaultParams par = cycle::BuildFaultParams(
      qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
   cycle::PrestressConfig pc = cfg.prestress;
   pc.a_max = cycle::GlobalMaxA(par, mpi.GetComm());
   const cycle::Prestress pre = cycle::BuildPrestress(par, pc);
   auto make = [&](real_t scale)
   {
      return std::make_unique<cycle::DynPhase>(pmesh, cfg, qd.Table(), par, pre, scale,
                                               lts.on ? &lts : nullptr, qd.Domain().FrameField(),
                                               nullptr, qd.PmlBandPtr(), qd.PmlWallSpeedFactor());
   };
   std::unique_ptr<cycle::DynPhase> dynp = make(1.0);
   r.dt_cfl = dynp->DtCFL();
   if (mode == "glob" && dt_glob > 0.0)
   {
      // global stepping at the LTS run's dt0 (PlanLts bounds the band's
      // kappa_e from above by its vertex values, so dt0 <= the global rule's)
      MFEM_VERIFY(dt_glob <= r.dt_cfl * (1.0 + 1e-12), "U-P2: dt0 = " << dt_glob
                  << " s is above the global step " << r.dt_cfl << " s");
      dynp.reset();
      dynp = make(dt_glob / r.dt_cfl);
   }
   cycle::DynPhase &dyn = *dynp;
   NitscheWaveOperator &W = dyn.Op();
   const int N = qd.Table().Size();
   Vector sk(2 * N), d, V, psi, us(W.Height());
   sk = 0.0;
   us = 0.0;
   cycle::bridge::ToNitsche(dyn.Frames(), sk, par.V_init, pre.psi0, d, V, psi);
   dyn.Enter(us, d, V, psi, pre.sigma0, 0.0);
   r.dt = dyn.Dt();
   r.dt0 = lts.on ? lts.dt0 : dyn.Dt();
   r.sync = dyn.StepDt();
   real_t peak = 0.0;
   dyn.SetFaultStepHook([&](real_t)
   {
      const Vector &Vf = W.SlipRate();
      for (int j = 0; j < Vf.Size() / 3; j++)
      {
         peak = std::max(peak, std::hypot(Vf(3 * j), Vf(3 * j + 1), Vf(3 * j + 2)));
      }
   });
   const long nsteps = std::lround(T_want / r.sync);
   MFEM_VERIFY(std::abs(nsteps * r.sync - T_want) <= 1e-9 * T_want, "U-P2 " << mode << ": T = "
               << T_want << " s is not a whole number of its steps of " << r.sync << " s");
   StopWatch sw;
   sw.Start();
   for (long k = 0; k < nsteps; k++)
   {
      dyn.Step();
      if (mpi.IsRoot() && std::getenv("PL3_TRACE") && (k + 1) % std::max(1L, nsteps / 20) == 0)
      {
         std::printf("    %s step %ld of %ld, tau %.3f s, %.0f s wall, rank-0 peak V so far %.4f, "
                     "now %.4f m/s\n", mode.c_str(), k + 1, nsteps, dyn.Tau(), sw.RealTime(),
                     peak, W.SlipRate().Normlinf());
         std::fflush(stdout);
      }
   }
   sw.Stop();
   r.seconds = sw.RealTime();
   r.T = dyn.Tau();
   r.peak_V = GMax(peak, mpi.GetComm());
   // fault rows
   std::vector<double> fl(9 * static_cast<size_t>(N));
   const Vector &X = W.FaultPositions(), &D = W.Slip(), &Vf = W.SlipRate();
   for (int j = 0; j < N; j++)
      for (int c = 0; c < 3; c++)
      {
         fl[9 * j + c] = X(3 * j + c);
         fl[9 * j + 3 + c] = D(3 * j + c);
         fl[9 * j + 6 + c] = Vf(3 * j + c);
      }
   GatherRows(mpi.GetComm(), fl, r.fault);
   // volume rows: velocity at the element centroids
   ParGridFunction gv(&W.Elasticity().GetFESpace());
   gv = dyn.Ud();
   std::vector<double> vl(7 * static_cast<size_t>(pmesh.GetNE()));
   Vector x(3), v(3);
   const PmlBand &band = *qd.PmlBandPtr();
   for (int e = 0; e < pmesh.GetNE(); e++)
   {
      ElementTransformation *Tr = pmesh.GetElementTransformation(e);
      const IntegrationPoint &ip = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetIntPoint(&ip);
      Tr->Transform(ip, x);
      gv.GetVectorValue(e, ip, v);
      for (int c = 0; c < 3; c++) { vl[7 * e + c] = x(c); vl[7 * e + 3 + c] = v(c); }
      vl[7 * e + 6] = band.InBand(x.GetData()) ? 1.0 : 0.0;
   }
   GatherRows(mpi.GetComm(), vl, r.vol);
   if (mpi.IsRoot())
   {
      std::printf("  %-5s: dt %.6f s, dt0 %.6f s, sync %.6f s, T = %.4f s, peak V %.6f m/s, "
                  "%.1f s wall\n", mode.c_str(), r.dt, r.dt0, r.sync, r.T, r.peak_V, r.seconds);
      if (lts.on) { std::cout << "         LTS: " << r.report << std::endl; }
   }
   return r;
}

using Key = std::array<long, 3>;
Key KeyOf(const double *x) { return {std::lround(x[0] * 100), std::lround(x[1] * 100), std::lround(x[2] * 100)}; }

/// Per point (coordinates to 1 cm): the mean of the rows' values (several
/// fault DOFs sit on one vertex).
std::map<Key, std::array<double, 6>> MeanRows(const std::vector<double> &rows, int width, int nval)
{
   std::map<Key, std::array<double, 7>> acc;
   for (size_t i = 0; i + width <= rows.size(); i += width)
   {
      auto &a = acc[KeyOf(&rows[i])];
      for (int c = 0; c < nval; c++) { a[c] += rows[i + 3 + c]; }
      a[6] += 1.0;
   }
   std::map<Key, std::array<double, 6>> out;
   for (auto &kv : acc)
   {
      std::array<double, 6> m{};
      for (int c = 0; c < nval; c++) { m[c] = kv.second[c] / kv.second[6]; }
      out[kv.first] = m;
   }
   return out;
}

void TestP2(MPIContext &mpi, const char *deck, int n_sync)
{
   if (mpi.IsRoot()) { std::cout << "\n=== U-P2 LTS with the band against global stepping, " << deck << " ===" << std::endl; }
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   // the rate2 run's sync interval sets T
   PhaseResult probe;
   {
      cycle::CycleConfig c2 = cfg;
      c2.lts.mode = "rate2";
      std::unique_ptr<Mesh> sm = cycle::LoadScaledMesh(c2.shared.mesh.path, c2.mesh_scale);
      cycle::LtsLayout l = cycle::DynPhase::PlanLts(*sm, c2, mpi.GetComm());
      probe.sync = std::ldexp(l.dt0, l.num_clusters - 1);
   }
   const real_t T = n_sync * probe.sync;
   if (mpi.IsRoot()) { std::printf("  T = %d sync intervals of %.6f s = %.4f s\n", n_sync, probe.sync, T); }
   PhaseResult lt = RunPhase(mpi, cfg, "rate2", T);
   if (std::getenv("PL3_RATE2_ONLY")) { return; }
   PhaseResult gl = RunPhase(mpi, cfg, "glob", T, lt.dt0);
   PhaseResult sa = RunPhase(mpi, cfg, "same", T);
   if (!mpi.IsRoot()) { return; }
   std::ostringstream os;
   // 1
   os << "U-P2 the band is in one cluster (" << lt.pml_cluster << " of " << lt.num_clusters
      << "), the three runs step at the same dt0 (" << gl.dt << ", " << sa.dt0 << ", " << lt.dt0
      << " s) and dt0 is within 10 % below the global rule's step " << gl.dt_cfl << " s";
   Check(lt.pml_cluster >= 0 && lt.num_clusters > 1 && std::abs(gl.dt / lt.dt0 - 1.0) <= 1e-14
         && sa.dt0 == lt.dt0 && lt.dt0 <= gl.dt_cfl * (1.0 + 1e-12) && lt.dt0 >= 0.9 * gl.dt_cfl,
         os.str());
   // maps
   auto fg = MeanRows(gl.fault, 9, 6), fs = MeanRows(sa.fault, 9, 6), fl = MeanRows(lt.fault, 9, 6);
   auto vg = MeanRows(gl.vol, 7, 3), vs = MeanRows(sa.vol, 7, 3), vl = MeanRows(lt.vol, 7, 3);
   std::map<Key, bool> inband;
   for (size_t i = 0; i + 7 <= gl.vol.size(); i += 7) { inband[KeyOf(&gl.vol[i])] = gl.vol[i + 6] > 0.5; }
   // 2
   double vband = 0.0, vall = 0.0;
   for (auto &kv : vg)
   {
      const double m = std::hypot(kv.second[0], kv.second[1], kv.second[2]);
      vall = std::max(vall, m);
      if (inband[kv.first]) { vband = std::max(vband, m); }
   }
   os.str("");
   os << "U-P2 the waves are inside the band at T = " << gl.T << " s: |v| up to " << vband
      << " m/s there (" << vall << " m/s in the whole mesh)";
   Check(vband > 1e-4, os.str());
   // 3 and 4
   auto diff = [](const std::map<Key, std::array<double, 6>> &a,
                  const std::map<Key, std::array<double, 6>> &b, int c0, int c1,
                  const std::map<Key, bool> *only, double &amax, bool &same_pts)
   {
      double dmax = 0.0;
      amax = 0.0;
      same_pts = (a.size() == b.size());
      for (auto &kv : a)
      {
         if (only && !only->at(kv.first)) { continue; }
         auto it = b.find(kv.first);
         if (it == b.end()) { same_pts = false; continue; }
         double d2 = 0.0, a2 = 0.0;
         for (int c = c0; c < c1; c++)
         {
            d2 += (kv.second[c] - it->second[c]) * (kv.second[c] - it->second[c]);
            a2 += kv.second[c] * kv.second[c];
         }
         dmax = std::max(dmax, std::sqrt(d2));
         amax = std::max(amax, std::sqrt(a2));
      }
      return dmax;
   };
   double smax = 0.0, vfmax = 0.0, vvmax = 0.0, bmax = 0.0;
   bool p1 = true, p2 = true, p3 = true, p4 = true;
   const double ds_same = diff(fg, fs, 0, 3, nullptr, smax, p1);
   const double dv_same = diff(fg, fs, 3, 6, nullptr, vfmax, p2);
   const double dvol_same = diff(vg, vs, 0, 3, nullptr, vvmax, p3);
   os.str("");
   os << "U-P2 same_step = global stepping: slip within 1e-9 m (" << ds_same << "), slip rate "
      << "within 1e-9 m/s (" << dv_same << "), velocity at " << vg.size()
      << " centroids within 1e-10 m/s (" << dvol_same << "); " << fg.size() << " fault points";
   Check(p1 && p2 && p3 && sa.T == gl.T && ds_same <= 1e-9 && dv_same <= 1e-9 && dvol_same <= 1e-10, os.str());
   const double ds_lts = diff(fg, fl, 0, 3, nullptr, smax, p1);
   const double db_lts = diff(vg, vl, 0, 3, &inband, bmax, p4);
   os.str("");
   os << "U-P2 rate2 against global stepping: slip within 0.1 % of the largest (" << ds_lts
      << " m of " << smax << " m = " << 100 * ds_lts / smax << " %)";
   Check(p1 && std::abs(lt.T - gl.T) <= 1e-9 * gl.T && ds_lts <= 1e-3 * smax, os.str());
   os.str("");
   os << "U-P2 rate2: peak slip rate " << lt.peak_V << " m/s against " << gl.peak_V
      << " (within 0.1 %: " << 100 * std::abs(lt.peak_V - gl.peak_V) / gl.peak_V << " %)";
   Check(std::abs(lt.peak_V - gl.peak_V) <= 1e-3 * gl.peak_V, os.str());
   os.str("");
   os << "U-P2 rate2: the band's velocity within 5 % of its largest (" << db_lts << " m/s of "
      << bmax << " m/s = " << 100 * db_lts / bmax << " %)";
   Check(p4 && db_lts <= 0.05 * bmax, os.str());
   if (const char *dir = std::getenv("PL3_DUMP"))
   {
      // x y z |slip| glob, rate2 per fault point; x y z |v| glob, rate2, in band per centroid
      std::ofstream ff(std::string(dir) + "/p2_fault.txt"), fv(std::string(dir) + "/p2_volume.txt");
      ff << std::setprecision(10);
      fv << std::setprecision(10);
      for (auto &kv : fg)
      {
         const auto &b = fl.at(kv.first);
         ff << kv.first[0] / 100.0 << " " << kv.first[1] / 100.0 << " " << kv.first[2] / 100.0 << " "
            << std::hypot(kv.second[0], kv.second[1], kv.second[2]) << " "
            << std::hypot(b[0], b[1], b[2]) << "\n";
      }
      for (auto &kv : vg)
      {
         const auto &b = vl.at(kv.first);
         fv << kv.first[0] / 100.0 << " " << kv.first[1] / 100.0 << " " << kv.first[2] / 100.0 << " "
            << std::hypot(kv.second[0], kv.second[1], kv.second[2]) << " "
            << std::hypot(b[0], b[1], b[2]) << " " << (inband[kv.first] ? 1 : 0) << "\n";
      }
   }
   std::printf("  wall time: global %.1f s, same_step %.1f s, rate2 %.1f s (global / rate2 = %.2f)\n",
               gl.seconds, sa.seconds, lt.seconds, gl.seconds / lt.seconds);
}

// ------------------------------------------------------------------ U-G4
/// U-G4 K (P1.10): the matrix-free product with the band's penalty factor
/// against the assembled K_p (ElasticityDomainOperator::SetPmlBand,
/// penalty_only) on a random field, and against the product without a band
/// (the factor must matter).
void TestG4K(const char *mesh_file, int p)
{
   const real_t X = 90e3, L = 30e3, mu = 3.2e10, rho = 2670.0;
   if (g_rank == 0)
   {
      std::cout << "\n=== U-G4 K: matrix-free K_p with the band against the assembled one, "
                << mesh_file << ", p = " << p << " ===" << std::endl;
   }
   Mesh sm(mesh_file, 1, 1);
   SplitWalls(sm, 5);
   ParMesh pm(MPI_COMM_WORLD, sm);
   BoundaryConfig bc;
   bc.fault_attr = 99;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1, 2, 7};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
   NitscheWaveOperator op(pm, p, MaterialField::MakeConstant(0.6 * mu, mu, rho), bc);
   NitscheWaveOperator plain(pm, p, MaterialField::MakeConstant(0.6 * mu, mu, rho), bc);
   PmlBandParams prm;
   prm.thickness = L;
   prm.power = 2;
   prm.R0 = 1e-3;
   prm.c_ref = std::sqrt(2.6 * mu / rho);
   prm.f_alpha = 0.01;
   PmlWall wp, wm;
   wp.n = {1.0, 0.0, 0.0};
   wp.h = X;
   wp.attr = 5;
   wm.n = {-1.0, 0.0, 0.0};
   wm.h = X;
   wm.attr = 5;
   op.SetPml(std::make_shared<PmlBand>(3, std::vector<PmlWall>{wp, wm}, prm));
   const real_t kmax = op.Pml()->KappaMax();
   op.SetMatrixFree(true);
   plain.SetMatrixFree(true);
   const int n = op.Height();
   Vector x(n), Kmf, Kas, K0;
   x.Randomize(4711 + g_rank);
   op.StiffnessMult(x, Kmf);
   op.Elasticity().ApplyStiffness(x, Kas);
   plain.StiffnessMult(x, K0);
   Kmf.HostRead();
   K0.HostRead();
   real_t d = 0.0, m = 0.0, d0 = 0.0;
   for (int i = 0; i < n; i++)
   {
      d = std::max(d, std::abs(Kmf(i) - Kas(i)));
      m = std::max(m, std::abs(Kas(i)));
      d0 = std::max(d0, std::abs(K0(i) - Kas(i)));
   }
   d = GMax(d, pm.GetComm());
   m = GMax(m, pm.GetComm());
   d0 = GMax(d0, pm.GetComm());
   std::ostringstream os;
   os << "U-G4 K p = " << p << ": matrix-free K_p x against the assembled K_p x: " << d / m
      << " of the largest entry (limit 1e-13); kappa_max " << kmax << "; without the band's "
      "factor the product is off by " << d0 / m;
   Check(op.MatrixFreeActive() && d <= 1e-13 * m && kmax > 1.0 && d0 > 1e-6 * m, os.str());
}

/// U-G4 W (P1.10): the plane wave into the band (U-P1's run, the walls
/// behind the band natural) on the device against the host: U and U' of
/// every DOF at three times (the pulse's centre a quarter and a half of the
/// way into the band, and at the wall, where the run ends), |dU| / A and
/// |dU'| / (c A / w) <= tol_v.  What the band returns is U-P1's check.
/// Cases: the column of tets (p = 2, Ricker) and T-pml (p = 1, Gaussian).
void TestG4W(const char *mesh_file, G4Opt g4, real_t col_h, real_t width)
{
   const real_t X = 90e3, L = 30e3, fa = 0.01, R0 = 1e-3;
   for (int k = 0; k < 2; k++)
   {
      G4Opt o = g4;
      o.name = k == 0 ? "G4Wcol" : "G4Wmesh";
      if (!o.only.empty() && o.only != (k == 0 ? "col" : "mesh")) { continue; }
      const int p = k == 0 ? 2 : 1;
      if (g_rank == 0)
      {
         std::cout << "\n=== U-G4 W plane P wave into the band, device " << o.device << ", "
                   << (k == 0 ? "column of tets" : mesh_file) << ", p = " << p << " ===" << std::endl;
      }
      WaveResult r;
      if (k == 0)
      {
         const int nx = static_cast<int>(std::lround(2 * X / col_h));
         Mesh sm = Mesh::MakeCartesian3D(nx, 1, 1, Element::TETRAHEDRON, 2 * X, col_h, col_h);
         for (int v = 0; v < sm.GetNV(); v++) { sm.GetVertex(v)[0] -= X; }
         for (int b = 0; b < sm.GetNBE(); b++) { sm.SetBdrAttribute(b, 5); }
         sm.SetAttributes();
         SplitWalls(sm, 5);
         r = RunWave(sm, p, o.band, X, L, width, fa, R0, 0.5, true, &o);
      }
      else
      {
         Mesh sm(mesh_file, 1, 1);
         SplitWalls(sm, 5);
         r = RunWave(sm, p, o.band, X, L, width, fa, R0, 0.5, false, &o);
      }
      PrintWave("band", r);
      std::ostringstream os;
      const bool dev = o.device != "cpu";
      os << "U-G4 W " << o.name << " (" << r.steps << " steps, kappa_max " << r.kappa_max
         << "): the band's memory steps "
         << (!o.band ? "nowhere (no band)" : dev ? "in device kernels" : "on the host")
         << "; the wave is inside the band at the samples (max |u| there " << o.in_band << " A)";
      Check(r.finite && o.sizes_ok && o.on_device == (dev && o.band) && o.in_band > 0.1, os.str());
      if (dev)
      {
         os.str("");
         os << "U-G4 W " << o.name << ": device against host at three times: U within " << o.du
            << " of A, U' within " << o.dv << " of c A / w (limits " << o.tol << " and " << o.tol_v << "; U' at the three "
            "times " << o.snap_dv[0] << ", " << o.snap_dv[1] << ", " << o.snap_dv[2] << ")";
         Check(o.sizes_ok && o.du <= o.tol && o.dv <= o.tol_v, os.str());
      }
   }
}

}  // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   g_rank = Mpi::WorldRank();
   const char *only = "";
   const char *mesh_file = "ucerf3d/mesh/T-pml.msh";
   int order = 2;
   real_t width = 12e3, col_h = 1.5e3, col_fa = 0.01;
   const char *deck = "tests/ucerf3d/configs/T-pml_lts_band.toml";
   int n_sync = 60;
   const char *dev_name = "cpu", *g4_write = "", *g4_compare = "";
   real_t g4_tol = 1e-12, g4_tol_v = 1e-11;
   const char *g4_only = "";
   int g4_band = 1;
   OptionsParser args(argc, argv);
   args.AddOption(&g4_tol_v, "-g4tv", "--g4-tol-v", "U-G4 W: limit for U'.");
   args.AddOption(&g4_only, "-g4c", "--g4-case", "U-G4 W: col or mesh (default both).");
   args.AddOption(&g4_band, "-g4b", "--g4-band", "U-G4 W: 0 runs without the band (reference).");
   args.AddOption(&dev_name, "-dev", "--device", "U-G4: cpu (writes) or debug / cuda (compares).");
   args.AddOption(&g4_write, "-wr", "--write", "U-G4: file prefix the host run writes.");
   args.AddOption(&g4_compare, "-cmp", "--compare", "U-G4: file prefix the device run compares.");
   args.AddOption(&g4_tol, "-g4t", "--g4-tol", "U-G4 W: limit of the device against the host.");
   args.AddOption(&only, "-o", "--only", "Comma-separated parts: P1, P2.");
   args.AddOption(&deck, "-d", "--deck", "The U-P2 deck.");
   args.AddOption(&col_h, "-ch", "--column-h", "U-P1 column: element size (m).");
   args.AddOption(&col_fa, "-cf", "--column-f-alpha", "U-P1 column: f_alpha (Hz).");
   args.AddOption(&n_sync, "-n", "--sync", "U-P2: the run's length in sync intervals of rate2.");
   args.AddOption(&mesh_file, "-m", "--mesh", "The T-pml mesh.");
   args.AddOption(&order, "-p", "--order", "Polynomial order of U-P1.");
   args.AddOption(&width, "-w", "--width", "Pulse width w of U-P1 (m).");
   args.Parse();
   if (!args.Good()) { if (g_rank == 0) { args.PrintUsage(std::cout); } return 1; }
   {
      std::stringstream ss(only);
      std::string tok;
      while (std::getline(ss, tok, ',')) { if (!tok.empty()) { g_only.insert(tok); } }
   }
   if (g_rank == 0)
   {
      std::cout << "test_pml_lts_3d (UCERF 3-D plan rev5, P1.4) on " << Mpi::WorldSize()
                << " rank(s)" << std::endl;
   }
   // U-G4: the mfem device of this process (none for the host runs)
   std::unique_ptr<Device> device;
   if (std::string(dev_name) != "cpu")
   {
      device.reset(new Device(dev_name));
#if defined(HYPRE_USING_GPU) && (MFEM_HYPRE_VERSION >= 23100)
      HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);   // as seas_nitsche_dyn_driver
      HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
#endif
      if (g_rank == 0) { device->Print(); }
   }

   if (Want("P1")) { TestP1(mesh_file, order, width); }
   if (Want("P1C")) { TestP1C(order, col_h, width, col_fa); }
   if (g_only.count("G4K")) { for (int p : {1, 2}) { TestG4K(mesh_file, p); } }
   if (g_only.count("G4W"))
   {
      G4Opt g4;
      g4.device = dev_name;
      g4.write = g4_write;
      g4.compare = g4_compare;
      g4.tol = g4_tol;
      g4.tol_v = g4_tol_v;
      g4.only = g4_only;
      g4.band = (g4_band != 0);
      MFEM_VERIFY((g4.device == "cpu") != g4.write.empty() && g4.write.empty() != g4.compare.empty(),
                  "U-G4 W: --device cpu --write PREFIX, or --device debug|cuda --compare PREFIX");
      TestG4W(mesh_file, g4, 3e3, width);
   }
   if (Want("P2"))
   {
#ifdef MFEM_USE_PETSC
      const cycle::CycleConfig c0 = cycle::CycleConfig::Load(deck);
      MFEMInitializePetsc(&argc, &argv, c0.qd.petsc_options.c_str(), NULL);
      TestP2(mpi, deck, n_sync);
      MFEMFinalizePetsc();
#else
      if (g_rank == 0) { std::cerr << "U-P2 needs PETSc.\n"; }
      num_failed++;
#endif
   }
   if (g_rank == 0)
   {
      std::cout << (num_tests - num_failed) << "/" << num_tests << " tests passed" << std::endl;
   }
   return num_failed == 0 ? 0 : 1;
}
