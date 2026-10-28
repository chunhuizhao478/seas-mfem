// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_2d_absorbing.cpp — 2D in-plane design rev 3, section 6.4: the dynamic
// regime's absorbing outer faces ([cycle.absorbing], cycle/absorbing_faces.hpp,
// hook H3 of the fully-dynamic-cycle P5) on a 2-D BP3 mesh, before the BP3-FD
// runs (Chunhui Zhao, 2026-10-27: "merge the absorbing boundary condition for
// 2d, verify before it runs BP3-FD").  It follows the 3-D BP5 test
// tests/cycle/test_cycle_absorbing.cpp check for check.
//
// The deck (default config/fd_cycle/bp3_fd_d60_thrust.toml) runs on its 1 km
// mesh (bp3_dNN_1km.msh) at p = 1.  BP3's outer boundary: the surface y = 0
// and the bottom y = y_min are traction free (attribute 1), the two sides
// (parallel to the fault: bp3/mesh/bp3.geo) are clamped (attribute 5) and carry
// the plate motion, as does the creep
// segment below W_f = 40 km down-dip (attribute 5, interior).  One QD balance
// at t_s = 3e9 s is handed over (DynPhase::Enter) to four dynamic regimes
// built in turn from the deck: reflecting (faces = "off", the reference),
// "outer", "bottom" and "all".  Accelerations are measured against a_K = max
// |M^-1 K u_s| of the reflecting W.
//
//   loaded_state       The fault state: dip slip (the kernel dip slot; the
//                      strike slot of a 2-D BP3 run is 0) S (1 - w / 2) +
//                      0.1 m w f_1(x), w = 1 - (x_d / W_f)^2 above W_f and 0
//                      below (x_d: the distance from the surface trace), so the
//                      slip is S where the fault meets the creep segment; psi0 +
//                      0.01 f_2(x).  S = +-Vp t_s: the sign for which a uniform
//                      slip S with psi0 (the two sides moved with the walls)
//                      leaves the balance's max V as at t = 0 to 1e-6; exactly
//                      one sign must.
//   retag              "outer": every re-tagged traction-free face lies on the
//                      bottom (y = y_min) and every kept one on the surface
//                      (y = 0); every wall face on a side (the line from a
//                      surface corner to the bottom corner below it); the
//                      creep segment keeps attribute 5; the QD mesh's boundary
//                      attributes unchanged; W absorbs on both new attributes.
//   frozen_load        b = K_abs u_s - G d_s - F_c(t_s) is not zero on the
//                      wall elements and zero elsewhere; its net force in x and
//                      in y is zero up to the QD solve residual: |sum b_c| /
//                      sum |b_c| <= 1e-6.
//   switch_equilibrium The first acceleration A^0 of the "outer" W, with b:
//                      off the wall elements equal to the reflecting W's, on
//                      them 0, both to 1e-12 a_K.  Control: without b, A^0 on
//                      the wall elements differs from the reflecting one by more
//                      than 0.5 max |M^-1 b| (the static load the walls drop).
//   first_steps        n steps of each: slip, slip rate and psi as in the
//                      reflecting run to 1e-12 (no wave has reached a boundary).
//   bottom             "bottom": the bottom faces only, no frozen load, A^0
//                      equal to the reflecting one to 1e-13 a_K, the slip rate
//                      and U after n steps to 1e-12.
//   all                "all" (every outer face absorbs, no free surface; the
//                      unit test is test_absorbing_all.cpp): the surface faces
//                      absorb too, the walls are "outer"'s; the frozen load
//                      equals "outer"'s to 1e-13 of its size; A^0 off the wall
//                      elements equal to the reflecting one and 0 on them, to
//                      1e-12 a_K.  The fault meets the absorbing surface, so
//                      the first steps differ from the reflecting run (printed).
//
// Usage (from miniapps/seas; PETSc needed):
//   mpirun -np N ./seas_test_2d_absorbing [--deck config/fd_cycle/bp3_fd_d60_thrust.toml]
//          [--order 1] [--dyn-steps 20]
// Exit code 0 = PASSED.

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

struct Checker
{
   const MPIContext &mpi;
   int failures = 0;

   void Report(const std::string &name, real_t m, real_t tol, const std::string &what)
   {
      const bool ok = (m <= tol);
      if (!ok) { failures++; }
      if (mpi.IsRoot())
      {
         std::cout << "  " << (ok ? "ok  " : "FAIL") << " " << std::left << std::setw(46)
                   << name << std::right << " " << what << " " << std::scientific
                   << std::setprecision(3) << m << " (tolerance " << tol << ")\n";
      }
   }
   void True(const std::string &name, bool ok_local, const std::string &detail)
   {
      const bool ok = mpi.GlobalMinInt(ok_local ? 1 : 0) == 1;
      if (!ok) { failures++; }
      if (mpi.IsRoot())
      {
         std::cout << "  " << (ok ? "ok  " : "FAIL") << " " << std::left << std::setw(46)
                   << name << std::right << " " << detail << "\n";
      }
   }
   real_t MaxAbs(const Vector &x) const
   {
      real_t m = 0.0;
      for (int i = 0; i < x.Size(); i++) { m = std::max(m, std::abs(x(i))); }
      return mpi.GlobalMax(m);
   }
   real_t DiffInf(const Vector &x, const Vector &y) const
   {
      MFEM_VERIFY(x.Size() == y.Size(), "DiffInf: sizes " << x.Size() << " vs " << y.Size());
      real_t d = 0.0;
      for (int i = 0; i < x.Size(); i++) { d = std::max(d, std::abs(x(i) - y(i))); }
      return mpi.GlobalMax(d);
   }
   real_t MaxAbsOn(const Vector &x, const std::vector<char> &mask, char want) const
   {
      real_t m = 0.0;
      for (int i = 0; i < x.Size(); i++)
      {
         if (mask[i] == want) { m = std::max(m, std::abs(x(i))); }
      }
      return mpi.GlobalMax(m);
   }
   real_t DiffInfOn(const Vector &x, const Vector &y, const std::vector<char> &mask,
                    char want) const
   {
      MFEM_VERIFY(x.Size() == y.Size(), "DiffInfOn: sizes " << x.Size() << " vs " << y.Size());
      real_t d = 0.0;
      for (int i = 0; i < x.Size(); i++)
      {
         if (mask[i] == want) { d = std::max(d, std::abs(x(i) - y(i))); }
      }
      return mpi.GlobalMax(d);
   }
   real_t RelInf(const Vector &x, const Vector &y) const
   {
      const real_t d = DiffInf(x, y), s = MaxAbs(y);
      return (s > 0.0) ? d / s : d;
   }
};

/// A smooth deterministic field of position, partition independent.
real_t Field(const real_t *x, int c)
{
   return std::sin(1.3e-4 * x[0] + 0.7 * c) * std::cos(2.1e-4 * x[1] - 0.3 * c)
          + 0.25 * std::sin(3.7e-5 * x[1] + 1.1 * c);
}

/// The loaded fault state (see loaded_state): dip slip S (1 - f w) + dip w f_1,
/// strike slip 0, psi0 + dpsi f_2; w = 1 - (x_d / W_f)^2 above W_f, 0 below.
void LoadedFaultState(const cycle::QDPhase &qd, const cycle::Prestress &pre, real_t S,
                      real_t f, real_t dip, real_t dpsi, Vector &s, Vector &psi)
{
   const real_t Wf = 40.0e3;
   const Vector &X = qd.Table().coords;
   const int N = qd.Table().Size();
   s.SetSize(2 * N);
   psi.SetSize(N);
   for (int i = 0; i < N; i++)
   {
      const real_t xd = std::hypot(X(3 * i), X(3 * i + 1)) / Wf;
      const real_t w = std::max(0.0, 1.0 - xd * xd);
      s(2 * i) = S * (1.0 - f * w) + dip * w * Field(&X(3 * i), 1);   // kernel pairs: (dip, strike)
      s(2 * i + 1) = 0.0;
      psi(i) = pre.psi0(i) + dpsi * Field(&X(3 * i), 2);
   }
}

struct Balance
{
   Vector u, sk, Vqd, psi, sn;
   real_t t = 0.0, Vmax = 0.0;
};

struct Handover
{
   Vector u, d, V, psi, sn;
   real_t t = 0.0;
};

Handover ToNitsche(const Balance &b, const cycle::DynPhase &dyn)
{
   Handover h;
   cycle::bridge::ToNitsche(dyn.Frames(), b.sk, b.Vqd, b.psi, h.d, h.V, h.psi);
   h.u = b.u;
   h.sn = b.sn;
   h.t = b.t;
   return h;
}

struct Run
{
   Vector A0, U, slip, V, psi;
   Vector Minv_b;      ///< M^-1 b ("outer" only)
   real_t tau = 0.0;   ///< W's clock after the steps (time since Enter)
};

/// Enter, A^0 = W.Mult(u_s, 0), n steps.
Run Phase(cycle::DynPhase &dyn, const Balance &b, int n)
{
   Run r;
   const Handover h = ToNitsche(b, dyn);
   dyn.Enter(h.u, h.d, h.V, h.psi, h.sn, h.t);
   Vector zero(h.u.Size());
   zero = 0.0;
   dyn.Op().Mult(h.u, zero, r.A0);
   if (dyn.FrozenWallLoad().Size() > 0)
   {
      dyn.Op().Mass().Apply(dyn.FrozenWallLoad(), r.Minv_b);
   }
   for (int k = 0; k < n; k++) { dyn.Step(); }
   r.U = dyn.U();
   r.slip = dyn.Slip();
   r.V = dyn.Op().SlipRate();
   r.psi = dyn.Psi();
   r.tau = dyn.Op().GetTime();
   return r;
}

/// 1 on the local DOFs of the elements with a boundary face of attribute
/// `attr` (DG: the element's own DOFs), 0 elsewhere.
std::vector<char> ElementDofMask(NitscheWaveOperator &W, int attr)
{
   ParMesh &m = W.GetMesh();
   const ParFiniteElementSpace &fes = W.Elasticity().GetFESpace();
   std::vector<char> on(fes.GetVSize(), 0);
   Array<int> vd;
   for (int be = 0; be < m.GetNBE(); be++)
   {
      if (m.GetBdrAttribute(be) != attr) { continue; }
      int e = -1, info = 0;
      m.GetBdrElementAdjacentElement(be, e, info);
      fes.GetElementVDofs(e, vd);
      for (int i = 0; i < vd.Size(); i++) { on[vd[i] >= 0 ? vd[i] : -1 - vd[i]] = 1; }
   }
   return on;
}

std::vector<int> BdrAttributes(const ParMesh &m)
{
   std::vector<int> a(m.GetNBE());
   for (int be = 0; be < m.GetNBE(); be++) { a[be] = m.GetBdrAttribute(be); }
   return a;
}

/// The midpoint of boundary element `be` (a segment).
void BdrCenter(const ParMesh &m, int be, real_t c[2])
{
   Array<int> v;
   m.GetBdrElementVertices(be, v);
   c[0] = c[1] = 0.0;
   for (int k = 0; k < v.Size(); k++)
   {
      for (int d = 0; d < 2; d++) { c[d] += m.GetVertex(v[k])[d] / v.Size(); }
   }
}

/// The 1 km mesh of a deck's 25 m mesh: bp3_dNN_25m.msh -> bp3_dNN_1km.msh.
std::string CoarseMesh(const std::string &path)
{
   const std::string fine = "_25m.msh";
   MFEM_VERIFY(path.size() > fine.size() && path.compare(path.size() - fine.size(),
                                                          fine.size(), fine) == 0,
               "test_2d_absorbing: the deck's mesh " << path << " is not a _25m.msh");
   return path.substr(0, path.size() - fine.size()) + "_1km.msh";
}

} // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "config/fd_cycle/bp3_fd_d60_thrust.toml";
   int order = 1, n_dyn = 20;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "BP3-FD deck (run on its 1 km mesh).");
   args.AddOption(&order, "-o", "--order", "Polynomial order.");
   args.AddOption(&n_dyn, "-nd", "--dyn-steps", "Dynamic steps per run.");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.IsRoot()) { args.PrintUsage(std::cout); }
      return 1;
   }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_2d_absorbing needs PETSc (the QD TS).\n"; }
   return 2;
#else
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   MFEM_VERIFY(cfg.switching, deck << ": [cycle] switching is false; the test needs the "
               "dynamic regime's sections");
   MFEM_VERIFY(!cfg.lts.On(), deck << ": the test compares global stepping runs");
   MFEM_VERIFY(n_dyn >= 1, "--dyn-steps must be >= 1");
   cfg.shared.mesh.path = CoarseMesh(cfg.shared.mesh.path);
   cfg.shared.mesh.order = order;
   cfg.absorbing.faces = "off";
   cycle::CycleConfig cfg_outer = cfg, cfg_bottom = cfg, cfg_all = cfg;
   cfg_outer.absorbing.faces = "outer";
   cfg_bottom.absorbing.faces = "bottom";
   cfg_all.absorbing.faces = "all";
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   int failures = 0;
   {
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      cfg.CheckDimension(smesh->Dimension());
      MFEM_VERIFY(smesh->Dimension() == 2, cfg.shared.mesh.path << " is not a 2-D mesh");
      ParMesh pmesh(mpi.GetComm(), *smesh);
      smesh.reset();
      if (mpi.IsRoot())
      {
         std::cout << "test_2d_absorbing: " << deck << ", " << cfg.shared.mesh.path << " x "
                   << cfg.mesh_scale << ", p = " << order << ", ranks " << mpi.Size()
                   << ", commit " << SEAS_CYCLE_COMMIT << "\n";
      }
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
      qd.Setup(par, pre);
      auto frames = qd.Domain().FrameField();   // one orientation for both regimes
      Checker c{mpi};

      // loaded_state: the slip's sign and size, then the QD balance at t_s.
      const real_t ts = 3.0e9;
      real_t Sp = 0.0;
      {
         Vector s0, psi0;
         LoadedFaultState(qd, pre, 0.0, 0.0, 0.0, 0.0, s0, psi0);
         qd.SetFaultState(s0, psi0);
         const real_t V0 = qd.Solve(0.0);
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(6) << "\nloaded_state: max V "
                      << V0 << " m/s at t = 0 (slip 0)\n";
         }
         int matches = 0;
         for (const real_t sg : {1.0, -1.0})
         {
            const real_t S = sg * cfg.plate_rate * ts;
            LoadedFaultState(qd, pre, S, 0.0, 0.0, 0.0, s0, psi0);
            qd.SetFaultState(s0, psi0);
            const real_t Vu = qd.Solve(ts), rel = std::abs(Vu / V0 - 1.0);
            if (mpi.IsRoot())
            {
               std::cout << "  uniform dip slip " << S << " m at t_s: max V " << Vu
                         << " m/s, |ratio - 1| " << rel << "\n";
            }
            if (rel <= 1e-6) { matches++; Sp = S; }
         }
         c.True("loaded_state: one sign of uniform slip is stress free", matches == 1,
                "|max V ratio - 1| <= 1e-6 for exactly one sign");
         if (matches != 1) { Sp = cfg.plate_rate * ts; }
      }
      Balance bal;
      {
         Vector s0, psi0;
         LoadedFaultState(qd, pre, Sp, 0.5, 0.1, 0.01, s0, psi0);
         qd.SetFaultState(s0, psi0);
         bal.Vmax = qd.Solve(ts);
         qd.Slip(bal.sk);
         qd.Psi(bal.psi);
         qd.SigmaN(bal.sn);
         bal.Vqd = qd.SlipRate();
         bal.u = qd.Displacement();
         bal.t = ts;
      }
      const real_t mus = c.MaxAbs(bal.u);
      if (mpi.IsRoot())
      {
         std::cout << std::scientific << std::setprecision(3) << "\nQD balance at t_s = " << ts
                   << " s: sides at +-" << 0.5 * cfg.plate_rate * ts << " m, max V "
                   << bal.Vmax << " m/s, max |u_s| " << mus << " m\n";
      }
      const std::vector<int> qd_attrs = BdrAttributes(pmesh);

      // the box: y_min (the bottom), and its corners on the surface (xs) and on
      // the bottom (xb), left [0] and right [1], over every rank; the sides run
      // from (xs[k], 0) to (xb[k], y_min)
      real_t ymin = 0.0, xabs = 0.0;
      for (int v = 0; v < pmesh.GetNV(); v++)
      {
         ymin = std::min(ymin, pmesh.GetVertex(v)[1]);
         xabs = std::max(xabs, std::abs(pmesh.GetVertex(v)[0]));
      }
      ymin = mpi.GlobalMin(ymin);
      xabs = mpi.GlobalMax(xabs);
      const real_t tol = 1e-9 * std::max(-ymin, xabs);
      real_t xs[2] = {xabs, -xabs}, xb[2] = {xabs, -xabs};
      for (int v = 0; v < pmesh.GetNV(); v++)
      {
         const real_t *X = pmesh.GetVertex(v);
         real_t *cx = (std::abs(X[1]) <= tol) ? xs
                      : (std::abs(X[1] - ymin) <= tol) ? xb : nullptr;
         if (cx) { cx[0] = std::min(cx[0], X[0]); cx[1] = std::max(cx[1], X[0]); }
      }
      for (int k = 0; k < 2; k++)
      {
         xs[k] = (k == 0) ? mpi.GlobalMin(xs[k]) : mpi.GlobalMax(xs[k]);
         xb[k] = (k == 0) ? mpi.GlobalMin(xb[k]) : mpi.GlobalMax(xb[k]);
      }
      auto on_side = [&](const real_t *x)
      {
         const real_t s = x[1] / ymin;   // 0 at the surface, 1 at the bottom
         for (int k = 0; k < 2; k++)
         {
            if (std::abs(x[0] - (xs[k] + s * (xb[k] - xs[k]))) <= tol) { return true; }
         }
         return false;
      };

      // ---- reflecting (the reference) ----
      Run refl;
      real_t aK = 0.0;
      {
         cycle::DynPhase dyn(pmesh, cfg, qd.Table(), par, pre, 1.0, nullptr, frames);
         refl = Phase(dyn, bal, n_dyn);
         c.True("off: no frozen load", dyn.FrozenWallLoad().Size() == 0, "");
         Vector Ku, MKu;
         dyn.Op().StiffnessMult(bal.u, Ku);
         dyn.Op().Mass().Apply(Ku, MKu);
         aK = c.MaxAbs(MKu);
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(3)
                      << "  reflecting W: a_K = max |M^-1 K u_s| = " << aK << " m/s^2\n";
         }
      }

      // ---- "outer" ----
      if (mpi.IsRoot()) { std::cout << "\nretag, frozen_load, switch_equilibrium (\"outer\")\n"; }
      Run outer;
      cycle::OuterFaceTags tags_outer;
      Vector b_outer;
      {
         cycle::DynPhase dyn(pmesh, cfg_outer, qd.Table(), par, pre, 1.0, nullptr, frames);
         const cycle::OuterFaceTags &t = dyn.AbsorbingFaces();
         tags_outer = t;
         if (mpi.IsRoot())
         {
            std::cout << "  re-tagged: " << t.natural << " traction-free faces -> attribute "
                      << t.natural_to << ", " << t.dirichlet << " wall faces -> attribute "
                      << t.dirichlet_to << "; " << t.free_surface << " free-surface faces kept\n";
         }
         c.True("retag: faces in each class", t.natural > 0 && t.dirichlet > 0
                && t.free_surface > 0, "bottom, side and surface faces all found");
         // where each class lies, on W's mesh
         ParMesh &wm = dyn.Op().GetMesh();
         long long n[5] = {0, 0, 0, 0, 0};   // bottom, sides, surface, creep, misplaced
         for (int be = 0; be < wm.GetNBE(); be++)
         {
            const int a = wm.GetBdrAttribute(be);
            const bool ext = wm.GetFaceInformation(wm.GetBdrElementFaceIndex(be)).IsBoundary();
            real_t x[2];
            BdrCenter(wm, be, x);
            if (a == t.natural_to)
            {
               if (ext && std::abs(x[1] - ymin) <= tol) { n[0]++; } else { n[4]++; }
            }
            else if (a == t.dirichlet_to)
            {
               if (ext && on_side(x)) { n[1]++; } else { n[4]++; }
            }
            else if (a == 1)
            {
               if (ext && std::abs(x[1]) <= tol) { n[2]++; } else { n[4]++; }
            }
            else if (a == 5)
            {
               if (!ext) { n[3]++; } else { n[4]++; }   // the creep segment, interior
            }
         }
         MPI_Allreduce(MPI_IN_PLACE, n, 5, MPI_LONG_LONG, MPI_SUM, mpi.GetComm());
         if (mpi.IsRoot())
         {
            std::cout << "  on W's mesh: " << n[0] << " absorbing faces at y = " << ymin
                      << " m, " << n[1] << " absorbing walls on the sides (surface corners x = "
                      << xs[0] << ", " << xs[1] << " m, bottom corners x = " << xb[0] << ", "
                      << xb[1] << " m), "
                      << n[2] << " traction-free faces at y = 0, " << n[3]
                      << " creep-segment faces (interior, attribute 5), " << n[4]
                      << " elsewhere\n";
         }
         c.True("retag: bottom, sides, surface where they belong",
                n[4] == 0 && n[0] == t.natural && n[1] == t.dirichlet && n[2] == t.free_surface
                && n[3] > 0, "no face out of place");
         c.True("retag: the QD mesh is unchanged", BdrAttributes(pmesh) == qd_attrs,
                "every boundary attribute of the QD regime's mesh");
         const std::set<int> &abs = dyn.Boundaries().absorbing_attrs;
         c.True("retag: W absorbs on both new attributes",
                abs.count(t.natural_to) == 1 && abs.count(t.dirichlet_to) == 1
                && abs.size() == 2, "");

         outer = Phase(dyn, bal, 0);   // A^0 with b, no step yet

         // frozen_load
         NitscheWaveOperator &W = dyn.Op();
         const std::vector<char> wall = ElementDofMask(W, t.dirichlet_to);
         const Vector &b = dyn.FrozenWallLoad();
         b_outer = b;
         const real_t bw = c.MaxAbsOn(b, wall, 1), bo = c.MaxAbsOn(b, wall, 0);
         const ParFiniteElementSpace &fes = W.Elasticity().GetFESpace();
         const int vdim = fes.GetVDim();
         c.True("frozen_load: two displacement components", vdim == 2, "vdim = "
                + std::to_string(vdim));
         real_t fsum[4] = {0.0, 0.0, 0.0, 0.0};   // sum b_c, sum |b_c|
         for (int i = 0; i < fes.GetNDofs(); i++)
         {
            for (int d = 0; d < 2; d++)
            {
               const real_t v = b(fes.DofToVDof(i, d));
               fsum[d] += v;
               fsum[2 + d] += std::abs(v);
            }
         }
         MPI_Allreduce(MPI_IN_PLACE, fsum, 4, MPITypeMap<real_t>::mpi_type, MPI_SUM,
                       mpi.GetComm());
         real_t net = 0.0;
         for (int d = 0; d < 2; d++)
         {
            net = std::max(net, fsum[2 + d] > 0.0 ? std::abs(fsum[d]) / fsum[2 + d] : 0.0);
         }
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(3)
                      << "  frozen load: max |b| on the wall elements " << bw << " N/m, elsewhere "
                      << bo << " N/m; net force (x, y) " << fsum[0] << ", " << fsum[1]
                      << " N/m against sum |b_c| " << fsum[2] << ", " << fsum[3] << " N/m\n";
         }
         c.True("frozen_load: the walls carry a load", bw > 0.0, "max |b| on the wall elements > 0");
         c.Report("frozen_load: b off the wall elements", bo, 0.0, "max |b| (N/m)");
         c.Report("frozen_load: net force of b", net, 1e-6, "max_c |sum b_c| / sum |b_c|");

         // switch_equilibrium
         const real_t mb = c.MaxAbs(outer.Minv_b);
         c.Report("switch_equilibrium: A0 off the wall elements",
                  c.DiffInfOn(outer.A0, refl.A0, wall, 0) / aK, 1e-12,
                  "max |A0_outer - A0_refl| / a_K");
         c.Report("switch_equilibrium: A0 on the wall elements",
                  c.MaxAbsOn(outer.A0, wall, 1) / aK, 1e-12, "max |A0_outer| / a_K");
         const real_t rw = c.MaxAbsOn(refl.A0, wall, 1) / aK;
         Vector none, a_nob, zero(bal.u.Size());
         zero = 0.0;
         W.SetConstantLoad(none);
         W.Mult(bal.u, zero, a_nob);
         const real_t ctl = c.DiffInfOn(a_nob, refl.A0, wall, 1) / mb;
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(3)
                      << "  reflecting A0 on the wall elements (M^-1 r, left out of b): " << rw
                      << " a_K; max |M^-1 b| = " << mb / aK << " a_K\n"
                      << "  control, without b: max |A0 - A0_refl| on the wall elements = "
                      << ctl << " max |M^-1 b|\n";
         }
         c.True("control: without b the walls unload", ctl > 0.5,
                "max |A0 - A0_refl| / max |M^-1 b| > 0.5 on the wall elements");

         // first_steps: enter again (b set again), n steps
         outer = Phase(dyn, bal, n_dyn);
      }
      if (mpi.IsRoot()) { std::cout << "\nfirst_steps (" << n_dyn << " steps)\n"; }
      c.Report("first_steps outer: slip", c.RelInf(outer.slip, refl.slip), 1e-12, "max rel diff");
      c.Report("first_steps outer: slip rate", c.RelInf(outer.V, refl.V), 1e-12, "max rel diff");
      c.Report("first_steps outer: psi", c.RelInf(outer.psi, refl.psi), 1e-12, "max rel diff");
      {
         Vector du = refl.U;
         du -= bal.u;
         const real_t dU = c.DiffInf(outer.U, refl.U), mU = c.MaxAbs(du);
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(3) << "  U: max |U_outer - U_refl| "
                      << dU << " m against max |U_refl - u_s| " << mU
                      << " m (printed only: the reflecting walls move Vp tau / 2 = "
                      << 0.5 * cfg.plate_rate * refl.tau << " m in tau = " << refl.tau
                      << " s)\n";
         }
      }

      // ---- "bottom" ----
      if (mpi.IsRoot()) { std::cout << "\nbottom\n"; }
      {
         cycle::DynPhase dyn(pmesh, cfg_bottom, qd.Table(), par, pre, 1.0, nullptr, frames);
         const cycle::OuterFaceTags &t = dyn.AbsorbingFaces();
         c.True("bottom: no wall faces, no frozen load", t.dirichlet == 0 && t.natural > 0,
                std::to_string(t.natural) + " traction-free faces absorb");
         const Run bot = Phase(dyn, bal, n_dyn);
         c.True("bottom: no frozen load after Enter", dyn.FrozenWallLoad().Size() == 0, "");
         c.Report("bottom: A0 vs reflecting", c.DiffInf(bot.A0, refl.A0) / aK, 1e-13,
                  "max |dA0| / a_K");
         c.Report("bottom: slip rate", c.RelInf(bot.V, refl.V), 1e-12, "max rel diff");
         c.Report("bottom: U", c.RelInf(bot.U, refl.U), 1e-12, "max rel diff");
      }

      // ---- "all" ----
      if (mpi.IsRoot()) { std::cout << "\nall\n"; }
      {
         cycle::DynPhase dyn(pmesh, cfg_all, qd.Table(), par, pre, 1.0, nullptr, frames);
         const cycle::OuterFaceTags &t = dyn.AbsorbingFaces();
         if (mpi.IsRoot())
         {
            std::cout << "  re-tagged: " << t.natural << " traction-free faces -> attribute "
                      << t.natural_to << ", " << t.dirichlet << " wall faces -> attribute "
                      << t.dirichlet_to << "; " << t.free_surface << " free-surface faces kept\n";
         }
         c.True("all: the surface absorbs too, the walls as outer",
                t.free_surface == 0 && t.natural == tags_outer.natural + tags_outer.free_surface
                && t.dirichlet == tags_outer.dirichlet,
                std::to_string(t.natural) + " traction-free faces, " + std::to_string(t.dirichlet)
                + " walls");
         const std::set<int> &abs = dyn.Boundaries().absorbing_attrs;
         c.True("all: W absorbs on both new attributes",
                abs.count(t.natural_to) == 1 && abs.count(t.dirichlet_to) == 1
                && abs.size() == 2, "");
         Run all = Phase(dyn, bal, 0);
         NitscheWaveOperator &W = dyn.Op();
         const std::vector<char> wall = ElementDofMask(W, t.dirichlet_to);
         c.Report("all: frozen load as outer's", c.RelInf(dyn.FrozenWallLoad(), b_outer), 1e-13,
                  "max |b_all - b_outer| / max |b_outer|");
         c.Report("all: A0 off the wall elements", c.DiffInfOn(all.A0, refl.A0, wall, 0) / aK,
                  1e-12, "max |A0_all - A0_refl| / a_K");
         c.Report("all: A0 on the wall elements", c.MaxAbsOn(all.A0, wall, 1) / aK, 1e-12,
                  "max |A0_all| / a_K");
         all = Phase(dyn, bal, n_dyn);
         const real_t ds = c.RelInf(all.slip, refl.slip), dv = c.RelInf(all.V, refl.V);
         if (mpi.IsRoot())
         {
            std::cout << std::scientific << std::setprecision(3) << "  after " << n_dyn << " steps (printed only: the fault meets the "
                      << "absorbing surface): slip max rel diff " << ds << ", slip rate " << dv
                      << " against the reflecting run\n";
         }
      }

      failures = c.failures;
      if (mpi.IsRoot())
      {
         std::cout << "\ntest_2d_absorbing: " << (failures ? "FAILED" : "PASSED") << " ("
                   << failures << " failing checks)\n";
      }
   }
   MFEMFinalizePetsc();
   return failures ? 1 : 0;
#endif
}
