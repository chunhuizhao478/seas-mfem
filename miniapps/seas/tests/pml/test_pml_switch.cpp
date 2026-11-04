// tests/pml/test_pml_switch.cpp - PML plan rev2, P2: the switch into the
// dynamic phase with the band (DynPhase with [cycle.pml], Sections 3.1 and
// 7.3), on the simple-case deck D0 (straight fault, +-240 km box, walls
// y = +-240 km attribute 5 loaded, sides x = +-240 km attribute 1 natural).
//
// The deck is run with [cycle.absorbing] off, Kelvin-Voigt off, and the band
// on every outer face ({5, 1}, the P1 checks' walls), wall speed "auto".  One
// QD balance at t_s with a smooth slip on the fault is handed to DynPhase:
//   loading    W's wall data at tau = 0 equals the QD regime's at t_s (the
//              same D3 wall factor in both), to rounding
//   rest       A^0 = W.Mult(u_s, 0) on the band elements' DOFs is zero to the
//              QD solve's residual: max |A^0| there / a_K <= 1e-9, a_K =
//              max |M^-1 K_p u_s|.  Control: without b_band (the static
//              band load) the same DOFs see > 1e-3 a_K
//   dt         DynPhase's step is the plain CFL step / sqrt(kappa_max)
// P3: --as-is runs a P3 deck's own band (the Mac, test-pml-decks);
// --parse-only only reads it.  --kind sponge (the comparator: no QD band, wall factor 1, no b_band,
// the plain step) and --kv-eta 0.1 (Kelvin-Voigt as the simple-case decks).
//   steps      n dynamic steps finite; the band DOFs' max |U - u_s| printed
//
// Usage (from miniapps/seas; PETSc needed):
//   mpirun -np N ./seas_test_pml_switch [--deck simple_cases/configs/D0_straight.toml]
//          [--dyn-steps 200]

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"
#include "../../common/mpi_context.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

namespace
{
int fails = 0;
void Check(MPIContext &mpi, bool ok, const char *what, real_t v)
{
   if (mpi.IsRoot()) { std::printf("  [%s] %s  (%.4e)\n", ok ? "PASS" : "FAIL", what, v); }
   if (!ok) { fails++; }
}
real_t GMax(MPIContext &mpi, real_t v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, mpi.GetComm());
   return v;
}
}  // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "simple_cases/configs/D0_straight.toml";
   int n_dyn = 200;
   const char *kind = "cfs";
   real_t kv = 0.0;
   bool as_is = false, parse_only = false;
   OptionsParser args(argc, argv);
   args.AddOption(&as_is, "-as", "--as-is", "-no-as", "--no-as-is",
                  "Use the deck's own [cycle.pml], [cycle.absorbing] and kv_eta (the P3 decks).");
   args.AddOption(&parse_only, "-po", "--parse-only", "-no-po", "--no-parse-only",
                  "Only read the deck and print its [cycle.pml] (no solve).");
   args.AddOption(&kind, "-k", "--kind", "[cycle.pml] kind: cfs or sponge.");
   args.AddOption(&kv, "-kv", "--kv-eta", "Kelvin-Voigt eta_bar of the dynamic phase (0: off).");
   args.AddOption(&deck, "-d", "--deck", "2-D switching deck with loaded walls (attribute 5).");
   args.AddOption(&n_dyn, "-n", "--dyn-steps", "Dynamic steps after the switch.");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.IsRoot()) { args.PrintUsage(std::cout); }
      return 1;
   }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_pml_switch needs PETSc.\n"; }
   return 2;
#else
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   MFEM_VERIFY(cfg.switching, deck << ": a switching deck is needed");
   if (as_is)
   {
      MFEM_VERIFY(cfg.pml.On(), deck << ": --as-is needs [cycle.pml] on");
      if (mpi.IsRoot())
      {
         std::printf("  deck %s: [cycle.pml] kind %s, walls %zu attrs, L %.0f m, N %d, R0 %g, "
                     "f_alpha %g Hz, c %.2f m/s, wall factor %s; absorbing \"%s\", kv_eta %g, "
                     "max_dynamic_phases %d\n", deck, cfg.pml.kind.c_str(),
                     cfg.pml.wall_attrs.size(), cfg.pml.thickness, cfg.pml.power, cfg.pml.R0,
                     cfg.pml.f_alpha, cfg.pml.c_ref,
                     cfg.pml.AutoWallSpeed() ? "auto" : std::to_string(cfg.pml.wall_speed_factor).c_str(),
                     cfg.absorbing.faces.c_str(), cfg.dyn.kv_eta,
                     cfg.switch_cfg.max_dynamic_phases);
      }
      if (parse_only) { return 0; }
   }
   else
   {
   cfg.absorbing.faces = "off";
   cfg.dyn.kv_eta = kv;
   cfg.pml.on = true;
   cfg.pml.kind = kind;
   cfg.pml.wall_attrs = {5, 1};
   cfg.pml.thickness = 100e3;
   cfg.pml.power = 2;
   cfg.pml.R0 = 1e-3;
   cfg.pml.f_alpha = 0.005;
   {
      const auto &m = cfg.shared.material_fallback;
      cfg.pml.c_ref = std::sqrt((m.lambda + 2.0 * m.mu) / m.rho);
   }
   cfg.pml.wall_speed_factor = cfg.pml.Sponge() ? 1.0 : 0.0;   // "auto" for the CFS-PML
   }
   const bool sponge = cfg.pml.Sponge();
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   {
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      ParMesh pmesh(mpi.GetComm(), *smesh);
      smesh.reset();
      if (mpi.IsRoot())
      {
         std::cout << "test_pml_switch: " << deck << ", p = " << cfg.shared.mesh.order
                   << ", ranks " << mpi.Size() << ", kind " << cfg.pml.kind << ", kv_eta "
                   << kv << "\n";
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

      // the QD balance at t_s: a smooth slip on the fault (both kernel slots;
      // 2-D carries one), psi0
      const real_t ts = 1.0e8;
      const int N = qd.Table().Size();
      Vector sk(2 * N), psi(N);
      for (int i = 0; i < N; i++)
      {
         const real_t x = qd.Table().coords(3 * i) / 15e3;
         const real_t w = std::max(0.0, 1.0 - x * x);
         sk(2 * i) = sk(2 * i + 1) = 0.05 * w;
         psi(i) = pre.psi0(i);
      }
      qd.SetFaultState(sk, psi);
      const real_t Vmax = qd.Solve(ts);
      Vector skq, psiq, sn;
      qd.Slip(skq);
      qd.Psi(psiq);
      qd.SigmaN(sn);
      const Vector Vqd = qd.SlipRate();
      const Vector us = qd.Displacement();
      if (mpi.IsRoot())
      {
         std::printf("  QD balance at t_s = %.3e s: max V %.3e m/s, wall factor %.4f\n", ts,
                     Vmax, qd.PmlWallSpeedFactor());
      }

      cycle::DynPhase dyn(pmesh, cfg, qd.Table(), par, pre, 1.0, nullptr,
                          qd.Domain().FrameField(), qd.SubstrateShape(), qd.PmlBandPtr(),
                          qd.PmlWallSpeedFactor());
      Vector d, V, psin;
      cycle::bridge::ToNitsche(dyn.Frames(), skq, Vqd, psiq, d, V, psin);
      dyn.Enter(us, d, V, psin, sn, ts);
      NitscheWaveOperator &W = dyn.Op();

      // loading: a point on the wall y = +-ymax
      {
         real_t ymax = 0.0;
         for (int v = 0; v < pmesh.GetNV(); v++) { ymax = std::max(ymax, pmesh.GetVertex(v)[1]); }
         ymax = GMax(mpi, ymax);
         Vector x(2), gq(2), gw(2);
         x(0) = 12.3e3;
         x(1) = ymax;
         qd.Domain().EvalDirichlet(5, x, ts, gq);
         W.Elasticity().EvalDirichlet(5, x, 0.0, gw);
         const real_t rel = std::hypot(gq(0) - gw(0), gq(1) - gw(1)) / std::hypot(gq(0), gq(1));
         Check(mpi, rel < 1e-14, "loading: W's wall data at tau = 0 = the QD regime's at t_s", rel);
      }

      // rest: A^0 on the band elements' DOFs
      const PmlBand &band = *qd.PmlBandPtr();
      const ParFiniteElementSpace &fes = W.Elasticity().GetFESpace();
      std::vector<char> inb(fes.GetVSize(), 0);
      {
         Array<int> vd;
         for (int e = 0; e < pmesh.GetNE(); e++)
         {
            if (!PmlBandElement(pmesh, e, band)) { continue; }
            fes.GetElementVDofs(e, vd);
            for (int k : vd) { inb[k] = 1; }
         }
      }
      Vector zero(us.Size()), A0, Ku, MKu;
      zero = 0.0;
      W.Mult(us, zero, A0);
      W.Elasticity().ApplyStiffness(us, Ku);
      W.Mass().Apply(Ku, MKu);
      const real_t aK = GMax(mpi, MKu.Normlinf());
      real_t a_band = 0.0;
      for (int i = 0; i < A0.Size(); i++) { if (inb[i]) { a_band = std::max(a_band, std::abs(A0(i))); } }
      a_band = GMax(mpi, a_band);
      // control: the same without b_band (the CFS-PML only; the sponge has none)
      real_t c_band = 0.0;
      if (!sponge)
      {
         Vector bb(W.Pml()->BandLoad()), Mbb;
         W.Mass().Apply(bb, Mbb);
         for (int i = 0; i < A0.Size(); i++)
         { if (inb[i]) { c_band = std::max(c_band, std::abs(A0(i) + Mbb(i))); } }
      }
      c_band = GMax(mpi, c_band);
      if (mpi.IsRoot())
      {
         std::printf("  a_K = max |M^-1 K_p u_s| = %.4e m/s^2; band DOFs: max |A^0| %.4e, "
                     "without b_band %.4e\n", aK, a_band, c_band);
      }
      Check(mpi, a_band <= 1e-9 * aK, "rest: A^0 on the band DOFs <= 1e-9 a_K", a_band / aK);
      if (!sponge)
      {
         Check(mpi, c_band > 1e-3 * aK, "rest control: without b_band > 1e-3 a_K", c_band / aK);
      }

      // dt
      {
         const real_t ceig = NitscheStableConstant(pmesh.Dimension(), cfg.shared.mesh.order);
         const real_t plain = ComputeMaxDt(W.Sample(), pmesh, cfg.shared.mesh.order,
                                           cfg.dyn.cfl * ceig);
         const real_t want = plain / std::sqrt(dyn.PmlKappaMax());
         if (mpi.IsRoot())
         {
            std::printf("  dt %.4e s = plain %.4e s / sqrt(kappa_max %.3f)\n", dyn.Dt(), plain,
                        dyn.PmlKappaMax());
         }
         Check(mpi, std::abs(dyn.Dt() / want - 1.0) < 1e-14, sponge ? "dt: the plain CFL step "
               "(the sponge leaves K alone, kappa_max = 1)" : "dt: plain CFL / sqrt(kappa_max)",
               std::abs(dyn.Dt() / want - 1.0));
      }

      // steps
      for (int k = 0; k < n_dyn; k++) { dyn.Step(); }
      const Vector &U = dyn.U();
      real_t du_b = 0.0, du = 0.0;
      bool fin = true;
      for (int i = 0; i < U.Size(); i++)
      {
         if (!std::isfinite(U(i))) { fin = false; }
         du = std::max(du, std::abs(U(i) - us(i)));
         if (inb[i]) { du_b = std::max(du_b, std::abs(U(i) - us(i))); }
      }
      du = GMax(mpi, du);
      du_b = GMax(mpi, du_b);
      int fin_all = fin ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &fin_all, 1, MPI_INT, MPI_MIN, mpi.GetComm());
      const real_t vmax_dyn = GMax(mpi, W.SlipRate().Normlinf());
      if (mpi.IsRoot())
      {
         std::printf("  %d steps (tau = %.3f s): max |U - u_s| %.4e m, on the band DOFs %.4e m, "
                     "max V %.4e m/s\n", n_dyn, W.GetTime(), du, du_b, vmax_dyn);
      }
      Check(mpi, fin_all == 1, "steps: finite", 0.0);
   }
   if (mpi.IsRoot()) { std::cout << (fails ? "FAILED: " : "ALL PASS: ") << fails << "\n"; }
   MFEMFinalizePetsc();
   return fails ? 1 : 0;
#endif
}
