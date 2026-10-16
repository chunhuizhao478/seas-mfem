// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_cycle_params_parity.cpp — fully dynamic BP5 cycle, P1 test
// test_cycle_params_parity (plan rev 3).
//
// The cycle driver builds its fault parameters from deck values
// (BuildFaultParams, BuildPrestress, hooks H1a and H1b).  This test checks,
// on every owned fault DOF, that they equal the analytic BP5 parameters that
// seas_bp5_full's gold run used (BP5Params, through the legacy FaultGeometry
// constructor on the same domain operator):
//
//   gold deck (ic "gold"):
//     ref = FaultGeometry(dom, BP5Params(), mpi)        # BP5Params at each DOF
//     a, Dc, V_init (2 components), tau0 (2 components) vs ref: rel. 1e-12
//     b, f0, V0, sigma_n, eta vs BP5Params' scalars:     rel. 1e-12
//     the QD geometry after H1a + H1b holds the same arrays as ref
//     psi0: kernel (PreInit + Init at zero traction) vs BuildPrestress
//   scec_fd deck (ic "scec_fd"):
//     ref with BP5Params V_nuc = 0.03 (SCEC V_i)
//     a, Dc, V_init vs ref; tau0 vs ref's tau0 without eta |V| (eq. 22)
//     psi0: kernel (SCEC psi, eq. 18) vs BuildPrestress
// The nucleation-patch DOFs and those on its edge are counted; both must be
// present (the edge is where the 1 mm inclusion tolerance decides).
//
// Zero traction is the operator's initial traction: zero slip and zero
// boundary displacement at t = 0 (the driver prints the same psi0 check
// after the real initial solve).  No domain solve runs here: the operator
// assembles and factors on its first solve.
//
// Usage (from miniapps/seas):
//   mpirun -np N ./seas_test_cycle_params_parity
//          [--deck config/fd_cycle/bp5_gold.toml]
//          [--deck-fd config/fd_cycle/bp5_scec_fd.toml]
// Exit code 0 = PASSED.

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{

constexpr real_t kTol = 1e-12;

/// Relative difference with a floor, so exact zeros compare exactly.
real_t RelDiff(real_t x, real_t y)
{
   const real_t d = std::abs(x - y);
   const real_t s = std::max(std::abs(x), std::abs(y));
   return (s > 0.0) ? d / s : d;
}

struct Checker
{
   const MPIContext &mpi;
   int failures = 0;

   /// Max relative difference of two equal-sized vectors, over all ranks.
   void Vec(const std::string &name, const Vector &x, const Vector &y)
   {
      int bad_size = (x.Size() != y.Size()) ? 1 : 0;
      bad_size = mpi.GlobalMaxInt(bad_size);
      if (bad_size)
      {
         if (mpi.IsRoot()) { std::cout << "  FAIL " << name << ": size mismatch\n"; }
         failures++;
         return;
      }
      real_t m = 0.0;
      for (int i = 0; i < x.Size(); i++) { m = std::max(m, RelDiff(x(i), y(i))); }
      Report(name, mpi.GlobalMax(m));
   }

   /// Every entry of x against the scalar y.
   void Scalar(const std::string &name, const Vector &x, real_t y)
   {
      real_t m = 0.0;
      for (int i = 0; i < x.Size(); i++) { m = std::max(m, RelDiff(x(i), y)); }
      Report(name, mpi.GlobalMax(m));
   }

   void Report(const std::string &name, real_t m)
   {
      const bool ok = (m <= kTol);
      if (!ok) { failures++; }
      if (mpi.IsRoot())
      {
         std::cout << "  " << (ok ? "ok  " : "FAIL") << " " << std::left
                   << std::setw(34) << name << std::right << " max rel diff "
                   << std::scientific << std::setprecision(3) << m << "\n";
      }
   }
};

/// The kernel's initial psi on the QD geometry after hooks H1a and H1b, as
/// QDPhase::Setup builds the kernel, at zero traction.
Vector KernelPsi0(cycle::QDPhase &qd, const cycle::CycleConfig &cfg,
                  const cycle::FaultParams &par, const cycle::Prestress &pre,
                  MPIContext &mpi)
{
   const BP5Params &seed = qd.Seed();
   Vector dir(2);
   dir(0) = 0.0;
   dir(1) = 1.0;
   qd.Geometry().SetRateStatePerDOF(par.rs, dir);
   qd.Geometry().SetInitialFieldsPerDOF(pre.tau0, par.V_init, pre.sigma0);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = seed.V0;
   fc.f0 = seed.f0;
   fc.b = seed.b;
   fc.Dc = seed.L0;
   DieterichRuinaFriction friction(fc);
   AgingLawPsi aging(seed.b, seed.V0, seed.f0);
   cycle::QDPhase::FaultOp fault(&qd.Geometry(), &friction, &aging, seed, &mpi);
   fault.SetScecPsiInit(!cfg.IsGold());

   Vector state(fault.StateSize());
   Vector traction(fault.TractionSize());
   traction = 0.0;
   fault.PreInit(state);
   fault.Init(traction, state);

   const int N = fault.NumNodes();
   Vector psi(N);
   for (int i = 0; i < N; i++)
   {
      psi(i) = state(i * cycle::QDPhase::FaultOp::StatePerNode +
                     cycle::QDPhase::FaultOp::PsiIndex);
   }
   return psi;
}

/// One IC set: deck parameters against BP5Params `bp5` on the same domain.
int CheckSet(const std::string &deck, ParMesh &pmesh, MPIContext &mpi,
             const BP5Params &bp5)
{
   const cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   const bool gold = cfg.IsGold();
   if (mpi.IsRoot())
   {
      std::cout << "\n" << deck << "  (ic " << cfg.ic << ", reference BP5Params V_nuc = "
                << bp5.V_nuc << ", nucleation_eps = " << bp5.nucleation_eps
                << ", delta_tau_factor = " << bp5.delta_tau_factor << ")\n";
   }

   cycle::QDPhase qd(pmesh, mpi, cfg);
   const spatial::SpatialFrictionConfig &s = cfg.shared;
   const MaterialField material = MaterialField::MakeConstant(
      s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
   const spatial::PorePressureSpec no_pp;
   const cycle::FaultParams par = cycle::BuildFaultParams(
      qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
   const cycle::Prestress pre = cycle::BuildPrestress(par, cfg.prestress);

   // The reference: seas_bp5_full's geometry on the same domain operator.
   FaultGeometry<ParMesh> ref(qd.Domain(), bp5, &mpi);
   const int N = ref.NumFaultDOFs();
   MFEM_VERIFY(N == qd.Table().Size(), "reference and QD geometry differ in size: "
               << N << " vs " << qd.Table().Size());

   // Nucleation patch and its edge, from the reference's (x2, x3).
   int in_patch = 0, on_edge = 0;
   {
      BP5Params exact = bp5;
      exact.nucleation_eps = 0.0;
      const real_t half_l = bp5.l_vw / 2.0;
      const real_t lo2 = -half_l, hi2 = -half_l + bp5.w_nuc;
      const real_t lo3 = bp5.hs + bp5.ht, hi3 = bp5.hs + bp5.ht + bp5.H;
      for (int i = 0; i < N; i++)
      {
         const real_t x2 = ref.GetCoordsX2()(i), x3 = ref.GetCoordsX3()(i);
         if (!bp5.IsNucleationZone(x2, x3)) { continue; }
         in_patch++;
         const real_t e = 1e-6 * bp5.w_nuc;
         if (std::abs(x2 - lo2) < e || std::abs(x2 - hi2) < e ||
             std::abs(x3 - lo3) < e || std::abs(x3 - hi3) < e)
         {
            on_edge++;
         }
      }
   }
   in_patch = mpi.GlobalSumInt(in_patch);
   on_edge = mpi.GlobalSumInt(on_edge);
   if (mpi.IsRoot())
   {
      std::cout << "  fault DOFs " << ref.NumGlobalFaultDOFs() << ", in the nucleation "
                << "patch " << in_patch << ", on its edge " << on_edge << "\n";
   }

   Checker c{mpi};
   if (in_patch == 0 || on_edge == 0)
   {
      if (mpi.IsRoot()) { std::cout << "  FAIL no DOF in the patch or on its edge\n"; }
      c.failures++;
   }

   // Deck parameters against BP5Params.
   c.Vec("a", par.rs.a, ref.GetAValues());
   c.Vec("Dc", par.rs.Dc, ref.GetDcValues());
   c.Vec("V_init (dip, strike)", par.V_init, ref.GetVInit());
   c.Vec("eta", par.rs.eta, ref.GetEtaValues());
   c.Scalar("b", par.rs.b, bp5.b);
   c.Scalar("f0", par.rs.f_0, bp5.f0);
   c.Scalar("V0", par.rs.V_0, bp5.V0);
   c.Scalar("sigma_n", par.rs.sigma_n_eff, bp5.sigma_n);
   c.Scalar("sigma0", pre.sigma0, bp5.sigma_n);
   c.Scalar("psi0 (BuildPrestress)", pre.psi0, bp5.f0 + bp5.b * std::log(bp5.V0 / bp5.Vp));
   if (gold)
   {
      c.Vec("tau0 (dip, strike)", pre.tau0, ref.GetTauPre());
   }
   else
   {
      // SCEC eq. 22: the reference's tau0 without eta |V|, same direction.
      Vector want(2 * N);
      const Vector &tp = ref.GetTauPre();
      const Vector &vi = ref.GetVInit();
      for (int i = 0; i < N; i++)
      {
         const real_t v0 = vi(2 * i), v1 = vi(2 * i + 1);
         const real_t vabs = std::sqrt(v0 * v0 + v1 * v1);
         const real_t ev = bp5.eta() * vabs;
         want(2 * i) = tp(2 * i) + ev * v0 / vabs;
         want(2 * i + 1) = tp(2 * i + 1) + ev * v1 / vabs;
      }
      c.Vec("tau0 (dip, strike), no eta|V|", pre.tau0, want);
   }

   // Hooks H1a + H1b: what the kernel caches, then the kernel's psi0.
   const Vector psi = KernelPsi0(qd, cfg, par, pre, mpi);
   const FaultGeometry<ParMesh> &g = qd.Geometry();
   c.Vec("H1a geometry a", g.GetAValues(), ref.GetAValues());
   c.Vec("H1a geometry Dc", g.GetDcValues(), ref.GetDcValues());
   c.Vec("H1a geometry eta", g.GetEtaValues(), ref.GetEtaValues());
   c.Vec("H1b geometry V_init", g.GetVInit(), par.V_init);
   c.Vec("H1b geometry tau_pre", g.GetTauPre(), pre.tau0);
   c.Scalar("H1b geometry sigma_n", g.sigma_n_per_dof(), bp5.sigma_n);
   c.Vec("psi0 kernel vs BuildPrestress", psi, pre.psi0);

   if (mpi.IsRoot())
   {
      std::cout << "  " << (c.failures ? "FAILED" : "PASSED") << " (" << cfg.ic << ")\n";
   }
   return c.failures;
}

} // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "config/fd_cycle/bp5_gold.toml";
   const char *deck_fd = "config/fd_cycle/bp5_scec_fd.toml";
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "Gold deck (ic = \"gold\").");
   args.AddOption(&deck_fd, "-dfd", "--deck-fd", "scec_fd deck (ic = \"scec_fd\").");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.IsRoot()) { args.PrintUsage(std::cout); }
      return 1;
   }

   const cycle::CycleConfig gold_cfg = cycle::CycleConfig::Load(deck);
   const cycle::CycleConfig fd_cfg = cycle::CycleConfig::Load(deck_fd);
   MFEM_VERIFY(gold_cfg.IsGold(), deck << " is not a gold deck");
   MFEM_VERIFY(!fd_cfg.IsGold(), deck_fd << " is not a scec_fd deck");
   MFEM_VERIFY(gold_cfg.shared.mesh.path == fd_cfg.shared.mesh.path &&
               gold_cfg.mesh_scale == fd_cfg.mesh_scale,
               "the two decks must use the same mesh");

   std::unique_ptr<Mesh> smesh =
      cycle::LoadScaledMesh(gold_cfg.shared.mesh.path, gold_cfg.mesh_scale);
   ParMesh pmesh(mpi.GetComm(), *smesh);
   smesh.reset();
   if (mpi.IsRoot())
   {
      std::cout << "test_cycle_params_parity: " << gold_cfg.shared.mesh.path << " x "
                << gold_cfg.mesh_scale << ", p = " << gold_cfg.shared.mesh.order
                << ", ranks " << mpi.Size() << ", tolerance " << kTol << "\n";
   }

   // seas_bp5_full's gold parameters are BP5Params' defaults (Tandem:
   // V_nuc = 0.01, nucleation_eps = 1 mm, delta_tau_factor = 0).
   BP5Params bp5_gold;
   BP5Params bp5_fd;
   bp5_fd.V_nuc = 0.03;
   int failures = CheckSet(deck, pmesh, mpi, bp5_gold);
   failures += CheckSet(deck_fd, pmesh, mpi, bp5_fd);

   if (mpi.IsRoot())
   {
      std::cout << "\ntest_cycle_params_parity: " << (failures ? "FAILED" : "PASSED")
                << " (" << failures << " failing checks)\n";
   }
   return failures ? 1 : 0;
}
