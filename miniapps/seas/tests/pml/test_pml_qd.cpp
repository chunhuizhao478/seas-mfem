// tests/pml/test_pml_qd.cpp - PML plan rev2, P1 at the QD operator level
// (T1 and the D3 loading check, as static solves; no time stepping).
//
// The deck (default simple_cases/configs/Q0_straight.toml: straight fault
// x in [-15, 15] km on y = 0, +-240 km box, walls y = +-240 km attribute 5
// loaded, sides x = +-240 km attribute 1 natural) is built three times:
//   plain   no band
//   band    [cycle.pml] on all four walls, wall speed factor 1
//   mapped  no band, on the mesh whose band vertices are moved to x~(x)
//           (Proposition 1: the same static problem as "band")
//   auto    the band with wall_speed_factor = "auto" (decision D3)
// Each one solves two static QD problems with the deck's operator:
//   load    zero slip, the loading at t = 1 s (the fault stressing rate)
//   slip    the plate at rest (t = 0), 1 m of uniform strike slip
// and the fault shear traction is compared at the owned fault DOFs:
//   T1      band against mapped: |dtau| / |tau| small (the mesh carries the
//           curved map x~ with straight elements, so it is a discretization
//           difference, printed with the threshold of the review)
//   D3      band / plain and auto / plain, the fault-averaged stressing rate,
//           printed against the estimate W / (W + 2 Delta L) of (4.5)
//   slip    band and plain: the slip response must hardly change (the band
//           is 225 km from the fault): relative difference printed
//
// Usage (from miniapps/seas; PETSc needed):
//   mpirun -np N ./seas_test_pml_qd [--deck ...] [--thickness-km 100] [--f-alpha 0.005]

#include "mfem.hpp"
#include "miniapps/seas/cycle/qd_phase.hpp"
#include "miniapps/seas/cycle/scaled_mesh.hpp"
#include "miniapps/seas/common/mpi_context.hpp"

#include <cmath>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

namespace
{

struct Result
{
   Vector tau_load, tau_slip;   // fault shear traction (strike slot), owned DOFs
};

// x~ of Proposition 1 at vertex x
void XTilde(const PmlBand &band, const real_t *x, real_t *xt)
{
   const int D = band.Dim();
   for (int i = 0; i < D; i++) { xt[i] = x[i]; }
   PmlPoint p = band.At(x);
   const real_t L = band.Params().thickness;
   const int N = band.Params().power;
   for (int k = 0; k < p.nact; k++)
   {
      const real_t Da = band.D0() * L / (N + 1) * std::pow(p.xi[k] / L, N + 1);
      for (int i = 0; i < D; i++) { xt[i] += Da / band.Alpha() * band.Walls()[p.wall[k]].n[i]; }
   }
}

PmlBandParams BandParams(const cycle::PmlConfig &P)
{
   PmlBandParams prm;
   prm.thickness = P.thickness;
   prm.power = P.power;
   prm.R0 = P.R0;
   prm.f_alpha = P.f_alpha;
   prm.c_ref = P.c_ref;
   return prm;
}

Result Solve(MPIContext &mpi, const cycle::CycleConfig &cfg, const PmlBand *map_band)
{
   std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
   if (map_band)
   {
      for (int v = 0; v < smesh->GetNV(); v++)
      {
         real_t xt[3];
         real_t *x = smesh->GetVertex(v);
         XTilde(*map_band, x, xt);
         for (int i = 0; i < smesh->Dimension(); i++) { x[i] = xt[i]; }
      }
   }
   ParMesh pmesh(mpi.GetComm(), *smesh);
   smesh.reset();
   cycle::QDPhase qd(pmesh, mpi, cfg);
   auto &dom = qd.Domain();
   const int nf = dom.GetNumFaultDOFs();
   const int nsc = dom.NumSlipComponents();
   ParGridFunction u(&dom.GetFESpace());
   Vector slip(nsc * nf), trac, own;
   Result r;
   for (int c = 0; c < 2; c++)
   {
      slip = 0.0;
      real_t t = 1.0;
      if (c == 1)
      {
         t = 0.0;
         for (int i = 0; i < nsc * nf; i++) { slip(i) = 1.0; }   // both slots; 2-D uses one
      }
      u = 0.0;
      dom.Solve(t, slip, u);
      dom.ComputeTraction(u, slip, trac);
      dom.RestrictToOwnedFault(trac, own, nsc);
      Vector tau(own.Size() / nsc);
      // the shear slot carrying the response: the one of larger norm
      real_t n0 = 0.0, n1 = 0.0;
      for (int i = 0; i < tau.Size(); i++)
      {
         n0 += own(i * nsc) * own(i * nsc);
         n1 += own(i * nsc + 1) * own(i * nsc + 1);
      }
      MPI_Allreduce(MPI_IN_PLACE, &n0, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, mpi.GetComm());
      MPI_Allreduce(MPI_IN_PLACE, &n1, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, mpi.GetComm());
      const int slot = (n1 > n0) ? 1 : 0;
      for (int i = 0; i < tau.Size(); i++) { tau(i) = own(i * nsc + slot); }
      (c == 0 ? r.tau_load : r.tau_slip) = tau;
   }
   return r;
}

real_t GMax(MPIContext &mpi, real_t v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, mpi.GetComm());
   return v;
}
real_t GSum(MPIContext &mpi, real_t v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, mpi.GetComm());
   return v;
}

// max |a - b| / max |b|
real_t RelDiff(MPIContext &mpi, const Vector &a, const Vector &b)
{
   real_t d = 0.0, n = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      d = std::max(d, std::abs(a(i) - b(i)));
      n = std::max(n, std::abs(b(i)));
   }
   return GMax(mpi, d) / GMax(mpi, n);
}

// mean(a) / mean(b) over the owned DOFs
real_t MeanRatio(MPIContext &mpi, const Vector &a, const Vector &b)
{
   return GSum(mpi, a.Sum()) / GSum(mpi, b.Sum());
}

}  // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "simple_cases/configs/Q0_straight.toml";
   real_t L_km = 100.0, f_alpha = 0.005;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "2-D QD deck with loaded walls (attribute 5).");
   args.AddOption(&L_km, "-L", "--thickness-km", "Band thickness [km].");
   args.AddOption(&f_alpha, "-fa", "--f-alpha", "Shift frequency [Hz].");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.IsRoot()) { args.PrintUsage(std::cout); }
      return 1;
   }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_pml_qd needs PETSc.\n"; }
   return 2;
#else
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   MFEM_VERIFY(!cfg.switching, deck << ": a QD-only deck is needed");
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   int fails = 0;
   {
      cycle::CycleConfig cb = cfg;
      cb.pml.on = true;
      cb.pml.wall_attrs = {5, 1};
      cb.pml.thickness = 1e3 * L_km;
      cb.pml.f_alpha = f_alpha;
      const real_t lam = cfg.shared.material_fallback.lambda, mu = cfg.shared.material_fallback.mu;
      cb.pml.c_ref = std::sqrt((lam + 2.0 * mu) / cfg.shared.material_fallback.rho);
      cb.pml.wall_speed_factor = 1.0;
      cycle::CycleConfig ca = cb;
      ca.pml.wall_speed_factor = 0.0;   // "auto"

      // the band on the serial mesh, for the mapped mesh and the estimate
      std::unique_ptr<Mesh> sm = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      PmlBand band(*sm, cb.pml.wall_attrs, BandParams(cb.pml));
      const real_t m[3] = {cfg.loading.plane_normal[0], cfg.loading.plane_normal[1], 0.0};
      const real_t W = band.OppositeWallDistance(m);
      const real_t est = W / (W + 2.0 * band.DeltaL());
      sm.reset();
      if (mpi.IsRoot())
      {
         std::cout << "test_pml_qd: " << deck << ", ranks " << mpi.Size() << ", L = " << L_km
                   << " km, f_alpha = " << f_alpha << " Hz, s0 wall = " << band.S0Wall()
                   << ", Delta L = " << band.DeltaL() / 1e3 << " km, W = " << W / 1e3
                   << " km\n";
      }
      const Result plain = Solve(mpi, cfg, nullptr);
      const Result bnd = Solve(mpi, cb, nullptr);
      const Result mapped = Solve(mpi, cfg, &band);
      const Result aut = Solve(mpi, ca, nullptr);

      auto check = [&](bool ok, const char *what, real_t v)
      {
         if (mpi.IsRoot()) { std::printf("  [%s] %s  (%.4e)\n", ok ? "PASS" : "FAIL", what, v); }
         if (!ok) { fails++; }
      };
      const real_t t1_load = RelDiff(mpi, bnd.tau_load, mapped.tau_load);
      const real_t t1_slip = RelDiff(mpi, bnd.tau_slip, mapped.tau_slip);
      check(t1_load < 0.05, "T1 load: band vs mapped mesh, max |dtau| / max |tau| < 5%", t1_load);
      check(t1_slip < 0.01, "T1 slip: band vs mapped mesh, max |dtau| / max |tau| < 1%", t1_slip);
      const real_t r_band = MeanRatio(mpi, bnd.tau_load, plain.tau_load);
      const real_t r_map = MeanRatio(mpi, mapped.tau_load, plain.tau_load);
      const real_t r_auto = MeanRatio(mpi, aut.tau_load, plain.tau_load);
      if (mpi.IsRoot())
      {
         std::printf("  stressing rate (fault mean) / plain: band %.4f, mapped %.4f, "
                     "estimate W/(W+2 dL) %.4f; auto %.4f (factor %.4f)\n",
                     r_band, r_map, est, r_auto, 1.0 / est);
      }
      check(r_band < 1.0, "the band lowers the stressing rate", r_band);
      check(std::abs(r_auto / r_band - 1.0 / est) < 1e-6 * (1.0 / est),
            "auto = band x (W + 2 dL) / W (the loading is linear in the wall speed)",
            r_auto / r_band);
      const real_t slip_d = RelDiff(mpi, bnd.tau_slip, plain.tau_slip);
      check(slip_d < 0.01, "slip response with the band vs plain < 1% (band far away)", slip_d);
      if (mpi.IsRoot())
      {
         std::printf("  D3 check: auto / plain = %.4f (1 = the plain box's fault loading)\n",
                     r_auto);
      }
   }
   if (mpi.IsRoot()) { std::cout << (fails ? "FAILED: " : "ALL PASS: ") << fails << "\n"; }
   MFEMFinalizePetsc();
   return fails ? 1 : 0;
#endif
}
