// tests/ucerf3d/test_tensor_prestress_3d.cpp - UCERF 3-D plan rev5, P1.5:
// the constant stress tensor on a dipping fault (test U-S1).
//
// The deck tests/ucerf3d/configs/T-dip_tensor.toml on T-dip (a planar fault
// dipping 58 degrees to +y, strike along x), through the cycle's own path:
// CycleConfig::Load, QDPhase (built, not set up), BuildFaultParams,
// QDPhase::ProjectTensor, BuildTensorPrestress and DynPhase (hook H2), as
// tests/2d_inplane/test_tensor_prestress.cpp does on its rotated box.
//
// Analytic values, with the dip delta = 58 degrees and the oriented normal
// (ref_normal (0, -1, 0)):
//    n   = (0, -sin delta, -cos delta)
//    e_s = (1, 0, 0)                          along strike
//    e_d = (0, cos delta, -sin delta)         down dip
//    T   = S n,  sigma_n = -n . T,  shear = T - (n . T) n
// Checks at every owned fault DOF, 1e-12 relative:
//   S1  sigma0 = -n . S n and |tau0| = |shear|
//   S2  the kernel's two shear slots: |tau0_1| = |e_d . T| (dip) and
//       |tau0_2| = |e_s . T| (strike), each with one sign on the whole fault
//   S3  V_init = -|V_init| tau0 / |tau0| with |V_init| the deck's, and psi0 =
//       a ln[(2 V0 / V) sinh(|tau0| / (sigma0 a))] (1e-12 absolute)
//   S4  the dynamic phase's background traction (hook H2: FrameMap::
//       TractionToGlobal(tau0, sigma0)) is S n, and the frame's normal is n
//
// Usage (from miniapps/seas; PETSc needed):
//   mpirun -np N ./seas_test_tensor_prestress_3d [--deck tests/ucerf3d/configs/T-dip_tensor.toml]

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
#include <sstream>

using namespace mfem;
using namespace mfem::seas;

namespace
{
int num_tests = 0, num_failed = 0;
void Check(MPIContext &mpi, bool ok, const std::string &what)
{
   num_tests++;
   if (!ok) { num_failed++; }
   if (mpi.IsRoot()) { std::cout << (ok ? "  PASSED: " : "  FAILED: ") << what << std::endl; }
}
real_t GMax(MPIContext &mpi, real_t v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, mpi.GetComm());
   return v;
}
real_t Dot3(const real_t *a, const real_t *b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
}  // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "tests/ucerf3d/configs/T-dip_tensor.toml";
   real_t dip_deg = 58.0;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "The U-S1 deck.");
   args.AddOption(&dip_deg, "-dip", "--dip", "Dip of the deck's fault (degrees, toward +y).");
   args.Parse();
   if (!args.Good()) { if (mpi.IsRoot()) { args.PrintUsage(std::cout); } return 1; }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_tensor_prestress_3d needs PETSc.\n"; }
   return 2;
#else
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   {
      if (mpi.IsRoot())
      {
         std::cout << "test_tensor_prestress_3d (UCERF 3-D plan rev5, P1.5, U-S1) on "
                   << mpi.Size() << " rank(s): " << deck << ", IC " << cfg.IcLabel() << std::endl;
      }
      MFEM_VERIFY(cfg.TensorPrestress(), deck << ": [stress] kind = \"constant_tensor\" is needed");
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      ParMesh pmesh(mpi.GetComm(), *smesh);
      smesh.reset();
      cycle::QDPhase qd(pmesh, mpi, cfg);
      const spatial::SpatialFrictionConfig &s = cfg.shared;
      const MaterialField material = MaterialField::MakeConstant(
         s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
      const spatial::PorePressureSpec no_pp;
      cycle::FaultParams par = cycle::BuildFaultParams(
         qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
      Vector vinit_abs(par.rs.V_init);   // the deck's |V_init| per DOF
      Vector K;
      qd.ProjectTensor(cfg.stress_tensor, K);
      const cycle::Prestress pre = cycle::BuildTensorPrestress(par, K, cfg.prestress, qd.Table());
      const int N = qd.Table().Size();
      long long ng = N;
      MPI_Allreduce(MPI_IN_PLACE, &ng, 1, MPI_LONG_LONG, MPI_SUM, mpi.GetComm());

      // analytic
      const real_t del = dip_deg * M_PI / 180.0, sd = std::sin(del), cd = std::cos(del);
      const real_t n[3] = {0.0, -sd, -cd}, es[3] = {1.0, 0.0, 0.0}, ed[3] = {0.0, cd, -sd};
      const std::array<real_t, 9> &S = cfg.stress_tensor;
      real_t T[3], sh[3];
      for (int i = 0; i < 3; i++) { T[i] = S[3 * i] * n[0] + S[3 * i + 1] * n[1] + S[3 * i + 2] * n[2]; }
      const real_t Tn = Dot3(T, n), sn_a = -Tn;
      for (int i = 0; i < 3; i++) { sh[i] = T[i] - Tn * n[i]; }
      const real_t tau_a = std::sqrt(Dot3(sh, sh)), td_a = Dot3(ed, T), ts_a = Dot3(es, T);
      if (mpi.IsRoot())
      {
         std::printf("  analytic: n = (0, %.6f, %.6f); T = S n = (%.6f, %.6f, %.6f) MPa; sigma_n = "
                     "%.6f MPa, |shear| = %.6f MPa (down dip %.6f, along strike %.6f)\n", n[1],
                     n[2], T[0] / 1e6, T[1] / 1e6, T[2] / 1e6, sn_a / 1e6, tau_a / 1e6, td_a / 1e6,
                     ts_a / 1e6);
      }

      // S1 - S3
      const real_t V0 = s.rate_state->V_0_default;
      real_t e_sn = 0, e_tau = 0, e_d = 0, e_s = 0, e_v = 0, e_dir = 0, e_psi = 0;
      real_t sgn1[2] = {0, 0}, sgn2[2] = {0, 0};   // (count of +, count of -)
      for (int i = 0; i < N; i++)
      {
         const real_t t1 = pre.tau0(2 * i), t2 = pre.tau0(2 * i + 1), sn = pre.sigma0(i);
         const real_t tau = std::hypot(t1, t2), a = par.rs.a(i);
         e_sn = std::max(e_sn, std::abs(sn - sn_a) / sn_a);
         e_tau = std::max(e_tau, std::abs(tau - tau_a) / tau_a);
         e_d = std::max(e_d, std::abs(std::abs(t1) - std::abs(td_a)) / tau_a);
         e_s = std::max(e_s, std::abs(std::abs(t2) - std::abs(ts_a)) / tau_a);
         sgn1[t1 > 0 ? 0 : 1] += 1;
         sgn2[t2 > 0 ? 0 : 1] += 1;
         const real_t v1 = par.V_init(2 * i), v2 = par.V_init(2 * i + 1), va = std::hypot(v1, v2);
         e_v = std::max(e_v, std::abs(va - vinit_abs(i)) / vinit_abs(i));
         e_dir = std::max(e_dir, std::hypot(v1 + va * t1 / tau, v2 + va * t2 / tau) / va);
         const real_t psi = a * std::log(2.0 * V0 / va * std::sinh(tau / (sn * a)));
         e_psi = std::max(e_psi, std::abs(pre.psi0(i) - psi));
      }
      MPI_Allreduce(MPI_IN_PLACE, sgn1, 2, MPITypeMap<real_t>::mpi_type, MPI_SUM, mpi.GetComm());
      MPI_Allreduce(MPI_IN_PLACE, sgn2, 2, MPITypeMap<real_t>::mpi_type, MPI_SUM, mpi.GetComm());
      e_sn = GMax(mpi, e_sn); e_tau = GMax(mpi, e_tau); e_d = GMax(mpi, e_d); e_s = GMax(mpi, e_s);
      e_v = GMax(mpi, e_v); e_dir = GMax(mpi, e_dir); e_psi = GMax(mpi, e_psi);
      std::ostringstream os;
      os << "U-S1 S1: at " << ng << " fault DOFs sigma0 = -n.S.n (" << e_sn << ") and |tau0| = "
         << "|S n - (n.S n) n| (" << e_tau << "), 1e-12 relative";
      Check(mpi, ng > 0 && e_sn <= 1e-12 && e_tau <= 1e-12, os.str());
      os.str("");
      os << "U-S1 S2: the dip slot is |e_d . S n| (" << e_d << ") and the strike slot |e_s . S n| ("
         << e_s << "), 1e-12; one sign each on the whole fault (dip +" << sgn1[0] << " / -"
         << sgn1[1] << ", strike +" << sgn2[0] << " / -" << sgn2[1] << ")";
      Check(mpi, e_d <= 1e-12 && e_s <= 1e-12 && (sgn1[0] == 0 || sgn1[1] == 0)
            && (sgn2[0] == 0 || sgn2[1] == 0), os.str());
      os.str("");
      os << "U-S1 S3: V_init = -|V_init| tau0 / |tau0| (size " << e_v << ", direction " << e_dir
         << ", 1e-12) and psi0 is the balance at V_init (" << e_psi << " absolute)";
      Check(mpi, e_v <= 1e-12 && e_dir <= 1e-12 && e_psi <= 1e-12, os.str());

      // S4: the dynamic phase's background traction
      cycle::DynPhase dyn(pmesh, cfg, qd.Table(), par, pre, 1.0, nullptr, qd.Domain().FrameField());
      Vector Tg;
      dyn.Frames().TractionToGlobal(pre.tau0, pre.sigma0, Tg);
      real_t e_T = 0, e_n = 0;
      const real_t Tabs = std::sqrt(Dot3(T, T));
      const Vector &fr = dyn.Op().FaultFrames();
      for (int i = 0; i < N; i++)
      {
         real_t d2 = 0, n2 = 0;
         for (int c = 0; c < 3; c++)
         {
            d2 += (Tg(3 * i + c) - T[c]) * (Tg(3 * i + c) - T[c]);
            n2 += (fr(9 * i + c) - n[c]) * (fr(9 * i + c) - n[c]);
         }
         e_T = std::max(e_T, std::sqrt(d2) / Tabs);
         e_n = std::max(e_n, std::sqrt(n2));
      }
      e_T = GMax(mpi, e_T);
      e_n = GMax(mpi, e_n);
      os.str("");
      os << "U-S1 S4: the dynamic phase's background traction is S n (" << e_T
         << " relative) and its fault normal is (0, -sin 58, -cos 58) (" << e_n << "), 1e-12";
      Check(mpi, dyn.Op().NumOwnedFaultDOFs() == N && e_T <= 1e-12 && e_n <= 1e-12, os.str());
   }
   MFEMFinalizePetsc();
   if (mpi.IsRoot()) { std::cout << (num_tests - num_failed) << "/" << num_tests << " tests passed" << std::endl; }
   return num_failed == 0 ? 0 : 1;
#endif
}
