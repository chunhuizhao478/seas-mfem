// tests/ucerf3d/test_device_cycle_3d.cpp - UCERF 3-D plan rev5, P1.7: the
// cycle's dynamic phase on the mfem device against the host (test U-G1).
//
// The deck tests/ucerf3d/configs/T-plane_device_step.toml through the
// cycle's own path: CycleConfig::Load, QDPhase (built, not set up),
// BuildFaultParams, BuildPrestress, DynPhase, bridge::ToNitsche, Enter (the
// switch-in), --steps steps and the hand-back state (slip, slip rate, state
// variable at the owned fault DOFs, and the displacement and velocity
// fields).  [cycle].device_step is set from the device: off on --device cpu
// (the host step), on with any other device (the matrix-free product and
// the device fault).
//
// One mfem device per process, so the comparison takes two runs:
//   --device cpu   --write FILE     the host step's hand-back state
//   --device debug --compare FILE   the device step's against it
// Checks of the compare run, relative to the largest value of each field.
// Limit 1e-12 (Chunhui Zhao 2026-11-26; the plan said 1e-13): the two steps
// order their sums differently, and the round-off grows with the steps
// (slip rate, 2 ranks: 5e-15 at 30 steps, 3e-14 at 100, 1.1e-13 at 300,
// 3e-13 at 1000).
//   G1  the device step ran (matrix-free product, device step, device fault)
//   G2  the same fault points, and the event ran (peak slip rate)
//   G3  slip, slip rate and state at every fault DOF
//   G4  the largest |u| and |v| of the volume fields
//
// Usage (from miniapps/seas; PETSc needed):
//   mpirun -np N ./seas_test_device_cycle_3d --device cpu --write /tmp/ug1.bin
//   mpirun -np N ./seas_test_device_cycle_3d --device debug --compare /tmp/ug1.bin

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"
#include "../../common/mpi_context.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

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
const int NR = 10;   // x (3), slip (3), slip rate (3), psi
}  // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "tests/ucerf3d/configs/T-plane_device_step.toml";
   const char *device_str = "cpu", *write = "", *compare = "";
   int steps = 300;
   real_t tol = 1e-12;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "The U-G1 deck.");
   args.AddOption(&device_str, "-dev", "--device", "mfem device: cpu (host step), debug, cuda.");
   args.AddOption(&steps, "-n", "--steps", "Steps of the dynamic phase.");
   args.AddOption(&write, "-w", "--write", "Write the hand-back state to this file.");
   args.AddOption(&compare, "-c", "--compare", "Compare the hand-back state with this file.");
   args.AddOption(&tol, "-tol", "--tol", "Tolerance of the comparison (relative to the largest).");
   args.Parse();
   if (!args.Good()) { if (mpi.IsRoot()) { args.PrintUsage(std::cout); } return 1; }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_device_cycle_3d needs PETSc.\n"; }
   return 2;
#else
   Device device(device_str);
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   cfg.device_step = std::string(device_str) != "cpu";
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   {
      if (mpi.IsRoot())
      {
         std::cout << "test_device_cycle_3d (UCERF 3-D plan rev5, P1.7, U-G1) on " << mpi.Size()
                   << " rank(s): " << deck << ", [cycle].device_step = "
                   << (cfg.device_step ? "true" : "false") << "\n";
         device.Print();
      }
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      ParMesh pmesh(mpi.GetComm(), *smesh);
      smesh.reset();
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
      cycle::DynPhase dyn(pmesh, cfg, qd.Table(), par, pre, 1.0, nullptr, qd.Domain().FrameField());
      NitscheWaveOperator &W = dyn.Op();
      const int N = qd.Table().Size();
      // the switch-in
      Vector sk(2 * N), d, V, psi, us(W.Height());
      sk = 0.0;
      us = 0.0;
      cycle::bridge::ToNitsche(dyn.Frames(), sk, par.V_init, pre.psi0, d, V, psi);
      dyn.Enter(us, d, V, psi, pre.sigma0, 0.0);
      real_t peak = 0.0;
      StopWatch sw;
      sw.Start();
      for (int k = 0; k < steps; k++)
      {
         dyn.Step();
         if ((k + 1) % 50 == 0 || k + 1 == steps)
         {
            // an output read of the fault state (copied from the device when read)
            const Vector &Vf = W.SlipRate();
            Vf.HostRead();
            real_t m = 0.0;
            for (int j = 0; j < N; j++)
            {
               m = std::max(m, std::hypot(Vf(3 * j), Vf(3 * j + 1), Vf(3 * j + 2)));
            }
            m = GMax(mpi, m);
            peak = std::max(peak, m);
            if (mpi.IsRoot())
            {
               std::printf("    step %d, tau %.4f s, max |V| %.6f m/s, %.1f s wall\n", k + 1,
                           dyn.Tau(), m, sw.RealTime());
               std::fflush(stdout);
            }
         }
      }
      // the hand-back state
      const Vector &X = W.FaultPositions(), &D = W.Slip(), &Vf = W.SlipRate(), &Ps = W.State();
      X.HostRead(); D.HostRead(); Vf.HostRead(); Ps.HostRead();
      std::vector<double> loc(NR * static_cast<size_t>(N)), all;
      for (int j = 0; j < N; j++)
      {
         for (int c = 0; c < 3; c++)
         {
            loc[NR * j + c] = X(3 * j + c);
            loc[NR * j + 3 + c] = D(3 * j + c);
            loc[NR * j + 6 + c] = Vf(3 * j + c);
         }
         loc[NR * j + 9] = Ps(j);
      }
      GatherRows(mpi.GetComm(), loc, all);
      dyn.U().HostRead();
      dyn.Ud().HostRead();
      const real_t umax = GMax(mpi, dyn.U().Normlinf()), vmax = GMax(mpi, dyn.Ud().Normlinf());
      const bool mf = W.MatrixFreeActive(), ds = W.DeviceStepActive(), df = W.DeviceFaultActive();
      if (mpi.IsRoot())
      {
         std::printf("  %d steps of %.6f s to tau = %.4f s; peak |V| %.6f m/s; max |u| %.6e m, "
                     "max |v| %.6e m/s; matrix-free %d, device step %d, device fault %d\n",
                     steps, dyn.Dt(), dyn.Tau(), peak, umax, vmax, mf, ds, df);
      }
      // rows sorted by position (the same order on any device and rank count)
      const size_t M = all.size() / NR;
      std::vector<std::array<double, NR>> rows(M);
      for (size_t i = 0; i < M; i++) { std::copy_n(all.begin() + NR * i, NR, rows[i].begin()); }
      std::sort(rows.begin(), rows.end());
      if (mpi.IsRoot() && write[0])
      {
         std::ofstream f(write, std::ios::binary);
         const double head[5] = {static_cast<double>(M), static_cast<double>(steps), dyn.Tau(),
                                 umax, vmax};
         f.write(reinterpret_cast<const char *>(head), sizeof(head));
         for (const auto &r : rows) { f.write(reinterpret_cast<const char *>(r.data()), sizeof(double) * NR); }
         MFEM_VERIFY(f.good(), "cannot write " << write);
         std::cout << "  wrote " << M << " fault rows to " << write << std::endl;
      }
      if (compare[0])
      {
         std::ostringstream os;
         os << "U-G1 G1: the device step ran: matrix-free product " << mf << ", device step "
            << ds << ", device fault " << df << " (device \"" << device_str << "\")";
         Check(mpi, cfg.device_step && mf && ds && df, os.str());
         if (mpi.IsRoot())
         {
            std::ifstream f(compare, std::ios::binary);
            double head[5] = {0, 0, 0, 0, 0};
            f.read(reinterpret_cast<char *>(head), sizeof(head));
            MFEM_VERIFY(f.good(), "cannot read " << compare);
            const size_t Mr = static_cast<size_t>(head[0]);
            std::vector<std::array<double, NR>> ref(Mr);
            for (auto &r : ref) { f.read(reinterpret_cast<char *>(r.data()), sizeof(double) * NR); }
            MFEM_VERIFY(f.good(), "short file " << compare);
            double ex = 0, ed = 0, ev = 0, ep = 0, md = 0, mv = 0, mp = 0;
            const bool same = (Mr == M) && static_cast<int>(head[1]) == steps && head[2] == dyn.Tau();
            for (size_t i = 0; same && i < M; i++)
            {
               for (int c = 0; c < 3; c++)
               {
                  ex = std::max(ex, std::abs(rows[i][c] - ref[i][c]));
                  ed = std::max(ed, std::abs(rows[i][3 + c] - ref[i][3 + c]));
                  ev = std::max(ev, std::abs(rows[i][6 + c] - ref[i][6 + c]));
                  md = std::max(md, std::abs(ref[i][3 + c]));
                  mv = std::max(mv, std::abs(ref[i][6 + c]));
               }
               ep = std::max(ep, std::abs(rows[i][9] - ref[i][9]));
               mp = std::max(mp, std::abs(ref[i][9]));
            }
            os.str("");
            os << "U-G1 G2: the same " << M << " fault points (" << ex << " m), " << steps
               << " steps to tau = " << dyn.Tau() << " s, and the event ran (peak |V| " << peak
               << " m/s > 1)";
            Check(mpi, same && M > 0 && ex == 0.0 && peak > 1.0, os.str());
            os.str("");
            os << "U-G1 G3: device against host at the hand-back: slip " << ed << " m of " << md
               << " (" << ed / md << "), slip rate " << ev << " m/s of " << mv << " ("
               << ev / mv << "), state " << ep << " of " << mp << " (" << ep / mp
               << "); relative limit " << tol;
            Check(mpi, same && ed <= tol * md && ev <= tol * mv && ep <= tol * mp, os.str());
            os.str("");
            const double eu = std::abs(umax - head[3]) / head[3], evv = std::abs(vmax - head[4]) / head[4];
            os << "U-G1 G4: the largest |u| (" << eu << ") and |v| (" << evv
               << ") of the volume fields, relative, limit " << tol;
            Check(mpi, eu <= tol && evv <= tol, os.str());
         }
         else { num_tests += 3; }
      }
   }
   MFEMFinalizePetsc();
   if (mpi.IsRoot() && compare[0])
   {
      std::cout << (num_tests - num_failed) << "/" << num_tests << " tests passed" << std::endl;
   }
   return num_failed == 0 ? 0 : 1;
#endif
}
