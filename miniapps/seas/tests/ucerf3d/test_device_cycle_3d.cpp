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
//   G5  (U-G2) the floor's counts: DOF-steps at the floor and distinct DOFs
//
// U-G2 (P1.8): the same run with the quench, the relaxation r and the floor
// of the law's normal stress, which the device step has since P1.8:
//   --quench-from K --quench-to L   SetQuench after step K, ClearQuench after
//                                   step L (-1: never), rate --quench-d
//   --relax 1                       r = --relax-lo Pa at x < 0 and --relax-hi
//                                   Pa at x >= 0 (DynPhase::SetNormalStressRelaxation)
//   --floor F                       [cycle].sigma_n_floor_pa = F
// With sigma_n = 25 MPa, r = 1 MPa / 0.2 MPa and the floor at 24.5 MPa the
// law's normal stress is below the floor at x < 0 and above it at x >= 0.
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
   int steps = 300, quench_from = -1, quench_to = -1;
   bool relax = false;
   real_t tol = 1e-12, min_peak = 1.0, quench_d = 0.5, relax_lo = 1.0e6, relax_hi = 0.2e6, floor = 0.0;
   const char *name = "U-G1";
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "The U-G1 deck.");
   args.AddOption(&device_str, "-dev", "--device", "mfem device: cpu (host step), debug, cuda.");
   args.AddOption(&steps, "-n", "--steps", "Steps of the dynamic phase.");
   args.AddOption(&write, "-w", "--write", "Write the hand-back state to this file.");
   args.AddOption(&compare, "-c", "--compare", "Compare the hand-back state with this file.");
   args.AddOption(&tol, "-tol", "--tol", "Tolerance of the comparison (relative to the largest).");
   args.AddOption(&min_peak, "-mp", "--min-peak", "G2: the peak slip rate must exceed this [m/s].");
   args.AddOption(&name, "-name", "--name", "Test name in the check lines (U-G1, U-G2 ...).");
   args.AddOption(&quench_from, "-qf", "--quench-from", "SetQuench after this step (-1: no quench).");
   args.AddOption(&quench_to, "-qt", "--quench-to", "ClearQuench after this step (-1: never).");
   args.AddOption(&quench_d, "-qd", "--quench-d", "Quench rate d [1/s].");
   args.AddOption(&relax, "-relax", "--relax", "-no-relax", "--no-relax",
                  "Relaxation r: --relax-lo at x < 0, --relax-hi at x >= 0.");
   args.AddOption(&relax_lo, "-rlo", "--relax-lo", "r at x < 0 [Pa].");
   args.AddOption(&relax_hi, "-rhi", "--relax-hi", "r at x >= 0 [Pa].");
   args.AddOption(&floor, "-floor", "--floor", "[cycle].sigma_n_floor_pa [Pa] (0: off).");
   args.Parse();
   if (!args.Good()) { if (mpi.IsRoot()) { args.PrintUsage(std::cout); } return 1; }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_device_cycle_3d needs PETSc.\n"; }
   return 2;
#else
   Device device(device_str);
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   cfg.device_step = std::string(device_str) != "cpu";
   if (quench_from >= 0)
   {
      cfg.quench.on = true;
      cfg.quench.d = quench_d;
   }
   if (floor > 0.0) { cfg.sigma_n_floor = floor; }
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   {
      if (mpi.IsRoot())
      {
         std::cout << "test_device_cycle_3d (UCERF 3-D plan rev5, P1.7 and P1.8, " << name
                   << ") on " << mpi.Size() << " rank(s): " << deck << ", [cycle].device_step = "
                   << (cfg.device_step ? "true" : "false") << "; quench after step "
                   << quench_from << " to " << quench_to << " (d = " << quench_d
                   << " 1/s), relaxation " << (relax ? 1 : 0) << " (" << relax_lo << ", "
                   << relax_hi << " Pa), floor " << floor << " Pa\n";
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
      // U-G2: the relaxation r per owned fault DOF, by position
      Vector rr;
      if (relax)
      {
         const Vector &Xf = W.FaultPositions();
         Xf.HostRead();
         rr.SetSize(N);
         for (int j = 0; j < N; j++) { rr(j) = (Xf(3 * j) < 0.0) ? relax_lo : relax_hi; }
         dyn.SetNormalStressRelaxation(rr);
      }
      // the switch-in
      Vector sk(2 * N), d, V, psi, us(W.Height());
      sk = 0.0;
      us = 0.0;
      cycle::bridge::ToNitsche(dyn.Frames(), sk, par.V_init, pre.psi0, d, V, psi);
      dyn.Enter(us, d, V, psi, pre.sigma0, 0.0);
      real_t peak = 0.0;
      StopWatch sw;
      sw.Start();
      int quench_steps = 0;
      for (int k = 0; k < steps; k++)
      {
         if (quench_from >= 0 && k == quench_from) { dyn.SetQuench(); }
         if (quench_to >= 0 && k == quench_to && dyn.Quenched()) { dyn.ClearQuench(); }
         if (dyn.Quenched()) { quench_steps++; }
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
      long fhits = 0, fdofs = 0;
      dyn.FloorCount(fhits, fdofs);
      if (mpi.IsRoot())
      {
         std::printf("  %d steps of %.6f s to tau = %.4f s; peak |V| %.6f m/s; max |u| %.6e m, "
                     "max |v| %.6e m/s; matrix-free %d, device step %d, device fault %d; "
                     "%d quenched steps; floor: %ld DOF-steps, %ld DOFs\n",
                     steps, dyn.Dt(), dyn.Tau(), peak, umax, vmax, mf, ds, df, quench_steps,
                     fhits, fdofs);
      }
      // rows sorted by position (the same order on any device and rank count)
      const size_t M = all.size() / NR;
      std::vector<std::array<double, NR>> rows(M);
      for (size_t i = 0; i < M; i++) { std::copy_n(all.begin() + NR * i, NR, rows[i].begin()); }
      std::sort(rows.begin(), rows.end());
      if (mpi.IsRoot() && write[0])
      {
         std::ofstream f(write, std::ios::binary);
         const double head[8] = {static_cast<double>(M), static_cast<double>(steps), dyn.Tau(),
                                 umax, vmax, static_cast<double>(fhits),
                                 static_cast<double>(fdofs), static_cast<double>(quench_steps)};
         f.write(reinterpret_cast<const char *>(head), sizeof(head));
         for (const auto &r : rows) { f.write(reinterpret_cast<const char *>(r.data()), sizeof(double) * NR); }
         MFEM_VERIFY(f.good(), "cannot write " << write);
         std::cout << "  wrote " << M << " fault rows to " << write << std::endl;
      }
      if (compare[0])
      {
         std::ostringstream os;
         os << name << " G1: the device step ran: matrix-free product " << mf << ", device step "
            << ds << ", device fault " << df << " (device \"" << device_str << "\")";
         Check(mpi, cfg.device_step && mf && ds && df, os.str());
         if (mpi.IsRoot())
         {
            std::ifstream f(compare, std::ios::binary);
            double head[8] = {0, 0, 0, 0, 0, 0, 0, 0};
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
            os << name << " G2: the same " << M << " fault points (" << ex << " m), " << steps
               << " steps to tau = " << dyn.Tau() << " s, and the event ran (peak |V| " << peak
               << " m/s > " << min_peak << ")";
            Check(mpi, same && M > 0 && ex == 0.0 && peak > min_peak, os.str());
            os.str("");
            os << name << " G3: device against host at the hand-back: slip " << ed << " m of " << md
               << " (" << ed / md << "), slip rate " << ev << " m/s of " << mv << " ("
               << ev / mv << "), state " << ep << " of " << mp << " (" << ep / mp
               << "); relative limit " << tol;
            Check(mpi, same && ed <= tol * md && ev <= tol * mv && ep <= tol * mp, os.str());
            os.str("");
            const double eu = std::abs(umax - head[3]) / head[3], evv = std::abs(vmax - head[4]) / head[4];
            os << name << " G4: the largest |u| (" << eu << ") and |v| (" << evv
               << ") of the volume fields, relative, limit " << tol;
            Check(mpi, eu <= tol && evv <= tol, os.str());
            os.str("");
            os << name << " G5: the floor's counts " << fhits << " DOF-steps and " << fdofs
               << " DOFs (host " << static_cast<long>(head[5]) << ", "
               << static_cast<long>(head[6]) << "), " << quench_steps
               << " quenched steps (host " << static_cast<long>(head[7]) << ")";
            Check(mpi, static_cast<long>(head[5]) == fhits && static_cast<long>(head[6]) == fdofs
                  && static_cast<long>(head[7]) == quench_steps
                  && (floor == 0.0 || !relax || fdofs > 0), os.str());
         }
         else { num_tests += 4; }
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
