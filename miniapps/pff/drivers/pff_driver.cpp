// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// pff_driver.cpp — explicit dynamic phase-field fracture (phases C1-C5).
//
//   mpirun -np N ./pff_driver -c spe_hftc/configs/section_homog_eta0p45_absorb.toml
//                             -fd [-tf END_TIME] [-o DIR] [-pi ITERS] [-pv INTERVAL]
//
// Phase C1 (this version): elastodynamics with the damage frozen at d = 0
// (-fd/--frozen-damage is required; damage evolution is phase C2). One step:
// central difference with lumped mass (solver/central_difference.hpp), split
// internal force (kernels/split_force_kernels.hpp), Lysmer dashpots and the
// ElkPulseLoadExperiment borehole pressure (domain/explicit_elasticity_operator.hpp).
//
// MOOSE parity (farms@cdms pulsepower/spe_hftc/code_{2d,3d}/elasticity_absorbing.i):
//   time windows of Tw = pulse_duration: zeta = zeta_pulse for the first t_free
//   of each window, zeta_quench after (zeta_fn); the step cap is dt_max before
//   t_free and dt_quiet after (dt_limit); pulse onsets k Tw and the t_free
//   switches are hit exactly (pulse_sync). The stable step is
//   cfl * CriticalStep(omega_max, zeta) (solver/cfl.hpp), omega_max^2 from a
//   power iteration on M^{-1} K at start.
//   One node at [mesh] corner is pinned in all components (fix_corner_*),
//   [boundaries] fixed_z holds u_z = 0 (3D).
//
// Output <dir>/<name>_c1.csv per step at t^n: t, dt, p_wall (Peak on the
// borehole times Shape(t)), v_max (nodal |v^n|), kinetic (lumped 1/2 v M v),
// elastic (int psi), external_work, damping_work, kv_work (trapezoidal in
// u), energy_error = kinetic + elastic - external_work + damping_work +
// kv_work. MOOSE's FarmsExternalWork is a left-endpoint rule in time with the
// consistent-mass kinetic energy, so the energies compare to O(dt) only.

#include "mfem.hpp"
#include "../config/pff_config_parser.hpp"
#include "../io/gmsh_physical_names.hpp"
#include "../domain/explicit_elasticity_operator.hpp"
#include "../solver/central_difference.hpp"
#include "../solver/cfl.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>

using namespace mfem;
using namespace mfem::pff;

namespace
{

Array<int> Marker(const std::vector<int> &attrs, int max_attr)
{
   Array<int> m(max_attr);
   m = 0;
   for (int a : attrs)
   {
      MFEM_VERIFY(a >= 1 && a <= max_attr, "boundary attribute " << a << " is not in the mesh");
      m[a - 1] = 1;
   }
   return m;
}

/// True dofs (all components) of the vertex nearest to x0, on the ranks that own them.
void NearestVertexTrueDofs(ParFiniteElementSpace &vfes, const std::array<real_t, 3> &x0,
                           Array<int> &tdofs, real_t &dist_out)
{
   ParMesh &mesh = *vfes.GetParMesh();
   const int dim = mesh.Dimension();
   real_t best = std::numeric_limits<real_t>::max();
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      const real_t *x = mesh.GetVertex(v);
      real_t s = 0.0;
      for (int k = 0; k < dim; k++) { s += (x[k] - x0[k]) * (x[k] - x0[k]); }
      best = std::min(best, s);
   }
   real_t gbest = 0.0;
   MPI_Allreduce(&best, &gbest, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, vfes.GetComm());
   dist_out = std::sqrt(gbest);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      const real_t *x = mesh.GetVertex(v);
      real_t s = 0.0;
      for (int k = 0; k < dim; k++) { s += (x[k] - x0[k]) * (x[k] - x0[k]); }
      if (s > gbest) { continue; }
      Array<int> vd;
      vfes.GetVertexVDofs(v, vd);
      for (int ld : vd)
      {
         const int t = vfes.GetLocalTDofNumber(ld);
         if (t >= 0) { tdofs.Append(t); }
      }
   }
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;

   const char *config_file = "";
   const char *output_override = "";
   real_t end_time_override = -1.0;
   bool frozen_damage = false;
   int power_its = 60;
   int pv_interval = 0;
   OptionsParser args(argc, argv);
   args.AddOption(&config_file, "-c", "--config", "Case TOML file.");
   args.AddOption(&output_override, "-o", "--output-dir", "Override [output] directory.");
   args.AddOption(&end_time_override, "-tf", "--end-time", "Override [time] end_time.");
   args.AddOption(&frozen_damage, "-fd", "--frozen-damage", "-nfd", "--no-frozen-damage",
                  "Hold d = 0 (phase C1). Required until damage evolution (C2) exists.");
   args.AddOption(&power_its, "-pi", "--power-iterations",
                  "Power iterations for omega_max^2 of M^-1 K.");
   args.AddOption(&pv_interval, "-pv", "--paraview-interval",
                  "Write ParaView output every N steps (0: none).");
   args.Parse();
   if (!args.Good() || std::string(config_file).empty())
   {
      if (Mpi::Root()) { args.PrintUsage(mfem::out); }
      return 1;
   }

   CaseConfig c = CaseConfigParser::ParseFile(config_file);
   if (std::string(output_override).size()) { c.output_dir = output_override; }
   if (end_time_override > 0.0) { c.end_time = end_time_override; }

   if (!frozen_damage)
   {
      MFEM_ABORT("pff_driver: damage evolution is phase C2 and not implemented yet; "
                 "run with -fd (d frozen at 0)");
   }
   MFEM_VERIFY(c.model == CaseModel::Solid, "pff_driver: hydromech is phase C4");
   MFEM_VERIFY(c.confinement_pressure == 0.0 && c.confinement_bdr.empty(),
               "pff_driver: confinement is not implemented yet (C2)");
   MFEM_VERIFY(c.initial_damage.empty(),
               "pff_driver: seeded initial damage needs the damage solve (C2)");
   if (c.zeta_quench > 0.0 && c.quench != QuenchMode::Stiffness)
   {
      MFEM_ABORT("pff_driver: quench modes 'mass' and 'implicit' are phase C3");
   }

   // Mesh (path relative to the config file).
   const std::filesystem::path cfg_dir = std::filesystem::path(config_file).parent_path();
   const std::string mesh_path = (cfg_dir / c.mesh_file).string();
   PhysicalNameMap names = ReadPhysicalNames(mesh_path);
   Mesh serial(mesh_path.c_str(), 1, 1);
   MFEM_VERIFY(serial.Dimension() == c.dim, "mesh dimension " << serial.Dimension()
               << " != [mesh] dim " << c.dim);
   for (int r = 0; r < c.uniform_refinements; r++) { serial.UniformRefinement(); }
   const int dim = serial.Dimension();
   const int max_attr = serial.bdr_attributes.Max();
   ParMesh pmesh(comm, serial);
   serial.Clear();

   H1_FECollection h1(1, dim);
   ParFiniteElementSpace Vfes(&pmesh, &h1, dim, Ordering::byNODES);
   ParFiniteElementSpace Sfes(&pmesh, &h1);

   // Fixed dofs: the corner node (all components) and fixed_z faces.
   Array<int> fixed;
   real_t corner_dist = 0.0;
   NearestVertexTrueDofs(Vfes, c.corner, fixed, corner_dist);
   if (!c.fixed_z_bdr.empty())
   {
      MFEM_VERIFY(dim == 3, "fixed_z needs a 3D mesh");
      Array<int> list;
      Vfes.GetEssentialTrueDofs(Marker(ResolveAttributes(c.fixed_z_bdr, names, dim - 1), max_attr),
                                list, 2);
      fixed.Append(list);
   }
   fixed.Sort();
   fixed.Unique();

   ExplicitElasticityParams ep;
   ep.K = c.K();
   ep.G = c.G();
   ep.density = c.density;
   ep.decomposition = c.decomposition;
   ep.deg = c.pf.deg;
   ep.pressure_marker = Marker(ResolveAttributes(c.pressure_bdr, names, dim - 1), max_attr);
   ep.absorbing_marker = Marker(ResolveAttributes(c.absorbing_bdr, names, dim - 1), max_attr);
   PulseLoad pulse(c.pulse);
   ExplicitElasticityOperator op(Vfes, Sfes, ep, pulse, fixed);

   const real_t omega2 = op.MaxEigenvalue(power_its);
   const real_t omega = std::sqrt(omega2);
   const real_t Tw = c.pulse.pulse_duration;
   const long long global_ne = pmesh.GetGlobalNE();
   const HYPRE_BigInt global_nodes = Sfes.GlobalTrueVSize();
   if (Mpi::Root())
   {
      mfem::out << "pff_driver (C1, d frozen at 0): " << c.name << "\n  source: " << c.source
                << "\n  mesh " << mesh_path << ": " << global_ne << " elements, "
                << global_nodes << " nodes"
                << "\n  c_p = " << ep.Cp() << " m/s, c_s = " << ep.Cs() << " m/s"
                << "\n  corner pin at distance " << corner_dist << " from the [mesh] corner"
                << "\n  omega_max = " << omega << " rad/s (" << power_its
                << " power iterations), dt_crit(zeta=0) = " << CriticalStep(omega, 0.0)
                << " s\n  peak pressure on the borehole " << op.PeakOnPressureBoundary()
                << " Pa, end time " << c.end_time << " s\n";
   }

   // State.
   CentralDifferenceState S;
   const int n = Vfes.GetTrueVSize();
   S.u.SetSize(n); S.u = 0.0;
   S.v_half.SetSize(n); S.v_half = 0.0;
   S.F_prev.SetSize(n); S.F_prev = 0.0;
   Vector D(Sfes.GetTrueVSize());
   D = 0.0;
   CentralDifference cd(op);
   CentralDifferenceStepInfo info;

   std::ofstream csv;
   if (Mpi::Root())
   {
      std::filesystem::create_directories(c.output_dir);
      const std::string path = c.output_dir + "/" + c.name + "_c1.csv";
      csv.open(path);
      MFEM_VERIFY(csv.good(), "cannot write " << path);
      csv << "step,t,dt,p_wall,v_max,kinetic,elastic,external_work,damping_work,kv_work,"
          "energy_error,shadow_energy,shadow_input,shadow_dissipated,shadow_error\n"
          << std::setprecision(12);
   }
   std::unique_ptr<ParaViewDataCollection> pv;
   ParGridFunction vel(&Vfes), disp(&Vfes);
   if (pv_interval > 0)
   {
      pv.reset(new ParaViewDataCollection(c.name + "_c1", &pmesh));
      pv->SetPrefixPath(c.output_dir);
      pv->RegisterField("velocity", &vel);
      pv->RegisterField("displacement", &disp);
   }

   // Work accounting (trapezoidal in u).
   Vector u_prev(S.u), fext_prev(n), fdamp_prev(n), fkv_prev(n);
   fext_prev = 0.0; fdamp_prev = 0.0; fkv_prev = 0.0;
   real_t W_ext = 0.0, W_damp = 0.0, W_kv = 0.0;
   // Discrete (shadow) energy of the central difference scheme, conserved up to
   // round-off for linear K and constant dt:
   //   H^{n-1/2} = 1/2 v^{n-1/2} M v^{n-1/2} + 1/2 u^n . F_int(u^{n-1}),
   //   H^{n+1/2} - H^{n-1/2} = (f_ext^n - f_damp^n - f_kv^n) . (u^{n+1} - u^{n-1}) / 2.
   // shadow_error^n = H^{n-1/2} - H^{1/2} - sum_{k=1}^{n-1} (...) (u^{k+1} - u^{k-1}) / 2.
   Vector u_nm1(S.u);
   real_t H_ref = 0.0, Ws_in = 0.0, Ws_out = 0.0, shadow_err_report = 0.0;
   const real_t eps_t = 1e-12 * Tw;

   auto window = [&](real_t t) { return std::fmod(t + eps_t, Tw); };
   auto zeta_at = [&](real_t t) { return window(t) < c.t_free ? c.zeta_pulse : c.zeta_quench; };
   auto cap_at = [&](real_t t) { return window(t) < c.t_free ? c.dt_max : c.dt_quiet; };
   auto next_event = [&](real_t t)
   {
      const real_t k = std::floor((t + eps_t) / Tw);
      real_t ev = (window(t) < c.t_free) ? k * Tw + c.t_free : (k + 1.0) * Tw;
      return std::min(ev, c.end_time);
   };

   tic_toc.Clear();
   tic_toc.Start();
   while (S.t < c.end_time - eps_t)
   {
      const real_t zeta = zeta_at(S.t);
      real_t dt = c.cfl * CriticalStep(omega, zeta);
      if (cap_at(S.t) > 0.0) { dt = std::min(dt, cap_at(S.t)); }
      const real_t ev = next_event(S.t);
      if (S.t + dt > ev - 1e-6 * dt) { dt = ev - S.t; }
      const real_t s = pulse.Shape(S.t);
      const real_t t_n = S.t;
      const long step_n = S.step;
      real_t H = 0.0;
      {
         real_t loc = 0.5 * InnerProduct(S.u, S.F_prev), glob = 0.0;
         MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
         H = op.KineticEnergy(S.v_half) + glob;
      }
      if (step_n == 1) { H_ref = H; }
      const Vector u_n(S.u);
      cd.Step(S, D, s, zeta, dt, info);
      if (step_n >= 1)
      {
         real_t li = 0.0, lo = 0.0;
         for (int i = 0; i < n; i++)
         {
            const real_t du2 = 0.5 * (S.u(i) - u_nm1(i));
            li += info.f_ext(i) * du2;
            lo += (info.f_damp(i) + info.f_kv(i)) * du2;
         }
         real_t loc2[2] = {li, lo}, glob2[2];
         MPI_Allreduce(loc2, glob2, 2, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
         // The step-n increment belongs to H^{n+1/2}: report the error before adding it.
         shadow_err_report = (step_n >= 1) ? H - H_ref - Ws_in + Ws_out : 0.0;
         Ws_in += glob2[0];
         Ws_out += glob2[1];
      }
      u_nm1 = u_n;

      // Works up to t^n (info refers to t^n).
      real_t dWe = 0.0, dWd = 0.0, dWk = 0.0;
      {
         // u^n = S.u - dt v^{n+1/2}
         Vector un(S.u);
         un.Add(-dt, S.v_half);
         for (int i = 0; i < n; i++)
         {
            const real_t d_u = un(i) - u_prev(i);
            dWe += 0.5 * (fext_prev(i) + info.f_ext(i)) * d_u;
            dWd += 0.5 * (fdamp_prev(i) + info.f_damp(i)) * d_u;
            dWk += 0.5 * (fkv_prev(i) + info.f_kv(i)) * d_u;
         }
         u_prev = un;
      }
      real_t loc[3] = {dWe, dWd, dWk}, glob[3];
      MPI_Allreduce(loc, glob, 3, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
      W_ext += glob[0]; W_damp += glob[1]; W_kv += glob[2];
      fext_prev = info.f_ext; fdamp_prev = info.f_damp; fkv_prev = info.f_kv;

      const real_t ke = op.KineticEnergy(info.v_n);
      const real_t vmax = op.MaxNodalSpeed(info.v_n);
      const real_t err = ke + info.psi - W_ext + W_damp + W_kv;
      if (Mpi::Root())
      {
         csv << step_n << ',' << t_n << ',' << dt << ',' << op.PeakOnPressureBoundary() * s << ','
             << vmax << ',' << ke << ',' << info.psi << ',' << W_ext << ',' << W_damp << ','
             << W_kv << ',' << err << ',' << H << ',' << Ws_in << ',' << Ws_out << ','
             << shadow_err_report << '\n';
         if (step_n % 200 == 0)
         {
            mfem::out << "step " << step_n << " t = " << t_n << " dt = " << dt << " v_max = "
                      << vmax << " KE = " << ke << " EE = " << info.psi << " W_ext = " << W_ext
                      << " W_damp = " << W_damp << " err = " << err << " shadow err = "
                      << shadow_err_report << "\n";
         }
      }
      if (pv && step_n % pv_interval == 0)
      {
         vel.SetFromTrueDofs(info.v_n);
         Vector un(S.u);
         un.Add(-dt, S.v_half);
         disp.SetFromTrueDofs(un);
         pv->SetCycle(static_cast<int>(step_n));
         pv->SetTime(t_n);
         pv->Save();
      }
   }
   tic_toc.Stop();
   if (Mpi::Root())
   {
      mfem::out << "done: " << S.step << " steps to t = " << S.t << " in " << tic_toc.RealTime()
                << " s (" << tic_toc.RealTime() / std::max<long>(S.step, 1) << " s/step)\n";
   }
   return 0;
}
