// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// pff_driver.cpp — explicit dynamic phase-field fracture (phases C1-C5).
//
//   mpirun -np N ./pff_driver -c spe_hftc/configs/section_homog_eta0p45_absorb.toml
//                             [-fd] [-tf END_TIME] [-o DIR] [-pi ITERS] [-pv INTERVAL] [-pm PARTS_DIR]
//   -pm: read the mesh as per-rank parts from pff_partition (large 3D meshes),
//        instead of loading the serial mesh on every rank.
//
// Phase C2 (this version): explicit elastodynamics plus one bound-constrained
// damage solve per step (below); -fd/--frozen-damage holds d = 0 (phase C1,
// output <name>_c1.csv, otherwise <name>_c2.csv). One step:
// central difference with consistent or lumped mass ([time] mass; solver/central_difference.hpp), split
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
// C2 columns: crack_energy (Gc/(c0 l) int alpha(d) + 1/2 d . kappa K d; its change since t^1, part
// of energy_error), d_max, fracture_extent and extent_090 (largest xy
// distance from [pulse] center of a node with d >= [metrics] d_threshold and
// d >= 0.9, minus r_hole), damage_qp_its (0: skipped step). The shadow_*
// columns assume a fixed d and are diagnostics of C1 only.
// quench = "mass" (C3, a pff choice approved by Andy 2026-11-22, not in the
// MOOSE inputs): after t_free of each window the stiffness quench is replaced
// by mass-proportional damping C = alpha_mass M (lumped), which leaves the
// stable step unchanged; mass_damping_work (part of energy_error) is its work.
// [confinement] (cmame): constant pressure on [boundaries] confinement from
// t = 0 (MOOSE Pressure BC static_pressure_outer, no preload); its work is
// confinement_work (external_work is the pulse only), part of energy_error.
//
// Phase C4, [case] model = "hydromech" (farms@cdms pulsepower/cmame_revision/
// 2d_hydromech/benchmark/elasticity_E1d25.i + fracture_E1d25.i): u-p Biot with
// the damage-dependent properties of multiphysics/poroelastic. One step
// (PLAN.md, single pass): u^{n+1} explicitly with the pore-pressure force
// -int alpha(d^n) p^n div(phi); p^{n+1} by backward Euler with eps_v(u^{n+1})
// (pore_pressure_operator.hpp); the history H = max(H, <2 psi+ + c_p p^2>_+)
// at the quadrature points (ElkPorousFlowHistoryEnergyEnhanced), whose element
// mean (the CONSTANT MONOMIAL psie_active_enhanced MOOSE transfers) drives
// d^{n+1} as 1/2 H. MOOSE integrates u with implicit Newmark (beta 1/4) at
// dt = 1e-8 and iterates the damage coupling (fixed_point_max_its = 10); pff
// is explicit with the step from the undrained stiffness bound
// (UndrainedStiffnessRatio). [confinement] initial_state = "static": u^0 is
// the static equilibrium under the confinement and the initial pore pressure
// (MOOSE static_solve.i read by SolutionIC), d = 0.
// Extra columns (hydromech): fluid_energy (1/2 int p^2/M), fluid_dissipated
// (int int grad p . k/mu grad p dt), fluid_boundary_work (work of the drained
// boundary on the fluid), p_min, p_max, pp_its. energy_error then is the
// change since t^1 of kinetic + elastic + crack + fluid_energy +
// fluid_dissipated - the works (the seeded damage and the preload energy are
// in the t^1 reference).

#include "mfem.hpp"
#include "../config/pff_config_parser.hpp"
#include "../io/gmsh_physical_names.hpp"
#include "../domain/explicit_elasticity_operator.hpp"
#include "../solver/central_difference.hpp"
#include "../solver/cfl.hpp"
#include "../domain/qs_damage_operator.hpp"
#include "../solver/bound_constrained_qp.hpp"
#include "../phasefield/initial_damage.hpp"
#include "../multiphysics/poroelastic/pore_pressure_operator.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
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

/// Inactive-set residual norm at x of the bound-constrained damage problem,
/// with the active set of BoundConstrainedQP (and of PETSc vinewtonrsls).
real_t InactiveResidualNorm(MPI_Comm comm, const Vector &x, const Vector &r, const Vector &lo,
                            const Vector &hi)
{
   real_t s = 0.0;
   for (int i = 0; i < x.Size(); i++)
   {
      const bool active = (x(i) <= lo(i) && r(i) > 0.0) || (x(i) >= hi(i) && r(i) < 0.0);
      if (!active) { s += r(i) * r(i); }
   }
   real_t g = 0.0;
   MPI_Allreduce(&s, &g, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   return std::sqrt(g);
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
   const char *mass_override = "";
   int power_its = 60;
   int pv_interval = 0;
   const char *parts_dir = "";
   OptionsParser args(argc, argv);
   args.AddOption(&config_file, "-c", "--config", "Case TOML file.");
   args.AddOption(&output_override, "-o", "--output-dir", "Override [output] directory.");
   args.AddOption(&end_time_override, "-tf", "--end-time", "Override [time] end_time.");
   args.AddOption(&frozen_damage, "-fd", "--frozen-damage", "-nfd", "--no-frozen-damage",
                  "Hold d = 0 (phase C1: elastodynamics only).");
   args.AddOption(&power_its, "-pi", "--power-iterations",
                  "Power iterations for omega_max^2 of M^-1 K.");
   args.AddOption(&mass_override, "-mass", "--mass", "Override [time] mass: consistent or lumped.");
   args.AddOption(&pv_interval, "-pv", "--paraview-interval",
                  "Write ParaView output every N steps (0: none).");
   args.AddOption(&parts_dir, "-pm", "--partitioned-mesh",
                  "Read the [mesh] file as per-rank parts written by pff_partition (one per MPI rank).");
   args.Parse();
   if (!args.Good() || std::string(config_file).empty())
   {
      if (Mpi::Root()) { args.PrintUsage(mfem::out); }
      return 1;
   }

   CaseConfig c = CaseConfigParser::ParseFile(config_file);
   if (std::string(output_override).size()) { c.output_dir = output_override; }
   if (end_time_override > 0.0) { c.end_time = end_time_override; }
   if (std::string(mass_override) == "consistent") { c.consistent_mass = true; }
   else if (std::string(mass_override) == "lumped") { c.consistent_mass = false; }
   else { MFEM_VERIFY(std::string(mass_override).empty(), "-mass must be consistent or lumped"); }

   const bool hydromech = c.model == CaseModel::Hydromech;
   MFEM_VERIFY(!hydromech || (c.zeta_pulse == 0.0 && c.zeta_quench == 0.0),
               "pff_driver: hydromech has no stiffness damping (zeta_pulse = zeta_quench = 0)");
   MFEM_VERIFY(hydromech || c.drained_bdr.empty(), "pff_driver: [boundaries] drained needs model = hydromech");
   MFEM_VERIFY((c.confinement_pressure == 0.0) == c.confinement_bdr.empty(),
               "pff_driver: [confinement] pressure and [boundaries] confinement go together");
   MFEM_VERIFY(frozen_damage || (c.pf.driving == DrivingEnergy::HalfHistoryEnhanced) == hydromech,
               "pff_driver: driving = half_history_enhanced goes with model = hydromech, active with solid");
   MFEM_VERIFY(frozen_damage || c.coupling == CouplingScheme::SinglePass,
               "pff_driver: [coupling] fixed_point is not implemented for the explicit step");
   MFEM_VERIFY(c.quench != QuenchMode::Implicit, "pff_driver: quench = 'implicit' is not implemented");
   MFEM_VERIFY(c.quench != QuenchMode::Mass || c.zeta_quench == 0.0,
               "pff_driver: quench = 'mass' uses alpha_mass; set zeta_quench = 0");

   // Mesh (path relative to the config file).
   const std::filesystem::path cfg_dir = std::filesystem::path(config_file).parent_path();
   const std::string mesh_path = (cfg_dir / c.mesh_file).string();
   PhysicalNameMap names;
   std::unique_ptr<ParMesh> pmesh_ptr;
   int max_attr = 0;
   if (std::string(parts_dir).empty())
   {
      // every rank loads the serial mesh, ParMesh partitions it (METIS)
      names = ReadPhysicalNames(mesh_path);
      Mesh serial(mesh_path.c_str(), 1, 1);
      MFEM_VERIFY(serial.Dimension() == c.dim, "mesh dimension " << serial.Dimension()
                  << " != [mesh] dim " << c.dim);
      for (int r = 0; r < c.uniform_refinements; r++) { serial.UniformRefinement(); }
      max_attr = serial.bdr_attributes.Max();
      pmesh_ptr.reset(new ParMesh(comm, serial));
   }
   else
   {
      // each rank reads its own part (pff_partition, same METIS partitioning)
      namespace fs = std::filesystem;
      const fs::path pdir(parts_dir);
      std::ifstream info(pdir / "info.txt");
      MFEM_VERIFY(info, "pff_driver: no " << (pdir / "info.txt") << " (run pff_partition)");
      std::map<std::string, std::string> kv;
      for (std::string k, v; info >> k && std::getline(info >> std::ws, v);) { kv[k] = v; }
      MFEM_VERIFY(std::stoi(kv["np"]) == Mpi::WorldSize(), "pff_driver: " << pdir << " has "
                  << kv["np"] << " parts, the run has " << Mpi::WorldSize() << " ranks");
      MFEM_VERIFY(fs::path(kv["source"]).filename() == fs::path(mesh_path).filename(),
                  "pff_driver: parts in " << pdir << " are of " << kv["source"]
                  << ", the config's [mesh] file is " << mesh_path);
      names = ReadPhysicalNames((pdir / "physical_names.msh").string());
      max_attr = std::stoi(kv["max_bdr_attribute"]);
      std::ifstream part(pdir / MakeParFilename("mesh.", Mpi::WorldRank()));
      MFEM_VERIFY(part, "pff_driver: cannot open part " << Mpi::WorldRank() << " in " << pdir);
      pmesh_ptr.reset(new ParMesh(comm, part));
      MFEM_VERIFY(pmesh_ptr->Dimension() == c.dim, "mesh dimension " << pmesh_ptr->Dimension()
                  << " != [mesh] dim " << c.dim);
      for (int r = 0; r < c.uniform_refinements; r++) { pmesh_ptr->UniformRefinement(); }
   }
   ParMesh &pmesh = *pmesh_ptr;
   const int dim = pmesh.Dimension();

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
   for (int comp = 0; comp < 2; comp++)
   {
      const std::vector<std::string> &bdr = comp == 0 ? c.fixed_x_bdr : c.fixed_y_bdr;
      if (bdr.empty()) { continue; }
      Array<int> list;
      Vfes.GetEssentialTrueDofs(Marker(ResolveAttributes(bdr, names, dim - 1), max_attr), list, comp);
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
   ep.consistent_mass = c.consistent_mass;
   ep.pressure_marker = Marker(ResolveAttributes(c.pressure_bdr, names, dim - 1), max_attr);
   ep.absorbing_marker = Marker(ResolveAttributes(c.absorbing_bdr, names, dim - 1), max_attr);
   ep.confinement_marker = Marker(ResolveAttributes(c.confinement_bdr, names, dim - 1), max_attr);
   ep.confinement_pressure = c.confinement_pressure;
   PulseLoad pulse(c.pulse);
   ExplicitElasticityOperator op(Vfes, Sfes, ep, pulse, fixed);

   const real_t omega2_drained = op.MaxEigenvalue(power_its);
   // hydromech: the undrained stiffness bounds the explicit step
   const real_t undrained_ratio = hydromech ? UndrainedStiffnessRatio(c.poro, c.G(), dim) : 1.0;
   const real_t omega2 = omega2_drained * undrained_ratio;
   const real_t omega = std::sqrt(omega2);
   const real_t Tw = c.pulse.pulse_duration;
   const long long global_ne = pmesh.GetGlobalNE();
   const HYPRE_BigInt global_nodes = Sfes.GlobalTrueVSize();
   if (Mpi::Root())
   {
      mfem::out << "pff_driver (" << (frozen_damage ? "C1, d frozen at 0" : "C2, damage")
                << "): " << c.name << "\n  source: " << c.source
                << "\n  mesh " << mesh_path << ": " << global_ne << " elements, "
                << global_nodes << " nodes"
                << "\n  c_p = " << ep.Cp() << " m/s, c_s = " << ep.Cs() << " m/s"
                << "\n  corner pin at distance " << corner_dist << " from the [mesh] corner"
                << "\n  omega_max = " << omega << " rad/s (" << power_its
                << " power iterations), dt_crit(zeta=0) = " << CriticalStep(omega, 0.0)
                << " s\n  peak pressure on the borehole " << op.PeakOnPressureBoundary()
                << " Pa, end time " << c.end_time << " s\n";
      if (hydromech)
      {
         const PoroPoint q0 = EvalPoroPoint(c.poro, 0.0);
         mfem::out << "  hydromech: alpha(0) = " << q0.alpha << ", M(0) = " << 1.0 / q0.inv_M
                   << " Pa, k0/mu = " << q0.mobility << " m^2/(Pa s), undrained stiffness ratio "
                   << undrained_ratio << " (in omega_max above)\n";
      }
   }

   // State.
   CentralDifferenceState S;
   const int n = Vfes.GetTrueVSize();
   S.u.SetSize(n); S.u = 0.0;
   S.v_half.SetSize(n); S.v_half = 0.0;
   S.F_prev.SetSize(n); S.F_prev = 0.0;
   Vector D(Sfes.GetTrueVSize());
   D = 0.0;

   // Damage (phase C2): farms@cdms pulsepower/spe_hftc/code_{2d,3d}/fracture.i.
   // One damage solve per explicit step (single pass), with psie_active of
   // u^{n+1} as the CONSTANT MONOMIAL driving field. Bounds: farms
   // VariableConstantIrreversibleBounds (bound_value = initial_damage_aux at
   // the first step, d_old after) and ConstantBounds 1, i.e. lo = max(d_old,
   // seed) with d^0 = 0. Tolerances of fracture.i: nl_rel_tol 1e-8,
   // nl_abs_tol 1e-10 with automatic_scaling (compute_scaling_once, Jacobian
   // diagonal): the absolute tolerance is on s |r| with s = 1/max_i A_ii of
   // the first solve. A step whose start point d_old already meets the
   // absolute tolerance is skipped, as PETSc stops at iteration 0 there.
   const IntegrationRule &dir = IntRules.Get(pmesh.GetTypicalElementGeometry(),
                                             ep.quadrature_order);
   L2_FECollection l2(0, dim);
   ParFiniteElementSpace Lfes(&pmesh, &l2);
   ParGridFunction Dgf(&Lfes);
   Dgf = 0.0;
   std::unique_ptr<QsDamageOperator> damage;
   BoundConstrainedQP qp;
   qp.rel_tol = 1e-8;
   real_t damage_abs_tol = -1.0;   // set at the first damage step
   const real_t fracture_nl_abs_tol = 1e-10;
   Vector seed, lo(D.Size()), hi(D.Size()), rd(D.Size()), dx(D.Size()), psi_el;
   hi = 1.0;
   Vector m_lump(D.Size());   // int phi_i
   std::unique_ptr<HypreParMatrix> Mcons;
   Vector node_r(D.Size());   // distance of each node from the pulse center in the xy plane
   long damage_solves = 0, damage_qp_its = 0;
   if (!frozen_damage)
   {
      damage.reset(new QsDamageOperator(Sfes, c.pf, &dir));
      SeedInitialDamage(Sfes, c.initial_damage, seed);
      ParLinearForm one(&Sfes);
      ConstantCoefficient onec(1.0);
      auto *li = new DomainLFIntegrator(onec);
      li->SetIntRule(&dir);
      one.AddDomainIntegrator(li);
      one.Assemble();
      one.ParallelAssemble(m_lump);
      ParBilinearForm mf(&Sfes);
      auto *mi = new MassIntegrator();
      mi->SetIntRule(&dir);
      mf.AddDomainIntegrator(mi);
      mf.Assemble();
      mf.Finalize();
      Mcons.reset(mf.ParallelAssemble());
      ParGridFunction rg(&Sfes);
      const real_t ctr[2] = {c.pulse.center[0], c.pulse.center[1]};
      FunctionCoefficient rf([&ctr](const Vector &x) { return std::hypot(x(0) - ctr[0], x(1) - ctr[1]); });
      rg.ProjectCoefficient(rf);
      rg.GetTrueDofs(node_r);
      if (Mpi::Root())
      {
         mfem::out << "  damage: " << (c.pf.crack.model == CrackModel::AT1 ? "AT1" : "AT2")
                   << ", Gc = " << c.pf.Gc << " J/m^2, l = " << c.pf.l << " m, kappa = "
                   << c.pf.DiffusionCoefficient() << ", onset psi+ = "
                   << c.pf.OnsetDrivingEnergy() << " J/m^3\n";
      }
   }
   // Hydromech state: pore pressure (true dofs), drained dofs, history at the points.
   std::unique_ptr<PorePressureOperator> poro;
   Vector Pp, drained_peak, Hq, f_poro;
   Array<int> drained;
   PorePressureOperator::StepInfo pinfo;
   real_t Ef_n = 0.0, W_fdiss = 0.0, W_fbdry = 0.0;
   long pp_its_total = 0;
   const real_t c_hist = hydromech ? HistoryPressureCoefficient(c.poro) : 0.0;
   if (hydromech)
   {
      if (!c.drained_bdr.empty())
      {
         Sfes.GetEssentialTrueDofs(Marker(ResolveAttributes(c.drained_bdr, names, dim - 1), max_attr),
                                   drained);
      }
      poro.reset(new PorePressureOperator(Vfes, Sfes, c.poro, drained, dir,
                                          c.poro_amg ? PoroPreconditioner::AMG : PoroPreconditioner::Jacobi));
      Pp.SetSize(Sfes.GetTrueVSize());
      Pp = c.initial_pore_pressure;
      // drained value: 0, or the pulse pressure Peak(x) Shape(t) (FunctionDirichletBC pp)
      drained_peak.SetSize(Sfes.GetTrueVSize());
      drained_peak = 0.0;
      if (c.drained_pulse)
      {
         ParGridFunction pk(&Sfes);
         FunctionCoefficient pf([&pulse](const Vector &x)
         {
            real_t X[3] = {0.0, 0.0, 0.0};
            for (int k = 0; k < x.Size(); k++) { X[k] = x(k); }
            return pulse.Peak(X);
         });
         pk.ProjectCoefficient(pf);
         pk.GetTrueDofs(drained_peak);
      }
      for (int t : drained) { Pp(t) = drained_peak(t) * pulse.Shape(0.0); }
      Hq.SetSize(pmesh.GetNE() * dir.GetNPoints());
      Hq = 0.0;
      f_poro.SetSize(n);
      real_t rate = 0.0;
      poro->Energies(D, Pp, Ef_n, rate);
      if (c.static_initial_state)
      {
         // K u0 = b_c - F_p(p0) with the fixed dofs held, d = 0 (linear): CG on F_int.
         struct KOp : public Operator
         {
            const ExplicitElasticityOperator &o; const Vector &D0;
            KOp(const ExplicitElasticityOperator &o_, const Vector &D0_, int n_) : Operator(n_), o(o_), D0(D0_) {}
            void Mult(const Vector &x, Vector &y) const override
            {
               Vector xx(x);
               for (int t : o.FixedTrueDofs()) { xx(t) = 0.0; }
               o.InternalForce(xx, D0, y);
               for (int t : o.FixedTrueDofs()) { y(t) = x(t); }
            }
         } Kop(op, D, n);
         Vector rhs(op.ConfinementVector());
         f_poro = 0.0;
         poro->AddCouplingForce(D, Pp, f_poro);
         rhs -= f_poro;
         for (int t : fixed) { rhs(t) = 0.0; }
         CGSolver cg(comm);
         cg.SetOperator(Kop);
         cg.SetRelTol(1e-12);
         cg.SetAbsTol(0.0);
         cg.SetMaxIter(200000);
         cg.SetPrintLevel(0);
         S.u = 0.0;
         cg.Mult(rhs, S.u);
         MFEM_VERIFY(cg.GetConverged(), "static initial state: CG did not converge ("
                     << cg.GetNumIterations() << " iterations, |r| = " << cg.GetFinalNorm() << ")");
         real_t psi0 = 0.0;
         Vector F0;
         op.InternalForce(S.u, D, F0, &psi0);
         S.F_prev = F0;
         if (Mpi::Root())
         {
            mfem::out << "  static initial state: CG " << cg.GetNumIterations()
                      << " iterations, elastic energy " << psi0 << "\n";
         }
      }
   }
   else
   {
      MFEM_VERIFY(!c.static_initial_state,
                  "pff_driver: [confinement] initial_state = 'static' is wired for hydromech only");
   }
   // Crack energy Gc/(c0 l) int alpha(d) + 1/2 d . kappa K d of the damage d.
   auto crack_energy = [&](const Vector &d)
   {
      if (frozen_damage) { return 0.0; }
      Vector Kd(d.Size());
      damage->DiffusionMatrix().Mult(d, Kd);
      real_t a = 0.0;
      if (c.pf.crack.model == CrackModel::AT1) { a = InnerProduct(comm, m_lump, d); }
      else
      {
         Vector Md(d.Size());
         Mcons->Mult(d, Md);
         a = InnerProduct(comm, d, Md);
      }
      return c.pf.Gc / (c.pf.crack.c0() * c.pf.l) * a + 0.5 * InnerProduct(comm, d, Kd);
   };
   // Largest xy distance from the pulse center of a node with d >= thr, minus r_hole (0 if none).
   auto extent = [&](const Vector &d, real_t thr)
   {
      real_t loc = 0.0, glob = 0.0;
      for (int i = 0; i < d.Size(); i++) { if (d(i) >= thr) { loc = std::max(loc, node_r(i) - c.r_hole); } }
      MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
      return glob;
   };
   CentralDifference cd(op);
   CentralDifferenceStepInfo info;

   std::ofstream csv;
   if (Mpi::Root())
   {
      std::filesystem::create_directories(c.output_dir);
      const std::string path = c.output_dir + "/" + c.name + (frozen_damage ? "_c1.csv" : "_c2.csv");
      csv.open(path);
      MFEM_VERIFY(csv.good(), "cannot write " << path);
      csv << "step,t,dt,p_wall,v_max,kinetic,elastic,external_work,damping_work,kv_work,"
          "energy_error,shadow_energy,shadow_input,shadow_dissipated,shadow_error,"
          "crack_energy,d_max,fracture_extent,extent_090,damage_qp_its,mass_damping_work,confinement_work"
          << (hydromech ? ",fluid_energy,fluid_dissipated,fluid_boundary_work,p_min,p_max,pp_its" : "")
          << "\n" << std::setprecision(12);
   }
   std::unique_ptr<ParaViewDataCollection> pv;
   ParGridFunction vel(&Vfes), disp(&Vfes), dout(&Sfes), ppout(&Sfes);
   if (pv_interval > 0)
   {
      pv.reset(new ParaViewDataCollection(c.name + "_c1", &pmesh));
      pv->SetPrefixPath(c.output_dir);
      pv->RegisterField("velocity", &vel);
      pv->RegisterField("displacement", &disp);
      pv->RegisterField("d", &dout);
      if (hydromech) { pv->RegisterField("pp", &ppout); }
   }

   // Work accounting (trapezoidal in u).
   Vector u_prev(S.u), fext_prev(n), fdamp_prev(n), fkv_prev(n), fmass_prev(n);
   fext_prev = 0.0; fdamp_prev = 0.0; fkv_prev = 0.0; fmass_prev = 0.0;
   const Vector &f_conf = op.ConfinementVector();   // constant in time
   real_t W_ext = 0.0, W_damp = 0.0, W_kv = 0.0, W_mass = 0.0, W_conf = 0.0;
   real_t Ec_first = 0.0, err_ref = 0.0;
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
   auto alpha_at = [&](real_t t)
   { return (c.quench == QuenchMode::Mass && window(t) >= c.t_free) ? c.alpha_mass : 0.0; };
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
      const Vector d_n(D);
      // hydromech: pore-pressure force of p^n at d^n; the values at t^n for the output
      const real_t Ef_out = Ef_n, Wfd_out = W_fdiss, Wfb_out = W_fbdry;
      real_t p_rng[2] = {0.0, 0.0};
      Vector p_n;
      if (hydromech)
      {
         f_poro = 0.0;
         poro->AddCouplingForce(D, Pp, f_poro);
         real_t loc[2] = {Pp.Size() ? -Pp.Min() : -std::numeric_limits<real_t>::max(),
                          Pp.Size() ? Pp.Max() : -std::numeric_limits<real_t>::max()};
         MPI_Allreduce(loc, p_rng, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
         p_rng[0] = -p_rng[0];
         if (pv && step_n % pv_interval == 0) { p_n = Pp; }
      }
      cd.Step(S, D, s, zeta, dt, info, alpha_at(S.t), hydromech ? &f_poro : nullptr);
      if (hydromech)
      {
         // p^{n+1} with eps_v(u^{n+1}) - eps_v(u^n), properties at d^n
         for (int t : drained) { Pp(t) = drained_peak(t) * pulse.Shape(S.t); }
         poro->Step(D, S.u, u_n, dt, Pp, pinfo);
         pp_its_total += pinfo.iterations;
         Ef_n = pinfo.fluid_energy;
         W_fdiss += dt * pinfo.dissipation_rate;
         W_fbdry += dt * pinfo.boundary_power;
      }

      // Damage d^{n+1} from psie_active(u^{n+1}).
      int qp_its_step = 0;
      if (!frozen_damage)
      {
         op.Force().ElementActiveEnergy(S.u, psi_el);
         psi_el.HostRead();
         if (hydromech)
         {
            // H at the quadrature points with p^{n+1} there; its element mean drives d
            // (psi+ is the element value: constant on P1 simplices)
            ParGridFunction pg(&Sfes);
            pg.SetFromTrueDofs(Pp);
            const int nq = dir.GetNPoints();
            for (int e = 0; e < pmesh.GetNE(); e++)
            {
               ElementTransformation &T = *pmesh.GetElementTransformation(e);
               real_t sum = 0.0, vol = 0.0;
               for (int q = 0; q < nq; q++)
               {
                  const IntegrationPoint &ip = dir.IntPoint(q);
                  T.SetIntPoint(&ip);
                  const real_t pq = pg.GetValue(T, ip);
                  real_t &H = Hq(e * nq + q);
                  H = std::max(std::max(2.0 * psi_el(e) + c_hist * pq * pq, real_t(0.0)), H);
                  const real_t w = ip.weight * T.Weight();
                  sum += w * H;
                  vol += w;
               }
               Dgf(e) = c.pf.DrivingFrom(sum / vol);
            }
         }
         else
         {
            for (int e = 0; e < pmesh.GetNE(); e++) { Dgf(e) = c.pf.DrivingFrom(psi_el(e)); }
         }
         for (int i = 0; i < D.Size(); i++) { lo(i) = std::max(D(i), seed(i)); }
         dx = lo;   // the projected start point d_old (SNESVI projects onto the bounds)
         bool assembled = false;
         if (damage_abs_tol < 0.0)
         {
            damage->Assemble(Dgf);
            assembled = true;
            Vector diag;
            damage->A().GetDiag(diag);
            real_t loc = diag.Normlinf(), glob = 0.0;
            MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
            damage_abs_tol = fracture_nl_abs_tol * glob;
            if (Mpi::Root())
            {
               mfem::out << "  damage automatic scaling factor 1/max A_ii = " << 1.0 / glob
                         << ", absolute tolerance on |r| = " << damage_abs_tol << "\n";
            }
         }
         damage->Residual(dx, Dgf, rd);
         if (InactiveResidualNorm(comm, dx, rd, lo, hi) > damage_abs_tol)
         {
            if (!assembled) { damage->Assemble(Dgf); }
            qp.abs_tol = damage_abs_tol;
            BoundConstrainedQP::Stats st;
            if (!qp.Solve(damage->A(), damage->b(), lo, hi, dx, st))
            {
               MFEM_ABORT("step " << step_n << " (t = " << S.t << "): damage QP did not converge, |r_I| = "
                          << st.r << " (|r_I,0| = " << st.r0 << ")");
            }
            qp_its_step = std::max(st.iterations, 1);
            damage_solves++;
            damage_qp_its += st.iterations;
         }
         D = dx;
      }
      if (step_n >= 1)
      {
         real_t li = 0.0, lo = 0.0;
         for (int i = 0; i < n; i++)
         {
            const real_t du2 = 0.5 * (S.u(i) - u_nm1(i));
            li += info.f_ext(i) * du2;
            lo += (info.f_damp(i) + info.f_kv(i) + info.f_mass(i)) * du2;
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
      real_t dWe = 0.0, dWd = 0.0, dWk = 0.0, dWm = 0.0, dWc = 0.0;
      {
         // u^n = S.u - dt v^{n+1/2}
         Vector un(S.u);
         un.Add(-dt, S.v_half);
         for (int i = 0; i < n; i++)
         {
            const real_t d_u = un(i) - u_prev(i);
            // fext_prev and info.f_ext hold pulse + confinement; book them apart
            dWe += (0.5 * (fext_prev(i) + info.f_ext(i)) - f_conf(i)) * d_u;
            dWc += f_conf(i) * d_u;
            dWd += 0.5 * (fdamp_prev(i) + info.f_damp(i)) * d_u;
            dWk += 0.5 * (fkv_prev(i) + info.f_kv(i)) * d_u;
            dWm += 0.5 * (fmass_prev(i) + info.f_mass(i)) * d_u;
         }
         u_prev = un;
      }
      real_t loc[5] = {dWe, dWd, dWk, dWm, dWc}, glob[5];
      MPI_Allreduce(loc, glob, 5, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
      W_ext += glob[0]; W_damp += glob[1]; W_kv += glob[2]; W_mass += glob[3]; W_conf += glob[4];
      fext_prev = info.f_ext; fdamp_prev = info.f_damp; fkv_prev = info.f_kv;
      fmass_prev = info.f_mass;

      const real_t ke = op.KineticEnergy(info.v_n);
      const real_t vmax = op.MaxNodalSpeed(info.v_n);
      // info and the works refer to t^n, where the damage is d^n.
      const real_t Ec = crack_energy(d_n);
      // The seeded initial damage (lo = seed at the first solve) is not part of
      // the budget: MOOSE subtracts dissipated_energy_first_step; here Ec at t^1.
      if (step_n == 1) { Ec_first = Ec; }
      real_t err = ke + info.psi + (Ec - Ec_first) - W_ext - W_conf + W_damp + W_kv + W_mass;
      if (hydromech)
      {
         err += Ef_out + Wfd_out - Wfb_out;
         if (step_n <= 1) { err_ref = err; }
         err -= err_ref;
      }
      real_t dmax = 0.0, ext = 0.0, ext90 = 0.0;
      if (!frozen_damage)
      {
         real_t loc = d_n.Size() ? d_n.Max() : 0.0;
         MPI_Allreduce(&loc, &dmax, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
         ext = extent(d_n, c.d_threshold);
         ext90 = extent(d_n, 0.9);
      }
      if (Mpi::Root())
      {
         csv << step_n << ',' << t_n << ',' << dt << ',' << op.PeakOnPressureBoundary() * s << ','
             << vmax << ',' << ke << ',' << info.psi << ',' << W_ext << ',' << W_damp << ','
             << W_kv << ',' << err << ',' << H << ',' << Ws_in << ',' << Ws_out << ','
             << shadow_err_report << ',' << Ec << ',' << dmax << ',' << ext << ',' << ext90
             << ',' << qp_its_step << ',' << W_mass << ',' << W_conf;
         if (hydromech)
         {
            csv << ',' << Ef_out << ',' << Wfd_out << ',' << Wfb_out << ',' << p_rng[0] << ','
                << p_rng[1] << ',' << pinfo.iterations;
         }
         csv << '\n';
         if (step_n % 200 == 0)
         {
            mfem::out << "step " << step_n << " t = " << t_n << " dt = " << dt << " v_max = "
                      << vmax << " KE = " << ke << " EE = " << info.psi << " W_ext = " << W_ext
                      << " W_damp = " << W_damp << " err = " << err << " shadow err = "
                      << shadow_err_report;
            if (!frozen_damage)
            {
               mfem::out << " d_max = " << dmax << " extent = " << ext << " crack = " << Ec;
            }
            if (hydromech)
            {
               mfem::out << " p = [" << p_rng[0] << ", " << p_rng[1] << "] E_f = " << Ef_out;
            }
            mfem::out << "\n";
         }
      }
      if (pv && step_n % pv_interval == 0)
      {
         vel.SetFromTrueDofs(info.v_n);
         Vector un(S.u);
         un.Add(-dt, S.v_half);
         disp.SetFromTrueDofs(un);
         dout.SetFromTrueDofs(d_n);
         if (hydromech) { ppout.SetFromTrueDofs(p_n); }
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
      if (hydromech)
      {
         mfem::out << "  pore-pressure PCG iterations " << pp_its_total << "\n";
      }
      if (!frozen_damage)
      {
         mfem::out << "  damage QP solves " << damage_solves << " (" << damage_qp_its
                   << " iterations), skipped steps " << S.step - damage_solves << "\n";
      }
   }
   return 0;
}
