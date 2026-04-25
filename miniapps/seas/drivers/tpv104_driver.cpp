// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV104 benchmark driver — 3D dynamic rupture on a vertical
// strike-slip fault with the slip-law-with-SRW (FVW / SCEC FL=103)
// friction law.
//
// Fully parallel MPI driver using ParMesh and WaveOperator<ParMesh>.
// Duplicates drivers/tpv102_driver.cpp's ADER time-stepping branch
// with plan §4.10 Step 9 substitutions:
//   - TPV104Params material / friction constants (not TPV102Params).
//   - InitializeFaultDOFs_TPV104 (per-QP a(x,z), ψ_ini via FVW inversion).
//   - Per-QP V_w side-channel via PopulateVwSideChannel_TPV104.
//   - ApplyNucleationIncremental_TPV104 (cumulative) per macro-step —
//     plan §3.9 directive, replaces TPV102's overwrite pattern.
//   - UpdateStateAnalyticSlipLawSRW (FVW analytic ψ) per macro-step —
//     plan §3.3 + §3.12 directives, replaces TPV102's forward-Euler
//     on AgingLawPsi.
//   - TPV104StationWriter / TPV104SurfaceStationWriter (SCEC trace
//     column order per plan §4.4).
//
// Usage:
//   ibrun ./seas_tpv104_driver \
//         --mesh tpv104/mesh/tpv104_200m.msh \
//         --tfinal 3.0 --ader-order 2 \
//         --friction-solver newton-stable \
//         --output-dir tpv104/results

#include "mfem.hpp"
#include "../dynamic/wave_state.hpp"
#include "../dynamic/wave_operator.hpp"
#include "../dynamic/fault_face_flux.hpp"
#include "../dynamic/friction_solver.hpp"
#include "../dynamic/tpv104_setup.hpp"
#include "../dynamic/tpv104_nucleation.hpp"
#include "../dynamic/tpv104_substep_iterator.hpp"
#include "../config/tpv104_params.hpp"
#include "../domain/boundary_config.hpp"
#include "../friction/slip_law_srw_psi.hpp"
#include "../dynamic/seas_diag_rank.hpp"
#include "../io/paraview_output.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <sys/stat.h>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

#ifdef SEAS_DIAG_FAULT_FLUX
// Global rank cache declared in dynamic/seas_diag_rank.hpp — one
// definition for the whole link.  Same pattern as tpv102_driver.cpp.
namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT
#endif

// --------------------------------------------------------------------------
// Small CLI parsing helpers (same convention as tpv102_driver.cpp so
// launch scripts can share argument shapes).
// --------------------------------------------------------------------------
static std::string GetStringArg(int argc, char *argv[], const char *flag,
                                const std::string &default_val)
{
   for (int i = 1; i < argc - 1; ++i)
   {
      if (std::string(argv[i]) == flag) { return argv[i + 1]; }
   }
   return default_val;
}

static real_t GetRealArg(int argc, char *argv[], const char *flag,
                         real_t default_val)
{
   const std::string val = GetStringArg(argc, argv, flag, "");
   if (val.empty()) { return default_val; }
   return std::stod(val);
}

static int GetIntArg(int argc, char *argv[], const char *flag,
                     int default_val)
{
   const std::string val = GetStringArg(argc, argv, flag, "");
   if (val.empty()) { return default_val; }
   return std::stoi(val);
}

static bool HasFlag(int argc, char *argv[], const char *flag)
{
   for (int i = 1; i < argc; ++i)
   {
      if (std::string(argv[i]) == flag) { return true; }
   }
   return false;
}

// Map the --friction-solver CLI name to FrictionSolver::Method (plan §4.10.X).
//
// R7-001/R7-005 note: on the current driver path this value is kept only
// for future iterator wiring.  The production time loop runs Brent via
// wave.AdvanceADER -> FaultFaceFlux::EvaluateADERTotal (see the honest
// banner below); the returned Method is not routed through that call.
// Nevertheless, MapSolver is kept strict so that (a) `newton` (bare) is
// the canonical shorthand for the stable-asinh variant and (b) unknown
// strings abort loudly rather than silently defaulting — both become
// load-bearing the moment the iterator is wired into the time loop.
static FrictionSolver::Method MapSolver(const std::string &s)
{
   if (s == "newton-stable" || s == "newton")
   {
      return FrictionSolver::Method::NewtonRaphsonStable;
   }
   if (s == "brent")        { return FrictionSolver::Method::Brent; }
   if (s == "newton-legacy"){ return FrictionSolver::Method::NewtonRaphson; }
   if (s == "hybrid")       { return FrictionSolver::Method::HybridNRBisection; }
   MFEM_ABORT("--friction-solver: unknown value '" << s
              << "'.  Accepted: newton-stable | newton | brent | "
              << "newton-legacy | hybrid.");
}

static std::string SolverBanner(const std::string &s)
{
   if (s == "newton-stable" || s == "newton")
   {
      return "Newton-Raphson (stable-asinh, plan §4.10 Step 5)";
   }
   if (s == "brent")        { return "Brent (log10-V, Tandem-verified; legacy dispatch)"; }
   if (s == "newton-legacy"){ return "Newton-Raphson (legacy, MFEM-native μ)"; }
   if (s == "hybrid")       { return "Hybrid NR+Bisection (legacy MFEM μ)"; }
   return std::string("unknown (") + s + ")";
}

// R8-001: route the dispatched-classification tags through an enum +
// `switch` rather than hardcoded constants.  Today — under R7-001
// option (b) — GetDispatched*(cli) unconditionally returns the
// Brent / OneShot / SlipSRW enum value regardless of CLI, because the
// driver's time loop ignores the CLI solver / iterator flags.  When
// R7-001 option (a) wires the iterator in, updating the Get*
// functions to branch on `cli` is the single required edit; TagOf() /
// BannerOf() stay correct because they're derived from the enum.
// Forgetting to extend those switches for a new enum value is a
// compiler warning (-Wswitch), preventing hardcoded-tag drift.
enum class DispatchedSolver { Brent, NewtonRaphsonStable, NewtonRaphsonLegacy, Hybrid };
enum class DispatchedIterator { OneShot, SubStep };
enum class DispatchedLaw { SlipSRW, Aging };

static DispatchedSolver GetDispatchedSolver(const std::string &/*friction_solver_cli*/)
{
   // R7-001 option (b): wave.AdvanceADER -> EvaluateADERTotal hard-codes
   // the default Method::Brent argument; the CLI value does not reach
   // the solver dispatch.  Replace with a CLI-aware switch when the
   // iterator is wired (option a).
   return DispatchedSolver::Brent;
}
static DispatchedIterator GetDispatchedIterator(const std::string &/*fault_iterator_cli*/)
{
   return DispatchedIterator::OneShot;
}
static DispatchedLaw GetDispatchedLaw(const std::string &/*fric_law_cli*/)
{
   return DispatchedLaw::SlipSRW;
}

static const char *TagOf(DispatchedSolver s)
{
   switch (s)
   {
      case DispatchedSolver::Brent:                return "brent";
      case DispatchedSolver::NewtonRaphsonStable:  return "newton-stable";
      case DispatchedSolver::NewtonRaphsonLegacy:  return "newton-legacy";
      case DispatchedSolver::Hybrid:               return "hybrid";
   }
   return "unknown";
}
static const char *TagOf(DispatchedIterator i)
{
   switch (i)
   {
      case DispatchedIterator::OneShot: return "oneshot";
      case DispatchedIterator::SubStep: return "substep";
   }
   return "unknown";
}
static const char *TagOf(DispatchedLaw l)
{
   switch (l)
   {
      case DispatchedLaw::SlipSRW: return "slip-srw";
      case DispatchedLaw::Aging:   return "aging";
   }
   return "unknown";
}

static std::string BannerOf(DispatchedSolver s)
{
   switch (s)
   {
      case DispatchedSolver::Brent:
         return "Brent (hard-coded via EvaluateADERTotal; "
                "--friction-solver flag IGNORED)";
      case DispatchedSolver::NewtonRaphsonStable:
         return "Newton-Raphson (stable-asinh, plan §4.10 Step 5)";
      case DispatchedSolver::NewtonRaphsonLegacy:
         return "Newton-Raphson (legacy, MFEM-native μ)";
      case DispatchedSolver::Hybrid:
         return "Hybrid NR+Bisection (legacy MFEM μ)";
   }
   return "unknown";
}
static std::string BannerOf(DispatchedIterator i)
{
   switch (i)
   {
      case DispatchedIterator::OneShot:
         return "one-shot (sub-step iterator NOT wired; "
                "--fault-iterator flag IGNORED — R7-001)";
      case DispatchedIterator::SubStep:
         return "sub-step (Tpv104SubStepIterator, plan §4.10 Step 7)";
   }
   return "unknown";
}
static std::string BannerOf(DispatchedLaw l)
{
   switch (l)
   {
      case DispatchedLaw::SlipSRW:
         return "slip-SRW (ψ-space, macro-step analytic "
                "UpdateStateAnalyticSlipLawSRW)";
      case DispatchedLaw::Aging:
         return "aging (ψ-space, forward-Euler)";
   }
   return "unknown";
}

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank, nprocs;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);
#else
   int rank = 0, nprocs = 1;
#endif

#ifdef SEAS_DIAG_FAULT_FLUX
   mfem::seas::g_seas_my_rank = rank;
#endif

   if (rank == 0)
   {
#ifdef SEAS_DIAG_FAULT_FLUX
      const char *diag_fault_flux = "SEAS_DIAG_FAULT_FLUX = ON";
#else
      const char *diag_fault_flux = "SEAS_DIAG_FAULT_FLUX = OFF";
#endif
#ifdef SEAS_DIAG_TPV104_STATE
      const char *diag_tpv104 = "SEAS_DIAG_TPV104_STATE = ON";
#else
      const char *diag_tpv104 = "SEAS_DIAG_TPV104_STATE = OFF";
#endif
      std::fprintf(stderr, "[BUILD] %s\n", diag_fault_flux);
      std::fprintf(stderr, "[BUILD] %s\n", diag_tpv104);
      std::ofstream binfo("build_info.txt");
      if (binfo.is_open())
      {
         binfo << "[BUILD] " << diag_fault_flux << "\n";
         binfo << "[BUILD] " << diag_tpv104 << "\n";
         binfo.close();
      }
   }

   // -----------------------------------------------------------------------
   // Parse command-line arguments
   // -----------------------------------------------------------------------
   std::string mesh_file = GetStringArg(argc, argv, "--mesh",
                                        "tpv104/mesh/tpv104_200m.msh");
   real_t mesh_scale     = GetRealArg(argc, argv, "--mesh-scale", 1.0);
   int order             = GetIntArg(argc, argv, "--order", 1);
   std::string bc_mode   = GetStringArg(argc, argv, "--bc-mode", "absorbing");
   real_t tfinal         = GetRealArg(argc, argv, "--tfinal", TPV104Params::t_final);
   std::string output_dir    = GetStringArg(argc, argv, "--output-dir", "tpv104/results");
   std::string output_prefix = GetStringArg(argc, argv, "--output-prefix", "tpv104");
   real_t cfl_factor     = GetRealArg(argc, argv, "--cfl", 0.5);
   int bc_free           = GetIntArg(argc, argv, "--bc-free", 1);
   int bc_fault          = GetIntArg(argc, argv, "--bc-fault", 3);
   int bc_absorb         = GetIntArg(argc, argv, "--bc-absorb", 5);
   int ader_order        = GetIntArg(argc, argv, "--ader-order", 2);
   real_t dt_override    = GetRealArg(argc, argv, "--dt", 0.0);
   real_t output_dt      = GetRealArg(argc, argv, "--output-dt", 0.01);
   bool disable_nucleation = HasFlag(argc, argv, "--disable-nucleation");
   bool debug_qnorm      = HasFlag(argc, argv, "--debug-qnorm");
   bool dry_run          = HasFlag(argc, argv, "--dry-run");
   bool verify_dispatch  = HasFlag(argc, argv, "--verify-dispatch");

   // ParaView output controls — mirror tpv102_driver.cpp + BP5 conventions:
   //   --paraview              : enable PVD/VTU output, interval matches --output-dt
   //   --paraview-every N      : write every N steps
   //   --paraview-dt X         : write every X seconds (overrides step interval)
   //   --paraview-bulk-dt X    : enable a SECOND collection in ParaView_bulk/
   //                             with velocity + sigma_yy/sigma_xy/sigma_xz at
   //                             coarser cadence (typical: 0.05 s)
   //   --pv-low-order          : linear tets only (~40x smaller volume output)
   //   --no-domain-pv          : suppress fault-schedule volume save (fault-
   //                             surface PVD/VTU still written; bulk collection
   //                             unaffected)
   bool use_paraview = HasFlag(argc, argv, "--paraview");
   bool pv_low_order = HasFlag(argc, argv, "--pv-low-order");
   bool pv_no_domain = HasFlag(argc, argv, "--no-domain-pv");
   int  paraview_step_interval = GetIntArg(argc, argv, "--paraview-every", 0);
   real_t paraview_dt_flag     = GetRealArg(argc, argv, "--paraview-dt", 0.0);
   real_t paraview_bulk_dt     = GetRealArg(argc, argv, "--paraview-bulk-dt", 0.0);
   if (paraview_step_interval > 0 || paraview_dt_flag > 0.0
       || paraview_bulk_dt > 0.0)
   {
      use_paraview = true;
   }
   // `--dry-run` is a shortcut for "no mesh, no time-stepping,
   // just print banner + verify wiring compiles/runs".  Used by
   // test_tpv104_smoke.cpp and the banner-check sbatch on Frontera.
   // `--verify-dispatch` (R7-004) additionally emits machine-readable
   // [dispatch] lines describing what the production time loop ACTUALLY
   // runs — the tri-consistency check between banner claims and runtime
   // dispatch.  Under R7-001 option (b) the dispatch is always
   // Brent/one-shot/slip-SRW, independent of --friction-solver /
   // --fault-iterator / --fric-law values.

   // TPV104-specific CLI (plan §4.10 Step 9).
   // R7-001/R7-003/R7-006: these flags are accepted so smoke tests and
   // sbatch scripts can exercise banner parity, but on the current
   // driver path (one-shot wave.AdvanceADER) they have NO effect on the
   // dispatched solver / iterator / friction law.  Runtime is:
   //   - friction solver: Brent (hard-coded via EvaluateADERTotal)
   //   - fault iterator : one-shot (Tpv104SubStepIterator not wired
   //                      — requires exposing per-sub-step I± from
   //                      wave_operator.inl; that file is on the
   //                      extreme-care no-touch list)
   //   - friction law   : slip-SRW via the free function
   //                      UpdateStateAnalyticSlipLawSRW (per-macro-step)
   // The banner below mirrors this disclosure exactly.
   std::string friction_solver =
      GetStringArg(argc, argv, "--friction-solver", "newton-stable");
   std::string fric_law =
      GetStringArg(argc, argv, "--fric-law", "slip-srw");
   std::string fault_iterator =
      GetStringArg(argc, argv, "--fault-iterator", "substep");

   // ader-order is accepted verbatim; wave.AdvanceADER clamps/validates
   // internally.  Allowing it through avoids false warnings when the
   // smoke test exercises --ader-order 5 for banner-text verification.
   if (ader_order < 1) { ader_order = 2; }

   if (rank == 0)
   {
      std::cout << "========================================\n";
      std::cout << "SCEC TPV104 Dynamic Rupture Simulation\n";
      std::cout << "========================================\n";
      std::cout << "Ranks: " << nprocs << "\n";
      std::cout << "Mesh: " << mesh_file << "\n";
      std::cout << "Scale: " << mesh_scale << "\n";
      std::cout << "Order: " << order << "\n";
      std::cout << "BC: " << bc_mode << "\n";
      std::cout << "Tfinal: " << tfinal << " s\n";
      std::cout << "Output: " << output_dir << "/" << output_prefix << "\n";
      std::cout << "CFL: " << cfl_factor << "\n";
      std::cout << "BC attrs: free=" << bc_free
                << ", fault=" << bc_fault
                << ", absorb=" << bc_absorb << "\n";
      // R7-001 option (b): the banner describes what the driver
      // ACTUALLY runs, not what the CLI requested.  The CLI values
      // (friction_solver, fault_iterator) are echoed on a separate
      // "CLI parsed (banner-only)" line so log parsers can still see
      // them, but the four load-bearing lines are the dispatch truth
      // derived from GetDispatched* / BannerOf (R8-001).
      const DispatchedSolver   actual_solver = GetDispatchedSolver(friction_solver);
      const DispatchedIterator actual_iter   = GetDispatchedIterator(fault_iterator);
      const DispatchedLaw      actual_law    = GetDispatchedLaw(fric_law);

      std::cout << "Time integrator: ADER-O" << ader_order
                << " (one-shot via wave.AdvanceADER)\n";
      std::cout << "Fault iterator: " << BannerOf(actual_iter) << "\n";
      std::cout << "Friction solver: " << BannerOf(actual_solver) << "\n";
      std::cout << "Friction law: " << BannerOf(actual_law) << "\n";
      std::cout << "Nucleation: "
                << (disable_nucleation ? "DISABLED"
                                       : "enabled (TPV104, macro-step "
                                         "incremental telescoping)")
                << "\n";
      std::cout << "CLI parsed (banner-only, not dispatched): "
                << "friction_solver=" << friction_solver
                << ", fault_iterator=" << fault_iterator
                << ", fric_law=" << fric_law << "\n";
      std::cout << "========================================\n\n";
      mkdir(output_dir.c_str(), 0755);
   }

   // R8-004: per-rank machine-readable tri-consistency lines.  Emit
   // OUTSIDE the rank==0 guard so multi-rank Frontera logs can self-
   // verify that every rank dispatches the same solver.  Each line is
   // prefixed with `rank=<N>` so interleaved output is still parseable.
   if (verify_dispatch)
   {
      const DispatchedSolver   actual_solver = GetDispatchedSolver(friction_solver);
      const DispatchedIterator actual_iter   = GetDispatchedIterator(fault_iterator);
      const DispatchedLaw      actual_law    = GetDispatchedLaw(fric_law);
      std::cout << "[dispatch] rank=" << rank
                << " friction_solver_actual=" << TagOf(actual_solver) << "\n";
      std::cout << "[dispatch] rank=" << rank
                << " fault_iterator_actual=" << TagOf(actual_iter) << "\n";
      std::cout << "[dispatch] rank=" << rank
                << " friction_law_actual=" << TagOf(actual_law) << "\n";
   }

#ifdef MFEM_USE_MPI
   MPI_Barrier(comm);
#endif

   // R7-005: validate --friction-solver eagerly so typos abort before
   // any simulation work (including --dry-run).
   //
   // R8-002: preserve the Method value in a named local (`method`)
   // rather than discarding MapSolver's return.  Under R7-001 option
   // (b) the value is not routed through the solver dispatch, so it
   // is cast to void here.  When R7-001 option (a) lands, the
   // required edit is a single site: remove `(void)method;` and pass
   // `method` into `Tpv104SubStepIterator::Advance`.  Keeping the
   // named local ensures grep / IDE reference finds the linkage point.
   const FrictionSolver::Method method = MapSolver(friction_solver);
   (void)method;  // R7-001 option (b): not routed through the time loop.

   // --dry-run: no mesh, no simulation.  Print a canonical end-of-run
   // line that test_tpv104_smoke.cpp scrapes ("[dry-run] OK.").  Used
   // as the cheapest possible Frontera startup verification.
   if (dry_run)
   {
      if (rank == 0)
      {
         std::cout << "[tpv104_driver] --dry-run: banner printed, "
                   << "no mesh / no simulation.\n";
         std::cout << "[dry-run] OK.\n";
      }
      Tpv104SubStepIterator::CloseAllProbeFiles();
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 0;
   }

   // -----------------------------------------------------------------------
   // 1. Load serial mesh, partition to ParMesh
   // -----------------------------------------------------------------------
   Mesh serial_mesh(mesh_file.c_str(), 1, 1);
   MFEM_VERIFY(serial_mesh.Dimension() == 3, "TPV104 requires 3D mesh");
   if (mesh_scale != 1.0)
   {
      serial_mesh.SetCurvature(1, false, 3, Ordering::byVDIM);
      Vector &nodes = *serial_mesh.GetNodes();
      nodes *= mesh_scale;
   }

#ifdef MFEM_USE_MPI
   ParMesh pmesh(comm, serial_mesh);
   using MeshT = ParMesh;
#else
   Mesh &pmesh = serial_mesh;
   using MeshT = Mesh;
#endif

   int ne_local = pmesh.GetNE();
   int ne_global = ne_local;
#ifdef MFEM_USE_MPI
   MPI_Allreduce(&ne_local, &ne_global, 1, MPI_INT, MPI_SUM, comm);
#endif
   if (rank == 0)
   {
      std::cout << "Mesh: " << ne_global << " elements, "
                << nprocs << " ranks\n";
   }

   // -----------------------------------------------------------------------
   // 2. Boundary conditions
   // -----------------------------------------------------------------------
   BoundaryConfig bc;
   bc.natural_attrs   = {bc_free};
   bc.fault_attr      = bc_fault;
   bc.absorbing_attrs = {bc_absorb};

   {
      int n_free = 0, n_fault = 0, n_absorb = 0;
      for (int b = 0; b < pmesh.GetNBE(); b++)
      {
         int attr = pmesh.GetBdrAttribute(b);
         if (bc.natural_attrs.count(attr))   { ++n_free; }
         if (attr == bc.fault_attr)          { ++n_fault; }
         if (bc.absorbing_attrs.count(attr)) { ++n_absorb; }
      }
      int n_free_g = n_free, n_fault_g = n_fault, n_absorb_g = n_absorb;
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&n_free,   &n_free_g,   1, MPI_INT, MPI_SUM, comm);
      MPI_Allreduce(&n_fault,  &n_fault_g,  1, MPI_INT, MPI_SUM, comm);
      MPI_Allreduce(&n_absorb, &n_absorb_g, 1, MPI_INT, MPI_SUM, comm);
#endif
      if (rank == 0)
      {
         std::cout << "BC faces — free: " << n_free_g
                   << ", fault: " << n_fault_g
                   << ", absorb: " << n_absorb_g << "\n";
      }
      MFEM_VERIFY(n_fault_g > 0,
                  "No fault faces with attr=" << bc.fault_attr
                  << " found.  Check --bc-fault or mesh Physical Surface tags.");
   }

   // -----------------------------------------------------------------------
   // 3. WaveOperator
   // -----------------------------------------------------------------------
   WaveOperator<MeshT> wave(pmesh, order,
                            TPV104Params::lambda, TPV104Params::mu,
                            TPV104Params::rho, bc);

   int ndof_total  = wave.GetScalarNDof();
   int ndof_global = ndof_total;
#ifdef MFEM_USE_MPI
   MPI_Allreduce(&ndof_total, &ndof_global, 1, MPI_INT, MPI_SUM, comm);
#endif
   if (rank == 0)
   {
      std::cout << "DOFs per component (global): " << ndof_global << "\n";
      std::cout << "Total DOFs (global): " << NUM_STATE * ndof_global << "\n";
   }

   // -----------------------------------------------------------------------
   // 4. Time step
   // -----------------------------------------------------------------------
   real_t cfl    = cfl_factor / (3.0 * (2.0 * order + 1.0));
   real_t dt_cfl = wave.ComputeMaxDt(cfl);
   real_t dt     = (dt_override > 0.0) ? dt_override : dt_cfl;
   int nsteps    = (tfinal > 0.0) ? static_cast<int>(std::ceil(tfinal / dt)) : 0;
   if (rank == 0)
   {
      std::cout << "CFL: " << cfl << ", dt_cfl = " << dt_cfl << " s\n";
      std::cout << "dt: " << dt << " s, steps: " << nsteps << "\n\n";
   }

   // -----------------------------------------------------------------------
   // 5. Fault DOF data (local partition) — TPV104-specific init
   // -----------------------------------------------------------------------
   L2_FECollection fec(order, 3, BasisType::GaussLobatto);
#ifdef MFEM_USE_MPI
   ParFiniteElementSpace fes(&pmesh, &fec);
#else
   FiniteElementSpace fes(&pmesh, &fec);
#endif

   const Array<int> &fault_int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &fault_shr_faces = wave.GetFaultSharedFaces();

   int nqp_per_face = 0;
   if (fault_int_faces.Size() > 0)
   {
      FaceElementTransformations *ftr0 =
         pmesh.GetInteriorFaceTransformations(fault_int_faces[0]);
      MFEM_VERIFY(ftr0, "fault interior face has null transformation");
      nqp_per_face = IntRules.Get(ftr0->GetGeometryType(), 2 * order).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   else if (fault_shr_faces.Size() > 0)
   {
      FaceElementTransformations *ftr0 =
         pmesh.GetSharedFaceTransformations(fault_shr_faces[0]);
      MFEM_VERIFY(ftr0, "fault shared face has null transformation");
      nqp_per_face = IntRules.Get(ftr0->GetGeometryType(), 2 * order).GetNPoints();
   }
   {
      int local_nqp = nqp_per_face;
      MPI_Allreduce(&local_nqp, &nqp_per_face, 1, MPI_INT, MPI_MAX, comm);
   }
#endif

   const int num_fault_local  = fault_int_faces.Size() * nqp_per_face;
   const int num_shared_fault = fault_shr_faces.Size() * nqp_per_face;
   const int num_fault_total  = num_fault_local + num_shared_fault;

   std::vector<Vector> fault_coords;
   fault_coords.reserve(num_fault_total);

   auto push_qps = [&](FaceElementTransformations *ftr)
   {
      const IntegrationRule &ir = IntRules.Get(
         ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); ++q)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   };
   for (int i = 0; i < fault_int_faces.Size(); ++i)
   {
      push_qps(pmesh.GetInteriorFaceTransformations(fault_int_faces[i]));
   }
#ifdef MFEM_USE_MPI
   for (int i = 0; i < fault_shr_faces.Size(); ++i)
   {
      push_qps(pmesh.GetSharedFaceTransformations(fault_shr_faces[i]));
   }
#endif

   int num_fault_global = num_fault_total;
#ifdef MFEM_USE_MPI
   MPI_Allreduce(&num_fault_total, &num_fault_global, 1, MPI_INT, MPI_SUM, comm);
#endif
   if (rank == 0)
   {
      std::cout << "Fault QPs (global): " << num_fault_global
                << " (local: " << num_fault_local
                << ", shared: " << num_shared_fault << ")\n";
   }

#ifdef SEAS_DIAG_TPV104_FAULT_BASIS
   // D1 instrumentation companion: dump fault QP coordinates indexed
   // by dof_idx (= the same index keyed by [diag-flip] in
   // wave_operator.inl).  Join via:
   //   awk '/diag-coords/{c[$2]=$0} /diag-flip/{f[$2]=$0} END{for(k in f)print f[k]" "c[k]}'
   for (int i = 0; i < num_fault_total; ++i)
   {
      const Vector &xqp = fault_coords[i];
      std::fprintf(stderr,
         "[diag-coords] dof=%d (x,y,z)=(%+9.1f,%+9.1f,%+9.1f) rank=%d\n",
         i, xqp(0), xqp(1), xqp(2), rank);
   }
#endif

   std::vector<DOFData> dof_data;
   std::vector<real_t>  V_w;
   if (num_fault_total > 0)
   {
      // TPV104 fluctuation-Q init: Q = 0, pre-stress lives in DOFData.
      InitializeFaultDOFs_TPV104(dof_data, num_fault_total, fault_coords);
      PopulateVwSideChannel_TPV104(V_w, fault_coords);
   }

   // Resolve the rank that owns the hypocenter QP (closest local fault QP
   // to (hypo_along_strike, -hypo_down_dip) in x/z), then tag that DOF
   // with diag_print = true so the C-1 EVAL / C-1n NORMAL probes in
   // dynamic/fault_face_flux.cpp emit one line per Evaluate call instead
   // of flooding stderr from every QP.  Mirrors the TPV102 pattern at
   // tpv102_driver.cpp:540-588.
   int hypo_rank = 0;
   int hypo_dof_local = -1;
   {
      real_t local_min_dist2 = std::numeric_limits<real_t>::max();
      for (int i = 0; i < num_fault_total; i++)
      {
         real_t dx = fault_coords[i](0) - TPV104Params::hypo_along_strike;
         real_t dz = std::abs(fault_coords[i](2)) - TPV104Params::hypo_down_dip;
         real_t d2 = dx*dx + dz*dz;
         if (d2 < local_min_dist2)
         {
            local_min_dist2 = d2;
            hypo_dof_local = i;
         }
      }
#ifdef MFEM_USE_MPI
      struct MinDist { double d; int r; };
      static_assert(offsetof(MinDist, d) == 0,
                    "MinDist.d must be at offset 0 for MPI_DOUBLE_INT");
      static_assert(offsetof(MinDist, r) == sizeof(double),
                    "MinDist.r must follow d with no padding for MPI_DOUBLE_INT");
      MinDist in{static_cast<double>(local_min_dist2), rank};
      MinDist out{};
      MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
      hypo_rank = out.r;
#endif
   }

#ifdef SEAS_DIAG_FAULT_FLUX
   // Only the rank that won the MINLOC sets diag_print = true on its
   // closest hypocenter DOF; at most one DOF is flagged globally.
   // Probe block: dynamic/fault_face_flux.cpp:259-303 (C-1 EVAL +
   // C-1n NORMAL).  Single-rank stderr printf, no MPI calls — no
   // deadlock at any rank count.
   //
   // C-2 BULK-PROBE EXTENSION: also identify the fault face containing
   // the diag DOF and the two adjacent tets, plus the 6 non-fault
   // interior faces of those tets.  Pushed to WaveOperator via setters
   // so the C-2A (Mult per-call) and C-2B (per non-fault face) probes
   // can fire without MPI on this single rank.  Other ranks have the
   // setters left at default -1 → silent.
   int diag_face_idx_for_dof = -1;
   int diag_elem_plus  = -1, diag_elem_minus = -1;
   int diag_face_dof_plus = -1, diag_face_dof_minus = -1;
   std::vector<int> diag_nonfault_faces;
   if (rank == hypo_rank && hypo_dof_local >= 0 && num_fault_total > 0)
   {
      dof_data[hypo_dof_local].diag_print = true;
      const Vector &hpos = fault_coords[hypo_dof_local];
      std::fprintf(stderr,
         "[diag] rank %d tagging hypo DOF %d at (%.1f, %.1f, %.1f)\n",
         rank, hypo_dof_local, hpos(0), hpos(1), hpos(2));

      // Map hypo_dof_local back to (interior fault face index, q-of-face).
      // Layout: dof_data[i] lives on interior face i / nqp_per_face,
      // QP index = i % nqp_per_face.  Shared-fault DOFs follow at
      // [num_fault_local, num_fault_total) — currently not used as
      // diag target (the hypo nucleus lives on interior faces).
      if (hypo_dof_local < num_fault_local)
      {
         const int fi = hypo_dof_local / nqp_per_face;     // index into fault_int_faces
         const int qi = hypo_dof_local % nqp_per_face;
         if (fi >= 0 && fi < fault_int_faces.Size())
         {
            diag_face_idx_for_dof = fault_int_faces[fi];
            FaceElementTransformations *ftr =
               pmesh.GetInteriorFaceTransformations(diag_face_idx_for_dof);
            if (ftr)
            {
               const int e1 = ftr->Elem1No;
               const int e2 = ftr->Elem2No;
               // Determine canonical + / - via FaultBasis sign_flipped at
               // this QP (matches wave_operator.inl:1226 logic:
               // elem1_on_plus = !sign_flipped).
               bool elem1_on_plus = true;
               const FaultBasis *fb = wave.GetFaultBasis();
               if (fb)
               {
                  const int fb_idx =
                     wave.LookupInteriorFaultBasisIndex(diag_face_idx_for_dof);
                  if (fb_idx >= 0 && fb_idx < fb->NumFaces())
                  {
                     const FaultBasisData &bd = fb->GetBasis(fb_idx);
                     if (qi < static_cast<int>(bd.qp_data.size()))
                     {
                        elem1_on_plus = !bd.qp_data[qi].sign_flipped;
                     }
                  }
               }
               diag_elem_plus  = elem1_on_plus ? e1 : e2;
               diag_elem_minus = elem1_on_plus ? e2 : e1;

               // Find the local DOF index on each tet that is closest
               // to the hypocenter QP physical position.  ndof_per_el =
               // (order+1)(order+2)(order+3)/6 for tet; for order=1
               // there are 4 DOFs per element, one at each vertex.
               auto closest_dof_idx = [&](int e) -> int
               {
                  const FiniteElement *fe = wave.GetFESpace().GetFE(e);
                  ElementTransformation *Tr =
                     wave.GetFESpace().GetElementTransformation(e);
                  const IntegrationRule &nodes = fe->GetNodes();
                  int best = -1;
                  real_t best_d2 = std::numeric_limits<real_t>::max();
                  for (int k = 0; k < nodes.GetNPoints(); k++)
                  {
                     Vector phys(3);
                     Tr->Transform(nodes.IntPoint(k), phys);
                     const real_t dx = phys(0) - hpos(0);
                     const real_t dy = phys(1) - hpos(1);
                     const real_t dz = phys(2) - hpos(2);
                     const real_t d2 = dx*dx + dy*dy + dz*dz;
                     if (d2 < best_d2) { best_d2 = d2; best = k; }
                  }
                  return best;
               };
               diag_face_dof_plus  = closest_dof_idx(diag_elem_plus);
               diag_face_dof_minus = closest_dof_idx(diag_elem_minus);

               // Collect the non-fault interior faces of the two diag
               // tets.  For tets, GetElementFaces returns 4 faces; we
               // skip the fault face itself and any boundary faces
               // (which are still printed elsewhere).
               Array<int> faces_p, ori_p, faces_m, ori_m;
               pmesh.GetElementFaces(diag_elem_plus,  faces_p, ori_p);
               pmesh.GetElementFaces(diag_elem_minus, faces_m, ori_m);
               auto add_nonfault = [&](const Array<int> &fs)
               {
                  for (int k = 0; k < fs.Size(); k++)
                  {
                     if (fs[k] == diag_face_idx_for_dof) { continue; }
                     diag_nonfault_faces.push_back(fs[k]);
                  }
               };
               add_nonfault(faces_p);
               add_nonfault(faces_m);

               std::fprintf(stderr,
                  "[diag-c2] rank %d face=%d e+=%d e-=%d "
                  "face_dof+=%d face_dof-=%d nonfault_faces=[",
                  rank, diag_face_idx_for_dof,
                  diag_elem_plus, diag_elem_minus,
                  diag_face_dof_plus, diag_face_dof_minus);
               for (size_t k = 0; k < diag_nonfault_faces.size(); k++)
               {
                  std::fprintf(stderr, "%s%d",
                               k == 0 ? "" : ",", diag_nonfault_faces[k]);
               }
               std::fprintf(stderr, "]\n");
            }
         }
      }
   }
   wave.SetDiagBulkElems(diag_elem_plus, diag_elem_minus);
   wave.SetDiagBulkFaceDofs(diag_face_dof_plus, diag_face_dof_minus);
   wave.SetDiagNonFaultFaces(diag_nonfault_faces);
#endif

   FaultFaceFlux fault_flux(TPV104Params::rho, TPV104Params::cp,
                            TPV104Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   // Zero Q_bg — fluctuation-Q dispatch.
   {
      real_t Q_bg[NUM_STATE] = {0};
      wave.SetAbsorbingBackground(Q_bg);
   }

   // R7-002: no SlipLawSRWPsi instance is constructed here.  The ψ
   // update is invoked via the free function
   // UpdateStateAnalyticSlipLawSRW(...) in the time loop below with
   // per-QP (V_w, a), bypassing the base-virtual dispatch entirely.
   // SlipLawSRWPsi::SetProductionMode()'s R-001 safety net is a guard
   // against silent base-virtual fallthrough, which cannot occur on
   // this path.  When the sub-step iterator is eventually wired into
   // the time loop (R7-001 option a), a SlipLawSRWPsi instance with
   // SetProductionMode() should be constructed and bound here.
   // -----------------------------------------------------------------------
   // 6. Initialize Q = 0 (fluctuation-Q).
   // -----------------------------------------------------------------------
   Vector Q(NUM_STATE * ndof_total);
   InitializeState_TPV104(Q, ndof_total);

   MFEM_VERIFY(wave.GetAbsorbingBackground() != nullptr,
               "tpv104_driver: SetAbsorbingBackground(Q_bg=0) not called.");

   if (rank == 0)
   {
      std::cout << "Background: τ_strike = "
                << TPV104Params::tau_ini / 1e6 << " MPa, "
                << "σ_n = " << TPV104Params::sigma_n / 1e6 << " MPa, "
                << "ψ_ini(a_in) ≈ 0.5636, V_ini = "
                << TPV104Params::V_ini << " m/s\n\n";
   }

   // -----------------------------------------------------------------------
   // 7. Station output
   // -----------------------------------------------------------------------
   auto stations = DefaultStations_TPV104();
   TPV104StationWriter station_writer;
#ifdef MFEM_USE_MPI
   station_writer.Open(output_dir, output_prefix, stations,
                       fault_coords, num_fault_local, comm);
#else
   station_writer.Open(output_dir, output_prefix, stations,
                       fault_coords, num_fault_local);
#endif
   station_writer.WriteStep(0.0, dof_data);

   auto surface_stations = DefaultSurfaceStations_TPV104();
   TPV104SurfaceStationWriter surface_writer;
   surface_writer.Open(output_dir, output_prefix, surface_stations,
                       pmesh, fes);
   surface_writer.WriteStep(0.0, Q);

   // -----------------------------------------------------------------------
   // 7b. ParaView output (mirrors tpv102_driver.cpp / BP5 seas::ParaViewOutput
   //     pattern).  Two collections:
   //       - pv_out (output_dir/ParaView): velocity + mpi_rank volume +
   //         fault-surface PVD/VTU (slip, slip_rate, traction dip+strike,
   //         psi, sigma_n, plus static a, Dc, x2, x3) at the fault schedule.
   //       - pv_bulk_out (output_dir/ParaView_bulk): velocity + sigma_yy +
   //         sigma_xy + sigma_xz + mpi_rank at --paraview-bulk-dt cadence.
   //
   // R-801 / BP5 component convention enforced project-wide: comp 0 = dip,
   // comp 1 = strike.  TPV104 is pure strike-slip so the strike channel
   // carries the rupture; dip stays near zero.
   //
   // Under R7-001 option (b) ADER one-shot, there is no "stage-4" snapshot
   // distinct from the time-averaged DOFData — write the same values into
   // both the averaged and the _k4 buffers so the ParaView Calculator
   // delta `<field>_k4 - <field>` reads zero everywhere (consistent with
   // what tpv102_driver.cpp does on the ADER path).
   // -----------------------------------------------------------------------
   using PvFES = typename seas::GFType<MeshT>::FESType;
   using PvGF  = typename seas::GFType<MeshT>::type;

   std::unique_ptr<seas::ParaViewOutput<MeshT>> pv_out;
   std::unique_ptr<L2_FECollection> pv_vel_fec, pv_rank_fec;
   std::unique_ptr<PvFES> pv_vel_fes, pv_rank_fes;
   std::unique_ptr<PvGF>  pv_vel_gf,  pv_rank_gf;

   std::unique_ptr<seas::ParaViewOutput<MeshT>> pv_bulk_out;
   std::unique_ptr<L2_FECollection> pv_bulk_sigma_fec;
   std::unique_ptr<PvFES> pv_bulk_sigma_fes;
   std::unique_ptr<PvGF>  pv_bulk_syy_gf, pv_bulk_sxy_gf, pv_bulk_sxz_gf;

   Vector pv_local_slip, pv_local_slip_rate, pv_local_traction;
   Vector pv_local_state, pv_local_normal_stress;
   Vector pv_local_a, pv_local_Dc, pv_local_x2, pv_local_x3;
   Vector pv_local_slip_rate_k4, pv_local_traction_k4, pv_local_normal_stress_k4;

   if (use_paraview)
   {
      if (rank == 0) { mkdir((output_dir + "/ParaView").c_str(), 0755); }
#ifdef MFEM_USE_MPI
      MPI_Barrier(comm);
#endif
      pv_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
         output_dir + "/ParaView", pmesh, order);

      if (pv_low_order)
      {
         pv_out->SetHighOrderOutput(false);
         pv_out->SetLevelsOfDetail(1);
      }

      // Velocity: 3-component vector L2, byNODES so memcpy from Q's
      // [VX..VZ] block is a single contiguous copy.
      pv_vel_fec = std::make_unique<L2_FECollection>(order, 3, BasisType::GaussLobatto);
      pv_vel_fes = std::make_unique<PvFES>(&pmesh, pv_vel_fec.get(),
                                           3, Ordering::byNODES);
      pv_vel_gf  = std::make_unique<PvGF>(pv_vel_fes.get());
      *pv_vel_gf = 0.0;
      pv_out->RegisterDomainField("velocity", pv_vel_gf.get());

      // MPI rank: L2 p=0 (one value per element).
      pv_rank_fec = std::make_unique<L2_FECollection>(0, 3);
      pv_rank_fes = std::make_unique<PvFES>(&pmesh, pv_rank_fec.get());
      pv_rank_gf  = std::make_unique<PvGF>(pv_rank_fes.get());
      *pv_rank_gf = static_cast<real_t>(rank);
      pv_out->RegisterDomainField("mpi_rank", pv_rank_gf.get());

      // Fault L2-p0 fields keyed off the wave operator's canonical face
      // lists (interior + shared) so ParaView indexing matches DOFData.
      pv_out->InitFaultOutputBP5(fault_int_faces, fault_shr_faces,
                                 nqp_per_face);

      pv_local_slip.SetSize(2 * num_fault_total);
      pv_local_slip_rate.SetSize(2 * num_fault_total);
      pv_local_traction.SetSize(2 * num_fault_total);
      pv_local_state.SetSize(num_fault_total);
      pv_local_normal_stress.SetSize(num_fault_total);
      pv_local_slip_rate_k4.SetSize(2 * num_fault_total);
      pv_local_traction_k4.SetSize(2 * num_fault_total);
      pv_local_normal_stress_k4.SetSize(num_fault_total);

      pv_local_a.SetSize(num_fault_total);
      pv_local_Dc.SetSize(num_fault_total);
      pv_local_x2.SetSize(num_fault_total);
      pv_local_x3.SetSize(num_fault_total);
      for (int i = 0; i < num_fault_total; i++)
      {
         pv_local_a(i)  = dof_data[i].a;
         pv_local_Dc(i) = dof_data[i].Dc;
         pv_local_x2(i) = fault_coords[i](0);
         pv_local_x3(i) = fault_coords[i](2);
      }
      pv_out->SetFaultParamsBP5(pv_local_a, pv_local_Dc,
                                pv_local_x2, pv_local_x3);

      // Schedule (CLI > step-interval > output_dt step interval):
      const int output_interval_for_pv =
         std::max(1, static_cast<int>(output_dt / dt));
      if (paraview_step_interval > 0)
      {
         pv_out->output_every_n_steps = paraview_step_interval;
      }
      else if (paraview_dt_flag > 0.0)
      {
         pv_out->fixed_dt = paraview_dt_flag;
      }
      else
      {
         pv_out->output_every_n_steps = output_interval_for_pv;
      }

      if (rank == 0)
      {
         std::cout << "ParaView output: ON (prefix="
                   << output_dir << "/ParaView)\n";
         if (pv_no_domain)
         {
            std::cout << "  Mode: fault-surface PVD only (--no-domain-pv)\n";
         }
         if (paraview_step_interval > 0)
         {
            std::cout << "  Interval: every " << paraview_step_interval
                      << " steps (--paraview-every)\n";
         }
         else if (paraview_dt_flag > 0.0)
         {
            std::cout << "  Interval: every " << paraview_dt_flag
                      << " s (--paraview-dt)\n";
         }
         else
         {
            std::cout << "  Interval: every " << output_interval_for_pv
                      << " steps (matches --output-dt=" << output_dt
                      << " s)\n";
         }
      }

      if (paraview_bulk_dt > 0.0)
      {
         pv_bulk_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
            output_dir + "/ParaView_bulk", pmesh, order);
         if (pv_low_order)
         {
            pv_bulk_out->SetHighOrderOutput(false);
            pv_bulk_out->SetLevelsOfDetail(1);
         }
         pv_bulk_out->RegisterDomainField("velocity", pv_vel_gf.get());
         pv_bulk_out->RegisterDomainField("mpi_rank", pv_rank_gf.get());

         pv_bulk_sigma_fec = std::make_unique<L2_FECollection>(
            order, 3, BasisType::GaussLobatto);
         pv_bulk_sigma_fes = std::make_unique<PvFES>(&pmesh,
                                                     pv_bulk_sigma_fec.get());
         pv_bulk_syy_gf = std::make_unique<PvGF>(pv_bulk_sigma_fes.get());
         pv_bulk_sxy_gf = std::make_unique<PvGF>(pv_bulk_sigma_fes.get());
         pv_bulk_sxz_gf = std::make_unique<PvGF>(pv_bulk_sigma_fes.get());
         *pv_bulk_syy_gf = 0.0;
         *pv_bulk_sxy_gf = 0.0;
         *pv_bulk_sxz_gf = 0.0;
         pv_bulk_out->RegisterDomainField("sigma_yy", pv_bulk_syy_gf.get());
         pv_bulk_out->RegisterDomainField("sigma_xy", pv_bulk_sxy_gf.get());
         pv_bulk_out->RegisterDomainField("sigma_xz", pv_bulk_sxz_gf.get());

         pv_bulk_out->fixed_dt = paraview_bulk_dt;

         if (rank == 0)
         {
            std::cout << "  Bulk collection: ON (prefix="
                      << output_dir << "/ParaView_bulk, every "
                      << paraview_bulk_dt
                      << " s; fields: velocity, sigma_yy, sigma_xy, "
                      << "sigma_xz, mpi_rank)\n";
         }
      }
   }

   // paraview_write: MPI-collective snapshot writer.  V_max must already be
   // globally reduced; step_num/time identify the frame.
   auto paraview_write = [&](int step_num, real_t time, real_t V_max)
   {
      if (!pv_out) { return; }

      const bool fault_wants =
         pv_out->PeekShouldWrite(step_num, time, V_max);
      const bool bulk_wants = pv_bulk_out &&
         pv_bulk_out->PeekShouldWrite(step_num, time, V_max);
      if (!fault_wants && !bulk_wants) { return; }

      if (!pv_no_domain || bulk_wants)
      {
         std::memcpy(pv_vel_gf->GetData(),
                     Q.GetData() + VX * ndof_total,
                     3 * ndof_total * sizeof(real_t));
      }

      if (bulk_wants)
      {
         std::memcpy(pv_bulk_syy_gf->GetData(),
                     Q.GetData() + SYY * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_sxy_gf->GetData(),
                     Q.GetData() + SXY * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_sxz_gf->GetData(),
                     Q.GetData() + SXZ * ndof_total,
                     ndof_total * sizeof(real_t));
         pv_bulk_out->ForceSave(step_num, time);
      }

      if (!fault_wants) { return; }

      for (int i = 0; i < num_fault_total; i++)
      {
         const DOFData &d = dof_data[i];
         pv_local_slip(2*i + 0)      = d.slip1;       // dip
         pv_local_slip(2*i + 1)      = d.slip2;       // strike
         pv_local_slip_rate(2*i + 0) = d.V1;
         pv_local_slip_rate(2*i + 1) = d.V2;
         pv_local_traction(2*i + 0)  = d.tau1_corr;
         pv_local_traction(2*i + 1)  = d.tau2_corr;
         pv_local_state(i)           = d.psi;
         pv_local_normal_stress(i)   = d.sigma_n_corr;

         // ADER one-shot: no stage-4 distinct from averaged DOFData (R7-001
         // option b).  Mirror tpv102 ADER path: write the same values into
         // _k4 buffers so the VTU delta reads exactly zero.
         pv_local_slip_rate_k4(2*i + 0) = d.V1;
         pv_local_slip_rate_k4(2*i + 1) = d.V2;
         pv_local_traction_k4(2*i + 0)  = d.tau1_corr;
         pv_local_traction_k4(2*i + 1)  = d.tau2_corr;
         pv_local_normal_stress_k4(i)   = d.sigma_n_corr;
      }

      if (pv_no_domain)
      {
         pv_out->CommitSchedule(time);
      }
      else
      {
         pv_out->UpdateFaultFieldsBP5(pv_local_slip, pv_local_slip_rate,
                                      pv_local_traction, pv_local_state,
                                      pv_local_normal_stress);
         pv_out->ForceSave(step_num, time);
      }

      pv_out->WriteFaultSurfaceVTU(
         output_dir, step_num, time, rank, nprocs,
         pv_local_slip, pv_local_slip_rate, pv_local_traction,
         pv_local_state, pv_local_normal_stress,
         pv_local_a, pv_local_Dc, pv_local_x2, pv_local_x3,
         pv_local_slip_rate_k4, pv_local_traction_k4,
         pv_local_normal_stress_k4);
   };

   // Initial snapshot at t=0 (V_max = V_ini since fault is quasi-static).
   paraview_write(0, 0.0, TPV104Params::V_ini);

   // If tfinal == 0: init-only run, stations at t=0 already written.
   // Skip the time loop and go straight to summary.  This is the
   // Phase-2 P2_D gate — the init sbatch verifies ψ_ini at every
   // station from the single t=0 row.
   if (nsteps == 0)
   {
      if (rank == 0)
      {
         std::cout << "[tpv104_driver] tfinal = 0 — init-only run. "
                   << "Station t=0 row written, exiting.\n";
      }
      station_writer.Flush();
      station_writer.Close();
      surface_writer.Flush();
      surface_writer.Close();
      Tpv104SubStepIterator::CloseAllProbeFiles();
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 0;
   }

   // -----------------------------------------------------------------------
   // 8. ADER time stepping (plan §4.10 Step 9: ADER-only, no RK4 path)
   // -----------------------------------------------------------------------
   int output_interval = std::max(1, static_cast<int>(output_dt / dt));

   real_t t = 0.0;
   real_t V_max_global = 0.0;
   Vector Q_new(Q.Size());
   std::vector<real_t> psi_n(num_fault_total);

   if (rank == 0)
   {
      std::cout << "Starting ADER-O(" << ader_order << ") time loop...\n";
   }

   for (int step = 0; step < nsteps; ++step)
   {
      real_t dt_step = std::min(dt, tfinal - t);
      if (dt_step <= 0.0) { break; }

      // Save ψ at step start for the analytic update below.
      for (int i = 0; i < num_fault_total; ++i)
      {
         psi_n[i] = dof_data[i].psi;
      }

      // Nucleation (§3.9 directive): cumulative per-macro-step increment.
      // ApplyNucleationIncremental_TPV104 adds  Δτ · smoothStepIncrement
      // to `tau2_nuc` at every call, so over [0, T_nuc] the channel
      // telescopes to the full perturbation.
      if (!disable_nucleation && num_fault_total > 0)
      {
         ApplyNucleationIncremental_TPV104(dof_data, fault_coords,
                                           t + dt_step, dt_step);
      }

      // One-shot ADER predictor-corrector.  AdvanceADER runs the bulk
      // wave update + the fault-face Riemann solve (through
      // FaultFaceFlux::EvaluateADERTotal); the solve reads DOFData.psi
      // and DOFData.tau*_nuc as set above, and writes V1/V2/slip_rate/
      // tau*_corr/sigma_n_corr back onto DOFData.
      wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
      Q.Swap(Q_new);
      t += dt_step;

#ifdef SEAS_DIAG_FAULT_FLUX
      // C-2C BULK-DELTA: per-macro-step delta of bulk Q at the diag
      // DOFs.  Captures the NET change Q^{n+1} - Q^n at the two diag
      // tets' face DOFs after one full ADER step (predictor + corrector
      // including fault Riemann + bulk wave op).  Lets the post-run
      // analyzer correlate per-step asymmetry growth against C-1n.
      if (diag_elem_plus >= 0 && diag_elem_minus >= 0 &&
          diag_face_dof_plus >= 0 && diag_face_dof_minus >= 0)
      {
         static std::vector<real_t> prev_Qp(NUM_STATE, 0.0);
         static std::vector<real_t> prev_Qm(NUM_STATE, 0.0);
         static bool have_prev = false;

         const real_t *Qd = Q.GetData();
         const int dof_off_p = diag_elem_plus  * wave.GetNDof();
         const int dof_off_m = diag_elem_minus * wave.GetNDof();
         real_t Qp[NUM_STATE], Qm[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Qp[c] = Qd[c*ndof_total + dof_off_p + diag_face_dof_plus];
            Qm[c] = Qd[c*ndof_total + dof_off_m + diag_face_dof_minus];
         }
         real_t dQp[NUM_STATE], dQm[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            dQp[c] = have_prev ? (Qp[c] - prev_Qp[c]) : 0.0;
            dQm[c] = have_prev ? (Qm[c] - prev_Qm[c]) : 0.0;
            prev_Qp[c] = Qp[c];
            prev_Qm[c] = Qm[c];
         }
         have_prev = true;
         std::fprintf(stderr,
            "[C-2C BULK-DELTA] rank=%d t=%.4e  "
            "Q+_SXX=%+.4e Q-_SXX=%+.4e diff_SXX=%+.4e  "
            "dQ+_SXX=%+.4e dQ-_SXX=%+.4e d_diff_SXX=%+.4e  "
            "Q+_SYY=%+.4e Q-_SYY=%+.4e diff_SYY=%+.4e  "
            "Q+_SZZ=%+.4e Q-_SZZ=%+.4e diff_SZZ=%+.4e\n",
            g_seas_my_rank, t,
            Qp[SXX], Qm[SXX], Qp[SXX] - Qm[SXX],
            dQp[SXX], dQm[SXX], dQp[SXX] - dQm[SXX],
            Qp[SYY], Qm[SYY], Qp[SYY] - Qm[SYY],
            Qp[SZZ], Qm[SZZ], Qp[SZZ] - Qm[SZZ]);
      }
#endif

      // ψ update + slip accumulation per plan §3.3 / §4.10 Step 9.
      // Replaces TPV102's forward-Euler on AgingLawPsi with the FVW
      // analytic exponential-relaxation update.  V_w is per-QP; `a` is
      // per-QP and already stored in DOFData (set by InitializeFaultDOFs).
      //
      // R7-007 disclosure: plan §3.12 mandates PER-SUB-STEP ψ integration
      // ("follow exactly SeisSol did").  On this driver path ψ is
      // integrated ONCE per macro-step with V = dof_data[i].slip_rate
      // (ADER-averaged V over the whole dt_step).  The per-sub-step
      // cadence requires wiring Tpv104SubStepIterator into the time
      // loop (blocked by R7-001 option a), which in turn requires
      // exposing per-sub-step I± from wave_operator.inl (extreme-care
      // no-touch).  Under a rapidly-changing V the macro-step ψ deviates
      // from the per-sub-step result by O(dt_macro²); Phase-3 probe 2
      // against SeisSol will quantify the gap.
      for (int i = 0; i < num_fault_total; ++i)
      {
         dof_data[i].psi = UpdateStateAnalyticSlipLawSRW(
            psi_n[i],
            dof_data[i].slip_rate,
            dof_data[i].Dc,
            dt_step,
            V_w[i],
            dof_data[i].a,
            TPV104Params::b,
            TPV104Params::V0,
            TPV104Params::f0,
            TPV104Params::f_w);
         dof_data[i].slip1 += dof_data[i].V1 * dt_step;
         dof_data[i].slip2 += dof_data[i].V2 * dt_step;
      }

      // V_max tracking.
      real_t V_max_local = 0.0;
      for (int i = 0; i < num_fault_total; ++i)
      {
         V_max_local = std::max(V_max_local, dof_data[i].slip_rate);
      }
#ifdef MFEM_USE_MPI
      real_t V_max_step;
      MPI_Allreduce(&V_max_local, &V_max_step, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#else
      real_t V_max_step = V_max_local;
#endif
      V_max_global = std::max(V_max_global, V_max_step);

      // Output cadence.
      if (step % output_interval == 0 || step == nsteps - 1)
      {
         station_writer.WriteStep(t, dof_data);
         surface_writer.WriteStep(t, Q);
         if (rank == 0)
         {
            std::cout << "Step " << step << "/" << nsteps
                      << ", t = " << t << " s"
                      << ", V_max = " << V_max_step << " m/s\n";
         }
      }

      // ParaView snapshot — has its own internal schedule (matches
      // --paraview-dt / --paraview-every / output_dt), so we always
      // call and let pv_out / pv_bulk_out's PeekShouldWrite gate the I/O.
      paraview_write(step + 1, t, V_max_step);

      if (step == 0)
      {
         wave.VerifySharedFaultDOFDataConsistency();
      }

      // NaN tripwire.
      {
         real_t local_nan = std::isnan(Q.Norml2()) ? 1.0 : 0.0;
         real_t global_nan = local_nan;
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&local_nan, &global_nan, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
         if (global_nan > 0.0)
         {
            std::cerr << "ERROR: NaN detected at step " << step
                      << ", t = " << t << " s (rank " << rank << ")\n";
            Tpv104SubStepIterator::CloseAllProbeFiles();
#ifdef MFEM_USE_MPI
            MPI_Finalize();
#endif
            return 1;
         }
      }

      if (debug_qnorm && step % output_interval == 0 && rank == 0)
      {
         std::cout << "  [qnorm] ||Q||_2 = " << Q.Norml2() << "\n";
      }
   }

   // -----------------------------------------------------------------------
   // 9. Summary
   // -----------------------------------------------------------------------
   station_writer.Flush();
   station_writer.Close();
   surface_writer.Flush();
   surface_writer.Close();
   Tpv104SubStepIterator::CloseAllProbeFiles();

   if (rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "TPV104 run complete.\n";
      std::cout << "  tfinal  = " << tfinal << " s\n";
      std::cout << "  steps   = " << nsteps << "\n";
      std::cout << "  V_max   = " << V_max_global << " m/s\n";
      std::cout << "  output  = " << output_dir << "/\n";
      std::cout << "========================================\n";
   }

#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return 0;
}
