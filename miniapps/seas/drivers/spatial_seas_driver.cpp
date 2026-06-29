// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// seas_spatial_seas_driver — Phase 0 of
// document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md
//
// Quasi-dynamic (QD) SEAS driver for spatially-varying material, pre-stress,
// and rate-and-state friction.  QD = quasi-static elasticity (no inertia /
// wave propagation) + radiation-damping term (eta * V) on the fault, advanced
// by an adaptive Dormand-Prince RK45 integrator.  This is the long-timescale
// (interseismic + coseismic) sibling of the single-event dynamic-rupture
// driver `spatial_dyn_driver.cpp`; the two share the spatial config schema
// (spatial/code/spatial_friction.hpp) but NOTHING from the dynamic/ wave stack.
//
// ---------------------------------------------------------------------------
// PHASE 0 SKELETON.  This is the empty-but-wired entry point: it initialises
// MPI, parses the CLI + TOML config, echoes a rank-0 banner + derived
// parameters, and exits on --dry-run BEFORE any operator construction.  The
// elasticity/fault operators, the RK45 time loop, and the ParaView/checkpoint
// I/O land in Phases 2-9.  A non-dry-run invocation therefore prints a clear
// "not yet implemented" notice and exits 0.
//
// CLI surface (Phase 0).  Generic scaffolding mirrored from
// spatial_dyn_driver.cpp:137-200 (CLI helpers) and :504-756 (config merge):
//     --config PATH.toml      (required)
//     --dry-run               parse + echo, then exit before construction
//     --mesh PATH             override [mesh].path
//     --tfinal SECONDS        override [time].tfinal
//     --restart PREFIX        checkpoint prefix to restart from (wired later)
//     --checkpoint-every N    checkpoint cadence in accepted steps
//     --output-dir DIR        override [output].output_dir
//     --paraview-volume MODE  vtu | off  (QD uses VTU/PVD, never vtkhdf)
//     --paraview-fault  MODE  vtu | off
//     --paraview-volume-dt S  volume snapshot cadence [s]
//     --paraview-fault-dt  S  fault  snapshot cadence [s]
// The dynamic-only flags (--ader-order, --mixed-flux, --pml*,
// --time-integrator) are intentionally absent: QD has no inertial CFL, no
// flux-mode choice, and a single (RK45) integrator.
//
// CLI override convention (matches the spatial driver): an empty-string flag
// or a negative numeric sentinel means "keep the TOML value"; any value given
// on the CLI overrides the TOML (CLI wins, applied after the TOML load).
//
// ---------------------------------------------------------------------------
// Build (from a normal checkout, miniapps/seas):
//     conda activate mfem-dev
//     make seas_spatial_seas_driver
//
// Build-from-worktree caveat (project memory).  Building under
// .claude/worktrees/<name>/miniapps/seas needs the MFEM build/lib/include
// dirs pointed at the MAIN checkout, because the worktree root has no
// libmfem.a (Makefile's `MFEM_DIR ?= ../..` would otherwise fail to link):
//     MAIN=/Users/<you>/projects/seas-mfem-spatial-dyn-driver
//     make MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN \
//          MFEM_LIB_DIR=$MAIN seas_spatial_seas_driver
// The gitignored extern/toml11 must also be symlinked from the main checkout
// (else the TOML parse path silently drops out / fails to find toml.hpp).

#include "mfem.hpp"

#include "../spatial/code/spatial_friction.hpp"
// Phase 2b: heterogeneous material.  MaterialField + DepthProfile1DMaterial come
// via spatial_friction.hpp -> dynamic/heterogeneous_material.hpp (a field
// abstraction, explicitly reusable per the plan — NOT wave physics).  The
// velocity sidecar bundle lives in the spatial scaffolding.
#include "../spatial/code/spatial_velocity.hpp"

// Phase 2: the elasticity domain operator (header-only template; transitively
// pulls in BoundaryConfig + MakeBP5DirichletFunc, DGMethod, LinearElastic,
// DomainConfig, SolverType, BCMode).  seas_driver links with only its own
// object, so this adds no link dependency to the QD driver.
#include "../domain/elasticity_operator.hpp"

// Phase 3: fault geometry + spatial per-DOF rate-state wiring (all header-only).
#include "../fault/fault_geometry.hpp"
#include "../config/bp5_params.hpp"
#include "../common/mpi_context.hpp"

// Phase 5: the quasi-dynamic coupling stack (the BP5/QD solver stack — NOT the
// dynamic-rupture wave code, which the plan forbids).  All header-only; the
// templated ComputeParams<StressSource> body + the stress-source primitives.
#include "../fault/rate_state_fault.hpp"
#include "../fault/fault_geometry_safs_templated.inl"
#include "../friction/dieterich_ruina.hpp"
#include "../friction/state_evolution.hpp"
#include "../friction/slip_law_srw_psi.hpp"
#include "../solver/seas_operator.hpp"
#include "../solver/time_stepper.hpp"
#include "../spatial/code/spatial_stress.hpp"

// Phase 6: QD-native checkpoint/restart (header-only; NO dynamic/ dep — R-009).
#include "../io/seas_qd_checkpoint.hpp"
// Phase 6: BP5 SCEC station/probe output (header-only; pulls bp5_benchmark_output
// for Station + DefaultStations).  Used only for BP5-parity runs.
#include "../io/bp5_parallel_output.hpp"

#ifdef SEAS_USE_TOML
#include <toml.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

// --------------------------------------------------------------------------
// Small CLI parsing helpers (same convention as spatial_dyn_driver.cpp:137).
// --------------------------------------------------------------------------
namespace
{

std::string GetStringArg(int argc, char *argv[], const char *flag,
                         const std::string &default_val)
{
   for (int i = 1; i < argc - 1; ++i)
   {
      if (std::string(argv[i]) == flag) { return argv[i + 1]; }
   }
   return default_val;
}

real_t GetRealArg(int argc, char *argv[], const char *flag, real_t default_val)
{
   const std::string val = GetStringArg(argc, argv, flag, "");
   if (val.empty()) { return default_val; }
   return std::stod(val);
}

int GetIntArg(int argc, char *argv[], const char *flag, int default_val)
{
   const std::string val = GetStringArg(argc, argv, flag, "");
   if (val.empty()) { return default_val; }
   return std::stoi(val);
}

bool HasFlag(int argc, char *argv[], const char *flag)
{
   for (int i = 1; i < argc; ++i)
   {
      if (std::string(argv[i]) == flag) { return true; }
   }
   return false;
}

const char *LawName(spatial::FrictionLawKind law)
{
   return (law == spatial::FrictionLawKind::SlipWeakening)
          ? "slip_weakening" : "rate_state";
}

}  // namespace

// =========================================================================
// main
// =========================================================================

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);
#else
   int rank = 0, nprocs = 1;
#endif

   // -----------------------------------------------------------------
   // 1.  CLI parse — config-driven; CLI overrides are merged after the
   //     TOML load (later wins).
   // -----------------------------------------------------------------
   const std::string config_path = GetStringArg(argc, argv, "--config", "");
   if (config_path.empty())
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: --config PATH.toml is required.\n"
                   << "  Usage: seas_spatial_seas_driver --config "
                   << "<file>.toml [OPTIONS]\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 2;
   }

   const bool dry_run = HasFlag(argc, argv, "--dry-run");

   const std::string cli_mesh       = GetStringArg(argc, argv, "--mesh", "");
   // Phase 6 checkpoint/restart: --restart PREFIX restores (t,dt,step)+state;
   // --checkpoint-every N writes a checkpoint every N accepted steps.
   const std::string restart_prefix = GetStringArg(argc, argv, "--restart", "");
   const std::string cli_output_dir =
      GetStringArg(argc, argv, "--output-dir", "");
   const std::string cli_pv_volume =
      GetStringArg(argc, argv, "--paraview-volume", "");
   const std::string cli_pv_fault =
      GetStringArg(argc, argv, "--paraview-fault", "");

   const real_t cli_tfinal      = GetRealArg(argc, argv, "--tfinal", -1.0);
   const real_t cli_pv_vol_dt   = GetRealArg(argc, argv, "--paraview-volume-dt", -1.0);
   const real_t cli_pv_fault_dt = GetRealArg(argc, argv, "--paraview-fault-dt", -1.0);
   const int    cli_checkpoint_every =
      GetIntArg(argc, argv, "--checkpoint-every", -1);

   // -----------------------------------------------------------------
   // 2.  TOML load — owns all defaults; CLI then overrides.
   // -----------------------------------------------------------------
   // R-001: LoadSpatialFrictionConfig reports a missing/unreadable file with
   // MFEM_ABORT (not a C++ exception), so the catch below cannot turn that
   // common case into the documented "exit 2".  Pre-check readability here
   // (QD-driver-local; the shared parser is untouched) so a bad --config path
   // gives a clean rank-0 message and exit 2.  A readable-but-malformed TOML
   // still aborts inside the shared parser — matching spatial_dyn_driver.
   {
      std::ifstream probe(config_path);
      if (!probe.good())
      {
         if (rank == 0)
         {
            std::cerr << "ERROR: --config '" << config_path
                      << "' is not readable.\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 2;
      }
   }

   spatial::SpatialFrictionConfig cfg;
   try
   {
      cfg = spatial::LoadSpatialFrictionConfig(config_path);
   }
   catch (...)
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: failed to load TOML config '" << config_path
                   << "'.\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 2;
   }

   // Merge CLI overrides into cfg (later wins; sentinels keep the TOML value).
   if (!cli_mesh.empty())       { cfg.mesh.path = cli_mesh; }
   if (cli_tfinal > 0.0)        { cfg.time.tfinal = cli_tfinal; }
   if (!cli_output_dir.empty()) { cfg.output.output_dir = cli_output_dir; }
   if (!cli_pv_volume.empty())  { cfg.output.paraview_volume = cli_pv_volume; }
   if (!cli_pv_fault.empty())   { cfg.output.paraview_fault = cli_pv_fault; }
   if (cli_pv_vol_dt > 0.0)     { cfg.output.paraview_volume_dt = cli_pv_vol_dt; }
   if (cli_pv_fault_dt > 0.0)   { cfg.output.paraview_fault_dt = cli_pv_fault_dt; }
   if (cli_checkpoint_every > 0)
   {
      cfg.output.checkpoint_every_steps = cli_checkpoint_every;
   }

   // Phase 6 restart safety: refuse if --output-dir resolves to the SAME
   // directory as the parent of --restart (a restart would overwrite the
   // checkpoint it is reading).  Mirrors spatial_dyn_driver.cpp:852-894.
   if (!restart_prefix.empty())
   {
      namespace fs = std::filesystem;
      try
      {
         fs::path restart_dir = fs::path(restart_prefix).parent_path();
         if (restart_dir.empty()) { restart_dir = "."; }
         if (fs::weakly_canonical(restart_dir)
             == fs::weakly_canonical(fs::path(cfg.output.output_dir)))
         {
            if (rank == 0)
            {
               std::cerr << "ERROR: --output-dir (" << cfg.output.output_dir
                         << ") resolves to the SAME directory as the parent of "
                         << "--restart (" << restart_dir.string() << "). Pick a "
                         << "DIFFERENT --output-dir for the restarted run.\n";
            }
#ifdef MFEM_USE_MPI
            MPI_Finalize();
#endif
            return 3;
         }
      }
      catch (const fs::filesystem_error &e)
      {
         if (rank == 0)
         {
            std::cerr << "ERROR: failed to canonicalise --restart / "
                      << "--output-dir paths: " << e.what() << "\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 3;
      }
   }

   // R-702: the BP5 prestress (tau0_vec) uses bp5_params.a(x2,x3) internally, so
   // the rate-state a MUST also come from bp5_params.  Refuse a half-set config
   // (one bp5_analytic flag without the other) here — a clear error instead of a
   // confusing SetInitialCondition equilibrium failure (R-007).
   {
      const bool rs_bp5 = cfg.rate_state.has_value()
                          && cfg.rate_state->bp5_analytic;
      MFEM_VERIFY(rs_bp5 == cfg.stress.bp5_analytic,
                  "spatial_seas: [friction.rate_state].bp5_analytic ("
                  << rs_bp5 << ") and [stress].bp5_analytic ("
                  << cfg.stress.bp5_analytic << ") must be set TOGETHER — the "
                  "BP5 prestress uses bp5_params.a(x2,x3), so the rate-state a "
                  "must too (R-007).  Mixing resolver a with bp5_analytic "
                  "prestress breaks equilibrium.");
   }

   // Strong rate weakening (state_evolution="slip_law_strong_rate_weakening")
   // IS wired into the QD fault operator (SlipLawSRWPsi via the per-DOF _SRW
   // route; see §5b below + PLAN_spatial_seas_srw_qd_2026-06-20.md).
   // But bp5_analytic prestress bypasses the resolver (it never fills per-DOF
   // V_w) and is BP5 aging-law only, so it cannot be combined with SRW — reject
   // that contradiction here with a clear message (else the operator ctor aborts
   // later on a missing-V_w, which is harder to diagnose).
   if (cfg.rate_state.has_value()
       && cfg.rate_state->bp5_analytic
       && cfg.rate_state->state_evolution
              == spatial::StateEvolutionKind::SlipLawStrongRateWeakening)
   {
      MFEM_ABORT("spatial_seas: [friction.rate_state].bp5_analytic=true is BP5 "
                 "aging-law only and bypasses the resolver (no per-DOF V_w), so "
                 "state_evolution=\"slip_law_strong_rate_weakening\" is "
                 "unsupported with it.  Use a resolver stress source "
                 "(constant_tensor / fault_local_prestress / depth_proportional) "
                 "for an SRW run, or set state_evolution=\"aging_law\".");
   }

   // -----------------------------------------------------------------
   // 3.  Rank-0 banner + derived-parameter echo.
   // -----------------------------------------------------------------
   if (rank == 0)
   {
      // Derived: QD runs span many seismic cycles; report tfinal in years
      // too (365.25-day Julian year) so the long-timescale scale is obvious.
      const real_t seconds_per_year = 365.25 * 24.0 * 3600.0;
      const real_t tfinal_years = cfg.time.tfinal / seconds_per_year;

      // R-002: the shared OutputSpec default for the ParaView modes is "hdf5",
      // but the QD driver writes VTU/PVD (plan I/O scope 2026-06-03) — Phase 6
      // forces vtu.  Annotate "hdf5" in the banner so an operator is not
      // misled into expecting .vtkhdf output from the QD driver.
      auto pv_mode_echo = [](const std::string &m)
      {
         return (m == "hdf5") ? std::string("hdf5 (Phase 6 will force vtu)") : m;
      };
      const std::string pv_volume_eff = pv_mode_echo(cfg.output.paraview_volume);
      const std::string pv_fault_eff  = pv_mode_echo(cfg.output.paraview_fault);

      std::cout << "================================================\n"
                << "seas_spatial_seas_driver — Phase 2 (mesh + elasticity op)\n"
                << "  quasi-dynamic SEAS (quasi-static + radiation damping)\n"
                << "================================================\n"
                << "config:           " << config_path << "\n"
                << "schema version:   " << cfg.schema_version << "\n"
                << "description:      " << cfg.description << "\n"
                << "mesh:             " << cfg.mesh.path << "\n"
                << "fe order:         " << cfg.mesh.order << "\n"
                << "law:              " << LawName(cfg.law) << "\n"
                << "tfinal:           " << cfg.time.tfinal << " s\n"
                << "t_initial:        " << cfg.time.t_initial << " s\n";
      if (cfg.time.dt_initial < 0.0)
      {
         std::cout << "dt_initial:       auto\n";
      }
      else
      {
         std::cout << "dt_initial:       " << cfg.time.dt_initial << " s\n";
      }
      std::cout << "dt_max:           " << cfg.time.dt_max << " s\n"
                << "output dir:       " << cfg.output.output_dir << "\n"
                << "paraview volume:  " << pv_volume_eff << "\n"
                << "paraview fault:   " << pv_fault_eff << "\n"
                << "checkpoint every: " << cfg.output.checkpoint_every_steps
                << " steps\n"
                << "restart:          "
                << (restart_prefix.empty() ? "(none)" : restart_prefix) << "\n"
                << "dry-run:          " << (dry_run ? "yes" : "no") << "\n"
                << "ranks:            " << nprocs << "\n"
                << "------------------------------------------------\n"
                << "[derived] tfinal: " << tfinal_years << " yr\n"
                << "================================================\n";
   }

   // -----------------------------------------------------------------
   // 3b. Dynamic-only key warning (Phase 1; R-001 / plan §C).  Printed
   //     AFTER the banner (R-105) so the config echo provides context.
   //     The QD driver ignores the dynamic-rupture [numerics] knobs
   //     (ader_order / mixed_flux / interior_flux / use_pml /
   //     time_integrator).  This warning lives HERE — never inside the
   //     shared ParseSpatialFrictionConfigString — so the dynamic driver
   //     does not warn on its own valid [numerics] keys.  Detection
   //     re-inspects the raw TOML for KEY PRESENCE (the parsed struct
   //     cannot distinguish "present" from "default"); the shared parser
   //     is left completely untouched.  Emitted ONCE on rank 0.
   // -----------------------------------------------------------------
#ifdef SEAS_USE_TOML
   if (rank == 0)
   {
      try
      {
         const toml::value root = toml::parse(config_path);
         if (root.contains("numerics"))
         {
            const auto &num = root.at("numerics");
            const char *dyn_keys[] = { "ader_order", "mixed_flux",
                                       "interior_flux", "use_pml",
                                       "time_integrator" };
            std::vector<std::string> present;
            for (const char *k : dyn_keys)
            {
               if (num.contains(k)) { present.emplace_back(k); }
            }
            if (!present.empty())
            {
               std::cerr << "WARNING: [numerics] key(s)";
               for (const auto &k : present) { std::cerr << " '" << k << "'"; }
               std::cerr << " are ignored by spatial_seas (quasi-dynamic) "
                         << "driver.\n";
            }
         }
      }
      catch (...)
      {
         // Re-parse is best-effort: LoadSpatialFrictionConfig already
         // succeeded above, so a failure here is not fatal — skip the
         // warning rather than abort.
      }
   }
#endif

   // -----------------------------------------------------------------
   // 4.  Mesh, boundary, material, elasticity domain operator (Phase 2,
   //     homogeneous-first).  Construction runs for BOTH --dry-run and a
   //     normal run; --dry-run then prints the fault-DOF counts and exits
   //     (section 5), before the (not-yet-implemented) time loop.
   // -----------------------------------------------------------------
#ifndef MFEM_USE_MPI
#  error "spatial_seas_driver requires MFEM_USE_MPI=YES."
#endif

   // 4.0 Material is built in 4.4 below (Phase 2b).  The Phase-2 fail-fast
   //     heterogeneous abort was REMOVED here: heterogeneous material (depth
   //     profile / velocity sidecar) is now supported via the IP coefficient
   //     operator ctor (Phase 2b Stage 1b).  Building the sidecar MaterialField
   //     needs the mesh, so it is deferred to 4.4 (after the mesh load); a
   //     sidecar load failure there aborts loudly (no silent fallback).

   // 4.1 Mesh (reuse spatial_dyn_driver.cpp:899-911).
   Mesh smesh(cfg.mesh.path.c_str(), 1, 1);
   MFEM_VERIFY(smesh.Dimension() == 3,
               "spatial_seas_driver: only 3D meshes supported; got dim="
               << smesh.Dimension());
   ParMesh pmesh(comm, smesh);
   smesh.Clear();
   pmesh.SetCurvature(cfg.mesh.order);

   // 4.2 Plate-loading rate Vp: [time].plate_rate_vp (>0), else the
   //     rate-state initial slip rate, else the BP5 default 1e-9 m/s.
   real_t plate_rate = cfg.time.plate_rate_vp;
   if (plate_rate <= 0.0)
   {
      plate_rate = cfg.rate_state.has_value()
                   ? cfg.rate_state->V_init_default : 1.0e-9;
   }

   // 4.3 BoundaryConfig.  QD needs Dirichlet plate loading on the far-field
   //     walls — the operator aborts on a fault-without-Dirichlet mesh
   //     (elasticity_operator_setup.inl:86).  Source the attrs from
   //     [boundary] with the BP5 convention (fault=3, natural=1, dirichlet=5)
   //     as the fallback, and the BP5 far-field loading function (SAF
   //     loading is Phase 8 — see the plan's open questions).
   BoundaryConfig bc;
   bc.fault_attr = (cfg.boundary.fault_attr > 0) ? cfg.boundary.fault_attr : 3;
   bc.natural_attrs = !cfg.boundary.natural_attrs.empty()
                      ? std::set<int>(cfg.boundary.natural_attrs.begin(),
                                      cfg.boundary.natural_attrs.end())
                      : std::set<int>{1};
   // R-201: Dirichlet plate loading is MANDATORY for the QD elasticity solve.
   //   - [boundary] present (fault_attr > 0) ⇒ require a non-empty
   //     dirichlet_attrs, with a clear QD-specific abort (not a downstream
   //     "Dirichlet attr 5 not in mesh").
   //   - [boundary] absent ⇒ fall back to the BP5 convention (attr 5).
   if (cfg.boundary.fault_attr > 0)
   {
      MFEM_VERIFY(!cfg.boundary.dirichlet_attrs.empty(),
                  "spatial_seas_driver: [boundary] is present but "
                  "dirichlet_attrs is empty — the quasi-dynamic elasticity "
                  "problem requires far-field Dirichlet plate-loading walls "
                  "(set [boundary].dirichlet_attrs).");
      bc.dirichlet_attrs = std::set<int>(cfg.boundary.dirichlet_attrs.begin(),
                                         cfg.boundary.dirichlet_attrs.end());
   }
   else
   {
      bc.dirichlet_attrs = std::set<int>{5};  // BP5 default (no [boundary] block)
   }
   // R-206: absorbing_attrs is a dynamic-wave concept; the quasi-static
   //   elasticity operator does not consume it (unmarked faces default to
   //   natural / traction-free).  Forwarded for completeness; empty for BP5.
   bc.absorbing_attrs = std::set<int>(cfg.boundary.absorbing_attrs.begin(),
                                      cfg.boundary.absorbing_attrs.end());

   // Far-field plate-loading function (Phase 8 / REVIEW R-001).
   //   "bp5"            : u_X = sgn(absolute y)*Vp*t/2 (BP5; default, byte-
   //                      identical) — valid ONLY for an origin-centred fault.
   //   "saf_recenter_y" : recentre the antisymmetry on the mesh y mid-plane y0
   //                      so a UTM/non-origin SAF box actually loads (on such a
   //                      mesh "bp5" degenerates to a UNIFORM +Vp/2 — no
   //                      differential plate motion).  All ranks share the
   //                      config value, so the collective reduce is balanced.
   if (cfg.boundary.plate_loading == "saf_recenter_y")
   {
      // Box bbox (curved-aware, R-004): diagnostics + fallback split plane.
      // GetBoundingBox is per-rank, so reduce across ranks (all ranks take this
      // config-driven branch, so the collectives are balanced).
      Vector pmin, pmax;
      pmesh.GetBoundingBox(pmin, pmax, /*ref=*/2);
      real_t ylo = pmin(1), yhi = pmax(1);
      MPI_Allreduce(MPI_IN_PLACE, &ylo, 1,
                    mfem::MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
      MPI_Allreduce(MPI_IN_PLACE, &yhi, 1,
                    mfem::MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
      const real_t y_box = 0.5 * (ylo + yhi);

      // R-005: split the plate-loading antisymmetry at the FAULT (where the two
      // plates meet), not the geometric box centre — correct even for a box that
      // is not symmetric about the fault.  Average y over the fault boundary
      // faces (attr == fault_attr; the SAFS fault is node-duplicated -> boundary
      // elements, cf. domain/seas_boundary_tags.hpp).  Fall back to the box
      // mid-plane + WARN if the mesh exposes no fault bdr faces.
      real_t fy_sum = 0.0;
      long   fy_cnt = 0;
      Array<int> fverts;
      for (int be = 0; be < pmesh.GetNBE(); ++be)
      {
         if (pmesh.GetBdrAttribute(be) != bc.fault_attr) { continue; }
         pmesh.GetBdrElementVertices(be, fverts);
         if (fverts.Size() == 0) { continue; }
         real_t yc = 0.0;
         for (int v = 0; v < fverts.Size(); ++v)
         { yc += pmesh.GetVertex(fverts[v])[1]; }
         fy_sum += yc / fverts.Size();
         fy_cnt += 1;
      }
      MPI_Allreduce(MPI_IN_PLACE, &fy_sum, 1,
                    mfem::MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
      MPI_Allreduce(MPI_IN_PLACE, &fy_cnt, 1, MPI_LONG, MPI_SUM, comm);

      real_t y0 = y_box;
      const char *y0_src = "box mid-plane (fallback)";
      if (fy_cnt > 0)
      {
         y0 = fy_sum / static_cast<real_t>(fy_cnt);   // fault y-centroid
         y0_src = "fault centroid";
         const real_t half_w = 0.5 * std::max(yhi - ylo, real_t(1.0));
         if (rank == 0 && std::abs(y0 - y_box) > 0.05 * half_w)
         {
            std::cerr << "[spatial_seas] NOTE: fault y-centroid (" << y0
                      << ") is off the box mid-plane (" << y_box << ") by "
                      << std::abs(y0 - y_box) << " m (>5% of the box half-width "
                      << half_w << ") — splitting the loading at the fault.\n";
         }
      }
      else if (rank == 0)
      {
         std::cerr << "[spatial_seas] WARNING: saf_recenter_y found no fault bdr "
                   << "faces (attr " << bc.fault_attr << "); splitting at the box "
                   << "mid-plane y0=" << y_box << " instead.\n";
      }
      bc.default_dirichlet_func = MakeSAFDirichletFunc(plate_rate, y0);
      if (rank == 0)
      {
         std::cout << "[spatial_seas] plate loading: saf_recenter_y, y0=" << y0
                   << " m (" << y0_src << "; box y in [" << ylo << ", " << yhi
                   << "]); per-side +/-Vp/2 = +/-" << 0.5 * plate_rate
                   << " m/s about the fault.\n";
      }
   }
   else
   {
      bc.default_dirichlet_func = MakeBP5DirichletFunc(plate_rate);
   }

   // 4.4 Material (Phase 2b).  Build the MaterialField from config:
   //   - [material].kind="constant" (default)  → Mode::Constant (LinearElastic).
   //   - [material].kind="depth_profile_1d"    → Mode::Coefficient (per-qp).
   //   - [velocity].use_sidecar=true           → Mode::Coefficient (sidecar).
   //   The heterogeneous (Mode::Coefficient) field is consumed per quadrature
   //   point by the DGMethod::IP integrators via the operator's coefficient ctor
   //   (Stage 1b) AND by ResolveRateState's auto-η (material.EvalAt, Phase 3).
   //   The owning wrappers (vel_bundle / depth_profile_wrapper) are declared
   //   HERE so they outlive domain_ptr (4.7) — the operator stores NON-OWNING
   //   Coefficient handles into them (PLAN §Phase 2b Edge Cases).
   const real_t mat_lambda = cfg.material_fallback.lambda;
   const real_t mat_mu     = cfg.material_fallback.mu;
   LinearElastic le(mat_lambda, mat_mu);   // used only on the constant path

   std::unique_ptr<spatial::SpatialVelocityBundle> vel_bundle;
   std::unique_ptr<DepthProfile1DMaterial>         depth_profile_wrapper;
   MaterialField material = MaterialField::MakeConstant(
      mat_lambda, mat_mu, cfg.material_fallback.rho);

   // No --no-sidecar-material CLI on this driver (unlike spatial_dyn_driver):
   // the TOML gate alone decides.  depth_profile_1d takes precedence over the
   // sidecar.  R-502: a sidecar_hdf5 material with use_sidecar=false is a
   // contradictory config — abort rather than silently fall back to constant.
   MFEM_VERIFY(!(cfg.material.kind == spatial::MaterialKind::SidecarHDF5
                 && !cfg.velocity.use_sidecar),
               "spatial_seas_driver: [material].kind=\"sidecar_hdf5\" requires "
               "[velocity].use_sidecar=true (got use_sidecar=false).");
   // R-501: honor an EXPLICIT kind=constant — do NOT attempt a sidecar just
   // because use_sidecar defaults to true.  After this gate, sidecar_requested
   // is true only for an explicit sidecar_hdf5 kind with use_sidecar.
   const bool sidecar_requested =
      cfg.velocity.use_sidecar
      && cfg.material.kind != spatial::MaterialKind::DepthProfile1D
      && cfg.material.kind != spatial::MaterialKind::Constant;

   if (cfg.material.kind == spatial::MaterialKind::DepthProfile1D)
   {
      depth_profile_wrapper = MakeDepthProfile1DMaterial(
         cfg.material.profile_layers, cfg.material.depth_axis);
      material = depth_profile_wrapper->field;
      if (rank == 0)
      {
         std::cout << "[spatial_seas] material: depth_profile_1d (axis='"
                   << cfg.material.depth_axis << "', "
                   << cfg.material.profile_layers.size() << " layers)\n";
      }
   }
   else if (sidecar_requested)
   {
      try
      {
         vel_bundle = std::make_unique<spatial::SpatialVelocityBundle>(
            spatial::LoadSpatialVelocityBundle(cfg.velocity, pmesh));
         material = vel_bundle->MakeMaterialField();
         if (rank == 0)
         {
            std::cout << "[spatial_seas] material: velocity sidecar loaded "
                      << "(Mode::Coefficient)\n";
         }
      }
      catch (const std::exception &e)
      {
         if (rank == 0)
         {
            std::cerr << "ERROR: spatial_seas_driver failed to load velocity "
                      << "sidecar: " << e.what() << "\n"
                      << "  Set [velocity].use_sidecar=false (or "
                      << "[material].kind=\"constant\") to use the "
                      << "[material_constant_fallback] block.\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 4;
      }
   }
   else if (rank == 0)
   {
      std::cout << "[spatial_seas] material: constant (lambda=" << mat_lambda
                << " mu=" << mat_mu << " rho=" << cfg.material_fallback.rho
                << ")\n";
   }

   // 4.5 DomainConfig from the QD [solver] knobs (Phase 1).
   DomainConfig domain_config;
   domain_config.check_residual = cfg.solver.residual_check;
   domain_config.blr_tol        = cfg.solver.blr_tol;
   // Phase 4: Krylov/AMG knobs for the CG_AMG / GMRES_AMG paths.
   domain_config.ksp_rtol        = cfg.solver.ksp_rtol;
   domain_config.ksp_atol        = cfg.solver.ksp_atol;
   domain_config.ksp_maxit       = cfg.solver.ksp_maxit;
   domain_config.amg_print_level = cfg.solver.amg_print_level;
   domain_config.amg_elasticity_options = cfg.solver.amg_elasticity_options;
   domain_config.amg_relax_type         = cfg.solver.amg_relax_type;
   domain_config.amg_aggressive_levels  = cfg.solver.amg_aggressive_levels;

   // 4.6 Fault dimensions.  The SpatialFrictionConfig schema has no
   //     fault-depth/length fields, and Wf_/lf_ do NOT enter stiffness
   //     assembly (stored-only getters + one BC-time guard,
   //     elasticity_operator_assembly.inl:1077).  Use the BP5 reference
   //     dimensions as documented placeholders; a config-sourced
   //     [fault_geometry] dimension block is deferred.
   const real_t fault_depth  = 40.0e3;   // BP5 Wf
   const real_t fault_length = 100.0e3;  // BP5 lf

   // 4.7 Construct.  dg_method = IP (the BP5 convention; the schema carries no
   //     dg_method field).  Constant material → the LinearElastic scalar ctor;
   //     heterogeneous (Mode::Coefficient) → the per-qp coefficient ctor
   //     (Phase 2b Stage 1b; IP-only).  Deferred construction via a unique_ptr
   //     so the ctor is selected by material.mode; `domain` aliases it so all
   //     downstream wiring (4b/5/6) is unchanged.  Construction triggers the
   //     stiffness assembly + the configured solver setup.
   std::unique_ptr<ElasticityDomainOperator<ParMesh>> domain_ptr;
   if (material.mode == MaterialField::Mode::Constant)
   {
      domain_ptr = std::make_unique<ElasticityDomainOperator<ParMesh>>(
         pmesh, cfg.mesh.order, le,
         plate_rate, fault_depth, fault_length,
         bc, DGMethod::IP, spatial::ParseQDSolverType(cfg.solver.type),
         domain_config);
   }
   else
   {
      // Heterogeneous: must be Mode::Coefficient (per-qp).  Mode::GridFunction
      // would require projecting onto the elasticity space (different FE space)
      // — refuse rather than silently mis-sample (PLAN §Phase 2b Edge Cases).
      MFEM_VERIFY(material.mode == MaterialField::Mode::Coefficient,
                  "spatial_seas_driver: heterogeneous material must be "
                  "MaterialField::Mode::Coefficient (per-qp); Mode::GridFunction "
                  "is not supported by the quasi-dynamic IP operator. Use the "
                  "depth-profile or velocity-sidecar Coefficient path.");
      MFEM_VERIFY(material.lambda_coef != nullptr && material.mu_coef != nullptr,
                  "spatial_seas_driver: Mode::Coefficient MaterialField has a "
                  "null lambda_coef/mu_coef handle.");
      domain_ptr = std::make_unique<ElasticityDomainOperator<ParMesh>>(
         pmesh, cfg.mesh.order, *material.lambda_coef, *material.mu_coef,
         plate_rate, fault_depth, fault_length,
         bc, DGMethod::IP, spatial::ParseQDSolverType(cfg.solver.type),
         domain_config);
   }
   ElasticityDomainOperator<ParMesh> &domain = *domain_ptr;

   if (rank == 0)
   {
      std::cout << "[spatial_seas] elasticity domain operator constructed "
                << "(rank 0 local counts):\n"
                << "  GetNumFaultDOFs:       " << domain.GetNumFaultDOFs() << "\n"
                << "  GetNumOwnedFaultDOFs:  " << domain.GetNumOwnedFaultDOFs() << "\n"
                << "  GetNumFaultFaces:      " << domain.GetNumFaultFaces() << "\n"
                << "  GetNbfPerFace:         " << domain.GetNbfPerFace() << "\n";
   }

   // -----------------------------------------------------------------
   // 4b. Fault geometry + spatial per-DOF rate-state (Phase 3).  NORMATIVE
   //     ORDER (R-001/R-002/R-003):  geom(compute_bp5_params=false) →
   //     owned-order resolver tables → ResolveRateState (TOTAL σ_n) →
   //     SetRateStatePerDOF.  The stress source (geom.ComputeParams*) and the
   //     RateStateFaultOperator ctor are DEFERRED to Phase 5.  Only the
   //     rate-state law is wired (the QD driver is rate-state).
   // -----------------------------------------------------------------
   int rc = 0;   // process exit code (set non-zero only on the no-rate-state error path)
   if (cfg.rate_state.has_value())
   {
      MPIContext mpi(comm);                              // owns_comm=false
      BP5Params seed;                                    // seeds c_s/μ → η scaling only
      FaultGeometry<ParMesh> geom(domain, seed, &mpi,
                                  /*compute_bp5_params=*/false);

      // (1) OWNED-order per-DOF tables.  GetFaultDOFToElem/Attr are owned;
      //     GetFaultDOFCoords3D is LOCAL and MUST be RestrictToOwnedFault'd
      //     before pairing with them (R-002).
      Array<int> dof_to_elem, dof_to_attr;
      domain.GetFaultDOFToElem(dof_to_elem);
      domain.GetFaultDOFToAttr(dof_to_attr);
      Vector local_coords_3d, dof_coords_3d;
      domain.GetFaultDOFCoords3D(local_coords_3d);
      domain.RestrictToOwnedFault(local_coords_3d, dof_coords_3d, /*comps=*/3);

      const int n_owned = domain.GetNumOwnedFaultDOFs();
      MFEM_VERIFY(dof_to_elem.Size() == n_owned
                  && dof_to_attr.Size() == n_owned
                  && dof_coords_3d.Size() == 3 * n_owned,
                  "spatial_seas: Phase 3 owned-table size mismatch (elem="
                  << dof_to_elem.Size() << ", attr=" << dof_to_attr.Size()
                  << ", coords/3=" << dof_coords_3d.Size() / 3
                  << ", owned=" << n_owned << ")");

      // (2) TOTAL effective-normal-stress per owned DOF — sourced
      //     INDEPENDENTLY of geom (the stress source runs in Phase 5).  Empty
      //     for uniform σ_n (BP5 parity uses [friction.rate_state].
      //     sigma_n_default, which is the EFFECTIVE σ_n).
      Vector sigma_n_total_owned;   // size 0 ⇒ resolver uses sigma_n_default

      // (3) Resolve per-DOF rate-state.  LAST arg is TOTAL σ_n (the resolver
      //     forms sigma_n_eff = sigma_n_total − P_p).  With an empty
      //     sigma_n_total the resolver uses sigma_n_default (already
      //     EFFECTIVE), so a ZERO PorePressureSpec is passed to avoid a
      //     second P_p subtraction (R-003; matches spatial_dyn_driver).
      spatial::RateStatePerDOFParams rs;
      if (cfg.rate_state->bp5_analytic)
      {
         // Phase 7 (R-007): BP5-native rate-state.  Fill rs DIRECTLY from
         // bp5_params at the owned fault coords, EXACTLY mirroring
         // FaultGeometry::ComputeBP5Params (a_of_x2_x3, Dc_of_x2_x3, V_init_vec,
         // eta()), so the rate-state a and the bp5_analytic prestress are
         // consistent and the run byte-matches seas_driver.cpp.  The spatial-rule
         // resolver cannot express BP5's a(x2,x3) smooth transition (Box rules
         // are hard; boxcar_taper is rejected), so this bypasses it.
         const Vector &cx2 = geom.GetCoordsX2();
         const Vector &cx3 = geom.GetCoordsX3();
         const int Nrs = cx2.Size();
         MFEM_VERIFY(cx3.Size() == Nrs,
                     "spatial_seas: fault coord size mismatch (x2=" << Nrs
                     << ", x3=" << cx3.Size() << ").");
         rs.a.SetSize(Nrs);      rs.b.SetSize(Nrs);       rs.Dc.SetSize(Nrs);
         rs.eta.SetSize(Nrs);    rs.V_init.SetSize(Nrs);  rs.sigma_n_eff.SetSize(Nrs);
         rs.f_0.SetSize(Nrs);    rs.V_0.SetSize(Nrs);
         const real_t eta_bp5 = seed.eta();
         for (int i = 0; i < Nrs; ++i)
         {
            rs.a(i)           = seed.a_of_x2_x3(cx2(i), cx3(i));
            rs.Dc(i)          = seed.Dc_of_x2_x3(cx2(i), cx3(i));
            rs.b(i)           = seed.b;
            rs.eta(i)         = eta_bp5;
            real_t Vi[2];
            seed.V_init_vec(cx2(i), cx3(i), Vi);
            // R-703: the magnitude + init_vel_dir=(0,1) reconstruction below
            // (in SetRateStatePerDOF) only reproduces bp5_params.V_init_vec when
            // V_init is pure +strike (true for BP5).  Refuse silently dropping a
            // SIGNIFICANT dip/negative-strike component.  The dip slot carries the
            // V_zero placeholder (=1e-20, ~11 orders below V_init=1e-9), which the
            // native BP5 path also keeps; use a RELATIVE tolerance (dip negligible
            // vs strike) so the reconstruction error stays ~1e-12, not an absolute
            // 1e-30 that erroneously rejects the legitimate placeholder.
            MFEM_VERIFY(std::abs(Vi[0]) <= 1e-6 * std::abs(Vi[1]) && Vi[1] >= 0.0,
                        "spatial_seas: bp5_analytic rate-state assumes pure "
                        "+strike V_init (dip << strike); got Vi=(" << Vi[0] << ", "
                        << Vi[1] << ") at fault DOF " << i << ".");
            rs.V_init(i)      = std::sqrt(Vi[0] * Vi[0] + Vi[1] * Vi[1]);
            rs.sigma_n_eff(i) = seed.sigma_n;
            rs.f_0(i)         = seed.f0;
            rs.V_0(i)         = seed.V0;
         }
      }
      else
      {
         spatial::SpatialFrictionResolver resolver;
         spatial::PorePressureSpec pp;  // default-empty for the uniform-σ_n path
         rs = resolver.ResolveRateState(
            *cfg.rate_state, dof_coords_3d, dof_to_elem, dof_to_attr,
            material, pmesh, pp, sigma_n_total_owned);
      }

      // Edge case: the QD fault operator bakes b/f0/V0 as SCALARS
      // (DieterichRuinaFriction/AgingLawPsi); refuse per-DOF variation rather
      // than silently averaging (see Phase 3b).
      //
      // R-301: the uniformity check must be MPI-GLOBAL, not per-rank.  A field
      // can be uniform on each rank yet vary globally (e.g. a region-partitioned
      // b whose region boundary aligns with the mesh partition), so reduce the
      // local extrema to global BEFORE comparing.  assert_uniform is called
      // unconditionally below on EVERY rank inside this collective rate_state
      // block, so MPIContext::GlobalMin/Max (each an MPI_Allreduce) are reached
      // by all ranks — do NOT add a per-rank early return (it would desync the
      // collective and deadlock).  Empty local view ⇒ contributes the identity
      // (max for min, lowest for max); if no rank has any fault DOF, hi < lo and
      // the check is skipped.
      auto assert_uniform = [&](const Vector &v, const char *name)
      {
         real_t lo = std::numeric_limits<real_t>::max();
         real_t hi = std::numeric_limits<real_t>::lowest();
         for (int i = 0; i < v.Size(); ++i)
         {
            lo = std::min(lo, v(i));
            hi = std::max(hi, v(i));
         }
         lo = mpi.GlobalMin(lo);
         hi = mpi.GlobalMax(hi);
         if (hi < lo) { return; }   // no fault DOFs on any rank
         const real_t mean = 0.5 * (lo + hi);
         MFEM_VERIFY(hi - lo <= 1e-12 * (std::abs(mean) + 1.0),
                     "spatial_seas: per-DOF " << name << " is not uniform "
                     "(global range " << (hi - lo) << ") — per-DOF b/f0/V0 is "
                     "not yet supported by the QD fault operator (see Phase 3b).");
      };
      assert_uniform(rs.b,   "b");
      assert_uniform(rs.f_0, "f_0");
      assert_uniform(rs.V_0, "V_0");

      // Fallback guard for SAFS-mode (degenerate fault normals would corrupt
      // the per-DOF τ_pre/σ_n mapping).
      MFEM_VERIFY(geom.NumZeroNormalFallbacks() == 0
                  && geom.NumT1Fallbacks() == 0,
                  "spatial_seas: FaultGeometry reported degenerate fault "
                  "normals (zero-normal=" << geom.NumZeroNormalFallbacks()
                  << ", t1-fallback=" << geom.NumT1Fallbacks()
                  << ") — refusing to map per-DOF rate-state.");

      // (4) Fill geom's per-DOF a/Dc/eta/V_init.  Runs BEFORE the fault-op
      //     ctor and MUST NOT set params_computed_ (R-001).  init_vel_dir is
      //     the prescribed loading direction (dip,strike); BP5 is pure strike.
      Vector init_vel_dir(2); init_vel_dir(0) = 0.0; init_vel_dir(1) = 1.0;
      // R-S1: the resolver fills rs.V_w to size N for EVERY rate-state config,
      // but the QD operator reads per-DOF V_w only for slip_law_srw.  Clear it
      // for non-SRW so geom's V_w stays EMPTY — preserving the operator's
      // missing-V_w sentinel + the "aging leaves V_w empty" invariant.  The
      // bp5_analytic path never sets rs.V_w, so this is a no-op there.
      if (!(cfg.rate_state.has_value() &&
            cfg.rate_state->state_evolution ==
               spatial::StateEvolutionKind::SlipLawStrongRateWeakening))
      { rs.V_w.SetSize(0); }
      geom.SetRateStatePerDOF(rs, init_vel_dir);

      if (rank == 0)
      {
         std::cout << "[spatial_seas] Phase 3: per-DOF rate-state resolved + "
                   << "applied (owned fault DOFs this rank: " << n_owned
                   << "); geom.HasParams()="
                   << (geom.HasParams() ? "true" : "false")
                   << " (must be false until ComputeParams* in Phase 5).\n";
      }

      // ---------------------------------------------------------------
      // 5.  Dry-run: stop after Phase-2/3 construction, BEFORE the stress
      //     source / fault-op / SetInitialCondition solve / loop.  (The
      //     single MPI_Finalize is at function scope, AFTER geom/mpi/fault_op
      //     are destroyed — calling it here would free MPI comms in their
      //     dtors after finalize.)
      // ---------------------------------------------------------------
      if (dry_run)
      {
         if (rank == 0)
         {
            std::cout << "[spatial_seas] --dry-run: elasticity operator + "
                      << "per-DOF rate-state constructed, exiting before the "
                      << "time loop.\n";
         }
      }
      else
      {
      // ===============================================================
      // 5b. Phase 5 — SEAS coupling + adaptive RK45.  Construction order
      //     is NORMATIVE (R-001): geom(false) → SetRateStatePerDOF (above)
      //     → fault-op ctor → ComputeParams* (stress source) → SetSAFSMode.
      // ===============================================================
      // Scalar friction constants from the Phase-3 asserted-uniform per-DOF
      // vectors (fields are V_0/f_0, Vectors → element 0).  Per-DOF Dc still
      // flows through geom's Dc_values_; fc.Dc is the scalar fallback (R-005).
      // A rank that owns ZERO fault DOFs (normal at high rank counts — the planar
      // fault does not reach every volume partition; cf. seas_driver.cpp:431-435,
      // which always sources fc from config and so runs at 8N×400r) has empty rs.*
      // here.  On a fault-owning rank read the (assert_uniform-checked, globally
      // uniform) per-DOF element 0 exactly as before; on a fault-free rank fall
      // back to the config uniform defaults (== rs.*(0) for every shipped config).
      // Per-DOF Dc still flows through geom's Dc_values_; fc.Dc is the scalar
      // fallback (R-005).
      DieterichRuinaFriction::Constants fc;
      if (rs.b.Size() > 0)
      {
         fc.V0 = rs.V_0(0); fc.f0 = rs.f_0(0); fc.b = rs.b(0); fc.Dc = rs.Dc(0);
      }
      else
      {
         MFEM_VERIFY(cfg.rate_state.has_value(),
                     "spatial_seas: [friction.rate_state] required to build the "
                     "Dieterich-Ruina friction constants on a fault-free rank.");
         fc.V0 = cfg.rate_state->V_0_default;
         fc.f0 = cfg.rate_state->f_0_default;
         fc.b  = cfg.rate_state->b_default;
         fc.Dc = cfg.rate_state->Dc_default;
      }
      DieterichRuinaFriction friction(fc);

      // State-evolution law: aging-law (default) or STRONG rate weakening
      // (TPV104 SlipLawSRWPsi).  Both objects live in this scope so they
      // outlive fault_op; we pass the selected one.  For SRW we flip on
      // production mode (the base virtuals throw, so the operator MUST use the
      // per-DOF _SRW route — which it does via StateRate_/StateSteady_, reading
      // per-DOF V_w/a from geom).  The scalar a here is only the law's fallback
      // (per-DOF a reaches it through geom); b/f0/V0 are the Phase-3 asserted-
      // uniform values, so fc.b/V0/f0 are exact.
      const bool use_srw =
         cfg.rate_state.has_value() &&
         cfg.rate_state->state_evolution ==
            spatial::StateEvolutionKind::SlipLawStrongRateWeakening;
      // SRW weakening params, read only on the SRW path (use_srw ⟹ has_value()).
      // The non-SRW fallback values are inert: `srw` is constructed unconditionally
      // (it must outlive fault_op) but is never selected unless use_srw.
      const real_t srw_fw = use_srw ? cfg.rate_state->f_w_default : 0.1;
      const real_t srw_vw = use_srw ? cfg.rate_state->V_w_default : 0.1;

      AgingLawPsi   aging(fc.b, fc.V0, fc.f0);
      // rs.a(0): srw is constructed UNCONDITIONALLY (it must outlive fault_op even
      // when inert), so guard the element-0 read — a fault-free rank has empty rs.a.
      // srw only ever DRIVES DOFs when use_srw, and then over 0 DOFs on such a rank,
      // so any finite a is inert; use the config a_default fallback.
      const real_t srw_a = (rs.a.Size() > 0)
                           ? rs.a(0)
                           : (cfg.rate_state.has_value()
                              ? cfg.rate_state->a_default : real_t(0.010));
      SlipLawSRWPsi srw(srw_a, fc.b, fc.V0, fc.f0, srw_fw, srw_vw);
      StateEvolution *evo = &aging;
      if (use_srw)
      {
         srw.SetProductionMode();   // base virtuals throw -> per-DOF _SRW only
         evo = &srw;
         if (rank == 0)
         {
            std::cout << "[spatial_seas] state law: SlipLawSRWPsi (strong rate "
                         "weakening; f_w=" << srw_fw
                      << ", per-DOF V_w/a from geom)\n";
         }
      }

      // (A) fault operator BEFORE the stress source (R-001: ctor asserts
      //     !geom.HasParams() and caches Dc/V_init from the Phase-3 fill).
      RateStateFaultOperator<ParMesh, 2> fault_op(&geom, &friction, evo,
                                                  seed, &mpi);

      // (B) stress source → geom's per-DOF tau_pre_/sigma_n_per_dof_ (+ sets
      //     params_computed_).  The σ_n here MUST match the TOTAL σ_n that fed
      //     ResolveRateState in Phase 3 (R-003): for BP5 parity set
      //     [friction.rate_state].sigma_n_default to BP5's effective σ_n; for
      //     fault_local_prestress set [stress].sigma_n_pa to the same value.
      if (cfg.stress.bp5_analytic)
      {
         // BP5-parity per-DOF prestress (R-002/R-007): embed bp5_params.tau0_vec
         // + effective σ_n into a Cauchy tensor on the constant planar-BP5 fault
         // basis (col 0 of geom.fault_dof_basis(); ± a global sign flip is
         // projection-invariant — verified in test_spatial_seas_bp5_analytic).
         // R-606: bp5_analytic takes precedence over [stress].kind — warn once
         // if a non-default kind was also set (it is silently ignored).
         if (rank == 0 && cfg.stress.kind != spatial::StressSourceKind::ConstantTensor)
         {
            std::cerr << "[spatial_seas] WARNING: [stress].bp5_analytic=true "
                      << "overrides [stress].kind (the kind is ignored).\n";
         }
         // Col-0 fault-basis axes for the single planar-BP5 Cauchy tensor.  A rank
         // owning ZERO fault DOFs (normal at high rank counts — the planar fault
         // does not reach every volume partition) has a 9x0 basis: use dummy
         // orthonormal axes (UNUSED — geom.ComputeParams early-returns when
         // num_fault_dofs_==0) so the source object still constructs.  A
         // fault-owning rank reads col 0 + runs the planar check exactly as before.
         real_t nrm[3] = { 1.0, 0.0, 0.0 };
         real_t t1v[3] = { 0.0, 1.0, 0.0 };
         real_t t2v[3] = { 0.0, 0.0, 1.0 };
         if (geom.NumFaultDOFs() > 0)
         {
            const DenseMatrix &Bsis = geom.fault_dof_basis();
            MFEM_VERIFY(Bsis.Height() == 9 && Bsis.Width() == geom.NumFaultDOFs(),
                        "spatial_seas: unexpected fault_dof_basis shape.");
            nrm[0] = Bsis(0,0); nrm[1] = Bsis(1,0); nrm[2] = Bsis(2,0);
            t1v[0] = Bsis(3,0); t1v[1] = Bsis(4,0); t1v[2] = Bsis(5,0);
            t2v[0] = Bsis(6,0); t2v[1] = Bsis(7,0); t2v[2] = Bsis(8,0);
            // R-604: bp5_analytic applies this single col-0 basis to ALL DOFs, so
            // it is valid only for a PLANAR fault.  Assert the columns share the
            // normal direction (up to sign) instead of silently mis-projecting on
            // a non-planar fault.
            const int ncheck = std::min(geom.NumFaultDOFs(), 16);
            for (int j = 1; j < ncheck; ++j)
            {
               const real_t dot = nrm[0]*Bsis(0,j) + nrm[1]*Bsis(1,j) + nrm[2]*Bsis(2,j);
               MFEM_VERIFY(std::abs(std::abs(dot) - 1.0) < 1e-9,
                           "spatial_seas: [stress].bp5_analytic requires a PLANAR "
                           "fault — column " << j << " normal differs from column 0 "
                           "(|n_j·n_0|=" << std::abs(dot) << ").");
            }
         }
         spatial::Bp5AnalyticStressSource bp5src(seed, nrm, t1v, t2v);
         geom.ComputeParams(bp5src);   // P_p = 0 (BP5 σ_n is effective)
      }
      else if (cfg.stress.kind == spatial::StressSourceKind::FaultLocalPrestress)
      {
         // R-003 / R-603: the uniform fault-local σ_n MUST match the TOTAL σ_n
         // the resolver used in Phase 3 (rs.sigma_n_eff), or the 4-phase init
         // can't equilibrate.  Guard here with a clear message instead of a
         // confusing downstream SetInitialCondition equilibrium failure.
         real_t sn_err = 0.0;
         for (int i = 0; i < rs.sigma_n_eff.Size(); ++i)
         {
            sn_err = std::max(sn_err,
                              std::abs(rs.sigma_n_eff(i) - cfg.stress.sigma_n_pa));
         }
         MFEM_VERIFY(sn_err <= 1e-6 * (std::abs(cfg.stress.sigma_n_pa) + 1.0),
                     "spatial_seas: [stress].sigma_n_pa (" << cfg.stress.sigma_n_pa
                     << ") must equal the rate-state effective σ_n the resolver "
                     "used (max abs diff " << sn_err << "; set "
                     "[friction.rate_state].sigma_n_default to the same value "
                     "— plan R-003).");
         geom.ComputeParamsFaultLocal(cfg.stress.tau_strike_pa,
                                      cfg.stress.tau_dip_pa,
                                      cfg.stress.sigma_n_pa, /*P_p=*/0.0);
      }
      else if (cfg.stress.kind == spatial::StressSourceKind::ConstantTensor)
      {
         // SAFS no-flip projection sign convention (matches spatial_dyn_driver
         // + test_constant_tensor_sign): NEGATE the xy off-diagonal so a
         // right-lateral-POSITIVE sigma_xy INPUT maps to a right-lateral
         // on-fault tau_strike.  Pass the config pore pressure (the tensor σ_n
         // is TOTAL; ComputeParams subtracts P_p for the effective σ_n).
         spatial::ConstantTensorStressSource src(
            cfg.stress.sigma_xx_pa, cfg.stress.sigma_yy_pa, cfg.stress.sigma_zz_pa,
            -cfg.stress.sigma_xy_pa, cfg.stress.sigma_yz_pa, cfg.stress.sigma_xz_pa);
         geom.ComputeParams(src,
                            cfg.stress.pore_pressure.P_p_pa,
                            cfg.stress.pore_pressure.P_p_grad_pa_per_m,
                            cfg.stress.pore_pressure.min_sigma_n_pa);
      }
      else if (cfg.stress.kind ==
               spatial::StressSourceKind::DepthProportionalToShearModulus)
      {
         // Phase 8 (SAF): depth-proportional pre-stress
         //   sigma(x,y,z) = sigma_*_per_mu * mu(x,y,z) / mu_ref.
         // mu(x,y,z) comes from the current material (Constant mu_const or
         // depth_profile_1d via eval_at_xyz; sidecar has no coordinate-only mu).
         // Mirrors spatial_dyn_driver including the -sigma_xy SAFS sign flip.
         spatial::DepthProportionalToShearModulusStressSource::MuAtFn mu_at_xyz;
         if (material.mode == MaterialField::Mode::Constant)
         {
            const real_t mu_const = material.mu_const;
            mu_at_xyz = [mu_const](real_t, real_t, real_t) { return mu_const; };
         }
         else if (depth_profile_wrapper != nullptr)
         {
            MFEM_VERIFY(static_cast<bool>(depth_profile_wrapper->eval_at_xyz),
                        "spatial_seas: depth-profile material has no eval_at_xyz "
                        "callback (heterogeneous_material wiring lost?).");
            auto &eval = depth_profile_wrapper->eval_at_xyz;
            mu_at_xyz = [&eval](real_t x, real_t y, real_t z) -> real_t
            {
               real_t lam, mu, rho;
               eval(x, y, z, lam, mu, rho);
               return mu;
            };
         }
         else
         {
            MFEM_ABORT("spatial_seas: [stress].kind=\"depth_proportional\" "
                       "requires [material] kind=\"constant\" or "
                       "\"depth_profile_1d\" (no coordinate-only mu lookup for "
                       "sidecar_hdf5).");
         }
         const auto &dp = cfg.stress.depth_proportional;
         spatial::DepthProportionalToShearModulusStressSource src(
            dp.sigma_xx_per_mu, dp.sigma_yy_per_mu, dp.sigma_zz_per_mu,
            -dp.sigma_xy_per_mu, dp.sigma_yz_per_mu, dp.sigma_xz_per_mu,
            dp.mu_ref_pa, std::move(mu_at_xyz));
         geom.ComputeParams(src,
                            cfg.stress.pore_pressure.P_p_pa,
                            cfg.stress.pore_pressure.P_p_grad_pa_per_m,
                            cfg.stress.pore_pressure.min_sigma_n_pa);
      }
      else   // SidecarHDF5: per-DOF tau0/σ_n from a precomputed CSM HDF5 sidecar.
      {
         spatial::ApplyCsmStressSidecar(cfg.stress, geom);
      }

      // (C) point the fault op at geom's per-DOF τ_pre / effective σ_n.
      fault_op.SetSAFSMode(true, &geom.GetTauPre(), &geom.sigma_n_per_dof());

      PBP5SEASOp seas_op(&domain, &fault_op, &mpi);
      seas_op.SetElasticSigmaN(true);

      // Either restart from a checkpoint or run the 4-phase init.
      Vector state(fault_op.StateSize());
      real_t t = 0.0;
      int step = 0;
      real_t dt_seed = 0.0;
      const std::string ckpt_prefix = cfg.output.output_dir + "/spatial_seas";
      if (!restart_prefix.empty())
      {
         // Restart: restore the evolved fault state + (t, dt, step) and SKIP the
         // 4-phase init (the checkpointed state is already past equilibrium).
         SpatialSeasCheckpoint ck;
         const bool ok = ReadSpatialSeasCheckpoint(restart_prefix, rank, ck,
                                                   "spatial_seas");
         MFEM_VERIFY(ok, "spatial_seas: failed to read checkpoint at prefix '"
                     << restart_prefix << "' (rank " << rank
                     << ") — missing file, magic, or driver_tag mismatch.");
         MFEM_VERIFY(ck.state.Size() == state.Size(),
                     "spatial_seas: checkpoint state size " << ck.state.Size()
                     << " != fault StateSize " << state.Size()
                     << " (mesh / partition / rank-count mismatch on restart).");
         // R-602: the static params were re-derived above (SetRateStatePerDOF +
         // ComputeParams* ran deterministically from config).  Validate them
         // against the checkpoint so a CHANGED config between the original run
         // and the restart aborts instead of silently continuing the old state
         // under new friction params.
         auto check_param = [&](const char *nm, const Vector &re, const Vector &ckv)
         {
            MFEM_VERIFY(re.Size() == ckv.Size(),
                        "spatial_seas restart: " << nm << " size "
                        << re.Size() << " != checkpoint " << ckv.Size());
            real_t e = 0.0;
            for (int i = 0; i < re.Size(); ++i)
            { e = std::max(e, std::abs(re(i) - ckv(i))); }
            MFEM_VERIFY(e <= 1e-9 * (re.Normlinf() + 1.0),
                        "spatial_seas restart: " << nm << " differs from the "
                        "checkpoint (config drift, max abs diff " << e << ") — "
                        "the restart config must match the original run.");
         };
         check_param("a",       geom.GetAValues(),      ck.a);
         check_param("Dc",      geom.GetDcValues(),      ck.Dc);
         check_param("eta",     geom.GetEtaValues(),     ck.eta);
         check_param("V_init",  geom.GetVInit(),         ck.V_init);
         check_param("sigma_n", geom.sigma_n_per_dof(),  ck.sigma_n);
         check_param("tau_pre", geom.GetTauPre(),        ck.tau_pre);
         state   = ck.state;
         t       = ck.t;
         step    = ck.step;
         dt_seed = ck.dt;
         if (rank == 0)
         {
            std::cout << "[spatial_seas] restart from '" << restart_prefix
                      << "': t=" << t << ", step=" << step << ", dt=" << dt_seed
                      << " (skipping SetInitialCondition).\n";
         }
      }
      else
      {
         // 4-phase init; aborts internally if the equilibrium residual > 1e-6.
         seas_op.SetInitialCondition(state);
         // rs.V_init(0): a fault-free rank has empty rs.V_init — fall back to the
         // config V_init_default.  dt_seed is only the initial RK45 step guess
         // (adapted immediately, then reconciled collectively across ranks).
         const real_t v_init_seed = (rs.V_init.Size() > 0)
                                    ? rs.V_init(0)
                                    : (cfg.rate_state.has_value()
                                       ? cfg.rate_state->V_init_default : real_t(1e-9));
         dt_seed = (cfg.time.dt_init > 0.0)
                   ? cfg.time.dt_init
                   : 0.01 * fc.Dc / std::max(v_init_seed, real_t(1e-15));
         // GetMaxSlipRate() is COLLECTIVE (MPI_Allreduce over the global V_max):
         // call it on ALL ranks, then print on rank 0.  Calling it inside
         // `if (rank == 0)` deadlocks at np>1 (rank 0 blocks in the Allreduce
         // while the other ranks reach the MPI_Barrier below).
         const real_t init_vmax = seas_op.GetMaxSlipRate();
         if (rank == 0)
         {
            std::cout << "[spatial_seas] SetInitialCondition done: V_max = "
                      << init_vmax << " m/s\n";
         }
      }

      // Adaptive Dormand-Prince RK45 (mirrors seas_driver.cpp: SetDt seeds the
      // trial step, Step returns the accepted dt; the operator is Step's first
      // argument — R-004).
      DormandPrinceRK45 ode;
      // MUST initialise the integrator before the first Step() — Init allocates
      // the 7 RK stage vectors (k_[i].SetSize(op.Width())) that Step passes to
      // seas_op.Mult() as the output `rate`.  Without it the stage vectors are
      // size 0 (null data) and the first Mult -> ComputeRHS write segfaults.
      ode.Init(seas_op);
      ode.SetAbsTol(cfg.time.rk45_atol);
      ode.SetRelTol(cfg.time.rk45_rtol);
      ode.SetDtMax(cfg.time.dt_max_years * BP5Params::seconds_per_year);
      ode.SetMPIContext(&mpi);
      ode.SetDt(dt_seed);

      // Checkpoint writer: flat state + (t,dt,step) + geom per-DOF static params.
      const int ckpt_every = cfg.output.checkpoint_every_steps;
      auto write_checkpoint = [&](real_t t_now, real_t dt_now, int step_now)
      {
         SpatialSeasCheckpoint ck;
         ck.t = t_now; ck.dt = dt_now; ck.step = step_now;
         ck.n_owned = geom.GetAValues().Size();
         ck.state   = state;
         ck.a       = geom.GetAValues();
         ck.Dc      = geom.GetDcValues();
         ck.eta     = geom.GetEtaValues();
         ck.V_init  = geom.GetVInit();
         ck.sigma_n = geom.sigma_n_per_dof();
         ck.tau_pre = geom.GetTauPre();
         WriteSpatialSeasCheckpoint(ckpt_prefix, rank, ck);
      };

      // R-605: ensure the output directory exists (rank 0 creates it, then a
      // barrier) when checkpointing OR BP5 station output is on.
      if (ckpt_every > 0 || cfg.output.bp5_stations)
      {
         if (rank == 0)
         {
            std::error_code ec;
            std::filesystem::create_directories(cfg.output.output_dir, ec);
         }
#ifdef MFEM_USE_MPI
         MPI_Barrier(comm);
#endif
      }

      // BP5 SCEC station/probe output (Phase 6; gated by [output].bp5_stations,
      // BP5-parity runs only).  Mirrors seas_driver.cpp:378-662 so the traces
      // byte-match the golden BP5 driver (confirm parity on Frontera).  Loop
      // writes rely on Dormand-Prince FSAL (the accepted step's last stage IS
      // the accepted state, so GetTraction()/GetMaxSlipRate() are accepted-state
      // values, matching the golden driver); the initial write re-solves once to
      // refresh them (a restart restored `state` without a Mult).
      std::unique_ptr<ParallelBP5BenchmarkOutput> bench_out;
      if (cfg.output.bp5_stations)
      {
         const std::string full_prefix = cfg.output.output_dir + "/spatial_seas";
         auto stations = BP5BenchmarkOutput<ParMesh>::DefaultStations();
         Vector local_x2 = geom.GetCoordsX2();
         Vector local_x3 = geom.GetCoordsX3();
         const int N_local = geom.NumLocalFaultDOFs();
         Vector tp_dip(N_local), tp_strike(N_local);
         const Vector &tp = geom.GetTauPre();   // [dip, strike] per DOF
         for (int i = 0; i < N_local; ++i)
         {
            tp_dip(i)    = tp(2 * i);
            tp_strike(i) = tp(2 * i + 1);
         }
         if (rank == 0)
         {
            for (const auto &st : stations)
            { std::remove((full_prefix + "_" + st.name + ".txt").c_str()); }
            std::remove((full_prefix + "_global.txt").c_str());
         }
#ifdef MFEM_USE_MPI
         MPI_Barrier(comm);
#endif
         bench_out = std::make_unique<ParallelBP5BenchmarkOutput>(
            full_prefix, seed, stations, geom, mpi, local_x2, local_x3,
            tp_dip, tp_strike, domain.GetNbfPerFace());
         Vector rate_scratch(state.Size());
         seas_op.Mult(state, rate_scratch);   // refresh V/traction from `state`
         bench_out->ForceWrite(t, state, fault_op, seas_op.GetTraction(),
                               seas_op.GetMaxSlipRate());
         bench_out->Flush();
      }

      const int max_steps = GetIntArg(argc, argv, "--max-steps", -1);
      while (t < cfg.time.tfinal && (max_steps < 0 || step < max_steps))
      {
         if (t + ode.GetDt() > cfg.time.tfinal)
         {
            ode.SetDt(cfg.time.tfinal - t);
         }
         real_t dt_taken;
         if (!ode.Step(seas_op, state, t, dt_taken)) { continue; }
         ++step;
         const real_t V_max = seas_op.GetMaxSlipRate();
         MFEM_VERIFY(std::isfinite(V_max),
                     "spatial_seas: non-finite V_max at step " << step
                     << " (t=" << t << ") — aborting.");
         if (rank == 0 && (step <= 5 || step % 10 == 0))
         {
            std::cout << "[spatial_seas] step " << step << "  t=" << t
                      << "  dt=" << dt_taken << "  V_max=" << V_max << "\n";
         }
         // BP5 station/probe write at the adaptive output cadence (FSAL ⇒
         // GetTraction()/V_max are accepted-state values).
         if (bench_out
             && bench_out->Write(t, state, fault_op, seas_op.GetTraction(), V_max))
         {
            bench_out->Flush();
         }
         if (ckpt_every > 0 && step % ckpt_every == 0)
         {
            // R-601: persist the NEXT trial dt (ode.GetDt() after the accepted
            // step), NOT the just-taken dt_taken — a restart seeds the RK45 with
            // this, and the un-restarted run's next step uses ode.GetDt().  The
            // final checkpoint below already does this.
            write_checkpoint(t, ode.GetDt(), step);
         }
      }
      // Final checkpoint (so a restart can continue from the last step even if
      // it was not on a --checkpoint-every boundary).
      if (ckpt_every > 0 && step > 0) { write_checkpoint(t, ode.GetDt(), step); }
      // Final forced station write + close (flush the last accepted state).
      if (bench_out)
      {
         bench_out->ForceWrite(t, state, fault_op, seas_op.GetTraction(),
                               seas_op.GetMaxSlipRate());
         bench_out->Close();
      }
      if (rank == 0)
      {
         std::cout << "[spatial_seas] time loop done: " << step << " steps, t="
                   << t << " s (tfinal=" << cfg.time.tfinal << " s).\n";
      }
      }   // end else (Phase 5 run path)
   }      // end if (cfg.rate_state.has_value()) — geom/mpi/fault_op destruct here
   else
   {
      // No [friction.rate_state]: the quasi-dynamic driver is rate-state, so
      // there is nothing to time-step.  Construct-only (dry-run) succeeds; a
      // real run errors out (do NOT silently no-op — CLAUDE.md).
      if (rank == 0)
      {
         if (dry_run)
         {
            std::cout << "[spatial_seas] --dry-run: elasticity operator "
                      << "constructed (no [friction.rate_state] block).\n";
         }
         else
         {
            std::cerr << "spatial_seas_driver: the quasi-dynamic driver requires "
                      << "a [friction.rate_state] block (the QD law is "
                      << "rate-state).\n";
         }
      }
      if (!dry_run) { rc = 1; }
   }

   // Single MPI_Finalize at function scope, AFTER all block-local MPI-owning
   // objects (geom, MPIContext, fault_op, seas_op) have destructed.
#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return rc;
}
