// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV102 benchmark driver: 3D dynamic rupture on a vertical strike-slip fault.
//
// Fully parallel MPI driver using ParMesh and WaveOperator<ParMesh>.
// Connects WaveOperator (Phase 1) + FaultFaceFlux (Phase 3) + first-order ABC
// with RK4 time stepping and TPV102-specific initialization.
//
// Usage:
//   ibrun ./seas_tpv102_driver --mesh tpv102/mesh/tpv102_fine.msh \
//       --mesh-scale 1000 --order 2 --bc-mode absorbing \
//       --tfinal 12.0 --output-dir results/ --output-prefix tpv102
//
// Reference: SCEC TPV101/102 benchmark specification.

#include "mfem.hpp"
#include "../dynamic/wave_state.hpp"
#include "../dynamic/wave_operator.hpp"
#include "../dynamic/fault_face_flux.hpp"
#include "../dynamic/pml_layer.hpp"
#include "../dynamic/seas_dynamic_operator.hpp"
#include "../dynamic/tpv102_setup.hpp"
#include "../dynamic/tpv102_setup_total.hpp"
#include "../config/tpv102_params.hpp"
#include "../domain/boundary_config.hpp"
#include "../io/paraview_output.hpp"
#include "../dynamic/seas_diag_rank.hpp"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <cmath>
#include <vector>
#include <sys/stat.h>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

#ifdef SEAS_DIAG_FAULT_FLUX
// v9.0.0 §0.5.2 preamble: definition of the global rank cache declared in
// dynamic/seas_diag_rank.hpp.  Written once at MPI init inside main();
// read-only thereafter.  Defined here (driver TU) so it has exactly one
// definition across the whole link.
namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT
#endif

static std::string GetStringArg(int argc, char *argv[], const char *flag,
                                const std::string &default_val)
{
   for (int i = 1; i < argc - 1; i++)
   {
      if (std::string(argv[i]) == flag) { return argv[i+1]; }
   }
   return default_val;
}

static real_t GetRealArg(int argc, char *argv[], const char *flag,
                         real_t default_val)
{
   std::string val = GetStringArg(argc, argv, flag, "");
   if (val.empty()) { return default_val; }
   return std::stod(val);
}

static int GetIntArg(int argc, char *argv[], const char *flag, int default_val)
{
   std::string val = GetStringArg(argc, argv, flag, "");
   if (val.empty()) { return default_val; }
   return std::stoi(val);
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
   // v9.0.0 §0.5.2: seed the global rank cache used by C-1/C-2/C-3 DIAG.
   mfem::seas::g_seas_my_rank = rank;
#endif

   // R-501 + R-606: diagnostic-flag startup banner.  Print on stderr AND
   // write to build_info.txt (rank 0 only).  Every Phase 2+ diagnostic
   // analysis must verify the banner shows the expected flag state
   // before interpreting DIAG output: MFEM's config.mk can silently drop
   // -D flags passed via CXXFLAGS+=..., and Frontera stderr may be
   // interleaved across ranks.  The file is rank-0-only and deterministic.
   if (rank == 0)
   {
#ifdef SEAS_DIAG_FAULT_FLUX
      const char *diag_fault_flux = "SEAS_DIAG_FAULT_FLUX = ON";
#else
      const char *diag_fault_flux = "SEAS_DIAG_FAULT_FLUX = OFF";
#endif
#ifdef SEAS_DIAG_GHOST_EXCHANGE
      const char *diag_ghost = "SEAS_DIAG_GHOST_EXCHANGE = ON";
#else
      const char *diag_ghost = "SEAS_DIAG_GHOST_EXCHANGE = OFF";
#endif
      std::fprintf(stderr, "[BUILD] %s\n", diag_fault_flux);
      std::fprintf(stderr, "[BUILD] %s\n", diag_ghost);

      std::ofstream binfo("build_info.txt");
      if (binfo.is_open())
      {
         binfo << "[BUILD] " << diag_fault_flux << "\n";
         binfo << "[BUILD] " << diag_ghost << "\n";
         binfo.close();
      }
      else
      {
         std::fprintf(stderr, "[WARNING] could not open build_info.txt for "
                              "banner (Phase 2 should fall back to stderr)\n");
      }
   }

   // -----------------------------------------------------------------------
   // Parse command-line arguments
   // -----------------------------------------------------------------------
   std::string mesh_file = GetStringArg(argc, argv, "--mesh",
                                        "tpv102/mesh/tpv102_coarse.msh");
   real_t mesh_scale = GetRealArg(argc, argv, "--mesh-scale", 1.0);
   // Round-7 R-001: default DG polynomial order is p=1 so the default
   // `--ader-order=2` satisfies the ADER plan's `O >= p + 1`
   // consistent-order rule (plan §Numerical constraints).  Higher p
   // requires a matching bump in --ader-order; at --ader-order >= 3 the
   // nucleation time-integration is only 2nd-order-accurate under the
   // current mid-step scheme, so land the deferred plan Phase 7 §4
   // higher-order nucleation before flipping this back to p=2.
   int order = GetIntArg(argc, argv, "--order", 1);
   std::string bc_mode = GetStringArg(argc, argv, "--bc-mode", "absorbing");
   real_t tfinal = GetRealArg(argc, argv, "--tfinal", TPV102Params::t_final);
   std::string output_dir = GetStringArg(argc, argv, "--output-dir", "tpv102/results");
   std::string output_prefix = GetStringArg(argc, argv, "--output-prefix", "tpv102");
   real_t cfl_factor = GetRealArg(argc, argv, "--cfl", 0.5);
   int bc_free = GetIntArg(argc, argv, "--bc-free", 1);
   int bc_fault = GetIntArg(argc, argv, "--bc-fault", 3);
   int bc_absorb = GetIntArg(argc, argv, "--bc-absorb", 5);

   // ADER I-05 Phase 7 (round-6: purpose change #1): time integrator
   // selection.  ADER is now the default; RK4 is kept only as an
   // alternative for byte-compat regression runs.
   // --time-integrator {ader|rk4}  : default ader.
   // --ader-order N                : ADER order in {2, 3, 4}, default 2.
   //                                 Ignored when --time-integrator=rk4.
   std::string time_integrator =
      GetStringArg(argc, argv, "--time-integrator", "ader");
   int ader_order = GetIntArg(argc, argv, "--ader-order", 2);
   // Normalise + validate.
   for (auto &ch : time_integrator) { ch = std::tolower(ch); }
   // R-005 (round-6): `use_ader` must be MUTABLE so the unknown-flag
   // fallback branch below can re-sync it with the post-fallback value
   // of `time_integrator` — otherwise a typo like `--time-integrator=ader3`
   // would warn "falling back to 'ader'" yet silently take the RK4
   // branch.
   bool use_ader = (time_integrator == "ader");
   if (!use_ader && time_integrator != "rk4")
   {
      if (rank == 0)
      {
         std::cerr << "[WARNING] --time-integrator '" << time_integrator
                   << "' unrecognised — falling back to 'ader'.\n";
      }
      time_integrator = "ader";
      use_ader = true;
   }
   if (use_ader && (ader_order < 2 || ader_order > 4))
   {
      if (rank == 0)
      {
         std::cerr << "[WARNING] --ader-order " << ader_order
                   << " out of {2,3,4} — clamping to 2.\n";
      }
      ader_order = 2;
   }

   // ParaView output controls (mirrors BP5 --paraview* flags).
   // --paraview              : enable PVD/VTU output, interval matches --output-dt
   // --paraview-every N      : write every N steps
   // --paraview-dt X         : write every X seconds (overrides step interval)
   // --pv-low-order          : disable high-order output + levels-of-detail=1
   //                           (linear tets only — ~40x smaller volume output)
   // --no-domain-pv          : suppress volume-mesh PVD (ParaView/Cycle*/*.vtu);
   //                           fault-surface PVD + VTUs are still written on the
   //                           same schedule.  Saves tremendous disk on large runs.
   // --debug-qnorm           : print per-rank ||Q||_inf at every station output
   //                           cycle.  Diagnostic for tpv102_debug_v1.md H1 —
   //                           checks whether bulk wave energy crosses rank
   //                           partition seams.  Rank 0 prints:
   //                             global {min, max, mean} of ||Q||_inf across ranks
   //                             per-rank ||Q||_inf for a small sampled set
   //                           Cost: one MPI_Gather per output cycle, negligible.
   bool use_paraview = false;
   bool pv_low_order = false;
   bool pv_no_domain = false;
   bool debug_qnorm  = false;
   int  paraview_step_interval = 0;
   real_t paraview_dt_flag = 0.0;
   // Bulk (volume) ParaView cadence.  <= 0 disables the bulk collection
   // entirely (default).  Positive values enable a SECOND ParaView
   // collection written at the specified seconds interval (to
   // `ParaView_bulk/`), carrying velocity, sigma_yy, sigma_xy, sigma_xz,
   // mpi_rank.  Independent of --paraview-dt / --paraview-every (those
   // set the fault-surface schedule), and independent of --no-domain-pv
   // (which suppresses only the fault-schedule domain save on pv_out).
   real_t paraview_bulk_dt = 0.0;
   for (int i = 1; i < argc; i++)
   {
      std::string a = argv[i];
      if (a == "--paraview") { use_paraview = true; }
      else if (a == "--pv-low-order") { pv_low_order = true; }
      else if (a == "--no-domain-pv") { pv_no_domain = true; }
      else if (a == "--debug-qnorm") { debug_qnorm = true; }
      else if (a == "--paraview-every" && i + 1 < argc)
      {
         use_paraview = true;
         paraview_step_interval = std::atoi(argv[++i]);
      }
      else if (a == "--paraview-dt" && i + 1 < argc)
      {
         use_paraview = true;
         paraview_dt_flag = std::atof(argv[++i]);
      }
      else if (a == "--paraview-bulk-dt" && i + 1 < argc)
      {
         use_paraview = true;
         paraview_bulk_dt = std::atof(argv[++i]);
      }
   }

   if (rank == 0)
   {
      std::cout << "========================================\n";
      std::cout << "SCEC TPV102 Dynamic Rupture Simulation\n";
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
      std::cout << "========================================\n\n";

      mkdir(output_dir.c_str(), 0755);
   }

#ifdef MFEM_USE_MPI
   MPI_Barrier(comm);
#endif

   // -----------------------------------------------------------------------
   // 1. Load serial mesh on all ranks, partition to ParMesh
   // -----------------------------------------------------------------------
   Mesh serial_mesh(mesh_file.c_str(), 1, 1);
   MFEM_VERIFY(serial_mesh.Dimension() == 3, "TPV102 requires 3D mesh");

   if (mesh_scale != 1.0)
   {
      serial_mesh.SetCurvature(1, false, 3, Ordering::byVDIM);
      Vector &nodes = *serial_mesh.GetNodes();
      nodes *= mesh_scale;
   }

   // R-010: Validate domain size after scaling
   {
      Vector bbox_min, bbox_max;
      serial_mesh.GetBoundingBox(bbox_min, bbox_max);
      real_t domain_x = bbox_max(0) - bbox_min(0);
      if (rank == 0 && (domain_x < 1e3 || domain_x > 1e6))
      {
         std::cerr << "WARNING: Domain X-extent = " << domain_x
                   << " m (after scale=" << mesh_scale
                   << "). Expected ~60000-120000 m for TPV102. "
                   << "Check --mesh-scale.\n";
      }
   }

#ifdef MFEM_USE_MPI
   ParMesh pmesh(comm, serial_mesh);
   // Use ParMesh for the wave operator
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
      std::cout << "Mesh: " << ne_global << " elements total, "
                << nprocs << " ranks\n";
   }

   // -----------------------------------------------------------------------
   // 2. Boundary conditions (Tandem convention)
   // -----------------------------------------------------------------------
   BoundaryConfig bc;
   bc.natural_attrs = {bc_free};
   bc.fault_attr = bc_fault;
   bc.absorbing_attrs = {bc_absorb};

   // R-012: Validate boundary attributes exist in the mesh
   {
      int n_free = 0, n_fault = 0, n_absorb = 0;
      for (int b = 0; b < pmesh.GetNBE(); b++)
      {
         int attr = pmesh.GetBdrAttribute(b);
         if (bc.natural_attrs.count(attr)) { n_free++; }
         if (attr == bc.fault_attr) { n_fault++; }
         if (bc.absorbing_attrs.count(attr)) { n_absorb++; }
      }
      int n_free_g = n_free, n_fault_g = n_fault, n_absorb_g = n_absorb;
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&n_free, &n_free_g, 1, MPI_INT, MPI_SUM, comm);
      MPI_Allreduce(&n_fault, &n_fault_g, 1, MPI_INT, MPI_SUM, comm);
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
                  << " found in mesh. Check --bc-fault flag or mesh Physical Surface tags.");
      if (n_free_g == 0 && rank == 0)
      {
         std::cerr << "WARNING: No free-surface faces found (attr="
                   << bc_free << "). Check --bc-free flag.\n";
      }
   }

   // -----------------------------------------------------------------------
   // 3. Construct WaveOperator (parallel)
   // -----------------------------------------------------------------------
   WaveOperator<MeshT> wave(pmesh, order,
                            TPV102Params::lambda, TPV102Params::mu,
                            TPV102Params::rho, bc);

   int ndof_total = wave.GetScalarNDof();
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
   // 4. Time step: --dt overrides; otherwise derive from CFL/h_min/cp.
   //    (h_min already reduced across ranks inside the WaveOperator ctor.)
   // -----------------------------------------------------------------------
   real_t dt_override = GetRealArg(argc, argv, "--dt", 0.0);
   real_t cfl = cfl_factor / (3.0 * (2.0 * order + 1.0));
   real_t dt_cfl = wave.ComputeMaxDt(cfl);
   real_t dt = (dt_override > 0.0) ? dt_override : dt_cfl;

   int nsteps = static_cast<int>(std::ceil(tfinal / dt));
   if (rank == 0)
   {
      std::cout << "CFL: " << cfl << ", dt_cfl = " << dt_cfl << " s\n";
      if (dt_override > 0.0)
      {
         std::cout << "dt (override, --dt): " << dt << " s";
         if (dt > dt_cfl)
         {
            std::cout << "  [WARNING: dt_override > dt_cfl; may be unstable]";
         }
         std::cout << "\n";
      }
      else
      {
         std::cout << "dt: " << dt << " s\n";
      }
      std::cout << "Steps: " << nsteps << "\n\n";
   }

   // -----------------------------------------------------------------------
   // 5. Fault DOF data (local partition)
   // -----------------------------------------------------------------------
   L2_FECollection fec(order, 3, BasisType::GaussLobatto);
#ifdef MFEM_USE_MPI
   ParFiniteElementSpace fes(&pmesh, &fec);
#else
   FiniteElementSpace fes(&pmesh, &fec);
#endif

   // Fault-face geometry lists are owned by the wave operator (single source
   // of truth, BP5 pattern).  We iterate them in order — local-interior
   // faces first, then shared — to build fault_coords, matching the DOFData
   // layout that SetFaultDOFData imposes.
   const Array<int> &fault_int_faces   = wave.GetFaultInteriorFaces();
   const Array<int> &fault_shr_faces   = wave.GetFaultSharedFaces();

   // nqp_per_face is derived once from the first local fault face (all
   // fault faces share the same geometry type, so the quadrature rule
   // has the same point count).
   int nqp_per_face = 0;
   if (fault_int_faces.Size() > 0)
   {
      FaceElementTransformations *ftr0 =
         pmesh.GetInteriorFaceTransformations(fault_int_faces[0]);
      MFEM_VERIFY(ftr0,
                  "wave.GetFaultInteriorFaces() returned a face without "
                  "interior transformation — wave operator invariant violated");
      nqp_per_face = IntRules.Get(ftr0->GetGeometryType(), 2*order).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   else if (fault_shr_faces.Size() > 0)
   {
      FaceElementTransformations *ftr0 =
         pmesh.GetSharedFaceTransformations(fault_shr_faces[0]);
      MFEM_VERIFY(ftr0, "shared fault face has null transformation");
      nqp_per_face = IntRules.Get(ftr0->GetGeometryType(), 2*order).GetNPoints();
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
         ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   };

   for (int i = 0; i < fault_int_faces.Size(); i++)
   {
      push_qps(pmesh.GetInteriorFaceTransformations(fault_int_faces[i]));
   }
#ifdef MFEM_USE_MPI
   for (int i = 0; i < fault_shr_faces.Size(); i++)
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

   std::vector<DOFData> dof_data;
   if (num_fault_total > 0)
   {
      InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
      // Round-6 purpose change #2 (total-Q only): zero the DOFData
      // pre-stress fields so the wave operator's EvaluateTotal /
      // EvaluateADERTotal dispatch does not double-count pre-stress
      // (pre-stress lives in bulk Q via InitializeStateTotal below).
      ZeroDOFDataPreStressTotal(dof_data, num_fault_total);
   }

   FaultFaceFlux fault_flux(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   // Round-6 purpose change #2: supply the bulk background state for the
   // total-Q-aware BC variants (AbsorbingTotal / FreeSurface*Total / PML
   // damping toward Q_bg).  Q_bg mirrors the TPV102 pre-stress tensor
   // in global coordinates (same rotation as InitializeStateTotal).
   {
      real_t Q_bg[NUM_STATE] = {0};
      Q_bg[SYY] =  TPV102Params::sigma_n;
      Q_bg[SXY] = -TPV102Params::tau_ini;
      wave.SetAbsorbingBackground(Q_bg);
   }

   // Nucleation under total-Q now uses the persistent-prestress channel:
   // ApplyNucleationTotalPrestress overwrites DOFData.tau2_nuc per call;
   // FaultFaceFlux::EvaluateTotal adds it to the trial traction so the
   // requested dtau is re-imposed at every Riemann solve.  No bulk-Q
   // injection — see debug_document/tpv102_debug_document/
   //   tpv102_nucleation_code_review_2026-04-22.md  (the bug analysis)
   //   tpv102_nucleation_code_fix_2026-04-22.md     (this fix)
   // The legacy FaultQPNodalMap / FaultQPNucleationState /
   // BuildFaultQPNodalMap machinery in tpv102_setup_total.hpp is no
   // longer used by the production driver; it remains for the
   // R-002 / R-009 / R-I06-003 unit tests that verify per-call
   // arithmetic.  No fault-QP-to-nodal map needed here.

   // R-002 fix: resolve the rank that owns the hypocenter QP (closest local
   // fault QP to (hypo_along_strike, -hypo_down_dip) in x/z).  Used only by
   // --debug-qnorm to print a per-rank ||Q||_inf watch list.
   //
   // R-105 fix: use a named struct with static_asserts so the MPI_DOUBLE_INT
   // layout assumption fails loudly at compile time if it is ever broken
   // (e.g. by a compiler with unusual padding of {double, int}).
   int hypo_rank = 0;
   int hypo_dof_local = -1;  // local DOF index of the closest hypo QP, or -1
   {
      real_t local_min_dist2 = std::numeric_limits<real_t>::max();
      for (int i = 0; i < num_fault_total; i++)
      {
         real_t dx = fault_coords[i](0) - TPV102Params::hypo_along_strike;
         real_t dz = std::abs(fault_coords[i](2)) - TPV102Params::hypo_down_dip;
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
   // v9.0.0 §0.5: tag the hypocenter DOF on the owning rank so the C-1/C-2/
   // C-3 printf blocks emit a single line per RK4 Mult instead of flooding
   // stderr from every fault QP.  Only the rank that won the MINLOC above
   // sets diag_print=true; at most one DOF is flagged globally.
   if (rank == hypo_rank && hypo_dof_local >= 0 && num_fault_total > 0)
   {
      dof_data[hypo_dof_local].diag_print = true;
      const Vector &hpos = fault_coords[hypo_dof_local];
      std::fprintf(stderr,
         "[diag] rank %d tagging hypo DOF %d at (%.1f, %.1f, %.1f)\n",
         rank, hypo_dof_local, hpos(0), hpos(1), hpos(2));
   }
#endif

   // -----------------------------------------------------------------------
   // 6. Initialize state Q = TPV102 pre-stress tensor at every DOF
   //    (Round-6 purpose change #2: total-Q only — bulk Q carries the
   //    pre-stress, nucleation writes deltas into Q[SXY] at fault QPs,
   //    the fault face flux dispatches `EvaluateTotal` / `EvaluateADERTotal`).
   // -----------------------------------------------------------------------
   Vector Q;
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);

   // Round-6 purpose change #2: total-Q is the only supported driver
   // mode.  Fail loud if SetAbsorbingBackground was never called — the
   // wave operator's non-fault BC branches + fault dispatch fall
   // through to fluctuation semantics in that case and would silently
   // radiate the pre-stress / disable the fault.  This is the
   // driver-level hardening analogue of R-004's wave-operator-level
   // MFEM_VERIFY proposal: keeping the wave-operator else branches
   // alive (they support BP5 + the non-TPV102 unit tests) while
   // making sure the TPV102 driver cannot accidentally dispatch
   // through them.
   MFEM_VERIFY(wave.GetAbsorbingBackground() != nullptr,
               "tpv102_driver: total-Q migration requires "
               "wave.SetAbsorbingBackground() to have been called "
               "before any Mult / AdvanceADER; none detected.");

   if (rank == 0)
   {
      std::cout << "State representation: total-Q (pre-stress baked into Q)\n"
                << "Background: tau_strike = "
                << TPV102Params::tau_ini / 1e6
                << " MPa, sigma_n = "
                << TPV102Params::sigma_n / 1e6 << " MPa\n\n";
   }

   // -----------------------------------------------------------------------
   // 7. Station output (rank 0 only for fault stations)
   // -----------------------------------------------------------------------
   auto stations = DefaultStations();
   TPV102StationWriter station_writer;
   // R-004 fix: MPI ownership resolution — only the rank with the globally
   // nearest DOF opens and writes each station file.
#ifdef MFEM_USE_MPI
   station_writer.Open(output_dir, output_prefix, stations,
                       fault_coords, num_fault_local, comm);
#else
   station_writer.Open(output_dir, output_prefix, stations,
                       fault_coords, num_fault_local);
#endif
   station_writer.WriteStep(0.0, dof_data);

   // R-006 fix: Wire surface station writer into driver.
   // Uses FindPoints() for proper element containment + reference coords.
   auto surface_stations = DefaultSurfaceStations();
   TPV102SurfaceStationWriter surface_writer;
   surface_writer.Open(output_dir, output_prefix, surface_stations, pmesh, fes);
   surface_writer.WriteStep(0.0, Q);

   // -----------------------------------------------------------------------
   // 7b. ParaView output (mirrors BP5 seas::ParaViewOutput pattern):
   //   - Domain velocity (3-component L2 vector GridFunction built from Q)
   //   - MPI rank (L2 p=0 scalar, for partition visualization)
   //   - Fault L2-p0 fields: slip/slip_rate/traction (dip,strike), psi, sigma_n
   //   - Fault friction parameter fields: param_a, param_Dc, fault_x2, fault_x3
   //   - Fault-surface VTU/PVD (proper triangle geometry, per-DOF field values)
   //
   // Component-to-name mapping (BP5 convention, enforced project-wide by
   // R-801 Option A: comp 0 = dip = tangent1, comp 1 = strike = tangent2):
   //   DOFData.V1/slip1/tau1_corr (dip,    mode-III in TPV102) -> dip channel
   //   DOFData.V2/slip2/tau2_corr (strike, mode-II  in TPV102) -> strike channel
   // TPV102 is pure strike-slip so the dip channel stays near zero and the
   // interesting rupture physics lives in the strike channel.
   // -----------------------------------------------------------------------
   using PvFES = typename seas::GFType<MeshT>::FESType;
   using PvGF  = typename seas::GFType<MeshT>::type;

   std::unique_ptr<seas::ParaViewOutput<MeshT>> pv_out;
   std::unique_ptr<L2_FECollection> pv_vel_fec, pv_rank_fec;
   std::unique_ptr<PvFES> pv_vel_fes, pv_rank_fes;
   std::unique_ptr<PvGF>  pv_vel_gf,  pv_rank_gf;

   // Second ParaView collection for coarse bulk output at --paraview-bulk-dt.
   // Populated when paraview_bulk_dt > 0, independent of --no-domain-pv:
   // the bulk collection writes to `ParaView_bulk/` (separate directory),
   // so it does not conflict with --no-domain-pv which only suppresses the
   // fault-schedule pv_out volume save.  Carries velocity (shared GF with
   // pv_out — allocated regardless of pv_no_domain), three stress
   // components (sigma_yy, sigma_xy, sigma_xz — the ones most informative
   // for fault loading / mode-II radiation), and mpi_rank.
   std::unique_ptr<seas::ParaViewOutput<MeshT>> pv_bulk_out;
   std::unique_ptr<L2_FECollection> pv_bulk_sigma_fec;
   std::unique_ptr<PvFES> pv_bulk_sigma_fes;
   std::unique_ptr<PvGF>  pv_bulk_syy_gf, pv_bulk_sxy_gf, pv_bulk_sxz_gf;

   // Scratch vectors passed to UpdateFaultFieldsBP5 / WriteFaultSurfaceVTU.
   Vector pv_local_slip, pv_local_slip_rate, pv_local_traction;
   Vector pv_local_state, pv_local_normal_stress;
   Vector pv_local_a, pv_local_Dc, pv_local_x2, pv_local_x3;
   // R-V92-E02 diagnostic: stage-4 (pre-RK4-averaging) DOFData snapshot.
   // Populated from the stage-4 kN buffers and shipped to the fault-surface
   // VTU alongside the averaged fields.  Plan §19 uses the VTU delta
   // `normal_stress_k4 - normal_stress` to discriminate H-V92-K from H-V92-G.
   Vector pv_local_slip_rate_k4, pv_local_traction_k4, pv_local_normal_stress_k4;

   // output_interval = step-count matching --output-dt, reused by station
   // writers, console logging, and (as a default) ParaView output.
   real_t output_dt = GetRealArg(argc, argv, "--output-dt", 0.01);
   int output_interval = std::max(1, static_cast<int>(output_dt / dt));

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
         // Write linear tets only: disables order-p curved geometry and
         // per-element sub-refinement, cutting per-cycle size ~40x.
         pv_out->SetHighOrderOutput(false);
         pv_out->SetLevelsOfDetail(1);
      }

      // Velocity: register a 3-component vector L2 GridFunction.  The
      // scalar-component DOF layout of byNODES (v_c[c * ndof_total + i])
      // matches Q's velocity block (Q[(VX+c) * ndof_total + i]), so one
      // memcpy of 3 * ndof_total * sizeof(real_t) populates the field.
      pv_vel_fec = std::make_unique<L2_FECollection>(order, 3, BasisType::GaussLobatto);
      pv_vel_fes = std::make_unique<PvFES>(&pmesh, pv_vel_fec.get(),
                                           3, Ordering::byNODES);
      pv_vel_gf  = std::make_unique<PvGF>(pv_vel_fes.get());
      *pv_vel_gf = 0.0;
      pv_out->RegisterDomainField("velocity", pv_vel_gf.get());

      // MPI rank: L2 p=0 (one value per element) — same as BP5.
      pv_rank_fec = std::make_unique<L2_FECollection>(0, 3);
      pv_rank_fes = std::make_unique<PvFES>(&pmesh, pv_rank_fec.get());
      pv_rank_gf  = std::make_unique<PvGF>(pv_rank_fes.get());
      *pv_rank_gf = static_cast<real_t>(rank);
      pv_out->RegisterDomainField("mpi_rank", pv_rank_gf.get());

      // Fault L2-p0 field registration. Forward the wave operator's
      // canonical face lists so the ParaView indexing lines up with
      // fault_coords / dof_data exactly (see BP5 pattern).
      pv_out->InitFaultOutputBP5(fault_int_faces, fault_shr_faces,
                                 nqp_per_face);

      // Allocate dynamic-field scratch vectors (vdim=2 for slip/rate/traction).
      pv_local_slip.SetSize(2 * num_fault_total);
      pv_local_slip_rate.SetSize(2 * num_fault_total);
      pv_local_traction.SetSize(2 * num_fault_total);
      pv_local_state.SetSize(num_fault_total);
      pv_local_normal_stress.SetSize(num_fault_total);

      // R-V92-E02 stage-4 buffers (same shape as the averaged ones).
      pv_local_slip_rate_k4.SetSize(2 * num_fault_total);
      pv_local_traction_k4.SetSize(2 * num_fault_total);
      pv_local_normal_stress_k4.SetSize(num_fault_total);

      // Static friction parameters / fault coordinates (one entry per QP).
      // x2 = along-strike, x3 = z (depth, negative below free surface) to
      // match BP5's fault_x2/fault_x3 convention.
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

      // Scheduling: CLI flags take precedence; default to step interval that
      // matches --output-dt so PV frames line up with station output.
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
         pv_out->output_every_n_steps = output_interval;
      }

      if (rank == 0)
      {
         std::cout << "ParaView output: ON (prefix="
                   << output_dir << "/ParaView)\n";
         if (pv_no_domain)
         {
            std::cout << "  Mode: fault-surface PVD only (--no-domain-pv)\n";
         }
         if (debug_qnorm)
         {
            std::cout << "  Diagnostic: --debug-qnorm ON "
                         "(per-rank ||Q||_inf each output cycle)\n";
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
            std::cout << "  Interval: every " << output_interval
                      << " steps (matches --output-dt=" << output_dt
                      << " s)\n";
         }
      }

      // Optional bulk collection: velocity + sigma_yy + sigma_xy +
      // sigma_xz + mpi_rank written at a separate --paraview-bulk-dt
      // cadence to `ParaView_bulk/`.  Independent of --no-domain-pv —
      // that flag suppresses pv_out's fault-schedule volume save,
      // while this collection writes its own PVD/VTU series.
      if (paraview_bulk_dt > 0.0)
      {
         pv_bulk_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
            output_dir + "/ParaView_bulk", pmesh, order);
         if (pv_low_order)
         {
            pv_bulk_out->SetHighOrderOutput(false);
            pv_bulk_out->SetLevelsOfDetail(1);
         }

         // Reuse velocity + mpi_rank GFs from pv_out — they are copied
         // into pv_vel_gf / pv_rank_gf at every paraview_write call.
         pv_bulk_out->RegisterDomainField("velocity", pv_vel_gf.get());
         pv_bulk_out->RegisterDomainField("mpi_rank", pv_rank_gf.get());

         // Scalar L2 order-`order` FESpace for the three stress
         // components extracted from Q's stress block.
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

      // R-007 / R-104 fix: single-shot schedule gate.  PeekShouldWrite is a
      // const read that does not advance last_write_time_.  The schedule is
      // committed below (ForceSave advances it for the !pv_no_domain path;
      // CommitSchedule advances it explicitly for the pv_no_domain path),
      // so the gate and the advance can never disagree on V_max or the
      // current last_write_time_.
      //
      // The bulk collection has its OWN schedule (via
      // `pv_bulk_out->fixed_dt = paraview_bulk_dt`), so we check it
      // independently.  If neither collection wants to fire this step,
      // early-return.  Otherwise fall through, pack the fields needed by
      // whichever collection IS firing, and dispatch both.
      const bool fault_wants =
         pv_out->PeekShouldWrite(step_num, time, V_max);
      const bool bulk_wants = pv_bulk_out &&
         pv_bulk_out->PeekShouldWrite(step_num, time, V_max);
      if (!fault_wants && !bulk_wants) { return; }

      if (!pv_no_domain || bulk_wants)
      {
         // Copy Q's velocity block (VX..VZ, length 3*ndof_total) into vel_gf.
         // byNODES ordering of the vector FES matches Q's component-major layout.
         // Needed by either the fault-schedule domain save (pv_out with
         // --no-domain-pv OFF) or the bulk schedule (pv_bulk_out).
         std::memcpy(pv_vel_gf->GetData(),
                     Q.GetData() + VX * ndof_total,
                     3 * ndof_total * sizeof(real_t));
      }

      // Pack bulk stress-component GFs from Q's SYY / SXY / SXZ blocks.
      // byNODES ordering of the scalar FES matches Q's component-major
      // layout for each of SYY/SXY/SXZ, so one memcpy per component
      // populates the corresponding GF.
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
      }

      // Bulk collection save — advances its own schedule (independent of
      // pv_out's schedule).  Uses the velocity GF packed above plus the
      // three sigma GFs packed above (and the mpi_rank GF which is
      // populated once at setup).  No fault-DOFData packing is needed.
      if (bulk_wants)
      {
         pv_bulk_out->ForceSave(step_num, time);
      }

      // Everything below is fault-schedule-only (fault DOFData packing,
      // fault-surface domain save, fault-surface VTU).  If the fault
      // schedule did not fire this step, nothing else to do.
      if (!fault_wants) { return; }

      // Pack fault fields from DOFData.  R-801 Option A: DOFData.V1/slip1/
      // tau1_corr are the DIP-aligned components and DOFData.V2/slip2/
      // tau2_corr are the STRIKE-aligned components, project-wide (BP5
      // convention).  BP5's ParaView writer takes comp 0 = dip, comp 1 =
      // strike, so this map is now the identity — no swap needed.
      //
      // Pre-R-801 this code swapped components because the SOURCE convention
      // was (t1 = strike, t2 = dip) on interior fault QPs via
      // GodunovFlux::BuildFrame.  Under Option A the source is BP5-canonical
      // everywhere, so the swap is removed (it would now double-invert).
      for (int i = 0; i < num_fault_total; i++)
      {
         const DOFData &d = dof_data[i];
         pv_local_slip(2*i + 0)      = d.slip1;       // dip
         pv_local_slip(2*i + 1)      = d.slip2;       // strike
         pv_local_slip_rate(2*i + 0) = d.V1;          // dip rate
         pv_local_slip_rate(2*i + 1) = d.V2;          // strike rate
         pv_local_traction(2*i + 0)  = d.tau1_corr;   // dip traction
         pv_local_traction(2*i + 1)  = d.tau2_corr;   // strike traction
         pv_local_state(i)           = d.psi;
         pv_local_normal_stress(i)   = d.sigma_n_corr;
      }
      // R-V92-E02 + R-V92-I01 fix: pv_local_*_k4 must carry the stage-4
      // DOFData snapshot captured by the main RK4 loop (line ~919) BEFORE
      // the averaging block overwrites DOFData.  The gate `step_num == 0`
      // restricts the "default = averaged" initialisation to the t=0
      // snapshot only — there is no stage-4 at t=0, so `<field>_k4 -
      // <field> == 0` by construction there, which is correct.  On every
      // subsequent paraview_write call the main-loop write is preserved
      // and the RK4-averaging residual is visible as the delta.  Without
      // this gate the lambda stomped the main-loop write at every output
      // step and the entire R-V92-E02 discriminator was dead.
      if (step_num == 0 && pv_local_normal_stress_k4.Size() == num_fault_total)
      {
         pv_local_slip_rate_k4 = pv_local_slip_rate;
         pv_local_traction_k4  = pv_local_traction;
         pv_local_normal_stress_k4 = pv_local_normal_stress;
      }

      if (pv_no_domain)
      {
         // Fault-surface PVD only; advance the schedule ourselves.
         pv_out->CommitSchedule(time);
      }
      else
      {
         pv_out->UpdateFaultFieldsBP5(pv_local_slip, pv_local_slip_rate,
                                      pv_local_traction, pv_local_state,
                                      pv_local_normal_stress);
         pv_out->ForceSave(step_num, time);  // advances last_write_time_ internally
      }

      // Fault-surface VTU (proper triangle geometry, no L2-p0 scatter) —
      // reached iff PeekShouldWrite returned true above.  R-V92-E02:
      // pass stage-4 buffers so the writer emits `<field>_k4` CellData
      // alongside the averaged `<field>` CellData.  ParaView Calculator
      // "normal_stress_k4 - normal_stress" then quantifies the RK4-
      // averaging residual per triangle.
      pv_out->WriteFaultSurfaceVTU(
         output_dir, step_num, time, rank, nprocs,
         pv_local_slip, pv_local_slip_rate, pv_local_traction,
         pv_local_state, pv_local_normal_stress,
         pv_local_a, pv_local_Dc, pv_local_x2, pv_local_x3,
         pv_local_slip_rate_k4, pv_local_traction_k4,
         pv_local_normal_stress_k4);
   };

   // Initial snapshot at t=0 (V_max = V_ini since the fault is quasi-static).
   paraview_write(0, 0.0, TPV102Params::V_ini);

   // -----------------------------------------------------------------------
   // 8. RK4 time stepping
   // -----------------------------------------------------------------------
   Vector k1(Q.Size()), k2(Q.Size()), k3(Q.Size()), k4(Q.Size());
   Vector Q_tmp(Q.Size());

   real_t t = 0.0;
   real_t V_max_global = 0.0;

   // Pre-allocate RK4 sub-stage storage (reused across steps).
   std::vector<real_t> psi_n(num_fault_total);
   std::vector<real_t> sr_k1(num_fault_total), sr_k2(num_fault_total);
   std::vector<real_t> sr_k3(num_fault_total), sr_k4(num_fault_total);
   std::vector<real_t> V1_k1(num_fault_total), V1_k2(num_fault_total);
   std::vector<real_t> V1_k3(num_fault_total), V1_k4(num_fault_total);
   std::vector<real_t> V2_k1(num_fault_total), V2_k2(num_fault_total);
   std::vector<real_t> V2_k3(num_fault_total), V2_k4(num_fault_total);
   // R-001 fix: Stage-wise corrected tractions for RK4-weighted output
   std::vector<real_t> t1c_k1(num_fault_total), t1c_k2(num_fault_total);
   std::vector<real_t> t1c_k3(num_fault_total), t1c_k4(num_fault_total);
   std::vector<real_t> t2c_k1(num_fault_total), t2c_k2(num_fault_total);
   std::vector<real_t> t2c_k3(num_fault_total), t2c_k4(num_fault_total);
   std::vector<real_t> snc_k1(num_fault_total), snc_k2(num_fault_total);
   std::vector<real_t> snc_k3(num_fault_total), snc_k4(num_fault_total);
   // v9.2.0 plan Step 5 F01+F02 fix: per-stage dpsi/dt buffers so psi
   // integrates inside the RK4 state (classical coupled RK4 on (Q, psi))
   // instead of the pre-fix operator-split pattern (stage-wise analytic
   // updates of psi using single-stage V).  The old pattern was O(dt^2)
   // in the (Q, psi) coupling; coupled RK4 is O(dt^4) to match the Q
   // integration order.
   //
   // Stability note (REVIEW R-V92-H06): UpdateStateAnalytic was
   // unconditionally stable in the constant-V limit via the closed-
   // form theta transform.  Explicit RK4 on dpsi/dt = (b V0 / Dc) *
   // (exp((f0-psi)/b) - V/V0) is only conditionally stable:
   //   local timescale ~ Dc / (b V0 exp((f0-psi)/b))
   //   stability requires dt << that timescale.
   // TPV102 during an event: psi ∈ [0.40, 0.85]; with b=0.012, V0=1e-6,
   // Dc=0.14 the timescale floor is ~0.6 s, and CFL dt ~ 2 ms gives
   // dt/tau ~ 3e-3 — safely stable.  Future SEAS cycle simulations
   // that drive psi below ~0.35 could enter an unstable regime; the
   // psi-floor warning printf below fires if psi drops unexpectedly
   // low.  Revisit if broadening past TPV102.
   std::vector<real_t> psi_k1(num_fault_total), psi_k2(num_fault_total);
   std::vector<real_t> psi_k3(num_fault_total), psi_k4(num_fault_total);
   AgingLawPsi aging_law(TPV102Params::b, TPV102Params::V0, TPV102Params::f0);
   // One-shot stability-envelope tripwire: fires at most once across
   // the whole run if any psi drops below 0.3.  Rank 0 only.
   bool psi_stability_warned = false;

   if (rank == 0)
   {
      std::cout << "Starting time stepping...\n";
      if (use_ader)
      {
         std::cout << "Time integrator: ADER-O(" << ader_order
                   << ")  (default)\n";
      }
      else
      {
         std::cout << "Time integrator: RK4  (alternative; default is ADER)\n";
      }
   }

   // Round-6 R-004: hoist the ADER step buffer out of the loop and use
   // Vector::Swap for an O(1) exchange each step (vs the previous
   // Q = Q_new deep copy).  Allocation is one-shot; the Swap-based
   // exchange also makes Q and Q_new refer to the same pair of
   // buffers throughout the run.
   Vector Q_new(Q.Size());

   for (int step = 0; step < nsteps; step++)
   {
      real_t dt_step = std::min(dt, tfinal - t);
      if (dt_step <= 0.0) { break; }

      // ================================================================
      // ADER I-05 Phase 7: ADER time-stepping branch.  RK4 path below is
      // unchanged (this branch is invoked only when --time-integrator=ader).
      // ================================================================
      if (use_ader)
      {
         // Save psi at step start for the sub-step psi update.
         for (int i = 0; i < num_fault_total; i++)
         {
            psi_n[i] = dof_data[i].psi;
         }

         // Nucleation at the stage midpoint for 2nd-order accuracy.
         // Higher-order nucleation time-integration is deferred (plan
         // Phase 7 §4 "Nucleation timing for higher order").
         // Persistent-prestress channel: overwrite DOFData.tau2_nuc;
         // EvaluateTotal additivity makes this a sustained driver
         // (re-imposed at every Riemann solve, not radiated through Q).
         if (num_fault_total > 0)
         {
            ApplyNucleationTotalPrestress(dof_data, fault_coords,
                                          t + dt_step / 2.0);
         }

         // One-shot ADER predictor-corrector.  AdvanceADER fills the
         // DOFData {V1, V2, slip_rate, tau*_corr, sigma_n_corr} with
         // the friction-solve output on Q̄ = I/dt (= time-averaged bulk
         // Q over [t, t+dt]).  NOTE: this is NOT exactly the time-
         // average of the friction-solve trajectory — the solve is
         // nonlinear — but it matches the one-shot equivalent to
         // O(dt²) per plan §Phase 5 §4 (R-007).
         // Q_new is hoisted above the step loop (round-6 R-004).  Q.Swap
         // exchanges the Vector backing buffers in O(1); no deep copy.
         wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
         Q.Swap(Q_new);

         t += dt_step;

         // Forward-Euler psi update using the time-averaged slip rate.
         // 1st-order in dt but matches the plan's pseudocode
         // (`UpdateStateAnalytic`-style closed-form from
         // psi_n and V_avg); the bulk scheme's 2nd-order accuracy is
         // preserved because psi is a scalar state that contributes only
         // to friction via a Lipschitz-smooth law.  Slip accumulation
         // uses the same time-averaged V directly.
         for (int i = 0; i < num_fault_total; i++)
         {
            const real_t dpsi = aging_law.Rate(dof_data[i].slip_rate,
                                               psi_n[i],
                                               dof_data[i].Dc);
            dof_data[i].psi = psi_n[i] + dt_step * dpsi;

            dof_data[i].slip1 += dof_data[i].V1 * dt_step;
            dof_data[i].slip2 += dof_data[i].V2 * dt_step;

            // psi stability tripwire (mirrors RK4 path).
            if (!psi_stability_warned && rank == 0 && dof_data[i].psi < 0.3)
            {
               std::fprintf(stderr,
                  "[WARNING] psi dropped to %.3f at fault QP %d (t=%.3f s). "
                  "ADER path uses forward-Euler psi update; switch to an "
                  "implicit psi solver if the stability margin shrinks.\n",
                  dof_data[i].psi, i, t);
               psi_stability_warned = true;
            }
         }

         // Under ADER, the DOFData values after AdvanceADER ARE the
         // one-shot endpoint state — there is no "stage-4" snapshot to
         // distinguish from an RK4-weighted average.  Populate the _k4
         // fields with the same values so `<field>_k4 − <field>` = 0 in
         // the VTU (plan Phase 7 §3: "documents that ADER eliminates
         // H-V92-K by formulation").
         if (use_paraview && pv_local_normal_stress_k4.Size() == num_fault_total)
         {
            for (int i = 0; i < num_fault_total; i++)
            {
               pv_local_slip_rate_k4(2*i + 0) = dof_data[i].V1;
               pv_local_slip_rate_k4(2*i + 1) = dof_data[i].V2;
               pv_local_traction_k4(2*i + 0)  = dof_data[i].tau1_corr;
               pv_local_traction_k4(2*i + 1)  = dof_data[i].tau2_corr;
               pv_local_normal_stress_k4(i)   = dof_data[i].sigma_n_corr;
            }
         }

         // V_max tracking — under ADER we only have the time-averaged
         // slip rate per step, so track its peak (vs the per-stage max
         // the RK4 path tracks).  This is a weaker tripwire under ADER
         // but still surfaces unbounded growth.
         real_t V_max_local = 0.0;
         for (int i = 0; i < num_fault_total; i++)
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

         // Output (mirrors the RK4 loop's station/surface write cadence).
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

         // R-004 (round-5): step-0 shared-fault DOFData consistency check
         // must also run under ADER so an MPI cross-rank regression
         // doesn't go silent.  Mirrors the RK4 path's call below.
         if (step == 0)
         {
            wave.VerifySharedFaultDOFDataConsistency();
         }

         // R-002 (round-5): ParaView output cadence must match the RK4
         // path — paraview_write is MPI-collective, so every rank must
         // call it every step.  The pv_local_*_k4 buffers were populated
         // above (so `<field>_k4 - <field>` = 0 under ADER by
         // construction, per plan Phase 7 §5).
         paraview_write(step + 1, t, V_max_step);

         // R-003 (round-5): NaN tripwire — AdvanceADER can produce NaN
         // under CFL-violating dt or a divergent friction solve; mirror
         // the RK4 path's detection and loud exit so a bad run is
         // caught immediately rather than propagated.
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
                         << ", t = " << t << " s (rank " << rank
                         << ", ADER path)\n";
#ifdef MFEM_USE_MPI
               MPI_Finalize();
#endif
               return 1;
            }
         }

         // All per-step epilogue actions handled in the ADER branch.
         // Skip the RK4 bookkeeping below and go back to the top of
         // the step loop.
         continue;
      }

      // ================================================================
      // RK4 path (unchanged from the default TPV102 driver).
      // ================================================================

      // Save psi at step start for sub-stage restoration.
      for (int i = 0; i < num_fault_total; i++)
      {
         psi_n[i] = dof_data[i].psi;
      }

      // RK4 stage 1: at time t, psi = psi_n
      // R-003 fix: nucleation evaluated at stage time.
      // Persistent-prestress channel: overwrite DOFData.tau2_nuc; the
      // EvaluateTotal trial-traction addition makes this a sustained
      // driver re-imposed at every wave-operator Riemann solve.
      if (num_fault_total > 0)
      {
         ApplyNucleationTotalPrestress(dof_data, fault_coords, t);
      }
      wave.Mult(Q, k1);
      for (int i = 0; i < num_fault_total; i++)
      {
         V1_k1[i] = dof_data[i].V1; V2_k1[i] = dof_data[i].V2;
         sr_k1[i] = dof_data[i].slip_rate;
         t1c_k1[i] = dof_data[i].tau1_corr; t2c_k1[i] = dof_data[i].tau2_corr;
         snc_k1[i] = dof_data[i].sigma_n_corr;
         // F01+F02: psi_k1 = dpsi/dt | (V = sr_k1, psi = psi_n, Dc).
         // Advance psi to stage 2 input: psi_n + (dt/2) * psi_k1.
         psi_k1[i] = aging_law.Rate(sr_k1[i], psi_n[i], dof_data[i].Dc);
         dof_data[i].psi = psi_n[i] + 0.5 * dt_step * psi_k1[i];
      }

      // RK4 stage 2: at time t + dt/2, psi = psi_n + (dt/2) * psi_k1
      if (num_fault_total > 0)
      {
         ApplyNucleationTotalPrestress(dof_data, fault_coords,
                                       t + dt_step / 2.0);
      }
      add(Q, dt_step / 2.0, k1, Q_tmp);
      wave.Mult(Q_tmp, k2);
      for (int i = 0; i < num_fault_total; i++)
      {
         V1_k2[i] = dof_data[i].V1; V2_k2[i] = dof_data[i].V2;
         sr_k2[i] = dof_data[i].slip_rate;
         t1c_k2[i] = dof_data[i].tau1_corr; t2c_k2[i] = dof_data[i].tau2_corr;
         snc_k2[i] = dof_data[i].sigma_n_corr;
         // F01+F02: psi_k2 = dpsi/dt at stage-2 state (V = sr_k2, psi = psi_n + dt/2 * psi_k1
         // which is the current dof_data[i].psi set at end of stage 1).
         // Advance psi to stage 3 input: psi_n + (dt/2) * psi_k2.
         psi_k2[i] = aging_law.Rate(sr_k2[i], dof_data[i].psi, dof_data[i].Dc);
         dof_data[i].psi = psi_n[i] + 0.5 * dt_step * psi_k2[i];
      }

      // RK4 stage 3: at time t + dt/2, psi = psi_n + (dt/2) * psi_k2
      add(Q, dt_step / 2.0, k2, Q_tmp);
      wave.Mult(Q_tmp, k3);
      for (int i = 0; i < num_fault_total; i++)
      {
         V1_k3[i] = dof_data[i].V1; V2_k3[i] = dof_data[i].V2;
         sr_k3[i] = dof_data[i].slip_rate;
         t1c_k3[i] = dof_data[i].tau1_corr; t2c_k3[i] = dof_data[i].tau2_corr;
         snc_k3[i] = dof_data[i].sigma_n_corr;
         // F01+F02: psi_k3 = dpsi/dt at stage-3 state (V = sr_k3, psi = psi_n + dt/2 * psi_k2).
         // Advance psi to stage 4 input: psi_n + dt * psi_k3.
         psi_k3[i] = aging_law.Rate(sr_k3[i], dof_data[i].psi, dof_data[i].Dc);
         dof_data[i].psi = psi_n[i] + dt_step * psi_k3[i];
      }

      // RK4 stage 4: at time t + dt, psi = psi_n + dt * psi_k3
      if (num_fault_total > 0)
      {
         ApplyNucleationTotalPrestress(dof_data, fault_coords, t + dt_step);
      }
      add(Q, dt_step, k3, Q_tmp);
      wave.Mult(Q_tmp, k4);
      for (int i = 0; i < num_fault_total; i++)
      {
         V1_k4[i] = dof_data[i].V1; V2_k4[i] = dof_data[i].V2;
         sr_k4[i] = dof_data[i].slip_rate;
         t1c_k4[i] = dof_data[i].tau1_corr; t2c_k4[i] = dof_data[i].tau2_corr;
         snc_k4[i] = dof_data[i].sigma_n_corr;
         // F01+F02: psi_k4 = dpsi/dt at stage-4 state (V = sr_k4, psi = psi_n + dt * psi_k3).
         // No further advance; the final RK4 combination lives in the post-stage block.
         psi_k4[i] = aging_law.Rate(sr_k4[i], dof_data[i].psi, dof_data[i].Dc);
      }

      // Update Q with RK4 weights
      for (int i = 0; i < Q.Size(); i++)
      {
         Q[i] += dt_step / 6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
      }

      t += dt_step;

      // R-V92-E02 (plan §19): capture stage-4 DOFData snapshot BEFORE the
      // RK4-averaging block below overwrites DOFData.V1/V2/tau*_corr/
      // sigma_n_corr with the averaged values.  paraview_write packs these
      // into `<field>_k4` CellData alongside the averaged `<field>` CellData
      // so the RK4-averaging residual `<field>_k4 - <field>` is directly
      // inspectable in ParaView — it is the H-V92-K vs H-V92-G discriminator.
      if (use_paraview && pv_local_normal_stress_k4.Size() == num_fault_total)
      {
         for (int i = 0; i < num_fault_total; i++)
         {
            pv_local_slip_rate_k4(2*i + 0)   = V1_k4[i];
            pv_local_slip_rate_k4(2*i + 1)   = V2_k4[i];
            pv_local_traction_k4(2*i + 0)    = t1c_k4[i];
            pv_local_traction_k4(2*i + 1)    = t2c_k4[i];
            pv_local_normal_stress_k4(i)     = snc_k4[i];
         }
      }

      // v9.2.0 Step 5 F01+F02: classical coupled RK4 on (Q, psi).
      // psi_k1..psi_k4 are dpsi/dt samples at the four RK4 stage states
      // (V, psi) — V from the just-finished Mult, psi from the stage input
      // set at the end of the previous stage block.  The combination below
      // is the O(dt^4) Butcher-tableau weighted sum, replacing the pre-fix
      // O(dt^2) operator-split pattern (analytic update with averaged V).
      // All other fields (V/slip/tau/sigma_n) retain their RK4-weighted
      // averages from the stage buffers — slip IS the RK4 integral of V
      // when V_avg uses (1+2+2+1)/6 weights (since dslip/dt = V), so
      // `slip += V_avg * dt` below is the same thing.
      for (int i = 0; i < num_fault_total; i++)
      {
         dof_data[i].psi = psi_n[i] + dt_step / 6.0 *
                           (psi_k1[i] + 2.0*psi_k2[i] + 2.0*psi_k3[i] + psi_k4[i]);
         // R-V92-H06 stability tripwire: explicit RK4 on aging law is
         // conditionally stable; safe margin fails around psi < 0.3.
         // Warn once on rank 0 to surface if a new scenario enters the
         // unstable regime without changing the integrator.
         if (!psi_stability_warned && rank == 0 && dof_data[i].psi < 0.3)
         {
            std::fprintf(stderr,
               "[WARNING] psi dropped to %.3f at fault QP %d (t=%.3f s). "
               "Explicit RK4 aging-law stability margin shrinks rapidly "
               "below psi=0.35 (see driver comment near aging_law ctor). "
               "If this run is a cycle simulation, switch to an implicit "
               "psi solver.\n",
               dof_data[i].psi, i, t);
            psi_stability_warned = true;
         }
         // R-V92-K01 (round-9): slip IS the RK4 integral of V when
         // V_avg uses the Butcher (1,2,2,1)/6 weights (since dslip/dt
         // = V), so `slip += V_avg * dt` below is the exact O(dt^4)
         // integral.  This is the ONLY place the Simpson-mean V_avg is
         // USED as a value in its own right — for the slip integral.
         const real_t V1_avg = (V1_k1[i] + 2*V1_k2[i] + 2*V1_k3[i] + V1_k4[i]) / 6.0;
         const real_t V2_avg = (V2_k1[i] + 2*V2_k2[i] + 2*V2_k3[i] + V2_k4[i]) / 6.0;
         dof_data[i].slip1 += V1_avg * dt_step;
         dof_data[i].slip2 += V2_avg * dt_step;
      }

      // ---------------------------------------------------------------
      // R-V92-K01 (round-9 REVIEW fix): endpoint re-evaluation of fault
      // observables at Q(t+dt).
      //
      // Prior behaviour: dof_data.{V1, V2, slip_rate, tau1_corr, tau2_corr,
      // sigma_n_corr} were overwritten with the RK4-Butcher-weighted
      // Simpson mean of stage-i snapshots, i.e. the TIME-AVERAGED values
      // over [t_n, t_n+dt] which approximate midpoint (t_n + dt/2).
      //
      // Problem: the station writer labels the sample time as t = t_{n+1}
      // (post-increment above), but the values correspond to t_n + dt/2 —
      // a half-step phase lag.  Semantically wrong for instantaneous
      // observables regardless of magnitude.
      //
      // Fix: invoke wave.Mult(Q, k_unused) one more time on the ENDPOINT
      // state Q(t_{n+1}).  This triggers FaultFaceFlux::Evaluate on the
      // endpoint Q, which overwrites dof_data.{tau1_corr, tau2_corr,
      // sigma_n_corr, V1, V2, slip_rate} with values self-consistent
      // with the Q used for output.  The slip accumulators already hold
      // the correct Simpson-integral of V (above) — do NOT re-accumulate.
      //
      // Cost: one extra Mult per step.  Not wrapped under `if output step`
      // because every subsequent RK4 step's stage-1 would do the same
      // Mult anyway; this just shifts that first Mult from "stage 1 of
      // step n+1" to "endpoint re-eval of step n", reusing the result.
      // To avoid the overhead we could cache k_endpoint and reuse it as
      // stage-1 k1 of the next step (FSAL-style), but that complicates
      // nucleation timing — deferred.
      // ---------------------------------------------------------------
      {
         Vector k_endpoint(Q.Size());
         // Persistent-prestress channel: re-impose tau2_nuc at the new
         // step time so the endpoint Mult sees the correct nucleation
         // amplitude.  Idempotent overwrite — safe to call after
         // stage-4 even though stage-4 already set tau2_nuc(t+dt).
         if (num_fault_total > 0)
         {
            ApplyNucleationTotalPrestress(dof_data, fault_coords, t);
         }
         wave.Mult(Q, k_endpoint);
         // dof_data[i].{tau1_corr, tau2_corr, sigma_n_corr, V1, V2,
         //              slip_rate} are now the endpoint-at-t values
         // produced by the stage Evaluate on Q(t+dt).
      }

      // R-006: Track peak V_max across all RK4 stages, not just the average
      real_t V_max_local = 0.0;
      for (int i = 0; i < num_fault_total; i++)
      {
         V_max_local = std::max(V_max_local,
            std::max({sr_k1[i], sr_k2[i], sr_k3[i], sr_k4[i]}));
      }
#ifdef MFEM_USE_MPI
      real_t V_max_step;
      MPI_Allreduce(&V_max_local, &V_max_step, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#else
      real_t V_max_step = V_max_local;
#endif
      V_max_global = std::max(V_max_global, V_max_step);

#ifdef SEAS_DIAG_FAULT_FLUX
      {
         // C-4 BULK: v9.0.0 §0.5 checkpoint — global max|Q[VX]|, max|Q[SXY]|
         // post-RK4 update.  Deadlock-safe by construction:
         //   * `c4_fire` depends only on (step, t) — both lockstep-identical
         //     across ranks — so every rank evaluates the same boolean.
         //   * MPI_Allreduce is called UNCONDITIONALLY on every rank
         //     whenever c4_fire is true; no rank-local gate on the
         //     collective.
         //   * Only the final fprintf is guarded by `rank == 0`, AFTER the
         //     collectives complete.
         // Throttle: every 100 steps routinely + every 10 steps inside the
         // 1.0 s <= t < 1.4 s breakaway window.  step=0 fires unconditionally
         // so the 2-rank deadlock unit-test sees at least one C-4 line.
         const bool c4_fire = (step % 100 == 0) ||
                              (step % 10 == 0 && t >= 1.0 && t < 1.4);
         if (c4_fire)
         {
            real_t q_vx_local  = 0.0;
            real_t q_sxy_local = 0.0;
            for (int i = 0; i < ndof_total; i++)
            {
               q_vx_local  = std::max(q_vx_local,
                                      std::abs(Q[VX  * ndof_total + i]));
               q_sxy_local = std::max(q_sxy_local,
                                      std::abs(Q[SXY * ndof_total + i]));
            }
            real_t q_vx_global  = q_vx_local;
            real_t q_sxy_global = q_sxy_local;
#ifdef MFEM_USE_MPI
            MPI_Allreduce(&q_vx_local,  &q_vx_global,  1,
                          MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
            MPI_Allreduce(&q_sxy_local, &q_sxy_global, 1,
                          MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
            if (rank == 0)
            {
               std::fprintf(stderr,
                  "[C-4 BULK] step=%d  t=%.4f  V_max=%.3e  "
                  "max|Q[VX]|=%.3e m/s  max|Q[SXY]|=%.3e Pa\n",
                  step, t, V_max_step, q_vx_global, q_sxy_global);
            }
         }
      }
#endif


      // R-101 fix: after the first RK4 step, verify that shared-fault
      // DOFData entries agree bit-identically across the two ranks that
      // own each shared fault QP.  The (+,-) canonicalisation in R-001
      // depends on an MFEM invariant (identical face normals on both
      // ranks) that was previously unverified; this call makes the
      // invariant failure mode loud and fatal instead of silent drift.
      if (step == 0)
      {
         wave.VerifySharedFaultDOFDataConsistency();
      }

      // Output
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

         // --debug-qnorm diagnostic (tpv102_debug_v1 H1): print per-rank
         // ||Q||_inf.  If all ranks other than the hypocenter's stay at 0,
         // bulk wave energy is not crossing partition seams.
         if (debug_qnorm)
         {
            real_t qn_local = Q.Normlinf();
#ifdef MFEM_USE_MPI
            // Gather all ranks' norms onto rank 0 for a compact summary.
            // R-303 fix: MPI datatype must match real_t at compile time.
            // Hardcoding MPI_DOUBLE silently corrupts qn_all on
            // MFEM_USE_SINGLE builds (where real_t = float = 4 bytes).
            std::vector<real_t> qn_all;
            if (rank == 0) { qn_all.resize(nprocs); }
            MPI_Gather(&qn_local, 1, MPITypeMap<real_t>::mpi_type,
                       rank == 0 ? qn_all.data() : nullptr, 1,
                       MPITypeMap<real_t>::mpi_type,
                       0, comm);
            if (rank == 0)
            {
               real_t qmin = qn_all[0], qmax = qn_all[0], qsum = 0.0;
               int nzero = 0;
               for (int r = 0; r < nprocs; r++)
               {
                  qmin = std::min(qmin, qn_all[r]);
                  qmax = std::max(qmax, qn_all[r]);
                  qsum += qn_all[r];
                  if (qn_all[r] == 0.0) { nzero++; }
               }
               std::cout << "  [qnorm] min=" << qmin
                         << " max=" << qmax
                         << " mean=" << (qsum / nprocs)
                         << " #ranks_with_||Q||=0: " << nzero
                         << "/" << nprocs << "\n";
               // R-002 / R-106 fix: per-rank breakout for H1 diagnosis.
               // Always show bilateral neighbours of the hypocenter rank
               // plus rank 0 and nprocs-1; this survives the edge case
               // where hypo_rank is at an endpoint of the rank range
               // (previously collapsed the watch list to a single entry).
               std::vector<int> watch = { hypo_rank };
               for (int off : {1, 4})
               {
                  if (hypo_rank + off <  nprocs) { watch.push_back(hypo_rank + off); }
                  if (hypo_rank - off >= 0)      { watch.push_back(hypo_rank - off); }
               }
               watch.push_back(nprocs - 1);
               watch.push_back(0);
               // De-dup while preserving order (small N so linear is fine).
               std::vector<int> watch_unique;
               for (int r : watch)
               {
                  bool dup = false;
                  for (int u : watch_unique) { if (u == r) { dup = true; break; } }
                  if (!dup) { watch_unique.push_back(r); }
               }
               // R-402 fix: format qnorm in scientific notation so the
               // post-run RESULT.txt regex operates on a stable format.
               // Default operator<< switches between fixed and scientific
               // by value magnitude, so a legitimate small ||Q||_inf like
               // 0.0001 would otherwise print as "0.0001" and trip the
               // dead-rank detector.
               std::ios::fmtflags prev_flags = std::cout.flags();
               std::streamsize    prev_prec  = std::cout.precision();
               std::cout << std::scientific << std::setprecision(3);
               std::cout << "  [qnorm:watch]";
               for (int r : watch_unique)
               {
                  std::cout << " r" << r << "=" << qn_all[r];
               }
               std::cout << " (hypo_rank=" << hypo_rank << ")\n";
               std::cout.flags(prev_flags);
               std::cout.precision(prev_prec);
            }
#else
            if (rank == 0)
            {
               std::cout << "  [qnorm] ||Q||_inf = " << qn_local << "\n";
            }
#endif
         }
      }

      // ParaView: MPI-collective; must be called every step (even if the
      // schedule gate inside Save() skips the actual write) so every rank
      // stays in lockstep.  V_max_step is already globally reduced.
      paraview_write(step + 1, t, V_max_step);

      // R-003: NaN detection reduced across all ranks (prevents deadlock)
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
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 1;
      }
   }

   // -----------------------------------------------------------------------
   // 9. Summary
   // -----------------------------------------------------------------------
   station_writer.Flush();
   surface_writer.Flush();

   if (rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "Simulation complete.\n";
      std::cout << "Final time: " << t << " s\n";
      std::cout << "Steps: " << nsteps << "\n";
      std::cout << "Max slip rate: " << V_max_global << " m/s\n";
      std::cout << "Ranks: " << nprocs << "\n";
      std::cout << "Output: " << output_dir << "/\n";
      std::cout << "========================================\n";
   }

#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return 0;
}
