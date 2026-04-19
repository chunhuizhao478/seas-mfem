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
#include "../config/tpv102_params.hpp"
#include "../domain/boundary_config.hpp"
#include "../io/paraview_output.hpp"

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

   // -----------------------------------------------------------------------
   // Parse command-line arguments
   // -----------------------------------------------------------------------
   std::string mesh_file = GetStringArg(argc, argv, "--mesh",
                                        "tpv102/mesh/tpv102_coarse.msh");
   real_t mesh_scale = GetRealArg(argc, argv, "--mesh-scale", 1.0);
   int order = GetIntArg(argc, argv, "--order", 2);
   std::string bc_mode = GetStringArg(argc, argv, "--bc-mode", "absorbing");
   real_t tfinal = GetRealArg(argc, argv, "--tfinal", TPV102Params::t_final);
   std::string output_dir = GetStringArg(argc, argv, "--output-dir", "tpv102/results");
   std::string output_prefix = GetStringArg(argc, argv, "--output-prefix", "tpv102");
   real_t cfl_factor = GetRealArg(argc, argv, "--cfl", 0.5);
   int bc_free = GetIntArg(argc, argv, "--bc-free", 1);
   int bc_fault = GetIntArg(argc, argv, "--bc-fault", 3);
   int bc_absorb = GetIntArg(argc, argv, "--bc-absorb", 5);

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
   }

   FaultFaceFlux fault_flux(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   // R-002 fix: resolve the rank that owns the hypocenter QP (closest local
   // fault QP to (hypo_along_strike, -hypo_down_dip) in x/z).  Used only by
   // --debug-qnorm to print a per-rank ||Q||_inf watch list.
   //
   // R-105 fix: use a named struct with static_asserts so the MPI_DOUBLE_INT
   // layout assumption fails loudly at compile time if it is ever broken
   // (e.g. by a compiler with unusual padding of {double, int}).
   int hypo_rank = 0;
   {
      real_t local_min_dist2 = std::numeric_limits<real_t>::max();
      for (int i = 0; i < num_fault_total; i++)
      {
         real_t dx = fault_coords[i](0) - TPV102Params::hypo_along_strike;
         real_t dz = std::abs(fault_coords[i](2)) - TPV102Params::hypo_down_dip;
         real_t d2 = dx*dx + dz*dz;
         if (d2 < local_min_dist2) { local_min_dist2 = d2; }
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

   // -----------------------------------------------------------------------
   // 6. Initialize state Q = 0 (perturbation field)
   // -----------------------------------------------------------------------
   Vector Q;
   InitializeState(Q, ndof_total);

   if (rank == 0)
   {
      std::cout << "Q = 0 (perturbation). Background: tau_strike = "
                << TPV102Params::tau_ini / 1e6
                << " MPa (BP5 convention: stored in tau2_0), sigma_n = "
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

   // Scratch vectors passed to UpdateFaultFieldsBP5 / WriteFaultSurfaceVTU.
   Vector pv_local_slip, pv_local_slip_rate, pv_local_traction;
   Vector pv_local_state, pv_local_normal_stress;
   Vector pv_local_a, pv_local_Dc, pv_local_x2, pv_local_x3;

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
      if (!pv_out->PeekShouldWrite(step_num, time, V_max)) { return; }

      if (!pv_no_domain)
      {
         // Copy Q's velocity block (VX..VZ, length 3*ndof_total) into vel_gf.
         // byNODES ordering of the vector FES matches Q's component-major layout.
         std::memcpy(pv_vel_gf->GetData(),
                     Q.GetData() + VX * ndof_total,
                     3 * ndof_total * sizeof(real_t));
      }

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
      // reached iff PeekShouldWrite returned true above.
      pv_out->WriteFaultSurfaceVTU(
         output_dir, step_num, time, rank, nprocs,
         pv_local_slip, pv_local_slip_rate, pv_local_traction,
         pv_local_state, pv_local_normal_stress,
         pv_local_a, pv_local_Dc, pv_local_x2, pv_local_x3);
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

   if (rank == 0) { std::cout << "Starting time stepping...\n"; }

   for (int step = 0; step < nsteps; step++)
   {
      real_t dt_step = std::min(dt, tfinal - t);
      if (dt_step <= 0.0) { break; }

      // Save psi at step start for sub-stage restoration.
      for (int i = 0; i < num_fault_total; i++)
      {
         psi_n[i] = dof_data[i].psi;
      }

      // RK4 stage 1: at time t, psi = psi_n
      // R-003 fix: nucleation evaluated at stage time
      if (num_fault_total > 0) { ApplyNucleation(dof_data, num_fault_total, fault_coords, t); }
      wave.Mult(Q, k1);
      for (int i = 0; i < num_fault_total; i++)
      {
         V1_k1[i] = dof_data[i].V1; V2_k1[i] = dof_data[i].V2;
         sr_k1[i] = dof_data[i].slip_rate;
         t1c_k1[i] = dof_data[i].tau1_corr; t2c_k1[i] = dof_data[i].tau2_corr;
         snc_k1[i] = dof_data[i].sigma_n_corr;
         // Advance psi to t + dt/2 for stage 2
         dof_data[i].psi = UpdateStateAnalytic(psi_n[i], sr_k1[i], dof_data[i].Dc,
            dt_step / 2.0, TPV102Params::f0, TPV102Params::b, TPV102Params::V0);
      }

      // RK4 stage 2: at time t + dt/2, psi advanced by dt/2 using k1 rate
      if (num_fault_total > 0) { ApplyNucleation(dof_data, num_fault_total, fault_coords, t + dt_step / 2.0); }
      add(Q, dt_step / 2.0, k1, Q_tmp);
      wave.Mult(Q_tmp, k2);
      for (int i = 0; i < num_fault_total; i++)
      {
         V1_k2[i] = dof_data[i].V1; V2_k2[i] = dof_data[i].V2;
         sr_k2[i] = dof_data[i].slip_rate;
         t1c_k2[i] = dof_data[i].tau1_corr; t2c_k2[i] = dof_data[i].tau2_corr;
         snc_k2[i] = dof_data[i].sigma_n_corr;
         // Re-advance psi from psi_n using k2 rate for dt/2
         dof_data[i].psi = UpdateStateAnalytic(psi_n[i], sr_k2[i], dof_data[i].Dc,
            dt_step / 2.0, TPV102Params::f0, TPV102Params::b, TPV102Params::V0);
      }

      // RK4 stage 3: at time t + dt/2, psi advanced by dt/2 using k2 rate
      add(Q, dt_step / 2.0, k2, Q_tmp);
      wave.Mult(Q_tmp, k3);
      for (int i = 0; i < num_fault_total; i++)
      {
         V1_k3[i] = dof_data[i].V1; V2_k3[i] = dof_data[i].V2;
         sr_k3[i] = dof_data[i].slip_rate;
         t1c_k3[i] = dof_data[i].tau1_corr; t2c_k3[i] = dof_data[i].tau2_corr;
         snc_k3[i] = dof_data[i].sigma_n_corr;
         // Advance psi from psi_n using k3 rate for full dt
         dof_data[i].psi = UpdateStateAnalytic(psi_n[i], sr_k3[i], dof_data[i].Dc,
            dt_step, TPV102Params::f0, TPV102Params::b, TPV102Params::V0);
      }

      // RK4 stage 4: at time t + dt, psi advanced by dt using k3 rate
      if (num_fault_total > 0) { ApplyNucleation(dof_data, num_fault_total, fault_coords, t + dt_step); }
      add(Q, dt_step, k3, Q_tmp);
      wave.Mult(Q_tmp, k4);
      for (int i = 0; i < num_fault_total; i++)
      {
         V1_k4[i] = dof_data[i].V1; V2_k4[i] = dof_data[i].V2;
         sr_k4[i] = dof_data[i].slip_rate;
         t1c_k4[i] = dof_data[i].tau1_corr; t2c_k4[i] = dof_data[i].tau2_corr;
         snc_k4[i] = dof_data[i].sigma_n_corr;
      }

      // Update Q with RK4 weights
      for (int i = 0; i < Q.Size(); i++)
      {
         Q[i] += dt_step / 6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
      }

      t += dt_step;

      // State update using RK4-weighted average slip rate.
      // NOTE: This gives O(dt^2) coupling accuracy for the wave+friction system.
      // The analytic update with averaged V loses RK4's higher-order corrections.
      // For CFL-limited dt ~ 0.02 ms on 100m meshes, cumulative psi error over
      // 12 s is ~ (dt)^2 * nsteps ~ 3e-4, negligible vs spatial O(h) error.
      // Full O(dt^4) coupling requires integrating psi inside the RK4 state vector.
      for (int i = 0; i < num_fault_total; i++)
      {
         real_t sr_avg = (sr_k1[i] + 2.0*sr_k2[i] + 2.0*sr_k3[i] + sr_k4[i]) / 6.0;
         dof_data[i].psi = UpdateStateAnalytic(psi_n[i], sr_avg, dof_data[i].Dc,
            dt_step, TPV102Params::f0, TPV102Params::b, TPV102Params::V0);
         dof_data[i].V1 = (V1_k1[i] + 2*V1_k2[i] + 2*V1_k3[i] + V1_k4[i]) / 6.0;
         dof_data[i].V2 = (V2_k1[i] + 2*V2_k2[i] + 2*V2_k3[i] + V2_k4[i]) / 6.0;
         // R-005: Recompute slip_rate from averaged V components for output consistency
         dof_data[i].slip_rate = std::sqrt(dof_data[i].V1 * dof_data[i].V1
                                         + dof_data[i].V2 * dof_data[i].V2);
         dof_data[i].slip1 += dof_data[i].V1 * dt_step;
         dof_data[i].slip2 += dof_data[i].V2 * dt_step;
         // R-001 fix: RK4-weighted corrected tractions for consistent station output
         dof_data[i].tau1_corr = (t1c_k1[i] + 2*t1c_k2[i] + 2*t1c_k3[i] + t1c_k4[i]) / 6.0;
         dof_data[i].tau2_corr = (t2c_k1[i] + 2*t2c_k2[i] + 2*t2c_k3[i] + t2c_k4[i]) / 6.0;
         dof_data[i].sigma_n_corr = (snc_k1[i] + 2*snc_k2[i] + 2*snc_k3[i] + snc_k4[i]) / 6.0;
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
