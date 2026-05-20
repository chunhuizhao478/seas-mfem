// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// seas_spatial_dyn_driver — Phase 4 of
// safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md (rev-3).
//
// Single-event dynamic-rupture driver for spatially-varying material,
// pre-stress, and friction inputs.  Reuses the TPV205 LSW closed-form
// machinery (Tpv205SubStepIterator + AdvanceADERWithSubStep + the new
// LSW_ForcedRupture dispatch from Phase H) on a SAFS / TPV26-27-style
// problem driven by a single TOML config.
//
// Deviation block (vs Phase 4 spec).  Documented up front for the
// adversarial code review (REVIEW.md).
//
// D-1.  WaveOperator(MaterialField, BoundaryConfig) ctor (Phase H.1)
//       is NOT yet wired (see dynamic/wave_operator.hpp:215-233,
//       SetGodunovFluxPool still aborts).  Until Phase H lands the
//       per-element flux dispatch, this driver constructs WaveOperator
//       through the scalar-material ctor only:
//
//         * MaterialField::Mode::Constant  ⇒ extract (lambda, mu, rho)
//           and call the scalar ctor (byte-equivalent to the
//           heterogeneous ctor on a Constant input).
//         * MaterialField::Mode::Coefficient (the SAFS sidecar path) ⇒
//           ABORT with a clear "Phase H not yet wired" message so the
//           caller cannot accidentally run sidecar-input physics on
//           the scalar code path.
//
//       --no-sidecar-material always uses the Constant path.  The CSM
//       stress sidecar path (Phase 3) is unaffected — pre-stress is a
//       *fault-DOF* quantity and runs through FaultGeometry, not the
//       WaveOperator material.
//
// D-2.  spatial::ComputePerDOFCoordsAndBasisFromWave (Phase 5 helper
//       was deferred — see dynamic/spatial_setup.hpp:14-22).  The
//       driver therefore mirrors the TPV205 driver's inline fault-DOF
//       walk to build per-DOF coordinates / basis / dof_to_elem from
//       wave.GetFaultInteriorFaces() + GetFaultSharedFaces().
//
// D-3.  Write/ReadTpv104Checkpoint actual signature is
//       (prefix, t, dt, step, Q, dof_data, mpi/(rank,size,comm),
//        driver_tag).  The plan's docstring referenced a hypothetical
//       6-extra-ParaView-state form; the actual on-disk schema is
//       smaller.  This driver uses the real raw-MPI overload with
//       driver_tag = "spatial_dyn".

#include "mfem.hpp"

#include "../dynamic/wave_state.hpp"
#include "../dynamic/wave_operator.hpp"
#include "../dynamic/fault_face_flux.hpp"
#include "../dynamic/friction_solver.hpp"
#include "../dynamic/tpv205_friction.hpp"
#include "../dynamic/tpv205_substep_iterator.hpp"
#include "../dynamic/tpv102_substep_iterator.hpp"
#include "../dynamic/tpv104_substep_iterator.hpp"
// REVIEW R-009: SCEC-trace station writer parity with tpv205_driver.cpp.
#include "../dynamic/tpv205_setup.hpp"
// REVIEW R-007: SCEC-trace station writers for TPV102 (rate-and-state
// aging law) and TPV104 (slip-law strong rate weakening).  Headers
// expose `TPV102StationWriter` + `DefaultStations()` (TPV102) and
// `TPV104StationWriter` + `DefaultStations_TPV104()`.  The unsuffixed
// `DefaultStations()` and `InitializeFaultDOFs()` from tpv102_setup.hpp
// do not collide with any other symbol in this translation unit.
#include "../dynamic/tpv102_setup.hpp"
#include "../dynamic/tpv104_setup.hpp"
#include "../friction/state_evolution.hpp"
#include "../friction/slip_law_srw_psi.hpp"
#include "../dynamic/heterogeneous_material.hpp"
#include "../dynamic/spatial_setup.hpp"
#include "../dynamic/seas_diag_rank.hpp"

#include "../domain/boundary_config.hpp"

#include "../fault/fault_basis.hpp"
#include "../fault/fault_geometry.hpp"
#include "../fault/fault_geometry_safs_templated.inl"

#include "../config/bp5_params.hpp"

#include "../common/mpi_context.hpp"

#include "../io/paraview_output.hpp"
#include "../io/tpv104_checkpoint.hpp"
#include "../io/data_field_3d.hpp"
#include "../io/stress_field_3d.hpp"
#include "../io/material_coefficients.hpp"

#include "../dynamic/spatial_nucleation.hpp"
#include "../dynamic/spatial_print_derived.hpp"

#include "../spatial/code/spatial_friction.hpp"
#include "../spatial/code/spatial_velocity.hpp"
#include "../spatial/code/spatial_stress.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

#ifdef SEAS_DIAG_FAULT_FLUX
namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT
#endif

// Precondition P-2 build-time guard (plan §Preconditions, rev-3).
// MaterialField::EvalAt must accept (elem, T, ip, lam, mu, rho) by ref.
static_assert(
   std::is_invocable_v<decltype(&mfem::seas::MaterialField::EvalAt),
                       const mfem::seas::MaterialField*, int,
                       mfem::ElementTransformation&,
                       const mfem::IntegrationPoint&,
                       mfem::real_t&, mfem::real_t&, mfem::real_t&>,
   "Precondition P-2 not met: MaterialField::EvalAt signature drift.");

// --------------------------------------------------------------------------
// Small CLI parsing helpers (same convention as tpv104_driver.cpp).
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

// Map TOML mixed_flux string -> WaveOperator enum.
MixedFluxMode ParseMixedFlux(const std::string &s)
{
   if (s == "none")            { return MixedFluxMode::None; }
   if (s == "adjacent")        { return MixedFluxMode::Adjacent; }
   if (s == "all_continuous")  { return MixedFluxMode::AllContinuous; }
   MFEM_ABORT("spatial_dyn_driver: --mixed-flux: unknown value '" << s
              << "'.  Accepted: none | adjacent | all_continuous.");
}

ParaViewOutput<ParMesh>::VolumeOutputMode ParseVolumeMode(
   const std::string &s, bool &enabled)
{
   if (s == "off") { enabled = false; return ParaViewOutput<ParMesh>::DefaultVolumeOutputMode(); }
   enabled = true;
   if (s == "vtu")  { return ParaViewOutput<ParMesh>::VolumeOutputMode::Vtu; }
   if (s == "hdf5") { return ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5; }
   MFEM_ABORT("spatial_dyn_driver: paraview_volume: unknown value '"
              << s << "'.  Accepted: hdf5 | vtu | off.");
}

ParaViewOutput<ParMesh>::FaultOutputMode ParseFaultMode(
   const std::string &s, bool &enabled)
{
   if (s == "off") { enabled = false; return ParaViewOutput<ParMesh>::FaultOutputMode::Vtu; }
   enabled = true;
   if (s == "vtu")  { return ParaViewOutput<ParMesh>::FaultOutputMode::Vtu; }
   if (s == "hdf5") { return ParaViewOutput<ParMesh>::FaultOutputMode::Hdf5; }
   MFEM_ABORT("spatial_dyn_driver: paraview_fault: unknown value '"
              << s << "'.  Accepted: hdf5 | vtu | off.");
}

// Compute the count of fault QPs per face using the same probe logic as
// drivers/tpv205_driver.cpp:1265-1284: peek at the first interior or
// shared fault face on this rank, then MPI_Allreduce(MAX) so every
// rank agrees.
int ProbeNbfPerFace(ParMesh &pmesh, int order,
                    const Array<int> &fault_int_faces,
                    const Array<int> &fault_shr_faces,
                    MPI_Comm comm)
{
   int nbf = 0;
   if (fault_int_faces.Size() > 0)
   {
      FaceElementTransformations *ftr =
         pmesh.GetInteriorFaceTransformations(fault_int_faces[0]);
      MFEM_VERIFY(ftr, "spatial_dyn: fault interior face has null FTR");
      nbf = IntRules.Get(ftr->GetGeometryType(), 2 * order).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   else if (fault_shr_faces.Size() > 0)
   {
      FaceElementTransformations *ftr =
         pmesh.GetSharedFaceTransformations(fault_shr_faces[0]);
      MFEM_VERIFY(ftr, "spatial_dyn: fault shared face has null FTR");
      nbf = IntRules.Get(ftr->GetGeometryType(), 2 * order).GetNPoints();
   }
   {
      int local_nbf = nbf;
      MPI_Allreduce(&local_nbf, &nbf, 1, MPI_INT, MPI_MAX, comm);
   }
#else
   (void)comm;
#endif
   return nbf;
}

// Walk the wave operator's fault face lists and produce the per-DOF
// arrays the new FaultGeometry ctor / spatial_setup.hpp expect.
//
// Outputs (all owned-fault, no double-counting; interior first, shared
// appended in fault_shared_faces order):
//   - phys_coords[i] = (x, y, z) at fault QP i
//   - dof_basis(9, i) = [n; t1=dip; t2=strike] at QP i
//   - dof_to_elem[i] = bulk Elem1 owning QP i
//   - dof_to_attr[i] = mesh boundary attribute at QP i (fault_attr for
//     interior fault faces; shared faces report the bdr attr if any,
//     else fault_attr)
//   - dof_ips[i] = reference-element IntegrationPoint at QP i
//     (FaultGeometry::fault_dof_ip cache; consumed by
//     spatial_setup.hpp's IP-aware overload).
void BuildPerDOFFaultTables(ParMesh &pmesh,
                            int order,
                            int fault_attr,
                            const FaultBasis &fbasis,
                            const Array<int> &fault_int_faces,
                            const Array<int> &fault_shr_faces,
                            int nbf_per_face,
                            std::vector<Vector> &phys_coords,
                            Vector &dof_coords_3d,
                            DenseMatrix &dof_basis,
                            Array<int> &dof_to_elem,
                            Array<int> &dof_to_attr,
                            std::vector<IntegrationPoint> &dof_ips)
{
   const int num_int  = fault_int_faces.Size() * nbf_per_face;
   const int num_shr  = fault_shr_faces.Size() * nbf_per_face;
   const int N        = num_int + num_shr;

   phys_coords.clear();
   phys_coords.reserve(N);
   dof_coords_3d.SetSize(3 * N);
   dof_basis.SetSize(9, N);
   dof_to_elem.SetSize(N);
   dof_to_attr.SetSize(N);
   dof_ips.assign(N, IntegrationPoint());

   auto write_dof = [&](int dof_idx, FaceElementTransformations *ftr,
                        const IntegrationPoint &ip,
                        const FaultBasisData &bdata,
                        int attr)
   {
      ftr->SetAllIntPoints(&ip);
      Vector phys(3);
      ftr->Face->Transform(ip, phys);
      phys_coords.push_back(phys);

      for (int d = 0; d < 3; ++d)
      {
         dof_coords_3d(3 * dof_idx + d) = phys(d);
         dof_basis(0 + d, dof_idx) = bdata.normal[d];
         dof_basis(3 + d, dof_idx) = bdata.tangent1[d];
         dof_basis(6 + d, dof_idx) = bdata.tangent2[d];
      }
      dof_to_elem[dof_idx] = ftr->Elem1No;
      dof_to_attr[dof_idx] = attr;
      dof_ips[dof_idx]     = ftr->GetElement1IntPoint();
   };

   int dof_idx = 0;

   for (int fi = 0; fi < fault_int_faces.Size(); ++fi)
   {
      const int face = fault_int_faces[fi];
      FaceElementTransformations *ftr =
         pmesh.GetInteriorFaceTransformations(face);
      MFEM_VERIFY(ftr, "spatial_dyn: interior fault face has null FTR");
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2 * order);
      MFEM_VERIFY(ir.GetNPoints() == nbf_per_face,
                  "spatial_dyn: interior fault face has " << ir.GetNPoints()
                  << " QPs, expected " << nbf_per_face);
      const FaultBasisData &bdata = fbasis.GetBasis(fi);
      for (int q = 0; q < nbf_per_face; ++q)
      {
         write_dof(dof_idx++, ftr, ir.IntPoint(q), bdata, fault_attr);
      }
   }
#ifdef MFEM_USE_MPI
   for (int si = 0; si < fault_shr_faces.Size(); ++si)
   {
      const int sf = fault_shr_faces[si];
      FaceElementTransformations *ftr =
         pmesh.GetSharedFaceTransformations(sf);
      MFEM_VERIFY(ftr, "spatial_dyn: shared fault face has null FTR");
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2 * order);
      MFEM_VERIFY(ir.GetNPoints() == nbf_per_face,
                  "spatial_dyn: shared fault face has " << ir.GetNPoints()
                  << " QPs, expected " << nbf_per_face);
      const FaultBasisData &bdata =
         fbasis.GetBasis(fault_int_faces.Size() + si);
      for (int q = 0; q < nbf_per_face; ++q)
      {
         write_dof(dof_idx++, ftr, ir.IntPoint(q), bdata, fault_attr);
      }
   }
#else
   (void)fault_shr_faces;
#endif
   MFEM_VERIFY(dof_idx == N,
               "spatial_dyn: per-DOF table walk wrote " << dof_idx
               << " entries, expected " << N);
}

// ADER macro-step driver: predictor in the bulk, per-QP fault solve
// per sub-step, corrector via wave.AdvanceADER with the side-channel
// I_imp.  Friction-law agnostic: the iterator-specific call (LSW closed
// form / RS aging law / RS slip-law SRW) is wrapped by the
// `do_iterate` lambda so this helper does not care which iterator is
// driving the per-QP solve.
//
// `set_substeps(deltaT_scaled, weights)` configures the wrapped
// iterator's quadrature each call (mirrors what was done inline before
// the refactor — the iterator owns its own quadrature state).
// `do_iterate(Q_pointwise_plus, Q_pointwise_minus, dt_step,
//             t_step_start, I_imp_plus, I_imp_minus, nuc_cb)`
// performs the iterator-specific advance.
void AdvanceADERWithSubStep_Spatial(
   WaveOperator<ParMesh> &wave,
   const std::function<void(const std::vector<real_t>&,
                            const std::vector<real_t>&)> &set_substeps,
   const std::function<void(const std::vector<std::vector<real_t>>&,
                            const std::vector<std::vector<real_t>>&,
                            real_t, real_t, real_t*, real_t*,
                            const std::function<void(real_t, real_t)>&)>
                            &do_iterate,
   const std::vector<real_t> &configured_deltaT,
   const std::vector<real_t> &configured_weights,
   const Vector &Q,
   real_t dt_step,
   int ader_order,
   real_t t_step_start,
   Vector &Q_new,
   const std::function<void(real_t, real_t)> &nuc_callback)
{
   MFEM_VERIFY(dt_step > 0.0,
               "AdvanceADERWithSubStep_Spatial: dt_step must be > 0, got "
               << dt_step);
   MFEM_VERIFY(ader_order >= 2 && ader_order <= 4,
               "AdvanceADERWithSubStep_Spatial: ader_order must be in "
               "{2,3,4}, got " << ader_order);

   const int O = static_cast<int>(configured_deltaT.size());
   MFEM_VERIFY(O >= 1,
               "AdvanceADERWithSubStep_Spatial: configured_deltaT empty");
   MFEM_VERIFY(static_cast<int>(configured_weights.size()) == O,
               "AdvanceADERWithSubStep_Spatial: configured_weights / "
               "configured_deltaT size mismatch");
   const real_t configured_sum =
      std::accumulate(configured_deltaT.begin(), configured_deltaT.end(),
                      static_cast<real_t>(0));
   MFEM_VERIFY(configured_sum > 0.0,
               "AdvanceADERWithSubStep_Spatial: Σ configured_deltaT = "
               << configured_sum << " ≤ 0");

   const real_t dt_scale = dt_step / configured_sum;
   std::vector<real_t> deltaT_scaled(O);
   for (int o = 0; o < O; ++o)
   {
      deltaT_scaled[o] = configured_deltaT[o] * dt_scale;
   }
   set_substeps(deltaT_scaled, configured_weights);

   // Sub-step MIDPOINT nodes on [0, dt_step].
   std::vector<real_t> tau_nodes(O);
   real_t acc = 0.0;
   for (int o = 0; o < O; ++o)
   {
      tau_nodes[o] = acc + 0.5 * deltaT_scaled[o];
      acc += deltaT_scaled[o];
   }

   std::vector<Vector> Q_per_node;
   wave.ComputeADERSubStepStates(Q, dt_step, ader_order, tau_nodes,
                                 Q_per_node);
   MFEM_VERIFY(static_cast<int>(Q_per_node.size()) == O,
               "AdvanceADERWithSubStep_Spatial: ComputeADERSubStepStates "
               "returned " << Q_per_node.size() << " nodes, expected " << O);

   const int n_total_fault_qps = wave.GetNumTotalFaultQPs();
   std::vector<std::vector<real_t>> Q_pointwise_plus(O), Q_pointwise_minus(O);
   for (int o = 0; o < O; ++o)
   {
      wave.EvaluateBulkAtFaultQPsCanonical(Q_per_node[o],
                                           Q_pointwise_plus[o],
                                           Q_pointwise_minus[o]);
   }

   const size_t n_words =
      static_cast<size_t>(NUM_STATE)
      * static_cast<size_t>(n_total_fault_qps);
   std::vector<real_t> I_imp_plus_flat(n_words, 0.0);
   std::vector<real_t> I_imp_minus_flat(n_words, 0.0);

   if (n_total_fault_qps > 0)
   {
      do_iterate(Q_pointwise_plus, Q_pointwise_minus,
                 dt_step, t_step_start,
                 I_imp_plus_flat.data(),
                 I_imp_minus_flat.data(),
                 nuc_callback);
   }

   struct ImposedGuard
   {
      WaveOperator<ParMesh> &w_;
      explicit ImposedGuard(WaveOperator<ParMesh> &w) : w_(w) {}
      ~ImposedGuard() { w_.ResetSubStepFaultImposedStates(); }
   };
   wave.SetSubStepFaultImposedStates(
      n_total_fault_qps > 0 ? I_imp_plus_flat.data()  : nullptr,
      n_total_fault_qps > 0 ? I_imp_minus_flat.data() : nullptr,
      n_total_fault_qps);
   ImposedGuard guard(wave);

   wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
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

#ifdef SEAS_DIAG_FAULT_FLUX
   mfem::seas::g_seas_my_rank = rank;
#endif

   // -----------------------------------------------------------------
   // 1.  CLI parse — config-driven; CLI overrides are merged after
   //     TOML load (later wins).
   // -----------------------------------------------------------------
   const std::string config_path =
      GetStringArg(argc, argv, "--config", "");
   if (config_path.empty())
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: --config PATH.toml is required.\n"
                   << "  Usage: seas_spatial_dyn_driver --config "
                   << "<file>.toml [OPTIONS]\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 2;
   }

   const bool dry_run         = HasFlag(argc, argv, "--dry-run");
   const bool verify_dispatch = HasFlag(argc, argv, "--verify-dispatch");
   const bool no_sidecar_material =
      HasFlag(argc, argv, "--no-sidecar-material");
   const bool print_derived = HasFlag(argc, argv, "--print-derived");

   const std::string cli_mesh        = GetStringArg(argc, argv, "--mesh", "");
   const std::string cli_vel_model   =
      GetStringArg(argc, argv, "--velocity-model", "");
   const std::string cli_vel_path    =
      GetStringArg(argc, argv, "--override-velocity-path", "");
   const std::string cli_stress_kind =
      GetStringArg(argc, argv, "--stress-kind", "");
   const std::string cli_stress_sidecar =
      GetStringArg(argc, argv, "--stress-sidecar", "");

   const real_t cli_tfinal     = GetRealArg(argc, argv, "--tfinal", -1.0);
   const real_t cli_cfl        = GetRealArg(argc, argv, "--cfl", -1.0);
   const int    cli_ader_order = GetIntArg(argc, argv, "--ader-order", -1);
   const std::string cli_mixed_flux =
      GetStringArg(argc, argv, "--mixed-flux", "");
   const bool   cli_pml        = HasFlag(argc, argv, "--pml");

   const std::string cli_output_dir =
      GetStringArg(argc, argv, "--output-dir", "");
   const std::string cli_pv_volume = GetStringArg(argc, argv, "--paraview-volume", "");
   const std::string cli_pv_bulk   = GetStringArg(argc, argv, "--paraview-bulk", "");
   const std::string cli_pv_fault  = GetStringArg(argc, argv, "--paraview-fault", "");
   const real_t cli_pv_vol_dt   = GetRealArg(argc, argv, "--paraview-volume-dt", -1.0);
   const real_t cli_pv_bulk_dt  = GetRealArg(argc, argv, "--paraview-bulk-dt", -1.0);
   const real_t cli_pv_fault_dt = GetRealArg(argc, argv, "--paraview-fault-dt", -1.0);
   const real_t cli_pv_vol_zfp  = GetRealArg(argc, argv, "--paraview-volume-zfp-tol", -1.0);
   const real_t cli_pv_bulk_zfp = GetRealArg(argc, argv, "--paraview-bulk-zfp-tol", -1.0);
   const real_t cli_pv_fault_zfp= GetRealArg(argc, argv, "--paraview-fault-zfp-tol", -1.0);
   const int    cli_pv_max_snap = GetIntArg(argc, argv, "--paraview-max-snapshots", -1);

   // Parity Phase 1: extended ParaView CLI flag surface.
   const bool   cli_pv                   = HasFlag(argc, argv, "--paraview");
   const real_t cli_pv_dt                = GetRealArg(argc, argv, "--paraview-dt", -1.0);
   const int    cli_pv_every             = GetIntArg (argc, argv, "--paraview-every", 0);
   const bool   cli_pv_force_fault_vtu   = HasFlag(argc, argv, "--paraview-fault-vtu");
   const bool   cli_pv_force_fault_hdf5  = HasFlag(argc, argv, "--paraview-fault-hdf5");
   const bool   cli_pv_legacy_ascii      = HasFlag(argc, argv, "--paraview-fault-legacy-ascii");
   const int    cli_pv_fault_deflate     = GetIntArg (argc, argv, "--paraview-fault-deflate-level", -1);
   const int    cli_pv_bulk_deflate      = GetIntArg (argc, argv, "--paraview-bulk-deflate-level",  -1);
   const bool   cli_pv_force_vol_vtu     = HasFlag(argc, argv, "--paraview-volume-vtu");
   const bool   cli_pv_force_vol_hdf5    = HasFlag(argc, argv, "--paraview-volume-hdf5");
   const int    cli_pv_volume_deflate    = GetIntArg (argc, argv, "--paraview-volume-deflate-level", -1);
   const real_t cli_pv_coseismic_dt      = GetRealArg(argc, argv, "--paraview-coseismic-dt",    -1.0);
   const real_t cli_pv_nucleation_dt     = GetRealArg(argc, argv, "--paraview-nucleation-dt",   -1.0);
   const real_t cli_pv_interseismic_dt   = GetRealArg(argc, argv, "--paraview-interseismic-dt", -1.0);

   const std::string restart_prefix =
      GetStringArg(argc, argv, "--restart", "");
   const int cli_checkpoint_every =
      GetIntArg(argc, argv, "--checkpoint-every", -1);

   // -----------------------------------------------------------------
   // 2.  TOML load — owns all defaults; CLI then overrides.
   // -----------------------------------------------------------------
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

   // Merge CLI overrides into cfg (later wins).
   if (!cli_mesh.empty())            { cfg.mesh.path = cli_mesh; }
   if (!cli_vel_model.empty())
   {
      if (cli_vel_model == "cvmh")
      { cfg.velocity.model = spatial::VelocityModel::CVMH; }
      else if (cli_vel_model == "cvm_s4.26.m01")
      { cfg.velocity.model = spatial::VelocityModel::CVMS_4_26_M01; }
      else if (cli_vel_model == "multiscale_statewise")
      { cfg.velocity.model = spatial::VelocityModel::MultiscaleStatewise; }
      else
      {
         MFEM_ABORT("--velocity-model: unknown value '" << cli_vel_model
                    << "'.  Accepted: cvmh | cvm_s4.26.m01 | "
                    << "multiscale_statewise.");
      }
   }
   if (!cli_vel_path.empty())        { cfg.velocity.override_path = cli_vel_path; }
   if (!cli_stress_kind.empty())
   {
      if (cli_stress_kind == "constant_tensor")
      { cfg.stress.kind = spatial::StressSourceKind::ConstantTensor; }
      else if (cli_stress_kind == "sidecar_hdf5")
      { cfg.stress.kind = spatial::StressSourceKind::SidecarHDF5; }
      else
      {
         MFEM_ABORT("--stress-kind: unknown value '" << cli_stress_kind
                    << "'.  Accepted: constant_tensor | sidecar_hdf5.");
      }
   }
   if (!cli_stress_sidecar.empty())  { cfg.stress.sidecar_path = cli_stress_sidecar; }
   if (cli_tfinal > 0.0)             { cfg.time.tfinal = cli_tfinal; }
   if (cli_cfl > 0.0)                { cfg.numerics.cfl = cli_cfl; }
   if (cli_ader_order > 0)           { cfg.numerics.ader_order = cli_ader_order; }
   if (!cli_mixed_flux.empty())      { cfg.numerics.mixed_flux = cli_mixed_flux; }
   if (cli_pml)                      { cfg.numerics.use_pml = true; }
   if (!cli_output_dir.empty())      { cfg.output.output_dir = cli_output_dir; }
   if (!cli_pv_volume.empty())       { cfg.output.paraview_volume = cli_pv_volume; }
   if (!cli_pv_bulk.empty())         { cfg.output.paraview_bulk   = cli_pv_bulk; }
   if (!cli_pv_fault.empty())        { cfg.output.paraview_fault  = cli_pv_fault; }
   if (cli_pv_vol_dt    > 0.0)       { cfg.output.paraview_volume_dt    = cli_pv_vol_dt; }
   if (cli_pv_bulk_dt   > 0.0)       { cfg.output.paraview_bulk_dt      = cli_pv_bulk_dt; }
   if (cli_pv_vol_zfp   > 0.0)       { cfg.output.paraview_volume_zfp_tol = cli_pv_vol_zfp; }
   if (cli_pv_bulk_zfp  > 0.0)       { cfg.output.paraview_bulk_zfp_tol  = cli_pv_bulk_zfp; }
   if (cli_pv_fault_zfp > 0.0)       { cfg.output.paraview_fault_zfp_tol = cli_pv_fault_zfp; }
   if (cli_pv_max_snap  > 0)         { cfg.output.max_snapshots          = cli_pv_max_snap; }
   if (cli_checkpoint_every > 0)     { cfg.output.checkpoint_every_steps = cli_checkpoint_every; }

   // Parity Phase 1: extended ParaView CLI overrides.  Order matters:
   // apply the umbrella `--paraview-dt` FIRST so that the more-specific
   // `--paraview-fault-dt` below overrides it when both are given (R-007
   // fix: specific flag wins over umbrella).
   if (cli_pv)                            { cfg.output.paraview_enabled = true; }
   if (cli_pv_dt              > 0.0)      { cfg.output.paraview_fault_dt    = cli_pv_dt;
                                            cfg.output.paraview_enabled    = true; }
   if (cli_pv_fault_dt  > 0.0)            { cfg.output.paraview_fault_dt   = cli_pv_fault_dt; }
   if (cli_pv_every           > 0)        { cfg.output.paraview_every_steps = cli_pv_every; }
   if (cli_pv_fault_deflate  >= 0)        { cfg.output.paraview_fault_deflate_level   = cli_pv_fault_deflate; }
   if (cli_pv_bulk_deflate   >= 0)        { cfg.output.paraview_bulk_deflate_level    = cli_pv_bulk_deflate; }
   if (cli_pv_volume_deflate >= 0)        { cfg.output.paraview_volume_deflate_level  = cli_pv_volume_deflate; }
   if (cli_pv_coseismic_dt    > 0.0)      { cfg.output.paraview_coseismic_dt    = cli_pv_coseismic_dt; }
   if (cli_pv_nucleation_dt   > 0.0)      { cfg.output.paraview_nucleation_dt   = cli_pv_nucleation_dt; }
   if (cli_pv_interseismic_dt > 0.0)      { cfg.output.paraview_interseismic_dt = cli_pv_interseismic_dt; }
   if (cli_pv_force_fault_vtu)            { cfg.output.paraview_fault  = "vtu";  }
   if (cli_pv_force_fault_hdf5)           { cfg.output.paraview_fault  = "hdf5"; }
   if (cli_pv_force_vol_vtu)              { cfg.output.paraview_volume = "vtu";  }
   if (cli_pv_force_vol_hdf5)             { cfg.output.paraview_volume = "hdf5"; }
   if (cli_pv_legacy_ascii)               { cfg.output.paraview_fault_legacy_ascii = true;
                                            cfg.output.paraview_fault = "vtu"; }

   // Parity Phase 1 build guards.  Fire BEFORE any pv_out construction
   // so a bad flag combination fails at parse time, not after the run.
   if (cli_pv_force_fault_vtu && cli_pv_force_fault_hdf5)
   {
      MFEM_ABORT("--paraview-fault-vtu and --paraview-fault-hdf5 are "
                 "mutually exclusive.");
   }
   if (cli_pv_legacy_ascii && cli_pv_force_fault_hdf5)
   {
      MFEM_ABORT("--paraview-fault-legacy-ascii implies --paraview-fault-vtu; "
                 "it cannot be combined with --paraview-fault-hdf5.");
   }
   if (cli_pv_force_vol_vtu && cli_pv_force_vol_hdf5)
   {
      MFEM_ABORT("--paraview-volume-vtu and --paraview-volume-hdf5 are "
                 "mutually exclusive.");
   }
   // ZFP + deflate combos are mutually exclusive per collection.
   if (cli_pv_vol_zfp > 0.0 && cli_pv_volume_deflate >= 0)
   {
      MFEM_ABORT("--paraview-volume-zfp-tol and --paraview-volume-deflate-level "
                 "are mutually exclusive.");
   }
   if (cli_pv_bulk_zfp > 0.0 && cli_pv_bulk_deflate >= 0)
   {
      MFEM_ABORT("--paraview-bulk-zfp-tol and --paraview-bulk-deflate-level "
                 "are mutually exclusive.");
   }
   if (cli_pv_fault_zfp > 0.0 && cli_pv_fault_deflate >= 0)
   {
      MFEM_ABORT("--paraview-fault-zfp-tol and --paraview-fault-deflate-level "
                 "are mutually exclusive.");
   }
#ifndef MFEM_USE_HDF5
   if (cli_pv_force_fault_hdf5 || cli_pv_force_vol_hdf5
       || cli_pv_vol_zfp > 0.0 || cli_pv_bulk_zfp > 0.0 || cli_pv_fault_zfp > 0.0
       || cli_pv_fault_deflate >= 0 || cli_pv_bulk_deflate >= 0
       || cli_pv_volume_deflate >= 0
       || cfg.output.paraview_volume == "hdf5"
       || cfg.output.paraview_bulk   == "hdf5"
       || cfg.output.paraview_fault  == "hdf5")
   {
      MFEM_ABORT("--paraview-*-hdf5 / --paraview-*-zfp-tol / "
                 "--paraview-*-deflate-level require an MFEM build with "
                 "MFEM_USE_HDF5=YES.  Re-run with --paraview-*-vtu or "
                 "switch the per-collection mode to \"vtu\" / \"off\".");
   }
#endif
#ifndef MFEM_USE_H5Z_ZFP
   if (cli_pv_vol_zfp > 0.0 || cli_pv_bulk_zfp > 0.0 || cli_pv_fault_zfp > 0.0)
   {
      MFEM_ABORT("--paraview-*-zfp-tol requires an MFEM build with "
                 "MFEM_USE_H5Z_ZFP=YES.  Use --paraview-*-deflate-level "
                 "for lossless compression instead.");
   }
#endif

   MFEM_VERIFY(cfg.law == spatial::FrictionLawKind::SlipWeakening ||
               cfg.law == spatial::FrictionLawKind::RateState,
               "spatial_dyn_driver: cfg.law must be SlipWeakening or "
               "RateState; got " << static_cast<int>(cfg.law));
   const bool is_lsw =
      (cfg.law == spatial::FrictionLawKind::SlipWeakening);
   if (!is_lsw)
   {
      MFEM_VERIFY(cfg.rate_state.has_value(),
                  "spatial_dyn_driver: [meta].law=\"rate_state\" but the "
                  "[friction.rate_state] block is absent in TOML.");
   }
   const bool rs_use_srw = !is_lsw
      && cfg.rate_state->state_evolution
         == spatial::StateEvolutionKind::SlipLawStrongRateWeakening;

   if (rank == 0)
   {
      std::cout << "================================================\n"
                << "seas_spatial_dyn_driver — Phase 4 (rev-3)\n"
                << "================================================\n"
                << "config:           " << config_path << "\n"
                << "mesh:             " << cfg.mesh.path << "\n"
                << "fe order:         " << cfg.mesh.order << "\n"
                << "law:              "
                << (is_lsw ? "slip_weakening"
                           : (rs_use_srw ? "rate_state (slip_law_srw)"
                                         : "rate_state (aging_law)")) << "\n"
                << "stress kind:      "
                << (cfg.stress.kind == spatial::StressSourceKind::ConstantTensor
                    ? "constant_tensor"
                    : cfg.stress.kind ==
                      spatial::StressSourceKind::ConstantTensorWithPatches
                      ? "constant_tensor_with_patches"
                      : cfg.stress.kind ==
                        spatial::StressSourceKind::DepthProportionalToShearModulus
                        ? "depth_proportional"
                        : "sidecar_hdf5") << "\n"
                << "tfinal:           " << cfg.time.tfinal << " s\n"
                << "cfl:              " << cfg.numerics.cfl << "\n"
                << "ader order:       " << cfg.numerics.ader_order << "\n"
                << "mixed flux:       " << cfg.numerics.mixed_flux << "\n"
                << "use pml:          " << (cfg.numerics.use_pml ? "yes" : "no")
                << "\n"
                << "nucleation:       "
                << (cfg.nucleation.enabled
                    ? (std::string(spatial::NucleationKindToString(cfg.nucleation.kind))
                       + " (enabled)")
                    : std::string("DISABLED"))
                << "\n"
                << "no-sidecar mat:   " << (no_sidecar_material ? "yes" : "no")
                << "\n"
                << "dry-run:          " << (dry_run ? "yes" : "no") << "\n"
                << "ranks:            " << nprocs << "\n"
                << "================================================\n";
   }

   // -----------------------------------------------------------------
   // 3.  Restart / output-dir safety check (mirrors tpv104_driver.cpp).
   // -----------------------------------------------------------------
   if (!restart_prefix.empty())
   {
      namespace fs = std::filesystem;
      try
      {
         const fs::path restart_path(restart_prefix);
         fs::path restart_dir_path = restart_path.parent_path();
         if (restart_dir_path.empty()) { restart_dir_path = "."; }

         const fs::path restart_canonical =
            fs::weakly_canonical(restart_dir_path);
         const fs::path output_canonical =
            fs::weakly_canonical(fs::path(cfg.output.output_dir));

         if (restart_canonical == output_canonical)
         {
            if (rank == 0)
            {
               std::cerr << "ERROR: --output-dir (" << cfg.output.output_dir
                         << ") resolves to the SAME directory as the "
                         << "parent of --restart (" << restart_dir_path.string()
                         << ").  Pick a DIFFERENT --output-dir for the "
                         << "restarted run.\n";
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

   // -----------------------------------------------------------------
   // 4.  Load mesh.
   //
   // REVIEW R-006: explicit preflight check.  MFEM's Mesh ctor emits a
   // raw `MFEM_ABORT("could not open Gmsh file ...")` that doesn't tell
   // the user how to regenerate it.  Surface a clearer error first so
   // a fresh checkout against a `.geo`-only mesh directory points the
   // user at the right `gmsh -format msh22` command.  (Per CLAUDE.md
   // "Known limitation — Gmsh .msh format", MFEM requires v2.2 format.)
   if (!std::filesystem::exists(cfg.mesh.path))
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: mesh file '" << cfg.mesh.path
                   << "' does not exist.  Generate it with:\n"
                   << "    gmsh -format msh22 -3 <input>.geo -o "
                   << cfg.mesh.path << "\n"
                   << "  (Gmsh v2.2 — required by MFEM; see CLAUDE.md "
                   << "\"Known limitation — Gmsh .msh format\".)\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Barrier(comm);
      MPI_Finalize();
#endif
      return 4;
   }
   Mesh smesh(cfg.mesh.path.c_str(), 1, 1);
   const int dim = smesh.Dimension();
   MFEM_VERIFY(dim == 3,
               "spatial_dyn_driver: only 3D meshes supported; got dim="
               << dim);

#ifdef MFEM_USE_MPI
   ParMesh pmesh(comm, smesh);
#else
#  error "spatial_dyn_driver requires MFEM_USE_MPI=YES."
#endif
   smesh.Clear();
   pmesh.SetCurvature(cfg.mesh.order);

   // -----------------------------------------------------------------
   // 5.  BoundaryConfig — Phase R.3 step 2: read from cfg.boundary.
   //
   // Defaults in `spatial::BoundarySpec` (fault=101, natural=[102],
   // absorbing=[103,104]) match the SAFS .geo so existing SAFS configs
   // continue to behave unchanged (R.3 step 11 backward-compat).
   // TPV205 + TPV31 configs override these via the new [boundary] TOML
   // block.
   // -----------------------------------------------------------------
   BoundaryConfig bc;
   bc.fault_attr      = cfg.boundary.fault_attr;
   bc.natural_attrs   = std::set<int>(cfg.boundary.natural_attrs.begin(),
                                      cfg.boundary.natural_attrs.end());
   bc.absorbing_attrs = std::set<int>(cfg.boundary.absorbing_attrs.begin(),
                                      cfg.boundary.absorbing_attrs.end());
   MFEM_VERIFY(bc.fault_attr > 0,
               "spatial_dyn_driver: [boundary].fault_attr must be > 0; got "
               << bc.fault_attr);

   // -----------------------------------------------------------------
   // 6.  Material — Phase R.3 step 6 + step 1: dispatch on
   //     cfg.material.kind.  Phase R.2 made the heterogeneous
   //     `WaveOperator(MaterialField, BoundaryConfig)` ctor work for
   //     all three modes at INTERIOR faces; the D-1 abort is removed.
   //
   //     Material kinds:
   //       * constant         — Mode::Constant from [material_constant_fallback]
   //       * depth_profile_1d — Mode::Coefficient from MakeDepthProfile1DMaterial
   //       * sidecar_hdf5     — Mode::GridFunction or Mode::Coefficient
   //                            (legacy SpatialVelocityBundle path)
   //
   //     Legacy CLI `--no-sidecar-material` forces "constant" regardless
   //     of `cfg.material.kind` so existing scripts keep working.
   // -----------------------------------------------------------------
   const real_t mat_lambda = cfg.material_fallback.lambda;
   const real_t mat_mu     = cfg.material_fallback.mu;
   const real_t mat_rho    = cfg.material_fallback.rho;

   std::unique_ptr<spatial::SpatialVelocityBundle> vel_bundle;
   std::unique_ptr<DepthProfile1DMaterial>         depth_profile_wrapper;
   MaterialField material = MaterialField::MakeConstant(mat_lambda,
                                                        mat_mu, mat_rho);

   const bool force_constant = no_sidecar_material
                               || cfg.material.kind == spatial::MaterialKind::Constant;

   if (force_constant)
   {
      if (rank == 0)
      {
         std::cout << "[material] kind=constant from "
                   << "[material_constant_fallback]: "
                   << "lambda=" << mat_lambda
                   << " mu=" << mat_mu
                   << " rho=" << mat_rho << "\n";
      }
   }
   else if (cfg.material.kind == spatial::MaterialKind::DepthProfile1D)
   {
      depth_profile_wrapper = MakeDepthProfile1DMaterial(
         cfg.material.profile_layers, cfg.material.depth_axis);
      material = depth_profile_wrapper->field;
      if (rank == 0)
      {
         std::cout << "[material] kind=depth_profile_1d (axis='"
                   << cfg.material.depth_axis << "', "
                   << cfg.material.profile_layers.size()
                   << " layers)\n";
         for (std::size_t i = 0; i < cfg.material.profile_layers.size(); ++i)
         {
            const auto &L = cfg.material.profile_layers[i];
            std::cout << "  layer " << i
                      << ": depth=[" << L.depth_top_m << ", "
                      << L.depth_bot_m << "] m  "
                      << "vp=" << L.vp_ms << " vs=" << L.vs_ms
                      << " rho=" << L.rho_kgm3
                      << " interp=" << L.interp << "\n";
         }
      }
   }
   else  // SidecarHDF5
   {
      try
      {
         vel_bundle = std::make_unique<spatial::SpatialVelocityBundle>(
            spatial::LoadSpatialVelocityBundle(cfg.velocity, pmesh));
         material = vel_bundle->MakeMaterialField();
         if (rank == 0)
         {
            std::cout << "[material] kind=sidecar_hdf5: "
                      << spatial::ResolveSpatialVelocitySidecarPath(cfg.velocity)
                      << "\n";
         }
      }
      catch (const std::exception &e)
      {
         if (rank == 0)
         {
            std::cerr << "ERROR: failed to load velocity sidecar: "
                      << e.what() << "\n"
                      << "  Re-run with --no-sidecar-material to use the "
                      << "[material_constant_fallback] block.\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 4;
      }
   }

   // -----------------------------------------------------------------
   // 7.  Construct WaveOperator.
   //
   //     Default (`[numerics] interior_flux = "bimaterial"`): use the
   //     heterogeneous (MaterialField, BoundaryConfig) ctor (Phase R.2
   //     / R.4 step 2).  Routes every interior face through
   //     `BimaterialFlux::ApplyPerFaceFlux` and exercises the
   //     bi-material precomputation on Mode::Constant input — drift
   //     from the scalar ctor is ~1e-12 relative (verified by
   //     T-PHASEH-SCALAR-PARITY).
   //
   //     REVIEW R-006 opt-in (`[numerics] interior_flux = "scalar"`):
   //     use the scalar ctor `WaveOperator(mesh, order, λ, μ, ρ, bc)`
   //     for byte parity with the native TPV102/TPV104/TPV205 drivers.
   //     Only valid when material is Mode::Constant (homogeneous).
   // -----------------------------------------------------------------
   const bool use_scalar_ctor = (cfg.numerics.interior_flux == "scalar");
   if (use_scalar_ctor)
   {
      MFEM_VERIFY(material.mode == MaterialField::Mode::Constant,
                  "spatial_dyn_driver: [numerics].interior_flux = \"scalar\" "
                  "requires homogeneous material (Mode::Constant); current "
                  "material kind cannot be reduced to scalar (λ, μ, ρ).  "
                  "Either set material.kind = \"constant\" or remove the "
                  "interior_flux opt-in to use the bimaterial path.  "
                  "(REVIEW R-006)");
      if (rank == 0)
      {
         std::cout << "[wave] interior_flux=scalar: using scalar "
                   << "WaveOperator(λ, μ, ρ) ctor for byte parity with "
                   << "native TPV102/TPV104/TPV205 (REVIEW R-006).\n";
      }
   }
   std::unique_ptr<WaveOperator<ParMesh>> wave_ptr =
      use_scalar_ctor
         ? std::make_unique<WaveOperator<ParMesh>>(
              pmesh, cfg.mesh.order,
              material.lambda_const, material.mu_const, material.rho_const,
              bc)
         : std::make_unique<WaveOperator<ParMesh>>(
              pmesh, cfg.mesh.order, material, bc);
   WaveOperator<ParMesh> &wave = *wave_ptr;

   // Representative material values for downstream consumers that
   // expect a scalar (reflection-time warning, FaultFaceFlux seed
   // impedance, PrintDerivedAndCheck mu_bulk).  For Mode::Constant
   // these ARE the actual constants; for non-constant modes they
   // fall back to the [material_constant_fallback] block which is
   // documented as the per-DOF override gets re-seeded by
   // InitializeFaultDOFs_Spatial (the seed values matter only for
   // initial diagnostics).
   const real_t seed_lambda =
      (material.mode == MaterialField::Mode::Constant)
      ? material.lambda_const : mat_lambda;
   const real_t seed_mu =
      (material.mode == MaterialField::Mode::Constant)
      ? material.mu_const     : mat_mu;
   const real_t seed_rho =
      (material.mode == MaterialField::Mode::Constant)
      ? material.rho_const    : mat_rho;

   // R-107 reflection-time warning: compute min_box_dim / cp_max from
   // mesh bounding box + scalar material.
   {
      const real_t cp = std::sqrt((seed_lambda
                                   + 2.0 * seed_mu)
                                   / seed_rho);
      Vector lo(3), hi(3);
      pmesh.GetBoundingBox(lo, hi, 1);
      real_t min_box_dim_local = std::numeric_limits<real_t>::infinity();
      for (int d = 0; d < 3; ++d)
      {
         min_box_dim_local = std::min(min_box_dim_local, hi(d) - lo(d));
      }
      real_t min_box_dim = min_box_dim_local;
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&min_box_dim_local, &min_box_dim, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
#endif
      const real_t t_reflect = (cp > 0.0) ? min_box_dim / cp : 0.0;
      if (rank == 0 && t_reflect < cfg.time.tfinal && !cfg.numerics.use_pml)
      {
         std::cout << "[spatial_dyn] WARNING: tfinal (" << cfg.time.tfinal
                   << " s) exceeds min_box_dim / cp_max ("
                   << t_reflect << " s).  Reflected waves will "
                   << "contaminate the rupture after this time.  Add "
                   << "--pml or shorten tfinal.\n";
      }
   }

   // Dispatch the wave operator's per-face friction law on `cfg.law`:
   //   - LSW: closed-form linear slip-weakening
   //   - RateState: Brent/Newton-Raphson per-QP via fault_flux EvaluateADER
   wave.SetFaultFrictionLaw(is_lsw ? FaultFrictionLaw::LSW
                                   : FaultFrictionLaw::RateAndState);

   // REVIEW R-007: `SetMixedFluxMode` is intentionally DEFERRED to
   // AFTER `SetFaultFlux` + `SetFaultDOFData` + `SetAbsorbingBackground`
   // (the R-1205 required call order, see `CLAUDE.md` and the native
   // TPV205 driver's wiring at `tpv205_driver.cpp:1606-1633`).  The
   // setter's internal cross-checks (e.g., bc.fault_attr > 0 for
   // Adjacent mode) assume the fault wiring is already in place; calling
   // it earlier would silently no-op for "none" but abort with a
   // misleading message for any non-trivial value.  See call site
   // further below (after SetAbsorbingBackground).

   wave.SetTime(cfg.time.t_initial);

   // -----------------------------------------------------------------
   // 8.  Build per-DOF fault tables (deviation D-2: inline walk).
   // -----------------------------------------------------------------
   const Array<int> &fault_int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &fault_shr_faces = wave.GetFaultSharedFaces();

   // MFEM_USE_MPI is unconditional in this driver (#error at L682
   // requires it), so the comm arg can be passed without the
   // #ifdef-in-argument-list dance (R-607 round-6).
   const int nbf_per_face = ProbeNbfPerFace(pmesh, cfg.mesh.order,
                                            fault_int_faces,
                                            fault_shr_faces, comm);

   const int num_fault_local  = fault_int_faces.Size() * nbf_per_face;
   const int num_shared_fault = fault_shr_faces.Size() * nbf_per_face;
   const int num_fault_total  = num_fault_local + num_shared_fault;

   // REVIEW R-003 — basis-consistency fix: the wave operator builds its
   // OWN `FaultBasis` internally at `dynamic/wave_operator.inl:349-350`
   // with hard-coded ref_normal = (0, -1, 0), up = (0, 0, 1) (the
   // Tandem convention used by `FaultFaceFlux::ComputeTrialTraction`).
   // Previously this driver built a SEPARATE FaultBasis from the TOML's
   // `[fault_geometry] ref_normal / up`, fed it to
   // `FaultGeometry::ComputeParams<StressSource>`, and produced
   // tau_pre in a basis that could disagree with the wave operator's
   // internal one.  Specifically, since `strike = up × n_ref`, flipping
   // the sign of ref_normal flips the strike basis vector — so a TOML
   // with ref_normal = (0, +1, 0) (the SAFS schema default!) produced
   // tau_pre.strike with the OPPOSITE sign of the runtime trial
   // traction, sending rupture backwards.
   //
   // Fix: reuse `wave.GetFaultBasis()` directly so external (pre-stress
   // projection) and internal (trial traction) live in the SAME frame
   // by construction.  The TOML's `[fault_geometry].ref_normal / up`
   // fields are now used ONLY for an early-startup check that the user
   // is not inadvertently asking for a frame that disagrees with the
   // wave operator's hard-coded one — this catches the SAFS-default-
   // (0,+1,0) foot-gun before any physics runs.  TPV31's (0,0,-1)
   // ref_normal also trips this check and surfaces the requirement
   // that the wave operator gain a parametric ref_normal before TPV31
   // is supported by this driver.
   Vector ref_normal(3);
   ref_normal(0) = cfg.fault_geometry.ref_normal[0];
   ref_normal(1) = cfg.fault_geometry.ref_normal[1];
   ref_normal(2) = cfg.fault_geometry.ref_normal[2];
   Vector up_vec(3);
   up_vec(0) = cfg.fault_geometry.up[0];
   up_vec(1) = cfg.fault_geometry.up[1];
   up_vec(2) = cfg.fault_geometry.up[2];

   {
      constexpr real_t kEps = 1e-12;
      const bool n_matches =
         std::abs(ref_normal(0) -  0.0) < kEps &&
         std::abs(ref_normal(1) - -1.0) < kEps &&
         std::abs(ref_normal(2) -  0.0) < kEps;
      const bool up_matches =
         std::abs(up_vec(0) - 0.0) < kEps &&
         std::abs(up_vec(1) - 0.0) < kEps &&
         std::abs(up_vec(2) - 1.0) < kEps;
      MFEM_VERIFY(n_matches && up_matches,
                  "spatial_dyn_driver: [fault_geometry].ref_normal / up = ("
                  << ref_normal(0) << "," << ref_normal(1) << ","
                  << ref_normal(2) << ") / (" << up_vec(0) << ","
                  << up_vec(1) << "," << up_vec(2) << ") does NOT match "
                  "the wave operator's internal Tandem convention (0,-1,0) "
                  "/ (0,0,1) at dynamic/wave_operator.inl:349-350.  Set "
                  "[fault_geometry] ref_normal = [0.0, -1.0, 0.0] and "
                  "up = [0.0, 0.0, 1.0] in the TOML so the external "
                  "Cauchy-projection basis matches the runtime trial-"
                  "traction basis.  (Drivers running with the SAFS schema "
                  "default (0,+1,0) need to update — see REVIEW R-003.)");
   }

   // GetFaultBasis() returns nullptr on ranks whose partition contains
   // zero local fault boundary elements (METIS routinely produces this
   // on 200+-rank runs of TPV205-scale meshes — the fault occupies a
   // thin slab around y=0, so most ranks own none of it).  All other
   // drivers (tpv102/tpv104/tpv205) treat the result as a pointer and
   // guard with `if (fb)`; the lone unguarded dereference here was
   // segfaulting at offset 0x4 (FaultBasis::num_faces_) on those ranks.
   // Fall back to a default-constructed empty FaultBasis so the verify
   // and BuildPerDOFFaultTables (which both no-op on empty inputs) see
   // consistent zero-sized state instead of UB.
   const FaultBasis *fbasis_ptr = wave.GetFaultBasis();
   static const FaultBasis empty_fault_basis_{};
   const FaultBasis &fbasis = fbasis_ptr ? *fbasis_ptr : empty_fault_basis_;
   MFEM_VERIFY(fbasis.NumFaces() >=
               fault_int_faces.Size() + fault_shr_faces.Size(),
               "spatial_dyn_driver: wave.GetFaultBasis() has "
               << fbasis.NumFaces() << " faces but driver expects "
               << (fault_int_faces.Size() + fault_shr_faces.Size())
               << " (interior + shared); the wave operator's internal "
               "FaultBasis is missing AppendSharedFaces?");

   std::vector<Vector>          fault_coords;
   Vector                       dof_coords_3d;
   DenseMatrix                  dof_basis;
   Array<int>                   dof_to_elem;
   Array<int>                   dof_to_attr;
   std::vector<IntegrationPoint> dof_ips;
   BuildPerDOFFaultTables(pmesh, cfg.mesh.order, bc.fault_attr,
                          fbasis, fault_int_faces, fault_shr_faces,
                          nbf_per_face,
                          fault_coords, dof_coords_3d, dof_basis,
                          dof_to_elem, dof_to_attr, dof_ips);

   int num_fault_global = num_fault_total;
   int num_shared_global = num_shared_fault;
#ifdef MFEM_USE_MPI
   MPI_Allreduce(&num_fault_total,  &num_fault_global,  1, MPI_INT,
                 MPI_SUM, comm);
   MPI_Allreduce(&num_shared_fault, &num_shared_global, 1, MPI_INT,
                 MPI_SUM, comm);
#endif
   if (rank == 0)
   {
      std::cout << "[fault] QPs per face = " << nbf_per_face
                << ", num_fault_global = " << num_fault_global
                << " (local = " << num_fault_local
                << ", shared = " << num_shared_fault << ")\n";
   }
   // REVIEW R-005: hard-abort if the configured fault_attr matches zero
   // mesh faces.  Mirrors the native TPV205/TPV102/TPV104 drivers
   // (drivers/tpv205_driver.cpp:1216-1218).  Without this, a typo'd
   // [boundary].fault_attr produces a silent useless run that completes
   // with no fault DOFs, no station traces, and no ParaView fault file.
   MFEM_VERIFY(num_fault_global > 0,
               "spatial_dyn_driver: no fault faces with attr="
               << bc.fault_attr << " found in mesh '"
               << cfg.mesh.path << "'.  Check [boundary].fault_attr in "
               "the TOML and the mesh's Physical Surface tags.");

   // -----------------------------------------------------------------
   // 9.  FaultGeometry via the new BP5 ctor (Phase 5a).  Seed BP5Params
   //     with its in-class defaults so ComputeBP5Params (skipped in
   //     the new ctor) doesn't divide by zero anywhere.
   // -----------------------------------------------------------------
   BP5Params bp5_seed;
   // Use the (argc, argv) MPIContext ctor (matches the seas_driver.cpp:190
   // reference pattern).  That ctor sets `comm_ = MPI_COMM_WORLD` directly
   // and leaves `owns_comm_ = false`, so the destructor is a no-op.
   // Avoids the `MPI_Comm_free after MPI_FINALIZE was invoked` abort that
   // the comm-dup ctor would trigger when `mpi_ctx`'s destructor fires
   // on the way out of main, AFTER `MPI_Finalize` has returned.
   // `Mpi::Init` inside the ctor is guarded by `Mpi::IsInitialized()`, so
   // the raw `MPI_Init` we already called above is the actual initialiser
   // and the ctor only re-runs `Hypre::Init()` (idempotent).
   MPIContext mpi_ctx(&argc, &argv);
   FaultGeometry<ParMesh> geom(bp5_seed,
                               dof_coords_3d, dof_basis, dof_to_elem,
                               nbf_per_face, &mpi_ctx, dof_ips);

   // -----------------------------------------------------------------
   // 10. Apply stress source (Phase 3 or 3b).
   // -----------------------------------------------------------------
   if (cfg.stress.kind == spatial::StressSourceKind::ConstantTensor)
   {
      spatial::ConstantTensorStressSource src(cfg.stress.sigma_xx_pa,
                                              cfg.stress.sigma_yy_pa,
                                              cfg.stress.sigma_zz_pa,
                                              cfg.stress.sigma_xy_pa,
                                              cfg.stress.sigma_yz_pa,
                                              cfg.stress.sigma_xz_pa);
      geom.ComputeParams(src,
                             cfg.stress.pore_pressure.P_p_pa,
                             cfg.stress.pore_pressure.P_p_grad_pa_per_m,
                             cfg.stress.pore_pressure.min_sigma_n_pa);
   }
   else if (cfg.stress.kind ==
            spatial::StressSourceKind::ConstantTensorWithPatches)
   {
      // TPV205-style: background constant tensor + N static square
      // patches that override per-component shear at t = 0 (NOT a
      // time-dependent perturbation — the patches are baked into
      // tau_pre_ via ComputeParams<StressSource>).
      spatial::ConstantTensorWithPatchesStressSource src(
         cfg.stress.sigma_xx_pa, cfg.stress.sigma_yy_pa,
         cfg.stress.sigma_zz_pa, cfg.stress.sigma_xy_pa,
         cfg.stress.sigma_yz_pa, cfg.stress.sigma_xz_pa,
         cfg.stress.patches);
      geom.ComputeParams(src,
                             cfg.stress.pore_pressure.P_p_pa,
                             cfg.stress.pore_pressure.P_p_grad_pa_per_m,
                             cfg.stress.pore_pressure.min_sigma_n_pa);
   }
   else if (cfg.stress.kind ==
            spatial::StressSourceKind::DepthProportionalToShearModulus)
   {
      // REVIEW R-003: TPV31-style depth-proportional pre-stress.
      //   σ(x, y, z) = sigma_*_per_mu · μ(x, y, z) / μ_ref
      // The `mu_at_xyz` callback comes from the current MaterialField.
      // Today the only material kind producing coordinate-only μ is
      // `depth_profile_1d` (via `DepthProfile1DMaterial::eval_at_xyz`).
      // For `constant` material, μ is a fixed scalar; for `sidecar_hdf5`
      // a coordinate-only lookup would require a point locator we have
      // not yet built — abort with a clear message in that case.
      spatial::DepthProportionalToShearModulusStressSource::MuAtFn mu_at_xyz;
      if (material.mode == MaterialField::Mode::Constant)
      {
         const real_t mu_const = material.mu_const;
         mu_at_xyz = [mu_const](real_t /*x*/, real_t /*y*/, real_t /*z*/)
                     { return mu_const; };
      }
      else if (depth_profile_wrapper != nullptr)
      {
         MFEM_VERIFY(static_cast<bool>(depth_profile_wrapper->eval_at_xyz),
                     "spatial_dyn_driver: depth-profile material has no "
                     "eval_at_xyz callback (heterogeneous_material.cpp "
                     "wiring lost?).");
         auto& eval = depth_profile_wrapper->eval_at_xyz;
         mu_at_xyz = [&eval](real_t x, real_t y, real_t z) -> real_t {
            real_t lam, mu, rho;
            eval(x, y, z, lam, mu, rho);
            return mu;
         };
      }
      else
      {
         MFEM_ABORT("spatial_dyn_driver: [stress] kind = "
                    "\"depth_proportional\" requires either [material] "
                    "kind = \"constant\" or \"depth_profile_1d\".  "
                    "Coordinate-only μ lookup for sidecar_hdf5 is not "
                    "yet implemented.");
      }
      const auto& dp = cfg.stress.depth_proportional;
      spatial::DepthProportionalToShearModulusStressSource src(
         dp.sigma_xx_per_mu, dp.sigma_yy_per_mu, dp.sigma_zz_per_mu,
         dp.sigma_xy_per_mu, dp.sigma_yz_per_mu, dp.sigma_xz_per_mu,
         dp.mu_ref_pa, std::move(mu_at_xyz));
      geom.ComputeParams(src,
                             cfg.stress.pore_pressure.P_p_pa,
                             cfg.stress.pore_pressure.P_p_grad_pa_per_m,
                             cfg.stress.pore_pressure.min_sigma_n_pa);
   }
   else
   {
      MFEM_VERIFY(cfg.stress.kind == spatial::StressSourceKind::SidecarHDF5,
                  "spatial_dyn_driver: unhandled stress kind "
                  << static_cast<int>(cfg.stress.kind)
                  << ".  Valid kinds: constant_tensor, "
                  "constant_tensor_with_patches, depth_proportional, "
                  "sidecar_hdf5.");
      spatial::ApplyCsmStressSidecar(cfg.stress, geom);
   }
   MFEM_VERIFY(geom.HasParams(),
               "spatial_dyn_driver: stress source projection failed");

   // -----------------------------------------------------------------
   // 11. Resolve per-DOF friction parameters (LSW or RateState).
   // -----------------------------------------------------------------
   spatial::SpatialFrictionResolver resolver;
   spatial::SlipWeakeningPerDOFParams lsw;
   spatial::RateStatePerDOFParams    rs;
   if (is_lsw)
   {
      MFEM_VERIFY(cfg.slip_weakening.has_value(),
                  "spatial_dyn_driver: [meta].law=slip_weakening but the "
                  "[friction.slip_weakening] block is absent in TOML.");
      lsw = resolver.ResolveSlipWeakening(*cfg.slip_weakening,
                                          dof_coords_3d, dof_to_attr);
   }
   else
   {
      rs = resolver.ResolveRateState(*cfg.rate_state,
                                     dof_coords_3d, dof_to_elem,
                                     dof_to_attr, material, pmesh,
                                     cfg.stress.pore_pressure,
                                     geom.sigma_n_per_dof());
   }

   // -----------------------------------------------------------------
   // 12. Phase N: resolve the active nucleation kind.  Four kinds:
   //       `gradual_overstress`              (Gaussian, smoothStep ramp)
   //       `square_overstress`               (rectangular, smoothStep ramp)
   //       `instantaneous_overstress_circular` (TPV31: circular cosine
   //                                            taper, INSTANTANEOUS,
   //                                            per-DOF µ-scaling)
   //       `gradual_overstress_compact_circular`
   //                                         (SCEC TPV101/102/104 spec:
   //                                          F(r)=exp(r²/(r²-R²)) +
   //                                          smoothStep ramp)
   //     Exactly one resolver returns non-empty per-DOF amplitudes; the
   //     others return zero-sized Vectors.  When `enabled == false`, ALL
   //     return zero-sized.
   // -----------------------------------------------------------------
   const bool gradual_active = cfg.nucleation.enabled
      && cfg.nucleation.kind == spatial::NucleationKind::GradualOverstress;
   const bool square_active  = cfg.nucleation.enabled
      && cfg.nucleation.kind == spatial::NucleationKind::SquareOverstress;
   const bool instant_active = cfg.nucleation.enabled
      && cfg.nucleation.kind
         == spatial::NucleationKind::InstantaneousOverstressCircular;
   const bool compact_circular_active = cfg.nucleation.enabled
      && cfg.nucleation.kind
         == spatial::NucleationKind::GradualOverstressCompactCircular;

   const spatial::GradualOverstressPerDOFParams nuc_params =
      spatial::ResolveGradualOverstress(
         cfg.nucleation.gradual_overstress,
         gradual_active,
         dof_coords_3d,
         dof_basis);

   const spatial::SquareOverstressPerDOFParams sq_nuc_params =
      spatial::ResolveSquareOverstress(
         cfg.nucleation.square_overstress,
         square_active,
         dof_coords_3d);

   // Per-DOF µ for the instantaneous_overstress_circular µ-scaling.
   // Empty (size 0) when this kind is inactive so the resolver
   // early-returns.
   Vector mu_per_fault_dof;
   if (instant_active && num_fault_total > 0)
   {
      mu_per_fault_dof.SetSize(num_fault_total);
      for (int i = 0; i < num_fault_total; ++i)
      {
         const int elem = dof_to_elem[i];
         mfem::IsoparametricTransformation Tr;
         pmesh.GetElementTransformation(elem, &Tr);
         Tr.SetIntPoint(&dof_ips[i]);
         real_t lam_i, mu_i, rho_i;
         material.EvalAt(elem, Tr, dof_ips[i], lam_i, mu_i, rho_i);
         mu_per_fault_dof(i) = mu_i;
      }
   }
   const spatial::InstantaneousOverstressCircularPerDOFParams ic_nuc_params =
      spatial::ResolveInstantaneousOverstressCircular(
         cfg.nucleation.instantaneous_overstress_circular,
         instant_active,
         dof_coords_3d,
         mu_per_fault_dof);

   const spatial::GradualOverstressCompactCircularPerDOFParams cc_nuc_params =
      spatial::ResolveGradualOverstressCompactCircular(
         cfg.nucleation.gradual_overstress_compact_circular,
         compact_circular_active,
         dof_coords_3d,
         dof_basis);

   // R-008: warn when a non-trivial fraction of fault DOFs live on
   // shared faces.  The wave operator's shared-face EvaluateADER_LSW
   // call reads DOFData::tau{1,2}_nuc AFTER the Phase N per-substep
   // iterator has accumulated the full smoothStep increment for the
   // macrostep, so those DOFs see the perturbation as an end-of-
   // macrostep step rather than a smooth ramp (1st-order time-
   // accuracy degradation).  Quantify and warn so the user can
   // tighten dt or accept the trade-off.
   if (cfg.nucleation.enabled && num_shared_global > 0 && rank == 0)
   {
      const real_t shared_frac = (num_fault_global > 0)
         ? (static_cast<real_t>(num_shared_global)
            / static_cast<real_t>(num_fault_global))
         : 0.0;
      std::cout << "[spatial_dyn] WARNING: " << num_shared_global
                << " of " << num_fault_global << " fault DOFs ("
                << (100.0 * shared_frac) << "%) live on shared faces "
                << "and will see the gradual_overstress perturbation as "
                << "an end-of-macrostep step rather than a smooth ramp.  "
                << "Tighten dt (smaller macrostep) to reduce the "
                << "1st-order time-accuracy error on shared faces.\n";
   }

   // -----------------------------------------------------------------
   // 13. Construct FaultFaceFlux with seed scalar impedances; the
   //     per-DOF impedances are overwritten by InitializeFaultDOFs_Spatial.
   // -----------------------------------------------------------------
   const real_t cp_seed = std::sqrt((seed_lambda
                                     + 2.0 * seed_mu)
                                     / seed_rho);
   const real_t cs_seed = std::sqrt(seed_mu / seed_rho);
   FaultFaceFlux fault_flux(seed_rho, cp_seed, cs_seed);
   wave.SetFaultFlux(&fault_flux);

   // -----------------------------------------------------------------
   // 14. Initialise per-DOF DOFData via the new free function
   //     (Phase 5c, IP-aware overload).
   //
   //     R-003 (Phase N): InitializeFaultDOFs_Spatial requires
   //     T_forced_s / t0_decay_s vectors of size num_fault_total (size-
   //     validated at spatial_setup.hpp:194-200).  Under the LSW
   //     dispatch (not LSW_ForcedRupture) these are NEVER read — supply
   //     the "never forced" sentinel so the size validator passes.  The
   //     dummy values are inert by construction.
   // -----------------------------------------------------------------
   Vector dummy_T_forced(num_fault_total);  dummy_T_forced = 1.0e9;
   Vector dummy_t0_decay(num_fault_total);  dummy_t0_decay = 0.0;

   std::vector<DOFData> dof_data;
   // Per-DOF V_w side-channel (FVW / SCEC FL=103 only).  Sized to
   // num_fault_total when state_evolution == slip_law_srw; left empty
   // otherwise.  Consumed by `Tpv104SubStepIterator::AdvanceWith...`.
   std::vector<real_t> Vw_per_dof;
   if (num_fault_total > 0)
   {
      if (is_lsw)
      {
         spatial::InitializeFaultDOFs_Spatial<ParMesh>(
            dof_data, num_fault_total, dof_to_elem, material, pmesh,
            lsw, geom.GetTauPre(), geom.sigma_n_per_dof(),
            dummy_T_forced, dummy_t0_decay,
            dof_ips);
      }
      else
      {
         // REVIEW R-008: IP-aware overload so per-DOF impedances are
         // evaluated at the actual fault QP, not at the bulk element
         // centroid.  Required for any RS config with non-constant
         // material (e.g., depth_profile_1d); harmless and consistent
         // for constant material.
         spatial::InitializeFaultDOFs_Spatial_RS<ParMesh>(
            dof_data, num_fault_total, dof_to_elem, material, pmesh,
            rs, geom.GetTauPre(), geom.sigma_n_per_dof(), dof_ips);

         // Override the at-rest seeding for rate-state:
         //   - psi from steady-state inversion using the friction
         //     solver's own DOFData inputs:
         //         tau = sigma_n0 * a * asinh((V/(2 V_0)) exp(psi/a))
         //     ⇒  psi = a * ln( (2 V_0 / V) sinh(tau/(sigma_n0 a)) )
         //     evaluated via the numerically stable
         //         log(x · sinh(c)) = |c| + log((x/2)·-sign(c)·expm1(-2|c|))
         //     identity (TPV102/TPV104 `ComputeInitialPsi*` helper).
         //   - slip_rate / V{1,2} from V_init magnitude and the pre-stress
         //     direction (V parallel to tau_pre, CLAUDE.md "Slip rate
         //     direction").
         //
         // R-006 (tpv102_tpv104_review.md): the friction-equation inputs
         // `a`, `sigma_n0`, `tau1_0`, `tau2_0` are read from `dof_data[i]`
         // (the values the runtime solver actually consults), NOT from
         // the resolver/geometry side-arrays.  `V_0` and `V_init` stay
         // on the resolver side because they have no DOFData mirror
         // (the spatial RS init does not duplicate them).  This keeps
         // the at-rest equilibrium `|tau| = sigma_n0 * a * asinh(...)`
         // self-consistent against whatever the underlying
         // `InitializeFaultDOFs_Spatial_RS` populated.
         for (int i = 0; i < num_fault_total; ++i)
         {
            DOFData &d = dof_data[i];
            const real_t a_i      = d.a;
            const real_t V0_i     = rs.V_0(i);
            const real_t V_init_i = std::max(rs.V_init(i),
                                             static_cast<real_t>(1.0e-300));
            const real_t sn_i     = std::abs(d.sigma_n0);
            const real_t tau1_pre = d.tau1_0;
            const real_t tau2_pre = d.tau2_0;
            const real_t tau_abs  = std::sqrt(tau1_pre * tau1_pre
                                              + tau2_pre * tau2_pre);

            // ψ from the stable-asinh inversion of the RS friction
            // coefficient (TPV102/TPV104 `ComputeInitialPsi*`).
            const real_t arg_c = tau_abs / (sn_i * a_i);
            const real_t x     = 2.0 * V0_i / V_init_i;
            const real_t sign_c =
               (arg_c >= 0.0) ? static_cast<real_t>(1.0)
                              : static_cast<real_t>(-1.0);
            const real_t absC  = std::abs(arg_c);
            d.psi = a_i * (absC + std::log(x / 2.0 * -sign_c
                                           * std::expm1(-2.0 * absC)));

            // Slip rate direction PARALLEL to tau_pre (R-801 BP5 frame:
            // tangent1 = dip, tangent2 = strike).  CLAUDE.md "Slip rate
            // direction": V_vec = (V_abs / tau_abs) · tau_vec.
            d.slip_rate = V_init_i;
            if (tau_abs > 0.0)
            {
               d.V1 = (V_init_i / tau_abs) * tau1_pre;
               d.V2 = (V_init_i / tau_abs) * tau2_pre;
            }
            else
            {
               d.V1 = 0.0;
               d.V2 = V_init_i;   // arbitrary tangential direction
            }
            d.slip1 = 0.0;
            d.slip2 = 0.0;
         }

         if (rs_use_srw)
         {
            MFEM_VERIFY(rs.V_w.Size() == num_fault_total,
                        "spatial_dyn_driver: rate_state.V_w size mismatch "
                        "(expected " << num_fault_total << ", got "
                        << rs.V_w.Size() << ")");
            Vw_per_dof.resize(num_fault_total);
            for (int i = 0; i < num_fault_total; ++i)
            { Vw_per_dof[i] = rs.V_w(i); }
         }
      }
   }

   // Instantaneous overstress (TPV31): write the full per-DOF Δτ into
   // DOFData::tau{1,2}_nuc ONCE, right after init.  No per-sub-step
   // ramp — the spec is `t = 0+`.  This is a no-op when the resolver
   // returned zero-sized Vectors (i.e. instant_active == false).
   spatial::ApplyInstantaneousOverstressCircular(dof_data, ic_nuc_params);

   wave.SetFaultDOFData(&dof_data, nbf_per_face);

   // Fluctuation-Q dispatch (matches TPV205): Q_bg = 0.
   {
      real_t Q_bg[NUM_STATE] = {0};
      wave.SetAbsorbingBackground(Q_bg);
   }

   // REVIEW R-007: SetMixedFluxMode at the R-1205-mandated point in the
   // setter sequence — AFTER ctor + SetFaultFlux + SetFaultDOFData +
   // SetAbsorbingBackground.  Mirrors `tpv205_driver.cpp:1606-1633`.
   wave.SetMixedFluxMode(ParseMixedFlux(cfg.numerics.mixed_flux));

   // -----------------------------------------------------------------
   // 15. CFL / Δt and derived numbers.  ComputeMaxDt is the scalar-
   //     material implementation (deviation D-1).
   //
   // REVIEW R-004: apply the DG safety factor 1 / (3*(2p+1)) when the
   // TOML opts in via `[numerics] cfl_safety = "dg"`.  This matches the
   // native TPV102/TPV104/TPV205 drivers (`cfl_factor / (3*(2*order+1))`)
   // and is required for byte parity with their outputs.  Default
   // `cfl_safety = "raw"` keeps existing SAFS configs unchanged.
   // -----------------------------------------------------------------
   real_t cfl_for_max_dt = cfg.numerics.cfl;
   if (cfg.numerics.cfl_safety == "dg")
   {
      cfl_for_max_dt /= (3.0 * (2.0 * cfg.mesh.order + 1.0));
      if (rank == 0)
      {
         std::cout << "[time] cfl_safety=dg: scaled cfl from "
                   << cfg.numerics.cfl << " to " << cfl_for_max_dt
                   << " (= cfl / (3*(2p+1)) with p=" << cfg.mesh.order
                   << ") to match native TPV/BP5 driver convention.\n";
      }
   }
   const real_t dt_cfl = wave.ComputeMaxDt(cfl_for_max_dt);
   real_t dt = (cfg.time.dt_initial > 0.0)
                ? cfg.time.dt_initial : dt_cfl;
   if (cfg.time.dt_max > 0.0 && dt > cfg.time.dt_max)
   {
      if (rank == 0)
      {
         std::cout << "[time] tightened by user dt_max: dt_cfl = "
                   << dt_cfl << " s -> dt = " << cfg.time.dt_max << " s\n";
      }
      dt = cfg.time.dt_max;
   }
   const int nsteps = (cfg.time.tfinal > 0.0)
                       ? static_cast<int>(std::ceil(cfg.time.tfinal / dt)) : 0;
   if (rank == 0)
   {
      std::cout << "[time] dt_cfl = " << dt_cfl << " s\n"
                << "[time] dt     = " << dt     << " s\n"
                << "[time] nsteps = " << nsteps << "\n";
   }

   if (print_derived && is_lsw)
   {
      // Compute mesh h_min for the L_nuc / h_min ratio.
      real_t h_min_local = std::numeric_limits<real_t>::infinity();
      for (int e = 0; e < pmesh.GetNE(); ++e)
      {
         h_min_local = std::min(h_min_local, pmesh.GetElementSize(e, 0));
      }
      real_t h_min_global = h_min_local;
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&h_min_local, &h_min_global, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
#endif
      spatial::PrintDerivedConfig pd_cfg;
      pd_cfg.enabled                = true;
      pd_cfg.abort_on_failure       = true;
      pd_cfg.outside_safety_factor  = 3.0;
      (void)spatial::PrintDerivedAndCheck(
         pd_cfg, lsw,
         geom.GetTauPre(), geom.sigma_n_per_dof(),
         dof_coords_3d,
         cfg.nucleation, nuc_params,
         cfg.stress,
         /*mu_bulk=*/seed_mu,
         cp_seed, cs_seed,
         h_min_global,
         dt_cfl, cfg.time.tfinal,
         num_fault_global,
         geom.NumZeroNormalFallbacks()
#ifdef MFEM_USE_MPI
         , comm
#endif
         , rank);
   }

   // -----------------------------------------------------------------
   // 16. Dry-run exit.
   // -----------------------------------------------------------------
   if (dry_run)
   {
      if (rank == 0)
      {
         std::cout << "[spatial_dyn] --dry-run: construction complete, "
                   << "exiting.\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 0;
   }

   // -----------------------------------------------------------------
   // 17. ParaView output wiring (parity Phases 1-6).
   //
   //     Parity Phase 3 master gate: if paraview_enabled == false AND
   //     every per-collection mode is "off", construct NOTHING.  The
   //     CLI `--paraview` / `--paraview-dt` overrides flip paraview_
   //     enabled on automatically (see merge block above).
   // -----------------------------------------------------------------
   // Parity Phase 4/5: wave-state component ordering assumptions.
   static_assert(SXX == 0 && SYY == 1 && SZZ == 2
                 && SXY == 3 && SYZ == 4 && SXZ == 5,
                 "Phase 4 sigma memcpy assumes wave_state.hpp ordering.");
   static_assert(VX == 6 && VY == 7 && VZ == 8,
                 "Phase 5 velocity memcpy assumes contiguous VX/VY/VZ; "
                 "wave_state.hpp ordering changed.");

   if (rank == 0)
   {
      std::filesystem::create_directories(cfg.output.output_dir);
   }
#ifdef MFEM_USE_MPI
   MPI_Barrier(comm);
#endif

   bool volume_pv_enabled = true;
   bool fault_pv_enabled  = true;
   bool bulk_pv_enabled   = true;
   auto volume_mode = ParseVolumeMode(cfg.output.paraview_volume,
                                      volume_pv_enabled);
   auto fault_mode  = ParseFaultMode(cfg.output.paraview_fault,
                                     fault_pv_enabled);
   {
      bool dummy;
      (void)ParseVolumeMode(cfg.output.paraview_bulk, dummy);
      bulk_pv_enabled = dummy;
   }

   // Parity Phase 1 §4 / Phase 3 §1: paraview_enabled is the HARD master
   // gate.  When false, force every per-collection enabled flag to false
   // regardless of the per-collection mode strings.  The CLI `--paraview`
   // merge above flips paraview_enabled to true when the user opts in.
   if (!cfg.output.paraview_enabled)
   {
      volume_pv_enabled = false;
      fault_pv_enabled  = false;
      bulk_pv_enabled   = false;
   }

   const bool any_pv_requested = volume_pv_enabled
                                || fault_pv_enabled
                                || bulk_pv_enabled;
   const bool primary_pv_active = any_pv_requested
                                  && (volume_pv_enabled || fault_pv_enabled);
   const bool bulk_pv_active    = any_pv_requested && bulk_pv_enabled;

   std::unique_ptr<seas::ParaViewOutput<ParMesh>> pv_out;
   std::unique_ptr<seas::ParaViewOutput<ParMesh>> pv_bulk_out;

   // Parity Phase 5: volume velocity + mpi_rank fields (lifetime of pv_out).
   std::unique_ptr<L2_FECollection> pv_vel_fec, pv_rank_fec;
   std::unique_ptr<ParFiniteElementSpace> pv_vel_fes, pv_rank_fes;
   std::unique_ptr<ParGridFunction>       pv_vel_gf, pv_rank_gf;

   // Parity Phase 4: secondary collection stress GFs (lifetime of pv_bulk_out).
   std::unique_ptr<L2_FECollection>       pv_bulk_sigma_fec;
   std::unique_ptr<ParFiniteElementSpace> pv_bulk_sigma_fes;
   std::unique_ptr<ParGridFunction>
      pv_bulk_sxx_gf, pv_bulk_syy_gf, pv_bulk_szz_gf,
      pv_bulk_sxy_gf, pv_bulk_syz_gf, pv_bulk_sxz_gf;

   // Parity Phase 6: SAFS fault static-parameter arrays (lifetime of pv_out).
   Vector pv_lsw_mu_s, pv_lsw_mu_d, pv_lsw_d_c;
   Vector pv_nuc_amplitude, pv_nuc_radial;
   Vector pv_sig_n_init, pv_tau1_init, pv_tau2_init;

   if (primary_pv_active)
   {
      pv_out = std::make_unique<seas::ParaViewOutput<ParMesh>>(
                  cfg.output.output_dir, pmesh, cfg.mesh.order,
                  "volume", volume_mode);
      pv_out->SetVolumePVDt(cfg.output.paraview_volume_dt);
      pv_out->SetVolumeSaveEnabled(volume_pv_enabled);
      if (cfg.output.paraview_every_steps > 0)
      {
         pv_out->output_every_n_steps = cfg.output.paraview_every_steps;
      }

      pv_out->fixed_dt = cfg.output.paraview_fault_dt;
      pv_out->GetSchedule().max_total_snapshots = cfg.output.max_snapshots;

      // Parity Phase 6: regime-adaptive cadence from TOML/CLI.
      pv_out->SetTotalRunTime(cfg.time.tfinal);
      if (cfg.output.paraview_coseismic_dt    > 0.0)
      { pv_out->GetSchedule().dt_coseismic    = cfg.output.paraview_coseismic_dt; }
      if (cfg.output.paraview_nucleation_dt   > 0.0)
      { pv_out->GetSchedule().dt_nucleation   = cfg.output.paraview_nucleation_dt; }
      if (cfg.output.paraview_interseismic_dt > 0.0)
      { pv_out->GetSchedule().dt_interseismic = cfg.output.paraview_interseismic_dt; }

#ifdef MFEM_USE_HDF5
      if (cfg.output.paraview_volume_zfp_tol > 0.0
          && volume_mode == ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5)
      {
         pv_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            cfg.output.paraview_volume_zfp_tol);
      }
      else if (cfg.output.paraview_volume_deflate_level >= 0
               && volume_mode == ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5)
      {
         pv_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::Deflate,
            static_cast<double>(cfg.output.paraview_volume_deflate_level));
      }
#endif

      pv_out->SetFaultOutputMode(fault_mode);
#ifdef MFEM_USE_HDF5
      if (fault_pv_enabled
          && fault_mode == ParaViewOutput<ParMesh>::FaultOutputMode::Hdf5
          && cfg.output.paraview_fault_zfp_tol > 0.0)
      {
         pv_out->SetFaultHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            cfg.output.paraview_fault_zfp_tol);
      }
#endif

      pv_out->InitFaultOutputBP5(fault_int_faces, fault_shr_faces,
                                 nbf_per_face);

      // Parity Phase 5: register velocity (3-component L2 p=order) and
      // mpi_rank (L2 p=0).  Only when the volume mode is opted in.
      if (volume_pv_enabled)
      {
         pv_vel_fec = std::make_unique<L2_FECollection>(cfg.mesh.order, 3,
                                                        BasisType::GaussLobatto);
         pv_vel_fes = std::make_unique<ParFiniteElementSpace>(
                         &pmesh, pv_vel_fec.get(), 3, Ordering::byNODES);
         pv_vel_gf  = std::make_unique<ParGridFunction>(pv_vel_fes.get());
         *pv_vel_gf = 0.0;
         pv_out->RegisterDomainField("velocity", pv_vel_gf.get());

         pv_rank_fec = std::make_unique<L2_FECollection>(0, 3);
         pv_rank_fes = std::make_unique<ParFiniteElementSpace>(
                          &pmesh, pv_rank_fec.get());
         pv_rank_gf  = std::make_unique<ParGridFunction>(pv_rank_fes.get());
         *pv_rank_gf = static_cast<real_t>(rank);
         pv_out->RegisterDomainField("mpi_rank", pv_rank_gf.get());
      }

      // Parity Phase 6: SAFS fault static fields.
      if (num_fault_total > 0)
      {
         pv_lsw_mu_s     .SetSize(num_fault_total);
         pv_lsw_mu_d     .SetSize(num_fault_total);
         pv_lsw_d_c      .SetSize(num_fault_total);
         pv_nuc_amplitude.SetSize(num_fault_total);
         pv_nuc_radial   .SetSize(num_fault_total);
         pv_sig_n_init   .SetSize(num_fault_total);
         pv_tau1_init    .SetSize(num_fault_total);
         pv_tau2_init    .SetSize(num_fault_total);
         for (int i = 0; i < num_fault_total; ++i)
         {
            const DOFData &d = dof_data[i];
            // LSW slots only carry meaning for the LSW path; for RS the
            // values are zero per `copy_lsw_and_forced_rupture_fields`
            // -> `_RS` init.
            pv_lsw_mu_s(i) = d.lsw_mu_s;
            pv_lsw_mu_d(i) = d.lsw_mu_d;
            pv_lsw_d_c (i) = d.lsw_d_c;
            // Sum contributions from all nucleation kinds; exactly one
            // resolver returned non-empty Vectors above.
            real_t ad = 0.0, as = 0.0, rad = 0.0;
            if (nuc_params.amplitude_dip.Size()    == num_fault_total)
            { ad = nuc_params.amplitude_dip(i); }
            if (nuc_params.amplitude_strike.Size() == num_fault_total)
            { as = nuc_params.amplitude_strike(i); }
            if (nuc_params.radial.Size()           == num_fault_total)
            { rad = nuc_params.radial(i); }
            if (sq_nuc_params.amplitude_dip.Size()    == num_fault_total)
            { ad += sq_nuc_params.amplitude_dip(i); }
            if (sq_nuc_params.amplitude_strike.Size() == num_fault_total)
            { as += sq_nuc_params.amplitude_strike(i); }
            if (sq_nuc_params.radial.Size()           == num_fault_total)
            { rad += sq_nuc_params.radial(i); }
            if (ic_nuc_params.amplitude_dip.Size()    == num_fault_total)
            { ad += ic_nuc_params.amplitude_dip(i); }
            if (ic_nuc_params.amplitude_strike.Size() == num_fault_total)
            { as += ic_nuc_params.amplitude_strike(i); }
            if (ic_nuc_params.radial.Size()           == num_fault_total)
            { rad += ic_nuc_params.radial(i); }
            if (cc_nuc_params.amplitude_dip.Size()    == num_fault_total)
            { ad += cc_nuc_params.amplitude_dip(i); }
            if (cc_nuc_params.amplitude_strike.Size() == num_fault_total)
            { as += cc_nuc_params.amplitude_strike(i); }
            if (cc_nuc_params.radial.Size()           == num_fault_total)
            { rad += cc_nuc_params.radial(i); }
            pv_nuc_amplitude(i) = std::sqrt(ad*ad + as*as);
            pv_nuc_radial   (i) = rad;
            pv_sig_n_init   (i) = d.sigma_n_corr;
            pv_tau1_init    (i) = d.tau1_corr;
            pv_tau2_init    (i) = d.tau2_corr;
         }
      }
      // Parity Phase 6: publish all 8 SAFS static fields via the new
      // SetFaultParamsSpatial method on seas::ParaViewOutput (added in
      // io/paraview_output.hpp).  Each field lands in fault.vtkhdf
      // under its proper SAFS-semantic name (NOT TPV104's a/Dc/x2/x3).
      // SetFaultParamsSpatial is safe to call on a rank with zero
      // local fault DOFs.
      pv_out->SetFaultParamsSpatial({
         {"lsw_mu_s",          &pv_lsw_mu_s},
         {"lsw_mu_d",          &pv_lsw_mu_d},
         {"lsw_d_c",           &pv_lsw_d_c},
         {"nuc_amplitude",     &pv_nuc_amplitude},
         {"nuc_radial_factor", &pv_nuc_radial},
         {"sigma_n_init",      &pv_sig_n_init},
         {"tau1_init",         &pv_tau1_init},
         {"tau2_init",         &pv_tau2_init},
      });
   }

   // Parity Phase 4: secondary `pv_bulk_out` collection (6 stress
   // components → ParaView_bulk/stress.vtkhdf).
   if (bulk_pv_active)
   {
      const std::string bulk_dir = cfg.output.output_dir + "/ParaView_bulk";
      if (rank == 0)
      {
         std::filesystem::create_directories(bulk_dir);
      }
#ifdef MFEM_USE_MPI
      MPI_Barrier(comm);
#endif
      bool bulk_mode_enabled = bulk_pv_enabled;
      auto bulk_mode = ParseVolumeMode(cfg.output.paraview_bulk,
                                       bulk_mode_enabled);
      pv_bulk_out = std::make_unique<seas::ParaViewOutput<ParMesh>>(
                       bulk_dir, pmesh, cfg.mesh.order,
                       "stress", bulk_mode);
      pv_bulk_out->fixed_dt = cfg.output.paraview_bulk_dt;

      pv_bulk_sigma_fec = std::make_unique<L2_FECollection>(
                            cfg.mesh.order, 3, BasisType::GaussLobatto);
      pv_bulk_sigma_fes = std::make_unique<ParFiniteElementSpace>(
                            &pmesh, pv_bulk_sigma_fec.get());
      pv_bulk_sxx_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      pv_bulk_syy_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      pv_bulk_szz_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      pv_bulk_sxy_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      pv_bulk_syz_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      pv_bulk_sxz_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      *pv_bulk_sxx_gf = 0.0; *pv_bulk_syy_gf = 0.0; *pv_bulk_szz_gf = 0.0;
      *pv_bulk_sxy_gf = 0.0; *pv_bulk_syz_gf = 0.0; *pv_bulk_sxz_gf = 0.0;
      pv_bulk_out->RegisterDomainField("sigma_xx", pv_bulk_sxx_gf.get());
      pv_bulk_out->RegisterDomainField("sigma_yy", pv_bulk_syy_gf.get());
      pv_bulk_out->RegisterDomainField("sigma_zz", pv_bulk_szz_gf.get());
      pv_bulk_out->RegisterDomainField("sigma_xy", pv_bulk_sxy_gf.get());
      pv_bulk_out->RegisterDomainField("sigma_yz", pv_bulk_syz_gf.get());
      pv_bulk_out->RegisterDomainField("sigma_xz", pv_bulk_sxz_gf.get());

#ifdef MFEM_USE_HDF5
      if (cfg.output.paraview_bulk_zfp_tol > 0.0
          && bulk_mode == ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5)
      {
         pv_bulk_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            cfg.output.paraview_bulk_zfp_tol);
      }
      else if (cfg.output.paraview_bulk_deflate_level >= 0
               && bulk_mode == ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5)
      {
         pv_bulk_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::Deflate,
            static_cast<double>(cfg.output.paraview_bulk_deflate_level));
      }
#endif
   }

   // Parity Phase 3: rank-0 banner.
   if (rank == 0)
   {
      if (any_pv_requested)
      {
         std::cout << "ParaView output: ON (prefix=" << cfg.output.output_dir
                   << ")\n"
                   << "  Volume:   " << cfg.output.paraview_volume
                   << " (dt=" << cfg.output.paraview_volume_dt << " s)\n"
                   << "  Bulk:     " << cfg.output.paraview_bulk
                   << " (dt=" << cfg.output.paraview_bulk_dt   << " s)\n"
                   << "  Fault:    " << cfg.output.paraview_fault
                   << " (dt=" << cfg.output.paraview_fault_dt  << " s)\n";
         if (cfg.output.paraview_volume_zfp_tol > 0.0)
         {
            std::cout << "  Volume ZFP tol: "
                      << cfg.output.paraview_volume_zfp_tol << "\n";
         }
         if (cfg.output.paraview_bulk_zfp_tol > 0.0)
         {
            std::cout << "  Bulk   ZFP tol: "
                      << cfg.output.paraview_bulk_zfp_tol   << "\n";
         }
         if (cfg.output.paraview_fault_zfp_tol > 0.0)
         {
            std::cout << "  Fault  ZFP tol: "
                      << cfg.output.paraview_fault_zfp_tol  << "\n";
         }
         if (cfg.output.paraview_coseismic_dt > 0.0
             || cfg.output.paraview_nucleation_dt > 0.0
             || cfg.output.paraview_interseismic_dt > 0.0)
         {
            std::cout << "  Regime-adaptive cadence active.\n";
         }
         if (cfg.output.paraview_enabled
             && !primary_pv_active && !bulk_pv_active)
         {
            std::cout << "  WARNING: paraview_enabled = true but every "
                      << "per-collection mode is 'off' — no ParaView files "
                      << "will be written.\n";
         }
      }
      else
      {
         std::cout << "ParaView output: OFF\n";
      }
   }

   // ParaView GFs are allocated for the time loop only; sized after
   // pv_out is up.  No-op when num_fault_total = 0 on this rank.
   Vector pv_local_slip, pv_local_slip_rate, pv_local_traction;
   Vector pv_local_state, pv_local_normal_stress;
   if (pv_out)
   {
      pv_local_slip.SetSize(2 * num_fault_total);
      pv_local_slip_rate.SetSize(2 * num_fault_total);
      pv_local_traction.SetSize(2 * num_fault_total);
      pv_local_state.SetSize(num_fault_total);
      pv_local_normal_stress.SetSize(num_fault_total);
   }

   // -----------------------------------------------------------------
   // 18. Restart (R-105 corrected schema: per-rank Q + dof_data).
   // -----------------------------------------------------------------
   const int ndof_total = wave.GetScalarNDof();
   Vector Q(NUM_STATE * ndof_total);
   Q = 0.0;

   real_t t      = cfg.time.t_initial;
   int    step0  = 0;
   real_t dt_now = dt;
   if (!restart_prefix.empty())
   {
      std::string driver_tag;
      const bool ok = ReadTpv104Checkpoint(
         restart_prefix, t, dt_now, step0, Q, NUM_STATE * ndof_total,
         dof_data, rank, nprocs
#ifdef MFEM_USE_MPI
         , comm
#endif
         , &driver_tag);
      if (!ok)
      {
         if (rank == 0)
         {
            std::cerr << "ERROR: restart prefix '" << restart_prefix
                      << "' has no per-rank file for rank " << rank
                      << ".\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 5;
      }
      if (!driver_tag.empty() && driver_tag != "spatial_dyn")
      {
         if (rank == 0)
         {
            std::cerr << "ERROR: refusing restart: checkpoint driver_tag '"
                      << driver_tag << "' != 'spatial_dyn'.\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 5;
      }
      if (driver_tag.empty() && rank == 0)
      {
         std::cout << "[restart] WARNING: restoring pre-tag checkpoint "
                   << "into spatial_dyn — proceed only if you know the "
                   << "checkpoint is compatible.\n";
      }
      // SetTime so the LSW_ForcedRupture guard sees a fresh t.
      wave.SetTime(t);
   }

   // -----------------------------------------------------------------
   // 19. Sub-step iterator — friction-law-dependent.  LSW uses the
   //     TPV205 closed-form integrator; rate-state uses TPV102 (aging
   //     law) or TPV104 (slip-law-strong-rate-weakening) iterators.
   //     Each iterator's callback-aware `AdvanceWithSubStepStates`
   //     drives the per-sub-step nucleation hook below.
   //
   // REVIEW R-005: `[numerics] fault_iterator` selects the sub-step
   // quadrature:
   //   "one-shot" — O = 1 (native TPV205 default; byte-parity).
   //   "substep"  — O = ader_order (per-sub-step ADER quadrature).
   // -----------------------------------------------------------------
   const bool substep_quadrature =
      (cfg.numerics.fault_iterator == "substep");
   const int substep_O = substep_quadrature
      ? std::max(1, cfg.numerics.ader_order)
      : 1;
   const std::vector<real_t> substep_deltaT(substep_O,
                                            dt / static_cast<real_t>(substep_O));
   const std::vector<real_t> substep_weights(substep_O,
                                             1.0 / static_cast<real_t>(substep_O));
   if (rank == 0)
   {
      std::cout << "[time] fault iterator quadrature: O = " << substep_O
                << " (" << cfg.numerics.fault_iterator << ")\n";
   }

   std::unique_ptr<Tpv205SubStepIterator>  iter_lsw;
   std::unique_ptr<AgingLawPsi>            rs_aging_state;
   std::unique_ptr<Tpv102SubStepIterator>  iter_rs_aging;
   std::unique_ptr<SlipLawSRWPsi>          rs_srw_state;
   std::unique_ptr<Tpv104SubStepIterator>  iter_rs_srw;

   if (is_lsw)
   {
      iter_lsw = std::make_unique<Tpv205SubStepIterator>(fault_flux);
      iter_lsw->SetSubSteps(substep_deltaT, substep_weights);
   }
   else if (!rs_use_srw)
   {
      // AgingLawPsi takes (b, V0, f0).  These are GLOBAL scalars (the
      // resolver's per-DOF arrays are read elsewhere); seed from the
      // block defaults so the global b/V0/f0 used inside the state
      // analytic update match the per-DOF resolver output.  For TPV102
      // these defaults already match the spec.
      rs_aging_state = std::make_unique<AgingLawPsi>(
         cfg.rate_state->b_default,
         cfg.rate_state->V_0_default,
         cfg.rate_state->f_0_default);
      iter_rs_aging = std::make_unique<Tpv102SubStepIterator>(
         fault_flux, *rs_aging_state);
      iter_rs_aging->SetSubSteps(substep_deltaT, substep_weights);
   }
   else
   {
      // SlipLawSRWPsi takes (a, b, V0, f0, muW, V_w_default).  The
      // `a` and `V_w_default` here are FALLBACKS — the per-QP friction
      // pipeline reads `d.a` from DOFData and `V_w[i]` from the
      // `Vw_per_dof` side-channel; these scalars are consumed only by
      // the base-class virtuals (which the production-mode guard
      // disables — see SlipLawSRWPsi::SetProductionMode).
      rs_srw_state = std::make_unique<SlipLawSRWPsi>(
         cfg.rate_state->a_default,
         cfg.rate_state->b_default,
         cfg.rate_state->V_0_default,
         cfg.rate_state->f_0_default,
         cfg.rate_state->f_w_default,
         cfg.rate_state->V_w_default);
      iter_rs_srw = std::make_unique<Tpv104SubStepIterator>(
         fault_flux, *rs_srw_state);
      iter_rs_srw->SetSubSteps(substep_deltaT, substep_weights);
   }

   // Per-sub-step nucleation accumulator hook.  Closes over the
   // four `*_nuc_params` resolver outputs + cfg; fires BEFORE each
   // sub-step's friction solve.
   auto nuc_cb = [&nuc_params, &sq_nuc_params, &cc_nuc_params, &dof_data, &cfg]
                 (real_t t_sub_end, real_t dt_sub)
   {
      if (!cfg.nucleation.enabled) { return; }
      switch (cfg.nucleation.kind)
      {
      case spatial::NucleationKind::GradualOverstress:
         spatial::ApplyGradualOverstressIncrement(
            dof_data, nuc_params,
            cfg.nucleation.gradual_overstress.T_nuc_s,
            t_sub_end, dt_sub);
         break;
      case spatial::NucleationKind::SquareOverstress:
         spatial::ApplySquareOverstressIncrement(
            dof_data, sq_nuc_params,
            cfg.nucleation.square_overstress.T_nuc_s,
            t_sub_end, dt_sub);
         break;
      case spatial::NucleationKind::GradualOverstressCompactCircular:
         spatial::ApplyGradualOverstressCompactCircularIncrement(
            dof_data, cc_nuc_params,
            cfg.nucleation.gradual_overstress_compact_circular.T_nuc_s,
            t_sub_end, dt_sub);
         break;
      case spatial::NucleationKind::InstantaneousOverstressCircular:
         // Applied ONCE at init; per-sub-step hook is a no-op.
         break;
      }
   };

   // ParaView snapshot writer (Parity Phases 1-6).  Updates the 5 BP5
   // fault projection GFs + writes the primary / bulk collections.
   // Splits writes into "fault collection wants" (primary pv_out) and
   // "bulk wants" (pv_bulk_out) so each can fire on its own cadence.
   auto paraview_write = [&](int step_num, real_t time, real_t V_max)
   {
      const bool fault_wants = pv_out
                              && pv_out->PeekShouldWrite(step_num, time, V_max);
      const bool bulk_wants  = pv_bulk_out
                              && pv_bulk_out->PeekShouldWrite(step_num, time, V_max);
      if (!fault_wants && !bulk_wants) { return; }

      // Parity Phase 4: publish 6 stress components to pv_bulk_out.
      if (bulk_wants)
      {
         const real_t* Q_data = Q.GetData();
         std::memcpy(pv_bulk_sxx_gf->GetData(), Q_data + SXX * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_syy_gf->GetData(), Q_data + SYY * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_szz_gf->GetData(), Q_data + SZZ * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_sxy_gf->GetData(), Q_data + SXY * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_syz_gf->GetData(), Q_data + SYZ * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_sxz_gf->GetData(), Q_data + SXZ * ndof_total,
                     ndof_total * sizeof(real_t));
         pv_bulk_out->ForceSave(step_num, time);
         pv_bulk_out->CommitSchedule(time, V_max);
      }

      if (!fault_wants) { return; }

      for (int i = 0; i < num_fault_total; ++i)
      {
         const DOFData &d = dof_data[i];
         pv_local_slip(2 * i + 0)      = d.slip1;
         pv_local_slip(2 * i + 1)      = d.slip2;
         pv_local_slip_rate(2 * i + 0) = d.V1;
         pv_local_slip_rate(2 * i + 1) = d.V2;
         pv_local_traction(2 * i + 0)  = d.tau1_corr;
         pv_local_traction(2 * i + 1)  = d.tau2_corr;
         if (is_lsw)
         {
            const real_t delta_norm = std::sqrt(d.slip1 * d.slip1
                                                + d.slip2 * d.slip2);
            pv_local_state(i) = mfem::seas::LSWFrictionCoefficient_TPV205(
                                   delta_norm,
                                   d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c);
         }
         else
         {
            // Rate-state: publish the state variable ψ (the spatial
            // driver does not maintain a per-DOF Dc in DOFData for the
            // LSW path, so the LSW friction coefficient is undefined —
            // ψ is the natural diagnostic).
            pv_local_state(i) = d.psi;
         }
         (void)time;
         pv_local_normal_stress(i)     = d.sigma_n_corr;
      }

      // Parity Phase 5: refresh the volume velocity field directly from
      // Q (byNODES so VX/VY/VZ are contiguous 3*ndof block).
      if (pv_out->GetVolumeSaveEnabled() && pv_vel_gf)
      {
         std::memcpy(pv_vel_gf->GetData(),
                     Q.GetData() + VX * ndof_total,
                     3 * ndof_total * sizeof(real_t));
      }

      if (pv_out->GetVolumeSaveEnabled())
      {
         pv_out->UpdateFaultFieldsBP5(pv_local_slip, pv_local_slip_rate,
                                      pv_local_traction, pv_local_state,
                                      pv_local_normal_stress);
         pv_out->ForceSave(step_num, time);
         pv_out->CommitSchedule(time, V_max);
      }
      else
      {
         pv_out->CommitSchedule(time, V_max);
      }
   };

   if (restart_prefix.empty())
   {
      paraview_write(0, cfg.time.t_initial, 0.0);
   }

   // -----------------------------------------------------------------
   // 19b. SCEC TPV205 station-trace writer (REVIEW R-009 parity with
   //      drivers/tpv205_driver.cpp:1851-1875).  Active only when the
   //      TOML's [problem].tag == "tpv205" so non-TPV205 SAFS / TPV31
   //      runs see no behaviour change.  The 16 station files match
   //      the SCEC TPV5 benchmark-trace filename convention so
   //      downstream comparison scripts can consume the spatial
   //      driver's output identically to the native driver's.
   // -----------------------------------------------------------------
   TPV205StationWriter tpv205_station_writer;
   const bool tpv205_stations_active = (cfg.problem.tag == "tpv205");
   if (tpv205_stations_active)
   {
      const std::vector<TPV205Station> stations = DefaultStations_TPV205();
#ifdef MFEM_USE_MPI
      tpv205_station_writer.Open(cfg.output.output_dir,
                                 /*prefix=*/"tpv205",
                                 stations, fault_coords,
                                 num_fault_local, comm);
#else
      tpv205_station_writer.Open(cfg.output.output_dir,
                                 /*prefix=*/"tpv205",
                                 stations, fault_coords, num_fault_local);
#endif
      if (restart_prefix.empty())
      {
         tpv205_station_writer.WriteStep(cfg.time.t_initial, dof_data);
      }
      if (rank == 0)
      {
         std::cout << "[stations] TPV205 station writer active ("
                   << stations.size() << " stations, prefix=tpv205_)\n";
      }
   }

   // REVIEW R-007: TPV102 SCEC station traces — active only when
   // `[problem].tag == "tpv102"`.  Mirrors the TPV205 wiring above.
   TPV102StationWriter tpv102_station_writer;
   const bool tpv102_stations_active = (cfg.problem.tag == "tpv102");
   if (tpv102_stations_active)
   {
      const std::vector<TPV102Station> stations = DefaultStations();
#ifdef MFEM_USE_MPI
      tpv102_station_writer.Open(cfg.output.output_dir,
                                 /*prefix=*/"tpv102",
                                 stations, fault_coords,
                                 num_fault_local, comm);
#else
      tpv102_station_writer.Open(cfg.output.output_dir,
                                 /*prefix=*/"tpv102",
                                 stations, fault_coords, num_fault_local);
#endif
      if (restart_prefix.empty())
      {
         tpv102_station_writer.WriteStep(cfg.time.t_initial, dof_data);
      }
      if (rank == 0)
      {
         std::cout << "[stations] TPV102 station writer active ("
                   << stations.size() << " stations, prefix=tpv102_)\n";
      }
   }

   // REVIEW R-007: TPV104 SCEC station traces — active only when
   // `[problem].tag == "tpv104"`.
   TPV104StationWriter tpv104_station_writer;
   const bool tpv104_stations_active = (cfg.problem.tag == "tpv104");
   if (tpv104_stations_active)
   {
      const std::vector<TPV104Station> stations = DefaultStations_TPV104();
#ifdef MFEM_USE_MPI
      tpv104_station_writer.Open(cfg.output.output_dir,
                                 /*prefix=*/"tpv104",
                                 stations, fault_coords,
                                 num_fault_local, comm);
#else
      tpv104_station_writer.Open(cfg.output.output_dir,
                                 /*prefix=*/"tpv104",
                                 stations, fault_coords, num_fault_local);
#endif
      if (restart_prefix.empty())
      {
         tpv104_station_writer.WriteStep(cfg.time.t_initial, dof_data);
      }
      if (rank == 0)
      {
         std::cout << "[stations] TPV104 station writer active ("
                   << stations.size() << " stations, prefix=tpv104_)\n";
      }
   }

   // -----------------------------------------------------------------
   // 20. Time loop.
   // -----------------------------------------------------------------
   Vector Q_new(Q.Size());
   real_t V_max_global = 0.0;
   // R-603 round-6: track the last completed step so the final
   // checkpoint at L1262 records the actual step the loop reached,
   // not `nsteps` unconditionally.
   int last_completed_step = step0;

   // Iterator-specific lambdas (one set captures whichever iterator
   // pointer the cfg.law branch above instantiated).  `set_substeps`
   // forwards to the iterator's SetSubSteps; `do_iterate` forwards to
   // its callback-aware `AdvanceWithSubStepStates`.
   auto set_substeps =
      [&iter_lsw, &iter_rs_aging, &iter_rs_srw]
      (const std::vector<real_t> &deltaT,
       const std::vector<real_t> &weights)
   {
      if (iter_lsw)      { iter_lsw->SetSubSteps(deltaT, weights); }
      else if (iter_rs_aging) { iter_rs_aging->SetSubSteps(deltaT, weights); }
      else if (iter_rs_srw)   { iter_rs_srw->SetSubSteps(deltaT, weights); }
   };
   auto do_iterate =
      [&iter_lsw, &iter_rs_aging, &iter_rs_srw,
       &dof_data, &fault_coords, &Vw_per_dof]
      (const std::vector<std::vector<real_t>> &Q_pointwise_plus,
       const std::vector<std::vector<real_t>> &Q_pointwise_minus,
       real_t dt_step, real_t t_step_start,
       real_t *I_imp_plus, real_t *I_imp_minus,
       const std::function<void(real_t, real_t)> &nuc_callback_inner)
   {
      if (iter_lsw)
      {
         iter_lsw->AdvanceWithSubStepStates(
            dof_data, fault_coords,
            Q_pointwise_plus, Q_pointwise_minus,
            dt_step, t_step_start,
            I_imp_plus, I_imp_minus,
            nuc_callback_inner);
      }
      else if (iter_rs_aging)
      {
         iter_rs_aging->AdvanceWithSubStepStates(
            dof_data, fault_coords,
            Q_pointwise_plus, Q_pointwise_minus,
            dt_step, t_step_start,
            I_imp_plus, I_imp_minus,
            nuc_callback_inner);
      }
      else if (iter_rs_srw)
      {
         iter_rs_srw->AdvanceWithSubStepStates(
            dof_data, fault_coords, Vw_per_dof,
            Q_pointwise_plus, Q_pointwise_minus,
            dt_step, t_step_start,
            I_imp_plus, I_imp_minus,
            nuc_callback_inner);
      }
   };

   for (int step = step0; step < nsteps; ++step)
   {
      const real_t dt_step = std::min(dt_now, cfg.time.tfinal - t);
      if (dt_step <= 0.0) { break; }
      wave.SetTime(t);

      AdvanceADERWithSubStep_Spatial(wave, set_substeps, do_iterate,
                                     substep_deltaT, substep_weights,
                                     Q, dt_step,
                                     cfg.numerics.ader_order, t, Q_new,
                                     nuc_cb);
      Q.Swap(Q_new);
      t += dt_step;
      last_completed_step = step + 1;

      // REVIEW R-010: NaN tripwire — match the native TPV205 driver's
      // per-step check at `drivers/tpv205_driver.cpp:2447-2464`.  Per
      // `CLAUDE.md` "What Constitutes a Regression" #3, dt going to
      // zero or NaN is a hard regression signature; without this
      // tripwire the simulation would run to completion writing
      // garbage.
      {
         real_t local_nan = std::isnan(Q.Norml2()) ? 1.0 : 0.0;
         real_t global_nan = local_nan;
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&local_nan, &global_nan, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
         if (global_nan > 0.0)
         {
            if (rank == 0)
            {
               std::cerr << "ERROR: NaN detected at step " << step
                         << ", t = " << t << " s (REVIEW R-010).\n";
            }
            if (tpv205_stations_active)
            {
               tpv205_station_writer.Flush();
               tpv205_station_writer.Close();
            }
            if (tpv102_stations_active)
            {
               tpv102_station_writer.Flush();
               tpv102_station_writer.Close();
            }
            if (tpv104_stations_active)
            {
               tpv104_station_writer.Flush();
               tpv104_station_writer.Close();
            }
#ifdef MFEM_USE_MPI
            MPI_Finalize();
#endif
            return 1;
         }
      }

      real_t V_max_local = 0.0;
      for (int i = 0; i < num_fault_total; ++i)
      {
         V_max_local = std::max(V_max_local, dof_data[i].slip_rate);
      }
      real_t V_max_step = V_max_local;
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&V_max_local, &V_max_step, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
      V_max_global = std::max(V_max_global, V_max_step);

      paraview_write(step + 1, t, V_max_step);

      // REVIEW R-009: per-step SCEC trace write (active only when
      // problem.tag == "tpv205").  Cadence is per-step here (matches
      // the spatial driver's paraview_fault_dt = 0.001s default in
      // TPV205 configs); the native TPV205 driver downsamples via
      // `output_dt = 0.01s` but per-step is a strict superset.
      if (tpv205_stations_active)
      {
         tpv205_station_writer.WriteStep(t, dof_data);
      }
      // REVIEW R-007: per-step SCEC trace write for TPV102 / TPV104.
      if (tpv102_stations_active)
      {
         tpv102_station_writer.WriteStep(t, dof_data);
      }
      if (tpv104_stations_active)
      {
         tpv104_station_writer.WriteStep(t, dof_data);
      }

      if (cfg.output.checkpoint_every_steps > 0
          && (step + 1) % cfg.output.checkpoint_every_steps == 0)
      {
         const std::string prefix = cfg.output.output_dir + "/"
                                    + cfg.output.restart_prefix;
         WriteTpv104Checkpoint(prefix, t, dt_now, step + 1, Q, dof_data,
                               rank, nprocs
#ifdef MFEM_USE_MPI
                               , comm
#endif
                               , "spatial_dyn");
      }

      if (rank == 0 && (step % 100 == 0 || step == nsteps - 1))
      {
         std::cout << "step " << step << "/" << nsteps
                   << "  t = " << t << " s"
                   << "  V_max = " << V_max_step << " m/s\n";
      }
   }

   // -----------------------------------------------------------------
   // 21. Final checkpoint.
   // -----------------------------------------------------------------
   if (cfg.output.checkpoint_every_steps > 0)
   {
      const std::string prefix = cfg.output.output_dir + "/"
                                 + cfg.output.restart_prefix;
      // R-603 round-6: use last_completed_step instead of nsteps so the
      // checkpoint reports the actual step the loop reached.
      WriteTpv104Checkpoint(prefix, t, dt_now, last_completed_step,
                            Q, dof_data,
                            rank, nprocs
#ifdef MFEM_USE_MPI
                            , comm
#endif
                            , "spatial_dyn");
   }

   // REVIEW R-007 / R-009: flush + close SCEC trace files for whichever
   // benchmark is active.  No-op for the inactive writers.
   if (tpv205_stations_active)
   {
      tpv205_station_writer.Flush();
      tpv205_station_writer.Close();
   }
   if (tpv102_stations_active)
   {
      tpv102_station_writer.Flush();
      tpv102_station_writer.Close();
   }
   if (tpv104_stations_active)
   {
      tpv104_station_writer.Flush();
      tpv104_station_writer.Close();
   }

   if (rank == 0)
   {
      std::cout << "[spatial_dyn] done.  Final t = " << t
                << " s, V_max_global = " << V_max_global << " m/s\n";
   }

#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return 0;
}
