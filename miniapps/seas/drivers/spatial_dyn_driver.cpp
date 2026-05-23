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

   // R-005: per-QP basis for curvilinear-face fidelity.  When the caller
   // populated `bdata.qp_data` (via `FaultBasis::ComputeQPBasis` /
   // `ComputeQPBasisShared`) the per-QP frame at QP `q` is used; otherwise
   // we fall back to the centroid frame (planar-face accuracy, byte-
   // compatible with the pre-R-005 inline walk).
   auto qp_normal = [](const FaultBasisData &bdata, int q, int d) -> real_t
   {
      return (static_cast<int>(bdata.qp_data.size()) > q)
                ? bdata.qp_data[q].normal[d]   : bdata.normal[d];
   };
   auto qp_t1 = [](const FaultBasisData &bdata, int q, int d) -> real_t
   {
      return (static_cast<int>(bdata.qp_data.size()) > q)
                ? bdata.qp_data[q].tangent1[d] : bdata.tangent1[d];
   };
   auto qp_t2 = [](const FaultBasisData &bdata, int q, int d) -> real_t
   {
      return (static_cast<int>(bdata.qp_data.size()) > q)
                ? bdata.qp_data[q].tangent2[d] : bdata.tangent2[d];
   };

   auto write_dof = [&](int dof_idx, FaceElementTransformations *ftr,
                        const IntegrationPoint &ip,
                        const FaultBasisData &bdata,
                        int q,
                        int attr)
   {
      ftr->SetAllIntPoints(&ip);
      Vector phys(3);
      ftr->Face->Transform(ip, phys);
      phys_coords.push_back(phys);

      for (int d = 0; d < 3; ++d)
      {
         dof_coords_3d(3 * dof_idx + d) = phys(d);
         dof_basis(0 + d, dof_idx) = qp_normal(bdata, q, d);
         dof_basis(3 + d, dof_idx) = qp_t1    (bdata, q, d);
         dof_basis(6 + d, dof_idx) = qp_t2    (bdata, q, d);
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
         write_dof(dof_idx++, ftr, ir.IntPoint(q), bdata, q, fault_attr);
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
         write_dof(dof_idx++, ftr, ir.IntPoint(q), bdata, q, fault_attr);
      }
   }
#else
   (void)fault_shr_faces;
#endif
   MFEM_VERIFY(dof_idx == N,
               "spatial_dyn: per-DOF table walk wrote " << dof_idx
               << " entries, expected " << N);
}

// TPV205-style ADER macro-step driver: predictor in the bulk, LSW
// closed-form per sub-step at fault QPs, corrector via wave.AdvanceADER
// with the side-channel I_imp.  Mirrors AdvanceADERWithSubStep in
// drivers/tpv205_driver.cpp (1:1 except no SEAS_DIAG hooks).
//
// Phase N: the trailing `nuc_callback` arg is forwarded to the
// per-sub-step callback overload of
// `Tpv205SubStepIterator::AdvanceWithSubStepStates`; the callback fires
// ONCE per ADER sub-step BEFORE the per-QP friction pipeline.  Pass
// `[](real_t, real_t){}` to opt out (no nucleation perturbation).
void AdvanceADERWithSubStep_Spatial(
   WaveOperator<ParMesh> &wave,
   Tpv205SubStepIterator &iterator,
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
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

   const std::vector<real_t> configured_deltaT  = iterator.GetDeltaT();
   const std::vector<real_t> configured_weights = iterator.GetTimeWeights();
   const int O = static_cast<int>(configured_deltaT.size());
   MFEM_VERIFY(O >= 1,
               "AdvanceADERWithSubStep_Spatial: iterator deltaT empty; "
               "SetSubSteps must be called first.");
   const real_t configured_sum =
      std::accumulate(configured_deltaT.begin(), configured_deltaT.end(),
                      static_cast<real_t>(0));
   MFEM_VERIFY(configured_sum > 0.0,
               "AdvanceADERWithSubStep_Spatial: Σ deltaT = "
               << configured_sum << " ≤ 0");

   const real_t dt_scale = dt_step / configured_sum;
   std::vector<real_t> deltaT_scaled(O);
   for (int o = 0; o < O; ++o)
   {
      deltaT_scaled[o] = configured_deltaT[o] * dt_scale;
   }
   iterator.SetSubSteps(deltaT_scaled, configured_weights);
   const std::vector<real_t> &deltaT = iterator.GetDeltaT();

   // Sub-step MIDPOINT nodes on [0, dt_step].
   std::vector<real_t> tau_nodes(O);
   real_t acc = 0.0;
   for (int o = 0; o < O; ++o)
   {
      tau_nodes[o] = acc + 0.5 * deltaT[o];
      acc += deltaT[o];
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
      iterator.AdvanceWithSubStepStates(dof_data, fault_coords,
                                        Q_pointwise_plus,
                                        Q_pointwise_minus,
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
   MFEM_VERIFY(is_lsw,
               "spatial_dyn_driver: only [meta].law = \"slip_weakening\" "
               "is supported in this commit; rate_state path is a "
               "deferred follow-up (plan §Phase 4 Edge Cases).");

   if (rank == 0)
   {
      std::cout << "================================================\n"
                << "seas_spatial_dyn_driver — Phase 4 (rev-3)\n"
                << "================================================\n"
                << "config:           " << config_path << "\n"
                << "mesh:             " << cfg.mesh.path << "\n"
                << "fe order:         " << cfg.mesh.order << "\n"
                << "law:              "
                << (is_lsw ? "slip_weakening" : "rate_state") << "\n"
                << "stress kind:      "
                << (cfg.stress.kind == spatial::StressSourceKind::ConstantTensor
                    ? "constant_tensor" : "sidecar_hdf5") << "\n"
                << "tfinal:           " << cfg.time.tfinal << " s\n"
                << "cfl:              " << cfg.numerics.cfl << "\n"
                << "ader order:       " << cfg.numerics.ader_order << "\n"
                << "mixed flux:       " << cfg.numerics.mixed_flux << "\n"
                << "use pml:          " << (cfg.numerics.use_pml ? "yes" : "no")
                << "\n"
                << "nucleation:       "
                << (cfg.nucleation.enabled
                    ? "gradual_overstress (enabled)"
                    : "DISABLED")
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
   // -----------------------------------------------------------------
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
   // 5.  BoundaryConfig (SAFS .geo: fault=101, top free=102, bottom +
   //     sides absorbing=103+104).
   // -----------------------------------------------------------------
   BoundaryConfig bc;
   bc.fault_attr = 101;
   bc.natural_attrs   = {102};
   bc.absorbing_attrs = {103, 104};

   // -----------------------------------------------------------------
   // 6.  Material — pick the right WaveOperator ctor (deviation D-1).
   // -----------------------------------------------------------------
   real_t mat_lambda = cfg.material_fallback.lambda;
   real_t mat_mu     = cfg.material_fallback.mu;
   real_t mat_rho    = cfg.material_fallback.rho;

   std::unique_ptr<spatial::SpatialVelocityBundle> vel_bundle;
   MaterialField material = MaterialField::MakeConstant(mat_lambda,
                                                        mat_mu, mat_rho);

   // Material-sidecar gating: TOML [velocity].use_sidecar AND no CLI
   // override.  The CLI flag always wins (it can force the constant
   // path even when the TOML asks for a sidecar load).
   const bool sidecar_requested =
      cfg.velocity.use_sidecar && !no_sidecar_material;

   // R-010 / D-1: heterogeneous WaveOperator(MaterialField) ctor +
   // per-element flux dispatch is not yet wired (see
   // dynamic/wave_operator.hpp:215).  A successful sidecar load would
   // produce a non-Constant MaterialField that the scalar-material
   // ctor below cannot consume.  Abort BEFORE the multi-minute HDF5
   // read so the user does not burn a Frontera dev-queue allocation
   // discovering this after the fact.
   if (sidecar_requested)
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: spatial_dyn_driver Phase H gap (D-1): the "
                   << "heterogeneous WaveOperator(MaterialField) ctor + "
                   << "per-element flux dispatch is NOT yet wired (see "
                   << "dynamic/wave_operator.hpp:215).  Re-run with "
                   << "--no-sidecar-material (or set "
                   << "[velocity].use_sidecar = false in the TOML) to "
                   << "use the [material_constant_fallback] block.\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 4;
   }

   if (sidecar_requested)
   {
      try
      {
         vel_bundle = std::make_unique<spatial::SpatialVelocityBundle>(
            spatial::LoadSpatialVelocityBundle(cfg.velocity, pmesh));
         material = vel_bundle->MakeMaterialField();
         if (rank == 0)
         {
            std::cout << "[material] loaded sidecar bundle: "
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
                      << "  Re-run with --no-sidecar-material (or set "
                      << "[velocity].use_sidecar = false in the TOML) "
                      << "to use the [material_constant_fallback] block.\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 4;
      }
   }
   else
   {
      if (rank == 0)
      {
         const char *gate_source =
            (!cfg.velocity.use_sidecar && no_sidecar_material) ? "TOML+CLI"
            : (!cfg.velocity.use_sidecar)                      ? "TOML"
            :                                                    "CLI";
         std::cout << "[material] using [material_constant_fallback] "
                   << "(gated by " << gate_source << "): "
                   << "lambda=" << mat_lambda
                   << " mu=" << mat_mu
                   << " rho=" << mat_rho << "\n";
      }
   }

   // Deviation D-1: WaveOperator(MaterialField) ctor is not yet wired.
   // Refuse to silently downgrade Mode::Coefficient input to scalar.
   MFEM_VERIFY(material.mode == MaterialField::Mode::Constant,
               "spatial_dyn_driver Phase H gap: the heterogeneous "
               "WaveOperator(MaterialField) ctor + per-element flux "
               "dispatch (plan §Phase H.1/H.2) is NOT yet wired (see "
               "dynamic/wave_operator.hpp:215 — SetGodunovFluxPool still "
               "aborts).  This driver supports only "
               "MaterialField::Mode::Constant via the existing scalar-"
               "material ctor.  Re-run with --no-sidecar-material to "
               "force the [material_constant_fallback] path.");

   // -----------------------------------------------------------------
   // 7.  Construct WaveOperator (scalar-material; deviation D-1).
   // -----------------------------------------------------------------
   WaveOperator<ParMesh> wave(pmesh, cfg.mesh.order,
                              material.lambda_const,
                              material.mu_const,
                              material.rho_const,
                              bc);

   // R-107 reflection-time warning: compute min_box_dim / cp_max from
   // mesh bounding box + scalar material.
   {
      const real_t cp = std::sqrt((material.lambda_const
                                   + 2.0 * material.mu_const)
                                   / material.rho_const);
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

   // Phase N: the spatial driver supports exactly one nucleation kind
   // (`gradual_overstress`) — the friction law is always plain LSW.
   // The gradual_overstress accumulator writes time-domain perturbations
   // into DOFData::tau{1,2}_nuc; the LSW solver consumes them via
   // s.tau{1,2}_total = tau{1,2}_0 + tau{1,2}_nuc + trial.  The obsolete
   // LSW_ForcedRupture dispatch arm + f_2(t) per-DOF friction reduction
   // are NOT used.  Native TPV* drivers continue to set
   // FaultFrictionLaw::LSW_ForcedRupture verbatim.
   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
   wave.SetMixedFluxMode(ParseMixedFlux(cfg.numerics.mixed_flux));

   // Prove the mixed-flux mode is NOT a silent no-op: report the GLOBAL
   // count of central-flux faces actually populated by
   // BuildCentralFluxFaceSet_.  The startup banner only echoes the config
   // STRING ("mixed flux: adjacent"); this prints what the operator really
   // built.  none -> 0 (upwind everywhere); adjacent / all_continuous -> > 0
   // iff fault-adjacent non-fault interior faces exist on this mesh.
   // (Rank-summed; a shared non-fault face is counted on each rank that
   // owns it, so seam faces are slightly over-counted — fine as an
   // "is it active" signal.  Compare none vs adjacent: 0 vs N.)
   {
      long long central_local =
         static_cast<long long>(wave.GetCentralFluxFaceSet().size());
      long long central_global = central_local;
#ifdef MFEM_USE_MPI
      MPI_Reduce(&central_local, &central_global, 1, MPI_LONG_LONG,
                 MPI_SUM, 0, comm);
#endif
      if (rank == 0)
      {
         std::cout << "[mixed-flux] mode = " << cfg.numerics.mixed_flux
                   << ", central-flux faces (rank-summed) = "
                   << central_global
                   << (central_global == 0
                       ? "  (upwind everywhere — mode is a NO-OP)"
                       : "  (central flux ACTIVE on these faces)")
                   << "\n";
      }
   }

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

   // Reconstruct FaultBasis on this rank (interior + shared appended).
   // IMPORTANT: this basis must match WaveOperator's runtime fault basis
   // exactly, otherwise tau_pre / sigma_n / nucleation are resolved in one
   // signed local frame and consumed by the friction solve in another.
   // WaveOperator uses the BP5 / Tandem convention ref_normal = -y and
   // up = +z (wave_operator.inl ctor path).  Keep the same convention here.
   Vector ref_normal(3);
   ref_normal(0) = 0.0;  ref_normal(1) = -1.0; ref_normal(2) = 0.0;
   Vector up_vec(3);
   up_vec(0)     = 0.0;  up_vec(1)     = 0.0;  up_vec(2)     = 1.0;

   FaultBasis fbasis;
   fbasis.Compute(pmesh, fault_int_faces, ref_normal, up_vec);
#ifdef MFEM_USE_MPI
   fbasis.AppendSharedFaces(pmesh, fault_shr_faces, ref_normal, up_vec);
#endif

   // R-005: populate per-QP basis vectors so curvilinear-face fidelity is
   // not lost when `BuildPerDOFFaultTables` stamps `dof_basis(*, q)`.
   // On a rank with no fault faces these calls are no-ops; on a rank with
   // only shared faces the IR geometry is probed from the first shared
   // face transformation.  `BuildPerDOFFaultTables` falls back to the
   // centroid frame if `qp_data` is empty (e.g., planar TPV-style mesh
   // where the probe IR doesn't match the per-face IR exactly).
   {
      FaceElementTransformations *ftr_probe = nullptr;
      if (fault_int_faces.Size() > 0)
      {
         ftr_probe = pmesh.GetInteriorFaceTransformations(fault_int_faces[0]);
      }
#ifdef MFEM_USE_MPI
      else if (fault_shr_faces.Size() > 0)
      {
         ftr_probe = pmesh.GetSharedFaceTransformations(fault_shr_faces[0]);
      }
#endif
      if (ftr_probe)
      {
         const IntegrationRule &qp_ir =
            IntRules.Get(ftr_probe->GetGeometryType(), 2 * cfg.mesh.order);
         fbasis.ComputeQPBasis(pmesh, fault_int_faces,
                               ref_normal, up_vec, qp_ir);
#ifdef MFEM_USE_MPI
         fbasis.ComputeQPBasisShared(pmesh, fault_shr_faces,
                                     ref_normal, up_vec, qp_ir,
                                     fault_int_faces.Size());
#endif
      }
   }

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
      geom.ComputeSAFSParams(src,
                             cfg.stress.pore_pressure.P_p_pa,
                             cfg.stress.pore_pressure.P_p_grad_pa_per_m,
                             cfg.stress.pore_pressure.min_sigma_n_pa);
   }
   else
   {
      spatial::ApplyCsmStressSidecar(cfg.stress, geom);
   }
   MFEM_VERIFY(geom.HasSAFSParams(),
               "spatial_dyn_driver: stress source projection failed");

   // -----------------------------------------------------------------
   // 11. Resolve per-DOF LSW parameters (Phase 1).
   // -----------------------------------------------------------------
   MFEM_VERIFY(cfg.slip_weakening.has_value(),
               "spatial_dyn_driver: [meta].law=slip_weakening but the "
               "[friction.slip_weakening] block is absent in TOML.");
   spatial::SpatialFrictionResolver resolver;
   const spatial::SlipWeakeningPerDOFParams lsw =
      resolver.ResolveSlipWeakening(*cfg.slip_weakening,
                                    dof_coords_3d, dof_to_attr);

   // -----------------------------------------------------------------
   // 12. Phase N: resolve the single nucleation kind
   //     (`gradual_overstress`).  Per-DOF amplitude_dip(i) /
   //     amplitude_strike(i) = F(r_i) · Δτ; per-DOF radial(i) = F(r_i).
   //     When `cfg.nucleation.enabled == false`, the resolver returns
   //     three zero-sized Vectors — the per-sub-step accumulator
   //     early-returns and the simulation runs with no nucleation.
   // -----------------------------------------------------------------
   const spatial::GradualOverstressPerDOFParams nuc_params =
      spatial::ResolveGradualOverstress(
         cfg.nucleation.gradual_overstress,
         cfg.nucleation.enabled,
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
   const real_t cp_seed = std::sqrt((material.lambda_const
                                     + 2.0 * material.mu_const)
                                     / material.rho_const);
   const real_t cs_seed = std::sqrt(material.mu_const / material.rho_const);
   FaultFaceFlux fault_flux(material.rho_const, cp_seed, cs_seed);
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
   if (num_fault_total > 0)
   {
      spatial::InitializeFaultDOFs_Spatial<ParMesh>(
         dof_data, num_fault_total, dof_to_elem, material, pmesh,
         lsw, geom.GetTauPre(), geom.sigma_n_per_dof(),
         dummy_T_forced, dummy_t0_decay,
         dof_ips);
   }
   wave.SetFaultDOFData(&dof_data, nbf_per_face);

   // Fluctuation-Q dispatch (matches TPV205): Q_bg = 0.
   {
      real_t Q_bg[NUM_STATE] = {0};
      wave.SetAbsorbingBackground(Q_bg);
   }

   // -----------------------------------------------------------------
   // 15. CFL / Δt and derived numbers.  ComputeMaxDt is the scalar-
   //     material implementation (deviation D-1).
   //
   //     The ADER-DG stable step carries a 1/(2N+1) spatial-order factor
   //     (N = cfg.mesh.order) plus a safety margin; ComputeMaxDt applies
   //     NEITHER (only its internal mixed-flux factor).  Mirror the
   //     verified tpv205_driver.cpp:1242 conversion so cfg.numerics.cfl
   //     is a TPV205-style safety knob, not the raw CFL number.  Passing
   //     0.5 raw stepped at 9x TPV205's dt -> above the N=1 stability
   //     boundary -> nucleation-end velocity blow-up (job 7743554).
   // -----------------------------------------------------------------
   const real_t dt_cfl = wave.ComputeMaxDt(
      cfg.numerics.cfl / (3.0 * (2.0 * cfg.mesh.order + 1.0)));
   real_t dt = (cfg.time.dt_initial > 0.0)
                ? cfg.time.dt_initial : dt_cfl;
   // R-006: warn if the user override exceeds the explicit CFL bound.
   // The ADER substep iterator is conditionally stable in dt; silently
   // crossing CFL produces nucleation-like blow-up with no error
   // pointing at the misconfiguration.
   if (cfg.time.dt_initial > 0.0 && cfg.time.dt_initial > dt_cfl)
   {
      if (rank == 0)
      {
         std::cerr << "[spatial_dyn] WARNING: [time].dt_initial = "
                   << cfg.time.dt_initial << " s exceeds CFL bound "
                   << dt_cfl << " s (cfl = " << cfg.numerics.cfl
                   << ").  ADER will be unconditionally unstable.  "
                   << "Tighten dt_initial or raise cfl.\n";
      }
   }
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

   if (print_derived)
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
         /*mu_bulk=*/material.mu_const,
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

      // R-003: opt back into legacy behaviour where the 8 L2-p0 fault
      // projection GFs are registered with the primary volume collection.
      // Without this the dynamic slip / slip_rate / traction / state
      // channels never reach `volume.vtkhdf` (only the static SAFS
      // params registered via `SetFaultParamsSpatial` do).  Must be
      // called BEFORE `InitFaultOutputBP5` (the registration happens
      // inside that call).
      pv_out->SetRegisterFaultProjectionsInVolumePV(true);

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
            pv_lsw_mu_s(i) = d.lsw_mu_s;
            pv_lsw_mu_d(i) = d.lsw_mu_d;
            pv_lsw_d_c (i) = d.lsw_d_c;
            const real_t ad = (nuc_params.amplitude_dip.Size()    == num_fault_total)
                              ? nuc_params.amplitude_dip(i)    : 0.0;
            const real_t as = (nuc_params.amplitude_strike.Size() == num_fault_total)
                              ? nuc_params.amplitude_strike(i) : 0.0;
            pv_nuc_amplitude(i) = std::sqrt(ad*ad + as*as);
            pv_nuc_radial   (i) = (nuc_params.radial.Size() == num_fault_total)
                                  ? nuc_params.radial(i)   : 0.0;
            pv_sig_n_init   (i) = d.sigma_n_corr;
            pv_tau1_init    (i) = d.tau1_corr;
            pv_tau2_init    (i) = d.tau2_corr;
         }
      }
      // Parity Phase 6: publish all 8 SAFS static fields via
      // SetFaultParamsSpatial on seas::ParaViewOutput.  Implementation
      // registers each field with `pv_dc_` (the PRIMARY volume PV
      // collection), so these fields land in `volume.vtkhdf` under
      // their SAFS-semantic names (NOT TPV104's a/Dc/x2/x3).  The
      // companion fault artefact (`fault.vtkhdf`) is populated
      // separately by `WriteFaultSurfaceVTU` (R-001) and carries the
      // dynamic slip / slip_rate / traction / state channels.
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
   // R-001: scratch buffers consumed by WriteFaultSurfaceVTU (ADER
   // one-shot: _k4 copies of the live values so the VTU delta channel
   // reads zero, mirroring drivers/tpv205_driver.cpp:2252-2259).  The
   // four `local_a / local_Dc / local_x2 / local_x3` slots in the
   // TPV-shaped signature stay empty — SAFS-semantic static params are
   // already registered with `pv_dc_` via `SetFaultParamsSpatial` and
   // appear in `volume.vtkhdf`.
   Vector pv_local_slip_rate_k4, pv_local_traction_k4;
   Vector pv_local_normal_stress_k4;
   const Vector pv_empty_param;  // unused TPV-shaped param slots
   if (pv_out)
   {
      pv_local_slip.SetSize(2 * num_fault_total);
      pv_local_slip_rate.SetSize(2 * num_fault_total);
      pv_local_traction.SetSize(2 * num_fault_total);
      pv_local_state.SetSize(num_fault_total);
      pv_local_normal_stress.SetSize(num_fault_total);
      pv_local_slip_rate_k4.SetSize(2 * num_fault_total);
      pv_local_traction_k4.SetSize(2 * num_fault_total);
      pv_local_normal_stress_k4.SetSize(num_fault_total);
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
   // 19. Sub-step iterator (TPV205 LSW closed form).  Phase N: only the
   //     plain LSW path is used in this driver — the forced-rupture
   //     mode toggle is gone with the dispatch flip above.
   // -----------------------------------------------------------------
   Tpv205SubStepIterator substep_iterator(fault_flux);
   {
      const int O = std::max(1, cfg.numerics.ader_order);
      std::vector<real_t> deltaT(O, dt / static_cast<real_t>(O));
      std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
      substep_iterator.SetSubSteps(deltaT, weights);
   }

   // Phase N: per-sub-step gradual_overstress accumulator hook.  Closes
   // over nuc_params + dof_data + cfg; fires inside
   // Tpv205SubStepIterator::AdvanceWithSubStepStates (callback overload)
   // BEFORE each sub-step's per-QP friction solve so that the perturbed
   // tau{1,2}_nuc is visible to s.tau{1,2}_total in the LSW path.
   auto nuc_cb = [&nuc_params, &dof_data, &cfg]
                 (real_t t_sub_end, real_t dt_sub)
   {
      if (!cfg.nucleation.enabled) { return; }
      spatial::ApplyGradualOverstressIncrement(
         dof_data, nuc_params,
         cfg.nucleation.gradual_overstress.T_nuc_s,
         t_sub_end, dt_sub);
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
         const real_t delta_norm = std::sqrt(d.slip1 * d.slip1
                                             + d.slip2 * d.slip2);
         // Phase N: spatial driver uses plain LSW.
         pv_local_state(i) = mfem::seas::LSWFrictionCoefficient_TPV205(
                                delta_norm,
                                d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c);
         (void)time;
         pv_local_normal_stress(i)     = d.sigma_n_corr;

         // R-001: ADER one-shot — no stage-4 distinct from averaged
         // DOFData (mirror of tpv205_driver.cpp:2252-2259).  Writing
         // the same values into the `_k4` buffers makes the VTU delta
         // channel read exactly zero.
         pv_local_slip_rate_k4(2 * i + 0) = d.V1;
         pv_local_slip_rate_k4(2 * i + 1) = d.V2;
         pv_local_traction_k4(2 * i + 0)  = d.tau1_corr;
         pv_local_traction_k4(2 * i + 1)  = d.tau2_corr;
         pv_local_normal_stress_k4(i)     = d.sigma_n_corr;
      }

      // Parity Phase 5: refresh the volume velocity field directly from
      // Q (byNODES so VX/VY/VZ are contiguous 3*ndof block).
      if (pv_out->GetVolumeSaveEnabled() && pv_vel_gf)
      {
         std::memcpy(pv_vel_gf->GetData(),
                     Q.GetData() + VX * ndof_total,
                     3 * ndof_total * sizeof(real_t));
      }

      // R-002: refresh fault projections unconditionally — they feed
      // both the volume PV (when registered via R-003) AND the fault
      // VTU/VTKHDF writer (R-001).  Only the volume save itself
      // (`ForceSave` → `pv_dc_->Save()`) is gated on the user toggle.
      pv_out->UpdateFaultFieldsBP5(pv_local_slip, pv_local_slip_rate,
                                   pv_local_traction, pv_local_state,
                                   pv_local_normal_stress);
      if (pv_out->GetVolumeSaveEnabled())
      {
         pv_out->ForceSave(step_num, time);
      }
      pv_out->CommitSchedule(time, V_max);

      // R-001: write the fault-surface VTU/VTKHDF artefact (the file
      // `verify_spatial_dyn_smoke_safs.py` opens and the plan
      // §Acceptance Criteria L491 requires).  Independent of volume
      // save — fault is its own collection.  The four TPV-shaped
      // param slots (`local_a / local_Dc / local_x2 / local_x3`) stay
      // empty: SAFS-semantic static params are already published to
      // `volume.vtkhdf` via `SetFaultParamsSpatial`.
      pv_out->WriteFaultSurfaceVTU(
         cfg.output.output_dir, step_num, time, rank, nprocs,
         pv_local_slip, pv_local_slip_rate, pv_local_traction,
         pv_local_state, pv_local_normal_stress,
         pv_empty_param, pv_empty_param, pv_empty_param, pv_empty_param,
         pv_local_slip_rate_k4, pv_local_traction_k4,
         pv_local_normal_stress_k4);
   };

   if (restart_prefix.empty())
   {
      paraview_write(0, cfg.time.t_initial, 0.0);
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
   for (int step = step0; step < nsteps; ++step)
   {
      const real_t dt_step = std::min(dt_now, cfg.time.tfinal - t);
      if (dt_step <= 0.0) { break; }
      wave.SetTime(t);

      AdvanceADERWithSubStep_Spatial(wave, substep_iterator, dof_data,
                                     fault_coords, Q, dt_step,
                                     cfg.numerics.ader_order, t, Q_new,
                                     nuc_cb);
      Q.Swap(Q_new);
      t += dt_step;
      last_completed_step = step + 1;

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

      // DIAG (env-gated, debug only): rupture-area proxy + onset localizer.
      // SEAS_DIAG_BLOWUP=1 enables.  Prints (a) global count of fault DOFs
      // with V>0.5 m/s (rupturing-area proxy: growing => propagating front,
      // ~constant => only the forced patch) + global max slip, and (b) at
      // onset (V_max>10) the owning rank dumps the argmax DOF's full state
      // so we can see WHAT runs away first (sigma_n? shear? slip?).
      if (std::getenv("SEAS_DIAG_BLOWUP")
          && (step % 100 == 0 || V_max_step > 10.0))
      {
         long long n_rup_local = 0;
         real_t    maxslip_local = 0.0;
         int       argmax_local = -1;
         real_t    vloc = -1.0;
         for (int i = 0; i < num_fault_total; ++i)
         {
            const DOFData &d = dof_data[i];
            if (d.slip_rate > 0.5) { n_rup_local++; }
            const real_t s = std::sqrt(d.slip1 * d.slip1 + d.slip2 * d.slip2);
            if (s > maxslip_local) { maxslip_local = s; }
            if (d.slip_rate > vloc) { vloc = d.slip_rate; argmax_local = i; }
         }
         long long n_rup_g = n_rup_local;
         real_t    maxslip_g = maxslip_local;
         struct { double v; int r; } in_{V_max_local, rank}, out_{V_max_local, rank};
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&n_rup_local, &n_rup_g, 1, MPI_LONG_LONG, MPI_SUM, comm);
         MPI_Allreduce(&maxslip_local, &maxslip_g, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
         MPI_Allreduce(&in_, &out_, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
#endif
         if (rank == 0)
         {
            std::cout << "[DIAG] step " << step << " t=" << t
                      << " V_max=" << V_max_step
                      << " n_rupturing(V>0.5)=" << n_rup_g
                      << " max_slip=" << maxslip_g << " m\n";
         }
         if (V_max_step > 10.0 && rank == out_.r && argmax_local >= 0)
         {
            const DOFData &d = dof_data[argmax_local];
            std::cout << "[DIAG-ONSET] rank " << rank << " dof " << argmax_local
                      << " xyz=(" << dof_coords_3d(3 * argmax_local) << ","
                      << dof_coords_3d(3 * argmax_local + 1) << ","
                      << dof_coords_3d(3 * argmax_local + 2) << ")"
                      << " V=" << d.slip_rate << " V1=" << d.V1
                      << " V2=" << d.V2 << " slip1=" << d.slip1
                      << " slip2=" << d.slip2
                      << " sigma_n_corr=" << d.sigma_n_corr
                      << " tau1_corr=" << d.tau1_corr
                      << " tau2_corr=" << d.tau2_corr
                      << " tau1_nuc=" << d.tau1_nuc
                      << " tau2_nuc=" << d.tau2_nuc << "\n";
         }
      }

      paraview_write(step + 1, t, V_max_step);

      // ------------------------------------------------------------------
      // MPI shared-fault consistency tripwire.  Follows tpv104_driver.cpp
      // (:2636-2639) and tpv205_driver.cpp (:2442-2445), both of which call
      // this once at step 0; the spatial driver had DROPPED it entirely.
      //
      // Extended past TPV104's step-0-only call to fire through the whole
      // gradual_overstress nucleation window: on the CURVILINEAR SAFS fault
      // the per-QP dip/strike frame (FaultBasis::ComputeOrientedFrame, via
      // sign(n_raw . ref_normal) with ref_normal=(0,-1,0)) is ill-conditioned
      // where the local normal is ~perpendicular to ref_normal, so the two
      // ranks owning a shared face can pick OPPOSITE strike directions.  The
      // fixed-sign nucleation increment (+F.dtau into tau2_nuc) then forces
      // the two sides of a shared face in opposite physical directions.  The
      // routine pairs same-physical-QP records across ranks and ABORTS with a
      // field/QP diagnostic if the 8 evolved fields differ > tol (so the
      // first offending shared QP inside the nucleation patch is pinpointed).
      // TPV104's planar fault has |n_raw . ref_normal| = 1 everywhere, so it
      // never flips and always passes this check.
      if (step == 0
          || (cfg.nucleation.enabled
              && t <= cfg.nucleation.gradual_overstress.T_nuc_s
              && (step % 100 == 0)))
      {
         // SEAS_R101_NONFATAL (env-gated, OFF by default => zero behaviour
         // change for TPV/BP5/production SAFS): downgrade this shared-fault
         // consistency tripwire from a hard MFEM_ABORT to a per-check rank-0
         // diagnostic.  Lets a run with a KNOWN cross-rank DOFData divergence
         // (the curvilinear-fault strike-frame flip) proceed to completion so
         // we can test whether a candidate fix (e.g. the fault-normal
         // canonicalization change in fault_basis.hpp) keeps worst_rel bounded
         // AND avoids the V_max blow-up.  The routine still runs every MPI
         // collective; abort_on_fail=false just suppresses the abort and writes
         // worst_rel / worst_field instead (wave_operator.inl:5826).
         const bool r101_nonfatal =
            (std::getenv("SEAS_R101_NONFATAL") != nullptr);
         double r101_worst_rel = 0.0;
         int    r101_worst_field = -1;
         wave.VerifySharedFaultDOFDataConsistency(
            /*tol=*/1e-10, &r101_worst_rel, &r101_worst_field,
            /*abort_on_fail=*/!r101_nonfatal);
         if (r101_nonfatal && rank == 0)
         {
            std::cout << "  [R-101 nonfatal] step " << step << " t=" << t
                      << "s  worst_rel=" << r101_worst_rel
                      << "  worst_field=" << r101_worst_field
                      << (r101_worst_rel > 1e-10 ? "  (DIVERGED)"
                                                 : "  (ok)")
                      << std::endl;
         }
      }

      // NaN tripwire — follows tpv104_driver.cpp:2641-2659 /
      // tpv205_driver.cpp:2447-2464 (also dropped by the spatial driver).
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
#ifdef MFEM_USE_MPI
            MPI_Finalize();
#endif
            return 1;
         }
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

      // R-008: also print on the step that exhausts `tfinal` so the
      // "final-step" log fires even when the loop exits via the early
      // `dt_step <= 0.0` break (when `nsteps = ceil(tfinal/dt)` over-
      // shoots, the last advance lands on a step != nsteps - 1).
      const bool reached_tfinal =
         (t + std::numeric_limits<real_t>::epsilon() >= cfg.time.tfinal);
      if (rank == 0 && (step % 100 == 0
                        || step == nsteps - 1
                        || reached_tfinal))
      {
         std::cout << "step " << step << "/" << nsteps
                   << "  t = " << t << " s"
                   << "  V_max = " << V_max_step << " m/s\n";
      }
   }

   // -----------------------------------------------------------------
   // 21. Final checkpoint.
   // -----------------------------------------------------------------
   // R-011: skip when the last in-loop checkpoint already covered
   // `last_completed_step` — the two writes would be byte-identical
   // and waste Lustre metadata ops on production.
   const bool already_checkpointed_final =
      (cfg.output.checkpoint_every_steps > 0)
      && (last_completed_step > step0)
      && (last_completed_step % cfg.output.checkpoint_every_steps == 0);
   if (cfg.output.checkpoint_every_steps > 0 && !already_checkpointed_final)
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
   else if (already_checkpointed_final && rank == 0)
   {
      std::cout << "[checkpoint] last in-loop checkpoint already covers "
                << "step " << last_completed_step
                << "; skipping redundant final write.\n";
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
