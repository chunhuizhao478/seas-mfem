// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_PARAVIEW_OUTPUT_HPP
#define MFEM_SEAS_PARAVIEW_OUTPUT_HPP

#include "mfem.hpp"
#include "../config/bp5_params.hpp"
#include "fault_vtu_binary.hpp"   // vtu::LocalFaultPack, GatherFaultPackToRoot,
                                   // WriteFaultPackVTU
#ifdef MFEM_USE_HDF5
// R-401: pull `mesh/vtkhdf.hpp` in directly so `MFEM_PARALLEL_HDF5` is
// visible to `DefaultVolumeOutputMode()`.  Without this include the
// macro is never set in this translation unit and the ParMesh branch
// silently fell back to per-rank VTU on every Frontera run.
#include "mesh/vtkhdf.hpp"
#include "fault_vtkhdf_writer.hpp" // vtkhdf::FaultHDFState, WriteFaultPackHdf
#endif

#include <string>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <sys/stat.h>
#include <vector>
#include <type_traits>
#include <set>

namespace mfem
{
namespace seas
{

// Helper to select GridFunction vs ParGridFunction based on MeshType.
template <typename MeshType> struct GFType
{ using type = GridFunction; using FESType = FiniteElementSpace; };

#ifdef MFEM_USE_MPI
template <> struct GFType<ParMesh>
{ using type = ParGridFunction; using FESType = ParFiniteElementSpace; };
#endif

/// @brief Manages ParaView-compatible visualization output for SEAS simulations.
///
/// Outputs domain fields (displacement) directly and fault fields (slip,
/// slip rate, traction, state variable) as piecewise-constant (L2 p=0)
/// GridFunctions on the domain mesh.  Fault values are mapped from owned
/// fault DOF vectors to the two volume elements adjacent to each fault face.
///
/// Output frequency is adaptive, keyed to V_max (slip rate).  Defaults:
///   - coseismic  (V > 1e-3 m/s)  : every 0.01 s
///   - nucleation (V > 1e-6)       : every 1.0 s
///   - interseismic                : every 1 year
/// Thresholds and intervals are runtime-configurable via GetSchedule()
/// (see AdaptiveSchedule).  A hysteresis_factor > 1 suppresses regime
/// chattering near the thresholds.
///
/// @tparam MeshType  Mesh (serial) or ParMesh (parallel)
template <typename MeshType = Mesh>
class ParaViewOutput
{
   using GF  = typename GFType<MeshType>::type;
   using FES = typename GFType<MeshType>::FESType;

   static constexpr real_t kOutputTimeTolerance = 0.99;

public:
   /// Selects the fault-surface VTU writer back end.
   ///
   ///   `Vtu`  - Phase 1 single per-cycle binary VTU on rank 0
   ///            (gather-to-rank-0; all ranks must call collectively).
   ///   `Hdf5` - Phase 2 single-file VTKHDF for the entire run.  Requires
   ///            `MFEM_USE_HDF5=YES` at build time.  Phase 1 leaves this
   ///            enum value defined for forward compatibility but
   ///            `WriteFaultSurfaceVTU` aborts if it is selected before
   ///            Phase 2 ships.
   ///
   /// Default is `Vtu`.  Phase 2b will switch the default to `Hdf5` on
   /// builds that define `MFEM_USE_HDF5`.
   enum class FaultOutputMode { Vtu, Hdf5 };

   /// Selects the volume-PV writer back end (Phase 6).
   ///
   ///   `Vtu`  - legacy `mfem::ParaViewDataCollection` (per-rank
   ///            ASCII / binary VTU + PVTU + PVD).  ParaView 5.4+.
   ///   `Hdf5` - `mfem::ParaViewHDFDataCollection` (single
   ///            `<prefix>/<collection_name>.vtkhdf` per run).
   ///            Requires `MFEM_USE_HDF5=YES` AND, for `ParMesh`,
   ///            `MFEM_PARALLEL_HDF5=YES` (the parallel HDF5 build).
   ///            ParaView 5.11+.
   ///
   /// Default is `DefaultVolumeOutputMode()` — picks `Hdf5` when the
   /// build can drive it (HDF5 + parallel-HDF5 for `ParMesh`; HDF5
   /// alone for serial `Mesh`), `Vtu` otherwise.
   enum class VolumeOutputMode { Vtu, Hdf5 };

   /// R-304 / Phase 6.1: compile-time selector that picks `Hdf5` on
   /// any build that can drive it.  R-401: the previous "fall back to
   /// Vtu when MFEM_PARALLEL_HDF5 is missing" branch silently disabled
   /// Phase 6 on every Frontera ParMesh run; that fallback has been
   /// replaced with a hard compile error.  Users who genuinely need
   /// the legacy per-rank VTU path must pass `VolumeOutputMode::Vtu`
   /// explicitly to the constructor.
   static constexpr VolumeOutputMode DefaultVolumeOutputMode()
   {
#if defined(MFEM_USE_HDF5)
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
#  if defined(MFEM_PARALLEL_HDF5)
         return VolumeOutputMode::Hdf5;
#  else
         // No silent VTU fallback (R-401).  Build needs MFEM with
         // parallel HDF5 enabled (HDF5_HAVE_PARALLEL=ON) to use the
         // default volume mode on a ParMesh.  Pass
         // VolumeOutputMode::Vtu explicitly if you really want the
         // legacy per-rank path on a serial-HDF5 build.
         static_assert(!std::is_same_v<MeshType, ParMesh>,
                       "ParaViewOutput<ParMesh> default volume mode "
                       "requires MFEM with MFEM_PARALLEL_HDF5 (parallel "
                       "HDF5).  Rebuild MFEM with HDF5_HAVE_PARALLEL "
                       "enabled, or pass VolumeOutputMode::Vtu "
                       "explicitly to the constructor.");
#  endif
      }
      else
      {
         return VolumeOutputMode::Hdf5;
      }
#else
      return VolumeOutputMode::Vtu;
#endif
   }

   /// Step-based output interval (write every N steps).  Set to 0 to
   /// use the adaptive V_max or fixed-dt schedule.
   int output_every_n_steps = 0;

   /// Fixed time interval between writes (seconds).  Set to 0 to use
   /// the adaptive V_max schedule.  Takes precedence over V_max schedule
   /// but not over step-based interval.
   real_t fixed_dt = 0.0;

   /// Runtime-configurable adaptive-schedule parameters.  Defaults match
   /// the legacy hardcoded OutputInterval thresholds/intervals and have
   /// hysteresis_factor = 1.0 (no hysteresis), so existing callers see
   /// bit-identical behavior.
   ///
   /// Regime indices: 0 = interseismic, 1 = nucleation, 2 = coseismic.
   struct AdaptiveSchedule
   {
      // V thresholds (m/s).  Coseismic enter when V exceeds v_coseismic;
      // nucleation enter when V exceeds v_nucleation.
      real_t v_coseismic       = 1e-3;
      real_t v_nucleation      = 1e-6;
      // Hysteresis factor (>= 1).  Leaving a regime requires V to drop
      // below (threshold / hysteresis_factor).  1.0 disables hysteresis.
      real_t hysteresis_factor = 1.0;
      // Output intervals (seconds) used in each regime.
      real_t dt_coseismic      = 0.01;
      real_t dt_nucleation     = 1.0;
      real_t dt_interseismic   = 1.0 * BP5Params::seconds_per_year;

      /// Phase 3: Soft cap on the total number of fault snapshots in a
      /// run.  0 = uncapped (default, preserves pre-Phase-3 behaviour).
      /// >0 = inflate `dt_interseismic` so the projected remaining
      /// snapshots fit under the cap.  The cap is a SOFT bound on
      /// interseismic-regime writes only — coseismic and nucleation
      /// regimes write at their physically-meaningful cadences and may
      /// push the final count slightly above K (by at most one write
      /// per event, see plan §Phase 3 step 1).  Plan §Phase 3 step 3
      /// recommends 5000 for production BP5.
      int max_total_snapshots = 0;

      /// Phase 3 R-002: latched on the first call to
      /// `RecomputeIntervalForCap` that hits the budget-exhausted branch
      /// (`remaining_budget <= 0`).  Used to print the plan-mandated
      /// one-line rank-0 warning exactly once per run.  `mutable`
      /// because the warning fires inside a `const` method.
      mutable bool cap_exhausted_warned_ = false;

      /// Return the output interval (seconds) for a given regime.
      /// V is accepted for signature symmetry with NextRegime but is
      /// unused — the cadence is stable within a regime.
      real_t Interval(real_t /*V*/, int regime) const
      {
         switch (regime)
         {
            case 2:  return dt_coseismic;
            case 1:  return dt_nucleation;
            default: return dt_interseismic;
         }
      }

      /// State machine mapping (V, prev_regime) -> new regime.
      ///
      /// Entry comparisons are STRICT > (not >=) to preserve byte
      /// compatibility with the legacy OutputInterval which used
      /// V_max > 1e-3 / V_max > 1e-6.  Exit comparisons are strict <,
      /// so an exact-threshold V stays in its current regime.
      ///
      /// Non-finite V (NaN/Inf) is treated as coseismic, so we keep
      /// dense output during a blowup.
      int NextRegime(real_t V, int prev) const
      {
         const real_t V_co_enter = v_coseismic;
         const real_t V_co_exit  = v_coseismic  / hysteresis_factor;
         const real_t V_nu_enter = v_nucleation;
         const real_t V_nu_exit  = v_nucleation / hysteresis_factor;

         if (!std::isfinite(V)) { return 2; }

         switch (prev)
         {
            case 0: // interseismic
               if (V > V_co_enter) { return 2; }
               if (V > V_nu_enter) { return 1; }
               return 0;
            case 1: // nucleation
               if (V > V_co_enter) { return 2; }
               if (V < V_nu_exit)  { return 0; }
               return 1;
            case 2: // coseismic
               if (V < V_co_exit)
               {
                  // Drop-past-threshold semantics (intentional):
                  // the target uses the STRICT ENTRY threshold, so a V
                  // landing in [V_nu_exit, V_nu_enter] goes straight to
                  // regime 0 (interseismic) — we do NOT transfer any
                  // accumulated nucleation-hysteresis credit from the
                  // outbound path.  Rationale: after a completed
                  // rupture, V is usually already deep in interseismic
                  // territory; returning to nucleation only when V is
                  // unambiguously above the entry threshold avoids
                  // over-sampling post-rupture ring-down.
                  if (V > V_nu_enter) { return 1; }
                  return 0;
               }
               return 2;
            default: return 0;
         }
      }

      /// Phase 3: cap-aware interseismic interval.
      ///
      /// Returns an interval (seconds) that, when used as the
      /// interseismic dt going forward, keeps the total snapshot count
      /// for the run under @a max_total_snapshots.  Coseismic and
      /// nucleation regimes are unaffected — their cadences are
      /// physically meaningful and the cap can only inflate
      /// interseismic spacing.
      ///
      /// Algorithm (plan §Phase 3 step 2):
      ///   remaining_budget = K - k
      ///   if remaining_budget <= 0:
      ///       dt_inter_new = T - t       // emit at most one more
      ///   else:
      ///       dt_inter_new = max(dt_inter, (T - t) / remaining_budget)
      ///
      /// @param time_to_end       Seconds remaining until the run ends.
      ///                          Negative or zero ⇒ no further writes.
      /// @param snapshots_so_far  Snapshots already written this run.
      /// @return                  The (possibly inflated) dt to use for
      ///                          subsequent interseismic writes.
      ///                          Returns dt_interseismic unchanged if
      ///                          max_total_snapshots == 0.
      real_t RecomputeIntervalForCap(real_t time_to_end,
                                     int snapshots_so_far) const
      {
         if (max_total_snapshots <= 0) { return dt_interseismic; }
         const int remaining_budget = max_total_snapshots - snapshots_so_far;
         if (remaining_budget <= 0)
         {
            // R-002: plan §Phase 3 edge case 1 — print one rank-0
            // warning the first time we hit the exhausted branch, so a
            // user who reads the run log sees that coseismic events
            // pushed past the cap.  `mfem::out` is rank-0-only on
            // parallel builds and a no-op sink on serial.
            if (!cap_exhausted_warned_)
            {
               mfem::out << "ParaViewOutput: max_total_snapshots="
                         << max_total_snapshots << " exceeded by "
                         "coseismic events (snapshots_so_far="
                         << snapshots_so_far << ").  Spreading "
                         "remaining " << time_to_end << " s into one "
                         "final interseismic write.\n";
               cap_exhausted_warned_ = true;
            }
            // R-005: when `time_to_end > 0` we honour the plan's "emit
            // at most one more" by returning T - t (the writer will
            // fire once when `time` accumulates that much past
            // `last_write_time_`).  When the simulation overruns
            // `total_run_time_` so `time_to_end <= 0`, return a hard
            // ceiling so further writes are LOCKED OUT rather than
            // re-enabling at user dt.  R-109: the caller multiplies
            // this by `kOutputTimeTolerance = 0.99` and compares to
            // `time - last_write_time_`; on `MFEM_USE_SINGLE=YES` builds
            // `numeric_limits<float>::max() * 0.99` is precariously
            // close to the float overflow threshold.  Halving the
            // sentinel keeps the multiply finite on float while still
            // being trivially larger than any plausible `time -
            // last_write_time_`.
            return time_to_end > 0.0
                   ? time_to_end
                   : (std::numeric_limits<real_t>::max() / real_t(2));
         }
         const real_t projected = (time_to_end > 0.0)
                                  ? (time_to_end / real_t(remaining_budget))
                                  : dt_interseismic;
         return std::max(dt_interseismic, projected);
      }

      /// One-shot invariant check.  Call once after applying any CLI
      /// overrides (NOT inside NextRegime / Interval — those are
      /// hot-path).  Aborts with a readable message on bad inputs.
      void Validate() const
      {
         MFEM_VERIFY(hysteresis_factor >= 1.0,
                     "AdaptiveSchedule: hysteresis_factor must be >= 1.0, got "
                     << hysteresis_factor);
         MFEM_VERIFY(v_coseismic > v_nucleation && v_nucleation > 0.0,
                     "AdaptiveSchedule: require v_coseismic (" << v_coseismic
                     << ") > v_nucleation (" << v_nucleation << ") > 0");
         MFEM_VERIFY(dt_coseismic > 0.0 && dt_nucleation > 0.0 &&
                     dt_interseismic > 0.0,
                     "AdaptiveSchedule: all dt_* must be positive");
         MFEM_VERIFY(max_total_snapshots >= 0,
                     "AdaptiveSchedule: max_total_snapshots must be >= 0 "
                     "(0 means uncapped), got " << max_total_snapshots);
      }
   };

   /// Mutable accessor so drivers can override thresholds/intervals at
   /// configuration time.  Call Validate() after tweaking.
   AdaptiveSchedule &GetSchedule() { return adaptive_; }
   const AdaptiveSchedule &GetSchedule() const { return adaptive_; }

   /// Phase 3: total simulation time (seconds) for the current run.
   /// Used by the snapshot-cap mechanism to compute `time_to_end`.
   /// Drivers must call this once before the time-stepping loop (e.g.
   /// `pv->SetTotalRunTime(tfinal)`).  When `max_total_snapshots > 0`
   /// AND `total_run_time_ <= 0`, the next call to a cap-aware schedule
   /// path aborts with `MFEM_VERIFY` rather than silently disabling the
   /// cap (R-001).
   ///
   /// Plan §Phase 3 follow-up TODO: wire `--paraview-max-snapshots N`
   /// + `--paraview-coseismic-dt`, `--paraview-nucleation-dt`,
   /// `--paraview-interseismic-dt` flags into the four drivers'
   /// argument parsers.  The C++ API is at:
   ///    pv->GetSchedule().max_total_snapshots = N;
   ///    pv->GetSchedule().dt_coseismic = X;        (etc.)
   ///    pv->SetTotalRunTime(tfinal);
   /// Blocked locally by the pre-existing MPIContext::GetComm /
   /// FaultScatter source skew that breaks the full-driver build.
   void SetTotalRunTime(real_t tfinal) { total_run_time_ = tfinal; }
   real_t GetTotalRunTime() const { return total_run_time_; }

   /// Phase 3: read-only counter of fault snapshots committed by Save /
   /// ShouldWrite / ForceSave so far this run.
   int GetTotalSnapshotsWritten() const { return total_snapshots_written_; }

   /// @brief Construct the ParaView output manager.
   ///
   /// R-303: new `collection_name` and `mode` parameters appended at
   /// the end with defaults so the existing 3-arg call sites keep
   /// compiling unchanged.  R-305: default collection name = "volume"
   /// for cross-driver uniformity (mirrors the fault path's "fault"
   /// — formerly "fault_surface", renamed in
   /// PLAN_split_bulk_solutions_2026-05-12).
   ///
   /// @param prefix          Output directory.  The volume writer
   ///                        emits either per-rank VTU files into
   ///                        `<prefix>/` (Vtu mode) or a single
   ///                        `<prefix>/<collection_name>.vtkhdf`
   ///                        (Hdf5 mode).
   /// @param mesh            Mesh / ParMesh to attach to the volume
   ///                        writer.  Mid-run mesh swap is NOT
   ///                        supported on the Hdf5 path (the file
   ///                        is keyed on the first mesh).
   /// @param order           Polynomial order for `SetLevelsOfDetail`.
   /// @param collection_name Volume-side collection name (basename
   ///                        of the .vtkhdf file in Hdf5 mode).
   /// @param mode            `VolumeOutputMode::Hdf5` or
   ///                        `VolumeOutputMode::Vtu`; default is
   ///                        `DefaultVolumeOutputMode()`.
   ///
   /// Re-running with the same `output_dir` overwrites any existing
   /// `.vtkhdf` (Hdf5 mode) — restart-mode is a planned follow-up.
   ParaViewOutput(const std::string &prefix,
                  MeshType &mesh,
                  int order,
                  const std::string &collection_name = "volume",
                  VolumeOutputMode mode = DefaultVolumeOutputMode())
      : mesh_(mesh),
        order_(order),
        last_write_time_(-1e30),
        volume_output_mode_(mode)
   {
      if (volume_output_mode_ == VolumeOutputMode::Hdf5)
      {
#ifdef MFEM_USE_HDF5
         pv_dc_ = std::make_unique<ParaViewHDFDataCollection>(
                     collection_name, &mesh);
#else
         MFEM_ABORT("VolumeOutputMode::Hdf5 selected on a build "
                    "without MFEM_USE_HDF5=YES");
#endif
      }
      else
      {
         pv_dc_ = std::make_unique<ParaViewDataCollection>(
                     collection_name, &mesh);
      }
      pv_dc_->SetPrefixPath(prefix);
      pv_dc_->SetDataFormat(VTKFormat::BINARY);
      pv_dc_->SetHighOrderOutput(true);
      pv_dc_->SetLevelsOfDetail(order);
   }

   /// R-302: getter only — no setter.  The mode is locked at
   /// construction (a setter would be unreachable from external
   /// callers anyway).
   VolumeOutputMode GetVolumeOutputMode() const { return volume_output_mode_; }

   // ---------------------------------------------------------------
   //  Domain field registration
   // ---------------------------------------------------------------

   /// Register a domain-mesh GridFunction for output (non-owning).
   ///
   /// R-311: forwards through the abstract `pv_dc_` (a
   /// `ParaViewDataCollectionBase*`).  `RegisterField` is virtual on
   /// the GRANDPARENT `mfem::DataCollection`; both
   /// `ParaViewDataCollection` and `ParaViewHDFDataCollection`
   /// inherit its body (direct insert into `field_map`).  Virtual
   /// dispatch reaches the same implementation regardless of which
   /// concrete subclass is active.
   void RegisterDomainField(const std::string &name, GF *gf)
   {
      pv_dc_->RegisterField(name, gf);
   }

   // ---------------------------------------------------------------
   //  Fault-surface VTU field filter
   // ---------------------------------------------------------------

   /// @brief Restrict which CellData fields WriteFaultSurfaceVTU emits.
   ///
   /// Default (empty set) emits ALL 12 standard fields plus the 5 `_k4`
   /// fields when stage-4 buffers are supplied, matching pre-feature
   /// behaviour.  When `fields` is non-empty, WriteFaultSurfaceVTU
   /// emits only the CellData arrays whose names appear in `fields`
   /// (comparison is case-sensitive, exact match).  Unknown names in
   /// `fields` are silently ignored — the resulting VTU simply has no
   /// such field.  Points / Cells / connectivity are always written.
   ///
   /// Use cases:
   ///  - Large-scale BP5 runs where only `slip_rate_strike` is needed
   ///    for visualisation ⇒ ~12× per-file size reduction vs the full
   ///    12-field default, without breaking backward compatibility.
   ///  - Debugging workflows where a single field is inspected.
   ///
   /// The filter does NOT affect `UpdateFaultFieldsBP5` (volume-PVD
   /// L2-p0 fields) — those remain fully populated.
   void SetFaultVTUFields(const std::set<std::string> &fields)
   { fault_vtu_fields_ = fields; }

   /// @brief Current filter set; empty ⇒ all fields.
   const std::set<std::string> &GetFaultVTUFields() const
   { return fault_vtu_fields_; }

   /// Select the fault-surface output back end.  See `FaultOutputMode`.
   void SetFaultOutputMode(FaultOutputMode mode) { output_mode_ = mode; }
   FaultOutputMode GetFaultOutputMode() const { return output_mode_; }

   /// Toggle the legacy per-rank ASCII fault-surface VTU writer.  When
   /// enabled, `WriteFaultSurfaceVTU` reverts to the pre-Phase-1 layout:
   /// one ASCII `fault_surface_r{r}_c{c}.vtu` per rank per cycle plus a
   /// per-cycle `fault_surface_c{c}.pvtu` index.  Useful for debugging
   /// and for tests whose parsers depend on the ASCII regex layout.
   /// Default `false`.  Has no effect when `GetFaultOutputMode() ==
   /// Hdf5` (HDF5 always emits a single file).
   void SetLegacyAsciiVTU(bool enable) { legacy_ascii_ = enable; }
   bool GetLegacyAsciiVTU() const { return legacy_ascii_; }

#ifdef MFEM_USE_HDF5
   /// @brief R-307: detect H5Z-ZFP plugin availability at runtime.
   ///
   /// Plan §Phase 2d.3 edge case mandates an actionable error message
   /// (rather than HDF5's default "required filter is not registered"
   /// error stack) when the user enables ZFP without a discoverable
   /// plugin.  We check `HDF5_PLUGIN_PATH` for `libh5zzfp.{so,dylib}`
   /// and abort with a one-line MFEM_VERIFY message naming the env
   /// var.  This is implemented at the seas wrapper (NOT in MFEM) so
   /// MFEM remains generic.
   static void ProbeH5ZZfpPluginOrAbort()
   {
      static bool probed = false;
      if (probed) { return; }
      probed = true;
      const char *plugin_path = std::getenv("HDF5_PLUGIN_PATH");
      auto exists = [](const std::string &p) -> bool
      {
         struct stat st;
         return ::stat(p.c_str(), &st) == 0;
      };
      // Search candidates: HDF5_PLUGIN_PATH (may be a colon-separated
      // list per HDF5 docs) plus the conventional install location at
      // ${HDF5_PLUGIN_PATH}/libh5zzfp.{so,dylib}.
      bool found = false;
      auto try_dir = [&](const std::string &dir)
      {
         if (dir.empty()) { return; }
         for (const char *name : {"libh5zzfp.dylib", "libh5zzfp.so",
                                  "libh5zzfp.1.dylib", "libh5zzfp.1.so"})
         {
            if (exists(dir + "/" + name)) { found = true; return; }
         }
      };
      if (plugin_path)
      {
         std::string acc;
         for (const char *p = plugin_path; ; ++p)
         {
            if (*p == ':' || *p == '\0')
            {
               try_dir(acc);
               if (found || *p == '\0') { break; }
               acc.clear();
            }
            else { acc += *p; }
         }
      }
      MFEM_VERIFY(found,
                  "H5Z-ZFP plugin not found at runtime.  Set "
                  "HDF5_PLUGIN_PATH to the directory containing "
                  "libh5zzfp.{so,dylib} and re-run.  Current value: "
                  << (plugin_path ? plugin_path : "<unset>"));
   }

   /// @brief Phase 2d.3 — choose the chunk filter for the fault VTKHDF
   /// output.
   ///
   /// `Deflate` (default) reproduces Phase 2b lossless behaviour.
   /// `ZfpAccuracy` switches FP datasets to LLNL ZFP (filter id 32013)
   /// in absolute-error accuracy mode at @a param tolerance; integer
   /// connectivity stays lossless via the deflate fallback in
   /// `VTKHDF::EnsureDataset`.
   ///
   /// Must be called BEFORE the first `WriteFaultSurfaceVTU` invocation
   /// in HDF5 mode (the data collection is constructed lazily on the
   /// first save).  Calling it after the first save updates the active
   /// dc but does not retroactively change already-written timesteps.
   ///
   /// @param alg    Filter algorithm.
   /// @param param  Deflate level (0..9) when alg == Deflate, or ZFP
   ///               accuracy tolerance (must be > 0) when alg ==
   ///               ZfpAccuracy.
   void SetFaultHDFCompression(
      ParaViewHDFDataCollection::HDFCompression alg,
      double param)
   {
      fault_hdf_state_.compression_alg   = alg;
      fault_hdf_state_.compression_param = param;
      if (alg == ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy)
      {
         MFEM_VERIFY(param > 0.0,
                     "SetFaultHDFCompression: ZfpAccuracy tolerance must "
                     "be > 0, got " << param);
         // R-307 / plan §Phase 2d.3 edge case "HDF5 plugin not on
         // HDF5_PLUGIN_PATH at run time": probe the env-var + library
         // BEFORE the first Save() so the user gets an actionable
         // message instead of the cryptic HDF5 "required filter is
         // not registered" error stack.  We do this at the seas-
         // wrapper level (NOT in MFEM) so MFEM stays generic.
         ProbeH5ZZfpPluginOrAbort();
      }
      if (fault_hdf_state_.initialised)
      {
         fault_hdf_state_.dc->SetHDFCompression(alg, param);
      }
   }

   /// @brief Phase 2d.3 / Phase 6.2 — choose the chunk filter for the
   /// volume PV path.
   ///
   /// `Deflate` (the default at level 6) reproduces the legacy VTKHDF
   /// lossless behaviour.  `ZfpAccuracy` switches FP datasets
   /// (velocity / sigma / displacement) to LLNL ZFP (filter id 32013)
   /// in absolute-error accuracy mode at @a param tolerance; integer
   /// connectivity stays lossless via the deflate fallback in
   /// `VTKHDF::EnsureDataset`.
   ///
   /// Phase 6.2: when the active volume writer is
   /// `ParaViewHDFDataCollection` (Hdf5 mode), the request actually
   /// applies.  When the active writer is `ParaViewDataCollection`
   /// (legacy Vtu mode — selected explicitly via `--paraview-volume-vtu`
   /// or implicitly on non-HDF5 builds), the call is a one-time
   /// rank-0 warning + no-op (the VTU writer cannot attach HDF5 chunk
   /// filters).
   ///
   /// Must be called BEFORE the first `pv_dc_->Save()` for the filter
   /// to attach to the dataset's filter chain (the Phase 2d.2 MFEM
   /// patch consults `algorithm` only at first-write time per
   /// `mesh/vtkhdf.cpp:EnsureDataset`).
   ///
   /// @param alg    Filter algorithm.
   /// @param param  Deflate level (0..9) when alg == Deflate, or ZFP
   ///               accuracy tolerance (must be > 0) when alg ==
   ///               ZfpAccuracy.
   void SetVolumeHDFCompression(
      ParaViewHDFDataCollection::HDFCompression alg,
      double param)
   {
      if (alg == ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy)
      {
         MFEM_VERIFY(param > 0.0,
                     "SetVolumeHDFCompression: ZfpAccuracy tolerance must "
                     "be > 0, got " << param);
      }
      auto *hdf_dc = GetHDFDC();
      if (hdf_dc)
      {
         // Phase 6.2: real apply on the active HDF collection.
         hdf_dc->SetHDFCompression(alg, param);
         if (alg == ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy)
         {
            // R-307: probe HDF5_PLUGIN_PATH at SetVolumeHDFCompression
            // time so the user gets an actionable error if the H5Z-ZFP
            // plugin is not discoverable, instead of cryptic HDF5
            // stack on first Save().
            ProbeH5ZZfpPluginOrAbort();
         }
         return;
      }
      // VTU active — emit a one-time warn-and-ignore.  R-316: warn for
      // every non-default request.  Only Deflate at the legacy default
      // level (6) is a true no-op; any other Deflate level cannot be
      // honoured by the VTU bulk path either.
      const bool nontrivial =
         (alg == ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy) ||
         (alg == ParaViewHDFDataCollection::HDFCompression::Deflate
          && static_cast<int>(param) != 6);
      if (nontrivial && !volume_hdf_compression_warned_)
      {
         mfem::out
            << "ParaViewOutput::SetVolumeHDFCompression: requested filter "
               "selector ignored — the active volume writer is per-rank "
               "VTU (`mfem::ParaViewDataCollection`), not VTKHDF.  The "
               "request is not applied; bulk output remains lossless.  "
               "Switch to VolumeOutputMode::Hdf5 (the default on "
               "MFEM_USE_HDF5=YES builds) to apply the filter, or pass "
               "--paraview-volume-hdf5 from the driver CLI.\n";
         volume_hdf_compression_warned_ = true;
      }
   }
#endif

   /// @brief Phase 4 — gate the fault-driven volume-PV save.
   ///
   /// When `enabled` is `false`, `ForceSaveImpl` skips the underlying
   /// `pv_.Save()` so the per-cycle volume PVD/VTU is not emitted, but
   /// the fault-surface PVD/VTU bookkeeping continues to fire (the
   /// snapshot counter still bumps; `WriteFaultSurfaceVTU` still
   /// dispatches).  This is the production default for BP5 quasi-
   /// dynamic runs where the volume save is large but rarely opened.
   void SetVolumeSaveEnabled(bool enabled) { emit_volume_save_ = enabled; }
   bool GetVolumeSaveEnabled() const { return emit_volume_save_; }

   /// @brief Phase 4 — set an independent cadence for the volume PV save.
   ///
   /// When @a dt > 0, `pv_.Save()` only fires when the time since the
   /// last volume write meets or exceeds @a dt (subject to the standard
   /// `kOutputTimeTolerance` slack).  When @a dt <= 0 (default), the
   /// volume save fires at every `ForceSaveImpl` call (matching the
   /// fault cadence).  R-315: a positive @a dt implicitly enables the
   /// volume save, overriding any prior `SetVolumeSaveEnabled(false)`,
   /// per plan §Phase 4 edge case "the explicit dt wins".
   void SetVolumePVDt(real_t dt)
   {
      volume_pv_dt_ = dt;
      if (dt > 0.0) { emit_volume_save_ = true; }
   }
   real_t GetVolumePVDt() const { return volume_pv_dt_; }

   /// @brief PLAN_split_bulk_solutions_2026-05-12: opt-in registration
   /// of the 12 L2-p0 fault projection fields with the primary
   /// (kinematics) volume PV collection.
   ///
   /// Default `false` — `InitFaultOutputBP5` still allocates the 12
   /// GridFunctions (so accessors and direct reads keep working) but
   /// does NOT register them with `pv_dc_`, so they are excluded from
   /// `kinematics.vtkhdf`.  Pass `true` BEFORE calling
   /// `InitFaultOutputBP5` to opt back in to the pre-2026-05-12
   /// behaviour where the projections appear alongside velocity /
   /// displacement on the volume mesh.  Rationale: the projections
   /// duplicate the per-DOF data already in `fault.vtkhdf` and were
   /// the largest contributor to `volume.vtkhdf` file size.
   void SetRegisterFaultProjectionsInVolumePV(bool enable)
   { register_fault_projections_in_volume_pv_ = enable; }
   bool GetRegisterFaultProjectionsInVolumePV() const
   { return register_fault_projections_in_volume_pv_; }

   // ---------------------------------------------------------------
   //  BP5 fault field output (2-component tangential + 1 state)
   // ---------------------------------------------------------------

   /// @brief Initialize L2-p0 fault fields on the domain mesh for BP5.
   ///
   /// Creates 7 scalar GridFunctions (slip_dip, slip_strike, slip_rate_dip,
   /// slip_rate_strike, traction_dip, traction_strike, state_variable) as
   /// piecewise-constant on the volume mesh.  Only elements adjacent to
   /// fault faces carry nonzero values.
   ///
   /// @param fault_interior_faces  Interior fault face indices (local faces)
   /// @param fault_shared_faces    Shared fault face indices (for ParMesh)
   /// @param nbf_per_face          Basis functions per fault face (e.g. 3 for p=1 tri)
   void InitFaultOutputBP5(const Array<int> &fault_interior_faces,
                           const Array<int> &fault_shared_faces,
                           int nbf_per_face)
   {
      nbf_per_face_ = nbf_per_face;
      int dim = mesh_.Dimension();

      // L2 p=0 on domain mesh — one DOF per element
      fault_fec_ = std::make_unique<L2_FECollection>(0, dim);
      fault_fes_ = MakeFES(mesh_, fault_fec_.get());

      // Allocate 7 scalar fields
      auto make_gf = [&]() {
         auto gf = std::make_unique<GF>(fault_fes_.get());
         *gf = 0.0;
         return gf;
      };
      fault_slip_dip_        = make_gf();
      fault_slip_strike_     = make_gf();
      fault_slip_rate_dip_   = make_gf();
      fault_slip_rate_strike_= make_gf();
      fault_trac_dip_        = make_gf();
      fault_trac_strike_     = make_gf();
      fault_state_           = make_gf();
      fault_normal_stress_   = make_gf();

      // Build fault face → (elem1, elem2) mapping for interior faces.
      // The driver list may include orphan / shared-as-bdr faces (kept to
      // preserve fi*nbf index alignment with fault_coords); for those the
      // interior-transformation lookup returns nullptr and we record -1 so
      // Set/UpdateFaultFieldsBP5 can skip the L2-p0 scatter.
      int n_int = fault_interior_faces.Size();
      fault_face_elem1_.resize(n_int);
      fault_face_elem2_.resize(n_int);
      for (int i = 0; i < n_int; i++)
      {
         int face = fault_interior_faces[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr)
         {
            fault_face_elem1_[i] = FTr->Elem1No;
            fault_face_elem2_[i] = FTr->Elem2No;
         }
         else
         {
            fault_face_elem1_[i] = -1;
            fault_face_elem2_[i] = -1;
         }
      }
      n_interior_fault_faces_ = n_int;

      // Shared faces: only Elem1 is local (Elem2 is on neighbor rank).
      // Map to local Elem1 only.
      n_shared_fault_faces_ = 0;
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
#ifdef MFEM_USE_MPI
         n_shared_fault_faces_ = fault_shared_faces.Size();
         fault_shared_elem1_.resize(n_shared_fault_faces_);
         for (int i = 0; i < n_shared_fault_faces_; i++)
         {
            int sf = fault_shared_faces[i];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            fault_shared_elem1_[i] = FTr->Elem1No;
         }
#endif
      }

      // Friction parameter / coordinate fields (static)
      fault_param_a_   = make_gf();
      fault_param_Dc_  = make_gf();
      fault_coord_x2_  = make_gf();
      fault_coord_x3_  = make_gf();

      // PLAN_split_bulk_solutions_2026-05-12: register the 12 L2-p0
      // projection fields with the kinematics collection ONLY when the
      // opt-in toggle `register_fault_projections_in_volume_pv_` is
      // set.  Default OFF — these duplicate the data in fault.vtkhdf
      // and used to be the largest contributor to volume.vtkhdf file
      // size.  The GFs are still allocated above so accessors and
      // direct reads keep working; only the registration with `pv_dc_`
      // (which controls what lands in the .vtkhdf) is gated.
      if (register_fault_projections_in_volume_pv_)
      {
         pv_dc_->RegisterField("slip_dip",         fault_slip_dip_.get());
         pv_dc_->RegisterField("slip_strike",      fault_slip_strike_.get());
         pv_dc_->RegisterField("slip_rate_dip",    fault_slip_rate_dip_.get());
         pv_dc_->RegisterField("slip_rate_strike", fault_slip_rate_strike_.get());
         pv_dc_->RegisterField("traction_dip",     fault_trac_dip_.get());
         pv_dc_->RegisterField("traction_strike",  fault_trac_strike_.get());
         pv_dc_->RegisterField("state_variable",   fault_state_.get());
         pv_dc_->RegisterField("normal_stress",    fault_normal_stress_.get());
         pv_dc_->RegisterField("param_a",          fault_param_a_.get());
         pv_dc_->RegisterField("param_Dc",         fault_param_Dc_.get());
         pv_dc_->RegisterField("fault_x2",         fault_coord_x2_.get());
         pv_dc_->RegisterField("fault_x3",         fault_coord_x3_.get());
      }

      fault_interior_faces_ = fault_interior_faces;
      fault_shared_faces_ = fault_shared_faces;

      has_fault_output_ = true;
   }

   /// @brief Set static friction parameters and coordinates on fault elements.
   ///
   /// Call once after InitFaultOutputBP5. Vectors are in local (all faces)
   /// layout with 1 component per DOF.
   void SetFaultParamsBP5(const Vector &local_a,
                          const Vector &local_Dc,
                          const Vector &local_x2,
                          const Vector &local_x3)
   {
      if (!has_fault_output_) { return; }

      *fault_param_a_  = 0.0;
      *fault_param_Dc_ = 0.0;
      *fault_coord_x2_ = 0.0;
      *fault_coord_x3_ = 0.0;

      const int nbf = nbf_per_face_;
      const int n_int = n_interior_fault_faces_;

      for (int fi = 0; fi < n_int; fi++)
      {
         int base = fi * nbf;
         real_t avg_a = 0, avg_Dc = 0, avg_x2 = 0, avg_x3 = 0;
         for (int k = 0; k < nbf; k++)
         {
            avg_a  += local_a(base + k);
            avg_Dc += local_Dc(base + k);
            avg_x2 += local_x2(base + k);
            avg_x3 += local_x3(base + k);
         }
         real_t inv = 1.0 / nbf;
         avg_a *= inv; avg_Dc *= inv; avg_x2 *= inv; avg_x3 *= inv;

         int e1 = fault_face_elem1_[fi];
         int e2 = fault_face_elem2_[fi];
         // Skip orphan / shared-as-bdr entries kept only for index alignment.
         if (e1 < 0) { continue; }
         (*fault_param_a_)(e1)  = avg_a;
         (*fault_param_Dc_)(e1) = avg_Dc;
         (*fault_coord_x2_)(e1) = avg_x2;
         (*fault_coord_x3_)(e1) = avg_x3;
         if (e2 >= 0)
         {
            (*fault_param_a_)(e2)  = avg_a;
            (*fault_param_Dc_)(e2) = avg_Dc;
            (*fault_coord_x2_)(e2) = avg_x2;
            (*fault_coord_x3_)(e2) = avg_x3;
         }
      }

      for (int si = 0; si < n_shared_fault_faces_; si++)
      {
         int fi = n_int + si;
         int base = fi * nbf;
         real_t avg_a = 0, avg_Dc = 0, avg_x2 = 0, avg_x3 = 0;
         for (int k = 0; k < nbf; k++)
         {
            avg_a  += local_a(base + k);
            avg_Dc += local_Dc(base + k);
            avg_x2 += local_x2(base + k);
            avg_x3 += local_x3(base + k);
         }
         real_t inv = 1.0 / nbf;
         avg_a *= inv; avg_Dc *= inv; avg_x2 *= inv; avg_x3 *= inv;

         int e1 = fault_shared_elem1_[si];
         (*fault_param_a_)(e1)  = avg_a;
         (*fault_param_Dc_)(e1) = avg_Dc;
         (*fault_coord_x2_)(e1) = avg_x2;
         (*fault_coord_x3_)(e1) = avg_x3;
      }
   }

   /// @brief Update all fault L2-p0 fields from owned-DOF vectors.
   ///
   /// Maps fault DOF values to the adjacent volume elements.  Each fault
   /// face has `nbf_per_face` DOFs; we average them for the element value.
   /// Vectors are "local" (all fault faces on this rank, interior+shared),
   /// with 2-component layout: [comp0_dof0, comp0_dof1, ..., comp1_dof0, ...]
   ///
   /// @param local_slip       Slip vector (2 * num_local_fault_dofs)
   /// @param local_slip_rate  Slip rate vector (2 * num_local_fault_dofs)
   /// @param local_traction   Traction vector (2 * num_local_fault_dofs)
   /// @param local_state      State (psi) vector (1 * num_local_fault_dofs)
   /// @param local_normal_stress  Normal stress vector (1 * num_local_fault_dofs),
   ///                             may be empty if elastic sigma_n is disabled
   void UpdateFaultFieldsBP5(const Vector &local_slip,
                             const Vector &local_slip_rate,
                             const Vector &local_traction,
                             const Vector &local_state,
                             const Vector &local_normal_stress = Vector())
   {
      if (!has_fault_output_) { return; }

      // Zero all fields (most elements won't have fault data)
      *fault_slip_dip_ = 0.0;        *fault_slip_strike_ = 0.0;
      *fault_slip_rate_dip_ = 0.0;   *fault_slip_rate_strike_ = 0.0;
      *fault_trac_dip_ = 0.0;        *fault_trac_strike_ = 0.0;
      *fault_state_ = 0.0;
      *fault_normal_stress_ = 0.0;
      const bool has_normal = (local_normal_stress.Size() > 0);

      const int nbf = nbf_per_face_;
      const int n_int = n_interior_fault_faces_;

      // Interior faces: assign to both elem1 and elem2
      for (int fi = 0; fi < n_int; fi++)
      {
         int base = fi * nbf;  // local fault DOF offset
         real_t s_d = 0, s_s = 0, sr_d = 0, sr_s = 0;
         real_t tr_d = 0, tr_s = 0, psi = 0, sn = 0;

         for (int k = 0; k < nbf; k++)
         {
            int dof = base + k;
            s_d  += local_slip(2 * dof);
            s_s  += local_slip(2 * dof + 1);
            sr_d += local_slip_rate(2 * dof);
            sr_s += local_slip_rate(2 * dof + 1);
            tr_d += local_traction(2 * dof);
            tr_s += local_traction(2 * dof + 1);
            psi  += local_state(dof);
            if (has_normal) { sn += local_normal_stress(dof); }
         }
         real_t inv = 1.0 / nbf;
         s_d *= inv; s_s *= inv; sr_d *= inv; sr_s *= inv;
         tr_d *= inv; tr_s *= inv; psi *= inv; sn *= inv;

         int e1 = fault_face_elem1_[fi];
         int e2 = fault_face_elem2_[fi];
         // Skip orphan / shared-as-bdr entries kept only for index alignment.
         if (e1 < 0) { continue; }
         (*fault_slip_dip_)(e1)         = s_d;
         (*fault_slip_strike_)(e1)      = s_s;
         (*fault_slip_rate_dip_)(e1)    = sr_d;
         (*fault_slip_rate_strike_)(e1) = sr_s;
         (*fault_trac_dip_)(e1)         = tr_d;
         (*fault_trac_strike_)(e1)      = tr_s;
         (*fault_state_)(e1)            = psi;
         if (has_normal) { (*fault_normal_stress_)(e1) = sn; }
         if (e2 >= 0)
         {
            (*fault_slip_dip_)(e2)         = s_d;
            (*fault_slip_strike_)(e2)      = s_s;
            (*fault_slip_rate_dip_)(e2)    = sr_d;
            (*fault_slip_rate_strike_)(e2) = sr_s;
            (*fault_trac_dip_)(e2)         = tr_d;
            (*fault_trac_strike_)(e2)      = tr_s;
            (*fault_state_)(e2)            = psi;
            if (has_normal) { (*fault_normal_stress_)(e2) = sn; }
         }
      }

      // Shared faces: only elem1 is local
      for (int si = 0; si < n_shared_fault_faces_; si++)
      {
         int fi = n_int + si;  // local fault face index (after interior)
         int base = fi * nbf;
         real_t s_d = 0, s_s = 0, sr_d = 0, sr_s = 0;
         real_t tr_d = 0, tr_s = 0, psi = 0, sn = 0;

         for (int k = 0; k < nbf; k++)
         {
            int dof = base + k;
            s_d  += local_slip(2 * dof);
            s_s  += local_slip(2 * dof + 1);
            sr_d += local_slip_rate(2 * dof);
            sr_s += local_slip_rate(2 * dof + 1);
            tr_d += local_traction(2 * dof);
            tr_s += local_traction(2 * dof + 1);
            psi  += local_state(dof);
            if (has_normal) { sn += local_normal_stress(dof); }
         }
         real_t inv = 1.0 / nbf;
         s_d *= inv; s_s *= inv; sr_d *= inv; sr_s *= inv;
         tr_d *= inv; tr_s *= inv; psi *= inv; sn *= inv;

         int e1 = fault_shared_elem1_[si];
         (*fault_slip_dip_)(e1)         = s_d;
         (*fault_slip_strike_)(e1)      = s_s;
         (*fault_slip_rate_dip_)(e1)    = sr_d;
         (*fault_slip_rate_strike_)(e1) = sr_s;
         (*fault_trac_dip_)(e1)         = tr_d;
         (*fault_trac_strike_)(e1)      = tr_s;
         (*fault_state_)(e1)            = psi;
         if (has_normal) { (*fault_normal_stress_)(e1) = sn; }
      }
   }

   // ---------------------------------------------------------------
   //  Fault surface VTU output (proper 2D face geometry)
   // ---------------------------------------------------------------

   /// @brief Write fault surface as a triangle VTU, one averaged value per cell.
   ///
   /// Each rank writes fault_surface_r{rank}_c{cycle}.vtu containing the
   /// fault-face triangles and one CellData value per triangle.  Each cell
   /// value is the arithmetic mean of the DG field evaluated at the face's
   /// `nbf_per_face_` quadrature/nodal points.  Averaging (rather than
   /// per-vertex interpolation from QP-at-vertex coordinates) avoids the
   /// R-001 / H-V91-A4 speckle artefact that arose when QP values were
   /// written to reference-triangle vertex positions directly (plan
   /// v9.1.0 §2.3).  It also side-steps the BR2 latent bug R-005 where
   /// the old hardcoded `k<3` loop crossed face boundaries for `nbf=1`.
   /// Rank 0 also writes a .pvtu index.
   ///
   /// @param prefix   Output directory path
   /// @param cycle    Time step number
   /// @param time     Simulation time
   /// @param rank     MPI rank
   /// @param nranks   Total number of MPI ranks
   /// @param local_slip_rate_k4  (optional, R-V92-E02 diagnostic) stage-4
   ///                            DOFData.V1/V2 snapshot, 2*num_fault_total.
   ///                            If `Size()==0` the `*_k4` fields are skipped.
   /// @param local_traction_k4   (optional, R-V92-E02 diagnostic) stage-4
   ///                            DOFData.tau1_corr/tau2_corr snapshot.
   /// @param local_normal_stress_k4  (optional, R-V92-E02 diagnostic) stage-4
   ///                                DOFData.sigma_n_corr snapshot.
   ///
   /// When the three `_k4` vectors are supplied (all non-empty), five
   /// extra CellData fields are emitted alongside the averaged ones:
   ///   slip_rate_dip_k4, slip_rate_strike_k4,
   ///   traction_dip_k4, traction_strike_k4,
   ///   normal_stress_k4.
   /// ParaView can then compute `sigma_n_k4 - normal_stress` to quantify
   /// the RK4-averaging discrepancy; a non-zero per-face delta confirms
   /// that the end-of-pipeline DOFData the station writer reports is a
   /// stage-averaged value distinct from the stage-4 friction-solve
   /// result (R-V92-E02 / plan §19 H-V92-K discriminator).
   void WriteFaultSurfaceVTU(
      const std::string &prefix,
      int cycle, real_t time, int rank, int nranks,
      const Vector &local_slip,
      const Vector &local_slip_rate,
      const Vector &local_traction,
      const Vector &local_state,
      const Vector &local_normal_stress,
      const Vector &local_a,
      const Vector &local_Dc,
      const Vector &local_x2,
      const Vector &local_x3,
      const Vector &local_slip_rate_k4 = Vector(),
      const Vector &local_traction_k4 = Vector(),
      const Vector &local_normal_stress_k4 = Vector())
   {
      // The early-return is uniform across ranks because `has_fault_output_`
      // is set in `InitFaultOutputBP5`, which all ranks call collectively.
      // This keeps the gather collective from deadlocking (Risk #1 in plan).
      if (!has_fault_output_) { return; }
      const int nbf = nbf_per_face_;
      const int n_int = n_interior_fault_faces_;
      const int n_shared = n_shared_fault_faces_;
      const bool has_normal = (local_normal_stress.Size() > 0);

      // R-007 (v9.1.0 rev 3): hard-fail on nbf<=0 in both Debug and Release.
      // The per-face average below divides by nbf; nbf=0 would produce
      // inv_nbf = 1/0 = +Inf and 0*Inf = NaN in every CellData entry.
      // InitFaultOutputBP5 assigns nbf_per_face_ = caller-supplied value
      // with no sanity check, so guard here at the call boundary.
      MFEM_VERIFY(nbf > 0,
                  "WriteFaultSurfaceVTU: nbf_per_face_ must be > 0, got "
                  << nbf);

#ifdef MFEM_USE_MPI
      // R-005: Defensive uniformity check.  The legacy branch below uses
      // no MPI; the binary branch enters a collective gather.  If a buggy
      // driver sets `legacy_ascii_` or `output_mode_` non-uniformly across
      // ranks, those branches diverge and the gather collective deadlocks.
      // Convert that silent hang into a clear abort here.  R-112: the
      // uniformity invariant is set up once at driver-init time and
      // never changes mid-run, so cache the verdict and skip the
      // collective on subsequent calls (an np=800 × 5000-cycle BP5 run
      // would otherwise pay 4M extra MPI_Bcast calls).
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
         if (!output_mode_uniformity_checked_)
         {
            int local_state[2] = { legacy_ascii_ ? 1 : 0,
                                   static_cast<int>(output_mode_) };
            int root_state[2]  = { local_state[0], local_state[1] };
            MPI_Bcast(root_state, 2, MPI_INT, 0, mesh_.GetComm());
            MFEM_VERIFY(local_state[0] == root_state[0] &&
                        local_state[1] == root_state[1],
                        "ParaViewOutput: legacy_ascii_ / output_mode_ differ "
                        "across ranks (rank " << rank << " has {"
                        << local_state[0] << "," << local_state[1]
                        << "}, rank 0 has {" << root_state[0] << ","
                        << root_state[1] << "}); call SetLegacyAsciiVTU / "
                        "SetFaultOutputMode uniformly on every rank.");
            output_mode_uniformity_checked_ = true;
         }
      }
#endif

      // R-002: Take the legacy ASCII back-compat path BEFORE any new MPI
      // collective so legacy mode remains a faithful preserve of the
      // pre-Phase-1 behaviour with no new MPI work.  The four pre-Phase-1
      // unit tests opt into this via SetLegacyAsciiVTU(true) so their
      // ASCII regex parsers keep working.
      if (legacy_ascii_)
      {
         WriteFaultSurfaceVTULegacyAscii(prefix, cycle, time, rank, nranks,
                                         local_slip, local_slip_rate,
                                         local_traction, local_state,
                                         local_normal_stress,
                                         local_a, local_Dc, local_x2, local_x3,
                                         local_slip_rate_k4,
                                         local_traction_k4,
                                         local_normal_stress_k4);
         return;
      }

      // -- HDF5 mode is dispatched AFTER the pack is built (it shares the
      //    per-face loop with the Vtu path).  On non-HDF5 builds the Hdf5
      //    enum value is unreachable because the constructor default
      //    falls back to Vtu.  Driver-side parse-time rejection of
      //    --paraview-fault-hdf5 on non-HDF5 builds (plan §Phase 2b
      //    "Files to Modify") guarantees external callers cannot reach
      //    Hdf5 without HDF5; the abort below is the last-resort guard
      //    for direct C++ misuse on non-HDF5 builds.
#ifndef MFEM_USE_HDF5
      if (output_mode_ == FaultOutputMode::Hdf5)
      {
         MFEM_ABORT("ParaViewOutput::WriteFaultSurfaceVTU: "
                    "FaultOutputMode::Hdf5 was selected but this build "
                    "of MFEM does not define MFEM_USE_HDF5.  Rebuild "
                    "with MFEM_USE_HDF5=YES, or call "
                    "SetFaultOutputMode(FaultOutputMode::Vtu).");
      }
#endif

      // R-V92-E02 diagnostic: stage-4 (non-averaged) DOFData snapshot.
      // Populated only if all three `_k4` vectors are non-empty on THIS
      // rank (`has_k4_local`).  Across ranks, `has_k4` is the OR-reduce
      // so every rank packs the same field-name list — required by the
      // gather's hash check (fault_vtu_binary.hpp:GatherFaultPackToRoot).
      const bool has_k4_local = (local_slip_rate_k4.Size() > 0
                                 && local_traction_k4.Size() > 0
                                 && local_normal_stress_k4.Size() > 0);
      bool has_k4 = has_k4_local;
#ifdef MFEM_USE_MPI
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
         int g = has_k4_local ? 1 : 0;
         MPI_Allreduce(MPI_IN_PLACE, &g, 1, MPI_INT, MPI_LOR,
                       mesh_.GetComm());
         has_k4 = (g != 0);
      }
#endif

      // -- Build the per-rank LocalFaultPack from the existing per-face
      //    average loop.  The pack carries the per-cell (one value per
      //    triangle) means; geometry uses the reference-triangle corners
      //    transformed to physical space (matches the legacy writer).
      // Per-cell field values
      std::vector<std::array<double,3>> vertices;
      std::vector<std::array<int,3>> triangles;
      std::vector<double> c_sd, c_ss, c_srd, c_srs, c_td, c_ts, c_psi, c_sn;
      std::vector<double> c_a, c_Dc, c_x2, c_x3;
      std::vector<double> c_srd_k4, c_srs_k4, c_td_k4, c_ts_k4, c_sn_k4;

      auto process_face = [&](int fi, int face_mesh_idx, bool is_shared)
      {
         FaceElementTransformations *FTr = nullptr;
         if (!is_shared)
         {
            int face = fault_interior_faces_[fi];
            FTr = mesh_.GetInteriorFaceTransformations(face);
         }
         else
         {
            if constexpr (std::is_same_v<MeshType, ParMesh>)
            {
#ifdef MFEM_USE_MPI
               int sf = fi;
               FTr = mesh_.GetSharedFaceTransformations(sf);
#endif
            }
         }
         if (!FTr) { return; }

         // Reference triangle vertices: (0,0), (1,0), (0,1).  Kept for
         // geometry; the per-cell data array below holds the single
         // averaged value for this triangle.
         int base_vert = static_cast<int>(vertices.size());
         const double ref_tri[3][2] = {{0,0}, {1,0}, {0,1}};
         for (int v = 0; v < 3; v++)
         {
            IntegrationPoint ip;
            ip.Set2(ref_tri[v][0], ref_tri[v][1]);
            Vector coords(3);
            FTr->Face->Transform(ip, coords);
            vertices.push_back({coords(0), coords(1), coords(2)});
         }
         triangles.push_back({base_vert, base_vert+1, base_vert+2});

         // Face-averaged cell-data values (over this face's nbf DOFs).
         // Loop bound `k<nbf` handles every DG method:
         //   BR2 (nbf=1): single centroid value, no mixing (R-005 fix).
         //   IP  (nbf=3): mean of three vertex values.
         //   TPV102 wave operator (nbf=3): mean of three interior QPs.
         //   Higher order (nbf=6, ...): mean of all DOFs (R-002 fix).
         // (nbf>0 enforced at function entry via MFEM_VERIFY — R-007.)
         const int base = face_mesh_idx * nbf;
         double a_sd = 0.0, a_ss = 0.0, a_srd = 0.0, a_srs = 0.0;
         double a_td = 0.0, a_ts = 0.0, a_psi = 0.0, a_sn = 0.0;
         double a_a  = 0.0, a_Dc = 0.0, a_x2  = 0.0, a_x3 = 0.0;
         double a_srd_k4 = 0.0, a_srs_k4 = 0.0;
         double a_td_k4  = 0.0, a_ts_k4  = 0.0, a_sn_k4 = 0.0;
         for (int k = 0; k < nbf; k++)
         {
            const int d = base + k;
            a_sd  += local_slip(2*d);
            a_ss  += local_slip(2*d+1);
            a_srd += local_slip_rate(2*d);
            a_srs += local_slip_rate(2*d+1);
            a_td  += local_traction(2*d);
            a_ts  += local_traction(2*d+1);
            a_psi += local_state(d);
            a_sn  += has_normal ? local_normal_stress(d) : 0.0;
            a_a   += local_a.Size()  > 0 ? local_a(d)  : 0.0;
            a_Dc  += local_Dc.Size() > 0 ? local_Dc(d) : 0.0;
            a_x2  += local_x2.Size() > 0 ? local_x2(d) : 0.0;
            a_x3  += local_x3.Size() > 0 ? local_x3(d) : 0.0;
            // R-006: gate the k4 reads on `has_k4_local`, NOT global
            // `has_k4`.  A rank with empty `_k4` Vectors but a global
            // has_k4=true (because rank 0 supplied them) would otherwise
            // OOB-read here.  When local is empty we leave the
            // accumulators at 0.0 and the post-loop push (gated by
            // global `has_k4`) records 0.0 for those k4 fields.
            if (has_k4_local)
            {
               a_srd_k4 += local_slip_rate_k4(2*d);
               a_srs_k4 += local_slip_rate_k4(2*d+1);
               a_td_k4  += local_traction_k4(2*d);
               a_ts_k4  += local_traction_k4(2*d+1);
               a_sn_k4  += local_normal_stress_k4(d);
            }
         }
         const double inv_nbf = 1.0 / static_cast<double>(nbf);
         c_sd.push_back(a_sd  * inv_nbf);
         c_ss.push_back(a_ss  * inv_nbf);
         c_srd.push_back(a_srd * inv_nbf);
         c_srs.push_back(a_srs * inv_nbf);
         c_td.push_back(a_td  * inv_nbf);
         c_ts.push_back(a_ts  * inv_nbf);
         c_psi.push_back(a_psi * inv_nbf);
         c_sn.push_back(a_sn  * inv_nbf);
         c_a.push_back(a_a   * inv_nbf);
         c_Dc.push_back(a_Dc  * inv_nbf);
         c_x2.push_back(a_x2  * inv_nbf);
         c_x3.push_back(a_x3  * inv_nbf);
         if (has_k4)
         {
            // Push uses global has_k4 so the field shape is identical
            // across ranks (the gather's hash check requires this).
            // Ranks where has_k4_local was false push 0.0 (accumulators
            // were untouched by the loop).
            c_srd_k4.push_back(a_srd_k4 * inv_nbf);
            c_srs_k4.push_back(a_srs_k4 * inv_nbf);
            c_td_k4.push_back (a_td_k4  * inv_nbf);
            c_ts_k4.push_back (a_ts_k4  * inv_nbf);
            c_sn_k4.push_back (a_sn_k4  * inv_nbf);
         }
      };

      // Interior faces
      for (int i = 0; i < n_int; i++)
      {
         process_face(i, i, false);
      }
      // Shared faces
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
#ifdef MFEM_USE_MPI
         for (int i = 0; i < n_shared; i++)
         {
            int sf_idx = fault_shared_faces_[i];
            process_face(sf_idx, n_int + i, true);
         }
#endif
      }

      // -- Field filter: same allow-list semantics as legacy.  Empty
      //    filter ⇒ emit all 12 (or 17 with _k4) fields.  Pruning at
      //    pack-build time means filtered arrays are not gathered, sent,
      //    or written.
      const bool filter_active = !fault_vtu_fields_.empty();
      auto want_field = [&](const char *name) -> bool
      {
         if (!filter_active) { return true; }
         return fault_vtu_fields_.count(std::string(name)) > 0;
      };

      // -- Assemble the LocalFaultPack.  Field order must be IDENTICAL
      //    on every rank (the gather hashes the concatenated field-name
      //    list).  The filter is uniform across ranks (it lives on a
      //    member of `*this`, set via SetFaultVTUFields), so order +
      //    membership are uniform.
      vtu::LocalFaultPack pack;
      pack.vertices  = std::move(vertices);
      pack.triangles = std::move(triangles);
      auto add_field = [&](const char *name, std::vector<double> &vals)
      {
         if (!want_field(name)) { return; }
         pack.field_names.emplace_back(name);
         pack.field_arrays.emplace_back(std::move(vals));
      };
      add_field("slip_dip",         c_sd);
      add_field("slip_strike",      c_ss);
      add_field("slip_rate_dip",    c_srd);
      add_field("slip_rate_strike", c_srs);
      add_field("traction_dip",     c_td);
      add_field("traction_strike",  c_ts);
      add_field("state_variable",   c_psi);
      add_field("normal_stress",    c_sn);
      add_field("param_a",          c_a);
      add_field("param_Dc",         c_Dc);
      add_field("fault_x2",         c_x2);
      add_field("fault_x3",         c_x3);
      if (has_k4)
      {
         add_field("slip_rate_dip_k4",    c_srd_k4);
         add_field("slip_rate_strike_k4", c_srs_k4);
         add_field("traction_dip_k4",     c_td_k4);
         add_field("traction_strike_k4",  c_ts_k4);
         add_field("normal_stress_k4",    c_sn_k4);
      }

      // -- Compute the comm once; both paths need it.
#ifdef MFEM_USE_MPI
      MPI_Comm comm = MPI_COMM_NULL;
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
         comm = mesh_.GetComm();
      }
#endif

      // -- Phase 2b HDF5 dispatch.  WriteFaultPackHdf does its own
      //    collective gather to rank 0 (reuses Phase 1's gather helper)
      //    and emits a single `<prefix>/fault.vtkhdf` (renamed from
      //    fault_surface.vtkhdf in PLAN_split_bulk_solutions_2026-05-12)
      //    rather than per-cycle VTUs.  No PVD bookkeeping in this mode.
      if (output_mode_ == FaultOutputMode::Hdf5)
      {
#ifdef MFEM_USE_HDF5
         // Filter selector (Deflate vs ZfpAccuracy) is configured on
         // `fault_hdf_state_` via `SetFaultHDFCompression` BEFORE the
         // first call (R-107).
         vtkhdf::WriteFaultPackHdf(fault_hdf_state_, prefix, cycle, time,
                                   pack, rank, nranks
#ifdef MFEM_USE_MPI
                                   , comm
#endif
                                  );
         return;
#else
         // Unreachable — guarded at function entry above on non-HDF5 builds.
         MFEM_ABORT("FaultOutputMode::Hdf5 requires MFEM_USE_HDF5=YES.");
#endif
      }

      // -- Vtu path (Phase 1): gather to rank 0 then emit one binary VTU.
      vtu::GatheredFaultPack gathered;
#ifdef MFEM_USE_MPI
      gathered = vtu::GatherFaultPackToRoot(pack, rank, nranks, comm);
#else
      (void)nranks;
      gathered = vtu::GatherFaultPackToRoot(pack, rank, nranks);
#endif

      // -- Rank 0: mkdir, write single binary VTU, append PVD entry.
      const std::string fault_dir = prefix + "/FaultSurface";
      if (rank == 0)
      {
         ::mkdir(fault_dir.c_str(), 0755);  // ignore EEXIST
         const std::string vtu_rel = "fault_surface_c"
                                   + std::to_string(cycle) + ".vtu";
         const std::string vtu_path = fault_dir + "/" + vtu_rel;
         vtu::WriteFaultPackVTU(vtu_path, gathered, VTKFormat::BINARY,
                                /*compression_level=*/0);
         fault_pvd_entries_.push_back({time, vtu_rel});
         WriteFaultPVD(fault_dir);
      }
      // R-003: No trailing barrier — only rank 0 does I/O in the binary
      // path (gather already synchronised) and `fault_dir` is not touched
      // by other ranks afterward.  Saves one collective per cycle.
   }

   // Stored shared face indices for surface VTU output
   Array<int> fault_shared_faces_;
   Array<int> fault_interior_faces_;

   // ---------------------------------------------------------------
   //  Save methods
   // ---------------------------------------------------------------

   /// Save based on adaptive V_max schedule.  All ranks must call
   /// collectively (ParaViewDataCollection::Save is MPI-collective).
   /// @param V_max  Global maximum slip rate (must be globally reduced
   ///               BEFORE calling so all ranks see the same value).
   bool Save(int cycle, real_t time, real_t V_max)
   {
      if (output_every_n_steps > 0)
      {
         // Step-based: write only at multiples of the interval
         if (cycle % output_every_n_steps == 0)
         {
            last_v_max_ = V_max;
            return ForceSaveImpl(cycle, time);
         }
         return false;
      }
      // Time-based: fixed dt or adaptive V_max schedule (with hysteresis).
      const int new_regime = adaptive_.NextRegime(V_max, current_regime_);
      real_t dt_out = (fixed_dt > 0.0)
                      ? fixed_dt
                      : SnapshotCapAwareInterval(V_max, new_regime, time);
      if (time - last_write_time_ < dt_out * kOutputTimeTolerance)
      {
         return false;
      }
      current_regime_ = new_regime;
      last_v_max_     = V_max;
      return ForceSaveImpl(cycle, time);
   }

   /// Schedule check without writing the volume PVD.  Returns true on the
   /// cycles/times when Save() would write, and advances last_write_time_
   /// (and current_regime_, last_v_max_) so subsequent scheduling stays
   /// consistent.  Use this when only the fault-surface VTU is wanted and
   /// the volume mesh+fields are suppressed to save disk space.
   bool ShouldWrite(int cycle, real_t time, real_t V_max)
   {
      if (output_every_n_steps > 0)
      {
         if (cycle % output_every_n_steps == 0)
         {
            last_write_time_ = time;
            last_v_max_      = V_max;
            // R-106: only bump if this cycle hasn't already been
            // committed by a paired Save() / ForceSave() / CommitSchedule
            // earlier in the same step.
            CommitOnceAtCycle(cycle);
            return true;
         }
         return false;
      }
      const int new_regime = adaptive_.NextRegime(V_max, current_regime_);
      real_t dt_out = (fixed_dt > 0.0)
                      ? fixed_dt
                      : SnapshotCapAwareInterval(V_max, new_regime, time);
      if (time - last_write_time_ < dt_out * kOutputTimeTolerance)
      {
         return false;
      }
      current_regime_  = new_regime;
      last_write_time_ = time;
      last_v_max_      = V_max;
      CommitOnceAtCycle(cycle);
      return true;
   }

   /// Read-only schedule check: returns true iff Save() (or ShouldWrite())
   /// would write at (cycle, time, V_max), without mutating last_write_time_,
   /// current_regime_, or last_v_max_.  Use this to gate expensive per-step
   /// packing work: then call the two-arg CommitSchedule(time, V_max) to
   /// advance state after the writes are committed.  The single-arg shim
   /// CommitSchedule(time) uses the last V_max recorded by Save/ShouldWrite
   /// (NOT by Peek — Peek is strictly read-only) and is correct for any
   /// call site that uses the default hysteresis_factor = 1.0 (where the
   /// regime is a stateless function of V).  Fault-only callers that set
   /// hysteresis_factor > 1 MUST use the two-arg CommitSchedule.
   bool PeekShouldWrite(int cycle, real_t time, real_t V_max) const
   {
      if (output_every_n_steps > 0)
      {
         return (cycle % output_every_n_steps == 0);
      }
      // Compute the prospective regime locally WITHOUT mutating state.
      const int prospective_regime = adaptive_.NextRegime(V_max, current_regime_);
      real_t dt_out = (fixed_dt > 0.0)
                      ? fixed_dt
                      : SnapshotCapAwareInterval(V_max, prospective_regime,
                                                 time);
      return (time - last_write_time_ >= dt_out * kOutputTimeTolerance);
   }

   /// Advance `last_write_time_` AND `current_regime_` AND `last_v_max_`
   /// in one atomic step after PeekShouldWrite has confirmed a write.
   /// Use this two-arg form from callers that skip Save (e.g. the
   /// fault-only path that only writes the fault-surface VTU) so the
   /// regime state machine keeps advancing — otherwise
   /// `current_regime_` would be pinned at 0 and hysteresis would
   /// silently no-op.
   void CommitSchedule(real_t time, real_t V_max)
   {
      // R-105: replace the original float-equality dedup with a
      // monotonic-step proxy.  CommitSchedule has no `cycle` argument,
      // so we synthesise one as `last_committed_cycle_ + 1`.  When a
      // paired Save / ShouldWrite / ForceSaveImpl has already committed
      // at the current step, `last_committed_cycle_` was just bumped,
      // and the increment-then-compare in CommitOnceAtCycle dedups
      // (because the synthesised cycle == last_committed_cycle_ + 1
      // exists, but only if the prior commit happened in this same
      // logical step — which we recognise by an FP-tolerant check on
      // `last_write_time_`).  Concretely: if `time` is within a tight
      // relative tolerance of `last_write_time_` (the most recent
      // recorded write time), treat this as the same step and dedup.
      // The tolerance is 8 ULPs in absolute units of `time`, which is
      // far below any meaningful simulation dt.
      const real_t tol =
         8.0 * std::numeric_limits<real_t>::epsilon()
             * std::max(std::abs(time), std::abs(last_write_time_));
      const bool same_step_as_last_commit =
         (last_committed_cycle_ != std::numeric_limits<int>::min())
         && (std::abs(time - last_write_time_) <= tol);
      current_regime_  = adaptive_.NextRegime(V_max, current_regime_);
      last_write_time_ = time;
      last_v_max_      = V_max;
      if (!same_step_as_last_commit)
      {
         CommitOnceAtCycle(last_committed_cycle_ + 1);
      }
   }

   /// Back-compat single-arg shim preserved for existing call sites
   /// (TPV102 driver and legacy BP5 paths).  Delegates to the two-arg
   /// form using the last V_max recorded by Save/ShouldWrite; this is
   /// correct for callers that either (a) use the default
   /// hysteresis_factor = 1.0 (regime is stateless in V), or (b) call
   /// Save/ShouldWrite just before this.  Callers that want hysteresis
   /// to advance through a PeekShouldWrite-gated path MUST call the
   /// two-arg overload instead.
   void CommitSchedule(real_t time) { CommitSchedule(time, last_v_max_); }

   /// Force a save at the current state.
   void ForceSave(int cycle, real_t time)
   {
      ForceSaveImpl(cycle, time);
   }

   /// Adaptive output interval based on maximum slip rate (legacy API,
   /// still used by test_io.cpp and any caller that wants the default
   /// thresholds/intervals without constructing a ParaViewOutput).  For
   /// runtime-configurable thresholds, mutate GetSchedule() instead.
   static real_t OutputInterval(real_t V_max)
   {
      AdaptiveSchedule s;  // defaults = legacy behavior
      return s.Interval(V_max, s.NextRegime(V_max, 0));
   }

   void SetDataFormat(VTKFormat fmt) { pv_dc_->SetDataFormat(fmt); }
   void SetHighOrderOutput(bool enable) { pv_dc_->SetHighOrderOutput(enable); }
   void SetLevelsOfDetail(int lod) { pv_dc_->SetLevelsOfDetail(lod); }
   bool HasFaultOutput() const { return has_fault_output_; }

private:
   MeshType &mesh_;
   int order_;
   real_t last_write_time_;
   /// Phase 6.1: abstract handle replaces the concrete
   /// `ParaViewDataCollection` member.  Either a `ParaViewDataCollection`
   /// (Vtu mode) or a `ParaViewHDFDataCollection` (Hdf5 mode), chosen
   /// at construction.  Virtual dispatch through the base for
   /// `RegisterField`, `SetCycle`, `SetTime`, `SetPrefixPath`,
   /// `SetDataFormat`, `SetHighOrderOutput`, `SetLevelsOfDetail`,
   /// `Save` keeps the call sites mode-agnostic.
   std::unique_ptr<ParaViewDataCollectionBase> pv_dc_;
   VolumeOutputMode volume_output_mode_;

#ifdef MFEM_USE_HDF5
   /// Phase 6.1 helper — return the active HDF collection if one is
   /// in use, else nullptr.  Used by `SetVolumeHDFCompression`
   /// (Phase 6.2) to dispatch HDF-only setters.
   ParaViewHDFDataCollection *GetHDFDC()
   {
      return dynamic_cast<ParaViewHDFDataCollection*>(pv_dc_.get());
   }
#endif

   // Adaptive-schedule state.  adaptive_ holds the tunable thresholds
   // and intervals (GetSchedule() exposes it for CLI overrides);
   // current_regime_ threads hysteresis across writes; last_v_max_
   // backs the single-arg CommitSchedule(time) shim.
   AdaptiveSchedule adaptive_;
   int    current_regime_ = 0;   // 0=interseismic, 1=nucleation, 2=coseismic
   real_t last_v_max_     = 0.0;

   // Phase 3 snapshot-cap state.  `total_run_time_` is set by drivers via
   // SetTotalRunTime; `total_snapshots_written_` is bumped via the
   // `CommitOnceAtCycle` helper on every committed save (Save /
   // ShouldWrite / CommitSchedule / ForceSaveImpl).
   real_t total_run_time_           = 0.0;
   int    total_snapshots_written_  = 0;

   // R-105 / R-106 dedup key.  All four committing paths gate the
   // counter bump on `cycle != last_committed_cycle_` via
   // CommitOnceAtCycle; this replaces the earlier float-equality on
   // (time, V_max) (fragile to 1-ULP FP noise) with an integer key
   // that is unambiguously the same step or not.  CommitSchedule
   // (which lacks a cycle argument) uses `commit_schedule_step_` as a
   // proxy: each PeekShouldWrite-true → CommitSchedule call sequence
   // is treated as advancing the step, but a Save/ShouldWrite that
   // ran at the same wallclock cycle prevents the counter bump.
   int last_committed_cycle_ = std::numeric_limits<int>::min();

   // Phase 4 volume-decouple state.  `emit_volume_save_` toggles whether
   // pv_.Save() fires inside ForceSaveImpl; `volume_pv_dt_` is an
   // optional cadence override (independent of the fault schedule);
   // `last_volume_write_time_` tracks when the volume save last fired
   // for the dt-cadence check.
   bool   emit_volume_save_       = true;
   real_t volume_pv_dt_           = 0.0;
   real_t last_volume_write_time_ = -1e30;
#ifdef MFEM_USE_HDF5
   // R-101 / Phase 2d.3: latched after the first SetVolumeHDFCompression
   // call that requests a non-default filter, so the explanatory rank-0
   // warning fires exactly once per run.
   mutable bool volume_hdf_compression_warned_ = false;
#endif

   // Fault L2-p0 output
   bool has_fault_output_ = false;
   int  nbf_per_face_ = 1;
   int  n_interior_fault_faces_ = 0;
   int  n_shared_fault_faces_ = 0;

   // Fault-surface VTU field filter.  Empty ⇒ emit all fields
   // (backward-compatible default).  Non-empty ⇒ emit only the
   // named CellData arrays.  See SetFaultVTUFields.
   std::set<std::string> fault_vtu_fields_;

   // PLAN_split_bulk_solutions_2026-05-12: opt-in for registering the
   // 12 L2-p0 fault projections with the primary (kinematics) volume
   // PV collection.  Default OFF.  See
   // SetRegisterFaultProjectionsInVolumePV.
   bool register_fault_projections_in_volume_pv_ = false;

   std::unique_ptr<FiniteElementCollection> fault_fec_;
   std::unique_ptr<FES> fault_fes_;

   std::unique_ptr<GF> fault_slip_dip_;
   std::unique_ptr<GF> fault_slip_strike_;
   std::unique_ptr<GF> fault_slip_rate_dip_;
   std::unique_ptr<GF> fault_slip_rate_strike_;
   std::unique_ptr<GF> fault_trac_dip_;
   std::unique_ptr<GF> fault_trac_strike_;
   std::unique_ptr<GF> fault_state_;
   std::unique_ptr<GF> fault_normal_stress_;

   // Friction parameter and coordinate fields (static, written once)
   std::unique_ptr<GF> fault_param_a_;
   std::unique_ptr<GF> fault_param_Dc_;
   std::unique_ptr<GF> fault_coord_x2_;
   std::unique_ptr<GF> fault_coord_x3_;

   // Interior fault face → element mapping
   std::vector<int> fault_face_elem1_;
   std::vector<int> fault_face_elem2_;
   // Shared fault face → local element mapping (Elem1 only)
   std::vector<int> fault_shared_elem1_;

   // Fault-surface output mode + PVD time-series entries (rank 0 only).
   // The PVD entry's `file` field points at either a per-cycle VTU
   // (Phase 1 binary path: `fault_surface_c{c}.vtu`) or a per-cycle
   // PVTU index (legacy ASCII path: `fault_surface_c{c}.pvtu`); the
   // PVD writer is agnostic to which.  The Hdf5 path bypasses the PVD
   // entirely (single .vtkhdf file at the prefix root).
   //
   // Default per plan §Phase 2b.Interfaces: tracks MFEM_USE_HDF5 so an
   // HDF5-enabled build defaults to the production Phase 2 path.
   FaultOutputMode  output_mode_  =
#ifdef MFEM_USE_HDF5
       FaultOutputMode::Hdf5;
#else
       FaultOutputMode::Vtu;
#endif
   bool             legacy_ascii_ = false;
   // R-112: latched on the first WriteFaultSurfaceVTU call after the
   // uniformity check passes; the check only depends on driver-init-time
   // state (output_mode_, legacy_ascii_) which is invariant mid-run.
   mutable bool     output_mode_uniformity_checked_ = false;
#ifdef MFEM_USE_HDF5
   // Persistent state for the rank-0 single-file VTKHDF writer.  Built
   // lazily on the first call to WriteFaultSurfaceVTU when output_mode_
   // == Hdf5.  Reused across all subsequent saves; the fault geometry
   // is fixed, only the per-cell field values change between cycles.
   mutable vtkhdf::FaultHDFState fault_hdf_state_;
#endif
   // Filename (relative to the FaultSurface directory) referenced by one
   // PVD `<DataSet>` entry.  Phase 1 binary path stores a `.vtu` filename
   // here; legacy ASCII path stores a `.pvtu` filename.  WriteFaultPVD
   // is agnostic to which.
   struct PVDEntry { real_t time; std::string file; };
   std::vector<PVDEntry> fault_pvd_entries_;

   /// R-105 / R-106: idempotent commit at @a cycle.  Increments
   /// `total_snapshots_written_` only when @a cycle differs from the
   /// last committed cycle.  This collapses a cycle-N → Save → ShouldWrite
   /// → CommitSchedule pattern (which previously triple-bumped the
   /// counter under R-106) into a single increment, AND removes the
   /// brittle float-equality dedup in the prior CommitSchedule (R-105).
   /// Returns true iff the counter was incremented.
   bool CommitOnceAtCycle(int cycle)
   {
      if (cycle == last_committed_cycle_) { return false; }
      last_committed_cycle_ = cycle;
      ++total_snapshots_written_;
      return true;
   }

   bool ForceSaveImpl(int cycle, real_t time)
   {
      // Phase 4: gate the volume PV save on the user toggles.  Three
      // cases:
      //   (1) emit_volume_save_=false AND volume_pv_dt_<=0:
      //       skip pv_.Save() entirely (production BP5 fault-only mode).
      //   (2) emit_volume_save_=true  AND volume_pv_dt_>0:
      //       fire only when (time - last_volume_write_time_) crosses
      //       the user dt, decoupling volume cadence from the fault
      //       schedule.
      //   (3) emit_volume_save_=true  AND volume_pv_dt_<=0 (default):
      //       fire on every fault commit (legacy behaviour).
      // The fault snapshot counter (`total_snapshots_written_`) bumps
      // unconditionally — the cap is on FAULT writes, not volume saves.
      if (emit_volume_save_)
      {
         bool volume_should_write = true;
         if (volume_pv_dt_ > 0.0)
         {
            volume_should_write =
               (time - last_volume_write_time_)
               >= volume_pv_dt_ * kOutputTimeTolerance;
         }
         if (volume_should_write)
         {
            // R-706: defensive — `pv_dc_` is always allocated (HDF or
            // VTU) by the constructor.  Catching a null dispatch here
            // converts a silent corruption into a loud abort if any
            // future refactor accidentally clears the unique_ptr.
            MFEM_VERIFY(pv_dc_,
                        "ParaViewOutput::ForceSaveImpl: pv_dc_ is null "
                        "but emit_volume_save_ is true.  This indicates "
                        "a constructor-time invariant violation.");
            // Phase 6.1: virtual dispatch through `pv_dc_` so
            // VolumeOutputMode::Vtu (legacy per-rank VTU) and
            // VolumeOutputMode::Hdf5 (single-file VTKHDF) share
            // the same call site.
            pv_dc_->SetCycle(cycle);
            pv_dc_->SetTime(time);
            pv_dc_->Save();
            last_volume_write_time_ = time;
         }
      }
      last_write_time_ = time;
      CommitOnceAtCycle(cycle);
      return true;
   }

   /// Phase 3: cap-aware variant of `adaptive_.Interval(V_max, regime)`.
   /// In the interseismic regime AND when `max_total_snapshots > 0`,
   /// inflates the interval per `RecomputeIntervalForCap`.  Coseismic
   /// and nucleation regimes return their own (uncapped) cadence.
   real_t SnapshotCapAwareInterval(real_t V_max, int regime,
                                    real_t time) const
   {
      const real_t base = adaptive_.Interval(V_max, regime);
      if (regime != 0 || adaptive_.max_total_snapshots <= 0)
      {
         return base;
      }
      // R-001: the cap requires a known total run time to project the
      // remaining budget over.  Drivers must call SetTotalRunTime(tfinal)
      // BEFORE the time-stepping loop.  Without it the cap silently
      // becomes a no-op — convert that into a clear abort so the user
      // notices instead of getting uncapped behaviour.
      MFEM_VERIFY(total_run_time_ > 0.0,
                  "ParaViewOutput: max_total_snapshots="
                  << adaptive_.max_total_snapshots
                  << " is set but SetTotalRunTime(tfinal) was not called "
                  "(total_run_time_ = " << total_run_time_ << ").  The cap "
                  "cannot project remaining writes without a known end time.");
      const real_t time_to_end = total_run_time_ - time;
      return adaptive_.RecomputeIntervalForCap(time_to_end,
                                               total_snapshots_written_);
   }

   /// Legacy per-rank ASCII fault-surface VTU writer.  Active only when
   /// `legacy_ascii_ == true` (set via `SetLegacyAsciiVTU(true)`).
   /// Layout is the pre-Phase-1 scheme: every rank writes
   /// `fault_surface_r{rank}_c{cycle}.vtu` (ASCII), rank 0 writes
   /// `fault_surface_c{cycle}.pvtu` indexing all per-rank pieces, and
   /// rank 0 appends a PVD entry pointing at the PVTU.  Preserved
   /// verbatim from the pre-Phase-1 implementation so the four ASCII-
   /// regex-based unit tests (continuity / k4 / field_filter /
   /// adjacent_triangle_uniformity) keep working without parser
   /// changes.  Drop after Phase 4 stabilises.
   void WriteFaultSurfaceVTULegacyAscii(
      const std::string &prefix,
      int cycle, real_t time, int rank, int nranks,
      const Vector &local_slip,
      const Vector &local_slip_rate,
      const Vector &local_traction,
      const Vector &local_state,
      const Vector &local_normal_stress,
      const Vector &local_a,
      const Vector &local_Dc,
      const Vector &local_x2,
      const Vector &local_x3,
      const Vector &local_slip_rate_k4,
      const Vector &local_traction_k4,
      const Vector &local_normal_stress_k4)
   {
      const int nbf = nbf_per_face_;
      const int n_int = n_interior_fault_faces_;
      const int n_shared = n_shared_fault_faces_;
      const bool has_normal = (local_normal_stress.Size() > 0);
      const bool has_k4 = (local_slip_rate_k4.Size() > 0
                           && local_traction_k4.Size() > 0
                           && local_normal_stress_k4.Size() > 0);

      std::vector<std::array<double,3>> vertices;
      std::vector<std::array<int,3>> triangles;
      std::vector<double> c_sd, c_ss, c_srd, c_srs, c_td, c_ts, c_psi, c_sn;
      std::vector<double> c_a, c_Dc, c_x2, c_x3;
      std::vector<double> c_srd_k4, c_srs_k4, c_td_k4, c_ts_k4, c_sn_k4;

      auto process_face = [&](int fi, int face_mesh_idx, bool is_shared)
      {
         FaceElementTransformations *FTr = nullptr;
         if (!is_shared)
         {
            int face = fault_interior_faces_[fi];
            FTr = mesh_.GetInteriorFaceTransformations(face);
         }
         else
         {
            if constexpr (std::is_same_v<MeshType, ParMesh>)
            {
#ifdef MFEM_USE_MPI
               int sf = fi;
               FTr = mesh_.GetSharedFaceTransformations(sf);
#endif
            }
         }
         if (!FTr) { return; }

         int base_vert = static_cast<int>(vertices.size());
         const double ref_tri[3][2] = {{0,0}, {1,0}, {0,1}};
         for (int v = 0; v < 3; v++)
         {
            IntegrationPoint ip;
            ip.Set2(ref_tri[v][0], ref_tri[v][1]);
            Vector coords(3);
            FTr->Face->Transform(ip, coords);
            vertices.push_back({coords(0), coords(1), coords(2)});
         }
         triangles.push_back({base_vert, base_vert+1, base_vert+2});

         const int base = face_mesh_idx * nbf;
         double a_sd = 0.0, a_ss = 0.0, a_srd = 0.0, a_srs = 0.0;
         double a_td = 0.0, a_ts = 0.0, a_psi = 0.0, a_sn = 0.0;
         double a_a  = 0.0, a_Dc = 0.0, a_x2  = 0.0, a_x3 = 0.0;
         double a_srd_k4 = 0.0, a_srs_k4 = 0.0;
         double a_td_k4  = 0.0, a_ts_k4  = 0.0, a_sn_k4 = 0.0;
         for (int k = 0; k < nbf; k++)
         {
            const int d = base + k;
            a_sd  += local_slip(2*d);
            a_ss  += local_slip(2*d+1);
            a_srd += local_slip_rate(2*d);
            a_srs += local_slip_rate(2*d+1);
            a_td  += local_traction(2*d);
            a_ts  += local_traction(2*d+1);
            a_psi += local_state(d);
            a_sn  += has_normal ? local_normal_stress(d) : 0.0;
            a_a   += local_a.Size()  > 0 ? local_a(d)  : 0.0;
            a_Dc  += local_Dc.Size() > 0 ? local_Dc(d) : 0.0;
            a_x2  += local_x2.Size() > 0 ? local_x2(d) : 0.0;
            a_x3  += local_x3.Size() > 0 ? local_x3(d) : 0.0;
            if (has_k4)
            {
               a_srd_k4 += local_slip_rate_k4(2*d);
               a_srs_k4 += local_slip_rate_k4(2*d+1);
               a_td_k4  += local_traction_k4(2*d);
               a_ts_k4  += local_traction_k4(2*d+1);
               a_sn_k4  += local_normal_stress_k4(d);
            }
         }
         const double inv_nbf = 1.0 / static_cast<double>(nbf);
         c_sd.push_back(a_sd  * inv_nbf);
         c_ss.push_back(a_ss  * inv_nbf);
         c_srd.push_back(a_srd * inv_nbf);
         c_srs.push_back(a_srs * inv_nbf);
         c_td.push_back(a_td  * inv_nbf);
         c_ts.push_back(a_ts  * inv_nbf);
         c_psi.push_back(a_psi * inv_nbf);
         c_sn.push_back(a_sn  * inv_nbf);
         c_a.push_back(a_a   * inv_nbf);
         c_Dc.push_back(a_Dc  * inv_nbf);
         c_x2.push_back(a_x2  * inv_nbf);
         c_x3.push_back(a_x3  * inv_nbf);
         if (has_k4)
         {
            c_srd_k4.push_back(a_srd_k4 * inv_nbf);
            c_srs_k4.push_back(a_srs_k4 * inv_nbf);
            c_td_k4.push_back (a_td_k4  * inv_nbf);
            c_ts_k4.push_back (a_ts_k4  * inv_nbf);
            c_sn_k4.push_back (a_sn_k4  * inv_nbf);
         }
      };

      for (int i = 0; i < n_int; i++) { process_face(i, i, false); }
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
#ifdef MFEM_USE_MPI
         for (int i = 0; i < n_shared; i++)
         {
            int sf_idx = fault_shared_faces_[i];
            process_face(sf_idx, n_int + i, true);
         }
#endif
      }

      const std::string fault_dir = prefix + "/FaultSurface";
      if (rank == 0) { ::mkdir(fault_dir.c_str(), 0755); }
#ifdef MFEM_USE_MPI
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
         MPI_Barrier(mesh_.GetComm());
      }
#endif
      const std::string vtu_name = fault_dir + "/fault_surface_r"
                                 + std::to_string(rank) + "_c"
                                 + std::to_string(cycle) + ".vtu";
      std::ofstream vtu(vtu_name);
      vtu << std::setprecision(10);
      const int npts   = static_cast<int>(vertices.size());
      const int ncells = static_cast<int>(triangles.size());
      vtu << "<?xml version=\"1.0\"?>\n";
      vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
      vtu << "<UnstructuredGrid>\n";
      vtu << "<Piece NumberOfPoints=\"" << npts
          << "\" NumberOfCells=\"" << ncells << "\">\n";
      vtu << "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" "
             "format=\"ascii\">\n";
      for (auto &v : vertices)
      {
         vtu << v[0] << " " << v[1] << " " << v[2] << "\n";
      }
      vtu << "</DataArray></Points>\n";
      vtu << "<Cells>\n";
      vtu << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
      for (auto &t : triangles)
      {
         vtu << t[0] << " " << t[1] << " " << t[2] << "\n";
      }
      vtu << "</DataArray>\n";
      vtu << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
      for (int i = 0; i < ncells; i++) { vtu << (i+1)*3 << "\n"; }
      vtu << "</DataArray>\n";
      vtu << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
      for (int i = 0; i < ncells; i++) { vtu << 5 << "\n"; }
      vtu << "</DataArray>\n";
      vtu << "</Cells>\n";

      const bool filter_active = !fault_vtu_fields_.empty();
      auto want_field = [&](const char *name) -> bool
      {
         if (!filter_active) { return true; }
         return fault_vtu_fields_.count(std::string(name)) > 0;
      };
      vtu << "<CellData>\n";
      auto write_field = [&](const char *name, const std::vector<double> &vals)
      {
         if (!want_field(name)) { return; }
         vtu << "<DataArray type=\"Float64\" Name=\"" << name
             << "\" format=\"ascii\">\n";
         for (double v : vals) { vtu << v << "\n"; }
         vtu << "</DataArray>\n";
      };
      write_field("slip_dip",         c_sd);
      write_field("slip_strike",      c_ss);
      write_field("slip_rate_dip",    c_srd);
      write_field("slip_rate_strike", c_srs);
      write_field("traction_dip",     c_td);
      write_field("traction_strike",  c_ts);
      write_field("state_variable",   c_psi);
      write_field("normal_stress",    c_sn);
      write_field("param_a",          c_a);
      write_field("param_Dc",         c_Dc);
      write_field("fault_x2",         c_x2);
      write_field("fault_x3",         c_x3);
      if (has_k4)
      {
         write_field("slip_rate_dip_k4",    c_srd_k4);
         write_field("slip_rate_strike_k4", c_srs_k4);
         write_field("traction_dip_k4",     c_td_k4);
         write_field("traction_strike_k4",  c_ts_k4);
         write_field("normal_stress_k4",    c_sn_k4);
      }
      vtu << "</CellData>\n";
      vtu << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
      vtu.close();

      if (rank == 0)
      {
         const std::string pvtu_rel = "fault_surface_c"
                                    + std::to_string(cycle) + ".pvtu";
         const std::string pvtu_name = fault_dir + "/" + pvtu_rel;
         std::ofstream pvtu(pvtu_name);
         pvtu << "<?xml version=\"1.0\"?>\n";
         pvtu << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">\n";
         pvtu << "<PUnstructuredGrid GhostLevel=\"0\">\n";
         pvtu << "<PPoints><PDataArray type=\"Float64\" "
                 "NumberOfComponents=\"3\"/></PPoints>\n";
         pvtu << "<PCellData>\n";
         const char *fields[] = {
            "slip_dip","slip_strike","slip_rate_dip","slip_rate_strike",
            "traction_dip","traction_strike","state_variable","normal_stress",
            "param_a","param_Dc","fault_x2","fault_x3"
         };
         for (auto f : fields)
         {
            if (!want_field(f)) { continue; }
            pvtu << "<PDataArray type=\"Float64\" Name=\"" << f << "\"/>\n";
         }
         if (has_k4)
         {
            const char *fields_k4[] = {
               "slip_rate_dip_k4","slip_rate_strike_k4",
               "traction_dip_k4","traction_strike_k4",
               "normal_stress_k4"
            };
            for (auto f : fields_k4)
            {
               if (!want_field(f)) { continue; }
               pvtu << "<PDataArray type=\"Float64\" Name=\"" << f
                    << "\"/>\n";
            }
         }
         pvtu << "</PCellData>\n";
         for (int r = 0; r < nranks; r++)
         {
            pvtu << "<Piece Source=\"fault_surface_r" << r
                 << "_c" << cycle << ".vtu\"/>\n";
         }
         pvtu << "</PUnstructuredGrid>\n</VTKFile>\n";
         pvtu.close();

         fault_pvd_entries_.push_back({time, pvtu_rel});
         WriteFaultPVD(fault_dir);
      }
   }

   /// Write (or overwrite) the fault surface PVD file with all entries so far.
   void WriteFaultPVD(const std::string &fault_dir)
   {
      std::string pvd_name = fault_dir + "/fault_surface.pvd";
      std::ofstream pvd(pvd_name, std::ios::trunc);
      pvd << std::setprecision(17);
      pvd << "<?xml version=\"1.0\"?>\n";
      pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
      pvd << "<Collection>\n";
      for (const auto &e : fault_pvd_entries_)
      {
         pvd << "<DataSet timestep=\"" << e.time
             << "\" file=\"" << e.file << "\"/>\n";
      }
      pvd << "</Collection>\n</VTKFile>\n";
      pvd.close();
   }

   // Factory for FES: serial vs parallel
   static std::unique_ptr<FES> MakeFES(MeshType &mesh,
                                       FiniteElementCollection *fec)
   {
      if constexpr (std::is_same_v<MeshType, ParMesh>)
      {
#ifdef MFEM_USE_MPI
         return std::make_unique<ParFiniteElementSpace>(&mesh, fec);
#else
         MFEM_ABORT("ParMesh requires MFEM_USE_MPI");
         return nullptr;
#endif
      }
      else
      {
         return std::make_unique<FiniteElementSpace>(&mesh, fec);
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PARAVIEW_OUTPUT_HPP
