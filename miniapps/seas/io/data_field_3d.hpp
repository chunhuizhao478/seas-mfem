// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// data_field_3d.hpp — schema-v1 HDF5 sidecar reader and trilinear
// interpolator for the data-projection feature.
//
// See miniapps/seas/document/features_dev/data_projection_schema_v1.md
// for the file-format specification and
// data_projection_feature_plan_v2.md for the consumer architecture.
//
// Behaviour highlights (Phase 3 contract):
//   * Loads a single named field from a schema-v1 HDF5 sidecar at
//     construction; verifies schema_version == "data_projection_v1",
//     crs == "EPSG:32611", z_positive == "elevation", axes monotone,
//     dataset shape == (Nx, Ny, Nz), forbids NaN, requires every cell
//     in [min_value, max_value].
//   * Evaluate(x, y, z) does trilinear interpolation; ANY out-of-bbox
//     query is a hard MFEM_ABORT (interpolation-only contract per
//     feature plan v2 §C-1 / §C-2).
//   * ContainsBBox(...) is the upfront mesh-vs-data containment gate
//     used by FieldProjector::Project (Phase 4).

#ifndef MFEM_SEAS_DATA_FIELD_3D_HPP
#define MFEM_SEAS_DATA_FIELD_3D_HPP

#include "mfem.hpp"

#include <array>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// Out-of-bbox query policy for ``DataField3D::Evaluate``.
///   - Abort (default, schema-v1 interpolation-only contract): any query
///     outside the data bbox is a hard MFEM_ABORT.
///   - Clamp (opt-in): the query coordinates are clamped to the nearest
///     data-axis edge before interpolation, so a mesh that extends beyond
///     the sidecar data hull (e.g. a large far-field absorbing box) yields
///     the edge value instead of aborting.  This is the ASAGI nearest-edge
///     hold behaviour SeisSol uses, reproduced here for SAFS meshes whose
///     far field exceeds the CVM data coverage.  In-bbox queries are
///     bit-identical to Abort (clamp only affects out-of-hull points).
enum class OOBPolicy : int
{
   Abort = 0,
   Clamp = 1
};

/// In-sidecar interpolation scheme for ``DataField3D::Evaluate``.
///
///   - Trilinear   (default, schema-v1 contract): 8-voxel stencil.
///                 O(h^2) accuracy on smooth fields.  Bit-exact at
///                 voxel corners.
///   - CatmullRom  (opt-in via SetInterpMode): tensor-product 64-voxel
///                 cardinal-cubic Catmull-Rom.  C^1 continuous,
///                 O(h^4) on smooth fields, no overshoot at C^1 kinks.
///                 Auto-falls back to trilinear within one voxel of
///                 the sidecar bbox where the 4-neighbour stencil is
///                 not available.  Bit-exact at voxel corners (cubic
///                 Catmull-Rom is interpolating).
enum class InterpMode : int
{
   Trilinear  = 0,
   CatmullRom = 1
};

class DataField3D
{
public:
   /// Load a single named field from a schema-v1 HDF5 sidecar.
   ///
   /// Aborts with a precise MFEM_ABORT message on:
   ///   - file open / dataset / attribute failure
   ///   - schema_version != "data_projection_v1"
   ///   - crs            != "EPSG:32611"
   ///   - z_positive     != "elevation"
   ///   - any axis non-monotone
   ///   - dataset shape != (Nx, Ny, Nz)
   ///   - any NaN in dataset
   ///   - any cell outside [min_value, max_value]
   ///   - oob not in {OOBPolicy::Abort, OOBPolicy::Clamp}.
   DataField3D(const std::string& sidecar_path,
               const std::string& field_name,
               OOBPolicy oob = OOBPolicy::Abort);

#ifdef MFEM_USE_MPI
   /// MPI-3 shared-memory overload
   /// (PLAN_sidecar_mpi_shared_memory_2026-07-17.md §4).
   ///
   /// ``node_comm`` must be a NODE-LOCAL communicator (all member ranks
   /// share physical memory; build it with
   /// ``MPI_Comm_split_type(..., MPI_COMM_TYPE_SHARED, ...)``).  The
   /// grid payload is then held in ONE ``MPI_Win_allocate_shared``
   /// window per node: node-local rank 0 performs the HDF5 read and the
   /// NaN/range validation; every other member maps the same physical
   /// pages via ``MPI_Win_shared_query``.  The tiny per-rank members
   /// (axes, bbox, scalars) stay replicated on every rank.
   ///
   /// Passing ``MPI_COMM_NULL`` selects the classic per-rank owning
   /// path, byte-identical to the 3-argument constructor.
   ///
   /// COLLECTIVE over ``node_comm`` (window allocation + publish
   /// barrier): every member rank must construct the same field from
   /// the same sidecar in the same order.  A validation failure on the
   /// loading rank aborts the WHOLE job (MFEM_ABORT calls MPI_Abort on
   /// MPI builds), so peers cannot hang at the publish barrier.
   DataField3D(const std::string& sidecar_path,
               const std::string& field_name,
               OOBPolicy oob,
               MPI_Comm node_comm);
#endif

   /// Frees the shared-memory window when this instance owns one.  If
   /// MPI has already been finalized (this driver calls MPI_Finalize
   /// before main's scope unwinds — see spatial_dyn_driver.cpp's
   /// MPIContext note), the free is skipped: the OS reclaims the
   /// shared segment at process exit, and calling MPI_Win_free after
   /// MPI_Finalize would abort.  SPMD control flow makes the
   /// skip-vs-free decision consistent across node-local ranks.
   ~DataField3D();

   /// Non-copyable: the instance may own an MPI shared-memory window
   /// handle, which cannot be duplicated.  (No in-repo code copied
   /// DataField3D before this change; StressField3D constructs its six
   /// members in place.)
   DataField3D(const DataField3D&) = delete;
   DataField3D& operator=(const DataField3D&) = delete;

   /// Evaluate at (x, y, z) in canonical CRS (UTM 11 N, m) using the
   /// currently-selected ``InterpMode`` (default Trilinear).
   /// Aborts on out-of-bbox query.
   real_t Evaluate(real_t x, real_t y, real_t z) const;

   /// Switch the in-sidecar interpolation scheme used by ``Evaluate``.
   /// Default is ``InterpMode::Trilinear`` (matches the pre-existing
   /// schema-v1 behaviour exactly; safe to leave unchanged).
   void SetInterpMode(InterpMode mode) { interp_mode_ = mode; }
   InterpMode GetInterpMode() const    { return interp_mode_; }

   /// Inclusive bbox in canonical CRS:
   /// {xmin, xmax, ymin, ymax, zmin, zmax}.
   const std::array<real_t, 6>& BBox() const { return bbox_; }

   /// True iff the supplied bbox is fully contained (with optional
   /// epsilon).  Used as the FieldProjector::Project pre-flight gate.
   bool ContainsBBox(real_t xmin, real_t xmax,
                     real_t ymin, real_t ymax,
                     real_t zmin, real_t zmax,
                     real_t eps = 0.0) const;

   const std::string& FieldName() const { return field_name_; }
   const std::string& Units()     const { return units_; }
   const std::string& Crs()       const { return crs_; }
   real_t MinValue()              const { return min_value_; }
   real_t MaxValue()              const { return max_value_; }

   /// Axis sizes for diagnostics / tests.
   int NumX() const { return static_cast<int>(x_.size()); }
   int NumY() const { return static_cast<int>(y_.size()); }
   int NumZ() const { return static_cast<int>(z_.size()); }

   /// Test-only accessor for find_index_().
   int FindIndexX(real_t v) const { return find_index_(x_, v); }
   int FindIndexY(real_t v) const { return find_index_(y_, v); }
   int FindIndexZ(real_t v) const { return find_index_(z_, v); }

private:
   /// Shared load body for both constructors: metadata + axes on every
   /// rank, then the field payload into either the per-rank owned
   /// vector (node_comm_ == MPI_COMM_NULL / serial build) or the
   /// node-shared MPI-3 window (loading rank fills + validates, peers
   /// map).  Sets ``data_`` / ``data_len_`` in both modes.
   void load_(const std::string& sidecar_path,
              const std::string& field_name);

   /// Trilinear (8-voxel) evaluation; the original Phase-3 path,
   /// unchanged.
   real_t evaluate_trilinear_(real_t x, real_t y, real_t z) const;

   /// Catmull-Rom tensor-product cubic (64-voxel) evaluation.  Falls
   /// back to ``evaluate_trilinear_`` when the query lies within one
   /// voxel of the sidecar bbox on any axis (no -1 / +2 neighbour
   /// available).
   real_t evaluate_catmull_rom_(real_t x, real_t y, real_t z) const;

   /// Return i with x_[i] <= v <= x_[i+1] for v ∈ [x_[0], x_[N-1]].
   /// Aborts if v is out of range.  N must be >= 2.
   int find_index_(const std::vector<real_t>& axis, real_t v) const;

   /// Index helper: flat index into data_ for 3-D coords (i, j, k).
   /// Uses (Nx, Ny, Nz) row-major layout (HDF5 C order).
   inline int flat_index_(int i, int j, int k) const
   {
      return (i * static_cast<int>(y_.size()) + j) *
             static_cast<int>(z_.size()) + k;
   }

   std::vector<real_t> x_, y_, z_;

   /// Grid payload, (Nx, Ny, Nz) row-major flat.  ``data_`` is the ONE
   /// read pointer every evaluator indexes; it aims at either
   ///   - ``data_owned_.data()``     (per-rank owning mode, the
   ///     pre-shared-memory behaviour, and the only mode on serial
   ///     builds), or
   ///   - the node-shared MPI-3 window base (shared mode; one physical
   ///     copy per node, mapped read-only by every node-local rank).
   /// Keeping the name ``data_`` leaves the trilinear / Catmull-Rom
   /// read code byte-for-byte unchanged.
   std::vector<real_t> data_owned_;
   const real_t* data_ = nullptr;
   size_t data_len_ = 0;
#ifdef MFEM_USE_MPI
   /// Node-local communicator used ONLY during load_ (not owned, not
   /// freed here); MPI_COMM_NULL selects the per-rank owning mode.
   MPI_Comm node_comm_ = MPI_COMM_NULL;
   /// Shared-memory window owning the payload in shared mode;
   /// MPI_WIN_NULL in owning mode.  Freed (collectively) in ~DataField3D
   /// unless MPI has already been finalized.
   MPI_Win shared_win_ = MPI_WIN_NULL;
#endif
   std::array<real_t, 6> bbox_;
   std::string field_name_;
   std::string units_;
   std::string crs_;
   real_t min_value_;
   real_t max_value_;
   OOBPolicy oob_policy_;
   InterpMode interp_mode_ = InterpMode::Trilinear;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DATA_FIELD_3D_HPP
