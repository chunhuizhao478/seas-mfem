// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// stress_field_3d.hpp — six-component (symmetric Cauchy) schema-v1
// HDF5 sidecar reader, Phase 6 §1 of PLAN_onfaultstress.md.
//
// Wraps six DataField3D instances in the schema-v1 canonical order
// (sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz) and
// returns the full symmetric tensor as a 3x3 mfem::DenseMatrix.
//
// Sign convention: this reader is a strict pass-through.  The bulk
// path's single-source sign flip from H&Z continuum-mechanics
// (compression negative) to SEAS-internal (compression positive)
// happens once in the Python preprocessor (Phase 3
// `bulk_stress_tensor_field`) at the source; Phase 5 (sidecar
// writer) does the MPa->Pa unit conversion only; this C++ reader
// returns the on-disk values unchanged (R-501 / R-502 contract).
//
// See PLAN_onfaultstress.md Phase 6 §1 (lines 1587-1645) for the
// detailed contract.
//
// Loaded fields per the schema-v1 spec at
// miniapps/seas/document/features_dev/data_projection_schema_v1.md:
//   sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz
// (units Pa, compression POSITIVE SEAS convention).

#ifndef MFEM_SEAS_STRESS_FIELD_3D_HPP
#define MFEM_SEAS_STRESS_FIELD_3D_HPP

#include "mfem.hpp"

#include "data_field_3d.hpp"

#include <array>
#include <string>

namespace mfem
{
namespace seas
{

/// Six-component symmetric Cauchy stress field reader.
///
/// Constructs six DataField3D instances from the same schema-v1
/// sidecar, one per component, in the canonical order
/// {sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz}.
/// Each component reader independently enforces the schema-v1
/// invariants (DataField3D ctor checks).
///
/// `Evaluate(x, y, z)` returns the full 3x3 symmetric tensor in Pa
/// (compression positive, SEAS internal convention).
class StressField3D
{
public:
   /// Load all six tensor components from a schema-v1 sidecar.
   ///
   /// Aborts via the underlying DataField3D ctor on:
   ///   - sidecar open / dataset / attribute failure on any of the
   ///     six fields,
   ///   - schema_version != "data_projection_v1", crs != "EPSG:32611",
   ///     z_positive != "elevation", non-monotone axes,
   ///   - any NaN or any cell outside the per-field [min, max].
   ///
   /// Additionally, asserts that all six components share an
   /// identical grid (axis sizes and bbox) — the Phase 5 writer
   /// guarantees this by construction; the assertion catches a
   /// hand-modified sidecar.
   StressField3D(const std::string& sidecar_path,
                 OOBPolicy oob = OOBPolicy::Abort);

#ifdef MFEM_USE_MPI
   /// MPI-3 shared-memory overload
   /// (PLAN_sidecar_mpi_shared_memory_2026-07-17.md §5 Phase 2):
   /// forwards ``node_comm`` to all six DataField3D component ctors, so
   /// the six grids live in six node-shared windows (one physical copy
   /// per node each) instead of six copies per rank.
   /// ``MPI_COMM_NULL`` selects the classic per-rank path.  COLLECTIVE
   /// over ``node_comm`` (see the DataField3D overload's contract).
   /// The windows are freed when this instance is destroyed — for the
   /// transient stress projection that is the end of apply_csm_impl,
   /// exactly the early-free the plan calls for.
   StressField3D(const std::string& sidecar_path,
                 OOBPolicy oob,
                 MPI_Comm node_comm);
#endif

   /// Evaluate the symmetric Cauchy stress tensor at (x, y, z) in
   /// canonical CRS (UTM 11 N, m), returned as a 3x3 dense matrix
   /// in Pa (compression POSITIVE, SEAS internal convention).
   ///
   /// Aborts on out-of-bbox query (interpolation-only contract
   /// inherited from DataField3D).
   mfem::DenseMatrix Evaluate(real_t x, real_t y, real_t z) const;

   /// Const access to the underlying component readers.
   /// component_index is interpreted in schema-v1 canonical order:
   ///   0 -> sigma_xx, 1 -> sigma_yy, 2 -> sigma_zz,
   ///   3 -> sigma_xy, 4 -> sigma_yz, 5 -> sigma_xz.
   /// Aborts on out-of-range index.
   const DataField3D& Field(int component_index) const;

   /// Inclusive bbox in canonical CRS:
   /// {xmin, xmax, ymin, ymax, zmin, zmax}.
   ///
   /// By construction (single Phase 5 grid), every component shares
   /// the same bbox; the ctor asserts this and stores the shared
   /// value. The intersection-of-six (per plan §1 line 1644-1645)
   /// equals the shared bbox in the well-formed case.
   const std::array<real_t, 6>& BBox() const { return bbox_; }

   /// True iff the supplied bbox is fully contained in every
   /// component's bbox.  Since all six share the same grid, this
   /// reduces to a single ContainsBBox check; provided as a
   /// convenience for callers that want the StressField3D-level
   /// gate analogous to FieldProjector::Project's pre-flight.
   bool ContainsBBox(real_t xmin, real_t xmax,
                     real_t ymin, real_t ymax,
                     real_t zmin, real_t zmax,
                     real_t eps = 0.0) const;

   /// Propagate the InterpMode to all six underlying readers in
   /// lock-step.  Default is InterpMode::Trilinear.
   void SetInterpMode(InterpMode mode);
   InterpMode GetInterpMode() const { return sigma_xx_.GetInterpMode(); }

   /// Number of components (always 6).
   static constexpr int NumComponents() { return 6; }

   /// Canonical field name for component_index (0..5).
   /// Returns the schema-v1 string used by the sidecar.
   static const std::string& ComponentName(int component_index);

private:
   /// Shared ctor tail: pin bbox_ (intersection of the six component
   /// bboxes) and run AssertConsistentGrid_.  Called by both ctors.
   void init_after_load_();

   /// Verify that all six underlying readers agree on grid axes /
   /// bbox / value range.  Aborts via MFEM_ABORT on mismatch.
   void AssertConsistentGrid_() const;

   // Six component readers, in schema-v1 canonical order.
   DataField3D sigma_xx_;
   DataField3D sigma_yy_;
   DataField3D sigma_zz_;
   DataField3D sigma_xy_;
   DataField3D sigma_yz_;
   DataField3D sigma_xz_;

   // Shared bbox (verified identical across all six in the ctor).
   std::array<real_t, 6> bbox_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_STRESS_FIELD_3D_HPP
