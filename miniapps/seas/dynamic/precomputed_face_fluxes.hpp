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

#ifndef MFEM_SEAS_PRECOMPUTED_FACE_FLUXES_HPP
#define MFEM_SEAS_PRECOMPUTED_FACE_FLUXES_HPP

#include "mfem.hpp"

#include "godunov_flux.hpp"
#include "wave_state.hpp"
#include "../domain/boundary_config.hpp"

// Forward-declare FreeSurfaceBCMode (defined in wave_operator.hpp) to
// avoid a circular include: wave_operator.hpp needs PrecomputedFaceFluxes
// as a WaveOperator member (TPV102 Phase 2a §6.3), and this header used
// to include wave_operator.hpp for the FreeSurfaceBCMode enum.  The
// enum-class forward declaration with explicit underlying type is
// sufficient for the Init parameter (C++11 §7.2/3).

#include <cstdint>
#include <set>
#include <unordered_map>
#include <vector>

namespace mfem
{
namespace seas
{

// Forward declaration (full definition in wave_operator.hpp).
enum class FreeSurfaceBCMode : int;

/// @brief BC classification for the precomputed per-cell per-side flux.
///
/// Per plan §5.1.2 / R-009.  SeisSol FaceType coverage:
///   Regular             -> Interior
///   FreeSurface         -> FreeSurfaceGodunov (mandatory for SeisSol parity)
///   DynamicRupture      -> delegated to FaultFaceFlux (untouched)
///   Outflow             -> NOT implemented (MFEM's Absorbing is not
///                          semantically equivalent; see plan §5.2.2)
///   FreeSurfaceGravity  -> out of scope
///   Dirichlet           -> out of scope
///   Periodic            -> out of scope
///   Analytical          -> out of scope
enum class FacePrecomputedBC
{
   Interior            = 0,
   FreeSurfaceGodunov  = 1,
   FreeSurfaceGamma    = 2,
   Absorbing           = 3
};

/// @brief Topology-based precomputed per-cell per-side flux rotation.
///
/// Implements Phase 1 of the "SeisSol-inspired topology-based precomputed
/// face-rotation" plan (2026-04-23).  Computes, at Init time, a 9x9 GLOBAL
/// frame matrix `nApNm1` (and `nAmNm1` for interior faces) for every
/// (element, local_side) whose mesh face is non-fault.  At run time,
/// `AddInteriorFaceRhs` / `AddBoundaryFaceRhs` apply the precomputed
/// matrices to the time-integrated state at each QP and emit the per-DOF
/// rhs contribution `-= w * shape(i) * F_h`.
///
/// For interior faces, the precomputed matrices reproduce
/// `GodunovFlux::Interior(n_outward_for_this_cell, I_self, I_nbr, F_h)` to
/// ULP; for boundary faces they reproduce `GodunovFlux::FreeSurfaceGodunovTotal`
/// / `FreeSurfaceTotal` / `AbsorbingTotal` on the per-cell outward normal.
/// See plan §2.4 and §2.5 for the equivalence contract.
class PrecomputedFaceFluxes
{
public:
   /// Per (element, local_side) stored entry.  Interior faces have TWO
   /// entries (one per adjacent element); boundary faces have ONE.
   struct FaceEntry
   {
      int    element       = -1;
      int    local_side    = -1;
      int    face_idx      = -1;
      int    neighbor_elem = -1;   ///< -1 if boundary
      FacePrecomputedBC bc = FacePrecomputedBC::Interior;
      real_t surface_area  = 0.0;
      real_t cell_volume   = 0.0;
      real_t normal[3]     = {0.0, 0.0, 0.0};
      // R4-003: for INTERIOR entries tangent1/tangent2 are load-bearing --
      // they are the topology-based frame used to build nApNm1 / nAmNm1.
      // For BOUNDARY entries they are DIAGNOSTIC-ONLY -- the probe-based
      // construction captures F_h via GodunovFlux::BuildFrame(n)'s own
      // internal tangents, which differ from these topology tangents.
      // Do NOT apply BuildRotation(n, tangent1, tangent2) to a boundary
      // entry's nApNm1 -- that would double-rotate and corrupt the
      // result.  Rotation-invariance of isotropic elastic BCs around `n`
      // makes the discrepancy harmless if the matrix is used as-is.
      real_t tangent1[3]   = {0.0, 0.0, 0.0};
      real_t tangent2[3]   = {0.0, 0.0, 0.0};
      DenseMatrix nApNm1;          ///< 9x9 GLOBAL frame
      DenseMatrix nAmNm1;          ///< 9x9 GLOBAL frame; zero for boundary

      // === Arm 3d: topology-based state-sampling tables =================
      // Phase 3 STOP findings (Arm 2b) proved MFEM's `CalcOrtho`-based face
      // normals and MFEM's `Loc1`/`Loc2` transforms are NOT D4-covariant
      // across orbit-pair diagonal faces: stored face-vertex ordering is
      // based on global coordinates, not per-cell outward orientation.
      //
      // Arm 3d fixes this by replacing the runtime `Loc1.Transform(ip, ip1)`
      // / `fe->CalcShape(ip1, shape1)` path with a PRECOMPUTED topology-
      // based shape table.  `shape_self[q][i]` is φ_i evaluated at the
      // reference point obtained by mapping QP q through the FACE's
      // topology-coordinate system (barycentric in the face triangle) to
      // the element's reference triangle using topology-consistent vertex
      // correspondence (by COORDINATE matching, not MFEM's stored
      // face-vertex ordering).  Same construction for the neighbor side.
      //
      // Invariant: on orbit-pair FaceEntries (face_a, face_b) under
      // y-mirror, `shape_self_a[q][i] == shape_self_b[q][i]` to ULP for
      // every (q, i) because the topology-based map is coordinate-
      // invariant and the MFEM reference vertex ordering is also
      // canonical (TETRAHEDRON reference vertices are fixed by MFEM).
      int nqp = 0;                       ///< face QP count
      std::vector<real_t> w_qp;          ///< topology weight per QP
      std::vector<std::vector<real_t>> shape_self;  ///< [nqp][ndof]
      std::vector<std::vector<real_t>> shape_nbr;   ///< [nqp][ndof]; empty for boundary
      int dof_offset_self = -1;          ///< caller's DOF offset
      int dof_offset_nbr  = -1;          ///< neighbor's DOF offset; -1 for boundary
      int ndof_per_el     = 0;
   };

   /// Build per-(elem, side) entries for every non-fault non-shared face.
   ///
   /// @param[in] mesh                 Mesh; non-const because
   ///                                 GetFaceElementTransformations is non-const.
   /// @param[in] fes                  Finite element space (for DOF sizing).
   /// @param[in] bc                   Boundary configuration (natural_attrs,
   ///                                 absorbing_attrs, fault_attr).
   /// @param[in] flux                 GodunovFlux (provides A^+ / A^- and
   ///                                 the probe entry points).
   /// @param[in] fs_mode              FreeSurfaceBCMode to use when
   ///                                 classifying natural_attrs faces.
   /// @param[in] fault_face_set       Face indices that must be skipped
   ///                                 (routed through FaultFaceFlux).
   /// @param[in] face_bdr_attr        Precomputed table sized mesh.GetNumFaces().
   /// @param[in] shared_mesh_face_set Shared (ParMesh) face indices,
   ///                                 skipped in Phase 1.
   /// @param[in] build_arm3d_tables   When true, populate the per-entry
   ///                                 `shape_self` / `shape_nbr` / `w_qp`
   ///                                 topology-based state-sampling
   ///                                 tables used by
   ///                                 AddInteriorFaceRhsFull /
   ///                                 AddBoundaryFaceRhsFull.  Defaults
   ///                                 OFF.  REVIEW.md R-003 (v9.5.0):
   ///                                 gated because the Full path shows
   ///                                 3.9e-01 rel drift vs runtime on
   ///                                 plane-wave inputs (see
   ///                                 wave_operator.inl:2101 note) and
   ///                                 the phase3 STOP freeze policy
   ///                                 forbids landing new flux-layer
   ///                                 mechanism on production callers
   ///                                 until Arm 1 justifies it.  The
   ///                                 Arm 3d unit tests opt in by
   ///                                 passing `true` explicitly.
   void Init(Mesh &mesh,
             const FiniteElementSpace &fes,
             const BoundaryConfig &bc,
             const GodunovFlux &flux,
             FreeSurfaceBCMode fs_mode,
             const std::set<int> &fault_face_set,
             const std::vector<int> &face_bdr_attr,
             const std::set<int> &shared_mesh_face_set,
             bool build_arm3d_tables = false);

#ifdef MFEM_USE_MPI
   /// TPV102 Phase 2b (§7.2): extend the precomputed tables with one
   /// FaceEntry per shared non-fault face owned locally.  Must be called
   /// AFTER `Init(...)` (which sets `initialized_ = true`); the shared
   /// pass appends to `entries_` and to `face_elem_to_entry_` keyed by
   /// `(mesh_face_idx, Elem1No)`.  `Elem1No` is the local owning cell's
   /// id and is always `>= 0`; the remote neighbor's id is the
   /// face-neighbor pseudo-id that MFEM returns as negative, which this
   /// rank does not own and does not store (peer rank builds its own
   /// entry independently).
   ///
   /// @param[in] pmesh           Parallel mesh.
   /// @param[in] fes             Finite element space (for DOF sizing).
   /// @param[in] flux            GodunovFlux (provides A^+ / A^-).
   /// @param[in] fault_face_set  Mesh face indices routed through
   ///                            FaultFaceFlux (skipped here).
   void InitSharedFaces(ParMesh &pmesh,
                        const FiniteElementSpace &fes,
                        const GodunovFlux &flux,
                        const std::set<int> &fault_face_set);
#endif

   /// Arm 3d: topology-based self-contained interior face accumulator.
   /// Reads Q directly from the global state vector using the precomputed
   /// topology shape tables (shape_self / shape_nbr) and writes the
   /// caller's per-DOF contribution.  D4-covariant by construction.
   ///
   /// Caller passes the global Q array and ndof_total; the function
   /// handles QP loop, state assembly, matvec, and accumulation.
   void AddInteriorFaceRhsFull(int face_idx,
                               int caller_elem,
                               const real_t *Q_data,
                               int ndof_total,
                               Vector &rhs) const;

   /// Arm 3d: topology-based self-contained boundary face accumulator.
   /// Reads Q from the caller's cell DOFs using the topology shape table
   /// and accumulates `-w · shape · (nApNm1 · I_self)` into rhs.
   void AddBoundaryFaceRhsFull(int face_idx,
                               const real_t *Q_data,
                               const real_t *bulk_bg_scaled,
                               int ndof_total,
                               Vector &rhs) const;

   /// Accumulate per-DOF rhs contribution from an interior face at one QP.
   ///
   /// Resolves the (face_idx, caller_elem) entry, computes
   ///   F_h = nApNm1 * I_self + nAmNm1 * I_nbr
   /// and subtracts `w * shape1(i) * F_h[c]` from rhs at the caller's DOFs.
   void AddInteriorFaceRhs(int face_idx,
                           int caller_elem,
                           const real_t *I_self, const real_t *I_nbr,
                           real_t w,
                           const real_t *shape1,
                           int ndof, int dof_offset1,
                           int ndof_total,
                           Vector &rhs) const;

   /// Accumulate per-DOF rhs contribution from a boundary face at one QP.
   ///
   /// Phase 1 precondition: `bulk_bg_scaled` must be bit-zero (the probe
   /// construction requires Q_bg = 0 for linearity; see plan §5.2.1).
   /// Any nonzero entry aborts.  Resolves the (face_idx, Elem1No) entry,
   /// computes F_h = nApNm1 * I_self, and subtracts
   /// `w * shape1(i) * F_h[c]` at the caller's DOFs.
   void AddBoundaryFaceRhs(int face_idx,
                           const real_t *I_self,
                           const real_t *bulk_bg_scaled,
                           real_t w,
                           const real_t *shape1,
                           int ndof, int dof_offset1,
                           int ndof_total,
                           Vector &rhs) const;

   bool IsInitialized() const { return initialized_; }

   /// Test-visibility accessor keyed by (face_idx, caller_elem_id).
   /// Aborts on miss.
   const FaceEntry &GetEntry(int face_idx, int caller_elem_id) const;

   std::size_t NumEntries() const { return entries_.size(); }

   /// Test-only: check whether an entry exists for (face_idx, caller_elem_id).
   bool HasEntry(int face_idx, int caller_elem_id) const;

   /// Test-only introspection: expose the face-keyed bookkeeping.
   const std::vector<FaceEntry> &GetEntries() const { return entries_; }

   /// Encoded (face_idx, caller_elem_id) key.  R3-006: uint32_t cast
   /// prevents sign-extension from negative caller_elem_id corrupting
   /// the upper 32 bits when bitwise-or'ed with the shifted face_idx.
   static long long EncodeKey(int face_idx, int caller_elem_id)
   {
      return (static_cast<long long>(face_idx) << 32) |
             static_cast<long long>(
                static_cast<std::uint32_t>(caller_elem_id));
   }

   // --- Internal helpers (exposed for targeted unit tests) --------------

   /// Topology-based outward normal + tangents + face area + cell volume
   /// for one (elem, local_side) tet face.  Algorithm per plan §5.1.4.
   static void ComputeCellFaceFrame(const Mesh &mesh,
                                    int elem, int local_side,
                                    real_t normal[3],
                                    real_t tangent1[3],
                                    real_t tangent2[3],
                                    real_t &surface_area,
                                    real_t &cell_volume);

   /// Build interior-face matrices in the GLOBAL frame:
   ///   nApNm1 = T * A_plus  * Tinv
   ///   nAmNm1 = T * A_minus * Tinv
   /// where (T, Tinv) are built from the supplied topology frame.
   static void BuildInteriorMatrices(const real_t normal[3],
                                     const real_t tangent1[3],
                                     const real_t tangent2[3],
                                     const GodunovFlux &flux,
                                     DenseMatrix &nApNm1,
                                     DenseMatrix &nAmNm1);

   /// Build boundary-face matrix by probing MFEM's runtime BC flux on 9
   /// basis states with bulk_bg = 0.  R4-003: (t1, t2) removed from
   /// signature -- dead parameters in the probe construction.
   static void BuildBoundaryMatricesGodunov(const real_t normal[3],
                                            const GodunovFlux &flux,
                                            DenseMatrix &nApNm1);
   static void BuildBoundaryMatricesGamma(const real_t normal[3],
                                          const GodunovFlux &flux,
                                          DenseMatrix &nApNm1);
   static void BuildBoundaryMatricesAbsorbing(const real_t normal[3],
                                              const GodunovFlux &flux,
                                              DenseMatrix &nApNm1);

   /// Asserts bit-equality between the hard-coded FACE2NODES_MFEM[4][3]
   /// and mfem::Geometry::Constants<TETRAHEDRON>::FaceVert[4][3].
   /// Guards against MFEM tet-table drift across releases.
   static void AssertMFEMFace2NodesTable();

   /// MFEM reference-tet face-to-vertex table, verified against
   /// `mfem::Geometry::Constants<TETRAHEDRON>::FaceVert` in
   /// `mfem/fem/geom.cpp:987-988`.
   static constexpr int FACE2NODES_MFEM[4][3] = {
      /* face 0 */ {1, 2, 3},
      /* face 1 */ {0, 3, 2},
      /* face 2 */ {0, 1, 3},
      /* face 3 */ {0, 2, 1}
   };

private:
   std::vector<FaceEntry> entries_;

   /// (face_idx, caller_elem_id) -> entry index.  See EncodeKey.
   std::unordered_map<long long, int> face_elem_to_entry_;

   /// face_idx -> entry index, boundary-only.  AddBoundaryFaceRhs does
   /// not know the owning element at the call site, so we keep a
   /// parallel O(1) lookup table populated at Init time.
   std::unordered_map<int, int> face_to_bdr_entry_;

   bool initialized_ = false;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PRECOMPUTED_FACE_FLUXES_HPP
