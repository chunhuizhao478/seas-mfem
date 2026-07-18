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

#ifndef MFEM_SEAS_BIMATERIAL_WAVE_OPERATOR_HPP
#define MFEM_SEAS_BIMATERIAL_WAVE_OPERATOR_HPP

#include "wave_operator.hpp"
#include "godunov_flux_bimaterial.hpp"  // BimaterialFlux per-face dispatch
#include "godunov_flux_pool.hpp"        // per-element GodunovFlux pool
#include "heterogeneous_material.hpp"   // MaterialField

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Heterogeneous (bi-material / `interior_flux="matrix"`) DG wave
/// operator (Phase 13).
///
/// `WaveOperator<MeshType>` is the SCALAR (homogeneous Godunov) operator — the
/// exact class the standalone `tpv{102,104,205}_driver.cpp` / `seas_bp5_full`
/// and the `interior_flux="scalar"` branch of `spatial_dyn_driver.cpp`
/// instantiate.  `BimaterialWaveOperator<MeshType>` is a SEPARATE subclass that
/// holds the per-element flux pool, the per-face bi-material flux matrices, and
/// the cross-rank neighbour-material map, and **overrides** the handful of
/// material-/flux-dependent virtual hooks so that EVERY material access goes
/// through per-element material.  Because the scalar placeholder logic lives
/// nowhere in this class, a missed material site cannot silently run the
/// `(1,1,1)` sentinel — the inherited `flux_` is poisoned with `(1,1,1)` (a
/// valid-but-unphysical triple) so any stray bare-`flux_` use becomes a hard
/// parity failure on the fault-bearing C-6 test rather than silent wrong
/// physics.
///
/// The ~7000-line DG/ADER skeleton (`Mult`, `AdvanceADER`, `ComputeADER*`,
/// `ComputeVolumeRHS`, `ComputeFaceFluxRHS`, …) lives ONCE in `WaveOperator`
/// and is inherited unchanged; this class injects per-element material behavior
/// only through the overridden hooks.
///
/// Reference: hrs-ref `wave_operator.inl` (the single-class form Phase 9
/// ported); Phase 13 splits it into two classes so the placeholder-leak class
/// is impossible by construction.
template <typename MeshType = Mesh>
class BimaterialWaveOperator : public WaveOperator<MeshType>
{
   using Base = WaveOperator<MeshType>;

   // Phase 13: the moved method bodies (ctor, builders, overrides) reference
   // these protected base members by their unqualified names.  In a class
   // template derived from a dependent base, unqualified lookup does NOT find
   // base members, so we pull them in with using-declarations (equivalent to
   // `this->member`) to keep the moved bodies byte-for-byte as they were in
   // `WaveOperator`.
   using Base::mesh_;
   using Base::ne_;
   using Base::ndof_per_el_;
   using Base::ndof_total_;
   using Base::bc_;
   using Base::fault_interior_faces_;
   using Base::fault_shared_faces_;      // IMPL-5: shared-fault disjointness assert
   using Base::face_bdr_attr_;
   using Base::shared_face_bdr_attr_;
   using Base::mixed_flux_mode_;
   using Base::central_flux_face_set_;   // Phase 2: central-flux build iterates it
   using Base::mf_on_;                   // Phase 3: per-face central dispatch gate
   using Base::cfl_rk_aware_;            // P3-4: interim central+ADER abort guard
   using Base::mixed_flux_contrast_tol_; // Part A: central-flux contrast guard tol

public:
   /// @brief Construct the heterogeneous operator from a `MaterialField`
   /// (Mode::Constant, Mode::Coefficient, or Mode::GridFunction).
   ///
   /// (Cross-rank Phase 1) GridFunction mode is now accepted: every material
   /// read routes through the uniform `MaterialAtLocal_`/`MaterialAtNbr_`
   /// accessors, whose GridFunction branch reads the source GFs (`GetValue`
   /// locally, shape×`FaceNbrData` across a seam) — so the per-element flux
   /// pool has a well-defined centroid evaluation in all three modes.
   ///
   /// Delegates to the SCALAR base ctor with the **valid sentinel triple
   /// `(lambda, mu, rho) = (1, 1, 1)`** (NOT NaN — `GodunovFlux` asserts
   /// `mu>0 && rho>0 && lambda+2*mu>0`, and `NaN>0` is false, so a NaN seed
   /// would abort EVERY matrix run).  On this object `flux_` is dead: every
   /// material access routes through the overridden `FluxForElem_(e)` →
   /// `owned_flux_pool_->At(e)`.  The `(1,1,1)` value is deliberately
   /// unphysical so any stray bare-`flux_` site is numerically visible on the
   /// fault-bearing parity test (the tripwire).  After base construction it
   /// builds the per-element pool, the per-element CFL cache, the cross-rank
   /// neighbour map, and the per-face bi-material flux matrices.
   BimaterialWaveOperator(MeshType &mesh, int order,
                          const MaterialField &material,
                          const BoundaryConfig &bc);

   ~BimaterialWaveOperator() override = default;

   /// True on this class (overrides the base's `false`).
   bool UsesGodunovFluxPool() const override
   { return owned_flux_pool_.get() != nullptr; }

   /// Per-element CFL data (production: also consumed by the LTS Phase-0
   /// `--lts-report` clustering, plus `test_phaseh_wave_operator_constant_parity`
   /// and `test_wave_operator`).  Override the base uniform-material accessors
   /// with this class's true per-element heterogeneous arrays, so a driver
   /// holding a base `WaveOperator&` sees the matrix-path values via dispatch.
   const std::vector<std::array<real_t, 3>> &GetPerElementMaterial() const override
   { return per_elem_lmr_; }
   const std::vector<real_t> &GetPerElementCflLength() const override
   { return per_elem_h_; }
   const std::unordered_map<int, std::array<real_t, 3>> &
   GetSharedFaceNeighbourMaterial() const
   { return shared_face_neighbour_material_; }
   const std::vector<std::array<std::array<mfem::DenseMatrix, 2>, 2>> &
   GetPerFaceBimaterialFlux() const
   { return per_face_bimaterial_flux_; }
   /// Phase 2 (mixed-flux): side-symmetric per-face central-flux matrices.
   /// `[mesh_face][0]` = ½·A_self (multiplier of Q_self), `[1]` = ½·A_nbr.
   /// Only central-set faces are present (empty unless mixed flux is enabled).
   const std::unordered_map<int, std::array<mfem::DenseMatrix, 2>> &
   GetPerFaceCentralFlux() const
   { return per_face_central_flux_; }
   std::size_t GetPhaserDispatchCount() const { return phaser_dispatch_count_; }
   void ResetPhaserDispatchCount() const { phaser_dispatch_count_ = 0; }

   /// (Cross-rank Phase 1) Number of MATERIAL `ParGridFunction::ExchangeFaceNbrData`
   /// collectives performed at construction: 3 in Mode::GridFunction (one per
   /// source GF), 0 in Mode::Constant / Mode::Coefficient (globally evaluable —
   /// NO material MPI).  Consumed by the `accessor_coeff_no_material_mpi` gate.
   std::size_t GetMaterialGfExchangeCount() const
   { return material_gf_exchange_count_; }

   /// (Cross-rank Phase 3) Number of SHARED (cross-rank) central faces this rank
   /// reclassified central -> upwind in the last `BuildPerFaceCentralFluxMatrices_`
   /// (the contrast guard).  Per-rank: a seam face shared by two ranks is counted
   /// by BOTH (the guard predicate is symmetric, so both agree).  0 when the guard
   /// is off (tol<0) or no shared corridor face has a strong contrast.
   std::size_t GetNReclassShared() const { return n_reclass_shared_; }

   /// CFL via the per-element `per_elem_h_` / `per_elem_lmr_` walk +
   /// `MPI_Allreduce(MIN)` (overrides the scalar `h_min_ / flux_.GetCp()`).
   real_t ComputeMaxDt(real_t cfl) const override;

   /// (PLAN Phase 3 — lifts R-003) Enable mixed flux on the matrix path:
   /// delegate to the base (validate, rank-consistency mode check,
   /// `BuildCentralFluxFaceSet_`, set `mf_on_`/`mixed_flux_mode_`), then
   /// precompute the per-face central matrices (`BuildPerFaceCentralFluxMatrices_`).
   /// The bi-material CENTRAL flux `½(A_self·Q_self + A_nbr·Q_nbr)` is then
   /// dispatched per-face on the central-set faces alongside the bi-material
   /// Godunov upwind elsewhere (drdg3d `get_flux` structure).
   void SetMixedFluxMode(MixedFluxMode m) override;

   /// (Part B / B2) Attach the fault flux AND affirm per-side-A conversion.  The
   /// matrix operator converts the per-side imposed state to the bulk flux with
   /// per-side A (FluxForElem_(elem_plus/elem_minus) at the fault conversion sites
   /// wave_operator.inl:3048 (Mult), :4301 (ADER), :3715 (shared)), so the
   /// bi-material-fault homogeneity guards in FaultFaceFlux are safe to relax.
   /// Verified: B0 (seas_test_bimaterial_fault_riemann +
   /// debug_document/tpv6_debug_document/bimaterial_fault_verification_2026-06-06.md).
   void SetFaultFlux(FaultFaceFlux *ff) override
   {
      Base::SetFaultFlux(ff);
      if (ff != nullptr) { ff->SetPerSideFluxApplied(true); }
   }

   /// (Part B / B1) Overwrite each fault DOF's per-side impedances with the material
   /// just inside each side of the fault (eps-offset rule, R-001/R-102).  Local
   /// 2-sided fault faces only; SHARED (cross-rank) fault faces keep the
   /// driver-seeded single-material values (R-101 peer eps-offset exchange deferred).
   void AssignFaultSidePerMaterialImpedances(
      std::vector<DOFData> &dof_data) const override;

   /// DEPRECATED (Cross-rank Phase 5) — NO BEHAVIORAL EFFECT.  Formerly affirmed
   /// the material is seam-continuous to permit the central (mixed) flux on a
   /// SHARED face despite the local-side neighbour-material stub (R-004).  The
   /// cross-rank exchange (Phase 2) now reads the TRUE peer material, so the
   /// central build no longer needs (or consults) this affirmation — the abort it
   /// used to gate is gone.  The setter is retained ONLY for `[material].
   /// seam_continuous` config back-compat (the flag is stored but never read);
   /// strong-contrast safety is handled by the contrast guard
   /// (`mixed_flux_contrast_tol >= 0`), not this flag.
   void SetSeamContinuous(bool v) { seam_continuous_ = v; }

   /// Precomputed (scalar Godunov) face fluxes are incompatible with the
   /// bi-material interior-face Riemann solve (REVIEW R-006): `Init` would
   /// build per-face flux tables from the inherited `(1,1,1)` sentinel `flux_`
   /// AND the precomputed dispatch bypasses the per-face bi-material
   /// `InteriorFaceFlux_`.  Abort on enable (mirrors `SetMixedFluxMode`); a
   /// no-op disable is allowed so generic teardown paths stay valid.
   void UsePrecomputedFaceFluxes(bool enable) override
   {
      MFEM_VERIFY(!enable,
                  "BimaterialWaveOperator::UsePrecomputedFaceFluxes: "
                  "precomputed face fluxes are scalar-only; the heterogeneous "
                  "interior_flux=\"matrix\" path already replaces the "
                  "interior-face flux with the bi-material Riemann solve.  Use "
                  "interior_flux=\"scalar\" (WaveOperator) for precomputed "
                  "face fluxes.");
   }

protected:
   /// Per-element material for the volume Jacobian, the boundary-face flux,
   /// and the fault imposed-state flux.
   const GodunovFlux &FluxForElem_(int e) const override
   { return owned_flux_pool_->At(e); }

   /// Interior non-fault flux for a fully-local (2-sided) face — per-face
   /// bi-material matrices (side 0 = Elem1, side 1 = Elem2).
   void InteriorFaceFlux_(int mesh_face, const real_t *Q_self,
                          const real_t *Q_nbr, const real_t *nor,
                          real_t *F_h_e1, real_t *F_h_e2) const override;

   /// Interior non-fault flux for a SHARED face — local (Elem1, side 0) only.
   void SharedInteriorFaceFlux_(int mesh_face, const real_t *Q_self,
                                const real_t *Q_nbr, const real_t *nor,
                                real_t *F_h) const override;

   /// ADER CK-recursion element Jacobian apply via per-element star matrices.
   void ApplyElementJacobian_(int dir, const Vector &X, Vector &Y,
                              real_t sign) const override;

   /// (LTS Phase 2) Element-subset twin of the above — per-element star matrices
   /// restricted to `elems[0..n)`.  MUST override the base (which would apply the
   /// dead (1,1,1) sentinel `flux_` star matrix on the matrix path).  Bit-
   /// identical to `ApplyElementJacobian_` over the full list.
   void ApplyElementJacobianElems_(int dir, const Vector &X, Vector &Y,
                                   real_t sign, const int *elems,
                                   int n) const override;

   // -----------------------------------------------------------------------
   // (Cross-rank Phase 1) ONE uniform material accessor, read by the flux
   // pool, the shared-face neighbour material, and (Phase 4) the per-side
   // fault — uniformly, locally AND across a partition seam.  Its ONLY
   // internal branch is the material REPRESENTATION (Constant / Coefficient /
   // GridFunction), never the problem.  No projection anywhere.  `protected`
   // so the np=2 gate (`tests/parallel/test_bimaterial_seam_material_np2.cpp`)
   // can expose them through a `TestableBimat` `using`-shim.
   // -----------------------------------------------------------------------

   /// Material `(lambda, mu, rho)` at a LOCAL element + reference IP.
   ///   - Constant / Coefficient: byte-exact with `MaterialField::EvalAt` at the
   ///     CALLER-supplied element transform `T` (the legacy read).
   ///   - GridFunction: `gf->GetValue(elem, ip)` per component (source GF; no
   ///     projection; `T` unused).
   /// (P1-001) The caller passes its OWN `ElementTransformation &T` (e.g. the
   /// flux-pool loop's per-element transform, or a held `ftr->Elem1`) so the
   /// accessor NEVER touches `Mesh::GetElementTransformation`'s shared scratch —
   /// which `ParMesh::GetSharedFaceTransformations` aliases as `ftr->Elem1`.  This
   /// keeps the accessor composable with a live `FaceElementTransformations`
   /// (the Phase 2/4 shared loops) and matches the existing fault-loop pattern
   /// `material_->EvalAt(Elem1No, *ftr->Elem1, ip, ...)`.
   /// Asserts `material_` set and the result is physical (rho>0, lambda+2mu>0).
   void MaterialAtLocal_(int elem, mfem::ElementTransformation &T,
                         const mfem::IntegrationPoint &ip,
                         real_t &lam, real_t &mu, real_t &rho) const;

   /// Material `(lambda, mu, rho)` on the PEER (face-neighbour, `ftr->Elem2`)
   /// side of a SHARED face, at the peer reference IP `ip_peer`.
   ///   - Constant / Coefficient: same Eval at the peer's face-neighbour
   ///     element transform (`ftr->Elem2`) — globally evaluable, NO material MPI.
   ///   - GridFunction: shape×`FaceNbrData` interpolation of the source GFs'
   ///     face-neighbour ghost layer (exchanged once at construction).
   /// Asserts `ftr`/`ftr->Elem2` non-null and the result is physical.
   void MaterialAtNbr_(mfem::FaceElementTransformations *ftr,
                       const mfem::IntegrationPoint &ip_peer,
                       real_t &lam, real_t &mu, real_t &rho) const;

private:
   /// Walk every local element, evaluate the material at the reference centroid
   /// via `MaterialAtLocal_` (reads `material_`), fill `per_elem_lmr_` /
   /// `per_elem_h_`, and `Build()` the pool.  (Cross-rank Phase 1) Routed through
   /// the accessor so the pool is well-defined in all three material modes
   /// (GridFunction included); the Constant/Coefficient result is byte-exact with
   /// the legacy `material.EvalAt`.  (P1-006) No `MaterialField` parameter — the
   /// member `material_` is the single source of truth (set in the ctor first).
   void BuildGodunovFluxPool_();

   /// (Cross-rank Phase 1) One-time face-neighbour exchange setup (parallel
   /// only, collective on ALL ranks): refresh the mesh geometry peer transforms
   /// (idempotent — the base ctor already did it) and, in Mode::GridFunction
   /// ONLY, `ExchangeFaceNbrData` the three source material GFs so
   /// `MaterialAtNbr_` can read their face-neighbour ghost layer.  Bumps
   /// `material_gf_exchange_count_` once per GF exchanged (0 for Constant /
   /// Coefficient).  Must run BEFORE `BuildGodunovFluxPool_` /
   /// `ExchangeBiMaterialNeighbours_`.
   void SetupMaterialFaceNbrExchange_();

   /// Populate `shared_face_neighbour_material_` with each shared face's
   /// neighbour-side material (LOCAL-side stub; correct for seam-continuous /
   /// depth-only materials — see the body comment).
   void ExchangeBiMaterialNeighbours_();

   /// Precompute per-(face, side) bi-material flux matrices into
   /// `per_face_bimaterial_flux_`.
   void BuildPerFaceBimaterialFluxMatrices_();

   /// (PLAN Phase 2) Centralised per-face operand derivation for a fully-local
   /// 2-sided interior face: returns `(flux_self=Elem1, flux_nbr=Elem2, centroid
   /// unit normal)` so the Godunov and central per-face builds derive
   /// byte-identical operands/orientation.  Returns false for 1-sided / shared /
   /// non-fault-boundary faces (the caller `continue`s); does NOT itself skip
   /// fault faces (the caller does).  The returned pointers alias the owned flux
   /// pool (stable for the operator's lifetime).  `const`: `mesh_` is a reference
   /// member, so the non-const `Mesh::GetFaceElementTransformations` is callable
   /// here.  (IMPL-1) No per-side swap parameter — every caller uses
   /// self=Elem1/nbr=Elem2 (the upwind build swaps internally via
   /// `compose_side`; the central build is side-symmetric).
   bool ResolveFaceFluxOperands_(int mesh_face,
                                 const GodunovFlux *&flux_self,
                                 const GodunovFlux *&flux_nbr,
                                 real_t nor_out[3]) const;

   /// (IMPL-3) Centroid unit normal for a face transformation — the single
   /// source of the GetCenter -> SetAllIntPoints -> CalcOrtho -> normalize
   /// sequence used by the Godunov per-face build (Pass 2), the interior operand
   /// derivation, and the central shared-face build (so the Godunov and central
   /// normals cannot silently desync).
   void CentroidUnitNormal_(FaceElementTransformations *ftr,
                            real_t nor_out[3]) const;

   /// (PLAN Phase 2, BUG-3) Precompute the side-symmetric per-face CENTRAL
   /// flux matrices (`½·A_self`, `½·A_nbr`) for the faces in
   /// `central_flux_face_set_`, into `per_face_central_flux_`.  Called by
   /// `SetMixedFluxMode` (Phase 3); a no-op when the central set is empty
   /// (`mixed_flux=none`).
   void BuildPerFaceCentralFluxMatrices_();

   /// Per-element variant of the file-local `ApplyJacobianPerDOF`: applies
   /// `FluxForElem_(e).GetReferenceStarMatrix(dir)` to element `e`'s DOFs only.
   void ApplyJacobianPerElementDOF_(int dir, const Vector &X, Vector &Y,
                                    real_t sign) const;

   /// Non-owning pointer to the MaterialField (caller owns it; must outlive).
   const MaterialField *material_ = nullptr;

   /// (Cross-rank Phase 1) Count of MATERIAL GridFunction face-nbr exchanges
   /// performed at construction (3 in GridFunction mode, 0 otherwise).  Set once
   /// by `SetupMaterialFaceNbrExchange_`; read by `GetMaterialGfExchangeCount`.
   std::size_t material_gf_exchange_count_ = 0;

   /// (Cross-rank Phase 3) Per-rank count of SHARED central faces reclassified to
   /// upwind by the contrast guard in the last `BuildPerFaceCentralFluxMatrices_`.
   /// Set by that method; read by `GetNReclassShared`.
   std::size_t n_reclass_shared_ = 0;

   /// Owned per-element flux cache (the bi-material pool).
   std::unique_ptr<GodunovFluxPool> owned_flux_pool_;

   /// Per-element (lambda, mu, rho) cache (one entry per local element).
   std::vector<std::array<real_t, 3>> per_elem_lmr_;

   /// Per-element CFL length cache (inscribed diameter for tets).
   std::vector<real_t> per_elem_h_;

   /// Bi-material shared-face neighbour material map (local-side stub).
   std::unordered_map<int, std::array<real_t, 3>> shared_face_neighbour_material_;

   /// Per-(mesh face index, side) precomputed bi-material flux matrices.
   /// `[face][side][0]` = fluxLocal, `[face][side][1]` = fluxNeighbor.
   std::vector<std::array<std::array<mfem::DenseMatrix, 2>, 2>>
      per_face_bimaterial_flux_;

   /// (PLAN Phase 2, BUG-3) Side-symmetric per-face CENTRAL flux matrices.
   /// `per_face_central_flux_[mesh_face] = { ½·A_self, ½·A_nbr }` — ONE pair
   /// per face (NO per-side swap; the dispatch deposits the single-valued
   /// F* = ½A_self·Q_self + ½A_nbr·Q_nbr identically to both sides).  Only
   /// central-set faces present; empty unless mixed flux is enabled.
   std::unordered_map<int, std::array<mfem::DenseMatrix, 2>>
      per_face_central_flux_;

   /// DEPRECATED (Cross-rank Phase 5) — stored but NEVER READ.  Formerly gated the
   /// central build on `Mode::Coefficient` SHARED faces (R-004 stub era); the
   /// cross-rank exchange (Phase 2) removed its only consumer.  Kept for
   /// `[material].seam_continuous` config back-compat.  Default false.
   bool seam_continuous_ = false;

   /// Diagnostic: count of bi-material per-face flux applies.  `mutable`
   /// because the dispatch runs inside const Mult/ADER paths.
   mutable std::size_t phaser_dispatch_count_ = 0;
};

// Implementation in bimaterial_wave_operator.inl
#include "bimaterial_wave_operator.inl"

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BIMATERIAL_WAVE_OPERATOR_HPP
