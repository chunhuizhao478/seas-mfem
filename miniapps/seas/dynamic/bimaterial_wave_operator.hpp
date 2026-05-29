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
   using Base::face_bdr_attr_;
   using Base::shared_face_bdr_attr_;
   using Base::mixed_flux_mode_;

public:
   /// @brief Construct the heterogeneous operator from a `MaterialField`
   /// (Mode::Constant or Mode::Coefficient; GridFunction rejected).
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

   /// Test-only accessors (moved verbatim from `WaveOperator`; consumed by
   /// `test_phaseh_wave_operator_constant_parity` and `test_wave_operator`).
   const std::vector<std::array<real_t, 3>> &GetPerElementMaterial() const
   { return per_elem_lmr_; }
   const std::vector<real_t> &GetPerElementCflLength() const
   { return per_elem_h_; }
   const std::unordered_map<int, std::array<real_t, 3>> &
   GetSharedFaceNeighbourMaterial() const
   { return shared_face_neighbour_material_; }
   const std::vector<std::array<std::array<mfem::DenseMatrix, 2>, 2>> &
   GetPerFaceBimaterialFlux() const
   { return per_face_bimaterial_flux_; }
   std::size_t GetPhaserDispatchCount() const { return phaser_dispatch_count_; }
   void ResetPhaserDispatchCount() const { phaser_dispatch_count_ = 0; }

   /// CFL via the per-element `per_elem_h_` / `per_elem_lmr_` walk +
   /// `MPI_Allreduce(MIN)` (overrides the scalar `h_min_ / flux_.GetCp()`).
   real_t ComputeMaxDt(real_t cfl) const override;

   /// Mixed flux is scalar-only (R-003): abort on any non-None mode.
   void SetMixedFluxMode(MixedFluxMode m) override;

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

private:
   /// Walk every local element, evaluate `material.EvalAt` at the reference
   /// centroid, fill `per_elem_lmr_` / `per_elem_h_`, and `Build()` the pool.
   void BuildGodunovFluxPool_(const MaterialField &material);

   /// Populate `shared_face_neighbour_material_` with each shared face's
   /// neighbour-side material (LOCAL-side stub; correct for seam-continuous /
   /// depth-only materials — see the body comment).
   void ExchangeBiMaterialNeighbours_();

   /// Precompute per-(face, side) bi-material flux matrices into
   /// `per_face_bimaterial_flux_`.
   void BuildPerFaceBimaterialFluxMatrices_();

   /// Per-element variant of the file-local `ApplyJacobianPerDOF`: applies
   /// `FluxForElem_(e).GetReferenceStarMatrix(dir)` to element `e`'s DOFs only.
   void ApplyJacobianPerElementDOF_(int dir, const Vector &X, Vector &Y,
                                    real_t sign) const;

   /// Non-owning pointer to the MaterialField (caller owns it; must outlive).
   const MaterialField *material_ = nullptr;

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

   /// Diagnostic: count of bi-material per-face flux applies.  `mutable`
   /// because the dispatch runs inside const Mult/ADER paths.
   mutable std::size_t phaser_dispatch_count_ = 0;
};

// Implementation in bimaterial_wave_operator.inl
#include "bimaterial_wave_operator.inl"

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BIMATERIAL_WAVE_OPERATOR_HPP
