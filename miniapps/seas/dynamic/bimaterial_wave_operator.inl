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

// ===========================================================================
// Implementation of BimaterialWaveOperator<MeshType>.
//
// Included at the END of bimaterial_wave_operator.hpp, INSIDE
// `namespace mfem { namespace seas { ... } }`.  All bodies are moved verbatim
// from the single-class form Phase 9 ported into wave_operator.inl; Phase 13
// splits them into this subclass so the scalar `WaveOperator` is bi-material-
// free and the `(1,1,1)` placeholder leak is impossible by construction.
// ===========================================================================

// ---------------------------------------------------------------------------
// Phase 13 ctor — heterogeneous-material.  Delegates to the SCALAR base ctor
// with the valid sentinel triple (1,1,1) (NOT NaN — GodunovFlux asserts
// mu>0/rho>0/lambda+2*mu>0, and NaN>0 is false, so a NaN seed would abort
// every matrix run).  On this object the inherited `flux_` is dead: every
// material access routes through the overridden FluxForElem_(e) ->
// owned_flux_pool_->At(e).  The unphysical (1,1,1) value makes any stray bare-
// `flux_` site numerically visible on the fault-bearing C-6 parity test.
// ---------------------------------------------------------------------------
template <typename MeshType>
BimaterialWaveOperator<MeshType>::BimaterialWaveOperator(
   MeshType &mesh, int order, const MaterialField &material,
   const BoundaryConfig &bc)
   : Base(mesh, order, real_t(1.0), real_t(1.0), real_t(1.0), bc)
{
   MFEM_VERIFY(material.mode == MaterialField::Mode::Constant
               || material.mode == MaterialField::Mode::Coefficient,
               "BimaterialWaveOperator(MaterialField): material.mode must be "
               "Constant or Coefficient; got Mode::GridFunction ("
               << static_cast<int>(material.mode) << ").  GridFunction "
               "mode is consumed via At(elem, dof, ...), which has no "
               "well-defined centroid evaluation needed by the per-"
               "element flux pool.");

   material_ = &material;

   BuildGodunovFluxPool_(material);
   ExchangeBiMaterialNeighbours_();
   BuildPerFaceBimaterialFluxMatrices_();
}

// ---------------------------------------------------------------------------
// Phase H.1 helper — fill per_elem_lmr_ / per_elem_h_, then Build()
// the owned flux pool.
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::BuildGodunovFluxPool_(
   const MaterialField &material)
{
   per_elem_lmr_.assign(static_cast<size_t>(ne_), std::array<real_t, 3>{0, 0, 0});
   per_elem_h_.assign(static_cast<size_t>(ne_), real_t(0));

   for (int e = 0; e < ne_; ++e)
   {
      ElementTransformation *T = mesh_.GetElementTransformation(e);
      const Geometry::Type   gtype = mesh_.GetElementBaseGeometry(e);
      const IntegrationPoint &ip   = Geometries.GetCenter(gtype);
      real_t lam, mu, rho;
      material.EvalAt(e, *T, ip, lam, mu, rho);
      MFEM_VERIFY(rho > 0.0,
                  "WaveOperator(MaterialField): rho (" << rho
                  << ") must be > 0 at element " << e);
      MFEM_VERIFY(lam + 2.0 * mu > 0.0,
                  "WaveOperator(MaterialField): (lambda + 2*mu) ("
                  << (lam + 2.0 * mu) << ") must be > 0 at element "
                  << e << " (lambda=" << lam << ", mu=" << mu << ").");
      per_elem_lmr_[e] = {lam, mu, rho};

      // Mirror the scalar ctor's per-element CFL-length formula
      // (inscribed diameter for tets, vol^{1/dim} for non-tets).
      // Byte-exact equality with the scalar ctor's `h_min_` loop on
      // Mode::Constant input.
      const real_t vol = mesh_.GetElementVolume(e);
      real_t h;
      if (gtype == Geometry::TETRAHEDRON)
      {
         Array<int> vert;
         mesh_.GetElementVertices(e, vert);
         Vector v0(mesh_.GetVertex(vert[0]), 3);
         Vector v1(mesh_.GetVertex(vert[1]), 3);
         Vector v2(mesh_.GetVertex(vert[2]), 3);
         Vector v3(mesh_.GetVertex(vert[3]), 3);
         auto tri_area = [](const Vector &a, const Vector &b,
                            const Vector &c) -> real_t {
            real_t e1[3] = {b(0)-a(0), b(1)-a(1), b(2)-a(2)};
            real_t e2[3] = {c(0)-a(0), c(1)-a(1), c(2)-a(2)};
            real_t cx = e1[1]*e2[2] - e1[2]*e2[1];
            real_t cy = e1[2]*e2[0] - e1[0]*e2[2];
            real_t cz = e1[0]*e2[1] - e1[1]*e2[0];
            return 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
         };
         const real_t A_total = tri_area(v0, v1, v2) + tri_area(v0, v2, v3)
                              + tri_area(v0, v3, v1) + tri_area(v1, v2, v3);
         h = (A_total > 0) ? 6.0 * vol / A_total
                           : std::pow(vol, 1.0 / 3.0);
      }
      else
      {
         h = std::pow(vol, 1.0 / mesh_.Dimension());
      }
      per_elem_h_[e] = h;
   }

   owned_flux_pool_ = std::make_unique<GodunovFluxPool>();
   owned_flux_pool_->Build(ne_, per_elem_lmr_, /*dedup_sig_figs=*/6);
}

// ---------------------------------------------------------------------------
// Phase H.4 helper — populate shared_face_neighbour_material_.
//
// REVIEW R-004 (stale-comment correction): Mode::Coefficient is now REACHABLE
// (Phase 9 wired the matrix path + removed the ctor abort), so the old "only
// Mode::Constant is reachable" contract is no longer true.  This body is still
// a LOCAL-SIDE STUB — it stores the local element's own material as the
// neighbour's (no MPI exchange).  That is correct ONLY when the neighbour's
// material equals the local material at the seam:
//   - Mode::Constant: always (every element shares the constants).
//   - Mode::Coefficient that is seam-continuous (e.g. depth-only, TPV31): the
//     centroid material agrees across the seam, so the stub is correct.
//   - Mode::Coefficient with LATERAL variation across a partition seam: WRONG
//     — the genuine peer-rank neighbour material is needed.  A real
//     MPI_Allgatherv exchange (key shared faces, pair, store the peer's
//     per_elem_lmr_) is unimplemented; the het ctor WARNs (below) for the
//     Coefficient + shared-faces case so a parallel laterally-heterogeneous
//     run does not silently use the wrong seam material.
//   - On serial Mesh this is a no-op (no shared faces).
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::ExchangeBiMaterialNeighbours_()
{
   shared_face_neighbour_material_.clear();
   if constexpr (!IsParallelMesh<MeshType>::value)
   {
      return;
   }
#ifdef MFEM_USE_MPI
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      auto &pmesh = static_cast<ParMesh &>(mesh_);
      const int n_shared = pmesh.GetNSharedFaces();
      for (int sf = 0; sf < n_shared; ++sf)
      {
         FaceElementTransformations *ftr =
            pmesh.GetSharedFaceTransformations(sf);
         if (!ftr) { continue; }
         const int local_elem = ftr->Elem1No;
         MFEM_ASSERT(local_elem >= 0 && local_elem < ne_,
                     "ExchangeBiMaterialNeighbours_: shared face " << sf
                     << " has Elem1No=" << local_elem
                     << " outside [0, " << ne_ << ").");
         shared_face_neighbour_material_[sf] = per_elem_lmr_[local_elem];
      }

      // REVIEW R-004: warn once (rank 0) if a genuinely heterogeneous
      // (Coefficient) material is used in parallel — the local-side stub
      // above uses the WRONG neighbour material at a partition seam where the
      // material varies laterally (it is correct only for seam-continuous /
      // depth-only materials like TPV31).  The real cross-rank exchange is
      // unimplemented.
      if (material_ != nullptr
          && material_->mode == MaterialField::Mode::Coefficient
          && n_shared > 0)
      {
         int rank = 0;
         MPI_Comm_rank(pmesh.GetComm(), &rank);
         if (rank == 0)
         {
            mfem::out << "[wave_operator] WARNING: ExchangeBiMaterialNeighbours_"
                         " is a LOCAL-SIDE stub — a parallel run with a "
                         "laterally-varying (Coefficient) material will use the "
                         "WRONG neighbour material at partition seams (correct "
                         "only for depth-only / seam-continuous materials like "
                         "TPV31).  The real cross-rank exchange is not yet "
                         "implemented.\n";
         }
      }
   }
#endif
}

// ---------------------------------------------------------------------------
// Phase 9 (Stage B) — per-element CK Jacobian apply (ported from hrs-ref).
// Loops elements, fetches FluxForElem_(e).GetReferenceStarMatrix(dir), and
// applies it ONLY to element e's ndof_per_el_ DOFs.  Bit-identical to the
// scalar ApplyJacobianPerDOF on Mode::Constant (same A on every element).
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::ApplyJacobianPerElementDOF_(
   int dir, const Vector &X, Vector &Y, real_t sign) const
{
   MFEM_ASSERT(dir >= 0 && dir < 3,
               "ApplyJacobianPerElementDOF_: dir must be in {0,1,2}, got "
               << dir);
   MFEM_ASSERT(X.Size() == NUM_STATE * ndof_total_,
               "ApplyJacobianPerElementDOF_: X size mismatch");
   MFEM_ASSERT(Y.Size() == NUM_STATE * ndof_total_,
               "ApplyJacobianPerElementDOF_: Y size mismatch");

   const real_t *Xd = X.GetData();
   real_t *Yd = Y.GetData();

   for (int e = 0; e < ne_; ++e)
   {
      const DenseMatrix &A = FluxForElem_(e).GetReferenceStarMatrix(dir);
      const int base = e * ndof_per_el_;
      for (int c = 0; c < NUM_STATE; ++c)
      {
         real_t *Yc = Yd + c * ndof_total_ + base;
         for (int cp = 0; cp < NUM_STATE; ++cp)
         {
            const real_t a = A(c, cp);
            if (a == 0.0) { continue; }
            const real_t w = sign * a;
            const real_t *Xcp = Xd + cp * ndof_total_ + base;
            for (int i = 0; i < ndof_per_el_; ++i)
            {
               Yc[i] += w * Xcp[i];
            }
         }
      }
   }
}

// ---------------------------------------------------------------------------
// Phase 2 (mixed-flux): centroid unit normal for a face transformation.
// (IMPL-3) Single source of the GetCenter -> SetAllIntPoints -> CalcOrtho ->
// normalize sequence, shared by the Godunov per-face build (Pass 2), the
// interior operand derivation (ResolveFaceFluxOperands_), and the central
// shared-face build — so the Godunov and central per-face normals cannot
// silently desync.
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::CentroidUnitNormal_(
   FaceElementTransformations *ftr, real_t nor_out[3]) const
{
   const Geometry::Type gtype = ftr->GetGeometryType();
   const IntegrationPoint &ip0 = Geometries.GetCenter(gtype);
   ftr->SetAllIntPoints(&ip0);
   Vector nor_vec(3);
   CalcOrtho(ftr->Face->Jacobian(), nor_vec);
   const real_t nor_len = nor_vec.Norml2();
   MFEM_VERIFY(nor_len > 0,
               "CentroidUnitNormal_: face has zero Jacobian normal at centroid.");
   nor_vec /= nor_len;
   nor_out[0] = nor_vec(0);
   nor_out[1] = nor_vec(1);
   nor_out[2] = nor_vec(2);
}

// ---------------------------------------------------------------------------
// Phase 2 (mixed-flux + bi-material): centralised per-face operand derivation.
// Returns (flux_self=Elem1, flux_nbr=Elem2, centroid unit normal) for a
// fully-local 2-sided interior face so the Godunov per-face build (Pass 1) and
// the central per-face build (BuildPerFaceCentralFluxMatrices_) derive IDENTICAL
// operands/orientation.  Returns false for 1-sided / shared / non-fault-boundary
// faces (the caller continues); does NOT itself skip fault faces (the caller
// does).  The returned pointers alias the owned flux pool (stable for the
// operator's lifetime).  (IMPL-1) No per-side swap parameter — every caller uses
// self=Elem1/nbr=Elem2; the upwind build swaps internally via compose_side, and
// the central build is side-symmetric.
// ---------------------------------------------------------------------------
template <typename MeshType>
bool BimaterialWaveOperator<MeshType>::ResolveFaceFluxOperands_(
   int mesh_face,
   const GodunovFlux *&flux_self, const GodunovFlux *&flux_nbr,
   real_t nor_out[3]) const
{
   FaceElementTransformations *ftr =
      mesh_.GetFaceElementTransformations(mesh_face);
   if (!ftr) { return false; }

   const int e1 = ftr->Elem1No;
   const int e2 = ftr->Elem2No;
   if (e2 < 0) { return false; }                          // 1-sided / shared
   if (face_bdr_attr_[mesh_face] != 0) { return false; }  // non-fault boundary

   CentroidUnitNormal_(ftr, nor_out);
   flux_self = &owned_flux_pool_->At(e1);
   flux_nbr  = &owned_flux_pool_->At(e2);
   return true;
}

// ---------------------------------------------------------------------------
// Phase 9 (Stage B) — precompute per-(face, side) bi-material flux matrices
// (ported from hrs-ref).  Skips fault and boundary faces.  Fully-local
// interior faces populate BOTH sides; shared faces populate ONLY side=0.
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::BuildPerFaceBimaterialFluxMatrices_()
{
   MFEM_VERIFY(owned_flux_pool_,
               "BuildPerFaceBimaterialFluxMatrices_: owned_flux_pool_ "
               "must be built first (call BuildGodunovFluxPool_).");

   const int n_faces = mesh_.GetNumFaces();
   per_face_bimaterial_flux_.assign(
      static_cast<size_t>(n_faces),
      std::array<std::array<DenseMatrix, 2>, 2>{});

   // Build the union of fault face sets so we skip fault faces without
   // depending on per-face boundary-attribute lookups (which may be 0
   // on a 2-sided interior fault face that lives entirely within the
   // mesh).
   std::set<int> fault_face_idx_set;
   for (int i = 0; i < fault_interior_faces_.Size(); ++i)
   {
      fault_face_idx_set.insert(fault_interior_faces_[i]);
   }

   // The per-(face, side) matrices share a single bi-material Riemann
   // state Q* per face.  Q* is defined by (n̂, L=Elem1, R=Elem2) and does
   // NOT depend on whose POV we compute from.  Both Elem1 and Elem2 apply
   // their OWN Jacobian A_self to the SAME Q*; side 0 has A_e1 baked in,
   // side 1 has A_e2.  The runtime apply on BOTH sides consumes
   // (Q_self=Q_e1, Q_nbr=Q_e2) — NOT swapped.  In the homogeneous limit
   // A_e1 = A_e2 so the two sides match `flux_.Interior(...)` to LU rounding.
   auto build_face_matrices = [&](const real_t *nor_unit,
                                  const GodunovFlux &flux_L,
                                  const GodunovFlux &flux_R,
                                  std::array<DenseMatrix, 2> &out_side0,
                                  std::array<DenseMatrix, 2> &out_side1)
   {
      const real_t n2 = nor_unit[0]*nor_unit[0] + nor_unit[1]*nor_unit[1]
                      + nor_unit[2]*nor_unit[2];
      MFEM_VERIFY(std::abs(n2 - 1.0) < 1e-10,
                  "BuildPerFaceBimaterialFluxMatrices_: face normal not "
                  "unit (|nor|^2 = " << n2 << ").");

      real_t t1[3], t2[3];
      GodunovFlux::BuildFrame(nor_unit, t1, t2);
      DenseMatrix T(NUM_STATE, NUM_STATE);
      DenseMatrix Tinv(NUM_STATE, NUM_STATE);
      GodunovFlux::BuildRotation(nor_unit, t1, t2, T);
      GodunovFlux::BuildRotationInverse(nor_unit, t1, t2, Tinv);

      DenseMatrix qGodL_FL, qGodN_FL;
      BimaterialFlux::BuildGodunovStateFaceLocal(
         flux_L.GetLambda(), flux_L.GetMu(), flux_L.GetRho(),
         flux_R.GetLambda(), flux_R.GetMu(), flux_R.GetRho(),
         qGodL_FL, qGodN_FL);

      auto compose_side = [&](const DenseMatrix &A_self_FL,
                              std::array<DenseMatrix, 2> &out_pair)
      {
         DenseMatrix tmp1(NUM_STATE, NUM_STATE);
         DenseMatrix tmp2(NUM_STATE, NUM_STATE);

         mfem::Mult(A_self_FL, qGodL_FL, tmp1);
         mfem::Mult(T, tmp1, tmp2);
         out_pair[0].SetSize(NUM_STATE, NUM_STATE);
         mfem::Mult(tmp2, Tinv, out_pair[0]);

         mfem::Mult(A_self_FL, qGodN_FL, tmp1);
         mfem::Mult(T, tmp1, tmp2);
         out_pair[1].SetSize(NUM_STATE, NUM_STATE);
         mfem::Mult(tmp2, Tinv, out_pair[1]);
      };

      compose_side(flux_L.GetAx(), out_side0);   // Elem1's POV (A_self = A_L)
      compose_side(flux_R.GetAx(), out_side1);   // Elem2's POV (A_self = A_R)
   };

   auto build_shared_face_side0 = [&](const real_t *nor_unit,
                                      const GodunovFlux &flux_local,
                                      const GodunovFlux &flux_nbr,
                                      std::array<DenseMatrix, 2> &out_pair)
   {
      BimaterialFlux::BuildPerFaceFluxMatricesGlobal(
         nor_unit, flux_local, flux_nbr,
         out_pair[0],   // fluxLocal (multiplier of Q_local)
         out_pair[1]);  // fluxNeighbor (multiplier of Q_nbr)
   };

   std::size_t n_int_built = 0;   // fully-local 2-sided interior faces
   std::size_t n_shr_built = 0;   // shared (ParMesh) interior faces

   // --- Pass 1: fully-local interior faces ---------------------------------
   for (int f = 0; f < n_faces; ++f)
   {
      if (fault_face_idx_set.count(f) > 0) { continue; }

      // Centralised operand derivation (self=Elem1, nbr=Elem2); returns false
      // (and we skip) for 1-sided / shared / non-fault-boundary faces.
      // Byte-identical to the previous inline derivation: same centroid normal
      // and the same owned-pool flux references.
      const GodunovFlux *flux_e1 = nullptr;
      const GodunovFlux *flux_e2 = nullptr;
      real_t nor[3];
      if (!ResolveFaceFluxOperands_(f, flux_e1, flux_e2, nor))
      {
         continue;
      }

      build_face_matrices(nor, *flux_e1, *flux_e2,
                          per_face_bimaterial_flux_[f][0],
                          per_face_bimaterial_flux_[f][1]);

      ++n_int_built;
   }

   // --- Pass 2: shared (cross-rank) faces, ParMesh only --------------------
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto &pmesh = static_cast<ParMesh &>(mesh_);
      const int n_shared = pmesh.GetNSharedFaces();
      for (int sf = 0; sf < n_shared; ++sf)
      {
         FaceElementTransformations *ftr =
            pmesh.GetSharedFaceTransformations(sf);
         if (!ftr) { continue; }

         const int mesh_face_idx = pmesh.GetSharedFace(sf);
         MFEM_VERIFY(mesh_face_idx >= 0 && mesh_face_idx < n_faces,
                     "BuildPerFaceBimaterialFluxMatrices_: shared face "
                     << sf << " -> mesh face " << mesh_face_idx
                     << " out of [0, " << n_faces << ").");

         const bool sf_fault =
            (sf < static_cast<int>(shared_face_bdr_attr_.size()))
            && (shared_face_bdr_attr_[sf] == bc_.fault_attr)
            && (bc_.fault_attr > 0);
         if (sf_fault) { continue; }

         const int local_elem = ftr->Elem1No;
         MFEM_ASSERT(local_elem >= 0 && local_elem < ne_,
                     "BuildPerFaceBimaterialFluxMatrices_: shared face "
                     << sf << " Elem1No=" << local_elem
                     << " outside [0, " << ne_ << ").");

         auto nbr_it = shared_face_neighbour_material_.find(sf);
         MFEM_VERIFY(nbr_it != shared_face_neighbour_material_.end(),
                     "BuildPerFaceBimaterialFluxMatrices_: shared face "
                     << sf << " missing entry in "
                     "shared_face_neighbour_material_.");
         const auto &lmr_nbr = nbr_it->second;
         GodunovFlux flux_nbr(lmr_nbr[0], lmr_nbr[1], lmr_nbr[2]);

         real_t nor[3];
         CentroidUnitNormal_(ftr, nor);   // (IMPL-3) shared centroid-normal helper

         const GodunovFlux &flux_local = owned_flux_pool_->At(local_elem);

         build_shared_face_side0(nor, flux_local, flux_nbr,
                                 per_face_bimaterial_flux_[mesh_face_idx][0]);

         ++n_shr_built;
      }
#endif
   }

   const std::size_t per_face_per_side_bytes =
      2 * NUM_STATE * NUM_STATE * sizeof(real_t);
   const std::size_t total_bytes =
      (2 * n_int_built + n_shr_built) * per_face_per_side_bytes;

   int print_rank = 0;
#ifdef MFEM_USE_MPI
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      auto &pmesh = static_cast<ParMesh &>(mesh_);
      MPI_Comm_rank(pmesh.GetComm(), &print_rank);
   }
#endif
   if (print_rank == 0)
   {
      mfem::out << "[wave_operator] BimaterialFlux precomputation:\n"
                << "  interior faces processed = " << n_int_built
                << "   (2 sides each)\n"
                << "  shared faces processed   = " << n_shr_built
                << "   (1 side each, Elem1 = local)\n"
                << "  total bytes              = " << total_bytes
                << " (per-rank)\n";

      constexpr std::size_t kSoftCapBytes = 4ULL * 1024 * 1024 * 1024;
      if (total_bytes > kSoftCapBytes)
      {
         mfem::out << "[wave_operator] WARNING: per-rank BimaterialFlux "
                   << "precomputation = " << total_bytes
                   << " bytes (> " << kSoftCapBytes
                   << " soft cap).  Partition more aggressively or "
                   << "enable a runtime-recompute fallback.\n";
      }
   }
}

// ---------------------------------------------------------------------------
// Phase 2 (mixed-flux + bi-material): precompute the side-symmetric per-face
// CENTRAL flux matrices for the fault-adjacent faces in central_flux_face_set_.
//
// Storage (PLAN BUG-3): ONE pair per face, per_face_central_flux_[mesh_face] =
// { centralE1 = 0.5*A_e1, centralE2 = 0.5*A_e2 }, NO per-side swap.  The Phase-3
// dispatch deposits the single-valued F* = centralE1.Q_e1 + centralE2.Q_e2
// identically to both sides.  Built from the SAME per-side operands as the
// Godunov matrices: ResolveFaceFluxOperands_ for interior faces; the sf-keyed
// shared-face path (mirrors BuildPerFaceBimaterialFluxMatrices_ Pass 2) for
// shared faces.
//
// NOT called from the ctor: SetMixedFluxMode (Phase 3) calls this after the
// base populates central_flux_face_set_.  With mixed_flux=none the set is
// empty so this is a no-op (the byte-exact mixed_flux=none path).
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::BuildPerFaceCentralFluxMatrices_()
{
   MFEM_VERIFY(owned_flux_pool_,
               "BuildPerFaceCentralFluxMatrices_: owned_flux_pool_ must be "
               "built first (call BuildGodunovFluxPool_).");

   per_face_central_flux_.clear();

   // Fault-face set for the debug disjointness assert (central_flux_face_set_
   // must exclude fault faces — the base BuildCentralFluxFaceSet_ guarantees it,
   // R-1100/R-1404).  (IMPL-5) Include SHARED fault faces too (mirroring the
   // base's fault-set construction at wave_operator.inl:1699-1710), so the
   // assert can catch a shared fault face wrongly present in the central set.
   std::set<int> fault_face_idx_set;
   for (int i = 0; i < fault_interior_faces_.Size(); ++i)
   {
      fault_face_idx_set.insert(fault_interior_faces_[i]);
   }
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto &pmesh_fault = static_cast<ParMesh &>(mesh_);
      for (int i = 0; i < fault_shared_faces_.Size(); ++i)
      {
         fault_face_idx_set.insert(
            pmesh_fault.GetSharedFace(fault_shared_faces_[i]));
      }
#endif
   }

   std::size_t n_int_built = 0;   // fully-local 2-sided central faces
   std::size_t n_shr_built = 0;   // shared (cross-rank) central faces

   // --- Interior (fully-local 2-sided) central faces -----------------------
   for (int mesh_face : central_flux_face_set_)
   {
      MFEM_ASSERT(fault_face_idx_set.count(mesh_face) == 0,
                  "BuildPerFaceCentralFluxMatrices_: central_flux_face_set_ "
                  "must exclude fault faces (got fault face " << mesh_face << ").");

      // self=Elem1, nbr=Elem2.  Returns false for shared / boundary faces
      // (shared central faces are handled by the sf loop below).
      const GodunovFlux *flux_e1 = nullptr;
      const GodunovFlux *flux_e2 = nullptr;
      real_t nor[3];
      if (!ResolveFaceFluxOperands_(mesh_face, flux_e1, flux_e2, nor))
      {
         continue;
      }
      auto &pair = per_face_central_flux_[mesh_face];
      BimaterialFlux::BuildPerFaceCentralMatricesGlobal(
         nor, *flux_e1, *flux_e2, pair[0], pair[1]);
      ++n_int_built;
   }

   // --- Shared (cross-rank) central faces, ParMesh only --------------------
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto &pmesh = static_cast<ParMesh &>(mesh_);
      const int n_faces = mesh_.GetNumFaces();   // (IMPL-7) used only here (parallel)
      const int n_shared = pmesh.GetNSharedFaces();
      for (int sf = 0; sf < n_shared; ++sf)
      {
         const int mesh_face_idx = pmesh.GetSharedFace(sf);
         // (IMPL-4) Fail loud on a corrupt sf->mesh_face map, matching the
         // Godunov build's Pass 2 — an out-of-range index is an invariant
         // violation, not an expected skip.
         MFEM_VERIFY(mesh_face_idx >= 0 && mesh_face_idx < n_faces,
                     "BuildPerFaceCentralFluxMatrices_: shared face " << sf
                     << " -> mesh face " << mesh_face_idx << " out of [0, "
                     << n_faces << ").");

         // Only central-set shared faces (skips non-central seams + faults,
         // which the base excludes from central_flux_face_set_).
         const bool is_central = central_flux_face_set_.count(mesh_face_idx) > 0;

         FaceElementTransformations *ftr =
            pmesh.GetSharedFaceTransformations(sf);
         // (P3-2) A null transform is a benign skip for a NON-central shared
         // face, but a central-set shared face MUST have a valid transform: a
         // silent `continue` here would leave a central-set face with NO entry
         // in per_face_central_flux_, tripping the exact-equality IMPL-8 assert
         // (per_face_central_flux_.size() == central_flux_face_set_.size()) at
         // SetMixedFluxMode time and aborting a VALID run.  Fail loud instead so
         // the invariant violation (a central face with no transform) is named
         // at its source rather than mis-attributed downstream.
         if (!ftr)
         {
            MFEM_VERIFY(!is_central,
                        "BuildPerFaceCentralFluxMatrices_: shared face " << sf
                        << " -> mesh face " << mesh_face_idx << " is in "
                        "central_flux_face_set_ but GetSharedFaceTransformations("
                        << sf << ") returned null — a central (mixed-flux) shared "
                        "face must have a valid face transform.  A silent skip "
                        "here would leave the central matrix unbuilt and trip the "
                        "IMPL-8 size-equality assert.");
            continue;
         }

         if (!is_central) { continue; }

         // Defensive fault skip (same guard as the Godunov build's Pass 2).
         const bool sf_fault =
            (sf < static_cast<int>(shared_face_bdr_attr_.size()))
            && (shared_face_bdr_attr_[sf] == bc_.fault_attr)
            && (bc_.fault_attr > 0);
         if (sf_fault) { continue; }

         // (PLAN BUG-6/BUG-10/BUG-15/BUG-16) Lateral-heterogeneity guard.  The
         // cross-rank neighbour material is a LOCAL-side stub
         // (ExchangeBiMaterialNeighbours_ stores local-as-neighbour, R-004);
         // for a non-Constant material on a partition seam the stub's A_nbr may
         // be wrong, and there is no communication-free way to detect lateral
         // variation (the field has no arbitrary-point evaluator the operator
         // can reach; an along-normal probe conflates depth with lateral
         // variation).  So a central shared face under a non-Constant material
         // requires the user to AFFIRM seam-continuity via
         // [material].seam_continuous=true.  This guard fires ONLY here, on a
         // central-set SHARED face, so mixed_flux=none (empty set) is unaffected.
         MFEM_VERIFY(material_->mode == MaterialField::Mode::Constant
                     || seam_continuous_,
                     "BimaterialWaveOperator::BuildPerFaceCentralFluxMatrices_: "
                     "central (mixed) flux on a fault-adjacent SHARED face "
                     "(shared face " << sf << " -> mesh face " << mesh_face_idx
                     << ") with a non-Constant [material] requires "
                     "[material].seam_continuous=true.  The cross-rank "
                     "neighbour-material exchange is a local-side stub (R-004): "
                     "it uses the LOCAL element's material as the neighbour's, "
                     "which is correct only for seam-continuous (e.g. depth-only) "
                     "materials.  Set seam_continuous=true to affirm the material "
                     "is seam-continuous, or wait for the cross-rank "
                     "MPI_Allgatherv exchange (deferred).");

         const int local_elem = ftr->Elem1No;
         MFEM_ASSERT(local_elem >= 0 && local_elem < ne_,
                     "BuildPerFaceCentralFluxMatrices_: shared face " << sf
                     << " Elem1No=" << local_elem << " outside [0, " << ne_
                     << ").");

         auto nbr_it = shared_face_neighbour_material_.find(sf);
         MFEM_VERIFY(nbr_it != shared_face_neighbour_material_.end(),
                     "BuildPerFaceCentralFluxMatrices_: shared face " << sf
                     << " missing entry in shared_face_neighbour_material_.");
         const auto &lmr_nbr = nbr_it->second;
         GodunovFlux flux_nbr(lmr_nbr[0], lmr_nbr[1], lmr_nbr[2]);

         real_t nor[3];
         CentroidUnitNormal_(ftr, nor);   // (IMPL-3) shared centroid-normal helper

         const GodunovFlux &flux_local = owned_flux_pool_->At(local_elem);

         // Side-0 pair (0.5*A_local, 0.5*A_nbr) — same self/nbr ordering as the
         // Godunov shared build, so the Phase-3 SharedInteriorFaceFlux_ deposit
         // 0.5*A_local.Q_self + 0.5*A_nbr.Q_nbr is consistent.
         auto &pair = per_face_central_flux_[mesh_face_idx];
         BimaterialFlux::BuildPerFaceCentralMatricesGlobal(
            nor, flux_local, flux_nbr, pair[0], pair[1]);
         ++n_shr_built;
      }
#endif
   }

   // Account the central-matrix bytes (mirror the Godunov precompute log + the
   // 4 GiB per-rank soft cap).  Each stored face holds ONE pair (2 matrices).
   const std::size_t per_face_bytes =
      2 * NUM_STATE * NUM_STATE * sizeof(real_t);
   const std::size_t total_bytes = (n_int_built + n_shr_built) * per_face_bytes;

   int print_rank = 0;
#ifdef MFEM_USE_MPI
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      auto &pmesh = static_cast<ParMesh &>(mesh_);
      MPI_Comm_rank(pmesh.GetComm(), &print_rank);
   }
#endif
   if (print_rank == 0)
   {
      mfem::out << "[wave_operator] BimaterialFlux CENTRAL precomputation:\n"
                << "  interior central faces = " << n_int_built << "\n"
                << "  shared central faces   = " << n_shr_built
                << "   (1 side each, Elem1 = local)\n"
                << "  total bytes            = " << total_bytes
                << " (per-rank)\n";
      constexpr std::size_t kSoftCapBytes = 4ULL * 1024 * 1024 * 1024;
      if (total_bytes > kSoftCapBytes)
      {
         mfem::out << "[wave_operator] WARNING: per-rank BimaterialFlux CENTRAL "
                   << "precomputation = " << total_bytes << " bytes (> "
                   << kSoftCapBytes << " soft cap).  Partition more "
                   << "aggressively or enable a runtime-recompute fallback.\n";
      }
   }
}

// ---------------------------------------------------------------------------
// Phase 13 hook overrides — material-/flux-dependent sites use per-element
// material (the per-face bi-material matrices / the flux pool) at every site.
// ---------------------------------------------------------------------------

// Interior NON-fault flux, fully-local (2-sided) face.  Each side applies its
// OWN A_self to the SAME per-face bi-material Riemann state; both consume
// (Q_self=Q_e1, Q_nbr=Q_e2) — NOT swapped.  `nor` is unused (the per-face
// matrices already bake in the centroid normal).
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::InteriorFaceFlux_(
   int mesh_face, const real_t *Q_self, const real_t *Q_nbr,
   const real_t * /*nor*/, real_t *F_h_e1, real_t *F_h_e2) const
{
   // (Phase 3, BUG-3) Mixed-flux CENTRAL branch on fault-adjacent faces: the
   // central flux is SINGLE-VALUED — compute F* once and copy it to side 2
   // (F_h_e2 = F_h_e1), mirroring the scalar reference.  The per-face pair is
   // (½A_e1, ½A_e2) (NO swap), and the runtime passes the UNSWAPPED
   // (Q_self=Q_e1, Q_nbr=Q_e2), so F* = ½A_e1·Q_e1 + ½A_e2·Q_e2 is identical on
   // both sides for ANY materials.  `.at()` aborts loudly on a missing central
   // matrix for a central-set face (never a silent zero).
   if (mf_on_ && central_flux_face_set_.count(mesh_face) > 0)
   {
      const auto &c = per_face_central_flux_.at(mesh_face);
      BimaterialFlux::ApplyPerFaceFlux(c[0], c[1], Q_self, Q_nbr, F_h_e1);
      for (int comp = 0; comp < NUM_STATE; ++comp) { F_h_e2[comp] = F_h_e1[comp]; }
      phaser_dispatch_count_ += 2;
      return;
   }

   const auto &mat_e1_local = per_face_bimaterial_flux_[mesh_face][0][0];
   const auto &mat_e1_nbr   = per_face_bimaterial_flux_[mesh_face][0][1];
   const auto &mat_e2_local = per_face_bimaterial_flux_[mesh_face][1][0];
   const auto &mat_e2_nbr   = per_face_bimaterial_flux_[mesh_face][1][1];
   BimaterialFlux::ApplyPerFaceFlux(mat_e1_local, mat_e1_nbr,
                                    Q_self, Q_nbr, F_h_e1);
   BimaterialFlux::ApplyPerFaceFlux(mat_e2_local, mat_e2_nbr,
                                    Q_self, Q_nbr, F_h_e2);
   phaser_dispatch_count_ += 2;
}

// Interior NON-fault flux, SHARED (cross-rank) face — local element is Elem1
// by MFEM convention; only side=0 of per_face_bimaterial_flux_ is populated.
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::SharedInteriorFaceFlux_(
   int mesh_face, const real_t *Q_self, const real_t *Q_nbr,
   const real_t * /*nor*/, real_t *F_h) const
{
   // (Phase 3, BUG-3) Mixed-flux CENTRAL branch (shared face, side 0 only): a
   // single single-valued F* = ½A_local·Q_self + ½A_nbr·Q_nbr.
   if (mf_on_ && central_flux_face_set_.count(mesh_face) > 0)
   {
      const auto &c = per_face_central_flux_.at(mesh_face);
      BimaterialFlux::ApplyPerFaceFlux(c[0], c[1], Q_self, Q_nbr, F_h);
      phaser_dispatch_count_ += 1;
      return;
   }

   const auto &mat_local = per_face_bimaterial_flux_[mesh_face][0][0];
   const auto &mat_nbr   = per_face_bimaterial_flux_[mesh_face][0][1];
   BimaterialFlux::ApplyPerFaceFlux(mat_local, mat_nbr, Q_self, Q_nbr, F_h);
   phaser_dispatch_count_ += 1;
}

// ADER CK-recursion element Jacobian apply via per-element star matrices.
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::ApplyElementJacobian_(
   int dir, const Vector &X, Vector &Y, real_t sign) const
{
   ApplyJacobianPerElementDOF_(dir, X, Y, sign);
}

// ---------------------------------------------------------------------------
// Phase 13 ComputeMaxDt override — per-element heterogeneous CFL walk.
// (The `if (!per_elem_h_.empty())` guard + scalar `h_min_/flux_.GetCp()`
// fallback of the single-class form are gone: on this class per_elem_h_ is
// always populated and `flux_` is the (1,1,1) poison, so the per-element walk
// + MPI_Allreduce(MIN) is the ONLY correct path.)
// ---------------------------------------------------------------------------
template <typename MeshType>
real_t BimaterialWaveOperator<MeshType>::ComputeMaxDt(real_t cfl) const
{
   MFEM_PERF_SCOPE("seas::BimaterialWaveOperator::ComputeMaxDt");

   // (Phase 4) Central+ADER guard — mirrors WaveOperator::ComputeMaxDt's guard
   // (same condition + message) so NEITHER operator can run central flux under
   // ADER.  The bi-material CENTRAL flux is non-dissipative (the average of the
   // two physical fluxes), so it is UNSTABLE under ADER's R(z)=1+z+z^2/2 (which
   // excludes the imaginary axis the central flux populates) and REQUIRES an RK
   // integrator (BUILD section 5.4/6.3; cfl_rk_aware_ is set true only by the
   // spatial driver's RK branch via SetCflRkAware).  Does NOT fire for
   // mixed_flux=none (the parity / TPV31-ADER byte-exact contract).
   MFEM_VERIFY(mixed_flux_mode_ == MixedFluxMode::None || cfl_rk_aware_,
               "BimaterialWaveOperator::ComputeMaxDt: mixed (central) flux on "
               "the matrix path is non-dissipative and UNSTABLE under ADER; it "
               "requires an RK integrator (set cfl_rk_aware_ via SetCflRkAware "
               "from the RK driver branch).  Use time_integrator=\"rk4\"/"
               "\"rk45\" with interior_flux=\"matrix\" + mixed_flux="
               "\"adjacent\"/\"all_continuous\", or mixed_flux=\"none\" "
               "under ADER.");

   // (Phase 4, BUG-2) Single source of truth: the SAME factor as the scalar
   // path (cfl_rk_aware_-gated 0.6/0.7, NOT the old divergent 0.9/0.4).  The
   // local divergent switch is deleted; both operators now route through
   // WaveOperator::MixedFluxCflFactor_().
   const real_t cfl_mixed_flux_factor = this->MixedFluxCflFactor_();

   MFEM_ASSERT(per_elem_h_.size() == static_cast<size_t>(ne_)
               && per_elem_lmr_.size() == static_cast<size_t>(ne_),
               "ComputeMaxDt: per_elem_h_ / per_elem_lmr_ sizes ("
               << per_elem_h_.size() << ", " << per_elem_lmr_.size()
               << ") must equal ne_ (" << ne_ << ").");
   real_t local_dt_min = std::numeric_limits<real_t>::infinity();
   for (int e = 0; e < ne_; ++e)
   {
      const real_t lam = per_elem_lmr_[e][0];
      const real_t mu  = per_elem_lmr_[e][1];
      const real_t rho = per_elem_lmr_[e][2];
      const real_t cp_e = std::sqrt((lam + 2.0 * mu) / rho);
      const real_t h_e  = per_elem_h_[e];
      const real_t dt_e = cfl_mixed_flux_factor * cfl * h_e / cp_e;
      local_dt_min = std::min(local_dt_min, dt_e);
   }
   real_t dt_global = local_dt_min;
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&local_dt_min, &dt_global, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN,
                    static_cast<ParMesh &>(mesh_).GetComm());
#endif
   }
   return dt_global;
}

// ---------------------------------------------------------------------------
// Phase 3 SetMixedFluxMode override — LIFTS the R-003 matrix x mixed-flux
// mutual exclusion (PLAN BUG-1).  The bi-material CENTRAL flux
// (½A_self·Q_self + ½A_nbr·Q_nbr, BuildPerFaceCentralMatricesGlobal) is now
// dispatched per-face on the central-set faces alongside the bi-material
// Godunov upwind elsewhere — exactly drdg3d's get_flux structure (central on
// fluxtype==1 faces, upwind otherwise).  Lifting R-003 is justified by that
// oracle and is the purpose of this feature (CLAUDE.md "never revert a prior
// fix without justification").
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::SetMixedFluxMode(MixedFluxMode m)
{
   // Base does: validate the mode, the rank-consistency Allgather/Allreduce
   // mode check, the `Adjacent => fault_attr>0` verify, the
   // `use_precomputed_face_fluxes_` mutual-exclusion (false here, defensive),
   // BuildCentralFluxFaceSet_ (the ONLY populator of central_flux_face_set_),
   // and sets mf_on_ / mixed_flux_mode_.  It MUST run BEFORE the central
   // precompute so the set is populated (IMPL-8 lifecycle).
   WaveOperator<MeshType>::SetMixedFluxMode(m);

   // Precompute the per-face central matrices for central_flux_face_set_
   // (clears first => re-entrant / idempotent on a mode change).
   BuildPerFaceCentralFluxMatrices_();

   // (IMPL-8 + Phase-2 acceptance) Lifecycle / no-silent-no-op guard: every
   // face in central_flux_face_set_ is either a 2-sided local interior face
   // (built by the interior arm) or a shared face (built by the sf arm), so the
   // map size MUST equal the set size.  A mismatch means the set was not
   // populated before the build (ordering bug) or a central face was silently
   // skipped.  Holds trivially for mixed_flux=none (both empty).
   MFEM_VERIFY(per_face_central_flux_.size() == central_flux_face_set_.size(),
               "BimaterialWaveOperator::SetMixedFluxMode: per_face_central_flux_ "
               "size (" << per_face_central_flux_.size() << ") != "
               "central_flux_face_set_ size (" << central_flux_face_set_.size()
               << ") after BuildPerFaceCentralFluxMatrices_ — central set not "
               "populated before the build, or a central face was skipped.");
}
