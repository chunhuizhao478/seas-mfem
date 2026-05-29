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

   auto compute_centroid_unit_normal = [](FaceElementTransformations *ftr,
                                          real_t out_n[3])
   {
      const Geometry::Type gtype = ftr->GetGeometryType();
      const IntegrationPoint &ip0 = Geometries.GetCenter(gtype);
      ftr->SetAllIntPoints(&ip0);
      Vector nor_vec(3);
      CalcOrtho(ftr->Face->Jacobian(), nor_vec);
      const real_t nor_len = nor_vec.Norml2();
      MFEM_VERIFY(nor_len > 0,
                  "BuildPerFaceBimaterialFluxMatrices_: face has zero "
                  "Jacobian normal at centroid.");
      nor_vec /= nor_len;
      out_n[0] = nor_vec(0);
      out_n[1] = nor_vec(1);
      out_n[2] = nor_vec(2);
   };

   std::size_t n_int_built = 0;   // fully-local 2-sided interior faces
   std::size_t n_shr_built = 0;   // shared (ParMesh) interior faces

   // --- Pass 1: fully-local interior faces ---------------------------------
   for (int f = 0; f < n_faces; ++f)
   {
      if (fault_face_idx_set.count(f) > 0) { continue; }

      FaceElementTransformations *ftr =
         mesh_.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;

      if (e2 < 0) { continue; }                  // 1-sided (boundary / seam)
      if (face_bdr_attr_[f] != 0) { continue; }  // non-fault boundary face

      real_t nor[3];
      compute_centroid_unit_normal(ftr, nor);

      const GodunovFlux &flux_e1 = owned_flux_pool_->At(e1);
      const GodunovFlux &flux_e2 = owned_flux_pool_->At(e2);

      build_face_matrices(nor, flux_e1, flux_e2,
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
         compute_centroid_unit_normal(ftr, nor);

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
   // Mixed flux is disabled on the matrix path (SetMixedFluxMode aborts on
   // non-None), so mixed_flux_mode_ is always None here and the factor is 1;
   // the switch is kept for parity with the scalar base implementation.
   real_t cfl_mixed_flux_factor = 1.0;
   switch (mixed_flux_mode_)
   {
      case MixedFluxMode::None:          cfl_mixed_flux_factor = 1.0; break;
      case MixedFluxMode::Adjacent:      cfl_mixed_flux_factor = 0.9; break;
      case MixedFluxMode::AllContinuous: cfl_mixed_flux_factor = 0.4; break;
      default:
         MFEM_ABORT("BimaterialWaveOperator::ComputeMaxDt: unknown "
                    "MixedFluxMode " << static_cast<int>(mixed_flux_mode_)
                    << " (R-1600 fall-through guard).");
   }

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
// Phase 13 SetMixedFluxMode override — mixed flux is scalar-only (R-003,
// now structural).  The bi-material Riemann solve already replaces the
// interior-face flux, so any non-None mode would be silently ignored.
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::SetMixedFluxMode(MixedFluxMode m)
{
   MFEM_VERIFY(m == MixedFluxMode::None,
               "BimaterialWaveOperator::SetMixedFluxMode: mixed flux is "
               "incompatible with the heterogeneous (bimaterial) "
               "interior_flux=\"matrix\" path; the bimaterial Riemann solve "
               "already replaces the interior-face flux.  Use "
               "interior_flux=\"scalar\" (WaveOperator) for mixed flux.");
   // m == None is the constructed default; nothing to do.
}
