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
   const BoundaryConfig &bc,
   const std::vector<int> *lts_cluster_id)
   : Base(mesh, order, real_t(1.0), real_t(1.0), real_t(1.0), bc, lts_cluster_id)
{
   // (Cross-rank Phase 1) All three material representations are accepted:
   // every material read routes through MaterialAtLocal_/MaterialAtNbr_, whose
   // GridFunction branch reads the source GFs natively (GetValue locally,
   // shape×FaceNbrData across a seam) — so the per-element flux pool's centroid
   // evaluation is well defined in every mode.  The enum is a fixed 3-value set;
   // assert it is one of them (guards against an uninitialised/garbage mode).
   MFEM_VERIFY(material.mode == MaterialField::Mode::Constant
               || material.mode == MaterialField::Mode::Coefficient
               || material.mode == MaterialField::Mode::GridFunction,
               "BimaterialWaveOperator(MaterialField): unrecognised material.mode ("
               << static_cast<int>(material.mode) << ").");

   material_ = &material;

   // Face-neighbour exchange (geometry + GridFunction material GFs) BEFORE the
   // pool / neighbour builds, so MaterialAtNbr_ peer reads are valid (R-6 ctor
   // order; R-5 collective symmetry — runs on ALL ranks).
   SetupMaterialFaceNbrExchange_();

   BuildGodunovFluxPool_();
   ExchangeBiMaterialNeighbours_();
   BuildPerFaceBimaterialFluxMatrices_();
}

// ---------------------------------------------------------------------------
// (Cross-rank Phase 1) One-time face-neighbour exchange setup.  Parallel only;
// collective on ALL ranks (no per-rank-conditional exchange — R-5/R-209).
//   - Geometry: pmesh.ExchangeFaceNbrData() (idempotent — the base WaveOperator
//     ctor already called it; ParMesh guards on have_face_nbr_data).  Required
//     for GetSharedFaceTransformations(sf)->Elem2 (the peer ElementTransformation
//     MaterialAtNbr_ reads in Coefficient mode).
//   - Material MPI: ONLY Mode::GridFunction needs it — the peer read interpolates
//     the source GFs' face-neighbour ghost layer (FaceNbrData()).  Constant /
//     Coefficient are globally evaluable, so they exchange NOTHING (asserted by
//     the accessor_coeff_no_material_mpi gate via material_gf_exchange_count_==0).
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::SetupMaterialFaceNbrExchange_()
{
   material_gf_exchange_count_ = 0;
   if constexpr (!IsParallelMesh<MeshType>::value)
   {
      return;   // serial: no peers; MaterialAtNbr_ is never called
   }
#ifdef MFEM_USE_MPI
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      auto &pmesh = static_cast<ParMesh &>(mesh_);
      pmesh.ExchangeFaceNbrData();   // idempotent; ensures peer transforms exist

      if (material_ != nullptr
          && material_->mode == MaterialField::Mode::GridFunction)
      {
         MFEM_VERIFY(material_->lambda_gf && material_->mu_gf && material_->rho_gf,
                     "SetupMaterialFaceNbrExchange_: GridFunction mode requires "
                     "all three source GFs to be non-null.");
         // R-402 / P1-003: material GFs MUST be scalar (vdim==1) — GetValue reads
         // a single scalar component (vdim>1 would read a per-component subvector).
         // The order MUST equal the state/wave order p (the plan's "scalar
         // L2(order p)" contract; a different-order GF samples the material at a
         // mismatched polynomial resolution).  Checked per GF (each reads through
         // its OWN ParFESpace — P1-002/P1-004: MaterialAtNbr_ uses the library
         // GetValue(Elem2No, ip) per GF, so the three GFs need NOT share a space).
         const int state_order = this->GetFESpace().GetMaxElementOrder();
         auto check_gf = [&](const mfem::ParGridFunction *gf, const char *name)
         {
            mfem::ParFiniteElementSpace *sp = gf->ParFESpace();
            MFEM_VERIFY(sp->GetVDim() == 1,
                        "SetupMaterialFaceNbrExchange_: material GridFunction '"
                        << name << "' must be scalar (vdim==1); got vdim="
                        << sp->GetVDim() << ".");
            MFEM_VERIFY(sp->GetMaxElementOrder() == state_order,
                        "SetupMaterialFaceNbrExchange_: material GridFunction '"
                        << name << "' order " << sp->GetMaxElementOrder()
                        << " != state order " << state_order << ".");
         };
         check_gf(material_->lambda_gf.get(), "lambda");
         check_gf(material_->mu_gf.get(),     "mu");
         check_gf(material_->rho_gf.get(),    "rho");
         // shared_ptr operator-> yields a non-const ParGridFunction* even through
         // the const MaterialField*, so these non-const exchanges are well-formed.
         material_->lambda_gf->ExchangeFaceNbrData(); ++material_gf_exchange_count_;
         material_->mu_gf->ExchangeFaceNbrData();     ++material_gf_exchange_count_;
         material_->rho_gf->ExchangeFaceNbrData();    ++material_gf_exchange_count_;
      }
   }
#endif
}

// ---------------------------------------------------------------------------
// (Cross-rank Phase 1) Uniform material accessors — local and peer.  The ONLY
// branch is the material representation; the consumers never ask "which problem?".
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::MaterialAtLocal_(
   int elem, ElementTransformation &T, const IntegrationPoint &ip,
   real_t &lam, real_t &mu, real_t &rho) const
{
   MFEM_ASSERT(material_ != nullptr, "MaterialAtLocal_: material_ not set.");
   if (material_->mode == MaterialField::Mode::GridFunction)
   {
      MFEM_ASSERT(material_->lambda_gf && material_->mu_gf && material_->rho_gf,
                  "MaterialAtLocal_: GridFunction mode requires three GFs.");
      lam = material_->lambda_gf->GetValue(elem, ip);
      mu  = material_->mu_gf->GetValue(elem, ip);
      rho = material_->rho_gf->GetValue(elem, ip);
   }
   else
   {
      // Constant + Coefficient: byte-exact with the legacy read — this IS
      // MaterialField::EvalAt at the CALLER-supplied element transform.  (P1-001)
      // No internal GetElementTransformation: the accessor never touches the mesh
      // shared scratch, so it composes with a live FaceElementTransformations.
      // (P2-005) Set the int point before EvalAt: byte-exact for the explicit-ip
      // FunctionCoefficients used today, but required for any coefficient that
      // reads T.GetIntPoint() (matches MaxCpInElement; closes the silent
      // wrong-material footgun for an analytic-CVM coefficient).
      T.SetIntPoint(&ip);
      material_->EvalAt(elem, T, ip, lam, mu, rho);
   }
   MFEM_ASSERT(rho > 0.0 && (lam + 2.0 * mu) > 0.0,
               "MaterialAtLocal_: non-physical material at element " << elem
               << " (rho=" << rho << ", lambda+2mu=" << (lam + 2.0 * mu) << ").");
}

template <typename MeshType>
void BimaterialWaveOperator<MeshType>::MaterialAtNbr_(
   FaceElementTransformations *ftr, const IntegrationPoint &ip_peer,
   real_t &lam, real_t &mu, real_t &rho) const
{
   MFEM_ASSERT(material_ != nullptr, "MaterialAtNbr_: material_ not set.");
   MFEM_ASSERT(ftr != nullptr && ftr->Elem2 != nullptr,
               "MaterialAtNbr_: null face transform or peer (Elem2) transform.");
   if (material_->mode == MaterialField::Mode::GridFunction)
   {
#ifdef MFEM_USE_MPI
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         MFEM_ASSERT(material_->lambda_gf && material_->mu_gf && material_->rho_gf,
                     "MaterialAtNbr_: GridFunction mode requires three GFs.");
         // Shared-face peer index: Elem2No = ne_ + nbr_idx (ParMesh convention).
         const int nbr_idx = ftr->Elem2No - ne_;
         MFEM_VERIFY(nbr_idx >= 0,
                     "MaterialAtNbr_: Elem2No=" << ftr->Elem2No
                     << " (- ne_=" << ne_ << " => " << nbr_idx << " < 0) is not a "
                     "face-neighbour; MaterialAtNbr_ is for SHARED faces only.");
         // (P1-002/P1-003/P1-004) Read each GF through its OWN space via the MFEM
         // library face-neighbour path: ParGridFunction::GetValue(i, ip) treats
         // i >= GetNE() (here Elem2No) as a face-nbr element — it does the
         // GetFaceNbrElementVDofs + GetFaceNbrFE + face_nbr_data interpolation,
         // applies the proper DofTransformation (InvTransformPrimal) and map type.
         // The three GFs therefore need NOT share one ParFiniteElementSpace.
         lam = material_->lambda_gf->GetValue(ftr->Elem2No, ip_peer);
         mu  = material_->mu_gf->GetValue(ftr->Elem2No, ip_peer);
         rho = material_->rho_gf->GetValue(ftr->Elem2No, ip_peer);
      }
      else
#endif
      {
         MFEM_ABORT("MaterialAtNbr_: GridFunction peer read requires a parallel "
                    "mesh build (MFEM_USE_MPI + ParMesh).");
      }
   }
   else
   {
      // Constant + Coefficient: globally evaluable — the SAME Eval at the peer's
      // face-neighbour element transform (ftr->Elem2).  No material MPI.  EvalAt
      // ignores the element index argument in these modes.  (P2-005) Set the peer
      // int point before EvalAt (the peer transform's stored ip may be stale from
      // the last GetSharedFaceTransformations/SetAllIntPoints): byte-exact for the
      // explicit-ip FunctionCoefficients, required for a GetIntPoint-reading one.
      ftr->Elem2->SetIntPoint(&ip_peer);
      material_->EvalAt(ftr->Elem2No, *ftr->Elem2, ip_peer, lam, mu, rho);
   }
   MFEM_ASSERT(rho > 0.0 && (lam + 2.0 * mu) > 0.0,
               "MaterialAtNbr_: non-physical peer material (rho=" << rho
               << ", lambda+2mu=" << (lam + 2.0 * mu) << ").");
}

// ---------------------------------------------------------------------------
// Phase H.1 helper — fill per_elem_lmr_ / per_elem_h_, then Build()
// the owned flux pool.
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::BuildGodunovFluxPool_()
{
   per_elem_lmr_.assign(static_cast<size_t>(ne_), std::array<real_t, 3>{0, 0, 0});
   per_elem_h_.assign(static_cast<size_t>(ne_), real_t(0));

   for (int e = 0; e < ne_; ++e)
   {
      ElementTransformation *T = mesh_.GetElementTransformation(e);
      const Geometry::Type   gtype = mesh_.GetElementBaseGeometry(e);
      const IntegrationPoint &ip   = Geometries.GetCenter(gtype);
      real_t lam, mu, rho;
      // (Cross-rank Phase 1) Route through MaterialAtLocal_ so the pool is well
      // defined in ALL three modes (GridFunction included).  Byte-exact with the
      // legacy `material.EvalAt(e, *T, ip, ...)` for Constant/Coefficient: the
      // accessor's local branch IS that same EvalAt at the same element transform
      // (the loop's own `T`, not a held face transform) and centroid IP.  This
      // loop holds no FaceElementTransformations, so reusing the mesh scratch via
      // GetElementTransformation here is safe (P1-001 concerns held-ftr callers).
      MaterialAtLocal_(e, *T, ip, lam, mu, rho);
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
// (Cross-rank Phase 2) Populate shared_face_neighbour_material_ with each shared
// face's TRUE peer (face-neighbour) material.
//
// Replaces the former LOCAL-SIDE STUB (R-004), which stored the local element's
// own material as the neighbour's — correct only for seam-continuous (depth-only)
// materials, and the reason the central build needed the seam_continuous
// affirmation.  We now read the genuine peer material through the uniform
// accessor MaterialAtNbr_ at the PEER element's reference centroid:
//   - Mode::Constant: the constant (byte-exact with the stub).
//   - Mode::Coefficient: coeff->Eval at ftr->Elem2's centroid — globally
//     evaluable, NO material MPI.  Byte-exact with the stub wherever the material
//     is seam-continuous (peer centroid material == local centroid material);
//     CORRECT (was wrong) wherever it varies laterally across the seam.
//   - Mode::GridFunction: the source GFs' face-neighbour ghost layer (exchanged
//     once at construction), via ParGridFunction::GetValue(Elem2No, centroid).
// On serial Mesh this is a no-op (no shared faces).
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
         MFEM_ASSERT(ftr->Elem2 != nullptr,
                     "ExchangeBiMaterialNeighbours_: shared face " << sf
                     << " has no peer (Elem2) transform; ExchangeFaceNbrData "
                     "(SetupMaterialFaceNbrExchange_) must run first.");
         // TRUE peer material at the peer element's reference centroid.
         const IntegrationPoint &cpeer =
            Geometries.GetCenter(ftr->Elem2->GetGeometryType());
         real_t lam, mu, rho;
         MaterialAtNbr_(ftr, cpeer, lam, mu, rho);
         shared_face_neighbour_material_[sf] = {lam, mu, rho};
      }

      // (P2-002) Informational rank-0 banner (NOT a warning — this is the correct
      // behavior).  A non-Constant material's seam neighbour is now the TRUE peer
      // (was a local-side stub); at a depth-/lateral-varying seam this re-baselines
      // the bulk shared Riemann vs the pre-2026-06 binary, so a parallel run may
      // differ at seams from an old gold by an O(grad·h) correction.
      if (material_ != nullptr
          && material_->mode != MaterialField::Mode::Constant
          && n_shared > 0)
      {
         int rank = 0;
         MPI_Comm_rank(pmesh.GetComm(), &rank);
         if (rank == 0)
         {
            mfem::out << "[wave_operator] cross-rank seam material = TRUE peer "
                         "(MaterialAtNbr_); bulk shared Riemann re-baselined at "
                         "material-varying seams vs the pre-2026-06 local-side "
                         "stub.\n";
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

   // (Unified bi-material plan, Part A) Central-flux CONTRAST GUARD.  When
   // mixed_flux_contrast_tol_ >= 0, a fault-adjacent corridor face whose two
   // elements span a STRONG impedance contrast is DROPPED from
   // central_flux_face_set_ so it dispatches the (dissipative) bi-material upwind
   // (per_face_bimaterial_flux_, built for every non-fault interior face in the
   // ctor) instead of the non-dissipative central flux.  tol < 0 (default) => the
   // whole guard is inert => byte-exact (no reclassification, no histogram, no new
   // MPI call).  The single `guard_on` flag is rank-uniform (tol is a config/CLI
   // scalar set identically on all ranks before SetMixedFluxMode), so the
   // MPI_Reduce histogram below is collective-safe (R-006).
   const real_t contrast_tol = mixed_flux_contrast_tol_;
   const bool   guard_on     = (contrast_tol >= 0.0);

   // (P2-003) The seam-continuity ABORT is gone (Phase 2).  When the guard is ALSO
   // off (default tol<0) and the material is non-Constant, a strong impedance
   // contrast on a central corridor face — local OR shared — builds a
   // non-dissipative central flux silently (no abort, no reclassify).  Emit a
   // one-time rank-0 WARNING so the absent safety net is visible until the contrast
   // guard is enabled (mixed_flux_contrast_tol >= 0).  (Consistent with the local
   // 2-sided central path, which has always built central across contrasts when the
   // guard is off; the warning now covers shared faces too.)
   if (!guard_on && material_ != nullptr
       && material_->mode != MaterialField::Mode::Constant
       && !central_flux_face_set_.empty())
   {
      int rank = 0;
#ifdef MFEM_USE_MPI
      if constexpr (IsParallelMesh<MeshType>::value)
      { MPI_Comm_rank(static_cast<ParMesh &>(mesh_).GetComm(), &rank); }
#endif
      if (rank == 0)
      {
         mfem::out << "[wave_operator] WARNING: mixed (central) flux on a "
                      "non-Constant material with the contrast guard DISABLED "
                      "(mixed_flux_contrast_tol < 0): a strong impedance contrast "
                      "on a corridor face builds a NON-dissipative central flux.  "
                      "Set mixed_flux_contrast_tol >= 0 to reclassify "
                      "strong-contrast faces to upwind.\n";
      }
   }

   std::vector<int> reclassified;       // mesh faces moved central -> upwind
   // Histogram of the contrast metric over the PRE-filter LOCAL (2-sided) corridor
   // faces only: those are disjoint across ranks, so MPI_Reduce(SUM) has NO double-
   // count (R-005).  Shared faces are excluded from the HISTOGRAM to avoid the MPI
   // double-count across the seam (both ranks see the same shared face); their
   // neighbour material is now the TRUE peer (Cross-rank Phase 2), so a shared face
   // CAN reclassify — the Phase-3 tally counts those separately (see the shared loop).
   long long hist[5] = {0, 0, 0, 0, 0};   // bins: [0,1) [1,5) [5,10) [10,20) [>=20] %
   long long n_local_hist   = 0;          // local corridor faces histogrammed
   long long n_reclass_local = 0;         // LOCAL faces reclassified (disjoint across
                                          // ranks => no MPI double-count; R-002)
   long long n_reclass_shared = 0;        // (Phase 3) SHARED faces reclassified on
                                          // THIS rank (a seam face is seen by both
                                          // adjacent ranks => the MPI_SUM below
                                          // counts it once per adjacent rank; R-005)

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
      if (guard_on)
      {
         const real_t contrast =
            BimaterialFlux::ContrastValue(*flux_e1, *flux_e2);
         if (contrast >= 0.0)            // histogram valid (non-acoustic) faces
         {
            const real_t pct = 100.0 * contrast;
            const int bin = (pct < 1.0)  ? 0 : (pct < 5.0)  ? 1
                          : (pct < 10.0) ? 2 : (pct < 20.0) ? 3 : 4;
            ++hist[bin];
            ++n_local_hist;
         }
         if (BimaterialFlux::IsStrongContrast(*flux_e1, *flux_e2, contrast_tol))
         {
            // Drop this face from the central corridor -> it dispatches the
            // bi-material upwind (per_face_bimaterial_flux_).  Defer the erase
            // (we are iterating central_flux_face_set_); skip the central build.
            reclassified.push_back(mesh_face);
            ++n_reclass_local;
            continue;
         }
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

         // (Cross-rank Phase 2) The former seam-continuity AFFIRMATION guard
         // (MFEM_VERIFY mode==Constant || seam_continuous_) is REMOVED: the
         // cross-rank neighbour material is now the TRUE peer material
         // (ExchangeBiMaterialNeighbours_ reads it via MaterialAtNbr_), so a
         // central (mixed) flux on a fault-adjacent SHARED face is correct for
         // ANY material (Constant, depth-only, lateral/CVM) without the user
         // affirming seam-continuity.  Strong lateral contrast on such a face is
         // handled below by the contrast guard (reclassify to upwind).
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

         if (guard_on &&
             BimaterialFlux::IsStrongContrast(flux_local, flux_nbr, contrast_tol))
         {
            // Strong contrast on a SHARED corridor face -> dispatch upwind.
            // (Cross-rank Phase 2) flux_nbr is now the TRUE peer material
            // (ExchangeBiMaterialNeighbours_ via MaterialAtNbr_), so for a real
            // lateral contrast across the seam this fires correctly; it stays
            // inert (flux_nbr == flux_local) for seam-continuous materials.
            // Shared faces are NOT histogrammed (avoids the MPI double-count,
            // R-005); the Phase-3 tally counts the per-rank reclassification.
            reclassified.push_back(mesh_face_idx);
            ++n_reclass_shared;
            continue;
         }

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

   // (Part A) Apply the contrast-guard reclassification: drop the strong-contrast
   // faces from central_flux_face_set_ so the dispatch (InteriorFaceFlux_ /
   // SharedInteriorFaceFlux_, which test central_flux_face_set_.count) routes them
   // to the existing bi-material upwind.  R-002: assert the upwind fallback exists
   // for each reclassified face BEFORE erasing (the ctor builds
   // per_face_bimaterial_flux_ for every non-fault interior face; a missing one
   // would be a silent garbage dispatch).  This keeps the IMPL-8 size invariant
   // (per_face_central_flux_.size() == central_flux_face_set_.size()) at
   // SetMixedFluxMode: reclassified faces are in NEITHER container.
   for (int f : reclassified)
   {
      MFEM_VERIFY(f >= 0
                  && f < static_cast<int>(per_face_bimaterial_flux_.size())
                  && per_face_bimaterial_flux_[f][0][0].Height() == NUM_STATE
                  && per_face_bimaterial_flux_[f][0][0].Width()  == NUM_STATE,
                  "BuildPerFaceCentralFluxMatrices_: contrast-guard reclassified "
                  "face " << f << " to upwind, but its bi-material upwind matrices "
                  "(per_face_bimaterial_flux_) are not built — no fallback flux to "
                  "dispatch.  The ctor must build the upwind matrices for every "
                  "non-fault interior face.");
      central_flux_face_set_.erase(f);
   }

   // (Phase 3) Publish this rank's shared reclassify count (0 if the guard is off,
   // since n_reclass_shared is only incremented in the guarded shared branch).
   n_reclass_shared_ = static_cast<std::size_t>(n_reclass_shared);

   // (Part A) Rank-0 contrast histogram + reclassification count.  ENTIRELY gated
   // on guard_on (rank-uniform), so tol<0 adds no statements and no collective
   // call (byte-exact); when on, the MPI_Reduce is called by all ranks (R-006).
   if (guard_on)
   {
      long long hist_g[5];
      for (int b = 0; b < 5; ++b) { hist_g[b] = hist[b]; }
      long long n_local_g  = n_local_hist;
      // R-002: LOCAL reclassified are disjoint across ranks => MPI_SUM is exact.
      // (Phase 3) SHARED reclassified are counted per-rank: a seam face is seen by
      // both adjacent ranks, so MPI_SUM(n_reclass_shared) counts each such face
      // ONCE PER ADJACENT RANK (2 for a 2-rank seam) — reported separately and
      // labelled, NOT folded into the exact local count (R-005).
      long long n_reclass  = n_reclass_local;
      long long n_reclass_g = n_reclass;
      long long n_reclass_shared_g = n_reclass_shared;
      int hrank = 0;
#ifdef MFEM_USE_MPI
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         auto &pmesh = static_cast<ParMesh &>(mesh_);
         MPI_Comm comm = pmesh.GetComm();
         MPI_Reduce(hist, hist_g, 5, MPI_LONG_LONG, MPI_SUM, 0, comm);
         MPI_Reduce(&n_local_hist, &n_local_g, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
         MPI_Reduce(&n_reclass, &n_reclass_g, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
         MPI_Reduce(&n_reclass_shared, &n_reclass_shared_g, 1, MPI_LONG_LONG,
                    MPI_SUM, 0, comm);
         MPI_Comm_rank(comm, &hrank);
      }
#endif
      // R-003: the MPI_Reduce above is UNCONDITIONAL under guard_on (collective-safe,
      // rank-uniform); only the rank-0 PRINT is gated on a nonzero global total, so a
      // matrix + mixed_flux=none + tol>=0 run does not emit a noise "0/0" line.
      if (hrank == 0 && (n_local_g > 0 || n_reclass_g > 0 || n_reclass_shared_g > 0))
      {
         mfem::out << "[wave_operator] central-flux CONTRAST GUARD (tol="
                   << contrast_tol << "):\n"
                   << "  local corridor faces histogrammed = " << n_local_g << "\n"
                   << "  contrast bins %% [0,1) [1,5) [5,10) [10,20) [>=20] = "
                   << hist_g[0] << " " << hist_g[1] << " " << hist_g[2] << " "
                   << hist_g[3] << " " << hist_g[4] << "\n"
                   << "  reclassified central -> upwind: local = " << n_reclass_g
                   << ", shared = " << n_reclass_shared_g
                   << " (shared = per-rank sum; a 2-rank seam face counts twice)\n";
      }
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

// (LTS Phase 2) Element-subset per-element-star Jacobian apply.  Same
// per-element FluxForElem_(e).GetReferenceStarMatrix(dir) and the same
// (c,cp,i) accumulation order as ApplyJacobianPerElementDOF_, restricted to the
// element list, so over the full list it is bit-identical to the whole-vector
// override above (the bimaterial single-cluster == GTS gate).
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::ApplyElementJacobianElems_(
   int dir, const Vector &X, Vector &Y, real_t sign, const int *elems, int n) const
{
   MFEM_ASSERT(dir >= 0 && dir < 3, "ApplyElementJacobianElems_: bad dir");
   MFEM_ASSERT(X.Size() == NUM_STATE * this->ndof_total_,
               "ApplyElementJacobianElems_: X size mismatch");
   MFEM_ASSERT(Y.Size() == NUM_STATE * this->ndof_total_,
               "ApplyElementJacobianElems_: Y size mismatch");
   const real_t *Xd = X.GetData();
   real_t *Yd = Y.GetData();
   for (int ei = 0; ei < n; ++ei)
   {
      const int e = elems[ei];
      const DenseMatrix &A = FluxForElem_(e).GetReferenceStarMatrix(dir);
      const int base = e * this->ndof_per_el_;
      for (int c = 0; c < NUM_STATE; ++c)
      {
         real_t *Yc = Yd + c * this->ndof_total_ + base;
         for (int cp = 0; cp < NUM_STATE; ++cp)
         {
            const real_t a = A(c, cp);
            if (a == 0.0) { continue; }
            const real_t w = sign * a;
            const real_t *Xcp = Xd + cp * this->ndof_total_ + base;
            for (int i = 0; i < this->ndof_per_el_; ++i) { Yc[i] += w * Xcp[i]; }
         }
      }
   }
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

// ---------------------------------------------------------------------------
// (Part B / B1) Per-side fault impedance assignment via the EPS-OFFSET rule.
//
// For each LOCAL 2-sided interior fault face, evaluate the material just INSIDE
// each side (the fault QP physical point displaced a small eps PERPENDICULAR to
// the fault into that element) and write Zp_plus/Zp_minus, Zs_plus/Zs_minus,
// eta_p/eta_s into the matching DOFData slot (fault_face_dof_offset_[f] + q).
//
// Why perpendicular (not toward the centroid, not at the face point):
//   - a fault-symmetric material (depth-only, TPV31): +/- eps along the fault
//     normal leaves the depth unchanged => both sides evaluate to the SAME
//     material => Zp_plus == Zp_minus (byte-exact; no spurious contrast on an
//     asymmetric mesh, the R-001 trap that centroid evaluation falls into);
//   - an across-fault material (halfspace, TPV6): +eps -> near side, -eps -> far
//     side => the correct per-side contrast (evaluating exactly at the face point
//     is ambiguous for a sign(.) coefficient).
// (Cross-rank Phase 4) SHARED (cross-rank) fault faces are now handled too: the
// LOCAL side is read with MaterialAtLocal_ and the PEER side across the seam with
// MaterialAtNbr_, BOTH at the eps-offset fault QP (in2 = -in1 ⇒ the SAME fault QP,
// opposite normal ⇒ exact symmetry on ANY mesh for a depth profile; correct
// per-side contrast for a halfspace; GridFunction supported).  Replaces the former
// "leave shared faces at the driver-seeded single-material values" stub (R-101).
// ---------------------------------------------------------------------------
template <typename MeshType>
void BimaterialWaveOperator<MeshType>::AssignFaultSidePerMaterialImpedances(
   std::vector<DOFData> &dof_data) const
{
   if (material_ == nullptr) { return; }
   // (Cross-rank Phase 4) All three modes are handled: the per-side reads route
   // through MaterialAtLocal_/MaterialAtNbr_, whose GridFunction branch uses
   // GetValue (no abort).  The former GridFunction early-return is removed.
   // Protected base members (the public getters are SEAS_TEST_INTERNAL-only):
   // interior_fault_elem1_on_plus_ is computed at ctor time; fault_face_dof_offset_
   // is populated by SetFaultDOFData (both done before the driver calls this).
   const std::vector<bool> &elem1_on_plus = this->interior_fault_elem1_on_plus_;
   const std::map<int, int> &dof_off_map = this->fault_face_dof_offset_;
   const int qdeg = this->FaultFaceQuadDegree();
   const int nbf  = this->GetNbfPerFace();
   const real_t eps_ref = 1.0e-3;   // small reference-space perpendicular offset

   // Perturb a face-QP reference IP a small eps along the (signed) physical fault
   // normal mapped into the element's reference frame (J^{-1}); affine (straight)
   // elements transform the offset point exactly, so no in-element clamp is needed.
   auto offset_ip = [&](ElementTransformation &T, const IntegrationPoint &ip,
                        const real_t nin[3]) -> IntegrationPoint
   {
      T.SetIntPoint(&ip);
      DenseMatrixInverse Jinv(T.Jacobian());
      Vector np(3), nr(T.GetDimension());
      np(0) = nin[0]; np(1) = nin[1]; np(2) = nin[2];
      Jinv.Mult(np, nr);
      const real_t rl = nr.Norml2();
      if (rl > 0.0) { nr /= rl; }
      IntegrationPoint out = ip;
      out.x = ip.x + eps_ref * nr(0);
      if (nr.Size() > 1) { out.y = ip.y + eps_ref * nr(1); }
      if (nr.Size() > 2) { out.z = ip.z + eps_ref * nr(2); }
      return out;
   };

   for (int fi = 0; fi < fault_interior_faces_.Size(); ++fi)
   {
      const int f = fault_interior_faces_[fi];
      FaceElementTransformations *ftr =
         const_cast<MeshType &>(mesh_).GetInteriorFaceTransformations(f);
      if (!ftr || ftr->Elem2No < 0) { continue; }   // need a 2-sided local face
      const int fb_idx = this->LookupInteriorFaultBasisIndex(f);
      if (fb_idx < 0 || fb_idx >= static_cast<int>(elem1_on_plus.size()))
      { continue; }
      const bool e1_plus = elem1_on_plus[fb_idx];
      auto off_it = dof_off_map.find(f);
      if (off_it == dof_off_map.end()) { continue; }
      const int dof_off = off_it->second;

      // Element centroids (physical) — used to orient the into-element direction
      // robustly (independent of CalcOrtho's face-normal sign convention).
      Vector c1(3), c2(3);
      ElementTransformation &T1c = *ftr->Elem1;
      ElementTransformation &T2c = *ftr->Elem2;
      T1c.Transform(Geometries.GetCenter(T1c.GetGeometryType()), c1);
      T2c.Transform(Geometries.GetCenter(T2c.GetGeometryType()), c2);

      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), qdeg);
      MFEM_VERIFY(ir.GetNPoints() == nbf,
                  "AssignFaultSidePerMaterialImpedances: fault face " << f
                  << " has " << ir.GetNPoints() << " QPs, expected nbf_per_face="
                  << nbf << " — integration rule does not match the DOFData "
                  "layout (would mis-index the per-side impedances).");

      for (int q = 0; q < nbf; ++q)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector xf(3);
         ftr->Face->Transform(ip, xf);
         Vector nvec(3);
         CalcOrtho(ftr->Face->Jacobian(), nvec);
         const real_t nlen = nvec.Norml2();
         if (nlen > 0.0) { nvec /= nlen; }

         // Into-Elem1 = sign((c1 - xf).n) * n ; Into-Elem2 = the opposite.
         real_t d1 = 0.0;
         for (int d = 0; d < 3; ++d) { d1 += (c1(d) - xf(d)) * nvec(d); }
         const real_t s1 = (d1 >= 0.0) ? 1.0 : -1.0;
         const real_t in1[3] = { s1 * nvec(0), s1 * nvec(1), s1 * nvec(2) };
         const real_t in2[3] = { -in1[0], -in1[1], -in1[2] };

         const IntegrationPoint ip1 =
            offset_ip(*ftr->Elem1, ftr->GetElement1IntPoint(), in1);
         const IntegrationPoint ip2 =
            offset_ip(*ftr->Elem2, ftr->GetElement2IntPoint(), in2);

         real_t lam1, mu1, rho1, lam2, mu2, rho2;
         // (Cross-rank Phase 4) Both sides are LOCAL for an interior fault face;
         // route through MaterialAtLocal_ (byte-exact with EvalAt for
         // Constant/Coefficient; GridFunction now supported via GetValue).
         MaterialAtLocal_(ftr->Elem1No, *ftr->Elem1, ip1, lam1, mu1, rho1);
         MaterialAtLocal_(ftr->Elem2No, *ftr->Elem2, ip2, lam2, mu2, rho2);
         MFEM_VERIFY(mu1 > 0.0 && rho1 > 0.0 && mu2 > 0.0 && rho2 > 0.0,
                     "AssignFaultSidePerMaterialImpedances: non-physical "
                     "material at fault face " << f << " QP " << q);

         const real_t Zp1 = std::sqrt((lam1 + 2.0 * mu1) * rho1);
         const real_t Zs1 = std::sqrt(mu1 * rho1);
         const real_t Zp2 = std::sqrt((lam2 + 2.0 * mu2) * rho2);
         const real_t Zs2 = std::sqrt(mu2 * rho2);

         const int idx = dof_off + q;
         if (idx < 0 || idx >= static_cast<int>(dof_data.size())) { continue; }
         DOFData &d = dof_data[idx];
         if (e1_plus)
         { d.Zp_plus = Zp1; d.Zs_plus = Zs1; d.Zp_minus = Zp2; d.Zs_minus = Zs2; }
         else
         { d.Zp_plus = Zp2; d.Zs_plus = Zs2; d.Zp_minus = Zp1; d.Zs_minus = Zs1; }
         d.eta_p = d.Zp_plus * d.Zp_minus / (d.Zp_plus + d.Zp_minus);
         d.eta_s = d.Zs_plus * d.Zs_minus / (d.Zs_plus + d.Zs_minus);
      }
   }

   // --- (Cross-rank Phase 4) SHARED (cross-rank) fault faces ----------------
   // Same eps-offset rule as the interior loop; the LOCAL side (Elem1) reads via
   // MaterialAtLocal_ and the PEER side (Elem2 = face-neighbour) via
   // MaterialAtNbr_.  in2 = -in1 ⇒ both sides at the SAME fault QP, offset ±eps
   // along the fault normal ⇒ EXACT symmetry for a depth profile on ANY mesh
   // (R-001/R-302), correct per-side contrast for a halfspace, GridFunction
   // supported.  Index spaces (R-204): shared_fault_elem1_on_plus_ is indexed by
   // POSITION si in fault_shared_faces_; shared_fault_dof_offset_ is keyed by the
   // RAW shared-face index sf (mirrors EvaluateBulkAtFaultQPsCanonical).
#ifdef MFEM_USE_MPI
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      auto &pmesh = static_cast<ParMesh &>(mesh_);
      const std::vector<bool> &sh_e1_plus = this->shared_fault_elem1_on_plus_;
      const std::map<int, int> &sh_dof_off = this->shared_fault_dof_offset_;
      for (int si = 0; si < fault_shared_faces_.Size(); ++si)
      {
         const int sf = fault_shared_faces_[si];
         FaceElementTransformations *ftr = pmesh.GetSharedFaceTransformations(sf);
         if (!ftr || ftr->Elem2 == nullptr) { continue; }
         if (si >= static_cast<int>(sh_e1_plus.size())) { continue; }
         const bool e1_plus = sh_e1_plus[si];
         auto off_it = sh_dof_off.find(sf);
         if (off_it == sh_dof_off.end()) { continue; }
         const int dof_off = off_it->second;

         // Physical centroids to orient the into-element direction (Elem1 local,
         // Elem2 the face-neighbour element transform — valid after the ctor's
         // ExchangeFaceNbrData).
         Vector c1(3), c2(3);
         ftr->Elem1->Transform(
            Geometries.GetCenter(ftr->Elem1->GetGeometryType()), c1);
         ftr->Elem2->Transform(
            Geometries.GetCenter(ftr->Elem2->GetGeometryType()), c2);

         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), qdeg);
         MFEM_VERIFY(ir.GetNPoints() == nbf,
                     "AssignFaultSidePerMaterialImpedances (shared): fault face sf="
                     << sf << " has " << ir.GetNPoints() << " QPs, expected "
                     "nbf_per_face=" << nbf << ".");

         for (int q = 0; q < nbf; ++q)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);

            Vector xf(3);
            ftr->Face->Transform(ip, xf);
            Vector nvec(3);
            CalcOrtho(ftr->Face->Jacobian(), nvec);
            const real_t nlen = nvec.Norml2();
            if (nlen > 0.0) { nvec /= nlen; }

            real_t d1 = 0.0;
            for (int d = 0; d < 3; ++d) { d1 += (c1(d) - xf(d)) * nvec(d); }
            const real_t s1 = (d1 >= 0.0) ? 1.0 : -1.0;
            const real_t in1[3] = { s1 * nvec(0), s1 * nvec(1), s1 * nvec(2) };
            const real_t in2[3] = { -in1[0], -in1[1], -in1[2] };

            const IntegrationPoint ip1 =
               offset_ip(*ftr->Elem1, ftr->GetElement1IntPoint(), in1);
            const IntegrationPoint ip2 =
               offset_ip(*ftr->Elem2, ftr->GetElement2IntPoint(), in2);

            real_t lam1, mu1, rho1, lam2, mu2, rho2;
            MaterialAtLocal_(ftr->Elem1No, *ftr->Elem1, ip1, lam1, mu1, rho1);
            MaterialAtNbr_(ftr, ip2, lam2, mu2, rho2);
            MFEM_VERIFY(mu1 > 0.0 && rho1 > 0.0 && mu2 > 0.0 && rho2 > 0.0,
                        "AssignFaultSidePerMaterialImpedances (shared): "
                        "non-physical material at fault sf=" << sf << " QP " << q);

            const real_t Zp1 = std::sqrt((lam1 + 2.0 * mu1) * rho1);
            const real_t Zs1 = std::sqrt(mu1 * rho1);
            const real_t Zp2 = std::sqrt((lam2 + 2.0 * mu2) * rho2);
            const real_t Zs2 = std::sqrt(mu2 * rho2);

            const int idx = dof_off + q;
            if (idx < 0 || idx >= static_cast<int>(dof_data.size())) { continue; }
            DOFData &d = dof_data[idx];
            if (e1_plus)
            { d.Zp_plus = Zp1; d.Zs_plus = Zs1; d.Zp_minus = Zp2; d.Zs_minus = Zs2; }
            else
            { d.Zp_plus = Zp2; d.Zs_plus = Zs2; d.Zp_minus = Zp1; d.Zs_minus = Zs1; }
            d.eta_p = d.Zp_plus * d.Zp_minus / (d.Zp_plus + d.Zp_minus);
            d.eta_s = d.Zs_plus * d.Zs_minus / (d.Zs_plus + d.Zs_minus);
         }
      }
   }
#endif
}
