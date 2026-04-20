// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// WaveOperator template implementation (included from wave_operator.hpp).

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------
template <typename MeshType>
WaveOperator<MeshType>::WaveOperator(MeshType &mesh, int order,
                                     real_t lambda, real_t mu, real_t rho,
                                     const BoundaryConfig &bc)
   : TimeDependentOperator(0),
     mesh_(mesh),
     order_(order),
     flux_(lambda, mu, rho),
     bc_(bc)
{
   int dim = mesh_.Dimension();
   MFEM_VERIFY(dim == 3, "WaveOperator requires 3D mesh, got dim=" << dim);

   fec_ = std::make_unique<L2_FECollection>(order, dim, BasisType::GaussLobatto);

   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      fes_ = std::make_unique<FESpaceType>(&mesh_, fec_.get());
#endif
   }
   else
   {
      fes_ = std::make_unique<FESpaceType>(&mesh_, fec_.get());
   }

   ne_ = mesh_.GetNE();
   ndof_per_el_ = fes_->GetFE(0)->GetDof();
   ndof_total_ = ne_ * ndof_per_el_;

   height = width = NUM_STATE * ndof_total_;

   flux_.BuildJacobian(0, Ax_);
   flux_.BuildJacobian(1, Ay_);
   flux_.BuildJacobian(2, Az_);

   AssembleElementMassInverse();

   // CFL length scale: inscribed diameter.
   // Hex: Vol^{1/3} (exact for cubes).
   // Tet: 6*Vol / total_face_area (exact inscribed diameter for any tet).
   //   Face areas computed from vertex cross products, no MFEM face API needed.
   h_min_ = std::numeric_limits<real_t>::max();
   for (int e = 0; e < ne_; e++)
   {
      real_t vol = mesh_.GetElementVolume(e);
      Geometry::Type geom = mesh_.GetElementGeometry(e);
      real_t h;
      if (geom == Geometry::TETRAHEDRON)
      {
         // Get the 4 tet vertices (non-owning Vector views into mesh storage;
         // safe because mesh data is stable during constructor).
         Array<int> vert;
         mesh_.GetElementVertices(e, vert);
         Vector v0(mesh_.GetVertex(vert[0]), 3);
         Vector v1(mesh_.GetVertex(vert[1]), 3);
         Vector v2(mesh_.GetVertex(vert[2]), 3);
         Vector v3(mesh_.GetVertex(vert[3]), 3);

         // Compute 4 face areas via cross products of edge vectors
         // Face 0: v1-v0, v2-v0  Face 1: v2-v0, v3-v0
         // Face 2: v3-v0, v1-v0  Face 3: v2-v1, v3-v1
         auto tri_area = [](const Vector &a, const Vector &b,
                            const Vector &c) -> real_t {
            real_t e1[3] = {b(0)-a(0), b(1)-a(1), b(2)-a(2)};
            real_t e2[3] = {c(0)-a(0), c(1)-a(1), c(2)-a(2)};
            real_t cx = e1[1]*e2[2] - e1[2]*e2[1];
            real_t cy = e1[2]*e2[0] - e1[0]*e2[2];
            real_t cz = e1[0]*e2[1] - e1[1]*e2[0];
            return 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
         };
         real_t A_total = tri_area(v0, v1, v2) + tri_area(v0, v2, v3)
                        + tri_area(v0, v3, v1) + tri_area(v1, v2, v3);
         // Inscribed diameter = 6*Vol / A_total
         h = (A_total > 0) ? 6.0 * vol / A_total : std::pow(vol, 1.0/3.0);
      }
      else
      {
         h = std::pow(vol, 1.0 / mesh_.Dimension());
      }
      h_min_ = std::min(h_min_, h);
   }

   // In parallel, reduce h_min across all ranks
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      real_t global_h_min;
      MPI_Allreduce(&h_min_, &global_h_min, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN,
                    mesh_.GetComm());
      h_min_ = global_h_min;
#endif
   }

   // Build face → boundary attribute map
   face_bdr_attr_.assign(mesh_.GetNumFaces(), 0);
   for (int b = 0; b < mesh_.GetNBE(); b++)
   {
      int face_idx = mesh_.GetBdrElementFaceIndex(b);
      face_bdr_attr_[face_idx] = mesh_.GetBdrAttribute(b);
   }

   num_fault_dofs_ = 0;

   // R-005 fix: Initialize face neighbor data structures ONCE in constructor.
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
      if (pfes)
      {
         auto &pmesh = const_cast<ParMesh &>(
            static_cast<const ParMesh &>(mesh_));
         pmesh.ExchangeFaceNbrData();
         pfes->ExchangeFaceNbrData();
         ghost_gf_ = std::make_unique<ParGridFunction>(pfes);
         ghost_initialized_ = true;
         // R-109 fix: cache MPI rank once here so
         // ComputeSharedFaceFluxRHS doesn't re-query it on every RK4 stage.
         my_rank_ = pmesh.GetMyRank();

         // R-001 fix: Identify shared faces that are fault faces.
         // A shared face is a fault face if it corresponds to a boundary
         // element with fault_attr on either rank. We check the local
         // boundary elements and match by face geometry (vertex keys).
         int n_shared = pmesh.GetNSharedFaces();
         shared_face_bdr_attr_.assign(n_shared, 0);

         // Global vertex IDs are the only stable key across ranks: local
         // vertex numbering is rank-private, so a shared face's local
         // verts on rank 0 and rank 1 are different integers.  MFEM's
         // GetGlobalVertexIndices assigns a consistent global ID per vertex.
         Array<HYPRE_BigInt> gvi;
         pmesh.GetGlobalVertexIndices(gvi);

         auto make_global_key = [&](const Array<int> &verts)
         {
            std::array<HYPRE_BigInt, 4> key = {0, 0, 0, 0};
            for (int v = 0; v < std::min(verts.Size(), 4); v++)
            {
               key[v] = gvi[verts[v]];
            }
            std::sort(key.begin(), key.begin() + verts.Size());
            return key;
         };

         // Build a set of local fault face *global* vertex keys.  When
         // the serial mesh tagged an interior face as a fault boundary
         // element (the TPV102 pattern via BooleanFragments, or the 2-tet
         // inline test), ParMesh partitioning keeps that BE on only ONE
         // of the two ranks that share the face.  So the key set must be
         // merged across ranks, otherwise the non-BE-owner rank would
         // silently misclassify its shared fault face as non-fault and
         // VerifySharedFaultDOFDataConsistency would see an unpaired
         // DOFData record (the R-305 abort condition).
         std::set<std::array<HYPRE_BigInt, 4>> global_fault_keys;
         for (int b = 0; b < mesh_.GetNBE(); b++)
         {
            if (mesh_.GetBdrAttribute(b) != bc_.fault_attr) { continue; }
            Array<int> verts;
            mesh_.GetBdrElementVertices(b, verts);
            global_fault_keys.insert(make_global_key(verts));
         }

         // MPI Allgatherv of local fault keys so every rank sees the
         // union — cheap (O(N_fault) per rank, no per-RK4-stage cost).
         {
            std::vector<HYPRE_BigInt> local_flat;
            local_flat.reserve(global_fault_keys.size() * 4);
            for (const auto &k : global_fault_keys)
            {
               for (int i = 0; i < 4; i++) { local_flat.push_back(k[i]); }
            }
            const int my_size = static_cast<int>(local_flat.size());
            int nprocs_ctor = 0;
            MPI_Comm_size(pmesh.GetComm(), &nprocs_ctor);
            std::vector<int> sizes(nprocs_ctor), displs(nprocs_ctor);
            MPI_Allgather(&my_size, 1, MPI_INT, sizes.data(), 1, MPI_INT,
                          pmesh.GetComm());
            int total = 0;
            for (int r = 0; r < nprocs_ctor; r++)
            { displs[r] = total; total += sizes[r]; }
            std::vector<HYPRE_BigInt> all_flat(total);
            MPI_Allgatherv(local_flat.data(), my_size,
                           MPITypeMap<HYPRE_BigInt>::mpi_type,
                           all_flat.data(), sizes.data(), displs.data(),
                           MPITypeMap<HYPRE_BigInt>::mpi_type,
                           pmesh.GetComm());
            for (int i = 0; i < total; i += 4)
            {
               std::array<HYPRE_BigInt, 4> key = {
                  all_flat[i], all_flat[i+1], all_flat[i+2], all_flat[i+3]
               };
               global_fault_keys.insert(key);
            }
         }

         // R-002 fix: Detect shared fault faces by matching each shared
         // face's global vertex key against the MPI-merged fault key set.
         //
         // We also record every shared face's mesh-face index in
         // shared_mesh_face_set_ so ComputeFaceFluxRHS can tell a real
         // domain-boundary face (Elem2No==-1 and NOT in the set) from
         // a partition-seam face (Elem2No==-1 BUT in the set) — the
         // latter are handled by ComputeSharedFaceFluxRHS and must
         // NOT be dispatched by the local-face loop.
         //
         // R-701 fix: peer-rank resolution and the SharedFacePeer table
         // are no longer needed.  The canonical-frame path in
         // ComputeSharedFaceFluxRHS reconstructs the pre-Step-5 frame
         // (identical on both ranks) from FaultBasis's sign_flipped bit,
         // so there is no owner/non-owner Evaluate split.  Both ranks
         // run Evaluate locally on bit-identical inputs.
         shared_mesh_face_set_.clear();
         for (int sf = 0; sf < n_shared; sf++)
         {
            FaceElementTransformations *ftr =
               pmesh.GetSharedFaceTransformations(sf);
            if (!ftr) { continue; }
            int face_idx = pmesh.GetSharedFace(sf);
            shared_mesh_face_set_.insert(face_idx);
            Array<int> verts;
            mesh_.GetFaceVertices(face_idx, verts);
            std::array<HYPRE_BigInt, 4> key = make_global_key(verts);
            shared_face_bdr_attr_[sf] = global_fault_keys.count(key)
                                      ? bc_.fault_attr : 0;
         }
      }
#endif
   }

   // Canonical fault-face geometry lists.  Same filter as
   // FaultBoundaryData::FindFaultInteriorFaces / FindFaultSharedFaces used
   // by BP5: interior list includes only truly 2-sided fault faces on this
   // rank (GetInteriorFaceTransformations non-null), shared list carries
   // the partition-seam counterparts separately.  Downstream code must
   // walk these two lists in order to keep fault_coords / DOFData / ParaView
   // geometry in lockstep.
   fault_interior_faces_.SetSize(0);
   if (bc_.fault_attr > 0)
   {
      for (int b = 0; b < mesh_.GetNBE(); b++)
      {
         if (mesh_.GetBdrAttribute(b) != bc_.fault_attr) { continue; }
         int face_idx = mesh_.GetBdrElementFaceIndex(b);
         if (mesh_.GetInteriorFaceTransformations(face_idx) == nullptr)
         { continue; }
         fault_interior_faces_.Append(face_idx);
      }
   }

   fault_shared_faces_.SetSize(0);
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto &pmesh = static_cast<const ParMesh &>(mesh_);
      for (int sf = 0; sf < pmesh.GetNSharedFaces(); sf++)
      {
         if (sf < static_cast<int>(shared_face_bdr_attr_.size()) &&
             shared_face_bdr_attr_[sf] == bc_.fault_attr)
         {
            fault_shared_faces_.Append(sf);
         }
      }
#endif
   }

   // R-801 fix: build mesh_face_idx → FaultBasis index map for interior
   // fault faces.  FaultBasis stores per-face data as [interior_faces...,
   // shared_faces...], so the position of face f in fault_interior_faces_
   // is directly usable as `fault_basis_->GetBasis(pos)`.  The
   // ComputeFaceFluxRHS interior-fault branch uses this to reconstruct the
   // BP5 canonical frame (same convention as the shared-fault branch) —
   // unifying the (t1, t2) = (dip, strike) semantics across every fault
   // code path so DOFData.V1/V2/tau1_corr/tau2_corr/slip1/slip2 carry the
   // SAME physical meaning on interior and shared fault QPs.
   fault_interior_face_to_basis_idx_.clear();
   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      fault_interior_face_to_basis_idx_[fault_interior_faces_[i]] = i;
   }

   // R-701 fix: populate BP5's FaultBasis for every fault face (interior
   // and shared).  The canonical pre-Step-5 frame — reconstructed below
   // in the flux routines from `FaultBasisQPData::sign_flipped` — is
   // identical on both ranks sharing a fault face.  This is what makes
   // the v5 R-501 owner-broadcast (and the entire v1–v5 (+,-) swap saga)
   // unnecessary: both ranks run `FaultFaceFlux::Evaluate` locally on
   // bit-identical inputs and produce bit-identical DOFData updates.
   //
   // Uses the same conventions as BP5's elasticity_operator_setup.inl:
   //   ref_normal = (0, -1, 0)  (Tandem convention — fault at y=0, normal -y)
   //   up         = (0,  0, 1)  (z-axis up).
   //
   // nbf_per_face_ is derived from face geometry here (not deferred to
   // SetFaultDOFData) so the per-QP FaultBasis can be populated at ctor
   // time.  The SetFaultDOFData verify still catches the case where the
   // driver configures a different nqp_per_face than the ctor computed.
   if (bc_.fault_attr > 0 &&
       (fault_interior_faces_.Size() > 0 || fault_shared_faces_.Size() > 0))
   {
      Vector ref_normal(3); ref_normal = 0.0; ref_normal(1) = -1.0;
      Vector up(3);         up = 0.0;         up(2) = 1.0;

      fault_basis_ = std::make_unique<FaultBasis>();
      // Always call Compute (even with an empty list) so FaultBasis's dim_
      // is initialised.  AppendSharedFaces / ComputeQPBasisShared use dim_
      // to size the `n_raw` vector passed to CalcOrtho; skipping Compute
      // on ranks with zero interior-fault faces but nonzero shared-fault
      // faces crashes in CalcOrtho.
      fault_basis_->Compute(static_cast<Mesh &>(mesh_),
                            fault_interior_faces_, ref_normal, up);
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         if (fault_shared_faces_.Size() > 0)
         {
            fault_basis_->AppendSharedFaces(
               static_cast<ParMesh &>(mesh_),
               fault_shared_faces_, ref_normal, up);
         }
#endif
      }

      // Determine face geometry + per-QP count from the first available
      // fault face.  Ranks may have only shared, only interior, or both.
      Geometry::Type face_geom = Geometry::TRIANGLE;
      if (fault_interior_faces_.Size() > 0)
      {
         auto *ftr0 = mesh_.GetInteriorFaceTransformations(
            fault_interior_faces_[0]);
         if (ftr0) { face_geom = ftr0->GetGeometryType(); }
      }
      else if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         if (fault_shared_faces_.Size() > 0)
         {
            auto &pmesh = static_cast<ParMesh &>(mesh_);
            auto *ftr0 = pmesh.GetSharedFaceTransformations(
               fault_shared_faces_[0]);
            if (ftr0) { face_geom = ftr0->GetGeometryType(); }
         }
#endif
      }
      const IntegrationRule &face_ir = IntRules.Get(face_geom, 2*order_);
      nbf_per_face_ = face_ir.GetNPoints();

      if (fault_interior_faces_.Size() > 0)
      {
         fault_basis_->ComputeQPBasis(static_cast<Mesh &>(mesh_),
                                      fault_interior_faces_,
                                      ref_normal, up, face_ir);
      }
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         if (fault_shared_faces_.Size() > 0)
         {
            fault_basis_->ComputeQPBasisShared(
               static_cast<ParMesh &>(mesh_),
               fault_shared_faces_, ref_normal, up, face_ir,
               fault_interior_faces_.Size());
         }
#endif
      }

      // R-701 swap-flag: determine per shared-fault face whether this rank's
      // Elem1 sits on the canonical "+" side (side FROM WHICH ref_normal
      // points AWAY — origin of the arrow).  We cannot use FaultBasis's
      // `sign_flipped` for the swap because MFEM's shared-face CalcOrtho
      // can produce the SAME raw nor on both ranks (the Elem1-outward
      // convention only holds strictly for serial interior faces; for
      // shared faces, depending on mesh topology, CalcOrtho may orient
      // globally rather than per-rank).  So we derive `elem1_on_plus_`
      // from geometry: compare Elem1 centroid vs face centroid along
      // ref_normal.  Exactly one of the two ranks sharing a fault face
      // will have elem1_on_plus_ == true; the other will have false.
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         shared_fault_elem1_on_plus_.assign(fault_shared_faces_.Size(), false);
         auto &pmesh = static_cast<ParMesh &>(mesh_);
         for (int sf_idx = 0; sf_idx < fault_shared_faces_.Size(); sf_idx++)
         {
            int sf = fault_shared_faces_[sf_idx];
            FaceElementTransformations *ftr =
               pmesh.GetSharedFaceTransformations(sf);
            if (!ftr) { continue; }

            // Face centroid.
            const IntegrationPoint &ip_center =
               Geometries.GetCenter(ftr->GetGeometryType());
            ftr->Face->SetIntPoint(&ip_center);
            Vector face_c(3);
            ftr->Face->Transform(ip_center, face_c);

            // Elem1 centroid (local element; vertex-average is sufficient
            // for tet/hex).
            Array<int> e1_verts;
            mesh_.GetElementVertices(ftr->Elem1No, e1_verts);
            Vector elem1_c(3);
            elem1_c = 0.0;
            for (int v = 0; v < e1_verts.Size(); v++)
            {
               const real_t *vp = mesh_.GetVertex(e1_verts[v]);
               for (int d = 0; d < 3; d++) { elem1_c(d) += vp[d]; }
            }
            if (e1_verts.Size() > 0)
            {
               elem1_c /= static_cast<real_t>(e1_verts.Size());
            }

            // Projections along ref_normal: canonical arrow points in
            // ref_normal direction.  "+ side" = origin half-space =
            // opposite to ref_normal direction.  A point with SMALLER
            // (more negative) projection onto ref_normal is on the
            // +side.  So elem1_on_plus = (elem1_proj < face_proj).
            real_t face_proj  = 0.0, elem1_proj = 0.0;
            for (int d = 0; d < 3; d++)
            {
               face_proj  += face_c(d)  * ref_normal(d);
               elem1_proj += elem1_c(d) * ref_normal(d);
            }
            shared_fault_elem1_on_plus_[sf_idx] = (elem1_proj < face_proj);
         }
#endif
      }
   }
}

template <typename MeshType>
WaveOperator<MeshType>::~WaveOperator() = default;

// ---------------------------------------------------------------------------
// Mult: dQ/dt = M^{-1} * (-Face + Vol)   [Eq. (3)]
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::Mult(const Vector &Q, Vector &dQdt) const
{
   MFEM_VERIFY(Q.Size() == height,
               "Q size mismatch: " << Q.Size() << " vs " << height);

   dQdt.SetSize(height);
   dQdt = 0.0;

   ComputeVolumeRHS(Q, dQdt);
   ComputeFaceFluxRHS(Q, dQdt);

   // Shared face flux for parallel (R-002 fix)
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      ComputeSharedFaceFluxRHS(Q, dQdt);
   }

   if (pml_layer_)
   {
      ApplyPMLDamping(Q, dQdt);
   }

   ApplyMassInverse(dQdt);
}

// ---------------------------------------------------------------------------
// Volume integral: Eq. (2c)
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ComputeVolumeRHS(const Vector &Q, Vector &rhs) const
{
   const real_t *Q_data = Q.GetData();

   for (int e = 0; e < ne_; e++)
   {
      const FiniteElement *fe = fes_->GetFE(e);
      ElementTransformation *Tr = fes_->GetElementTransformation(e);
      int ndof = fe->GetDof();

      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2*order_);
      int nqp = ir.GetNPoints();

      int dof_offset = e * ndof_per_el_;

      Vector shape(ndof);
      DenseMatrix dshape(ndof, 3);

      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetIntPoint(&ip);
         real_t w = ip.weight * Tr->Weight();

         fe->CalcShape(ip, shape);
         fe->CalcPhysDShape(*Tr, dshape);

         real_t Q_qp[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_qp[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               Q_qp[c] += shape(i) * Q_data[c * ndof_total_ + dof_offset + i];
            }
         }

         real_t F[3][NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            F[0][c] = 0.0; F[1][c] = 0.0; F[2][c] = 0.0;
            for (int k = 0; k < NUM_STATE; k++)
            {
               F[0][c] += Ax_(c, k) * Q_qp[k];
               F[1][c] += Ay_(c, k) * Q_qp[k];
               F[2][c] += Az_(c, k) * Q_qp[k];
            }
         }

         for (int c = 0; c < NUM_STATE; c++)
         {
            for (int i = 0; i < ndof; i++)
            {
               real_t val = 0.0;
               for (int j = 0; j < 3; j++)
               {
                  val += dshape(i, j) * F[j][c];
               }
               rhs[c * ndof_total_ + dof_offset + i] += w * val;
            }
         }
      }
   }
}

// ---------------------------------------------------------------------------
// Face flux: Eq. (2d) — local faces only
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ComputeFaceFluxRHS(const Vector &Q, Vector &rhs) const
{
   const real_t *Q_data = Q.GetData();

   for (int f = 0; f < mesh_.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh_.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      int e1 = ftr->Elem1No;
      int e2 = ftr->Elem2No;

      int bdr_attr = face_bdr_attr_[f];

      // Skip ALL faces that belong to another rank (partition seam).
      // MFEM leaves Elem2No==-1 on such faces, which makes them look
      // like true boundary faces even when they carry a bdr attribute
      // (e.g. a fault face at a partition seam is tagged with
      // fault_attr AND is a shared face).  Those are dispatched by
      // ComputeSharedFaceFluxRHS; dispatching them here as well would
      // both double-count the flux and mis-classify a fault face as
      // an outflow boundary.
      if (e2 < 0 && shared_mesh_face_set_.count(f) > 0) { continue; }

      bool is_boundary = (e2 < 0) && (bdr_attr > 0);

      // Regular interior shared face (not relevant in serial)
      if (e2 < 0 && bdr_attr == 0) { continue; }

      const FiniteElement *fe1 = fes_->GetFE(e1);
      int ndof = fe1->GetDof();
      int dof_offset1 = e1 * ndof_per_el_;

      const IntegrationRule &ir = IntRules.Get(
         ftr->GetGeometryType(), 2*order_);
      int nqp = ir.GetNPoints();

      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0) { nor_vec /= nor_len; }

         real_t w = ip.weight * nor_len;
         real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         real_t Q_self[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_self[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               Q_self[c] += shape1(i) * Q_data[c * ndof_total_ + dof_offset1 + i];
            }
         }

         real_t F_h[NUM_STATE];

         if (is_boundary)
         {
            FaceBC bc_type = ClassifyBoundaryFace(bdr_attr);

            switch (bc_type)
            {
               case FaceBC::Absorbing:
                  flux_.Absorbing(nor, Q_self, F_h);
                  break;
               case FaceBC::FreeSurface:
                  flux_.FreeSurface(nor, Q_self, F_h);
                  break;
               case FaceBC::Fault:
                  // A fault face must have an element on both sides.  Partition-seam
                  // fault faces are skipped above via shared_mesh_face_set_ — reaching
                  // this branch means the face is a true domain-boundary face (no
                  // neighbor exists anywhere, not on any rank) that the mesh file
                  // tagged with fault_attr.  This is a mesh configuration error:
                  // fault faces are interior physics and cannot be closed by an
                  // absorbing (or any) boundary condition.  Fail loudly rather
                  // than silently degrade to absorbing.
                  {
                     Array<int> fv;
                     const_cast<MeshType&>(mesh_).GetFaceVertices(f, fv);
                     double cx=0, cy=0, cz=0;
                     for (int iv = 0; iv < fv.Size(); iv++)
                     {
                        const double *x = mesh_.GetVertex(fv[iv]);
                        cx += x[0]; cy += x[1]; cz += x[2];
                     }
                     if (fv.Size() > 0) { cx/=fv.Size(); cy/=fv.Size(); cz/=fv.Size(); }
                     MFEM_ABORT("Fault face " << f << " centroid ("
                        << cx << ", " << cy << ", " << cz
                        << ") has no neighbor element on any rank.  Fault faces must "
                        "be interior (2-sided).  Re-check the mesh: every face tagged "
                        "with fault_attr must be an embedded surface inside the volume, "
                        "not a face on the domain boundary.");
                  }
                  break;
               default:
                  flux_.Absorbing(nor, Q_self, F_h);
                  break;
            }

            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  rhs[c * ndof_total_ + dof_offset1 + i] -= w * shape1(i) * F_h[c];
               }
            }
         }
         else
         {
            const FiniteElement *fe2 = fes_->GetFE(e2);
            int dof_offset2 = e2 * ndof_per_el_;

            IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2);
            Vector shape2(ndof);
            fe2->CalcShape(ip2, shape2);

            real_t Q_nbr[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               Q_nbr[c] = 0.0;
               for (int i = 0; i < ndof; i++)
               {
                  Q_nbr[c] += shape2(i) * Q_data[c * ndof_total_ + dof_offset2 + i];
               }
            }

            // Check if this interior face is a fault face (R-002 fix).
            // BooleanFragments meshes have fault at Y=0 as interior faces
            // with bdr_attr = fault_attr. Detect and mark for future
            // FaultFaceFlux dispatch.
            bool is_fault = (bdr_attr == bc_.fault_attr) && (bc_.fault_attr > 0);

            if (is_fault && fault_flux_ && fault_dof_data_)
            {
               // Fault face: dispatch to FaultFaceFlux with tracked DOFData.

               // R-001 fix: Look up the tracked DOFData for this face QP.
               auto it = fault_face_dof_offset_.find(f);
               int dof_idx = -1;
               if (it != fault_face_dof_offset_.end())
               {
                  dof_idx = it->second + q;
               }

               // R-801 fix: look up this interior fault face in FaultBasis
               // so we can use the same BP5 canonical (dip, strike) frame
               // convention as ComputeSharedFaceFluxRHS.  Without this
               // unification, interior-fault QPs store DOFData.V1 =
               // along-strike (from GodunovFlux::BuildFrame's t1=x) while
               // shared-fault QPs store DOFData.V1 = along-dip (from BP5's
               // tangent1 = dip), corrupting the TPV102 pure-strike-slip
               // initialisation and mixing the physical components.
               const int fb_idx = LookupInteriorFaultBasisIndex(f);
               const bool have_basis = (fault_basis_ && fb_idx >= 0 &&
                                        fb_idx < fault_basis_->NumFaces());
               const FaultBasisQPData *qpd_ptr = nullptr;
               if (have_basis)
               {
                  const FaultBasisData &bd = fault_basis_->GetBasis(fb_idx);
                  if (q < static_cast<int>(bd.qp_data.size()))
                  {
                     qpd_ptr = &bd.qp_data[q];
                  }
               }

               if (dof_idx >= 0 &&
                   dof_idx < static_cast<int>(fault_dof_data_->size()) &&
                   qpd_ptr != nullptr)
               {
                  DOFData &fdata = (*fault_dof_data_)[dof_idx];

                  // 1. Reconstruct canonical (pre-Step-5) BP5 frame from
                  // (stored_basis, sign_flipped) — ref-normal-aligned.
                  // BP5 convention: can_t1 = dip, can_t2 = strike.
                  const FaultBasisQPData &qpd = *qpd_ptr;
                  real_t can_n[3], can_t1[3], can_t2[3];
                  for (int d = 0; d < 3; d++)
                  {
                     can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]
                                                  :  qpd.normal[d];
                     can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d]
                                                  :  qpd.tangent1[d];
                     can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d]
                                                  :  qpd.tangent2[d];
                  }

                  DenseMatrix T_can(NUM_STATE), Tinv_can(NUM_STATE);
                  GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);
                  GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2,
                                                    Tinv_can);

                  // For interior faces MFEM's CalcOrtho returns Elem1-outward
                  // normal; sign_flipped = true iff that is anti-aligned with
                  // ref_normal.  can_n points along ref_normal (+→−).  So
                  // Elem1 sits on the canonical + side iff its outward lies
                  // along can_n iff sign_flipped == false.
                  const bool elem1_on_plus = !qpd.sign_flipped;

                  // Rotate self/nbr into canonical frame.
                  real_t Q_self_can[NUM_STATE], Q_nbr_can[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     Q_self_can[c] = 0.0; Q_nbr_can[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        Q_self_can[c] += Tinv_can(c, k) * Q_self[k];
                        Q_nbr_can[c]  += Tinv_can(c, k) * Q_nbr[k];
                     }
                  }

                  const real_t *Q_plus_local  = elem1_on_plus ? Q_self_can
                                                              : Q_nbr_can;
                  const real_t *Q_minus_local = elem1_on_plus ? Q_nbr_can
                                                              : Q_self_can;

                  // 2. Evaluate: trial traction → friction solve → imposed states
                  real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
                  fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
                                        Q_imp_plus, Q_imp_minus);

                  // 3. Rotate imposed states back to global via T_can.
                  real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     Q_imp_plus_g[c] = 0.0; Q_imp_minus_g[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        Q_imp_plus_g[c]  += T_can(c, k) * Q_imp_plus[k];
                        Q_imp_minus_g[c] += T_can(c, k) * Q_imp_minus[k];
                     }
                  }

                  // Godunov flux from imposed states.  Standard DG
                  // accumulation: Elem1 -= F, Elem2 += F.  MFEM's `nor` is
                  // Elem1-outward on interior faces, so this element's
                  // Q_self comes from the (+) side iff elem1_on_plus is
                  // true.  Route Q_imp_plus_g/Q_imp_minus_g to self/nbr
                  // slots accordingly; the Godunov identity
                  // `F(L,R,+n) = -F(R,L,-n)` then makes the Elem1/Elem2
                  // accumulation conservative.
                  const real_t *Q_self_imp = elem1_on_plus ? Q_imp_plus_g
                                                           : Q_imp_minus_g;
                  const real_t *Q_nbr_imp  = elem1_on_plus ? Q_imp_minus_g
                                                           : Q_imp_plus_g;
                  real_t F_h_total[NUM_STATE];
                  flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h_total);

                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     for (int i = 0; i < ndof; i++)
                     {
                        rhs[c * ndof_total_ + dof_offset1 + i] -=
                           w * shape1(i) * F_h_total[c];
                     }
                  }
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     for (int i = 0; i < ndof; i++)
                     {
                        rhs[c * ndof_total_ + dof_offset2 + i] +=
                           w * shape2(i) * F_h_total[c];
                     }
                  }
               }
               else
               {
                  // Fault face but no DOFData mapping (or FaultBasis not
                  // populated) — fall back to welded interior flux.
                  flux_.Interior(nor, Q_self, Q_nbr, F_h);
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     for (int i = 0; i < ndof; i++)
                     {
                        rhs[c * ndof_total_ + dof_offset1 + i] -=
                           w * shape1(i) * F_h[c];
                        rhs[c * ndof_total_ + dof_offset2 + i] +=
                           w * shape2(i) * F_h[c];
                     }
                  }
               }
            }
            else
            {
               // Regular interior face: standard Godunov flux
               flux_.Interior(nor, Q_self, Q_nbr, F_h);

               for (int c = 0; c < NUM_STATE; c++)
               {
                  for (int i = 0; i < ndof; i++)
                  {
                     rhs[c * ndof_total_ + dof_offset1 + i] -=
                        w * shape1(i) * F_h[c];
                     rhs[c * ndof_total_ + dof_offset2 + i] +=
                        w * shape2(i) * F_h[c];
                  }
               }
            }
         }
      }
   }
}

// ---------------------------------------------------------------------------
// Shared face flux: R-002 fix — wave propagation across rank boundaries
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ComputeSharedFaceFluxRHS(
   const Vector &Q, Vector &rhs) const
{
   if constexpr (!IsParallelMesh<MeshType>::value)
   {
      return;  // No shared faces in serial
   }
   else
   {
#ifdef MFEM_USE_MPI
      auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
      MFEM_VERIFY(pfes, "FESpace must be ParFiniteElementSpace for ParMesh");

      auto &pmesh = static_cast<const ParMesh &>(mesh_);
      int n_shared = pmesh.GetNSharedFaces();
      if (n_shared == 0) { return; }

      // R-001 fix: canonicalize (+,-) side at shared fault faces.  Both ranks
      // share the same physical QP but each stores an independent DOFData
      // entry.  Without the swap, Evaluate is called with opposite (+,-)
      // arguments on the two ranks — the velocity-jump term in Eq. 7b flips
      // sign while the stress term does not, so tau1_trial (and every
      // downstream quantity) diverges and the two DOFData entries drift
      // apart step by step.  The lower rank ID is the canonical "+" owner;
      // the peer swaps (Q_plus, Q_minus) and the Q_imp output buffers so
      // both ranks drive Evaluate with identical inputs and end up with
      // bit-identical DOFData updates.
      // R-109 fix: MPI rank is cached in the ctor to avoid per-call queries.

      // Exchange ghost Q data: one component at a time via persistent
      // ParGridFunction (R-003 fix: initialized once in constructor).
      const real_t *Q_data = Q.GetData();
      MFEM_VERIFY(ghost_gf_, "Ghost GF not initialized");

      ParGridFunction &q_gf = *ghost_gf_;
      std::vector<Vector> nbr_data(NUM_STATE);

      for (int c = 0; c < NUM_STATE; c++)
      {
         for (int i = 0; i < ndof_total_; i++)
         {
            q_gf[i] = Q_data[c * ndof_total_ + i];
         }
         q_gf.ExchangeFaceNbrData();
         // R-005 fix: force own-storage deep copy. Vector::operator= can take
         // a shallow branch; without this the per-component ghost data may
         // all alias q_gf.FaceNbrData(), which is overwritten by the NEXT
         // ExchangeFaceNbrData call — every nbr_data[c] would then end up
         // pointing at the last component's ghost buffer (silent aliasing).
         const Vector &src = q_gf.FaceNbrData();
         nbr_data[c].SetSize(src.Size());
         std::memcpy(nbr_data[c].GetData(), src.GetData(),
                     src.Size() * sizeof(real_t));
         // R-302 Part B (refined by R-405): MFEM_VERIFY in Release to catch
         // a future refactor that replaces SetSize+memcpy with a shallow
         // `operator=` (the H2 regression).  Every component uses the same
         // allocator pathway, so a single check on the first iteration
         // gives the same protection without per-Mult-per-component branch
         // overhead in the hot loop.
         if (c == 0)
         {
            MFEM_VERIFY(nbr_data[c].GetData() != q_gf.FaceNbrData().GetData(),
                        "nbr_data[0] aliases q_gf.FaceNbrData() — H2 "
                        "regression.  The deep-copy invariant is load-"
                        "bearing for cross-rank bulk wave propagation and "
                        "must hold in Release.");
         }
      }

      // R-701 fix: canonical-fault-frame refactor.
      //
      // Before R-701, shared-fault face processing required owner broadcast
      // of post-Evaluate state (v5 R-501) because each rank built its own
      // fault-local frame via `GodunovFlux::BuildFrame(rank_local_nor, ...)`
      // and the two ranks' rotation matrices disagreed → Evaluate inputs
      // disagreed → DOFData drifted.  Under R-701, the frame comes from
      // BP5's `FaultBasis` (populated in the ctor), which exposes a per-QP
      // `sign_flipped` bit.  Reconstructing the pre-Step-5 canonical frame
      // from `(stored_basis, sign_flipped)` gives a (normal, tangent1,
      // tangent2) triple that is BIT-IDENTICAL on both ranks (it is the
      // ref-normal-aligned frame both ranks compute in BP5's Step 3/4
      // before the Step-5 negation flips the stored copy).
      //
      // Consequence: both ranks feed `Evaluate` with bit-identical arguments
      // and produce bit-identical DOFData updates locally.  No MPI
      // broadcast, no owner/non-owner split, no peer-rank table.  The
      // ~270 lines of v5 R-501 machinery are deleted.
      //
      // Per-rank self/nbr labeling still differs (rank A's Elem1 is on
      // canonical-+ side, rank B's on canonical-- side).  The swap logic
      // comes from `qpd.sign_flipped`: sign_flipped == true identifies the
      // rank whose Elem1 is on the canonical-- side (rank B).  Flux
      // assembly uses each rank's own MFEM `nor` with its own
      // (Q_self_imp, Q_nbr_imp) assignment; the Godunov conservation
      // identity `F(L, R, +n) = -F(R, L, -n)` then gives consistent
      // accumulation on both ranks without an explicit sign flip.
      const int nfs = fault_shared_faces_.Size();
      const bool fault_active = (fault_flux_ && fault_dof_data_
                                 && nfs > 0 && nbf_per_face_ > 0
                                 && fault_basis_);

      // Map each shared-face index sf to its FaultBasis entry.  BP5 lays
      // out `basis_[]` as (interior-fault faces, then shared-fault faces);
      // the shared half starts at offset `fault_interior_faces_.Size()`.
      std::vector<int> sf_to_basis_idx(n_shared, -1);
      if (fault_active)
      {
         for (int sf_idx = 0; sf_idx < nfs; sf_idx++)
         {
            int sf = fault_shared_faces_[sf_idx];
            if (sf >= 0 && sf < n_shared)
            {
               sf_to_basis_idx[sf] = fault_interior_faces_.Size() + sf_idx;
            }
         }
      }

      for (int sf = 0; sf < n_shared; sf++)
      {
         FaceElementTransformations *ftr =
            const_cast<ParMesh &>(pmesh).GetSharedFaceTransformations(sf);
         if (!ftr) { continue; }

         int e1 = ftr->Elem1No;
         const FiniteElement *fe1 = fes_->GetFE(e1);
         int ndof = fe1->GetDof();
         int dof_offset1 = e1 * ndof_per_el_;

         // Ghost element: index in face neighbor block
         int nbr_idx = ftr->Elem2No - ne_;
         const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
         int ndof2 = fe2->GetDof();

         const IntegrationRule &ir = IntRules.Get(
            ftr->GetGeometryType(), 2*order_);

         const bool sf_fault =
            fault_active && (sf < static_cast<int>(shared_face_bdr_attr_.size()))
            && (shared_face_bdr_attr_[sf] == bc_.fault_attr)
            && (bc_.fault_attr > 0)
            && (sf_to_basis_idx[sf] >= 0);
         const int basis_idx = sf_fault ? sf_to_basis_idx[sf] : -1;

         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);

            Vector nor_vec(3);
            CalcOrtho(ftr->Face->Jacobian(), nor_vec);
            real_t nor_len = nor_vec.Norml2();
            if (nor_len > 0) { nor_vec /= nor_len; }

            real_t w = ip.weight * nor_len;
            real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

            // Evaluate Q on local element (Elem1)
            IntegrationPoint ip1;
            ftr->Loc1.Transform(ip, ip1);
            Vector shape1(ndof);
            fe1->CalcShape(ip1, shape1);

            real_t Q_self[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               Q_self[c] = 0.0;
               for (int i = 0; i < ndof; i++)
               {
                  Q_self[c] += shape1(i) * Q_data[c * ndof_total_ + dof_offset1 + i];
               }
            }

            // Evaluate Q on ghost element (Elem2) from exchanged data
            IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2);
            Vector shape2(ndof2);
            fe2->CalcShape(ip2, shape2);

            MFEM_ASSERT(ndof2 == ndof_per_el_,
                        "Mixed element types: ghost ndof " << ndof2
                        << " != ndof_per_el " << ndof_per_el_);

            real_t Q_nbr[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               Q_nbr[c] = 0.0;
               for (int i = 0; i < ndof2; i++)
               {
                  Q_nbr[c] += shape2(i) * nbr_data[c][nbr_idx * ndof_per_el_ + i];
               }
            }

            real_t F_h[NUM_STATE];

            if (sf_fault)
            {
               auto it = shared_fault_dof_offset_.find(sf);
               int dof_idx = (it != shared_fault_dof_offset_.end())
                           ? it->second + q : -1;

               if (dof_idx >= 0 &&
                   dof_idx < static_cast<int>(fault_dof_data_->size()))
               {
                  // Look up the rank's canonical-side flag for this
                  // shared fault face (computed once in the ctor from
                  // Elem1 geometry).  sf_idx in fault_shared_faces_:
                  const int sf_idx_in_fault = basis_idx
                                              - fault_interior_faces_.Size();
                  const bool elem1_on_plus =
                     (sf_idx_in_fault >= 0 && sf_idx_in_fault <
                      static_cast<int>(shared_fault_elem1_on_plus_.size()))
                     ? shared_fault_elem1_on_plus_[sf_idx_in_fault] : false;

                  const FaultBasisData &bd = fault_basis_->GetBasis(basis_idx);
                  MFEM_ASSERT(q < static_cast<int>(bd.qp_data.size()),
                              "FaultBasis::qp_data not populated for shared "
                              "fault face — ComputeQPBasisShared missed "
                              "this face");
                  const FaultBasisQPData &qpd = bd.qp_data[q];

                  // Reconstruct canonical (pre-Step-5) frame: aligned with
                  // ref_normal on BOTH ranks → bit-identical across ranks
                  // regardless of whether MFEM's CalcOrtho gave identical
                  // or opposite raw normals.
                  real_t can_n[3], can_t1[3], can_t2[3];
                  for (int d = 0; d < 3; d++)
                  {
                     can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]
                                                  :  qpd.normal[d];
                     can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d]
                                                  :  qpd.tangent1[d];
                     can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d]
                                                  :  qpd.tangent2[d];
                  }

                  DenseMatrix T_can(NUM_STATE), Tinv_can(NUM_STATE);
                  GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);
                  GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2,
                                                    Tinv_can);

                  // Rotate self/nbr Q into canonical frame.
                  real_t Q_self_can[NUM_STATE], Q_nbr_can[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     Q_self_can[c] = 0.0; Q_nbr_can[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        Q_self_can[c] += Tinv_can(c, k) * Q_self[k];
                        Q_nbr_can[c]  += Tinv_can(c, k) * Q_nbr[k];
                     }
                  }

                  // Swap based on geometric `elem1_on_plus_`: exactly
                  // one of the two ranks has elem1_on_plus=true, the
                  // other false.  Both ranks feed Evaluate with the
                  // same (Q_plus, Q_minus) = (Q at canonical-+ side,
                  // Q at canonical-- side) — bit-identical across
                  // ranks, because Q_self_can ⊕ Q_nbr_can across ranks
                  // exhaust the same two physical values (one rank's
                  // self = other rank's nbr by MFEM ghost exchange).
                  const real_t *Q_plus_local  = elem1_on_plus ? Q_self_can
                                                              : Q_nbr_can;
                  const real_t *Q_minus_local = elem1_on_plus ? Q_nbr_can
                                                              : Q_self_can;

                  DOFData &fdata = (*fault_dof_data_)[dof_idx];
                  real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
                  fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
                                        Q_imp_plus, Q_imp_minus);

                  // Rotate imposed states back to global via T_can (same
                  // on both ranks).  Q_imp_plus_g, Q_imp_minus_g are
                  // bit-identical across ranks.
                  real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     Q_imp_plus_g[c] = 0.0; Q_imp_minus_g[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        Q_imp_plus_g[c]  += T_can(c, k) * Q_imp_plus[k];
                        Q_imp_minus_g[c] += T_can(c, k) * Q_imp_minus[k];
                     }
                  }

                  // R-802 fix: use canonical normal `can_n` (bit-identical
                  // on both ranks, by construction) for the Godunov flux so
                  // both ranks compute the SAME F_h regardless of whether
                  // MFEM's CalcOrtho gives `nor_A = -nor_B` (classical) or
                  // `nor_A = nor_B` (observed empirically on the inline
                  // 2-tet mesh in v6 fix report).  The accumulation sign is
                  // then gated on `elem1_on_plus`:
                  //   Elem1 on canonical "+" side:  rhs[Elem1] -= F_can*w
                  //     (flux leaves Elem1 into the "-" half-space).
                  //   Elem1 on canonical "-" side:  rhs[Elem1] += F_can*w
                  //     (flux enters Elem1 from the "+" half-space).
                  // Pre-R-802 this code fed `nor` + self/nbr to Interior and
                  // relied on the Godunov identity F(L,R,+n) = -F(R,L,-n) to
                  // produce conservation — which is true only when
                  // nor_A = -nor_B.  When MFEM returns identical normals on
                  // the two ranks (v6 observation) the identity does not
                  // apply and bulk momentum leaks across the shared fault.
                  flux_.Interior(can_n, Q_imp_plus_g, Q_imp_minus_g, F_h);
                  const real_t accum_sign = elem1_on_plus ? +1.0 : -1.0;
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     for (int i = 0; i < ndof; i++)
                     {
                        rhs[c * ndof_total_ + dof_offset1 + i] -=
                           accum_sign * w * shape1(i) * F_h[c];
                     }
                  }
                  continue;  // Fault accumulation done; skip the generic
                             // -= F_h path below (which assumes rank-local
                             // `nor` + standard Godunov orientation).
               }
               else
               {
                  flux_.Interior(nor, Q_self, Q_nbr, F_h);
               }
            }
            else
            {
               flux_.Interior(nor, Q_self, Q_nbr, F_h);
            }

            // Accumulate into local element only (no ghost writes).
            // Non-fault shared faces: standard Godunov conservation holds
            // because both ranks use rank-local `nor` (opposite) on
            // bit-identical Q_self/Q_nbr inputs (each rank's `Q_self` is
            // the other rank's `Q_nbr`, and ghost exchange enforces
            // agreement).  For the fault path we use a separate
            // canonical-normal assembly + explicit accum_sign gated on
            // elem1_on_plus (see R-802 block above and its `continue`).
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  rhs[c * ndof_total_ + dof_offset1 + i] -= w * shape1(i) * F_h[c];
               }
            }
         }
      }
#endif
   }
}

// ---------------------------------------------------------------------------
// Apply per-element inverse mass matrix
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ApplyMassInverse(Vector &dQdt) const
{
   Vector elem_rhs(ndof_per_el_);
   Vector elem_result(ndof_per_el_);

   for (int e = 0; e < ne_; e++)
   {
      int dof_offset = e * ndof_per_el_;
      const DenseMatrix &Minv = elem_mass_inv_[e];

      for (int c = 0; c < NUM_STATE; c++)
      {
         for (int i = 0; i < ndof_per_el_; i++)
         {
            elem_rhs(i) = dQdt[c * ndof_total_ + dof_offset + i];
         }

         Minv.Mult(elem_rhs, elem_result);

         for (int i = 0; i < ndof_per_el_; i++)
         {
            dQdt[c * ndof_total_ + dof_offset + i] = elem_result(i);
         }
      }
   }
}

// ---------------------------------------------------------------------------
// Assemble per-element inverse mass matrices
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::AssembleElementMassInverse()
{
   elem_mass_inv_.resize(ne_);

   for (int e = 0; e < ne_; e++)
   {
      const FiniteElement *fe = fes_->GetFE(e);
      ElementTransformation *Tr = fes_->GetElementTransformation(e);
      int ndof = fe->GetDof();

      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2*order_);

      DenseMatrix M(ndof);
      M = 0.0;

      Vector shape(ndof);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetIntPoint(&ip);
         real_t w = ip.weight * Tr->Weight();

         fe->CalcShape(ip, shape);
         AddMult_a_VVt(w, shape, M);
      }

      elem_mass_inv_[e].SetSize(ndof);
      DenseMatrixInverse M_solver(M);
      M_solver.GetInverseMatrix(elem_mass_inv_[e]);
   }
}

// ---------------------------------------------------------------------------
// CFL time step
// ---------------------------------------------------------------------------
template <typename MeshType>
real_t WaveOperator<MeshType>::ComputeMaxDt(real_t cfl) const
{
   return cfl * h_min_ / flux_.GetCp();
}

// ---------------------------------------------------------------------------
// PML damping: rhs -= d(x) * D * Q  (Eq. 16)
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ApplyPMLDamping(const Vector &Q, Vector &rhs) const
{
   const real_t *Q_data = Q.GetData();

   for (int e = 0; e < ne_; e++)
   {
      const FiniteElement *fe = fes_->GetFE(e);
      ElementTransformation *Tr = fes_->GetElementTransformation(e);
      int ndof = fe->GetDof();
      int dof_offset = e * ndof_per_el_;

      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2*order_);

      Vector shape(ndof);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetIntPoint(&ip);
         real_t w = ip.weight * Tr->Weight();

         fe->CalcShape(ip, shape);

         Vector phys(3);
         Tr->Transform(ip, phys);

         real_t dx, dy, dz;
         pml_layer_->ComputeDamping(phys(0), phys(1), phys(2), dx, dy, dz);

         if (dx <= 0.0 && dy <= 0.0 && dz <= 0.0) { continue; }

         real_t Q_qp[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_qp[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               Q_qp[c] += shape(i) * Q_data[c * ndof_total_ + dof_offset + i];
            }
         }

         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t d_c = dx * PMLLayer::Dx[c]
                       + dy * PMLLayer::Dy[c]
                       + dz * PMLLayer::Dz[c];
            if (d_c > 0.0)
            {
               for (int i = 0; i < ndof; i++)
               {
                  rhs[c * ndof_total_ + dof_offset + i] -= w * d_c * shape(i) * Q_qp[c];
               }
            }
         }
      }
   }
}

// ---------------------------------------------------------------------------
// R-101 fix: verify shared-fault DOFData consistency across ranks.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::VerifySharedFaultDOFDataConsistency(
   real_t tol) const
{
   if constexpr (!IsParallelMesh<MeshType>::value)
   {
      return;  // Serial: nothing to compare.
   }
   else
   {
#ifdef MFEM_USE_MPI
      auto &pmesh = const_cast<ParMesh &>(static_cast<const ParMesh &>(mesh_));
      MPI_Comm comm = pmesh.GetComm();
      int nprocs;
      MPI_Comm_size(comm, &nprocs);

      // Short-circuit globally only: every rank must still participate in
      // the MPI_Allgather below even if it has nothing to contribute, or
      // ranks with shared fault faces would deadlock waiting on ranks
      // without any.
      int any_shared_local = (fault_dof_data_ &&
                              fault_shared_faces_.Size() > 0) ? 1 : 0;
      int any_shared_global = 0;
      MPI_Allreduce(&any_shared_local, &any_shared_global, 1, MPI_INT,
                    MPI_MAX, comm);
      if (!any_shared_global) { return; }

      // Pack per-QP records: centroid + every DOFData field that
      // FaultFaceFlux::Evaluate or the driver-side RK4-averaging step
      // writes.  One record = 3 centroid coords + 8 mutable fields +
      // 1 emitting-rank tag = 12 doubles.  Each rank contributes one
      // record per shared fault QP it owns; the same physical QP appears
      // on both ranks that share the fault face.  Any field whose values
      // disagree between the two owners is flagged as an R-101
      // consistency failure.  The emitting-rank tag is used in the
      // diagnostic print for unpaired entries so we can identify which
      // rank contributed a singleton.
      constexpr int REC = 12;
      constexpr int FIELD_BASE = 3;
      constexpr int NUM_FIELDS = 8;
      constexpr int RANK_OFFSET = 11;   // FIELD_BASE + NUM_FIELDS
      static const char *FIELD_NAMES[NUM_FIELDS] = {
         "tau1_corr", "tau2_corr", "sigma_n_corr",
         "V1", "V2", "psi", "slip1", "slip2"
      };
      std::vector<double> local_data;

      for (int sf_idx = 0; sf_idx < fault_shared_faces_.Size(); sf_idx++)
      {
         int sf = fault_shared_faces_[sf_idx];
         FaceElementTransformations *ftr = pmesh.GetSharedFaceTransformations(sf);
         if (!ftr) { continue; }

         auto it = shared_fault_dof_offset_.find(sf);
         if (it == shared_fault_dof_offset_.end()) { continue; }
         int dof_base = it->second;

         const IntegrationRule &ir =
            IntRules.Get(ftr->GetGeometryType(), 2*order_);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            Vector phys(3);
            ftr->Face->Transform(ip, phys);

            int idx = dof_base + q;
            if (idx < 0 ||
                idx >= static_cast<int>(fault_dof_data_->size())) { continue; }
            const DOFData &d = (*fault_dof_data_)[idx];

            local_data.push_back(phys(0));
            local_data.push_back(phys(1));
            local_data.push_back(phys(2));
            local_data.push_back(d.tau1_corr);
            local_data.push_back(d.tau2_corr);
            local_data.push_back(d.sigma_n_corr);
            local_data.push_back(d.V1);
            local_data.push_back(d.V2);
            local_data.push_back(d.psi);
            local_data.push_back(d.slip1);
            local_data.push_back(d.slip2);
            local_data.push_back(static_cast<double>(my_rank_));
         }
      }

      // Allgatherv.  Use int sizes/displacements; total record count
      // fits easily for a TPV102-scale (< 1e5 shared fault QPs globally).
      int my_size = static_cast<int>(local_data.size());
      std::vector<int> sizes(nprocs), displs(nprocs);
      MPI_Allgather(&my_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, comm);
      int total = 0;
      for (int r = 0; r < nprocs; r++) { displs[r] = total; total += sizes[r]; }
      std::vector<double> all_data(total);
      MPI_Allgatherv(local_data.data(), my_size, MPI_DOUBLE,
                     all_data.data(), sizes.data(), displs.data(),
                     MPI_DOUBLE, comm);

      int n_entries = total / REC;
      if (n_entries == 0) { return; }

      // Sort entries lexicographically by exact (cx, cy, cz).  Both ranks
      // that own a shared fault QP compute the centroid by evaluating the
      // same MFEM Face transformation at the same reference IntegrationPoint,
      // so matching pairs agree bit-for-bit.  Using exact equality (a) avoids
      // the strict-weak-ordering UB of any tolerance-based comparator, and
      // (b) exposes topological mismatches (QPs owned by only one rank) as
      // unpaired entries instead of silently fusing them into 3-element
      // "groups" via transitive tolerance.
      std::vector<int> idx(n_entries);
      std::iota(idx.begin(), idx.end(), 0);
      std::sort(idx.begin(), idx.end(), [&](int a, int b)
      {
         for (int k = 0; k < 3; k++)
         {
            double va = all_data[a*REC + k], vb = all_data[b*REC + k];
            if (va != vb) { return va < vb; }
         }
         return false;
      });

      // Track both absolute and relative diff so the diagnostic can report
      // the physical units while thresholding against a scale-relative
      // tolerance (R-502).
      double max_rel_diff = 0.0;
      double max_diff_abs = 0.0;
      double max_diff_scale = 1.0;
      int    max_diff_field = -1;
      int    max_diff_entry = -1;
      int    n_pairs = 0;
      int    n_unpaired = 0;

      auto same_centroid = [&](int a, int b)
      {
         // R-503 fix: hybrid tolerance — `scale * DBL_EPSILON` alone
         // collapses to ~1e-28 at coordinates near zero (fault plane at
         // y=0 is the TPV102 case), so two ranks with 1-ULP-of-nonzero
         // drift in y get classified as different centroids and trigger
         // a spurious R-305 unpair-abort.
         //
         // v8.0.0 Phase 1 fix (job 7666171 diagnosis): the prior 1e-9 floor
         // was too tight.  At 200 m production mesh with 50-rank partition,
         // MFEM's shared-face `ftr->Face->Transform(ip, phys)` produces FP
         // drift > 1e-9 between rank pairs for some QP configurations
         // (vertex coordinate replication through face_nbr data is not
         // bit-identical to the local-face path through all aggregation
         // orders).  Two orphaned entries were observed with centroids
         // matching to 10 printed digits yet failing `same_centroid`.
         //
         // Raising the floor to 1e-6 m keeps us safely above any observed
         // inter-rank FP drift while remaining 8 orders of magnitude below
         // the smallest realistic fault-mesh element (min h ≈ 100 m on
         // production BP5/TPV102 meshes).  Two distinct fault faces have
         // centroids separated by O(h/2) ≈ 100 m, so 1e-6 m cannot cause
         // false grouping of unrelated faces.
         const double abs_floor = 1e-6;
         for (int k = 0; k < 3; k++)
         {
            double va = all_data[a*REC + k], vb = all_data[b*REC + k];
            double scale = std::max(std::abs(va), std::abs(vb));
            double tol_k = std::max(scale
                                    * std::numeric_limits<double>::epsilon(),
                                    abs_floor);
            if (std::abs(va - vb) > tol_k) { return false; }
         }
         return true;
      };

      // R-101 diagnostic enhancement (v8.0.0 Phase 1): collect unpaired
      // entries' centroid + emitting rank so the abort message identifies
      // WHICH QPs are orphaned.  Needed to distinguish (a) a rank
      // asymmetrically classifying one face as fault (both-side symmetric
      // orphan → 2 singletons of DIFFERENT centroids) from (b) a
      // transitive centroid collision (3-way group collapsing to one).
      // Capped at 32 entries so the diagnostic stays bounded.
      struct UnpairedEntry
      {
         double cx, cy, cz;
         int    rank;
         int    group_size;
      };
      std::vector<UnpairedEntry> unpaired;
      constexpr int MAX_UNPAIRED_REPORT = 32;

      int i = 0;
      while (i < n_entries)
      {
         int j = i + 1;
         while (j < n_entries && same_centroid(idx[i], idx[j])) { j++; }
         int group_size = j - i;
         if (group_size == 2)
         {
            int a = idx[i], b = idx[i+1];
            // Guard: both entries of a valid pair must come from DIFFERENT
            // ranks (same rank emitting twice for the same centroid is a
            // ctor bug — dof_offset double-mapping or similar).
            int rank_a = static_cast<int>(all_data[a*REC + RANK_OFFSET]);
            int rank_b = static_cast<int>(all_data[b*REC + RANK_OFFSET]);
            if (rank_a == rank_b)
            {
               if (static_cast<int>(unpaired.size()) < MAX_UNPAIRED_REPORT)
               {
                  UnpairedEntry e = {all_data[a*REC+0], all_data[a*REC+1],
                                     all_data[a*REC+2], rank_a,
                                     /*group_size=*/2};
                  unpaired.push_back(e);
               }
               n_unpaired += 2;  // both counted as anomalous
               i = j;
               continue;
            }
            for (int k = FIELD_BASE; k < FIELD_BASE + NUM_FIELDS; k++)
            {
               double va = all_data[a*REC + k];
               double vb = all_data[b*REC + k];
               double diff = std::abs(va - vb);
               // R-502: scale-relative tolerance.  Fields span ~1e-12
               // (slip rate at nucleation) to ~1e+8 Pa (normal stress);
               // a fixed absolute tol either fails on large-magnitude
               // fields (sub-ULP reassociation) or permits gross drift
               // on small ones.  `max(|va|, |vb|, 1.0)` gives a
               // physically meaningful scale that (a) never divides by
               // zero, (b) is O(1) for dimensionless fields like psi,
               // and (c) matches the field magnitude for large
               // dimensioned fields.
               double field_scale = std::max({std::abs(va), std::abs(vb),
                                              1.0});
               double rel_diff = diff / field_scale;
               if (rel_diff > max_rel_diff)
               {
                  max_rel_diff = rel_diff;
                  max_diff_abs = diff;
                  max_diff_scale = field_scale;
                  max_diff_field = k;
                  max_diff_entry = a;
               }
            }
            n_pairs++;
         }
         else
         {
            // group_size == 1 (singleton — only one rank claimed this QP),
            // or group_size >= 3 (centroid-match over-aggregated due to
            // transitive tolerance — less likely at the 1e-9 floor, but
            // still reported so the two failure modes are distinguishable
            // in the log).
            for (int m = i; m < j; m++)
            {
               int e = idx[m];
               if (static_cast<int>(unpaired.size()) < MAX_UNPAIRED_REPORT)
               {
                  UnpairedEntry u = {
                     all_data[e*REC+0], all_data[e*REC+1], all_data[e*REC+2],
                     static_cast<int>(all_data[e*REC + RANK_OFFSET]),
                     group_size
                  };
                  unpaired.push_back(u);
               }
            }
            n_unpaired += group_size;
         }
         i = j;
      }

      // Only rank 0 prints a summary; all ranks cooperate on the abort.
      // A shared fault QP must appear on exactly 2 ranks.  n_unpaired > 0
      // indicates a mesh-partitioning pathology (or a peer-rank resolution
      // bug that MFEM_VERIFY at the call site of Evaluate did not catch);
      // treat it as a hard diagnostic failure.  `tol` is interpreted as
      // a RELATIVE tolerance (R-502): a field with |field| ≈ 1e8 is
      // allowed to differ by ~tol × 1e8 across ranks.
      int fail_local = (max_rel_diff > tol || n_unpaired > 0) ? 1 : 0;
      int fail_global = 0;
      MPI_Allreduce(&fail_local, &fail_global, 1, MPI_INT, MPI_MAX, comm);

      if (fail_global)
      {
         if (n_unpaired > 0 && max_rel_diff <= tol)
         {
            // Rank 0 prints the detailed unpaired-entry report directly
            // to stderr with an explicit flush BEFORE calling MFEM_ABORT.
            // MFEM_ABORT → MPI_Abort can swallow long cerr-buffered
            // messages under SLURM/ibrun (job 7666151 saw only the
            // short non-rank-0 messages); fprintf+fflush is the
            // low-level escape hatch.  Every rank then calls a short
            // MFEM_ABORT for cooperative shutdown.
            if (my_rank_ == 0)
            {
               std::fprintf(stderr,
                  "\n[R-101 rank-0 detail]\n"
                  "  n_entries=%d, n_pairs=%d, n_unpaired=%d, "
                  "reported=%d%s\n",
                  n_entries, n_pairs, n_unpaired,
                  static_cast<int>(unpaired.size()),
                  (static_cast<int>(unpaired.size()) >= MAX_UNPAIRED_REPORT
                   ? " (truncated)" : ""));
               for (size_t u = 0; u < unpaired.size(); u++)
               {
                  // %.17e preserves full double precision so any
                  // sub-printable-digit FP drift that causes pairing
                  // failure is visible in the log (v8.0.0 Phase 1: job
                  // 7666171 printed "identical" coords at %.10e that
                  // actually differed by > 1e-9 and caused the failure).
                  std::fprintf(stderr,
                     "  [UNPAIRED] centroid=(%.17e, %.17e, %.17e) "
                     "rank=%d group_size=%d\n",
                     unpaired[u].cx, unpaired[u].cy, unpaired[u].cz,
                     unpaired[u].rank, unpaired[u].group_size);
               }
               std::fprintf(stderr,
                  "  Likely mesh-partitioning pathology (one rank "
                  "classifies a shared face as fault, peer does not; "
                  "global-vertex-key mismatch) OR tolerance collapse.\n\n");
               // v8.0.0 Phase 1 self-check (job 7666233 followup): the
               // abs_floor=1e-6 fix should have paired two entries whose
               // centroids agree to ~1 ULP.  If they still orphan, this
               // block distinguishes (a) same_centroid-logic failure vs
               // (b) sort-ordering pathology (a third entry sorts between
               // them and is NOT same_centroid with the first).
               if (unpaired.size() == 2)
               {
                  int idx_a = -1, idx_b = -1;
                  for (int e = 0; e < n_entries; e++)
                  {
                     int er = static_cast<int>(all_data[e*REC + RANK_OFFSET]);
                     double ex = all_data[e*REC + 0];
                     double ey = all_data[e*REC + 1];
                     double ez = all_data[e*REC + 2];
                     if (idx_a < 0 &&
                         ex == unpaired[0].cx && ey == unpaired[0].cy &&
                         ez == unpaired[0].cz && er == unpaired[0].rank)
                     { idx_a = e; }
                     if (idx_b < 0 &&
                         ex == unpaired[1].cx && ey == unpaired[1].cy &&
                         ez == unpaired[1].cz && er == unpaired[1].rank)
                     { idx_b = e; }
                  }
                  if (idx_a >= 0 && idx_b >= 0)
                  {
                     bool sc = same_centroid(idx_a, idx_b);
                     double dx = std::abs(all_data[idx_a*REC+0]
                                        - all_data[idx_b*REC+0]);
                     double dy = std::abs(all_data[idx_a*REC+1]
                                        - all_data[idx_b*REC+1]);
                     double dz = std::abs(all_data[idx_a*REC+2]
                                        - all_data[idx_b*REC+2]);
                     std::fprintf(stderr,
                        "  [self-check] same_centroid(orphan_0, orphan_1) "
                        "= %s\n"
                        "  [self-check]   dx=%.3e, dy=%.3e, dz=%.3e; "
                        "abs_floor=1.0e-6\n",
                        sc ? "TRUE" : "FALSE", dx, dy, dz);

                     int pos_a = -1, pos_b = -1;
                     for (int p = 0; p < n_entries; p++)
                     {
                        if (idx[p] == idx_a) { pos_a = p; }
                        if (idx[p] == idx_b) { pos_b = p; }
                     }
                     std::fprintf(stderr,
                        "  [self-check]   sort positions: orphan_0 at %d, "
                        "orphan_1 at %d (distance %d)\n",
                        pos_a, pos_b, std::abs(pos_a - pos_b));

                     // If non-adjacent, dump every entry between them so
                     // we can see what the "intruder" record(s) look like.
                     // Cap the window at 16 to keep the log bounded.
                     if (pos_a >= 0 && pos_b >= 0 &&
                         std::abs(pos_a - pos_b) != 1)
                     {
                        int lo = std::min(pos_a, pos_b);
                        int hi = std::max(pos_a, pos_b);
                        int span = hi - lo + 1;
                        int cap = std::min(span, 16);
                        for (int p = lo; p < lo + cap; p++)
                        {
                           int e = idx[p];
                           int er = static_cast<int>(
                              all_data[e*REC + RANK_OFFSET]);
                           std::fprintf(stderr,
                              "  [self-check]   pos=%d: "
                              "cx=%.17e cy=%.17e cz=%.17e rank=%d\n",
                              p,
                              all_data[e*REC+0],
                              all_data[e*REC+1],
                              all_data[e*REC+2], er);
                        }
                        if (span > cap)
                        {
                           std::fprintf(stderr,
                              "  [self-check]   ... (%d entries truncated)\n",
                              span - cap);
                        }
                     }
                  }
                  else
                  {
                     std::fprintf(stderr,
                        "  [self-check] unable to re-locate orphan data "
                        "indices (idx_a=%d, idx_b=%d)\n", idx_a, idx_b);
                  }
               }
               std::fflush(stderr);
            }
            MFEM_ABORT("R-101 shared-fault DOFData: " << n_unpaired
                       << " unpaired entries (every shared QP should have "
                       "exactly 2 ranks).  See rank-0 [R-101 rank-0 detail] "
                       "lines above this abort trace.");
         }
         double cx = all_data[max_diff_entry*REC + 0];
         double cy = all_data[max_diff_entry*REC + 1];
         double cz = all_data[max_diff_entry*REC + 2];
         const int field_off = max_diff_field - FIELD_BASE;
         const char *field_name =
            (field_off >= 0 && field_off < NUM_FIELDS)
               ? FIELD_NAMES[field_off] : "<unknown>";
         // R-504: generic diagnostic wording.  Prior to R-501 this blamed
         // "R-001's (+,-) canonicalisation"; with the owner-broadcast fix
         // in place the underlying R-001 path no longer runs, so any
         // future trip is from a different bug (missed field in the
         // broadcast, R-501 MPI exchange error, stale auth_state, ...).
         MFEM_ABORT("R-101 shared-fault DOFData consistency FAILED.  "
                    "Field '" << field_name << "' at centroid ("
                    << cx << ", " << cy << ", " << cz
                    << ") differs by " << max_diff_abs
                    << " across the two ranks sharing the face (scale="
                    << max_diff_scale << ", rel_diff=" << max_rel_diff
                    << ", rel_tol=" << tol << ").  The two ranks' DOFData "
                    "diverged — check that the R-501 owner-broadcast in "
                    "ComputeSharedFaceFluxRHS covers every mutable field "
                    "written by FaultFaceFlux::Evaluate and the driver's "
                    "RK4 averaging step, and that the Q_imp exchange is "
                    "packed/unpacked in a consistent order on both ranks.");
      }

      if (my_rank_ == 0)
      {
         std::cout << "  [R-101 check] shared-fault DOFData consistency OK: "
                   << n_pairs << " pairs matched, " << n_unpaired
                   << " unpaired entries, max_rel_diff=" << max_rel_diff
                   << " (max_abs_diff=" << max_diff_abs
                   << ", scale=" << max_diff_scale
                   << ", rel_tol=" << tol << ")\n";
      }
#else
      (void)tol;
#endif
   }
}

// ---------------------------------------------------------------------------
// Classify boundary face by attribute
// ---------------------------------------------------------------------------
template <typename MeshType>
typename WaveOperator<MeshType>::FaceBC
WaveOperator<MeshType>::ClassifyBoundaryFace(int bdr_attr) const
{
   if (bc_.absorbing_attrs.count(bdr_attr)) { return FaceBC::Absorbing; }
   if (bc_.natural_attrs.count(bdr_attr))   { return FaceBC::FreeSurface; }
   if (bdr_attr == bc_.fault_attr)          { return FaceBC::Fault; }
   return FaceBC::Absorbing;
}
