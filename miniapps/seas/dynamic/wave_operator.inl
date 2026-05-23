// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// WaveOperator template implementation (included from wave_operator.hpp).
// NOTE: seas_diag_rank.hpp must be included at file scope — include it from
// wave_operator.hpp before the `namespace mfem::seas` block opens, NOT here
// (wave_operator.hpp #includes this .inl from inside `namespace mfem::seas`,
// so any `#include` here would be nested).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
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
         // R-1510: STATIC-MESH assumption.  These two exchanges build the
         // face-neighbour topology graph used by every per-step ghost
         // exchange downstream (Mult, AdvanceADER, ComputeADERSharedFaceFluxRHS,
         // EvaluateBulkAtFaultQPsCanonical).  Re-calling them is unnecessary
         // for a static partition.  An AMR / moving-mesh extension MUST
         // re-call BOTH (in this order: ParMesh first, then
         // ParFiniteElementSpace) before any subsequent macro-step or the
         // per-step `q_gf.ExchangeFaceNbrData()` calls will exchange values
         // through a stale topology.  TPV104 production uses a static mesh
         // so the once-only call is correct here.
         pmesh.ExchangeFaceNbrData();
         pfes->ExchangeFaceNbrData();
         ghost_gf_ = std::make_unique<ParGridFunction>(pfes);
         ghost_initialized_ = true;

         // R-1601: build the vdim=NUM_STATE byNODES batched ghost exchange
         // PFE space + PGF used by EvaluateBulkAtFaultQPsCanonical to
         // amortise the 9·O per-substep collectives down to O.  Static-
         // mesh assumption (R-1510) means we exchange topology once here
         // and reuse it for every macro-step.  byNODES ordering matches
         // the wave operator's component-major Q layout
         // (`Q[c * ndof_total + i]`), so packing/unpacking are simple
         // strided copies (no per-DOF stride conversion).
         pfes_full_state_ = std::make_unique<ParFiniteElementSpace>(
            &pmesh, fec_.get(), NUM_STATE, Ordering::byNODES);
         pfes_full_state_->ExchangeFaceNbrData();
         ghost_gf_full_state_ = std::make_unique<ParGridFunction>(
            pfes_full_state_.get());
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

      // R-101: per-interior-fault-face elem1_on_plus computed from
      // geometry (element centroid vs face centroid projected onto
      // ref_normal), identical recipe to the shared-fault block below.
      // Replaces the per-QP `!qpd.sign_flipped` derivation that was
      // FP-sensitive on near-axis-aligned faces (the per-QP CalcOrtho
      // normal can flip sign across QPs of the same face when the
      // dot product with ref_normal is near zero, producing per-QP
      // bimodal flux orientation and breaking y-mirror invariance).
      // The per-face flag is FP-stable: the centroid-projection margin
      // is O(dx/2), far above any FP noise.
      interior_fault_elem1_on_plus_.assign(fault_interior_faces_.Size(),
                                           false);
      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         int f = fault_interior_faces_[i];
         FaceElementTransformations *ftr =
            mesh_.GetInteriorFaceTransformations(f);
         if (!ftr) { continue; }

         // Face centroid.
         const IntegrationPoint &ip_center =
            Geometries.GetCenter(ftr->GetGeometryType());
         ftr->Face->SetIntPoint(&ip_center);
         Vector face_c(3);
         ftr->Face->Transform(ip_center, face_c);

         // Elem1 centroid (vertex-average; sufficient for tet/hex).
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

         // R-SAFS (curvilinear-fault side fix): project the centroid offset
         // onto the CANONICAL FACE NORMAL, not the hardcoded ref_normal.
         // Projecting onto ref_normal mislabels the +/- side once the fault
         // tilts enough that ref_normal's tangential component dominates the
         // (asymmetric) centroid offset — the margin |elem1_proj - face_proj|
         // is O(|n·ref_normal|) and flips sign at a geometry-set tilt angle,
         // breaking the "exactly one side is +" invariant on shared faces
         // (tilted-fault serial-vs-parallel test).  The true face normal gives
         // a robust O(element-size) margin at EVERY fault-trace angle.  For
         // planar TPV/BP5 faults the canonical normal == ref_normal (up to a
         // positive |J_F| scale), so the comparison sign — hence the stored
         // bool — is unchanged (byte-exact).  Convention preserved: a point
         // with SMALLER projection onto the canonical normal is on the "+"
         // side.
         Vector fn(3);
         CalcOrtho(ftr->Face->Jacobian(), fn);
         if (FaultBasis::NormalNeedsFlipToCanonical(fn, ref_normal, 3))
         {
            fn.Neg();
         }
         real_t face_proj = 0.0, elem1_proj = 0.0;
         for (int d = 0; d < 3; d++)
         {
            face_proj  += face_c(d)  * fn(d);
            elem1_proj += elem1_c(d) * fn(d);
         }
         interior_fault_elem1_on_plus_[i] = (elem1_proj < face_proj);
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

            // R-SAFS (curvilinear-fault side fix): project the centroid offset
            // onto the CANONICAL FACE NORMAL, not the hardcoded ref_normal —
            // see the interior-fault block above for the full rationale.  On a
            // shared face CalcOrtho gives OPPOSITE raw normals on the two
            // ranks (±n_phys); NormalNeedsFlipToCanonical maps both to the SAME
            // canonical direction (a pure function of the physical normal
            // line), so projecting each rank's local elem1 centroid offset
            // onto it yields exactly ONE rank with elem1_on_plus == true at
            // every fault-trace angle (the invariant the old ref_normal
            // projection broke past a geometry-set tilt).  Convention
            // preserved: SMALLER projection onto the canonical normal => "+"
            // side.  Planar TPV/BP5: canonical normal == ref_normal => bool
            // unchanged (byte-exact).
            Vector fn(3);
            CalcOrtho(ftr->Face->Jacobian(), fn);
            if (FaultBasis::NormalNeedsFlipToCanonical(fn, ref_normal, 3))
            {
               fn.Neg();
            }
            real_t face_proj  = 0.0, elem1_proj = 0.0;
            for (int d = 0; d < 3; d++)
            {
               face_proj  += face_c(d)  * fn(d);
               elem1_proj += elem1_c(d) * fn(d);
            }
            shared_fault_elem1_on_plus_[sf_idx] = (elem1_proj < face_proj);

#ifdef SEAS_DIAG_CENTROID_MARGIN
            // v9.2.0 §17 (rev-3d, post-REVIEW R-V92-C04): dump the
            // per-face FP margin of the `elem1_on_plus` comparison so
            // we can quantify the production-mesh FP-fragility risk
            // without running any time stepping.  Margin = |elem1_proj
            // − face_proj|; a small margin means the sign of the
            // comparison is fragile (centroid close to the fault plane
            // in the ref_normal direction).  On the 4-km Cartesian
            // fixture every margin is O(10³ m); on Gmsh-duplicated
            // production meshes near x = 0 the margin can be O(ε_FP).
            const real_t margin = std::abs(elem1_proj - face_proj);
            std::fprintf(stderr,
               "[CENTROID-MARGIN] sf_idx=%d elem1_proj=%.17e face_proj=%.17e "
               "margin=%.6e elem1_on_plus=%d face_c=(%.3e,%.3e,%.3e)\n",
               sf_idx, static_cast<double>(elem1_proj),
               static_cast<double>(face_proj),
               static_cast<double>(margin),
               shared_fault_elem1_on_plus_[sf_idx] ? 1 : 0,
               static_cast<double>(face_c(0)),
               static_cast<double>(face_c(1)),
               static_cast<double>(face_c(2)));
            std::fflush(stderr);
#endif
         }
#endif
      }
   }
}

template <typename MeshType>
WaveOperator<MeshType>::~WaveOperator() = default;

// ---------------------------------------------------------------------------
// Phase H.1 (Stage 1) — heterogeneous-material ctor.
//
// Delegates to the scalar ctor above with the constants extracted from
// the MaterialField, so the resulting `flux_` member is byte-identical
// to the scalar-ctor path on Mode::Constant input.  Then builds the
// per-element flux pool (H.1), the per-element CFL length cache
// (H.3 input), and the bi-material shared-face neighbour map stub
// (H.4 Stage 1 in-place fill; real cross-rank exchange in Stage 2).
//
// Mode::Coefficient input is ABORTED unconditionally in Stage 1: the
// per-element flux dispatch (H.2) at every flux_.X() call site in
// wave_operator.inl is the Stage 2 deliverable.  Until then the
// scalar `flux_` is still consulted by every Mult/AdvanceADER call;
// dispatching it on a per-element-Eval seed would silently produce
// wrong physics on heterogeneous input.  The abort here is the only
// gate that prevents that.
//
// In Mode::Constant the new ctor's observable output is BYTE-IDENTICAL
// to the scalar ctor's, proven by
// `test_phaseh_wave_operator_constant_parity` (np=1 + np=4).
// ---------------------------------------------------------------------------
template <typename MeshType>
WaveOperator<MeshType>::WaveOperator(MeshType &mesh, int order,
                                     const MaterialField &material,
                                     const BoundaryConfig &bc)
   : WaveOperator(mesh, order,
                  material.mode == MaterialField::Mode::Constant
                     ? material.lambda_const : real_t(1.0),
                  material.mode == MaterialField::Mode::Constant
                     ? material.mu_const     : real_t(1.0),
                  material.mode == MaterialField::Mode::Constant
                     ? material.rho_const    : real_t(1.0),
                  bc)
{
   MFEM_VERIFY(material.mode == MaterialField::Mode::Constant
               || material.mode == MaterialField::Mode::Coefficient,
               "WaveOperator(MaterialField): material.mode must be "
               "Constant or Coefficient; got Mode::GridFunction ("
               << static_cast<int>(material.mode) << ").  GridFunction "
               "mode is consumed via At(elem, dof, ...), which has no "
               "well-defined centroid evaluation needed by the per-"
               "element flux pool.");

   MFEM_VERIFY(material.mode == MaterialField::Mode::Constant,
               "WaveOperator(MaterialField) Phase H Stage 1 gap: "
               "Mode::Coefficient input requires the per-element flux "
               "dispatch in wave_operator.inl (plan Phase H.2), which "
               "is NOT yet wired in this commit.  Constructing on a "
               "Coefficient material would silently fall back to the "
               "scalar `flux_` member built from the placeholder seed "
               "(1.0, 1.0, 1.0) inside the delegating ctor above — "
               "i.e., wrong physics on every Mult call.  Re-run with "
               "MaterialField::MakeConstant(...) until Stage 2 lands.");

   material_ = &material;

   BuildGodunovFluxPool_(material);
   ExchangeBiMaterialNeighbours_();
}

// ---------------------------------------------------------------------------
// Phase H.1 helper — fill per_elem_lmr_ / per_elem_h_, then Build()
// the owned flux pool.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::BuildGodunovFluxPool_(
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
// Stage 1 contract:
//   - Mode::Constant input is the only reachable path (the ctor
//     MFEM_VERIFY above aborts Coefficient mode).
//   - In Mode::Constant every neighbour's per-element material is, by
//     construction, equal to the local material.  We therefore fill
//     the map with the LOCAL-side (lambda, mu, rho) for every shared
//     face on this rank — no MPI is performed.
//   - On serial Mesh this is a no-op (no shared faces).
//
// Stage 2 will replace this body with the real MPI_Allgatherv exchange.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ExchangeBiMaterialNeighbours_()
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
   }
#endif
}

// ---------------------------------------------------------------------------
// TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23
// Phase 2a (§6.3 + §6.3a): opt-in precomputed-flux switch with late
// `fault_face_set_` population (R4-001 + R5-002 FIXES).
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::UsePrecomputedFaceFluxes(bool enable)
{
   // R-1203 cross-check: precomputed-flux and mixed-flux are mutually
   // exclusive.  The precomputed flux cache bakes upwind into its tables;
   // mixed-flux at any non-None mode would be silently ignored if the
   // precomputed-flux dispatch fires first.  Abort if the OTHER setter
   // is in a non-default state.  Symmetric guard in SetMixedFluxMode.
   if (enable && mixed_flux_mode_ != MixedFluxMode::None)
   {
      MFEM_ABORT("UsePrecomputedFaceFluxes(true): mutually exclusive "
                 "with mixed-flux mode (currently "
                 << static_cast<int>(mixed_flux_mode_)
                 << ").  Disable mixed flux first via "
                 "SetMixedFluxMode(MixedFluxMode::None).");
   }

   if (enable)
   {
      // === §6.3a (R5-002 FIX) — late population of fault_face_set_ ===
      //
      // fault_interior_faces_ and fault_shared_faces_ are populated by
      // the WaveOperator ctor (see wave_operator.inl construction block
      // above: `fault_interior_faces_.Append(face_idx)` loop and the
      // parallel `fault_shared_faces_.Append(sf)` loop).  Populating
      // fault_face_set_ HERE — not in the ctor — is robust against
      // future refactors that might move the fault-face array population
      // to a post-ctor helper, and it keeps the initialization at one
      // unambiguous call site.
      //
      // R-002 (v9.5.0): empty local fault arrays are legitimate on MPI
      // ranks whose partition contains no fault-adjacent elements.  The
      // former MFEM_VERIFY(fault_interior_faces_ + fault_shared_faces_
      // > 0) aborted those ranks on any call to UsePrecomputedFaceFluxes
      // (true), even though the condition is a normal consequence of
      // domain decomposition on a fault that touches only a minority of
      // ranks.  The remaining logic below is already a no-op for empty
      // arrays (the fault_face_set_ stays empty), so we simply drop the
      // guard.  A globally-empty-fault / bc_.fault_attr > 0 mis-config
      // is not detectable here without an MPI_Allreduce; that check
      // belongs to the mesh/config layer.
      if (bc_.fault_attr > 0)
      {
         fault_face_set_.clear();
         for (int i = 0; i < fault_interior_faces_.Size(); i++)
         {
            fault_face_set_.insert(fault_interior_faces_[i]);
         }
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            // fault_shared_faces_ stores shared-face indices (sf), not
            // mesh face indices; we need the mesh face index for the
            // set.  ParMesh::GetSharedFace translates sf → mesh face.
            if constexpr (IsParallelMesh<MeshType>::value)
            {
#ifdef MFEM_USE_MPI
               auto &pmesh = static_cast<const ParMesh &>(mesh_);
               int mesh_face_idx = const_cast<ParMesh &>(pmesh).GetSharedFace(
                  fault_shared_faces_[i]);
               fault_face_set_.insert(mesh_face_idx);
#endif
            }
         }
      }
      else
      {
         fault_face_set_.clear();
      }
      // === end §6.3a population =====================================

      if (!precomputed_face_fluxes_.IsInitialized())
      {
         // R-002 (v9.5.0): the companion fault_face_set_ emptiness guard
         // is the same class of false-positive on fault-free MPI ranks.
         // Empty set on a rank whose partition contains no fault slice
         // is legitimate; only the face_bdr_attr_ sizing check remains
         // as a meaningful ctor-bookkeeping tripwire.
         MFEM_VERIFY(face_bdr_attr_.size() ==
                     static_cast<size_t>(mesh_.GetNumFaces()),
                     "WaveOperator::UsePrecomputedFaceFluxes: "
                     "face_bdr_attr_ not sized to mesh.GetNumFaces(); "
                     "WaveOperator ctor bookkeeping is incomplete.");

         precomputed_face_fluxes_.Init(
            mesh_,                                         // Mesh &     (R3-008)
            *fes_,                                         // const FES&
            bc_,                                           // const BoundaryConfig&
            flux_,                                         // const GodunovFlux&
            free_surface_bc_mode_,                         // FreeSurfaceBCMode
            fault_face_set_,                               // R4-002 / R5-002
            face_bdr_attr_,                                // R3-003
            shared_mesh_face_set_,                         // R3-002
            arm3d_tables_enabled_);                        // R-003 v9.5.0 (Arm 3d gate)

         // Phase 2b (§7.2): extend Init with a ParMesh-aware shared-face
         // pass.  One FaceEntry per shared non-fault face per rank, keyed
         // by (face_idx, Elem1No).  No-op on serial builds.
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            auto *pmesh_ptr = dynamic_cast<ParMesh *>(&mesh_);
            if (pmesh_ptr)
            {
               precomputed_face_fluxes_.InitSharedFaces(
                  *pmesh_ptr,
                  *fes_,
                  flux_,
                  fault_face_set_);
            }
#endif
         }
      }
   }
   use_precomputed_face_fluxes_ = enable;
}

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

#ifdef SEAS_DIAG_FAULT_FLUX
   // C-2A BULK-DOF: per-Mult-call dump of bulk Q (input) and dQdt
   // (pre-mass-inverse RHS) at the hypocenter face's two adjacent tets.
   // Reads the two driver-tagged DOFs (one per side) directly from the
   // bulk Q vector — no face-QP interpolation, no fault-Riemann logic.
   // Single-rank stderr printf (other ranks have diag_elem_*_=-1).
   if (diag_elem_plus_ >= 0 && diag_elem_minus_ >= 0 &&
       diag_face_dof_plus_ >= 0 && diag_face_dof_minus_ >= 0)
   {
      const real_t *Qd = Q.GetData();
      const real_t *Rd = dQdt.GetData();
      const int dof_off_p = diag_elem_plus_  * ndof_per_el_;
      const int dof_off_m = diag_elem_minus_ * ndof_per_el_;
      const int idx_p = diag_face_dof_plus_;
      const int idx_m = diag_face_dof_minus_;

      // Bounds check: silent on bogus indices (defensive).
      if (idx_p < ndof_per_el_ && idx_m < ndof_per_el_)
      {
         real_t Qp[NUM_STATE], Qm[NUM_STATE], Rp[NUM_STATE], Rm[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Qp[c] = Qd[c*ndof_total_ + dof_off_p + idx_p];
            Qm[c] = Qd[c*ndof_total_ + dof_off_m + idx_m];
            Rp[c] = Rd[c*ndof_total_ + dof_off_p + idx_p];
            Rm[c] = Rd[c*ndof_total_ + dof_off_m + idx_m];
         }
         std::fprintf(stderr,
            "[C-2A BULK-DOF Q+]  rank=%d e=%d  "
            "SXX=%+.4e SYY=%+.4e SZZ=%+.4e SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
            "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
            g_seas_my_rank, diag_elem_plus_,
            Qp[SXX], Qp[SYY], Qp[SZZ], Qp[SXY], Qp[SYZ], Qp[SXZ],
            Qp[VX],  Qp[VY],  Qp[VZ]);
         std::fprintf(stderr,
            "[C-2A BULK-DOF Q-]  rank=%d e=%d  "
            "SXX=%+.4e SYY=%+.4e SZZ=%+.4e SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
            "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
            g_seas_my_rank, diag_elem_minus_,
            Qm[SXX], Qm[SYY], Qm[SZZ], Qm[SXY], Qm[SYZ], Qm[SXZ],
            Qm[VX],  Qm[VY],  Qm[VZ]);
         std::fprintf(stderr,
            "[C-2A BULK-DOF k+]  rank=%d e=%d  "
            "SXX=%+.4e SYY=%+.4e SZZ=%+.4e SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
            "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
            g_seas_my_rank, diag_elem_plus_,
            Rp[SXX], Rp[SYY], Rp[SZZ], Rp[SXY], Rp[SYZ], Rp[SXZ],
            Rp[VX],  Rp[VY],  Rp[VZ]);
         std::fprintf(stderr,
            "[C-2A BULK-DOF k-]  rank=%d e=%d  "
            "SXX=%+.4e SYY=%+.4e SZZ=%+.4e SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
            "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
            g_seas_my_rank, diag_elem_minus_,
            Rm[SXX], Rm[SYY], Rm[SZZ], Rm[SXY], Rm[SYZ], Rm[SXZ],
            Rm[VX],  Rm[VY],  Rm[VZ]);
      }
   }
#endif

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
// ApplySpatialDerivative: element-local L2-projected ∂_{x_dir} Q (ADER Phase 1)
// ---------------------------------------------------------------------------
// Computes dQ_dxdir[c,i] = (M_e^{-1} K_d^e Q_c)[i] per element, per component,
// with NO flux coupling.  Matches ComputeVolumeRHS's quadrature rule
// (IntRules.Get(geom, 2*order_)) so the projection is consistent with the
// volume integral semantics.  Building block for the Cauchy-Kovalevskaya
// recursion (Phase 3).
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ApplySpatialDerivative(int dir,
                                                    const Vector &Q,
                                                    Vector &dQ_dxdir) const
{
   MFEM_VERIFY(dir >= 0 && dir < 3,
               "ApplySpatialDerivative: dir must be in {0,1,2}, got " << dir);
   MFEM_VERIFY(Q.Size() == NUM_STATE * ndof_total_,
               "ApplySpatialDerivative: Q size "
               << Q.Size() << " != NUM_STATE * ndof_total_ = "
               << NUM_STATE * ndof_total_);

   dQ_dxdir.SetSize(NUM_STATE * ndof_total_);
   dQ_dxdir = 0.0;

   if (ne_ == 0) { return; }

   const real_t *Q_data = Q.GetData();
   real_t *dQ_data = dQ_dxdir.GetData();

   // Scratch: K_d^e * Q_c accumulator, shape [NUM_STATE, ndof].  Reused per
   // element.  Size follows ndof_per_el_ which is constant across elements
   // for a homogeneous-order L2 space.
   DenseMatrix KdQ(NUM_STATE, ndof_per_el_);

   for (int e = 0; e < ne_; e++)
   {
      const FiniteElement *fe = fes_->GetFE(e);
      ElementTransformation *Tr = fes_->GetElementTransformation(e);
      int ndof = fe->GetDof();

      // R-001 guard: this routine (and the whole WaveOperator) assumes a
      // homogeneous DG space — every element has the same `GetDof()`.  The
      // element loop writes into `KdQ(c, i < ndof)` which is sized to
      // `ndof_per_el_`; a heterogeneous mesh (mixed element types or
      // variable-order L2) would silently overflow that row.  Fail loud so
      // a future extension cannot introduce that corruption without
      // updating the scratch sizing.
      // R-1508: promoted from MFEM_ASSERT (Debug-only) to MFEM_VERIFY so
      // a Release build still fails loud rather than silently overflowing
      // the per-element scratch row on a heterogeneous mesh.
      MFEM_VERIFY(ndof == ndof_per_el_,
                  "ApplySpatialDerivative assumes homogeneous elements: "
                  "elem=" << e << " has ndof=" << ndof
                  << " but ndof_per_el_=" << ndof_per_el_);

      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2*order_);
      int nqp = ir.GetNPoints();

      int dof_offset = e * ndof_per_el_;

      Vector shape(ndof);
      DenseMatrix dshape(ndof, 3);

      KdQ = 0.0;

      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr->SetIntPoint(&ip);
         real_t w = ip.weight * Tr->Weight();

         fe->CalcShape(ip, shape);
         fe->CalcPhysDShape(*Tr, dshape);

         // Evaluate ∂_dir Q_c at this QP: dQdir_qp[c] = Σ_j dshape(j,dir) * Q_c(j)
         real_t dQdir_qp[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t s = 0.0;
            const real_t *Qc = Q_data + c * ndof_total_ + dof_offset;
            for (int j = 0; j < ndof; j++)
            {
               s += dshape(j, dir) * Qc[j];
            }
            dQdir_qp[c] = s;
         }

         // Accumulate (K_d Q_c)[i] += w * φ_i(q) * ∂_dir Q_c(q)
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t wd = w * dQdir_qp[c];
            for (int i = 0; i < ndof; i++)
            {
               KdQ(c, i) += shape(i) * wd;
            }
         }
      }

      // Apply element mass inverse: dQ_c[i] = Σ_j M_e^{-1}(i,j) * (K_d Q_c)[j]
      const DenseMatrix &M_inv = elem_mass_inv_[e];
      for (int c = 0; c < NUM_STATE; c++)
      {
         real_t *dQc = dQ_data + c * ndof_total_ + dof_offset;
         for (int i = 0; i < ndof; i++)
         {
            real_t val = 0.0;
            for (int j = 0; j < ndof; j++)
            {
               val += M_inv(i, j) * KdQ(c, j);
            }
            dQc[i] = val;
         }
      }
   }
}

// ---------------------------------------------------------------------------
// ADER I-05 Phase 3: helper — per-DOF 9x9 Jacobian application with sign.
// ---------------------------------------------------------------------------
// Computes Y[c, dof] += sign * Σ_{c'} A(c, c') * X[c', dof] for every DOF.
// Used to implement L(D) = -Σ_d A_d · ∂_d D on the time-integrated state.
// Layout: X/Y are Vectors of size NUM_STATE * ndof_total, component-major
// (X[c * ndof_total + dof]).
//
// R-004: `static inline` (file-local linkage) is preferred to an
// anonymous namespace here because this .inl is #include'd from the
// header — an anonymous namespace would otherwise emit a separate
// internal-linkage copy in every TU that transitively sees wave_operator.hpp.
// `static inline` expresses the same intent more idiomatically.
// ---------------------------------------------------------------------------
static inline void ApplyJacobianPerDOF(const DenseMatrix &A,
                                       const Vector &X, Vector &Y,
                                       int ndof_total, real_t sign)
{
   MFEM_ASSERT(A.Height() == NUM_STATE && A.Width() == NUM_STATE,
               "ApplyJacobianPerDOF: A must be NUM_STATE x NUM_STATE");
   MFEM_ASSERT(X.Size() == NUM_STATE * ndof_total,
               "ApplyJacobianPerDOF: X size mismatch");
   MFEM_ASSERT(Y.Size() == NUM_STATE * ndof_total,
               "ApplyJacobianPerDOF: Y size mismatch");

   const real_t *Xd = X.GetData();
   real_t *Yd = Y.GetData();

   // Row-major walk over components.  For each row c, accumulate
   // Σ_c' A(c,c') * X[c', dof] into Y[c, dof] at every DOF.
   for (int c = 0; c < NUM_STATE; c++)
   {
      real_t *Yc = Yd + c * ndof_total;
      for (int cp = 0; cp < NUM_STATE; cp++)
      {
         const real_t a = A(c, cp);
         if (a == 0.0) { continue; }
         const real_t w = sign * a;
         const real_t *Xcp = Xd + cp * ndof_total;
         for (int i = 0; i < ndof_total; i++)
         {
            Yc[i] += w * Xcp[i];
         }
      }
   }
}

// ---------------------------------------------------------------------------
// ADER I-05 Phase 3: Cauchy-Kovalevskaya time-integrated state predictor.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ComputeADERTimeIntegrated(
   const Vector &Q, real_t dt, int order, Vector &I) const
{
   MFEM_VERIFY(order >= 2 && order <= 4,
               "ComputeADERTimeIntegrated: order must be in {2,3,4}, got "
               << order);
   MFEM_VERIFY(Q.Size() == NUM_STATE * ndof_total_,
               "ComputeADERTimeIntegrated: Q size "
               << Q.Size() << " != NUM_STATE * ndof_total_ = "
               << NUM_STATE * ndof_total_);
   // R-002 guard: if the caller aliases Q and I, the I = 0.0 below would
   // zero the input before D_curr is copied, silently yielding I = 0.
   MFEM_VERIFY(&Q != &I,
               "ComputeADERTimeIntegrated: Q and I must be distinct Vectors "
               "(aliasing would zero the input before the CK copy)");

   I.SetSize(NUM_STATE * ndof_total_);
   I = 0.0;

   if (dt <= 0.0 || ndof_total_ == 0) { return; }

   // R-1501: lazy-init mutable scratch buffers (ck_D_curr_buf_,
   // ck_D_next_buf_, ck_dQ_dxd_buf_) reused across macro-steps.
   // See header docstring for the cost rationale.  SetSize is a no-op
   // when the size already matches; the first call sizes once.  Note:
   // `D_curr = Q` (assignment) writes a fresh copy of Q into the
   // member buffer on every call, so prior-call residual values are
   // overwritten before the first read.  Same applies to `D_next` /
   // `dQ_dxd`: D_next is zeroed at L1010 before each iteration's
   // accumulation, and dQ_dxd is fully overwritten by
   // ApplySpatialDerivative at L1013 before any read.
   const int N = NUM_STATE * ndof_total_;
   if (ck_D_curr_buf_.Size() != N) { ck_D_curr_buf_.SetSize(N); }
   if (ck_D_next_buf_.Size() != N) { ck_D_next_buf_.SetSize(N); }
   if (ck_dQ_dxd_buf_.Size() != N) { ck_dQ_dxd_buf_.SetSize(N); }
   Vector &D_curr = ck_D_curr_buf_;
   Vector &D_next = ck_D_next_buf_;
   Vector &dQ_dxd = ck_dQ_dxd_buf_;
   D_curr = Q;                                   // D(0) = Q

   // k = 0 contribution: I += dt * D(0)
   real_t fac = dt;
   I.Add(fac, D_curr);

   for (int k = 0; k < order - 1; k++)
   {
      // D_next = L(D_curr) = -Σ_d A_d · ∂_{x_d} D_curr
      D_next = 0.0;
      for (int d = 0; d < 3; d++)
      {
         ApplySpatialDerivative(d, D_curr, dQ_dxd);
         const DenseMatrix &A_d = flux_.GetReferenceStarMatrix(d);
         ApplyJacobianPerDOF(A_d, dQ_dxd, D_next, ndof_total_, /*sign=*/-1.0);
      }

      // Advance factorial factor: fac *= dt / (k+2).
      //
      // R-1503 cross-reference: `ComputeADERSubStepStates` below uses
      // `denom = k + 1` for the analogous CK update; the two parameterise
      // the recursion at different phases:
      //   - Here `fac` is initialised to `dt` at L1004 (so before the loop,
      //     fac == dt^{k=0+1}/0! · 1 represents the k=0 term coefficient
      //     dt^{k+1}/(k+1)!).  At loop iteration k we are computing D(k+1),
      //     which contributes dt^{k+2}/(k+2)! · D(k+1) to I; the update is
      //     therefore `fac *= dt / (k+2)`.
      //   - `ComputeADERSubStepStates` initialises `fac[o] = 1.0` at L1093
      //     (its k=0 term `(τ^0 / 0!) · D(0) = D(0)` is pre-added at
      //     L1086-1089 with no scaling).  At its iteration k it computes
      //     D(k+1) and accumulates `(τ^{k+1}/(k+1)!) · D(k+1)`, so its
      //     update is `fac[o] *= τ / (k+1)`.
      // Both are mathematically equivalent (different k-indexing, same
      // Taylor coefficient).  DO NOT "harmonise" by changing one to match
      // the other without re-deriving the recursion phase — see R-1503.
      fac *= dt / static_cast<real_t>(k + 2);
      I.Add(fac, D_next);

      // Swap buffers: D_curr <- D_next for the next iteration.
      mfem::Swap(D_curr, D_next);
   }
}

// ---------------------------------------------------------------------------
// R-602: ComputeADERSubStepStates — pointwise Q at sub-step nodes via
//        Cauchy-Kovalevskaya Taylor expansion.
//
// Same recursion D(k+1) = -Σ_d A_d · ∂_{x_d} D(k) as
// ComputeADERTimeIntegrated; instead of integrating with weight
// dt^{k+1}/(k+1)!, this evaluates the Taylor polynomial
//   Q(τ) = Σ_{k=0}^{O-1} (τ^k / k!) · D(k)
// at each supplied τ in tau_nodes.
//
// Used by Tpv104SubStepIterator::AdvanceWithSubStepStates (and the driver's
// AdvanceADERWithSubStep helper) to match SeisSol's per-substep
// `qInterpolated[o]` semantics.  The Gauss-Lobatto nodes on [0, dt] are
// the standard ADER-O sub-step quadrature.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ComputeADERSubStepStates(
   const Vector &Q, real_t dt, int order,
   const std::vector<real_t> &tau_nodes,
   std::vector<Vector> &Q_per_node) const
{
   MFEM_VERIFY(order >= 2 && order <= 4,
               "ComputeADERSubStepStates: order must be in {2,3,4}, got "
               << order);
   MFEM_VERIFY(dt > 0.0,
               "ComputeADERSubStepStates: dt must be > 0, got " << dt);
   MFEM_VERIFY(Q.Size() == NUM_STATE * ndof_total_,
               "ComputeADERSubStepStates: Q size " << Q.Size()
               << " != NUM_STATE * ndof_total_ = "
               << NUM_STATE * ndof_total_);
   const int O_nodes = static_cast<int>(tau_nodes.size());
   MFEM_VERIFY(O_nodes >= 1,
               "ComputeADERSubStepStates: tau_nodes must be non-empty");
   for (int o = 0; o < O_nodes; o++)
   {
      MFEM_VERIFY(std::isfinite(tau_nodes[o]),
                  "ComputeADERSubStepStates: tau_nodes[" << o
                  << "] must be finite, got " << tau_nodes[o]);
      MFEM_VERIFY(tau_nodes[o] >= 0.0 && tau_nodes[o] <= dt,
                  "ComputeADERSubStepStates: tau_nodes[" << o
                  << "] = " << tau_nodes[o]
                  << " is out of [0, dt] = [0, " << dt << "]");
   }

   Q_per_node.resize(O_nodes);
   for (int o = 0; o < O_nodes; o++)
   {
      Q_per_node[o].SetSize(NUM_STATE * ndof_total_);
      Q_per_node[o] = 0.0;
   }

   if (ndof_total_ == 0) { return; }

   // R-1501: lazy-init separate mutable scratch buffers for the substep
   // predictor (`ck_substep_*_buf_`) — distinct from `ck_*_buf_` used by
   // ComputeADERTimeIntegrated so the two predictors are independent.
   // SetSize is a no-op on resized-once buffers; D_curr is overwritten by
   // `= Q`, D_next is zeroed at L1156, dQ_dxd is fully overwritten by
   // ApplySpatialDerivative — no stale-read hazards.
   const int N = NUM_STATE * ndof_total_;
   if (ck_substep_D_curr_buf_.Size() != N)
   { ck_substep_D_curr_buf_.SetSize(N); }
   if (ck_substep_D_next_buf_.Size() != N)
   { ck_substep_D_next_buf_.SetSize(N); }
   if (ck_substep_dQ_dxd_buf_.Size() != N)
   { ck_substep_dQ_dxd_buf_.SetSize(N); }
   Vector &D_curr = ck_substep_D_curr_buf_;
   Vector &D_next = ck_substep_D_next_buf_;
   Vector &dQ_dxd = ck_substep_dQ_dxd_buf_;
   D_curr = Q;                                    // D(0) = Q

   // k = 0 contribution: Q(τ) += (τ^0 / 0!) * D(0) = D(0).
   for (int o = 0; o < O_nodes; o++)
   {
      Q_per_node[o].Add(1.0, D_curr);
   }

   // Recursion factor f(k) = τ^k / k!, evolved per node as
   //   f(k+1) = f(k) * τ / (k+1)
   std::vector<real_t> fac(O_nodes, 1.0);

   for (int k = 0; k < order - 1; k++)
   {
      // D_next = L(D_curr) = -Σ_d A_d · ∂_{x_d} D_curr.  Same kernel as
      // ComputeADERTimeIntegrated; element-local, no cross-element coupling.
      D_next = 0.0;
      for (int d = 0; d < 3; d++)
      {
         ApplySpatialDerivative(d, D_curr, dQ_dxd);
         const DenseMatrix &A_d = flux_.GetReferenceStarMatrix(d);
         ApplyJacobianPerDOF(A_d, dQ_dxd, D_next, ndof_total_, /*sign=*/-1.0);
      }

      // Update factorial factors and accumulate D(k+1) into each node.
      //
      // R-1503 cross-reference: `ComputeADERTimeIntegrated` above uses
      // `fac *= dt / (k+2)`.  The two recursions parameterise the factorial
      // at different phases:
      //   - Here `fac[o]` was initialised to 1.0 at L1093 and the k=0
      //     contribution `(τ^0/0!) · D(0) = D(0)` was pre-added at
      //     L1086-1089 WITHOUT scaling.  At iteration k we are forming
      //     D(k+1) and accumulating `(τ^{k+1}/(k+1)!) · D(k+1)`, so the
      //     update is `fac[o] *= τ / (k+1)` — i.e. `denom = k + 1`.
      //   - `ComputeADERTimeIntegrated` initialised `fac = dt` at L1004
      //     (its k=0 contribution was pre-added with `dt` already folded in).
      //     There the iteration computes D(k+1) and contributes
      //     `dt^{k+2}/(k+2)! · D(k+1)`, hence `fac *= dt / (k+2)`.
      // Both are mathematically equivalent.  DO NOT "harmonise" the two
      // denominators without re-deriving the recursion phase — see R-1503.
      const real_t denom = static_cast<real_t>(k + 1);
      for (int o = 0; o < O_nodes; o++)
      {
         fac[o] *= tau_nodes[o] / denom;
         Q_per_node[o].Add(fac[o], D_next);
      }

      // Swap buffers: D_curr <- D_next for the next iteration.
      mfem::Swap(D_curr, D_next);
   }
}

// ---------------------------------------------------------------------------
// R-602/R-603: SubStep iterator side-channel setters.  See header docstring.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::SetSubStepFaultImposedStates(
   const real_t *I_imp_plus_flat, const real_t *I_imp_minus_flat,
   int n_total_fault_qps) const
{
   MFEM_VERIFY(n_total_fault_qps >= 0,
               "SetSubStepFaultImposedStates: n_total_fault_qps must be "
               ">= 0, got " << n_total_fault_qps);
   if (n_total_fault_qps > 0)
   {
      MFEM_VERIFY(I_imp_plus_flat != nullptr && I_imp_minus_flat != nullptr,
                  "SetSubStepFaultImposedStates: both pointers must be "
                  "non-null when n_total_fault_qps > 0");
   }
   substep_I_imp_plus_flat_  = I_imp_plus_flat;
   substep_I_imp_minus_flat_ = I_imp_minus_flat;
   substep_n_total_fault_qps_ = n_total_fault_qps;
}

template <typename MeshType>
void WaveOperator<MeshType>::ResetSubStepFaultImposedStates() const
{
   substep_I_imp_plus_flat_  = nullptr;
   substep_I_imp_minus_flat_ = nullptr;
   substep_n_total_fault_qps_ = 0;
}

// ---------------------------------------------------------------------------
// Round-11 Mixed-Flux mode setter (Zhang et al. 2023).  Default `None` ==
// pre-Mixed-Flux behavior.  Mode flip rebuilds central_flux_face_set_.
// R-1203 cross-check with precomputed-flux state; symmetric guard in
// UsePrecomputedFaceFluxes.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::SetMixedFluxMode(MixedFluxMode m)
{
   // R-1205: enforce collective-call discipline.  All ranks must invoke
   // SetMixedFluxMode with the SAME mode in the same order; otherwise
   // BuildCentralFluxFaceSet_ deadlocks at MPI_Allgather (Adjacent runs
   // collectives, None / AllContinuous skip them).  Allreduce-MIN/MAX
   // over the mode value catches mismatched callers cheaply at the
   // setter rather than at a downstream collective.
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto &pmesh_consensus = static_cast<ParMesh &>(mesh_);
      const int my_mode = static_cast<int>(m);
      int max_mode = 0, min_mode = 0;
      MPI_Allreduce(&my_mode, &max_mode, 1, MPI_INT, MPI_MAX,
                    pmesh_consensus.GetComm());
      MPI_Allreduce(&my_mode, &min_mode, 1, MPI_INT, MPI_MIN,
                    pmesh_consensus.GetComm());
      MFEM_VERIFY(max_mode == min_mode,
                  "SetMixedFluxMode: ranks disagree on the mode value "
                  "(max=" << max_mode << ", min=" << min_mode
                  << ").  Every rank must call SetMixedFluxMode with "
                  "the SAME MixedFluxMode in the same call order; "
                  "otherwise the post-walk Allgatherv deadlocks.");
#endif
   }

   // R-1203 cross-check.
   if (m != MixedFluxMode::None && use_precomputed_face_fluxes_)
   {
      MFEM_ABORT("SetMixedFluxMode("
                 << (m == MixedFluxMode::Adjacent ? "Adjacent"
                                                  : "AllContinuous")
                 << "): mutually exclusive with precomputed face flux "
                 "(UsePrecomputedFaceFluxes is currently true).  "
                 "Disable precomputed flux first.");
   }

   // R-1205: Adjacent mode requires fault attribute and a populated fault
   // face list.  AllContinuous works even without a fault.
   if (m == MixedFluxMode::Adjacent)
   {
      MFEM_VERIFY(bc_.fault_attr > 0,
                  "SetMixedFluxMode(Adjacent): bc_.fault_attr = "
                  << bc_.fault_attr << " (≤ 0).  Adjacent mode requires "
                  "a fault attribute; use AllContinuous or None for "
                  "no-fault meshes.");
   }

   mixed_flux_mode_ = m;
   mf_on_ = (m != MixedFluxMode::None);   // R-1208: keep cached flag in sync.
   BuildCentralFluxFaceSet_();
}

// ---------------------------------------------------------------------------
// BuildCentralFluxFaceSet_ — populate central_flux_face_set_ per the
// configured mixed-flux mode.  Algorithm follows Zhang 2023 §3.2 / Fig 5.
// R-1207: explicit "is interior" predicate that handles MPI shared seams.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::BuildCentralFluxFaceSet_()
{
   central_flux_face_set_.clear();
   if (mixed_flux_mode_ == MixedFluxMode::None) { return; }

   // R-1207 membership predicate factored as a lambda for reuse.
   auto is_interior_face = [&](int f) -> bool
   {
      FaceElementTransformations *ftr =
         mesh_.GetFaceElementTransformations(f);
      const bool two_sided_interior =
         (ftr != nullptr && ftr->Elem2No >= 0);
      const bool shared_seam =
         (shared_mesh_face_set_.count(f) > 0);
      // R-1404: exclude faces carrying a NON-fault boundary attribute
      // (free-surface, absorbing).  Such faces fall through to the BC
      // branch only when Elem2No < 0; if a 2-sided face happens to also
      // carry a non-zero non-fault bdr_attr (custom mesh with internal
      // boundary layer), the BC branch's `is_boundary = (e2 < 0) &&
      // (bdr_attr > 0)` skips them and the mixed-flux dispatch would
      // wrongly fire Central instead of the user-intended BC.  Defensive:
      // require absence of non-fault boundary attribute for membership
      // in the central set.
      const bool nonfault_bc =
         (f < static_cast<int>(face_bdr_attr_.size()) &&
          face_bdr_attr_[f] > 0 &&
          face_bdr_attr_[f] != bc_.fault_attr);
      return (two_sided_interior || shared_seam) && !nonfault_bc;
   };
   // R-1100 fix: build a mesh-face-index keyed fault-face set that is
   // rank-symmetric.  `face_bdr_attr_` is populated only on the BE-owner
   // rank, so on a shared fault face the non-BE-owner's `face_bdr_attr_`
   // entry is 0 — `is_fault_face` would mis-classify the shared fault
   // face as non-fault and the local Adjacent walk would wrongly insert
   // it into central_flux_face_set_, breaking the documented invariant
   // `central_flux_face_set_ ∩ fault_faces == ∅` on the non-BE-owner.
   //
   // Reuse the canonical rank-symmetric lists `fault_interior_faces_`
   // and `fault_shared_faces_` (both populated correctly at ctor R-002
   // stage via the merged `shared_face_bdr_attr_` exchange).
   // R-1206: build the set unconditionally, then enforce ctor-state
   // consistency: a non-empty fault-face list requires a positive
   // fault attribute.  Lambdas no longer need the redundant
   // `bc_.fault_attr > 0` short-circuit.
   std::unordered_set<int> fault_mesh_face_idx_set;
   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      fault_mesh_face_idx_set.insert(fault_interior_faces_[i]);
   }
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto &pmesh_lookup = static_cast<ParMesh &>(mesh_);
      for (int sf_i = 0; sf_i < fault_shared_faces_.Size(); sf_i++)
      {
         const int sf = fault_shared_faces_[sf_i];
         const int f = pmesh_lookup.GetSharedFace(sf);
         fault_mesh_face_idx_set.insert(f);
      }
#endif
   }
   MFEM_VERIFY(fault_mesh_face_idx_set.empty() || bc_.fault_attr > 0,
               "BuildCentralFluxFaceSet_: fault-face lists are non-empty "
               "(set size = " << fault_mesh_face_idx_set.size()
               << ") but bc_.fault_attr = " << bc_.fault_attr
               << " (≤ 0).  Inconsistent ctor state — the fault face "
               "lists must be empty when no fault attribute is set "
               "(R-1206 invariant).");
   auto is_fault_face = [&](int f) -> bool
   {
      // R-1206 + R-1401: rely on set membership AS THE TRUTH SOURCE,
      // but assert agreement with the local `face_bdr_attr_` table for
      // BE-owner ranks (debug-only).  If a future ctor refactor
      // decouples `fault_mesh_face_idx_set` from `face_bdr_attr_`, the
      // dispatch-time fault check (which uses face_bdr_attr_) and this
      // set-build check would diverge — central could be dispatched on
      // a face that the fault path also activates, double-counting.
      const bool by_set = (fault_mesh_face_idx_set.count(f) > 0);
#ifndef NDEBUG
      // The bdr-attr table is rank-asymmetric (BE-owner only) on
      // shared fault faces; only assert the implication "by_attr →
      // by_set" (the strict direction), not equivalence.
      const bool by_attr =
         (bc_.fault_attr > 0 &&
          f < static_cast<int>(face_bdr_attr_.size()) &&
          face_bdr_attr_[f] == bc_.fault_attr);
      MFEM_ASSERT(!by_attr || by_set,
                  "R-1401: face " << f << " has fault bdr_attr but is "
                  "NOT in fault_mesh_face_idx_set — ctor refactor has "
                  "decoupled the two; mixed-flux dispatch would "
                  "double-dispatch on this face.");
#endif
      return by_set;
   };

   if (mixed_flux_mode_ == MixedFluxMode::AllContinuous)
   {
      // Mixed-Flux 1 (Zhang Fig 4a): central on every interior non-fault.
      for (int f = 0; f < mesh_.GetNumFaces(); f++)
      {
         if (is_interior_face(f) && !is_fault_face(f))
         {
            central_flux_face_set_.insert(f);
         }
      }
      return;
   }

   // Mixed-Flux 2 (Zhang Fig 4b): central only on faces immediately
   // adjacent to a fault element.

   // Step 1: build E_fault_adj from interior + shared fault faces.
   std::unordered_set<int> E_fault_adj;
   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      const int f = fault_interior_faces_[i];
      FaceElementTransformations *ftr =
         mesh_.GetInteriorFaceTransformations(f);
      if (!ftr) { continue; }
      if (ftr->Elem1No >= 0) { E_fault_adj.insert(ftr->Elem1No); }
      if (ftr->Elem2No >= 0) { E_fault_adj.insert(ftr->Elem2No); }
   }
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto &pmesh = static_cast<ParMesh &>(mesh_);
      for (int sf_i = 0; sf_i < fault_shared_faces_.Size(); sf_i++)
      {
         const int sf = fault_shared_faces_[sf_i];
         FaceElementTransformations *ftr =
            pmesh.GetSharedFaceTransformations(sf);
         if (!ftr) { continue; }
         if (ftr->Elem1No >= 0) { E_fault_adj.insert(ftr->Elem1No); }
         // Elem2 of a shared face is on a remote rank (Elem2No < 0 with
         // ghost-block index); only this rank's Elem1 is added.
      }
#endif
   }

   // Step 2: for each element in E_fault_adj, walk its faces; insert any
   // interior non-fault face into central_flux_face_set_.
   Array<int> elem_faces, elem_orient;
   for (int e : E_fault_adj)
   {
      mesh_.GetElementFaces(e, elem_faces, elem_orient);
      for (int j = 0; j < elem_faces.Size(); j++)
      {
         const int f = elem_faces[j];
         if (is_interior_face(f) && !is_fault_face(f))
         {
            central_flux_face_set_.insert(f);
         }
      }
   }

   // R-001 fix (post-local-walk MPI exchange): the local walk above only
   // inserts shared non-fault faces whose LOCAL element is in E_fault_adj.
   // A shared non-fault face whose REMOTE element is fault-adjacent on
   // the peer rank (but local element isn't) would be missed, causing
   // the two ranks to dispatch DIFFERENT fluxes on the same physical
   // face → conservation broken across the rank seam (a 0.5·|A_n|·jump
   // imbalance per ill-dispatched face per macro-step).
   //
   // Allgatherv the global vertex keys of every shared non-fault face
   // this rank has inserted; every rank then re-walks its shared faces
   // and inserts any whose key appears in the global merged set.  Same
   // pattern as the ctor's `global_fault_keys` exchange (lines 159-219).
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto &pmesh = static_cast<ParMesh &>(mesh_);
      const int n_shared = pmesh.GetNSharedFaces();
      Array<HYPRE_BigInt> gvi;
      pmesh.GetGlobalVertexIndices(gvi);
      auto make_global_key = [&](const Array<int> &verts)
      {
         // R-1106: guard against silent truncation if a future mesh has
         // faces with > 4 vertices (impossible for tet/hex but possible
         // for polygonal elements).
         MFEM_VERIFY(verts.Size() <= 4,
                     "BuildCentralFluxFaceSet_::make_global_key: face has "
                     << verts.Size() << " vertices > 4; the "
                     "std::array<HYPRE_BigInt,4> key bucket would "
                     "silently truncate.  Extend the bucket size or use "
                     "std::vector<HYPRE_BigInt> if higher-vertex faces "
                     "are introduced.");
         // R-1207 (deferred): for a hybrid mesh with both triangular
         // and quadrilateral faces, the key bucket size (4) is wide
         // enough to hold either, but a 3-vert triangle key
         // {0, v1, v2, v3} (sorted) could in principle alias against a
         // quad key whose sorted form happens to start with 0.  TPV104
         // production uses tet meshes only (all faces triangular), so
         // the alias is impossible by construction.  When a hybrid
         // mesh fixture is added, prepend `verts.Size()` as a
         // disambiguator slot per the R-1207 suggested fix.
         std::array<HYPRE_BigInt, 4> key = {0, 0, 0, 0};
         for (int v = 0; v < std::min(verts.Size(), 4); v++)
         {
            key[v] = gvi[verts[v]];
         }
         std::sort(key.begin(), key.begin() + verts.Size());
         return key;
      };

      // Collect this rank's "shared face inserted into central set" keys.
      std::set<std::array<HYPRE_BigInt, 4>> local_keys;
      for (int sf = 0; sf < n_shared; sf++)
      {
         const int f = pmesh.GetSharedFace(sf);
         if (central_flux_face_set_.count(f) > 0)
         {
            Array<int> verts;
            mesh_.GetFaceVertices(f, verts);
            local_keys.insert(make_global_key(verts));
         }
      }

      // Allgatherv local_keys → global_keys (collective; must run on all
      // ranks, including those with my_size = 0).
      std::vector<HYPRE_BigInt> local_flat;
      local_flat.reserve(local_keys.size() * 4);
      for (const auto &k : local_keys)
      {
         for (int i = 0; i < 4; i++) { local_flat.push_back(k[i]); }
      }
      const int my_size = static_cast<int>(local_flat.size());
      int nprocs_loc = 0;
      MPI_Comm_size(pmesh.GetComm(), &nprocs_loc);
      std::vector<int> sizes(nprocs_loc), displs(nprocs_loc);
      MPI_Allgather(&my_size, 1, MPI_INT, sizes.data(), 1, MPI_INT,
                    pmesh.GetComm());
      int total = 0;
      for (int r = 0; r < nprocs_loc; r++)
      { displs[r] = total; total += sizes[r]; }
      std::vector<HYPRE_BigInt> all_flat(total);
      MPI_Allgatherv(local_flat.data(), my_size,
                     MPITypeMap<HYPRE_BigInt>::mpi_type,
                     all_flat.data(), sizes.data(), displs.data(),
                     MPITypeMap<HYPRE_BigInt>::mpi_type,
                     pmesh.GetComm());
      std::set<std::array<HYPRE_BigInt, 4>> global_keys;
      for (int i = 0; i < total; i += 4)
      {
         std::array<HYPRE_BigInt, 4> key = {
            all_flat[i], all_flat[i+1], all_flat[i+2], all_flat[i+3]
         };
         global_keys.insert(key);
      }

      // Insert any shared non-fault face whose key is in the merged
      // global set (idempotent for keys this rank already contributed).
      for (int sf = 0; sf < n_shared; sf++)
      {
         const int f = pmesh.GetSharedFace(sf);
         if (is_fault_face(f)) { continue; }
         Array<int> verts;
         mesh_.GetFaceVertices(f, verts);
         if (global_keys.count(make_global_key(verts)) > 0)
         {
            central_flux_face_set_.insert(f);
         }
      }
#endif
   }
}

// ---------------------------------------------------------------------------
// EvaluateBulkAtFaultQPsCanonical — read bulk Q at every fault QP (interior
// + shared, R-1003), rotate to the canonical fault-local frame, route Elem1's
// evaluation into the canonical-+/− output bucket, and pack into flat per-QP
// arrays.
//
// Layout:
//   Q_*_flat[ dof_idx * NUM_STATE + c ]
// where
//   dof_idx = fault_face_dof_offset_[f]   + q   for interior fault face f
//   dof_idx = shared_fault_dof_offset_[sf] + q   for shared   fault face sf
// (the latter map already includes the GetNumLocalFaultQPs() base offset by
// construction; see SetFaultDOFData in wave_operator.hpp).  This is the same
// indexing used by the fault branches of ComputeADERFaceFluxRHS and
// ComputeADERSharedFaceFluxRHS, so the iterator's I_imp accumulator and the
// flux's substep gate share a single absolute index.
//
// Reuses the canonical frame from FaultBasis (sign_flipped reconstruction)
// and the per-face elem1_on_plus flag (R-101 for interior; geometric for
// shared) so the rotation and side-labeling are bit-identical to the
// production fault flux path on each branch.  The interior loop uses
// `!interior_fault_elem1_on_plus_[fi]` for frame negation (R-101); the
// shared loop uses `qpd.sign_flipped` to mirror ComputeADERSharedFaceFluxRHS
// (R-1305 / R-801 documents this convention split as a separate latent
// risk; matching the existing flux convention is required so the iterator's
// canonical Q matches what the flux's substep gate will rotate back).
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::EvaluateBulkAtFaultQPsCanonical(
   const Vector &Q_bulk,
   std::vector<real_t> &Q_plus_flat,
   std::vector<real_t> &Q_minus_flat) const
{
   MFEM_VERIFY(Q_bulk.Size() == NUM_STATE * ndof_total_,
               "EvaluateBulkAtFaultQPsCanonical: Q_bulk size "
               << Q_bulk.Size() << " != NUM_STATE * ndof_total_ = "
               << NUM_STATE * ndof_total_);

   // R-1003 §4: size to TOTAL fault QPs (interior + shared).  Interior
   // entries occupy [0, GetNumLocalFaultQPs()); shared entries occupy
   // [GetNumLocalFaultQPs(), GetNumTotalFaultQPs()).  At np=1 the shared
   // count is 0 so this is byte-identical to the pre-R-1003 sizing.
   const int n_local_qps = GetNumLocalFaultQPs();
   const int n_total_qps = GetNumTotalFaultQPs();
   const size_t expect_words =
      static_cast<size_t>(NUM_STATE) * static_cast<size_t>(n_total_qps);
   Q_plus_flat.assign(expect_words, 0.0);
   Q_minus_flat.assign(expect_words, 0.0);

   // R-1600: at np>1 the parallel block below performs a per-substep
   // PAIRWISE collective (`q_gf.ExchangeFaceNbrData`) that EVERY rank
   // with any face neighbours must participate in, even ranks with
   // zero fault QPs (and zero `fault_basis_` — the ctor only allocates
   // fault_basis_ on ranks with at least one local fault face, see
   // wave_operator.inl ctor L332-333).  Pre-R-1600 this site short-
   // circuited on `n_total_qps == 0` (or equivalently `!fault_basis_`),
   // causing fault-adjacent ranks to hang in MPI_Wait at np>1 (the
   // production 13.5-min spin).
   //
   // On serial builds the parallel block is `if constexpr` skipped at
   // compile time, so a no-fault-QPs rank can return early without harm.
   // The interior and shared loops below are bounded by
   // `fault_interior_faces_.Size()` and `fault_shared_faces_.Size()`
   // respectively — both are 0 when `fault_basis_` is null, so the loop
   // bodies (which DO dereference `fault_basis_`) never run on those
   // ranks.  No early return on `!fault_basis_` is needed.
   if (n_total_qps == 0 && !IsParallelMesh<MeshType>::value) { return; }

   const real_t *Q_data = Q_bulk.GetData();

   for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
   {
      const int f = fault_interior_faces_[fi];
      FaceElementTransformations *ftr =
         mesh_.GetInteriorFaceTransformations(f);
      // R-1004: ftr should NEVER be null on a face in fault_interior_faces_
      // (the ctor populates this list filtering on
      // GetInteriorFaceTransformations(f) != nullptr).  Silently skipping
      // here would leave the corresponding entries of Q_*_flat at zero,
      // producing wrong friction inputs for the substep iterator.  Abort
      // loudly instead.
      MFEM_VERIFY(ftr != nullptr,
                  "EvaluateBulkAtFaultQPsCanonical: face " << f
                  << " (fault_interior_faces_[" << fi << "]) has no "
                  "InteriorFaceTransformations.  This violates the "
                  "invariant established at WaveOperator ctor; the "
                  "fault face list has been corrupted.");

      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      const FiniteElement *fe1 = fes_->GetFE(e1);
      const FiniteElement *fe2 = fes_->GetFE(e2);
      const int ndof = fe1->GetDof();
      MFEM_ASSERT(ndof == ndof_per_el_,
                  "Heterogeneous DOF counts not supported in "
                  "EvaluateBulkAtFaultQPsCanonical.");
      const int dof_offset1 = e1 * ndof_per_el_;
      const int dof_offset2 = e2 * ndof_per_el_;

      const IntegrationRule &ir = IntRules.Get(
         ftr->GetGeometryType(), 2 * order_);
      const int nqp = ir.GetNPoints();
      MFEM_VERIFY(nqp == nbf_per_face_,
                  "EvaluateBulkAtFaultQPsCanonical: face nqp "
                  << nqp << " != nbf_per_face_ " << nbf_per_face_);

      auto it_off = fault_face_dof_offset_.find(f);
      MFEM_VERIFY(it_off != fault_face_dof_offset_.end(),
                  "EvaluateBulkAtFaultQPsCanonical: fault_face_dof_offset_"
                  " missing entry for fault face " << f);
      const int base_dof_idx = it_off->second;

      // R-101: per-face geometric side label (matches ComputeADERFaceFluxRHS).
      MFEM_VERIFY(fi >= 0 &&
                  fi < static_cast<int>(interior_fault_elem1_on_plus_.size()),
                  "EvaluateBulkAtFaultQPsCanonical: "
                  "interior_fault_elem1_on_plus_ size mismatch");
      const bool elem1_on_plus = interior_fault_elem1_on_plus_[fi];

      // FaultBasis index for this interior fault face.
      const int fb_idx = LookupInteriorFaultBasisIndex(f);
      MFEM_VERIFY(fb_idx >= 0 && fb_idx < fault_basis_->NumFaces(),
                  "EvaluateBulkAtFaultQPsCanonical: invalid FaultBasis "
                  "index for face " << f);
      const FaultBasisData &bd = fault_basis_->GetBasis(fb_idx);

      for (int q = 0; q < nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         IntegrationPoint ip1, ip2;
         ftr->Loc1.Transform(ip, ip1);
         ftr->Loc2.Transform(ip, ip2);
         Vector shape1(ndof), shape2(ndof);
         fe1->CalcShape(ip1, shape1);
         fe2->CalcShape(ip2, shape2);

         // Bulk Q at the QP from each side, in the global frame.
         real_t Q_self[NUM_STATE], Q_nbr[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t s_self = 0.0, s_nbr = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               s_self += shape1(i) * Q_data[c * ndof_total_ + dof_offset1 + i];
               s_nbr  += shape2(i) * Q_data[c * ndof_total_ + dof_offset2 + i];
            }
            Q_self[c] = s_self;
            Q_nbr[c]  = s_nbr;
         }

         // Reconstruct the canonical frame from FaultBasis (sign-corrected).
         MFEM_VERIFY(q < static_cast<int>(bd.qp_data.size()),
                     "EvaluateBulkAtFaultQPsCanonical: qp_data missing q="
                     << q << " for face " << f);
         const FaultBasisQPData &qpd = bd.qp_data[q];
         real_t can_n[3], can_t1[3], can_t2[3];
         const bool should_negate_frame = !elem1_on_plus;
         for (int d = 0; d < 3; d++)
         {
            can_n[d]  = should_negate_frame ? -qpd.normal[d]
                                            :  qpd.normal[d];
            can_t1[d] = should_negate_frame ? -qpd.tangent1[d]
                                            :  qpd.tangent1[d];
            can_t2[d] = should_negate_frame ? -qpd.tangent2[d]
                                            :  qpd.tangent2[d];
         }
         DenseMatrix Tinv_can(NUM_STATE);
         GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv_can);

         real_t Q_self_can[NUM_STATE], Q_nbr_can[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t s_self = 0.0, s_nbr = 0.0;
            for (int k = 0; k < NUM_STATE; k++)
            {
               s_self += Tinv_can(c, k) * Q_self[k];
               s_nbr  += Tinv_can(c, k) * Q_nbr[k];
            }
            Q_self_can[c] = s_self;
            Q_nbr_can[c]  = s_nbr;
         }

         // Route Elem1's evaluation into the +/− bucket per elem1_on_plus.
         const real_t *src_plus  = elem1_on_plus ? Q_self_can : Q_nbr_can;
         const real_t *src_minus = elem1_on_plus ? Q_nbr_can  : Q_self_can;

         const int dof_idx = base_dof_idx + q;
         real_t *dst_p = Q_plus_flat.data()  + dof_idx * NUM_STATE;
         real_t *dst_m = Q_minus_flat.data() + dof_idx * NUM_STATE;
         for (int c = 0; c < NUM_STATE; c++)
         {
            dst_p[c] = src_plus[c];
            dst_m[c] = src_minus[c];
         }
      }
   }

   // R-1003 §3: shared-fault loop.  Populates the
   // [GetNumLocalFaultQPs(), GetNumTotalFaultQPs()) slice of Q_*_flat from
   // self-side (locally owned) Q + neighbour-side (face-nbr ghost layer) Q,
   // matching the canonical-frame rotation in ComputeADERSharedFaceFluxRHS
   // (lines ~3998-4034).  At np=1 fault_shared_faces_ is empty and this
   // block is a no-op (no MPI, no ghost exchange).  At np>1 we perform one
   // per-component ExchangeFaceNbrData per call (NUM_STATE collectives) so
   // the substep predictor's per-sub-step Q lands on neighbour ranks; this
   // is R-1003 §1 (the per-substep ghost coverage required by R-1303 — the
   // CK recursion in ComputeADERSubStepStates is element-local and does
   // not populate ghost cells on its own).  Cost: O macro-step calls ×
   // NUM_STATE collectives = 9·O exchanges per macro-step.
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      // R-1600: gate the per-substep COLLECTIVE on TOTAL shared faces
      // (`pmesh.GetNSharedFaces()`), NOT on the fault-only shared count
      // (`fault_shared_faces_.Size()`).  `q_gf.ExchangeFaceNbrData()` is
      // a pairwise MPI exchange among the rank's face-neighbour set: a
      // rank that has shared NON-FAULT faces with a fault-adjacent peer
      // MUST still call the exchange, otherwise its peers' MPI_Irecv
      // never gets paired and they hang in MPI_Wait at 100% CPU.  Pre-
      // R-1600 this site gated on `fault_shared_faces_.Size() > 0`,
      // producing the production-blocking 13.5-min hang at np=10
      // observed on the symmirror 1000m mesh (every interior rank with
      // shared non-fault faces but no fault faces of its own would
      // skip, deadlocking the fault-adjacent ranks).  The mirror site
      // ComputeADERSharedFaceFluxRHS uses the same `GetNSharedFaces()`
      // gate (`wave_operator.inl:4167`) and was never affected.
      auto &pmesh = const_cast<ParMesh &>(
         static_cast<const ParMesh &>(mesh_));
      const int n_shared_total = pmesh.GetNSharedFaces();
      const int n_shared_faces = fault_shared_faces_.Size();
      if (n_shared_total > 0)
      {
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         MFEM_VERIFY(pfes,
                     "EvaluateBulkAtFaultQPsCanonical: FESpace must be "
                     "ParFiniteElementSpace for ParMesh.");

         // R-1601: batched per-substep ghost exchange via the
         // vdim=NUM_STATE byNODES `ghost_gf_full_state_` ParGridFunction
         // built once in the ctor (see wave_operator.hpp:755-770 for
         // rationale).  Pre-R-1601 this site looped 9 sequential
         // single-component `q_gf.ExchangeFaceNbrData()` calls (~9·O
         // collectives per macro-step).  Now the byNODES vdim=NUM_STATE
         // PGF carries all 9 components contiguously and one
         // ExchangeFaceNbrData call covers them all → O collectives per
         // macro-step → saves 8·O.  R-1600 contract is preserved: every
         // rank with `pmesh.GetNSharedFaces() > 0` calls the exchange
         // exactly once per invocation; the per-fault-face loop below
         // remains gated on `fault_shared_faces_.Size() > 0`.
         MFEM_VERIFY(ghost_gf_full_state_,
                     "EvaluateBulkAtFaultQPsCanonical: ghost_gf_full_state_ "
                     "not initialised (ctor's ParMesh path should have "
                     "built it; R-1601).");
         ParGridFunction &q_gf_full = *ghost_gf_full_state_;

         // Pack Q (component-major) into the byNODES vdim=NUM_STATE PGF.
         // Both layouts are component-major so the copy is contiguous
         // per component; total size NUM_STATE * ndof_total_ matches Q.
         {
            real_t *q_full_data = q_gf_full.GetData();
            std::memcpy(q_full_data, Q_data,
                        static_cast<size_t>(NUM_STATE) *
                        static_cast<size_t>(ndof_total_) * sizeof(real_t));
         }

         q_gf_full.ExchangeFaceNbrData();   // 1 collective for all NUM_STATE components.

         // Unpack FaceNbrData into the per-component nbr_data[c] arrays
         // expected by the per-fault-face loop below.  byNODES layout in
         // FaceNbrData: `src[c * n_face_nbr_dofs + i]`.  We deep-copy
         // (rather than alias) because subsequent `q_gf_full` reads in
         // the loop body could in principle invalidate the storage —
         // matches the pre-R-1601 deep-copy semantics.
         const Vector &src = q_gf_full.FaceNbrData();
         MFEM_VERIFY(src.Size() % NUM_STATE == 0,
                     "EvaluateBulkAtFaultQPsCanonical: FaceNbrData size "
                     << src.Size() << " is not a multiple of NUM_STATE = "
                     << NUM_STATE);
         const int n_face_nbr_dofs = src.Size() / NUM_STATE;

         std::vector<Vector> nbr_data(NUM_STATE);
         for (int c = 0; c < NUM_STATE; c++)
         {
            nbr_data[c].SetSize(n_face_nbr_dofs);
            std::memcpy(nbr_data[c].GetData(),
                        src.GetData() + c * n_face_nbr_dofs,
                        static_cast<size_t>(n_face_nbr_dofs) * sizeof(real_t));
         }

         // R-1600: a rank with shared faces but NO fault-shared faces
         // has now satisfied the pairwise-collective contract for its
         // peers.  No per-fault-face work to do here — fall through to
         // the function tail.  The output buffer is already zeroed at
         // the function entry (Q_*_flat.assign(expect_words, 0.0) above).
         for (int sf_idx = 0; sf_idx < n_shared_faces; sf_idx++)
         {
            const int sf = fault_shared_faces_[sf_idx];
            // R-1603 (defer-with-rationale): `pmesh.GetSharedFaceTransformations`
            // returns a pointer to a `FaceElementTransformations` object held
            // in the ParMesh's INTERNAL storage; that object is reused across
            // calls (the next call with a different `sf` overwrites its
            // members).  A naïve "cache the pointer per sf" therefore
            // aliases stale state.  A correct cache would deep-copy the FET
            // (Loc1/Loc2 transforms, geometry type, Elem1No/Elem2No) — non-
            // trivial because FET members include MFEM-internal pointers.
            // At TPV104 production scale (O × n_shared_faces per macro-step
            // ≈ 4 × ~24 = 96 calls per macro-step), the cost is dominated by
            // R-1601's collectives, not by these constructions.  Defer the
            // FET cache until profiling shows it dominates after the R-1601
            // batched-exchange landing.
            FaceElementTransformations *ftr =
               pmesh.GetSharedFaceTransformations(sf);
            // R-1004-style invariant: the ctor only adds entries to
            // fault_shared_faces_ when GetSharedFaceTransformations
            // succeeds, so a null here means the list has been corrupted.
            MFEM_VERIFY(ftr != nullptr,
                        "EvaluateBulkAtFaultQPsCanonical: shared face sf="
                        << sf << " (fault_shared_faces_[" << sf_idx
                        << "]) has no SharedFaceTransformations.");

            const int e1 = ftr->Elem1No;
            const FiniteElement *fe1 = fes_->GetFE(e1);
            const int ndof = fe1->GetDof();
            MFEM_ASSERT(ndof == ndof_per_el_,
                        "Heterogeneous DOF counts not supported in "
                        "EvaluateBulkAtFaultQPsCanonical (shared self).");
            const int dof_offset1 = e1 * ndof_per_el_;

            const int nbr_idx = ftr->Elem2No - ne_;
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
            const int ndof2 = fe2->GetDof();
            MFEM_ASSERT(ndof2 == ndof_per_el_,
                        "Heterogeneous DOF counts not supported in "
                        "EvaluateBulkAtFaultQPsCanonical (shared ghost).");

            const IntegrationRule &ir = IntRules.Get(
               ftr->GetGeometryType(), 2 * order_);
            const int nqp = ir.GetNPoints();
            MFEM_VERIFY(nqp == nbf_per_face_,
                        "EvaluateBulkAtFaultQPsCanonical (shared): face "
                        "nqp " << nqp << " != nbf_per_face_ "
                        << nbf_per_face_);

            auto it_off = shared_fault_dof_offset_.find(sf);
            MFEM_VERIFY(it_off != shared_fault_dof_offset_.end(),
                        "EvaluateBulkAtFaultQPsCanonical: "
                        "shared_fault_dof_offset_ missing entry for "
                        "shared fault face sf=" << sf);
            const int base_dof_idx = it_off->second;
            // Sanity: shared offsets must lie in the upper slice.
            MFEM_ASSERT(base_dof_idx >= n_local_qps &&
                        base_dof_idx + nqp <= n_total_qps,
                        "EvaluateBulkAtFaultQPsCanonical: shared offset "
                        "out of [n_local_qps, n_total_qps).");

            MFEM_VERIFY(sf_idx >= 0 &&
                        sf_idx < static_cast<int>(
                                    shared_fault_elem1_on_plus_.size()),
                        "EvaluateBulkAtFaultQPsCanonical: "
                        "shared_fault_elem1_on_plus_ size mismatch.");
            const bool elem1_on_plus =
               shared_fault_elem1_on_plus_[sf_idx];

            // FaultBasis indexes interior faces first, then shared.
            const int fb_idx = fault_interior_faces_.Size() + sf_idx;
            MFEM_VERIFY(fb_idx >= 0 && fb_idx < fault_basis_->NumFaces(),
                        "EvaluateBulkAtFaultQPsCanonical: invalid "
                        "FaultBasis index " << fb_idx
                        << " for shared face sf=" << sf);
            const FaultBasisData &bd = fault_basis_->GetBasis(fb_idx);

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               ftr->SetAllIntPoints(&ip);

               IntegrationPoint ip1, ip2;
               ftr->Loc1.Transform(ip, ip1);
               ftr->Loc2.Transform(ip, ip2);
               Vector shape1(ndof), shape2(ndof2);
               fe1->CalcShape(ip1, shape1);
               fe2->CalcShape(ip2, shape2);

               // Self side: locally-owned Q via dof_offset1.  Neighbour
               // side: face-nbr ghost layer via nbr_idx * ndof_per_el_.
               // Same indexing as ComputeADERSharedFaceFluxRHS:3954-3976.
               real_t Q_self[NUM_STATE], Q_nbr[NUM_STATE];
               for (int c = 0; c < NUM_STATE; c++)
               {
                  real_t s_self = 0.0, s_nbr = 0.0;
                  for (int i = 0; i < ndof; i++)
                  {
                     s_self += shape1(i)
                               * Q_data[c * ndof_total_ + dof_offset1 + i];
                  }
                  for (int i = 0; i < ndof2; i++)
                  {
                     s_nbr += shape2(i)
                              * nbr_data[c][nbr_idx * ndof_per_el_ + i];
                  }
                  Q_self[c] = s_self;
                  Q_nbr[c]  = s_nbr;
               }

               MFEM_VERIFY(q < static_cast<int>(bd.qp_data.size()),
                           "EvaluateBulkAtFaultQPsCanonical (shared): "
                           "qp_data missing q=" << q << " for sf=" << sf);
               const FaultBasisQPData &qpd = bd.qp_data[q];

               // R-1305 / R-801 caveat: the existing shared-fault ADER
               // branch (ComputeADERSharedFaceFluxRHS:4007-4012) gates
               // frame negation on `qpd.sign_flipped`, while the interior
               // ADER branch uses the R-101 geometric flag
               // `!interior_fault_elem1_on_plus_`.  Mirror the shared-
               // flux convention here so the iterator's canonical Q
               // matches the rotation T_can the substep gate uses to
               // map I_imp_± back to global on the corresponding shared
               // QPs (R-1003 invariant: same Tinv_can on both ends).
               // sign_flipped is RANK-INDEPENDENT (set at FaultBasis
               // ctor from face geometry); !elem1_on_plus is RANK-
               // DEPENDENT (rank A's Elem1 = rank B's Elem2 → opposite
               // canonical frames on the same physical face).  Using
               // sign_flipped keeps the canonical frame consistent
               // across the partition seam.
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

               DenseMatrix Tinv_can(NUM_STATE);
               GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2,
                                                 Tinv_can);

               real_t Q_self_can[NUM_STATE], Q_nbr_can[NUM_STATE];
               for (int c = 0; c < NUM_STATE; c++)
               {
                  real_t s_self = 0.0, s_nbr = 0.0;
                  for (int k = 0; k < NUM_STATE; k++)
                  {
                     s_self += Tinv_can(c, k) * Q_self[k];
                     s_nbr  += Tinv_can(c, k) * Q_nbr[k];
                  }
                  Q_self_can[c] = s_self;
                  Q_nbr_can[c]  = s_nbr;
               }

               // Side-routing uses the geometric flag (matches
               // ComputeADERSharedFaceFluxRHS:4031-4034 for the +/− label).
               const real_t *src_plus  = elem1_on_plus ? Q_self_can
                                                       : Q_nbr_can;
               const real_t *src_minus = elem1_on_plus ? Q_nbr_can
                                                       : Q_self_can;

               const int dof_idx = base_dof_idx + q;
               real_t *dst_p = Q_plus_flat.data()  + dof_idx * NUM_STATE;
               real_t *dst_m = Q_minus_flat.data() + dof_idx * NUM_STATE;
               for (int c = 0; c < NUM_STATE; c++)
               {
                  dst_p[c] = src_plus[c];
                  dst_m[c] = src_minus[c];
               }
            }
         }
      }
#endif
   }
}

// ---------------------------------------------------------------------------
// ADER I-05 Phase 4: volume-integral contribution on the time-integrated state.
// ---------------------------------------------------------------------------
// Forwards to ComputeVolumeRHS, which already (a) uses the strong-form DG
// volume integrand the plan specifies and (b) ADDS to rhs rather than
// overwriting.  Defining the ADER entry point as a thin wrapper keeps the
// two routines guaranteed-consistent (same kernel, same quadrature, same
// signed Jacobians) without duplicating 60 LOC of inner loops.  The name
// is kept in the public API so Phase 6's `AdvanceADER` can call it
// explicitly — signalling intent to readers of the corrector.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ComputeADERVolumeUpdate(const Vector &I,
                                                     Vector &rhs) const
{
   MFEM_VERIFY(I.Size() == NUM_STATE * ndof_total_,
               "ComputeADERVolumeUpdate: I size "
               << I.Size() << " != NUM_STATE * ndof_total_ = "
               << NUM_STATE * ndof_total_);
   if (rhs.Size() == 0)
   {
      rhs.SetSize(NUM_STATE * ndof_total_);
      rhs = 0.0;
   }
   else
   {
      MFEM_VERIFY(rhs.Size() == NUM_STATE * ndof_total_,
                  "ComputeADERVolumeUpdate: rhs size "
                  << rhs.Size() << " != NUM_STATE * ndof_total_ = "
                  << NUM_STATE * ndof_total_);
   }

   ComputeVolumeRHS(I, rhs);
}

// ---------------------------------------------------------------------------
// Face flux: Eq. (2d) — local faces only
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ComputeFaceFluxRHS(const Vector &Q, Vector &rhs) const
{
   const real_t *Q_data = Q.GetData();

   // Round-11 Mixed-Flux dispatch flag (R-1206 short-circuit): hoisted
   // ONCE per Mult call so the per-face dispatch can short-circuit the
   // unordered_set lookup when mode == None.  Bit-identical to pre-
   // Mixed-Flux path AND zero per-face overhead in default mode.
   const bool mf_on = mf_on_;  // R-1208: cached member, kept in sync by SetMixedFluxMode

   // R-204: hoist the R-002 fault-bookkeeping guard out of the per-QP
   // loop.  If bc_.fault_attr > 0 the mesh has a fault, and the
   // interior-face loop below can dispatch to the fault branch at any
   // face whose bdr_attr matches — in that case fault_flux_ and
   // fault_dof_data_ MUST be wired.  One check per Mult call replaces
   // billions of per-QP checks with identical semantics.
   if (bc_.fault_attr > 0)
   {
      MFEM_VERIFY(fault_flux_ && fault_dof_data_,
                  "WaveOperator::ComputeFaceFluxRHS: bc_.fault_attr="
                  << bc_.fault_attr
                  << " > 0 (fault configured) but fault_flux_="
                  << (void*)fault_flux_
                  << ", fault_dof_data_=" << (void*)fault_dof_data_
                  << "; ctor did not populate fault bookkeeping.  "
                  "v9.0.0 Pelties-9 per-side flux requires both; "
                  "welded-flux fallback is no longer physical.");
   }

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
            // R-001 (v9.5.0): honor `use_precomputed_face_fluxes_` in the
            // RK4 Mult path, symmetric with ComputeADERFaceFluxRHS
            // (wave_operator.inl:2101).  Pre-fix the flag was consulted
            // only by the ADER dispatch, so drivers that mixed Mult + the
            // ADER switch silently got mixed dispatches and the Arm 1
            // boundary-lift probe was a vacuous no-op.  `bulk_bg_zero`
            // matches the linearity contract the Phase 1 probe used to
            // construct the precomputed-boundary tables (Q_bg = 0).
            if (use_precomputed_face_fluxes_)
            {
               real_t bulk_bg_zero[NUM_STATE] = {0.0};
               precomputed_face_fluxes_.AddBoundaryFaceRhs(
                  f, Q_self, bulk_bg_zero, w, shape1.GetData(),
                  ndof, dof_offset1, ndof_total_, rhs);
               continue;
            }

            FaceBC bc_type = ClassifyBoundaryFace(bdr_attr);

            switch (bc_type)
            {
               case FaceBC::Absorbing:
                  // v9.4.0 (REVIEW R-007): updated stale "Total-Q only"
                  // wording.  Callers must supply a bulk background via
                  // `WaveOperator::SetAbsorbingBackground(Q_bg)`; under
                  // fluctuation-Q dispatch, Q_bg = 0 is the natural and
                  // correct choice (AbsorbingTotal(Q_self, 0)
                  // = Interior(Q_self, 0) = Absorbing(Q_self)).
                  MFEM_VERIFY(has_bulk_bg_,
                              "wave.Mult(): SetAbsorbingBackground(Q_bg) "
                              "must be called before Absorbing BC dispatch "
                              "(Q_bg = 0 is valid under fluctuation-Q).");
                  flux_.AbsorbingTotal(nor, Q_self, bulk_bg_, F_h);
                  break;
               case FaceBC::FreeSurface:
                  // v9.4.0 (REVIEW R-007): updated stale "Total-Q only"
                  // wording.  The mode flag selects gamma-mirror vs
                  // Godunov-projection; both variants honor the bulk
                  // background (Q_bg = 0 under fluctuation-Q).
                  MFEM_VERIFY(has_bulk_bg_,
                              "wave.Mult(): SetAbsorbingBackground(Q_bg) "
                              "must be called before FreeSurface BC "
                              "dispatch (Q_bg = 0 is valid).");
                  if (free_surface_bc_mode_ == FreeSurfaceBCMode::Godunov)
                  {
                     flux_.FreeSurfaceGodunovTotal(nor, Q_self,
                                                   bulk_bg_, F_h);
                  }
                  else
                  {
                     flux_.FreeSurfaceTotal(nor, Q_self, bulk_bg_, F_h);
                  }
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
                  // Default fallback mirrors the FaceBC::Absorbing branch.
                  MFEM_VERIFY(has_bulk_bg_,
                              "wave.Mult(): SetAbsorbingBackground(Q_bg) "
                              "must be called before the default BC "
                              "dispatch (Q_bg = 0 is valid).");
                  flux_.AbsorbingTotal(nor, Q_self, bulk_bg_, F_h);
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

            if (is_fault)
            {
               // R-204: hoisted check at the top of ComputeFaceFluxRHS
               // already verified fault_flux_ && fault_dof_data_ when
               // bc_.fault_attr > 0.  Keep a debug-only assertion here
               // as a tripwire for a future refactor that bypasses the
               // hoisted check.
               MFEM_ASSERT(fault_flux_ && fault_dof_data_,
                           "R-204: bookkeeping hoisted check should have "
                           "fired at the top of ComputeFaceFluxRHS.  If "
                           "execution reaches here with null fault_flux_, "
                           "the hoisted guard was removed or bypassed.");
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

                  // R-101: per-face geometric flag (computed at ctor
                  // from element/face centroid projection on ref_normal),
                  // FP-stable across all QPs of this face.  Replaces the
                  // per-QP `!qpd.sign_flipped` derivation, which on
                  // near-axis-aligned faces could flip sign between QPs
                  // due to FP noise in the CalcOrtho dot product and
                  // produced y-mirror non-invariance of the per-side
                  // assembly.
                  MFEM_ASSERT(fb_idx >= 0 &&
                              fb_idx < static_cast<int>(
                                 interior_fault_elem1_on_plus_.size()),
                              "R-101: interior_fault_elem1_on_plus_ "
                              "missing entry for fb_idx=" << fb_idx);
                  const bool elem1_on_plus =
                     interior_fault_elem1_on_plus_[fb_idx];

                  // 1. Reconstruct canonical (pre-Step-5) BP5 frame.
                  // BP5 convention: can_t1 = dip, can_t2 = strike.
                  // Use the per-face flag (R-101) to decide whether to
                  // negate qpd.{normal,tangent1,tangent2} so can_n
                  // points along ref_normal.  Equivalent in exact
                  // arithmetic to `qpd.sign_flipped ? -x : x` but
                  // FP-stable per face (no per-QP sign bimodality).
                  const FaultBasisQPData &qpd = *qpd_ptr;
                  const bool should_negate_frame = !elem1_on_plus;
                  real_t can_n[3], can_t1[3], can_t2[3];
                  for (int d = 0; d < 3; d++)
                  {
                     can_n[d]  = should_negate_frame ? -qpd.normal[d]
                                                     :  qpd.normal[d];
                     can_t1[d] = should_negate_frame ? -qpd.tangent1[d]
                                                     :  qpd.tangent1[d];
                     can_t2[d] = should_negate_frame ? -qpd.tangent2[d]
                                                     :  qpd.tangent2[d];
                  }

                  DenseMatrix T_can(NUM_STATE), Tinv_can(NUM_STATE);
                  GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);
                  GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2,
                                                    Tinv_can);

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

                  // 2. v9.4.0 Commit 3: fluctuation-Q dispatch.  Bulk Q
                  //    carries only the dynamic fluctuation (Q = 0 at
                  //    rest); static pre-stress and nucleation live in
                  //    DOFData (tau*_0, tau*_nuc).  FaultFaceFlux::
                  //    Evaluate sums tau*_total = tau*_0 + tau*_nuc +
                  //    tau*_trial(Q) internally (v9.4.0 Commit 1).
                  //    Fault dispatch does NOT read Q_bg (REVIEW R-003);
                  //    Q_bg is consumed only by Absorbing/FreeSurface/
                  //    PML BC branches, which carry their own guards.
                  real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
                  fault_flux_->Evaluate(fdata,
                                        Q_plus_local, Q_minus_local,
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

#ifdef SEAS_DIAG_FAULT_FLUX
                  // C-1s STRESS-ROT (interior-fault path): rotation-
                  // pipeline trace at the hypocenter QP.  Records the
                  // 9-component Q on each canonical side at four
                  // pipeline stages so the C-1n NORMAL stress-channel
                  // leak (debug doc 2026-04-25) can be back-traced:
                  //   GLOB-IN:  bulk Q the wave op reads (global)
                  //   LOC-IN:   after Tinv*Q (fault-local; Evaluate input)
                  //   LOC-IMP:  after BuildImposedState (local)
                  //   GLOB-OUT: after T*Q_imp (global; injected to bulk)
                  // BASIS line dumps can_n/can_t1/can_t2 + sign-flip so
                  // we can verify the rotation matrix per-QP.  Gated by
                  // data.diag_print → only fires at the tagged hypo QP.
                  if (dof_idx >= 0 &&
                      dof_idx < static_cast<int>(fault_dof_data_->size()) &&
                      (*fault_dof_data_)[dof_idx].diag_print)
                  {
                     std::fprintf(stderr,
                        "[C-1s INT BASIS] rank=%d  "
                        "can_n=(%+.4e,%+.4e,%+.4e)  "
                        "can_t1=(%+.4e,%+.4e,%+.4e)  "
                        "can_t2=(%+.4e,%+.4e,%+.4e)  "
                        "sign_flipped=%d  elem1_on_plus=%d\n",
                        g_seas_my_rank,
                        can_n[0], can_n[1], can_n[2],
                        can_t1[0], can_t1[1], can_t1[2],
                        can_t2[0], can_t2[1], can_t2[2],
                        qpd.sign_flipped ? 1 : 0,
                        elem1_on_plus  ? 1 : 0);

                     auto _c1s_print = [&](const char *tag, const real_t *q)
                     {
                        std::fprintf(stderr,
                           "[C-1s INT %s] rank=%d  "
                           "SXX=%+.4e SYY=%+.4e SZZ=%+.4e  "
                           "SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
                           "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
                           tag, g_seas_my_rank,
                           q[SXX], q[SYY], q[SZZ],
                           q[SXY], q[SYZ], q[SXZ],
                           q[VX],  q[VY],  q[VZ]);
                     };

                     // Self/Nbr → canonical ± mapping for global-frame Q.
                     const real_t *Q_glob_p = elem1_on_plus ? Q_self : Q_nbr;
                     const real_t *Q_glob_m = elem1_on_plus ? Q_nbr  : Q_self;

                     _c1s_print("GLOB-IN+",  Q_glob_p);
                     _c1s_print("GLOB-IN-",  Q_glob_m);
                     _c1s_print("LOC-IN+",   Q_plus_local);
                     _c1s_print("LOC-IN-",   Q_minus_local);
                     _c1s_print("LOC-IMP+",  Q_imp_plus);
                     _c1s_print("LOC-IMP-",  Q_imp_minus);
                     _c1s_print("GLOB-OUT+", Q_imp_plus_g);
                     _c1s_print("GLOB-OUT-", Q_imp_minus_g);
                  }
#endif

                  // v9.0.0 Pelties-9 per-side flux (plan §14.2).  Each
                  // side's bulk rhs gets its OWN imposed-state flux
                  // `A_{can_n} . Q_imp_side` in the global frame,
                  // realised via the identity
                  //   flux_.Interior(n, Q, Q) = T . (A_x^+ + A_x^-) . T^{-1} . Q
                  //                           = T . A_x . T^{-1} . Q
                  //                           = A_n . Q
                  // Using `can_n` (canonical, rank-invariant) rather
                  // than MFEM's local `nor` removes the L/R routing
                  // step and matches the shared-fault convention.
                  real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
                  flux_.Interior(can_n, Q_imp_plus_g,  Q_imp_plus_g,
                                 F_h_plus);
                  flux_.Interior(can_n, Q_imp_minus_g, Q_imp_minus_g,
                                 F_h_minus);

#ifdef SEAS_DIAG_FAULT_FLUX
                  // C-2 FLUX: print Elem1's own-side flux in the GLOBAL
                  // frame (the flux the DG rhs update actually sees for
                  // this element).
                  const real_t *F_h_elem1 = elem1_on_plus ? F_h_plus
                                                          : F_h_minus;
                  if (dof_idx >= 0 &&
                      dof_idx < static_cast<int>(fault_dof_data_->size()) &&
                      (*fault_dof_data_)[dof_idx].diag_print)
                  {
                     real_t F_v_mag = std::sqrt(
                        F_h_elem1[VX]*F_h_elem1[VX]
                      + F_h_elem1[VY]*F_h_elem1[VY]
                      + F_h_elem1[VZ]*F_h_elem1[VZ]);
                     real_t F_s_max = 0.0;
                     for (int c = 0; c < 6; c++)
                     {
                        F_s_max = std::max(F_s_max, std::abs(F_h_elem1[c]));
                     }
                     std::fprintf(stderr,
                        "[C-2 FLUX] rank=%d  dof=%d  side=%c  "
                        "|F_v|=%.3e m2/s2  max|F_stress|=%.3e Pa*m/s  "
                        "F_h[VX]=%+.3e  F_h[SXY]=%+.3e  F_h[VY]=%+.3e  "
                        "F_h[SXZ]=%+.3e\n",
                        g_seas_my_rank, dof_idx,
                        elem1_on_plus ? '+' : '-',
                        F_v_mag, F_s_max,
                        F_h_elem1[VX], F_h_elem1[SXY],
                        F_h_elem1[VY], F_h_elem1[SXZ]);
                  }
#endif

#ifdef SEAS_DIAG_FAULT_FLUX
                  // C-3 RHS: bracket the Elem1 accumulation to measure
                  // the delta injected into rhs[VX, elem1, dof0].
                  const bool _c3_diag =
                     (dof_idx >= 0 &&
                      dof_idx < static_cast<int>(fault_dof_data_->size()) &&
                      (*fault_dof_data_)[dof_idx].diag_print);
                  const int  _c3_probe =
                     VX * ndof_total_ + dof_offset1 + 0;
                  const real_t _c3_pre = _c3_diag ? rhs[_c3_probe] : 0.0;
#endif
                  // Per-side DG assembly.  Each element's rhs gets its
                  // own side's flux, signed by its outward normal
                  // relative to can_n:
                  //   plus-side  elem (outward = +can_n): rhs -= w*shape*F_h_plus
                  //   minus-side elem (outward = -can_n): rhs += w*shape*F_h_minus
                  if (elem1_on_plus)
                  {
                     for (int c = 0; c < NUM_STATE; c++)
                     {
                        for (int i = 0; i < ndof; i++)
                        {
                           rhs[c * ndof_total_ + dof_offset1 + i] -=
                              w * shape1(i) * F_h_plus[c];
                           rhs[c * ndof_total_ + dof_offset2 + i] +=
                              w * shape2(i) * F_h_minus[c];
                        }
                     }
                  }
                  else
                  {
                     for (int c = 0; c < NUM_STATE; c++)
                     {
                        for (int i = 0; i < ndof; i++)
                        {
                           rhs[c * ndof_total_ + dof_offset1 + i] +=
                              w * shape1(i) * F_h_minus[c];
                           rhs[c * ndof_total_ + dof_offset2 + i] -=
                              w * shape2(i) * F_h_plus[c];
                        }
                     }
                  }

#ifdef SEAS_DIAG_FAULT_FLUX
                  if (_c3_diag)
                  {
                     const real_t _c3_post = rhs[_c3_probe];
                     std::fprintf(stderr,
                        "[C-3 RHS]  rank=%d  dof=%d  side=%c  "
                        "rhs[VX,elem1,dof0] pre=%+.3e post=%+.3e  "
                        "delta=%+.3e  w=%.3e  shape1(0)=%.3e  "
                        "F_h_elem1[VX]=%+.3e\n",
                        g_seas_my_rank, dof_idx,
                        elem1_on_plus ? '+' : '-',
                        _c3_pre, _c3_post, _c3_post - _c3_pre, w,
                        shape1(0), F_h_elem1[VX]);
                  }
#endif
               }
               else
               {
                  // v9.0.0 Pelties-9 (plan §14.2 + R-F02): a fault face
                  // without a valid DOFData / FaultBasis mapping is a
                  // ctor-population bug.  The previous welded-flux
                  // fallback silently bypassed Pelties-9 eq. (9) and
                  // under-radiated the affected QP.  Force the bug to
                  // surface instead of masking it.
                  MFEM_ABORT("interior fault face f=" << f
                             << " has no FaultBasis/DOFData mapping"
                             << " (fb_idx=" << fb_idx
                             << ", dof_idx=" << dof_idx
                             << ", have_basis=" << have_basis << "). "
                             << "v9.0.0 Pelties-9 per-side flux requires"
                             << " a valid mapping; welded-flux fallback"
                             << " is no longer physical.");
               }
            }
            else
            {
               // R-001 (v9.5.0): honor `use_precomputed_face_fluxes_` in
               // the RK4 Mult interior non-fault path.  Symmetric with
               // ComputeADERFaceFluxRHS (wave_operator.inl:2117); both
               // sides are deposited via AddInteriorFaceRhs, keyed by
               // (face, Elem1) and (face, Elem2) respectively.
               if (use_precomputed_face_fluxes_)
               {
                  precomputed_face_fluxes_.AddInteriorFaceRhs(
                     f, /*caller_elem=*/e1, Q_self, Q_nbr, w,
                     shape1.GetData(), ndof, dof_offset1, ndof_total_,
                     rhs);
                  precomputed_face_fluxes_.AddInteriorFaceRhs(
                     f, /*caller_elem=*/e2, Q_nbr, Q_self, w,
                     shape2.GetData(), ndof, dof_offset2, ndof_total_,
                     rhs);
               }
               else
               {
                  // Round-11 Mixed-Flux dispatch (R-1206 short-circuit):
                  // when mixed-flux mode is None (default), the mf_on
                  // boolean is false and `count(f)` is never evaluated —
                  // bit-identical AND zero-cost vs pre-Mixed-Flux path.
                  if (mf_on && central_flux_face_set_.count(f) > 0)
                  {
                     flux_.Central(nor, Q_self, Q_nbr, F_h);
                  }
                  else
                  {
                     // Regular interior face: standard Godunov flux
                     flux_.Interior(nor, Q_self, Q_nbr, F_h);
                  }

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

#ifdef SEAS_DIAG_FAULT_FLUX
                  // C-2B NONFAULT-FACE: per-flux dump on the 6 non-fault
                  // interior faces of the hypocenter's two adjacent tets.
                  // Captures Q_self, Q_nbr, F_h and the per-DOF
                  // contribution at the diag DOF (when one of the two
                  // diag elements participates as Elem1 or Elem2).
                  // Fires only on the diag rank (others have empty list)
                  // and only on the FIRST QP per face per Mult call (q==0,
                  // user-decision: per-Mult, not per-stage / per-QP).
                  if (q == 0 && !diag_nonfault_faces_.empty())
                  {
                     bool match = false;
                     for (size_t kk = 0;
                          kk < diag_nonfault_faces_.size(); kk++)
                     {
                        if (diag_nonfault_faces_[kk] == f) { match = true; break; }
                     }
                     if (match)
                     {
                        // Determine which side (+ or -) this face touches.
                        const char *side_p = (e1 == diag_elem_plus_  ||
                                              e2 == diag_elem_plus_)  ? "P" : "";
                        const char *side_m = (e1 == diag_elem_minus_ ||
                                              e2 == diag_elem_minus_) ? "M" : "";
                        std::fprintf(stderr,
                           "[C-2B NONFAULT-FACE] rank=%d face=%d e1=%d e2=%d "
                           "side=%s%s nor=(%+.4e,%+.4e,%+.4e) w=%.4e\n",
                           g_seas_my_rank, f, e1, e2, side_p, side_m,
                           nor[0], nor[1], nor[2], w);
                        std::fprintf(stderr,
                           "[C-2B NONFAULT-FACE Q_self] rank=%d face=%d "
                           "SXX=%+.4e SYY=%+.4e SZZ=%+.4e "
                           "SXY=%+.4e SYZ=%+.4e SXZ=%+.4e "
                           "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
                           g_seas_my_rank, f,
                           Q_self[SXX], Q_self[SYY], Q_self[SZZ],
                           Q_self[SXY], Q_self[SYZ], Q_self[SXZ],
                           Q_self[VX],  Q_self[VY],  Q_self[VZ]);
                        std::fprintf(stderr,
                           "[C-2B NONFAULT-FACE Q_nbr]  rank=%d face=%d "
                           "SXX=%+.4e SYY=%+.4e SZZ=%+.4e "
                           "SXY=%+.4e SYZ=%+.4e SXZ=%+.4e "
                           "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
                           g_seas_my_rank, f,
                           Q_nbr[SXX], Q_nbr[SYY], Q_nbr[SZZ],
                           Q_nbr[SXY], Q_nbr[SYZ], Q_nbr[SXZ],
                           Q_nbr[VX],  Q_nbr[VY],  Q_nbr[VZ]);
                        std::fprintf(stderr,
                           "[C-2B NONFAULT-FACE F_h]    rank=%d face=%d "
                           "SXX=%+.4e SYY=%+.4e SZZ=%+.4e "
                           "SXY=%+.4e SYZ=%+.4e SXZ=%+.4e "
                           "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
                           g_seas_my_rank, f,
                           F_h[SXX], F_h[SYY], F_h[SZZ],
                           F_h[SXY], F_h[SYZ], F_h[SXZ],
                           F_h[VX],  F_h[VY],  F_h[VZ]);
                     }
                  }
#endif
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

      // Round-11 Mixed-Flux short-circuit (R-1206): hoist once per call.
      const bool mf_on = mf_on_;  // R-1208: cached member, kept in sync by SetMixedFluxMode

      auto &pmesh = static_cast<const ParMesh &>(mesh_);
      int n_shared = pmesh.GetNSharedFaces();
      if (n_shared == 0) { return; }

      // R-204 (symmetric with ComputeFaceFluxRHS): hoist the fault
      // bookkeeping guard out of the per-face / per-QP loops.  If the
      // mesh has a fault (bc_.fault_attr > 0) then fault_flux_ and
      // fault_dof_data_ MUST be wired — the per-QP `fault_active`
      // short-circuit at line ~1063 would otherwise silently skip the
      // fault path on this rank, under-radiating the same way the
      // pre-R-002 welded-flux fallback did.
      if (bc_.fault_attr > 0)
      {
         MFEM_VERIFY(fault_flux_ && fault_dof_data_,
                     "WaveOperator::ComputeSharedFaceFluxRHS: "
                     "bc_.fault_attr=" << bc_.fault_attr
                     << " > 0 (fault configured) but fault_flux_="
                     << (void*)fault_flux_
                     << ", fault_dof_data_=" << (void*)fault_dof_data_
                     << "; ctor did not populate fault bookkeeping.  "
                     "v9.0.0 Pelties-9 per-side flux requires both on "
                     "every rank that owns a shared fault face.");
      }

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

         // R-001 (v9.5.0): mesh face index for the precomputed-table
         // lookup.  `sf` is the shared-face index; the precomputed-face
         // entry key is (mesh_face_idx, Elem1No).  Symmetric with
         // ComputeADERSharedFaceFluxRHS (wave_operator.inl:2308).
         const int mesh_face_idx =
            const_cast<ParMesh &>(pmesh).GetSharedFace(sf);

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
                  // R-003: post-v9.0.0 per-side flux, the failure modes for
                  // a wrong `elem1_on_plus` are NOT symmetric (the pre-fix
                  // welded flux partially masked them via cancellation).
                  // A silent fallback to `false` would simultaneously
                  // mis-select Q_imp_± and flip assemble_sign.  Fail loud.
                  MFEM_VERIFY(sf_idx_in_fault >= 0 &&
                              sf_idx_in_fault <
                              static_cast<int>(shared_fault_elem1_on_plus_.size()),
                              "shared_fault_elem1_on_plus_ missing entry for "
                              "sf_idx_in_fault=" << sf_idx_in_fault
                              << " (array size="
                              << shared_fault_elem1_on_plus_.size()
                              << "; fault_shared_faces_ size="
                              << fault_shared_faces_.Size() << "); "
                              "ctor did not populate. v9.0.0 Pelties-9 "
                              "per-side flux cannot default-route silently.");
                  const bool elem1_on_plus =
                     shared_fault_elem1_on_plus_[sf_idx_in_fault];

                  const FaultBasisData &bd = fault_basis_->GetBasis(basis_idx);
                  // R-202: use MFEM_VERIFY (not MFEM_ASSERT).  MFEM_ASSERT
                  // is a no-op in release builds; the next line would then
                  // read bd.qp_data[q] out-of-bounds — undefined behaviour.
                  // qpd.normal/tangent1/tangent2 are load-bearing for the
                  // per-side flux and cross-rank consistency, so garbage
                  // here silently corrupts the result on one rank.
                  MFEM_VERIFY(q < static_cast<int>(bd.qp_data.size()),
                              "FaultBasis::qp_data not populated for shared "
                              "fault face — ComputeQPBasisShared missed "
                              "this face (basis_idx=" << basis_idx
                              << ", q=" << q
                              << ", qp_data.size()=" << bd.qp_data.size()
                              << ").  v9.0.0 Pelties-9 per-side flux reads "
                              "qpd.normal/tangent1/tangent2 from this entry; "
                              "a wrong can_n corrupts the flux.");
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
                  // v9.4.0 Commit 3: fluctuation-Q dispatch.  REVIEW
                  // R-003: no has_bulk_bg_ guard — fault dispatch does
                  // not consume Q_bg.
                  fault_flux_->Evaluate(fdata,
                                        Q_plus_local, Q_minus_local,
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

#ifdef SEAS_DIAG_FAULT_FLUX
                  // C-1s STRESS-ROT (shared-fault path): same trace as
                  // the interior-fault block above, but tagged "SHR" so
                  // we can tell whether the hypocenter QP sits on a
                  // rank-internal interior fault face or a rank-boundary
                  // shared fault face.  Both paths share the rotation
                  // pipeline; only one fires per QP per macro step.
                  if (dof_idx >= 0 &&
                      dof_idx < static_cast<int>(fault_dof_data_->size()) &&
                      (*fault_dof_data_)[dof_idx].diag_print)
                  {
                     std::fprintf(stderr,
                        "[C-1s SHR BASIS] rank=%d  "
                        "can_n=(%+.4e,%+.4e,%+.4e)  "
                        "can_t1=(%+.4e,%+.4e,%+.4e)  "
                        "can_t2=(%+.4e,%+.4e,%+.4e)  "
                        "sign_flipped=%d  elem1_on_plus=%d\n",
                        g_seas_my_rank,
                        can_n[0], can_n[1], can_n[2],
                        can_t1[0], can_t1[1], can_t1[2],
                        can_t2[0], can_t2[1], can_t2[2],
                        qpd.sign_flipped ? 1 : 0,
                        elem1_on_plus  ? 1 : 0);

                     auto _c1s_print = [&](const char *tag, const real_t *q)
                     {
                        std::fprintf(stderr,
                           "[C-1s SHR %s] rank=%d  "
                           "SXX=%+.4e SYY=%+.4e SZZ=%+.4e  "
                           "SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
                           "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
                           tag, g_seas_my_rank,
                           q[SXX], q[SYY], q[SZZ],
                           q[SXY], q[SYZ], q[SXZ],
                           q[VX],  q[VY],  q[VZ]);
                     };

                     const real_t *Q_glob_p = elem1_on_plus ? Q_self : Q_nbr;
                     const real_t *Q_glob_m = elem1_on_plus ? Q_nbr  : Q_self;

                     _c1s_print("GLOB-IN+",  Q_glob_p);
                     _c1s_print("GLOB-IN-",  Q_glob_m);
                     _c1s_print("LOC-IN+",   Q_plus_local);
                     _c1s_print("LOC-IN-",   Q_minus_local);
                     _c1s_print("LOC-IMP+",  Q_imp_plus);
                     _c1s_print("LOC-IMP-",  Q_imp_minus);
                     _c1s_print("GLOB-OUT+", Q_imp_plus_g);
                     _c1s_print("GLOB-OUT-", Q_imp_minus_g);
                  }
#endif

                  // R-802 fix: use canonical normal `can_n` (bit-identical
                  // on both ranks, by construction) for the Godunov flux so
                  // both ranks compute the SAME F_h regardless of whether
                  // MFEM's CalcOrtho gives `nor_A = -nor_B` (classical) or
                  // `nor_A = nor_B` (observed empirically on the inline
                  // v9.0.0 Pelties-9 per-side flux (plan §14.3).  This
                  // rank owns Elem1 only; it assembles that side's own
                  // `A_{can_n} . Q_imp_side` contribution.  The paired
                  // rank (with elem1_on_plus flipped) assembles the
                  // opposite-side contribution separately with the
                  // opposite sign.
                  //
                  // Sign convention:
                  //   plus-side  elem1 (outward = +can_n) :  rhs -= w*shape*F
                  //   minus-side elem1 (outward = -can_n) :  rhs += w*shape*F
                  real_t F_h_side[NUM_STATE];
                  const real_t *Q_imp_side = elem1_on_plus
                                             ? Q_imp_plus_g
                                             : Q_imp_minus_g;
                  flux_.Interior(can_n, Q_imp_side, Q_imp_side, F_h_side);

                  const real_t assemble_sign = elem1_on_plus ? -1.0 : +1.0;
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     for (int i = 0; i < ndof; i++)
                     {
                        rhs[c * ndof_total_ + dof_offset1 + i] +=
                           assemble_sign * w * shape1(i) * F_h_side[c];
                     }
                  }
                  continue;  // Fault accumulation done; skip the generic
                             // -= F_h path below (which assumes rank-local
                             // `nor` + standard Godunov orientation).
               }
               else
               {
                  // v9.0.0 Pelties-9 (plan §14.3 + R-F02): see §14.2's
                  // matching interior-fault fallback for rationale.
                  // Welded-flux fallback at a SHARED fault face is
                  // equally unacceptable post-fix; hard-abort.
                  MFEM_ABORT("shared fault face sf=" << sf
                             << " has no DOFData mapping (dof_idx="
                             << dof_idx << "). v9.0.0 Pelties-9 per-side"
                             << " flux requires a valid mapping;"
                             << " welded-flux fallback is no longer"
                             << " physical.");
               }
            }
            else if (use_precomputed_face_fluxes_)
            {
               // R-001 (v9.5.0): one-sided precomputed interior dispatch
               // for shared non-fault faces, symmetric with
               // ComputeADERSharedFaceFluxRHS (wave_operator.inl:2459).
               // Uses the LEGACY AddInteriorFaceRhs path because the
               // neighbor DOFs live in the ghost buffer; the Arm 3d
               // Full-path's dof_offset_nbr assumes local Q indexing.
               precomputed_face_fluxes_.AddInteriorFaceRhs(
                  mesh_face_idx, /*caller_elem=*/e1,
                  Q_self, Q_nbr, w, shape1.GetData(),
                  ndof, dof_offset1, ndof_total_, rhs);
               continue;
            }
            else
            {
               // Non-fault shared face: standard welded Godunov flux.
               //
               // I-04 invariant (addresses REVIEW.md R-I04-005): a shared
               // face in MFEM's ParMesh model always has elements on BOTH
               // sides and therefore carries no boundary attribute.  A
               // free-surface face is boundary-attributed and must live
               // on a non-shared face; it is dispatched exclusively by
               // ComputeFaceFluxRHS's FaceBC::FreeSurface branch, which
               // honours `free_surface_bc_mode_`.  If a future driver or
               // mesh convention breaks this invariant (shared face that
               // also carries a free-surface attribute), the BC-mode flag
               // will be silently ignored here — add a symmetric branch
               // or a loud MFEM_VERIFY at that time.
               //
               // Round-11 Mixed-Flux dispatch (R-1206 short-circuit + R-1207
               // shared-seam membership): when mf_on, check the central
               // set keyed by global mesh face index.
               if (mf_on && central_flux_face_set_.count(mesh_face_idx) > 0)
               {
                  flux_.Central(nor, Q_self, Q_nbr, F_h);
               }
               else
               {
                  flux_.Interior(nor, Q_self, Q_nbr, F_h);
               }
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
// ADER I-05 Phase 6: ComputeADERFaceFluxRHS — local faces.
// ---------------------------------------------------------------------------
// Duplicates ComputeFaceFluxRHS's control flow with the following diffs:
//   * Input is `I` (time-integrated state) and `dt`; not `Q`.
//   * Non-fault branches apply the standard Godunov flux treating `I` as
//     a linear state (F_h(I) = dt · F_h(Q_avg) by linearity).  For
//     background-aware BCs (AbsorbingTotal, FreeSurface*Total), the
//     background is scaled by `dt` so the resulting flux equals
//     `dt · F_h(Q_avg, Q_bg)` — the correctly time-integrated boundary
//     flux.
//   * Fault faces dispatch to `fault_flux_->EvaluateADER` (v9.4.0
//     Commit 3: fluctuation-Q dispatch) which handles the 1/dt → dt
//     scaling internally (see Phase 5).
//
// The duplication (vs. refactoring ComputeFaceFluxRHS into a shared core)
// is intentional: `wave_operator.inl` is on the CLAUDE.md "Files Requiring
// Extreme Care" list and the RK4 path must remain byte-identical.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ComputeADERFaceFluxRHS(const Vector &I,
                                                    real_t dt,
                                                    Vector &rhs) const
{
   MFEM_VERIFY(dt > 0.0,
               "ComputeADERFaceFluxRHS: dt must be > 0, got " << dt);
   MFEM_VERIFY(I.Size() == NUM_STATE * ndof_total_,
               "ComputeADERFaceFluxRHS: I size mismatch");

   // v9.4.0: Q_bg feeds Absorbing / FreeSurface / PML BC branches and
   // the PML damp-toward-bg loop.  Q_bg = 0 is valid under fluctuation-Q
   // dispatch, but SetAbsorbingBackground() must still be called so
   // downstream reads of bulk_bg_ are from an initialised buffer rather
   // than the default-ctor zero-fill (defensive gate, kept for explicit
   // contract).  REVIEW R-007: updated stale "Total-Q only" wording.
   //
   // R-1505: COLLECTIVE consensus check.  A rank-local MFEM_VERIFY would
   // abort one rank while the others proceeded into the next collective
   // (q_gf.ExchangeFaceNbrData below in the shared corrector), causing a
   // deadlock instead of a fail-loud abort across all ranks.  Mirror the
   // SetMixedFluxMode pattern (wave_operator.inl:1210-1226): Allreduce-MIN
   // of has_bulk_bg_ and abort everywhere if ANY rank lacks it.
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto &pmesh_consensus = static_cast<const ParMesh &>(mesh_);
      const int my_has = has_bulk_bg_ ? 1 : 0;
      int min_has = 0;
      MPI_Allreduce(&my_has, &min_has, 1, MPI_INT, MPI_MIN,
                    pmesh_consensus.GetComm());
      MFEM_VERIFY(min_has == 1,
                  "ComputeADERFaceFluxRHS: SetAbsorbingBackground(Q_bg) was "
                  "NOT called on every rank (min has_bulk_bg_=" << min_has
                  << ").  A rank-local missing call would deadlock at the "
                  "shared corrector's ExchangeFaceNbrData; the collective "
                  "check fails loud everywhere instead.  See R-1505.");
#endif
   }
   else
   {
      MFEM_VERIFY(has_bulk_bg_,
                  "wave.AdvanceADER() requires SetAbsorbingBackground(Q_bg) "
                  "to have been called (Q_bg = 0 is valid under "
                  "fluctuation-Q dispatch).");
   }

   const real_t *I_data = I.GetData();

   // Time-integrated background: bulk_bg_scaled[c] = dt · bulk_bg_[c].
   real_t bulk_bg_scaled[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      bulk_bg_scaled[c] = dt * bulk_bg_[c];
   }

   // Round-11 Mixed-Flux short-circuit (R-1206): hoist once per call.
   const bool mf_on = mf_on_;  // R-1208: cached member, kept in sync by SetMixedFluxMode

   if (bc_.fault_attr > 0)
   {
      MFEM_VERIFY(fault_flux_ && fault_dof_data_,
                  "WaveOperator::ComputeADERFaceFluxRHS: bc_.fault_attr="
                  << bc_.fault_attr << " > 0 but fault_flux_ or "
                  "fault_dof_data_ is null.");
   }

   // Round-12 Patch 2: per-face stage-averaging hook (test-only).  When
   // SEAS_TEST_EVAL_FACE_AVG_STAGE is set to "trial" / "theta" / "vabs"
   // / "tcorr", the fault branch below batches all QPs on a face,
   // averages the selected `FaultFaceFlux::EvalStageState` field across
   // the face, then completes the downstream chain.  Unset (or
   // "") → None → existing per-QP path is used verbatim, byte-identical
   // to the pre-patch code.  This hook MUST live in the face-face loop
   // (not inside `FaultFaceFlux::Evaluate`) because averaging requires
   // simultaneous access to every QP on a face — `Evaluate` only sees
   // one.
   enum class FaultEvalStageAvgMode { None, Trial, Theta, Vabs, Tcorr };
   FaultEvalStageAvgMode eval_avg_mode = FaultEvalStageAvgMode::None;
   {
      const char *env = std::getenv("SEAS_TEST_EVAL_FACE_AVG_STAGE");
      if (env && env[0] != '\0')
      {
         if      (std::strcmp(env, "trial") == 0) { eval_avg_mode = FaultEvalStageAvgMode::Trial; }
         else if (std::strcmp(env, "theta") == 0) { eval_avg_mode = FaultEvalStageAvgMode::Theta; }
         else if (std::strcmp(env, "vabs")  == 0) { eval_avg_mode = FaultEvalStageAvgMode::Vabs;  }
         else if (std::strcmp(env, "tcorr") == 0) { eval_avg_mode = FaultEvalStageAvgMode::Tcorr; }
         else
         {
            MFEM_VERIFY(false,
                        "ComputeADERFaceFluxRHS: unknown "
                        "SEAS_TEST_EVAL_FACE_AVG_STAGE=\"" << env
                        << "\" (expected one of: trial, theta, vabs, tcorr).");
         }
      }
   }

   // R-1203: SEAS_TEST_NONFAULT_BOTH_SYM is mathematically meaningless.
   // For ANY conservative flux (Interior, Central) the anti-symmetry
   // identity F(+n,L,R) + F(-n,R,L) = 0 holds exactly (Round-2 Gate-3
   // verified Central to 1e-8 relative; Interior to 1e-15).  So
   // `0.5*(F(nor, L, R) + F(-nor, R, L))` ZEROES every non-fault face's
   // contribution — a complete physics break — and silently masks
   // any Mixed-Flux dispatch error at the central faces.  Hook is
   // permanently disabled.  Setting the env var is now an immediate
   // abort: the operator stops the run rather than letting the user
   // believe they are auditing physics.
   const bool nonfault_both_sym = false;
   if (std::getenv("SEAS_TEST_NONFAULT_BOTH_SYM"))
   {
      MFEM_ABORT("SEAS_TEST_NONFAULT_BOTH_SYM is permanently disabled "
                 "(R-1203).  The (n,L,R) -> (-n,R,L) symmetrization is "
                 "identically zero for any conservative flux including "
                 "Central, so this hook ZEROED the bulk wave equation "
                 "rather than 'auditing' it.  Drop the env var.  If a "
                 "rotation-pipeline diagnostic is genuinely needed, use "
                 "the (n)->(-n) variant WITHOUT the L<->R swap "
                 "(it tests rotation parity instead of conservation).");
   }

   for (int f = 0; f < mesh_.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh_.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      int e1 = ftr->Elem1No;
      int e2 = ftr->Elem2No;

      int bdr_attr = face_bdr_attr_[f];

      if (e2 < 0 && shared_mesh_face_set_.count(f) > 0) { continue; }

      bool is_boundary = (e2 < 0) && (bdr_attr > 0);

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

         real_t I_self[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               I_self[c] += shape1(i) * I_data[c * ndof_total_ + dof_offset1 + i];
            }
         }

         real_t F_h[NUM_STATE];

         // === §6.2.hoist (R5-001 FIX) ================================
         // The pre-patch ComputeADERFaceFluxRHS declared dof_offset2,
         // shape2, I_nbr, and is_fault INSIDE its `else (interior)`
         // branch.  §6.2's new dispatch ordering must branch on
         // is_fault as the OUTERMOST predicate, so we hoist these four
         // declarations to outer scope.  The `!is_boundary` guard
         // ensures boundary faces never evaluate e2 * ndof_per_el_
         // with e2 < 0 (which would produce a negative offset and
         // out-of-bounds reads).
         int    dof_offset2 = -1;
         Vector shape2;
         real_t I_nbr[NUM_STATE] = {0};
         bool   is_fault = false;

         if (!is_boundary)
         {
            const FiniteElement *fe2 = fes_->GetFE(e2);
            dof_offset2 = e2 * ndof_per_el_;
            shape2.SetSize(ndof);
            IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2);
            fe2->CalcShape(ip2, shape2);
            for (int c = 0; c < NUM_STATE; c++)
            {
               I_nbr[c] = 0.0;
               for (int i = 0; i < ndof; i++)
               {
                  I_nbr[c] += shape2(i) * I_data[c * ndof_total_
                                                  + dof_offset2 + i];
               }
            }
            is_fault = (bdr_attr == bc_.fault_attr) && (bc_.fault_attr > 0);
         }
         // is_fault remains `false` on boundary faces by construction of
         // the `!is_boundary` guard.  A 1-sided face tagged with
         // bc_.fault_attr is a mesh-configuration error caught below in
         // the boundary dispatch (ClassifyBoundaryFace → FaceBC::Fault →
         // MFEM_ABORT).
         // === end §6.2.hoist =========================================

         if (is_fault)
         {
            // R2-008 FIX — fault-branch body lifted verbatim from the
            // pre-patch ComputeADERFaceFluxRHS interior-else branch.
            // Every line is load-bearing (BP5 canonical frame +
            // FaultFaceFlux::EvaluateADER + per-side accumulation).
            MFEM_ASSERT(fault_flux_ && fault_dof_data_,
                        "R-204-ader: bookkeeping hoisted check should "
                        "have fired at the top of ComputeADERFaceFluxRHS.");

            if (eval_avg_mode == FaultEvalStageAvgMode::None)
            {
               auto it = fault_face_dof_offset_.find(f);
               int dof_idx = -1;
               if (it != fault_face_dof_offset_.end())
               {
                  dof_idx = it->second + q;
               }

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

                  // R-101: per-face geometric flag (FP-stable across
                  // all QPs of this face); see ctor population around
                  // line ~395 and the corresponding Mult-path block
                  // around line ~1290.
                  MFEM_ASSERT(fb_idx >= 0 &&
                              fb_idx < static_cast<int>(
                                 interior_fault_elem1_on_plus_.size()),
                              "R-101: interior_fault_elem1_on_plus_ "
                              "missing entry for fb_idx=" << fb_idx);
                  const bool elem1_on_plus =
                     interior_fault_elem1_on_plus_[fb_idx];

                  const FaultBasisQPData &qpd = *qpd_ptr;
                  const bool should_negate_frame = !elem1_on_plus;
                  real_t can_n[3], can_t1[3], can_t2[3];
                  for (int d = 0; d < 3; d++)
                  {
                     can_n[d]  = should_negate_frame ? -qpd.normal[d]
                                                     :  qpd.normal[d];
                     can_t1[d] = should_negate_frame ? -qpd.tangent1[d]
                                                     :  qpd.tangent1[d];
                     can_t2[d] = should_negate_frame ? -qpd.tangent2[d]
                                                     :  qpd.tangent2[d];
                  }

                  DenseMatrix T_can(NUM_STATE), Tinv_can(NUM_STATE);
                  GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);
                  GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2,
                                                    Tinv_can);

#ifdef SEAS_DIAG_TPV104_FAULT_BASIS
                  // D1 instrumentation (TPV104 σ_n perturbation diagnostic).
                  // For each fault QP, print sign_flipped and the resulting
                  // canonical basis once.  A planar TPV104 fault SHOULD have
                  // sign_flipped uniform across all QPs.  Bimodal distribution
                  // = the per-QP rotation-asymmetry bug we suspect.  Build with
                  //   make CXXFLAGS_USER='-DSEAS_DIAG_TPV104_FAULT_BASIS' ...
                  static std::set<int> seen_dofs;
                  if (seen_dofs.insert(dof_idx).second)
                  {
                     std::fprintf(stderr,
                        "[diag-flip] dof=%d sf=%d "
                        "raw_n=(%+.6f,%+.6f,%+.6f) "
                        "raw_t1=(%+.6f,%+.6f,%+.6f) "
                        "raw_t2=(%+.6f,%+.6f,%+.6f) "
                        "can_n=(%+.6f,%+.6f,%+.6f) "
                        "elem1_on_plus=%d\n",
                        dof_idx, qpd.sign_flipped ? 1 : 0,
                        qpd.normal[0], qpd.normal[1], qpd.normal[2],
                        qpd.tangent1[0], qpd.tangent1[1], qpd.tangent1[2],
                        qpd.tangent2[0], qpd.tangent2[1], qpd.tangent2[2],
                        can_n[0], can_n[1], can_n[2],
                        elem1_on_plus ? 1 : 0);
                  }
#endif

                  real_t I_self_can[NUM_STATE], I_nbr_can[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     I_self_can[c] = 0.0; I_nbr_can[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        I_self_can[c] += Tinv_can(c, k) * I_self[k];
                        I_nbr_can[c]  += Tinv_can(c, k) * I_nbr[k];
                     }
                  }

                  const real_t *I_plus_local  = elem1_on_plus ? I_self_can
                                                              : I_nbr_can;
                  const real_t *I_minus_local = elem1_on_plus ? I_nbr_can
                                                              : I_self_can;

                  real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
                  // R-602/R-603 substep dispatch:  if the driver has
                  // pre-computed per-substep imposed states via
                  // Tpv104SubStepIterator::AdvanceWithSubStepStates and
                  // installed them via SetSubStepFaultImposedStates, consume
                  // them here in lieu of running EvaluateADER inline.  Layout:
                  //   substep_I_imp_*_flat_[dof_idx * NUM_STATE + c]
                  // Both arrays are in the canonical fault-local frame —
                  // SAME frame as EvaluateADER's outputs — so the downstream
                  // T_can rotation back to global is unchanged.  When the
                  // pointers are null (default), the inline EvaluateADER
                  // path runs unchanged (bit-identical to pre-change).
                  if (substep_I_imp_plus_flat_ != nullptr &&
                      substep_I_imp_minus_flat_ != nullptr &&
                      dof_idx >= 0 &&
                      dof_idx < substep_n_total_fault_qps_)
                  {
                     const real_t *src_p =
                        substep_I_imp_plus_flat_ + dof_idx * NUM_STATE;
                     const real_t *src_m =
                        substep_I_imp_minus_flat_ + dof_idx * NUM_STATE;
                     for (int c = 0; c < NUM_STATE; c++)
                     {
                        I_imp_plus[c]  = src_p[c];
                        I_imp_minus[c] = src_m[c];
                     }
                  }
                  else
                  {
                     // v9.4.0 Commit 3: fluctuation-Q ADER dispatch;
                     // has_bulk_bg_ already asserted at the top of
                     // ComputeADERFaceFluxRHS (Q_bg = 0 is valid).
                     // REVIEW R-016 + Phase H.6 (rev-3): route LSW callers
                     // (TPV205) through EvaluateADER_LSW; LSW_ForcedRupture
                     // callers (SAFS dynamic-rupture driver, D-4 nucleation)
                     // through the time-dependent variant; the default
                     // RateAndState path is byte-identical to pre-change.
                     if (fault_friction_law_ ==
                         FaultFrictionLaw::LSW_ForcedRupture)
                     {
                        // R-401 / R-502 guard hoisted out of the per-DOF
                        // loop: time_was_set_ is per-WaveOperator-call.
                        WaveOperator<MeshType>::VerifyForcedRuptureTimeReady(
                           time_was_set_, fdata.T_forced_rupture);
                        fault_flux_->EvaluateADER_LSW_ForcedRupture(
                           fdata,
                           I_plus_local, I_minus_local,
                           dt,
                           GetTime(),
                           I_imp_plus, I_imp_minus);
                     }
                     else if (fault_friction_law_ == FaultFrictionLaw::LSW)
                     {
                        fault_flux_->EvaluateADER_LSW(
                           fdata,
                           I_plus_local, I_minus_local,
                           dt,
                           I_imp_plus, I_imp_minus);
                     }
                     else
                     {
                        fault_flux_->EvaluateADER(fdata,
                                                  I_plus_local, I_minus_local,
                                                  dt,
                                                  I_imp_plus, I_imp_minus);
                     }
                  }

                  real_t I_imp_plus_g[NUM_STATE], I_imp_minus_g[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     I_imp_plus_g[c] = 0.0; I_imp_minus_g[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        I_imp_plus_g[c]  += T_can(c, k) * I_imp_plus[k];
                        I_imp_minus_g[c] += T_can(c, k) * I_imp_minus[k];
                     }
                  }

#ifdef SEAS_DIAG_FAULT_FLUX
                  // C-1s STRESS-ROT (ADER interior-fault path): production
                  // dispatch goes through EvaluateADER here, NOT the
                  // per-stage Evaluate path probed by C-1s INT/SHR above
                  // (Frontera 7677529 confirmed 0 INT/SHR lines, 7024
                  // C-1n NORMAL → production path is this ADER one).
                  // Variables are time-integrated I_* — print Q_avg = I/dt
                  // to match C-1n NORMAL probe semantics (which sees Q_avg
                  // inside Evaluate via EvaluateADER's I/dt division).
                  if (dof_idx >= 0 &&
                      dof_idx < static_cast<int>(fault_dof_data_->size()) &&
                      (*fault_dof_data_)[dof_idx].diag_print &&
                      dt > 0)
                  {
                     const real_t inv_dt = 1.0 / dt;

                     std::fprintf(stderr,
                        "[C-1s ADER-INT BASIS] rank=%d  "
                        "can_n=(%+.4e,%+.4e,%+.4e)  "
                        "can_t1=(%+.4e,%+.4e,%+.4e)  "
                        "can_t2=(%+.4e,%+.4e,%+.4e)  "
                        "sign_flipped=%d  elem1_on_plus=%d  dt=%.4e\n",
                        g_seas_my_rank,
                        can_n[0], can_n[1], can_n[2],
                        can_t1[0], can_t1[1], can_t1[2],
                        can_t2[0], can_t2[1], can_t2[2],
                        qpd.sign_flipped ? 1 : 0,
                        elem1_on_plus  ? 1 : 0,
                        dt);

                     auto _c1s_print = [&](const char *tag, const real_t *I)
                     {
                        // Q_avg = I / dt — matches Evaluate's view.
                        std::fprintf(stderr,
                           "[C-1s ADER-INT %s] rank=%d  "
                           "SXX=%+.4e SYY=%+.4e SZZ=%+.4e  "
                           "SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
                           "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
                           tag, g_seas_my_rank,
                           I[SXX]*inv_dt, I[SYY]*inv_dt, I[SZZ]*inv_dt,
                           I[SXY]*inv_dt, I[SYZ]*inv_dt, I[SXZ]*inv_dt,
                           I[VX]*inv_dt,  I[VY]*inv_dt,  I[VZ]*inv_dt);
                     };

                     // I_self/I_nbr are GLOBAL frame (read at line ~2167).
                     // Map to canonical ± via elem1_on_plus.
                     const real_t *I_glob_p = elem1_on_plus ? I_self : I_nbr;
                     const real_t *I_glob_m = elem1_on_plus ? I_nbr  : I_self;

                     _c1s_print("GLOB-IN+",  I_glob_p);
                     _c1s_print("GLOB-IN-",  I_glob_m);
                     _c1s_print("LOC-IN+",   I_plus_local);
                     _c1s_print("LOC-IN-",   I_minus_local);
                     _c1s_print("LOC-IMP+",  I_imp_plus);
                     _c1s_print("LOC-IMP-",  I_imp_minus);
                     _c1s_print("GLOB-OUT+", I_imp_plus_g);
                     _c1s_print("GLOB-OUT-", I_imp_minus_g);
                  }
#endif

                  real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
                  flux_.Interior(can_n, I_imp_plus_g,  I_imp_plus_g,
                                 F_h_plus);
                  flux_.Interior(can_n, I_imp_minus_g, I_imp_minus_g,
                                 F_h_minus);

                  if (elem1_on_plus)
                  {
                     for (int c = 0; c < NUM_STATE; c++)
                     {
                        for (int i = 0; i < ndof; i++)
                        {
                           rhs[c * ndof_total_ + dof_offset1 + i] -=
                              w * shape1(i) * F_h_plus[c];
                           rhs[c * ndof_total_ + dof_offset2 + i] +=
                              w * shape2(i) * F_h_minus[c];
                        }
                     }
                  }
                  else
                  {
                     for (int c = 0; c < NUM_STATE; c++)
                     {
                        for (int i = 0; i < ndof; i++)
                        {
                           rhs[c * ndof_total_ + dof_offset1 + i] +=
                              w * shape1(i) * F_h_minus[c];
                           rhs[c * ndof_total_ + dof_offset2 + i] -=
                              w * shape2(i) * F_h_plus[c];
                        }
                     }
                  }
               }
               else
               {
                  MFEM_ABORT("ComputeADERFaceFluxRHS: interior fault face f="
                             << f << " has no FaultBasis/DOFData mapping.");
               }
            }  // end if (mode == None)
            else
            {
               // Round-12 Patch 2: face-averaging path.  Process at
               // q == 0; for q > 0 on this fault face, skip the rest
               // of the outer-q body (the face was fully batched at
               // q = 0 already).
               if (q != 0) { continue; }

               auto it = fault_face_dof_offset_.find(f);
               MFEM_VERIFY(it != fault_face_dof_offset_.end(),
                           "ComputeADERFaceFluxRHS (avg): fault face "
                           << f << " has no DOFData offset.");
               const int face_dof_base = it->second;

               const int fb_idx = LookupInteriorFaultBasisIndex(f);
               MFEM_VERIFY(fault_basis_ && fb_idx >= 0 &&
                           fb_idx < fault_basis_->NumFaces(),
                           "ComputeADERFaceFluxRHS (avg): fault face "
                           << f << " has no FaultBasis mapping.");
               const FaultBasisData &bd = fault_basis_->GetBasis(fb_idx);
               MFEM_VERIFY(nqp == static_cast<int>(bd.qp_data.size()),
                           "ComputeADERFaceFluxRHS (avg): face nqp="
                           << nqp << " != basis qp_data size="
                           << bd.qp_data.size());

               const FiniteElement *fe2_avg = fes_->GetFE(e2);
               const int face_nqp = nqp;
               const real_t inv_dt = 1.0 / dt;

               std::vector<EvalStageState> states(face_nqp);
               std::vector<real_t>   w_per_qp(face_nqp);
               std::vector<Vector>   shape1_per_qp(face_nqp);
               std::vector<Vector>   shape2_per_qp(face_nqp);
               std::vector<std::array<real_t, 3>> can_n_per_qp(face_nqp);
               std::vector<DenseMatrix> T_can_per_qp(face_nqp);
               std::vector<bool>     elem1_on_plus_per_qp(face_nqp);
               std::vector<int>      dof_idx_per_qp(face_nqp);
               std::vector<std::array<real_t, NUM_STATE>> I_plus_loc_per_qp(face_nqp);
               std::vector<std::array<real_t, NUM_STATE>> I_minus_loc_per_qp(face_nqp);

               // Pass 1: reconstruct everything per QP and compute the
               // full stage chain (ComputeStageState).
               for (int qq = 0; qq < face_nqp; qq++)
               {
                  const IntegrationPoint &ip_qq = ir.IntPoint(qq);
                  ftr->SetAllIntPoints(&ip_qq);

                  Vector nor_vec_qq(3);
                  CalcOrtho(ftr->Face->Jacobian(), nor_vec_qq);
                  const real_t nor_len_qq = nor_vec_qq.Norml2();
                  if (nor_len_qq > 0) { nor_vec_qq /= nor_len_qq; }
                  w_per_qp[qq] = ip_qq.weight * nor_len_qq;

                  IntegrationPoint ip1_qq;
                  ftr->Loc1.Transform(ip_qq, ip1_qq);
                  shape1_per_qp[qq].SetSize(ndof);
                  fe1->CalcShape(ip1_qq, shape1_per_qp[qq]);

                  IntegrationPoint ip2_qq;
                  ftr->Loc2.Transform(ip_qq, ip2_qq);
                  shape2_per_qp[qq].SetSize(ndof);
                  fe2_avg->CalcShape(ip2_qq, shape2_per_qp[qq]);

                  real_t I_self_qq[NUM_STATE], I_nbr_qq[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     I_self_qq[c] = 0.0; I_nbr_qq[c] = 0.0;
                     for (int i = 0; i < ndof; i++)
                     {
                        I_self_qq[c] += shape1_per_qp[qq](i)
                           * I_data[c * ndof_total_ + dof_offset1 + i];
                        I_nbr_qq[c]  += shape2_per_qp[qq](i)
                           * I_data[c * ndof_total_ + dof_offset2 + i];
                     }
                  }

                  dof_idx_per_qp[qq] = face_dof_base + qq;
                  MFEM_VERIFY(dof_idx_per_qp[qq] >= 0 &&
                              dof_idx_per_qp[qq] <
                                 static_cast<int>(fault_dof_data_->size()),
                              "ComputeADERFaceFluxRHS (avg): dof_idx out of range");

                  // R-101: per-face geometric flag, FP-stable across QPs.
                  MFEM_ASSERT(fb_idx >= 0 &&
                              fb_idx < static_cast<int>(
                                 interior_fault_elem1_on_plus_.size()),
                              "R-101: interior_fault_elem1_on_plus_ "
                              "missing entry for fb_idx=" << fb_idx);
                  elem1_on_plus_per_qp[qq] =
                     interior_fault_elem1_on_plus_[fb_idx];

                  const FaultBasisQPData &qpd = bd.qp_data[qq];
                  const bool should_negate_frame_qq =
                     !elem1_on_plus_per_qp[qq];
                  real_t can_t1[3], can_t2[3];
                  for (int d = 0; d < 3; d++)
                  {
                     can_n_per_qp[qq][d] = should_negate_frame_qq
                                              ? -qpd.normal[d]
                                              :  qpd.normal[d];
                     can_t1[d] = should_negate_frame_qq ? -qpd.tangent1[d]
                                                        :  qpd.tangent1[d];
                     can_t2[d] = should_negate_frame_qq ? -qpd.tangent2[d]
                                                        :  qpd.tangent2[d];
                  }

                  T_can_per_qp[qq].SetSize(NUM_STATE);
                  DenseMatrix Tinv_can_qq(NUM_STATE);
                  GodunovFlux::BuildRotation(can_n_per_qp[qq].data(),
                                              can_t1, can_t2,
                                              T_can_per_qp[qq]);
                  GodunovFlux::BuildRotationInverse(can_n_per_qp[qq].data(),
                                                     can_t1, can_t2,
                                                     Tinv_can_qq);

                  real_t I_self_can_qq[NUM_STATE], I_nbr_can_qq[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     I_self_can_qq[c] = 0.0; I_nbr_can_qq[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        I_self_can_qq[c] += Tinv_can_qq(c, k) * I_self_qq[k];
                        I_nbr_can_qq[c]  += Tinv_can_qq(c, k) * I_nbr_qq[k];
                     }
                  }

                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     if (elem1_on_plus_per_qp[qq])
                     {
                        I_plus_loc_per_qp[qq][c]  = I_self_can_qq[c];
                        I_minus_loc_per_qp[qq][c] = I_nbr_can_qq[c];
                     }
                     else
                     {
                        I_plus_loc_per_qp[qq][c]  = I_nbr_can_qq[c];
                        I_minus_loc_per_qp[qq][c] = I_self_can_qq[c];
                     }
                  }

                  real_t Q_avg_plus_qq[NUM_STATE];
                  real_t Q_avg_minus_qq[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     Q_avg_plus_qq[c]  = I_plus_loc_per_qp[qq][c]  * inv_dt;
                     Q_avg_minus_qq[c] = I_minus_loc_per_qp[qq][c] * inv_dt;
                  }

                  const DOFData &fdata_const =
                     (*fault_dof_data_)[dof_idx_per_qp[qq]];
                  fault_flux_->ComputeStageState(fdata_const,
                                                  Q_avg_plus_qq,
                                                  Q_avg_minus_qq,
                                                  states[qq]);
               }

               // Pass 2: average the selected stage field and run the
               // appropriate completion helper.
               const real_t inv_nqp = 1.0 / static_cast<real_t>(face_nqp);
               if (eval_avg_mode == FaultEvalStageAvgMode::Trial)
               {
                  real_t sn_avg = 0.0, t1_avg = 0.0, t2_avg = 0.0;
                  for (int qq = 0; qq < face_nqp; qq++)
                  {
                     sn_avg += states[qq].sigma_n_trial;
                     t1_avg += states[qq].tau1_trial;
                     t2_avg += states[qq].tau2_trial;
                  }
                  sn_avg *= inv_nqp; t1_avg *= inv_nqp; t2_avg *= inv_nqp;
                  for (int qq = 0; qq < face_nqp; qq++)
                  {
                     states[qq].sigma_n_trial = sn_avg;
                     states[qq].tau1_trial    = t1_avg;
                     states[qq].tau2_trial    = t2_avg;
                     const DOFData &fdata_c =
                        (*fault_dof_data_)[dof_idx_per_qp[qq]];
                     fault_flux_->CompleteFromTrial(fdata_c, states[qq]);
                  }
               }
               else if (eval_avg_mode == FaultEvalStageAvgMode::Theta)
               {
                  real_t Theta_avg = 0.0;
                  for (int qq = 0; qq < face_nqp; qq++)
                  {
                     Theta_avg += states[qq].Theta;
                  }
                  Theta_avg *= inv_nqp;
                  for (int qq = 0; qq < face_nqp; qq++)
                  {
                     states[qq].Theta = Theta_avg;
                     const DOFData &fdata_c =
                        (*fault_dof_data_)[dof_idx_per_qp[qq]];
                     fault_flux_->CompleteFromTheta(fdata_c, states[qq]);
                  }
               }
               else if (eval_avg_mode == FaultEvalStageAvgMode::Vabs)
               {
                  real_t Vabs_avg = 0.0;
                  for (int qq = 0; qq < face_nqp; qq++)
                  {
                     Vabs_avg += states[qq].V_abs;
                  }
                  Vabs_avg *= inv_nqp;
                  for (int qq = 0; qq < face_nqp; qq++)
                  {
                     states[qq].V_abs = Vabs_avg;
                     const DOFData &fdata_c =
                        (*fault_dof_data_)[dof_idx_per_qp[qq]];
                     fault_flux_->CompleteFromVabs(fdata_c, states[qq]);
                  }
               }
               else if (eval_avg_mode == FaultEvalStageAvgMode::Tcorr)
               {
                  real_t sn_avg = 0.0, t1_avg = 0.0, t2_avg = 0.0;
                  for (int qq = 0; qq < face_nqp; qq++)
                  {
                     sn_avg += states[qq].sigma_n_corr;
                     t1_avg += states[qq].tau1_corr;
                     t2_avg += states[qq].tau2_corr;
                  }
                  sn_avg *= inv_nqp; t1_avg *= inv_nqp; t2_avg *= inv_nqp;
                  for (int qq = 0; qq < face_nqp; qq++)
                  {
                     states[qq].sigma_n_corr = sn_avg;
                     states[qq].tau1_corr    = t1_avg;
                     states[qq].tau2_corr    = t2_avg;
                     // Tcorr: no downstream recompute — per plan spec.
                  }
               }

               // Pass 3: build imposed states, write back DOFData,
               // compute per-side flux, deposit into rhs.
               for (int qq = 0; qq < face_nqp; qq++)
               {
                  DOFData &fdata_qq =
                     (*fault_dof_data_)[dof_idx_per_qp[qq]];

                  real_t Q_avg_plus_qq[NUM_STATE];
                  real_t Q_avg_minus_qq[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     Q_avg_plus_qq[c]  = I_plus_loc_per_qp[qq][c]  * inv_dt;
                     Q_avg_minus_qq[c] = I_minus_loc_per_qp[qq][c] * inv_dt;
                  }

                  real_t Q_imp_plus_qq[NUM_STATE];
                  real_t Q_imp_minus_qq[NUM_STATE];
                  fault_flux_->BuildImposedState(fdata_qq, states[qq],
                                                  Q_avg_plus_qq,
                                                  Q_avg_minus_qq,
                                                  Q_imp_plus_qq,
                                                  Q_imp_minus_qq);

                  real_t I_imp_plus_qq[NUM_STATE];
                  real_t I_imp_minus_qq[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     I_imp_plus_qq[c]  = Q_imp_plus_qq[c]  * dt;
                     I_imp_minus_qq[c] = Q_imp_minus_qq[c] * dt;
                  }

                  fault_flux_->WriteBackState(fdata_qq, states[qq]);

                  real_t I_imp_plus_g_qq[NUM_STATE];
                  real_t I_imp_minus_g_qq[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     I_imp_plus_g_qq[c]  = 0.0;
                     I_imp_minus_g_qq[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        I_imp_plus_g_qq[c]  +=
                           T_can_per_qp[qq](c, k) * I_imp_plus_qq[k];
                        I_imp_minus_g_qq[c] +=
                           T_can_per_qp[qq](c, k) * I_imp_minus_qq[k];
                     }
                  }

                  real_t F_h_plus_qq[NUM_STATE];
                  real_t F_h_minus_qq[NUM_STATE];
                  flux_.Interior(can_n_per_qp[qq].data(),
                                  I_imp_plus_g_qq, I_imp_plus_g_qq,
                                  F_h_plus_qq);
                  flux_.Interior(can_n_per_qp[qq].data(),
                                  I_imp_minus_g_qq, I_imp_minus_g_qq,
                                  F_h_minus_qq);

                  const Vector &sh1_qq = shape1_per_qp[qq];
                  const Vector &sh2_qq = shape2_per_qp[qq];
                  const real_t   w_qq  = w_per_qp[qq];

                  if (elem1_on_plus_per_qp[qq])
                  {
                     for (int c = 0; c < NUM_STATE; c++)
                     {
                        for (int i = 0; i < ndof; i++)
                        {
                           rhs[c * ndof_total_ + dof_offset1 + i] -=
                              w_qq * sh1_qq(i) * F_h_plus_qq[c];
                           rhs[c * ndof_total_ + dof_offset2 + i] +=
                              w_qq * sh2_qq(i) * F_h_minus_qq[c];
                        }
                     }
                  }
                  else
                  {
                     for (int c = 0; c < NUM_STATE; c++)
                     {
                        for (int i = 0; i < ndof; i++)
                        {
                           rhs[c * ndof_total_ + dof_offset1 + i] +=
                              w_qq * sh1_qq(i) * F_h_minus_qq[c];
                           rhs[c * ndof_total_ + dof_offset2 + i] -=
                              w_qq * sh2_qq(i) * F_h_plus_qq[c];
                        }
                     }
                  }
               }
            }  // end avg-mode path
         }
         else if (use_precomputed_face_fluxes_ && is_boundary)
         {
            // Phase 2a (§6.2): boundary dispatch through precomputed
            // table.  bulk_bg_scaled must be bit-zero (Phase 1 probe
            // construction required Q_bg = 0 for linearity).
            // NOTE: Arm 3d Full-path (AddBoundaryFaceRhsFull) is
            // available but NOT wired into this dispatch yet — it showed
            // 3.9e-01 rel drift vs runtime on plane-wave inputs in
            // P_SWITCH_ON_LINEAR_EQ despite reducing orbit drift by 9+
            // orders of magnitude on constant-Q.  Root cause under
            // investigation; keep legacy path until Full matches runtime
            // bit-for-bit on non-orbit-asymmetric inputs.
            precomputed_face_fluxes_.AddBoundaryFaceRhs(
               f, I_self, bulk_bg_scaled, w, shape1.GetData(),
               ndof, dof_offset1, ndof_total_, rhs);
         }
         else if (use_precomputed_face_fluxes_)   // interior non-fault
         {
            // Phase 2a (§6.2 R2-003): legacy AddInteriorFaceRhs path.
            // Arm 3d Full-path available via direct call
            // (AddInteriorFaceRhsFull) — used by the Phase 3 audit with
            // SEAS_TEST_USE_ARM3D=1 as an opt-in diagnostic, but not
            // wired here pending plane-wave equivalence investigation.
            precomputed_face_fluxes_.AddInteriorFaceRhs(
               f, /*caller_elem=*/e1, I_self, I_nbr, w, shape1.GetData(),
               ndof, dof_offset1, ndof_total_, rhs);
            precomputed_face_fluxes_.AddInteriorFaceRhs(
               f, /*caller_elem=*/e2, I_nbr, I_self, w, shape2.GetData(),
               ndof, dof_offset2, ndof_total_, rhs);
         }
         else if (is_boundary)
         {
            // Runtime boundary dispatch — lifted VERBATIM from the
            // pre-patch `if (is_boundary) { switch(bc_type) {...} ... }`
            // block of ComputeADERFaceFluxRHS (R2-008 discipline:
            // preserved unchanged when the flag is off).  Round-13C
            // Patch 1 wraps the dispatch with an optional n↔-n
            // symmetrization (SEAS_TEST_NONFAULT_BOTH_SYM=1); unset
            // env → byte-identical.
            FaceBC bc_type = ClassifyBoundaryFace(bdr_attr);

            auto compute_bc_flux = [&](const real_t *nvec, real_t *F_out)
            {
               switch (bc_type)
               {
                  case FaceBC::Absorbing:
                     flux_.AbsorbingTotal(nvec, I_self, bulk_bg_scaled, F_out);
                     break;
                  case FaceBC::FreeSurface:
                     if (free_surface_bc_mode_ == FreeSurfaceBCMode::Godunov)
                     {
                        flux_.FreeSurfaceGodunovTotal(nvec, I_self,
                                                      bulk_bg_scaled, F_out);
                     }
                     else
                     {
                        flux_.FreeSurfaceTotal(nvec, I_self,
                                               bulk_bg_scaled, F_out);
                     }
                     break;
                  case FaceBC::Fault:
                     MFEM_ABORT("ComputeADERFaceFluxRHS: fault face " << f
                                << " is 1-sided (no neighbor).  Fault faces "
                                "must be 2-sided interior faces.");
                     break;
                  default:
                     flux_.AbsorbingTotal(nvec, I_self, bulk_bg_scaled, F_out);
                     break;
               }
            };

            // R-1411: nonfault_both_sym always false post-R-1203 abort;
            // the dead else-branch (0.5*(F(n)+F(-n))) was identically
            // zero for any conservative flux and has been removed.
            compute_bc_flux(nor, F_h);

            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  rhs[c * ndof_total_ + dof_offset1 + i] -=
                     w * shape1(i) * F_h[c];
               }
            }
         }
         else
         {
            // Runtime interior dispatch — lifted VERBATIM from the
            // pre-patch `else { flux_.Interior(nor, I_self, I_nbr, F_h); ... }`
            // block (non-fault interior fallback; R2-008 discipline).
            // Round-13C Patch 1 wraps with optional n↔-n,L↔R symmetrization
            // (SEAS_TEST_NONFAULT_BOTH_SYM=1); unset → byte-identical.
            //
            // Round-11 Mixed-Flux dispatch (R-1206): pick Central or
            // Interior based on the central-flux face set.  Both
            // symmetrization branches honor the mode.
            const bool use_central_here =
               (mf_on && central_flux_face_set_.count(f) > 0);
            auto interior_or_central = [&](const real_t *n_arg,
                                           const real_t *Q_a,
                                           const real_t *Q_b,
                                           real_t *F_out)
            {
               if (use_central_here) { flux_.Central(n_arg, Q_a, Q_b, F_out); }
               else                  { flux_.Interior(n_arg, Q_a, Q_b, F_out); }
            };

            // R-1411: nonfault_both_sym dead else-branch removed.
            interior_or_central(nor, I_self, I_nbr, F_h);

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

#ifdef SEAS_DIAG_FAULT_FLUX
            // C-2B NONFAULT-FACE (ADER path): per-AdvanceADER-call dump
            // on the 6 non-fault interior faces of the hypocenter's two
            // adjacent tets.  Captures I_self, I_nbr (time-integrated
            // states on both sides) and F_h (Godunov flux).  Fires
            // FIRST QP per face per AdvanceADER call (q==0,
            // user-decision: per-Mult, not per-stage / per-QP).
            if (q == 0 && !diag_nonfault_faces_.empty())
            {
               bool match = false;
               for (size_t kk = 0;
                    kk < diag_nonfault_faces_.size(); kk++)
               {
                  if (diag_nonfault_faces_[kk] == f) { match = true; break; }
               }
               if (match)
               {
                  const char *side_p = (e1 == diag_elem_plus_  ||
                                        e2 == diag_elem_plus_)  ? "P" : "";
                  const char *side_m = (e1 == diag_elem_minus_ ||
                                        e2 == diag_elem_minus_) ? "M" : "";
                  std::fprintf(stderr,
                     "[C-2B NONFAULT-FACE] rank=%d face=%d e1=%d e2=%d "
                     "side=%s%s nor=(%+.4e,%+.4e,%+.4e) w=%.4e\n",
                     g_seas_my_rank, f, e1, e2, side_p, side_m,
                     nor[0], nor[1], nor[2], w);
                  std::fprintf(stderr,
                     "[C-2B NONFAULT-FACE I_self] rank=%d face=%d "
                     "SXX=%+.4e SYY=%+.4e SZZ=%+.4e "
                     "SXY=%+.4e SYZ=%+.4e SXZ=%+.4e "
                     "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
                     g_seas_my_rank, f,
                     I_self[SXX], I_self[SYY], I_self[SZZ],
                     I_self[SXY], I_self[SYZ], I_self[SXZ],
                     I_self[VX],  I_self[VY],  I_self[VZ]);
                  std::fprintf(stderr,
                     "[C-2B NONFAULT-FACE I_nbr]  rank=%d face=%d "
                     "SXX=%+.4e SYY=%+.4e SZZ=%+.4e "
                     "SXY=%+.4e SYZ=%+.4e SXZ=%+.4e "
                     "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
                     g_seas_my_rank, f,
                     I_nbr[SXX], I_nbr[SYY], I_nbr[SZZ],
                     I_nbr[SXY], I_nbr[SYZ], I_nbr[SXZ],
                     I_nbr[VX],  I_nbr[VY],  I_nbr[VZ]);
                  std::fprintf(stderr,
                     "[C-2B NONFAULT-FACE F_h]    rank=%d face=%d "
                     "SXX=%+.4e SYY=%+.4e SZZ=%+.4e "
                     "SXY=%+.4e SYZ=%+.4e SXZ=%+.4e "
                     "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
                     g_seas_my_rank, f,
                     F_h[SXX], F_h[SYY], F_h[SZZ],
                     F_h[SXY], F_h[SYZ], F_h[SXZ],
                     F_h[VX],  F_h[VY],  F_h[VZ]);
               }
            }
#endif
         }
      }
   }
}

// ---------------------------------------------------------------------------
// ADER I-05 Phase 6: ComputeADERSharedFaceFluxRHS — partition-seam faces.
// ---------------------------------------------------------------------------
// Parallel counterpart of ComputeADERFaceFluxRHS.  Serial builds return
// immediately.  Structure mirrors ComputeSharedFaceFluxRHS with the same
// I/dt/EvaluateADER[Total] substitutions as the local-face routine.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ComputeADERSharedFaceFluxRHS(const Vector &I,
                                                          real_t dt,
                                                          Vector &rhs) const
{
   if constexpr (!IsParallelMesh<MeshType>::value)
   {
      return;
   }
   else
   {
#ifdef MFEM_USE_MPI
      MFEM_VERIFY(dt > 0.0,
                  "ComputeADERSharedFaceFluxRHS: dt must be > 0");
      MFEM_VERIFY(I.Size() == NUM_STATE * ndof_total_,
                  "ComputeADERSharedFaceFluxRHS: I size mismatch");

      // Round-11 Mixed-Flux short-circuit (R-1206): hoist once per call.
      const bool mf_on = mf_on_;  // R-1208: cached member, kept in sync by SetMixedFluxMode

      auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
      MFEM_VERIFY(pfes, "FESpace must be ParFiniteElementSpace for ParMesh");

      auto &pmesh = static_cast<const ParMesh &>(mesh_);
      int n_shared = pmesh.GetNSharedFaces();
      if (n_shared == 0) { return; }

      // v9.4.0: see ComputeADERFaceFluxRHS top-level note.  REVIEW R-007:
      // updated stale "Total-Q only" wording.
      //
      // R-1505: COLLECTIVE consensus check (shared corrector mirror).  See
      // the matching block in ComputeADERFaceFluxRHS for the rationale.
      // Only the shared corrector's ExchangeFaceNbrData below would
      // actually deadlock on a rank-local abort, but checking here too
      // costs ~one Allreduce per AdvanceADER call (negligible) and keeps
      // the two corrector entry points self-documenting about their
      // collective contract.
      {
#ifdef MFEM_USE_MPI
         auto &pmesh_consensus = static_cast<const ParMesh &>(mesh_);
         const int my_has = has_bulk_bg_ ? 1 : 0;
         int min_has = 0;
         MPI_Allreduce(&my_has, &min_has, 1, MPI_INT, MPI_MIN,
                       pmesh_consensus.GetComm());
         MFEM_VERIFY(min_has == 1,
                     "ComputeADERSharedFaceFluxRHS: SetAbsorbingBackground"
                     "(Q_bg) was NOT called on every rank (min has_bulk_bg_="
                     << min_has << ").  A rank-local abort here would "
                     "deadlock the next q_gf.ExchangeFaceNbrData; the "
                     "collective check fails loud everywhere instead.  "
                     "See R-1505.");
#else
         MFEM_VERIFY(has_bulk_bg_,
                     "wave.AdvanceADER() requires SetAbsorbingBackground("
                     "Q_bg) to have been called (Q_bg = 0 is valid under "
                     "fluctuation-Q dispatch).");
#endif
      }

      if (bc_.fault_attr > 0)
      {
         MFEM_VERIFY(fault_flux_ && fault_dof_data_,
                     "WaveOperator::ComputeADERSharedFaceFluxRHS: "
                     "bc_.fault_attr=" << bc_.fault_attr
                     << " > 0 but fault bookkeeping unset.");
      }

      // R-1203: SEAS_TEST_NONFAULT_BOTH_SYM is permanently disabled.
      // The local-face hook in ComputeADERFaceFluxRHS aborts on the env
      // var; the shared-face mirror does the same here for parity.
      // See R-1203 description in ComputeADERFaceFluxRHS for the
      // anti-symmetry argument.
      const bool nonfault_both_sym = false;
      if (std::getenv("SEAS_TEST_NONFAULT_BOTH_SYM"))
      {
         MFEM_ABORT("SEAS_TEST_NONFAULT_BOTH_SYM is permanently disabled "
                    "(R-1203, shared-face mirror).  See the local-face "
                    "abort message in ComputeADERFaceFluxRHS.");
      }

      const real_t *I_data = I.GetData();
      MFEM_VERIFY(ghost_gf_, "Ghost GF not initialized");

      ParGridFunction &q_gf = *ghost_gf_;
      std::vector<Vector> nbr_data(NUM_STATE);

      for (int c = 0; c < NUM_STATE; c++)
      {
         for (int i = 0; i < ndof_total_; i++)
         {
            q_gf[i] = I_data[c * ndof_total_ + i];
         }
         q_gf.ExchangeFaceNbrData();
         const Vector &src = q_gf.FaceNbrData();
         nbr_data[c].SetSize(src.Size());
         std::memcpy(nbr_data[c].GetData(), src.GetData(),
                     src.Size() * sizeof(real_t));
      }

      const int nfs = fault_shared_faces_.Size();
      const bool fault_active = (fault_flux_ && fault_dof_data_
                                 && nfs > 0 && nbf_per_face_ > 0
                                 && fault_basis_);

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

         // Phase 2b (§7.2): mesh face index for the precomputed-table
         // lookup.  `sf` is the shared-face index; the key for
         // AddInteriorFaceRhs is (mesh_face_idx, Elem1No).
         const int mesh_face_idx =
            const_cast<ParMesh &>(pmesh).GetSharedFace(sf);

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

            IntegrationPoint ip1;
            ftr->Loc1.Transform(ip, ip1);
            Vector shape1(ndof);
            fe1->CalcShape(ip1, shape1);

            real_t I_self[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               I_self[c] = 0.0;
               for (int i = 0; i < ndof; i++)
               {
                  I_self[c] += shape1(i) * I_data[c * ndof_total_ + dof_offset1 + i];
               }
            }

            IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2);
            Vector shape2(ndof2);
            fe2->CalcShape(ip2, shape2);

            // R-1508: promoted from MFEM_ASSERT (Debug-only) to MFEM_VERIFY
            // so a Release build still fails loud on a heterogeneous ghost
            // layer rather than silently producing wrong shape evaluations.
            MFEM_VERIFY(ndof2 == ndof_per_el_,
                        "ComputeADERSharedFaceFluxRHS: heterogeneous ghost "
                        "element (ndof2=" << ndof2 << " vs ndof_per_el_="
                        << ndof_per_el_ << ") not supported.");

            real_t I_nbr[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               I_nbr[c] = 0.0;
               for (int i = 0; i < ndof2; i++)
               {
                  I_nbr[c] += shape2(i) * nbr_data[c][nbr_idx * ndof_per_el_ + i];
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
                  const int sf_idx_in_fault = basis_idx
                                              - fault_interior_faces_.Size();
                  MFEM_VERIFY(sf_idx_in_fault >= 0 &&
                              sf_idx_in_fault <
                              static_cast<int>(shared_fault_elem1_on_plus_.size()),
                              "shared_fault_elem1_on_plus_ missing entry");
                  const bool elem1_on_plus =
                     shared_fault_elem1_on_plus_[sf_idx_in_fault];

                  const FaultBasisData &bd = fault_basis_->GetBasis(basis_idx);
                  MFEM_VERIFY(q < static_cast<int>(bd.qp_data.size()),
                              "FaultBasis::qp_data not populated for shared "
                              "fault face");
                  const FaultBasisQPData &qpd = bd.qp_data[q];

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

                  real_t I_self_can[NUM_STATE], I_nbr_can[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     I_self_can[c] = 0.0; I_nbr_can[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        I_self_can[c] += Tinv_can(c, k) * I_self[k];
                        I_nbr_can[c]  += Tinv_can(c, k) * I_nbr[k];
                     }
                  }

                  const real_t *I_plus_local  = elem1_on_plus ? I_self_can
                                                              : I_nbr_can;
                  const real_t *I_minus_local = elem1_on_plus ? I_nbr_can
                                                              : I_self_can;

                  DOFData &fdata = (*fault_dof_data_)[dof_idx];
                  real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
                  // R-1601 (np>1 substep stability fix): the shared-fault
                  // substep dispatch was producing rank-dependent
                  // canonical-frame mismatches that overflowed Q within
                  // ~7 macro-steps at np=10 (tau=4e28 fed into Brent →
                  // SIGABRT; reproduced on symmirror_1000m.msh).  The
                  // root cause is that R-101 `elem1_on_plus` is rank-
                  // local (rank A's Elem1 = rank B's Elem2 → opposite
                  // canonical frames on the same physical face) while
                  // `qpd.sign_flipped` is rank-independent.  Until the
                  // iterator's per-shared-QP physics is reconciled with
                  // the corrector's frame convention, fall back to
                  // inline EvaluateADER on SHARED QPs only — INTERIOR
                  // QPs continue to use the substep buffer.  This keeps
                  // substep semantics where they are well-tested and
                  // restores np>1 stability at the cost of mixed-mode
                  // dispatch on shared faces.  See
                  // SUBSTEP_NP_GT_1_HANG_REVIEW.md (R-1600 + R-1601).
                  //
                  // SHARED FALLBACK: this branch always runs the inline
                  // ADER closure regardless of substep_I_imp_*_flat_
                  // (R-1600/R-1601 frame-mismatch on shared QPs).
                  // REVIEW R-016 + Phase H.6 (rev-3): dispatch on the
                  // friction-law tag.  RateAndState keeps the Brent path
                  // byte-identical (TPV102/TPV104); LSW (TPV205) runs the
                  // closed-form solver; LSW_ForcedRupture (SAFS spatial
                  // dyn-driver, D-4 nucleation) runs the time-dependent
                  // variant.  Without LSW dispatch, TPV205 shared-fault
                  // QPs at np > 1 silently consume LSW values via Brent
                  // and stall the rupture front.
                  if (fault_friction_law_ == FaultFrictionLaw::LSW_ForcedRupture)
                  {
                     WaveOperator<MeshType>::VerifyForcedRuptureTimeReady(
                        time_was_set_, fdata.T_forced_rupture);
                     fault_flux_->EvaluateADER_LSW_ForcedRupture(
                        fdata,
                        I_plus_local, I_minus_local,
                        dt,
                        GetTime(),
                        I_imp_plus, I_imp_minus);
                  }
                  else if (fault_friction_law_ == FaultFrictionLaw::LSW)
                  {
                     fault_flux_->EvaluateADER_LSW(
                        fdata,
                        I_plus_local, I_minus_local,
                        dt,
                        I_imp_plus, I_imp_minus);
                  }
                  else
                  {
                     fault_flux_->EvaluateADER(fdata,
                                               I_plus_local, I_minus_local,
                                               dt,
                                               I_imp_plus, I_imp_minus);
                  }
                  (void)substep_I_imp_plus_flat_;
                  (void)substep_I_imp_minus_flat_;
                  (void)substep_n_total_fault_qps_;

                  // R-DIP2: cross-rank divergence diagnostic (RUNTIME env-gated,
                  // zero overhead when off).  At the target shared fault QP
                  // (env SEAS_DIAG_XRANK_QP="x,y,z"; default = the Frontera Dc2
                  // diverging QP) dump, per rank per sub-step, the decomposition
                  // INPUTS (canonical dip tangent + rotated bulk shear I_*_can)
                  // and OUTPUTS (DOFData dip/strike traction + slip-rate).
                  // Diff the two ranks' lines to localise the seed:
                  //   can_t1 differs            => H-C (curvilinear frame not
                  //                                bit-identical across ranks)
                  //   frame same, I_*_can differ => bulk-Q FP seed (H-B
                  //                                substrate; self vs nbr tells
                  //                                local-bulk vs ghost)
                  //   tau1_0 / tau1_nuc differ   => H-A (prestress/nuc gap)
                  //   then V1 growth over steps  => the weak-channel loop gain.
                  {
                     static const bool diag_on = []{
                        const char *e = std::getenv("SEAS_DIAG_XRANK");
                        return e && e[0] && !(e[0]=='0' && e[1]=='\0');
                     }();
                     if (diag_on && dof_idx >= 0 &&
                         dof_idx < static_cast<int>(fault_dof_data_->size()))
                     {
                        static const std::array<double,3> tgt = []{
                           std::array<double,3> t = {6.0751786666666670e+05,
                                                     3.7063591666666665e+06,
                                                     -4.5430145000000002e+03};
                           if (const char *e = std::getenv("SEAS_DIAG_XRANK_QP"))
                           { std::sscanf(e, "%lf,%lf,%lf", &t[0], &t[1], &t[2]); }
                           return t;
                        }();
                        static const double rad = []{
                           const char *e = std::getenv("SEAS_DIAG_XRANK_R");
                           return e ? std::atof(e) : 250.0;
                        }();
                        static const double tmin = []{
                           const char *e = std::getenv("SEAS_DIAG_XRANK_TMIN");
                           return e ? std::atof(e) : 0.0;
                        }();
                        Vector _phys(3);
                        ftr->Face->Transform(ip, _phys);
                        const double _dx = _phys(0)-tgt[0];
                        const double _dy = _phys(1)-tgt[1];
                        const double _dz = _phys(2)-tgt[2];
                        if (_dx*_dx + _dy*_dy + _dz*_dz <= rad*rad &&
                            GetTime() >= tmin)
                        {
                           std::fprintf(stderr,
                              "[XRANK] rank=%d t=%.6e c=(%.1f,%.1f,%.1f) "
                              "can_t1=(%+.10e,%+.10e,%+.10e) "
                              "Iself_xy=%+.10e Iself_xz=%+.10e "
                              "Inbr_xy=%+.10e Inbr_xz=%+.10e "
                              "tau1_0=%+.10e tau1_nuc=%+.10e "
                              "tau1_corr=%+.10e tau2_corr=%+.10e "
                              "V1=%+.10e V2=%+.10e "
                              "sign_flipped=%d elem1_on_plus=%d\n",
                              my_rank_, GetTime(),
                              _phys(0), _phys(1), _phys(2),
                              can_t1[0], can_t1[1], can_t1[2],
                              I_self_can[SXY], I_self_can[SXZ],
                              I_nbr_can[SXY],  I_nbr_can[SXZ],
                              fdata.tau1_0, fdata.tau1_nuc,
                              fdata.tau1_corr, fdata.tau2_corr,
                              fdata.V1, fdata.V2,
                              qpd.sign_flipped ? 1 : 0,
                              elem1_on_plus ? 1 : 0);
                           std::fflush(stderr);
                        }
                     }
                  }

                  real_t I_imp_plus_g[NUM_STATE], I_imp_minus_g[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     I_imp_plus_g[c] = 0.0; I_imp_minus_g[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        I_imp_plus_g[c]  += T_can(c, k) * I_imp_plus[k];
                        I_imp_minus_g[c] += T_can(c, k) * I_imp_minus[k];
                     }
                  }

#ifdef SEAS_DIAG_FAULT_FLUX
                  // C-1s STRESS-ROT (ADER shared-fault path): same trace
                  // as ADER-INT above, tagged ADER-SHR for distinguishing
                  // when the hypocenter QP sits on a rank-boundary
                  // shared face.  Print Q_avg = I/dt.
                  if (dof_idx >= 0 &&
                      dof_idx < static_cast<int>(fault_dof_data_->size()) &&
                      (*fault_dof_data_)[dof_idx].diag_print &&
                      dt > 0)
                  {
                     const real_t inv_dt = 1.0 / dt;

                     std::fprintf(stderr,
                        "[C-1s ADER-SHR BASIS] rank=%d  "
                        "can_n=(%+.4e,%+.4e,%+.4e)  "
                        "can_t1=(%+.4e,%+.4e,%+.4e)  "
                        "can_t2=(%+.4e,%+.4e,%+.4e)  "
                        "sign_flipped=%d  elem1_on_plus=%d  dt=%.4e\n",
                        g_seas_my_rank,
                        can_n[0], can_n[1], can_n[2],
                        can_t1[0], can_t1[1], can_t1[2],
                        can_t2[0], can_t2[1], can_t2[2],
                        qpd.sign_flipped ? 1 : 0,
                        elem1_on_plus  ? 1 : 0,
                        dt);

                     auto _c1s_print = [&](const char *tag, const real_t *I)
                     {
                        std::fprintf(stderr,
                           "[C-1s ADER-SHR %s] rank=%d  "
                           "SXX=%+.4e SYY=%+.4e SZZ=%+.4e  "
                           "SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
                           "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
                           tag, g_seas_my_rank,
                           I[SXX]*inv_dt, I[SYY]*inv_dt, I[SZZ]*inv_dt,
                           I[SXY]*inv_dt, I[SYZ]*inv_dt, I[SXZ]*inv_dt,
                           I[VX]*inv_dt,  I[VY]*inv_dt,  I[VZ]*inv_dt);
                     };

                     const real_t *I_glob_p = elem1_on_plus ? I_self : I_nbr;
                     const real_t *I_glob_m = elem1_on_plus ? I_nbr  : I_self;

                     _c1s_print("GLOB-IN+",  I_glob_p);
                     _c1s_print("GLOB-IN-",  I_glob_m);
                     _c1s_print("LOC-IN+",   I_plus_local);
                     _c1s_print("LOC-IN-",   I_minus_local);
                     _c1s_print("LOC-IMP+",  I_imp_plus);
                     _c1s_print("LOC-IMP-",  I_imp_minus);
                     _c1s_print("GLOB-OUT+", I_imp_plus_g);
                     _c1s_print("GLOB-OUT-", I_imp_minus_g);
                  }
#endif

                  real_t F_h_side[NUM_STATE];
                  const real_t *I_imp_side = elem1_on_plus
                                             ? I_imp_plus_g
                                             : I_imp_minus_g;
                  flux_.Interior(can_n, I_imp_side, I_imp_side, F_h_side);

                  const real_t assemble_sign = elem1_on_plus ? -1.0 : +1.0;
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     for (int i = 0; i < ndof; i++)
                     {
                        rhs[c * ndof_total_ + dof_offset1 + i] +=
                           assemble_sign * w * shape1(i) * F_h_side[c];
                     }
                  }
                  continue;
               }
               else
               {
                  MFEM_ABORT("ComputeADERSharedFaceFluxRHS: shared fault "
                             "face sf=" << sf << " missing DOFData.");
               }
            }
            else if (use_precomputed_face_fluxes_)
            {
               // Phase 2b (§7.2): one-sided precomputed interior dispatch
               // for shared non-fault faces.  Uses the LEGACY
               // AddInteriorFaceRhs path rather than Arm 3d's Full path,
               // because the Full path reads I_nbr from the LOCAL Q array
               // at dof_offset_nbr, and on shared faces the neighbor DOFs
               // live in the ghost buffer (pfes->GetFaceNbrFE /
               // q_gf.FaceNbrData) at a different offset.  The legacy
               // path already handles the ghost-buffer sampling via
               // MFEM's shape2 on the face-nbr FE.  Arm 3d's shared-face
               // extension (Full_Shared with Q_local + Q_ghost) is
               // deferred — not required for the Phase 3 audit (which
               // runs serial).
               precomputed_face_fluxes_.AddInteriorFaceRhs(
                  mesh_face_idx, /*caller_elem=*/e1,
                  I_self, I_nbr, w, shape1.GetData(),
                  ndof, dof_offset1, ndof_total_, rhs);
               continue;
            }
            else
            {
               // Round-14C: optional n↔-n symmetrization for shared
               // non-fault faces, parallel to the Round-13C hook in
               // ComputeADERFaceFluxRHS.  Env unset → verbatim
               // pre-patch path.
               //
               // Round-11 Mixed-Flux dispatch (R-1206 + R-1207): central
               // flux on shared non-fault faces if their global mesh face
               // index is in central_flux_face_set_ (the set construction
               // honors shared seams via shared_mesh_face_set_).
               const bool use_central_here_shared =
                  (mf_on && central_flux_face_set_.count(mesh_face_idx) > 0);
               auto interior_or_central_shared = [&](const real_t *n_arg,
                                                     const real_t *Q_a,
                                                     const real_t *Q_b,
                                                     real_t *F_out)
               {
                  if (use_central_here_shared)
                     { flux_.Central(n_arg, Q_a, Q_b, F_out); }
                  else
                     { flux_.Interior(n_arg, Q_a, Q_b, F_out); }
               };

               // R-1411: nonfault_both_sym dead else-branch removed
               // (shared-face mirror of the local-face cleanup).
               interior_or_central_shared(nor, I_self, I_nbr, F_h);
            }

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
// ADER I-05 Phase 6: AdvanceADER — one-step predictor-corrector.
// ---------------------------------------------------------------------------
// Flow:
//   1. I = ∫_0^dt Q(t+τ) dτ  via ComputeADERTimeIntegrated (Phase 3).
//   2. rhs = ∫ ∇φ · A · I dV − ∫ φ · F_h(I) dA − [PML · I if enabled].
//   3. rhs *= M^{-1}.
//   4. Q_new = Q + rhs.
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::AdvanceADER(const Vector &Q, real_t dt,
                                         int order, Vector &Q_new) const
{
   MFEM_VERIFY(dt > 0.0,
               "AdvanceADER: dt must be > 0, got " << dt);
   MFEM_VERIFY(order >= 2 && order <= 4,
               "AdvanceADER: order must be in {2,3,4}, got " << order);
   MFEM_VERIFY(Q.Size() == NUM_STATE * ndof_total_,
               "AdvanceADER: Q size mismatch");
   MFEM_VERIFY(&Q != &Q_new,
               "AdvanceADER: Q and Q_new must be distinct Vectors");

   // R-1506: hoist Q_new sizing to the top of the routine so any future
   // path that wants to pre-stage state into Q_new (e.g., a one-shot
   // corrector-only entry, an inline driver-side rk-stage tee, or a
   // PML-on-Q_new variant) finds Q_new pre-sized rather than zero-sized.
   // Pre-R-1506 the SetSize ran AFTER all corrector work — benign today
   // because the PML branch writes only `rhs`, but the contract was
   // fragile.  No semantic change at the current call sites.
   Q_new.SetSize(NUM_STATE * ndof_total_);

   // R-1501: lazy-init member scratch buffers (ader_I_buf_, ader_rhs_buf_)
   // instead of stack-allocating fresh Vectors on every macro-step.  See
   // header docstring for the cost rationale (~32 TB malloc traffic
   // saved on a 60 s TPV104 run).  SetSize() is a no-op when the size
   // already matches; the first call sizes once.
   const int N = NUM_STATE * ndof_total_;
   if (ader_I_buf_.Size()   != N) { ader_I_buf_.SetSize(N); }
   if (ader_rhs_buf_.Size() != N) { ader_rhs_buf_.SetSize(N); }
   Vector &I   = ader_I_buf_;
   Vector &rhs = ader_rhs_buf_;

   // 1. CK predictor.
   ComputeADERTimeIntegrated(Q, dt, order, I);

   // 2. Corrector RHS accumulation.
   rhs = 0.0;

   ComputeADERVolumeUpdate(I, rhs);
   ComputeADERFaceFluxRHS(I, dt, rhs);

   if constexpr (IsParallelMesh<MeshType>::value)
   {
      ComputeADERSharedFaceFluxRHS(I, dt, rhs);
   }

   if (pml_layer_)
   {
      // PML damps (I - dt·Q_bg) toward zero.  Under total-Q, the background
      // is also time-integrated, so the damping target is dt·Q_bg.
      // Re-use ApplyPMLDamping on I: it internally subtracts bulk_bg_ from
      // Q_qp; we need it to subtract dt·bulk_bg_ instead.  Implement inline
      // to avoid fighting the Mult-centric signature.
      const real_t *I_data_pml = I.GetData();
      for (int e = 0; e < ne_; e++)
      {
         const FiniteElement *fe = fes_->GetFE(e);
         ElementTransformation *Tr = fes_->GetElementTransformation(e);
         int ndof = fe->GetDof();
         int dof_offset = e * ndof_per_el_;

         const IntegrationRule &ir = IntRules.Get(
            fe->GetGeomType(), 2*order_);

         Vector shape(ndof);
         for (int qp = 0; qp < ir.GetNPoints(); qp++)
         {
            const IntegrationPoint &ip = ir.IntPoint(qp);
            Tr->SetIntPoint(&ip);
            real_t w = ip.weight * Tr->Weight();

            fe->CalcShape(ip, shape);
            Vector phys(3);
            Tr->Transform(ip, phys);

            real_t dx, dy, dz;
            pml_layer_->ComputeDamping(phys(0), phys(1), phys(2), dx, dy, dz);
            if (dx <= 0.0 && dy <= 0.0 && dz <= 0.0) { continue; }

            real_t I_qp[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               I_qp[c] = 0.0;
               for (int i = 0; i < ndof; i++)
               {
                  I_qp[c] += shape(i) * I_data_pml[c * ndof_total_ + dof_offset + i];
               }
            }
            // Round-6 R-002: total-Q only; has_bulk_bg_ already asserted
            // in ComputeADERFaceFluxRHS earlier in AdvanceADER.
            for (int c = 0; c < NUM_STATE; c++)
            {
               I_qp[c] -= dt * bulk_bg_[c];
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
                     rhs[c * ndof_total_ + dof_offset + i] -=
                        w * d_c * shape(i) * I_qp[c];
                  }
               }
            }
         }
      }
   }

#ifdef SEAS_DIAG_FAULT_FLUX
   // C-2A BULK-DOF (ADER path): per-AdvanceADER-call dump of bulk Q
   // (input) and the corrector rhs at the hypocenter face's two
   // adjacent tets' diag DOFs.  rhs here is the time-INTEGRATED
   // corrector residual (volume + non-fault flux + fault flux,
   // pre-mass-inverse, pre-Q+rhs accumulation).
   if (diag_elem_plus_ >= 0 && diag_elem_minus_ >= 0 &&
       diag_face_dof_plus_ >= 0 && diag_face_dof_minus_ >= 0)
   {
      const real_t *Qd = Q.GetData();
      const real_t *Rd = rhs.GetData();
      const int dof_off_p = diag_elem_plus_  * ndof_per_el_;
      const int dof_off_m = diag_elem_minus_ * ndof_per_el_;
      const int idx_p = diag_face_dof_plus_;
      const int idx_m = diag_face_dof_minus_;
      if (idx_p < ndof_per_el_ && idx_m < ndof_per_el_)
      {
         real_t Qp[NUM_STATE], Qm[NUM_STATE], Rp[NUM_STATE], Rm[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Qp[c] = Qd[c*ndof_total_ + dof_off_p + idx_p];
            Qm[c] = Qd[c*ndof_total_ + dof_off_m + idx_m];
            Rp[c] = Rd[c*ndof_total_ + dof_off_p + idx_p];
            Rm[c] = Rd[c*ndof_total_ + dof_off_m + idx_m];
         }
         std::fprintf(stderr,
            "[C-2A BULK-DOF Q+]  rank=%d e=%d  "
            "SXX=%+.4e SYY=%+.4e SZZ=%+.4e SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
            "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
            g_seas_my_rank, diag_elem_plus_,
            Qp[SXX], Qp[SYY], Qp[SZZ], Qp[SXY], Qp[SYZ], Qp[SXZ],
            Qp[VX],  Qp[VY],  Qp[VZ]);
         std::fprintf(stderr,
            "[C-2A BULK-DOF Q-]  rank=%d e=%d  "
            "SXX=%+.4e SYY=%+.4e SZZ=%+.4e SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
            "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
            g_seas_my_rank, diag_elem_minus_,
            Qm[SXX], Qm[SYY], Qm[SZZ], Qm[SXY], Qm[SYZ], Qm[SXZ],
            Qm[VX],  Qm[VY],  Qm[VZ]);
         std::fprintf(stderr,
            "[C-2A BULK-DOF k+]  rank=%d e=%d  "
            "SXX=%+.4e SYY=%+.4e SZZ=%+.4e SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
            "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
            g_seas_my_rank, diag_elem_plus_,
            Rp[SXX], Rp[SYY], Rp[SZZ], Rp[SXY], Rp[SYZ], Rp[SXZ],
            Rp[VX],  Rp[VY],  Rp[VZ]);
         std::fprintf(stderr,
            "[C-2A BULK-DOF k-]  rank=%d e=%d  "
            "SXX=%+.4e SYY=%+.4e SZZ=%+.4e SXY=%+.4e SYZ=%+.4e SXZ=%+.4e  "
            "VX=%+.4e VY=%+.4e VZ=%+.4e\n",
            g_seas_my_rank, diag_elem_minus_,
            Rm[SXX], Rm[SYY], Rm[SZZ], Rm[SXY], Rm[SYZ], Rm[SXZ],
            Rm[VX],  Rm[VY],  Rm[VZ]);
      }
   }
#endif

   // 3. rhs *= M^{-1}.
   ApplyMassInverse(rhs);

   // 4. Q_new = Q + rhs.  (Q_new sized at top of routine, R-1506.)
   add(Q, 1.0, rhs, Q_new);
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
   // R-1403 + R-1502: central flux is non-dissipative; the CFL stability
   // factor depends on the FRACTION of faces using central vs upwind,
   // not just whether mixed-flux is engaged.  Zhang 2023 §3.3 cites
   // CFL=0.3 for mixed-flux runs vs CFL=0.5 for pure upwind; the ratio
   // that applies to a given mode depends on its central-face density:
   //
   //   Adjacent      (~5-10% central): operator dominantly upwind →
   //                                    CFL near pure upwind.  0.9×
   //                                    is a slight guard band.
   //   AllContinuous (~95% central):    operator nearly non-dissipative
   //                                    → explicit RK4/ADER stability
   //                                    significantly tighter.  0.4×
   //                                    is close to Zhang's 0.3-equiv.
   //
   // These factors are interim placeholders pending a multi-step
   // stability calibration on the production fixture.  Drivers that
   // calibrate CFL externally can compensate via the `cfl` argument.
   real_t cfl_mixed_flux_factor = 1.0;
   switch (mixed_flux_mode_)
   {
      case MixedFluxMode::None:          cfl_mixed_flux_factor = 1.0; break;
      case MixedFluxMode::Adjacent:      cfl_mixed_flux_factor = 0.9; break;
      case MixedFluxMode::AllContinuous: cfl_mixed_flux_factor = 0.4; break;
      default:
         // R-1600: any future MixedFluxMode value reaches this default
         // and ABORTS — silent fall-through (with cfl_factor=1.0, i.e.,
         // full upwind CFL) would destabilize multi-step production.
         // Adding a new mode REQUIRES adding the corresponding factor
         // here.
         MFEM_ABORT("ComputeMaxDt: unknown MixedFluxMode "
                    << static_cast<int>(mixed_flux_mode_)
                    << " (R-1600 fall-through guard).  Adding a new "
                    "MixedFluxMode REQUIRES adding the corresponding "
                    "CFL factor in this switch.");
   }

   // Phase H.3 (Stage 1): heterogeneous-CFL path.  When the
   // heterogeneous ctor populated `per_elem_h_` / `per_elem_lmr_`,
   // walk every owned element to compute
   //   dt_e = cfl_factor * cfl * h_e / c_p,e
   // and `MPI_Allreduce(MIN)` across ranks.
   //
   // For Mode::Constant input (Stage 1's only supported mode) every
   // element shares the same (lambda, mu, rho), so c_p,e is constant
   // and the local min over `cfl_factor * cfl * h_e / c_p` is exactly
   // `cfl_factor * cfl * (min_e h_e) / c_p`.  The MPI MIN reduction
   // over these per-rank locals matches the legacy path's reduction
   // on `h_min_` (which was reduced over the same per-element `h_e`
   // set in the scalar ctor).  So this branch is BYTE-IDENTICAL to
   // the scalar formula below — proven by
   // `test_phaseh_wave_operator_constant_parity` Test C-3.
   if (!per_elem_h_.empty())
   {
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

   return cfl_mixed_flux_factor * cfl * h_min_ / flux_.GetCp();
}

// ---------------------------------------------------------------------------
// PML damping: rhs -= d(x) * D * (Q - Q_bg)  (Eq. 16)
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::ApplyPMLDamping(const Vector &Q, Vector &rhs) const
{
   // Round-6 R-002: Total-Q only.  PML damps toward the bulk background
   // Q_bg supplied via SetAbsorbingBackground; Q_bg=0 is valid
   // (reproduces the old damping-toward-zero behaviour) but
   // SetAbsorbingBackground must have been called.
   MFEM_VERIFY(has_bulk_bg_,
               "Total-Q only: SetAbsorbingBackground(Q_bg) must be called "
               "before enabling PML and running wave.Mult().");

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

         // Round-6 R-002: Total-Q only.  Damp the FLUCTUATION
         // (Q - Q_bg) toward zero, not Q itself — damping Q toward zero
         // would erode the pre-stress tensor in the PML region and
         // create a stress-gradient artefact at the PML/interior
         // interface.  has_bulk_bg_ is asserted in ComputeFaceFluxRHS
         // upstream of this call in Mult.
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_qp[c] -= bulk_bg_[c];
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
   real_t tol, double *worst_rel_out, int *worst_field_out,
   bool abort_on_fail) const
{
   // R-DIP1: default the monitoring out-params so the serial / no-shared-face
   // short-circuits below report "no divergence" rather than leaving them
   // uninitialised.
   if (worst_rel_out)   { *worst_rel_out = 0.0; }
   if (worst_field_out) { *worst_field_out = -1; }

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

      // Pack per-QP records with BP5-style face-vertex-key for cross-rank
      // matching.  One record (16 doubles):
      //
      //   [0..2]  face_key (sorted triple of HYPRE_BigInt global vertex IDs)
      //           encoded as double — bit-exact for |v| < 2^53 (holds for any
      //           realistic mesh; TPV102-200m has ~500K vertices).
      //   [3]     qp_idx within the face (0, 1, or 2 for triangles at order=1).
      //   [4..6]  centroid cx, cy, cz — diagnostic only, not used for pairing.
      //   [7..14] 8 mutable DOFData fields that
      //           FaultFaceFlux::Evaluate / the RK4 averaging step writes.
      //   [15]    emitting rank.
      //
      // v8.0.0 Phase 1 / job 7666233 root cause: the previous centroid-based
      // pairing (sort by (cx,cy,cz) with a tolerance) is sensitive to
      // sub-ULP FP drift between the two ranks' `ftr->Face->Transform` paths
      // for the same shared face, as MFEM's shared-face code does not
      // guarantee bit-identical physical coordinates across ranks for every
      // QP.  Replacing the key with the sorted global-vertex-ID triple
      // (the same approach BP5's `ElasticityOperator` uses at scale —
      // `MakeFaceKey` in `domain/elasticity_operator.hpp`, mirrored locally
      // in `dynamic/shared_fault_key.hpp`) makes pairing integer-exact and
      // immune to FP drift.  The (key, qp_idx) composite key pairs the same
      // physical QP on both ranks deterministically: IntRules is deterministic
      // given the same geometry type, so rank A's qp_idx=k and rank B's
      // qp_idx=k are the same reference point on the same canonical face.
      constexpr int REC = 16;
      constexpr int KEY_OFFSET = 0;
      constexpr int NUM_KEY_COMPONENTS = 3;
      constexpr int QP_IDX_OFFSET = 3;
      constexpr int CENTROID_OFFSET = 4;
      constexpr int FIELD_BASE = 7;
      constexpr int NUM_FIELDS = 8;
      constexpr int RANK_OFFSET = 15;
      static const char *FIELD_NAMES[NUM_FIELDS] = {
         "tau1_corr", "tau2_corr", "sigma_n_corr",
         "V1", "V2", "psi", "slip1", "slip2"
      };

      // Need the global vertex index table to build face keys.
      Array<HYPRE_BigInt> gvi;
      pmesh.GetGlobalVertexIndices(gvi);

      std::vector<double> local_data;

      for (int sf_idx = 0; sf_idx < fault_shared_faces_.Size(); sf_idx++)
      {
         int sf = fault_shared_faces_[sf_idx];
         FaceElementTransformations *ftr = pmesh.GetSharedFaceTransformations(sf);
         if (!ftr) { continue; }

         auto it = shared_fault_dof_offset_.find(sf);
         if (it == shared_fault_dof_offset_.end()) { continue; }
         int dof_base = it->second;

         int local_face = pmesh.GetSharedFace(sf);
         dynamic::FaceVertexKey key =
            dynamic::MakeFaceKey(local_face, gvi, pmesh);

         const IntegrationRule &ir =
            IntRules.Get(ftr->GetGeometryType(), 2*order_);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            Vector phys(3);
            ftr->Face->Transform(ip, phys);

            int didx = dof_base + q;
            if (didx < 0 ||
                didx >= static_cast<int>(fault_dof_data_->size())) { continue; }
            const DOFData &d = (*fault_dof_data_)[didx];

            local_data.push_back(static_cast<double>(key.v[0]));
            local_data.push_back(static_cast<double>(key.v[1]));
            local_data.push_back(static_cast<double>(key.v[2]));
            local_data.push_back(static_cast<double>(q));
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

      // Pairing strategy (v8.0.0 Phase 2 refinement, np=14 DIAG result):
      //
      // Primary key: face vertex triple (integer, exact) — groups the 6
      //   entries for one shared fault face (3 QPs × 2 owning ranks) into
      //   one contiguous block.  No tolerance, no FP pathology.
      //
      // Secondary key within a face group: physical centroid (with ~1e-6 m
      //   tolerance) — this is necessary because MFEM's shared-face
      //   orientation can flip between the two ranks owning the face
      //   (observed in DIAG output: rank 5 had `sign_flipped=1`, rank 6
      //   `sign_flipped=0`).  When that happens, rank A's qp_idx=k and
      //   rank B's qp_idx=k are at DIFFERENT physical points on the face,
      //   and pairing by qp_idx alone compares mismatched physical QPs
      //   (whose psi/V differ by the spatial-gradient over face_size/3 —
      //   the spurious "2.19e-2 psi drift" of Step 1.11).
      //
      // Within a face, the 3 QPs are spatially well-separated (~ face_size /
      //   3 ≈ 300 m on the 1000 m mesh), while the two ranks' views of the
      //   SAME physical QP agree to ~1 ULP (1e-12 at magnitude 1e4).  A
      //   1e-6 m tolerance on the centroid secondary key therefore
      //   unambiguously identifies same-physical-QP pairs without risk of
      //   false grouping.  qp_idx is ignored for pairing but kept in the
      //   record + diagnostic output so orientation-flip cases are visible.
      auto key_component = [&](int entry, int k) -> HYPRE_BigInt
      {
         return static_cast<HYPRE_BigInt>(all_data[entry*REC + k]);
      };

      std::vector<int> idx(n_entries);
      std::iota(idx.begin(), idx.end(), 0);
      std::sort(idx.begin(), idx.end(), [&](int a, int b)
      {
         // Primary: face_key (3 int64 components, exact).
         for (int k = 0; k < NUM_KEY_COMPONENTS; k++)
         {
            HYPRE_BigInt va = key_component(a, k);
            HYPRE_BigInt vb = key_component(b, k);
            if (va != vb) { return va < vb; }
         }
         // Secondary: physical centroid lex (tie-breaker within a face).
         for (int k = 0; k < 3; k++)
         {
            double va = all_data[a*REC + CENTROID_OFFSET + k];
            double vb = all_data[b*REC + CENTROID_OFFSET + k];
            if (va != vb) { return va < vb; }
         }
         return false;
      });

      auto same_face_qp = [&](int a, int b)
      {
         // Must share face (exact integer match)...
         for (int k = 0; k < NUM_KEY_COMPONENTS; k++)
         {
            if (key_component(a, k) != key_component(b, k)) { return false; }
         }
         // ...AND have centroids agreeing to 1e-6 m (safe within-face tol).
         constexpr double centroid_abs_floor = 1e-6;
         for (int k = 0; k < 3; k++)
         {
            double va = all_data[a*REC + CENTROID_OFFSET + k];
            double vb = all_data[b*REC + CENTROID_OFFSET + k];
            double scale = std::max(std::abs(va), std::abs(vb));
            double tol_k = std::max(scale
                                    * std::numeric_limits<double>::epsilon(),
                                    centroid_abs_floor);
            if (std::abs(va - vb) > tol_k) { return false; }
         }
         return true;
      };

      double max_rel_diff = 0.0;
      double max_diff_abs = 0.0;
      double max_diff_scale = 1.0;
      int    max_diff_field = -1;
      int    max_diff_entry = -1;
      int    n_pairs = 0;
      int    n_unpaired = 0;

      // Unpaired entries still reported so a genuine classification
      // pathology (one rank doesn't see this face as fault; global-key
      // mismatch) stays visible in the abort trace.  With integer-exact
      // pairing, tolerance-induced false singletons can no longer happen.
      struct UnpairedEntry
      {
         HYPRE_BigInt key[3];
         int          qp_idx;
         double       cx, cy, cz;
         int          rank;
         int          group_size;
      };
      std::vector<UnpairedEntry> unpaired;
      constexpr int MAX_UNPAIRED_REPORT = 32;

      auto push_unpaired = [&](int e, int group_size)
      {
         if (static_cast<int>(unpaired.size()) >= MAX_UNPAIRED_REPORT) { return; }
         UnpairedEntry u;
         u.key[0] = key_component(e, 0);
         u.key[1] = key_component(e, 1);
         u.key[2] = key_component(e, 2);
         u.qp_idx = static_cast<int>(key_component(e, QP_IDX_OFFSET));
         u.cx = all_data[e*REC + CENTROID_OFFSET + 0];
         u.cy = all_data[e*REC + CENTROID_OFFSET + 1];
         u.cz = all_data[e*REC + CENTROID_OFFSET + 2];
         u.rank = static_cast<int>(all_data[e*REC + RANK_OFFSET]);
         u.group_size = group_size;
         unpaired.push_back(u);
      };

      // v8.0.0 Phase 3 (job 7666323 diagnosis): the sequential-adjacency
      // grouping loop is fragile to sort-order accidents.  On the 50-rank
      // 200 m topology, one shared fault face has two vertices V0 and V1
      // with bit-identical x — this makes rank 6's q=1 (at P1) and q=2
      // (at P0) compute to bit-identical cx while rank 3's q=0 (also at
      // P0) has a 1-ULP smaller cx.  The sort then places rank 3's P0
      // entry ALONE (smaller cx) and the three rank-6 entries + rank-3's
      // P1 together at the larger cx.  `same_face_qp` between two
      // physical-P0 entries is TRUE, but they are non-adjacent, so the
      // sequential grouping loop misses the pair.
      //
      // Fix: within a face group, do all-pairs nearest-neighbor matching
      // by centroid (O(n²) per face, n=6 typical, cheap).  Robust to any
      // sort artifact, including orientation-flip + shared-vertex-x
      // combined edge cases.  `same_face_qp` tolerance check becomes a
      // PAIR TOLERANCE (must hold for a claimed pair), not a GROUP
      // predicate.
      auto centroid_max_dist = [&](int a, int b) -> double
      {
         double dx = std::abs(all_data[a*REC + CENTROID_OFFSET + 0]
                            - all_data[b*REC + CENTROID_OFFSET + 0]);
         double dy = std::abs(all_data[a*REC + CENTROID_OFFSET + 1]
                            - all_data[b*REC + CENTROID_OFFSET + 1]);
         double dz = std::abs(all_data[a*REC + CENTROID_OFFSET + 2]
                            - all_data[b*REC + CENTROID_OFFSET + 2]);
         return std::max({dx, dy, dz});
      };

      auto same_face_key = [&](int a, int b)
      {
         for (int k = 0; k < NUM_KEY_COMPONENTS; k++)
         {
            if (key_component(a, k) != key_component(b, k)) { return false; }
         }
         return true;
      };

      constexpr double PAIR_TOL_M = 1e-6;

      int i = 0;
      while (i < n_entries)
      {
         // Find end of face-key group (contiguous after the face_key
         // primary sort).
         int j = i + 1;
         while (j < n_entries && same_face_key(idx[i], idx[j])) { j++; }
         int n_in_group = j - i;

         // All-pairs nearest-neighbor matching within this face group.
         // Repeatedly pick the closest unmatched pair whose ranks differ
         // and whose centroid separation is below PAIR_TOL_M; pair them,
         // compare fields, mark matched.  Stop when no more valid pairs.
         std::vector<bool> matched(n_in_group, false);
         while (true)
         {
            double best_d = std::numeric_limits<double>::max();
            int best_a = -1, best_b = -1;
            for (int a = 0; a < n_in_group; a++)
            {
               if (matched[a]) { continue; }
               int ea = idx[i + a];
               int rank_a = static_cast<int>(all_data[ea*REC + RANK_OFFSET]);
               for (int b = a + 1; b < n_in_group; b++)
               {
                  if (matched[b]) { continue; }
                  int eb = idx[i + b];
                  int rank_b = static_cast<int>(
                     all_data[eb*REC + RANK_OFFSET]);
                  if (rank_a == rank_b)
                  {
                     // Same-rank entries are not a valid cross-rank pair.
                     continue;
                  }
                  double d = centroid_max_dist(ea, eb);
                  if (d < best_d)
                  {
                     best_d = d;
                     best_a = a;
                     best_b = b;
                  }
               }
            }
            if (best_a < 0 || best_d > PAIR_TOL_M) { break; }

            matched[best_a] = true;
            matched[best_b] = true;
            int ea = idx[i + best_a];
            int eb = idx[i + best_b];
            for (int k = FIELD_BASE; k < FIELD_BASE + NUM_FIELDS; k++)
            {
               double va = all_data[ea*REC + k];
               double vb = all_data[eb*REC + k];
               double diff = std::abs(va - vb);
               double field_scale = std::max({std::abs(va), std::abs(vb),
                                              1.0});
               double rel_diff = diff / field_scale;
               if (rel_diff > max_rel_diff)
               {
                  max_rel_diff = rel_diff;
                  max_diff_abs = diff;
                  max_diff_scale = field_scale;
                  max_diff_field = k;
                  max_diff_entry = ea;
               }
            }
            n_pairs++;
         }

         // Unmatched entries in this face group are singletons.  Genuine
         // classification asymmetry (one rank doesn't classify this face
         // as fault) or same-rank double emission (ctor dof_offset bug)
         // are both reported here.
         for (int a = 0; a < n_in_group; a++)
         {
            if (!matched[a])
            {
               push_unpaired(idx[i + a], 1);
               n_unpaired++;
            }
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

      // R-DIP1: report the worst pairing diff via the out-params (every rank
      // holds the full allgathered record set, so max_rel_diff/max_diff_field
      // are identical across ranks).  Set on BOTH success and failure.
      if (worst_rel_out)   { *worst_rel_out = max_rel_diff; }
      if (worst_field_out) { *worst_field_out = max_diff_field; }

      // R-DIP1: non-aborting diagnostic mode (tests) — return on mismatch
      // without the verbose dump / MFEM_ABORT.  Collective-safe: every rank
      // computed the same `fail_global`, so all return together.
      if (fail_global && !abort_on_fail) { return; }

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
                  // Print the face-vertex-key (integer, bit-exact across
                  // ranks) + qp_idx + physical centroid (diagnostic only).
                  // With integer-exact pairing, a group_size != 2 means
                  // either (a) only one rank classified this face as fault
                  // (classification asymmetry — likely a global-fault-key
                  // allgather bug) or (b) the same rank emitted twice for
                  // the same (face, qp) (dof-offset double mapping).
                  std::fprintf(stderr,
                     "  [UNPAIRED] face_key=(%lld, %lld, %lld) qp_idx=%d "
                     "centroid=(%.17e, %.17e, %.17e) rank=%d group_size=%d\n",
                     static_cast<long long>(unpaired[u].key[0]),
                     static_cast<long long>(unpaired[u].key[1]),
                     static_cast<long long>(unpaired[u].key[2]),
                     unpaired[u].qp_idx,
                     unpaired[u].cx, unpaired[u].cy, unpaired[u].cz,
                     unpaired[u].rank, unpaired[u].group_size);
               }
               std::fprintf(stderr,
                  "  Likely fault-classification asymmetry (one rank's "
                  "ctor did not add this face-key to fault_shared_faces_ — "
                  "check global_fault_keys allgather) OR ctor "
                  "dof-offset double mapping.\n\n");

               // v8.0.0 Phase 2 follow-up diagnostic: when entries remain
               // unpaired under the Phase 2 (face_key, centroid) predicate,
               // dump ALL entries that share a face_key with any orphan.
               // This reveals the full 6-entry picture for the affected
               // face so we can see which QP pairs succeeded and why the
               // orphaned pair fell through.  Also prints
               // `same_face_qp(orphan_i, orphan_j)` for every orphan
               // pairing to distinguish (a) predicate returning FALSE
               // unexpectedly from (b) algorithm grouping them apart
               // despite predicate TRUE (e.g., intruder entry breaks
               // adjacency).
               std::fprintf(stderr,
                  "  [R-101 face-context dump for unpaired entries]\n");
               std::set<std::array<HYPRE_BigInt,3>> unpaired_keys;
               for (size_t u = 0; u < unpaired.size(); u++)
               {
                  unpaired_keys.insert({unpaired[u].key[0],
                                        unpaired[u].key[1],
                                        unpaired[u].key[2]});
               }
               for (const auto &k : unpaired_keys)
               {
                  std::fprintf(stderr,
                     "    face_key=(%lld, %lld, %lld) all entries in sort "
                     "order:\n",
                     static_cast<long long>(k[0]),
                     static_cast<long long>(k[1]),
                     static_cast<long long>(k[2]));
                  for (int p = 0; p < n_entries; p++)
                  {
                     int e = idx[p];
                     HYPRE_BigInt k0 = key_component(e, 0);
                     HYPRE_BigInt k1 = key_component(e, 1);
                     HYPRE_BigInt k2 = key_component(e, 2);
                     if (k0 != k[0] || k1 != k[1] || k2 != k[2]) { continue; }
                     int er = static_cast<int>(all_data[e*REC + RANK_OFFSET]);
                     int qidx = static_cast<int>(
                        key_component(e, QP_IDX_OFFSET));
                     std::fprintf(stderr,
                        "      sort_pos=%d qp_idx=%d rank=%d "
                        "centroid=(%.17e, %.17e, %.17e)\n",
                        p, qidx, er,
                        all_data[e*REC + CENTROID_OFFSET + 0],
                        all_data[e*REC + CENTROID_OFFSET + 1],
                        all_data[e*REC + CENTROID_OFFSET + 2]);
                  }
               }
               // Direct same_face_qp() check between every orphan pair that
               // shares a face_key — this pinpoints whether the predicate
               // is the cause or the sort/group order is.
               if (unpaired.size() >= 2)
               {
                  std::fprintf(stderr,
                     "  [R-101 same_face_qp self-check between orphans]\n");
                  for (size_t u = 0; u < unpaired.size(); u++)
                  {
                     for (size_t v = u + 1; v < unpaired.size(); v++)
                     {
                        if (unpaired[u].key[0] != unpaired[v].key[0] ||
                            unpaired[u].key[1] != unpaired[v].key[1] ||
                            unpaired[u].key[2] != unpaired[v].key[2])
                        { continue; }
                        // Find the data indices by linear search.
                        int idx_u = -1, idx_v = -1;
                        for (int e = 0; e < n_entries; e++)
                        {
                           int er = static_cast<int>(
                              all_data[e*REC + RANK_OFFSET]);
                           if (idx_u < 0 && er == unpaired[u].rank &&
                               all_data[e*REC + CENTROID_OFFSET + 0] ==
                                  unpaired[u].cx)
                           { idx_u = e; }
                           if (idx_v < 0 && er == unpaired[v].rank &&
                               all_data[e*REC + CENTROID_OFFSET + 0] ==
                                  unpaired[v].cx)
                           { idx_v = e; }
                        }
                        if (idx_u < 0 || idx_v < 0)
                        {
                           std::fprintf(stderr,
                              "    could not re-locate orphan pair u=%zu "
                              "v=%zu\n", u, v);
                           continue;
                        }
                        bool sc = same_face_qp(idx_u, idx_v);
                        double dx = std::abs(all_data[idx_u*REC+CENTROID_OFFSET+0]
                                            - all_data[idx_v*REC+CENTROID_OFFSET+0]);
                        double dy = std::abs(all_data[idx_u*REC+CENTROID_OFFSET+1]
                                            - all_data[idx_v*REC+CENTROID_OFFSET+1]);
                        double dz = std::abs(all_data[idx_u*REC+CENTROID_OFFSET+2]
                                            - all_data[idx_v*REC+CENTROID_OFFSET+2]);
                        std::fprintf(stderr,
                           "    same_face_qp(orphan u=%zu rank=%d, "
                           "v=%zu rank=%d) = %s "
                           "(dx=%.3e dy=%.3e dz=%.3e)\n",
                           u, unpaired[u].rank, v, unpaired[v].rank,
                           sc ? "TRUE" : "FALSE", dx, dy, dz);
                     }
                  }
               }
               std::fflush(stderr);
            }
            MFEM_ABORT("R-101 shared-fault DOFData: " << n_unpaired
                       << " unpaired entries (every shared QP should have "
                       "exactly 2 ranks).  See rank-0 [R-101 rank-0 detail] "
                       "lines above this abort trace.");
         }
         double cx = all_data[max_diff_entry*REC + CENTROID_OFFSET + 0];
         double cy = all_data[max_diff_entry*REC + CENTROID_OFFSET + 1];
         double cz = all_data[max_diff_entry*REC + CENTROID_OFFSET + 2];
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
