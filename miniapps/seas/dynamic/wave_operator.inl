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
// TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23
// Phase 2a (§6.3 + §6.3a): opt-in precomputed-flux switch with late
// `fault_face_set_` population (R4-001 + R5-002 FIXES).
// ---------------------------------------------------------------------------
template <typename MeshType>
void WaveOperator<MeshType>::UsePrecomputedFaceFluxes(bool enable)
{
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
      MFEM_ASSERT(ndof == ndof_per_el_,
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

   // Ping-pong buffers for the recursion.
   Vector D_curr(Q);                             // D(0) = Q
   Vector D_next(NUM_STATE * ndof_total_);
   Vector dQ_dxd(NUM_STATE * ndof_total_);

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

      // Advance factorial factor: fac *= dt / (k+2)
      fac *= dt / static_cast<real_t>(k + 2);
      I.Add(fac, D_next);

      // Swap buffers: D_curr <- D_next for the next iteration.
      mfem::Swap(D_curr, D_next);
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
   MFEM_VERIFY(has_bulk_bg_,
               "wave.AdvanceADER() requires SetAbsorbingBackground(Q_bg) "
               "to have been called (Q_bg = 0 is valid under "
               "fluctuation-Q dispatch).");

   const real_t *I_data = I.GetData();

   // Time-integrated background: bulk_bg_scaled[c] = dt · bulk_bg_[c].
   real_t bulk_bg_scaled[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++)
   {
      bulk_bg_scaled[c] = dt * bulk_bg_[c];
   }

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

   // Round-13C Patch 1: test-only n ↔ -n symmetrization of the ADER
   // local non-fault branches (both interior non-fault and boundary).
   // When `SEAS_TEST_NONFAULT_BOTH_SYM=1`, the per-QP flux is replaced
   // by 0.5 * (F(nor) + F(-nor)) at every non-fault face.  Mirrors
   // `ComputeInteriorFhSym` / `ComputeBoundaryFhSym` in
   // `test_adjacent_triangle_fault_first_step_audit.cpp`.  Fault faces
   // and shared faces are unaffected in this round.  Unset env →
   // byte-identical to the pre-patch path.
   const bool nonfault_both_sym = []()
   {
      const char *env = std::getenv("SEAS_TEST_NONFAULT_BOTH_SYM");
      return env && env[0] == '1';
   }();

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

                  const bool elem1_on_plus = !qpd.sign_flipped;

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
                  // v9.4.0 Commit 3: fluctuation-Q ADER dispatch;
                  // has_bulk_bg_ already asserted at the top of
                  // ComputeADERFaceFluxRHS (Q_bg = 0 is valid).
                  fault_flux_->EvaluateADER(fdata,
                                            I_plus_local, I_minus_local,
                                            dt,
                                            I_imp_plus, I_imp_minus);

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

                  const FaultBasisQPData &qpd = bd.qp_data[qq];
                  real_t can_t1[3], can_t2[3];
                  for (int d = 0; d < 3; d++)
                  {
                     can_n_per_qp[qq][d] = qpd.sign_flipped ? -qpd.normal[d]
                                                            :  qpd.normal[d];
                     can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d]
                                                  :  qpd.tangent1[d];
                     can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d]
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

                  elem1_on_plus_per_qp[qq] = !qpd.sign_flipped;

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

            if (!nonfault_both_sym)
            {
               compute_bc_flux(nor, F_h);
            }
            else
            {
               real_t F_pos[NUM_STATE], F_neg[NUM_STATE];
               real_t nor_neg[3] = {-nor[0], -nor[1], -nor[2]};
               compute_bc_flux(nor,     F_pos);
               compute_bc_flux(nor_neg, F_neg);
               for (int c = 0; c < NUM_STATE; c++)
               {
                  F_h[c] = 0.5 * (F_pos[c] + F_neg[c]);
               }
            }

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
            if (!nonfault_both_sym)
            {
               flux_.Interior(nor, I_self, I_nbr, F_h);
            }
            else
            {
               real_t F_pos[NUM_STATE], F_neg[NUM_STATE];
               real_t nor_neg[3] = {-nor[0], -nor[1], -nor[2]};
               flux_.Interior(nor,     I_self, I_nbr, F_pos);
               flux_.Interior(nor_neg, I_nbr,  I_self, F_neg);
               for (int c = 0; c < NUM_STATE; c++)
               {
                  F_h[c] = 0.5 * (F_pos[c] + F_neg[c]);
               }
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

      auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
      MFEM_VERIFY(pfes, "FESpace must be ParFiniteElementSpace for ParMesh");

      auto &pmesh = static_cast<const ParMesh &>(mesh_);
      int n_shared = pmesh.GetNSharedFaces();
      if (n_shared == 0) { return; }

      // v9.4.0: see ComputeADERFaceFluxRHS top-level note.  REVIEW R-007:
      // updated stale "Total-Q only" wording.
      MFEM_VERIFY(has_bulk_bg_,
                  "wave.AdvanceADER() requires SetAbsorbingBackground("
                  "Q_bg) to have been called (Q_bg = 0 is valid under "
                  "fluctuation-Q dispatch).");

      if (bc_.fault_attr > 0)
      {
         MFEM_VERIFY(fault_flux_ && fault_dof_data_,
                     "WaveOperator::ComputeADERSharedFaceFluxRHS: "
                     "bc_.fault_attr=" << bc_.fault_attr
                     << " > 0 but fault bookkeeping unset.");
      }

      // Round-14C: reuse SEAS_TEST_NONFAULT_BOTH_SYM for MPI/Frontera
      // parity with the local-face hook added in Round-13C Patch 1.
      // Applies only to shared NON-FAULT faces (the runtime else
      // branch below); shared fault faces and the precomputed path
      // are untouched.  Env unset → byte-identical to pre-patch.
      const bool nonfault_both_sym = []()
      {
         const char *env = std::getenv("SEAS_TEST_NONFAULT_BOTH_SYM");
         return env && env[0] == '1';
      }();

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

            MFEM_ASSERT(ndof2 == ndof_per_el_,
                        "Mixed element types in ghost");

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
                  // v9.4.0 Commit 3: fluctuation-Q ADER dispatch;
                  // has_bulk_bg_ asserted at the top of
                  // ComputeADERSharedFaceFluxRHS (Q_bg = 0 is valid).
                  fault_flux_->EvaluateADER(fdata,
                                            I_plus_local, I_minus_local,
                                            dt,
                                            I_imp_plus, I_imp_minus);

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
               if (!nonfault_both_sym)
               {
                  flux_.Interior(nor, I_self, I_nbr, F_h);
               }
               else
               {
                  real_t F_pos[NUM_STATE], F_neg[NUM_STATE];
                  real_t nor_neg[3] = {-nor[0], -nor[1], -nor[2]};
                  flux_.Interior(nor,     I_self, I_nbr, F_pos);
                  flux_.Interior(nor_neg, I_nbr,  I_self, F_neg);
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     F_h[c] = 0.5 * (F_pos[c] + F_neg[c]);
                  }
               }
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

   // 1. CK predictor.
   Vector I(NUM_STATE * ndof_total_);
   ComputeADERTimeIntegrated(Q, dt, order, I);

   // 2. Corrector RHS accumulation.
   Vector rhs(NUM_STATE * ndof_total_);
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

   // 3. rhs *= M^{-1}.
   ApplyMassInverse(rhs);

   // 4. Q_new = Q + rhs.
   Q_new.SetSize(NUM_STATE * ndof_total_);
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
   return cfl * h_min_ / flux_.GetCp();
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
