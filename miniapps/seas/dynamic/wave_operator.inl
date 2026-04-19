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
         shared_mesh_face_set_.clear();
         // R-001 fix: resolve peer rank per shared face via face_nbr_elements_offset.
         // R-107 fix: resolved/peer_rank fields are independent (a false
         // `resolved` means "ctor couldn't resolve this face", distinct from
         // "peer rank happens to be 0").
         shared_face_peer_.assign(n_shared, SharedFacePeer{});
         const Array<int> &fn_offs = pmesh.face_nbr_elements_offset;
         const int num_fn = pmesh.GetNFaceNeighbors();
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

            // Resolve peer rank: map the ghost element to its face neighbor
            // group, then look up that group's remote MPI rank.
            int nbr_elem_idx = ftr->Elem2No - ne_;
            int fn = 0;
            while (fn < num_fn && nbr_elem_idx >= fn_offs[fn + 1]) { fn++; }
            if (fn < num_fn)
            {
               shared_face_peer_[sf].resolved  = true;
               shared_face_peer_[sf].peer_rank = pmesh.GetFaceNbrRank(fn);
            }
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

               if (dof_idx >= 0 &&
                   dof_idx < static_cast<int>(fault_dof_data_->size()))
               {
                  DOFData &fdata = (*fault_dof_data_)[dof_idx];

                  // 1. Rotate Q± to fault-local coordinates
                  real_t t1[3], t2[3];
                  GodunovFlux::BuildFrame(nor, t1, t2);
                  DenseMatrix T(NUM_STATE), Tinv(NUM_STATE);
                  GodunovFlux::BuildRotation(nor, t1, t2, T);
                  GodunovFlux::BuildRotationInverse(nor, t1, t2, Tinv);

                  real_t Q_plus_local[NUM_STATE], Q_minus_local[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     Q_plus_local[c] = 0.0;
                     Q_minus_local[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        Q_plus_local[c]  += Tinv(c, k) * Q_self[k];
                        Q_minus_local[c] += Tinv(c, k) * Q_nbr[k];
                     }
                  }

                  // 2. Evaluate: trial traction → friction solve → imposed states
                  real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
                  fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
                                        Q_imp_plus, Q_imp_minus);

                  // 3. Rotate imposed states back to global
                  real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     Q_imp_plus_g[c] = 0.0;
                     Q_imp_minus_g[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        Q_imp_plus_g[c]  += T(c, k) * Q_imp_plus[k];
                        Q_imp_minus_g[c] += T(c, k) * Q_imp_minus[k];
                     }
                  }

                  // Godunov flux from imposed states: F = A_n^+ Q^{+,imp} + A_n^- Q^{-,imp}.
                  // Standard DG accumulation (Elem1 -= F, Elem2 += F) ensures conservation.
                  real_t F_h_total[NUM_STATE];
                  flux_.Interior(nor, Q_imp_plus_g, Q_imp_minus_g, F_h_total);

                  // Standard accumulation: Elem1 -= F, Elem2 += F
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
                  // Fault face but no DOFData mapping — fall back to welded
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

            // R-001 Part C: Check if shared face is a fault face
            bool is_fault = false;
            if (sf < static_cast<int>(shared_face_bdr_attr_.size()))
            {
               is_fault = (shared_face_bdr_attr_[sf] == bc_.fault_attr)
                        && (bc_.fault_attr > 0);
            }

            real_t F_h[NUM_STATE];

            if (is_fault && fault_flux_ && fault_dof_data_)
            {
               // Shared fault face: same pipeline as local fault faces
               auto it = shared_fault_dof_offset_.find(sf);
               int dof_idx = (it != shared_fault_dof_offset_.end())
                           ? it->second + q : -1;

               if (dof_idx >= 0 &&
                   dof_idx < static_cast<int>(fault_dof_data_->size()))
               {
                  DOFData &fdata = (*fault_dof_data_)[dof_idx];

                  real_t t1[3], t2[3];
                  GodunovFlux::BuildFrame(nor, t1, t2);
                  DenseMatrix T(NUM_STATE), Tinv(NUM_STATE);
                  GodunovFlux::BuildRotation(nor, t1, t2, T);
                  GodunovFlux::BuildRotationInverse(nor, t1, t2, Tinv);

                  real_t Q_plus_local[NUM_STATE], Q_minus_local[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     Q_plus_local[c] = 0.0;
                     Q_minus_local[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        Q_plus_local[c]  += Tinv(c, k) * Q_self[k];
                        Q_minus_local[c] += Tinv(c, k) * Q_nbr[k];
                     }
                  }

                  real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
                  // R-001 fix: peer with lower rank ID is the canonical "+" owner.
                  // R-102 fix: hard-abort on unresolved peer rank instead of
                  // silently falling back to the pre-R-001 buggy path.  Both
                  // ranks sharing a fault face would fail resolution together
                  // and both would then skip the swap, re-introducing the H3
                  // DOFData-drift bug with no warning.
                  MFEM_VERIFY(sf < static_cast<int>(shared_face_peer_.size()) &&
                              shared_face_peer_[sf].resolved,
                              "shared_face_peer_[" << sf << "] was not resolved "
                              "in the WaveOperator ctor.  Cannot canonicalise "
                              "the (+,-) side of a shared fault face — check "
                              "that face_nbr_elements_offset is populated "
                              "(pmesh.ExchangeFaceNbrData has run) and that "
                              "GetFaceNbrRank returned a valid rank.");
                  const int peer_rank = shared_face_peer_[sf].peer_rank;
                  MFEM_VERIFY(peer_rank != my_rank_,
                              "shared face " << sf << " peer_rank == my_rank_ "
                              "— ParMesh invariant violated");
                  const bool owner = (my_rank_ < peer_rank);
                  if (owner)
                  {
                     fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
                                           Q_imp_plus, Q_imp_minus);
                  }
                  else
                  {
                     // Non-owner: swap so Evaluate sees the same (Q+, Q-) as
                     // the owner.  Swap output buffers so that, in caller
                     // vocabulary, Q_imp_plus still means the e1/self-side
                     // imposed state and Q_imp_minus still means the
                     // e2/nbr-side imposed state (matches the downstream
                     // flux_.Interior convention).
                     fault_flux_->Evaluate(fdata, Q_minus_local, Q_plus_local,
                                           Q_imp_minus, Q_imp_plus);
                  }

                  real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     Q_imp_plus_g[c] = 0.0;
                     Q_imp_minus_g[c] = 0.0;
                     for (int k = 0; k < NUM_STATE; k++)
                     {
                        Q_imp_plus_g[c]  += T(c, k) * Q_imp_plus[k];
                        Q_imp_minus_g[c] += T(c, k) * Q_imp_minus[k];
                     }
                  }

                  flux_.Interior(nor, Q_imp_plus_g, Q_imp_minus_g, F_h);
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

            // Accumulate into local element only (no ghost writes)
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
      // writes.  One record = 3 centroid coords + 8 mutable fields = 11
      // doubles.  Each rank contributes one record per shared fault QP
      // it owns; the same physical QP appears on both ranks that share
      // the fault face.  Any field whose values disagree between the
      // two owners is flagged as an R-101 consistency failure.
      constexpr int REC = 11;
      constexpr int FIELD_BASE = 3;
      constexpr int NUM_FIELDS = 8;
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

      double max_diff = 0.0;
      int    max_diff_field = -1;
      int    max_diff_entry = -1;
      int    n_pairs = 0;
      int    n_unpaired = 0;

      auto same_centroid = [&](int a, int b)
      {
         // R-404 fix: sub-ULP tolerance (1 ULP * max(|va|, |vb|)).  The
         // exact compare imported by R-304 assumed `ftr->Face->Transform(ip,
         // phys)` produces bit-identical output on both ranks that share a
         // face — which holds for the current MFEM build but is not
         // documented as a contract.  The scale-relative ULP floor keeps
         // strict weak ordering (threshold bounded by the larger magnitude,
         // never an absolute constant) while accommodating a future
         // compiler-reassociation round-off that would otherwise flip
         // matched pairs into R-305 unpaired-aborts.
         for (int k = 0; k < 3; k++)
         {
            double va = all_data[a*REC + k], vb = all_data[b*REC + k];
            double scale = std::max(std::abs(va), std::abs(vb));
            double ulp = scale * std::numeric_limits<double>::epsilon();
            if (std::abs(va - vb) > ulp) { return false; }
         }
         return true;
      };

      int i = 0;
      while (i < n_entries)
      {
         int j = i + 1;
         while (j < n_entries && same_centroid(idx[i], idx[j])) { j++; }
         int group_size = j - i;
         if (group_size == 2)
         {
            int a = idx[i], b = idx[i+1];
            for (int k = FIELD_BASE; k < REC; k++)
            {
               double diff = std::abs(all_data[a*REC + k] - all_data[b*REC + k]);
               if (diff > max_diff)
               {
                  max_diff = diff;
                  max_diff_field = k;
                  max_diff_entry = a;
               }
            }
            n_pairs++;
         }
         else if (group_size != 1)
         {
            // A shared face is owned by exactly 2 ranks; anything else means
            // the centroid-match collapsed unrelated faces, which would
            // already be a logic error in this diagnostic.
            n_unpaired += group_size;
         }
         else
         {
            n_unpaired++;
         }
         i = j;
      }

      // Only rank 0 prints a summary; all ranks cooperate on the abort.
      // A shared fault QP must appear on exactly 2 ranks.  n_unpaired > 0
      // indicates a mesh-partitioning pathology (or a peer-rank resolution
      // bug that MFEM_VERIFY at the call site of Evaluate did not catch);
      // treat it as a hard diagnostic failure.
      int fail_local = (max_diff > tol || n_unpaired > 0) ? 1 : 0;
      int fail_global = 0;
      MPI_Allreduce(&fail_local, &fail_global, 1, MPI_INT, MPI_MAX, comm);

      if (fail_global)
      {
         if (n_unpaired > 0 && max_diff <= tol)
         {
            MFEM_ABORT("R-101 shared-fault DOFData: " << n_unpaired
                       << " unpaired entries (every shared QP should have "
                       "exactly 2 ranks).  Likely mesh-partitioning "
                       "pathology or peer-rank resolution failure.");
         }
         double cx = all_data[max_diff_entry*REC + 0];
         double cy = all_data[max_diff_entry*REC + 1];
         double cz = all_data[max_diff_entry*REC + 2];
         const int field_off = max_diff_field - FIELD_BASE;
         const char *field_name =
            (field_off >= 0 && field_off < NUM_FIELDS)
               ? FIELD_NAMES[field_off] : "<unknown>";
         MFEM_ABORT("R-101 shared-fault DOFData consistency FAILED.  "
                    "Field '" << field_name << "' at centroid ("
                    << cx << ", " << cy << ", " << cz
                    << ") differs by " << max_diff
                    << " across the two ranks sharing the face (tol="
                    << tol << ").  R-001's (+,-) canonicalisation is not "
                    "sufficient under the current MFEM face-normal "
                    "convention; the fix must be extended (e.g. by having "
                    "the owner rank broadcast its DOFData to the non-owner "
                    "after Evaluate).");
      }

      if (my_rank_ == 0)
      {
         std::cout << "  [R-101 check] shared-fault DOFData consistency OK: "
                   << n_pairs << " pairs matched, " << n_unpaired
                   << " unpaired entries, max_diff=" << max_diff
                   << " (tol=" << tol << ")\n";
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
