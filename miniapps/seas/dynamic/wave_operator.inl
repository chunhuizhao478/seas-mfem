// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// WaveOperator template implementation (included from wave_operator.hpp).

#include <cmath>
#include <limits>

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
      MPI_Allreduce(&h_min_, &global_h_min, 1, MPI_DOUBLE, MPI_MIN,
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

         // R-001 fix: Identify shared faces that are fault faces.
         // A shared face is a fault face if it corresponds to a boundary
         // element with fault_attr on either rank. We check the local
         // boundary elements and match by face geometry (vertex keys).
         int n_shared = pmesh.GetNSharedFaces();
         shared_face_bdr_attr_.assign(n_shared, 0);

         // Build a set of local fault face vertex keys
         std::set<std::array<int, 4>> local_fault_keys;
         for (int b = 0; b < mesh_.GetNBE(); b++)
         {
            if (mesh_.GetBdrAttribute(b) != bc_.fault_attr) { continue; }
            Array<int> verts;
            mesh_.GetBdrElementVertices(b, verts);
            std::array<int, 4> key = {0, 0, 0, 0};
            for (int v = 0; v < std::min(verts.Size(), 4); v++)
            {
               key[v] = verts[v];
            }
            std::sort(key.begin(), key.begin() + verts.Size());
            local_fault_keys.insert(key);
         }

         // R-002 fix: Detect shared fault faces using LOCAL vertex-key
         // matching only (no MPI_Allreduce — n_shared differs per rank).
         // Both ranks sharing a fault face independently detect it via
         // their respective boundary elements.
         //
         // We also record every shared face's mesh-face index in
         // shared_mesh_face_set_ so ComputeFaceFluxRHS can tell a real
         // domain-boundary face (Elem2No==-1 and NOT in the set) from
         // a partition-seam face (Elem2No==-1 BUT in the set) — the
         // latter are handled by ComputeSharedFaceFluxRHS and must
         // NOT be dispatched by the local-face loop.
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
            std::array<int, 4> key = {0, 0, 0, 0};
            for (int v = 0; v < std::min(verts.Size(), 4); v++)
            {
               key[v] = verts[v];
            }
            std::sort(key.begin(), key.begin() + verts.Size());
            shared_face_bdr_attr_[sf] = local_fault_keys.count(key)
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
         nbr_data[c] = q_gf.FaceNbrData();
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
                  fault_flux_->Evaluate(fdata, Q_plus_local, Q_minus_local,
                                        Q_imp_plus, Q_imp_minus);

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
