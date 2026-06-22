// elasticity_operator_traction.inl
// Original range: Lines 4899-6448 of elasticity_operator.hpp (pre-decomposition)
// This file is included from elasticity_operator.hpp — do not include directly.

// ============================================================================
// Out-of-class implementation of remaining methods
// ============================================================================

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GetFaultDepths(Vector &depths) const
{
   if (!fault_depths_computed_)
   {
      fault_depths_.SetSize(num_fault_dofs_);
      fault_depths_ = 0.0;

      // Nodal rule for per-DOF coordinate evaluation.
      // At nbf=1 (BR2 / order-0 face space): single centroid point.
      // At nbf>1 (IP nodal fault space): face nodes on the reference triangle.
      const IntegrationRule &nir = face_quad_->GetNodalRule();
      int nbf = nbf_per_face_;

      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         int face = fault_interior_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         for (int kk = 0; kk < nbf; kk++)
         {
            const IntegrationPoint &nip = nir.IntPoint(kk);
            FTr->SetAllIntPoints(&nip);
            Vector coords(3);
            FTr->Elem1->Transform(FTr->GetElement1IntPoint(), coords);

            // Depth: -Z (Z is negative downward in Tandem, depth is positive)
            fault_depths_(i * nbf + kk) = -coords(2);
         }
      }

      // Shared fault faces (parallel only)
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            int face_idx = fault_interior_faces_.Size() + i;
            for (int kk = 0; kk < nbf; kk++)
            {
               const IntegrationPoint &nip = nir.IntPoint(kk);
               FTr->SetAllIntPoints(&nip);
               Vector coords(3);
               FTr->Elem1->Transform(FTr->GetElement1IntPoint(), coords);

               fault_depths_(face_idx * nbf + kk) = -coords(2);
            }
         }
#endif
      }

      fault_depths_computed_ = true;
   }

   depths = fault_depths_;
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GetFaultCoords2D(
   Vector &coords_x2, Vector &coords_x3) const
{
   if (!fault_coords_computed_)
   {
      fault_x2_.SetSize(num_fault_dofs_);
      fault_x3_.SetSize(num_fault_dofs_);
      fault_x2_ = 0.0;
      fault_x3_ = 0.0;

      // Nodal rule for per-DOF coordinate evaluation
      const IntegrationRule &nir = face_quad_->GetNodalRule();
      int nbf = nbf_per_face_;

      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         int face = fault_interior_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         for (int kk = 0; kk < nbf; kk++)
         {
            const IntegrationPoint &nip = nir.IntPoint(kk);

            // Use Elem1->Transform for robust coordinate evaluation
            // (Face->Transform can disagree on BP5/tet path — v59 Section 6)
            FTr->SetAllIntPoints(&nip);
            Vector coords(3);
            FTr->Elem1->Transform(FTr->GetElement1IntPoint(), coords);

            // Tandem: X=along-strike=coords(0), depth=-Z=-coords(2)
            fault_x2_(i * nbf + kk) = coords(0);
            fault_x3_(i * nbf + kk) = -coords(2);
         }
      }

      // Shared fault faces (parallel only)
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            int face_idx = fault_interior_faces_.Size() + i;
            for (int kk = 0; kk < nbf; kk++)
            {
               const IntegrationPoint &nip = nir.IntPoint(kk);
               FTr->SetAllIntPoints(&nip);
               Vector coords(3);
               FTr->Elem1->Transform(FTr->GetElement1IntPoint(), coords);

               fault_x2_(face_idx * nbf + kk) = coords(0);
               fault_x3_(face_idx * nbf + kk) = -coords(2);
            }
         }
#endif
      }

      fault_coords_computed_ = true;
   }

   coords_x2 = fault_x2_;
   coords_x3 = fault_x3_;
}

namespace
{

// R-104: per-face helpers shared between the interior- and shared-face
// loops of GetFaultDOFCoords3D / GetFaultDOFBasis.  Keeping the per-DOF
// computation in one place avoids the maintenance hazard of changing
// it in one loop and forgetting the other.

inline void WriteFaceDOFCoords3D_(
   mfem::FaceElementTransformations *FTr,
   int face_dof_offset, int nbf,
   const mfem::IntegrationRule &nir,
   mfem::Vector &dof_coords_3d)
{
   for (int kk = 0; kk < nbf; kk++)
   {
      const mfem::IntegrationPoint &nip = nir.IntPoint(kk);
      FTr->SetAllIntPoints(&nip);
      mfem::Vector coords(3);
      FTr->Elem1->Transform(FTr->GetElement1IntPoint(), coords);
      const int dof_idx = face_dof_offset + kk;
      dof_coords_3d(3 * dof_idx + 0) = coords(0);
      dof_coords_3d(3 * dof_idx + 1) = coords(1);
      dof_coords_3d(3 * dof_idx + 2) = coords(2);
   }
}

inline void WriteFaceDOFBasis_(
   mfem::FaceElementTransformations *FTr,
   int face_dof_offset, int nbf, int dim,
   const mfem::IntegrationRule &nir,
   const mfem::Vector &ref_normal,
   const mfem::Vector &up,
   const char *face_kind,    // "face" or "shared face" — for error msg
   int face_id_for_error,
   mfem::DenseMatrix &dof_basis)
{
   for (int kk = 0; kk < nbf; kk++)
   {
      const mfem::IntegrationPoint &nip = nir.IntPoint(kk);
      FTr->Face->SetIntPoint(&nip);

      mfem::Vector n_raw(dim);
      mfem::CalcOrtho(FTr->Face->Jacobian(), n_raw);

      mfem::real_t normal[3], t1[3], t2[3];
      bool sign_flipped;
      mfem::real_t nl;
      mfem::seas::FaultBasis::ComputeOrientedFrame(
         n_raw, dim, ref_normal, up,
         normal, t1, t2, sign_flipped, nl);
      MFEM_VERIFY(nl > 0.0,
                  "GetFaultDOFBasis: zero-length " << face_kind
                  << " normal at " << face_kind << " "
                  << face_id_for_error << " DOF " << kk);

      const int dof_idx = face_dof_offset + kk;
      for (int d = 0; d < 3; d++)
      {
         dof_basis(d,     dof_idx) = normal[d];
         dof_basis(3 + d, dof_idx) = t1[d];
         dof_basis(6 + d, dof_idx) = t2[d];
      }
   }
}

} // anonymous namespace

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GetFaultDOFCoords3D(
   Vector &dof_coords_3d) const
{
   // Phase 6.A — per-DOF (x_i, y_i, z_i) global UTM coordinates for the
   // SAFS sidecar lookup. Mirrors GetFaultCoords2D's iteration pattern;
   // writes all three components without the GetFaultCoords2D depth flip.
   dof_coords_3d.SetSize(3 * num_fault_dofs_);
   dof_coords_3d = 0.0;

   if (num_fault_dofs_ == 0) { return; }

   const IntegrationRule &nir = face_quad_->GetNodalRule();
   const int nbf = nbf_per_face_;

   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      int face = fault_interior_faces_[i];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);
      if (FTr == nullptr) { continue; }
      WriteFaceDOFCoords3D_(FTr, i * nbf, nbf, nir, dof_coords_3d);
   }

   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         int sf = fault_shared_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(sf);
         if (FTr == nullptr) { continue; }
         const int face_idx = fault_interior_faces_.Size() + i;
         WriteFaceDOFCoords3D_(FTr, face_idx * nbf, nbf, nir,
                               dof_coords_3d);
      }
#endif
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GetFaultDOFBasis(
   DenseMatrix &dof_basis) const
{
   // Phase 6.A — per-DOF (n, t1, t2) basis via FaultBasis::ComputeOrientedFrame
   // evaluated at the same nodal reference points used for the per-DOF
   // coordinate accessors. Sign convention identical to fault_basis_
   // (sign flip baked in, Tandem convention).
   dof_basis.SetSize(9, num_fault_dofs_);
   dof_basis = 0.0;

   if (num_fault_dofs_ == 0) { return; }

   const int dim = mesh_.Dimension();
   MFEM_VERIFY(dim == 3,
               "GetFaultDOFBasis: 3-D mesh required for per-DOF "
               "(n, t1, t2) accessor (got dim = " << dim << ")");

   const IntegrationRule &nir = face_quad_->GetNodalRule();
   const int nbf = nbf_per_face_;

   // Use the same up vector as the fault_basis_ initialiser
   // (elasticity_operator_setup.inl:735-737) so per-DOF and per-face
   // bases agree at constant-Jacobian faces.
   Vector up(3);
   up = 0.0;
   up(2) = 1.0;

   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      int face = fault_interior_faces_[i];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);
      if (FTr == nullptr) { continue; }
      WriteFaceDOFBasis_(FTr, i * nbf, nbf, dim, nir, ref_normal_, up,
                         "face", face, dof_basis);
   }

   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         int sf = fault_shared_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(sf);
         if (FTr == nullptr) { continue; }
         const int face_idx = fault_interior_faces_.Size() + i;
         WriteFaceDOFBasis_(FTr, face_idx * nbf, nbf, dim, nir,
                            ref_normal_, up,
                            "shared face", sf, dof_basis);
      }
#endif
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::RestrictToOwnedFault(
   const Vector &local_data, Vector &owned_data, int comps_per_dof) const
{
   MFEM_VERIFY(comps_per_dof > 0, "comps_per_dof must be positive");
   MFEM_VERIFY(local_data.Size() == comps_per_dof * num_fault_dofs_,
               "Local fault vector size mismatch: got " << local_data.Size()
               << ", expected " << comps_per_dof * num_fault_dofs_);

   owned_data.SetSize(comps_per_dof * num_owned_fault_dofs_);
   for (int owned_dof = 0; owned_dof < num_owned_fault_dofs_; owned_dof++)
   {
      const int local_dof = owned_fault_dof_to_local_dof_[owned_dof];
      for (int c = 0; c < comps_per_dof; c++)
      {
         owned_data(comps_per_dof * owned_dof + c) =
            local_data(comps_per_dof * local_dof + c);
      }
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ExpandOwnedToLocalFault(
   const Vector &owned_data, Vector &local_data, int comps_per_dof) const
{
   MFEM_VERIFY(comps_per_dof > 0, "comps_per_dof must be positive");
   MFEM_VERIFY(owned_data.Size() == comps_per_dof * num_owned_fault_dofs_,
               "Owned fault vector size mismatch: got " << owned_data.Size()
               << ", expected " << comps_per_dof * num_owned_fault_dofs_);

   local_data.SetSize(comps_per_dof * num_fault_dofs_);
   local_data = 0.0;

   for (int owned_dof = 0; owned_dof < num_owned_fault_dofs_; owned_dof++)
   {
      const int local_dof = owned_fault_dof_to_local_dof_[owned_dof];
      for (int c = 0; c < comps_per_dof; c++)
      {
         local_data(comps_per_dof * local_dof + c) =
            owned_data(comps_per_dof * owned_dof + c);
      }
   }

   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      if (shared_fault_comm_blocks_.empty()) { return; }

      // Use FaultScatter to handle all point-to-point communication.
      FaultScatter scatter(shared_fault_comm_blocks_, mesh_.GetComm());
      scatter.BeginScatter(owned_data, comps_per_dof, nbf_per_face_);
      scatter.WaitScatter();

      // Unpack received data with canonical→local DOF permutation.
      const int block_size = comps_per_dof * nbf_per_face_;
      for (int bi = 0; bi < scatter.NumBlocks(); bi++)
      {
         const auto &block = scatter.GetBlock(bi);
         const Vector &recv = scatter.GetRecvBuffer(bi);
         if (recv.Size() == 0) { continue; }

         for (int j = 0; j < static_cast<int>(block.recv_local_faces.size()); j++)
         {
            const int local_face = block.recv_local_faces[j];
            // Received data is in canonical order (sorted by global vertex ID).
            // Apply this rank's canonical→local permutation so the values
            // land at the correct MFEM-local DOF positions.
            for (int canonical_kk = 0; canonical_kk < nbf_per_face_; canonical_kk++)
            {
               const int mfem_kk =
                  canonical_to_local_perm_[local_face * nbf_per_face_ + canonical_kk];
               const int local_dof = local_face * nbf_per_face_ + mfem_kk;
               for (int c = 0; c < comps_per_dof; c++)
               {
                  local_data(comps_per_dof * local_dof + c) =
                     recv(j * block_size + canonical_kk * comps_per_dof + c);
               }
            }
         }
      }
#endif
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GetFaultDOFToElem(
   Array<int> &dof_to_elem) const
{
   // LOCAL-order Elem1No (same face walk as GetFaultDOFCoords3D), then
   // restricted to the OWNED set via owned_fault_dof_to_local_dof_.
   Array<int> local(num_fault_dofs_);
   local = -1;
   const int nbf = nbf_per_face_;
   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(fault_interior_faces_[i]);
      if (FTr == nullptr) { continue; }
      for (int kk = 0; kk < nbf; kk++) { local[i * nbf + kk] = FTr->Elem1No; }
   }
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(fault_shared_faces_[i]);
         if (FTr == nullptr) { continue; }
         const int face_idx = fault_interior_faces_.Size() + i;
         for (int kk = 0; kk < nbf; kk++)
         {
            local[face_idx * nbf + kk] = FTr->Elem1No;
         }
      }
#endif
   }
   dof_to_elem.SetSize(num_owned_fault_dofs_);
   for (int od = 0; od < num_owned_fault_dofs_; od++)
   {
      dof_to_elem[od] = local[owned_fault_dof_to_local_dof_[od]];
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GetFaultDOFToAttr(
   Array<int> &dof_to_attr) const
{
   // Every owned fault DOF lies on the fault interface, so its attribute is
   // the fault boundary attribute (matching spatial_dyn_driver's dof_to_attr).
   dof_to_attr.SetSize(num_owned_fault_dofs_);
   dof_to_attr = bdr_config_.fault_attr;
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GetFaultDOFIntegrationPoints(
   std::vector<IntegrationPoint> &dof_ips) const
{
   // LOCAL-order reference IPs in Elem1's frame (same nodal rule + face walk
   // as GetFaultDOFCoords3D), then restricted to the OWNED set.
   std::vector<IntegrationPoint> local(num_fault_dofs_);
   if (num_fault_dofs_ == 0) { dof_ips.clear(); return; }

   const IntegrationRule &nir = face_quad_->GetNodalRule();
   const int nbf = nbf_per_face_;
   for (int i = 0; i < fault_interior_faces_.Size(); i++)
   {
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(fault_interior_faces_[i]);
      if (FTr == nullptr) { continue; }
      for (int kk = 0; kk < nbf; kk++)
      {
         const IntegrationPoint &nip = nir.IntPoint(kk);
         FTr->SetAllIntPoints(&nip);
         local[i * nbf + kk] = FTr->GetElement1IntPoint();
      }
   }
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(fault_shared_faces_[i]);
         if (FTr == nullptr) { continue; }
         const int face_idx = fault_interior_faces_.Size() + i;
         for (int kk = 0; kk < nbf; kk++)
         {
            const IntegrationPoint &nip = nir.IntPoint(kk);
            FTr->SetAllIntPoints(&nip);
            local[face_idx * nbf + kk] = FTr->GetElement1IntPoint();
         }
      }
#endif
   }
   dof_ips.resize(num_owned_fault_dofs_);
   for (int od = 0; od < num_owned_fault_dofs_; od++)
   {
      dof_ips[od] = local[owned_fault_dof_to_local_dof_[od]];
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::Solve(
   real_t time, const Vector &slip_bc, GridFuncType &displacement)
{
   if (!stiffness_assembled_)
   {
      AssembleStiffness();
   }

   // Build RHS
   LinFormType b(fes_.get());
   b.Assemble();
   Vector &rhs = b;

   const bool debug_first_step = (method_ == DGMethod::IP) &&
                                 DebugEnabledForTime(time);
   debug_time_ = time;

   // Production assembly path: same on ALL ranks (MPI-safe).
   // Debug dumps happen AFTER the solve, gated by rank — no separate
   // assembly path that could cause MPI divergence.
   if (method_ == DGMethod::IP)
   {
      AssembleSlipContributionIP(rhs, slip_bc);
      AssembleSlipContributionIPShared(rhs, slip_bc,
                                       fault_interior_faces_.Size());
   }
   else
   {
      AssembleSlipContributionBR2(rhs, slip_bc);
      AssembleSlipContributionBR2Shared(rhs, slip_bc,
                                        fault_interior_faces_.Size());
   }

   // Snapshot slip-only RHS before adding Dirichlet.
   // Taken on ALL ranks when debug is enabled (needed for MPI-collective norms).
   Vector rhs_slip_snapshot;
   const bool debug_norm = first_step_debug_.enabled && (time > 0.0);
   if (debug_first_step || debug_norm)
   {
      rhs_slip_snapshot = rhs;  // copy before Dirichlet is added
   }

   // Add Dirichlet loading
   AssembleDirichletLoading(rhs, time);

   X_ = 0.0;
   B_ = rhs;

   solver_->Mult(B_, X_);

   // Post-solve residual check: ||K*x - b|| / ||b|| (global norms).  Phase 4:
   // mandatory on the AMG paths (cfg.solver.residual_check defaults true) — a
   // silently non-converged Krylov solve corrupts traction -> friction ->
   // blowup.  Threshold tied to the Krylov tolerance (10*ksp_rtol); warn LOUDLY
   // with the iteration count (iterative solvers only; -1 for direct).
   if (check_residual_)
   {
      Vector R_(B_.Size());
      const real_t res_threshold = 10.0 * ksp_rtol_;
      auto *itsolver = dynamic_cast<IterativeSolver *>(solver_.get());
      const int  ksp_iters = itsolver ? itsolver->GetNumIterations() : -1;
      const bool ksp_conv  = itsolver ? itsolver->GetConverged() : true;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         cached_Ah_.As<HypreParMatrix>()->Mult(X_, R_);
         R_ -= B_;
         // Use global norms (MPI AllReduce) for meaningful relative residual
         real_t local_res2 = R_ * R_;
         real_t local_rhs2 = B_ * B_;
         real_t global_res2 = 0.0, global_rhs2 = 0.0;
         MPI_Allreduce(&local_res2, &global_res2, 1, MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
         MPI_Allreduce(&local_rhs2, &global_rhs2, 1, MPI_DOUBLE, MPI_SUM, mesh_.GetComm());
         real_t global_res = std::sqrt(global_res2);
         real_t global_rhs = std::sqrt(global_rhs2);
         real_t rel_res = (global_rhs > 0.0) ? global_res / global_rhs : global_res;
         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         if (rel_res > res_threshold && rank == 0)
         {
            mfem::err << "RESIDUAL WARNING: global ||K*x-b||/||b|| = " << rel_res
                      << " > " << res_threshold << " (10*ksp_rtol); Krylov iters="
                      << ksp_iters
                      << (ksp_conv ? "" : " (DID NOT CONVERGE)")
                      << " — non-converged solve corrupts traction.\n";
         }
#endif
      }
      else
      {
         cached_a_->SpMat().Mult(X_, R_);
         R_ -= B_;
         real_t res_norm = R_.Norml2();
         real_t rhs_norm = B_.Norml2();
         real_t rel_res = (rhs_norm > 0.0) ? res_norm / rhs_norm : res_norm;
         if (rel_res > res_threshold)
         {
            mfem::err << "RESIDUAL WARNING: ||K*x-b||/||b|| = " << rel_res
                      << " > " << res_threshold << " (10*ksp_rtol); Krylov iters="
                      << ksp_iters
                      << (ksp_conv ? "" : " (DID NOT CONVERGE)") << "\n";
         }
      }
   }

   // Check for NaN/Inf in solution
   {
      real_t u_max = X_.Normlinf();
      MFEM_VERIFY(std::isfinite(u_max),
         "Domain solve produced NaN/Inf in displacement (||u||_inf = "
         << u_max << ")");
   }

   // Check convergence and log solver info when RHS is large
   {
      auto *iter_solver = dynamic_cast<IterativeSolver*>(solver_.get());
      if (iter_solver)
      {
         // A stalled-but-finite iterative solve silently corrupts the elastic
         // displacement -> traction -> friction -> state, so it must fail fast
         // rather than poison the time loop.  (MUMPS / direct solvers are not
         // IterativeSolver, so this guard does not apply to them.)
         MFEM_VERIFY(iter_solver->GetConverged(),
            "Iterative domain solve did NOT converge after "
            << iter_solver->GetNumIterations() << " iterations (final norm = "
            << iter_solver->GetFinalNorm() << ").  This corrupts traction and "
            "the friction coupling.  Improve the AMG preconditioner, raise "
            "max_iter, or switch to a direct solver (solver_type=\"mumps\").");
      }
   }

   displacement = X_;

   // Debug dumps: rank-local file I/O only, no MPI collectives.
   // Placed AFTER displacement = X_ so the GridFunction has the current solution.
   if (debug_first_step)
   {
      // Compute Dirichlet-only RHS by subtraction (no extra assembly)
      Vector rhs_dir(rhs.Size());
      subtract(rhs, rhs_slip_snapshot, rhs_dir);

      DebugDumpElementVector("rhs_slip", rhs_slip_snapshot);
      DebugDumpElementVector("rhs_dirichlet", rhs_dir);
      DebugDumpElementVector("rhs_total", rhs);
      DebugDumpElementVector("u", X_);

      // Dump fault jumps and traction using the solved displacement.
      // Interior faces only — no ExchangeFaceNbrData (MPI-safe).
      DebugDumpFaultJumpsLocal(displacement, slip_bc);
      DebugDumpFaultTractionLocal(displacement, slip_bc);
      DebugDumpKContributions();
      first_step_debug_done_ = true;
   }

   // Global norm comparison — OUTSIDE rank-gated block.
   // All ranks must participate in MPI_Reduce (no deadlock).
   // Condition uses only rank-independent flags (enabled + time).
   // Prints at EVERY Mult() call during first step to trace RK45 stages.
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      static int norm_count = 0;
      if (norm_count < 10 && first_step_debug_.enabled && time > 0.0)
      {
         norm_count++;
         int myrank;
         MPI_Comm_rank(mesh_.GetComm(), &myrank);
         if (myrank == 0)
         {
            mfem::out << "  [NORM] Stage " << norm_count
                      << " at time = " << std::setprecision(17) << time << "\n";
         }

         auto print_global_norm = [&](const char *label, const Vector &v) {
            double local_n1 = 0.0, local_n2sq = 0.0, local_ninf = 0.0;
            for (int i = 0; i < v.Size(); i++)
            {
               double a = std::abs(v(i));
               local_n1 += a;
               local_n2sq += a * a;
               if (a > local_ninf) { local_ninf = a; }
            }
            double g1, g2sq, ginf;
            MPI_Reduce(&local_n1, &g1, 1, MPI_DOUBLE, MPI_SUM, 0, mesh_.GetComm());
            MPI_Reduce(&local_n2sq, &g2sq, 1, MPI_DOUBLE, MPI_SUM, 0, mesh_.GetComm());
            MPI_Reduce(&local_ninf, &ginf, 1, MPI_DOUBLE, MPI_MAX, 0, mesh_.GetComm());
            if (myrank == 0)
            {
               mfem::out << std::setprecision(15);
               mfem::out << "  [NORM] ||" << label << "||_1   = " << g1 << "\n";
               mfem::out << "  [NORM] ||" << label << "||_2   = "
                         << std::sqrt(g2sq) << "\n";
               mfem::out << "  [NORM] ||" << label << "||_inf = " << ginf << "\n";
            }
         };

         print_global_norm("b_slip", rhs_slip_snapshot);
         print_global_norm("b_total", rhs);
         print_global_norm("u", X_);
      }
#endif
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTraction(
   const GridFuncType &displacement,
   const Vector &slip_bc,
   Vector &traction,
   Vector *normal_traction)
{
   ComputeTractionImpl(displacement, slip_bc, traction, normal_traction,
                       nullptr, nullptr, nullptr);
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::AssembleSlipOnlyRHS(
   Vector &rhs, const Vector &slip_bc) const
{
   rhs.SetSize(fes_->GetVSize());
   rhs = 0.0;

   if (method_ == DGMethod::IP)
   {
      AssembleSlipContributionIP(rhs, slip_bc);
      AssembleSlipContributionIPShared(rhs, slip_bc,
                                       fault_interior_faces_.Size());
   }
   else
   {
      AssembleSlipContributionBR2(rhs, slip_bc);
      AssembleSlipContributionBR2Shared(rhs, slip_bc,
                                        fault_interior_faces_.Size());
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::AssembleDirichletOnlyRHS(
   Vector &rhs, real_t time) const
{
   rhs.SetSize(fes_->GetVSize());
   rhs = 0.0;
   AssembleDirichletLoading(rhs, time);
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTractionComponents(
   const GridFuncType &displacement,
   const Vector &slip_bc,
   Vector &traction,
   Vector &traction_stress,
   Vector &traction_correction,
   Vector *normal_traction)
{
   ComputeTractionImpl(displacement, slip_bc, traction, normal_traction,
                       &traction_stress, &traction_correction, nullptr);
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTractionDiagnostics(
   const GridFuncType &displacement,
   const Vector &slip_bc,
   Vector &traction,
   Vector &traction_stress,
   Vector &traction_correction,
   Vector &jump_residual,
   Vector *normal_traction,
   Vector *normal_stress,
   Vector *normal_correction)
{
   ComputeTractionImpl(displacement, slip_bc, traction, normal_traction,
                       &traction_stress, &traction_correction, &jump_residual,
                       normal_stress, normal_correction);
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTractionImpl(
   const GridFuncType &displacement,
   const Vector &slip_bc,
   Vector &traction,
   Vector *normal_traction,
   Vector *traction_stress_out,
   Vector *traction_correction_out,
   Vector *jump_residual_out,
   Vector *normal_stress_out,
   Vector *normal_correction_out)
{
   int dim = 3;
   traction.SetSize(2 * num_fault_dofs_);
   traction = 0.0;

    if (traction_stress_out)
    {
       traction_stress_out->SetSize(2 * num_fault_dofs_);
       *traction_stress_out = 0.0;
    }
    if (traction_correction_out)
    {
       traction_correction_out->SetSize(2 * num_fault_dofs_);
       *traction_correction_out = 0.0;
    }
    if (jump_residual_out)
    {
       jump_residual_out->SetSize(2 * num_fault_dofs_);
       *jump_residual_out = 0.0;
    }
   if (normal_stress_out)
   {
      normal_stress_out->SetSize(num_fault_dofs_);
      *normal_stress_out = 0.0;
   }
   if (normal_correction_out)
   {
      normal_correction_out->SetSize(num_fault_dofs_);
      *normal_correction_out = 0.0;
   }

   // v51: Elastic normal traction for sigma_n feedback
   if (normal_traction)
   {
      normal_traction->SetSize(num_fault_dofs_);
      *normal_traction = 0.0;
   }

   for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
   {
      int face = fault_interior_faces_[fi];
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(face);
      if (FTr == nullptr) { continue; }

      // ============================================================
      // v35: Per-quadrature-point traction evaluation with face averaging
      //
      // Previous code evaluated {σ·n̂} at face centroid only and used
      // face-averaged shapes for BR2 correction. This is exact for p=1
      // but wrong for p≥2 (stress varies across face, lifted function
      // has spatial variation). Now we evaluate everything at each face
      // quadrature point and face-average the result.
      // ============================================================

      // Face centroid IP (used for diagnostic coordinate output only)
      const IntegrationPoint &ip =
         Geometries.GetCenter(FTr->GetGeometryType());

      // Get DOFs
      Array<int> vdofs1, vdofs2;
      fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
      fes_->GetElementVDofs(FTr->Elem2No, vdofs2);

      const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
      const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
      int ndof1 = fe1->GetDof();
      int ndof2 = fe2->GetDof();

      // Get displacement values
      Vector u1_all(vdofs1.Size()), u2_all(vdofs2.Size());
      displacement.GetSubVector(vdofs1, u1_all);
      displacement.GetSubVector(vdofs2, u2_all);

      // Fault basis
      const auto &basis = fault_basis_.GetBasis(fi);

      // Element Jacobian inverses (constant for linear tets/hexes)
      // Must set an integration point first so Jacobian() is valid.
      FTr->SetAllIntPoints(&ip);
      DenseMatrix Jinv1(dim), Jinv2(dim);
      CalcInverse(FTr->Elem1->Jacobian(), Jinv1);
      CalcInverse(FTr->Elem2->Jacobian(), Jinv2);

      // Face quadrature rule
      int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
      int quad_order_trac = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
      const IntegrationRule &ir_trac = IntRules.Get(
         FTr->GetGeometryType(), quad_order_trac);
      int nqp = ir_trac.GetNPoints();

      // Face-averaged accumulators (used for BR2 path and diagnostics)
      real_t T_global[3] = {0.0, 0.0, 0.0};
      real_t T_stress[3] = {0.0, 0.0, 0.0};
      real_t correction[3] = {0.0, 0.0, 0.0};
      real_t residual_avg[3] = {0.0, 0.0, 0.0};
      real_t sum_wq = 0.0;

      if (method_ == DGMethod::IP)
      {
         // Tandem-style traction recovery using combined integrator.
         // ComputeTractionAtQuadPoints (per-qp geometry, same penalty as K)
         // + ProjectTractionToFaultDOFs (nl_q-weighted L2 projection)
         // When decomposition diagnostics are requested, uses the decomposed
         // variant to also output stress, correction, and jump residual.
         int nbf = nbf_per_face_;
         bool need_decomp = traction_stress_out || traction_correction_out ||
                            jump_residual_out ||
                            normal_stress_out || normal_correction_out;

         // Build sign-corrected slip at quad points (Tandem evaluate_slip).
         // 1. Collect tangential slip components (dip, strike) per DOF
         // 2. Interpolate to quad points using face basis functions e_q
         // 3. Embed at each QP using per-QP tangent frame
         // This matches Tandem's tensor contraction:
         //   slip_q[p,q] = e_q[l,q] * fault_basis_q[p,o,q] * slip[l,n] * copy_slip[n,o]
         int qo_slip = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
         const IntegrationRule &ir_slip = IntRules.Get(
            FTr->GetGeometryType(), qo_slip);
         int nqp_slip = ir_slip.GetNPoints();
         Vector delta_u_quad_t;

         if (!basis.qp_data.empty())
         {
            // Per-QP tangent embedding (matches Tandem for any nbf)
            // Step 1: tangential slip components per DOF
            Vector slip_tang(2 * nbf);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf + kk;
               slip_tang(0 * nbf + kk) = slip_bc(2 * dof_idx);      // dip
               slip_tang(1 * nbf + kk) = slip_bc(2 * dof_idx + 1);  // strike
            }
            // Step 2: interpolate tangential components to QPs
            Vector slip_tang_q;
            face_quad_->InterpolateToQuadPoints(2, slip_tang, slip_tang_q);
            // Step 3: embed at each QP using per-QP tangent frame
            delta_u_quad_t.SetSize(dim * nqp_slip);
            for (int q = 0; q < nqp_slip; q++)
            {
               real_t sl_q[2] = {slip_tang_q(q), slip_tang_q(nqp_slip + q)};
               real_t du[3];
               fault_basis_.EmbedSlipQP(fi, q, sl_q, du);
               for (int c = 0; c < dim; c++)
                  delta_u_quad_t(c * nqp_slip + q) = du[c];
            }
         }
         else
         {
            // Fallback: centroid tangents (exact for flat faces, used by BR2)
            Vector delta_u_nodal_t(dim * nbf);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf + kk;
               real_t sl[2] = {slip_bc(2 * dof_idx), slip_bc(2 * dof_idx + 1)};
               real_t du[3];
               fault_basis_.EmbedSlip(fi, sl, du);
               for (int c = 0; c < dim; c++)
                  delta_u_nodal_t(c * nbf + kk) = du[c];
            }
            face_quad_->InterpolateToQuadPoints(dim, delta_u_nodal_t, delta_u_quad_t);
         }

         // Step 1: Traction at quad points (Tandem compute_traction)
         DGElasticityIPCombinedIntegrator trac_integ(
            lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
         Vector T_quad_new, nl_q_vec;
         Vector T_stress_quad_dec, T_corr_quad_dec, R_quad_dec;
         if (need_decomp)
         {
            trac_integ.ComputeTractionAtQuadPointsDecomposed(
               *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad_t,
               T_quad_new,
               traction_stress_out ? &T_stress_quad_dec : nullptr,
               traction_correction_out ? &T_corr_quad_dec : nullptr,
               jump_residual_out ? &R_quad_dec : nullptr,
               nullptr, &nl_q_vec);
         }
         else
         {
            trac_integ.ComputeTractionAtQuadPoints(
               *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad_t,
               T_quad_new, nullptr, &nl_q_vec);
         }
         int nqp_new = T_quad_new.Size() / dim;

         // Step 2: Project to fault DOFs (Tandem evaluate_traction)
         int quad_order_new = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
         const IntegrationRule &ir_new = IntRules.Get(
            FTr->GetGeometryType(), quad_order_new);
         const DenseMatrix &e_q = face_quad_->BasisAtQuadPoints();

         const auto *qpd = basis.qp_data.empty() ? nullptr : &basis.qp_data;

         // Unified projection: project all local components through one
         // L2 solve (matching Tandem's single-kernel approach).
         // When normal_traction is requested, project 3 components
         // (normal, dip, strike); otherwise project 2 (dip, strike).
         int ncomp_proj = normal_traction ? 3 : 2;
         real_t basis_vecs[3][3];
         if (normal_traction)
         {
            for (int d = 0; d < 3; d++)
            {
               basis_vecs[0][d] = basis.normal[d];
               basis_vecs[1][d] = basis.tangent1[d];
               basis_vecs[2][d] = basis.tangent2[d];
            }
         }
         else
         {
            for (int d = 0; d < 3; d++)
            {
               basis_vecs[0][d] = basis.tangent1[d];
               basis_vecs[1][d] = basis.tangent2[d];
               basis_vecs[2][d] = 0.0;
            }
         }

         Vector trac_local_new;
         DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
            dim, ncomp_proj, T_quad_new, nl_q_vec, ir_new, nbf, e_q,
            basis_vecs, trac_local_new, qpd);

         // Store per-DOF shear traction
         int tang_offset = normal_traction ? 1 : 0;
         for (int kk = 0; kk < nbf; kk++)
         {
            int dof_idx = fi * nbf_per_face_ + kk;
            traction(2 * dof_idx)     = trac_local_new((tang_offset + 0) * nbf + kk);
            traction(2 * dof_idx + 1) = trac_local_new((tang_offset + 1) * nbf + kk);
         }

         // Normal stress: extracted from the unified 3-component projection
         // above (Tandem-style: normal + dip + strike through one L2 solve).
         // Sign convention: positive in compression = -(T · n̂).
         if (normal_traction)
         {
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*normal_traction)(dof_idx) = -trac_local_new(0 * nbf + kk);
            }
         }

         // Project decomposition components to fault DOFs (if requested)
         if (traction_stress_out)
         {
            Vector stress_local;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, 2, T_stress_quad_dec, nl_q_vec, ir_new, nbf, e_q,
               basis_vecs + tang_offset, stress_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*traction_stress_out)(2 * dof_idx)     = stress_local(0 * nbf + kk);
               (*traction_stress_out)(2 * dof_idx + 1) = stress_local(1 * nbf + kk);
            }
         }
         if (traction_correction_out)
         {
            Vector corr_local;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, 2, T_corr_quad_dec, nl_q_vec, ir_new, nbf, e_q,
               basis_vecs + tang_offset, corr_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*traction_correction_out)(2 * dof_idx)     = corr_local(0 * nbf + kk);
               (*traction_correction_out)(2 * dof_idx + 1) = corr_local(1 * nbf + kk);
            }
         }
         if (jump_residual_out)
         {
            Vector res_local;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, 2, R_quad_dec, nl_q_vec, ir_new, nbf, e_q,
               basis_vecs + tang_offset, res_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*jump_residual_out)(2 * dof_idx)     = res_local(0 * nbf + kk);
               (*jump_residual_out)(2 * dof_idx + 1) = res_local(1 * nbf + kk);
            }
         }

         // Normal decomposition: project stress and correction onto normal
         // Same sign convention as total normal traction: -(T · n̂)
         if (normal_stress_out && normal_traction)
         {
            Vector ns_local;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, 1, T_stress_quad_dec, nl_q_vec, ir_new, nbf, e_q,
               basis_vecs, ns_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*normal_stress_out)(dof_idx) = -ns_local(kk);
            }
         }
         if (normal_correction_out && normal_traction)
         {
            Vector nc_local;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, 1, T_corr_quad_dec, nl_q_vec, ir_new, nbf, e_q,
               basis_vecs, nc_local, qpd);
            for (int kk = 0; kk < nbf; kk++)
            {
               int dof_idx = fi * nbf_per_face_ + kk;
               (*normal_correction_out)(dof_idx) = -nc_local(kk);
            }
         }

      }
      else  // BR2
      {
         // BR2 always uses nbf=1, so face-index == DOF-index
         real_t slip_local_br2[2] = {slip_bc(2 * fi), slip_bc(2 * fi + 1)};
         real_t delta_u[3];
         fault_basis_.EmbedSlip(fi, slip_local_br2, delta_u);

         if (!mass_inv_computed_) {
            PrecomputeMassInverse();
            // Re-fetch FTr: PrecomputeMassInverse calls GetElementTransformation
            // for every element, which may invalidate MFEM's cached face transforms.
            FTr = mesh_.GetInteriorFaceTransformations(face);
         }

         Geometry::Type geom = mesh_.GetElementGeometry(FTr->Elem1No);
         real_t br2_penalty = (geom == Geometry::TETRAHEDRON)
                                  ? real_t(dim + 1) : real_t(2 * dim);

         const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
         const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];

         // Precompute shapes, normals, and face integrals for BR2 lifting
         DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
         Vector nor_arr(dim * nqp), w_arr(nqp);

         DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
         face_int1 = 0.0;
         face_int2 = 0.0;

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &fip = ir_trac.IntPoint(q);
            FTr->SetAllIntPoints(&fip);
            const IntegrationPoint &eip1_q = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2_q = FTr->GetElement2IntPoint();

            Vector s1q(shapes1.GetColumn(q), ndof1);
            Vector s2q(shapes2.GetColumn(q), ndof2);
            fe1->CalcShape(eip1_q, s1q);
            fe2->CalcShape(eip2_q, s2q);

            Vector nor_q(&nor_arr[q * dim], dim);
            CalcOrtho(FTr->Jacobian(), nor_q);
            w_arr[q] = fip.weight;

            // Per-quadrature-point displacement jump for lifting
            real_t jump_q[3] = {0.0, 0.0, 0.0};
            for (int c = 0; c < dim; c++)
            {
               real_t u1q = 0.0, u2q = 0.0;
               for (int k = 0; k < ndof1; k++)
                  u1q += s1q(k) * u1_all(c * ndof1 + k);
               for (int k = 0; k < ndof2; k++)
                  u2q += s2q(k) * u2_all(c * ndof2 + k);
               jump_q[c] = (u1q - u2q) - delta_u[c];
            }

            for (int u = 0; u < dim; u++)
               for (int s = 0; s < dim; s++)
               {
                  for (int m = 0; m < ndof1; m++)
                     face_int1(u * dim + s, m) +=
                        w_arr[q] * s1q(m) * jump_q[u] * nor_q(s);
                  for (int m = 0; m < ndof2; m++)
                     face_int2(u * dim + s, m) +=
                        w_arr[q] * s2q(m) * jump_q[u] * nor_q(s);
               }
         }

         // Lift: f_lifted = 0.5 * Minv * face_int
         DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
         MultABt(face_int1, Minv1, f_lifted1);
         f_lifted1 *= 0.5;
         MultABt(face_int2, Minv2, f_lifted2);
         f_lifted2 *= 0.5;

         // Pass 2: evaluate traction at each quadrature point
         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &fip = ir_trac.IntPoint(q);
            FTr->SetAllIntPoints(&fip);
            const IntegrationPoint &eip1_q = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2_q = FTr->GetElement2IntPoint();
            real_t wq = w_arr[q];

            // {σ·n̂} at quadrature point q
            DenseMatrix dshape1_ref(ndof1, dim), dshape2_ref(ndof2, dim);
            fe1->CalcDShape(eip1_q, dshape1_ref);
            fe2->CalcDShape(eip2_q, dshape2_ref);
            DenseMatrix dshape1_phys(ndof1, dim), dshape2_phys(ndof2, dim);
            Mult(dshape1_ref, Jinv1, dshape1_phys);
            Mult(dshape2_ref, Jinv2, dshape2_phys);

            DenseMatrix grad1(dim, dim), grad2(dim, dim);
            grad1 = 0.0; grad2 = 0.0;
            for (int c = 0; c < dim; c++)
               for (int d = 0; d < dim; d++)
               {
                  for (int k = 0; k < ndof1; k++)
                     grad1(c, d) += dshape1_phys(k, d)
                                    * u1_all(c * ndof1 + k);
                  for (int k = 0; k < ndof2; k++)
                     grad2(c, d) += dshape2_phys(k, d)
                                    * u2_all(c * ndof2 + k);
               }

            real_t T_stress_q[3] = {0.0, 0.0, 0.0};
            for (int ci = 0; ci < dim; ci++)
            {
               for (int cj = 0; cj < dim; cj++)
               {
                  real_t ag = 0.5 * (grad1(ci, cj) + grad2(ci, cj));
                  real_t ag_t = 0.5 * (grad1(cj, ci) + grad2(cj, ci));
                  real_t eps_ij = 0.5 * (ag + ag_t);
                  real_t tr_contrib = (ci == cj)
                     ? lambda_val_ * (0.5*((grad1(0,0)+grad2(0,0))
                        + (grad1(1,1)+grad2(1,1))
                        + (grad1(2,2)+grad2(2,2)))) : 0.0;
                  real_t stress_ij = tr_contrib + 2.0 * mu_val_ * eps_ij;
                  T_stress_q[ci] += stress_ij * basis.normal[cj];
               }
            }

            // BR2 correction at quadrature point q:
            // Evaluate lifted function using per-point shapes (not avg_shapes)
            real_t correction_q[3] = {0.0, 0.0, 0.0};
            for (int i = 0; i < dim; i++)
            {
               real_t sum = 0.0;
               for (int u = 0; u < dim; u++)
                  for (int s = 0; s < dim; s++)
                  {
                     real_t tn = lambda_val_
                           * (u == s ? 1.0 : 0.0) * basis.normal[i]
                        + mu_val_
                           * ((i == u ? 1.0 : 0.0) * basis.normal[s]
                              + (i == s ? 1.0 : 0.0) * basis.normal[u]);
                     real_t eval1 = 0.0, eval2 = 0.0;
                     for (int m = 0; m < ndof1; m++)
                        eval1 += shapes1(m, q)
                                 * f_lifted1(u * dim + s, m);
                     for (int m = 0; m < ndof2; m++)
                        eval2 += shapes2(m, q)
                                 * f_lifted2(u * dim + s, m);
                     sum += tn * (eval1 + eval2);
                  }
               correction_q[i] = br2_penalty * 0.5 * sum;
            }

            // Accumulate face average
            for (int c = 0; c < dim; c++)
            {
               T_stress[c] += wq * T_stress_q[c];
               correction[c] += wq * correction_q[c];
            }
            sum_wq += wq;
         }
      }

      // BR2 path: normalize face average and store (BR2 uses nbf=1 always)
      if (method_ != DGMethod::IP)
      {
         if (sum_wq > 0.0)
         {
            for (int c = 0; c < dim; c++)
            {
               T_stress[c] /= sum_wq;
               correction[c] /= sum_wq;
               T_global[c] = T_stress[c] - correction[c];
            }
         }

         // Project to local frame: (tau_dip, tau_strike)
         // BR2 always has nbf_per_face_=1, so dof_idx = fi
         real_t tau_local[2];
         fault_basis_.ProjectTraction(fi, T_global, tau_local);
         traction(2 * fi)     = tau_local[0];
         traction(2 * fi + 1) = tau_local[1];
         if (traction_stress_out || traction_correction_out)
         {
            real_t tau_stress_local[2], tau_corr_local[2];
            fault_basis_.ProjectTraction(fi, T_stress, tau_stress_local);
            real_t corr_neg[3] = {-correction[0], -correction[1], -correction[2]};
            fault_basis_.ProjectTraction(fi, corr_neg, tau_corr_local);
            if (traction_stress_out)
            {
               (*traction_stress_out)(2 * fi) = tau_stress_local[0];
               (*traction_stress_out)(2 * fi + 1) = tau_stress_local[1];
            }
            if (traction_correction_out)
            {
               (*traction_correction_out)(2 * fi) = tau_corr_local[0];
               (*traction_correction_out)(2 * fi + 1) = tau_corr_local[1];
            }
         }
         if (jump_residual_out)
         {
            real_t res_local[2];
            fault_basis_.ProjectTraction(fi, residual_avg, res_local);
            (*jump_residual_out)(2 * fi) = res_local[0];
            (*jump_residual_out)(2 * fi + 1) = res_local[1];
         }
         // v51: elastic normal traction for sigma_n feedback
         if (normal_traction)
         {
            (*normal_traction)(fi) = fault_basis_.NormalStress(fi, T_global);
         }
      }
   }

   // Shared fault faces (parallel only)
   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
      int rank;
      MPI_Comm_rank(mesh_.GetComm(), &rank);
      ParGridFunction par_u(pfes);
      par_u = displacement;
      pfes->ExchangeFaceNbrData();
      par_u.ExchangeFaceNbrData();

      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         int sf = fault_shared_faces_[i];
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(sf);
         if (FTr == nullptr) { continue; }

         int trac_idx = fault_interior_faces_.Size() + i;

         // v35: Per-quadrature-point traction (same approach as interior faces)
         const IntegrationPoint &ip =
            Geometries.GetCenter(FTr->GetGeometryType());

         // Elem1 (local)
         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         int ndof1 = fe1->GetDof();

         Array<int> vdofs1;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         Vector u1_all(vdofs1.Size());
         par_u.GetSubVector(vdofs1, u1_all);

         // Elem2 (face-neighbor)
         int nbr_idx = FTr->Elem2No - mesh_.GetNE();
         const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
         int ndof2 = fe2->GetDof();

         Array<int> vdofs2;
         pfes->GetFaceNbrElementVDofs(nbr_idx, vdofs2);
         const Vector &nbr_data = par_u.FaceNbrData();
         Vector u2_all(vdofs2.Size());
         for (int j = 0; j < vdofs2.Size(); j++)
         {
            u2_all(j) = nbr_data(vdofs2[j]);
         }

         // Fault basis
         const auto &basis = fault_basis_.GetBasis(trac_idx);

         // Element Jacobian inverses (constant for linear tets)
         DenseMatrix Jinv1(dim), Jinv2(dim);
         CalcInverse(FTr->Elem1->Jacobian(), Jinv1);
         CalcInverse(FTr->Elem2->Jacobian(), Jinv2);

         // Face quadrature rule
         int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
         int quad_order_trac_sh = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
         const IntegrationRule &ir_trac = IntRules.Get(
            FTr->GetGeometryType(), quad_order_trac_sh);
         int nqp = ir_trac.GetNPoints();

         // Face-averaged accumulators (used for BR2 path and diagnostics)
         real_t T_global[3] = {0.0, 0.0, 0.0};
         real_t T_stress[3] = {0.0, 0.0, 0.0};
         real_t correction[3] = {0.0, 0.0, 0.0};
         real_t residual_avg[3] = {0.0, 0.0, 0.0};
         real_t sum_wq = 0.0;


         if (method_ == DGMethod::IP)
         {
            // Tandem-style traction for shared faces (same as interior).
            // When decomposition is requested, uses decomposed variant.
            int nbf_sh = nbf_per_face_;
            bool need_decomp_sh = traction_stress_out ||
                                  traction_correction_out || jump_residual_out ||
                                  normal_stress_out || normal_correction_out;
            const auto &basis_sh = fault_basis_.GetBasis(trac_idx);

            // Build sign-corrected slip at quad points (shared faces).
            // Same Tandem-style per-QP tangent embedding as interior faces.
            int quad_order_sh_slip = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
            const IntegrationRule &ir_sh_slip = IntRules.Get(
               FTr->GetGeometryType(), quad_order_sh_slip);
            int nqp_sh_slip = ir_sh_slip.GetNPoints();
            Vector delta_u_quad_sh;

            if (!basis_sh.qp_data.empty())
            {
               // Per-QP tangent embedding (Tandem evaluate_slip)
               Vector slip_tang_sh(2 * nbf_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_sh + kk;
                  slip_tang_sh(0 * nbf_sh + kk) = slip_bc(2 * dof_idx);      // dip
                  slip_tang_sh(1 * nbf_sh + kk) = slip_bc(2 * dof_idx + 1);  // strike
               }
               Vector slip_tang_q_sh;
               face_quad_->InterpolateToQuadPoints(2, slip_tang_sh, slip_tang_q_sh);
               delta_u_quad_sh.SetSize(dim * nqp_sh_slip);
               for (int q = 0; q < nqp_sh_slip; q++)
               {
                  real_t sl_q[2] = {slip_tang_q_sh(q),
                                    slip_tang_q_sh(nqp_sh_slip + q)};
                  real_t du[3];
                  fault_basis_.EmbedSlipQP(trac_idx, q, sl_q, du);
                  for (int c = 0; c < dim; c++)
                     delta_u_quad_sh(c * nqp_sh_slip + q) = du[c];
               }
            }
            else
            {
               // Fallback: centroid tangents (BR2 path)
               Vector delta_u_nodal_sh(dim * nbf_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_sh + kk;
                  real_t sl[2] = {slip_bc(2 * dof_idx), slip_bc(2 * dof_idx + 1)};
                  real_t du[3];
                  fault_basis_.EmbedSlip(trac_idx, sl, du);
                  for (int c = 0; c < dim; c++)
                     delta_u_nodal_sh(c * nbf_sh + kk) = du[c];
               }
               face_quad_->InterpolateToQuadPoints(dim, delta_u_nodal_sh, delta_u_quad_sh);
            }

            DGElasticityIPCombinedIntegrator trac_integ_sh(
               lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
            Vector T_quad_sh, nl_q_sh;
            Vector T_stress_quad_sh, T_corr_quad_sh, R_quad_sh;
            if (need_decomp_sh)
            {
               trac_integ_sh.ComputeTractionAtQuadPointsDecomposed(
                  *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad_sh,
                  T_quad_sh,
                  traction_stress_out ? &T_stress_quad_sh : nullptr,
                  traction_correction_out ? &T_corr_quad_sh : nullptr,
                  jump_residual_out ? &R_quad_sh : nullptr,
                  nullptr, &nl_q_sh);
            }
            else
            {
               trac_integ_sh.ComputeTractionAtQuadPoints(
                  *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad_sh,
                  T_quad_sh, nullptr, &nl_q_sh);
            }
            int nqp_sh = T_quad_sh.Size() / dim;

            int quad_order_sh2 = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
            const IntegrationRule &ir_sh2 = IntRules.Get(
               FTr->GetGeometryType(), quad_order_sh2);
            const DenseMatrix &e_q_sh = face_quad_->BasisAtQuadPoints();

            const auto *qpd_sh = basis_sh.qp_data.empty() ? nullptr : &basis_sh.qp_data;

            // Unified projection for shared faces (same as interior)
            int ncomp_sh = normal_traction ? 3 : 2;
            real_t basis_vecs_sh[3][3];
            if (normal_traction)
            {
               for (int d = 0; d < 3; d++)
               {
                  basis_vecs_sh[0][d] = basis_sh.normal[d];
                  basis_vecs_sh[1][d] = basis_sh.tangent1[d];
                  basis_vecs_sh[2][d] = basis_sh.tangent2[d];
               }
            }
            else
            {
               for (int d = 0; d < 3; d++)
               {
                  basis_vecs_sh[0][d] = basis_sh.tangent1[d];
                  basis_vecs_sh[1][d] = basis_sh.tangent2[d];
                  basis_vecs_sh[2][d] = 0.0;
               }
            }

            Vector trac_local_sh;
            DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
               dim, ncomp_sh, T_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
               basis_vecs_sh, trac_local_sh, qpd_sh);

            int tang_off_sh = normal_traction ? 1 : 0;
            for (int kk = 0; kk < nbf_sh; kk++)
            {
               int dof_idx = trac_idx * nbf_per_face_ + kk;
               traction(2 * dof_idx)     = trac_local_sh((tang_off_sh + 0) * nbf_sh + kk);
               traction(2 * dof_idx + 1) = trac_local_sh((tang_off_sh + 1) * nbf_sh + kk);
            }
            if (normal_traction)
            {
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*normal_traction)(dof_idx) = -trac_local_sh(0 * nbf_sh + kk);
               }
            }

            // Project decomposition to fault DOFs (shared faces)
            if (traction_stress_out)
            {
               Vector stress_local_sh;
               DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
                  dim, 2, T_stress_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
                  basis_vecs_sh + tang_off_sh, stress_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*traction_stress_out)(2 * dof_idx)     = stress_local_sh(0 * nbf_sh + kk);
                  (*traction_stress_out)(2 * dof_idx + 1) = stress_local_sh(1 * nbf_sh + kk);
               }
            }
            if (traction_correction_out)
            {
               Vector corr_local_sh;
               DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
                  dim, 2, T_corr_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
                  basis_vecs_sh + tang_off_sh, corr_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*traction_correction_out)(2 * dof_idx)     = corr_local_sh(0 * nbf_sh + kk);
                  (*traction_correction_out)(2 * dof_idx + 1) = corr_local_sh(1 * nbf_sh + kk);
               }
            }
            if (jump_residual_out)
            {
               Vector res_local_sh;
               DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
                  dim, 2, R_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
                  basis_vecs_sh + tang_off_sh, res_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*jump_residual_out)(2 * dof_idx)     = res_local_sh(0 * nbf_sh + kk);
                  (*jump_residual_out)(2 * dof_idx + 1) = res_local_sh(1 * nbf_sh + kk);
               }
            }

            // Normal decomposition for shared faces
            if (normal_stress_out && normal_traction)
            {
               Vector ns_local_sh;
               DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
                  dim, 1, T_stress_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
                  basis_vecs_sh, ns_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*normal_stress_out)(dof_idx) = -ns_local_sh(kk);
               }
            }
            if (normal_correction_out && normal_traction)
            {
               Vector nc_local_sh;
               DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
                  dim, 1, T_corr_quad_sh, nl_q_sh, ir_sh2, nbf_sh, e_q_sh,
                  basis_vecs_sh, nc_local_sh, qpd_sh);
               for (int kk = 0; kk < nbf_sh; kk++)
               {
                  int dof_idx = trac_idx * nbf_per_face_ + kk;
                  (*normal_correction_out)(dof_idx) = -nc_local_sh(kk);
               }
            }

         }
         else  // BR2
         {
            // BR2 always uses nbf=1, so face-index == DOF-index
            real_t slip_local_br2[2] = {slip_bc(2 * trac_idx),
                                        slip_bc(2 * trac_idx + 1)};
            real_t delta_u[3];
            fault_basis_.EmbedSlip(trac_idx, slip_local_br2, delta_u);

            if (!mass_inv_computed_) {
               PrecomputeMassInverse();
               // Re-fetch FTr: PrecomputeMassInverse invalidates cached transforms.
               FTr = mesh_.GetSharedFaceTransformations(sf);
            }

            Geometry::Type geom = mesh_.GetElementGeometry(FTr->Elem1No);
            real_t br2_penalty = (geom == Geometry::TETRAHEDRON)
                                     ? real_t(dim + 1) : real_t(2 * dim);

            const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
            const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];

            DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
            Vector nor_arr(dim * nqp), w_arr(nqp);

            DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
            face_int1 = 0.0;
            face_int2 = 0.0;

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &fip = ir_trac.IntPoint(q);
               FTr->SetAllIntPoints(&fip);
               const IntegrationPoint &eip1_q = FTr->GetElement1IntPoint();
               const IntegrationPoint &eip2_q = FTr->GetElement2IntPoint();

               Vector s1q(shapes1.GetColumn(q), ndof1);
               Vector s2q(shapes2.GetColumn(q), ndof2);
               fe1->CalcShape(eip1_q, s1q);
               fe2->CalcShape(eip2_q, s2q);

               Vector nor_q(&nor_arr[q * dim], dim);
               CalcOrtho(FTr->Jacobian(), nor_q);
               w_arr[q] = fip.weight;

               real_t jump_q[3] = {0.0, 0.0, 0.0};
               for (int c = 0; c < dim; c++)
               {
                  real_t u1q = 0.0, u2q = 0.0;
                  for (int k = 0; k < ndof1; k++)
                     u1q += s1q(k) * u1_all(c * ndof1 + k);
                  for (int k = 0; k < ndof2; k++)
                     u2q += s2q(k) * u2_all(c * ndof2 + k);
                  jump_q[c] = (u1q - u2q) - delta_u[c];
               }

               for (int u = 0; u < dim; u++)
                  for (int s = 0; s < dim; s++)
                  {
                     for (int m = 0; m < ndof1; m++)
                        face_int1(u * dim + s, m) +=
                           w_arr[q] * s1q(m) * jump_q[u] * nor_q(s);
                     for (int m = 0; m < ndof2; m++)
                        face_int2(u * dim + s, m) +=
                           w_arr[q] * s2q(m) * jump_q[u] * nor_q(s);
                  }
            }

            DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
            MultABt(face_int1, Minv1, f_lifted1);
            f_lifted1 *= 0.5;
            MultABt(face_int2, Minv2, f_lifted2);
            f_lifted2 *= 0.5;

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &fip = ir_trac.IntPoint(q);
               FTr->SetAllIntPoints(&fip);
               const IntegrationPoint &eip1_q = FTr->GetElement1IntPoint();
               const IntegrationPoint &eip2_q = FTr->GetElement2IntPoint();
               real_t wq = w_arr[q];

               DenseMatrix dshape1_ref(ndof1, dim), dshape2_ref(ndof2, dim);
               fe1->CalcDShape(eip1_q, dshape1_ref);
               fe2->CalcDShape(eip2_q, dshape2_ref);
               DenseMatrix dshape1_phys(ndof1, dim), dshape2_phys(ndof2, dim);
               Mult(dshape1_ref, Jinv1, dshape1_phys);
               Mult(dshape2_ref, Jinv2, dshape2_phys);

               DenseMatrix grad1(dim, dim), grad2(dim, dim);
               grad1 = 0.0; grad2 = 0.0;
               for (int c = 0; c < dim; c++)
                  for (int d = 0; d < dim; d++)
                  {
                     for (int k = 0; k < ndof1; k++)
                        grad1(c, d) += dshape1_phys(k, d)
                                       * u1_all(c * ndof1 + k);
                     for (int k = 0; k < ndof2; k++)
                        grad2(c, d) += dshape2_phys(k, d)
                                       * u2_all(c * ndof2 + k);
                  }

               real_t T_stress_q[3] = {0.0, 0.0, 0.0};
               for (int ci = 0; ci < dim; ci++)
                  for (int cj = 0; cj < dim; cj++)
                  {
                     real_t ag = 0.5 * (grad1(ci, cj) + grad2(ci, cj));
                     real_t ag_t = 0.5 * (grad1(cj, ci) + grad2(cj, ci));
                     real_t eps_ij = 0.5 * (ag + ag_t);
                     real_t tr_contrib = (ci == cj)
                        ? lambda_val_ * (0.5*((grad1(0,0)+grad2(0,0))
                           + (grad1(1,1)+grad2(1,1))
                           + (grad1(2,2)+grad2(2,2)))) : 0.0;
                     real_t stress_ij = tr_contrib + 2.0 * mu_val_ * eps_ij;
                     T_stress_q[ci] += stress_ij * basis.normal[cj];
                  }

               real_t correction_q[3] = {0.0, 0.0, 0.0};
               for (int ci = 0; ci < dim; ci++)
               {
                  real_t sum = 0.0;
                  for (int u = 0; u < dim; u++)
                     for (int s = 0; s < dim; s++)
                     {
                        real_t tn = lambda_val_
                              * (u == s ? 1.0 : 0.0) * basis.normal[ci]
                           + mu_val_
                              * ((ci == u ? 1.0 : 0.0) * basis.normal[s]
                                 + (ci == s ? 1.0 : 0.0) * basis.normal[u]);
                        real_t eval1 = 0.0, eval2 = 0.0;
                        for (int m = 0; m < ndof1; m++)
                           eval1 += shapes1(m, q)
                                    * f_lifted1(u * dim + s, m);
                        for (int m = 0; m < ndof2; m++)
                           eval2 += shapes2(m, q)
                                    * f_lifted2(u * dim + s, m);
                        sum += tn * (eval1 + eval2);
                     }
                  correction_q[ci] = br2_penalty * 0.5 * sum;
               }

               for (int c = 0; c < dim; c++)
               {
                  T_stress[c] += wq * T_stress_q[c];
                  correction[c] += wq * correction_q[c];
               }
               sum_wq += wq;
            }
         }

         // BR2 path: normalize face average and store (BR2 uses nbf=1 always)
         if (method_ != DGMethod::IP)
         {
            if (sum_wq > 0.0)
            {
               for (int c = 0; c < dim; c++)
               {
                  T_stress[c] /= sum_wq;
                  correction[c] /= sum_wq;
                  residual_avg[c] /= sum_wq;
                  T_global[c] = T_stress[c] - correction[c];
               }
            }

            // BR2 always has nbf_per_face_=1, so dof_idx = trac_idx
            real_t tau_local[2];
            fault_basis_.ProjectTraction(trac_idx, T_global, tau_local);
            traction(2 * trac_idx)     = tau_local[0];
            traction(2 * trac_idx + 1) = tau_local[1];
            if (traction_stress_out || traction_correction_out)
            {
               real_t tau_stress_local[2], tau_corr_local[2];
               fault_basis_.ProjectTraction(trac_idx, T_stress, tau_stress_local);
               real_t corr_neg[3] = {-correction[0], -correction[1], -correction[2]};
               fault_basis_.ProjectTraction(trac_idx, corr_neg, tau_corr_local);
               if (traction_stress_out)
               {
                  (*traction_stress_out)(2 * trac_idx) = tau_stress_local[0];
                  (*traction_stress_out)(2 * trac_idx + 1) = tau_stress_local[1];
               }
               if (traction_correction_out)
               {
                  (*traction_correction_out)(2 * trac_idx) = tau_corr_local[0];
                  (*traction_correction_out)(2 * trac_idx + 1) = tau_corr_local[1];
               }
            }
            if (jump_residual_out)
            {
               real_t res_local[2];
               fault_basis_.ProjectTraction(trac_idx, residual_avg, res_local);
               (*jump_residual_out)(2 * trac_idx) = res_local[0];
               (*jump_residual_out)(2 * trac_idx + 1) = res_local[1];
            }
            // v51: elastic normal traction for sigma_n feedback
            if (normal_traction)
            {
               (*normal_traction)(trac_idx) =
                  fault_basis_.NormalStress(trac_idx, T_global);
            }
         }
      }
#endif
   }

   // Absolute traction sanity guard.  The post-solve relative residual check
   // (||K x-b||/||b||) CANNOT catch a near-null-space-contaminated solution: a
   // contaminated x = x_true + c*v (K v ~ eps*v, eps tiny) keeps ||K x-b|| small
   // while c blows up, so the relative residual passes but the traction (which
   // depends on grad(x), amplified by the DG penalty) is enormous.  A loaded CG
   // solve that STAGNATES on the ill-conditioned DG operator (e.g. with
   // aggressive coarsening) lands exactly here and would otherwise feed garbage
   // (|tau| ~ 1e18, sigma_n_eff ~ 1e17) silently into the friction solver as a
   // bare NaN.  Abort here, at the source, with diagnostics.  Physical QD
   // traction is O(sigma_n) ~ 1e8 Pa; 1e12 (~1e4x) is a safe non-physical ceiling
   // (above any legitimate imposed-slip transient, far below garbage ~1e18).
#ifdef MFEM_USE_MPI
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      int rank;
      MPI_Comm_rank(mesh_.GetComm(), &rank);
      const int num_interior_faces = fault_interior_faces_.Size();
      real_t local_max_tau = 0.0;
      for (int i = 0; i < num_fault_dofs_; i++)
      {
         real_t tau_mag = std::sqrt(traction(2*i)*traction(2*i) +
                                    traction(2*i+1)*traction(2*i+1));
         if (std::isfinite(tau_mag)) { local_max_tau = std::max(local_max_tau, tau_mag); }
         else                        { local_max_tau = 1e300; } // sentinel >> 1e9 for NaN/Inf
         if (tau_mag > 1e12 || std::isnan(tau_mag))
         {
            // Convert DOF index to face index (nbf_per_face_ DOFs per face)
            int face_i = i / nbf_per_face_;
            // Get face coordinates for diagnostics
            Vector face_center(3);
            face_center = 0.0;
            if (face_i < num_interior_faces)
            {
               int face_idx = fault_interior_faces_[face_i];
               FaceElementTransformations *FTr =
                  mesh_.GetInteriorFaceTransformations(face_idx);
               if (FTr)
               {
                  const IntegrationPoint &ip =
                     Geometries.GetCenter(FTr->GetGeometryType());
                  FTr->Face->SetIntPoint(&ip);
                  FTr->Face->Transform(ip, face_center);
               }
            }
            else
            {
               int shared_idx = face_i - num_interior_faces;
               int sf = fault_shared_faces_[shared_idx];
               FaceElementTransformations *FTr =
                  mesh_.GetSharedFaceTransformations(sf);
               if (FTr)
               {
                  const IntegrationPoint &ip =
                     Geometries.GetCenter(FTr->GetGeometryType());
                  FTr->Face->SetIntPoint(&ip);
                  FTr->Face->Transform(ip, face_center);
               }
            }
            mfem::out << "[Rank " << rank << "] TRACTION BLOWUP: DOF " << i
                      << " (face " << face_i << "/"
                      << (face_i < num_interior_faces ? "interior" : "shared")
                      << ")"
                      << " tau_mag=" << tau_mag
                      << " tau=(" << traction(2*i) << ","
                      << traction(2*i+1) << ")"
                      << " at x=(" << face_center(0) << ","
                      << face_center(1) << "," << face_center(2) << ")"
                      << " slip=(" << slip_bc(2*i) << ","
                      << slip_bc(2*i+1) << ")\n";
         }
      }
      // Collective abort if ANY rank saw a non-physical traction.  All ranks
      // reduce so the MFEM_VERIFY fires consistently (avoids a one-rank hang).
      real_t global_max_tau = local_max_tau;
      MPI_Allreduce(MPI_IN_PLACE, &global_max_tau, 1, MPI_DOUBLE,
                    MPI_MAX, mesh_.GetComm());
      MFEM_VERIFY(global_max_tau <= 1e12,
                  "Domain solve produced non-physical fault traction (max |tau| = "
                  << global_max_tau << " Pa >> O(sigma_n) ~ 1e8).  This is a "
                  "near-null-space-contaminated / non-converged elasticity solve "
                  "(the relative residual check cannot detect it).  Use a stronger "
                  "preconditioner: set [solver].amg_aggressive_levels=0 (default), "
                  "raise ksp_maxit, or use a direct solver (mumps).");
   }
#endif

}
