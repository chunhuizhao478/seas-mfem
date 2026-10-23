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

// elasticity_operator_qd.inl — GPU plan G2: the per-stage face work of the
// quasi-dynamic solve as operators (qd_stage_operators.hpp).
//
// BuildQDStageOperators probes the face-local functions the face loops use
// (BuildSlipAtQuadPointsFace + AssembleSlipFaceRHS, AssembleBoundaryFaceRHS,
// ComputeIPFaceTraction) with unit inputs; each response is one column of
// that face's block.  The face sets, quadrature rules, Dirichlet points and
// skeleton signs are the loops' own.  Everything is linear, so the matrices
// reproduce the loops up to round-off, without the loops' small-slip skip.
//
// Index spaces:
//   volume   true DOFs of the DG space (local DOFs = true DOFs)
//   slip     owned fault DOFs x [dip, strike]; global index 2*g + c
//   normal   owned fault DOFs; global index g
//   g(t)     this rank's Dirichlet points x 3 components; 3*p + c + offset
// g = global owned index of a fault DOF.  Every local copy of a fault DOF
// (shared faces have one on each side) learns its owner's g by expanding a
// vector holding g + 1 on the owned DOFs with ExpandOwnedToLocalFault, so the
// canonical face-DOF permutation is the one the loops use.

// ---------------------------------------------------------------------------
// Face-local IP traction (shared with ComputeTractionImpl)
// ---------------------------------------------------------------------------

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeIPFaceTraction(
   int fi, FaceElementTransformations &FTr,
   const FiniteElement &fe1, const FiniteElement &fe2,
   const Vector &u1_all, const Vector &u2_all,
   const real_t *face_slip, bool with_normal, Vector &trac_local) const
{
   const int dim = 3;
   const int nbf = nbf_per_face_;
   const auto &basis = fault_basis_.GetBasis(fi);

   // Sign-corrected slip at quad points (Tandem evaluate_slip)
   Vector delta_u_quad;
   BuildSlipAtQuadPointsFace(fi, face_slip, delta_u_quad);

   // Traction at quad points (Tandem compute_traction)
   DGElasticityIPCombinedIntegrator trac_integ(
      *lambda_ptr_, *mu_ptr_, dim, epsilon_, penalty_factor_);
   Vector T_quad, nl_q;
   trac_integ.ComputeTractionAtQuadPoints(fe1, fe2, FTr, u1_all, u2_all,
                                          delta_u_quad, T_quad, nullptr,
                                          &nl_q);

   // Projection to fault DOFs (Tandem evaluate_traction): normal, dip and
   // strike through one L2 solve, or dip and strike only.
   const int quad_order = 2 * std::max(fe1.GetOrder(), fe2.GetOrder()) + 1;
   const IntegrationRule &ir = IntRules.Get(FTr.GetGeometryType(), quad_order);
   const DenseMatrix &e_q = face_quad_->BasisAtQuadPoints();
   const auto *qpd = basis.qp_data.empty() ? nullptr : &basis.qp_data;
   real_t basis_vecs[3][3];
   for (int d = 0; d < 3; d++)
   {
      if (with_normal)
      {
         basis_vecs[0][d] = basis.normal[d];
         basis_vecs[1][d] = basis.tangent1[d];
         basis_vecs[2][d] = basis.tangent2[d];
      }
      else
      {
         basis_vecs[0][d] = basis.tangent1[d];
         basis_vecs[1][d] = basis.tangent2[d];
         basis_vecs[2][d] = 0.0;
      }
   }
   DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
      dim, with_normal ? 3 : 2, T_quad, nl_q, ir, nbf, e_q, basis_vecs,
      trac_local, qpd);
}

#ifdef MFEM_USE_MPI

// ---------------------------------------------------------------------------
// Build
// ---------------------------------------------------------------------------

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::BuildQDStageOperators(
   bool with_normal)
{
   if constexpr (!IsParallelMesh<MeshType>::value)
   {
      MFEM_ABORT("--qd-operators matrix needs a parallel mesh (with_normal = "
                 << with_normal << ")");
   }
   else
   {
      MFEM_VERIFY(method_ == DGMethod::IP,
                  "--qd-operators matrix: IP method only (BR2 is not "
                  "supported)");
      MPI_Comm comm = mesh_.GetComm();
      const double t_start = MPI_Wtime();

      auto *pfes = dynamic_cast<ParFiniteElementSpace *>(fes_.get());
      MFEM_VERIFY(pfes != nullptr, "QD stage operators: no parallel space");
      MFEM_VERIFY(pfes->GetVSize() == pfes->GetTrueVSize(),
                  "QD stage operators: DG space expected (local DOFs = "
                  "true DOFs)");
      pfes->ExchangeFaceNbrData();
      const HYPRE_BigInt *nbr_gdof = pfes->GetFaceNbrGlobalDofMap();

      const int dim = 3;
      const int nbf = nbf_per_face_;
      const int n_int = fault_interior_faces_.Size();
      const int n_slots = n_int + fault_shared_faces_.Size();
      const int nvol = pfes->GetTrueVSize();
      const int nown = num_owned_fault_dofs_;

      // ---- Index spaces ----
      Array<HYPRE_BigInt> vol_part, slip_part, norm_part, g_part;
      HYPRE_BigInt vol_first, fault_first, slip_first, g_first;
      const HYPRE_BigInt nvol_glob =
         MakeQDPartition(comm, nvol, vol_part, vol_first);
      const HYPRE_BigInt nfault_glob =
         MakeQDPartition(comm, nown, norm_part, fault_first);
      const HYPRE_BigInt nslip_glob =
         MakeQDPartition(comm, 2 * nown, slip_part, slip_first);
      MFEM_VERIFY(vol_first == pfes->GetMyTDofOffset() &&
                  nvol_glob == pfes->GlobalTrueVSize(),
                  "QD stage operators: volume partition differs from the "
                  "space's true-DOF partition");

      // Owner's global index of every local fault DOF copy
      Vector owned_id(nown), local_id;
      for (int i = 0; i < nown; i++)
      {
         owned_id(i) = static_cast<real_t>(fault_first + i + 1);
      }
      ExpandOwnedToLocalFault(owned_id, local_id, 1);
      Array<HYPRE_BigInt> fault_gid(num_fault_dofs_);
      for (int d = 0; d < num_fault_dofs_; d++)
      {
         const real_t v = local_id(d);
         MFEM_VERIFY(v >= 1.0 && v == std::floor(v) &&
                     v <= static_cast<real_t>(nfault_glob),
                     "QD stage operators: local fault DOF " << d
                     << " got no owner index (" << v << ")");
         fault_gid[d] = static_cast<HYPRE_BigInt>(v) - 1;
      }
      // Owned index of each local fault DOF (-1: a copy owned elsewhere)
      Array<int> local_to_owned(num_fault_dofs_);
      local_to_owned = -1;
      for (int i = 0; i < nown; i++)
      {
         const int ld = owned_fault_dof_to_local_dof_[i];
         local_to_owned[ld] = i;
         MFEM_VERIFY(fault_gid[ld] == fault_first + i,
                     "QD stage operators: owned fault DOF " << i
                     << " expands to global index " << fault_gid[ld]);
      }

      DGElasticityIPCombinedIntegrator integ(
         *lambda_ptr_, *mu_ptr_, dim, epsilon_, penalty_factor_);
      QDTriplets tB, tD, tCu, tCd, tNu, tNd;
      Array<int> vdofs1, vdofs2;
      Vector face_slip(2 * nbf), delta, elvec1, elvec2, trac_local;

      // Row (or column) of a local vdof, with the sign of a negative vdof
      auto vdof_index = [](int v, real_t &sign)
      {
         if (v >= 0) { sign = 1.0; return v; }
         sign = -1.0;
         return -1 - v;
      };
      auto add_rhs = [&](QDTriplets &t, const Array<int> &vdofs,
                         const Vector &ev, HYPRE_BigInt col)
      {
         for (int j = 0; j < vdofs.Size(); j++)
         {
            real_t sign;
            const int row = vdof_index(vdofs[j], sign);
            t.Add(row, col, sign * ev(j));
         }
      };

      // ---- Fault faces: B, C_u, C_δ, N_u, N_δ ----
      for (int s = 0; s < n_slots; s++)
      {
         const bool shared = (s >= n_int);
         FaceElementTransformations *FTr = shared
            ? mesh_.GetSharedFaceTransformations(fault_shared_faces_[s - n_int])
            : mesh_.GetInteriorFaceTransformations(fault_interior_faces_[s]);
         if (FTr == nullptr) { continue; }

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = nullptr;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         if (shared)
         {
            const int nbr = FTr->Elem2No - mesh_.GetNE();
            fe2 = pfes->GetFaceNbrFE(nbr);
            pfes->GetFaceNbrElementVDofs(nbr, vdofs2);
         }
         else
         {
            fe2 = scalar_fes_->GetFE(FTr->Elem2No);
            fes_->GetElementVDofs(FTr->Elem2No, vdofs2);
         }

         // B: unit slip -> right-hand side of this rank's elements (shared
         // faces: Elem1 only, as AssembleSlipContributionIPShared)
         for (int kk = 0; kk < nbf; kk++)
         {
            for (int c = 0; c < 2; c++)
            {
               face_slip = 0.0;
               face_slip(2 * kk + c) = 1.0;
               BuildSlipAtQuadPointsFace(s, face_slip.GetData(), delta);
               integ.AssembleSlipFaceRHS(*fe1, *fe2, *FTr, delta,
                                         elvec1, elvec2);
               const HYPRE_BigInt col = 2 * fault_gid[s * nbf + kk] + c;
               add_rhs(tB, vdofs1, elvec1, col);
               if (!shared) { add_rhs(tB, vdofs2, elvec2, col); }
            }
         }

         // C, N: rows only where this rank owns the face
         const bool owned = (local_to_owned[s * nbf] >= 0);
         for (int kk = 1; kk < nbf; kk++)
         {
            MFEM_VERIFY((local_to_owned[s * nbf + kk] >= 0) == owned,
                        "QD stage operators: fault face slot " << s
                        << " is partly owned");
         }
         if (!owned) { continue; }

         const int tang_offset = with_normal ? 1 : 0;
         auto add_traction = [&](QDTriplets &tC, QDTriplets &tN,
                                 HYPRE_BigInt col, real_t sign)
         {
            for (int kk = 0; kk < nbf; kk++)
            {
               const int o = local_to_owned[s * nbf + kk];
               for (int c = 0; c < 2; c++)
               {
                  tC.Add(2 * o + c, col,
                         sign * trac_local((tang_offset + c) * nbf + kk));
               }
               if (with_normal)
               {
                  tN.Add(o, col, -sign * trac_local(0 * nbf + kk));
               }
            }
         };

         Vector u1(vdofs1.Size()), u2(vdofs2.Size());
         face_slip = 0.0;
         for (int j = 0; j < vdofs1.Size(); j++)
         {
            u1 = 0.0;
            u2 = 0.0;
            u1(j) = 1.0;
            ComputeIPFaceTraction(s, *FTr, *fe1, *fe2, u1, u2,
                                  face_slip.GetData(), with_normal,
                                  trac_local);
            real_t sign;
            const int ld = vdof_index(vdofs1[j], sign);
            add_traction(tCu, tNu, pfes->GetGlobalTDofNumber(ld), sign);
         }
         for (int j = 0; j < vdofs2.Size(); j++)
         {
            u1 = 0.0;
            u2 = 0.0;
            u2(j) = 1.0;
            ComputeIPFaceTraction(s, *FTr, *fe1, *fe2, u1, u2,
                                  face_slip.GetData(), with_normal,
                                  trac_local);
            if (shared)
            {
               // The loops read FaceNbrData(vdofs2[j]) directly.
               add_traction(tCu, tNu, nbr_gdof[vdofs2[j]], 1.0);
            }
            else
            {
               real_t sign;
               const int ld = vdof_index(vdofs2[j], sign);
               add_traction(tCu, tNu, pfes->GetGlobalTDofNumber(ld), sign);
            }
         }
         u1 = 0.0;
         u2 = 0.0;
         for (int kk = 0; kk < nbf; kk++)
         {
            for (int c = 0; c < 2; c++)
            {
               face_slip = 0.0;
               face_slip(2 * kk + c) = 1.0;
               ComputeIPFaceTraction(s, *FTr, *fe1, *fe2, u1, u2,
                                     face_slip.GetData(), with_normal,
                                     trac_local);
               add_traction(tCd, tNd, 2 * fault_gid[s * nbf + kk] + c, 1.0);
            }
         }
      }

      // ---- Dirichlet faces: D and the points g(t) is evaluated at ----
      // Same faces, rules, points and skeleton signs as
      // AssembleDirichletLoading; the sign goes into D, g holds the raw
      // callback values.  Columns are local (3p + c) until the offset of
      // this rank's points is known.
      std::vector<int> pt_attr;
      std::vector<real_t> pt_xyz;
      Vector u_D_3d, phys(dim);
      auto add_point = [&](FaceElementTransformations *Tr,
                           const IntegrationPoint &ipq, int attr)
      {
         Tr->SetAllIntPoints(&ipq);
         Tr->Elem1->Transform(Tr->GetElement1IntPoint(), phys);
         pt_attr.push_back(attr);
         for (int d = 0; d < dim; d++) { pt_xyz.push_back(phys(d)); }
      };

      // (1) Boundary faces (boundary-face term, Elem1)
      Array<int> dir_bdr_faces, dir_bdr_attrs;
      CollectDirichletBoundaryFaces(dir_bdr_faces, dir_bdr_attrs);
      for (int bi = 0; bi < dir_bdr_faces.Size(); bi++)
      {
         FaceElementTransformations *FTr =
            mesh_.GetFaceElementTransformations(dir_bdr_faces[bi]);
         if (FTr == nullptr) { continue; }
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         const FiniteElement *fe = scalar_fes_->GetFE(FTr->Elem1No);
         const IntegrationRule &ir_dir =
            IntRules.Get(FTr->FaceGeom, 2 * fe->GetOrder() + 1);
         const int nq = ir_dir.GetNPoints();
         const int p0 = static_cast<int>(pt_attr.size());
         for (int q = 0; q < nq; q++)
         {
            add_point(FTr, ir_dir.IntPoint(q), dir_bdr_attrs[bi]);
         }
         u_D_3d.SetSize(dim * nq);
         for (int q = 0; q < nq; q++)
         {
            for (int c = 0; c < dim; c++)
            {
               u_D_3d = 0.0;
               u_D_3d(c * nq + q) = 1.0;
               integ.AssembleBoundaryFaceRHS(*fe, *FTr, u_D_3d, elvec1);
               add_rhs(tD, vdofs1, elvec1, 3 * (p0 + q) + c);
            }
         }
      }

      // (2) Interior and (3) shared skeleton faces (slip-face term with the
      // per-point orientation sign; shared faces: Elem1 only)
      const int n_dir_int = dirichlet_interior_faces_.Size();
      const int n_dir = n_dir_int + dirichlet_shared_faces_.Size();
      Vector sign_q;
      for (int k = 0; k < n_dir; k++)
      {
         const bool shared = (k >= n_dir_int);
         FaceElementTransformations *FTr = shared
            ? mesh_.GetSharedFaceTransformations(
                 dirichlet_shared_faces_[k - n_dir_int])
            : mesh_.GetInteriorFaceTransformations(
                 dirichlet_interior_faces_[k]);
         if (FTr == nullptr) { continue; }
         const int attr = shared ? dirichlet_shared_attrs_[k - n_dir_int]
                                 : dirichlet_interior_attrs_[k];
         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = nullptr;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         if (shared)
         {
            fe2 = pfes->GetFaceNbrFE(FTr->Elem2No - mesh_.GetNE());
         }
         else
         {
            fe2 = scalar_fes_->GetFE(FTr->Elem2No);
            fes_->GetElementVDofs(FTr->Elem2No, vdofs2);
         }
         const IntegrationRule &ir_dir = IntRules.Get(
            FTr->FaceGeom, 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1);
         const int nq = ir_dir.GetNPoints();
         const int p0 = static_cast<int>(pt_attr.size());
         sign_q.SetSize(nq);
         for (int q = 0; q < nq; q++)
         {
            add_point(FTr, ir_dir.IntPoint(q), attr);
            sign_q(q) = ComputeSkeletonDirichletSign(FTr);
         }
         u_D_3d.SetSize(dim * nq);
         for (int q = 0; q < nq; q++)
         {
            for (int c = 0; c < dim; c++)
            {
               u_D_3d = 0.0;
               u_D_3d(c * nq + q) = sign_q(q);
               integ.AssembleSlipFaceRHS(*fe1, *fe2, *FTr, u_D_3d,
                                         elvec1, elvec2);
               add_rhs(tD, vdofs1, elvec1, 3 * (p0 + q) + c);
               if (!shared) { add_rhs(tD, vdofs2, elvec2, 3 * (p0 + q) + c); }
            }
         }
      }

      const int ng = 3 * static_cast<int>(pt_attr.size());
      const HYPRE_BigInt ng_glob = MakeQDPartition(comm, ng, g_part, g_first);
      tD.OffsetColumns(g_first);

      // ---- Matrices ----
      auto ops = std::unique_ptr<QDStageOperators>(new QDStageOperators);
      ops->with_normal = with_normal;
      ops->B = tB.Build(comm, nvol, nvol_glob, vol_part, nslip_glob, slip_part);
      ops->D = tD.Build(comm, nvol, nvol_glob, vol_part, ng_glob, g_part);
      ops->Cu = tCu.Build(comm, 2 * nown, nslip_glob, slip_part,
                          nvol_glob, vol_part);
      ops->Cd = tCd.Build(comm, 2 * nown, nslip_glob, slip_part,
                          nslip_glob, slip_part);
      if (with_normal)
      {
         ops->Nu = tNu.Build(comm, nown, nfault_glob, norm_part,
                             nvol_glob, vol_part);
         ops->Nd = tNd.Build(comm, nown, nfault_glob, norm_part,
                             nslip_glob, slip_part);
      }
      ops->g_attr.SetSize(static_cast<int>(pt_attr.size()));
      for (int p = 0; p < ops->g_attr.Size(); p++) { ops->g_attr[p] = pt_attr[p]; }
      ops->g_xyz.SetSize(static_cast<int>(pt_xyz.size()));
      for (int i = 0; i < ops->g_xyz.Size(); i++) { ops->g_xyz(i) = pt_xyz[i]; }
      ops->g.SetSize(ng);
      ops->g.UseDevice(true);   // D·g reads it on the device (GPU plan G3)
      ops->g = 0.0;
      ops->build_seconds = MPI_Wtime() - t_start;
      qd_ops_ = std::move(ops);
   }
}

// ---------------------------------------------------------------------------
// Stage
// ---------------------------------------------------------------------------

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::EvaluateQDDirichletValues(
   real_t time) const
{
   QDStageOperators &ops = *qd_ops_;
   const double t0 = MPI_Wtime();
   // The callback runs on the host; on a device, the previous stage left g
   // there, so mark the host copy as the one being written.
   real_t *g = ops.g.HostWrite();
   Vector u_D_val(3);
   for (int p = 0; p < ops.NumPoints(); p++)
   {
      const Vector x(ops.g_xyz.GetData() + 3 * p, 3);
      EvalDirichletFunc(ops.g_attr[p], x, time, u_D_val);
      for (int c = 0; c < 3; c++) { g[3 * p + c] = u_D_val(c); }
   }
   ops.g_eval_seconds += MPI_Wtime() - t0;
   ops.g_evals++;
   if (Device::Allows(Backend::DEVICE_MASK))
   {
      // Copy g to the device here, so the copy is timed on its own rather
      // than inside D·g (on a GPU, after the kernels queued before it).
      const bool gpu = Device::Allows(Backend::CUDA_MASK | Backend::HIP_MASK);
      if (gpu) { MFEM_DEVICE_SYNC; }
      const double t1 = MPI_Wtime();
      {
         // The one planned host-to-device copy in a stage (GPU plan G3a);
         // its own NVTX range, so the G4 trace check can tell it apart.
         NvtxRange range("seas:g_copy");
         ops.g.Read();
         if (gpu) { MFEM_DEVICE_SYNC; }
      }
      ops.g_copy_seconds += MPI_Wtime() - t1;
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::QDStageRHS(
   real_t time, const Vector &slip_owned, Vector &b) const
{
   MFEM_VERIFY(qd_ops_, "QD stage operators: BuildQDStageOperators() was "
               "not called");
   const QDStageOperators &ops = *qd_ops_;
   b.SetSize(ops.B->Height());
   {
      ScopedStageTimer timer(stage_timers_, StageTimers::SlipRHS);
      ops.B->Mult(slip_owned, b);
   }
   {
      ScopedStageTimer timer(stage_timers_, StageTimers::DirichletRHS);
      EvaluateQDDirichletValues(time);
      ops.D->Mult(1.0, ops.g, 1.0, b);
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::QDStageTraction(
   const Vector &u, const Vector &slip_owned, Vector &traction_owned,
   Vector *normal_owned) const
{
   MFEM_VERIFY(qd_ops_, "QD stage operators: BuildQDStageOperators() was "
               "not called");
   const QDStageOperators &ops = *qd_ops_;
   MFEM_VERIFY(normal_owned == nullptr || ops.with_normal,
               "QD stage operators: normal traction requested, but N_u and "
               "N_delta were not built (BuildQDStageOperators(false))");
   traction_owned.SetSize(ops.Cu->Height());
   ops.Cu->Mult(u, traction_owned);
   ops.Cd->Mult(1.0, slip_owned, 1.0, traction_owned);
   if (normal_owned)
   {
      normal_owned->SetSize(ops.Nu->Height());
      ops.Nu->Mult(u, *normal_owned);
      ops.Nd->Mult(1.0, slip_owned, 1.0, *normal_owned);
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::SolveStage(
   real_t time, const Vector &slip_owned, GridFuncType &displacement,
   Vector &traction_owned, Vector *normal_owned)
{
   MFEM_VERIFY(!assemble_matrix_only_,
               "SolveStage() after AssembleMatrixOnly(): no solver was set up");
   if (!stiffness_assembled_)
   {
      ScopedStageTimer timer(stage_timers_, StageTimers::MatrixSetup);
      AssembleStiffness();
   }
   QDStageRHS(time, slip_owned, B_);
   SolveLoadedRHS(displacement);
   ScopedStageTimer timer(stage_timers_, StageTimers::Traction);
   QDStageTraction(X_, slip_owned, traction_owned, normal_owned);
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::PrintQDStageEvalSummary(
   long long accepted_steps, long long rejected_steps, std::ostream &os) const
{
   MFEM_VERIFY(qd_ops_, "QD stage operators: BuildQDStageOperators() was "
               "not called");
   MPI_Comm comm = mesh_.GetComm();
   const QDStageOperators &ops = *qd_ops_;
   double sec_max = 0.0, pts = ops.NumPoints(), pts_max = 0.0, copy_max = 0.0;
   MPI_Reduce(&ops.g_eval_seconds, &sec_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
   MPI_Reduce(&pts, &pts_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
   MPI_Reduce(&ops.g_copy_seconds, &copy_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
   // Host loading (callback + copy) as a share of the stage time, per rank
   // (GPU plan G3a threshold: 5%); needs --stage-timers.
   const double stage_sec =
      stage_timers_ ? stage_timers_->Seconds(StageTimers::StageTotal) : 0.0;
   const double share = stage_sec > 0.0
                        ? (ops.g_eval_seconds + ops.g_copy_seconds) / stage_sec
                        : -1.0;
   double share_max = -1.0;
   MPI_Reduce(&share, &share_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
   int rank;
   MPI_Comm_rank(comm, &rank);
   if (rank != 0) { return; }
   const bool on_device = Device::Allows(Backend::DEVICE_MASK);
   const long long attempted = accepted_steps + rejected_steps;
   std::ostringstream o;
   o << "[QD-OPERATORS] g(t): " << ops.g_evals << " evaluations over "
     << accepted_steps << " accepted + " << rejected_steps
     << " rejected steps (" << std::fixed << std::setprecision(2)
     << (attempted > 0 ? double(ops.g_evals) / attempted : 0.0)
     << " per attempted step); callback " << std::scientific
     << std::setprecision(4) << sec_max << " s total, "
     << (ops.g_evals > 0 ? 1e3 * sec_max / ops.g_evals : 0.0)
     << " ms per evaluation (slowest rank, at most "
     << static_cast<long long>(pts_max) << " points = "
     << static_cast<long long>(3 * pts_max) << " values per rank); ";
   if (on_device)
   {
      o << "host-to-device copy " << copy_max << " s total, "
        << (ops.g_evals > 0 ? 1e3 * copy_max / ops.g_evals : 0.0)
        << " ms per copy (slowest rank)\n";
   }
   else
   {
      o << "no host-to-device copy on the CPU\n";
   }
   if (share_max >= 0.0)
   {
      o << "[QD-OPERATORS] host loading (g(t) callback + copy): "
        << std::fixed << std::setprecision(3) << 100.0 * share_max
        << "% of the stage time (largest share over ranks; G3a threshold "
        << "5%)\n";
   }
   os << o.str();
   os.flush();
}

// ---------------------------------------------------------------------------
// Green's function (QD speed item 3, --qd-green)
// ---------------------------------------------------------------------------

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::QDSolvePathTraction(
   real_t time, const Vector &slip_owned, Vector &T, Vector *Tn)
{
   StageTimers *timers = stage_timers_;
   stage_timers_ = nullptr;
   Vector b, x;
   QDStageRHS(time, slip_owned, b);
   x.SetSize(b.Size());
   x = 0.0;
   solver_->Mult(b, x);
   QDStageTraction(x, slip_owned, T, Tn);
   stage_timers_ = timers;
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::BuildQDGreen(
   const QDGreenConfig &config)
{
   if constexpr (!IsParallelMesh<MeshType>::value)
   {
      MFEM_ABORT("--qd-green needs a parallel mesh (block = "
                 << config.block << ")");
   }
   else
   {
      MFEM_VERIFY(qd_ops_, "--qd-green needs the QD stage operators: call "
                  "BuildQDStageOperators() first (--qd-operators matrix)");
      MFEM_VERIFY(!Device::Allows(Backend::DEVICE_MASK),
                  "--qd-green runs on the CPU only (no --device)");
      MFEM_VERIFY(config.compress == 0.0 ||
                  !Device::Allows(Backend::DEVICE_MASK),
                  "--qd-green-compress runs on the CPU only (no --device)");
      MFEM_VERIFY(config.block >= 1, "--qd-green-block must be >= 1, got "
                  << config.block);
      if (!stiffness_assembled_)
      {
         ScopedStageTimer timer(stage_timers_, StageTimers::MatrixSetup);
         AssembleStiffness();
      }
      MFEM_VERIFY(solver_, "--qd-green: no linear solver was set up");
      MPI_Comm comm = mesh_.GetComm();
      int rank, nranks;
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &nranks);
      QDStageOperators &ops = *qd_ops_;
      const bool with_normal = ops.with_normal;
      const int nown = num_owned_fault_dofs_;
      Array<HYPRE_BigInt> part;
      HYPRE_BigInt fault_first = 0;
      const HYPRE_BigInt nf_glob = MakeQDPartition(comm, nown, part,
                                                   fault_first);
      const HYPRE_BigInt nslip_glob = 2 * nf_glob;
      const HYPRE_BigInt slip_first = 2 * fault_first;
      const double t_start = MPI_Wtime();

      // The loading must have one time shape, g(t) = t·ĝ: g(0) = 0 and
      // g(2) = 2·g(1).
      Vector g_hat;
      {
         EvaluateQDDirichletValues(0.0);
         real_t g0 = ops.g.Normlinf();
         EvaluateQDDirichletValues(1.0);
         g_hat = ops.g;
         EvaluateQDDirichletValues(2.0);
         real_t dev = 0.0, gmax = g_hat.Normlinf();
         for (int i = 0; i < g_hat.Size(); i++)
         {
            dev = std::max(dev, std::abs(ops.g(i) - 2.0 * g_hat(i)));
         }
         real_t loc[3] = {g0, dev, gmax}, glob[3];
         MPI_Allreduce(loc, glob, 3, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                       comm);
         MFEM_VERIFY(glob[0] == 0.0 && glob[1] <= 1e-12 * glob[2],
                     "--qd-green needs Dirichlet loading g(t) = t*g_hat: "
                     "max|g(0)| = " << glob[0] << ", max|g(2) - 2 g(1)| = "
                     << glob[1] << " (max|g(1)| = " << glob[2]
                     << ").  Run without --qd-green.");
      }

      // QD speed item 4 (--qd-green-compress): G as an H-matrix
      const bool compress = config.compress > 0.0;
      if (rank == 0)
      {
         const double gib = double(with_normal ? 3 : 2) * double(nf_glob) *
                            double(nslip_glob) * sizeof(real_t) / 1073741824.0;
         if (compress)
         {
            mfem::out << "[QD-GREEN] dense G would be " << std::fixed
                      << std::setprecision(2) << gib << " GiB; kept as an "
                      "H-matrix (dense rows only while compressing)\n";
         }
         else
         {
            mfem::out << "[QD-GREEN] allocating G: " << std::fixed
                      << std::setprecision(2) << gib << " GiB in all, about "
                      << gib / nranks << " GiB per rank (rows split over all "
                      << nranks << " ranks)\n";
         }
         mfem::out.flush();
      }
      MFEM_VERIFY(config.hcache_path.empty() || compress,
                  "--qd-green-hcache needs --qd-green-compress");
      const double mem_node = config.compress_mem_gib * 1073741824.0;
      auto green = std::unique_ptr<QDGreenOperator>(new QDGreenOperator(
         comm, nf_glob, fault_first, nown, with_normal, !compress));
      if (!config.cache_path.empty() || compress)
      {
         // Position key of each owned fault DOF (node, face centroid), so
         // the cache file can be read on another partition.
         MFEM_VERIFY(mesh_.SpaceDimension() == 3, "--qd-green-cache and "
                     "--qd-green-compress need a 3-D mesh (they key fault "
                     "DOFs by position)");
         Vector xyz, keys(6 * nown);
         GetFaultDOFCoords3D(xyz);
         const int nbf = nbf_per_face_;
         for (int o = 0; o < nown; o++)
         {
            const int ld = owned_fault_dof_to_local_dof_[o];
            const int f0 = (ld / nbf) * nbf;
            for (int d = 0; d < 3; d++)
            {
               real_t c = 0.0;
               for (int k = 0; k < nbf; k++) { c += xyz(3 * (f0 + k) + d); }
               keys(6 * o + d) = xyz(3 * ld + d);
               keys(6 * o + 3 + d) = c / nbf;
            }
         }
         green->SetKeys(keys);
      }
      std::unique_ptr<QDHMatrix> hmat;
      std::string hinfo;
      bool h_loaded = false;
      // QD speed item C: the bytes layout is planned from the H file's
      // block table, so it applies to a read (a compress run reads its
      // file back, below).
      MFEM_VERIFY(config.layout == "bytes" || config.layout == "count",
                  "--qd-green-layout: bytes or count, got '" << config.layout
                  << "'");
      const bool by_bytes = compress && config.layout == "bytes" &&
                            !config.hcache_path.empty();
      std::vector<real_t> keys_all;
      // Read the H file: header, the bytes layout if asked for, then data
      auto read_h = [&](QDHMatrix &h, std::string &hi)
      {
         if (!h.ReadHeader(config.hcache_path, config.compress, hi))
         {
            return false;
         }
         if (by_bytes) { green->SetRowLayout(h.PlanLayout()); }
         h.ReadData(config.hcache_path, mem_node, green->LoadingRowsWrite());
         return true;
      };
      if (compress)
      {
         green->GatherKeysAll(keys_all);
         hmat.reset(new QDHMatrix(comm, nf_glob, with_normal ? 3 : 2,
                                  green->RowLo(), keys_all, config.leaf,
                                  config.eta));
         if (!config.hcache_path.empty())
         {
            h_loaded = read_h(*hmat, hinfo);
         }
         if (!h_loaded) { green->AllocateDense(); }
      }
      std::string info;
      bool loaded = false;
      if (!h_loaded && !config.cache_path.empty())
      {
         loaded = green->ReadFile(config.cache_path, info);
      }
      double t_solve = 0.0, t_cols = 0.0;
      if (!loaded && !h_loaded)
      {
         const int m_max = config.block;
         std::vector<Vector> E(m_max), R(m_max), X(m_max), Ts(m_max),
             Tn(with_normal ? m_max : 0);
         for (int i = 0; i < m_max; i++)
         {
            E[i].SetSize(2 * nown);
            R[i].SetSize(ops.B->Height());
            X[i].SetSize(ops.B->Height());
         }
         for (HYPRE_BigInt j0 = 0; j0 < nslip_glob; j0 += m_max)
         {
            const int m = static_cast<int>(
               std::min<HYPRE_BigInt>(m_max, nslip_glob - j0));
            Array<const Vector *> Rp(m);
            Array<Vector *> Xp(m);
            for (int i = 0; i < m; i++)
            {
               E[i] = 0.0;
               const HYPRE_BigInt jl = j0 + i - slip_first;
               if (jl >= 0 && jl < 2 * nown) { E[i](static_cast<int>(jl)) = 1.0; }
               ops.B->Mult(E[i], R[i]);
               X[i] = 0.0;
               Rp[i] = &R[i];
               Xp[i] = &X[i];
            }
            const double t0 = MPI_Wtime();
            solver_->ArrayMult(Rp, Xp);
            const double t1 = MPI_Wtime();
            for (int i = 0; i < m; i++)
            {
               QDStageTraction(X[i], E[i], Ts[i],
                               with_normal ? &Tn[i] : nullptr);
            }
            green->SetColumns(j0, m, Ts, Tn);
            t_solve += t1 - t0;
            t_cols += MPI_Wtime() - t1;
         }
         // T_D: the loading ĝ with zero slip
         Vector zero(2 * nown), b(ops.B->Height()), x(ops.B->Height());
         zero = 0.0;
         ops.D->Mult(g_hat, b);
         x = 0.0;
         solver_->Mult(b, x);
         Vector TDs, TDn;
         QDStageTraction(x, zero, TDs, with_normal ? &TDn : nullptr);
         green->SetLoading(TDs, TDn);
      }
      if (h_loaded) { green->SetCompressed(std::move(hmat)); }
      MFEM_VERIFY(green->Complete(), "--qd-green: G is incomplete after "
                  << (loaded || h_loaded ? "reading" : "building") << " it");
      const double t_built = MPI_Wtime() - t_start;
      if (!loaded && !h_loaded && !config.cache_path.empty())
      {
         std::ostringstream o;
         o << "N_f " << nf_glob << ", order " << order_ << ", penalty factor "
           << penalty_factor_ << ", " << (with_normal ? 3 : 2)
           << " components, built on " << nranks << " ranks";
         green->WriteFile(config.cache_path, o.str());
      }
      std::string hsummary;
      if (compress && !h_loaded)
      {
         const bool keep = !config.hcache_path.empty();
         hmat->Compress(green->DenseRows(), config.compress, mem_node, keep);
         if (keep)
         {
            std::ostringstream o;
            o << "N_f " << nf_glob << ", order " << order_
              << ", penalty factor " << penalty_factor_ << ", "
              << (with_normal ? 3 : 2) << " components, compressed on "
              << nranks << " ranks from G "
              << (loaded ? "read from " + config.cache_path
                  : std::string("built by solves"));
            hmat->WriteFile(config.hcache_path, o.str(), green->LoadingRows());
         }
         hsummary = hmat->Summary("compressed");
         if (by_bytes)
         {
            // Item C: drop this operator and the dense rows, and read the
            // file just written back in the bytes layout (qd-c-design-rev1
            // §3.1); this also reads every block of it back whole.
            hmat.reset();
            green->FreeDense();
            hmat.reset(new QDHMatrix(comm, nf_glob, with_normal ? 3 : 2,
                                     green->RowLo(), keys_all, config.leaf,
                                     config.eta));
            MFEM_VERIFY(read_h(*hmat, hinfo), "--qd-green-hcache: cannot "
                        "reopen '" << config.hcache_path << "' after "
                        "writing it");
            hsummary += hmat->Summary("read back from " + config.hcache_path);
         }
         else if (config.layout == "bytes")
         {
            hsummary += "[QD-HMAT] layout count: --qd-green-layout bytes is "
                        "read from an H file, and this run writes none "
                        "(no --qd-green-hcache)\n";
         }
         green->SetCompressed(std::move(hmat));
      }
      else if (h_loaded)
      {
         hsummary = green->HMatrix()->Summary("read from " +
                                              config.hcache_path);
      }
      const double t_ready = MPI_Wtime() - t_start;
      qd_green_ = std::move(green);

      // Check against the solve path: slip of order 1 m varying along the
      // fault and t = 100 yr, so both G·δ and t·T_D matter.
      const real_t t_v = 100.0 * 365.25 * 86400.0;
      Vector dv(2 * nown);
      for (int i = 0; i < 2 * nown; i++)
      {
         const real_t j = static_cast<real_t>(slip_first + i);
         dv(i) = std::sin(0.37 * j + 0.1) + 0.5 * std::cos(0.011 * j);
      }
      Vector Tg, Tng, Ts_, Tns;
      qd_green_->Apply(t_v, dv, Tg, with_normal ? &Tng : nullptr, nullptr);
      QDSolvePathTraction(t_v, dv, Ts_, with_normal ? &Tns : nullptr);
      auto rel = [&](const Vector &a, const Vector &b)
      {
         real_t loc[2] = {0.0, 0.0}, glob[2];
         for (int i = 0; i < a.Size(); i++)
         {
            loc[0] = std::max(loc[0], std::abs(a(i) - b(i)));
            loc[1] = std::max(loc[1], std::abs(b(i)));
         }
         MPI_Allreduce(loc, glob, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                       comm);
         return glob[1] > 0.0 ? glob[0] / glob[1] : glob[0];
      };
      const real_t rel_s = rel(Tg, Ts_);
      const real_t rel_n = with_normal ? rel(Tng, Tns) : 0.0;

      // The build's g(t) calls are not stage evaluations.
      ops.g_evals = 0;
      ops.g_eval_seconds = 0.0;
      ops.g_copy_seconds = 0.0;

      double bytes = qd_green_->LocalBytes(), bmax = 0.0, btot = 0.0;
      MPI_Reduce(&bytes, &bmax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
      MPI_Reduce(&bytes, &btot, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
      double tm[4] = {t_built, t_solve, t_cols, t_ready}, tmax[4];
      MPI_Reduce(tm, tmax, 4, MPI_DOUBLE, MPI_MAX, 0, comm);
      if (rank == 0)
      {
         std::ostringstream o;
         o << "[QD-GREEN] G " << qd_green_->NComp() * nf_glob << " x "
           << nslip_glob << " (N_f " << nf_glob << ", "
           << (with_normal ? "shear + normal" : "shear") << " rows), "
           << std::fixed << std::setprecision(2) << btot / 1073741824.0
           << " GiB in all" << (compress ? " as an H-matrix" : "") << ", "
           << bmax / 1073741824.0 << " GiB largest rank\n";
         if (h_loaded)
         {
            o << "[QD-GREEN] H-matrix read from " << config.hcache_path
              << " (" << hinfo << ")\n";
         }
         else if (loaded)
         {
            o << "[QD-GREEN] read from " << config.cache_path << " in "
              << std::setprecision(1) << tmax[0] << " s (" << info << ")\n";
         }
         else
         {
            o << "[QD-GREEN] built in " << std::setprecision(1) << tmax[0]
              << " s (blocked solves " << tmax[1] << " s, traction and "
              << "redistribution " << tmax[2] << " s; " << config.block
              << " columns per solve)";
            if (!config.cache_path.empty())
            {
               o << ", written to " << config.cache_path;
            }
            o << "\n";
         }
         if (compress)
         {
            o << hsummary << "[QD-GREEN] compressed G ready "
              << std::fixed << std::setprecision(1) << tmax[3]
              << " s after the start of the build\n";
         }
         o << "[QD-GREEN] check against the solve path (t = 100 yr, test "
           << "slip): max rel. difference shear " << std::scientific
           << std::setprecision(3) << rel_s;
         if (with_normal) { o << ", normal " << rel_n; }
         o << " (tolerance " << config.tol << ")\n";
         mfem::out << o.str();
         mfem::out.flush();
      }
      MFEM_VERIFY(rel_s <= config.tol && rel_n <= config.tol,
                  "--qd-green: G differs from the solve path by " << rel_s
                  << " (shear), " << rel_n << " (normal), above --qd-green-tol "
                  << config.tol << (h_loaded ? ".  The H-matrix file "
                  "belongs to another model: remove it or pass another "
                  "--qd-green-hcache." : loaded ? ".  The cache file "
                  "belongs to another model: remove it or pass another "
                  "--qd-green-cache." : "."));
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GreenStage(
   real_t time, const Vector &slip_owned, Vector &traction_owned,
   Vector *normal_owned)
{
   MFEM_VERIFY(qd_green_, "--qd-green: BuildQDGreen() was not called");
   qd_green_->Apply(time, slip_owned, traction_owned, normal_owned,
                    stage_timers_);
}

#endif // MFEM_USE_MPI
