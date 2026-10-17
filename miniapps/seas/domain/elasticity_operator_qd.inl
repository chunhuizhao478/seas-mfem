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
      lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
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
         lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
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
      ops.g.Read();
      if (gpu) { MFEM_DEVICE_SYNC; }
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
   // Friction reads T and T_n on the host; on a device, copy them back.
   traction_owned.HostRead();
   if (normal_owned) { normal_owned->HostRead(); }
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

#endif // MFEM_USE_MPI
