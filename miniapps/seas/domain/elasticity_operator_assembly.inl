// elasticity_operator_assembly.inl
// Original range: Lines 3054-4896 of elasticity_operator.hpp (pre-decomposition)
// This file is included from elasticity_operator.hpp — do not include directly.

   // ========================================================================
   // Stiffness assembly
   // ========================================================================

   void AssembleStiffness() const
   {
      // In parallel, exchange face neighbor data before assembly.
      // ParBilinearForm::AssembleSharedFaces needs this to access
      // face neighbor elements and DOFs.
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         mesh_.ExchangeFaceNbrData();
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         if (pfes) { pfes->ExchangeFaceNbrData(); }
#endif
      }

      cached_a_ = std::make_unique<BilinFormType>(fes_.get());

      // v57: Enable full neighbor-block assembly for DG shared faces.
      // With keep_nbr_block=false (MFEM default), shared face contributions
      // only assemble Elem1 rows into a rectangular local matrix. This causes
      // the global HypreParMatrix to have partition-dependent NNZ and slight
      // numerical differences in the stiffness matrix, producing ~6%
      // displacement errors across different MPI rank counts.
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         cached_a_->KeepNbrBlock(true);
      }

      // Volume term: ∫_Ω σ(u):ε(v) dV
      cached_a_->AddDomainIntegrator(
         new ElasticityIntegrator(lambda_coeff_, mu_coeff_));

      if (method_ == DGMethod::BR2)
      {
         if (!mass_inv_computed_)
         {
            PrecomputeMassInverse();
         }
         // BR2 interior face integrator with elasticity tensor coupling
         cached_a_->AddInteriorFaceIntegrator(
            new DGElasticityBR2Integrator(
               lambda_coeff_, mu_coeff_, epsilon_, elem_mass_inv_, 3));

         // BR2 boundary face integrator for Dirichlet walls
         if (dirichlet_bdr_marker_.Size() > 0)
         {
            cached_a_->AddBdrFaceIntegrator(
               new DGElasticityBR2BoundaryIntegrator(
                  lambda_coeff_, mu_coeff_, epsilon_, elem_mass_inv_, 3),
               dirichlet_bdr_marker_);
         }
      }
      else  // IP
      {
         // v55: Combined integrator matching Tandem's assembleSurface kernel.
         // All three DG terms (consistency + symmetry + penalty) computed in
         // one pass with the same traction operator, same quadrature (2p+1),
         // and same penalty formula. K, fault slip b, Dirichlet b, and traction
         // recovery all use this same integrator.
         cached_a_->AddInteriorFaceIntegrator(
            new DGElasticityIPCombinedIntegrator(
               lambda_coeff_, mu_coeff_, 3, epsilon_, penalty_factor_));

         if (dirichlet_bdr_marker_.Size() > 0)
         {
            cached_a_->AddBdrFaceIntegrator(
               new DGElasticityIPCombinedIntegrator(
                  lambda_coeff_, mu_coeff_, 3, epsilon_, penalty_factor_),
               dirichlet_bdr_marker_);
         }
      }

      // v57: Assemble with skip_zeros=0 to ensure partition-independent NNZ.
      // With skip_zeros=1 (default), shared vs interior faces produce
      // slightly different face matrices (due to elem1/elem2 swap changing
      // floating-point order), causing O(ε) entries to be zero-skipped
      // on one partition but not another → different K matrix.
      cached_a_->Assemble(0);
      cached_a_->Finalize();

      // K assembly diagnostics: count faces/elements processed.
      // All ranks participate in MPI_Allreduce — MPI-safe.
      {
         long long local_volume = mesh_.GetNE();

         long long local_interior_face = 0;
         for (int f = 0; f < mesh_.GetNumFaces(); f++)
         {
            if (mesh_.GetInteriorFaceTransformations(f)) { local_interior_face++; }
         }

         // Count TRUE one-sided boundary faces (not interior) that K processes.
         // A boundary element is a true boundary face only if its underlying
         // mesh face does NOT have two local parent elements.
         // Use FaceIsInterior (checks Elem2No >= 0) instead of
         // GetBdrFaceTransformations (which can deadlock on ParMesh).
         long long local_bdr_face_true = 0;   // attr-5 true boundary faces
         long long local_bdr_all_true = 0;    // all-attr true boundary faces
         long long local_bdr_elem_total = mesh_.GetNBE();  // all boundary elements
         for (int be = 0; be < mesh_.GetNBE(); be++)
         {
            int face_idx, face_info_val;
            mesh_.GetBdrElementFace(be, &face_idx, &face_info_val);
            // Skip if the underlying face is interior (two local parents)
            // or shared (FaceIsTrueInterior includes shared faces)
            if (mesh_.FaceIsInterior(face_idx)) { continue; }
            local_bdr_all_true++;
            int attr = mesh_.GetBdrAttribute(be);
            if (dirichlet_bdr_marker_.Size() > 0 &&
                dirichlet_bdr_marker_[attr - 1] == 1)
            {
               local_bdr_face_true++;
            }
         }

         long long local_shared = 0;
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            local_shared = mesh_.GetNSharedFaces();
#endif
         }

         long long counts[5] = {local_volume, local_interior_face, local_shared,
                                local_bdr_face_true, local_bdr_all_true};
         long long global[5] = {0};

         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Allreduce(counts, global, 5, MPI_LONG_LONG, MPI_SUM,
                          mesh_.GetComm());
#endif
         }
         else
         {
            for (int i = 0; i < 5; i++) { global[i] = counts[i]; }
         }

         // Count true one-sided boundary faces per attribute.
         // Uses FaceIsTrueInterior to skip interior and shared faces
         // (same filter as the main boundary count above).
         int local_max_attr = mesh_.GetNBE() > 0 ? mesh_.bdr_attributes.Max() : 0;
         int max_attr = local_max_attr;
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Allreduce(&local_max_attr, &max_attr, 1, MPI_INT, MPI_MAX,
                          mesh_.GetComm());
#endif
         }
         std::vector<long long> local_bdr_by_attr(max_attr + 1, 0);
         for (int be = 0; be < mesh_.GetNBE(); be++)
         {
            int face_idx2, face_info_val2;
            mesh_.GetBdrElementFace(be, &face_idx2, &face_info_val2);
            if (mesh_.FaceIsInterior(face_idx2)) { continue; }
            int attr = mesh_.GetBdrAttribute(be);
            if (attr >= 1 && attr <= max_attr) { local_bdr_by_attr[attr]++; }
         }
         std::vector<long long> global_bdr_by_attr(max_attr + 1, 0);
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Allreduce(local_bdr_by_attr.data(), global_bdr_by_attr.data(),
                          max_attr + 1, MPI_LONG_LONG, MPI_SUM, mesh_.GetComm());
#endif
         }
         else
         {
            global_bdr_by_attr = local_bdr_by_attr;
         }

         if (DebugRank() == 0)
         {
            mfem::out << "  [K-DIAG] Volume elements:            " << global[0] << "\n";
            mfem::out << "  [K-DIAG] Interior faces (K):         " << global[1] << "\n";
            mfem::out << "  [K-DIAG] Shared faces (K):           " << global[2]
                      << " (each face counted by both ranks)\n";
            mfem::out << "  [K-DIAG] Bdr faces w/ K (Dirichlet): " << global[3]
                      << " (true one-sided, attr 5)\n";
            mfem::out << "  [K-DIAG] Bdr faces (all attrs):      " << global[4]
                      << " (true one-sided, any attr)\n";
            mfem::out << "  [K-DIAG] Bdr elements total:         " << local_bdr_elem_total
                      << " (local rank, includes interior)\n";
            for (int a = 1; a <= max_attr; a++)
            {
               if (global_bdr_by_attr[a] > 0)
               {
                  mfem::out << "  [K-DIAG]   attr " << a << ": "
                            << global_bdr_by_attr[a] << " bdr faces"
                            << (dirichlet_bdr_marker_.Size() >= a &&
                                dirichlet_bdr_marker_[a-1] == 1
                                ? " (Dirichlet K)" : " (NO K integrator)")
                            << "\n";
               }
            }
            mfem::out << "  [K-DIAG] Total K contributions:  "
                      << global[0] + global[1] + global[2] + global[3]
                      << " (vol + interior + shared + bdr_dir)\n";
         }
      }

      // Set up solver operator
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         cached_Ah_.SetType(Operator::Hypre_ParCSR);
         cached_a_->ParallelAssemble(cached_Ah_);

         // One-shot matrix analysis (enabled via --matrix-stats).  Runs
         // before solver factorization so it covers the same matrix the
         // solver sees, and is cheap (O(global_nnz) for one Transpose +
         // Add).
         if (matrix_stats_enabled_)
         {
            auto *K = cached_Ah_.As<HypreParMatrix>();
            const int dofs_per_elem = fes_->GetTypicalFE()->GetDof();
            const int vdim = fes_->GetVDim();
            MatrixStats stats =
               ComputeMatrixStats(*K, dofs_per_elem, vdim, mesh_.GetComm(),
                                  /*do_symmetry=*/true);
            int my_rank = 0;
            MPI_Comm_rank(mesh_.GetComm(), &my_rank);
            if (my_rank == 0)
            {
               PrintMatrixStatsSummary(stats);
               if (!matrix_stats_path_.empty())
               {
                  WriteMatrixStatsJson(stats, matrix_stats_path_);
                  mfem::out << "  [MATRIX-STATS] JSON written to "
                            << matrix_stats_path_ << "\n";
               }
            }
         }

         // K·v diagnostic: apply global K to all-ones vector,
         // dump result at target elements for cross-code comparison.
         if (first_step_debug_.enabled)
         {
            auto *K = cached_Ah_.As<HypreParMatrix>();
            int local_size = K->Height();
            HypreParVector ones(K->GetComm(), K->GetGlobalNumRows(),
                                K->GetRowStarts());
            ones = 1.0;
            HypreParVector Kv(K->GetComm(), K->GetGlobalNumRows(),
                              K->GetRowStarts());
            K->Mult(ones, Kv);

            // Dump K·1 for target elements
            if (DebugRank() == first_step_debug_.target_rank ||
                first_step_debug_.target_rank < 0)
            {
               std::string kv_path = DebugFilePath("first_step_Kv");
               std::ofstream out(kv_path, std::ios::trunc);
               out << "elem,component,local_dof,vdof,Kv_value\n";
               out << std::setprecision(17);

               const auto &target_elems = DebugTargetLocalElements();
               for (int e : target_elems)
               {
                  Array<int> vdofs;
                  fes_->GetElementVDofs(e, vdofs);
                  for (int j = 0; j < vdofs.Size(); j++)
                  {
                     int gj = vdofs[j];
                     int sign = 1;
                     if (gj < 0) { gj = -1 - gj; sign = -1; }
                     int comp = j / (vdofs.Size() / 3);
                     int ldof = j % (vdofs.Size() / 3);
                     double val = (gj < local_size) ? Kv(gj) * sign : 0.0;
                     out << e << "," << comp << "," << ldof << ","
                         << vdofs[j] << "," << val << "\n";
                  }
               }
               mfem::out << "  [K-DIAG] K·1 dumped to " << kv_path
                         << " (" << target_elems.size() << " elements)\n";
            }

            // Print global norms (collective, rank 0 prints)
            double local_norm1 = 0.0, local_norm2sq = 0.0, local_norminf = 0.0;
            for (int i = 0; i < local_size; i++)
            {
               double v = std::abs(Kv(i));
               local_norm1 += v;
               local_norm2sq += v * v;
               if (v > local_norminf) { local_norminf = v; }
            }
            double global_norm1, global_norm2sq, global_norminf;
            MPI_Reduce(&local_norm1, &global_norm1, 1, MPI_DOUBLE, MPI_SUM,
                        0, mesh_.GetComm());
            MPI_Reduce(&local_norm2sq, &global_norm2sq, 1, MPI_DOUBLE, MPI_SUM,
                        0, mesh_.GetComm());
            MPI_Reduce(&local_norminf, &global_norminf, 1, MPI_DOUBLE, MPI_MAX,
                        0, mesh_.GetComm());
            int myrank;
            MPI_Comm_rank(mesh_.GetComm(), &myrank);
            if (myrank == 0)
            {
               mfem::out << std::setprecision(15);
               mfem::out << "  [K-DIAG] ||K·1||_1   = " << global_norm1 << "\n";
               mfem::out << "  [K-DIAG] ||K·1||_2   = "
                         << std::sqrt(global_norm2sq) << "\n";
               mfem::out << "  [K-DIAG] ||K·1||_inf = " << global_norminf << "\n";
            }
         }

#ifdef MFEM_USE_MUMPS
         const int mumps_print_level = matrix_stats_enabled_ ? 2 : 1;
         if (solver_type_ == SolverType::MUMPS)
         {
            auto *mumps = new MUMPSSolver(mesh_.GetComm());
            mumps->SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
            mumps->SetPrintLevel(mumps_print_level);
            mumps->SetOperator(*cached_Ah_.As<HypreParMatrix>());
            solver_.reset(mumps);
         }
         else if (solver_type_ == SolverType::MUMPS_BLR)
         {
            auto *mumps = new MUMPSSolver(mesh_.GetComm());
            mumps->SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
            mumps->SetPrintLevel(mumps_print_level);
            mumps->SetBLRTol(blr_tol_);
            mumps->SetOperator(*cached_Ah_.As<HypreParMatrix>());
            solver_.reset(mumps);
         }
         else
#endif
#ifdef MFEM_USE_SUPERLU
         if (solver_type_ == SolverType::SUPERLU)
         {
            cached_superlu_mat_.reset(
               new SuperLURowLocMatrix(*cached_Ah_.As<HypreParMatrix>()));
            auto *superlu = new SuperLUSolver(mesh_.GetComm());
            superlu->SetOperator(*cached_superlu_mat_);
            superlu->SetColumnPermutation(superlu::PARMETIS);
            superlu->SetIterativeRefine(superlu::SLU_DOUBLE);
            superlu->SetSymmetricPattern(true);
            superlu->SetPrintStatistics(true);
            solver_.reset(superlu);
         }
         else
#endif
#ifdef MFEM_USE_STRUMPACK
         if (solver_type_ == SolverType::STRUMPACK)
         {
            cached_strumpack_mat_.reset(
               new STRUMPACKRowLocMatrix(*cached_Ah_.As<HypreParMatrix>()));
            auto *strumpack = new STRUMPACKSolver(mesh_.GetComm());
            strumpack->SetOperator(*cached_strumpack_mat_);
            strumpack->SetCompression(strumpack::CompressionType::BLR);
            strumpack->SetCompressionRelTol(1e-10);
            strumpack->SetReorderingStrategy(strumpack::ReorderingStrategy::METIS);
            strumpack->SetKrylovSolver(strumpack::KrylovSolver::DIRECT);
            strumpack->SetPrintFactorStatistics(true);
            strumpack->SetPrintSolveStatistics(true);
            solver_.reset(strumpack);
         }
         else
#endif
         if (solver_type_ == SolverType::GMRES_BlockILU)
         {
            auto *gmres = new GMRESSolver(mesh_.GetComm());
            gmres->SetRelTol(1e-10);
            gmres->SetAbsTol(0.0);
            gmres->SetMaxIter(1000);
            gmres->SetKDim(50);
            gmres->SetPrintLevel(0);

            // BlockILU with element-sized blocks — natural for DG
            int block_size = fes_->GetTypicalFE()->GetDof() * 3;  // vdim=3
            cached_prec_.reset(new BlockILU(block_size,
               BlockILU::Reordering::MINIMUM_DISCARDED_FILL));

            gmres->SetPreconditioner(*cached_prec_);
            gmres->SetOperator(*cached_Ah_.As<HypreParMatrix>());
            solver_.reset(gmres);
         }
         else if (solver_type_ == SolverType::GMRES_AMG)
         {
            // GMRES + BoomerAMG: for DG at high p where:
            //   - MUMPS runs out of memory (dense DG factorization)
            //   - GMRES+BlockILU doesn't converge (local preconditioner too weak)
            //   - CG+AMG fails because AMG is non-SPD for DG matrices
            // GMRES tolerates non-SPD preconditioners while still benefiting
            // from AMG's multilevel global coarse-grid correction.
            // Phase 4: built ONCE here (guarded by stiffness_assembled_) and
            // reused for every solver_->Mult; tolerances from DomainConfig.
            MFEM_PERF_SCOPE("seas::elasticity::SetupSolver::GMRES_AMG");
            auto *gmres = new GMRESSolver(mesh_.GetComm());
            gmres->SetRelTol(ksp_rtol_);
            gmres->SetAbsTol(ksp_atol_);
            gmres->SetMaxIter(ksp_maxit_);
            gmres->SetKDim(100);
            gmres->SetPrintLevel(amg_print_level_);

            auto *amg = new HypreBoomerAMG(*cached_Ah_.As<HypreParMatrix>());
            // DG vector space is Ordering::byNODES (elasticity_operator_setup.inl);
            // map the 3 vector components to AMG DOF-functions accordingly.
            MFEM_VERIFY(fes_->GetOrdering() == Ordering::byNODES,
                        "GMRES_AMG: SetSystemsOptions assumes byNODES ordering "
                        "(elasticity DG space is byNODES)");
            amg->SetSystemsOptions(3, /*order_bynodes=*/true);
            // Skip SetElasticityOptions for the non-CG path — DG sparsity breaks
            // the CFEM rigid-body near-null-space assumption (R-006).
            // Same DG perf tuning as CG_AMG (parallel smoother + aggressive
            // coarsening; see the CG_AMG path).
            if (amg_relax_type_ >= 0) { amg->SetRelaxType(amg_relax_type_); }
            if (amg_aggressive_levels_ > 0)
            {
               HYPRE_BoomerAMGSetAggNumLevels(
                  static_cast<HYPRE_Solver>(*amg), amg_aggressive_levels_);
            }
            amg->SetPrintLevel(amg_print_level_);
            cached_prec_.reset(amg);

            gmres->SetPreconditioner(*cached_prec_);
            gmres->SetOperator(*cached_Ah_.As<HypreParMatrix>());
            solver_.reset(gmres);
         }
         else
         {
            // CG + BoomerAMG (SPD).  Phase 4: built ONCE here, reused per solve;
            // tolerances from DomainConfig.
            MFEM_PERF_SCOPE("seas::elasticity::SetupSolver::CG_AMG");
            auto *cg = new CGSolver(mesh_.GetComm());
            cg->SetRelTol(ksp_rtol_);
            cg->SetAbsTol(ksp_atol_);
            cg->SetMaxIter(ksp_maxit_);
            cg->SetPrintLevel(amg_print_level_);

            auto *amg = new HypreBoomerAMG(*cached_Ah_.As<HypreParMatrix>());
            MFEM_VERIFY(fes_->GetOrdering() == Ordering::byNODES,
                        "CG_AMG: SetSystemsOptions assumes byNODES ordering "
                        "(elasticity DG space is byNODES)");
            amg->SetSystemsOptions(3, /*order_bynodes=*/true);
            // R-006: rigid-body near-null-space modes for the SPD CG path.
            // Unproven for DG (MFEM ex17p does not use it) and bloats the AMG
            // hierarchy there (slow V-cycles); gated by amg_elasticity_options_
            // ([solver].amg_elasticity_options, default true for back-compat).
            // Set false on DG to use SetSystemsOptions only (the GMRES_AMG path,
            // which the iterative-vs-direct test verifies converges).  The
            // mandatory post-solve residual check guards a non-converged solve.
            if (amg_elasticity_options_)
            {
               amg->SetElasticityOptions(
                  dynamic_cast<ParFiniteElementSpace*>(fes_.get()));
            }
            // PERF (DG): the default hierarchy has operator complexity ~4.5 ->
            // expensive V-cycles, and the default hybrid-GS smoother (relax type
            // 3) couples sequentially across ranks -> np>1 stalls.  Config knobs
            // [solver].amg_relax_type (default 8 = l1-symmetric-GS, fully
            // parallel + SPD-robust) and amg_aggressive_levels (default 1 =
            // distance-2 coarsening on the top level -> complexity ~4.5->2.1,
            // ~2x cheaper V-cycles).  Measured ~10x faster init on the BP5 mesh;
            // iterative-vs-direct verifies it still converges to 1e-7.
            if (amg_relax_type_ >= 0) { amg->SetRelaxType(amg_relax_type_); }
            if (amg_aggressive_levels_ > 0)
            {
               HYPRE_BoomerAMGSetAggNumLevels(
                  static_cast<HYPRE_Solver>(*amg), amg_aggressive_levels_);
            }
            amg->SetPrintLevel(amg_print_level_);
            cached_prec_.reset(amg);

            cg->SetPreconditioner(*cached_prec_);
            cg->SetOperator(*cached_Ah_.As<HypreParMatrix>());
            solver_.reset(cg);
         }
#endif
      }
      else
      {
         SparseMatrix &A = cached_a_->SpMat();
         serial_prec_ = std::make_unique<GSSmoother>(A);
         static_cast<CGSolver*>(solver_.get())->SetPreconditioner(*serial_prec_);
         solver_->SetOperator(A);
      }

      stiffness_assembled_ = true;
      ++num_stiffness_assemblies_;   // Phase 4 (R-402): reuse counter
   }

   // ========================================================================
   // Slip BC assembly — IP method
   // ========================================================================

   void AssembleSlipContributionIP(Vector &rhs, const Vector &slip_bc) const
   {
      int dim = 3;
      int nbf = nbf_per_face_;

      // v55: Use the combined integrator's AssembleSlipFaceRHS to guarantee
      // K-b consistency. Same traction operator, quadrature, and penalty as K.
      DGElasticityIPCombinedIntegrator slip_integrator(
         lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         // Check if any slip is non-zero
         bool all_zero = true;
         for (int kk = 0; kk < nbf && all_zero; kk++)
         {
            int dof_idx = fi * nbf + kk;
            if (std::abs(slip_bc(2 * dof_idx)) >= 1e-15 ||
                std::abs(slip_bc(2 * dof_idx + 1)) >= 1e-15)
            { all_zero = false; }
         }
         if (all_zero) { continue; }

         // Build slip at quad points (Tandem evaluate_slip)
         Vector delta_u_quad;
         BuildSlipAtQuadPoints(fi, slip_bc, delta_u_quad);

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);

         // Assemble slip RHS using the combined integrator
         Vector elvec1, elvec2;
         slip_integrator.AssembleSlipFaceRHS(
            *fe1, *fe2, *FTr, delta_u_quad, elvec1, elvec2);
         DebugDumpFaceData(fi, face, "fault_interior", FTr, nullptr,
                           &delta_u_quad, elvec1, &elvec2);

         // Scatter into global RHS
         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);

         for (int j = 0; j < vdofs1.Size(); j++)
         {
            int gj = vdofs1[j];
            if (gj >= 0) { rhs(gj) += elvec1(j); }
            else { rhs(-1 - gj) -= elvec1(j); }
         }
         for (int j = 0; j < vdofs2.Size(); j++)
         {
            int gj = vdofs2[j];
            if (gj >= 0) { rhs(gj) += elvec2(j); }
            else { rhs(-1 - gj) -= elvec2(j); }
         }
      }
   }

   // ========================================================================
   // Slip BC assembly — BR2 method
   // ========================================================================

   void AssembleSlipContributionBR2(Vector &rhs, const Vector &slip_bc) const
   {
      int dim = 3;
      if (!mass_inv_computed_) { PrecomputeMassInverse(); }

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         if (FTr == nullptr) { continue; }

         // BR2 penalty: num faces per element (geometry-dependent)
         Geometry::Type geom = mesh_.GetElementGeometry(FTr->Elem1No);
         real_t penalty = (geom == Geometry::TETRAHEDRON)
                              ? real_t(dim + 1) : real_t(2 * dim);

         // Embed local slip to global frame
         real_t slip_local[2] = {slip_bc(2 * fi), slip_bc(2 * fi + 1)};
         real_t delta_u[3];
         fault_basis_.EmbedSlip(fi, slip_local, delta_u);

         real_t max_delta = std::max({std::abs(delta_u[0]),
                                       std::abs(delta_u[1]),
                                       std::abs(delta_u[2])});
         if (max_delta < 1e-15) { continue; }

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
         int ndof1 = fe1->GetDof();
         int ndof2 = fe2->GetDof();

         const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
         const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];

         int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
         int quad_order_br2 = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, quad_order_br2);
         int nqp = ir.GetNPoints();

         // Precompute shapes and normals
         DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
         Vector nor_all(dim * nqp);
         Vector w_all(nqp);
         Vector slip_sign_all(nqp);

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            FTr->SetAllIntPoints(&ip);
            const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

            Vector s1q(shapes1.GetColumn(q), ndof1);
            Vector s2q(shapes2.GetColumn(q), ndof2);
            fe1->CalcShape(eip1, s1q);
            fe2->CalcShape(eip2, s2q);

            Vector nor_q(&nor_all[q * dim], dim);
            CalcOrtho(FTr->Jacobian(), nor_q);

            w_all[q] = ip.weight;
            slip_sign_all[q] = (nor_q(1) > 0) ? 1.0 : -1.0;
         }

         // ===================================================================
         // Compute BR2 lifted slip with elasticity tensor coupling
         //
         // Step 1: face_int[elem][m,s] = sum_q E[elem][m,q] * delta_u_signed[s,q] * w[q]
         //   where delta_u_signed[s,q] = sum_c slip_sign * delta_u[c] * ... no, simpler:
         //   Actually for the RHS, the "source" is the imposed slip delta_u.
         //   The lifting formula is:
         //     f_lifted[elem][u,m,s] = 0.5 * Minv[m,o] * sum_q E[o,q] * delta_u_signed[u] * n[s,q] * w[q]
         //   Then f_lifted_q[i,q] = 0.5 * sum_{elem,s,m} test_normal_{iu,sq} * E[elem][m,q] * f_lifted[elem][u,m,s]
         //
         // This follows exactly the scalar BR2 pattern from antiplane but with
         // delta_u being 3-component and test_normal coupling.
         // ===================================================================

         // Compute face integrals for each element
         // face_int[elem][u*dim+s, m] = sum_q E[m,q] * sign_q * delta_u[u] * n[s,q] * w[q]
         DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
         face_int1 = 0.0;
         face_int2 = 0.0;

         for (int q = 0; q < nqp; q++)
         {
            real_t wq = w_all[q];
            real_t sign_q = slip_sign_all[q];

            for (int u = 0; u < dim; u++)
            {
               real_t delta_u_signed = sign_q * delta_u[u];
               for (int s = 0; s < dim; s++)
               {
                  real_t n_s = nor_all[q * dim + s];
                  real_t factor = delta_u_signed * n_s * wq;

                  for (int m = 0; m < ndof1; m++)
                  {
                     face_int1(u * dim + s, m) += shapes1(m, q) * factor;
                  }
                  for (int m = 0; m < ndof2; m++)
                  {
                     face_int2(u * dim + s, m) += shapes2(m, q) * factor;
                  }
               }
            }
         }

         // f_lifted[elem][u,s,m] = 0.5 * Minv[m,o] * face_int[elem][u*dim+s, o]
         DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
         MultABt(face_int1, Minv1, f_lifted1);
         f_lifted1 *= 0.5;
         MultABt(face_int2, Minv2, f_lifted2);
         f_lifted2 *= 0.5;

         // f_lifted_q[i,q] = 0.5 * sum_{elem,u,s,m} test_normal_{iu,s} * E[elem][m,q] * f_lifted[elem][u,s,m]
         DenseMatrix f_lifted_q(dim, nqp);
         f_lifted_q = 0.0;

         for (int q = 0; q < nqp; q++)
         {
            Vector n_q(dim);
            for (int d = 0; d < dim; d++) { n_q(d) = nor_all[q * dim + d]; }

            for (int i = 0; i < dim; i++)
            {
               real_t sum = 0.0;
               for (int u = 0; u < dim; u++)
               {
                  for (int s = 0; s < dim; s++)
                  {
                     real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * n_q(i)
                        + mu_val_ * ((i == u ? 1.0 : 0.0) * n_q(s)
                                     + (i == s ? 1.0 : 0.0) * n_q(u));

                     // Eval lifted: sum_m E[elem][m,q] * f_lifted[elem][u*dim+s, m]
                     real_t eval1 = 0.0, eval2 = 0.0;
                     for (int m = 0; m < ndof1; m++)
                     {
                        eval1 += shapes1(m, q) * f_lifted1(u * dim + s, m);
                     }
                     for (int m = 0; m < ndof2; m++)
                     {
                        eval2 += shapes2(m, q) * f_lifted2(u * dim + s, m);
                     }
                     sum += tn * (eval1 + eval2);
                  }
               }
               f_lifted_q(i, q) = 0.5 * sum;
            }
         }

         // ===================================================================
         // Assemble RHS:
         // b_k^i += c1 * sum_u [σ_i(φ_k e_i)·n]_u * delta_u_signed[u] * w_elem
         //        ± penalty * φ_k * f_lifted_q[i,q] * w
         // ===================================================================

         Vector elvec1(vdofs1.Size()), elvec2(vdofs2.Size());
         elvec1 = 0.0;
         elvec2 = 0.0;

         real_t c1 = epsilon_ * 0.5;  // symmetry data: −∫ g^F · {{σ(v)·n}} ds

         for (int q = 0; q < nqp; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            FTr->SetAllIntPoints(&ip);
            const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
            const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

            real_t wq = w_all[q];
            real_t sign_q = slip_sign_all[q];

            Vector nor_q(dim);
            for (int d = 0; d < dim; d++) { nor_q[d] = nor_all[q * dim + d]; }

            // Gradients
            DenseMatrix dshape1_ref(ndof1, dim), dshape2_ref(ndof2, dim);
            fe1->CalcDShape(eip1, dshape1_ref);
            fe2->CalcDShape(eip2, dshape2_ref);

            DenseMatrix adjJ1(dim), adjJ2(dim);
            CalcAdjugate(FTr->Elem1->Jacobian(), adjJ1);
            CalcAdjugate(FTr->Elem2->Jacobian(), adjJ2);

            DenseMatrix dshape1_adj(ndof1, dim), dshape2_adj(ndof2, dim);
            Mult(dshape1_ref, adjJ1, dshape1_adj);
            Mult(dshape2_ref, adjJ2, dshape2_adj);

            real_t detJ1 = FTr->Elem1->Weight();
            real_t detJ2 = FTr->Elem2->Weight();

            // Element 1 symmetry data + BR2 lifting data
            for (int k = 0; k < ndof1; k++)
            {
               real_t grad_dot_n = 0.0;
               for (int d = 0; d < dim; d++)
               {
                  grad_dot_n += dshape1_adj(k, d) * nor_q(d);
               }

               for (int i = 0; i < dim; i++)
               {
                  real_t sym_val = 0.0;
                  for (int u = 0; u < dim; u++)
                  {
                     real_t trac_iu = lambda_val_ * dshape1_adj(k, i) * nor_q(u)
                        + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                     + dshape1_adj(k, u) * nor_q(i));
                     sym_val += trac_iu * sign_q * delta_u[u];
                  }

                  int idx = i * ndof1 + k;  // byNODES ordering
                  elvec1(idx) += c1 * wq * sym_val / detJ1;
                  elvec1(idx) += penalty * wq * shapes1(k, q) * f_lifted_q(i, q);
               }
            }

            // Element 2 symmetry data + BR2 lifting data
            for (int k = 0; k < ndof2; k++)
            {
               real_t grad_dot_n = 0.0;
               for (int d = 0; d < dim; d++)
               {
                  grad_dot_n += dshape2_adj(k, d) * nor_q(d);
               }

               for (int i = 0; i < dim; i++)
               {
                  real_t sym_val = 0.0;
                  for (int u = 0; u < dim; u++)
                  {
                     real_t trac_iu = lambda_val_ * dshape2_adj(k, i) * nor_q(u)
                        + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                     + dshape2_adj(k, u) * nor_q(i));
                     sym_val += trac_iu * sign_q * delta_u[u];
                  }

                  int idx = i * ndof2 + k;
                  elvec2(idx) += c1 * wq * sym_val / detJ2;
                  elvec2(idx) += (-penalty) * wq * shapes2(k, q) * f_lifted_q(i, q);
               }
            }
         }

         // Add to global RHS
         for (int j = 0; j < vdofs1.Size(); j++)
         {
            int gj = vdofs1[j];
            if (gj >= 0) { rhs(gj) += elvec1(j); }
            else { rhs(-1 - gj) -= elvec1(j); }
         }
         for (int j = 0; j < vdofs2.Size(); j++)
         {
            int gj = vdofs2[j];
            if (gj >= 0) { rhs(gj) += elvec2(j); }
            else { rhs(-1 - gj) -= elvec2(j); }
         }
      }
   }

   // ========================================================================
   // Shared Face Slip Assembly — IP method (Parallel)
   // ========================================================================

   void AssembleSlipContributionIPShared(Vector &rhs, const Vector &slip_bc,
                                         int interior_face_count) const
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int dim = 3;
         int nbf = nbf_per_face_;

         // v55: Use combined integrator for K-b consistency (shared faces)
         DGElasticityIPCombinedIntegrator slip_integrator(
            lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            int slip_idx = interior_face_count + i;
            bool all_zero = true;
            for (int kk = 0; kk < nbf && all_zero; kk++)
            {
               int dof_idx = slip_idx * nbf + kk;
               if (std::abs(slip_bc(2 * dof_idx)) >= 1e-15 ||
                   std::abs(slip_bc(2 * dof_idx + 1)) >= 1e-15)
               { all_zero = false; }
            }
            if (all_zero) { continue; }

            Vector delta_u_quad;
            BuildSlipAtQuadPoints(slip_idx, slip_bc, delta_u_quad);

            const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
            auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
            int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);

            // Assemble using combined integrator — only use elvec1 (local elem)
            Vector elvec1, elvec2;
            slip_integrator.AssembleSlipFaceRHS(
               *fe1, *fe2, *FTr, delta_u_quad, elvec1, elvec2);
            DebugDumpFaceData(slip_idx, mesh_.GetSharedFace(sf), "fault_shared",
                              FTr, nullptr, &delta_u_quad, elvec1, &elvec2);

            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
            for (int j = 0; j < vdofs1.Size(); j++)
            {
               int gj = vdofs1[j];
               if (gj >= 0) { rhs(gj) += elvec1(j); }
               else { rhs(-1 - gj) -= elvec1(j); }
            }
         }
#endif
      }
   }

   // ========================================================================
   // Shared Face Slip Assembly — BR2 method (Parallel)
   // ========================================================================

   void AssembleSlipContributionBR2Shared(Vector &rhs, const Vector &slip_bc,
                                           int interior_face_count) const
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int dim = 3;
         if (!mass_inv_computed_) { PrecomputeMassInverse(); }

         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());

         for (int fi = 0; fi < fault_shared_faces_.Size(); fi++)
         {
            int sf = fault_shared_faces_[fi];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            int slip_idx = interior_face_count + fi;
            real_t slip_local[2] = {slip_bc(2 * slip_idx),
                                    slip_bc(2 * slip_idx + 1)};
            real_t delta_u[3];
            fault_basis_.EmbedSlip(slip_idx, slip_local, delta_u);

            real_t max_delta = std::max({std::abs(delta_u[0]),
                                          std::abs(delta_u[1]),
                                          std::abs(delta_u[2])});
            if (max_delta < 1e-15) { continue; }

            // BR2 penalty
            Geometry::Type geom = mesh_.GetElementGeometry(FTr->Elem1No);
            real_t penalty = (geom == Geometry::TETRAHEDRON)
                                 ? real_t(dim + 1) : real_t(2 * dim);

            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);

            const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
            int ndof1 = fe1->GetDof();

            int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
            int ndof2 = fe2->GetDof();

            const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
            const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];

            int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
            int quad_order_br2_sh = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
            const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, quad_order_br2_sh);
            int nqp = ir.GetNPoints();

            // Precompute shapes and normals
            DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
            Vector nor_all(dim * nqp);
            Vector w_all(nqp);
            Vector slip_sign_all(nqp);

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();
               const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

               Vector s1q(shapes1.GetColumn(q), ndof1);
               Vector s2q(shapes2.GetColumn(q), ndof2);
               fe1->CalcShape(eip1, s1q);
               fe2->CalcShape(eip2, s2q);

               Vector nor_q(&nor_all[q * dim], dim);
               CalcOrtho(FTr->Jacobian(), nor_q);

               w_all[q] = ip.weight;
               slip_sign_all[q] = (nor_q(1) > 0) ? 1.0 : -1.0;
            }

            // Two-sided BR2 lifting (both elements contribute to lifting)
            DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
            face_int1 = 0.0;
            face_int2 = 0.0;

            for (int q = 0; q < nqp; q++)
            {
               real_t wq = w_all[q];
               real_t sign_q = slip_sign_all[q];

               for (int u = 0; u < dim; u++)
               {
                  real_t delta_u_signed = sign_q * delta_u[u];
                  for (int s = 0; s < dim; s++)
                  {
                     real_t n_s = nor_all[q * dim + s];
                     real_t factor = delta_u_signed * n_s * wq;

                     for (int m = 0; m < ndof1; m++)
                     {
                        face_int1(u * dim + s, m) += shapes1(m, q) * factor;
                     }
                     for (int m = 0; m < ndof2; m++)
                     {
                        face_int2(u * dim + s, m) += shapes2(m, q) * factor;
                     }
                  }
               }
            }

            DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
            MultABt(face_int1, Minv1, f_lifted1);
            f_lifted1 *= 0.5;
            MultABt(face_int2, Minv2, f_lifted2);
            f_lifted2 *= 0.5;

            // f_lifted_q[i,q] with elasticity tensor coupling
            DenseMatrix f_lifted_q(dim, nqp);
            f_lifted_q = 0.0;

            for (int q = 0; q < nqp; q++)
            {
               Vector n_q(dim);
               for (int d = 0; d < dim; d++) { n_q(d) = nor_all[q * dim + d]; }

               for (int ci = 0; ci < dim; ci++)
               {
                  real_t sum = 0.0;
                  for (int u = 0; u < dim; u++)
                  {
                     for (int s = 0; s < dim; s++)
                     {
                        real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * n_q(ci)
                           + mu_val_ * ((ci == u ? 1.0 : 0.0) * n_q(s)
                                        + (ci == s ? 1.0 : 0.0) * n_q(u));

                        real_t eval1 = 0.0, eval2 = 0.0;
                        for (int m = 0; m < ndof1; m++)
                        {
                           eval1 += shapes1(m, q) * f_lifted1(u * dim + s, m);
                        }
                        for (int m = 0; m < ndof2; m++)
                        {
                           eval2 += shapes2(m, q) * f_lifted2(u * dim + s, m);
                        }
                        sum += tn * (eval1 + eval2);
                     }
                  }
                  f_lifted_q(ci, q) = 0.5 * sum;
               }
            }

            // Assemble RHS (Elem1 only — local element)
            Vector elvec1(vdofs1.Size());
            elvec1 = 0.0;

            real_t c1 = epsilon_ * 0.5;  // symmetry data: −∫ g^F · {{σ(v)·n}} ds

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();

               real_t wq = w_all[q];
               real_t sign_q = slip_sign_all[q];

               Vector nor_q(dim);
               for (int d = 0; d < dim; d++) { nor_q[d] = nor_all[q * dim + d]; }

               DenseMatrix dshape1_ref(ndof1, dim);
               fe1->CalcDShape(eip1, dshape1_ref);

               DenseMatrix adjJ1(dim);
               CalcAdjugate(FTr->Elem1->Jacobian(), adjJ1);

               DenseMatrix dshape1_adj(ndof1, dim);
               Mult(dshape1_ref, adjJ1, dshape1_adj);

               real_t detJ1 = FTr->Elem1->Weight();

               for (int k = 0; k < ndof1; k++)
               {
                  real_t grad_dot_n = 0.0;
                  for (int d = 0; d < dim; d++)
                  {
                     grad_dot_n += dshape1_adj(k, d) * nor_q(d);
                  }

                  for (int ci = 0; ci < dim; ci++)
                  {
                     real_t sym_val = 0.0;
                     for (int u = 0; u < dim; u++)
                     {
                        real_t trac_iu = lambda_val_ * dshape1_adj(k, ci) * nor_q(u)
                           + mu_val_ * ((ci == u ? 1.0 : 0.0) * grad_dot_n
                                        + dshape1_adj(k, u) * nor_q(ci));
                        sym_val += trac_iu * sign_q * delta_u[u];
                     }

                     int idx = ci * ndof1 + k;
                     elvec1(idx) += c1 * wq * sym_val / detJ1;
                     elvec1(idx) += penalty * wq * shapes1(k, q) * f_lifted_q(ci, q);
                  }
               }
            }

            for (int j = 0; j < vdofs1.Size(); j++)
            {
               int gj = vdofs1[j];
               if (gj >= 0) { rhs(gj) += elvec1(j); }
               else { rhs(-1 - gj) -= elvec1(j); }
            }
         }
#endif
      }
   }

   // ========================================================================
   // Dirichlet loading assembly
   // ========================================================================

   void AssembleDirichletLoading(Vector &rhs, real_t time) const
   {
      // BP5 Dirichlet loading: u_D = (sgn(y)*Vp*t/2, 0, 0) on far-field,
      //                        u_D = (Vp*t, 0, 0) on y=0 skeleton faces.
      // Applied on attr-5 faces only (controlled by BCMode).
      //
      // DG Dirichlet BC contribution:
      //   b[k,i] += c0 * [σ(φ_k e_i)·n]_u * u_D_u * (1/detJ)
      //           + penalty * φ_k * u_D_i

      if (std::abs(time * Vp_) < 1e-30) { return; }

      int dim = 3;

      // Build set of face indices already handled as interior/shared Dirichlet.
      // These must NOT also be processed as boundary faces (would double-load
      // with incompatible DG formulas: boundary uses c1=epsilon, skeleton uses
      // c1=0.5*epsilon).
      std::set<int> dir_interior_set;
      for (int fi = 0; fi < dirichlet_interior_faces_.Size(); fi++)
      {
         dir_interior_set.insert(dirichlet_interior_faces_[fi]);
      }
      std::set<int> dir_shared_set;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int fi = 0; fi < dirichlet_shared_faces_.Size(); fi++)
         {
            dir_shared_set.insert(mesh_.GetSharedFace(dirichlet_shared_faces_[fi]));
         }
#endif
      }

      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (dirichlet_bdr_marker_[attr - 1] != 1) { continue; }

         // Get face transformation
         int face_idx, face_info;
         mesh_.GetBdrElementFace(be, &face_idx, &face_info);

         // Skip faces already handled as interior or shared Dirichlet
         if (dir_interior_set.count(face_idx) > 0) { continue; }
         if (dir_shared_set.count(face_idx) > 0) { continue; }

         FaceElementTransformations *FTr =
            mesh_.GetFaceElementTransformations(face_idx);
         if (FTr == nullptr) { continue; }

         Array<int> vdofs;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs);
         const FiniteElement *fe = scalar_fes_->GetFE(FTr->Elem1No);
         int ndof = fe->GetDof();

         int face_order = fe->GetOrder();
         int quad_order_bdr = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, quad_order_bdr);

         Vector elvec(vdofs.Size());
         elvec = 0.0;

         if (method_ == DGMethod::IP)
         {
            // v55: Use combined integrator's AssembleBoundaryFaceRHS
            // (same formula as K's boundary face, Tandem's rhs_boundary)
            DGElasticityIPCombinedIntegrator dir_integ(
               lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

            int quad_order_dir = 2 * fe->GetOrder() + 1;
            const IntegrationRule &ir_dir = IntRules.Get(
               FTr->FaceGeom, quad_order_dir);
            int nq_dir = ir_dir.GetNPoints();

            // v55: Evaluate Tandem's boundary(x,y,z,t) at each quad point
            Vector u_D_3d(dim * nq_dir);
            u_D_3d = 0.0;
            Vector phys_y_qp(nq_dir);
            Vector u_D_val(dim);
            for (int q = 0; q < nq_dir; q++)
            {
               const IntegrationPoint &ipq = ir_dir.IntPoint(q);
               FTr->SetAllIntPoints(&ipq);
               Vector phys(dim);
               FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
               phys_y_qp(q) = phys(1);
               EvalDirichletFunc(attr, phys, time, u_D_val);
               for (int c = 0; c < dim; c++)
               {
                  u_D_3d(c * nq_dir + q) = u_D_val(c);
               }
            }

            Vector elvec_dir;
            dir_integ.AssembleBoundaryFaceRHS(*fe, *FTr, u_D_3d, elvec_dir);
            DebugDumpFaceData(-1, face_idx, "dirichlet_boundary", FTr,
                              &phys_y_qp, &u_D_3d, elvec_dir, nullptr);

            for (int j = 0; j < elvec_dir.Size(); j++)
               elvec(j) += elvec_dir(j);
         }
         else  // BR2
         {
            // BR2: compute u_D from face centroid
            const IntegrationPoint &ip_c = Geometries.GetCenter(FTr->GetGeometryType());
            FTr->Face->SetIntPoint(&ip_c);
            Vector fc_br2(dim);
            FTr->Face->Transform(ip_c, fc_br2);
            Vector u_D_br2;
            EvalDirichletFunc(attr, fc_br2, time, u_D_br2);
            real_t u_D[3] = {u_D_br2(0), u_D_br2(1), u_D_br2(2)};

            const DenseMatrix &Minv = elem_mass_inv_[FTr->Elem1No];
            int nqp = ir.GetNPoints();

            // Precompute shapes and normals
            DenseMatrix shapes(ndof, nqp);
            Vector nor_arr(dim * nqp), w_arr(nqp);

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip = FTr->GetElement1IntPoint();

               Vector sq(shapes.GetColumn(q), ndof);
               fe->CalcShape(eip, sq);

               Vector nq(&nor_arr[q * dim], dim);
               CalcOrtho(FTr->Jacobian(), nq);
               w_arr[q] = ip.weight;
            }

            // Compute lifted Dirichlet: similar to slip but for boundary
            DenseMatrix face_int(dim * dim, ndof);
            face_int = 0.0;
            for (int q = 0; q < nqp; q++)
            {
               real_t wq = w_arr[q];
               for (int u = 0; u < dim; u++)
               {
                  for (int s = 0; s < dim; s++)
                  {
                     real_t n_s = nor_arr[q * dim + s];
                     real_t factor = u_D[u] * n_s * wq;
                     for (int m = 0; m < ndof; m++)
                     {
                        face_int(u * dim + s, m) += shapes(m, q) * factor;
                     }
                  }
               }
            }

            DenseMatrix f_lifted(dim * dim, ndof);
            MultABt(face_int, Minv, f_lifted);
            // No 0.5 for boundary (full factor)

            // f_lifted_q[i,q]
            DenseMatrix fl_q(dim, nqp);
            fl_q = 0.0;
            for (int q = 0; q < nqp; q++)
            {
               Vector n_q(dim);
               for (int d = 0; d < dim; d++) { n_q(d) = nor_arr[q * dim + d]; }

               for (int i = 0; i < dim; i++)
               {
                  real_t sum = 0.0;
                  for (int u = 0; u < dim; u++)
                  {
                     for (int s = 0; s < dim; s++)
                     {
                        real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * n_q(i)
                           + mu_val_ * ((i == u ? 1.0 : 0.0) * n_q(s)
                                        + (i == s ? 1.0 : 0.0) * n_q(u));
                        real_t ev = 0.0;
                        for (int m = 0; m < ndof; m++)
                        {
                           ev += shapes(m, q) * f_lifted(u * dim + s, m);
                        }
                        sum += tn * ev;
                     }
                  }
                  fl_q(i, q) = sum;
               }
            }

            // BR2 penalty: num faces per element (geometry-dependent)
            Geometry::Type br2_geom = mesh_.GetElementGeometry(FTr->Elem1No);
            real_t br2_pen = (br2_geom == Geometry::TETRAHEDRON)
                                 ? real_t(dim + 1) : real_t(2 * dim);

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip = FTr->GetElement1IntPoint();

               real_t wq = w_arr[q];
               Vector nor_q(dim);
               for (int d = 0; d < dim; d++) { nor_q[d] = nor_arr[q * dim + d]; }

               DenseMatrix dshape_ref(ndof, dim);
               fe->CalcDShape(eip, dshape_ref);
               DenseMatrix adjJ(dim);
               CalcAdjugate(FTr->Elem1->Jacobian(), adjJ);
               DenseMatrix dshape_adj(ndof, dim);
               Mult(dshape_ref, adjJ, dshape_adj);
               real_t detJ = FTr->Elem1->Weight();

               for (int k = 0; k < ndof; k++)
               {
                  real_t grad_dot_n = 0.0;
                  for (int d = 0; d < dim; d++)
                  {
                     grad_dot_n += dshape_adj(k, d) * nor_q(d);
                  }

                  for (int i = 0; i < dim; i++)
                  {
                     real_t sym_val = 0.0;
                     for (int u = 0; u < dim; u++)
                     {
                        real_t trac = lambda_val_ * dshape_adj(k, i) * nor_q(u)
                           + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                        + dshape_adj(k, u) * nor_q(i));
                        sym_val += trac * u_D[u];
                     }

                     int idx = i * ndof + k;
                     elvec(idx) += epsilon_ * wq * sym_val / detJ;
                     elvec(idx) += br2_pen * wq * shapes(k, q) * fl_q(i, q);
                  }
               }
            }
         }

         for (int j = 0; j < vdofs.Size(); j++)
         {
            int gj = vdofs[j];
            if (gj >= 0) { rhs(gj) += elvec(j); }
            else { rhs(-1 - gj) -= elvec(j); }
         }
      }

      // ---------------------------------------------------------------
      // Interior faces with Dirichlet BC (skeleton pattern).
      // In Tandem's BP5 mesh, Physical Surface(5) includes Y=0 faces
      // outside the fault. These are interior faces in MFEM and need
      // the DG skeleton RHS contribution from both elements.
      // Tandem's bp5.lua boundary(): at Y=0, neither y>1 nor y<-1
      // triggers, so u_D = (Vp*t, 0, 0) — the full plate velocity.
      // This represents a locked fault (no relative slip).
      // ---------------------------------------------------------------
      for (int fi = 0; fi < dirichlet_interior_faces_.Size(); fi++)
      {
         int f = dirichlet_interior_faces_[fi];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(f);
         if (FTr == nullptr) { continue; }

         // Evaluate bp5.lua boundary(x,y,z,t) at each quad point:
         //   y > 1000 m (1 km):  u_D = (Vp*t/2, 0, 0)
         //   y < -1000 m:        u_D = (-Vp*t/2, 0, 0)
         //   |y| <= 1000 m:      u_D = (Vp*t, 0, 0)
         //
         // Orientation sign applied per QP via ComputeSkeletonDirichletSign,
         // matching Tandem DGCurvilinearCommon.h:92-99.

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
         int ndof1 = fe1->GetDof();
         int ndof2 = fe2->GetDof();

         int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
         int quad_order_dir = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, quad_order_dir);

         Vector elvec1(vdofs1.Size()), elvec2(vdofs2.Size());
         elvec1 = 0.0;
         elvec2 = 0.0;

         if (method_ == DGMethod::IP)
         {
            DGElasticityIPCombinedIntegrator dir_integ(
               lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

            int quad_order_dir = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
            const IntegrationRule &ir_dir = IntRules.Get(
               FTr->FaceGeom, quad_order_dir);
            int nq_dir = ir_dir.GetNPoints();

            // Evaluate boundary(x,y,z,t) at each QP with per-QP orientation
            // sign, matching Tandem DGCurvilinearCommon.h:92-99.
            Vector u_D_3d(dim * nq_dir);
            u_D_3d = 0.0;
            Vector phys_y_qp(nq_dir);
            Vector u_D_val(dim);
            for (int q = 0; q < nq_dir; q++)
            {
               const IntegrationPoint &ipq = ir_dir.IntPoint(q);
               FTr->SetAllIntPoints(&ipq);
               Vector phys(dim);
               FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
               phys_y_qp(q) = phys(1);
               EvalDirichletFunc(dirichlet_interior_attrs_[fi], phys, time, u_D_val);
               // Per-QP orientation sign (Tandem DGCurvilinearCommon.h:97-98)
               real_t dir_sign = ComputeSkeletonDirichletSign(FTr);
               for (int c = 0; c < dim; c++)
               {
                  u_D_3d(c * nq_dir + q) = dir_sign * u_D_val(c);
               }
            }

            Vector ev1, ev2;
            dir_integ.AssembleSlipFaceRHS(*fe1, *fe2, *FTr, u_D_3d, ev1, ev2);
            DebugDumpFaceData(-1, f, "dirichlet_interior", FTr,
                              &phys_y_qp, &u_D_3d, ev1, &ev2);

            for (int j = 0; j < ev1.Size(); j++) { elvec1(j) += ev1(j); }
            for (int j = 0; j < ev2.Size(); j++) { elvec2(j) += ev2(j); }
         }
         else  // BR2
         {
            // BR2 path: compute u_D_int from face centroid (legacy)
            const IntegrationPoint &ip_c = Geometries.GetCenter(FTr->GetGeometryType());
            FTr->Face->SetIntPoint(&ip_c);
            Vector fc_br2(dim);
            FTr->Face->Transform(ip_c, fc_br2);
            Vector u_D_br2;
            EvalDirichletFunc(dirichlet_interior_attrs_[fi], fc_br2, time, u_D_br2);
            // Orientation sign for BR2 (face-constant, affine faces)
            const IntegrationPoint &ip_s = Geometries.GetCenter(FTr->GetGeometryType());
            FTr->SetAllIntPoints(&ip_s);
            real_t dir_sign_br2 = ComputeSkeletonDirichletSign(FTr);
            real_t u_D_int[3] = {dir_sign_br2 * u_D_br2(0),
                                 dir_sign_br2 * u_D_br2(1),
                                 dir_sign_br2 * u_D_br2(2)};

            const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
            const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];
            int nqp = ir.GetNPoints();

            // Precompute shapes and normals
            DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
            Vector nor_arr(dim * nqp), w_arr(nqp);

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);

               Vector sq1(shapes1.GetColumn(q), ndof1);
               fe1->CalcShape(FTr->GetElement1IntPoint(), sq1);

               Vector sq2(shapes2.GetColumn(q), ndof2);
               fe2->CalcShape(FTr->GetElement2IntPoint(), sq2);

               Vector nq(&nor_arr[q * dim], dim);
               CalcOrtho(FTr->Jacobian(), nq);
               w_arr[q] = ip.weight;
            }

            // Compute lifted Dirichlet for skeleton:
            // face_int_k[u*dim+s, m] = sum_q shape_k_m(q) * u_D[u] * n[s] * w[q] * 0.5
            //
            // CRITICAL: Both face_int1 and face_int2 use the SAME sign (+=).
            // This matches the bilinear form (DGElasticityBR2Integrator) and the
            // slip assembly (AssembleSlipContributionBR2), which both use:
            //   IntFace_11 += shapes1 * factor
            //   IntFace_21 += shapes2 * factor  (same sign)
            //
            // The u_D_int is a prescribed JUMP (u1-u2), not a per-element value.
            // The BR2 lifting of the jump into both elements uses the same
            // reference normal. The -sign for elem2's penalty is handled later
            // via "elvec2 -= penalty * ..." (not here in face_int).
            //
            // BUG FIX (v37): Previously face_int2 used -= (opposite sign), which:
            //   - In v34 (no cross-element): pushed elem2 in WRONG direction
            //     (resisting correct loading, causing delayed recurrence)
            //   - In v36 (cross-element + 0.5): caused eval1+eval2 cancellation
            //     → zero penalty → fault lockup
            DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
            face_int1 = 0.0;
            face_int2 = 0.0;
            for (int q = 0; q < nqp; q++)
            {
               real_t wq = w_arr[q];
               for (int u = 0; u < dim; u++)
               {
                  for (int s = 0; s < dim; s++)
                  {
                     real_t n_s = nor_arr[q * dim + s];
                     real_t factor = u_D_int[u] * n_s * wq * 0.5;
                     for (int m = 0; m < ndof1; m++)
                     {
                        face_int1(u * dim + s, m) += shapes1(m, q) * factor;
                     }
                     for (int m = 0; m < ndof2; m++)
                     {
                        face_int2(u * dim + s, m) += shapes2(m, q) * factor;
                     }
                  }
               }
            }

            DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
            MultABt(face_int1, Minv1, f_lifted1);
            MultABt(face_int2, Minv2, f_lifted2);

            // BR2 penalty
            Geometry::Type br2_geom = mesh_.GetElementGeometry(FTr->Elem1No);
            real_t br2_pen = (br2_geom == Geometry::TETRAHEDRON)
                                 ? real_t(dim + 1) : real_t(2 * dim);

            // Combined BR2 lifted function evaluation at each quadrature point.
            // Uses cross-element lifting: f_lifted_q = C:n:(r_1 + r_2)
            // where r_k = Minv_k * face_int_k is the lift into element k.
            // The 0.5 factor matches the bilinear form's L_q = 0.5 * sum
            // structure (face_int already contains the {ψ} average factor 0.5).
            DenseMatrix f_lifted_q(dim, nqp);
            f_lifted_q = 0.0;
            for (int q = 0; q < nqp; q++)
            {
               Vector n_q(dim);
               for (int d = 0; d < dim; d++) { n_q(d) = nor_arr[q * dim + d]; }

               for (int i = 0; i < dim; i++)
               {
                  real_t sum = 0.0;
                  for (int u = 0; u < dim; u++)
                  {
                     for (int s = 0; s < dim; s++)
                     {
                        real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * n_q(i)
                           + mu_val_ * ((i == u ? 1.0 : 0.0) * n_q(s)
                                        + (i == s ? 1.0 : 0.0) * n_q(u));

                        // Cross-element evaluation: sum over BOTH elements
                        real_t eval1 = 0.0, eval2 = 0.0;
                        for (int m = 0; m < ndof1; m++)
                        {
                           eval1 += shapes1(m, q) * f_lifted1(u * dim + s, m);
                        }
                        for (int m = 0; m < ndof2; m++)
                        {
                           eval2 += shapes2(m, q) * f_lifted2(u * dim + s, m);
                        }
                        sum += tn * (eval1 + eval2);
                     }
                  }
                  f_lifted_q(i, q) = 0.5 * sum;
               }
            }

            // Elem1 contribution
            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();

               real_t wq = w_arr[q];
               Vector nor_q(dim);
               for (int d = 0; d < dim; d++) { nor_q[d] = nor_arr[q * dim + d]; }

               DenseMatrix dshape_ref(ndof1, dim);
               fe1->CalcDShape(eip1, dshape_ref);
               DenseMatrix adjJ(dim);
               CalcAdjugate(FTr->Elem1->Jacobian(), adjJ);
               DenseMatrix dshape_adj(ndof1, dim);
               Mult(dshape_ref, adjJ, dshape_adj);
               real_t detJ = FTr->Elem1->Weight();

               for (int k = 0; k < ndof1; k++)
               {
                  real_t grad_dot_n = 0.0;
                  for (int d = 0; d < dim; d++)
                  {
                     grad_dot_n += dshape_adj(k, d) * nor_q(d);
                  }

                  for (int i = 0; i < dim; i++)
                  {
                     real_t sym_val = 0.0;
                     for (int u = 0; u < dim; u++)
                     {
                        real_t trac = lambda_val_ * dshape_adj(k, i) * nor_q(u)
                           + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                        + dshape_adj(k, u) * nor_q(i));
                        sym_val += trac * u_D_int[u];
                     }

                     int idx = i * ndof1 + k;
                     elvec1(idx) += epsilon_ * wq * sym_val / (2.0 * detJ);
                     elvec1(idx) += br2_pen * wq * shapes1(k, q) * f_lifted_q(i, q);
                  }
               }
            }

            // Elem2 contribution (opposite penalty sign)
            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               FTr->SetAllIntPoints(&ip);
               const IntegrationPoint &eip2 = FTr->GetElement2IntPoint();

               real_t wq = w_arr[q];
               Vector nor_q(dim);
               for (int d = 0; d < dim; d++) { nor_q[d] = nor_arr[q * dim + d]; }

               DenseMatrix dshape_ref(ndof2, dim);
               fe2->CalcDShape(eip2, dshape_ref);
               DenseMatrix adjJ(dim);
               CalcAdjugate(FTr->Elem2->Jacobian(), adjJ);
               DenseMatrix dshape_adj(ndof2, dim);
               Mult(dshape_ref, adjJ, dshape_adj);
               real_t detJ = FTr->Elem2->Weight();

               for (int k = 0; k < ndof2; k++)
               {
                  real_t grad_dot_n = 0.0;
                  for (int d = 0; d < dim; d++)
                  {
                     grad_dot_n += dshape_adj(k, d) * nor_q(d);
                  }

                  for (int i = 0; i < dim; i++)
                  {
                     real_t sym_val = 0.0;
                     for (int u = 0; u < dim; u++)
                     {
                        real_t trac = lambda_val_ * dshape_adj(k, i) * nor_q(u)
                           + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                        + dshape_adj(k, u) * nor_q(i));
                        sym_val += trac * u_D_int[u];
                     }

                     int idx = i * ndof2 + k;
                     elvec2(idx) += epsilon_ * wq * sym_val / (2.0 * detJ);
                     elvec2(idx) -= br2_pen * wq * shapes2(k, q) * f_lifted_q(i, q);
                  }
               }
            }
         }

         // Scatter to global RHS
         for (int j = 0; j < vdofs1.Size(); j++)
         {
            int gj = vdofs1[j];
            if (gj >= 0) { rhs(gj) += elvec1(j); }
            else { rhs(-1 - gj) -= elvec1(j); }
         }
         for (int j = 0; j < vdofs2.Size(); j++)
         {
            int gj = vdofs2[j];
            if (gj >= 0) { rhs(gj) += elvec2(j); }
            else { rhs(-1 - gj) -= elvec2(j); }
         }
      }

      // ---------------------------------------------------------------
      // Shared faces with Dirichlet BC (parallel only).
      // Same skeleton pattern as interior Dirichlet faces, but only
      // elem1 is local — elem2 is on a neighboring rank.
      // Each rank contributes its local elem1 side; the neighboring rank
      // handles that same face's elem2 side as ITS elem1.
      // ---------------------------------------------------------------
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());

         for (int fi = 0; fi < dirichlet_shared_faces_.Size(); fi++)
         {
            int sf = dirichlet_shared_faces_[fi];
            FaceElementTransformations *FTr =
               mesh_.GetSharedFaceTransformations(sf);
            if (FTr == nullptr) { continue; }

            // v55: Per-quad-point boundary evaluation (done below in IP block)

            // Only elem1 is local
            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);

            const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
            int ndof1 = fe1->GetDof();

            // Elem2 (face-neighbor)
            int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
            int ndof2 = fe2->GetDof();

            int face_order = std::max(fe1->GetOrder(), fe2->GetOrder());
            int quad_order_dir_sh = match_quad_order_ ? (2 * face_order) : (2 * face_order + 1);
            const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, quad_order_dir_sh);

            Vector elvec1(vdofs1.Size());
            elvec1 = 0.0;

            if (method_ == DGMethod::IP)
            {
               DGElasticityIPCombinedIntegrator dir_integ(
                  lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

               int quad_order_dir = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
               const IntegrationRule &ir_dir = IntRules.Get(
                  FTr->FaceGeom, quad_order_dir);
               int nq_dir = ir_dir.GetNPoints();

               // Per-QP boundary evaluation with orientation sign
               Vector u_D_3d(dim * nq_dir);
               u_D_3d = 0.0;
               Vector phys_y_qp(nq_dir);
               Vector u_D_val(dim);
               for (int q = 0; q < nq_dir; q++)
               {
                  const IntegrationPoint &ipq = ir_dir.IntPoint(q);
                  FTr->SetAllIntPoints(&ipq);
                  Vector phys(dim);
                  FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
                  phys_y_qp(q) = phys(1);
                  EvalDirichletFunc(dirichlet_shared_attrs_[fi], phys, time, u_D_val);
                  real_t dir_sign = ComputeSkeletonDirichletSign(FTr);
                  for (int c = 0; c < dim; c++)
                  {
                     u_D_3d(c * nq_dir + q) = dir_sign * u_D_val(c);
                  }
               }

               Vector ev1, ev2;
               dir_integ.AssembleSlipFaceRHS(*fe1, *fe2, *FTr, u_D_3d, ev1, ev2);
               DebugDumpFaceData(-1, mesh_.GetSharedFace(sf), "dirichlet_shared",
                                 FTr, &phys_y_qp, &u_D_3d, ev1, &ev2);
               for (int j = 0; j < ev1.Size(); j++) { elvec1(j) += ev1(j); }
               // ev2 goes to neighbor rank, not used here
            }
            else  // BR2
            {
               if (!mass_inv_computed_) { PrecomputeMassInverse(); }

               // BR2: compute u_D_int from face centroid (legacy)
               const IntegrationPoint &ip_c = Geometries.GetCenter(FTr->GetGeometryType());
               FTr->Face->SetIntPoint(&ip_c);
               Vector fc_br2(dim);
               FTr->Face->Transform(ip_c, fc_br2);
               Vector u_D_br2;
               EvalDirichletFunc(dirichlet_shared_attrs_[fi], fc_br2, time, u_D_br2);
               // Orientation sign for BR2 (face-constant, affine faces)
               FTr->SetAllIntPoints(&ip_c);
               real_t dir_sign_br2 = ComputeSkeletonDirichletSign(FTr);
               real_t u_D_int[3] = {dir_sign_br2 * u_D_br2(0),
                                    dir_sign_br2 * u_D_br2(1),
                                    dir_sign_br2 * u_D_br2(2)};

               const DenseMatrix &Minv1 = elem_mass_inv_[FTr->Elem1No];
               const DenseMatrix &Minv2 = elem_mass_inv_[FTr->Elem2No];
               int nqp = ir.GetNPoints();

               // Precompute shapes and normals
               DenseMatrix shapes1(ndof1, nqp), shapes2(ndof2, nqp);
               Vector nor_arr(dim * nqp), w_arr(nqp);

               for (int q = 0; q < nqp; q++)
               {
                  const IntegrationPoint &ip = ir.IntPoint(q);
                  FTr->SetAllIntPoints(&ip);

                  Vector sq1(shapes1.GetColumn(q), ndof1);
                  fe1->CalcShape(FTr->GetElement1IntPoint(), sq1);

                  Vector sq2(shapes2.GetColumn(q), ndof2);
                  fe2->CalcShape(FTr->GetElement2IntPoint(), sq2);

                  Vector nq(&nor_arr[q * dim], dim);
                  CalcOrtho(FTr->Jacobian(), nq);
                  w_arr[q] = ip.weight;
               }

               // Compute lifted Dirichlet for skeleton (both sides)
               DenseMatrix face_int1(dim * dim, ndof1), face_int2(dim * dim, ndof2);
               face_int1 = 0.0;
               face_int2 = 0.0;
               for (int q = 0; q < nqp; q++)
               {
                  real_t wq = w_arr[q];
                  for (int u = 0; u < dim; u++)
                  {
                     for (int s = 0; s < dim; s++)
                     {
                        real_t n_s = nor_arr[q * dim + s];
                        real_t factor = u_D_int[u] * n_s * wq * 0.5;
                        for (int m = 0; m < ndof1; m++)
                        {
                           face_int1(u * dim + s, m) += shapes1(m, q) * factor;
                        }
                        for (int m = 0; m < ndof2; m++)
                        {
                           face_int2(u * dim + s, m) += shapes2(m, q) * factor;
                        }
                     }
                  }
               }

               DenseMatrix f_lifted1(dim * dim, ndof1), f_lifted2(dim * dim, ndof2);
               MultABt(face_int1, Minv1, f_lifted1);
               MultABt(face_int2, Minv2, f_lifted2);

               // BR2 penalty
               Geometry::Type br2_geom = mesh_.GetElementGeometry(FTr->Elem1No);
               real_t br2_pen = (br2_geom == Geometry::TETRAHEDRON)
                                    ? real_t(dim + 1) : real_t(2 * dim);

               // Cross-element lifted function evaluation
               DenseMatrix f_lifted_q(dim, nqp);
               f_lifted_q = 0.0;
               for (int q = 0; q < nqp; q++)
               {
                  Vector n_q(dim);
                  for (int d = 0; d < dim; d++) { n_q(d) = nor_arr[q * dim + d]; }

                  for (int i = 0; i < dim; i++)
                  {
                     real_t sum = 0.0;
                     for (int u = 0; u < dim; u++)
                     {
                        for (int s = 0; s < dim; s++)
                        {
                           real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * n_q(i)
                              + mu_val_ * ((i == u ? 1.0 : 0.0) * n_q(s)
                                           + (i == s ? 1.0 : 0.0) * n_q(u));

                           real_t eval1 = 0.0, eval2 = 0.0;
                           for (int m = 0; m < ndof1; m++)
                           {
                              eval1 += shapes1(m, q) * f_lifted1(u * dim + s, m);
                           }
                           for (int m = 0; m < ndof2; m++)
                           {
                              eval2 += shapes2(m, q) * f_lifted2(u * dim + s, m);
                           }
                           sum += tn * (eval1 + eval2);
                        }
                     }
                     f_lifted_q(i, q) = 0.5 * sum;
                  }
               }

               // Elem1 contribution only
               for (int q = 0; q < nqp; q++)
               {
                  const IntegrationPoint &ip = ir.IntPoint(q);
                  FTr->SetAllIntPoints(&ip);
                  const IntegrationPoint &eip1 = FTr->GetElement1IntPoint();

                  real_t wq = w_arr[q];
                  Vector nor_q(dim);
                  for (int d = 0; d < dim; d++) { nor_q[d] = nor_arr[q * dim + d]; }

                  DenseMatrix dshape_ref(ndof1, dim);
                  fe1->CalcDShape(eip1, dshape_ref);
                  DenseMatrix adjJ(dim);
                  CalcAdjugate(FTr->Elem1->Jacobian(), adjJ);
                  DenseMatrix dshape_adj(ndof1, dim);
                  Mult(dshape_ref, adjJ, dshape_adj);
                  real_t detJ = FTr->Elem1->Weight();

                  for (int k = 0; k < ndof1; k++)
                  {
                     real_t grad_dot_n = 0.0;
                     for (int d = 0; d < dim; d++)
                     {
                        grad_dot_n += dshape_adj(k, d) * nor_q(d);
                     }

                     for (int i = 0; i < dim; i++)
                     {
                        real_t sym_val = 0.0;
                        for (int u = 0; u < dim; u++)
                        {
                           real_t trac = lambda_val_ * dshape_adj(k, i) * nor_q(u)
                              + mu_val_ * ((i == u ? 1.0 : 0.0) * grad_dot_n
                                           + dshape_adj(k, u) * nor_q(i));
                           sym_val += trac * u_D_int[u];
                        }

                        int idx = i * ndof1 + k;
                        elvec1(idx) += epsilon_ * wq * sym_val / (2.0 * detJ);
                        elvec1(idx) += br2_pen * wq * shapes1(k, q) * f_lifted_q(i, q);
                     }
                  }
               }
            }

            // Scatter elem1 to global RHS
            for (int j = 0; j < vdofs1.Size(); j++)
            {
               int gj = vdofs1[j];
               if (gj >= 0) { rhs(gj) += elvec1(j); }
               else { rhs(-1 - gj) -= elvec1(j); }
            }
         }
#endif
      }
   }
