// elasticity_operator_verify.inl
// Original range: Lines 276-723 of elasticity_operator.hpp (pre-decomposition)
// This file is included from elasticity_operator.hpp — do not include directly.

   // ========================================================================
   // Production-mesh verification diagnostics (triggered by --verify)
   // ========================================================================

   /// Verify ghost DOF communication on the production mesh.
   /// Sets owned DOFs to f(x2,x3) = 7*sin(x2) + 3*cos(x3), expands to
   /// ghost DOFs via MPI, and verifies all local DOFs match expected values
   /// computed from independently-evaluated local coordinates.
   /// Returns global max error across all ranks.
   real_t VerifyGhostDOFCommunication() const
   {
      real_t max_err = 0.0;
      int n_mismatch = 0;

      // Ranks with fault DOFs run the full test; ranks without
      // participate in the collective MPI_Allreduce with zero values.
      if (num_fault_dofs_ > 0)
      {
         Vector local_x2, local_x3;
         GetFaultCoords2D(local_x2, local_x3);

         Vector owned_x2, owned_x3;
         RestrictToOwnedFault(local_x2, owned_x2);
         RestrictToOwnedFault(local_x3, owned_x3);

         Vector owned_data(2 * num_owned_fault_dofs_);
         for (int i = 0; i < num_owned_fault_dofs_; i++)
         {
            owned_data(2 * i)     = 7.0 * std::sin(owned_x2(i))
                                   + 3.0 * std::cos(owned_x3(i));
            owned_data(2 * i + 1) = 2.0 * owned_x2(i) - 5.0 * owned_x3(i);
         }

         Vector local_data;
         ExpandOwnedToLocalFault(owned_data, local_data, 2);

         Vector expected(2 * num_fault_dofs_);
         for (int i = 0; i < num_fault_dofs_; i++)
         {
            expected(2 * i)     = 7.0 * std::sin(local_x2(i))
                                 + 3.0 * std::cos(local_x3(i));
            expected(2 * i + 1) = 2.0 * local_x2(i) - 5.0 * local_x3(i);
         }

         for (int i = 0; i < num_fault_dofs_; i++)
         {
            real_t err = std::max(
               std::abs(local_data(2 * i)     - expected(2 * i)),
               std::abs(local_data(2 * i + 1) - expected(2 * i + 1)));
            max_err = std::max(max_err, err);
            if (err > 1e-10) { n_mismatch++; }
         }
      }

      // All ranks participate in the collective reduction
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         real_t global_max = 0.0;
         MPI_Allreduce(&max_err, &global_max, 1, MPI_DOUBLE, MPI_MAX,
                       mesh_.GetComm());
         int global_mismatch = 0;
         MPI_Allreduce(&n_mismatch, &global_mismatch, 1, MPI_INT, MPI_SUM,
                       mesh_.GetComm());
         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         if (rank == 0)
         {
            mfem::out << "  [VERIFY] Ghost DOF communication: max_err="
                      << std::scientific << std::setprecision(6) << global_max
                      << ", mismatches=" << global_mismatch
                      << " -> " << (global_max < 1e-10 ? "PASS" : "FAIL")
                      << "\n";
         }
         return global_max;
#endif
      }
      return max_err;
   }

   /// Diagnose the Dirichlet skip-set mechanism in AssembleDirichletLoading.
   /// Counts how many attr=5 boundary elements are: processed by boundary
   /// loop, skipped by dir_interior_set, skipped by dir_shared_set, or
   /// skipped because FTr==null. Mismatches between serial and parallel
   /// indicate the skip set misses some shared Dirichlet faces.
   void VerifyDirichletSkipSets() const
   {
      int n_attr5 = 0, n_skip_interior = 0, n_skip_shared = 0;
      int n_skip_null = 0, n_processed = 0;

      // Rebuild the same skip sets as AssembleDirichletLoading
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
            dir_shared_set.insert(
               mesh_.GetSharedFace(dirichlet_shared_faces_[fi]));
         }
#endif
      }

      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (dirichlet_bdr_marker_[attr - 1] != 1) { continue; }
         n_attr5++;

         int face_idx, face_info;
         mesh_.GetBdrElementFace(be, &face_idx, &face_info);

         if (dir_interior_set.count(face_idx) > 0)
         {
            n_skip_interior++;
            continue;
         }
         if (dir_shared_set.count(face_idx) > 0)
         {
            n_skip_shared++;
            continue;
         }

         FaceElementTransformations *FTr =
            mesh_.GetFaceElementTransformations(face_idx);
         if (FTr == nullptr) { n_skip_null++; continue; }

         n_processed++;
      }

      int n_dir_interior = dirichlet_interior_faces_.Size();
      int n_dir_shared = dirichlet_shared_faces_.Size();

      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         // Global sums
         int locals[6] = {n_attr5, n_skip_interior, n_skip_shared,
                          n_skip_null, n_processed, n_dir_shared};
         int globals[6] = {0};
         MPI_Allreduce(locals, globals, 6, MPI_INT, MPI_SUM,
                       mesh_.GetComm());
         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         if (rank == 0)
         {
            mfem::out << "  [VERIFY] Dirichlet skip-set audit:\n"
                      << "    attr=5 bdr elems:      " << globals[0] << "\n"
                      << "    skipped (interior):    " << globals[1] << "\n"
                      << "    skipped (shared set):  " << globals[2] << "\n"
                      << "    skipped (FTr null):    " << globals[3] << "\n"
                      << "    processed (boundary):  " << globals[4] << "\n"
                      << "    dirichlet_shared list: " << globals[5] << "\n"
                      << "    expected: processed + interior_skip + shared_skip"
                      << " + null_skip = attr5\n"
                      << "    actual:   " << globals[4] << " + " << globals[1]
                      << " + " << globals[2] << " + " << globals[3]
                      << " = " << (globals[4]+globals[1]+globals[2]+globals[3])
                      << (globals[4]+globals[1]+globals[2]+globals[3] == globals[0]
                          ? " -> OK" : " -> MISMATCH") << "\n";
            // Check if shared set size matches skip count
            if (globals[2] != globals[5])
            {
               mfem::out << "    WARNING: shared_skip ("
                         << globals[2] << ") != dirichlet_shared list ("
                         << globals[5] << ") -> some shared faces NOT skipped"
                         << " from boundary loop!\n";
            }
         }
#endif
      }
      else
      {
         mfem::out << "  [VERIFY] Dirichlet skip-set audit (serial):\n"
                   << "    attr=5 bdr elems: " << n_attr5
                   << ", skip_interior: " << n_skip_interior
                   << ", processed: " << n_processed << "\n";
      }
   }

   /// Per-face diagnostic for shared Dirichlet faces.
   /// For each shared Dirichlet face, dumps face normal, dir_sign, and
   /// ev1/ev2 norms from both ranks. Gathers via Allgather and compares
   /// rank A's ||ev2|| against rank B's ||ev1|| for the same physical face.
   void VerifySharedDirichletPerFace(real_t time) const
   {
      if constexpr (!IsParallelMesh<MeshType>::value) { return; }
#ifdef MFEM_USE_MPI
      if (dirichlet_shared_faces_.Size() == 0 &&
          mesh_.GetNRanks() > 1)
      {
         // Some ranks may have 0 shared Dirichlet faces — still participate
      }

      int dim = 3;
      int n_sh = dirichlet_shared_faces_.Size();

      // Get global vertex IDs for face keying
      Array<HYPRE_BigInt> gvert;
      mesh_.GetGlobalVertexIndices(gvert);

      // Pack: [gv0, gv1, gv2, nor0, nor1, nor2, dir_sign,
      //        ev1_norm, ev2_norm, rank] = 10 doubles per face
      const int ENTRY = 10;
      std::vector<double> local_pack(ENTRY * n_sh, 0.0);

      DGElasticityIPCombinedIntegrator dir_integ(
         lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
      auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());

      for (int fi = 0; fi < n_sh; fi++)
      {
         int sf = dirichlet_shared_faces_[fi];
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(sf);
         if (!FTr) { continue; }

         // Face vertex key
         int lf = mesh_.GetSharedFace(sf);
         Array<int> verts;
         mesh_.GetFaceVertices(lf, verts);
         HYPRE_BigInt gv[3] = {0, 0, 0};
         for (int j = 0; j < 3 && j < verts.Size(); j++)
            gv[j] = gvert[verts[j]];
         if (gv[0] > gv[1]) { std::swap(gv[0], gv[1]); }
         if (gv[1] > gv[2]) { std::swap(gv[1], gv[2]); }
         if (gv[0] > gv[1]) { std::swap(gv[0], gv[1]); }

         // Normal and dir_sign at face centroid
         const IntegrationPoint &ip0 =
            Geometries.GetCenter(FTr->GetGeometryType());
         FTr->SetAllIntPoints(&ip0);
         Vector nor(3);
         CalcOrtho(FTr->Jacobian(), nor);
         real_t dir_sign = ComputeSkeletonDirichletSign(FTr);

         // Compute ev1, ev2 with the same formula as the assembly
         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         int nbr_idx = FTr->Elem2No - mesh_.GetNE();
         const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);

         int qo = 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1;
         const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, qo);
         int nq = ir.GetNPoints();

         real_t Vh = Vp_ * time;
         // y=0 faces: no scaling (|y| < 1000)
         Vector u_D_3d(dim * nq);
         u_D_3d = 0.0;
         for (int q = 0; q < nq; q++)
         {
            FTr->SetAllIntPoints(&ir.IntPoint(q));
            real_t ds = ComputeSkeletonDirichletSign(FTr);
            u_D_3d(0 * nq + q) = ds * Vh;
         }

         Vector ev1, ev2;
         dir_integ.AssembleSlipFaceRHS(*fe1, *fe2, *FTr, u_D_3d, ev1, ev2);

         int b = ENTRY * fi;
         local_pack[b + 0] = static_cast<double>(gv[0]);
         local_pack[b + 1] = static_cast<double>(gv[1]);
         local_pack[b + 2] = static_cast<double>(gv[2]);
         local_pack[b + 3] = nor(0);
         local_pack[b + 4] = nor(1);
         local_pack[b + 5] = nor(2);
         local_pack[b + 6] = dir_sign;
         local_pack[b + 7] = ev1.Norml2();
         local_pack[b + 8] = ev2.Norml2();
         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         local_pack[b + 9] = static_cast<double>(rank);
      }

      // Allgather
      int lc = ENTRY * n_sh;
      int nranks = mesh_.GetNRanks();
      std::vector<int> rc(nranks), rd(nranks);
      MPI_Allgather(&lc, 1, MPI_INT, rc.data(), 1, MPI_INT,
                    mesh_.GetComm());
      int total = 0;
      for (int r = 0; r < nranks; r++)
      {
         rd[r] = total;
         total += rc[r];
      }
      std::vector<double> all(total);
      MPI_Allgatherv(local_pack.data(), lc, MPI_DOUBLE,
                     all.data(), rc.data(), rd.data(),
                     MPI_DOUBLE, mesh_.GetComm());

      // Group by face key, compare pairs
      struct VKey3 {
         int64_t v[3];
         bool operator<(const VKey3 &o) const {
            if (v[0]!=o.v[0]) return v[0]<o.v[0];
            if (v[1]!=o.v[1]) return v[1]<o.v[1];
            return v[2]<o.v[2];
         }
      };
      std::map<VKey3, std::vector<int>> key_entries;
      int n_entries = total / ENTRY;
      for (int e = 0; e < n_entries; e++)
      {
         VKey3 key;
         key.v[0] = static_cast<int64_t>(all[ENTRY*e + 0]);
         key.v[1] = static_cast<int64_t>(all[ENTRY*e + 1]);
         key.v[2] = static_cast<int64_t>(all[ENTRY*e + 2]);
         key_entries[key].push_back(e);
      }

      int rank = 0;
      MPI_Comm_rank(mesh_.GetComm(), &rank);
      if (rank == 0)
      {
         mfem::out << "  [VERIFY] Shared Dirichlet per-face diagnostic"
                   << " (t=" << time << "):\n";
         int n_pairs = 0, n_sign_same = 0, n_ev_mismatch = 0;
         for (auto &kv : key_entries)
         {
            auto &entries = kv.second;
            if (entries.size() != 2) { continue; }
            n_pairs++;
            int b0 = ENTRY * entries[0], b1 = ENTRY * entries[1];

            real_t ds0 = all[b0+6], ds1 = all[b1+6];
            real_t ev1_0 = all[b0+7], ev2_0 = all[b0+8];
            real_t ev1_1 = all[b1+7], ev2_1 = all[b1+8];
            int r0 = static_cast<int>(all[b0+9]);
            int r1 = static_cast<int>(all[b1+9]);

            bool sign_same = (ds0 * ds1 > 0);
            if (sign_same) { n_sign_same++; }

            // rank A's ev2 should match rank B's ev1
            real_t ev_scale = std::max({ev1_0, ev2_0, ev1_1, ev2_1, 1e-30});
            real_t cross_err = std::max(
               std::abs(ev2_0 - ev1_1) / ev_scale,
               std::abs(ev1_0 - ev2_1) / ev_scale);
            if (cross_err > 0.01) { n_ev_mismatch++; }

            auto &k = kv.first;
            mfem::out << "    face(" << k.v[0] << "," << k.v[1] << ","
                      << k.v[2] << "):"
                      << " r" << r0 << " ds=" << std::setw(2) << ds0
                      << " nor_y=" << std::setprecision(4) << all[b0+4]
                      << " |ev1|=" << std::setprecision(6) << ev1_0
                      << " |ev2|=" << ev2_0
                      << "  |  r" << r1 << " ds=" << std::setw(2) << ds1
                      << " nor_y=" << std::setprecision(4) << all[b1+4]
                      << " |ev1|=" << std::setprecision(6) << ev1_1
                      << " |ev2|=" << ev2_1
                      << (sign_same ? " SIGN_SAME!" : "")
                      << (cross_err > 0.01 ? " EV_MISMATCH!" : "")
                      << "\n";
         }
         mfem::out << "    summary: " << n_pairs << " pairs, "
                   << n_sign_same << " with SAME dir_sign, "
                   << n_ev_mismatch << " with ev cross-mismatch\n";
         if (n_sign_same > 0)
         {
            mfem::out << "    -> dir_sign does NOT flip on " << n_sign_same
                      << " faces! This is the bug.\n";
         }
      }
#endif
   }

   /// Verify serial-parallel Dirichlet loading consistency.
   /// At the given time, assembles the full RHS (slip + Dirichlet) and
   /// reports the global RHS norm. When run at 1 rank vs N ranks, the
   /// norms should match. Logs b_slip and b_dirichlet norms separately.
   void VerifyRHSNorms(real_t time, const Vector &slip_bc) const
   {
      if (!stiffness_assembled_) { return; }

      // Assemble slip RHS
      Vector rhs_slip(fes_->GetTrueVSize());
      rhs_slip = 0.0;
      if (method_ == DGMethod::IP)
      {
         AssembleSlipContributionIP(rhs_slip, slip_bc);
         AssembleSlipContributionIPShared(rhs_slip, slip_bc,
                                          fault_interior_faces_.Size());
      }
      else
      {
         AssembleSlipContributionBR2(rhs_slip, slip_bc);
         AssembleSlipContributionBR2Shared(rhs_slip, slip_bc,
                                           fault_interior_faces_.Size());
      }
      real_t slip_norm = rhs_slip.Norml2();

      // Assemble Dirichlet RHS
      Vector rhs_dir(fes_->GetTrueVSize());
      rhs_dir = 0.0;
      // Copy the slip RHS to get the total, then subtract to isolate Dirichlet
      Vector rhs_total = rhs_slip;
      AssembleDirichletLoading(rhs_total, time);
      for (int i = 0; i < rhs_dir.Size(); i++)
      {
         rhs_dir(i) = rhs_total(i) - rhs_slip(i);
      }
      real_t dir_norm = rhs_dir.Norml2();
      real_t total_norm = rhs_total.Norml2();

      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         // Use global L2 norms
         real_t local_s2 = slip_norm * slip_norm;
         real_t local_d2 = dir_norm * dir_norm;
         real_t local_t2 = total_norm * total_norm;
         real_t global_s2 = 0, global_d2 = 0, global_t2 = 0;
         MPI_Allreduce(&local_s2, &global_s2, 1, MPI_DOUBLE, MPI_SUM,
                       mesh_.GetComm());
         MPI_Allreduce(&local_d2, &global_d2, 1, MPI_DOUBLE, MPI_SUM,
                       mesh_.GetComm());
         MPI_Allreduce(&local_t2, &global_t2, 1, MPI_DOUBLE, MPI_SUM,
                       mesh_.GetComm());
         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         if (rank == 0)
         {
            mfem::out << "  [VERIFY] RHS norms at t=" << time << "s:"
                      << " ||b_slip||=" << std::scientific << std::setprecision(12)
                      << std::sqrt(global_s2)
                      << " ||b_dir||=" << std::sqrt(global_d2)
                      << " ||b_total||=" << std::sqrt(global_t2)
                      << "\n";
         }
#endif
      }
      else
      {
         mfem::out << "  [VERIFY] RHS norms at t=" << time << "s:"
                   << " ||b_slip||=" << std::scientific << std::setprecision(12)
                   << slip_norm
                   << " ||b_dir||=" << dir_norm
                   << " ||b_total||=" << total_norm
                   << "\n";
      }
   }
