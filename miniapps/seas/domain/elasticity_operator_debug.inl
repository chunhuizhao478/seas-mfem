// elasticity_operator_debug.inl
// Original range: Lines 758-1605 of elasticity_operator.hpp (pre-decomposition)
// This file is included from elasticity_operator.hpp — do not include directly.

   void ComputeTractionImpl(const GridFuncType &displacement,
                            const Vector &slip_bc,
                            Vector &traction,
                            Vector *normal_traction,
                            Vector *traction_stress_out,
                            Vector *traction_correction_out,
                            Vector *jump_residual_out,
                            Vector *normal_stress_out = nullptr,
                            Vector *normal_correction_out = nullptr);

   int DebugRank() const
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         return rank;
#endif
      }
      return 0;
   }

   bool DebugEnabledForTime(real_t time) const
   {
      if (!first_step_debug_.enabled || first_step_debug_done_) { return false; }
      // Target: dump at first accepted step (last RK45 stage).
      // With dt_init=0.02, the last stage is at t=0.02.
      // Threshold skips intermediate RK45 stages (t < 0.019).
      if (time < 0.019) { return false; }
      return (first_step_debug_.target_rank < 0 ||
              DebugRank() == first_step_debug_.target_rank);
   }

   const char *DebugPhaseName() const
   {
      switch (debug_phase_)
      {
         case DebugAssemblePhase::Slip: return "slip";
         case DebugAssemblePhase::Dirichlet: return "dirichlet";
         default: return "none";
      }
   }

   std::string DebugFilePath(const std::string &stem) const
   {
      std::ostringstream oss;
      oss << first_step_debug_.output_dir;
      if (!first_step_debug_.output_dir.empty() &&
          first_step_debug_.output_dir.back() != '/')
      {
         oss << "/";
      }
      oss << stem << "_r" << DebugRank() << ".csv";
      return oss.str();
   }

   const std::set<int> &DebugTargetLocalElements() const
   {
      if (debug_target_elems_cached_) { return debug_target_elems_; }
      debug_target_elems_.clear();
      const int num_local = mesh_.GetNE();
      const int interior_count = fault_interior_faces_.Size();
      for (int fault_idx : first_step_debug_.target_fault_faces)
      {
         if (fault_idx < 0) { continue; }
         FaceElementTransformations *FTr = nullptr;
         if (fault_idx < interior_count)
         {
            int face = fault_interior_faces_[fault_idx];
            FTr = mesh_.GetInteriorFaceTransformations(face);
         }
         else if constexpr (IsParallelMesh<MeshType>::value)
         {
            int sh_idx = fault_idx - interior_count;
            if (sh_idx >= 0 && sh_idx < fault_shared_faces_.Size())
            {
               FTr = mesh_.GetSharedFaceTransformations(fault_shared_faces_[sh_idx]);
            }
         }
         if (!FTr) { continue; }
         if (FTr->Elem1No >= 0 && FTr->Elem1No < num_local) { debug_target_elems_.insert(FTr->Elem1No); }
         if (FTr->Elem2No >= 0 && FTr->Elem2No < num_local) { debug_target_elems_.insert(FTr->Elem2No); }
      }
      debug_target_elems_cached_ = true;
      return debug_target_elems_;
   }

   bool DebugShouldDumpFace(int logical_fault_idx,
                            FaceElementTransformations *FTr) const
   {
      if (!DebugEnabledForTime(debug_time_) || FTr == nullptr) { return false; }
      if (logical_fault_idx >= 0 &&
          first_step_debug_.target_fault_faces.count(logical_fault_idx) > 0)
      {
         return true;
      }

      const auto &target_elems = DebugTargetLocalElements();
      return target_elems.count(FTr->Elem1No) > 0 ||
             target_elems.count(FTr->Elem2No) > 0;
   }

   void DebugDumpFaceData(int logical_fault_idx,
                          int mesh_face_idx,
                          const char *face_kind,
                          FaceElementTransformations *FTr,
                          const Vector *phys_y_qp,
                          const Vector *input_qp,
                          const Vector &elvec1,
                          const Vector *elvec2) const
   {
      if (!DebugShouldDumpFace(logical_fault_idx, FTr)) { return; }

      // Truncate on first write to discard stale data from previous runs
      std::ofstream out(DebugFilePath("first_step_face_rhs"),
                        debug_face_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_face_header_written_)
      {
         out << "time,phase,face_kind,fault_idx,mesh_face,elem1,elem2,record,index0,index1,value\n";
         debug_face_header_written_ = true;
      }

      const int elem1 = FTr ? FTr->Elem1No : -1;
      const int elem2 = FTr ? FTr->Elem2No : -1;
      auto write_row = [&](const char *record, int i0, int i1, real_t value)
      {
         out << std::setprecision(17) << debug_time_ << ","
             << DebugPhaseName() << ","
             << face_kind << ","
             << logical_fault_idx << ","
             << mesh_face_idx << ","
             << elem1 << ","
             << elem2 << ","
             << record << ","
             << i0 << ","
             << i1 << ","
             << value << "\n";
      };

      if (phys_y_qp)
      {
         for (int q = 0; q < phys_y_qp->Size(); q++)
         {
            write_row("phys_y", q, -1, (*phys_y_qp)(q));
         }
      }
      if (input_qp)
      {
         const int nq = input_qp->Size() / 3;
         for (int q = 0; q < nq; q++)
         {
            for (int c = 0; c < 3; c++)
            {
               write_row("input_qp", q, c, (*input_qp)(c * nq + q));
            }
         }
      }
      for (int j = 0; j < elvec1.Size(); j++)
      {
         write_row("elvec1", j, -1, elvec1(j));
      }
      if (elvec2)
      {
         for (int j = 0; j < elvec2->Size(); j++)
         {
            write_row("elvec2", j, -1, (*elvec2)(j));
         }
      }
   }

   void DebugDumpElementVector(const char *quantity, const Vector &vec) const
   {
      if (!DebugEnabledForTime(debug_time_)) { return; }

      std::ofstream out(DebugFilePath("first_step_elem_data"),
                        debug_elem_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_elem_header_written_)
      {
         out << "time,quantity,elem,component,local_dof,vdof,value\n";
         debug_elem_header_written_ = true;
      }

      const auto &target_elems = DebugTargetLocalElements();
      for (int elem : target_elems)
      {
         Array<int> vdofs;
         fes_->GetElementVDofs(elem, vdofs);
         const int ndof = scalar_fes_->GetFE(elem)->GetDof();
         for (int j = 0; j < vdofs.Size(); j++)
         {
            const int vdof = vdofs[j];
            const int lid = (vdof >= 0) ? vdof : (-1 - vdof);
            const real_t value = (vdof >= 0) ? vec(lid) : -vec(lid);
            out << std::setprecision(17) << debug_time_ << ","
                << quantity << ","
                << elem << ","
                << (j / ndof) << ","
                << (j % ndof) << ","
                << vdof << ","
                << value << "\n";
         }
      }
   }

   /// Dump per-element and per-face K matrix contributions for target elements.
   /// Each contribution is dumped separately so we can compare against Tandem's
   /// assemble_volume / assemble_skeleton / assemble_boundary piece by piece.
   /// MPI-safe: only local elements and local interior faces.
   void DebugDumpKContributions() const
   {
      if (!DebugEnabledForTime(debug_time_) || method_ != DGMethod::IP) { return; }

      const int dim = 3;
      const auto &target_elems = DebugTargetLocalElements();
      if (target_elems.empty()) { return; }

      auto write_mat = [](std::ofstream &os, const char *source, int id,
                          int e1, int e2, const char *block,
                          const DenseMatrix &M)
      {
         for (int i = 0; i < M.Height(); i++)
         {
            for (int j = 0; j < M.Width(); j++)
            {
               real_t v = M(i, j);
               if (std::abs(v) > 1e-30)
               {
                  os << source << "," << id << "," << e1 << "," << e2
                     << "," << block << "," << i << "," << j << ","
                     << std::setprecision(17) << v << "\n";
               }
            }
         }
      };

      // 1. Volume contributions: ElasticityIntegrator on target elements
      {
         std::ofstream out(DebugFilePath("first_step_K_volume"), std::ios::trunc);
         if (out)
         {
            out << "source,elem,elem1,elem2,block,row,col,value\n";
            ElasticityIntegrator vol_integ(lambda_coeff_, mu_coeff_);
            for (int elem : target_elems)
            {
               const FiniteElement *fe = scalar_fes_->GetFE(elem);
               ElementTransformation *eltrans =
                  mesh_.GetElementTransformation(elem);
               DenseMatrix K_vol;
               vol_integ.AssembleElementMatrix(*fe, *eltrans, K_vol);
               write_mat(out, "volume", elem, elem, -1, "A00", K_vol);
            }
         }
      }

      // 2. Interior face contributions: DGElasticityIPCombinedIntegrator
      {
         std::ofstream out(DebugFilePath("first_step_K_skeleton"),
                           std::ios::trunc);
         if (out)
         {
            out << "source,face,elem1,elem2,block,row,col,value\n";
            DGElasticityIPCombinedIntegrator face_integ(
               lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);

            int nfaces = mesh_.GetNumFaces();
            for (int f = 0; f < nfaces; f++)
            {
               FaceElementTransformations *FTr =
                  mesh_.GetInteriorFaceTransformations(f);
               if (!FTr) { continue; }

               bool e1_target = target_elems.count(FTr->Elem1No) > 0;
               bool e2_target = (FTr->Elem2No >= 0 &&
                                 target_elems.count(FTr->Elem2No) > 0);
               if (!e1_target && !e2_target) { continue; }

               const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
               const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
               DenseMatrix K_face;
               face_integ.AssembleFaceMatrix(*fe1, *fe2, *FTr, K_face);

               int n1 = fe1->GetDof() * dim;
               int n2 = fe2->GetDof() * dim;
               DenseMatrix A00(n1, n1), A01(n1, n2), A10(n2, n1), A11(n2, n2);
               K_face.GetSubMatrix(0, n1, 0, n1, A00);
               K_face.GetSubMatrix(0, n1, n1, n1 + n2, A01);
               K_face.GetSubMatrix(n1, n1 + n2, 0, n1, A10);
               K_face.GetSubMatrix(n1, n1 + n2, n1, n1 + n2, A11);

               write_mat(out, "interior_face", f, FTr->Elem1No, FTr->Elem2No,
                         "A00", A00);
               write_mat(out, "interior_face", f, FTr->Elem1No, FTr->Elem2No,
                         "A01", A01);
               write_mat(out, "interior_face", f, FTr->Elem1No, FTr->Elem2No,
                         "A10", A10);
               write_mat(out, "interior_face", f, FTr->Elem1No, FTr->Elem2No,
                         "A11", A11);
            }

            // Boundary face contributions (if any target element touches one)
            for (int be = 0; be < mesh_.GetNBE(); be++)
            {
               int attr = mesh_.GetBdrAttribute(be);
               if (dirichlet_bdr_marker_.Size() > 0 &&
                   dirichlet_bdr_marker_[attr - 1] != 1) { continue; }

               int face_idx, face_info_val;
               mesh_.GetBdrElementFace(be, &face_idx, &face_info_val);
               FaceElementTransformations *FTr =
                  mesh_.GetFaceElementTransformations(face_idx);
               if (!FTr || target_elems.count(FTr->Elem1No) == 0)
               { continue; }

               const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
               DenseMatrix K_bdr;
               face_integ.AssembleFaceMatrix(*fe1, *fe1, *FTr, K_bdr);
               write_mat(out, "boundary_face", face_idx, FTr->Elem1No, -1,
                         "A00", K_bdr);
            }
         }
      }
   }

   /// Local-only variant: dumps interior fault faces only (no MPI exchange).
   /// Safe to call on a single rank without deadlocking.
   void DebugDumpFaultJumpsLocal(const GridFuncType &displacement,
                                 const Vector &slip_bc) const
   {
      if (!DebugEnabledForTime(debug_time_) || method_ != DGMethod::IP) { return; }

      std::ofstream out(DebugFilePath("first_step_face_jump"),
                        debug_jump_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_jump_header_written_)
      {
         out << "time,face_kind,fault_idx,mesh_face,elem1,elem2,record,q,component,value\n";
         debug_jump_header_written_ = true;
      }

      auto write_row = [&](const char *face_kind, int logical_fault_idx,
                           int mesh_face_idx, int elem1, int elem2,
                           const char *record, int q, int c, real_t value)
      {
         out << std::setprecision(17) << debug_time_ << ","
             << face_kind << "," << logical_fault_idx << ","
             << mesh_face_idx << "," << elem1 << "," << elem2 << ","
             << record << "," << q << "," << c << "," << value << "\n";
      };

      const int dim = 3;

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(face);
         if (!DebugShouldDumpFace(fi, FTr)) { continue; }

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);
         const int ndof1 = fe1->GetDof();
         const int ndof2 = fe2->GetDof();
         const int nq = face_quad_->NumQuadPoints();

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);
         Vector u1_all(vdofs1.Size()), u2_all(vdofs2.Size());
         displacement.GetSubVector(vdofs1, u1_all);
         displacement.GetSubVector(vdofs2, u2_all);

         Vector delta_u_quad;
         BuildSlipAtQuadPoints(fi, slip_bc, delta_u_quad);

         const IntegrationRule &ir = face_quad_->GetQuadRule();
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &fip = ir.IntPoint(q);
            FTr->SetAllIntPoints(&fip);

            Vector phys(dim);
            FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_x", q, -1, phys(0));
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_y", q, -1, phys(1));
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_z", q, -1, phys(2));

            Vector s1(ndof1), s2(ndof2);
            fe1->CalcShape(FTr->GetElement1IntPoint(), s1);
            fe2->CalcShape(FTr->GetElement2IntPoint(), s2);
            for (int c = 0; c < dim; c++)
            {
               real_t u1q = 0.0, u2q = 0.0;
               for (int k = 0; k < ndof1; k++) { u1q += s1(k) * u1_all(c * ndof1 + k); }
               for (int k = 0; k < ndof2; k++) { u2q += s2(k) * u2_all(c * ndof2 + k); }
               const real_t slipq = delta_u_quad(c * nq + q);
               write_row("fault_interior", fi, face,
                         FTr->Elem1No, FTr->Elem2No, "u1_q", q, c, u1q);
               write_row("fault_interior", fi, face,
                         FTr->Elem1No, FTr->Elem2No, "u2_q", q, c, u2q);
               write_row("fault_interior", fi, face,
                         FTr->Elem1No, FTr->Elem2No, "slip_q", q, c, slipq);
               write_row("fault_interior", fi, face,
                         FTr->Elem1No, FTr->Elem2No, "jump_minus_slip", q, c,
                         (u1q - u2q) - slipq);
            }
         }
      }
   }

   /// Local-only variant: dumps interior fault face traction (no MPI exchange).
   /// Local-only variant: computes traction per-face inline for interior
   /// faces only. No ComputeTraction call (which uses MPI for shared faces).
   void DebugDumpFaultTractionLocal(const GridFuncType &displacement,
                                     const Vector &slip_bc)
   {
      if (!DebugEnabledForTime(debug_time_) || method_ != DGMethod::IP) { return; }

      std::ofstream out(DebugFilePath("first_step_face_trac"),
                        debug_trac_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_trac_header_written_)
      {
         out << "time,face_kind,fault_idx,mesh_face,elem1,elem2,record,index0,index1,value\n";
         debug_trac_header_written_ = true;
      }

      auto write_row = [&](const char *face_kind, int logical_fault_idx,
                           int mesh_face_idx, int elem1, int elem2,
                           const char *record, int i0, int i1, real_t value)
      {
         out << std::setprecision(17) << debug_time_ << ","
             << face_kind << "," << logical_fault_idx << ","
             << mesh_face_idx << "," << elem1 << "," << elem2 << ","
             << record << "," << i0 << "," << i1 << "," << value << "\n";
      };

      const int dim = 3;
      const int nbf = nbf_per_face_;

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(face);
         if (!DebugShouldDumpFace(fi, FTr)) { continue; }

         const FiniteElement *fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         const FiniteElement *fe2 = scalar_fes_->GetFE(FTr->Elem2No);

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);
         Vector u1_all(vdofs1.Size()), u2_all(vdofs2.Size());
         displacement.GetSubVector(vdofs1, u1_all);
         displacement.GetSubVector(vdofs2, u2_all);

         Vector delta_u_quad;
         BuildSlipAtQuadPoints(fi, slip_bc, delta_u_quad);

         DGElasticityIPCombinedIntegrator trac_integ(
            lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
         Vector T_quad, nl_q;
         trac_integ.ComputeTractionAtQuadPoints(
            *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad, T_quad, nullptr, &nl_q);

         const int nq = T_quad.Size() / dim;
         const IntegrationRule &ir = IntRules.Get(
            FTr->GetGeometryType(), 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1);
         for (int q = 0; q < nq; q++)
         {
            FTr->SetAllIntPoints(&ir.IntPoint(q));
            Vector phys(dim);
            FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_x", q, -1, phys(0));
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_y", q, -1, phys(1));
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "phys_z", q, -1, phys(2));
            write_row("fault_interior", fi, face,
                      FTr->Elem1No, FTr->Elem2No, "nl_q", q, -1, nl_q(q));
            for (int c = 0; c < dim; c++)
            {
               write_row("fault_interior", fi, face,
                         FTr->Elem1No, FTr->Elem2No, "traction_q", q, c,
                         T_quad(c * nq + q));
            }
         }

         // Note: projected DOF-level traction (traction_dip/strike/normal)
         // is omitted here because it requires ComputeTraction which uses MPI.
         // The per-QP traction_q values above are sufficient for cross-code
         // comparison — project offline if needed.
      }
   }

   void DebugDumpFaultJumps(const GridFuncType &displacement,
                            const Vector &slip_bc) const
   {
      if (!DebugEnabledForTime(debug_time_) || method_ != DGMethod::IP) { return; }

      std::ofstream out(DebugFilePath("first_step_face_jump"),
                        debug_jump_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_jump_header_written_)
      {
         out << "time,face_kind,fault_idx,mesh_face,elem1,elem2,record,q,component,value\n";
         debug_jump_header_written_ = true;
      }

      auto write_row = [&](const char *face_kind, int logical_fault_idx,
                           int mesh_face_idx, int elem1, int elem2,
                           const char *record, int q, int c, real_t value)
      {
         out << std::setprecision(17) << debug_time_ << ","
             << face_kind << ","
             << logical_fault_idx << ","
             << mesh_face_idx << ","
             << elem1 << ","
             << elem2 << ","
             << record << ","
             << q << ","
             << c << ","
             << value << "\n";
      };

      auto dump_face = [&](int logical_fault_idx, int mesh_face_idx,
                           const char *face_kind, FaceElementTransformations *FTr,
                           const FiniteElement *fe1, const Vector &u1_all,
                           const FiniteElement *fe2, const Vector &u2_all)
      {
         if (!DebugShouldDumpFace(logical_fault_idx, FTr)) { return; }

         const int dim = 3;
         const int ndof1 = fe1->GetDof();
         const int ndof2 = fe2->GetDof();
         const int nq = face_quad_->NumQuadPoints();

         Vector delta_u_quad;
         BuildSlipAtQuadPoints(logical_fault_idx, slip_bc, delta_u_quad);

         const IntegrationRule &ir = face_quad_->GetQuadRule();
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &fip = ir.IntPoint(q);
            FTr->SetAllIntPoints(&fip);

            // Physical coordinates at this QP (for cross-code matching)
            Vector phys(dim);
            FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_x", q, -1, phys(0));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_y", q, -1, phys(1));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_z", q, -1, phys(2));

            Vector s1(ndof1), s2(ndof2);
            fe1->CalcShape(FTr->GetElement1IntPoint(), s1);
            fe2->CalcShape(FTr->GetElement2IntPoint(), s2);
            for (int c = 0; c < dim; c++)
            {
               real_t u1q = 0.0, u2q = 0.0;
               for (int k = 0; k < ndof1; k++) { u1q += s1(k) * u1_all(c * ndof1 + k); }
               for (int k = 0; k < ndof2; k++) { u2q += s2(k) * u2_all(c * ndof2 + k); }
               const real_t slipq = delta_u_quad(c * nq + q);
               const real_t jump_minus_slip = (u1q - u2q) - slipq;
               write_row(face_kind, logical_fault_idx, mesh_face_idx,
                         FTr->Elem1No, FTr->Elem2No, "u1_q", q, c, u1q);
               write_row(face_kind, logical_fault_idx, mesh_face_idx,
                         FTr->Elem1No, FTr->Elem2No, "u2_q", q, c, u2q);
               write_row(face_kind, logical_fault_idx, mesh_face_idx,
                         FTr->Elem1No, FTr->Elem2No, "slip_q", q, c, slipq);
               write_row(face_kind, logical_fault_idx, mesh_face_idx,
                         FTr->Elem1No, FTr->Elem2No,
                         "jump_minus_slip", q, c, jump_minus_slip);
            }
         }
      };

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(face);
         if (!DebugShouldDumpFace(fi, FTr)) { continue; }

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);
         Vector u1_all(vdofs1.Size()), u2_all(vdofs2.Size());
         displacement.GetSubVector(vdofs1, u1_all);
         displacement.GetSubVector(vdofs2, u2_all);
         dump_face(fi, face, "fault_interior", FTr,
                   scalar_fes_->GetFE(FTr->Elem1No), u1_all,
                   scalar_fes_->GetFE(FTr->Elem2No), u2_all);
      }

      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         if (!pfes) { return; }

         ParGridFunction par_u(pfes);
         par_u = displacement;
         pfes->ExchangeFaceNbrData();
         par_u.ExchangeFaceNbrData();

         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            const int logical_fault_idx = fault_interior_faces_.Size() + i;
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr = mesh_.GetSharedFaceTransformations(sf);
            if (!DebugShouldDumpFace(logical_fault_idx, FTr)) { continue; }

            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
            Vector u1_all(vdofs1.Size());
            par_u.GetSubVector(vdofs1, u1_all);

            const int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
            Array<int> vdofs2;
            pfes->GetFaceNbrElementVDofs(nbr_idx, vdofs2);
            const Vector &nbr_data = par_u.FaceNbrData();
            Vector u2_all(vdofs2.Size());
            for (int j = 0; j < vdofs2.Size(); j++)
            {
               u2_all(j) = nbr_data(vdofs2[j]);
            }

            dump_face(logical_fault_idx, mesh_.GetSharedFace(sf), "fault_shared", FTr,
                      scalar_fes_->GetFE(FTr->Elem1No), u1_all, fe2, u2_all);
         }
#endif
      }
   }

   void DebugDumpFaultTraction(const GridFuncType &displacement,
                               const Vector &slip_bc)
   {
      if (!DebugEnabledForTime(debug_time_) || method_ != DGMethod::IP) { return; }

      Vector traction, normal_traction;
      ComputeTraction(displacement, slip_bc, traction, &normal_traction);

      std::ofstream out(DebugFilePath("first_step_face_trac"),
                        debug_trac_header_written_ ? std::ios::app : std::ios::trunc);
      if (!out) { return; }
      if (!debug_trac_header_written_)
      {
         out << "time,face_kind,fault_idx,mesh_face,elem1,elem2,record,index0,index1,value\n";
         debug_trac_header_written_ = true;
      }

      auto write_row = [&](const char *face_kind, int logical_fault_idx,
                           int mesh_face_idx, int elem1, int elem2,
                           const char *record, int i0, int i1, real_t value)
      {
         out << std::setprecision(17) << debug_time_ << ","
             << face_kind << ","
             << logical_fault_idx << ","
             << mesh_face_idx << ","
             << elem1 << ","
             << elem2 << ","
             << record << ","
             << i0 << ","
             << i1 << ","
             << value << "\n";
      };

      auto dump_face = [&](int logical_fault_idx, int mesh_face_idx,
                           const char *face_kind, FaceElementTransformations *FTr,
                           const FiniteElement *fe1, const Vector &u1_all,
                           const FiniteElement *fe2, const Vector &u2_all)
      {
         if (!DebugShouldDumpFace(logical_fault_idx, FTr)) { return; }

         const int dim = 3;
         const int nbf = nbf_per_face_;
         const int nq = face_quad_->NumQuadPoints();

         Vector delta_u_quad;
         BuildSlipAtQuadPoints(logical_fault_idx, slip_bc, delta_u_quad);

         DGElasticityIPCombinedIntegrator trac_integ(
            lambda_coeff_, mu_coeff_, dim, epsilon_, penalty_factor_);
         Vector T_quad, nl_q;
         trac_integ.ComputeTractionAtQuadPoints(
            *fe1, *fe2, *FTr, u1_all, u2_all, delta_u_quad, T_quad, nullptr, &nl_q);

         // Compute physical coordinates at each QP for cross-code matching
         const IntegrationRule &ir_trac = IntRules.Get(
            FTr->GetGeometryType(), 2 * std::max(fe1->GetOrder(), fe2->GetOrder()) + 1);
         for (int q = 0; q < nq; q++)
         {
            FTr->SetAllIntPoints(&ir_trac.IntPoint(q));
            Vector phys(dim);
            FTr->Elem1->Transform(FTr->GetElement1IntPoint(), phys);
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_x", q, -1, phys(0));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_y", q, -1, phys(1));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "phys_z", q, -1, phys(2));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "nl_q", q, -1, nl_q(q));
            for (int c = 0; c < dim; c++)
            {
               write_row(face_kind, logical_fault_idx, mesh_face_idx,
                         FTr->Elem1No, FTr->Elem2No, "traction_q", q, c,
                         T_quad(c * nq + q));
            }
         }

         for (int kk = 0; kk < nbf; kk++)
         {
            const int dof_idx = logical_fault_idx * nbf + kk;
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "traction_dip", kk, -1,
                      traction(2 * dof_idx));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "traction_strike", kk, -1,
                      traction(2 * dof_idx + 1));
            write_row(face_kind, logical_fault_idx, mesh_face_idx,
                      FTr->Elem1No, FTr->Elem2No, "normal_traction", kk, -1,
                      normal_traction(dof_idx));
         }
      };

      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
      {
         int face = fault_interior_faces_[fi];
         FaceElementTransformations *FTr = mesh_.GetInteriorFaceTransformations(face);
         if (!DebugShouldDumpFace(fi, FTr)) { continue; }

         Array<int> vdofs1, vdofs2;
         fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
         fes_->GetElementVDofs(FTr->Elem2No, vdofs2);
         Vector u1_all(vdofs1.Size()), u2_all(vdofs2.Size());
         displacement.GetSubVector(vdofs1, u1_all);
         displacement.GetSubVector(vdofs2, u2_all);
         dump_face(fi, face, "fault_interior", FTr,
                   scalar_fes_->GetFE(FTr->Elem1No), u1_all,
                   scalar_fes_->GetFE(FTr->Elem2No), u2_all);
      }

      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         if (!pfes) { return; }

         ParGridFunction par_u(pfes);
         par_u = displacement;
         pfes->ExchangeFaceNbrData();
         par_u.ExchangeFaceNbrData();

         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            const int logical_fault_idx = fault_interior_faces_.Size() + i;
            int sf = fault_shared_faces_[i];
            FaceElementTransformations *FTr = mesh_.GetSharedFaceTransformations(sf);
            if (!DebugShouldDumpFace(logical_fault_idx, FTr)) { continue; }

            Array<int> vdofs1;
            fes_->GetElementVDofs(FTr->Elem1No, vdofs1);
            Vector u1_all(vdofs1.Size());
            par_u.GetSubVector(vdofs1, u1_all);

            const int nbr_idx = FTr->Elem2No - mesh_.GetNE();
            const FiniteElement *fe2 = pfes->GetFaceNbrFE(nbr_idx);
            Array<int> vdofs2;
            pfes->GetFaceNbrElementVDofs(nbr_idx, vdofs2);
            const Vector &nbr_data = par_u.FaceNbrData();
            Vector u2_all(vdofs2.Size());
            for (int j = 0; j < vdofs2.Size(); j++)
            {
               u2_all(j) = nbr_data(vdofs2[j]);
            }

            dump_face(logical_fault_idx, mesh_.GetSharedFace(sf), "fault_shared", FTr,
                      scalar_fes_->GetFE(FTr->Elem1No), u1_all, fe2, u2_all);
         }
#endif
      }
   }

   /// Build 3D slip at quadrature points for a single fault face.
   ///
   /// Shared by production assembly (AssembleSlipContributionIP) and debug
   /// diagnostics (DebugDumpFaultJumps, DebugDumpFaultTraction).  Single
   /// source of truth for the interpolation + embedding + sign chain.
   ///
   /// @param logical_fault_idx  Face index in the combined interior+shared list
   /// @param slip_bc            Full local slip vector [2 * num_fault_dofs]
   /// @param[out] delta_u_quad  3D slip at QPs [dim * nq], sign-corrected
   void BuildSlipAtQuadPoints(int logical_fault_idx, const Vector &slip_bc,
                              Vector &delta_u_quad) const
   {
      const int dim = 3;
      const int nbf = nbf_per_face_;
      const auto &basis = fault_basis_.GetBasis(logical_fault_idx);
      // Tandem convention: sign is baked into the basis vectors.
      // No separate sign factor needed.

      if (!basis.qp_data.empty())
      {
         Vector slip_tang(2 * nbf);
         for (int kk = 0; kk < nbf; kk++)
         {
            int dof_idx = logical_fault_idx * nbf + kk;
            slip_tang(0 * nbf + kk) = slip_bc(2 * dof_idx);
            slip_tang(1 * nbf + kk) = slip_bc(2 * dof_idx + 1);
         }
         Vector slip_tang_q;
         face_quad_->InterpolateToQuadPoints(2, slip_tang, slip_tang_q);
         int nqp = slip_tang_q.Size() / 2;
         delta_u_quad.SetSize(dim * nqp);
         for (int q = 0; q < nqp; q++)
         {
            real_t sl_q[2] = {slip_tang_q(q), slip_tang_q(nqp + q)};
            real_t du[3];
            fault_basis_.EmbedSlipQP(logical_fault_idx, q, sl_q, du);
            for (int c = 0; c < dim; c++)
               delta_u_quad(c * nqp + q) = du[c];
         }
      }
      else
      {
         Vector delta_u_nodal(dim * nbf);
         for (int kk = 0; kk < nbf; kk++)
         {
            int dof_idx = logical_fault_idx * nbf + kk;
            real_t slip_local[2] = {slip_bc(2 * dof_idx),
                                    slip_bc(2 * dof_idx + 1)};
            real_t du[3];
            fault_basis_.EmbedSlip(logical_fault_idx, slip_local, du);
            for (int c = 0; c < dim; c++)
               delta_u_nodal(c * nbf + kk) = du[c];
         }
         face_quad_->InterpolateToQuadPoints(dim, delta_u_nodal, delta_u_quad);
      }
   }
