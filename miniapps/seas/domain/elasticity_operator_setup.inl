// elasticity_operator_setup.inl
// Original range: Lines 1757-3052 of elasticity_operator.hpp (pre-decomposition)
// This file is included from elasticity_operator.hpp — do not include directly.

   // ========================================================================
   // Setup methods
   // ========================================================================

   void SetupFESpace()
   {
      fec_ = std::make_unique<DG_FECollection>(order_, 3, BasisType::GaussLobatto);
      // Vector DG space with 3 components, byNODES ordering
      fes_ = std::make_unique<FESpaceType>(&mesh_, fec_.get(), 3, Ordering::byNODES);
      // Scalar DG space for BR2 mass matrix inverses
      scalar_fes_ = std::make_unique<FiniteElementSpace>(&mesh_, fec_.get());
   }

   void SetupBoundaryMarkers()
   {
      int num_bdr = mesh_.bdr_attributes.Size() > 0 ? mesh_.bdr_attributes.Max() : 0;
      dirichlet_bdr_marker_.SetSize(num_bdr);
      dirichlet_bdr_marker_ = 0;

      // Always use BoundaryConfig (populated by both constructors)
      for (int attr : bdr_config_.dirichlet_attrs)
      {
         MFEM_VERIFY(attr >= 1 && attr <= num_bdr,
                     "BoundaryConfig: Dirichlet attr " << attr
                     << " not in mesh (max attr = " << num_bdr << ")");
         dirichlet_bdr_marker_[attr - 1] = 1;
      }
   }

   /// Build the facet BC tables from mesh tags. Single source of truth.
   ///
   /// Mirrors Tandem's facet-BC model: import tags → propagate exactly →
   /// store per-facet class → consume per-facet class.
   ///
   /// 1. Scan boundary elements with attr 3 (fault) and attr 5 (Dirichlet)
   /// 2. Mark interior faces directly in face_bc_
   /// 3. Collect shared-face canonical vertex keys, Allgatherv across ranks
   /// 4. Mark shared faces in shared_face_bc_
   /// 5. Validate: attr 3 and attr 5 must exist; fault ∩ Dirichlet = ∅
   /// 6. Derive legacy face arrays for downstream assembly
   void BuildFacetBCTables()
   {
      const int num_faces = mesh_.GetNumFaces();
      face_bc_.assign(num_faces, FacetBC::None);
      face_bc_attr_.assign(num_faces, 0);

      int num_shared = 0;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         num_shared = mesh_.GetNSharedFaces();
#endif
      }
      shared_face_bc_.assign(num_shared, FacetBC::None);
      shared_face_bc_attr_.assign(num_shared, 0);

      // ---- Resolve BC attr numbers ----
      const int fault_attr = bdr_config_.fault_attr;
      const auto &dir_attrs = bdr_config_.dirichlet_attrs;

      // ---- Check global attr existence ----
      int local_has_fault = 0, local_has_dir = 0;
      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (attr == fault_attr) { local_has_fault = 1; }
         if (dir_attrs.count(attr)) { local_has_dir = 1; }
      }
      int global_has_fault = local_has_fault;
      int global_has_dir = local_has_dir;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         MPI_Allreduce(MPI_IN_PLACE, &global_has_fault, 1, MPI_INT,
                       MPI_MAX, mesh_.GetComm());
         MPI_Allreduce(MPI_IN_PLACE, &global_has_dir, 1, MPI_INT,
                       MPI_MAX, mesh_.GetComm());
#endif
      }
      bool has_fault = (global_has_fault > 0);
      bool has_dirichlet = (global_has_dir > 0);
      if (has_fault)
      {
         MFEM_VERIFY(has_dirichlet,
            "ERROR: Mesh has fault faces (attr " << fault_attr
            << ") but no Dirichlet faces. Requires both.");
      }
      if (!has_fault && !has_dirichlet) { return; }

      // ---- Build reverse map: local face → shared face index ----
      std::unordered_map<int, int> lface_to_sface;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int sf = 0; sf < num_shared; sf++)
         {
            lface_to_sface[mesh_.GetSharedFace(sf)] = sf;
         }
#endif
      }

      // ---- Phase 1: Scan boundary elements, mark local faces ----
      // Collect shared-face tags locally (only 1 rank has the bdr element)
      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         bool is_fault = (attr == fault_attr);
         bool is_dir = dir_attrs.count(attr) > 0;
         if (!is_fault && !is_dir) { continue; }
         FacetBC bc = is_fault ? FacetBC::Fault : FacetBC::Dirichlet;
         int face_idx = mesh_.GetBdrElementFaceIndex(be);

         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face_idx);
         if (FTr != nullptr)
         {
            // Interior face: mark in face_bc_ table
            MFEM_VERIFY(face_bc_[face_idx] == FacetBC::None ||
                        face_bc_[face_idx] == bc,
               "ERROR: Interior face " << face_idx << " tagged as both "
               "fault (attr " << fault_attr << ") and Dirichlet.");
            face_bc_[face_idx] = bc;
            face_bc_attr_[face_idx] = attr;

            if (bc == FacetBC::Fault)
            {
               int e1 = FTr->Elem1No;
               int e2 = FTr->Elem2No;
               long key = (long)std::min(e1, e2) * mesh_.GetNE()
                          + std::max(e1, e2);
               fault_face_keys_.insert(key);
            }
         }
         else if constexpr (IsParallelMesh<MeshType>::value)
         {
            auto it = lface_to_sface.find(face_idx);
            if (it != lface_to_sface.end())
            {
               int sf = it->second;
               MFEM_VERIFY(shared_face_bc_[sf] == FacetBC::None ||
                           shared_face_bc_[sf] == bc,
                  "ERROR: Shared face " << sf << " tagged as both "
                  "attr 3 and attr 5.");
               shared_face_bc_[sf] = bc;
               shared_face_bc_attr_[sf] = attr;
            }
         }
      }

      // ---- Phase 2: Propagate shared-face tags across MPI ranks ----
      // A shared face's boundary element exists on only ONE rank. The
      // other rank must discover it via canonical vertex key exchange.
      // We propagate fault and Dirichlet tags in a single Allgatherv,
      // packing (key, bc_class) per face.
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         Array<HYPRE_BigInt> gvert;
         mesh_.GetGlobalVertexIndices(gvert);

         // Pack: 5 values per locally-tagged shared face (3 key + bc_class + attr)
         std::vector<HYPRE_BigInt> local_flat;
         for (int sf = 0; sf < num_shared; sf++)
         {
            if (shared_face_bc_[sf] == FacetBC::None) { continue; }
            int lf = mesh_.GetSharedFace(sf);
            FaceVertexKey fk = MakeFaceKey(lf, gvert);
            local_flat.push_back(fk.v[0]);
            local_flat.push_back(fk.v[1]);
            local_flat.push_back(fk.v[2]);
            local_flat.push_back(static_cast<HYPRE_BigInt>(shared_face_bc_[sf]));
            local_flat.push_back(static_cast<HYPRE_BigInt>(shared_face_bc_attr_[sf]));
         }

         int lc = static_cast<int>(local_flat.size());
         int nranks = 1;
         MPI_Comm_size(mesh_.GetComm(), &nranks);
         std::vector<int> rc(nranks), dp(nranks);
         MPI_Allgather(&lc, 1, MPI_INT, rc.data(), 1, MPI_INT,
                       mesh_.GetComm());
         int tot = 0;
         for (int r = 0; r < nranks; r++) { dp[r] = tot; tot += rc[r]; }
         std::vector<HYPRE_BigInt> all(tot);
         MPI_Allgatherv(local_flat.data(), lc, HYPRE_MPI_BIG_INT,
                        all.data(), rc.data(), dp.data(),
                        HYPRE_MPI_BIG_INT, mesh_.GetComm());

         // Build global key → (bc, attr) map
         std::map<FaceVertexKey, std::pair<FacetBC, int>> global_shared_bc;
         for (int i = 0; i < tot; i += 5)
         {
            FaceVertexKey k;
            k.v[0] = all[i]; k.v[1] = all[i+1]; k.v[2] = all[i+2];
            FacetBC bc = static_cast<FacetBC>(all[i+3]);
            int bc_attr = static_cast<int>(all[i+4]);
            auto [it, inserted] = global_shared_bc.emplace(
               k, std::make_pair(bc, bc_attr));
            MFEM_VERIFY(inserted || (it->second.first == bc &&
                                     it->second.second == bc_attr),
               "ERROR: Shared face key (" << k.v[0] << "," << k.v[1]
               << "," << k.v[2] << ") has conflicting tags from "
               "different ranks (bc=" << static_cast<int>(bc)
               << " attr=" << bc_attr << ").");
         }

         // Fill untagged local shared faces from the global map
         for (int sf = 0; sf < num_shared; sf++)
         {
            if (shared_face_bc_[sf] != FacetBC::None) { continue; }
            int lf = mesh_.GetSharedFace(sf);
            FaceVertexKey key = MakeFaceKey(lf, gvert);
            auto it = global_shared_bc.find(key);
            if (it != global_shared_bc.end())
            {
               shared_face_bc_[sf] = it->second.first;
               shared_face_bc_attr_[sf] = it->second.second;
            }
         }
#endif
      }

      // ---- Phase 3: Derive legacy face arrays from BC tables ----
      fault_interior_faces_.SetSize(0);
      fault_shared_faces_.SetSize(0);
      dirichlet_interior_faces_.SetSize(0);
      dirichlet_shared_faces_.SetSize(0);
      dirichlet_interior_attrs_.SetSize(0);
      dirichlet_shared_attrs_.SetSize(0);

      for (int f = 0; f < num_faces; f++)
      {
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(f);
         if (!FTr) { continue; }
         if (face_bc_[f] == FacetBC::Fault)
         {
            fault_interior_faces_.Append(f);
         }
         else if (face_bc_[f] == FacetBC::Dirichlet)
         {
            dirichlet_interior_faces_.Append(f);
            dirichlet_interior_attrs_.Append(face_bc_attr_[f]);
         }
      }
      for (int sf = 0; sf < num_shared; sf++)
      {
         if (shared_face_bc_[sf] == FacetBC::Fault)
         {
            fault_shared_faces_.Append(sf);
         }
         else if (shared_face_bc_[sf] == FacetBC::Dirichlet)
         {
            dirichlet_shared_faces_.Append(sf);
            dirichlet_shared_attrs_.Append(shared_face_bc_attr_[sf]);
         }
      }

      // ---- Phase 4: Assert no None-classified face in any assembly array ----
      for (int i = 0; i < fault_interior_faces_.Size(); i++)
      {
         MFEM_VERIFY(face_bc_[fault_interior_faces_[i]] == FacetBC::Fault,
            "ERROR: fault_interior_faces_[" << i << "] = "
            << fault_interior_faces_[i] << " has BC="
            << static_cast<int>(face_bc_[fault_interior_faces_[i]])
            << ", expected Fault.");
      }
      for (int i = 0; i < dirichlet_interior_faces_.Size(); i++)
      {
         MFEM_VERIFY(face_bc_[dirichlet_interior_faces_[i]] == FacetBC::Dirichlet,
            "ERROR: dirichlet_interior_faces_[" << i << "] = "
            << dirichlet_interior_faces_[i] << " has BC="
            << static_cast<int>(face_bc_[dirichlet_interior_faces_[i]])
            << ", expected Dirichlet.");
      }
      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         MFEM_VERIFY(shared_face_bc_[fault_shared_faces_[i]] == FacetBC::Fault,
            "ERROR: fault_shared_faces_[" << i << "] = "
            << fault_shared_faces_[i] << " has BC="
            << static_cast<int>(shared_face_bc_[fault_shared_faces_[i]])
            << ", expected Fault.");
      }
      for (int i = 0; i < dirichlet_shared_faces_.Size(); i++)
      {
         MFEM_VERIFY(shared_face_bc_[dirichlet_shared_faces_[i]] == FacetBC::Dirichlet,
            "ERROR: dirichlet_shared_faces_[" << i << "] = "
            << dirichlet_shared_faces_[i] << " has BC="
            << static_cast<int>(shared_face_bc_[dirichlet_shared_faces_[i]])
            << ", expected Dirichlet.");
      }

      // ---- Phase 4b: Audit y=0 face classification ----
      // Detect interior faces on y=0 that are FacetBC::None (unclassified).
      // These faces get DG penalty enforcing zero jump, conflicting with
      // adjacent fault slip or Dirichlet displacement.
      {
         int rank = 0;
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Comm_rank(mesh_.GetComm(), &rank);
#endif
         }

         int local_y0_none = 0, local_y0_fault = 0, local_y0_dir = 0;
         int local_y0_none_shared = 0;
         std::vector<std::array<double, 4>> none_faces;  // x, y, z, face_idx

         // Interior faces
         for (int f = 0; f < num_faces; f++)
         {
            auto *FTr = mesh_.GetInteriorFaceTransformations(f);
            if (!FTr) { continue; }
            const IntegrationPoint &ip =
               Geometries.GetCenter(FTr->GetGeometryType());
            FTr->SetAllIntPoints(&ip);
            Vector fc(3);
            FTr->Face->Transform(ip, fc);
            if (std::abs(fc(1)) > 1.0) { continue; }  // not on y=0

            if (face_bc_[f] == FacetBC::None)
            {
               local_y0_none++;
               none_faces.push_back({fc(0), fc(1), fc(2),
                                     static_cast<double>(f)});
            }
            else if (face_bc_[f] == FacetBC::Fault) { local_y0_fault++; }
            else if (face_bc_[f] == FacetBC::Dirichlet) { local_y0_dir++; }
         }

         // Shared faces — count by classification only (avoid
         // GetSharedFaceTransformations which can fail for boundary faces).
         // Use GetSharedFace to get the local face index, then compute
         // centroid via the local face geometry.
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            for (int sf = 0; sf < num_shared; sf++)
            {
               if (shared_face_bc_[sf] != FacetBC::None) { continue; }

               // Get centroid of shared face via local face index
               int lf = mesh_.GetSharedFace(sf);
               auto *face_tr = mesh_.GetFaceTransformation(lf);
               if (!face_tr) { local_y0_none_shared++; continue; }
               const IntegrationPoint &ip =
                  Geometries.GetCenter(face_tr->GetGeometryType());
               face_tr->SetIntPoint(&ip);
               Vector fc(3);
               face_tr->Transform(ip, fc);
               if (std::abs(fc(1)) > 1.0) { continue; }

               local_y0_none_shared++;
               none_faces.push_back({fc(0), fc(1), fc(2),
                                     static_cast<double>(lf)});
            }
#endif
         }

         // Reduce counts
         int global_counts[4] = {local_y0_fault, local_y0_dir,
                                  local_y0_none, local_y0_none_shared};
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            MPI_Allreduce(MPI_IN_PLACE, global_counts, 4, MPI_INT,
                          MPI_SUM, mesh_.GetComm());
#endif
         }

         if (rank == 0)
         {
            std::cout << "\n=== Y=0 Face Classification Audit ===\n"
                      << "  Fault (attr=3):      " << global_counts[0] << "\n"
                      << "  Dirichlet (attr=5):  " << global_counts[1] << "\n"
                      << "  NONE (interior):     " << global_counts[2] << "\n"
                      << "  NONE (shared):       " << global_counts[3] << "\n";
            if (global_counts[2] + global_counts[3] > 0)
            {
               std::cout << "  *** WARNING: " << global_counts[2] + global_counts[3]
                         << " y=0 faces have NO BC! These enforce zero jump "
                         << "via DG penalty, conflicting with fault slip.\n";
            }
            else
            {
               std::cout << "  OK: all y=0 interior faces classified.\n";
            }
            std::cout << "=====================================\n\n";
         }

         // Dump unclassified face locations to CSV
         if (!none_faces.empty())
         {
            std::string fname = "face_audit_r" + std::to_string(rank) + ".csv";
            std::ofstream ofs(fname);
            ofs << "x,y,z,face_idx,bc\n";
            for (auto &nf : none_faces)
            {
               ofs << nf[0] << "," << nf[1] << ","
                   << nf[2] << "," << static_cast<int>(nf[3])
                   << ",None\n";
            }
            // Also dump fault faces near the boundary for context
            for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
            {
               int f = fault_interior_faces_[fi];
               auto *FTr = mesh_.GetInteriorFaceTransformations(f);
               if (!FTr) { continue; }
               const IntegrationPoint &ip =
                  Geometries.GetCenter(FTr->GetGeometryType());
               FTr->SetAllIntPoints(&ip);
               Vector fc(3);
               FTr->Face->Transform(ip, fc);
               ofs << fc(0) << "," << fc(1) << ","
                   << fc(2) << "," << f << ",Fault\n";
            }
            for (int di = 0; di < dirichlet_interior_faces_.Size(); di++)
            {
               int f = dirichlet_interior_faces_[di];
               auto *FTr = mesh_.GetInteriorFaceTransformations(f);
               if (!FTr) { continue; }
               const IntegrationPoint &ip =
                  Geometries.GetCenter(FTr->GetGeometryType());
               FTr->SetAllIntPoints(&ip);
               Vector fc(3);
               FTr->Face->Transform(ip, fc);
               ofs << fc(0) << "," << fc(1) << ","
                   << fc(2) << "," << f << ",Dirichlet\n";
            }
         }
      }

      // ---- Phase 5: Exact canonical-key validation ----
      ValidateFacetBCTables();
   }

   /// Allgather a local set of FaceVertexKeys and return the global union.
   std::set<FaceVertexKey> AllgatherKeys(
      const std::set<FaceVertexKey> &local_keys) const
   {
      std::set<FaceVertexKey> global_keys = local_keys;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int lc = static_cast<int>(local_keys.size());
         int nranks = 1;
         MPI_Comm_size(mesh_.GetComm(), &nranks);
         std::vector<int> rc(nranks), dp(nranks);
         MPI_Allgather(&lc, 1, MPI_INT, rc.data(), 1, MPI_INT,
                       mesh_.GetComm());
         int tot = 0;
         for (int r = 0; r < nranks; r++) { dp[r] = tot; tot += rc[r]; }
         std::vector<HYPRE_BigInt> lf(3 * lc);
         int idx = 0;
         for (auto &k : local_keys)
         {
            lf[3*idx] = k.v[0]; lf[3*idx+1] = k.v[1]; lf[3*idx+2] = k.v[2];
            idx++;
         }
         std::vector<int> rc3(nranks), dp3(nranks);
         for (int r = 0; r < nranks; r++)
         { rc3[r] = 3*rc[r]; dp3[r] = 3*dp[r]; }
         std::vector<HYPRE_BigInt> af(3 * tot);
         MPI_Allgatherv(lf.data(), 3*lc, HYPRE_MPI_BIG_INT,
                        af.data(), rc3.data(), dp3.data(),
                        HYPRE_MPI_BIG_INT, mesh_.GetComm());
         for (int i = 0; i < tot; i++)
         {
            FaceVertexKey k;
            k.v[0] = af[3*i]; k.v[1] = af[3*i+1]; k.v[2] = af[3*i+2];
            global_keys.insert(k);
         }
#endif
      }
      return global_keys;
   }

   /// Exact canonical-key validation of recovered face sets.
   void ValidateFacetBCTables() const
   {
      Array<HYPRE_BigInt> gvert;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         mesh_.GetGlobalVertexIndices(gvert);
#endif
      }
      else
      {
         gvert.SetSize(mesh_.GetNV());
         for (int i = 0; i < mesh_.GetNV(); i++) { gvert[i] = i; }
      }

      // Build tagged key sets from boundary elements (ground truth).
      const int fault_attr = bdr_config_.fault_attr;
      const auto &dir_attrs = bdr_config_.dirichlet_attrs;

      std::set<FaceVertexKey> local_tagged_fault;
      std::set<FaceVertexKey> local_tagged_dir_recoverable;  // interior/shared only

      // R-006: build a single lface_to_sface map once instead of
      // scanning shared faces inside the per-boundary-element loop.
      std::unordered_map<int, int> lf2sf;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int sf = 0; sf < mesh_.GetNSharedFaces(); sf++)
         {
            lf2sf[mesh_.GetSharedFace(sf)] = sf;
         }
#endif
      }

      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         int attr = mesh_.GetBdrAttribute(be);
         if (attr != fault_attr && !dir_attrs.count(attr)) { continue; }
         int face_idx = mesh_.GetBdrElementFaceIndex(be);
         FaceVertexKey key = MakeFaceKey(face_idx, gvert);

         if (attr == fault_attr)
         {
            // Fault attrs MUST be on 2-sided interior or shared faces
            // (split-mesh convention).  A 1-sided outer-boundary fault
            // tag is almost always a mesh-construction bug — keep the
            // strict tagging so the bidirectional equality check
            // surfaces the issue loudly.  Tests that need an unused
            // outer attribute should use a non-conflicting number
            // (e.g., 7), not fault_attr (REVIEW R-001).
            local_tagged_fault.insert(key);
         }
         else
         {
            // Dirichlet attrs can legitimately be on 1-sided outer
            // faces (far-field BC on +Lx, etc.).  Only the interior /
            // shared subset is "recoverable" through the face_bc_
            // machinery; 1-sided cases are handled elsewhere.  Filter
            // to that subset before adding to the recoverable set.
            const bool is_interior =
               (mesh_.GetInteriorFaceTransformations(face_idx) != nullptr);
            bool is_shared = false;
            if constexpr (IsParallelMesh<MeshType>::value)
            {
#ifdef MFEM_USE_MPI
               if (!is_interior && lf2sf.find(face_idx) != lf2sf.end())
               {
                  is_shared = true;
               }
#endif
            }
            if (!is_interior && !is_shared) { continue; }
            local_tagged_dir_recoverable.insert(key);
         }
      }
      std::set<FaceVertexKey> global_tagged_fault =
         AllgatherKeys(local_tagged_fault);
      std::set<FaceVertexKey> global_tagged_dir_recoverable =
         AllgatherKeys(local_tagged_dir_recoverable);

      // Build recovered key sets from the derived face arrays
      std::set<FaceVertexKey> local_recovered_fault, local_recovered_dir;
      for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
         local_recovered_fault.insert(MakeFaceKey(fault_interior_faces_[fi], gvert));
      for (int fi = 0; fi < dirichlet_interior_faces_.Size(); fi++)
         local_recovered_dir.insert(MakeFaceKey(dirichlet_interior_faces_[fi], gvert));
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
            local_recovered_fault.insert(
               MakeFaceKey(mesh_.GetSharedFace(fault_shared_faces_[i]), gvert));
         for (int i = 0; i < dirichlet_shared_faces_.Size(); i++)
            local_recovered_dir.insert(
               MakeFaceKey(mesh_.GetSharedFace(dirichlet_shared_faces_[i]), gvert));
#endif
      }
      std::set<FaceVertexKey> global_recovered_fault =
         AllgatherKeys(local_recovered_fault);
      std::set<FaceVertexKey> global_recovered_dir =
         AllgatherKeys(local_recovered_dir);

      // Exact fault key equality: tagged == recovered (bidirectional)
      for (auto &k : global_tagged_fault)
      {
         MFEM_VERIFY(global_recovered_fault.count(k) > 0,
            "ERROR: Tagged attr-3 face key (" << k.v[0] << "," << k.v[1]
            << "," << k.v[2] << ") was not recovered as a fault face.");
      }
      for (auto &k : global_recovered_fault)
      {
         MFEM_VERIFY(global_tagged_fault.count(k) > 0,
            "ERROR: Recovered fault face key (" << k.v[0] << "," << k.v[1]
            << "," << k.v[2] << ") has no matching attr-3 boundary element.");
      }

      // Exact Dirichlet key equality for recoverable faces (bidirectional)
      for (auto &k : global_tagged_dir_recoverable)
      {
         MFEM_VERIFY(global_recovered_dir.count(k) > 0,
            "ERROR: Tagged attr-5 interior/shared face key (" << k.v[0]
            << "," << k.v[1] << "," << k.v[2]
            << ") was not recovered as a Dirichlet face.");
      }
      for (auto &k : global_recovered_dir)
      {
         MFEM_VERIFY(global_tagged_dir_recoverable.count(k) > 0,
            "ERROR: Recovered Dirichlet face key (" << k.v[0] << ","
            << k.v[1] << "," << k.v[2]
            << ") has no matching recoverable attr-5 boundary element.");
      }

      // Fault ∩ Dirichlet = ∅
      for (auto &k : global_recovered_fault)
      {
         MFEM_VERIFY(!global_recovered_dir.count(k),
            "ERROR: Face key (" << k.v[0] << "," << k.v[1] << ","
            << k.v[2] << ") classified as both fault AND Dirichlet.");
      }

      // Summary
      bool is_root = true;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int rank; MPI_Comm_rank(mesh_.GetComm(), &rank);
         is_root = (rank == 0);
#endif
      }
      if (is_root)
      {
         mfem::out << "  Tag validation: "
                   << global_tagged_fault.size() << " fault keys (exact), "
                   << global_recovered_dir.size() << " Dirichlet keys (exact, "
                   << global_tagged_dir_recoverable.size() << " recoverable)"
                   << " (OK)\n";
         mfem::out << "  Fault faces: "
                   << fault_interior_faces_.Size() << " interior + "
                   << fault_shared_faces_.Size() << " shared (this rank)\n";
         mfem::out << "  Dirichlet faces: "
                   << dirichlet_interior_faces_.Size() << " interior + "
                   << dirichlet_shared_faces_.Size() << " shared (this rank)\n";
      }
   }

   void SetupFaultInfo()
   {
      // Build facet BC tables (single source of truth) and derive
      // legacy face arrays. Includes exact canonical-key validation.
      BuildFacetBCTables();

      num_fault_faces_ = fault_interior_faces_.Size() + fault_shared_faces_.Size();

      // Exchange face-neighbor data for shared face assembly
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         mesh_.ExchangeFaceNbrData();
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         if (pfes) { pfes->ExchangeFaceNbrData(); }
#endif
      }

      // Multi-DOF fault quadrature
      // IP: use the same nodal triangle order as the volume space, matching
      // Tandem's fault discretization even at p=1.
      // BR2: keep the legacy face-averaged path (nbf=1).
      int face_fe_order = (method_ == DGMethod::IP) ? order_ : 0;
      face_quad_ = std::make_unique<FaceQuadrature>(face_fe_order,
                                                     std::max(order_, 1),
                                                     Geometry::TRIANGLE,
                                                     face_basis_type_);
      nbf_per_face_ = face_quad_->NumBasisFunctions();
      num_fault_dofs_ = num_fault_faces_ * nbf_per_face_;
      BuildOwnedFaultLayout();

      // Diagnostic: report owned vs local fault DOFs across all ranks
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int rank = 0, nranks = 1;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         MPI_Comm_size(mesh_.GetComm(), &nranks);
         int local_shared = fault_shared_faces_.Size();
         int global_shared = 0;
         MPI_Reduce(&local_shared, &global_shared, 1, MPI_INT, MPI_SUM, 0,
                    mesh_.GetComm());
         int global_owned = 0;
         MPI_Reduce(&num_owned_fault_dofs_, &global_owned, 1, MPI_INT, MPI_SUM,
                    0, mesh_.GetComm());
         int global_local = 0;
         MPI_Reduce(&num_fault_dofs_, &global_local, 1, MPI_INT, MPI_SUM,
                    0, mesh_.GetComm());
         // Check how many faces have non-identity permutation
         int local_permuted = 0;
         for (int fi = 0; fi < num_fault_faces_; fi++)
         {
            for (int k = 0; k < nbf_per_face_; k++)
            {
               if (canonical_to_local_perm_[fi * nbf_per_face_ + k] != k)
               {
                  local_permuted++;
                  break;
               }
            }
         }
         int global_permuted = 0;
         MPI_Reduce(&local_permuted, &global_permuted, 1, MPI_INT, MPI_SUM, 0,
                    mesh_.GetComm());
         if (rank == 0)
         {
            mfem::out << "  Owned fault layout: "
                      << "global_shared_faces=" << global_shared
                      << ", global_owned_dofs=" << global_owned
                      << ", global_local_dofs=" << global_local
                      << ", faces_with_nontrivial_perm=" << global_permuted
                      << "\n";
         }
#endif
      }

      fault_dofs_.SetSize(num_fault_dofs_);
      for (int i = 0; i < num_fault_dofs_; i++)
      {
         fault_dofs_[i] = i;
      }

      // ref_normal_ is needed by ComputeSkeletonDirichletSign for ALL
      // ranks that have shared Dirichlet faces, even ranks with zero fault
      // DOFs. Set it unconditionally BEFORE the fault DOF guard.
      {
         Vector ref_normal(3);
         ref_normal = 0.0;
         ref_normal(1) = -1.0;  // Y = fault-normal, pointing -Y
         ref_normal_ = ref_normal;
      }

      // Compute FaultBasis for coordinate transforms
      if (num_fault_dofs_ > 0)
      {
         // Tandem coordinate system:
         //   X = along-strike, Y = fault-normal, Z = depth (negative down)
         //   Fault at Y = 0
         //   ref_normal = (0, -1, 0) matches Tandem's convention
         //   Up = (0, 0, 1)
         const Vector &ref_normal = ref_normal_;

         Vector up(3);
         up = 0.0;
         up(2) = 1.0;  // +Z = upward

         fault_basis_.Compute(mesh_, fault_interior_faces_, ref_normal, up);

         // Append shared faces to FaultBasis
         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            fault_basis_.AppendSharedFaces(mesh_, fault_shared_faces_,
                                           ref_normal, up);
#endif
         }

         // Per-quad-point basis (matching Tandem AdapterBase::prepare).
         // Uses the same 2p+1 face quadrature order as the combined integrator.
         // Detect face geometry from the first available fault face (interior
         // or shared), so ranks with only shared faces get the right rule.
         const int face_quad_order = 2 * order_ + 1;
         Geometry::Type face_geom = Geometry::TRIANGLE;  // default for tets
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
               auto *ftr0 = mesh_.GetSharedFaceTransformations(
                  fault_shared_faces_[0]);
               if (ftr0) { face_geom = ftr0->GetGeometryType(); }
            }
#endif
         }
         const IntegrationRule &face_ir =
            IntRules.Get(face_geom, face_quad_order);

         fault_basis_.ComputeQPBasis(mesh_, fault_interior_faces_,
                                      ref_normal, up, face_ir);

         if constexpr (IsParallelMesh<MeshType>::value)
         {
#ifdef MFEM_USE_MPI
            fault_basis_.ComputeQPBasisShared(mesh_, fault_shared_faces_,
                                               ref_normal, up, face_ir,
                                               fault_interior_faces_.Size());
#endif
         }
      }
   }

   /// Startup face audit (validation-only). Verifies shared-face classification
   /// agreement between neighboring MPI ranks. No coordinate-based logic —
   /// only validates that the tag-recovered sets are self-consistent.
   void RunStartupFaceAudit()
   {
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         int rank = 0, nranks = 1;
         MPI_Comm_rank(mesh_.GetComm(), &rank);
         MPI_Comm_size(mesh_.GetComm(), &nranks);

         // Build classification lookup for shared faces:
         // 0=unclassified, 1=fault, 2=dirichlet
         std::map<int, int> shared_cls;
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
            shared_cls[fault_shared_faces_[i]] = 1;
         for (int i = 0; i < dirichlet_shared_faces_.Size(); i++)
            shared_cls[dirichlet_shared_faces_[i]] = 2;

         // Pack: 4 values per classified shared face (3 key + cls)
         Array<HYPRE_BigInt> gvert;
         mesh_.GetGlobalVertexIndices(gvert);

         std::vector<HYPRE_BigInt> local_flat;
         for (auto &[sf, cls] : shared_cls)
         {
            int lf = mesh_.GetSharedFace(sf);
            FaceVertexKey fk = MakeFaceKey(lf, gvert);
            local_flat.push_back(fk.v[0]);
            local_flat.push_back(fk.v[1]);
            local_flat.push_back(fk.v[2]);
            local_flat.push_back(cls);
         }

         int lc = static_cast<int>(local_flat.size());
         std::vector<int> counts(nranks), displs(nranks);
         MPI_Allgather(&lc, 1, MPI_INT, counts.data(), 1, MPI_INT,
                        mesh_.GetComm());
         int total = 0;
         for (int r = 0; r < nranks; r++)
         {
            displs[r] = total;
            total += counts[r];
         }
         std::vector<HYPRE_BigInt> all(total);
         MPI_Allgatherv(local_flat.data(), lc, HYPRE_MPI_BIG_INT,
                         all.data(), counts.data(), displs.data(),
                         HYPRE_MPI_BIG_INT, mesh_.GetComm());

         // Check: each key must appear exactly twice (once per rank sharing
         // the face), and both ranks must agree on the classification.
         std::map<FaceVertexKey, std::vector<int>> kmap;
         for (int i = 0; i < total; i += 4)
         {
            FaceVertexKey k;
            k.v[0] = all[i]; k.v[1] = all[i+1]; k.v[2] = all[i+2];
            kmap[k].push_back(static_cast<int>(all[i+3]));
         }
         int mismatch = 0, single_rank = 0;
         for (auto &[k, entries] : kmap)
         {
            if (entries.size() < 2)
            {
               // Face seen from only one rank — valid in some MFEM
               // partitioning layouts where only the owning rank
               // contributes the shared face to the allgather.
               single_rank++;
               continue;
            }
            // Check that all ranks contributing this face agree
            for (size_t e = 1; e < entries.size(); e++)
            {
               if (entries[e] != entries[0])
               {
                  mismatch++;
                  if (rank == 0 && mismatch <= 10)
                  {
                     mfem::out << "  AUDIT MISMATCH: key=("
                               << k.v[0] << "," << k.v[1] << "," << k.v[2]
                               << ") cls=" << entries[0] << " vs cls="
                               << entries[e] << "\n";
                  }
                  break;
               }
            }
         }
         MFEM_VERIFY(mismatch == 0,
            "ERROR: " << mismatch << " shared faces have classification "
            "mismatch between neighboring ranks.");
         if (rank == 0)
         {
            int paired = static_cast<int>(kmap.size()) - single_rank;
            mfem::out << "  AUDIT: shared-face agreement OK ("
                      << paired << " paired, "
                      << single_rank << " single-rank)\n";
         }
#endif
      }
   }

   void BuildOwnedFaultLayout()
   {
      owned_fault_face_to_local_face_.SetSize(0);
      full_fault_face_to_owned_face_.SetSize(num_fault_faces_);
      full_fault_face_to_owned_face_ = -1;
      owned_fault_dof_to_local_dof_.SetSize(0);
      shared_fault_comm_blocks_.clear();

      const int num_interior = fault_interior_faces_.Size();

      for (int face_idx = 0; face_idx < num_interior; face_idx++)
      {
         full_fault_face_to_owned_face_[face_idx] =
            owned_fault_face_to_local_face_.Size();
         owned_fault_face_to_local_face_.Append(face_idx);
      }

      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         Array<HYPRE_BigInt> gvert;
         mesh_.GetGlobalVertexIndices(gvert);

         int rank = 0;
         MPI_Comm_rank(mesh_.GetComm(), &rank);

         // Build canonical vertex-key for each local shared fault face
         std::vector<FaceVertexKey> local_shared_keys(fault_shared_faces_.Size());
         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            const int sf = fault_shared_faces_[i];
            const int local_face = mesh_.GetSharedFace(sf);
            local_shared_keys[i] = MakeFaceKey(local_face, gvert);
         }

         // Allgather: collect all shared fault face keys + originating rank
         int comm_size = 1;
         MPI_Comm_size(mesh_.GetComm(), &comm_size);
         const int local_shared_count = fault_shared_faces_.Size();
         std::vector<int> recv_counts(comm_size, 0), recv_displs(comm_size, 0);
         MPI_Allgather(&local_shared_count, 1, MPI_INT,
                       recv_counts.data(), 1, MPI_INT, mesh_.GetComm());

         int total_shared = 0;
         for (int r = 0; r < comm_size; r++)
         {
            recv_displs[r] = total_shared;
            total_shared += recv_counts[r];
         }

         // Pack as flat array: 3 HYPRE_BigInt per face (sorted vertex IDs)
         std::vector<HYPRE_BigInt> local_flat(3 * local_shared_count);
         for (int i = 0; i < local_shared_count; i++)
         {
            local_flat[3*i]   = local_shared_keys[i].v[0];
            local_flat[3*i+1] = local_shared_keys[i].v[1];
            local_flat[3*i+2] = local_shared_keys[i].v[2];
         }
         std::vector<int> recv3(comm_size), disp3(comm_size);
         for (int r = 0; r < comm_size; r++)
         {
            recv3[r] = 3 * recv_counts[r];
            disp3[r] = 3 * recv_displs[r];
         }
         std::vector<HYPRE_BigInt> all_flat(3 * total_shared);
         MPI_Allgatherv(local_flat.data(), 3 * local_shared_count,
                        HYPRE_MPI_BIG_INT,
                        all_flat.data(), recv3.data(), disp3.data(),
                        HYPRE_MPI_BIG_INT, mesh_.GetComm());

         // Build map: face key → list of ranks that have it
         std::map<FaceVertexKey, std::vector<int>> key_ranks;
         for (int r = 0; r < comm_size; r++)
         {
            for (int j = 0; j < recv_counts[r]; j++)
            {
               int idx = recv_displs[r] + j;
               FaceVertexKey k;
               k.v[0] = all_flat[3*idx];
               k.v[1] = all_flat[3*idx+1];
               k.v[2] = all_flat[3*idx+2];
               key_ranks[k].push_back(r);
            }
         }
         for (auto &kv : key_ranks)
         {
            auto &ranks = kv.second;
            std::sort(ranks.begin(), ranks.end());
            ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
         }

         std::unordered_map<int, std::vector<SharedFaultFaceBlock>> send_by_rank;
         std::unordered_map<int, std::vector<SharedFaultFaceBlock>> recv_by_rank;

         for (int i = 0; i < fault_shared_faces_.Size(); i++)
         {
            const int sf = fault_shared_faces_[i];
            const int local_face = mesh_.GetSharedFace(sf);
            const int face_idx = num_interior + i;
            const FaceVertexKey key = local_shared_keys[i];
            auto key_it = key_ranks.find(key);
            MFEM_VERIFY(key_it != key_ranks.end() && !key_it->second.empty(),
                        "Missing shared fault ownership info for face key");
            const auto &sharing_ranks = key_it->second;
            const int owner_rank = sharing_ranks.front();

            if (rank == owner_rank)
            {
               const int owned_face =
                  owned_fault_face_to_local_face_.Size();
               full_fault_face_to_owned_face_[face_idx] = owned_face;
               owned_fault_face_to_local_face_.Append(face_idx);
               for (int other_rank : sharing_ranks)
               {
                  if (other_rank == rank) { continue; }
                  send_by_rank[other_rank].push_back({key, owned_face});
               }
            }
            else
            {
               recv_by_rank[owner_rank].push_back({key, face_idx});
            }
         }

         std::set<int> neighbors;
         for (const auto &kv : send_by_rank) { neighbors.insert(kv.first); }
         for (const auto &kv : recv_by_rank) { neighbors.insert(kv.first); }

         for (int neighbor : neighbors)
         {
            auto &send_faces = send_by_rank[neighbor];
            auto &recv_faces = recv_by_rank[neighbor];

            std::sort(send_faces.begin(), send_faces.end(),
                      [](const SharedFaultFaceBlock &a,
                         const SharedFaultFaceBlock &b)
                      {
                         return a.key < b.key;
                      });
            std::sort(recv_faces.begin(), recv_faces.end(),
                      [](const SharedFaultFaceBlock &a,
                         const SharedFaultFaceBlock &b)
                      {
                         return a.key < b.key;
                      });

            SharedFaultCommBlock block;
            block.neighbor_rank = neighbor;
            block.send_owned_faces.reserve(send_faces.size());
            block.recv_local_faces.reserve(recv_faces.size());
            for (const auto &entry : send_faces)
            {
               block.send_owned_faces.push_back(entry.face_idx);
            }
            for (const auto &entry : recv_faces)
            {
               block.recv_local_faces.push_back(entry.face_idx);
            }
            shared_fault_comm_blocks_.push_back(std::move(block));
         }
#endif
      }

      num_owned_fault_faces_ = owned_fault_face_to_local_face_.Size();
      num_owned_fault_dofs_ = num_owned_fault_faces_ * nbf_per_face_;

      // ------------------------------------------------------------------
      // Build per-face canonical DOF permutation (Tandem sorted-simplex).
      // For each fault face, sort the face vertices by global vertex ID.
      // canonical_to_local_perm_[face][k] = MFEM local DOF index for
      // the k-th vertex in sorted-global-ID order.
      // For p=1 triangles: DOF k = vertex k, so vertex perm = DOF perm.
      // ------------------------------------------------------------------
      canonical_to_local_perm_.SetSize(num_fault_faces_ * nbf_per_face_);

      // Get global vertex IDs (parallel) or use local indices (serial).
      // Use int64_t for sorting to avoid HYPRE dependency in serial builds.
      std::vector<int64_t> global_vert_ids;
      if constexpr (IsParallelMesh<MeshType>::value)
      {
#ifdef MFEM_USE_MPI
         Array<HYPRE_BigInt> gvi;
         mesh_.GetGlobalVertexIndices(gvi);
         global_vert_ids.resize(gvi.Size());
         for (int i = 0; i < gvi.Size(); i++)
         {
            global_vert_ids[i] = static_cast<int64_t>(gvi[i]);
         }
#endif
      }
      else
      {
         global_vert_ids.resize(mesh_.GetNV());
         for (int i = 0; i < mesh_.GetNV(); i++)
         {
            global_vert_ids[i] = static_cast<int64_t>(i);
         }
      }

      for (int fi = 0; fi < num_fault_faces_; fi++)
      {
         // Get local face index in the mesh
         int local_face_idx;
         if (fi < num_interior)
         {
            local_face_idx = fault_interior_faces_[fi];
         }
         else
         {
            const int si = fi - num_interior;
            if constexpr (IsParallelMesh<MeshType>::value)
            {
               local_face_idx = mesh_.GetSharedFace(fault_shared_faces_[si]);
            }
            else
            {
               MFEM_ABORT("Shared fault face in serial mesh");
               local_face_idx = -1;
            }
         }

         // Get face vertices (local indices)
         Array<int> vert;
         mesh_.GetFaceVertices(local_face_idx, vert);
         const int nv = vert.Size();

         // Build (global_id, mfem_local_index) pairs and sort by global ID
         // (Tandem sorted-simplex convention)
         std::vector<std::pair<int64_t, int>> gid_idx(nv);
         for (int k = 0; k < nv; k++)
         {
            gid_idx[k] = {global_vert_ids[vert[k]], k};
         }
         std::sort(gid_idx.begin(), gid_idx.end());

         // canonical_to_local_perm_[fi * nbf + canonical_k] = mfem_local_k
         if (nbf_per_face_ == 1)
         {
            // BR2: single centroid DOF per face, no vertex permutation
            canonical_to_local_perm_[fi] = 0;
         }
         else
         {
            // IP (nbf == nv for p=1): reorder by sorted global vertex IDs
            for (int k = 0; k < nv && k < nbf_per_face_; k++)
            {
               canonical_to_local_perm_[fi * nbf_per_face_ + k] =
                  gid_idx[k].second;
            }
            // For higher-order DOFs beyond vertices (p>=2): identity
            for (int k = nv; k < nbf_per_face_; k++)
            {
               canonical_to_local_perm_[fi * nbf_per_face_ + k] = k;
            }
         }
      }

      // Build owned DOF → local DOF map (with permutation applied)
      owned_fault_dof_to_local_dof_.SetSize(num_owned_fault_dofs_);
      for (int owned_face = 0; owned_face < num_owned_fault_faces_; owned_face++)
      {
         const int local_face = owned_fault_face_to_local_face_[owned_face];
         for (int kk = 0; kk < nbf_per_face_; kk++)
         {
            const int mfem_kk =
               canonical_to_local_perm_[local_face * nbf_per_face_ + kk];
            owned_fault_dof_to_local_dof_[owned_face * nbf_per_face_ + kk] =
               local_face * nbf_per_face_ + mfem_kk;
         }
      }
   }

   /// Facet BC table lookup for interior faces.
   FacetBC GetFaceBC(int face_idx) const
   {
      return (face_idx >= 0 && face_idx < static_cast<int>(face_bc_.size()))
             ? face_bc_[face_idx] : FacetBC::None;
   }

   /// Facet BC table lookup for shared faces.
   FacetBC GetSharedFaceBC(int shared_face) const
   {
      return (shared_face >= 0 &&
              shared_face < static_cast<int>(shared_face_bc_.size()))
             ? shared_face_bc_[shared_face] : FacetBC::None;
   }

   void SetupSolver()
   {
      X_.SetSize(fes_->GetTrueVSize());
      B_.SetSize(fes_->GetTrueVSize());

      if constexpr (IsParallelMesh<MeshType>::value)
      {
         // Parallel solver created in AssembleStiffness
      }
      else
      {
         auto *cg = new CGSolver();
         cg->SetRelTol(1e-12);
         cg->SetAbsTol(0.0);
         cg->SetMaxIter(10000);
         cg->SetPrintLevel(0);
         solver_.reset(cg);
      }
   }

   // ========================================================================
   // Mass inverse precomputation (for BR2)
   // ========================================================================

   void PrecomputeMassInverse() const
   {
      if (mass_inv_computed_) { return; }

      int ne = mesh_.GetNE();
      elem_mass_inv_.resize(ne);

      for (int e = 0; e < ne; e++)
      {
         // Use scalar FE space for mass matrix (same for all components)
         const FiniteElement *fe = scalar_fes_->GetFE(e);
         ElementTransformation *T = mesh_.GetElementTransformation(e);

         int ndof = fe->GetDof();
         DenseMatrix M(ndof);

         const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(),
                                                   2 * fe->GetOrder());
         Vector shape(ndof);
         M = 0.0;

         for (int j = 0; j < ir.GetNPoints(); j++)
         {
            const IntegrationPoint &ip = ir.IntPoint(j);
            T->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);

            real_t w = ip.weight * T->Weight();
            for (int k = 0; k < ndof; k++)
            {
               for (int l = 0; l < ndof; l++)
               {
                  M(k, l) += w * shape(k) * shape(l);
               }
            }
         }

         elem_mass_inv_[e].SetSize(ndof);
         DenseMatrixInverse M_inv(M);
         M_inv.GetInverseMatrix(elem_mass_inv_[e]);
      }

#ifdef MFEM_USE_MPI
      // In parallel, also compute mass inverses for face-neighbor elements.
      // These are needed by DGElasticityBR2Integrator::AssembleFaceMatrix when
      // ParBilinearForm::AssembleSharedFaces passes a shared face with
      // Trans.Elem2No >= mesh.GetNE().
      if constexpr (IsParallelMesh<MeshType>::value)
      {
         auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
         MFEM_VERIFY(pfes, "Expected ParFiniteElementSpace in parallel");
         pfes->ExchangeFaceNbrData();

         ParMesh *pmesh = pfes->GetParMesh();
         int nel_nbr = pmesh->GetNFaceNeighborElements();
         elem_mass_inv_.resize(ne + nel_nbr);

         for (int i = 0; i < nel_nbr; i++)
         {
            const FiniteElement *fe = pfes->GetFaceNbrFE(i);
            ElementTransformation *T =
               pmesh->GetFaceNbrElementTransformation(i);

            int ndof = fe->GetDof();
            DenseMatrix M(ndof);

            const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(),
                                                      2 * fe->GetOrder());
            Vector shape(ndof);
            M = 0.0;

            for (int j = 0; j < ir.GetNPoints(); j++)
            {
               const IntegrationPoint &ip = ir.IntPoint(j);
               T->SetIntPoint(&ip);
               fe->CalcShape(ip, shape);

               real_t w = ip.weight * T->Weight();
               for (int k = 0; k < ndof; k++)
               {
                  for (int l = 0; l < ndof; l++)
                  {
                     M(k, l) += w * shape(k) * shape(l);
                  }
               }
            }

            elem_mass_inv_[ne + i].SetSize(ndof);
            DenseMatrixInverse M_inv(M);
            M_inv.GetInverseMatrix(elem_mass_inv_[ne + i]);
         }
      }
#endif

      mass_inv_computed_ = true;
   }
