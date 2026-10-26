// elasticity_operator_fault_ops.inl
// Profiling plan item 3 (Chunhui Zhao, 2026-10-17): precomputed per-fault-face
// operators for the IP traction recovery and the slip load.  This file is
// included from elasticity_operator.hpp -- do not include directly.
//
// ComputeTractionImpl and AssembleSlipContributionIP(Shared) rebuild, every
// call and for every fault face, the face transformations, the slip at the
// quadrature points, the traction at the quadrature points and the
// nl_q-weighted L2 projection onto the face's fault DOFs.  All of it is
// linear in (u on Elem1, u on Elem2, slip of the face) and depends only on
// the fixed mesh and material, so BuildFaultFaceOps() evaluates the same
// per-face kernels once per input column (unit vectors) and stores
//
//   A_f  (3 nbf x (n1 + n2 + 2 nbf)):  [dip; strike; normal] traction at the
//        face's fault DOFs, in ComputeTraction's layout (normal positive in
//        compression), from [u_1, u_2, slip (dip, strike interleaved)];
//   G_f  ((n1 + n2) x 2 nbf), shared faces n1 rows:  the slip-load vectors
//        elvec1 (and elvec2) of AssembleSlipFaceRHS.
//
// Per step, ComputeTractionFast is then one small dense product per face and
// AddSlipLoadFast one per face, with no allocation and no whole-vector copy.
// Shared fault faces need u on the neighbour's element: instead of the full
// face-neighbour exchange (every shared face of the partition), only the
// elements adjacent to shared FAULT faces travel, between the ranks that
// own them (FaultNbrExchange_).  The kernels are the production ones, so the
// result equals ComputeTraction / AssembleSlipOnlyRHS to round-off (checked
// by tests/unit/bench_nitsche_kernels.cpp).  One difference: the old slip
// load skips a face whose slip is below 1e-15 m at every DOF; the matrix
// form adds that face's (round-off sized) contribution.

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::FaultFaceTractionKernel_(
   int idx, FaceElementTransformations &FTr, const FiniteElement &fe1,
   const FiniteElement &fe2, const Vector &u1, const Vector &u2,
   const Vector &slip_bc, Vector &out) const
{
   // Same steps as the IP branch of ComputeTractionImpl with a normal
   // traction requested and no decomposition.
   const int dim = 3, nbf = nbf_per_face_;
   const auto &basis = fault_basis_.GetBasis(idx);
   FTr.SetAllIntPoints(&Geometries.GetCenter(FTr.GetGeometryType()));
   Vector du_q;
   BuildSlipAtQuadPoints(idx, slip_bc, du_q);
   DGElasticityIPCombinedIntegrator integ(*lambda_ptr_, *mu_ptr_, dim, epsilon_,
                                          penalty_factor_);
   Vector T_q, nl_q;
   integ.ComputeTractionAtQuadPoints(fe1, fe2, FTr, u1, u2, du_q, T_q, nullptr, &nl_q);
   const int qo = 2 * std::max(fe1.GetOrder(), fe2.GetOrder()) + 1;
   const IntegrationRule &ir = IntRules.Get(FTr.GetGeometryType(), qo);
   real_t bv[3][3];
   for (int d = 0; d < 3; d++)
   {
      bv[0][d] = basis.normal[d];
      bv[1][d] = basis.tangent1[d];
      bv[2][d] = basis.tangent2[d];
   }
   Vector tl;
   DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
      dim, 3, T_q, nl_q, ir, nbf, face_quad_->BasisAtQuadPoints(), bv, tl,
      basis.qp_data.empty() ? nullptr : &basis.qp_data);
   out.SetSize(3 * nbf);
   for (int kk = 0; kk < nbf; kk++)
   {
      out(kk) = tl(1 * nbf + kk);            // dip
      out(nbf + kk) = tl(2 * nbf + kk);      // strike
      out(2 * nbf + kk) = -tl(0 * nbf + kk); // normal, positive in compression
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::BuildFaultFaceOps()
{
   Require3D_("BuildFaultFaceOps");
   MFEM_VERIFY(method_ == DGMethod::IP, "BuildFaultFaceOps: IP method only");
   const int nfi = fault_interior_faces_.Size(), nfs = fault_shared_faces_.Size();
   const int nbf = nbf_per_face_, nf = nfi + nfs;
   MFEM_VERIFY(num_fault_dofs_ == nf * nbf, "BuildFaultFaceOps: " << num_fault_dofs_
               << " fault DOFs != " << nf << " faces x " << nbf);
#ifdef MFEM_USE_MPI
   ParFiniteElementSpace *pfes = nullptr;
#endif
#ifdef MFEM_USE_MPI
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      mesh_.ExchangeFaceNbrData();
      pfes = dynamic_cast<ParFiniteElementSpace *>(fes_.get());
      pfes->ExchangeFaceNbrData();
   }
#endif
   ffo_.nf = nf;
   ffo_.nfi = nfi;
   ffo_.off_vd.assign(nf + 1, 0);
   ffo_.off_A.assign(nf + 1, 0);
   ffo_.off_G.assign(nf + 1, 0);
   ffo_.n1.assign(nf, 0);
   ffo_.n2.assign(nf, 0);
   ffo_.vd.clear();
   ffo_.A.clear();
   ffo_.G.clear();
   ffo_.nbr_face_elem.assign(nfs, -1);

   Vector probe(2 * num_fault_dofs_);
   probe = 0.0;
   DGElasticityIPCombinedIntegrator slip_integ(*lambda_ptr_, *mu_ptr_, 3, epsilon_,
                                               penalty_factor_);
   Vector out, ev1, ev2, du_q;
   for (int f = 0; f < nf; f++)
   {
      const bool shared = (f >= nfi);
      FaceElementTransformations *FTr = nullptr;
      const FiniteElement *fe1 = nullptr, *fe2 = nullptr;
      Array<int> vd1, vd2;
      if (!shared)
      {
         FTr = mesh_.GetInteriorFaceTransformations(fault_interior_faces_[f]);
         MFEM_VERIFY(FTr, "BuildFaultFaceOps: no transformation for fault face " << f);
         fe1 = scalar_fes_->GetFE(FTr->Elem1No);
         fe2 = scalar_fes_->GetFE(FTr->Elem2No);
         fes_->GetElementVDofs(FTr->Elem1No, vd1);
         fes_->GetElementVDofs(FTr->Elem2No, vd2);
      }
      else
      {
#ifdef MFEM_USE_MPI
         if constexpr (IsParallelMesh<MeshType>::value)
         {
            FTr = mesh_.GetSharedFaceTransformations(fault_shared_faces_[f - nfi]);
            MFEM_VERIFY(FTr, "BuildFaultFaceOps: no transformation for shared fault face "
                        << f - nfi);
            const int nbr = FTr->Elem2No - mesh_.GetNE();
            fe1 = scalar_fes_->GetFE(FTr->Elem1No);
            fe2 = pfes->GetFaceNbrFE(nbr);
            fes_->GetElementVDofs(FTr->Elem1No, vd1);
            ffo_.nbr_face_elem[f - nfi] = nbr;
            vd2.SetSize(3 * fe2->GetDof());   // values come from the exchange
         }
#endif
      }
      for (int j = 0; j < vd1.Size(); j++)
      {
         MFEM_VERIFY(vd1[j] >= 0, "BuildFaultFaceOps: signed vdof (not DG?)");
      }
      for (int j = 0; j < vd2.Size() && !shared; j++)
      {
         MFEM_VERIFY(vd2[j] >= 0, "BuildFaultFaceOps: signed vdof (not DG?)");
      }
      const int n1 = vd1.Size(), n2 = vd2.Size(), nc = n1 + n2 + 2 * nbf;
      ffo_.n1[f] = n1;
      ffo_.n2[f] = n2;
      for (int j = 0; j < n1; j++) { ffo_.vd.push_back(vd1[j]); }
      for (int j = 0; j < n2; j++) { ffo_.vd.push_back(shared ? -1 : vd2[j]); }
      ffo_.off_vd[f + 1] = static_cast<int>(ffo_.vd.size());

      // A_f, column by column (row-major storage, 3 nbf rows)
      const std::size_t a0 = ffo_.A.size();
      ffo_.A.resize(a0 + static_cast<std::size_t>(3 * nbf) * nc, 0.0);
      Vector u1(n1), u2(n2);
      for (int c = 0; c < nc; c++)
      {
         u1 = 0.0;
         u2 = 0.0;
         int sj = -1;
         if (c < n1) { u1(c) = 1.0; }
         else if (c < n1 + n2) { u2(c - n1) = 1.0; }
         else { sj = 2 * f * nbf + (c - n1 - n2); probe(sj) = 1.0; }
         FaultFaceTractionKernel_(f, *FTr, *fe1, *fe2, u1, u2, probe, out);
         if (sj >= 0) { probe(sj) = 0.0; }
         for (int r = 0; r < 3 * nbf; r++) { ffo_.A[a0 + static_cast<std::size_t>(r) * nc + c] = out(r); }
      }
      ffo_.off_A[f + 1] = static_cast<int>(ffo_.A.size());

      // G_f: slip-load vectors, column by column (row-major, ng rows)
      const int ng = shared ? n1 : n1 + n2;
      const std::size_t g0 = ffo_.G.size();
      ffo_.G.resize(g0 + static_cast<std::size_t>(ng) * 2 * nbf, 0.0);
      for (int c = 0; c < 2 * nbf; c++)
      {
         const int sj = 2 * f * nbf + c;
         probe(sj) = 1.0;
         BuildSlipAtQuadPoints(f, probe, du_q);
         probe(sj) = 0.0;
         slip_integ.AssembleSlipFaceRHS(*fe1, *fe2, *FTr, du_q, ev1, ev2);
         for (int r = 0; r < n1; r++) { ffo_.G[g0 + static_cast<std::size_t>(r) * 2 * nbf + c] = ev1(r); }
         if (!shared)
         {
            for (int r = 0; r < n2; r++)
            {
               ffo_.G[g0 + static_cast<std::size_t>(n1 + r) * 2 * nbf + c] = ev2(r);
            }
         }
      }
      ffo_.off_G[f + 1] = static_cast<int>(ffo_.G.size());
   }
   BuildFaultNbrExchange_();
   ffo_.built = true;
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::BuildFaultNbrExchange_()
{
   ffo_.recv_ranks.clear();
   ffo_.recv_off.assign(1, 0);
   ffo_.send_ranks.clear();
   ffo_.send_off.assign(1, 0);
   ffo_.send_elems.clear();
   ffo_.nbr_off.assign(fault_shared_faces_.Size(), -1);
#ifdef MFEM_USE_MPI
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      MPI_Comm comm = mesh_.GetComm();
      int nranks = 1;
      MPI_Comm_size(comm, &nranks);
      const int nfs = fault_shared_faces_.Size();
      const int vdim = 3 * scalar_fes_->GetFE(0)->GetDof();
      // positions (in the owner's send_face_nbr_elements row) requested
      // from each neighbour rank, in the order of the shared fault faces
      std::map<int, std::vector<int>> want;       // rank -> positions
      std::map<int, std::vector<int>> want_face;  // rank -> shared face index
      const Array<int> &off = mesh_.face_nbr_elements_offset;
      for (int i = 0; i < nfs; i++)
      {
         const int nbr = ffo_.nbr_face_elem[i];
         int fn = 0;
         while (fn + 1 < off.Size() && off[fn + 1] <= nbr) { fn++; }
         const int r = mesh_.GetFaceNbrRank(fn);
         want[r].push_back(nbr - off[fn]);
         want_face[r].push_back(i);
      }
      std::vector<int> nsend_req(nranks, 0), nrecv_req(nranks, 0);
      for (auto &kv : want) { nsend_req[kv.first] = static_cast<int>(kv.second.size()); }
      MPI_Alltoall(nsend_req.data(), 1, MPI_INT, nrecv_req.data(), 1, MPI_INT, comm);
      const int tag = 7311;
      std::vector<MPI_Request> req;
      std::map<int, std::vector<int>> asked;   // rank -> positions it asks of me
      for (int r = 0; r < nranks; r++)
      {
         if (nrecv_req[r] > 0)
         {
            asked[r].resize(nrecv_req[r]);
            req.emplace_back();
            MPI_Irecv(asked[r].data(), nrecv_req[r], MPI_INT, r, tag, comm, &req.back());
         }
      }
      for (auto &kv : want)
      {
         req.emplace_back();
         MPI_Isend(kv.second.data(), static_cast<int>(kv.second.size()), MPI_INT, kv.first,
                   tag, comm, &req.back());
      }
      MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
      // what I receive: per rank, its elements in my request order
      int roff = 0;
      for (auto &kv : want)
      {
         ffo_.recv_ranks.push_back(kv.first);
         const auto &faces = want_face[kv.first];
         for (std::size_t k = 0; k < faces.size(); k++)
         {
            ffo_.nbr_off[faces[k]] = roff + static_cast<int>(k) * vdim;
         }
         roff += static_cast<int>(faces.size()) * vdim;
         ffo_.recv_off.push_back(roff);
      }
      ffo_.recv_buf.assign(roff, 0.0);
      // what I send: per asking rank, my local elements for its positions
      const Table &sfe = mesh_.send_face_nbr_elements;
      for (auto &kv : asked)
      {
         int fn = -1;
         for (int g = 0; g < mesh_.GetNFaceNeighbors(); g++)
         {
            if (mesh_.GetFaceNbrRank(g) == kv.first) { fn = g; }
         }
         MFEM_VERIFY(fn >= 0, "BuildFaultFaceOps: rank " << kv.first
                     << " asks for fault-face neighbour elements but is not a face neighbour");
         const int *row = sfe.GetRow(fn);
         ffo_.send_ranks.push_back(kv.first);
         for (int p : kv.second)
         {
            MFEM_VERIFY(p >= 0 && p < sfe.RowSize(fn), "BuildFaultFaceOps: bad position");
            ffo_.send_elems.push_back(row[p]);
         }
         ffo_.send_off.push_back(static_cast<int>(ffo_.send_elems.size()));
      }
      ffo_.send_buf.assign(ffo_.send_elems.size() * vdim, 0.0);
      ffo_.send_vd.clear();
      Array<int> vd;
      for (int e : ffo_.send_elems)
      {
         fes_->GetElementVDofs(e, vd);
         for (int j = 0; j < vd.Size(); j++) { ffo_.send_vd.push_back(vd[j]); }
      }
      for (int i = 0; i < nfs; i++)
      {
         MFEM_VERIFY(ffo_.nbr_off[i] >= 0, "BuildFaultFaceOps: shared fault face " << i
                     << " has no neighbour data slot");
      }
   }
#endif
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::FaultNbrExchange_(const Vector &u) const
{
#ifdef MFEM_USE_MPI
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      if (ffo_.send_ranks.empty() && ffo_.recv_ranks.empty()) { return; }
      MPI_Comm comm = mesh_.GetComm();
      const int tag = 7312;
      const int vdim = ffo_.send_elems.empty() ? 0
                       : static_cast<int>(ffo_.send_vd.size() / ffo_.send_elems.size());
      for (std::size_t i = 0; i < ffo_.send_vd.size(); i++) { ffo_.send_buf[i] = u(ffo_.send_vd[i]); }
      ffo_.reqs.resize(ffo_.recv_ranks.size() + ffo_.send_ranks.size());
      int q = 0;
      for (std::size_t k = 0; k < ffo_.recv_ranks.size(); k++)
      {
         MPI_Irecv(ffo_.recv_buf.data() + ffo_.recv_off[k], ffo_.recv_off[k + 1] - ffo_.recv_off[k],
                   MPITypeMap<real_t>::mpi_type, ffo_.recv_ranks[k], tag, comm, &ffo_.reqs[q++]);
      }
      for (std::size_t k = 0; k < ffo_.send_ranks.size(); k++)
      {
         const int b = ffo_.send_off[k] * vdim, n = (ffo_.send_off[k + 1] - ffo_.send_off[k]) * vdim;
         MPI_Isend(ffo_.send_buf.data() + b, n, MPITypeMap<real_t>::mpi_type, ffo_.send_ranks[k],
                   tag, comm, &ffo_.reqs[q++]);
      }
      MPI_Waitall(q, ffo_.reqs.data(), MPI_STATUSES_IGNORE);
   }
#else
   (void)u;
#endif
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::ComputeTractionFast(
   const Vector &u, const Vector &slip_bc, Vector &traction, Vector &normal_traction) const
{
   MFEM_VERIFY(ffo_.built, "ComputeTractionFast: call BuildFaultFaceOps first");
   const int nbf = nbf_per_face_;
   traction.SetSize(2 * num_fault_dofs_);
   normal_traction.SetSize(num_fault_dofs_);
   FaultNbrExchange_(u);
   const real_t *pu = u.GetData(), *ps = slip_bc.GetData();
   real_t *pt = traction.GetData(), *pn = normal_traction.GetData();
   std::vector<real_t> &x = ffo_.x;
   for (int f = 0; f < ffo_.nf; f++)
   {
      const int n1 = ffo_.n1[f], n2 = ffo_.n2[f], nc = n1 + n2 + 2 * nbf;
      x.resize(nc);
      const int *vd = &ffo_.vd[ffo_.off_vd[f]];
      for (int j = 0; j < n1; j++) { x[j] = pu[vd[j]]; }
      if (f < ffo_.nfi)
      {
         for (int j = 0; j < n2; j++) { x[n1 + j] = pu[vd[n1 + j]]; }
      }
      else
      {
         const real_t *nb = ffo_.recv_buf.data() + ffo_.nbr_off[f - ffo_.nfi];
         for (int j = 0; j < n2; j++) { x[n1 + j] = nb[j]; }
      }
      for (int j = 0; j < 2 * nbf; j++) { x[n1 + n2 + j] = ps[2 * f * nbf + j]; }
      const real_t *A = &ffo_.A[ffo_.off_A[f]];
      for (int kk = 0; kk < nbf; kk++)
      {
         real_t s[3] = {0.0, 0.0, 0.0};
         for (int r = 0; r < 3; r++)
         {
            const real_t *a = A + static_cast<std::size_t>(r * nbf + kk) * nc;
            real_t t = 0.0;
            for (int c = 0; c < nc; c++) { t += a[c] * x[c]; }
            s[r] = t;
         }
         const int k = f * nbf + kk;
         pt[2 * k] = s[0];
         pt[2 * k + 1] = s[1];
         pn[k] = s[2];
      }
   }
}

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::AddSlipLoadFast(const Vector &slip_bc,
                                                         Vector &rhs) const
{
   MFEM_VERIFY(ffo_.built, "AddSlipLoadFast: call BuildFaultFaceOps first");
   const int nbf = nbf_per_face_, m = 2 * nbf;
   const real_t *ps = slip_bc.GetData();
   real_t *pr = rhs.GetData();
   for (int f = 0; f < ffo_.nf; f++)
   {
      const real_t *s = ps + 2 * f * nbf;
      const int ng = (f < ffo_.nfi) ? ffo_.n1[f] + ffo_.n2[f] : ffo_.n1[f];
      const int *vd = &ffo_.vd[ffo_.off_vd[f]];
      const real_t *G = &ffo_.G[ffo_.off_G[f]];
      for (int r = 0; r < ng; r++)
      {
         const real_t *g = G + static_cast<std::size_t>(r) * m;
         real_t t = 0.0;
         for (int c = 0; c < m; c++) { t += g[c] * s[c]; }
         pr[vd[r]] += t;
      }
   }
}
