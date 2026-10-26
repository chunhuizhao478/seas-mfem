// elasticity_operator_penalty.inl
// C-04 (Nitsche dynamic rupture plan rev 8.2): per-fault-face interior
// penalty sigma_F.  This file is included from elasticity_operator.hpp --
// do not include directly.
//
// Walks fault_interior_faces_ then fault_shared_faces_ exactly as
// ComputeTractionImpl does, and at every point of the rule
// DGElasticityIPCombinedIntegrator::ComputeTractionAtQuadPoints uses
// (order 2 max(p1, p2) + 1) calls the integrator's public GetPenalty with the
// same per-side inputs that function passes to ComputePenalty: the two
// element Jacobian determinants, the per-side (lambda, mu) evaluated through
// lambda_ptr_/mu_ptr_ (C-09) and the length of the unnormalised face normal.
//
// Scope (Chunhui Zhao, 2026-10-14).

template <typename MeshType>
void ElasticityDomainOperator<MeshType>::GetFaultFacePenalties(
   Vector &sigma_F_face, Vector &spread) const
{
   Require3D_("GetFaultFacePenalties");
   MFEM_VERIFY(method_ == DGMethod::IP,
               "GetFaultFacePenalties: IP method only (sigma_F is the "
               "interior-penalty coefficient)");
   const int dim = 3;
   const int nfaces = fault_interior_faces_.Size() + fault_shared_faces_.Size();
   sigma_F_face.SetSize(nfaces);
   spread.SetSize(nfaces);
   sigma_F_face = 0.0;
   spread = 0.0;

   DGElasticityIPCombinedIntegrator integ(*lambda_ptr_, *mu_ptr_, dim,
                                          epsilon_, penalty_factor_);

   // sigma_F at every traction quadrature point of one face -> (mean, spread)
   auto face_values = [&](FaceElementTransformations &Tr,
                          const FiniteElement &fe1, const FiniteElement &fe2,
                          int f)
   {
      const int qo = 2 * std::max(fe1.GetOrder(), fe2.GetOrder()) + 1;
      const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), qo);
      const int nq = ir.GetNPoints();
      std::vector<real_t> s(nq);
      Vector nor(dim);
      real_t sum = 0.0;
      for (int q = 0; q < nq; q++)
      {
         Tr.SetAllIntPoints(&ir.IntPoint(q));
         const IntegrationPoint &eip1 = Tr.GetElement1IntPoint();
         const IntegrationPoint &eip2 = Tr.GetElement2IntPoint();
         CalcOrtho(Tr.Jacobian(), nor);
         const real_t nl = nor.Norml2();
         const real_t lam1 = lambda_ptr_->Eval(*Tr.Elem1, eip1);
         const real_t mu1 = mu_ptr_->Eval(*Tr.Elem1, eip1);
         const real_t lam2 = lambda_ptr_->Eval(*Tr.Elem2, eip2);
         const real_t mu2 = mu_ptr_->Eval(*Tr.Elem2, eip2);
         s[q] = integ.GetPenalty(fe1, fe2, Tr.Elem1->Weight(), Tr.Elem2->Weight(),
                                 lam1, mu1, lam2, mu2, nl);
         sum += s[q];
      }
      const real_t mean = sum / nq;
      real_t sp = 0.0;
      for (int q = 0; q < nq; q++) { sp = std::max(sp, std::abs(s[q] - mean)); }
      sigma_F_face(f) = mean;
      spread(f) = sp / mean;
   };

   for (int fi = 0; fi < fault_interior_faces_.Size(); fi++)
   {
      FaceElementTransformations *FTr =
         mesh_.GetInteriorFaceTransformations(fault_interior_faces_[fi]);
      MFEM_VERIFY(FTr != nullptr, "GetFaultFacePenalties: no transformation "
                  "for fault face " << fault_interior_faces_[fi]);
      face_values(*FTr, *scalar_fes_->GetFE(FTr->Elem1No),
                  *scalar_fes_->GetFE(FTr->Elem2No), fi);
   }

   if constexpr (IsParallelMesh<MeshType>::value)
   {
#ifdef MFEM_USE_MPI
      auto *pfes = dynamic_cast<ParFiniteElementSpace*>(fes_.get());
      pfes->ExchangeFaceNbrData();   // collective: every rank calls it
      for (int i = 0; i < fault_shared_faces_.Size(); i++)
      {
         FaceElementTransformations *FTr =
            mesh_.GetSharedFaceTransformations(fault_shared_faces_[i]);
         MFEM_VERIFY(FTr != nullptr, "GetFaultFacePenalties: no transformation "
                     "for shared fault face " << fault_shared_faces_[i]);
         const int nbr_idx = FTr->Elem2No - mesh_.GetNE();
         face_values(*FTr, *scalar_fes_->GetFE(FTr->Elem1No),
                     *pfes->GetFaceNbrFE(nbr_idx),
                     fault_interior_faces_.Size() + i);
      }
#endif
   }
}
