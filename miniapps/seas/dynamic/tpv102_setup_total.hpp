// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.3.0 Phase 4 (I-06 part B): total-stress companion to
// tpv102_setup.hpp.  Under the Phase 4 migration, TPV102's wave state Q
// carries the background pre-stress tensor in every bulk DOF (in global
// Cartesian coordinates), so the fault-face flux dispatches to
// FaultFaceFlux::EvaluateTotal on already-total inputs and DOFData
// pre-stress fields (sigma_n0, tau1_0, tau2_0) are zeroed at init to
// avoid double-counting.
//
// This file MUST NOT modify tpv102_setup.hpp — BP5 uses tpv102_setup.hpp's
// DOFData-initialization utilities, and the CLAUDE.md rule "TPV102 fix:
// reuse BP5 patterns, never modify BP5 source" applies.

#ifndef MFEM_SEAS_TPV102_SETUP_TOTAL_HPP
#define MFEM_SEAS_TPV102_SETUP_TOTAL_HPP

#include "mfem.hpp"
#include "wave_state.hpp"
#include "fault_face_flux.hpp"
#include "../common/seas_types.hpp"
#include "../config/tpv102_params.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Initialize the bulk Q vector with the TPV102 pre-stress tensor
///        in global coordinates (I-06 Phase 4).
///
/// For TPV102 (homogeneous half-space under a vertical strike-slip fault
/// at y=0, BP5 canonical frame n=(0,-1,0), t1=(0,0,-1), t2=(+1,0,0)),
/// the rotation of the canonical-local tensor
///     sigma_nn_local = +sigma_n0   (compression)
///     sigma_nt2_local = +tau_ini   (along-strike shear)
///     (all other components zero)
/// into GLOBAL Cartesian coordinates gives:
///     sigma_yy_global = +sigma_n0  (compression on the y=const plane,
///                                   matches MFEM compression-positive)
///     sigma_xy_global = -tau_ini   (sign flip from canonical-local
///                                   sigma_nt2 because t2 = +x, n = -y
///                                   under Voigt rotation)
/// All other tensor components vanish in exact arithmetic.  See the
/// rotation check in the review-incorporation plan §R-001.
///
/// Under this layout, ComputeTrialTraction on TOTAL Q produces the
/// TOTAL trial traction directly (eta_p * 2 * sigma_nn_local / Zp =
/// sigma_n0 under homogeneous eta_p = Zp/2), so EvaluateTotal can
/// skip the Evaluate path's sigma_n0 reconstruction.
///
/// @param[out] Q            State vector to resize + fill.
/// @param[in]  ndof_total   Total scalar DOFs per component.
/// @param[in]  sigma_n0     TPV102 normal pre-stress [Pa] (positive = compression).
/// @param[in]  tau_ini      TPV102 initial along-strike shear [Pa].
inline void InitializeStateTotal(Vector &Q, int ndof_total,
                                 real_t sigma_n0, real_t tau_ini)
{
   MFEM_VERIFY(ndof_total > 0,
               "InitializeStateTotal: ndof_total must be positive, got "
               << ndof_total);

   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;

   for (int i = 0; i < ndof_total; i++)
   {
      Q[SYY * ndof_total + i] =  sigma_n0;   // R-001: +sigma_n0 (compression)
      Q[SXY * ndof_total + i] = -tau_ini;    // R-001: -tau_ini  (canonical-
                                             //                   frame sign flip)
   }
}

/// @brief Zero the pre-stress fields of every DOFData entry (I-06 Phase 4).
///
/// Under the v9.3.0 total-stress migration, the driver bakes the physical
/// pre-stress into bulk Q via InitializeStateTotal.  DOFData's
/// sigma_n0 / tau1_0 / tau2_0 MUST be zeroed so EvaluateTotal's
/// ComputeTrialTraction (which reads Q) does not double-count pre-stress
/// via DOFData.  All OTHER DOFData fields are left alone.  In particular
/// InitializeFaultDOFs pre-seeds several fields with physically
/// meaningful t=0 values so the first ParaView / station snapshot
/// matches user expectations; each of these is OUTPUT-only for
/// EvaluateTotal (or INPUT but must survive the zeroing):
///
///   - sigma_n_corr / tau_{1,2}_corr : seeded with the physical pre-
///     stress.  Overwritten on the first EvaluateTotal call.  OUTPUT.
///   - slip_rate / V1 / V2           : seeded (V_ini, 0, V_ini) —
///     TPV102 is pure strike-slip under BP5 canonical convention.
///     Also overwritten on the first EvaluateTotal call.  OUTPUT.
///   - slip1 / slip2                 : accumulated slip, seeded 0.
///     Updated by the driver's RK4 Simpson-integral block.  OUTPUT.
///   - psi                           : seeded at steady-state for t=0
///     equilibrium; treated as INPUT state variable by EvaluateTotal
///     (read but not written; integrated by the driver's coupled-RK4-
///     on-psi block).  INPUT + driver-integrated.
///   - Zp_{plus,minus}, Zs_{plus,minus}, eta_{p,s}, a, Dc : material
///     parameters.  INPUT.
///
/// A refactor that "notices" the V_ini / psi seeds and adds them to
/// this function would zero the t=0 snapshot (for V_ini) or cripple
/// the equilibrium solve (for psi).  Do NOT extend beyond the three
/// pre-stress fields without revisiting the station-output contract.
///
/// @param[in,out] dof_data  Array of DOFData already filled by
///                          InitializeFaultDOFs.
/// @param[in]     ndof      Number of fault DOFs to clear.
inline void ZeroDOFDataPreStressTotal(std::vector<DOFData> &dof_data, int ndof)
{
   MFEM_VERIFY(static_cast<int>(dof_data.size()) >= ndof,
               "ZeroDOFDataPreStressTotal: dof_data size "
               << dof_data.size() << " < ndof " << ndof);

   for (int i = 0; i < ndof; i++)
   {
      dof_data[i].sigma_n0   = 0.0;
      dof_data[i].tau1_0     = 0.0;
      dof_data[i].tau2_0     = 0.0;
      // Persistent-prestress channel: under total-Q the nucleation
      // amplitude lives in tau*_nuc and is overwritten per step by
      // ApplyNucleationTotalPrestress.  Reset here so a re-init that
      // calls InitializeFaultDOFs + ZeroDOFDataPreStressTotal in
      // sequence (e.g. checkpoint restart, test fixture rebuild) does
      // NOT carry stale nucleation state from a prior run.
      dof_data[i].sigma_n_nuc = 0.0;
      dof_data[i].tau1_nuc    = 0.0;
      dof_data[i].tau2_nuc    = 0.0;
   }
}

/// @brief Per-fault-QP nodal-injection map (I-06 Phase 5).
///
/// Under the v9.3.0 total-stress migration, nucleation perturbs Q[SXY] at
/// the fault QPs (see the R-002 sign note below) instead of writing to
/// DOFData.tau2_0.  For a GaussLobatto interpolatory basis, injecting at
/// a fault QP reduces to updating the single element DOF whose shape
/// function peaks at that QP (the Lagrange node closest to the QP in
/// each element's reference coordinates).
///
/// Per R-009 (review-incorporation plan): both the + and - sides of the
/// fault must be updated at each QP — the trial traction in
/// EvaluateTotal averages Q+[SXY]/Zs + Q-[SXY]/Zs across the two sides,
/// so a one-sided injection would halve the perturbation.
///
/// Partition-seam fault faces (shared between MPI ranks) only have one
/// locally-owned side per rank; the other side's DOF update is performed
/// on the peer rank.  Entries set to -1 in this map signal "non-local,
/// skip".
struct FaultQPNodalMap
{
   /// Per fault QP (indexed in the same order as the driver's
   /// fault_coords: interior QPs first, then shared QPs):
   ///   qp_to_elem_{plus,minus}[i]  = element index, or -1 if non-local
   ///   qp_to_dof_{plus,minus}[i]   = local DOF index within that element,
   ///                                 or -1 if non-local
   ///   qp_shape_max_{plus,minus}[i] = shape value of the picked nodal
   ///                                 DOF at the QP (addresses R-I06-003:
   ///                                 used to scale the ApplyNucleation
   ///                                 delta so the shape-weighted QP sum
   ///                                 still equals the requested dtau
   ///                                 when the face QP does NOT coincide
   ///                                 with a GaussLobatto Lagrange node).
   ///                                 0.0 for non-local sides.
   std::vector<int>    qp_to_elem_plus;
   std::vector<int>    qp_to_dof_plus;
   std::vector<int>    qp_to_elem_minus;
   std::vector<int>    qp_to_dof_minus;
   std::vector<real_t> qp_shape_max_plus;
   std::vector<real_t> qp_shape_max_minus;

   /// Sanity: after Build(), every entry must be either (valid, valid)
   /// on the local rank or (-1, -1) for the non-local side.
   bool Consistent() const
   {
      const int n = static_cast<int>(qp_to_elem_plus.size());
      if (static_cast<int>(qp_to_dof_plus.size())       != n ||
          static_cast<int>(qp_to_elem_minus.size())     != n ||
          static_cast<int>(qp_to_dof_minus.size())      != n ||
          static_cast<int>(qp_shape_max_plus.size())    != n ||
          static_cast<int>(qp_shape_max_minus.size())   != n)
      {
         return false;
      }
      for (int i = 0; i < n; i++)
      {
         if ((qp_to_elem_plus[i]  < 0) != (qp_to_dof_plus[i]  < 0)) { return false; }
         if ((qp_to_elem_minus[i] < 0) != (qp_to_dof_minus[i] < 0)) { return false; }
      }
      return true;
   }
};

/// @brief Pick the interpolatory-nodal DOF of element @p elem_id at a
/// face integration point (I-06 Phase 5, R-009).
///
/// For a GaussLobatto basis, the element Lagrange node closest to the QP
/// has shape-function value ~1 and all others ~0 — so injecting a field
/// at the QP reduces to a single-DOF update at that node.  This helper
/// evaluates all element shape functions at @p ip_ref and returns the
/// DOF index with the largest value, together with the value itself for
/// a caller-side interpolatory check.
///
/// @param[in]  fe        Element finite element (L2 + GaussLobatto).
/// @param[in]  ip_ref    Reference-coordinate integration point in @p fe.
/// @param[out] shape_max On return, the max shape value at ip_ref.
/// @return Local DOF index with shape_max (or -1 on empty element).
inline int PickNearestNodalDOF(const FiniteElement &fe,
                               const IntegrationPoint &ip_ref,
                               real_t &shape_max)
{
   const int ndof = fe.GetDof();
   if (ndof <= 0) { shape_max = 0.0; return -1; }
   Vector shape(ndof);
   fe.CalcShape(ip_ref, shape);
   int best = 0;
   real_t bestv = shape(0);
   for (int d = 1; d < ndof; d++)
   {
      if (shape(d) > bestv) { bestv = shape(d); best = d; }
   }
   shape_max = bestv;
   return best;
}

/// @brief Verify that no (elem, dof) pair repeats across fault QPs on
/// the same side (I-06 round 3, R-002).
///
/// If two fault QPs pick the same nearest-nodal DOF on the same side,
/// `ApplyNucleationTotal`'s cumulative add-delta at that DOF double-
/// counts the nucleation amplitude on every colliding QP (each call
/// adds `-delta/shape_max` independently).  Abort at setup rather
/// than silently double-count during the RK4 loop.
///
/// For standard MFEM L2 + GaussLobatto at order ≥ 2 on a tet with
/// Dunavant face QPs, each face has nqp_per_face ≈ nqp_per_face_DOFs
/// and the nearest-nodal Voronoi cells of distinct Lagrange nodes
/// generally contain distinct QPs (no collision).  Lower-order meshes
/// or tensor-product face quadratures may collide.
inline void VerifyFaultQPNodalMapUnique(const FaultQPNodalMap &map)
{
   auto check_unique = [](const std::vector<int> &elems,
                          const std::vector<int> &dofs,
                          const char *label)
   {
      // Pack (elem, dof) into a single 64-bit key; skip entries where
      // elem < 0 (non-local sides).  elem < 2^44 and dof < 2^20 fit
      // safely; TPV102 meshes have elem counts well below 2^44.
      const int n = static_cast<int>(elems.size());
      std::vector<std::uint64_t> keys;
      keys.reserve(n);
      for (int i = 0; i < n; i++)
      {
         const int e = elems[i];
         const int d = dofs[i];
         if (e < 0 || d < 0) { continue; }
         const std::uint64_t key =
            (static_cast<std::uint64_t>(e) << 20) |
            static_cast<std::uint64_t>(d);
         keys.push_back(key);
      }
      std::sort(keys.begin(), keys.end());
      for (std::size_t k = 1; k < keys.size(); k++)
      {
         MFEM_VERIFY(keys[k] != keys[k-1],
                     "VerifyFaultQPNodalMapUnique: collision on " << label
                     << " side — two fault QPs share (elem="
                     << (keys[k] >> 20) << ", dof="
                     << (keys[k] & ((std::uint64_t(1)<<20) - 1))
                     << ").  ApplyNucleationTotal's cumulative add-"
                     "delta would double-count nucleation amplitude "
                     "at this DOF.  The face quadrature / basis "
                     "combination places two QPs in the "
                     "\"nearest-node\" Voronoi cell of one Lagrange "
                     "node; switch to a richer face quadrature or "
                     "the deferred L2-projection nucleation path.");
      }
   };
   check_unique(map.qp_to_elem_plus,  map.qp_to_dof_plus,  "+");
   check_unique(map.qp_to_elem_minus, map.qp_to_dof_minus, "-");
}

/// @brief Build the FaultQPNodalMap for TPV102 (I-06 Phase 5).
///
/// Iterates the driver's fault-face list in the SAME order as
/// fault_coords (interior QPs first, then shared QPs).  For each QP:
///   - Interior fault face: both Elem1 and Elem2 are local; populate
///     both + and - side entries based on which element lies on the
///     canonical "+" side.  The discriminator projects
///     (elem_centroid - face_centroid) onto -ref_normal: positive
///     dot-product means the element is on the + side (opposite to
///     where ref_normal points).  For TPV102 with ref_normal = (0,-1,0)
///     this reduces to "centroid_y > face_y_midpoint ⇒ + side".  Caller
///     may override ref_normal for dipping / non-planar fault
///     geometries (R-I06-008).
///   - Shared fault face: only Elem1 is local; populate only the
///     corresponding side (+ or -); the peer rank populates the other.
///
/// Verifies at each QP that the picked DOF has shape_max > kShapeMaxFloor
/// (ApplyNucleationTotal's 1/shape_max scaling handles the rest).
///
/// @param[out] map           Populated on return.
/// @param[in]  mesh          The mesh (Mesh or ParMesh).
/// @param[in]  fes           The L2 FESpace (order / basis read here).
/// @param[in]  fault_int_faces Interior fault face indices (wave operator).
/// @param[in]  fault_shr_faces Shared fault face indices (ParMesh only).
/// @param[in]  order         DG polynomial order (for integration rule).
/// @param[in]  nqp_per_face  Number of QPs per fault face (integration rule).
/// @param[in]  ref_normal    Optional canonical fault normal (3 comps).
///                           Default = (0, -1, 0) matches BP5 / TPV102.
///                           Any driver with a non-y=0 planar fault
///                           (e.g., dipping fault) must pass its own
///                           ref_normal.
template <typename MeshType>
inline void BuildFaultQPNodalMap(FaultQPNodalMap &map,
                                 MeshType &mesh,
                                 FESpaceForMesh<MeshType> &fes,
                                 const Array<int> &fault_int_faces,
                                 const Array<int> &fault_shr_faces,
                                 int order, int nqp_per_face,
                                 const real_t *ref_normal = nullptr)
{
   (void)fes;
   // Shape-max floor for the picked nearest-nodal DOF at each fault QP.
   // Addresses round 2 REVIEW.md R-I06-006 (the previous 0.7 threshold
   // was unjustified and empirically aborted the driver at setup).
   //
   // Empirical shape_max (measured on a 2-tet fixture, MFEM L2 +
   // BasisType::GaussLobatto, default IntRules.Get(TRIANGLE, 2*order)):
   //    order=1, Dunavant 3-pt :  shape_max = 0.667 at every QP
   //    order=2, Dunavant 6-pt :  shape_max in {0.518, 0.796}
   //    order=3, Dunavant 12-pt:  shape_max in {0.396, 0.775, 0.841}
   //
   // Because ApplyNucleationTotal scales the injected delta by
   // 1/shape_max (REVIEW.md R-I06-003 fix), the only requirement this
   // floor enforces is that the picked DOF has POSITIVE shape at the
   // QP — i.e., the basis is nodal-interpolatory in the sense that
   // SOME Lagrange node is the "nearest" one with non-negligible
   // shape.  A modal / hierarchical basis would have shape_max ~ 0 at
   // mid-face QPs; a GaussLobatto L2 basis always has shape_max > 0.3
   // at every Gauss / Dunavant QP observed above.
   //
   // Keeping the floor at 0.1 (kShapeMaxFloor below) still catches a
   // non-nodal basis without rejecting any of the observed orders.  Do
   // NOT raise this back to 0.7 without implementing the plan's
   // deferred L2-projection path — the scaling in ApplyNucleationTotal
   // delivers the correct QP-effective amplitude as long as
   // shape_max > 0.
   const real_t kShapeMaxFloor = 0.1;

   const int n_int = fault_int_faces.Size();
   const int n_shr = fault_shr_faces.Size();
   const int n_qp_total = (n_int + n_shr) * nqp_per_face;

   map.qp_to_elem_plus   .assign(n_qp_total, -1);
   map.qp_to_dof_plus    .assign(n_qp_total, -1);
   map.qp_to_elem_minus  .assign(n_qp_total, -1);
   map.qp_to_dof_minus   .assign(n_qp_total, -1);
   map.qp_shape_max_plus .assign(n_qp_total, 0.0);
   map.qp_shape_max_minus.assign(n_qp_total, 0.0);

   // R-I06-008: use a ref-normal-aware side discriminator.  An element
   // is on the canonical "+" side when (elem_centroid - face_centroid)
   // projects onto -ref_normal with a positive value (i.e., the element
   // lies opposite to the direction ref_normal points).  For TPV102
   // (ref_normal = (0,-1,0), fault at y=0) this reduces to the original
   // "centroid_y > 0" check.  Caller may override ref_normal for
   // dipping or non-planar fault geometries.
   const real_t default_ref_n[3] = {0.0, -1.0, 0.0};
   const real_t *ref_n = ref_normal ? ref_normal : default_ref_n;

   auto elem_on_plus_side = [&mesh, ref_n](int elem_id, int face_id) -> bool
   {
      Array<int> ev;
      mesh.GetElementVertices(elem_id, ev);
      real_t ce[3] = {0.0, 0.0, 0.0};
      for (int v = 0; v < ev.Size(); v++)
      {
         const real_t *vpos = mesh.GetVertex(ev[v]);
         for (int d = 0; d < 3; d++) { ce[d] += vpos[d]; }
      }
      const int nv = std::max(1, ev.Size());
      for (int d = 0; d < 3; d++) { ce[d] /= nv; }

      Array<int> fv;
      mesh.GetFaceVertices(face_id, fv);
      real_t cf[3] = {0.0, 0.0, 0.0};
      for (int v = 0; v < fv.Size(); v++)
      {
         const real_t *vpos = mesh.GetVertex(fv[v]);
         for (int d = 0; d < 3; d++) { cf[d] += vpos[d]; }
      }
      const int nf = std::max(1, fv.Size());
      for (int d = 0; d < 3; d++) { cf[d] /= nf; }

      const real_t dot = (ce[0] - cf[0]) * (-ref_n[0])
                       + (ce[1] - cf[1]) * (-ref_n[1])
                       + (ce[2] - cf[2]) * (-ref_n[2]);
      return dot > 0.0;
   };

   int qp_offset = 0;
   // Interior fault faces
   for (int fi = 0; fi < n_int; fi++)
   {
      const int face = fault_int_faces[fi];
      auto *ftr = mesh.GetInteriorFaceTransformations(face);
      MFEM_VERIFY(ftr, "BuildFaultQPNodalMap: interior fault face " << face
                  << " has no transformation");

      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      MFEM_VERIFY(e1 >= 0 && e2 >= 0,
                  "BuildFaultQPNodalMap: interior fault face " << face
                  << " missing an element (Elem1=" << e1 << ", Elem2=" << e2
                  << ")");

      const bool e1_is_plus = elem_on_plus_side(e1, face);
      const int elem_plus  = e1_is_plus ? e1 : e2;
      const int elem_minus = e1_is_plus ? e2 : e1;

      const FiniteElement *fe_plus  = fes.GetFE(elem_plus);
      const FiniteElement *fe_minus = fes.GetFE(elem_minus);

      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      MFEM_VERIFY(ir.GetNPoints() == nqp_per_face,
                  "BuildFaultQPNodalMap: nqp_per_face mismatch on face " << face);

      for (int q = 0; q < nqp_per_face; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         IntegrationPoint ip_plus, ip_minus;
         if (e1_is_plus)
         {
            ftr->Loc1.Transform(ip, ip_plus);
            ftr->Loc2.Transform(ip, ip_minus);
         }
         else
         {
            ftr->Loc2.Transform(ip, ip_plus);
            ftr->Loc1.Transform(ip, ip_minus);
         }

         real_t sh_p, sh_m;
         const int d_plus  = PickNearestNodalDOF(*fe_plus,  ip_plus,  sh_p);
         const int d_minus = PickNearestNodalDOF(*fe_minus, ip_minus, sh_m);
         MFEM_VERIFY(sh_p > kShapeMaxFloor,
                     "BuildFaultQPNodalMap: + side shape_max=" << sh_p
                     << " at or below " << kShapeMaxFloor
                     << " at QP " << q << " of face " << face
                     << ".  The ApplyNucleationTotal 1/shape_max scaling "
                     "would diverge.  Require GaussLobatto (or any nodal-"
                     "interpolatory) basis; modal / hierarchical bases "
                     "need the deferred L2-projection path (plan Phase 5 "
                     "Risk Assessment).");
         MFEM_VERIFY(sh_m > kShapeMaxFloor,
                     "BuildFaultQPNodalMap: - side shape_max=" << sh_m
                     << " at or below " << kShapeMaxFloor
                     << " at QP " << q << " of face " << face);

         const int qp_idx = qp_offset + q;
         map.qp_to_elem_plus    [qp_idx] = elem_plus;
         map.qp_to_dof_plus     [qp_idx] = d_plus;
         map.qp_to_elem_minus   [qp_idx] = elem_minus;
         map.qp_to_dof_minus    [qp_idx] = d_minus;
         map.qp_shape_max_plus  [qp_idx] = sh_p;
         map.qp_shape_max_minus [qp_idx] = sh_m;
      }
      qp_offset += nqp_per_face;
   }

#ifdef MFEM_USE_MPI
   // Shared fault faces (ParMesh only) — Elem1 is local, Elem2 is ghost.
   if constexpr (IsParallelMesh<MeshType>::value)
   {
      auto &pmesh = static_cast<ParMesh &>(mesh);
      for (int si = 0; si < n_shr; si++)
      {
         const int sf = fault_shr_faces[si];
         auto *ftr = pmesh.GetSharedFaceTransformations(sf);
         MFEM_VERIFY(ftr, "BuildFaultQPNodalMap: shared fault face " << sf
                     << " has no transformation");

         const int e1 = ftr->Elem1No;
         MFEM_VERIFY(e1 >= 0,
                     "BuildFaultQPNodalMap: shared fault face " << sf
                     << " has no local Elem1");

         const bool e1_is_plus = elem_on_plus_side(e1, sf);
         const FiniteElement *fe1 = fes.GetFE(e1);

         const IntegrationRule &ir =
            IntRules.Get(ftr->GetGeometryType(), 2*order);
         MFEM_VERIFY(ir.GetNPoints() == nqp_per_face,
                     "BuildFaultQPNodalMap: nqp_per_face mismatch on shared "
                     "face " << sf);

         for (int q = 0; q < nqp_per_face; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            IntegrationPoint ip1;
            ftr->Loc1.Transform(ip, ip1);
            real_t sh1;
            const int d1 = PickNearestNodalDOF(*fe1, ip1, sh1);
            MFEM_VERIFY(sh1 > kShapeMaxFloor,
                        "BuildFaultQPNodalMap: shared-face shape_max=" << sh1
                        << " at or below " << kShapeMaxFloor
                        << " at shared face " << sf << " QP " << q);

            const int qp_idx = qp_offset + q;
            if (e1_is_plus)
            {
               map.qp_to_elem_plus   [qp_idx] = e1;
               map.qp_to_dof_plus    [qp_idx] = d1;
               map.qp_shape_max_plus [qp_idx] = sh1;
               // minus-side update happens on the peer rank
            }
            else
            {
               map.qp_to_elem_minus   [qp_idx] = e1;
               map.qp_to_dof_minus    [qp_idx] = d1;
               map.qp_shape_max_minus [qp_idx] = sh1;
            }
         }
         qp_offset += nqp_per_face;
      }
   }
#else
   (void)fault_shr_faces;
#endif

   MFEM_VERIFY(qp_offset == n_qp_total,
               "BuildFaultQPNodalMap: populated " << qp_offset
               << " QPs but expected " << n_qp_total);

   // R-I06-round3 R-002: verify that no (elem, dof) pair appears on
   // the + side of more than one QP, and likewise for the - side.
   // See VerifyFaultQPNodalMapUnique for the rationale and abort message.
   VerifyFaultQPNodalMapUnique(map);
}

/// @brief Mutable state accompanying a FaultQPNodalMap: the nucleation
/// amplitude (dtau) most recently applied to each fault QP's nearest
/// nodal DOF on the + and - sides.  Driver owns; ApplyNucleationTotal
/// consumes and mutates.  Addresses REVIEW.md R-I06-004 (the overwrite
/// → add-delta migration).
struct FaultQPNucleationState
{
   /// dtau already applied to Q[SXY, DOF_plus[i]] / Q[SXY, DOF_minus[i]].
   /// Initial 0.0 (InitializeStateTotal lays down Q[SXY] = -tau_ini at
   /// every DOF, which corresponds to dtau=0).  Updated by
   /// ApplyNucleationTotal on every call.  Size = number of fault QPs;
   /// entries for non-local sides are also tracked (kept at 0.0 — their
   /// peer rank tracks the applied dtau there).
   std::vector<real_t> dtau_applied_plus;
   std::vector<real_t> dtau_applied_minus;

   /// Reset to all-zero.  Size to n_qp.  Call once at driver init,
   /// immediately after InitializeStateTotal.
   void Reset(int n_qp)
   {
      dtau_applied_plus .assign(n_qp, 0.0);
      dtau_applied_minus.assign(n_qp, 0.0);
   }
};

// v9.4.0 R-006 (review): bulk-Q `ApplyNucleationTotal` removed.
// Rationale: it poked nucleation as a single-DOF point source into
// Q[SXY] at the fault QP, then relied on the wave operator's
// shape-weighted QP reconstruction to surface it.  The point source
// radiates in O(h/cp), diluting the requested nucleation amplitude
// by 1e4-1e5x by the time the friction solver reads it.  Replaced
// by `ApplyNucleationPrestress` (DOFData persistent channel, below)
// and the supporting FaultQPNodalMap / FaultQPNucleationState
// infrastructure above is retained only for legacy-test access;
// under v9.4.0 production dispatch nothing reads it.

/// @brief Persistent-prestress nucleation channel for TPV102
///        (production driver's only nucleation entry point after v9.4.0).
///
/// OVERWRITES the per-DOF nucleation prestress each call.
/// `FaultFaceFlux::Evaluate` (v9.4.0 Commit 1) reads `data.tau2_nuc` and
/// adds it to the friction-input traction, so the requested `dtau` is
/// re-imposed at every Riemann solve and never gets diluted by wave
/// radiation.
///
/// Why this replaces the bulk-Q `ApplyNucleationTotal` (DEPRECATED):
///   The bulk-Q version pokes a single nodal DOF in `Q[SXY]` with
///   `-delta/shape_max`.  The wave operator at `wave_operator.inl:905-913`
///   reads Q via `Q_self[c] = sum_i shape1(i) * Q[..., dof_i]`, treating
///   the poke as a point source.  For a 3D DG wave system that point
///   source radiates outward in `O(h/cp) ≈ 3×10⁻⁵ s` on a 200 m mesh —
///   far shorter than the 1 s nucleation ramp.  The shape-weighted QP
///   amplitude effective in `ComputeTrialTraction` ends up ~10⁴-10⁵×
///   below the requested `nuc_dtau = 25 MPa` (back-of-envelope: ramp
///   rate × residence time = 25 MPa/s × 30 µs ≈ 750 Pa).  See
///   `debug_document/tpv102_debug_document/tpv102_nucleation_code_review_2026-04-22.md`
///   and the companion fix doc `tpv102_nucleation_code_fix_2026-04-22.md`.
///
/// The persistent-prestress channel sidesteps the wave-radiation problem
/// entirely: nucleation does not enter bulk Q at all.  At every
/// Evaluate call the friction solver sees
///     tau2_total = tau2_0 + data.tau2_nuc + tau2_trial(Q)
/// where `tau2_nuc` is overwritten per call from this function.  Bulk Q
/// at the fault DOFs evolves naturally under the slip-driven Riemann
/// imposed-state coupling.
///
/// BP5 convention (tangent2 = strike): TPV102's along-strike pre-stress
/// lives in component 2.  This function therefore writes to `tau2_nuc`,
/// leaving `tau1_nuc` and `sigma_n_nuc` at their default 0.
///
/// @param[in,out] dof_data      Fault DOFData (`tau2_nuc` is the only
///                              field mutated).  Sized to one entry per
///                              fault QP (matches `InitializeFaultDOFs`).
/// @param[in]     fault_coords  Per-QP physical (x, y, z); `x` =
///                              along-strike, `|z|` = down-dip.
/// @param[in]     t             Simulation time [s].  The per-QP
///                              nucleation amplitude is
///                              `NucleationPerturbation(x, |z|, t)`.
inline void ApplyNucleationPrestress(std::vector<DOFData> &dof_data,
                                     const std::vector<Vector> &fault_coords,
                                     real_t t)
{
   const int n = static_cast<int>(dof_data.size());
   MFEM_VERIFY(static_cast<int>(fault_coords.size()) >= n,
               "ApplyNucleationPrestress: fault_coords size "
               << fault_coords.size() << " < dof_data size " << n);
   for (int i = 0; i < n; i++)
   {
      const real_t along_strike = fault_coords[i](0);
      const real_t down_dip     = std::abs(fault_coords[i](2));
      const real_t dtau         = NucleationPerturbation(along_strike,
                                                          down_dip, t);
      dof_data[i].tau2_nuc = dtau;
      // tau1_nuc / sigma_n_nuc intentionally left at their default 0
      // (TPV102 is pure strike-slip, no normal-stress nucleation).
   }
}

/// @brief DEPRECATED alias preserved for out-of-tree callers.
///
/// v9.4.0 renames `ApplyNucleationTotalPrestress` to
/// `ApplyNucleationPrestress` (the function no longer has anything
/// "total"-specific to it — it writes a persistent channel read by the
/// fluctuation-path `Evaluate`).  In-tree callers migrate in Commit 2;
/// this alias remains solely so that any external driver that never got
/// updated still compiles.  Prefer the new name for all new code.
inline void ApplyNucleationTotalPrestress(std::vector<DOFData> &dof_data,
                                          const std::vector<Vector> &fault_coords,
                                          real_t t)
{
   ApplyNucleationPrestress(dof_data, fault_coords, t);
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV102_SETUP_TOTAL_HPP
