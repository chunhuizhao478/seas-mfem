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

#include "precomputed_face_fluxes.hpp"

// Full definition of FreeSurfaceBCMode (forward-declared in the header
// to break a circular include chain).  Only the .cpp needs the full
// enum values because the header parameter declaration is satisfied by
// the underlying-type-specified forward declaration alone.
#include "wave_operator.hpp"

#include <cmath>
#include <cstring>

namespace mfem
{
namespace seas
{

namespace
{

inline void Cross3(const real_t a[3], const real_t b[3], real_t out[3])
{
   out[0] = a[1]*b[2] - a[2]*b[1];
   out[1] = a[2]*b[0] - a[0]*b[2];
   out[2] = a[0]*b[1] - a[1]*b[0];
}

inline real_t Dot3(const real_t a[3], const real_t b[3])
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline real_t Norm3(const real_t a[3])
{
   return std::sqrt(Dot3(a, a));
}

inline void Normalize3(real_t a[3])
{
   real_t n = Norm3(a);
   MFEM_VERIFY(n > 0.0,
               "PrecomputedFaceFluxes: degenerate (zero-length) vector; "
               "topology frame construction cannot continue.");
   a[0] /= n;
   a[1] /= n;
   a[2] /= n;
}

// Dimensional identity check is commented in BuildInteriorMatrices.
// Per plan §5.1.6:
//   w * nApNm1 * I_self  ==  ip.weight * nor_len * T * A_plus * Tinv * I_self
//                        ==  MFEM's `w * flux_.Interior(nor_unit, ...)_self`
// so `nor_len` is NOT baked into `nApNm1`; the integration weight is
// entirely at the runtime call site.  See plan §5.1.6.

// |det(3x3)| for three rows given as real_t[3].
inline real_t Det3x3(const real_t r0[3], const real_t r1[3], const real_t r2[3])
{
   return r0[0] * (r1[1]*r2[2] - r1[2]*r2[1])
        - r0[1] * (r1[0]*r2[2] - r1[2]*r2[0])
        + r0[2] * (r1[0]*r2[1] - r1[1]*r2[0]);
}

} // anonymous namespace

// --------------------------------------------------------------------------
// AssertMFEMFace2NodesTable
// --------------------------------------------------------------------------
void PrecomputedFaceFluxes::AssertMFEMFace2NodesTable()
{
   using Tet = mfem::Geometry::Constants<mfem::Geometry::TETRAHEDRON>;
   for (int s = 0; s < 4; s++)
   {
      for (int j = 0; j < 3; j++)
      {
         MFEM_VERIFY(FACE2NODES_MFEM[s][j] == Tet::FaceVert[s][j],
                     "PrecomputedFaceFluxes: MFEM tet face-vertex table "
                     "drift at side " << s << " vertex " << j
                     << " (MFEM=" << Tet::FaceVert[s][j]
                     << " vs hard-coded=" << FACE2NODES_MFEM[s][j] << ")");
      }
   }
}

// --------------------------------------------------------------------------
// ComputeCellFaceFrame  (plan §5.1.4)
// --------------------------------------------------------------------------
void PrecomputedFaceFluxes::ComputeCellFaceFrame(
   const Mesh &mesh,
   int elem, int local_side,
   real_t normal[3],
   real_t tangent1[3], real_t tangent2[3],
   real_t &surface_area,
   real_t &cell_volume)
{
   Array<int> ev;
   mesh.GetElementVertices(elem, ev);
   MFEM_VERIFY(ev.Size() == 4,
               "PrecomputedFaceFluxes: expected 4 element vertices (tet); "
               "got " << ev.Size());
   MFEM_VERIFY(local_side >= 0 && local_side < 4,
               "PrecomputedFaceFluxes: local_side out of range "
               "[0, 4): " << local_side);

   real_t p_elem[4][3];
   for (int i = 0; i < 4; i++)
   {
      const real_t *v = mesh.GetVertex(ev[i]);
      p_elem[i][0] = v[0];
      p_elem[i][1] = v[1];
      p_elem[i][2] = v[2];
   }

   int f0 = FACE2NODES_MFEM[local_side][0];
   int f1 = FACE2NODES_MFEM[local_side][1];
   int f2 = FACE2NODES_MFEM[local_side][2];

   real_t p0[3] = {p_elem[f0][0], p_elem[f0][1], p_elem[f0][2]};
   real_t p1[3] = {p_elem[f1][0], p_elem[f1][1], p_elem[f1][2]};
   real_t p2[3] = {p_elem[f2][0], p_elem[f2][1], p_elem[f2][2]};

   real_t ab[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
   real_t ac[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};

   real_t normal_raw[3];
   Cross3(ab, ac, normal_raw);

   real_t c_elem[3] = {
      (p_elem[0][0] + p_elem[1][0] + p_elem[2][0] + p_elem[3][0]) / 4.0,
      (p_elem[0][1] + p_elem[1][1] + p_elem[2][1] + p_elem[3][1]) / 4.0,
      (p_elem[0][2] + p_elem[1][2] + p_elem[2][2] + p_elem[3][2]) / 4.0
   };
   real_t c_face[3] = {
      (p0[0] + p1[0] + p2[0]) / 3.0,
      (p0[1] + p1[1] + p2[1]) / 3.0,
      (p0[2] + p1[2] + p2[2]) / 3.0
   };
   real_t outward[3] = { c_face[0]-c_elem[0],
                         c_face[1]-c_elem[1],
                         c_face[2]-c_elem[2] };
   if (Dot3(normal_raw, outward) < 0.0)
   {
      // Swap fn[0] <-> fn[1] (i.e. swap p0 and p1) and recompute.
      real_t tmp[3] = {p0[0], p0[1], p0[2]};
      p0[0] = p1[0]; p0[1] = p1[1]; p0[2] = p1[2];
      p1[0] = tmp[0]; p1[1] = tmp[1]; p1[2] = tmp[2];
      ab[0] = p1[0]-p0[0]; ab[1] = p1[1]-p0[1]; ab[2] = p1[2]-p0[2];
      ac[0] = p2[0]-p0[0]; ac[1] = p2[1]-p0[1]; ac[2] = p2[2]-p0[2];
      Cross3(ab, ac, normal_raw);
      MFEM_VERIFY(Dot3(normal_raw, outward) > 0.0,
                  "PrecomputedFaceFluxes::ComputeCellFaceFrame: "
                  "outward orientation not recovered after vertex swap "
                  "(elem=" << elem << " side=" << local_side << ")");
   }

   real_t tangent1_raw[3] = {ab[0], ab[1], ab[2]};
   real_t tangent2_raw[3];
   Cross3(normal_raw, tangent1_raw, tangent2_raw);

   real_t n_len = Norm3(normal_raw);
   surface_area = 0.5 * n_len;

   normal[0] = normal_raw[0]; normal[1] = normal_raw[1]; normal[2] = normal_raw[2];
   Normalize3(normal);
   tangent1[0] = tangent1_raw[0]; tangent1[1] = tangent1_raw[1]; tangent1[2] = tangent1_raw[2];
   Normalize3(tangent1);
   tangent2[0] = tangent2_raw[0]; tangent2[1] = tangent2_raw[1]; tangent2[2] = tangent2_raw[2];
   Normalize3(tangent2);

   // Cell volume: rows of the coordinate-difference matrix, NOT vertex
   // indices.  Unit-tet sanity: |det| / 6 == 1/6 to 1e-14.
   real_t dX0[3] = {p_elem[1][0]-p_elem[0][0],
                    p_elem[1][1]-p_elem[0][1],
                    p_elem[1][2]-p_elem[0][2]};
   real_t dX1[3] = {p_elem[2][0]-p_elem[0][0],
                    p_elem[2][1]-p_elem[0][1],
                    p_elem[2][2]-p_elem[0][2]};
   real_t dX2[3] = {p_elem[3][0]-p_elem[0][0],
                    p_elem[3][1]-p_elem[0][1],
                    p_elem[3][2]-p_elem[0][2]};
   cell_volume = std::abs(Det3x3(dX0, dX1, dX2)) / 6.0;
}

// --------------------------------------------------------------------------
// BuildInteriorMatrices  (plan §5.1.5)
// --------------------------------------------------------------------------
void PrecomputedFaceFluxes::BuildInteriorMatrices(
   const real_t normal[3],
   const real_t tangent1[3],
   const real_t tangent2[3],
   const GodunovFlux &flux,
   DenseMatrix &nApNm1,
   DenseMatrix &nAmNm1)
{
   DenseMatrix T(NUM_STATE, NUM_STATE);
   DenseMatrix Tinv(NUM_STATE, NUM_STATE);
   GodunovFlux::BuildRotation(normal, tangent1, tangent2, T);
   GodunovFlux::BuildRotationInverse(normal, tangent1, tangent2, Tinv);

   const DenseMatrix &A_plus  = flux.GetAxPlus();
   const DenseMatrix &A_minus = flux.GetAxMinus();

   DenseMatrix tmp(NUM_STATE, NUM_STATE);

   nApNm1.SetSize(NUM_STATE, NUM_STATE);
   Mult(A_plus, Tinv, tmp);
   Mult(T, tmp, nApNm1);

   nAmNm1.SetSize(NUM_STATE, NUM_STATE);
   Mult(A_minus, Tinv, tmp);
   Mult(T, tmp, nAmNm1);
}

// --------------------------------------------------------------------------
// Boundary matrix probes  (plan §5.2.1 / §5.2.2 / §5.2.3)
// --------------------------------------------------------------------------
static void ProbeBoundary9x9(
   void (*probe)(const GodunovFlux&, const real_t*, const real_t*,
                 const real_t*, real_t*),
   const real_t normal[3],
   const GodunovFlux &flux,
   DenseMatrix &nApNm1)
{
   real_t bulk_bg_zero[NUM_STATE] = {0.0};
   nApNm1.SetSize(NUM_STATE, NUM_STATE);
   for (int j = 0; j < NUM_STATE; j++)
   {
      real_t I_self[NUM_STATE] = {0.0};
      I_self[j] = 1.0;
      real_t F_h[NUM_STATE] = {0.0};
      probe(flux, normal, I_self, bulk_bg_zero, F_h);
      for (int i = 0; i < NUM_STATE; i++) { nApNm1(i, j) = F_h[i]; }
   }
}

static void ProbeGodunovTotal(const GodunovFlux &flux,
                              const real_t *n, const real_t *I_self,
                              const real_t *Q_bg, real_t *F_h)
{
   flux.FreeSurfaceGodunovTotal(n, I_self, Q_bg, F_h);
}

static void ProbeFreeSurfaceTotal(const GodunovFlux &flux,
                                  const real_t *n, const real_t *I_self,
                                  const real_t *Q_bg, real_t *F_h)
{
   flux.FreeSurfaceTotal(n, I_self, Q_bg, F_h);
}

static void ProbeAbsorbingTotal(const GodunovFlux &flux,
                                const real_t *n, const real_t *I_self,
                                const real_t *Q_bg, real_t *F_h)
{
   flux.AbsorbingTotal(n, I_self, Q_bg, F_h);
}

void PrecomputedFaceFluxes::BuildBoundaryMatricesGodunov(
   const real_t normal[3], const GodunovFlux &flux, DenseMatrix &nApNm1)
{
   ProbeBoundary9x9(&ProbeGodunovTotal, normal, flux, nApNm1);
}

void PrecomputedFaceFluxes::BuildBoundaryMatricesGamma(
   const real_t normal[3], const GodunovFlux &flux, DenseMatrix &nApNm1)
{
   ProbeBoundary9x9(&ProbeFreeSurfaceTotal, normal, flux, nApNm1);
}

void PrecomputedFaceFluxes::BuildBoundaryMatricesAbsorbing(
   const real_t normal[3], const GodunovFlux &flux, DenseMatrix &nApNm1)
{
   ProbeBoundary9x9(&ProbeAbsorbingTotal, normal, flux, nApNm1);
}

// --------------------------------------------------------------------------
// Init  (plan §5.3, post-R3/R4 fixes)
// --------------------------------------------------------------------------
void PrecomputedFaceFluxes::Init(
   Mesh &mesh,
   const FiniteElementSpace &fes,
   const BoundaryConfig &bc,
   const GodunovFlux &flux,
   FreeSurfaceBCMode fs_mode,
   const std::set<int> &fault_face_set,
   const std::vector<int> &face_bdr_attr,
   const std::set<int> &shared_mesh_face_set,
   bool build_arm3d_tables)
{
   (void)fes;   // sizing is not stored here; AddInteriorFaceRhs passes ndof_total
   AssertMFEMFace2NodesTable();

   MFEM_VERIFY(face_bdr_attr.size() == static_cast<std::size_t>(mesh.GetNumFaces()),
               "PrecomputedFaceFluxes::Init: face_bdr_attr sized "
               << face_bdr_attr.size() << " but mesh has "
               << mesh.GetNumFaces() << " faces.");

   const int ne = mesh.GetNE();
   entries_.clear();
   face_elem_to_entry_.clear();
   face_to_bdr_entry_.clear();

   // Phase 1 precondition: GodunovFlux carries a SINGLE homogeneous
   // (lambda, mu, rho).  Bimaterial faces cannot arise with this design,
   // so the plan §2.4 "homog_ok" guard is STRUCTURAL (enforced at the
   // type level), not a runtime check.  A naive `homog_ok(flux.GetZp(),
   // flux.GetZp())` comparison would be a tautology (R-001).  When (and
   // if) per-element material support is added, reinstate a per-side
   // check that mirrors fault_face_flux.cpp:95-110's pattern using the
   // extended material lookup (homog_ok(Zp_plus, Zp_minus) with the two
   // sides' distinct impedances).  Until then we avoid a false-positive
   // MFEM_VERIFY here.
   (void)flux;  // acknowledged — no per-cell impedance to compare

   for (int e = 0; e < ne; e++)
   {
      Array<int> el_faces, el_orients;
      mesh.GetElementFaces(e, el_faces, el_orients);
      MFEM_VERIFY(el_faces.Size() == 4,
                  "PrecomputedFaceFluxes::Init: element " << e
                  << " has " << el_faces.Size()
                  << " faces (expected 4; tet-only Phase 1 precondition).");

      for (int s = 0; s < 4; s++)
      {
         const int face_idx = el_faces[s];

         if (shared_mesh_face_set.count(face_idx) > 0) { continue; }
         if (fault_face_set.count(face_idx) > 0)       { continue; }

         FaceElementTransformations *ftr =
            mesh.GetFaceElementTransformations(face_idx);
         if (!ftr) { continue; }

         const int e2_ftr = ftr->Elem2No;
         const int bdr_attr = face_bdr_attr[face_idx];
         const bool is_boundary = (e2_ftr < 0) && (bdr_attr > 0);

         // Non-boundary with no neighbor is an orphan face; skip.
         if (e2_ftr < 0 && bdr_attr == 0) { continue; }

         // R3-009 precondition: face Jacobian constant across QPs (p=1).
         {
            const IntegrationRule &ir = IntRules.Get(
               ftr->GetGeometryType(), 2);
            const int npts = ir.GetNPoints();
            MFEM_VERIFY(npts >= 2,
                        "PrecomputedFaceFluxes::Init: face has fewer than "
                        "two quadrature points; cannot verify constant "
                        "Jacobian.");
            Vector nor_first(mesh.SpaceDimension());
            Vector nor_last(mesh.SpaceDimension());
            {
               const IntegrationPoint &ip0 = ir.IntPoint(0);
               ftr->SetAllIntPoints(&ip0);
               CalcOrtho(ftr->Face->Jacobian(), nor_first);
            }
            {
               const IntegrationPoint &ipl = ir.IntPoint(npts-1);
               ftr->SetAllIntPoints(&ipl);
               CalcOrtho(ftr->Face->Jacobian(), nor_last);
            }
            const real_t len0 = nor_first.Norml2();
            const real_t lenN = nor_last.Norml2();
            MFEM_VERIFY(len0 > 0.0,
                        "PrecomputedFaceFluxes::Init: zero-length face "
                        "Jacobian on face " << face_idx);
            MFEM_VERIFY(std::abs(len0 - lenN) <= 1e-12 * len0,
                        "PrecomputedFaceFluxes requires p=1 isoparametric "
                        "(affine) mesh; face " << face_idx << " has "
                        "QP-varying Jacobian (face is curved): |len_first - "
                        "len_last| = " << std::abs(len0 - lenN)
                        << " vs len_first=" << len0);
         }

         FaceEntry fe;
         fe.element    = e;
         fe.local_side = s;
         fe.face_idx   = face_idx;

         ComputeCellFaceFrame(mesh, e, s,
                              fe.normal, fe.tangent1, fe.tangent2,
                              fe.surface_area, fe.cell_volume);

         // === Arm 3d: shape-table precomputation via MFEM's Loc1/Loc2 ====
         // The orbit-asymmetry that Arm 2b diagnosed is in the FACE NORMAL
         // from CalcOrtho, not in the face-ref-to-elem-ref QP mapping.
         // Loc1.Transform / Loc2.Transform are orbit-covariant on their
         // own — they map face ref IP to elem ref IP via MFEM's internal
         // orientation-aware tables.  Using Loc1/Loc2 here guarantees
         // shape values match MFEM's runtime sampling EXACTLY on every
         // mesh, while the precomputed nApNm1/nAmNm1 matrices come from
         // topology-based (n_outward, t1, t2) — the SINGLE place the
         // orbit asymmetry was injected.
         //
         // This is the correct reading of Phase 3 STOP findings (Arm 1 +
         // Arm 2 + Arm 2b + Arm 3d):
         //   * MFEM's kernels + Loc1/Loc2 are D4-covariant at sampling.
         //   * CalcOrtho's stored face normal is NOT D4-covariant on
         //     orbit-pair diagonal faces.
         //   * The fix is: replace CalcOrtho-based normal with
         //     topology-based normal (done in BuildInteriorMatrices);
         //     keep Loc1/Loc2-based QP sampling (done here).
         //
         // R-003 (v9.5.0): gated behind `build_arm3d_tables`.  The Full
         // path this table feeds is not wired into production dispatch
         // (see wave_operator.inl:2101 note) pending plane-wave
         // equivalence investigation; leaving the table resident in
         // every FaceEntry wastes memory proportional to
         // n_entries × nqp × ndof_per_el × 2 on production runs that
         // will never read it.  Arm 3d unit tests opt in explicitly.
         if (build_arm3d_tables)
         {
            const FiniteElement *fe_self = fes.GetFE(e);
            const int ndof_self = fe_self->GetDof();

            const IntegrationRule &ir = IntRules.Get(
               ftr->GetGeometryType(), 2);
            const int nqp = ir.GetNPoints();

            fe.nqp = nqp;
            fe.w_qp.resize(nqp);
            fe.shape_self.assign(nqp, std::vector<real_t>(ndof_self, 0.0));
            fe.dof_offset_self = e * ndof_self;
            fe.ndof_per_el     = ndof_self;

            int nbr_elem_id = -1;
            if (!is_boundary)
            {
               nbr_elem_id = (ftr->Elem1No == e) ? ftr->Elem2No : ftr->Elem1No;
               fe.shape_nbr.assign(nqp, std::vector<real_t>(ndof_self, 0.0));
               fe.dof_offset_nbr = nbr_elem_id * ndof_self;
            }

            Vector shape_eval(ndof_self);

            // Which MFEM Loc to use for "self" depends on whether e is
            // Elem1 or Elem2 of this face.  MFEM sets Elem1 to the cell
            // whose outward-normal alignment is canonical; we use Loc1
            // when e == Elem1No, else Loc2.
            const bool self_is_elem1 = (ftr->Elem1No == e);

            for (int q = 0; q < nqp; q++)
            {
               const IntegrationPoint &ipq = ir.IntPoint(q);
               ftr->SetAllIntPoints(&ipq);

               // nor_len: use MFEM's face-Jacobian magnitude via CalcOrtho
               // exactly as the runtime does.  For affine faces this
               // equals 2 * surface_area_topology to ULP.
               Vector nor_vec(mesh.SpaceDimension());
               CalcOrtho(ftr->Face->Jacobian(), nor_vec);
               fe.w_qp[q] = ipq.weight * nor_vec.Norml2();

               IntegrationPoint ip_self, ip_nbr;
               if (self_is_elem1)
               { ftr->Loc1.Transform(ipq, ip_self); }
               else
               { ftr->Loc2.Transform(ipq, ip_self); }
               fe_self->CalcShape(ip_self, shape_eval);
               for (int i = 0; i < ndof_self; i++)
               { fe.shape_self[q][i] = shape_eval(i); }

               if (!is_boundary)
               {
                  if (self_is_elem1)
                  { ftr->Loc2.Transform(ipq, ip_nbr); }
                  else
                  { ftr->Loc1.Transform(ipq, ip_nbr); }
                  const FiniteElement *fe_nbr = fes.GetFE(nbr_elem_id);
                  fe_nbr->CalcShape(ip_nbr, shape_eval);
                  for (int i = 0; i < ndof_self; i++)
                  { fe.shape_nbr[q][i] = shape_eval(i); }
               }
            }
         }
         // === end Arm 3d ==================================================

         if (is_boundary)
         {
            FacePrecomputedBC bc_class;
            if (bc.natural_attrs.count(bdr_attr) > 0)
            {
               bc_class = (fs_mode == FreeSurfaceBCMode::Godunov)
                          ? FacePrecomputedBC::FreeSurfaceGodunov
                          : FacePrecomputedBC::FreeSurfaceGamma;
            }
            else if (bc.absorbing_attrs.count(bdr_attr) > 0)
            {
               bc_class = FacePrecomputedBC::Absorbing;
            }
            else if (bdr_attr == bc.fault_attr)
            {
               MFEM_ABORT("PrecomputedFaceFluxes::Init: face " << face_idx
                          << " is 1-sided with fault_attr; should have "
                          "been skipped via fault_face_set.");
            }
            else
            {
               MFEM_ABORT("PrecomputedFaceFluxes::Init: face " << face_idx
                          << " has unclassified boundary attribute "
                          << bdr_attr
                          << " (not in natural_attrs, absorbing_attrs, "
                          "or == fault_attr).");
            }

            fe.bc = bc_class;
            fe.neighbor_elem = -1;

            switch (bc_class)
            {
               case FacePrecomputedBC::FreeSurfaceGodunov:
                  BuildBoundaryMatricesGodunov(fe.normal, flux, fe.nApNm1);
                  break;
               case FacePrecomputedBC::FreeSurfaceGamma:
                  BuildBoundaryMatricesGamma(fe.normal, flux, fe.nApNm1);
                  break;
               case FacePrecomputedBC::Absorbing:
                  BuildBoundaryMatricesAbsorbing(fe.normal, flux, fe.nApNm1);
                  break;
               default:
                  MFEM_ABORT("PrecomputedFaceFluxes::Init: unexpected "
                             "FacePrecomputedBC classification.");
            }

            fe.nAmNm1.SetSize(NUM_STATE, NUM_STATE);
            fe.nAmNm1 = 0.0;
         }
         else
         {
            // Interior non-fault face.  Bimaterial guard per plan §2.4
            // is STRUCTURAL (single GodunovFlux, single material) and
            // documented once at the top of Init; no per-face runtime
            // check is possible at this layer (R-001).
            fe.bc = FacePrecomputedBC::Interior;
            fe.neighbor_elem =
               (ftr->Elem1No == e) ? ftr->Elem2No : ftr->Elem1No;

            BuildInteriorMatrices(fe.normal, fe.tangent1, fe.tangent2,
                                  flux, fe.nApNm1, fe.nAmNm1);
         }

         entries_.push_back(std::move(fe));
      }
   }

   // Build face_elem_to_entry_ after entries_ is stable.
   for (std::size_t i = 0; i < entries_.size(); i++)
   {
      const FaceEntry &fe = entries_[i];
      MFEM_VERIFY(fe.element >= 0,
                  "PrecomputedFaceFluxes::Init: Phase 1 requires "
                  "non-negative caller_elem ids.");
      long long key = EncodeKey(fe.face_idx, fe.element);
      auto res = face_elem_to_entry_.emplace(key, static_cast<int>(i));
      MFEM_VERIFY(res.second,
                  "PrecomputedFaceFluxes::Init: duplicate (face_idx, "
                  "element) entry at face " << fe.face_idx
                  << " element " << fe.element);

      if (fe.bc != FacePrecomputedBC::Interior)
      {
         auto bres = face_to_bdr_entry_.emplace(fe.face_idx,
                                                 static_cast<int>(i));
         MFEM_VERIFY(bres.second,
                     "PrecomputedFaceFluxes::Init: duplicate boundary "
                     "entry at face " << fe.face_idx);
      }
   }

   initialized_ = true;
}

#ifdef MFEM_USE_MPI
// --------------------------------------------------------------------------
// InitSharedFaces  (plan §7.2 — Phase 2b shared-face pass)
// --------------------------------------------------------------------------
// One FaceEntry per shared non-fault face per rank, keyed by
// (mesh_face_idx, Elem1No).  Elem1No is always >= 0 (the local owning
// cell's id); the remote neighbor's pseudo-id is never used as a key
// because this rank does not own those DOFs.  The peer rank stores its
// own side's entry independently.
void PrecomputedFaceFluxes::InitSharedFaces(
   ParMesh &pmesh,
   const FiniteElementSpace &fes,
   const GodunovFlux &flux,
   const std::set<int> &fault_face_set)
{
   (void)fes;
   MFEM_VERIFY(initialized_,
               "PrecomputedFaceFluxes::InitSharedFaces: call Init(...) "
               "first; shared-face pass APPENDS to Phase 1 tables.");

   const int n_shared = pmesh.GetNSharedFaces();
   for (int sf = 0; sf < n_shared; sf++)
   {
      const int mesh_face_idx = pmesh.GetSharedFace(sf);

      // Skip shared-fault faces; they are routed through FaultFaceFlux.
      if (fault_face_set.count(mesh_face_idx) > 0) { continue; }

      FaceElementTransformations *ftr = pmesh.GetSharedFaceTransformations(sf);
      if (!ftr) { continue; }

      const int e1 = ftr->Elem1No;
      MFEM_VERIFY(e1 >= 0,
                  "PrecomputedFaceFluxes::InitSharedFaces: shared face "
                  << sf << " (mesh face " << mesh_face_idx << ") has "
                  "Elem1No=" << e1 << " < 0; local owning cell must have "
                  "a non-negative id.");

      // Determine which local side of Elem1 this face is.
      Array<int> el_faces, el_orients;
      pmesh.GetElementFaces(e1, el_faces, el_orients);
      int local_side = -1;
      for (int s = 0; s < el_faces.Size(); s++)
      {
         if (el_faces[s] == mesh_face_idx) { local_side = s; break; }
      }
      MFEM_VERIFY(local_side >= 0 && local_side < 4,
                  "PrecomputedFaceFluxes::InitSharedFaces: local side for "
                  "shared face mesh_face_idx=" << mesh_face_idx
                  << " on element " << e1 << " not found in element faces.");

      // R3-009 precondition: face Jacobian constant across QPs (p=1).
      {
         const IntegrationRule &ir = IntRules.Get(
            ftr->GetGeometryType(), 2);
         const int npts = ir.GetNPoints();
         MFEM_VERIFY(npts >= 2,
                     "PrecomputedFaceFluxes::InitSharedFaces: face has "
                     "fewer than two quadrature points.");
         Vector nor_first(pmesh.SpaceDimension());
         Vector nor_last(pmesh.SpaceDimension());
         {
            const IntegrationPoint &ip0 = ir.IntPoint(0);
            ftr->SetAllIntPoints(&ip0);
            CalcOrtho(ftr->Face->Jacobian(), nor_first);
         }
         {
            const IntegrationPoint &ipl = ir.IntPoint(npts-1);
            ftr->SetAllIntPoints(&ipl);
            CalcOrtho(ftr->Face->Jacobian(), nor_last);
         }
         const real_t len0 = nor_first.Norml2();
         const real_t lenN = nor_last.Norml2();
         MFEM_VERIFY(len0 > 0.0,
                     "PrecomputedFaceFluxes::InitSharedFaces: zero-length "
                     "face Jacobian on shared face sf=" << sf);
         MFEM_VERIFY(std::abs(len0 - lenN) <= 1e-12 * len0,
                     "PrecomputedFaceFluxes::InitSharedFaces: p=1 "
                     "precondition violated on shared face sf=" << sf);
      }

      FaceEntry fe;
      fe.element       = e1;
      fe.local_side    = local_side;
      fe.face_idx      = mesh_face_idx;
      fe.bc            = FacePrecomputedBC::Interior;
      fe.neighbor_elem = -1;   // remote pseudo-id not stored

      ComputeCellFaceFrame(pmesh, e1, local_side,
                           fe.normal, fe.tangent1, fe.tangent2,
                           fe.surface_area, fe.cell_volume);

      BuildInteriorMatrices(fe.normal, fe.tangent1, fe.tangent2,
                            flux, fe.nApNm1, fe.nAmNm1);

      const int idx = static_cast<int>(entries_.size());
      entries_.push_back(std::move(fe));

      long long key = EncodeKey(mesh_face_idx, e1);
      auto res = face_elem_to_entry_.emplace(key, idx);
      MFEM_VERIFY(res.second,
                  "PrecomputedFaceFluxes::InitSharedFaces: duplicate "
                  "(face_idx, Elem1No) entry at shared face "
                  << mesh_face_idx << " element " << e1
                  << "; Phase 1 Init should have skipped shared faces.");
   }
}
#endif // MFEM_USE_MPI

// --------------------------------------------------------------------------
// GetEntry / HasEntry
// --------------------------------------------------------------------------
const PrecomputedFaceFluxes::FaceEntry &
PrecomputedFaceFluxes::GetEntry(int face_idx, int caller_elem_id) const
{
   MFEM_VERIFY(initialized_,
               "PrecomputedFaceFluxes::GetEntry: Init has not run.");
   long long key = EncodeKey(face_idx, caller_elem_id);
   auto it = face_elem_to_entry_.find(key);
   MFEM_VERIFY(it != face_elem_to_entry_.end(),
               "PrecomputedFaceFluxes::GetEntry: no entry for "
               "(face_idx=" << face_idx
               << ", caller_elem_id=" << caller_elem_id << ")");
   return entries_[it->second];
}

bool PrecomputedFaceFluxes::HasEntry(int face_idx, int caller_elem_id) const
{
   if (!initialized_) { return false; }
   long long key = EncodeKey(face_idx, caller_elem_id);
   return face_elem_to_entry_.find(key) != face_elem_to_entry_.end();
}

// --------------------------------------------------------------------------
// AddInteriorFaceRhs  (plan §5.1.5 accumulation; §6.2 call site)
// --------------------------------------------------------------------------
void PrecomputedFaceFluxes::AddInteriorFaceRhs(
   int face_idx,
   int caller_elem,
   const real_t *I_self, const real_t *I_nbr,
   real_t w,
   const real_t *shape1,
   int ndof, int dof_offset1,
   int ndof_total,
   Vector &rhs) const
{
   MFEM_ASSERT(initialized_,
               "PrecomputedFaceFluxes::AddInteriorFaceRhs: Init not run.");
   MFEM_ASSERT(I_self != nullptr && I_nbr != nullptr && shape1 != nullptr,
               "PrecomputedFaceFluxes::AddInteriorFaceRhs: null buffer.");
   MFEM_ASSERT(dof_offset1 >= 0 && ndof > 0 && ndof_total >= dof_offset1 + ndof,
               "PrecomputedFaceFluxes::AddInteriorFaceRhs: invalid DOF "
               "offset/sizing (dof_offset1=" << dof_offset1
               << " ndof=" << ndof << " ndof_total=" << ndof_total << ").");

   const FaceEntry &fe = GetEntry(face_idx, caller_elem);
   MFEM_ASSERT(fe.bc == FacePrecomputedBC::Interior,
               "PrecomputedFaceFluxes::AddInteriorFaceRhs: entry is not "
               "interior (face_idx=" << face_idx
               << ", caller_elem=" << caller_elem << ").");

   real_t F_h[NUM_STATE] = {0.0};
   for (int c = 0; c < NUM_STATE; c++)
   {
      real_t acc = 0.0;
      for (int k = 0; k < NUM_STATE; k++)
      {
         acc += fe.nApNm1(c, k) * I_self[k]
              + fe.nAmNm1(c, k) * I_nbr[k];
      }
      F_h[c] = acc;
   }

   real_t *rhs_data = rhs.GetData();
   for (int c = 0; c < NUM_STATE; c++)
   {
      const real_t scale = w * F_h[c];
      real_t *row = rhs_data + c * ndof_total + dof_offset1;
      for (int i = 0; i < ndof; i++)
      {
         row[i] -= scale * shape1[i];
      }
   }
}

// --------------------------------------------------------------------------
// AddBoundaryFaceRhs  (plan §5.2.1 / §5.2.4 P_BFACE_NONZERO_BG_ABORT)
// --------------------------------------------------------------------------
void PrecomputedFaceFluxes::AddBoundaryFaceRhs(
   int face_idx,
   const real_t *I_self,
   const real_t *bulk_bg_scaled,
   real_t w,
   const real_t *shape1,
   int ndof, int dof_offset1,
   int ndof_total,
   Vector &rhs) const
{
   MFEM_ASSERT(initialized_,
               "PrecomputedFaceFluxes::AddBoundaryFaceRhs: Init not run.");
   MFEM_ASSERT(I_self != nullptr && shape1 != nullptr,
               "PrecomputedFaceFluxes::AddBoundaryFaceRhs: null buffer.");
   MFEM_ASSERT(dof_offset1 >= 0 && ndof > 0 && ndof_total >= dof_offset1 + ndof,
               "PrecomputedFaceFluxes::AddBoundaryFaceRhs: invalid DOF "
               "offset/sizing.");

   // Phase 1 probe-linearity precondition: bulk_bg must be bit-zero.
   if (bulk_bg_scaled != nullptr)
   {
      for (int c = 0; c < NUM_STATE; c++)
      {
         MFEM_VERIFY(std::abs(bulk_bg_scaled[c]) <= 1e-15,
                     "PrecomputedFaceFluxes::AddBoundaryFaceRhs: "
                     "Phase 1 requires bulk_bg = 0 (probe precondition); "
                     "got |bulk_bg_scaled[" << c << "]| = "
                     << std::abs(bulk_bg_scaled[c]) << ".");
      }
   }

   auto it = face_to_bdr_entry_.find(face_idx);
   MFEM_ASSERT(it != face_to_bdr_entry_.end(),
               "PrecomputedFaceFluxes::AddBoundaryFaceRhs: no boundary "
               "entry for face_idx=" << face_idx);
   const FaceEntry &fe = entries_[it->second];

   real_t F_h[NUM_STATE] = {0.0};
   for (int c = 0; c < NUM_STATE; c++)
   {
      real_t acc = 0.0;
      for (int k = 0; k < NUM_STATE; k++)
      {
         acc += fe.nApNm1(c, k) * I_self[k];
      }
      F_h[c] = acc;
   }

   real_t *rhs_data = rhs.GetData();
   for (int c = 0; c < NUM_STATE; c++)
   {
      const real_t scale = w * F_h[c];
      real_t *row = rhs_data + c * ndof_total + dof_offset1;
      for (int i = 0; i < ndof; i++)
      {
         row[i] -= scale * shape1[i];
      }
   }
}

// --------------------------------------------------------------------------
// Arm 3d: AddInteriorFaceRhsFull  — topology-based state sampling.
// --------------------------------------------------------------------------
void PrecomputedFaceFluxes::AddInteriorFaceRhsFull(
   int face_idx,
   int caller_elem,
   const real_t *Q_data,
   int ndof_total,
   Vector &rhs) const
{
   MFEM_ASSERT(initialized_,
               "PrecomputedFaceFluxes::AddInteriorFaceRhsFull: Init not run.");
   MFEM_ASSERT(Q_data != nullptr,
               "PrecomputedFaceFluxes::AddInteriorFaceRhsFull: null Q.");

   const FaceEntry &fe = GetEntry(face_idx, caller_elem);
   MFEM_ASSERT(fe.bc == FacePrecomputedBC::Interior,
               "AddInteriorFaceRhsFull: entry is not interior.");
   MFEM_ASSERT(fe.nqp > 0 && fe.ndof_per_el > 0,
               "AddInteriorFaceRhsFull: topology shape tables not built.");
   MFEM_ASSERT(fe.dof_offset_nbr >= 0,
               "AddInteriorFaceRhsFull: neighbor DOF offset invalid.");

   const int ndof = fe.ndof_per_el;
   real_t *rhs_data = rhs.GetData();

   for (int q = 0; q < fe.nqp; q++)
   {
      // Assemble I_self, I_nbr at QP q using topology shape tables.
      real_t I_self[NUM_STATE] = {0.0};
      real_t I_nbr [NUM_STATE] = {0.0};
      const std::vector<real_t> &ss = fe.shape_self[q];
      const std::vector<real_t> &sn = fe.shape_nbr [q];
      for (int c = 0; c < NUM_STATE; c++)
      {
         real_t acc_s = 0.0, acc_n = 0.0;
         for (int i = 0; i < ndof; i++)
         {
            acc_s += ss[i] * Q_data[c * ndof_total + fe.dof_offset_self + i];
            acc_n += sn[i] * Q_data[c * ndof_total + fe.dof_offset_nbr  + i];
         }
         I_self[c] = acc_s;
         I_nbr [c] = acc_n;
      }

      // F_h = nApNm1 * I_self + nAmNm1 * I_nbr  (global frame)
      real_t F_h[NUM_STATE] = {0.0};
      for (int c = 0; c < NUM_STATE; c++)
      {
         real_t acc = 0.0;
         for (int k = 0; k < NUM_STATE; k++)
         {
            acc += fe.nApNm1(c, k) * I_self[k]
                 + fe.nAmNm1(c, k) * I_nbr [k];
         }
         F_h[c] = acc;
      }

      // Accumulate -= w * shape_self(i) * F_h[c] at caller's DOFs.
      const real_t w = fe.w_qp[q];
      for (int c = 0; c < NUM_STATE; c++)
      {
         real_t *row = rhs_data + c * ndof_total + fe.dof_offset_self;
         const real_t scale = w * F_h[c];
         for (int i = 0; i < ndof; i++)
         {
            row[i] -= scale * ss[i];
         }
      }
   }
}

// --------------------------------------------------------------------------
// Arm 3d: AddBoundaryFaceRhsFull  — topology-based self-sample.
// --------------------------------------------------------------------------
void PrecomputedFaceFluxes::AddBoundaryFaceRhsFull(
   int face_idx,
   const real_t *Q_data,
   const real_t *bulk_bg_scaled,
   int ndof_total,
   Vector &rhs) const
{
   MFEM_ASSERT(initialized_,
               "PrecomputedFaceFluxes::AddBoundaryFaceRhsFull: Init not run.");
   MFEM_ASSERT(Q_data != nullptr,
               "PrecomputedFaceFluxes::AddBoundaryFaceRhsFull: null Q.");

   // Phase 1 precondition: bulk_bg must be bit-zero for probe linearity.
   if (bulk_bg_scaled != nullptr)
   {
      for (int c = 0; c < NUM_STATE; c++)
      {
         MFEM_VERIFY(std::abs(bulk_bg_scaled[c]) <= 1e-15,
                     "AddBoundaryFaceRhsFull: Phase 1 requires bulk_bg=0; "
                     "got |bulk_bg_scaled[" << c << "]|="
                     << std::abs(bulk_bg_scaled[c]));
      }
   }

   auto it = face_to_bdr_entry_.find(face_idx);
   MFEM_ASSERT(it != face_to_bdr_entry_.end(),
               "AddBoundaryFaceRhsFull: no boundary entry for face_idx="
               << face_idx);
   const FaceEntry &fe = entries_[it->second];
   MFEM_ASSERT(fe.nqp > 0 && fe.ndof_per_el > 0,
               "AddBoundaryFaceRhsFull: topology shape tables not built.");

   const int ndof = fe.ndof_per_el;
   real_t *rhs_data = rhs.GetData();

   for (int q = 0; q < fe.nqp; q++)
   {
      real_t I_self[NUM_STATE] = {0.0};
      const std::vector<real_t> &ss = fe.shape_self[q];
      for (int c = 0; c < NUM_STATE; c++)
      {
         real_t acc = 0.0;
         for (int i = 0; i < ndof; i++)
         {
            acc += ss[i] * Q_data[c * ndof_total + fe.dof_offset_self + i];
         }
         I_self[c] = acc;
      }

      real_t F_h[NUM_STATE] = {0.0};
      for (int c = 0; c < NUM_STATE; c++)
      {
         real_t acc = 0.0;
         for (int k = 0; k < NUM_STATE; k++)
         {
            acc += fe.nApNm1(c, k) * I_self[k];
         }
         F_h[c] = acc;
      }

      const real_t w = fe.w_qp[q];
      for (int c = 0; c < NUM_STATE; c++)
      {
         real_t *row = rhs_data + c * ndof_total + fe.dof_offset_self;
         const real_t scale = w * F_h[c];
         for (int i = 0; i < ndof; i++)
         {
            row[i] -= scale * ss[i];
         }
      }
   }
}

} // namespace seas
} // namespace mfem
