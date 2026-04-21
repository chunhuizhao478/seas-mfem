// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-443211.
// All rights reserved. See file COPYRIGHT for details.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_FAULT_BASIS_HPP
#define MFEM_SEAS_FAULT_BASIS_HPP

#include "mfem.hpp"
#include <cmath>
#include <vector>

namespace mfem
{
namespace seas
{

/// Per-quad-point data on a fault face (matching Tandem AdapterBase).
/// Sign convention: all vectors already include the sign_flipped negation
/// (Tandem convention).  Callers use them directly, no sign factor needed.
struct FaultBasisQPData
{
   real_t normal[3];       ///< Unit normal (with sign baked in)
   real_t tangent1[3];     ///< First tangent / dip (with sign baked in)
   real_t tangent2[3];     ///< Second tangent / strike (with sign baked in)
   real_t nl;              ///< Normal length |n_raw| at this quad point
   bool sign_flipped;      ///< Diagnostic only; sign already baked into basis vectors
};

/// Per-face summary (for backward compatibility and fast access).
/// Sign convention: all vectors already include the sign_flipped negation
/// (Tandem convention).  Callers use them directly, no sign factor needed.
struct FaultBasisData
{
   real_t normal[3];    ///< Unit normal at centroid (with sign baked in)
   real_t tangent1[3];  ///< First tangent / dip at centroid (with sign baked in)
   real_t tangent2[3];  ///< Second tangent / strike at centroid (with sign baked in)
   bool sign_flipped;   ///< Diagnostic only; sign already baked into basis vectors
   std::vector<FaultBasisQPData> qp_data;  ///< Per-quad-point data (Tandem convention)
};

/// Computes and stores per-face local coordinate frames on a fault surface.
///
/// Algorithm (following Tandem Curvilinear::facetBasis):
///   1. Get raw face normal from mesh geometry (CalcOrtho)
///   2. Orient: if dot(n_raw, ref_normal) < 0, flip
///   3. Normalize: n = n_raw / |n_raw|
///   4. strike = normalize(up x n)
///   5. dip = strike x n  (unit since strike _|_ n)
///   6. tangent1 = dip, tangent2 = strike  (Tandem convention)
///
/// For 2D (dim=2): Single tangent computed as 90-degree rotation of normal,
/// with sign chosen so that "up" direction is consistent.
///
/// Reference: Tandem src/geometry/Curvilinear.cpp:facetBasis()
class FaultBasis
{
public:
   /// Compute basis for all fault faces.
   ///
   /// One basis per fault face (constant normal per planar face).
   /// For DG methods, each face has independent DOFs sharing one basis.
   ///
   /// @param mesh        The mesh (serial)
   /// @param fault_faces Interior face indices forming the fault
   /// @param ref_normal  Reference normal for consistent orientation
   /// @param up          Reference "up" vector for strike/dip decomposition
   void Compute(Mesh &mesh,
                const Array<int> &fault_faces,
                const Vector &ref_normal,
                const Vector &up)
   {
      dim_ = mesh.Dimension();
      num_faces_ = fault_faces.Size();
      basis_.resize(num_faces_);

      for (int i = 0; i < num_faces_; i++)
      {
         int face_idx = fault_faces[i];

         auto *ftr = mesh.GetInteriorFaceTransformations(face_idx);
         MFEM_VERIFY(ftr != nullptr,
                     "Face " << face_idx << " is not an interior face");

         const IntegrationPoint &ip =
            Geometries.GetCenter(ftr->GetGeometryType());
         ftr->Face->SetIntPoint(&ip);

         Vector n_raw(dim_);
         CalcOrtho(ftr->Face->Jacobian(), n_raw);

         real_t nl;
         ComputeOrientedFrame(n_raw, dim_, ref_normal, up,
                              basis_[i].normal, basis_[i].tangent1,
                              basis_[i].tangent2, basis_[i].sign_flipped, nl);
         MFEM_VERIFY(nl > 0.0, "Zero-length face normal");
      }
   }

   /// Populate per-quad-point basis data for all faces (Tandem AdapterBase style).
   ///
   /// Must be called AFTER Compute(). Evaluates normal, tangent1, tangent2,
   /// nl, sign_flipped at each quad point on each face.
   ///
   /// @param mesh The mesh
   /// @param fault_faces Interior face indices
   /// @param ref_normal Reference normal
   /// @param up Up vector
   /// @param ir Quadrature rule for face integration
   void ComputeQPBasis(Mesh &mesh,
                        const Array<int> &fault_faces,
                        const Vector &ref_normal,
                        const Vector &up,
                        const IntegrationRule &ir)
   {
      int nq = ir.GetNPoints();
      for (int i = 0; i < fault_faces.Size(); i++)
      {
         int face_idx = fault_faces[i];
         auto *ftr = mesh.GetInteriorFaceTransformations(face_idx);
         if (!ftr) { continue; }

         basis_[i].qp_data.resize(nq);
         for (int q = 0; q < nq; q++)
         {
            ftr->Face->SetIntPoint(&ir.IntPoint(q));
            Vector n_raw(dim_);
            CalcOrtho(ftr->Face->Jacobian(), n_raw);

            auto &qpd = basis_[i].qp_data[q];
            ComputeOrientedFrame(n_raw, dim_, ref_normal, up,
                                 qpd.normal, qpd.tangent1, qpd.tangent2,
                                 qpd.sign_flipped, qpd.nl);
         }
      }
   }

   /// Same as ComputeQPBasis but for shared faces (parallel).
   template <typename PMeshType>
   void ComputeQPBasisShared(PMeshType &mesh,
                              const Array<int> &shared_faces,
                              const Vector &ref_normal,
                              const Vector &up,
                              const IntegrationRule &ir,
                              int interior_face_count)
   {
      int nq = ir.GetNPoints();
      for (int i = 0; i < shared_faces.Size(); i++)
      {
         int sf = shared_faces[i];
         auto *ftr = mesh.GetSharedFaceTransformations(sf);
         if (!ftr) { continue; }

         int bi = interior_face_count + i;
         basis_[bi].qp_data.resize(nq);
         for (int q = 0; q < nq; q++)
         {
            ftr->Face->SetIntPoint(&ir.IntPoint(q));
            Vector n_raw(dim_);
            CalcOrtho(ftr->Face->Jacobian(), n_raw);

            auto &qpd = basis_[bi].qp_data[q];
            ComputeOrientedFrame(n_raw, dim_, ref_normal, up,
                                 qpd.normal, qpd.tangent1, qpd.tangent2,
                                 qpd.sign_flipped, qpd.nl);
         }
      }
   }

   /// Append basis data for shared fault faces (parallel).
   ///
   /// Uses GetSharedFaceTransformations instead of GetInteriorFaceTransformations.
   /// Appends to existing basis_ array computed by Compute().
   template <typename MeshType>
   void AppendSharedFaces(MeshType &mesh,
                          const Array<int> &shared_faces,
                          const Vector &ref_normal,
                          const Vector &up)
   {
#ifdef MFEM_USE_MPI
      int ns = shared_faces.Size();
      if (ns == 0) { return; }

      int old_count = num_faces_;
      num_faces_ += ns;
      basis_.resize(num_faces_);

      for (int i = 0; i < ns; i++)
      {
         int sf = shared_faces[i];
         auto *ftr = mesh.GetSharedFaceTransformations(sf);
         MFEM_VERIFY(ftr != nullptr,
                     "Shared face " << sf << " returned null transformation");

         const IntegrationPoint &ip =
            Geometries.GetCenter(ftr->GetGeometryType());
         ftr->Face->SetIntPoint(&ip);

         Vector n_raw(dim_);
         CalcOrtho(ftr->Face->Jacobian(), n_raw);

         int bi = old_count + i;
         real_t nl;
         ComputeOrientedFrame(n_raw, dim_, ref_normal, up,
                              basis_[bi].normal, basis_[bi].tangent1,
                              basis_[bi].tangent2, basis_[bi].sign_flipped, nl);
         MFEM_VERIFY(nl > 0.0, "Zero-length shared face normal");
      }
#endif
   }

   /// Number of fault faces with computed basis.
   int NumFaces() const { return num_faces_; }

   /// Dimension of the mesh (2 or 3).
   int Dimension() const { return dim_; }

   /// Number of tangential components (1 for 2D, 2 for 3D).
   int NumTangentComponents() const { return dim_ - 1; }

   /// Get basis data for face fi.
   const FaultBasisData &GetBasis(int fi) const
   {
      MFEM_ASSERT(fi >= 0 && fi < num_faces_,
                  "FaultBasis::GetBasis: face index " << fi
                  << " out of range [0, " << num_faces_ << ")");
      return basis_[fi];
   }

   /// Project global traction vector to fault-local tangential components.
   ///
   /// For 3D: tau_local[0] = traction . tangent1 (dip)
   ///         tau_local[1] = traction . tangent2 (strike)
   /// For 2D: tau_local[0] = traction . tangent1
   ///
   /// NOTE: tangent vectors already contain the sign_flipped negation
   /// (Tandem convention).  No additional sign correction is needed.
   void ProjectTraction(int fi, const real_t *traction_global,
                        real_t *tau_local) const
   {
      MFEM_ASSERT(fi >= 0 && fi < num_faces_,
                  "FaultBasis::ProjectTraction: face index " << fi
                  << " out of range [0, " << num_faces_ << ")");
      const auto &b = basis_[fi];
      tau_local[0] = 0.0;
      for (int d = 0; d < dim_; d++)
      {
         tau_local[0] += traction_global[d] * b.tangent1[d];
      }
      if (dim_ == 3)
      {
         tau_local[1] = 0.0;
         for (int d = 0; d < dim_; d++)
         {
            tau_local[1] += traction_global[d] * b.tangent2[d];
         }
      }
   }

   /// Embed fault-local slip into global displacement jump vector.
   ///
   /// For 3D: delta_u = slip[0]*tangent1 + slip[1]*tangent2
   /// For 2D: delta_u = slip[0]*tangent1
   void EmbedSlip(int fi, const real_t *slip_local,
                  real_t *delta_u_global) const
   {
      MFEM_ASSERT(fi >= 0 && fi < num_faces_,
                  "FaultBasis::EmbedSlip: face index " << fi
                  << " out of range [0, " << num_faces_ << ")");
      const auto &b = basis_[fi];
      for (int d = 0; d < dim_; d++)
      {
         delta_u_global[d] = slip_local[0] * b.tangent1[d];
      }
      if (dim_ == 3)
      {
         for (int d = 0; d < dim_; d++)
         {
            delta_u_global[d] += slip_local[1] * b.tangent2[d];
         }
      }
   }

   /// Per-quad-point embed: uses qp_data tangent vectors instead of centroid.
   ///
   /// For 3D: delta_u = slip[0]*tangent1(q) + slip[1]*tangent2(q)
   void EmbedSlipQP(int fi, int q, const real_t *slip_local,
                    real_t *delta_u_global) const
   {
      MFEM_ASSERT(fi >= 0 && fi < num_faces_,
                  "FaultBasis::EmbedSlipQP: face index out of range");
      const auto &b = basis_[fi];
      if (q >= 0 && q < static_cast<int>(b.qp_data.size()))
      {
         const auto &qd = b.qp_data[q];
         for (int d = 0; d < dim_; d++)
         {
            delta_u_global[d] = slip_local[0] * qd.tangent1[d];
         }
         if (dim_ == 3)
         {
            for (int d = 0; d < dim_; d++)
            {
               delta_u_global[d] += slip_local[1] * qd.tangent2[d];
            }
         }
      }
      else
      {
         // Fallback to centroid basis if qp_data not populated
         EmbedSlip(fi, slip_local, delta_u_global);
      }
   }

   /// Get the normal stress component: sigma_n = -traction . normal.
   /// Positive in compression.
   real_t NormalStress(int fi, const real_t *traction_global) const
   {
      MFEM_ASSERT(fi >= 0 && fi < num_faces_,
                  "FaultBasis::NormalStress: face index " << fi
                  << " out of range [0, " << num_faces_ << ")");
      const auto &b = basis_[fi];
      real_t dot = 0.0;
      for (int d = 0; d < dim_; d++)
      {
         dot += traction_global[d] * b.normal[d];
      }
      return -dot;
   }

private:
   int dim_ = 0;
   int num_faces_ = 0;
   std::vector<FaultBasisData> basis_;

public:
   /// Compute oriented fault frame (normal, dip, strike) from a raw face
   /// normal, replicating Tandem's AdapterBase::prepare() convention exactly.
   ///
   /// Exposed as public (v9.1.0) so `test_fault_basis_dip_strike_symmetry.cpp`
   /// can exercise the R-003 orthonormality contract on 1000 perturbed-
   /// normal samples without constructing a mesh per sample.  The method is
   /// pure (static, no class-state dependency); public promotion has no
   /// semantic consequence.
   ///
   /// Tandem's algorithm (AdapterBase.cpp:62-84, Curvilinear.cpp:260-293):
   ///   1. Record nl = |n_raw|
   ///   2. sign_flipped = dot(n_raw, ref_normal) < 0
   ///   3. If sign_flipped: flip n_raw to ref-aligned  (line 72-73)
   ///   4. Compute tangent frame from ref-aligned normal:
   ///        strike = normalize(up × n_ref)   (facetBasis:287)
   ///        dip    = normalize(strike × n_ref) (facetBasis:291)
   ///   5. If sign_flipped: negate ALL stored vectors  (lines 78-84)
   ///
   /// The stored basis already contains the sign.  Callers use it directly
   /// — NO separate sign factor is needed.
   ///
   /// @param n_raw       Raw face normal (modified in-place)
   /// @param dim         Spatial dimension (2 or 3)
   /// @param ref_normal  Reference normal for orientation
   /// @param up          Reference up vector
   /// @param[out] normal      Unit normal with sign (3 components)
   /// @param[out] tangent1    Dip tangent with sign (3 components)
   /// @param[out] tangent2    Strike tangent with sign (3 components)
   /// @param[out] sign_flipped  True if raw normal was anti-aligned with ref
   /// @param[out] nl          Raw normal length (before normalization)
   static void ComputeOrientedFrame(Vector &n_raw, int dim,
                                     const Vector &ref_normal,
                                     const Vector &up,
                                     real_t normal[3],
                                     real_t tangent1[3],
                                     real_t tangent2[3],
                                     bool &sign_flipped,
                                     real_t &nl)
   {
      nl = n_raw.Norml2();

      // Step 1-2: Determine orientation relative to reference normal
      real_t dot = 0.0;
      for (int d = 0; d < dim; d++) { dot += n_raw(d) * ref_normal(d); }
      sign_flipped = (dot < 0.0);

      // Step 3: Flip to ref-aligned (Tandem AdapterBase.cpp:72-73)
      if (sign_flipped) { n_raw.Neg(); }

      // Normalize
      if (nl > 0.0) { n_raw /= nl; }

      // Zero-initialize output
      for (int d = 0; d < 3; d++)
      {
         normal[d] = 0.0;
         tangent1[d] = 0.0;
         tangent2[d] = 0.0;
      }
      for (int d = 0; d < dim; d++) { normal[d] = n_raw(d); }

      // Step 4: Compute tangent frame from ref-aligned normal
      if (dim == 3)
      {
         // strike = normalize(up × n_ref) (Tandem facetBasis:287)
         real_t s[3];
         s[0] = up(1) * n_raw(2) - up(2) * n_raw(1);
         s[1] = up(2) * n_raw(0) - up(0) * n_raw(2);
         s[2] = up(0) * n_raw(1) - up(1) * n_raw(0);
         real_t s_len = std::sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]);
         MFEM_VERIFY(s_len > 1e-12,
                     "ComputeOrientedFrame: up vector is nearly collinear with "
                     "face normal (|up x n| = " << s_len << "). Cannot form "
                     "strike/dip tangent frame.");
         s[0] /= s_len; s[1] /= s_len; s[2] /= s_len;

         // dip = strike × n_ref (Tandem facetBasis:291)
         real_t dv[3];
         dv[0] = s[1] * n_raw(2) - s[2] * n_raw(1);
         dv[1] = s[2] * n_raw(0) - s[0] * n_raw(2);
         dv[2] = s[0] * n_raw(1) - s[1] * n_raw(0);

         // R-003 (v9.1.0, plan §3.4 H-V91-B1): normalize dip explicitly.
         // |strike|=1 and strike⟂n guarantee |dip|=1 in exact arithmetic,
         // but FP rounding in the strike normalize above leaks into this
         // cross product so |dip| drifts ~5 ULP from 1.  Downstream
         // BuildRotation / BuildRotationInverse and the per-side flux
         // consumers assume orthonormality; defensive normalize.
         const real_t d_len = std::sqrt(dv[0]*dv[0] + dv[1]*dv[1]
                                        + dv[2]*dv[2]);
         MFEM_VERIFY(d_len > 1e-12,
                     "ComputeOrientedFrame: degenerate dip vector |dv|="
                     << d_len << ".  strike/ref_normal nearly parallel?");
         const real_t d_inv = 1.0 / d_len;
         dv[0] *= d_inv; dv[1] *= d_inv; dv[2] *= d_inv;

         for (int d = 0; d < 3; d++) { tangent1[d] = dv[d]; tangent2[d] = s[d]; }
      }
      else // dim == 2
      {
         real_t cross = up(0) * n_raw(1) - up(1) * n_raw(0);
         real_t sv = (cross >= 0.0) ? 1.0 : -1.0;
         tangent1[0] = -sv * n_raw(1);
         tangent1[1] =  sv * n_raw(0);
      }

      // Step 5: Negate ALL if sign_flipped (Tandem AdapterBase.cpp:78-84)
      // The stored basis now contains the sign — callers use it directly.
      //
      // NOTE (v61 analysis): For shared faces in parallel, both ranks
      // compute basis independently. CalcOrtho gives opposite normals on
      // each rank → one rank has sign_flipped=true, the other false →
      // basis vectors differ by a global sign: tangent_A = -tangent_B.
      //
      // This is CORRECT because the DG face normal also flips between
      // ranks (nor_B = -nor_A, proven by K matching to 15 digits).
      // The embedded displacement jump flips accordingly (f_q_B = -f_q_A),
      // and the assembly produces elvec1_B = elvec2_A (verified by
      // test_serial_parallel_displacement_match to rel_err < 1e-12).
      // See debug document v61 Section 8.9 for the full proof.
      if (sign_flipped)
      {
         for (int d = 0; d < 3; d++)
         {
            normal[d] = -normal[d];
            tangent1[d] = -tangent1[d];
            tangent2[d] = -tangent2[d];
         }
      }
   }

private:
   // R-008 (v9.1.0 rev 3): explicit closing `private:` section.  The
   // v9.1.0 R-003 patch opened a mid-class `public:` block above the
   // static `ComputeOrientedFrame` (so the test can call it directly).
   // Without this trailing section, any member added after the method
   // would silently inherit `public:` visibility.  Intentionally empty;
   // add future private members here.
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_BASIS_HPP
