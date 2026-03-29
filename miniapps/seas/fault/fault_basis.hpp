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

/// Per-face orthonormal basis data on fault surface.
struct FaultBasisData
{
   real_t normal[3];    ///< Unit outward normal (oriented to ref_normal)
   real_t tangent1[3];  ///< First tangent (dip direction for vertical fault)
   real_t tangent2[3];  ///< Second tangent (strike direction for vertical fault)
   bool sign_flipped;   ///< True if mesh normal was flipped to align with ref_normal
                        ///< (Tandem AdapterBase convention). When true, the DG
                        ///< compute_traction uses the OPPOSITE normal, so the
                        ///< traction projection must negate the fault basis to
                        ///< compensate (double negation → correct result).
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

         // Get face element transformations
         auto *ftr = mesh.GetInteriorFaceTransformations(face_idx);
         MFEM_VERIFY(ftr != nullptr,
                     "Face " << face_idx << " is not an interior face");

         // Compute normal at face centroid
         const IntegrationPoint &ip =
            Geometries.GetCenter(ftr->GetGeometryType());
         ftr->Face->SetIntPoint(&ip);

         // Get face Jacobian and compute (unnormalized) normal
         const DenseMatrix &J = ftr->Face->Jacobian();
         Vector n_raw(dim_);
         CalcOrtho(J, n_raw);

         // Orient with reference normal (matching Tandem AdapterBase::prepare)
         real_t dot = 0.0;
         for (int d = 0; d < dim_; d++) { dot += n_raw(d) * ref_normal(d); }
         basis_[i].sign_flipped = (dot < 0.0);
         if (basis_[i].sign_flipped) { n_raw.Neg(); }

         // Normalize
         real_t n_len = n_raw.Norml2();
         MFEM_VERIFY(n_len > 0.0, "Zero-length face normal");
         n_raw /= n_len;

         // Zero-initialize all components
         for (int d = 0; d < 3; d++)
         {
            basis_[i].normal[d] = 0.0;
            basis_[i].tangent1[d] = 0.0;
            basis_[i].tangent2[d] = 0.0;
         }

         // Store normal
         for (int d = 0; d < dim_; d++)
         {
            basis_[i].normal[d] = n_raw(d);
         }

         if (dim_ == 3)
         {
            // Following Tandem: strike = normalize(up x n), dip = strike x n
            real_t s[3];
            s[0] = up(1) * n_raw(2) - up(2) * n_raw(1);
            s[1] = up(2) * n_raw(0) - up(0) * n_raw(2);
            s[2] = up(0) * n_raw(1) - up(1) * n_raw(0);

            real_t s_len = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
            MFEM_VERIFY(s_len > 1e-12,
                        "Up vector and normal are nearly collinear");
            s[0] /= s_len;
            s[1] /= s_len;
            s[2] /= s_len;

            // dip = strike x n (no negation)
            // Matches Tandem's Curvilinear.cpp:facetBasis().
            // With ref_normal=(0,-1,0), up=(0,0,1):
            //   n=(0,-1,0), strike=(1,0,0), dip=(0,0,-1) = downward.
            real_t d_vec[3];
            d_vec[0] = s[1] * n_raw(2) - s[2] * n_raw(1);
            d_vec[1] = s[2] * n_raw(0) - s[0] * n_raw(2);
            d_vec[2] = s[0] * n_raw(1) - s[1] * n_raw(0);

            // tangent1 = dip, tangent2 = strike (Tandem convention)
            for (int d = 0; d < 3; d++)
            {
               basis_[i].tangent1[d] = d_vec[d];
               basis_[i].tangent2[d] = s[d];
            }
         }
         else // dim_ == 2
         {
            // Tandem 2D convention: tangent perpendicular to normal.
            // The sign is chosen from the cross product of "up" with n so
            // that the tangent points in the "depth" direction when "up"
            // points upward in the mesh (e.g., up=(0,-1) with depth
            // positive downward gives tangent=(0,1), matching the existing
            // antiplane convention).
            real_t cross = up(0) * n_raw(1) - up(1) * n_raw(0);
            MFEM_VERIFY(std::abs(cross) > 1e-12,
                        "2D: Up vector and normal are nearly collinear");
            real_t sign_val = (cross >= 0.0) ? 1.0 : -1.0;
            basis_[i].tangent1[0] = -sign_val * n_raw(1);
            basis_[i].tangent1[1] =  sign_val * n_raw(0);
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

         const DenseMatrix &J = ftr->Face->Jacobian();
         Vector n_raw(dim_);
         CalcOrtho(J, n_raw);

         real_t dot = 0.0;
         for (int d = 0; d < dim_; d++) { dot += n_raw(d) * ref_normal(d); }
         int bi = old_count + i;
         basis_[bi].sign_flipped = (dot < 0.0);
         if (basis_[bi].sign_flipped) { n_raw.Neg(); }

         real_t n_len = n_raw.Norml2();
         MFEM_VERIFY(n_len > 0.0, "Zero-length shared face normal");
         n_raw /= n_len;

         for (int d = 0; d < 3; d++)
         {
            basis_[bi].normal[d] = 0.0;
            basis_[bi].tangent1[d] = 0.0;
            basis_[bi].tangent2[d] = 0.0;
         }

         for (int d = 0; d < dim_; d++)
         {
            basis_[bi].normal[d] = n_raw(d);
         }

         if (dim_ == 3)
         {
            real_t s[3];
            s[0] = up(1) * n_raw(2) - up(2) * n_raw(1);
            s[1] = up(2) * n_raw(0) - up(0) * n_raw(2);
            s[2] = up(0) * n_raw(1) - up(1) * n_raw(0);

            real_t s_len = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
            MFEM_VERIFY(s_len > 1e-12,
                        "Up vector and normal are nearly collinear (shared)");
            s[0] /= s_len;
            s[1] /= s_len;
            s[2] /= s_len;

            // dip = strike x n (no negation, matches Tandem)
            real_t d_vec[3];
            d_vec[0] = s[1] * n_raw(2) - s[2] * n_raw(1);
            d_vec[1] = s[2] * n_raw(0) - s[0] * n_raw(2);
            d_vec[2] = s[0] * n_raw(1) - s[1] * n_raw(0);

            for (int d = 0; d < 3; d++)
            {
               basis_[bi].tangent1[d] = d_vec[d];
               basis_[bi].tangent2[d] = s[d];
            }
         }
         else // dim_ == 2
         {
            real_t cross = up(0) * n_raw(1) - up(1) * n_raw(0);
            MFEM_VERIFY(std::abs(cross) > 1e-12,
                        "2D: Up and normal nearly collinear (shared)");
            real_t sign_val = (cross >= 0.0) ? 1.0 : -1.0;
            basis_[bi].tangent1[0] = -sign_val * n_raw(1);
            basis_[bi].tangent1[1] =  sign_val * n_raw(0);
         }
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
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_BASIS_HPP
