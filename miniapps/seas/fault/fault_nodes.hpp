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

#ifndef MFEM_SEAS_FAULT_NODES_HPP
#define MFEM_SEAS_FAULT_NODES_HPP

#include "mfem.hpp"
#include "../common/seas_types.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Manages per-face fault DOFs following Tandem's approach.
///
/// Each fault face has 2 independent DOFs (P1 basis on the face),
/// NOT shared between adjacent faces. This avoids the alternating
/// vertex mode instability that occurs with shared vertex DOFs.
///
/// For N fault faces:
///   - Total DOFs = 2*N (two per face)
///   - Face f owns DOFs {2*f, 2*f+1}
///   - DOF 2*f corresponds to face vertex at ip.x=0 (top)
///   - DOF 2*f+1 corresponds to face vertex at ip.x=1 (bottom)
///
/// Each face stores a precomputed 2x2 local mass matrix inverse
/// for L2-projected traction computation:
///   M_f = (h_f/6) * [[2,1],[1,2]]
///   M_f^{-1} = (2/h_f) * [[2,-1],[-1,2]]
///
/// Reference: Tandem AdapterBase.cpp lines 37-99
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class FaultNodes
{
public:
   using FESpaceType = FESpaceForMesh<MeshType>;

   /// @brief Construct fault nodes from mesh and face lists.
   ///
   /// @param mesh The computational mesh
   /// @param fes The DG finite element space
   /// @param order Polynomial order of the DG discretization
   /// @param interior_faces Array of MFEM interior face indices on the fault
   /// @param shared_faces Array of MFEM shared face indices on the fault
   /// @param vertex_tol Unused (kept for API compatibility)
   FaultNodes(MeshType &mesh, FESpaceType &fes, int order,
              const Array<int> &interior_faces,
              const Array<int> &shared_faces,
              real_t vertex_tol = 0.1)
      : mesh_(mesh), fes_(fes), order_(order),
        num_interior_faces_(interior_faces.Size()),
        num_shared_faces_(shared_faces.Size()),
        nodes_per_face_(order + 1)
   {
      int num_faces = num_interior_faces_ + num_shared_faces_;

      // Store face indices (interior first, then shared)
      mesh_face_indices_.SetSize(num_faces);
      is_interior_face_.SetSize(num_faces);
      for (int i = 0; i < num_interior_faces_; i++)
      {
         mesh_face_indices_[i] = interior_faces[i];
         is_interior_face_[i] = true;
      }
      for (int i = 0; i < num_shared_faces_; i++)
      {
         mesh_face_indices_[num_interior_faces_ + i] = shared_faces[i];
         is_interior_face_[num_interior_faces_ + i] = false;
      }

      // Extract vertex positions, build per-face DOF mapping, compute mass matrices
      ExtractFaceVertices();
      BuildPerFaceDOFs();
      ComputeFaceLengths();
      ComputeFaceMassMatrices();
   }

   /// @brief Number of fault DOFs (= 2 * NumFaces, per-face DOFs).
   int NumNodes() const { return num_nodes_; }

   /// @brief Total number of fault faces (interior + shared).
   int NumFaces() const { return num_interior_faces_ + num_shared_faces_; }

   /// @brief Number of interior fault faces.
   int NumInteriorFaces() const { return num_interior_faces_; }

   /// @brief Number of shared fault faces (parallel only).
   int NumSharedFaces() const { return num_shared_faces_; }

   /// @brief Number of nodes per face (= order + 1).
   int NodesPerFace() const { return nodes_per_face_; }

   /// @brief Get depths (z-coordinates) of each DOF.
   ///
   /// Returns 2*NumFaces() depths. DOFs at shared physical vertices
   /// between adjacent faces will have duplicate depth values.
   const Vector &GetDepths() const { return node_depths_; }

   /// @brief For face f (local index), return array of global DOF indices.
   ///
   /// Returns {2*f, 2*f+1} for order 1.
   /// DOF 2*f is at ip.x=0 (top vertex), DOF 2*f+1 is at ip.x=1 (bottom).
   ///
   /// @param f Local face index [0, NumFaces())
   const Array<int> &FaceToNodes(int f) const
   {
      MFEM_ASSERT(f >= 0 && f < NumFaces(),
                  "Face index out of range: " << f);
      return face_to_nodes_[f];
   }

   /// @brief Number of faces sharing a given node (always 1 for per-face DOFs).
   int NodeFaceCount(int n) const
   {
      MFEM_ASSERT(n >= 0 && n < num_nodes_,
                  "Node index out of range: " << n);
      return 1;  // Per-face DOFs are never shared
   }

   /// @brief Check if local face index f is an interior face (vs shared).
   bool IsInteriorFace(int f) const
   {
      MFEM_ASSERT(f >= 0 && f < NumFaces(),
                  "Face index out of range: " << f);
      return is_interior_face_[f];
   }

   /// @brief Get the MFEM face index for a given local fault face.
   int GetMeshFaceIndex(int f) const
   {
      MFEM_ASSERT(f >= 0 && f < NumFaces(),
                  "Face index out of range: " << f);
      return mesh_face_indices_[f];
   }

   /// @brief Get the length of a fault face.
   real_t FaceLength(int f) const
   {
      MFEM_ASSERT(f >= 0 && f < NumFaces(),
                  "Face index out of range: " << f);
      return face_lengths_(f);
   }

   /// @brief Get the precomputed 2x2 local mass matrix inverse for face f.
   ///
   /// M_f^{-1} = (2/h_f) * [[2,-1],[-1,2]] for P1 basis on face of length h_f.
   const DenseMatrix &FaceMassInverse(int f) const
   {
      MFEM_ASSERT(f >= 0 && f < NumFaces(),
                  "Face index out of range: " << f);
      return face_mass_inv_[f];
   }

   /// @brief Print fault node information for diagnostics.
   void Print(std::ostream &os = mfem::out) const
   {
      os << "FaultNodes (per-face DOFs):\n";
      os << "  Polynomial order: " << order_ << "\n";
      os << "  Nodes per face: " << nodes_per_face_ << "\n";
      os << "  Interior faces: " << num_interior_faces_ << "\n";
      os << "  Shared faces: " << num_shared_faces_ << "\n";
      os << "  Total faces: " << NumFaces() << "\n";
      os << "  Total DOFs: " << num_nodes_ << "\n";
      if (num_nodes_ > 0)
      {
         real_t z_min = node_depths_.Min();
         real_t z_max = node_depths_.Max();
         os << "  Depth range: [" << z_max / 1000.0 << ", "
            << z_min / 1000.0 << "] km\n";
      }
   }

private:
   MeshType &mesh_;
   FESpaceType &fes_;
   int order_;
   int num_interior_faces_;
   int num_shared_faces_;
   int nodes_per_face_;

   // Face data
   Array<int> mesh_face_indices_;  // MFEM face index for each local face
   Array<bool> is_interior_face_;  // true if interior, false if shared

   // Raw vertex depths (2 per face for order 1)
   std::vector<real_t> raw_vertex_depths_;

   // Per-face DOF data
   int num_nodes_;                // = 2 * NumFaces()
   Vector node_depths_;           // z-coordinate of each DOF

   // Face geometry
   Vector face_lengths_;          // length of each fault face [meters]

   // Per-face local mass matrix inverses (2x2 each)
   std::vector<DenseMatrix> face_mass_inv_;

   // Topology mappings
   std::vector<Array<int>> face_to_nodes_;  // face f → {2*f, 2*f+1}

   /// @brief Extract vertex positions from each fault face.
   void ExtractFaceVertices()
   {
      int num_faces = NumFaces();
      int dim = mesh_.Dimension();

      raw_vertex_depths_.resize(num_faces * nodes_per_face_);

      for (int f = 0; f < num_faces; f++)
      {
         FaceElementTransformations *FTr = nullptr;

         if (is_interior_face_[f])
         {
            FTr = mesh_.GetInteriorFaceTransformations(mesh_face_indices_[f]);
         }
         else
         {
            if constexpr (IsParallelMesh<MeshType>::value)
            {
#ifdef MFEM_USE_MPI
               FTr = mesh_.GetSharedFaceTransformations(mesh_face_indices_[f]);
#endif
            }
         }

         if (FTr == nullptr)
         {
            for (int j = 0; j < nodes_per_face_; j++)
            {
               raw_vertex_depths_[f * nodes_per_face_ + j] = 0.0;
            }
            continue;
         }

         // Extract endpoint positions on the face
         for (int j = 0; j < nodes_per_face_; j++)
         {
            IntegrationPoint ip;
            if (nodes_per_face_ == 1)
            {
               ip.x = 0.5;
            }
            else
            {
               ip.x = static_cast<real_t>(j) / (nodes_per_face_ - 1);
            }

            FTr->SetAllIntPoints(&ip);
            Vector coords(dim);
            FTr->Face->Transform(ip, coords);

            raw_vertex_depths_[f * nodes_per_face_ + j] = coords(1);
         }
      }
   }

   /// @brief Build per-face DOF mapping (no deduplication).
   ///
   /// Face f owns DOFs {2*f, 2*f+1}. Each DOF gets the z-coordinate
   /// of its corresponding face vertex.
   void BuildPerFaceDOFs()
   {
      int num_faces = NumFaces();
      num_nodes_ = nodes_per_face_ * num_faces;  // 2*N for order 1

      // Assign depths from raw vertex data
      node_depths_.SetSize(num_nodes_);
      for (int f = 0; f < num_faces; f++)
      {
         for (int j = 0; j < nodes_per_face_; j++)
         {
            int dof = f * nodes_per_face_ + j;
            node_depths_(dof) = raw_vertex_depths_[dof];
         }
      }

      // Build trivial face-to-node mapping: face f → {2*f, 2*f+1}
      face_to_nodes_.resize(num_faces);
      for (int f = 0; f < num_faces; f++)
      {
         face_to_nodes_[f].SetSize(nodes_per_face_);
         for (int j = 0; j < nodes_per_face_; j++)
         {
            face_to_nodes_[f][j] = f * nodes_per_face_ + j;
         }
      }
   }

   /// @brief Compute face lengths from vertex depths.
   void ComputeFaceLengths()
   {
      int nf = NumFaces();
      face_lengths_.SetSize(nf);
      for (int f = 0; f < nf; f++)
      {
         real_t z0 = node_depths_(face_to_nodes_[f][0]);
         real_t z1 = node_depths_(face_to_nodes_[f][1]);
         face_lengths_(f) = std::abs(z1 - z0);
      }
   }

   /// @brief Precompute per-face local mass matrix inverses.
   ///
   /// For P1 basis on face of length h:
   ///   M = (h/6) * [[2,1],[1,2]]
   ///   M^{-1} = (2/h) * [[2,-1],[-1,2]]
   void ComputeFaceMassMatrices()
   {
      int nf = NumFaces();
      face_mass_inv_.resize(nf);

      for (int f = 0; f < nf; f++)
      {
         real_t h = face_lengths_(f);
         MFEM_ASSERT(h > 0, "Zero-length face " << f);

         real_t coeff = 2.0 / h;
         face_mass_inv_[f].SetSize(nodes_per_face_);

         // M^{-1} = (2/h) * [[2,-1],[-1,2]]
         face_mass_inv_[f](0, 0) =  2.0 * coeff;
         face_mass_inv_[f](0, 1) = -1.0 * coeff;
         face_mass_inv_[f](1, 0) = -1.0 * coeff;
         face_mass_inv_[f](1, 1) =  2.0 * coeff;
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_NODES_HPP
