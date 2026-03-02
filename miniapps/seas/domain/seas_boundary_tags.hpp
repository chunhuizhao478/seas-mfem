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

#ifndef MFEM_SEAS_BOUNDARY_TAGS_HPP
#define MFEM_SEAS_BOUNDARY_TAGS_HPP

#include "mfem.hpp"

namespace mfem
{
namespace seas
{

/// @brief Standard boundary attribute tags for SEAS meshes.
///
/// These tags correspond to Gmsh Physical Curve (2D) or Physical Surface (3D)
/// IDs used in .geo files. They follow Tandem's convention:
///   - Tags 1-4: outer boundaries
///   - Tag 5: fault interior interface
///
/// Existing BP2BoundaryAttributes (tags 1-4) remains for backward compatibility.
struct SEASBoundaryTags
{
   static constexpr int FARFIELD_LEFT  = 1;  ///< x = -Lx boundary
   static constexpr int FARFIELD_RIGHT = 2;  ///< x = +Lx boundary
   static constexpr int FREE_SURFACE   = 3;  ///< z = 0 (top) boundary
   static constexpr int BOTTOM         = 4;  ///< z = -Lz (bottom) boundary
   static constexpr int FAULT          = 5;  ///< Fault interior interface
};

/// @brief Utility for finding fault faces from Gmsh Physical Group tags.
///
/// When a Gmsh .geo file tags the fault interior edges/surfaces with a
/// Physical Curve/Surface, MFEM loads them as boundary elements with
/// the corresponding attribute. This class scans those boundary elements
/// to find the interior face indices that constitute the fault.
///
/// Usage:
/// @code
///   Array<int> fault_faces = FaultBoundaryData::FindFaultInteriorFaces(
///       mesh, SEASBoundaryTags::FAULT);
/// @endcode
class FaultBoundaryData
{
public:
   /// @brief Find interior faces tagged as fault in the mesh.
   ///
   /// Scans all boundary elements for the given fault_tag attribute.
   /// For each matching boundary element, checks whether the corresponding
   /// face is an interior face (has two adjacent volume elements).
   /// Returns the array of interior face indices that lie on the fault.
   ///
   /// @param mesh The mesh to scan (non-const due to MFEM's GetInteriorFaceTransformations)
   /// @param fault_tag The Gmsh Physical Curve/Surface tag for the fault
   /// @return Array of interior face indices on the fault
   static Array<int> FindFaultInteriorFaces(Mesh &mesh, int fault_tag)
   {
      Array<int> fault_faces;

      // Check if the mesh has any boundary attributes matching the tag
      if (mesh.bdr_attributes.Size() == 0 ||
          fault_tag > mesh.bdr_attributes.Max())
      {
         return fault_faces;
      }

      for (int be = 0; be < mesh.GetNBE(); be++)
      {
         if (mesh.GetBdrAttribute(be) != fault_tag) { continue; }

         // Get the face index for this boundary element
         int face, ori;
         mesh.GetBdrElementFace(be, &face, &ori);

         // Check if this face is an interior face (two adjacent elements)
         // In MFEM, GetInteriorFaceTransformations returns nullptr for
         // true boundary faces, and a valid pointer for interior faces.
         if (mesh.GetInteriorFaceTransformations(face) != nullptr)
         {
            fault_faces.Append(face);
         }
      }

      return fault_faces;
   }

#ifdef MFEM_USE_MPI
   /// @brief Find shared faces on the fault for parallel meshes.
   ///
   /// In parallel, fault faces at MPI partition boundaries become "shared
   /// faces" rather than interior faces. MFEM does not directly expose
   /// boundary attributes on shared faces, so this method uses a
   /// coordinate-based fallback: a shared face is on the fault if its
   /// centroid has |x| < coord_tol.
   ///
   /// This is acceptable for axis-aligned faults (BP1/BP2/BP5 at x=0).
   /// For non-axis-aligned faults, a more general approach (e.g.,
   /// communicating boundary element tags across partitions) would be
   /// needed.
   ///
   /// @param pmesh The parallel mesh
   /// @param fault_tag The fault tag (unused for shared faces, reserved
   ///                  for future tag-based shared face detection)
   /// @param coord_tol Tolerance for coordinate-based shared face detection
   /// @return Array of shared face indices on the fault
   static Array<int> FindFaultSharedFaces(
      ParMesh &pmesh, int fault_tag, real_t coord_tol = 1e-6)
   {
      // Suppress unused parameter warning — fault_tag reserved for future use
      (void)fault_tag;

      Array<int> fault_shared_faces;

      for (int sf = 0; sf < pmesh.GetNSharedFaces(); sf++)
      {
         auto *FTr = pmesh.GetSharedFaceTransformations(sf);
         if (FTr == nullptr) { continue; }

         // Evaluate face centroid
         IntegrationPoint ip;
         ip.x = 0.5;
         if (pmesh.Dimension() == 3) { ip.y = 0.5; }

         FTr->SetAllIntPoints(&ip);
         Vector coords(pmesh.Dimension());
         FTr->Face->Transform(ip, coords);

         // Fault is at x = 0 (first coordinate)
         if (std::abs(coords(0)) < coord_tol)
         {
            fault_shared_faces.Append(sf);
         }
      }

      return fault_shared_faces;
   }
#endif
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BOUNDARY_TAGS_HPP
