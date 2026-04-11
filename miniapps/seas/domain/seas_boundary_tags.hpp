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

/// @brief Standard boundary attribute tags for SEAS meshes (Tandem convention).
///
/// These tags match Tandem's Gmsh Physical Surface IDs:
///   - Tag 1 = Natural (top Z=0 + bottom Z=Z0), zero traction
///   - Tag 3 = Fault (Y=0 interior interface)
///   - Tag 5 = Dirichlet (far-field vertical faces), plate loading
struct SEASBoundaryTags
{
   static constexpr int NATURAL   = 1;  ///< Top (Z=0) + bottom (Z=Z0), zero traction
   static constexpr int FAULT     = 3;  ///< Fault interior (Y=0)
   static constexpr int DIRICHLET = 5;  ///< Far-field vertical faces, plate loading
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
   /// centroid has |Y| < coord_tol (Tandem convention: fault at Y=0).
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

         // Fault is at Y = 0 (second coordinate, Tandem convention)
         if (std::abs(coords(1)) < coord_tol)
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
