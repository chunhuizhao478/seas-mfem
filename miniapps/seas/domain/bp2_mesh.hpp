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
// !!! what's the relation between this file and antiplane_operator.hpp

#ifndef MFEM_SEAS_BP2_MESH_HPP
#define MFEM_SEAS_BP2_MESH_HPP

#include "mfem.hpp"
#include "../config/bp2_params.hpp"

#include <memory>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

namespace mfem
{
namespace seas
{

/// @brief Boundary attribute indices for the FULL DOMAIN BP2 mesh
///
/// Full domain: x ∈ [-Lx, +Lx], z ∈ [-Lz, 0]
/// The fault at x = 0 is an INTERIOR interface, not a boundary.
///
/// Z-coordinate convention (geophysics standard):
/// - z = 0: free surface (Earth's surface)
/// - z < 0: depth (negative values indicate greater depth)
/// - z = -Lz: bottom of computational domain
///
/// All outer boundaries have Natural BC (zero traction).
/// Tectonic loading is applied directly on the fault at z < -Wf
/// by prescribing slip rate = Vp in the fault operator.
struct BP2BoundaryAttributes
{
   static constexpr int FARFIELD_LEFT = 1;   ///< Left far-field (x = -Lx) - Natural BC
   static constexpr int FARFIELD_RIGHT = 2;  ///< Right far-field (x = +Lx) - Natural BC
   static constexpr int FREE_SURFACE = 3;    ///< Free surface (z = 0) - Natural BC
   static constexpr int BOTTOM = 4;          ///< Bottom (z = -Lz) - Natural BC

   // Note: The fault is NOT a boundary in full domain - it's interior faces at x = 0
};

/// @brief BP2 mesh generator for the FULL DOMAIN 2D antiplane shear problem
///
/// Creates a 2D rectangular mesh for the full domain:
///   x ∈ [-Lx, +Lx], z ∈ [-Lz, 0]
///
/// Z-coordinate convention:
///   z = 0: free surface (Earth's surface)
///   z < 0: depth (negative values going into the Earth)
///   z = -Lz: bottom boundary
///
/// The fault at x = 0 is an INTERIOR interface where slip is imposed
/// as a jump condition: [[u]] = u(0⁺) - u(0⁻) = slip
///
/// Boundary attributes (for OUTER boundaries only):
///   1 = Far-field left (x = -Lx) - Natural BC: zero traction
///   2 = Far-field right (x = +Lx) - Natural BC: zero traction
///   3 = Free surface (z = 0) - Natural BC: zero traction
///   4 = Bottom (z = -Lz) - Natural BC: zero traction
///
/// Tectonic loading: prescribed V=Vp on fault faces below Wf (z < -Wf)
/// Both sides of the fault are explicitly modeled
class BP2MeshGenerator
{
public:
   /// @brief Parameters for mesh generation
   struct Parameters
   {
      real_t Lx = 100.0e3;   ///< Domain half-width [m] (100 km), full domain is [-Lx, +Lx]
      real_t Lz = 100.0e3;   ///< Domain depth [m] (100 km)
      real_t Wf = 40.0e3;    ///< Rate-state fault depth [m] (40 km)

      int nx = 100;          ///< Number of elements in x-direction (for half domain, 2*nx total)
      int nz = 100;          ///< Number of elements in z-direction

      /// @brief Mesh grading parameters for refinement near fault
      real_t grading_x = 1.0; ///< Grading in x (1.0 = uniform)
      real_t grading_z = 1.0; ///< Grading in z (1.0 = uniform)
   };

   /// @brief Create a BP2 full domain mesh with default parameters
   static std::unique_ptr<Mesh> Create()
   {
      return Create(Parameters());
   }

   /// @brief Create a BP2 full domain mesh with specified parameters
   ///
   /// Creates mesh for x ∈ [-Lx, +Lx], z ∈ [-Lz, 0]
   /// Ensures that x = 0 aligns with element boundaries (for clean fault interface)
   static std::unique_ptr<Mesh> Create(const Parameters &params)
   {
      // Create mesh for [0, 2*Lx] × [0, Lz] initially
      // Use 2*nx elements in x to ensure symmetric mesh about x=Lx (which becomes x=0)
      auto mesh = std::make_unique<Mesh>(
         Mesh::MakeCartesian2D(2 * params.nx, params.nz,
                               Element::QUADRILATERAL,
                               true,  // generate_edges
                               2.0 * params.Lx, params.Lz));

      // Shift coordinates:
      // x: [0, 2*Lx] → [-Lx, +Lx]
      // z: [0, Lz] → [-Lz, 0]
      ShiftMeshCoordinates(*mesh, params.Lx, params.Lz);

      // Set boundary attributes for the outer boundaries
      SetBoundaryAttributes(*mesh, params);

      return mesh;
   }

   /// @brief Create a refined BP2 mesh with grading near fault
   static std::unique_ptr<Mesh> CreateGraded(const Parameters &params)
   {
      auto mesh = Create(params);

      if (params.grading_x != 1.0 || params.grading_z != 1.0)
      {
         ApplyGrading(*mesh, params);
      }

      return mesh;
   }

   /// @brief Create a simple test mesh for unit testing
   ///
   /// Creates a small 8x4 mesh (4 on each side of fault) for quick tests
   /// Domain: x ∈ [-10km, +10km], z ∈ [-10km, 0]
   static std::unique_ptr<Mesh> CreateTestMesh()
   {
      Parameters params;
      params.Lx = 10.0e3;  // 10 km half-width, full domain [-10km, +10km]
      params.Lz = 10.0e3;  // 10 km depth
      params.Wf = 10.0e3;  // Entire depth is fault
      params.nx = 4;       // 4 elements per side, 8 total in x
      params.nz = 4;
      return Create(params);
   }

   /// @brief Create a mesh for MMS (Method of Manufactured Solutions) testing
   ///
   /// Creates a mesh suitable for convergence studies
   /// Domain: x ∈ [-1, +1], z ∈ [-1, 0]
   static std::unique_ptr<Mesh> CreateMMSMesh(int refinement_level)
   {
      Parameters params;
      params.Lx = 1.0;  // Unit half-width, full domain [-1, +1]
      params.Lz = 1.0;
      params.Wf = 1.0;

      // Increase resolution with refinement level
      params.nx = 4 * (1 << refinement_level);  // 4, 8, 16, 32, ... per side
      params.nz = params.nx;

      return Create(params);
   }

   /// @brief Load a Gmsh .msh mesh file and scale coordinates.
   ///
   /// The .geo file uses km; this method scales to meters (scale=1000).
   /// Physical Curve IDs in the .geo file must match BP2BoundaryAttributes.
   ///
   /// @param filename Path to .msh file (Gmsh format 2.2)
   /// @param scale    Coordinate scale factor (default 1000 = km→m)
   static std::unique_ptr<Mesh> LoadGmshMesh(const std::string &filename,
                                              real_t scale = 1000.0)
   {
      auto mesh = std::make_unique<Mesh>(filename.c_str(), 1, 1, true);

      // Scale coordinates (e.g., km → m)
      if (std::abs(scale - 1.0) > 1e-14)
      {
         for (int i = 0; i < mesh->GetNV(); i++)
         {
            real_t *v = mesh->GetVertex(i);
            for (int d = 0; d < mesh->SpaceDimension(); d++)
            {
               v[d] *= scale;
            }
         }
      }

      return mesh;
   }

   /// @brief Save mesh to VTK format for ParaView visualization
   static void SaveVTK(Mesh &mesh, const std::string &filename)
   {
      std::ofstream vtk_file(filename);
      if (!vtk_file.is_open())
      {
         MFEM_WARNING("Could not open VTK file: " << filename);
         return;
      }
      mesh.PrintVTK(vtk_file);
      vtk_file.close();
   }

   /// @brief Get the x-coordinate of a point and determine if it's on the fault
   ///
   /// @param x x-coordinate
   /// @param z z-coordinate (z=0 at surface, z<0 at depth)
   /// @param Wf fault depth (positive value, fault extends from z=0 to z=-Wf)
   /// @param tol tolerance for x=0 check
   /// @return true if point is on the fault (x≈0, z > -Wf)
   static bool IsOnFault(real_t x, real_t z, real_t Wf, real_t tol = 1e-10)
   {
      // Fault is at x=0, extending from z=0 (surface) to z=-Wf (depth)
      return std::abs(x) < tol && z > -Wf - tol;
   }

private:
   /// @brief Shift mesh coordinates to center domain
   ///
   /// Shifts x from [0, 2*Lx] to [-Lx, +Lx]
   /// Shifts z from [0, Lz] to [-Lz, 0]
   static void ShiftMeshCoordinates(Mesh &mesh, real_t Lx, real_t Lz)
   {
      for (int i = 0; i < mesh.GetNV(); i++)
      {
         real_t *v = mesh.GetVertex(i);
         v[0] -= Lx;  // Shift from [0, 2*Lx] to [-Lx, +Lx]
         v[1] -= Lz;  // Shift from [0, Lz] to [-Lz, 0]
      }
   }

   /// @brief Set boundary attributes for the full domain mesh
   ///
   /// Only the 4 OUTER boundaries get attributes.
   /// Interior faces at x=0 (the fault) are NOT boundaries.
   ///
   /// After coordinate shift, domain is x ∈ [-Lx, +Lx], z ∈ [-Lz, 0]
   static void SetBoundaryAttributes(Mesh &mesh, const Parameters &params)
   {
      for (int be = 0; be < mesh.GetNBE(); be++)
      {
         // Get the center of this boundary element
         Array<int> vertices;
         mesh.GetBdrElementVertices(be, vertices);

         // Compute centroid of boundary element
         real_t x_center = 0.0;
         real_t z_center = 0.0;
         for (int v = 0; v < vertices.Size(); v++)
         {
            const real_t *coords = mesh.GetVertex(vertices[v]);
            x_center += coords[0];
            z_center += coords[1];
         }
         x_center /= vertices.Size();
         z_center /= vertices.Size();

         // Determine boundary attribute based on position
         const real_t tol = 1e-10 * std::max(params.Lx, params.Lz);

         int attr = 0;

         // Check outer boundaries (full domain: x ∈ [-Lx, +Lx], z ∈ [-Lz, 0])
         if (std::abs(x_center + params.Lx) < tol)
         {
            // Left boundary (x = -Lx): Natural BC (zero traction at far-field)
            attr = BP2BoundaryAttributes::FARFIELD_LEFT;
         }
         else if (std::abs(x_center - params.Lx) < tol)
         {
            // Right boundary (x = +Lx): Natural BC (zero traction at far-field)
            attr = BP2BoundaryAttributes::FARFIELD_RIGHT;
         }
         else if (std::abs(z_center) < tol)
         {
            // Top boundary (z = 0): Free surface with Natural BC
            attr = BP2BoundaryAttributes::FREE_SURFACE;
         }
         else if (std::abs(z_center + params.Lz) < tol)
         {
            // Bottom boundary (z = -Lz): Dirichlet BC with position-dependent plate loading
            attr = BP2BoundaryAttributes::BOTTOM;
         }
         else
         {
            // Should not happen for rectangular mesh
            MFEM_WARNING("Unexpected boundary element position: ("
                         << x_center << ", " << z_center << ")");
            attr = BP2BoundaryAttributes::FARFIELD_LEFT;  // Default
         }

         mesh.SetBdrAttribute(be, attr);
      }

      // Update boundary attributes array
      mesh.SetAttributes();
   }

   /// @brief Apply geometric grading to concentrate mesh near fault (x=0) and surface (z=0)
   ///
   /// After coordinate shift, domain is x ∈ [-Lx, +Lx], z ∈ [-Lz, 0]
   static void ApplyGrading(Mesh &mesh, const Parameters &params)
   {
      // Get mesh nodes
      GridFunction *nodes = mesh.GetNodes();
      if (!nodes)
      {
         mesh.SetCurvature(1);
         nodes = mesh.GetNodes();
      }

      const int dim = mesh.Dimension();
      const int nnodes = nodes->Size() / dim;

      for (int i = 0; i < nnodes; i++)
      {
         real_t x = (*nodes)(i * dim);
         real_t z = (*nodes)(i * dim + 1);

         // Normalize to [-1, 1] for x, [-1, 0] for z
         real_t xi = x / params.Lx;           // xi ∈ [-1, +1]
         real_t zeta = z / params.Lz;         // zeta ∈ [-1, 0] (z is negative)

         // Apply grading in x (concentrate near x=0, i.e., the fault)
         // Using sinh-based grading: preserves [-1, 1] range
         if (params.grading_x != 1.0)
         {
            real_t alpha_x = params.grading_x;
            xi = std::sinh(alpha_x * xi) / std::sinh(alpha_x);
         }

         // Apply grading in z (concentrate near z=0, the surface/fault zone)
         // zeta is in [-1, 0], sinh grading concentrates elements near zeta=0
         if (params.grading_z != 1.0)
         {
            real_t alpha_z = params.grading_z;
            zeta = std::sinh(alpha_z * zeta) / std::sinh(alpha_z);
         }

         // Transform back to physical coordinates
         (*nodes)(i * dim) = xi * params.Lx;
         (*nodes)(i * dim + 1) = zeta * params.Lz;
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BP2_MESH_HPP
