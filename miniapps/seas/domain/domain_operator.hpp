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

#ifndef MFEM_SEAS_DOMAIN_OPERATOR_HPP
#define MFEM_SEAS_DOMAIN_OPERATOR_HPP

#include "mfem.hpp"
#include "../common/seas_types.hpp"

namespace mfem
{
namespace seas
{

// Forward declaration — only a pointer is returned, no include needed.
class FaultBasis;

/// @brief Abstract base class for domain operators
///
/// This class defines the interface for solving the domain PDE
/// (e.g., Laplace equation for antiplane shear) with fault slip
/// boundary conditions, and computing traction on the fault.
///
/// The template parameter MeshType allows the same interface to
/// work with both serial (Mesh) and parallel (ParMesh) execution.
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class DomainOperator
{
public:
   // Type aliases derived from mesh type
   using FESpaceType = FESpaceForMesh<MeshType>;
   using GridFuncType = GridFunctionForMesh<MeshType>;
   using BilinFormType = BilinearFormForMesh<MeshType>;
   using LinFormType = LinearFormForMesh<MeshType>;

   virtual ~DomainOperator() = default;

   /// @brief Get number of displacement components
   ///
   /// Returns 1 for antiplane (scalar), 2 for plane strain, 3 for 3D elasticity.
   virtual int NumComponents() const = 0;

   /// @brief Get the spatial dimension of the problem
   ///
   /// Returns 2 for 2D problems (antiplane, plane strain), 3 for 3D.
   virtual int Dimension() const = 0;

   /// @brief Solve domain problem with given fault slip boundary condition
   ///
   /// Given the current time and slip values on the fault, solve the
   /// domain equilibrium equation (e.g., ∇²u = 0 for antiplane) with:
   /// - Fault boundary: u = slip/2 on x=0⁺, u = -slip/2 on x=0⁻
   /// - Far-field: u = ±Vp*t/2 at x → ±∞
   /// - Free surface: σ·n = 0 at z = 0
   ///
   /// @param[in] time Current simulation time
   /// @param[in] slip_bc Slip values at fault DOFs
   /// @param[out] displacement Solution displacement field
   virtual void Solve(real_t time, const Vector &slip_bc,
                      GridFuncType &displacement) = 0;

   /// @brief Compute traction at fault from displacement field
   ///
   /// Computes the full DG numerical flux traction:
   ///   τ = μ * {∂u/∂x} + penalty * ([[u]] - δ)
   /// where {·} is the average, [[·]] is the jump, and δ is the slip BC.
   ///
   /// @param[in] displacement Current displacement field
   /// @param[in] slip_bc Slip boundary condition at fault DOFs
   /// @param[out] traction Computed traction at fault DOFs
   virtual void ComputeTraction(const GridFuncType &displacement,
                                const Vector &slip_bc,
                                Vector &traction) = 0;

   /// @brief Get reference to finite element space
   virtual FESpaceType &GetFESpace() = 0;

   /// @brief Get const reference to finite element space
   virtual const FESpaceType &GetFESpace() const = 0;

   /// @brief Get the underlying mesh
   virtual MeshType &GetMesh() = 0;

   /// @brief Get const reference to the mesh
   virtual const MeshType &GetMesh() const = 0;

   /// @brief Get shear modulus
   virtual real_t GetShearModulus() const = 0;

   /// @brief Get the number of DOFs on the fault boundary
   virtual int GetNumFaultDOFs() const = 0;

   /// @brief Get the depth (z-coordinate) at each fault DOF
   ///
   /// This is needed for applying depth-dependent friction parameters.
   virtual void GetFaultDepths(Vector &depths) const = 0;

   /// @brief Get fault DOF indices in the global finite element space
   ///
   /// These indices can be used to apply boundary conditions.
   virtual const Array<int> &GetFaultDOFs() const = 0;

   /// @brief Get number of slip components per fault DOF
   ///
   /// Returns 1 for antiplane (scalar slip), 2 for 3D elasticity (strike + dip).
   virtual int NumSlipComponents() const { return 1; }

   /// @brief Get 2D fault coordinates (along-strike and depth) at each fault DOF
   ///
   /// For 2D problems (antiplane), the along-strike coordinate is zero and
   /// the depth coordinate matches GetFaultDepths(). For 3D problems (BP5),
   /// this returns the (x2, x3) coordinates of each fault DOF on the 2D
   /// fault surface.
   ///
   /// @param[out] coords_x2 Along-strike coordinates [NumFaultDOFs]
   /// @param[out] coords_x3 Depth coordinates [NumFaultDOFs]
   virtual void GetFaultCoords2D(Vector &coords_x2, Vector &coords_x3) const
   {
      Vector depths;
      GetFaultDepths(depths);
      coords_x2.SetSize(depths.Size());
      coords_x2 = 0.0;
      coords_x3 = depths;
   }

   /// @brief Get the per-face fault basis. Returns nullptr by default.
   ///
   /// Override in 3D elasticity operators to provide the global-to-local
   /// coordinate transformation on the fault surface.
   virtual const FaultBasis *GetFaultBasis() const { return nullptr; }

   /// @brief Get off-fault displacement at specified spatial points
   ///
   /// Evaluates the current displacement field at given 3D points.
   /// Required by BP5 for off-fault station output. Default: no-op.
   ///
   /// @param[in] points Vector of spatial coordinates (each is a dim-Vector)
   /// @param[out] displacements Displacement values at each point
   virtual void GetOffFaultDisplacement(
      const std::vector<Vector> &points, Vector &displacements) const
   {
      displacements.SetSize(0);
   }

#ifdef MFEM_USE_MPI
   /// Get MPI communicator (parallel only)
   virtual MPI_Comm GetComm() const { return MPI_COMM_WORLD; }
#endif
};

// Convenience type aliases
using SerialDomainOperator = DomainOperator<Mesh>;
#ifdef MFEM_USE_MPI
using ParallelDomainOperator = DomainOperator<ParMesh>;
#endif

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DOMAIN_OPERATOR_HPP
