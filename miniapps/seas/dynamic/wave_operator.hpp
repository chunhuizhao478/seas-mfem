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

#ifndef MFEM_SEAS_WAVE_OPERATOR_HPP
#define MFEM_SEAS_WAVE_OPERATOR_HPP

#include "mfem.hpp"
#include "wave_state.hpp"
#include "godunov_flux.hpp"
#include "pml_layer.hpp"
#include "../domain/boundary_config.hpp"
#include "../fault/fault_basis.hpp"

#include <memory>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief DG wave operator for the 3D velocity-stress elastic wave equation.
///
/// Solves dQ/dt + A dQ/dx + B dQ/dy + C dQ/dz = 0 using DG with
/// Godunov upwind flux. Inherits only TimeDependentOperator (R-001 fix),
/// NOT DomainOperator.
///
/// Owns its own mesh, L2 FE space, per-element inverse mass matrix,
/// and Godunov flux. Shared fault-coupling code (FaultBasis) is accessed
/// via composition.
///
/// The state vector Q has 9 components per DOF:
///   [sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz, v_x, v_y, v_z]
/// stored as Q[c * ndof_total + local_dof] (component-major).
///
/// Reference: Dumbser & Kaser (2006), de la Puente et al. (2009).
class WaveOperator : public TimeDependentOperator
{
public:
   /// @brief Construct a WaveOperator on the given mesh.
   ///
   /// @param[in] mesh  Mesh (serial or parallel). WaveOperator does NOT own this.
   /// @param[in] order  Polynomial order for L2 DG space.
   /// @param[in] lambda  First Lame parameter [Pa].
   /// @param[in] mu  Shear modulus [Pa].
   /// @param[in] rho  Density [kg/m^3].
   /// @param[in] bc  Boundary configuration (which attrs are absorbing/free/etc).
   WaveOperator(Mesh &mesh, int order,
                real_t lambda, real_t mu, real_t rho,
                const BoundaryConfig &bc);

   /// Destructor.
   ~WaveOperator() override;

   /// @brief Compute dQ/dt = (M^{-1}) * (-Face + Vol).
   ///
   /// Implements Eq. (3) from the plan: the full semi-discrete ODE RHS.
   /// This is the main entry point called by the time integrator.
   void Mult(const Vector &Q, Vector &dQdt) const override;

   /// @name Accessors
   ///@{
   int GetOrder() const { return order_; }
   int GetNDof() const { return ndof_per_el_; }
   int NumElements() const { return ne_; }

   const GodunovFlux &GetFlux() const { return flux_; }
   const FiniteElementSpace &GetFESpace() const { return *fes_; }

   /// Number of DOFs per scalar component (total across all elements).
   int GetScalarNDof() const { return ndof_total_; }

   /// Compute CFL-limited time step: dt = cfl * h_min / c_p.
   real_t ComputeMaxDt(real_t cfl) const;

   /// Get the FaultBasis (may be null if no fault faces found).
   const FaultBasis *GetFaultBasis() const { return fault_basis_.get(); }

   /// Get number of fault DOFs per component.
   int GetNumFaultDOFs() const { return num_fault_dofs_; }

   /// Get per-element inverse mass matrix (for testing).
   const DenseMatrix &GetElementMassInverse(int e) const { return elem_mass_inv_[e]; }

   /// Set PML layer (non-owning). Pass nullptr to disable.
   void SetPML(PMLLayer *pml) { pml_layer_ = pml; }
   const PMLLayer *GetPML() const { return pml_layer_; }
   ///@}

private:
   Mesh &mesh_;
   int order_;
   int ndof_per_el_;  ///< DOFs per element per component
   int ne_;           ///< Number of elements
   int ndof_total_;   ///< Total DOFs per component = ne_ * ndof_per_el_

   std::unique_ptr<L2_FECollection> fec_;
   std::unique_ptr<FiniteElementSpace> fes_;

   GodunovFlux flux_;
   BoundaryConfig bc_;

   /// Precomputed Jacobian matrices for x, y, z directions (constant for homogeneous).
   DenseMatrix Ax_, Ay_, Az_;

   /// Per-element inverse mass matrix (dense ndof x ndof).
   std::vector<DenseMatrix> elem_mass_inv_;

   /// Fault basis (constructed if mesh has fault faces).
   std::unique_ptr<FaultBasis> fault_basis_;
   int num_fault_dofs_ = 0;

   /// Minimum inscribed element diameter (for CFL).
   real_t h_min_;

   /// Per-face boundary attribute (0 = interior/shared face).
   std::vector<int> face_bdr_attr_;

   /// Optional PML layer (non-owning, nullptr if disabled).
   PMLLayer *pml_layer_ = nullptr;

   /// @brief Compute volume integral contribution to RHS.
   /// Implements Eq. (2c): for each element, accumulate
   ///   rhs[c*ndof+i] += w * dshape(i,j) * F[j][c]
   void ComputeVolumeRHS(const Vector &Q, Vector &rhs) const;

   /// @brief Compute face flux contribution to RHS.
   /// Implements Eq. (2d): for each face, compute Godunov flux and
   /// subtract from both elements' RHS.
   void ComputeFaceFluxRHS(const Vector &Q, Vector &rhs) const;

   /// @brief Apply per-element inverse mass matrix to the RHS.
   /// dQdt[c*ndof+i] = sum_j M_inv[i][j] * rhs[c*ndof+j], per element.
   void ApplyMassInverse(Vector &dQdt) const;

   /// @brief Assemble and store per-element inverse mass matrices.
   void AssembleElementMassInverse();

   /// @brief Apply PML damping: rhs -= d(x) * D * Q (Eq. 16).
   /// Per-component damping with directional splitting (R-004 fix).
   void ApplyPMLDamping(const Vector &Q, Vector &rhs) const;

   /// Classify a boundary face: returns "absorbing", "free", or "interior".
   enum class FaceBC { Interior, Absorbing, FreeSurface, Fault };
   FaceBC ClassifyBoundaryFace(int bdr_attr) const;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_WAVE_OPERATOR_HPP
