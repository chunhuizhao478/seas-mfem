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

#ifndef MFEM_PFF_SOLVER_HPP
#define MFEM_PFF_SOLVER_HPP

#include "mfem.hpp"
#include "materials/pff_material.hpp"
#include "integrators/degraded_elasticity_integrator.hpp"
#include "operators/elasticity_operator.hpp"
#include "operators/damage_operator.hpp"

#include <memory>
#include <cmath>
#include <iomanip>

namespace mfem
{

namespace pff
{

/** @brief Coefficient that returns g(d) * base_value, where g is degradation.
 *
 * g(d) = (1-d)^p * (1-eta) + eta
 *
 * This is used to degrade elastic moduli based on the damage field.
 */
class DegradedCoefficient : public Coefficient
{
public:
   DegradedCoefficient(GridFunction *d_gf, real_t base_value, real_t eta, real_t p = 2.0)
      : d_gf_(d_gf), base_value_(base_value), eta_(eta), p_(p) {}

   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
      real_t d_val = d_gf_->GetValue(T, ip);
      d_val = std::max(0.0, std::min(1.0, d_val));
      real_t g = std::pow(1.0 - d_val, p_) * (1.0 - eta_) + eta_;
      return g * base_value_;
   }

private:
   GridFunction *d_gf_;
   real_t base_value_;
   real_t eta_;
   real_t p_;
};

/** @brief Staggered solver for phase field fracture.
 *
 * Implements the staggered (alternate minimization) scheme:
 * 1. Given d^k, solve elasticity for u^{k+1}
 * 2. Given u^{k+1}, update strain energy history H
 * 3. Solve damage equation for d^{k+1}
 * 4. Check convergence and repeat
 *
 * The solver uses Newton's method with GMRES and AMG preconditioner
 * for scalability to large parallel problems.
 */
class PFFSolver
{
public:
   /** @brief Construct PFF solver.
    *
    * @param[in] mesh Parallel mesh
    * @param[in] mat Material parameters
    * @param[in] order Finite element order (default: 1)
    */
   PFFSolver(ParMesh &mesh, const PFFMaterialParameters &mat, int order = 1)
      : pmesh_(mesh), mat_(mat), order_(order), dim_(mesh.Dimension()),
        t_(0.0), prescribed_disp_(nullptr)
   {
      // Create finite element collections
      u_fec_ = std::make_unique<H1_FECollection>(order_, dim_);
      d_fec_ = std::make_unique<H1_FECollection>(order_, dim_);
      psi_fec_ = std::make_unique<L2_FECollection>(order_ - 1, dim_);

      // Create finite element spaces
      u_fes_ = std::make_unique<ParFiniteElementSpace>(&pmesh_, u_fec_.get(), dim_);
      d_fes_ = std::make_unique<ParFiniteElementSpace>(&pmesh_, d_fec_.get());
      psi_fes_ = std::make_unique<ParFiniteElementSpace>(&pmesh_, psi_fec_.get());

      // Create solution GridFunctions
      u_gf_ = std::make_unique<ParGridFunction>(u_fes_.get());
      d_gf_ = std::make_unique<ParGridFunction>(d_fes_.get());
      d_old_gf_ = std::make_unique<ParGridFunction>(d_fes_.get());
      psi_gf_ = std::make_unique<ParGridFunction>(psi_fes_.get());

      // Initialize to zero
      *u_gf_ = 0.0;
      *d_gf_ = 0.0;
      *d_old_gf_ = 0.0;
      *psi_gf_ = 0.0;

      // Default essential boundary (none)
      ess_bdr_u_.SetSize(pmesh_.bdr_attributes.Max());
      ess_bdr_u_ = 0;

      // Note: Operators and solvers are set up lazily when first needed
      // This avoids issues with uninitialized boundary conditions
      operators_initialized_ = false;
      solvers_initialized_ = false;

      // Print info - must call global functions on all ranks before printing
      long long global_ne = pmesh_.GetGlobalNE();
      HYPRE_BigInt u_dofs = u_fes_->GlobalTrueVSize();
      HYPRE_BigInt d_dofs = d_fes_->GlobalTrueVSize();
      if (Mpi::Root())
      {
         mfem::out << "PFFSolver initialized:\n"
                   << "  Mesh elements: " << global_ne << "\n"
                   << "  Displacement DOFs: " << u_dofs << "\n"
                   << "  Damage DOFs: " << d_dofs << "\n"
                   << "  FE order: " << order_ << "\n";
      }
   }

   ~PFFSolver() = default;

   // =========================================================================
   // Boundary Condition Configuration
   // =========================================================================

   /** @brief Set essential boundary attributes for displacement.
    *
    * @param[in] bdr_attr Array marking essential boundary attributes (1 = essential)
    */
   void SetElasticityEssentialBdr(const Array<int> &bdr_attr)
   {
      MFEM_VERIFY(bdr_attr.Size() == pmesh_.bdr_attributes.Max(),
                  "Boundary attribute array size mismatch");
      ess_bdr_u_ = bdr_attr;

      // Rebuild operators with new BCs
      SetupOperators();
   }

   /** @brief Set prescribed displacement coefficient.
    *
    * The coefficient should implement SetTime() for time-dependent loading.
    *
    * @param[in] u_bc Vector coefficient for prescribed displacement
    */
   void SetPrescribedDisplacement(VectorCoefficient &u_bc)
   {
      prescribed_disp_ = &u_bc;
   }

   /** @brief Configure simple tension boundary conditions.
    *
    * Bottom (attribute 1): fixed in y
    * Top (attribute 3): prescribed displacement in y
    * Left corner: fixed in x (to prevent rigid body motion)
    */
   void ConfigureSimpleTensionBCs()
   {
      // Essential boundaries: bottom and top
      ess_bdr_u_ = 0;
      if (pmesh_.bdr_attributes.Max() >= 1) { ess_bdr_u_[0] = 1; }  // bottom
      if (pmesh_.bdr_attributes.Max() >= 3) { ess_bdr_u_[2] = 1; }  // top

      // Get essential true DOFs
      u_fes_->GetEssentialTrueDofs(ess_bdr_u_, ess_tdof_u_);

      // Create component-wise markers for selective BC application
      const int tsize = u_fes_->GetTrueVSize();
      // Number of nodes (MFEM uses byNODES ordering by default)
      const int num_nodes = u_fes_->GetTrueVSize() / dim_;

      marker_top_y_.SetSize(tsize);
      marker_bot_y_.SetSize(tsize);
      marker_bot_x_.SetSize(tsize);
      marker_top_y_ = 0;
      marker_bot_y_ = 0;
      marker_bot_x_ = 0;

      // Identify DOFs on boundaries by coordinate
      // Project coordinate function to get reliable node positions
      VectorFunctionCoefficient coord_coeff(dim_,
         [](const Vector &x, Vector &v) { v = x; });
      ParGridFunction coords(u_fes_.get());
      coords.ProjectCoefficient(coord_coeff);

      Vector coords_true;
      coords.GetTrueDofs(coords_true);

      // Find domain y-bounds
      // Note: MFEM uses byNODES ordering: DOFs are [x0,y0,x1,y1,...] for 2D
      real_t y_min = std::numeric_limits<real_t>::max();
      real_t y_max = std::numeric_limits<real_t>::lowest();

      for (int i = 0; i < num_nodes; i++)
      {
         // byNODES ordering: y-component is at index i*dim + 1
         real_t y = coords_true(i * dim_ + 1);
         y_min = std::min(y_min, y);
         y_max = std::max(y_max, y);
      }

      // MPI reduce to get global bounds
      real_t global_y_min, global_y_max;
      MPI_Allreduce(&y_min, &global_y_min, 1, MPI_DOUBLE, MPI_MIN, pmesh_.GetComm());
      MPI_Allreduce(&y_max, &global_y_max, 1, MPI_DOUBLE, MPI_MAX, pmesh_.GetComm());

      const real_t tol = 1e-10 * (global_y_max - global_y_min);

      // Mark DOFs using byNODES ordering
      for (int i = 0; i < num_nodes; i++)
      {
         // byNODES ordering: for node i, x-DOF is at i*dim, y-DOF is at i*dim+1
         real_t y = coords_true(i * dim_ + 1);

         // Bottom boundary: fix both x and y components to prevent rigid body motion
         if (std::abs(y - global_y_min) < tol)
         {
            marker_bot_y_[i * dim_ + 1] = 1;  // Fix y on bottom
            marker_bot_x_[i * dim_] = 1;      // Also fix x on bottom
         }

         // Top y-component DOFs (prescribed displacement)
         if (std::abs(y - global_y_max) < tol)
         {
            marker_top_y_[i * dim_ + 1] = 1;
         }
      }

      // Combine markers into essential DOF list
      ess_tdof_u_.SetSize(0);
      for (int i = 0; i < tsize; i++)
      {
         if (marker_bot_y_[i] || marker_bot_x_[i] || marker_top_y_[i])
         {
            ess_tdof_u_.Append(i);
         }
      }

      // Reset initialization flags so operators will be rebuilt with new BCs
      operators_initialized_ = false;
      solvers_initialized_ = false;
      // Reset preconditioners so they get rebuilt
      amg_u_.reset();
      amg_d_.reset();

      if (Mpi::Root())
      {
         mfem::out << "Simple tension BCs configured:\n"
                   << "  Bottom y-DOFs fixed: " << marker_bot_y_.Sum() << "\n"
                   << "  Bottom x-DOFs fixed: " << marker_bot_x_.Sum() << "\n"
                   << "  Top y-DOFs prescribed: " << marker_top_y_.Sum() << "\n";
      }
   }

   /** @brief Configure Mode 1 fracture boundary conditions (Raccoon tutorial).
    *
    * This reproduces the boundary conditions from the Raccoon mode1_brittle_fracture
    * tutorial:
    * - Top boundary (y = y_max): u_y = prescribed (tension), u_x = 0 (fixed)
    * - Bottom right (y = y_min, x >= x_mid): u_y = 0 (support)
    * - Bottom left (y = y_min, x < x_mid): FREE (represents pre-crack)
    *
    * @param[in] x_crack_tip X-coordinate of crack tip (default: 0.5 for unit domain)
    */
   void ConfigureMode1FractureBCs(real_t x_crack_tip = 0.5)
   {
      // Get mesh bounds
      Vector bb_min, bb_max;
      pmesh_.GetBoundingBox(bb_min, bb_max);
      real_t global_x_min = bb_min(0), global_x_max = bb_max(0);
      real_t global_y_min = bb_min(1), global_y_max = bb_max(1);

      // Use reasonable tolerance
      const real_t tol = 1e-8 * std::max(global_x_max - global_x_min,
                                          global_y_max - global_y_min);

      // Scale crack tip coordinate to actual domain if given as fraction
      real_t crack_x = x_crack_tip;
      if (x_crack_tip <= 1.0 && global_x_max > 1.0)
      {
         crack_x = global_x_min + x_crack_tip * (global_x_max - global_x_min);
      }

      // Use MFEM's coordinate coefficient to mark DOFs
      // This approach iterates through local DOFs and marks based on coordinates

      const int ndof = u_fes_->GetNDofs();  // Number of scalar DOFs locally
      const int vdim = u_fes_->GetVDim();   // Should be dim_

      // Initialize markers (local DOF size = ndof * vdim)
      const int local_size = ndof * vdim;
      marker_top_y_.SetSize(local_size);
      marker_top_x_.SetSize(local_size);
      marker_bot_y_.SetSize(local_size);
      marker_bot_x_.SetSize(local_size);
      marker_top_y_ = 0;
      marker_top_x_ = 0;
      marker_bot_y_ = 0;
      marker_bot_x_ = 0;

      int num_top = 0, num_bot = 0;

      // Iterate over mesh vertices (for order 1, DOFs are at vertices)
      const int nv = pmesh_.GetNV();
      Array<int> vdofs;

      for (int v = 0; v < nv; v++)
      {
         const real_t *vc = pmesh_.GetVertex(v);
         real_t x = vc[0];
         real_t y = vc[1];

         // Find the DOFs at this vertex
         // For H1 with vdim components, vertex v maps to DOFs based on ordering
         // With byNODES: DOF for component c at vertex v is v*vdim + c

         int dof_x = v * vdim;      // x-component DOF
         int dof_y = v * vdim + 1;  // y-component DOF

         // Check if within local DOF range
         if (dof_x >= local_size || dof_y >= local_size) { continue; }

         // Top boundary: fix x, prescribe y
         if (std::abs(y - global_y_max) < tol)
         {
            marker_top_x_[dof_x] = 1;
            marker_top_y_[dof_y] = 1;
            num_top++;
         }

         // Bottom right (noncrack): fix y only
         if (std::abs(y - global_y_min) < tol && x >= crack_x - tol)
         {
            marker_bot_y_[dof_y] = 1;
            num_bot++;
         }
      }

      // Convert local DOF markers to true DOF list
      ess_tdof_u_.SetSize(0);
      Array<int> ldof_marker(local_size);
      ldof_marker = 0;
      for (int i = 0; i < local_size; i++)
      {
         if (marker_top_x_[i] || marker_top_y_[i] || marker_bot_y_[i])
         {
            ldof_marker[i] = 1;
         }
      }

      // Use GetEssentialTrueDofs-style conversion
      // Actually, we can directly convert local DOF indices to true DOF indices
      for (int i = 0; i < local_size; i++)
      {
         if (ldof_marker[i])
         {
            int tdof = u_fes_->GetLocalTDofNumber(i);
            if (tdof >= 0)  // Only add if this is a true DOF on this process
            {
               ess_tdof_u_.Append(tdof);
            }
         }
      }

      // Sort and remove duplicates
      ess_tdof_u_.Sort();
      ess_tdof_u_.Unique();

      // Reset initialization flags
      operators_initialized_ = false;
      solvers_initialized_ = false;
      amg_u_.reset();
      amg_d_.reset();

      // Global count
      int global_top, global_bot;
      MPI_Allreduce(&num_top, &global_top, 1, MPI_INT, MPI_SUM, pmesh_.GetComm());
      MPI_Allreduce(&num_bot, &global_bot, 1, MPI_INT, MPI_SUM, pmesh_.GetComm());

      if (Mpi::Root())
      {
         mfem::out << "Mode 1 fracture BCs configured:\n"
                   << "  Domain: [" << global_x_min << ", " << global_x_max
                   << "] x [" << global_y_min << ", " << global_y_max << "]\n"
                   << "  Crack tip x-coordinate: " << crack_x << "\n"
                   << "  Top boundary nodes (u_x=0, u_y=prescribed): " << global_top << "\n"
                   << "  Bottom right nodes (u_y=0, support): " << global_bot << "\n"
                   << "  Total essential true DOFs: " << ess_tdof_u_.Size() << "\n";
      }
   }

   /** @brief Configure Mode 2 (shear) fracture boundary conditions.
    *
    * This reproduces the boundary conditions from the Raccoon mode2_brittle_fracture
    * tutorial:
    * - Domain: [0,1] x [-0.5, 0.5] with horizontal crack at y=0, x < crack_tip
    * - Top boundary (y = y_max): u_x = prescribed (shear), u_y = 0 (fixed)
    * - Bottom boundary (y = y_min): u_x = 0 (fixed), u_y = 0 (fixed)
    *
    * Note: The pre-existing crack is represented by the mesh geometry (nodes not
    * connected at y=0 for x < crack_tip). For a conforming mesh without crack,
    * the phase field will evolve to represent the crack.
    */
   void ConfigureMode2FractureBCs()
   {
      // Use MFEM's boundary attribute approach which properly handles parallel DOFs.
      // For Mode 2 shear:
      // - Top (attribute 3): u_x = prescribed, u_y = 0
      // - Bottom (attribute 1): u_x = 0, u_y = 0
      //
      // Key insight: After mesh refinement, DOF numbering changes. We must use
      // boundary attributes to identify essential DOFs, NOT coordinate-based markers.

      int max_bdr_attr = pmesh_.bdr_attributes.Max();

      // Get essential TRUE DOFs for top boundary (attribute 3)
      Array<int> top_bdr(max_bdr_attr);
      top_bdr = 0;
      if (max_bdr_attr >= 3) { top_bdr[2] = 1; }
      u_fes_->GetEssentialTrueDofs(top_bdr, ess_tdof_top_);

      // Get essential TRUE DOFs for bottom boundary (attribute 1)
      Array<int> bot_bdr(max_bdr_attr);
      bot_bdr = 0;
      if (max_bdr_attr >= 1) { bot_bdr[0] = 1; }
      u_fes_->GetEssentialTrueDofs(bot_bdr, ess_tdof_bot_);

      // Combined essential DOFs for FormLinearSystem
      ess_bdr_u_.SetSize(max_bdr_attr);
      ess_bdr_u_ = 0;
      if (max_bdr_attr >= 1) { ess_bdr_u_[0] = 1; }  // bottom
      if (max_bdr_attr >= 3) { ess_bdr_u_[2] = 1; }  // top
      u_fes_->GetEssentialTrueDofs(ess_bdr_u_, ess_tdof_u_);

      // For byNODES ordering, identify x-component DOFs in the top boundary list
      // True DOF layout: [x-comp DOFs (0 to N-1)] [y-comp DOFs (N to 2N-1)]
      // where N = TrueVSize() / vdim
      const int vdim = u_fes_->GetVDim();
      const int scalar_tdof_size = u_fes_->GetTrueVSize() / vdim;

      // Build list of top boundary x-component TRUE DOFs
      ess_tdof_top_x_.SetSize(0);
      for (int i = 0; i < ess_tdof_top_.Size(); i++)
      {
         int tdof = ess_tdof_top_[i];
         // For byNODES: tdof < scalar_tdof_size means x-component
         if (tdof < scalar_tdof_size)
         {
            ess_tdof_top_x_.Append(tdof);
         }
      }

      // Reset initialization flags
      operators_initialized_ = false;
      solvers_initialized_ = false;
      amg_u_.reset();
      amg_d_.reset();

      // Count for reporting - also count x vs y components
      int num_top_x = ess_tdof_top_x_.Size();
      int num_top = ess_tdof_top_.Size();
      int num_bot = ess_tdof_bot_.Size();

      // Count x and y components in bottom boundary
      int num_bot_x = 0, num_bot_y = 0;
      for (int i = 0; i < ess_tdof_bot_.Size(); i++)
      {
         int tdof = ess_tdof_bot_[i];
         if (tdof < scalar_tdof_size) { num_bot_x++; }
         else { num_bot_y++; }
      }

      // Global counts - all processes must participate in collective operations
      int global_top_x, global_top, global_bot, global_bot_x, global_bot_y;
      MPI_Allreduce(&num_top_x, &global_top_x, 1, MPI_INT, MPI_SUM, pmesh_.GetComm());
      MPI_Allreduce(&num_top, &global_top, 1, MPI_INT, MPI_SUM, pmesh_.GetComm());
      MPI_Allreduce(&num_bot, &global_bot, 1, MPI_INT, MPI_SUM, pmesh_.GetComm());
      MPI_Allreduce(&num_bot_x, &global_bot_x, 1, MPI_INT, MPI_SUM, pmesh_.GetComm());
      MPI_Allreduce(&num_bot_y, &global_bot_y, 1, MPI_INT, MPI_SUM, pmesh_.GetComm());

      // Get bounding box on all processes (may involve MPI communication)
      Vector bb_min, bb_max;
      pmesh_.GetBoundingBox(bb_min, bb_max);

      if (Mpi::Root())
      {
         mfem::out << "Mode 2 fracture BCs configured:\n"
                   << "  Domain: [" << bb_min(0) << ", " << bb_max(0)
                   << "] x [" << bb_min(1) << ", " << bb_max(1) << "]\n"
                   << "  Top boundary true DOFs: " << global_top
                   << " (x-component: " << global_top_x << ")\n"
                   << "  Bottom boundary true DOFs: " << global_bot
                   << " (x: " << global_bot_x << ", y: " << global_bot_y << ")\n"
                   << "  Total essential true DOFs: " << ess_tdof_u_.Size() << "\n"
                   << "  Scalar true DOF size: " << scalar_tdof_size << "\n";
      }
   }

   // =========================================================================
   // Time Stepping
   // =========================================================================

   /** @brief Solve one time step using staggered iteration (MOOSE-style Picard).
    *
    * Implements MOOSE-style Picard iteration between elasticity and damage:
    * - Alternates between solving elasticity (given d) and damage (given u)
    * - Tracks both |u_new - u_old| and |d_new - d_old| for convergence
    * - Controlled by fixed_point_max_its, fixed_point_rel_tol, fixed_point_abs_tol
    * - accept_on_max_fixed_point_iteration controls whether to accept unconverged solution
    *
    * @param[in] dt Time step size
    * @param[out] num_iters Optional pointer to receive number of iterations
    * @return true if converged, false if max iterations reached without convergence
    */
   bool SolveTimeStep(real_t dt, int *num_iters = nullptr)
   {
      // Ensure operators are initialized before solving
      EnsureInitialized();

      // Update prescribed displacement time
      if (prescribed_disp_)
      {
         prescribed_disp_->SetTime(t_ + dt);
      }

      // Store u and d from previous staggered iteration for fixed-point convergence
      // (MOOSE checks both u and d change for convergence)
      Vector u_prev_iter;
      u_gf_->GetTrueDofs(u_prev_iter);
      Vector d_prev_iter;
      d_gf_->GetTrueDofs(d_prev_iter);

      int iter = 0;
      bool converged = false;
      real_t fp_res_u = 0.0, fp_res_d = 0.0;
      real_t fp_rel_u = 0.0, fp_rel_d = 0.0;

      while (iter < max_fp_iter_ && !converged)
      {
         // =====================================================================
         // Elasticity sub-problem (MOOSE-style output)
         // =====================================================================
         if (Mpi::Root() && print_level_ > 0)
         {
            mfem::out << " Elasticity system solve:\n";
         }

         // Compute initial residual before solve
         real_t u_res_initial = ComputeElasticityResidual();
         if (Mpi::Root() && print_level_ > 0)
         {
            mfem::out << "  * Nonlinear |R| = " << std::scientific
                      << std::setprecision(6) << u_res_initial
                      << " (before solve)" << std::defaultfloat << "\n";
         }

         // Solve elasticity
         SolveElasticity();

         // Compute final residual after solve
         real_t u_res_final = ComputeElasticityResidual();
         if (Mpi::Root() && print_level_ > 0)
         {
            mfem::out << "  0 Nonlinear |R| = " << std::scientific
                      << std::setprecision(6) << u_res_final
                      << std::defaultfloat << "\n";
         }

         // =====================================================================
         // Update strain energy history
         // =====================================================================
         ComputeAndUpdateStrainEnergy();

         // =====================================================================
         // Damage sub-problem (MOOSE-style output)
         // =====================================================================
         if (Mpi::Root() && print_level_ > 0)
         {
            mfem::out << " Damage system solve:\n";
         }

         // Compute initial residual before solve
         real_t d_res_initial = ComputeDamageResidual();
         if (Mpi::Root() && print_level_ > 0)
         {
            mfem::out << "  * Nonlinear |R| = " << std::scientific
                      << std::setprecision(6) << d_res_initial
                      << " (before solve)" << std::defaultfloat << "\n";
         }

         // Solve damage
         SolveDamage();

         // Apply irreversibility (enforce d >= d_old from previous time step)
         Vector d_true;
         d_gf_->GetTrueDofs(d_true);
         Vector d_old_timestep;
         d_old_gf_->GetTrueDofs(d_old_timestep);
         damage_op_->ApplyIrreversibility(d_true, d_old_timestep);
         d_gf_->SetFromTrueDofs(d_true);

         // Compute final residual after solve
         real_t d_res_final = ComputeDamageResidual();
         if (Mpi::Root() && print_level_ > 0)
         {
            mfem::out << "  0 Nonlinear |R| = " << std::scientific
                      << std::setprecision(6) << d_res_final
                      << std::defaultfloat << "\n";
         }

         // =====================================================================
         // Fixed-point convergence check (MOOSE-style: check BOTH u and d)
         // =====================================================================

         // Get current u
         Vector u_true;
         u_gf_->GetTrueDofs(u_true);

         // Compute |u_new - u_old| for fixed-point convergence
         Vector diff_u(u_true);
         diff_u -= u_prev_iter;
         real_t local_norm_sq_u = diff_u * diff_u;
         real_t global_norm_sq_u;
         MPI_Allreduce(&local_norm_sq_u, &global_norm_sq_u, 1, MPI_DOUBLE, MPI_SUM,
                       pmesh_.GetComm());
         fp_res_u = std::sqrt(global_norm_sq_u);

         // Compute |u_new| for relative convergence
         real_t u_norm_local = u_true * u_true;
         real_t u_norm_global;
         MPI_Allreduce(&u_norm_local, &u_norm_global, 1, MPI_DOUBLE, MPI_SUM,
                       pmesh_.GetComm());
         real_t u_norm = std::sqrt(u_norm_global);
         fp_rel_u = fp_res_u / (u_norm + 1e-14);

         // Compute |d_new - d_old| for fixed-point convergence
         Vector diff_d(d_true);
         diff_d -= d_prev_iter;
         real_t local_norm_sq_d = diff_d * diff_d;
         real_t global_norm_sq_d;
         MPI_Allreduce(&local_norm_sq_d, &global_norm_sq_d, 1, MPI_DOUBLE, MPI_SUM,
                       pmesh_.GetComm());
         fp_res_d = std::sqrt(global_norm_sq_d);

         // Compute |d_new| for relative convergence
         real_t d_norm_local = d_true * d_true;
         real_t d_norm_global;
         MPI_Allreduce(&d_norm_local, &d_norm_global, 1, MPI_DOUBLE, MPI_SUM,
                       pmesh_.GetComm());
         real_t d_norm = std::sqrt(d_norm_global);
         fp_rel_d = fp_res_d / (d_norm + 1e-14);

         // Combined fixed-point residual (MOOSE uses max of both)
         real_t fp_res = std::max(fp_res_u, fp_res_d);
         real_t fp_rel = std::max(fp_rel_u, fp_rel_d);

         // MOOSE-style convergence: either absolute or relative (BOTH u and d must converge)
         bool u_converged = (fp_res_u < fp_abs_tol_) || (fp_rel_u < fp_rel_tol_);
         bool d_converged = (fp_res_d < fp_abs_tol_) || (fp_rel_d < fp_rel_tol_);
         converged = u_converged && d_converged;

         if (Mpi::Root() && print_level_ > 0)
         {
            mfem::out << " Fixed point (u): |delta_u| = " << std::scientific
                      << std::setprecision(6) << fp_res_u
                      << ", |delta_u|/|u| = " << fp_rel_u << std::defaultfloat;
            if (u_converged) { mfem::out << " [CONVERGED]"; }
            mfem::out << "\n";

            mfem::out << " Fixed point (d): |delta_d| = " << std::scientific
                      << std::setprecision(6) << fp_res_d
                      << ", |delta_d|/|d| = " << fp_rel_d << std::defaultfloat;
            if (d_converged) { mfem::out << " [CONVERGED]"; }
            mfem::out << "\n";

            if (converged)
            {
               mfem::out << " Fixed point iteration CONVERGED at iteration " << (iter + 1) << "\n";
            }
            mfem::out << "\n";
         }

         // Update for next iteration
         u_prev_iter = u_true;
         d_prev_iter = d_true;
         iter++;
      }

      // Handle non-convergence (like MOOSE's accept_on_max_fixed_point_iteration)
      if (!converged)
      {
         if (Mpi::Root() && print_level_ >= 0)
         {
            mfem::out << " WARNING: Fixed-point iteration did not converge in "
                      << max_fp_iter_ << " iterations\n";
            mfem::out << "   Final |delta_u| = " << std::scientific << fp_res_u
                      << ", |delta_u|/|u| = " << fp_rel_u << "\n";
            mfem::out << "   Final |delta_d| = " << fp_res_d
                      << ", |delta_d|/|d| = " << fp_rel_d << std::defaultfloat << "\n";

            if (!accept_on_max_fp_iter_)
            {
               mfem::out << " Solution NOT accepted (accept_on_max_fixed_point_iteration=false)\n";
            }
            else
            {
               mfem::out << " Solution ACCEPTED despite non-convergence "
                         << "(accept_on_max_fixed_point_iteration=true)\n";
            }
         }
      }

      if (num_iters) { *num_iters = iter; }
      return converged || accept_on_max_fp_iter_;
   }

   /** @brief Backward-compatible version that returns iteration count.
    *
    * @param[in] dt Time step size
    * @return Number of staggered iterations
    */
   int SolveTimeStepIterations(real_t dt)
   {
      int iters = 0;
      SolveTimeStep(dt, &iters);
      return iters;
   }

   /** @brief Compute the true nonlinear residual for elasticity: R = K*u - f.
    *
    * This is the actual PDE residual, not the linear solver residual.
    * For elasticity: R = ∫ σ(u,d):ε(v) dx - ∫ f·v dx - ∫ t·v ds
    *
    * @return Global L2 norm of the residual
    */
   real_t ComputeElasticityResidual()
   {
      // Use DegradedElasticityIntegrator with SPECTRAL DECOMPOSITION
      // to compute true PDE residual: R = ∫ σ(u,d) : ε(v) dΩ
      // where σ = g(d)*σ⁺ + σ⁻ (only tensile stress is degraded)
      ParNonlinearForm nlf(u_fes_.get());
      nlf.AddDomainIntegrator(new DegradedElasticityIntegrator(mat_, *d_gf_));

      // Get current displacement as true DOF vector
      Vector u_true(u_fes_->GetTrueVSize());
      u_gf_->GetTrueDofs(u_true);

      // Compute residual: R(u) = internal forces
      Vector R(u_fes_->GetTrueVSize());
      nlf.Mult(u_true, R);

      // Zero out essential DOF contributions (they are satisfied exactly)
      for (int i = 0; i < ess_tdof_u_.Size(); i++)
      {
         R(ess_tdof_u_[i]) = 0.0;
      }

      // Compute global L2 norm
      real_t local_norm_sq = R * R;
      real_t global_norm_sq;
      MPI_Allreduce(&local_norm_sq, &global_norm_sq, 1, MPI_DOUBLE, MPI_SUM,
                    pmesh_.GetComm());

      return std::sqrt(global_norm_sq);
   }

   /** @brief Compute the true nonlinear residual for damage: R = A*d - b.
    *
    * This is the actual PDE residual for the damage equation.
    * R = ∫ (Gc*l/c0)*∇d·∇v dx + ∫ (2*Gc/(c0*l))*d*v dx - ∫ (-g'(d)*H)*v dx
    *
    * @return Global L2 norm of the residual
    */
   real_t ComputeDamageResidual()
   {
      // Assemble A and b for damage equation
      // Use direct ParallelAssemble to get true PDE residual
      ParBilinearForm a(d_fes_.get());
      ParLinearForm b(d_fes_.get());

      // Diffusion coefficient: Gc * l / c0
      real_t diff_coeff = mat_.Gc * mat_.l / mat_.c0;
      ConstantCoefficient kappa(diff_coeff);

      // Reaction coefficient: 2 * Gc / (c0 * l)
      real_t react_coeff = 2.0 * mat_.Gc / (mat_.c0 * mat_.l);
      ConstantCoefficient sigma(react_coeff);

      a.AddDomainIntegrator(new DiffusionIntegrator(kappa));
      a.AddDomainIntegrator(new MassIntegrator(sigma));
      a.Assemble();
      a.Finalize();

      // Source term: -g'(d)*H
      DamageSourceCoefficient src_coeff(d_gf_.get(),
                                        &damage_op_->GetStrainEnergyHistory(),
                                        mat_.eta);
      b.AddDomainIntegrator(new DomainLFIntegrator(src_coeff));
      b.Assemble();

      // Get parallel matrix A (no BC elimination for damage)
      OperatorPtr A_ptr;
      A_ptr.Reset(a.ParallelAssemble());
      HypreParMatrix *A = A_ptr.As<HypreParMatrix>();

      // Get parallel RHS vector B
      Vector B(d_fes_->GetTrueVSize());
      b.ParallelAssemble(B);

      // Get current damage as true DOF vector
      Vector d_true(d_fes_->GetTrueVSize());
      d_gf_->GetTrueDofs(d_true);

      // Compute residual: R = A*d - B
      Vector R(d_fes_->GetTrueVSize());
      A->Mult(d_true, R);
      R -= B;

      // Compute global L2 norm using parallel inner product
      HypreParVector R_hyp(A->GetComm(), A->GetGlobalNumRows(),
                           R.GetData(), A->GetRowStarts());
      real_t norm_sq = InnerProduct(R_hyp, R_hyp);

      return std::sqrt(norm_sq);
   }

   /** @brief Advance time after a converged step.
    *
    * @param[in] dt Time step size
    */
   void AdvanceTime(real_t dt)
   {
      t_ += dt;
      *d_old_gf_ = *d_gf_;
   }

   // =========================================================================
   // Accessors
   // =========================================================================

   /// Get current displacement field
   const ParGridFunction &GetDisplacement() const { return *u_gf_; }
   ParGridFunction &GetDisplacement() { return *u_gf_; }

   /// Get current damage field
   const ParGridFunction &GetDamage() const { return *d_gf_; }
   ParGridFunction &GetDamage() { return *d_gf_; }

   /// Get strain energy field
   const ParGridFunction &GetStrainEnergy() const { return *psi_gf_; }

   /// Get strain energy history field (H = max(ψ⁺) over time)
   const ParGridFunction &GetStrainEnergyHistory() const
   {
      return damage_op_->GetStrainEnergyHistory();
   }

   /// Get displacement FE space
   const ParFiniteElementSpace &GetDisplacementFES() const { return *u_fes_; }

   /// Get damage FE space
   const ParFiniteElementSpace &GetDamageFES() const { return *d_fes_; }

   /// Get current time
   real_t GetTime() const { return t_; }

   /// Get boundary markers (for debugging)
   const Array<int> &GetMaskTopY() const { return marker_top_y_; }
   const Array<int> &GetMaskBotY() const { return marker_bot_y_; }
   const Array<int> &GetMaskBotX() const { return marker_bot_x_; }

   /** @brief Compute reaction force on the top boundary using geometry detection.
    *
    * Computes R = K*u (internal force) and sums the specified component
    * over the DOFs at y = y_max.
    *
    * For mode 2 shear: component=0 (x) gives shear reaction.
    * For mode 1 tension: component=1 (y) gives tensile reaction.
    *
    * @param[in] component Displacement component (0=x, 1=y)
    * @return Total reaction force on top boundary (summed over all processes)
    */
   real_t ComputeReactionForceOnTop(int component)
   {
      EnsureInitialized();

      // Assemble the stiffness matrix (without BC elimination)
      ParBilinearForm a(u_fes_.get());

      DegradedCoefficient lambda_coeff(d_gf_.get(), mat_.lambda, mat_.eta, mat_.p);
      DegradedCoefficient mu_coeff(d_gf_.get(), mat_.mu, mat_.eta, mat_.p);

      a.AddDomainIntegrator(new ElasticityIntegrator(lambda_coeff, mu_coeff));
      a.Assemble();
      a.Finalize();

      // Get the parallel matrix
      HypreParMatrix *K = a.ParallelAssemble();

      // Get displacement as true DOF vector
      Vector u_true;
      u_gf_->GetTrueDofs(u_true);

      // Compute internal force: f_int = K * u
      Vector f_int(u_true.Size());
      K->Mult(u_true, f_int);

      delete K;

      // Get domain bounds
      Vector bb_min, bb_max;
      pmesh_.GetBoundingBox(bb_min, bb_max);
      real_t y_max = bb_max(1);
      real_t tol = 1e-6 * (bb_max(1) - bb_min(1));

      // Sum reaction forces at top boundary (y = y_max) for the given component
      // Use geometry-based detection instead of boundary markers
      real_t local_reaction = 0.0;

      // Get node coordinates
      int ndof = u_fes_->GetNDofs();  // Scalar DOFs per component
      ParGridFunction coords(u_fes_.get());
      VectorFunctionCoefficient coord_coeff(dim_,
         [](const Vector &x, Vector &v) { v = x; });
      coords.ProjectCoefficient(coord_coeff);

      // Iterate over scalar DOFs and find those at y_max
      for (int i = 0; i < ndof; i++)
      {
         // Get y-coordinate of this DOF (stored in second component)
         real_t y_coord = coords(ndof + i);

         // Check if this DOF is on the top boundary
         if (std::abs(y_coord - y_max) < tol)
         {
            // Get the true DOF index for the specified component
            // For vector FE space: DOF i of component c is at local index c*ndof + i
            int vdof = component * ndof + i;

            // Convert to true DOF
            int tdof = u_fes_->GetLocalTDofNumber(vdof);
            if (tdof >= 0)  // Only process locally owned DOFs
            {
               local_reaction += f_int(tdof);
            }
         }
      }

      // Global sum across all processes
      real_t global_reaction;
      MPI_Allreduce(&local_reaction, &global_reaction, 1, MPI_DOUBLE, MPI_SUM,
                    pmesh_.GetComm());

      return global_reaction;
   }

   /** @brief Compute x-direction reaction force on top boundary (for mode 2 shear).
    * @return Shear reaction force Fx on top boundary
    */
   real_t ComputeShearReactionForce()
   {
      return ComputeReactionForceOnTop(0);  // x-component on top
   }

   /** @brief Compute y-direction reaction force on top boundary (for mode 1 tension).
    * @return Tensile reaction force Fy on top boundary
    */
   real_t ComputeTensileReactionForce()
   {
      return ComputeReactionForceOnTop(1);  // y-component on top
   }

   // =========================================================================
   // Solver Parameters
   // =========================================================================

   /// Set maximum staggered iterations
   void SetMaxIterations(int max_it) { max_fp_iter_ = max_it; }

   /// Set staggered convergence tolerance (relative)
   void SetRelativeTolerance(real_t tol) { fp_rel_tol_ = tol; }

   /// Set staggered convergence tolerance (absolute)
   void SetAbsoluteTolerance(real_t tol) { fp_abs_tol_ = tol; }

   /// Set Newton solver tolerances
   void SetNewtonRelTol(real_t tol)
   {
      newton_u_->SetRelTol(tol);
      newton_d_->SetRelTol(tol);
   }

   void SetNewtonAbsTol(real_t tol)
   {
      newton_u_->SetAbsTol(tol);
      newton_d_->SetAbsTol(tol);
   }

   void SetNewtonMaxIter(int max_it)
   {
      newton_u_->SetMaxIter(max_it);
      newton_d_->SetMaxIter(max_it);
   }

   /// Set print level (0 = silent, 1 = summary, 2 = verbose)
   void SetPrintLevel(int level) { print_level_ = level; }

   /// Set whether to accept solution when max fixed-point iterations reached
   /// (Like MOOSE's accept_on_max_fixed_point_iteration)
   void SetAcceptOnMaxFixedPointIteration(bool accept) { accept_on_max_fp_iter_ = accept; }

private:
   // =========================================================================
   // Internal Methods
   // =========================================================================

   void SetupOperators()
   {
      // Create damage operator (used for strain energy history tracking)
      damage_op_ = std::make_unique<DamageOperator>(*d_fes_, mat_);
   }

   void SetupNewtonSolvers()
   {
      // Newton solvers are no longer used - using direct linear solves instead
      // This function is kept for interface compatibility
   }

   void EnsureInitialized() const
   {
      if (!operators_initialized_)
      {
         const_cast<PFFSolver*>(this)->SetupOperators();
         operators_initialized_ = true;
      }
      if (!solvers_initialized_)
      {
         const_cast<PFFSolver*>(this)->SetupNewtonSolvers();
         solvers_initialized_ = true;
      }
   }

   /** @brief Solve elasticity sub-problem with SPECTRAL DECOMPOSITION.
    *
    * Uses Newton iteration because spectral decomposition makes the
    * problem nonlinear (stress depends on sign of strain eigenvalues).
    * This matches Raccoon's approach: only tensile stress is degraded.
    *
    * @return Final residual norm from the Newton solver
    */
   real_t SolveElasticity()
   {
      EnsureInitialized();

      // Build prescribed displacement vector for BCs
      Vector u_prescribed(u_fes_->GetTrueVSize());
      u_prescribed = 0.0;

      if (prescribed_disp_)
      {
         real_t shear_val = 0.0;

         Vector bb_min, bb_max;
         pmesh_.GetBoundingBox(bb_min, bb_max);
         real_t y_max = bb_max(1);
         real_t tol = 1e-6 * (y_max - bb_min(1));

         for (int be = 0; be < pmesh_.GetNBE(); be++)
         {
            ElementTransformation *bdr_tr = pmesh_.GetBdrElementTransformation(be);
            if (!bdr_tr) { continue; }

            IntegrationPoint bdr_center;
            bdr_center.x = 0.5;
            bdr_center.y = 0.0;
            bdr_center.z = 0.0;
            bdr_center.weight = 1.0;
            bdr_tr->SetIntPoint(&bdr_center);

            Vector bdr_pt(dim_);
            bdr_tr->Transform(bdr_center, bdr_pt);

            if (std::abs(bdr_pt(1) - y_max) < tol)
            {
               Vector V(dim_);
               prescribed_disp_->Eval(V, *bdr_tr, bdr_center);
               shear_val = V(0);
               break;
            }
         }

         real_t global_shear_val;
         MPI_Allreduce(&shear_val, &global_shear_val, 1, MPI_DOUBLE, MPI_MAX,
                       pmesh_.GetComm());
         shear_val = global_shear_val;

         for (int i = 0; i < ess_tdof_top_x_.Size(); i++)
         {
            int tdof = ess_tdof_top_x_[i];
            u_prescribed(tdof) = shear_val;
         }
      }

      // Create nonlinear form with DegradedElasticityIntegrator
      // This uses SPECTRAL DECOMPOSITION: σ = g(d)*σ⁺ + σ⁻
      // Only tensile stress is degraded, matching Raccoon
      ParNonlinearForm nlf(u_fes_.get());
      nlf.AddDomainIntegrator(new DegradedElasticityIntegrator(mat_, *d_gf_));
      nlf.SetEssentialTrueDofs(ess_tdof_u_);

      // Set initial guess from previous solution
      Vector u_true(u_fes_->GetTrueVSize());
      u_gf_->GetTrueDofs(u_true);

      // Apply prescribed displacements to initial guess
      for (int i = 0; i < ess_tdof_u_.Size(); i++)
      {
         u_true(ess_tdof_u_[i]) = u_prescribed(ess_tdof_u_[i]);
      }

      // Newton iteration
      const int max_newton_iter = 10;
      const real_t newton_rel_tol = 1e-8;
      const real_t newton_abs_tol = 1e-12;

      Vector r(u_fes_->GetTrueVSize());
      Vector du(u_fes_->GetTrueVSize());

      real_t res_norm = 0.0;
      real_t res_norm_0 = 0.0;

      for (int newton_iter = 0; newton_iter < max_newton_iter; newton_iter++)
      {
         // Compute residual: r = -R(u)
         nlf.Mult(u_true, r);
         r.Neg();

         // Enforce essential BCs on residual: r_i = u_prescribed_i - u_i for ess DOFs
         for (int i = 0; i < ess_tdof_u_.Size(); i++)
         {
            int tdof = ess_tdof_u_[i];
            r(tdof) = u_prescribed(tdof) - u_true(tdof);
         }

         // Compute residual norm
         real_t local_norm_sq = r * r;
         real_t global_norm_sq;
         MPI_Allreduce(&local_norm_sq, &global_norm_sq, 1, MPI_DOUBLE, MPI_SUM,
                       pmesh_.GetComm());
         res_norm = std::sqrt(global_norm_sq);

         if (newton_iter == 0) { res_norm_0 = res_norm; }

         // Check convergence
         if (res_norm < newton_abs_tol ||
             (res_norm_0 > 0 && res_norm / res_norm_0 < newton_rel_tol))
         {
            break;
         }

         // Get Jacobian K = dR/du
         Operator &K_op = nlf.GetGradient(u_true);
         HypreParMatrix *K = dynamic_cast<HypreParMatrix*>(&K_op);
         MFEM_VERIFY(K != nullptr, "Failed to get Jacobian as HypreParMatrix");

         // Eliminate essential DOFs from Jacobian
         HypreParMatrix *K_elim = K->EliminateRowsCols(ess_tdof_u_);
         delete K_elim;

         // Solve K * du = r with AMG preconditioned GMRES
         HypreBoomerAMG amg(*K);
         amg.SetPrintLevel(0);
         amg.SetElasticityOptions(u_fes_.get());

         GMRESSolver gmres(pmesh_.GetComm());
         gmres.SetRelTol(1e-10);
         gmres.SetAbsTol(1e-14);
         gmres.SetMaxIter(200);
         gmres.SetPrintLevel(0);
         gmres.SetPreconditioner(amg);
         gmres.SetOperator(*K);

         du = 0.0;
         gmres.Mult(r, du);

         // Update solution: u = u + du
         u_true += du;
      }

      // Set solution back to grid function
      u_gf_->SetFromTrueDofs(u_true);

      return res_norm;
   }

   /** @brief Compute eigenvalues of 2D symmetric strain tensor.
    *
    * For strain tensor:
    *   ε = [ε₁₁  ε₁₂]
    *       [ε₁₂  ε₂₂]
    *
    * Eigenvalues: λ₁,₂ = (ε₁₁+ε₂₂)/2 ± sqrt[((ε₁₁-ε₂₂)/2)² + ε₁₂²]
    *
    * @param[in] strain Symmetric 2x2 strain matrix
    * @param[out] eig_vals Vector of 2 eigenvalues (sorted: λ₁ ≤ λ₂)
    */
   void ComputeStrainEigenvalues2D(const DenseMatrix &strain, Vector &eig_vals) const
   {
      eig_vals.SetSize(2);

      real_t e11 = strain(0, 0);
      real_t e22 = strain(1, 1);
      real_t e12 = strain(0, 1);

      real_t trace = e11 + e22;
      real_t diff = e11 - e22;
      real_t disc = std::sqrt(0.25 * diff * diff + e12 * e12);

      // Eigenvalues (λ₁ ≤ λ₂)
      eig_vals(0) = 0.5 * trace - disc;
      eig_vals(1) = 0.5 * trace + disc;
   }

   /** @brief Compute spectral decomposition strain energy (tensile part only).
    *
    * Spectral decomposition splits strain energy into tensile (positive eigenvalues)
    * and compressive (negative eigenvalues) parts. Only tensile part drives damage.
    *
    * For AT2 model:
    *   ψ⁺ = 0.5 * λ * <tr(ε)>₊² + μ * ε⁺:ε⁺
    *
    * where:
    *   <x>₊ = max(x, 0)
    *   ε⁺ = Σᵢ <λᵢ>₊ nᵢ⊗nᵢ (positive part of strain)
    *   ε⁺:ε⁺ = Σᵢ <λᵢ>₊²
    *
    * @param[in] strain Symmetric strain matrix
    * @param[in] dim Spatial dimension (2 or 3)
    * @return Tensile strain energy density ψ⁺
    */
   real_t ComputeSpectralStrainEnergy(const DenseMatrix &strain, int dim) const
   {
      if (dim == 2)
      {
         Vector eig_vals;
         ComputeStrainEigenvalues2D(strain, eig_vals);

         // Positive eigenvalues
         real_t lam1_pos = std::max(eig_vals(0), 0.0);
         real_t lam2_pos = std::max(eig_vals(1), 0.0);

         // Trace terms
         real_t tr_eps = eig_vals(0) + eig_vals(1);
         real_t tr_eps_pos = std::max(tr_eps, 0.0);

         // ε⁺:ε⁺ = λ₁⁺² + λ₂⁺²
         real_t eps_pos_sq = lam1_pos * lam1_pos + lam2_pos * lam2_pos;

         // ψ⁺ = 0.5*λ*<tr(ε)>₊² + μ*ε⁺:ε⁺
         real_t psi_pos = 0.5 * mat_.lambda * tr_eps_pos * tr_eps_pos +
                          mat_.mu * eps_pos_sq;

         return psi_pos;
      }
      else // dim == 3
      {
         // For 3D, use MFEM's eigenvalue solver
         DenseMatrix strain_copy(strain);
         Vector eig_vals(3);
         DenseMatrix eig_vecs(3, 3);
         strain_copy.Eigenvalues(eig_vals, eig_vecs);

         // Positive eigenvalues
         real_t lam1_pos = std::max(eig_vals(0), 0.0);
         real_t lam2_pos = std::max(eig_vals(1), 0.0);
         real_t lam3_pos = std::max(eig_vals(2), 0.0);

         // Trace terms
         real_t tr_eps = eig_vals(0) + eig_vals(1) + eig_vals(2);
         real_t tr_eps_pos = std::max(tr_eps, 0.0);

         // ε⁺:ε⁺
         real_t eps_pos_sq = lam1_pos * lam1_pos + lam2_pos * lam2_pos +
                             lam3_pos * lam3_pos;

         // ψ⁺ = 0.5*λ*<tr(ε)>₊² + μ*ε⁺:ε⁺
         real_t psi_pos = 0.5 * mat_.lambda * tr_eps_pos * tr_eps_pos +
                          mat_.mu * eps_pos_sq;

         return psi_pos;
      }
   }

   void ComputeAndUpdateStrainEnergy()
   {
      // Compute strain energy density from displacement field using
      // SPECTRAL decomposition (only tensile part drives damage).
      //
      // This matches Raccoon's decomposition = SPECTRAL option.

      const int NE = u_fes_->GetNE();
      *psi_gf_ = 0.0;

      DenseMatrix grad_u, strain;
      Vector elfun;
      Array<int> vdofs, psi_dofs;

      for (int el = 0; el < NE; el++)
      {
         const FiniteElement &u_fe = *u_fes_->GetFE(el);
         ElementTransformation &Tr = *u_fes_->GetElementTransformation(el);

         u_fes_->GetElementVDofs(el, vdofs);
         u_gf_->GetSubVector(vdofs, elfun);

         int nd = u_fe.GetDof();
         int dim = u_fe.GetDim();

         // Compute strain at element center
         const IntegrationRule &ir = IntRules.Get(u_fe.GetGeomType(), 2);
         real_t psi_avg = 0.0;
         real_t vol = 0.0;

         DenseMatrix dshape(nd, dim), dshapedxt(nd, dim);

         for (int i = 0; i < ir.GetNPoints(); i++)
         {
            const IntegrationPoint &ip = ir.IntPoint(i);
            Tr.SetIntPoint(&ip);

            u_fe.CalcDShape(ip, dshape);
            Mult(dshape, Tr.InverseJacobian(), dshapedxt);

            // Compute gradient of u
            grad_u.SetSize(dim, dim);
            grad_u = 0.0;
            for (int j = 0; j < dim; j++)
            {
               for (int k = 0; k < dim; k++)
               {
                  for (int l = 0; l < nd; l++)
                  {
                     grad_u(j, k) += elfun(l + j * nd) * dshapedxt(l, k);
                  }
               }
            }

            // Compute symmetric strain tensor
            strain.SetSize(dim, dim);
            for (int j = 0; j < dim; j++)
            {
               for (int k = 0; k < dim; k++)
               {
                  strain(j, k) = 0.5 * (grad_u(j, k) + grad_u(k, j));
               }
            }

            // Compute tensile strain energy using spectral decomposition
            real_t psi = ComputeSpectralStrainEnergy(strain, dim);

            real_t w = ip.weight * Tr.Weight();
            psi_avg += w * psi;
            vol += w;
         }

         psi_avg /= (vol > 0 ? vol : 1.0);

         // Set element value in psi_gf
         psi_fes_->GetElementDofs(el, psi_dofs);
         for (int i = 0; i < psi_dofs.Size(); i++)
         {
            (*psi_gf_)(psi_dofs[i]) = psi_avg;
         }
      }

      // Update history in damage operator
      damage_op_->UpdateStrainEnergyHistory(*psi_gf_);
   }

   /** @brief Coefficient for damage source term: -g'(d)*H = 2*(1-d)*(1-η)*H.
    *
    * This is the correct driving force for damage evolution. It naturally
    * reduces to zero as d→1 (fully damaged), preventing over-driving in
    * already damaged regions.
    */
   class DamageSourceCoefficient : public Coefficient
   {
   public:
      DamageSourceCoefficient(GridFunction *d_gf, GridFunction *H_gf,
                              real_t eta)
         : d_gf_(d_gf), H_gf_(H_gf), eta_(eta) {}

      real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
         real_t d_val = d_gf_->GetValue(T, ip);
         d_val = std::max(0.0, std::min(1.0, d_val));
         real_t H_val = H_gf_->GetValue(T, ip);

         // -g'(d) = 2*(1-d)*(1-η) for AT2 model with p=2
         // Source = -g'(d) * H = 2*(1-d)*(1-η)*H
         real_t neg_g_prime = 2.0 * (1.0 - d_val) * (1.0 - eta_);
         return neg_g_prime * H_val;
      }

   private:
      GridFunction *d_gf_;
      GridFunction *H_gf_;
      real_t eta_;
   };

   /** @brief Solve damage sub-problem with BOUNDED Newton (VI solver).
    *
    * Enforces bounds during solve (like Raccoon's vinewtonrsls):
    *   - Lower bound: d >= d_old (irreversibility)
    *   - Upper bound: d <= 1
    *
    * Uses projected Newton iteration. After each linear solve, the solution
    * is projected onto the feasible set [d_old, 1]. This is a simplified
    * semi-smooth Newton approach that works well for the damage problem.
    *
    * @return Final residual norm
    */
   real_t SolveDamage()
   {
      EnsureInitialized();

      // Get bounds for damage: [d_old, 1]
      Vector d_lower(d_fes_->GetTrueVSize());
      Vector d_upper(d_fes_->GetTrueVSize());
      d_old_gf_->GetTrueDofs(d_lower);  // Lower bound = d from previous time step
      d_upper = 1.0;                     // Upper bound = 1

      // Get current damage as initial guess
      Vector d_true(d_fes_->GetTrueVSize());
      d_gf_->GetTrueDofs(d_true);

      // Project initial guess onto bounds
      for (int i = 0; i < d_true.Size(); i++)
      {
         d_true(i) = std::max(d_lower(i), std::min(d_upper(i), d_true(i)));
      }

      // Bounded Newton iteration parameters
      const int max_vi_iter = 20;
      const real_t vi_rel_tol = 1e-8;
      const real_t vi_abs_tol = 1e-12;

      real_t res_norm = 0.0;
      real_t res_norm_0 = 0.0;

      for (int vi_iter = 0; vi_iter < max_vi_iter; vi_iter++)
      {
         // Update grid function from true DOFs for coefficient evaluation
         d_gf_->SetFromTrueDofs(d_true);

         // Assemble the damage system
         ParBilinearForm a(d_fes_.get());
         ParLinearForm b(d_fes_.get());

         // Diffusion coefficient: Gc * l / c0
         real_t diff_coeff = mat_.Gc * mat_.l / mat_.c0;
         ConstantCoefficient kappa(diff_coeff);

         // Reaction coefficient: 2 * Gc / (c0 * l)
         real_t react_coeff = 2.0 * mat_.Gc / (mat_.c0 * mat_.l);
         ConstantCoefficient sigma(react_coeff);

         a.AddDomainIntegrator(new DiffusionIntegrator(kappa));
         a.AddDomainIntegrator(new MassIntegrator(sigma));
         a.Assemble();
         a.Finalize();

         // Source term: -g'(d)*H = 2*(1-d)*(1-η)*H
         DamageSourceCoefficient src_coeff(d_gf_.get(),
                                           &damage_op_->GetStrainEnergyHistory(),
                                           mat_.eta);
         b.AddDomainIntegrator(new DomainLFIntegrator(src_coeff));
         b.Assemble();

         // Get parallel matrix and RHS
         OperatorPtr A_ptr;
         A_ptr.Reset(a.ParallelAssemble());
         HypreParMatrix *A = A_ptr.As<HypreParMatrix>();

         Vector B(d_fes_->GetTrueVSize());
         b.ParallelAssemble(B);

         // Compute residual: r = A*d - B
         Vector r(d_fes_->GetTrueVSize());
         A->Mult(d_true, r);
         r -= B;

         // Apply complementarity conditions for VI:
         // Residual r = A*d - B:
         //   r > 0 means d is too high, solution wants to DECREASE d
         //   r < 0 means d is too low, solution wants to INCREASE d
         // At LOWER bound: if r > 0 (wants to decrease), can't → active constraint
         // At UPPER bound: if r < 0 (wants to increase), can't → active constraint
         real_t bound_tol = 1e-12;

         for (int i = 0; i < d_true.Size(); i++)
         {
            bool at_lower = (d_true(i) - d_lower(i)) < bound_tol;
            bool at_upper = (d_upper(i) - d_true(i)) < bound_tol;

            if (at_lower && r(i) > 0.0)
            {
               r(i) = 0.0;  // Can't decrease below lower bound
            }
            else if (at_upper && r(i) < 0.0)
            {
               r(i) = 0.0;  // Can't increase above upper bound
            }
         }

         // Compute residual norm
         real_t local_norm_sq = r * r;
         real_t global_norm_sq;
         MPI_Allreduce(&local_norm_sq, &global_norm_sq, 1, MPI_DOUBLE, MPI_SUM,
                       pmesh_.GetComm());
         res_norm = std::sqrt(global_norm_sq);

         if (vi_iter == 0) { res_norm_0 = res_norm; }

         // Check convergence
         if (res_norm < vi_abs_tol ||
             (res_norm_0 > 0 && res_norm / res_norm_0 < vi_rel_tol))
         {
            break;
         }

         // Solve for new damage: A * d_new = B
         // (This is a direct solve, not incremental)
         HypreBoomerAMG amg(*A);
         amg.SetPrintLevel(0);

         CGSolver cg(pmesh_.GetComm());
         cg.SetRelTol(1e-10);
         cg.SetAbsTol(1e-14);
         cg.SetMaxIter(500);
         cg.SetPrintLevel(0);
         cg.SetPreconditioner(amg);
         cg.SetOperator(*A);

         Vector d_new(d_fes_->GetTrueVSize());
         d_new = 0.0;
         cg.Mult(B, d_new);

         // Project onto bounds [d_lower, d_upper]
         for (int i = 0; i < d_new.Size(); i++)
         {
            d_new(i) = std::max(d_lower(i), std::min(d_upper(i), d_new(i)));
         }

         // Update d_true
         d_true = d_new;
      }

      // Set final solution
      d_gf_->SetFromTrueDofs(d_true);

      return res_norm;
   }

   bool CheckConvergence(const Vector &d_new, const Vector &d_old) const
   {
      Vector diff(d_new);
      diff -= d_old;

      real_t norm_diff = diff.Norml2();
      real_t norm_d = d_new.Norml2();

      // Global reduction for parallel
      real_t global_norm_diff, global_norm_d;
      MPI_Allreduce(&norm_diff, &global_norm_diff, 1, MPI_DOUBLE, MPI_SUM,
                    pmesh_.GetComm());
      MPI_Allreduce(&norm_d, &global_norm_d, 1, MPI_DOUBLE, MPI_SUM,
                    pmesh_.GetComm());

      global_norm_diff = std::sqrt(global_norm_diff);
      global_norm_d = std::sqrt(global_norm_d);

      real_t rel_change = global_norm_diff / (global_norm_d + 1e-14);

      return (rel_change < fp_rel_tol_) || (global_norm_diff < fp_abs_tol_);
   }

   // =========================================================================
   // Member Variables
   // =========================================================================

   // Mesh and material
   ParMesh &pmesh_;
   const PFFMaterialParameters &mat_;
   int order_;
   int dim_;

   // Finite element collections and spaces
   std::unique_ptr<H1_FECollection> u_fec_;
   std::unique_ptr<H1_FECollection> d_fec_;
   std::unique_ptr<L2_FECollection> psi_fec_;

   std::unique_ptr<ParFiniteElementSpace> u_fes_;
   std::unique_ptr<ParFiniteElementSpace> d_fes_;
   std::unique_ptr<ParFiniteElementSpace> psi_fes_;

   // Solution fields
   std::unique_ptr<ParGridFunction> u_gf_;
   std::unique_ptr<ParGridFunction> d_gf_;
   std::unique_ptr<ParGridFunction> d_old_gf_;
   std::unique_ptr<ParGridFunction> psi_gf_;

   // Operators
   std::unique_ptr<ElasticityOperator> elasticity_op_;
   std::unique_ptr<DamageOperator> damage_op_;

   // Newton solvers and linear solvers
   std::unique_ptr<NewtonSolver> newton_u_;
   std::unique_ptr<NewtonSolver> newton_d_;
   std::unique_ptr<GMRESSolver> gmres_u_;
   std::unique_ptr<CGSolver> cg_d_;
   std::unique_ptr<HypreBoomerAMG> amg_u_;
   std::unique_ptr<HypreBoomerAMG> amg_d_;

   // Boundary conditions
   Array<int> ess_bdr_u_;
   Array<int> ess_tdof_u_;       // Combined essential true DOFs (top + bottom)
   Array<int> ess_tdof_top_;     // Essential true DOFs on top boundary (all components)
   Array<int> ess_tdof_bot_;     // Essential true DOFs on bottom boundary (all components)
   Array<int> ess_tdof_top_x_;   // Essential true DOFs for x-component on top boundary

   // Local DOF markers for boundary conditions (used by ConfigureSimpleTensionBCs and ConfigureMode1FractureBCs)
   Array<int> marker_top_y_;     // Marker for y-component DOFs on top boundary
   Array<int> marker_top_x_;     // Marker for x-component DOFs on top boundary
   Array<int> marker_bot_y_;     // Marker for y-component DOFs on bottom boundary
   Array<int> marker_bot_x_;     // Marker for x-component DOFs on bottom boundary

   VectorCoefficient *prescribed_disp_;

   // Solver parameters (MOOSE-style fixed-point iteration)
   int max_fp_iter_ = 50;
   real_t fp_rel_tol_ = 1e-6;
   real_t fp_abs_tol_ = 1e-10;
   int print_level_ = 1;
   bool accept_on_max_fp_iter_ = false;  // Like MOOSE's accept_on_max_fixed_point_iteration

   // Lazy initialization flags
   mutable bool operators_initialized_;
   mutable bool solvers_initialized_;

   // Time
   real_t t_;
};

} // namespace pff

} // namespace mfem

#endif // MFEM_PFF_SOLVER_HPP
