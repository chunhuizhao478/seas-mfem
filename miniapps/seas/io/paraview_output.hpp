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

#ifndef MFEM_SEAS_PARAVIEW_OUTPUT_HPP
#define MFEM_SEAS_PARAVIEW_OUTPUT_HPP

#include "mfem.hpp"
#include "../config/bp2_params.hpp"

#include <string>
#include <cmath>
#include <memory>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Manages ParaView-compatible visualization output for SEAS simulations.
///
/// Outputs all fields (displacement and fault fields) in a single
/// ParaViewDataCollection on the 2D domain mesh. Fault fields (slip, slip
/// rate, shear stress, state variable) are stored as piecewise-constant (L2
/// p=0) GridFunctions on the domain mesh: elements adjacent to the fault
/// carry fault values, all other elements are zero.
///
/// Output frequency controlled by step interval (every N time steps).
///
/// @tparam MeshType Either Mesh for serial or ParMesh for parallel
template <typename MeshType = Mesh>
class ParaViewOutput
{
public:
   /// @brief Output interval in time steps (write every N steps).
   /// Can be modified directly to change output frequency.
   int output_every_n_steps = 100;

   /// @brief Construct the ParaView output manager.
   ///
   /// @param prefix Output directory prefix (e.g., "ParaView/bp2")
   /// @param mesh The computational mesh
   /// @param order Finite element order (for levels of detail)
   ParaViewOutput(const std::string &prefix,
                  MeshType &mesh,
                  int order)
      : mesh_(mesh),
        order_(order),
        last_write_time_(-1e30),
        pv_(prefix, &mesh)
   {
      pv_.SetPrefixPath(prefix);
      pv_.SetDataFormat(VTKFormat::BINARY);
      pv_.SetHighOrderOutput(true);
      pv_.SetLevelsOfDetail(order);
   }

   /// @brief Register a domain field for output.
   ///
   /// @param name Field name (e.g., "displacement")
   /// @param gf Pointer to the GridFunction (must remain valid)
   void RegisterDomainField(const std::string &name, GridFunction *gf)
   {
      pv_.RegisterField(name, gf);
   }

   /// @brief Initialize fault output on the domain mesh.
   ///
   /// Creates L2 p=0 GridFunctions on the 2D domain mesh for fault fields.
   /// Elements adjacent to fault faces carry fault values; others are zero.
   ///
   /// @param fault_interior_faces Array of interior face indices on the fault
   void InitFaultOutput(const Array<int> &fault_interior_faces)
   {
      int nf = fault_interior_faces.Size();
      if (nf == 0) { return; }

      int dim = mesh_.Dimension();

      // L2 p=0 on the domain mesh: one DOF per element
      fault_fec_ = std::make_unique<L2_FECollection>(0, dim);
      fault_fes_ = std::make_unique<FiniteElementSpace>(
         &mesh_, fault_fec_.get());

      // Allocate grid functions
      fault_slip_ = std::make_unique<GridFunction>(fault_fes_.get());
      fault_slip_rate_ = std::make_unique<GridFunction>(fault_fes_.get());
      fault_stress_ = std::make_unique<GridFunction>(fault_fes_.get());
      fault_state_ = std::make_unique<GridFunction>(fault_fes_.get());

      *fault_slip_ = 0.0;
      *fault_slip_rate_ = 0.0;
      *fault_stress_ = 0.0;
      *fault_state_ = 0.0;

      // Build fault face -> element mapping
      // For each fault face, store the two adjacent element indices
      fault_face_elem1_.resize(nf);
      fault_face_elem2_.resize(nf);
      for (int i = 0; i < nf; i++)
      {
         int face = fault_interior_faces[i];
         FaceElementTransformations *FTr =
            mesh_.GetInteriorFaceTransformations(face);
         fault_face_elem1_[i] = FTr->Elem1No;
         fault_face_elem2_[i] = FTr->Elem2No;
      }

      // Register fault fields in the same collection
      pv_.RegisterField("slip", fault_slip_.get());
      pv_.RegisterField("slip_rate", fault_slip_rate_.get());
      pv_.RegisterField("shear_stress", fault_stress_.get());
      pv_.RegisterField("state_variable", fault_state_.get());

      has_fault_output_ = true;
   }

   /// @brief Update fault field data before saving.
   ///
   /// Maps fault DOF values to adjacent domain elements. Each fault face
   /// has (order+1) DOFs; we average them for the adjacent element value.
   ///
   /// @param slip Slip values at fault DOFs [m]
   /// @param slip_rate Slip rate at fault DOFs [m/s]
   /// @param stress Shear stress at fault DOFs [Pa]
   /// @param state State variable (theta) at fault DOFs [s]
   /// @param dofs_per_face Number of DOFs per fault face (1 = midpoint)
   void UpdateFaultFields(const Vector &slip,
                          const Vector &slip_rate,
                          const Vector &stress,
                          const Vector &state,
                          int dofs_per_face = 1)
   {
      if (!has_fault_output_) { return; }

      // Zero out all elements first
      *fault_slip_ = 0.0;
      *fault_slip_rate_ = 0.0;
      *fault_stress_ = 0.0;
      *fault_state_ = 0.0;

      int nf = static_cast<int>(fault_face_elem1_.size());

      // For each fault face, average the (dofs_per_face) DOFs and
      // assign to both adjacent elements
      for (int i = 0; i < nf; i++)
      {
         // Average DOFs for this face
         real_t avg_slip = 0.0, avg_rate = 0.0;
         real_t avg_stress = 0.0, avg_state = 0.0;
         int start = i * dofs_per_face;
         int count = 0;
         for (int j = 0; j < dofs_per_face && start + j < slip.Size(); j++)
         {
            avg_slip += slip(start + j);
            avg_rate += slip_rate(start + j);
            avg_stress += stress(start + j);
            avg_state += state(start + j);
            count++;
         }
         if (count > 0)
         {
            avg_slip /= count;
            avg_rate /= count;
            avg_stress /= count;
            avg_state /= count;
         }

         int e1 = fault_face_elem1_[i];
         int e2 = fault_face_elem2_[i];

         // L2 p=0: one DOF per element, DOF index = element index
         (*fault_slip_)(e1) = avg_slip;
         (*fault_slip_)(e2) = avg_slip;
         (*fault_slip_rate_)(e1) = avg_rate;
         (*fault_slip_rate_)(e2) = avg_rate;
         (*fault_stress_)(e1) = avg_stress / 1e6;  // Pa to MPa
         (*fault_stress_)(e2) = avg_stress / 1e6;
         (*fault_state_)(e1) = avg_state;
         (*fault_state_)(e2) = avg_state;
      }
   }

   /// @brief Write output if step count matches the interval.
   ///
   /// @param cycle Time step index
   /// @param time Physical time [s]
   /// @return true if data was written
   bool Save(int cycle, real_t time)
   {
      if (cycle % output_every_n_steps != 0) { return false; }

      pv_.SetCycle(cycle);
      pv_.SetTime(time);
      pv_.Save();

      last_write_time_ = time;
      return true;
   }

   /// @brief Write output based on adaptive time schedule (legacy).
   ///
   /// @param cycle Time step index
   /// @param time Physical time [s]
   /// @param V_max Maximum slip rate [m/s]
   /// @return true if data was written
   bool Save(int cycle, real_t time, real_t V_max)
   {
      real_t dt_out = OutputInterval(V_max);
      if (time - last_write_time_ < dt_out * 0.99) { return false; }

      return ForceSaveImpl(cycle, time);
   }

   /// @brief Force a save at the current state.
   void ForceSave(int cycle, real_t time)
   {
      ForceSaveImpl(cycle, time);
   }

   /// @brief Compute adaptive output interval for ParaView.
   static real_t OutputInterval(real_t V_max)
   {
      if (V_max > 1e-3)
      {
         return 0.01;  // Coseismic: every 0.01 s
      }
      else if (V_max > 1e-6)
      {
         return 1.0;   // Nucleation: every 1.0 s
      }
      else
      {
         return 1.0 * BP2Params::seconds_per_year;  // Interseismic: every 1 year
      }
   }

   /// Set data format (default: VTKFormat::BINARY)
   void SetDataFormat(VTKFormat fmt) { pv_.SetDataFormat(fmt); }

   /// Enable/disable high-order output (default: true)
   void SetHighOrderOutput(bool enable) { pv_.SetHighOrderOutput(enable); }

   /// @brief Get the number of fault face-element pairs mapped.
   int GetNumFaultFaces() const
   {
      return static_cast<int>(fault_face_elem1_.size());
   }

   /// @brief Check if fault output is initialized.
   bool HasFaultOutput() const { return has_fault_output_; }

private:
   MeshType &mesh_;
   int order_;
   real_t last_write_time_;
   ParaViewDataCollection pv_;

   // Fault output on domain mesh (L2 p=0)
   bool has_fault_output_ = false;
   std::unique_ptr<FiniteElementCollection> fault_fec_;
   std::unique_ptr<FiniteElementSpace> fault_fes_;
   std::unique_ptr<GridFunction> fault_slip_;
   std::unique_ptr<GridFunction> fault_slip_rate_;
   std::unique_ptr<GridFunction> fault_stress_;
   std::unique_ptr<GridFunction> fault_state_;

   // Fault face -> element mapping
   std::vector<int> fault_face_elem1_;
   std::vector<int> fault_face_elem2_;

   bool ForceSaveImpl(int cycle, real_t time)
   {
      pv_.SetCycle(cycle);
      pv_.SetTime(time);
      pv_.Save();

      last_write_time_ = time;
      return true;
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PARAVIEW_OUTPUT_HPP
