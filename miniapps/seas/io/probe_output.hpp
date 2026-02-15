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

#ifndef MFEM_SEAS_PROBE_OUTPUT_HPP
#define MFEM_SEAS_PROBE_OUTPUT_HPP

#include "mfem.hpp"
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

namespace mfem
{
namespace seas
{

/// @brief Manages time series output at a single probe location.
///
/// Writes columnar ASCII data to a file with a header comment.
/// Each call to WriteStep() appends one line.
class ProbeOutput
{
public:
   /// @brief Construct a probe output file.
   ///
   /// @param filename Output file path
   /// @param column_names Names for the header comment
   /// @param filename Output file path
   /// @param column_names Names for the header comment
   /// @param description Optional descriptive first line (e.g., probe location)
   ProbeOutput(const std::string &filename,
               const std::vector<std::string> &column_names,
               const std::string &description = "")
      : num_columns_(static_cast<int>(column_names.size()))
   {
      file_.open(filename, std::ios::out);
      MFEM_VERIFY(file_.is_open(), "Cannot open probe output file: " << filename);

      // Write optional description header (e.g., probe location)
      if (!description.empty())
      {
         file_ << "# " << description << "\n";
      }

      // Write column header
      file_ << "# Columns:";
      for (size_t i = 0; i < column_names.size(); i++)
      {
         file_ << " " << column_names[i];
         if (i + 1 < column_names.size()) { file_ << ","; }
      }
      file_ << "\n";
      file_ << std::scientific << std::uppercase << std::setprecision(12);
   }

   ~ProbeOutput()
   {
      if (file_.is_open()) { file_.close(); }
   }

   // Non-copyable
   ProbeOutput(const ProbeOutput &) = delete;
   ProbeOutput &operator=(const ProbeOutput &) = delete;

   // Movable
   ProbeOutput(ProbeOutput &&other) noexcept
      : file_(std::move(other.file_)),
        num_columns_(other.num_columns_) {}

   /// @brief Write a single time step row.
   ///
   /// @param values Vector of column values (size must match num_columns)
   void WriteStep(const std::vector<real_t> &values)
   {
      MFEM_ASSERT(static_cast<int>(values.size()) == num_columns_,
                  "WriteStep: expected " << num_columns_ << " columns, got "
                  << values.size());
      for (int i = 0; i < num_columns_; i++)
      {
         if (i > 0) { file_ << " "; }
         file_ << values[i];
      }
      file_ << "\n";
   }

   /// Flush output buffer to disk.
   void Flush()
   {
      if (file_.is_open()) { file_.flush(); }
   }

   /// Close the file.
   void Close()
   {
      if (file_.is_open()) { file_.close(); }
   }

   /// Check if the file is open and healthy.
   bool Good() const { return file_.is_open() && file_.good(); }

private:
   std::ofstream file_;
   int num_columns_;
};

/// @brief Interpolates fault data to probe locations.
///
/// Given a set of probe depths and fault DOF depths, precomputes
/// interpolation weights for efficient repeated evaluation.
/// Uses linear interpolation between the two nearest DOFs.
class ProbeInterpolator
{
public:
   /// @brief Construct the interpolator.
   ///
   /// @param fault_depths Z-coordinates of fault DOFs (from domain operator)
   /// @param probe_depths Z-coordinates of probe locations (negative = depth)
   ProbeInterpolator(const Vector &fault_depths,
                     const std::vector<real_t> &probe_depths)
      : num_probes_(static_cast<int>(probe_depths.size())),
        num_dofs_(fault_depths.Size())
   {
      idx_lo_.resize(num_probes_);
      idx_hi_.resize(num_probes_);
      weight_.resize(num_probes_);

      for (int p = 0; p < num_probes_; p++)
      {
         ComputeWeights(fault_depths, probe_depths[p], p);
      }
   }

   /// @brief Interpolate a scalar field to all probe locations.
   ///
   /// @param field Scalar field at fault DOFs [num_dofs_]
   /// @param probe_values Output values at probe locations [num_probes_]
   void Interpolate(const Vector &field, Vector &probe_values) const
   {
      MFEM_ASSERT(field.Size() == num_dofs_,
                  "Field size mismatch");
      probe_values.SetSize(num_probes_);

      for (int p = 0; p < num_probes_; p++)
      {
         int lo = idx_lo_[p];
         int hi = idx_hi_[p];
         real_t w = weight_[p];

         if (lo == hi || lo < 0 || hi < 0)
         {
            // Nearest-neighbor fallback
            int idx = (lo >= 0) ? lo : hi;
            idx = std::max(0, std::min(idx, num_dofs_ - 1));
            probe_values(p) = field(idx);
         }
         else
         {
            probe_values(p) = (1.0 - w) * field(lo) + w * field(hi);
         }
      }
   }

   /// @brief Interpolate a single value from a scalar field.
   ///
   /// @param field Scalar field at fault DOFs
   /// @param probe_idx Probe index
   /// @return Interpolated value
   real_t InterpolateOne(const Vector &field, int probe_idx) const
   {
      MFEM_ASSERT(probe_idx >= 0 && probe_idx < num_probes_,
                  "Probe index out of range");

      int lo = idx_lo_[probe_idx];
      int hi = idx_hi_[probe_idx];
      real_t w = weight_[probe_idx];

      if (lo == hi || lo < 0 || hi < 0)
      {
         int idx = (lo >= 0) ? lo : hi;
         idx = std::max(0, std::min(idx, num_dofs_ - 1));
         return field(idx);
      }
      return (1.0 - w) * field(lo) + w * field(hi);
   }

   /// Number of probe locations.
   int NumProbes() const { return num_probes_; }

private:
   int num_probes_;
   int num_dofs_;
   std::vector<int> idx_lo_;     ///< Lower bracket DOF index
   std::vector<int> idx_hi_;     ///< Upper bracket DOF index
   std::vector<real_t> weight_;  ///< Interpolation weight (0=lo, 1=hi)

   /// Compute interpolation weights for one probe.
   void ComputeWeights(const Vector &fault_depths, real_t z_probe, int p)
   {
      // Find the two DOFs bracketing z_probe
      int best_lo = -1, best_hi = -1;
      real_t z_lo = -std::numeric_limits<real_t>::max();
      real_t z_hi = std::numeric_limits<real_t>::max();
      int closest = 0;
      real_t min_dist = std::abs(fault_depths(0) - z_probe);

      for (int i = 0; i < num_dofs_; i++)
      {
         real_t zi = fault_depths(i);
         real_t dist = std::abs(zi - z_probe);
         if (dist < min_dist)
         {
            min_dist = dist;
            closest = i;
         }
         if (zi <= z_probe && zi > z_lo)
         {
            z_lo = zi;
            best_lo = i;
         }
         if (zi >= z_probe && zi < z_hi)
         {
            z_hi = zi;
            best_hi = i;
         }
      }

      if (best_lo >= 0 && best_hi >= 0 && best_lo != best_hi)
      {
         real_t dz = z_hi - z_lo;
         idx_lo_[p] = best_lo;
         idx_hi_[p] = best_hi;
         weight_[p] = (std::abs(dz) > 1e-30) ? (z_probe - z_lo) / dz : 0.0;
      }
      else
      {
         // Exact match or extrapolation: use closest
         idx_lo_[p] = closest;
         idx_hi_[p] = closest;
         weight_[p] = 0.0;
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PROBE_OUTPUT_HPP
