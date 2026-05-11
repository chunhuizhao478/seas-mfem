// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// data_field_3d.cpp — implementation of DataField3D (Phase 3 of the
// data-projection feature, plan v2).

#include "data_field_3d.hpp"

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>

namespace mfem
{
namespace seas
{

namespace
{

// ----------------------------------------------------------------------
// Local HDF5 helpers (C API).  Each helper that opens an HDF5 handle
// closes it before returning to keep the caller code straightforward
// even on error paths.
// ----------------------------------------------------------------------

std::string read_string_attr_(hid_t loc, const char* name)
{
   if (H5Aexists(loc, name) <= 0)
   {
      MFEM_ABORT("DataField3D: required string attribute '" << name
                 << "' is missing");
   }
   hid_t aid = H5Aopen(loc, name, H5P_DEFAULT);
   if (aid < 0)
   {
      MFEM_ABORT("DataField3D: failed to open attribute '" << name << "'");
   }
   hid_t atype = H5Aget_type(aid);
   if (atype < 0)
   {
      H5Aclose(aid);
      MFEM_ABORT("DataField3D: failed to get type of attribute '" << name
                 << "'");
   }
   H5T_class_t cls = H5Tget_class(atype);
   if (cls != H5T_STRING)
   {
      H5Tclose(atype);
      H5Aclose(aid);
      MFEM_ABORT("DataField3D: attribute '" << name << "' is not a string");
   }
   const bool is_variable_len = H5Tis_variable_str(atype) > 0;
   std::string out;
   if (is_variable_len)
   {
      // h5py defaults to variable-length strings.  Use H5Tcopy(atype)
      // for the read buffer so the encoding (ASCII vs UTF-8) matches
      // the on-disk attribute exactly; the original code used
      // H5Tcopy(H5T_C_S1) which fails the HDF5 conversion-path lookup
      // on h5py-written attrs that carry an explicit UTF-8 cset.
      hid_t mem_type = H5Tcopy(atype);
      char* tmp = nullptr;
      hid_t aspace = H5Aget_space(aid);
      herr_t status = H5Aread(aid, mem_type, &tmp);
      if (status < 0 || tmp == nullptr)
      {
         H5Sclose(aspace);
         H5Tclose(mem_type);
         H5Tclose(atype);
         H5Aclose(aid);
         MFEM_ABORT("DataField3D: failed to read variable-length string "
                    "attribute '" << name << "'");
      }
      out = tmp;
      H5Treclaim(mem_type, aspace, H5P_DEFAULT, &tmp);
      H5Sclose(aspace);
      H5Tclose(mem_type);
   }
   else
   {
      size_t size = H5Tget_size(atype);
      std::vector<char> buf(size + 1, '\0');
      herr_t status = H5Aread(aid, atype, buf.data());
      if (status < 0)
      {
         H5Tclose(atype);
         H5Aclose(aid);
         MFEM_ABORT("DataField3D: failed to read fixed-length string "
                    "attribute '" << name << "'");
      }
      out.assign(buf.data());
      // Strip trailing NULs.
      while (!out.empty() && out.back() == '\0')
      {
         out.pop_back();
      }
   }
   H5Tclose(atype);
   H5Aclose(aid);
   return out;
}

double read_double_attr_(hid_t loc, const char* name)
{
   if (H5Aexists(loc, name) <= 0)
   {
      MFEM_ABORT("DataField3D: required attribute '" << name
                 << "' is missing");
   }
   hid_t aid = H5Aopen(loc, name, H5P_DEFAULT);
   if (aid < 0)
   {
      MFEM_ABORT("DataField3D: failed to open attribute '" << name << "'");
   }
   double value = 0.0;
   herr_t status = H5Aread(aid, H5T_NATIVE_DOUBLE, &value);
   H5Aclose(aid);
   if (status < 0)
   {
      MFEM_ABORT("DataField3D: failed to read scalar attribute '" << name
                 << "'");
   }
   return value;
}

std::vector<double> read_axis_(hid_t file, const char* path,
                               const char* axis_label)
{
   hid_t did = H5Dopen2(file, path, H5P_DEFAULT);
   if (did < 0)
   {
      MFEM_ABORT("DataField3D: missing axis dataset '" << path << "'");
   }
   hid_t sid = H5Dget_space(did);
   if (sid < 0)
   {
      H5Dclose(did);
      MFEM_ABORT("DataField3D: failed to get dataspace for '" << path
                 << "'");
   }
   const int ndims = H5Sget_simple_extent_ndims(sid);
   if (ndims != 1)
   {
      H5Sclose(sid);
      H5Dclose(did);
      MFEM_ABORT("DataField3D: axis '" << axis_label << "' must be 1-D; "
                 "got ndims=" << ndims);
   }
   hsize_t dim = 0;
   H5Sget_simple_extent_dims(sid, &dim, nullptr);
   std::vector<double> data(static_cast<size_t>(dim));
   herr_t status = H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, data.data());
   H5Sclose(sid);
   H5Dclose(did);
   if (status < 0)
   {
      MFEM_ABORT("DataField3D: failed to read axis '" << path << "'");
   }
   if (data.size() < 1)
   {
      MFEM_ABORT("DataField3D: axis '" << axis_label << "' is empty");
   }
   for (size_t k = 1; k < data.size(); ++k)
   {
      if (!(data[k] > data[k - 1]))
      {
         MFEM_ABORT("DataField3D: axis '" << axis_label
                    << "' is not strictly monotone increasing at index "
                    << k << ": " << data[k - 1] << " -> " << data[k]);
      }
   }
   return data;
}

} // anonymous namespace


// ----------------------------------------------------------------------
// DataField3D ctor
// ----------------------------------------------------------------------

DataField3D::DataField3D(const std::string& sidecar_path,
                         const std::string& field_name,
                         OOBPolicy oob)
   : field_name_(field_name)
   , min_value_(0.0)
   , max_value_(0.0)
   , oob_policy_(oob)
{
   if (oob != OOBPolicy::Abort)
   {
      MFEM_ABORT("DataField3D: only OOBPolicy::Abort is supported in "
                 "schema-v1 (interpolation-only contract); got policy "
                 "value " << static_cast<int>(oob));
   }

   hid_t file = H5Fopen(sidecar_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0)
   {
      MFEM_ABORT("DataField3D: cannot open sidecar '" << sidecar_path << "'");
   }

   // Schema gates.
   const std::string schema = read_string_attr_(file, "schema_version");
   if (schema != "data_projection_v1")
   {
      H5Fclose(file);
      MFEM_ABORT("DataField3D: unexpected schema_version='" << schema
                 << "' in '" << sidecar_path
                 << "', expected 'data_projection_v1'");
   }
   crs_ = read_string_attr_(file, "crs");
   if (crs_ != "EPSG:32611")
   {
      H5Fclose(file);
      MFEM_ABORT("DataField3D: only crs='EPSG:32611' supported, got '"
                 << crs_ << "' in '" << sidecar_path << "'");
   }
   const std::string units_root = read_string_attr_(file, "units");
   if (units_root != "m")
   {
      H5Fclose(file);
      MFEM_ABORT("DataField3D: only units='m' supported, got '"
                 << units_root << "'");
   }
   const std::string z_positive = read_string_attr_(file, "z_positive");
   if (z_positive != "elevation")
   {
      H5Fclose(file);
      MFEM_ABORT("DataField3D: only z_positive='elevation' supported, "
                 "got '" << z_positive << "'");
   }

   // Axes.
   {
      auto x = read_axis_(file, "/grid/x", "x");
      auto y = read_axis_(file, "/grid/y", "y");
      auto z = read_axis_(file, "/grid/z", "z");
      x_.assign(x.begin(), x.end());
      y_.assign(y.begin(), y.end());
      z_.assign(z.begin(), z.end());
   }
   bbox_ = { x_.front(), x_.back(),
             y_.front(), y_.back(),
             z_.front(), z_.back() };

   // Field dataset.
   const std::string ds_path = std::string("/fields/") + field_name;
   hid_t did = H5Dopen2(file, ds_path.c_str(), H5P_DEFAULT);
   if (did < 0)
   {
      // Enumerate available fields for the abort message.
      std::ostringstream oss;
      oss << "DataField3D: field '" << field_name << "' not found in '"
          << sidecar_path << "'.";
      hid_t gid = H5Gopen2(file, "/fields", H5P_DEFAULT);
      if (gid >= 0)
      {
         oss << "  Available fields:";
         hsize_t nobj = 0;
         H5Gget_num_objs(gid, &nobj);
         for (hsize_t i = 0; i < nobj; ++i)
         {
            char name[256];
            ssize_t len = H5Gget_objname_by_idx(gid, i, name, sizeof(name));
            if (len > 0) { oss << ' ' << name; }
         }
         H5Gclose(gid);
      }
      H5Fclose(file);
      MFEM_ABORT(oss.str());
   }

   // Per-field attrs.
   units_ = read_string_attr_(did, "units");
   min_value_ = static_cast<real_t>(read_double_attr_(did, "min_value"));
   max_value_ = static_cast<real_t>(read_double_attr_(did, "max_value"));
   if (!(min_value_ < max_value_))
   {
      H5Dclose(did);
      H5Fclose(file);
      MFEM_ABORT("DataField3D: field '" << field_name
                 << "' has min_value (" << min_value_ << ") >= max_value ("
                 << max_value_ << ") — schema violation");
   }

   // Shape check.
   hid_t sid = H5Dget_space(did);
   const int ndims = H5Sget_simple_extent_ndims(sid);
   if (ndims != 3)
   {
      H5Sclose(sid);
      H5Dclose(did);
      H5Fclose(file);
      MFEM_ABORT("DataField3D: field '" << field_name
                 << "' must be 3-D; got ndims=" << ndims);
   }
   hsize_t dims[3] = {0, 0, 0};
   H5Sget_simple_extent_dims(sid, dims, nullptr);
   const size_t Nx = x_.size();
   const size_t Ny = y_.size();
   const size_t Nz = z_.size();
   if (dims[0] != Nx || dims[1] != Ny || dims[2] != Nz)
   {
      H5Sclose(sid);
      H5Dclose(did);
      H5Fclose(file);
      MFEM_ABORT("DataField3D: field '" << field_name << "' shape ["
                 << dims[0] << ", " << dims[1] << ", " << dims[2]
                 << "] does not match grid (" << Nx << ", " << Ny
                 << ", " << Nz << ")");
   }
   H5Sclose(sid);

   // Read data.
   data_.assign(Nx * Ny * Nz, 0.0);
   {
      std::vector<double> tmp(Nx * Ny * Nz);
      herr_t status = H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                              H5P_DEFAULT, tmp.data());
      if (status < 0)
      {
         H5Dclose(did);
         H5Fclose(file);
         MFEM_ABORT("DataField3D: failed to read field '" << field_name
                    << "' data");
      }
      for (size_t i = 0; i < tmp.size(); ++i)
      {
         data_[i] = static_cast<real_t>(tmp[i]);
      }
   }

   H5Dclose(did);
   H5Fclose(file);

   // Post-load value sanity check (defense-in-depth vs. writer guards).
   for (size_t i = 0; i < data_.size(); ++i)
   {
      const real_t v = data_[i];
      if (std::isnan(v))
      {
         // Decode i back to (i_x, j_y, k_z) for diagnostics.
         const size_t k = i % Nz;
         const size_t j = (i / Nz) % Ny;
         const size_t ii = i / (Ny * Nz);
         MFEM_ABORT("DataField3D: NaN in field '" << field_name
                    << "' at flat index " << i << " (i, j, k) = (" << ii
                    << ", " << j << ", " << k << "). v1 schema forbids "
                    "NaN; regenerate the sidecar.");
      }
      if (v < min_value_ || v > max_value_)
      {
         const size_t k = i % Nz;
         const size_t j = (i / Nz) % Ny;
         const size_t ii = i / (Ny * Nz);
         MFEM_ABORT("DataField3D: field '" << field_name
                    << "' value " << v << " out of declared range ["
                    << min_value_ << ", " << max_value_
                    << "] at (i, j, k) = (" << ii << ", " << j << ", "
                    << k << "). Regenerate the source dataset.");
      }
   }
}


// ----------------------------------------------------------------------
// Public methods
// ----------------------------------------------------------------------

bool DataField3D::ContainsBBox(real_t xmin, real_t xmax,
                               real_t ymin, real_t ymax,
                               real_t zmin, real_t zmax,
                               real_t eps) const
{
   return (xmin >= bbox_[0] - eps) && (xmax <= bbox_[1] + eps)
          && (ymin >= bbox_[2] - eps) && (ymax <= bbox_[3] + eps)
          && (zmin >= bbox_[4] - eps) && (zmax <= bbox_[5] + eps);
}

real_t DataField3D::Evaluate(real_t x, real_t y, real_t z) const
{
   switch (interp_mode_)
   {
      case InterpMode::Trilinear:
         return evaluate_trilinear_(x, y, z);
      case InterpMode::CatmullRom:
         return evaluate_catmull_rom_(x, y, z);
   }
   // Should be unreachable; the enum has only the two values above.
   MFEM_ABORT("DataField3D::Evaluate: unknown InterpMode value "
              << static_cast<int>(interp_mode_));
   return real_t{0};
}

real_t DataField3D::evaluate_trilinear_(real_t x, real_t y, real_t z) const
{
   const int i = find_index_(x_, x);
   const int j = find_index_(y_, y);
   const int k = find_index_(z_, z);

   const real_t x0 = x_[i],   x1 = x_[i + 1];
   const real_t y0 = y_[j],   y1 = y_[j + 1];
   const real_t z0 = z_[k],   z1 = z_[k + 1];

   const real_t u = (x - x0) / (x1 - x0);
   const real_t v = (y - y0) / (y1 - y0);
   const real_t w = (z - z0) / (z1 - z0);

   const real_t f000 = data_[flat_index_(i,   j,   k  )];
   const real_t f100 = data_[flat_index_(i+1, j,   k  )];
   const real_t f010 = data_[flat_index_(i,   j+1, k  )];
   const real_t f110 = data_[flat_index_(i+1, j+1, k  )];
   const real_t f001 = data_[flat_index_(i,   j,   k+1)];
   const real_t f101 = data_[flat_index_(i+1, j,   k+1)];
   const real_t f011 = data_[flat_index_(i,   j+1, k+1)];
   const real_t f111 = data_[flat_index_(i+1, j+1, k+1)];

   const real_t one_u = 1.0 - u;
   const real_t one_v = 1.0 - v;
   const real_t one_w = 1.0 - w;

   return (one_u * one_v * one_w) * f000
        + (    u * one_v * one_w) * f100
        + (one_u *     v * one_w) * f010
        + (    u *     v * one_w) * f110
        + (one_u * one_v *     w) * f001
        + (    u * one_v *     w) * f101
        + (one_u *     v *     w) * f011
        + (    u *     v *     w) * f111;
}


// ----------------------------------------------------------------------
// Catmull-Rom tricubic
// ----------------------------------------------------------------------
//
// 1-D Catmull-Rom basis on the unit interval t in [0, 1] interpolating
// between p1 (at t=0) and p2 (at t=1), with p0 / p3 as the surrounding
// "tangent" samples (one-voxel-out on either side).  Standard cardinal-
// cubic blending function (same form as Bartels / Beatty / Barsky):
//
//   b0(t) = -0.5 t + t^2 - 0.5 t^3       (weight on p0 = the i-1 sample)
//   b1(t) =  1   - 2.5 t^2 + 1.5 t^3     (weight on p1 = the i   sample)
//   b2(t) =  0.5 t + 2 t^2 - 1.5 t^3     (weight on p2 = the i+1 sample)
//   b3(t) = -0.5 t^2 + 0.5 t^3           (weight on p3 = the i+2 sample)
//
// At t == 0 (the i sample) the basis collapses to {0, 1, 0, 0}; at
// t == 1 to {0, 0, 1, 0}.  Hence the spline is interpolating
// (passes through every voxel value).  C^1 continuity follows from
// matching b'_(0..3)(0) and b'_(0..3)(1) at adjacent intervals.
//
// The 3-D evaluator does the tensor-product reduction in three stages:
//   (1) Along z: 16 cubic interpolations -> 16 z-row values
//   (2) Along y:  4 cubic interpolations -> 4 (y, z=t_z) values
//   (3) Along x:  1 cubic interpolation  -> the final scalar
// Total: 21 cubic 1-D evaluations, 16 + 4 + 1 = 21 calls to a 4-tap
// inner product.  The 64 source voxel values are each touched exactly
// once.
//
// Edge handling: when the query cell at (i, j, k) on any axis does
// not have BOTH a (-1) and a (+2) neighbour (i.e. cell at the leading
// or trailing axis edge), we fall back to the ``evaluate_trilinear_``
// path for the whole query rather than asymmetric clamping or
// reflection — bounded fallback preserves the schema-v1 invariant
// that the result is always an interpolation of real sidecar values.
// ----------------------------------------------------------------------

namespace
{

// 1-D Catmull-Rom: 4 weights at parameter ``t`` (in [0, 1]).
inline void catmull_rom_weights_(real_t t,
                                 real_t& b0, real_t& b1,
                                 real_t& b2, real_t& b3)
{
   const real_t t2 = t * t;
   const real_t t3 = t2 * t;
   b0 = real_t{-0.5} * t + t2 - real_t{0.5} * t3;
   b1 = real_t{1.0} - real_t{2.5} * t2 + real_t{1.5} * t3;
   b2 = real_t{0.5} * t + real_t{2.0} * t2 - real_t{1.5} * t3;
   b3 = real_t{-0.5} * t2 + real_t{0.5} * t3;
}

inline real_t cubic_dot_(const real_t b0, const real_t b1,
                         const real_t b2, const real_t b3,
                         const real_t p0, const real_t p1,
                         const real_t p2, const real_t p3)
{
   return b0 * p0 + b1 * p1 + b2 * p2 + b3 * p3;
}

} // anonymous namespace


real_t DataField3D::evaluate_catmull_rom_(real_t x, real_t y, real_t z) const
{
   const int i = find_index_(x_, x);
   const int j = find_index_(y_, y);
   const int k = find_index_(z_, z);

   const int Nx = static_cast<int>(x_.size());
   const int Ny = static_cast<int>(y_.size());
   const int Nz = static_cast<int>(z_.size());

   // Bounded fallback: we need one extra neighbour on either side of
   // the [i, i+1] cell on every axis.  If any axis is at its leading
   // or trailing cell, defer to trilinear (still correct, no new
   // schema-v1 violation).
   if (i < 1 || i > Nx - 3 ||
       j < 1 || j > Ny - 3 ||
       k < 1 || k > Nz - 3)
   {
      return evaluate_trilinear_(x, y, z);
   }

   // Local Catmull-Rom requires UNIFORM spacing along each axis to
   // collapse to the standard cardinal-cubic basis above.  For
   // non-uniform spacing the basis would need tangents
   // (p2 - p0) / (axis[i+1] - axis[i-1]) etc. — schema-v1 sidecars
   // produced by the canonical CVM-H pipeline are uniform in x and y
   // and PIECEWISE-uniform in z; for cells where the local spacing
   // departs from uniform by > 1 % on any axis, also fall back to
   // trilinear for safety.
   const real_t dxm1 = x_[i]     - x_[i - 1];
   const real_t dx0  = x_[i + 1] - x_[i];
   const real_t dxp1 = x_[i + 2] - x_[i + 1];
   const real_t dym1 = y_[j]     - y_[j - 1];
   const real_t dy0  = y_[j + 1] - y_[j];
   const real_t dyp1 = y_[j + 2] - y_[j + 1];
   const real_t dzm1 = z_[k]     - z_[k - 1];
   const real_t dz0  = z_[k + 1] - z_[k];
   const real_t dzp1 = z_[k + 2] - z_[k + 1];
   const real_t tol = real_t{0.01};
   auto nearly_uniform = [tol](real_t a, real_t b, real_t c)
   {
      const real_t mid = b;
      const real_t lo = std::min({a, b, c});
      const real_t hi = std::max({a, b, c});
      return (hi - lo) <= tol * mid;
   };
   if (!nearly_uniform(dxm1, dx0, dxp1) ||
       !nearly_uniform(dym1, dy0, dyp1) ||
       !nearly_uniform(dzm1, dz0, dzp1))
   {
      return evaluate_trilinear_(x, y, z);
   }

   const real_t tx = (x - x_[i]) / dx0;
   const real_t ty = (y - y_[j]) / dy0;
   const real_t tz = (z - z_[k]) / dz0;

   real_t bx0, bx1, bx2, bx3;
   real_t by0, by1, by2, by3;
   real_t bz0, bz1, bz2, bz3;
   catmull_rom_weights_(tx, bx0, bx1, bx2, bx3);
   catmull_rom_weights_(ty, by0, by1, by2, by3);
   catmull_rom_weights_(tz, bz0, bz1, bz2, bz3);

   // Stage 1: along z.  For each of the 16 (i_off, j_off) pairs
   // collapse the 4 z-samples into a single value.
   real_t plane_yz[4][4];
   for (int ix = 0; ix < 4; ++ix)
   {
      const int ii = i - 1 + ix;
      for (int iy = 0; iy < 4; ++iy)
      {
         const int jj = j - 1 + iy;
         const real_t pz0 = data_[flat_index_(ii, jj, k - 1)];
         const real_t pz1 = data_[flat_index_(ii, jj, k    )];
         const real_t pz2 = data_[flat_index_(ii, jj, k + 1)];
         const real_t pz3 = data_[flat_index_(ii, jj, k + 2)];
         plane_yz[ix][iy] =
            cubic_dot_(bz0, bz1, bz2, bz3, pz0, pz1, pz2, pz3);
      }
   }

   // Stage 2: along y.  Collapse each (ix, *) row of 4 y-samples.
   real_t row_x[4];
   for (int ix = 0; ix < 4; ++ix)
   {
      row_x[ix] = cubic_dot_(by0, by1, by2, by3,
                             plane_yz[ix][0], plane_yz[ix][1],
                             plane_yz[ix][2], plane_yz[ix][3]);
   }

   // Stage 3: along x.  Final cubic.
   return cubic_dot_(bx0, bx1, bx2, bx3,
                     row_x[0], row_x[1], row_x[2], row_x[3]);
}


// ----------------------------------------------------------------------
// Private
// ----------------------------------------------------------------------

int DataField3D::find_index_(const std::vector<real_t>& axis,
                             real_t v) const
{
   const int N = static_cast<int>(axis.size());
   if (N < 2)
   {
      MFEM_ABORT("DataField3D::find_index_: axis has fewer than 2 "
                 "samples (size=" << N << "); a single-cell axis cannot "
                 "support trilinear evaluation");
   }
   if (v < axis.front() || v > axis.back())
   {
      MFEM_ABORT("DataField3D::Evaluate: query value " << v
                 << " is OUTSIDE field '" << field_name_
                 << "' axis bbox [" << axis.front() << ", " << axis.back()
                 << "]. v1 schema enforces interpolation-only; either "
                 "shrink the mesh or regenerate the sidecar.");
   }
   // upper_bound returns first element > v.  i = (it - begin) - 1 gives
   // the leftmost index of the cell containing v, except at v == axis[0]
   // where upper_bound returns axis.begin()+1 already → it-begin = 1 →
   // i = 0.  At v == axis[N-1], upper_bound returns axis.end() → i = N-1
   // which would index axis[N] (overrun); clamp to N-2.
   auto it = std::upper_bound(axis.begin(), axis.end(), v);
   int i = static_cast<int>(it - axis.begin()) - 1;
   if (i < 0)            { i = 0; }
   if (i >= N - 1)       { i = N - 2; }
   return i;
}

} // namespace seas
} // namespace mfem
