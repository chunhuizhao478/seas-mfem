// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// grid_field_2d.hpp — 2-D gridded input fields for the 2.5D cycle
// (heterogeneity plan rev4, Step B1; fully-dynamic-fault-opening-hetergeneity).
//
// Reads named fields of a "seas_grid2d_v1" NetCDF-4 file through the HDF5 C
// API (a NetCDF-4 file is an HDF5 file; decision D7 = a, no NetCDF library):
//
//   datasets  x (Nx), y (Ny)          float64, strictly ascending, metres in the
//                                      mesh frame (the mesh's own coordinates)
//             <field> (Ny, Nx)         float64, row-major (y slowest), no NaN
//   root attribute  schema = "seas_grid2d_v1"
//
// The writer is miniapps/seas/ucerf2d/hetero/build_hetero_grid.py.
//
// Evaluate(x, y) is bilinear in the cell [x_i, x_{i+1}] x [y_j, y_{j+1}]
// that holds the point, written in the nested lerp form
//
//    f_0 = f(i, j)   + tx (f(i+1, j)   - f(i, j))
//    f_1 = f(i, j+1) + tx (f(i+1, j+1) - f(i, j+1))
//    f   = f_0       + ty (f_1 - f_0)            (t = 1 returns the far value)
//
// so that a constant field returns its value bit for bit (gate U4 of the
// plan: a uniform file reproduces the constant-material path exactly), a
// grid node returns its stored value exactly, and a linear field is
// reproduced to round-off (gate U1).
//
// Outside the grid: OOBPolicy::Abort (default) aborts with the point and the
// grid's box; OOBPolicy::Clamp clamps the point to the box first.
//
// Every rank reads the whole file (the UCERF file is ~1.6 M nodes x 7 fields).

#ifndef MFEM_SEAS_GRID_FIELD_2D_HPP
#define MFEM_SEAS_GRID_FIELD_2D_HPP

#include "mfem.hpp"
#include "data_field_3d.hpp"   // OOBPolicy

#include <hdf5.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

class GridField2D
{
public:
   /// Load @a names from @a path.  Aborts on a missing file, a missing or
   /// wrong "schema" attribute, a missing dataset, a non-ascending or
   /// too-short axis, a field of the wrong shape, or a NaN/Inf value.
   GridField2D(const std::string &path, const std::vector<std::string> &names,
               OOBPolicy oob = OOBPolicy::Abort)
      : path_(path), names_(names), oob_(oob)
   {
      MFEM_VERIFY(!names_.empty(), "GridField2D: no field names for " << path_);
      hid_t f = H5Fopen(path_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      MFEM_VERIFY(f >= 0, "GridField2D: cannot open " << path_);
      const std::string schema = ReadStringAttr_(f, "schema");
      if (schema != "seas_grid2d_v1")
      {
         H5Fclose(f);
         MFEM_ABORT("GridField2D: " << path_ << " has schema \"" << schema
                    << "\", expected \"seas_grid2d_v1\"");
      }
      x_ = Read1D_(f, "x");
      y_ = Read1D_(f, "y");
      CheckAxis_(x_, "x");
      CheckAxis_(y_, "y");
      data_.resize(names_.size());
      for (size_t k = 0; k < names_.size(); k++) { data_[k] = Read2D_(f, names_[k]); }
      H5Fclose(f);
   }

   int NumFields() const { return static_cast<int>(names_.size()); }
   const std::vector<std::string> &Names() const { return names_; }
   int Index(const std::string &name) const
   {
      for (size_t k = 0; k < names_.size(); k++) { if (names_[k] == name) { return static_cast<int>(k); } }
      MFEM_ABORT("GridField2D: no field \"" << name << "\" loaded from " << path_);
      return -1;
   }
   int NumX() const { return static_cast<int>(x_.size()); }
   int NumY() const { return static_cast<int>(y_.size()); }
   /// {xmin, xmax, ymin, ymax}
   std::array<real_t, 4> BBox() const { return {x_.front(), x_.back(), y_.front(), y_.back()}; }
   bool Contains(real_t x, real_t y) const
   { return x >= x_.front() && x <= x_.back() && y >= y_.front() && y <= y_.back(); }
   /// The stored value at node (i, j) (i along x).
   real_t Node(int field, int i, int j) const { return data_[field][static_cast<size_t>(j) * x_.size() + i]; }

   /// Bilinear value of field @a field at (x, y) (file comment).
   real_t Evaluate(int field, real_t x, real_t y) const
   {
      int i, j; real_t tx, ty;
      Locate_(x, y, i, j, tx, ty);
      return Lerp2_(data_[field], i, j, tx, ty);
   }
   /// All loaded fields at (x, y), in load order, into @a out (size NumFields()).
   void EvaluateAll(real_t x, real_t y, real_t *out) const
   {
      int i, j; real_t tx, ty;
      Locate_(x, y, i, j, tx, ty);
      for (size_t k = 0; k < data_.size(); k++) { out[k] = Lerp2_(data_[k], i, j, tx, ty); }
   }

private:
   std::string path_;
   std::vector<std::string> names_;
   OOBPolicy oob_;
   std::vector<real_t> x_, y_;
   std::vector<std::vector<real_t>> data_;

   real_t Lerp2_(const std::vector<real_t> &d, int i, int j, real_t tx, real_t ty) const
   {
      const size_t nx = x_.size();
      const real_t f00 = d[j * nx + i], f10 = d[j * nx + i + 1];
      const real_t f01 = d[(j + 1) * nx + i], f11 = d[(j + 1) * nx + i + 1];
      const real_t f0 = Lerp_(f00, f10, tx);
      const real_t f1 = Lerp_(f01, f11, tx);
      return Lerp_(f0, f1, ty);
   }

   /// a + t (b - a), with t == 1 returning b exactly (a point on the last
   /// grid line of an axis sits in the last cell with t = 1).
   static real_t Lerp_(real_t a, real_t b, real_t t) { return t == 1.0 ? b : a + t * (b - a); }

   void Locate_(real_t x, real_t y, int &i, int &j, real_t &tx, real_t &ty) const
   {
      if (!Contains(x, y))
      {
         if (oob_ == OOBPolicy::Clamp)
         {
            x = std::min(std::max(x, x_.front()), x_.back());
            y = std::min(std::max(y, y_.front()), y_.back());
         }
         else
         {
            MFEM_ABORT("GridField2D: point (" << x << ", " << y << ") is outside the grid of "
                       << path_ << ": x [" << x_.front() << ", " << x_.back() << "], y ["
                       << y_.front() << ", " << y_.back() << "] (outside = \"abort\")");
         }
      }
      i = Interval_(x_, x); j = Interval_(y_, y);
      tx = (x - x_[i]) / (x_[i + 1] - x_[i]);
      ty = (y - y_[j]) / (y_[j + 1] - y_[j]);
   }

   /// i with a[i] <= v <= a[i+1]; v on the last node maps to the last cell.
   static int Interval_(const std::vector<real_t> &a, real_t v)
   {
      const int n = static_cast<int>(a.size());
      int i = static_cast<int>(std::upper_bound(a.begin(), a.end(), v) - a.begin()) - 1;
      return std::min(std::max(i, 0), n - 2);
   }

   void CheckAxis_(const std::vector<real_t> &a, const char *name) const
   {
      MFEM_VERIFY(a.size() >= 2, "GridField2D: axis " << name << " of " << path_
                  << " has " << a.size() << " values (need >= 2)");
      for (size_t i = 0; i + 1 < a.size(); i++)
      {
         MFEM_VERIFY(std::isfinite(a[i]) && a[i + 1] > a[i], "GridField2D: axis " << name
                     << " of " << path_ << " is not strictly ascending at index " << i);
      }
   }

   std::string ReadStringAttr_(hid_t f, const char *name) const
   {
      MFEM_VERIFY(H5Aexists(f, name) > 0, "GridField2D: " << path_
                  << " has no root attribute \"" << name << "\"");
      hid_t a = H5Aopen(f, name, H5P_DEFAULT);
      hid_t t = H5Aget_type(a);
      std::string s;
      if (H5Tget_class(t) == H5T_STRING)
      {
         if (H5Tis_variable_str(t) > 0)
         {
            char *buf = nullptr;
            hid_t mt = H5Tcopy(H5T_C_S1); H5Tset_size(mt, H5T_VARIABLE);
            H5Aread(a, mt, &buf);
            if (buf) { s = buf; H5free_memory(buf); }
            H5Tclose(mt);
         }
         else
         {
            const size_t n = H5Tget_size(t);
            std::vector<char> buf(n + 1, '\0');
            H5Aread(a, t, buf.data());
            s.assign(buf.data());
         }
      }
      H5Tclose(t); H5Aclose(a);
      return s;
   }

   std::vector<real_t> Read1D_(hid_t f, const char *name) const
   {
      hid_t d = H5Dopen2(f, name, H5P_DEFAULT);
      MFEM_VERIFY(d >= 0, "GridField2D: " << path_ << " has no dataset \"" << name << "\"");
      hid_t sp = H5Dget_space(d);
      MFEM_VERIFY(H5Sget_simple_extent_ndims(sp) == 1, "GridField2D: dataset \"" << name
                  << "\" of " << path_ << " is not 1-D");
      hsize_t n = 0; H5Sget_simple_extent_dims(sp, &n, nullptr);
      std::vector<double> v(n);
      MFEM_VERIFY(H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >= 0,
                  "GridField2D: read of \"" << name << "\" from " << path_ << " failed");
      H5Sclose(sp); H5Dclose(d);
      return std::vector<real_t>(v.begin(), v.end());
   }

   std::vector<real_t> Read2D_(hid_t f, const std::string &name) const
   {
      hid_t d = H5Dopen2(f, name.c_str(), H5P_DEFAULT);
      MFEM_VERIFY(d >= 0, "GridField2D: " << path_ << " has no field \"" << name << "\"");
      hid_t sp = H5Dget_space(d);
      MFEM_VERIFY(H5Sget_simple_extent_ndims(sp) == 2, "GridField2D: field \"" << name
                  << "\" of " << path_ << " is not 2-D");
      hsize_t dims[2]; H5Sget_simple_extent_dims(sp, dims, nullptr);
      MFEM_VERIFY(dims[0] == y_.size() && dims[1] == x_.size(), "GridField2D: field \""
                  << name << "\" of " << path_ << " is " << dims[0] << " x " << dims[1]
                  << ", expected (Ny, Nx) = " << y_.size() << " x " << x_.size());
      std::vector<double> v(dims[0] * dims[1]);
      MFEM_VERIFY(H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >= 0,
                  "GridField2D: read of \"" << name << "\" from " << path_ << " failed");
      H5Sclose(sp); H5Dclose(d);
      for (size_t k = 0; k < v.size(); k++)
      {
         MFEM_VERIFY(std::isfinite(v[k]), "GridField2D: field \"" << name << "\" of " << path_
                     << " has a non-finite value at flat index " << k);
      }
      return std::vector<real_t>(v.begin(), v.end());
   }
};

}  // namespace seas
}  // namespace mfem

#endif // MFEM_SEAS_GRID_FIELD_2D_HPP
