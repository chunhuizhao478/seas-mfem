// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// snapshots.hpp — fully dynamic BP5 cycle: whole-field snapshots of the
// fault and of the free surface ([cycle.snapshots], cycle_config.hpp), one
// file per kind for the whole run:
//   <dir>/<prefix>_fault_snapshots.vtkhdf     on every owned fault DOF:
//       slip_strike, slip_dip [m], slip_rate_strike, slip_rate_dip [m/s],
//       state (theta) [s], normal_stress (sigma_n_eff, compression
//       positive) [MPa], shear_stress_strike, shear_stress_dip [MPa]
//   <dir>/<prefix>_surface_snapshots.vtkhdf   on the free surface:
//       displacement, velocity (global XYZ) [m, m/s]
// The fault values are the station files' at each DOF, term by term
// (fd_cycle_output.hpp): slip and V in kernel components, shear stress
// -(tau_a + tau_b + eta V) / 1e6, theta, from the regime's FillStationFields.
// The surface values are the displacement field evaluated at the points of
// each selected boundary face (its vertices at p = 1; p > 1: the face
// refined p times).  The velocity is U' in the dynamic regime; in the QD
// regime, which has no velocity, it is the backward difference of the
// displacement over the last accepted QD step (0 at a run's first snapshot,
// and from the switch-back balance after a dynamic phase).
//
// Format: VTKHDF 2.0, a transient UnstructuredGrid whose mesh is written
// once (ParaView 5.12 or later reads it; a time step per snapshot).  The
// fault grid holds one triangle per fault face at p = 1 (the face's three
// DOFs) and one vertex cell per DOF at p > 1.  Points are per face (the
// fault fields are discontinuous), in rank order then owned order, so a
// restart needs the same rank count (the checkpoint needs it too).  Besides
// /VTKHDF, the group /fd_cycle holds per snapshot: regime (0 QD, 1
// dynamic), phase (dynamic phases entered so far), reason (SnapReason) and
// step (QD TS step or dynamic step in its phase), and the file's
// description, commit and units as attributes.
//
// The root gathers each snapshot (MPI_Gatherv) and writes it with the HDF5
// C API (serial calls on the root; the HDF5 1.10 API of the cloud build).
// Every snapshot ends with H5Fflush.  A job killed while a snapshot is being
// written can still leave the file unreadable; the station files are text.
//
// Schedule (every rank decides the same way, from t and tau alone):
//   QD:       mark k = floor(t / dt_qd); a snapshot when k > the last QD
//             mark written (so t = 0 is one), or at a switch to dynamic
//   dynamic:  mark k = floor(tau / dt_dyn), counted from the phase start
//             (whose time already has the switch snapshot)
//   one at the end of each dynamic phase and at the end of the run; when a
//   snapshot already has that time, it takes that reason instead.
// Restart: the checkpoint keeps each file's snapshot count, last QD mark and
// last time (GetState); a run restarted from it trims the file to that
// count and continues.  A checkpoint written with snapshots off cannot
// continue a snapshot file: the restart refuses.
//
// No solver header (rule M4): fault fields arrive through the regime
// adapters' FillStationFields, the displacement as a true-DOF vector of the
// QD regime's DG space.

#ifndef MFEM_SEAS_CYCLE_SNAPSHOTS_HPP
#define MFEM_SEAS_CYCLE_SNAPSHOTS_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "fault_fields.hpp"
#include "fd_cycle_output.hpp"

#include <hdf5.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// Why a snapshot was written (the /fd_cycle/reason dataset).
enum class SnapReason : int
{
   Cadence = 0,        ///< a dt_qd or dt_dyn mark
   Start = 1,          ///< t = 0 of a run that starts dynamic
   ToDynamic = 2,      ///< the QD step at which the run switches to dynamic
   EndOfDynamic = 3,   ///< the end of a dynamic phase
   End = 4             ///< the end of the run
};

/// A transient VTKHDF UnstructuredGrid with a fixed mesh and appended point
/// data, written by one process (the root).
class VtkhdfSeries
{
public:
   struct Field
   {
      std::string name;
      int ncomp;   ///< 1 (scalar) or 3 (vector)
   };

   /// A new file (an existing one is replaced) with the grid and no step.
   VtkhdfSeries(const std::string &path, const std::vector<double> &points,
                const std::vector<long long> &connectivity,
                const std::vector<long long> &offsets, const std::vector<unsigned char> &types,
                const std::vector<Field> &fields,
                const std::vector<std::pair<std::string, std::string>> &attrs)
      : path_(path), fields_(fields)
   {
      npts_ = static_cast<long long>(points.size() / 3);
      MFEM_VERIFY(npts_ > 0 && static_cast<long long>(points.size()) == 3 * npts_,
                  "VtkhdfSeries: " << points.size() << " point coordinates for " << path);
      const long long ncells = static_cast<long long>(types.size());
      MFEM_VERIFY(static_cast<long long>(offsets.size()) == ncells + 1 &&
                  offsets.back() == static_cast<long long>(connectivity.size()),
                  "VtkhdfSeries: offsets do not match the connectivity for " << path);
      file_ = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      Check(file_, "H5Fcreate");
      vtk_ = H5Gcreate2(file_, "VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      Check(vtk_, "H5Gcreate2 VTKHDF");
      const long long version[2] = {2, 0};
      AttrInt64(vtk_, "Version", version, 2);
      AttrString(vtk_, "Type", "UnstructuredGrid");

      const long long np = npts_, nc = ncells;
      const long long nconn = static_cast<long long>(connectivity.size());
      Fixed(vtk_, "NumberOfPoints", H5T_NATIVE_LLONG, H5T_STD_I64LE, {1}, &np);
      Fixed(vtk_, "NumberOfCells", H5T_NATIVE_LLONG, H5T_STD_I64LE, {1}, &nc);
      Fixed(vtk_, "NumberOfConnectivityIds", H5T_NATIVE_LLONG, H5T_STD_I64LE, {1}, &nconn);
      Fixed(vtk_, "Points", H5T_NATIVE_DOUBLE, H5T_IEEE_F64LE,
            {static_cast<hsize_t>(npts_), 3}, points.data());
      Fixed(vtk_, "Connectivity", H5T_NATIVE_LLONG, H5T_STD_I64LE,
            {static_cast<hsize_t>(nconn)}, connectivity.data());
      Fixed(vtk_, "Offsets", H5T_NATIVE_LLONG, H5T_STD_I64LE,
            {static_cast<hsize_t>(ncells + 1)}, offsets.data());
      Fixed(vtk_, "Types", H5T_NATIVE_UCHAR, H5T_STD_U8LE, {static_cast<hsize_t>(ncells)},
            types.data());

      const hid_t pd = H5Gcreate2(vtk_, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      Check(pd, "H5Gcreate2 PointData");
      for (const Field &f : fields_)
      {
         Close(Growing(pd, f.name, H5T_IEEE_F64LE, f.ncomp, static_cast<hsize_t>(npts_)));
      }
      Close(pd, true);
      const hid_t st = H5Gcreate2(vtk_, "Steps", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      Check(st, "H5Gcreate2 Steps");
      const long long zero = 0;
      AttrInt64(st, "NSteps", &zero, 1);
      Close(Growing(st, "Values", H5T_IEEE_F64LE, 1, 1024));
      for (const char *n : kStepOffsets) { Close(Growing(st, n, H5T_STD_I64LE, 1, 1024)); }
      const hid_t pdo = H5Gcreate2(st, "PointDataOffsets", H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT);
      Check(pdo, "H5Gcreate2 Steps/PointDataOffsets");
      for (const Field &f : fields_) { Close(Growing(pdo, f.name, H5T_STD_I64LE, 1, 1024)); }
      Close(pdo, true);
      Close(st, true);

      const hid_t fc = H5Gcreate2(file_, "fd_cycle", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      Check(fc, "H5Gcreate2 fd_cycle");
      for (const auto &a : attrs) { AttrString(fc, a.first.c_str(), a.second); }
      for (const char *n : kMeta) { Close(Growing(fc, n, H5T_STD_I64LE, 1, 1024)); }
      Close(fc, true);
      Flush();
   }

   /// An existing file from a restart, trimmed to its first `keep` steps;
   /// `points` must be the file's (same rank count and partition).
   VtkhdfSeries(const std::string &path, const std::vector<double> &points,
                const std::vector<Field> &fields, long long keep)
      : path_(path), fields_(fields)
   {
      file_ = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      Check(file_, "H5Fopen (restart)");
      vtk_ = H5Gopen2(file_, "VTKHDF", H5P_DEFAULT);
      Check(vtk_, "H5Gopen2 VTKHDF");
      long long np = 0;
      ReadFixed(vtk_, "NumberOfPoints", H5T_NATIVE_LLONG, &np, 1);
      npts_ = np;
      MFEM_VERIFY(3 * npts_ == static_cast<long long>(points.size()),
                  "snapshots restart: " << path << " has " << npts_ << " points, this run "
                  << points.size() / 3 << " (same rank count and order needed)");
      std::vector<double> have(points.size());
      ReadFixed(vtk_, "Points", H5T_NATIVE_DOUBLE, have.data(), have.size());
      MFEM_VERIFY(have == points, "snapshots restart: the points of " << path << " differ "
                  "from this run's (another partition or mesh)");
      nsteps_ = Rows("Steps/Values");
      MFEM_VERIFY(nsteps_ >= keep, "snapshots restart: " << path << " has " << nsteps_
                  << " snapshots, the checkpoint recorded " << keep);
      for (const Field &f : fields_)
      {
         MFEM_VERIFY(Rows("PointData/" + f.name) >= keep * npts_,
                     "snapshots restart: " << path << ": PointData/" << f.name
                     << " is shorter than " << keep << " snapshots");
         Trim("PointData/" + f.name, keep * npts_);
         Trim("Steps/PointDataOffsets/" + f.name, keep);
      }
      Trim("Steps/Values", keep);
      for (const char *n : kStepOffsets) { Trim(std::string("Steps/") + n, keep); }
      for (const char *n : kMeta) { Trim(std::string("/fd_cycle/") + n, keep); }
      nsteps_ = keep;
      SetNSteps();
      Flush();
   }

   ~VtkhdfSeries()
   {
      if (vtk_ >= 0) { H5Gclose(vtk_); }
      if (file_ >= 0) { H5Fclose(file_); }
   }
   VtkhdfSeries(const VtkhdfSeries &) = delete;
   VtkhdfSeries &operator=(const VtkhdfSeries &) = delete;

   /// One step: `values[i]` holds ncomp_i x NumPoints() doubles of field i
   /// (point-major for vectors).
   void Append(double t, const long long meta[4], const std::vector<const double *> &values)
   {
      MFEM_VERIFY(values.size() == fields_.size(), "VtkhdfSeries::Append: "
                  << values.size() << " fields for " << fields_.size());
      for (size_t i = 0; i < fields_.size(); i++)
      {
         AppendRows("PointData/" + fields_[i].name, H5T_NATIVE_DOUBLE, fields_[i].ncomp,
                    static_cast<hsize_t>(npts_), values[i]);
         const long long off = nsteps_ * npts_;
         AppendRows("Steps/PointDataOffsets/" + fields_[i].name, H5T_NATIVE_LLONG, 1, 1, &off);
      }
      AppendRows("Steps/Values", H5T_NATIVE_DOUBLE, 1, 1, &t);
      // PartOffsets, NumberOfParts, PointOffsets, CellOffsets,
      // ConnectivityIdOffsets: the one part of the fixed grid.
      const long long so[5] = {0, 1, 0, 0, 0};
      for (int k = 0; k < 5; k++)
      {
         AppendRows(std::string("Steps/") + kStepOffsets[k], H5T_NATIVE_LLONG, 1, 1, &so[k]);
      }
      for (int k = 0; k < 4; k++)
      {
         AppendRows(std::string("/fd_cycle/") + kMeta[k], H5T_NATIVE_LLONG, 1, 1, &meta[k]);
      }
      nsteps_++;
      SetNSteps();
      Flush();
   }

   /// The last step's /fd_cycle/reason becomes `reason` (a switch or an end
   /// that falls on a step already written).
   void SetLastReason(long long reason)
   {
      MFEM_VERIFY(nsteps_ > 0, "VtkhdfSeries::SetLastReason: no step in " << path_);
      const hid_t d = Open("/fd_cycle/reason");
      const hid_t sp = H5Dget_space(d);
      const hsize_t start = static_cast<hsize_t>(nsteps_ - 1), count = 1;
      Check(H5Sselect_hyperslab(sp, H5S_SELECT_SET, &start, nullptr, &count, nullptr),
            "H5Sselect_hyperslab");
      const hid_t ms = H5Screate_simple(1, &count, nullptr);
      Check(H5Dwrite(d, H5T_NATIVE_LLONG, ms, sp, H5P_DEFAULT, &reason), "H5Dwrite (reason)");
      H5Sclose(ms);
      H5Sclose(sp);
      H5Dclose(d);
      Flush();
   }

   long long Steps() const { return nsteps_; }
   long long NumPoints() const { return npts_; }
   void Flush() { Check(H5Fflush(file_, H5F_SCOPE_GLOBAL), "H5Fflush"); }

private:
   static constexpr const char *kStepOffsets[5] = {"PartOffsets", "NumberOfParts",
                                                   "PointOffsets", "CellOffsets",
                                                   "ConnectivityIdOffsets"};
   static constexpr const char *kMeta[4] = {"regime", "phase", "reason", "step"};

   void Check(long long rc, const char *what) const
   {
      MFEM_VERIFY(rc >= 0, "snapshots: " << what << " failed for " << path_);
   }
   void Close(hid_t id, bool group = false) const
   {
      Check(group ? H5Gclose(id) : H5Dclose(id), group ? "H5Gclose" : "H5Dclose");
   }

   void AttrInt64(hid_t loc, const char *name, const long long *v, hsize_t n) const
   {
      const hid_t sp = H5Screate_simple(1, &n, nullptr);
      Check(sp, "H5Screate_simple");
      hid_t a = -1;
      if (H5Aexists(loc, name) > 0) { a = H5Aopen(loc, name, H5P_DEFAULT); }
      else { a = H5Acreate2(loc, name, H5T_STD_I64LE, sp, H5P_DEFAULT, H5P_DEFAULT); }
      Check(a, "H5Acreate2/H5Aopen");
      Check(H5Awrite(a, H5T_NATIVE_LLONG, v), "H5Awrite");
      H5Aclose(a);
      H5Sclose(sp);
   }

   /// A fixed-length string without a terminator (the VTKHDF reader in
   /// ParaView trips over a null; mesh/vtkhdf.cpp writes Type the same way).
   void AttrString(hid_t loc, const char *name, const std::string &s) const
   {
      const hid_t ty = H5Tcopy(H5T_C_S1);
      Check(ty, "H5Tcopy");
      Check(H5Tset_size(ty, std::max<size_t>(s.size(), 1)), "H5Tset_size");
      Check(H5Tset_strpad(ty, H5T_STR_NULLPAD), "H5Tset_strpad");
      const hid_t sp = H5Screate(H5S_SCALAR);
      const hid_t a = H5Acreate2(loc, name, ty, sp, H5P_DEFAULT, H5P_DEFAULT);
      Check(a, "H5Acreate2 (string)");
      const std::string v = s.empty() ? std::string(" ") : s;
      Check(H5Awrite(a, ty, v.data()), "H5Awrite (string)");
      H5Aclose(a);
      H5Sclose(sp);
      H5Tclose(ty);
   }

   void Fixed(hid_t loc, const char *name, hid_t memtype, hid_t filetype,
              std::vector<hsize_t> dims, const void *data) const
   {
      const hid_t sp = H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr);
      Check(sp, "H5Screate_simple");
      const hid_t d = H5Dcreate2(loc, name, filetype, sp, H5P_DEFAULT, H5P_DEFAULT,
                                 H5P_DEFAULT);
      Check(d, "H5Dcreate2");
      Check(H5Dwrite(d, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "H5Dwrite");
      H5Dclose(d);
      H5Sclose(sp);
   }

   void ReadFixed(hid_t loc, const char *name, hid_t memtype, void *data, size_t n) const
   {
      const hid_t d = H5Dopen2(loc, name, H5P_DEFAULT);
      Check(d, "H5Dopen2");
      const hid_t sp = H5Dget_space(d);
      MFEM_VERIFY(H5Sget_simple_extent_npoints(sp) == static_cast<hssize_t>(n),
                  "snapshots restart: " << path_ << ": " << name << " has "
                  << H5Sget_simple_extent_npoints(sp) << " values, expected " << n);
      Check(H5Dread(d, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "H5Dread");
      H5Sclose(sp);
      H5Dclose(d);
   }

   /// An empty dataset that grows along its first dimension: (0) or (0, 3).
   hid_t Growing(hid_t loc, const std::string &name, hid_t filetype, int ncomp,
                 hsize_t chunk_rows) const
   {
      const int rank = (ncomp == 1) ? 1 : 2;
      hsize_t dims[2] = {0, static_cast<hsize_t>(ncomp)};
      hsize_t maxd[2] = {H5S_UNLIMITED, static_cast<hsize_t>(ncomp)};
      hsize_t chunk[2] = {chunk_rows, static_cast<hsize_t>(ncomp)};
      const hid_t sp = H5Screate_simple(rank, dims, maxd);
      Check(sp, "H5Screate_simple");
      const hid_t pl = H5Pcreate(H5P_DATASET_CREATE);
      Check(H5Pset_chunk(pl, rank, chunk), "H5Pset_chunk");
      const hid_t d = H5Dcreate2(loc, name.c_str(), filetype, sp, H5P_DEFAULT, pl,
                                 H5P_DEFAULT);
      Check(d, "H5Dcreate2 (growing)");
      H5Pclose(pl);
      H5Sclose(sp);
      return d;
   }

   long long Rows(const std::string &name) const
   {
      const hid_t d = H5Dopen2(vtk_, name.c_str(), H5P_DEFAULT);
      Check(d, "H5Dopen2");
      const hid_t sp = H5Dget_space(d);
      hsize_t dims[2] = {0, 0};
      Check(H5Sget_simple_extent_dims(sp, dims, nullptr), "H5Sget_simple_extent_dims");
      H5Sclose(sp);
      H5Dclose(d);
      return static_cast<long long>(dims[0]);
   }

   /// Opens a dataset under /VTKHDF, or by an absolute path.
   hid_t Open(const std::string &name) const
   {
      const hid_t d = (name[0] == '/') ? H5Dopen2(file_, name.c_str(), H5P_DEFAULT)
                                       : H5Dopen2(vtk_, name.c_str(), H5P_DEFAULT);
      Check(d, "H5Dopen2");
      return d;
   }

   void Trim(const std::string &name, long long rows) const
   {
      const hid_t d = Open(name);
      const hid_t sp = H5Dget_space(d);
      hsize_t dims[2] = {0, 0};
      const int rank = H5Sget_simple_extent_dims(sp, dims, nullptr);
      Check(rank, "H5Sget_simple_extent_dims");
      H5Sclose(sp);
      MFEM_VERIFY(static_cast<long long>(dims[0]) >= rows, "snapshots restart: "
                  << path_ << ": " << name << " has " << dims[0] << " rows, fewer than "
                  << rows);
      dims[0] = static_cast<hsize_t>(rows);
      Check(H5Dset_extent(d, dims), "H5Dset_extent (trim)");
      H5Dclose(d);
   }

   void AppendRows(const std::string &name, hid_t memtype, int ncomp, hsize_t nrows,
                   const void *data) const
   {
      const hid_t d = Open(name);
      hid_t sp = H5Dget_space(d);
      hsize_t dims[2] = {0, static_cast<hsize_t>(ncomp)};
      const int rank = H5Sget_simple_extent_dims(sp, dims, nullptr);
      Check(rank, "H5Sget_simple_extent_dims");
      H5Sclose(sp);
      hsize_t start[2] = {dims[0], 0};
      hsize_t count[2] = {nrows, static_cast<hsize_t>(ncomp)};
      dims[0] += nrows;
      Check(H5Dset_extent(d, dims), "H5Dset_extent");
      sp = H5Dget_space(d);
      Check(H5Sselect_hyperslab(sp, H5S_SELECT_SET, start, nullptr, count, nullptr),
            "H5Sselect_hyperslab");
      const hid_t ms = H5Screate_simple(rank, count, nullptr);
      Check(H5Dwrite(d, memtype, ms, sp, H5P_DEFAULT, data), "H5Dwrite (append)");
      H5Sclose(ms);
      H5Sclose(sp);
      H5Dclose(d);
   }

   void SetNSteps() const
   {
      const hid_t st = H5Gopen2(vtk_, "Steps", H5P_DEFAULT);
      Check(st, "H5Gopen2 Steps");
      const long long n = nsteps_;
      AttrInt64(st, "NSteps", &n, 1);
      H5Gclose(st);
   }

   std::string path_;
   std::vector<Field> fields_;
   hid_t file_ = -1, vtk_ = -1;
   long long npts_ = 0;
   long long nsteps_ = 0;
};

/// The run's snapshot files and their schedule.  Every member function is
/// collective; only the root opens files.
class Snapshots
{
public:
   using FieldsFn = std::function<void(StationFields &)>;

   static constexpr int kFaultFields = 8;

   /// `fes`: the QD regime's displacement space (DG, vdim 3); `table` and
   /// `nbf`: the owned fault DOFs, grouped by face (nbf per face).  With
   /// `restart` (a checkpoint's GetState), the files are continued.
   Snapshots(const CycleConfig &cfg, MPI_Comm comm, ParMesh &pmesh,
             const ParFiniteElementSpace &fes, const FaultDOFTable &table, int nbf,
             const std::string &commit, const std::vector<SnapshotFileState> *restart)
      : cfg_(cfg.snapshots), comm_(comm), fes_(fes), order_(cfg.shared.mesh.order)
   {
      MPI_Comm_rank(comm_, &rank_);
      MPI_Comm_size(comm_, &size_);
      const std::string base = cfg.output.dir + "/" + cfg.output.prefix;
      if (cfg_.fault.enabled)
      {
         BuildFault_(table, nbf);
         fault_.path = base + "_fault_snapshots.vtkhdf";
         fault_.cfg = cfg_.fault;
         fault_.fields = {{"slip_strike", 1}, {"slip_dip", 1}, {"slip_rate_strike", 1},
                          {"slip_rate_dip", 1}, {"state", 1}, {"normal_stress", 1},
                          {"shear_stress_strike", 1}, {"shear_stress_dip", 1}};
         Open_(fault_, "fault", restart, commit,
               "fd_cycle fault snapshots: slip [m], slip rate [m/s] (kernel strike and dip "
               "components, as the station files), state theta [s], normal stress "
               "sigma_n_eff (compression positive) [MPa], shear stress -(tau + eta V) "
               "[MPa]; one point per fault DOF of each fault face",
               "slip m, slip_rate m/s, state s, normal_stress MPa, shear_stress MPa, time s");
      }
      if (cfg_.surface.enabled)
      {
         BuildSurface_(pmesh);
         surf_.path = base + "_surface_snapshots.vtkhdf";
         surf_.cfg = cfg_.surface;
         surf_.fields = {{"displacement", 3}, {"velocity", 3}};
         Open_(surf_, "surface", restart, commit,
               "fd_cycle free-surface snapshots: displacement [m] and velocity [m/s] in "
               "global XYZ at the points of the selected boundary faces; QD velocity: "
               "backward difference over the last QD step",
               "displacement m, velocity m/s, time s");
      }
   }

   bool Any() const { return cfg_.Any(); }

   /// After each monitored QD step at (t, step) with the QD fields, and at the
   /// start of a run that starts dynamic (reason Start).  `u`: the QD
   /// displacement of that state.  reason != Cadence writes regardless of
   /// the marks.
   void QD(real_t t, long long step, int phase, SnapReason reason, const FieldsFn &fill,
           const Vector &u)
   {
      if (!Any()) { return; }
      Vector usurf, vsurf;
      if (cfg_.surface.enabled)
      {
         SampleSurface_(u, usurf);
         vsurf.SetSize(usurf.Size());
         if (have_prev_ && t > t_prev_)
         {
            for (int i = 0; i < usurf.Size(); i++)
            {
               vsurf(i) = (usurf(i) - u_prev_(i)) / (t - t_prev_);
            }
            v_last_ = vsurf;
         }
         else if (have_prev_ && t == t_prev_ && v_last_.Size() == usurf.Size())
         {
            vsurf = v_last_;   // the same state again (the end of the run)
         }
         else { vsurf = 0.0; }
      }
      const long long mark_f = Mark_(t, cfg_.fault.dt_qd);
      const long long mark_s = Mark_(t, cfg_.surface.dt_qd);
      const bool forced = (reason != SnapReason::Cadence);
      if (cfg_.fault.enabled && t > fault_.st.last_t &&
          (forced || mark_f > fault_.st.qd_mark))
      {
         WriteFault_(t, 0, phase, reason, step, fill);
         fault_.st.qd_mark = std::max(fault_.st.qd_mark, mark_f);
      }
      else if (cfg_.fault.enabled && forced) { Relabel_(fault_, t, reason); }
      if (cfg_.surface.enabled && t > surf_.st.last_t &&
          (forced || mark_s > surf_.st.qd_mark))
      {
         WriteSurface_(t, 0, phase, reason, step, usurf, vsurf);
         surf_.st.qd_mark = std::max(surf_.st.qd_mark, mark_s);
      }
      else if (cfg_.surface.enabled && forced) { Relabel_(surf_, t, reason); }
      if (cfg_.surface.enabled && (!have_prev_ || t > t_prev_))
      {
         u_prev_ = usurf;
         t_prev_ = t;
         have_prev_ = true;
      }
   }

   /// The QD velocity history restarts from the displacement u at t: after
   /// the switch back (the QD balance) and after a restart.
   void ResetQDHistory(real_t t, const Vector &u)
   {
      if (!cfg_.surface.enabled) { return; }
      SampleSurface_(u, u_prev_);
      t_prev_ = t;
      have_prev_ = true;
      v_last_.SetSize(0);
   }

   /// A dynamic phase starts (its start time already has a snapshot).
   void EnterDynamic()
   {
      dyn_mark_f_ = 0;
      dyn_mark_s_ = 0;
   }

   /// After each dynamic step: t, the phase time tau, the dynamic fields,
   /// U and U'.  end = true: the phase's last snapshot (EndOfDynamic).
   void Dyn(real_t t, real_t tau, long long step, int phase, const FieldsFn &fill,
            const Vector &U, const Vector &Ud, bool end = false)
   {
      if (!Any()) { return; }
      const long long mf = Mark_(tau, cfg_.fault.dt_dyn);
      const long long ms = Mark_(tau, cfg_.surface.dt_dyn);
      const SnapReason r = end ? SnapReason::EndOfDynamic : SnapReason::Cadence;
      if (cfg_.fault.enabled && t > fault_.st.last_t && (end || mf > dyn_mark_f_))
      {
         WriteFault_(t, 1, phase, r, step, fill);
         dyn_mark_f_ = std::max(dyn_mark_f_, mf);
      }
      else if (cfg_.fault.enabled && end) { Relabel_(fault_, t, r); }
      if (cfg_.surface.enabled && t > surf_.st.last_t && (end || ms > dyn_mark_s_))
      {
         Vector us, vs;
         SampleSurface_(U, us);
         SampleSurface_(Ud, vs);
         WriteSurface_(t, 1, phase, r, step, us, vs);
         dyn_mark_s_ = std::max(dyn_mark_s_, ms);
      }
      else if (cfg_.surface.enabled && end) { Relabel_(surf_, t, r); }
   }

   /// For the checkpoint (the same on every rank).
   std::vector<SnapshotFileState> GetState() const
   {
      std::vector<SnapshotFileState> s;
      if (cfg_.fault.enabled) { s.push_back(fault_.st); }
      if (cfg_.surface.enabled) { s.push_back(surf_.st); }
      return s;
   }

   long long FaultSteps() const { return fault_.st.steps; }
   long long SurfaceSteps() const { return surf_.st.steps; }
   long long FaultPoints() const { return fault_.npts; }
   long long SurfacePoints() const { return surf_.npts; }
   const std::string &FaultPath() const { return fault_.path; }
   const std::string &SurfacePath() const { return surf_.path; }

private:
   struct Kind
   {
      std::string path;
      SnapshotKindConfig cfg;
      std::vector<VtkhdfSeries::Field> fields;
      SnapshotFileState st;
      std::unique_ptr<VtkhdfSeries> file;   ///< root only
      // The grid, gathered on the root.
      std::vector<double> points;
      std::vector<long long> conn, offsets;
      std::vector<unsigned char> types;
      long long npts = 0;
      // Gather layout (values per rank).
      std::vector<int> counts, displs;
      int nloc = 0;
   };

   static long long Mark_(real_t t, real_t dt)
   {
      return static_cast<long long>(std::floor(t / dt * (1.0 + 1e-12) + 1e-9));
   }

   /// The gather layout for `nloc` points per rank and the points on the root.
   void GatherPoints_(Kind &k, const std::vector<double> &xyz_loc)
   {
      k.nloc = static_cast<int>(xyz_loc.size() / 3);
      k.counts.assign(size_, 0);
      k.displs.assign(size_, 0);
      MPI_Gather(&k.nloc, 1, MPI_INT, k.counts.data(), 1, MPI_INT, 0, comm_);
      long long total = 0;
      if (rank_ == 0)
      {
         for (int r = 0; r < size_; r++) { k.displs[r] = static_cast<int>(total); total += k.counts[r]; }
      }
      MPI_Bcast(&total, 1, MPI_LONG_LONG, 0, comm_);
      k.npts = total;
      std::vector<int> c3(size_), d3(size_);
      for (int r = 0; r < size_; r++) { c3[r] = 3 * k.counts[r]; d3[r] = 3 * k.displs[r]; }
      if (rank_ == 0) { k.points.resize(static_cast<size_t>(3 * total)); }
      MPI_Gatherv(xyz_loc.data(), 3 * k.nloc, MPI_DOUBLE, k.points.data(), c3.data(),
                  d3.data(), MPI_DOUBLE, 0, comm_);
   }

   /// Gathers `nf` values per local point (point-major) into field-major
   /// arrays on the root: out[f][j] (vectors: 3 consecutive values).
   void GatherFields_(const Kind &k, const std::vector<double> &loc, const std::vector<int> &ncomp,
                      std::vector<std::vector<double>> &out)
   {
      int nv = 0;
      for (int c : ncomp) { nv += c; }
      std::vector<int> cn(size_), dn(size_);
      for (int r = 0; r < size_; r++) { cn[r] = nv * k.counts[r]; dn[r] = nv * k.displs[r]; }
      std::vector<double> all;
      if (rank_ == 0) { all.resize(static_cast<size_t>(nv * k.npts)); }
      MPI_Gatherv(loc.data(), nv * k.nloc, MPI_DOUBLE, all.data(), cn.data(), dn.data(),
                  MPI_DOUBLE, 0, comm_);
      if (rank_ != 0) { return; }
      out.assign(ncomp.size(), {});
      for (size_t f = 0; f < ncomp.size(); f++)
      {
         out[f].resize(static_cast<size_t>(ncomp[f] * k.npts));
      }
      for (long long j = 0; j < k.npts; j++)
      {
         int off = 0;
         for (size_t f = 0; f < ncomp.size(); f++)
         {
            for (int c = 0; c < ncomp[f]; c++)
            {
               out[f][static_cast<size_t>(ncomp[f] * j + c)] = all[static_cast<size_t>(nv * j + off + c)];
            }
            off += ncomp[f];
         }
      }
   }

   void Open_(Kind &k, const std::string &kind, const std::vector<SnapshotFileState> *restart,
              const std::string &commit, const std::string &description,
              const std::string &units)
   {
      k.st.kind = kind;
      if (restart)
      {
         const SnapshotFileState *found = nullptr;
         for (const SnapshotFileState &s : *restart) { if (s.kind == kind) { found = &s; } }
         MFEM_VERIFY(found, "snapshots restart: the checkpoint has no " << kind << " snapshot "
                     "record (written with [cycle.snapshots." << kind << "] off or by an "
                     "older commit); start a new run or turn it off");
         k.st = *found;
         if (rank_ == 0)
         {
            k.file = std::make_unique<VtkhdfSeries>(k.path, k.points, k.fields, k.st.steps);
         }
         return;
      }
      if (rank_ == 0)
      {
         k.file = std::make_unique<VtkhdfSeries>(
            k.path, k.points, k.conn, k.offsets, k.types, k.fields,
            std::vector<std::pair<std::string, std::string>>{
               {"description", description}, {"commit", commit}, {"units", units},
               {"regime", "0 QD, 1 dynamic"},
               {"reason", "0 cadence, 1 start, 2 switch to dynamic, 3 end of a dynamic "
                          "phase, 4 end of the run"},
               {"cadence", "dt_qd " + Num_(k.cfg.dt_qd) + " s from t = 0, dt_dyn " +
                           Num_(k.cfg.dt_dyn) + " s from each phase start"}});
      }
   }

   static std::string Num_(real_t v)
   {
      std::ostringstream o;
      o << std::setprecision(12) << v;
      return o.str();
   }

   /// p = 1: one triangle per fault face (its nbf = 3 consecutive owned
   /// DOFs, which share the bulk element); p > 1: one vertex per DOF.
   void BuildFault_(const FaultDOFTable &table, int nbf)
   {
      const int N = table.Size();
      std::vector<double> xyz(static_cast<size_t>(3 * N));
      for (int i = 0; i < 3 * N; i++) { xyz[static_cast<size_t>(i)] = table.coords(i); }
      const bool tri = (order_ == 1);
      if (tri)
      {
         MFEM_VERIFY(nbf == 3 && N % 3 == 0, "snapshots: p = 1 expects 3 DOFs per fault "
                     "face; got nbf = " << nbf << ", " << N << " owned DOFs");
         for (int f = 0; f < N / 3; f++)
         {
            const int e = table.elem[3 * f];
            MFEM_VERIFY(table.elem[3 * f + 1] == e && table.elem[3 * f + 2] == e,
                        "snapshots: owned fault DOFs " << 3 * f << ".." << 3 * f + 2
                        << " are not one face (different bulk elements)");
            const double *a = &xyz[9 * f], *b = a + 3, *c = a + 6;
            const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                 u[0] * v[1] - u[1] * v[0]};
            MFEM_VERIFY(n[0] * n[0] + n[1] * n[1] + n[2] * n[2] > 0.0,
                        "snapshots: fault face of owned DOFs " << 3 * f << ".." << 3 * f + 2
                        << " has zero area");
         }
      }
      GatherPoints_(fault_, xyz);
      if (rank_ != 0) { return; }
      const long long np = fault_.npts;
      fault_.conn.resize(static_cast<size_t>(np));
      for (long long j = 0; j < np; j++) { fault_.conn[static_cast<size_t>(j)] = j; }
      const int per = tri ? 3 : 1;
      const long long nc = np / per;
      fault_.offsets.resize(static_cast<size_t>(nc + 1));
      for (long long c = 0; c <= nc; c++) { fault_.offsets[static_cast<size_t>(c)] = per * c; }
      fault_.types.assign(static_cast<size_t>(nc), tri ? 5 : 1);   // VTK_TRIANGLE, VTK_VERTEX
   }

   /// The selected boundary faces: listed attribute, all vertices on the
   /// plane.  Each face refined `order_` times (p = 1: its vertices).
   void BuildSurface_(ParMesh &pmesh)
   {
      MFEM_VERIFY(fes_.GetVDim() == 3 && pmesh.Dimension() == 3,
                  "snapshots: the surface needs a 3-D mesh and a vector displacement space");
      const std::array<real_t, 3> &n = cfg_.surface_normal;
      std::vector<double> xyz;
      std::vector<long long> conn;
      long long nloc_cells = 0;
      Vector x(3);
      for (int be = 0; be < pmesh.GetNBE(); be++)
      {
         const int attr = pmesh.GetBdrAttribute(be);
         if (std::find(cfg_.surface_attrs.begin(), cfg_.surface_attrs.end(), attr) ==
             cfg_.surface_attrs.end()) { continue; }
         Array<int> v;
         pmesh.GetBdrElementVertices(be, v);
         bool on = true;
         for (int k = 0; k < v.Size() && on; k++)
         {
            const real_t *p = pmesh.GetVertex(v[k]);
            on = std::abs(n[0] * p[0] + n[1] * p[1] + n[2] * p[2] - cfg_.surface_offset)
                 <= cfg_.surface_tol;
         }
         if (!on) { continue; }
         FaceElementTransformations *T = pmesh.GetBdrFaceTransformations(be);
         MFEM_VERIFY(T, "snapshots: boundary element " << be << " has no face transformation");
         MFEM_VERIFY(T->GetGeometryType() == Geometry::TRIANGLE,
                     "snapshots: the surface handles triangular faces only");
         RefinedGeometry *rg = GlobGeometryRefiner.Refine(Geometry::TRIANGLE, order_, 1);
         const long long first = static_cast<long long>(xyz.size() / 3);
         for (int i = 0; i < rg->RefPts.GetNPoints(); i++)
         {
            const IntegrationPoint &fip = rg->RefPts.IntPoint(i);
            T->SetAllIntPoints(&fip);
            T->Transform(fip, x);
            for (int d = 0; d < 3; d++) { xyz.push_back(x(d)); }
            elem_.push_back(T->Elem1No);
            eip_.push_back(T->GetElement1IntPoint());
         }
         for (int i = 0; i < rg->RefGeoms.Size(); i++) { conn.push_back(first + rg->RefGeoms[i]); }
         nloc_cells += rg->RefGeoms.Size() / 3;
      }
      GatherPoints_(surf_, xyz);
      long long total_cells = 0;
      MPI_Allreduce(&nloc_cells, &total_cells, 1, MPI_LONG_LONG, MPI_SUM, comm_);
      MFEM_VERIFY(total_cells > 0, "snapshots: no boundary face with an attribute in "
                  "[cycle.snapshots.surface].attrs lies on the plane");
      // Connectivity: gather each rank's, shifted by its first point.
      int nc = static_cast<int>(conn.size());
      std::vector<int> cc(size_), cd(size_, 0);
      MPI_Gather(&nc, 1, MPI_INT, cc.data(), 1, MPI_INT, 0, comm_);
      if (rank_ == 0)
      {
         for (int r = 1; r < size_; r++) { cd[r] = cd[r - 1] + cc[r - 1]; }
         surf_.conn.resize(static_cast<size_t>(cd[size_ - 1] + cc[size_ - 1]));
      }
      MPI_Gatherv(conn.data(), nc, MPI_LONG_LONG, surf_.conn.data(), cc.data(), cd.data(),
                  MPI_LONG_LONG, 0, comm_);
      if (rank_ != 0) { return; }
      for (int r = 0; r < size_; r++)
      {
         for (int i = cd[r]; i < cd[r] + cc[r]; i++) { surf_.conn[static_cast<size_t>(i)] += surf_.displs[r]; }
      }
      const long long ncell = static_cast<long long>(surf_.conn.size() / 3);
      surf_.offsets.resize(static_cast<size_t>(ncell + 1));
      for (long long c = 0; c <= ncell; c++) { surf_.offsets[static_cast<size_t>(c)] = 3 * c; }
      surf_.types.assign(static_cast<size_t>(ncell), 5);   // VTK_TRIANGLE
   }

   /// The displacement-space vector `w` at the surface points (3 per point).
   void SampleSurface_(const Vector &w, Vector &out)
   {
      MFEM_VERIFY(w.Size() == fes_.GetVSize(), "snapshots: displacement vector size "
                  << w.Size() << " != " << fes_.GetVSize() << " (the QD regime's DG space)");
      ParGridFunction gf(const_cast<ParFiniteElementSpace *>(&fes_), const_cast<real_t *>(w.GetData()));
      out.SetSize(3 * static_cast<int>(elem_.size()));
      Vector val(3);
      for (size_t i = 0; i < elem_.size(); i++)
      {
         gf.GetVectorValue(elem_[i], eip_[i], val);
         for (int d = 0; d < 3; d++) { out(3 * static_cast<int>(i) + d) = val(d); }
      }
   }

   void WriteFault_(real_t t, int regime, int phase, SnapReason reason, long long step,
                    const FieldsFn &fill)
   {
      StationFields f;
      fill(f);
      const int N = fault_.nloc;
      MFEM_VERIFY(f.slip_strike.Size() == N && f.sigma_n.Size() == N && f.theta.Size() == N,
                  "snapshots: the regime filled " << f.slip_strike.Size() << " fault values "
                  "(sigma_n " << f.sigma_n.Size() << ") for " << N << " owned DOFs");
      std::vector<double> loc(static_cast<size_t>(kFaultFields * N));
      for (int j = 0; j < N; j++)
      {
         double *r = &loc[static_cast<size_t>(kFaultFields * j)];
         r[0] = f.slip_strike(j);
         r[1] = f.slip_dip(j);
         r[2] = f.V_strike(j);
         r[3] = f.V_dip(j);
         r[4] = f.theta(j);
         r[5] = f.sigma_n(j) / 1e6;
         r[6] = -(f.tau_a_strike(j) + f.tau_b_strike(j) + f.eta * f.V_strike(j)) / 1e6;
         r[7] = -(f.tau_a_dip(j) + f.tau_b_dip(j) + f.eta * f.V_dip(j)) / 1e6;
      }
      std::vector<std::vector<double>> all;
      GatherFields_(fault_, loc, std::vector<int>(kFaultFields, 1), all);
      Append_(fault_, t, regime, phase, reason, step, all);
   }

   void WriteSurface_(real_t t, int regime, int phase, SnapReason reason, long long step,
                      const Vector &u, const Vector &v)
   {
      const int n = surf_.nloc;
      std::vector<double> loc(static_cast<size_t>(6 * n));
      for (int j = 0; j < n; j++)
      {
         for (int d = 0; d < 3; d++)
         {
            loc[static_cast<size_t>(6 * j + d)] = u(3 * j + d);
            loc[static_cast<size_t>(6 * j + 3 + d)] = v(3 * j + d);
         }
      }
      std::vector<std::vector<double>> all;
      GatherFields_(surf_, loc, {3, 3}, all);
      Append_(surf_, t, regime, phase, reason, step, all);
   }

   /// A switch or an end at the time of the last snapshot: that snapshot
   /// takes its reason.
   void Relabel_(Kind &k, real_t t, SnapReason reason)
   {
      if (k.st.steps == 0 || t != k.st.last_t) { return; }
      if (rank_ == 0) { k.file->SetLastReason(static_cast<long long>(reason)); }
   }

   void Append_(Kind &k, real_t t, int regime, int phase, SnapReason reason, long long step,
                const std::vector<std::vector<double>> &all)
   {
      if (rank_ == 0)
      {
         std::vector<const double *> p;
         for (const auto &a : all) { p.push_back(a.data()); }
         const long long meta[4] = {regime, phase, static_cast<long long>(reason), step};
         k.file->Append(static_cast<double>(t), meta, p);
      }
      k.st.steps++;
      k.st.last_t = t;
   }

   SnapshotConfig cfg_;
   MPI_Comm comm_;
   const ParFiniteElementSpace &fes_;
   int order_ = 1;
   int rank_ = 0, size_ = 1;
   Kind fault_, surf_;
   // Surface points: local element and its reference point.
   std::vector<int> elem_;
   std::vector<IntegrationPoint> eip_;
   // QD velocity history of the surface samples.
   Vector u_prev_, v_last_;
   real_t t_prev_ = 0.0;
   bool have_prev_ = false;
   long long dyn_mark_f_ = 0, dyn_mark_s_ = 0;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_SNAPSHOTS_HPP
