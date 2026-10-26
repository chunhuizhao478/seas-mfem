// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/fault_station_file.hpp — on-fault stations read from a file the
// deck names (benchmark plan rev 1.3, section 4.5; Chunhui Zhao, 2026-10-16).
//
// Station file: plain CSV.  Lines starting with '#' are comments; the first
// other line is the header "name,x,y,z"; every following line is one
// station: name, X, Y, Z in metres in the global frame of the mesh.
// scripts/make_fault_stations.py writes the files for TPV102, TPV205 and
// TPV104; a user can write their own.  No station position lives in C++.
//
// Snapping (SnapFaultStations): each station takes the nearest fault node by
// 3-D distance over the owned fault DOF positions of every rank, so every
// station sits on the fault.  Ties go to the lowest local index on a rank
// and then to the lowest rank (the rule of TPV102StationWriter, R-005).  The
// caller logs the requested point, the snapped node and the distance, and
// aborts when the distance exceeds the deck's snap tolerance.
//
// Writer (FaultStationWriter): one file per station,
// <prefix>_station_<name>.dat, written by the rank that owns the snapped
// DOF, with TPV102StationWriter's columns and number format
//     t slip1 slip2 V1 V2 tau1 tau2 sigma_n state
// (scientific, 10 digits, flushed every write).  The caller supplies the
// state column per law: log10 theta for the aging law, psi for the slip law
// with strong rate weakening, the slip path for slip weakening.
// SetExtraColumns (2D in-plane design rev 3, phase D: T and p of thermal
// pressurization) appends named columns after the state column, from a
// value per snapped DOF; without it every byte is as before.
//
// dynamic/tpv102_stations.hpp, which the native TPV102 driver and the SAFS
// driver use, is not edited.

#ifndef MFEM_SEAS_FAULT_STATION_FILE_HPP
#define MFEM_SEAS_FAULT_STATION_FILE_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"   // DOFData

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// One requested station: a name (used in the output file name) and a point
/// in the global frame of the mesh (m).
struct FaultStation
{
   std::string name;
   real_t x[3] = {0.0, 0.0, 0.0};
};

namespace fault_station_detail
{
inline std::string Trim(const std::string &s)
{
   size_t a = 0, b = s.size();
   while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) { a++; }
   while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) { b--; }
   return s.substr(a, b - a);
}

inline std::vector<std::string> SplitComma(const std::string &s)
{
   std::vector<std::string> f;
   std::string cur;
   std::istringstream in(s);
   while (std::getline(in, cur, ',')) { f.push_back(Trim(cur)); }
   if (!s.empty() && s.back() == ',') { f.push_back(""); }
   return f;
}
} // namespace fault_station_detail

/// Read a station file (format in the header comment).  Aborts, naming the
/// file and line, on a missing header, a line without four fields, a number
/// that does not parse or is not finite, an empty name, a name with
/// whitespace or '/', a duplicate name, or a file with no station.
inline std::vector<FaultStation> ReadFaultStationFile(const std::string &path)
{
   namespace fd = fault_station_detail;
   std::ifstream in(path);
   MFEM_VERIFY(in.good(), "station file '" << path << "' cannot be opened");
   std::vector<FaultStation> st;
   std::set<std::string> names;
   std::string line;
   int ln = 0;
   bool header = false;
   while (std::getline(in, line))
   {
      ln++;
      const std::string t = fd::Trim(line);
      if (t.empty() || t[0] == '#') { continue; }
      const std::vector<std::string> f = fd::SplitComma(t);
      MFEM_VERIFY(f.size() == 4, "station file '" << path << "', line " << ln << ": expected 4 "
                  "comma-separated fields (name,x,y,z); got " << f.size() << ": '" << t << "'");
      if (!header)
      {
         std::string h[4];
         for (int i = 0; i < 4; i++)
         {
            for (char c : f[i]) { h[i] += static_cast<char>(std::tolower(static_cast<unsigned char>(c))); }
         }
         MFEM_VERIFY(h[0] == "name" && h[1] == "x" && h[2] == "y" && h[3] == "z",
                     "station file '" << path << "', line " << ln << ": expected the header "
                     "line 'name,x,y,z'; got '" << t << "'");
         header = true;
         continue;
      }
      FaultStation s;
      s.name = f[0];
      MFEM_VERIFY(!s.name.empty(), "station file '" << path << "', line " << ln << ": empty name");
      for (char c : s.name)
      {
         MFEM_VERIFY(!std::isspace(static_cast<unsigned char>(c)) && c != '/',
                     "station file '" << path << "', line " << ln << ": station name '"
                     << s.name << "' contains whitespace or '/' (it becomes part of a file name)");
      }
      MFEM_VERIFY(names.insert(s.name).second, "station file '" << path << "', line " << ln
                  << ": duplicate station name '" << s.name << "'");
      for (int c = 0; c < 3; c++)
      {
         const char *b = f[c + 1].c_str();
         char *e = nullptr;
         const double v = std::strtod(b, &e);
         MFEM_VERIFY(!f[c + 1].empty() && e && *e == '\0' && std::isfinite(v),
                     "station file '" << path << "', line " << ln << ": station '" << s.name
                     << "': coordinate " << "xyz"[c] << " = '" << f[c + 1]
                     << "' is not a finite number");
         s.x[c] = static_cast<real_t>(v);
      }
      st.push_back(s);
   }
   MFEM_VERIFY(header, "station file '" << path << "': no header line 'name,x,y,z'");
   MFEM_VERIFY(!st.empty(), "station file '" << path << "': no stations");
   return st;
}

/// 3-D distance from a station to a point.
inline real_t FaultStationDistance(const FaultStation &s, const real_t *x)
{
   const real_t dx = x[0] - s.x[0], dy = x[1] - s.x[1], dz = x[2] - s.x[2];
   return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/// Where a station landed.  `node` and `distance` are the same on every
/// rank; `dof` is the owned DOF index on the owning rank and -1 elsewhere.
struct FaultStationSnap
{
   int dof = -1;
   bool owner = false;
   real_t node[3] = {0.0, 0.0, 0.0};
   real_t distance = 0.0;
};

/// Snap every station to its nearest fault node (3-D distance) over the
/// owned fault DOF positions `X` (3 per DOF) of all ranks.  Collective.
inline std::vector<FaultStationSnap> SnapFaultStations(const std::vector<FaultStation> &st,
                                                       const Vector &X, MPI_Comm comm)
{
   const int ns = static_cast<int>(st.size()), no = X.Size() / 3;
   int rank = 0, nranks = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nranks);
   std::vector<real_t> dloc(ns, std::numeric_limits<real_t>::max()), dmin(ns);
   std::vector<int> jloc(ns, -1);
   for (int s = 0; s < ns; s++)
   {
      for (int j = 0; j < no; j++)
      {
         const real_t d = FaultStationDistance(st[s], &X(3 * j));
         if (d < dloc[s]) { dloc[s] = d; jloc[s] = j; }
      }
   }
   MPI_Allreduce(dloc.data(), dmin.data(), ns, MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
   std::vector<FaultStationSnap> snap(ns);
   for (int s = 0; s < ns; s++)
   {
      MFEM_VERIFY(dmin[s] < std::numeric_limits<real_t>::max(),
                  "station " << st[s].name << ": the mesh has no fault DOF");
      // lowest rank among the equidistant ones writes the station (R-005)
      const bool cand = jloc[s] >= 0 && std::abs(dloc[s] - dmin[s]) < 1e-10;
      const int mine = cand ? rank : nranks;
      int win = nranks;
      MPI_Allreduce(&mine, &win, 1, MPI_INT, MPI_MIN, comm);
      real_t buf[4] = {0.0, 0.0, 0.0, 0.0};
      if (rank == win)
      {
         snap[s].owner = true;
         snap[s].dof = jloc[s];
         for (int c = 0; c < 3; c++) { buf[c] = X(3 * jloc[s] + c); }
         buf[3] = dloc[s];
      }
      MPI_Bcast(buf, 4, MPITypeMap<real_t>::mpi_type, win, comm);
      for (int c = 0; c < 3; c++) { snap[s].node[c] = buf[c]; }
      snap[s].distance = buf[3];
   }
   return snap;
}

/// Per-station time-series files (format in the header comment).
class FaultStationWriter
{
public:
   using StateFn = std::function<real_t(const DOFData &)>;
   /// Value of extra column c at owned fault DOF `dof`.
   using ExtraFn = std::function<real_t(int dof, int c)>;

   /// Extra columns after the state column (header); call before Open.
   void SetExtraColumns(const std::vector<std::string> &labels, ExtraFn value)
   {
      extra_labels_ = labels;
      extra_ = std::move(value);
   }

   /// Opens the files of the stations this rank owns.  `state_label` names
   /// the ninth column; `state` computes it from a DOF's data.
   void Open(const std::string &dir, const std::string &prefix,
             const std::vector<FaultStation> &st, const std::vector<FaultStationSnap> &snap,
             const std::string &state_label, StateFn state)
   {
      st_ = st;
      snap_ = snap;
      state_ = std::move(state);
      files_.resize(st.size());
      for (size_t s = 0; s < st.size(); s++)
      {
         if (!snap[s].owner) { continue; }
         const std::string fname = dir + "/" + prefix + "_station_" + st[s].name + ".dat";
         files_[s].open(fname);
         MFEM_VERIFY(files_[s].is_open(), "station file " << fname << " cannot be written");
         files_[s] << std::setprecision(12) << "# fault station " << st[s].name << ": requested ("
                   << st[s].x[0] << ", " << st[s].x[1] << ", " << st[s].x[2]
                   << ") m, snapped to the fault node (" << snap[s].node[0] << ", "
                   << snap[s].node[1] << ", " << snap[s].node[2] << ") m, distance "
                   << snap[s].distance << " m\n"
                   << "# Columns: time slip1 slip2 V1 V2 tau1 tau2 sigma_n " << state_label;
         for (const std::string &l : extra_labels_) { files_[s] << " " << l; }
         files_[s] << "\n"
                   << "# 1 = dip, 2 = strike; frame and signs in " << prefix
                   << "_station_conventions.txt\n";
      }
   }

   /// One row per owned station: the data of its snapped DOF.
   void WriteStep(real_t t, const std::vector<DOFData> &dd)
   {
      for (size_t s = 0; s < st_.size(); s++)
      {
         if (!snap_[s].owner || !files_[s].is_open()) { continue; }
         const DOFData &d = dd[snap_[s].dof];
         files_[s] << std::scientific << std::setprecision(10)
                   << t << " "
                   << d.slip1 << " "
                   << d.slip2 << " "
                   << d.V1 << " "
                   << d.V2 << " "
                   << d.tau1_corr << " "
                   << d.tau2_corr << " "
                   << d.sigma_n_corr << " "
                   << state_(d);
         for (size_t c = 0; c < extra_labels_.size(); c++)
         {
            files_[s] << " " << extra_(snap_[s].dof, static_cast<int>(c));
         }
         files_[s] << "\n";
         // flushed every write, so a job killed at its wall-time limit keeps its rows
         files_[s].flush();
      }
   }

   void Close()
   {
      for (auto &f : files_) { if (f.is_open()) { f.close(); } }
   }

private:
   std::vector<FaultStation> st_;
   std::vector<FaultStationSnap> snap_;
   StateFn state_;
   std::vector<std::string> extra_labels_;
   ExtraFn extra_;
   std::vector<std::ofstream> files_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_STATION_FILE_HPP
