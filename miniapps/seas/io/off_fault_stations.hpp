// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// off_fault_stations.hpp — off-fault (body) station time series for the
// Nitsche dynamic-rupture driver (rheology plan rev 1, phase R1; SCEC
// TPV26/27 description v13, Part 9).
//
// Each station is a point in the volume (on the free surface for TPV26/27).
// ParMesh::FindPoints places it in one element on one rank; that rank
// evaluates the displacement u and the velocity v there at every write and
// owns the station's file <prefix>_offfault_<name>.dat:
//
//    t  h-disp  h-vel  v-disp  v-vel  n-disp  n-vel
//
// in the station frame of the deck: h along strike shat = up x nhat_ref,
// v downward (-up), n toward the far side (-nhat_ref), nhat_ref the
// [fault_geometry] ref_normal.  For TPV26/27 (up = Z, nhat_ref = -Y) that
// is h = u_X, v = -u_Z, n = u_Y, the SCEC columns and signs of Part 9.  u
// is the displacement from the initial state (the stress-change method).
// The DG field is evaluated in the element FindPoints picked; on an
// element face that is one of the two traces.
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_OFF_FAULT_STATIONS_HPP
#define MFEM_SEAS_OFF_FAULT_STATIONS_HPP

#include "mfem.hpp"

#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

struct OffFaultStation
{
   std::string name;
   real_t x[3];
};

/// Reads "name,x,y,z" lines ('#' starts a comment; the header line
/// "name,x,y,z" is skipped).  Aborts with the line on a malformed entry.
inline std::vector<OffFaultStation> ReadOffFaultStations(const std::string &path)
{
   std::ifstream in(path);
   MFEM_VERIFY(in.is_open(), "off-fault stations: cannot open " << path);
   std::vector<OffFaultStation> out;
   std::string line;
   int ln = 0;
   while (std::getline(in, line))
   {
      ln++;
      if (line.empty() || line[0] == '#' || line.rfind("name,", 0) == 0) { continue; }
      std::stringstream ss(line);
      OffFaultStation s;
      std::string f[4];
      for (int k = 0; k < 4; k++)
      {
         MFEM_VERIFY(std::getline(ss, f[k], ','), "off-fault stations: " << path << ":" << ln
                     << ": expected name,x,y,z");
      }
      s.name = f[0];
      for (int k = 0; k < 3; k++) { s.x[k] = std::stod(f[k + 1]); }
      out.push_back(s);
   }
   MFEM_VERIFY(!out.empty(), "off-fault stations: " << path << " lists no station");
   return out;
}

class OffFaultStationWriter
{
public:
   /// strike, up, far: unit vectors of the station frame (header).
   OffFaultStationWriter(ParFiniteElementSpace &fes, const std::vector<OffFaultStation> &st,
                         const real_t strike[3], const real_t up[3], const real_t far[3])
      : fes_(fes), st_(st)
   {
      for (int i = 0; i < 3; i++) { s_[i] = strike[i]; d_[i] = -up[i]; n_[i] = far[i]; }
      ParMesh &pm = *fes.GetParMesh();
      const int np = static_cast<int>(st.size());
      DenseMatrix pts(3, np);
      for (int k = 0; k < np; k++)
      {
         for (int i = 0; i < 3; i++) { pts(i, k) = st[k].x[i]; }
      }
      pm.FindPoints(pts, elem_, ips_, false);
      for (int k = 0; k < np; k++)
      {
         int found = elem_[k] >= 0 ? 1 : (elem_[k] == -2 ? 1 : 0);
         MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, pm.GetComm());
         MFEM_VERIFY(found, "off-fault station " << st[k].name << " at (" << st[k].x[0] << ", "
                     << st[k].x[1] << ", " << st[k].x[2] << ") is not in the mesh");
      }
   }

   /// Opens the files of the stations this rank owns.
   void Open(const std::string &dir, const std::string &prefix, const std::string &problem)
   {
      files_.resize(st_.size());
      for (size_t k = 0; k < st_.size(); k++)
      {
         if (elem_[k] < 0) { continue; }
         files_[k].reset(new std::ofstream(dir + "/" + prefix + "_offfault_" + st_[k].name + ".dat"));
         std::ofstream &f = *files_[k];
         MFEM_VERIFY(f.is_open(), "off-fault station " << st_[k].name << ": cannot open its file");
         auto z = [](real_t a) { return a + 0.0; };   // prints -0 as 0
         f << "# problem=" << problem << "\n# code=seas-mfem Nitsche dynamic rupture\n"
           << "# location= (" << z(st_[k].x[0]) << ", " << z(st_[k].x[1]) << ", " << z(st_[k].x[2])
           << ") m, global frame\n"
           << "# h along strike (" << z(s_[0]) << ", " << z(s_[1]) << ", " << z(s_[2]) << "), v down ("
           << z(d_[0]) << ", " << z(d_[1]) << ", " << z(d_[2]) << "), n toward the far side ("
           << z(n_[0]) << ", " << z(n_[1]) << ", " << z(n_[2]) << ")\n"
           << "# Column #1 = Time (s)\n# Column #2 = horizontal displacement (m)\n"
           << "# Column #3 = horizontal velocity (m/s)\n# Column #4 = vertical displacement (m)\n"
           << "# Column #5 = vertical velocity (m/s)\n# Column #6 = normal displacement (m)\n"
           << "# Column #7 = normal velocity (m/s)\n"
           << "t h-disp h-vel v-disp v-vel n-disp n-vel\n";
      }
   }

   /// One line per owned station at time t; u, v are the displacement and
   /// velocity vectors of fes (size VSize).
   void WriteStep(real_t t, const Vector &u, const Vector &v)
   {
      ParGridFunction ug(&fes_, const_cast<real_t *>(u.GetData()));
      ParGridFunction vg(&fes_, const_cast<real_t *>(v.GetData()));
      Vector uu(3), vv(3);
      char buf[256];
      for (size_t k = 0; k < st_.size(); k++)
      {
         if (elem_[k] < 0) { continue; }
         ug.GetVectorValue(elem_[k], ips_[k], uu);
         vg.GetVectorValue(elem_[k], ips_[k], vv);
         auto dot = [](const real_t a[3], const Vector &b)
         { return a[0] * b(0) + a[1] * b(1) + a[2] * b(2) + 0.0; };   // + 0.0: no -0 in the files
         std::snprintf(buf, sizeof(buf), "%.12e %.9e %.9e %.9e %.9e %.9e %.9e\n", (double)t,
                       (double)dot(s_, uu), (double)dot(s_, vv), (double)dot(d_, uu),
                       (double)dot(d_, vv), (double)dot(n_, uu), (double)dot(n_, vv));
         *files_[k] << buf;
      }
   }

   void Close()
   {
      for (auto &f : files_) { if (f) { f->close(); } }
   }

   /// Number of stations this rank owns.
   int NumOwned() const
   {
      int n = 0;
      for (int e : elem_) { n += e >= 0; }
      return n;
   }

private:
   ParFiniteElementSpace &fes_;
   std::vector<OffFaultStation> st_;
   real_t s_[3], d_[3], n_[3];
   Array<int> elem_;
   Array<IntegrationPoint> ips_;
   std::vector<std::unique_ptr<std::ofstream>> files_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_OFF_FAULT_STATIONS_HPP
