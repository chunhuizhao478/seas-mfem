// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// REVIEW.md R-001 regression (PLAN_sidecar_mpi_shared_memory §7): a
// shared-mode DataField3D destroyed AFTER MPI_Finalize must not crash.
//
// This is the real lifetime of the driver's long-lived readers:
// spatial_dyn_driver calls MPI_Finalize() before main's locals unwind
// (see its MPIContext comment), so `vel_bundle`'s three windows and the
// friction readers reach ~DataField3D after finalize.  The dtor must
// detect MPI_Finalized() and skip the collective unlock/free (the OS
// reclaims the shm segment at process exit).
//
// Run single-rank:  mpirun -np 1 ./seas_test_data_field_3d_shared_mem_finalize
// PASS = prints "dtor-after-finalize OK" and exits 0.  A regression
// (dtor calling MPI_Win_unlock_all / MPI_Win_free after finalize)
// aborts inside MPI and exits nonzero.

#include "mfem.hpp"

#include "../../io/data_field_3d.hpp"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas;

namespace
{

void write_string_attr(hid_t loc, const char* name, const std::string& v)
{
   hid_t aspace = H5Screate(H5S_SCALAR);
   hid_t atype = H5Tcopy(H5T_C_S1);
   H5Tset_size(atype, v.size() + 1);
   H5Tset_strpad(atype, H5T_STR_NULLTERM);
   hid_t aid = H5Acreate2(loc, name, atype, aspace, H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(aid, atype, v.c_str());
   H5Aclose(aid);
   H5Tclose(atype);
   H5Sclose(aspace);
}

void write_double_attr(hid_t loc, const char* name, double v)
{
   hid_t aspace = H5Screate(H5S_SCALAR);
   hid_t aid = H5Acreate2(loc, name, H5T_NATIVE_DOUBLE, aspace,
                          H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(aid, H5T_NATIVE_DOUBLE, &v);
   H5Aclose(aid);
   H5Sclose(aspace);
}

void write_axis(hid_t group, const char* name, const std::vector<double>& v)
{
   hsize_t dim = v.size();
   hid_t sid = H5Screate_simple(1, &dim, nullptr);
   hid_t did = H5Dcreate2(group, name, H5T_NATIVE_DOUBLE, sid,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
   H5Dclose(did);
   H5Sclose(sid);
}

/// Minimal valid schema-v1 fixture: 2x2x2 grid, field "F" == 2.0.
void write_min_fixture(const std::string& path)
{
   hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                          H5P_DEFAULT);
   write_string_attr(file, "schema_version", "data_projection_v1");
   write_string_attr(file, "crs",            "EPSG:32611");
   write_string_attr(file, "units",          "m");
   write_string_attr(file, "z_positive",     "elevation");
   write_string_attr(file, "created_at",     "1970-01-01T00:00:00Z");
   hid_t grid = H5Gcreate2(file, "/grid", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
   write_axis(grid, "x", {0.0, 1.0});
   write_axis(grid, "y", {0.0, 1.0});
   write_axis(grid, "z", {0.0, 1.0});
   H5Gclose(grid);
   hid_t flds = H5Gcreate2(file, "/fields", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
   std::vector<double> d(8, 2.0);
   hsize_t shape[3] = {2, 2, 2};
   hid_t sid = H5Screate_simple(3, shape, nullptr);
   hid_t did = H5Dcreate2(flds, "F", H5T_NATIVE_DOUBLE, sid,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, d.data());
   write_string_attr(did, "units", "unit");
   write_double_attr(did, "min_value", 0.0);
   write_double_attr(did, "max_value", 5.0);
   H5Dclose(did);
   H5Sclose(sid);
   H5Gclose(flds);
   H5Fclose(file);
}

} // anonymous namespace


int main(int argc, char* argv[])
{
   MPI_Init(&argc, &argv);

   const char* tmpdir = std::getenv("TMPDIR");
   if (!tmpdir || !*tmpdir) { tmpdir = "/tmp"; }
   const std::string path = std::string(tmpdir)
                            + "/seas_test_df3d_sm_finalize_"
                            + std::to_string(::getpid()) + ".h5";
   ::unlink(path.c_str());
   write_min_fixture(path);

   // Heap-allocate so destruction can be deferred past MPI_Finalize —
   // exactly what happens to the driver's main-scope readers.
   auto* fld = new DataField3D(path, "F", OOBPolicy::Abort, MPI_COMM_SELF);
   const real_t v = fld->Evaluate(real_t(0.5), real_t(0.5), real_t(0.5));
   if (v != real_t(2.0))
   {
      std::fprintf(stderr,
                   "FAILED: shared-mode Evaluate returned %g, expected 2\n",
                   static_cast<double>(v));
      ::unlink(path.c_str());
      delete fld;
      MPI_Finalize();
      return 1;
   }

   ::unlink(path.c_str());
   MPI_Finalize();

   // R-001: the dtor must detect MPI_Finalized() and skip the collective
   // MPI_Win_unlock_all / MPI_Win_free.  A regression aborts right here.
   delete fld;

   std::printf("dtor-after-finalize OK\n");
   return 0;
}
