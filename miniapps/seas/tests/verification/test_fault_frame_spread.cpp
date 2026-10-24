// Nitsche dynamic rupture, fault frames compressed on the other ranks
// (fault-writer note rev1, W0/W1, Chunhui Zhao 2026-10-24).
//
//   FS-1  VTKHDF::FirstAppendChunk: a first append up to 0.5 MB keeps its
//         length (at least 1024 rows); a longer one gets its largest divisor
//         up to 0.5 MB (TPV205 200 m: 314736 values, chunks of 52456), or
//         0.5 MB when no divisor reaches 1/16 of it.
//   FS-2  FaultHDFStream with SetSpread(true) and SetSpread(false) write the
//         same fault.vtkhdf: every dataset with the same type, shape,
//         chunk length, filters and values; every chunk of every dataset
//         with the same stored bytes and filter mask; every attribute with
//         the same bytes.  Cases: frames of several chunks per field
//         (deflate 1 and 6), one chunk per field with fewer chunks than
//         helpers, ranks holding uneven parts (one holding none), random
//         bits that deflate makes larger than they are, and a frame with no
//         chunk length that divides it (the spread falls back to rank 0,
//         SpreadActive false).  On one rank the spread is off.
//   FS-3  The spread is on exactly when it can be: more than one rank,
//         deflate, a dividing chunk length; with no filter it is off.
//   FS-4  The comparison sees one value changed by one part in 10^12 in the
//         last frame, in the values and, with the values left out, in the
//         stored chunks.
//
// Run on 1-4 ranks (make test-fault-frame-spread, which then removes the
// files).  Output files go to test_fault_frame_spread_out/ in the working
// directory.
//
// --compare a.vtkhdf b.vtkhdf [--values-only] (one rank) compares two fault
// files the same way and exits 0 when they agree: attributes, dataset types
// and shapes, values (blocks of 64 MB), and the stored bytes of every chunk
// with the chunk lengths and filters; --values-only leaves out the chunks,
// the chunk lengths and the filters (files written before W0 have other
// chunk lengths).

#include "mfem.hpp"
#include "../../io/fault_vtkhdf_writer.hpp"

#include <hdf5.h>
#if !H5_VERSION_GE(1, 12, 0)
#error "test_fault_frame_spread reads the files with the HDF5 1.12 object API (H5Ovisit3)"
#endif

#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <sys/stat.h>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0, g_size = 1;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const std::string kOut = "test_fault_frame_spread_out";

// ---------------------------------------------------------------------------
// Two files compared: attributes and dataset shapes from a listing, then each
// dataset's values and stored chunks block by block (a Vista fault file holds
// 15 GB of values)

struct Listing
{
   std::map<std::string, std::string> small;   ///< attributes and shapes -> bytes
   std::vector<std::string> datasets;
   std::string error;
   bool values_only = false;
};

static std::string Size(const std::string &s) { return std::to_string(s.size()) + " bytes"; }

static void Attributes(hid_t obj, const std::string &path, Listing &c)
{
   H5O_info2_t info;
   H5Oget_info3(obj, &info, H5O_INFO_NUM_ATTRS);
   for (hsize_t i = 0; i < info.num_attrs; i++)
   {
      const hid_t a = H5Aopen_by_idx(obj, ".", H5_INDEX_NAME, H5_ITER_INC, i, H5P_DEFAULT,
                                     H5P_DEFAULT);
      char name[256] = {0};
      H5Aget_name(a, sizeof(name), name);
      const hid_t t = H5Aget_type(a), s = H5Aget_space(a);
      const hssize_t n = H5Sget_simple_extent_npoints(s);
      std::string b(static_cast<std::size_t>(n) * H5Tget_size(t), '\0');
      if (H5Tis_variable_str(t) > 0) { c.error = path + "@" + name + " is a variable-length string"; }
      else if (!b.empty()) { H5Aread(a, t, &b[0]); }
      c.small["attr " + path + "@" + name] = b;
      H5Sclose(s);
      H5Tclose(t);
      H5Aclose(a);
   }
}

/// type, dims, and (unless values_only) chunk length and filters
static std::string Shape(hid_t d, bool values_only)
{
   const hid_t t = H5Dget_type(d), s = H5Dget_space(d), dcpl = H5Dget_create_plist(d);
   const int nd = H5Sget_simple_extent_ndims(s);
   std::vector<hsize_t> dims(nd);
   H5Sget_simple_extent_dims(s, dims.data(), nullptr);
   std::string shape = "type " + std::to_string(H5Tget_class(t)) + " size " +
                       std::to_string(H5Tget_size(t)) + " dims";
   for (hsize_t x : dims) { shape += " " + std::to_string(x); }
   if (!values_only && H5Pget_layout(dcpl) == H5D_CHUNKED)
   {
      std::vector<hsize_t> ch(nd);
      H5Pget_chunk(dcpl, nd, ch.data());
      shape += " chunk";
      for (hsize_t x : ch) { shape += " " + std::to_string(x); }
   }
   const int nf = values_only ? 0 : H5Pget_nfilters(dcpl);
   for (int i = 0; i < nf; i++)
   {
      unsigned flags = 0, cd[8] = {0};
      size_t ncd = 8;
      char fname[64];
      const H5Z_filter_t id = H5Pget_filter2(dcpl, i, &flags, &ncd, cd, sizeof(fname), fname,
                                             nullptr);
      shape += " filter " + std::to_string(id) + " flags " + std::to_string(flags);
      for (size_t j = 0; j < ncd && j < 8; j++) { shape += " " + std::to_string(cd[j]); }
   }
   H5Pclose(dcpl);
   H5Sclose(s);
   H5Tclose(t);
   return shape;
}

static herr_t Visit(hid_t root, const char *name, const H5O_info2_t *info, void *data)
{
   Listing &c = *static_cast<Listing *>(data);
   const std::string path = name;
   const hid_t o = H5Oopen(root, name, H5P_DEFAULT);
   Attributes(o, path, c);
   if (info->type == H5O_TYPE_DATASET)
   {
      c.small["shape " + path] = Shape(o, c.values_only);
      c.datasets.push_back(path);
   }
   H5Oclose(o);
   return 0;
}

static Listing List(hid_t f, bool values_only)
{
   Listing c;
   c.values_only = values_only;
   H5Ovisit3(f, H5_INDEX_NAME, H5_ITER_INC, Visit, &c, H5O_INFO_BASIC);
   return c;
}

/// "" when the values (blocks of rows) and, unless values_only, the stored
/// chunks of dataset path agree; counts the chunks and value bytes compared
static std::string CompareDataset(hid_t fa, hid_t fb, const std::string &path, bool values_only,
                                  bool chunks_only, long long &chunks, long long &bytes)
{
   const hid_t da = H5Dopen2(fa, path.c_str(), H5P_DEFAULT), db = H5Dopen2(fb, path.c_str(), H5P_DEFAULT);
   const hid_t t = H5Dget_type(da), s = H5Dget_space(da);
   const int nd = H5Sget_simple_extent_ndims(s);
   std::vector<hsize_t> dims(nd);
   H5Sget_simple_extent_dims(s, dims.data(), nullptr);
   std::size_t row = H5Tget_size(t);
   for (int i = 1; i < nd; i++) { row *= dims[i]; }
   const hsize_t rows = nd > 0 ? dims[0] : 1;
   const hsize_t block = std::max<hsize_t>(1, (64u << 20) / std::max<std::size_t>(row, 1));
   std::string why;
   std::vector<char> va, vb;
   for (hsize_t r0 = 0; r0 < rows && why.empty() && nd > 0 && !chunks_only; r0 += block)
   {
      std::vector<hsize_t> start(nd, 0), count = dims;
      start[0] = r0;
      count[0] = std::min(block, rows - r0);
      va.resize(count[0] * row);
      vb.resize(count[0] * row);
      const hid_t ms = H5Screate_simple(nd, count.data(), nullptr);
      const hid_t sa = H5Dget_space(da), sb = H5Dget_space(db);
      H5Sselect_hyperslab(sa, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr);
      H5Sselect_hyperslab(sb, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr);
      const herr_t ea = H5Dread(da, t, ms, sa, H5P_DEFAULT, va.data());
      const herr_t eb = H5Dread(db, t, ms, sb, H5P_DEFAULT, vb.data());
      H5Sclose(sb);
      H5Sclose(sa);
      H5Sclose(ms);
      if (ea < 0 || eb < 0) { why = "values " + path + ": read failed"; }
      else if (std::memcmp(va.data(), vb.data(), va.size()) != 0)
      {
         why = "values " + path + " rows " + std::to_string(r0) + ".." +
               std::to_string(r0 + count[0] - 1) + " differ";
      }
      bytes += static_cast<long long>(va.size());
   }
   const hid_t dcpl = H5Dget_create_plist(da);
   if (why.empty() && !values_only && H5Pget_layout(dcpl) == H5D_CHUNKED)
   {
      hsize_t na = 0, nb = 0;
      H5Dget_num_chunks(da, H5S_ALL, &na);
      H5Dget_num_chunks(db, H5S_ALL, &nb);
      if (na != nb) { why = "chunks " + path + ": " + std::to_string(na) + " vs " + std::to_string(nb); }
      std::string ra, rb;
      for (hsize_t k = 0; k < na && why.empty(); k++)
      {
         std::vector<hsize_t> off(nd);
         unsigned ma = 0, mb = 0;
         haddr_t aa = 0, ab = 0;
         hsize_t za = 0, zb = 0;
         H5Dget_chunk_info(da, H5S_ALL, k, off.data(), &ma, &aa, &za);
         std::string key = "chunk " + path;
         for (hsize_t x : off) { key += " " + std::to_string(x); }
         if (H5Dget_chunk_info_by_coord(db, off.data(), &mb, &ab, &zb) < 0 || ab == HADDR_UNDEF)
         {
            why = key + " only in the first file";
            break;
         }
         ra.assign(static_cast<std::size_t>(za), '\0');
         rb.assign(static_cast<std::size_t>(zb), '\0');
         uint32_t fa_ = 0, fb_ = 0;
         if (za > 0) { H5Dread_chunk(da, H5P_DEFAULT, off.data(), &fa_, &ra[0]); }
         if (zb > 0) { H5Dread_chunk(db, H5P_DEFAULT, off.data(), &fb_, &rb[0]); }
         if (ma != mb || fa_ != fb_ || ra != rb)
         {
            why = key + " differs (" + Size(ra) + " vs " + Size(rb) + ", masks " +
                  std::to_string(ma) + " vs " + std::to_string(mb) + ")";
         }
         chunks++;
      }
   }
   H5Pclose(dcpl);
   H5Sclose(s);
   H5Tclose(t);
   H5Dclose(db);
   H5Dclose(da);
   return why;
}

struct Comparison
{
   std::string why;   ///< "" = the same
   long long chunks = 0, bytes = 0;
   std::size_t datasets = 0;
};

/// file a against file b: attributes, shapes, then values and chunks
/// (chunks_only: the test's check that the chunks alone see a change)
static Comparison CompareFiles(const std::string &a, const std::string &b, bool values_only,
                               bool chunks_only = false)
{
   Comparison r;
   const hid_t fa = H5Fopen(a.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   const hid_t fb = fa < 0 ? -1 : H5Fopen(b.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (fa < 0 || fb < 0)
   {
      r.why = "cannot open " + (fa < 0 ? a : b);
      if (fa >= 0) { H5Fclose(fa); }
      return r;
   }
   const Listing la = List(fa, values_only), lb = List(fb, values_only);
   if (!la.error.empty() || !lb.error.empty()) { r.why = la.error.empty() ? lb.error : la.error; }
   for (const auto &kv : la.small)
   {
      if (!r.why.empty()) { break; }
      const auto it = lb.small.find(kv.first);
      if (it == lb.small.end()) { r.why = kv.first + " only in the first file"; }
      else if (it->second != kv.second)
      {
         r.why = kv.first + " differs (" + Size(kv.second) + " vs " + Size(it->second) + ")";
      }
   }
   for (const auto &kv : lb.small)
   {
      if (!r.why.empty()) { break; }
      if (la.small.find(kv.first) == la.small.end()) { r.why = kv.first + " only in the second file"; }
   }
   for (const std::string &d : la.datasets)
   {
      if (!r.why.empty()) { break; }
      r.why = CompareDataset(fa, fb, d, values_only, chunks_only, r.chunks, r.bytes);
   }
   r.datasets = la.datasets.size();
   H5Fclose(fb);
   H5Fclose(fa);
   return r;
}

// ---------------------------------------------------------------------------
// A fault of ncell triangles split over the ranks; fields change per frame

struct Case
{
   std::string name;
   int ncell = 0;
   int nfields = 0, nconst = 0, frames = 0;
   double level = 1.0;   ///< deflate level; < 0 = no filter
   bool empty_rank = false;   ///< the last rank holds no cells
   bool noise = false;   ///< changing fields of random bits: deflate output larger than its input
};

static double Noise(unsigned long long x)
{
   // splitmix64
   x += 0x9e3779b97f4a7c15ULL;
   x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
   x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
   x ^= x >> 31;
   double d;
   std::memcpy(&d, &x, sizeof(d));
   return d;
}

/// cells [c0, c1) of this rank
static void MyCells(const Case &k, int &c0, int &c1)
{
   // parts proportional to rank + 1 (uneven), or the last rank none
   const int nparts = k.empty_rank && g_size > 1 ? g_size - 1 : g_size;
   long long wsum = 0, w0 = 0;
   for (int r = 0; r < nparts; r++) { wsum += r + 1; if (r < g_rank) { w0 += r + 1; } }
   if (g_rank >= nparts) { c0 = c1 = k.ncell; return; }
   c0 = static_cast<int>(k.ncell * w0 / wsum);
   c1 = static_cast<int>(k.ncell * (w0 + g_rank + 1) / wsum);
}

struct Written
{
   bool active = false;
   int helpers = 0;
   long long chunk = 0;
};

static Written WriteCase(const Case &k, bool spread, const std::string &dir,
                         bool perturb = false)
{
   int c0 = 0, c1 = 0;
   MyCells(k, c0, c1);
   vtu::LocalFaultPack pack;
   for (int g = c0; g < c1; g++)
   {
      const int b = static_cast<int>(pack.vertices.size());
      const double x = 10.0 * (g % 1000), z = -10.0 * (g / 1000);
      pack.vertices.push_back({x, 0.0, z});
      pack.vertices.push_back({x + 10.0, 0.0, z});
      pack.vertices.push_back({x, 0.0, z - 10.0});
      pack.triangles.push_back({b, b + 1, b + 2});
   }
   std::vector<char> constant(k.nfields, 0);
   for (int j = 0; j < k.nfields; j++)
   {
      // names out of alphabetical order, as the driver's
      pack.field_names.push_back(std::string(1, static_cast<char>('z' - j)) + "_field");
      pack.field_arrays.emplace_back(c1 - c0, 0.0);
      if (j >= k.nfields - k.nconst) { constant[j] = 1; }
   }
   vtkhdf::FaultHDFStream stream;
   stream.SetCompression(ParaViewHDFDataCollection::HDFCompression::Deflate, k.level);
   stream.SetSpread(spread);
   for (int f = 0; f < k.frames; f++)
   {
      for (int j = 0; j < k.nfields; j++)
      {
         for (int g = c0; g < c1; g++)
         {
            // smooth in g with a rupture front (zeros ahead of it), as slip rates
            const double front = 0.25 * k.ncell * (f + 1) / k.frames;
            double v = constant[j] ? 0.008 + 1e-6 * std::cos(0.01 * g + j)
                       : (g < front ? 1e6 * std::sin(1e-4 * g + 0.3 * f + j) + j : 0.0);
            if (j == 0 && g % 7 == 0) { v = 1.0 / (1.0 + g); }
            if (k.noise && !constant[j])
            {
               v = Noise(static_cast<unsigned long long>(g) + 1000003ULL * f + 1000000007ULL * j);
            }
            if (perturb && j == 0 && f == k.frames - 1 && g == k.ncell / 8 + 1) { v *= 1.0 + 1e-12; }
            pack.field_arrays[j][g - c0] = v;
         }
      }
      stream.Write(dir, 0.25 * f, pack, constant, g_rank, g_size, MPI_COMM_WORLD);
   }
   Written w;
   w.active = stream.SpreadActive();
   w.helpers = stream.SpreadHelpers();
   w.chunk = stream.ChunkValues();
   return w;
}

/// FS-2 for one case: the spread and the rank-0 files hold the same content
static void RunCase(const Case &k, bool expect_active, long long expect_chunk)
{
   const std::string a = kOut + "/" + k.name + "_spread", b = kOut + "/" + k.name + "_rank0";
   const Written ws = WriteCase(k, true, a);
   const Written wr = WriteCase(k, false, b);
   MPI_Barrier(MPI_COMM_WORLD);
   const bool on = expect_active && g_size > 1;
   TEST_ASSERT(ws.active == on && !wr.active, k.name << ": spread " << (ws.active ? "on" : "off")
               << " (expected " << (on ? "on" : "off") << "), rank0 " << (wr.active ? "on" : "off"));
   if (on)
   {
      const int ntask = (k.nfields - k.nconst) * static_cast<int>(3LL * k.ncell / expect_chunk);
      TEST_ASSERT(ws.chunk == expect_chunk && ws.helpers == std::min(ntask, g_size - 1),
                  k.name << ": chunks of " << ws.chunk << " values (expected " << expect_chunk
                  << "), " << ws.helpers << " helpers (expected " << std::min(ntask, g_size - 1)
                  << ")");
   }
   int ok = 1;
   std::string why;
   long long nch = 0;
   if (g_rank == 0)
   {
      const Comparison c = CompareFiles(a + "/fault.vtkhdf", b + "/fault.vtkhdf", false);
      why = c.why;
      ok = why.empty() ? 1 : 0;
      nch = c.chunks;
   }
   MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
   TEST_ASSERT(ok == 1, k.name << ": the two files hold the same datasets, chunks ("
               << nch << ") and attributes" << (ok ? "" : ": " + why));
}

/// FS-4: a file with one value changed differs from the spread file in a chunk
static void RunPerturbed(const Case &k)
{
   const std::string a = kOut + "/" + k.name + "_spread", b = kOut + "/" + k.name + "_changed";
   WriteCase(k, true, a);
   WriteCase(k, false, b, true);
   MPI_Barrier(MPI_COMM_WORLD);
   int ok = 1;
   std::string why;
   if (g_rank == 0)
   {
      why = CompareFiles(a + "/fault.vtkhdf", b + "/fault.vtkhdf", false).why;
      // the stored chunks see it too, with the values left out
      const std::string wc = CompareFiles(a + "/fault.vtkhdf", b + "/fault.vtkhdf", false, true).why;
      ok = (why.rfind("values ", 0) == 0 && wc.rfind("chunk ", 0) == 0) ? 1 : 0;
      why += "; chunks alone: " + wc;
   }
   MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
   TEST_ASSERT(ok == 1, "FS-4 " << k.name << ": one value changed by 1e-12 is found: "
               << (why.empty() ? std::string("no difference") : why));
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   g_rank = Mpi::WorldRank();
   g_size = Mpi::WorldSize();
   for (int i = 1; i < argc; i++)
   {
      if (std::string(argv[i]) != "--compare") { continue; }
      MFEM_VERIFY(i + 2 < argc && g_size == 1, "usage: seas_test_fault_frame_spread --compare "
                  "a.vtkhdf b.vtkhdf [--values-only], on one rank");
      bool values_only = false;
      for (int j = 1; j < argc; j++) { values_only = values_only || std::string(argv[j]) == "--values-only"; }
      const auto t0 = std::chrono::steady_clock::now();
      const Comparison c = CompareFiles(argv[i + 1], argv[i + 2], values_only);
      const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      std::cout << (c.why.empty() ? "SAME" : "DIFFERENT") << ": " << argv[i + 1] << " vs "
                << argv[i + 2] << (values_only ? " (values only)" : "") << ": " << c.datasets
                << " datasets, " << c.bytes << " value bytes, " << c.chunks << " chunks compared in "
                << sec << " s" << (c.why.empty() ? "" : "; first difference: " + c.why) << "\n";
      return c.why.empty() ? 0 : 1;
   }
   if (g_rank == 0)
   {
      std::cout << "test_fault_frame_spread on " << g_size << " rank(s)\n";
      ::mkdir(kOut.c_str(), 0775);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   // FS-1
   {
      const hsize_t cmax = 65536;
      TEST_ASSERT(VTKHDF::FirstAppendChunk(314736, cmax) == 52456,
                  "FS-1 314736 values (TPV205 200 m): chunks of 52456");
      TEST_ASSERT(VTKHDF::FirstAppendChunk(6000, cmax) == 6000
                  && VTKHDF::FirstAppendChunk(100, cmax) == 1024
                  && VTKHDF::FirstAppendChunk(65536, cmax) == 65536,
                  "FS-1 up to 0.5 MB: the append's length, at least 1024");
      TEST_ASSERT(VTKHDF::FirstAppendChunk(3 * 65537, cmax) == cmax
                  && VTKHDF::FirstAppendChunk(131072, cmax) == 65536
                  && VTKHDF::FirstAppendChunk(78684, cmax) == 39342,
                  "FS-1 longer: the largest divisor, or 65536 when none reaches 4096");
   }

   // FS-2 / FS-3
   // TPV205 200 m frame length, 6 chunks per field, deflate 1 (the decks)
   RunCase({"tpv205_200m_d1", 104912, 6, 2, 4, 1.0, false}, true, 52456);
   // deflate 6 (the writer's default), uneven parts with one rank empty
   RunCase({"d6_empty_rank", 104912, 4, 1, 3, 6.0, true}, true, 52456);
   // one chunk per field, one changing field: fewer chunks than helpers
   RunCase({"one_chunk", 2000, 3, 2, 3, 1.0, false}, true, 6000);
   // random bits: chunks that deflate makes larger, stored as HDF5 stores them
   RunCase({"noise", 104912, 3, 1, 2, 1.0, false, true}, true, 52456);
   // no dividing chunk length (3 x 65537 values): rank 0 alone
   RunCase({"no_divisor", 65537, 3, 1, 2, 1.0, false}, false, 0);
   // no filter: rank 0 alone
   RunCase({"no_filter", 2000, 3, 1, 2, -1.0, false}, false, 0);
   // FS-4
   RunPerturbed({"perturbed", 104912, 3, 1, 2, 1.0, false});

   if (g_rank == 0)
   {
      std::cout << "\n" << num_passed << "/" << num_tests << " passed, " << num_failed
                << " failed\n";
   }
   return (num_failed == 0) ? 0 : 1;
}
