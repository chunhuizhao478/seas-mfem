// Nitsche dynamic rupture, station output (speed-next O2, Chunhui Zhao
// 2026-10-20): the buffered station rows and the members-only copy average.
//
//   O2-1  StationRowFile: with flush every n the file on disk is, after each
//         row, the header and the first n*floor(rows/n) rows of the
//         flush-every-1 file; Flush writes everything; the closed files are
//         identical.
//   O2-2  FaultStationWriter: flush every 100 rows and every row give the
//         same bytes; before Close the file ends on a whole row.
//   O2-3  StationCopyAverage: the members-only version (sub, the driver's
//         default) and the world version give the same DOF data and the
//         same envelope files bit for bit, with the copies spread over all
//         ranks, kept off rank 0 (the root is not rank 0), and on the last
//         rank only (one member); the average and the envelope match a
//         serial reference; DOFs outside Dofs() are not touched.
//   O2-4  Signals (a child process, forked before MPI_Init): the first
//         SIGTERM / SIGINT sets the flag and the process keeps running; the
//         second ends it with that signal; SIGUSR1 never ends it; flush
//         every 0 rows is refused.
//
// Run on 1-4 ranks (make test-station-output, which then removes the
// files).  Output files go to test_station_output_out/ in the working
// directory.

#include "mfem.hpp"
#include "../../dynamic/fault_station_file.hpp"
#include "../../dynamic/station_copy_average.hpp"

#include <algorithm>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0, g_size = 1;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const std::string kOut = "test_station_output_out";

static int ISum(int v) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); return v; }

static std::string ReadFile(const std::string &name)
{
   std::ifstream f(name, std::ios::binary);
   std::ostringstream s;
   s << f.rdbuf();
   return s.str();
}

// ---------------------------------------------------------------------------
// O2-4: signals, in a child forked before MPI_Init (results: 1 = pass)

static int ForkCheck(void (*body)(int), int expect_signal, int expect_exit)
{
   int fd[2];
   if (pipe(fd) != 0) { return 0; }
   const pid_t pid = fork();
   if (pid < 0) { return 0; }
   if (pid == 0)
   {
      close(fd[0]);
      body(fd[1]);
      _exit(3);   // body returned: no signal ended it
   }
   close(fd[1]);
   char buf[8] = {0};
   const ssize_t n = read(fd[0], buf, sizeof(buf));
   close(fd[0]);
   int status = 0;
   waitpid(pid, &status, 0);
   const bool alive_after_first = (n > 0 && buf[0] == '1');
   if (expect_signal > 0)
   {
      return alive_after_first && WIFSIGNALED(status) && WTERMSIG(status) == expect_signal;
   }
   return alive_after_first && WIFEXITED(status) && WEXITSTATUS(status) == expect_exit;
}

static void Mark(int fd) { const char c = '1'; if (write(fd, &c, 1) != 1) { _exit(4); } }

static void TwoSignals(int fd, int sig)
{
   InstallStationSignalHandlers();
   std::raise(sig);
   if (StationSignal() != sig) { _exit(5); }
   Mark(fd);              // alive after the first signal, flag set
   std::raise(sig);       // ends the process with sig
}
static void TwoTerm(int fd) { TwoSignals(fd, SIGTERM); }
static void TwoInt(int fd) { TwoSignals(fd, SIGINT); }
static void Usr1(int fd)
{
   InstallStationSignalHandlers();
   std::raise(SIGUSR1);
   std::raise(SIGUSR1);
   std::raise(SIGUSR1);
   if (StationSignal() != SIGUSR1) { _exit(5); }
   Mark(fd);
   _exit(0);
}
static void UsrThenTerm(int fd)
{
   // the flag holds the first signal; a later SIGTERM ends the process
   InstallStationSignalHandlers();
   std::raise(SIGUSR1);
   if (StationSignal() != SIGUSR1) { _exit(5); }
   Mark(fd);
   std::raise(SIGTERM);
}
static void FlushZero(int fd)
{
   Mark(fd);
   const int null = open("/dev/null", O_WRONLY);   // the expected abort message
   if (null >= 0) { dup2(null, 2); }
   StationRowFile f;
   f.SetFlushEvery(0);    // MFEM_VERIFY: aborts
   _exit(0);
}

struct SignalResults { int term = 0, intr = 0, usr1 = 0, usr_term = 0, zero = 0; };

static SignalResults RunSignalChecks()
{
   SignalResults r;
   r.term = ForkCheck(TwoTerm, SIGTERM, 0);
   r.intr = ForkCheck(TwoInt, SIGINT, 0);
   r.usr1 = ForkCheck(Usr1, 0, 0);
   r.usr_term = ForkCheck(UsrThenTerm, SIGTERM, 0);
   // abort() when MPI is not initialized
   r.zero = ForkCheck(FlushZero, SIGABRT, 0);
   return r;
}

static void ReportSignals(const SignalResults &r)
{
   if (g_rank == 0) { std::cout << "\n--- O2-4: signals (child processes) ---\n"; }
   TEST_ASSERT(ISum(r.term) == g_size, "first SIGTERM sets the flag and the process runs on; the second ends it");
   TEST_ASSERT(ISum(r.intr) == g_size, "first SIGINT sets the flag and the process runs on; the second ends it");
   TEST_ASSERT(ISum(r.usr1) == g_size, "SIGUSR1 three times: flag set, the process exits normally");
   TEST_ASSERT(ISum(r.usr_term) == g_size, "SIGUSR1 then SIGTERM: the SIGTERM ends the process");
   TEST_ASSERT(ISum(r.zero) == g_size, "flush every 0 rows is refused (abort)");
}

// ---------------------------------------------------------------------------
// O2-1: StationRowFile

static std::string Row(int i)
{
   std::ostringstream s;
   s << std::scientific << std::setprecision(10) << 0.0125 * i << " " << std::sin(0.7 * i) << " "
     << -1.0e6 * std::cos(0.3 * i) << "\n";
   return s.str();
}

static void TestRowFile()
{
   if (g_rank == 0) { std::cout << "\n--- O2-1: StationRowFile ---\n"; }
   const std::string a = kOut + "/rowfile_every1_r" + std::to_string(g_rank) + ".dat";
   const std::string b = kOut + "/rowfile_every7_r" + std::to_string(g_rank) + ".dat";
   const std::string head = "# header line\n# columns t x y\n";
   int prefix_ok = 1, every1_ok = 1, flush_ok = 1, switch_ok = 1;
   {
      StationRowFile fa, fb;
      fa.Open(a, 1);
      fb.Open(b, 7);
      fa.Text() << head;
      fb.Text() << head;
      std::string all = head;   // what the files hold in the end
      const int n1 = 23, n2 = 40;
      for (int i = 1; i <= n2; i++)
      {
         const std::string r = Row(i);
         fa.Text() << r; fa.EndRow();
         fb.Text() << r; fb.EndRow();
         all += r;
         // every 1: everything so far is on disk
         if (ReadFile(a) != all) { every1_ok = 0; }
         if (i <= n1)
         {
            // every 7: the header and the first 7*floor(i/7) rows
            std::string want;
            const int k = 7 * (i / 7);
            if (k > 0) { want = head; for (int j = 1; j <= k; j++) { want += Row(j); } }
            if (ReadFile(b) != want) { prefix_ok = 0; }
         }
         if (i == n1)
         {
            fb.Flush();
            if (ReadFile(b) != all) { flush_ok = 0; }
            fb.SetFlushEvery(1);
         }
         if (i > n1 && ReadFile(b) != all) { switch_ok = 0; }
      }
   }
   TEST_ASSERT(ISum(every1_ok) == g_size, "flush every row: the file holds every row after each row");
   TEST_ASSERT(ISum(prefix_ok) == g_size, "flush every 7 rows: the file holds the header and 7*floor(rows/7) whole rows");
   TEST_ASSERT(ISum(flush_ok) == g_size, "Flush writes the rows held in memory");
   TEST_ASSERT(ISum(switch_ok) == g_size, "after SetFlushEvery(1) every row reaches the file");
   const int same = ReadFile(a) == ReadFile(b) && !ReadFile(a).empty();
   TEST_ASSERT(ISum(same) == g_size, "closed files: flush every 1 and every 7 are byte-identical");
}

// ---------------------------------------------------------------------------
// O2-2: FaultStationWriter

static DOFData Data(int g, int step)
{
   DOFData d;
   const real_t x = 0.37 * g + 0.011 * step;
   d.slip1 = 1.0e-3 * std::sin(x);
   d.slip2 = 2.5 * std::cos(1.3 * x) + 0.01 * step;
   d.V1 = 1.0e-4 * std::sin(2.1 * x + 0.4);
   d.V2 = 1.7 * std::cos(0.9 * x) * std::exp(-0.001 * step);
   d.tau1_corr = 1.0e5 * std::sin(0.5 * x + 1.0);
   d.tau2_corr = 7.0e7 + 3.0e6 * std::cos(x);
   d.sigma_n_corr = 1.2e8 + 1.0e6 * std::sin(0.2 * x);
   d.psi = 0.6 + 0.05 * std::cos(3.0 * x);
   d.slip_rate = std::hypot(d.V1, d.V2);
   return d;
}

static void TestWriter()
{
   if (g_rank == 0) { std::cout << "\n--- O2-2: FaultStationWriter ---\n"; }
   const int ns = 5, no = 8, nsteps = 250;
   std::vector<FaultStation> st(ns);
   std::vector<FaultStationSnap> snap(ns);
   for (int s = 0; s < ns; s++)
   {
      st[s].name = "st" + std::to_string(s);
      st[s].x[0] = 1000.0 * s; st[s].x[1] = 0.0; st[s].x[2] = -500.0 * s;
      snap[s].owner = (s % g_size == g_rank);
      snap[s].dof = snap[s].owner ? (3 * s + 1) % no : -1;
      for (int c = 0; c < 3; c++) { snap[s].node[c] = st[s].x[c] + 0.5; }
      snap[s].distance = std::sqrt(0.75);
   }
   auto state = [](const DOFData &d) { return std::log10(d.psi); };
   // same prefix (the header names it), two folders
   const std::string d1 = kOut + "/writer_every1", d100 = kOut + "/writer_every100";
   if (g_rank == 0) { mkdir(d1.c_str(), 0755); mkdir(d100.c_str(), 0755); }
   MPI_Barrier(MPI_COMM_WORLD);
   FaultStationWriter w1, w100;
   w1.Open(d1, "wr", st, snap, "log10_theta", state, 1);
   w100.Open(d100, "wr", st, snap, "log10_theta", state, 100);
   std::vector<DOFData> dd(no);
   for (int k = 0; k < nsteps; k++)
   {
      for (int j = 0; j < no; j++) { dd[j] = Data(100 * g_rank + j, k); }
      w1.WriteStep(1.0e-3 * k, dd);
      w100.WriteStep(1.0e-3 * k, dd);
   }
   int whole = 1, rows200 = 1;
   for (int s = 0; s < ns; s++)
   {
      if (!snap[s].owner) { continue; }
      const std::string f1 = d1 + "/wr_station_" + st[s].name + ".dat";
      const std::string f100 = d100 + "/wr_station_" + st[s].name + ".dat";
      const std::string t1 = ReadFile(f1), t100 = ReadFile(f100);
      if (t100.empty() || t100.back() != '\n' || t1.compare(0, t100.size(), t100) != 0) { whole = 0; }
      int rows = 0;
      std::istringstream in(t100);
      for (std::string l; std::getline(in, l);) { if (!l.empty() && l[0] != '#') { rows++; } }
      if (rows != 200) { rows200 = 0; }
   }
   TEST_ASSERT(ISum(whole) == g_size, "before Close: the every-100 file is a prefix of the every-row file ending on a whole row");
   TEST_ASSERT(ISum(rows200) == g_size, "before Close: 250 rows written, 200 on disk with flush every 100");
   w1.Close();
   w100.Close();
   int same = 1, full = 1;
   for (int s = 0; s < ns; s++)
   {
      if (!snap[s].owner) { continue; }
      const std::string t1 = ReadFile(d1 + "/wr_station_" + st[s].name + ".dat");
      const std::string t100 = ReadFile(d100 + "/wr_station_" + st[s].name + ".dat");
      if (t1 != t100) { same = 0; }
      int rows = 0;
      std::istringstream in(t1);
      for (std::string l; std::getline(in, l);) { if (!l.empty() && l[0] != '#') { rows++; } }
      if (rows != nsteps) { full = 0; }
   }
   TEST_ASSERT(ISum(same) == g_size, "after Close: flush every 100 and every row are byte-identical");
   TEST_ASSERT(ISum(full) == g_size, "after Close: every file holds all 250 rows");
   std::vector<int> want;
   for (int s = 0; s < ns; s++) { if (snap[s].owner) { want.push_back(snap[s].dof); } }
   TEST_ASSERT(ISum(w1.Dofs() == want) == g_size, "Dofs() lists the snapped DOFs of the owned stations");
}

// ---------------------------------------------------------------------------
// O2-3: StationCopyAverage, members only vs world

// A synthetic fault: per station a node with k copies (DG DOFs at the same
// point up to 1e-9 m, each with its own face centroid and area; two copies
// share a centroid, as the two sides of a face do), and filler DOFs far from
// every station.  Every rank builds the whole list and keeps its own part.
struct SynthDof { int g; int station; real_t x[3]; real_t c[3]; real_t area; };

static std::vector<SynthDof> Synth(std::vector<FaultStation> &st)
{
   const real_t node[3][3] = {{1000.0, 0.0, -2000.0}, {-4500.0, 0.0, -7500.0}, {12000.0, 0.0, -500.0}};
   const int kc[3] = {6, 4, 9};
   st.assign(3, {});
   std::vector<SynthDof> all;
   int g = 0;
   for (int s = 0; s < 3; s++)
   {
      st[s].name = "cs" + std::to_string(s);
      st[s].x[0] = node[s][0] + 0.25; st[s].x[1] = node[s][1] + 0.1; st[s].x[2] = node[s][2] - 0.2;
      for (int k = 0; k < kc[s]; k++)
      {
         SynthDof d;
         d.g = g++;
         d.station = s;
         for (int c = 0; c < 3; c++) { d.x[c] = node[s][c]; }
         d.x[0] += 1.0e-9 * (k % 3);   // round-off scatter, inside kTie
         d.x[2] -= 1.0e-10 * k;
         const int kk = (k == 1) ? 0 : k;   // copies 0 and 1 share a centroid
         const real_t a = 0.9 * kk + 0.3 * s;
         d.c[0] = node[s][0] + 60.0 * std::cos(a);
         d.c[1] = 0.0;
         d.c[2] = node[s][2] + 60.0 * std::sin(a);
         d.area = 4000.0 + 731.0 * k + 97.0 * s;
         all.push_back(d);
      }
   }
   for (int f = 0; f < 40; f++)
   {
      SynthDof d;
      d.g = g++;
      d.station = -1;
      d.x[0] = -20000.0 + 1000.0 * f + 333.0; d.x[1] = 0.0; d.x[2] = -150.0 - 400.0 * (f % 7);
      for (int c = 0; c < 3; c++) { d.c[c] = d.x[c] + 20.0; }
      d.area = 5000.0 + 10.0 * f;
      all.push_back(d);
   }
   return all;
}

// The rank that holds DOF d: "spread" over all ranks, "skip0" with the
// copies off rank 0, "last" with every copy on the last rank.
static int Owner(const SynthDof &d, const std::string &layout)
{
   if (d.station < 0 || g_size == 1) { return (7 * d.g + 3) % g_size; }
   if (layout == "skip0") { return 1 + (5 * d.g + 1) % (g_size - 1); }
   if (layout == "last") { return g_size - 1; }
   return (7 * d.g + 3) % g_size;
}

static std::string EnvelopeRow(real_t t, int nc, const real_t *lo, const real_t *hi)
{
   // StationCopyAverage::WriteEnvelope_'s format: sigma_n, tau1, tau2, V1, V2
   std::ostringstream s;
   s << std::scientific << std::setprecision(10) << t << " " << nc << " "
     << lo[6] << " " << hi[6] << " " << lo[4] << " " << hi[4] << " " << lo[5] << " "
     << hi[5] << " " << lo[2] << " " << hi[2] << " " << lo[3] << " " << hi[3] << "\n";
   return s.str();
}

static void Fields(const DOFData &q, real_t *v)
{
   const real_t f[8] = {q.slip1, q.slip2, q.V1, q.V2, q.tau1_corr, q.tau2_corr, q.sigma_n_corr, q.psi};
   for (int i = 0; i < 8; i++) { v[i] = f[i]; }
}

static bool SameBits(const DOFData &a, const DOFData &b)
{
   real_t va[8], vb[8];
   Fields(a, va);
   Fields(b, vb);
   return std::memcmp(va, vb, sizeof(va)) == 0
          && std::memcmp(&a.slip_rate, &b.slip_rate, sizeof(real_t)) == 0;
}

static void TestCopyAverage(const std::string &layout)
{
   if (g_rank == 0) { std::cout << "\n--- O2-3: StationCopyAverage, layout " << layout << " ---\n"; }
   std::vector<FaultStation> st;
   const std::vector<SynthDof> all = Synth(st);
   const int ns = static_cast<int>(st.size());
   // this rank's DOFs, in a scrambled local order
   std::vector<SynthDof> mine;
   for (const SynthDof &d : all) { if (Owner(d, layout) == g_rank) { mine.push_back(d); } }
   std::sort(mine.begin(), mine.end(), [](const SynthDof &a, const SynthDof &b) {
      return (2654435761u * static_cast<unsigned>(a.g)) % 1000003u
             < (2654435761u * static_cast<unsigned>(b.g)) % 1000003u;
   });
   const int no = static_cast<int>(mine.size());
   Vector X(3 * no), ca(4 * no);
   for (int j = 0; j < no; j++)
   {
      for (int c = 0; c < 3; c++) { X(3 * j + c) = mine[j].x[c]; ca(4 * j + c) = mine[j].c[c]; }
      ca(4 * j + 3) = mine[j].area;
   }
   const std::string dw = kOut + "/copies_world_" + layout, ds = kOut + "/copies_sub_" + layout;
   if (g_rank == 0) { mkdir(dw.c_str(), 0755); mkdir(ds.c_str(), 0755); }
   MPI_Barrier(MPI_COMM_WORLD);
   StationCopyAverage world(ca, X, st, dw, "syn", false, 1);
   StationCopyAverage sub(ca, X, st, ds, "syn", true, 3);

   // members: the ranks that hold a station copy; the root is the lowest
   int member = 0;
   for (const SynthDof &d : mine) { member = member || d.station >= 0; }
   std::vector<int> mem(g_size);
   MPI_Allgather(&member, 1, MPI_INT, mem.data(), 1, MPI_INT, MPI_COMM_WORLD);
   int nmem = 0, lowest = -1;
   for (int r = 0; r < g_size; r++) { if (mem[r]) { nmem++; if (lowest < 0) { lowest = r; } } }
   const int members_ok = sub.Members() == (member ? nmem : 0) && world.Members() == g_size;
   TEST_ASSERT(ISum(members_ok) == g_size, "Members(): " << nmem << " of " << g_size
               << " rank(s) hold copies (world: all ranks)");
   const int root_ok = sub.Root() == (g_rank == lowest) && world.Root() == (g_rank == 0);
   TEST_ASSERT(ISum(root_ok) == g_size, "Root(): the envelope is written on rank " << lowest
               << " (world: rank 0)");
   std::vector<int> want;
   for (int j = 0; j < no; j++) { if (mine[j].station >= 0) { want.push_back(j); } }
   TEST_ASSERT(ISum(sub.Dofs() == want && world.Dofs() == want) == g_size,
               "Dofs(): exactly the local station copies");

   // serial reference: per station the copies in centroid order
   std::vector<std::vector<int>> ref(ns);
   for (size_t i = 0; i < all.size(); i++) { if (all[i].station >= 0) { ref[all[i].station].push_back(i); } }
   for (int s = 0; s < ns; s++)
   {
      std::stable_sort(ref[s].begin(), ref[s].end(), [&](int a, int b) {
         return std::lexicographical_compare(all[a].c, all[a].c + 3, all[b].c, all[b].c + 3);
      });
   }
   const int nsteps = 25;
   int bits_ok = 1, untouched_ok = 1, avg_ok = 1, rate_ok = 1;
   real_t max_rel = 0.0;
   std::vector<std::string> env_ref(ns);
   for (int k = 0; k < nsteps; k++)
   {
      const real_t t = 0.004 * k;
      std::vector<DOFData> d0(no), dw_(no), ds_(no);
      for (int j = 0; j < no; j++) { d0[j] = Data(mine[j].g, k); }
      dw_ = d0;
      ds_ = d0;
      world.Apply(t, dw_);
      sub.Apply(t, ds_);
      for (int j = 0; j < no; j++)
      {
         if (!SameBits(dw_[j], ds_[j])) { bits_ok = 0; }
         if (mine[j].station < 0 && !SameBits(ds_[j], d0[j])) { untouched_ok = 0; }
      }
      for (int s = 0; s < ns; s++)
      {
         real_t wsum = 0.0, avg[8] = {0}, mag[8] = {0}, lo[8], hi[8], v[8];
         for (int f = 0; f < 8; f++) { lo[f] = std::numeric_limits<real_t>::max(); hi[f] = -lo[f]; }
         for (int i : ref[s])
         {
            Fields(Data(all[i].g, k), v);
            wsum += all[i].area;
            for (int f = 0; f < 8; f++)
            {
               avg[f] += all[i].area * v[f];
               mag[f] += all[i].area * std::abs(v[f]);
               lo[f] = std::min(lo[f], v[f]);
               hi[f] = std::max(hi[f], v[f]);
            }
         }
         for (int f = 0; f < 8; f++) { avg[f] /= wsum; mag[f] /= wsum; }
         env_ref[s] += EnvelopeRow(t, static_cast<int>(ref[s].size()), lo, hi);
         for (int j = 0; j < no; j++)
         {
            if (mine[j].station != s) { continue; }
            Fields(ds_[j], v);
            for (int f = 0; f < 8; f++)
            {
               // relative to the weighted mean magnitude (a mean near 0 cancels)
               const real_t rel = std::abs(v[f] - avg[f]) / std::max(mag[f], real_t(1e-300));
               max_rel = std::max(max_rel, rel);
               if (rel > 1e-14) { avg_ok = 0; }
            }
            if (ds_[j].slip_rate != std::hypot(ds_[j].V1, ds_[j].V2)) { rate_ok = 0; }
         }
      }
   }
   MPI_Allreduce(MPI_IN_PLACE, &max_rel, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, MPI_COMM_WORLD);
   TEST_ASSERT(ISum(bits_ok) == g_size, "members only == world, bit for bit, on every DOF over "
               << nsteps << " steps");
   TEST_ASSERT(ISum(untouched_ok) == g_size, "DOFs away from the stations are not changed");
   TEST_ASSERT(ISum(avg_ok) == g_size, "station DOFs hold the area-weighted average (max rel. diff. vs serial "
               << std::scientific << std::setprecision(2) << max_rel << std::defaultfloat << ")");
   TEST_ASSERT(ISum(rate_ok) == g_size, "slip_rate = |(V1, V2)| of the average");

   // envelope files: flush both, compare on rank 0
   world.Flush();
   sub.Flush();
   MPI_Barrier(MPI_COMM_WORLD);
   int env_same = 1, env_ref_ok = 1;
   if (g_rank == 0)
   {
      for (int s = 0; s < ns; s++)
      {
         const std::string n = "/syn_station_" + st[s].name + "_copies.dat";
         const std::string tw = ReadFile(dw + n), ts = ReadFile(ds + n);
         if (tw.empty() || tw != ts) { env_same = 0; }
         // the rows after the two header lines == the serial reference
         size_t p = 0;
         for (int h = 0; h < 2 && p != std::string::npos; h++) { p = tw.find('\n', p); if (p != std::string::npos) { p++; } }
         if (p == std::string::npos || tw.substr(p) != env_ref[s]) { env_ref_ok = 0; }
      }
   }
   MPI_Bcast(&env_same, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&env_ref_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
   TEST_ASSERT(env_same, "envelope files: members only (flush every 3) == world (every row), byte for byte");
   TEST_ASSERT(env_ref_ok, "envelope rows (copies, min, max) == the serial reference");
}

int main(int argc, char *argv[])
{
   // O2-4 forks children; fork before MPI_Init
   const SignalResults sig = RunSignalChecks();

   Mpi::Init(argc, argv);
   g_rank = Mpi::WorldRank();
   g_size = Mpi::WorldSize();
   if (g_rank == 0)
   {
      std::cout << "=== Station output (speed-next O2) on " << g_size << " rank(s) ===\n";
      mkdir(kOut.c_str(), 0755);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   TestRowFile();
   TestWriter();
   for (const char *layout : {"spread", "skip0", "last"}) { TestCopyAverage(layout); }
   ReportSignals(sig);

   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                << " passed, " << num_failed << " failed ===\n";
   }
   return (num_failed == 0) ? 0 : 1;
}
