/* gemv_bw.c - QD speed item 2: where the ranks sit, and the memory bandwidth
   they get.

     mpicc -O3 -mcpu=native gemv_bw.c -o gemv_bw
     ibrun ./gemv_bw GB_PER_NODE [fp32]

   1. Layout.  Every rank reports its node, the CPUs it may run on
      (sched_getaffinity) and the socket of the first of them
      (/sys/devices/system/cpu/cpuN/topology/physical_package_id, numbered
      0, 1, ... in the order of the node's distinct ids: Vista's Grace
      nodes report ids such as 56 and 2384); rank 0
      prints the ranks per node and socket and whether ranks are bound
      (one rank's CPU set smaller than the node).  "LAYOUT balanced" means
      every node holds the same number of ranks and, on every node, the
      sockets differ by at most one rank; otherwise "LAYOUT UNBALANCED"
      and exit code 3, before the bandwidth test.  Unbound ranks are
      counted on the socket of their first allowed CPU, so for them the
      socket counts mean nothing and only the ranks per node are checked;
      the first line says which.
   2. Bandwidth.  Each rank applies its own row block of a dense matrix
      (GB_PER_NODE x nodes in total, 20,000 columns) to a vector, best of
      20 repeats between barriers: memory traffic / time, the rate a
      bandwidth-bound kernel like the triangular solves can reach with
      this layout. */
#define _GNU_SOURCE
#include <mpi.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#define REPS 20
#define HOSTLEN 64

static int cpu_package(int cpu) {
  char path[128];
  snprintf(path, sizeof path,
           "/sys/devices/system/cpu/cpu%d/topology/physical_package_id", cpu);
  FILE *f = fopen(path, "r");
  int pkg = -1;
  if (f) { if (fscanf(f, "%d", &pkg) != 1) pkg = -1; fclose(f); }
  return pkg;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int r, p; MPI_Comm_rank(MPI_COMM_WORLD, &r); MPI_Comm_size(MPI_COMM_WORLD, &p);
  double gb_node = argc > 1 ? atof(argv[1]) : 8.0;
  int f32 = argc > 2 && !strcmp(argv[2], "fp32");

  /* ---- 1. Layout ---- */
  char host[HOSTLEN] = {0};
  gethostname(host, HOSTLEN - 1);
  cpu_set_t set; CPU_ZERO(&set);
  sched_getaffinity(0, sizeof set, &set);
  /* CPUs allowed, first allowed CPU, its socket, sockets on the node */
  int info[4] = {CPU_COUNT(&set), -1, -1, 1};
  for (int c = 0; c < CPU_SETSIZE; c++) if (CPU_ISSET(c, &set)) { info[1] = c; break; }
  long ncpu_node = sysconf(_SC_NPROCESSORS_ONLN);
  /* the node's distinct package ids, ascending (at most 8) */
  int ids[8], nid = 0;
  for (int c = 0; c < ncpu_node; c++) {
    int pk = cpu_package(c), k = 0;
    if (pk < 0) continue;
    while (k < nid && ids[k] < pk) k++;
    if (k < nid && ids[k] == pk) continue;
    if (nid == 8) continue;
    memmove(ids + k + 1, ids + k, (nid - k) * sizeof(int));
    ids[k] = pk; nid++;
  }
  if (nid > 0) info[3] = nid;
  if (info[1] >= 0) {
    int pk = cpu_package(info[1]);
    for (int k = 0; k < nid; k++) if (ids[k] == pk) { info[2] = k; break; }
  }

  char *hosts = r ? NULL : malloc((size_t)p * HOSTLEN);
  int *infos = r ? NULL : malloc((size_t)p * 4 * sizeof(int));
  MPI_Gather(host, HOSTLEN, MPI_CHAR, hosts, HOSTLEN, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Gather(info, 4, MPI_INT, infos, 4, MPI_INT, 0, MPI_COMM_WORLD);

  int bad = 0, nnodes = 0;
  if (!r) {
    /* nodes in order of first appearance */
    int *node = malloc(p * sizeof(int));
    for (int i = 0; i < p; i++) {
      node[i] = -1;
      for (int j = 0; j < i; j++)
        if (!strncmp(hosts + (size_t)i * HOSTLEN, hosts + (size_t)j * HOSTLEN, HOSTLEN)) { node[i] = node[j]; break; }
      if (node[i] < 0) node[i] = nnodes++;
    }
    int bound = 1, minw = 1 << 30, maxw = 0;
    for (int i = 0; i < p; i++) {
      int w = infos[4 * i];
      if (w < minw) minw = w;
      if (w > maxw) maxw = w;
      if (w >= ncpu_node) bound = 0;
    }
    printf("LAYOUT %d ranks on %d node(s); CPUs per rank %d..%d of %ld; ranks %s\n",
           p, nnodes, minw, maxw, ncpu_node, bound ? "bound" : "NOT bound (socket counts are first allowed CPU only)");
    int first_count = -1;
    for (int n = 0; n < nnodes; n++) {
      int cnt[8] = {0}, other = 0, tot = 0;
      for (int i = 0; i < p; i++) if (node[i] == n) {
        tot++;
        int pk = infos[4 * i + 2];
        if (pk >= 0 && pk < 8) cnt[pk]++; else other++;
      }
      int npk = 0, lo = 1 << 30, hi = 0;
      for (int i = 0; i < p; i++)
        if (node[i] == n && infos[4 * i + 3] > npk) npk = infos[4 * i + 3];
      for (int k = 0; k < npk; k++) { if (cnt[k] < lo) lo = cnt[k]; if (cnt[k] > hi) hi = cnt[k]; }
      int i0 = 0; while (node[i0] != n) i0++;
      printf("LAYOUT node %d (%s): %d ranks;", n, hosts + (size_t)i0 * HOSTLEN, tot);
      for (int k = 0; k < npk; k++) printf(" socket %d: %d", k, cnt[k]);
      if (other) printf(" unknown socket: %d", other);
      printf("\n");
      if (first_count < 0) first_count = tot;
      if (tot != first_count || (bound && (hi - lo > 1 || other))) bad = 1;
    }
    printf("LAYOUT %s\n", bad ? "UNBALANCED" : "balanced");
    fflush(stdout);
    free(node); free(hosts); free(infos);
  }
  MPI_Bcast(&bad, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nnodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (bad) { MPI_Finalize(); return 3; }

  /* ---- 2. Bandwidth ---- */
  size_t es = f32 ? 4 : 8, n = 20000;                   /* columns */
  size_t m = (size_t)(gb_node * nnodes * 1e9 / p / es / n);  /* rows per rank */
  void *A = malloc(m * n * es); double *x = malloc(n * 8), *y = malloc(m * 8);
  float *xf = malloc(n * 4);
  for (size_t j = 0; j < n; j++) { x[j] = 1.0 / (j + 1); xf[j] = (float)x[j]; }
  for (size_t i = 0; i < m * n; i++) {
    if (f32) ((float *)A)[i] = (float)(i % 97) * 1e-3f; else ((double *)A)[i] = (i % 97) * 1e-3;
  }
  double best = 1e30;
  for (int rep = 0; rep < REPS; rep++) {
    MPI_Barrier(MPI_COMM_WORLD); double t0 = MPI_Wtime();
    for (size_t i = 0; i < m; i++) {
      double s = 0.0;
      if (f32) { const float *a = (const float *)A + i * n; float sf = 0.f;
                 for (size_t j = 0; j < n; j++) sf += a[j] * xf[j];
                 s = sf; }
      else { const double *a = (const double *)A + i * n;
             for (size_t j = 0; j < n; j++) s += a[j] * x[j]; }
      y[i] = s;
    }
    MPI_Barrier(MPI_COMM_WORLD); double t = MPI_Wtime() - t0;
    if (t < best) best = t;
  }
  double tot = (double)m * n * es * p;
  if (!r) printf("GEMV %s ranks %d nodes %d matrix %.2f GB best %.4f s -> %.1f GB/s total, %.1f GB/s per node (y[0]=%g)\n",
                 f32 ? "fp32" : "fp64", p, nnodes, tot / 1e9, best, tot / best / 1e9,
                 tot / best / 1e9 / nnodes, y[0]);
  free(A); free(x); free(y); free(xf);
  MPI_Finalize(); return 0;
}
