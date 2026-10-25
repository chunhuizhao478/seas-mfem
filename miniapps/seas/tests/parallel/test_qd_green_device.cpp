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

// test_qd_green_device.cpp — GPU plan G3c acceptance: the device product of
// QDGreenOperator (domain/qd_green_operator.hpp, UseDevice) against its host
// product.
//
//   mpirun -np N ./seas_test_qd_green_device [--device cuda] [--nf 2000]
//        [--no-normal] [--uneven] [--gpu-mpi] [--block 64] [--repeat 5]
//        [--seed 1] [--tol 1e-13]
//
// G ((2|3)·nf rows × 2·nf columns) and T_D get random entries in [-1, 1],
// set through SetColumns and SetLoading as BuildQDGreen sets them, from
// owned tractions.  The fault DOFs are split evenly over the ranks, or with
// --uneven so that rank 0 owns half of them and the last rank none.  Three
// cases (t, slip) go through Apply() on the host, then UseDevice(--gpu-mpi)
// and Apply() again on the device:
//   zero    t = 0, slip = 0;
//   load    t = 1e9 s, slip = 0 (only t·T_D);
//   mixed   t = 3.15e9 s, slip random in [-1, 1].
// Checks, per case:
//   1. max |T_dev - T_host| / max |T_host| <= tol, shear and normal (zero:
//      T_dev = 0 exactly);
//   2. a second device Apply() gives the same bits;
//   3. the output vectors are written on the device: the host copy is
//      stale until HostRead() (checked only on the debug device, where a
//      host access without HostRead() faults).
// Printed: the differences, and the time of one Apply() (mean of --repeat)
// on the host and on the device with the bandwidth it implies (bytes of G
// per second).  Exit code 0 when every check passes.

#include "mfem.hpp"
#include "../../domain/qd_green_operator.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <iostream>
#include <sstream>
#include <random>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int myrank = 0;
int failures = 0;

void Check(bool ok, const std::string &what)
{
   if (myrank == 0)
   {
      std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << "\n";
   }
   if (!ok) { failures++; }
}

void DeviceSync()
{
   if (Device::Allows(Backend::CUDA_MASK | Backend::HIP_MASK))
   {
      MFEM_DEVICE_SYNC;
   }
}

/// Max over all ranks of max_i |a_i - b_i| and of max_i |b_i|.
void MaxDiff(const Vector &a, const Vector &b, double &diff, double &ref)
{
   const real_t *pa = a.HostRead();
   const real_t *pb = b.HostRead();
   double loc[2] = {0.0, 0.0};
   for (int i = 0; i < a.Size(); i++)
   {
      loc[0] = std::max(loc[0], double(std::abs(pa[i] - pb[i])));
      loc[1] = std::max(loc[1], double(std::abs(pb[i])));
   }
   double glob[2];
   MPI_Allreduce(loc, glob, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   diff = glob[0];
   ref = glob[1];
}

bool SameBits(const Vector &a, const Vector &b)
{
   const real_t *pa = a.HostRead();
   const real_t *pb = b.HostRead();
   int same = 1;
   for (int i = 0; i < a.Size(); i++)
   {
      if (std::memcmp(&pa[i], &pb[i], sizeof(real_t)) != 0) { same = 0; break; }
   }
   int all = 0;
   MPI_Allreduce(&same, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   return all == 1;
}

struct Case
{
   std::string name;
   real_t t = 0.0;
   Vector slip;
   Vector T_host, Tn_host;
};

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   myrank = Mpi::WorldRank();
   const int nranks = Mpi::WorldSize();
   MPI_Comm comm = MPI_COMM_WORLD;

   int nf_opt = 2000;
   bool with_normal = true;
   bool uneven = false;
   bool gpu_mpi = false;
   int block = 64;
   int repeat = 5;
   int seed = 1;
   real_t tol = 1e-13;
#ifdef MFEM_USE_CUDA
   std::string device_str = "cuda";
#else
   std::string device_str = "debug";
#endif
   OptionsParser args(argc, argv);
   args.AddOption(&device_str, "-d", "--device", "MFEM device string.");
   args.AddOption(&nf_opt, "-n", "--nf", "Fault DOFs in all.");
   args.AddOption(&with_normal, "-normal", "--normal", "-no-normal",
                  "--no-normal", "Normal rows (elastic sigma_n).");
   args.AddOption(&uneven, "-u", "--uneven", "-no-u", "--no-uneven",
                  "Rank 0 owns half the fault DOFs, the last rank none.");
   args.AddOption(&gpu_mpi, "-g", "--gpu-mpi", "-no-g", "--no-gpu-mpi",
                  "Device pointers to MPI (CUDA-aware MPI).");
   args.AddOption(&block, "-b", "--block", "Columns per SetColumns call.");
   args.AddOption(&repeat, "-r", "--repeat", "Apply() calls timed per side.");
   args.AddOption(&seed, "-s", "--seed", "Random seed.");
   args.AddOption(&tol, "-t", "--tol", "Relative tolerance.");
   args.Parse();
   if (!args.Good())
   {
      if (myrank == 0) { args.PrintUsage(std::cout); }
      return 1;
   }

   const long long nf = nf_opt;

   // One GPU per rank on a node, as seas_bp5_full chooses it.
   int device_id = 0;
#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
   if (device_str.find("cuda") != std::string::npos ||
       device_str.find("hip") != std::string::npos)
   {
      const int num_gpu = Device::GetDeviceCount();
      if (num_gpu > 1)
      {
         MPI_Comm node_comm;
         MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                             &node_comm);
         int node_rank = 0;
         MPI_Comm_rank(node_comm, &node_rank);
         MPI_Comm_free(&node_comm);
         device_id = node_rank % num_gpu;
      }
   }
#endif
   Device device(device_str, device_id);
   if (myrank == 0) { device.Print(); }
   MFEM_VERIFY(Device::Allows(Backend::DEVICE_MASK), "--device must name a "
               "device backend (cuda or debug), got '" << device_str << "'");

   // ---- Ownership ----
   std::vector<long long> nown_all(nranks, 0);
   if (uneven && nranks > 1)
   {
      nown_all[0] = nf / 2;
      const long long rest = nf - nown_all[0];
      for (int r = 1; r < nranks - 1; r++)
      {
         nown_all[r] = rest / (nranks - 2) +
                       (r - 1 < rest % (nranks - 2) ? 1 : 0);
      }
      if (nranks == 2) { nown_all[0] = nf; }
   }
   else
   {
      for (int r = 0; r < nranks; r++)
      {
         nown_all[r] = nf / nranks + (r < nf % nranks ? 1 : 0);
      }
   }
   long long first = 0;
   for (int r = 0; r < myrank; r++) { first += nown_all[r]; }
   const int nown = static_cast<int>(nown_all[myrank]);
   const int ncomp = with_normal ? 3 : 2;
   if (myrank == 0)
   {
      std::cout << "QD Green device test: N_f " << nf << ", " << ncomp
                << " components, " << nranks << " ranks (owned:";
      for (int r = 0; r < nranks; r++) { std::cout << " " << nown_all[r]; }
      std::cout << "), G " << std::fixed << std::setprecision(3)
                << double(ncomp) * nf * 2 * nf * sizeof(real_t) / 1073741824.0
                << " GiB in all" << std::defaultfloat << "\n";
   }

   // ---- G and T_D from random owned tractions ----
   QDGreenOperator green(comm, nf, first, nown, with_normal);
   std::mt19937_64 rng(static_cast<unsigned long long>(seed) * 1000003ull +
                       static_cast<unsigned long long>(myrank));
   std::uniform_real_distribution<double> uni(-1.0, 1.0);
   const auto t_fill = std::chrono::steady_clock::now();
   {
      const long long ncols = 2 * nf;
      std::vector<Vector> Ts(block), Tn(with_normal ? block : 0);
      for (long long j0 = 0; j0 < ncols; j0 += block)
      {
         const int m = static_cast<int>(std::min<long long>(block,
                                                            ncols - j0));
         for (int i = 0; i < m; i++)
         {
            Ts[i].SetSize(2 * nown);
            for (int k = 0; k < 2 * nown; k++) { Ts[i](k) = uni(rng); }
            if (with_normal)
            {
               Tn[i].SetSize(nown);
               for (int k = 0; k < nown; k++) { Tn[i](k) = uni(rng); }
            }
         }
         green.SetColumns(j0, m, Ts, Tn);
      }
      Vector TDs(2 * nown), TDn(with_normal ? nown : 0);
      for (int k = 0; k < TDs.Size(); k++) { TDs(k) = 1e-9 * uni(rng); }
      for (int k = 0; k < TDn.Size(); k++) { TDn(k) = 1e-9 * uni(rng); }
      green.SetLoading(TDs, TDn);
   }
   Check(green.Complete(), "every column and T_D set");
   const double fill_s = std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t_fill).count();

   // ---- Host product ----
   std::vector<Case> cases(3);
   cases[0].name = "zero";
   cases[0].t = 0.0;
   cases[1].name = "load";
   cases[1].t = 1e9;
   cases[2].name = "mixed";
   cases[2].t = 3.15e9;
   for (Case &c : cases)
   {
      c.slip.SetSize(2 * nown);
      c.slip = 0.0;
   }
   for (int k = 0; k < 2 * nown; k++) { cases[2].slip(k) = uni(rng); }
   for (Case &c : cases)
   {
      green.Apply(c.t, c.slip, c.T_host, with_normal ? &c.Tn_host : nullptr,
                  nullptr);
   }
   double host_s = 0.0;
   {
      Vector T, Tn;
      MPI_Barrier(comm);
      const auto t0 = std::chrono::steady_clock::now();
      for (int r = 0; r < repeat; r++)
      {
         green.Apply(cases[2].t, cases[2].slip, T,
                     with_normal ? &Tn : nullptr, nullptr);
      }
      MPI_Barrier(comm);
      host_s = std::chrono::duration<double>(
                  std::chrono::steady_clock::now() - t0).count() / repeat;
   }

   // ---- Device product ----
   const auto t_up = std::chrono::steady_clock::now();
   green.UseDevice(gpu_mpi);
   DeviceSync();
   const double upload_s = std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - t_up).count();
   if (myrank == 0)
   {
      std::cout << "Device product: " << green.ProductName()
                << (nranks > 1 ? (gpu_mpi ? ", MPI on device buffers"
                                          : ", MPI staged through the host")
                               : "")
                << "\n";
   }
   for (Case &c : cases)
   {
      for (Vector *v : {&c.slip}) { v->UseDevice(true); v->Read(); }
      Vector T, Tn, T2, Tn2;
      T.UseDevice(true);
      Tn.UseDevice(true);
      T2.UseDevice(true);
      Tn2.UseDevice(true);
      green.Apply(c.t, c.slip, T, with_normal ? &Tn : nullptr, nullptr);
      green.Apply(c.t, c.slip, T2, with_normal ? &Tn2 : nullptr, nullptr);
      DeviceSync();
      double ds = 0.0, rs = 0.0, dn = 0.0, rn = 0.0;
      MaxDiff(T, c.T_host, ds, rs);
      if (with_normal) { MaxDiff(Tn, c.Tn_host, dn, rn); }
      std::ostringstream o;
      o << std::scientific << std::setprecision(3);
      if (c.name == "zero")
      {
         o << c.name << ": device T = 0 exactly (max |T| shear " << ds;
         if (with_normal) { o << ", normal " << dn; }
         o << ")";
         Check(ds == 0.0 && dn == 0.0, o.str());
      }
      else
      {
         const double es = rs > 0.0 ? ds / rs : ds;
         const double en = rn > 0.0 ? dn / rn : dn;
         o << c.name << ": device vs host max rel. difference shear " << es;
         if (with_normal) { o << ", normal " << en; }
         o << " (tol " << tol << ")";
         Check(es <= tol && en <= tol, o.str());
      }
      const bool same = SameBits(T, T2) && (!with_normal || SameBits(Tn, Tn2));
      Check(same, c.name + ": a second device Apply() gives the same bits");
   }
   double dev_s = 0.0;
   {
      Vector T, Tn;
      T.UseDevice(true);
      Tn.UseDevice(true);
      green.Apply(cases[2].t, cases[2].slip, T, with_normal ? &Tn : nullptr,
                  nullptr);  // warm-up
      DeviceSync();
      MPI_Barrier(comm);
      const auto t0 = std::chrono::steady_clock::now();
      for (int r = 0; r < repeat; r++)
      {
         green.Apply(cases[2].t, cases[2].slip, T,
                     with_normal ? &Tn : nullptr, nullptr);
      }
      DeviceSync();
      MPI_Barrier(comm);
      dev_s = std::chrono::duration<double>(
                 std::chrono::steady_clock::now() - t0).count() / repeat;
   }

   double bytes = green.LocalBytes(), bmax = 0.0;
   MPI_Reduce(&bytes, &bmax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
   if (myrank == 0)
   {
      std::cout << std::fixed << std::setprecision(3)
                << "Times: fill " << fill_s << " s, move to device "
                << upload_s << " s\n"
                << "One Apply(): host " << host_s * 1e3 << " ms ("
                << bmax / host_s / 1e9 << " GB/s of G on the largest rank), "
                << "device " << dev_s * 1e3 << " ms ("
                << bmax / dev_s / 1e9 << " GB/s)\n";
      std::cout << (failures == 0 ? "ALL PASS" : "FAILED") << " ("
                << failures << " failed)\n";
   }
   return failures == 0 ? 0 : 1;
}
