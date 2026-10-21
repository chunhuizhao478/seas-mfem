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

// test_qd_stage_operators.cpp — GPU plan G2, acceptance tests 1-4: the
// per-stage operators of qd_stage_operators.hpp against the face loops they
// replace, on a BP5 mesh (SIPG, triangle fault faces).
//
//   mpirun -np N ./seas_test_qd_stage_operators [--mesh FILE]
//          [--mesh-scale S] [--order P] [--seed K] [--partition metis|fault]
//
//   1. 20 random physical-scale (slip, t, u): b = B·δ + D·g(t) within 1e-13
//      and T, T_n within 1e-12 of the loops (relative, max norm).
//   2. Sub-threshold slip (one face, then every face), which the loops skip:
//      ‖b_matrix − b_loops‖∞ ≤ 1e-13·‖b_loops‖∞ + ‖B‖∞·max|δ_skipped|,
//      ‖B‖∞ the exact largest absolute row sum of the assembled B.
//   3. Shared fault faces: a displacement on each single DOF of the
//      elements next to a fault face gives the same T and T_n through C_u,
//      N_u as through the loops (1e-12), including the element across a
//      partition boundary, with the face owner on either side.  A fault
//      face has one element on each side, so DOF j of every element on one
//      side is displaced at once: one displaced DOF per face per batch.
//      METIS seldom cuts along the fault on a few ranks, so --partition
//      fault puts the partition boundary on the fault plane (even N); test
//      3 then fails unless owned shared faces sit on both sides and every
//      one of them responds.
//   4. D·g(t) against AssembleDirichletLoading at t = 0, just below and
//      above the old cutoff |t·Vp| = 1e-30, at stage times across a run,
//      and for an offset loading g0 + t²h (nonzero at t = 0, not linear in
//      t): ‖D·g − b_loops‖∞ ≤ 1e-13·‖b_loops‖∞ + 1e-13·‖D‖∞·‖g‖∞, and a
//      zero loading must give exactly zero.
//
// Exit code 0 when every check passes.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../constitutive/linear_elastic.hpp"
#include "../../config/bp5_params.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

MPI_Comm comm = MPI_COMM_WORLD;
int myrank = 0;
int failures = 0;

real_t GlobalMaxAbs(const Vector &v)
{
   real_t a = 0.0;
   for (int i = 0; i < v.Size(); i++) { a = std::max(a, std::abs(v(i))); }
   MPI_Allreduce(MPI_IN_PLACE, &a, 1, MFEM_MPI_REAL_T, MPI_MAX, comm);
   return a;
}

real_t GlobalMaxAbsDiff(const Vector &a, const Vector &b)
{
   MFEM_VERIFY(a.Size() == b.Size(), "size mismatch " << a.Size() << " vs "
               << b.Size());
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++) { m = std::max(m, std::abs(a(i) - b(i))); }
   MPI_Allreduce(MPI_IN_PLACE, &m, 1, MFEM_MPI_REAL_T, MPI_MAX, comm);
   return m;
}

/// Exact largest absolute row sum over all ranks (diag + offd blocks).
real_t MaxAbsRowSum(const HypreParMatrix &A)
{
   hypre_ParCSRMatrix *h = A;
   hypre_CSRMatrix *blocks[2] = {hypre_ParCSRMatrixDiag(h),
                                 hypre_ParCSRMatrixOffd(h)};
   std::vector<real_t> rows(A.Height(), 0.0);
   for (hypre_CSRMatrix *M : blocks)
   {
      const HYPRE_Int *I = hypre_CSRMatrixI(M);
      const HYPRE_Real *V = hypre_CSRMatrixData(M);
      if (I == nullptr || hypre_CSRMatrixNumNonzeros(M) == 0) { continue; }
      for (int r = 0; r < A.Height(); r++)
      {
         for (HYPRE_Int k = I[r]; k < I[r + 1]; k++) { rows[r] += std::abs(V[k]); }
      }
   }
   real_t m = 0.0;
   for (real_t s : rows) { m = std::max(m, s); }
   MPI_Allreduce(MPI_IN_PLACE, &m, 1, MFEM_MPI_REAL_T, MPI_MAX, comm);
   return m;
}

void Check(const std::string &name, real_t err, real_t bound)
{
   const bool ok = (err <= bound);
   if (!ok) { failures++; }
   if (myrank == 0)
   {
      std::printf("  %-58s %10.3e  <= %10.3e  %s\n", name.c_str(), err, bound,
                  ok ? "PASS" : "FAIL");
   }
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   myrank = Mpi::WorldRank();
   const int nranks = Mpi::WorldSize();

   std::string mesh_file = "bp5/mesh/reference/bp5_tandem_exact.msh";
   real_t mesh_scale = 1000.0;
   int order = 1;
   int seed = 2026;
   std::string partition = "metis";
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Gmsh BP5 mesh.");
   args.AddOption(&mesh_scale, "-s", "--mesh-scale", "Coordinate scale.");
   args.AddOption(&order, "-o", "--order", "Polynomial order.");
   args.AddOption(&seed, "-r", "--seed", "Random seed.");
   args.AddOption(&partition, "-P", "--partition",
                  "metis, or fault: the partition boundary runs along the "
                  "fault (y = 0), so test 3 has shared fault faces, with "
                  "owners on both sides.");
   args.Parse();
   if (!args.Good())
   {
      if (myrank == 0) { args.PrintUsage(std::cout); }
      return 1;
   }

   // ---- Mesh and operator (as seas_bp5_full, with a switchable loading) ----
   Mesh smesh(mesh_file.c_str(), 1, 1);
   if (mesh_scale != 1.0)
   {
      for (int i = 0; i < smesh.GetNV(); i++)
      {
         real_t *v = smesh.GetVertex(i);
         for (int d = 0; d < 3; d++) { v[d] *= mesh_scale; }
      }
      smesh.SetAttributes();
   }
   MFEM_VERIFY(partition == "metis" || partition == "fault",
               "--partition must be metis or fault");
   Array<int> part;   // empty: METIS
   if (partition == "fault" && nranks > 1)
   {
      // Strips along strike (x) over the fault's extent; in each strip
      // the two sides of the fault go to two ranks.  The lower rank owns
      // a shared fault face, so the owner's side alternates strip by
      // strip (with 2 ranks: two strips, each rank on both sides).
      MFEM_VERIFY(nranks % 2 == 0, "--partition fault needs an even "
                  "number of ranks");
      real_t xmin = std::numeric_limits<real_t>::max(), xmax = -xmin;
      Array<int> bv;
      for (int be = 0; be < smesh.GetNBE(); be++)
      {
         if (smesh.GetBdrAttribute(be) != 3) { continue; }
         smesh.GetBdrElementVertices(be, bv);
         for (int v : bv)
         {
            xmin = std::min(xmin, smesh.GetVertex(v)[0]);
            xmax = std::max(xmax, smesh.GetVertex(v)[0]);
         }
      }
      const int nstrips = std::max(2, nranks / 2);
      part.SetSize(smesh.GetNE());
      Vector c(3);
      for (int e = 0; e < smesh.GetNE(); e++)
      {
         smesh.GetElementCenter(e, c);
         const int side = c(1) > 0.0 ? 1 : 0;
         int strip = static_cast<int>(nstrips * (c(0) - xmin) / (xmax - xmin));
         strip = std::min(std::max(strip, 0), nstrips - 1);
         const int flip = side ^ (strip & 1);
         part[e] = (nranks == 2) ? flip : 2 * strip + flip;
      }
   }
   ParMesh pmesh(comm, smesh, part.Size() ? part.GetData() : nullptr);
   smesh.Clear();

   BP5Params params;
   LinearElastic model(params.lambda(), params.mu());
   int loading = 0;   // 0: BP5 far-field loading, 1: offset g0 + t^2 h
   const DirichletFunc bp5_loading = MakeBP5DirichletFunc(params.Vp);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func =
      [&loading, bp5_loading](const Vector &x, real_t t, Vector &u)
   {
      if (loading == 0) { bp5_loading(x, t, u); return; }
      const real_t s = 1e-5;   // 1/m
      u.SetSize(3);
      u(0) = 0.3 + 0.1 * std::sin(s * x(0))
             + t * t * 2e-19 * (1.0 + 0.5 * std::cos(s * x(1)));
      u(1) = -0.2 + 0.05 * std::cos(s * x(2))
             + t * t * 1e-19 * std::sin(s * x(0));
      u(2) = 0.1 + t * t * 3e-19;
   };
   ElasticityDomainOperator<ParMesh> domain(
      pmesh, order, model, params.Vp, params.Wf, params.lf, bc,
      DGMethod::IP, SolverType::MUMPS);
   auto &pfes = domain.GetFESpace();
   const int nown = domain.GetNumOwnedFaultDOFs();
   const int nbf = domain.GetNbfPerFace();

   domain.BuildQDStageOperators(true);
   const QDStageOperators &ops = *domain.GetQDStageOperators();
   const long long global_ne = pmesh.GetGlobalNE();   // collective
   if (myrank == 0)
   {
      std::cout << "G2 stage-operator tests: " << mesh_file << ", p = "
                << order << ", " << nranks << " ranks, "
                << global_ne << " elements\n";
   }
   ops.PrintSummary(comm, std::cout);
   const real_t B_inf = MaxAbsRowSum(*ops.B);
   const real_t D_inf = MaxAbsRowSum(*ops.D);
   if (myrank == 0)
   {
      std::printf("  ||B||inf = %.6e   ||D||inf = %.6e\n", B_inf, D_inf);
   }

   // ---- Loops reference ----
   auto loops_b = [&](const Vector &slip_owned, real_t t, Vector &b)
   {
      Vector local, bs, bd;
      domain.ExpandOwnedToLocalFault(slip_owned, local, 2);
      domain.AssembleSlipOnlyRHS(bs, local);
      domain.AssembleDirichletOnlyRHS(bd, t);
      b = bs;
      b += bd;
   };
   auto loops_T = [&](const Vector &u, const Vector &slip_owned, Vector &T,
                      Vector &Tn)
   {
      ParGridFunction ugf(&pfes);
      ugf = u;
      Vector local, Tl, Tnl;
      domain.ExpandOwnedToLocalFault(slip_owned, local, 2);
      domain.ComputeTraction(ugf, local, Tl, &Tnl);
      domain.RestrictToOwnedFault(Tl, T, 2);
      domain.RestrictToOwnedFault(Tnl, Tn, 1);
   };

   const real_t yr = BP5Params::seconds_per_year;
   std::mt19937_64 rng_local(seed + 7919 * myrank);  // owned values
   std::mt19937_64 rng_global(seed);                  // same on every rank
   std::uniform_real_distribution<real_t> U(-1.0, 1.0);
   const int nvol = pfes.GetTrueVSize();
   Vector slip(2 * nown), u(nvol), bm, bl, Tm, Tnm, Tl, Tnl;

   // ---- 1. Random physical-scale inputs ----
   if (myrank == 0) { std::cout << "1. random (slip, t, u), 20 draws\n"; }
   {
      real_t eb = 0.0, eT = 0.0, eTn = 0.0;
      for (int k = 0; k < 20; k++)
      {
         for (int i = 0; i < slip.Size(); i++) { slip(i) = 10.0 * U(rng_local); }
         for (int i = 0; i < nvol; i++) { u(i) = U(rng_local); }
         const real_t t = 1500.0 * yr * 0.5 * (1.0 + U(rng_global));
         loops_b(slip, t, bl);
         domain.QDStageRHS(t, slip, bm);
         eb = std::max(eb, GlobalMaxAbsDiff(bm, bl) / GlobalMaxAbs(bl));
         // Traction at MPa scale (the loops print every |T| > 1 GPa)
         u *= 1e-2;
         slip *= 1e-3;
         loops_T(u, slip, Tl, Tnl);
         domain.QDStageTraction(u, slip, Tm, &Tnm);
         eT = std::max(eT, GlobalMaxAbsDiff(Tm, Tl) / GlobalMaxAbs(Tl));
         eTn = std::max(eTn, GlobalMaxAbsDiff(Tnm, Tnl) / GlobalMaxAbs(Tnl));
      }
      Check("b = B slip + D g(t), worst relative", eb, 1e-13);
      Check("T = C_u u + C_d slip, worst relative", eT, 1e-12);
      Check("T_n = N_u u + N_d slip, worst relative", eTn, 1e-12);
   }

   // ---- 2. Sub-threshold slip ----
   if (myrank == 0) { std::cout << "2. sub-threshold slip (loops skip it)\n"; }
   for (int every = 0; every <= 1; every++)
   {
      real_t skipped = 0.0;
      for (int i = 0; i < slip.Size(); i++)
      {
         slip(i) = every ? 0.9e-15 * U(rng_local) : 10.0 * U(rng_local);
         if (every) { skipped = std::max(skipped, std::abs(slip(i))); }
      }
      if (!every && myrank == 0 && nown >= nbf)
      {
         // the first owned fault face (owned DOFs 0 .. nbf-1)
         for (int i = 0; i < 2 * nbf; i++)
         {
            slip(i) = 0.9e-15 * U(rng_local);
            skipped = std::max(skipped, std::abs(slip(i)));
         }
      }
      MPI_Allreduce(MPI_IN_PLACE, &skipped, 1, MFEM_MPI_REAL_T, MPI_MAX, comm);
      const real_t t = 1500.0 * yr * 0.5 * (1.0 + U(rng_global));
      loops_b(slip, t, bl);
      domain.QDStageRHS(t, slip, bm);
      const real_t err = GlobalMaxAbsDiff(bm, bl);
      const real_t bound = 1e-13 * GlobalMaxAbs(bl) + B_inf * skipped;
      Check(every ? "every face below 1e-15: |b_m - b_l|inf"
                  : "one face below 1e-15: |b_m - b_l|inf", err, bound);
   }

   // ---- 3. Shared fault faces ----
   if (myrank == 0) { std::cout << "3. shared fault faces, displaced DOFs\n"; }
   {
      // The elements next to a fault face, by side of the fault (sign of y
      // at the centroid).  Every fault face has one element on each side,
      // so displacing local DOF j of every element on one side at once
      // puts exactly one displaced DOF into each face's traction: per face
      // this is the one-DOF column test, in 2 x (DOFs per element) batches.
      const Array<int> &shared = domain.GetFaultSharedFaces();
      const Array<int> &interior = domain.GetFaultInteriorFaces();
      const Array<int> &owned_map = domain.GetOwnedFaultFaceMap();
      const int n_int = interior.Size();
      std::set<int> side_elems[2];
      Vector c(3);
      auto side_of = [&](int e)
      {
         pmesh.GetElementCenter(e, c);
         return c(1) > 0.0 ? 1 : 0;
      };
      for (int i = 0; i < n_int; i++)
      {
         FaceElementTransformations *FTr =
            pmesh.GetInteriorFaceTransformations(interior[i]);
         if (FTr == nullptr) { continue; }
         side_elems[side_of(FTr->Elem1No)].insert(FTr->Elem1No);
         side_elems[side_of(FTr->Elem2No)].insert(FTr->Elem2No);
      }
      // owned shared faces, by the side of the owner's element
      std::vector<int> owned_shared;   // owned face indices
      int counts[3] = {0, 0, 0};       // owned shared: total, y < 0, y > 0
      for (int i = 0; i < shared.Size(); i++)
      {
         FaceElementTransformations *FTr =
            pmesh.GetSharedFaceTransformations(shared[i]);
         if (FTr == nullptr) { continue; }
         side_elems[side_of(FTr->Elem1No)].insert(FTr->Elem1No);
      }
      for (int o = 0; o < owned_map.Size(); o++)
      {
         if (owned_map[o] < n_int) { continue; }
         FaceElementTransformations *FTr =
            pmesh.GetSharedFaceTransformations(shared[owned_map[o] - n_int]);
         owned_shared.push_back(o);
         counts[0]++;
         counts[1 + side_of(FTr->Elem1No)]++;
      }
      MPI_Allreduce(MPI_IN_PLACE, counts, 3, MPI_INT, MPI_SUM, comm);
      if (myrank == 0)
      {
         std::printf("  shared fault faces: %d (owner's element on y < 0: "
                     "%d, on y > 0: %d)\n", counts[0], counts[1], counts[2]);
      }
      if (partition == "fault" && nranks > 1)
      {
         Check("fault partition: owners on both sides (0 = yes)",
               (counts[1] > 0 && counts[2] > 0) ? 0.0 : 1.0, 0.0);
      }

      const int nd = pfes.GetFE(0)->GetDof() * pfes.GetVDim();
      Vector zero_slip(2 * nown);
      zero_slip = 0.0;
      std::vector<char> hit(owned_shared.size(), 0);
      Array<int> vd;
      const char *side_name[2] = {"y < 0", "y > 0"};
      for (int side = 0; side < 2; side++)
      {
         real_t eT = 0.0, eTn = 0.0;
         for (int j = 0; j < nd; j++)
         {
            u = 0.0;
            for (int e : side_elems[side])
            {
               pfes.GetElementVDofs(e, vd);
               const int d = vd[j] >= 0 ? vd[j] : -1 - vd[j];
               u(d) = 1e-3;
            }
            loops_T(u, zero_slip, Tl, Tnl);
            domain.QDStageTraction(u, zero_slip, Tm, &Tnm);
            const real_t sT = GlobalMaxAbs(Tl), sTn = GlobalMaxAbs(Tnl);
            const real_t dT = GlobalMaxAbsDiff(Tm, Tl);
            const real_t dTn = GlobalMaxAbsDiff(Tnm, Tnl);
            // relative to the response; a zero response must stay zero
            eT = std::max(eT, sT > 0.0 ? dT / sT : (dT > 0.0 ? 1.0 : 0.0));
            eTn = std::max(eTn, sTn > 0.0 ? dTn / sTn : (dTn > 0.0 ? 1.0 : 0.0));
            for (size_t k = 0; k < owned_shared.size(); k++)
            {
               for (int q = 0; q < 2 * nbf; q++)
               {
                  if (Tl(2 * nbf * owned_shared[k] + q) != 0.0) { hit[k] = 1; }
               }
            }
         }
         Check(std::string("T,   DOFs on ") + side_name[side] +
               " displaced, worst relative", eT, 1e-12);
         Check(std::string("T_n, DOFs on ") + side_name[side] +
               " displaced, worst relative", eTn, 1e-12);
      }
      int missed = 0;
      for (char h : hit) { missed += h ? 0 : 1; }
      MPI_Allreduce(MPI_IN_PLACE, &missed, 1, MPI_INT, MPI_SUM, comm);
      if (myrank == 0)
      {
         std::printf("  %d batches (%d DOFs per element, 2 sides)\n", 2 * nd,
                     nd);
      }
      Check("shared fault faces with no traction response", missed, 0.0);
   }

   // ---- 4. D·g(t) against the Dirichlet loops ----
   if (myrank == 0) { std::cout << "4. D g(t) vs AssembleDirichletLoading\n"; }
   {
      Vector zero_slip(2 * nown);
      zero_slip = 0.0;
      const real_t tc = 1e-30 / params.Vp;   // old cutoff |t Vp| = 1e-30
      struct Case { int loading; real_t t; const char *what; };
      std::vector<Case> cases = {
         {0, 0.0, "BP5, t = 0"},
         {0, tc * (1.0 - 1e-6), "BP5, t just below the old cutoff"},
         {0, tc * (1.0 + 1e-6), "BP5, t just above the old cutoff"},
         {0, 1.0, "BP5, t = 1 s"},
         {0, yr, "BP5, t = 1 yr"},
         {0, 225.29 * yr, "BP5, t = 225.29 yr"},
         {0, 1500.0 * yr, "BP5, t = 1500 yr"},
         {1, 0.0, "offset g0 + t^2 h, t = 0"},
         {1, tc, "offset g0 + t^2 h, t at the old cutoff"},
         {1, 1e7, "offset g0 + t^2 h, t = 1e7 s"},
         {1, 100.0 * yr, "offset g0 + t^2 h, t = 100 yr"},
         {1, 1500.0 * yr, "offset g0 + t^2 h, t = 1500 yr"}};
      for (int k = 0; k < 3; k++)
      {
         cases.push_back({0, 1500.0 * yr * 0.5 * (1.0 + U(rng_global)),
                          "BP5, random stage time"});
      }
      for (const Case &c : cases)
      {
         loading = c.loading;
         domain.QDStageRHS(c.t, zero_slip, bm);
         domain.AssembleDirichletOnlyRHS(bl, c.t);
         const real_t err = GlobalMaxAbsDiff(bm, bl);
         const real_t bl_inf = GlobalMaxAbs(bl);
         const real_t g_inf = GlobalMaxAbs(ops.g);
         if (c.loading == 0 && c.t == 0.0)
         {
            // zero loading: both exactly zero
            Check(std::string(c.what) + ": |b_m|inf + |b_l|inf",
                  GlobalMaxAbs(bm) + bl_inf, 0.0);
            continue;
         }
         Check(c.what, err, 1e-13 * bl_inf + 1e-13 * D_inf * g_inf);
      }
      loading = 0;
   }

   if (myrank == 0)
   {
      std::cout << (failures == 0 ? "ALL PASSED" : "FAILURES: ")
                << (failures == 0 ? "" : std::to_string(failures)) << "\n";
   }
   return failures == 0 ? 0 : 1;
}
