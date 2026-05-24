// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// R-004 (REVIEW_speckle_seissol_drdg3d_rootcause_2026-05-24.md) — guards the
// neighbour-side read of the R-1601 batched ghost exchange used by the
// per-sub-step predictor (`WaveOperator::EvaluateBulkAtFaultQPsCanonical`).
//
// Background: the predictor exchanges the full state through a vdim=NUM_STATE,
// Ordering::byNODES ParGridFunction (`ghost_gf_full_state_`) in ONE batched
// `ExchangeFaceNbrData` (R-1601, saves 8·O collectives/macro-step vs the macro
// path's NUM_STATE per-component scalar exchanges).  The ORIGINAL unpack assumed
// the resulting `FaceNbrData()` is component-major,
//        nbr(c, nbr_elem, i)  ==  src[c * n_face_nbr_dofs + nbr_elem*ndof + i]
// (wave_operator.inl:2050-2058, 2170).  That assumption is FALSE — MFEM does NOT
// lay out a byNODES vector space's `FaceNbrData()` as contiguous component slabs,
// so the predictor read scrambled neighbour components (the SAFS shared-fault
// runaway seed).
//
// The CORRECT, layout-agnostic read is `GetFaceNbrElementVDofs(nbr_elem, vdofs)`
// → `FaceNbrData()(vdofs[k])` (the canonical MFEM pattern, cf.
// elasticity_operator_debug.inl:629-638).  For a byNODES vector space the vdof
// list is component-major, so component c / local-dof i lives at
//        vdofs[c * ndof_e + i].
//
// This test fills the local true-dofs with a per-(rank,component,dof)
// fingerprint, runs BOTH the byNODES batched exchange and the per-component
// scalar exchange (the macro path's verified-correct mechanism) on the SAME
// data, and for every shared face-neighbour element asserts:
//
//   src_full[ vdofs_full[c*ndof_e + i] ]   ==   scalar_nbr[c][ vdofs_scalar[i] ]
//   ^ FIX: GetFaceNbrElementVDofs read           ^ ground truth (scalar exchange)
//
// It also reports how many reads the ORIGINAL flat-slab formula
//   src_full[ c*n_fn + vdofs_scalar[i] ]
// gets wrong (informational — that is the bug the fix removes).
//
// PASS ⇒ the layout-agnostic GetFaceNbrElementVDofs read (now used by
// EvaluateBulkAtFaultQPsCanonical) recovers the neighbour state bit-exactly.
//
// Run: mpirun -np 2 ./seas_test_ghost_exchange_bynodes_vs_scalar

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"   // mfem::seas::NUM_STATE

#include <cmath>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const int my_rank = Mpi::WorldRank();
   const int nprocs  = Mpi::WorldSize();

   if (my_rank == 0)
   {
      std::cout << "=== test_ghost_exchange_bynodes_vs_scalar (R-004) ===\n"
                << "    NUM_STATE = " << NUM_STATE << ", nprocs = " << nprocs
                << "\n";
   }

   // Only meaningful at np >= 2 (np = 1 has no face-neighbour DOFs — the same
   // condition under which the SAFS shared-fault runaway cannot fire).
   if (nprocs < 2)
   {
      if (my_rank == 0)
      {
         std::cout << "  SKIP: needs >= 2 MPI ranks (run with mpirun -np 2).\n"
                   << "=== Summary: SKIPPED ===\n";
      }
      return 0;
   }

   // ---- Build the two spaces exactly as the WaveOperator ctor does ----------
   // (wave_operator.inl:149,160-164): scalar `pfes` (vdim=1, macro path) and the
   // vdim=NUM_STATE byNODES space (predictor path).
   Mesh serial_mesh = Mesh::MakeCartesian3D(4, 2, 2, Element::HEXAHEDRON);
   ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);
   serial_mesh.Clear();

   const int dim   = pmesh.Dimension();
   const int order = 2;
   L2_FECollection fec(order, dim);
   ParFiniteElementSpace pfes(&pmesh, &fec);                                // vdim=1
   ParFiniteElementSpace pfes_full(&pmesh, &fec, NUM_STATE, Ordering::byNODES);

   pmesh.ExchangeFaceNbrData();
   pfes.ExchangeFaceNbrData();
   pfes_full.ExchangeFaceNbrData();

   const int ndof_total = pfes.GetNDofs();

   // ---- Per-(rank,component,dof) fingerprint (exactly representable) ---------
   auto encode = [](int rank, int c, int i) -> real_t
   {
      return 1.0e6 * static_cast<real_t>(rank)
             + 1.0e3 * static_cast<real_t>(c)
             + static_cast<real_t>(i);
   };
   std::vector<real_t> Q(static_cast<size_t>(NUM_STATE) * ndof_total);
   for (int c = 0; c < NUM_STATE; c++)
      for (int i = 0; i < ndof_total; i++)
         Q[static_cast<size_t>(c) * ndof_total + i] = encode(my_rank, c, i);

   // ---- (A) byNODES batched exchange (predictor path, R-1601) ----------------
   ParGridFunction gf_full(&pfes_full);
   MFEM_VERIFY(gf_full.Size() ==
               static_cast<int>(static_cast<size_t>(NUM_STATE) * ndof_total),
               "byNODES PGF local size mismatch with NUM_STATE*ndof_total");
   std::memcpy(gf_full.GetData(), Q.data(),
               static_cast<size_t>(NUM_STATE) * ndof_total * sizeof(real_t));
   gf_full.ExchangeFaceNbrData();
   const Vector &src_full = gf_full.FaceNbrData();

   // ---- (B) per-component scalar exchange (macro path, verified-correct) -----
   ParGridFunction gf_scalar(&pfes);
   std::vector<Vector> scalar_nbr(NUM_STATE);
   for (int c = 0; c < NUM_STATE; c++)
   {
      for (int i = 0; i < ndof_total; i++)
      {
         gf_scalar(i) = Q[static_cast<size_t>(c) * ndof_total + i];
      }
      gf_scalar.ExchangeFaceNbrData();
      const Vector &src_c = gf_scalar.FaceNbrData();
      scalar_nbr[c].SetSize(src_c.Size());
      scalar_nbr[c] = src_c;   // deep copy (next iteration overwrites src_c)
   }
   const int n_fn = scalar_nbr[0].Size();   // scalar face-nbr DOF count (slab stride)

   // ---- Element-level comparison over every shared face-neighbour element ----
   const int n_fn_elems = pmesh.GetNFaceNeighborElements();
   int local_fail = 0, fix_mismatch = 0, old_mismatch = 0;
   real_t worst_fix = 0.0;
   long long compared = 0;

   for (int e = 0; e < n_fn_elems; e++)
   {
      Array<int> sv, fv;
      pfes.GetFaceNbrElementVDofs(e, sv);        // scalar: sv[i],  size ndof_e
      pfes_full.GetFaceNbrElementVDofs(e, fv);   // vector: fv[c*ndof_e + i] (byNODES)
      const int ndof_e = sv.Size();
      if (fv.Size() != NUM_STATE * ndof_e)
      {
         local_fail++;
         std::cout << "  [rank " << my_rank << "] FAIL: vector vdofs size "
                   << fv.Size() << " != NUM_STATE*ndof_e (" << NUM_STATE
                   << "*" << ndof_e << ")\n";
         continue;
      }
      for (int c = 0; c < NUM_STATE; c++)
      {
         for (int i = 0; i < ndof_e; i++)
         {
            const real_t truth = scalar_nbr[c][sv[i]];          // ground truth
            const real_t fixed = src_full[fv[c * ndof_e + i]];  // FIX (vdof read)
            const real_t old   = src_full[c * n_fn + sv[i]];    // ORIGINAL (buggy)
            compared++;
            if (fixed != truth)
            {
               fix_mismatch++;
               worst_fix = std::max(worst_fix, std::abs(fixed - truth));
               if (fix_mismatch <= 5)
               {
                  std::cout << "  [rank " << my_rank << "] FIX MISMATCH e=" << e
                            << " c=" << c << " i=" << i << " : vdof-read="
                            << fixed << " truth=" << truth << "\n";
               }
            }
            if (old != truth) { old_mismatch++; }
         }
      }
   }
   if (fix_mismatch > 0) { local_fail++; }

   // Require at least one rank to have compared something (else vacuous).
   int local_has = (compared > 0) ? 1 : 0, global_has = 0;
   MPI_Allreduce(&local_has, &global_has, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   int global_fail = 0;
   MPI_Allreduce(&local_fail, &global_fail, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   int global_old = 0;
   MPI_Allreduce(&old_mismatch, &global_old, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   std::cout << "  [rank " << my_rank << "] face_nbr_elems=" << n_fn_elems
             << " compared=" << compared << " fix_mismatch=" << fix_mismatch
             << " old_formula_mismatch=" << old_mismatch
             << " worst_fix=" << worst_fix << "\n";

   if (my_rank == 0)
   {
      if (global_has == 0)
      {
         std::cout << "  FAIL: no rank had face-neighbour elements — partition "
                      "produced no shared faces; test is vacuous.\n"
                   << "=== Summary: 0/1 passed, 1 failed ===\n";
         return 1;
      }
      std::cout << "  (informational) ORIGINAL flat-slab formula src[c*n_fn+j] "
                   "got " << global_old << " reads wrong — that was R-004.\n";
      if (global_fail == 0)
      {
         std::cout << "  PASSED: GetFaceNbrElementVDofs read of the byNODES "
                      "batched exchange is bit-identical to the per-component "
                      "scalar exchange at every (element, component, dof).\n"
                   << "  ⇒ R-004 fix verified: the predictor reads the correct "
                      "neighbour state.\n"
                   << "=== Summary: 1/1 passed, 0 failed ===\n";
      }
      else
      {
         std::cout << "  FAILED: GetFaceNbrElementVDofs read disagrees with the "
                      "scalar exchange — the byNODES vdof mapping assumed here "
                      "(vdofs[c*ndof_e+i]) is wrong.\n"
                   << "=== Summary: 0/1 passed, 1 failed ===\n";
      }
   }

   return (global_fail == 0 && global_has > 0) ? 0 : 1;
}
