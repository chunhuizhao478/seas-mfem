// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_fault_frame_bit_census.cpp — Phase 0 of
// document/PLAN_unify_interior_shared_fault_substep_2026-07-09.md
//
// READ-ONLY DIAGNOSTIC.  This is the plan's STOP GATE: it decides whether the
// plan's stated root cause is real before any behavior changes.
//
// Background.  The canonical fault frame `can_n` is built from `qpd.normal`
// with a sign factor, and the codebase uses TWO DIFFERENT BITS for that factor:
//
//   interior fault sites  (wave_operator.inl:2221, 2972, 4217):
//       can_n = (!elem1_on_plus) ? -qpd.normal : qpd.normal
//   shared   fault sites  (wave_operator.inl:2513, 3646, 5252):
//       can_n = ( qpd.sign_flipped) ? -qpd.normal : qpd.normal
//
// The plan hypothesised that these two rules disagree (the "frame bit" defect),
// and that R-1601 disabled sub-stepping on shared faces because of it.
//
// TWO MEASUREMENTS, because the hypothesis has two separable parts:
//
//   C1 (per rank): count QPs where `sign_flipped != !elem1_on_plus`.
//       Analytic prediction: ZERO.  MFEM's face normal points away from Elem1,
//       and `elem1_on_plus` is derived by projecting (elem1_c - face_c) onto the
//       canonicalised normal, so the two bits are the SAME bit computed two ways
//       (see the plan's Phase-1 requirement 0).  If C1 == 0 everywhere, the
//       "two different rules" framing is FALSIFIED and the plan must be revised.
//
//   C2 (cross rank, np>1): for the SAME physical shared fault face, do the two
//       owning ranks reconstruct the SAME `can_n`, or OPPOSITE ones?
//       This is what R-1601's comment actually claims
//       ("rank A's Elem1 = rank B's Elem2 -> opposite canonical frames").
//       Measured for BOTH bit choices, so we learn whether swapping the bit
//       (the plan's Phase 1) would change anything at all.
//
// Robust to per-rank QP ordering: C2 compares the FACE-AVERAGED can_n direction
// (mean over the face's QPs, then normalised), so it does not depend on how the
// two ranks order their quadrature points.
//
// Fixture: the 2-tet shared-fault mesh + y-sign partition lifted verbatim from
// tests/parallel/test_bimaterial_seam_fault_np2.cpp, so the fault face is
// INTERIOR at np=1 and SHARED at np=2 — exactly the comparison we need.
//
// This test changes NOTHING in dynamic/ or drivers/.  It only reads.

#include "mfem.hpp"

#include "../../dynamic/wave_operator.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../fault/fault_basis.hpp"
#include "shared_fault_fixtures.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using mfem::seas::test_fixtures::BuildTwoTetSharedFault;
using mfem::seas::test_fixtures::PartitionByYSign;

namespace
{

int g_rank = 0, g_nprocs = 1;
int g_num_tests = 0, g_num_passed = 0, g_num_failed = 0;

#define CENSUS_ASSERT(cond, msg)                                              \
   do {                                                                       \
      ++g_num_tests;                                                          \
      if (cond) { ++g_num_passed;                                             \
         if (g_rank == 0) { std::cout << "  PASSED: " << msg << "\n"; } }     \
      else { ++g_num_failed;                                                  \
         std::cerr << "  FAILED [rank " << g_rank << ", line " << __LINE__    \
                   << "]: " << msg << "\n"; }                                 \
   } while (0)

BoundaryConfig FaultBC()
{
   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};
   return bc;
}

// Normalise in place; returns false if the vector is degenerate.
bool Normalize3(real_t v[3])
{
   const real_t n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
   if (!(n > 1e-300)) { return false; }
   v[0] /= n; v[1] /= n; v[2] /= n;
   return true;
}

real_t Dot3(const real_t a[3], const real_t b[3])
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

// Census result for one rank.
struct Census
{
   int  n_interior_faces = 0;
   int  n_shared_faces   = 0;
   int  n_qps_examined   = 0;
   int  c1_disagreements = 0;      ///< count(sign_flipped != !elem1_on_plus)
   // Face-averaged canonical normal of the FIRST shared fault face, built each way.
   real_t can_n_from_sign_flipped[3]  = {0, 0, 0};
   real_t can_n_from_elem1_on_plus[3] = {0, 0, 0};
   bool   have_shared = false;
};

// Walk every fault QP, comparing the two frame bits and recording the two
// candidate canonical normals on the first shared face.
Census RunCensus(WaveOperator<ParMesh> &op)
{
   Census c;
   const FaultBasis *fb = op.GetFaultBasis();
   MFEM_VERIFY(fb != nullptr, "census: WaveOperator has no FaultBasis");

   const Array<int> &int_faces = op.GetFaultInteriorFaces();
   const Array<int> &shr_faces = op.GetFaultSharedFaces();
   const std::vector<bool> &e1p_int = op.GetInteriorFaultElem1OnPlus();
   const std::vector<bool> &e1p_shr = op.GetSharedFaultElem1OnPlus();
   const std::map<int, int> &face2basis = op.GetFaultInteriorFaceToBasisIdx();

   // R-003: FaultBasis is built INTERIOR-FIRST -- wave_operator.inl:358 calls
   // fault_basis_->Compute(mesh, fault_interior_faces_, ...) and then appends
   // the shared faces.  The `int_faces.Size() + si` index below depends on that
   // layout, and this fixture never exercises the offset (it has either 0
   // interior or 0 shared faces), so assert the convention explicitly.
   MFEM_VERIFY(fb->NumFaces() == int_faces.Size() + shr_faces.Size(),
               "census: FaultBasis face count (" << fb->NumFaces()
               << ") != interior (" << int_faces.Size() << ") + shared ("
               << shr_faces.Size() << ") -- the interior-first basis-index "
               "convention this census relies on has changed.");

   c.n_interior_faces = int_faces.Size();
   c.n_shared_faces   = shr_faces.Size();

   // ---- interior fault faces -------------------------------------------
   for (int i = 0; i < int_faces.Size(); ++i)
   {
      auto it = face2basis.find(int_faces[i]);
      MFEM_VERIFY(it != face2basis.end(), "census: interior face has no basis idx");
      const FaultBasisData &bd = fb->GetBasis(it->second);
      const bool e1p = e1p_int[i];
      for (size_t q = 0; q < bd.qp_data.size(); ++q)
      {
         const FaultBasisQPData &qpd = bd.qp_data[q];
         ++c.n_qps_examined;
         if (qpd.sign_flipped != (!e1p)) { ++c.c1_disagreements; }
      }
   }

   // ---- shared fault faces ---------------------------------------------
   // Basis-index convention (BuildPerDOFFaultTables, spatial_dyn_driver.cpp:372):
   // shared face `si` uses basis index `int_faces.Size() + si`.
   for (int si = 0; si < shr_faces.Size(); ++si)
   {
      const int basis_idx = int_faces.Size() + si;
      const FaultBasisData &bd = fb->GetBasis(basis_idx);
      const bool e1p = e1p_shr[si];

      real_t acc_sf[3] = {0, 0, 0};
      real_t acc_e1[3] = {0, 0, 0};
      for (size_t q = 0; q < bd.qp_data.size(); ++q)
      {
         const FaultBasisQPData &qpd = bd.qp_data[q];
         ++c.n_qps_examined;
         if (qpd.sign_flipped != (!e1p)) { ++c.c1_disagreements; }

         const real_t s_sf = qpd.sign_flipped ? -1.0 : 1.0;   // shared-site rule
         const real_t s_e1 = (!e1p)           ? -1.0 : 1.0;   // interior-site rule
         for (int d = 0; d < 3; ++d)
         {
            acc_sf[d] += s_sf * qpd.normal[d];
            acc_e1[d] += s_e1 * qpd.normal[d];
         }
      }
      if (si == 0 && !bd.qp_data.empty())
      {
         for (int d = 0; d < 3; ++d)
         {
            c.can_n_from_sign_flipped[d]  = acc_sf[d];
            c.can_n_from_elem1_on_plus[d] = acc_e1[d];
         }
         // R-002: evaluate BOTH.  `&&` would short-circuit and leave the second
         // vector unnormalized on a degenerate face -- precisely the pathology
         // this census exists to detect.
         const bool ok_sf = Normalize3(c.can_n_from_sign_flipped);
         const bool ok_e1 = Normalize3(c.can_n_from_elem1_on_plus);
         c.have_shared = ok_sf && ok_e1;
      }
   }
   return c;
}

void PrintCensus(const Census &c)
{
   std::printf("[rank %d] fault faces: interior=%d shared=%d | QPs=%d | "
               "C1 disagreements (sign_flipped != !elem1_on_plus) = %d\n",
               g_rank, c.n_interior_faces, c.n_shared_faces,
               c.n_qps_examined, c.c1_disagreements);
   if (c.have_shared)
   {
      std::printf("[rank %d] shared face 0: can_n(sign_flipped) = "
                  "(%+.6f, %+.6f, %+.6f)   can_n(!elem1_on_plus) = "
                  "(%+.6f, %+.6f, %+.6f)\n",
                  g_rank,
                  c.can_n_from_sign_flipped[0], c.can_n_from_sign_flipped[1],
                  c.can_n_from_sign_flipped[2],
                  c.can_n_from_elem1_on_plus[0], c.can_n_from_elem1_on_plus[1],
                  c.can_n_from_elem1_on_plus[2]);
   }
   std::fflush(stdout);
}

}  // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank    = Mpi::WorldRank();
   g_nprocs  = Mpi::WorldSize();

   // R-007: C2 is only defined for the 2-rank seam.  At np>2 the run would
   // silently skip it and exit 0 -- a vacuous pass.
   MFEM_VERIFY(g_nprocs == 1 || g_nprocs == 2,
               "census fixture is defined for np=1 (interior fault face) or "
               "np=2 (shared fault face); got np=" << g_nprocs);

   if (g_rank == 0)
   {
      std::cout << "\n=== Phase 0 census: fault canonical-frame bits (np="
                << g_nprocs << ") ===\n";
   }

   ParMesh pmesh = PartitionByYSign(BuildTwoTetSharedFault());
   const BoundaryConfig bc = FaultBC();
   const int order = 1;
   WaveOperator<ParMesh> op(pmesh, order,
                            /*lambda=*/32.04e9, /*mu=*/32.04e9, /*rho=*/2670.0,
                            bc);

   const Census c = RunCensus(op);
   for (int r = 0; r < g_nprocs; ++r)
   {
      if (r == g_rank) { PrintCensus(c); }
      MPI_Barrier(MPI_COMM_WORLD);
   }

   // ---------------------------------------------------------------------
   // C1 — do the two bits ever disagree?
   //
   // MEASURED RESULT (2026-07-10): they DO, and necessarily so.
   //   * at np=1 the single fault face is INTERIOR and the bits agree (0);
   //   * at np=2 the face is SHARED and the bits disagree on EXACTLY ONE rank
   //     (the rank whose local Elem1 sits on the canonical minus side).
   // This is the direct consequence of C2 below: `sign_flipped` is
   // rank-independent while `elem1_on_plus` is rank-local, so on one of the two
   // ranks they must differ.  The bits are NOT interchangeable.
   // ---------------------------------------------------------------------
   int c1_global = 0;
   MPI_Allreduce(&c.c1_disagreements, &c1_global, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "\n[C1] global sign_flipped != !elem1_on_plus disagreements: "
                << c1_global << "\n";
   }
   if (g_nprocs == 1)
   {
      CENSUS_ASSERT(c1_global == 0,
                    "C1: on an INTERIOR fault face the two bits agree "
                    "(both elements are rank-local, so no ambiguity)");
   }
   else
   {
      CENSUS_ASSERT(c1_global > 0,
                    "C1: on a SHARED fault face the two bits DISAGREE on one "
                    "rank -> they are NOT interchangeable");
   }

   // ---------------------------------------------------------------------
   // C2 — at np=2, do the two ranks build the SAME or OPPOSITE canonical
   // normal for the same physical shared fault face?  Measured both ways.
   // ---------------------------------------------------------------------
   if (g_nprocs == 2)
   {
      // Exactly one shared fault face is expected on each rank.
      CENSUS_ASSERT(c.n_shared_faces == 1 && c.n_interior_faces == 0,
                    "np=2 fixture: each rank sees exactly 1 SHARED fault face "
                    "and 0 interior ones");
      CENSUS_ASSERT(c.have_shared, "np=2: shared-face canonical normals built");

      // R-002: never run the dot-product assertions on a degenerate (zero /
      // unnormalized) face-averaged normal -- they would fail for the wrong
      // reason and mask the real diagnosis.
      if (!c.have_shared)
      {
         if (g_rank == 0)
         {
            std::cerr << "  [C2] SKIPPED: degenerate face-averaged canonical "
                         "normal on the shared fault face.\n";
         }
      }
      else
      {

      real_t send[6] = { c.can_n_from_sign_flipped[0],
                         c.can_n_from_sign_flipped[1],
                         c.can_n_from_sign_flipped[2],
                         c.can_n_from_elem1_on_plus[0],
                         c.can_n_from_elem1_on_plus[1],
                         c.can_n_from_elem1_on_plus[2] };
      real_t recv[12] = {0};
      MPI_Allgather(send, 6, MPITypeMap<real_t>::mpi_type,
                    recv, 6, MPITypeMap<real_t>::mpi_type, MPI_COMM_WORLD);

      const real_t *sf0 = recv + 0, *e10 = recv + 3;
      const real_t *sf1 = recv + 6, *e11 = recv + 9;
      const real_t dot_sf = Dot3(sf0, sf1);
      const real_t dot_e1 = Dot3(e10, e11);

      if (g_rank == 0)
      {
         std::cout << "\n[C2] cross-rank canonical-normal agreement on the "
                      "shared fault face:\n";
         std::printf("     dot(can_n_rank0, can_n_rank1) using sign_flipped   "
                     "= %+.12f\n", dot_sf);
         std::printf("     dot(can_n_rank0, can_n_rank1) using !elem1_on_plus "
                     "= %+.12f\n", dot_e1);
         std::cout << "     (+1 => ranks agree; -1 => OPPOSITE frames, which is "
                      "what R-1601's comment claims)\n";
      }

      // MEASURED RESULT (2026-07-10).  These assertions now pin the TRUTH the
      // census discovered, and stand as a regression guard against anyone
      // "unifying" the frame bit in the future:
      //
      //   sign_flipped   -> dot = +1  : RANK-INDEPENDENT.  Correct for shared.
      //   !elem1_on_plus -> dot = -1  : RANK-LOCAL.  Using it on a shared face
      //                                 would give the two ranks OPPOSITE
      //                                 canonical frames -> antiparallel
      //                                 traction -> the 4e28 blow-up that
      //                                 R-1601 was written to prevent.
      //
      // Conclusion: the six frame sites are NOT a bug.  Each face class uses
      // the only bit that is well-defined for it.
      CENSUS_ASSERT(std::abs(dot_sf - 1.0) < 1e-10,
                    "C2: sign_flipped gives the SAME canonical frame on both "
                    "ranks (rank-independent) -> correct for shared faces");
      CENSUS_ASSERT(std::abs(dot_e1 + 1.0) < 1e-10,
                    "C2: !elem1_on_plus gives OPPOSITE canonical frames across "
                    "ranks (rank-local) -> must NOT be used on shared faces");
      }  // if (c.have_shared)
   }

   // R-005: pass/test counts are per-rank; reduce them so the summary is not
   // self-contradictory when a failure lives on only one rank (as C1 does).
   int failed_global = 0, passed_global = 0, tests_global = 0;
   MPI_Allreduce(&g_num_failed, &failed_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_passed, &passed_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_tests,  &tests_global,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "  Phase 0 census: " << passed_global << " / " << tests_global
                << " passed, " << failed_global << " failed\n";
      std::cout << "========================================\n";
   }
   return (failed_global == 0) ? 0 : 1;
}
