// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_fault_contact.cpp — fault opening version B, phase
// B1 (plan 2026-11-05, Section 5.2 and 9: the opening as the third kernel
// slip component), operator-level checks on the miniature BP3 layout of
// test_2d_fault_layout.cpp (ShearBox, dip 60, p = 1 and 2):
//
//   B1a  g = 0 equivalence: the stage operators built with slip_comps = 3
//        give, for the same tangential slip, the same assembled stage RHS
//        as the build with slip_comps = 2 — exactly (the g entries of the
//        slip are 0.0, so every extra assembly term adds an exact zero).
//        The solved stage traction then agrees to within the measured
//        instance-to-instance noise of the direct solver: each operator
//        owns its own MUMPS factorization, and two factorizations of the
//        same matrix are not bit-reproducible (alignment-dependent BLAS
//        kernels), so a second slip_comps = 2 operator provides the noise
//        baseline, with a 1e-12 relative floor.  (Observed on np = 2,
//        p = 1: RHS exactly equal, u differing by 2 ulp between solver
//        instances, amplified to ~3e-8 Pa by the Nitsche penalty.)
//   B1b  the dense Green's function built on the slip_comps = 3 operators
//        (its own internal check against the solve path runs the dense
//        column <-> stage slip index mapping) equals the slip_comps = 2
//        G stage for the same slip, within the same instance-noise bound.
//   B1c  the opening column: a unit PHYSICAL opening (faces apart; the
//        embedding is -g x basis normal, orientation measured by B1f) at
//        one fault DOF, solved through the stage (t = 0, no loading),
//        RAISES the elastic normal traction at that DOF (compression
//        positive: inserting an opening between the Dirichlet walls loads
//        the fault, so the diagonal of G^{sigma g} is positive), and its
//        out-of-plane shear row is exactly zero (plane strain).
//   B1d  energy reciprocity: with u_e = K^{-1} B e the solve of one unit
//        column, u_i^T (B e_j) = u_j^T (B e_i) for a slip column i and an
//        opening column j, because K is symmetric (SIPG) — exact up to the
//        direct solver.  (The ENTRYWISE identity G^{sigma delta}(j, i) =
//        G^{tau g}(i, j) does NOT hold: the traction projection weights by
//        the face mass, so reciprocity holds in the energy pairing, not
//        entry by entry.  T2 of the plan measures the entrywise asymmetry
//        of the gathered G_OO block where the active set uses it.)
//   B1e  a uniform opening of the whole fault loads its interior: between
//        Dirichlet walls, prescribing the faces apart inserts material, so
//        the elastic normal-traction response is positive (more
//        compression) at every DOF in the middle 60% of the fault, and in
//        the mean.  (A tip row can unload — crack mechanics, not a bug.)
//   B1f  the kinematic sign: the solved uniform-gap displacement moves the
//        two sides of every fault face APART along the geometric normal
//        ((u_A - u_B) . n > 0 with side A along +n) — the measurement that
//        fixes the embedding orientation (and the s of plan eq. 4.3).
//
// Usage (from miniapps/seas; PETSc/MUMPS needed):
//   mpirun -np N ./seas_test_fault_contact

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_frame_field.hpp"
#include "../../constitutive/linear_elastic.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kLambda = 32.04e9, kMu = 32.04e9;
static const real_t kLs = 20e3, kD = 20e3, kWf = 12e3;
static const int kNs = 8, kNy = 10;

static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static real_t RSum(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_SUM, c); return v; }

// The BP3 layout in miniature (test_2d_fault_layout.cpp): a grid (s, y)
// mapped by x = s - y cot(dip); the grid line s = 0 is the fault, attr 3
// for depth -y < Wf, attr 5 below; top attr 1, other sides attr 5.
static Mesh *ShearBox(real_t dip_deg)
{
   const real_t c = (dip_deg == 90.0) ? 0.0 : 1.0 / std::tan(dip_deg * M_PI / 180.0);
   const int nvx = kNs + 1, nvy = kNy + 1;
   Mesh *m = new Mesh(2, nvx * nvy, 2 * kNs * kNy, 0, 2);
   for (int j = 0; j < nvy; j++)
      for (int i = 0; i < nvx; i++)
      {
         const real_t s = -kLs + 2.0 * kLs * i / kNs, y = -kD + kD * j / kNy;
         m->AddVertex(s - y * c, y);
      }
   auto v = [&](int i, int j) { return i + nvx * j; };
   for (int j = 0; j < kNy; j++)
      for (int i = 0; i < kNs; i++)
      {
         m->AddTriangle(v(i, j), v(i + 1, j), v(i + 1, j + 1), 1);
         m->AddTriangle(v(i, j), v(i + 1, j + 1), v(i, j + 1), 1);
      }
   for (int i = 0; i < kNs; i++)
   {
      m->AddBdrSegment(v(i, 0), v(i + 1, 0), 5);
      m->AddBdrSegment(v(i, kNy), v(i + 1, kNy), 1);
   }
   for (int j = 0; j < kNy; j++)
   {
      m->AddBdrSegment(v(0, j), v(0, j + 1), 5);
      m->AddBdrSegment(v(kNs, j), v(kNs, j + 1), 5);
      const real_t ymid = -kD + kD * (j + 0.5) / kNy;
      m->AddBdrSegment(v(kNs / 2, j), v(kNs / 2, j + 1), (-ymid < kWf) ? 3 : 5);
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

static std::shared_ptr<const FaultFrameConfig> FrameConfig(const real_t ref[3],
                                                           const real_t up[3])
{
   auto c = std::make_shared<FaultFrameConfig>();
   for (int d = 0; d < 3; d++)
   {
      c->orientation.ref_normal[d] = ref[d];
      c->orientation.up[d] = up[d];
      c->up[d] = up[d];
   }
   c->orientation.fault_attrs = {3};
   return c;
}

// The plate loading of the miniature (test_2d_fault_layout.cpp Q1): linear
// in t, so g(0) = 0 and the Green build accepts it.
static void PlateLoading2D(const Vector &x, real_t t, Vector &u)
{
   const real_t n[2] = {std::sqrt(3.0) / 2.0, 0.5}, sd[2] = {-0.5, std::sqrt(3.0) / 2.0};
   const real_t Vp = 1e-9, d = x(0) * n[0] + x(1) * n[1];
   const real_t Vh = (d > 1.0) ? 0.5 * Vp * t : (d < -1.0) ? -0.5 * Vp * t : Vp * t;
   u.SetSize(2);
   u(0) = Vh * sd[0];
   u(1) = Vh * sd[1];
}

using Op = ElasticityDomainOperator<ParMesh>;

static std::unique_ptr<Op> MakeOp(ParMesh &pm, int order)
{
   static const real_t ref[3] = {1, 0, 0}, up[3] = {0, 1, 0};
   LinearElastic *model = new LinearElastic(kLambda, kMu);   // leaked: test lifetime
   DomainConfig dc;
   FaultFrameSource src;
   src.config = FrameConfig(ref, up);
   dc.fault_frames = src;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = PlateLoading2D;
   return std::make_unique<Op>(pm, order, *model, 0.0, kWf, 2.0 * kLs, bc,
                               DGMethod::IP, SolverType::MUMPS, dc);
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   g_rank = Mpi::WorldRank();

   const real_t t = 100.0 * 365.25 * 86400.0;
   for (int p : {1, 2})
   {
      if (g_rank == 0) { std::cout << "\n=== B1: opening slip component, p = " << p << " ===\n"; }
      std::unique_ptr<Mesh> ser(ShearBox(60.0));
      ParMesh pm2(comm, *ser), pm3(comm, *ser);
      auto op2 = MakeOp(pm2, p);
      auto op3 = MakeOp(pm3, p);
      op2->BuildQDStageOperators(true);        // today's build (slip_comps 2)
      op3->BuildQDStageOperators(true, 3);     // B1: with the opening column
      const int nown = op2->GetNumOwnedFaultDOFs();
      MFEM_VERIFY(op3->GetNumOwnedFaultDOFs() == nown, "owned DOF mismatch");
      const int slot = op2->FrameField()->InPlaneSlot();   // 0 (dip) here
      const int out = 1 - slot;
      long long first = 0, mine = nown;
      MPI_Exscan(&mine, &first, 1, MPI_LONG_LONG, MPI_SUM, comm);
      if (g_rank == 0) { first = 0; }
      long long nf = first + nown;
      MPI_Allreduce(MPI_IN_PLACE, &nf, 1, MPI_LONG_LONG, MPI_MAX, comm);

      // ---- B1a: the same tangential slip through both builds ----
      Vector slip2(2 * nown), slip3(3 * nown);
      slip3 = 0.0;
      for (int i = 0; i < nown; i++)
      {
         const real_t a = 1.0 + 0.5 * std::sin(0.37 * i + 0.1 + g_rank);
         slip2(2 * i + slot) = a;
         slip2(2 * i + out) = 0.0;
         slip3(3 * i + slot) = a;
      }
      // The assembled stage RHS: exact, no solver in the loop.
      Vector b2a, b3a;
      op2->QDStageRHS(t, slip2, b2a);
      op3->QDStageRHS(t, slip3, b3a);
      real_t db = 0.0;
      MFEM_VERIFY(b2a.Size() == b3a.Size(), "stage RHS size mismatch");
      for (int i = 0; i < b2a.Size(); i++) { db = std::max(db, std::abs(b2a(i) - b3a(i))); }
      db = RMax(db, comm);
      TEST_ASSERT(db == 0.0, "B1a p = " << p << ": slip_comps 3 stage RHS = "
                  "slip_comps 2 stage RHS for g = 0, exactly (" << db << ")");
      // The solved stage traction: bounded by the direct solver's
      // instance-to-instance noise, measured by a second slip_comps = 2
      // operator built the same way.
      ParMesh pm2b(comm, *ser);
      auto op2b = MakeOp(pm2b, p);
      op2b->BuildQDStageOperators(true);
      ParGridFunction u2(&op2->GetFESpace()), u3(&op3->GetFESpace()),
                      u2b(&op2b->GetFESpace());
      Vector T2, Tn2, T3, Tn3, T2b, Tn2b;
      op2->SolveStage(t, slip2, u2, T2, &Tn2);
      op3->SolveStage(t, slip3, u3, T3, &Tn3);
      op2b->SolveStage(t, slip2, u2b, T2b, &Tn2b);
      real_t ds = 0.0, dn = 0.0, ns = 0.0, nn = 0.0, Tsc = 0.0;
      for (int i = 0; i < 2 * nown; i++)
      {
         ds = std::max(ds, std::abs(T2(i) - T3(i)));
         ns = std::max(ns, std::abs(T2(i) - T2b(i)));
         Tsc = std::max(Tsc, std::abs(T2(i)));
      }
      for (int i = 0; i < nown; i++)
      {
         dn = std::max(dn, std::abs(Tn2(i) - Tn3(i)));
         nn = std::max(nn, std::abs(Tn2(i) - Tn2b(i)));
         Tsc = std::max(Tsc, std::abs(Tn2(i)));
      }
      ds = RMax(ds, comm);
      dn = RMax(dn, comm);
      ns = RMax(ns, comm);
      nn = RMax(nn, comm);
      Tsc = RMax(Tsc, comm);
      const real_t noise = std::max(std::max(ns, nn), 1e-12 * Tsc);
      TEST_ASSERT(ds <= noise && dn <= noise, "B1a p = " << p << ": slip_comps 3 "
                  "stage traction = slip_comps 2 for g = 0 within solver instance "
                  "noise (shear " << ds << ", normal " << dn << ", measured noise "
                  << std::max(ns, nn) << ", bound " << noise << ", scale " << Tsc << " Pa)");

      // ---- B1b: the dense G built on the slip_comps = 3 operators ----
      // BuildQDGreen's own check (1e-8 against the solve path) runs the
      // dense-column -> stage-entry mapping; then the stage product must
      // equal the slip_comps = 2 build's exactly.
      op2->BuildQDGreen(QDGreenConfig{});
      op3->BuildQDGreen(QDGreenConfig{});
      Vector Tg2, Tng2, Tg3, Tng3;
      op2->GreenStage(t, slip2, Tg2, &Tng2);
      op3->GreenStage(t, slip2, Tg3, &Tng3);   // G keeps the 2-per-DOF slip space
      ds = 0.0;
      dn = 0.0;
      for (int i = 0; i < 2 * nown; i++) { ds = std::max(ds, std::abs(Tg2(i) - Tg3(i))); }
      for (int i = 0; i < nown; i++) { dn = std::max(dn, std::abs(Tng2(i) - Tng3(i))); }
      ds = RMax(ds, comm);
      dn = RMax(dn, comm);
      // G's columns come from each operator's own solver instance, so the
      // bound is the B1a instance noise (measured 0.0 here on every
      // np/p combination so far — the blocked multi-RHS solves reproduce).
      TEST_ASSERT(ds <= noise && dn <= noise, "B1b p = " << p << ": G of the "
                  "slip_comps 3 build = G of the slip_comps 2 build within solver "
                  "instance noise (shear " << ds << ", normal " << dn
                  << ", bound " << noise << ")");

      // ---- B1c: the opening column at one DOF (t = 0: no loading) ----
      const long long jg = nf / 2, ig = nf / 3;
      auto unit = [&](long long g, int c, Vector &e)
      {
         e.SetSize(3 * nown);
         e = 0.0;
         if (g >= first && g < first + nown) { e(3 * static_cast<int>(g - first) + c) = 1.0; }
      };
      auto at = [&](const Vector &v, long long g, int stride, int c)
      {
         real_t val = 0.0;
         if (g >= first && g < first + nown)
         {
            val = v(stride * static_cast<int>(g - first) + c);
         }
         return RSum(val, comm);
      };
      Vector eg, Tcol, Tncol;
      ParGridFunction ug(&op3->GetFESpace());
      unit(jg, 2, eg);
      op3->SolveStage(0.0, eg, ug, Tcol, &Tncol);
      const real_t gdiag = at(Tncol, jg, 1, 0);
      const real_t gout = RMax(std::abs(at(Tcol, jg, 2, out)), comm);
      real_t outrow = 0.0;
      for (int i = 0; i < nown; i++) { outrow = std::max(outrow, std::abs(Tcol(2 * i + out))); }
      outrow = RMax(outrow, comm);
      TEST_ASSERT(gdiag > 0.0, "B1c p = " << p << ": a unit physical opening raises "
                  "the elastic normal traction at its own DOF (G^{sigma g} diagonal "
                  << gdiag << " Pa/m)");
      TEST_ASSERT(outrow == 0.0, "B1c p = " << p << ": the opening column's out-of-plane "
                  "shear rows are exactly zero (" << outrow << "; probe slot value " << gout << ")");

      // ---- B1d: energy reciprocity u_i^T B e_j = u_j^T B e_i ----
      // u_e = K^{-1} B e (t = 0: no loading), so u_i^T B e_j = u_i^T K u_j,
      // symmetric because K is SIPG-symmetric.  Exact up to the direct
      // solver, and it runs the g column of B through the same energy
      // pairing as the slip columns.
      Vector ei, Ti, Tni, bi, bj;
      ParGridFunction ui(&op3->GetFESpace());
      unit(ig, slot, ei);
      op3->SolveStage(0.0, ei, ui, Ti, &Tni);
      op3->QDStageRHS(0.0, ei, bi);
      op3->QDStageRHS(0.0, eg, bj);
      auto pdot = [&](const Vector &a, const Vector &b)
      {
         real_t d = 0.0;
         for (int i = 0; i < a.Size(); i++) { d += a(i) * b(i); }
         return RSum(d, comm);
      };
      const real_t w_ij = pdot(ui, bj);   // u_i^T K u_j
      const real_t w_ji = pdot(ug, bi);   // u_j^T K u_i
      const real_t wscale = std::max(std::abs(w_ij), std::abs(w_ji));
      const real_t wrel = wscale > 0.0 ? std::abs(w_ij - w_ji) / wscale : 0.0;
      TEST_ASSERT(wscale > 0.0 && wrel <= 1e-9, "B1d p = " << p << ": energy "
                  "reciprocity of a slip and an opening column, u_i^T B e_j = "
                  "u_j^T B e_i (rel. diff " << wrel << ", value " << w_ij << ")");

      // ---- B1e: a uniform opening relieves compression on the interior ----
      Vector eu(3 * nown), Tu, Tnu;
      ParGridFunction uu(&op3->GetFESpace());
      eu = 0.0;
      for (int i = 0; i < nown; i++) { eu(3 * i + 2) = 1.0e-3; }
      op3->SolveStage(0.0, eu, uu, Tu, &Tnu);
      Vector xyz_loc, xyz_own;
      op3->GetFaultDOFCoords3D(xyz_loc);
      op3->RestrictToOwnedFault(xyz_loc, xyz_own, 3);
      real_t tn_mid = 1e300, tn_mean = 0.0;
      for (int i = 0; i < nown; i++)
      {
         tn_mean += Tnu(i);
         const real_t depth = -xyz_own(3 * i + 1);   // y <= 0 on the fault
         if (depth > 0.2 * kWf && depth < 0.8 * kWf)
         {
            tn_mid = std::min(tn_mid, Tnu(i));
         }
      }
      tn_mid = -RMax(-tn_mid, comm);
      tn_mean = RSum(tn_mean, comm) / static_cast<real_t>(nf);
      TEST_ASSERT(tn_mid > 0.0 && tn_mean > 0.0, "B1e p = " << p << ": uniform 1 mm "
                  "opening loads the middle 60% of the fault (interior min "
                  << tn_mid << " Pa, mean " << tn_mean << " Pa)");

      // ---- B1f: the KINEMATIC direction of the g column (B2 sign s) ----
      // uu solves the uniform +1 mm gap.  At the centre of each local
      // interior fault face, the jump (u_A - u_B) . n_hat, with side A the
      // element whose centroid lies along +n_hat of the face centre and
      // n_hat = (sin dip, cos dip, 0) the geometric fault normal: POSITIVE
      // = +g moves the faces APART, the physical opening.  (Measured
      // 2026-11-05: the raw +basis-normal jump moves them TOGETHER by
      // exactly g, which is why EmbedOpening embeds -g x normal.)
      {
         const real_t rad = 60.0 * M_PI / 180.0;
         const real_t nh[2] = {std::sin(rad), std::cos(rad)};
         real_t sep_min = 1e300, sep_max = -1e300;
         const auto &ff = op3->GetFaultInteriorFaces();
         for (int k = 0; k < ff.Size(); k++)
         {
            auto *FTr = pm3.GetInteriorFaceTransformations(ff[k]);
            if (!FTr) { continue; }
            const IntegrationPoint &ipc = Geometries.GetCenter(FTr->GetGeometryType());
            FTr->SetAllIntPoints(&ipc);
            Vector u1v, u2v, fc(2), c1(2), c2(2);
            FTr->Face->Transform(ipc, fc);
            uu.GetVectorValue(*FTr->Elem1, FTr->GetElement1IntPoint(), u1v);
            uu.GetVectorValue(*FTr->Elem2, FTr->GetElement2IntPoint(), u2v);
            pm3.GetElementCenter(FTr->Elem1No, c1);
            pm3.GetElementCenter(FTr->Elem2No, c2);
            const real_t side1 = (c1(0) - fc(0)) * nh[0] + (c1(1) - fc(1)) * nh[1];
            const real_t side2 = (c2(0) - fc(0)) * nh[0] + (c2(1) - fc(1)) * nh[1];
            MFEM_VERIFY(side1 * side2 < 0.0, "B1f: sides not straddling");
            const real_t ju = ((u1v(0) - u2v(0)) * nh[0] + (u1v(1) - u2v(1)) * nh[1])
                              * (side1 > 0.0 ? 1.0 : -1.0);
            sep_min = std::min(sep_min, ju);
            sep_max = std::max(sep_max, ju);
         }
         sep_max = RMax(sep_max, comm);
         sep_min = -RMax(-sep_min, comm);
         if (g_rank == 0)
         {
            std::cout << "  B1f p = " << p << ": face separation (u_A - u_B) . n per"
                      " +1 mm gap in [" << sep_min << ", " << sep_max << "] m ("
                      << (sep_min > 0.0 ? "+g = faces APART"
                          : (sep_max < 0.0 ? "+g = faces TOGETHER" : "MIXED")) << ")\n";
         }
         TEST_ASSERT(sep_min > 0.0, "B1f p = " << p << ": a positive gap moves the "
                     "two sides of every fault face APART (kinematic sign of the g "
                     "column: separation in [" << sep_min << ", " << sep_max << "] m)");
      }
   }

   real_t failed = num_failed;
   failed = RMax(failed, comm);
   if (g_rank == 0)
   {
      std::cout << "\n" << (failed > 0 ? "FAILED" : "PASSED") << " ("
                << num_tests << " checks, " << static_cast<int>(failed)
                << " failed)\n";
   }
   return failed > 0 ? 1 : 0;
}
