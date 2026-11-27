// tests/ucerf3d/test_loading_3d.cpp - UCERF 3-D plan rev5, P1.6: the loading
// options A and B on a 3-D mesh (tests U-L1, U-L2).
//
// Option A ([cycle.loading] creep_attrs): the walls move by -+ Vp t / 2 and
// every creep face (prescribed slip, attribute 6) slips at Vp along the plate
// direction projected onto the face (cycle/loading.hpp, MakeCreepFaceLoading).
// Option B: the same deck without creep_attrs, walls only.
//
// The deck tests/ucerf3d/configs/T-bend_loading_A.toml on T-bend: the fault
// and the creep surface below it are two vertical planes, y = 0 for x < 0
// and turned by 20 degrees about z for x > 0.  With ref_normal (0, -1, 0), m
// = (0, 1, 0) and s = (1, 0, 0) the analytic values are
//    west:  n = (0, -1, 0),            jump  Vp t (1, 0, 0)
//    east:  n = (sin 20, -cos 20, 0),  jump  Vp t (cos 20, sin 20, 0)
// (jump = u(side opposite to n) - u(side n points to) = u(y > plane) -
// u(y < plane): the side m points to moves by +s against the other).
//
// U-L1 (option A)
//   A1  the creep census: every creep face's oriented normal is the
//       analytic one (1e-12) and the creep faces number the mesh's
//   A2  MakeCreepFaceLoading at every creep face: the analytic jump, of
//       size Vp t and in the face's plane (1e-12)
//   A3  the QD operator's Dirichlet load vector equals that of a reference
//       operator whose creep attribute takes a function of position (the
//       analytic jump by the sign of x), 1e-12 relative
//   A4  the same through the QD stage operators (D g(t))
//   A5  the dynamic phase's operator gives the same load at t0 + tau
//   A6  a static solve with the fault locked: the displacement jump across
//       the creep faces is the prescribed one (imposed weakly: the mean
//       within --jump-tol of Vp t along the projected direction)
// U-L2 (option B)
//   B1  no prescribed-slip face in the census
//   B2  the load vector equals that of an operator with the walls only and
//       the BP5 function (MakeBP5DirichletFunc), 1e-12, and differs from
//       option A's
//   B3  the static solve has no jump across the attribute-6 faces
//
// Usage (from miniapps/seas; PETSc needed):
//   mpirun -np N ./seas_test_loading_3d [--deck tests/ucerf3d/configs/T-bend_loading_A.toml]
// --split --no-solve: the mesh is cut along the fault and creep planes, so
// their faces are shared between ranks (the default partition shares none
// on T-bend); the solves A6 and B3 are left out (they measure the jump
// where both elements are on one rank).
// L3D_DUMP=DIR writes the creep-face jumps of A6 and B3 (rank files).

#include "mfem.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/loading.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../common/mpi_context.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace mfem;
using namespace mfem::seas;

namespace
{
int num_tests = 0, num_failed = 0;
void Check(MPIContext &mpi, bool ok, const std::string &what)
{
   num_tests++;
   if (!ok) { num_failed++; }
   if (mpi.IsRoot()) { std::cout << (ok ? "  PASSED: " : "  FAILED: ") << what << std::endl; }
}
real_t GMax(MPIContext &mpi, real_t v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, mpi.GetComm());
   return v;
}
real_t GSum(MPIContext &mpi, real_t v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, mpi.GetComm());
   return v;
}
real_t Dot3(const real_t *a, const real_t *b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

typedef ElasticityDomainOperator<ParMesh> DomOp;

/// The largest |a - b| over the largest |b|, over the ranks.
real_t RelDiff(MPIContext &mpi, const Vector &a, const Vector &b)
{
   MFEM_VERIFY(a.Size() == b.Size(), "RelDiff: sizes " << a.Size() << ", " << b.Size());
   real_t d = 0.0, m = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      d = std::max(d, std::abs(a(i) - b(i)));
      m = std::max(m, std::abs(b(i)));
   }
   d = GMax(mpi, d);
   m = GMax(mpi, m);
   return m > 0.0 ? d / m : (d > 0.0 ? 1.0 : 0.0);
}

/// The displacement jump u(side opposite to n) - u(side n points to) at the
/// centroid of every face with boundary attribute `attr` whose two elements
/// are on this rank; n: nrm(x) (unit).  Rows: x y z, jump (3), n (3).
std::vector<std::array<real_t, 9>> FaceJumps(
   ParMesh &pm, const ParGridFunction &u, int attr,
   const std::function<void(int face, const Vector &x, real_t *n)> &nrm)
{
   std::vector<std::array<real_t, 9>> rows;
   Vector xc(3), v1(3), v2(3), c1(3);
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      if (pm.GetBdrAttribute(be) != attr) { continue; }
      const int f = pm.GetBdrElementFaceIndex(be);
      int e1, e2;
      pm.GetFaceElements(f, &e1, &e2);
      if (e2 < 0) { continue; }   // the other element is on another rank
      ElementTransformation *FT = pm.GetFaceTransformation(f);
      FT->Transform(Geometries.GetCenter(pm.GetFaceGeometry(f)), xc);
      IntegrationPoint ip;
      pm.GetElementTransformation(e1)->TransformBack(xc, ip);
      u.GetVectorValue(e1, ip, v1);
      pm.GetElementTransformation(e1)->Transform(Geometries.GetCenter(pm.GetElementGeometry(e1)), c1);
      pm.GetElementTransformation(e2)->TransformBack(xc, ip);
      u.GetVectorValue(e2, ip, v2);
      real_t n[3];
      nrm(f, xc, n);
      // e1 on the side n points to?
      const real_t side = (c1(0) - xc(0)) * n[0] + (c1(1) - xc(1)) * n[1] + (c1(2) - xc(2)) * n[2];
      std::array<real_t, 9> r;
      for (int i = 0; i < 3; i++)
      {
         r[i] = xc(i);
         r[3 + i] = (side > 0.0) ? v2(i) - v1(i) : v1(i) - v2(i);
         r[6 + i] = n[i];
      }
      rows.push_back(r);
   }
   return rows;
}

void Dump(const char *dir, const char *name, int rank, const std::vector<std::array<real_t, 9>> &rows)
{
   std::ostringstream fn;
   fn << dir << "/" << name << "." << rank << ".txt";
   std::ofstream f(fn.str());
   f << std::setprecision(10);
   for (const auto &r : rows)
   {
      for (int i = 0; i < 9; i++) { f << r[i] << (i == 8 ? "\n" : " "); }
   }
}
}  // namespace

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   const char *deck = "tests/ucerf3d/configs/T-bend_loading_A.toml";
   real_t bend_deg = 20.0, jump_tol = 0.1, t = 3.0e9;
   int creep_attr = 6;
   bool split = false, solve = true;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "The option A deck.");
   args.AddOption(&bend_deg, "-bend", "--bend", "Turn of the plane x > 0 about z (degrees).");
   args.AddOption(&jump_tol, "-jt", "--jump-tol", "A6: tolerance of the mean jump over Vp t.");
   args.AddOption(&t, "-t", "--time", "The loading time of the checks (s).");
   args.AddOption(&split, "-sp", "--split", "-nsp", "--no-split",
                  "Partition along the fault and creep planes (their faces become shared).");
   args.AddOption(&solve, "-so", "--solve", "-nso", "--no-solve",
                  "Run the static solves A6 and B3 (they see faces with both elements on "
                  "one rank only).");
   args.Parse();
   if (!args.Good()) { if (mpi.IsRoot()) { args.PrintUsage(std::cout); } return 1; }
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cerr << "test_loading_3d needs PETSc.\n"; }
   return 2;
#else
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
   MFEMInitializePetsc(&argc, &argv, cfg.qd.petsc_options.c_str(), NULL);
   {
      const char *dump = std::getenv("L3D_DUMP");
      const cycle::LoadingConfig &L = cfg.loading;
      MFEM_VERIFY(L.Plate() && L.creep_attrs == std::vector<int>{creep_attr},
                  deck << ": [cycle.loading] kind = \"plate\" with creep_attrs = [6] is needed");
      const real_t Vp = cfg.plate_rate, Vt = Vp * t;
      const real_t th = bend_deg * M_PI / 180.0, cb = std::cos(th), sb = std::sin(th);
      const real_t n_w[3] = {0.0, -1.0, 0.0}, n_e[3] = {sb, -cb, 0.0};
      const real_t j_w[3] = {Vt, 0.0, 0.0}, j_e[3] = {Vt * cb, Vt * sb, 0.0};
      if (mpi.IsRoot())
      {
         std::cout << "test_loading_3d (UCERF 3-D plan rev5, P1.6, U-L1 and U-L2) on " << mpi.Size()
                   << " rank(s)" << (split ? ", split along the fault and creep planes" : "")
                   << ": " << deck << "\n  Vp = " << Vp << " m/s, t = " << t
                   << " s, Vp t = " << Vt << " m; bend " << bend_deg << " deg" << std::endl;
      }
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      long long n_creep_mesh = 0;
      for (int be = 0; be < smesh->GetNBE(); be++)
      {
         if (smesh->GetBdrAttribute(be) == creep_attr) { n_creep_mesh++; }
      }
      // --split: cut the mesh along the fault and creep planes (and at the
      // bend), so that their faces are shared between ranks
      Array<int> part;
      if (split)
      {
         part.SetSize(smesh->GetNE());
         Vector c(3);
         for (int e = 0; e < smesh->GetNE(); e++)
         {
            smesh->GetElementTransformation(e)->Transform(
               Geometries.GetCenter(smesh->GetElementGeometry(e)), c);
            const int east = c(0) > 0.0 ? 1 : 0;
            const int north = (c(1) - (east ? c(0) * sb / cb : 0.0)) > 0.0 ? 1 : 0;
            part[e] = (north + 2 * east) % mpi.Size();
         }
      }
      ParMesh pmesh(mpi.GetComm(), *smesh, split ? part.GetData() : nullptr);
      smesh.reset();
      const spatial::SpatialFrictionConfig &s = cfg.shared;

      // ================= U-L1: option A =================
      if (mpi.IsRoot()) { std::cout << "\n=== U-L1 option A: walls and creep faces ===" << std::endl; }
      Vector bA;
      LinearElastic model(s.material_fallback.lambda, s.material_fallback.mu);
      std::ostringstream os;
      {   // option A's operators (and their factors) end with this block
      cycle::QDPhase qd(pmesh, mpi, cfg);
      DomOp &dom = qd.Domain();
      std::shared_ptr<const FaultFrameField> frames = dom.FrameField();
      MFEM_VERIFY(frames, "the QD operator has no fault frame field");

      // A1, A2: the census normals and the face function
      const DirichletFaceFunc creep = cycle::MakeCreepFaceLoading(L, Vp);
      real_t e_n = 0.0, e_j = 0.0, e_size = 0.0, e_plane = 0.0;
      long long n_faces = 0, n_west = 0, n_east = 0;
      {
         Vector xc(3), u(3);
         for (int be = 0; be < pmesh.GetNBE(); be++)
         {
            if (pmesh.GetBdrAttribute(be) != creep_attr) { continue; }
            const int f = pmesh.GetBdrElementFaceIndex(be);
            MFEM_VERIFY(frames->IsCreepFace(f), "an attribute-6 face is not in the creep census");
            pmesh.GetFaceTransformation(f)->Transform(
               Geometries.GetCenter(pmesh.GetFaceGeometry(f)), xc);
            const bool west = xc(0) < 0.0;
            const real_t *na = west ? n_w : n_e, *ja = west ? j_w : j_e;
            const real_t *n = frames->CreepNormal(f);
            creep(xc, t, n, u);
            real_t dn = 0.0, dj = 0.0;
            for (int i = 0; i < 3; i++)
            {
               dn += (n[i] - na[i]) * (n[i] - na[i]);
               dj += (u(i) - ja[i]) * (u(i) - ja[i]);
            }
            e_n = std::max(e_n, std::sqrt(dn));
            e_j = std::max(e_j, std::sqrt(dj) / Vt);
            e_size = std::max(e_size, std::abs(u.Norml2() - Vt) / Vt);
            e_plane = std::max(e_plane, std::abs(u(0) * n[0] + u(1) * n[1] + u(2) * n[2]) / Vt);
            n_faces++;
            (west ? n_west : n_east)++;
         }
      }
      e_n = GMax(mpi, e_n); e_j = GMax(mpi, e_j); e_size = GMax(mpi, e_size);
      e_plane = GMax(mpi, e_plane);
      const long long census = frames->Report().creep.faces;
      {
         long long cnt[3] = {n_faces, n_west, n_east};
         MPI_Allreduce(MPI_IN_PLACE, cnt, 3, MPI_LONG_LONG, MPI_SUM, mpi.GetComm());
         n_faces = cnt[0]; n_west = cnt[1]; n_east = cnt[2];
      }
      os << "U-L1 A1: the census holds " << census << " creep faces (mesh: " << n_creep_mesh
         << "); normals are (0, -1, 0) west and (sin 20, -cos 20, 0) east of the bend ("
         << e_n << ", 1e-12; " << n_west << " west and " << n_east
         << " east face visits over the ranks)";
      Check(mpi, census == n_creep_mesh && n_west > 0 && n_east > 0 && e_n <= 1e-12, os.str());
      os.str("");
      os << "U-L1 A2: the face function gives Vp t (1, 0, 0) west and Vp t (cos 20, sin 20, 0) "
         << "east (" << e_j << "), of size Vp t (" << e_size << ") and in the face (" << e_plane
         << "), 1e-12";
      Check(mpi, e_j <= 1e-12 && e_size <= 1e-12 && e_plane <= 1e-12, os.str());

      // A3: the load vector against a reference operator (creep by position)
      DomainConfig dcfg;
      dcfg.penalty_factor = cfg.penalty_factor;
      dcfg.fault_frames = FaultFrameSource{frames, nullptr};
      BoundaryConfig rbc;
      rbc.dirichlet_attrs = {5, creep_attr};
      rbc.natural_attrs = {1};
      rbc.fault_attr = cfg.fault_attr;
      rbc.default_dirichlet_func = cycle::MakeCycleLoading(cfg);
      rbc.dirichlet_funcs[creep_attr] = [Vp, cb, sb](const Vector &x, real_t tt, Vector &u)
      {
         u.SetSize(3);
         u = 0.0;
         if (x(0) < 0.0) { u(0) = Vp * tt; }
         else { u(0) = Vp * tt * cb; u(1) = Vp * tt * sb; }
      };
      DomOp ref(pmesh, s.mesh.order, model, Vp, 0.0, 0.0, rbc, DGMethod::IP, SolverType::MUMPS,
                dcfg);
      Vector bR;
      dom.AssembleDirichletOnlyRHS(bA, t);
      ref.AssembleDirichletOnlyRHS(bR, t);
      const real_t e_b = RelDiff(mpi, bA, bR);
      os.str("");
      os << "U-L1 A3: the QD operator's Dirichlet load equals the reference operator's (creep "
         << "jump by position), " << e_b << " of the largest entry " << GMax(mpi, bR.Normlinf())
         << " (1e-12)";
      Check(mpi, GMax(mpi, bR.Normlinf()) > 0.0 && e_b <= 1e-12, os.str());

      // A4: through the stage operators
      {
         dom.BuildQDStageOperators(false);
         ref.BuildQDStageOperators(false);
         Vector zero(2 * dom.GetNumOwnedFaultDOFs()), sA, sR;
         zero = 0.0;
         dom.QDStageRHS(t, zero, sA);
         ref.QDStageRHS(t, zero, sR);
         const real_t e_s = RelDiff(mpi, sA, sR);
         os.str("");
         os << "U-L1 A4: the stage operators' D g(t) equals the reference's (" << e_s
            << " of the largest " << GMax(mpi, sR.Normlinf()) << ", 1e-12)";
         Check(mpi, GMax(mpi, sR.Normlinf()) > 0.0 && e_s <= 1e-12, os.str());
      }

      // A5: the dynamic phase's operator, at t0 + tau with t0 = 0
      {
         const MaterialField material = MaterialField::MakeConstant(
            s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
         const spatial::PorePressureSpec no_pp;
         const cycle::FaultParams par = cycle::BuildFaultParams(
            qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
         cycle::PrestressConfig pc = cfg.prestress;
         pc.a_max = cycle::GlobalMaxA(par, mpi.GetComm());
         const cycle::Prestress pre = cycle::BuildPrestress(par, pc);
         cycle::DynPhase dyn(pmesh, cfg, qd.Table(), par, pre, 1.0, nullptr, frames);
         Vector bD;
         dyn.Op().Elasticity().AssembleDirichletOnlyRHS(bD, t);
         const real_t e_d = RelDiff(mpi, bD, bR);
         os.str("");
         os << "U-L1 A5: the dynamic phase's operator gives the same load (" << e_d
            << ", 1e-12)";
         Check(mpi, e_d <= 1e-12, os.str());
      }

      // A6: the static solve with the fault locked
      if (solve)
      {
         ParGridFunction uA(&dom.GetFESpace());
         Vector zero(2 * dom.GetNumFaultDOFs());
         zero = 0.0;
         dom.Solve(t, zero, uA);
         const auto rows = FaceJumps(pmesh, uA, creep_attr,
                                     [&](int f, const Vector &, real_t *n)
         { const real_t *c = frames->CreepNormal(f); for (int i = 0; i < 3; i++) { n[i] = c[i]; } });
         if (dump) { Dump(dump, "jumps_A", mpi.Rank(), rows); }
         // along the analytic direction, and across it
         real_t sum = 0.0, cnt = 0.0, lo = 1e300, hi = -1e300, off = 0.0;
         for (const auto &r : rows)
         {
            const real_t *ja = (r[0] < 0.0) ? j_w : j_e;
            const real_t al = (r[3] * ja[0] + r[4] * ja[1] + r[5] * ja[2]) / (Vt * Vt);
            const real_t jn2 = r[3] * r[3] + r[4] * r[4] + r[5] * r[5];
            sum += al; cnt += 1.0;
            lo = std::min(lo, al); hi = std::max(hi, al);
            off = std::max(off, std::sqrt(std::max(0.0, jn2 / (Vt * Vt) - al * al)));
         }
         sum = GSum(mpi, sum); cnt = GSum(mpi, cnt);
         lo = -GMax(mpi, -lo); hi = GMax(mpi, hi); off = GMax(mpi, off);
         const real_t mean = sum / cnt;
         os.str("");
         os << "U-L1 A6: static solve, fault locked: the jump across " << cnt
            << " creep faces along the projected plate direction is " << mean
            << " Vp t on average (range " << lo << " to " << hi << "; largest part across it "
            << off << " Vp t); mean within " << jump_tol << " of 1";
         Check(mpi, cnt > 0 && std::abs(mean - 1.0) <= jump_tol, os.str());
      }
      }

      // ================= U-L2: option B =================
      if (mpi.IsRoot()) { std::cout << "\n=== U-L2 option B: walls only ===" << std::endl; }
      cycle::CycleConfig cfgB = cfg;
      cfgB.loading.creep_attrs.clear();
      cycle::QDPhase qdB(pmesh, mpi, cfgB);
      DomOp &domB = qdB.Domain();
      std::shared_ptr<const FaultFrameField> framesB = domB.FrameField();
      const long long censusB = framesB->Report().creep.faces;
      os.str("");
      os << "U-L2 B1: option B has " << censusB << " prescribed-slip faces (0)";
      Check(mpi, censusB == 0 && !framesB->HasCreep(), os.str());
      {
         DomainConfig dB;
         dB.penalty_factor = cfg.penalty_factor;
         dB.fault_frames = FaultFrameSource{framesB, nullptr};
         BoundaryConfig wbc;
         wbc.dirichlet_attrs = {5};
         wbc.natural_attrs = {1};
         wbc.fault_attr = cfg.fault_attr;
         wbc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);
         DomOp refB(pmesh, s.mesh.order, model, Vp, 0.0, 0.0, wbc, DGMethod::IP,
                    SolverType::MUMPS, dB);
         Vector bB, bW;
         domB.AssembleDirichletOnlyRHS(bB, t);
         refB.AssembleDirichletOnlyRHS(bW, t);
         const real_t e_w = RelDiff(mpi, bB, bW), dAB = RelDiff(mpi, bA, bB);
         os.str("");
         os << "U-L2 B2: option B's load equals the walls-only operator's with the BP5 function ("
            << e_w << ", 1e-12) and differs from option A's (" << dAB << " of its largest entry)";
         Check(mpi, GMax(mpi, bW.Normlinf()) > 0.0 && e_w <= 1e-12 && dAB > 1e-3, os.str());
      }
      if (solve)
      {
         ParGridFunction uB(&domB.GetFESpace());
         Vector zero(2 * domB.GetNumFaultDOFs());
         zero = 0.0;
         domB.Solve(t, zero, uB);
         const auto rows = FaceJumps(pmesh, uB, creep_attr,
                                     [&](int, const Vector &x, real_t *n)
         { const real_t *c = (x(0) < 0.0) ? n_w : n_e; for (int i = 0; i < 3; i++) { n[i] = c[i]; } });
         if (dump) { Dump(dump, "jumps_B", mpi.Rank(), rows); }
         real_t big = 0.0, cnt = 0.0;
         for (const auto &r : rows)
         {
            big = std::max(big, std::sqrt(r[3] * r[3] + r[4] * r[4] + r[5] * r[5]) / Vt);
            cnt += 1.0;
         }
         big = GMax(mpi, big); cnt = GSum(mpi, cnt);
         os.str("");
         os << "U-L2 B3: static solve, fault locked: the largest jump across " << cnt
            << " attribute-6 faces is " << big << " Vp t (interior faces of the DG solution: "
            << "below " << jump_tol << ")";
         Check(mpi, cnt > 0 && big <= jump_tol, os.str());
      }
   }
   MFEMFinalizePetsc();
   if (mpi.IsRoot()) { std::cout << (num_tests - num_failed) << "/" << num_tests << " tests passed" << std::endl; }
   return num_failed == 0 ? 0 : 1;
#endif
}
