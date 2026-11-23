// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_cycle_hetero.cpp — unit tests of the heterogeneity inputs of the
// cycle (heterogeneity plan rev4; branch fully-dynamic-fault-opening-hetergeneity).
//
//   make -f cycle/cycle.mk seas_test_cycle_hetero SEAS_FS_LIB= && mpirun -np 1 ./seas_test_cycle_hetero
//
// T1 parser: [velocity] kind = "grid" and [stress] kind = "grid_tensor" with
//    their file / sampling / outside keys parse; every malformed form is refused (forked child).
// T2 GridFieldCoefficient: the value at a quadrature point is the grid at the
//    point's physical coordinates (bit for bit), and the grid is exact for a
//    linear field.
// T3 HeterogeneityConfig: CycleConfig::Load fills it from the deck; Material()
//    and StressAt() read the file; the restart tag names the file hash.
// T4 per-element sampling (SampleElementMaterial): StiffestPoint and CellAverage
//    against an independent hand computation; a uniform material returns its
//    values bit for bit in every mode.
// T5 rho-weighted NitscheMass: constant rho gives the old mass bit for bit;
//    rho(x) linear gives MFEM's MassIntegrator with the same rule to 1e-13.
// T6 volume stiffness with GridFieldCoefficient lambda(x), mu(x): for the
//    constant-strain field u = (x, 0), u^T K_vol u = int (lambda + 2 mu) dx,
//    against the exact integral of the linear fields (1e-12).

#include "mfem.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/grid_material.hpp"
#include "../../domain/nitsche_material_sample.hpp"
#include "../../domain/nitsche_mass.hpp"
#include "../../solver/nitsche_lts_run.hpp"

#include <hdf5.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
#define CHECK(cond, msg) do { num_tests++; if (!(cond)) { num_failed++; \
   std::cerr << "FAILED: " << msg << " (line " << __LINE__ << ")\n"; } \
   else { std::cout << "  PASSED: " << msg << "\n"; } } while (0)

// grid: lambda = L0 + Lx x + Ly y, mu = M0 + Mx x + My y, rho = R0 + Rx x
static const double L0 = 3.0e10, Lx = 2.0e4, Ly = -1.0e4;
static const double M0 = 3.2e10, Mx = 1.5e4, My = 2.5e4;
static const double R0 = 2700.0, Rx = 1.0e-3;
static double Lam(double x, double y) { return L0 + Lx * x + Ly * y; }
static double Mu(double x, double y) { return M0 + Mx * x + My * y; }
static double Rho(double x, double) { return R0 + Rx * x; }
// tensor (Pa, tension +): s_xx = -60e6 + 10 x, s_yy = -40e6 - 5 y, s_xy = 3e6 + 2 x, s_zz = -50e6
static double Sxx(double x, double) { return -60e6 + 10.0 * x; }
static double Syy(double, double y) { return -40e6 - 5.0 * y; }
static double Sxy(double x, double) { return 3e6 + 2.0 * x; }

static void WriteGrid(const std::string &path, bool uniform)
{
   std::vector<double> x, y;
   for (int i = 0; i <= 40; i++) { x.push_back(-2000.0 + 100.0 * i); }
   for (int j = 0; j <= 30; j++) { y.push_back(-1500.0 + 100.0 * j); }
   const size_t nx = x.size(), ny = y.size();
   std::vector<std::pair<std::string, std::vector<double>>> F = {
      {"lambda", {}}, {"mu", {}}, {"rho", {}}, {"s_xx", {}}, {"s_yy", {}}, {"s_xy", {}}, {"s_zz", {}}};
   for (auto &f : F) { f.second.resize(nx * ny); }
   for (size_t j = 0; j < ny; j++)
      for (size_t i = 0; i < nx; i++)
      {
         const double X = uniform ? 0.0 : x[i], Y = uniform ? 0.0 : y[j];
         const size_t k = j * nx + i;
         F[0].second[k] = Lam(X, Y); F[1].second[k] = Mu(X, Y); F[2].second[k] = Rho(X, Y);
         F[3].second[k] = Sxx(X, Y); F[4].second[k] = Syy(X, Y); F[5].second[k] = Sxy(X, Y);
         F[6].second[k] = -50e6;
      }
   hid_t f = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   auto w1 = [&](const char *name, const std::vector<double> &v)
   {
      hsize_t n = v.size(); hid_t sp = H5Screate_simple(1, &n, nullptr);
      hid_t d = H5Dcreate2(f, name, H5T_IEEE_F64LE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()); H5Dclose(d); H5Sclose(sp);
   };
   w1("x", x); w1("y", y);
   for (const auto &fv : F)
   {
      hsize_t dims[2] = {ny, nx}; hid_t sp = H5Screate_simple(2, dims, nullptr);
      hid_t d = H5Dcreate2(f, fv.first.c_str(), H5T_IEEE_F64LE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fv.second.data());
      H5Dclose(d); H5Sclose(sp);
   }
   const std::string schema = "seas_grid2d_v1";
   hid_t st = H5Tcopy(H5T_C_S1); H5Tset_size(st, schema.size());
   hid_t asp = H5Screate(H5S_SCALAR);
   hid_t a = H5Acreate2(f, "schema", st, asp, H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(a, st, schema.c_str()); H5Aclose(a); H5Sclose(asp); H5Tclose(st);
   H5Fclose(f);
}

static bool AbortsInChild(const std::function<void()> &fn)
{
   std::cout.flush(); std::cerr.flush();
   pid_t pid = fork();
   if (pid == 0) { freopen("/dev/null", "w", stderr); freopen("/dev/null", "w", stdout); fn(); _exit(0); }
   int status = 0; waitpid(pid, &status, 0);
   return !(WIFEXITED(status) && WEXITSTATUS(status) == 0);
}

static std::string ReadFile(const std::string &p)
{
   std::ifstream f(p); std::stringstream s; s << f.rdbuf(); return s.str();
}

/// The mapview tensor deck with its [velocity] / [stress] sections replaced.
static std::string MakeDeck(const std::string &base, const std::string &velocity, const std::string &stress,
                            const std::string &extra)
{
   std::string d = ReadFile(base);
   auto replace_section = [&d](const std::string &head, const std::string &body)
   {
      const size_t i0 = d.find("\n" + head + "\n");
      MFEM_VERIFY(i0 != std::string::npos, "test deck: no " << head);
      size_t i1 = d.find("\n[", i0 + 1);
      const size_t ic = d.find("\n#", i0 + 1);
      if (ic != std::string::npos && ic < i1) { i1 = ic; }
      d = d.substr(0, i0) + "\n" + head + "\n" + body + d.substr(i1);
   };
   replace_section("[velocity]", velocity);
   replace_section("[stress]", stress);
   return d + "\n" + extra;
}

int main(int argc, char **argv)
{
   const std::string grid = "/tmp/seas_test_cycle_hetero_grid.h5";
   const std::string grid_u = "/tmp/seas_test_cycle_hetero_uniform.h5";
   WriteGrid(grid, false);
   WriteGrid(grid_u, true);
   const std::string base = "config/fd_cycle/mapview_rs_tensor_fd.toml";
   const std::string src = "";
   const std::string vel = "kind = \"grid\"\nfile = \"" + grid + "\"\nsampling = \"quadrature\"\n";
   const std::string str = "kind = \"grid_tensor\"\nfile = \"" + grid + "\"\n";
   const std::string gfile = "file = \"" + grid + "\"\n";
   auto write_deck = [](const std::string &text) { const std::string p = "/tmp/seas_test_cycle_hetero_deck.toml";
                                                   std::ofstream(p) << text; return p; };

   std::cout << "T1 parser\n";
   {
      const spatial::SpatialFrictionConfig s =
         spatial::LoadSpatialFrictionConfig(write_deck(MakeDeck(base, vel, str, src)));
      CHECK(s.velocity.kind == "grid" && s.velocity.file == grid && s.velocity.sampling == "quadrature"
            && s.velocity.outside == "abort" && !s.velocity.use_sidecar,
            "[velocity] kind = grid parsed (use_sidecar forced off, outside defaults to abort)");
      CHECK(s.stress.kind == spatial::StressSourceKind::GridTensor && s.stress.grid_file == grid
            && s.stress.grid_outside == "abort", "[stress] kind = grid_tensor parsed");
      const spatial::SpatialFrictionConfig s0 = spatial::LoadSpatialFrictionConfig(base);
      CHECK(s0.velocity.kind == "constant" && s0.stress.kind == spatial::StressSourceKind::ConstantTensor
            && s0.velocity.file.empty(), "an existing deck parses as before");
      auto refused = [&](const std::string &v, const std::string &st, const std::string &ex)
      { const std::string p = write_deck(MakeDeck(base, v, st, ex));
        return AbortsInChild([p] { cycle::CycleConfig::Load(p); }); };
      CHECK(refused("kind = \"grid\"\n", str, src), "grid without file refused");
      CHECK(refused("kind = \"grid\"\n" + gfile + "use_sidecar = true\n", str, src), "grid with use_sidecar = true refused");
      CHECK(refused("kind = \"grid\"\n" + gfile + "sampling = \"mean\"\n", str, src), "unknown sampling refused");
      CHECK(refused("use_sidecar = false\nsampling = \"quadrature\"\n", str, src), "sampling without kind = grid refused");
      CHECK(refused("kind = \"voxels\"\n", str, src), "unknown velocity kind refused");
      CHECK(refused(vel, "kind = \"grid_tensor\"\n" + gfile + "sigma_xx_pa = 1.0\n", src), "grid_tensor with sigma keys refused");
      CHECK(refused(vel, "kind = \"grid_tensor\"\n", src), "grid_tensor without file refused");
      CHECK(refused("kind = \"grid\"\nfile = \"/tmp/no_such_grid.nc\"\n", str, src), "missing file refused");
      CHECK(refused(vel + "outside = \"wrap\"\n", str, src), "unknown outside policy refused");
   }

   // MPI after T1: the refusal checks fork, and a forked MPI_Abort could reach the parent
   Mpi::Init(argc, argv);
   std::cout << "T2/T3 coefficient and HeterogeneityConfig\n";
   const cycle::CycleConfig cfg = cycle::CycleConfig::Load(write_deck(MakeDeck(base, vel, str, src)));
   CHECK(cfg.hetero.velocity && cfg.hetero.stress && cfg.hetero.Pointwise() && cfg.TensorPrestress(),
         "CycleConfig::Load: grid velocity and tensor on");
   CHECK(cfg.IcLabel().find("+gridvel:quadrature:" + cfg.hetero.vel_hash) != std::string::npos
         && cfg.IcLabel().find("+gridtensor:" + cfg.hetero.str_hash) != std::string::npos,
         "restart tag names the sampling and the file hash");
   Mesh smesh = Mesh::MakeCartesian2D(8, 6, Element::TRIANGLE, true, 1600.0, 1200.0);
   smesh.Transform([](const Vector &X, Vector &Y) { Y(0) = X(0) - 900.0; Y(1) = X(1) - 700.0; });
   ParMesh pmesh(MPI_COMM_WORLD, smesh);
   const MaterialField mq = cfg.hetero.Material();
   {
      bool same = true; double err = 0.0;
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
         ElementTransformation *T = pmesh.GetElementTransformation(e);
         const IntegrationRule &ir = IntRules.Get(pmesh.GetElementBaseGeometry(e), 4);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            T->SetIntPoint(&ip);
            Vector X(2); T->Transform(ip, X);
            real_t l, m, r;
            mq.EvalAt(e, *T, ip, l, m, r);
            real_t gl, gm, gr;
            const real_t x3[3] = {X(0), X(1), 0.0};
            cfg.hetero.MaterialAt(x3, gl, gm, gr);
            same &= (l == gl && m == gm && r == gr);
            err = std::max({err, std::abs(l / Lam(X(0), X(1)) - 1), std::abs(m / Mu(X(0), X(1)) - 1),
                            std::abs(r / Rho(X(0), X(1)) - 1)});
            real_t s[9];
            cfg.hetero.StressAt(x3, s);
            err = std::max({err, std::abs(s[0] / Sxx(X(0), X(1)) - 1), std::abs(s[4] / Syy(X(0), X(1)) - 1),
                            std::abs(s[1] / Sxy(X(0), X(1)) - 1)});
            same &= (s[1] == s[3] && s[2] == 0.0 && s[5] == 0.0 && s[8] == -50e6);
         }
      }
      std::cout << "    max relative error of linear fields: " << err << "\n";
      CHECK(same, "coefficient == grid at the physical point; tensor symmetric, s_zz from the file");
      CHECK(err < 1e-12, "linear lambda, mu, rho and tensor reproduced to 1e-12");
   }

   std::cout << "T4 per-element sampling\n";
   {
      const int order = 1;
      bool stiff_ok = true, avg_ok = true, cen_ok = true;
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
         ElementTransformation *T = pmesh.GetElementTransformation(e);
         const Geometry::Type g = pmesh.GetElementBaseGeometry(e);
         const IntegrationRule &ir = IntRules.Get(g, 2 * order);
         // hand: stiffest point and SeisSol average
         real_t best = -1, bl = 0, bm = 0, br = 0, ws = 0, rm = 0, mi = 0, ve = 0;
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            Vector X(2); T->Transform(ir.IntPoint(q), X);
            const double l = Lam(X(0), X(1)), m = Mu(X(0), X(1)), r = Rho(X(0), X(1));
            if ((l + 2 * m) / r > best) { best = (l + 2 * m) / r; bl = l; bm = m; br = r; }
            const double w = ir.IntPoint(q).weight;
            ws += w; rm += w * r; mi += w / m; ve += w * l / (2 * m * (3 * l + 2 * m));
         }
         rm /= ws; mi /= ws; ve /= ws;
         const double am = 1.0 / mi, al = 4 * am * am * ve / (1 - 6 * am * ve);
         MaterialField ms = mq; ms.sampling = MaterialField::Sampling::StiffestPoint;
         real_t l, m, r;
         SampleElementMaterial(ms, *T, e, g, order, l, m, r);
         stiff_ok &= std::abs(l / bl - 1) < 1e-12 && std::abs(m / bm - 1) < 1e-12 && std::abs(r / br - 1) < 1e-12;
         ms.sampling = MaterialField::Sampling::CellAverage;
         SampleElementMaterial(ms, *T, e, g, order, l, m, r);
         avg_ok &= std::abs(l / al - 1) < 1e-11 && std::abs(m / am - 1) < 1e-12 && std::abs(r / rm - 1) < 1e-12;
         ms.sampling = MaterialField::Sampling::Centroid;
         SampleElementMaterial(ms, *T, e, g, order, l, m, r);
         Vector C(2); T->Transform(Geometries.GetCenter(g), C);
         cen_ok &= std::abs(m / Mu(C(0), C(1)) - 1) < 1e-12;
      }
      CHECK(stiff_ok, "StiffestPoint == the hand search over IntRules(geom, 2p)");
      CHECK(avg_ok, "CellAverage == SeisSol's rule computed by hand");
      CHECK(cen_ok, "Centroid == the value at the centroid");
      // uniform: every mode returns the stored values bit for bit
      cycle::HeterogeneityConfig hu;
      hu.velocity = true; hu.vel_path = grid_u; hu.Load();
      bool exact = true;
      for (auto smode : {MaterialField::Sampling::Centroid, MaterialField::Sampling::StiffestPoint,
                         MaterialField::Sampling::CellAverage})
      {
         MaterialField mu_ = hu.Material(); mu_.sampling = smode;
         for (int e = 0; e < pmesh.GetNE(); e++)
         {
            real_t l, m, r;
            SampleElementMaterial(mu_, *pmesh.GetElementTransformation(e), e, pmesh.GetElementBaseGeometry(e), 1, l, m, r);
            exact &= (l == L0 && m == M0 && r == R0);
         }
      }
      CHECK(exact, "uniform material: every sampling mode returns the stored values bit for bit");
   }

   std::cout << "T5 rho-weighted mass\n";
   {
      const int order = 1;
      NitscheMaterialSample sample(pmesh, mq, order);
      Vector rho_e; RhoPerElement(sample, rho_e);
      NitscheMass M0m(pmesh, order, rho_e);
      // constant rho per element as a pointwise coefficient: bit for bit the old mass
      class ElemConst : public Coefficient
      {
      public:
         explicit ElemConst(const Vector &v) : v_(v) {}
         real_t Eval(ElementTransformation &T, const IntegrationPoint &) override { return v_(T.ElementNo); }
      private:
         const Vector &v_;
      } rc(rho_e);
      NitscheMass M1(pmesh, order, rho_e, nullptr, 0.0, nullptr, &rc);
      bool same = true;
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
         const DenseMatrix &A = M0m.ElementMassInverse(e), &B = M1.ElementMassInverse(e);
         for (int i = 0; i < A.Height(); i++)
            for (int j = 0; j < A.Width(); j++) { same &= A(i, j) == B(i, j); }
      }
      CHECK(same, "constant rho per element: the pointwise mass equals the old one bit for bit");
      // rho(x) linear: (1/rho_e) Minv_e == (int rho phi phi)^-1 by MFEM's MassIntegrator, same rule
      NitscheMass M2(pmesh, order, rho_e, nullptr, 0.0, nullptr, mq.rho_coef);
      L2_FECollection fec(order, 2, BasisType::GaussLobatto);
      double err = 0.0;
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
         const FiniteElement &fe = *fec.FiniteElementForGeometry(pmesh.GetElementBaseGeometry(e));
         MassIntegrator mi(*mq.rho_coef);
         mi.SetIntRule(&IntRules.Get(pmesh.GetElementBaseGeometry(e), 2 * order));
         DenseMatrix Me; mi.AssembleElementMatrix(fe, *pmesh.GetElementTransformation(e), Me);
         DenseMatrixInverse inv(Me); DenseMatrix Mi; inv.GetInverseMatrix(Mi);
         const DenseMatrix &B = M2.ElementMassInverse(e);
         for (int i = 0; i < Mi.Height(); i++)
            for (int j = 0; j < Mi.Width(); j++)
            {
               err = std::max(err, std::abs(B(i, j) / rho_e(e) - Mi(i, j)) / Mi.MaxMaxNorm());
            }
      }
      std::cout << "    max relative difference: " << err << "\n";
      CHECK(err < 1e-13, "rho(x) linear: the pointwise mass equals MassIntegrator(rho) with the same rule");
   }

   std::cout << "T6 volume stiffness with lambda(x), mu(x)\n";
   {
      H1_FECollection fec(1, 2);
      ParFiniteElementSpace fes(&pmesh, &fec, 2);
      ParBilinearForm a(&fes);
      a.AddDomainIntegrator(new ElasticityIntegrator(*mq.lambda_coef, *mq.mu_coef));
      a.Assemble(); a.Finalize();
      ParGridFunction u(&fes);
      VectorFunctionCoefficient uc(2, [](const Vector &X, Vector &U) { U(0) = X(0); U(1) = 0.0; });
      u.ProjectCoefficient(uc);
      ParGridFunction Ku(&fes); a.Mult(u, Ku);
      const double e_h = InnerProduct(MPI_COMM_WORLD, u, Ku);
      // exact: int (lambda + 2 mu) over [-900, 700] x [-700, 500] for linear fields = (L + 2M) at the centre x area
      const double xc = -100.0, yc = -100.0, area = 1600.0 * 1200.0;
      const double e_x = (Lam(xc, yc) + 2.0 * Mu(xc, yc)) * area;
      std::cout << "    u^T K u = " << e_h << ", exact " << e_x << ", rel " << std::abs(e_h / e_x - 1) << "\n";
      CHECK(std::abs(e_h / e_x - 1) < 1e-12, "u^T K_vol u = int (lambda + 2 mu) dx for u = (x, 0)");
   }

   std::cout << "T7 local time stepping blocks: the volume moduli per quadrature point\n";
   {
      // NitscheLtsBlocks::Build with the operator's volume moduli (H4: the
      // LTS run took the per-element moduli for K, so K_LTS != K_QD and the
      // hand-back jumped).  The face terms take lam_e, mu_e in both builds,
      // so y(vol) - y(no vol) = (K_vol(lambda(x), mu(x)) - K_vol(lam_e, mu_e)) u.
      L2_FECollection fec(1, 2);
      ParFiniteElementSpace fes(&pmesh, &fec, 2);
      pmesh.ExchangeFaceNbrData();
      ConstantCoefficient lam_e(3.0e10), mu_e(3.2e10);
      Array<int> dir(pmesh.bdr_attributes.Max()); dir = 0;
      std::vector<int> cl(pmesh.GetNE(), 0), cl_nbr(pmesh.GetNFaceNeighborElements(), 0);
      ParGridFunction u(&fes);
      VectorFunctionCoefficient uc(2, [](const Vector &X, Vector &U)
      { U(0) = 1e-4 * X(0) + 3e-5 * X(1) + 1e-8 * X(0) * X(1); U(1) = -2e-5 * X(0) + 5e-9 * X(1) * X(1); });
      u.ProjectCoefficient(uc);
      auto apply = [&](Coefficient *vl, Coefficient *vm)
      {
         NitscheLtsBlocks b;
         b.Build(fes, lam_e, mu_e, nullptr, nullptr, dir, 10.0, cl, cl_nbr, 1, nullptr, nullptr, vl, vm);
         Vector y(fes.GetVSize()); y = 0.0;
         b.MultInterior(0, u, nullptr, y);
         b.AddShared(0, y);
         return y;
      };
      Vector y_pw = apply(mq.lambda_coef, mq.mu_coef), y_e = apply(nullptr, nullptr);
      auto kvol = [&](Coefficient &l, Coefficient &m)
      {
         ParBilinearForm a(&fes);
         a.AddDomainIntegrator(new ElasticityIntegrator(l, m));
         a.Assemble(); a.Finalize();
         Vector y(fes.GetVSize()); a.Mult(u, y); return y;
      };
      Vector d_blk(y_pw); d_blk -= y_e;
      Vector d_ref = kvol(*mq.lambda_coef, *mq.mu_coef); d_ref -= kvol(lam_e, mu_e);
      Vector r(d_blk); r -= d_ref;
      // round-off scale: |K u| (the blocks sum in their own order)
      const double rel = r.Normlinf() / y_pw.Normlinf();
      std::cout << "    |K u|_inf " << y_pw.Normlinf() << ", |dK u|_inf " << d_ref.Normlinf()
                << ", blocks - assembled " << r.Normlinf() << " (rel to |K u| " << rel << ")\n";
      CHECK(d_ref.Normlinf() > 1e-4 * y_pw.Normlinf() && rel < 1e-12,
            "LTS blocks with the volume moduli: K_vol per quadrature point, as AssembleStiffness");
      Vector y_e2 = apply(nullptr, nullptr); y_e2 -= y_e;
      CHECK(y_e2.Normlinf() == 0.0, "LTS blocks without volume moduli: unchanged (lam_e, mu_e)");
   }

   std::remove(grid.c_str()); std::remove(grid_u.c_str());
   std::cout << (num_failed ? "FAIL" : "PASS") << ": " << num_tests - num_failed << " / " << num_tests << "\n";
   return num_failed ? 1 : 0;
}
