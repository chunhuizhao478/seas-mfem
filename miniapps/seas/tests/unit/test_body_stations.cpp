// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_body_stations.cpp — the off-fault station writer of the Nitsche
// dynamic-rupture driver (dynamic/body_station_writer.hpp; TPV34 plan rev 2,
// step S3; Chunhui Zhao, 2026-11-12).
//
//   T1  Barycentric coordinates of points in and outside a tetrahedron.
//   T2  Copies: one element for a point inside an element or on the outer
//       boundary, two on a shared face, more on a mesh vertex.
//   T3  Values, order 1 (linear field) and order 2 (quadratic field): every
//       row of every station file equals the analytic h, v, n displacement
//       and velocity (the mean over the copies of a continuous field is the
//       field), with rotated axes.
//   T4  Rows: 7 steps with a flush every 3 rows give 7 rows after Close.
//
// Run on 1, 2 and 3 ranks (make test-body-stations).

#include "mfem.hpp"
#include "../../dynamic/body_station_writer.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int g_pass = 0, g_fail = 0, g_order = 1;

void Check(bool ok, const std::string &what)
{
   if (Mpi::WorldRank() != 0) { return; }
   (ok ? g_pass : g_fail)++;
   std::cout << (ok ? "  PASSED: " : "  FAILED: ") << what << "\n";
}

/// The test displacement: linear for order 1, plus a quadratic part for 2.
void Disp(const Vector &x, Vector &u)
{
   const real_t X = x(0) / 1e3, Y = x(1) / 1e3, Z = x(2) / 1e3;
   u(0) = 0.3 + 0.2 * X - 0.1 * Y + 0.05 * Z;
   u(1) = -0.1 + 0.07 * X + 0.3 * Y - 0.2 * Z;
   u(2) = 0.4 - 0.15 * X + 0.02 * Y + 0.11 * Z;
   if (g_order >= 2)
   {
      u(0) += 0.03 * X * Y - 0.02 * Z * Z;
      u(1) += 0.01 * X * X + 0.04 * Y * Z;
      u(2) += -0.05 * X * Z + 0.02 * Y * Y;
   }
}

void Vel(const Vector &x, Vector &v)
{
   Vector u(3);
   Disp(x, u);
   v(0) = 2.0 * u(1) - 0.5;
   v(1) = -1.5 * u(2) + 0.25;
   v(2) = 0.7 * u(0) + 1.0;
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const int rank = Mpi::WorldRank(), nranks = Mpi::WorldSize();

   // T1
   {
      const real_t V[4][3] = {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}, {0, 0, 2}};
      const real_t Pin[3] = {0.5, 0.5, 0.5}, Pout[3] = {1.5, 1.5, 0.5};
      real_t l[4];
      bool ok = BodyStationWriter::Barycentric(V, Pin, l);
      Check(ok && std::fabs(l[0] - 0.25) < 1e-14 && std::fabs(l[1] - 0.25) < 1e-14
            && std::fabs(l[2] - 0.25) < 1e-14 && std::fabs(l[3] - 0.25) < 1e-14,
            "T1 barycentric coordinates inside");
      ok = BodyStationWriter::Barycentric(V, Pout, l);
      Check(ok && l[0] < -0.1, "T1 a point outside has a negative coordinate");
   }

   // stations: inside an element, on a mesh vertex, on the top boundary,
   // on a face shared by two elements (the plane x = 1000 between two cells)
   std::vector<FaultStation> st(4);
   st[0].name = "inside";   st[0].x[0] = 1234.5; st[0].x[1] = 2345.6; st[0].x[2] = 987.6;
   st[1].name = "vertex";   st[1].x[0] = 1000.0; st[1].x[1] = 1000.0; st[1].x[2] = 1000.0;
   st[2].name = "top";      st[2].x[0] = 1500.3; st[2].x[1] = 2100.7; st[2].x[2] = 4000.0;
   st[3].name = "face";     st[3].x[0] = 1000.0; st[3].x[1] = 1300.0; st[3].x[2] = 1750.0;
   // rotated axes (unit, orthogonal)
   const real_t c = std::cos(0.3), s = std::sin(0.3);
   const real_t axes[3][3] = {{c, s, 0.0}, {0.0, 0.0, -1.0}, {-s, c, 0.0}};

   char tmpl[] = "/tmp/seas_body_st_XXXXXX";
   std::string dir = "";
   if (rank == 0) { dir = mkdtemp(tmpl); }
   {
      int n = static_cast<int>(dir.size());
      MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
      dir.resize(n);
      MPI_Bcast(&dir[0], n, MPI_CHAR, 0, MPI_COMM_WORLD);
   }

   for (int order = 1; order <= 2; order++)
   {
      g_order = order;
      Mesh mesh = Mesh::MakeCartesian3D(4, 4, 4, Element::TETRAHEDRON, 4000.0, 4000.0, 4000.0);
      ParMesh pm(MPI_COMM_WORLD, mesh);
      L2_FECollection fec(order, 3);
      ParFiniteElementSpace fes(&pm, &fec, 3);
      ParGridFunction u0(&fes), v0(&fes);
      VectorFunctionCoefficient uc(3, Disp), vc(3, Vel);
      u0.ProjectCoefficient(uc);
      v0.ProjectCoefficient(vc);

      const std::string prefix = "p" + std::to_string(order);
      const int nsteps = 7;
      BodyStationWriter w(fes, st, axes, dir, prefix, 3);
      if (order == 1)
      {
         Check(w.NumCopies(0) == 1, "T2 inside an element: 1 copy");
         Check(w.NumCopies(1) > 2, "T2 on a mesh vertex: " + std::to_string(w.NumCopies(1))
               + " copies");
         Check(w.NumCopies(2) == 1, "T2 on the outer boundary: 1 copy");
         Check(w.NumCopies(3) == 2, "T2 on a shared face: 2 copies (got "
               + std::to_string(w.NumCopies(3)) + ")");
         int loc = w.LocalCopies(), tot = 0;
         MPI_Allreduce(&loc, &tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         Check(tot == w.NumCopies(0) + w.NumCopies(1) + w.NumCopies(2) + w.NumCopies(3),
               "T2 the ranks' copies add up");
      }
      Vector U(u0.Size()), Ud(v0.Size());
      for (int k = 0; k < nsteps; k++)
      {
         const real_t t = 0.1 * k;
         for (int i = 0; i < U.Size(); i++) { U(i) = (1.0 + t) * u0(i); Ud(i) = (1.0 - 2.0 * t) * v0(i); }
         w.WriteStep(t, U, Ud);
      }
      w.Close();

      if (rank == 0)
      {
         real_t worst = 0.0;
         bool rows_ok = true;
         for (int q = 0; q < 4; q++)
         {
            std::ifstream f(dir + "/" + prefix + "_offfault_" + st[q].name + ".dat");
            std::string line;
            int k = 0;
            Vector x(3), u(3), v(3);
            for (int d = 0; d < 3; d++) { x(d) = st[q].x[d]; }
            Disp(x, u);
            Vel(x, v);
            while (std::getline(f, line))
            {
               if (line.empty() || line[0] == '#') { continue; }
               std::istringstream is(line);
               real_t t, val[6];
               is >> t;
               for (int j = 0; j < 6; j++) { is >> val[j]; }
               if (std::fabs(t - 0.1 * k) > 1e-9) { rows_ok = false; }
               for (int a = 0; a < 3; a++)
               {
                  const real_t eu = (1.0 + t) * (u(0) * axes[a][0] + u(1) * axes[a][1] + u(2) * axes[a][2]);
                  const real_t ev = (1.0 - 2.0 * t) * (v(0) * axes[a][0] + v(1) * axes[a][1] + v(2) * axes[a][2]);
                  worst = std::max(worst, std::fabs(val[2 * a] - eu));
                  worst = std::max(worst, std::fabs(val[2 * a + 1] - ev));
               }
               k++;
            }
            if (k != nsteps) { rows_ok = false; }
            std::remove((dir + "/" + prefix + "_offfault_" + st[q].name + ".dat").c_str());
         }
         char buf[64];
         std::snprintf(buf, sizeof(buf), "%.3e", worst);
         Check(worst < 1e-9, "T3 order " + std::to_string(order) + ": station values = the field"
               " (max abs err " + buf + ")");
         Check(rows_ok, "T4 order " + std::to_string(order) + ": 7 rows with their times");
      }
   }
   if (rank == 0)
   {
      rmdir(dir.c_str());
      std::cout << "test_body_stations (" << nranks << " rank(s)): " << g_pass << " passed, "
                << g_fail << " failed\n";
   }
   return g_fail == 0 ? 0 : 1;
}
