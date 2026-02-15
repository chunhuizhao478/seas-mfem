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

// Phase 7: Parallel MMS convergence test for antiplane operator
// Run with: mpirun -np 2 ./seas_mms_antiplane_parallel

#include "mfem.hpp"
#include "../../common/mpi_context.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../domain/bp2_mesh.hpp"

#include <iostream>
#include <cmath>

using namespace mfem;
using namespace mfem::seas;

/// MMS exact solution: u(x,z) = sin(a*(x+Lx)) * exp(a*z)
/// where a = pi/(2*Lx). Satisfies Laplace equation when Lx = Lz.
class MMSSolution : public Coefficient
{
public:
   real_t Lx, Lz;
   MMSSolution(real_t Lx_, real_t Lz_) : Lx(Lx_), Lz(Lz_) {}

   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
      Vector x(2);
      T.Transform(ip, x);
      real_t a = M_PI / (2.0 * Lx);
      real_t xi = x(0) + Lx;
      return std::sin(a * xi) * std::exp(a * x(1));
   }
};

int main(int argc, char *argv[])
{
   MPIContext ctx(&argc, &argv);

   if (ctx.IsRoot())
   {
      std::cout << "===============================================" << std::endl;
      std::cout << "   Parallel MMS Convergence Test" << std::endl;
      std::cout << "   (np=" << ctx.Size() << ")" << std::endl;
      std::cout << "===============================================" << std::endl;
      std::cout << "\n  Level  h         L2 error     rate" << std::endl;
   }

   real_t prev_error = 0.0;
   real_t prev_h = 0.0;
   bool all_passed = true;

   for (int level = 0; level < 4; level++)
   {
      // Create serial mesh, then partition
      auto serial_mesh = BP2MeshGenerator::CreateMMSMesh(level);
      ParMesh pmesh(ctx.GetComm(), *serial_mesh);

      real_t Lx = 1.0;
      real_t Lz = 1.0;
      real_t mu = 1.0;
      real_t Vp = 0.0;
      real_t Wf = 0.0;  // No fault for MMS

      AntiplaneDomainOperator<ParMesh> op(pmesh, 1, mu, Vp, Wf);

      MMSSolution u_exact(Lx, Lz);
      ParGridFunction u_h(&op.GetFESpace());
      op.SolveMMS(u_exact, u_h);

      // Compute L2 error (parallel-aware)
      real_t l2_error = u_h.ComputeL2Error(u_exact);

      real_t h = 1.0 / (4 * (1 << level));

      if (ctx.IsRoot())
      {
         if (level > 0 && prev_error > 1e-14)
         {
            real_t rate = std::log(prev_error / l2_error) / std::log(prev_h / h);
            std::cout << "    " << level << "     " << h
                      << "   " << l2_error
                      << "   " << rate << std::endl;

            if (rate < 1.5)
            {
               std::cerr << "  FAIL: convergence rate " << rate
                         << " < 1.5" << std::endl;
               all_passed = false;
            }
         }
         else
         {
            std::cout << "    " << level << "     " << h
                      << "   " << l2_error << "   ---" << std::endl;
         }
      }

      prev_error = l2_error;
      prev_h = h;
   }

   // ====== BR2 MMS Convergence Test ======
   if (ctx.IsRoot())
   {
      std::cout << "\n  --- BR2 Method ---" << std::endl;
      std::cout << "  Level  h         L2 error     rate" << std::endl;
   }

   real_t br2_prev_error = 0.0;
   real_t br2_prev_h = 0.0;

   for (int level = 0; level < 4; level++)
   {
      auto serial_mesh = BP2MeshGenerator::CreateMMSMesh(level);
      ParMesh pmesh(ctx.GetComm(), *serial_mesh);

      real_t Lx = 1.0;
      real_t Lz = 1.0;
      real_t mu = 1.0;
      real_t Vp = 0.0;
      real_t Wf = 0.0;

      AntiplaneDomainOperator<ParMesh> op(pmesh, 1, mu, Vp, Wf,
                                           DGMethod::BR2);

      MMSSolution u_exact(Lx, Lz);
      ParGridFunction u_h(&op.GetFESpace());
      op.SolveMMS(u_exact, u_h);

      real_t l2_error = u_h.ComputeL2Error(u_exact);
      real_t h = 1.0 / (4 * (1 << level));

      if (ctx.IsRoot())
      {
         if (level > 0 && br2_prev_error > 1e-14)
         {
            real_t rate = std::log(br2_prev_error / l2_error) /
                          std::log(br2_prev_h / h);
            std::cout << "    " << level << "     " << h
                      << "   " << l2_error
                      << "   " << rate << std::endl;

            if (rate < 1.5)
            {
               std::cerr << "  FAIL: BR2 convergence rate " << rate
                         << " < 1.5" << std::endl;
               all_passed = false;
            }
         }
         else
         {
            std::cout << "    " << level << "     " << h
                      << "   " << l2_error << "   ---" << std::endl;
         }
      }

      br2_prev_error = l2_error;
      br2_prev_h = h;
   }

   if (ctx.IsRoot())
   {
      std::cout << "\nParallel MMS test: "
                << (all_passed ? "PASSED" : "FAILED") << std::endl;
   }

   return all_passed ? 0 : 1;
}
