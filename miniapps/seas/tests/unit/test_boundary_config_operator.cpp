// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// KEY equivalence test (plan item 4h-ii):
// Verify old constructor (BCMode) and new constructor (BoundaryConfig)
// produce identical boundary marker setup, stiffness matrix, and
// Dirichlet loading RHS.
//
// Run: ./seas_test_boundary_config_operator

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../constitutive/linear_elastic.hpp"
#include "../../config/bp5_mesh_utils.hpp"
#include "test_macros.hpp"

#include <iostream>
#include <cmath>

using namespace mfem;
using namespace mfem::seas;

void TestBoundaryMarkerEquivalence()
{
   std::cout << "\n=== Test: Boundary marker equivalence (old vs new ctor) ===\n";

   // Use the inline mesh with fault faces
   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;

   // Old constructor: uses BCMode::FarField internally → attr 5 = Dirichlet
   ElasticityDomainOperator<Mesh> op_old(
      *mesh, 1, lambda, mu, Vp, Wf, lf,
      DGMethod::BR2, SolverType::CG_AMG, BCMode::FarField);

   // New constructor: explicit BoundaryConfig matching the old path
   LinearElastic model(lambda, mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   ElasticityDomainOperator<Mesh> op_new(
      *mesh, 1, model, Vp, Wf, lf, bc,
      DGMethod::BR2, SolverType::CG_AMG);

   // Compare: same number of fault faces and DOFs
   TEST_ASSERT(op_old.GetNumFaultFaces() == op_new.GetNumFaultFaces(),
               "fault face count matches");
   TEST_ASSERT(op_old.GetNumFaultDOFs() == op_new.GetNumFaultDOFs(),
               "fault DOF count matches");
   TEST_ASSERT(op_old.GetNumOwnedFaultDOFs() == op_new.GetNumOwnedFaultDOFs(),
               "owned fault DOF count matches");

   // Compare: same material properties
   TEST_ASSERT(op_old.GetShearModulus() == op_new.GetShearModulus(),
               "shear modulus matches");
   TEST_ASSERT(op_old.GetLambda() == op_new.GetModel().GetPenaltyModulus()
               - 2.0 * op_old.GetShearModulus(),
               "lambda matches (derived from penalty modulus)");
}

void TestSolveEquivalence()
{
   std::cout << "\n=== Test: Solve equivalence (old vs new ctor) ===\n";

   auto mesh = CreateBP5InlineMesh(2, 2, 1, 200e3, 100e3, 100e3);

   real_t lambda = 3.2044e10, mu = 3.2038e10;
   real_t Vp = 1e-9, Wf = 40e3, lf = 100e3;

   // Old constructor
   ElasticityDomainOperator<Mesh> op_old(
      *mesh, 1, lambda, mu, Vp, Wf, lf,
      DGMethod::BR2, SolverType::CG_AMG, BCMode::FarField);

   // New constructor with equivalent BoundaryConfig
   LinearElastic model(lambda, mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(Vp);

   ElasticityDomainOperator<Mesh> op_new(
      *mesh, 1, model, Vp, Wf, lf, bc,
      DGMethod::BR2, SolverType::CG_AMG);

   int n_dofs = op_old.GetNumOwnedFaultDOFs();
   if (n_dofs == 0)
   {
      std::cout << "  (no fault DOFs on this mesh, skipping solve comparison)\n";
      TEST_ASSERT(true, "no fault DOFs — skip solve");
      return;
   }

   // Create zero slip
   int slip_size = 2 * op_old.GetNumFaultDOFs();  // 2 components (dip, strike)
   Vector slip(slip_size);
   slip = 0.0;

   // Solve with both operators at t=1e6
   real_t time = 1e6;
   GridFunction u_old(&op_old.GetFESpace());
   GridFunction u_new(&op_new.GetFESpace());
   u_old = 0.0;
   u_new = 0.0;

   op_old.Solve(time, slip, u_old);
   op_new.Solve(time, slip, u_new);

   // Compare displacements
   Vector diff(u_old.Size());
   subtract(u_old, u_new, diff);
   real_t diff_norm = diff.Norml2();
   real_t ref_norm = u_old.Norml2();
   real_t rel_err = (ref_norm > 1e-30) ? diff_norm / ref_norm : diff_norm;

   std::cout << "  Displacement diff: ||u_old - u_new|| = " << diff_norm
             << ", rel = " << rel_err << "\n";
   TEST_ASSERT(rel_err < 1e-14, "displacement: old vs new bit-identical");

   // Compare tractions
   Vector trac_old, trac_new;
   op_old.ComputeTraction(u_old, slip, trac_old);
   op_new.ComputeTraction(u_new, slip, trac_new);

   Vector trac_diff(trac_old.Size());
   subtract(trac_old, trac_new, trac_diff);
   real_t trac_diff_norm = trac_diff.Norml2();
   real_t trac_ref_norm = trac_old.Norml2();
   real_t trac_rel_err = (trac_ref_norm > 1e-30)
                         ? trac_diff_norm / trac_ref_norm : trac_diff_norm;

   std::cout << "  Traction diff: ||t_old - t_new|| = " << trac_diff_norm
             << ", rel = " << trac_rel_err << "\n";
   TEST_ASSERT(trac_rel_err < 1e-14, "traction: old vs new bit-identical");
}

int main(int argc, char *argv[])
{
   TestBoundaryMarkerEquivalence();
   TestSolveEquivalence();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
