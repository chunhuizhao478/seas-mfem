// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// Unit tests for BoundaryConfig, DirichletFunc, and MakeBP5DirichletFunc.
// Run: ./seas_test_boundary_config

#include "mfem.hpp"
#include "../../domain/boundary_config.hpp"
#include "test_macros.hpp"

#include <cmath>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

void TestBoundaryConfigDefaults()
{
   std::cout << "\n=== Test: BoundaryConfig defaults ===\n";
   BoundaryConfig bc;
   TEST_ASSERT(bc.dirichlet_attrs.empty(), "default dirichlet_attrs empty");
   TEST_ASSERT(bc.natural_attrs.empty(), "default natural_attrs empty");
   TEST_ASSERT(bc.fault_attr == 3, "default fault_attr == 3");
   TEST_ASSERT(bc.dirichlet_funcs.empty(), "default dirichlet_funcs empty");
   TEST_ASSERT(!bc.default_dirichlet_func, "default func is null");
}

void TestBoundaryConfigExplicitAttrs()
{
   std::cout << "\n=== Test: BoundaryConfig with explicit attrs ===\n";
   BoundaryConfig bc;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.fault_attr = 3;

   TEST_ASSERT(bc.dirichlet_attrs.count(5) == 1, "has Dirichlet attr 5");
   TEST_ASSERT(bc.natural_attrs.count(1) == 1, "has Natural attr 1");
   TEST_ASSERT(bc.dirichlet_attrs.count(1) == 0, "1 is not Dirichlet");
   TEST_ASSERT(bc.fault_attr == 3, "fault attr is 3");
}

void TestBoundaryConfigMultipleAttrs()
{
   std::cout << "\n=== Test: BoundaryConfig with multiple attrs ===\n";
   BoundaryConfig bc;
   bc.dirichlet_attrs = {1, 2, 3, 4};
   bc.natural_attrs = {5, 6};
   bc.fault_attr = 7;

   TEST_ASSERT(bc.dirichlet_attrs.size() == 4, "4 Dirichlet attrs");
   TEST_ASSERT(bc.natural_attrs.size() == 2, "2 Natural attrs");
   TEST_ASSERT(bc.dirichlet_attrs.count(3) == 1, "attr 3 is Dirichlet");
   TEST_ASSERT(bc.natural_attrs.count(6) == 1, "attr 6 is Natural");
}

void TestMakeBP5DirichletFunc()
{
   std::cout << "\n=== Test: MakeBP5DirichletFunc ===\n";
   real_t Vp = 1e-9;
   auto func = MakeBP5DirichletFunc(Vp);

   // Test y > 1000: u_D = (Vp*t/2, 0, 0)
   {
      Vector x(3); x(0) = 0; x(1) = 5000; x(2) = -10000;
      Vector u(3);
      real_t t = 1e9;  // 1 billion seconds
      func(x, t, u);
      real_t expected = Vp * t * 0.5;
      TEST_NEAR(u(0), expected, 1e-20, "y>1000: u_x = Vp*t/2");
      TEST_NEAR(u(1), 0.0, 1e-30, "y>1000: u_y = 0");
      TEST_NEAR(u(2), 0.0, 1e-30, "y>1000: u_z = 0");
   }

   // Test y < -1000: u_D = (-Vp*t/2, 0, 0)
   {
      Vector x(3); x(0) = 0; x(1) = -5000; x(2) = -10000;
      Vector u(3);
      real_t t = 1e9;
      func(x, t, u);
      real_t expected = -Vp * t * 0.5;
      TEST_NEAR(u(0), expected, 1e-20, "y<-1000: u_x = -Vp*t/2");
   }

   // Test |y| <= 1000: u_D = (Vp*t, 0, 0)
   {
      Vector x(3); x(0) = 0; x(1) = 0; x(2) = -10000;
      Vector u(3);
      real_t t = 1e9;
      func(x, t, u);
      real_t expected = Vp * t;
      TEST_NEAR(u(0), expected, 1e-20, "|y|<=1000: u_x = Vp*t");
   }

   // Test t = 0: u_D = (0, 0, 0)
   {
      Vector x(3); x(0) = 0; x(1) = 5000; x(2) = -10000;
      Vector u(3);
      func(x, 0.0, u);
      TEST_NEAR(u(0), 0.0, 1e-30, "t=0: u_x = 0");
   }
}

void TestPerAttrDirichletFunc()
{
   std::cout << "\n=== Test: Per-attr DirichletFunc dispatch ===\n";
   BoundaryConfig bc;
   bc.dirichlet_attrs = {5, 7};

   // Default function: constant (1, 0, 0)
   bc.default_dirichlet_func = [](const Vector &x, real_t t, Vector &u)
   {
      u.SetSize(3);
      u = 0.0;
      u(0) = 1.0;
   };

   // Per-attr function for attr 7: constant (0, 2, 0)
   bc.dirichlet_funcs[7] = [](const Vector &x, real_t t, Vector &u)
   {
      u.SetSize(3);
      u = 0.0;
      u(1) = 2.0;
   };

   Vector x(3); x = 0.0;
   Vector u(3);

   // Attr 5: should use default → (1, 0, 0)
   auto it5 = bc.dirichlet_funcs.find(5);
   if (it5 != bc.dirichlet_funcs.end())
   {
      it5->second(x, 0.0, u);
   }
   else if (bc.default_dirichlet_func)
   {
      bc.default_dirichlet_func(x, 0.0, u);
   }
   TEST_NEAR(u(0), 1.0, 1e-14, "attr 5 uses default: u_x = 1");
   TEST_NEAR(u(1), 0.0, 1e-14, "attr 5 uses default: u_y = 0");

   // Attr 7: should use per-attr → (0, 2, 0)
   auto it7 = bc.dirichlet_funcs.find(7);
   TEST_ASSERT(it7 != bc.dirichlet_funcs.end(), "attr 7 has per-attr func");
   it7->second(x, 0.0, u);
   TEST_NEAR(u(0), 0.0, 1e-14, "attr 7 per-attr: u_x = 0");
   TEST_NEAR(u(1), 2.0, 1e-14, "attr 7 per-attr: u_y = 2");
}

void TestDirichletFuncCapturesByValue()
{
   std::cout << "\n=== Test: DirichletFunc captures by value ===\n";
   real_t Vp = 1e-9;
   auto func = MakeBP5DirichletFunc(Vp);

   // Changing Vp after capture should NOT affect the lambda
   Vp = 999.0;

   Vector x(3); x(0) = 0; x(1) = 5000; x(2) = -10000;
   Vector u(3);
   func(x, 1e9, u);
   real_t expected = 1e-9 * 1e9 * 0.5;  // original Vp, not 999
   TEST_NEAR(u(0), expected, 1e-20, "lambda captures Vp by value");
}

int main(int argc, char *argv[])
{
   TestBoundaryConfigDefaults();
   TestBoundaryConfigExplicitAttrs();
   TestBoundaryConfigMultipleAttrs();
   TestMakeBP5DirichletFunc();
   TestPerAttrDirichletFunc();
   TestDirichletFuncCapturesByValue();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? 0 : 1;
}
