// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// REVIEW.md R-003 regression: verify that the base-class default
// `DomainOperator::ComputeTractionDiagnostics` sizes the diagnostic
// outputs (`traction_stress`, `traction_correction`, `jump_residual`)
// to match `traction.Size()` rather than leaving them empty.
//
// Pre-fix the default did `traction_stress.SetSize(0)`, which left
// downstream consumers like `FaceTraceLogger::RecordMult` reading
// past the end of an empty Vector when an Antiplane operator was
// combined with an active face tracer.

#include "mfem.hpp"
#include "../../domain/antiplane_operator.hpp"
#include "../../domain/bp2_mesh.hpp"
#include "../../config/bp2_params.hpp"

#include <iostream>
#include <memory>

using namespace mfem;
using namespace mfem::seas;

static int num_tests_passed = 0;
static int num_tests_failed = 0;

#define TEST_ASSERT(condition, message) do { \
   if (!(condition)) { \
      std::cerr << "FAILED: " << message << " (line " << __LINE__ << ")\n"; \
      num_tests_failed++; \
   } else { \
      std::cout << "  PASSED: " << message << "\n"; \
      num_tests_passed++; \
   } } while (0)

int main(int, char**)
{
   std::cout << "[R-003] Antiplane ComputeTractionDiagnostics default sizing\n";

   BP2Params params;
   BP2MeshGenerator::Parameters mp;
   mp.Lx = 10.0e3; mp.Lz = 10.0e3; mp.nx = 4; mp.nz = 4;
   auto mesh = BP2MeshGenerator::Create(mp);

   const int order = 1;
   AntiplaneDomainOperator<Mesh> domain(*mesh, order, params.mu(),
                                        params.Vp, mp.Lz);
   const int nf = domain.GetNumFaultDOFs();
   TEST_ASSERT(nf > 0, "BP2 mesh produces non-zero fault DOFs");

   // Solve once to get a non-trivial displacement.
   DG_FECollection fec(order, 2, BasisType::GaussLobatto);
   FiniteElementSpace fes(mesh.get(), &fec);
   GridFunction u(&fes); u = 0.0;
   Vector slip(nf); slip = 0.0;
   domain.Solve(0.0, slip, u);

   // Pre-size the outputs differently so we can confirm the default
   // resizes them (not just leaves them at the caller's old size).
   Vector traction(7);              traction = 7.0;
   Vector traction_stress(99);      traction_stress = 99.0;
   Vector traction_correction(99);  traction_correction = 99.0;
   Vector jump_residual(99);        jump_residual = 99.0;
   Vector normal_traction(99);      normal_traction = 99.0;
   Vector normal_stress(99);        normal_stress = 99.0;
   Vector normal_correction(99);    normal_correction = 99.0;

   // Call through a base-class reference so we exercise the virtual
   // default (Antiplane inherits it; ElasticityDomainOperator would
   // dispatch its own override).
   DomainOperator<Mesh> &base = domain;
   base.ComputeTractionDiagnostics(u, slip, traction,
                                   traction_stress,
                                   traction_correction,
                                   jump_residual,
                                   &normal_traction,
                                   &normal_stress,
                                   &normal_correction);

   // R-003 contract: every diagnostic vector must be sized to match
   // `traction` (so face-tracer downstream can iterate without OOB).
   TEST_ASSERT(traction.Size() == nf,
               "traction sized to NumFaultDOFs after default ComputeTractionDiagnostics");
   TEST_ASSERT(traction_stress.Size() == traction.Size(),
               "traction_stress sized to match traction (was 99 before)");
   TEST_ASSERT(traction_correction.Size() == traction.Size(),
               "traction_correction sized to match traction (was 99 before)");
   TEST_ASSERT(jump_residual.Size() == traction.Size(),
               "jump_residual sized to match traction (was 99 before)");

   // Diagnostics should be zero (the default fills with zeros).
   real_t max_abs_stress = 0.0, max_abs_corr = 0.0, max_abs_jr = 0.0;
   for (int i = 0; i < traction_stress.Size(); i++)
   {
      max_abs_stress = std::max(max_abs_stress, std::abs(traction_stress(i)));
      max_abs_corr   = std::max(max_abs_corr,   std::abs(traction_correction(i)));
      max_abs_jr     = std::max(max_abs_jr,     std::abs(jump_residual(i)));
   }
   TEST_ASSERT(max_abs_stress == 0.0,
               "traction_stress entries are all zero in default");
   TEST_ASSERT(max_abs_corr == 0.0,
               "traction_correction entries are all zero in default");
   TEST_ASSERT(max_abs_jr == 0.0,
               "jump_residual entries are all zero in default");

   // Optional normal_* outputs stay empty per the contract.
   TEST_ASSERT(normal_stress.Size() == 0,
               "normal_stress stays empty (caller guarded by elastic_sigma_n_)");
   TEST_ASSERT(normal_correction.Size() == 0,
               "normal_correction stays empty");

   std::cout << "========================================\n";
   std::cout << "R-003: " << num_tests_passed << " / "
             << (num_tests_passed + num_tests_failed)
             << " passed, " << num_tests_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_tests_failed == 0) ? 0 : 1;
}
