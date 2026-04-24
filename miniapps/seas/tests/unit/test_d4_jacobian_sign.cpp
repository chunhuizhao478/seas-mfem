// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// R3-R003 diagnostic: check whether BuildD4Mesh's orient=-1 tets have
// negative Jacobian determinants.  If they do, the DG volume kernel's
// integrand sees different-sign Jacobians on y-mirror-paired tets, and
// the per-DOF rhs picks up a sign flip that makes orbit-paired values
// equal-opposite instead of equal.  This is NOT a kernel bug — it is
// a fixture orientation issue masquerading as "antisymmetric" signature.

#include "mfem.hpp"
#include "../../dynamic/d4_tet_mesh.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

int main()
{
   constexpr real_t kL = 1000.0;
   std::cout << "=== D4 fixture Jacobian sign audit ===\n";
   Mesh mesh = BuildD4Mesh(false, kL);

   int n_pos = 0, n_neg = 0, n_zero = 0;
   int n_pos_lower = 0, n_neg_lower = 0, n_pos_upper = 0, n_neg_upper = 0;
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      ElementTransformation *Tr = mesh.GetElementTransformation(e);
      IntegrationPoint ip;  ip.Init(0);  // reference-tet center is fine
      ip.x = 0.25; ip.y = 0.25; ip.z = 0.25; ip.weight = 1.0;
      Tr->SetIntPoint(&ip);
      const real_t detJ = Tr->Jacobian().Det();
      Array<int> ev; mesh.GetElementVertices(e, ev);
      real_t cy = 0.0;
      for (int i = 0; i < ev.Size(); i++)
      {
         cy += mesh.GetVertex(ev[i])[1];
      }
      cy /= ev.Size();
      const bool is_upper = (cy > 0.5 * kL);
      if (detJ > 0)      { n_pos++; if (is_upper) n_pos_upper++; else n_pos_lower++; }
      else if (detJ < 0) { n_neg++; if (is_upper) n_neg_upper++; else n_neg_lower++; }
      else               { n_zero++; }
   }
   std::cout << "  total elements    : " << mesh.GetNE() << "\n";
   std::cout << "  positive det(J)   : " << n_pos
             << "   (" << n_pos_lower << " lower, " << n_pos_upper << " upper)\n";
   std::cout << "  negative det(J)   : " << n_neg
             << "   (" << n_neg_lower << " lower, " << n_neg_upper << " upper)\n";
   std::cout << "  zero det(J)       : " << n_zero << "\n";

   bool orientation_uniform = (n_neg == 0) || (n_pos == 0);
   std::cout << "\n  D4 fixture orientation uniform: "
             << (orientation_uniform ? "PASS" : "FAIL")
             << "   (mixed-sign Jacobians break y-mirror symmetry of\n"
                "    ∫ shape_i ∂Q/∂x dV even on slot-equivariant meshes)\n";
   return orientation_uniform ? 0 : 1;
}
