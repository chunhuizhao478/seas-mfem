// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// initial_damage.hpp — the seeded lower bound of the damage (MOOSE
// initial_damage_aux, a ConstantAux on SubdomainBoundingBoxGenerator blocks;
// farms@cdms pulsepower/cmame_revision/*/fracture*.i).

#ifndef MFEM_PFF_INITIAL_DAMAGE_HPP
#define MFEM_PFF_INITIAL_DAMAGE_HPP

#include "mfem.hpp"
#include "../config/pff_config.hpp"

#include <algorithm>
#include <vector>

namespace mfem
{
namespace pff
{

/// Node-interior seeding of initial_damage (MOOSE initial_damage_aux): a
/// node gets a box value only if every element containing it has its vertex
/// average in the box (SubdomainBoundingBoxGenerator block); MOOSE writes the
/// value on the box block and then 0 on the other blocks, so nodes on the box
/// boundary end at 0 (PLAN.md, parity details). Several boxes: the maximum.
inline void SeedInitialDamage(ParFiniteElementSpace &sfes, const std::vector<DamageBox> &boxes,
                       Vector &seed)
{
   ParMesh &mesh = *sfes.GetParMesh();
   const int dim = mesh.Dimension();
   seed.SetSize(sfes.GetTrueVSize());
   seed = 0.0;
   const HypreParMatrix *P = sfes.Dof_TrueDof_Matrix();
   for (const DamageBox &b : boxes)
   {
      Vector in(sfes.GetVSize()), out(sfes.GetVSize());
      in = 0.0; out = 0.0;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         Array<int> verts;
         mesh.GetElementVertices(e, verts);
         real_t xc[3] = {0.0, 0.0, 0.0};
         for (int v : verts)
         {
            const real_t *x = mesh.GetVertex(v);
            for (int k = 0; k < dim; k++) { xc[k] += x[k] / verts.Size(); }
         }
         bool inside = true;
         for (int k = 0; k < dim; k++) { inside = inside && xc[k] >= b.lo[k] && xc[k] <= b.hi[k]; }
         Array<int> dofs;
         sfes.GetElementDofs(e, dofs);
         for (int ld : dofs) { (inside ? in : out)(ld) += 1.0; }
      }
      Vector tin(sfes.GetTrueVSize()), tout(sfes.GetTrueVSize());
      P->MultTranspose(in, tin);
      P->MultTranspose(out, tout);
      for (int i = 0; i < seed.Size(); i++)
      {
         if (tin(i) > 0.0 && tout(i) == 0.0) { seed(i) = std::max(seed(i), b.value); }
      }
   }
}

} // namespace pff
} // namespace mfem

#endif
