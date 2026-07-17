// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// stress_field_3d.cpp — implementation of the six-component Cauchy
// stress reader (Phase 6 §1 of PLAN_onfaultstress.md).

#include "stress_field_3d.hpp"

#include "mfem.hpp"

#include <algorithm>
#include <array>
#include <sstream>
#include <string>

namespace mfem
{
namespace seas
{

namespace
{

// Schema-v1 canonical component names, indexed 0..5.
const std::array<std::string, 6>& component_names_()
{
   static const std::array<std::string, 6> names = {
      "sigma_xx", "sigma_yy", "sigma_zz",
      "sigma_xy", "sigma_yz", "sigma_xz"
   };
   return names;
}

} // namespace


const std::string& StressField3D::ComponentName(int component_index)
{
   MFEM_VERIFY(component_index >= 0 && component_index < 6,
               "StressField3D::ComponentName: component_index "
               << component_index << " out of range [0, 6)");
   return component_names_()[component_index];
}


StressField3D::StressField3D(const std::string& sidecar_path,
                             OOBPolicy oob)
   : sigma_xx_(sidecar_path, "sigma_xx", oob)
   , sigma_yy_(sidecar_path, "sigma_yy", oob)
   , sigma_zz_(sidecar_path, "sigma_zz", oob)
   , sigma_xy_(sidecar_path, "sigma_xy", oob)
   , sigma_yz_(sidecar_path, "sigma_yz", oob)
   , sigma_xz_(sidecar_path, "sigma_xz", oob)
{
   init_after_load_();
}

#ifdef MFEM_USE_MPI
StressField3D::StressField3D(const std::string& sidecar_path,
                             OOBPolicy oob,
                             MPI_Comm node_comm)
   : sigma_xx_(sidecar_path, "sigma_xx", oob, node_comm)
   , sigma_yy_(sidecar_path, "sigma_yy", oob, node_comm)
   , sigma_zz_(sidecar_path, "sigma_zz", oob, node_comm)
   , sigma_xy_(sidecar_path, "sigma_xy", oob, node_comm)
   , sigma_yz_(sidecar_path, "sigma_yz", oob, node_comm)
   , sigma_xz_(sidecar_path, "sigma_xz", oob, node_comm)
{
   init_after_load_();
}
#endif

void StressField3D::init_after_load_()
{
   // Pin bbox_ as the strict intersection of all six component
   // bboxes per plan §1 line 1644-1645 (R-903).  Under the
   // schema-v1 contract the six are bit-exact identical (single
   // /grid/* per file), so this reduces to the shared value;
   // computing the intersection explicitly is defensive against
   // future schema relaxation that would allow per-component grids.
   const DataField3D* comps[6] = {
      &sigma_xx_, &sigma_yy_, &sigma_zz_,
      &sigma_xy_, &sigma_yz_, &sigma_xz_
   };
   bbox_ = comps[0]->BBox();
   for (int c = 1; c < 6; ++c)
   {
      const auto& bc = comps[c]->BBox();
      // Axis-aligned intersection: (max of mins, min of maxes).
      bbox_[0] = std::max(bbox_[0], bc[0]);
      bbox_[1] = std::min(bbox_[1], bc[1]);
      bbox_[2] = std::max(bbox_[2], bc[2]);
      bbox_[3] = std::min(bbox_[3], bc[3]);
      bbox_[4] = std::max(bbox_[4], bc[4]);
      bbox_[5] = std::min(bbox_[5], bc[5]);
   }

   AssertConsistentGrid_();
}


void StressField3D::AssertConsistentGrid_() const
{
   const DataField3D* components[6] = {
      &sigma_xx_, &sigma_yy_, &sigma_zz_,
      &sigma_xy_, &sigma_yz_, &sigma_xz_
   };

   const int nx0 = sigma_xx_.NumX();
   const int ny0 = sigma_xx_.NumY();
   const int nz0 = sigma_xx_.NumZ();
   const std::array<real_t, 6>& bbox0 = sigma_xx_.BBox();

   for (int c = 1; c < 6; ++c)
   {
      if (components[c]->NumX() != nx0 ||
          components[c]->NumY() != ny0 ||
          components[c]->NumZ() != nz0)
      {
         std::ostringstream oss;
         oss << "StressField3D: component '"
             << components[c]->FieldName()
             << "' has axis sizes (" << components[c]->NumX() << ", "
             << components[c]->NumY() << ", " << components[c]->NumZ()
             << ") which differ from sigma_xx (" << nx0 << ", "
             << ny0 << ", " << nz0
             << ").  All six components must share the schema-v1 grid.";
         MFEM_ABORT(oss.str());
      }

      const auto& bbox_c = components[c]->BBox();
      for (int k = 0; k < 6; ++k)
      {
         // Exact equality is the spec: the Phase 5 writer emits the
         // same (x, y, z) axes for every component.
         if (bbox_c[k] != bbox0[k])
         {
            std::ostringstream oss;
            oss << "StressField3D: component '"
                << components[c]->FieldName()
                << "' bbox differs from sigma_xx at index " << k
                << " (component=" << bbox_c[k] << ", sigma_xx="
                << bbox0[k] << ").  All six components must share "
                << "the schema-v1 grid.";
            MFEM_ABORT(oss.str());
         }
      }
   }
}


mfem::DenseMatrix StressField3D::Evaluate(real_t x, real_t y, real_t z) const
{
   mfem::DenseMatrix S(3);

   const real_t xx = sigma_xx_.Evaluate(x, y, z);
   const real_t yy = sigma_yy_.Evaluate(x, y, z);
   const real_t zz = sigma_zz_.Evaluate(x, y, z);
   const real_t xy = sigma_xy_.Evaluate(x, y, z);
   const real_t yz = sigma_yz_.Evaluate(x, y, z);
   const real_t xz = sigma_xz_.Evaluate(x, y, z);

   S(0, 0) = xx; S(1, 1) = yy; S(2, 2) = zz;
   S(0, 1) = xy; S(1, 0) = xy;
   S(0, 2) = xz; S(2, 0) = xz;
   S(1, 2) = yz; S(2, 1) = yz;

   return S;
}


const DataField3D& StressField3D::Field(int component_index) const
{
   switch (component_index)
   {
      case 0: return sigma_xx_;
      case 1: return sigma_yy_;
      case 2: return sigma_zz_;
      case 3: return sigma_xy_;
      case 4: return sigma_yz_;
      case 5: return sigma_xz_;
      default:
         MFEM_ABORT("StressField3D::Field: component_index "
                    << component_index << " out of range [0, 6)");
   }
   // Unreachable; suppress compiler warning.
   return sigma_xx_;
}


bool StressField3D::ContainsBBox(real_t xmin, real_t xmax,
                                 real_t ymin, real_t ymax,
                                 real_t zmin, real_t zmax,
                                 real_t eps) const
{
   // The shared-grid invariant means every component agrees, so we
   // gate via the first one.  The assertion in the ctor guarantees
   // this is safe.
   return sigma_xx_.ContainsBBox(xmin, xmax, ymin, ymax,
                                 zmin, zmax, eps);
}


void StressField3D::SetInterpMode(InterpMode mode)
{
   sigma_xx_.SetInterpMode(mode);
   sigma_yy_.SetInterpMode(mode);
   sigma_zz_.SetInterpMode(mode);
   sigma_xy_.SetInterpMode(mode);
   sigma_yz_.SetInterpMode(mode);
   sigma_xz_.SetInterpMode(mode);
}

} // namespace seas
} // namespace mfem
