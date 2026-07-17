// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// spatial_velocity.cpp — Phase 2 implementation.

#include "spatial_velocity.hpp"

#include "../../io/field_coefficient.hpp"   // FieldProjector::AbortContainmentFailure

#include <fstream>
#include <string>

namespace mfem
{
namespace seas
{
namespace spatial
{

namespace
{

const char* model_dir(VelocityModel m)
{
   switch (m)
   {
   case VelocityModel::CVMH:                return "cvmh";
   case VelocityModel::CVMS_4_26_M01:       return "cvm_s4.26.m01";
   case VelocityModel::MultiscaleStatewise: return "multiscale_statewise_cvm";
   }
   MFEM_ABORT("VelocityModel enum out of range");
   return "";
}

bool file_exists(const std::string& path)
{
   std::ifstream f(path);
   return f.good();
}

template <typename MeshT>
void compute_bbox(MeshT& mesh,
                  real_t& xmin, real_t& xmax,
                  real_t& ymin, real_t& ymax,
                  real_t& zmin, real_t& zmax)
{
   // GetBoundingBox is defined on both Mesh and ParMesh; for ParMesh it
   // returns the global bbox via MPI_Allreduce internally (per MFEM
   // convention).
   mfem::Vector lo, hi;
   mesh.GetBoundingBox(lo, hi);
   MFEM_VERIFY(lo.Size() >= 3 && hi.Size() >= 3,
               "compute_bbox: mesh dimension < 3");
   xmin = lo(0); xmax = hi(0);
   ymin = lo(1); ymax = hi(1);
   zmin = lo(2); zmax = hi(2);
}

template <typename MeshT>
SpatialVelocityBundle load_impl(const VelocitySpec& spec, MeshT& mesh
#ifdef MFEM_USE_MPI
                                , MPI_Comm node_comm = MPI_COMM_NULL
#endif
                                )
{
   const std::string path = ResolveSpatialVelocitySidecarPath(spec);
   MFEM_VERIFY(file_exists(path),
               "LoadSpatialVelocityBundle: velocity sidecar file not "
               "found at '" << path << "'.  Check `[velocity].model` and "
               "`[velocity].dataset_root` (or `override_path`).");

   // Opt-in ASAGI-style far-field clamp: when the mesh extends beyond the
   // velocity data hull (e.g. safv4_deep's large absorbing box), edge-clamp
   // out-of-hull material queries instead of aborting.  Default = Abort
   // preserves the strict interpolation-only contract (TPV/BP5 unchanged).
   const OOBPolicy oob = spec.far_field_clamp ? OOBPolicy::Clamp
                                              : OOBPolicy::Abort;

   SpatialVelocityBundle b;
#ifdef MFEM_USE_MPI
   if (node_comm != MPI_COMM_NULL)
   {
      // MPI-3 shared-memory mode: one physical copy of each grid per
      // node (PLAN_sidecar_mpi_shared_memory_2026-07-17.md §4).
      b.vp_field  = std::make_unique<DataField3D>(path, "Vp",      oob,
                                                  node_comm);
      b.vs_field  = std::make_unique<DataField3D>(path, "Vs",      oob,
                                                  node_comm);
      b.rho_field = std::make_unique<DataField3D>(path, "density", oob,
                                                  node_comm);
   }
   else
#endif
   {
      b.vp_field  = std::make_unique<DataField3D>(path, "Vp",      oob);
      b.vs_field  = std::make_unique<DataField3D>(path, "Vs",      oob);
      b.rho_field = std::make_unique<DataField3D>(path, "density", oob);
   }

   real_t mxmin, mxmax, mymin, mymax, mzmin, mzmax;
   compute_bbox(mesh, mxmin, mxmax, mymin, mymax, mzmin, mzmax);

   if (spec.far_field_clamp)
   {
      // Far-field clamp ON: the ContainsBBox gate is intentionally skipped
      // (an out-of-hull mesh is the whole point); report the breach so the
      // edge-clamp is visible in the run log, mirroring SeisSol/ASAGI.
      const std::array<real_t, 6>& vb = b.vp_field->BBox();
      mfem::out << "LoadSpatialVelocityBundle: [velocity].far_field_clamp = "
                   "true -> ASAGI-style nearest-edge hold for mesh extent "
                   "outside the velocity data hull (ContainsBBox gate skipped).\n"
                << "  mesh bbox x[" << mxmin << "," << mxmax << "] y["
                << mymin << "," << mymax << "] z[" << mzmin << "," << mzmax
                << "]\n  data bbox x[" << vb[0] << "," << vb[1] << "] y["
                << vb[2] << "," << vb[3] << "] z[" << vb[4] << "," << vb[5]
                << "] (Vp; Vs/density share the grid)\n";
   }
   else
   {
      // ContainsBBox on all three fields.  Delegate to the shared
      // FieldProjector::AbortContainmentFailure helper so the error
      // message stays in sync with the rest of the SEAS I/O layer (R-012;
      // helper exposed as public on field_coefficient.hpp).
      auto check = [&](const DataField3D& field)
      {
         const bool inside = field.ContainsBBox(mxmin, mxmax,
                                                mymin, mymax,
                                                mzmin, mzmax);
         if (!inside)
         {
            FieldProjector::AbortContainmentFailure(
               field, mxmin, mxmax, mymin, mymax, mzmin, mzmax);
         }
      };
      check(*b.vp_field);
      check(*b.vs_field);
      check(*b.rho_field);
   }

   b.lambda_coef = std::make_unique<LambdaFromSidecar>(
                      *b.vp_field, *b.vs_field, *b.rho_field);
   b.mu_coef     = std::make_unique<MuFromSidecar>(
                      *b.vs_field, *b.rho_field);
   b.rho_coef    = std::make_unique<RhoFromSidecar>(*b.rho_field);

   return b;
}

}  // namespace

std::string ResolveSpatialVelocitySidecarPath(const VelocitySpec& spec)
{
   if (!spec.override_path.empty()) { return spec.override_path; }
   MFEM_VERIFY(!spec.dataset_root.empty(),
               "ResolveSpatialVelocitySidecarPath: dataset_root must be "
               "non-empty when override_path is empty.");
   std::string path = spec.dataset_root;
   if (!path.empty() && path.back() != '/') { path.push_back('/'); }
   path += "velocity/results/";
   path += model_dir(spec.model);
   path += "/velocity_safs.h5";
   return path;
}

SpatialVelocityBundle LoadSpatialVelocityBundle(const VelocitySpec& spec,
                                                ParMesh&            pmesh)
{
   return load_impl<ParMesh>(spec, pmesh);
}

#ifdef MFEM_USE_MPI
SpatialVelocityBundle LoadSpatialVelocityBundle(const VelocitySpec& spec,
                                                ParMesh&            pmesh,
                                                MPI_Comm            node_comm)
{
   return load_impl<ParMesh>(spec, pmesh, node_comm);
}
#endif

SpatialVelocityBundle LoadSpatialVelocityBundle(const VelocitySpec& spec,
                                                mfem::Mesh&         mesh)
{
   return load_impl<mfem::Mesh>(spec, mesh);
}

}  // namespace spatial
}  // namespace seas
}  // namespace mfem
