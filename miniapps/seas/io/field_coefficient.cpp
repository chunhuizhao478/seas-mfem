// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// field_coefficient.cpp — implementation of FieldProjector.

#include "field_coefficient.hpp"

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>

namespace mfem
{
namespace seas
{

std::atomic<int> FieldProjector::call_count_{0};


// ----------------------------------------------------------------------
// Mesh bbox helpers
// ----------------------------------------------------------------------

void FieldProjector::ComputeMeshBBoxSerial(const mfem::Mesh& mesh,
                                           real_t& xmin, real_t& xmax,
                                           real_t& ymin, real_t& ymax,
                                           real_t& zmin, real_t& zmax)
{
   const real_t inf = std::numeric_limits<real_t>::infinity();
   xmin = ymin = zmin =  inf;
   xmax = ymax = zmax = -inf;
   const int nv = mesh.GetNV();
   for (int v = 0; v < nv; ++v)
   {
      const real_t* p = mesh.GetVertex(v);
      xmin = std::min(xmin, p[0]); xmax = std::max(xmax, p[0]);
      ymin = std::min(ymin, p[1]); ymax = std::max(ymax, p[1]);
      zmin = std::min(zmin, p[2]); zmax = std::max(zmax, p[2]);
   }
}

void FieldProjector::ComputeMeshBBoxParallel(const mfem::ParMesh& mesh,
                                             real_t& xmin, real_t& xmax,
                                             real_t& ymin, real_t& ymax,
                                             real_t& zmin, real_t& zmax)
{
   // REVIEW R-013: do not use `ParMesh::GetBoundingBox` here.  That
   // method is non-const and would force a `const_cast` at the call
   // site, breaking the const promise on `mesh` and creating a
   // maintenance hazard (a future change to `Mesh::GetBoundingBox`
   // that mutates observable state would silently corrupt callers).
   //
   // Instead, compute the local bbox via the existing serial helper
   // (which only reads vertex coords) and do two clean MPI_Allreduce
   // calls (one MIN, one MAX) over packed 3-component buffers.  This
   // is both const-correct and avoids the previous hand-rolled
   // packing trick that mixed MIN/MAX into a single MAX-reduce.
   ComputeMeshBBoxSerial(mesh, xmin, xmax, ymin, ymax, zmin, zmax);

   real_t local_min[3] = { xmin, ymin, zmin };
   real_t local_max[3] = { xmax, ymax, zmax };
   real_t global_min[3];
   real_t global_max[3];
   MPI_Allreduce(local_min, global_min, 3,
                 MPITypeMap<real_t>::mpi_type, MPI_MIN, mesh.GetComm());
   MPI_Allreduce(local_max, global_max, 3,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, mesh.GetComm());
   xmin = global_min[0]; xmax = global_max[0];
   ymin = global_min[1]; ymax = global_max[1];
   zmin = global_min[2]; zmax = global_max[2];

   // R-009 (review round-2): if every rank had zero local vertices,
   // the per-rank serial bbox is (+inf, -inf, ...) and the global
   // reductions produce a degenerate bbox.  `ContainsBBox` would
   // then trivially pass (+inf >= anything is TRUE, -inf <= anything
   // is TRUE), masking the real problem.  Abort with a clear message
   // instead.
   if (xmin > xmax || ymin > ymax || zmin > zmax)
   {
      MFEM_ABORT(
         "FieldProjector::ComputeMeshBBoxParallel: global mesh bbox "
         "is degenerate.  Global bounds: "
         "x=[" << xmin << ", " << xmax << "], "
         "y=[" << ymin << ", " << ymax << "], "
         "z=[" << zmin << ", " << zmax << "].  "
         "(A degenerate axis has min > max.)  This typically indicates "
         "the mesh has zero elements on every MPI rank — check that the "
         "mesh file is non-empty and that the ParMesh partitioning did "
         "not silently drop all elements.");
   }
}


// ----------------------------------------------------------------------
// GF min/max helpers
// ----------------------------------------------------------------------

void FieldProjector::ComputeGridFunctionMinMaxSerial(
   const mfem::GridFunction& gf, real_t& lo, real_t& hi)
{
   lo =  std::numeric_limits<real_t>::infinity();
   hi = -std::numeric_limits<real_t>::infinity();
   const int n = gf.Size();
   for (int i = 0; i < n; ++i)
   {
      const real_t v = gf(i);
      if (v < lo) { lo = v; }
      if (v > hi) { hi = v; }
   }
}

void FieldProjector::ComputeGridFunctionMinMaxParallel(
   const mfem::ParGridFunction& gf, real_t& lo, real_t& hi)
{
   ComputeGridFunctionMinMaxSerial(gf, lo, hi);
   real_t local[2] = {-lo, hi};
   real_t global[2];
   MPI_Allreduce(local, global, 2, MPITypeMap<real_t>::mpi_type,
                 MPI_MAX, gf.ParFESpace()->GetComm());
   lo = -global[0]; hi = global[1];
}


// ----------------------------------------------------------------------
// Abort helpers
// ----------------------------------------------------------------------

void FieldProjector::AbortContainmentFailure(
   const DataField3D& field,
   real_t mxmin, real_t mxmax,
   real_t mymin, real_t mymax,
   real_t mzmin, real_t mzmax)
{
   const auto& fb = field.BBox();
   std::ostringstream oss;
   oss << "FieldProjector::Project: mesh bbox is NOT contained in "
       << "field '" << field.FieldName() << "' data bbox.\n"
       << "  mesh bbox  (UTM 11 N, m): x=[" << mxmin << ", " << mxmax
       << "] y=[" << mymin << ", " << mymax << "] z=[" << mzmin << ", "
       << mzmax << "]\n"
       << "  data bbox  (UTM 11 N, m): x=[" << fb[0] << ", " << fb[1]
       << "] y=[" << fb[2] << ", " << fb[3] << "] z=[" << fb[4] << ", "
       << fb[5] << "]\n"
       << "v1 schema enforces interpolation-only; either shrink the "
       << "mesh padding or regenerate the sidecar over a wider region.";
   MFEM_ABORT(oss.str());
}

void FieldProjector::AbortRangeFailure(
   const std::string& field_name,
   real_t observed_lo, real_t observed_hi,
   real_t declared_lo, real_t declared_hi)
{
   std::ostringstream oss;
   oss << "FieldProjector::Project: projected field '" << field_name
       << "' violates the declared sanity bounds.\n"
       << "  observed   = [" << observed_lo << ", " << observed_hi << "]\n"
       << "  declared   = [" << declared_lo << ", " << declared_hi << "]\n"
       << "Either widen the bounds in the sidecar (and reload) or "
       << "regenerate the source dataset.";
   MFEM_ABORT(oss.str());
}


// ----------------------------------------------------------------------
// Project / ProjectSerial
// ----------------------------------------------------------------------

FieldProjector::Result FieldProjector::ProjectSerial(
   const DataField3D&            field,
   mfem::FiniteElementSpace&     target_fes,
   real_t scale, real_t offset)
{
   call_count_.fetch_add(1);

   real_t mxmin, mxmax, mymin, mymax, mzmin, mzmax;
   const mfem::Mesh* mesh = target_fes.GetMesh();
   MFEM_VERIFY(mesh != nullptr,
               "FieldProjector: target_fes has null mesh");
   ComputeMeshBBoxSerial(*mesh, mxmin, mxmax, mymin, mymax, mzmin, mzmax);
   if (!field.ContainsBBox(mxmin, mxmax, mymin, mymax, mzmin, mzmax))
   {
      AbortContainmentFailure(field,
                              mxmin, mxmax, mymin, mymax, mzmin, mzmax);
   }

   FieldCoefficient coef(field, scale, offset);
   auto gf = std::make_shared<mfem::GridFunction>(&target_fes);
   gf->ProjectCoefficient(coef);

   real_t lo, hi;
   ComputeGridFunctionMinMaxSerial(*gf, lo, hi);

   const real_t expected_lo = scale * field.MinValue() + offset;
   const real_t expected_hi = scale * field.MaxValue() + offset;
   const real_t lo_band = std::min(expected_lo, expected_hi);
   const real_t hi_band = std::max(expected_lo, expected_hi);
   if (lo < lo_band || hi > hi_band)
   {
      AbortRangeFailure(field.FieldName(), lo, hi, lo_band, hi_band);
   }

   Result r;
   r.gf = gf;
   r.min_value = lo;
   r.max_value = hi;
   return r;
}


FieldProjector::Result FieldProjector::Project(
   const DataField3D&             field,
   mfem::ParFiniteElementSpace&   target_fes,
   real_t scale, real_t offset)
{
   call_count_.fetch_add(1);

   const mfem::ParMesh* pmesh = target_fes.GetParMesh();
   MFEM_VERIFY(pmesh != nullptr,
               "FieldProjector::Project: target_fes has null ParMesh");

   real_t mxmin, mxmax, mymin, mymax, mzmin, mzmax;
   ComputeMeshBBoxParallel(*pmesh, mxmin, mxmax, mymin, mymax, mzmin, mzmax);
   if (!field.ContainsBBox(mxmin, mxmax, mymin, mymax, mzmin, mzmax))
   {
      AbortContainmentFailure(field,
                              mxmin, mxmax, mymin, mymax, mzmin, mzmax);
   }

   FieldCoefficient coef(field, scale, offset);
   auto gf = std::make_shared<mfem::ParGridFunction>(&target_fes);
   gf->ProjectCoefficient(coef);

   real_t lo, hi;
   ComputeGridFunctionMinMaxParallel(*gf, lo, hi);

   const real_t expected_lo = scale * field.MinValue() + offset;
   const real_t expected_hi = scale * field.MaxValue() + offset;
   const real_t lo_band = std::min(expected_lo, expected_hi);
   const real_t hi_band = std::max(expected_lo, expected_hi);

   // Catmull-Rom is a non-monotonic cubic; it can over/undershoot
   // by ~ 5% of the local range at C^0 jumps in the source (e.g.
   // basin/basement Vs contrast).  Relax the per-source range check
   // by 5 % of the band width when CatmullRom is active so the
   // interp-only contract isn't violated by a known property of the
   // chosen scheme.  Trilinear keeps the original strict check.
   real_t lo_check = lo_band;
   real_t hi_check = hi_band;
   if (field.GetInterpMode() == InterpMode::CatmullRom)
   {
      const real_t slack = real_t{0.05} * (hi_band - lo_band);
      lo_check = lo_band - slack;
      hi_check = hi_band + slack;
   }
   if (lo < lo_check || hi > hi_check)
   {
      AbortRangeFailure(field.FieldName(), lo, hi, lo_check, hi_check);
   }

   Result r;
   r.gf = gf;
   r.min_value = lo;
   r.max_value = hi;
   return r;
}


// ----------------------------------------------------------------------
// ProjectVelocity
// ----------------------------------------------------------------------

FieldProjector::VelocityFields FieldProjector::ProjectVelocity(
   const std::string&             sidecar_path,
   mfem::ParFiniteElementSpace&   target_fes,
   real_t lambda_min_pa, real_t lambda_max_pa,
   real_t mu_min_pa,     real_t mu_max_pa,
   InterpMode interp)
{
   DataField3D vp(sidecar_path, "Vp");
   DataField3D vs(sidecar_path, "Vs");
   DataField3D rho(sidecar_path, "density");
   vp.SetInterpMode(interp);
   vs.SetInterpMode(interp);
   rho.SetInterpMode(interp);

   Result vp_r = Project(vp, target_fes);
   Result vs_r = Project(vs, target_fes);
   Result rho_r = Project(rho, target_fes);

   auto rho_pgf = std::dynamic_pointer_cast<mfem::ParGridFunction>(rho_r.gf);
   auto vp_pgf  = std::dynamic_pointer_cast<mfem::ParGridFunction>(vp_r.gf);
   auto vs_pgf  = std::dynamic_pointer_cast<mfem::ParGridFunction>(vs_r.gf);
   MFEM_VERIFY(rho_pgf && vp_pgf && vs_pgf,
               "FieldProjector::ProjectVelocity: parallel project did "
               "not return ParGridFunction (internal error)");
   const int n = rho_pgf->Size();
   MFEM_VERIFY(vp_pgf->Size() == n && vs_pgf->Size() == n,
               "FieldProjector::ProjectVelocity: Vp/Vs/rho GFs have "
               "different sizes (Vp=" << vp_pgf->Size() << ", Vs="
               << vs_pgf->Size() << ", rho=" << n << ")");

   auto mu = std::make_shared<mfem::ParGridFunction>(&target_fes);
   auto la = std::make_shared<mfem::ParGridFunction>(&target_fes);
   for (int i = 0; i < n; ++i)
   {
      const real_t r  = (*rho_pgf)(i);
      const real_t cp = (*vp_pgf)(i);
      const real_t cs = (*vs_pgf)(i);
      const real_t mu_i = r * cs * cs;
      const real_t la_i = r * cp * cp - 2.0 * mu_i;
      (*mu)(i) = mu_i;
      (*la)(i) = la_i;
   }

   real_t mu_lo, mu_hi, la_lo, la_hi;
   ComputeGridFunctionMinMaxParallel(*mu, mu_lo, mu_hi);
   ComputeGridFunctionMinMaxParallel(*la, la_lo, la_hi);

   // 5% slack for CatmullRom over/undershoot (see Project() above).
   real_t mu_lo_check = mu_min_pa, mu_hi_check = mu_max_pa;
   real_t la_lo_check = lambda_min_pa, la_hi_check = lambda_max_pa;
   if (interp == InterpMode::CatmullRom)
   {
      const real_t mu_slack = real_t{0.05} * (mu_max_pa - mu_min_pa);
      const real_t la_slack = real_t{0.05} * (lambda_max_pa - lambda_min_pa);
      mu_lo_check = mu_min_pa - mu_slack;
      mu_hi_check = mu_max_pa + mu_slack;
      la_lo_check = lambda_min_pa - la_slack;
      la_hi_check = lambda_max_pa + la_slack;
   }
   if (mu_lo < mu_lo_check || mu_hi > mu_hi_check)
   {
      AbortRangeFailure("derived mu", mu_lo, mu_hi,
                        mu_lo_check, mu_hi_check);
   }
   if (la_lo < la_lo_check || la_hi > la_hi_check)
   {
      AbortRangeFailure("derived lambda", la_lo, la_hi,
                        la_lo_check, la_hi_check);
   }

   VelocityFields vf;
   vf.rho     = rho_pgf;
   vf.lambda  = la;
   vf.mu      = mu;
   vf.vp      = vp_pgf;
   vf.vs      = vs_pgf;
   vf.min_vp  = vp_r.min_value;  vf.max_vp  = vp_r.max_value;
   vf.min_vs  = vs_r.min_value;  vf.max_vs  = vs_r.max_value;
   vf.min_rho = rho_r.min_value; vf.max_rho = rho_r.max_value;
   vf.min_lambda = la_lo; vf.max_lambda = la_hi;
   vf.min_mu     = mu_lo; vf.max_mu     = mu_hi;
   return vf;
}

} // namespace seas
} // namespace mfem
