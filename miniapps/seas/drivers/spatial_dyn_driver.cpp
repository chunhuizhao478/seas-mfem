// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// seas_spatial_dyn_driver — Phase 4 of
// safs/project_7.0_alternative/document/spatial_dynamic_rupture_plan.md (rev-3).
//
// Single-event dynamic-rupture driver for spatially-varying material,
// pre-stress, and friction inputs.  Reuses the TPV205 LSW closed-form
// machinery (Tpv205SubStepIterator + AdvanceADERWithSubStep + the new
// LSW_ForcedRupture dispatch from Phase H) on a SAFS / TPV26-27-style
// problem driven by a single TOML config.
//
// Deviation block (vs Phase 4 spec).  Documented up front for the
// adversarial code review (REVIEW.md).
//
// D-1.  WaveOperator(MaterialField, BoundaryConfig) ctor (Phase H.1)
//       is NOT yet wired (see dynamic/wave_operator.hpp:215-233,
//       SetGodunovFluxPool still aborts).  Until Phase H lands the
//       per-element flux dispatch, this driver constructs WaveOperator
//       through the scalar-material ctor only:
//
//         * MaterialField::Mode::Constant  ⇒ extract (lambda, mu, rho)
//           and call the scalar ctor (byte-equivalent to the
//           heterogeneous ctor on a Constant input).
//         * MaterialField::Mode::Coefficient (the SAFS sidecar path) ⇒
//           ABORT with a clear "Phase H not yet wired" message so the
//           caller cannot accidentally run sidecar-input physics on
//           the scalar code path.
//
//       --no-sidecar-material always uses the Constant path.  The CSM
//       stress sidecar path (Phase 3) is unaffected — pre-stress is a
//       *fault-DOF* quantity and runs through FaultGeometry, not the
//       WaveOperator material.
//
// D-2.  spatial::ComputePerDOFCoordsAndBasisFromWave (Phase 5 helper
//       was deferred — see dynamic/spatial_setup.hpp:14-22).  The
//       driver therefore mirrors the TPV205 driver's inline fault-DOF
//       walk to build per-DOF coordinates / basis / dof_to_elem from
//       wave.GetFaultInteriorFaces() + GetFaultSharedFaces().
//
// D-3.  Write/ReadTpv104Checkpoint actual signature is
//       (prefix, t, dt, step, Q, dof_data, mpi/(rank,size,comm),
//        driver_tag).  The plan's docstring referenced a hypothetical
//       6-extra-ParaView-state form; the actual on-disk schema is
//       smaller.  This driver uses the real raw-MPI overload with
//       driver_tag = "spatial_dyn".

#include "mfem.hpp"

#include "../dynamic/wave_state.hpp"
#include "../dynamic/wave_operator.hpp"
#include "../dynamic/bimaterial_wave_operator.hpp"  // Phase 13: matrix (bimaterial) operator
#include "../dynamic/lts_clustering.hpp"             // LTS Phase 0: --lts-report clustering
#include "../dynamic/lts_layout.hpp"                 // LTS Phase 1: run-side layout
#include "../dynamic/lts_partition.hpp"              // LTS Phase 4 Step 0: cluster-weighted METIS partition
#include "../dynamic/lts_bulk_stepper.hpp"           // LTS Phase 2: sync-interval stepper
#include "../dynamic/lts_fault_stepper.hpp"          // LTS Phase 3 (P3-4): fault-aware sync-interval stepper
#include <cstdint>
#include "../dynamic/fault_face_flux.hpp"
#include "../dynamic/friction_solver.hpp"
#include "../dynamic/tpv205_friction.hpp"
#include "../dynamic/fault_state_channel.hpp"        // R-024: shared state-channel value
#include "../dynamic/tpv205_substep_iterator.hpp"
#include "../dynamic/friction_iterator.hpp"           // Phase 2: IFrictionIterator + adapters
#include "../dynamic/friction_iterator_factory.hpp"   // Phase 3: MakeFrictionIterator
#include "../dynamic/fault_resample.hpp"              // Phase 3: BuildFaultResampleMatrix
#include "../dynamic/heterogeneous_material.hpp"
#include "../dynamic/spatial_setup.hpp"
#include "../dynamic/seas_diag_rank.hpp"

#include "../domain/boundary_config.hpp"

#include "../fault/fault_basis.hpp"
#include "../fault/fault_geometry.hpp"
#include "../fault/fault_geometry_safs_templated.inl"

#include "../config/bp5_params.hpp"

#include "../common/mpi_context.hpp"

#include "../io/paraview_output.hpp"
#include "../io/free_surface_output.hpp"
#include "../io/tpv104_checkpoint.hpp"
#include "../io/data_field_3d.hpp"
#include "../io/stress_field_3d.hpp"
#include "../io/material_coefficients.hpp"

#include "../dynamic/spatial_nucleation.hpp"
#include "../dynamic/nucleation_factory.hpp"   // Phase 7: MakeNucleation
#include "../dynamic/rk_time_stepper.hpp"      // Phase 14: RK4/RK45 coupled stepper
#include "../dynamic/tpv31_stations.hpp"       // Phase 10: TPV31 SCEC station writer
#include "../dynamic/tpv102_stations.hpp"      // SCEC TPV102 station writer (rate-state)
#include "../dynamic/tpv104_stations.hpp"      // SCEC TPV104 station writer (rate-state)
#include "../dynamic/tpv205_stations.hpp"      // SCEC TPV205 station writer (LSW)
#include "../dynamic/tpv6_stations.hpp"        // Part C: TPV6/7 per-side on-fault station writer
#include "../dynamic/fault_locality_partition.hpp"  // keep fault element-pairs co-resident (np>1)
#include "../dynamic/spatial_print_derived.hpp"

#include "../spatial/code/spatial_friction.hpp"
#include "../spatial/code/spatial_velocity.hpp"
#include "../spatial/code/spatial_stress.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <functional>
#include <memory>
#include <type_traits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

#ifdef SEAS_DIAG_FAULT_FLUX
namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT
#endif

// Precondition P-2 build-time guard (plan §Preconditions, rev-3).
// MaterialField::EvalAt must accept (elem, T, ip, lam, mu, rho) by ref.
static_assert(
   std::is_invocable_v<decltype(&mfem::seas::MaterialField::EvalAt),
                       const mfem::seas::MaterialField*, int,
                       mfem::ElementTransformation&,
                       const mfem::IntegrationPoint&,
                       mfem::real_t&, mfem::real_t&, mfem::real_t&>,
   "Precondition P-2 not met: MaterialField::EvalAt signature drift.");

// --------------------------------------------------------------------------
// Small CLI parsing helpers (same convention as tpv104_driver.cpp).
// --------------------------------------------------------------------------
namespace
{

std::string GetStringArg(int argc, char *argv[], const char *flag,
                         const std::string &default_val)
{
   for (int i = 1; i < argc - 1; ++i)
   {
      if (std::string(argv[i]) == flag) { return argv[i + 1]; }
   }
   return default_val;
}

real_t GetRealArg(int argc, char *argv[], const char *flag, real_t default_val)
{
   const std::string val = GetStringArg(argc, argv, flag, "");
   if (val.empty()) { return default_val; }
   return std::stod(val);
}

int GetIntArg(int argc, char *argv[], const char *flag, int default_val)
{
   const std::string val = GetStringArg(argc, argv, flag, "");
   if (val.empty()) { return default_val; }
   return std::stoi(val);
}

bool HasFlag(int argc, char *argv[], const char *flag)
{
   for (int i = 1; i < argc; ++i)
   {
      if (std::string(argv[i]) == flag) { return true; }
   }
   return false;
}

// (LTS Phase 0/1) Clustering inputs derived from a CONSTRUCTED wave operator +
// its ParMesh: per-local-element CFL dt (only the ratios matter to the cluster
// assignment; the absolute scale is recovered by the caller when it wants
// seconds), the local interior fault-face element pairs (forced diff 0), and
// every local face as an LtsFaceSpec for BuildLtsLayout.  Shared by
// --lts-report and the run-path layout wiring so the dt / fault-face definition
// is single-sourced.
struct LtsMeshInputs
{
   std::vector<double>                  dt_e;         // per local element (h_e / c_p,e)
   std::vector<std::pair<int, int>>     fault_pairs;  // local fault-face elem pairs
   std::vector<mfem::seas::LtsFaceSpec> faces;        // every local face
   bool per_elem_material = false;                    // true on the matrix path
};

LtsMeshInputs BuildLtsMeshInputs(WaveOperator<ParMesh> &wave, ParMesh &pmesh)
{
   LtsMeshInputs in;
   const int neL = pmesh.GetNE();
   const std::vector<real_t> &hvec = wave.GetPerElementCflLength();
   const std::vector<std::array<real_t, 3>> &lmr = wave.GetPerElementMaterial();
   MFEM_VERIFY(static_cast<int>(hvec.size()) == neL,
               "BuildLtsMeshInputs: per-element CFL length size " << hvec.size()
               << " != local NE " << neL << " (operator not built?)");
   in.per_elem_material = !lmr.empty();
   MFEM_VERIFY(!in.per_elem_material || static_cast<int>(lmr.size()) == neL,
               "BuildLtsMeshInputs: per-element material size " << lmr.size()
               << " != local NE " << neL);

   in.dt_e.resize(static_cast<std::size_t>(neL));
   for (int e = 0; e < neL; ++e)
   {
      double cp = 1.0;  // scalar path: uniform c_p cancels out of cluster ratios
      if (in.per_elem_material)
      {
         const double lam = lmr[e][0], mu = lmr[e][1], rho = lmr[e][2];
         cp = std::sqrt((lam + 2.0 * mu) / rho);
      }
      in.dt_e[static_cast<std::size_t>(e)] =
         static_cast<double>(hvec[e]) / cp;
   }

   // Fault faces: mark the wave operator's LOCAL interior fault-face indices so
   // BuildLtsLayout tracks them as diff-0 and BuildLtsClustering forces the pair
   // into one cluster.
   std::vector<char> is_fault_face(pmesh.GetNumFaces(), 0);
   const Array<int> &fif = wave.GetFaultInteriorFaces();
   for (int i = 0; i < fif.Size(); ++i)
   {
      const int f = fif[i];
      if (f >= 0 && f < pmesh.GetNumFaces()) { is_fault_face[f] = 1; }
      int e1 = -1, e2 = -1;
      pmesh.GetFaceElements(f, &e1, &e2);
      if (e1 >= 0 && e2 >= 0) { in.fault_pairs.emplace_back(e1, e2); }
   }

   const int nfaces = pmesh.GetNumFaces();
   in.faces.reserve(static_cast<std::size_t>(nfaces));
   for (int f = 0; f < nfaces; ++f)
   {
      int e1 = -1, e2 = -1;
      pmesh.GetFaceElements(f, &e1, &e2);
      if (e1 < 0) { continue; }   // face with no local element (skip)
      mfem::seas::LtsFaceSpec fs;
      fs.face_id  = f;
      fs.elem1    = e1;
      fs.elem2    = e2;   // < 0 for boundary / rank seam
      fs.is_fault = (is_fault_face[f] != 0);
      in.faces.push_back(fs);
   }
   return in;
}

// (LTS Phase 4 Step 0) Operator-INDEPENDENT clustering inputs, sourced from the
// mesh geometry + the MaterialField directly — so the cluster ids can be
// computed BEFORE the wave operator is constructed and passed into its ctor for
// the fault-QP reorder (Phase 3, P-006).  Produces the SAME per-element dt_e
// binning as BuildLtsMeshInputs above:
//   dt_e = h_e / c_p,e,  h_e = insphere diameter (tets) / vol^{1/dim},
//   c_p,e = sqrt((lambda + 2 mu) / rho) evaluated at the element centroid.
// For Mode::Constant (the scalar-flux path) c_p is uniform, so — like the
// operator-sourced builder's cp=1 — it leaves the cluster binning unchanged
// (the assignment is invariant to a uniform scale of dt_e); for Mode::Coefficient
// (the matrix/bimaterial path) c_p is the true per-element value, matching the
// operator's per_elem_lmr_ EvalAt-at-centroid.  `fault_attr <= 0` => no fault
// faces.  Runs on the ParMesh at np==1 (where local == global element order);
// the deterministic serial-mesh version is the np>1 half of this step.
// Takes a base `Mesh&` so it serves BOTH the local ParMesh (np=1 reorder
// activation) AND the SERIAL mesh (np>1 serial clustering, Phase 4 Step 0).
LtsMeshInputs BuildLtsMeshInputsFromMaterial(Mesh &pmesh,
                                             const MaterialField &material,
                                             int fault_attr)
{
   LtsMeshInputs in;
   const int neL = pmesh.GetNE();
   in.per_elem_material = (material.mode != MaterialField::Mode::Constant);
   in.dt_e.resize(static_cast<std::size_t>(neL));

   for (int e = 0; e < neL; ++e)
   {
      const real_t vol = pmesh.GetElementVolume(e);
      const Geometry::Type gtype = pmesh.GetElementBaseGeometry(e);
      real_t h;
      if (gtype == Geometry::TETRAHEDRON)
      {
         Array<int> vert; pmesh.GetElementVertices(e, vert);
         Vector v0(pmesh.GetVertex(vert[0]), 3), v1(pmesh.GetVertex(vert[1]), 3);
         Vector v2(pmesh.GetVertex(vert[2]), 3), v3(pmesh.GetVertex(vert[3]), 3);
         auto tri_area = [](const Vector &a, const Vector &b,
                            const Vector &c) -> real_t {
            real_t e1[3] = {b(0)-a(0), b(1)-a(1), b(2)-a(2)};
            real_t e2[3] = {c(0)-a(0), c(1)-a(1), c(2)-a(2)};
            real_t cx = e1[1]*e2[2] - e1[2]*e2[1];
            real_t cy = e1[2]*e2[0] - e1[0]*e2[2];
            real_t cz = e1[0]*e2[1] - e1[1]*e2[0];
            return 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
         };
         const real_t A = tri_area(v0,v1,v2) + tri_area(v0,v2,v3)
                        + tri_area(v0,v3,v1) + tri_area(v1,v2,v3);
         h = (A > 0) ? 6.0 * vol / A : std::pow(vol, 1.0/3.0);
      }
      else
      {
         h = std::pow(vol, 1.0 / pmesh.Dimension());
      }

      double cp = 1.0;
      if (in.per_elem_material)
      {
         // REVIEW A-3: this LTS path evaluates the material as a coordinate
         // COEFFICIENT (EvalAt); Mode::GridFunction would abort inside EvalAt.
         // The spatial driver only ever builds Constant/Coefficient material, so
         // guard loudly rather than aborting deep in EvalAt.
         MFEM_VERIFY(material.mode != MaterialField::Mode::GridFunction,
                     "BuildLtsMeshInputsFromMaterial: Mode::GridFunction is not "
                     "supported on the pre-operator LTS clustering path (it needs "
                     "a coordinate-evaluable coefficient).");
         ElementTransformation *T = pmesh.GetElementTransformation(e);
         const IntegrationPoint &ip = Geometries.GetCenter(gtype);
         // REVIEW A-2: SET the integration point before EvalAt, mirroring the
         // operator's MaterialAtLocal_ (P2-005) — a coefficient that reads
         // T.GetIntPoint() would otherwise see a stale point and yield a
         // different c_p (hence a different dt_e binning than the operator).
         T->SetIntPoint(&ip);
         real_t lam = 0, mu = 0, rho = 0;
         material.EvalAt(e, *T, ip, lam, mu, rho);
         MFEM_VERIFY(rho > 0.0 && (lam + 2.0*mu) > 0.0,
                     "BuildLtsMeshInputsFromMaterial: non-physical material at "
                     "element " << e << " (rho=" << rho << ", lambda+2mu="
                     << (lam + 2.0*mu) << ").");
         cp = std::sqrt((static_cast<double>(lam) + 2.0*mu) / rho);
      }
      in.dt_e[static_cast<std::size_t>(e)] = static_cast<double>(h) / cp;
   }

   // Fault faces from the mesh boundary (fault_attr); both elements of each
   // 2-sided fault face form a fault pair (forced diff-0 by BuildLtsClustering).
   std::vector<char> is_fault_face(pmesh.GetNumFaces(), 0);
   if (fault_attr > 0)
   {
      for (int b = 0; b < pmesh.GetNBE(); ++b)
      {
         if (pmesh.GetBdrAttribute(b) != fault_attr) { continue; }
         const int f = pmesh.GetBdrElementFaceIndex(b);
         if (f < 0 || f >= pmesh.GetNumFaces()) { continue; }
         int e1 = -1, e2 = -1;
         pmesh.GetFaceElements(f, &e1, &e2);
         if (e1 >= 0 && e2 >= 0)
         { is_fault_face[f] = 1; in.fault_pairs.emplace_back(e1, e2); }
      }
   }
   const int nfaces = pmesh.GetNumFaces();
   in.faces.reserve(static_cast<std::size_t>(nfaces));
   for (int f = 0; f < nfaces; ++f)
   {
      int e1 = -1, e2 = -1;
      pmesh.GetFaceElements(f, &e1, &e2);
      if (e1 < 0) { continue; }
      mfem::seas::LtsFaceSpec fs;
      fs.face_id  = f;
      fs.elem1    = e1;
      fs.elem2    = e2;
      fs.is_fault = (is_fault_face[f] != 0);
      in.faces.push_back(fs);
   }
   return in;
}

// (LTS Phase 4 Step 0) Build a symmetric, self-loop-free CSR element-adjacency
// (for METIS) from a mesh's element-to-element table.
void BuildLtsCsrFromMesh(Mesh &mesh, std::vector<int> &xadj,
                         std::vector<int> &adjncy)
{
   const Table &e2e = mesh.ElementToElementTable();
   const int n = e2e.Size();
   const int *I = e2e.GetI();
   const int *J = e2e.GetJ();
   xadj.assign(static_cast<std::size_t>(n) + 1, 0);
   adjncy.clear();
   adjncy.reserve(static_cast<std::size_t>(I[n]));
   for (int e = 0; e < n; ++e)
   {
      for (int k = I[e]; k < I[e + 1]; ++k)
      {
         const int nb = J[k];
         if (nb != e && nb >= 0) { adjncy.push_back(nb); }   // drop self-loops
      }
      xadj[static_cast<std::size_t>(e) + 1] = static_cast<int>(adjncy.size());
   }
}

// (LTS Phase 4 Step 0) Populate an LtsClusteringOptions from the [numerics]
// config (single-sourced; used by every clustering call site).
LtsClusteringOptions MakeLtsClusteringOptions(int nc_cap, double merge_loss_tol,
                                              const std::string &wiggle)
{
   LtsClusteringOptions opt;
   opt.nc_cap         = nc_cap;
   opt.merge_loss_tol = merge_loss_tol;
   if (wiggle == "scan") { opt.wiggle_scan = true; }
   else if (wiggle == "off") { opt.wiggle_scan = false; opt.lambda_fixed = 1.0; }
   else
   {
      opt.wiggle_scan  = false;
      opt.lambda_fixed = std::strtod(wiggle.c_str(), nullptr);
   }
   return opt;
}

// Map TOML mixed_flux string -> WaveOperator enum.
MixedFluxMode ParseMixedFlux(const std::string &s)
{
   if (s == "none")            { return MixedFluxMode::None; }
   if (s == "adjacent")        { return MixedFluxMode::Adjacent; }
   if (s == "all_continuous")  { return MixedFluxMode::AllContinuous; }
   MFEM_ABORT("spatial_dyn_driver: --mixed-flux: unknown value '" << s
              << "'.  Accepted: none | adjacent | all_continuous.");
}

ParaViewOutput<ParMesh>::VolumeOutputMode ParseVolumeMode(
   const std::string &s, bool &enabled)
{
   if (s == "off") { enabled = false; return ParaViewOutput<ParMesh>::DefaultVolumeOutputMode(); }
   enabled = true;
   if (s == "vtu")  { return ParaViewOutput<ParMesh>::VolumeOutputMode::Vtu; }
   if (s == "hdf5") { return ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5; }
   MFEM_ABORT("spatial_dyn_driver: paraview_volume: unknown value '"
              << s << "'.  Accepted: hdf5 | vtu | off.");
}

ParaViewOutput<ParMesh>::FaultOutputMode ParseFaultMode(
   const std::string &s, bool &enabled)
{
   if (s == "off") { enabled = false; return ParaViewOutput<ParMesh>::FaultOutputMode::Vtu; }
   enabled = true;
   if (s == "vtu")  { return ParaViewOutput<ParMesh>::FaultOutputMode::Vtu; }
   if (s == "hdf5") { return ParaViewOutput<ParMesh>::FaultOutputMode::Hdf5; }
   MFEM_ABORT("spatial_dyn_driver: paraview_fault: unknown value '"
              << s << "'.  Accepted: hdf5 | vtu | off.");
}

// Phase 3: free-surface slice mode (mirror ParseVolumeMode/ParseFaultMode).
seas::FreeSurfaceOutput::Mode ParseFreeSurfaceMode(
   const std::string &s, bool &enabled)
{
   if (s == "off") { enabled = false; return seas::FreeSurfaceOutput::Mode::Vtu; }
   enabled = true;
   if (s == "vtu")  { return seas::FreeSurfaceOutput::Mode::Vtu; }
   if (s == "hdf5") { return seas::FreeSurfaceOutput::Mode::Hdf5; }
   MFEM_ABORT("spatial_dyn_driver: paraview_free_surface: unknown value '"
              << s << "'.  Accepted: hdf5 | vtu | off.");
}

// Compute the count of fault QPs per face for a given fault-face quadrature
// exactness degree.  Peek at the first interior or shared fault face on this
// rank, then MPI_Allreduce(MAX) so every rank agrees.
//
// Phase 4: callers pass `wave.FaultFaceQuadDegree()` (= 2*order with
// over-integration off, so this is byte-identical to the prior `2*order`
// probe; = 2*(order+k) with `--fault-overint k`).
int ProbeNbfPerFace(ParMesh &pmesh, int fault_quad_degree,
                    const Array<int> &fault_int_faces,
                    const Array<int> &fault_shr_faces,
                    MPI_Comm comm)
{
   int nbf = 0;
   if (fault_int_faces.Size() > 0)
   {
      FaceElementTransformations *ftr =
         pmesh.GetInteriorFaceTransformations(fault_int_faces[0]);
      MFEM_VERIFY(ftr, "spatial_dyn: fault interior face has null FTR");
      nbf = IntRules.Get(ftr->GetGeometryType(), fault_quad_degree).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   else if (fault_shr_faces.Size() > 0)
   {
      FaceElementTransformations *ftr =
         pmesh.GetSharedFaceTransformations(fault_shr_faces[0]);
      MFEM_VERIFY(ftr, "spatial_dyn: fault shared face has null FTR");
      nbf = IntRules.Get(ftr->GetGeometryType(), fault_quad_degree).GetNPoints();
   }
   {
      int local_nbf = nbf;
      MPI_Allreduce(&local_nbf, &nbf, 1, MPI_INT, MPI_MAX, comm);
   }
#else
   (void)comm;
#endif
   return nbf;
}

// Walk the wave operator's fault face lists and produce the per-DOF
// arrays the new FaultGeometry ctor / spatial_setup.hpp expect.
//
// Outputs (all owned-fault, no double-counting; interior first, shared
// appended in fault_shared_faces order):
//   - phys_coords[i] = (x, y, z) at fault QP i
//   - dof_basis(9, i) = [n; t1=dip; t2=strike] at QP i
//   - dof_to_elem[i] = bulk Elem1 owning QP i
//   - dof_to_attr[i] = mesh boundary attribute at QP i (fault_attr for
//     interior fault faces; shared faces report the bdr attr if any,
//     else fault_attr)
//   - dof_ips[i] = reference-element IntegrationPoint at QP i
//     (FaultGeometry::fault_dof_ip cache; consumed by
//     spatial_setup.hpp's IP-aware overload).
void BuildPerDOFFaultTables(ParMesh &pmesh,
                            int fault_quad_degree,
                            int fault_attr,
                            const FaultBasis &fbasis,
                            const Array<int> &fault_int_faces,
                            const Array<int> &fault_shr_faces,
                            int nbf_per_face,
                            std::vector<Vector> &phys_coords,
                            Vector &dof_coords_3d,
                            DenseMatrix &dof_basis,
                            Array<int> &dof_to_elem,
                            Array<int> &dof_to_attr,
                            std::vector<IntegrationPoint> &dof_ips)
{
   const int num_int  = fault_int_faces.Size() * nbf_per_face;
   const int num_shr  = fault_shr_faces.Size() * nbf_per_face;
   const int N        = num_int + num_shr;

   phys_coords.clear();
   phys_coords.reserve(N);
   dof_coords_3d.SetSize(3 * N);
   dof_basis.SetSize(9, N);
   dof_to_elem.SetSize(N);
   dof_to_attr.SetSize(N);
   dof_ips.assign(N, IntegrationPoint());

   // R-005: per-QP basis for curvilinear-face fidelity.  When the caller
   // populated `bdata.qp_data` (via `FaultBasis::ComputeQPBasis` /
   // `ComputeQPBasisShared`) the per-QP frame at QP `q` is used; otherwise
   // we fall back to the centroid frame (planar-face accuracy, byte-
   // compatible with the pre-R-005 inline walk).
   auto qp_normal = [](const FaultBasisData &bdata, int q, int d) -> real_t
   {
      return (static_cast<int>(bdata.qp_data.size()) > q)
                ? bdata.qp_data[q].normal[d]   : bdata.normal[d];
   };
   auto qp_t1 = [](const FaultBasisData &bdata, int q, int d) -> real_t
   {
      return (static_cast<int>(bdata.qp_data.size()) > q)
                ? bdata.qp_data[q].tangent1[d] : bdata.tangent1[d];
   };
   auto qp_t2 = [](const FaultBasisData &bdata, int q, int d) -> real_t
   {
      return (static_cast<int>(bdata.qp_data.size()) > q)
                ? bdata.qp_data[q].tangent2[d] : bdata.tangent2[d];
   };

   auto write_dof = [&](int dof_idx, FaceElementTransformations *ftr,
                        const IntegrationPoint &ip,
                        const FaultBasisData &bdata,
                        int q,
                        int attr)
   {
      ftr->SetAllIntPoints(&ip);
      Vector phys(3);
      ftr->Face->Transform(ip, phys);
      phys_coords.push_back(phys);

      for (int d = 0; d < 3; ++d)
      {
         dof_coords_3d(3 * dof_idx + d) = phys(d);
         dof_basis(0 + d, dof_idx) = qp_normal(bdata, q, d);
         dof_basis(3 + d, dof_idx) = qp_t1    (bdata, q, d);
         dof_basis(6 + d, dof_idx) = qp_t2    (bdata, q, d);
      }
      dof_to_elem[dof_idx] = ftr->Elem1No;
      dof_to_attr[dof_idx] = attr;
      dof_ips[dof_idx]     = ftr->GetElement1IntPoint();
   };

   int dof_idx = 0;

   for (int fi = 0; fi < fault_int_faces.Size(); ++fi)
   {
      const int face = fault_int_faces[fi];
      FaceElementTransformations *ftr =
         pmesh.GetInteriorFaceTransformations(face);
      MFEM_VERIFY(ftr, "spatial_dyn: interior fault face has null FTR");
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), fault_quad_degree);
      MFEM_VERIFY(ir.GetNPoints() == nbf_per_face,
                  "spatial_dyn: interior fault face has " << ir.GetNPoints()
                  << " QPs, expected " << nbf_per_face);
      const FaultBasisData &bdata = fbasis.GetBasis(fi);
      for (int q = 0; q < nbf_per_face; ++q)
      {
         write_dof(dof_idx++, ftr, ir.IntPoint(q), bdata, q, fault_attr);
      }
   }
#ifdef MFEM_USE_MPI
   for (int si = 0; si < fault_shr_faces.Size(); ++si)
   {
      const int sf = fault_shr_faces[si];
      FaceElementTransformations *ftr =
         pmesh.GetSharedFaceTransformations(sf);
      MFEM_VERIFY(ftr, "spatial_dyn: shared fault face has null FTR");
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), fault_quad_degree);
      MFEM_VERIFY(ir.GetNPoints() == nbf_per_face,
                  "spatial_dyn: shared fault face has " << ir.GetNPoints()
                  << " QPs, expected " << nbf_per_face);
      const FaultBasisData &bdata =
         fbasis.GetBasis(fault_int_faces.Size() + si);
      for (int q = 0; q < nbf_per_face; ++q)
      {
         write_dof(dof_idx++, ftr, ir.IntPoint(q), bdata, q, fault_attr);
      }
   }
#else
   (void)fault_shr_faces;
#endif
   MFEM_VERIFY(dof_idx == N,
               "spatial_dyn: per-DOF table walk wrote " << dof_idx
               << " entries, expected " << N);
}

// TPV205-style ADER macro-step driver: predictor in the bulk, LSW
// closed-form per sub-step at fault QPs, corrector via wave.AdvanceADER
// with the side-channel I_imp.  Mirrors AdvanceADERWithSubStep in
// drivers/tpv205_driver.cpp (1:1 except no SEAS_DIAG hooks).
//
// Phase N: the trailing `nuc_callback` arg is forwarded to the
// per-sub-step `IFrictionIterator::Advance` callback; the callback fires
// ONCE per ADER sub-step BEFORE the per-QP friction pipeline.  Pass
// `[](real_t, real_t){}` to opt out (no nucleation perturbation).
//
// Phase 2 (R-010): the iterator is the runtime-dispatch IFrictionIterator
// strategy (LSW or rate-and-state adapter), not the concrete
// Tpv205SubStepIterator; the body calls only interface methods.
void AdvanceADERWithSubStep_Spatial(
   WaveOperator<ParMesh> &wave,
   IFrictionIterator &iterator,
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   const Vector &Q,
   real_t dt_step,
   int ader_order,
   real_t t_step_start,
   Vector &Q_new,
   const std::function<void(real_t, real_t)> &nuc_callback,
   bool use_shared_ck)
{
   MFEM_PERF_SCOPE("seas::spatial_dyn::AdvanceADERWithSubStep");
   MFEM_VERIFY(dt_step > 0.0,
               "AdvanceADERWithSubStep_Spatial: dt_step must be > 0, got "
               << dt_step);
   MFEM_VERIFY(ader_order >= 2 && ader_order <= 4,
               "AdvanceADERWithSubStep_Spatial: ader_order must be in "
               "{2,3,4}, got " << ader_order);

   const std::vector<real_t> configured_deltaT  = iterator.GetDeltaT();
   const std::vector<real_t> configured_weights = iterator.GetTimeWeights();
   const int O = static_cast<int>(configured_deltaT.size());
   MFEM_VERIFY(O >= 1,
               "AdvanceADERWithSubStep_Spatial: iterator deltaT empty; "
               "SetSubSteps must be called first.");
   const real_t configured_sum =
      std::accumulate(configured_deltaT.begin(), configured_deltaT.end(),
                      static_cast<real_t>(0));
   MFEM_VERIFY(configured_sum > 0.0,
               "AdvanceADERWithSubStep_Spatial: Σ deltaT = "
               << configured_sum << " ≤ 0");

   const real_t dt_scale = dt_step / configured_sum;
   std::vector<real_t> deltaT_scaled(O);
   for (int o = 0; o < O; ++o)
   {
      deltaT_scaled[o] = configured_deltaT[o] * dt_scale;
   }
   iterator.SetSubSteps(deltaT_scaled, configured_weights);
   const std::vector<real_t> &deltaT = iterator.GetDeltaT();

   // Sub-step MIDPOINT nodes on [0, dt_step].
   std::vector<real_t> tau_nodes(O);
   real_t acc = 0.0;
   for (int o = 0; o < O; ++o)
   {
      tau_nodes[o] = acc + 0.5 * deltaT[o];
      acc += deltaT[o];
   }

   std::vector<Vector> Q_per_node;
   // Lever 2: when enabled, ONE CK recursion produces both Q_per_node and the
   // time integral shared_I (passed to AdvanceADER below), removing the second
   // recursion.  Runs on ALL ranks (use_shared_ck is the same CLI value
   // everywhere) so Q_per_node still feeds the per-substep ExchangeFaceNbrData
   // collective in EvaluateBulkAtFaultQPsCanonical — REVIEW R-001.  Default
   // (use_shared_ck=false) is the original two-recursion path, byte-identical.
   Vector shared_I;
   if (use_shared_ck)
   {
      wave.ComputeADERSubStepStatesAndIntegral(Q, dt_step, ader_order, tau_nodes,
                                               Q_per_node, shared_I);
   }
   else
   {
      wave.ComputeADERSubStepStates(Q, dt_step, ader_order, tau_nodes,
                                    Q_per_node);
   }
   MFEM_VERIFY(static_cast<int>(Q_per_node.size()) == O,
               "AdvanceADERWithSubStep_Spatial: ComputeADERSubStepStates "
               "returned " << Q_per_node.size() << " nodes, expected " << O);

   const int n_total_fault_qps = wave.GetNumTotalFaultQPs();
   std::vector<std::vector<real_t>> Q_pointwise_plus(O), Q_pointwise_minus(O);
   for (int o = 0; o < O; ++o)
   {
      wave.EvaluateBulkAtFaultQPsCanonical(Q_per_node[o],
                                           Q_pointwise_plus[o],
                                           Q_pointwise_minus[o]);
   }

   const size_t n_words =
      static_cast<size_t>(NUM_STATE)
      * static_cast<size_t>(n_total_fault_qps);
   std::vector<real_t> I_imp_plus_flat(n_words, 0.0);
   std::vector<real_t> I_imp_minus_flat(n_words, 0.0);

   if (n_total_fault_qps > 0)
   {
      MFEM_PERF_BEGIN("seas::spatial_dyn::friction_substep");
      iterator.Advance(dof_data, fault_coords,
                       Q_pointwise_plus,
                       Q_pointwise_minus,
                       dt_step, t_step_start,
                       I_imp_plus_flat.data(),
                       I_imp_minus_flat.data(),
                       nuc_callback);
      MFEM_PERF_END("seas::spatial_dyn::friction_substep");
   }

   struct ImposedGuard
   {
      WaveOperator<ParMesh> &w_;
      explicit ImposedGuard(WaveOperator<ParMesh> &w) : w_(w) {}
      ~ImposedGuard() { w_.ResetSubStepFaultImposedStates(); }
   };
   wave.SetSubStepFaultImposedStates(
      n_total_fault_qps > 0 ? I_imp_plus_flat.data()  : nullptr,
      n_total_fault_qps > 0 ? I_imp_minus_flat.data() : nullptr,
      n_total_fault_qps);
   ImposedGuard guard(wave);

   wave.AdvanceADER(Q, dt_step, ader_order, Q_new,
                    use_shared_ck ? &shared_I : nullptr);
}

// -------------------------------------------------------------------------
// Phase 3 of PLAN_thermal_case2_mixedflux_port_2026-07-08.md.
//
// Open the `[friction.rate_state.sidecar]` HDF5 readers and type-erase them
// into the `RateStateSidecarFields` callables the resolver consumes.  The
// `DataField3D` readers are held by `shared_ptr` CAPTURED BY VALUE inside each
// lambda, so they outlive the resolver and the returned struct owns them
// transitively -- no dangling reference if the caller drops its own handle.
//
// `far_field_clamp` maps to `OOBPolicy::Clamp` (ASAGI nearest-edge hold): the
// SCEC CTM hull covers the SAFS fault but not the far-field absorbing box.
// With `false` an out-of-hull query is a hard abort (the schema-v1 default).
//
// Returns nullptr when the block is absent, which `SpatialFrictionResolver`
// accepts as "no sidecar" -- so the call site needs no branch.
// -------------------------------------------------------------------------
std::shared_ptr<const seas::spatial::RateStateSidecarFields>
MakeRateStateSidecarFields(const seas::spatial::FrictionSidecarSpec &spec
#ifdef MFEM_USE_MPI
                           // PLAN_sidecar_mpi_shared_memory Phase 0: pass a
                           // node-local comm to hold each field grid in ONE
                           // MPI-3 shared window per node; MPI_COMM_NULL =
                           // the classic per-rank replication.
                           , MPI_Comm node_comm = MPI_COMM_NULL
#endif
                           )
{
   using seas::spatial::RateStateSidecarFields;
   if (!spec.enabled) { return nullptr; }

   const seas::OOBPolicy oob = spec.far_field_clamp
                               ? seas::OOBPolicy::Clamp
                               : seas::OOBPolicy::Abort;

   // One reader per field.  Each ctor validates the schema-v1 invariants
   // (schema_version / crs / z_positive / monotone axes / no NaN / bounds)
   // and aborts with a precise message on a bad file, so a missing dataset
   // fails HERE, at setup, not mid-time-loop.
   auto open_reader = [&spec, oob
#ifdef MFEM_USE_MPI
                       , node_comm
#endif
                      ](const std::string &field)
   {
#ifdef MFEM_USE_MPI
      if (node_comm != MPI_COMM_NULL)
      {
         return std::make_shared<seas::DataField3D>(spec.path, field, oob,
                                                    node_comm);
      }
#endif
      return std::make_shared<seas::DataField3D>(spec.path, field, oob);
   };
   // Type-erase a reader into the resolver's callable.  The shared_ptr is
   // CAPTURED BY VALUE, so the reader outlives every handle the caller drops.
   auto bind = [](std::shared_ptr<seas::DataField3D> reader)
   {
      return RateStateSidecarFields::FieldFn(
                [reader](real_t x, real_t y, real_t z)
                { return reader->Evaluate(x, y, z); });
   };

   auto fields = std::make_shared<RateStateSidecarFields>();

   // R-001: record the data hull so `--print-derived` can refuse to sample a
   // point outside it.  Every field of one schema-v1 sidecar shares the same
   // /grid, so the `a` reader's bbox describes them all (StressField3D asserts
   // this invariant for its own six components).
   auto a_reader = open_reader(spec.a_field);
   fields->bbox               = a_reader->BBox();
   fields->clamps_out_of_hull = spec.far_field_clamp;
   fields->a                  = bind(a_reader);

   fields->V_w = bind(open_reader(spec.V_w_field));
   if (!spec.b_field.empty())  { fields->b  = bind(open_reader(spec.b_field)); }
   if (!spec.Dc_field.empty()) { fields->Dc = bind(open_reader(spec.Dc_field)); }
   return fields;
}

// -------------------------------------------------------------------------
// `--print-derived` diagnostic for the friction sidecar (Phase 3 acceptance).
//
// Reports, over ALL fault DOFs on ALL ranks: min / max / mean of the resolved
// `a` and `V_w`, the velocity-weakening DOF count (a - b < 0), and the sampled
// sidecar values at the nucleation centre (the hypocentre for this problem).
// This is the online counterpart of `build_pref_friction_sidecar.py`'s own G2
// (range) and G4 (hypocentre anchor) gates -- if the two disagree, the mesh,
// the sidecar, or the CRS is wrong.
//
// Collective: every rank must call it.  Prints on rank 0 only.
// -------------------------------------------------------------------------
void PrintFrictionSidecarSummary(
   const seas::spatial::RateStatePerDOFParams &rs,
   const seas::spatial::RateStateSidecarFields &fields,
   const seas::spatial::HypocenterSpec &hypo,
   const seas::spatial::NucleationSpec &nuc,
#ifdef MFEM_USE_MPI
   MPI_Comm comm,
#endif
   int rank)
{
   const int N = rs.a.Size();

   real_t a_min = std::numeric_limits<real_t>::infinity();
   real_t a_max = -std::numeric_limits<real_t>::infinity();
   real_t v_min = std::numeric_limits<real_t>::infinity();
   real_t v_max = -std::numeric_limits<real_t>::infinity();
   real_t a_sum = 0.0, v_sum = 0.0;
   long   n_vw  = 0;   // a - b < 0  (velocity-weakening)

   for (int i = 0; i < N; ++i)
   {
      a_min = std::min(a_min, rs.a(i));
      a_max = std::max(a_max, rs.a(i));
      v_min = std::min(v_min, rs.V_w(i));
      v_max = std::max(v_max, rs.V_w(i));
      a_sum += rs.a(i);
      v_sum += rs.V_w(i);
      if (rs.a(i) - rs.b(i) < 0.0) { ++n_vw; }
   }

   long n_tot = N;
#ifdef MFEM_USE_MPI
   auto all_min = [comm](real_t v)
   { real_t o; MPI_Allreduce(&v, &o, 1, MPITypeMap<real_t>::mpi_type,
                             MPI_MIN, comm); return o; };
   auto all_max = [comm](real_t v)
   { real_t o; MPI_Allreduce(&v, &o, 1, MPITypeMap<real_t>::mpi_type,
                             MPI_MAX, comm); return o; };
   auto all_sum_r = [comm](real_t v)
   { real_t o; MPI_Allreduce(&v, &o, 1, MPITypeMap<real_t>::mpi_type,
                             MPI_SUM, comm); return o; };
   auto all_sum_l = [comm](long v)
   { long o; MPI_Allreduce(&v, &o, 1, MPI_LONG, MPI_SUM, comm); return o; };

   a_min = all_min(a_min);  a_max = all_max(a_max);
   v_min = all_min(v_min);  v_max = all_max(v_max);
   a_sum = all_sum_r(a_sum); v_sum = all_sum_r(v_sum);
   n_vw  = all_sum_l(n_vw);  n_tot = all_sum_l(n_tot);
#endif

   if (rank != 0) { return; }

   // A rank with zero fault DOFs leaves min/max at +-inf; a GLOBAL n_tot of 0
   // means the fault carries no DOFs at all, which is already fatal upstream.
   MFEM_VERIFY(n_tot > 0,
               "PrintFrictionSidecarSummary: zero fault DOFs globally.");

   std::cout << "\n[print-derived] friction sidecar (a, V_w)\n"
             << "  fault DOFs (global) : " << n_tot << "\n"
             << "  a      min/mean/max : " << a_min << " / "
             << (a_sum / static_cast<real_t>(n_tot)) << " / " << a_max << "\n"
             << "  V_w    min/mean/max : " << v_min << " / "
             << (v_sum / static_cast<real_t>(n_tot)) << " / " << v_max
             << " m/s\n"
             << "  velocity-weakening  : " << n_vw << " / " << n_tot
             << " DOFs with a - b < 0  ("
             << (100.0 * static_cast<real_t>(n_vw)
                 / static_cast<real_t>(n_tot)) << " %)\n";

   // Sample the sidecar at the nucleation centre when nucleation is enabled
   // (that IS the hypocentre for the SAFS decks); otherwise fall back to the
   // [hypocenter] block.  Both are plain coordinate lookups into the same
   // field the resolver used, so they cross-check the offline G4 gate.
   real_t hx = hypo.x_m, hy = hypo.y_m, hz = hypo.z_m;
   const char *origin = "[hypocenter]";
   bool have_point = (hx != 0.0 || hy != 0.0 || hz != 0.0);
   if (nuc.enabled)
   {
      // R-007: `have_point` is cleared here and re-set INSIDE each case.  No
      // `default:` label, so `-Wswitch` still flags a newly-added NucleationKind
      // at compile time; and if one slips through anyway, the diagnostic
      // degrades to "skipped" rather than silently sampling the [hypocenter]
      // default of (0, 0, 0).
      have_point = false;
      switch (nuc.kind)
      {
         case seas::spatial::NucleationKind::GradualOverstress:
            hx = nuc.gradual_overstress.center_x_m;
            hy = nuc.gradual_overstress.center_y_m;
            hz = nuc.gradual_overstress.center_z_m;
            origin = "[nucleation.gradual_overstress] centre";
            have_point = true;
            break;
         case seas::spatial::NucleationKind::GradualOverstressCompactCircular:
            hx = nuc.compact_circular.center_x_m;
            hy = nuc.compact_circular.center_y_m;
            hz = nuc.compact_circular.center_z_m;
            origin = "[nucleation.gradual_overstress_compact_circular] centre";
            have_point = true;
            break;
         case seas::spatial::NucleationKind::InstantaneousOverstressCircular:
            hx = nuc.instantaneous_circular.center_x_m;
            hy = nuc.instantaneous_circular.center_y_m;
            hz = nuc.instantaneous_circular.center_z_m;
            origin = "[nucleation.instantaneous_overstress_circular] centre";
            have_point = true;
            break;
      }
   }

   // Skip the point sample when neither a nucleation centre nor an explicit
   // [hypocenter] was given: the (0,0,0) default is outside any UTM 11 N hull.
   if (!have_point)
   {
      std::cout << "  sidecar point sample : skipped (no [nucleation] centre "
                   "and no [hypocenter])\n" << std::flush;
      return;
   }

   // R-001: NEVER evaluate outside the data hull.
   //   * Under OOBPolicy::Clamp the read returns a nearest-edge value that is
   //     indistinguishable from a real sample.  For the SCEC CTM CASE2 field the
   //     edge values are EXACTLY the hypocentre anchors (a = 0.015, V_w = 0.05)
   //     this diagnostic exists to verify, so an out-of-hull point would print a
   //     perfect "pass" — the gate could never fail.
   //   * Under OOBPolicy::Abort the read aborts, killing the run from inside a
   //     diagnostic.
   // Report the miss instead; it is far more informative than either.
   if (!fields.InHull(hx, hy, hz))
   {
      std::cout << "  sidecar at " << origin << " (" << hx << ", " << hy << ", "
                << hz << "): OUTSIDE the sidecar hull "
                << "x[" << fields.bbox[0] << ", " << fields.bbox[1] << "] "
                << "y[" << fields.bbox[2] << ", " << fields.bbox[3] << "] "
                << "z[" << fields.bbox[4] << ", " << fields.bbox[5] << "]\n"
                << "      NOT SAMPLED — an out-of-hull read would "
                << (fields.clamps_out_of_hull
                    ? "silently clamp to an edge value (a false PASS)."
                    : "abort the run.")
                << "\n" << std::flush;
      return;
   }

   std::cout << "  sidecar at " << origin << " (" << hx << ", " << hy << ", "
             << hz << ")  [in hull]:\n"
             << "      a   = " << fields.a(hx, hy, hz) << "\n"
             << "      V_w = " << fields.V_w(hx, hy, hz) << " m/s\n";
   if (fields.b)  { std::cout << "      b   = " << fields.b(hx, hy, hz) << "\n"; }
   if (fields.Dc) { std::cout << "      Dc  = " << fields.Dc(hx, hy, hz) << " m\n"; }
   std::cout << std::flush;
}

}  // namespace

// =========================================================================
// main
// =========================================================================

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);
#else
   int rank = 0, nprocs = 1;
#endif

#ifdef SEAS_DIAG_FAULT_FLUX
   mfem::seas::g_seas_my_rank = rank;
#endif

   // -----------------------------------------------------------------
   // 1.  CLI parse — config-driven; CLI overrides are merged after
   //     TOML load (later wins).
   // -----------------------------------------------------------------
   const std::string config_path =
      GetStringArg(argc, argv, "--config", "");
   if (config_path.empty())
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: --config PATH.toml is required.\n"
                   << "  Usage: seas_spatial_dyn_driver --config "
                   << "<file>.toml [OPTIONS]\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 2;
   }

   // (LTS Phase 0) --lts-report: cluster the mesh with our own per-element dt
   // and print the RAW (SeisSol-comparable) + PRODUCTION histograms and the
   // predicted element-update speedups, then exit.  IMPLIES --dry-run (report
   // and quit before any stepping).
   const bool lts_report      = HasFlag(argc, argv, "--lts-report");
   // LTS Track-A A1: merge the per-correct seam exchanges into ONE collective
   // per buffer per tick (125 -> 64 rounds/sync at Nc=6).  Opt-in; bitwise
   // identical by construction (document/comm_dev/DESIGN_a1_per_tick_merge_2026-07-24.md).
   // Bulk (fault-free) LTS path only -- the fault-interleave stepper is not
   // wired yet and the flag is refused there rather than silently ignored.
   const bool lts_merge_exch = HasFlag(argc, argv, "--lts-merge-exchanges");
   const bool dry_run         = HasFlag(argc, argv, "--dry-run") || lts_report;
   const bool verify_dispatch = HasFlag(argc, argv, "--verify-dispatch");
   const bool no_sidecar_material =
      HasFlag(argc, argv, "--no-sidecar-material");
   const bool print_derived = HasFlag(argc, argv, "--print-derived");
   // PLAN_sidecar_mpi_shared_memory_2026-07-17.md Phase 0: opt-in MPI-3
   // shared-memory sidecar windows — ONE physical copy of each gridded
   // sidecar (material / stress / friction) per NODE instead of one per
   // rank.  Default OFF = the classic per-rank replication, byte-identical
   // to the pre-flag behaviour.
   const bool sidecar_shared_mem = HasFlag(argc, argv, "--sidecar-shared-mem");

   const std::string cli_mesh        = GetStringArg(argc, argv, "--mesh", "");
   const std::string cli_vel_model   =
      GetStringArg(argc, argv, "--velocity-model", "");
   const std::string cli_vel_path    =
      GetStringArg(argc, argv, "--override-velocity-path", "");
   const std::string cli_stress_kind =
      GetStringArg(argc, argv, "--stress-kind", "");
   const std::string cli_stress_sidecar =
      GetStringArg(argc, argv, "--stress-sidecar", "");

   const real_t cli_tfinal     = GetRealArg(argc, argv, "--tfinal", -1.0);
   const real_t cli_cfl        = GetRealArg(argc, argv, "--cfl", -1.0);
   // Lever A (2026-07-01): CLI override for the extra DG safety margin
   // [numerics].cfl_dg_safety (default 3.0 = today; 1.0 = SeisSol-equivalent,
   // dt×3).  Sweep knob — makes a cfl_dg_safety A/B one flag, no config dupes.
   // The Courant number `cfl` is a SEPARATE knob and stays ≤ 1.
   const real_t cli_cfl_dg_safety = GetRealArg(argc, argv, "--cfl-dg-safety", -1.0);
   const int    cli_ader_order = GetIntArg(argc, argv, "--ader-order", -1);
   // Phase 4: --fault-overint K — fault-flux over-integration factor (0 = off,
   // byte-exact).  Applied via WaveOperator::SetFaultOverint after the operator
   // is built (below).
   const int    cli_fault_overint = GetIntArg(argc, argv, "--fault-overint", 0);
   // Phase 3: --fault-resample — secular resample of the rate-state Δψ increment
   // (degree-N L2 projection per fault face).  Default off ⇒ byte-exact.  A
   // no-op without over-integration at a unisolvent rule (R = I); rate-state
   // only (LSW resample is not implemented — see the guard below).
   const bool   cli_fault_resample = HasFlag(argc, argv, "--fault-resample");
   // Phase 14: --time-integrator ader|rk4|rk45 (mirror --ader-order; empty ⇒
   // keep the TOML/default).  Validated in the config-override block below.
   const std::string cli_time_integrator =
      GetStringArg(argc, argv, "--time-integrator", "");
   const std::string cli_mixed_flux =
      GetStringArg(argc, argv, "--mixed-flux", "");
   // (Unified bi-material plan, Part A) central-flux corridor contrast guard.
   // Sentinel -1e30 = "not given on the CLI" (keep TOML/default); any value >= that
   // sentinel was given and overrides the config (CLI wins).  A negative given value
   // explicitly DISABLES the guard (so the CLI can override an enabling TOML).
   constexpr real_t kContrastTolUnset = -1.0e30;
   const real_t cli_contrast_tol =
      GetRealArg(argc, argv, "--mixed-flux-contrast-tol", kContrastTolUnset);
   const bool   cli_pml        = HasFlag(argc, argv, "--pml");
   // Phase 12.2: PML overrides (CLI wins over [numerics]); sentinels mean
   // "not given on the CLI" so the TOML / struct default is kept.
   const real_t cli_pml_thickness = GetRealArg(argc, argv, "--pml-thickness", -1.0);
   const real_t cli_pml_target_R  = GetRealArg(argc, argv, "--pml-target-R",  -1.0);
   const int    cli_pml_cells     = GetIntArg (argc, argv, "--pml-cells",     -1);
   const int    cli_pml_damp_bot  = GetIntArg (argc, argv, "--pml-damp-bottom", -1);
   const int    cli_pml_damp_top  = GetIntArg (argc, argv, "--pml-damp-top",    -1);

   const std::string cli_output_dir =
      GetStringArg(argc, argv, "--output-dir", "");
   const std::string cli_pv_volume = GetStringArg(argc, argv, "--paraview-volume", "");
   const std::string cli_pv_bulk   = GetStringArg(argc, argv, "--paraview-bulk", "");
   const std::string cli_pv_fault  = GetStringArg(argc, argv, "--paraview-fault", "");
   const real_t cli_pv_vol_dt   = GetRealArg(argc, argv, "--paraview-volume-dt", -1.0);
   const real_t cli_pv_bulk_dt  = GetRealArg(argc, argv, "--paraview-bulk-dt", -1.0);
   const real_t cli_pv_fault_dt = GetRealArg(argc, argv, "--paraview-fault-dt", -1.0);
   const real_t cli_pv_vol_zfp  = GetRealArg(argc, argv, "--paraview-volume-zfp-tol", -1.0);
   const real_t cli_pv_bulk_zfp = GetRealArg(argc, argv, "--paraview-bulk-zfp-tol", -1.0);
   const real_t cli_pv_fault_zfp= GetRealArg(argc, argv, "--paraview-fault-zfp-tol", -1.0);
   const int    cli_pv_max_snap = GetIntArg(argc, argv, "--paraview-max-snapshots", -1);

   // Parity Phase 1: extended ParaView CLI flag surface.
   const bool   cli_pv                   = HasFlag(argc, argv, "--paraview");
   const real_t cli_pv_dt                = GetRealArg(argc, argv, "--paraview-dt", -1.0);
   const int    cli_pv_every             = GetIntArg (argc, argv, "--paraview-every", 0);
   const bool   cli_pv_force_fault_vtu   = HasFlag(argc, argv, "--paraview-fault-vtu");
   const bool   cli_pv_force_fault_hdf5  = HasFlag(argc, argv, "--paraview-fault-hdf5");
   const bool   cli_pv_legacy_ascii      = HasFlag(argc, argv, "--paraview-fault-legacy-ascii");
   const int    cli_pv_fault_deflate     = GetIntArg (argc, argv, "--paraview-fault-deflate-level", -1);
   const int    cli_pv_bulk_deflate      = GetIntArg (argc, argv, "--paraview-bulk-deflate-level",  -1);
   const bool   cli_pv_force_vol_vtu     = HasFlag(argc, argv, "--paraview-volume-vtu");
   const bool   cli_pv_force_vol_hdf5    = HasFlag(argc, argv, "--paraview-volume-hdf5");
   const bool   cli_pv_force_bulk_vtu    = HasFlag(argc, argv, "--paraview-bulk-vtu");
   const bool   cli_pv_force_bulk_hdf5   = HasFlag(argc, argv, "--paraview-bulk-hdf5");
   const int    cli_pv_volume_deflate    = GetIntArg (argc, argv, "--paraview-volume-deflate-level", -1);
   const real_t cli_pv_coseismic_dt      = GetRealArg(argc, argv, "--paraview-coseismic-dt",    -1.0);
   const real_t cli_pv_nucleation_dt     = GetRealArg(argc, argv, "--paraview-nucleation-dt",   -1.0);
   const real_t cli_pv_interseismic_dt   = GetRealArg(argc, argv, "--paraview-interseismic-dt", -1.0);
   // Phase 3: free-surface slice CLI overrides.
   const std::string cli_pv_free_surface = GetStringArg(argc, argv, "--paraview-free-surface", "");
   const real_t cli_pv_free_surface_dt   = GetRealArg(argc, argv, "--paraview-free-surface-dt", -1.0);

   const std::string restart_prefix =
      GetStringArg(argc, argv, "--restart", "");
   const int cli_checkpoint_every =
      GetIntArg(argc, argv, "--checkpoint-every", -1);

   // -----------------------------------------------------------------
   // 2.  TOML load — owns all defaults; CLI then overrides.
   // -----------------------------------------------------------------
   spatial::SpatialFrictionConfig cfg;
   try
   {
      cfg = spatial::LoadSpatialFrictionConfig(config_path);
   }
   catch (...)
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: failed to load TOML config '" << config_path
                   << "'.\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 2;
   }

   // Merge CLI overrides into cfg (later wins).
   if (!cli_mesh.empty())            { cfg.mesh.path = cli_mesh; }
   if (!cli_vel_model.empty())
   {
      if (cli_vel_model == "cvmh")
      { cfg.velocity.model = spatial::VelocityModel::CVMH; }
      else if (cli_vel_model == "cvm_s4.26.m01")
      { cfg.velocity.model = spatial::VelocityModel::CVMS_4_26_M01; }
      else if (cli_vel_model == "multiscale_statewise")
      { cfg.velocity.model = spatial::VelocityModel::MultiscaleStatewise; }
      else
      {
         MFEM_ABORT("--velocity-model: unknown value '" << cli_vel_model
                    << "'.  Accepted: cvmh | cvm_s4.26.m01 | "
                    << "multiscale_statewise.");
      }
   }
   if (!cli_vel_path.empty())        { cfg.velocity.override_path = cli_vel_path; }
   if (!cli_stress_kind.empty())
   {
      if (cli_stress_kind == "constant_tensor")
      { cfg.stress.kind = spatial::StressSourceKind::ConstantTensor; }
      else if (cli_stress_kind == "sidecar_hdf5")
      { cfg.stress.kind = spatial::StressSourceKind::SidecarHDF5; }
      else
      {
         MFEM_ABORT("--stress-kind: unknown value '" << cli_stress_kind
                    << "'.  Accepted: constant_tensor | sidecar_hdf5.");
      }
   }
   if (!cli_stress_sidecar.empty())  { cfg.stress.sidecar_path = cli_stress_sidecar; }
   if (cli_tfinal > 0.0)             { cfg.time.tfinal = cli_tfinal; }
   if (cli_cfl > 0.0)                { cfg.numerics.cfl = cli_cfl; }
   if (cli_cfl_dg_safety > 0.0)
   {
      // Re-validate here: the parser's >=1.0 guard ran before this override.
      MFEM_VERIFY(cli_cfl_dg_safety >= 1.0,
                  "--cfl-dg-safety must be >= 1.0 (1.0 = SeisSol's order-factor "
                  "floor 1/(2N+1); below it the explicit DG step is unstable); "
                  "got " << cli_cfl_dg_safety);
      cfg.numerics.cfl_dg_safety = cli_cfl_dg_safety;
   }
   if (cli_ader_order > 0)           { cfg.numerics.ader_order = cli_ader_order; }
   if (!cli_time_integrator.empty())
   {
      if      (cli_time_integrator == "ader")
      { cfg.numerics.time_integrator = spatial::TimeIntegratorKind::ADER; }
      else if (cli_time_integrator == "rk4")
      { cfg.numerics.time_integrator = spatial::TimeIntegratorKind::RK4; }
      else if (cli_time_integrator == "rk45")
      { cfg.numerics.time_integrator = spatial::TimeIntegratorKind::RK45; }
      else
      {
         MFEM_ABORT("--time-integrator: unknown value '" << cli_time_integrator
                    << "'.  Accepted: ader | rk4 | rk45.");
      }
   }
   if (!cli_mixed_flux.empty())      { cfg.numerics.mixed_flux = cli_mixed_flux; }
   if (cli_contrast_tol > kContrastTolUnset)
   {
      cfg.numerics.mixed_flux_contrast_tol = cli_contrast_tol;
   }
   if (cli_pml)                      { cfg.numerics.use_pml = true; }
   if (cli_pml_thickness > 0.0)      { cfg.numerics.pml_thickness_m = cli_pml_thickness; }
   if (cli_pml_target_R  > 0.0)      { cfg.numerics.pml_target_R    = cli_pml_target_R; }
   if (cli_pml_cells     > 0)        { cfg.numerics.pml_cells       = cli_pml_cells; }
   if (cli_pml_damp_bot  >= 0)       { cfg.numerics.pml_damp_bottom = (cli_pml_damp_bot != 0); }
   if (cli_pml_damp_top  >= 0)       { cfg.numerics.pml_damp_top    = (cli_pml_damp_top != 0); }
   if (!cli_output_dir.empty())      { cfg.output.output_dir = cli_output_dir; }
   if (!cli_pv_volume.empty())       { cfg.output.paraview_volume = cli_pv_volume; }
   if (!cli_pv_bulk.empty())         { cfg.output.paraview_bulk   = cli_pv_bulk; }
   if (!cli_pv_fault.empty())        { cfg.output.paraview_fault  = cli_pv_fault; }
   if (cli_pv_vol_dt    > 0.0)       { cfg.output.paraview_volume_dt    = cli_pv_vol_dt; }
   if (cli_pv_bulk_dt   > 0.0)       { cfg.output.paraview_bulk_dt      = cli_pv_bulk_dt; }
   if (cli_pv_vol_zfp   > 0.0)       { cfg.output.paraview_volume_zfp_tol = cli_pv_vol_zfp; }
   if (cli_pv_bulk_zfp  > 0.0)       { cfg.output.paraview_bulk_zfp_tol  = cli_pv_bulk_zfp; }
   if (cli_pv_fault_zfp > 0.0)       { cfg.output.paraview_fault_zfp_tol = cli_pv_fault_zfp; }
   if (cli_pv_max_snap  > 0)         { cfg.output.max_snapshots          = cli_pv_max_snap; }
   if (cli_checkpoint_every > 0)     { cfg.output.checkpoint_every_steps = cli_checkpoint_every; }

   // Parity Phase 1: extended ParaView CLI overrides.  Order matters:
   // apply the umbrella `--paraview-dt` FIRST so that the more-specific
   // `--paraview-fault-dt` below overrides it when both are given (R-007
   // fix: specific flag wins over umbrella).
   if (cli_pv)                            { cfg.output.paraview_enabled = true; }
   if (cli_pv_dt              > 0.0)      { cfg.output.paraview_fault_dt    = cli_pv_dt;
                                            cfg.output.paraview_enabled    = true; }
   if (cli_pv_fault_dt  > 0.0)            { cfg.output.paraview_fault_dt   = cli_pv_fault_dt; }
   if (cli_pv_every           > 0)        { cfg.output.paraview_every_steps = cli_pv_every; }
   if (cli_pv_fault_deflate  >= 0)        { cfg.output.paraview_fault_deflate_level   = cli_pv_fault_deflate; }
   if (cli_pv_bulk_deflate   >= 0)        { cfg.output.paraview_bulk_deflate_level    = cli_pv_bulk_deflate; }
   if (cli_pv_volume_deflate >= 0)        { cfg.output.paraview_volume_deflate_level  = cli_pv_volume_deflate; }
   if (cli_pv_coseismic_dt    > 0.0)      { cfg.output.paraview_coseismic_dt    = cli_pv_coseismic_dt; }
   if (cli_pv_nucleation_dt   > 0.0)      { cfg.output.paraview_nucleation_dt   = cli_pv_nucleation_dt; }
   if (cli_pv_interseismic_dt > 0.0)      { cfg.output.paraview_interseismic_dt = cli_pv_interseismic_dt; }
   // Phase 3: free-surface slice overrides (specific value/dt win over TOML).
   if (!cli_pv_free_surface.empty())      { cfg.output.paraview_free_surface = cli_pv_free_surface; }
   if (cli_pv_free_surface_dt > 0.0)      { cfg.output.paraview_free_surface_dt = cli_pv_free_surface_dt; }
   if (cli_pv_force_fault_vtu)            { cfg.output.paraview_fault  = "vtu";  }
   if (cli_pv_force_fault_hdf5)           { cfg.output.paraview_fault  = "hdf5"; }
   if (cli_pv_force_vol_vtu)              { cfg.output.paraview_volume = "vtu";  }
   if (cli_pv_force_vol_hdf5)             { cfg.output.paraview_volume = "hdf5"; }
   if (cli_pv_force_bulk_vtu)             { cfg.output.paraview_bulk   = "vtu";  }
   if (cli_pv_force_bulk_hdf5)            { cfg.output.paraview_bulk   = "hdf5"; }
   if (cli_pv_legacy_ascii)               { cfg.output.paraview_fault_legacy_ascii = true;
                                            cfg.output.paraview_fault = "vtu"; }

   // Parity Phase 1 build guards.  Fire BEFORE any pv_out construction
   // so a bad flag combination fails at parse time, not after the run.
   if (cli_pv_force_fault_vtu && cli_pv_force_fault_hdf5)
   {
      MFEM_ABORT("--paraview-fault-vtu and --paraview-fault-hdf5 are "
                 "mutually exclusive.");
   }
   if (cli_pv_legacy_ascii && cli_pv_force_fault_hdf5)
   {
      MFEM_ABORT("--paraview-fault-legacy-ascii implies --paraview-fault-vtu; "
                 "it cannot be combined with --paraview-fault-hdf5.");
   }
   if (cli_pv_force_vol_vtu && cli_pv_force_vol_hdf5)
   {
      MFEM_ABORT("--paraview-volume-vtu and --paraview-volume-hdf5 are "
                 "mutually exclusive.");
   }
   if (cli_pv_force_bulk_vtu && cli_pv_force_bulk_hdf5)
   {
      MFEM_ABORT("--paraview-bulk-vtu and --paraview-bulk-hdf5 are "
                 "mutually exclusive.");
   }
   // ZFP + deflate combos are mutually exclusive per collection.
   if (cli_pv_vol_zfp > 0.0 && cli_pv_volume_deflate >= 0)
   {
      MFEM_ABORT("--paraview-volume-zfp-tol and --paraview-volume-deflate-level "
                 "are mutually exclusive.");
   }
   if (cli_pv_bulk_zfp > 0.0 && cli_pv_bulk_deflate >= 0)
   {
      MFEM_ABORT("--paraview-bulk-zfp-tol and --paraview-bulk-deflate-level "
                 "are mutually exclusive.");
   }
   if (cli_pv_fault_zfp > 0.0 && cli_pv_fault_deflate >= 0)
   {
      MFEM_ABORT("--paraview-fault-zfp-tol and --paraview-fault-deflate-level "
                 "are mutually exclusive.");
   }
#ifndef MFEM_USE_HDF5
   if (cli_pv_force_fault_hdf5 || cli_pv_force_vol_hdf5 || cli_pv_force_bulk_hdf5
       || cli_pv_vol_zfp > 0.0 || cli_pv_bulk_zfp > 0.0 || cli_pv_fault_zfp > 0.0
       || cli_pv_fault_deflate >= 0 || cli_pv_bulk_deflate >= 0
       || cli_pv_volume_deflate >= 0
       || cfg.output.paraview_volume == "hdf5"
       || cfg.output.paraview_bulk   == "hdf5"
       || cfg.output.paraview_fault  == "hdf5")
   {
      MFEM_ABORT("--paraview-*-hdf5 / --paraview-*-zfp-tol / "
                 "--paraview-*-deflate-level require an MFEM build with "
                 "MFEM_USE_HDF5=YES.  Re-run with --paraview-*-vtu or "
                 "switch the per-collection mode to \"vtu\" / \"off\".");
   }
#endif
#ifndef MFEM_USE_H5Z_ZFP
   if (cli_pv_vol_zfp > 0.0 || cli_pv_bulk_zfp > 0.0 || cli_pv_fault_zfp > 0.0)
   {
      MFEM_ABORT("--paraview-*-zfp-tol requires an MFEM build with "
                 "MFEM_USE_H5Z_ZFP=YES.  Use --paraview-*-deflate-level "
                 "for lossless compression instead.");
   }
#endif

   MFEM_VERIFY(cfg.law == spatial::FrictionLawKind::SlipWeakening ||
               cfg.law == spatial::FrictionLawKind::RateState,
               "spatial_dyn_driver: cfg.law must be SlipWeakening or "
               "RateState; got " << static_cast<int>(cfg.law));
   const bool is_lsw =
      (cfg.law == spatial::FrictionLawKind::SlipWeakening);
   // Phase 3: both slip_weakening and rate_state are wired (the RS branch
   // below resolves RS params, seeds equilibrium psi, and selects the aging
   // iterator via MakeFrictionIterator).

   // Phase 14: the explicit-RK time integrator is an additive, CLI-gated
   // alternative to ADER (default).  Its preconditions:
   //   (a) friction law — BOTH rate_state and slip_weakening (LSW) are now
   //       supported on the RK/Mult path.  RS uses the instantaneous `Evaluate`
   //       + the coupled-RK-on-(Q,ψ,slip) stepper `AdvanceRKCoupled_Spatial`;
   //       LSW uses the instantaneous `EvaluateLSW` (the dt→0 limit of
   //       EvaluateADER_LSW) + the coupled-RK-on-(Q,slip) stepper
   //       `AdvanceRKCoupledLSW_Spatial`.  The ONLY LSW sub-case still rejected
   //       is forced rupture (TPV26/27): is_lsw maps to FaultFrictionLaw::LSW
   //       (never LSW_ForcedRupture) at the SetFaultFrictionLaw call below, and
   //       the Mult-path dispatch aborts on LSW_ForcedRupture, so there is no
   //       reachable forced-rupture RK path here today;
   //   (b) (Phase 5, BUG-4) interior flux — central/mixed flux is now SUPPORTED
   //       on the matrix (bi-material) path too (PLAN_mixed_flux_hetero_riemann.md
   //       Phases 1-4), but central flux is non-dissipative and UNSTABLE under
   //       ADER, so `matrix + mixed_flux != none` REQUIRES RK (guard below,
   //       mirroring WaveOperator::ComputeMaxDt).  `matrix + mixed_flux == none`
   //       (bi-material Godunov upwind, dissipative) runs under either
   //       integrator; `scalar` is unaffected here (scalar mixed+ADER is caught
   //       downstream by ComputeMaxDt, P4-1).
   const bool is_rk =
      (cfg.numerics.time_integrator != spatial::TimeIntegratorKind::ADER);
   const bool matrix_mixed =
      (cfg.numerics.interior_flux == spatial::InteriorFlux::Matrix
       && cfg.numerics.mixed_flux != "none");
   // G2 (BUG-4): central flux is non-dissipative ⇒ unstable under ADER.  Use
   // the shared pure predicate (== the old `matrix_mixed && !is_rk`) so the
   // EXACT driver guard decision is table-testable without a mesh (Tests
   // 5.1/5.2 via spatial::MatrixMixedFluxUnderAder).
   MFEM_VERIFY(!spatial::MatrixMixedFluxUnderAder(cfg),
               "spatial_dyn_driver: interior_flux=\"matrix\" + mixed_flux=\""
               << cfg.numerics.mixed_flux << "\" uses the non-dissipative "
               "central flux, which is UNSTABLE under ADER; it requires an RK "
               "integrator.  Set --time-integrator rk4|rk45, or use "
               "mixed_flux=\"none\" under ADER.");

   // (LTS Phase 1) Clustered LTS couples clusters through the ADER Taylor
   // predictor, so it is ADER-only.  This is the post-CLI-merge companion to the
   // parser's config-level guard (which rejects lts + mixed_flux): the final
   // integrator is only known after the --time-integrator override is applied.
   MFEM_VERIFY(!(cfg.numerics.LtsEnabled() && is_rk),
               "spatial_dyn_driver: [numerics].lts=\"" << cfg.numerics.lts
               << "\" requires the ADER integrator (clustered LTS couples "
               "clusters through the ADER Taylor predictor); got "
               "time_integrator=" << (is_rk ? "rk4|rk45" : "ader")
               << ".  Use --time-integrator ader or set lts=\"off\".");

   // (LTS Phase 3, REVIEW R-006) The secular fault-resample projects the
   // per-macro-step Δψ / Δ|slip| over WHOLE per-face QP blocks of the GLOBAL
   // fault vector; the per-cluster LTS friction sweep advances one contiguous
   // QP range at a time, so the two are incompatible in v1.  Reject the
   // combination at parse time (a clean setup error) rather than letting the
   // range Advance abort mid-sweep from deep inside the fault loop.
   MFEM_VERIFY(!(cfg.numerics.LtsEnabled() && cli_fault_resample),
               "spatial_dyn_driver: [numerics].lts=\"" << cfg.numerics.lts
               << "\" is incompatible with --fault-resample (the secular "
               "resample projects over whole per-face QP blocks of the global "
               "fault vector, which the per-cluster LTS range sweep does not "
               "support in v1).  Drop --fault-resample or set lts=\"off\".");

   // σ_n strength floor banner string (sliver-blowup plan 2026-05-26):
   // "DISABLED" for the negative sentinel, else the value in MPa.
   std::string sigma_n_floor_banner;
   if (cfg.sigma_n_strength_floor_pa < 0.0)
   {
      sigma_n_floor_banner = "DISABLED";
   }
   else
   {
      std::ostringstream oss;
      oss << (cfg.sigma_n_strength_floor_pa / 1.0e6) << " MPa";
      sigma_n_floor_banner = oss.str();
   }

   if (rank == 0)
   {
      std::cout << "================================================\n"
                << "seas_spatial_dyn_driver — Phase 4 (rev-3)\n"
                << "================================================\n"
                << "config:           " << config_path << "\n"
                << "mesh:             " << cfg.mesh.path << "\n"
                << "fe order:         " << cfg.mesh.order << "\n"
                << "law:              "
                << (is_lsw ? "slip_weakening" : "rate_state") << "\n"
                << "stress kind:      "
                << (cfg.stress.kind == spatial::StressSourceKind::ConstantTensor
                    ? "constant_tensor"
                    : cfg.stress.kind == spatial::StressSourceKind::FaultLocalPrestress
                      ? "fault_local_prestress"
                    : cfg.stress.kind ==
                      spatial::StressSourceKind::DepthProportionalToShearModulus
                      ? "depth_proportional" : "sidecar_hdf5") << "\n"
                << "tfinal:           " << cfg.time.tfinal << " s\n"
                << "cfl:              " << cfg.numerics.cfl << "\n"
                << "ader order:       " << cfg.numerics.ader_order << "\n"
                << "time integrator:  "
                << (cfg.numerics.time_integrator
                       == spatial::TimeIntegratorKind::ADER ? "ader"
                    : cfg.numerics.time_integrator
                       == spatial::TimeIntegratorKind::RK4  ? "rk4"
                                                            : "rk45") << "\n"
                << "mixed flux:       " << cfg.numerics.mixed_flux << "\n"
                << "mixed flux contrast tol: "
                << cfg.numerics.mixed_flux_contrast_tol
                << (cfg.numerics.mixed_flux_contrast_tol < 0.0
                    ? " (disabled)" : "") << "\n"
                << "use pml:          " << (cfg.numerics.use_pml ? "yes" : "no")
                << "\n"
                << "sigma_n strength floor: " << sigma_n_floor_banner << "\n"
                << "nucleation:       "
                << (cfg.nucleation.enabled
                    ? "gradual_overstress (enabled)"
                    : "DISABLED")
                << "\n"
                << "no-sidecar mat:   " << (no_sidecar_material ? "yes" : "no")
                << "\n"
                << "dry-run:          " << (dry_run ? "yes" : "no") << "\n"
                << "ranks:            " << nprocs << "\n"
                << "================================================\n";
   }

   // -----------------------------------------------------------------
   // 3.  Restart / output-dir safety check (mirrors tpv104_driver.cpp).
   // -----------------------------------------------------------------
   if (!restart_prefix.empty())
   {
      namespace fs = std::filesystem;
      try
      {
         const fs::path restart_path(restart_prefix);
         fs::path restart_dir_path = restart_path.parent_path();
         if (restart_dir_path.empty()) { restart_dir_path = "."; }

         const fs::path restart_canonical =
            fs::weakly_canonical(restart_dir_path);
         const fs::path output_canonical =
            fs::weakly_canonical(fs::path(cfg.output.output_dir));

         if (restart_canonical == output_canonical)
         {
            if (rank == 0)
            {
               std::cerr << "ERROR: --output-dir (" << cfg.output.output_dir
                         << ") resolves to the SAME directory as the "
                         << "parent of --restart (" << restart_dir_path.string()
                         << ").  Pick a DIFFERENT --output-dir for the "
                         << "restarted run.\n";
            }
#ifdef MFEM_USE_MPI
            MPI_Finalize();
#endif
            return 3;
         }
      }
      catch (const fs::filesystem_error &e)
      {
         if (rank == 0)
         {
            std::cerr << "ERROR: failed to canonicalise --restart / "
                      << "--output-dir paths: " << e.what() << "\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 3;
      }
   }

   // -----------------------------------------------------------------
   // 4.  Load mesh.
   // -----------------------------------------------------------------
   Mesh smesh(cfg.mesh.path.c_str(), 1, 1);
   const int dim = smesh.Dimension();
   MFEM_VERIFY(dim == 3,
               "spatial_dyn_driver: only 3D meshes supported; got dim="
               << dim);

#ifdef MFEM_USE_MPI
   // Fault-locality partition (opt-in via --partition-fault-locality; default OFF
   // => byte-exact for existing spatial runs).  Forces both elements of every
   // fault face to be co-resident on one rank so the y=0 fault/material interface
   // is NOT cut across ranks.  OPTIONAL for a BI-MATERIAL fault (TPV6/7) at np>1:
   // the cross-rank neighbour-material exchange now reads the TRUE peer element's
   // material (merged from fix-cross-rank-material-exchange), so a ParMETIS-cut
   // fault is handled correctly — the per-side fault material and the mixed-flux
   // seam are right on SHARED fault faces without this partition.  Enabling it is
   // still useful: it keeps every fault face LOCAL, removing the per-face
   // cross-rank exchange on the fault altogether (this is what the native TPV
   // drivers do, dynamic/fault_locality_partition.hpp).  smesh must stay alive
   // until the ParMesh ctor; fl_part must outlive it (MFEM may store the pointer).
   Array<int> fl_part;
   int *part_data = nullptr;
   if (HasFlag(argc, argv, "--partition-fault-locality"))
   {
      const int fault_attr_part = (cfg.boundary.fault_attr > 0)
                                  ? cfg.boundary.fault_attr : 101;
      Array<int> fault_faces;
      seas::FindFaultFaceIndices(smesh, fault_attr_part, fault_faces);
      int n_relocated = 0;
      seas::BuildFaultLocalityPartitioning(smesh, fault_faces, nprocs,
                                           fl_part, &n_relocated);
      const int violations =
         seas::VerifyFaultLocality(smesh, fault_faces, fl_part);
      MFEM_VERIFY(violations == 0,
                  "spatial_dyn: fault-locality partition has " << violations
                  << " violations (partitioning logic bug).");
      if (rank == 0)
      {
         std::cout << "[partition] fault-locality ENABLED (fault_attr="
                   << fault_attr_part << "): " << fault_faces.Size()
                   << " fault faces, " << n_relocated << " elements relocated, "
                   << violations << " violations\n";
      }
      part_data = fl_part.GetData();
   }

   // === (LTS Phase 4 Step 0) Serial-mesh clustering + LTS-aware partition ====
   // Rank 0 clusters the SERIAL mesh (material evaluated on smesh BEFORE the
   // ParMesh — the "material-before-ParMesh" requirement; all spatial-driver
   // material modes are coordinate-evaluable) and broadcasts the rank-count-
   // INDEPENDENT cluster ids.  At np>1 it also builds the cluster-weighted LTS
   // partition (P-009) and injects it as `part_data` (overriding fault-locality,
   // which the LTS partition already enforces via its fault-lock).  Gated
   // lts != "off": for lts="off" nothing here runs (byte-identical).  These
   // outlive the ParMesh ctor (serial_cluster + lts_part_vec are reused after it
   // to map serial->local cluster ids and to re-base the checkpoint hash).
   std::vector<int> serial_cluster;    // per SERIAL-element cluster id
   LtsClustering    serial_cl;
   const int        serial_ne = smesh.GetNE();
   std::vector<int> lts_part_vec;      // LTS partition (int*, must outlive ctor)
   bool             serial_lts_ready = false;
   if (cfg.numerics.LtsEnabled())
   {
      // REVIEW F-1b: the serial->local element map (below the ParMesh) relies on
      // MFEM assigning local elements in ascending SERIAL order — a CONFORMING-
      // mesh guarantee; a nonconforming/AMR serial mesh would silently mis-map.
      MFEM_VERIFY(!smesh.Nonconforming(),
                  "spatial_dyn: LTS serial clustering requires a conforming serial "
                  "mesh (the nonconforming serial->local element map is not "
                  "handled).");
      // REVIEW F-1 (CRITICAL): use the EFFECTIVE fault attribute.  SAFS decks omit
      // the [boundary] block => cfg.boundary.fault_attr = -1, and the driver falls
      // back to 101 (into bc.fault_attr) only LATER.  Passing the raw -1 here empties
      // in.fault_pairs => BuildLtsClustering never forces a fault face's two elements
      // into one cluster (D-1) AND the partition fault-lock is a no-op (D-2) — both
      // silently broken on every SAFS mesh (invisible to TPV tests, which set fault=3).
      const int fault_attr_eff =
         (cfg.boundary.fault_attr > 0) ? cfg.boundary.fault_attr : 101;

      serial_cluster.assign(static_cast<std::size_t>(serial_ne), 0);
      int nc = 1; double lam = 1.0, dtb = 0.0, modeled = 0.0;
      if (nprocs > 1) { lts_part_vec.assign(static_cast<std::size_t>(serial_ne), 0); }
      if (rank == 0)
      {
       try
       {
         // Serial material on smesh (wrappers local to this scope; smat borrows
         // them and is consumed by BuildLtsMeshInputsFromMaterial below).
         std::unique_ptr<DepthProfile1DMaterial>          s_depth;
         std::unique_ptr<HalfspaceAcrossFaultMaterial>    s_half;
         std::unique_ptr<spatial::SpatialVelocityBundle>  s_vel;
         MaterialField smat = MaterialField::MakeConstant(
            cfg.material_fallback.lambda, cfg.material_fallback.mu,
            cfg.material_fallback.rho);
         const bool sc_sidecar =
            cfg.velocity.use_sidecar && !no_sidecar_material
            && cfg.material.kind != spatial::MaterialKind::DepthProfile1D;
         if (cfg.material.kind == spatial::MaterialKind::DepthProfile1D)
         {
            s_depth = MakeDepthProfile1DMaterial(cfg.material.profile_layers,
                                                 cfg.material.depth_axis);
            smat = s_depth->field;
         }
         else if (cfg.material.kind == spatial::MaterialKind::HalfspaceAcrossFault)
         {
            const auto &hs = cfg.material.halfspace;
            s_half = MakeHalfspaceAcrossFaultMaterial(
               hs.vp_near, hs.vs_near, hs.rho_near,
               hs.vp_far,  hs.vs_far,  hs.rho_far, hs.x0, hs.normal);
            smat = s_half->field;
         }
         else if (sc_sidecar)
         {
            s_vel = std::make_unique<spatial::SpatialVelocityBundle>(
               spatial::LoadSpatialVelocityBundle(cfg.velocity, smesh));
            smat = s_vel->MakeMaterialField();
         }

         LtsMeshInputs in =
            BuildLtsMeshInputsFromMaterial(smesh, smat, fault_attr_eff);
         const LtsClusteringOptions opt = MakeLtsClusteringOptions(
            cfg.numerics.lts_nc_cap, cfg.numerics.lts_merge_loss_tol,
            cfg.numerics.lts_wiggle);
         serial_cl = BuildLtsClustering(in.dt_e, smesh.ElementToElementTable(),
                                        in.fault_pairs, opt);
         serial_cluster = serial_cl.cluster;
         nc = serial_cl.num_clusters; lam = serial_cl.lambda;
         dtb = serial_cl.dt_base; modeled = serial_cl.modeled_cost;

         if (nprocs > 1)
         {
            std::vector<int> xadj, adjncy;
            BuildLtsCsrFromMesh(smesh, xadj, adjncy);
            LtsPartitionOptions popt;
            try
            {
               BuildLtsAwarePartition(serial_ne, xadj, adjncy, serial_cluster,
                                      nc, nprocs, nullptr, in.fault_pairs, popt,
                                      lts_part_vec);
            }
            catch (const std::exception &ex)
            {
               // REVIEW F-4: fall back to a RECONSTRUCTABLE explicit partition
               // (fault-locality = MFEM GeneratePartitioning + fault-lock), NOT
               // MFEM's internal default — the serial->local map below can only
               // reproduce an explicit part_data, and this preserves D-2.
               std::cerr << "[lts] LTS-aware partition failed (" << ex.what()
                         << "); falling back to the fault-locality partition.\n";
               Array<int> ff, fbp;
               seas::FindFaultFaceIndices(smesh, fault_attr_eff, ff);
               seas::BuildFaultLocalityPartitioning(smesh, ff, nprocs, fbp);
               MFEM_VERIFY(fbp.Size() == serial_ne,
                           "spatial_dyn: fallback partition size mismatch.");
               lts_part_vec.assign(fbp.GetData(), fbp.GetData() + serial_ne);
            }

            // Phase-5 realization-fraction attribution: log the achieved
            // per-cluster / fault-work imbalance of the injected partition.
            {
               const auto pib = seas::ComputeLtsPartitionImbalance(
                  serial_cluster, nc, lts_part_vec, nprocs, nullptr,
                  in.fault_pairs);
               std::cout << "[lts-partition] imbalance: overall=" << pib.overall
                         << "  worst-cluster c" << pib.worst_cluster << "="
                         << pib.worst_cluster_imbalance << "  fault=" << pib.fault
                         << "\n[lts-partition]   per-cluster:";
               for (int c = 0; c < nc; ++c)
               { std::cout << " c" << c << "=" << pib.per_cluster[c]; }
               std::cout << "\n";
            }
         }
       }
       catch (const std::exception &ex)
       {
          // REVIEW F-2: a throw here (e.g. bad_alloc from the ~1.3 GB sidecar)
          // would otherwise terminate rank 0 while peers wait at MPI_Bcast below
          // (hang).  Abort the whole job loudly instead.
          std::cerr << "[lts] serial clustering failed on rank 0: " << ex.what()
                    << "\n";
          MPI_Abort(comm, 1);
       }
      }
#ifdef MFEM_USE_MPI
      if (nprocs > 1)
      {
         MPI_Bcast(serial_cluster.data(), serial_ne, MPI_INT, 0, comm);
         int nc_b = nc; MPI_Bcast(&nc_b, 1, MPI_INT, 0, comm); nc = nc_b;
         double d3[3] = {lam, dtb, modeled};
         MPI_Bcast(d3, 3, MPI_DOUBLE, 0, comm);
         lam = d3[0]; dtb = d3[1]; modeled = d3[2];
         // The LTS (or fault-locality-fallback) partition is ALWAYS filled at
         // np>1, so it is always injected as part_data (the serial->local map
         // below requires an explicit partition).
         if (static_cast<int>(lts_part_vec.size()) != serial_ne)
         { lts_part_vec.assign(static_cast<std::size_t>(serial_ne), 0); }
         MPI_Bcast(lts_part_vec.data(), serial_ne, MPI_INT, 0, comm);
         part_data = lts_part_vec.data();   // override fault-locality
         if (rank == 0)
         {
            std::cout << "[lts] LTS-aware partition injected (Nc=" << nc
                      << ", np=" << nprocs << ")\n";
         }
      }
#endif
      serial_cl.cluster = serial_cluster;
      serial_cl.num_clusters = nc; serial_cl.lambda = lam;
      serial_cl.dt_base = dtb; serial_cl.modeled_cost = modeled;
      serial_lts_ready = true;
   }

   ParMesh pmesh(comm, smesh, part_data);
#else
#  error "spatial_dyn_driver requires MFEM_USE_MPI=YES."
#endif
   smesh.Clear();
   pmesh.SetCurvature(cfg.mesh.order);

   // -----------------------------------------------------------------
   // 5.  BoundaryConfig.
   //     Phase 8: honor the [boundary] block from the config when present
   //     (TPV meshes tag the fault as Physical Surface 3, free surface 1,
   //     absorbing 5 — NOT the SAFS .geo convention).  When [boundary] is
   //     ABSENT the parser leaves cfg.boundary.fault_attr = -1 (struct
   //     default) with empty natural/absorbing lists; fall back to the SAFS
   //     .geo convention (fault=101, top free=102, bottom+sides
   //     absorbing=103/104) so existing SAFS configs are byte-unchanged.
   // -----------------------------------------------------------------
   BoundaryConfig bc;
   if (cfg.boundary.fault_attr > 0)
   {
      bc.fault_attr      = cfg.boundary.fault_attr;
      bc.natural_attrs   = std::set<int>(cfg.boundary.natural_attrs.begin(),
                                         cfg.boundary.natural_attrs.end());
      bc.absorbing_attrs = std::set<int>(cfg.boundary.absorbing_attrs.begin(),
                                         cfg.boundary.absorbing_attrs.end());
   }
   else
   {
      bc.fault_attr      = 101;
      bc.natural_attrs   = {102};
      bc.absorbing_attrs = {103, 104};
   }

   // -----------------------------------------------------------------
   // 5b. Free-surface slice attributes (Phase 3).  Resolve against the
   //     RESOLVED `bc.natural_attrs` (a std::set<int> that already includes
   //     the SAFS {102} fallback) — NOT the raw `cfg.boundary.natural_attrs`,
   //     which is empty for every SAFS config (no [boundary] block) and would
   //     silently disable the slice on its primary target (R-001).  An
   //     explicit [output].paraview_free_surface_attrs override wins.
   // -----------------------------------------------------------------
   Array<int> fs_attrs;
   if (!cfg.output.paraview_free_surface_attrs.empty())
   {
      for (int a : cfg.output.paraview_free_surface_attrs) { fs_attrs.Append(a); }
   }
   else
   {
      for (int a : bc.natural_attrs) { fs_attrs.Append(a); }
   }

   // -----------------------------------------------------------------
   // 6.  Material — pick the right WaveOperator ctor (deviation D-1).
   // -----------------------------------------------------------------
   real_t mat_lambda = cfg.material_fallback.lambda;
   real_t mat_mu     = cfg.material_fallback.mu;
   real_t mat_rho    = cfg.material_fallback.rho;

#ifdef MFEM_USE_MPI
   // PLAN_sidecar_mpi_shared_memory_2026-07-17.md Phase 0: node-local
   // communicator for the MPI-3 shared-memory sidecar windows
   // (--sidecar-shared-mem).  MPI_COMM_NULL when the flag is off — every
   // loader then takes its classic per-rank path.  RAII guard: the comm
   // is freed on destruction UNLESS MPI_Finalize already ran (this
   // driver finalizes before main's locals unwind — see the MPIContext
   // note below; MPI_Comm_free after finalize aborts, and MPI_Finalize
   // itself reclaims the communicator).  The windows created on this
   // comm do NOT need it alive to be freed (MPI keeps an internal group
   // reference), so guard-vs-window destruction order is immaterial.
   struct NodeCommGuard
   {
      MPI_Comm comm = MPI_COMM_NULL;
      ~NodeCommGuard()
      {
         if (comm != MPI_COMM_NULL)
         {
            int finalized = 0;
            MPI_Finalized(&finalized);
            if (!finalized) { MPI_Comm_free(&comm); }
         }
      }
   } sidecar_node_comm;
   if (sidecar_shared_mem)
   {
      // key = world rank: node-local rank order mirrors MPI_COMM_WORLD,
      // so "node-local rank 0" (the loading rank) is deterministic.
      MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, /*key=*/rank,
                          MPI_INFO_NULL, &sidecar_node_comm.comm);
      int node_size = 0;
      MPI_Comm_size(sidecar_node_comm.comm, &node_size);
      if (rank == 0)
      {
         std::cout << "[sidecar] --sidecar-shared-mem: MPI-3 shared-memory "
                   << "windows ON (rank 0's node: " << node_size
                   << " rank(s) share ONE copy of each gridded sidecar)\n";
      }
   }
#endif

   // REVIEW R-008 (lifetime): `vel_bundle` (owns the Coefficient objects the
   // sidecar MaterialField points at) and `material` MUST outlive `wave_ptr`
   // — the matrix-path het ctor stores a non-owning `material_ = &material`,
   // and a Coefficient `material` borrows the bundle's coefficients.  They are
   // declared here, before `wave_ptr` (constructed ~120 lines below), so C++
   // destroys them AFTER `wave_ptr`; keep this ordering.
   std::unique_ptr<spatial::SpatialVelocityBundle> vel_bundle;
   // Phase 10 (TPV31): owns the three FunctionCoefficients that the
   // depth-profile MaterialField points at; MUST outlive `wave_ptr` and
   // `material` (same lifetime contract as `vel_bundle`).
   std::unique_ptr<DepthProfile1DMaterial> depth_profile_wrapper;
   // (Part B / B3) across-fault halfspace wrapper; same lifetime contract as
   // depth_profile_wrapper (the MaterialField borrows its FunctionCoefficients).
   std::unique_ptr<HalfspaceAcrossFaultMaterial> halfspace_wrapper;
   MaterialField material = MaterialField::MakeConstant(mat_lambda,
                                                        mat_mu, mat_rho);

   // Material-sidecar gating: TOML [velocity].use_sidecar AND no CLI
   // override.  The CLI flag always wins (it can force the constant
   // path even when the TOML asks for a sidecar load).
   const bool sidecar_requested =
      cfg.velocity.use_sidecar && !no_sidecar_material
      && cfg.material.kind != spatial::MaterialKind::DepthProfile1D;

   // REVIEW R-002: the stale "Phase H gap (D-1)" abort that used to reject
   // every sidecar_requested run was REMOVED here — Phase 9 (Stage B) wired
   // the heterogeneous WaveOperator(MaterialField) ctor + per-element flux
   // dispatch, so a non-Constant (Mode::Coefficient) sidecar material is now a
   // valid input on the interior_flux="matrix" path.  The remaining guards
   // cover both paths: the scalar path forces Mode::Constant (the
   // material.mode guard below), and the matrix branch forces non-Constant.
   if (cfg.material.kind == spatial::MaterialKind::DepthProfile1D)
   {
      // Phase 10 (TPV31): build the depth-profile (Mode::Coefficient)
      // MaterialField.  `material` borrows the wrapper's FunctionCoefficients
      // (the wrapper outlives both `material` and `wave_ptr`).  The matrix
      // (bimaterial) interior-flux path consumes it.
      depth_profile_wrapper = MakeDepthProfile1DMaterial(
         cfg.material.profile_layers, cfg.material.depth_axis);
      material = depth_profile_wrapper->field;
      if (rank == 0)
      {
         std::cout << "[material] kind=depth_profile_1d (axis='"
                   << cfg.material.depth_axis << "', "
                   << cfg.material.profile_layers.size()
                   << " layers)\n";
         for (std::size_t i = 0; i < cfg.material.profile_layers.size(); ++i)
         {
            const auto &L = cfg.material.profile_layers[i];
            std::cout << "  layer " << i
                      << ": depth=[" << L.depth_top_m << ", "
                      << L.depth_bot_m << "] m  vp=" << L.vp_ms
                      << " vs=" << L.vs_ms << " rho=" << L.rho_kgm3
                      << " interp=" << L.interp << "\n";
         }
      }
   }
   else if (cfg.material.kind == spatial::MaterialKind::HalfspaceAcrossFault)
   {
      // (Part B / B3, TPV6) build the across-fault halfspace (Mode::Coefficient)
      // MaterialField.  `material` borrows the wrapper's FunctionCoefficients
      // (the wrapper outlives both `material` and `wave_ptr`).  The matrix
      // (bimaterial) interior-flux path + B1 per-side fault assignment consume it.
      const auto &hs = cfg.material.halfspace;
      halfspace_wrapper = MakeHalfspaceAcrossFaultMaterial(
         hs.vp_near, hs.vs_near, hs.rho_near,
         hs.vp_far,  hs.vs_far,  hs.rho_far,
         hs.x0, hs.normal);
      material = halfspace_wrapper->field;
      if (rank == 0)
      {
         std::cout << "[material] kind=halfspace_across_fault\n"
                   << "  near (sign((x-x0).n)>=0): vp=" << hs.vp_near
                   << " vs=" << hs.vs_near << " rho=" << hs.rho_near << "\n"
                   << "  far                     : vp=" << hs.vp_far
                   << " vs=" << hs.vs_far << " rho=" << hs.rho_far << "\n"
                   << "  plane x0=(" << hs.x0[0] << "," << hs.x0[1] << ","
                   << hs.x0[2] << ") n=(" << hs.normal[0] << "," << hs.normal[1]
                   << "," << hs.normal[2] << ")\n";
      }
   }
   else if (sidecar_requested)
   {
      try
      {
#ifdef MFEM_USE_MPI
         vel_bundle = std::make_unique<spatial::SpatialVelocityBundle>(
            spatial::LoadSpatialVelocityBundle(cfg.velocity, pmesh,
                                               sidecar_node_comm.comm));
#else
         vel_bundle = std::make_unique<spatial::SpatialVelocityBundle>(
            spatial::LoadSpatialVelocityBundle(cfg.velocity, pmesh));
#endif
         material = vel_bundle->MakeMaterialField();
         if (rank == 0)
         {
            std::cout << "[material] loaded sidecar bundle: "
                      << spatial::ResolveSpatialVelocitySidecarPath(cfg.velocity)
                      << "\n";
         }
      }
      catch (const std::exception &e)
      {
         if (rank == 0)
         {
            std::cerr << "ERROR: failed to load velocity sidecar: "
                      << e.what() << "\n"
                      << "  Re-run with --no-sidecar-material (or set "
                      << "[velocity].use_sidecar = false in the TOML) "
                      << "to use the [material_constant_fallback] block.\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 4;
      }
   }
   else
   {
      if (rank == 0)
      {
         const char *gate_source =
            (!cfg.velocity.use_sidecar && no_sidecar_material) ? "TOML+CLI"
            : (!cfg.velocity.use_sidecar)                      ? "TOML"
            :                                                    "CLI";
         std::cout << "[material] using [material_constant_fallback] "
                   << "(gated by " << gate_source << "): "
                   << "lambda=" << mat_lambda
                   << " mu=" << mat_mu
                   << " rho=" << mat_rho << "\n";
      }
   }

   // Phase 9 (Stage B): the heterogeneous WaveOperator(MaterialField) ctor is
   // now wired.  On the SCALAR interior-flux path the scalar ctor still
   // consumes only the constant (lambda,mu,rho), so a non-Constant material
   // there would be silently downgraded — keep the guard for scalar.  On the
   // MATRIX path the het ctor consumes the MaterialField directly (a
   // non-Constant material is required and checked at the matrix branch below).
   MFEM_VERIFY(material.mode == MaterialField::Mode::Constant
               || cfg.numerics.interior_flux == spatial::InteriorFlux::Matrix,
               "spatial_dyn_driver: interior_flux=\"scalar\" requires a "
               "MaterialField::Mode::Constant material (the scalar WaveOperator "
               "ctor consumes only constant lambda/mu/rho).  Re-run with "
               "--no-sidecar-material to force the "
               "[material_constant_fallback] path, or set "
               "interior_flux=\"matrix\".");

   // -----------------------------------------------------------------
   // 7.  Construct WaveOperator — scalar vs matrix interior flux (Phase 9).
   // -----------------------------------------------------------------
   // REVIEW R-003: this driver always sub-steps (O = ader_order via the
   // IFrictionIterator below); [numerics].fault_iterator="one-shot" is not
   // implemented.  Reject it rather than silently sub-stepping under a
   // misleading label.  (Placed before the --dry-run exit so a dry-run
   // catches it.)
   MFEM_VERIFY(spatial::FaultIteratorSupported(cfg),
               "spatial_dyn_driver: [numerics].fault_iterator=\"one-shot\" is "
               "not implemented (the spatial driver always sub-steps with "
               "O=ader_order); use fault_iterator=\"substep\".");

   // -----------------------------------------------------------------
   // (LTS Phase 4 Step 0 / Phase 3 P-006) Map the SERIAL cluster ids (computed
   // before the ParMesh) to LOCAL element order, so they can be passed into the
   // operator ctor to reorder the fault-face list (P-006) and drive the layout.
   // MFEM's ParMesh ctor assigns local elements in ascending serial-element order
   // among those with partitioning[e]==rank (pmesh.cpp), so walking the serial
   // elements in order and picking part[e]==rank reproduces the local order.  At
   // np=1 (part all-0 / null) local == serial — identical to the operator-sourced
   // per-pmesh clustering it replaces.  For lts="off" nothing runs (byte-exact).
   // `lts_cl` carries the SERIAL cluster ids so the checkpoint layout hash is
   // rank-count-independent; `lts_cluster_id` / `lts_layout` are LOCAL.
   // -----------------------------------------------------------------
   std::vector<int> lts_cluster_id;    // per-LOCAL-element cluster id; empty => no reorder
   LtsClustering    lts_cl;            // = serial_cl (SERIAL ids; num_clusters/lambda)
   LtsLayout        lts_layout;
   LtsGlobalMeta    lts_meta;
   bool             lts_layout_ready = false;
   if (cfg.numerics.LtsEnabled() && serial_lts_ready)
   {
      lts_cluster_id.reserve(static_cast<std::size_t>(pmesh.GetNE()));
      for (int e = 0; e < serial_ne; ++e)
      {
         const int r = (part_data != nullptr) ? part_data[e] : 0;
         if (r == rank) { lts_cluster_id.push_back(serial_cluster[e]); }
      }
      MFEM_VERIFY(static_cast<int>(lts_cluster_id.size()) == pmesh.GetNE(),
                  "spatial_dyn: serial->local cluster map size "
                  << lts_cluster_id.size() << " != local NE " << pmesh.GetNE()
                  << " (partition / element-ordering mismatch).");
      // Local layout (faces / fault pairs on pmesh) from the LOCAL cluster ids.
      LtsMeshInputs lin =
         BuildLtsMeshInputsFromMaterial(pmesh, material, bc.fault_attr);
      lts_cl     = serial_cl;   // SERIAL ids + num_clusters/lambda/dt_base (hash)
      lts_layout = BuildLtsLayout(lts_cluster_id, serial_cl.num_clusters, lin.faces);
      // LTS Phase 4a: reduce the LOCAL per-cluster counts into the GLOBAL counts
      // that gate collectives (P-007).  At np=1 the sum hook is a no-op (identity).
      lts_meta   = ReduceGlobalMeta(
         lts_layout,
         (nprocs > 1)
            ? std::function<void(std::vector<long long>&)>(
                 [comm](std::vector<long long> &buf)
                 {
                    MPI_Allreduce(MPI_IN_PLACE, buf.data(),
                                  static_cast<int>(buf.size()),
                                  MPI_LONG_LONG, MPI_SUM, comm);
                 })
            : nullptr);
      // LTS Phase 4b: the batched exchange makes the collective count independent
      // of NUM_STATE (seam I = 1 batched call/correct; provider D(k) = ader_order
      // batched calls/predict), so the former `lts_meta.num_state == NUM_STATE`
      // pin is now vacuous and was removed (REVIEW p4b LOW) — num_state no longer
      // feeds n_collectives.
      lts_layout_ready = true;
   }
   const std::vector<int> *lts_cluster_ptr =
      lts_layout_ready ? &lts_cluster_id : nullptr;

   std::unique_ptr<WaveOperator<ParMesh>> wave_ptr;
   if (cfg.numerics.interior_flux == spatial::InteriorFlux::Scalar)
   {
      // Scalar (homogeneous Godunov) path — byte-identical to pre-Phase-9.
      // `lts_cluster_ptr` (P-006) is nullptr unless lts != "off", so lts="off"
      // is byte-identical (no reorder).
      wave_ptr = std::make_unique<WaveOperator<ParMesh>>(
                    pmesh, cfg.mesh.order,
                    material.lambda_const,
                    material.mu_const,
                    material.rho_const,
                    bc, lts_cluster_ptr);
   }
   else  // spatial::InteriorFlux::Matrix
   {
      // Heterogeneous (bimaterial) Riemann path (Phase 9 Stage B): the het
      // ctor consumes the MaterialField directly and builds owned_flux_pool_
      // + the per-face bimaterial flux matrices.  A Constant material here
      // would collapse to scalar and is forbidden by the parser (matrix
      // requires material.kind != Constant); guard loudly so a mis-built
      // homogeneous material cannot silently run as "matrix".
      //
      // Phase 10 (TPV31): [material].kind=\"depth_profile_1d\" is now wired
      // above (depth_profile_wrapper → Mode::Coefficient `material`).  The
      // matrix path requires a non-Constant material; depth_profile_1d and a
      // CVM velocity sidecar both satisfy it.
      MFEM_VERIFY(material.mode != MaterialField::Mode::Constant,
                  "spatial_dyn_driver: interior_flux=\"matrix\" requires a "
                  "non-Constant MaterialField, but the constructed material is "
                  "Mode::Constant.  Use [material].kind=\"depth_profile_1d\" "
                  "(with [[material_profile.layer]]) or a CVM velocity sidecar "
                  "([velocity].use_sidecar=true, --no-sidecar-material OFF).");
      // Phase 13: the matrix path constructs the separate
      // BimaterialWaveOperator<ParMesh> (held as a base WaveOperator<ParMesh>
      // unique_ptr).  All 18 `wave.*` driver calls dispatch through the base;
      // the material-dependent ones (FluxForElem_/InteriorFaceFlux_/
      // ApplyElementJacobian_/ComputeMaxDt/SetMixedFluxMode) resolve virtually
      // to the per-element bimaterial overrides.
      wave_ptr = std::make_unique<BimaterialWaveOperator<ParMesh>>(
                    pmesh, cfg.mesh.order, material, bc, lts_cluster_ptr);
      // (Cross-rank Phase 5) DEPRECATED: [material].seam_continuous now has NO
      // effect — the cross-rank exchange reads the TRUE peer material, so the
      // central build no longer gates on this affirmation (the abort it guarded is
      // gone).  The call is kept for config back-compat (the flag is stored, never
      // read); strong-contrast safety is the contrast guard (mixed_flux_contrast_tol).
      static_cast<BimaterialWaveOperator<ParMesh>&>(*wave_ptr)
         .SetSeamContinuous(cfg.material.seam_continuous);
   }
   WaveOperator<ParMesh> &wave = *wave_ptr;

   // -----------------------------------------------------------------
   // (LTS Phase 1) Compute + carry the clustered-LTS layout on every run.
   //   Gated on lts != "off": for lts="off" NOTHING here runs, so the driver is
   //   byte-identical to the pre-LTS global-time-stepping (GTS) path.  Phase 1
   //   still STEPS GLOBALLY — the layout is built and logged (proving the run
   //   path constructs it without error), but not yet consumed; Phase 2 promotes
   //   it to a persistent member and drives the multi-cluster tick loop from it.
   //
   //   STAGED (Phase 1b): the DETERMINISTIC serial-mesh clustering (rank 0,
   //   pre-ParMesh, so cluster ids are independent of the rank count) + the
   //   LTS-aware METIS partition + the fault-QP reorder are deferred.  They
   //   depend on building the material BEFORE the ParMesh (it is built after)
   //   plus a serial<->local element map, and their acceptance gate (rate2
   //   stations byte-identical under the QP permutation) is only validatable on
   //   the production meshes (Frontera/Expanse).  Until then the per-rank
   //   clustering here is exact only at np=1 (the go/no-go + Phase-0 report
   //   configuration): at np>1 a rank's local mesh can be disconnected and the
   //   maxdiff fixpoint is per-rank, so the build is skipped with a clear log.
   // -----------------------------------------------------------------
   if (cfg.numerics.LtsEnabled())
   {
      if (lts_layout_ready)
      {
         // Serial-mesh clustering built the layout (single-sourced with the P-006
         // fault-face reorder in the operator).  REVIEW F-6: at np>1 this is NO
         // LONGER "byte-identical to lts=off" — the reorder + LTS partition are
         // ACTIVE; only the multi-cluster stepping stays np=1 (the MPI stepper is
         // Phase 4, the fault-half interleave Phase 3), so np>1 currently steps
         // GTS but with the Step-0 machinery active (np>1 parity is Frontera).
         const LtsClustering &cl     = lts_cl;
         const LtsLayout      &layout = lts_layout;

         if (rank == 0)
         {
            const auto hist = cl.cells_per_cluster();
            std::cout << "[lts] ENABLED (mode=" << cfg.numerics.lts
                      << ", np=" << nprocs << "): Nc = " << cl.num_clusters
                      << ", lambda = " << cl.lambda
                      << ", predicted harmonic speedup "
                      << HarmonicUpdateSpeedup(cl) << "x"
                      << (wave.FaultFacesReordered()
                          ? "; fault-face reorder ACTIVE (cluster-contiguous, "
                            "P-006)" : "")
                      << (nprocs > 1
                          ? "; LTS partition ACTIVE (Step 0; stepping GTS — MPI "
                            "stepping is Phase 4, np>1 parity Frontera-validated)"
                          : "")
                      << "\n[lts]   histogram (cells per cluster c=0..):";
            for (long long h : hist) { std::cout << ' ' << h; }
            std::cout << "\n[lts]   layout (rank 0): " << layout.provider_elems.size()
                      << " provider elems, " << layout.consumer_owner_elems.size()
                      << " consumer-owner elems, " << layout.num_owned_faces()
                      << " owned faces, "
                      << layout.local_fault_faces_per_cluster.size()
                      << "-cluster fault map\n";
         }
      }
      else if (rank == 0)
      {
         std::cout << "[lts] ENABLED (mode=" << cfg.numerics.lts
                   << ") but the serial clustering was not built; running GTS "
                      "(byte-identical to lts=off).\n";
      }
   }

   // Lever 1 (ADER hot-path optimization) opt-in.  --deriv-cache precomputes the
   // per-element D_d^e = M_e^{-1} K_d^e and switches ApplySpatialDerivative to a
   // dense mat-vec (the dominant ~2/3 of step time; 4-5x on the macro-step
   // aggregate locally).  DEFAULT (flag absent) stays DerivMode::OnTheFly, so
   // every TPV*/BP5 run is byte-identical; with the flag the result changes at
   // round-off only (REVIEW R-002 — machine-eps, not bit-exact).  Works for both
   // the scalar and bimaterial (matrix) operators (the cache is geometry-only).
   // R-004: SetDerivMode aborts fail-loud if the cache exceeds the per-rank
   // budget (default 1 GiB; p2 ~18 MB/rank, p3 ~72 MB/rank for TPV31).
   const bool use_deriv_cache = HasFlag(argc, argv, "--deriv-cache");
   if (use_deriv_cache)
   {
      wave.SetDerivMode(DerivMode::Cached);
      if (rank == 0)
      {
         std::cout << "[deriv] DerivMode::Cached enabled "
                      "(precomputed D_d^e = M^-1 K_d; --deriv-cache)\n";
      }
   }

   // Opt 2026-06-24 (ADER hot-path).  --face-cache precomputes per-(interior
   // non-fault face, QP) {normal, weight, shape1, shape2} so ComputeADERFaceFluxRHS
   // reuses them instead of recomputing GetFaceElementTransformations/CalcOrtho/
   // CalcShape every macro-step.  ≤1e-12 (NOT bit-exact); default OFF.  Mutually
   // exclusive with precomputed fluxes (the SAFS driver never enables those).
   const bool use_face_cache = HasFlag(argc, argv, "--face-cache");
   if (use_face_cache)
   {
      wave.SetUseFaceCache(true);
      if (rank == 0)
      {
         std::cout << "[deriv] face geometry cache enabled "
                      "(precomputed interior-face normal/shape; --face-cache)\n";
      }
   }

   // Lever 2 (ADER hot-path optimization) opt-in.  --shared-ck-recursion makes
   // the substep dispatch compute the Cauchy-Kovalevskaya recursion ONCE per
   // macro-step (producing both the substep nodal states and the time integral)
   // instead of twice, removing the second ApplySpatialDerivative subtree.
   // Byte-identical to the default two-recursion path at a fixed DerivMode
   // (REVIEW R-002/R-003); default OFF preserves the current behaviour.
   const bool use_shared_ck = HasFlag(argc, argv, "--shared-ck-recursion");
   if (use_shared_ck && rank == 0)
   {
      std::cout << "[deriv] shared CK recursion enabled "
                   "(one predictor recursion/step; --shared-ck-recursion)\n";
   }

   // R-107 reflection-time warning + Phase 12.2 PML geometry.
   //
   // GetBoundingBox is per-rank LOCAL, so reduce to the GLOBAL box; both the
   // reflection warning's min_box_dim and the PML shell placement (built far
   // below, after SetAbsorbingBackground) need the global extent, and on a
   // partitioned ParMesh no single rank owns the whole box.  These three
   // values are hoisted to outer scope so the PML construction site can reuse
   // them without recomputing.
   Vector pml_box_lo(3), pml_box_hi(3);
   pmesh.GetBoundingBox(pml_box_lo, pml_box_hi, 1);
#ifdef MFEM_USE_MPI
   {
      real_t loc[3], glob[3];
      for (int d = 0; d < 3; ++d) { loc[d] = pml_box_lo(d); }
      MPI_Allreduce(loc, glob, 3, MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
      for (int d = 0; d < 3; ++d) { pml_box_lo(d) = glob[d]; }
      for (int d = 0; d < 3; ++d) { loc[d] = pml_box_hi(d); }
      MPI_Allreduce(loc, glob, 3, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
      for (int d = 0; d < 3; ++d) { pml_box_hi(d) = glob[d]; }
   }
#endif
   real_t pml_cp = -1.0;   // <0 sentinel = "unknown" (non-Constant material)

   // Phase 13 (REVIEW R-002): the reflection warning (and PML) are valid ONLY
   // on the scalar (Constant) path.  On the matrix path `material` is
   // Mode::Coefficient, whose lambda_const/mu_const/rho_const are 0
   // (MakeCoefficient sets only the Coefficient pointers), so cp = sqrt(0/0)
   // = NaN and the warning would fire bogusly ("cp_max (0 s)").  Skip it for
   // non-Constant materials; a per-element cp_max estimate is a follow-up.
   if (material.mode == MaterialField::Mode::Constant)
   {
      pml_cp = std::sqrt((material.lambda_const
                          + 2.0 * material.mu_const)
                         / material.rho_const);
      real_t min_box_dim = std::numeric_limits<real_t>::infinity();
      for (int d = 0; d < 3; ++d)
      {
         min_box_dim = std::min(min_box_dim, pml_box_hi(d) - pml_box_lo(d));
      }
      const real_t t_reflect = (pml_cp > 0.0) ? min_box_dim / pml_cp : 0.0;
      if (rank == 0 && t_reflect < cfg.time.tfinal && !cfg.numerics.use_pml)
      {
         std::cout << "[spatial_dyn] WARNING: tfinal (" << cfg.time.tfinal
                   << " s) exceeds min_box_dim / cp_max ("
                   << t_reflect << " s).  Reflected waves will "
                   << "contaminate the rupture after this time.  Add "
                   << "--pml or shorten tfinal.\n";
      }
   }

   // Phase N: the spatial driver supports exactly one nucleation kind
   // (`gradual_overstress`) — the friction law is always plain LSW.
   // The gradual_overstress accumulator writes time-domain perturbations
   // into DOFData::tau{1,2}_nuc; the LSW solver consumes them via
   // s.tau{1,2}_total = tau{1,2}_0 + tau{1,2}_nuc + trial.  The obsolete
   // LSW_ForcedRupture dispatch arm + f_2(t) per-DOF friction reduction
   // are NOT used.  Native TPV* drivers continue to set
   // FaultFrictionLaw::LSW_ForcedRupture verbatim.
   wave.SetFaultFrictionLaw(is_lsw ? FaultFrictionLaw::LSW
                                   : FaultFrictionLaw::RateAndState);
   // R-020: at np>1 the shared (rank-seam) fault QPs run inline EvaluateADER at
   // macro dt with END-OF-STEP psi (1st-order at ader_order>=2; the R-1601
   // fallback), and the plan-mandated np=2 psi-consistency gate (R-004) is not
   // yet validated.  Warn so a parallel RS result is not mistaken for validated.
   if (!is_lsw && nprocs > 1 && rank == 0)
   {
      std::cout << "[spatial_dyn] WARNING: rate_state at np>1 (" << nprocs
                << " ranks) — shared (rank-seam) fault QPs use end-of-step psi "
                   "(1st-order) and the np=2 psi-consistency gate (plan R-004) "
                   "is not yet validated.  Treat parallel RS results as "
                   "PRELIMINARY.\n";
   }
   // (Unified bi-material plan, Part A) Set the central-flux contrast-guard
   // tolerance BEFORE SetMixedFluxMode: the bi-material SetMixedFluxMode override
   // rebuilds central_flux_face_set_ + the per-face central matrices via
   // BuildPerFaceCentralFluxMatrices_, which reads mixed_flux_contrast_tol_ to drop
   // strong-contrast corridor faces.  Order matters (R-006).  tol < 0 = disabled.
   wave.SetMixedFluxContrastTol(cfg.numerics.mixed_flux_contrast_tol);
   wave.SetMixedFluxMode(ParseMixedFlux(cfg.numerics.mixed_flux));

   // (Phase 5, req 4) Headline banner for the combination this feature enables:
   // the bi-material (matrix) operator dispatching central flux per-face.
   if (rank == 0 && matrix_mixed)
   {
      std::cout << "[mixed-flux] matrix (bi-material) + "
                << cfg.numerics.mixed_flux << " central flux enabled "
                << "(seam_continuous=" << (cfg.material.seam_continuous
                                           ? "true" : "false")
                << "; DEPRECATED, no effect — cross-rank seam material uses the "
                   "TRUE peer)\n";
   }

   // Prove the mixed-flux mode is NOT a silent no-op: report the GLOBAL
   // count of central-flux faces actually populated by
   // BuildCentralFluxFaceSet_.  The startup banner only echoes the config
   // STRING ("mixed flux: adjacent"); this prints what the operator really
   // built.  none -> 0 (upwind everywhere); adjacent / all_continuous -> > 0
   // iff fault-adjacent non-fault interior faces exist on this mesh.
   // (Rank-summed; a shared non-fault face is counted on each rank that
   // owns it, so seam faces are slightly over-counted — fine as an
   // "is it active" signal.  Compare none vs adjacent: 0 vs N.)
   {
      long long central_local =
         static_cast<long long>(wave.GetCentralFluxFaceSet().size());
      long long central_global = central_local;
#ifdef MFEM_USE_MPI
      MPI_Reduce(&central_local, &central_global, 1, MPI_LONG_LONG,
                 MPI_SUM, 0, comm);
#endif
      if (rank == 0)
      {
         std::cout << "[mixed-flux] mode = " << cfg.numerics.mixed_flux
                   << ", central-flux faces (rank-summed) = "
                   << central_global
                   << (central_global == 0
                       ? "  (upwind everywhere — mode is a NO-OP)"
                       : "  (central flux ACTIVE on these faces)")
                   << "\n";
      }
   }

   // Phase 4 (fault-dealiasing §6): fault-flux over-integration.
   // --fault-overint K (default 0 = off = byte-exact) raises the FAULT-face
   // quadrature to degree 2*(order+K) for the friction solve + flux assembly,
   // decoupled from the bulk 2*order rule.  Placed AFTER SetMixedFluxMode so
   // SetFaultOverint's guard sees the real mixed-flux mode, and BEFORE the
   // fault-table setup below — it grows nbf_per_face_, which ProbeNbfPerFace /
   // the FaultBasis QP probe / BuildPerDOFFaultTables / SetFaultDOFData all
   // read via wave.FaultFaceQuadDegree() / wave.GetNbfPerFace().  Honoured on
   // both the scalar and matrix (BimaterialWaveOperator) paths (the fault-flux
   // routines live in the base WaveOperator).
   MFEM_VERIFY(cli_fault_overint >= 0,
               "--fault-overint: factor K must be >= 0, got "
               << cli_fault_overint);
   if (cli_fault_overint > 0)
   {
      wave.SetFaultOverint(cli_fault_overint);
      if (rank == 0)
      {
         // NB: do NOT print wave.GetNbfPerFace() here — it is the PER-RANK
         // local fault-QP count, which is 0 on a rank that owns no fault faces
         // (e.g. rank 0 in most partitions).  The authoritative global per-face
         // count is the Allreduce'd "[fault] QPs per face = N" line printed
         // after ProbeNbfPerFace below.
         std::cout << "[fault] over-integration ON (--fault-overint "
                   << cli_fault_overint << "): fault-face quad degree "
                   << wave.FaultFaceQuadDegree() << " vs baseline "
                   << 2 * cfg.mesh.order << " (per-face QP count reported "
                   "below as '[fault] QPs per face').\n";
      }
   }

   wave.SetTime(cfg.time.t_initial);

   // -----------------------------------------------------------------
   // 8.  Build per-DOF fault tables (deviation D-2: inline walk).
   // -----------------------------------------------------------------
   const Array<int> &fault_int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &fault_shr_faces = wave.GetFaultSharedFaces();

   // MFEM_USE_MPI is unconditional in this driver (#error at L682
   // requires it), so the comm arg can be passed without the
   // #ifdef-in-argument-list dance (R-607 round-6).
   const int nbf_per_face = ProbeNbfPerFace(pmesh, wave.FaultFaceQuadDegree(),
                                            fault_int_faces,
                                            fault_shr_faces, comm);

   const int num_fault_local  = fault_int_faces.Size() * nbf_per_face;
   const int num_shared_fault = fault_shr_faces.Size() * nbf_per_face;
   const int num_fault_total  = num_fault_local + num_shared_fault;

   // Reconstruct FaultBasis on this rank (interior + shared appended).
   // IMPORTANT: this basis must match WaveOperator's runtime fault basis
   // exactly, otherwise tau_pre / sigma_n / nucleation are resolved in one
   // signed local frame and consumed by the friction solve in another.
   // WaveOperator uses the BP5 / Tandem convention ref_normal = -y and
   // up = +z (wave_operator.inl ctor path).  Keep the same convention here.
   Vector ref_normal(3);
   ref_normal(0) = 0.0;  ref_normal(1) = -1.0; ref_normal(2) = 0.0;
   Vector up_vec(3);
   up_vec(0)     = 0.0;  up_vec(1)     = 0.0;  up_vec(2)     = 1.0;

   FaultBasis fbasis;
   fbasis.Compute(pmesh, fault_int_faces, ref_normal, up_vec);
#ifdef MFEM_USE_MPI
   fbasis.AppendSharedFaces(pmesh, fault_shr_faces, ref_normal, up_vec);
#endif

   // R-005: populate per-QP basis vectors so curvilinear-face fidelity is
   // not lost when `BuildPerDOFFaultTables` stamps `dof_basis(*, q)`.
   // On a rank with no fault faces these calls are no-ops; on a rank with
   // only shared faces the IR geometry is probed from the first shared
   // face transformation.  `BuildPerDOFFaultTables` falls back to the
   // centroid frame if `qp_data` is empty (e.g., planar TPV-style mesh
   // where the probe IR doesn't match the per-face IR exactly).
   {
      FaceElementTransformations *ftr_probe = nullptr;
      if (fault_int_faces.Size() > 0)
      {
         ftr_probe = pmesh.GetInteriorFaceTransformations(fault_int_faces[0]);
      }
#ifdef MFEM_USE_MPI
      else if (fault_shr_faces.Size() > 0)
      {
         ftr_probe = pmesh.GetSharedFaceTransformations(fault_shr_faces[0]);
      }
#endif
      if (ftr_probe)
      {
         // Phase 4: match WaveOperator's internal fault basis exactly — use the
         // same (possibly over-integrated) fault-face quadrature degree.
         const IntegrationRule &qp_ir =
            IntRules.Get(ftr_probe->GetGeometryType(),
                         wave.FaultFaceQuadDegree());
         fbasis.ComputeQPBasis(pmesh, fault_int_faces,
                               ref_normal, up_vec, qp_ir);
#ifdef MFEM_USE_MPI
         fbasis.ComputeQPBasisShared(pmesh, fault_shr_faces,
                                     ref_normal, up_vec, qp_ir,
                                     fault_int_faces.Size());
#endif
      }
   }

   std::vector<Vector>          fault_coords;
   Vector                       dof_coords_3d;
   DenseMatrix                  dof_basis;
   Array<int>                   dof_to_elem;
   Array<int>                   dof_to_attr;
   std::vector<IntegrationPoint> dof_ips;
   BuildPerDOFFaultTables(pmesh, wave.FaultFaceQuadDegree(), bc.fault_attr,
                          fbasis, fault_int_faces, fault_shr_faces,
                          nbf_per_face,
                          fault_coords, dof_coords_3d, dof_basis,
                          dof_to_elem, dof_to_attr, dof_ips);

   int num_fault_global = num_fault_total;
   int num_shared_global = num_shared_fault;
#ifdef MFEM_USE_MPI
   MPI_Allreduce(&num_fault_total,  &num_fault_global,  1, MPI_INT,
                 MPI_SUM, comm);
   MPI_Allreduce(&num_shared_fault, &num_shared_global, 1, MPI_INT,
                 MPI_SUM, comm);
#endif
   if (rank == 0)
   {
      std::cout << "[fault] QPs per face = " << nbf_per_face
                << ", num_fault_global = " << num_fault_global
                << " (local = " << num_fault_local
                << ", shared = " << num_shared_fault << ")\n";
   }

   // Phase 3: build the per-face resample projector R once (degree-`order` L2
   // projection on the fault-face over-integration GP set).  Affine faces share
   // one reference R (the |J_F| scale cancels), so it is built from the
   // reference rule on the fault-face geometry at wave.FaultFaceQuadDegree() —
   // the SAME rule that defines the fault QPs (nbf_per_face).  Built on every
   // rank (consistent geometry + rule ⇒ identical R, required for shared-face
   // consistency); harmlessly unused on a rank with no fault QPs.
   DenseMatrix fault_resample_R;
   if (cli_fault_resample)
   {
      Geometry::Type face_geom = Geometry::TRIANGLE;
      FaceElementTransformations *probe_ftr = nullptr;
      if (fault_int_faces.Size() > 0)
      {
         probe_ftr = pmesh.GetInteriorFaceTransformations(fault_int_faces[0]);
         if (probe_ftr) { face_geom = probe_ftr->GetGeometryType(); }
      }
#ifdef MFEM_USE_MPI
      else if (fault_shr_faces.Size() > 0)
      {
         probe_ftr = pmesh.GetSharedFaceTransformations(fault_shr_faces[0]);
         if (probe_ftr) { face_geom = probe_ftr->GetGeometryType(); }
      }
#endif
      const IntegrationRule &rs_ir =
         IntRules.Get(face_geom, wave.FaultFaceQuadDegree());
      MFEM_VERIFY(rs_ir.GetNPoints() == nbf_per_face,
                  "spatial_dyn: resample rule has " << rs_ir.GetNPoints()
                  << " QPs but nbf_per_face = " << nbf_per_face
                  << " — the resample rule must match the fault QP rule.");

      // R-003: the single reference R assumes AFFINE fault faces (|J_F| constant
      // over the face ⇒ the geometric scale cancels in the projector,
      // dynamic/fault_resample.hpp).  On a curved (isoparametric) face |J_F|
      // varies within the face and the reference-measure R is NOT the
      // physical-L2(dA) projector of plan Eq. (4.2).  Verify on the probe face
      // by sampling |J_F| at the rule's QPs (all current targets use straight-
      // sided tets ⇒ flat faces ⇒ this passes trivially).
      if (probe_ftr)
      {
         real_t jmin = std::numeric_limits<real_t>::max(), jmax = 0.0;
         for (int q = 0; q < rs_ir.GetNPoints(); ++q)
         {
            probe_ftr->SetAllIntPoints(&rs_ir.IntPoint(q));
            const real_t jw = probe_ftr->Face->Weight();   // |J_F| at this QP
            jmin = std::min(jmin, jw);
            jmax = std::max(jmax, jw);
         }
         MFEM_VERIFY(jmax - jmin <= 1e-10 * jmax,
                     "spatial_dyn: --fault-resample requires AFFINE (straight-"
                     "sided) fault faces; |J_F| varies by "
                     << (jmax - jmin) / jmax << " on the probe face, so the "
                     "single reference resample R is invalid (curved fault face). "
                     "Build R per face with per-QP |J_F| weighting to support it.");
      }

      BuildFaultResampleMatrix(face_geom, cfg.mesh.order, rs_ir, fault_resample_R);
      if (rank == 0)
      {
         std::cout << "[fault] resample ON (--fault-resample): degree-"
                   << cfg.mesh.order << " L2 projector, " << nbf_per_face
                   << " QPs/face"
                   << (cli_fault_overint > 0
                       ? ""
                       : "  (INACTIVE without --fault-overint: a true no-op)")
                   << "\n";
      }
   }

   // -----------------------------------------------------------------
   // 9.  FaultGeometry via the new BP5 ctor (Phase 5a).  Seed BP5Params
   //     with its in-class defaults so ComputeBP5Params (skipped in
   //     the new ctor) doesn't divide by zero anywhere.
   // -----------------------------------------------------------------
   BP5Params bp5_seed;
   // Use the (argc, argv) MPIContext ctor (matches the seas_driver.cpp:190
   // reference pattern).  That ctor sets `comm_ = MPI_COMM_WORLD` directly
   // and leaves `owns_comm_ = false`, so the destructor is a no-op.
   // Avoids the `MPI_Comm_free after MPI_FINALIZE was invoked` abort that
   // the comm-dup ctor would trigger when `mpi_ctx`'s destructor fires
   // on the way out of main, AFTER `MPI_Finalize` has returned.
   // `Mpi::Init` inside the ctor is guarded by `Mpi::IsInitialized()`, so
   // the raw `MPI_Init` we already called above is the actual initialiser
   // and the ctor only re-runs `Hypre::Init()` (idempotent).
   MPIContext mpi_ctx(&argc, &argv);
   FaultGeometry<ParMesh> geom(bp5_seed,
                               dof_coords_3d, dof_basis, dof_to_elem,
                               nbf_per_face, &mpi_ctx, dof_ips);

   // -----------------------------------------------------------------
   // 10. Apply stress source (Phase 3 or 3b).
   // -----------------------------------------------------------------
   if (cfg.stress.kind == spatial::StressSourceKind::ConstantTensor)
   {
      // D3.1: the TOML stores stress in the right-lateral-POSITIVE
      // convention (tension-positive Cauchy shear sigma_xy^tens), while
      // the regional StressSource is the SEAS compression-positive Cauchy
      // tensor.  For the canonical y=0 vertical strike-slip fault the
      // rotation that takes the input to the stored tensor inverts only
      // the xy off-diagonal; normals (compression-positive) and the dip
      // yz shear keep their sign.  Negating sigma_xy_pa at construction
      // makes a positive (right-lateral) input map to a positive
      // (right-lateral) on-fault tau_strike, with the stored Cauchy tensor
      // and every on-fault quantity bit-identical to the pre-D3.1 configs
      // (which stored the already-negated value).  R-007: this negation
      // and the sign flip of sigma_xy_pa in all affected configs land
      // together; the parametrized golden test_constant_tensor_sign guards
      // every flipped config so a half-applied change fails loudly.
      spatial::ConstantTensorStressSource src(cfg.stress.sigma_xx_pa,
                                              cfg.stress.sigma_yy_pa,
                                              cfg.stress.sigma_zz_pa,
                                              -cfg.stress.sigma_xy_pa,
                                              cfg.stress.sigma_yz_pa,
                                              cfg.stress.sigma_xz_pa);
      geom.ComputeParams(src,
                             cfg.stress.pore_pressure.P_p_pa,
                             cfg.stress.pore_pressure.P_p_grad_pa_per_m,
                             cfg.stress.pore_pressure.min_sigma_n_pa);
   }
   else if (cfg.stress.kind == spatial::StressSourceKind::FaultLocalPrestress)
   {
      // D3.2 (Phase 6 req 2): fault-local background pre-stress, seeded
      // DIRECTLY (no Cauchy projection).  tau2_0=tau_strike (right-lateral
      // POSITIVE), tau1_0=tau_dip, sigma_n0=sigma_n_pa-P_p (compression
      // POSITIVE).  P_p here is the uniform scalar — the fault-local path
      // has no depth-gradient input (constant background only).
      geom.ComputeParamsFaultLocal(cfg.stress.tau_strike_pa,
                                   cfg.stress.tau_dip_pa,
                                   cfg.stress.sigma_n_pa,
                                   cfg.stress.pore_pressure.P_p_pa);
      // Optional rectangular tau_strike patches (last-match-wins over the
      // patch list).  Only the strike slot is overridden; the background
      // dip / sigma_n seeded above are untouched.
      if (!cfg.stress.fault_local_patches.empty())
      {
         const auto &patches = cfg.stress.fault_local_patches;
         geom.ApplyStrikePreStressOverride(
            [&patches](real_t x, real_t y, real_t z,
                       real_t &tau_strike_out) -> bool
            {
               bool matched = false;
               for (const auto &p : patches)   // last-match-wins
               {
                  if (p.inside(x, y, z))
                  {
                     tau_strike_out = p.tau_strike_pa;
                     matched = true;
                  }
               }
               return matched;
            });
      }
   }
   else if (cfg.stress.kind ==
            spatial::StressSourceKind::DepthProportionalToShearModulus)
   {
      // Phase 10 (TPV31): depth-proportional pre-stress
      //   sigma(x,y,z) = sigma_*_per_mu * mu(x,y,z) / mu_ref
      // mu(x,y,z) comes from the current material.  Today the coordinate-only
      // mu is produced by `depth_profile_1d` (DepthProfile1DMaterial::
      // eval_at_xyz) or a Constant material; sidecar_hdf5 has no coordinate-
      // only lookup yet.
      spatial::DepthProportionalToShearModulusStressSource::MuAtFn mu_at_xyz;
      if (material.mode == MaterialField::Mode::Constant)
      {
         const real_t mu_const = material.mu_const;
         mu_at_xyz = [mu_const](real_t /*x*/, real_t /*y*/, real_t /*z*/)
                     { return mu_const; };
      }
      else if (depth_profile_wrapper != nullptr)
      {
         MFEM_VERIFY(static_cast<bool>(depth_profile_wrapper->eval_at_xyz),
                     "spatial_dyn_driver: depth-profile material has no "
                     "eval_at_xyz callback (heterogeneous_material wiring lost?).");
         auto &eval = depth_profile_wrapper->eval_at_xyz;
         mu_at_xyz = [&eval](real_t x, real_t y, real_t z) -> real_t
         {
            real_t lam, mu, rho;
            eval(x, y, z, lam, mu, rho);
            return mu;
         };
      }
      else
      {
         MFEM_ABORT("spatial_dyn_driver: [stress] kind=\"depth_proportional\" "
                    "requires [material] kind=\"constant\" or "
                    "\"depth_profile_1d\".  Coordinate-only mu lookup for "
                    "sidecar_hdf5 is not yet implemented.");
      }
      const auto &dp = cfg.stress.depth_proportional;
      // D3.1 sign convention (SAME as the constant_tensor path at the negation
      // above): the TOML stores stress in the right-lateral-POSITIVE convention,
      // but safs uses the no-flip on-fault projection (fault_geometry_safs_
      // templated.inl:104: tau2 = +t2·(S·n) = -S_xy for the canonical y=0
      // fault).  To map a positive (right-lateral) sigma_xy INPUT to a positive
      // (right-lateral) on-fault tau_strike, negate the xy off-diagonal at
      // construction — exactly as ConstantTensorStressSource is built with
      // -cfg.stress.sigma_xy_pa above.  Normals (compression-positive) and the
      // dip yz/xz shears (rotation coefficient +1) keep their sign.  Without
      // this negation TPV31 would project to a LEFT-lateral background that
      // opposes the (positive) instantaneous nucleation overstress.
      spatial::DepthProportionalToShearModulusStressSource src(
         dp.sigma_xx_per_mu, dp.sigma_yy_per_mu, dp.sigma_zz_per_mu,
         -dp.sigma_xy_per_mu, dp.sigma_yz_per_mu, dp.sigma_xz_per_mu,
         dp.mu_ref_pa, std::move(mu_at_xyz));
      geom.ComputeParams(src,
                         cfg.stress.pore_pressure.P_p_pa,
                         cfg.stress.pore_pressure.P_p_grad_pa_per_m,
                         cfg.stress.pore_pressure.min_sigma_n_pa);
   }
   else
   {
#ifdef MFEM_USE_MPI
      spatial::ApplyCsmStressSidecar(cfg.stress, geom,
                                     sidecar_node_comm.comm);
#else
      spatial::ApplyCsmStressSidecar(cfg.stress, geom);
#endif
   }
   MFEM_VERIFY(geom.HasParams(),
               "spatial_dyn_driver: stress source projection failed");

   // -----------------------------------------------------------------
   // 11. Resolve per-DOF friction parameters (Phase 1/3).  Both result
   //     structs live at outer scope so the DOF-init below sees whichever
   //     the law selected; exactly one is filled.
   // -----------------------------------------------------------------
   // Phase 3 (thermal port): open the optional 3-D friction sidecar BEFORE the
   // resolver so its readers are validated at setup.  Null when the
   // [friction.rate_state.sidecar] block is absent -- and null is exactly what
   // the resolver's legacy path expects, so no branch is needed here.
   // LSW configs carry no [friction.rate_state] block at all.
   std::shared_ptr<const spatial::RateStateSidecarFields> rs_sidecar;
   if (!is_lsw && cfg.rate_state.has_value())
   {
#ifdef MFEM_USE_MPI
      rs_sidecar = MakeRateStateSidecarFields(cfg.rate_state->sidecar,
                                              sidecar_node_comm.comm);
#else
      rs_sidecar = MakeRateStateSidecarFields(cfg.rate_state->sidecar);
#endif
      if (rs_sidecar && rank == 0)
      {
         std::cout << "[spatial_dyn] friction sidecar: "
                   << cfg.rate_state->sidecar.path << " (a='"
                   << cfg.rate_state->sidecar.a_field << "', V_w='"
                   << cfg.rate_state->sidecar.V_w_field << "', "
                   << (cfg.rate_state->sidecar.far_field_clamp
                       ? "clamp" : "abort")
                   << " out-of-hull)\n";
      }
   }

   spatial::SpatialFrictionResolver     resolver(rs_sidecar);
   spatial::SlipWeakeningPerDOFParams   lsw;  // filled iff is_lsw
   spatial::RateStatePerDOFParams       rs;   // filled iff !is_lsw
   if (is_lsw)
   {
      // R-002: the slip_weakening-block check must be gated under is_lsw —
      // the parser forbids a [friction.slip_weakening] block when
      // law="rate_state", so an RS config has cfg.slip_weakening == nullopt
      // and an un-gated guard would abort the SAFS-RS run during setup.
      MFEM_VERIFY(cfg.slip_weakening.has_value(),
                  "spatial_dyn_driver: [meta].law=slip_weakening but the "
                  "[friction.slip_weakening] block is absent in TOML.");
      lsw = resolver.ResolveSlipWeakening(*cfg.slip_weakening,
                                          dof_coords_3d, dof_to_attr);
   }
   else
   {
      MFEM_VERIFY(cfg.rate_state.has_value(),
                  "spatial_dyn_driver: [meta].law=rate_state but the "
                  "[friction.rate_state] block is absent in TOML.");
      // R-003: geom.sigma_n_per_dof() is ALREADY effective (sigma_n - P_p,
      // from ProjectFaultPreStress).  ResolveRateState's last argument is
      // sigma_n_total and it subtracts pp internally to fill rs.sigma_n_eff,
      // so pass a ZERO PorePressureSpec{} here — passing
      // cfg.stress.pore_pressure would double-subtract P_p.
      rs = resolver.ResolveRateState(
         *cfg.rate_state, dof_coords_3d, dof_to_elem, dof_to_attr,
         material, pmesh,
         spatial::PorePressureSpec{},
         geom.sigma_n_per_dof());
   }

   // -----------------------------------------------------------------
   // 12. Phase 7: build the nucleation method behind INucleationMethod.
   //     MakeNucleation dispatches on cfg.nucleation.kind (absent block ->
   //     StaticOverstress).  The SAFS+RS Gaussian path is routed through
   //     GaussianGradualOverstress, which wraps the EXACT same
   //     ResolveGradualOverstress + ApplyGradualOverstressIncrement as the
   //     pre-Phase-7 inline call — byte-identical.  `nuc->ApplyIncrement`
   //     replaces the per-sub-step nuc_cb; `nuc->ApplyOnce` (after init/restart)
   //     seeds the one-shot instantaneous patch (no-op for Gaussian/static).
   // -----------------------------------------------------------------
   // Phase 10 (TPV31 spec p. 7): per-point shear-modulus lookup so the
   // instantaneous-circular nucleation amplitude is scaled by mu(depth)/mu_ref
   // (the SAME depth-dependent mu that scales the depth-proportional background
   // stress).  Empty unless a coordinate-only mu is available; MakeNucleation
   // forwards it only to the instantaneous resolver and only when its
   // mu_ref_pa > 0, so every other config/path is byte-unchanged.
   std::function<real_t(real_t, real_t, real_t)> nuc_mu_at_xyz;
   if (material.mode == MaterialField::Mode::Constant)
   {
      const real_t mu_const = material.mu_const;
      nuc_mu_at_xyz = [mu_const](real_t, real_t, real_t) { return mu_const; };
   }
   else if (depth_profile_wrapper != nullptr &&
            static_cast<bool>(depth_profile_wrapper->eval_at_xyz))
   {
      auto &eval = depth_profile_wrapper->eval_at_xyz;
      nuc_mu_at_xyz = [&eval](real_t x, real_t y, real_t z) -> real_t
      {
         real_t lam, mu, rho;
         eval(x, y, z, lam, mu, rho);
         return mu;
      };
   }
   std::unique_ptr<INucleationMethod> nuc =
      MakeNucleation(cfg, dof_coords_3d, dof_basis, nuc_mu_at_xyz);

   // The PrintDerivedAndCheck* diagnostics + the static ParaView nucleation
   // fields read the resolved per-DOF amplitudes/radial.  Source them from the
   // active method regardless of kind (single resolve; identical data to the
   // old inline nuc_params for the Gaussian path).  R-002: the compact-circular
   // (TPV102/104) and instantaneous-circular (TPV31) kinds resolve a real
   // strike amplitude / radial factor; map those into the diagnostic struct so
   // the ParaView nuc_amplitude / nuc_radial_factor fields and the derived
   // banner are not silently zeroed for the non-Gaussian kinds.  The
   // instantaneous kind has no radial field (one-shot patch), so only
   // amplitude_strike is populated.  Static / disabled kinds leave nuc_params
   // zero-sized — those consumers already guard on `.Size() == num_fault_total`.
   spatial::GradualOverstressPerDOFParams nuc_params;
   if (auto *g = dynamic_cast<GaussianGradualOverstress *>(nuc.get()))
   {
      nuc_params = g->Params();
   }
   else if (auto *c = dynamic_cast<CompactCircularGradualOverstress *>(nuc.get()))
   {
      nuc_params.amplitude_strike = c->Params().amplitude_strike;
      nuc_params.radial           = c->Params().radial;
   }
   else if (auto *inst =
               dynamic_cast<InstantaneousOverstressCircular *>(nuc.get()))
   {
      nuc_params.amplitude_strike = inst->Params().amplitude_strike;
   }

   // R-008: warn when a non-trivial fraction of fault DOFs live on
   // shared faces.  The wave operator's shared-face EvaluateADER_LSW
   // call reads DOFData::tau{1,2}_nuc AFTER the Phase N per-substep
   // iterator has accumulated the full smoothStep increment for the
   // macrostep, so those DOFs see the perturbation as an end-of-
   // macrostep step rather than a smooth ramp (1st-order time-
   // accuracy degradation).  Quantify and warn so the user can
   // tighten dt or accept the trade-off.
   if (cfg.nucleation.enabled && num_shared_global > 0 && rank == 0)
   {
      const real_t shared_frac = (num_fault_global > 0)
         ? (static_cast<real_t>(num_shared_global)
            / static_cast<real_t>(num_fault_global))
         : 0.0;
      std::cout << "[spatial_dyn] WARNING: " << num_shared_global
                << " of " << num_fault_global << " fault DOFs ("
                << (100.0 * shared_frac) << "%) live on shared faces "
                << "and will see the gradual_overstress perturbation as "
                << "an end-of-macrostep step rather than a smooth ramp.  "
                << "Tighten dt (smaller macrostep) to reduce the "
                << "1st-order time-accuracy error on shared faces.\n";
   }

   // -----------------------------------------------------------------
   // 13. Construct FaultFaceFlux with seed scalar impedances; the
   //     per-DOF impedances are overwritten by InitializeFaultDOFs_Spatial.
   // -----------------------------------------------------------------
   // The seed impedances are immediately overwritten per-DOF by
   // InitializeFaultDOFs_Spatial (and FaultFaceFlux's scalar rho_/cp_/cs_/Zp_/
   // Zs_ members are write-only — never read in the physics).  On the matrix
   // path `material` is Mode::Coefficient, so material.*_const == 0 and the
   // naive sqrt((0 + 0) / 0) would seed NaN (REVIEW R-007).  Derive a finite
   // representative seed from element 0 (the same EvalAt pattern the bi-material
   // flux-pool builder uses); the scalar (Constant) path keeps the *_const
   // values byte-identically.
   real_t lam_seed = material.lambda_const;
   real_t mu_seed  = material.mu_const;
   real_t rho_seed = material.rho_const;
   if (material.mode != MaterialField::Mode::Constant && pmesh.GetNE() > 0)
   {
      ElementTransformation *T0 = pmesh.GetElementTransformation(0);
      const IntegrationPoint &ip0 =
         Geometries.GetCenter(pmesh.GetElementBaseGeometry(0));
      material.EvalAt(0, *T0, ip0, lam_seed, mu_seed, rho_seed);
   }
   const real_t cp_seed = std::sqrt((lam_seed + 2.0 * mu_seed) / rho_seed);
   const real_t cs_seed = std::sqrt(mu_seed / rho_seed);
   FaultFaceFlux fault_flux(rho_seed, cp_seed, cs_seed);
   // σ_n strength floor (sliver-blowup plan 2026-05-26).  Sentinel < 0 ⇒
   // disabled (each law keeps its exact current strength expression ⇒
   // byte-exact for the TPV/BP5 regressions); >= 0 floors the σ_n that
   // enters the shear strength only.
   fault_flux.SetSigmaNStrengthFloor(cfg.sigma_n_strength_floor_pa);
   wave.SetFaultFlux(&fault_flux);

   // -----------------------------------------------------------------
   // 14. Initialise per-DOF DOFData via the new free function
   //     (Phase 5c, IP-aware overload).
   //
   //     R-003 (Phase N): InitializeFaultDOFs_Spatial requires
   //     T_forced_s / t0_decay_s vectors of size num_fault_total (size-
   //     validated at spatial_setup.hpp:194-200).  Under the LSW
   //     dispatch (not LSW_ForcedRupture) these are NEVER read — supply
   //     the "never forced" sentinel so the size validator passes.  The
   //     dummy values are inert by construction.
   // -----------------------------------------------------------------
   Vector dummy_T_forced(num_fault_total);  dummy_T_forced = 1.0e9;
   Vector dummy_t0_decay(num_fault_total);  dummy_t0_decay = 0.0;

   std::vector<DOFData> dof_data;
   if (num_fault_total > 0)
   {
      if (is_lsw)
      {
         spatial::InitializeFaultDOFs_Spatial<ParMesh>(
            dof_data, num_fault_total, dof_to_elem, material, pmesh,
            lsw, geom.GetTauPre(), geom.sigma_n_per_dof(),
            dummy_T_forced, dummy_t0_decay,
            dof_ips);
      }
      else
      {
         // RS: InitializeFaultDOFs_Spatial_RS has no IP-aware overload
         // (it uses the element centroid).  SAFS material is homogeneous,
         // so centroid == IP — correct here; an IP-aware RS overload is a
         // heterogeneous-material follow-up (not needed for SAFS).
         spatial::InitializeFaultDOFs_Spatial_RS<ParMesh>(
            dof_data, num_fault_total, dof_to_elem, material, pmesh,
            rs, geom.GetTauPre(), geom.sigma_n_per_dof());
         // Phase-1 helper: overwrite the psi=0 stub with the equilibrium
         // state variable for the resolved pre-stress.
         spatial::SeedEquilibriumPsi_RS(dof_data, rs, *cfg.rate_state);
      }
   }
   wave.SetFaultDOFData(&dof_data, nbf_per_face);

   // (Part B / B1) Overwrite per-side fault impedances with the material just
   // inside each side (eps-offset rule).  No-op on the scalar operator and on a
   // fault-symmetric material (Zp_plus == Zp_minus to round-off => byte-exact);
   // gives the bi-material-fault contrast for an across-fault material (TPV6).
   // MUST run AFTER SetFaultDOFData (needs fault_face_dof_offset_) and after
   // InitializeFaultDOFs_Spatial (which seeded the single-material default).
   wave.AssignFaultSidePerMaterialImpedances(dof_data);

   // Fluctuation-Q dispatch (matches TPV205): Q_bg = 0.
   {
      real_t Q_bg[NUM_STATE] = {0};
      wave.SetAbsorbingBackground(Q_bg);
   }

   // -----------------------------------------------------------------
   // 14b. Phase 12.2 — optional absorbing PML on the box walls.
   //
   //   Gated on cfg.numerics.use_pml; when off, NO PMLLayer is built and the
   //   run is byte-identical to the pre-PML driver.  Constructed AFTER
   //   SetAbsorbingBackground (the total-Q background is a hard prerequisite:
   //   the ADER corrector damps the fluctuation Q - Q_bg, not Q itself) and
   //   held in a unique_ptr at driver scope so it outlives the time loop.  PML
   //   is integrator-agnostic — both the ADER corrector and the RK Mult path
   //   call ApplyPMLDamping internally once wave.SetPML(...) is set.  On a
   //   restart PML is stateless: bbox/cp recompute identically.
   // -----------------------------------------------------------------
   std::unique_ptr<PMLLayer> pml_layer;
   if (cfg.numerics.use_pml)
   {
      MFEM_VERIFY(wave.GetAbsorbingBackground() != nullptr,
                  "PML requires a total-Q background; call "
                  "SetAbsorbingBackground(Q_bg) before SetPML.");
      MFEM_VERIFY(pml_cp > 0.0,
                  "PML needs a scalar (Constant) material to derive c_p; the "
                  "matrix interior-flux path exposes no single c_p.  Run the "
                  "SAFS PML on interior_flux=\"scalar\".");

      // (1) Thickness: explicit pml_thickness_m wins; else pml_cells*lc_far.
      real_t L_pml = cfg.numerics.pml_thickness_m;
      if (L_pml <= 0.0)
      {
         MFEM_VERIFY(cfg.mesh.lc_far_m > 0.0,
                     "PML thickness must be derived (pml_cells * "
                     "[mesh].lc_far_m) but [mesh].lc_far_m is unset.  Set it, "
                     "or pass --pml-thickness / [numerics].pml_thickness_m.  "
                     "(Refusing to hardcode a cell size in C++.)");
         L_pml = cfg.numerics.pml_cells * cfg.mesh.lc_far_m;
      }
      MFEM_VERIFY(L_pml > 0.0, "PML thickness resolved to a non-positive value");

      // (2) Half-face mask: damp x±, y±, and z_min if pml_damp_bottom; the
      //     free surface at z=z_max stays undamped unless pml_damp_top.
      int face_mask = PMLLayer::FaceXLo | PMLLayer::FaceXHi
                    | PMLLayer::FaceYLo | PMLLayer::FaceYHi;
      if (cfg.numerics.pml_damp_bottom) { face_mask |= PMLLayer::FaceZLo; }
      if (cfg.numerics.pml_damp_top)    { face_mask |= PMLLayer::FaceZHi; }
      if (cfg.numerics.pml_damp_top && rank == 0)
      {
         std::cout << "[spatial_dyn] WARNING: pml_damp_top=true damps the free "
                      "surface at z=z_max — unphysical for a half-space SAFS "
                      "run.  Leave it false unless you know why.\n";
      }

      // (3) Fault->PML clearance guard.  dof_coords_3d holds this rank's fault
      //     DOF coords (3*i+d); reduce to the GLOBAL fault bbox and require
      //     every DAMPED inner edge to clear it (else the PML would arrest
      //     slip).  Empty-fault ranks contribute the infinite sentinels.
      const real_t INF = std::numeric_limits<real_t>::infinity();
      real_t f_lo[3] = { INF, INF, INF }, f_hi[3] = { -INF, -INF, -INF };
      const int n_fault_dofs = dof_coords_3d.Size() / 3;
      for (int i = 0; i < n_fault_dofs; ++i)
      {
         for (int d = 0; d < 3; ++d)
         {
            const real_t v = dof_coords_3d(3 * i + d);
            f_lo[d] = std::min(f_lo[d], v);
            f_hi[d] = std::max(f_hi[d], v);
         }
      }
#ifdef MFEM_USE_MPI
      {
         real_t g[3];
         MPI_Allreduce(f_lo, g, 3, MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
         for (int d = 0; d < 3; ++d) { f_lo[d] = g[d]; }
         MPI_Allreduce(f_hi, g, 3, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
         for (int d = 0; d < 3; ++d) { f_hi[d] = g[d]; }
      }
#endif
      const bool have_fault = std::isfinite(f_lo[0]) && (f_hi[0] >= f_lo[0]);

      // Per-damped-face clearance from the fault to the PML inner edge.
      struct FaceClear { const char *name; int bit; real_t clearance; };
      const FaceClear faces[6] = {
         {"x_lo", PMLLayer::FaceXLo, f_lo[0] - (pml_box_lo(0) + L_pml)},
         {"x_hi", PMLLayer::FaceXHi, (pml_box_hi(0) - L_pml) - f_hi[0]},
         {"y_lo", PMLLayer::FaceYLo, f_lo[1] - (pml_box_lo(1) + L_pml)},
         {"y_hi", PMLLayer::FaceYHi, (pml_box_hi(1) - L_pml) - f_hi[1]},
         {"z_lo", PMLLayer::FaceZLo, f_lo[2] - (pml_box_lo(2) + L_pml)},
         {"z_hi", PMLLayer::FaceZHi, (pml_box_hi(2) - L_pml) - f_hi[2]},
      };
      if (have_fault)
      {
         for (const auto &fc : faces)
         {
            if (!(face_mask & fc.bit)) { continue; }
            MFEM_VERIFY(fc.clearance > 0.0,
                        "PML shell reaches the fault on the " << fc.name
                        << " wall (clearance " << fc.clearance << " m <= 0). "
                        "The PML would arrest slip.  Reduce pml_cells / "
                        "pml_thickness_m, disable that face, or remesh with a "
                        "larger buffer.");
         }
      }
      else if (rank == 0)
      {
         std::cout << "[spatial_dyn] WARNING: PML clearance check skipped — no "
                      "fault DOFs found globally.\n";
      }

      // (4) Construct + wire.  dirs=7 is ignored (an explicit mask is passed);
      //     pml_target_R is the EFFECTIVE reflection (see pml_layer.cpp Eq.17).
      pml_layer = std::make_unique<PMLLayer>(
                     pml_box_lo, pml_box_hi, L_pml, pml_cp,
                     cfg.numerics.pml_target_R, 7, face_mask);
      wave.SetPML(pml_layer.get());

      // (5) Banner — greppable "PML: ACTIVE" token for the sbatch guard.
      if (rank == 0)
      {
         auto on = [&](int bit) { return (face_mask & bit) ? "on" : "off"; };
         // R-002: pml_target_R is the INPUT R_0; the cubic profile uses the
         // n=2 prefactor (pml_layer.cpp Eq.17), so the REALIZED reflection is
         // R_eff = R_0^(3/4).  Print both so the knob is not read as R_eff.
         std::cout << "[spatial_dyn] PML: ACTIVE (L=" << L_pml << " m"
                   << ", R_input=" << cfg.numerics.pml_target_R
                   << " (R_eff~=" << std::pow(cfg.numerics.pml_target_R, 0.75)
                   << ")"
                   << ", d_max=" << pml_layer->GetDmax() << " 1/s, faces="
                   << "x_lo:" << on(PMLLayer::FaceXLo)
                   << " x_hi:" << on(PMLLayer::FaceXHi)
                   << " y_lo:" << on(PMLLayer::FaceYLo)
                   << " y_hi:" << on(PMLLayer::FaceYHi)
                   << " z_lo:" << on(PMLLayer::FaceZLo)
                   << " z_hi:" << on(PMLLayer::FaceZHi) << ")\n";
         std::cout << "             box=[" << pml_box_lo(0) << "," << pml_box_hi(0)
                   << "]x[" << pml_box_lo(1) << "," << pml_box_hi(1)
                   << "]x[" << pml_box_lo(2) << "," << pml_box_hi(2) << "] m\n";
         std::cout << "             inner edges: x[" << pml_box_lo(0) + L_pml
                   << "," << pml_box_hi(0) - L_pml << "] y["
                   << pml_box_lo(1) + L_pml << "," << pml_box_hi(1) - L_pml
                   << "] z[" << pml_box_lo(2) + L_pml << ","
                   << pml_box_hi(2) - L_pml << "]\n";
         if (have_fault)
         {
            std::cout << "             fault->PML clearances [m]:";
            for (const auto &fc : faces)
            {
               if (face_mask & fc.bit)
               { std::cout << " " << fc.name << "=" << fc.clearance; }
            }
            std::cout << "\n";
         }
         if (!(face_mask & PMLLayer::FaceZHi))
         {
            std::cout << "             free surface z=" << pml_box_hi(2)
                      << " UNDAMPED (half-space).\n";
         }
      }
   }

   // -----------------------------------------------------------------
   // 15. CFL / Δt and derived numbers.  ComputeMaxDt is the scalar-
   //     material implementation (deviation D-1).
   //
   //     The ADER-DG stable step carries a 1/(2N+1) spatial-order factor
   //     (N = cfg.mesh.order) plus a safety margin; ComputeMaxDt applies
   //     NEITHER (only its internal mixed-flux factor).  Mirror the
   //     verified tpv205_driver.cpp:1242 conversion so cfg.numerics.cfl
   //     is a TPV205-style safety knob, not the raw CFL number.  Passing
   //     0.5 raw stepped at 9x TPV205's dt -> above the N=1 stability
   //     boundary -> nucleation-end velocity blow-up (job 7743554).
   //
   //     REVIEW R-002 / plan D2: the DG factor is now opt-in via
   //     [numerics].cfl_safety.  CflSafetyFactor returns 1/(3·(2N+1)) for
   //     "dg" (byte-identical to the previous unconditional hardcode and to
   //     the gold) and 1.0 for "raw" (experimental escape hatch).
   // -----------------------------------------------------------------
   // Phase 14.4: the RK path replaces the ADER 1/(3(2N+1)) de-rating with the
   // RK imaginary-axis stability bound (spatial::RkCflFactor) and flips the
   // operator's CFL switch to the RK-calibrated mixed-flux factors
   // (SetCflRkAware).  For mixed_flux=none the ADER branch is byte-unchanged
   // (CflSafetyFactor + the default cfl_rk_aware_=false), so `--time-integrator
   // ader` keeps the pre-Phase-14 dt to the bit.
   // (Phase 4, P4-1) For mixed_flux != none, ComputeMaxDt now ABORTS under ADER
   // (central flux is non-dissipative => unstable under ADER's stability
   // region) — central/mixed flux REQUIRES an RK integrator.  A bare ADER run
   // of a mixed_flux=adjacent config (e.g. the spec-default tpv102/205_spatial
   // TOMLs, which real jobs CLI-override to rk45 or --mixed-flux none)
   // therefore fails loud with an actionable message rather than running an
   // unstable scheme.
   if (is_rk) { wave.SetCflRkAware(true); }
   const real_t dt_cfl =
      is_rk
      ? wave.ComputeMaxDt(cfg.numerics.cfl * spatial::RkCflFactor(cfg))
      : wave.ComputeMaxDt(cfg.numerics.cfl * spatial::CflSafetyFactor(cfg));
   real_t dt = (cfg.time.dt_initial > 0.0)
                ? cfg.time.dt_initial : dt_cfl;
   // R-006: warn if the user override exceeds the explicit CFL bound.
   // The ADER substep iterator is conditionally stable in dt; silently
   // crossing CFL produces nucleation-like blow-up with no error
   // pointing at the misconfiguration.
   if (cfg.time.dt_initial > 0.0 && cfg.time.dt_initial > dt_cfl)
   {
      if (rank == 0)
      {
         std::cerr << "[spatial_dyn] WARNING: [time].dt_initial = "
                   << cfg.time.dt_initial << " s exceeds CFL bound "
                   << dt_cfl << " s (cfl = " << cfg.numerics.cfl
                   << ").  ADER will be unconditionally unstable.  "
                   << "Tighten dt_initial or raise cfl.\n";
      }
   }
   if (cfg.time.dt_max > 0.0 && dt > cfg.time.dt_max)
   {
      if (rank == 0)
      {
         std::cout << "[time] tightened by user dt_max: dt_cfl = "
                   << dt_cfl << " s -> dt = " << cfg.time.dt_max << " s\n";
      }
      dt = cfg.time.dt_max;
   }
   const int nsteps = (cfg.time.tfinal > 0.0)
                       ? static_cast<int>(std::ceil(cfg.time.tfinal / dt)) : 0;
   if (rank == 0)
   {
      std::cout << "[time] dt_cfl = " << dt_cfl << " s\n"
                << "[time] dt     = " << dt     << " s\n"
                << "[time] nsteps = " << nsteps << "\n";
      // Lever A (2026-06-30): make the CFL safety breakdown legible so a
      // cfl_dg_safety sweep is auditable in the log.  Dg/ADER path only.
      if (!is_rk && cfg.numerics.cfl_safety == spatial::CflSafety::Dg)
      {
         // dt scales as 1/cfl_dg_safety, so the cfl_dg_safety=1.0 (SeisSol)
         // step is the current dt_cfl times the current cfl_dg_safety.
         const real_t seissol_equiv = dt_cfl * cfg.numerics.cfl_dg_safety;
         std::cout << "[time] cfl = " << cfg.numerics.cfl
                   << " (Courant number, <= 1); cfl_dg_safety = "
                   << cfg.numerics.cfl_dg_safety
                   << " (extra DG margin; 1.0 = SeisSol-equivalent)\n"
                   << "[time]   SeisSol-equivalent dt_cfl (cfl_dg_safety=1.0) = "
                   << seissol_equiv << " s\n";
      }
   }

   if (print_derived)
   {
      // Compute mesh h_min for the L_nuc / h_min ratio.
      real_t h_min_local = std::numeric_limits<real_t>::infinity();
      for (int e = 0; e < pmesh.GetNE(); ++e)
      {
         h_min_local = std::min(h_min_local, pmesh.GetElementSize(e, 0));
      }
      real_t h_min_global = h_min_local;
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&h_min_local, &h_min_global, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
#endif
      spatial::PrintDerivedConfig pd_cfg;
      pd_cfg.enabled                = true;
      pd_cfg.abort_on_failure       = true;
      pd_cfg.outside_safety_factor  = 3.0;

      // Phase 11b/c: depth-profile summary (informational).  The L_nuc / f_ss
      // ranges printed by PrintDerivedAndCheckRS already reflect the depth-
      // varying a/b through rs.a(i)/rs.b(i); this block just echoes the source
      // CSVs, the coverage vs the mesh fault depth, and the VW->VS transition.
      if (!is_lsw && cfg.rate_state && cfg.rate_state->depth_profile.enabled)
      {
         const auto& dp = cfg.rate_state->depth_profile;
         real_t fault_depth_max_local = 0.0;
         for (int i = 0; i < dof_coords_3d.Size() / 3; ++i)
         {
            const real_t z = dof_coords_3d(3 * i + 2);
            fault_depth_max_local =
               std::max(fault_depth_max_local, std::max(real_t(0.0), -z));
         }
         real_t fault_depth_max = fault_depth_max_local;
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&fault_depth_max_local, &fault_depth_max, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
         // R-005 (Phase 11 review): reduce the RESOLVED per-DOF a/b range across
         // ranks so the summary reports the actual seeded a/b, not just the CSV
         // depth extents.  Collective — runs on every rank before the rank-0
         // print (empty-fault ranks contribute +inf/-inf, harmless for MIN/MAX).
         real_t a_min_l =  std::numeric_limits<real_t>::infinity();
         real_t a_max_l = -std::numeric_limits<real_t>::infinity();
         real_t b_min_l =  std::numeric_limits<real_t>::infinity();
         real_t b_max_l = -std::numeric_limits<real_t>::infinity();
         for (int i = 0; i < rs.a.Size(); ++i)
         {
            a_min_l = std::min(a_min_l, rs.a(i));
            a_max_l = std::max(a_max_l, rs.a(i));
            b_min_l = std::min(b_min_l, rs.b(i));
            b_max_l = std::max(b_max_l, rs.b(i));
         }
         real_t a_min = a_min_l, a_max = a_max_l, b_min = b_min_l, b_max = b_max_l;
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&a_min_l, &a_min, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
         MPI_Allreduce(&a_max_l, &a_max, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
         MPI_Allreduce(&b_min_l, &b_min, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
         MPI_Allreduce(&b_max_l, &b_max, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
         if (rank == 0)
         {
            const auto& acv = dp.profile.a_of_depth;
            const auto& amb = dp.profile.amb_of_depth;
            std::cout << "[derived] depth-profile ENABLED (Phase 11b):\n"
                      << "[derived]   param_a_csv         = " << dp.param_a_csv << "\n"
                      << "[derived]   param_a_minus_b_csv = " << dp.param_a_minus_b_csv << "\n"
                      << "[derived]   depth scale         = " << dp.depth_to_m
                      << " m per CSV depth unit\n"
                      << "[derived]   a   depth range = [" << acv.x.front() << ", "
                      << acv.x.back() << "] m\n"
                      << "[derived]   a-b depth range = [" << amb.x.front() << ", "
                      << amb.x.back() << "] m\n"
                      << "[derived]   fault max depth (mesh) = " << fault_depth_max
                      << " m  (profile flat-clamped beyond its sampled range)\n";
            // R-005 (Phase 11 review): the CSV sample tables + the resolved a/b
            // range over the fault DOFs (the plan §11c summary requirement).
            std::cout << "[derived]   a(depth_m) knots    :";
            for (std::size_t k = 0; k < acv.x.size(); ++k)
            {
               std::cout << " (" << acv.x[k] << "->" << acv.y[k] << ")";
            }
            std::cout << "\n[derived]   (a-b)(depth_m) knots:";
            for (std::size_t k = 0; k < amb.x.size(); ++k)
            {
               std::cout << " (" << amb.x[k] << "->" << amb.y[k] << ")";
            }
            std::cout << "\n[derived]   resolved a over fault DOFs in ["
                      << a_min << ", " << a_max << "]\n"
                      << "[derived]   resolved b over fault DOFs in ["
                      << b_min << ", " << b_max << "]\n";
            // R-029: a(z) and (a-b)(z) flat-clamp on their OWN depth grids.  If
            // the fault reaches deeper than the shallower CSV's last sample,
            // b = a - (a-b) there mixes a flat-clamped (constant) curve with a
            // still-varying one — a modeling choice worth flagging.  Mesh-aware
            // so it stays silent when the fault is shallower than both CSVs
            // (the shipped config: fault ~16.5 km, CSVs to 52.7/60 km).
            const real_t shallower_max = std::min(acv.x.back(), amb.x.back());
            if (acv.x.back() != amb.x.back() && fault_depth_max > shallower_max)
            {
               const char* clamped = (acv.x.back() < amb.x.back()) ? "a" : "a-b";
               const char* varying = (acv.x.back() < amb.x.back()) ? "a-b" : "a";
               std::cout << "[derived]   WARNING (R-029): the two CSV depth ranges "
                            "differ and the fault (" << fault_depth_max
                         << " m) is deeper than the shallower CSV (" << shallower_max
                         << " m); below that depth b = a-(a-b) mixes a flat-clamped "
                         << clamped << "(z) with a varying " << varying
                         << "(z) — extend " << clamped
                         << "'s CSV to match if unintended.\n";
            }
            bool found = false;
            for (std::size_t k = 1; k < amb.x.size(); ++k)
            {
               const bool straddles =
                  (amb.y[k - 1] < 0.0) != (amb.y[k] < 0.0);
               if (straddles)
               {
                  const real_t zc = amb.x[k - 1]
                     + (0.0 - amb.y[k - 1]) / (amb.y[k] - amb.y[k - 1])
                       * (amb.x[k] - amb.x[k - 1]);
                  std::cout << "[derived]   VW->VS transition (a-b=0) at depth = "
                            << zc << " m";
                  if (zc > fault_depth_max)
                  {
                     std::cout << "  (BELOW the fault bottom: the meshed fault is "
                                  "entirely velocity-weakening)";
                  }
                  std::cout << "\n";
                  found = true;
                  break;
               }
            }
            if (!found)
            {
               std::cout << "[derived]   (a-b) does not cross 0 within the "
                            "sampled range — fault is entirely "
                         << (amb.y.front() >= 0.0 ? "velocity-strengthening"
                                                  : "velocity-weakening")
                         << "\n";
            }
         }
      }
      // Phase 3 (thermal port): the sidecar-sourced a / V_w summary.  Collective
      // (MPI_Allreduce inside), so it must sit OUTSIDE any rank==0 branch and
      // every rank must reach it.  `rs` is only filled on the !is_lsw path.
      if (!is_lsw && rs_sidecar)
      {
         PrintFrictionSidecarSummary(rs, *rs_sidecar, cfg.hypocenter,
                                     cfg.nucleation
#ifdef MFEM_USE_MPI
                                     , comm
#endif
                                     , rank);
      }

      // PLAN DEVIATION (documented): PrintDerivedAndCheck is LSW-only (it reads
      // lsw.mu_s/mu_d/d_c) and would abort on an RS run; the RS overload
      // PrintDerivedAndCheckRS prints RS-grounded derived quantities (L_nuc =
      // mu*Dc/((b-a)*sigma_n), steady-state friction f_ss, RS nucleation gate).
      if (is_lsw)
      {
         (void)spatial::PrintDerivedAndCheck(
            pd_cfg, lsw,
            geom.GetTauPre(), geom.sigma_n_per_dof(),
            dof_coords_3d,
            cfg.nucleation, nuc_params,
            cfg.stress,
            /*mu_bulk=*/material.mu_const,
            cp_seed, cs_seed,
            h_min_global,
            dt_cfl, cfg.time.tfinal,
            num_fault_global,
            geom.NumZeroNormalFallbacks()
#ifdef MFEM_USE_MPI
            , comm
#endif
            , rank);
      }
      else
      {
         (void)spatial::PrintDerivedAndCheckRS(
            pd_cfg, rs,
            geom.GetTauPre(), geom.sigma_n_per_dof(),
            dof_coords_3d,
            cfg.nucleation, nuc_params,
            cfg.stress,
            /*mu_bulk=*/material.mu_const,
            cp_seed, cs_seed,
            h_min_global,
            dt_cfl, cfg.time.tfinal,
            num_fault_global,
            geom.NumZeroNormalFallbacks()
#ifdef MFEM_USE_MPI
            , comm
#endif
            , rank);
      }
   }

   // --verify-dispatch: print the selected friction / iterator / nucleation
   // dispatch so a --dry-run can confirm the law without running the time
   // loop.  Diagnostic only (gated on the flag; no behaviour change off).
   if (verify_dispatch && rank == 0)
   {
      std::cout << "[verify-dispatch] friction law  : "
                << (is_lsw ? "LSW (slip-weakening)"
                           : "RateAndState (aging)") << "\n"
                << "[verify-dispatch] friction iter : "
                << (is_lsw
                       ? "LswFrictionIterator (Tpv205 closed-form)"
                       : "RateStateAgingFrictionIterator (Tpv102 aging, Brent)")
                << "\n"
                << "[verify-dispatch] nucleation    : "
                << (cfg.nucleation.enabled ? "gradual_overstress (enabled)"
                                           : "none (disabled)") << "\n"
                << "[verify-dispatch] interior flux : "
                << (cfg.numerics.interior_flux == spatial::InteriorFlux::Matrix
                       ? "matrix (heterogeneous bimaterial Riemann)\n"
                       : "scalar (homogeneous Godunov)\n");
   }

   // (Part C / B1, R-001) Per-side fault impedance summary: confirms whether the
   // fault is BI-MATERIAL (Zp_plus != Zp_minus, e.g. TPV6/7) or symmetric
   // (Zp_plus == Zp_minus, e.g. TPV31/TPV205).  Scans the constructed dof_data
   // (after AssignFaultSidePerMaterialImpedances).  GLOBAL: all ranks accumulate
   // their local fault DOFs and MPI_Reduce to rank 0 (verify_dispatch is
   // rank-uniform), so the banner is correct under MPI (not rank-0-local).
   if (verify_dispatch)
   {
      long long n_bimat_l = 0, n_dof_l = static_cast<long long>(dof_data.size());
      real_t zpmn_l = std::numeric_limits<real_t>::max(), zpmx_l = 0.0;
      real_t zmmn_l = std::numeric_limits<real_t>::max(), zmmx_l = 0.0;
      for (const auto &d : dof_data)
      {
         zpmn_l = std::min(zpmn_l, d.Zp_plus);  zpmx_l = std::max(zpmx_l, d.Zp_plus);
         zmmn_l = std::min(zmmn_l, d.Zp_minus); zmmx_l = std::max(zmmx_l, d.Zp_minus);
         const real_t denom = std::max(real_t(1.0),
                                       std::max(std::abs(d.Zp_plus), std::abs(d.Zp_minus)));
         if (std::abs(d.Zp_plus - d.Zp_minus) > 1.0e-9 * denom) { ++n_bimat_l; }
      }
      long long n_bimat = n_bimat_l, n_dof = n_dof_l;
      real_t zpmn = zpmn_l, zpmx = zpmx_l, zmmn = zmmn_l, zmmx = zmmx_l;
#ifdef MFEM_USE_MPI
      MPI_Reduce(&n_bimat_l, &n_bimat, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
      MPI_Reduce(&n_dof_l,   &n_dof,   1, MPI_LONG_LONG, MPI_SUM, 0, comm);
      MPI_Reduce(&zpmn_l, &zpmn, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
      MPI_Reduce(&zpmx_l, &zpmx, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
      MPI_Reduce(&zmmn_l, &zmmn, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
      MPI_Reduce(&zmmx_l, &zmmx, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
#endif
      if (rank == 0 && n_dof > 0)
      {
         std::cout << "[verify-dispatch] fault per-side : "
                   << (n_bimat > 0 ? "BI-MATERIAL fault" : "symmetric fault")
                   << " (" << n_bimat << "/" << n_dof
                   << " DOFs with Zp_plus != Zp_minus); Zp_plus in ["
                   << zpmn << ", " << zpmx << "], Zp_minus in ["
                   << zmmn << ", " << zmmx << "]\n";
      }
   }

   // -----------------------------------------------------------------
   // 15b. LTS cluster report (--lts-report; report-and-exit, runs just
   //      before the dry-run exit which --lts-report implies).  Clusters
   //      THIS rank's local mesh with our own per-element CFL dt and prints
   //      the RAW (SeisSol-comparable) and PRODUCTION histograms plus the
   //      predicted element-update speedups.  Run at np=1 for the GLOBAL
   //      histogram (SeisSol's published counts are global).
   // -----------------------------------------------------------------
   if (lts_report)
   {
      // REVIEW D-1: the report clusters THIS rank's LOCAL mesh, but at np>1 a
      // METIS subdomain can be disconnected — an all-coarse island with no
      // adjacency path to the rank's finest element bins to cluster>=2, leaving a
      // gap in the used ids, and BuildLtsClustering's contiguity assert (which
      // holds only for a connected mesh) would MFEM_ABORT with a misleading
      // "maxdiff bug".  So run the clustering ONLY at np==1 (where the single
      // SAFS body is connected and the histogram is also the GLOBAL, SeisSol-
      // comparable one); at np>1 tell the user to rerun serially and skip.  This
      // mirrors the run-path block's nprocs==1 guard.
      if (nprocs > 1)
      {
         if (rank == 0)
         {
            std::cout << "[lts-report] SKIPPED: running on " << nprocs
                      << " ranks.  The report clusters each rank's LOCAL mesh, "
                         "which at np>1 can be disconnected (a graceful global "
                         "histogram needs the Phase-1b serial clustering).  "
                         "Rerun with np=1 to reproduce SeisSol's GLOBAL cluster "
                         "counts.\n";
         }
      }
      else
      {
         const int neL = pmesh.GetNE();
         LtsMeshInputs in = BuildLtsMeshInputs(wave, pmesh);
         const bool per_elem_mat = in.per_elem_material;

         // Recover the true per-element dt scale (seconds) from the already-
         // computed dt_cfl so the report reflects the run's real timestep WITHOUT
         // the (private) MixedFluxCflFactor_() constant.  BuildLtsMeshInputs
         // returns the scale-free ratios r_e = h_e / c_p,e; multiply by
         // dt_cfl / min(r_e) (exact at np==1, where r_min is the global min).
         std::vector<double> dt_e = in.dt_e;
         double r_min = std::numeric_limits<double>::infinity();
         for (double r : dt_e) { r_min = std::min(r_min, r); }
         const double Kscale = static_cast<double>(dt_cfl) / r_min;
         for (double &v : dt_e) { v *= Kscale; }

         const std::vector<std::pair<int, int>> &fault_pairs = in.fault_pairs;
         const Table &e2e = pmesh.ElementToElementTable();
         const int nc_cap_cli = GetIntArg(argc, argv, "--lts-nc-cap", 6);

         auto print_report = [&](const char *tag, const LtsClustering &cl)
         {
            const auto hist = cl.cells_per_cluster();
            std::cout << "[lts-report] " << tag << ": Nc = " << cl.num_clusters
                      << ", lambda = " << cl.lambda
                      << ", dt_base = " << cl.dt_base << " s\n"
                      << "[lts-report]   histogram (cells per cluster c=0..):";
            for (long long h : hist) { std::cout << ' ' << h; }
            std::cout << "\n[lts-report]   element-update speedup: harmonic "
                      << HarmonicUpdateSpeedup(cl)
                      << "x (honest workload ratio), arithmetic-mean "
                      << ArithmeticMeanSpeedup(cl)
                      << "x (SeisSol-printed form)\n";
         };

         // RAW mode: lambda = 1, no wiggle, no merge — SeisSol's exact algorithm.
         LtsClusteringOptions raw_opt;
         raw_opt.wiggle_scan  = false;
         raw_opt.lambda_fixed = 1.0;
         raw_opt.nc_cap       = 0;   // <= 0 disables the auto-merge
         const LtsClustering raw =
            BuildLtsClustering(dt_e, e2e, fault_pairs, raw_opt);

         // PRODUCTION mode: lambda-scan + Nc-cap auto-merge (the run defaults).
         LtsClusteringOptions prod_opt;
         prod_opt.wiggle_scan    = true;
         prod_opt.nc_cap         = nc_cap_cli;
         prod_opt.merge_loss_tol = 0.05;
         const LtsClustering prod =
            BuildLtsClustering(dt_e, e2e, fault_pairs, prod_opt);

         if (rank == 0)
         {
            std::cout << "[lts-report] local NE = " << neL
                      << ", fault-face pairs = " << fault_pairs.size()
                      << ", material path = "
                      << (per_elem_mat ? "matrix (per-element c_p)"
                                       : "scalar (uniform c_p)")
                      << ", dt_cfl = " << dt_cfl << " s\n";
            print_report("RAW (SeisSol algorithm: lambda=1, no merge)", raw);
            print_report("PRODUCTION (lambda-scan + Nc-cap merge)", prod);

            // --lts-report-nparts N: build a HYPOTHETICAL N-rank LTS-aware
            // partition of THIS (serial) mesh and report its load imbalance —
            // the Phase-5 realization-fraction predictor, measured on one core
            // WITHOUT launching N ranks.  Uses the PRODUCTION clustering (what a
            // real run steps) and the same METIS path as the run.
            const int report_nparts =
               GetIntArg(argc, argv, "--lts-report-nparts", 0);
            if (report_nparts > 1)
            {
               std::vector<int> xadj, adjncy;
               xadj.reserve(static_cast<std::size_t>(neL) + 1);
               xadj.push_back(0);
               const int *I = e2e.GetI();
               const int *J = e2e.GetJ();
               for (int e = 0; e < neL; ++e)
               {
                  for (int k = I[e]; k < I[e + 1]; ++k)
                  {
                     const int nb = J[k];
                     if (nb >= 0 && nb != e) { adjncy.push_back(nb); }
                  }
                  xadj.push_back(static_cast<int>(adjncy.size()));
               }
               std::vector<int> part_hyp;
               seas::LtsPartitionOptions popt;
               try
               {
                  seas::BuildLtsAwarePartition(neL, xadj, adjncy, prod.cluster,
                                               prod.num_clusters, report_nparts,
                                               nullptr, fault_pairs, popt,
                                               part_hyp);
                  const int viol = seas::VerifyLtsPartitionFaultLocality(
                     fault_pairs, part_hyp);
                  const auto pib = seas::ComputeLtsPartitionImbalance(
                     prod.cluster, prod.num_clusters, part_hyp, report_nparts,
                     nullptr, fault_pairs);
                  std::cout << "[lts-report] HYPOTHETICAL " << report_nparts
                            << "-rank LTS partition (production clustering):\n"
                            << "[lts-report]   imbalance overall=" << pib.overall
                            << "  worst-cluster c" << pib.worst_cluster << "="
                            << pib.worst_cluster_imbalance
                            << "  fault=" << pib.fault
                            << "  fault-lock-violations=" << viol << "\n"
                            << "[lts-report]   per-cluster imbalance:";
                  for (int c = 0; c < prod.num_clusters; ++c)
                  { std::cout << " c" << c << "=" << pib.per_cluster[c]; }
                  std::cout << "\n";
               }
               catch (const std::exception &ex)
               {
                  std::cout << "[lts-report] partition for " << report_nparts
                            << " ranks failed: " << ex.what() << "\n";
               }
            }
         }
      }
   }

   // -----------------------------------------------------------------
   // 16. Dry-run exit.
   // -----------------------------------------------------------------
   if (dry_run)
   {
      if (rank == 0)
      {
         std::cout << "[spatial_dyn] --dry-run: construction complete, "
                   << "exiting.\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 0;
   }

   // -----------------------------------------------------------------
   // 17. ParaView output wiring (parity Phases 1-6).
   //
   //     Parity Phase 3 master gate: if paraview_enabled == false AND
   //     every per-collection mode is "off", construct NOTHING.  The
   //     CLI `--paraview` / `--paraview-dt` overrides flip paraview_
   //     enabled on automatically (see merge block above).
   // -----------------------------------------------------------------
   // Parity Phase 4/5: wave-state component ordering assumptions.
   static_assert(SXX == 0 && SYY == 1 && SZZ == 2
                 && SXY == 3 && SYZ == 4 && SXZ == 5,
                 "Phase 4 sigma memcpy assumes wave_state.hpp ordering.");
   static_assert(VX == 6 && VY == 7 && VZ == 8,
                 "Phase 5 velocity memcpy assumes contiguous VX/VY/VZ; "
                 "wave_state.hpp ordering changed.");

   if (rank == 0)
   {
      std::filesystem::create_directories(cfg.output.output_dir);
   }
#ifdef MFEM_USE_MPI
   MPI_Barrier(comm);
#endif

   bool volume_pv_enabled = true;
   bool fault_pv_enabled  = true;
   bool bulk_pv_enabled   = true;
   auto volume_mode = ParseVolumeMode(cfg.output.paraview_volume,
                                      volume_pv_enabled);
   auto fault_mode  = ParseFaultMode(cfg.output.paraview_fault,
                                     fault_pv_enabled);
   {
      bool dummy;
      (void)ParseVolumeMode(cfg.output.paraview_bulk, dummy);
      bulk_pv_enabled = dummy;
   }

   // Parity Phase 1 §4 / Phase 3 §1: paraview_enabled is the HARD master
   // gate.  When false, force every per-collection enabled flag to false
   // regardless of the per-collection mode strings.  The CLI `--paraview`
   // merge above flips paraview_enabled to true when the user opts in.
   if (!cfg.output.paraview_enabled)
   {
      volume_pv_enabled = false;
      fault_pv_enabled  = false;
      bulk_pv_enabled   = false;
   }

   const bool any_pv_requested = volume_pv_enabled
                                || fault_pv_enabled
                                || bulk_pv_enabled;
   const bool primary_pv_active = any_pv_requested
                                  && (volume_pv_enabled || fault_pv_enabled);
   const bool bulk_pv_active    = any_pv_requested && bulk_pv_enabled;

   std::unique_ptr<seas::ParaViewOutput<ParMesh>> pv_out;
   std::unique_ptr<seas::ParaViewOutput<ParMesh>> pv_bulk_out;

   // Parity Phase 5: volume velocity + mpi_rank fields (lifetime of pv_out).
   std::unique_ptr<L2_FECollection> pv_vel_fec, pv_rank_fec;
   std::unique_ptr<ParFiniteElementSpace> pv_vel_fes, pv_rank_fes;
   std::unique_ptr<ParGridFunction>       pv_vel_gf, pv_rank_gf;

   // Parity Phase 4: secondary collection stress GFs (lifetime of pv_bulk_out).
   std::unique_ptr<L2_FECollection>       pv_bulk_sigma_fec;
   std::unique_ptr<ParFiniteElementSpace> pv_bulk_sigma_fes;
   std::unique_ptr<ParGridFunction>
      pv_bulk_sxx_gf, pv_bulk_syy_gf, pv_bulk_szz_gf,
      pv_bulk_sxy_gf, pv_bulk_syz_gf, pv_bulk_sxz_gf;

   // Parity Phase 6: SAFS fault static-parameter arrays (lifetime of pv_out).
   Vector pv_lsw_mu_s, pv_lsw_mu_d, pv_lsw_d_c;
   Vector pv_nuc_amplitude, pv_nuc_radial;
   Vector pv_sig_n_init, pv_tau1_init, pv_tau2_init;

   if (primary_pv_active)
   {
      pv_out = std::make_unique<seas::ParaViewOutput<ParMesh>>(
                  cfg.output.output_dir, pmesh, cfg.mesh.order,
                  "volume", volume_mode);
      pv_out->SetVolumePVDt(cfg.output.paraview_volume_dt);
      pv_out->SetVolumeSaveEnabled(volume_pv_enabled);
      if (cfg.output.paraview_every_steps > 0)
      {
         pv_out->output_every_n_steps = cfg.output.paraview_every_steps;
      }

      pv_out->fixed_dt = cfg.output.paraview_fault_dt;
      pv_out->GetSchedule().max_total_snapshots = cfg.output.max_snapshots;

      // Parity Phase 6: regime-adaptive cadence from TOML/CLI.
      pv_out->SetTotalRunTime(cfg.time.tfinal);
      if (cfg.output.paraview_coseismic_dt    > 0.0)
      { pv_out->GetSchedule().dt_coseismic    = cfg.output.paraview_coseismic_dt; }
      if (cfg.output.paraview_nucleation_dt   > 0.0)
      { pv_out->GetSchedule().dt_nucleation   = cfg.output.paraview_nucleation_dt; }
      if (cfg.output.paraview_interseismic_dt > 0.0)
      { pv_out->GetSchedule().dt_interseismic = cfg.output.paraview_interseismic_dt; }

#ifdef MFEM_USE_HDF5
      // Only set HDF5 chunk compression (and probe the H5Z-ZFP plugin) when the
      // volume collection is ACTUALLY being written.  ParseVolumeMode("off")
      // leaves volume_mode at the build default (Hdf5) while disabling saving,
      // so without the volume_pv_enabled guard a paraview_volume="off" config
      // that still carries a stale paraview_volume_zfp_tol would call
      // SetVolumeHDFCompression(ZfpAccuracy) -> ProbeH5ZZfpPluginOrAbort and
      // abort a run that emits no volume HDF5 at all.  Mirrors the
      // fault_pv_enabled guard on the fault path below.
      if (volume_pv_enabled
          && cfg.output.paraview_volume_zfp_tol > 0.0
          && volume_mode == ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5)
      {
         pv_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            cfg.output.paraview_volume_zfp_tol);
      }
      else if (volume_pv_enabled
               && cfg.output.paraview_volume_deflate_level >= 0
               && volume_mode == ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5)
      {
         pv_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::Deflate,
            static_cast<double>(cfg.output.paraview_volume_deflate_level));
      }
#endif

      pv_out->SetFaultOutputMode(fault_mode);
#ifdef MFEM_USE_HDF5
      if (fault_pv_enabled
          && fault_mode == ParaViewOutput<ParMesh>::FaultOutputMode::Hdf5
          && cfg.output.paraview_fault_zfp_tol > 0.0)
      {
         pv_out->SetFaultHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            cfg.output.paraview_fault_zfp_tol);
      }
#endif

      // R-003: opt back into legacy behaviour where the 8 L2-p0 fault
      // projection GFs are registered with the primary volume collection.
      // Without this the dynamic slip / slip_rate / traction / state
      // channels never reach `volume.vtkhdf` (only the static SAFS
      // params registered via `SetFaultParamsSpatial` do).  Must be
      // called BEFORE `InitFaultOutputBP5` (the registration happens
      // inside that call).
      pv_out->SetRegisterFaultProjectionsInVolumePV(true);

      pv_out->InitFaultOutputBP5(fault_int_faces, fault_shr_faces,
                                 nbf_per_face);

      // Parity Phase 5: register velocity (3-component L2 p=order) and
      // mpi_rank (L2 p=0).  Only when the volume mode is opted in.
      if (volume_pv_enabled)
      {
         pv_vel_fec = std::make_unique<L2_FECollection>(cfg.mesh.order, 3,
                                                        BasisType::GaussLobatto);
         pv_vel_fes = std::make_unique<ParFiniteElementSpace>(
                         &pmesh, pv_vel_fec.get(), 3, Ordering::byNODES);
         pv_vel_gf  = std::make_unique<ParGridFunction>(pv_vel_fes.get());
         *pv_vel_gf = 0.0;
         pv_out->RegisterDomainField("velocity", pv_vel_gf.get());

         pv_rank_fec = std::make_unique<L2_FECollection>(0, 3);
         pv_rank_fes = std::make_unique<ParFiniteElementSpace>(
                          &pmesh, pv_rank_fec.get());
         pv_rank_gf  = std::make_unique<ParGridFunction>(pv_rank_fes.get());
         *pv_rank_gf = static_cast<real_t>(rank);
         pv_out->RegisterDomainField("mpi_rank", pv_rank_gf.get());
      }

      // Parity Phase 6: SAFS fault static fields.
      if (num_fault_total > 0)
      {
         pv_lsw_mu_s     .SetSize(num_fault_total);
         pv_lsw_mu_d     .SetSize(num_fault_total);
         pv_lsw_d_c      .SetSize(num_fault_total);
         pv_nuc_amplitude.SetSize(num_fault_total);
         pv_nuc_radial   .SetSize(num_fault_total);
         pv_sig_n_init   .SetSize(num_fault_total);
         pv_tau1_init    .SetSize(num_fault_total);
         pv_tau2_init    .SetSize(num_fault_total);
         for (int i = 0; i < num_fault_total; ++i)
         {
            const DOFData &d = dof_data[i];
            pv_lsw_mu_s(i) = d.lsw_mu_s;
            pv_lsw_mu_d(i) = d.lsw_mu_d;
            pv_lsw_d_c (i) = d.lsw_d_c;
            const real_t ad = (nuc_params.amplitude_dip.Size()    == num_fault_total)
                              ? nuc_params.amplitude_dip(i)    : 0.0;
            const real_t as = (nuc_params.amplitude_strike.Size() == num_fault_total)
                              ? nuc_params.amplitude_strike(i) : 0.0;
            pv_nuc_amplitude(i) = std::sqrt(ad*ad + as*as);
            pv_nuc_radial   (i) = (nuc_params.radial.Size() == num_fault_total)
                                  ? nuc_params.radial(i)   : 0.0;
            pv_sig_n_init   (i) = d.sigma_n_corr;
            pv_tau1_init    (i) = d.tau1_corr;
            pv_tau2_init    (i) = d.tau2_corr;
         }
      }
      // Parity Phase 6: publish all 8 SAFS static fields via
      // SetFaultParamsSpatial on seas::ParaViewOutput.  Implementation
      // registers each field with `pv_dc_` (the PRIMARY volume PV
      // collection), so these fields land in `volume.vtkhdf` under
      // their SAFS-semantic names (NOT TPV104's a/Dc/x2/x3).  The
      // companion fault artefact (`fault.vtkhdf`) is populated
      // separately by `WriteFaultSurfaceVTU` (R-001) and carries the
      // dynamic slip / slip_rate / traction / state channels.
      // SetFaultParamsSpatial is safe to call on a rank with zero
      // local fault DOFs.
      pv_out->SetFaultParamsSpatial({
         {"lsw_mu_s",          &pv_lsw_mu_s},
         {"lsw_mu_d",          &pv_lsw_mu_d},
         {"lsw_d_c",           &pv_lsw_d_c},
         {"nuc_amplitude",     &pv_nuc_amplitude},
         {"nuc_radial_factor", &pv_nuc_radial},
         {"sigma_n_init",      &pv_sig_n_init},
         {"tau1_init",         &pv_tau1_init},
         {"tau2_init",         &pv_tau2_init},
      });
   }

   // Parity Phase 4: secondary `pv_bulk_out` collection (6 stress
   // components → ParaView_bulk/stress.vtkhdf).
   if (bulk_pv_active)
   {
      const std::string bulk_dir = cfg.output.output_dir + "/ParaView_bulk";
      if (rank == 0)
      {
         std::filesystem::create_directories(bulk_dir);
      }
#ifdef MFEM_USE_MPI
      MPI_Barrier(comm);
#endif
      bool bulk_mode_enabled = bulk_pv_enabled;
      auto bulk_mode = ParseVolumeMode(cfg.output.paraview_bulk,
                                       bulk_mode_enabled);
      pv_bulk_out = std::make_unique<seas::ParaViewOutput<ParMesh>>(
                       bulk_dir, pmesh, cfg.mesh.order,
                       "stress", bulk_mode);
      pv_bulk_out->fixed_dt = cfg.output.paraview_bulk_dt;

      pv_bulk_sigma_fec = std::make_unique<L2_FECollection>(
                            cfg.mesh.order, 3, BasisType::GaussLobatto);
      pv_bulk_sigma_fes = std::make_unique<ParFiniteElementSpace>(
                            &pmesh, pv_bulk_sigma_fec.get());
      pv_bulk_sxx_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      pv_bulk_syy_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      pv_bulk_szz_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      pv_bulk_sxy_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      pv_bulk_syz_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      pv_bulk_sxz_gf = std::make_unique<ParGridFunction>(pv_bulk_sigma_fes.get());
      *pv_bulk_sxx_gf = 0.0; *pv_bulk_syy_gf = 0.0; *pv_bulk_szz_gf = 0.0;
      *pv_bulk_sxy_gf = 0.0; *pv_bulk_syz_gf = 0.0; *pv_bulk_sxz_gf = 0.0;
      pv_bulk_out->RegisterDomainField("sigma_xx", pv_bulk_sxx_gf.get());
      pv_bulk_out->RegisterDomainField("sigma_yy", pv_bulk_syy_gf.get());
      pv_bulk_out->RegisterDomainField("sigma_zz", pv_bulk_szz_gf.get());
      pv_bulk_out->RegisterDomainField("sigma_xy", pv_bulk_sxy_gf.get());
      pv_bulk_out->RegisterDomainField("sigma_yz", pv_bulk_syz_gf.get());
      pv_bulk_out->RegisterDomainField("sigma_xz", pv_bulk_sxz_gf.get());

#ifdef MFEM_USE_HDF5
      if (cfg.output.paraview_bulk_zfp_tol > 0.0
          && bulk_mode == ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5)
      {
         pv_bulk_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            cfg.output.paraview_bulk_zfp_tol);
      }
      else if (cfg.output.paraview_bulk_deflate_level >= 0
               && bulk_mode == ParaViewOutput<ParMesh>::VolumeOutputMode::Hdf5)
      {
         pv_bulk_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::Deflate,
            static_cast<double>(cfg.output.paraview_bulk_deflate_level));
      }
#endif
   }

   // Phase 3: free-surface slice writer (default-ON; independent of the
   // any_pv_requested master gate).  Constructed HERE — after the volume/bulk
   // collections, BEFORE the banner (which reports it) — because the ctor needs
   // only pmesh / fs_attrs / order / rank / dt, NOT ndof_total (R-004).
   bool fs_enabled = false;
   auto fs_mode = ParseFreeSurfaceMode(cfg.output.paraview_free_surface, fs_enabled);
   std::unique_ptr<seas::FreeSurfaceOutput> fs_out;
   if (fs_enabled)
   {
      if (fs_attrs.Size() == 0)
      {
         if (rank == 0)
         {
            std::cerr << "[free-surface] WARNING: paraview_free_surface on but "
                         "no free-surface attributes ([boundary].natural_attrs "
                         "/ [output].paraview_free_surface_attrs both empty); "
                         "slice disabled.\n";
         }
      }
      else
      {
         const std::string fs_dir =
            cfg.output.output_dir + "/ParaView_free_surface";
         if (rank == 0) { std::filesystem::create_directories(fs_dir); }
#ifdef MFEM_USE_MPI
         MPI_Barrier(comm);
#endif
         // Collective: make_unique<FreeSurfaceOutput> calls CreateFromBoundary,
         // which is MPI-collective — ALL ranks must reach it.  The guards
         // (fs_enabled, fs_attrs.Size()) are rank-invariant, so they do.
         fs_out = std::make_unique<seas::FreeSurfaceOutput>(
                     fs_dir, pmesh, fs_attrs, cfg.mesh.order, rank,
                     cfg.output.paraview_free_surface_dt, "free_surface",
                     fs_mode);
         if (fs_out->GlobalNE() == 0)
         {
            if (rank == 0)
            {
               std::cerr << "[free-surface] WARNING: 0 surface elements matched "
                            "the resolved free-surface attributes; slice "
                            "DISABLED. Check [boundary].natural_attrs / "
                            "[output].paraview_free_surface_attrs.\n";
            }
            // R-303: GlobalNE() is MPI_Allreduce'd, so "== 0" is identical on
            // every rank; all ranks drop fs_out together (collective-safe — no
            // half-disabled state, no per-step empty Save for the whole run).
            fs_out.reset();
         }
      }
   }

   // Parity Phase 3: rank-0 banner.
   if (rank == 0)
   {
      if (any_pv_requested)
      {
         std::cout << "ParaView output: ON (prefix=" << cfg.output.output_dir
                   << ")\n"
                   << "  Volume:   " << cfg.output.paraview_volume
                   << " (dt=" << cfg.output.paraview_volume_dt << " s)\n"
                   << "  Bulk:     " << cfg.output.paraview_bulk
                   << " (dt=" << cfg.output.paraview_bulk_dt   << " s)\n"
                   << "  Fault:    " << cfg.output.paraview_fault
                   << " (dt=" << cfg.output.paraview_fault_dt  << " s)\n";
         if (cfg.output.paraview_volume_zfp_tol > 0.0)
         {
            std::cout << "  Volume ZFP tol: "
                      << cfg.output.paraview_volume_zfp_tol << "\n";
         }
         if (cfg.output.paraview_bulk_zfp_tol > 0.0)
         {
            std::cout << "  Bulk   ZFP tol: "
                      << cfg.output.paraview_bulk_zfp_tol   << "\n";
         }
         if (cfg.output.paraview_fault_zfp_tol > 0.0)
         {
            std::cout << "  Fault  ZFP tol: "
                      << cfg.output.paraview_fault_zfp_tol  << "\n";
         }
         if (cfg.output.paraview_coseismic_dt > 0.0
             || cfg.output.paraview_nucleation_dt > 0.0
             || cfg.output.paraview_interseismic_dt > 0.0)
         {
            std::cout << "  Regime-adaptive cadence active.\n";
         }
         if (cfg.output.paraview_enabled
             && !primary_pv_active && !bulk_pv_active)
         {
            std::cout << "  WARNING: paraview_enabled = true but every "
                      << "per-collection mode is 'off' — no ParaView files "
                      << "will be written.\n";
         }
      }
      else
      {
         std::cout << "ParaView output: OFF\n";
      }
      // Phase 3: the free-surface slice is independent of the master gate, so
      // report it AFTER the any_pv_requested if/else — even when volume/bulk/
      // fault are all off (any_pv_requested == false), the slice still writes
      // and must be banner-visible (R-006).
      if (fs_out)
      {
         // R-302: echo the RESOLVED attributes so silent over-inclusion (e.g. a
         // bottom zero-traction face tagged natural) is visible in the banner.
         std::cout << "  FreeSurf: " << cfg.output.paraview_free_surface
                   << " (dt=" << cfg.output.paraview_free_surface_dt << " s, "
                   << fs_out->GlobalNE() << " surf elems, attrs {";
         for (int i = 0; i < fs_attrs.Size(); ++i)
         {
            std::cout << (i ? "," : "") << fs_attrs[i];
         }
         std::cout << "})\n";
      }
   }

   // ParaView GFs are allocated for the time loop only; sized after
   // pv_out is up.  No-op when num_fault_total = 0 on this rank.
   Vector pv_local_slip, pv_local_slip_rate, pv_local_traction;
   Vector pv_local_state, pv_local_normal_stress;
   // R-001: scratch buffers consumed by WriteFaultSurfaceVTU (ADER
   // one-shot: _k4 copies of the live values so the VTU delta channel
   // reads zero, mirroring drivers/tpv205_driver.cpp:2252-2259).  The
   // four `local_a / local_Dc / local_x2 / local_x3` slots in the
   // TPV-shaped signature stay empty — SAFS-semantic static params are
   // already registered with `pv_dc_` via `SetFaultParamsSpatial` and
   // appear in `volume.vtkhdf`.
   Vector pv_local_slip_rate_k4, pv_local_traction_k4;
   Vector pv_local_normal_stress_k4;
   const Vector pv_empty_param;  // unused TPV-shaped param slots
   if (pv_out)
   {
      pv_local_slip.SetSize(2 * num_fault_total);
      pv_local_slip_rate.SetSize(2 * num_fault_total);
      pv_local_traction.SetSize(2 * num_fault_total);
      pv_local_state.SetSize(num_fault_total);
      pv_local_normal_stress.SetSize(num_fault_total);
      pv_local_slip_rate_k4.SetSize(2 * num_fault_total);
      pv_local_traction_k4.SetSize(2 * num_fault_total);
      pv_local_normal_stress_k4.SetSize(num_fault_total);
   }

   // -----------------------------------------------------------------
   // 18. Restart (R-105 corrected schema: per-rank Q + dof_data).
   // -----------------------------------------------------------------
   const int ndof_total = wave.GetScalarNDof();
   Vector Q(NUM_STATE * ndof_total);
   Q = 0.0;

   real_t t      = cfg.time.t_initial;
   int    step0  = 0;
   real_t dt_now = dt;
   // (Checkpoint V2) Recompute the current clustering's layout hash — the same
   // (rate, Nc, cluster ids, dt_base=lambda*dt_cfl, lambda) the LTS loop uses —
   // for the restart layout-match check.  Deterministic; the LTS loop rebuilds
   // it identically.  Only meaningful when lts is enabled (fault-free bulk).
   auto compute_lts_layout_hash = [&]() -> std::uint64_t
   {
      // Reuse the pre-operator clustering (single-sourced with the reorder), so
      // the restart hash matches the write-time hash exactly.
      MFEM_VERIFY(lts_layout_ready,
                  "spatial_dyn: LTS layout hash requested but the pre-operator "
                  "clustering was not built (internal error).");
      return LtsLayoutHash(/*rate=*/2, lts_cl.num_clusters, lts_cl.cluster,
                           lts_cl.lambda * dt_cfl, lts_cl.lambda);
   };

   // (P-006) Per-QP canonical permutation of the (reordered) fault stack, used
   // to serialize checkpoint dof_data in CANONICAL (layout-independent) order.
   // Empty (=> nullptr) when the reorder is off, so the write/read are unchanged.
   const std::vector<int> qp_canon_perm =
      wave.FaultFacesReordered() ? wave.GetFaultQpCanonicalPerm()
                                 : std::vector<int>{};
   const std::vector<int> *qp_canon_ptr =
      qp_canon_perm.empty() ? nullptr : &qp_canon_perm;

   // (P-006 canary) A GTS-path checkpoint under the fault-face reorder is written
   // as V2 (layout hash + canonical perm) so the on-disk dof_data stays
   // layout-independent and LtsCheckpointCheck refuses a cross-mode (lts on<->off)
   // restart; without the reorder it is the usual V1 write (byte-identical).
   auto write_gts_checkpoint =
      [&](const std::string &prefix, real_t tt, real_t dtt, int stp)
   {
      if (wave.FaultFacesReordered())
      {
         mfem::seas::internal::WriteTpv104CheckpointV2Impl(
            prefix, tt, dtt, stp, /*lts_mode=*/1, compute_lts_layout_hash(),
            Q, dof_data, rank, nprocs, "spatial_dyn", qp_canon_ptr);
      }
      else
      {
         WriteTpv104Checkpoint(prefix, tt, dtt, stp, Q, dof_data, rank, nprocs
#ifdef MFEM_USE_MPI
                               , comm
#endif
                               , "spatial_dyn");
      }
   };

   if (!restart_prefix.empty())
   {
      std::string driver_tag;
      const int ckpt_version = mfem::seas::internal::PeekTpv104CheckpointVersion(restart_prefix, rank);
      if (ckpt_version == 2)
      {
         // V2 (LTS) checkpoint: refuse under lts=off; else validate the hash.
         MFEM_VERIFY(cfg.numerics.LtsEnabled(),
                     "spatial_dyn: refusing restart — a V2 (LTS) checkpoint "
                     "cannot resume an lts=\"off\" (GTS) run (different scheme).");
         int lts_mode = 0; std::uint64_t stored_hash = 0;
         const bool ok = mfem::seas::internal::ReadTpv104CheckpointV2Impl(
            restart_prefix, t, dt_now, step0, lts_mode, stored_hash, Q,
            NUM_STATE * ndof_total, dof_data, rank, nprocs, &driver_tag,
            qp_canon_ptr);
         MFEM_VERIFY(ok, "spatial_dyn: V2 checkpoint read failed for rank " << rank);
         const std::uint64_t computed_hash = compute_lts_layout_hash();
         const LtsCheckpointDecision decision =
            LtsCheckpointCheck(2, true, stored_hash, computed_hash);
         MFEM_VERIFY(decision == LtsCheckpointDecision::Accept,
                     "spatial_dyn: refusing V2 restart — layout hash mismatch "
                     "(stored 0x" << std::hex << stored_hash << " != computed 0x"
                     << computed_hash << std::dec << "); mesh/clustering/config "
                     "changed since the checkpoint was written.");
      }
      else
      {
      // (P-006) A V1 checkpoint holds CANONICAL (non-reordered) dof_data; it
      // cannot be read directly into a run whose fault stack was reordered
      // (lts != "off").  Refuse loudly rather than silently mis-scattering.
      // Restart the canary from its own V2 checkpoint, or start it fresh.
      MFEM_VERIFY(!qp_canon_ptr,
                  "spatial_dyn: refusing restart — a V1 (non-reordered) "
                  "checkpoint cannot resume a run with the fault-face reorder "
                  "active (lts != \"off\"); its dof_data is in canonical order. "
                  "Restart from a V2 checkpoint or start the run fresh.");
      const bool ok = ReadTpv104Checkpoint(
         restart_prefix, t, dt_now, step0, Q, NUM_STATE * ndof_total,
         dof_data, rank, nprocs
#ifdef MFEM_USE_MPI
         , comm
#endif
         , &driver_tag);
      if (!ok)
      {
         if (rank == 0)
         {
            std::cerr << "ERROR: restart prefix '" << restart_prefix
                      << "' has no per-rank file for rank " << rank
                      << ".\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 5;
      }
      // (Checkpoint V2) An lts != off run must NOT resume a V1 (GTS) checkpoint.
      if (cfg.numerics.LtsEnabled())
      {
         const LtsCheckpointDecision decision =
            LtsCheckpointCheck(/*file_version=*/1, /*lts_enabled=*/true, 0, 0);
         MFEM_VERIFY(decision == LtsCheckpointDecision::Accept,
                     "spatial_dyn: refusing restart — lts=\"" << cfg.numerics.lts
                     << "\" cannot resume a V1 (GTS) checkpoint (different "
                     "trajectory).  LTS runs write V2 checkpoints.");
      }
      }   // end else (V1 path)
      if (!driver_tag.empty() && driver_tag != "spatial_dyn")
      {
         if (rank == 0)
         {
            std::cerr << "ERROR: refusing restart: checkpoint driver_tag '"
                      << driver_tag << "' != 'spatial_dyn'.\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 5;
      }
      if (driver_tag.empty() && rank == 0)
      {
         std::cout << "[restart] WARNING: restoring pre-tag checkpoint "
                   << "into spatial_dyn — proceed only if you know the "
                   << "checkpoint is compatible.\n";
      }
      // SetTime so the LSW_ForcedRupture guard sees a fresh t.
      wave.SetTime(t);
   }

   // Phase 7: seed the one-shot instantaneous-overstress patch ONCE, after the
   // fault DOFs are initialized.  No-op for the per-sub-step kinds (Gaussian /
   // compact-circular) and for StaticOverstress.  Skipped on restart — the
   // checkpoint already restored the seeded tau2_nuc into dof_data, so
   // re-seeding would double-count the patch.
   if (restart_prefix.empty())
   {
      nuc->ApplyOnce(dof_data);
   }

   // -----------------------------------------------------------------
   // 19. Sub-step iterator (Phase 2/3): dispatch through the
   //     IFrictionIterator strategy.  The LswFrictionIterator is a
   //     transparent forwarder over the same Tpv205SubStepIterator, so an
   //     LSW run is byte-identical to the pre-Phase-2 driver.
   // -----------------------------------------------------------------
   // Phase 3: dispatch the friction iterator through the factory.  LSW ->
   // LswFrictionIterator (byte-identical forwarder); rate_state -> the
   // aging RateStateAgingFrictionIterator.  rs is passed only for the RS
   // path (reserved for future per-DOF wiring; the aging adapter reads the
   // scalar RateStateBlock).
   auto friction_iterator =
      MakeFrictionIterator(cfg, fault_flux, is_lsw ? nullptr : &rs);
   IFrictionIterator &substep_iterator = *friction_iterator;
   // R-023: keep IFrictionIterator::WaveOpLaw() live on the production path as
   // a single-source-of-truth cross-check — the wave-op law was set from is_lsw
   // (:975) ~900 lines earlier (before the iterator exists), so assert the
   // iterator the factory built agrees with it rather than leaving WaveOpLaw()
   // dead surface area.
   MFEM_VERIFY(substep_iterator.WaveOpLaw() ==
               (is_lsw ? FaultFrictionLaw::LSW : FaultFrictionLaw::RateAndState),
               "spatial_dyn_driver: friction iterator WaveOpLaw() ("
               << static_cast<int>(substep_iterator.WaveOpLaw())
               << ") disagrees with the wave-operator fault law set from "
                  "is_lsw — friction dispatch is inconsistent.");
   {
      const int O = std::max(1, cfg.numerics.ader_order);
      std::vector<real_t> deltaT(O, dt / static_cast<real_t>(O));
      std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
      substep_iterator.SetSubSteps(deltaT, weights);
   }
   // PLAN_speckle_slip_runaway (diag): tell the iterator how many fault QPs
   // are LOCAL/interior so the env-gated [SLIP] trace can tag each QP shared
   // vs interior (dof_data is laid out interior [0,n_local) then shared).
   substep_iterator.SetDiagNumLocalFaultQPs(wave.GetNumLocalFaultQPs());

   // Phase 3: wire the secular slip/state resample onto the iterator.  Both laws
   // are now supported (R-001 fix, 2026-06-04):
   //   - rate-state (TPV102/104): resamples the per-macro-step Δψ increment.
   //   - LSW (TPV205/TPV31): resamples the per-macro-step accumulated-slip-
   //     MAGNITUDE increment (path length Σ|V|·dt, plan §4.2/§6 Phase 3) and
   //     rescales the directional slip (slip1,slip2) to the dealiased magnitude,
   //     preserving direction.  τ_corr / the slip-rate stay from the un-resampled
   //     friction solve (plan §6 Phase 3 — "do not rebuild τ_corr from a
   //     resampled quantity").
   // R-002: the resample is the SECULAR layer applied ON TOP OF over-integration
   // (plan §4.2: "resample is a no-op without over-integration; over-integration
   // is the prerequisite").  Gate it on over-integration being ON, NOT on the
   // per-face head-count: on a TRIANGLE face the minimal 2*order rule is
   // over-determined for p>=3 (#QP > #DOF ⇒ R != I), so keying off (#QP != #DOF)
   // would make --fault-resample ALONE act at p>=3, violating the §Acceptance
   // "resample alone ⇒ byte-exact" contract.  With over-integration off,
   // resample_active is false at EVERY order ⇒ a TRUE no-op (byte-exact).
   const bool resample_active = cli_fault_resample && (cli_fault_overint > 0);

   // R-004 (plan §4.2 / §Phase 4 step 2): the resample must be OFF across a
   // genuine fault-normal material contrast (SeisSol's BiMaterialFault —
   // "resampling introduces artificial oscillations").  The spatial setup sets
   // Zp_plus==Zp_minus / Zs_plus==Zs_minus per DOF (the matrix path is
   // depth-heterogeneity, equal across the fault at each QP), so this never
   // trips today; it guards a future genuine-contrast configuration.
   if (resample_active)
   {
      int local_contrast = 0;
      for (const DOFData &d : dof_data)
      {
         if (d.Zp_plus != d.Zp_minus || d.Zs_plus != d.Zs_minus)
         { local_contrast = 1; break; }
      }
      int global_contrast = local_contrast;
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&local_contrast, &global_contrast, 1, MPI_INT, MPI_MAX, comm);
#endif
      MFEM_VERIFY(!global_contrast,
                  "--fault-resample is invalid across a genuine fault-normal "
                  "material contrast (Zp_plus != Zp_minus or Zs_plus != Zs_minus): "
                  "resampling introduces artificial oscillations (SeisSol "
                  "BiMaterialFault, plan §4.2).  Disable --fault-resample for "
                  "that configuration.");
   }

   substep_iterator.SetFaultResample(
      resample_active ? &fault_resample_R : nullptr,
      nbf_per_face, resample_active);

   // Phase 7: per-sub-step nucleation hook, dispatched through the
   // INucleationMethod strategy.  Fires once per ADER sub-step BEFORE the
   // per-QP friction solve so the perturbed tau{1,2}_nuc is visible to
   // s.tau{1,2}_total.  For the Gaussian path ApplyIncrement is the same
   // ApplyGradualOverstressIncrement call as the pre-Phase-7 lambda (the
   // cfg.nucleation.enabled guard is subsumed: a disabled config yields a
   // StaticOverstress whose ApplyIncrement is a no-op).
   // R-005: gate on IsPerSubStep() so the strategy's per-sub-step intent is
   // honoured explicitly.  For kinds that do no per-sub-step work (Static,
   // InstantaneousOverstressCircular) we pass an explicit no-op rather than
   // calling their (already no-op) ApplyIncrement every sub-step.  Note the
   // iterator REJECTS an empty std::function (friction_substep_iterator.hpp:108),
   // so the opt-out must be a no-op lambda, not a default-constructed function.
   std::function<void(real_t, real_t)> nuc_cb;
   if (nuc->IsPerSubStep())
   {
      nuc_cb = [&nuc, &dof_data](real_t t_sub_end, real_t dt_sub)
      {
         nuc->ApplyIncrement(dof_data, t_sub_end, dt_sub);
      };
   }
   else
   {
      nuc_cb = [](real_t, real_t) {};
   }

   // Phase 14: select the RK Butcher tableau ONCE (only consumed on the RK
   // branch of the time loop).  ADER ignores it.  Both laws use the SAME
   // tableaus; only the rate_state RK path reads the [friction.rate_state]
   // block (for the PsiRate ψ coupling), so that requirement is gated under
   // !is_lsw — the LSW RK path (AdvanceRKCoupledLSW_Spatial) integrates slip,
   // not ψ, and needs no rate_state block.
   RKTableau rk_tab;
   if (is_rk)
   {
      if (!is_lsw)
      {
         MFEM_VERIFY(cfg.rate_state.has_value(),
                     "spatial_dyn_driver: --time-integrator rk4|rk45 with "
                     "[meta].law=\"rate_state\" requires a [friction.rate_state] "
                     "block (the RK ψ coupling reads it).");
      }
      rk_tab = (cfg.numerics.time_integrator
                == spatial::TimeIntegratorKind::RK45)
               ? MakeDormandPrinceRK45Tableau()
               : MakeRK4Tableau();
      ValidateTableau(rk_tab);
      if (rank == 0)
      {
         std::cout << "[time-integrator] RK tableau = " << rk_tab.name
                   << " (" << rk_tab.stages << " stages)"
                   << (is_lsw ? "  [LSW coupled-RK on (Q, slip)]"
                              : "  [rate-state coupled-RK on (Q, ψ, slip)]")
                   << "\n";
      }
   }

   // ParaView snapshot writer (Parity Phases 1-6).  Updates the 5 BP5
   // fault projection GFs + writes the primary / bulk collections.
   // Splits writes into "fault collection wants" (primary pv_out) and
   // "bulk wants" (pv_bulk_out) so each can fire on its own cadence.
   auto paraview_write = [&](int step_num, real_t time, real_t V_max)
   {
      const bool fault_wants = pv_out
                              && pv_out->PeekShouldWrite(step_num, time, V_max);
      const bool bulk_wants  = pv_bulk_out
                              && pv_bulk_out->PeekShouldWrite(step_num, time, V_max);
      const bool fs_wants    = fs_out && fs_out->ShouldWrite(time);
      if (!fault_wants && !bulk_wants && !fs_wants) { return; }

      // Free surface: independent of fault/bulk AND of the
      // `if (!fault_wants) return;` guard below — must be written HERE, before
      // the bulk block and that guard, NOT at the end of the lambda (R-002).
      // Reuses the same Q velocity block the volume path memcpys below.
      if (fs_wants)
      {
         fs_out->UpdateVelocity(Q.GetData() + VX * ndof_total, ndof_total);
         fs_out->Save(step_num, time);
      }

      // Parity Phase 4: publish 6 stress components to pv_bulk_out.
      if (bulk_wants)
      {
         const real_t* Q_data = Q.GetData();
         std::memcpy(pv_bulk_sxx_gf->GetData(), Q_data + SXX * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_syy_gf->GetData(), Q_data + SYY * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_szz_gf->GetData(), Q_data + SZZ * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_sxy_gf->GetData(), Q_data + SXY * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_syz_gf->GetData(), Q_data + SYZ * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_sxz_gf->GetData(), Q_data + SXZ * ndof_total,
                     ndof_total * sizeof(real_t));
         pv_bulk_out->ForceSave(step_num, time);
         pv_bulk_out->CommitSchedule(time, V_max);
      }

      if (!fault_wants) { return; }

      for (int i = 0; i < num_fault_total; ++i)
      {
         const DOFData &d = dof_data[i];
         pv_local_slip(2 * i + 0)      = d.slip1;
         pv_local_slip(2 * i + 1)      = d.slip2;
         pv_local_slip_rate(2 * i + 0) = d.V1;
         pv_local_slip_rate(2 * i + 1) = d.V2;
         pv_local_traction(2 * i + 0)  = d.tau1_corr;
         pv_local_traction(2 * i + 1)  = d.tau2_corr;
         // R-005/R-024: LSW writes the friction coefficient; RS writes the
         // state variable psi.  InitializeFaultDOFs_Spatial_RS never sets the
         // lsw_* fields (they stay 0), so the LSW formula would emit 0/NaN for
         // an RS run.  The selection lives in the shared FaultStateChannelValue
         // (dynamic/fault_state_channel.hpp) so this writer and the R-005 test
         // exercise the SAME code path.
         pv_local_state(i) = mfem::seas::FaultStateChannelValue(is_lsw, d);
         (void)time;
         pv_local_normal_stress(i)     = d.sigma_n_corr;

         // R-001: ADER one-shot — no stage-4 distinct from averaged
         // DOFData (mirror of tpv205_driver.cpp:2252-2259).  Writing
         // the same values into the `_k4` buffers makes the VTU delta
         // channel read exactly zero.
         pv_local_slip_rate_k4(2 * i + 0) = d.V1;
         pv_local_slip_rate_k4(2 * i + 1) = d.V2;
         pv_local_traction_k4(2 * i + 0)  = d.tau1_corr;
         pv_local_traction_k4(2 * i + 1)  = d.tau2_corr;
         pv_local_normal_stress_k4(i)     = d.sigma_n_corr;
      }

      // Parity Phase 5: refresh the volume velocity field directly from
      // Q (byNODES so VX/VY/VZ are contiguous 3*ndof block).
      if (pv_out->GetVolumeSaveEnabled() && pv_vel_gf)
      {
         std::memcpy(pv_vel_gf->GetData(),
                     Q.GetData() + VX * ndof_total,
                     3 * ndof_total * sizeof(real_t));
      }

      // R-002: refresh fault projections unconditionally — they feed
      // both the volume PV (when registered via R-003) AND the fault
      // VTU/VTKHDF writer (R-001).  Only the volume save itself
      // (`ForceSave` → `pv_dc_->Save()`) is gated on the user toggle.
      pv_out->UpdateFaultFieldsBP5(pv_local_slip, pv_local_slip_rate,
                                   pv_local_traction, pv_local_state,
                                   pv_local_normal_stress);
      if (pv_out->GetVolumeSaveEnabled())
      {
         pv_out->ForceSave(step_num, time);
      }
      pv_out->CommitSchedule(time, V_max);

      // R-001: write the fault-surface VTU/VTKHDF artefact (the file
      // `verify_spatial_dyn_smoke_safs.py` opens and the plan
      // §Acceptance Criteria L491 requires).  Independent of volume
      // save — fault is its own collection.  The four TPV-shaped
      // param slots (`local_a / local_Dc / local_x2 / local_x3`) stay
      // empty: SAFS-semantic static params are already published to
      // `volume.vtkhdf` via `SetFaultParamsSpatial`.
      pv_out->WriteFaultSurfaceVTU(
         cfg.output.output_dir, step_num, time, rank, nprocs,
         pv_local_slip, pv_local_slip_rate, pv_local_traction,
         pv_local_state, pv_local_normal_stress,
         pv_empty_param, pv_empty_param, pv_empty_param, pv_empty_param,
         pv_local_slip_rate_k4, pv_local_traction_k4,
         pv_local_normal_stress_k4);
   };

   if (restart_prefix.empty())
   {
      paraview_write(0, cfg.time.t_initial, 0.0);
   }

   // SCEC on-fault station traces (one `.dat` per station) — the verification
   // artifact consumed by tpv{31,102,104,205}/visualize_results.py.  The
   // spatial driver otherwise emits only ParaView/HDF5 fault output.  Each
   // problem reuses the SAME benchmark-format writer the native driver uses
   // (extracted to the lean dynamic/tpv*_stations.hpp), selected by the
   // [meta].tag, so the columns are byte-for-byte the SCEC layout:
   //   tpv31  -> TPV31StationWriter   (LSW μ_eff, 30-station grid)
   //   tpv102 -> TPV102StationWriter  (rate-state, log10_theta)
   //   tpv104 -> TPV104StationWriter  (rate-state, psi)
   //   tpv205 -> TPV205StationWriter  (LSW μ_eff, 16-station grid)
   // `fault_coords` is sized num_fault_total with the same (interior-then-
   // shared) DOF ordering as `dof_data`; each writer resolves cross-rank
   // station ownership internally.  The four writers are distinct types that
   // share one interface (Open(...,comm) / WriteStep(t,dof_data) / Close()),
   // so WriteStep/Close are type-erased into std::function closures over a
   // shared owner — at most one is wired (no [meta].tag match -> none, e.g.
   // SAFS, which has no SCEC station grid).
   std::shared_ptr<void> stations_owner;
   std::function<void(real_t, const std::vector<DOFData> &)> stations_write;
   std::function<void()> stations_close;
   {
      auto wire_stations = [&](auto *writer_tag, const auto &stations)
      {
         using WriterT = std::remove_pointer_t<decltype(writer_tag)>;
         auto w = std::make_shared<WriterT>();
#ifdef MFEM_USE_MPI
         w->Open(cfg.output.output_dir, cfg.problem.tag, stations,
                 fault_coords, num_fault_total, comm);
#else
         w->Open(cfg.output.output_dir, cfg.problem.tag, stations,
                 fault_coords, num_fault_total);
#endif
         if (restart_prefix.empty())
         {
            w->WriteStep(cfg.time.t_initial, dof_data);
         }
         stations_owner = w;
         stations_write = [w](real_t tt, const std::vector<DOFData> &dd)
         { w->WriteStep(tt, dd); };
         stations_close = [w]() { w->Close(); };
      };

      const std::string &tag = cfg.problem.tag;
      if      (tag == "tpv31")  { wire_stations((TPV31StationWriter  *)nullptr, DefaultStations_TPV31()); }
      else if (tag == "tpv102") { wire_stations((TPV102StationWriter *)nullptr, DefaultStations()); }
      else if (tag == "tpv104") { wire_stations((TPV104StationWriter *)nullptr, DefaultStations_TPV104()); }
      else if (tag == "tpv205") { wire_stations((TPV205StationWriter *)nullptr, DefaultStations_TPV205()); }
      else if (tag == "tpv6")   { wire_stations((TPV6StationWriter   *)nullptr, DefaultStations_TPV6()); }
   }

   // -----------------------------------------------------------------
   // 20. Time loop.
   // -----------------------------------------------------------------
   Vector Q_new(Q.Size());
   real_t V_max_global = 0.0;
   // R-603 round-6: track the last completed step so the final
   // checkpoint at L1262 records the actual step the loop reached,
   // not `nsteps` unconditionally.
   int last_completed_step = step0;

   // -----------------------------------------------------------------
   // 20a. (LTS Phase 2) Clustered sync-interval loop — FAULT-FREE bulk only.
   //   Entirely gated on lts != "off", so lts="off" runs the unchanged GTS
   //   per-step loop below (byte-identical).  The production SAFS driver
   //   interleaves fault sub-stepping with the bulk corrector every step, so
   //   fault-LTS is Phase 3; this path REQUIRES a fault-free config and np=1
   //   (the serial-mesh clustering for np>1 determinism is Phase 1b).  The
   //   stepping engine (predictor/corrector/multi-rate coupling) is unit-
   //   validated (test_lts_predictor / test_lts_layout); this is its wiring.
   // -----------------------------------------------------------------
   // (Phase 3 canary) A fault run under lts != "off" REORDERS the fault-face
   // list into cluster-contiguous order at operator construction (P-006, done
   // above via lts_cluster_ptr) but STEPS GTS — the fault-half LTS interleave is
   // Phase 3.  So the LTS multi-cluster sync loop runs ONLY for the fault-free
   // bulk problem (Phase 2); a fault + lts run falls through to the GTS loop
   // below with the reorder ACTIVE (the still-GTS reorder canary, P-006 gate).
   // LTS Phase 4a: the fault-free bulk stepper now runs at np>1 (the MPI
   // ghost-field seam exchange).  The gate MUST use the GLOBAL fault count
   // (`num_fault_global`), not the per-rank `num_fault_total`: a fault mesh
   // partitioned so some ranks own no fault faces would otherwise send those
   // ranks into the bulk sync loop while fault-bearing ranks take a different
   // path (the REVIEW A-1 divergence).  With the global count, EVERY rank agrees:
   // a genuinely fault-free mesh ⇒ all ranks run the bulk stepper; any fault ⇒
   // all ranks fall to the fault path (np=1 interleave / GTS).
   const bool lts_stepping =
      cfg.numerics.LtsEnabled() && num_fault_global == 0;
   // LTS Phase 3 (P3-4): the fault-half INTERLEAVE — each cluster advances its
   // fault QPs at its own rate.  OPT-IN via SEAS_LTS_FAULT_INTERLEAVE (default OFF
   // ⇒ the still-GTS reorder canary is preserved, and the default fault+lts science
   // is unchanged until the Frontera physics gate + the D-5 default flip).  Requires
   // np=1 (D-2: no shared fault faces) + the active P-006 reorder (cluster-contiguous
   // fault-QP ranges).
   // REVIEW P3-4 C1: parse the VALUE (not mere presence) — this flag changes the
   // science, so `SEAS_LTS_FAULT_INTERLEAVE=0` (or empty / "false") must DISABLE it,
   // not silently enable it (unlike the presence-based SEAS_DIAG_* debug gates).
   const char *fault_interleave_env = std::getenv("SEAS_LTS_FAULT_INTERLEAVE");
   const bool lts_fault_interleave_optin =
      fault_interleave_env != nullptr && fault_interleave_env[0] != '\0'
      && std::strcmp(fault_interleave_env, "0") != 0
      && std::strcmp(fault_interleave_env, "false") != 0;
   // LTS Phase 4a: the fault interleave now runs at np>1.  D-2 (fault-locality)
   // keeps every fault face rank-INTERIOR, so the fault half stays local; the BULK
   // half's rank seams are handled by the same seam machinery as the bulk stepper.
   // The gate MUST be PURELY rank-uniform (all three terms are the same on every
   // rank) — a collective branch decision must not hinge on per-rank state.
   // REVIEW p4-fault MODERATE: an earlier `(FaultFacesReordered() ||
   // num_fault_total==0)` clause could diverge (a rank owning only SHARED fault
   // faces has num_fault_total>0 but an empty interior reorder => false, while
   // clean ranks get true => deadlock).  The P-006 reorder + D-2 are instead
   // enforced fail-loud on ALL ranks by the fault stepper's ctor asserts
   // (GetNumSharedFaultQPs()==0 + per-cluster QP contiguity) and the
   // num_shared_global==0 check just below.
   const bool lts_fault_stepping =
      cfg.numerics.LtsEnabled() && num_fault_global > 0
      && lts_fault_interleave_optin;
   if (cfg.numerics.LtsEnabled() && num_fault_global > 0 && rank == 0)
   {
      if (lts_fault_stepping)
      {
         std::cout << "[lts] fault + lts=\"" << cfg.numerics.lts
                   << "\": fault-half LTS INTERLEAVE ACTIVE (P3-4; np" << nprocs
                   << ") — each cluster advances its fault QPs at dt_base*2^c.  "
                      "Fault faces are rank-interior (D-2); bulk rank seams use the "
                      "seam machinery.  Multi-cluster physics fidelity is Frontera-"
                      "gated; single-cluster == GTS (machine-eps) is the local gate.\n";
      }
      else
      {
         std::cout << "[lts] fault + lts=\"" << cfg.numerics.lts
                   << "\": fault-face reorder "
                   << (wave.FaultFacesReordered() ? "ACTIVE" : "inactive")
                   << " (cluster-contiguous, P-006) but stepping GTS — the fault-half "
                      "LTS interleave (P3-4) is OPT-IN via SEAS_LTS_FAULT_INTERLEAVE "
                      "(this is the still-GTS reorder canary).\n";
      }
   }
   bool lts_v2_checkpoint_written = false;
   if (lts_stepping)
   {
      // LTS Phase 4a: np>1 is enabled for the fault-free bulk stepper (the MPI
      // ghost-field seam exchange); serial clustering (Step 0) gives rank-count-
      // independent cluster ids, so the per-rank stepping is deterministic.
      // REVIEW OUT-2: the per-sync volume writer keys off the adaptive
      // slip-rate schedule, but a fault-free bulk run has no slip rate, so it
      // would pin the slowest (interseismic) regime and under-sample.  Warn so
      // the user sets a fixed time-based volume cadence (--paraview-volume-pv-dt
      // / paraview_fault_dt); free-surface + fixed-dt bulk-stress output are
      // unaffected (own schedules).
      if (rank == 0 && pv_out)
      {
         std::cout << "[lts] NOTE: per-sync volume ParaView uses a TIME-based "
                      "cadence (no slip-rate under fault-free bulk).  Set a fixed "
                      "volume dt (paraview_fault_dt / --paraview-volume-pv-dt); "
                      "an adaptive (slip-rate) schedule under-samples here.\n";
      }

      // Reuse the pre-operator clustering (single-sourced with the fault-face
      // reorder above; the pre-operator (pmesh, material) inputs give the SAME
      // cluster ids BuildLtsMeshInputs(wave, pmesh) would — identical dt_e
      // binning), so the reorder and the checkpoint layout_hash cannot disagree.
      MFEM_VERIFY(lts_layout_ready,
                  "spatial_dyn: LTS stepping requested but the pre-operator "
                  "clustering was not built (internal error).");
      const LtsClustering &cl     = lts_cl;
      const LtsLayout      &layout = lts_layout;
      const LtsGlobalMeta  &meta   = lts_meta;
      LtsClusteringOptions opt;   // default rate=2 for the layout hash below

      // Fine-cluster dt in seconds = lambda * dt_cfl (dt_cfl is the min stable
      // dt over the mesh — the fine cluster's dt; cluster c steps at dt_base*2^c,
      // CFL-safe by the GAP-A1 clustering guarantee).
      const real_t dt_base = cl.lambda * dt_cfl;
      const real_t T_s = dt_base * static_cast<real_t>(1LL << (cl.num_clusters - 1));
      // The hash uses the SERIAL cluster ids (cl.cluster, rank-independent — so a
      // restart is layout-checked independent of rank count); the stepper uses the
      // LOCAL ids (lts_cluster_id).  LTS Phase 4a: this block now runs at np>1
      // too (the bulk seam exchange); the SERIAL-vs-LOCAL split is exactly what
      // makes it rank-count-independent (they coincide only at np=1).
      const std::uint64_t layout_hash =
         LtsLayoutHash(opt.rate, cl.num_clusters, cl.cluster, dt_base, cl.lambda);

      LtsBulkSyncStepper<ParMesh> stepper(wave, layout, lts_cluster_id, Q,
                                          cfg.numerics.ader_order, dt_base);
      // LTS Phase 4a Stage 2: enable the coarse provider-D(k) ghost exchange for
      // CROSS-CLUSTER (diff-1) rank seams — only meaningful at np>1 with >1
      // cluster.  Must be set consistently on the stepper AND in `meta` (which
      // feeds BuildTickTable's matched-collective count).  `meta` is a reference
      // to `lts_meta`, so this update is visible to the loop below.
      lts_meta.exchange_bulk_provider_dk = (nprocs > 1 && cl.num_clusters > 1);
      stepper.SetExchangeProviderDk(lts_meta.exchange_bulk_provider_dk);
      // A1 (opt-in): the tick table PREDICTS the collective count, so the meta
      // flag and the stepper flag must be set from the same value or the driver's
      // matched-collective VERIFY aborts (P-007).
      lts_meta.merge_tick_seam_exchanges = lts_merge_exch;
      stepper.SetMergeTickExchanges(lts_merge_exch);
      if (rank == 0)
      {
         std::cout << "[lts] STEPPING (rate2, fault-free bulk): Nc = "
                   << cl.num_clusters << ", dt_base = " << dt_base
                   << " s, T_sync = " << T_s << " s, predicted harmonic speedup "
                   << HarmonicUpdateSpeedup(cl) << "x, layout_hash = 0x"
                   << std::hex << layout_hash << std::dec << "\n";
      }

      real_t t_lts = t;
      // REVIEW P3-4 C2: seed the sync/cycle counter from the restored step0 (a
      // sync count for the LTS V2 checkpoint) so output cycles + logs CONTINUE
      // across --restart instead of restarting at 0 and colliding pre-restart
      // ParaView frames.  Fresh run ⇒ step0 == 0.
      int sync = step0;
      while (t_lts < cfg.time.tfinal - 1e-12 * T_s)
      {
         MFEM_PERF_SCOPE("seas::spatial_dyn::step");
         const real_t T_actual = std::min(T_s, cfg.time.tfinal - t_lts);
         if (T_actual <= 0.0) { break; }
         auto tab = BuildTickTable(cl.num_clusters, dt_base, T_actual,
                                   cfg.numerics.ader_order, meta);
         stepper.SetSyncInterval(t_lts, T_actual);
         wave.SetTime(t_lts);
         // LTS Phase 4a (P-007): the ghost seam exchanges made this sync must
         // equal the tick table's summed n_collectives on EVERY rank that OWNS
         // shared faces (matched collectives — a mismatch is the class of bug that
         // hangs np>=10).  A rank with no shared faces (np=1, or an isolated
         // subdomain at np>1) does zero exchanges — expect 0 there, else the
         // check would false-abort a correct run (REVIEW p4a MODERATE).
         wave.ResetGhostExchangeCount();
         long long sched_exchanges = 0;
         for (const auto &tk : tab) { sched_exchanges += tk.n_collectives; }
         const long long expected_exchanges =
            wave.HasSharedFaces() ? sched_exchanges : 0;
         RunSyncInterval(tab, stepper);
         MFEM_VERIFY(wave.GhostExchangeCount() == expected_exchanges,
                     "spatial_dyn: LTS ghost-exchange count "
                     << wave.GhostExchangeCount() << " != expected "
                     << expected_exchanges << " (tick-table n_collectives "
                     << sched_exchanges << ", has_shared="
                     << wave.HasSharedFaces() << ") at sync " << sync
                     << " (matched-collective violation, P-007).");
         MFEM_VERIFY(stepper.BuffersZero(),
                     "spatial_dyn: LTS accumulate buffers non-zero at sync point "
                     << sync << " (buffer lifecycle bug — invariant ii).");
         t_lts += T_actual;
         ++sync;
         // Per-sync bulk output: NaN tripwire + volume ParaView (fault-free, so
         // no station/fault output) + a V2 (LTS) checkpoint on the sync cadence.
         MFEM_VERIFY(Q.CheckFinite() == 0,
                     "spatial_dyn: LTS Q non-finite at sync " << sync
                     << ", t = " << t_lts << " s (NaN tripwire).");
         paraview_write(sync, t_lts, /*V_max=*/0.0);   // volume/free-surface (fault-free)
         if (cfg.output.checkpoint_every_steps > 0
             && (sync % cfg.output.checkpoint_every_steps == 0))
         {
            const std::string prefix = cfg.output.output_dir + "/" + cfg.output.restart_prefix;
            mfem::seas::internal::WriteTpv104CheckpointV2Impl(prefix, t_lts, dt_base, sync,
                                                  /*lts_mode=*/1, layout_hash, Q,
                                                  dof_data, rank, nprocs, "spatial_dyn");
         }
         if (rank == 0 && (sync % 100 == 0))
         { std::cout << "[lts] sync " << sync << ", t = " << t_lts << " s\n"; }
      }
      t = t_lts;
      last_completed_step = sync;
      // Final V2 checkpoint (unless the last sync already wrote one on cadence).
      if (cfg.output.checkpoint_every_steps > 0
          && !(sync > 0 && sync % cfg.output.checkpoint_every_steps == 0))
      {
         const std::string prefix = cfg.output.output_dir + "/" + cfg.output.restart_prefix;
         mfem::seas::internal::WriteTpv104CheckpointV2Impl(prefix, t, dt_base, sync,
                                               /*lts_mode=*/1, layout_hash, Q,
                                               dof_data, rank, nprocs, "spatial_dyn");
      }
      lts_v2_checkpoint_written = (cfg.output.checkpoint_every_steps > 0);
      if (rank == 0)
      {
         std::cout << "[lts] done: " << sync << " sync intervals, t = " << t
                   << " s (fault-free bulk; per-sync volume output + V2 "
                      "checkpoints written; fault/station output is Phase 3).\n";
      }
   }

   // ------------------------------------------------------------------
   // (LTS Phase 3, P3-4) Fault-half INTERLEAVE sync loop (np=1, opt-in).
   // ------------------------------------------------------------------
   // Mirrors the fault-free bulk sync loop above but drives the fault-aware
   // LtsFaultSyncStepper and writes per-sync fault + station output.  Each due
   // cluster advances its fault QPs over its cluster-contiguous global-QP range
   // at dt_base*2^c (range friction Advance + absolute range nucleation + the
   // cluster corrector's fault-face flux).  Single-cluster == GTS (machine-eps)
   // is the local gate; multi-cluster physics fidelity is Frontera-staged.
   if (lts_fault_stepping)
   {
      // LTS Phase 4a: np>1 enabled.  D-2 (fault-locality) keeps fault faces
      // rank-interior; the bulk half's rank seams use the seam machinery.
      // Enforce D-2 GLOBALLY + uniformly (num_shared_global is an Allreduce, so
      // this aborts on ALL ranks together if fault-locality was not applied — the
      // fault interleave has no cross-rank shared-fault path).
      MFEM_VERIFY(num_shared_global == 0,
                  "spatial_dyn: LTS fault interleave requires fault-locality (D-2): "
                  "expected zero SHARED fault faces globally, got " << num_shared_global
                  << ".  lts != off implies fault-locality (P-017) — check the "
                  "partition (cross-rank shared-fault stepping is not supported).");
      MFEM_VERIFY(lts_layout_ready,
                  "spatial_dyn: LTS fault interleave requested but the pre-operator "
                  "clustering was not built (internal error).");
      const LtsClustering &cl     = lts_cl;
      const LtsLayout      &layout = lts_layout;
      const LtsGlobalMeta  &meta   = lts_meta;
      LtsClusteringOptions opt;   // default rate=2 for the layout hash
      const real_t dt_base = cl.lambda * dt_cfl;
      const real_t T_s = dt_base * static_cast<real_t>(1LL << (cl.num_clusters - 1));
      const std::uint64_t layout_hash =
         LtsLayoutHash(opt.rate, cl.num_clusters, cl.cluster, dt_base, cl.lambda);

      LtsFaultSyncStepper<ParMesh> stepper(
         wave, layout, lts_cluster_id, Q, cfg.numerics.ader_order, dt_base,
         substep_iterator, dof_data, fault_coords, nuc.get());
      // LTS Phase 4a Stage 2: enable the bulk diff-1 seam D(k) exchange at np>1
      // with >1 cluster (consistent with the tick-table matched-collective count).
      lts_meta.exchange_bulk_provider_dk = (nprocs > 1 && cl.num_clusters > 1);
      stepper.SetExchangeProviderDk(lts_meta.exchange_bulk_provider_dk);
      // A1 is NOT wired into the fault-interleave stepper yet (it has no
      // BeforeCorrects override), so enabling the merge here would make the tick
      // table predict 64 rounds while the stepper still fires 125 -> a
      // matched-collective ABORT mid-run.  Refuse up front instead.
      MFEM_VERIFY(!lts_merge_exch,
                  "spatial_dyn: --lts-merge-exchanges is not supported on the LTS "
                  "fault-interleave path yet (bulk-only; see document/comm_dev/"
                  "DESIGN_a1_per_tick_merge_2026-07-24.md).  Re-run without it.");
      if (rank == 0)
      {
         std::cout << "[lts] STEPPING (rate2, fault interleave): Nc = "
                   << cl.num_clusters << ", dt_base = " << dt_base
                   << " s, T_sync = " << T_s << " s, predicted harmonic speedup "
                   << HarmonicUpdateSpeedup(cl) << "x, layout_hash = 0x"
                   << std::hex << layout_hash << std::dec << "\n";
      }
      // REVIEW P3-4 C3: at Nc>1 the fault + SCEC-station output fires once per
      // SYNC (period T_sync = dt_base*2^(Nc-1), the COARSEST cluster's step), so a
      // fast near-fault fine-cluster signal is under-sampled — the SCEC .dat is the
      // benchmark artifact, so this aliases the very signal being validated.  This
      // coarse-cadence sampling is acceptable for the Nc=1 local gate; sub-sync
      // fault/station sampling is a Phase-4 item.  Warn loudly for Nc>1 runs.
      if (rank == 0 && cl.num_clusters > 1)
      {
         std::cout << "[lts] WARNING: fault/station output is written once per SYNC "
                      "(coarse cadence T_sync = " << T_s << " s) — at Nc = "
                   << cl.num_clusters << " > 1 this UNDER-SAMPLES the near-fault "
                      "signal.  Multi-cluster fault-output fidelity + physics "
                      "acceptance are Frontera-staged (sub-sync sampling is Phase 4); "
                      "the local gate is Nc=1 == GTS.\n";
      }

      real_t t_lts = t;
      // REVIEW P3-4 C2: seed the sync/cycle counter from the restored step0 (a
      // sync count for the LTS V2 checkpoint) so output cycles + logs CONTINUE
      // across --restart instead of restarting at 0 and colliding pre-restart
      // ParaView frames.  Fresh run ⇒ step0 == 0.
      int sync = step0;
      while (t_lts < cfg.time.tfinal - 1e-12 * T_s)
      {
         MFEM_PERF_SCOPE("seas::spatial_dyn::step");
         const real_t T_actual = std::min(T_s, cfg.time.tfinal - t_lts);
         if (T_actual <= 0.0) { break; }
         // Per-sync reset of the honest sub-step |V| max (mirrors the GTS loop's
         // per-step reset); the friction range Advance takes running max over it.
         for (int i = 0; i < num_fault_total; ++i)
         { dof_data[i].slip_rate_substep_max = 0.0; }

         auto tab = BuildTickTable(cl.num_clusters, dt_base, T_actual,
                                   cfg.numerics.ader_order, meta);
         stepper.SetSyncInterval(t_lts, T_actual);
         wave.SetTime(t_lts);
         // LTS Phase 4a (P-007): matched-collective audit — the bulk seam ghost
         // exchanges this sync must equal the tick-table n_collectives on every
         // rank that owns shared faces (fault faces are rank-interior, D-2, so the
         // fault half adds no exchanges).
         wave.ResetGhostExchangeCount();
         long long sched_exchanges = 0;
         for (const auto &tk : tab) { sched_exchanges += tk.n_collectives; }
         const long long expected_exchanges =
            wave.HasSharedFaces() ? sched_exchanges : 0;
         RunSyncInterval(tab, stepper);
         MFEM_VERIFY(wave.GhostExchangeCount() == expected_exchanges,
                     "spatial_dyn: LTS fault-path ghost-exchange count "
                     << wave.GhostExchangeCount() << " != expected "
                     << expected_exchanges << " at sync " << sync
                     << " (matched-collective violation, P-007).");
         MFEM_VERIFY(stepper.BuffersZero(),
                     "spatial_dyn: LTS accumulate buffers non-zero at sync point "
                     << sync << " (buffer lifecycle bug — invariant ii).");
         MFEM_VERIFY(wave.SeamCoarseBuffersZero(),
                     "spatial_dyn: LTS coarse-side seam in-tray non-zero at sync "
                     << sync << " (Stage-2 buffer lifecycle bug).");
         t_lts += T_actual;
         ++sync;

         real_t V_max_local = 0.0;
         for (int i = 0; i < num_fault_total; ++i)
         {
            V_max_local = std::max(V_max_local,
               std::max(dof_data[i].slip_rate, dof_data[i].slip_rate_substep_max));
         }
         real_t V_max_step = V_max_local;
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&V_max_local, &V_max_step, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
         V_max_global = std::max(V_max_global, V_max_step);

         MFEM_VERIFY(Q.CheckFinite() == 0,
                     "spatial_dyn: LTS Q non-finite at sync " << sync
                     << ", t = " << t_lts << " s (NaN tripwire).");
         paraview_write(sync, t_lts, V_max_step);         // fault + volume/free-surface
         if (stations_write) { stations_write(t_lts, dof_data); }
         if (cfg.output.checkpoint_every_steps > 0
             && (sync % cfg.output.checkpoint_every_steps == 0))
         {
            const std::string prefix =
               cfg.output.output_dir + "/" + cfg.output.restart_prefix;
            mfem::seas::internal::WriteTpv104CheckpointV2Impl(
               prefix, t_lts, dt_base, sync, /*lts_mode=*/1, layout_hash, Q,
               dof_data, rank, nprocs, "spatial_dyn", qp_canon_ptr);
         }
         if (rank == 0 && (sync % 100 == 0))
         {
            std::cout << "[lts] sync " << sync << ", t = " << t_lts
                      << " s, V_max = " << V_max_step << " m/s\n";
         }
      }
      t = t_lts;
      last_completed_step = sync;
      if (cfg.output.checkpoint_every_steps > 0
          && !(sync > 0 && sync % cfg.output.checkpoint_every_steps == 0))
      {
         const std::string prefix =
            cfg.output.output_dir + "/" + cfg.output.restart_prefix;
         mfem::seas::internal::WriteTpv104CheckpointV2Impl(
            prefix, t, dt_base, sync, /*lts_mode=*/1, layout_hash, Q,
            dof_data, rank, nprocs, "spatial_dyn", qp_canon_ptr);
      }
      lts_v2_checkpoint_written = (cfg.output.checkpoint_every_steps > 0);
      if (rank == 0)
      {
         std::cout << "[lts] done: " << sync << " sync intervals (fault "
                      "interleave), t = " << t << " s, V_max = " << V_max_global
                   << " m/s.\n";
      }
   }

   for (int step = step0;
        step < nsteps && !lts_stepping && !lts_fault_stepping; ++step)
   {
      MFEM_PERF_SCOPE("seas::spatial_dyn::step");
      const real_t dt_step = std::min(dt_now, cfg.time.tfinal - t);
      if (dt_step <= 0.0) { break; }
      wave.SetTime(t);

      // PLAN_speckle_slip_runaway Phase 1 (R-004): reset the per-macro-step
      // honest |V| max before the sub-step solve takes running max over it
      // (iterator sub-steps + shared-fault macro solve).  Transient
      // diagnostic field only; does not feed back into state.
      for (int i = 0; i < num_fault_total; ++i)
      {
         dof_data[i].slip_rate_substep_max = 0.0;
      }
      // [DIAG-SIGN] NOTE (R-003): sigma_n_substep_min is intentionally NOT
      // reset here.  It accumulates the most-tensile sub-step normal traction
      // across the WHOLE diag interval and is reset only inside the diag block
      // after it has been reported, so a tensile transient on a step the diag
      // does not sample (step%100!=0, V<10) is still captured at the next
      // print.  slip_rate_substep_max keeps its per-macro-step reset above.

      // Phase 14: ONLY change to the stepping control is this branch; the
      // Q.Swap(Q_new) and t += dt_step below are shared with the ADER branch.
      if (is_rk && is_lsw)
      {
         // LSW coupled-RK stepper: drives wave.Mult directly, couples (Q, slip)
         // with the tableau weights (no ψ — LSW's μ(δ) depends on accumulated
         // slip, which is staged before each Mult), and applies absolute
         // nucleation at each stage time (no-op for TPV205).  Scalar-flux
         // (guarded at setup).  The per-stage max-|V| is reduced into
         // slip_rate_substep_max.  Mult dispatches to the slip-stateless
         // EvaluateLSW because fault_friction_law_ == LSW (set above).
         AdvanceRKCoupledLSW_Spatial(wave, dof_data, Q, dt_step, t, Q_new,
                                     rk_tab, nuc.get());
      }
      else if (is_rk)
      {
         // Rate-state RK4 / RK45 coupled stepper: drives wave.Mult directly,
         // couples (Q, ψ, slip) with the tableau weights, and applies §14.3
         // absolute nucleation at each stage time (nuc_cb / the substep
         // iterator are NOT used).  Scalar-flux (guarded at setup).  The
         // per-stage max-|V| is reduced into slip_rate_substep_max (§14.5).
         AdvanceRKCoupled_Spatial(wave, dof_data, *cfg.rate_state, rs,
                                  Q, dt_step, t, Q_new, rk_tab, nuc.get());
      }
      else
      {
         AdvanceADERWithSubStep_Spatial(wave, substep_iterator, dof_data,
                                        fault_coords, Q, dt_step,
                                        cfg.numerics.ader_order, t, Q_new,
                                        nuc_cb, use_shared_ck);
      }
      Q.Swap(Q_new);
      t += dt_step;
      last_completed_step = step + 1;

      real_t V_max_local = 0.0;
      for (int i = 0; i < num_fault_total; ++i)
      {
         V_max_local = std::max(V_max_local, dof_data[i].slip_rate);
      }
      real_t V_max_step = V_max_local;
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&V_max_local, &V_max_step, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
      V_max_global = std::max(V_max_global, V_max_step);

      // DIAG (env-gated, debug only): rupture-area proxy + onset localizer.
      // SEAS_DIAG_BLOWUP=1 enables.  Prints (a) global count of fault DOFs
      // with V>0.5 m/s (rupturing-area proxy: growing => propagating front,
      // ~constant => only the forced patch) + global max slip, and (b) at
      // onset (V_max>10) the owning rank dumps the argmax DOF's full state
      // so we can see WHAT runs away first (sigma_n? shear? slip?).
      if (std::getenv("SEAS_DIAG_BLOWUP")
          && (step % 100 == 0 || V_max_step > 10.0))
      {
         long long n_rup_local = 0;
         real_t    maxslip_local = 0.0;
         int       argmax_local = -1;
         real_t    vloc = -1.0;
         // Phase 1 (R-004): honest sub-step-aware |V| = max(slip_rate,
         // slip_rate_substep_max) so the print cannot under-report vs the
         // last-sub-step-only slip_rate (shared-QP / reconcile safe).
         real_t    vsub_local = 0.0;
         // [DIAG-SIGN] speckle instrumentation (sliver_blowup plan 2026-05-26):
         // track the MOST-TENSILE normal traction and how many fault DOFs have
         // the sigma_n strength floor ENGAGED (sigma_n_corr < floor).  Two
         // questions this answers directly: (a) "is the cap triggered at all?"
         // -> n_below_floor; (b) "where/when does sigma_n go tensile?" ->
         // sigma_n_min + its DOF (the speckle onset, which precedes the V>10
         // blowup the [DIAG-ONSET] line catches).  NB: data.sigma_n_corr ==
         // the flux's sigma_n_fric (= sigma_n_trial + sigma_n_nuc), so
         // (sigma_n_corr < floor) is EXACTLY fault_face_flux.cpp's engagement
         // test.  The floor only clamps the shear STRENGTH, not this imposed
         // normal traction (sigma_n_corr = sigma_n_trial), so tensile excursions
         // here still drive the bulk normal Riemann update ([[v_n]] opening).
         const real_t SN_UNSET = std::numeric_limits<real_t>::max();
         const real_t sn_floor = cfg.sigma_n_strength_floor_pa;  // < 0 = disabled
         real_t    sigman_min_local = SN_UNSET;    // end-of-step sigma_n_corr
         int       argmin_sn_local  = -1;
         long long n_tensile_local     = 0;   // sigma_n_corr < 0 (tensile)
         long long n_below_floor_local = 0;   // sigma_n_corr < floor (engaged)
         // SUB-STEP trackers (the key signal: a tensile transient that the
         // end-of-step sigma_n_corr above has already recovered from).
         real_t    ss_sigman_min_local = SN_UNSET;  // interval min over sub-steps
         long long n_ss_tensile_local     = 0;  // sigma_n_substep_min < 0
         long long n_ss_below_floor_local = 0;  // sigma_n_substep_min < floor
         for (int i = 0; i < num_fault_total; ++i)
         {
            const DOFData &d = dof_data[i];
            if (d.slip_rate > 0.5) { n_rup_local++; }
            const real_t s = std::sqrt(d.slip1 * d.slip1 + d.slip2 * d.slip2);
            if (s > maxslip_local) { maxslip_local = s; }
            if (d.slip_rate > vloc) { vloc = d.slip_rate; argmax_local = i; }
            vsub_local = std::max(vsub_local,
                                  std::max(d.slip_rate, d.slip_rate_substep_max));
            const real_t sn = d.sigma_n_corr;
            if (sn < sigman_min_local) { sigman_min_local = sn; argmin_sn_local = i; }
            if (sn < 0.0) { n_tensile_local++; }
            if (sn_floor >= 0.0 && sn < sn_floor) { n_below_floor_local++; }
            const real_t snss = d.sigma_n_substep_min;
            if (snss < SN_UNSET)   // a sub-step value was recorded this interval
            {
               if (snss < ss_sigman_min_local) { ss_sigman_min_local = snss; }
               if (snss < 0.0) { n_ss_tensile_local++; }
               if (sn_floor >= 0.0 && snss < sn_floor) { n_ss_below_floor_local++; }
            }
         }
         long long n_rup_g = n_rup_local;
         real_t    maxslip_g = maxslip_local;
         real_t    vsub_g = vsub_local;
         struct { double v; int r; } in_{V_max_local, rank}, out_{V_max_local, rank};
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&n_rup_local, &n_rup_g, 1, MPI_LONG_LONG, MPI_SUM, comm);
         MPI_Allreduce(&maxslip_local, &maxslip_g, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
         MPI_Allreduce(&in_, &out_, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
         MPI_Allreduce(&vsub_local, &vsub_g, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
         if (rank == 0)
         {
            std::cout << "[DIAG] step " << step << " t=" << t
                      << " V_max=" << V_max_step
                      << " V_substep_max=" << vsub_g
                      << " n_rupturing(V>0.5)=" << n_rup_g
                      << " max_slip=" << maxslip_g << " m\n";
         }
         if (V_max_step > 10.0 && rank == out_.r && argmax_local >= 0)
         {
            const DOFData &d = dof_data[argmax_local];
            std::cout << "[DIAG-ONSET] rank " << rank << " dof " << argmax_local
                      << " xyz=(" << dof_coords_3d(3 * argmax_local) << ","
                      << dof_coords_3d(3 * argmax_local + 1) << ","
                      << dof_coords_3d(3 * argmax_local + 2) << ")"
                      << " V=" << d.slip_rate << " V1=" << d.V1
                      << " V2=" << d.V2 << " slip1=" << d.slip1
                      << " slip2=" << d.slip2
                      << " sigma_n_corr=" << d.sigma_n_corr
                      << " tau1_corr=" << d.tau1_corr
                      << " tau2_corr=" << d.tau2_corr
                      << " tau1_nuc=" << d.tau1_nuc
                      << " tau2_nuc=" << d.tau2_nuc << "\n";
         }

         // [DIAG-SIGN]: tensile-normal-stress / floor-engagement tracker.
         // Two timescales: end-of-step (sigma_n_corr, what ParaView shows) AND
         // sub-step (sigma_n_substep_min).  A LARGE n_ss_below_floor with a
         // SMALL n_below_floor is the signature of a sub-step tensile transient
         // that recovered by output time — i.e. the cap engages on the sub-step
         // even though the output field is never below the floor (resolves the
         // "I checked those points, sigma_n is not below 10 MPa" puzzle).
         long long n_tensile_g = n_tensile_local;
         long long n_below_floor_g = n_below_floor_local;
         long long n_ss_tensile_g = n_ss_tensile_local;
         long long n_ss_below_floor_g = n_ss_below_floor_local;
         struct { double v; int r; } sn_in{sigman_min_local, rank},
                                     sn_out{sigman_min_local, rank},
                                     ss_in{ss_sigman_min_local, rank},
                                     ss_out{ss_sigman_min_local, rank};
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&n_tensile_local, &n_tensile_g, 1, MPI_LONG_LONG,
                       MPI_SUM, comm);
         MPI_Allreduce(&n_below_floor_local, &n_below_floor_g, 1, MPI_LONG_LONG,
                       MPI_SUM, comm);
         MPI_Allreduce(&n_ss_tensile_local, &n_ss_tensile_g, 1, MPI_LONG_LONG,
                       MPI_SUM, comm);
         MPI_Allreduce(&n_ss_below_floor_local, &n_ss_below_floor_g, 1,
                       MPI_LONG_LONG, MPI_SUM, comm);
         MPI_Allreduce(&sn_in, &sn_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
         MPI_Allreduce(&ss_in, &ss_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
#endif
         // R-006: once V>10 the diag block fires every step; print the
         // aggregate line only every 100 steps OR when something is actually
         // tensile, so the log does not bloat with all-compressive lines.
         // (step + n_ss_tensile_g are globally consistent, so all ranks agree
         // -> no divergent collective / deadlock and a consistent reset below.)
         const bool print_sign = (step % 100 == 0 || n_ss_tensile_g > 0);
         if (rank == 0 && print_sign)
         {
            std::cout << "[DIAG-SIGN] step " << step << " t=" << t
                      << " floor="
                      << (sn_floor >= 0.0 ? sn_floor / 1.0e6 : -1.0) << " MPa"
                      << (sn_floor < 0.0 ? " (DISABLED)" : "")
                      << " | end-of-step: sigma_n_min=" << sn_out.v
                      << " n_tensile=" << n_tensile_g
                      << " n_below_floor=" << n_below_floor_g
                      << " | SUB-STEP(interval): sigma_n_min=" << ss_out.v
                      << " n_tensile=" << n_ss_tensile_g
                      << " n_below_floor=" << n_ss_below_floor_g << "\n";
         }
         if (print_sign)
         {
            // R-001: EVERY rank dumps ALL its local DOFs whose interval sub-step
            // sigma_n is below the floor (or tensile, when the cap is off), so
            // EVERY concurrent speckle spot is localized -- not just the single
            // global-worst MINLOC DOF.  Capped per rank to bound log volume.
            const real_t sign_thr = (sn_floor >= 0.0) ? sn_floor : 0.0;
            int dumped = 0;
            for (int i = 0; i < num_fault_total && dumped < 32; ++i)
            {
               const DOFData &d = dof_data[i];
               if (d.sigma_n_substep_min < SN_UNSET &&
                   d.sigma_n_substep_min < sign_thr)
               {
                  std::cout << "[DIAG-SIGN-DOF] rank " << rank << " dof " << i
                            << " xyz=(" << dof_coords_3d(3 * i) << ","
                            << dof_coords_3d(3 * i + 1) << ","
                            << dof_coords_3d(3 * i + 2) << ")"
                            << " sigma_n_substep_min=" << d.sigma_n_substep_min
                            << " sigma_n_corr(end)=" << d.sigma_n_corr
                            << " V=" << d.slip_rate
                            << " V_substep_max=" << d.slip_rate_substep_max
                            << " tau1_corr=" << d.tau1_corr
                            << " tau2_corr=" << d.tau2_corr << "\n";
                  ++dumped;
               }
            }
            // R-003: the interval accumulator has now been reported; reset it on
            // EVERY rank (consistent gate) so the next interval starts fresh.
            for (int i = 0; i < num_fault_total; ++i)
            {
               dof_data[i].sigma_n_substep_min = SN_UNSET;
            }
         }
      }

      paraview_write(step + 1, t, V_max_step);

      // Append the post-step state to each SCEC station trace (every step;
      // each station is ~80 B/row).  No-op unless a [meta].tag-matched writer
      // was wired above (tpv31/102/104/205).
      if (stations_write) { stations_write(t, dof_data); }

      // ------------------------------------------------------------------
      // MPI shared-fault consistency tripwire.  Follows tpv104_driver.cpp
      // (:2636-2639) and tpv205_driver.cpp (:2442-2445), both of which call
      // this once at step 0; the spatial driver had DROPPED it entirely.
      //
      // Extended past TPV104's step-0-only call to fire through the whole
      // gradual_overstress nucleation window: on the CURVILINEAR SAFS fault
      // the per-QP dip/strike frame (FaultBasis::ComputeOrientedFrame, via
      // sign(n_raw . ref_normal) with ref_normal=(0,-1,0)) is ill-conditioned
      // where the local normal is ~perpendicular to ref_normal, so the two
      // ranks owning a shared face can pick OPPOSITE strike directions.  The
      // fixed-sign nucleation increment (+F.dtau into tau2_nuc) then forces
      // the two sides of a shared face in opposite physical directions.  The
      // routine pairs same-physical-QP records across ranks and ABORTS with a
      // field/QP diagnostic if the 8 evolved fields differ > tol (so the
      // first offending shared QP inside the nucleation patch is pinpointed).
      // TPV104's planar fault has |n_raw . ref_normal| = 1 everywhere, so it
      // never flips and always passes this check.
      // SEAS_R101_SKIP (env-gated): skip the shared-fault consistency tripwire
      // ENTIRELY.  Default (unset) runs it (TPV/BP5 byte-exact — those configs
      // never set it).  Set for the CURVILINEAR SAFS fault, where the per-QP
      // strike frame can flip between the two ranks owning a shared face and trip
      // this check during nucleation; only ~0.02% of fault DOFs are on shared
      // faces, so skipping it leaves the bulk rupture unaffected.  The env is
      // rank-uniform (exported before srun), so every rank skips together — no
      // rank-conditional collective.  (SEAS_R101_NONFATAL below is the milder
      // option: keep the check but warn instead of abort.)
      if (std::getenv("SEAS_R101_SKIP") == nullptr
          && (step == 0
              || (cfg.nucleation.enabled
                  && t <= cfg.nucleation.gradual_overstress.T_nuc_s
                  && (step % 100 == 0))))
      {
         // SEAS_R101_NONFATAL (env-gated, OFF by default => zero behaviour
         // change for TPV/BP5/production SAFS): downgrade this shared-fault
         // consistency tripwire from a hard MFEM_ABORT to a per-check rank-0
         // diagnostic.  Lets a run with a KNOWN cross-rank DOFData divergence
         // (the curvilinear-fault strike-frame flip) proceed to completion so
         // we can test whether a candidate fix (e.g. the fault-normal
         // canonicalization change in fault_basis.hpp) keeps worst_rel bounded
         // AND avoids the V_max blow-up.  The routine still runs every MPI
         // collective; abort_on_fail=false just suppresses the abort and writes
         // worst_rel / worst_field instead (wave_operator.inl:5826).
         const bool r101_nonfatal =
            (std::getenv("SEAS_R101_NONFATAL") != nullptr);
         double r101_worst_rel = 0.0;
         int    r101_worst_field = -1;
         wave.VerifySharedFaultDOFDataConsistency(
            /*tol=*/1e-10, &r101_worst_rel, &r101_worst_field,
            /*abort_on_fail=*/!r101_nonfatal);
         if (r101_nonfatal && rank == 0)
         {
            std::cout << "  [R-101 nonfatal] step " << step << " t=" << t
                      << "s  worst_rel=" << r101_worst_rel
                      << "  worst_field=" << r101_worst_field
                      << (r101_worst_rel > 1e-10 ? "  (DIVERGED)"
                                                 : "  (ok)")
                      << std::endl;
         }
      }

      // NaN tripwire — follows tpv104_driver.cpp:2641-2659 /
      // tpv205_driver.cpp:2447-2464 (also dropped by the spatial driver).
      {
         real_t local_nan = std::isnan(Q.Norml2()) ? 1.0 : 0.0;
         real_t global_nan = local_nan;
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&local_nan, &global_nan, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
         if (global_nan > 0.0)
         {
            std::cerr << "ERROR: NaN detected at step " << step
                      << ", t = " << t << " s (rank " << rank << ")\n";
#ifdef MFEM_USE_MPI
            MPI_Finalize();
#endif
            return 1;
         }
      }

      if (cfg.output.checkpoint_every_steps > 0
          && (step + 1) % cfg.output.checkpoint_every_steps == 0)
      {
         const std::string prefix = cfg.output.output_dir + "/"
                                    + cfg.output.restart_prefix;
         write_gts_checkpoint(prefix, t, dt_now, step + 1);
      }

      // R-008: also print on the step that exhausts `tfinal` so the
      // "final-step" log fires even when the loop exits via the early
      // `dt_step <= 0.0` break (when `nsteps = ceil(tfinal/dt)` over-
      // shoots, the last advance lands on a step != nsteps - 1).
      const bool reached_tfinal =
         (t + std::numeric_limits<real_t>::epsilon() >= cfg.time.tfinal);
      if (rank == 0 && (step % 100 == 0
                        || step == nsteps - 1
                        || reached_tfinal))
      {
         std::cout << "step " << step << "/" << nsteps
                   << "  t = " << t << " s"
                   << "  V_max = " << V_max_step << " m/s\n";
      }
   }

   // Flush + close the SCEC station traces (no-op unless a writer was wired).
   // Files are flushed every WriteStep, so this is belt-and-suspenders before
   // the (RAII) close at scope exit.
   if (stations_close) { stations_close(); }

   // -----------------------------------------------------------------
   // 21. Final checkpoint.
   // -----------------------------------------------------------------
   // R-011: skip when the last in-loop checkpoint already covered
   // `last_completed_step` — the two writes would be byte-identical
   // and waste Lustre metadata ops on production.
   // (LTS Phase 2, REVIEW DL-1) On the LTS path `last_completed_step` is a
   // sync-interval count, NOT a GTS step count, so it must not flow into the
   // step-based checkpoint modulus below.  The LTS loop already wrote V2
   // (layout-hash) checkpoints on the sync cadence + a final one; do NOT also run
   // the step-based V1 write here (that would emit a V1 file an LTS restart
   // refuses).  The GTS path (else) is byte-identical to before.
   if (lts_stepping || lts_fault_stepping)
   {
      // (LTS Phase 2/3) `last_completed_step` is a sync-interval count, not a GTS
      // step count; the LTS loop already wrote V2 (layout-hash) checkpoints on the
      // sync cadence + a final one, so the step-based V1 write below must NOT run
      // (an LTS restart refuses a V1 file).
      if (rank == 0 && lts_v2_checkpoint_written)
      {
         std::cout << "[checkpoint] LTS V2 checkpoint(s) written on the sync "
                      "cadence; skipping the step-based V1 final write.\n";
      }
   }
   else
   {
   const bool already_checkpointed_final =
      (cfg.output.checkpoint_every_steps > 0)
      && (last_completed_step > step0)
      && (last_completed_step % cfg.output.checkpoint_every_steps == 0);
   if (cfg.output.checkpoint_every_steps > 0 && !already_checkpointed_final)
   {
      const std::string prefix = cfg.output.output_dir + "/"
                                 + cfg.output.restart_prefix;
      // R-603 round-6: use last_completed_step instead of nsteps so the
      // checkpoint reports the actual step the loop reached.
      write_gts_checkpoint(prefix, t, dt_now, last_completed_step);
   }
   else if (already_checkpointed_final && rank == 0)
   {
      std::cout << "[checkpoint] last in-loop checkpoint already covers "
                << "step " << last_completed_step
                << "; skipping redundant final write.\n";
   }
   }   // end else (!lts_stepping) — REVIEW DL-1

   if (rank == 0)
   {
      std::cout << "[spatial_dyn] done.  Final t = " << t
                << " s, V_max_global = " << V_max_global << " m/s\n";
   }

#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return 0;
}
