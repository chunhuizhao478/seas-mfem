// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper bug: first-step audit on the 8-triangle adjacent-fault fixture.
//
// This test narrows the first observable symmetry break reported by
// `test_adjacent_triangle_fault_uniformity`:
//   1. replay the step-0 ADER fault solve on the 8-triangle fixture and
//      assert the fault-local outputs are uniform across all fault QPs;
//   2. manually assemble ONLY the fault-face contribution into the bulk and
//      check symmetry of the fault-adjacent tetrahedra before and after
//      element mass inverse;
//   3. verify that this manual fault-only increment matches the production
//      `AdvanceADER` first step from Q=0 exactly.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/precomputed_face_fluxes.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../dynamic/d4_tet_mesh.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <tuple>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   const double vv = (v), tt = (tol); \
   if (vv <= tt) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(6) << vv << ", tol " << tt << ")\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(6) \
                << vv << ", expected <= " << tt << ")\n"; \
   } \
} while (0)

namespace {

constexpr real_t kL = 1000.0;
constexpr real_t kDt = 5.0e-5;
constexpr int kOrder = 1;
constexpr int kAderOrder = 2;
int g_nx = 2;
int g_ny = 2;
int g_nz = 2;

// Round-11 Step 2: branch-isolation toggles for the non-fault face
// audit helpers.  When either flag is '1', the corresponding branch
// (boundary face accumulation, or interior non-fault face
// accumulation) is skipped in RunADERNonFaultFaceAudit,
// RunPrecomputedFluxLiftedAudit, and RunAnalyticTraceLiftAudit.
// Used to determine whether the first asymmetry in Gates 4 / 7 /
// 14 / 16 comes from boundary flux, interior non-fault flux, or
// neither.
bool EnvFlagSet(const char *name)
{
   const char *v = std::getenv(name);
   return v && v[0] == '1';
}

const char *kCompName[NUM_STATE] = {
   "SXX", "SYY", "SZZ", "SXY", "SYZ", "SXZ", "VX", "VY", "VZ"
};

enum class AffineFieldKind
{
   X,
   Z,
   XPlusZ
};

const char *AffineFieldName(AffineFieldKind kind)
{
   switch (kind)
   {
      case AffineFieldKind::X: return "x";
      case AffineFieldKind::Z: return "z";
      case AffineFieldKind::XPlusZ: return "x+z";
   }
   MFEM_ABORT("unknown affine field kind");
}

real_t EvaluateAffineField(AffineFieldKind kind, real_t x, real_t y, real_t z)
{
   (void)y;
   switch (kind)
   {
      case AffineFieldKind::X: return x / kL;
      case AffineFieldKind::Z: return z / kL;
      case AffineFieldKind::XPlusZ: return (x + z) / kL;
   }
   MFEM_ABORT("unknown affine field kind");
}

void FillQFromAffineField(const FiniteElementSpace &fes,
                          int component,
                          AffineFieldKind kind,
                          Vector &Q)
{
   const int ndof_total = fes.GetNDofs();
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
   for (int e = 0; e < fes.GetNE(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *tr = fes.GetElementTransformation(e);
      const IntegrationRule &nodes = fe->GetNodes();
      Array<int> edofs;
      fes.GetElementDofs(e, edofs);
      for (int j = 0; j < fe->GetDof(); j++)
      {
         Vector xj(3);
         tr->Transform(nodes.IntPoint(j), xj);
         Q(component * ndof_total + edofs[j]) =
            EvaluateAffineField(kind, xj(0), xj(1), xj(2));
      }
   }
}

Mesh BuildCartesianFaultMesh()
{
   // Phase 2D (2026-04-23): honor SEAS_TEST_FIXTURE env var.  When set
   // to "d4", swap the Kuhn Mesh::MakeCartesian3D(TETRAHEDRON) for the
   // slot-equivariant BuildD4Mesh (see dynamic/d4_tet_mesh.hpp).  This
   // lets the audit produce a Gate 14′ measurement on D4 for §C's
   // decision table.  D4 is only supported at the fixture's native
   // 2×2×2 size; any other dimensions fall back to Kuhn.
   const char *fixture_env = std::getenv("SEAS_TEST_FIXTURE");
   const bool use_d4 = (fixture_env != nullptr
                        && std::string(fixture_env) == "d4");
   if (use_d4)
   {
      if (g_nx == 2 && g_ny == 2 && g_nz == 2)
      {
         std::cout << "  [SEAS_TEST_FIXTURE=d4] using BuildD4Mesh(true).\n";
         Mesh d4 = BuildD4Mesh(/*add_fault=*/true, kL);
         // R4-R001: consolidated invariant gate.  Abort if broken.
         std::cout << "  [R4-R001] audit D4 fixture invariants:\n";
         AssertD4FixtureValid(d4, kL, /*abort_on_fail=*/true);
         return d4;
      }
      std::cout << "  [SEAS_TEST_FIXTURE=d4] requested but mesh dims "
                << g_nx << "x" << g_ny << "x" << g_nz
                << " != 2x2x2; falling back to Kuhn.\n";
   }

   Mesh mesh = Mesh::MakeCartesian3D(g_nx, g_ny, g_nz, Element::TETRAHEDRON,
                                     kL, kL, kL, false);
   mesh.FinalizeTopology();
   mesh.Finalize();

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0.0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy - 0.5 * kL) < 1e-8)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }

   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

void ReadMeshDimsFromEnv()
{
   auto read_dim = [](const char *name, int fallback) {
      const char *env = std::getenv(name);
      if (!env) { return fallback; }
      const int parsed = std::atoi(env);
      if (parsed < 2 || (parsed % 2) != 0) { return fallback; }
      return parsed;
   };
   g_nx = read_dim("SEAS_TEST_FAULT_NX", g_nx);
   g_ny = read_dim("SEAS_TEST_FAULT_NY", g_ny);
   g_nz = read_dim("SEAS_TEST_FAULT_NZ", g_nz);
}

int SetupFault(WaveOperator<Mesh> &wave, Mesh &mesh, int order,
               std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
               std::vector<Vector> &fault_coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   MFEM_VERIFY(shr_faces.Size() == 0, "serial audit fixture should have no shared fault faces");

   int nqp = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault FTR null");
      nqp = IntRules.Get(ftr->GetGeometryType(), 2 * order).GetNPoints();
   }

   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }

   const int n_fault = int_faces.Size() * nqp;
   if (n_fault > 0)
   {
      InitializeFaultDOFs(dof_data, n_fault, fault_coords);
      for (int i = 0; i < n_fault; i++)
      {
         dof_data[i].tau2_nuc = TPV102Params::nuc_dtau;
      }
   }

   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);
   return n_fault;
}

struct FaultQPSample
{
   int face_idx = -1;
   int qp_idx = -1;
   std::array<real_t, NUM_STATE> i_plus_local{};
   std::array<real_t, NUM_STATE> i_minus_local{};
   std::array<real_t, NUM_STATE> q_imp_plus_local{};
   std::array<real_t, NUM_STATE> q_imp_minus_local{};
   std::array<real_t, NUM_STATE> f_h_plus{};
   std::array<real_t, NUM_STATE> f_h_minus{};
};

struct FaultOnlyAudit
{
   std::vector<FaultQPSample> samples;
   Vector rhs_before_minv;
   Vector q_after_minv;
};

struct NonFaultFaceAudit
{
   Vector boundary_rhs_before_minv;
   Vector interior_rhs_before_minv;
   Vector boundary_after_minv;
   Vector interior_after_minv;
};

// Round-13A: fault-face stage-mean audit.  Used to distinguish two
// residual-pepper hypotheses after Round 12 (AVG_TRIAL / AVG_TCORR
// only closed ~40% — the plan reinterprets this as ~40% within-face
// per-QP + ~60% remaining).  The remaining ~60% is either (A)
// between-face mean variation across mirror-related fault faces,
// or (B) asymmetry introduced when equal face means are deposited
// into bulk via shape/weight/orientation.
//
// FaultMeanMode::None   — baseline per-QP stage values.
// FaultMeanMode::Face   — each face's per-QP stage field replaced
//                         by that face's own mean (within-face
//                         uniform, between-face as baseline).
// FaultMeanMode::Global — every face's per-QP stage field replaced
//                         by a single mean over every fault QP
//                         (globally uniform).
//
// FaultMeanSource::Trial — operate on sigma_n_trial / tau1_trial /
//                          tau2_trial, then CompleteFromTrial.
// FaultMeanSource::Tcorr — operate on sigma_n_corr / tau1_corr /
//                          tau2_corr (no downstream recompute).
enum class FaultMeanMode   { None, Face, Global };
enum class FaultMeanSource { Trial, Tcorr };

const char *FaultMeanModeLabel(FaultMeanMode m)
{
   switch (m)
   {
      case FaultMeanMode::None:   return "none";
      case FaultMeanMode::Face:   return "face-mean";
      case FaultMeanMode::Global: return "global-mean";
   }
   return "?";
}

const char *FaultMeanSourceLabel(FaultMeanSource s)
{
   switch (s)
   {
      case FaultMeanSource::Trial: return "trial";
      case FaultMeanSource::Tcorr: return "tcorr";
   }
   return "?";
}

// Per-face aggregate used to detect between-face drift in face-mean
// quantities.  The F_h_* fields capture the mean of the PER-SIDE
// flux carried downstream into the rhs deposit; together with
// tau*_mean these let Gate 4b decide whether face means already
// disagree across orbit-mates.
struct FaultFaceMeanSample
{
   int face_idx = -1;
   bool elem1_on_plus = false;
   real_t cx = 0.0, cy = 0.0, cz = 0.0;

   real_t sigma_n_trial_mean = 0.0;
   real_t tau1_trial_mean    = 0.0;
   real_t tau2_trial_mean    = 0.0;

   real_t sigma_n_corr_mean = 0.0;
   real_t tau1_corr_mean    = 0.0;
   real_t tau2_corr_mean    = 0.0;

   real_t fh_plus_sxy_mean  = 0.0;
   real_t fh_plus_sxz_mean  = 0.0;
   real_t fh_minus_sxy_mean = 0.0;
   real_t fh_minus_sxz_mean = 0.0;
};

struct FaultFaceMeanAudit
{
   std::vector<FaultFaceMeanSample> faces;
   Vector rhs_before_minv;
   Vector rhs_after_minv;
};

// Round-13B: non-fault face seed audit.  Round-13A ruled out
// between-face fault-mean variation at step 1; Round 12 showed that
// within-face fault AVG_TRIAL / AVG_TCORR only closes ~40% of pepper.
// The remaining signal must originate in the STEP-0 non-fault face
// deposit that generates the asymmetric Q_step0 later steps see.
// Round 13B enumerates interior and boundary non-fault face flux
// values, optionally symmetrizes (Gate 5b / Gate 6c style) or
// orbit-averages per class, and measures the resulting
// fault-adjacent sorted-signature drift.
enum class NonFaultSeedMode
{
   Baseline,
   InteriorSym,
   BoundarySym,
   BothSym,
   InteriorFaceMean,
   BoundaryFaceMean,
   BothFaceMean,
   // Round-14A: interior-sym background + variable x-side boundary
   // treatment.  Used by Gate 4g to isolate where the remaining ~1 Pa
   // SXY floor lives after the dominant interior seed is suppressed.
   InteriorSym_XBoundaryRaw, // interior sym, all boundary raw
   InteriorSym_XBoundarySym, // interior sym, x=0 AND x=L boundary sym
   InteriorSym_XMinSym,      // interior sym, only x=0 boundary sym
   InteriorSym_XMaxSym       // interior sym, only x=L boundary sym
};

// Cube-boundary side classification for the 2x2x2 Kuhn/D4 fixture
// (kL = 1000).  None is used for interior non-fault faces.
enum class BoundarySide
{
   None,
   XMin, XMax,
   YMin, YMax,
   ZMin, ZMax
};

struct NonFaultFaceSample
{
   int face_idx = -1;
   int qp_idx   = -1;
   int e1       = -1;
   int e2       = -1;
   int dof_offset1 = -1;
   int dof_offset2 = -1;
   bool is_boundary = false;
   BoundarySide side = BoundarySide::None;

   real_t cx = 0.0, cy = 0.0, cz = 0.0;
   real_t w  = 0.0;
   std::array<real_t, 3> nor{{0, 0, 0}};

   Vector shape1;
   Vector shape2;   // empty on boundary faces

   std::array<real_t, NUM_STATE> fh_raw{};
   std::array<real_t, NUM_STATE> fh_sym{};
};

struct NonFaultSeedAudit
{
   std::vector<NonFaultFaceSample> faces;
   Vector rhs_before_minv;
   Vector rhs_after_minv;
};

// Round-15: x-side boundary trace dump audit.  Round 14 decided
// Branch B — the residual SXY floor on x-aligned outer sides is not
// an n↔-n orientation issue.  Round 15 dumps orbit-paired x-side
// boundary face inputs AND outputs (I_self, bulk_bg_scaled, F_h) so
// the next diagnosis can separate three hypotheses:
//   (a) I_self[SXY] already asymmetric → bug is upstream of BC dispatch
//   (b) I_self matches but F_h[SXY] differs → BC formula path
//   (c) orbit mates take different bdr_attr / bc_type → classification bug
//
// A z-side control family (BoundarySide::ZMin + ZMax) is collected in
// the same gate for contrast: if z-sides are ULP-clean, the residual
// is x-side-specific.
// WaveOperator::FaceBC is private, so the audit declares its own
// mirror enum with identical semantics to the production classes
// used in DispatchBoundaryBC.
enum class AuditFaceBC { Interior, Absorbing, FreeSurface, Fault };

const char *FaceBCLabel(AuditFaceBC bc)
{
   switch (bc)
   {
      case AuditFaceBC::Interior:    return "Interior";
      case AuditFaceBC::Absorbing:   return "Absorbing";
      case AuditFaceBC::FreeSurface: return "FreeSurface";
      case AuditFaceBC::Fault:       return "Fault";
   }
   return "?";
}

struct XSideBoundaryTraceSample
{
   int face_idx = -1;
   int elem     = -1;
   int dof_offset = -1;
   int bdr_attr = 0;
   int q        = -1;

   BoundarySide side = BoundarySide::None;
   AuditFaceBC  bc_type = AuditFaceBC::Absorbing;

   real_t cx = 0.0, cy = 0.0, cz = 0.0;
   real_t w  = 0.0;

   std::array<real_t, 3> nor{{0.0, 0.0, 0.0}};
   std::array<real_t, NUM_STATE> I_self{};
   std::array<real_t, NUM_STATE> bulk_bg_scaled{};
   std::array<real_t, NUM_STATE> F_h{};

   std::vector<real_t> shape1;
};

struct XSideBoundaryTraceAudit
{
   std::vector<XSideBoundaryTraceSample> x_faces;
   std::vector<XSideBoundaryTraceSample> z_faces;
};

// Round-16: boundary trace sampler comparison.  Gate 4k stores
// I_self[SXY] computed two ways from the same I_data:
//   offset: production-style `shape1(i) * I_data(SXY*ndof_total +
//           dof_offset1 + i)` with dof_offset1 = e1 * ndof.
//   edofs:  explicit-VDof style `shape1(i) * I_data(SXY*ndof_total
//           + edofs1[i])` via fes.GetElementDofs(e1, edofs1).
// If the two disagree in any QP, production boundary trace
// addressing is broken.  If they agree, the bug is upstream of
// boundary interpolation.
struct BoundaryTraceSamplerSample
{
   int face_idx   = -1;
   int elem       = -1;
   int q          = -1;
   int dof_offset = -1;

   BoundarySide side = BoundarySide::None;

   real_t cx = 0.0, cy = 0.0, cz = 0.0;
   std::array<real_t, 3> nor{{0.0, 0.0, 0.0}};

   std::vector<int>    edofs1;
   std::vector<real_t> shape1;

   real_t i_self_offset_sxy = 0.0;
   real_t i_self_edofs_sxy  = 0.0;
};

struct BoundaryTraceSamplerAudit
{
   std::vector<BoundaryTraceSamplerSample> x_faces;
   std::vector<BoundaryTraceSamplerSample> z_faces;
};

// Round-16 element-level I_data[SXY] audit.  Samples element DOFs
// directly via fes.GetElementDofs + I_data lookup, computes the
// element-local SXY mean, and groups by canonicalized centroid so
// orbit-mates can be compared.  Answers: is element I_data[SXY]
// already asymmetric BEFORE any face interpolation?
struct ElementIDataSample
{
   int elem = -1;
   real_t cx = 0.0, cy = 0.0, cz = 0.0;

   std::vector<int>    edofs;
   std::vector<real_t> i_sxy_dofs;
   real_t i_sxy_mean = 0.0;
};

struct ElementIDataAudit
{
   std::vector<ElementIDataSample> elems;
};

// Round-17: per-DOF and face-interp canonicalization audit.
// `CanonicalLocalDof` pairs a local DOF's reference-space position
// with its I_data[SXY] value.  `canonical_dofs` is the `raw_dofs`
// vector sorted by reference position so orbit-paired elements can
// be compared without any ordering ambiguity.  If `raw` is dirty
// but `canonical` is clean, the element data is permutation-
// equivalent and the bug is local DOF ordering / equivariance.
// If `canonical` is still dirty, the predictor injects higher-order
// asymmetry that is not just a reorder.
struct CanonicalLocalDof
{
   int local_idx = -1;
   real_t xr = 0.0, yr = 0.0, zr = 0.0;
   real_t value = 0.0;
};

struct ElementSxyDofSample
{
   int elem = -1;
   real_t cx = 0.0, cy = 0.0, cz = 0.0;

   std::vector<CanonicalLocalDof> raw_dofs;
   std::vector<CanonicalLocalDof> canonical_dofs;

   real_t raw_mean = 0.0;
};

struct ElementSxyDofAudit
{
   std::vector<ElementSxyDofSample> elems;
};

struct FaceInterpolationSample
{
   int face_idx = -1;
   int elem = -1;
   int q = -1;
   BoundarySide side = BoundarySide::None;

   real_t cx = 0.0, cy = 0.0, cz = 0.0;
   std::array<real_t, 3> nor{{0.0, 0.0, 0.0}};

   std::vector<CanonicalLocalDof> raw_dofs;
   std::vector<CanonicalLocalDof> canonical_dofs;
   std::vector<real_t> shape_raw;
   std::vector<real_t> shape_canonical;

   real_t interp_raw = 0.0;
   real_t interp_canonical = 0.0;
};

struct FaceInterpolationAudit
{
   std::vector<FaceInterpolationSample> x_faces;
   std::vector<FaceInterpolationSample> z_faces;
};

// Round-18A: face-local canonical interpolation probe.  Stores the
// same face QP's interpolation under three orderings:
//   raw            — current local-DOF order.
//   elem-canonical — Round-17 sort by reference position (3D).
//   face-canonical — sort by reference position restricted to the
//                    sampled face coordinates: (yr, zr) on x-faces,
//                    (xr, yr) on z-faces.  Targets the hypothesis
//                    that face sampling itself must be permutation-
//                    aware relative to face-local coords, not full
//                    element-local coords.
struct FaceLocalCanonicalSample
{
   int face_idx = -1;
   int elem = -1;
   int q = -1;
   BoundarySide side = BoundarySide::None;

   real_t cy = 0.0, cz = 0.0;
   std::array<real_t, 3> nor{{0.0, 0.0, 0.0}};

   std::vector<CanonicalLocalDof> elem_canonical_dofs;
   std::vector<CanonicalLocalDof> face_canonical_dofs;

   std::vector<real_t> shape_raw;
   std::vector<real_t> shape_elem_canonical;
   std::vector<real_t> shape_face_canonical;

   real_t interp_raw = 0.0;
   real_t interp_elem_canonical = 0.0;
   real_t interp_face_canonical = 0.0;
};

struct FaceLocalCanonicalAudit
{
   std::vector<FaceLocalCanonicalSample> x_faces;
   std::vector<FaceLocalCanonicalSample> z_faces;
};

const char *NonFaultSeedModeLabel(NonFaultSeedMode mode)
{
   switch (mode)
   {
      case NonFaultSeedMode::Baseline:                 return "baseline";
      case NonFaultSeedMode::InteriorSym:              return "interiorSym";
      case NonFaultSeedMode::BoundarySym:              return "boundarySym";
      case NonFaultSeedMode::BothSym:                  return "bothSym";
      case NonFaultSeedMode::InteriorFaceMean:         return "interiorMean";
      case NonFaultSeedMode::BoundaryFaceMean:         return "boundaryMean";
      case NonFaultSeedMode::BothFaceMean:             return "bothMean";
      case NonFaultSeedMode::InteriorSym_XBoundaryRaw: return "iSym+xBndRaw";
      case NonFaultSeedMode::InteriorSym_XBoundarySym: return "iSym+xBndSym";
      case NonFaultSeedMode::InteriorSym_XMinSym:      return "iSym+xMinSym";
      case NonFaultSeedMode::InteriorSym_XMaxSym:      return "iSym+xMaxSym";
   }
   return "?";
}

const char *BoundarySideLabel(BoundarySide side)
{
   switch (side)
   {
      case BoundarySide::None: return "interior";
      case BoundarySide::XMin: return "x=0";
      case BoundarySide::XMax: return "x=L";
      case BoundarySide::YMin: return "y=0";
      case BoundarySide::YMax: return "y=L";
      case BoundarySide::ZMin: return "z=0";
      case BoundarySide::ZMax: return "z=L";
   }
   return "?";
}

struct AffineTraceProbeResult
{
   real_t worst_self_err = 0.0;
   real_t worst_nbr_err = 0.0;
   real_t worst_jump_err = 0.0;
   real_t worst_bdry_err = 0.0;
   AffineFieldKind worst_self_field = AffineFieldKind::X;
   AffineFieldKind worst_nbr_field = AffineFieldKind::X;
   AffineFieldKind worst_jump_field = AffineFieldKind::X;
   AffineFieldKind worst_bdry_field = AffineFieldKind::X;
};

struct ElementMeanData
{
   int elem = -1;
   real_t cy = 0.0;
   Array<int> edofs;
   Vector mean_weights;
   real_t volume = 0.0;
};

struct AnalyticTraceLiftResult
{
   real_t worst_iface_drift = 0.0;
   real_t worst_bface_drift = 0.0;
   AffineFieldKind worst_iface_field = AffineFieldKind::X;
   AffineFieldKind worst_bface_field = AffineFieldKind::X;
};

void ResetFaultState(std::vector<DOFData> &dof_data,
                     const std::vector<Vector> &fault_coords)
{
   InitializeFaultDOFs(dof_data, static_cast<int>(dof_data.size()), fault_coords);
   for (auto &d : dof_data) { d.tau2_nuc = TPV102Params::nuc_dtau; }
}

template <typename MeshT>
void ApplyMassInverseManually(const WaveOperator<MeshT> &wave, Vector &rhs)
{
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   for (int e = 0; e < wave.NumElements(); e++)
   {
      Array<int> edofs;
      fes.GetElementDofs(e, edofs);
      Vector elem_rhs(edofs.Size());
      Vector elem_result(edofs.Size());
      const DenseMatrix &minv = wave.GetElementMassInverse(e);
      for (int c = 0; c < NUM_STATE; c++)
      {
         for (int i = 0; i < edofs.Size(); i++)
         {
            elem_rhs(i) = rhs(c * ndof_total + edofs[i]);
         }
         minv.Mult(elem_rhs, elem_result);
         for (int i = 0; i < edofs.Size(); i++)
         {
            rhs(c * ndof_total + edofs[i]) = elem_result(i);
         }
      }
   }
}

ElementMeanData BuildElementMeanData(WaveOperator<Mesh> &wave, Mesh &mesh, int elem)
{
   ElementMeanData data;
   data.elem = elem;

   Array<int> ev;
   mesh.GetElementVertices(elem, ev);
   for (int v = 0; v < ev.Size(); v++) { data.cy += mesh.GetVertex(ev[v])[1]; }
   data.cy /= ev.Size();

   const auto &fes = wave.GetFESpace();
   fes.GetElementDofs(elem, data.edofs);
   const FiniteElement *fe = fes.GetFE(elem);
   const int ndof = fe->GetDof();
   data.mean_weights.SetSize(ndof);
   data.mean_weights = 0.0;

   ElementTransformation *tr = mesh.GetElementTransformation(elem);
   const IntegrationRule &ir =
      IntRules.Get(mesh.GetElementBaseGeometry(elem), 2 * wave.GetOrder() + 2);
   Vector shape(ndof);
   for (int q = 0; q < ir.GetNPoints(); q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      tr->SetIntPoint(&ip);
      fe->CalcShape(ip, shape);
      const real_t w = ip.weight * tr->Weight();
      data.volume += w;
      for (int i = 0; i < ndof; i++)
      {
         data.mean_weights(i) += w * shape(i);
      }
   }
   return data;
}

std::vector<ElementMeanData> BuildFaultAdjacentElementData(WaveOperator<Mesh> &wave,
                                                           Mesh &mesh)
{
   std::vector<int> elems;
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      elems.push_back(ftr->Elem1No);
      elems.push_back(ftr->Elem2No);
   }
   std::sort(elems.begin(), elems.end());
   elems.erase(std::unique(elems.begin(), elems.end()), elems.end());

   std::vector<ElementMeanData> data;
   data.reserve(elems.size());
   for (int elem : elems)
   {
      data.push_back(BuildElementMeanData(wave, mesh, elem));
   }
   return data;
}

real_t ComputeElementMean(const Vector &Q, int comp, int ndof_total,
                          const ElementMeanData &elem_data)
{
   real_t sum = 0.0;
   for (int i = 0; i < elem_data.edofs.Size(); i++)
   {
      sum += elem_data.mean_weights(i) * Q(comp * ndof_total + elem_data.edofs[i]);
   }
   return sum / elem_data.volume;
}

real_t MaxRelativeSpread(const std::vector<FaultQPSample> &samples,
                         const std::array<real_t, NUM_STATE> FaultQPSample::*field,
                         int &worst_comp)
{
   real_t worst = 0.0;
   worst_comp = -1;
   for (int c = 0; c < NUM_STATE; c++)
   {
      real_t min_v = std::numeric_limits<real_t>::max();
      real_t max_v = std::numeric_limits<real_t>::lowest();
      for (const auto &sample : samples)
      {
         const real_t v = (sample.*field)[c];
         min_v = std::min(min_v, v);
         max_v = std::max(max_v, v);
      }
      const real_t mid = 0.5 * (max_v + min_v);
      const real_t spread = (max_v - min_v) /
         std::max(std::abs(mid), real_t(1.0));
      if (spread > worst)
      {
         worst = spread;
         worst_comp = c;
      }
   }
   return worst;
}

real_t MaxFaultStateSpread(const std::vector<DOFData> &dof_data,
                           int which, int &worst_comp)
{
   auto value_of = [&](const DOFData &d) -> real_t {
      switch (which)
      {
         case 0: return d.slip_rate;
         case 1: return d.tau1_corr;
         case 2: return d.tau2_corr;
         case 3: return d.sigma_n_corr;
         default: MFEM_ABORT("bad DOFData field index");
      }
   };

   real_t min_v = std::numeric_limits<real_t>::max();
   real_t max_v = std::numeric_limits<real_t>::lowest();
   for (const auto &d : dof_data)
   {
      const real_t v = value_of(d);
      min_v = std::min(min_v, v);
      max_v = std::max(max_v, v);
   }
   const real_t mid = 0.5 * (max_v + min_v);
   worst_comp = which;
   return (max_v - min_v) / std::max(std::abs(mid), real_t(1.0));
}

real_t CheckElementSideUniformity(const Vector &Q,
                                  const std::vector<ElementMeanData> &elems,
                                  int ndof_total,
                                  int comp,
                                  bool upper_side,
                                  real_t &mean_min,
                                  real_t &mean_max)
{
   mean_min = std::numeric_limits<real_t>::max();
   mean_max = std::numeric_limits<real_t>::lowest();
   for (const auto &elem : elems)
   {
      const bool is_upper = elem.cy > 0.5 * kL;
      if (is_upper != upper_side) { continue; }
      const real_t mean_val = ComputeElementMean(Q, comp, ndof_total, elem);
      mean_min = std::min(mean_min, mean_val);
      mean_max = std::max(mean_max, mean_val);
   }
   const real_t mid = 0.5 * (mean_max + mean_min);
   return (mean_max - mean_min) / std::max(std::abs(mid), real_t(1.0));
}

real_t MaxSortedSignatureDrift(const Vector &Q,
                               const std::vector<ElementMeanData> &elems,
                               int ndof_total,
                               int comp,
                               bool upper_side)
{
   std::vector<std::vector<real_t>> signatures;
   for (const auto &elem : elems)
   {
      const bool is_upper = elem.cy > 0.5 * kL;
      if (is_upper != upper_side) { continue; }
      std::vector<real_t> sig(elem.edofs.Size());
      for (int i = 0; i < elem.edofs.Size(); i++)
      {
         sig[i] = Q(comp * ndof_total + elem.edofs[i]);
      }
      std::sort(sig.begin(), sig.end());
      signatures.push_back(sig);
   }

   if (signatures.size() < 2) { return 0.0; }
   const auto &ref = signatures.front();
   real_t worst = 0.0;
   for (size_t s = 1; s < signatures.size(); s++)
   {
      for (size_t i = 0; i < ref.size(); i++)
      {
         const real_t scale = std::max({std::abs(ref[i]),
                                        std::abs(signatures[s][i]),
                                        real_t(1.0)});
         const real_t drift = std::abs(ref[i] - signatures[s][i]) / scale;
         worst = std::max(worst, drift);
      }
   }
   return worst;
}

// Round-13A: orbit bucketing for fault faces.  For TPV102 the fault
// lies at y = L/2 (fixed under y-mirror), so the relevant symmetries
// within the fault plane are x ↔ (L-x) and z ↔ (L-z).  The key
// canonicalizes both via min(c, L-c) so mirror-paired faces
// (e.g. (cx, cz) and (L-cx, cz)) share a bucket.  `side` = 0/1 keyed
// on cy is retained (always 0 or always 1 for a flat fault at y=L/2),
// per the plan's key shape, but falls out as a singleton for TPV102.
using FaultFaceOrbitKey = std::tuple<long long, long long, int>;

FaultFaceOrbitKey MakeFaultFaceOrbitKey(real_t cx, real_t cy, real_t cz)
{
   const real_t s = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v * s)); };
   const real_t cx_can = std::min(cx, kL - cx);
   const real_t cz_can = std::min(cz, kL - cz);
   const int side = (cy < 0.5 * kL) ? 0 : 1;
   return std::make_tuple(q(cx_can), q(cz_can), side);
}

real_t MaxFaultFaceMeanOrbitDrift(
   const std::vector<FaultFaceMeanSample> &faces,
   std::function<real_t(const FaultFaceMeanSample&)> getter,
   std::string &worst_label)
{
   std::map<FaultFaceOrbitKey, std::vector<real_t>> buckets;
   for (const auto &f : faces)
   {
      const auto key = MakeFaultFaceOrbitKey(f.cx, f.cy, f.cz);
      buckets[key].push_back(getter(f));
   }
   real_t worst = 0.0;
   worst_label.clear();
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      const real_t vmin = *std::min_element(kv.second.begin(), kv.second.end());
      const real_t vmax = *std::max_element(kv.second.begin(), kv.second.end());
      const real_t drift = vmax - vmin;
      if (drift > worst)
      {
         worst = drift;
         std::ostringstream oss;
         oss << "key(cx_q=" << std::get<0>(kv.first)
             << ", cz_q=" << std::get<1>(kv.first)
             << ", side=" << std::get<2>(kv.first)
             << ", n=" << kv.second.size() << ")";
         worst_label = oss.str();
      }
   }
   return worst;
}

// Round-13B: orbit bucketing for non-fault faces.  Boundary faces are
// bucketed by cube side (to prevent mixing x=0 and x=L, say, under
// canonicalized centroid); interior faces use side = 0 (None).
// Centroid is canonicalized via min(c, L-c) on all three axes, and
// the face normal is taken in absolute value so faces with opposite
// outward-normal orientation but otherwise-symmetric centroids share
// an orbit.
using NonFaultFaceOrbitKey =
   std::tuple<int, long long, long long, long long,
              long long, long long, long long>;

BoundarySide ClassifyBoundarySide(real_t cx, real_t cy, real_t cz)
{
   const real_t tol = 1e-8 * kL;
   if (std::abs(cx)        < tol) { return BoundarySide::XMin; }
   if (std::abs(cx - kL)   < tol) { return BoundarySide::XMax; }
   if (std::abs(cy)        < tol) { return BoundarySide::YMin; }
   if (std::abs(cy - kL)   < tol) { return BoundarySide::YMax; }
   if (std::abs(cz)        < tol) { return BoundarySide::ZMin; }
   if (std::abs(cz - kL)   < tol) { return BoundarySide::ZMax; }
   return BoundarySide::None;
}

NonFaultFaceOrbitKey MakeNonFaultFaceOrbitKey(const NonFaultFaceSample &s)
{
   const real_t scale = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v * scale)); };

   const real_t cx_can = std::min(s.cx, kL - s.cx);
   const real_t cy_can = std::min(s.cy, kL - s.cy);
   const real_t cz_can = std::min(s.cz, kL - s.cz);

   return std::make_tuple(
      static_cast<int>(s.side),
      q(cx_can), q(cy_can), q(cz_can),
      q(std::abs(s.nor[0])),
      q(std::abs(s.nor[1])),
      q(std::abs(s.nor[2])));
}

// Round-13C Patch 3: per-boundary-side orbit drift reducer for the
// Gate 4f SXY decomposition.  Filters faces to those with
// `is_boundary && sample.side == side`, buckets via the usual
// `MakeNonFaultFaceOrbitKey`, and returns the worst max-min drift.
real_t MaxBoundarySideOrbitDrift(
   const std::vector<NonFaultFaceSample> &faces,
   BoundarySide side,
   std::function<real_t(const NonFaultFaceSample&)> getter,
   std::string &worst_label)
{
   std::map<NonFaultFaceOrbitKey, std::vector<real_t>> buckets;
   int n_matched = 0;
   for (const auto &f : faces)
   {
      if (!f.is_boundary || f.side != side) { continue; }
      ++n_matched;
      const auto key = MakeNonFaultFaceOrbitKey(f);
      buckets[key].push_back(getter(f));
   }
   real_t worst = 0.0;
   worst_label.clear();
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      const real_t vmin = *std::min_element(kv.second.begin(), kv.second.end());
      const real_t vmax = *std::max_element(kv.second.begin(), kv.second.end());
      const real_t drift = vmax - vmin;
      if (drift > worst)
      {
         worst = drift;
         std::ostringstream oss;
         oss << "n=" << kv.second.size();
         worst_label = oss.str();
      }
   }
   if (worst_label.empty())
   {
      std::ostringstream oss;
      oss << "n_faces=" << n_matched << " (no mirror)";
      worst_label = oss.str();
   }
   return worst;
}

// Round-15 x-side trace orbit key: groups by (side, bdr_attr,
// cy, cz, |n|).  cy and cz are NOT canonicalized because we want
// to discover whether specific (cy, cz) orbit mates on the same
// x-side disagree in their dispatch or outputs.  Within one
// x-side (cx ∈ {0, L}), orbit mates share |cx|, |cz| up to
// face-splitting but the inner pair of triangles on one hex-face
// will land in the same bucket.
using XSideTraceOrbitKey =
   std::tuple<int, int, long long, long long,
              long long, long long, long long>;

XSideTraceOrbitKey MakeBoundaryTraceOrbitKey(const XSideBoundaryTraceSample &s)
{
   const real_t scale = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v * scale)); };

   const real_t cy_can = std::min(s.cy, kL - s.cy);
   const real_t cz_can = std::min(s.cz, kL - s.cz);
   return std::make_tuple(
      static_cast<int>(s.side),
      s.bdr_attr,
      q(cy_can), q(cz_can),
      q(std::abs(s.nor[0])),
      q(std::abs(s.nor[1])),
      q(std::abs(s.nor[2])));
}

template <typename Getter>
real_t MaxTraceOrbitDrift(
   const std::vector<XSideBoundaryTraceSample> &samples,
   Getter getter,
   std::string &worst_label,
   std::vector<XSideBoundaryTraceSample> *worst_bucket = nullptr)
{
   std::map<XSideTraceOrbitKey,
            std::vector<XSideBoundaryTraceSample>> buckets;
   for (const auto &s : samples)
   {
      buckets[MakeBoundaryTraceOrbitKey(s)].push_back(s);
   }
   real_t worst = 0.0;
   worst_label.clear();
   if (worst_bucket) { worst_bucket->clear(); }
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      real_t vmin = getter(kv.second.front());
      real_t vmax = vmin;
      for (const auto &s : kv.second)
      {
         const real_t v = getter(s);
         vmin = std::min(vmin, v);
         vmax = std::max(vmax, v);
      }
      const real_t drift = vmax - vmin;
      if (drift > worst)
      {
         worst = drift;
         std::ostringstream oss;
         oss << "bucket(side=" << std::get<0>(kv.first)
             << ", attr=" << std::get<1>(kv.first)
             << ", cy=" << std::get<2>(kv.first)
             << ", cz=" << std::get<3>(kv.first)
             << ", n=" << kv.second.size() << ")";
         worst_label = oss.str();
         if (worst_bucket) { *worst_bucket = kv.second; }
      }
   }
   return worst;
}

// Round-16: orbit keys / reducers / printers for boundary trace
// sampler comparison and element-level I_data[SXY].
using BoundarySamplerOrbitKey =
   std::tuple<int, long long, long long, long long, long long>;

BoundarySamplerOrbitKey MakeBoundarySamplerOrbitKey(
   const BoundaryTraceSamplerSample &s)
{
   const real_t scale = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v * scale)); };
   const real_t cy_can = std::min(s.cy, kL - s.cy);
   const real_t cz_can = std::min(s.cz, kL - s.cz);
   return std::make_tuple(
      static_cast<int>(s.side),
      q(cy_can), q(cz_can),
      q(std::abs(s.nor[0])),
      q(std::abs(s.nor[2])));
}

using ElementOrbitKey = std::tuple<long long, long long, long long>;

ElementOrbitKey MakeElementOrbitKey(const ElementIDataSample &s)
{
   const real_t scale = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v * scale)); };
   return std::make_tuple(
      q(std::min(s.cx, kL - s.cx)),
      q(std::min(s.cy, kL - s.cy)),
      q(std::min(s.cz, kL - s.cz)));
}

template <typename Getter>
real_t MaxSamplerOrbitDrift(
   const std::vector<BoundaryTraceSamplerSample> &samples,
   Getter getter,
   std::string &worst_label,
   std::vector<BoundaryTraceSamplerSample> *worst_bucket = nullptr)
{
   std::map<BoundarySamplerOrbitKey,
            std::vector<BoundaryTraceSamplerSample>> buckets;
   for (const auto &s : samples)
   {
      buckets[MakeBoundarySamplerOrbitKey(s)].push_back(s);
   }
   real_t worst = 0.0;
   worst_label.clear();
   if (worst_bucket) { worst_bucket->clear(); }
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      real_t vmin = getter(kv.second.front());
      real_t vmax = vmin;
      for (const auto &s : kv.second)
      {
         const real_t v = getter(s);
         vmin = std::min(vmin, v);
         vmax = std::max(vmax, v);
      }
      const real_t drift = vmax - vmin;
      if (drift > worst)
      {
         worst = drift;
         std::ostringstream oss;
         oss << "bucket(side=" << std::get<0>(kv.first)
             << ", cy=" << std::get<1>(kv.first)
             << ", cz=" << std::get<2>(kv.first)
             << ", n=" << kv.second.size() << ")";
         worst_label = oss.str();
         if (worst_bucket) { *worst_bucket = kv.second; }
      }
   }
   return worst;
}

template <typename Getter>
real_t MaxElementOrbitDrift(
   const std::vector<ElementIDataSample> &samples,
   Getter getter,
   std::string &worst_label,
   std::vector<ElementIDataSample> *worst_bucket = nullptr)
{
   std::map<ElementOrbitKey, std::vector<ElementIDataSample>> buckets;
   for (const auto &s : samples)
   {
      buckets[MakeElementOrbitKey(s)].push_back(s);
   }
   real_t worst = 0.0;
   worst_label.clear();
   if (worst_bucket) { worst_bucket->clear(); }
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      real_t vmin = getter(kv.second.front());
      real_t vmax = vmin;
      for (const auto &s : kv.second)
      {
         const real_t v = getter(s);
         vmin = std::min(vmin, v);
         vmax = std::max(vmax, v);
      }
      const real_t drift = vmax - vmin;
      if (drift > worst)
      {
         worst = drift;
         std::ostringstream oss;
         oss << "bucket(cx=" << std::get<0>(kv.first)
             << ", cy=" << std::get<1>(kv.first)
             << ", cz=" << std::get<2>(kv.first)
             << ", n=" << kv.second.size() << ")";
         worst_label = oss.str();
         if (worst_bucket) { *worst_bucket = kv.second; }
      }
   }
   return worst;
}

void PrintBoundarySamplerBucket(
   const std::vector<BoundaryTraceSamplerSample> &bucket,
   const char *label)
{
   if (bucket.empty()) { return; }
   std::cout << "      " << label << "  (bucket size=" << bucket.size() << ")\n";
   for (const auto &s : bucket)
   {
      const real_t delta = s.i_self_offset_sxy - s.i_self_edofs_sxy;
      std::cout << "        face=" << s.face_idx
                << " elem=" << s.elem
                << " side=" << BoundarySideLabel(s.side)
                << " q=" << s.q
                << std::fixed << std::setprecision(2)
                << " cy=" << s.cy << " cz=" << s.cz
                << " nor=(" << std::showpos << s.nor[0]
                << "," << s.nor[1] << "," << s.nor[2] << ")"
                << std::noshowpos
                << std::scientific << std::setprecision(3)
                << " offset=" << s.i_self_offset_sxy
                << " edofs=" << s.i_self_edofs_sxy
                << " delta=" << delta
                << "\n";
   }
}

// Round-17: reference-space DOF position extraction and canonical
// ordering helpers.
std::vector<std::array<real_t, 3>> GetElementRefDofPositions(
   const FiniteElement &fe)
{
   const IntegrationRule &nodes = fe.GetNodes();
   const int ndof = fe.GetDof();
   std::vector<std::array<real_t, 3>> out(ndof);
   for (int i = 0; i < ndof; i++)
   {
      const IntegrationPoint &ip = nodes.IntPoint(i);
      out[i][0] = ip.x;
      out[i][1] = ip.y;
      out[i][2] = ip.z;
   }
   return out;
}

using RefPosKey = std::tuple<long long, long long, long long>;

RefPosKey MakeRefPosKey(real_t xr, real_t yr, real_t zr)
{
   const real_t scale = 1e9;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v * scale)); };
   return std::make_tuple(q(xr), q(yr), q(zr));
}

std::vector<CanonicalLocalDof> BuildRawLocalDofs(
   const FiniteElement &fe,
   const Array<int> &edofs,
   const Vector &I_data,
   int ndof_total)
{
   const auto pos = GetElementRefDofPositions(fe);
   const int ndof = fe.GetDof();
   std::vector<CanonicalLocalDof> out(ndof);
   for (int i = 0; i < ndof; i++)
   {
      out[i].local_idx = i;
      out[i].xr = pos[i][0];
      out[i].yr = pos[i][1];
      out[i].zr = pos[i][2];
      out[i].value = I_data(SXY * ndof_total + edofs[i]);
   }
   return out;
}

std::vector<CanonicalLocalDof> BuildCanonicalLocalDofs(
   const FiniteElement &fe,
   const Array<int> &edofs,
   const Vector &I_data,
   int ndof_total)
{
   auto out = BuildRawLocalDofs(fe, edofs, I_data, ndof_total);
   std::sort(out.begin(), out.end(),
             [](const CanonicalLocalDof &a, const CanonicalLocalDof &b)
             {
                return MakeRefPosKey(a.xr, a.yr, a.zr) <
                       MakeRefPosKey(b.xr, b.yr, b.zr);
             });
   return out;
}

// Compare element samples within their shared orbit bucket (keyed by
// canonicalized centroid).  Sample-pair drift is the max |Δ| over
// equal-index DOF comparisons; bucket drift is the max over all
// sample pairs in the bucket.  Requires all samples in a bucket to
// have the same DOF count.
real_t MaxElementRawDofOrbitDrift(
   const std::vector<ElementSxyDofSample> &samples,
   std::string &worst_label,
   std::vector<ElementSxyDofSample> *worst_bucket = nullptr)
{
   std::map<ElementOrbitKey, std::vector<ElementSxyDofSample>> buckets;
   for (const auto &s : samples)
   {
      const auto key = std::make_tuple(
         static_cast<long long>(std::round(std::min(s.cx, kL - s.cx) * 1e6)),
         static_cast<long long>(std::round(std::min(s.cy, kL - s.cy) * 1e6)),
         static_cast<long long>(std::round(std::min(s.cz, kL - s.cz) * 1e6)));
      buckets[key].push_back(s);
   }
   real_t worst = 0.0;
   worst_label.clear();
   if (worst_bucket) { worst_bucket->clear(); }
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      const int ndof = static_cast<int>(kv.second.front().raw_dofs.size());
      bool same_ndof = true;
      for (const auto &s : kv.second)
      {
         if (static_cast<int>(s.raw_dofs.size()) != ndof) { same_ndof = false; break; }
      }
      if (!same_ndof) { continue; }
      real_t bucket_worst = 0.0;
      for (int i = 0; i < ndof; i++)
      {
         real_t vmin = kv.second.front().raw_dofs[i].value;
         real_t vmax = vmin;
         for (const auto &s : kv.second)
         {
            vmin = std::min(vmin, s.raw_dofs[i].value);
            vmax = std::max(vmax, s.raw_dofs[i].value);
         }
         bucket_worst = std::max(bucket_worst, vmax - vmin);
      }
      if (bucket_worst > worst)
      {
         worst = bucket_worst;
         std::ostringstream oss;
         oss << "bucket(cx=" << std::get<0>(kv.first)
             << ", cy=" << std::get<1>(kv.first)
             << ", cz=" << std::get<2>(kv.first)
             << ", n=" << kv.second.size() << ")";
         worst_label = oss.str();
         if (worst_bucket) { *worst_bucket = kv.second; }
      }
   }
   return worst;
}

real_t MaxElementCanonicalDofOrbitDrift(
   const std::vector<ElementSxyDofSample> &samples,
   std::string &worst_label,
   std::vector<ElementSxyDofSample> *worst_bucket = nullptr)
{
   std::map<ElementOrbitKey, std::vector<ElementSxyDofSample>> buckets;
   for (const auto &s : samples)
   {
      const auto key = std::make_tuple(
         static_cast<long long>(std::round(std::min(s.cx, kL - s.cx) * 1e6)),
         static_cast<long long>(std::round(std::min(s.cy, kL - s.cy) * 1e6)),
         static_cast<long long>(std::round(std::min(s.cz, kL - s.cz) * 1e6)));
      buckets[key].push_back(s);
   }
   real_t worst = 0.0;
   worst_label.clear();
   if (worst_bucket) { worst_bucket->clear(); }
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      const int ndof = static_cast<int>(kv.second.front().canonical_dofs.size());
      bool same_ndof = true;
      for (const auto &s : kv.second)
      {
         if (static_cast<int>(s.canonical_dofs.size()) != ndof) { same_ndof = false; break; }
      }
      if (!same_ndof) { continue; }
      real_t bucket_worst = 0.0;
      for (int i = 0; i < ndof; i++)
      {
         real_t vmin = kv.second.front().canonical_dofs[i].value;
         real_t vmax = vmin;
         for (const auto &s : kv.second)
         {
            vmin = std::min(vmin, s.canonical_dofs[i].value);
            vmax = std::max(vmax, s.canonical_dofs[i].value);
         }
         bucket_worst = std::max(bucket_worst, vmax - vmin);
      }
      if (bucket_worst > worst)
      {
         worst = bucket_worst;
         std::ostringstream oss;
         oss << "bucket(cx=" << std::get<0>(kv.first)
             << ", cy=" << std::get<1>(kv.first)
             << ", cz=" << std::get<2>(kv.first)
             << ", n=" << kv.second.size() << ")";
         worst_label = oss.str();
         if (worst_bucket) { *worst_bucket = kv.second; }
      }
   }
   return worst;
}

void PrintElementSxyDofBucket(
   const std::vector<ElementSxyDofSample> &bucket,
   const char *label, bool canonical)
{
   if (bucket.empty()) { return; }
   std::cout << "      " << label
             << "  (bucket size=" << bucket.size()
             << ", order=" << (canonical ? "canonical" : "raw") << ")\n";
   for (const auto &s : bucket)
   {
      const auto &dofs = canonical ? s.canonical_dofs : s.raw_dofs;
      std::cout << "        elem=" << s.elem
                << std::fixed << std::setprecision(2)
                << " c=(" << s.cx << "," << s.cy << "," << s.cz << ")"
                << std::scientific << std::setprecision(3)
                << " dofs[0..4)={";
      const int n_show = std::min<int>(4, static_cast<int>(dofs.size()));
      for (int i = 0; i < n_show; i++)
      {
         if (i) { std::cout << ","; }
         std::cout << dofs[i].value;
      }
      std::cout << "}\n";
   }
}

// Face-interpolation orbit key: group samples by (side, cy, cz, |n|, q)
// so orbit mates with the same reference-QP position on mirror-paired
// faces share a bucket.
using FaceInterpOrbitKey =
   std::tuple<int, long long, long long, long long, long long, int>;

FaceInterpOrbitKey MakeFaceInterpOrbitKey(const FaceInterpolationSample &s)
{
   const real_t scale = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v * scale)); };
   const real_t cy_can = std::min(s.cy, kL - s.cy);
   const real_t cz_can = std::min(s.cz, kL - s.cz);
   return std::make_tuple(
      static_cast<int>(s.side),
      q(cy_can), q(cz_can),
      q(std::abs(s.nor[0])),
      q(std::abs(s.nor[2])),
      s.q);
}

template <typename Getter>
real_t MaxFaceInterpolationOrbitDrift(
   const std::vector<FaceInterpolationSample> &samples,
   Getter getter,
   std::string &worst_label,
   std::vector<FaceInterpolationSample> *worst_bucket = nullptr)
{
   std::map<FaceInterpOrbitKey, std::vector<FaceInterpolationSample>> buckets;
   for (const auto &s : samples)
   {
      buckets[MakeFaceInterpOrbitKey(s)].push_back(s);
   }
   real_t worst = 0.0;
   worst_label.clear();
   if (worst_bucket) { worst_bucket->clear(); }
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      real_t vmin = getter(kv.second.front());
      real_t vmax = vmin;
      for (const auto &s : kv.second)
      {
         const real_t v = getter(s);
         vmin = std::min(vmin, v);
         vmax = std::max(vmax, v);
      }
      const real_t drift = vmax - vmin;
      if (drift > worst)
      {
         worst = drift;
         std::ostringstream oss;
         oss << "bucket(side=" << std::get<0>(kv.first)
             << ", cy=" << std::get<1>(kv.first)
             << ", cz=" << std::get<2>(kv.first)
             << ", q=" << std::get<5>(kv.first)
             << ", n=" << kv.second.size() << ")";
         worst_label = oss.str();
         if (worst_bucket) { *worst_bucket = kv.second; }
      }
   }
   return worst;
}

// Round-18A: face-local canonical DOF builder.  Sorts DOFs by their
// reference coordinates restricted to the face plane — (yr, zr) for
// x-normal faces, (xr, yr) for z-normal faces — and returns both the
// sorted DOFs and the permutation (canonical_k → raw_index).
std::vector<CanonicalLocalDof> BuildFaceLocalCanonicalLocalDofs(
   const FiniteElement &fe,
   const Array<int> &edofs,
   const Vector &I_data,
   int ndof_total,
   BoundarySide side,
   std::vector<int> *perm_out = nullptr)
{
   auto out = BuildRawLocalDofs(fe, edofs, I_data, ndof_total);
   const bool x_face = (side == BoundarySide::XMin || side == BoundarySide::XMax);
   const bool z_face = (side == BoundarySide::ZMin || side == BoundarySide::ZMax);

   auto key_of = [&](const CanonicalLocalDof &d) -> RefPosKey
   {
      if (x_face)      { return MakeRefPosKey(d.yr, d.zr, d.xr); }
      else if (z_face) { return MakeRefPosKey(d.xr, d.yr, d.zr); }
      return MakeRefPosKey(d.xr, d.yr, d.zr);
   };
   std::sort(out.begin(), out.end(),
             [&](const CanonicalLocalDof &a, const CanonicalLocalDof &b)
             { return key_of(a) < key_of(b); });
   if (perm_out)
   {
      perm_out->resize(out.size());
      for (size_t k = 0; k < out.size(); k++)
      { (*perm_out)[k] = out[k].local_idx; }
   }
   return out;
}

template <typename Getter>
real_t MaxFaceLocalCanonicalOrbitDrift(
   const std::vector<FaceLocalCanonicalSample> &samples,
   Getter getter,
   std::string &worst_label,
   std::vector<FaceLocalCanonicalSample> *worst_bucket = nullptr)
{
   // Reuse FaceInterpOrbitKey shape by packing the sample into a
   // temporary FaceInterpolationSample-like key.
   using KT = std::tuple<int, long long, long long, long long, int>;
   std::map<KT, std::vector<FaceLocalCanonicalSample>> buckets;
   for (const auto &s : samples)
   {
      const real_t scale = 1e6;
      auto q = [&](real_t v) { return static_cast<long long>(std::round(v * scale)); };
      const real_t cy_can = std::min(s.cy, kL - s.cy);
      const real_t cz_can = std::min(s.cz, kL - s.cz);
      KT key = std::make_tuple(
         static_cast<int>(s.side),
         q(cy_can), q(cz_can),
         q(std::abs(s.nor[0])),
         s.q);
      buckets[key].push_back(s);
   }
   real_t worst = 0.0;
   worst_label.clear();
   if (worst_bucket) { worst_bucket->clear(); }
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      real_t vmin = getter(kv.second.front());
      real_t vmax = vmin;
      for (const auto &s : kv.second)
      {
         const real_t v = getter(s);
         vmin = std::min(vmin, v);
         vmax = std::max(vmax, v);
      }
      const real_t drift = vmax - vmin;
      if (drift > worst)
      {
         worst = drift;
         std::ostringstream oss;
         oss << "bucket(side=" << std::get<0>(kv.first)
             << ", cy=" << std::get<1>(kv.first)
             << ", cz=" << std::get<2>(kv.first)
             << ", q=" << std::get<4>(kv.first)
             << ", n=" << kv.second.size() << ")";
         worst_label = oss.str();
         if (worst_bucket) { *worst_bucket = kv.second; }
      }
   }
   return worst;
}

void PrintFaceLocalCanonicalBucket(
   const std::vector<FaceLocalCanonicalSample> &bucket,
   const char *label, const char *which)
{
   if (bucket.empty()) { return; }
   std::cout << "      " << label << "  (bucket size=" << bucket.size()
             << ", which=" << which << ")\n";
   for (const auto &s : bucket)
   {
      std::cout << "        face=" << s.face_idx
                << " elem=" << s.elem
                << " side=" << BoundarySideLabel(s.side)
                << " q=" << s.q
                << std::fixed << std::setprecision(2)
                << " cy=" << s.cy << " cz=" << s.cz
                << std::scientific << std::setprecision(3)
                << " raw=" << s.interp_raw
                << " eCan=" << s.interp_elem_canonical
                << " fCan=" << s.interp_face_canonical
                << "\n";
   }
}

void PrintFaceInterpolationBucket(
   const std::vector<FaceInterpolationSample> &bucket,
   const char *label, bool canonical)
{
   if (bucket.empty()) { return; }
   std::cout << "      " << label
             << "  (bucket size=" << bucket.size()
             << ", order=" << (canonical ? "canonical" : "raw") << ")\n";
   for (const auto &s : bucket)
   {
      const auto &dofs  = canonical ? s.canonical_dofs : s.raw_dofs;
      const auto &shape = canonical ? s.shape_canonical : s.shape_raw;
      std::cout << "        face=" << s.face_idx
                << " elem=" << s.elem
                << " side=" << BoundarySideLabel(s.side)
                << " q=" << s.q
                << std::fixed << std::setprecision(2)
                << " cy=" << s.cy << " cz=" << s.cz
                << std::scientific << std::setprecision(3)
                << " interp=" << (canonical ? s.interp_canonical : s.interp_raw)
                << " dofs[0..4)={";
      const int n_show = std::min<int>(4, static_cast<int>(dofs.size()));
      for (int i = 0; i < n_show; i++)
      {
         if (i) { std::cout << ","; }
         std::cout << dofs[i].value;
      }
      std::cout << "}";
      std::cout << " shape[0..4)={";
      const int n_show_s = std::min<int>(4, static_cast<int>(shape.size()));
      for (int i = 0; i < n_show_s; i++)
      {
         if (i) { std::cout << ","; }
         std::cout << shape[i];
      }
      std::cout << "}\n";
   }
}

void PrintElementIDataBucket(
   const std::vector<ElementIDataSample> &bucket,
   const char *label)
{
   if (bucket.empty()) { return; }
   std::cout << "      " << label << "  (bucket size=" << bucket.size() << ")\n";
   for (const auto &s : bucket)
   {
      std::cout << "        elem=" << s.elem
                << std::fixed << std::setprecision(2)
                << " c=(" << s.cx << "," << s.cy << "," << s.cz << ")"
                << std::scientific << std::setprecision(3)
                << " mean=" << s.i_sxy_mean;
      const int n_show = std::min<int>(4, static_cast<int>(s.i_sxy_dofs.size()));
      std::cout << " dofs[0.." << n_show << ")={";
      for (int i = 0; i < n_show; i++)
      {
         if (i) { std::cout << ","; }
         std::cout << s.i_sxy_dofs[i];
      }
      std::cout << "}\n";
   }
}

void PrintTraceBucket(const std::vector<XSideBoundaryTraceSample> &bucket,
                      const char *label)
{
   if (bucket.empty()) { return; }
   std::cout << "      " << label << "  (bucket size=" << bucket.size() << ")\n";
   for (const auto &s : bucket)
   {
      std::cout << "        face=" << s.face_idx
                << " side=" << BoundarySideLabel(s.side)
                << " attr=" << s.bdr_attr
                << " bc=" << FaceBCLabel(s.bc_type)
                << " q=" << s.q
                << std::fixed << std::setprecision(2)
                << " cy=" << s.cy << " cz=" << s.cz
                << " nor=(" << std::showpos << s.nor[0]
                << "," << s.nor[1] << "," << s.nor[2] << ")"
                << std::noshowpos
                << std::scientific << std::setprecision(3)
                << " I_self[SXY]=" << s.I_self[SXY]
                << " bulk_bg[SXY]=" << s.bulk_bg_scaled[SXY]
                << " F_h[SXY]=" << s.F_h[SXY]
                << "\n";
   }
}

real_t MaxNonFaultFaceOrbitDrift(
   const std::vector<NonFaultFaceSample> &faces,
   std::function<real_t(const NonFaultFaceSample&)> getter,
   bool boundary_only, bool interior_only,
   std::string &worst_label)
{
   std::map<NonFaultFaceOrbitKey, std::vector<real_t>> buckets;
   for (const auto &f : faces)
   {
      if (boundary_only && !f.is_boundary)  { continue; }
      if (interior_only &&  f.is_boundary)  { continue; }
      const auto key = MakeNonFaultFaceOrbitKey(f);
      buckets[key].push_back(getter(f));
   }
   real_t worst = 0.0;
   worst_label.clear();
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      const real_t vmin = *std::min_element(kv.second.begin(), kv.second.end());
      const real_t vmax = *std::max_element(kv.second.begin(), kv.second.end());
      const real_t drift = vmax - vmin;
      if (drift > worst)
      {
         worst = drift;
         std::ostringstream oss;
         oss << "key(side=" << std::get<0>(kv.first)
             << ", cx=" << std::get<1>(kv.first)
             << ", cy=" << std::get<2>(kv.first)
             << ", cz=" << std::get<3>(kv.first)
             << ", n=" << kv.second.size() << ")";
         worst_label = oss.str();
      }
   }
   return worst;
}

// Round-13B flux builders — factor the Gate 5 / Gate 6 formulas into
// reusable helpers so the seed audit never re-derives them.  Matches:
//   - Interior Raw = flux.Interior(nor, I_self, I_nbr, F_h)
//   - Interior Sym = 0.5 * (Interior(n,L,R) + Interior(-n,R,L))
//   - Boundary Raw = existing RunADERNonFaultFaceAudit dispatch
//                    (absorbing → AbsorbingTotal; natural →
//                     FreeSurface[Godunov]Total; fallback → AbsorbingTotal)
//   - Boundary Sym = 0.5 * (Gamma(n) + Gamma(-n)) where Gamma is the
//                    same BC family as the raw path.
void ComputeInteriorFhRaw(const GodunovFlux &flux,
                           const real_t *nor,
                           const real_t *I_self,
                           const real_t *I_nbr,
                           real_t *F_h)
{
   flux.Interior(nor, I_self, I_nbr, F_h);
}

void ComputeInteriorFhSym(const GodunovFlux &flux,
                           const real_t *nor,
                           const real_t *I_self,
                           const real_t *I_nbr,
                           real_t *F_h)
{
   real_t F_fwd[NUM_STATE], F_rev[NUM_STATE];
   flux.Interior(nor, I_self, I_nbr, F_fwd);
   const real_t nor_rev[3] = {-nor[0], -nor[1], -nor[2]};
   flux.Interior(nor_rev, I_nbr, I_self, F_rev);
   for (int c = 0; c < NUM_STATE; c++)
   {
      F_h[c] = 0.5 * (F_fwd[c] + F_rev[c]);
   }
}

void DispatchBoundaryBC(const WaveOperator<Mesh> &wave,
                        const BoundaryConfig &bc,
                        const real_t *nor,
                        const real_t *I_self,
                        const real_t *bulk_bg_scaled,
                        real_t *F_h)
{
   // Mirrors the RunADERNonFaultFaceAudit dispatch; the fixture uses
   // bc.natural_attrs = {1}, so this routes via FreeSurface[Godunov]Total
   // in practice but absorbing / fallback branches are preserved so the
   // helper stays reusable.
   const int boundary_attr = 1;
   if (bc.absorbing_attrs.count(boundary_attr))
   {
      wave.GetFlux().AbsorbingTotal(nor, I_self, bulk_bg_scaled, F_h);
   }
   else if (bc.natural_attrs.count(boundary_attr))
   {
      if (wave.GetFreeSurfaceBCMode() == FreeSurfaceBCMode::Godunov)
      {
         wave.GetFlux().FreeSurfaceGodunovTotal(nor, I_self, bulk_bg_scaled, F_h);
      }
      else
      {
         wave.GetFlux().FreeSurfaceTotal(nor, I_self, bulk_bg_scaled, F_h);
      }
   }
   else
   {
      wave.GetFlux().AbsorbingTotal(nor, I_self, bulk_bg_scaled, F_h);
   }
}

void ComputeBoundaryFhRaw(const WaveOperator<Mesh> &wave,
                           const BoundaryConfig &bc,
                           const real_t *nor,
                           const real_t *I_self,
                           const real_t *bulk_bg_scaled,
                           real_t *F_h)
{
   DispatchBoundaryBC(wave, bc, nor, I_self, bulk_bg_scaled, F_h);
}

void ComputeBoundaryFhSym(const WaveOperator<Mesh> &wave,
                           const BoundaryConfig &bc,
                           const real_t *nor,
                           const real_t *I_self,
                           const real_t *bulk_bg_scaled,
                           real_t *F_h)
{
   real_t F_fwd[NUM_STATE], F_rev[NUM_STATE];
   DispatchBoundaryBC(wave, bc, nor, I_self, bulk_bg_scaled, F_fwd);
   const real_t nor_rev[3] = {-nor[0], -nor[1], -nor[2]};
   DispatchBoundaryBC(wave, bc, nor_rev, I_self, bulk_bg_scaled, F_rev);
   for (int c = 0; c < NUM_STATE; c++)
   {
      F_h[c] = 0.5 * (F_fwd[c] + F_rev[c]);
   }
}

FaultOnlyAudit RunFaultOnlyAudit(WaveOperator<Mesh> &wave, Mesh &mesh,
                                 FaultFaceFlux &ff,
                                 std::vector<DOFData> &dof_data, real_t dt)
{
   FaultOnlyAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();
   audit.rhs_before_minv.SetSize(size);
   audit.rhs_before_minv = 0.0;

   const FaultBasis *fault_basis = wave.GetFaultBasis();
   MFEM_VERIFY(fault_basis != nullptr, "fault basis must be populated");
   const auto &flux = wave.GetFlux();
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();

   real_t zero_state[NUM_STATE] = {0.0};

   for (int fi = 0; fi < int_faces.Size(); fi++)
   {
      const int face_idx = int_faces[fi];
      auto *ftr = mesh.GetInteriorFaceTransformations(face_idx);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2 * kOrder);

      Array<int> edofs1, edofs2;
      fes.GetElementDofs(ftr->Elem1No, edofs1);
      fes.GetElementDofs(ftr->Elem2No, edofs2);
      const FiniteElement *fe1 = fes.GetFE(ftr->Elem1No);
      const FiniteElement *fe2 = fes.GetFE(ftr->Elem2No);

      const int fb_idx = wave.LookupInteriorFaultBasisIndex(face_idx);
      MFEM_VERIFY(fb_idx >= 0, "missing FaultBasis entry for interior fault face");
      const FaultBasisData &bd = fault_basis->GetBasis(fb_idx);
      MFEM_VERIFY(static_cast<int>(bd.qp_data.size()) == ir.GetNPoints(),
                  "FaultBasis qp count mismatch");

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector n_raw(3);
         CalcOrtho(ftr->Face->Jacobian(), n_raw);
         const real_t w = ip.weight * n_raw.Norml2();

         IntegrationPoint ip1, ip2;
         ftr->Loc1.Transform(ip, ip1);
         ftr->Loc2.Transform(ip, ip2);
         Vector shape1(fe1->GetDof()), shape2(fe2->GetDof());
         fe1->CalcShape(ip1, shape1);
         fe2->CalcShape(ip2, shape2);

         const FaultBasisQPData &qpd = bd.qp_data[q];
         real_t can_n[3], can_t1[3], can_t2[3];
         for (int d = 0; d < 3; d++)
         {
            can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]   : qpd.normal[d];
            can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d] : qpd.tangent1[d];
            can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d] : qpd.tangent2[d];
         }

         DenseMatrix T_can(NUM_STATE), Tinv_can(NUM_STATE);
         GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);
         GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv_can);

         const bool elem1_on_plus = !qpd.sign_flipped;
         const int dof_idx = fi * ir.GetNPoints() + q;
         DOFData &fdata = dof_data[dof_idx];

         real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
         ff.EvaluateADER(fdata, zero_state, zero_state, dt,
                         I_imp_plus, I_imp_minus);

         FaultQPSample sample;
         sample.face_idx = face_idx;
         sample.qp_idx = q;
         for (int c = 0; c < NUM_STATE; c++)
         {
            sample.q_imp_plus_local[c] = I_imp_plus[c] / dt;
            sample.q_imp_minus_local[c] = I_imp_minus[c] / dt;
         }

         real_t I_imp_plus_g[NUM_STATE], I_imp_minus_g[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_imp_plus_g[c] = 0.0;
            I_imp_minus_g[c] = 0.0;
            for (int k = 0; k < NUM_STATE; k++)
            {
               I_imp_plus_g[c] += T_can(c, k) * I_imp_plus[k];
               I_imp_minus_g[c] += T_can(c, k) * I_imp_minus[k];
            }
         }

         real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
         flux.Interior(can_n, I_imp_plus_g, I_imp_plus_g, F_h_plus);
         flux.Interior(can_n, I_imp_minus_g, I_imp_minus_g, F_h_minus);

         for (int c = 0; c < NUM_STATE; c++)
         {
            sample.f_h_plus[c] = F_h_plus[c] / dt;
            sample.f_h_minus[c] = F_h_minus[c] / dt;
         }
         audit.samples.push_back(sample);

         if (elem1_on_plus)
         {
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < edofs1.Size(); i++)
               {
                  audit.rhs_before_minv(c * ndof_total + edofs1[i]) -=
                     w * shape1(i) * F_h_plus[c];
               }
               for (int i = 0; i < edofs2.Size(); i++)
               {
                  audit.rhs_before_minv(c * ndof_total + edofs2[i]) +=
                     w * shape2(i) * F_h_minus[c];
               }
            }
         }
         else
         {
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < edofs1.Size(); i++)
               {
                  audit.rhs_before_minv(c * ndof_total + edofs1[i]) +=
                     w * shape1(i) * F_h_minus[c];
               }
               for (int i = 0; i < edofs2.Size(); i++)
               {
                  audit.rhs_before_minv(c * ndof_total + edofs2[i]) -=
                     w * shape2(i) * F_h_plus[c];
               }
            }
         }
      }
   }

   audit.q_after_minv = audit.rhs_before_minv;
   ApplyMassInverseManually(wave, audit.q_after_minv);
   return audit;
}

FaultOnlyAudit RunADERTraceAudit(WaveOperator<Mesh> &wave, Mesh &mesh,
                                 FaultFaceFlux &ff,
                                 std::vector<DOFData> &dof_data,
                                 const Vector &I_data,
                                 real_t dt)
{
   FaultOnlyAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();
   audit.rhs_before_minv.SetSize(size);
   audit.rhs_before_minv = 0.0;

   const FaultBasis *fault_basis = wave.GetFaultBasis();
   MFEM_VERIFY(fault_basis != nullptr, "fault basis must be populated");
   const auto &flux = wave.GetFlux();
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();

   for (int fi = 0; fi < int_faces.Size(); fi++)
   {
      const int face_idx = int_faces[fi];
      auto *ftr = mesh.GetInteriorFaceTransformations(face_idx);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2 * kOrder);

      Array<int> edofs1, edofs2;
      fes.GetElementDofs(ftr->Elem1No, edofs1);
      fes.GetElementDofs(ftr->Elem2No, edofs2);
      const FiniteElement *fe1 = fes.GetFE(ftr->Elem1No);
      const FiniteElement *fe2 = fes.GetFE(ftr->Elem2No);

      const int fb_idx = wave.LookupInteriorFaultBasisIndex(face_idx);
      MFEM_VERIFY(fb_idx >= 0, "missing FaultBasis entry for interior fault face");
      const FaultBasisData &bd = fault_basis->GetBasis(fb_idx);
      MFEM_VERIFY(static_cast<int>(bd.qp_data.size()) == ir.GetNPoints(),
                  "FaultBasis qp count mismatch");

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector n_raw(3);
         CalcOrtho(ftr->Face->Jacobian(), n_raw);
         const real_t w = ip.weight * n_raw.Norml2();

         IntegrationPoint ip1, ip2;
         ftr->Loc1.Transform(ip, ip1);
         ftr->Loc2.Transform(ip, ip2);
         Vector shape1(fe1->GetDof()), shape2(fe2->GetDof());
         fe1->CalcShape(ip1, shape1);
         fe2->CalcShape(ip2, shape2);

         real_t I_self[NUM_STATE], I_nbr[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self[c] = 0.0;
            I_nbr[c] = 0.0;
            for (int i = 0; i < edofs1.Size(); i++)
            {
               I_self[c] += shape1(i) * I_data(c * ndof_total + edofs1[i]);
            }
            for (int i = 0; i < edofs2.Size(); i++)
            {
               I_nbr[c] += shape2(i) * I_data(c * ndof_total + edofs2[i]);
            }
         }

         const FaultBasisQPData &qpd = bd.qp_data[q];
         real_t can_n[3], can_t1[3], can_t2[3];
         for (int d = 0; d < 3; d++)
         {
            can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]   : qpd.normal[d];
            can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d] : qpd.tangent1[d];
            can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d] : qpd.tangent2[d];
         }

         DenseMatrix T_can(NUM_STATE), Tinv_can(NUM_STATE);
         GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);
         GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv_can);

         real_t I_self_can[NUM_STATE], I_nbr_can[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self_can[c] = 0.0;
            I_nbr_can[c] = 0.0;
            for (int k = 0; k < NUM_STATE; k++)
            {
               I_self_can[c] += Tinv_can(c, k) * I_self[k];
               I_nbr_can[c] += Tinv_can(c, k) * I_nbr[k];
            }
         }

         const bool elem1_on_plus = !qpd.sign_flipped;
         const real_t *I_plus_local = elem1_on_plus ? I_self_can : I_nbr_can;
         const real_t *I_minus_local = elem1_on_plus ? I_nbr_can : I_self_can;

         const int dof_idx = fi * ir.GetNPoints() + q;
         DOFData &fdata = dof_data[dof_idx];

         real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
         ff.EvaluateADER(fdata, I_plus_local, I_minus_local, dt,
                         I_imp_plus, I_imp_minus);

         FaultQPSample sample;
         sample.face_idx = face_idx;
         sample.qp_idx = q;
         for (int c = 0; c < NUM_STATE; c++)
         {
            sample.i_plus_local[c] = I_plus_local[c] / dt;
            sample.i_minus_local[c] = I_minus_local[c] / dt;
            sample.q_imp_plus_local[c] = I_imp_plus[c] / dt;
            sample.q_imp_minus_local[c] = I_imp_minus[c] / dt;
         }

         real_t I_imp_plus_g[NUM_STATE], I_imp_minus_g[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_imp_plus_g[c] = 0.0;
            I_imp_minus_g[c] = 0.0;
            for (int k = 0; k < NUM_STATE; k++)
            {
               I_imp_plus_g[c] += T_can(c, k) * I_imp_plus[k];
               I_imp_minus_g[c] += T_can(c, k) * I_imp_minus[k];
            }
         }

         real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
         flux.Interior(can_n, I_imp_plus_g, I_imp_plus_g, F_h_plus);
         flux.Interior(can_n, I_imp_minus_g, I_imp_minus_g, F_h_minus);

         for (int c = 0; c < NUM_STATE; c++)
         {
            sample.f_h_plus[c] = F_h_plus[c] / dt;
            sample.f_h_minus[c] = F_h_minus[c] / dt;
         }
         audit.samples.push_back(sample);

         if (elem1_on_plus)
         {
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < edofs1.Size(); i++)
               {
                  audit.rhs_before_minv(c * ndof_total + edofs1[i]) -=
                     w * shape1(i) * F_h_plus[c];
               }
               for (int i = 0; i < edofs2.Size(); i++)
               {
                  audit.rhs_before_minv(c * ndof_total + edofs2[i]) +=
                     w * shape2(i) * F_h_minus[c];
               }
            }
         }
         else
         {
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < edofs1.Size(); i++)
               {
                  audit.rhs_before_minv(c * ndof_total + edofs1[i]) +=
                     w * shape1(i) * F_h_minus[c];
               }
               for (int i = 0; i < edofs2.Size(); i++)
               {
                  audit.rhs_before_minv(c * ndof_total + edofs2[i]) -=
                     w * shape2(i) * F_h_plus[c];
               }
            }
         }
      }
   }

   audit.q_after_minv = audit.rhs_before_minv;
   ApplyMassInverseManually(wave, audit.q_after_minv);
   return audit;
}

// Round-13A: fault-only deposit audit with three mean-override modes.
// Structurally parallel to RunADERTraceAudit, but instead of calling
// `ff.EvaluateADER` per QP, this helper:
//   1. computes time-integrated I via ComputeADERTimeIntegrated(Q_step0),
//   2. gathers per-QP `EvalStageState` via `ff.ComputeStageState`,
//   3. (if Face or Global mode) overrides the selected stage field
//      with its face mean (Face) or global mean (Global),
//   4. runs the appropriate completion helper,
//   5. builds imposed Q± via `ff.BuildImposedState`,
//   6. rotates back to the global frame, computes F_h±, and deposits
//      into `audit.rhs_before_minv` with the same sign convention as
//      the production `ComputeADERFaceFluxRHS` fault branch.
// Does NOT call WriteBackState — successive calls with different
// modes must leave `dof_data` pristine for the caller's subsequent
// inspection.  Does NOT call wave.AdvanceADER.  Processes fault faces
// only.
FaultFaceMeanAudit RunFaultFaceMeanAudit(WaveOperator<Mesh> &wave,
                                         Mesh &mesh,
                                         FaultFaceFlux &ff,
                                         std::vector<DOFData> &dof_data,
                                         const Vector &Q_step0,
                                         real_t dt,
                                         FaultMeanMode mean_mode,
                                         FaultMeanSource mean_source)
{
   FaultFaceMeanAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();
   audit.rhs_before_minv.SetSize(size);
   audit.rhs_before_minv = 0.0;

   // Time-integrated state used by the ADER fault Riemann solve.
   Vector I_data(size);
   wave.ComputeADERTimeIntegrated(Q_step0, dt, kAderOrder, I_data);

   const FaultBasis *fault_basis = wave.GetFaultBasis();
   MFEM_VERIFY(fault_basis != nullptr, "fault basis must be populated");
   const auto &flux = wave.GetFlux();
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();

   // Per-QP scratch.  We keep everything we need for passes 2 and 3
   // in parallel std::vectors (per fault-face and per-QP flattened).
   struct PerQP
   {
      int face_list_idx = -1;   // index into int_faces
      int qp_idx        = -1;
      int face_mesh_idx = -1;   // mesh face index
      int dof_idx       = -1;
      int edof1_offset  = 0;
      int edof2_offset  = 0;
      Array<int> edofs1;
      Array<int> edofs2;
      Vector shape1;
      Vector shape2;
      real_t w = 0.0;
      std::array<real_t, 3> can_n{};
      DenseMatrix T_can;
      bool elem1_on_plus = false;
      std::array<real_t, NUM_STATE> I_plus_local{};
      std::array<real_t, NUM_STATE> I_minus_local{};
      EvalStageState state;
   };
   std::vector<PerQP> qps;
   audit.faces.clear();

   // Pass 1: gather per-QP data and compute the full stage chain.
   for (int fi = 0; fi < int_faces.Size(); fi++)
   {
      const int face_idx = int_faces[fi];
      auto *ftr = mesh.GetInteriorFaceTransformations(face_idx);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                               2 * kOrder);

      Array<int> edofs1, edofs2;
      fes.GetElementDofs(ftr->Elem1No, edofs1);
      fes.GetElementDofs(ftr->Elem2No, edofs2);
      const FiniteElement *fe1 = fes.GetFE(ftr->Elem1No);
      const FiniteElement *fe2 = fes.GetFE(ftr->Elem2No);

      const int fb_idx = wave.LookupInteriorFaultBasisIndex(face_idx);
      MFEM_VERIFY(fb_idx >= 0,
                  "missing FaultBasis entry for interior fault face");
      const FaultBasisData &bd = fault_basis->GetBasis(fb_idx);
      MFEM_VERIFY(static_cast<int>(bd.qp_data.size()) == ir.GetNPoints(),
                  "FaultBasis qp count mismatch");

      // Face centroid from the mesh (physical space).
      Array<int> fv;
      mesh.GetFaceVertices(face_idx, fv);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < fv.Size(); v++)
      {
         cx += mesh.GetVertex(fv[v])[0];
         cy += mesh.GetVertex(fv[v])[1];
         cz += mesh.GetVertex(fv[v])[2];
      }
      if (fv.Size() > 0)
      {
         cx /= fv.Size(); cy /= fv.Size(); cz /= fv.Size();
      }

      FaultFaceMeanSample face_sample;
      face_sample.face_idx = face_idx;
      face_sample.cx = cx; face_sample.cy = cy; face_sample.cz = cz;

      // Per-face accumulators for means.
      real_t sum_sn_trial = 0.0, sum_t1_trial = 0.0, sum_t2_trial = 0.0;
      real_t sum_sn_corr  = 0.0, sum_t1_corr  = 0.0, sum_t2_corr  = 0.0;
      const int face_nqp = ir.GetNPoints();

      // Index of the first per-QP slot for this face in `qps`.
      const int qps_face_base = static_cast<int>(qps.size());

      for (int q = 0; q < face_nqp; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector n_raw(3);
         CalcOrtho(ftr->Face->Jacobian(), n_raw);
         const real_t w = ip.weight * n_raw.Norml2();

         IntegrationPoint ip1, ip2;
         ftr->Loc1.Transform(ip, ip1);
         ftr->Loc2.Transform(ip, ip2);
         Vector shape1(fe1->GetDof()), shape2(fe2->GetDof());
         fe1->CalcShape(ip1, shape1);
         fe2->CalcShape(ip2, shape2);

         real_t I_self[NUM_STATE], I_nbr[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self[c] = 0.0; I_nbr[c] = 0.0;
            for (int i = 0; i < edofs1.Size(); i++)
            {
               I_self[c] += shape1(i) * I_data(c * ndof_total + edofs1[i]);
            }
            for (int i = 0; i < edofs2.Size(); i++)
            {
               I_nbr[c] += shape2(i) * I_data(c * ndof_total + edofs2[i]);
            }
         }

         const FaultBasisQPData &qpd = bd.qp_data[q];
         real_t can_n[3], can_t1[3], can_t2[3];
         for (int d = 0; d < 3; d++)
         {
            can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]   : qpd.normal[d];
            can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d] : qpd.tangent1[d];
            can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d] : qpd.tangent2[d];
         }

         DenseMatrix T_can(NUM_STATE), Tinv_can(NUM_STATE);
         GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);
         GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv_can);

         real_t I_self_can[NUM_STATE], I_nbr_can[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self_can[c] = 0.0; I_nbr_can[c] = 0.0;
            for (int k = 0; k < NUM_STATE; k++)
            {
               I_self_can[c] += Tinv_can(c, k) * I_self[k];
               I_nbr_can[c]  += Tinv_can(c, k) * I_nbr[k];
            }
         }

         const bool elem1_on_plus = !qpd.sign_flipped;

         PerQP p;
         p.face_list_idx  = fi;
         p.qp_idx         = q;
         p.face_mesh_idx  = face_idx;
         p.dof_idx        = fi * face_nqp + q;
         p.edofs1         = edofs1;
         p.edofs2         = edofs2;
         p.shape1         = shape1;
         p.shape2         = shape2;
         p.w              = w;
         for (int d = 0; d < 3; d++) { p.can_n[d] = can_n[d]; }
         p.T_can          = T_can;
         p.elem1_on_plus  = elem1_on_plus;

         for (int c = 0; c < NUM_STATE; c++)
         {
            const real_t i_plus  = elem1_on_plus ? I_self_can[c] : I_nbr_can[c];
            const real_t i_minus = elem1_on_plus ? I_nbr_can[c]  : I_self_can[c];
            p.I_plus_local[c]  = i_plus;
            p.I_minus_local[c] = i_minus;
         }

         // Compute stage state on the time-AVERAGED Q (I / dt).
         real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
         const real_t inv_dt = 1.0 / dt;
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_avg_plus[c]  = p.I_plus_local[c]  * inv_dt;
            Q_avg_minus[c] = p.I_minus_local[c] * inv_dt;
         }
         const DOFData &fdata_const = dof_data[p.dof_idx];
         ff.ComputeStageState(fdata_const, Q_avg_plus, Q_avg_minus, p.state);

         sum_sn_trial += p.state.sigma_n_trial;
         sum_t1_trial += p.state.tau1_trial;
         sum_t2_trial += p.state.tau2_trial;
         sum_sn_corr  += p.state.sigma_n_corr;
         sum_t1_corr  += p.state.tau1_corr;
         sum_t2_corr  += p.state.tau2_corr;

         qps.push_back(std::move(p));
      }

      const real_t inv_fnqp = (face_nqp > 0)
                              ? 1.0 / static_cast<real_t>(face_nqp)
                              : 0.0;
      face_sample.sigma_n_trial_mean = sum_sn_trial * inv_fnqp;
      face_sample.tau1_trial_mean    = sum_t1_trial * inv_fnqp;
      face_sample.tau2_trial_mean    = sum_t2_trial * inv_fnqp;
      face_sample.sigma_n_corr_mean  = sum_sn_corr  * inv_fnqp;
      face_sample.tau1_corr_mean     = sum_t1_corr  * inv_fnqp;
      face_sample.tau2_corr_mean     = sum_t2_corr  * inv_fnqp;
      // face_sample.elem1_on_plus and F_h fields filled in pass 3.
      face_sample.elem1_on_plus = qps[qps_face_base].elem1_on_plus;
      audit.faces.push_back(face_sample);
   }

   // Pass 2: compute mean-overrides per the selected mode.
   if (mean_mode != FaultMeanMode::None && !qps.empty())
   {
      if (mean_mode == FaultMeanMode::Face)
      {
         // Per-face means already computed and stored in audit.faces.
         // Overwrite each QP's selected stage field with its face mean.
         for (auto &p : qps)
         {
            const FaultFaceMeanSample &fs = audit.faces[p.face_list_idx];
            if (mean_source == FaultMeanSource::Trial)
            {
               p.state.sigma_n_trial = fs.sigma_n_trial_mean;
               p.state.tau1_trial    = fs.tau1_trial_mean;
               p.state.tau2_trial    = fs.tau2_trial_mean;
               ff.CompleteFromTrial(dof_data[p.dof_idx], p.state);
            }
            else
            {
               p.state.sigma_n_corr = fs.sigma_n_corr_mean;
               p.state.tau1_corr    = fs.tau1_corr_mean;
               p.state.tau2_corr    = fs.tau2_corr_mean;
               // Tcorr: no downstream recompute.
            }
         }
      }
      else   // Global
      {
         real_t gsum_sn_t = 0.0, gsum_t1_t = 0.0, gsum_t2_t = 0.0;
         real_t gsum_sn_c = 0.0, gsum_t1_c = 0.0, gsum_t2_c = 0.0;
         for (const auto &p : qps)
         {
            gsum_sn_t += p.state.sigma_n_trial;
            gsum_t1_t += p.state.tau1_trial;
            gsum_t2_t += p.state.tau2_trial;
            gsum_sn_c += p.state.sigma_n_corr;
            gsum_t1_c += p.state.tau1_corr;
            gsum_t2_c += p.state.tau2_corr;
         }
         const real_t inv_n = 1.0 / static_cast<real_t>(qps.size());
         const real_t g_sn_t = gsum_sn_t * inv_n;
         const real_t g_t1_t = gsum_t1_t * inv_n;
         const real_t g_t2_t = gsum_t2_t * inv_n;
         const real_t g_sn_c = gsum_sn_c * inv_n;
         const real_t g_t1_c = gsum_t1_c * inv_n;
         const real_t g_t2_c = gsum_t2_c * inv_n;
         for (auto &p : qps)
         {
            if (mean_source == FaultMeanSource::Trial)
            {
               p.state.sigma_n_trial = g_sn_t;
               p.state.tau1_trial    = g_t1_t;
               p.state.tau2_trial    = g_t2_t;
               ff.CompleteFromTrial(dof_data[p.dof_idx], p.state);
            }
            else
            {
               p.state.sigma_n_corr = g_sn_c;
               p.state.tau1_corr    = g_t1_c;
               p.state.tau2_corr    = g_t2_c;
            }
         }
      }
   }

   // Pass 3: build imposed states, compute F_h, deposit into rhs.
   // Accumulate per-face F_h means for Gate 4b reporting.
   std::vector<real_t> fh_plus_sxy_sum(audit.faces.size(), 0.0);
   std::vector<real_t> fh_plus_sxz_sum(audit.faces.size(), 0.0);
   std::vector<real_t> fh_minus_sxy_sum(audit.faces.size(), 0.0);
   std::vector<real_t> fh_minus_sxz_sum(audit.faces.size(), 0.0);
   std::vector<int>    fh_face_count(audit.faces.size(), 0);

   const real_t inv_dt = 1.0 / dt;

   for (const auto &p : qps)
   {
      real_t Q_avg_plus[NUM_STATE], Q_avg_minus[NUM_STATE];
      for (int c = 0; c < NUM_STATE; c++)
      {
         Q_avg_plus[c]  = p.I_plus_local[c]  * inv_dt;
         Q_avg_minus[c] = p.I_minus_local[c] * inv_dt;
      }

      real_t Q_imp_plus[NUM_STATE], Q_imp_minus[NUM_STATE];
      ff.BuildImposedState(dof_data[p.dof_idx], p.state,
                           Q_avg_plus, Q_avg_minus,
                           Q_imp_plus, Q_imp_minus);

      real_t I_imp_plus[NUM_STATE], I_imp_minus[NUM_STATE];
      for (int c = 0; c < NUM_STATE; c++)
      {
         I_imp_plus[c]  = Q_imp_plus[c]  * dt;
         I_imp_minus[c] = Q_imp_minus[c] * dt;
      }

      // Rotate imposed states back to the global frame.
      real_t I_imp_plus_g[NUM_STATE], I_imp_minus_g[NUM_STATE];
      for (int c = 0; c < NUM_STATE; c++)
      {
         I_imp_plus_g[c]  = 0.0;
         I_imp_minus_g[c] = 0.0;
         for (int k = 0; k < NUM_STATE; k++)
         {
            I_imp_plus_g[c]  += p.T_can(c, k) * I_imp_plus[k];
            I_imp_minus_g[c] += p.T_can(c, k) * I_imp_minus[k];
         }
      }

      real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
      flux.Interior(p.can_n.data(), I_imp_plus_g,  I_imp_plus_g,  F_h_plus);
      flux.Interior(p.can_n.data(), I_imp_minus_g, I_imp_minus_g, F_h_minus);

      // Per-face F_h means (for Gate 4b reporting).
      fh_plus_sxy_sum [p.face_list_idx] += F_h_plus [SXY];
      fh_plus_sxz_sum [p.face_list_idx] += F_h_plus [SXZ];
      fh_minus_sxy_sum[p.face_list_idx] += F_h_minus[SXY];
      fh_minus_sxz_sum[p.face_list_idx] += F_h_minus[SXZ];
      fh_face_count   [p.face_list_idx] += 1;

      if (p.elem1_on_plus)
      {
         for (int c = 0; c < NUM_STATE; c++)
         {
            for (int i = 0; i < p.edofs1.Size(); i++)
            {
               audit.rhs_before_minv(c * ndof_total + p.edofs1[i]) -=
                  p.w * p.shape1(i) * F_h_plus[c];
            }
            for (int i = 0; i < p.edofs2.Size(); i++)
            {
               audit.rhs_before_minv(c * ndof_total + p.edofs2[i]) +=
                  p.w * p.shape2(i) * F_h_minus[c];
            }
         }
      }
      else
      {
         for (int c = 0; c < NUM_STATE; c++)
         {
            for (int i = 0; i < p.edofs1.Size(); i++)
            {
               audit.rhs_before_minv(c * ndof_total + p.edofs1[i]) +=
                  p.w * p.shape1(i) * F_h_minus[c];
            }
            for (int i = 0; i < p.edofs2.Size(); i++)
            {
               audit.rhs_before_minv(c * ndof_total + p.edofs2[i]) -=
                  p.w * p.shape2(i) * F_h_plus[c];
            }
         }
      }
   }

   // Store per-face F_h means.
   for (size_t fi = 0; fi < audit.faces.size(); fi++)
   {
      const int n = std::max(fh_face_count[fi], 1);
      const real_t inv = 1.0 / static_cast<real_t>(n);
      audit.faces[fi].fh_plus_sxy_mean  = fh_plus_sxy_sum [fi] * inv;
      audit.faces[fi].fh_plus_sxz_mean  = fh_plus_sxz_sum [fi] * inv;
      audit.faces[fi].fh_minus_sxy_mean = fh_minus_sxy_sum[fi] * inv;
      audit.faces[fi].fh_minus_sxz_mean = fh_minus_sxz_sum[fi] * inv;
   }

   audit.rhs_after_minv = audit.rhs_before_minv;
   ApplyMassInverseManually(wave, audit.rhs_after_minv);
   return audit;
}

// Round-13B: second-step non-fault face seed audit.  Given Q_step0,
// compute I = ComputeADERTimeIntegrated(Q_step0, dt), iterate every
// non-fault face (interior + boundary), compute both the raw and
// symmetrized F_h per face, optionally orbit-average, then deposit
// into rhs with the same sign convention as production.  Returns the
// post-M^-1 rhs and the per-face samples.  Does not touch fault faces
// or fault DOFData.
NonFaultSeedAudit RunSecondStepNonFaultSeedAudit(WaveOperator<Mesh> &wave,
                                                 Mesh &mesh,
                                                 const BoundaryConfig &bc,
                                                 const Vector &Q_step0,
                                                 real_t dt,
                                                 NonFaultSeedMode mode)
{
   NonFaultSeedAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();
   audit.rhs_before_minv.SetSize(size);
   audit.rhs_before_minv = 0.0;

   Vector I_data(size);
   wave.ComputeADERTimeIntegrated(Q_step0, dt, kAderOrder, I_data);

   real_t bulk_bg_scaled[NUM_STATE] = {0.0};
   if (const real_t *bulk_bg = wave.GetAbsorbingBackground())
   {
      for (int c = 0; c < NUM_STATE; c++) { bulk_bg_scaled[c] = dt * bulk_bg[c]; }
   }

   const auto &flux = wave.GetFlux();

   // Pass 1: gather per-face samples.  Mirrors the structure of
   // RunADERNonFaultFaceAudit but computes BOTH fh_raw and fh_sym
   // unconditionally so Pass 2 can pick between them.
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;

      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < fv.Size(); v++)
      {
         cx += mesh.GetVertex(fv[v])[0];
         cy += mesh.GetVertex(fv[v])[1];
         cz += mesh.GetVertex(fv[v])[2];
      }
      if (fv.Size() > 0)
      {
         cx /= fv.Size(); cy /= fv.Size(); cz /= fv.Size();
      }

      const bool is_boundary = (e2 < 0);
      // Skip fault faces: interior face with centroid at y = L/2.
      const bool is_fault = (e2 >= 0) && (std::abs(cy - 0.5 * kL) < 1e-8);
      if (is_fault) { continue; }

      const FiniteElement *fe1 = fes.GetFE(e1);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = e1 * wave.GetNDof();
      int dof_offset2 = -1;
      const FiniteElement *fe2 = nullptr;
      if (!is_boundary)
      {
         fe2 = fes.GetFE(e2);
         dof_offset2 = e2 * wave.GetNDof();
      }
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                               2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }
         const real_t w = ip.weight * nor_len;
         real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         Vector shape2;
         if (!is_boundary)
         {
            shape2.SetSize(ndof);
            IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2);
            fe2->CalcShape(ip2, shape2);
         }

         real_t I_self[NUM_STATE], I_nbr[NUM_STATE] = {0.0};
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               I_self[c] += shape1(i) * I_data(c * ndof_total + dof_offset1 + i);
            }
         }
         if (!is_boundary)
         {
            for (int c = 0; c < NUM_STATE; c++)
            {
               I_nbr[c] = 0.0;
               for (int i = 0; i < ndof; i++)
               {
                  I_nbr[c] += shape2(i) * I_data(c * ndof_total + dof_offset2 + i);
               }
            }
         }

         NonFaultFaceSample s;
         s.face_idx    = f;
         s.qp_idx      = q;
         s.e1          = e1;
         s.e2          = e2;
         s.dof_offset1 = dof_offset1;
         s.dof_offset2 = dof_offset2;
         s.is_boundary = is_boundary;
         s.side        = is_boundary ? ClassifyBoundarySide(cx, cy, cz)
                                     : BoundarySide::None;
         s.cx = cx; s.cy = cy; s.cz = cz;
         s.w  = w;
         for (int d = 0; d < 3; d++) { s.nor[d] = nor[d]; }
         s.shape1 = shape1;
         s.shape2 = shape2;

         if (is_boundary)
         {
            real_t F_raw[NUM_STATE], F_sym[NUM_STATE];
            ComputeBoundaryFhRaw(wave, bc, nor, I_self, bulk_bg_scaled, F_raw);
            ComputeBoundaryFhSym(wave, bc, nor, I_self, bulk_bg_scaled, F_sym);
            for (int c = 0; c < NUM_STATE; c++)
            {
               s.fh_raw[c] = F_raw[c];
               s.fh_sym[c] = F_sym[c];
            }
         }
         else
         {
            real_t F_raw[NUM_STATE], F_sym[NUM_STATE];
            ComputeInteriorFhRaw(flux, nor, I_self, I_nbr, F_raw);
            ComputeInteriorFhSym(flux, nor, I_self, I_nbr, F_sym);
            for (int c = 0; c < NUM_STATE; c++)
            {
               s.fh_raw[c] = F_raw[c];
               s.fh_sym[c] = F_sym[c];
            }
         }

         audit.faces.push_back(std::move(s));
      }
   }

   // Pass 2: produce the per-sample F_h actually deposited, according
   // to the mode.
   const int nfaces = static_cast<int>(audit.faces.size());
   std::vector<std::array<real_t, NUM_STATE>> fh_use(nfaces);
   auto set_from = [&](int i, bool use_sym)
   {
      const auto &s = audit.faces[i];
      for (int c = 0; c < NUM_STATE; c++)
      {
         fh_use[i][c] = use_sym ? s.fh_sym[c] : s.fh_raw[c];
      }
   };
   // Default: raw.
   for (int i = 0; i < nfaces; i++) { set_from(i, /*use_sym=*/false); }

   // Symmetrization overrides.
   if (mode == NonFaultSeedMode::InteriorSym ||
       mode == NonFaultSeedMode::BothSym)
   {
      for (int i = 0; i < nfaces; i++)
      {
         if (!audit.faces[i].is_boundary) { set_from(i, true); }
      }
   }
   if (mode == NonFaultSeedMode::BoundarySym ||
       mode == NonFaultSeedMode::BothSym)
   {
      for (int i = 0; i < nfaces; i++)
      {
         if (audit.faces[i].is_boundary) { set_from(i, true); }
      }
   }

   // Face-mean orbit averaging overrides (applied after raw baseline;
   // sym + mean are not combined per plan §4).
   auto apply_orbit_mean = [&](bool boundary_only, bool interior_only)
   {
      std::map<NonFaultFaceOrbitKey,
               std::pair<std::array<real_t, NUM_STATE>, int>> bucket_sum;
      for (int i = 0; i < nfaces; i++)
      {
         const auto &s = audit.faces[i];
         if (boundary_only && !s.is_boundary) { continue; }
         if (interior_only &&  s.is_boundary) { continue; }
         const auto key = MakeNonFaultFaceOrbitKey(s);
         auto &entry = bucket_sum[key];
         for (int c = 0; c < NUM_STATE; c++) { entry.first[c] += fh_use[i][c]; }
         entry.second += 1;
      }
      std::map<NonFaultFaceOrbitKey, std::array<real_t, NUM_STATE>> bucket_mean;
      for (const auto &kv : bucket_sum)
      {
         std::array<real_t, NUM_STATE> m{};
         const real_t inv = 1.0 / static_cast<real_t>(kv.second.second);
         for (int c = 0; c < NUM_STATE; c++) { m[c] = kv.second.first[c] * inv; }
         bucket_mean[kv.first] = m;
      }
      for (int i = 0; i < nfaces; i++)
      {
         const auto &s = audit.faces[i];
         if (boundary_only && !s.is_boundary) { continue; }
         if (interior_only &&  s.is_boundary) { continue; }
         const auto key = MakeNonFaultFaceOrbitKey(s);
         auto it = bucket_mean.find(key);
         if (it == bucket_mean.end()) { continue; }
         for (int c = 0; c < NUM_STATE; c++) { fh_use[i][c] = it->second[c]; }
      }
   };

   if (mode == NonFaultSeedMode::InteriorFaceMean ||
       mode == NonFaultSeedMode::BothFaceMean)
   {
      apply_orbit_mean(/*boundary_only=*/false, /*interior_only=*/true);
   }
   if (mode == NonFaultSeedMode::BoundaryFaceMean ||
       mode == NonFaultSeedMode::BothFaceMean)
   {
      apply_orbit_mean(/*boundary_only=*/true, /*interior_only=*/false);
   }

   // Round-14A: interior-sym background + conditional x-boundary
   // treatment.  Deviation from plan §1.3's SelectNonFaultFh helper:
   // to preserve parallelism with the existing set_from / fh_use
   // infrastructure, the same effect is implemented inline here.
   // Each mode turns on interior sym across all interior non-fault
   // faces, then toggles boundary sym only on the named x-side(s).
   // Non-x boundaries (y, z, interior-class=None on boundary) remain
   // raw under these modes.
   if (mode == NonFaultSeedMode::InteriorSym_XBoundaryRaw ||
       mode == NonFaultSeedMode::InteriorSym_XBoundarySym ||
       mode == NonFaultSeedMode::InteriorSym_XMinSym ||
       mode == NonFaultSeedMode::InteriorSym_XMaxSym)
   {
      // Interior non-fault faces: sym.
      for (int i = 0; i < nfaces; i++)
      {
         if (!audit.faces[i].is_boundary) { set_from(i, true); }
      }
      // Boundary faces: sym only on the selected x-side(s).
      for (int i = 0; i < nfaces; i++)
      {
         const auto &s = audit.faces[i];
         if (!s.is_boundary) { continue; }
         bool apply_sym = false;
         if (mode == NonFaultSeedMode::InteriorSym_XBoundarySym)
         {
            apply_sym = (s.side == BoundarySide::XMin ||
                         s.side == BoundarySide::XMax);
         }
         else if (mode == NonFaultSeedMode::InteriorSym_XMinSym)
         {
            apply_sym = (s.side == BoundarySide::XMin);
         }
         else if (mode == NonFaultSeedMode::InteriorSym_XMaxSym)
         {
            apply_sym = (s.side == BoundarySide::XMax);
         }
         // InteriorSym_XBoundaryRaw: no boundary sym at all.
         if (apply_sym) { set_from(i, true); }
      }
   }

   // Pass 3: deposit using the production sign convention.  Boundary:
   // subtract from e1 only.  Interior: subtract from e1, add to e2.
   for (int i = 0; i < nfaces; i++)
   {
      const auto &s = audit.faces[i];
      const auto &F  = fh_use[i];
      const int ndof = s.shape1.Size();
      if (s.is_boundary)
      {
         for (int c = 0; c < NUM_STATE; c++)
         {
            for (int k = 0; k < ndof; k++)
            {
               audit.rhs_before_minv(c * ndof_total + s.dof_offset1 + k) -=
                  s.w * s.shape1(k) * F[c];
            }
         }
      }
      else
      {
         for (int c = 0; c < NUM_STATE; c++)
         {
            for (int k = 0; k < ndof; k++)
            {
               audit.rhs_before_minv(c * ndof_total + s.dof_offset1 + k) -=
                  s.w * s.shape1(k) * F[c];
               audit.rhs_before_minv(c * ndof_total + s.dof_offset2 + k) +=
                  s.w * s.shape2(k) * F[c];
            }
         }
      }
   }

   audit.rhs_after_minv = audit.rhs_before_minv;
   ApplyMassInverseManually(wave, audit.rhs_after_minv);
   return audit;
}

// Round-15 forced-kernel override: bypasses DispatchBoundaryBC and
// forces a specific boundary kernel on every x-side face.  Gate 4i
// is run three extra times with these overrides if the first gate
// shows I_self clean but F_h dirty on x-sides.
enum class BoundaryKernelOverride
{
   None,
   ForceAbsorbing,
   ForceFreeSurface,
   ForceFreeSurfaceGodunov
};

const char *BoundaryKernelOverrideLabel(BoundaryKernelOverride o)
{
   switch (o)
   {
      case BoundaryKernelOverride::None:                   return "dispatch";
      case BoundaryKernelOverride::ForceAbsorbing:         return "force_abs";
      case BoundaryKernelOverride::ForceFreeSurface:       return "force_fs_gamma";
      case BoundaryKernelOverride::ForceFreeSurfaceGodunov: return "force_fs_god";
   }
   return "?";
}

// Round-15 trace-only boundary audit.  Reuses the gather path from
// RunSecondStepNonFaultSeedAudit but filters to boundary faces on
// x=0 / x=L (diagnostic target) and z=0 / z=L (quiet control).
// Does NOT deposit and does NOT touch any production state.
XSideBoundaryTraceAudit RunXSideBoundaryTraceAudit(
   WaveOperator<Mesh> &wave,
   Mesh &mesh,
   const BoundaryConfig &bc,
   const Vector &Q_step0,
   real_t dt,
   BoundaryKernelOverride override_mode = BoundaryKernelOverride::None)
{
   XSideBoundaryTraceAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   Vector I_data(wave.Height());
   wave.ComputeADERTimeIntegrated(Q_step0, dt, kAderOrder, I_data);

   real_t bulk_bg_scaled_vec[NUM_STATE] = {0.0};
   if (const real_t *bulk_bg = wave.GetAbsorbingBackground())
   {
      for (int c = 0; c < NUM_STATE; c++) { bulk_bg_scaled_vec[c] = dt * bulk_bg[c]; }
   }

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      if (e2 >= 0) { continue; }   // boundary faces only

      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < fv.Size(); v++)
      {
         cx += mesh.GetVertex(fv[v])[0];
         cy += mesh.GetVertex(fv[v])[1];
         cz += mesh.GetVertex(fv[v])[2];
      }
      if (fv.Size() > 0)
      {
         cx /= fv.Size(); cy /= fv.Size(); cz /= fv.Size();
      }

      const BoundarySide side = ClassifyBoundarySide(cx, cy, cz);
      const bool want_x = (side == BoundarySide::XMin || side == BoundarySide::XMax);
      const bool want_z = (side == BoundarySide::ZMin || side == BoundarySide::ZMax);
      if (!want_x && !want_z) { continue; }

      // The audit file historically uses boundary_attr = 1 (all
      // exterior faces in this fixture; natural_attrs = {1}).
      // `WaveOperator::ClassifyBoundaryFace` is private, so the
      // classification is replicated locally to mirror the
      // DispatchBoundaryBC logic used throughout this file:
      //   absorbing_attrs → Absorbing
      //   natural_attrs   → FreeSurface
      //   else            → Absorbing (fallback)
      const int bdr_attr = 1;
      AuditFaceBC bc_type;
      if (bc.absorbing_attrs.count(bdr_attr))
      {
         bc_type = AuditFaceBC::Absorbing;
      }
      else if (bc.natural_attrs.count(bdr_attr))
      {
         bc_type = AuditFaceBC::FreeSurface;
      }
      else
      {
         bc_type = AuditFaceBC::Absorbing;
      }

      const FiniteElement *fe1 = fes.GetFE(e1);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = e1 * wave.GetNDof();
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                               2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }
         const real_t w = ip.weight * nor_len;
         real_t nor_arr[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         real_t I_self[NUM_STATE] = {0.0};
         for (int c = 0; c < NUM_STATE; c++)
         {
            for (int i = 0; i < ndof; i++)
            {
               I_self[c] += shape1(i) * I_data(c * ndof_total + dof_offset1 + i);
            }
         }

         real_t F_h[NUM_STATE] = {0.0};
         if (override_mode == BoundaryKernelOverride::None)
         {
            ComputeBoundaryFhRaw(wave, bc, nor_arr, I_self,
                                 bulk_bg_scaled_vec, F_h);
         }
         else if (override_mode == BoundaryKernelOverride::ForceAbsorbing)
         {
            wave.GetFlux().AbsorbingTotal(nor_arr, I_self,
                                          bulk_bg_scaled_vec, F_h);
         }
         else if (override_mode == BoundaryKernelOverride::ForceFreeSurface)
         {
            wave.GetFlux().FreeSurfaceTotal(nor_arr, I_self,
                                            bulk_bg_scaled_vec, F_h);
         }
         else   // ForceFreeSurfaceGodunov
         {
            wave.GetFlux().FreeSurfaceGodunovTotal(nor_arr, I_self,
                                                    bulk_bg_scaled_vec, F_h);
         }

         XSideBoundaryTraceSample s;
         s.face_idx   = f;
         s.elem       = e1;
         s.dof_offset = dof_offset1;
         s.bdr_attr   = bdr_attr;
         s.q          = q;
         s.side       = side;
         s.bc_type    = bc_type;
         s.cx = cx; s.cy = cy; s.cz = cz;
         s.w  = w;
         for (int d = 0; d < 3; d++) { s.nor[d] = nor_arr[d]; }
         for (int c = 0; c < NUM_STATE; c++)
         {
            s.I_self[c] = I_self[c];
            s.bulk_bg_scaled[c] = bulk_bg_scaled_vec[c];
            s.F_h[c] = F_h[c];
         }
         s.shape1.assign(shape1.GetData(), shape1.GetData() + shape1.Size());

         if (want_x) { audit.x_faces.push_back(std::move(s)); }
         else        { audit.z_faces.push_back(std::move(s)); }
      }
   }
   return audit;
}

// Round-16 Gate 4k source: sample I_self[SXY] two ways on the same
// I_data and same QP — via production contiguous offset style AND
// via explicit GetElementDofs.  Trace-only; no flux, no deposit.
BoundaryTraceSamplerAudit RunBoundaryTraceSamplerComparisonAudit(
   WaveOperator<Mesh> &wave,
   Mesh &mesh,
   const Vector &Q_step0,
   real_t dt)
{
   BoundaryTraceSamplerAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   Vector I_data(wave.Height());
   wave.ComputeADERTimeIntegrated(Q_step0, dt, kAderOrder, I_data);

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }
      if (ftr->Elem2No >= 0) { continue; }    // boundary only

      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < fv.Size(); v++)
      {
         cx += mesh.GetVertex(fv[v])[0];
         cy += mesh.GetVertex(fv[v])[1];
         cz += mesh.GetVertex(fv[v])[2];
      }
      if (fv.Size() > 0)
      {
         cx /= fv.Size(); cy /= fv.Size(); cz /= fv.Size();
      }

      const BoundarySide side = ClassifyBoundarySide(cx, cy, cz);
      const bool want_x = (side == BoundarySide::XMin || side == BoundarySide::XMax);
      const bool want_z = (side == BoundarySide::ZMin || side == BoundarySide::ZMax);
      if (!want_x && !want_z) { continue; }

      const int e1 = ftr->Elem1No;
      const FiniteElement *fe1 = fes.GetFE(e1);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = e1 * wave.GetNDof();

      Array<int> edofs1;
      fes.GetElementDofs(e1, edofs1);

      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                               2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         // Two samplings of the SXY block from the same I_data.
         real_t sxy_offset = 0.0;
         real_t sxy_edofs  = 0.0;
         for (int i = 0; i < ndof; i++)
         {
            sxy_offset += shape1(i) * I_data(SXY * ndof_total + dof_offset1 + i);
            sxy_edofs  += shape1(i) * I_data(SXY * ndof_total + edofs1[i]);
         }

         BoundaryTraceSamplerSample s;
         s.face_idx   = f;
         s.elem       = e1;
         s.q          = q;
         s.dof_offset = dof_offset1;
         s.side       = side;
         s.cx = cx; s.cy = cy; s.cz = cz;
         for (int d = 0; d < 3; d++) { s.nor[d] = nor_vec(d); }
         s.edofs1.assign(edofs1.begin(), edofs1.end());
         s.shape1.assign(shape1.GetData(), shape1.GetData() + shape1.Size());
         s.i_self_offset_sxy = sxy_offset;
         s.i_self_edofs_sxy  = sxy_edofs;

         if (want_x) { audit.x_faces.push_back(std::move(s)); }
         else        { audit.z_faces.push_back(std::move(s)); }
      }
   }
   return audit;
}

// Round-16 Gate 4l source: element-level I_data[SXY] coefficients
// per element DOF, with element-orbit grouping by canonicalized
// centroid.  Answers: is I_data[SXY] already orbit-asymmetric at
// the element level, before any face interpolation?
ElementIDataAudit RunElementIDataAudit(WaveOperator<Mesh> &wave,
                                        Mesh &mesh,
                                        const Vector &Q_step0,
                                        real_t dt)
{
   ElementIDataAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   Vector I_data(wave.Height());
   wave.ComputeADERTimeIntegrated(Q_step0, dt, kAderOrder, I_data);

   for (int e = 0; e < wave.NumElements(); e++)
   {
      ElementIDataSample s;
      s.elem = e;

      Array<int> ev;
      mesh.GetElementVertices(e, ev);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < ev.Size(); v++)
      {
         cx += mesh.GetVertex(ev[v])[0];
         cy += mesh.GetVertex(ev[v])[1];
         cz += mesh.GetVertex(ev[v])[2];
      }
      if (ev.Size() > 0)
      {
         cx /= ev.Size(); cy /= ev.Size(); cz /= ev.Size();
      }
      s.cx = cx; s.cy = cy; s.cz = cz;

      Array<int> edofs;
      fes.GetElementDofs(e, edofs);
      s.edofs.assign(edofs.begin(), edofs.end());
      s.i_sxy_dofs.resize(edofs.Size());
      real_t sum = 0.0;
      for (int i = 0; i < edofs.Size(); i++)
      {
         const real_t v = I_data(SXY * ndof_total + edofs[i]);
         s.i_sxy_dofs[i] = v;
         sum += v;
      }
      if (edofs.Size() > 0) { s.i_sxy_mean = sum / edofs.Size(); }

      audit.elems.push_back(std::move(s));
   }
   return audit;
}

// Round-17 element SXY per-DOF audit: for every element, build both
// raw-order and canonical (sorted-by-reference-position) views of
// the local I_data[SXY] DOFs.  Answers: are orbit-paired elements
// equivalent up to a DOF permutation?
ElementSxyDofAudit RunElementSxyDofAudit(WaveOperator<Mesh> &wave,
                                          Mesh &mesh,
                                          const Vector &Q_step0,
                                          real_t dt)
{
   ElementSxyDofAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   Vector I_data(wave.Height());
   wave.ComputeADERTimeIntegrated(Q_step0, dt, kAderOrder, I_data);

   for (int e = 0; e < wave.NumElements(); e++)
   {
      ElementSxyDofSample s;
      s.elem = e;

      Array<int> ev;
      mesh.GetElementVertices(e, ev);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < ev.Size(); v++)
      {
         cx += mesh.GetVertex(ev[v])[0];
         cy += mesh.GetVertex(ev[v])[1];
         cz += mesh.GetVertex(ev[v])[2];
      }
      if (ev.Size() > 0)
      {
         cx /= ev.Size(); cy /= ev.Size(); cz /= ev.Size();
      }
      s.cx = cx; s.cy = cy; s.cz = cz;

      Array<int> edofs;
      fes.GetElementDofs(e, edofs);
      const FiniteElement *fe = fes.GetFE(e);

      s.raw_dofs       = BuildRawLocalDofs(*fe, edofs, I_data, ndof_total);
      s.canonical_dofs = BuildCanonicalLocalDofs(*fe, edofs, I_data, ndof_total);
      real_t sum = 0.0;
      for (const auto &d : s.raw_dofs) { sum += d.value; }
      if (!s.raw_dofs.empty()) { s.raw_mean = sum / s.raw_dofs.size(); }

      audit.elems.push_back(std::move(s));
   }
   return audit;
}

// Round-17 face interpolation canonicalization audit: for every x/z
// boundary non-fault face QP, store both raw-order and canonical-
// order local DOF vectors AND the matching shape functions.  The
// canonical order applies the SAME permutation to shape so
// `Σ shape_canonical[i] * canonical_dofs[i].value` is a permutation-
// invariant reformulation of the raw contraction.  Diagnostic only.
FaceInterpolationAudit RunFaceInterpolationAudit(WaveOperator<Mesh> &wave,
                                                 Mesh &mesh,
                                                 const Vector &Q_step0,
                                                 real_t dt)
{
   FaceInterpolationAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   Vector I_data(wave.Height());
   wave.ComputeADERTimeIntegrated(Q_step0, dt, kAderOrder, I_data);

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }
      if (ftr->Elem2No >= 0) { continue; }   // boundary only

      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < fv.Size(); v++)
      {
         cx += mesh.GetVertex(fv[v])[0];
         cy += mesh.GetVertex(fv[v])[1];
         cz += mesh.GetVertex(fv[v])[2];
      }
      if (fv.Size() > 0)
      {
         cx /= fv.Size(); cy /= fv.Size(); cz /= fv.Size();
      }

      const BoundarySide side = ClassifyBoundarySide(cx, cy, cz);
      const bool want_x = (side == BoundarySide::XMin || side == BoundarySide::XMax);
      const bool want_z = (side == BoundarySide::ZMin || side == BoundarySide::ZMax);
      if (!want_x && !want_z) { continue; }

      const int e1 = ftr->Elem1No;
      const FiniteElement *fe1 = fes.GetFE(e1);
      const int ndof = fe1->GetDof();

      Array<int> edofs1;
      fes.GetElementDofs(e1, edofs1);

      // Precompute raw + canonical DOFs and the local permutation.
      std::vector<CanonicalLocalDof> raw_dofs =
         BuildRawLocalDofs(*fe1, edofs1, I_data, ndof_total);
      std::vector<CanonicalLocalDof> canonical_dofs =
         BuildCanonicalLocalDofs(*fe1, edofs1, I_data, ndof_total);
      // Permutation: canonical_dofs[k] came from raw_dofs[perm[k]].
      std::vector<int> perm(ndof);
      for (int k = 0; k < ndof; k++) { perm[k] = canonical_dofs[k].local_idx; }

      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                               2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         FaceInterpolationSample s;
         s.face_idx = f;
         s.elem     = e1;
         s.q        = q;
         s.side     = side;
         s.cx = cx; s.cy = cy; s.cz = cz;
         for (int d = 0; d < 3; d++) { s.nor[d] = nor_vec(d); }
         s.raw_dofs       = raw_dofs;
         s.canonical_dofs = canonical_dofs;
         s.shape_raw.assign(shape1.GetData(), shape1.GetData() + shape1.Size());
         s.shape_canonical.resize(ndof);
         for (int k = 0; k < ndof; k++)
         {
            s.shape_canonical[k] = shape1(perm[k]);
         }

         real_t interp_raw = 0.0, interp_can = 0.0;
         for (int i = 0; i < ndof; i++)
         {
            interp_raw += s.shape_raw[i] * raw_dofs[i].value;
            interp_can += s.shape_canonical[i] * canonical_dofs[i].value;
         }
         s.interp_raw       = interp_raw;
         s.interp_canonical = interp_can;

         if (want_x) { audit.x_faces.push_back(std::move(s)); }
         else        { audit.z_faces.push_back(std::move(s)); }
      }
   }
   return audit;
}

// Round-18A: face-local canonical interpolation audit.  For every
// x/z boundary non-fault face QP, build three orderings — raw,
// element-canonical (sorted by full 3D reference position), and
// face-canonical (sorted by reference position restricted to the
// face plane).  Apply the same permutation to both DOF values and
// shape functions in each case.  Trace-only; no flux, no deposit.
FaceLocalCanonicalAudit RunFaceLocalCanonicalInterpolationAudit(
   WaveOperator<Mesh> &wave,
   Mesh &mesh,
   const Vector &Q_step0,
   real_t dt)
{
   FaceLocalCanonicalAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   Vector I_data(wave.Height());
   wave.ComputeADERTimeIntegrated(Q_step0, dt, kAderOrder, I_data);

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }
      if (ftr->Elem2No >= 0) { continue; }   // boundary only

      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < fv.Size(); v++)
      {
         cx += mesh.GetVertex(fv[v])[0];
         cy += mesh.GetVertex(fv[v])[1];
         cz += mesh.GetVertex(fv[v])[2];
      }
      if (fv.Size() > 0)
      {
         cx /= fv.Size(); cy /= fv.Size(); cz /= fv.Size();
      }

      const BoundarySide side = ClassifyBoundarySide(cx, cy, cz);
      const bool want_x = (side == BoundarySide::XMin || side == BoundarySide::XMax);
      const bool want_z = (side == BoundarySide::ZMin || side == BoundarySide::ZMax);
      if (!want_x && !want_z) { continue; }

      const int e1 = ftr->Elem1No;
      const FiniteElement *fe1 = fes.GetFE(e1);
      const int ndof = fe1->GetDof();

      Array<int> edofs1;
      fes.GetElementDofs(e1, edofs1);

      // Raw DOFs and element-canonical DOFs (Round-17).
      std::vector<CanonicalLocalDof> raw_dofs =
         BuildRawLocalDofs(*fe1, edofs1, I_data, ndof_total);
      std::vector<int> perm_elem(ndof);
      std::vector<CanonicalLocalDof> elem_canonical = [&]() {
         auto d = BuildCanonicalLocalDofs(*fe1, edofs1, I_data, ndof_total);
         for (int k = 0; k < ndof; k++) { perm_elem[k] = d[k].local_idx; }
         return d;
      }();

      // Face-canonical DOFs (Round-18A).
      std::vector<int> perm_face;
      std::vector<CanonicalLocalDof> face_canonical =
         BuildFaceLocalCanonicalLocalDofs(*fe1, edofs1, I_data, ndof_total,
                                           side, &perm_face);

      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                               2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         FaceLocalCanonicalSample s;
         s.face_idx = f;
         s.elem     = e1;
         s.q        = q;
         s.side     = side;
         s.cy = cy; s.cz = cz;
         for (int d = 0; d < 3; d++) { s.nor[d] = nor_vec(d); }

         s.elem_canonical_dofs = elem_canonical;
         s.face_canonical_dofs = face_canonical;

         s.shape_raw.assign(shape1.GetData(), shape1.GetData() + shape1.Size());
         s.shape_elem_canonical.resize(ndof);
         s.shape_face_canonical.resize(ndof);
         for (int k = 0; k < ndof; k++)
         {
            s.shape_elem_canonical[k] = shape1(perm_elem[k]);
            s.shape_face_canonical[k] = shape1(perm_face[k]);
         }

         real_t v_raw = 0.0, v_eCan = 0.0, v_fCan = 0.0;
         for (int i = 0; i < ndof; i++)
         {
            v_raw  += s.shape_raw[i]            * raw_dofs[i].value;
            v_eCan += s.shape_elem_canonical[i] * elem_canonical[i].value;
            v_fCan += s.shape_face_canonical[i] * face_canonical[i].value;
         }
         s.interp_raw             = v_raw;
         s.interp_elem_canonical  = v_eCan;
         s.interp_face_canonical  = v_fCan;

         if (want_x) { audit.x_faces.push_back(std::move(s)); }
         else        { audit.z_faces.push_back(std::move(s)); }
      }
   }
   return audit;
}

NonFaultFaceAudit RunADERNonFaultFaceAudit(WaveOperator<Mesh> &wave, Mesh &mesh,
                                           const BoundaryConfig &bc,
                                           const Vector &I_data,
                                           real_t dt)
{
   NonFaultFaceAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();
   audit.boundary_rhs_before_minv.SetSize(size);
   audit.interior_rhs_before_minv.SetSize(size);
   audit.boundary_rhs_before_minv = 0.0;
   audit.interior_rhs_before_minv = 0.0;

   // Round-11 Step 2: branch-isolation toggles.
   const bool skip_boundary          = EnvFlagSet("SEAS_TEST_SKIP_BOUNDARY");
   const bool skip_interior_nonfault = EnvFlagSet("SEAS_TEST_SKIP_INTERIOR_NONFAULT");

   real_t bulk_bg_scaled[NUM_STATE] = {0.0};
   if (const real_t *bulk_bg = wave.GetAbsorbingBackground())
   {
      for (int c = 0; c < NUM_STATE; c++) { bulk_bg_scaled[c] = dt * bulk_bg[c]; }
   }

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cy = 0.0;
      for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
      cy /= std::max(fv.Size(), 1);
      const bool is_boundary = (e2 < 0);
      const bool is_fault = (e2 >= 0) && (std::abs(cy - 0.5 * kL) < 1e-8);

      // Branch-isolation skips must happen at the FACE level (before
      // the quadrature loop opens) so the boundary/interior branches
      // contribute nothing to the audit totals when disabled.
      if (is_boundary && skip_boundary) { continue; }
      if (!is_boundary && !is_fault && skip_interior_nonfault) { continue; }

      const FiniteElement *fe1 = fes.GetFE(e1);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = e1 * wave.GetNDof();
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }
         const real_t w = ip.weight * nor_len;
         real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         real_t I_self[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               I_self[c] += shape1(i) * I_data(c * ndof_total + dof_offset1 + i);
            }
         }

         if (is_boundary)
         {
            real_t F_h[NUM_STATE];
            const int boundary_attr = 1;
            if (bc.absorbing_attrs.count(boundary_attr))
            {
               wave.GetFlux().AbsorbingTotal(nor, I_self, bulk_bg_scaled, F_h);
            }
            else if (bc.natural_attrs.count(boundary_attr))
            {
               if (wave.GetFreeSurfaceBCMode() == FreeSurfaceBCMode::Godunov)
               {
                  wave.GetFlux().FreeSurfaceGodunovTotal(nor, I_self, bulk_bg_scaled, F_h);
               }
               else
               {
                  wave.GetFlux().FreeSurfaceTotal(nor, I_self, bulk_bg_scaled, F_h);
               }
            }
            else
            {
               wave.GetFlux().AbsorbingTotal(nor, I_self, bulk_bg_scaled, F_h);
            }

            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  audit.boundary_rhs_before_minv(c * ndof_total + dof_offset1 + i) -=
                     w * shape1(i) * F_h[c];
               }
            }
         }
         else
         {
            if (is_fault) { continue; }

            const FiniteElement *fe2 = fes.GetFE(e2);
            const int dof_offset2 = e2 * wave.GetNDof();
            IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2);
            Vector shape2(ndof);
            fe2->CalcShape(ip2, shape2);

            real_t I_nbr[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               I_nbr[c] = 0.0;
               for (int i = 0; i < ndof; i++)
               {
                  I_nbr[c] += shape2(i) * I_data(c * ndof_total + dof_offset2 + i);
               }
            }

            real_t F_h[NUM_STATE];
            wave.GetFlux().Interior(nor, I_self, I_nbr, F_h);
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  audit.interior_rhs_before_minv(c * ndof_total + dof_offset1 + i) -=
                     w * shape1(i) * F_h[c];
                  audit.interior_rhs_before_minv(c * ndof_total + dof_offset2 + i) +=
                     w * shape2(i) * F_h[c];
               }
            }
         }
      }
   }

   audit.boundary_after_minv = audit.boundary_rhs_before_minv;
   audit.interior_after_minv = audit.interior_rhs_before_minv;
   ApplyMassInverseManually(wave, audit.boundary_after_minv);
   ApplyMassInverseManually(wave, audit.interior_after_minv);
   return audit;
}

// Phase 3 (§8.2): precomputed-path variant of RunADERNonFaultFaceAudit.
// Mirrors the structure of the runtime-path audit above but dispatches
// every non-fault, non-shared face through PrecomputedFaceFluxes's
// public AddInteriorFaceRhs / AddBoundaryFaceRhs.  Used by Gate 7' and
// Gate 14' to measure the precomputed-path lifted-rhs drift directly.
//
// Precondition: wave.UsePrecomputedFaceFluxes(true) has been called so
// wave.GetPrecomputedFaceFluxes().IsInitialized() is true.
NonFaultFaceAudit
RunPrecomputedFluxLiftedAudit(WaveOperator<Mesh> &wave, Mesh &mesh,
                               const Vector &I_data, real_t dt)
{
   NonFaultFaceAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();
   audit.boundary_rhs_before_minv.SetSize(size);
   audit.interior_rhs_before_minv.SetSize(size);
   audit.boundary_rhs_before_minv = 0.0;
   audit.interior_rhs_before_minv = 0.0;

   // Round-11 Step 2: branch-isolation toggles.
   const bool skip_boundary          = EnvFlagSet("SEAS_TEST_SKIP_BOUNDARY");
   const bool skip_interior_nonfault = EnvFlagSet("SEAS_TEST_SKIP_INTERIOR_NONFAULT");

   const PrecomputedFaceFluxes &pff = wave.GetPrecomputedFaceFluxes();
   MFEM_VERIFY(pff.IsInitialized(),
               "RunPrecomputedFluxLiftedAudit: PrecomputedFaceFluxes::Init "
               "has not run; call wave.UsePrecomputedFaceFluxes(true) first.");

   real_t bulk_bg_scaled[NUM_STATE] = {0.0};
   if (const real_t *bulk_bg = wave.GetAbsorbingBackground())
   {
      for (int c = 0; c < NUM_STATE; c++) { bulk_bg_scaled[c] = dt * bulk_bg[c]; }
   }

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cy = 0.0;
      for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
      cy /= std::max(fv.Size(), 1);
      const bool is_boundary = (e2 < 0);
      const bool is_fault = (e2 >= 0) && (std::abs(cy - 0.5 * kL) < 1e-8);

      // Skip fault faces: they have no entry in the precomputed table
      // (fault_face_set excludes them) and are handled by FaultFaceFlux.
      if (is_fault) { continue; }

      // Round-11 Step 2 branch-isolation skips.
      if (is_boundary && skip_boundary) { continue; }
      if (!is_boundary && skip_interior_nonfault) { continue; }

      const FiniteElement *fe1 = fes.GetFE(e1);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = e1 * wave.GetNDof();
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         const real_t w = ip.weight * nor_len;

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         real_t I_self[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               I_self[c] += shape1(i) * I_data(c * ndof_total + dof_offset1 + i);
            }
         }

         if (is_boundary)
         {
            pff.AddBoundaryFaceRhs(f, I_self, bulk_bg_scaled, w,
                                   shape1.GetData(),
                                   ndof, dof_offset1, ndof_total,
                                   audit.boundary_rhs_before_minv);
         }
         else
         {
            const FiniteElement *fe2 = fes.GetFE(e2);
            const int dof_offset2 = e2 * wave.GetNDof();
            IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2);
            Vector shape2(ndof);
            fe2->CalcShape(ip2, shape2);

            real_t I_nbr[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               I_nbr[c] = 0.0;
               for (int i = 0; i < ndof; i++)
               {
                  I_nbr[c] += shape2(i) * I_data(c * ndof_total + dof_offset2 + i);
               }
            }

            // §6.2 two-sided dispatch: one AddInteriorFaceRhs per side,
            // matching the production ComputeADERFaceFluxRHS path.
            pff.AddInteriorFaceRhs(f, /*caller_elem=*/e1, I_self, I_nbr,
                                   w, shape1.GetData(),
                                   ndof, dof_offset1, ndof_total,
                                   audit.interior_rhs_before_minv);
            pff.AddInteriorFaceRhs(f, /*caller_elem=*/e2, I_nbr, I_self,
                                   w, shape2.GetData(),
                                   ndof, dof_offset2, ndof_total,
                                   audit.interior_rhs_before_minv);
         }
      }
   }

   audit.boundary_after_minv = audit.boundary_rhs_before_minv;
   audit.interior_after_minv = audit.interior_rhs_before_minv;
   ApplyMassInverseManually(wave, audit.boundary_after_minv);
   ApplyMassInverseManually(wave, audit.interior_after_minv);
   return audit;
}

// Element-level lifted probes: replay the non-fault branches with
// symmetrizing modifiers, then apply M^{-1} and return the result so that
// the caller can compare orbit sorted-signatures on fault-adjacent tets.
// Interior-branch mode controls the Interior(nor, L, R) call; boundary-
// branch mode controls the FreeSurface[Godunov]Total call.
enum class IfaceLiftMode { Raw, Symmetrized };
enum class BFaceLiftMode { GammaRaw, GammaSymmetrized };

NonFaultFaceAudit
RunADERNonFaultFaceAuditModed(WaveOperator<Mesh> &wave, Mesh &mesh,
                              const BoundaryConfig &bc,
                              const Vector &I_data, real_t dt,
                              IfaceLiftMode iface_mode,
                              BFaceLiftMode bface_mode)
{
   NonFaultFaceAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();
   audit.boundary_rhs_before_minv.SetSize(size);
   audit.interior_rhs_before_minv.SetSize(size);
   audit.boundary_rhs_before_minv = 0.0;
   audit.interior_rhs_before_minv = 0.0;

   real_t bulk_bg_scaled[NUM_STATE] = {0.0};
   if (const real_t *bulk_bg = wave.GetAbsorbingBackground())
   {
      for (int c = 0; c < NUM_STATE; c++) { bulk_bg_scaled[c] = dt * bulk_bg[c]; }
   }

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cy = 0.0;
      for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
      cy /= std::max(fv.Size(), 1);
      const bool is_boundary = (e2 < 0);
      const bool is_fault = (e2 >= 0) && (std::abs(cy - 0.5 * kL) < 1e-8);

      const FiniteElement *fe1 = fes.GetFE(e1);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = e1 * wave.GetNDof();
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                               2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }
         const real_t w = ip.weight * nor_len;
         real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         real_t I_self[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               I_self[c] += shape1(i) * I_data(c * ndof_total + dof_offset1 + i);
            }
         }

         if (is_boundary)
         {
            real_t F_h[NUM_STATE];
            auto call_gamma = [&](const real_t n[3], real_t *out_f) {
               if (bc.natural_attrs.count(1) > 0)
               {
                  wave.GetFlux().FreeSurfaceTotal(n, I_self,
                                                  bulk_bg_scaled, out_f);
               }
               else
               {
                  wave.GetFlux().AbsorbingTotal(n, I_self, bulk_bg_scaled,
                                                 out_f);
               }
            };
            if (bface_mode == BFaceLiftMode::GammaSymmetrized)
            {
               real_t F_fwd[NUM_STATE], F_rev[NUM_STATE];
               call_gamma(nor, F_fwd);
               const real_t nor_rev[3] = {-nor[0], -nor[1], -nor[2]};
               call_gamma(nor_rev, F_rev);
               for (int c = 0; c < NUM_STATE; c++)
               {
                  F_h[c] = 0.5 * (F_fwd[c] + F_rev[c]);
               }
            }
            else
            {
               call_gamma(nor, F_h);
            }
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  audit.boundary_rhs_before_minv(c * ndof_total + dof_offset1 + i) -=
                     w * shape1(i) * F_h[c];
               }
            }
         }
         else
         {
            if (is_fault) { continue; }

            const FiniteElement *fe2 = fes.GetFE(e2);
            const int dof_offset2 = e2 * wave.GetNDof();
            IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2);
            Vector shape2(ndof);
            fe2->CalcShape(ip2, shape2);

            real_t I_nbr[NUM_STATE];
            for (int c = 0; c < NUM_STATE; c++)
            {
               I_nbr[c] = 0.0;
               for (int i = 0; i < ndof; i++)
               {
                  I_nbr[c] += shape2(i) * I_data(c * ndof_total + dof_offset2 + i);
               }
            }

            real_t F_h[NUM_STATE];
            if (iface_mode == IfaceLiftMode::Symmetrized)
            {
               real_t F_fwd[NUM_STATE], F_rev[NUM_STATE];
               wave.GetFlux().Interior(nor, I_self, I_nbr, F_fwd);
               const real_t nor_rev[3] = {-nor[0], -nor[1], -nor[2]};
               wave.GetFlux().Interior(nor_rev, I_nbr, I_self, F_rev);
               for (int c = 0; c < NUM_STATE; c++)
               {
                  F_h[c] = 0.5 * (F_fwd[c] + F_rev[c]);
               }
            }
            else
            {
               wave.GetFlux().Interior(nor, I_self, I_nbr, F_h);
            }
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  audit.interior_rhs_before_minv(c * ndof_total + dof_offset1 + i) -=
                     w * shape1(i) * F_h[c];
                  audit.interior_rhs_before_minv(c * ndof_total + dof_offset2 + i) +=
                     w * shape2(i) * F_h[c];
               }
            }
         }
      }
   }

   audit.boundary_after_minv = audit.boundary_rhs_before_minv;
   audit.interior_after_minv = audit.interior_rhs_before_minv;
   ApplyMassInverseManually(wave, audit.boundary_after_minv);
   ApplyMassInverseManually(wave, audit.interior_after_minv);
   return audit;
}

// Per-face samples for the second-step non-fault orbit drill-down probes.
struct FaceOrbitSample
{
   int face_idx = -1;
   int qp_idx = -1;
   int e1 = -1;
   int e2 = -1;
   real_t cx = 0.0, cy = 0.0, cz = 0.0;
   real_t qx = 0.0, qy = 0.0, qz = 0.0;
   real_t nor[3] = {0.0, 0.0, 0.0};
   real_t nor_len = 0.0;
   real_t w = 0.0;
   std::array<real_t, NUM_STATE> I_self{};
   std::array<real_t, NUM_STATE> I_nbr{};
   std::array<real_t, NUM_STATE> F_h{};
};

// Orbit classification for step-2 non-fault face probes.  The production
// `test_adjacent_triangle_fault_uniformity` asserts orbit identity across
// the y=L/2 reflection; orbits are defined by mirroring `(x, y, z)` →
// `(x, L-y, z)`.  Two faces sit in the same orbit when their centroids
// and normal axes match after reflecting either face across y=L/2.
std::string OrbitKey(real_t cx, real_t cy, real_t cz, const real_t nor[3])
{
   constexpr real_t s = 1e3;          // micron-scale coordinate quantization
   constexpr real_t n_quant = 1e4;    // ULP-free normal quantization
   const real_t abs_yc = std::abs(cy - 0.5 * kL);
   auto r = [&](real_t v, real_t q) { return std::round(v * q); };
   // Also fold the sign of the normal under the y=L/2 reflection so that
   // a face with `nor = +y_hat` matches its mirror at `L - y` with
   // `nor = -y_hat`; same for the y-component of non-axis-aligned faults.
   const real_t ny_abs = std::abs(nor[1]);
   char buf[256];
   std::snprintf(buf, sizeof(buf),
                 "x=%g;y=%g;z=%g;nx=%g;ny=%g;nz=%g",
                 r(std::abs(cx - 0.5 * kL), s),
                 r(abs_yc, s),
                 r(std::abs(cz - 0.5 * kL), s),
                 r(std::abs(nor[0]), n_quant),
                 r(ny_abs, n_quant),
                 r(std::abs(nor[2]), n_quant));
   return std::string(buf);
}

// Max per-component F_h sorted-signature drift within each orbit, taken
// over all orbits.  Operates on a bare vector of face samples.
real_t MaxOrbitFhDrift(const std::vector<FaceOrbitSample> &samples,
                       int &worst_comp,
                       std::string &worst_orbit)
{
   worst_comp = -1;
   worst_orbit.clear();
   std::map<std::string, std::vector<const FaceOrbitSample *>> by_orbit;
   for (const auto &s : samples)
   {
      const std::string key = OrbitKey(s.cx, s.cy, s.cz, s.nor);
      by_orbit[key].push_back(&s);
   }
   real_t worst = 0.0;
   for (const auto &kv : by_orbit)
   {
      const auto &orbit = kv.second;
      if (orbit.size() < 2) { continue; }
      for (int c = 0; c < NUM_STATE; c++)
      {
         std::vector<real_t> vals;
         vals.reserve(orbit.size());
         for (const auto *p : orbit) { vals.push_back(p->F_h[c]); }
         std::sort(vals.begin(), vals.end());
         const real_t scale =
            std::max({std::abs(vals.front()), std::abs(vals.back()), real_t(1.0)});
         const real_t drift = (vals.back() - vals.front()) / scale;
         if (drift > worst)
         {
            worst = drift;
            worst_comp = c;
            worst_orbit = kv.first;
         }
      }
   }
   return worst;
}

enum class IfaceProbeMode
{
   Raw,                 // flux_.Interior(nor, I_self, I_nbr)
   Symmetrized,         // 0.5 * (Interior(nor, L, R) + Interior(-nor, R, L))
   CanonicalNormal      // replace nor with outward-positive-y canonical
};

enum class BFaceProbeMode
{
   GammaBaseline,       // FreeSurfaceTotal(nor, I_self, bg)  [matches prod]
   GodunovBaseline,     // FreeSurfaceGodunovTotal(nor, I_self, bg)
   GammaSymmetrized     // 0.5 * (Gamma(nor,...) + Gamma(-nor,...))
};

std::vector<FaceOrbitSample>
CollectInteriorFaceSamples(WaveOperator<Mesh> &wave, Mesh &mesh,
                           const Vector &I, real_t /*dt*/,
                           IfaceProbeMode mode)
{
   std::vector<FaceOrbitSample> out;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr || ftr->Elem2No < 0) { continue; }

      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < fv.Size(); v++)
      {
         cx += mesh.GetVertex(fv[v])[0];
         cy += mesh.GetVertex(fv[v])[1];
         cz += mesh.GetVertex(fv[v])[2];
      }
      cx /= fv.Size(); cy /= fv.Size(); cz /= fv.Size();
      if (std::abs(cy - 0.5 * kL) < 1e-8) { continue; }  // skip fault faces

      const FiniteElement *fe1 = fes.GetFE(ftr->Elem1No);
      const FiniteElement *fe2 = fes.GetFE(ftr->Elem2No);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = ftr->Elem1No * wave.GetNDof();
      const int dof_offset2 = ftr->Elem2No * wave.GetNDof();
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                               2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }
         const real_t w = ip.weight * nor_len;
         real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         if (mode == IfaceProbeMode::CanonicalNormal)
         {
            // Flip nor so that the dominant axis points in the +direction
            // (orbit-invariant across y=L/2 reflection only when the
            // non-y components are canonicalized too).
            int dom = 0;
            for (int d = 1; d < 3; d++)
            {
               if (std::abs(nor[d]) > std::abs(nor[dom])) { dom = d; }
            }
            if (nor[dom] < 0.0)
            {
               for (int d = 0; d < 3; d++) { nor[d] = -nor[d]; }
            }
         }

         IntegrationPoint ip1, ip2;
         ftr->Loc1.Transform(ip, ip1);
         ftr->Loc2.Transform(ip, ip2);
         Vector shape1(ndof), shape2(ndof);
         fe1->CalcShape(ip1, shape1);
         fe2->CalcShape(ip2, shape2);

         real_t I_self[NUM_STATE], I_nbr[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self[c] = 0.0;
            I_nbr[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               I_self[c] += shape1(i) * I(c * ndof_total + dof_offset1 + i);
               I_nbr[c]  += shape2(i) * I(c * ndof_total + dof_offset2 + i);
            }
         }

         real_t F_h[NUM_STATE];
         if (mode == IfaceProbeMode::Symmetrized)
         {
            real_t F_fwd[NUM_STATE], F_rev[NUM_STATE];
            wave.GetFlux().Interior(nor, I_self, I_nbr, F_fwd);
            real_t nor_rev[3] = {-nor[0], -nor[1], -nor[2]};
            wave.GetFlux().Interior(nor_rev, I_nbr, I_self, F_rev);
            for (int c = 0; c < NUM_STATE; c++)
            {
               F_h[c] = 0.5 * (F_fwd[c] + F_rev[c]);
            }
         }
         else
         {
            wave.GetFlux().Interior(nor, I_self, I_nbr, F_h);
         }

         FaceOrbitSample s;
         s.face_idx = f;
         s.qp_idx = q;
         s.e1 = ftr->Elem1No;
         s.e2 = ftr->Elem2No;
         s.cx = cx; s.cy = cy; s.cz = cz;
         s.qx = phys(0); s.qy = phys(1); s.qz = phys(2);
         for (int d = 0; d < 3; d++) { s.nor[d] = nor[d]; }
         s.nor_len = nor_len;
         s.w = w;
         for (int c = 0; c < NUM_STATE; c++)
         {
            s.I_self[c] = I_self[c];
            s.I_nbr[c]  = I_nbr[c];
            s.F_h[c]    = F_h[c];
         }
         out.push_back(s);
      }
   }
   return out;
}

std::vector<FaceOrbitSample>
CollectBoundaryFaceSamples(WaveOperator<Mesh> &wave, Mesh &mesh,
                           const BoundaryConfig &bc,
                           const Vector &I, real_t dt,
                           BFaceProbeMode mode)
{
   std::vector<FaceOrbitSample> out;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   real_t bulk_bg_scaled[NUM_STATE] = {0.0};
   if (const real_t *bulk_bg = wave.GetAbsorbingBackground())
   {
      for (int c = 0; c < NUM_STATE; c++) { bulk_bg_scaled[c] = dt * bulk_bg[c]; }
   }

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr || ftr->Elem2No >= 0) { continue; }  // interior or fault

      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int v = 0; v < fv.Size(); v++)
      {
         cx += mesh.GetVertex(fv[v])[0];
         cy += mesh.GetVertex(fv[v])[1];
         cz += mesh.GetVertex(fv[v])[2];
      }
      cx /= fv.Size(); cy /= fv.Size(); cz /= fv.Size();

      const FiniteElement *fe1 = fes.GetFE(ftr->Elem1No);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = ftr->Elem1No * wave.GetNDof();
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                               2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }
         const real_t w = ip.weight * nor_len;
         real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         real_t I_self[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               I_self[c] += shape1(i) * I(c * ndof_total + dof_offset1 + i);
            }
         }

         // Fixture uses bc.natural_attrs = {1}: all exterior faces map to
         // FreeSurface[Godunov]Total via ClassifyBoundaryFace.
         const int boundary_attr = 1;
         const bool is_natural = bc.natural_attrs.count(boundary_attr) > 0;

         real_t F_h[NUM_STATE];
         auto call_gamma = [&](const real_t n[3], real_t *out_f) {
            wave.GetFlux().FreeSurfaceTotal(n, I_self, bulk_bg_scaled, out_f);
         };
         auto call_godunov = [&](const real_t n[3], real_t *out_f) {
            wave.GetFlux().FreeSurfaceGodunovTotal(n, I_self, bulk_bg_scaled,
                                                   out_f);
         };
         if (!is_natural)
         {
            // AbsorbingTotal fallback (not exercised by this fixture but
            // kept so the helper is reusable).
            wave.GetFlux().AbsorbingTotal(nor, I_self, bulk_bg_scaled, F_h);
         }
         else if (mode == BFaceProbeMode::GodunovBaseline)
         {
            call_godunov(nor, F_h);
         }
         else if (mode == BFaceProbeMode::GammaSymmetrized)
         {
            real_t F_fwd[NUM_STATE], F_rev[NUM_STATE];
            call_gamma(nor, F_fwd);
            real_t nor_rev[3] = {-nor[0], -nor[1], -nor[2]};
            call_gamma(nor_rev, F_rev);
            for (int c = 0; c < NUM_STATE; c++)
            {
               F_h[c] = 0.5 * (F_fwd[c] + F_rev[c]);
            }
         }
         else
         {
            call_gamma(nor, F_h);
         }

         FaceOrbitSample s;
         s.face_idx = f;
         s.qp_idx = q;
         s.e1 = ftr->Elem1No;
         s.e2 = -1;
         s.cx = cx; s.cy = cy; s.cz = cz;
         s.qx = phys(0); s.qy = phys(1); s.qz = phys(2);
         for (int d = 0; d < 3; d++) { s.nor[d] = nor[d]; }
         s.nor_len = nor_len;
         s.w = w;
         for (int c = 0; c < NUM_STATE; c++)
         {
            s.I_self[c] = I_self[c];
            s.I_nbr[c]  = 0.0;
            s.F_h[c]    = F_h[c];
         }
         out.push_back(s);
      }
   }
   return out;
}

AffineTraceProbeResult RunAffineTraceProbe(WaveOperator<Mesh> &wave,
                                           Mesh &mesh,
                                           const BoundaryConfig &bc)
{
   (void)bc;
   AffineTraceProbeResult result;
   const auto &fes = wave.GetFESpace();
   for (AffineFieldKind field : {AffineFieldKind::X,
                                 AffineFieldKind::Z,
                                 AffineFieldKind::XPlusZ})
   {
      Vector Q_affine;
      FillQFromAffineField(fes, SXZ, field, Q_affine);

      const auto iface = CollectInteriorFaceSamples(wave, mesh, Q_affine, 0.0,
                                                    IfaceProbeMode::Raw);
      for (const auto &s : iface)
      {
         const real_t exact = EvaluateAffineField(field, s.qx, s.qy, s.qz);
         const real_t scale = std::max(std::abs(exact), real_t(1.0));
         const real_t self_err = std::abs(s.I_self[SXZ] - exact) / scale;
         const real_t nbr_err  = std::abs(s.I_nbr[SXZ]  - exact) / scale;
         const real_t jump_err = std::abs(s.I_self[SXZ] - s.I_nbr[SXZ]) / scale;
         if (self_err > result.worst_self_err)
         {
            result.worst_self_err = self_err;
            result.worst_self_field = field;
         }
         if (nbr_err > result.worst_nbr_err)
         {
            result.worst_nbr_err = nbr_err;
            result.worst_nbr_field = field;
         }
         if (jump_err > result.worst_jump_err)
         {
            result.worst_jump_err = jump_err;
            result.worst_jump_field = field;
         }
      }

      const auto bface = CollectBoundaryFaceSamples(wave, mesh, bc, Q_affine,
                                                    kDt,
                                                    BFaceProbeMode::GammaBaseline);
      for (const auto &s : bface)
      {
         const real_t exact = EvaluateAffineField(field, s.qx, s.qy, s.qz);
         const real_t scale = std::max(std::abs(exact), real_t(1.0));
         const real_t self_err = std::abs(s.I_self[SXZ] - exact) / scale;
         if (self_err > result.worst_bdry_err)
         {
            result.worst_bdry_err = self_err;
            result.worst_bdry_field = field;
         }
      }
   }
   return result;
}

NonFaultFaceAudit RunAnalyticTraceLiftAudit(WaveOperator<Mesh> &wave,
                                            Mesh &mesh,
                                            const BoundaryConfig &bc,
                                            real_t dt,
                                            int component,
                                            AffineFieldKind field)
{
   NonFaultFaceAudit audit;
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();
   audit.boundary_rhs_before_minv.SetSize(size);
   audit.interior_rhs_before_minv.SetSize(size);
   audit.boundary_rhs_before_minv = 0.0;
   audit.interior_rhs_before_minv = 0.0;

   // Round-11 Step 2: branch-isolation toggles.
   const bool skip_boundary          = EnvFlagSet("SEAS_TEST_SKIP_BOUNDARY");
   const bool skip_interior_nonfault = EnvFlagSet("SEAS_TEST_SKIP_INTERIOR_NONFAULT");

   real_t bulk_bg_scaled[NUM_STATE] = {0.0};
   if (const real_t *bulk_bg = wave.GetAbsorbingBackground())
   {
      for (int c = 0; c < NUM_STATE; c++) { bulk_bg_scaled[c] = dt * bulk_bg[c]; }
   }

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }

      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cy = 0.0;
      for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
      cy /= std::max(fv.Size(), 1);
      const bool is_boundary = (e2 < 0);
      const bool is_fault = (e2 >= 0) && (std::abs(cy - 0.5 * kL) < 1e-8);

      // Round-11 Step 2 branch-isolation skips.
      if (is_boundary && skip_boundary) { continue; }
      if (!is_boundary && !is_fault && skip_interior_nonfault) { continue; }

      const FiniteElement *fe1 = fes.GetFE(e1);
      const int ndof = fe1->GetDof();
      const int dof_offset1 = e1 * wave.GetNDof();
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                               2 * wave.GetOrder());

      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }
         const real_t w = ip.weight * nor_len;
         real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);

         real_t I_self[NUM_STATE] = {0.0};
         real_t I_nbr[NUM_STATE] = {0.0};
         const real_t exact = EvaluateAffineField(field, phys(0), phys(1), phys(2));
         I_self[component] = exact;
         I_nbr[component] = exact;

         if (is_boundary)
         {
            real_t F_h[NUM_STATE];
            if (bc.natural_attrs.count(1) > 0)
            {
               wave.GetFlux().FreeSurfaceTotal(nor, I_self, bulk_bg_scaled, F_h);
            }
            else
            {
               wave.GetFlux().AbsorbingTotal(nor, I_self, bulk_bg_scaled, F_h);
            }
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  audit.boundary_rhs_before_minv(c * ndof_total + dof_offset1 + i) -=
                     w * shape1(i) * F_h[c];
               }
            }
         }
         else
         {
            if (is_fault) { continue; }

            const FiniteElement *fe2 = fes.GetFE(e2);
            const int dof_offset2 = e2 * wave.GetNDof();
            IntegrationPoint ip2;
            ftr->Loc2.Transform(ip, ip2);
            Vector shape2(ndof);
            fe2->CalcShape(ip2, shape2);

            real_t F_h[NUM_STATE];
            wave.GetFlux().Interior(nor, I_self, I_nbr, F_h);
            for (int c = 0; c < NUM_STATE; c++)
            {
               for (int i = 0; i < ndof; i++)
               {
                  audit.interior_rhs_before_minv(c * ndof_total + dof_offset1 + i) -=
                     w * shape1(i) * F_h[c];
                  audit.interior_rhs_before_minv(c * ndof_total + dof_offset2 + i) +=
                     w * shape2(i) * F_h[c];
               }
            }
         }
      }
   }

   audit.boundary_after_minv = audit.boundary_rhs_before_minv;
   audit.interior_after_minv = audit.interior_rhs_before_minv;
   ApplyMassInverseManually(wave, audit.boundary_after_minv);
   ApplyMassInverseManually(wave, audit.interior_after_minv);
   return audit;
}

AnalyticTraceLiftResult RunAnalyticTraceLiftSweep(WaveOperator<Mesh> &wave,
                                                  Mesh &mesh,
                                                  const BoundaryConfig &bc,
                                                  const std::vector<ElementMeanData> &fault_adjacent)
{
   AnalyticTraceLiftResult result;
   const int ndof_total = wave.GetFESpace().GetNDofs();
   auto sig = [&](const Vector &v, int comp, bool upper) {
      return MaxSortedSignatureDrift(v, fault_adjacent, ndof_total, comp, upper);
   };

   for (AffineFieldKind field : {AffineFieldKind::X,
                                 AffineFieldKind::Z,
                                 AffineFieldKind::XPlusZ})
   {
      const auto audit = RunAnalyticTraceLiftAudit(wave, mesh, bc, kDt, SXZ, field);
      real_t iface_worst = 0.0;
      real_t bface_worst = 0.0;
      for (int comp : {SXY, SXZ, SXX, SYY, SZZ})
      {
         for (bool upper : {false, true})
         {
            iface_worst = std::max(iface_worst, sig(audit.interior_after_minv, comp, upper));
            bface_worst = std::max(bface_worst, sig(audit.boundary_after_minv, comp, upper));
         }
      }
      if (iface_worst > result.worst_iface_drift)
      {
         result.worst_iface_drift = iface_worst;
         result.worst_iface_field = field;
      }
      if (bface_worst > result.worst_bface_drift)
      {
         result.worst_bface_drift = bface_worst;
         result.worst_bface_field = field;
      }
   }
   return result;
}

// Report helper: print orbit-by-orbit F_h drift breakdown when a gate fails.
void ReportOrbitBreakdown(const std::string &label,
                          const std::vector<FaceOrbitSample> &samples)
{
   std::map<std::string, std::vector<const FaceOrbitSample *>> by_orbit;
   for (const auto &s : samples)
   {
      by_orbit[OrbitKey(s.cx, s.cy, s.cz, s.nor)].push_back(&s);
   }
   std::cout << "      orbits in " << label << " probe: "
             << by_orbit.size() << "  (total samples: "
             << samples.size() << ")\n";
   int shown = 0;
   for (const auto &kv : by_orbit)
   {
      if (kv.second.size() < 2) { continue; }
      real_t orbit_worst = 0.0;
      int worst_c = -1;
      for (int c = 0; c < NUM_STATE; c++)
      {
         std::vector<real_t> vals;
         for (const auto *p : kv.second) { vals.push_back(p->F_h[c]); }
         std::sort(vals.begin(), vals.end());
         const real_t scale = std::max({std::abs(vals.front()),
                                        std::abs(vals.back()),
                                        real_t(1.0)});
         const real_t drift = (vals.back() - vals.front()) / scale;
         if (drift > orbit_worst) { orbit_worst = drift; worst_c = c; }
      }
      if (orbit_worst > 1e-10 && shown < 4)
      {
         std::cout << "        orbit " << kv.first
                   << "  size=" << kv.second.size()
                   << "  worst F_h drift=" << std::scientific
                   << std::setprecision(3) << orbit_worst
                   << " on " << kCompName[worst_c] << "\n";
         ++shown;
      }
   }
}

} // namespace

int main()
{
   ReadMeshDimsFromEnv();
   std::cout << "\n=== TPV102 adjacent-triangle first-step audit ===\n";
   std::cout << "  mesh dims: " << g_nx << " x " << g_ny << " x " << g_nz
             << " hexes (override with SEAS_TEST_FAULT_NX/NY/NZ, even integers)\n";

   Mesh mesh = BuildCartesianFaultMesh();
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};

   real_t bulk_bg[NUM_STATE] = {0.0};

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                           TPV102Params::mu, TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bulk_bg);

   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   const int n_fault = SetupFault(wave, mesh, kOrder, dof_data, ff, fault_coords);
   std::cout << "  fault QPs: " << n_fault << "\n";
   MFEM_VERIFY(n_fault > 0, "audit fixture must have fault QPs");

   FaultOnlyAudit audit = RunFaultOnlyAudit(wave, mesh, ff, dof_data, kDt);

   std::cout << "\n-- Gate 1: step-0 fault-local uniformity --\n";
   int worst_idx = -1;
   const real_t slip_spread = MaxFaultStateSpread(dof_data, 0, worst_idx);
   const real_t tau1_spread = MaxFaultStateSpread(dof_data, 1, worst_idx);
   const real_t tau2_spread = MaxFaultStateSpread(dof_data, 2, worst_idx);
   const real_t sn_spread = MaxFaultStateSpread(dof_data, 3, worst_idx);
   int worst_comp = -1;
   const real_t qimp_plus_spread =
      MaxRelativeSpread(audit.samples, &FaultQPSample::q_imp_plus_local, worst_comp);
   const real_t qimp_minus_spread =
      MaxRelativeSpread(audit.samples, &FaultQPSample::q_imp_minus_local, worst_comp);
   const real_t fh_plus_spread =
      MaxRelativeSpread(audit.samples, &FaultQPSample::f_h_plus, worst_comp);
   const real_t fh_minus_spread =
      MaxRelativeSpread(audit.samples, &FaultQPSample::f_h_minus, worst_comp);

   TEST_LE(slip_spread, 1.0e-10, "fault slip_rate uniform across all QPs at step 0");
   TEST_LE(tau1_spread, 1.0e-10, "fault tau1_corr uniform across all QPs at step 0");
   TEST_LE(tau2_spread, 1.0e-10, "fault tau2_corr uniform across all QPs at step 0");
   TEST_LE(sn_spread, 1.0e-10, "fault sigma_n_corr uniform across all QPs at step 0");
   TEST_LE(qimp_plus_spread, 1.0e-10, "Q_imp_plus_local uniform across all fault QPs");
   TEST_LE(qimp_minus_spread, 1.0e-10, "Q_imp_minus_local uniform across all fault QPs");
   TEST_LE(fh_plus_spread, 1.0e-10, "F_h_plus uniform across all fault QPs");
   TEST_LE(fh_minus_spread, 1.0e-10, "F_h_minus uniform across all fault QPs");

   const auto fault_adjacent = BuildFaultAdjacentElementData(wave, mesh);
   const int ndof_total = wave.GetFESpace().GetNDofs();

   std::cout << "\n-- Gate 2: fault-face lift symmetry on adjacent bulk tets --\n";
   real_t min_lower = 0.0, max_lower = 0.0, min_upper = 0.0, max_upper = 0.0;
   const real_t rhs_sxy_lower = CheckElementSideUniformity(
      audit.rhs_before_minv, fault_adjacent, ndof_total, SXY, false, min_lower, max_lower);
   const real_t rhs_sxy_upper = CheckElementSideUniformity(
      audit.rhs_before_minv, fault_adjacent, ndof_total, SXY, true, min_upper, max_upper);
   const real_t rhs_sxz_lower = CheckElementSideUniformity(
      audit.rhs_before_minv, fault_adjacent, ndof_total, SXZ, false, min_lower, max_lower);
   const real_t rhs_sxz_upper = CheckElementSideUniformity(
      audit.rhs_before_minv, fault_adjacent, ndof_total, SXZ, true, min_upper, max_upper);
   const real_t q_sxy_lower = CheckElementSideUniformity(
      audit.q_after_minv, fault_adjacent, ndof_total, SXY, false, min_lower, max_lower);
   const real_t q_sxy_upper = CheckElementSideUniformity(
      audit.q_after_minv, fault_adjacent, ndof_total, SXY, true, min_upper, max_upper);
   const real_t q_sxz_lower = CheckElementSideUniformity(
      audit.q_after_minv, fault_adjacent, ndof_total, SXZ, false, min_lower, max_lower);
   const real_t q_sxz_upper = CheckElementSideUniformity(
      audit.q_after_minv, fault_adjacent, ndof_total, SXZ, true, min_upper, max_upper);

   TEST_LE(rhs_sxy_lower, 1.0e-10, "pre-Minv rhs SXY cell means uniform on lower-side fault-adjacent tets");
   TEST_LE(rhs_sxy_upper, 1.0e-10, "pre-Minv rhs SXY cell means uniform on upper-side fault-adjacent tets");
   TEST_LE(rhs_sxz_lower, 1.0e-10, "pre-Minv rhs SXZ cell means uniform on lower-side fault-adjacent tets");
   TEST_LE(rhs_sxz_upper, 1.0e-10, "pre-Minv rhs SXZ cell means uniform on upper-side fault-adjacent tets");
   TEST_LE(q_sxy_lower, 1.0e-10, "post-Minv SXY cell means uniform on lower-side fault-adjacent tets");
   TEST_LE(q_sxy_upper, 1.0e-10, "post-Minv SXY cell means uniform on upper-side fault-adjacent tets");
   TEST_LE(q_sxz_lower, 1.0e-10, "post-Minv SXZ cell means uniform on lower-side fault-adjacent tets");
   TEST_LE(q_sxz_upper, 1.0e-10, "post-Minv SXZ cell means uniform on upper-side fault-adjacent tets");

   std::cout << "\n-- Gate 3: manual fault-only increment equals production first ADER step --\n";
   ResetFaultState(dof_data, fault_coords);
   wave.SetFaultDOFData(&dof_data, wave.GetNbfPerFace());
   Vector Q0(wave.Height());
   Vector Q1(wave.Height());
   Q0 = 0.0;
   wave.AdvanceADER(Q0, kDt, kAderOrder, Q1);
   Q1 -= audit.q_after_minv;
   const real_t max_abs_diff = Q1.Normlinf();
   TEST_LE(max_abs_diff, 1.0e-11,
           "manual fault-face-only first-step increment matches production AdvanceADER(Q=0)");

   // Phase 3 (§8.2): opt-in to the precomputed-flux path for Gates 4+.
   // Env: SEAS_TEST_USE_PRECOMPUTED_FLUX=1 (or any nonzero value).
   // When set, `wave.UsePrecomputedFaceFluxes(true)` is called here —
   // BEFORE Gate 4 — so every subsequent `wave.AdvanceADER` call
   // dispatches non-fault faces through PrecomputedFaceFluxes per
   // §6.2.  Gates 1-3 above still ran on the runtime path (default);
   // they are the R-001 regression guard (fault branch untouched).
   bool use_precomputed_env = false;
   if (const char *s = std::getenv("SEAS_TEST_USE_PRECOMPUTED_FLUX"))
   {
      use_precomputed_env = (std::strlen(s) > 0 && std::strcmp(s, "0") != 0);
   }
   if (use_precomputed_env)
   {
      wave.UsePrecomputedFaceFluxes(true);
      std::cout << "\n[Phase 3 §8.2] SEAS_TEST_USE_PRECOMPUTED_FLUX=1 — "
                << "wave.UsePrecomputedFaceFluxes(true) fired before Gate 4.\n";
   }

   std::cout << "\n-- Gate 4: second-step bulk signature audit --\n";
   ResetFaultState(dof_data, fault_coords);
   wave.SetFaultDOFData(&dof_data, wave.GetNbfPerFace());
   Vector Q_step0(wave.Height());
   Vector Q_step1(wave.Height());
    Vector I_step1(wave.Height());
   Vector volume_rhs(wave.Height());
   Vector volume_update(wave.Height());
   Vector full_step1_increment(wave.Height());
   Vector nonfault_residual(wave.Height());
   Vector nonfault_face_manual(wave.Height());
   Q0 = 0.0;
   wave.AdvanceADER(Q0, kDt, kAderOrder, Q_step0);
   wave.ComputeADERTimeIntegrated(Q_step0, kDt, kAderOrder, I_step1);
   FaultOnlyAudit step1_face_audit =
      RunADERTraceAudit(wave, mesh, ff, dof_data, I_step1, kDt);
   NonFaultFaceAudit step1_nonfault_face_audit =
      RunADERNonFaultFaceAudit(wave, mesh, bc, I_step1, kDt);
   volume_rhs = 0.0;
   wave.ComputeADERVolumeUpdate(I_step1, volume_rhs);
   volume_update = volume_rhs;
   ApplyMassInverseManually(wave, volume_update);
   wave.AdvanceADER(Q_step0, kDt, kAderOrder, Q_step1);
   full_step1_increment = Q_step1;
   full_step1_increment -= Q_step0;
   nonfault_residual = full_step1_increment;
   nonfault_residual -= step1_face_audit.q_after_minv;
   nonfault_residual -= volume_update;
   nonfault_face_manual = step1_nonfault_face_audit.boundary_after_minv;
   nonfault_face_manual += step1_nonfault_face_audit.interior_after_minv;

   const real_t step1_slip_spread = MaxFaultStateSpread(dof_data, 0, worst_idx);
   const real_t step1_tau1_spread = MaxFaultStateSpread(dof_data, 1, worst_idx);
   const real_t step1_tau2_spread = MaxFaultStateSpread(dof_data, 2, worst_idx);
   const real_t step1_sn_spread = MaxFaultStateSpread(dof_data, 3, worst_idx);
   const real_t step1_i_plus_spread =
      MaxRelativeSpread(step1_face_audit.samples, &FaultQPSample::i_plus_local, worst_comp);
   const real_t step1_i_minus_spread =
      MaxRelativeSpread(step1_face_audit.samples, &FaultQPSample::i_minus_local, worst_comp);
   const real_t step1_qimp_plus_spread =
      MaxRelativeSpread(step1_face_audit.samples, &FaultQPSample::q_imp_plus_local, worst_comp);
   const real_t step1_qimp_minus_spread =
      MaxRelativeSpread(step1_face_audit.samples, &FaultQPSample::q_imp_minus_local, worst_comp);
   const real_t step1_fh_plus_spread =
      MaxRelativeSpread(step1_face_audit.samples, &FaultQPSample::f_h_plus, worst_comp);
   const real_t step1_fh_minus_spread =
      MaxRelativeSpread(step1_face_audit.samples, &FaultQPSample::f_h_minus, worst_comp);

   TEST_LE(step1_i_plus_spread, 1.0e-10, "second-step I_plus_local traces uniform across all fault QPs");
   TEST_LE(step1_i_minus_spread, 1.0e-10, "second-step I_minus_local traces uniform across all fault QPs");
   TEST_LE(step1_qimp_plus_spread, 1.0e-10, "second-step Q_imp_plus_local uniform across all fault QPs");
   TEST_LE(step1_qimp_minus_spread, 1.0e-10, "second-step Q_imp_minus_local uniform across all fault QPs");
   TEST_LE(step1_fh_plus_spread, 1.0e-10, "second-step F_h_plus uniform across all fault QPs");
   TEST_LE(step1_fh_minus_spread, 1.0e-10, "second-step F_h_minus uniform across all fault QPs");
   TEST_LE(step1_slip_spread, 1.0e-10, "fault slip_rate still uniform after the second ADER step");
   TEST_LE(step1_tau1_spread, 1.0e-10, "fault tau1_corr still uniform after the second ADER step");
   TEST_LE(step1_tau2_spread, 1.0e-10, "fault tau2_corr still uniform after the second ADER step");
   TEST_LE(step1_sn_spread, 1.0e-10, "fault sigma_n_corr still uniform after the second ADER step");

   const real_t step0_sxy_lower_sig =
      MaxSortedSignatureDrift(Q_step0, fault_adjacent, ndof_total, SXY, false);
   const real_t step0_sxy_upper_sig =
      MaxSortedSignatureDrift(Q_step0, fault_adjacent, ndof_total, SXY, true);
   const real_t step0_sxz_lower_sig =
      MaxSortedSignatureDrift(Q_step0, fault_adjacent, ndof_total, SXZ, false);
   const real_t step0_sxz_upper_sig =
      MaxSortedSignatureDrift(Q_step0, fault_adjacent, ndof_total, SXZ, true);
   const real_t step1_sxy_lower_sig =
      MaxSortedSignatureDrift(Q_step1, fault_adjacent, ndof_total, SXY, false);
   const real_t step1_sxy_upper_sig =
      MaxSortedSignatureDrift(Q_step1, fault_adjacent, ndof_total, SXY, true);
   const real_t step1_sxz_lower_sig =
      MaxSortedSignatureDrift(Q_step1, fault_adjacent, ndof_total, SXZ, false);
   const real_t step1_sxz_upper_sig =
      MaxSortedSignatureDrift(Q_step1, fault_adjacent, ndof_total, SXZ, true);
   const real_t step1_fault_only_sxy_lower_sig =
      MaxSortedSignatureDrift(step1_face_audit.q_after_minv, fault_adjacent, ndof_total, SXY, false);
   const real_t step1_fault_only_sxy_upper_sig =
      MaxSortedSignatureDrift(step1_face_audit.q_after_minv, fault_adjacent, ndof_total, SXY, true);
   const real_t step1_fault_only_sxz_lower_sig =
      MaxSortedSignatureDrift(step1_face_audit.q_after_minv, fault_adjacent, ndof_total, SXZ, false);
   const real_t step1_fault_only_sxz_upper_sig =
      MaxSortedSignatureDrift(step1_face_audit.q_after_minv, fault_adjacent, ndof_total, SXZ, true);
   const real_t step1_volume_only_sxy_lower_sig =
      MaxSortedSignatureDrift(volume_update, fault_adjacent, ndof_total, SXY, false);
   const real_t step1_volume_only_sxy_upper_sig =
      MaxSortedSignatureDrift(volume_update, fault_adjacent, ndof_total, SXY, true);
   const real_t step1_volume_only_sxz_lower_sig =
      MaxSortedSignatureDrift(volume_update, fault_adjacent, ndof_total, SXZ, false);
   const real_t step1_volume_only_sxz_upper_sig =
      MaxSortedSignatureDrift(volume_update, fault_adjacent, ndof_total, SXZ, true);
   const real_t step1_nonfault_resid_sxy_lower_sig =
      MaxSortedSignatureDrift(nonfault_residual, fault_adjacent, ndof_total, SXY, false);
   const real_t step1_nonfault_resid_sxy_upper_sig =
      MaxSortedSignatureDrift(nonfault_residual, fault_adjacent, ndof_total, SXY, true);
   const real_t step1_nonfault_resid_sxz_lower_sig =
      MaxSortedSignatureDrift(nonfault_residual, fault_adjacent, ndof_total, SXZ, false);
   const real_t step1_nonfault_resid_sxz_upper_sig =
      MaxSortedSignatureDrift(nonfault_residual, fault_adjacent, ndof_total, SXZ, true);
   const real_t step1_boundary_sxy_lower_sig =
      MaxSortedSignatureDrift(step1_nonfault_face_audit.boundary_after_minv, fault_adjacent, ndof_total, SXY, false);
   const real_t step1_boundary_sxy_upper_sig =
      MaxSortedSignatureDrift(step1_nonfault_face_audit.boundary_after_minv, fault_adjacent, ndof_total, SXY, true);
   const real_t step1_boundary_sxz_lower_sig =
      MaxSortedSignatureDrift(step1_nonfault_face_audit.boundary_after_minv, fault_adjacent, ndof_total, SXZ, false);
   const real_t step1_boundary_sxz_upper_sig =
      MaxSortedSignatureDrift(step1_nonfault_face_audit.boundary_after_minv, fault_adjacent, ndof_total, SXZ, true);
   const real_t step1_interior_nf_sxy_lower_sig =
      MaxSortedSignatureDrift(step1_nonfault_face_audit.interior_after_minv, fault_adjacent, ndof_total, SXY, false);
   const real_t step1_interior_nf_sxy_upper_sig =
      MaxSortedSignatureDrift(step1_nonfault_face_audit.interior_after_minv, fault_adjacent, ndof_total, SXY, true);
   const real_t step1_interior_nf_sxz_lower_sig =
      MaxSortedSignatureDrift(step1_nonfault_face_audit.interior_after_minv, fault_adjacent, ndof_total, SXZ, false);
   const real_t step1_interior_nf_sxz_upper_sig =
      MaxSortedSignatureDrift(step1_nonfault_face_audit.interior_after_minv, fault_adjacent, ndof_total, SXZ, true);
   const real_t step1_nonfault_manual_sxy_lower_sig =
      MaxSortedSignatureDrift(nonfault_face_manual, fault_adjacent, ndof_total, SXY, false);
   const real_t step1_nonfault_manual_sxy_upper_sig =
      MaxSortedSignatureDrift(nonfault_face_manual, fault_adjacent, ndof_total, SXY, true);
   const real_t step1_nonfault_manual_sxz_lower_sig =
      MaxSortedSignatureDrift(nonfault_face_manual, fault_adjacent, ndof_total, SXZ, false);
   const real_t step1_nonfault_manual_sxz_upper_sig =
      MaxSortedSignatureDrift(nonfault_face_manual, fault_adjacent, ndof_total, SXZ, true);
   Vector nonfault_compare = nonfault_residual;
   nonfault_compare -= nonfault_face_manual;
   const real_t nonfault_compare_max_abs = nonfault_compare.Normlinf();

   TEST_LE(step0_sxy_lower_sig, 1.0e-10, "step-0 lower-side sorted SXY signatures match across all fault-adjacent tets");
   TEST_LE(step0_sxy_upper_sig, 1.0e-10, "step-0 upper-side sorted SXY signatures match across all fault-adjacent tets");
   TEST_LE(step0_sxz_lower_sig, 1.0e-10, "step-0 lower-side sorted SXZ signatures match across all fault-adjacent tets");
   TEST_LE(step0_sxz_upper_sig, 1.0e-10, "step-0 upper-side sorted SXZ signatures match across all fault-adjacent tets");
   TEST_LE(step1_fault_only_sxy_lower_sig, 1.0e-10, "second-step fault-face-only lower-side sorted SXY signatures match");
   TEST_LE(step1_fault_only_sxy_upper_sig, 1.0e-10, "second-step fault-face-only upper-side sorted SXY signatures match");
   TEST_LE(step1_fault_only_sxz_lower_sig, 1.0e-10, "second-step fault-face-only lower-side sorted SXZ signatures match");
   TEST_LE(step1_fault_only_sxz_upper_sig, 1.0e-10, "second-step fault-face-only upper-side sorted SXZ signatures match");
   TEST_LE(step1_volume_only_sxy_lower_sig, 1.0e-10, "second-step volume-only lower-side sorted SXY signatures match");
   TEST_LE(step1_volume_only_sxy_upper_sig, 1.0e-10, "second-step volume-only upper-side sorted SXY signatures match");
   TEST_LE(step1_volume_only_sxz_lower_sig, 1.0e-10, "second-step volume-only lower-side sorted SXZ signatures match");
   TEST_LE(step1_volume_only_sxz_upper_sig, 1.0e-10, "second-step volume-only upper-side sorted SXZ signatures match");
   TEST_LE(step1_boundary_sxy_lower_sig, 1.0e-10, "second-step boundary-face lower-side sorted SXY signatures match");
   TEST_LE(step1_boundary_sxy_upper_sig, 1.0e-10, "second-step boundary-face upper-side sorted SXY signatures match");
   TEST_LE(step1_boundary_sxz_lower_sig, 1.0e-10, "second-step boundary-face lower-side sorted SXZ signatures match");
   TEST_LE(step1_boundary_sxz_upper_sig, 1.0e-10, "second-step boundary-face upper-side sorted SXZ signatures match");
   TEST_LE(step1_interior_nf_sxy_lower_sig, 1.0e-10, "second-step interior-nonfault lower-side sorted SXY signatures match");
   TEST_LE(step1_interior_nf_sxy_upper_sig, 1.0e-10, "second-step interior-nonfault upper-side sorted SXY signatures match");
   TEST_LE(step1_interior_nf_sxz_lower_sig, 1.0e-10, "second-step interior-nonfault lower-side sorted SXZ signatures match");
   TEST_LE(step1_interior_nf_sxz_upper_sig, 1.0e-10, "second-step interior-nonfault upper-side sorted SXZ signatures match");
   TEST_LE(step1_nonfault_manual_sxy_lower_sig, 1.0e-10, "second-step manual nonfault-face lower-side sorted SXY signatures match");
   TEST_LE(step1_nonfault_manual_sxy_upper_sig, 1.0e-10, "second-step manual nonfault-face upper-side sorted SXY signatures match");
   TEST_LE(step1_nonfault_manual_sxz_lower_sig, 1.0e-10, "second-step manual nonfault-face lower-side sorted SXZ signatures match");
   TEST_LE(step1_nonfault_manual_sxz_upper_sig, 1.0e-10, "second-step manual nonfault-face upper-side sorted SXZ signatures match");
   TEST_LE(step1_nonfault_resid_sxy_lower_sig, 1.0e-10, "second-step nonfault residual lower-side sorted SXY signatures match");
   TEST_LE(step1_nonfault_resid_sxy_upper_sig, 1.0e-10, "second-step nonfault residual upper-side sorted SXY signatures match");
   TEST_LE(step1_nonfault_resid_sxz_lower_sig, 1.0e-10, "second-step nonfault residual lower-side sorted SXZ signatures match");
   TEST_LE(step1_nonfault_resid_sxz_upper_sig, 1.0e-10, "second-step nonfault residual upper-side sorted SXZ signatures match");
   TEST_LE(nonfault_compare_max_abs, 1.0e-11, "manual nonfault-face assembly matches nonfault residual exactly");
   TEST_LE(step1_sxy_lower_sig, 1.0e-10, "step-1 lower-side sorted SXY signatures match across all fault-adjacent tets");
   TEST_LE(step1_sxy_upper_sig, 1.0e-10, "step-1 upper-side sorted SXY signatures match across all fault-adjacent tets");
   TEST_LE(step1_sxz_lower_sig, 1.0e-10, "step-1 lower-side sorted SXZ signatures match across all fault-adjacent tets");
   TEST_LE(step1_sxz_upper_sig, 1.0e-10, "step-1 upper-side sorted SXZ signatures match across all fault-adjacent tets");

   // ---------------------------------------------------------------
   // Round-13A Gate 4b / 4c: fault-face mean decomposition.
   // ---------------------------------------------------------------
   // Reset dof_data to the post-step-0 state and re-run the fault
   // Riemann solve in three mean-override modes on the SAME input
   // Q_step0.  RunFaultFaceMeanAudit never calls WriteBackState, so
   // dof_data's friction outputs are not disturbed between calls.
   //
   // Gate 4b answers: "do mirrored fault faces already disagree at the
   // face-mean level?"  If yes, the remaining ~60% residual from Round
   // 12 is primarily between-face mean variation (hypothesis A).
   //
   // Gate 4c answers: "does forcing global fault uniformity close the
   // fault-only deposit drift?"  If face-mean doesn't close but
   // global-mean does, hypothesis A is confirmed.  If neither closes,
   // the residual is deposit/topology (hypothesis B).

   // Reset dof_data to the post-step-0 state so the mean audit sees
   // the same psi / friction inputs the second ADER step actually uses.
   ResetFaultState(dof_data, fault_coords);
   wave.SetFaultDOFData(&dof_data, wave.GetNbfPerFace());
   Vector Q_mean_in(wave.Height());
   Vector Q_step0_for_mean(wave.Height());
   Q_mean_in = 0.0;
   wave.AdvanceADER(Q_mean_in, kDt, kAderOrder, Q_step0_for_mean);
   // `dof_data` now holds the step-0 friction outputs; psi is the
   // step-0-exit psi (same as what RunADERTraceAudit above used).

   // Start with Tcorr source (per plan §9 minimal matrix).
   FaultFaceMeanAudit fault_mean_baseline =
      RunFaultFaceMeanAudit(wave, mesh, ff, dof_data, Q_step0_for_mean, kDt,
                            FaultMeanMode::None, FaultMeanSource::Tcorr);
   FaultFaceMeanAudit fault_mean_face =
      RunFaultFaceMeanAudit(wave, mesh, ff, dof_data, Q_step0_for_mean, kDt,
                            FaultMeanMode::Face, FaultMeanSource::Tcorr);
   FaultFaceMeanAudit fault_mean_global =
      RunFaultFaceMeanAudit(wave, mesh, ff, dof_data, Q_step0_for_mean, kDt,
                            FaultMeanMode::Global, FaultMeanSource::Tcorr);

   std::cout << "\n-- Gate 4b: fault-face mean orbit audit --\n";
   {
      const auto &faces = fault_mean_baseline.faces;
      std::cout << "    (n fault faces = " << faces.size() << ")\n";
      auto report = [&](const char *label,
                        std::function<real_t(const FaultFaceMeanSample&)> g)
      {
         std::string bucket;
         const real_t drift = MaxFaultFaceMeanOrbitDrift(faces, g, bucket);
         std::cout << "    " << label << " face-mean orbit drift = "
                   << std::scientific << std::setprecision(3) << drift
                   << "  " << (bucket.empty() ? std::string("(no mirror pair)") : bucket)
                   << "\n";
      };
      report("sigma_n_trial ", [](const FaultFaceMeanSample &f) { return f.sigma_n_trial_mean; });
      report("tau1_trial    ", [](const FaultFaceMeanSample &f) { return f.tau1_trial_mean; });
      report("tau2_trial    ", [](const FaultFaceMeanSample &f) { return f.tau2_trial_mean; });
      report("sigma_n_corr  ", [](const FaultFaceMeanSample &f) { return f.sigma_n_corr_mean; });
      report("tau1_corr     ", [](const FaultFaceMeanSample &f) { return f.tau1_corr_mean; });
      report("tau2_corr     ", [](const FaultFaceMeanSample &f) { return f.tau2_corr_mean; });
      report("F_h_plus[SXY] ", [](const FaultFaceMeanSample &f) { return f.fh_plus_sxy_mean; });
      report("F_h_plus[SXZ] ", [](const FaultFaceMeanSample &f) { return f.fh_plus_sxz_mean; });
      report("F_h_minus[SXY]", [](const FaultFaceMeanSample &f) { return f.fh_minus_sxy_mean; });
      report("F_h_minus[SXZ]", [](const FaultFaceMeanSample &f) { return f.fh_minus_sxz_mean; });
   }

   std::cout << "\n-- Gate 4c: fault-only deposit mean-mode comparison --\n";
   {
      auto drift_row = [&](const char *label, const FaultFaceMeanAudit &a)
      {
         const real_t sxy_l = MaxSortedSignatureDrift(a.rhs_after_minv,
                                 fault_adjacent, ndof_total, SXY, false);
         const real_t sxy_u = MaxSortedSignatureDrift(a.rhs_after_minv,
                                 fault_adjacent, ndof_total, SXY, true);
         const real_t sxz_l = MaxSortedSignatureDrift(a.rhs_after_minv,
                                 fault_adjacent, ndof_total, SXZ, false);
         const real_t sxz_u = MaxSortedSignatureDrift(a.rhs_after_minv,
                                 fault_adjacent, ndof_total, SXZ, true);
         std::cout << "    " << label << ":"
                   << " SXY lower=" << std::scientific << std::setprecision(3) << sxy_l
                   << " upper=" << sxy_u
                   << "  SXZ lower=" << sxz_l
                   << " upper=" << sxz_u << "\n";
         return std::array<real_t, 4>{sxy_l, sxy_u, sxz_l, sxz_u};
      };
      const auto d_base   = drift_row("baseline   ", fault_mean_baseline);
      const auto d_face   = drift_row("face-mean  ", fault_mean_face);
      const auto d_global = drift_row("global-mean", fault_mean_global);

      // Assertions — monotonicity only (per plan §6, Gate 4c): each
      // successive mode should not make drift worse.  Small slack tol
      // allows for accumulated roundoff in the sum-then-average path.
      const real_t mono_tol = 1e-10;
      for (int i = 0; i < 4; i++)
      {
         const char *comp = (i == 0) ? "SXY lower"
                         : (i == 1) ? "SXY upper"
                         : (i == 2) ? "SXZ lower"
                                    : "SXZ upper";
         std::ostringstream msg1;
         msg1 << "Gate 4c: face-mean drift <= baseline + tol (" << comp << ")";
         std::ostringstream msg2;
         msg2 << "Gate 4c: global-mean drift <= face-mean + tol (" << comp << ")";
         TEST_LE(d_face[i]   - d_base[i], mono_tol, msg1.str().c_str());
         TEST_LE(d_global[i] - d_face[i], mono_tol, msg2.str().c_str());
      }
   }

   // ---------------------------------------------------------------
   // Round-13B Gate 4d / 4e: second-step non-fault seed audit.
   // ---------------------------------------------------------------
   // Round 12 left ~60% of pepper unexplained; Round 13A ruled out
   // between-face fault-mean variation at step 1.  Round 13B probes
   // the non-fault face flux (interior + boundary) that generates
   // Q_step0's asymmetry via seven mode-controlled audit runs.
   // Gate 4d prints per-class orbit drift of raw F_h; Gate 4e prints
   // fault-adjacent sorted-signature drift for each mode and asserts
   // monotonicity only.

   NonFaultSeedAudit nf_baseline = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt, NonFaultSeedMode::Baseline);
   NonFaultSeedAudit nf_isym = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt, NonFaultSeedMode::InteriorSym);
   NonFaultSeedAudit nf_bsym = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt, NonFaultSeedMode::BoundarySym);
   NonFaultSeedAudit nf_bothsym = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt, NonFaultSeedMode::BothSym);
   NonFaultSeedAudit nf_imean = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt, NonFaultSeedMode::InteriorFaceMean);
   NonFaultSeedAudit nf_bmean = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt, NonFaultSeedMode::BoundaryFaceMean);
   NonFaultSeedAudit nf_bothmean = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt, NonFaultSeedMode::BothFaceMean);

   std::cout << "\n-- Gate 4d: second-step non-fault face orbit audit --\n";
   {
      std::cout << "    (n non-fault faces = " << nf_baseline.faces.size() << ")\n";
      auto drift_on = [&](bool interior_only, bool boundary_only,
                          std::function<real_t(const NonFaultFaceSample&)> g)
      {
         std::string bucket;
         const real_t d = MaxNonFaultFaceOrbitDrift(nf_baseline.faces, g,
                                                    boundary_only, interior_only, bucket);
         return std::make_pair(d, bucket);
      };
      auto rep_interior = [&](const char *label,
                              std::function<real_t(const NonFaultFaceSample&)> g)
      {
         const auto d = drift_on(/*interior_only=*/true, /*boundary_only=*/false, g);
         std::cout << "    interior raw F_h orbit drift: " << label
                   << " = " << std::scientific << std::setprecision(3) << d.first
                   << "  " << (d.second.empty() ? std::string("(no mirror pair)") : d.second)
                   << "\n";
      };
      auto rep_boundary = [&](const char *label,
                              std::function<real_t(const NonFaultFaceSample&)> g)
      {
         const auto d = drift_on(/*interior_only=*/false, /*boundary_only=*/true, g);
         std::cout << "    boundary raw F_h orbit drift: " << label
                   << " = " << std::scientific << std::setprecision(3) << d.first
                   << "  " << (d.second.empty() ? std::string("(no mirror pair)") : d.second)
                   << "\n";
      };
      rep_interior("SXY", [](const NonFaultFaceSample &s) { return s.fh_raw[SXY]; });
      rep_interior("SXZ", [](const NonFaultFaceSample &s) { return s.fh_raw[SXZ]; });
      rep_boundary("SXY", [](const NonFaultFaceSample &s) { return s.fh_raw[SXY]; });
      rep_boundary("SXZ", [](const NonFaultFaceSample &s) { return s.fh_raw[SXZ]; });
   }

   std::cout << "\n-- Gate 4e: second-step non-fault seed comparison --\n";
   {
      auto four_drifts = [&](const NonFaultSeedAudit &a)
      {
         return std::array<real_t, 4>{
            MaxSortedSignatureDrift(a.rhs_after_minv, fault_adjacent, ndof_total, SXY, false),
            MaxSortedSignatureDrift(a.rhs_after_minv, fault_adjacent, ndof_total, SXY, true),
            MaxSortedSignatureDrift(a.rhs_after_minv, fault_adjacent, ndof_total, SXZ, false),
            MaxSortedSignatureDrift(a.rhs_after_minv, fault_adjacent, ndof_total, SXZ, true)
         };
      };
      auto max4 = [](const std::array<real_t, 4> &d)
      {
         return std::max({d[0], d[1], d[2], d[3]});
      };
      auto row = [&](const char *label, const std::array<real_t, 4> &d)
      {
         std::cout << "    " << label << ":"
                   << " SXY L=" << std::scientific << std::setprecision(3) << d[0]
                   << " U=" << d[1]
                   << "  SXZ L=" << d[2]
                   << " U=" << d[3] << "\n";
      };

      const auto d_base    = four_drifts(nf_baseline);
      const auto d_isym    = four_drifts(nf_isym);
      const auto d_bsym    = four_drifts(nf_bsym);
      const auto d_both    = four_drifts(nf_bothsym);
      const auto d_imean   = four_drifts(nf_imean);
      const auto d_bmean   = four_drifts(nf_bmean);
      const auto d_bothm   = four_drifts(nf_bothmean);

      row("baseline    ", d_base);
      row("interiorSym ", d_isym);
      row("boundarySym ", d_bsym);
      row("bothSym     ", d_both);
      row("interiorMean", d_imean);
      row("boundaryMean", d_bmean);
      row("bothMean    ", d_bothm);

      // Assertions — monotonicity only (per plan §8).  Compare the
      // max-over-four-components scalar for each pair.  Small slack
      // (1e-12) allows for accumulated roundoff in the mode passes.
      const real_t mono_tol = 1e-12;
      const real_t m_base  = max4(d_base);
      const real_t m_isym  = max4(d_isym);
      const real_t m_bsym  = max4(d_bsym);
      const real_t m_both  = max4(d_both);
      const real_t m_imean = max4(d_imean);
      const real_t m_bmean = max4(d_bmean);
      const real_t m_bothm = max4(d_bothm);
      TEST_LE(m_isym  - m_base,  mono_tol, "Gate 4e: interiorSym   <= baseline + tol");
      TEST_LE(m_bsym  - m_base,  mono_tol, "Gate 4e: boundarySym   <= baseline + tol");
      TEST_LE(m_both  - m_isym,  mono_tol, "Gate 4e: bothSym       <= interiorSym + tol");
      TEST_LE(m_both  - m_bsym,  mono_tol, "Gate 4e: bothSym       <= boundarySym + tol");
      TEST_LE(m_imean - m_base,  mono_tol, "Gate 4e: interiorMean  <= baseline + tol");
      TEST_LE(m_bmean - m_base,  mono_tol, "Gate 4e: boundaryMean  <= baseline + tol");
      TEST_LE(m_bothm - m_imean, mono_tol, "Gate 4e: bothMean      <= interiorMean + tol");
      TEST_LE(m_bothm - m_bmean, mono_tol, "Gate 4e: bothMean      <= boundaryMean + tol");
   }

   // ---------------------------------------------------------------
   // Round-14A Gate 4g: x-side boundary isolation.
   // ---------------------------------------------------------------
   // Suppress the dominant interior seed with InteriorSym, then test
   // four x-boundary treatments to isolate the residual SXY floor.
   NonFaultSeedAudit nf_ixraw = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt,
      NonFaultSeedMode::InteriorSym_XBoundaryRaw);
   NonFaultSeedAudit nf_ixboth = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt,
      NonFaultSeedMode::InteriorSym_XBoundarySym);
   NonFaultSeedAudit nf_ixmin = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt,
      NonFaultSeedMode::InteriorSym_XMinSym);
   NonFaultSeedAudit nf_ixmax = RunSecondStepNonFaultSeedAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt,
      NonFaultSeedMode::InteriorSym_XMaxSym);

   std::cout << "\n-- Gate 4g: x-side boundary isolation (interior-sym baseline) --\n";
   {
      auto four_drifts = [&](const NonFaultSeedAudit &a)
      {
         return std::array<real_t, 4>{
            MaxSortedSignatureDrift(a.rhs_after_minv, fault_adjacent, ndof_total, SXY, false),
            MaxSortedSignatureDrift(a.rhs_after_minv, fault_adjacent, ndof_total, SXY, true),
            MaxSortedSignatureDrift(a.rhs_after_minv, fault_adjacent, ndof_total, SXZ, false),
            MaxSortedSignatureDrift(a.rhs_after_minv, fault_adjacent, ndof_total, SXZ, true)
         };
      };
      auto max4 = [](const std::array<real_t, 4> &d)
      {
         return std::max({d[0], d[1], d[2], d[3]});
      };
      auto row = [&](const char *label, const std::array<real_t, 4> &d)
      {
         std::cout << "    " << label << ":"
                   << " SXY L=" << std::scientific << std::setprecision(3) << d[0]
                   << " U=" << d[1]
                   << "  SXZ L=" << d[2]
                   << " U=" << d[3] << "\n";
      };

      const auto d_xraw  = four_drifts(nf_ixraw);
      const auto d_xboth = four_drifts(nf_ixboth);
      const auto d_xmin  = four_drifts(nf_ixmin);
      const auto d_xmax  = four_drifts(nf_ixmax);

      row("iSym+xBndRaw", d_xraw);
      row("iSym+xBndSym", d_xboth);
      row("iSym+xMinSym", d_xmin);
      row("iSym+xMaxSym", d_xmax);

      // Weak monotonicity only: each selective-sym variant should
      // not WORSEN drift relative to the interior-only-sym baseline.
      const real_t mono_tol = 1e-12;
      const real_t m_xraw  = max4(d_xraw);
      const real_t m_xboth = max4(d_xboth);
      const real_t m_xmin  = max4(d_xmin);
      const real_t m_xmax  = max4(d_xmax);
      TEST_LE(m_xboth - m_xraw, mono_tol, "Gate 4g: xBndSym <= xBndRaw + tol");
      TEST_LE(m_xmin  - m_xraw, mono_tol, "Gate 4g: xMinSym <= xBndRaw + tol");
      TEST_LE(m_xmax  - m_xraw, mono_tol, "Gate 4g: xMaxSym <= xBndRaw + tol");
   }

   // ---------------------------------------------------------------
   // Round-14B Gate 4h: x-side boundary SXY orbit drift.
   // ---------------------------------------------------------------
   // Reuse the Round-13B boundary-mode samples.  For each of
   // baseline / boundarySym / boundaryMean, print the x=0 and x=L
   // orbit drifts side-by-side, plus face counts.  Diagnostic only.
   std::cout << "\n-- Gate 4h: x-side boundary SXY orbit drift --\n";
   {
      auto n_on_side = [&](const std::vector<NonFaultFaceSample> &faces,
                           BoundarySide side)
      {
         int n = 0;
         for (const auto &f : faces)
         { if (f.is_boundary && f.side == side) { ++n; } }
         return n;
      };
      const int n_x0 = n_on_side(nf_baseline.faces, BoundarySide::XMin);
      const int n_xL = n_on_side(nf_baseline.faces, BoundarySide::XMax);
      std::cout << "    (n x=0 faces = " << n_x0
                << ", n x=L faces = " << n_xL << ")\n";

      auto getter_raw = [](const NonFaultFaceSample &s) { return s.fh_raw[SXY]; };
      auto getter_sym = [](const NonFaultFaceSample &s) { return s.fh_sym[SXY]; };

      auto print_row = [&](const char *label,
                           std::function<real_t(const NonFaultFaceSample&)> g)
      {
         std::string l_x0, l_xL;
         const real_t d_x0 = MaxBoundarySideOrbitDrift(
            nf_baseline.faces, BoundarySide::XMin, g, l_x0);
         const real_t d_xL = MaxBoundarySideOrbitDrift(
            nf_baseline.faces, BoundarySide::XMax, g, l_xL);
         std::cout << "    " << label << ":"
                   << " x=0=" << std::scientific << std::setprecision(3) << d_x0
                   << " [" << l_x0 << "]"
                   << "   x=L=" << d_xL
                   << " [" << l_xL << "]\n";
      };
      print_row("baseline    ", getter_raw);
      print_row("boundarySym ", getter_sym);
      print_row("boundaryMean", getter_raw);
   }

   // ---------------------------------------------------------------
   // Round-15 Gate 4i: x-side boundary trace comparison.
   // ---------------------------------------------------------------
   XSideBoundaryTraceAudit xtrace = RunXSideBoundaryTraceAudit(
      wave, mesh, bc, Q_step0_for_mean, kDt);

   std::cout << "\n-- Gate 4i: x-side boundary trace comparison --\n";
   std::cout << "    (n x_faces samples = " << xtrace.x_faces.size()
             << ", n z_faces samples = " << xtrace.z_faces.size() << ")\n";
   {
      auto get_Iself_sxy = [](const XSideBoundaryTraceSample &s)
      { return s.I_self[SXY]; };
      auto get_bg_sxy = [](const XSideBoundaryTraceSample &s)
      { return s.bulk_bg_scaled[SXY]; };
      auto get_Fh_sxy = [](const XSideBoundaryTraceSample &s)
      { return s.F_h[SXY]; };

      std::vector<XSideBoundaryTraceSample> xbucket_Iself, xbucket_Fh;
      std::string xl_I, xl_bg, xl_F;
      const real_t x_I  = MaxTraceOrbitDrift(xtrace.x_faces, get_Iself_sxy, xl_I, &xbucket_Iself);
      const real_t x_bg = MaxTraceOrbitDrift(xtrace.x_faces, get_bg_sxy,    xl_bg);
      const real_t x_F  = MaxTraceOrbitDrift(xtrace.x_faces, get_Fh_sxy,    xl_F,  &xbucket_Fh);

      std::vector<XSideBoundaryTraceSample> zbucket_Fh;
      std::string zl_I, zl_bg, zl_F;
      const real_t z_I  = MaxTraceOrbitDrift(xtrace.z_faces, get_Iself_sxy, zl_I);
      const real_t z_bg = MaxTraceOrbitDrift(xtrace.z_faces, get_bg_sxy,    zl_bg);
      const real_t z_F  = MaxTraceOrbitDrift(xtrace.z_faces, get_Fh_sxy,    zl_F, &zbucket_Fh);

      std::cout << "    x-side worst drift: I_self[SXY]="
                << std::scientific << std::setprecision(3) << x_I
                << " bulk_bg[SXY]=" << x_bg
                << " F_h[SXY]=" << x_F << "\n";
      std::cout << "    z-side worst drift: I_self[SXY]=" << z_I
                << " bulk_bg[SXY]=" << z_bg
                << " F_h[SXY]=" << z_F << "\n";

      if (x_F > 1e-13)
      {
         std::cout << "    worst x-side F_h[SXY] bucket " << xl_F << ":\n";
         PrintTraceBucket(xbucket_Fh, "x F_h[SXY] worst");
      }
      if (x_I > 1e-13 && !xbucket_Iself.empty())
      {
         std::cout << "    worst x-side I_self[SXY] bucket " << xl_I << ":\n";
         PrintTraceBucket(xbucket_Iself, "x I_self[SXY] worst");
      }
      if (z_F > 1e-13 && !zbucket_Fh.empty())
      {
         std::cout << "    worst z-side F_h[SXY] bucket " << zl_F << ":\n";
         PrintTraceBucket(zbucket_Fh, "z F_h[SXY] worst");
      }

      // Gate 4j: dispatch-consistency within x-side orbit buckets.
      // Scan each x orbit bucket and verify all samples share the
      // same bdr_attr and bc_type.  If they differ, that is already
      // the root cause class.
      std::cout << "\n-- Gate 4j: x-side dispatch consistency --\n";
      std::map<XSideTraceOrbitKey,
               std::vector<XSideBoundaryTraceSample>> buckets;
      for (const auto &s : xtrace.x_faces)
      {
         buckets[MakeBoundaryTraceOrbitKey(s)].push_back(s);
      }
      int n_attr_mismatch = 0, n_bc_mismatch = 0, n_checked = 0;
      for (const auto &kv : buckets)
      {
         if (kv.second.size() < 2) { continue; }
         ++n_checked;
         const int attr0 = kv.second.front().bdr_attr;
         const AuditFaceBC bc0 = kv.second.front().bc_type;
         bool attr_ok = true, bc_ok = true;
         for (const auto &s : kv.second)
         {
            if (s.bdr_attr != attr0) { attr_ok = false; }
            if (s.bc_type  != bc0)   { bc_ok   = false; }
         }
         if (!attr_ok)
         {
            ++n_attr_mismatch;
            std::cout << "    attr mismatch in bucket(side="
                      << std::get<0>(kv.first) << ", attr(first)="
                      << attr0 << ")\n";
         }
         if (!bc_ok)
         {
            ++n_bc_mismatch;
            std::cout << "    bc_type mismatch in bucket(side="
                      << std::get<0>(kv.first) << ", bc(first)="
                      << FaceBCLabel(bc0) << ")\n";
         }
      }
      std::cout << "    checked " << n_checked
                << " x-side orbit buckets; attr mismatches=" << n_attr_mismatch
                << ", bc_type mismatches=" << n_bc_mismatch << "\n";
      TEST_LE(static_cast<real_t>(n_attr_mismatch), 0.0,
              "Gate 4j: x-side orbit mates share bdr_attr");
      TEST_LE(static_cast<real_t>(n_bc_mismatch), 0.0,
              "Gate 4j: x-side orbit mates share bc_type");
   }

   // Round-15 step 9: forced-kernel follow-up.  Gate 4i showed
   // I_self[SXY] is non-zero on BOTH x-sides and z-sides (~0.09 Pa),
   // but F_h[SXY] amplifies to 323 Pa on x-sides and 0 on z-sides.
   // Conclusion: the BC formula is the side-discriminating amplifier,
   // not the trace reconstruction.  Run the same audit forcing each
   // kernel to see which one(s) amplify the trace asymmetry on
   // x-sides.  Diagnostic only.
   std::cout << "\n-- Gate 4i-forced: x-side F_h drift under forced BC kernels --\n";
   {
      const BoundaryKernelOverride overrides[3] = {
         BoundaryKernelOverride::ForceAbsorbing,
         BoundaryKernelOverride::ForceFreeSurface,
         BoundaryKernelOverride::ForceFreeSurfaceGodunov
      };
      for (auto o : overrides)
      {
         XSideBoundaryTraceAudit xa = RunXSideBoundaryTraceAudit(
            wave, mesh, bc, Q_step0_for_mean, kDt, o);
         std::string xl_F, zl_F;
         const real_t x_F = MaxTraceOrbitDrift(xa.x_faces,
            [](const XSideBoundaryTraceSample &s) { return s.F_h[SXY]; }, xl_F);
         const real_t z_F = MaxTraceOrbitDrift(xa.z_faces,
            [](const XSideBoundaryTraceSample &s) { return s.F_h[SXY]; }, zl_F);
         std::cout << "    " << BoundaryKernelOverrideLabel(o)
                   << ": x F_h[SXY]="
                   << std::scientific << std::setprecision(3) << x_F
                   << "  z F_h[SXY]=" << z_F << "\n";
      }
   }

   // ---------------------------------------------------------------
   // Round-16 Gate 4k / 4l: I_data vs trace-sampling split.
   // ---------------------------------------------------------------
   // Round 15 established the BC formulas are clean amplifiers.  The
   // remaining question is whether the I_self[SXY] asymmetry arises
   // from (a) production-style contiguous `dof_offset + i` addressing
   // diverging from explicit `edofs[i]` lookups, or (b) clean trace
   // sampling applied to an I_data that is itself already dirty.
   BoundaryTraceSamplerAudit sampler_audit =
      RunBoundaryTraceSamplerComparisonAudit(wave, mesh, Q_step0_for_mean, kDt);
   ElementIDataAudit elem_i_audit =
      RunElementIDataAudit(wave, mesh, Q_step0_for_mean, kDt);

   std::cout << "\n-- Gate 4k: boundary sampler comparison --\n";
   std::cout << "    (n x samples = " << sampler_audit.x_faces.size()
             << ", n z samples = " << sampler_audit.z_faces.size() << ")\n";
   {
      auto get_offset = [](const BoundaryTraceSamplerSample &s)
      { return s.i_self_offset_sxy; };
      auto get_edofs  = [](const BoundaryTraceSamplerSample &s)
      { return s.i_self_edofs_sxy; };
      auto get_delta  = [](const BoundaryTraceSamplerSample &s)
      { return s.i_self_offset_sxy - s.i_self_edofs_sxy; };

      std::string lxo, lxe, lxd, lzo, lze, lzd;
      std::vector<BoundaryTraceSamplerSample> xb_delta, xb_offset, zb_delta;
      const real_t x_off   = MaxSamplerOrbitDrift(sampler_audit.x_faces, get_offset, lxo, &xb_offset);
      const real_t x_ed    = MaxSamplerOrbitDrift(sampler_audit.x_faces, get_edofs,  lxe);
      const real_t x_delta = MaxSamplerOrbitDrift(sampler_audit.x_faces, get_delta,  lxd, &xb_delta);
      const real_t z_off   = MaxSamplerOrbitDrift(sampler_audit.z_faces, get_offset, lzo);
      const real_t z_ed    = MaxSamplerOrbitDrift(sampler_audit.z_faces, get_edofs,  lze);
      const real_t z_delta = MaxSamplerOrbitDrift(sampler_audit.z_faces, get_delta,  lzd, &zb_delta);

      std::cout << "    x-side: offset=" << std::scientific << std::setprecision(3) << x_off
                << " edofs=" << x_ed
                << " delta(offset-edofs)=" << x_delta << "\n";
      std::cout << "    z-side: offset=" << z_off
                << " edofs=" << z_ed
                << " delta(offset-edofs)=" << z_delta << "\n";

      if (x_delta > 1e-13)
      {
         std::cout << "    worst x-side DELTA bucket " << lxd << ":\n";
         PrintBoundarySamplerBucket(xb_delta, "x DELTA worst");
      }
      if (z_delta > 1e-13)
      {
         std::cout << "    worst z-side DELTA bucket " << lzd << ":\n";
         PrintBoundarySamplerBucket(zb_delta, "z DELTA worst");
      }
      if (x_off > 1e-13 && !xb_offset.empty())
      {
         std::cout << "    worst x-side offset bucket " << lxo << ":\n";
         PrintBoundarySamplerBucket(xb_offset, "x offset worst");
      }

      // Weak assertion: delta must not exceed offset by more than the
      // roundoff tolerance (i.e. contiguous addressing can at worst
      // equal the full signal; it should never introduce drift that
      // exceeds the physical signal amplitude).  Hard enforcement of
      // offset == edofs is deferred — report both magnitudes.
      TEST_LE(x_delta - x_off, 1e-10, "Gate 4k: x delta <= offset + tol");
      TEST_LE(z_delta - z_off, 1e-10, "Gate 4k: z delta <= offset + tol");
   }

   std::cout << "\n-- Gate 4l: element I_data[SXY] orbit audit --\n";
   {
      auto get_mean = [](const ElementIDataSample &s) { return s.i_sxy_mean; };
      std::string lm;
      std::vector<ElementIDataSample> worst_bucket;
      const real_t d_mean = MaxElementOrbitDrift(elem_i_audit.elems,
                                                  get_mean, lm, &worst_bucket);
      std::cout << "    (n elements = " << elem_i_audit.elems.size() << ")\n";
      std::cout << "    element I_data[SXY] mean orbit drift = "
                << std::scientific << std::setprecision(3) << d_mean << "\n";
      if (d_mean > 1e-13 && !worst_bucket.empty())
      {
         std::cout << "    worst element-orbit bucket " << lm << ":\n";
         PrintElementIDataBucket(worst_bucket, "element SXY worst");
      }
   }

   // ---------------------------------------------------------------
   // Round-17 Gate 4m / 4n: canonical per-DOF + face-interp audit.
   // ---------------------------------------------------------------
   // Round 16 showed element MEANS are orbit-clean but per-QP face
   // traces are not.  Round 17 distinguishes two remaining causes:
   // (i) local DOF ordering is non-equivariant but the data itself is
   //     permutation-equivalent → canonical sort closes drift;
   // (ii) higher modes of I_data already differ → canonical sort
   //     still dirty.
   ElementSxyDofAudit elem_dof_audit =
      RunElementSxyDofAudit(wave, mesh, Q_step0_for_mean, kDt);
   FaceInterpolationAudit face_interp_audit =
      RunFaceInterpolationAudit(wave, mesh, Q_step0_for_mean, kDt);

   std::cout << "\n-- Gate 4m: element SXY per-DOF orbit audit --\n";
   {
      std::string lr, lc;
      std::vector<ElementSxyDofSample> bucket_raw, bucket_canon;
      const real_t d_raw = MaxElementRawDofOrbitDrift(elem_dof_audit.elems,
                                                       lr, &bucket_raw);
      const real_t d_can = MaxElementCanonicalDofOrbitDrift(elem_dof_audit.elems,
                                                             lc, &bucket_canon);
      std::cout << "    (n elements = " << elem_dof_audit.elems.size() << ")\n";
      std::cout << "    raw-order drift       = "
                << std::scientific << std::setprecision(3) << d_raw << "\n";
      std::cout << "    canonical-order drift = " << d_can << "\n";
      if (d_raw > 1e-13 && !bucket_raw.empty())
      {
         std::cout << "    worst raw-order bucket " << lr << ":\n";
         PrintElementSxyDofBucket(bucket_raw, "element SXY worst raw",
                                  /*canonical=*/false);
      }
      if (d_can > 1e-13 && !bucket_canon.empty())
      {
         std::cout << "    worst canonical-order bucket " << lc << ":\n";
         PrintElementSxyDofBucket(bucket_canon, "element SXY worst canonical",
                                  /*canonical=*/true);
      }
      // Weak monotonicity — canonical should never exceed raw by
      // more than roundoff (sorting cannot fabricate drift).
      TEST_LE(d_can - d_raw, 1e-10,
              "Gate 4m: canonical-order drift <= raw-order drift + tol");
   }

   std::cout << "\n-- Gate 4n: face interpolation canonicalization audit --\n";
   {
      auto get_raw = [](const FaceInterpolationSample &s) { return s.interp_raw; };
      auto get_can = [](const FaceInterpolationSample &s) { return s.interp_canonical; };

      std::string lxr, lxc, lzr, lzc;
      std::vector<FaceInterpolationSample> xbucket_raw, xbucket_can;
      const real_t x_raw = MaxFaceInterpolationOrbitDrift(face_interp_audit.x_faces,
                                                           get_raw, lxr, &xbucket_raw);
      const real_t x_can = MaxFaceInterpolationOrbitDrift(face_interp_audit.x_faces,
                                                           get_can, lxc, &xbucket_can);
      const real_t z_raw = MaxFaceInterpolationOrbitDrift(face_interp_audit.z_faces,
                                                           get_raw, lzr);
      const real_t z_can = MaxFaceInterpolationOrbitDrift(face_interp_audit.z_faces,
                                                           get_can, lzc);
      std::cout << "    (n x samples = " << face_interp_audit.x_faces.size()
                << ", n z samples = " << face_interp_audit.z_faces.size() << ")\n";
      std::cout << "    x-side raw interp drift       = "
                << std::scientific << std::setprecision(3) << x_raw << "\n";
      std::cout << "    x-side canonical interp drift = " << x_can << "\n";
      std::cout << "    z-side raw interp drift       = " << z_raw << "\n";
      std::cout << "    z-side canonical interp drift = " << z_can << "\n";
      if (x_raw > 1e-13 && !xbucket_raw.empty())
      {
         std::cout << "    worst x-side raw bucket " << lxr << ":\n";
         PrintFaceInterpolationBucket(xbucket_raw, "x raw worst",
                                      /*canonical=*/false);
      }
      if (x_can > 1e-13 && !xbucket_can.empty())
      {
         std::cout << "    worst x-side canonical bucket " << lxc << ":\n";
         PrintFaceInterpolationBucket(xbucket_can, "x canonical worst",
                                      /*canonical=*/true);
      }
      TEST_LE(x_can - x_raw, 1e-10,
              "Gate 4n: x canonical interp drift <= raw interp drift + tol");
      TEST_LE(z_can - z_raw, 1e-10,
              "Gate 4n: z canonical interp drift <= raw interp drift + tol");
   }

   // ---------------------------------------------------------------
   // Round-18A Gate 4o: face-local canonical interpolation probe.
   // ---------------------------------------------------------------
   // Round-17 ruled out (a) local DOF ordering-only (Gate 4m raw !=
   // canonical) and (b) a straightforward element-canonical closure
   // (Gate 4n raw == canonical).  This gate adds a third ordering:
   // face-local canonical, where DOFs are sorted by their reference
   // coordinates restricted to the sampled face plane.  If this
   // closes the x/z interp drift, the bug is specifically in
   // boundary face sampling's sensitivity to element-local DOF
   // labeling — the fix target is a topology-aware face sampling
   // path (AddBoundaryFaceRhsFull-style) rather than mesh renumbering.
   FaceLocalCanonicalAudit flc_audit =
      RunFaceLocalCanonicalInterpolationAudit(wave, mesh,
                                               Q_step0_for_mean, kDt);

   std::cout << "\n-- Gate 4o: face-local canonical interpolation --\n";
   {
      auto get_raw  = [](const FaceLocalCanonicalSample &s) { return s.interp_raw; };
      auto get_eCan = [](const FaceLocalCanonicalSample &s) { return s.interp_elem_canonical; };
      auto get_fCan = [](const FaceLocalCanonicalSample &s) { return s.interp_face_canonical; };

      std::string lxr, lxe, lxf, lzr, lze, lzf;
      std::vector<FaceLocalCanonicalSample> xbf;
      const real_t x_raw  = MaxFaceLocalCanonicalOrbitDrift(flc_audit.x_faces, get_raw,  lxr);
      const real_t x_eCan = MaxFaceLocalCanonicalOrbitDrift(flc_audit.x_faces, get_eCan, lxe);
      const real_t x_fCan = MaxFaceLocalCanonicalOrbitDrift(flc_audit.x_faces, get_fCan, lxf, &xbf);
      const real_t z_raw  = MaxFaceLocalCanonicalOrbitDrift(flc_audit.z_faces, get_raw,  lzr);
      const real_t z_eCan = MaxFaceLocalCanonicalOrbitDrift(flc_audit.z_faces, get_eCan, lze);
      const real_t z_fCan = MaxFaceLocalCanonicalOrbitDrift(flc_audit.z_faces, get_fCan, lzf);

      std::cout << "    (n x samples = " << flc_audit.x_faces.size()
                << ", n z samples = " << flc_audit.z_faces.size() << ")\n";
      std::cout << "    x-side raw=" << std::scientific << std::setprecision(3) << x_raw
                << " elem-canonical=" << x_eCan
                << " face-canonical=" << x_fCan << "\n";
      std::cout << "    z-side raw=" << z_raw
                << " elem-canonical=" << z_eCan
                << " face-canonical=" << z_fCan << "\n";

      if (x_fCan > 1e-13 && !xbf.empty())
      {
         std::cout << "    worst x-side face-canonical bucket " << lxf << ":\n";
         PrintFaceLocalCanonicalBucket(xbf, "x face-canonical worst",
                                       "face-canonical");
      }

      // Weak monotonicity: face-canonical must not exceed raw by
      // more than roundoff (face-local sort cannot fabricate drift).
      TEST_LE(x_fCan - x_raw, 1e-10,
              "Gate 4o: x face-canonical drift <= raw drift + tol");
      TEST_LE(z_fCan - z_raw, 1e-10,
              "Gate 4o: z face-canonical drift <= raw drift + tol");
   }

   // ---------------------------------------------------------------
   // Round-13C Gate 4f: boundary-side SXY decomposition.
   // ---------------------------------------------------------------
   // Round-13B showed SXY has a residual ~1.0 Pa floor that neither
   // boundarySym nor boundaryMean can fully eliminate.  This gate
   // prints the boundary-side orbit drift of the face-F_h SXY
   // component per cube side (x=0, x=L, y=0, y=L, z=0, z=L) for
   // three modes:
   //   baseline     — raw F_h per boundary face
   //   boundarySym  — 0.5·(F(nor)+F(-nor)) per boundary face
   //   boundaryMean — raw F_h (plan §3 getter: fh_raw); the
   //                  orbit-mean transformation zeros within-bucket
   //                  drift by construction, so this row is
   //                  identical to baseline and serves as the plan's
   //                  parallel-format template.
   //
   // Diagnostic only — no assertions.  The intent is to localize the
   // remaining SXY floor to one specific cube side family.
   std::cout << "\n-- Gate 4f: boundary-side SXY decomposition --\n";
   {
      const BoundarySide sides[6] = {
         BoundarySide::XMin, BoundarySide::XMax,
         BoundarySide::YMin, BoundarySide::YMax,
         BoundarySide::ZMin, BoundarySide::ZMax
      };
      auto getter_raw = [](const NonFaultFaceSample &s) { return s.fh_raw[SXY]; };
      auto getter_sym = [](const NonFaultFaceSample &s) { return s.fh_sym[SXY]; };

      auto print_row = [&](const char *label,
                           std::function<real_t(const NonFaultFaceSample&)> g)
      {
         std::cout << "    " << label << ":";
         for (int i = 0; i < 6; i++)
         {
            std::string bucket;
            const real_t d = MaxBoundarySideOrbitDrift(nf_baseline.faces,
                                                       sides[i], g, bucket);
            std::cout << "  " << BoundarySideLabel(sides[i]) << "="
                      << std::scientific << std::setprecision(3) << d;
         }
         std::cout << "\n";
      };
      print_row("baseline    ", getter_raw);
      print_row("boundarySym ", getter_sym);
      print_row("boundaryMean", getter_raw);
   }

   std::cout << "\n-- Gate 5: per-face orbit drift in interior non-fault branch --\n";
   {
      int worst_c_raw = -1, worst_c_sym = -1, worst_c_can = -1;
      std::string orbit_raw, orbit_sym, orbit_can;
      const auto raw = CollectInteriorFaceSamples(wave, mesh, I_step1, kDt,
                                                  IfaceProbeMode::Raw);
      const real_t drift_raw = MaxOrbitFhDrift(raw, worst_c_raw, orbit_raw);
      std::cout << "  5a raw Interior(nor,L,R): max orbit F_h drift = "
                << std::scientific << std::setprecision(3) << drift_raw
                << " on " << (worst_c_raw >= 0 ? kCompName[worst_c_raw] : "-")
                << "\n";
      if (drift_raw > 1e-10) { ReportOrbitBreakdown("iface raw", raw); }

      const auto sym = CollectInteriorFaceSamples(wave, mesh, I_step1, kDt,
                                                  IfaceProbeMode::Symmetrized);
      const real_t drift_sym = MaxOrbitFhDrift(sym, worst_c_sym, orbit_sym);
      std::cout << "  5b symmetrized 0.5*(Interior(n,L,R)+Interior(-n,R,L)):"
                << " max orbit F_h drift = "
                << std::scientific << std::setprecision(3) << drift_sym
                << " on " << (worst_c_sym >= 0 ? kCompName[worst_c_sym] : "-")
                << "\n";
      if (drift_sym > 1e-10) { ReportOrbitBreakdown("iface sym", sym); }

      const auto can = CollectInteriorFaceSamples(wave, mesh, I_step1, kDt,
                                                  IfaceProbeMode::CanonicalNormal);
      const real_t drift_can = MaxOrbitFhDrift(can, worst_c_can, orbit_can);
      std::cout << "  5c canonical-normal Interior(nor_can,L,R): max orbit "
                << "F_h drift = "
                << std::scientific << std::setprecision(3) << drift_can
                << " on " << (worst_c_can >= 0 ? kCompName[worst_c_can] : "-")
                << "\n";
      if (drift_can > 1e-10) { ReportOrbitBreakdown("iface can", can); }

      TEST_LE(drift_raw, 1.0e-10,
              "Gate 5a: per-face F_h orbit-uniform on interior non-fault (raw Interior)");
      TEST_LE(drift_sym, 1.0e-10,
              "Gate 5b: per-face F_h orbit-uniform when Interior symmetrized in (nor,L<->R)");
      TEST_LE(drift_can, 1.0e-10,
              "Gate 5c: per-face F_h orbit-uniform when normal is canonicalized outward");
   }

   std::cout << "\n-- Gate 6: per-face orbit drift in boundary branch --\n";
   {
      int worst_c_gam = -1, worst_c_god = -1, worst_c_gsy = -1;
      std::string orbit_gam, orbit_god, orbit_gsy;
      const auto gam = CollectBoundaryFaceSamples(wave, mesh, bc, I_step1, kDt,
                                                  BFaceProbeMode::GammaBaseline);
      const real_t drift_gam = MaxOrbitFhDrift(gam, worst_c_gam, orbit_gam);
      std::cout << "  6a Gamma baseline FreeSurfaceTotal: max orbit F_h drift = "
                << std::scientific << std::setprecision(3) << drift_gam
                << " on " << (worst_c_gam >= 0 ? kCompName[worst_c_gam] : "-")
                << "\n";
      if (drift_gam > 1e-10) { ReportOrbitBreakdown("bface gamma", gam); }

      const auto god = CollectBoundaryFaceSamples(wave, mesh, bc, I_step1, kDt,
                                                  BFaceProbeMode::GodunovBaseline);
      const real_t drift_god = MaxOrbitFhDrift(god, worst_c_god, orbit_god);
      std::cout << "  6b Godunov variant FreeSurfaceGodunovTotal: max orbit "
                << "F_h drift = "
                << std::scientific << std::setprecision(3) << drift_god
                << " on " << (worst_c_god >= 0 ? kCompName[worst_c_god] : "-")
                << "\n";
      if (drift_god > 1e-10) { ReportOrbitBreakdown("bface godunov", god); }

      const auto gsy = CollectBoundaryFaceSamples(wave, mesh, bc, I_step1, kDt,
                                                  BFaceProbeMode::GammaSymmetrized);
      const real_t drift_gsy = MaxOrbitFhDrift(gsy, worst_c_gsy, orbit_gsy);
      std::cout << "  6c Gamma symmetrized 0.5*(Gamma(n)+Gamma(-n)): max orbit "
                << "F_h drift = "
                << std::scientific << std::setprecision(3) << drift_gsy
                << " on " << (worst_c_gsy >= 0 ? kCompName[worst_c_gsy] : "-")
                << "\n";
      if (drift_gsy > 1e-10) { ReportOrbitBreakdown("bface gsym", gsy); }

      TEST_LE(drift_gam, 1.0e-10,
              "Gate 6a: per-face F_h orbit-uniform on boundary (Gamma baseline)");
      TEST_LE(drift_god, 1.0e-10,
              "Gate 6b: per-face F_h orbit-uniform on boundary (Godunov variant)");
      TEST_LE(drift_gsy, 1.0e-10,
              "Gate 6c: per-face F_h orbit-uniform when boundary flux sign-symmetrized");
   }

   std::cout << "\n-- Gate 7: element-level lifted orbit drift after branch fixes --\n";
   {
      auto drift_of = [&](const NonFaultFaceAudit &a) {
         real_t worst = 0.0;
         for (int comp : {SXY, SXZ})
         {
            for (bool upper : {false, true})
            {
               worst = std::max(worst, MaxSortedSignatureDrift(
                  a.interior_after_minv, fault_adjacent, ndof_total, comp, upper));
               worst = std::max(worst, MaxSortedSignatureDrift(
                  a.boundary_after_minv, fault_adjacent, ndof_total, comp, upper));
            }
         }
         return worst;
      };
      const auto raw = RunADERNonFaultFaceAuditModed(
         wave, mesh, bc, I_step1, kDt, IfaceLiftMode::Raw, BFaceLiftMode::GammaRaw);
      const auto isym = RunADERNonFaultFaceAuditModed(
         wave, mesh, bc, I_step1, kDt, IfaceLiftMode::Symmetrized,
         BFaceLiftMode::GammaRaw);
      const auto bsym = RunADERNonFaultFaceAuditModed(
         wave, mesh, bc, I_step1, kDt, IfaceLiftMode::Raw,
         BFaceLiftMode::GammaSymmetrized);
      const auto both = RunADERNonFaultFaceAuditModed(
         wave, mesh, bc, I_step1, kDt, IfaceLiftMode::Symmetrized,
         BFaceLiftMode::GammaSymmetrized);

      auto report = [&](const char *label, const NonFaultFaceAudit &a) {
         std::cout << "    " << label
                   << "  iface SXY lower="
                   << std::scientific << std::setprecision(3)
                   << MaxSortedSignatureDrift(a.interior_after_minv,
                                              fault_adjacent, ndof_total, SXY, false)
                   << "  iface SXZ lower="
                   << MaxSortedSignatureDrift(a.interior_after_minv,
                                              fault_adjacent, ndof_total, SXZ, false)
                   << "  bface SXY lower="
                   << MaxSortedSignatureDrift(a.boundary_after_minv,
                                              fault_adjacent, ndof_total, SXY, false)
                   << "  bface SXZ lower="
                   << MaxSortedSignatureDrift(a.boundary_after_minv,
                                              fault_adjacent, ndof_total, SXZ, false)
                   << "\n";
      };
      report("[raw                ]", raw);
      report("[iface symmetrized  ]", isym);
      report("[bface symmetrized  ]", bsym);
      report("[both symmetrized   ]", both);

      const real_t drift_raw_total  = drift_of(raw);
      const real_t drift_isym_total = drift_of(isym);
      const real_t drift_bsym_total = drift_of(bsym);
      const real_t drift_both_total = drift_of(both);
      std::cout << "  totals:"
                << "  raw="         << std::scientific << std::setprecision(3)
                << drift_raw_total
                << "  iface_sym="   << drift_isym_total
                << "  bface_sym="   << drift_bsym_total
                << "  both_sym="    << drift_both_total << "\n";

      TEST_LE(drift_raw_total,  1.0e+01, "Gate 7 raw baseline reproduces the failure (sanity)");
      TEST_LE(drift_isym_total, 1.0e-10, "Gate 7a: interior-only symmetrization zeroes the lifted drift");
      TEST_LE(drift_bsym_total, 1.0e-10, "Gate 7b: boundary-only symmetrization zeroes the lifted drift");
      TEST_LE(drift_both_total, 1.0e-10, "Gate 7c: both symmetrizations together zero the lifted drift");
   }

   // Phase 3 (§8.2): Gate 7′ — duplicate Gate 7's drift measurement
   // but dispatched through PrecomputedFaceFluxes::Add*FaceRhs (the
   // actual §6.2 production dispatch when use_precomputed_face_fluxes_
   // is on).  Runs only when SEAS_TEST_USE_PRECOMPUTED_FLUX=1.
   // Acceptance (plan §8.3):
   //   P3.3 — iface SXY and SXZ orbit drift both ≤ 1e-10
   //   P3.4 — bface SXY and SXZ orbit drift both ≤ 1e-10 (critical)
   if (use_precomputed_env)
   {
      std::cout << "\n-- Gate 7′ (§8.3 P3.3/P3.4): precomputed-path "
                << "lifted drift on I_step1 --\n";
      const auto prec = RunPrecomputedFluxLiftedAudit(wave, mesh, I_step1, kDt);
      real_t iface_sxy_lower = MaxSortedSignatureDrift(
         prec.interior_after_minv, fault_adjacent, ndof_total, SXY, false);
      real_t iface_sxy_upper = MaxSortedSignatureDrift(
         prec.interior_after_minv, fault_adjacent, ndof_total, SXY, true);
      real_t iface_sxz_lower = MaxSortedSignatureDrift(
         prec.interior_after_minv, fault_adjacent, ndof_total, SXZ, false);
      real_t iface_sxz_upper = MaxSortedSignatureDrift(
         prec.interior_after_minv, fault_adjacent, ndof_total, SXZ, true);
      real_t bface_sxy_lower = MaxSortedSignatureDrift(
         prec.boundary_after_minv, fault_adjacent, ndof_total, SXY, false);
      real_t bface_sxy_upper = MaxSortedSignatureDrift(
         prec.boundary_after_minv, fault_adjacent, ndof_total, SXY, true);
      real_t bface_sxz_lower = MaxSortedSignatureDrift(
         prec.boundary_after_minv, fault_adjacent, ndof_total, SXZ, false);
      real_t bface_sxz_upper = MaxSortedSignatureDrift(
         prec.boundary_after_minv, fault_adjacent, ndof_total, SXZ, true);
      std::cout << "    iface SXY lower=" << std::scientific
                << std::setprecision(3) << iface_sxy_lower
                << "  upper=" << iface_sxy_upper
                << "  iface SXZ lower=" << iface_sxz_lower
                << "  upper=" << iface_sxz_upper << "\n";
      std::cout << "    bface SXY lower=" << bface_sxy_lower
                << "  upper=" << bface_sxy_upper
                << "  bface SXZ lower=" << bface_sxz_lower
                << "  upper=" << bface_sxz_upper << "\n";
      TEST_LE(iface_sxy_lower, 1.0e-10,
              "Gate 7′ (P3.3): precomputed iface SXY lower orbit drift");
      TEST_LE(iface_sxy_upper, 1.0e-10,
              "Gate 7′ (P3.3): precomputed iface SXY upper orbit drift");
      TEST_LE(iface_sxz_lower, 1.0e-10,
              "Gate 7′ (P3.3): precomputed iface SXZ lower orbit drift");
      TEST_LE(iface_sxz_upper, 1.0e-10,
              "Gate 7′ (P3.3): precomputed iface SXZ upper orbit drift");
      TEST_LE(bface_sxy_lower, 1.0e-10,
              "Gate 7′ (P3.4): precomputed bface SXY lower orbit drift "
              "(CRITICAL — if > 1e-10 but SXZ is ULP, invoke §8.5 F1)");
      TEST_LE(bface_sxy_upper, 1.0e-10,
              "Gate 7′ (P3.4): precomputed bface SXY upper orbit drift "
              "(CRITICAL — if > 1e-10 but SXZ is ULP, invoke §8.5 F1)");
      TEST_LE(bface_sxz_lower, 1.0e-10,
              "Gate 7′ (P3.4): precomputed bface SXZ lower orbit drift");
      TEST_LE(bface_sxz_upper, 1.0e-10,
              "Gate 7′ (P3.4): precomputed bface SXZ upper orbit drift");
   }

   std::cout << "\n-- Gate 8: within-outer-side boundary F_h orbit breakdown --\n";
   {
      // Classify a boundary face by which of the six outer sides of the
      // cube it lies on.  Within a side, the orbit action is the 2D
      // reflection / D4 of the square's in-plane coordinates.  A bug
      // that lives in single-side orbits has a different mechanism from
      // one that lives across sides.
      auto side_of = [&](real_t cx, real_t cy, real_t cz) -> const char * {
         const real_t eps = 1e-6 * kL;
         if (std::abs(cx)       < eps) { return "x=0";  }
         if (std::abs(cx - kL)  < eps) { return "x=L";  }
         if (std::abs(cy)       < eps) { return "y=0";  }
         if (std::abs(cy - kL)  < eps) { return "y=L";  }
         if (std::abs(cz)       < eps) { return "z=0";  }
         if (std::abs(cz - kL)  < eps) { return "z=L";  }
         return "unknown";
      };
      auto inplane_key = [&](const char *side, real_t cx, real_t cy, real_t cz) {
         constexpr real_t s = 1e3;
         auto r = [&](real_t v) { return std::round(v * s); };
         char buf[128];
         if (std::string(side) == "y=0" || std::string(side) == "y=L")
         {
            std::snprintf(buf, sizeof(buf), "%s|%g|%g", side,
                          r(std::abs(cx - 0.5 * kL)),
                          r(std::abs(cz - 0.5 * kL)));
         }
         else if (std::string(side) == "x=0" || std::string(side) == "x=L")
         {
            std::snprintf(buf, sizeof(buf), "%s|%g|%g", side,
                          r(std::abs(cy - 0.5 * kL)),
                          r(std::abs(cz - 0.5 * kL)));
         }
         else
         {
            std::snprintf(buf, sizeof(buf), "%s|%g|%g", side,
                          r(std::abs(cx - 0.5 * kL)),
                          r(std::abs(cy - 0.5 * kL)));
         }
         return std::string(buf);
      };

      const auto samples = CollectBoundaryFaceSamples(
         wave, mesh, bc, I_step1, kDt, BFaceProbeMode::GammaBaseline);
      std::map<std::string, std::vector<const FaceOrbitSample *>> by_side_orbit;
      std::map<std::string, int> side_count;
      for (const auto &s : samples)
      {
         const char *side = side_of(s.cx, s.cy, s.cz);
         side_count[side]++;
         const std::string k = inplane_key(side, s.cx, s.cy, s.cz);
         by_side_orbit[k].push_back(&s);
      }

      std::cout << "    boundary-face side count:";
      for (const auto &kv : side_count)
      { std::cout << "  " << kv.first << "=" << kv.second; }
      std::cout << "\n";

      // Per-side per-component worst orbit drift
      std::map<std::string, std::array<real_t, NUM_STATE>> per_side_worst;
      for (const auto &kv : by_side_orbit)
      {
         if (kv.second.size() < 2) { continue; }
         const std::string full_key = kv.first;
         const std::string side_name = full_key.substr(0, full_key.find('|'));
         auto &side_worst = per_side_worst[side_name];
         for (int c = 0; c < NUM_STATE; c++)
         {
            real_t lo = std::numeric_limits<real_t>::max();
            real_t hi = std::numeric_limits<real_t>::lowest();
            for (const auto *p : kv.second)
            {
               lo = std::min(lo, p->F_h[c]);
               hi = std::max(hi, p->F_h[c]);
            }
            const real_t scale = std::max({std::abs(lo), std::abs(hi), real_t(1.0)});
            const real_t drift = (hi - lo) / scale;
            side_worst[c] = std::max(side_worst[c], drift);
         }
      }

      real_t worst_within_side = 0.0;
      std::cout << "    per-side worst F_h orbit drift (within side):\n";
      for (const auto &kv : per_side_worst)
      {
         std::cout << "      " << kv.first << ":";
         for (int c = 0; c < NUM_STATE; c++)
         {
            std::cout << "  " << kCompName[c] << "="
                      << std::scientific << std::setprecision(2) << kv.second[c];
            worst_within_side = std::max(worst_within_side, kv.second[c]);
         }
         std::cout << "\n";
      }
      TEST_LE(worst_within_side, 1.0e-10,
              "Gate 8: within-side F_h orbits uniform on each outer-cube side");
   }

   std::cout << "\n-- Gate 9: Gamma vs Godunov lifted bface drift comparison --\n";
   {
      const auto bsym = RunADERNonFaultFaceAuditModed(
         wave, mesh, bc, I_step1, kDt, IfaceLiftMode::Symmetrized,
         BFaceLiftMode::GammaSymmetrized);
      // Repeat with the Godunov BC mode on FreeSurface.  Toggle the
      // wave operator's mode, rerun the lifted audit, restore.
      const FreeSurfaceBCMode saved_mode = wave.GetFreeSurfaceBCMode();
      wave.SetFreeSurfaceBCMode(FreeSurfaceBCMode::Godunov);
      const auto god_raw = RunADERNonFaultFaceAuditModed(
         wave, mesh, bc, I_step1, kDt, IfaceLiftMode::Raw,
         BFaceLiftMode::GammaRaw);  // GammaRaw → Godunov (since mode switched)
      const auto god_both = RunADERNonFaultFaceAuditModed(
         wave, mesh, bc, I_step1, kDt, IfaceLiftMode::Symmetrized,
         BFaceLiftMode::GammaSymmetrized);
      wave.SetFreeSurfaceBCMode(saved_mode);

      auto sig = [&](const Vector &v, int comp, bool upper) {
         return MaxSortedSignatureDrift(v, fault_adjacent, ndof_total,
                                        comp, upper);
      };
      std::cout << "    Gamma + both sym   : bface SXY lower="
                << std::scientific << std::setprecision(3)
                << sig(bsym.boundary_after_minv, SXY, false)
                << "  bface SXY upper=" << sig(bsym.boundary_after_minv, SXY, true)
                << "  bface SXZ lower=" << sig(bsym.boundary_after_minv, SXZ, false)
                << "  bface SXZ upper=" << sig(bsym.boundary_after_minv, SXZ, true)
                << "\n";
      std::cout << "    Godunov raw        : bface SXY lower="
                << sig(god_raw.boundary_after_minv, SXY, false)
                << "  bface SXY upper=" << sig(god_raw.boundary_after_minv, SXY, true)
                << "  bface SXZ lower=" << sig(god_raw.boundary_after_minv, SXZ, false)
                << "  bface SXZ upper=" << sig(god_raw.boundary_after_minv, SXZ, true)
                << "\n";
      std::cout << "    Godunov + both sym : bface SXY lower="
                << sig(god_both.boundary_after_minv, SXY, false)
                << "  bface SXY upper=" << sig(god_both.boundary_after_minv, SXY, true)
                << "  bface SXZ lower=" << sig(god_both.boundary_after_minv, SXZ, false)
                << "  bface SXZ upper=" << sig(god_both.boundary_after_minv, SXZ, true)
                << "\n";

      real_t godunov_both_total = 0.0;
      for (int comp : {SXY, SXZ})
      {
         for (bool upper : {false, true})
         {
            godunov_both_total = std::max(godunov_both_total,
                                          sig(god_both.boundary_after_minv, comp, upper));
         }
      }
      TEST_LE(godunov_both_total, 1.0e-10,
              "Gate 9: Godunov BC mode + both symmetrizations zero the bface drift");
   }

   std::cout << "\n-- Gate 10: zero-input boundary sanity --\n";
   {
      // Replay the lifted non-fault audit with I = 0.  All F_h should be
      // zero; any nonzero drift under I=0 would indicate the audit itself
      // produces spurious output (e.g. the bulk_bg_scaled or a precomputed
      // flux matrix is leaking).
      Vector I_zero(wave.Height());
      I_zero = 0.0;
      const auto zero_lift = RunADERNonFaultFaceAuditModed(
         wave, mesh, bc, I_zero, kDt, IfaceLiftMode::Raw,
         BFaceLiftMode::GammaRaw);
      const real_t zero_iface_max = zero_lift.interior_after_minv.Normlinf();
      const real_t zero_bface_max = zero_lift.boundary_after_minv.Normlinf();
      std::cout << "    with I=0: interior rhs max="
                << std::scientific << std::setprecision(3) << zero_iface_max
                << "  boundary rhs max=" << zero_bface_max << "\n";
      TEST_LE(zero_iface_max, 1.0e-20,
              "Gate 10a: zero input gives zero interior non-fault contribution");
      TEST_LE(zero_bface_max, 1.0e-20,
              "Gate 10b: zero input gives zero boundary contribution");
   }

   std::cout << "\n-- Gate 11: constant-input boundary D4 covariance probe --\n";
   {
      // Build I with per-DOF constants so that every face-quadrature's
      // I_self is identical regardless of shape-ordering: set each state
      // component to a distinct nonzero constant across every DOF.
      // If per-face F_h within one outer side is NOT orbit-uniform here,
      // the bug is inside FreeSurfaceTotal (since I_self is literally the
      // same constant 9-vector at every face QP by construction).  If it
      // IS orbit-uniform here but not under I_step1, the drift comes from
      // per-face I_self variation driven by the Q_step1 pattern, not the
      // flux.
      Vector I_const(wave.Height());
      for (int c = 0; c < NUM_STATE; c++)
      {
         const real_t val = 1.0e6 * (c + 1);
         for (int i = 0; i < ndof_total; i++)
         {
            I_const(c * ndof_total + i) = val;
         }
      }

      auto dump_side_drift = [&](const std::vector<FaceOrbitSample> &samples,
                                 const char *label)
      {
         auto side_of = [&](real_t cx, real_t cy, real_t cz) -> const char * {
            const real_t eps = 1e-6 * kL;
            if (std::abs(cx)      < eps) { return "x=0"; }
            if (std::abs(cx - kL) < eps) { return "x=L"; }
            if (std::abs(cy)      < eps) { return "y=0"; }
            if (std::abs(cy - kL) < eps) { return "y=L"; }
            if (std::abs(cz)      < eps) { return "z=0"; }
            if (std::abs(cz - kL) < eps) { return "z=L"; }
            return "unknown";
         };
         std::map<std::string, std::array<real_t, NUM_STATE>> per_side_max;
         std::map<std::string, std::array<real_t, NUM_STATE>> per_side_min;
         for (const auto &s : samples)
         {
            const std::string side = side_of(s.cx, s.cy, s.cz);
            if (per_side_max.find(side) == per_side_max.end())
            {
               for (int c = 0; c < NUM_STATE; c++)
               {
                  per_side_max[side][c] = std::numeric_limits<real_t>::lowest();
                  per_side_min[side][c] = std::numeric_limits<real_t>::max();
               }
            }
            auto &mx = per_side_max[side];
            auto &mn = per_side_min[side];
            for (int c = 0; c < NUM_STATE; c++)
            {
               mx[c] = std::max(mx[c], s.F_h[c]);
               mn[c] = std::min(mn[c], s.F_h[c]);
            }
         }
         // Use RELATIVE spread (absolute spread / max|F_h|) to separate
         // ULP noise from real orbit drift.  Constant-I F_h values reach
         // O(1e17–1e18) so absolute spreads in the 30–130 range are
         // 1e-16 relative — pure floating-point roundoff.
         std::cout << "    " << label << " per-side relative F_h drift on constant I:\n";
         real_t worst = 0.0;
         for (const auto &kv : per_side_max)
         {
            const auto &mx = kv.second;
            const auto &mn = per_side_min.at(kv.first);
            std::cout << "      " << kv.first << ":";
            for (int c = 0; c < NUM_STATE; c++)
            {
               const real_t scale =
                  std::max({std::abs(mx[c]), std::abs(mn[c]), real_t(1.0)});
               const real_t rel = (mx[c] - mn[c]) / scale;
               std::cout << "  " << kCompName[c] << "="
                         << std::scientific << std::setprecision(2) << rel;
               worst = std::max(worst, rel);
            }
            std::cout << "\n";
         }
         return worst;
      };

      const auto samples_gamma = CollectBoundaryFaceSamples(
         wave, mesh, bc, I_const, kDt, BFaceProbeMode::GammaBaseline);
      const real_t gamma_worst = dump_side_drift(samples_gamma, "Gamma  ");

      const auto samples_god = CollectBoundaryFaceSamples(
         wave, mesh, bc, I_const, kDt, BFaceProbeMode::GodunovBaseline);
      const real_t god_worst = dump_side_drift(samples_god, "Godunov");

      TEST_LE(gamma_worst, 1.0e-13,
              "Gate 11a: constant-I F_h within-side orbit uniform under Gamma BC");
      TEST_LE(god_worst,   1.0e-13,
              "Gate 11b: constant-I F_h within-side orbit uniform under Godunov BC");
   }

   std::cout << "\n-- Gate 12: boundary-face normal orientation audit --\n";
   {
      // For every boundary face, compute nor via CalcOrtho(J_F) and test
      // whether it points outward (away from Elem1's centroid).  A
      // physically outward-oriented `nor` for a boundary face must
      // satisfy nor . (c_face - c_elem1) > 0, because Elem1 lies
      // "inside" the domain and the face is on the outer boundary.
      int total = 0;
      int inward = 0;
      std::map<std::string, int> inward_per_side;
      auto side_of = [&](real_t cx, real_t cy, real_t cz) -> const char * {
         const real_t eps = 1e-6 * kL;
         if (std::abs(cx)      < eps) { return "x=0"; }
         if (std::abs(cx - kL) < eps) { return "x=L"; }
         if (std::abs(cy)      < eps) { return "y=0"; }
         if (std::abs(cy - kL) < eps) { return "y=L"; }
         if (std::abs(cz)      < eps) { return "z=0"; }
         if (std::abs(cz - kL) < eps) { return "z=L"; }
         return "unknown";
      };
      for (int f = 0; f < mesh.GetNumFaces(); f++)
      {
         FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
         if (!ftr || ftr->Elem2No >= 0) { continue; }
         const int e1 = ftr->Elem1No;

         Array<int> fv;
         mesh.GetFaceVertices(f, fv);
         real_t cx_f = 0.0, cy_f = 0.0, cz_f = 0.0;
         for (int v = 0; v < fv.Size(); v++)
         {
            cx_f += mesh.GetVertex(fv[v])[0];
            cy_f += mesh.GetVertex(fv[v])[1];
            cz_f += mesh.GetVertex(fv[v])[2];
         }
         cx_f /= fv.Size(); cy_f /= fv.Size(); cz_f /= fv.Size();

         Array<int> ev;
         mesh.GetElementVertices(e1, ev);
         real_t cx_e = 0.0, cy_e = 0.0, cz_e = 0.0;
         for (int v = 0; v < ev.Size(); v++)
         {
            cx_e += mesh.GetVertex(ev[v])[0];
            cy_e += mesh.GetVertex(ev[v])[1];
            cz_e += mesh.GetVertex(ev[v])[2];
         }
         cx_e /= ev.Size(); cy_e /= ev.Size(); cz_e /= ev.Size();

         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                  2 * wave.GetOrder());
         const IntegrationPoint &ip = ir.IntPoint(0);
         ftr->SetAllIntPoints(&ip);
         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0.0) { nor_vec /= nor_len; }

         const real_t out_x = cx_f - cx_e;
         const real_t out_y = cy_f - cy_e;
         const real_t out_z = cz_f - cz_e;
         const real_t dot = nor_vec(0)*out_x + nor_vec(1)*out_y + nor_vec(2)*out_z;
         ++total;
         if (dot < 0.0)
         {
            ++inward;
            inward_per_side[side_of(cx_f, cy_f, cz_f)]++;
         }
      }

      std::cout << "    boundary-face inward-normal count: "
                << inward << " of " << total << " total\n";
      for (const auto &kv : inward_per_side)
      {
         std::cout << "      side " << kv.first << ": "
                   << kv.second << " inward-pointing normals\n";
      }
      TEST_LE(static_cast<real_t>(inward), 0.0,
              "Gate 12: every boundary-face nor from CalcOrtho points outward");
   }

   std::cout << "\n-- Gate 14: constant-I lifted element rhs orbit drift --\n";
   {
      // Same constant I used in Gate 11.  Compute the lifted interior
      // and boundary rhs under constant I, apply M^{-1}, and measure
      // element-level sorted-signature drift on fault-adjacent tets.
      //
      // If this PASSES (drift ~ ULP): the per-DOF `rhs ±= w*shape*F_h`
      //   accumulation under the Kuhn mesh's D4 is correct — the only
      //   source of the non-constant I_step1 drift is that Q_step1 is
      //   sorted-sig orbit-identical but not per-DOF orbit-identical,
      //   which is physically fine (the next step's FE operator then
      //   preserves orbit symmetry in EXPECTED REFLECTED form rather
      //   than bit-equal, and the original uniformity test's 1e-10
      //   threshold is too tight for this regime).
      // If this FAILS: the Kuhn mesh's element-local DOF ordering
      //   breaks D4 at the lift step itself — independent of Q — and
      //   the fix must restructure the assembly to be DOF-ordering
      //   invariant.
      Vector I_const(wave.Height());
      for (int c = 0; c < NUM_STATE; c++)
      {
         const real_t val = 1.0e6 * (c + 1);
         for (int i = 0; i < ndof_total; i++)
         {
            I_const(c * ndof_total + i) = val;
         }
      }
      const auto lifted = RunADERNonFaultFaceAuditModed(
         wave, mesh, bc, I_const, kDt, IfaceLiftMode::Symmetrized,
         BFaceLiftMode::GammaRaw);

      auto sig = [&](const Vector &v, int comp, bool upper) {
         return MaxSortedSignatureDrift(v, fault_adjacent, ndof_total,
                                        comp, upper);
      };
      real_t iface_worst = 0.0, bface_worst = 0.0;
      for (int comp : {SXY, SXZ, SXX, SYY, SZZ})
      {
         for (bool upper : {false, true})
         {
            iface_worst = std::max(iface_worst,
                                   sig(lifted.interior_after_minv, comp, upper));
            bface_worst = std::max(bface_worst,
                                   sig(lifted.boundary_after_minv, comp, upper));
         }
      }
      std::cout << "    constant-I lifted drift: iface worst="
                << std::scientific << std::setprecision(3) << iface_worst
                << "  bface worst=" << bface_worst << "\n";
      TEST_LE(iface_worst, 1.0e-10,
              "Gate 14a: constant-I lifted interior rhs orbit-uniform (element-level)");
      TEST_LE(bface_worst, 1.0e-10,
              "Gate 14b: constant-I lifted boundary rhs orbit-uniform (element-level)");
   }

   std::cout << "\n-- Gate 15: affine trace reconstruction probe --\n";
   {
      const auto probe = RunAffineTraceProbe(wave, mesh, bc);
      std::cout << "    worst interior self trace error : "
                << std::scientific << std::setprecision(3)
                << probe.worst_self_err
                << "  on field " << AffineFieldName(probe.worst_self_field)
                << "\n";
      std::cout << "    worst interior nbr  trace error : "
                << probe.worst_nbr_err
                << "  on field " << AffineFieldName(probe.worst_nbr_field)
                << "\n";
      std::cout << "    worst interior self-nbr jump    : "
                << probe.worst_jump_err
                << "  on field " << AffineFieldName(probe.worst_jump_field)
                << "\n";
      std::cout << "    worst boundary self trace error : "
                << probe.worst_bdry_err
                << "  on field " << AffineFieldName(probe.worst_bdry_field)
                << "\n";
      TEST_LE(probe.worst_self_err, 1.0e-12,
              "Gate 15a: interior self trace reproduces affine field");
      TEST_LE(probe.worst_nbr_err, 1.0e-12,
              "Gate 15b: interior neighbor trace reproduces affine field");
      TEST_LE(probe.worst_jump_err, 1.0e-12,
              "Gate 15c: interior self/nbr traces agree on affine field");
      TEST_LE(probe.worst_bdry_err, 1.0e-12,
              "Gate 15d: boundary self trace reproduces affine field");
   }

   std::cout << "\n-- Gate 16: analytic-trace lift probe --\n";
   {
      const auto analytic_lift =
         RunAnalyticTraceLiftSweep(wave, mesh, bc, fault_adjacent);
      std::cout << "    worst analytic-trace interior lifted drift : "
                << std::scientific << std::setprecision(3)
                << analytic_lift.worst_iface_drift
                << "  on field "
                << AffineFieldName(analytic_lift.worst_iface_field) << "\n";
      std::cout << "    worst analytic-trace boundary lifted drift : "
                << analytic_lift.worst_bface_drift
                << "  on field "
                << AffineFieldName(analytic_lift.worst_bface_field) << "\n";
      TEST_LE(analytic_lift.worst_iface_drift, 1.0e-10,
              "Gate 16a: analytic interior traces remain orbit-uniform after lift");
      TEST_LE(analytic_lift.worst_bface_drift, 1.0e-10,
              "Gate 16b: analytic boundary traces remain orbit-uniform after lift");
   }

   // Phase 3 (§8.2): Gate 14′ — duplicate Gate 14's constant-I lifted
   // drift measurement but through PrecomputedFaceFluxes's §6.2
   // dispatch.  Runs only under SEAS_TEST_USE_PRECOMPUTED_FLUX=1.
   // Acceptance (plan §8.3 P3.5): iface and bface worst drift ≤ 1e-14.
   if (use_precomputed_env)
   {
      std::cout << "\n-- Gate 14′ (§8.3 P3.5): constant-I lifted drift "
                << "under precomputed path --\n";
      Vector I_const(wave.Height());
      for (int c = 0; c < NUM_STATE; c++)
      {
         const real_t val = 1.0e6 * (c + 1);
         for (int i = 0; i < ndof_total; i++)
         {
            I_const(c * ndof_total + i) = val;
         }
      }
      const auto lifted_prec =
         RunPrecomputedFluxLiftedAudit(wave, mesh, I_const, kDt);

      auto sig = [&](const Vector &v, int comp, bool upper) {
         return MaxSortedSignatureDrift(v, fault_adjacent, ndof_total,
                                        comp, upper);
      };
      real_t iface_worst_prec = 0.0, bface_worst_prec = 0.0;
      for (int comp : {SXY, SXZ, SXX, SYY, SZZ})
      {
         for (bool upper : {false, true})
         {
            iface_worst_prec = std::max(iface_worst_prec,
               sig(lifted_prec.interior_after_minv, comp, upper));
            bface_worst_prec = std::max(bface_worst_prec,
               sig(lifted_prec.boundary_after_minv, comp, upper));
         }
      }
      std::cout << "    precomputed constant-I lifted drift: iface worst="
                << std::scientific << std::setprecision(3)
                << iface_worst_prec
                << "  bface worst=" << bface_worst_prec << "\n";

      // R-001 round 3 diagnostic (2026-04-23): per-orbit drift dump.
      // Interpretation: if a single orbit / single component dominates,
      // the drift is fixture-localized (potentially a metric artifact).
      // If multiple orbits contribute, the drift is distributed and
      // reflects a systemic kernel/dispatch asymmetry.
      {
         const char *comp_name[9] = {"SXX","SYY","SZZ","SXY","SYZ","SXZ",
                                      "VX","VY","VZ"};
         int n_upper = 0, n_lower = 0;
         for (const auto &e : fault_adjacent)
         {
            if (e.cy > 0.5 * kL) { n_upper++; } else { n_lower++; }
         }
         std::cout << "    [R-001 r3 diag] fault_adjacent: "
                   << n_upper << " upper, " << n_lower << " lower\n";
         for (int comp : {SXY, SXZ, SXX, SYY, SZZ, VX, VY, VZ})
         {
            for (bool upper : {false, true})
            {
               real_t iface_d = sig(lifted_prec.interior_after_minv,
                                    comp, upper);
               real_t bface_d = sig(lifted_prec.boundary_after_minv,
                                    comp, upper);
               if (iface_d > 1e-14 || bface_d > 1e-14)
               {
                  std::cout << "    [R-001 r3] comp=" << comp_name[comp]
                            << (upper ? " upper" : " lower")
                            << "  iface=" << std::scientific
                            << std::setprecision(3) << iface_d
                            << "  bface=" << bface_d << "\n";
               }
            }
         }
      }
      TEST_LE(iface_worst_prec, 1.0e-14,
              "Gate 14′ (P3.5): precomputed interior constant-I lifted "
              "rhs orbit-uniform (element-level)");
      TEST_LE(bface_worst_prec, 1.0e-14,
              "Gate 14′ (P3.5): precomputed boundary constant-I lifted "
              "rhs orbit-uniform (element-level)");
   }

   std::cout << "\n-- Gate 13: raw per-face nor / I_self / F_h dump on x=0 side --\n";
   {
      Vector I_const(wave.Height());
      for (int c = 0; c < NUM_STATE; c++)
      {
         const real_t val = 1.0e6 * (c + 1);
         for (int i = 0; i < ndof_total; i++)
         {
            I_const(c * ndof_total + i) = val;
         }
      }
      real_t bulk_bg_scaled[NUM_STATE] = {0.0};
      int shown = 0;
      const int max_show = 6;
      for (int f = 0; f < mesh.GetNumFaces() && shown < max_show; f++)
      {
         FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
         if (!ftr || ftr->Elem2No >= 0) { continue; }
         Array<int> fv; mesh.GetFaceVertices(f, fv);
         real_t cx = 0.0;
         for (int v = 0; v < fv.Size(); v++) { cx += mesh.GetVertex(fv[v])[0]; }
         cx /= fv.Size();
         if (std::abs(cx) > 1e-6 * kL) { continue; }  // x=0 only

         const FiniteElement *fe1 = wave.GetFESpace().GetFE(ftr->Elem1No);
         const int ndof = fe1->GetDof();
         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                  2 * wave.GetOrder());
         const IntegrationPoint &ip = ir.IntPoint(0);
         ftr->SetAllIntPoints(&ip);
         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         Vector nor_unit = nor_vec;
         if (nor_len > 0.0) { nor_unit /= nor_len; }
         IntegrationPoint ip1;
         ftr->Loc1.Transform(ip, ip1);
         Vector shape1(ndof);
         fe1->CalcShape(ip1, shape1);
         real_t shape_sum = 0.0;
         for (int i = 0; i < ndof; i++) { shape_sum += shape1(i); }

         const int dof_offset1 = ftr->Elem1No * wave.GetNDof();
         real_t I_self[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            I_self[c] = 0.0;
            for (int i = 0; i < ndof; i++)
            {
               I_self[c] += shape1(i) *
                            I_const(c * ndof_total + dof_offset1 + i);
            }
         }
         real_t nor_arr[3] = {nor_unit(0), nor_unit(1), nor_unit(2)};
         real_t F_h[NUM_STATE];
         wave.GetFlux().FreeSurfaceTotal(nor_arr, I_self, bulk_bg_scaled, F_h);

         std::cout << "    f=" << f << "  e1=" << ftr->Elem1No
                   << "  nor=(" << std::fixed << std::setprecision(4)
                   << nor_unit(0) << "," << nor_unit(1) << "," << nor_unit(2)
                   << ")  nor_len=" << std::scientific << std::setprecision(3)
                   << nor_len
                   << "  shape_sum=" << std::fixed << std::setprecision(6)
                   << shape_sum
                   << "  I_self[SXX]=" << std::scientific << std::setprecision(3)
                   << I_self[SXX]
                   << "  F_h[SXX]=" << F_h[SXX]
                   << "  F_h[SXY]=" << F_h[SXY]
                   << "  F_h[SXZ]=" << F_h[SXZ]
                   << "\n";
         ++shown;
      }
   }

   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}
