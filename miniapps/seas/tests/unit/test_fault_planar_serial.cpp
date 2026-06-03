// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 0 of the fault-dealiasing plan
// (document/fault_dealiasing_dev/seisol_overintegration_resample_speckle_2026-06-02.md):
// LOCAL REPRODUCTION HARNESS for the on-fault normal-stress speckle / drift.
//
// ============================================================================
// What this is
// ============================================================================
// A serial (MPI-capable), homogeneous, PLANAR-fault mini-rupture testbed large
// enough for a rupture front to propagate (a structured nx x nz x 2 cartesian
// tet slab with the fault on the mid-plane y = Ly/2 -> 2*nx*nz fault triangles;
// the 2-tet fixture used by the pepper-bug suite is far too small to develop the
// secular sigma_n drift).
//
// It reproduces, on a purely local problem at p1, the symptom the plan targets:
//   - max_F |sigma_n - sigma_n0|(t)  : the on-fault normal-stress excursion.
//                                      Analytically ZERO for a symmetric planar
//                                      strike-slip fault; in pure-upwind ADER it
//                                      grows secularly.  PRIMARY metric (stored
//                                      directly in DOFData::sigma_n_corr).
//   - max_F |[[v_n]]|(t)             : the genuine fault-normal velocity jump
//                                      |v_y(+) - v_y(-)| at the fault QPs (the
//                                      planar fault normal is +/- y_hat).  This
//                                      is the coherent leak that drives
//                                      delta sigma_n = eta_p * [[v_n]].
//
// It is the A/B testbed for Phases 1-3 (over-integration + resample).  The
// --fault-overint knob is LIVE (Phase 1): it raises the fault-flux quadrature
// degree to 2*(order+K) via WaveOperator::SetFaultOverint, decoupled from the
// bulk 2*order rule.  --fault-resample is still ACCEPTED but inert (Phases 2-3).
//
// ============================================================================
// Mesh variants (near-fault stencil symmetry)
// ============================================================================
// Both variants are the SAME nx x nz x 2 tet slab; they differ ONLY in the
// near-fault tet orientation across y = Ly/2:
//   - asymmetric (default): both y-layers use the same Kuhn 6-tet split
//        (translation-invariant) -> the +y / -y stencils are NOT mirror images.
//   - symmetric            : the upper y-layer's split is the y-mirror of the
//        lower (generalized from dynamic/d4_tet_mesh.hpp's orient flag) -> the
//        mesh is geometrically mirror-symmetric across the fault.
// Comparing the two probes the Zhang (2023) near-fault-upwind-dissipation
// sensitivity, which is a DISTINCT speckle source from the friction aliasing
// the plan's dealiasing cures (see plan Section 8).
//
// ============================================================================
// Usage
// ============================================================================
//   make seas_test_fault_planar_serial
//   ./seas_test_fault_planar_serial                       # serial baseline (p1)
//   ./seas_test_fault_planar_serial --order 1 --nx 8 --nz 8 --nsteps 8000
//   mpirun -np 8 ./seas_test_fault_planar_serial --nx 16 --nz 16  # bigger
//
// Exit 0 iff the acceptance checks pass.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_TRUE(cond, msg) do { \
   num_tests++; \
   const bool c = (cond); \
   if (g_rank == 0) { \
      if (c) { num_passed++; \
         std::cout << "  PASSED: " << msg << "\n"; } \
      else { num_failed++; \
         std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
   } \
} while (0)

namespace {

// ---------------------------------------------------------------------------
// Command-line options
// ---------------------------------------------------------------------------
struct Options {
   int    order   = 1;        // FE polynomial order (Phase 0 acceptance is p1)
   int    nx      = 8;        // in-plane cells (strike, x)
   int    nz      = 8;        // in-plane cells (dip, z)
   int    ny_half = 2;        // element layers per fault side; absorber ny_half*he
                              // away (REVIEW R-003: 1 was too close to the fault)
   real_t he      = 1000.0;   // element edge length [m]
   real_t dt      = 1.0e-4;   // time step [s]; CFL-calibrated for order=1, he=1 km.
                              // Auto-scaled with he/order in main() unless --dt is
                              // given (REVIEW R-001).
   bool   dt_user_set = false;// true once --dt is parsed (disables auto-scaling)
   int    nsteps  = 8000;     // steps (the secular drift accrues slowly vs CFL)
   int    print_every = 1000; // diagnostic print cadence
   bool   fault_overint = false;  // Phase 1 knob (inert in Phase 0)
   int    fault_overint_k = 0;
   bool   fault_resample = false; // Phase 2/3 knob (inert in Phase 0)
   bool   run_symmetry_probe = true; // run the asym-vs-sym comparison
   int    sym_nx = 6, sym_nz = 6;    // matched size for the comparison
   int    sym_nsteps = 4000;
};

void PrintUsage()
{
   std::cout <<
      "Usage: seas_test_fault_planar_serial [opts]\n"
      "  --order N        FE order (default 1; Phase 0 acceptance is p1)\n"
      "  --nx N --nz N    in-plane cells (default 8x8 -> 2*nx*nz fault faces)\n"
      "  --ny-half N      element layers per fault side (default 2; absorber N*he away)\n"
      "  --he L           element edge [m] (default 1000)\n"
      "  --dt T           time step [s] (default 1e-4 @ order1/he1km; auto-scaled\n"
      "                   with he/order unless set explicitly)\n"
      "  --nsteps N       number of ADER-2 steps (default 8000)\n"
      "  --print-every N  diagnostic cadence (default 1000)\n"
      "  --no-symmetry-probe   skip the asymmetric-vs-symmetric comparison\n"
      "  --fault-overint K     fault-flux over-integration factor (degree 2*(order+K))\n"
      "  --fault-resample      (Phase 2/3 knob; ACCEPTED but INERT)\n";
}

bool ParseArgs(int argc, char *argv[], Options &o)
{
   auto need = [&](int i) -> bool {
      if (i + 1 >= argc) {
         if (g_rank == 0) { std::cout << "ERROR: missing value for "
                                      << argv[i] << "\n"; }
         return false;
      }
      return true;
   };
   for (int i = 1; i < argc; i++) {
      std::string a = argv[i];
      if      (a == "--order"       && need(i)) { o.order  = std::atoi(argv[++i]); }
      else if (a == "--nx"          && need(i)) { o.nx     = std::atoi(argv[++i]); }
      else if (a == "--nz"          && need(i)) { o.nz     = std::atoi(argv[++i]); }
      else if (a == "--ny-half"     && need(i)) { o.ny_half = std::atoi(argv[++i]); }
      else if (a == "--he"          && need(i)) { o.he     = std::atof(argv[++i]); }
      else if (a == "--dt"          && need(i)) { o.dt     = std::atof(argv[++i]);
                                                  o.dt_user_set = true; }
      else if (a == "--nsteps"      && need(i)) { o.nsteps = std::atoi(argv[++i]); }
      else if (a == "--print-every" && need(i)) { o.print_every = std::atoi(argv[++i]); }
      else if (a == "--no-symmetry-probe")      { o.run_symmetry_probe = false; }
      else if (a == "--fault-overint" && need(i)) {
         o.fault_overint = true; o.fault_overint_k = std::atoi(argv[++i]);
      }
      else if (a == "--fault-resample")         { o.fault_resample = true; }
      else if (a == "--help" || a == "-h")      { if (g_rank==0) PrintUsage(); return false; }
      else {
         if (g_rank == 0) { std::cout << "ERROR: unknown arg '" << a << "'\n";
                            PrintUsage(); }
         return false;
      }
   }
   if (o.order < 1 || o.nx < 2 || o.nz < 2 || o.ny_half < 1 || o.he <= 0.0 ||
       o.dt <= 0.0 || o.nsteps < 1) {
      if (g_rank == 0) { std::cout << "ERROR: invalid option value(s)\n"; }
      return false;
   }
   return true;
}

// ---------------------------------------------------------------------------
// Cartesian tet fault slab, fault on the mid-plane y = Ly/2 (Ly = 2*ny_half*he,
// so the fault sits at y = ny_half*he with `ny_half` element layers per side and
// the absorbing boundary ny_half*he away from the fault — REVIEW R-003).
//
// Generalizes the add_hex 6-tet Kuhn split + orient flag from
// dynamic/d4_tet_mesh.hpp to an nx x nz in-plane tiling with 2*ny_half layers
// in y.
//   symmetric == false : every layer orient +1 (translation-invariant Kuhn)
//                        -> near-fault stencil is NOT mirror-symmetric.
//   symmetric == true  : lower half (iy < ny_half) orient +1, upper half
//                        orient -1 (the y-mirror dicing) -> mesh is geometrically
//                        mirror-symmetric across y = ny_half*he.  Because orient
//                        -1 is exactly the within-hex y-mirror of orient +1,
//                        reflecting any lower +1 hex across the fault plane lands
//                        a -1 hex in the matching upper slot, for any ny_half.
// Boundary attrs: fault triangles (y=Ly/2 interior) -> 3; all exterior -> 5
// (absorbing, so radiated waves leave cleanly).
// ---------------------------------------------------------------------------
Mesh BuildFaultSlab(int nx, int nz, real_t he, int ny_half, bool symmetric)
{
   const real_t Lx = nx * he, Lz = nz * he, Ly = 2.0 * ny_half * he;
   const int nyv = 2 * ny_half + 1;  // y-levels; fault plane at y = ny_half*he
   Mesh mesh(3, (nx + 1) * nyv * (nz + 1), 0, 0, 3);

   auto vid = [&](int ix, int iy, int iz) {
      return ix + (nx + 1) * (iy + nyv * iz);
   };
   for (int iz = 0; iz <= nz; iz++)
      for (int iy = 0; iy < nyv; iy++)
         for (int ix = 0; ix <= nx; ix++) {
            mesh.AddVertex(ix * he, iy * he, iz * he);
         }

   // 6-tet Kuhn split of the hex (ix,iy,iz)->(ix+1,iy+1,iz+1), orient lifted
   // verbatim from dynamic/d4_tet_mesh.hpp BuildD4Mesh::add_hex.
   auto add_hex = [&](int ix, int iy, int iz, int orient) {
      const int v000=vid(ix,iy,iz),  v100=vid(ix+1,iy,iz);
      const int v010=vid(ix,iy+1,iz),v110=vid(ix+1,iy+1,iz);
      const int v001=vid(ix,iy,iz+1),v101=vid(ix+1,iy,iz+1);
      const int v011=vid(ix,iy+1,iz+1),v111=vid(ix+1,iy+1,iz+1);
      auto v = [&](int k) -> int {
         if (orient > 0) {
            switch(k){case 0:return v000;case 1:return v100;case 2:return v010;case 3:return v110;
                      case 4:return v001;case 5:return v101;case 6:return v011;case 7:return v111;}
         } else {
            switch(k){case 0:return v010;case 1:return v110;case 2:return v000;case 3:return v100;
                      case 4:return v011;case 5:return v111;case 6:return v001;case 7:return v101;}
         }
         return -1;
      };
      int T[6][4] = {
         {v(0),v(1),v(3),v(7)},{v(0),v(2),v(3),v(7)},{v(0),v(2),v(6),v(7)},
         {v(0),v(4),v(6),v(7)},{v(0),v(4),v(5),v(7)},{v(0),v(1),v(5),v(7)}
      };
      for (int t = 0; t < 6; t++) { mesh.AddTet(T[t], 1); }
   };

   for (int iz = 0; iz < nz; iz++)
      for (int ix = 0; ix < nx; ix++)
         for (int iy = 0; iy < 2 * ny_half; iy++) {
            // Lower half (iy < ny_half) is the -y side; upper half the +y side.
            // symmetric -> upper half uses the y-mirror dicing (orient -1).
            const int orient = (iy < ny_half) ? +1 : (symmetric ? -1 : +1);
            add_hex(ix, iy, iz, orient);
         }

   mesh.FinalizeTopology();
   // AddTet orientations are not auto-fixed by Finalize (see d4_tet_mesh.hpp);
   // force positive Jacobians.  This is geometry-preserving (only swaps vertex
   // slots), so it does NOT break the geometric mirror symmetry.
   mesh.CheckElementOrientation(/*fix_it=*/true);
   mesh.Finalize();

   for (int f = 0; f < mesh.GetNumFaces(); f++) {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      real_t cy = 0.0;
      for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
      cy /= fv.Size();
      if (ftr && ftr->Elem2No >= 0) {
         if (std::abs(cy - 0.5 * Ly) < 1e-8) {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);   // FAULT
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 5);          // ABSORBING (exterior)
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   (void)Lx; (void)Lz;
   return mesh;
}

// Geometric y-mirror verification (for the symmetric variant): every vertex
// (x,y,z) must have a partner (x, Ly-y, z), and every element centroid likewise.
bool VerifyMirrorSymmetry(Mesh &mesh, real_t Ly, real_t tol)
{
   // Vertices.
   const int nv = mesh.GetNV();
   for (int i = 0; i < nv; i++) {
      const real_t *p = mesh.GetVertex(i);
      const real_t my = Ly - p[1];
      bool found = false;
      for (int j = 0; j < nv && !found; j++) {
         const real_t *q = mesh.GetVertex(j);
         if (std::abs(q[0]-p[0]) < tol && std::abs(q[1]-my) < tol &&
             std::abs(q[2]-p[2]) < tol) { found = true; }
      }
      if (!found) { return false; }
   }
   // Element centroids.
   const int ne = mesh.GetNE();
   std::vector<std::array<real_t,3>> cen(ne);
   for (int e = 0; e < ne; e++) {
      Array<int> ev; mesh.GetElementVertices(e, ev);
      real_t c[3] = {0,0,0};
      for (int v = 0; v < ev.Size(); v++) {
         const real_t *p = mesh.GetVertex(ev[v]);
         c[0]+=p[0]; c[1]+=p[1]; c[2]+=p[2];
      }
      cen[e] = {c[0]/ev.Size(), c[1]/ev.Size(), c[2]/ev.Size()};
   }
   for (int e = 0; e < ne; e++) {
      const real_t mx = cen[e][0], my = Ly - cen[e][1], mz = cen[e][2];
      bool found = false;
      for (int f = 0; f < ne && !found; f++) {
         if (std::abs(cen[f][0]-mx) < tol && std::abs(cen[f][1]-my) < tol &&
             std::abs(cen[f][2]-mz) < tol) { found = true; }
      }
      if (!found) { return false; }
   }
   return true;
}

// ---------------------------------------------------------------------------
// Fault setup (pattern from test_rupture_multistep_serial_vs_parallel.cpp):
// gather fault QP coords, init TPV102 rate-state DOFData, wire flux + state.
// Returns local fault QP count and fills fault_coords / dof_data / nqp_per_face.
// ---------------------------------------------------------------------------
template <typename MeshT>
int SetupFault(WaveOperator<MeshT> &wave, MeshT &mesh, int order,
               std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
               std::vector<Vector> &fault_coords, int &nqp_per_face)
{
   (void)order;  // fault quadrature degree now comes from wave.FaultFaceQuadDegree()
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();

   nqp_per_face = 0;
   if (int_faces.Size() > 0) {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault FTR null");
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), wave.FaultFaceQuadDegree()).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>) {
      if (nqp_per_face == 0 && shr_faces.Size() > 0) {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[0]);
         MFEM_VERIFY(ftr, "shared fault FTR null");
         nqp_per_face = IntRules.Get(ftr->GetGeometryType(), wave.FaultFaceQuadDegree()).GetNPoints();
      }
   }
#endif

   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++) {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), wave.FaultFaceQuadDegree());
      for (int q = 0; q < ir.GetNPoints(); q++) {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>) {
      for (int i = 0; i < shr_faces.Size(); i++) {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[i]);
         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), wave.FaultFaceQuadDegree());
         for (int q = 0; q < ir.GetNPoints(); q++) {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            Vector phys(3); ftr->Face->Transform(ip, phys);
            fault_coords.push_back(phys);
         }
      }
   }
#endif

   const int n_fault = (int_faces.Size() + shr_faces.Size()) * nqp_per_face;
   if (n_fault > 0) {
      InitializeFaultDOFs(dof_data, n_fault, fault_coords);
      ZeroDOFDataPreStressTotal(dof_data, n_fault);
   }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);
   return n_fault;
}

// Localized central overstress patch -> a rupture front that propagates
// outward (the front is the non-polynomial aliasing source, plan Section 3.1).
void ApplyPatchNucleation(std::vector<DOFData> &dof_data,
                          const std::vector<Vector> &fault_coords,
                          real_t cx, real_t cz, real_t r_nuc, real_t dtau)
{
   const int n = static_cast<int>(dof_data.size());
   for (int i = 0; i < n; i++) {
      const real_t dx = fault_coords[i](0) - cx;
      const real_t dz = fault_coords[i](2) - cz;
      const real_t r2 = dx*dx + dz*dz;
      dof_data[i].tau2_nuc = (r2 <= r_nuc * r_nuc) ? dtau : 0.0;
   }
}

// Evaluate a scalar bulk component (e.g. VY) of Q at an element reference IP.
template <typename FES>
real_t EvalComp(const FES &fes, const Vector &Q, int ndof_total, int comp,
                int e, const IntegrationPoint &eip)
{
   const FiniteElement *fe = fes.GetFE(e);
   Array<int> edofs; fes.GetElementDofs(e, edofs);
   Vector shape(fe->GetDof());
   fe->CalcShape(eip, shape);
   real_t val = 0.0;
   for (int j = 0; j < shape.Size(); j++) {
      val += shape[j] * Q(comp * ndof_total + edofs[j]);
   }
   return val;
}

// Genuine fault-normal velocity jump max_F |v_y(+) - v_y(-)| over interior
// fault faces (the planar fault normal is +/- y_hat).  Complete in serial;
// covers the interior subset under MPI (shared faces are skipped — the sigma_n
// metric still covers every QP via DOFData).
template <typename MeshT>
real_t MaxVnJumpInterior(WaveOperator<MeshT> &wave, MeshT &mesh,
                         const Vector &Q, int order)
{
   (void)order;  // fault QP rule now from wave.FaultFaceQuadDegree()
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   real_t worst = 0.0;
   for (int i = 0; i < int_faces.Size(); i++) {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      if (!ftr || ftr->Elem2No < 0) { continue; }
      // Determine +y / -y side by element centroid y.
      auto cy = [&](int e) {
         Array<int> ev; mesh.GetElementVertices(e, ev);
         real_t s = 0; for (int v=0;v<ev.Size();v++){ s += mesh.GetVertex(ev[v])[1]; }
         return s / ev.Size();
      };
      const real_t cy1 = cy(ftr->Elem1No), cy2 = cy(ftr->Elem2No);
      const int e_plus  = (cy1 > cy2) ? ftr->Elem1No : ftr->Elem2No;
      const int e_minus = (cy1 > cy2) ? ftr->Elem2No : ftr->Elem1No;
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), wave.FaultFaceQuadDegree());
      for (int q = 0; q < ir.GetNPoints(); q++) {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         const IntegrationPoint &eip1 = ftr->GetElement1IntPoint();
         const IntegrationPoint &eip2 = ftr->GetElement2IntPoint();
         const IntegrationPoint &eip_plus  = (e_plus  == ftr->Elem1No) ? eip1 : eip2;
         const IntegrationPoint &eip_minus = (e_minus == ftr->Elem1No) ? eip1 : eip2;
         const real_t vyp = EvalComp(fes, Q, ndof_total, VY, e_plus,  eip_plus);
         const real_t vym = EvalComp(fes, Q, ndof_total, VY, e_minus, eip_minus);
         worst = std::max(worst, std::abs(vyp - vym));
      }
   }
   return worst;
}

real_t GlobalMax(real_t local)
{
#ifdef MFEM_USE_MPI
   // The SEAS unit tests assume real_t == double (MPI_DOUBLE used throughout,
   // e.g. test_rupture_multistep_serial_vs_parallel.cpp).
   real_t g = local;
   MPI_Allreduce(&local, &g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   return g;
#else
   return local;
#endif
}

struct RunResult {
   real_t excursion_early = 0.0;   // max_F |sigma_n - sigma_n0| at first print
   real_t excursion_final = 0.0;   // ... at the last step
   real_t excursion_peak  = 0.0;   // ... peak over the run
   real_t vn_jump_final   = 0.0;   // max_F |[[v_n]]| at the last step
   real_t peak_slip_rate  = 0.0;   // confirms a rupture actually happened
};

// Run an nsteps ADER-2 rupture on a planar fault slab and record the
// sigma_n excursion / [[v_n]] time series.  Templated over Mesh / ParMesh.
template <typename MeshT>
RunResult RunRupture(MeshT &mesh, const Options &o, int nsteps,
                     const char *label, bool verbose)
{
   RunResult R;

   BoundaryConfig bc;
   bc.natural_attrs   = {};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {5};

   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;

   WaveOperator<MeshT> wave(mesh, o.order, TPV102Params::lambda,
                            TPV102Params::mu, TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bulk_bg);

   // CFL backstop (REVIEW R-001).  AdvanceADER performs no internal stability
   // check, and --he/--dt/--order are independent knobs.  ComputeMaxDt(1) =
   // h_min/cp (h_min = tet inscribed diameter ~ 0.41*he).  The default
   // dt=1e-4 @ order=1, he=1 km is very conservative (cfl ~ 1.5e-3 of
   // ComputeMaxDt(1); the ADER-2 stability bound is ~cfl 2), so a mis-sized
   // --he alone will not destabilize it — the original review overstated this.
   // main() still auto-scales dt with he/order so non-default runs keep that
   // margin (and higher orders get the 3/(2N+1) de-rating the p1-vs-p2
   // acceptance wants).  This guard catches only a grossly over-CFL explicit
   // --dt (or a hand-built mesh near the cfl~2 bound), aborting loud rather
   // than reporting an instability-driven sigma_n excursion as a reproduction.
   const real_t dt_cfl1 = wave.ComputeMaxDt(/*cfl=*/1.0);
   MFEM_VERIFY(o.dt <= 2.0 * dt_cfl1,
               "dt=" << o.dt << " s exceeds 2x the ADER-2 CFL reference "
               << dt_cfl1 << " s (= h_min/cp) for he=" << o.he
               << " m, order=" << o.order << " — near/over the cfl~2 "
               "stability bound.  The run would blow up and the sigma_n "
               "acceptance could pass off a numerical instability rather "
               "than the friction aliasing this harness must isolate.  Pass "
               "a smaller --dt, or omit --dt to use the auto-scaled value.");

   // Phase 1: enable fault-flux over-integration when requested.  When off, no
   // call is made, so the ctor's k=0 fault quadrature (degree 2*order) stands
   // ⇒ byte-identical fault QP set.  Must precede SetupFault so DOFData /
   // fault_coords are sized to the (possibly grown) per-face QP count.
   if (o.fault_overint) { wave.SetFaultOverint(o.fault_overint_k); }

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   int nqp_per_face = 0;
   const int n_fault = SetupFault(wave, mesh, o.order, dof_data, ff,
                                  fault_coords, nqp_per_face);

   int n_fault_global = n_fault;
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>) {
      MPI_Allreduce(&n_fault, &n_fault_global, 1, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD);
   }
#endif
   MFEM_VERIFY(n_fault_global > 0, "no fault QPs");

   // Central overstress patch (in-plane center, radius ~2 elements).
   const real_t cx = 0.5 * o.nx * o.he, cz = 0.5 * o.nz * o.he;
   const real_t r_nuc = 2.0 * o.he;
   if (n_fault > 0) {
      ApplyPatchNucleation(dof_data, fault_coords, cx, cz, r_nuc,
                           TPV102Params::nuc_dtau);
   }

   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();

   Vector Q(size);
   InitializeStateTotal(Q, ndof_total, TPV102Params::sigma_n,
                        TPV102Params::tau_ini);
   Vector Q_new(size);

   if (verbose && g_rank == 0) {
      std::cout << "\n-- Run [" << label << "]  order=" << o.order
                << "  mesh=" << o.nx << "x" << o.nz << "x2 tets"
                << "  fault QPs=" << n_fault_global
                << "  dt=" << std::scientific << std::setprecision(2) << o.dt
                << "  nsteps=" << nsteps << " --\n";
      std::cout << std::setw(8) << "step" << std::setw(10) << "t[s]"
                << std::setw(16) << "max|dSigN|[Pa]"
                << std::setw(16) << "max|[[v_n]]|"
                << std::setw(16) << "peakV[m/s]" << "\n";
   }

   bool early_recorded = false;
   for (int step = 0; step < nsteps; step++) {
      wave.AdvanceADER(Q, o.dt, /*ader_order=*/2, Q_new);
      Q.Swap(Q_new);

      // Per-step diagnostics.  NOTE: under the total-Q dispatch the background
      // normal stress lives in the bulk Q and DOFData::sigma_n0 is ZEROED by
      // ZeroDOFDataPreStressTotal; the excursion is therefore measured against
      // the PHYSICAL background TPV102Params::sigma_n (the constant the on-fault
      // sigma_n must hold), matching test_interior_fault_flux_path_asymmetric.cpp.
      const real_t sig_bg = TPV102Params::sigma_n;
      real_t loc_dsn = 0.0, loc_v = 0.0;
      for (int i = 0; i < n_fault; i++) {
         loc_dsn = std::max(loc_dsn,
                            std::abs(dof_data[i].sigma_n_corr - sig_bg));
         loc_v   = std::max(loc_v, std::abs(dof_data[i].slip_rate));
      }
      const real_t dsn = GlobalMax(loc_dsn);
      const real_t pkV = GlobalMax(loc_v);
      R.excursion_peak = std::max(R.excursion_peak, dsn);
      R.peak_slip_rate = std::max(R.peak_slip_rate, pkV);

      // "early" baseline for the secular-growth check: sampled at ~10% of the
      // run (NOT step 0, where the excursion is identically 0), so the growth
      // ratio final/early is a meaningful measure of accumulation.
      if (!early_recorded && (step + 1) >= std::max(1, nsteps / 10)) {
         R.excursion_early = dsn; early_recorded = true;
      }

      const bool do_print = verbose &&
         ((step % o.print_every) == 0 || step == nsteps - 1);
      if (do_print || step == nsteps - 1) {
         const real_t vnj = GlobalMax(MaxVnJumpInterior(wave, mesh, Q, o.order));
         R.excursion_final  = dsn;
         R.vn_jump_final    = vnj;
         if (do_print && g_rank == 0) {
            std::cout << std::setw(8) << step
                      << std::setw(10) << std::fixed << std::setprecision(4)
                      << (step + 1) * o.dt
                      << std::setw(16) << std::scientific << std::setprecision(4) << dsn
                      << std::setw(16) << vnj
                      << std::setw(16) << pkV << "\n";
         }
      }
   }
   return R;
}

} // anonymous namespace

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
#endif
   int nprocs = 1;
#ifdef MFEM_USE_MPI
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif

   Options o;
   if (!ParseArgs(argc, argv, o)) {
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 2;
   }

   // CFL-aware dt (REVIEW R-001).  The default dt (1e-4) is calibrated for the
   // validated config order=1, he=1 km (already deeply sub-CFL).  The stable
   // dt scales with element size (~h_min ∝ he) and with the DG order de-rating
   // 3/(2N+1); auto-scale so non-default --he / --order runs keep the same CFL
   // margin (the order de-rating in particular gives p2/p3 a proportionally
   // smaller dt, matching the plan's p1-vs-p2 acceptance).  An explicit --dt
   // is honoured verbatim (still bounded by the RunRupture CFL backstop).  At
   // the default order=1, he=1000 the scale factor is exactly 1.0 (byte-exact).
   if (!o.dt_user_set) {
      o.dt *= (o.he / 1000.0) * (3.0 / (2.0 * o.order + 1.0));
   }

   if (g_rank == 0) {
      std::cout << "\n=== Phase 0: planar-fault speckle/drift reproduction harness ==="
                << "\n  MPI ranks: " << nprocs
                << (nprocs == 1 ? " (serial Mesh — the Phase-0 acceptance path)"
                                : " (ParMesh — bigger reproduction)")
                << "\n";
      if (o.fault_overint) {
         std::cout << "  Phase 1: fault-flux OVER-INTEGRATION ON (k="
                   << o.fault_overint_k << "): fault quad degree "
                   << 2*(o.order + o.fault_overint_k) << " vs baseline "
                   << 2*o.order << ".\n";
      }
      if (o.fault_resample) {
         std::cout << "  NOTE: --fault-resample is ACCEPTED but INERT "
                   << "(wired in Phases 2-3).\n";
      }
   }

   // -----------------------------------------------------------------------
   // Part A — primary reproduction: large ASYMMETRIC slab, p1.
   // Acceptance: a growing sigma_n excursion (analytically zero).
   // -----------------------------------------------------------------------
   RunResult main_run;
   {
      Mesh smesh = BuildFaultSlab(o.nx, o.nz, o.he, o.ny_half, /*symmetric=*/false);
      if (nprocs == 1) {
         main_run = RunRupture<Mesh>(smesh, o, o.nsteps, "asymmetric/main", true);
      }
#ifdef MFEM_USE_MPI
      else {
         ParMesh pmesh(MPI_COMM_WORLD, smesh);
         main_run = RunRupture<ParMesh>(pmesh, o, o.nsteps, "asymmetric/main", true);
      }
#endif
   }

   if (g_rank == 0) {
      std::cout << "\n  [main] sigma_n excursion: early=" << std::scientific
                << std::setprecision(4) << main_run.excursion_early
                << "  final=" << main_run.excursion_final
                << "  peak="  << main_run.excursion_peak << " Pa\n"
                << "  [main] peak slip rate = " << main_run.peak_slip_rate
                << " m/s   final |[[v_n]]| = " << main_run.vn_jump_final
                << " m/s\n";
   }

   // A rupture must actually have happened (otherwise the test is vacuous).
   TEST_TRUE(main_run.peak_slip_rate > 1.0e-2,
             "a rupture front developed (peak slip rate > 1e-2 m/s)");
   // sigma_n must be a meaningful excursion: on a symmetric planar strike-slip
   // fault it is ANALYTICALLY ZERO, and at step 0 it is machine round-off
   // (|[[v_n]]| ~ 1e-17).  A > 1e5 Pa excursion is ~10 orders above that floor
   // and physically significant (it is a spurious leak from a constant field).
   TEST_TRUE(main_run.excursion_peak > 1.0e5,
             "on-fault sigma_n excursion reaches > 1e5 Pa — far above round-off, "
             "physically significant (sigma_n is analytically constant)");
   // ... and it must GROW secularly (one-way drift, not a bounded transient):
   // final > 3x the 10%-mark value.  This is the rectified accumulation of the
   // [[v_n]] leak (plan Section 3.4), the symptom the dealiasing must remove.
   TEST_TRUE(main_run.excursion_final > 3.0 * main_run.excursion_early &&
             main_run.excursion_early > 0.0,
             "sigma_n excursion GROWS secularly (final > 3x the 10%-mark value) "
             "at p1 — the rectified [[v_n]]-leak drift");

   // -----------------------------------------------------------------------
   // Part B — near-fault mesh-asymmetry sensitivity: asymmetric vs symmetric
   // at matched size.  Reports both; asserts symmetric is verified mirror and
   // that the two differ measurably.
   // -----------------------------------------------------------------------
   if (o.run_symmetry_probe) {
      Options os = o;
      os.nx = o.sym_nx; os.nz = o.sym_nz; os.print_every = o.sym_nsteps; // print only last

      Mesh asym = BuildFaultSlab(os.nx, os.nz, os.he, os.ny_half, /*symmetric=*/false);
      Mesh sym  = BuildFaultSlab(os.nx, os.nz, os.he, os.ny_half, /*symmetric=*/true);

      const real_t Ly = 2.0 * os.ny_half * os.he;
      const bool sym_ok = VerifyMirrorSymmetry(sym, Ly, 1e-6 * os.he);
      TEST_TRUE(sym_ok, "symmetric variant is geometrically y-mirror-symmetric");
      const bool asym_not = !VerifyMirrorSymmetry(asym, Ly, 1e-6 * os.he);
      TEST_TRUE(asym_not, "asymmetric variant is NOT y-mirror-symmetric (control)");

      RunResult ra, rs;
      if (nprocs == 1) {
         ra = RunRupture<Mesh>(asym, os, os.sym_nsteps, "asymmetric/probe", true);
         rs = RunRupture<Mesh>(sym,  os, os.sym_nsteps, "symmetric/probe",  true);
      }
#ifdef MFEM_USE_MPI
      else {
         ParMesh pa(MPI_COMM_WORLD, asym), ps(MPI_COMM_WORLD, sym);
         ra = RunRupture<ParMesh>(pa, os, os.sym_nsteps, "asymmetric/probe", true);
         rs = RunRupture<ParMesh>(ps, os, os.sym_nsteps, "symmetric/probe",  true);
      }
#endif
      if (g_rank == 0) {
         std::cout << "\n  [probe] asymmetric peak |dSigN| = " << std::scientific
                   << std::setprecision(4) << ra.excursion_peak
                   << " Pa   symmetric peak |dSigN| = " << rs.excursion_peak
                   << " Pa   ratio(asym/sym) = "
                   << (rs.excursion_peak > 0 ? ra.excursion_peak / rs.excursion_peak
                                             : 0.0) << "\n";
      }
      // Sensitivity exists: the two stencils produce measurably different drift
      // (we do NOT assert which is larger — friction aliasing also drives sym).
      const real_t rel_diff =
         std::abs(ra.excursion_peak - rs.excursion_peak) /
         std::max({ra.excursion_peak, rs.excursion_peak, real_t(1.0)});
      TEST_TRUE(rel_diff > 1.0e-3,
                "near-fault stencil symmetry measurably changes the drift "
                "(asym vs sym differ > 0.1%)");
   }

   if (g_rank == 0) {
      std::cout << "\n========================================\n"
                << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests << " tests\n"
                << "========================================\n";
   }

   int exit_code = (num_failed == 0) ? 0 : 1;
#ifdef MFEM_USE_MPI
   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
#endif
   return exit_code;
}
