// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// BP5 Parallel Smoke Test
//
// Tests the BP5 simulation pipeline without full ODE integration:
// 1. CreateBP5InlineMesh() with coarse resolution
// 2. ElasticityDomainOperator<ParMesh> with BR2
// 3. FaultGeometry<ParMesh> for BP5 parameters
// 4. RateStateFaultOperator<ParMesh, 2> with DR friction
// 5. SEAS operator: single Mult() call (= Solve + ComputeTraction + ComputeRHS)
// 6. Verify: no NaN/Inf, traction bounded, displacement bounded
//
// Usage: mpirun -np N ./seas_test_bp5_parallel_smoke

#include "mfem.hpp"
#include "../../solver/seas_operator.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp5_params.hpp"
#include "../../common/mpi_context.hpp"
#include "../../common/parallel_utils.hpp"

#include <iostream>
#include <cmath>
#include <memory>
#include <algorithm>
#include <array>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

std::unique_ptr<Mesh> LoadScaledBP5Mesh(const std::string &mesh_file,
                                        real_t mesh_scale)
{
   auto mesh = std::make_unique<Mesh>(mesh_file.c_str(), 1, 1);
   if (mesh_scale != 1.0)
   {
      for (int i = 0; i < mesh->GetNV(); i++)
      {
         real_t *v = mesh->GetVertex(i);
         v[0] *= mesh_scale;
         v[1] *= mesh_scale;
         v[2] *= mesh_scale;
      }
   }
   return mesh;
}

long long QuantizeCoord(real_t x, real_t scale = 1.0e9)
{
   return static_cast<long long>(std::llround(x * scale));
}

std::string MakeFaceLocalKey(const std::vector<std::pair<real_t, real_t>> &face_nodes,
                             real_t dof_x2, real_t dof_x3)
{
   std::array<std::pair<long long, long long>, 3> verts;
   for (int k = 0; k < 3; k++)
   {
      verts[k] = {QuantizeCoord(face_nodes[k].first),
                  QuantizeCoord(face_nodes[k].second)};
   }
   std::sort(verts.begin(), verts.end());

   std::ostringstream os;
   os << verts[0].first << ":" << verts[0].second << "|"
      << verts[1].first << ":" << verts[1].second << "|"
      << verts[2].first << ":" << verts[2].second << "|"
      << QuantizeCoord(dof_x2) << ":" << QuantizeCoord(dof_x3);
   return os.str();
}

struct BP5FieldEntry
{
   std::string key;
   real_t x2 = 0.0;
   real_t x3 = 0.0;
   real_t tau_dip = 0.0;
   real_t tau_strike = 0.0;
   real_t V_dip = 0.0;
   real_t V_strike = 0.0;
   real_t rate_dip = 0.0;
   real_t rate_strike = 0.0;
   real_t rate_psi = 0.0;
};

std::vector<BP5FieldEntry> BuildFaceLocalEntries(
   const Vector &x2, const Vector &x3,
   const Vector &traction, const Vector &slip_rate,
   const Vector &rate, int nbf_per_face)
{
   MFEM_VERIFY(nbf_per_face == 3,
               "This diagnostic currently assumes p=1 triangular fault faces");
   const int ndofs = x2.Size();
   MFEM_VERIFY(x3.Size() == ndofs, "Coordinate size mismatch");
   MFEM_VERIFY(traction.Size() == 2 * ndofs, "Traction size mismatch");
   MFEM_VERIFY(slip_rate.Size() == 2 * ndofs, "Slip-rate size mismatch");
   MFEM_VERIFY(rate.Size() == 3 * ndofs, "Rate size mismatch");
   MFEM_VERIFY(ndofs % nbf_per_face == 0, "DOF count not divisible by nbf");

   std::vector<BP5FieldEntry> entries;
   entries.reserve(ndofs);
   const int nfaces = ndofs / nbf_per_face;

   for (int f = 0; f < nfaces; f++)
   {
      std::vector<std::pair<real_t, real_t>> face_nodes(nbf_per_face);
      for (int k = 0; k < nbf_per_face; k++)
      {
         const int idx = f * nbf_per_face + k;
         face_nodes[k] = {x2(idx), x3(idx)};
      }

      for (int k = 0; k < nbf_per_face; k++)
      {
         const int idx = f * nbf_per_face + k;
         BP5FieldEntry e;
         e.key = MakeFaceLocalKey(face_nodes, x2(idx), x3(idx));
         e.x2 = x2(idx);
         e.x3 = x3(idx);
         e.tau_dip = traction(2 * idx);
         e.tau_strike = traction(2 * idx + 1);
         e.V_dip = slip_rate(2 * idx);
         e.V_strike = slip_rate(2 * idx + 1);
         e.rate_dip = rate(3 * idx);
         e.rate_strike = rate(3 * idx + 1);
         e.rate_psi = rate(3 * idx + 2);
         entries.push_back(e);
      }
   }

   std::sort(entries.begin(), entries.end(),
             [](const BP5FieldEntry &a, const BP5FieldEntry &b)
             {
                return a.key < b.key;
             });
   return entries;
}

std::vector<BP5FieldEntry> BuildDedupedParallelEntries(
   const Vector &x2, const Vector &x3,
   const Vector &traction, const Vector &slip_rate,
   const Vector &rate, int nbf_per_face)
{
   auto raw = BuildFaceLocalEntries(x2, x3, traction, slip_rate, rate, nbf_per_face);
   std::map<std::string, BP5FieldEntry> dedup;
   for (const auto &e : raw)
   {
      dedup.emplace(e.key, e);
   }

   std::vector<BP5FieldEntry> entries;
   entries.reserve(dedup.size());
   for (const auto &kv : dedup)
   {
      entries.push_back(kv.second);
   }
   return entries;
}

real_t MaxRelativeFieldDiff(const std::vector<BP5FieldEntry> &a,
                            const std::vector<BP5FieldEntry> &b,
                            std::string field_name)
{
   MFEM_VERIFY(a.size() == b.size(),
               "Entry count mismatch for field " << field_name);
   real_t max_rel = 0.0;
   for (size_t i = 0; i < a.size(); i++)
   {
      MFEM_VERIFY(a[i].key == b[i].key,
                  "Entry key mismatch at index " << i
                  << " for field " << field_name);

      real_t va = 0.0, vb = 0.0;
      if (field_name == "tau_dip") { va = a[i].tau_dip; vb = b[i].tau_dip; }
      else if (field_name == "tau_strike") { va = a[i].tau_strike; vb = b[i].tau_strike; }
      else if (field_name == "V_dip") { va = a[i].V_dip; vb = b[i].V_dip; }
      else if (field_name == "V_strike") { va = a[i].V_strike; vb = b[i].V_strike; }
      else if (field_name == "rate_dip") { va = a[i].rate_dip; vb = b[i].rate_dip; }
      else if (field_name == "rate_strike") { va = a[i].rate_strike; vb = b[i].rate_strike; }
      else if (field_name == "rate_psi") { va = a[i].rate_psi; vb = b[i].rate_psi; }
      else { MFEM_ABORT("Unknown field: " << field_name); }

      real_t denom = std::max({std::abs(va), std::abs(vb), 1.0e-14});
      max_rel = std::max(max_rel, std::abs(va - vb) / denom);
   }
   return max_rel;
}

real_t GetFieldValue(const BP5FieldEntry &e, const std::string &field_name)
{
   if (field_name == "tau_dip") { return e.tau_dip; }
   if (field_name == "tau_strike") { return e.tau_strike; }
   if (field_name == "V_dip") { return e.V_dip; }
   if (field_name == "V_strike") { return e.V_strike; }
   if (field_name == "rate_dip") { return e.rate_dip; }
   if (field_name == "rate_strike") { return e.rate_strike; }
   if (field_name == "rate_psi") { return e.rate_psi; }
   MFEM_ABORT("Unknown field: " << field_name);
   return 0.0;
}

real_t MaxDuplicateFieldSpread(const std::vector<BP5FieldEntry> &entries,
                               const std::string &field_name)
{
   real_t max_rel = 0.0;
   size_t i = 0;
   while (i < entries.size())
   {
      size_t j = i + 1;
      real_t vmin = GetFieldValue(entries[i], field_name);
      real_t vmax = vmin;
      while (j < entries.size() && entries[j].key == entries[i].key)
      {
         real_t v = GetFieldValue(entries[j], field_name);
         vmin = std::min(vmin, v);
         vmax = std::max(vmax, v);
         j++;
      }

      if (j - i > 1)
      {
         real_t denom = std::max({std::abs(vmin), std::abs(vmax), 1.0e-14});
         max_rel = std::max(max_rel, std::abs(vmax - vmin) / denom);
      }
      i = j;
   }
   return max_rel;
}

int CountDuplicateKeys(const std::vector<BP5FieldEntry> &entries)
{
   int duplicate_keys = 0;
   size_t i = 0;
   while (i < entries.size())
   {
      size_t j = i + 1;
      while (j < entries.size() && entries[j].key == entries[i].key)
      {
         j++;
      }
      if (j - i > 1)
      {
         duplicate_keys++;
      }
      i = j;
   }
   return duplicate_keys;
}

void ApplyDeterministicSlipPattern(const Vector &x2, const Vector &x3,
                                   const BP5Params &params,
                                   Vector &state)
{
   const int ndofs = x2.Size();
   MFEM_VERIFY(x3.Size() == ndofs, "Coordinate size mismatch");
   MFEM_VERIFY(state.Size() == 3 * ndofs, "State size mismatch for BP5");

   const real_t pi = M_PI;
   for (int i = 0; i < ndofs; i++)
   {
      const real_t sx = x2(i) / params.lf;
      const real_t sz = x3(i) / params.Wf;
      const real_t dip =
         0.01 * std::cos(pi * sx) * std::sin(0.5 * pi * sz);
      const real_t strike =
         0.05 * std::sin(pi * (sx + 0.5)) * std::sin(pi * sz);

      state(3 * i) = dip;
      state(3 * i + 1) = strike;
   }
}

void GatherArbitraryVectorToRoot(const Vector &local_data,
                                 Vector &global_data,
                                 MPIContext &mpi)
{
   std::vector<real_t> local_vec(local_data.Size());
   for (int i = 0; i < local_data.Size(); i++)
   {
      local_vec[i] = local_data(i);
   }

   std::vector<real_t> global_vec;
   GatherVectorToRoot(local_vec, global_vec, mpi.GetComm());

   if (mpi.IsRoot())
   {
      global_data.SetSize(static_cast<int>(global_vec.size()));
      for (int i = 0; i < global_data.Size(); i++)
      {
         global_data(i) = global_vec[static_cast<size_t>(i)];
      }
   }
   else
   {
      global_data.SetSize(0);
   }
}

struct ReferenceIPComparisonResult
{
   int init_duplicate_ok = 1;
   int init_compare_ok = 1;
   int pert_duplicate_ok = 1;
   int pert_compare_ok = 1;
};

ReferenceIPComparisonResult TestReferenceMeshIPSerialParallel(MPIContext &mpi)
{
   ReferenceIPComparisonResult result;

   if (mpi.IsRoot())
   {
      std::cout << "\nReference BP5 IP serial-vs-parallel check:\n";
   }

   const std::string mesh_file = "bp5/mesh/reference/bp5_tandem_coarse.msh";
   const real_t mesh_scale = 1000.0;
   const int order = 1;
   const DGMethod dg_method = DGMethod::IP;
   const real_t duplicate_tol = 1.0e-12;
   const real_t compare_tol = 1.0e-6;
   BP5Params params;

   std::vector<BP5FieldEntry> serial_init_entries;
   std::vector<BP5FieldEntry> serial_pert_entries;
   if (mpi.IsRoot())
   {
      auto serial_mesh = LoadScaledBP5Mesh(mesh_file, mesh_scale);

      ElasticityDomainOperator<Mesh> serial_domain(
         *serial_mesh, order, params.lambda(), params.mu(),
         params.Vp, params.Wf, params.lf, dg_method);
      FaultGeometry<Mesh> serial_geom(serial_domain, params);

      DieterichRuinaFriction::Constants fc;
      fc.V0 = params.V0;
      fc.f0 = params.f0;
      fc.b = params.b;
      fc.Dc = params.L0;
      DieterichRuinaFriction friction(fc);
      AgingLawPsi aging(params.b, params.V0, params.f0);

      RateStateFaultOperator<Mesh, 2> serial_fault(
         &serial_geom, &friction, &aging, params);
      BP5SEASOp serial_seas(&serial_domain, &serial_fault);

      Vector serial_state(serial_fault.StateSize());
      serial_seas.SetInitialCondition(serial_state);

      Vector serial_x2, serial_x3;
      serial_domain.GetFaultCoords2D(serial_x2, serial_x3);

      Vector serial_init_rate(serial_fault.StateSize());
      serial_fault.ComputeRHS(serial_seas.GetTraction(),
                              serial_state, serial_init_rate);
      serial_init_entries = BuildFaceLocalEntries(
         serial_x2, serial_x3, serial_seas.GetTraction(),
         serial_fault.GetSlipRate(), serial_init_rate,
         serial_domain.GetNbfPerFace());

      ApplyDeterministicSlipPattern(serial_x2, serial_x3, params, serial_state);
      serial_seas.SetTime(0.0);
      Vector serial_pert_rate(serial_fault.StateSize());
      serial_seas.Mult(serial_state, serial_pert_rate);
      serial_pert_entries = BuildFaceLocalEntries(
         serial_x2, serial_x3, serial_seas.GetTraction(),
         serial_fault.GetSlipRate(), serial_pert_rate,
         serial_domain.GetNbfPerFace());
   }

   auto parallel_mesh = LoadScaledBP5Mesh(mesh_file, mesh_scale);
   ParMesh pmesh(mpi.GetComm(), *parallel_mesh);
   parallel_mesh.reset();

   ElasticityDomainOperator<ParMesh> par_domain(
      pmesh, order, params.lambda(), params.mu(),
      params.Vp, params.Wf, params.lf, dg_method);
   FaultGeometry<ParMesh> par_geom(par_domain, params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;
   fc.f0 = params.f0;
   fc.b = params.b;
   fc.Dc = params.L0;
   DieterichRuinaFriction friction(fc);
   AgingLawPsi aging(params.b, params.V0, params.f0);

   RateStateFaultOperator<ParMesh, 2> par_fault(
      &par_geom, &friction, &aging, params, &mpi);
   PBP5SEASOp par_seas(&par_domain, &par_fault, &mpi);

   Vector par_state(par_fault.StateSize());
   par_seas.SetInitialCondition(par_state);

   Vector par_x2_local, par_x3_local;
   par_domain.GetFaultCoords2D(par_x2_local, par_x3_local);

   Vector par_init_rate_local(par_fault.StateSize());
   par_fault.ComputeRHS(par_seas.GetTraction(), par_state, par_init_rate_local);

   Vector par_x2_raw, par_x3_raw, par_trac_init_raw, par_slip_init_raw, par_rate_init_raw;
   par_geom.GatherToRoot(par_x2_local, par_x2_raw);
   par_geom.GatherToRoot(par_x3_local, par_x3_raw);
   GatherArbitraryVectorToRoot(par_seas.GetTraction(), par_trac_init_raw, mpi);
   GatherArbitraryVectorToRoot(par_fault.GetSlipRate(), par_slip_init_raw, mpi);
   GatherArbitraryVectorToRoot(par_init_rate_local, par_rate_init_raw, mpi);

   ApplyDeterministicSlipPattern(par_x2_local, par_x3_local, params, par_state);
   par_seas.SetTime(0.0);
   Vector par_pert_rate_local(par_fault.StateSize());
   par_seas.Mult(par_state, par_pert_rate_local);

   Vector par_trac_pert_raw, par_slip_pert_raw, par_rate_pert_raw;
   GatherArbitraryVectorToRoot(par_seas.GetTraction(), par_trac_pert_raw, mpi);
   GatherArbitraryVectorToRoot(par_fault.GetSlipRate(), par_slip_pert_raw, mpi);
   GatherArbitraryVectorToRoot(par_pert_rate_local, par_rate_pert_raw, mpi);

   if (mpi.IsRoot())
   {
      auto par_init_raw_entries = BuildFaceLocalEntries(
         par_x2_raw, par_x3_raw, par_trac_init_raw,
         par_slip_init_raw, par_rate_init_raw,
         par_domain.GetNbfPerFace());
      auto par_init_entries = BuildDedupedParallelEntries(
         par_x2_raw, par_x3_raw, par_trac_init_raw,
         par_slip_init_raw, par_rate_init_raw,
         par_domain.GetNbfPerFace());

      auto par_pert_raw_entries = BuildFaceLocalEntries(
         par_x2_raw, par_x3_raw, par_trac_pert_raw,
         par_slip_pert_raw, par_rate_pert_raw,
         par_domain.GetNbfPerFace());
      auto par_pert_entries = BuildDedupedParallelEntries(
         par_x2_raw, par_x3_raw, par_trac_pert_raw,
         par_slip_pert_raw, par_rate_pert_raw,
         par_domain.GetNbfPerFace());

      const int init_duplicate_keys = CountDuplicateKeys(par_init_raw_entries);
      const int pert_duplicate_keys = CountDuplicateKeys(par_pert_raw_entries);

      const real_t init_dup_tau = MaxDuplicateFieldSpread(par_init_raw_entries, "tau_strike");
      const real_t init_dup_v = MaxDuplicateFieldSpread(par_init_raw_entries, "V_strike");
      const real_t init_dup_rate = MaxDuplicateFieldSpread(par_init_raw_entries, "rate_strike");
      const real_t pert_dup_tau = MaxDuplicateFieldSpread(par_pert_raw_entries, "tau_strike");
      const real_t pert_dup_v = MaxDuplicateFieldSpread(par_pert_raw_entries, "V_strike");
      const real_t pert_dup_rate = MaxDuplicateFieldSpread(par_pert_raw_entries, "rate_strike");

      result.init_duplicate_ok =
         (init_dup_tau < duplicate_tol &&
          init_dup_v < duplicate_tol &&
          init_dup_rate < duplicate_tol) ? 1 : 0;
      result.pert_duplicate_ok =
         (pert_dup_tau < duplicate_tol &&
          pert_dup_v < duplicate_tol &&
          pert_dup_rate < duplicate_tol) ? 1 : 0;

      const bool init_size_ok = (par_init_entries.size() == serial_init_entries.size());
      const bool pert_size_ok = (par_pert_entries.size() == serial_pert_entries.size());

      real_t init_tau_rel = -1.0, init_v_rel = -1.0, init_rate_rel = -1.0;
      real_t pert_tau_rel = -1.0, pert_v_rel = -1.0, pert_rate_rel = -1.0;
      if (init_size_ok)
      {
         init_tau_rel = MaxRelativeFieldDiff(serial_init_entries, par_init_entries, "tau_strike");
         init_v_rel = MaxRelativeFieldDiff(serial_init_entries, par_init_entries, "V_strike");
         init_rate_rel = MaxRelativeFieldDiff(serial_init_entries, par_init_entries, "rate_strike");
      }
      if (pert_size_ok)
      {
         pert_tau_rel = MaxRelativeFieldDiff(serial_pert_entries, par_pert_entries, "tau_strike");
         pert_v_rel = MaxRelativeFieldDiff(serial_pert_entries, par_pert_entries, "V_strike");
         pert_rate_rel = MaxRelativeFieldDiff(serial_pert_entries, par_pert_entries, "rate_strike");
      }

      result.init_compare_ok =
         (init_size_ok &&
          init_tau_rel < compare_tol &&
          init_v_rel < compare_tol &&
          init_rate_rel < compare_tol) ? 1 : 0;
      result.pert_compare_ok =
         (pert_size_ok &&
          pert_tau_rel < compare_tol &&
          pert_v_rel < compare_tol &&
          pert_rate_rel < compare_tol) ? 1 : 0;

      std::cout << "  Init shared duplicates: " << init_duplicate_keys
                << ", max rel spread tau/V/rate = "
                << init_dup_tau << " / " << init_dup_v
                << " / " << init_dup_rate << "\n";
      std::cout << "  Init serial-parallel rel diff tau/V/rate = "
                << init_tau_rel << " / " << init_v_rel
                << " / " << init_rate_rel << "\n";
      std::cout << "  Pert shared duplicates: " << pert_duplicate_keys
                << ", max rel spread tau/V/rate = "
                << pert_dup_tau << " / " << pert_dup_v
                << " / " << pert_dup_rate << "\n";
      std::cout << "  Pert serial-parallel rel diff tau/V/rate = "
                << pert_tau_rel << " / " << pert_v_rel
                << " / " << pert_rate_rel << "\n";
   }

   BroadcastFromRoot(result.init_duplicate_ok, mpi.GetComm());
   BroadcastFromRoot(result.init_compare_ok, mpi.GetComm());
   BroadcastFromRoot(result.pert_duplicate_ok, mpi.GetComm());
   BroadcastFromRoot(result.pert_compare_ok, mpi.GetComm());

   return result;
}

}  // namespace

// ============================================================================
// Inline mesh creation (same as bp5_verification_full.cpp)
// ============================================================================

std::unique_ptr<Mesh> CreateBP5InlineMesh(
   int nx, int ny, int nz,
   real_t Lx, real_t Ly, real_t Lz)
{
   auto mesh = std::make_unique<Mesh>(
      Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz,
                            Element::HEXAHEDRON,
                            2.0 * Lx, 2.0 * Ly, Lz));

   for (int i = 0; i < mesh->GetNV(); i++)
   {
      real_t *v = mesh->GetVertex(i);
      v[0] -= Lx;
      v[1] -= Ly;
      v[2] -= Lz;  // Z ranges [-Lz, 0] (Tandem depth convention)
   }

   const real_t tol = 1e-6 * std::max({Lx, Ly, Lz});

   for (int i = 0; i < mesh->GetNBE(); i++)
   {
      Array<int> vertices;
      mesh->GetBdrElementVertices(i, vertices);

      real_t cx = 0.0, cy = 0.0, cz = 0.0;
      for (int j = 0; j < vertices.Size(); j++)
      {
         const real_t *v = mesh->GetVertex(vertices[j]);
         cx += v[0]; cy += v[1]; cz += v[2];
      }
      cx /= vertices.Size();
      cy /= vertices.Size();
      cz /= vertices.Size();

      int attr;
      // Tandem convention: top (z=0) and bottom (z=-Lz) → Natural (attr 1)
      // Far-field sides → Dirichlet (attr 5)
      if (std::abs(cz) < tol || std::abs(cz + Lz) < tol)
      {
         attr = 1;  // Natural (top/bottom)
      }
      else
      {
         attr = 5;  // Dirichlet (far-field)
      }

      mesh->SetBdrAttribute(i, attr);
   }

   mesh->SetAttributes();
   return mesh;
}

// ============================================================================
// Test macros
// ============================================================================

static int num_passed = 0;
static int num_failed = 0;

#define TEST_CHECK(ctx, name, cond) \
   do { \
      bool ok = (cond); \
      int local_ok = ok ? 1 : 0; \
      int global_ok = ctx.GlobalMinInt(local_ok); \
      if (ctx.IsRoot()) { \
         if (global_ok) { \
            std::cout << "  PASS: " << (name) << "\n"; \
            num_passed++; \
         } else { \
            std::cout << "  FAIL: " << (name) << "\n"; \
            num_failed++; \
         } \
      } \
   } while (0)

// ============================================================================
// Main test
// ============================================================================

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);

   if (mpi.IsRoot())
   {
      std::cout << "BP5 Parallel Smoke Test (NP=" << mpi.Size() << ")\n";
      std::cout << std::string(50, '=') << "\n";
   }

   // ========================================================================
   // Create coarse mesh
   // ========================================================================
   BP5Params params;

   // Use small mesh matching parallel_elasticity tests (which work reliably)
   // with Wf/lf sized to capture fault DOFs at these dimensions.
   int nx = 2, ny = 1, nz = 1;
   real_t Lx = 4.0, Ly = 2.0, Lz = 2.0;

   auto serial_mesh = CreateBP5InlineMesh(nx, ny, nz, Lx, Ly, Lz);

   if (mpi.IsRoot())
   {
      std::cout << "  Serial mesh: " << serial_mesh->GetNE()
                << " elements\n" << std::flush;
   }

   ParMesh pmesh(mpi.GetComm(), *serial_mesh);
   serial_mesh.reset();

   // GetGlobalNE() may be collective — call on all ranks
   long long global_ne = pmesh.GetGlobalNE();
   if (mpi.IsRoot())
   {
      std::cout << "  ParMesh: " << global_ne
                << " global elements\n";
   }

   // ========================================================================
   // Domain operator
   // ========================================================================
   int order = 1;
   // Use mesh-matching Wf/lf to ensure fault DOFs are detected on small mesh
   real_t Wf = Lz, lf_domain = 2.0 * Ly;
   ElasticityDomainOperator<ParMesh> domain(
      pmesh, order, params.lambda(), params.mu(),
      params.Vp, Wf, lf_domain, DGMethod::BR2);

   int local_fault_dofs = domain.GetNumFaultDOFs();
   int global_fault_dofs = mpi.GlobalSumInt(local_fault_dofs);

   if (mpi.IsRoot())
   {
      std::cout << "  Global fault DOFs: " << global_fault_dofs << "\n";
   }

   TEST_CHECK(mpi, "Fault DOFs detected", global_fault_dofs > 0);

   // ========================================================================
   // Fault components
   // ========================================================================
   FaultGeometry<ParMesh> fault_geom(domain, params, &mpi);

   DieterichRuinaFriction::Constants fc;
   fc.V0 = params.V0;
   fc.f0 = params.f0;
   fc.b = params.b;
   fc.Dc = params.L0;
   DieterichRuinaFriction friction(fc);
   AgingLawPsi aging(params.b, params.V0, params.f0);

   RateStateFaultOperator<ParMesh, 2> fault_op(
      &fault_geom, &friction, &aging, params, &mpi);

   // ========================================================================
   // SEAS operator + initial condition
   // ========================================================================
   PBP5SEASOp seas_op(&domain, &fault_op, &mpi);

   Vector state(fault_op.StateSize());
   seas_op.SetInitialCondition(state);

   real_t V_init = seas_op.GetMaxSlipRate();

   if (mpi.IsRoot())
   {
      std::cout << "  V_init = " << V_init << " m/s\n";
      std::cout << "  StateSize = " << fault_op.StateSize() << "\n";
   }

   TEST_CHECK(mpi, "Initial V_max finite",
              std::isfinite(V_init) && V_init >= 0.0);

   // ========================================================================
   // Check initial traction is bounded
   // ========================================================================
   {
      const Vector &trac = seas_op.GetTraction();
      bool trac_ok = true;
      real_t max_tau = 0.0;
      for (int i = 0; i < trac.Size(); i++)
      {
         if (!std::isfinite(trac(i)))
         {
            trac_ok = false;
            break;
         }
         max_tau = std::max(max_tau, std::abs(trac(i)));
      }
      real_t global_max_tau = mpi.GlobalMax(max_tau);

      TEST_CHECK(mpi, "Initial traction finite", trac_ok);
      TEST_CHECK(mpi, "Initial traction bounded (< 1 GPa)",
                 global_max_tau < 1e9);

      if (mpi.IsRoot())
      {
         std::cout << "  Max initial |tau| = "
                   << global_max_tau / 1e6 << " MPa\n";
      }
   }

   // ========================================================================
   // Single SEAS Mult() call — exercises Solve + ComputeTraction + ComputeRHS
   // This is the critical path that blows up on Frontera.
   // ========================================================================
   if (mpi.IsRoot())
   {
      std::cout << "\nSingle Mult() call:\n";
   }

   seas_op.SetTime(0.0);
   Vector rate(fault_op.StateSize());
   rate = 0.0;
   seas_op.Mult(state, rate);

   // Check rate vector for NaN/Inf
   {
      bool rate_ok = true;
      for (int i = 0; i < rate.Size(); i++)
      {
         if (!std::isfinite(rate(i))) { rate_ok = false; break; }
      }
      int local_ok = rate_ok ? 1 : 0;
      int global_ok = mpi.GlobalMinInt(local_ok);
      TEST_CHECK(mpi, "Rate vector finite after Mult()", global_ok == 1);
   }

   // Check traction after Mult()
   {
      const Vector &trac = seas_op.GetTraction();
      real_t max_tau = 0.0;
      bool trac_ok = true;
      for (int i = 0; i < trac.Size(); i++)
      {
         if (!std::isfinite(trac(i))) { trac_ok = false; break; }
         max_tau = std::max(max_tau, std::abs(trac(i)));
      }
      real_t global_max_tau = mpi.GlobalMax(max_tau);

      TEST_CHECK(mpi, "Post-Mult traction finite", trac_ok);
      TEST_CHECK(mpi, "Post-Mult traction bounded (< 1 GPa)",
                 global_max_tau < 1e9);

      if (mpi.IsRoot())
      {
         std::cout << "  Max |tau| after Mult = "
                   << global_max_tau / 1e6 << " MPa\n";
      }
   }

   // Check displacement after Mult()
   {
      const auto &u = seas_op.GetDisplacement();
      real_t u_max = 0.0;
      for (int i = 0; i < u.Size(); i++)
      {
         u_max = std::max(u_max, std::abs(u(i)));
      }
      real_t global_u_max = mpi.GlobalMax(u_max);

      TEST_CHECK(mpi, "Displacement bounded (< 1e6 m)",
                 global_u_max < 1e6);

      if (mpi.IsRoot())
      {
         std::cout << "  ||u||_inf = " << global_u_max << " m\n";
      }
   }

   // Check V_max after Mult()
   {
      real_t V_max = seas_op.GetMaxSlipRate();
      TEST_CHECK(mpi, "V_max finite after Mult()",
                 std::isfinite(V_max) && V_max >= 0.0);
      TEST_CHECK(mpi, "V_max bounded (< 1e4 m/s)", V_max < 1e4);

      if (mpi.IsRoot())
      {
         std::cout << "  V_max = " << V_max << " m/s\n";
      }
   }

   // ========================================================================
   // Explicit Euler step: state += dt * rate, then call Mult() again
   // Tests that a second evaluation doesn't blow up.
   // ========================================================================
   if (mpi.IsRoot())
   {
      std::cout << "\nExplicit Euler step (dt=0.1s):\n";
   }

   real_t dt_euler = 0.1;
   Vector state2(state);
   state2.Add(dt_euler, rate);

   seas_op.SetTime(dt_euler);
   Vector rate2(fault_op.StateSize());
   rate2 = 0.0;
   seas_op.Mult(state2, rate2);

   // Check after second Mult()
   {
      bool rate_ok = true;
      for (int i = 0; i < rate2.Size(); i++)
      {
         if (!std::isfinite(rate2(i))) { rate_ok = false; break; }
      }
      int local_ok = rate_ok ? 1 : 0;
      int global_ok = mpi.GlobalMinInt(local_ok);
      TEST_CHECK(mpi, "Rate2 finite after Euler step", global_ok == 1);

      const Vector &trac = seas_op.GetTraction();
      real_t max_tau = 0.0;
      for (int i = 0; i < trac.Size(); i++)
      {
         max_tau = std::max(max_tau, std::abs(trac(i)));
      }
      real_t global_max_tau = mpi.GlobalMax(max_tau);

      TEST_CHECK(mpi, "Post-Euler traction bounded (< 1 GPa)",
                 global_max_tau < 1e9);

      real_t V_max = seas_op.GetMaxSlipRate();
      TEST_CHECK(mpi, "Post-Euler V_max bounded (< 1e4 m/s)",
                 std::isfinite(V_max) && V_max < 1e4);

      if (mpi.IsRoot())
      {
         std::cout << "  Max |tau| = " << global_max_tau / 1e6
                   << " MPa, V_max = " << V_max << " m/s\n";
      }
   }

   // ========================================================================
   // Reference BP5 IP serial-vs-parallel comparison on the Tandem mesh
   // ========================================================================
   auto ref_ip_compare = TestReferenceMeshIPSerialParallel(mpi);
   TEST_CHECK(mpi, "BP5 IP init shared-face duplicates agree",
              ref_ip_compare.init_duplicate_ok == 1);
   TEST_CHECK(mpi, "BP5 IP init serial matches parallel",
              ref_ip_compare.init_compare_ok == 1);
   TEST_CHECK(mpi, "BP5 IP perturbed shared-face duplicates agree",
              ref_ip_compare.pert_duplicate_ok == 1);
   TEST_CHECK(mpi, "BP5 IP perturbed serial matches parallel",
              ref_ip_compare.pert_compare_ok == 1);

   // ========================================================================
   // Summary
   // ========================================================================
   if (mpi.IsRoot())
   {
      std::cout << "\n" << std::string(50, '=') << "\n";
      std::cout << "Results: " << num_passed << " passed, "
                << num_failed << " failed\n";
   }

   return num_failed;
}
