// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_fields.hpp — fully dynamic BP5 cycle (plan rev 3, Section 3.3 b):
// RectTaperField and BuildFaultParams.
//
// RectTaperField is a generic field on a fault: a rectangle in the fault
// plane, given in global XYZ by a centre, two unit in-plane axes and two
// half-lengths, with a linear taper of width w in the max norm.  It carries
// no benchmark code; the deck gives every number.
//
// BuildFaultParams turns the deck's rate-and-state block and its fields into
// per-DOF parameters on the owned fault DOFs.  a and Dc come from
// RectTaperFields through the shared resolver's field callables
// (spatial::RateStateSidecarFields), so the resolver's validation and its
// eta = "auto" path run unchanged.  The initial slip rate has two kernel
// components: strike from a RectTaperField, dip a constant (V_zero).
//
// No solver header (rule M4): the DOF table comes from the QD adapter.

#ifndef MFEM_SEAS_CYCLE_FAULT_FIELDS_HPP
#define MFEM_SEAS_CYCLE_FAULT_FIELDS_HPP

#include "mfem.hpp"
#include "../spatial/code/spatial_friction.hpp"
#include "../dynamic/heterogeneous_material.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// A rectangle in the fault plane with a linear taper (deck [cycle.fields.*]).
struct RectTaperSpec
{
   std::array<real_t, 3> center = {{0.0, 0.0, 0.0}};  ///< global XYZ [m]
   std::array<real_t, 3> axis1 = {{1.0, 0.0, 0.0}};   ///< unit, in the fault plane
   std::array<real_t, 3> axis2 = {{0.0, 0.0, 1.0}};   ///< unit, in the fault plane
   real_t half1 = 0.0;      ///< half-length along axis1 [m]
   real_t half2 = 0.0;      ///< half-length along axis2 [m]
   real_t taper = 0.0;      ///< taper width w [m]; 0 = sharp edge
   real_t eps = 0.0;        ///< inclusion tolerance of a sharp edge [m]
   real_t inside = 0.0;     ///< value inside the rectangle
   real_t outside = 0.0;    ///< value beyond the taper
   bool set = false;        ///< the deck gave this field

   /// Refuses non-unit or non-orthogonal axes, negative lengths.
   void Validate(const std::string &name) const
   {
      auto dot = [](const std::array<real_t, 3> &u, const std::array<real_t, 3> &v)
      { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
      MFEM_VERIFY(std::abs(dot(axis1, axis1) - 1.0) < 1e-12 &&
                  std::abs(dot(axis2, axis2) - 1.0) < 1e-12,
                  "[cycle.fields." << name << "] axis1 and axis2 must be unit vectors");
      MFEM_VERIFY(std::abs(dot(axis1, axis2)) < 1e-12,
                  "[cycle.fields." << name << "] axis1 and axis2 must be orthogonal");
      MFEM_VERIFY(half1 >= 0.0 && half2 >= 0.0 && taper >= 0.0 && eps >= 0.0,
                  "[cycle.fields." << name << "] half1, half2, taper and eps must be >= 0");
   }
};

/// The field value at global point x:
///   r1 = |(x - c).e1| - h1,  r2 = |(x - c).e2| - h2,  r = max(r1, r2)
///   w = 0:  inside if r <= eps, else outside            (sharp rectangle)
///   w > 0:  inside + (outside - inside) * clamp(r / w, 0, 1)
/// The clamped ends return inside and outside exactly.
inline real_t RectTaperField(const RectTaperSpec &f, const real_t x[3])
{
   const real_t d[3] = {x[0] - f.center[0], x[1] - f.center[1], x[2] - f.center[2]};
   const real_t r1 = std::abs(d[0] * f.axis1[0] + d[1] * f.axis1[1] + d[2] * f.axis1[2])
                     - f.half1;
   const real_t r2 = std::abs(d[0] * f.axis2[0] + d[1] * f.axis2[1] + d[2] * f.axis2[2])
                     - f.half2;
   const real_t r = std::max(r1, r2);
   if (f.taper == 0.0) { return (r <= f.eps) ? f.inside : f.outside; }
   const real_t s = r / f.taper;
   if (s <= 0.0) { return f.inside; }
   if (s >= 1.0) { return f.outside; }
   return f.inside + s * (f.outside - f.inside);
}

/// The owned fault DOFs in the QD operator's owned order: global XYZ (3 per
/// DOF), the bulk element of each DOF and its attribute.
struct FaultDOFTable
{
   Vector coords;
   Array<int> elem;
   Array<int> attr;
   int Size() const { return elem.Size(); }
};

/// Deck fields of BuildFaultParams.
struct FaultFieldsConfig
{
   RectTaperSpec a;           ///< direct effect a
   RectTaperSpec Dc;          ///< critical slip distance
   RectTaperSpec V_init;      ///< initial slip rate, kernel strike component
   real_t V_init_dip = 0.0;   ///< initial slip rate, kernel dip component
};

/// Per-DOF parameters on the owned fault DOFs.
struct FaultParams
{
   spatial::RateStatePerDOFParams rs;  ///< a, b, Dc, f_0, V_0, eta, sigma_n_eff;
                                       ///< rs.V_init = |V_init|
   Vector V_init;                      ///< [2N]: (dip, strike) per DOF, kernel frame
};

/// BuildFaultParams (plan 3.3 b):
///   fields.a, fields.Dc = RectTaperField(cfg.a), RectTaperField(cfg.Dc)
///   fields.V_w          = constant (required by the struct; aging ignores it)
///   rs = SpatialFrictionResolver(fields).ResolveRateState(block, table, ...)
///   V_init(dip, strike) = (cfg.V_init_dip, RectTaperField(cfg.V_init))
/// `block` is the deck's [friction.rate_state]; its spatial rules still apply
/// after the fields, as in the resolver.  The fields are mandatory here, so a
/// deck cannot fall back to the scalar defaults by omission.
inline FaultParams BuildFaultParams(const FaultDOFTable &table,
                                    const spatial::RateStateBlock &block,
                                    const FaultFieldsConfig &cfg,
                                    const MaterialField &material,
                                    ParMesh &pmesh,
                                    const spatial::PorePressureSpec &pp)
{
   MFEM_VERIFY(cfg.a.set && cfg.Dc.set && cfg.V_init.set,
               "BuildFaultParams: [cycle.fields.a], [cycle.fields.Dc] and "
               "[cycle.fields.V_init] are all required");
   cfg.a.Validate("a");
   cfg.Dc.Validate("Dc");
   cfg.V_init.Validate("V_init");
   MFEM_VERIFY(block.state_evolution == spatial::StateEvolutionKind::AgingLaw,
               "BuildFaultParams: the QD kernel of this branch runs the aging law only");
   MFEM_VERIFY(!block.depth_profile.enabled,
               "BuildFaultParams: [friction.rate_state] depth_profile and the "
               "cycle fields both seed a; use one");

   const int N = table.Size();
   MFEM_VERIFY(table.coords.Size() == 3 * N && table.attr.Size() == N,
               "BuildFaultParams: inconsistent DOF table");

   auto fields = std::make_shared<spatial::RateStateSidecarFields>();
   const RectTaperSpec fa = cfg.a, fdc = cfg.Dc;
   fields->a = [fa](real_t x, real_t y, real_t z)
   { const real_t p[3] = {x, y, z}; return RectTaperField(fa, p); };
   fields->Dc = [fdc](real_t x, real_t y, real_t z)
   { const real_t p[3] = {x, y, z}; return RectTaperField(fdc, p); };
   const real_t V_w = block.V_w_default;
   fields->V_w = [V_w](real_t, real_t, real_t) { return V_w; };
   // A field defined everywhere: no hull limits.
   const real_t big = std::numeric_limits<real_t>::max();
   fields->bbox = {{-big, big, -big, big, -big, big}};

   spatial::SpatialFrictionResolver resolver(fields);
   FaultParams out;
   out.rs = resolver.ResolveRateState(block, table.coords, table.elem, table.attr,
                                      material, pmesh, pp, Vector());

   out.V_init.SetSize(2 * N);
   for (int i = 0; i < N; i++)
   {
      const real_t *x = table.coords.GetData() + 3 * i;
      const real_t vs = RectTaperField(cfg.V_init, x);
      const real_t vd = cfg.V_init_dip;
      MFEM_VERIFY(vs > 0.0 || vd > 0.0, "BuildFaultParams: V_init = 0 at DOF " << i);
      out.V_init(2 * i) = vd;
      out.V_init(2 * i + 1) = vs;
      out.rs.V_init(i) = std::sqrt(vd * vd + vs * vs);
   }
   // The aging law reads no V_w; an empty V_w keeps the QD geometry's
   // "aging leaves V_w empty" invariant (fault_geometry.hpp SetRateStatePerDOF).
   out.rs.V_w.SetSize(0);
   return out;
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_FAULT_FIELDS_HPP
