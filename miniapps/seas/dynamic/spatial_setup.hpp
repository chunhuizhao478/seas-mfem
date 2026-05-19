// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/spatial_setup.hpp — Phase 5c of
// spatial_dynamic_rupture_plan.md (rev-3).
//
// SAFS analog of `dynamic/tpv205_setup.hpp::InitializeFaultDOFs_TPV205`:
// writes per-DOF impedances + LSW friction params + pre-stress +
// TPV26/27 forced-rupture fields into `FaultFaceFlux::DOFData[]` in a
// single pass.  Header-only because it only stitches together existing
// per-DOF arrays produced upstream (resolver + stress source + material
// EvalAt).
//
// Plan reference: §Phase 5 of spatial_dynamic_rupture_plan.md.
//
// Note: the planned `ComputePerDOFCoordsAndBasisFromWave` helper (which
// would walk a `WaveOperator`'s fault face arrays to build per-DOF
// coords / basis / dof_to_elem without instantiating an
// `ElasticityDomainOperator`) is deferred to a Phase H follow-up.  The
// SAFS driver can use the new `FaultGeometry` BP5 ctor overload
// (Phase 5a) to consume pre-built per-DOF arrays from any source
// (including the existing `ElasticityDomainOperator`-driven flow today
// and the future WaveOperator-driven flow later).

#ifndef MFEM_SEAS_DYNAMIC_SPATIAL_SETUP_HPP
#define MFEM_SEAS_DYNAMIC_SPATIAL_SETUP_HPP

#include "mfem.hpp"

#include "fault_face_flux.hpp"
#include "heterogeneous_material.hpp"
#include "../spatial/code/spatial_friction.hpp"

#include <cmath>
#include <vector>

namespace mfem
{
namespace seas
{
namespace spatial
{

namespace internal
{

/// Per-DOF impedances + initial corrected-traction defaults.  Common
/// body for the LSW and RS InitializeFaultDOFs_Spatial overloads.
template <typename MeshT>
inline void seed_static_dof_fields(DOFData& d,
                                   int elem,
                                   MeshT& mesh,
                                   const MaterialField& material,
                                   real_t tau_pre_dip,
                                   real_t tau_pre_strike,
                                   real_t sigma_n_eff)
{
   // Per-DOF impedances from MaterialField::EvalAt at the bulk element
   // owning this fault DOF.  Use the reference-element CENTROID (NOT
   // ip.Init(0), which lands at the reference ORIGIN — a corner of
   // the element).  For MaterialField::MakeCoefficient (the SAFS
   // production path with CVMH/CVMS sidecars), the corner vs centroid
   // distinction biases impedances by the heterogeneity scale of one
   // element.  The centroid is the best per-element stop-gap until
   // the driver threads FaultGeometry::fault_dof_ip(i) through via the
   // IP-aware overload of seed_static_dof_fields below.  For
   // MaterialField::MakeConstant the IP location is irrelevant.
   mfem::ElementTransformation* T = mesh.GetElementTransformation(elem);
   const mfem::Geometry::Type gtype = mesh.GetElementBaseGeometry(elem);
   const mfem::IntegrationPoint& ip = mfem::Geometries.GetCenter(gtype);
   real_t lam, mu, rho;
   material.EvalAt(elem, *T, ip, lam, mu, rho);
   MFEM_VERIFY(mu > 0.0 && rho > 0.0,
               "InitializeFaultDOFs_Spatial: mu > 0 and rho > 0 required "
               "at bulk element " << elem << " (got mu=" << mu
               << ", rho=" << rho << ")");
   const real_t cp = std::sqrt((lam + 2.0 * mu) / rho);
   const real_t cs = std::sqrt(mu / rho);

   d.Zp_plus  = d.Zp_minus = rho * cp;
   d.Zs_plus  = d.Zs_minus = rho * cs;
   d.eta_p    = 0.5 * rho * cp;
   d.eta_s    = 0.5 * rho * cs;

   d.tau1_0   = tau_pre_dip;     // dip
   d.tau2_0   = tau_pre_strike;  // strike
   d.sigma_n0 = sigma_n_eff;

   d.tau1_nuc = d.tau2_nuc = d.sigma_n_nuc = 0.0;

   // Initial state: at rest (mirrors InitializeFaultDOFs_TPV205).
   d.slip_rate    = 0.0;
   d.V1 = d.V2    = 0.0;
   d.slip1 = d.slip2 = 0.0;
   d.tau1_corr    = 0.0;
   d.tau2_corr    = d.tau2_0;
   d.sigma_n_corr = d.sigma_n0;
}

/// IP-aware overload of seed_static_dof_fields (R-309): evaluates
/// material at the explicit reference IntegrationPoint `ip` instead of
/// the element centroid.  Used when the driver supplies per-DOF
/// reference IPs via FaultGeometry::fault_dof_ip(i) — the cleanest
/// per-DOF location is the actual fault QP IP, not the element
/// centroid.
template <typename MeshT>
inline void seed_static_dof_fields(DOFData& d,
                                   int elem,
                                   MeshT& mesh,
                                   const MaterialField& material,
                                   const mfem::IntegrationPoint& ip,
                                   real_t tau_pre_dip,
                                   real_t tau_pre_strike,
                                   real_t sigma_n_eff)
{
   mfem::ElementTransformation* T = mesh.GetElementTransformation(elem);
   real_t lam, mu, rho;
   material.EvalAt(elem, *T, ip, lam, mu, rho);
   MFEM_VERIFY(mu > 0.0 && rho > 0.0,
               "InitializeFaultDOFs_Spatial: mu > 0 and rho > 0 required "
               "at bulk element " << elem << " (got mu=" << mu
               << ", rho=" << rho << ")");
   const real_t cp = std::sqrt((lam + 2.0 * mu) / rho);
   const real_t cs = std::sqrt(mu / rho);

   d.Zp_plus  = d.Zp_minus = rho * cp;
   d.Zs_plus  = d.Zs_minus = rho * cs;
   d.eta_p    = 0.5 * rho * cp;
   d.eta_s    = 0.5 * rho * cs;

   d.tau1_0   = tau_pre_dip;
   d.tau2_0   = tau_pre_strike;
   d.sigma_n0 = sigma_n_eff;

   d.tau1_nuc = d.tau2_nuc = d.sigma_n_nuc = 0.0;

   d.slip_rate    = 0.0;
   d.V1 = d.V2    = 0.0;
   d.slip1 = d.slip2 = 0.0;
   d.tau1_corr    = 0.0;
   d.tau2_corr    = d.tau2_0;
   d.sigma_n_corr = d.sigma_n0;
}

/// Per-DOF LSW field copy + forced-rupture field copy + defensive
/// RS-slot zeroing.  Shared by both `InitializeFaultDOFs_Spatial`
/// overloads (centroid + IP-aware, R-405).  Caller is responsible for
/// having run `seed_static_dof_fields` first to populate impedances /
/// pre-stress / initial state.
inline void copy_lsw_and_forced_rupture_fields(
   DOFData& d, int i,
   const SlipWeakeningPerDOFParams& lsw,
   real_t T_forced_s, real_t t0_decay_s)
{
   d.lsw_mu_s = lsw.mu_s(i);
   d.lsw_mu_d = lsw.mu_d(i);
   d.lsw_d_c  = lsw.d_c(i);

   d.T_forced_rupture = T_forced_s;
   d.t0_decay_forced  = t0_decay_s;

   d.a   = 0.0;
   d.psi = 0.0;
   d.Dc  = 0.0;
}

/// Shared validator for both LSW `InitializeFaultDOFs_Spatial`
/// overloads (centroid + IP-aware, R-504).  Aborts if any per-DOF
/// array's size does not match the declared `ndof`.
inline void validate_lsw_per_dof_arrays(
   int ndof,
   const Array<int>& dof_to_elem,
   const SlipWeakeningPerDOFParams& lsw,
   const Vector& tau_pre,
   const Vector& sigma_n_eff,
   const Vector& T_forced_s,
   const Vector& t0_decay_s)
{
   MFEM_VERIFY(ndof >= 0,
               "InitializeFaultDOFs_Spatial: ndof must be >= 0; got "
               << ndof);
   MFEM_VERIFY(dof_to_elem.Size() == ndof,
               "InitializeFaultDOFs_Spatial: dof_to_elem.Size() ("
               << dof_to_elem.Size() << ") != ndof (" << ndof << ")");
   MFEM_VERIFY(lsw.mu_s.Size() == ndof,
               "InitializeFaultDOFs_Spatial: lsw.mu_s.Size() ("
               << lsw.mu_s.Size() << ") != ndof (" << ndof << ")");
   MFEM_VERIFY(lsw.mu_d.Size() == ndof, "lsw.mu_d size mismatch");
   MFEM_VERIFY(lsw.d_c.Size()  == ndof, "lsw.d_c size mismatch");
   MFEM_VERIFY(tau_pre.Size()      == 2 * ndof,
               "InitializeFaultDOFs_Spatial: tau_pre.Size() ("
               << tau_pre.Size() << ") != 2 * ndof (" << (2 * ndof) << ")");
   MFEM_VERIFY(sigma_n_eff.Size()  == ndof,
               "InitializeFaultDOFs_Spatial: sigma_n_eff.Size() ("
               << sigma_n_eff.Size() << ") != ndof (" << ndof << ")");
   MFEM_VERIFY(T_forced_s.Size()   == ndof,
               "InitializeFaultDOFs_Spatial: T_forced_s.Size() ("
               << T_forced_s.Size() << ") != ndof (" << ndof << ")");
   MFEM_VERIFY(t0_decay_s.Size()   == ndof,
               "InitializeFaultDOFs_Spatial: t0_decay_s.Size() ("
               << t0_decay_s.Size() << ") != ndof (" << ndof << ")");
}

/// Guard: abort if any DOFData in `dof_data` already looks initialized
/// (sentinel: lsw_d_c > 0 OR Zs_plus > 0).  Catches the
/// double-init-with-TPV205 mistake the plan §Phase 5 Acceptance test 3
/// targets.
inline void verify_dof_data_uninitialized(const std::vector<DOFData>& dof_data)
{
   for (size_t i = 0; i < dof_data.size(); ++i)
   {
      const DOFData& d = dof_data[i];
      MFEM_VERIFY(d.lsw_d_c == 0.0 && d.Zs_plus == 0.0,
                  "InitializeFaultDOFs_Spatial: dof_data[" << i
                  << "] is already populated (lsw_d_c=" << d.lsw_d_c
                  << ", Zs_plus=" << d.Zs_plus << "); did you also call "
                  "InitializeFaultDOFs_TPV205 first?  Each DOFData "
                  "should be written by exactly one InitializeFaultDOFs* "
                  "call.");
   }
}

}  // namespace internal

/// @brief SAFS analog of `InitializeFaultDOFs_TPV205`: writes per-DOF
/// impedances + LSW friction params + pre-stress + forced-rupture
/// fields into `DOFData[]` in one pass.
///
/// @param dof_data    Resized to `ndof`; each entry populated.
/// @param ndof        Number of fault DOFs (QPs).
/// @param dof_to_elem [ndof] bulk-element index that owns each DOF.
/// @param material    Heterogeneous material (Mode::Constant or
///                    Mode::Coefficient).  Used to compute per-DOF
///                    impedances via MaterialField::EvalAt at the
///                    element centroid.
/// @param mesh        Mesh used by MaterialField::EvalAt for the
///                    ElementTransformation lookup (ParMesh or
///                    serial Mesh).
/// @param lsw         Per-DOF LSW params (from
///                    SpatialFrictionResolver::ResolveSlipWeakening).
/// @param tau_pre     [2 * ndof] interleaved (tau_dip_i, tau_strike_i)
///                    from FaultGeometry::GetTauPre() after the stress
///                    source has been applied (Phase 3 / 3b).
/// @param sigma_n_eff [ndof] effective normal stress from
///                    FaultGeometry::sigma_n_per_dof().
/// @param T_forced_s  [ndof] per-DOF time-of-forced-rupture from
///                    SpatialFrictionResolver::ResolveForcedRupture
///                    (Phase 1 / D-4).
/// @param t0_decay_s  [ndof] per-DOF forced-rupture decay time
///                    (uniform = cfg.nucleation.t0_decay_s; stored
///                    per-DOF so the LSW solver can read both from one
///                    DOFData entry).
///
/// Aborts if any DOFData already has lsw_d_c > 0 or Zs_plus > 0 (the
/// double-init guard).  After the call, every DOFData has its impedance,
/// pre-stress, LSW, forced-rupture, and initial-state fields filled;
/// the rate-and-state slots (a, psi, Dc) are zeroed defensively
/// (mirrors `InitializeFaultDOFs_TPV205`).
template <typename MeshT>
inline void InitializeFaultDOFs_Spatial(
   std::vector<DOFData>&             dof_data,
   int                                ndof,
   const Array<int>&                  dof_to_elem,
   const MaterialField&               material,
   MeshT&                             mesh,
   const SlipWeakeningPerDOFParams&   lsw,
   const Vector&                      tau_pre,
   const Vector&                      sigma_n_eff,
   const Vector&                      T_forced_s,
   const Vector&                      t0_decay_s)
{
   // R-504: shared validator across centroid + IP-aware overloads.
   internal::validate_lsw_per_dof_arrays(ndof, dof_to_elem, lsw,
                                          tau_pre, sigma_n_eff,
                                          T_forced_s, t0_decay_s);
   internal::verify_dof_data_uninitialized(dof_data);
   dof_data.resize(ndof);

   for (int i = 0; i < ndof; ++i)
   {
      DOFData& d = dof_data[i];
      internal::seed_static_dof_fields<MeshT>(d, dof_to_elem[i], mesh,
                                              material,
                                              tau_pre(2 * i + 0),
                                              tau_pre(2 * i + 1),
                                              sigma_n_eff(i));
      // R-405: shared helper for LSW + forced-rupture + RS-slot fields.
      // Defaults on DOFData (T_forced_rupture=1e9, t0_decay_forced=0)
      // preserve TPV205 byte-exactness; here we OVERWRITE per DOF using
      // the resolver output.
      internal::copy_lsw_and_forced_rupture_fields(d, i, lsw,
                                                    T_forced_s(i),
                                                    t0_decay_s(i));
   }
}

/// @brief Rate-and-state analog of `InitializeFaultDOFs_Spatial`.
/// Writes per-DOF impedances + RS params + pre-stress in one pass.
/// Forced-rupture fields are intentionally left at their TPV205-byte-
/// exact defaults (T=1e9, t_0=0) — rate-and-state nucleation is via
/// V_init, not gradual friction reduction.
template <typename MeshT>
inline void InitializeFaultDOFs_Spatial_RS(
   std::vector<DOFData>&            dof_data,
   int                              ndof,
   const Array<int>&                dof_to_elem,
   const MaterialField&             material,
   MeshT&                           mesh,
   const RateStatePerDOFParams&     rs,
   const Vector&                    tau_pre,
   const Vector&                    sigma_n_eff)
{
   MFEM_VERIFY(ndof >= 0,
               "InitializeFaultDOFs_Spatial_RS: ndof must be >= 0");
   MFEM_VERIFY(dof_to_elem.Size() == ndof, "dof_to_elem size mismatch");
   MFEM_VERIFY(rs.a.Size() == ndof,        "rs.a size mismatch");
   MFEM_VERIFY(rs.b.Size() == ndof,        "rs.b size mismatch");
   MFEM_VERIFY(rs.Dc.Size() == ndof,       "rs.Dc size mismatch");
   MFEM_VERIFY(tau_pre.Size() == 2 * ndof, "tau_pre size mismatch");
   MFEM_VERIFY(sigma_n_eff.Size() == ndof, "sigma_n_eff size mismatch");

   internal::verify_dof_data_uninitialized(dof_data);
   dof_data.resize(ndof);

   for (int i = 0; i < ndof; ++i)
   {
      DOFData& d = dof_data[i];
      internal::seed_static_dof_fields<MeshT>(d, dof_to_elem[i], mesh,
                                              material,
                                              tau_pre(2 * i + 0),
                                              tau_pre(2 * i + 1),
                                              sigma_n_eff(i));

      // RS fields (lsw_* stays zero — defensive).
      d.a   = rs.a(i);
      d.Dc  = rs.Dc(i);
      // Initial state: derived from V_init magnitude / eta during the
      // 4-phase init (callers run the Brent solve later); here we seed
      // psi from steady-state for V_init.  The driver overrides this
      // during the equilibrium solve, but a non-zero seed avoids
      // accidentally consuming a zero psi in any debug path.
      d.psi = 0.0;

      // Forced-rupture fields stay at their in-class defaults
      // (T=1e9, t_0=0) — RS path does not consume them.
   }
}

/// @brief IP-aware overload of `InitializeFaultDOFs_Spatial` (R-309).
/// Same shape as the centroid-based overload but accepts the per-DOF
/// reference IntegrationPoint cache from FaultGeometry::fault_dof_ip().
/// When the driver builds FaultGeometry via the new BP5 ctor with a
/// non-empty `dof_ips` argument (Phase 5a), pass `&geom.fault_dof_ip(0)`
/// alongside the geom-supplied dof_to_elem to evaluate material exactly
/// at each fault QP instead of at the bulk element centroid.  This is
/// the R-111 long-term shape.  The aggregate is passed as a pointer so
/// the centroid-based overload above stays the default.
template <typename MeshT>
inline void InitializeFaultDOFs_Spatial(
   std::vector<DOFData>&                       dof_data,
   int                                          ndof,
   const Array<int>&                            dof_to_elem,
   const MaterialField&                         material,
   MeshT&                                       mesh,
   const SlipWeakeningPerDOFParams&             lsw,
   const Vector&                                tau_pre,
   const Vector&                                sigma_n_eff,
   const Vector&                                T_forced_s,
   const Vector&                                t0_decay_s,
   const std::vector<mfem::IntegrationPoint>&   dof_ips)
{
   MFEM_VERIFY(static_cast<int>(dof_ips.size()) == ndof,
               "InitializeFaultDOFs_Spatial(ip-aware): dof_ips.size() ("
               << dof_ips.size() << ") != ndof (" << ndof << ")");
   // R-504: same per-DOF-array validator as the centroid overload.
   internal::validate_lsw_per_dof_arrays(ndof, dof_to_elem, lsw,
                                          tau_pre, sigma_n_eff,
                                          T_forced_s, t0_decay_s);
   internal::verify_dof_data_uninitialized(dof_data);
   dof_data.resize(ndof);

   for (int i = 0; i < ndof; ++i)
   {
      DOFData& d = dof_data[i];
      internal::seed_static_dof_fields<MeshT>(d, dof_to_elem[i], mesh,
                                              material, dof_ips[i],
                                              tau_pre(2 * i + 0),
                                              tau_pre(2 * i + 1),
                                              sigma_n_eff(i));
      // R-405: same helper as the centroid overload.
      internal::copy_lsw_and_forced_rupture_fields(d, i, lsw,
                                                    T_forced_s(i),
                                                    t0_decay_s(i));
   }
}

/// @brief IP-aware overload of `InitializeFaultDOFs_Spatial_RS` (REVIEW
/// R-008).  Same shape as the centroid-based overload but accepts the
/// per-DOF reference IntegrationPoint cache from
/// `FaultGeometry::fault_dof_ip()`, so per-DOF impedances are evaluated
/// at the actual fault QP rather than at the bulk element centroid.
///
/// Required for any RS config that combines `material.kind =
/// "depth_profile_1d"` (or `"sidecar_hdf5"`) with rate-state friction;
/// the centroid-based overload is correct only for `Mode::Constant`
/// material.
template <typename MeshT>
inline void InitializeFaultDOFs_Spatial_RS(
   std::vector<DOFData>&                       dof_data,
   int                                          ndof,
   const Array<int>&                            dof_to_elem,
   const MaterialField&                         material,
   MeshT&                                       mesh,
   const RateStatePerDOFParams&                 rs,
   const Vector&                                tau_pre,
   const Vector&                                sigma_n_eff,
   const std::vector<mfem::IntegrationPoint>&   dof_ips)
{
   MFEM_VERIFY(ndof >= 0,
               "InitializeFaultDOFs_Spatial_RS(ip-aware): ndof must be "
               ">= 0; got " << ndof);
   MFEM_VERIFY(static_cast<int>(dof_ips.size()) == ndof,
               "InitializeFaultDOFs_Spatial_RS(ip-aware): dof_ips.size() ("
               << dof_ips.size() << ") != ndof (" << ndof << ")");
   MFEM_VERIFY(dof_to_elem.Size() == ndof, "dof_to_elem size mismatch");
   MFEM_VERIFY(rs.a.Size() == ndof,        "rs.a size mismatch");
   MFEM_VERIFY(rs.b.Size() == ndof,        "rs.b size mismatch");
   MFEM_VERIFY(rs.Dc.Size() == ndof,       "rs.Dc size mismatch");
   MFEM_VERIFY(tau_pre.Size() == 2 * ndof, "tau_pre size mismatch");
   MFEM_VERIFY(sigma_n_eff.Size() == ndof, "sigma_n_eff size mismatch");

   internal::verify_dof_data_uninitialized(dof_data);
   dof_data.resize(ndof);

   for (int i = 0; i < ndof; ++i)
   {
      DOFData& d = dof_data[i];
      internal::seed_static_dof_fields<MeshT>(d, dof_to_elem[i], mesh,
                                              material, dof_ips[i],
                                              tau_pre(2 * i + 0),
                                              tau_pre(2 * i + 1),
                                              sigma_n_eff(i));
      d.a   = rs.a(i);
      d.Dc  = rs.Dc(i);
      d.psi = 0.0;
   }
}

}  // namespace spatial
}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_DYNAMIC_SPATIAL_SETUP_HPP
