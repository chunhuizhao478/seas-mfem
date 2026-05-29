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
// R-001: DieterichRuinaFriction (Constants / InitialStatePsi) is used by
// SeedEquilibriumPsi_RS below and is not pulled in by any of the includes
// above.  FrictionSolver (R-009 V0 guard) is visible transitively via
// fault_face_flux.hpp -> friction_solver.hpp.
#include "../friction/dieterich_ruina.hpp"

#include <cmath>
#include <cstddef>
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
   d.lsw_cohesion = lsw.cohesion(i);   // Phase 10 (TPV31): per-DOF C0 [Pa]

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
      d.b   = rs.b(i);   // Phase 11a: per-DOF state-evolution b (was scalar blk.b_default)
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

/// @brief Seed the per-DOF rate-and-state ψ to equilibrium for the
/// resolved pre-stress (overwrites the `d.psi = 0.0` stub left by
/// `InitializeFaultDOFs_Spatial_RS`).
///
/// For each DOF, ψ is chosen so the fault is in steady sliding at
/// `rs.V_init(i)` under the current total pre-stress
///   tau0 = |(tau1_0 + tau1_nuc, tau2_0 + tau2_nuc)|
/// i.e. it inverts the regularized friction balance
///   tau0 = σ_n0 · a · asinh[V_init/(2·V0) · exp(ψ/a)] + η · V_init
/// via `DieterichRuinaFriction::InitialStatePsi`.  When the effective
/// driving traction `tau_eff = tau0 − η·V_init` is non-positive (no
/// pre-stress shear, or radiation damping `η·V_init` that exceeds tau0),
/// `InitialStatePsi`'s `log[(2V0/V_init)·sinh(f/a)]` would go to inf/NaN,
/// so the locked steady-state ψ = f0 + b·ln(V0/V_init) is used instead
/// (R-014).
///
/// @param dof_data  Per-DOF state (already seeded by
///                  `InitializeFaultDOFs_Spatial_RS`); only `psi` is
///                  written.
/// @param rs        Per-DOF RS params: a, b, Dc, V_init are read (Phase 11a:
///                  b is now per-DOF `rs.b(ii)`, not the scalar blk.b_default).
/// @param blk       Global RS scalars: V_0, f_0 are read (b moved to `rs`).
///
/// R-001 (CRITICAL): the radiation-damping argument is `d.eta_s` (a
/// `DOFData` member set by `seed_static_dof_fields` at :82 == 0.5·sqrt(μρ),
/// exactly what `FaultFaceFlux` uses in the friction solve), NOT
/// `rs.eta_s` — `RateStatePerDOFParams` has only a per-DOF `eta` Vector
/// (spatial_friction.hpp:295) and no `eta_s` member, so `rs.eta_s` would
/// be a hard compile error.  For the homogeneous SAFS material
/// `d.eta_s == rs.eta(i)` at the same element centroid, so the seed's
/// damping matches the friction solve's.
///
/// R-009 (V0 single source of truth): the leading guard asserts the
/// config `blk.V_0_default` equals the solver's compile-time
/// `FrictionSolver::V0` (= 1e-6).  The force solve hard-codes
/// `FrictionSolver::V0` (fault_face_flux.cpp:224) while this seed and the
/// `AgingLawPsi` ψ-update use `blk.V_0_default`; if they disagree the
/// fault is NOT at V_init at t=0 and the friction is wrong throughout
/// with no error.  The deliverable uses the default V_0 = 1e-6 (the guard
/// passes); it is a tripwire for future tuning.  (Proper fix — thread the
/// resolved V0 into FrictionSolver — is a no-hardcoded-numbers follow-up
/// gated out of the byte-exact TPV oracle path.)
inline void SeedEquilibriumPsi_RS(
   std::vector<DOFData>&            dof_data,
   const RateStatePerDOFParams&     rs,   // a, b, Dc, V_init per-DOF
   const RateStateBlock&            blk)  // V_0, f_0 scalar globals (b now per-DOF in rs)
{
   // R-012: the loop reads rs.{V_init,Dc,a}(ii) for ii in [0, nd); mfem::
   // Vector::operator() is unchecked in release, so a too-short rs is an
   // OOB read.  Validate up front like InitializeFaultDOFs_Spatial_RS.
   const int nd = static_cast<int>(dof_data.size());
   MFEM_VERIFY(rs.a.Size() >= nd && rs.b.Size() >= nd && rs.Dc.Size() >= nd
               && rs.V_init.Size() >= nd,
               "SeedEquilibriumPsi_RS: rs vectors (a=" << rs.a.Size()
               << ", b=" << rs.b.Size() << ", Dc=" << rs.Dc.Size()
               << ", V_init=" << rs.V_init.Size()
               << ") must cover dof_data.size()=" << nd
               << "; call after InitializeFaultDOFs_Spatial_RS with the same rs.");

   // R-009: V_0 has two consumers — the config (this seed + AgingLawPsi)
   // and the hardcoded FrictionSolver::V0 used by the force solve.  They
   // MUST agree or the fault is not in equilibrium at t=0.
   MFEM_VERIFY(std::abs(blk.V_0_default - FrictionSolver::V0)
               <= 1e-30 + 1e-12 * FrictionSolver::V0,
               "SeedEquilibriumPsi_RS: [friction.rate_state].V_0 ("
               << blk.V_0_default << ") must equal FrictionSolver::V0 ("
               << FrictionSolver::V0 << "); the force solve hardcodes V0.");

   for (std::size_t i = 0; i < dof_data.size(); ++i)
   {
      const int ii = static_cast<int>(i);  // mfem::Vector accessors take int
      DOFData& d = dof_data[i];

      // R-028: tripwire that InitializeFaultDOFs_*_RS set the per-DOF d.b
      // (Phase 11a) before seeding.  DOFData.b defaults to NaN; a forgotten
      // init would otherwise reach the aging-law UpdateStateAnalytic where a
      // NaN/0 b silently mis-evolves psi (b=0: psi snaps to f0 for psi<f0).
      MFEM_VERIFY(std::isfinite(d.b) && d.b > 0.0,
                  "SeedEquilibriumPsi_RS: DOFData.b not set (= " << d.b
                  << ") at DOF " << ii << "; InitializeFaultDOFs_*_RS must set "
                  "d.b > 0 before seeding (Phase 11a per-DOF b).");

      // R-001 edge guard: else log()/InitialStatePsi -> inf/NaN psi.
      MFEM_VERIFY(rs.V_init(ii) > 0.0,
                  "SeedEquilibriumPsi_RS: V_init <= 0 at DOF " << ii
                  << " (ill-posed steady state)");

      // R-017: InitialStatePsi computes f = (tau0 - eta*V_init)/sigma_n0;
      // sigma_n0 == 0 -> f = inf -> psi = inf, sigma_n0 < 0 -> log(negative)
      // -> NaN, and InitialStatePsi's only check (MFEM_ASSERT log_arg>0) is
      // debug-only and passes the +inf case.  ResolveRateState enforces
      // sigma_n_eff > 0 upstream (spatial_friction.cpp:1073), so this guard
      // is a defensive tripwire consistent with the V_init guard above.
      MFEM_VERIFY(d.sigma_n0 > 0.0,
                  "SeedEquilibriumPsi_RS: sigma_n0 <= 0 at DOF " << ii
                  << " (" << d.sigma_n0 << "); effective normal stress must "
                  "be compressive (>0) or InitialStatePsi seeds inf/NaN psi.");

      const DieterichRuinaFriction fr(
         DieterichRuinaFriction::Constants{
            blk.V_0_default, blk.f_0_default, rs.b(ii), rs.Dc(ii)});

      const real_t tau0 = std::hypot(d.tau1_0 + d.tau1_nuc,
                                     d.tau2_0 + d.tau2_nuc);

      // R-014: InitialStatePsi computes f = (tau0 - eta*V_init)/sigma_n
      // and log[(2V0/V_init)*sinh(f/a)]; if the radiation damping eta*V_init
      // exceeds the prestress (tau_eff <= 0) the log argument goes negative
      // -> NaN psi.  Branch on the EFFECTIVE traction, not tau0, so the
      // damping-dominated case falls back to the locked steady state.
      // R-001: damping is d.eta_s (DOFData member), NEVER rs.eta_s
      // (RateStatePerDOFParams has no such member).
      const real_t tau_eff = tau0 - d.eta_s * rs.V_init(ii);
      d.psi = (tau_eff > 0.0)
         ? fr.InitialStatePsi(tau0, rs.V_init(ii),
                              d.sigma_n0, d.eta_s, rs.a(ii))
         : (blk.f_0_default
            + rs.b(ii) * std::log(blk.V_0_default / rs.V_init(ii)));
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

}  // namespace spatial
}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_DYNAMIC_SPATIAL_SETUP_HPP
