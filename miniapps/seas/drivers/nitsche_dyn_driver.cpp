// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_dyn_driver.cpp — Nitsche dynamic rupture, plan item N-07 (build
// sequence rev 5.2, phase P7): the driver.
//
// Reads the SAFS spatial deck through spatial::LoadSpatialFrictionConfig (the
// SAFS driver's own parser) plus the Nitsche-only keys the SAFS parser does
// not know ([numerics] penalty_factor, the fault-coordinate specification's
// [fault_geometry] keys, [diagnostics] energy_balance), builds the mesh, the
// material, the boundary configuration and the fault setup, constructs the
// NitscheWaveOperator (N-01) and advances it with MFEM's
// CentralDifferenceSolver (N-02).  Nothing benchmark-specific lives here: the
// law comes from a factory keyed on [meta] law, the friction parameters from
// the SAFS resolver, the prestress and the nucleation from the deck as global
// tensors (spec I1), the station list from the [problem] tag.
//
// Outputs (in --output-dir):
//   <prefix>_station_<name>.dat   TPV102StationWriter (column layout of the
//                                 existing runs) through a DOFData adapter;
//   <prefix>_station_conventions.txt  frame, pin and sign map (I3);
//   <prefix>_fault_r<rank>.csv    fault-field dump every --dump-every steps;
//   <prefix>_rupture_time_r<rank>.csv  first |V| > 1 mm/s per fault DOF;
//   stdout log every --log-every steps.
//
// Station columns (SCEC map, plan "The SCEC conventions, mapped"):
//   slip1 = delta.dhat (v-slip), slip2 = delta.shat (h-slip),
//   V1 = V.dhat, V2 = V.shat,
//   tau1 = -T.dhat, tau2 = -T.shat (SCEC shear tractions, Pa),
//   sigma_n = -T.nhat (compression positive, Pa),
// with T the total traction (Nitsche traction at the end of the step plus
// prestress plus nucleation) on the face with normal nhat.
//
// Scope (Chunhui Zhao, 2026-10-14): tetrahedral meshes; constant material
// (the sidecar / depth-profile / halfspace selections are P8);
// rate-and-state (slip weakening comes with TPV205).

#include "mfem.hpp"
#include "../domain/nitsche_wave_operator.hpp"
#include "../domain/nitsche_energy.hpp"
#include "../solver/nitsche_cfl.hpp"
#include "../friction/rate_state_closure.hpp"
#include "../fault/fault_frame.hpp"
#include "../spatial/code/spatial_friction.hpp"
#include "../dynamic/spatial_nucleation.hpp"
#include "../dynamic/tpv102_stations.hpp"

#include <toml.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
namespace ff = mfem::seas::fault_frame;

namespace
{

/// P3's measured stable constants C(p), p = 1..4 (dt = cfl C(p) h/(p^2 c_p)).
const real_t kCeig[5] = {0.0, 0.113, 0.212, 0.282, 0.332};

int g_rank = 0;

std::string ArgStr(int argc, char *argv[], const char *flag, const std::string &def)
{
   for (int i = 1; i + 1 < argc; i++) { if (std::string(argv[i]) == flag) { return argv[i + 1]; } }
   return def;
}
real_t ArgReal(int argc, char *argv[], const char *flag, real_t def)
{
   const std::string v = ArgStr(argc, argv, flag, "");
   return v.empty() ? def : std::atof(v.c_str());
}
int ArgInt(int argc, char *argv[], const char *flag, int def)
{
   const std::string v = ArgStr(argc, argv, flag, "");
   return v.empty() ? def : std::atoi(v.c_str());
}

const toml::value *DeckTable(const toml::value &root, const char *name)
{
   return root.contains(name) ? &root.at(name) : nullptr;
}
real_t TomlReal(const toml::value *t, const char *key, real_t def)
{
   if (!t || !t->contains(key)) { return def; }
   const auto &v = t->at(key);
   if (v.is_floating()) { return static_cast<real_t>(v.as_floating()); }
   if (v.is_integer()) { return static_cast<real_t>(v.as_integer()); }
   MFEM_ABORT("deck key '" << key << "' must be a number");
   return def;
}
std::string TomlStr(const toml::value *t, const char *key, const std::string &def)
{
   if (!t || !t->contains(key)) { return def; }
   return t->at(key).as_string();
}
bool TomlBool(const toml::value *t, const char *key, bool def)
{
   if (!t || !t->contains(key)) { return def; }
   return t->at(key).as_boolean();
}

/// Friction-law factory keyed on [meta] law (plan N-07, Generality).
struct ClosureBundle
{
   std::unique_ptr<DieterichRuinaFriction> fr;
   std::unique_ptr<IFaultClosure> closure;
   std::string name;
};

ClosureBundle MakeClosure(const spatial::SpatialFrictionConfig &cfg)
{
   ClosureBundle B;
   if (cfg.law == spatial::FrictionLawKind::RateState)
   {
      MFEM_VERIFY(cfg.rate_state.has_value(), "[meta] law = rate_state needs "
                  "[friction.rate_state]");
      const auto &rs = *cfg.rate_state;
      B.fr = std::make_unique<DieterichRuinaFriction>(
         DieterichRuinaFriction::Constants{rs.V_0_default, rs.f_0_default,
                                        rs.b_default, rs.Dc_default});
      RateStateClosure::StateLaw law = RateStateClosure::StateLaw::Aging;
      if (rs.state_evolution == spatial::StateEvolutionKind::SlipLawStrongRateWeakening)
      {
         MFEM_ABORT("Nitsche driver: state_evolution = slip-law strong rate "
                    "weakening (TPV104) has no IFaultClosure yet.");
      }
      B.closure = std::make_unique<RateStateClosure>(*B.fr, law);
      B.name = "rate_state (aging law)";
      return B;
   }
   MFEM_ABORT("Nitsche driver: [meta] law = slip_weakening has no "
              "IFaultClosure yet (it comes with TPV205).");
   return B;
}

/// Unit vector helpers.
void Normalize(real_t v[3])
{
   const real_t l = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
   MFEM_VERIFY(l > 0.0, "zero-length vector in [fault_geometry]");
   for (int i = 0; i < 3; i++) { v[i] /= l; }
}

void MakeDirs(const std::string &d)
{
   if (g_rank == 0) { ::mkdir(d.c_str(), 0755); }
   MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   const int nranks = Mpi::WorldSize();
   auto wall0 = std::chrono::steady_clock::now();

   const std::string cfg_path = ArgStr(argc, argv, "--config", "");
   if (cfg_path.empty())
   {
      if (g_rank == 0)
      {
         std::cerr << "usage: seas_nitsche_dyn_driver --config deck.toml [--mesh m.msh]"
                      " [--order p] [--tfinal s] [--output-dir d] [--prefix name]"
                      " [--station-every n] [--dump-every n] [--log-every n]\n";
      }
      return 1;
   }

   // ------------------------------------------------------------ deck
   spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(cfg_path);
   const toml::value root = toml::parse(cfg_path);
   const toml::value *t_num = DeckTable(root, "numerics");
   const toml::value *t_fg = DeckTable(root, "fault_geometry");
   const toml::value *t_diag = DeckTable(root, "diagnostics");

   const std::string mesh_path = ArgStr(argc, argv, "--mesh", cfg.mesh.path);
   const int p = ArgInt(argc, argv, "--order", cfg.mesh.order);
   const real_t tfinal = ArgReal(argc, argv, "--tfinal", cfg.time.tfinal);
   const std::string out_dir = ArgStr(argc, argv, "--output-dir", cfg.output.output_dir);
   const std::string prefix = ArgStr(argc, argv, "--prefix",
                                     cfg.problem.tag.empty() ? "nitsche" : cfg.problem.tag);
   const int station_every = ArgInt(argc, argv, "--station-every", 1);
   const int dump_every = ArgInt(argc, argv, "--dump-every", 0);
   const int log_every = ArgInt(argc, argv, "--log-every", 100);
   const real_t cfl = cfg.numerics.cfl;
   const real_t penalty_factor = TomlReal(t_num, "penalty_factor", 1.0);
   const bool energy_on = TomlBool(t_diag, "energy_balance", false);
   MFEM_VERIFY(p >= 1 && p <= 4, "order p must be in 1..4 (P3 constants); got " << p);
   MFEM_VERIFY(cfl > 0.0 && cfl <= 1.0, "[numerics] cfl (fraction of C(p)) must be in (0,1]");
   MakeDirs(out_dir);

   if (g_rank == 0)
   {
      std::cout << "=== seas_nitsche_dyn_driver (plan N-07) on " << nranks << " rank(s) ===\n"
                << "  deck " << cfg_path << "\n  " << cfg.description << "\n";
      // velocity-stress numerics keys: listed once, ignored
      const char *vs_keys[] = {"ader_order", "time_integrator", "mixed_flux",
                               "interior_flux", "fault_iterator", "use_pml",
                               "cfl_safety", "cfl_dg_safety", "lts"};
      std::string ign;
      for (const char *k : vs_keys) { if (t_num && t_num->contains(k)) { ign += std::string(" ") + k; } }
      if (!ign.empty()) { std::cout << "  ignored [numerics] keys (velocity-stress path):" << ign << "\n"; }
   }

   // ------------------------------------------------------------ mesh
   std::unique_ptr<Mesh> smesh = std::make_unique<Mesh>(mesh_path.c_str(), 1, 1);
   MFEM_VERIFY(smesh->Dimension() == 3, "mesh must be 3-D");
   for (int e = 0; e < smesh->GetNE(); e++)
   {
      MFEM_VERIFY(smesh->GetElementBaseGeometry(e) == Geometry::TETRAHEDRON,
                  "Nitsche driver: tetrahedral meshes only; element " << e << " is "
                  << Geometry::Name[smesh->GetElementBaseGeometry(e)]);
   }
   const long long ne_global = smesh->GetNE();
   ParMesh pmesh(MPI_COMM_WORLD, *smesh);
   smesh.reset();

   // ------------------------------------------------------------ material
   MFEM_VERIFY(cfg.material.kind == spatial::MaterialKind::Constant && !cfg.velocity.use_sidecar,
               "Nitsche driver (P7): only the constant material "
               "([velocity] use_sidecar = false, [material] kind constant) is "
               "wired; the sidecar / depth-profile / halfspace selections are P8.");
   MaterialField material = MaterialField::MakeConstant(
      cfg.material_fallback.lambda, cfg.material_fallback.mu, cfg.material_fallback.rho);

   // ------------------------------------------------------------ boundary
   MFEM_VERIFY(cfg.boundary.fault_attr > 0, "[boundary] fault_attr is required");
   BoundaryConfig bc;
   bc.fault_attr = cfg.boundary.fault_attr;
   for (int a : cfg.boundary.natural_attrs) { bc.natural_attrs.insert(a); }
   for (int a : cfg.boundary.absorbing_attrs) { bc.absorbing_attrs.insert(a); }

   // ------------------------------------------------------------ fault geometry
   NitscheFaultSetup fs;
   {
      const std::string orient = TomlStr(t_fg, "orientation", "propagated");
      const std::string pin_by = TomlStr(t_fg, "pin_by", "reference");
      const std::string tangent = TomlStr(t_fg, "tangent", "pointwise");
      const std::string junction = TomlStr(t_fg, "junction_policy", "refuse");
      MFEM_VERIFY(orient == "propagated", "[fault_geometry] orientation = '" << orient
                  << "': only \"propagated\" is implemented (N-08)");
      MFEM_VERIFY(tangent == "pointwise", "[fault_geometry] tangent = '" << tangent
                  << "': only \"pointwise\" is implemented (N-09; F8 transport is not)");
      MFEM_VERIFY(junction == "refuse", "[fault_geometry] junction_policy = '" << junction
                  << "': only \"refuse\" is implemented");
      auto &oc = fs.orientation;
      for (int i = 0; i < 3; i++)
      {
         oc.ref_normal[i] = cfg.fault_geometry.ref_normal[i];
         oc.up[i] = cfg.fault_geometry.up[i];
         fs.up[i] = cfg.fault_geometry.up[i];
      }
      if (pin_by == "reference") { oc.pin = FaultOrientationConfig::Pin::Reference; }
      else
      {
         MFEM_ABORT("[fault_geometry] pin_by = '" << pin_by << "': only \"reference\" "
                    "is wired in the driver (Pin::Point needs a pin_point key)");
      }
      oc.c_min = TomlReal(t_fg, "pin_min_margin", oc.c_min);
      oc.fault_attrs = {bc.fault_attr};
   }

   // ------------------------------------------------------------ closure
   ClosureBundle law = MakeClosure(cfg);
   fs.closure = law.closure.get();
   fs.sigma_min = cfg.stress.pore_pressure.min_sigma_n_pa;
   MFEM_VERIFY(cfg.stress.pore_pressure.P_p_pa == 0.0
               && cfg.stress.pore_pressure.P_p_grad_pa_per_m == 0.0,
               "Nitsche driver (P7): pore pressure is not wired (P_p must be 0)");

   // ------------------------------------------------------------ prestress (I1)
   MFEM_VERIFY(cfg.stress.kind == spatial::StressSourceKind::ConstantTensor,
               "Nitsche driver: [stress] kind must be \"constant_tensor\" (the global "
               "Cauchy tensor, tension positive; spec I1).  fault_local_prestress is "
               "the slot pattern I2 forbids.");
   const real_t s0[9] = {cfg.stress.sigma_xx_pa, cfg.stress.sigma_xy_pa, cfg.stress.sigma_xz_pa,
                         cfg.stress.sigma_xy_pa, cfg.stress.sigma_yy_pa, cfg.stress.sigma_yz_pa,
                         cfg.stress.sigma_xz_pa, cfg.stress.sigma_yz_pa, cfg.stress.sigma_zz_pa};
   fs.prestress = [s0](const real_t[3], real_t s[9]) { for (int i = 0; i < 9; i++) { s[i] = s0[i]; } };

   // pin frame (F1 with the reference normal): nhat, shat, dhat
   real_t pn[3], pu[3], ps[3], pd[3];
   for (int i = 0; i < 3; i++) { pn[i] = cfg.fault_geometry.ref_normal[i]; pu[i] = cfg.fault_geometry.up[i]; }
   Normalize(pn);
   ff::Cross(pu, pn, ps);
   Normalize(ps);
   ff::Cross(ps, pn, pd);

   // ------------------------------------------------------------ nucleation
   if (cfg.nucleation.enabled)
   {
      MFEM_VERIFY(cfg.nucleation.kind == spatial::NucleationKind::GradualOverstressCompactCircular,
                  "Nitsche driver (P7): only [nucleation] kind = "
                  "\"gradual_overstress_compact_circular\" is wired");
      const auto spec = cfg.nucleation.compact_circular;
      MFEM_VERIFY(spec.radius_m > 0.0 && spec.T_nuc_s > 0.0, "compact-circular nucleation: "
                  "radius_m and T_nuc_s must be > 0");
      // strike overstress delta_tau F(r) G(t) as the global tensor
      //   dsigma = -delta_tau F G (shat (x) nhat + nhat (x) shat),
      // so that -(dsigma nhat).shat = +delta_tau F G.  r is the in-plane
      // distance from the centre measured with (dhat, shat), as the SAFS
      // resolver ResolveGradualOverstressCompactCircular measures it.
      std::array<real_t, 3> c = {spec.center_x_m, spec.center_y_m, spec.center_z_m};
      std::array<real_t, 3> sv = {ps[0], ps[1], ps[2]}, nv = {pn[0], pn[1], pn[2]},
                            dv = {pd[0], pd[1], pd[2]};
      fs.nucleation_spatial = [spec, c, sv, nv, dv](const real_t x[3], real_t s[9]) {
         real_t e[3] = {x[0] - c[0], x[1] - c[1], x[2] - c[2]};
         const real_t ed = e[0] * dv[0] + e[1] * dv[1] + e[2] * dv[2];
         const real_t es = e[0] * sv[0] + e[1] * sv[1] + e[2] * sv[2];
         const real_t amp = -spec.delta_tau_pa
                            * spatial::CompactBellFactor(std::sqrt(ed * ed + es * es), spec.radius_m);
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) { s[3 * i + j] = amp * (sv[i] * nv[j] + nv[i] * sv[j]); }
      };
      const real_t T = spec.T_nuc_s;
      fs.nucleation_temporal = [T](real_t t) { return spatial::SmoothStep(t, T); };
   }

   // ------------------------------------------------------------ friction parameters
   spatial::RateStateBlock rs = *cfg.rate_state;
   rs.eta_auto = false;      // eta is quasi-dynamic radiation damping; the
   rs.eta_default = 0.0;     // Nitsche closure uses gamma = 1/2 sigma_F dt
   spatial::SpatialFrictionResolver resolver;
   const real_t b0 = rs.b_default, f00 = rs.f_0_default, V00 = rs.V_0_default;
   fs.parameters = [&](const real_t x[3], FaultPointState &st) {
      Vector xc(3);
      for (int i = 0; i < 3; i++) { xc(i) = x[i]; }
      Array<int> el(1), at(1);
      el[0] = 0;
      at[0] = bc.fault_attr;
      Vector sn0;
      spatial::RateStatePerDOFParams r = resolver.ResolveRateState(
         rs, xc, el, at, material, pmesh, cfg.stress.pore_pressure, sn0);
      MFEM_VERIFY(r.b(0) == b0 && r.f_0(0) == f00 && r.V_0(0) == V00,
                  "Nitsche driver: per-DOF b / f_0 / V_0 differ from the defaults at ("
                  << x[0] << ", " << x[1] << ", " << x[2] << "); the closure takes "
                  "them as constants (spatially varying b is not wired)");
      st.a = r.a(0);
      st.Dc = r.Dc(0);
   };

   // ------------------------------------------------------------ operator
   DomainConfig dcfg;
   dcfg.penalty_factor = penalty_factor;
   NitscheWaveOperator op(pmesh, p, material, bc, dcfg, &fs);
   MFEM_VERIFY(op.Orientation().GetReport().junction_edges == 0,
               "[fault_geometry] junction_policy = \"refuse\": the fault has "
               << op.Orientation().GetReport().junction_edges << " junction edges");
   const real_t dt_cfl = ComputeMaxDt(op.Sample(), pmesh, p, cfl * kCeig[p]);
   const real_t dt = std::min(dt_cfl, cfg.time.dt_max);
   op.SetTimeStep(dt);

   // initial fault state: d = 0, psi from the continuous friction condition
   // at the projected prestress, V = V_init along -X_tan (C-07)
   const int no = op.NumOwnedFaultDOFs();
   {
      Vector d(3 * no), V(3 * no), psi(no), S0;
      d = 0.0;
      op.BackgroundTraction(cfg.time.t_initial, S0);
      Vector xc(3);
      Array<int> el(1), at(1);
      el[0] = 0; at[0] = bc.fault_attr;
      for (int j = 0; j < no; j++)
      {
         const real_t *n = &op.FaultFrames()(9 * j), *dd = &op.FaultFrames()(9 * j + 3),
                       *ss = &op.FaultFrames()(9 * j + 6);
         real_t X[2], sig;
         ff::ToLocal(&S0(3 * j), 0.0, n, dd, ss, X[0], X[1], sig);
         const real_t xm = std::hypot(X[0], X[1]);
         for (int c = 0; c < 3; c++) { xc(c) = op.FaultPositions()(3 * j + c); }
         Vector sn0;
         spatial::RateStatePerDOFParams r = resolver.ResolveRateState(
            rs, xc, el, at, material, pmesh, cfg.stress.pore_pressure, sn0);
         const real_t Vi = r.V_init(0);
         FaultPointState st;
         st.a = op.FaultA()(j);
         st.Dc = op.FaultDc()(j);
         law.closure->SeedEquilibrium(xm, Vi, sig, st);
         psi(j) = st.psi;
         real_t Vl[2] = {0.0, 0.0};
         if (xm > 0.0) { Vl[0] = -Vi * X[0] / xm; Vl[1] = -Vi * X[1] / xm; }
         ff::Embed(Vl[0], Vl[1], dd, ss, &V(3 * j));
      }
      op.SetFaultState(d, V, psi);
      // initial fault state per owned DOF (generic; e.g. to check psi^0
      // against a benchmark's closed form outside the driver)
      std::ofstream ini(out_dir + "/" + prefix + "_initial_r" + std::to_string(g_rank) + ".csv");
      ini << std::setprecision(17) << "# x y z a Dc psi0 tau0_strike(-T.shat) sigma_n0 V0_strike\n";
      for (int j = 0; j < no; j++)
      {
         const real_t *n = &op.FaultFrames()(9 * j), *dd = &op.FaultFrames()(9 * j + 3),
                       *ss = &op.FaultFrames()(9 * j + 6);
         real_t X[2], sig;
         ff::ToLocal(&S0(3 * j), 0.0, n, dd, ss, X[0], X[1], sig);
         const real_t *x = &op.FaultPositions()(3 * j);
         ini << x[0] << " " << x[1] << " " << x[2] << " " << op.FaultA()(j) << " "
             << op.FaultDc()(j) << " " << psi(j) << " " << -X[1] << " " << sig << " "
             << ff::Dot(&V(3 * j), ss) << "\n";
      }
   }

   // ------------------------------------------------------------ nucleation check
   // The projected SCEC along-strike traction -(S.shat) at the fault DOF
   // nearest the deck's hypocentre against the analytic tau_0 + delta_tau
   // F(r) G(t) at that DOF's position (plan P7 decision 9).
   if (cfg.nucleation.enabled)
   {
      const auto &spec = cfg.nucleation.compact_circular;
      const real_t h[3] = {cfg.hypocenter.x_m, cfg.hypocenter.y_m, cfg.hypocenter.z_m};
      struct { double d; int r; } loc = {1e300, g_rank}, glob;
      int jh = -1;
      for (int j = 0; j < no; j++)
      {
         const real_t *x = &op.FaultPositions()(3 * j);
         const real_t dd = std::hypot(x[0] - h[0], x[1] - h[1], x[2] - h[2]);
         if (dd < loc.d) { loc.d = dd; jh = j; }
      }
      MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
      if (glob.r == g_rank && jh >= 0)
      {
         const real_t *x = &op.FaultPositions()(3 * jh);
         const real_t *ss = &op.FaultFrames()(9 * jh + 6), *dd = &op.FaultFrames()(9 * jh + 3);
         real_t e[3] = {x[0] - spec.center_x_m, x[1] - spec.center_y_m, x[2] - spec.center_z_m};
         const real_t r = std::hypot(ff::Dot(e, dd), ff::Dot(e, ss));
         const real_t F = spatial::CompactBellFactor(r, spec.radius_m);
         Vector S0;
         op.BackgroundTraction(0.0, S0);
         const real_t tau0 = -ff::Dot(&S0(3 * jh), ss);
         std::cout << std::setprecision(8) << "  nucleation check at the fault DOF nearest the hypocentre ("
                   << x[0] << ", " << x[1] << ", " << x[2] << "), " << glob.d << " m away, r = "
                   << r << " m, F(r) = " << F << ":\n     t [s]   projected -T.shat [MPa]   analytic "
                   "tau_0 + dtau F G [MPa]   difference [Pa]\n";
         for (real_t tc : {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5})
         {
            Vector S;
            op.BackgroundTraction(tc, S);
            const real_t proj = -ff::Dot(&S(3 * jh), ss);
            const real_t ana = tau0 + spec.delta_tau_pa * F * spatial::SmoothStep(tc, spec.T_nuc_s);
            std::cout << "     " << std::setw(5) << tc << "   " << std::setw(20) << proj / 1e6
                      << "   " << std::setw(26) << ana / 1e6 << "   " << proj - ana << "\n";
         }
         std::cout << std::flush;
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }

   // ------------------------------------------------------------ outputs
   // station adapter: owned DOF positions -> TPV102StationWriter
   std::vector<Vector> fcoords(no);
   for (int j = 0; j < no; j++)
   {
      fcoords[j].SetSize(3);
      for (int c = 0; c < 3; c++) { fcoords[j](c) = op.FaultPositions()(3 * j + c); }
   }
   std::unique_ptr<TPV102StationWriter> stations;
   if (cfg.problem.tag == "tpv102")
   {
      stations = std::make_unique<TPV102StationWriter>();
      stations->Open(out_dir, prefix, DefaultStations(), fcoords, no, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::ofstream h(out_dir + "/" + prefix + "_station_conventions.txt");
         h << "# Nitsche dynamic-rupture station conventions (spec I3)\n"
           << "# frame: nhat = pin reference normal (" << pn[0] << ", " << pn[1] << ", " << pn[2]
           << "), up = (" << pu[0] << ", " << pu[1] << ", " << pu[2] << ")\n"
           << "#        shat = up x nhat / |up x nhat| = (" << ps[0] << ", " << ps[1] << ", " << ps[2]
           << "), dhat = shat x nhat = (" << pd[0] << ", " << pd[1] << ", " << pd[2] << ")\n"
           << "# pin: orientation propagated, pinned by the reference normal; plus side = side nhat points into\n"
           << "# columns (TPV102StationWriter layout): t slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta\n"
           << "#   slip1 = delta.dhat (SCEC v-slip)      slip2 = delta.shat (SCEC h-slip)\n"
           << "#   V1    = V.dhat     (SCEC v-slip-rate) V2    = V.shat     (SCEC h-slip-rate)\n"
           << "#   tau1  = -T.dhat    (SCEC v-shear-stress, Pa)  tau2 = -T.shat (SCEC h-shear-stress, Pa)\n"
           << "#   sigma_n = -T.nhat (compression positive, Pa); T = total traction on the face with normal nhat\n"
           << "#   delta = u(+) - u(-), V = d delta/dt (SCEC slip definition); log10_theta from psi\n";
      }
   }
   std::vector<DOFData> dofdata(no);
   auto fill_dofdata = [&]() {
      const Vector &d = op.Slip(), &V = op.SlipRate(), &psi = op.State();
      Vector Xs = op.LastShear(), sg = op.LastSigmaEff();
      const bool have_x = Xs.Size() == 2 * no;
      if (!have_x)
      {
         // before the first step: no bulk traction yet, the total traction
         // is the projected background (prestress + nucleation at t)
         Vector S;
         op.BackgroundTraction(op.GetTime(), S);
         Xs.SetSize(2 * no);
         sg.SetSize(no);
         for (int j = 0; j < no; j++)
         {
            const real_t *n = &op.FaultFrames()(9 * j), *dd = &op.FaultFrames()(9 * j + 3),
                          *ss = &op.FaultFrames()(9 * j + 6);
            ff::ToLocal(&S(3 * j), 0.0, n, dd, ss, Xs(2 * j), Xs(2 * j + 1), sg(j));
         }
      }
      for (int j = 0; j < no; j++)
      {
         const real_t *dd = &op.FaultFrames()(9 * j + 3), *ss = &op.FaultFrames()(9 * j + 6);
         DOFData &q = dofdata[j];
         q.slip1 = ff::Dot(&d(3 * j), dd);
         q.slip2 = ff::Dot(&d(3 * j), ss);
         q.V1 = ff::Dot(&V(3 * j), dd);
         q.V2 = ff::Dot(&V(3 * j), ss);
         const real_t g = have_x ? 0.5 * op.FaultSigmaF()(j) * dt : 0.0;
         // endpoint traction tangential components in the frame: X* + gamma V
         const real_t X1 = Xs(2 * j) + g * q.V1;
         const real_t X2 = Xs(2 * j + 1) + g * q.V2;
         q.tau1_corr = -X1;
         q.tau2_corr = -X2;
         q.sigma_n_corr = sg(j);
         q.psi = psi(j);
         q.Dc = op.FaultDc()(j);
         q.slip_rate = std::hypot(q.V1, q.V2);
      }
   };
   std::unique_ptr<std::ofstream> dump;
   Vector cf;
   if (dump_every > 0)
   {
      dump = std::make_unique<std::ofstream>(out_dir + "/" + prefix + "_fault_r"
                                             + std::to_string(g_rank) + ".csv");
      *dump << std::setprecision(12);
      // face centroid per owned DOF (keys DG copies uniquely)
      ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
      const int nd = E.GetNumFaultDOFs(), nbf = E.GetNbfPerFace();
      Vector cl(3 * nd);
      for (int k = 0; k < nd; k++)
      {
         long long key[3];
         real_t xv[3][3];
         op.Orientation().FaceVertexData(k / nbf, key, xv);
         for (int c = 0; c < 3; c++) { cl(3 * k + c) = (xv[0][c] + xv[1][c] + xv[2][c]) / 3.0; }
      }
      E.RestrictToOwnedFault(cl, cf, 3);
   }
   std::vector<real_t> trupt(no, -1.0);

   std::unique_ptr<NitscheEnergyMonitor> mon;

   const long nsteps = std::lround(std::ceil((tfinal - cfg.time.t_initial) / dt - 1e-9));
   if (g_rank == 0)
   {
      long long nloc = op.Height();
      std::cout << "  mesh " << mesh_path << ": " << ne_global << " tets; p = " << p
                << "; law " << law.name << "; fault faces " << op.Orientation().GetReport().faces
                << ", shared fault faces " << op.NumSharedFaultFacesGlobal() << "\n"
                << "  dt = " << dt << " s (CFL walk: cfl " << cfl << " x C(p) " << kCeig[p]
                << " -> " << dt_cfl << "; dt_max " << cfg.time.dt_max << "), " << nsteps
                << " steps to t = " << tfinal << " s; penalty_factor " << penalty_factor
                << "; energy balance " << (energy_on ? "on" : "off") << "\n";
      (void)nloc;
   }

   // ------------------------------------------------------------ loop
   const int n = op.Height();
   Vector U(n), Ud(n);
   U = 0.0;
   Ud = 0.0;
   real_t t = cfg.time.t_initial;
   op.SetTime(t);
   if (energy_on)
   {
      mon = std::make_unique<NitscheEnergyMonitor>(op);
      Vector A0(n);
      op.Mult(U, Ud, A0);
      mon->Start(mon->Evaluate(U, Ud, A0, t));
   }
   CentralDifferenceSolver ode;
   ode.Init(op);
   if (stations) { fill_dofdata(); stations->WriteStep(t, dofdata); }
   real_t max_ratio_run = 0.0, max_V_run = 0.0, min_sig_run = 1e300, max_sig_run = -1e300;
   for (long s = 1; s <= nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
      const NitscheStepRecord &rec = op.LastStep();
      max_ratio_run = std::max(max_ratio_run, rec.max_V_dt_over_Dc);
      max_V_run = std::max(max_V_run, rec.max_V);
      min_sig_run = std::min(min_sig_run, rec.min_sigma_eff);
      max_sig_run = std::max(max_sig_run, rec.max_sigma_eff);
      MFEM_VERIFY(std::isfinite(rec.max_V), "non-finite slip rate at step " << s);
      const Vector &V = op.SlipRate();
      for (int j = 0; j < no; j++)
      {
         if (trupt[j] < 0.0 && std::sqrt(V(3 * j) * V(3 * j) + V(3 * j + 1) * V(3 * j + 1)
                                         + V(3 * j + 2) * V(3 * j + 2)) > 1e-3) { trupt[j] = t; }
      }
      real_t rho = 0.0;
      NitscheEnergyRecord er;
      if (mon)
      {
         er = mon->Evaluate(U, Ud, op.LastAcceleration(), t);
         rho = mon->Step(er);
      }
      if (stations && s % station_every == 0) { fill_dofdata(); stations->WriteStep(t, dofdata); }
      if (dump && s % dump_every == 0)
      {
         const Vector &d = op.Slip(), &psi = op.State();
         const Vector &Xs = op.LastShear(), &sg = op.LastSigmaEff();
         const Vector &X = op.FaultPositions();
         for (int j = 0; j < no; j++)
         {
            *dump << s << "," << t << "," << X(3 * j) << "," << X(3 * j + 1) << "," << X(3 * j + 2)
                  << "," << cf(3 * j) << "," << cf(3 * j + 1) << "," << cf(3 * j + 2)
                  << "," << d(3 * j) << "," << d(3 * j + 1) << "," << d(3 * j + 2)
                  << "," << V(3 * j) << "," << V(3 * j + 1) << "," << V(3 * j + 2)
                  << "," << psi(j) << "," << Xs(2 * j) << "," << Xs(2 * j + 1) << "," << sg(j) << "\n";
         }
      }
      if (g_rank == 0 && (s % log_every == 0 || s == nsteps))
      {
         const real_t wall = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - wall0).count();
         std::cout << std::setprecision(6) << "  step " << s << " t = " << t << "  max|V| = "
                   << rec.max_V << "  max|V|dt/Dc = " << rec.max_V_dt_over_Dc
                   << "  sigma_n_eff in [" << rec.min_sigma_eff << ", " << rec.max_sigma_eff << "]";
         if (mon) { std::cout << "  E_h = " << er.E_h << "  sum rho = " << mon->SumRho(); }
         std::cout << "  wall " << wall << " s\n" << std::flush;
         (void)rho;
      }
   }
   if (stations) { stations->Close(); }

   // rupture times
   {
      std::ofstream rt(out_dir + "/" + prefix + "_rupture_time_r" + std::to_string(g_rank) + ".csv");
      rt << std::setprecision(12) << "# x y z t_rupture (first |V| > 1e-3 m/s; -1 = never)\n";
      const Vector &X = op.FaultPositions();
      for (int j = 0; j < no; j++)
      {
         rt << X(3 * j) << " " << X(3 * j + 1) << " " << X(3 * j + 2) << " " << trupt[j] << "\n";
      }
   }
   const real_t wall = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - wall0).count();
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "=== done: " << nsteps << " steps, t = " << t
                << " s; max|V| " << max_V_run << " m/s; max|V|dt/Dc " << max_ratio_run
                << "; sigma_n_eff in [" << min_sig_run << ", " << max_sig_run << "] Pa; wall "
                << wall << " s ===\n";
   }
   return 0;
}
