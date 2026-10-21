// Poroelastic plan rev 3, phase P0 (steps P0.1-P0.3): the option, the hooks
// in the elastic operator, and the static pore pressure.
//
//   U0-3  Hook H1: Force(x, v, d) followed by the mass inverse gives the
//         acceleration of Mult and AccelerationWithSlip bit for bit, p = 1-3,
//         with and without absorbing faces and the Kelvin-Voigt layer.
//   U0-4  ParsePoroConfig: a valid table reads back exactly; every required
//         key missing, misspelt or of the wrong type, unknown keys, the
//         source and starting-pressure tables, --poroelastic 1 without the
//         table, the option off with a table present (warning), and the bulk
//         density check.  The plan's [multiphysics] poroelastic = true next
//         to a [multiphysics.poroelastic] table is not valid TOML (checked).
//   U0-5  Static P_p: the TPV102 box with sigma_yy = -120 MPa and P_p = 10 MPa
//         against sigma_yy = -110 MPa and P_p = 0, 500 steps: slip, slip rate,
//         state and displacement agree to 1e-10 relative (class R); without
//         P_p the two differ (control).  P_p = 0 set through the pressure path
//         is bit-identical to no P_p.  The plan asked the bulk to be bitwise;
//         it cannot be, because sigma_n' is (T + S).n - P_p in one case and
//         (T + S').n in the other, which differ by round-off, and the slip
//         rate feeds the bulk through G d.
//
// The driver-level checks (U0-4 abort text, U0-6 refusals) and U0-1/U0-2 (the
// off path, every existing test) are Makefile recipes; see test-poro-p0.
//
// Usage: mpirun -np 1|2 ./seas_test_poro_p0 [--u03] [--u04] [--u05]
//        (no flag = all)

#include "mfem.hpp"
#include "../../../domain/nitsche_wave_operator.hpp"
#include "../../../solver/nitsche_cfl.hpp"
#include "../../../friction/rate_state_closure.hpp"
#include "../../../config/tpv102_params.hpp"
#include "../poro_config.hpp"
#include "poro_fixtures.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static const real_t kCuse[4] = {0.0, 0.057, 0.11, 0.14};   // 0.5 C_eig (P3)

static BoundaryConfig TPV102BC(bool absorbing)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   if (absorbing) { bc.absorbing_attrs = {5}; }
   else
   {
      // a time-dependent far field, so Force exercises the Dirichlet load
      bc.dirichlet_attrs = {5};
      bc.default_dirichlet_func = [](const Vector &x, real_t t, Vector &u) {
         u.SetSize(3);
         u(0) = 1e-6 * t * x(1) / kL;
         u(1) = 0.0;
         u(2) = 2e-7 * t;
      };
   }
   return bc;
}

static MaterialField TPV102Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu,
                                      TPV102Params::rho);
}

/// TPV102 fault setup with sigma_yy = -sigma_yy_pa and an optional P_p.  The
/// TPV102 nucleation (3 km) is not resolved on this 15 km box, so a broad
/// gentle one (12 MPa, e-fold radius 20 km, TPV102's G(t)) makes the fault
/// slip at up to ~1e-6 m/s in the test window, where the slip rate depends
/// on sigma_n' through the rate-and-state law, without dynamic runaway.
static NitscheFaultSetup TPV102Setup(const IFaultClosure &cl, real_t sigma_yy_pa,
                                     bool with_pp, real_t pp)
{
   NitscheFaultSetup fs;
   fs.closure = &cl;
   fs.prestress = [sigma_yy_pa](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = TPV102Params::tau_ini;
      s[4] = -sigma_yy_pa;
   };
   fs.nucleation_spatial = [](const real_t x[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      const real_t dx = x[0] - TPV102Params::hypo_along_strike;
      const real_t dz = -x[2] - TPV102Params::hypo_down_dip;
      const real_t r2 = (dx * dx + dz * dz) / (20.0e3 * 20.0e3);
      s[1] = s[3] = 12.0e6 * std::exp(-r2);
   };
   fs.nucleation_temporal = [](real_t t) { return NucleationTemporal(t); };
   fs.parameters = [](const real_t x[3], FaultPointState &st) {
      st.a = ComputeA(x[0], -x[2]);
      st.Dc = TPV102Params::Dc;
   };
   fs.sigma_min = 0.0;
   if (with_pp) { fs.pore_pressure = [pp](const real_t[3]) { return pp; }; }
   return fs;
}

/// d = 0, V = V_ini along strike, psi in equilibrium with the projected
/// prestress and sigma_n' (as the driver seeds it).
static void SeedState(NitscheWaveOperator &op, const IFaultClosure &cl)
{
   const int no = op.NumOwnedFaultDOFs();
   Vector d(3 * no), V(3 * no), psi(no), S;
   d = 0.0;
   op.BackgroundTraction(0.0, S);
   for (int j = 0; j < no; j++)
   {
      const real_t *f = &op.FaultFrames()(9 * j);
      real_t X[2], sig;
      fault_frame::ToLocal(&S(3 * j), 0.0, f, f + 3, f + 6, X[0], X[1], sig);
      if (op.HasFaultPressure()) { sig -= op.FaultPressure()(j); }
      FaultPointState st = op.FaultParams(j);
      cl.SeedEquilibrium(std::hypot(X[0], X[1]), TPV102Params::V_ini, sig, st);
      psi(j) = st.psi;
      for (int c = 0; c < 3; c++) { V(3 * j + c) = TPV102Params::V_ini * f[6 + c]; }
   }
   op.SetFaultState(d, V, psi);
}

// ------------------------------------------------------------------- U0-3
static void TestU03(ParMesh &pm)
{
   if (g_rank == 0) { std::cout << "--- U0-3: hook H1, Force then the mass inverse\n"; }
   DieterichRuinaFriction fr{{TPV102Params::V0, TPV102Params::f0, TPV102Params::b,
                              TPV102Params::Dc}};
   RateStateClosure aging(fr, RateStateClosure::StateLaw::Aging);
   for (int p = 1; p <= 3; p++)
   {
      for (int absorb = 0; absorb <= 1; absorb++)
      {
         for (int kv = 0; kv <= 1; kv++)
         {
            NitscheFaultSetup fs = TPV102Setup(aging, TPV102Params::sigma_n, false, 0.0);
            NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(absorb), {}, &fs);
            const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
            op.SetTimeStep(dt);
            if (kv) { op.SetKelvinVoigt(0.1, NitscheWaveOperator::KVScale::Node, -1.0); }
            SeedState(op, aging);
            op.SetTime(0.3);
            const int n = op.Height(), no = op.NumOwnedFaultDOFs();
            Vector x(n), v(n), d(3 * no);
            x.Randomize(11 + g_rank);
            v.Randomize(23 + g_rank);
            d.Randomize(37 + g_rank);
            x *= 1e-3;
            v *= 1e-2;
            Vector a_ref(n), a_new(n), F;
            // Mult: the stored slip
            op.Mult(x, v, a_ref);
            op.Force(x, &v, &op.Slip(), F);
            op.Mass().Apply(F, a_new);
            const long b1 = LSum(BitDiff(a_ref, a_new), pm.GetComm());
            // AccelerationWithSlip: a given slip, the damped mass with absorbing faces
            op.AccelerationWithSlip(x, v, d, a_ref);
            const Vector *vp = (op.HasDamping() || kv) ? &v : nullptr;
            op.Force(x, vp, &d, F);
            if (op.HasDamping()) { op.DampedMass()->Apply(F, a_new); }
            else { op.Mass().Apply(F, a_new); }
            const long b2 = LSum(BitDiff(a_ref, a_new), pm.GetComm());
            std::ostringstream tag;
            tag << "U0-3 p = " << p << (absorb ? ", absorbing" : ", no absorbing")
                << (kv ? ", Kelvin-Voigt" : "");
            TEST_ASSERT(b1 == 0, tag.str() << ": Mult == Force + Mass().Apply bit for bit");
            TEST_ASSERT(b2 == 0, tag.str() << ": AccelerationWithSlip == Force + "
                        << (op.HasDamping() ? "DampedMass" : "Mass") << "().Apply bit for bit");
         }
      }
   }
}

// ------------------------------------------------------------------- U0-4
static const char *kValid =
   "[multiphysics.poroelastic]\n"
   "enabled = true\n"
   "porosity = 0.4\n"
   "permeability_m2 = 6.0e-13\n"
   "viscosity_pa_s = 1.0e-3\n"
   "tortuosity = 2.0\n"
   "rho_fluid = 1040.0\n"
   "rho_solid = 2080.0\n"
   "bulk_fluid_pa = 2.5e9\n"
   "bulk_solid_pa = 2.0e10\n"
   "fault_hydraulics = \"permeable\"\n"
   "background_p_pa = 0\n"
   "fluid_penalty = 1.0\n"
   "drag_gamma = \"fitted\"\n"
   "velocity_coupling = \"momentum\"\n";

static const char *kInitTable =
   "[multiphysics.poroelastic.initial_pressure]\n"
   "origin_m = [0.0, -1000.0, -7500.0]\n"
   "magnitude_pa = 5.0e7\n"
   "width_per_m2 = 1.0e-7\n";

static const char *kKeys[] = {"enabled", "porosity", "permeability_m2", "viscosity_pa_s",
                              "tortuosity", "rho_fluid", "rho_solid", "bulk_fluid_pa",
                              "bulk_solid_pa", "fault_hydraulics", "background_p_pa",
                              "fluid_penalty", "drag_gamma", "velocity_coupling"};

static poro::PoroParseResult Parse(const std::string &deck, int cli = -1)
{
   std::istringstream is(deck);
   const toml::value root = toml::parse(is, "u04.toml");
   return poro::ParsePoroConfig(root, cli);
}

/// kValid with the line of `key` removed or replaced by `line`.
static std::string Edit(const std::string &key, const std::string &line)
{
   std::istringstream is(kValid);
   std::ostringstream os;
   for (std::string l; std::getline(is, l);)
   {
      if (l.rfind(key + " =", 0) == 0) { if (!line.empty()) { os << line << "\n"; } }
      else { os << l << "\n"; }
   }
   return os.str();
}

static bool Has(const std::string &s, const std::string &sub) { return s.find(sub) != std::string::npos; }

static void TestU04()
{
   if (g_rank == 0) { std::cout << "--- U0-4: the parser\n"; }
   // a valid table reads back exactly
   {
      const poro::PoroParseResult r = Parse(std::string(kValid) + kInitTable);
      const poro::PoroConfig &c = r.cfg;
      TEST_ASSERT(r.error.empty() && c.enabled, "U0-4 valid tpv5-poro table parses (" << r.error << ")");
      TEST_ASSERT(c.porosity == 0.4 && c.permeability_m2 == 6.0e-13 && c.viscosity_pa_s == 1.0e-3
                  && c.tortuosity == 2.0 && c.rho_fluid == 1040.0 && c.rho_solid == 2080.0
                  && c.bulk_fluid_pa == 2.5e9 && c.bulk_solid_pa == 2.0e10
                  && c.fault_hydraulics == "permeable" && c.background_p_pa == 0.0
                  && c.fluid_penalty == 1.0 && c.drag_gamma == "fitted"
                  && c.velocity_coupling == "momentum",
                  "U0-4 every value reads back exactly (integer 0 accepted as a number)");
      TEST_ASSERT(c.initial_pressure.on && c.initial_pressure.origin_m[1] == -1000.0
                  && c.initial_pressure.origin_m[2] == -7500.0
                  && c.initial_pressure.magnitude_pa == 5.0e7
                  && c.initial_pressure.width_per_m2 == 1.0e-7 && !c.source.on,
                  "U0-4 the starting-pressure table reads back; no source");
      TEST_ASSERT(std::abs(c.BulkDensity() - 1664.0) <= 1e-12 * 1664.0,
                  "U0-4 bulk density (1 - phi) rho_s + phi rho_f = 1664 kg/m^3");
      TEST_ASSERT(poro::CheckBulkDensity(c, 1664.0).empty(), "U0-4 rho 1664 accepted");
      const std::string e = poro::CheckBulkDensity(c, 2670.0);
      TEST_ASSERT(Has(e, "rho = 2670") && Has(e, "1664"), "U0-4 rho mismatch refused: " << e);
   }
   // every required key missing, then of the wrong type
   int miss_ok = 0, type_ok = 0;
   const int nk = sizeof(kKeys) / sizeof(kKeys[0]);
   for (int k = 0; k < nk; k++)
   {
      const std::string key = kKeys[k];
      const poro::PoroParseResult r = Parse(Edit(key, ""));
      if (Has(r.error, key + ": missing")) { miss_ok++; }
      else if (g_rank == 0) { std::cout << "    missing " << key << " gave: '" << r.error << "'\n"; }
      const bool is_string = key == "fault_hydraulics" || key == "drag_gamma"
                             || key == "velocity_coupling";
      const std::string bad = key + (key == "enabled" ? " = 1" : (is_string ? " = 1.0" : " = \"x\""));
      const poro::PoroParseResult r2 = Parse(Edit(key, bad));
      if (Has(r2.error, key + ": must be")) { type_ok++; }
      else if (g_rank == 0) { std::cout << "    wrong type " << key << " gave: '" << r2.error << "'\n"; }
   }
   TEST_ASSERT(miss_ok == nk, "U0-4 each of the " << nk << " keys missing: error names it ("
               << miss_ok << "/" << nk << ")");
   TEST_ASSERT(type_ok == nk, "U0-4 each of the " << nk << " keys of the wrong type: error names it ("
               << type_ok << "/" << nk << ")");
   // misspelt and unknown keys
   {
      const poro::PoroParseResult r = Parse(Edit("porosity", "porosty = 0.4"));
      TEST_ASSERT(Has(r.error, "porosty: unknown key"), "U0-4 misspelt key named: " << r.error);
      const poro::PoroParseResult r2 = Parse(std::string(kValid) + "[multiphysics]\nthermal = 1\n");
      TEST_ASSERT(Has(r2.error, "[multiphysics] thermal: unknown key"),
                  "U0-4 unknown key in [multiphysics] named: " << r2.error);
      const poro::PoroParseResult r3 = Parse(std::string(kValid) + kInitTable + "radius_m = 3.0\n");
      TEST_ASSERT(Has(r3.error, "initial_pressure] radius_m: unknown key"),
                  "U0-4 unknown key in the starting-pressure table named: " << r3.error);
   }
   // values out of range and the choices
   {
      const poro::PoroParseResult r = Parse(Edit("porosity", "porosity = 1.0"));
      TEST_ASSERT(Has(r.error, "porosity: must be in (0, 1)"), "U0-4 porosity 1 refused: " << r.error);
      const poro::PoroParseResult r2 = Parse(Edit("fault_hydraulics", "fault_hydraulics = \"sealed\""));
      TEST_ASSERT(Has(r2.error, "arrives in step P6"), "U0-4 sealed fault refused until P6: " << r2.error);
      const poro::PoroParseResult r3 = Parse(Edit("drag_gamma", "drag_gamma = \"exact\""));
      TEST_ASSERT(Has(r3.error, "drag_gamma: \"exact\" is not one of"), "U0-4 unknown choice named: " << r3.error);
      const poro::PoroParseResult r4 = Parse(Edit("tortuosity", "tortuosity = 0.5"));
      TEST_ASSERT(Has(r4.error, "tortuosity: must be >= 1"), "U0-4 tortuosity < 1 refused: " << r4.error);
   }
   // the source table
   {
      const std::string g = "[multiphysics.poroelastic.source]\nshape = \"gaussian\"\n"
                            "origin_m = [0.0, 0.0, -5000.0]\nrate = 1.0e-6\n";
      const poro::PoroParseResult r1 = Parse(std::string(kValid) + g + "width_per_m2 = 1.0e-7\n"
                                             "time_function = \"smooth_step\"\nrise_time_s = 0.1\n");
      TEST_ASSERT(r1.error.empty() && r1.cfg.source.on && r1.cfg.source.rise_time_s == 0.1,
                  "U0-4 gaussian smooth-step source parses (" << r1.error << ")");
      const poro::PoroParseResult r2 = Parse(std::string(kValid) + g + "time_function = \"step\"\n");
      TEST_ASSERT(Has(r2.error, "width_per_m2: missing"), "U0-4 gaussian without width refused: " << r2.error);
      const poro::PoroParseResult r3 = Parse(std::string(kValid) + g + "width_per_m2 = 1.0e-7\n"
                                             "time_function = \"step\"\nrise_time_s = 0.1\n");
      TEST_ASSERT(Has(r3.error, "rise_time_s: applies to"), "U0-4 rise time with a step refused: " << r3.error);
      const poro::PoroParseResult r4 = Parse(std::string(kValid) + g + "width_per_m2 = 1.0e-7\n"
                                             "time_function = \"ramp\"\n");
      TEST_ASSERT(Has(r4.error, "rise_time_s: missing"), "U0-4 ramp without rise time refused: " << r4.error);
      const std::string w = "[multiphysics.poroelastic.source]\nshape = \"well\"\n"
                            "origin_m = [0.0, 0.0, -5000.0]\nrate = 1.0e-3\ntime_function = \"step\"\n";
      const poro::PoroParseResult r5 = Parse(std::string(kValid) + w);
      TEST_ASSERT(r5.error.empty() && r5.cfg.source.shape == "well", "U0-4 well source parses ("
                  << r5.error << ")");
      const poro::PoroParseResult r6 = Parse(std::string(kValid) + w + "width_per_m2 = 1.0e-7\n");
      TEST_ASSERT(Has(r6.error, "width_per_m2: applies to shape = \"gaussian\" only"),
                  "U0-4 width with a well refused: " << r6.error);
   }
   // on and off
   {
      const poro::PoroParseResult r1 = Parse("[numerics]\ncfl = 0.5\n", 1);
      TEST_ASSERT(Has(r1.error, "--poroelastic 1 needs the [multiphysics.poroelastic] table"),
                  "U0-4 --poroelastic 1 without the table refused: " << r1.error);
      const poro::PoroParseResult r2 = Parse("[numerics]\ncfl = 0.5\n");
      TEST_ASSERT(r2.error.empty() && !r2.cfg.enabled && r2.cfg.warning.empty(),
                  "U0-4 no table: off, no warning");
      const poro::PoroParseResult r3 = Parse(Edit("enabled", "enabled = false"));
      TEST_ASSERT(r3.error.empty() && !r3.cfg.enabled && Has(r3.cfg.warning, "is ignored"),
                  "U0-4 enabled = false with a table: off, warning '" << r3.cfg.warning << "'");
      const poro::PoroParseResult r4 = Parse(kValid, 0);
      TEST_ASSERT(r4.error.empty() && !r4.cfg.enabled && Has(r4.cfg.warning, "--poroelastic 0"),
                  "U0-4 --poroelastic 0 overrides enabled = true, with a warning");
      const poro::PoroParseResult r5 = Parse(Edit("enabled", "enabled = false"), 1);
      TEST_ASSERT(r5.error.empty() && r5.cfg.enabled, "U0-4 --poroelastic 1 overrides enabled = false");
      const poro::PoroParseResult r6 = Parse(Edit("enabled", "enabled = false")
                                             + "porosty = 3\n");
      TEST_ASSERT(r6.error.empty() && !r6.cfg.enabled,
                  "U0-4 with the option off the table is not checked (ignored, as warned)");
      bool threw = false;
      try
      {
         std::istringstream is("[multiphysics]\nporoelastic = true\n[multiphysics.poroelastic]\n"
                               "porosity = 0.4\n");
         (void)toml::parse(is, "plan_layout.toml");
      }
      catch (const std::exception &) { threw = true; }
      TEST_ASSERT(threw, "U0-4 the plan's layout ([multiphysics] poroelastic = true and a "
                  "[multiphysics.poroelastic] table) is invalid TOML, hence the key 'enabled'");
   }
}

// ------------------------------------------------------------------- U0-5
struct RunState { Vector U, d, V, psi, sig; };

static RunState RunTPV102(ParMesh &pm, int p, real_t sigma_yy, bool with_pp, real_t pp,
                          int nsteps)
{
   DieterichRuinaFriction fr{{TPV102Params::V0, TPV102Params::f0, TPV102Params::b,
                              TPV102Params::Dc}};
   RateStateClosure aging(fr, RateStateClosure::StateLaw::Aging);
   NitscheFaultSetup fs = TPV102Setup(aging, sigma_yy, with_pp, pp);
   NitscheWaveOperator op(pm, p, TPV102Material(), TPV102BC(true), {}, &fs);
   const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   op.SetTimeStep(dt);
   SeedState(op, aging);
   const int n = op.Height();
   Vector U(n), Ud(n);
   U = 0.0;
   Ud = 0.0;
   real_t t = 0.0;
   op.SetTime(t);
   CentralDifferenceSolver ode;
   ode.Init(op);
   for (int s = 0; s < nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
   }
   RunState r;
   r.U = U;
   r.d = op.Slip();
   r.V = op.SlipRate();
   r.psi = op.State();
   r.sig = op.LastSigmaEff();
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "    sigma_yy = -" << sigma_yy / 1e6 << " MPa, P_p "
                << (with_pp ? std::to_string(pp / 1e6) + " MPa" : std::string("unset"))
                << ": " << nsteps << " steps of " << dt << " s, max V " << op.LastStep().max_V
                << " m/s, sigma_n' in [" << op.LastStep().min_sigma_eff / 1e6 << ", "
                << op.LastStep().max_sigma_eff / 1e6 << "] MPa\n";
   }
   return r;
}

static void TestU05(ParMesh &pm, int nsteps)
{
   if (g_rank == 0) { std::cout << "--- U0-5: static pore pressure\n"; }
   const MPI_Comm c = pm.GetComm();
   for (int p = 1; p <= 2; p++)
   {
      const RunState A = RunTPV102(pm, p, 120e6, true, 10e6, nsteps);
      const RunState B = RunTPV102(pm, p, 110e6, false, 0.0, nsteps);
      const RunState Z = RunTPV102(pm, p, 110e6, true, 0.0, nsteps);
      const RunState N = RunTPV102(pm, p, 120e6, false, 0.0, nsteps);   // control
      const real_t ed = RelDiff(A.d, B.d, c), eV = RelDiff(A.V, B.V, c);
      const real_t ep = RelDiff(A.psi, B.psi, c), es = RelDiff(A.sig, B.sig, c);
      const real_t eU = RelDiff(A.U, B.U, c);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(3) << "    p = " << p << ": relative differences slip "
                   << ed << ", slip rate " << eV << ", state " << ep << ", sigma_n' " << es
                   << ", displacement " << eU << "\n";
      }
      const real_t vmax = RMax(B.V.Normlinf(), c);
      TEST_ASSERT(vmax > 1e-8, "U0-5 p = " << p << ": the fault slips in the test window, so "
                  "the law sees sigma_n' (max |V| = " << vmax << " m/s)");
      TEST_ASSERT(ed <= 1e-10 && eV <= 1e-10 && ep <= 1e-10 && es <= 1e-10,
                  "U0-5 p = " << p << ": P_p = 10 MPa with sigma_yy = -120 MPa matches "
                  "sigma_yy = -110 MPa on the fault to 1e-10");
      TEST_ASSERT(eU <= 1e-10, "U0-5 p = " << p << ": the displacement matches to 1e-10 "
                  "(it sees the stress only through the slip)");
      const real_t eN = RelDiff(N.V, B.V, c);
      TEST_ASSERT(eN > 1e-3, "U0-5 p = " << p << ": control: without P_p, sigma_yy = -120 MPa "
                  "differs from -110 MPa (slip rate " << eN << " relative), so the match above "
                  "is the subtraction's doing");
      const long bz = LSum(BitDiff(Z.U, B.U) + BitDiff(Z.d, B.d) + BitDiff(Z.V, B.V)
                           + BitDiff(Z.psi, B.psi), c);
      TEST_ASSERT(bz == 0, "U0-5 p = " << p << ": P_p = 0 through the pressure path is "
                  "bit-identical to no P_p");
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   bool u03 = false, u04 = false, u05 = false;
   int nsteps = 500;
   for (int i = 1; i < argc; i++)
   {
      const std::string a = argv[i];
      if (a == "--u03") { u03 = true; }
      else if (a == "--u04") { u04 = true; }
      else if (a == "--u05") { u05 = true; }
      else if (a == "--steps" && i + 1 < argc) { nsteps = std::atoi(argv[++i]); }
   }
   if (!u03 && !u04 && !u05) { u03 = u04 = u05 = true; }
   if (g_rank == 0)
   {
      std::cout << "=== Poroelastic P0 (plan rev 3) on " << Mpi::WorldSize() << " rank(s) ===\n";
   }
   if (u04 && g_rank == 0) { TestU04(); }   // pure parsing, rank 0
   std::unique_ptr<Mesh> m(TPV102Box(4, 4, 2));
   ParMesh pm(MPI_COMM_WORLD, *m);
   if (u03) { TestU03(pm); }
   if (u05) { TestU05(pm, nsteps); }
   int nf = num_failed;
   MPI_Allreduce(MPI_IN_PLACE, &nf, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "=== Summary: " << num_passed << "/" << num_tests << " passed, "
                << num_failed << " failed ===\n";
   }
   return nf == 0 ? 0 : 1;
}
