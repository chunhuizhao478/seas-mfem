// Poroelastic plan rev 3, phase P0 (steps P0.1-P0.2): the option and the
// hooks in the elastic operator.
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
//
// The driver-level checks (U0-4 abort text, U0-6 refusals) and U0-1/U0-2 (the
// off path, every existing test) are Makefile recipes; see test-poro-p0.
//
// Usage: mpirun -np 1|2 ./seas_test_poro_p0 [--u03] [--u04]
//        (no flag = all)

#include "mfem.hpp"
#include "../../../domain/nitsche_wave_operator.hpp"
#include "../../../solver/nitsche_cfl.hpp"
#include "../../../friction/rate_state_closure.hpp"
#include "../../../config/tpv102_params.hpp"
#include "../poro_config.hpp"

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

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kL = 30.0e3;
static const real_t kCuse[4] = {0.0, 0.057, 0.11, 0.14};   // 0.5 C_eig (P3)

static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static long LSum(long v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_LONG, MPI_SUM, c); return v; }

static long BitDiff(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return -1; }
   long n = 0;
   for (int i = 0; i < a.Size(); i++)
   {
      if (std::memcmp(&a(i), &b(i), sizeof(real_t)) != 0) { n++; }
   }
   return n;
}

/// max_i |a_i - b_i| / max_i |b_i| over all ranks.
static real_t RelDiff(const Vector &a, const Vector &b, MPI_Comm c)
{
   real_t num = 0.0, den = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      num = std::max(num, std::abs(a(i) - b(i)));
      den = std::max(den, std::abs(b(i)));
   }
   num = RMax(num, c);
   den = RMax(den, c);
   return den > 0.0 ? num / den : num;
}

// ---------------------------------------------------------------- fixture
// F-B of plan rev 3, section 6.1: the TPV102 inline box of
// test_nitsche_hook_spy.cpp (24 tets per hex, fault on Y = 0).

static Mesh *TPV102Box(int nx, int ny, int nz)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2 * kL, 2 * kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
   }
   const real_t tol = kL / nz * 0.1;
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      real_t zc = 0.0;
      for (int i = 0; i < verts.Size(); i++) { zc += mesh->GetVertex(verts[i])[2]; }
      zc /= verts.Size();
      mesh->SetBdrAttribute(b, (std::abs(zc) < tol) ? 1 : 5);
   }
   const real_t ytol = 1e-6 * kL;
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool on_plane = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         if (std::abs(mesh->GetVertex(fv[i])[1]) > ytol) { on_plane = false; }
      }
      if (on_plane) { tris.push_back({fv[0], fv[1], fv[2]}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   return mesh;
}

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

/// TPV102 fault setup with sigma_yy = -sigma_yy_pa.  The TPV102 nucleation
/// (3 km) is not resolved on this 15 km box; a broad gentle one (12 MPa,
/// e-fold radius 20 km, TPV102's G(t)) stands in for it.
static NitscheFaultSetup TPV102Setup(const IFaultClosure &cl, real_t sigma_yy_pa)
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
   return fs;
}

/// d = 0, V = V_ini along strike, psi in equilibrium with the projected
/// prestress (as the driver seeds it).
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
            NitscheFaultSetup fs = TPV102Setup(aging, TPV102Params::sigma_n);
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

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   bool u03 = false, u04 = false;
   for (int i = 1; i < argc; i++)
   {
      const std::string a = argv[i];
      if (a == "--u03") { u03 = true; }
      else if (a == "--u04") { u04 = true; }
   }
   if (!u03 && !u04) { u03 = u04 = true; }
   if (g_rank == 0)
   {
      std::cout << "=== Poroelastic P0 (plan rev 3) on " << Mpi::WorldSize() << " rank(s) ===\n";
   }
   if (u04 && g_rank == 0) { TestU04(); }   // pure parsing, rank 0
   std::unique_ptr<Mesh> m(TPV102Box(4, 4, 2));
   ParMesh pm(MPI_COMM_WORLD, *m);
   if (u03) { TestU03(pm); }
   int nf = num_failed;
   MPI_Allreduce(MPI_IN_PLACE, &nf, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "=== Summary: " << num_passed << "/" << num_tests << " passed, "
                << num_failed << " failed ===\n";
   }
   return nf == 0 ? 0 : 1;
}
