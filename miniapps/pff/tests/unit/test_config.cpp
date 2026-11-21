// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests of the case configuration: every case file in <model>/configs/ parses, the
// derived values match the MOOSE inputs, unknown keys are refused, and Gmsh
// physical names resolve to attributes.
//
// Run from miniapps/pff (the Makefile's `make test` does).

#include "mfem.hpp"
#include <map>
#include "../../config/pff_config_parser.hpp"
#include "../../config/qs_config_parser.hpp"
#include "../../domain/qs_mesh.hpp"
#include "../../io/gmsh_physical_names.hpp"
#include "../test_framework.hpp"

#include <sstream>

using namespace mfem;
using namespace mfem::pff;

static void TestSpe()
{
   std::cout << "\n[spe_hftc cube case]\n";
   const CaseConfig c = CaseConfigParser::ParseFile("spe_hftc/configs/cube_homog_eta0p45_absorb.toml");
   TEST_ASSERT(c.model == CaseModel::Solid, "model = solid");
   TEST_NEAR(c.pf.Gc, 2720.5714285714284, 1e-9, "Gc derived from ft = 69 MPa at l = 0.6 mm");
   TEST_NEAR(c.Cp(), 1943.0, 1.0, "c_p = 1943 m/s (formulation write-up)");
   TEST_NEAR(c.Cs(), 934.0, 1.0, "c_s = 934 m/s");
   TEST_NEAR(c.end_time, 20 * 6e-5, 1e-15, "end_time defaults to Np Tw");
   TEST_ASSERT(c.absorbing_bdr.size() == 4 && c.fixed_z_bdr.size() == 2, "boundary names read");
   TEST_ASSERT(c.pf.driving == DrivingEnergy::Active, "driving = active");
   const CaseConfig s = CaseConfigParser::ParseFile("spe_hftc/configs/section_homog_eta0p45_absorb.toml");
   TEST_ASSERT(s.dim == 2 && s.fixed_z_bdr.empty(), "2D section: dim 2, no end faces");
}

static void TestCmame()
{
   std::cout << "\n[cmame_revision cases]\n";
   const CaseConfig s2 = CaseConfigParser::ParseFile("cmame_revision/configs/2d_solid_benchmark.toml");
   TEST_ASSERT(s2.initial_damage.size() == 2, "two seeded slots");
   TEST_NEAR(s2.initial_damage[0].value, 0.9, 0.0, "seed value 0.9");
   // ft = sqrt(3/8 E Gc / l) = 96.8 MPa for E = 50 GPa, Gc = 100, l = 0.2 mm.
   // (The "#137 MPa" comment next to ft in the MOOSE inputs is stale: MOOSE
   // evaluates the fparse expression, which gives 96.8 MPa.)
   TEST_NEAR(PhaseFieldModel::AT1Strength(s2.E, s2.pf.Gc, s2.pf.l), 96.82e6, 0.01e6,
             "2D solid: AT1 strength 96.8 MPa");
   TEST_NEAR(s2.confinement_pressure, 1e6, 0.0, "1 MPa confinement");

   const CaseConfig h2 = CaseConfigParser::ParseFile("cmame_revision/configs/2d_hydromech_benchmark.toml");
   TEST_ASSERT(h2.model == CaseModel::Hydromech, "model = hydromech");
   TEST_ASSERT(h2.poro.perm_model == PermeabilityModel::DarcyPoiseuille, "darcy_poiseuille");
   TEST_REL(h2.poro.wc, 100.0 / std::sqrt(3.0 / 8.0 * 50e9 * 100.0 / 2e-4), 1e-12,
            "w_c defaults to Gc / ft (elasticity_E1d25.i)");
   TEST_REL(h2.poro.K, 50e9 / 3.0 / (1.0 - 0.6), 1e-14, "poro K = intact drained bulk modulus");
   TEST_ASSERT(h2.pf.driving == DrivingEnergy::HalfHistoryEnhanced, "driving = half_history_enhanced");

   const CaseConfig h3 = CaseConfigParser::ParseFile("cmame_revision/configs/3d_hydromech_em_5J.toml");
   TEST_ASSERT(h3.poro.perm_model == PermeabilityModel::NormalStrain, "normal_strain permeability");
   TEST_NEAR(h3.initial_pore_pressure, 0.0965e6, 0.0, "p0 = 0.0965 MPa");
   TEST_NEAR(h3.end_time, 50e-5, 1e-18, "end_time 50e-5");

   const CaseConfig s3 = CaseConfigParser::ParseFile("cmame_revision/configs/3d_solid_em_5J.toml");
   TEST_ASSERT(s3.absorbing_bdr.size() == 3 && s3.absorbing_bdr[0] == "4", "3D solid dashpots on 4 5 2");
}

static void TestStrictKeys()
{
   std::cout << "\n[unknown keys are refused]\n";
   const std::string base =
      "[case]\nname='x'\n[mesh]\nfile='m.msh'\n[material]\nE=1e9\nnu=0.25\ndensity=1000\n"
      "[phase_field]\nGc=1\nl=1e-3\n[time]\nend_time=1\n";
   bool ok = true;
   try { CaseConfigParser::ParseString(base); } catch (...) { ok = false; }
   TEST_ASSERT(ok, "a minimal case parses");
   bool threw = false;
   try { CaseConfigParser::ParseString(base + "[damping]\nzeta_pulsee=1e-8\n"); }
   catch (const std::invalid_argument &) { threw = true; }
   TEST_ASSERT(threw, "misspelt key throws");
   threw = false;
   try { CaseConfigParser::ParseString(base + "[pulses]\nEM=1\n"); }
   catch (const std::invalid_argument &) { threw = true; }
   TEST_ASSERT(threw, "unknown table throws");
   threw = false;
   try
   {
      CaseConfigParser::ParseString(
         "[case]\nname='x'\nmodel='hydromech'\n[mesh]\nfile='m.msh'\n[material]\nE=1e9\nnu=0.25\n"
         "density=1000\n[phase_field]\nGc=1\nl=1e-3\n[time]\nend_time=1\n");
   }
   catch (const std::invalid_argument &) { threw = true; }
   TEST_ASSERT(threw, "hydromech without [poro] throws");
   const CaseConfig ci = CaseConfigParser::ParseString(
                            "[case]\nname='x'\n[mesh]\nfile='m.msh'\n[material]\nE=1000000000\n"
                            "nu=0.25\ndensity=1000\n[phase_field]\nGc=1\nl=1e-3\n[time]\nend_time=1\n");
   TEST_NEAR(ci.E, 1e9, 0.0, "an integer TOML value is accepted for a real");
}

static void TestPhysicalNames()
{
   std::cout << "\n[Gmsh $PhysicalNames]\n";
   std::istringstream msh("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n$PhysicalNames\n4\n"
                          "2 3 \"borehole\"\n2 5 \"xmin\"\n3 9 \"acrylic\"\n1 1 \"edge\"\n"
                          "$EndPhysicalNames\n$Entities\n");
   const PhysicalNameMap m = ReadPhysicalNames(msh);
   TEST_ASSERT(m.size() == 4, "four names read");
   const std::vector<int> a = ResolveAttributes({"borehole", "xmin", "7"}, m, 2);
   TEST_ASSERT(a.size() == 3 && a[0] == 3 && a[1] == 5 && a[2] == 7, "names and numbers resolve");
   bool threw = false;
   try { ResolveAttributes({"acrylic"}, m, 2); } catch (const std::invalid_argument &) { threw = true; }
   TEST_ASSERT(threw, "a volume name is not a boundary");
}

static void TestRaccoonMode1()
{
   std::cout << "\n[raccoon mode1 case]\n";
   const QsCaseConfig c = QsCaseConfigParser::ParseFile(
                             "raccoon/mode1_brittle_fracture/configs/mode1.toml");
   TEST_NEAR(c.K(), 2.1e5 / 3.0 / (1.0 - 0.6), 1e-9, "K = E/3/(1-2 nu) (elasticity.i)");
   TEST_NEAR(c.G(), 2.1e5 / 2.0 / 1.3, 1e-9, "G = E/2/(1+nu)");
   TEST_ASSERT(c.pf.crack.model == CrackModel::AT2 && c.pf.crack.c0() == 2.0, "AT2, c0 = 2");
   TEST_ASSERT(c.decomposition == Decomposition::None && c.driving_element_average,
               "no split, element-average driving energy");
   TEST_ASSERT(c.dirichlet.size() == 3 && c.dirichlet[0].rate == 1.0, "three Dirichlet BCs, u_y = t");
   TEST_ASSERT(c.automatic_scaling, "automatic_scaling = true (elasticity.i Executioner)");
   BoundaryNames names;
   Mesh m = BuildQsMesh(c, names);
   // tutorials/mode1_brittle_fracture/gold/elasticity_out.e: 888 elements, 994 nodes
   TEST_ASSERT(m.GetNE() == 888 && m.GetNV() == 994, "refined mesh matches the gold (888, 994)");
   TEST_ASSERT(names.at("noncrack") == 5, "noncrack node set is attribute 5");
   bool threw = false;
   try { QsCaseConfigParser::ParseString("[case]\nname='x'\n[mesh]\nnx=1\nny=1\n[material]\nE=1\nnu=0\n"
                                          "[phase_field]\nGc=1\nl=1\n[time]\ndt=1\n[solver]\nfixed_point_max_it=3\n"); }
   catch (const std::invalid_argument &) { threw = true; }
   TEST_ASSERT(threw, "misspelt [solver] key throws");
}

static void TestRaccoonMode2()
{
   std::cout << "\n[raccoon mode2 case]\n";
   const QsCaseConfig c = QsCaseConfigParser::ParseFile(
                             "raccoon/mode2_brittle_fracture/configs/mode2.toml");
   TEST_ASSERT(c.mesh_type == "slit" && c.slit_x == 0.5, "slit mesh stitched for x >= 0.5");
   TEST_ASSERT(c.decomposition == Decomposition::Spectral, "spectral split (elasticity.i)");
   TEST_ASSERT(!c.dirichlet[0].preset && c.dirichlet[1].preset,
               "xdisp preset = false, the DirichletBCs keep MOOSE's default preset = true");
   BoundaryNames names;
   Mesh m = BuildQsMesh(c, names);
   // M-R2 and tutorials/mode2_brittle_fracture/gold/elasticity_out.e: 2826
   // elements, 3026 nodes, 28 coincident node pairs on y = 0, x < 0.5.
   int pairs = 0;
   std::map<std::pair<long long, long long>, int> seen;
   for (int v = 0; v < m.GetNV(); v++)
   {
      const real_t *x = m.GetVertex(v);
      const auto key = std::make_pair(std::llround(x[0] * 1e9), std::llround(x[1] * 1e9));
      if (seen[key]++ == 1) { pairs++; }
   }
   TEST_ASSERT(m.GetNE() == 2826 && m.GetNV() == 3026, "refined slit mesh matches MOOSE (2826, 3026)");
   TEST_ASSERT(pairs == 28, "28 duplicated node pairs along the slit");
   TEST_ASSERT(names.at("top_half_top") == 3 && names.at("bottom_half_bottom") == 5,
               "boundary_name_prefix names");
}

int main()
{
   TestSpe();
   TestCmame();
   TestStrictKeys();
   TestPhysicalNames();
   TestRaccoonMode1();
   TestRaccoonMode2();
   return TestSummary("test_config");
}
