// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_2d_deck.cpp — 2D in-plane design rev 3, phase B
// (part 2): the deck's vectors in 2D (section 5.5).  A 2D deck is a 3D deck
// whose vectors have two entries; the readers pad z = 0 and record the
// entries, and the check against the mesh dimension
// (CycleConfig::CheckDimension, spatial::CheckFaultGeometryDimension)
// refuses a vector of the other dimension.
//
//   K1  the first-event deck (3-D) loads as before: every vector key read
//       with 3 entries, [fault_geometry] ref_normal and up 3 entries, and
//       CheckDimension(3) passes
//   K2  its 2-D form (every vector [x, y]; ref_normal [0, -1]; up [0, 0, 1],
//       out of the plane): the same keys with 2 entries, z = 0 in every
//       padded vector, the in-plane entries as written, up_entries 3, and
//       CheckDimension(2) passes
//   K3  up in the plane ([0, 1] with ref_normal [1, 0], a vertical section)
//       and up = [0, 0, -1] also pass CheckDimension(2)
//   K4  the spatial loader alone (ParseSpatialFrictionConfigString): a
//       2-entry [fault_geometry] vector is read with z = 0 and normalized,
//       a 3-entry one as before, an absent one keeps the default with
//       entries 0
//
// Modes: (default) K1-K4;  --refuse N (N = 1..10) loads a deck and checks
// it against a dimension, which must abort (the makefile checks the exit
// status and the message):
//    1 3-D, ref_normal [nx, ny]                   "ref_normal has 2 entries, but the mesh is 3-D"
//    2 3-D, up [ux, uy]                           "up has 2 entries, but the mesh is 3-D"
//    3 3-D, a station x [x, y]                    "has 2 entries, but the mesh is 3-D"
//    4 2-D, plane_origin [x, y, z]                "has 3 entries, but the mesh is 2-D"
//    5 2-D, no up                                 "up is required on a 2-D mesh"
//    6 2-D, up [0, 0.6, 0.8]                      "must be [0, 0, +-1]"
//    7 2-D, ref_normal [nx, ny, nz]               "ref_normal has 3 entries, but the mesh is 2-D"
//    8 2-D, no ref_normal                         "ref_normal is required on a 2-D mesh"
//    9 a 4-entry [fault_geometry] vector          "must be a 2- or 3-element array"
//   10 a 4-entry cycle vector                     "must be an array of 3 numbers (2 on a 2-D mesh)"
// The make target also runs the two drivers on a 3-D mesh with a 2-entry
// ref_normal, which they must refuse right after reading the mesh.
//
//   make -f tests/2d_inplane/2d_inplane.mk test-2d-deck
//
// Scope (Chunhui Zhao, 2026-10-25).

#include "mfem.hpp"
#include "../../cycle/cycle_config.hpp"

#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const char *kDeck = "config/fd_cycle/bp5_fd_first_event.toml";
static const int kDeckVectors = 22;   // 3 fields x 3, 3 output plane, 10 stations

static std::string ReadFile(const std::string &p)
{
   std::ifstream f(p);
   MFEM_VERIFY(f, "cannot read " << p);
   std::stringstream ss;
   ss << f.rdbuf();
   return ss.str();
}

static std::string WriteDeck(const std::string &tag, const std::string &text)
{
   const char *tmp = std::getenv("TMPDIR");
   const std::string path = std::string(tmp ? tmp : "/tmp") + "/seas_test_2d_deck_" + tag
                            + "_" + std::to_string(getpid()) + "_" + std::to_string(g_rank) + ".toml";
   std::ofstream f(path);
   f << text;
   MFEM_VERIFY(f, "cannot write " << path);
   return path;
}

static cycle::CycleConfig LoadText(const std::string &tag, const std::string &text)
{
   const std::string p = WriteDeck(tag, text);
   cycle::CycleConfig c = cycle::CycleConfig::Load(p);
   std::remove(p.c_str());
   return c;
}

// Replaces the one line `old` of the text by `rep` (verified to occur once).
static std::string ReplaceLine(std::string text, const std::string &old, const std::string &rep)
{
   const size_t at = text.find(old);
   MFEM_VERIFY(at != std::string::npos && text.find(old, at + 1) == std::string::npos,
               "test_2d_deck: the deck line '" << old << "' moved");
   return text.replace(at, old.size(), rep);
}

static const std::string kRef3 = "ref_normal      = [0.0, -1.0, 0.0]";
static const std::string kUp3 = "up              = [0.0,  0.0, 1.0]";

// The deck's 2-D form: every 3-number array [a, b, c] -> [a, b]; then the
// axes that drop to zero get in-plane values (axis2 = [0, 1] across the
// fault line, plane_axis2 = [0, -1]) and up goes back to [0, 0, 1].
static std::string To2D(const std::string &deck)
{
   const std::regex v3(R"(\[\s*([^,\]\s]+)\s*,\s*([^,\]\s]+)\s*,\s*([^,\]\s]+)\s*\])");
   std::string t = std::regex_replace(deck, v3, "[$1, $2]");
   t = std::regex_replace(t, std::regex(R"(axis2   = \[0\.0, 0\.0\])"), "axis2   = [0.0, 1.0]");
   t = ReplaceLine(t, "plane_axis2  = [0.0, 0.0]", "plane_axis2  = [0.0, -1.0]");
   t = ReplaceLine(t, "up              = [0.0, 0.0]", kUp3);
   return t;
}

static bool AllEntries(const cycle::CycleConfig &c, int n)
{
   for (const cycle::DeckVector &v : c.deck_vectors) { if (v.entries != n) { return false; } }
   return true;
}

// ------------------------------------------------------------------ K1, K2
static void TestDeck3D2D()
{
   if (g_rank == 0) { std::cout << "\n=== K1 the 3-D deck, K2 its 2-D form ===\n"; }
   const std::string deck = ReadFile(kDeck);
   const cycle::CycleConfig c3 = cycle::CycleConfig::Load(kDeck);
   const spatial::FaultGeometrySpec &g3 = c3.shared.fault_geometry;
   TEST_ASSERT(static_cast<int>(c3.deck_vectors.size()) == kDeckVectors && AllEntries(c3, 3),
               "K1 " << c3.deck_vectors.size() << " vector keys read, all with 3 entries");
   TEST_ASSERT(g3.ref_normal_entries == 3 && g3.up_entries == 3,
               "K1 [fault_geometry] ref_normal and up: 3 entries");
   c3.CheckDimension(3);
   TEST_ASSERT(true, "K1 CheckDimension(3) passes");

   const std::string d2 = To2D(deck);
   const cycle::CycleConfig c2 = LoadText("k2", d2);
   const spatial::FaultGeometrySpec &g2 = c2.shared.fault_geometry;
   bool keys = c2.deck_vectors.size() == c3.deck_vectors.size();
   for (size_t i = 0; keys && i < c2.deck_vectors.size(); i++)
   {
      keys = (c2.deck_vectors[i].key == c3.deck_vectors[i].key);
   }
   TEST_ASSERT(keys && AllEntries(c2, 2), "K2 the same " << c2.deck_vectors.size()
               << " keys, all with 2 entries");
   // padded: z = 0, x and y as the 3-D deck wrote them (axis2 and
   // plane_axis2 set by To2D)
   bool pad = true;
   for (size_t i = 0; i < c3.output.stations.size(); i++)
   {
      const auto &a = c3.output.stations[i].x, &b = c2.output.stations[i].x;
      pad = pad && b[0] == a[0] && b[1] == a[1] && b[2] == 0.0;
   }
   for (const cycle::RectTaperSpec *f : {&c2.fields.a, &c2.fields.Dc, &c2.fields.V_init})
   {
      pad = pad && f->center[2] == 0.0 && f->axis1[2] == 0.0 && f->axis2[0] == 0.0
            && f->axis2[1] == 1.0 && f->axis2[2] == 0.0;
   }
   pad = pad && c2.output.plane_axis2[0] == 0.0 && c2.output.plane_axis2[1] == -1.0
         && c2.output.plane_axis2[2] == 0.0;
   TEST_ASSERT(pad, "K2 every padded vector has z = 0 and the in-plane entries as written");
   TEST_ASSERT(g2.ref_normal_entries == 2 && g2.up_entries == 3 && g2.ref_normal[0] == 0.0
               && g2.ref_normal[1] == -1.0 && g2.ref_normal[2] == 0.0 && g2.up[2] == 1.0,
               "K2 [fault_geometry] ref_normal [0, -1] (2 entries, padded), up [0, 0, 1] (3)");
   TEST_ASSERT(c2.dyn.ref_normal == g2.ref_normal && c2.dyn.up == g2.up,
               "K2 the dynamic regime's ref_normal and up are the padded vectors");
   c2.CheckDimension(2);
   TEST_ASSERT(true, "K2 CheckDimension(2) passes");
}

// ---------------------------------------------------------------------- K3
static void TestUpChoices()
{
   if (g_rank == 0) { std::cout << "\n=== K3 the two up choices in 2-D ===\n"; }
   const std::string d2 = To2D(ReadFile(kDeck));
   // a vertical section: ref_normal [1, 0], up [0, 1] (BP3).  The default
   // [cycle.loading] kind "bp5" needs ref_normal (0, -1, 0), so the deck
   // takes a plate loading across the line x = 0.
   std::string t = ReplaceLine(d2, "ref_normal      = [0.0, -1.0]", "ref_normal      = [1.0, 0.0]");
   t = ReplaceLine(t, kUp3, "up              = [0.0, 1.0]");
   t += "\n[cycle.loading]\nkind = \"plate\"\nplane_point = [0.0, 0.0]\n"
        "plane_normal = [1.0, 0.0]\nslip_dir = [0.0, -1.0]\nband = 1000.0\n";
   const cycle::CycleConfig cb = LoadText("k3b", t);
   const spatial::FaultGeometrySpec &gb = cb.shared.fault_geometry;
   TEST_ASSERT(gb.up_entries == 2 && gb.up[0] == 0.0 && gb.up[1] == 1.0 && gb.up[2] == 0.0
               && AllEntries(cb, 2) && cb.deck_vectors.size() == kDeckVectors + 3,
               "K3 up [0, 1]: 2 entries, padded (0, 1, 0); the 3 [cycle.loading] vectors with 2");
   cb.CheckDimension(2);
   TEST_ASSERT(true, "K3 up in the plane: CheckDimension(2) passes");
   const cycle::CycleConfig cm = LoadText("k3m", ReplaceLine(d2, kUp3,
                                          "up              = [0.0, 0.0, -1.0]"));
   cm.CheckDimension(2);
   TEST_ASSERT(cm.shared.fault_geometry.up[2] == -1.0,
               "K3 up [0, 0, -1]: CheckDimension(2) passes");
}

// ---------------------------------------------------------------------- K4
static void TestSpatialLoader()
{
   if (g_rank == 0) { std::cout << "\n=== K4 the spatial loader alone ===\n"; }
   const std::string deck = ReadFile(kDeck);
   const std::string t = ReplaceLine(ReplaceLine(deck, kRef3, "ref_normal = [3.0, 4.0]"),
                                     kUp3, "up = [0.0, 0.0, 2.0]");
   const spatial::SpatialFrictionConfig s = spatial::ParseSpatialFrictionConfigString(t);
   const spatial::FaultGeometrySpec &g = s.fault_geometry;
   TEST_ASSERT(g.ref_normal_entries == 2 && g.ref_normal[0] == 0.6 && g.ref_normal[1] == 0.8
               && g.ref_normal[2] == 0.0,
               "K4 ref_normal [3, 4]: 2 entries, normalized (0.6, 0.8, 0)");
   TEST_ASSERT(g.up_entries == 3 && g.up[0] == 0.0 && g.up[1] == 0.0 && g.up[2] == 1.0,
               "K4 up [0, 0, 2]: 3 entries, normalized (0, 0, 1)");
   std::string none = ReplaceLine(ReplaceLine(deck, kRef3, ""), kUp3, "");
   const spatial::FaultGeometrySpec g0 = spatial::ParseSpatialFrictionConfigString(none).fault_geometry;
   const spatial::FaultGeometrySpec gd;
   TEST_ASSERT(g0.ref_normal_entries == 0 && g0.up_entries == 0 && g0.ref_normal == gd.ref_normal
               && g0.up == gd.up, "K4 absent keys: entries 0, the defaults (0, -1, 0), (0, 0, 1)");
   spatial::CheckFaultGeometryDimension(g0, 3);
   TEST_ASSERT(true, "K4 absent keys pass CheckFaultGeometryDimension(3)");
}

// ---------------------------------------------------------------- refusals
static int Refuse(int n)
{
   const std::string deck = ReadFile(kDeck);
   const std::string d2 = To2D(deck);
   const std::string st3 = "x    = [0.0, 0.0, -22.0e3]";
   std::string text;
   int dim = 3;
   switch (n)
   {
      case 1: text = ReplaceLine(deck, kRef3, "ref_normal = [0.0, -1.0]"); break;
      case 2: text = ReplaceLine(deck, kUp3, "up = [1.0, 0.0]"); break;
      case 3: text = ReplaceLine(deck, st3, "x    = [0.0, -22.0e3]"); break;
      case 4: text = ReplaceLine(d2, "plane_origin = [0.0, 0.0]", "plane_origin = [0.0, 0.0, 0.0]");
              dim = 2; break;
      case 5: text = ReplaceLine(d2, kUp3, ""); dim = 2; break;
      case 6: text = ReplaceLine(d2, kUp3, "up = [0.0, 0.6, 0.8]"); dim = 2; break;
      case 7: text = ReplaceLine(d2, "ref_normal      = [0.0, -1.0]", kRef3); dim = 2; break;
      case 8: text = ReplaceLine(d2, "ref_normal      = [0.0, -1.0]", ""); dim = 2; break;
      case 9: text = ReplaceLine(deck, kRef3, "ref_normal = [0.0, -1.0, 0.0, 0.0]"); break;
      case 10: text = ReplaceLine(deck, st3, "x    = [0.0, 0.0, -22.0e3, 0.0]"); break;
      default: MFEM_ABORT("--refuse N: N = 1..10");
   }
   if (g_rank == 0) { std::cout << "--refuse " << n << ": a " << dim << "-D check\n"; }
   const cycle::CycleConfig c = LoadText("refuse" + std::to_string(n), text);
   c.CheckDimension(dim);   // must abort
   if (g_rank == 0) { std::cout << "--refuse " << n << ": the deck was accepted\n"; }
   return 0;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   g_rank = Mpi::WorldRank();
   int refuse = 0;
   for (int i = 1; i < argc; i++)
   {
      if (std::strcmp(argv[i], "--refuse") == 0 && i + 1 < argc) { refuse = std::atoi(argv[++i]); }
   }
   if (refuse) { return Refuse(refuse); }
   if (g_rank == 0)
   {
      std::cout << "test_2d_deck (2D in-plane design rev 3, phase B) on " << Mpi::WorldSize()
                << " rank(s)\n";
   }
   TestDeck3D2D();
   TestUpChoices();
   TestSpatialLoader();
   if (g_rank == 0)
   {
      std::cout << "\ntest_2d_deck: " << (num_tests - num_failed) << " / " << num_tests
                << " passed" << (num_failed ? " -- FAILED" : "") << "\n";
   }
   return num_failed ? 1 : 0;
}
