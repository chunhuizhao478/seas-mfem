// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests of the borehole pulse (loading/pulse_load.hpp) against the numbers
// of the spe_hftc inputs and formulation write-up (farms branch cdms).

#include "mfem.hpp"
#include "../../loading/pulse_load.hpp"
#include "../test_framework.hpp"

using namespace mfem;
using namespace mfem::pff;

static PulseLoadParams SpeParams(real_t eta)
{
   PulseLoadParams P;
   P.shape_alpha = 4.1032e5;
   P.shape_beta = 1.1529e7;
   P.rise_time = 3e-7;
   P.pulse_duration = 6e-5;
   P.EM = 0.144;
   P.gap = 7.17e-3;
   P.efficiency = eta;
   P.kappa = 0.35;
   P.m_exp = 0.25;
   P.base_factor = 8000;
   P.center[2] = 49.14e-3 / 3.0;
   P.num_pulses = 20;
   return P;
}

static void TestPeak()
{
   std::cout << "\n[peak pressure at the wall next to the gap, rb = 4.15 mm]\n";
   const real_t rb = 4.15e-3;
   const real_t etas[3] = {0.4, 0.45, 0.5};
   const real_t paper[3] = {206e6, 215e6, 223e6};   // formulation write-up, section 5
   for (int i = 0; i < 3; i++)
   {
      PulseLoad pl(SpeParams(etas[i]));
      const real_t x[3] = {rb, 0.0, pl.Params().center[2]};
      TEST_NEAR(pl.Peak(x), paper[i], 0.6e6, "Pp(eta) matches the write-up to 0.6 MPa");
   }
   PulseLoad pl(SpeParams(0.45));
   const real_t xc[3] = {rb, 0.0, pl.Params().center[2]};
   const real_t xe[3] = {rb, 0.0, pl.Params().center[2] + 7.17e-3/2 + 3e-3};
   TEST_REL(pl.DistanceMM(xe), 1000.0 * std::sqrt(rb*rb + 9e-6), 1e-12,
            "beyond the gap the distance includes the axial offset");
   TEST_ASSERT(pl.Peak(xe) < pl.Peak(xc), "peak decays away from the gap");
}

static void TestShape()
{
   std::cout << "\n[time history]\n";
   PulseLoad pl(SpeParams(0.45));
   TEST_NEAR(pl.Shape(3e-7), 1.0, 1e-12, "shape = 1 at the rise time");
   TEST_NEAR(pl.Shape(0.0), 0.0, 1e-15, "shape = 0 at each onset");
   TEST_NEAR(pl.Shape(6e-5 * 7 + 1.234e-6), pl.Shape(1.234e-6), 1e-9, "periodic with T_w");
   TEST_NEAR(pl.Shape(20 * 6e-5 + 1e-6), 0.0, 0.0, "zero after N_p windows");
   // 10 % of the peak near t_d = 6 us (the defining condition of the rates)
   TEST_NEAR(pl.Shape(6e-6), 0.1, 2e-3, "about 10 % at t_d = 6 us");
   PulseLoadParams P = SpeParams(0.45); P.peak_pressure = 1e8;
   PulseLoad po(P); const real_t x[3] = {1.0, 0.0, 0.0};
   TEST_NEAR(po.Value(3e-7, x), 1e8, 1e-3, "peak_pressure overrides the distance law");
}

static void TestRates()
{
   std::cout << "\n[rates from (t_r, t_d)]\n";
   real_t a, b;
   const bool ok = ShapeRatesFromRiseAndDuration(3e-7, 6e-6, a, b);
   TEST_ASSERT(ok, "root bracketed");
   TEST_REL(a, 4.1032e5, 2e-3, "a = 4.1032e5 1/s (elasticity_absorbing.i, t_dur = 6 us)");
   TEST_REL(b, 1.1529e7, 2e-3, "b = 1.1529e7 1/s");
   ShapeRatesFromRiseAndDuration(3e-7, 3e-6, a, b);
   TEST_REL(a, 8.9489e5, 2e-3, "t_dur = 3 us: a = 8.9489e5");
   TEST_REL(b, 8.3321e6, 2e-3, "t_dur = 3 us: b = 8.3321e6");
}

int main()
{
   TestPeak();
   TestShape();
   TestRates();
   return TestSummary("test_pulse_load");
}
