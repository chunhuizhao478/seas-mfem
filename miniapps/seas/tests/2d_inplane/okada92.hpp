// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/okada92.hpp — Okada (1992), "Internal deformation due to
// shear and tensile faults in a half-space", BSSA 82(2), 1018-1040: the
// displacement and its gradient at a point of an elastic half-space z <= 0
// due to uniform slip on a buried rectangular fault.  A line-by-line C++
// transcription of Okada's subroutine DC3D (coded by Y. Okada, Sep. 1991,
// revised May 2002; with UA, UB, UC, DCCON0, DCCON2), in double precision
// throughout (the distributed Fortran takes REAL*4 arguments).  Source: the
// DC3D.f that NIED distributes (www.bosai.go.jp/study/application/dc3d/,
// no licence stated there), as packaged in okada_wrapper 24.6.15 (PyPI,
// T. Ben Thompson, MIT licence for the wrapper).  Used by
// test_okada_dislocation.cpp, which also checks the transcription's physics
// (jump across the patch, free surface, equilibrium, traction continuity).
//
// Okada's coordinates: x along strike, z up (z <= 0 in the medium); the
// fault plane passes through (0, 0, -depth) and dips by `dip` degrees, its
// points are (xi, eta cos(dip), -depth + eta sin(dip)) with
// al1 <= xi <= al2 and aw1 <= eta <= aw2 (eta up-dip).  disl1, disl2,
// disl3: strike-, dip- and tensile dislocation, the motion of the hanging
// wall (the block above the plane) relative to the foot wall.  alpha =
// (lambda + mu) / (lambda + 2 mu).  Lengths in any one unit (Okada's
// singular-point tolerance is 1e-6 of it; the test uses km), displacements
// in the unit of the dislocation, gradients in dislocation unit per length.
//
// Scope (Chunhui Zhao, 2026-10-25).

#ifndef MFEM_SEAS_TESTS_2D_INPLANE_OKADA92_HPP
#define MFEM_SEAS_TESTS_2D_INPLANE_OKADA92_HPP

#include <cmath>

namespace okada92
{

/// DC3D's output: u[i] = u_i, du[j][i] = d u_i / d x_j (Okada's UXX, UYX,
/// UZX = du[0][0..2], UXY, UYY, UZY = du[1][0..2], UXZ, UYZ, UZZ =
/// du[2][0..2]); iret 0 normal, 1 singular (on a fault edge), 2 z > 0.
struct Result
{
   double u[3] = {0, 0, 0};
   double du[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
   int iret = 0;
};

namespace detail
{
static const double kF0 = 0.0, kF1 = 1.0, kF2 = 2.0, kF3 = 3.0;
static const double kPI2 = 6.283185307179586;
static const double kEPS = 1.0e-6;

// COMMON /C0/
struct C0 { double alp1, alp2, alp3, alp4, alp5, sd, cd, sdsd, cdcd, sdcd, s2d, c2d; };
// COMMON /C2/
struct C2
{
   double xi2, et2, q2, r, r2, r3, r5, y, d, tt, alx, ale, x11, y11, x32, y32,
          ey, ez, fy, fz, gy, gz, hy, hz;
};

inline void DCCON0(double alpha, double dip, C0 &c)
{
   c.alp1 = (kF1 - alpha) / kF2;
   c.alp2 = alpha / kF2;
   c.alp3 = (kF1 - alpha) / alpha;
   c.alp4 = kF1 - alpha;
   c.alp5 = alpha;
   const double p18 = kPI2 / 360.0;
   c.sd = std::sin(dip * p18);
   c.cd = std::cos(dip * p18);
   if (std::abs(c.cd) < kEPS)
   {
      c.cd = kF0;
      if (c.sd > kF0) { c.sd = kF1; }
      if (c.sd < kF0) { c.sd = -kF1; }
   }
   c.sdsd = c.sd * c.sd;
   c.cdcd = c.cd * c.cd;
   c.sdcd = c.sd * c.cd;
   c.s2d = kF2 * c.sdcd;
   c.c2d = c.cdcd - c.sdsd;
}

inline void DCCON2(double &xi, double &et, double &q, double sd, double cd, int kxi, int ket,
                   C2 &c)
{
   if (std::abs(xi) < kEPS) { xi = kF0; }
   if (std::abs(et) < kEPS) { et = kF0; }
   if (std::abs(q) < kEPS) { q = kF0; }
   c.xi2 = xi * xi;
   c.et2 = et * et;
   c.q2 = q * q;
   c.r2 = c.xi2 + c.et2 + c.q2;
   c.r = std::sqrt(c.r2);
   if (c.r == kF0) { return; }
   c.r3 = c.r * c.r2;
   c.r5 = c.r3 * c.r2;
   c.y = et * cd + q * sd;
   c.d = et * sd - q * cd;
   if (q == kF0) { c.tt = kF0; }
   else          { c.tt = std::atan(xi * et / (q * c.r)); }
   if (kxi == 1)
   {
      c.alx = -std::log(c.r - xi);
      c.x11 = kF0;
      c.x32 = kF0;
   }
   else
   {
      const double rxi = c.r + xi;
      c.alx = std::log(rxi);
      c.x11 = kF1 / (c.r * rxi);
      c.x32 = (c.r + rxi) * c.x11 * c.x11 / c.r;
   }
   if (ket == 1)
   {
      c.ale = -std::log(c.r - et);
      c.y11 = kF0;
      c.y32 = kF0;
   }
   else
   {
      const double ret = c.r + et;
      c.ale = std::log(ret);
      c.y11 = kF1 / (c.r * ret);
      c.y32 = (c.r + ret) * c.y11 * c.y11 / c.r;
   }
   c.ey = sd / c.r - c.y * q / c.r3;
   c.ez = cd / c.r + c.d * q / c.r3;
   c.fy = c.d / c.r3 + c.xi2 * c.y32 * sd;
   c.fz = c.y / c.r3 + c.xi2 * c.y32 * cd;
   c.gy = kF2 * c.x11 * sd - c.y * q * c.x32;
   c.gz = kF2 * c.x11 * cd + c.d * q * c.x32;
   c.hy = c.d * q * c.x32 + xi * q * c.y32 * sd;
   c.hz = c.y * q * c.x32 + xi * q * c.y32 * cd;
}

// Part A (u[0..11] = Okada's U(1..12)).
inline void UA(double xi, double et, double q, double disl1, double disl2, double disl3,
               const C0 &a, const C2 &c, double u[12])
{
   const double alp1 = a.alp1, alp2 = a.alp2, sd = a.sd, cd = a.cd;
   const double xi2 = c.xi2, q2 = c.q2, r = c.r, r3 = c.r3, y = c.y, d = c.d, tt = c.tt,
                alx = c.alx, ale = c.ale, x11 = c.x11, y11 = c.y11, y32 = c.y32, ey = c.ey,
                ez = c.ez, fy = c.fy, fz = c.fz, gy = c.gy, gz = c.gz, hy = c.hy, hz = c.hz;
   double du[12];
   for (int i = 0; i < 12; i++) { u[i] = kF0; }
   const double xy = xi * y11, qx = q * x11, qy = q * y11;
   if (disl1 != kF0)
   {
      du[0] = tt / kF2 + alp2 * xi * qy;
      du[1] = alp2 * q / r;
      du[2] = alp1 * ale - alp2 * q * qy;
      du[3] = -alp1 * qy - alp2 * xi2 * q * y32;
      du[4] = -alp2 * xi * q / r3;
      du[5] = alp1 * xy + alp2 * xi * q2 * y32;
      du[6] = alp1 * xy * sd + alp2 * xi * fy + d / kF2 * x11;
      du[7] = alp2 * ey;
      du[8] = alp1 * (cd / r + qy * sd) - alp2 * q * fy;
      du[9] = alp1 * xy * cd + alp2 * xi * fz + y / kF2 * x11;
      du[10] = alp2 * ez;
      du[11] = -alp1 * (sd / r - qy * cd) - alp2 * q * fz;
      for (int i = 0; i < 12; i++) { u[i] += disl1 / kPI2 * du[i]; }
   }
   if (disl2 != kF0)
   {
      du[0] = alp2 * q / r;
      du[1] = tt / kF2 + alp2 * et * qx;
      du[2] = alp1 * alx - alp2 * q * qx;
      du[3] = -alp2 * xi * q / r3;
      du[4] = -qy / kF2 - alp2 * et * q / r3;
      du[5] = alp1 / r + alp2 * q2 / r3;
      du[6] = alp2 * ey;
      du[7] = alp1 * d * x11 + xy / kF2 * sd + alp2 * et * gy;
      du[8] = alp1 * y * x11 - alp2 * q * gy;
      du[9] = alp2 * ez;
      du[10] = alp1 * y * x11 + xy / kF2 * cd + alp2 * et * gz;
      du[11] = -alp1 * d * x11 - alp2 * q * gz;
      for (int i = 0; i < 12; i++) { u[i] += disl2 / kPI2 * du[i]; }
   }
   if (disl3 != kF0)
   {
      du[0] = -alp1 * ale - alp2 * q * qy;
      du[1] = -alp1 * alx - alp2 * q * qx;
      du[2] = tt / kF2 - alp2 * (et * qx + xi * qy);
      du[3] = -alp1 * xy + alp2 * xi * q2 * y32;
      du[4] = -alp1 / r + alp2 * q2 / r3;
      du[5] = -alp1 * qy - alp2 * q * q2 * y32;
      du[6] = -alp1 * (cd / r + qy * sd) - alp2 * q * fy;
      du[7] = -alp1 * y * x11 - alp2 * q * gy;
      du[8] = alp1 * (d * x11 + xy * sd) + alp2 * q * hy;
      du[9] = alp1 * (sd / r - qy * cd) - alp2 * q * fz;
      du[10] = alp1 * d * x11 - alp2 * q * gz;
      du[11] = alp1 * (y * x11 + xy * cd) + alp2 * q * hz;
      for (int i = 0; i < 12; i++) { u[i] += disl3 / kPI2 * du[i]; }
   }
}

// Part B.
inline void UB(double xi, double et, double q, double disl1, double disl2, double disl3,
               const C0 &a, const C2 &c, double u[12])
{
   const double alp3 = a.alp3, sd = a.sd, cd = a.cd, sdsd = a.sdsd, cdcd = a.cdcd,
                sdcd = a.sdcd;
   const double xi2 = c.xi2, q2 = c.q2, r = c.r, r3 = c.r3, y = c.y, d = c.d, tt = c.tt,
                ale = c.ale, x11 = c.x11, y11 = c.y11, y32 = c.y32, ey = c.ey, ez = c.ez,
                fy = c.fy, fz = c.fz, gy = c.gy, gz = c.gz, hy = c.hy, hz = c.hz;
   double du[12];
   const double rd = r + d;
   const double d11 = kF1 / (r * rd);
   const double aj2 = xi * y / rd * d11;
   const double aj5 = -(d + y * y / rd) * d11;
   double ai3, ai4, ak1, ak3, aj3, aj6;
   if (cd != kF0)
   {
      if (xi == kF0)
      {
         ai4 = kF0;
      }
      else
      {
         const double x = std::sqrt(xi2 + q2);
         ai4 = kF1 / cdcd * (xi / rd * sdcd +
                             kF2 * std::atan((et * (x + q * cd) + x * (r + x) * sd) /
                                             (xi * (r + x) * cd)));
      }
      ai3 = (y * cd / rd - ale + sd * std::log(rd)) / cdcd;
      ak1 = xi * (d11 - y11 * sd) / cd;
      ak3 = (q * y11 - y * d11) / cd;
      aj3 = (ak1 - aj2 * sd) / cd;
      aj6 = (ak3 - aj5 * sd) / cd;
   }
   else
   {
      const double rd2 = rd * rd;
      ai3 = (et / rd + y * q / rd2 - ale) / kF2;
      ai4 = xi * y / rd2 / kF2;
      ak1 = xi * q / rd * d11;
      ak3 = sd / rd * (xi2 * d11 - kF1);
      aj3 = -xi / rd2 * (q2 * d11 - kF1 / kF2);
      aj6 = -y / rd2 * (xi2 * d11 - kF1 / kF2);
   }
   const double xy = xi * y11;
   const double ai1 = -xi / rd * cd - ai4 * sd;
   const double ai2 = std::log(rd) + ai3 * sd;
   const double ak2 = kF1 / r + ak3 * sd;
   const double ak4 = xy * cd - ak1 * sd;
   const double aj1 = aj5 * cd - aj6 * sd;
   const double aj4 = -xy - aj2 * cd + aj3 * sd;
   for (int i = 0; i < 12; i++) { u[i] = kF0; }
   const double qx = q * x11, qy = q * y11;
   if (disl1 != kF0)
   {
      du[0] = -xi * qy - tt - alp3 * ai1 * sd;
      du[1] = -q / r + alp3 * y / rd * sd;
      du[2] = q * qy - alp3 * ai2 * sd;
      du[3] = xi2 * q * y32 - alp3 * aj1 * sd;
      du[4] = xi * q / r3 - alp3 * aj2 * sd;
      du[5] = -xi * q2 * y32 - alp3 * aj3 * sd;
      du[6] = -xi * fy - d * x11 + alp3 * (xy + aj4) * sd;
      du[7] = -ey + alp3 * (kF1 / r + aj5) * sd;
      du[8] = q * fy - alp3 * (qy - aj6) * sd;
      du[9] = -xi * fz - y * x11 + alp3 * ak1 * sd;
      du[10] = -ez + alp3 * y * d11 * sd;
      du[11] = q * fz + alp3 * ak2 * sd;
      for (int i = 0; i < 12; i++) { u[i] += disl1 / kPI2 * du[i]; }
   }
   if (disl2 != kF0)
   {
      du[0] = -q / r + alp3 * ai3 * sdcd;
      du[1] = -et * qx - tt - alp3 * xi / rd * sdcd;
      du[2] = q * qx + alp3 * ai4 * sdcd;
      du[3] = xi * q / r3 + alp3 * aj4 * sdcd;
      du[4] = et * q / r3 + qy + alp3 * aj5 * sdcd;
      du[5] = -q2 / r3 + alp3 * aj6 * sdcd;
      du[6] = -ey + alp3 * aj1 * sdcd;
      du[7] = -et * gy - xy * sd + alp3 * aj2 * sdcd;
      du[8] = q * gy + alp3 * aj3 * sdcd;
      du[9] = -ez - alp3 * ak3 * sdcd;
      du[10] = -et * gz - xy * cd - alp3 * xi * d11 * sdcd;
      du[11] = q * gz - alp3 * ak4 * sdcd;
      for (int i = 0; i < 12; i++) { u[i] += disl2 / kPI2 * du[i]; }
   }
   if (disl3 != kF0)
   {
      du[0] = q * qy - alp3 * ai3 * sdsd;
      du[1] = q * qx + alp3 * xi / rd * sdsd;
      du[2] = et * qx + xi * qy - tt - alp3 * ai4 * sdsd;
      du[3] = -xi * q2 * y32 - alp3 * aj4 * sdsd;
      du[4] = -q2 / r3 - alp3 * aj5 * sdsd;
      du[5] = q * q2 * y32 - alp3 * aj6 * sdsd;
      du[6] = q * fy - alp3 * aj1 * sdsd;
      du[7] = q * gy - alp3 * aj2 * sdsd;
      du[8] = -q * hy - alp3 * aj3 * sdsd;
      du[9] = q * fz + alp3 * ak3 * sdsd;
      du[10] = q * gz + alp3 * xi * d11 * sdsd;
      du[11] = -q * hz + alp3 * ak4 * sdsd;
      for (int i = 0; i < 12; i++) { u[i] += disl3 / kPI2 * du[i]; }
   }
}

// Part C.
inline void UC(double xi, double et, double q, double z, double disl1, double disl2,
               double disl3, const C0 &a, const C2 &c, double u[12])
{
   const double alp4 = a.alp4, alp5 = a.alp5, sd = a.sd, cd = a.cd, sdsd = a.sdsd,
                cdcd = a.cdcd, sdcd = a.sdcd;
   const double xi2 = c.xi2, et2 = c.et2, q2 = c.q2, r = c.r, r2 = c.r2, r3 = c.r3,
                r5 = c.r5, y = c.y, d = c.d, x11 = c.x11, y11 = c.y11, x32 = c.x32,
                y32 = c.y32;
   double du[12];
   const double cc = d + z;
   const double x53 = (8.0 * r2 + 9.0 * r * xi + kF3 * xi2) * x11 * x11 * x11 / r2;
   const double y53 = (8.0 * r2 + 9.0 * r * et + kF3 * et2) * y11 * y11 * y11 / r2;
   const double h = q * cd - z;
   const double z32 = sd / r3 - h * y32;
   const double z53 = kF3 * sd / r5 - h * y53;
   const double y0 = y11 - xi2 * y32;
   const double z0 = z32 - xi2 * z53;
   const double ppy = cd / r3 + q * y32 * sd;
   const double ppz = sd / r3 - q * y32 * cd;
   const double qq = z * y32 + z32 + z0;
   const double qqy = kF3 * cc * d / r5 - qq * sd;
   const double qqz = kF3 * cc * y / r5 - qq * cd + q * y32;
   const double xy = xi * y11;
   const double qx = q * x11;
   const double qy = q * y11;
   const double qr = kF3 * q / r5;
   const double cqx = cc * q * x53;
   const double cdr = (cc + d) / r3;
   const double yy0 = y / r3 - y0 * cd;
   (void)qx; (void)cqx;   // set by Okada, unused in his formulas
   for (int i = 0; i < 12; i++) { u[i] = kF0; }
   if (disl1 != kF0)
   {
      du[0] = alp4 * xy * cd - alp5 * xi * q * z32;
      du[1] = alp4 * (cd / r + kF2 * qy * sd) - alp5 * cc * q / r3;
      du[2] = alp4 * qy * cd - alp5 * (cc * et / r3 - z * y11 + xi2 * z32);
      du[3] = alp4 * y0 * cd - alp5 * q * z0;
      du[4] = -alp4 * xi * (cd / r3 + kF2 * q * y32 * sd) + alp5 * cc * xi * qr;
      du[5] = -alp4 * xi * q * y32 * cd + alp5 * xi * (kF3 * cc * et / r5 - qq);
      du[6] = -alp4 * xi * ppy * cd - alp5 * xi * qqy;
      du[7] = alp4 * kF2 * (d / r3 - y0 * sd) * sd - y / r3 * cd -
              alp5 * (cdr * sd - et / r3 - cc * y * qr);
      du[8] = -alp4 * q / r3 + yy0 * sd + alp5 * (cdr * cd + cc * d * qr - (y0 * cd + q * z0) * sd);
      du[9] = alp4 * xi * ppz * cd - alp5 * xi * qqz;
      du[10] = alp4 * kF2 * (y / r3 - y0 * cd) * sd + d / r3 * cd - alp5 * (cdr * cd + cc * d * qr);
      du[11] = yy0 * cd - alp5 * (cdr * sd - cc * y * qr - y0 * sdsd + q * z0 * cd);
      for (int i = 0; i < 12; i++) { u[i] += disl1 / kPI2 * du[i]; }
   }
   if (disl2 != kF0)
   {
      du[0] = alp4 * cd / r - qy * sd - alp5 * cc * q / r3;
      du[1] = alp4 * y * x11 - alp5 * cc * et * q * x32;
      du[2] = -d * x11 - xy * sd - alp5 * cc * (x11 - q2 * x32);
      du[3] = -alp4 * xi / r3 * cd + alp5 * cc * xi * qr + xi * q * y32 * sd;
      du[4] = -alp4 * y / r3 + alp5 * cc * et * qr;
      du[5] = d / r3 - y0 * sd + alp5 * cc / r3 * (kF1 - kF3 * q2 / r2);
      du[6] = -alp4 * et / r3 + y0 * sdsd - alp5 * (cdr * sd - cc * y * qr);
      du[7] = alp4 * (x11 - y * y * x32) - alp5 * cc * ((d + kF2 * q * cd) * x32 - y * et * q * x53);
      du[8] = xi * ppy * sd + y * d * x32 + alp5 * cc * ((y + kF2 * q * sd) * x32 - y * q2 * x53);
      du[9] = -q / r3 + y0 * sdcd - alp5 * (cdr * cd + cc * d * qr);
      du[10] = alp4 * y * d * x32 - alp5 * cc * ((y - kF2 * q * sd) * x32 + d * et * q * x53);
      du[11] = -xi * ppz * sd + x11 - d * d * x32 - alp5 * cc * ((d - kF2 * q * cd) * x32 - d * q2 * x53);
      for (int i = 0; i < 12; i++) { u[i] += disl2 / kPI2 * du[i]; }
   }
   if (disl3 != kF0)
   {
      du[0] = -alp4 * (sd / r + qy * cd) - alp5 * (z * y11 - q2 * z32);
      du[1] = alp4 * kF2 * xy * sd + d * x11 - alp5 * cc * (x11 - q2 * x32);
      du[2] = alp4 * (y * x11 + xy * cd) + alp5 * q * (cc * et * x32 + xi * z32);
      du[3] = alp4 * xi / r3 * sd + xi * q * y32 * cd + alp5 * xi * (kF3 * cc * et / r5 - kF2 * z32 - z0);
      du[4] = alp4 * kF2 * y0 * sd - d / r3 + alp5 * cc / r3 * (kF1 - kF3 * q2 / r2);
      du[5] = -alp4 * yy0 - alp5 * (cc * et * qr - q * z0);
      du[6] = alp4 * (q / r3 + y0 * sdcd) + alp5 * (z / r3 * cd + cc * d * qr - q * z0 * sd);
      du[7] = -alp4 * kF2 * xi * ppy * sd - y * d * x32 +
              alp5 * cc * ((y + kF2 * q * sd) * x32 - y * q2 * x53);
      du[8] = -alp4 * (xi * ppy * cd - x11 + y * y * x32) +
              alp5 * (cc * ((d + kF2 * q * cd) * x32 - y * et * q * x53) + xi * qqy);
      du[9] = -et / r3 + y0 * cdcd - alp5 * (z / r3 * sd - cc * y * qr - y0 * sdsd + q * z0 * cd);
      du[10] = alp4 * kF2 * xi * ppz * sd - x11 + d * d * x32 -
               alp5 * cc * ((d - kF2 * q * cd) * x32 - d * q2 * x53);
      du[11] = alp4 * (xi * ppz * cd + y * d * x32) +
               alp5 * (cc * ((y - kF2 * q * sd) * x32 + d * et * q * x53) + xi * qqz);
      for (int i = 0; i < 12; i++) { u[i] += disl3 / kPI2 * du[i]; }
   }
}

// One source (real: image = false, d = depth + z; image: d = depth - z):
// the edge rejection, the negative-extension flags and the four corners.
inline bool Source(bool image, double y, double zz, double dd, const double xi_in[2],
                   double aw1, double aw2, double dd1, double dd2, double dd3, const C0 &a,
                   double U[12])
{
   const double sd = a.sd, cd = a.cd;
   double xi[2] = {xi_in[0], xi_in[1]};
   const double p = y * cd + dd * sd;
   double q = y * sd - dd * cd;
   double et[2] = {p - aw1, p - aw2};
   if (std::abs(q) < kEPS) { q = kF0; }
   if (std::abs(et[0]) < kEPS) { et[0] = kF0; }
   if (std::abs(et[1]) < kEPS) { et[1] = kF0; }
   // on a fault edge
   if (q == kF0 && ((xi[0] * xi[1] <= kF0 && et[0] * et[1] == kF0) ||
                    (et[0] * et[1] <= kF0 && xi[0] * xi[1] == kF0)))
   {
      return false;
   }
   // on the negative extension of a fault edge
   int kxi[2] = {0, 0}, ket[2] = {0, 0};
   const double r12 = std::sqrt(xi[0] * xi[0] + et[1] * et[1] + q * q);
   const double r21 = std::sqrt(xi[1] * xi[1] + et[0] * et[0] + q * q);
   const double r22 = std::sqrt(xi[1] * xi[1] + et[1] * et[1] + q * q);
   if (xi[0] < kF0 && r21 + xi[1] < kEPS) { kxi[0] = 1; }
   if (xi[0] < kF0 && r22 + xi[1] < kEPS) { kxi[1] = 1; }
   if (et[0] < kF0 && r12 + et[1] < kEPS) { ket[0] = 1; }
   if (et[0] < kF0 && r22 + et[1] < kEPS) { ket[1] = 1; }
   double du[12], dua[12], dub[12], duc[12];
   for (int k = 0; k < 2; k++)
   {
      for (int j = 0; j < 2; j++)
      {
         C2 c;
         // Okada: CALL DCCON2(XI(J),ET(K),Q,SD,CD,KXI(K),KET(J))
         DCCON2(xi[j], et[k], q, sd, cd, kxi[k], ket[j], c);
         UA(xi[j], et[k], q, dd1, dd2, dd3, a, c, dua);
         if (!image)
         {
            for (int i = 0; i < 12; i += 3)
            {
               du[i] = -dua[i];
               du[i + 1] = -dua[i + 1] * cd + dua[i + 2] * sd;
               du[i + 2] = -dua[i + 1] * sd - dua[i + 2] * cd;
               if (i < 9) { continue; }
               du[i] = -du[i];
               du[i + 1] = -du[i + 1];
               du[i + 2] = -du[i + 2];
            }
         }
         else
         {
            UB(xi[j], et[k], q, dd1, dd2, dd3, a, c, dub);
            UC(xi[j], et[k], q, zz, dd1, dd2, dd3, a, c, duc);
            for (int i = 0; i < 12; i += 3)
            {
               du[i] = dua[i] + dub[i] + zz * duc[i];
               du[i + 1] = (dua[i + 1] + dub[i + 1] + zz * duc[i + 1]) * cd -
                           (dua[i + 2] + dub[i + 2] + zz * duc[i + 2]) * sd;
               du[i + 2] = (dua[i + 1] + dub[i + 1] - zz * duc[i + 1]) * sd +
                           (dua[i + 2] + dub[i + 2] - zz * duc[i + 2]) * cd;
               if (i < 9) { continue; }
               du[9] = du[9] + duc[0];
               du[10] = du[10] + duc[1] * cd - duc[2] * sd;
               du[11] = du[11] - duc[1] * sd - duc[2] * cd;
            }
         }
         // Okada: IF(J+K.NE.3) add, IF(J+K.EQ.3) subtract (1-based J, K)
         for (int i = 0; i < 12; i++)
         {
            if (j + k != 1) { U[i] += du[i]; }
            else            { U[i] -= du[i]; }
         }
      }
   }
   return true;
}
} // namespace detail

/// Okada's DC3D (see the file comment).
inline Result DC3D(double alpha, double x, double y, double z, double depth, double dip,
                   double al1, double al2, double aw1, double aw2,
                   double disl1, double disl2, double disl3)
{
   using namespace detail;
   Result res;
   if (z > 0.0) { res.iret = 2; return res; }
   C0 a;
   DCCON0(alpha, dip, a);
   double xi[2] = {x - al1, x - al2};
   if (std::abs(xi[0]) < kEPS) { xi[0] = kF0; }
   if (std::abs(xi[1]) < kEPS) { xi[1] = kF0; }
   double U[12];
   for (int i = 0; i < 12; i++) { U[i] = kF0; }
   if (!Source(false, y, z, depth + z, xi, aw1, aw2, disl1, disl2, disl3, a, U) ||
       !Source(true, y, z, depth - z, xi, aw1, aw2, disl1, disl2, disl3, a, U))
   {
      res.iret = 1;
      return res;
   }
   for (int i = 0; i < 3; i++) { res.u[i] = U[i]; }
   for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++) { res.du[j][i] = U[3 + 3 * j + i]; }
   return res;
}

} // namespace okada92

#endif // MFEM_SEAS_TESTS_2D_INPLANE_OKADA92_HPP
