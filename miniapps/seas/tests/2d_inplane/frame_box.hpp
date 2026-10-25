// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/frame_box.hpp — shared by the 2D in-plane branch's tests
// (test_fault_frame_field.cpp, test_plate_loading.cpp): rigid rotations and
// the test box of 4 x 4 x 2 hexes split into 24 tets each, with a fault
// patch (attr 3) and a coplanar prescribed-slip band (attr 6) on y = 0, walls
// and bottom attr 5, top attr 1 (kind 0), or the L-shaped seam surface
// (kind 1).  Vertices rotated by R after the faces are tagged.
//
// Scope (Chunhui Zhao, 2026-10-25).

#ifndef MFEM_SEAS_TESTS_2D_INPLANE_FRAME_BOX_HPP
#define MFEM_SEAS_TESTS_2D_INPLANE_FRAME_BOX_HPP

#include "mfem.hpp"

#include <cmath>
#include <vector>

namespace frame_box
{
using mfem::real_t;
using mfem::Mesh;
using mfem::Array;

static const real_t kL = 10.0e3;   ///< box half-width [m]

// ---------------------------------------------------------------- rotations
struct Rot
{
   real_t R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
   void Apply(const real_t a[3], real_t b[3]) const
   {
      for (int i = 0; i < 3; i++) { b[i] = R[i][0] * a[0] + R[i][1] * a[1] + R[i][2] * a[2]; }
   }
   // Rodrigues: angle (degrees) about axis.
   static Rot Axis(real_t ax, real_t ay, real_t az, real_t deg)
   {
      const real_t l = std::sqrt(ax * ax + ay * ay + az * az);
      const real_t k[3] = {ax / l, ay / l, az / l};
      const real_t t = deg * M_PI / 180.0, c = std::cos(t), s = std::sin(t);
      Rot r;
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++)
         {
            r.R[i][j] = (i == j ? c : 0.0) + (1 - c) * k[i] * k[j];
         }
      r.R[0][1] -= s * k[2]; r.R[0][2] += s * k[1];
      r.R[1][0] += s * k[2]; r.R[1][2] -= s * k[0];
      r.R[2][0] -= s * k[1]; r.R[2][1] += s * k[0];
      return r;
   }
};

// ------------------------------------------------------------------- meshes
// kind 0: the coplanar box (fault + band on y = 0); kind 1: the L-shaped
// seam surface.  Vertices rotated by R after the faces are tagged.
inline Mesh *FrameBox(int kind, const Rot &rot)
{
   const int nx = 4, ny = 4, nz = (kind == 0) ? 2 : 1;
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(nx, ny, nz, 2 * kL, 2 * kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
   }
   const real_t tol = 1e-6 * kL;
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      real_t zc = 0.0;
      for (int i = 0; i < verts.Size(); i++) { zc += mesh->GetVertex(verts[i])[2]; }
      zc /= verts.Size();
      mesh->SetBdrAttribute(b, (std::abs(zc) < tol) ? 1 : 5);
   }
   struct Tri { int v[3]; int attr; };
   std::vector<Tri> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      real_t c[3] = {0, 0, 0};
      bool y0 = true, x0 = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         const real_t *x = mesh->GetVertex(fv[i]);
         for (int d = 0; d < 3; d++) { c[d] += x[d] / fv.Size(); }
         if (std::abs(x[1]) > tol) { y0 = false; }
         if (std::abs(x[0]) > tol) { x0 = false; }
      }
      int attr = 0;
      if (kind == 0 && y0)
      {
         attr = (std::abs(c[0]) < 0.5 * kL && c[2] > -0.5 * kL) ? 3 : 6;
      }
      else if (kind == 1)
      {
         if (y0 && c[0] > -0.5 * kL && c[0] < 0.0) { attr = 3; }
         if (x0 && c[1] > 0.0 && c[1] < 0.5 * kL) { attr = 6; }
      }
      if (attr) { tris.push_back({{fv[0], fv[1], fv[2]}, attr}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t.v[0], t.v[1], t.v[2], t.attr); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      const real_t x[3] = {c[0], c[1], c[2]};
      rot.Apply(x, c);
   }
   return mesh;
}


} // namespace frame_box

#endif // MFEM_SEAS_TESTS_2D_INPLANE_FRAME_BOX_HPP
