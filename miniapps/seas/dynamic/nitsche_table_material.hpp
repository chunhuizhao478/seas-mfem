// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_table_material.hpp — a material read from a 3-D table, for the
// Nitsche dynamic-rupture driver (TPV34 plan rev 2, step S2; Chunhui Zhao,
// 2026-11-12).
//
// The table is a data_projection_v1 HDF5 file (document/features_dev/
// data_projection_schema_v1.md) with the fields Vp (m/s), Vs (m/s) and
// density (kg/m3) on one rectilinear grid, as written by
// scripts/cvm/build_cvm_table.py or by the SAFS velocity converter.  It is
// read with DataField3D (trilinear interpolation; with a node communicator
// one copy of the grid per node).
//
//   TableFrame
//       maps a mesh point x to the table point q:
//           q_0 = origin[0] + x_axis[0] x_0 + y_axis[0] x_1
//           q_1 = origin[1] + x_axis[1] x_0 + y_axis[1] x_1
//           q_2 = x_2
//       The default is the identity (the table is in the mesh frame).  A
//       mesh in a fault-local frame reads a table in map coordinates by
//       giving the map position of the mesh origin and the map directions
//       of the mesh X and Y axes.
//   NitscheTableMaterial / MakeNitscheTableMaterial
//       owns the three fields and the coefficients
//           rho(x), mu(x) = rho Vs^2, lambda(x) = rho (Vp^2 - 2 Vs^2)
//       evaluated at q(x), the MaterialField over them, and eval_at_xyz.
//
// The material is smooth in the table: no jump surfaces.  A point outside
// the table aborts (OOBPolicy::Abort) or takes the edge value
// (OOBPolicy::Clamp), as the deck asks.

#ifndef MFEM_SEAS_NITSCHE_TABLE_MATERIAL_HPP
#define MFEM_SEAS_NITSCHE_TABLE_MATERIAL_HPP

#include "mfem.hpp"
#include "heterogeneous_material.hpp"
#include "../io/data_field_3d.hpp"

#include <memory>
#include <string>

namespace mfem
{
namespace seas
{

struct TableFrame
{
   real_t origin[2] = {0.0, 0.0};
   real_t x_axis[2] = {1.0, 0.0};
   real_t y_axis[2] = {0.0, 1.0};

   bool IsIdentity() const
   {
      return origin[0] == 0.0 && origin[1] == 0.0 && x_axis[0] == 1.0 && x_axis[1] == 0.0
             && y_axis[0] == 0.0 && y_axis[1] == 1.0;
   }
   void Map(const real_t x[3], real_t q[3]) const
   {
      q[0] = origin[0] + x_axis[0] * x[0] + y_axis[0] * x[1];
      q[1] = origin[1] + x_axis[1] * x[0] + y_axis[1] * x[1];
      q[2] = x[2];
   }
};

struct NitscheTableMaterial
{
   /// One of rho, mu, lambda at the table point of a mesh point.
   class Coef : public Coefficient
   {
   public:
      Coef(const NitscheTableMaterial &m, int which) : m_(m), which_(which) { }
      real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
      {
         real_t xd[3];
         Vector x(xd, 3);
         T.Transform(ip, x);
         real_t v[3];   // lambda, mu, rho
         m_.eval_at_xyz(xd, v[0], v[1], v[2]);
         return v[which_];
      }
   private:
      const NitscheTableMaterial &m_;
      const int which_;
   };

   std::unique_ptr<DataField3D> vp, vs, rho;
   TableFrame frame;
   std::unique_ptr<Coef> lambda_coef, mu_coef, rho_coef;
   MaterialField field;

   void eval_at_xyz(const real_t x[3], real_t &lambda_out, real_t &mu_out, real_t &rho_out) const
   {
      real_t q[3];
      frame.Map(x, q);
      const real_t r = rho->Evaluate(q[0], q[1], q[2]);
      const real_t a = vp->Evaluate(q[0], q[1], q[2]);
      const real_t b = vs->Evaluate(q[0], q[1], q[2]);
      rho_out = r;
      mu_out = r * b * b;
      lambda_out = r * (a * a - 2.0 * b * b);
   }

   /// True when the table holds the eight corners of the mesh box [lo, hi]
   /// (mapped by the frame), which bounds every mesh point since the map is
   /// affine and the table is a box.
   bool ContainsMeshBox(const real_t lo[3], const real_t hi[3]) const
   {
      const std::array<real_t, 6> &b = vp->BBox();
      for (int c = 0; c < 8; c++)
      {
         const real_t x[3] = {(c & 1) ? hi[0] : lo[0], (c & 2) ? hi[1] : lo[1],
                              (c & 4) ? hi[2] : lo[2]};
         real_t q[3];
         frame.Map(x, q);
         for (int d = 0; d < 3; d++)
         {
            const real_t tol = 1e-9 * (b[2 * d + 1] - b[2 * d]);
            if (q[d] < b[2 * d] - tol || q[d] > b[2 * d + 1] + tol) { return false; }
         }
      }
      return true;
   }

   NitscheTableMaterial() = default;
   NitscheTableMaterial(const NitscheTableMaterial &) = delete;
   NitscheTableMaterial &operator=(const NitscheTableMaterial &) = delete;
};

/// Reads Vp, Vs and density of `path`.  `node_comm` as DataField3D's: a
/// node-local communicator keeps one copy of each grid per node (collective
/// over it); MPI_COMM_NULL keeps one per rank.
inline std::unique_ptr<NitscheTableMaterial>
MakeNitscheTableMaterial(const std::string &path, const TableFrame &frame, bool clamp_outside
#ifdef MFEM_USE_MPI
                         , MPI_Comm node_comm = MPI_COMM_NULL
#endif
                        )
{
   auto m = std::make_unique<NitscheTableMaterial>();
   const OOBPolicy oob = clamp_outside ? OOBPolicy::Clamp : OOBPolicy::Abort;
   const char *names[3] = {"Vp", "Vs", "density"};
   std::unique_ptr<DataField3D> *dst[3] = {&m->vp, &m->vs, &m->rho};
   for (int i = 0; i < 3; i++)
   {
#ifdef MFEM_USE_MPI
      if (node_comm != MPI_COMM_NULL)
      {
         *dst[i] = std::make_unique<DataField3D>(path, names[i], oob, node_comm);
         continue;
      }
#endif
      *dst[i] = std::make_unique<DataField3D>(path, names[i], oob);
   }
   m->frame = frame;
   m->lambda_coef = std::make_unique<NitscheTableMaterial::Coef>(*m, 0);
   m->mu_coef = std::make_unique<NitscheTableMaterial::Coef>(*m, 1);
   m->rho_coef = std::make_unique<NitscheTableMaterial::Coef>(*m, 2);
   m->field = MaterialField::MakeCoefficient(m->lambda_coef.get(), m->mu_coef.get(),
                                             m->rho_coef.get());
   return m;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_TABLE_MATERIAL_HPP
