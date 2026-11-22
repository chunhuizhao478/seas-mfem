// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// grid_material.hpp — the 2.5D cycle's gridded velocity and initial stress
// (heterogeneity plan rev4, Steps B1-B3; branch
// fully-dynamic-fault-opening-hetergeneity).
//
// Deck (the shared [velocity] / [stress] sections; spatial_friction.cpp
// parses them, and drivers other than fd_cycle refuse these kinds):
//
//   [velocity]
//   kind     = "grid"          # lambda, mu, rho from the file (fields lambda, mu, rho)
//   file     = "ucerf2d/hetero/ucerf2d_z10km_cvms426_so2.nc"   # build_hetero_grid.py
//   sampling = "quadrature"    # or "cell_average" (SeisSol's homogenization), both regimes
//
//   [stress]
//   kind = "grid_tensor"       # the prestress tensor per fault DOF (s_xx, s_yy, s_xy, s_zz)
//   file = "ucerf2d/hetero/ucerf2d_z10km_cvms426_so2.nc"
//
// (optional in both: outside = "abort" (default) or "clamp")
//
// Without them ([velocity] kind "constant", [stress] another kind) every path
// of the cycle is as before, bit for bit.
//
// "quadrature" (D11): K and M of both regimes read the material at every
// quadrature point (QD: the ElasticityDomainOperator coefficients; dynamic:
// NitscheWaveOperator with MaterialField::pointwise, the rho-weighted
// NitscheMass); the per-element scalars (CFL, penalty scale, Kelvin-Voigt,
// damping rates, Q) take the element's stiffest quadrature point
// (MaterialField::Sampling::StiffestPoint).  "cell_average": one value per
// element by SeisSol's rule (MaterialField::Sampling::CellAverage), used by
// K and M of both regimes and by the per-element scalars.
//
// The file's coordinates are the mesh's ([cycle].mesh_scale must be 1).

#ifndef MFEM_SEAS_CYCLE_GRID_MATERIAL_HPP
#define MFEM_SEAS_CYCLE_GRID_MATERIAL_HPP

#include "mfem.hpp"
#include "../io/grid_field_2d.hpp"
#include "../dynamic/heterogeneous_material.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <sstream>
#include <string>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// One field of a GridField2D as an mfem::Coefficient at physical points.
class GridFieldCoefficient : public Coefficient
{
public:
   GridFieldCoefficient(const GridField2D &g, int field) : g_(g), k_(field) {}
   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
      real_t x[3] = {0.0, 0.0, 0.0};
      Vector X(x, T.GetSpaceDim());
      T.Transform(ip, X);
      return g_.Evaluate(k_, x[0], x[1]);
   }
private:
   const GridField2D &g_;
   int k_;
};

/// FNV-1a 64 of a file's bytes, for the restart key (a checkpoint refuses a
/// deck whose grid file differs).
inline std::string FileHash64(const std::string &path)
{
   FILE *f = std::fopen(path.c_str(), "rb");
   MFEM_VERIFY(f, "FileHash64: cannot open " << path);
   std::uint64_t h = 1469598103934665603ULL;
   unsigned char buf[1 << 16];
   size_t n;
   while ((n = std::fread(buf, 1, sizeof(buf), f)) > 0)
   {
      for (size_t i = 0; i < n; i++) { h ^= buf[i]; h *= 1099511628211ULL; }
   }
   std::fclose(f);
   std::ostringstream s; s << std::hex << h;
   return s.str();
}

/// [velocity] kind = "grid" / [stress] kind = "grid_tensor" and the loaded
/// grids (CycleConfig::Load fills it from the shared config).  Copies share
/// the grids.
struct HeterogeneityConfig
{
   bool velocity = false;
   bool stress = false;
   std::string sampling = "quadrature";
   std::string vel_path, vel_outside = "abort";   ///< [velocity] file, outside
   std::string str_path, str_outside = "abort";   ///< [stress] file, outside
   std::string vel_hash, str_hash;                ///< FileHash64, set by Load

   bool On() const { return velocity || stress; }
   bool Pointwise() const { return sampling == "quadrature"; }

   /// Read the fields the switches need (collective in effect: every rank
   /// reads the file).  Called by CycleConfig::Load.
   void Load()
   {
      MFEM_VERIFY(sampling == "quadrature" || sampling == "cell_average",
                  "[velocity].sampling = \"" << sampling
                  << "\"; use \"quadrature\" or \"cell_average\"");
      auto oob = [](const std::string &o) { return o == "clamp" ? OOBPolicy::Clamp : OOBPolicy::Abort; };
      if (velocity)
      {
         vel_ = std::make_shared<GridField2D>(vel_path, std::vector<std::string>{"lambda", "mu", "rho"},
                                              oob(vel_outside));
         lam_ = std::make_shared<GridFieldCoefficient>(*vel_, 0);
         mu_ = std::make_shared<GridFieldCoefficient>(*vel_, 1);
         rho_ = std::make_shared<GridFieldCoefficient>(*vel_, 2);
         vel_hash = FileHash64(vel_path);
      }
      if (stress)
      {
         str_ = std::make_shared<GridField2D>(str_path, std::vector<std::string>{"s_xx", "s_yy", "s_xy", "s_zz"},
                                              oob(str_outside));
         str_hash = (velocity && str_path == vel_path) ? vel_hash : FileHash64(str_path);
      }
   }

   /// The material of both regimes (velocity = true only).
   MaterialField Material() const
   {
      MFEM_VERIFY(velocity && lam_, "HeterogeneityConfig::Material without [velocity] kind = \"grid\"");
      MaterialField m = MaterialField::MakeCoefficient(lam_.get(), mu_.get(), rho_.get());
      m.pointwise = Pointwise();
      m.sampling = Pointwise() ? MaterialField::Sampling::StiffestPoint
                               : MaterialField::Sampling::CellAverage;
      return m;
   }
   /// lambda, mu, rho at a point (velocity = true only).
   void MaterialAt(const real_t x[3], real_t &lam, real_t &mu, real_t &rho) const
   {
      MFEM_VERIFY(velocity && vel_, "HeterogeneityConfig::MaterialAt without velocity");
      real_t v[3];
      vel_->EvaluateAll(x[0], x[1], v);
      lam = v[0]; mu = v[1]; rho = v[2];
   }
   Coefficient *Lambda() const { return lam_.get(); }
   Coefficient *Mu() const { return mu_.get(); }

   /// The prestress tensor at (x, y), row major, tension positive [Pa]
   /// (stress = true only); s_xz = s_yz = 0 (Andersonian, vertical faults).
   void StressAt(const real_t x[3], real_t s[9]) const
   {
      MFEM_VERIFY(stress && str_, "HeterogeneityConfig::StressAt without [stress] kind = \"grid_tensor\"");
      real_t v[4];
      str_->EvaluateAll(x[0], x[1], v);
      s[0] = v[0]; s[1] = v[2]; s[2] = 0.0;
      s[3] = v[2]; s[4] = v[1]; s[5] = 0.0;
      s[6] = 0.0;  s[7] = 0.0;  s[8] = v[3];
   }

   /// For the restart key and the run log.
   std::string Tag() const
   {
      std::string t;
      if (velocity) { t += "+gridvel:" + sampling + ":" + vel_hash; }
      if (stress) { t += "+gridtensor:" + str_hash; }
      return t;
   }

private:
   std::shared_ptr<GridField2D> vel_, str_;
   std::shared_ptr<GridFieldCoefficient> lam_, mu_, rho_;
};

}  // namespace cycle
}  // namespace seas
}  // namespace mfem

#endif // MFEM_SEAS_CYCLE_GRID_MATERIAL_HPP
