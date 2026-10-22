// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_source.hpp — poroelastic option, step P2.3 (plan rev 3).
//
// PoroSource is the fluid source of derivation rev 1.2, eq. (4.3), from the
// deck table [multiphysics.poroelastic.source] (PoroSourceConfig):
//
//   Q(x, t) = rate g(x) f(t),    V(x, t) = int_0^t Q dt = rate g(x) F(t)
//
//   shape = "gaussian"   g = exp(-width |x - x0|^2), rate = peak Q in 1/s
//   shape = "well"       g = delta(x - x_w), x_w = origin_m, rate = Q_w in m^3/s
//
// The time functions (SourceTimeFunction), with s = t / rise_time_s and
// T = rise_time_s; f = F = 0 for t < 0:
//
//   step          f = 1                         F = t
//   ramp          f = s            (s <= 1)     F = T s^2 / 2
//   smooth_step   f = 10 s^3 - 15 s^4 + 6 s^5   F = T s^4 (5/2 - 3 s + s^2)
//   ramp and smooth_step: f = 1 and F = t - T/2 for s > 1.
// Plan rev 3 asks only that each F be closed form; the quintic smooth step
// (zero slope and curvature at both ends) was chosen by Chunhui Zhao on
// 2026-10-21 in the review of step P2.3, over the cubic and over the
// TPV102/104 nucleation curve, whose integral is not closed form.
//
// Q is separable, so V at every quadrature point is rate g(x) times the one
// time factor F; PoroSource keeps that factor (VolumeFactor) for the time it
// was last advanced to.  Advance(t) moves it, by the closed form F(t) (the
// default) or by the trapezoid rule V^{n+1} = V^n + dt/2 (Q^n + Q^{n+1}),
// step 2b of plan section 5.2.
//
// Loads, eq. (6.6):
//   gaussian  InjectedVolume() is the coefficient rate g(x) VolumeFactor() for
//             PoroNitscheOperator, which assembles F_V with its face terms;
//   well      inside one element K, the load is M V_w(t) zeta(v, q)(x_w) =
//             M V_w (alpha div v + div q)(x_w) on K's rows only, no face terms.
//             WellLoad() is that vector for V_w = 1, [solid rows; fluid rows]
//             in the operator's layout; AddWellLoad adds rate F times it.
//             A well within 1e-10 (barycentric) of a face, edge or vertex of
//             an element, or outside the mesh, aborts: the point load would
//             belong to two elements, or to none.
// InitialPressure(ip) is p_init = magnitude exp(-width |x - x0|^2), the
// starting pressure of section 4.2, for PoroNitscheOperator::SetInitialPressure.
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only.

#ifndef MFEM_SEAS_PORO_SOURCE_HPP
#define MFEM_SEAS_PORO_SOURCE_HPP

#include "mfem.hpp"
#include "poro_config.hpp"
#include "poro_material.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{
namespace poro
{

/// f(t) and its closed-form integral F(t) for one deck time function.
class SourceTimeFunction
{
public:
   enum class Kind { kStep, kRamp, kSmoothStep };

   /// name = "step", "ramp" or "smooth_step"; rise > 0 for the last two.
   SourceTimeFunction(const std::string &name, real_t rise)
      : rise_(rise)
   {
      if (name == "step") { kind_ = Kind::kStep; }
      else if (name == "ramp") { kind_ = Kind::kRamp; }
      else if (name == "smooth_step") { kind_ = Kind::kSmoothStep; }
      else { MFEM_ABORT("SourceTimeFunction: unknown time function '" << name << "'."); }
      MFEM_VERIFY(kind_ == Kind::kStep || rise > 0.0, "SourceTimeFunction: "
                  << name << " needs rise_time_s > 0; got " << rise << ".");
   }

   Kind GetKind() const { return kind_; }
   real_t RiseTime() const { return rise_; }

   /// f(t), the rate's time factor.
   real_t Rate(real_t t) const
   {
      if (t < 0.0) { return 0.0; }
      if (kind_ == Kind::kStep || t >= rise_) { return 1.0; }
      const real_t s = t / rise_;
      if (kind_ == Kind::kRamp) { return s; }
      return s * s * s * (10.0 + s * (-15.0 + 6.0 * s));
   }

   /// F(t) = int_0^t f dt, closed form.
   real_t Integral(real_t t) const
   {
      if (t <= 0.0) { return 0.0; }
      if (kind_ == Kind::kStep) { return t; }
      if (t >= rise_) { return t - 0.5 * rise_; }
      const real_t s = t / rise_;
      if (kind_ == Kind::kRamp) { return 0.5 * rise_ * s * s; }
      return rise_ * s * s * s * s * (2.5 + s * (-3.0 + s));
   }

private:
   Kind kind_ = Kind::kStep;
   real_t rise_ = 0.0;
};

/// a exp(-width |x - x0|^2) at the physical point; a from a callback, so the
/// injected volume follows its source.
class GaussianCoefficient : public Coefficient
{
public:
   GaussianCoefficient(const real_t origin[3], real_t width, std::function<real_t()> amp)
      : width_(width), amp_(std::move(amp))
   {
      for (int d = 0; d < 3; d++) { x0_[d] = origin[d]; }
   }

   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
      T.Transform(ip, x_);
      real_t r2 = 0.0;
      for (int d = 0; d < x_.Size(); d++) { r2 += (x_(d) - x0_[d]) * (x_(d) - x0_[d]); }
      return amp_() * std::exp(-width_ * r2);
   }

private:
   real_t x0_[3];
   real_t width_;
   std::function<real_t()> amp_;
   Vector x_;
};

/// p_init = magnitude exp(-width |x - x0|^2) (derivation section 4.2).
inline std::unique_ptr<Coefficient> InitialPressure(const PoroInitialPressure &ip)
{
   MFEM_VERIFY(ip.on, "InitialPressure: [multiphysics.poroelastic.initial_pressure] "
               "is not set.");
   const real_t a = ip.magnitude_pa;
   return std::make_unique<GaussianCoefficient>(ip.origin_m, ip.width_per_m2,
                                                [a]() { return a; });
}

enum class SourceIntegration { kExact, kTrapezoid };

#ifdef MFEM_USE_MPI

class PoroSource
{
public:
   /// The barycentric distance below which a well counts as on a face.
   static constexpr real_t kWellFaceTol = 1e-10;

   /// cfg.on must be set.  A well needs `fes`, the operator's vector DG space
   /// (byNODES, PoroNitscheOperator::Elastic().Elasticity().GetFESpace()), and
   /// the Biot sample for alpha M and M of its element; a gaussian ignores them.
   PoroSource(const PoroSourceConfig &cfg,
              SourceIntegration integ = SourceIntegration::kExact,
              ParFiniteElementSpace *fes = nullptr,
              const PoroMaterialSample *sample = nullptr)
      : cfg_(cfg), time_(cfg.time_function, cfg.rise_time_s), integ_(integ)
   {
      MFEM_VERIFY(cfg.on, "PoroSource: [multiphysics.poroelastic.source] is not set.");
      MFEM_VERIFY(cfg.shape == "gaussian" || cfg.shape == "well",
                  "PoroSource: unknown shape '" << cfg.shape << "'.");
      if (IsWell())
      {
         MFEM_VERIFY(fes && sample, "PoroSource: a well needs the vector DG space "
                     "and the Biot sample.");
         LocateWell_(*fes, *sample);
      }
      else
      {
         volume_ = std::make_unique<GaussianCoefficient>(
                      cfg.origin_m, cfg.width_per_m2,
                      [this]() { return cfg_.rate * factor_; });
      }
   }

   // The volume coefficient reads this object: no copies.
   PoroSource(const PoroSource &) = delete;
   PoroSource &operator=(const PoroSource &) = delete;

   bool IsWell() const { return cfg_.shape == "well"; }
   const SourceTimeFunction &TimeFunction() const { return time_; }
   real_t Rate() const { return cfg_.rate; }

   /// The time the source was last advanced to (0 at construction), and
   /// V / (rate g) there.
   real_t Time() const { return t_; }
   real_t VolumeFactor() const { return factor_; }

   /// Move the source to t >= Time(): the closed form, or one trapezoid step.
   void Advance(real_t t)
   {
      MFEM_VERIFY(t >= t_, "PoroSource::Advance: t = " << t << " is before the "
                  "source's time " << t_ << ".");
      if (integ_ == SourceIntegration::kExact) { factor_ = time_.Integral(t); }
      else { factor_ += 0.5 * (t - t_) * (time_.Rate(t_) + time_.Rate(t)); }
      t_ = t;
   }

   /// Gaussian only: V(x, Time()) = rate g(x) VolumeFactor().
   Coefficient &InjectedVolume()
   {
      MFEM_VERIFY(!IsWell(), "PoroSource::InjectedVolume: a well has no volume "
                  "coefficient; use AddWellLoad.");
      return *volume_;
   }

   /// Well only: the local element holding x_w on this rank, or -1.
   int WellElement() const { return well_elem_; }
   /// Well only: the load for V_w = 1, [alpha M div v; M div q](x_w), size 2n.
   const Vector &WellLoad() const { return well_load_; }

   /// Well only: F += rate VolumeFactor() WellLoad(), on the well's rows.
   void AddWellLoad(Vector &F) const
   {
      MFEM_VERIFY(IsWell(), "PoroSource::AddWellLoad: the source is a gaussian.");
      MFEM_VERIFY(F.Size() == well_load_.Size(), "PoroSource::AddWellLoad: F has "
                  "size " << F.Size() << " != 2 n = " << well_load_.Size());
      if (well_elem_ < 0) { return; }
      const real_t s = cfg_.rate * factor_;
      for (int i : well_rows_) { F(i) += s * well_load_(i); }
   }

private:
   /// Find the element holding x_w (all ranks agree, or all abort) and build
   /// the unit load there.
   void LocateWell_(ParFiniteElementSpace &fes, const PoroMaterialSample &sample)
   {
      ParMesh &mesh = *fes.GetParMesh();
      MFEM_VERIFY(fes.GetOrdering() == Ordering::byNODES && fes.GetVDim() == 3,
                  "PoroSource: the well load needs the vector DG space (byNODES, "
                  "vdim 3).");
      const real_t *xw = cfg_.origin_m;
      int inside = -1, n_inside = 0, n_touch = 0;
      real_t lam[4];
      IntegrationPoint ip_in;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         MFEM_VERIFY(mesh.GetElementGeometry(e) == Geometry::TETRAHEDRON,
                     "PoroSource: a well needs a tetrahedral mesh.");
         Barycentric_(mesh, e, xw, lam);
         const real_t lmin = std::min(std::min(lam[0], lam[1]), std::min(lam[2], lam[3]));
         if (lmin > kWellFaceTol)
         {
            inside = e;
            n_inside++;
            ip_in.Set3(lam[1], lam[2], lam[3]);
         }
         else if (lmin >= -kWellFaceTol) { n_touch++; }
      }
      int counts[2] = {n_inside, n_touch}, total[2];
      MPI_Allreduce(counts, total, 2, MPI_INT, MPI_SUM, mesh.GetComm());
      MFEM_VERIFY(total[1] == 0, "PoroSource: the well at (" << xw[0] << ", " << xw[1]
                  << ", " << xw[2] << ") lies on a face, edge or vertex of the mesh "
                  "(a barycentric coordinate within " << kWellFaceTol << " of zero); "
                  "a point source must be inside one element (derivation (4.3)), so "
                  "move origin_m.");
      MFEM_VERIFY(total[0] == 1, "PoroSource: the well at (" << xw[0] << ", " << xw[1]
                  << ", " << xw[2] << ") is outside the mesh (" << total[0]
                  << " elements hold it).");

      const int n = fes.GetVSize();
      well_load_.SetSize(2 * n);
      well_load_ = 0.0;
      if (inside < 0) { return; }
      well_elem_ = inside;
      const FiniteElement &fe = *fes.GetFE(inside);
      ElementTransformation &T = *mesh.GetElementTransformation(inside);
      T.SetIntPoint(&ip_in);
      Vector x;
      T.Transform(ip_in, x);
      real_t miss = 0.0, size = 0.0;
      for (int d = 0; d < 3; d++)
      {
         miss = std::max(miss, std::abs(x(d) - xw[d]));
         size = std::max(size, std::abs(xw[d]));
      }
      MFEM_VERIFY(miss <= 1e-10 * std::max(size, real_t(1.0)),
                  "PoroSource: the well's reference point maps to a point "
                  << miss << " m away from it.");
      const int nd = fe.GetDof();
      DenseMatrix dshape(nd, 3);
      fe.CalcPhysDShape(T, dshape);
      Array<int> vdofs;
      fes.GetElementVDofs(inside, vdofs);
      const BiotConstants &c = sample.Constants(inside);
      const real_t aM = c.alpha * c.M;
      for (int d = 0; d < 3; d++)
      {
         for (int k = 0; k < nd; k++)
         {
            const int i = vdofs[d * nd + k];
            well_load_(i) = aM * dshape(k, d);
            well_load_(n + i) = c.M * dshape(k, d);
            well_rows_.push_back(i);
            well_rows_.push_back(n + i);
         }
      }
   }

   /// Barycentric coordinates of x in tetrahedron e from its four vertices.
   static void Barycentric_(const Mesh &mesh, int e, const real_t *x, real_t lam[4])
   {
      Array<int> v;
      mesh.GetElementVertices(e, v);
      const real_t *p0 = mesh.GetVertex(v[0]);
      DenseMatrix J(3, 3);
      for (int j = 0; j < 3; j++)
      {
         const real_t *pj = mesh.GetVertex(v[j + 1]);
         for (int d = 0; d < 3; d++) { J(d, j) = pj[d] - p0[d]; }
      }
      DenseMatrixInverse Jinv(J);
      Vector r(3), xi(3);
      for (int d = 0; d < 3; d++) { r(d) = x[d] - p0[d]; }
      Jinv.Mult(r, xi);
      lam[1] = xi(0);
      lam[2] = xi(1);
      lam[3] = xi(2);
      lam[0] = 1.0 - xi(0) - xi(1) - xi(2);
   }

   PoroSourceConfig cfg_;
   SourceTimeFunction time_;
   SourceIntegration integ_;
   real_t t_ = 0.0;
   real_t factor_ = 0.0;
   std::unique_ptr<GaussianCoefficient> volume_;
   int well_elem_ = -1;
   Vector well_load_;
   std::vector<int> well_rows_;
};

#endif // MFEM_USE_MPI

} // namespace poro
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PORO_SOURCE_HPP
