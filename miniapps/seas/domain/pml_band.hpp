// domain/pml_band.hpp - the PML band geometry and its zero-frequency tensor.
//
// PML design plan rev2 (Chunhui Zhao, 3 Nov 2026), Sections 1, 2.2, 2.4,
// 2.6, 2.7, 3.5 and 4.1; implementation step P1 (Section 7.2), decisions
// D1-D6 at their Section 7.7 defaults (Chunhui Zhao, 4 Nov 2026).
//
// What this class gives, at any point x of the mesh:
//   A(x)        the walls whose bands contain x (0..D of them)
//   xi_a(x)     depth into band a, (1.1): xi_a = L_a - dist(x, wall a)
//   d_a(xi)     damping profile, (2.13): d_a = d0_a (xi_a / L_a)^N,
//               d0_a = (N + 1) c ln(1 / R0) / (2 L_a)
//   s0_a        zero-frequency stretch, (2.15): s0_a = 1 + d_a / alpha
//   Lambda0^-1  (2.18) with s0:  I + sum_a (1/s0_a - 1) n_a n_a^T
//   S0          det Lambda0 = prod_a s0_a
//   kappa       penalty factor, (3.14): max_j prod_{a in A, a != j} s0_a
// and C0 : G, the band stress of a displacement gradient G (G_kl = d_l u_k),
// (4.2):  sigma = S0 * sigma~(G Lambda0^-1) Lambda0^-1, with
// sigma~(H) = lambda tr(H) I + mu (H + H^T).
//
// The same code serves 2-D (plane strain) and 3-D (decision D5): in 2-D the
// vectors have D = 2 components and no wall faces e_3 (Section 2.7).
//
// Walls are found from the boundary faces of the listed attributes: faces are
// grouped by their outward normal (one attribute may hold several walls, the
// UCERF attribute 5 holds two opposite ones) and each group must be flat.
// Walls that meet must be perpendicular (|n_a . n_b| <= 1e-6) and the bands of
// opposite walls must not overlap.

#ifndef MFEM_SEAS_PML_BAND_HPP
#define MFEM_SEAS_PML_BAND_HPP

#include "mfem.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// The [cycle.pml] numbers the band needs (Section 7.2 deck keys).
struct PmlBandParams
{
   real_t thickness = 0.0;   ///< L (m), the same for every wall
   int power = 2;            ///< N of the profile (2.13)
   real_t R0 = 1e-3;         ///< target normal-incidence round-trip amplitude
   real_t f_alpha = 0.0;     ///< shift frequency f_alpha (Hz), alpha = 2 pi f_alpha
   real_t c_ref = 0.0;       ///< c of (2.13) (m/s), c_p by default
};

/// One flat wall: outward unit normal n and plane n . x = h.
struct PmlWall
{
   std::array<real_t, 3> n{{0.0, 0.0, 0.0}};
   real_t h = 0.0;
   int attr = 0;             ///< boundary attribute the faces came from
};

/// The band at one point (Section 2.6, at zero frequency).
struct PmlPoint
{
   int nact = 0;                     ///< |A(x)|
   int wall[3] = {-1, -1, -1};       ///< the active walls (indices into Walls())
   real_t xi[3] = {0.0, 0.0, 0.0};   ///< depth into each active band (m)
   real_t d[3] = {0.0, 0.0, 0.0};    ///< d_a (1/s)
   real_t s0[3] = {1.0, 1.0, 1.0};   ///< s0_a
   real_t Linv[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
   real_t S0 = 1.0;
   real_t kappa = 1.0;
};

class PmlBand
{
public:
   /// Walls given directly (unit tests, and the mesh constructor below).
   PmlBand(int dim, const std::vector<PmlWall> &walls, const PmlBandParams &prm)
      : dim_(dim), walls_(walls), prm_(prm)
   {
      Setup_();
   }

   /// Walls found from the boundary faces of `wall_attrs` on `mesh`.  With
   /// MPI, `comm` is the mesh's communicator: the faces of every rank are
   /// merged, so every rank holds the same walls in the same order.
   PmlBand(const Mesh &mesh, const std::vector<int> &wall_attrs,
           const PmlBandParams &prm
#ifdef MFEM_USE_MPI
           , MPI_Comm comm = MPI_COMM_NULL
#endif
          )
      : dim_(mesh.Dimension()), prm_(prm)
   {
      std::vector<real_t> loc;   // per local group: n[3], hmin, hmax, attr, count
      FindLocalWalls_(mesh, wall_attrs, loc);
      std::vector<real_t> all = loc;
#ifdef MFEM_USE_MPI
      if (comm != MPI_COMM_NULL)
      {
         int np = 1;
         MPI_Comm_size(comm, &np);
         int nloc = static_cast<int>(loc.size());
         std::vector<int> counts(np), displs(np, 0);
         MPI_Allgather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
         for (int r = 1; r < np; r++) { displs[r] = displs[r - 1] + counts[r - 1]; }
         all.assign(displs[np - 1] + counts[np - 1], 0.0);
         MPI_Allgatherv(loc.data(), nloc, MPITypeMap<real_t>::mpi_type, all.data(),
                        counts.data(), displs.data(), MPITypeMap<real_t>::mpi_type, comm);
      }
#endif
      MergeGroups_(all, mesh);
      Setup_();
   }

   int Dim() const { return dim_; }
   const std::vector<PmlWall> &Walls() const { return walls_; }
   const PmlBandParams &Params() const { return prm_; }
   real_t Alpha() const { return alpha_; }
   real_t D0() const { return d0_; }
   /// s0 on the wall itself, 1 + d0 / alpha.
   real_t S0Wall() const { return 1.0 + d0_ / alpha_; }
   /// The static extra length of one band, (4.4): D(L) / alpha.
   real_t DeltaL() const { return d0_ * prm_.thickness / ((prm_.power + 1) * alpha_); }

   /// Distance between the two opposite walls whose normals are along +-m
   /// (|n . m| > 1 - 1e-6, m a unit vector of Dim() components); 0 when no
   /// such pair is in the band set.
   real_t OppositeWallDistance(const real_t *m) const
   {
      const int nw = static_cast<int>(walls_.size());
      for (int a = 0; a < nw; a++)
         for (int b = a + 1; b < nw; b++)
         {
            real_t ab = 0.0, am = 0.0;
            for (int i = 0; i < dim_; i++)
            {
               ab += walls_[a].n[i] * walls_[b].n[i];
               am += walls_[a].n[i] * m[i];
            }
            if (ab < -1.0 + 1e-6 && std::abs(am) > 1.0 - 1e-6)
            { return walls_[a].h + walls_[b].h; }
         }
      return 0.0;
   }

   /// Decision D3 "auto": (W + 2 Delta L) / W, W the distance between the
   /// opposite walls along m (Section 4.4, (4.5)).
   real_t WallSpeedFactor(const real_t *m) const
   {
      const real_t W = OppositeWallDistance(m);
      MFEM_VERIFY(W > 0.0, "PmlBand::WallSpeedFactor: no pair of opposite walls with "
                  "normals along the loading plane normal; give wall_speed_factor as a "
                  "number");
      return (W + 2.0 * DeltaL()) / W;
   }

   /// The band at x (x has Dim() components).
   PmlPoint At(const real_t *x) const
   {
      PmlPoint p;
      const real_t L = prm_.thickness;
      for (int a = 0; a < static_cast<int>(walls_.size()); a++)
      {
         const PmlWall &w = walls_[a];
         real_t nx = 0.0;
         for (int i = 0; i < dim_; i++) { nx += w.n[i] * x[i]; }
         real_t xi = L - (w.h - nx);          // (1.1)
         if (xi <= 0.0) { continue; }
         if (xi > L) { xi = L; }              // round-off outside the wall
         MFEM_VERIFY(p.nact < dim_, "PmlBand: a point lies in more than " << dim_
                     << " bands");
         const int k = p.nact++;
         p.wall[k] = a;
         p.xi[k] = xi;
         p.d[k] = d0_ * std::pow(xi / L, prm_.power);
         p.s0[k] = 1.0 + p.d[k] / alpha_;
      }
      // Lambda0^-1 = I + sum_a (1/s0_a - 1) n_a n_a^T ; S0 = prod_a s0_a  (2.18)
      for (int k = 0; k < p.nact; k++)
      {
         const PmlWall &w = walls_[p.wall[k]];
         const real_t c = 1.0 / p.s0[k] - 1.0;
         for (int i = 0; i < dim_; i++)
            for (int j = 0; j < dim_; j++) { p.Linv[i][j] += c * w.n[i] * w.n[j]; }
         p.S0 *= p.s0[k];
      }
      // kappa (3.14): over an orthonormal frame holding the active normals,
      // prod_{a != j} s0_a is S0 / s0_j, and S0 for a frame direction that is
      // no active normal (exists when |A| < D).
      if (p.nact < dim_) { p.kappa = p.S0; }
      else
      {
         real_t smin = p.s0[0];
         for (int k = 1; k < p.nact; k++) { smin = std::min(smin, p.s0[k]); }
         p.kappa = p.S0 / smin;
      }
      if (p.nact == 0) { p.kappa = 1.0; }
      return p;
   }

   PmlPoint At(const Vector &x) const { return At(x.GetData()); }

   /// True when some band contains x.
   bool InBand(const real_t *x) const
   {
      for (const PmlWall &w : walls_)
      {
         real_t nx = 0.0;
         for (int i = 0; i < dim_; i++) { nx += w.n[i] * x[i]; }
         if (prm_.thickness - (w.h - nx) > 0.0) { return true; }
      }
      return false;
   }

   /// sigma = C0 : G at a band point (4.2), by the three stages of Section 3.3
   /// at zero frequency.  G and sigma are Dim() x Dim(), row-major,
   /// G[k][l] = d_l u_k.  lam, mu: the Lame parameters at the point.
   void ApplyC0(const PmlPoint &p, real_t lam, real_t mu,
                const real_t G[3][3], real_t sig[3][3]) const
   {
      const int D = dim_;
      real_t H[3][3] = {{0.0}};                 // stage 1: G Lambda0^-1
      for (int k = 0; k < D; k++)
         for (int r = 0; r < D; r++)
            for (int l = 0; l < D; l++) { H[k][r] += G[k][l] * p.Linv[l][r]; }
      real_t tr = 0.0;
      for (int k = 0; k < D; k++) { tr += H[k][k]; }
      real_t st[3][3];                          // stage 2: Hooke's law
      for (int i = 0; i < D; i++)
         for (int q = 0; q < D; q++)
         { st[i][q] = mu * (H[i][q] + H[q][i]) + (i == q ? lam * tr : 0.0); }
      for (int i = 0; i < D; i++)               // stage 3: S0 sigma~ Lambda0^-1
         for (int j = 0; j < D; j++)
         {
            real_t v = 0.0;
            for (int q = 0; q < D; q++) { v += st[i][q] * p.Linv[q][j]; }
            sig[i][j] = p.S0 * v;
         }
   }

private:
   int dim_ = 0;
   std::vector<PmlWall> walls_;
   PmlBandParams prm_;
   real_t alpha_ = 0.0;
   real_t d0_ = 0.0;

   void Setup_()
   {
      MFEM_VERIFY(dim_ == 2 || dim_ == 3, "PmlBand: dimension " << dim_);
      MFEM_VERIFY(prm_.thickness > 0.0, "PmlBand: thickness must be > 0 m");
      MFEM_VERIFY(prm_.power >= 1, "PmlBand: power must be >= 1");
      MFEM_VERIFY(prm_.R0 > 0.0 && prm_.R0 < 1.0, "PmlBand: R0 must be in (0, 1)");
      MFEM_VERIFY(prm_.f_alpha > 0.0, "PmlBand: f_alpha must be > 0 Hz (CFS-PML, D1)");
      MFEM_VERIFY(prm_.c_ref > 0.0, "PmlBand: c_ref must be > 0 m/s");
      MFEM_VERIFY(!walls_.empty(), "PmlBand: no walls");
      alpha_ = 2.0 * M_PI * prm_.f_alpha;
      d0_ = (prm_.power + 1) * prm_.c_ref * std::log(1.0 / prm_.R0)
            / (2.0 * prm_.thickness);                                 // (2.13)
      const int nw = static_cast<int>(walls_.size());
      for (int a = 0; a < nw; a++)
      {
         real_t nn = 0.0;
         for (int i = 0; i < dim_; i++) { nn += walls_[a].n[i] * walls_[a].n[i]; }
         MFEM_VERIFY(std::abs(nn - 1.0) < 1e-10, "PmlBand: wall " << a
                     << " normal is not a unit vector");
         for (int i = dim_; i < 3; i++) { walls_[a].n[i] = 0.0; }
      }
      for (int a = 0; a < nw; a++)
         for (int b = a + 1; b < nw; b++)
         {
            real_t c = 0.0;
            for (int i = 0; i < dim_; i++) { c += walls_[a].n[i] * walls_[b].n[i]; }
            if (c < -1.0 + 1e-6)
            {
               // opposite walls: n_b = -n_a, planes n_a . x = h_a and = -h_b
               const real_t gap = walls_[a].h + walls_[b].h;
               MFEM_VERIFY(2.0 * prm_.thickness < gap, "PmlBand: the bands of "
                           "opposite walls " << a << " and " << b << " overlap ("
                           << "2 L = " << 2.0 * prm_.thickness << " m, walls "
                           << gap << " m apart)");
               continue;
            }
            MFEM_VERIFY(c < 1.0 - 1e-6, "PmlBand: walls " << a << " and " << b
                        << " have the same normal");
            MFEM_VERIFY(std::abs(c) <= 1e-6, "PmlBand: walls " << a << " and " << b
                        << " are not perpendicular (n_a . n_b = " << c << ")");
         }
   }

   static void Centroid_(const Mesh &mesh, const Array<int> &v, int dim, real_t *c)
   {
      for (int i = 0; i < dim; i++) { c[i] = 0.0; }
      for (int k = 0; k < v.Size(); k++)
      {
         const real_t *x = mesh.GetVertex(v[k]);
         for (int i = 0; i < dim; i++) { c[i] += x[i] / v.Size(); }
      }
   }

   /// Groups this rank's wall faces by outward normal.  Output per group:
   /// n0, n1, n2, hmin, hmax, attr, count (7 numbers).
   void FindLocalWalls_(const Mesh &mesh, const std::vector<int> &attrs,
                        std::vector<real_t> &out) const
   {
      const int D = dim_;
      out.clear();
      Array<int> v, ev;
      Vector nor(D);
      for (int be = 0; be < mesh.GetNBE(); be++)
      {
         const int attr = mesh.GetBdrAttribute(be);
         if (std::find(attrs.begin(), attrs.end(), attr) == attrs.end()) { continue; }
         int el = -1, info = 0;
         mesh.GetBdrElementAdjacentElement(be, el, info);
         ElementTransformation *T =
            const_cast<Mesh &>(mesh).GetBdrElementTransformation(be);
         const IntegrationPoint &c =
            Geometries.GetCenter(mesh.GetBdrElementGeometry(be));
         T->SetIntPoint(&c);
         CalcOrtho(T->Jacobian(), nor);
         nor /= nor.Norml2();
         mesh.GetBdrElementVertices(be, v);
         real_t fc[3], ec[3];
         Centroid_(mesh, v, D, fc);
         mesh.GetElementVertices(el, ev);
         Centroid_(mesh, ev, D, ec);
         real_t out_dot = 0.0;
         for (int i = 0; i < D; i++) { out_dot += nor(i) * (fc[i] - ec[i]); }
         if (out_dot < 0.0) { nor.Neg(); }
         // offsets of the face vertices along the normal
         real_t hmin = 1e300, hmax = -1e300;
         for (int k = 0; k < v.Size(); k++)
         {
            const real_t *x = mesh.GetVertex(v[k]);
            real_t h = 0.0;
            for (int i = 0; i < D; i++) { h += nor(i) * x[i]; }
            hmin = std::min(hmin, h);
            hmax = std::max(hmax, h);
         }
         int g = -1;
         for (int k = 0; k < static_cast<int>(out.size()) / 7; k++)
         {
            real_t dot = 0.0;
            for (int i = 0; i < D; i++) { dot += out[7 * k + i] * nor(i); }
            if (dot > 1.0 - 1e-6 && static_cast<int>(out[7 * k + 5]) == attr) { g = k; break; }
         }
         if (g < 0)
         {
            g = static_cast<int>(out.size()) / 7;
            out.insert(out.end(), {0.0, 0.0, 0.0, hmin, hmax, real_t(attr), 0.0});
            for (int i = 0; i < D; i++) { out[7 * g + i] = nor(i); }
         }
         out[7 * g + 3] = std::min(out[7 * g + 3], hmin);
         out[7 * g + 4] = std::max(out[7 * g + 4], hmax);
         out[7 * g + 6] += 1.0;
      }
   }

   /// Merges the groups of every rank into walls (same order on every rank:
   /// by attribute, then by normal, compared component by component).
   void MergeGroups_(const std::vector<real_t> &all, const Mesh &mesh)
   {
      const int D = dim_;
      struct G { std::array<real_t, 3> n; real_t hmin, hmax; int attr; real_t cnt; };
      std::vector<G> gs;
      for (size_t k = 0; k + 7 <= all.size(); k += 7)
      {
         G g{{{all[k], all[k + 1], all[k + 2]}}, all[k + 3], all[k + 4],
             static_cast<int>(all[k + 5]), all[k + 6]};
         bool merged = false;
         for (G &o : gs)
         {
            real_t dot = 0.0;
            for (int i = 0; i < D; i++) { dot += o.n[i] * g.n[i]; }
            if (o.attr == g.attr && dot > 1.0 - 1e-6)
            {
               o.hmin = std::min(o.hmin, g.hmin);
               o.hmax = std::max(o.hmax, g.hmax);
               o.cnt += g.cnt;
               merged = true;
               break;
            }
         }
         if (!merged) { gs.push_back(g); }
      }
      std::sort(gs.begin(), gs.end(), [D](const G &a, const G &b)
      {
         if (a.attr != b.attr) { return a.attr < b.attr; }
         for (int i = 0; i < D; i++)
         {
            if (std::abs(a.n[i] - b.n[i]) > 1e-6) { return a.n[i] < b.n[i]; }
         }
         return false;
      });
      // size of the mesh, for the flatness tolerance
      real_t ext = 0.0;
      {
         Vector lo, hi;
         const_cast<Mesh &>(mesh).GetBoundingBox(lo, hi, 1);
         for (int i = 0; i < D; i++) { ext = std::max(ext, hi(i) - lo(i)); }
      }
      walls_.clear();
      for (const G &g : gs)
      {
         MFEM_VERIFY(g.hmax - g.hmin <= 1e-6 * std::max(ext, real_t(1.0)),
                     "PmlBand: the faces of attribute " << g.attr << " with normal ("
                     << g.n[0] << ", " << g.n[1] << ", " << g.n[2] << ") are not one "
                     "flat wall (offsets " << g.hmin << " .. " << g.hmax << " m)");
         PmlWall w;
         real_t nn = 0.0;
         for (int i = 0; i < D; i++) { nn += g.n[i] * g.n[i]; }
         for (int i = 0; i < D; i++) { w.n[i] = g.n[i] / std::sqrt(nn); }
         w.h = 0.5 * (g.hmin + g.hmax);
         w.attr = g.attr;
         walls_.push_back(w);
      }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PML_BAND_HPP
