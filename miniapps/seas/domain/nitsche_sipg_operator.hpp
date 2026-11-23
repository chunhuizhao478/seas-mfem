// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_sipg_operator.hpp — GPU dynamic-rupture plan D2 (gpu-dyn-nitsche
// plan rev2 §4.3; profiling item 8b, lts-kernels-design-rev1 §4.3-4.4):
// y = K u + K_eta v without a matrix, from per-element geometry and
// reference tables, in mfem::forall kernels that run on the CPU and the GPU.
//
// K is the IP stiffness of ElasticityDomainOperator (DGMethod::IP): MFEM's
// ElasticityIntegrator(lambda, mu) on the elements, and
// DGElasticityIPCombinedIntegrator(lambda, mu, 3, eps, penalty_factor) on
// every interior face (fault faces included), every shared face (the rows of
// this rank's element, as KeepNbrBlock(true) gives) and every Dirichlet
// face.  K_eta is the same operator with the Kelvin-Voigt moduli beta_e
// lambda_e, beta_e mu_e.  Absorbing and free faces add nothing to K.
//
// In the element e's own terms, with n_e the unnormalized normal of a face
// pointing out of e (CalcOrtho, |n_e| the face Jacobian), u_h the neighbour
// across the face, z = u + beta v on each side and sigma(G) = lambda tr(G) I
// + mu (G + G^T), K u + K_eta v at row (node k, component c) of e is
//
//   volume   sum_q w_q |J_e| (sigma_e(grad z_e) : grad phi_k)_c
//   faces    sum_q w_q [ phi_k ( -(a_e sigma_e(grad z_e) + a_h sigma_h(grad z_h)) n_e
//                                + pen ju + pen_eta jv )_c
//                        + c1 (sigma_e(sym(J (x) n_e)) grad phi_k)_c ]
//
// with ju = u_e - u_h, jv = v_e - v_h, J = ju + beta_e jv (beta_e is the
// test side's: the symmetry term of K_eta uses the test element's moduli),
// and per face kind
//   interior and shared  a_e = a_h = 1/2, c1 = eps/2,
//                        pen = pf (p_e + p_h)/4 |n| (ComputePenalty)
//   Dirichlet            a_e = 1, a_h = 0, u_h = v_h = 0, c1 = eps,
//                        pen = pf p_e |n|
// with p_i = 4 (p (p + 2)/3) (3 |n| / |J_i|) (3 lambda_i + 2 mu_i)^2 / (2 mu_i)
// (the integrator's (dim + 1) c_N (dim nl / detJ) ratio), and pen_eta the
// same with the Kelvin-Voigt moduli (linear in beta).  The integrators'
// coefficients c0 = -1/2, +1/2 (consistency), c1 (symmetry) and c2 = +-pen
// (penalty) per block, read row by row, give exactly these terms
// (dg_elasticity_ip_combined_integrator.hpp:122-180).
//
// On affine tets the integrands have degree at most 2p on a face and 2(p-1)
// in a tet, so the triangle rule of order 2p and MFEM's tet rule of order
// 2(p-1) (ElasticityIntegrator's own) reproduce the assembled matrices to
// round-off; only the order of the sums differs.  The test
// tests/unit/test_nitsche_sipg_mf.cpp checks ||y - y_asm||_inf <= 1e-13
// ||y_asm||_inf (plan rev2 §4.3 acceptance).
//
// Each element computes only its own rows, reading the neighbour's values
// across each face, so no row is written by two threads and the result does
// not depend on scheduling (plan rev2 Appendix C); the face work is done
// once from each side.  The kernels:
//   A  every element: the volume term and its interior and Dirichlet sides
//      (y_e = ...), while the halo of the shared faces travels;
//   B  the elements with shared sides: those sides (y_e += ...), after the
//      halo arrived.
// Each kernel has two layouts (SetLayout) with the same terms, which agree
// to round-off (the team forms each side's sigma n separately):
//   element  a thread per element, its rows in registers (FaceSide); the
//            host default, and the device default at p = 1.  On a GH200 at
//            p = 2 it needs more than 255 registers and spills (nvcc
//            -Xptxas -v).
//   team     a team of threads per element, 2 to 8 elements per block
//            (TeamElement): the values staged in shared memory, a thread per
//            quadrature point and side, then a thread per row; 52 to 78
//            registers, no spills.  The device default at p = 2, with 2
//            elements per block.
// The device defaults are the faster layouts of the D2 microbenchmark on a
// GH200 (job 1084877, TPV205 1000 m): at p = 1 the element layout, at p = 2
// the team with 2 elements per block.
// The halo is u (and v) on the DOFs of the face-neighbour elements, packed
// together in one message per neighbour rank (MFEM's face-neighbour tables,
// as ParGridFunction::ExchangeFaceNbrData), device buffers when MPI is GPU
// aware.
//
// Data per element: J^-1 (9), |J|, lambda, mu, beta (16 slots); per element
// side: the neighbour, the two sides' reference tables, n_e, pen, pen_eta,
// a_e.  Reference tables (shared by all elements): shape values and gradients
// at the volume points and, for each local face and orientation of a tet
// (64 lf + ori, Mesh::GetLocalFaceTransformation), at the face points.
//
// Scope: tetrahedra (the Nitsche path), uniform order p = 1 or 2, an order-0
// L2 material (NitscheMaterialSample), Kelvin-Voigt moduli proportional to
// the elastic ones per element (lambda_eta / lambda = mu_eta / mu = beta_e,
// as SetKelvinVoigt builds them); the constructor refuses anything else.

#ifndef MFEM_SEAS_NITSCHE_SIPG_OPERATOR_HPP
#define MFEM_SEAS_NITSCHE_SIPG_OPERATOR_HPP

#include "mfem.hpp"
#include "../common/stage_timers.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

namespace sipg_mf
{

/// sigma(G) = lambda tr(G) I + mu (G + G^T), G[c][d] = d z_c / d x_d.
MFEM_HOST_DEVICE inline void Stress(const real_t G[3][3], real_t lam, real_t mu,
                                    real_t S[3][3])
{
   const real_t tr = lam * (G[0][0] + G[1][1] + G[2][2]);
   for (int c = 0; c < 3; c++)
   {
      for (int d = 0; d < 3; d++) { S[c][d] = mu * (G[c][d] + G[d][c]); }
      S[c][c] += tr;
   }
}

/// G = Gr J^-1: the physical gradient from the reference one,
/// Jinv[3 r + d] = (J^-1)_{rd} (MFEM's Mult(dshape, InverseJacobian)).
MFEM_HOST_DEVICE inline void ToPhysical(const real_t Gr[3][3], const real_t *Jinv,
                                        real_t G[3][3])
{
   for (int c = 0; c < 3; c++)
      for (int d = 0; d < 3; d++)
      {
         G[c][d] = Gr[c][0] * Jinv[d] + Gr[c][1] * Jinv[3 + d] + Gr[c][2] * Jinv[6 + d];
      }
}

/// Sr = W S J^-T: the stress contracted against reference gradients,
/// (S grad phi_k)_c = sum_r Sr[c][r] dphi_k/dxi_r / W.
MFEM_HOST_DEVICE inline void ToReference(const real_t S[3][3], const real_t *Jinv,
                                         real_t W, real_t Sr[3][3])
{
   for (int c = 0; c < 3; c++)
      for (int r = 0; r < 3; r++)
      {
         Sr[c][r] = W * (S[c][0] * Jinv[3 * r] + S[c][1] * Jinv[3 * r + 1]
                         + S[c][2] * Jinv[3 * r + 2]);
      }
}

/// Element nodes of a tet of order P.
template <int P> struct Tet
{
   static constexpr int ND = (P + 1) * (P + 2) * (P + 3) / 6;
   static constexpr int NQV = (P == 1) ? 1 : 4;   ///< tet rule of order 2(P-1)
   static constexpr int NQF = (P == 1) ? 3 : 6;   ///< triangle rule of order 2P
};

/// Reference gradient Gr[c][r] = sum_k (u + beta v)[c][k] dphi_k/dxi_r of an
/// element whose values are read through get(c, k, m) (m = 0: u, 1: v).
template <int P, int NM, typename Get>
MFEM_HOST_DEVICE inline void RefGradient(const Get &get, real_t beta, const real_t *dphi,
                                         real_t Gr[3][3])
{
   constexpr int ND = Tet<P>::ND;
   for (int c = 0; c < 3; c++)
      for (int r = 0; r < 3; r++) { Gr[c][r] = 0.0; }
   for (int k = 0; k < ND; k++)
   {
      const real_t d0 = dphi[3 * k], d1 = dphi[3 * k + 1], d2 = dphi[3 * k + 2];
      for (int c = 0; c < 3; c++)
      {
         const real_t z = (NM == 2) ? get(c, k, 0) + beta * get(c, k, 1) : get(c, k, 0);
         Gr[c][0] += z * d0;
         Gr[c][1] += z * d1;
         Gr[c][2] += z * d2;
      }
   }
}

/// Traces t[m][c] = sum_k phi_k x_m[c][k].
template <int P, int NM, typename Get>
MFEM_HOST_DEVICE inline void Trace(const Get &get, const real_t *phi, real_t t[2][3])
{
   constexpr int ND = Tet<P>::ND;
   for (int m = 0; m < NM; m++)
      for (int c = 0; c < 3; c++)
      {
         real_t s = 0.0;
         for (int k = 0; k < ND; k++) { s += phi[k] * get(c, k, m); }
         t[m][c] = s;
      }
}

/// a_e sigma_e n_e + a_h sigma_h n_e at face point q: the stress traces of
/// both sides (sigma of grad(u + beta v)), weighted, against the normal.
template <int P, int NM, typename GetE, typename GetH>
MFEM_HOST_DEVICE inline void AvgTraction(const GetE &gete, const real_t *ge, const GetH &geth,
                                         const real_t *gh, bool nbr, real_t ae, const real_t *n,
                                         const real_t *dphe, const real_t *dphh, real_t t[3])
{
   real_t Gr[3][3], G[3][3], Se[3][3], Sh[3][3];
   RefGradient<P, NM>(gete, ge[12], dphe, Gr);
   ToPhysical(Gr, ge, G);
   Stress(G, ge[10], ge[11], Se);
   if (nbr)
   {
      RefGradient<P, NM>(geth, gh[12], dphh, Gr);
      ToPhysical(Gr, gh, G);
      Stress(G, gh[10], gh[11], Sh);
   }
   const real_t ah = 1.0 - ae;
   for (int c = 0; c < 3; c++)
   {
      t[c] = 0.0;
      for (int d = 0; d < 3; d++) { t[c] += (nbr ? ae * Se[c][d] + ah * Sh[c][d] : Se[c][d]) * n[d]; }
   }
}

/// One face point of an element side: q (3) and Sr (3 x 3) such that the
/// point adds phi_k q_c + sum_r Sr[c][r] dphi_k/dxi_r to row (k, c) of e.
/// t is the weighted stress trace (AvgTraction) at the point, phe / phh the
/// two sides' shape values there, W the face weight.  sd: n_e, pen, pen_eta,
/// a_e of the side; ge / gh as in FaceSide.
template <int P, int NM, typename GetE, typename GetH>
MFEM_HOST_DEVICE inline void FacePoint(const GetE &gete, const real_t *ge, const GetH &geth,
                                       bool nbr, const real_t *sd, const real_t *phe,
                                       const real_t *phh, real_t W, real_t eps, const real_t t[3],
                                       real_t qv[3], real_t Sr[3][3])
{
   const real_t lam_e = ge[10], mu_e = ge[11], beta_e = ge[12];
   const real_t pen = sd[3], pen_eta = sd[4], ae = sd[5];
   real_t te[2][3], th[2][3], ju[3], jv[3], J[3];
   Trace<P, NM>(gete, phe, te);
   if (nbr) { Trace<P, NM>(geth, phh, th); }
   for (int c = 0; c < 3; c++)
   {
      ju[c] = nbr ? te[0][c] - th[0][c] : te[0][c];
      jv[c] = (NM == 2) ? (nbr ? te[1][c] - th[1][c] : te[1][c]) : 0.0;
      qv[c] = W * (-t[c] + pen * ju[c] + ((NM == 2) ? pen_eta * jv[c] : 0.0));
      J[c] = (NM == 2) ? ju[c] + beta_e * jv[c] : ju[c];
   }
   // c1 sigma_e(sym(J (x) n_e)), then against the reference gradients
   const real_t Jn = lam_e * (J[0] * sd[0] + J[1] * sd[1] + J[2] * sd[2]);
   real_t S[3][3];
   for (int c = 0; c < 3; c++)
   {
      for (int d = 0; d < 3; d++) { S[c][d] = mu_e * (J[c] * sd[d] + sd[c] * J[d]); }
      S[c][c] += Jn;
   }
   ToReference(S, ge, eps * ae * W, Sr);
}

/// One element side's rows, added to ye.  e: the element (values through
/// gete, geometry ge); h: the neighbour (geth, gh), absent for a Dirichlet
/// side (nbr = false).  sd: n_e[3], pen, pen_eta, a_e.  Tables: phi / dphi of
/// the two sides at the face points (table index * NQF * ND), weights fw.
/// At P = 1 the gradients are constant on a tet, so the stress traces are
/// computed once per side instead of once per face point.
template <int P, int NM, typename GetE, typename GetH>
MFEM_HOST_DEVICE inline void FaceSide(const GetE &gete, const real_t *ge, const GetH &geth,
                                      const real_t *gh, bool nbr, const real_t *sd,
                                      int tab_e, int tab_h, real_t eps, const real_t *fw,
                                      const real_t *fphi, const real_t *fdphi,
                                      real_t ye[3][Tet<P>::ND])
{
   constexpr int ND = Tet<P>::ND, NQF = Tet<P>::NQF;
   real_t t[3];
   if (P == 1)
   {
      AvgTraction<P, NM>(gete, ge, geth, gh, nbr, sd[5], sd, fdphi + tab_e * NQF * ND * 3,
                         fdphi + tab_h * NQF * ND * 3, t);
   }
   for (int q = 0; q < NQF; q++)
   {
      const real_t *phe = fphi + (tab_e * NQF + q) * ND;
      const real_t *dphe = fdphi + (tab_e * NQF + q) * ND * 3;
      if (P != 1)
      {
         AvgTraction<P, NM>(gete, ge, geth, gh, nbr, sd[5], sd, dphe,
                            fdphi + (tab_h * NQF + q) * ND * 3, t);
      }
      real_t qv[3], Sr[3][3];
      FacePoint<P, NM>(gete, ge, geth, nbr, sd, phe, fphi + (tab_h * NQF + q) * ND, fw[q], eps,
                       t, qv, Sr);
      for (int k = 0; k < ND; k++)
      {
         const real_t p = phe[k], d0 = dphe[3 * k], d1 = dphe[3 * k + 1], d2 = dphe[3 * k + 2];
         for (int c = 0; c < 3; c++)
         {
            ye[c][k] += p * qv[c] + Sr[c][0] * d0 + Sr[c][1] * d1 + Sr[c][2] * d2;
         }
      }
   }
}

/// Data layout constants and side kinds (side_i_[3 s]: the neighbour, or a
/// kind below).
constexpr int kGeo = 16;     ///< per element: J^-1, |J|, lambda, mu, beta, lam_eta, mu_eta
constexpr int kSide = 6;     ///< per side: n_e, pen, pen_eta, a_e
constexpr int kNone = -1;    ///< side without a term (free, absorbing)
constexpr int kDirichlet = -2;
constexpr int kGhost0 = -3;  ///< kGhost0 - nbr: face-neighbour element nbr

/// Device pointers and sizes of one product, for the team kernels.
struct KernelData
{
   const real_t *du, *dv, *halo, *geo, *sd, *vw, *vdphi, *fw, *fphi, *fdphi;
   const int *si, *gv;
   int nsc, ne;
   real_t eps;
   bool do_vol, do_faces;
};

/// Team layout sizes: NQT stress-trace points per side (one at P = 1, where
/// the gradients are constant on a tet), NPT quadrature points (volume, then
/// 4 sides x NQF), NV values per field (3 ND rows), TX threads per element.
template <int P> struct Team
{
   static constexpr int ND = Tet<P>::ND, NQV = Tet<P>::NQV, NQF = Tet<P>::NQF;
   static constexpr int NQT = (P == 1) ? 1 : NQF;
   static constexpr int NV = 3 * ND, NPT = NQV + 4 * NQF;
   static constexpr int N1 = NQV + 8 * NQT;   ///< phase 1 items
   static constexpr int TX = (N1 > NV ? N1 : NV) > 4 * NQF ? (N1 > NV ? N1 : NV) : 4 * NQF;
   /// Shared memory of a block of BZ elements (Team_'s arrays), and whether
   /// it fits the 48 KB of static shared memory of a CUDA block.
   static constexpr long Bytes(int NM, int BZ)
   {
      return long(sizeof(real_t)) * BZ * (5 * NM * NV + NPT * 12 + 8 * NQT * 3);
   }
   static constexpr bool Fits(int NM, int BZ) { return Bytes(NM, BZ) <= 48 * 1024; }
};

/// One element in the team layout (a team of threads per element, the
/// shared arrays below per element):
///   0  stage u (and v) of e and of its side neighbours (ue: [NM][3][ND],
///      uh: [4][NM][3][ND]);
///   1  a thread per volume point (its Sr in sp) and per (side, trace point,
///      side of the face): a sigma(grad z) n_e of that side in tr, so each
///      thread holds one gradient;
///   2  a thread per face point: its q and Sr (FacePoint) in sp;
///   3  a thread per row (c, k): the sum over the points.
/// GHOST = false (kernel A): the volume and the interior and Dirichlet
/// sides, y = ...; GHOST = true (kernel B): the shared sides, y += ....  The
/// terms are FaceSide's and kernel A's, summed in another order.
template <int P, int NM, bool GHOST>
MFEM_HOST_DEVICE inline void TeamElement(int e, const KernelData &a, real_t *dy, real_t *ue,
                                         real_t *uh, real_t (*sp)[12], real_t (*tr)[3])
{
   constexpr int ND = Tet<P>::ND, NQV = Tet<P>::NQV, NQF = Tet<P>::NQF;
   constexpr int NQT = Team<P>::NQT, NV = Team<P>::NV;
   const real_t *ge = a.geo + kGeo * e;
   const int *si = a.si + 12 * e;
   const real_t *sde = a.sd + 4 * kSide * e;
   auto active = [&](int s) {
      const int nbr = si[3 * s];
      return a.do_faces && (GHOST ? nbr <= kGhost0 : (nbr != kNone && nbr > kGhost0));
   };

   MFEM_FOREACH_THREAD(i, x, 5 * NM * NV)
   {
      const int blk = i / (NM * NV), j = i % (NM * NV), m = j / NV, ck = j % NV;
      const int c = ck / ND, k = ck % ND;
      const real_t *src = (m == 0) ? a.du : a.dv;
      if (blk == 0) { ue[j] = src[c * a.nsc + e * ND + k]; }
      else
      {
         const int s = blk - 1, nbr = si[3 * s];
         real_t val = 0.0;
         if (active(s) && nbr >= 0) { val = src[c * a.nsc + nbr * ND + k]; }
         else if (active(s) && nbr <= kGhost0)
         {
            val = a.halo[NM * a.gv[NV * (kGhost0 - nbr) + ck] + m];
         }
         uh[s * NM * NV + j] = val;
      }
   }
   MFEM_SYNC_THREAD;

   MFEM_FOREACH_THREAD(it, x, NQV + 8 * NQT)
   {
      real_t Gr[3][3], G[3][3], S[3][3];
      if (it < NQV)
      {
         real_t *o = sp[it];
         for (int i = 0; i < 12; i++) { o[i] = 0.0; }
         if (!GHOST && a.do_vol)
         {
            const real_t *dph = a.vdphi + it * ND * 3;
            auto get = [&](int c, int k, int m) { return ue[m * NV + c * ND + k]; };
            real_t Sr[3][3];
            RefGradient<P, NM>(get, ge[12], dph, Gr);
            ToPhysical(Gr, ge, G);
            Stress(G, ge[10], ge[11], S);
            ToReference(S, ge, a.vw[it] * ge[9], Sr);
            for (int c = 0; c < 3; c++)
               for (int r = 0; r < 3; r++) { o[3 + 3 * c + r] = Sr[c][r]; }
         }
      }
      else
      {
         // a_i sigma_i(grad z_i) n_e of side i (0: e, 1: the neighbour) of
         // face s at trace point q; a_e = 1 on a Dirichlet side, whose
         // neighbour half stays zero
         const int j = it - NQV, s = j / (2 * NQT), q = (j / 2) % NQT, h = j % 2;
         real_t *o = tr[j];
         o[0] = o[1] = o[2] = 0.0;
         if (active(s) && (h == 0 || si[3 * s] != kDirichlet))
         {
            const int nbr = si[3 * s];
            const real_t *sd = sde + kSide * s;
            const real_t *uv = (h == 0) ? ue : uh + s * NM * NV;
            const real_t *g = (h == 0) ? ge
                              : a.geo + kGeo * (nbr >= 0 ? nbr : a.ne + kGhost0 - nbr);
            auto get = [&](int c, int k, int m) { return uv[m * NV + c * ND + k]; };
            const int tab = si[3 * s + 1 + h];
            RefGradient<P, NM>(get, g[12], a.fdphi + (tab * NQF + q) * ND * 3, Gr);
            ToPhysical(Gr, g, G);
            Stress(G, g[10], g[11], S);
            const real_t w = (h == 0) ? sd[5] : 1.0 - sd[5];
            for (int c = 0; c < 3; c++)
            {
               o[c] = w * (S[c][0] * sd[0] + S[c][1] * sd[1] + S[c][2] * sd[2]);
            }
         }
      }
   }
   MFEM_SYNC_THREAD;

   MFEM_FOREACH_THREAD(j, x, 4 * NQF)
   {
      const int s = j / NQF, q = j % NQF;
      real_t *o = sp[NQV + j];
      for (int i = 0; i < 12; i++) { o[i] = 0.0; }
      if (active(s))
      {
         const real_t *uhs = uh + s * NM * NV;
         auto gete = [&](int c, int k, int m) { return ue[m * NV + c * ND + k]; };
         auto geth = [&](int c, int k, int m) { return uhs[m * NV + c * ND + k]; };
         const real_t *sd = sde + kSide * s;
         const int tab_e = si[3 * s + 1], tab_h = si[3 * s + 2];
         const real_t *t0 = tr[2 * (s * NQT + (P == 1 ? 0 : q))];
         const real_t t[3] = {t0[0] + t0[3], t0[1] + t0[4], t0[2] + t0[5]};
         real_t Sr[3][3];
         FacePoint<P, NM>(gete, ge, geth, si[3 * s] != kDirichlet, sd,
                          a.fphi + (tab_e * NQF + q) * ND, a.fphi + (tab_h * NQF + q) * ND,
                          a.fw[q], a.eps, t, o, Sr);
         for (int c = 0; c < 3; c++)
            for (int r = 0; r < 3; r++) { o[3 + 3 * c + r] = Sr[c][r]; }
      }
   }
   MFEM_SYNC_THREAD;

   MFEM_FOREACH_THREAD(row, x, NV)
   {
      const int c = row / ND, k = row % ND;
      real_t acc = 0.0;
      if (!GHOST && a.do_vol)
      {
         for (int q = 0; q < NQV; q++)
         {
            const real_t *d = a.vdphi + (q * ND + k) * 3;
            const real_t *o = sp[q];
            acc += o[3 + 3 * c] * d[0] + o[4 + 3 * c] * d[1] + o[5 + 3 * c] * d[2];
         }
      }
      for (int s = 0; s < 4; s++)
      {
         if (!active(s)) { continue; }
         const int tab_e = si[3 * s + 1];
         for (int q = 0; q < NQF; q++)
         {
            const real_t *o = sp[NQV + s * NQF + q];
            const real_t ph = a.fphi[(tab_e * NQF + q) * ND + k];
            const real_t *d = a.fdphi + ((tab_e * NQF + q) * ND + k) * 3;
            acc += ph * o[c] + o[3 + 3 * c] * d[0] + o[4 + 3 * c] * d[1] + o[5 + 3 * c] * d[2];
         }
      }
      if (GHOST) { dy[c * a.nsc + e * ND + k] += acc; }
      else { dy[c * a.nsc + e * ND + k] = acc; }
   }
}

} // namespace sipg_mf

/// y = K u (+ K_eta v) of the Nitsche IP stiffness, matrix-free (header).
class NitscheSipgOperator
{
public:
   /// What Mult computes: both passes (the product) or one of them (the
   /// microbenchmark times them apart; one pass alone is not a product).
   enum Pass { kVolume = 1, kFaces = 2, kAll = 3 };
   /// Kernel layout: a thread per element (registers), or a team of threads
   /// per element (shared memory, TeamElement).  kAuto: on a device the team
   /// at p = 2 and a thread per element at p = 1 (header), a thread per
   /// element on the host.  Same terms, sums in another order (round-off).
   enum Layout { kAuto = 0, kElement = 1, kTeam = 2 };

   /// fes: the DG space of K (ElasticityDomainOperator::GetFESpace; vdim 3,
   /// byNODES).  lambda, mu: the order-0 L2 material (NitscheMaterialSample)
   /// with face-neighbour data exchanged; kv_lambda, kv_mu: the Kelvin-Voigt
   /// moduli on the same space (both or neither).  dirichlet_marker: per
   /// boundary attribute, 1 = Dirichlet (K's boundary-face integrator).
   /// Collective (face-neighbour exchanges).
   NitscheSipgOperator(ParFiniteElementSpace &fes, const ParGridFunction &lambda,
                       const ParGridFunction &mu, const ParGridFunction *kv_lambda,
                       const ParGridFunction *kv_mu, const Array<int> &dirichlet_marker,
                       real_t penalty_factor, real_t epsilon = -1.0)
      : fes_(fes), eps_(epsilon), kv_(kv_lambda != nullptr)
   {
      MFEM_VERIFY((kv_lambda == nullptr) == (kv_mu == nullptr),
                  "NitscheSipgOperator: give both Kelvin-Voigt moduli or neither");
      ParMesh &pm = *fes.GetParMesh();
      MFEM_VERIFY(pm.Dimension() == 3 && pm.SpaceDimension() == 3,
                  "NitscheSipgOperator: 3D meshes only");
      MFEM_VERIFY(fes.GetVDim() == 3 && fes.GetOrdering() == Ordering::byNODES,
                  "NitscheSipgOperator: the space must be vdim 3, byNODES");
      ne_ = pm.GetNE();
      order_ = fes.GetMaxElementOrder();
      MFEM_VERIFY(order_ == 1 || order_ == 2, "NitscheSipgOperator: order " << order_
                  << "; the kernels cover p = 1 and 2");
      nd_ = (order_ + 1) * (order_ + 2) * (order_ + 3) / 6;
      nsc_ = fes.GetNDofs();
      MFEM_VERIFY(nsc_ == ne_ * nd_, "NitscheSipgOperator: " << nsc_ << " scalar DOFs for "
                  << ne_ << " elements of " << nd_ << " nodes; not a DG space");
      pm.ExchangeFaceNbrData();
      fes.ExchangeFaceNbrData();

      // the element DOFs are e nd + k in each component block (the kernels
      // index them so)
      Array<int> vd;
      for (int e = 0; e < ne_; e++)
      {
         MFEM_VERIFY(pm.GetElementBaseGeometry(e) == Geometry::TETRAHEDRON,
                     "NitscheSipgOperator: element " << e << " is not a tetrahedron");
         MFEM_VERIFY(fes.GetFE(e)->GetOrder() == order_,
                     "NitscheSipgOperator: element " << e << " has order "
                     << fes.GetFE(e)->GetOrder() << " != " << order_);
         fes.GetElementVDofs(e, vd);
         MFEM_VERIFY(vd.Size() == 3 * nd_, "NitscheSipgOperator: element " << e << " has "
                     << vd.Size() << " vector DOFs");
         for (int c = 0; c < 3; c++)
            for (int k = 0; k < nd_; k++)
            {
               MFEM_VERIFY(vd[c * nd_ + k] == c * nsc_ + e * nd_ + k,
                           "NitscheSipgOperator: element " << e << " DOF (" << c << ", " << k
                           << ") is " << vd[c * nd_ + k] << ", not " << c * nsc_ + e * nd_ + k);
            }
      }

      BuildTables_(pm, *fes.GetFE(0));
      const int nfn = pm.face_nbr_elements.Size();
      ng_ = nfn;
      geo_.SetSize(kGeo * (ne_ + ng_));
      geo_ = 0.0;
      gvdof_.SetSize(3 * nd_ * ng_);
      gvdof_ = -1;
      side_i_.SetSize(4 * 3 * ne_);
      side_d_.SetSize(4 * kSide * ne_);
      for (int i = 0; i < side_i_.Size(); i++) { side_i_[i] = kNone; }
      side_d_ = 0.0;

      // material, per element and per face-neighbour element
      const ParFiniteElementSpace &mfes = *lambda.ParFESpace();
      MFEM_VERIFY(mfes.GetVSize() == ne_ && mu.ParFESpace() == &mfes
                  && (!kv_ || (kv_lambda->ParFESpace() == &mfes && kv_mu->ParFESpace() == &mfes)),
                  "NitscheSipgOperator: the material must be order-0 L2 (one value per element), "
                  "all fields on one space");
      // read on the host below; with a device the fields and their
      // face-neighbour data may hold their valid copy there
      for (const ParGridFunction *g : {&lambda, &mu, kv_lambda, kv_mu})
      {
         if (g) { g->HostRead(); g->FaceNbrData().HostRead(); }
      }
      auto material = [&](int el, bool ghost, real_t m[4]) {
         // lambda, mu, lambda_eta, mu_eta; a face-neighbour element reads the
         // GFs' face-neighbour data (as GridFunctionCoefficient::Eval does)
         Array<int> md;
         if (ghost) { mfes.GetFaceNbrElementVDofs(el, md); }
         else { mfes.GetElementDofs(el, md); }
         MFEM_VERIFY(md.Size() == 1, "NitscheSipgOperator: material space is not order-0 L2");
         auto val = [&](const ParGridFunction &g) {
            return ghost ? g.FaceNbrData()(md[0]) : g(md[0]);
         };
         m[0] = val(lambda);
         m[1] = val(mu);
         m[2] = kv_ ? val(*kv_lambda) : 0.0;
         m[3] = kv_ ? val(*kv_mu) : 0.0;
      };
      auto fill_geo = [&](int slot, ElementTransformation &T, const real_t m[4]) {
         // affine tets: J^-1 and |J| at the centre hold everywhere
         T.SetIntPoint(&Geometries.GetCenter(Geometry::TETRAHEDRON));
         const DenseMatrix &Ji = T.InverseJacobian();
         real_t *g = &geo_[kGeo * slot];
         for (int r = 0; r < 3; r++)
            for (int d = 0; d < 3; d++) { g[3 * r + d] = Ji(r, d); }
         g[9] = T.Weight();
         g[10] = m[0];
         g[11] = m[1];
         MFEM_VERIFY(m[1] > 0.0, "NitscheSipgOperator: mu = " << m[1] << " <= 0 at element "
                     << slot << " (the IP penalty divides by mu)");
         real_t beta = 0.0;
         if (kv_)
         {
            // lambda_eta = beta lambda and mu_eta = beta mu (SetKelvinVoigt);
            // the kernels use z = u + beta v, which needs one beta per element
            beta = m[3] / m[1];
            MFEM_VERIFY(std::abs(m[2] - beta * m[0]) <= 1e-13 * std::abs(m[2]) + 1e-300,
                        "NitscheSipgOperator: Kelvin-Voigt moduli not proportional to the "
                        "elastic ones at element " << slot << " (lambda_eta / lambda = "
                        << m[2] / m[0] << ", mu_eta / mu = " << beta << ")");
         }
         g[12] = beta;
         g[13] = m[2];
         g[14] = m[3];
      };
      for (int e = 0; e < ne_; e++)
      {
         real_t m[4];
         material(e, false, m);
         fill_geo(e, *pm.GetElementTransformation(e), m);
      }

      // penalty of one face side: p_i of ComputePenalty, from the side's
      // moduli and |J|, at |n| = nl
      const real_t pf = penalty_factor;
      const int dim = 3;
      const real_t c_N = order_ * (order_ + dim - 1.0) / dim;
      auto p_side = [&](real_t lam, real_t mu, real_t detJ, real_t nl) {
         const real_t c1v = dim * lam + 2.0 * mu;
         const real_t ratio = c1v * c1v / (2.0 * mu);
         return (dim + 1) * c_N * (real_t(dim) * nl / detJ) * ratio;
      };
      auto set_side = [&](int e, int lf, int nbr, int tab_e, int tab_h, const Vector &nor,
                          real_t sign, real_t pen, real_t pen_eta, real_t ae) {
         const int s = 4 * e + lf;
         MFEM_VERIFY(side_i_[3 * s] == kNone, "NitscheSipgOperator: element " << e
                     << " local face " << lf << " met twice");
         side_i_[3 * s] = nbr;
         side_i_[3 * s + 1] = tab_e;
         side_i_[3 * s + 2] = tab_h;
         real_t *d = &side_d_[kSide * s];
         for (int c = 0; c < 3; c++) { d[c] = sign * nor(c); }
         d[3] = pen;
         d[4] = pen_eta;
         d[5] = ae;
      };
      auto tab = [&](int inf) {
         const int lf = inf / 64, ori = inf % 64;
         MFEM_VERIFY(lf >= 0 && lf < 4 && ori >= 0 && ori < 6,
                     "NitscheSipgOperator: face info " << inf << " out of range");
         return 6 * lf + ori;
      };
      Vector nor(3);
      const IntegrationPoint &fc = Geometries.GetCenter(Geometry::TRIANGLE);

      // interior faces (fault faces included)
      for (int f = 0; f < pm.GetNumFaces(); f++)
      {
         int e1, e2, inf1, inf2;
         pm.GetFaceElements(f, &e1, &e2);
         if (e1 < 0 || e2 < 0 || e1 == e2) { continue; }
         pm.GetFaceInfos(f, &inf1, &inf2);
         FaceElementTransformations *T = pm.GetInteriorFaceTransformations(f);
         MFEM_VERIFY(T, "NitscheSipgOperator: no transformation for interior face " << f);
         T->SetAllIntPoints(&fc);
         CalcOrtho(T->Jacobian(), nor);
         const real_t nl = nor.Norml2();
         const real_t d1 = T->Elem1->Weight(), d2 = T->Elem2->Weight();
         const real_t *g1 = &geo_[kGeo * e1], *g2 = &geo_[kGeo * e2];
         const real_t pen = pf * (p_side(g1[10], g1[11], d1, nl) + p_side(g2[10], g2[11], d2, nl))
                            / 4.0 * nl;
         const real_t pen_eta = kv_ ? pf * (p_side(g1[13], g1[14], d1, nl)
                                            + p_side(g2[13], g2[14], d2, nl)) / 4.0 * nl : 0.0;
         set_side(e1, inf1 / 64, e2, tab(inf1), tab(inf2), nor, 1.0, pen, pen_eta, 0.5);
         set_side(e2, inf2 / 64, e1, tab(inf2), tab(inf1), nor, -1.0, pen, pen_eta, 0.5);
         n_int_sides_ += 2;
      }

      // shared faces: this rank's element is Elem1; the neighbour is a
      // face-neighbour element, a ghost slot ne + nbr after the local ones
      for (int sf = 0; sf < pm.GetNSharedFaces(); sf++)
      {
         const int f = pm.GetSharedFace(sf);
         int e1, e2, inf1, inf2;
         pm.GetFaceElements(f, &e1, &e2);
         pm.GetFaceInfos(f, &inf1, &inf2);
         MFEM_VERIFY(e1 >= 0 && e2 < 0, "NitscheSipgOperator: shared face " << sf
                     << " has elements " << e1 << ", " << e2);
         const int nb = -1 - e2;
         FaceElementTransformations *T = pm.GetSharedFaceTransformations(sf);
         T->SetAllIntPoints(&fc);
         CalcOrtho(T->Jacobian(), nor);
         const real_t nl = nor.Norml2();
         if (gvdof_[3 * nd_ * nb] < 0)
         {
            real_t m[4];
            material(nb, true, m);
            fill_geo(ne_ + nb, *T->Elem2, m);
            fes.GetFaceNbrElementVDofs(nb, vd);
            MFEM_VERIFY(vd.Size() == 3 * nd_, "NitscheSipgOperator: face-neighbour element "
                        << nb << " has " << vd.Size() << " vector DOFs");
            for (int i = 0; i < 3 * nd_; i++)
            {
               MFEM_VERIFY(vd[i] >= 0, "NitscheSipgOperator: negative face-neighbour DOF");
               gvdof_[3 * nd_ * nb + i] = vd[i];
            }
         }
         const real_t d1 = T->Elem1->Weight(), d2 = T->Elem2->Weight();
         const real_t *g1 = &geo_[kGeo * e1], *g2 = &geo_[kGeo * (ne_ + nb)];
         const real_t pen = pf * (p_side(g1[10], g1[11], d1, nl) + p_side(g2[10], g2[11], d2, nl))
                            / 4.0 * nl;
         const real_t pen_eta = kv_ ? pf * (p_side(g1[13], g1[14], d1, nl)
                                            + p_side(g2[13], g2[14], d2, nl)) / 4.0 * nl : 0.0;
         set_side(e1, inf1 / 64, kGhost0 - nb, tab(inf1), tab(inf2), nor, 1.0, pen, pen_eta, 0.5);
         n_shared_sides_++;
      }

      // Dirichlet faces (true boundary faces with a marked attribute)
      for (int be = 0; be < pm.GetNBE(); be++)
      {
         const int attr = pm.GetBdrAttribute(be);
         if (attr < 1 || attr > dirichlet_marker.Size() || dirichlet_marker[attr - 1] == 0) { continue; }
         const int f = pm.GetBdrElementFaceIndex(be);
         int e1, e2, inf1, inf2;
         pm.GetFaceElements(f, &e1, &e2);
         pm.GetFaceInfos(f, &inf1, &inf2);
         if (e2 >= 0 || inf2 >= 0) { continue; }   // interior or shared: not a boundary face
         FaceElementTransformations *T = pm.GetBdrFaceTransformations(be);
         MFEM_VERIFY(T && T->Elem1No == e1, "NitscheSipgOperator: boundary element " << be
                     << " has no boundary-face transformation");
         T->SetAllIntPoints(&fc);
         CalcOrtho(T->Jacobian(), nor);
         const real_t nl = nor.Norml2();
         const real_t d1 = T->Elem1->Weight();
         const real_t *g1 = &geo_[kGeo * e1];
         const real_t pen = pf * p_side(g1[10], g1[11], d1, nl) * nl;
         const real_t pen_eta = kv_ ? pf * p_side(g1[13], g1[14], d1, nl) * nl : 0.0;
         set_side(e1, inf1 / 64, kDirichlet, tab(inf1), 0, nor, 1.0, pen, pen_eta, 1.0);
         n_dir_sides_++;
      }

      // the elements with shared sides (kernel B)
      std::vector<int> be_list;
      for (int e = 0; e < ne_; e++)
      {
         for (int s = 0; s < 4; s++)
         {
            if (side_i_[3 * (4 * e + s)] <= kGhost0) { be_list.push_back(e); break; }
         }
      }
      shared_elems_.SetSize(static_cast<int>(be_list.size()));
      for (int i = 0; i < shared_elems_.Size(); i++) { shared_elems_[i] = be_list[i]; }

      SetupHalo_(pm);
      for (auto *a : {&geo_, &side_d_}) { a->UseDevice(true); }
      for (auto *a : {&side_i_, &gvdof_, &shared_elems_}) { a->GetMemory().UseDevice(true); }
   }

   bool HasKelvinVoigt() const { return kv_; }
   int Order() const { return order_; }
   int NumElements() const { return ne_; }
   /// Element sides of this rank by kind: interior (both directions),
   /// shared, Dirichlet.
   long long NumInteriorSides() const { return n_int_sides_; }
   long long NumSharedSides() const { return n_shared_sides_; }
   long long NumDirichletSides() const { return n_dir_sides_; }
   /// Bytes of the per-element data, the side data and the halo index on
   /// this rank (the reference tables are shared and small).
   long long Bytes() const
   {
      return static_cast<long long>(sizeof(real_t)) * (geo_.Size() + side_d_.Size())
             + static_cast<long long>(sizeof(int)) * (side_i_.Size() + gvdof_.Size()
                                                     + shared_elems_.Size() + send_map_.Size());
   }
   /// Model flop count of one Mult on this rank (both passes, or one): the
   /// operations of the loops below, a multiply-add counted as two flops.
   double Flops(bool with_kv, int pass = kAll) const
   {
      const double nd = nd_, dz = with_kv ? 1.0 : 0.0, nm = with_kv ? 2.0 : 1.0;
      const double nqv = (order_ == 1) ? 1 : 4, nqf = (order_ == 1) ? 3 : 6;
      // gradient (z = u + beta v and 3 multiply-adds per node and component),
      // J^-1 (45), stress (24)
      const double grad = 2.0 * 3.0 * nd * (3.0 + dz) + 45.0 + 24.0;
      const double vol = nqv * (grad + 55.0 + 18.0 * nd);   // J^-T (55), residual
      // per side: the stress traces of both sides and the flux (45), per
      // point at p = 2, once per side at p = 1; per point: the value traces,
      // jumps, flux and symmetry stress (75), J^-T (54), residual (24 nd)
      auto side = [&](bool nbr) {
         const double ntr = nbr ? 2.0 : 1.0;
         const double stress = ntr * grad + 45.0;
         const double pt = ntr * 2.0 * nm * 3.0 * nd + 75.0 + 54.0 + 24.0 * nd;
         return (order_ == 1 ? stress : nqf * stress) + nqf * pt;
      };
      const double faces = side(true) * static_cast<double>(n_int_sides_ + n_shared_sides_)
                           + side(false) * static_cast<double>(n_dir_sides_);
      return ((pass & kVolume) ? vol * ne_ : 0.0) + ((pass & kFaces) ? faces : 0.0);
   }
   /// Time the halo wait (StageTimers::DynHaloWait; null: off).
   void SetStageTimers(StageTimers *timers) { timers_ = timers; }
   /// Whether the team layout takes team_batch elements per block for K u
   /// (with_kv false) or K u + K_eta v: 8 at p = 2 with Kelvin-Voigt needs
   /// 49.9 KB of shared memory, over the 48 KB of a CUDA block.
   bool TeamBatchFits(int team_batch, bool with_kv) const
   {
      const int nm = with_kv ? 2 : 1;
      return (order_ == 1) ? sipg_mf::Team<1>::Fits(nm, team_batch)
                           : sipg_mf::Team<2>::Fits(nm, team_batch);
   }
   /// Kernel layout, and the elements per thread block of the team layout
   /// (2, 4 or 8; see TeamBatchFits).
   void SetLayout(Layout layout, int team_batch = 2)
   {
      MFEM_VERIFY(team_batch == 2 || team_batch == 4 || team_batch == 8,
                  "NitscheSipgOperator::SetLayout: team batch " << team_batch
                  << " is not 2, 4 or 8");
      layout_ = layout;
      team_batch_ = team_batch;
   }

   /// y = K u.
   void Mult(const Vector &u, Vector &y, int pass = kAll) const
   {
      Dispatch_<1>(u, nullptr, y, pass);
   }
   /// y = K u + K_eta v (needs the Kelvin-Voigt moduli).
   void Mult(const Vector &u, const Vector &v, Vector &y, int pass = kAll) const
   {
      MFEM_VERIFY(kv_, "NitscheSipgOperator::Mult(u, v, y): built without Kelvin-Voigt moduli");
      Dispatch_<2>(u, &v, y, pass);
   }

private:
   static constexpr int kGeo = sipg_mf::kGeo, kSide = sipg_mf::kSide;
   static constexpr int kNone = sipg_mf::kNone, kDirichlet = sipg_mf::kDirichlet;
   static constexpr int kGhost0 = sipg_mf::kGhost0;

   void BuildTables_(ParMesh &pm, const FiniteElement &fe)
   {
      const int nd = nd_;
      const IntegrationRule &irv = IntRules.Get(Geometry::TETRAHEDRON, 2 * (order_ - 1));
      const IntegrationRule &irf = IntRules.Get(Geometry::TRIANGLE, 2 * order_);
      nqv_ = irv.GetNPoints();
      nqf_ = irf.GetNPoints();
      MFEM_VERIFY(nqv_ == (order_ == 1 ? 1 : 4) && nqf_ == (order_ == 1 ? 3 : 6),
                  "NitscheSipgOperator: MFEM's rules have " << nqv_ << " tet and " << nqf_
                  << " triangle points; the kernels expect " << (order_ == 1 ? 1 : 4) << " and "
                  << (order_ == 1 ? 3 : 6));
      DenseMatrix dsh(nd, 3);
      Vector sh(nd);
      vw_.SetSize(nqv_);
      vdphi_.SetSize(nqv_ * nd * 3);
      for (int q = 0; q < nqv_; q++)
      {
         const IntegrationPoint &ip = irv.IntPoint(q);
         vw_(q) = ip.weight;
         fe.CalcDShape(ip, dsh);
         for (int k = 0; k < nd; k++)
            for (int r = 0; r < 3; r++) { vdphi_((q * nd + k) * 3 + r) = dsh(k, r); }
      }
      fw_.SetSize(nqf_);
      fphi_.SetSize(24 * nqf_ * nd);
      fdphi_.SetSize(24 * nqf_ * nd * 3);
      for (int q = 0; q < nqf_; q++) { fw_(q) = irf.IntPoint(q).weight; }
      IntegrationPointTransformation loc;
      for (int lf = 0; lf < 4; lf++)
         for (int ori = 0; ori < 6; ori++)
         {
            // face reference point -> element reference point, as the face
            // transformations' Loc1 / Loc2 (pmesh.cpp, GetLocalFaceTransformation)
            pm.GetLocalFaceTransformation(Element::TRIANGLE, Element::TETRAHEDRON, loc.Transf,
                                          64 * lf + ori);
            const int t = 6 * lf + ori;
            for (int q = 0; q < nqf_; q++)
            {
               IntegrationPoint eip;
               loc.Transform(irf.IntPoint(q), eip);
               fe.CalcShape(eip, sh);
               fe.CalcDShape(eip, dsh);
               for (int k = 0; k < nd; k++)
               {
                  fphi_((t * nqf_ + q) * nd + k) = sh(k);
                  for (int r = 0; r < 3; r++) { fdphi_(((t * nqf_ + q) * nd + k) * 3 + r) = dsh(k, r); }
               }
            }
         }
      for (auto *a : {&vw_, &vdphi_, &fw_, &fphi_, &fdphi_}) { a->UseDevice(true); }
   }

   void SetupHalo_(ParMesh &pm)
   {
      nfn_ = pm.GetNFaceNeighbors();
      if (nfn_ == 0) { return; }
      comm_ = fes_.GetComm();
      const Table &st = fes_.send_face_nbr_ldof, &rt = fes_.face_nbr_ldof;
      MFEM_VERIFY(st.Size() == nfn_ && rt.Size() == nfn_,
                  "NitscheSipgOperator: face-neighbour tables of " << st.Size() << " and "
                  << rt.Size() << " rows for " << nfn_ << " neighbours");
      send_off_.assign(st.GetI(), st.GetI() + nfn_ + 1);
      recv_off_.assign(rt.GetI(), rt.GetI() + nfn_ + 1);
      nbr_rank_.resize(nfn_);
      for (int fn = 0; fn < nfn_; fn++) { nbr_rank_[fn] = pm.GetFaceNbrRank(fn); }
      MFEM_VERIFY(recv_off_[nfn_] == fes_.GetFaceNbrVSize(),
                  "NitscheSipgOperator: halo of " << recv_off_[nfn_] << " values, face-neighbour "
                  "space of " << fes_.GetFaceNbrVSize());
      send_map_.SetSize(send_off_[nfn_]);
      const int *J = st.GetJ();
      for (int i = 0; i < send_map_.Size(); i++)
      {
         MFEM_VERIFY(J[i] >= 0, "NitscheSipgOperator: negative send DOF (not an L2 space)");
         send_map_[i] = J[i];
      }
      send_map_.GetMemory().UseDevice(true);
      reqs_.resize(2 * nfn_);
   }

public:
   // The functions holding the kernels are public: nvcc refuses an extended
   // __host__ __device__ lambda inside a private or protected member.

   /// Pack u (and v) on the send DOFs and post the messages.
   template <int NM>
   void HaloBegin_(const real_t *du, const real_t *dv, bool use_dev) const
   {
      if (nfn_ == 0) { return; }
      const int ns = send_map_.Size();
      send_buf_.UseDevice(true);
      halo_.UseDevice(true);
      send_buf_.SetSize(NM * ns);
      halo_.SetSize(NM * recv_off_[nfn_]);
      const int *smap = send_map_.Read(use_dev);
      real_t *sb = send_buf_.Write(use_dev);
      forall_switch(use_dev, ns, [=] MFEM_HOST_DEVICE (int i)
      {
         sb[NM * i] = du[smap[i]];
         if (NM == 2) { sb[NM * i + 1] = dv[smap[i]]; }
      });
      const bool gpu_mpi = use_dev && Device::GetGPUAwareMPI();
      const real_t *sp = gpu_mpi ? send_buf_.Read() : send_buf_.HostRead();
      real_t *rp = gpu_mpi ? halo_.Write() : halo_.HostWrite();
      if (gpu_mpi) { MFEM_STREAM_SYNC; }
      const int tag = 7323;
      for (int fn = 0; fn < nfn_; fn++)
      {
         MPI_Irecv(rp + NM * recv_off_[fn], NM * (recv_off_[fn + 1] - recv_off_[fn]),
                   MPITypeMap<real_t>::mpi_type, nbr_rank_[fn], tag, comm_, &reqs_[fn]);
         MPI_Isend(sp + NM * send_off_[fn], NM * (send_off_[fn + 1] - send_off_[fn]),
                   MPITypeMap<real_t>::mpi_type, nbr_rank_[fn], tag, comm_, &reqs_[nfn_ + fn]);
      }
   }

   void HaloEnd_() const
   {
      if (nfn_ == 0) { return; }
      ScopedStageTimer wait(timers_, StageTimers::DynHaloWait);
      MPI_Waitall(2 * nfn_, reqs_.data(), MPI_STATUSES_IGNORE);
   }

   template <int NM>
   void Dispatch_(const Vector &u, const Vector *v, Vector &y, int pass) const
   {
      MFEM_VERIFY(u.Size() == 3 * nsc_ && (!v || v->Size() == 3 * nsc_),
                  "NitscheSipgOperator::Mult: vector size " << u.Size() << " != " << 3 * nsc_);
      y.SetSize(3 * nsc_);
      if (order_ == 1) { Apply_<1, NM>(u, v, y, pass); }
      else { Apply_<2, NM>(u, v, y, pass); }
   }

   template <int P, int NM>
   void Apply_(const Vector &u, const Vector *v, Vector &y, int pass) const
   {
      using namespace sipg_mf;
      constexpr int ND = Tet<P>::ND, NQV = Tet<P>::NQV;
      const bool use_dev = u.UseDevice() || y.UseDevice() || (v && v->UseDevice());
      const real_t *du = u.Read(use_dev);
      const real_t *dv = (NM == 2) ? v->Read(use_dev) : nullptr;
      if (pass & kFaces) { HaloBegin_<NM>(du, dv, use_dev); }
      real_t *dy = y.Write(use_dev);
      const real_t *geo = geo_.Read(use_dev), *sd = side_d_.Read(use_dev);
      const int *si = side_i_.Read(use_dev);
      const real_t *vw = vw_.Read(use_dev), *vdphi = vdphi_.Read(use_dev);
      const real_t *fw = fw_.Read(use_dev), *fphi = fphi_.Read(use_dev),
                    *fdphi = fdphi_.Read(use_dev);
      const int nsc = nsc_;
      const real_t eps = eps_;
      const bool do_vol = (pass & kVolume) != 0, do_faces = (pass & kFaces) != 0;
      const bool team = (layout_ == kTeam) || (layout_ == kAuto && use_dev && P == 2);
      if (team)
      {
         KernelData kd{du, dv, nullptr, geo, sd, vw, vdphi, fw, fphi, fdphi, si, nullptr,
                       nsc, ne_, eps, do_vol, do_faces};
         TeamBatch_<P, NM, false>(use_dev, ne_, nullptr, kd, dy);
         if (!do_faces || nfn_ == 0) { return; }
         HaloEnd_();
         kd.halo = halo_.Read(use_dev);
         kd.gv = gvdof_.Read(use_dev);
         TeamBatch_<P, NM, true>(use_dev, shared_elems_.Size(), shared_elems_.Read(use_dev), kd,
                                 dy);
         return;
      }

      // kernel A: volume and the interior and Dirichlet sides
      forall_switch(use_dev, ne_, [=] MFEM_HOST_DEVICE (int e)
      {
         const real_t *ge = geo + kGeo * e;
         auto gete = [&](int c, int k, int m) {
            return (m == 0 ? du : dv)[c * nsc + e * ND + k];
         };
         real_t ye[3][ND];
         for (int c = 0; c < 3; c++)
            for (int k = 0; k < ND; k++) { ye[c][k] = 0.0; }
         if (do_vol)
         {
            for (int q = 0; q < NQV; q++)
            {
               const real_t *dph = vdphi + q * ND * 3;
               real_t Gr[3][3], G[3][3], S[3][3], Sr[3][3];
               RefGradient<P, NM>(gete, ge[12], dph, Gr);
               ToPhysical(Gr, ge, G);
               Stress(G, ge[10], ge[11], S);
               ToReference(S, ge, vw[q] * ge[9], Sr);
               for (int k = 0; k < ND; k++)
               {
                  const real_t d0 = dph[3 * k], d1 = dph[3 * k + 1], d2 = dph[3 * k + 2];
                  for (int c = 0; c < 3; c++) { ye[c][k] += Sr[c][0] * d0 + Sr[c][1] * d1 + Sr[c][2] * d2; }
               }
            }
         }
         if (do_faces)
         {
            for (int s = 0; s < 4; s++)
            {
               const int *ss = si + 3 * (4 * e + s);
               const int nbr = ss[0];
               if (nbr == kNone || nbr <= kGhost0) { continue; }
               const bool has_nbr = (nbr >= 0);
               const int h = has_nbr ? nbr : e;
               auto geth = [&](int c, int k, int m) {
                  return (m == 0 ? du : dv)[c * nsc + h * ND + k];
               };
               FaceSide<P, NM>(gete, ge, geth, geo + kGeo * h, has_nbr, sd + kSide * (4 * e + s),
                               ss[1], ss[2], eps, fw, fphi, fdphi, ye);
            }
         }
         for (int c = 0; c < 3; c++)
            for (int k = 0; k < ND; k++) { dy[c * nsc + e * ND + k] = ye[c][k]; }
      });

      if (!do_faces || nfn_ == 0) { return; }
      HaloEnd_();
      // kernel B: the shared sides, from the halo
      const real_t *halo = halo_.Read(use_dev);
      const int *gv = gvdof_.Read(use_dev), *sel = shared_elems_.Read(use_dev);
      const int ne = ne_;
      forall_switch(use_dev, shared_elems_.Size(), [=] MFEM_HOST_DEVICE (int i)
      {
         const int e = sel[i];
         const real_t *ge = geo + kGeo * e;
         auto gete = [&](int c, int k, int m) {
            return (m == 0 ? du : dv)[c * nsc + e * ND + k];
         };
         real_t ye[3][ND];
         for (int c = 0; c < 3; c++)
            for (int k = 0; k < ND; k++) { ye[c][k] = 0.0; }
         for (int s = 0; s < 4; s++)
         {
            const int *ss = si + 3 * (4 * e + s);
            if (ss[0] > kGhost0) { continue; }
            const int g = kGhost0 - ss[0];
            const int *gd = gv + 3 * ND * g;
            auto geth = [&](int c, int k, int m) { return halo[NM * gd[c * ND + k] + m]; };
            FaceSide<P, NM>(gete, ge, geth, geo + kGeo * (ne + g), true, sd + kSide * (4 * e + s),
                            ss[1], ss[2], eps, fw, fphi, fdphi, ye);
         }
         for (int c = 0; c < 3; c++)
            for (int k = 0; k < ND; k++) { dy[c * nsc + e * ND + k] += ye[c][k]; }
      });
   }

   template <int P, int NM, bool GHOST>
   void TeamBatch_(bool use_dev, int n, const int *sel, const sipg_mf::KernelData &kd,
                   real_t *dy) const
   {
      switch (team_batch_)
      {
         case 2: Team_<P, NM, 2, GHOST>(use_dev, n, sel, kd, dy); break;
         case 4: Team_<P, NM, 4, GHOST>(use_dev, n, sel, kd, dy); break;
         default:
            // not instantiated where the block's arrays do not fit
            if constexpr (sipg_mf::Team<P>::Fits(NM, 8))
            {
               Team_<P, NM, 8, GHOST>(use_dev, n, sel, kd, dy);
            }
            else
            {
               MFEM_ABORT("NitscheSipgOperator: team batch 8 at p = " << P << (NM == 2 ? " with "
                          "Kelvin-Voigt" : "") << " needs " << sipg_mf::Team<P>::Bytes(NM, 8)
                          << " bytes of shared memory per block (48 KB max); use 2 or 4");
            }
            break;
      }
   }

   /// The team layout: TX threads per element, BZ elements per block.
   /// sel: the elements (kernel B), null for 0 .. n - 1 (kernel A).
   template <int P, int NM, int BZ, bool GHOST>
   void Team_(bool use_dev, int n, const int *sel, const sipg_mf::KernelData &kd,
              real_t *dy) const
   {
      using namespace sipg_mf;
      constexpr int NV = Team<P>::NV, NPT = Team<P>::NPT, TX = Team<P>::TX;
      const KernelData a = kd;
      ForallWrap<2>(use_dev, n, [=] MFEM_HOST_DEVICE (int i)
      {
         MFEM_SHARED real_t s_ue[BZ][NM * NV], s_uh[BZ][4 * NM * NV];
         MFEM_SHARED real_t s_sp[BZ][NPT][12], s_tr[BZ][8 * Team<P>::NQT][3];
         const int tz = MFEM_THREAD_ID(z);
         TeamElement<P, NM, GHOST>(sel ? sel[i] : i, a, dy, s_ue[tz], s_uh[tz], s_sp[tz],
                                   s_tr[tz]);
      }, TX, 1, BZ);
   }

private:
   ParFiniteElementSpace &fes_;
   real_t eps_ = -1.0;
   bool kv_ = false;
   int ne_ = 0, ng_ = 0, order_ = 0, nd_ = 0, nsc_ = 0, nqv_ = 0, nqf_ = 0;
   long long n_int_sides_ = 0, n_shared_sides_ = 0, n_dir_sides_ = 0;
   Vector geo_, side_d_;                 ///< per element (and ghost); per side
   Array<int> side_i_, gvdof_, shared_elems_;
   Vector vw_, vdphi_, fw_, fphi_, fdphi_;   ///< reference tables
   // halo
   int nfn_ = 0;
   MPI_Comm comm_ = MPI_COMM_NULL;
   std::vector<int> send_off_, recv_off_, nbr_rank_;
   Array<int> send_map_;
   mutable Vector send_buf_, halo_;
   mutable std::vector<MPI_Request> reqs_;
   StageTimers *timers_ = nullptr;
   Layout layout_ = kAuto;
   int team_batch_ = 2;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_SIPG_OPERATOR_HPP
