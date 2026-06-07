// INVESTIGATION PROBE (dip-drift, 2026-06-06): symmetry-equivariance of the
// Godunov fluxes under the two fault-plane symmetries.  A correct elastodynamic
// flux must commute with each mirror M:  M . Flux(n, Q[, Qbg]) == Flux(Mn, M.Q[, M.Qbg]).
// The dip slip is ANTI-symmetric in along-strike x => it breaks the STRIKE-REVERSAL
// (x->-x) symmetry.  Tests Interior (off-fault upwind), Central, and the
// FREE-SURFACE flux (drift peaks at the surface) under x->-x.  Diagnostic only.
#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include <cmath>
#include <cstdio>
using namespace mfem;
using namespace mfem::seas;

// STRIKE-REVERSAL x->-x: x-ODD state components flip: VX(6), SXY(3), SXZ(5).
static inline real_t sx(int c){ return (c==VX||c==SXY||c==SXZ) ? -1.0 : 1.0; }
static void Mx(const real_t *Q, real_t *MQ){ for(int c=0;c<NUM_STATE;c++) MQ[c]=sx(c)*Q[c]; }
static void norm3(real_t v[3]){real_t L=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);v[0]/=L;v[1]/=L;v[2]/=L;}
static real_t l1(const real_t*a,const real_t*b){real_t s=0;for(int c=0;c<NUM_STATE;c++)s+=std::abs(a[c]-b[c]);return s;}

int main()
{
   const real_t rho=2670.0, cp=6000.0, cs=3464.0;
   GodunovFlux flux(rho, cp, cs);

   real_t Qs[NUM_STATE], Qn[NUM_STATE], Qbg[NUM_STATE];
   for(int c=0;c<NUM_STATE;c++){ Qs[c]=1.0+0.37*c; Qn[c]=-0.5+0.21*c; Qbg[c]=0.13*c-0.4; }
   real_t MQs[NUM_STATE],MQn[NUM_STATE],MQbg[NUM_STATE]; Mx(Qs,MQs);Mx(Qn,MQn);Mx(Qbg,MQbg);

   std::printf("\n=== FLUX STRIKE-REVERSAL (x->-x) EQUIVARIANCE PROBE ===\n");
   std::printf("residual = | M.Flux(n,Q) - Flux(Mn, M.Q) |.  0 = symmetric.\n\n");

   // --- (A) off-fault interior + central on strike-tilted faces ---
   std::printf("(A) interior/central, face tilted along strike (n_x varies, n_z=0.3):\n");
   std::printf("    %6s | %13s %13s\n","n_x","UPWIND |R|","CENTRAL |R|");
   for (real_t nx : {0.0,0.25,0.5,0.75,1.0}) {
      real_t n[3]={nx,1.0,0.3}; norm3(n); real_t Mn[3]={-n[0],n[1],n[2]};
      real_t Fi[NUM_STATE],MFi[NUM_STATE],Fim[NUM_STATE]; flux.Interior(n,Qs,Qn,Fi); Mx(Fi,MFi); flux.Interior(Mn,MQs,MQn,Fim);
      real_t Fc[NUM_STATE],MFc[NUM_STATE],Fcm[NUM_STATE]; flux.Central (n,Qs,Qn,Fc); Mx(Fc,MFc); flux.Central (Mn,MQs,MQn,Fcm);
      std::printf("    %6.2f | %13.4e %13.4e\n", nx, l1(MFi,Fim), l1(MFc,Fcm));
   }

   // --- (B) FREE SURFACE (horizontal face, n=(0,0,1)); drift peaks here ---
   // n_x = 0 so Mn = n; test M.FS(n,Q,Qbg) vs FS(n, M.Q, M.Qbg).
   std::printf("\n(B) FREE-SURFACE flux (n=(0,0,1), background-aware gamma-mirror):\n");
   real_t nfs[3]={0.0,0.0,1.0};
   real_t Ft[NUM_STATE],MFt[NUM_STATE],Ftm[NUM_STATE];
   flux.FreeSurfaceGodunovTotal(nfs, Qs, Qbg, Ft);  Mx(Ft,MFt);
   flux.FreeSurfaceGodunovTotal(nfs, MQs, MQbg, Ftm);
   std::printf("    FreeSurfaceGodunovTotal |R| = %.4e   R[VZ]=%.4e  R[SYZ]=%.4e\n",
               l1(MFt,Ftm), MFt[VZ]-Ftm[VZ], MFt[SYZ]-Ftm[SYZ]);
   real_t Gt[NUM_STATE],MGt[NUM_STATE],Gtm[NUM_STATE];
   flux.FreeSurfaceTotal(nfs, Qs, Qbg, Gt);  Mx(Gt,MGt);
   flux.FreeSurfaceTotal(nfs, MQs, MQbg, Gtm);
   std::printf("    FreeSurfaceTotal        |R| = %.4e   R[VZ]=%.4e  R[SYZ]=%.4e\n",
               l1(MGt,Gtm), MGt[VZ]-Gtm[VZ], MGt[SYZ]-Gtm[SYZ]);

   std::printf("\nNonzero |R| => that flux breaks strike-reversal symmetry and can\n"
               "radiate the anti-symmetric-in-x dip seed.\n");
   return 0;
}
