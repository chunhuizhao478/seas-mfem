// INVESTIGATION PROBE (dip-drift debug, 2026-06-06): does any FAULT-FRAME error
// manufacture DIP-channel traction (local SXY = n.sigma.t1) from a pure
// strike-slip global stress?  Uses the REAL GodunovFlux::BuildRotationInverse.
// This is a diagnostic, not a functional change. Builds standalone like
// test_canonical_rotation_pure_strikeslip.cpp (links only godunov_flux.o).
#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include <cmath>
#include <cstdio>
using namespace mfem;
using namespace mfem::seas;

static void run(const char *label, const real_t n[3], const real_t t1[3],
                const real_t t2[3], const Vector &sig_glob)
{
   DenseMatrix Tinv;
   GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv);
   Vector loc(NUM_STATE); Tinv.Mult(sig_glob, loc);
   std::printf("  %-34s dip(SXY)=% .6e   strike(SXZ)=% .6e\n",
               label, loc(SXY), loc(SXZ));
}
static void unit3(real_t v[3]){real_t L=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);v[0]/=L;v[1]/=L;v[2]/=L;}
static void cross(const real_t a[3],const real_t b[3],real_t o[3]){
   o[0]=a[1]*b[2]-a[2]*b[1];o[1]=a[2]*b[0]-a[0]*b[2];o[2]=a[0]*b[1]-a[1]*b[0];}
// FaultBasis-style: strike=up x n, dip=strike x n; t1=dip,t2=strike; canon negate if n.ref<0
static void faultbasis(const real_t nraw[3], real_t n[3], real_t t1[3], real_t t2[3])
{
   const real_t up[3]={0,0,1}, ref[3]={0,-1,0};
   n[0]=nraw[0];n[1]=nraw[1];n[2]=nraw[2]; unit3(n);
   real_t strike[3],dip[3]; cross(up,n,strike); unit3(strike); cross(strike,n,dip); unit3(dip);
   for(int i=0;i<3;i++){t1[i]=dip[i];t2[i]=strike[i];}
   if(n[0]*ref[0]+n[1]*ref[1]+n[2]*ref[2]<0){for(int i=0;i<3;i++){n[i]=-n[i];t1[i]=-t1[i];t2[i]=-t2[i];}}
}

int main()
{
   const real_t tau=75.0e6;
   Vector sig(NUM_STATE); sig=0.0; sig(SXY)=tau;   // pure strike-slip global stress
   std::printf("\n=== Frame dip-leak probe (real BuildRotationInverse) ===\n");
   std::printf("Pure strike-slip sigma_xy=%.3e; ideal dip channel must be 0, strike=-tau\n", tau);

   real_t n[3]={0,-1,0}, t1[3]={0,0,-1}, t2[3]={1,0,0};
   run("[1] ideal frame", n,t1,t2, sig);

   { real_t a[3]={0, 1,0},b[3]={0,0,-1},c[3]={1,0,0}; run("[2] flip n", a,b,c, sig); }
   { real_t a[3]={0,-1,0},b[3]={0,0, 1},c[3]={1,0,0}; run("[2] flip t1(dip)", a,b,c, sig); }
   { real_t a[3]={0,-1,0},b[3]={0,0,-1},c[3]={-1,0,0}; run("[2] flip t2(strike)", a,b,c, sig); }

   for (real_t e : {1e-2, 1e-1}) {
      real_t raw[3]={e,-1,0}, fn[3],f1[3],f2[3]; faultbasis(raw,fn,f1,f2);
      char L[64]; std::snprintf(L,64,"[3] normal tilt x eps=%.0e",e); run(L,fn,f1,f2,sig);
      real_t raw2[3]={0,-1,e}; faultbasis(raw2,fn,f1,f2);
      std::snprintf(L,64,"[3] normal tilt z eps=%.0e",e); run(L,fn,f1,f2,sig);
   }
   for (real_t e : {1e-2, 1e-1}) {
      real_t t1i[3]={0,0,-1},t2i[3]={1,0,0};
      real_t a[3]={0,-1,0};
      real_t b[3]={std::cos(e)*t1i[0]+std::sin(e)*t2i[0],std::cos(e)*t1i[1]+std::sin(e)*t2i[1],std::cos(e)*t1i[2]+std::sin(e)*t2i[2]};
      real_t c[3]={-std::sin(e)*t1i[0]+std::cos(e)*t2i[0],-std::sin(e)*t1i[1]+std::cos(e)*t2i[1],-std::sin(e)*t1i[2]+std::cos(e)*t2i[2]};
      char L[80]; std::snprintf(L,80,"[4] in-plane rot eps=%.0e (=-eps*tau=%.2e)",e,-e*tau); run(L,a,b,c,sig);
   }
   std::printf("Ideal frame, global sigma carries a small DIP-COUPLED component:\n");
   { Vector s=sig; s(SYZ)+=1e-3*tau; run("[5] +sigma_yz=1e-3*tau", n,t1,t2, s); }
   { Vector s=sig; s(SXZ)+=1e-3*tau; run("[5] +sigma_xz=1e-3*tau", n,t1,t2, s); }
   std::printf("\nCONCLUSION: only [4] in-plane rotation (structurally prevented by\n"
               "FaultBasis up-pinning) and [5] a genuine global sigma_yz leak into dip.\n"
               "Sign flips [2] and normal tilts [3] do NOT. => frame is not the dip seed.\n");
   return 0;
}
