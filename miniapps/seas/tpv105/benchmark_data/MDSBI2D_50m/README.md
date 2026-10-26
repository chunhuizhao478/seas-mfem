# TPV105-2D reference: MDSBI, 50 m

SCEC TPV105-2D (rate-and-state slip law with strong rate weakening and thermal
pressurization, plane strain, a slice at 7.5 km depth) from the SCEC dynamic
rupture code verification uploads: Eric M. Dunham, MDSBI 4.1.6, element size
50 m, one sample every 0.005 s from 0 to 15 s (3001 rows), dated 2011-02-10 in
the headers.  Chunhui Zhao placed the files in the Mac clone in October 2026
(the download URL was not recorded); they are here byte for byte (sha256 below).

- On-fault stations `tpv105_mdsbi2d_x2_<x>_x3_7.5.txt`, x = -18, -15, -9, 0,
  9, 15, 18 km: `t h-slip h-slip-rate h-shear-stress v-slip v-slip-rate
  v-shear-stress n-stress psi temperature pressure` (m, m/s, MPa; n-stress =
  effective normal stress, compression positive; temperature in K; pore
  pressure in MPa).
- Checked against SPECFEM2D_50m (10-25 benchmark check): arrivals within
  0.03 s, final slip within 3.5 % except x = -18 km (1.255 against 1.100 m),
  peak temperature within 2.5 %, pore pressure within 1.2 MPa.

Read by `scripts/plot_tpv105_2d_stations.py` (the TPV105-2D gold's plots,
where it is the reference the misfits are measured against).

sha256:

```
0a5ef1fdbb9b2584f012d66a61886b2784e39e0b327be2561b76981635085a55  tpv105_mdsbi2d_x2_-15.0_x3_7.5.txt
55ede1ed463d5f236241e51503e21f4dd2e2fd7603a2f553b250c24473182218  tpv105_mdsbi2d_x2_-18.0_x3_7.5.txt
ccd94d2f9471ff84a7fdccb77a16269838617d6a6eaae4b1888a45d64937684e  tpv105_mdsbi2d_x2_-9.0_x3_7.5.txt
80191baa9f0b46df491382f58b0f9301448d7651a039800f601d9181328ff577  tpv105_mdsbi2d_x2_0_x3_7.5.txt
b3f034372106ade3f7a3a281e0a48d09b5fdfeaf5f963e61fd32f8acbd414e2e  tpv105_mdsbi2d_x2_15.0_x3_7.5.txt
34a929edbc8159a677a86d0eb17d6a8c60c56eadfe0714f6498d545f5438fe98  tpv105_mdsbi2d_x2_18.0_x3_7.5.txt
35c025ac20d8aaeaf82c070420016a3140bd6aaa2d761af4c61312ed345238aa  tpv105_mdsbi2d_x2_9.0_x3_7.5.txt
```
