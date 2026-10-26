# TPV205-2D reference: PyLith, 200 m

SCEC TPV205-2D (the plane-strain slice of TPV205 at 7.5 km depth) from the SCEC
dynamic rupture code verification uploads: Brad Aagaard, PyLith 1.7.0a
(scecdynrup branch), element size 200 m, one sample every 0.06 s from 0 to 12 s
(201 rows), dated 2012-02-01 in the headers.  Chunhui Zhao placed the files in
the Mac clone in October 2026 (the download URL was not recorded); they are
here byte for byte (sha256 below).

- On-fault stations `tpv205_pylith2d_x2_<x>_x3_7.5.txt`, x = -12, -7.5, -4.5,
  0, 4.5, 7.5, 12 km: `t h-slip h-slip-rate h-shear-stress v-slip v-slip-rate
  v-shear-stress` (m, m/s, MPa).  There is no normal-stress column.  The
  headers place the +-4.5 and +-7.5 km stations at +-4.4 and +-7.6 km (200 m
  mesh nodes).  The first row (t = 0) has shear stress 0; from t = 0.06 s it
  holds the initial stress.  At x = 0, h-slip is 4.5e-4 m at t = 0.
- Off-fault stations `tpv205_pylith2d_body_3.0_x2_+-12.0_x3_7.5.txt`, 3 km off
  the fault at x = +-12 km: `t h-disp h-vel v-disp v-vel n-disp n-vel`.
- PyLith is the only TPV205-2D code here; there is no second code to
  cross-check against.

Read by `scripts/plot_tpv205_2d_stations.py` (the TPV205-2D gold's plots).

sha256:

```
b1213840a1cef05faf45af27bf1c3591cd563b0a08541b5ba820a374002cc20d  tpv205_pylith2d_body_3.0_x2_-12.0_x3_7.5.txt
f1c9f647716f9d9e86f041bcdee7759fa93b67f5b9c8f236b686971a5c132063  tpv205_pylith2d_body_3.0_x2_12.0_x3_7.5.txt
cf4e51b11c28366c5bc289faad274b77bf74615d1286d84e883eecfa9be6bb33  tpv205_pylith2d_x2_-12.0_x3_7.5.txt
4933c17365ccff07d6d97d03f9df5ed59d249efbbc9ec82fb4cc1ae298e2f27d  tpv205_pylith2d_x2_-4.5_x3_7.5.txt
69bf8ee58d17df2194913b3837cba159d1aefbca964f0bdd0d4a2dd9222422b5  tpv205_pylith2d_x2_-7.5_x3_7.5.txt
482202cbb70de336a0a8038d04463cc5e72b053b360f34cb44f6e31538b7f0bc  tpv205_pylith2d_x2_0_x3_7.5.txt
d6a0c44fa104660c8b2de683dc1d66fa09711b2d922b29d186bd1ef3cef9ef41  tpv205_pylith2d_x2_12.0_x3_7.5.txt
2fc8314a47cd74f3edbc7866dee1d6c272c2d8351f34a5d62c6bf929fe841f2e  tpv205_pylith2d_x2_4.5_x3_7.5.txt
916f7ea1c5ae5079bd85cd57fe8fd2f99fc981aed2f18c776f5569bdfe869415  tpv205_pylith2d_x2_7.5_x3_7.5.txt
```
