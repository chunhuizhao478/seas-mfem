# TPV105-2D reference: SPECFEM2D (spectral elements), 50 m

SCEC TPV105-2D from the SCEC dynamic rupture code verification uploads:
Kaneko, Lapusta and Ampuero, SPECFEM3D (se) 0.3.1 in its 2D form, element size
50 m with 5 GLL nodes, one sample every 0.0025 s from 0.0025 s to 12.5 s (5000
rows; the benchmark asks for 15 s), dated 2011-02-20 in the headers.  Chunhui
Zhao placed the files in the Mac clone in October 2026 (the download URL was
not recorded); they are here byte for byte (sha256 below).

- On-fault stations `tpv105_specfem2d_x2_<x>_x3_7.5.txt` (header location
  faultst<x in 0.1 km>dp075), x = -18, -15, -9, 0, 9, 15, 18 km: the same 11
  columns as MDSBI2D_50m (`t h-slip h-slip-rate h-shear-stress v-slip
  v-slip-rate v-shear-stress n-stress psi temperature pressure`); the header
  lines start with a space before the #.
- Because the record ends at 12.5 s, the TPV105-2D comparisons use the window
  0 to 12.5 s.

Read by `scripts/plot_tpv105_2d_stations.py` (the TPV105-2D gold's plots).

sha256:

```
80c8480615faa3f30174d15de4f20a7d82f7a3ef7842ed784df8a9f293fbec9d  tpv105_specfem2d_x2_-15.0_x3_7.5.txt
bb2e5e11a36cc43433900311bc869cce15a2f9a41c26663d827d135d7da99117  tpv105_specfem2d_x2_-18.0_x3_7.5.txt
93e8591335d4d42b87919a0d2ce9e9b53f280202f3ca932e7aef745766a880a1  tpv105_specfem2d_x2_-9.0_x3_7.5.txt
430928f31f470b220027bd68f2a996440c8d7d5e86c4058a506ff04fd1305c77  tpv105_specfem2d_x2_0_x3_7.5.txt
0c4990ad855e2c690117da9b33890fb3813425122fc7092ec0b7039d5c6b6e1f  tpv105_specfem2d_x2_15.0_x3_7.5.txt
79f051ad30165f156ddaa1475d50d563cced08214a0ad39105e0fff95aa67a17  tpv105_specfem2d_x2_18.0_x3_7.5.txt
c52ce5c90df0142c8534fe57e3a8f25c6c8a32772355aba889d4f39e3428a028  tpv105_specfem2d_x2_9.0_x3_7.5.txt
```
