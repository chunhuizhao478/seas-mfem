#!/bin/zsh
# Rebuild every view carrying the per-ELEMENT resolution fields:
#   f_min_p3, f_min_p5, Vs, edge_min, edge_max, vs_over_dx  (+ bc on surfaces)
#
#   f_min = (p/4) * Vs/dx , dx = element MAX edge, Vs nearest-grid at the
#   barycentre from MUSCAL.  p3 -> 0.75*Vs/dx (gate 0.6667 for 0.5 Hz),
#   p5 -> 1.25*Vs/dx (gate 0.8 for 1 Hz).
#
# `volume` keeps geometry/connect IN the .puml.h5 and writes only a field
# sidecar, so the big meshes are not duplicated.
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
R=$A/meshing_shakeoutbox_gate/results
V=$A/meshing_shakeoutbox_gate/view
INT=$R/safalt_0d5Hz_p3_deep40km_shakeoutbox_muscal.puml.h5
HVY=$R/safalt_fb200_1Hz_p5_shakeoutbox_muscal.puml.h5
cd $V

echo "########## intermediate surface"; $PY $CODE/mesh_xdmf_fields.py --mesh $INT --mode surface --out $V/intermediate_surface.xdmf; echo "rc=$?"
echo "########## heavy fault";          $PY $CODE/mesh_xdmf_fields.py --mesh $HVY --mode fault   --out $V/heavy_fault.xdmf;          echo "rc=$?"
echo "########## heavy surface";        $PY $CODE/mesh_xdmf_fields.py --mesh $HVY --mode surface --out $V/heavy_surface.xdmf;        echo "rc=$?"
echo "########## intermediate volume";  $PY $CODE/mesh_xdmf_fields.py --mesh $INT --mode volume  --out $V/intermediate_volume.xdmf;  echo "rc=$?"
echo "########## heavy volume";         $PY $CODE/mesh_xdmf_fields.py --mesh $HVY --mode volume  --out $V/heavy_volume.xdmf;         echo "rc=$?"
rm -f $V/intermediate_full.xdmf $V/heavy_full.xdmf   # superseded by *_volume
echo ALL DONE
