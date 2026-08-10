#!/bin/zsh
# Generate the three ParaView views for every PREFERRED shakeout-domain mesh.
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
R=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_preferred
cd $CODE
for t in small intermediate heavy; do
  M=$R/meshing_shakeout_domain_$t/results/safpref_${t}_shakeout.puml.h5
  D=$R/meshing_shakeout_domain_$t/view
  mkdir -p $D
  echo "########## PREF $t ##########"
  $PY mesh_xdmf.py --mesh $M --mode wrap    --out $D/${t}_full.xdmf
  $PY mesh_xdmf.py --mesh $M --mode fault   --out $D/${t}_fault.xdmf
  $PY mesh_xdmf.py --mesh $M --mode surface --out $D/${t}_surface.xdmf
done
