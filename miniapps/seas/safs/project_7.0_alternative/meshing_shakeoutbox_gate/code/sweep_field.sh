#!/bin/zsh
# beta and hgrad are the two knobs that move the collar's price by 3x.
# beta = barycentre-depth / max-edge; MEASURED on the shipped collar as
# p1 0.083 / p10 0.124 / median 0.294.  The FINE branch is beta-independent
# (small cells land in the z=0 bin whatever beta is), but the COARSE branch is
# not, which is exactly why a uniform 2,500 m collar still fails 91,784 cells.
PY=/Users/chunhuizhao/miniforge/envs/pythonenv/bin/python
CODE=${0:A:h}
cd /Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative/meshing_shakeoutbox_gate/build_tmp
for gate in 0.8 0.6667; do
  for beta in 0.294 0.124; do
    for hgrad in 1.3 1.5 2.0; do
      echo "=================== gate $gate  beta $beta  hgrad $hgrad"
      $PY $CODE/field3d.py --gate $gate --beta $beta --hgrad $hgrad 2>&1 \
        | grep -E "measured build factor|z=0 |collar volume"
    done
  done
done
