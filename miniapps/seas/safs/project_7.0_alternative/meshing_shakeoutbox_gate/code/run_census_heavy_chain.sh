#!/bin/zsh
# The two big censuses, sequentially -- each holds ~0.6 GB of geometry plus its
# chunk, and each streams a ~6 GB file, so do not race them for page cache.
CODE=${0:A:h}
A=/Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/project_7.0_alternative
cd $A/meshing_shakeoutbox_gate/build_tmp
echo "### heavy shakeoutbox product"
/usr/bin/time -l $CODE/run_census.sh heavy > census_heavy.log 2>&1
echo "heavy rc=$?"
echo "### 1 Hz p5 parent (candidate heavy parent)"
/usr/bin/time -l $CODE/run_census.sh parent1hz > census_parent1hz.log 2>&1
echo "parent1hz rc=$?"
echo "### refine2 parent (current heavy parent)"
/usr/bin/time -l $CODE/run_census.sh parentref2 > census_parentref2.log 2>&1
echo "parentref2 rc=$?"
echo ALL DONE
