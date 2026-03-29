#!/bin/bash
# Setup directories and build Tandem for p2, p4, p6 on Frontera
# Run this on the Frontera login node

BASE=/scratch2/10024/zhaochun/tandem
SRC=$BASE/tandem

module load gcc/9.1.0
module load impi/19.0.9
module unload cmake/4.1.1 2>/dev/null
module load cmake/3.31.5

# ---- Build p2 ----
echo "=== Building Tandem p2 ==="
mkdir -p $BASE/build_3d_2p_staging && cd $BASE/build_3d_2p_staging
cmake $SRC \
  -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_PREFIX_PATH=$BASE/../petsc/build \
  -DDOMAIN_DIMENSION=3 -DPOLYNOMIAL_DEGREE=2
make -j20

# ---- Build p4 ----
echo "=== Building Tandem p4 ==="
mkdir -p $BASE/build_3d_4p_staging && cd $BASE/build_3d_4p_staging
cmake $SRC \
  -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_PREFIX_PATH=$BASE/../petsc/build \
  -DDOMAIN_DIMENSION=3 -DPOLYNOMIAL_DEGREE=4
make -j20

# ---- Build p6 ----
echo "=== Building Tandem p6 ==="
mkdir -p $BASE/build_3d_6p_staging && cd $BASE/build_3d_6p_staging
cmake $SRC \
  -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_PREFIX_PATH=$BASE/../petsc/build \
  -DDOMAIN_DIMENSION=3 -DPOLYNOMIAL_DEGREE=6
make -j20

# ---- Setup run directories ----
echo "=== Setting up run directories ==="

# p2 1000m
mkdir -p $BASE/bp5_1000m_p2
cp $SRC/examples/tandem/3d/bp5.toml $BASE/bp5_1000m_p2/
cp $SRC/examples/tandem/3d/bp5.lua  $BASE/bp5_1000m_p2/
echo "  bp5_1000m_p2: scp bp5_tandem_1000m.msh -> bp5.msh"

# p4 2500m
mkdir -p $BASE/bp5_2500m_p4
cp $SRC/examples/tandem/3d/bp5.toml $BASE/bp5_2500m_p4/
cp $SRC/examples/tandem/3d/bp5.lua  $BASE/bp5_2500m_p4/
echo "  bp5_2500m_p4: scp bp5_tandem_2500m.msh -> bp5.msh"

# p6 4000m
mkdir -p $BASE/bp5_4000m_p6
cp $SRC/examples/tandem/3d/bp5.toml $BASE/bp5_4000m_p6/
cp $SRC/examples/tandem/3d/bp5.lua  $BASE/bp5_4000m_p6/
echo "  bp5_4000m_p6: scp bp5_tandem_4000m.msh -> bp5.msh"

echo "=== Done. Now scp mesh files and submit jobs. ==="
