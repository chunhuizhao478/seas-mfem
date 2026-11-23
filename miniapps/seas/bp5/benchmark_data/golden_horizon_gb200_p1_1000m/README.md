# gpu_green_p1_1500yr_g1_job69011

Produced from commit `3f4053845fdf87ad833326fffecb339f3405fb96` (branch `gpu-qd-bp5`) of https://github.com/chunhuizhao478/seas-mfem:
https://github.com/chunhuizhao478/seas-mfem/commit/3f4053845fdf87ad833326fffecb339f3405fb96

The checkout had no uncommitted changes in tracked files, so these results come from that commit exactly.

The job wrote this README with `jobs/vista/gold_readme.sh`, which was written for Vista. The cluster, the
run line and the reproduce steps were corrected by hand on 2026-10-25, before the commit; see
"Corrections" at the end. Every other file in this folder is as the job wrote it.

## Job

- Cluster: Horizon (c105-035.horizon.tacc.utexas.edu), Slurm job 69011 (h_full_cycle), started 2026-10-25T01:36:13Z
- Partition debug, account EAR26024, 1 node(s), 1 task(s), one GB200 GPU
- Job script: `miniapps/seas/jobs/horizon/h_full_cycle.sbatch` at that commit. It sources
  `jobs/horizon/h_env.sh` and runs `jobs/vista/g5_gpu_full_cycle.sbatch`.
- Submit line: `sbatch --parsable --export=ALL,QD_PATH=green jobs/horizon/h_full_cycle.sbatch`
- Settings from the environment: ORDER=1, NP=1, TFINAL_YR=1500, CUDA_MODULE=cuda/13.1, G1_GOLDEN_ROOT=/scratch/10024/zhaochun/seas-gpu-goldens/g0_vista

## Inputs

- Run line (`jobs/vista/g5_gpu_full_cycle.sbatch` with QD_PATH=green, from `miniapps/seas`; run.log ends
  with the same options in PETSc's list of options it did not use):

```
ibrun -n 1 -o 0 ./seas_bp5_full --mesh bp5/mesh/reference/bp5_tandem_exact.msh \
   --mesh-scale 1000 --dg-method IP --tandem-time-stepping --solver mumps --petsc-ts \
   --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg \
   --order 1 --device cuda --solver cudss \
   --qd-operators matrix --qd-green --qd-green-cache /scratch/10024/zhaochun/seas-qd-green/bp5_prod_p1.bin \
   --stage-timers 5000 --tfinal 4.733640e+10 --checkpoint-interval 5000 \
   --output-dir bp5/results_vista_cycle/gpu_green_p1_1500yr_g1_job69011 --output-prefix bp5_cycle_p1
```

  The first four lines are `G0_V63_FLAGS` (`jobs/vista/g0_env.sh`), the deck of the CPU gold 1074099. The
  later `--solver cudss` replaces `--solver mumps`. `--qd-green` makes each stage one product with the dense
  Green's function G on the GPU (cuBLAS dgemv). cuDSS solves the displacement only when a checkpoint
  writes it (62 solves in the run; `stage_timers_final.txt`, solve calls minus stages).
- G: `bp5_prod_p1.bin`, 37,462,176,512 bytes (83,808 x 55,872 doubles), built on the CPU with MUMPS on 144
  Grace ranks by `jobs/horizon/h_green_cache.sbatch` (job 68897, same commit). The run read it in 214 s, moved
  it to the GPU and checked it against the solve path at t = 100 yr (`qd_green.txt`: max relative difference
  3.2e-15 shear, 5.3e-14 normal).

- Meshes at that commit (paths relative to `miniapps/seas`), sha256:

```
f413b957cf1c6d02d2c84b65fc05d1f99232301ed8083cb7d8533c9fa442c3b4  bp5/mesh/reference/bp5_tandem_exact.msh
34dd17748445d1982f867389023b3873e401380baf18a7a0c2c97591e3ff6121  bp5/mesh/reference/bp5_tandem_500m.msh
c2bffe326e3e7c8df0be73c97fed734225d8968b1f8eb776fbd9f7055bd23dbc  bp5/mesh/reference/bp5_tandem_250m.msh
```

## Build

- Binary: `miniapps/seas/seas_bp5_full`, sha256 9b54b311f2195794b20f3b1eefdab2dabaa5715c3242b96d71e6c87425ebc0df
- Built with `MODE=cuda bash build_horizon.sh` (PETSc and its packages, MFEM, then seas). The job rebuilt
  seas at its own commit (the last 'seas' stage below); PETSc and MFEM are from the build at d079573a.
  The build's "dirty files" also counts untracked files (job logs and result folders in the clone).
  What the build used:

```
# build_horizon.sh MODE=cuda stage 'fetch' on login1.horizon.tacc.utexas.edu at 2026-10-22T03:40:56Z
source: 4e4e52ccc9f87b753349a5eeedb99aeb0b104c4e (gpu-qd-bp5), dirty files: 0
modules: ucc/1.9.0 ucx/1.22.0 xalt/3.2.4 TACC gcc/15.3.0 cuda/13.3 openmpi/5.0.11 nvpl/26.9 cmake/4.4.0 
mpicc: gcc (GCC) 15.3.0
mpif90: GNU Fortran (GCC) 15.3.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
nvcc: Cuda compilation tools, release 13.3, V13.3.73
gpu arch: sm_100
cudss: /scratch/10024/zhaochun/seas-gpu-deps/cudss-0.7.1.6 (CUDSS_VERSION 0.7.1.6)
petsc: 3.24.6 at /scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda
# build_horizon.sh MODE=cuda stage 'deps' on c109-028.horizon.tacc.utexas.edu at 2026-10-23T16:39:31Z
source: 0a31efaafd0a68abf684b5d2111725f2c492615b (gpu-qd-bp5), dirty files: 6
modules: ucx/1.22.0 xalt/3.2.4 TACC gcc/15.3.0 cuda/13.1 ucc/1.9.0 openmpi/5.0.11 nvpl/26.9 cmake/4.4.0 
mpicc: gcc (GCC) 15.3.0
mpif90: GNU Fortran (GCC) 15.3.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
nvcc: Cuda compilation tools, release 13.1, V13.1.115
gpu arch: sm_100
cudss: /scratch/10024/zhaochun/seas-gpu-deps/cudss-0.7.1.6 (CUDSS_VERSION 0.7.1.6)
petsc: 3.24.6 at /scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda
petsc configure options:  PETSC_ARCH=arch-linux-cuda-opt --prefix=/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-debugging=0 COPTFLAGS=\"-O3 -mcpu=neoverse-v2\" CXXOPTFLAGS=\"-O3 -mcpu=neoverse-v2\" FOPTFLAGS=\"-O3 -mcpu=neoverse-v2\" --with-shared-libraries=1 --with-64-bit-indices=0 --with-x=0 --with-blaslapack-lib=\"-L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core\" --download-scalapack --download-mumps --download-metis --download-parmetis --download-ptscotch --download-hypre --with-packages-download-dir=/scratch/10024/zhaochun/seas-gpu-deps/downloads --with-cuda=1 --with-cuda-arch=100
mumps: (no version line)
hypre: (no version line)
metis: -DCMAKE_INSTALL_PREFIX=/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda -DCMAKE_INSTALL_NAME_DIR:STRING="/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/15.3.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/15.3.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DGKLIB_PATH=../GKlib -DGKRAND=1 -DSHARED=1 -DMATH_LIB="-lm"
parmetis: -DCMAKE_INSTALL_PREFIX=/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda -DCMAKE_INSTALL_NAME_DIR:STRING="/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/15.3.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/15.3.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DGKLIB_PATH=../headers -DMETIS_PATH=/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda -DSHARED=1
ptscotch: (no version line)
scalapack: -DCMAKE_INSTALL_PREFIX=/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda -DCMAKE_INSTALL_NAME_DIR:STRING="/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/15.3.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/15.3.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DLAPACK_LIBRARIES="-Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.3/cuda/13.1/lib64 -L/home1/apps/nvidia/Linux_aarch64/26.3/cuda/13.1/lib64 -lcudart -lnvtx3interop -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.3/math_libs/13.1/lib64 -L/home1/apps/nvidia/Linux_aarch64/26.3/math_libs/13.1/lib64 -lcufft -lcublas -lcusparse -lcusolver -lcurand -L/home1/apps/nvidia/Linux_aarch64/26.3/cuda/13.1/lib64/stubs -lcuda -lm -Wl,-rpath,/opt/apps/gcc15/cuda13/openmpi/5.0.11/lib -Wl,-rpath,/opt/apps/gcc15/cuda13/openmpi/5.0.11/lib -L/opt/apps/gcc15/cuda13/openmpi/5.0.11/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -Wl,-rpath,/opt/apps/gcc/15.3.0/lib/gcc/aarch64-unknown-linux-gnu/15.3.0 -Wl,-rpath,/opt/apps/gcc/15.3.0/lib/gcc/aarch64-unknown-linux-gnu/15.3.0 -L/opt/apps/gcc/15.3.0/lib/gcc/aarch64-unknown-linux-gnu/15.3.0 -Wl,-rpath,/opt/apps/gcc/15.3.0/lib64 -Wl,-rpath,/opt/apps/gcc/15.3.0/lib64 -L/opt/apps/gcc/15.3.0/lib64 -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.3/cuda/13.1/lib64/stubs -L/home1/apps/nvidia/Linux_aarch64/26.3/cuda/13.1/lib64/stubs -Wl,-rpath,/opt/apps/gcc/15.3.0/aarch64-unknown-linux-gnu/lib -Wl,-rpath,/opt/apps/gcc/15.3.0/aarch64-unknown-linux-gnu/lib -L/opt/apps/gcc/15.3.0/aarch64-unknown-linux-gnu/lib -Wl,-rpath,/opt/apps/gcc/15.3.0/lib -Wl,-rpath,/opt/apps/gcc/15.3.0/lib -L/opt/apps/gcc/15.3.0/lib -Wl,-rpath,/opt/apps/gcc15/cuda13/openmpi/5.0.11/lib -lgfortran -lm -lgcc_s" -DSCALAPACK_BUILD_TESTS=OFF
#define HYPRE_RELEASE_VERSION "3.0.0"
#define HYPRE_USING_CUDA 1
#define HYPRE_USING_CUDA_STREAMS 1
#define HYPRE_USING_GPU 1
# build_horizon.sh MODE=cuda stage 'deps' on c101-036.horizon.tacc.utexas.edu at 2026-10-24T06:05:48Z
source: d079573ae1fc6f0fec349af6a94db35fb08d18bc (gpu-qd-bp5), dirty files: 9
modules: ucx/1.22.0 xalt/3.2.4 TACC gcc/15.3.0 cuda/13.1 ucc/1.9.0 openmpi/5.0.11 nvpl/26.9 cmake/4.4.0 
mpicc: gcc (GCC) 15.3.0
mpif90: GNU Fortran (GCC) 15.3.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
nvcc: Cuda compilation tools, release 13.1, V13.1.115
gpu arch: sm_100, mfem -arch=compute_90 -code=sm_100
cudss: /scratch/10024/zhaochun/seas-gpu-deps/cudss-0.7.1.6 (CUDSS_VERSION 0.7.1.6)
petsc: 3.24.6 at /scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda
petsc configure options:  PETSC_ARCH=arch-linux-cuda-opt --prefix=/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-debugging=0 COPTFLAGS=\"-O3 -mcpu=neoverse-v2\" CXXOPTFLAGS=\"-O3 -mcpu=neoverse-v2\" FOPTFLAGS=\"-O3 -mcpu=neoverse-v2\" --with-shared-libraries=1 --with-64-bit-indices=0 --with-x=0 --with-blaslapack-lib=\"-L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core\" --download-scalapack --download-mumps --download-metis --download-parmetis --download-ptscotch --download-hypre --with-packages-download-dir=/scratch/10024/zhaochun/seas-gpu-deps/downloads --with-cuda=1 --with-cuda-arch=100
mumps: (no version line)
hypre: (no version line)
metis: -DCMAKE_INSTALL_PREFIX=/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda -DCMAKE_INSTALL_NAME_DIR:STRING="/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/15.3.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/15.3.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DGKLIB_PATH=../GKlib -DGKRAND=1 -DSHARED=1 -DMATH_LIB="-lm"
parmetis: -DCMAKE_INSTALL_PREFIX=/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda -DCMAKE_INSTALL_NAME_DIR:STRING="/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/15.3.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/15.3.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DGKLIB_PATH=../headers -DMETIS_PATH=/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda -DSHARED=1
ptscotch: (no version line)
scalapack: -DCMAKE_INSTALL_PREFIX=/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda -DCMAKE_INSTALL_NAME_DIR:STRING="/scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/15.3.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/15.3.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc15/cuda13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DLAPACK_LIBRARIES="-Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.3/cuda/13.1/lib64 -L/home1/apps/nvidia/Linux_aarch64/26.3/cuda/13.1/lib64 -lcudart -lnvtx3interop -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.3/math_libs/13.1/lib64 -L/home1/apps/nvidia/Linux_aarch64/26.3/math_libs/13.1/lib64 -lcufft -lcublas -lcusparse -lcusolver -lcurand -L/home1/apps/nvidia/Linux_aarch64/26.3/cuda/13.1/lib64/stubs -lcuda -lm -Wl,-rpath,/opt/apps/gcc15/cuda13/openmpi/5.0.11/lib -Wl,-rpath,/opt/apps/gcc15/cuda13/openmpi/5.0.11/lib -L/opt/apps/gcc15/cuda13/openmpi/5.0.11/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -Wl,-rpath,/opt/apps/gcc/15.3.0/lib/gcc/aarch64-unknown-linux-gnu/15.3.0 -Wl,-rpath,/opt/apps/gcc/15.3.0/lib/gcc/aarch64-unknown-linux-gnu/15.3.0 -L/opt/apps/gcc/15.3.0/lib/gcc/aarch64-unknown-linux-gnu/15.3.0 -Wl,-rpath,/opt/apps/gcc/15.3.0/lib64 -Wl,-rpath,/opt/apps/gcc/15.3.0/lib64 -L/opt/apps/gcc/15.3.0/lib64 -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.3/cuda/13.1/lib64/stubs -L/home1/apps/nvidia/Linux_aarch64/26.3/cuda/13.1/lib64/stubs -Wl,-rpath,/opt/apps/gcc/15.3.0/aarch64-unknown-linux-gnu/lib -Wl,-rpath,/opt/apps/gcc/15.3.0/aarch64-unknown-linux-gnu/lib -L/opt/apps/gcc/15.3.0/aarch64-unknown-linux-gnu/lib -Wl,-rpath,/opt/apps/gcc/15.3.0/lib -Wl,-rpath,/opt/apps/gcc/15.3.0/lib -L/opt/apps/gcc/15.3.0/lib -Wl,-rpath,/opt/apps/gcc15/cuda13/openmpi/5.0.11/lib -lgfortran -lm -lgcc_s" -DSCALAPACK_BUILD_TESTS=OFF
#define HYPRE_RELEASE_VERSION "3.0.0"
#define HYPRE_USING_CUDA 1
#define HYPRE_USING_CUDA_STREAMS 1
#define HYPRE_USING_GPU 1
# build_horizon.sh MODE=cuda stage 'mfem' on c101-036.horizon.tacc.utexas.edu at 2026-10-24T06:09:11Z
source: d079573ae1fc6f0fec349af6a94db35fb08d18bc (gpu-qd-bp5), dirty files: 9
modules: ucx/1.22.0 xalt/3.2.4 TACC gcc/15.3.0 cuda/13.1 ucc/1.9.0 openmpi/5.0.11 nvpl/26.9 cmake/4.4.0 
mpicc: gcc (GCC) 15.3.0
mpif90: GNU Fortran (GCC) 15.3.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
nvcc: Cuda compilation tools, release 13.1, V13.1.115
gpu arch: sm_100, mfem -arch=compute_90 -code=sm_100
cudss: /scratch/10024/zhaochun/seas-gpu-deps/cudss-0.7.1.6 (CUDSS_VERSION 0.7.1.6)
petsc: 3.24.6 at /scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda
MFEM_USE_CUDSS         = YES
MFEM_USE_MPI           = YES
MFEM_USE_METIS         = YES
MFEM_USE_LAPACK        = YES
MFEM_USE_MUMPS         = YES
MFEM_USE_HDF5          = NO
MFEM_USE_PETSC         = YES
MFEM_USE_CUDA          = YES
# build_horizon.sh MODE=cuda stage 'seas' on c101-036.horizon.tacc.utexas.edu at 2026-10-24T06:10:00Z
source: d079573ae1fc6f0fec349af6a94db35fb08d18bc (gpu-qd-bp5), dirty files: 9
modules: ucx/1.22.0 xalt/3.2.4 TACC gcc/15.3.0 cuda/13.1 ucc/1.9.0 openmpi/5.0.11 nvpl/26.9 cmake/4.4.0 
mpicc: gcc (GCC) 15.3.0
mpif90: GNU Fortran (GCC) 15.3.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
nvcc: Cuda compilation tools, release 13.1, V13.1.115
gpu arch: sm_100, mfem -arch=compute_90 -code=sm_100
cudss: /scratch/10024/zhaochun/seas-gpu-deps/cudss-0.7.1.6 (CUDSS_VERSION 0.7.1.6)
petsc: 3.24.6 at /scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda
-rwx------ 1 zhaochun G-821164 57590200 Oct 24 01:09 seas_bp5_full
-rwx------ 1 zhaochun G-821164 54668744 Oct 24 01:09 seas_test_elasticity_operator
-rwx------ 1 zhaochun G-821164 56948224 Oct 24 01:09 seas_test_g0_matrix_only
# build_horizon.sh MODE=cuda stage 'seas' on c101-020.horizon.tacc.utexas.edu at 2026-10-24T20:41:00Z
source: 3f4053845fdf87ad833326fffecb339f3405fb96 (gpu-qd-bp5), dirty files: 25
modules: ucx/1.22.0 xalt/3.2.4 TACC gcc/15.3.0 cuda/13.1 ucc/1.9.0 openmpi/5.0.11 nvpl/26.9 cmake/4.4.0 
mpicc: gcc (GCC) 15.3.0
mpif90: GNU Fortran (GCC) 15.3.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
nvcc: Cuda compilation tools, release 13.1, V13.1.115
gpu arch: sm_100, mfem -arch=compute_90 -code=sm_100
cudss: /scratch/10024/zhaochun/seas-gpu-deps/cudss-0.7.1.6 (CUDSS_VERSION 0.7.1.6)
petsc: 3.24.6 at /scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda
-rwx------ 1 zhaochun G-821164 57590200 Oct 24 15:41 seas_bp5_full
-rwx------ 1 zhaochun G-821164 54668744 Oct 24 15:40 seas_test_elasticity_operator
-rwx------ 1 zhaochun G-821164 56948224 Oct 24 15:40 seas_test_g0_matrix_only
# build_horizon.sh MODE=cuda stage 'seas' on c101-020.horizon.tacc.utexas.edu at 2026-10-24T21:07:41Z
source: 3f4053845fdf87ad833326fffecb339f3405fb96 (gpu-qd-bp5), dirty files: 26
modules: ucx/1.22.0 xalt/3.2.4 TACC gcc/15.3.0 cuda/13.1 ucc/1.9.0 openmpi/5.0.11 nvpl/26.9 cmake/4.4.0 
mpicc: gcc (GCC) 15.3.0
mpif90: GNU Fortran (GCC) 15.3.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
nvcc: Cuda compilation tools, release 13.1, V13.1.115
gpu arch: sm_100, mfem -arch=compute_90 -code=sm_100
cudss: /scratch/10024/zhaochun/seas-gpu-deps/cudss-0.7.1.6 (CUDSS_VERSION 0.7.1.6)
petsc: 3.24.6 at /scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda
-rwx------ 1 zhaochun G-821164 57590200 Oct 24 16:07 seas_bp5_full
-rwx------ 1 zhaochun G-821164 54668744 Oct 24 16:07 seas_test_elasticity_operator
-rwx------ 1 zhaochun G-821164 56948224 Oct 24 16:07 seas_test_g0_matrix_only
# build_horizon.sh MODE=cuda stage 'seas' on c101-020.horizon.tacc.utexas.edu at 2026-10-24T21:14:13Z
source: 3f4053845fdf87ad833326fffecb339f3405fb96 (gpu-qd-bp5), dirty files: 26
modules: ucx/1.22.0 xalt/3.2.4 TACC gcc/15.3.0 cuda/13.1 ucc/1.9.0 openmpi/5.0.11 nvpl/26.9 cmake/4.4.0 
mpicc: gcc (GCC) 15.3.0
mpif90: GNU Fortran (GCC) 15.3.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
nvcc: Cuda compilation tools, release 13.1, V13.1.115
gpu arch: sm_100, mfem -arch=compute_90 -code=sm_100
cudss: /scratch/10024/zhaochun/seas-gpu-deps/cudss-0.7.1.6 (CUDSS_VERSION 0.7.1.6)
petsc: 3.24.6 at /scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda
-rwx------ 1 zhaochun G-821164 57590200 Oct 24 16:14 seas_bp5_full
-rwx------ 1 zhaochun G-821164 54668744 Oct 24 16:14 seas_test_elasticity_operator
-rwx------ 1 zhaochun G-821164 56948224 Oct 24 16:14 seas_test_g0_matrix_only
# build_horizon.sh MODE=cuda stage 'seas' on c105-035.horizon.tacc.utexas.edu at 2026-10-25T01:36:11Z
source: 3f4053845fdf87ad833326fffecb339f3405fb96 (gpu-qd-bp5), dirty files: 29
modules: ucx/1.22.0 xalt/3.2.4 TACC gcc/15.3.0 cuda/13.1 ucc/1.9.0 openmpi/5.0.11 nvpl/26.9 cmake/4.4.0 
mpicc: gcc (GCC) 15.3.0
mpif90: GNU Fortran (GCC) 15.3.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
nvcc: Cuda compilation tools, release 13.1, V13.1.115
gpu arch: sm_100, mfem -arch=compute_90 -code=sm_100
cudss: /scratch/10024/zhaochun/seas-gpu-deps/cudss-0.7.1.6 (CUDSS_VERSION 0.7.1.6)
petsc: 3.24.6 at /scratch/10024/zhaochun/seas-gpu-deps/petsc-3.24.6-cuda
-rwx------ 1 zhaochun G-821164 57590200 Oct 24 20:36 seas_bp5_full
-rwx------ 1 zhaochun G-821164 54668744 Oct 24 20:36 seas_test_elasticity_operator
-rwx------ 1 zhaochun G-821164 56948224 Oct 24 20:35 seas_test_g0_matrix_only
```

## Reproduce

On Horizon, from a clone of branch `gpu-qd-bp5` (the build job builds only a clean checkout at its upstream
commit, so not a detached checkout of 3f40538). Between 3f40538 and the commit that adds this folder, the only
code change is the `NSYS` line of `jobs/horizon/h_env.sh` (74aef91c), which the cycle does not use.

```
git clone -b gpu-qd-bp5 https://github.com/chunhuizhao478/seas-mfem.git seas-mfem-gpu-cuda
cd seas-mfem-gpu-cuda/miniapps/seas
source jobs/horizon/h_env.sh && MODE=cuda STAGE=fetch bash ../../build_horizon.sh   # login node (internet)
QD_PATH=green SUBMIT=1 H_CHAIN="build green cycle" bash jobs/horizon/h_submit.sh   # build, G file, this run
```

## Speed

From `speed.txt` and `stage_timers_final.txt`: 1500 yr, 7 events, 313,950 accepted steps (1,968,895 RK
stages) in 10,729 s = 2.98 h, 0.0342 s/step with setup (0.0333 s/step in the time loop). The CPU gold
1074099 (Vista, 144 ranks, MUMPS) took 0.36 s/step, 31.6 h: 10.6x per step. The G product is 95.8 % of the
loop (5.09 ms per stage); friction 0.9 %, slip expansion 1.0 %, traction 0.2 %, output 0.35 %.

## Corrections

The job's README (sha256 7a808deb5c45c4304a78480042722000aa388ecd199043f326ecf0b5698127be) was changed in three places:

- Cluster: it said Vista; the job ran on Horizon.
- Inputs: it listed only `G0_V63_FLAGS`, which still hold `--solver mumps`. It now gives the whole run line.
- Reproduce: it gave the Vista build (`GPU_ARCH=90`, `jobs/vista/g1_cuda_build.sbatch`). It now gives the
  Horizon chain.

The Speed section and the G line under Inputs were added. MANIFEST.sha256 holds the new README's sha256 and
the sha256 of `compare_tandem_p1/README.md`; its other lines are as they came from Horizon.
