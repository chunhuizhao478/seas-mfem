# p1_1500yr_r144_job1074099

Produced from commit `c1fbc7bc069806afdbee053f32a4cad6c4ef101d` (branch `gpu-qd-bp5`) of https://github.com/chunhuizhao478/seas-mfem:
https://github.com/chunhuizhao478/seas-mfem/commit/c1fbc7bc069806afdbee053f32a4cad6c4ef101d

The checkout had no uncommitted changes in tracked files, so these results come from that commit exactly.

This README was written after the run, on 2026-10-18, from the job's provenance.txt, sacct and the commit; the run itself recorded only provenance.txt.

## Job

- Cluster: Vista (i615-091.vista.tacc.utexas.edu), Slurm job 1074099 (bp5_cycle), started 2026-10-16T12:54:19Z
- Partition gg, account ear20006, 1 node(s), 144 task(s)
- Job script: `miniapps/seas/jobs/vista/bp5_full_cycle.sbatch` at that commit
- Submit line: `sbatch jobs/vista/bp5_full_cycle.sbatch`
- Settings from the environment: (not recorded by this job; see the submit line)

## Inputs

- Flags shared by the G-phase runs (`G0_V63_FLAGS` in `miniapps/seas/jobs/vista/g0_env.sh`; the job script adds the mesh, order and per-run flags):

```
--mesh-scale 1000 --dg-method IP --tandem-time-stepping --solver mumps --petsc-ts --petsc-ts-options tests/verification/petsc_ts_rk45_tandem.cfg
```

- Meshes at that commit (paths relative to `miniapps/seas`), sha256:

```
f413b957cf1c6d02d2c84b65fc05d1f99232301ed8083cb7d8533c9fa442c3b4  bp5/mesh/reference/bp5_tandem_exact.msh
34dd17748445d1982f867389023b3873e401380baf18a7a0c2c97591e3ff6121  bp5/mesh/reference/bp5_tandem_500m.msh
c2bffe326e3e7c8df0be73c97fed734225d8968b1f8eb776fbd9f7055bd23dbc  bp5/mesh/reference/bp5_tandem_250m.msh
```

## Build

- Binary: `miniapps/seas/seas_bp5_full`, sha256 8b9bf3077644e919 (first 16 hex digits, from provenance.txt)
- Built with `MODE=cpu bash build_horizon.sh` (PETSc and its packages, MFEM, then seas). What the build used:

```
# build_horizon.sh MODE=cpu stage 'fetch' on login2.vista.tacc.utexas.edu at 2026-10-15T01:37:32Z
source: ecc6036e6547de6e69ff658956732593eb2c8f5e (gpu-qd-bp5), dirty files: 1
modules: ucc/1.9.0 ucx/1.22.0 xalt/3.1 TACC gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5 
mpicc: gcc (GCC) 13.2.0
mpif90: GNU Fortran (GCC) 13.2.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
petsc: 3.24.6 at /work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu
# build_horizon.sh MODE=cpu stage 'deps' on i614-052.vista.tacc.utexas.edu at 2026-10-15T03:07:07Z
source: ecc6036e6547de6e69ff658956732593eb2c8f5e (gpu-qd-bp5), dirty files: 4
modules: ucc/1.9.0 ucx/1.22.0 xalt/3.1 TACC gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5 
mpicc: gcc (GCC) 13.2.0
mpif90: GNU Fortran (GCC) 13.2.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
petsc: 3.24.6 at /work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu
petsc configure options:  --prefix=/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-debugging=0 COPTFLAGS=\"-O3 -mcpu=neoverse-v2\" CXXOPTFLAGS=\"-O3 -mcpu=neoverse-v2\" FOPTFLAGS=\"-O3 -mcpu=neoverse-v2\" --with-shared-libraries=1 --with-64-bit-indices=0 --with-x=0 --with-cuda=0 --with-blaslapack-lib=\"-L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core\" --download-scalapack --download-mumps --download-metis --download-parmetis --download-ptscotch --download-hypre --with-packages-download-dir=/work/10024/zhaochun/vista/seas-gpu-deps/downloads
mumps: (no version line)
hypre: (no version line)
metis: -DCMAKE_INSTALL_PREFIX=/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu -DCMAKE_INSTALL_NAME_DIR:STRING="/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/13.2.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/13.2.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DGKLIB_PATH=../GKlib -DGKRAND=1 -DSHARED=1 -DMATH_LIB="-lm"
parmetis: -DCMAKE_INSTALL_PREFIX=/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu -DCMAKE_INSTALL_NAME_DIR:STRING="/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/13.2.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/13.2.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DGKLIB_PATH=../headers -DMETIS_PATH=/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu -DSHARED=1
ptscotch: (no version line)
scalapack: -DCMAKE_INSTALL_PREFIX=/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu -DCMAKE_INSTALL_NAME_DIR:STRING="/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/13.2.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/13.2.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DLAPACK_LIBRARIES="-Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core -lm -Wl,-rpath,/opt/apps/gcc13/openmpi/5.0.11/lib -Wl,-rpath,/opt/apps/gcc13/openmpi/5.0.11/lib -L/opt/apps/gcc13/openmpi/5.0.11/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -Wl,-rpath,/opt/apps/gcc/13.2.0/lib/gcc/aarch64-unknown-linux-gnu/13.2.0 -Wl,-rpath,/opt/apps/gcc/13.2.0/lib/gcc/aarch64-unknown-linux-gnu/13.2.0 -L/opt/apps/gcc/13.2.0/lib/gcc/aarch64-unknown-linux-gnu/13.2.0 -Wl,-rpath,/opt/apps/gcc/13.2.0/lib64 -Wl,-rpath,/opt/apps/gcc/13.2.0/lib64 -L/opt/apps/gcc/13.2.0/lib64 -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -Wl,-rpath,/opt/apps/gcc/13.2.0/aarch64-unknown-linux-gnu/lib -Wl,-rpath,/opt/apps/gcc/13.2.0/aarch64-unknown-linux-gnu/lib -L/opt/apps/gcc/13.2.0/aarch64-unknown-linux-gnu/lib -Wl,-rpath,/opt/apps/gcc/13.2.0/lib -Wl,-rpath,/opt/apps/gcc/13.2.0/lib -L/opt/apps/gcc/13.2.0/lib -Wl,-rpath,/opt/apps/gcc13/openmpi/5.0.11/lib -lgfortran -lm -lgcc_s" -DSCALAPACK_BUILD_TESTS=OFF
#define HYPRE_RELEASE_VERSION "3.0.0"
# build_horizon.sh MODE=cpu stage 'mfem' on i614-052.vista.tacc.utexas.edu at 2026-10-15T03:07:42Z
source: ecc6036e6547de6e69ff658956732593eb2c8f5e (gpu-qd-bp5), dirty files: 4
modules: ucc/1.9.0 ucx/1.22.0 xalt/3.1 TACC gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5 
mpicc: gcc (GCC) 13.2.0
mpif90: GNU Fortran (GCC) 13.2.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
petsc: 3.24.6 at /work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu
MFEM_USE_MPI           = YES
MFEM_USE_METIS         = YES
MFEM_USE_LAPACK        = YES
MFEM_USE_MUMPS         = YES
MFEM_USE_HDF5          = NO
MFEM_USE_PETSC         = YES
MFEM_USE_CUDA          = NO
# build_horizon.sh MODE=cpu stage 'seas' on i614-052.vista.tacc.utexas.edu at 2026-10-15T03:08:05Z
source: ecc6036e6547de6e69ff658956732593eb2c8f5e (gpu-qd-bp5), dirty files: 4
modules: ucc/1.9.0 ucx/1.22.0 xalt/3.1 TACC gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5 
mpicc: gcc (GCC) 13.2.0
mpif90: GNU Fortran (GCC) 13.2.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
petsc: 3.24.6 at /work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu
-rwx------ 1 zhaochun G-821164 15954992 Oct 14 22:08 seas_bp5_full
-rwx------ 1 zhaochun G-821164 13238608 Oct 14 22:08 seas_test_elasticity_operator
-rwx------ 1 zhaochun G-821164 15657488 Oct 14 22:07 seas_test_g0_matrix_only
# build_horizon.sh MODE=cpu stage 'deps' on i617-131.vista.tacc.utexas.edu at 2026-10-15T03:40:22Z
source: a4e6f90e1c4840d0094479be78ced02cf75214a2 (gpu-qd-bp5), dirty files: 8
modules: ucc/1.9.0 ucx/1.22.0 xalt/3.1 TACC gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5 
mpicc: gcc (GCC) 13.2.0
mpif90: GNU Fortran (GCC) 13.2.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
petsc: 3.24.6 at /work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu
petsc configure options:  --prefix=/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-debugging=0 COPTFLAGS=\"-O3 -mcpu=neoverse-v2\" CXXOPTFLAGS=\"-O3 -mcpu=neoverse-v2\" FOPTFLAGS=\"-O3 -mcpu=neoverse-v2\" --with-shared-libraries=1 --with-64-bit-indices=0 --with-x=0 --with-cuda=0 --with-blaslapack-lib=\"-L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core\" --download-scalapack --download-mumps --download-metis --download-parmetis --download-ptscotch --download-hypre --with-packages-download-dir=/work/10024/zhaochun/vista/seas-gpu-deps/downloads
mumps: (no version line)
hypre: (no version line)
metis: -DCMAKE_INSTALL_PREFIX=/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu -DCMAKE_INSTALL_NAME_DIR:STRING="/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/13.2.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/13.2.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DGKLIB_PATH=../GKlib -DGKRAND=1 -DSHARED=1 -DMATH_LIB="-lm"
parmetis: -DCMAKE_INSTALL_PREFIX=/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu -DCMAKE_INSTALL_NAME_DIR:STRING="/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/13.2.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/13.2.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DGKLIB_PATH=../headers -DMETIS_PATH=/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu -DSHARED=1
ptscotch: (no version line)
scalapack: -DCMAKE_INSTALL_PREFIX=/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu -DCMAKE_INSTALL_NAME_DIR:STRING="/work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu/lib" -DCMAKE_INSTALL_LIBDIR:STRING="lib" -DCMAKE_VERBOSE_MAKEFILE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_AR="/opt/apps/gcc/13.2.0/bin/ar" -DCMAKE_C_COMPILER="mpicc" -DMPI_C_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicc" -DCMAKE_RANLIB=/opt/apps/gcc/13.2.0/bin/ranlib -DCMAKE_C_FLAGS:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_C_FLAGS_DEBUG:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_C_FLAGS_RELEASE:STRING="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -O3 -mcpu=neoverse-v2 -Wno-implicit-int -Wno-int-conversion -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -fno-common" -DCMAKE_CXX_COMPILER="mpicxx" -DMPI_CXX_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpicxx" -DCMAKE_CXX_FLAGS:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_DEBUG:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_CXX_FLAGS_RELEASE:STRING="-Wno-lto-type-mismatch -Wno-psabi -O3 -mcpu=neoverse-v2 -fPIC" -DCMAKE_Fortran_COMPILER="mpif90" -DMPI_Fortran_COMPILER="/opt/apps/gcc13/openmpi/5.0.11/bin/mpif90" -DCMAKE_Fortran_FLAGS:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_DEBUG:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O3 -mcpu=neoverse-v2 -fallow-argument-mismatch" -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_STATIC_LIBS:BOOL=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DLAPACK_LIBRARIES="-Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core -lm -Wl,-rpath,/opt/apps/gcc13/openmpi/5.0.11/lib -Wl,-rpath,/opt/apps/gcc13/openmpi/5.0.11/lib -L/opt/apps/gcc13/openmpi/5.0.11/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -Wl,-rpath,/opt/apps/gcc/13.2.0/lib/gcc/aarch64-unknown-linux-gnu/13.2.0 -Wl,-rpath,/opt/apps/gcc/13.2.0/lib/gcc/aarch64-unknown-linux-gnu/13.2.0 -L/opt/apps/gcc/13.2.0/lib/gcc/aarch64-unknown-linux-gnu/13.2.0 -Wl,-rpath,/opt/apps/gcc/13.2.0/lib64 -Wl,-rpath,/opt/apps/gcc/13.2.0/lib64 -L/opt/apps/gcc/13.2.0/lib64 -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -Wl,-rpath,/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -Wl,-rpath,/opt/apps/gcc/13.2.0/aarch64-unknown-linux-gnu/lib -Wl,-rpath,/opt/apps/gcc/13.2.0/aarch64-unknown-linux-gnu/lib -L/opt/apps/gcc/13.2.0/aarch64-unknown-linux-gnu/lib -Wl,-rpath,/opt/apps/gcc/13.2.0/lib -Wl,-rpath,/opt/apps/gcc/13.2.0/lib -L/opt/apps/gcc/13.2.0/lib -Wl,-rpath,/opt/apps/gcc13/openmpi/5.0.11/lib -lgfortran -lm -lgcc_s" -DSCALAPACK_BUILD_TESTS=OFF
#define HYPRE_RELEASE_VERSION "3.0.0"
# build_horizon.sh MODE=cpu stage 'mfem' on i617-131.vista.tacc.utexas.edu at 2026-10-15T03:40:58Z
source: a4e6f90e1c4840d0094479be78ced02cf75214a2 (gpu-qd-bp5), dirty files: 8
modules: ucc/1.9.0 ucx/1.22.0 xalt/3.1 TACC gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5 
mpicc: gcc (GCC) 13.2.0
mpif90: GNU Fortran (GCC) 13.2.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
petsc: 3.24.6 at /work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu
MFEM_USE_MPI           = YES
MFEM_USE_METIS         = YES
MFEM_USE_LAPACK        = YES
MFEM_USE_MUMPS         = YES
MFEM_USE_HDF5          = NO
MFEM_USE_PETSC         = YES
MFEM_USE_CUDA          = NO
# build_horizon.sh MODE=cpu stage 'seas' on i617-131.vista.tacc.utexas.edu at 2026-10-15T03:41:21Z
source: a4e6f90e1c4840d0094479be78ced02cf75214a2 (gpu-qd-bp5), dirty files: 8
modules: ucc/1.9.0 ucx/1.22.0 xalt/3.1 TACC gcc/13.2.0 openmpi/5.0.11 nvpl/26.9 cmake/3.31.5 
mpicc: gcc (GCC) 13.2.0
mpif90: GNU Fortran (GCC) 13.2.0
mpi: Open MPI: 5.0.11
blas/lapack: -L/home1/apps/nvidia/Linux_aarch64/26.9/math_libs/nvpl/lib -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core
cpu flags: -mcpu=neoverse-v2
petsc: 3.24.6 at /work/10024/zhaochun/vista/seas-gpu-deps/petsc-3.24.6-cpu
-rwx------ 1 zhaochun G-821164 15954992 Oct 14 22:41 seas_bp5_full
-rwx------ 1 zhaochun G-821164 13238632 Oct 14 22:41 seas_test_elasticity_operator
-rwx------ 1 zhaochun G-821164 15657488 Oct 14 22:41 seas_test_g0_matrix_only
```

## Reproduce

```
git clone https://github.com/chunhuizhao478/seas-mfem.git && cd seas-mfem
git checkout c1fbc7bc069806afdbee053f32a4cad6c4ef101d
STAGE=fetch bash build_horizon.sh                          # login node
cd miniapps/seas && sbatch jobs/vista/g0_build.sbatch       # deps, MFEM, seas
sbatch jobs/vista/bp5_full_cycle.sbatch   # from miniapps/seas, after the build
```
