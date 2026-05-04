# SAFS C++ tools — CGAL-based corefinement

This directory contains C++ tools that operate on the SAFS per-fault
STLs in support of `PLAN_cgal_corefine.md`. Phase 0 lands the build
infrastructure and a hello-world smoke binary
(`corefine_faults_hello`); Phase 1 will replace it with the real
`corefine_faults` tool.

## 1. Prerequisites

| Tool / library | Min version | Tested with | Notes |
|---|---|---|---|
| CGAL | 5.4 | **5.6.1** (pinned) | Header-only since 5.0; `protect_constraints(true)` semantics for `isotropic_remeshing` documented as exact preservation since 5.4. CGAL 6.0 introduces ABI changes and is **out of bounds** — enforced by `find_package(CGAL 5.4)` + version-guard in `CMakeLists.txt`. |
| GMP | 6.0 | 6.3.0 | Required by CGAL kernels. |
| MPFR | 4.0 | 4.2.1 | Required by CGAL kernels. |
| Boost | 1.66 | 1.84.0 | Headers only — CGAL uses Boost.Graph and Boost.Property_map. |
| CMake | 3.18 | 4.2.1 | Modern CMake `find_package(CGAL)` integration. |
| C++ toolchain | C++17 | clang (macOS), gcc 9.1 (Frontera) | Need `std::string_view`, `std::filesystem`. |

## 2. Laptop recipe (macOS / Linux, conda-forge)

```bash
# 1. Activate the same env that builds the rest of the SEAS C++ code.
conda activate mfem-dev

# 2. Install CGAL + Boost + GMP/MPFR.  The exact pin matches what was
#    tested against the corefinement_difference_remeshed CGAL example.
conda install -c conda-forge cgal=5.6.1 boost-cpp gmp mpfr

# 3. Configure and build.
cmake -S miniapps/seas/safs/tools \
      -B miniapps/seas/safs/tools/build \
      -DCMAKE_BUILD_TYPE=Release
cmake --build miniapps/seas/safs/tools/build

# 4. Smoke test.
./miniapps/seas/safs/tools/build/corefine_faults_hello
# Expected:  CGAL version: 5.6.1
```

### What the binary links against

CGAL 5.x is header-only; the produced binary therefore does NOT link
against a `libcgal.dylib` on macOS or `libcgal.so` on Linux. Boost on
conda-forge ships a header-only-usable distribution, so likewise no
`libboost_*` line appears.

GMP / MPFR appear in `otool -L` output **only** for binaries that
actually instantiate kernels using exact arithmetic (e.g. `EPECK`).
The Phase-0 hello binary uses `EPICK` (double-precision predicates,
inexact constructions) so its only dynamic dependencies are the C++
runtime and `libSystem`:

```text
$ otool -L tools/build/corefine_faults_hello
tools/build/corefine_faults_hello:
        @rpath/libc++.1.dylib (compatibility version 1.0.0, current version 1.0.0)
        /usr/lib/libSystem.B.dylib (...)
```

Phase 1's `corefine_faults` will pull `libgmp.10.dylib` and
`libmpfr.6.dylib` once it instantiates `corefine` code paths that
internally use exact predicates. Plan §Numerical constraints
explicitly avoids EPECK for Phases 1–3, so the dependency footprint
stays small.

### CGAL 5.6.1 + macOS Apple-Clang workaround

The conda-forge `cgal-cpp 5.6.1` package has a known compile error
on Apple Clang in `<CGAL/boost/graph/iterator.h>` — three of the
fault-mesh iterator classes' `operator bool()` overloads call
`this->base()` but don't actually inherit a `base()` method.

The build sidesteps this with an in-tree header shim at
`tools/cgal_patch/CGAL/boost/graph/iterator.h`. The shim replaces
the broken `this->base()` calls with the semantically-equivalent
`this->g != nullptr` (default-constructed iterators set `g=nullptr`).
`tools/CMakeLists.txt` puts the shim FIRST on the include path so it
shadows the conda env's CGAL header for our targets only — the conda
env stays untouched. See `tools/cgal_patch/README.md` for full
provenance and removal conditions (delete the shim if/when conda-
forge ships a fixed package, or the project upgrades to CGAL ≥ 6.0).

## 3. Frontera recipe — DOCUMENTED, NOT YET EXECUTED

> **Important:** Per project memory `feedback_frontera_approval.md`,
> Frontera commands (sbatch, ssh, allocation checks) require explicit
> user approval before being run. The recipe below is documentation
> only; do not execute without authorisation.

### a. Module list

Copy verbatim from a known-working seas-mfem sbatch — per project
memory `feedback_sbatch_modules.md`, do NOT trim modules even if
they look unrelated (`fftw3` in particular has caused link
failures in past Frontera builds):

```bash
module reset
module load intel/19.1.1 impi/19.0.9 boost/1.71 \
            gcc/9.1.0 cmake/3.24.2 fftw3/3.3.10 \
            gmp mpfr
export LD_LIBRARY_PATH=$TACC_GMP_LIB:$TACC_MPFR_LIB:$LD_LIBRARY_PATH
```

Verify the exact module versions by reading the live `seas_tpv102`
sbatch file the next time the user authorises a login session.

### b. CGAL fetch (header-only since 5.0 → no compile step)

```bash
cd $WORK
git clone --depth 1 -b v5.6.1 https://github.com/CGAL/cgal.git cgal-5.6.1
export CGAL_DIR=$WORK/cgal-5.6.1/lib/cmake/CGAL
```

### c. CMake build

```bash
cd $WORK/seas-mfem/miniapps/seas/safs/tools
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH="$CGAL_DIR;$TACC_BOOST_DIR"
cmake --build build
./build/corefine_faults_hello
```

### d. Acceptance for Frontera Phase 0

The `corefine_faults_hello` smoke binary builds and prints the CGAL
version on a Frontera **login node** before any compute-node work
begins.

### e. Failure mode — what to do if Frontera build fails

If the build dies with missing fftw3 / `LD_LIBRARY_PATH` / Boost
errors: copy the module list **VERBATIM** from the live
`seas_tpv102` sbatch (per `feedback_sbatch_modules.md`); do NOT
trim modules even if they look unrelated to the CGAL build. The
seas_tpv102 module list has been validated end-to-end and is the
canonical reference.

## 4. Smoke build target

```bash
cmake --build miniapps/seas/safs/tools/build --target corefine_faults_hello
```

Successful output ends with `Linking CXX executable corefine_faults_hello`
and the binary appears in `tools/build/`.

## 5. CMake version-guard rationale

The CMakeLists requires CGAL ≥ 5.4 and refuses CGAL ≥ 6.0:

- **5.4 minimum** — `protect_constraints(true)` on
  `CGAL::Polygon_mesh_processing::isotropic_remeshing` is documented
  as exact preservation since this version. The tool relies on this
  property in Phase 2.
- **5.6.1 tested** — matches the
  `corefinement_difference_remeshed` example referenced in
  `mesh/REVIEW.md` R-010, which is the basis for the Phase 1
  implementation.
- **6.0 forbidden** — major-version bump may have changed ABI or
  the `protect_constraints` semantics. A future ticket would need
  to validate against CGAL 6.x before lifting this guard.

Tested window: **`[5.4, 6.0)`**. Pinned install: **5.6.1**.
