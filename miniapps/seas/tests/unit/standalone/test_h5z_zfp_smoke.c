/* Phase 2d.1 standalone smoke test (PLAN_paraview_compaction_2026-04-28.md §2d.1.2).
 *
 * Vanilla HDF5 + H5Z-ZFP round-trip with no MFEM dependency.  Verifies:
 *   1. The H5Z-ZFP plugin loads from HDF5_PLUGIN_PATH (or compiled-in path).
 *   2. A 1M-element double dataset compresses with filter id 32013 in
 *      ACCURACY mode at tol=1e-3.
 *   3. After decompression max abs error <= 1e-3.
 *   4. Compression ratio > 4x vs the uncompressed equivalent (8 MB -> < 2 MB).
 *
 * Build:
 *   h5cc -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib \
 *        -o test_h5z_zfp_smoke test_h5z_zfp_smoke.c
 * Run:
 *   HDF5_PLUGIN_PATH=$CONDA_PREFIX/plugin ./test_h5z_zfp_smoke
 * Inspect:
 *   h5dump -p -A -d /vals /tmp/zfp_smoke.h5 | grep -i 'FILTERS\|ZFP\|32013'
 */
#include <hdf5.h>
#include "H5Zzfp_plugin.h"     /* H5Pset_zfp_accuracy_cdata macro */
#include "H5Zzfp_version.h"    /* H5Z_FILTER_ZFP, H5Z_ZFP_MODE_ACCURACY */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>

int main(void)
{
    const hsize_t N = 1024 * 1024;        /* 1M doubles = 8 MB raw */
    const double  TOL = 1e-3;             /* ZFP accuracy mode tolerance */
    const char   *PATH = "/tmp/zfp_smoke.h5";

    /* ----- Generate input data with smooth-but-nontrivial structure. */
    double *buf = (double *) malloc(N * sizeof(double));
    if (!buf) { fprintf(stderr, "malloc failed\n"); return 1; }
    for (hsize_t i = 0; i < N; ++i)
    {
        const double x = (double)i / (double)N;
        buf[i] = sin(50.0 * x) * exp(-3.0 * x) + 0.1 * x;
    }

    /* ----- Create file + dataset with chunking + ZFP filter. */
    hid_t f = H5Fcreate(PATH, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (f < 0) { fprintf(stderr, "FAIL: H5Fcreate\n"); return 1; }
    hid_t s = H5Screate_simple(1, &N, NULL);
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    const hsize_t chunk = 65536;
    H5Pset_chunk(dcpl, 1, &chunk);

    /* H5Z-ZFP cd_values via the macro (cd_nelmts == 4 in ACCURACY mode). */
    unsigned int cd_values[10] = {0};
    size_t cd_nelmts = 10;
    H5Pset_zfp_accuracy_cdata(TOL, cd_nelmts, cd_values);
    if (H5Pset_filter(dcpl, H5Z_FILTER_ZFP, H5Z_FLAG_MANDATORY,
                      cd_nelmts, cd_values) < 0)
    {
        fprintf(stderr, "FAIL: H5Pset_filter (filter %u not registered? "
                "set HDF5_PLUGIN_PATH)\n", H5Z_FILTER_ZFP);
        return 1;
    }

    hid_t d = H5Dcreate2(f, "vals", H5T_IEEE_F64LE, s,
                          H5P_DEFAULT, dcpl, H5P_DEFAULT);
    if (d < 0) { fprintf(stderr, "FAIL: H5Dcreate2\n"); return 1; }
    if (H5Dwrite(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                 H5P_DEFAULT, buf) < 0)
    {
        fprintf(stderr, "FAIL: H5Dwrite (filter chain failed at first chunk)\n");
        return 1;
    }
    H5Dclose(d);  H5Pclose(dcpl);  H5Sclose(s);  H5Fclose(f);

    /* ----- Reopen + read back; verify numerical accuracy. */
    f = H5Fopen(PATH, H5F_ACC_RDONLY, H5P_DEFAULT);
    d = H5Dopen2(f, "vals", H5P_DEFAULT);
    double *rb = (double *) malloc(N * sizeof(double));
    if (!rb) { fprintf(stderr, "malloc failed\n"); return 1; }
    H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rb);
    double max_err = 0.0;
    for (hsize_t i = 0; i < N; ++i)
    {
        double e = rb[i] - buf[i];
        if (e < 0) e = -e;
        if (e > max_err) max_err = e;
    }
    H5Dclose(d);  H5Fclose(f);

    /* ----- Verify compression ratio. */
    struct stat st;
    if (stat(PATH, &st) != 0)
    { fprintf(stderr, "FAIL: stat(%s)\n", PATH); return 1; }
    const double raw_bytes = (double)N * 8.0;
    const double file_bytes = (double)st.st_size;
    const double ratio = raw_bytes / file_bytes;

    /* ----- Report. */
    fprintf(stderr, "max abs err: %.3e (tol %.3e)\n", max_err, TOL);
    fprintf(stderr, "file size:   %.3f MB (raw %.3f MB)\n",
            file_bytes / 1.0e6, raw_bytes / 1.0e6);
    fprintf(stderr, "ratio:       %.2fx\n", ratio);

    int fail = 0;
    if (max_err > TOL)
    { fprintf(stderr, "FAIL: max abs err %.3e > tol %.3e\n", max_err, TOL); fail = 1; }
    if (ratio < 4.0)
    { fprintf(stderr, "FAIL: compression ratio %.2fx < 4x\n", ratio); fail = 1; }
    if (!fail) { fprintf(stderr, "PASS\n"); }

    free(buf); free(rb);
    return fail;
}
