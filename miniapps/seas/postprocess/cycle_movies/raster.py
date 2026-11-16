# Pixel grid <- triangles: for every pixel centre the triangle that holds it
# and the three barycentric weights, so that a field given at the triangle
# corners (the per-element points of a volume snapshot) is interpolated
# linearly inside each triangle.  Computed once per movie and reused for
# every frame: a frame is then one gather and one weighted sum.
import numpy as np


def locate(points, tris, x0, dx, nx, y0, dy, ny, pairs_per_chunk=2000000, tol=1e-9):
    """points (n, 2), tris (m, 3) node ids; pixel (i, j) has its centre at
    (x0 + (i + 1/2) dx, y0 + (j + 1/2) dy), i < nx, j < ny.

    Returns tri (ny * nx,) int, -1 where no triangle holds the centre, and
    w (ny * nx, 3), the barycentric weights (rows of pixels without a
    triangle are 0).  Row-major: index j * nx + i.

    Exact: every triangle is tested against every pixel centre inside its
    bounding box, so a triangle smaller than a pixel is found whenever a
    centre lies in it and no pixel is given a neighbour's values.
    """
    a, b, c = points[tris[:, 0]], points[tris[:, 1]], points[tris[:, 2]]
    lo = np.minimum(np.minimum(a, b), c)
    hi = np.maximum(np.maximum(a, b), c)
    i0 = np.clip(np.ceil((lo[:, 0] - x0) / dx - 0.5).astype(np.int64), 0, nx)
    i1 = np.clip(np.floor((hi[:, 0] - x0) / dx - 0.5).astype(np.int64), -1, nx - 1)
    j0 = np.clip(np.ceil((lo[:, 1] - y0) / dy - 0.5).astype(np.int64), 0, ny)
    j1 = np.clip(np.floor((hi[:, 1] - y0) / dy - 0.5).astype(np.int64), -1, ny - 1)
    wi = np.maximum(i1 - i0 + 1, 0)
    hj = np.maximum(j1 - j0 + 1, 0)
    npair = wi * hj
    cand = np.where(npair > 0)[0]
    tri = np.full(nx * ny, -1, dtype=np.int64)
    w = np.zeros((nx * ny, 3))
    best = np.full(nx * ny, -np.inf)
    cum = np.cumsum(npair[cand])
    start = 0
    while start < len(cand):
        base = cum[start - 1] if start else 0
        stop = int(np.searchsorted(cum, base + pairs_per_chunk, side="right"))
        stop = max(stop, start + 1)
        t = cand[start:stop]
        n = npair[t]
        rep = np.repeat(np.arange(len(t)), n)
        local = np.arange(int(n.sum())) - np.repeat(np.cumsum(n) - n, n)
        tt = t[rep]
        ii = i0[tt] + local % wi[tt]
        jj = j0[tt] + local // wi[tt]
        q = np.stack([x0 + (ii + 0.5) * dx, y0 + (jj + 0.5) * dy], axis=1)
        v0, v1, v2 = b[tt] - a[tt], c[tt] - a[tt], q - a[tt]
        den = v0[:, 0] * v1[:, 1] - v1[:, 0] * v0[:, 1]
        ok = den != 0
        den = np.where(ok, den, 1.0)
        l1 = (v2[:, 0] * v1[:, 1] - v1[:, 0] * v2[:, 1]) / den
        l2 = (v0[:, 0] * v2[:, 1] - v2[:, 0] * v0[:, 1]) / den
        bw = np.stack([1.0 - l1 - l2, l1, l2], axis=1)
        score = np.where(ok, bw.min(axis=1), -np.inf)
        pix = jj * nx + ii
        keep = score >= -tol
        # a centre on an edge is in two triangles: keep the one it is deepest in
        order = np.argsort(score[keep], kind="stable")
        pk, sk, tk, bk = pix[keep][order], score[keep][order], tt[keep][order], bw[keep][order]
        better = sk > best[pk]
        pk, sk, tk, bk = pk[better], sk[better], tk[better], bk[better]
        tri[pk], w[pk], best[pk] = tk, bk, sk  # later (deeper) entries win
        start = stop
    w = np.clip(w, 0.0, None)
    s = w.sum(axis=1, keepdims=True)
    w = np.where(s > 0, w / np.where(s > 0, s, 1.0), 0.0)
    return tri, w
