#!/usr/bin/env python3
"""
Compute peak ground-motion parameters (PGV, PGA, PGD) from SeisSol free-surface output.

Two free-surface sources are supported:

  receivers : SeisSol off-fault receiver files (safs-receiver-*.dat).
              Columns: Time, xx,yy,zz,xy,yz,xz, v1,v2,v3   (velocity in m/s).
              Fine time sampling (dt ~ 5 ms) -> accurate point PGA/PGV/PGD.
              Displacement is obtained by time-integrating velocity (cumulative trapezoid).

  surface   : SeisSol free-surface field (safs-surface.xdmf + safs-surface_cell/...).
              Per-cell v1,v2,v3 (m/s) AND u1,u2,u3 (m, displacement already integrated).
              Coarse time sampling (dt ~ 0.5 s) -> good for spatial PGV/PGD maps;
              PGA is unreliable (aliased) and is reported only as a rough estimate.

Conventions
  v1,u1 = x component, v2,u2 = y component, v3,u3 = z (vertical) component.
  "horizontal" magnitude = sqrt(c1^2 + c2^2);  "3d" = sqrt(c1^2+c2^2+c3^2);  "vertical" = |c3|.
  PGV [m/s], PGD [m], PGA reported in g (1 g = 9.80665 m/s^2) and in m/s^2.

Usage
  python compute_pgm.py receivers <output_dir> [--prefix safs] [--out pgm_receivers.csv]
                                  [--detrend-disp]
  python compute_pgm.py surface   <output_dir> [--prefix safs] [--out pgm_surface.vtu]
"""
import argparse, glob, os, re, sys
import numpy as np

G = 9.80665  # m/s^2


# ----------------------------------------------------------------------------- helpers
def peaks(c1, c2, c3):
    """Return dict of peak magnitudes (horizontal, vertical, 3d) of a 3-component series."""
    hor = np.sqrt(c1 * c1 + c2 * c2)
    mag = np.sqrt(c1 * c1 + c2 * c2 + c3 * c3)
    return {
        "h": float(np.max(hor)),
        "v": float(np.max(np.abs(c3))),
        "3d": float(np.max(mag)),
        "c1": float(np.max(np.abs(c1))),
        "c2": float(np.max(np.abs(c2))),
    }


def cumtrapz0(y, t):
    """Cumulative trapezoid integral of y(t), same length as y, starting at 0."""
    out = np.zeros_like(y)
    dt = np.diff(t)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * dt)
    return out


# ----------------------------------------------------------------------------- receivers
def parse_receiver(path):
    """Return (t, v[N,3], coord[3]) from a SeisSol receiver .dat file."""
    x = [None, None, None]
    with open(path) as fh:
        for line in fh:
            m = re.match(r"#\s*x([123])\s+([-+0-9.eE]+)", line)
            if m:
                x[int(m.group(1)) - 1] = float(m.group(2))
            if not line.startswith(("#", "TITLE", "VARIABLES", " ")) and line.strip():
                pass
    data = np.loadtxt(path, comments="#", skiprows=2)  # skip TITLE + VARIABLES lines
    t = data[:, 0]
    v = data[:, 7:10]  # v1,v2,v3
    return t, v, np.array(x, dtype=float)


def run_receivers(args):
    files = sorted(glob.glob(os.path.join(args.dir, f"{args.prefix}-receiver-*.dat")))
    if not files:
        sys.exit(f"no {args.prefix}-receiver-*.dat files in {args.dir}")

    rows = []
    header = ("receiver,x1,x2,x3,"
              "PGV_h,PGV_v,PGV_3d,"
              "PGA_h_g,PGA_v_g,PGA_3d_g,PGA_h_ms2,"
              "PGD_h,PGD_v,PGD_3d,"
              "perm_disp_h,perm_disp_v")
    print(f"{'recv':>5} {'x1':>12} {'x2':>14} | "
          f"{'PGV_h':>8} {'PGV_v':>8} | {'PGA_h[g]':>9} {'PGA_v[g]':>9} | "
          f"{'PGD_h':>8} {'PGD_v':>8} {'perm_h':>8}")
    for f in files:
        t, v, x = parse_receiver(f)
        rid = re.search(r"receiver-(\d+)", os.path.basename(f)).group(1)

        # PGV directly from velocity
        pv = peaks(v[:, 0], v[:, 1], v[:, 2])

        # PGA from time-derivative of velocity
        a = np.gradient(v, t, axis=0)
        pa = peaks(a[:, 0], a[:, 1], a[:, 2])

        # PGD from time-integral of velocity (displacement)
        d = np.column_stack([cumtrapz0(v[:, i], t) for i in range(3)])
        if args.detrend_disp:  # remove linear baseline (drops permanent offset!)
            for i in range(3):
                d[:, i] -= np.polyval(np.polyfit(t, d[:, i], 1), t)
        pd = peaks(d[:, 0], d[:, 1], d[:, 2])
        perm_h = float(np.hypot(d[-1, 0], d[-1, 1]))  # final (static) offset
        perm_v = float(abs(d[-1, 2]))

        rows.append(
            f"{rid},{x[0]:.3f},{x[1]:.3f},{x[2]:.3f},"
            f"{pv['h']:.6e},{pv['v']:.6e},{pv['3d']:.6e},"
            f"{pa['h']/G:.6e},{pa['v']/G:.6e},{pa['3d']/G:.6e},{pa['h']:.6e},"
            f"{pd['h']:.6e},{pd['v']:.6e},{pd['3d']:.6e},"
            f"{perm_h:.6e},{perm_v:.6e}")
        print(f"{rid:>5} {x[0]:>12.1f} {x[1]:>14.1f} | "
              f"{pv['h']:>8.3f} {pv['v']:>8.3f} | {pa['h']/G:>9.3f} {pa['v']/G:>9.3f} | "
              f"{pd['h']:>8.3f} {pd['v']:>8.3f} {perm_h:>8.3f}")

    out = args.out or os.path.join(args.dir, "pgm_receivers.csv")
    with open(out, "w") as fh:
        fh.write(header + "\n")
        fh.write("\n".join(rows) + "\n")
    print(f"\nwrote {out}  ({len(rows)} receivers)")
    if args.detrend_disp:
        print("NOTE: --detrend-disp removed the linear baseline -> PGD excludes the "
              "permanent tectonic offset.")


# ----------------------------------------------------------------------------- surface
def load_cell_field(d, name, nstep, ncell):
    arr = np.fromfile(os.path.join(d, name + ".bin"), dtype="<f8")
    if arr.size != nstep * ncell:
        sys.exit(f"{name}.bin has {arr.size} values, expected {nstep*ncell}")
    return arr.reshape(nstep, ncell)


def run_surface(args):
    xdmf = os.path.join(args.dir, f"{args.prefix}-surface.xdmf")
    if not os.path.exists(xdmf):
        sys.exit(f"missing {xdmf}")
    txt = open(xdmf).read()
    ncell = int(re.search(r'TopologyType="Triangle" NumberOfElements="(\d+)"', txt).group(1))
    nvert = int(re.search(r'GeometryType="XYZ" NumberOfElements="(\d+)"', txt).group(1))
    nstep = txt.count("<Time Value")
    times = [float(m) for m in re.findall(r'<Time Value="([^"]+)"', txt)]
    dt = times[1] - times[0] if len(times) > 1 else 0.0
    cell = os.path.join(args.dir, f"{args.prefix}-surface_cell", "mesh0")
    vert = os.path.join(args.dir, f"{args.prefix}-surface_vertex", "mesh0")
    print(f"cells={ncell} verts={nvert} steps={nstep} dt={dt} s")

    # peak velocity / displacement per cell over all time steps
    v = [load_cell_field(cell, f"v{i}", nstep, ncell) for i in (1, 2, 3)]
    u = [load_cell_field(cell, f"u{i}", nstep, ncell) for i in (1, 2, 3)]

    pgv_h = np.sqrt(v[0] ** 2 + v[1] ** 2).max(axis=0)
    pgv_3d = np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2).max(axis=0)
    pgv_v = np.abs(v[2]).max(axis=0)
    pgd_h = np.sqrt(u[0] ** 2 + u[1] ** 2).max(axis=0)
    pgd_3d = np.sqrt(u[0] ** 2 + u[1] ** 2 + u[2] ** 2).max(axis=0)
    pgd_v = np.abs(u[2]).max(axis=0)
    perm_h = np.hypot(u[0][-1], u[1][-1])  # final static offset map
    # coarse PGA from dt=0.5 s velocity (aliased - rough only)
    a = [np.gradient(v[i], dt, axis=0) for i in range(3)]
    pga_h_g = np.sqrt(a[0] ** 2 + a[1] ** 2).max(axis=0) / G

    geom = np.fromfile(os.path.join(vert, "geometry.bin"), dtype="<f8").reshape(nvert, 3)
    conn = np.fromfile(os.path.join(cell, "connect.bin"), dtype="<i8").reshape(ncell, 3)

    out = args.out or os.path.join(args.dir, "pgm_surface.vtu")
    fields = {
        "PGV_h": pgv_h, "PGV_v": pgv_v, "PGV_3d": pgv_3d,
        "PGD_h": pgd_h, "PGD_v": pgd_v, "PGD_3d": pgd_3d,
        "perm_disp_h": perm_h, "PGA_h_g_coarse": pga_h_g,
    }
    write_vtu(out, geom, conn, fields)
    print(f"wrote {out}")
    print(f"  PGV_h  max={pgv_h.max():.3f} m/s   PGD_h max={pgd_h.max():.3f} m   "
          f"perm_h max={perm_h.max():.3f} m")
    print("NOTE: surface dt=0.5 s -> PGA_h_g_coarse is aliased; use receivers for true PGA.")


def write_vtu(path, points, cells, cell_fields):
    """Minimal ASCII VTU (UnstructuredGrid, VTK triangle type=5) with cell data."""
    npn, ncl = len(points), len(cells)
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" '
                'byte_order="LittleEndian">\n  <UnstructuredGrid>\n')
        f.write(f'    <Piece NumberOfPoints="{npn}" NumberOfCells="{ncl}">\n')
        f.write('      <Points>\n        <DataArray type="Float64" NumberOfComponents="3" '
                'format="ascii">\n')
        f.write(" ".join(f"{p:.6e}" for p in points.ravel()) + "\n")
        f.write("        </DataArray>\n      </Points>\n      <Cells>\n")
        f.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        f.write(" ".join(map(str, cells.ravel())) + "\n        </DataArray>\n")
        f.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        f.write(" ".join(str(3 * (i + 1)) for i in range(ncl)) + "\n        </DataArray>\n")
        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        f.write(" ".join("5" for _ in range(ncl)) + "\n        </DataArray>\n      </Cells>\n")
        f.write("      <CellData>\n")
        for name, arr in cell_fields.items():
            f.write(f'        <DataArray type="Float64" Name="{name}" format="ascii">\n')
            f.write(" ".join(f"{x:.6e}" for x in arr) + "\n        </DataArray>\n")
        f.write("      </CellData>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n")


# ----------------------------------------------------------------------------- cli
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("receivers", help="point PGV/PGA/PGD from receiver .dat files")
    pr.add_argument("dir")
    pr.add_argument("--prefix", default="safs")
    pr.add_argument("--out", default=None)
    pr.add_argument("--detrend-disp", action="store_true",
                    help="remove linear baseline from displacement (drops permanent offset)")
    pr.set_defaults(func=run_receivers)
    ps = sub.add_parser("surface", help="spatial PGV/PGD maps from surface field -> VTU")
    ps.add_argument("dir")
    ps.add_argument("--prefix", default="safs")
    ps.add_argument("--out", default=None)
    ps.set_defaults(func=run_surface)
    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
