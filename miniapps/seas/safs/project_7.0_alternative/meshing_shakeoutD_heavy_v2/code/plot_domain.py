#!/usr/bin/env python3
"""plot_domain.py -- the ShakeOut-D MESH domain on the tiled topo basemap.

Draws exactly one thing: the 520 x 300 km rectangle the heavy mesh is being built
in, over the UTM 11N topo basemap already rendered for the domain comparison.

The rectangle is Olsen et al. (2009) ShakeOut-D, TRIMMED 80 km off its southeast
end.  The published 600 x 300 km box reaches into Sonora/Baja where MUSCAL -- the
velocity model the resolution gate is scored on -- has no data (1.25 % of columns,
lon -116.60..-113.94, lat 31.08..33.17).  Trimming to 520 km buys 100.0000 %
coverage at every depth and still leaves 61.5 km of absorbing margin beyond the
fault's southeast tip.
"""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

HERE = "/Users/chunhuizhao/Downloads/connect_with_bbp/shakeoutv1v2_comparison/domain_comparison"
MESH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Q = np.load(os.path.join(MESH, "build_tmp", "domain_corners_utm.npy"))   # W,N,E,S
W, N, E, S = Q
poly = np.array([W, N, E, S, W]) / 1e3

meta = json.load(open(os.path.join(HERE, "basemap_domains_utm11n_topo.json")))
ext = meta["extent_utm11n"]
img = mpimg.imread(os.path.join(HERE, "basemap_domains_utm11n_topo.png"))

L1 = np.linalg.norm(S - W) / 1e3
L2 = np.linalg.norm(N - W) / 1e3

plt.rcParams.update({"font.size": 11, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})
fig, ax = plt.subplots(figsize=(11, 10.5))
ax.imshow(img, extent=[ext[0]/1e3, ext[1]/1e3, ext[2]/1e3, ext[3]/1e3],
          origin="upper", interpolation="bilinear", zorder=0)
ax.plot(poly[:, 0], poly[:, 1], color="#D1495B", lw=3.2, zorder=5,
        label=f"ShakeOut-D mesh domain  {L1:.0f} x {L2:.0f} km x 80 km deep")
ax.fill(poly[:, 0], poly[:, 1], color="#D1495B", alpha=0.10, zorder=4)

# fault SURFACE TRACE -- every DR-triangulation edge at the free surface, so the
# multi-strand system and its branches show as they are.  A fitted polyline would
# merge strands; a footprint band is not a trace.
tr = os.path.join(MESH, "build_tmp", "fault_trace.npz")
if os.path.exists(tr):
    from matplotlib.collections import LineCollection
    seg = np.load(tr)["seg"] / 1e3
    ax.add_collection(LineCollection(seg, colors="#111111", linewidths=1.15, zorder=8))
    ax.plot([], [], color="#111111", lw=1.6,
            label="San Andreas fault system — surface trace")

from pyproj import Transformer
inv = Transformer.from_crs("EPSG:32611", "EPSG:4326", always_xy=True)
for nm, P in (("W", W), ("N", N), ("E", E), ("S", S)):
    lo, la = inv.transform(P[0], P[1])
    ax.plot(P[0]/1e3, P[1]/1e3, "o", ms=7, mfc="white", mec="#D1495B", mew=2.2, zorder=6)
    ax.annotate(f"{nm}\n{abs(lo):.4f}°W\n{la:.4f}°N", (P[0]/1e3, P[1]/1e3),
                textcoords="offset points", xytext=(9, 9), fontsize=9,
                zorder=7, bbox=dict(fc="white", ec="#D1495B", alpha=0.86, lw=1.0, pad=0.35))

ax.set_xlim(ext[0]/1e3, ext[1]/1e3); ax.set_ylim(ext[2]/1e3, ext[3]/1e3)

# Cities from the SINGLE canonical table used by every SAFS figure -- importing it
# rather than copying coordinates is the whole point of that module (a per-figure
# copy is exactly how two panels end up disagreeing about where Los Angeles is).
sys.path.insert(0, os.path.expanduser(
    "~/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/seisol_quakeworx/"
    "v3_under_construction/toolbox/combined_workflow/pgv_postprocess"))
try:
    import map_cities
    map_cities.verify_projection()
    drawn = map_cities.draw_cities(ax, fontsize=11, ms=6.5)
    print(f"  cities drawn: {', '.join(drawn)}")
except Exception as exc:
    print(f"  [cities] SKIPPED: {exc}")
ax.set_xlabel("UTM 11N easting (km)"); ax.set_ylabel("UTM 11N northing (km)")
ax.set_title(f"SAFS heavy mesh domain — ShakeOut-D trimmed to full MUSCAL coverage\n"
             f"{L1:.1f} × {L2:.1f} km, 80 km deep · {L1*L2:,.0f} km² · bearing 130°",
             fontsize=12.5)
ax.legend(loc="lower left", fontsize=10, framealpha=0.94)
ax.set_aspect("equal")
out = os.path.join(HERE, "shakeoutD_mesh_domain.png")
fig.tight_layout(); fig.savefig(out, dpi=150)
print(f"[plot] {out}")
print(f"  domain {L1:.3f} x {L2:.3f} km, area {L1*L2:,.0f} km2")
