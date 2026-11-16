# A line map behind the movies: coastline, state and country borders, lakes,
# a graticule and city names, drawn in the model's x, y (km).
#
# The region file (regions/*.json) says how the model frame sits on the Earth
# and what to show:
#   projection   kind "azimuthal_equidistant", lat0, lon0 (degrees), radius_km
#   window_km    [x0, x1, y0, y1] of the picture
#   lonlat_window [lon0, lon1, lat0, lat1]: map data outside it is dropped
#   graticule    lat_step, lon_step (degrees)
#   cities       [[name, lat, lon], ...]
# The map lines are Natural Earth 1:10m GeoJSON files (public domain), read
# from a folder filled once by fetch_natural_earth.sh (compute nodes have no
# internet).
import json
import os

import numpy as np

FILES = {"lakes": "ne_10m_lakes.geojson",
         "states": "ne_10m_admin_1_states_provinces_lines.geojson",
         "borders": "ne_10m_admin_0_boundary_lines_land.geojson",
         "coast": "ne_10m_coastline.geojson"}
STYLE = {"lakes": ("#7aa6c2", 0.6), "states": ("#999999", 0.7), "borders": ("#666666", 0.9),
         "coast": ("#33566b", 0.9)}
BACKGROUND = "#f4f1ea"


class Basemap:
    def __init__(self, region_file, data_dir):
        r = json.load(open(region_file))
        p = r["projection"]
        if p["kind"] != "azimuthal_equidistant":
            raise SystemExit("%s: projection kind %s is not known" % (region_file, p["kind"]))
        self.R = float(p["radius_km"])
        self.p0, self.l0 = np.radians(p["lat0"]), np.radians(p["lon0"])
        self.window = tuple(float(v) for v in r["window_km"])
        self.lonlat = tuple(float(v) for v in r["lonlat_window"])
        self.cities = r.get("cities", [])
        self.grat = r.get("graticule", {"lat_step": 1, "lon_step": 2})
        self.lines = {}
        for key, name in FILES.items():
            path = os.path.join(data_dir, name)
            if not os.path.exists(path):
                raise SystemExit("map data %s is missing: run fetch_natural_earth.sh %s on a "
                                 "node with internet" % (path, data_dir))
            self.lines[key] = self._lines(path)

    def project(self, lon, lat):
        """Degrees -> model x, y in km (azimuthal equidistant on a sphere)."""
        p = np.radians(np.asarray(lat, float))
        dl = np.radians(np.asarray(lon, float)) - self.l0
        c = np.arccos(np.clip(np.sin(self.p0) * np.sin(p) + np.cos(self.p0) * np.cos(p) * np.cos(dl), -1, 1))
        k = np.where(c == 0, 1.0, c / np.where(c == 0, 1.0, np.sin(c)))
        return (self.R * k * np.cos(p) * np.sin(dl),
                self.R * k * (np.cos(self.p0) * np.sin(p) - np.sin(self.p0) * np.cos(p) * np.cos(dl)))

    def _lines(self, path):
        lo0, lo1, la0, la1 = self.lonlat
        out = []

        def add(coords):
            a = np.asarray(coords, float)[:, :2]
            ok = (a[:, 0] >= lo0) & (a[:, 0] <= lo1) & (a[:, 1] >= la0) & (a[:, 1] <= la1)
            edge = np.flatnonzero(np.diff(np.r_[False, ok, False].astype(int)))
            for i, j in zip(edge[::2], edge[1::2]):  # runs of points inside the window
                if j - i >= 2:
                    out.append(np.column_stack(self.project(a[i:j, 0], a[i:j, 1])))

        for f in json.load(open(path))["features"]:
            g = f.get("geometry")
            if not g:
                continue
            t, c = g["type"], g["coordinates"]
            parts = {"LineString": [c], "MultiLineString": c, "Polygon": c,
                     "MultiPolygon": [ring for poly in c for ring in poly]}.get(t, [])
            for part in parts:
                add(part)
        return out

    def draw_under(self, ax, zorder=0):
        """Background colour and graticule: below the data."""
        from matplotlib.collections import LineCollection
        ax.set_facecolor(BACKGROUND)
        lo0, lo1, la0, la1 = self.lonlat
        g = []
        for la in np.arange(np.ceil(la0), la1 + 1e-9, self.grat["lat_step"]):
            lo = np.linspace(lo0, lo1, 200)
            g.append(np.column_stack(self.project(lo, np.full_like(lo, la))))
        for lo in np.arange(np.ceil(lo0 / self.grat["lon_step"]) * self.grat["lon_step"], lo1 + 1e-9,
                            self.grat["lon_step"]):
            la = np.linspace(la0, la1, 200)
            g.append(np.column_stack(self.project(np.full_like(la, lo), la)))
        ax.add_collection(LineCollection(g, colors="#b9b4a8", linewidths=0.5, linestyles=(0, (1, 2)),
                                         zorder=zorder + 0.1))

    def draw_over(self, ax, zorder=5, fontsize=7):
        """Map lines and city names: above the data."""
        import matplotlib.patheffects as pe
        from matplotlib.collections import LineCollection
        for key in ("lakes", "states", "borders", "coast"):
            colour, width = STYLE[key]
            ax.add_collection(LineCollection(self.lines[key], colors=colour, linewidths=width, zorder=zorder))
        for name, lat, lon in self.cities:
            x, y = self.project(lon, lat)
            ax.plot([x], [y], "o", ms=3, mfc="k", mec="white", mew=0.6, zorder=zorder + 1)
            ax.text(x + 5, y + 3, name, fontsize=fontsize, color="#222222", zorder=zorder + 1,
                    path_effects=[pe.withStroke(linewidth=2, foreground="white")])
        ax.text(0.988, 0.012, "map: Natural Earth", transform=ax.transAxes, fontsize=6, color="#8a857a",
                ha="right", va="bottom", zorder=zorder + 2)
