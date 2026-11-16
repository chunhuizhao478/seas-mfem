#!/bin/bash
# The four Natural Earth 1:10m line files of basemap.py (public domain,
# naturalearthdata.com; 38.5 MB) into the folder $1.  Run once on a node with
# internet (a login node); compute nodes read the folder.
set -e
DIR=${1:?usage: fetch_natural_earth.sh FOLDER}
BASE=https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson
mkdir -p "$DIR"
for f in ne_10m_coastline ne_10m_admin_1_states_provinces_lines ne_10m_admin_0_boundary_lines_land ne_10m_lakes; do
   [ -s "$DIR/$f.geojson" ] || curl -fsS -o "$DIR/$f.geojson" "$BASE/$f.geojson"
   echo "$DIR/$f.geojson $(stat -c %s "$DIR/$f.geojson") bytes"
done
