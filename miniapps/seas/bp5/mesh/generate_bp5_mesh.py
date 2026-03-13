#!/usr/bin/env python3
"""Generate BP5 mesh with region-conformal surfaces.

Usage:
    conda activate pythonenv
    python generate_bp5_mesh.py --res 1.0 -o bp5_v2_1000m.msh
    python generate_bp5_mesh.py --res 5.0 -o bp5_v2_5000m.msh
"""
import subprocess
import argparse
import os

def generate(res_f_km, output, geo_file="bp5_v2.geo"):
    """Generate BP5 mesh using Gmsh."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    geo_path = os.path.join(script_dir, geo_file)
    if not os.path.exists(geo_path):
        raise FileNotFoundError(f"Geo file not found: {geo_path}")
    cmd = ["gmsh", "-3", geo_path,
           "-setnumber", "res_f", str(res_f_km),
           "-o", output]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"Generated: {output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BP5 mesh")
    parser.add_argument("--res", type=float, default=1.0,
                        help="Fault resolution in km (default: 1.0)")
    parser.add_argument("-o", "--output", default="bp5_v2_1000m.msh",
                        help="Output mesh file (default: bp5_v2_1000m.msh)")
    parser.add_argument("--geo", default="bp5_v2.geo",
                        help="Geo file to use (default: bp5_v2.geo)")
    args = parser.parse_args()
    generate(args.res, args.output, args.geo)
