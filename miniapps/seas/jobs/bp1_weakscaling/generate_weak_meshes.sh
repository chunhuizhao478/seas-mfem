#!/bin/bash
# ============================================================================
# Generate weak scaling meshes for BP1
#
# This script runs LOCALLY (not on the cluster).
# Prerequisites: conda activate pythonenv
#
# The weak scaling formula: h = 0.30 * sqrt(4 / np)
# This keeps the number of elements per process approximately constant.
# ============================================================================

echo "=============================================="
echo " BP1 Weak Scaling Mesh Generation"
echo "=============================================="
echo ""
echo "Make sure you have activated the correct environment:"
echo "  conda activate pythonenv"
echo ""

# Navigate to the mesh directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MESH_DIR="${SCRIPT_DIR}/../../bp1/mesh"

if [ ! -d "$MESH_DIR" ]; then
    echo "ERROR: Mesh directory not found: $MESH_DIR"
    exit 1
fi

cd "$MESH_DIR"
echo "Working directory: $(pwd)"
echo ""

# Check that the geo file exists
if [ ! -f "bp1_selfsimilar.geo" ]; then
    echo "ERROR: bp1_selfsimilar.geo not found in $(pwd)"
    exit 1
fi

# Generate meshes for each np value
# h = 0.30 * sqrt(4/np) in km

echo "--- np=4, h=0.3000 km ---"
gmsh -2 bp1_selfsimilar.geo -o bp1_ss_weak_np4.msh -setnumber h 0.3000
echo "Element count for np4:"
gmsh -0 bp1_ss_weak_np4.msh 2>&1 | grep -i "element"
echo ""

echo "--- np=8, h=0.2121 km ---"
gmsh -2 bp1_selfsimilar.geo -o bp1_ss_weak_np8.msh -setnumber h 0.2121
echo "Element count for np8:"
gmsh -0 bp1_ss_weak_np8.msh 2>&1 | grep -i "element"
echo ""

echo "--- np=16, h=0.1500 km ---"
gmsh -2 bp1_selfsimilar.geo -o bp1_ss_weak_np16.msh -setnumber h 0.1500
echo "Element count for np16:"
gmsh -0 bp1_ss_weak_np16.msh 2>&1 | grep -i "element"
echo ""

echo "--- np=32, h=0.1061 km ---"
gmsh -2 bp1_selfsimilar.geo -o bp1_ss_weak_np32.msh -setnumber h 0.1061
echo "Element count for np32:"
gmsh -0 bp1_ss_weak_np32.msh 2>&1 | grep -i "element"
echo ""

echo "--- np=40, h=0.0949 km ---"
gmsh -2 bp1_selfsimilar.geo -o bp1_ss_weak_np40.msh -setnumber h 0.0949
echo "Element count for np40:"
gmsh -0 bp1_ss_weak_np40.msh 2>&1 | grep -i "element"
echo ""

echo "--- np=80, h=0.0671 km ---"
gmsh -2 bp1_selfsimilar.geo -o bp1_ss_weak_np80.msh -setnumber h 0.0671
echo "Element count for np80:"
gmsh -0 bp1_ss_weak_np80.msh 2>&1 | grep -i "element"
echo ""

echo "=============================================="
echo " Mesh generation complete!"
echo "=============================================="
echo ""
echo "Generated meshes in: $(pwd)"
ls -lh bp1_ss_weak_np*.msh 2>/dev/null
echo ""
echo "Remember to scp the meshes to the cluster, e.g.:"
echo "  scp bp1_ss_weak_np*.msh <user>@<cluster>:<path>/bp1/mesh/"
