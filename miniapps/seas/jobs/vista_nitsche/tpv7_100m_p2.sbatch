#!/bin/bash
#SBATCH -J tpv7_100m_p2
#SBATCH -o tpv7_100m_p2_%j.out
#SBATCH -e tpv7_100m_p2_%j.err
#SBATCH -p gg
#SBATCH -N 32
#SBATCH -n 4608
#SBATCH -t 03:00:00
#SBATCH -A EAR24034
# Nitsche dynamic rupture, tpv7 on the 100 m mesh, p = 2, 0 -> 12.0 s, Kelvin-Voigt eta_bar 0.1 +
# Prakash-Clifton as committed in tpv7/configs/tpv7_nitsche.toml (Chunhui Zhao, TPV6/TPV7/TPV31
# plan rev 1.1, 2026-10-19: 100 m, the node spacing the description asks for; layout
# revised 2026-10-25 for the rank 0 mesh setup).
# Submit from miniapps/seas/jobs/vista_nitsche/ AFTER build_verify_vista.sbatch passed at this commit.
# Mesh: tpv6/mesh/tpv6_nitsche_100m.msh, made on Vista by build_verify_vista.sbatch (TPV6 and TPV7
# share it).  The cloud gmsh 4.15.2 mesh of the same .geo (built-in kernel): 6,043,707 tets, min
# 6V/A 14.4 m on the near (6000 m/s) side, so dt is 0.469 x the TPV205 200 m dt at every p.
# 144 ranks per node: since the rank 0 mesh setup (trunk bc91ac91, merged at 746b9c4f) only rank 0
# reads the serial mesh.  Memory is now set by the first time step's operator assembly: the TPV31
# 200 m p = 2 run on 16 x 144 (Vista 1096105, 126.9 k tets per node) peaked at 51.8 GiB committed per
# node, about 21 GiB plus 0.25 MiB per tet, so 67-77 GiB here (188.9 k tets per node) of the 118.8 GiB
# commit limit (Vista compute nodes refuse allocations past it).
# Estimate: 199,000 steps to 12.0 s (TPV7 p = 1 took 93,447, Vista 1085892; TPV31 p = 2 took 2.13 x
# its p = 1 steps).  Step cost from TPV31 p = 2 (11.19 ms per step on the slowest rank) scaled by tets
# per node: 16.7 ms per step, 0.9 h, plus setup and fault output, about 1.1 h: 32 nodes x 0.33 SU =
# 12 SU.
# Requested: 32 nodes x 3 h x 0.33 SU = 32 SU at most.
PROB=tpv7
RES=100m
DECK=tpv7/configs/tpv7_nitsche.toml
MESH=tpv6/mesh/tpv6_nitsche_100m.msh
ORDER=2
TFINAL=12.0
source "$SLURM_SUBMIT_DIR/run_nitsche_vista.sh"
