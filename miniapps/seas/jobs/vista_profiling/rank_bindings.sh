#!/bin/bash
# S1 placement probe (slow-parts study, Chunhui Zhao 2026-10-25): one line per rank, "rank local-rank host
# cores", the cores this process may run on as the launcher left it (Cpus_allowed_list).
echo "${OMPI_COMM_WORLD_RANK:-?} ${OMPI_COMM_WORLD_LOCAL_RANK:-?} $(hostname -s) $(awk '/^Cpus_allowed_list/{print $2}' /proc/self/status)"
