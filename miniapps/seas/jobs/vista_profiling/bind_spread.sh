#!/bin/bash
# S1 launcher (slow-parts study, Chunhui Zhao 2026-10-25): binds this rank to core
# local rank x (cores of the node / ranks on the node), so 72 ranks on a 144-core Vista node take every
# second core and fill both Grace CPUs (36 each, cores 0-71 and 72-143); 144 ranks get core = local rank.
lr=${OMPI_COMM_WORLD_LOCAL_RANK:?bind_spread.sh: not started by Open MPI}
ls=${OMPI_COMM_WORLD_LOCAL_SIZE:?bind_spread.sh: not started by Open MPI}
exec taskset -c $(( lr * $(nproc --all) / ls )) "$@"
