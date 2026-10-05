#!/usr/bin/env bash
set -euo pipefail
# Build against the MPI implementation used by the cluster launcher.
output_dir="${1:?Usage: bash native/mpi/build.sh OUTPUT_DIRECTORY}"
mkdir -p "$output_dir"
mpicc -O2 -fPIC -shared "$(dirname "$0")/tomostar_mpi.c" -o "$output_dir/libtomostar_mpi.so"
