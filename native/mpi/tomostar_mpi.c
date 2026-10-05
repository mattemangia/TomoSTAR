/* Copyright 2026 Matteo Mangiagalli; SPDX-License-Identifier: Apache-2.0 */
#include <mpi.h>
#include <stdint.h>

int ts_mpi_init(void) {
    int provided;
    int status = MPI_Init_thread(0, 0, MPI_THREAD_FUNNELED, &provided);
    if (status != MPI_SUCCESS) return status;
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    return provided >= MPI_THREAD_FUNNELED ? MPI_SUCCESS : -1;
}
int ts_mpi_rank(int *rank) { return MPI_Comm_rank(MPI_COMM_WORLD, rank); }
int ts_mpi_size(int *size) { return MPI_Comm_size(MPI_COMM_WORLD, size); }
int ts_mpi_bcast_length(int64_t *length, int root) {
    return MPI_Bcast(length, 1, MPI_INT64_T, root, MPI_COMM_WORLD);
}
int ts_mpi_bcast_bytes(unsigned char *data, int count, int root) {
    return MPI_Bcast(data, count, MPI_BYTE, root, MPI_COMM_WORLD);
}
int ts_mpi_abort(int code) { return MPI_Abort(MPI_COMM_WORLD, code); }
int ts_mpi_finalize(void) { return MPI_Finalize(); }
