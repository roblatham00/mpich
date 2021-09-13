/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "ad_unify.h"
#include "ad_unify_common.h"

/* as with ADIOI_PVFS2_Flush, implement the resize operation in a scalable
 * manner. one process does the work, then broadcasts the result to everyone
 * else.  fortunately, this operation is defined to be collective */
void ADIOI_UNIFY_Resize(ADIO_File fd, ADIO_Offset size, int *error_code)
{
    int ret, rank;
    ADIOI_UNIFY_fs *unifyfs_blob = fd->fs_ptr;

    *error_code = MPI_SUCCESS;

    unifyfs_io_request truncate_req = {
        .gfid = unifyfs_blob->gfid,
        .offset = size,
        .op = UNIFYFS_IOREQ_OP_TRUNC,
    };

    MPI_Comm_rank(fd->comm, &rank);

    /* MPI-IO semantics treat conflicting MPI_File_set_size requests the
     * same as conflicting write requests. Thus, a resize from one
     * process does not have to be visible to the other processes until a
     * synchronization point is reached */

    if (rank == fd->hints->ranklist[0]) {
        ret = unifyfs_dispatch_io(unifyfs_blob->fshdl, 1, &truncate_req);
        if (ret != UNIFYFS_SUCCESS) {
            *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                               MPIR_ERR_RECOVERABLE,
                                               __func__, __LINE__, MPI_ERR_IO, "Error in Dispatch:",
                                               "%s", strerror(ret));
            goto exit;
        }
        ret = unifyfs_wait_io(unifyfs_blob->fshdl, 1, &truncate_req, 1);
        if (ret != UNIFYFS_SUCCESS) {
            *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                               MPIR_ERR_RECOVERABLE,
                                               __func__, __LINE__, MPI_ERR_IO, "Error in I/O Wait:",
                                               "%s", strerror(ret));
            goto exit;
        }
        if (truncate_req.result.error != UNIFYFS_SUCCESS) {
            *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                               MPIR_ERR_RECOVERABLE,
                                               __func__, __LINE__, MPI_ERR_IO,
                                               "Error in I/O operation:", "%s",
                                               strerror(truncate_req.result.error));
            ret = truncate_req.result.error;
            goto exit;
        }

    }
  exit:
    MPI_Bcast(&ret, 1, MPI_INT, fd->hints->ranklist[0], fd->comm);
    /* --BEGIN ERROR HANDLING-- */
    if (ret != 0) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           __func__, __LINE__,
                                           MPI_ERR_IO, "Error in Unify Truncate", 0);
        return;
    }
    /* --END ERROR HANDLING-- */
}
