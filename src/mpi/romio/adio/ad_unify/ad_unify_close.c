/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "ad_unify.h"
#include "ad_unify_common.h"
#include <assert.h>

/* What does it mean to "close" a stateless UnifyFS file?
 * - most straightforwrad: send a SYNC_META operation
 * - "Laminate" is not a great fit: "laminate" means not only "i'm done with
 *    this file" but also "i will never ever write this again"
 * - future work (?): the i/o operations never wait for completion and instead
 *   do all the waiting in sync or close.  Requires some more plumbing -- no
 *   other ROMIO driver has a concept of operations living across ADIO calls */
void ADIOI_UNIFY_Close(ADIO_File fd, int *error_code)
{
    int ret;
    ADIOI_UNIFY_fs *unifyfs_blob = (ADIOI_UNIFY_fs *) fd->fs_ptr;

    unifyfs_io_request sync_req = {
        .op = UNIFYFS_IOREQ_OP_SYNC_META,
        .gfid = unifyfs_blob->gfid,
    };
    ret = unifyfs_dispatch_io(unifyfs_blob->fshdl, 1, &sync_req);
    if (ret == UNIFYFS_SUCCESS) {
        ret = unifyfs_wait_io(unifyfs_blob->fshdl, 1, &sync_req, 1);
        *error_code = MPI_SUCCESS;
        if (ret != UNIFYFS_SUCCESS) {
            *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                               MPIR_ERR_RECOVERABLE,
                                               "ADIOI_UNIFY_Close", __LINE__,
                                               MPI_ERR_UNKNOWN, "Error in unify wait request", "%s",
                                               strerror(ret));
        }
        if (sync_req.result.error != 0) {
            *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                               MPIR_ERR_RECOVERABLE,
                                               "ADIOI_UNIFY_Close", __LINE__,
                                               MPI_ERR_UNKNOWN, "Error in unify sync operation",
                                               "%s", strerror(ret));
        }
    } else {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           "ADIOI_UNIFY_Close", __LINE__,
                                           MPI_ERR_UNKNOWN, "Error in unify dispatch operation",
                                           "%s", strerror(ret));
    }
    ADIOI_Free(unifyfs_blob);
    fd->fs_ptr = NULL;
}
