
#include "ad_unify.h"
#include "ad_unify_common.h"

void ADIOI_UNIFY_Fence(ADIO_File fd, int *error_code)
{
    int ret;
    *error_code = MPI_SUCCESS;
    ADIOI_UNIFY_fs *unifyfs_blob = fd->fs_ptr;

    /* Synchronize with file system to flush writes
     * to make data visible to other processes. */
    unifyfs_io_request sync_req = {
        .op = UNIFYFS_IOREQ_OP_SYNC_META,
        .gfid = unifyfs_blob->gfid,
        .result.error = 0,
    };
    ret = unifyfs_dispatch_io(unifyfs_blob->fshdl, 1, &sync_req);
    if (ret == UNIFYFS_SUCCESS) {
        ret = unifyfs_wait_io(unifyfs_blob->fshdl, 1, &sync_req, 1);
        if (ret == UNIFYFS_SUCCESS)
            if (sync_req.result.error != 0) {
                *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                                   MPIR_ERR_RECOVERABLE, __func__, __LINE__,
                                                   MPI_ERR_IO, "**io", "**io %s",
                                                   strerror(sync_req.result.error));
            }
    }

    /* Synchronize processes so that upon returning from fence,
     * all procs know that all other procs have completed their flush. */
    MPI_Barrier(fd->comm);

    return;
}
