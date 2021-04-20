/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "ad_unify.h"
#include "ad_unify_common.h"
#include <assert.h>

void ADIOI_UNIFY_Close(ADIO_File fd, int *error_code)
{
    int ret;
    ADIOI_UNIFY_fs *unifyfs_blob = (ADIOI_UNIFY_fs *) fd->fs_ptr;

    ret = unifyfs_laminate(unifyfs_blob->fshdl, fd->filename);
    if (ret != 0) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           "ADIOI_UNIFY_Close", __LINE__,
                                           MPI_ERR_UNKNOWN, "Error in unify_laminate", 0);
    } else {
        *error_code = MPI_SUCCESS;
    }

    ADIOI_Free(unifyfs_blob);

}
