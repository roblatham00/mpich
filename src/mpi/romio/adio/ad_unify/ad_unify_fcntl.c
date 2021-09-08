/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "ad_unify.h"
#include "ad_unify_common.h"
#include "adio.h"

void ADIOI_UNIFY_Fcntl(ADIO_File fd, int flag, ADIO_Fcntl_t * fcntl_struct, int *error_code)
{
    int ret;
    ADIOI_UNIFY_fs *unifyfs_blob = fd->fs_ptr;
    unifyfs_status u_status;

    switch (flag) {
        case ADIO_FCNTL_GET_FSIZE:
            ret = unifyfs_stat(unifyfs_blob->fshdl, unifyfs_blob->gfid, &u_status);
            if (ret != UNIFYFS_SUCCESS) {
                /* --BEGIN ERROR HANDLING-- */
                *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                                   MPIR_ERR_RECOVERABLE,
                                                   __func__, __LINE__,
                                                   MPI_ERR_IO, "Error in unify_stat", 0);
                /* --END ERROR HANDLING-- */
            } else {
                *error_code = MPI_SUCCESS;
            }
            fcntl_struct->fsize = u_status.global_file_size;
            return;

        case ADIO_FCNTL_SET_DISKSPACE:
            ADIOI_GEN_Prealloc(fd, fcntl_struct->diskspace, error_code);
            break;

            /* --BEGIN ERROR HANDLING-- */
        case ADIO_FCNTL_SET_ATOMICITY: /* atomic mode not supported */
        default:
            *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                               MPIR_ERR_RECOVERABLE,
                                               __func__, __LINE__,
                                               MPI_ERR_ARG, "**flag", "**flag %d", flag);
            /* --END ERROR HANDLING-- */
    }
}
