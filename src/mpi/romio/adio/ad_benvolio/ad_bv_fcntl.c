/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *
 *   Copyright (C) 2019 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_bv.h"
#include <bv.h>

void ADIOI_BV_Fcntl(ADIO_File fd, int flag, ADIO_Fcntl_t * fcntl_struct, int *error_code)
{
    ADIO_Offset fsize;

    switch (flag) {
        case ADIO_FCNTL_GET_FSIZE:
            fsize = bv_getsize(fd->fs_ptr, fd->filename);
            if (fsize < 0) {
                /* --BEGIN ERROR HANDLING-- */
                *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                                   MPIR_ERR_RECOVERABLE,
                                                   __func__, __LINE__,
                                                   fsize, "Error in bv_getsize", 0);
                /* --END ERROR HANDLING-- */
            } else {
                *error_code = MPI_SUCCESS;
            }
            fcntl_struct->fsize = fsize;
            return;
        case ADIO_FCNTL_SET_DISKSPACE:
            ADIOI_GEN_Prealloc(fd, fcntl_struct->diskspace, error_code);
            break;
        case ADIO_FCNTL_SET_ATOMICITY:
        default:
            /* --BEGIN ERROR HANDLING-- */
            *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                               MPIR_ERR_RECOVERABLE,
                                               __func__, __LINE__,
                                               MPI_ERR_ARG, "**flag", "**flag %d", flag);
            /* --END ERROR HANDLING-- */
    }
}
