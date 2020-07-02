/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *
 *   Copyright (C) 1997 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_bv.h"
#include "adio.h"
#include "ad_bv_common.h"
#include <bv.h>

void ADIOI_BV_Flush(ADIO_File fd, int *error_code)
{
    int err;
    static char myname[] = "ADIOI_BV_FLUSH";

    /* the deferred-open optimization may mean that a file has not been opened
     * on this processor */
    if (fd->is_open > 0) {
        err = bv_flush(fd->fs_ptr, fd->filename);
        /* --BEGIN ERROR HANDLING-- */
        if (err != 0) {
            *error_code = MPIO_Err_create_code(MPI_SUCCESS, MPIR_ERR_RECOVERABLE,
                                               myname, __LINE__, MPI_ERR_IO,
                                               "Benvolio error", 0);
            return;
        }
        /* --END ERROR HANDLING-- */
    }

    *error_code = MPI_SUCCESS;
}
