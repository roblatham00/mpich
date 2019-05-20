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

/* Benvolio clients do not have anything to flush/sync.  Instead we request a
 * flush from the remote end. Simliar to PVFS2, we can take advantage of
 * collective MPI_Flush and arrange it so only one process sends the flush
 * request to the remote side */
void ADIOI_BV_Flush(ADIO_File fd, int *error_code)
{
    int err, dummy = 0, dummy_in = 0, rank;

    static char myname[] = "ADIOI_BV_FLUSH";

    MPI_Comm_rank(fd->comm, &rank);
    /* We don't need a full barrier: just need to know everyone has entered
     * this collective (and finished their writes) */
    MPI_Reduce(&dummy_in, &dummy, 1, MPI_INT, MPI_SUM, fd->hints->ranklist[0], fd->comm);

    if (rank = fd->hints->ranklist[0]) {
        err = bv_flush(fd->fs_ptr, fd->filename);
    }
    MPI_Bcast(&err, 1, MPI_INT, fd->hints->ranklist[0], fd->comm);
    /* --BEGIN ERROR HANDLING-- */
    if (err != 0) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS, MPIR_ERR_RECOVERABLE,
                                           myname, __LINE__, MPI_ERR_IO, "Benvolio error", 0);
        return;
    }
    /* --END ERROR HANDLING-- */

    *error_code = MPI_SUCCESS;
}
