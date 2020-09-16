/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

#include "ad_bv.h"

#include "adio.h"

#include <bv.h>

/* setting a file size is collective, so as with other drivers implementing
 * scalable resize, we will have only one process do the resize and broadcast
 * result to everyone else */
void ADIOI_BV_Resize(ADIO_File fd, ADIO_Offset size, int *error_code)
{
    int err, rank;
    static char myname[] = "ADIOI_BV_RESIZE";
    MPI_Comm_rank(fd->comm, &rank);

    if (rank == fd->hints->ranklist[0]) {
        err = bv_setsize(fd->fs_ptr, fd->filename, size);
    }

    MPI_Bcast (&err, 1, MPI_INT, fd->hints->ranklist[0], fd->comm);

    /* --BEGIN ERROR HANDLING-- */
    if (err == -1) {
        *error_code = ADIOI_Err_create_code(myname, fd->filename, errno);
        return;
    }
    
    /* --END ERROR HANDLING-- */

    *error_code = MPI_SUCCESS;
}


