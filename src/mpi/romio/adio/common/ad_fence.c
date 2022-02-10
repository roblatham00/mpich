/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "adio.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

void ADIOI_GEN_Fence(ADIO_File fd, int *error_code)
{
    int err;
    static char myname[] = "ADIOI_GEN_FENCE";

    /* For the general (common) case, we assume that the underlying
     * file system implements a POSIX-compliant write,
     * meaning that written data becomes visible to other processes
     * immediately upon returning from a successful write call.
     *
     * For such file systems, there is no need to explicitly flush
     * data or to synchronize with the file system to expose the
     * data from those write calls. */

    /* Even though the data flush step is a NOP,
     * fence semantics still require that we synchronize processes. */
    MPI_Barrier(fd->comm);

    *error_code = MPI_SUCCESS;
}
