/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_iwrite.c,v 1.3 2002/10/24 17:01:04 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"

/**
 *
 *  Could improve this:
 *    add functions to writebuf to free up memory; in
 *    the start function, inform writebuf that it should try to free at least
 *    this many bytes; in the done function, write dump it to the writebuf;
 *    either free or not
 *
 *  Hmm; this is not needed; calling the progress function on the writebuf
 *  already tries to free mem; it would be enough to call the progress
 *  function in the start func
 */



/* ADIOI_LOGFS_IwriteContig()
 *
 * Implemented by immediately calling WriteContig()
 */
void ADIOI_LOGFS_IwriteContig(ADIO_File fd, const void *buf, int count,
                              MPI_Datatype datatype, int file_ptr_type,
                              ADIO_Offset offset, ADIO_Request * request, int
                              *error_code)
{
    ADIO_Status status;
    MPI_Offset len;
    MPI_Count typesize;
    int myrank, nprocs;

    MPI_Type_size_x(datatype, &typesize);
    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    len = count * typesize;
    ADIOI_LOGFS_WriteContig(fd, buf, count, datatype, file_ptr_type, offset, &status, error_code);

    MPIO_Completed_request_create(&fd, len, error_code, request);
}

void ADIOI_LOGFS_IwriteStrided(ADIO_File fd, const void *buf, int count,
                               MPI_Datatype datatype, int file_ptr_type,
                               ADIO_Offset offset, ADIO_Request * request, int
                               *error_code)
{
    ADIO_Status status;
    MPI_Count typesize;

    MPI_Type_size_x(datatype, &typesize);

    ADIOI_LOGFS_WriteStrided(fd, buf, count, datatype, file_ptr_type, offset, &status, error_code);

    MPIO_Completed_request_create(&fd, typesize * count, error_code, request);
}
