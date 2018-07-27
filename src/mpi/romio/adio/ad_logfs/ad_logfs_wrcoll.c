/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_wrcoll.c,v 1.4 2002/10/24 17:01:06 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"
#include "logfs.h"
#include <assert.h>

void ADIOI_LOGFS_WriteStridedColl(ADIO_File fd, const void *buf, int count,
                                  MPI_Datatype datatype, int file_ptr_type,
                                  ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    int myrank, nprocs;

    *error_code = MPI_SUCCESS;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    if (ADIO_INDIVIDUAL == file_ptr_type) {
        /* fp->ind is in bytes ignoring view */
        offset = fd->fp_ind - fd->disp;
        offset /= fd->etype_size;
    }

    /* now offset in etypes rel to displacement of view */
    logfs_writedata(fd->fs_ptr, buf, count, datatype, offset, 1);

    if (file_ptr_type == ADIO_INDIVIDUAL) {
        int size;
        int datasize;
        MPI_Aint extent;
        MPI_Type_extent(fd->filetype, &extent);
        MPI_Type_size(fd->filetype, &size);
        MPI_Type_size(datatype, &datasize);
        assert(!((datasize * count) % size));
        fd->fp_ind = offset + extent * ((datasize * count) / size);
    }

    if (status) {
        int bufsize, size;
        /* Don't set status if it isn't needed */
        MPI_Type_size(datatype, &size);
        bufsize = size * count;
        MPIR_Status_set_bytes(status, datatype, bufsize);
    }
}
