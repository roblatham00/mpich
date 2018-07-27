/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_write.c,v 1.5 2002/10/24 17:01:06 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"
#include "logfs.h"
#include <assert.h>

void ADIOI_LOGFS_WriteContig(ADIO_File fd, const void *buf, int count,
                             MPI_Datatype datatype, int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status * status, int
                             *error_code)
{
    int myrank, nprocs, datatype_size;
    /*int etype_extent; */

    *error_code = MPI_SUCCESS;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);
    MPI_Type_size(datatype, &datatype_size);
    /*   FPRINTF(stdout, "[%d/%d] ADIOI_LOGFS_WriteContig called on %s\n", myrank,
     * nprocs, fd->filename);
     * FPRINTF(stdout, "Filetype: fd->file_system: %i", fd->file_system);
     */

    if (file_ptr_type == ADIO_INDIVIDUAL)
        offset = fd->fp_ind;

    /* convert offset into view (including disp+etype) */
    offset -= fd->disp;
    offset /= fd->etype_size;

    logfs_writedata(fd->fs_ptr, buf, count, datatype, offset, 0);

    if (file_ptr_type == ADIO_INDIVIDUAL)
        fd->fp_ind += datatype_size * count;

    if (status) {
        int bufsize, size;
        /* Don't set status if it isn't needed */
        MPI_Type_size(datatype, &size);
        bufsize = size * count;
        MPIR_Status_set_bytes(status, datatype, bufsize);
    }
}

void ADIOI_LOGFS_WriteStrided(ADIO_File fd, const void *buf, int count,
                              MPI_Datatype datatype, int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    int myrank, nprocs;

    *error_code = MPI_SUCCESS;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);
    /*FPRINTF(stdout, "[%d/%d] ADIOI_LOGFS_WriteStrided called on %s\n",
     * myrank, nprocs, fd->filename);
     * FPRINTF(stdout, "Filetype: fd->file_system: %i", fd->file_system);
     */
    if (ADIO_INDIVIDUAL == file_ptr_type) {
        offset = fd->fp_ind - fd->disp;
        offset /= fd->etype_size;
    }

    /* now offset in etypes rel to displacement of view */
    logfs_writedata(fd->fs_ptr, buf, count, datatype, offset, 0);

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
