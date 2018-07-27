/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_read.c,v 1.6 2002/10/24 17:01:05 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"
#include "logfs.h"

/* here offset has been converted from etypes into bytes */
void ADIOI_LOGFS_ReadContig(ADIO_File fd, void *buf, int count,
                            MPI_Datatype datatype, int file_ptr_type,
                            ADIO_Offset offset, ADIO_Status * status, int
                            *error_code)
{
    int myrank, nprocs, datatype_size;
    int ret;

    *error_code = MPI_SUCCESS;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);
    MPI_Type_size(datatype, &datatype_size);

    if (file_ptr_type == ADIO_INDIVIDUAL)
        offset = fd->fp_ind;

    /* read independent */
    ret = logfs_readdata(fd->fs_ptr, buf, count, datatype, offset, 0, status);

    /* update indiv. filepointer (== in bytes) */
    if (file_ptr_type != ADIO_EXPLICIT_OFFSET) {
        offset = fd->fp_ind;
        fd->fp_ind += datatype_size * count;
        fd->fp_sys_posn = fd->fp_ind;
    }

    *error_code = ret;
}

void ADIOI_LOGFS_ReadStrided(ADIO_File fd, void *buf, int count,
                             MPI_Datatype datatype, int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status * status, int
                             *error_code)
{
    int myrank, nprocs;
    int ret;

    *error_code = MPI_SUCCESS;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    if (ADIO_INDIVIDUAL == file_ptr_type) {
        offset = fd->fp_ind - fd->disp;
        offset /= fd->etype_size;
    }

    /* now offset in etypes rel to displacement of view */
    /* Read data non collective */
    ret = logfs_readdata(fd->fs_ptr, buf, count, datatype, offset, 0, status);

    if (file_ptr_type == ADIO_INDIVIDUAL) {
        int size;
        int datasize;
        MPI_Aint extent;
        MPI_Type_extent(fd->filetype, &extent);
        MPI_Type_size(fd->filetype, &size);
        MPI_Type_size(datatype, &datasize);
        ADIOI_Assert(!((datasize * count) % size));
        fd->fp_ind = offset + extent * ((datasize * count) / size);
    }

    *error_code = ret;
}
