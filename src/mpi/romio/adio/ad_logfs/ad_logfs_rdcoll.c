/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_rdcoll.c,v 1.4 2002/10/24 17:01:05 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"
#include "logfs.h"

void ADIOI_LOGFS_ReadStridedColl(ADIO_File fd, void *buf, int count,
                                 MPI_Datatype datatype, int file_ptr_type,
                                 ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    int ret;

    *error_code = MPI_SUCCESS;

    if (ADIO_INDIVIDUAL == file_ptr_type) {
        /* fp->ind is in bytes ignoring view */
        offset = fd->fp_ind - fd->disp;
        offset /= fd->etype_size;
    }

    /* now offset in etypes rel to displacement of view */
    /* collective read */
    ret = logfs_readdata(fd->fs_ptr, buf, count, datatype, offset, 1, status);

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
