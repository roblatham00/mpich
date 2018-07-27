/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_iread.c,v 1.3 2002/10/24 17:01:04 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"

/* ADIOI_LOGFS_IreadContig()
 *
 * Implemented by immediately calling ReadContig()
 */
void ADIOI_LOGFS_IreadContig(ADIO_File fd, void *buf, int count,
                             MPI_Datatype datatype, int file_ptr_type,
                             ADIO_Offset offset, ADIO_Request * request, int
                             *error_code)
{
    MPI_Count typesize;

    MPI_Status status;

    MPI_Type_size_x(datatype, &typesize);
    *error_code = MPI_SUCCESS;

    ADIOI_LOGFS_ReadContig(fd, buf, count, datatype, file_ptr_type, offset, &status, error_code);
    MPIO_Completed_request_create(&fd, count * typesize, error_code, request);

}

void ADIOI_LOGFS_IreadStrided(ADIO_File fd, void *buf, int count,
                              MPI_Datatype datatype, int file_ptr_type,
                              ADIO_Offset offset, ADIO_Request * request, int
                              *error_code)
{
    MPI_Status status;
    MPI_Count typesize;
    MPI_Type_size_x(datatype, &typesize);

    ADIOI_LOGFS_ReadStrided(fd, buf, count, datatype, file_ptr_type, offset, &status, error_code);
    MPIO_Completed_request_create(&fd, count * typesize, error_code, request);
}
