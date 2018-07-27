/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_setsh.c,v 1.2 2002/10/24 17:01:05 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"

void ADIOI_LOGFS_Set_shared_fp(ADIO_File fd, ADIO_Offset offset, int *error_code)
{
    int myrank, nprocs;

    *error_code = MPI_SUCCESS;


    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    /* what to do here ??? */
    ADIOI_Assert(0);

}
