/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_resize.c,v 1.2 2002/10/24 17:01:05 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"
#include "logfs.h"

void ADIOI_LOGFS_Resize(ADIO_File fd, ADIO_Offset size, int *error_code)
{
    /* resize always works */
    *error_code = MPI_SUCCESS;

    logfs_resize(fd->fs_ptr, size);
}
