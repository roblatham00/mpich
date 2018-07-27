/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_delete.c,v 1.2 2002/10/24 17:01:03 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include <assert.h>

#include "ad_logfs.h"
#include "adioi.h"
#include "logfs.h"

/* delete can only be called in standalone mode
 * (unless we also modify the global delete to detect logfs files) */
void ADIOI_LOGFS_Delete(const char *filename, int *error_code)
{
    *error_code = MPI_SUCCESS;

    /* logfs_delete removes all the logfs ancilary files (.logfs, .meta,
     * .data) related to 'filename', but not the actual file */
    logfs_delete(filename);

    /* MPI_File_delete already stripped the prefix.  TODO: How do we avoid a
     * loop if we are using logfs as a data representation ? */
    *error_code = MPI_File_delete(filename, MPI_INFO_NULL);

}
