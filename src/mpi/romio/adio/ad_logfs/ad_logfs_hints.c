/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_logfs.hints.c,v 1.4 2002/10/24 17:01:04 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"
#include "logfs.h"

void ADIOI_LOGFS_SetInfo(ADIO_File fd, MPI_Info users_info, int *error_code)
{
    *error_code = MPI_SUCCESS;


    /* In standalone mode we have complete control over hints, so just
     * process the users' info structure and return;
     * However, somehow we still need to call the gen_setinfo
     * since otherwise there are segmentation faults in other parts
     * of the code that depend on certaing things to be set
     * (two-phase for example) */

    /* this modifies fd->info */
    ADIOI_GEN_SetInfo(fd, users_info, error_code);
    logfs_transfer_hints(users_info, fd->info);

    *error_code = MPI_SUCCESS;
    return;

}
