/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_wait.c,v 1.4 2002/10/24 17:01:05 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"

void ADIOI_LOGFS_ReadComplete(ADIO_Request * request, ADIO_Status * status, int
                              *error_code)
{
    return;
}

void ADIOI_LOGFS_WriteComplete(ADIO_Request * request, ADIO_Status * status, int
                               *error_code)
{
    return;
}
