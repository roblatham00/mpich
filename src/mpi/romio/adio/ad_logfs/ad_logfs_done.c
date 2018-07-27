/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_done.c,v 1.3 2002/10/24 17:01:03 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"

int ADIOI_LOGFS_ReadDone(ADIO_Request * request, ADIO_Status * status, int
                         *error_code)
{
    return 1;
}

int ADIOI_LOGFS_WriteDone(ADIO_Request * request, ADIO_Status * status, int
                          *error_code)
{
    return 1;
}
