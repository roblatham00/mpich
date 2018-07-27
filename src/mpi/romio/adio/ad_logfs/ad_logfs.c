/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_logfs.c,v 1.2 2002/10/24 17:01:03 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"

/* adioi.h has the ADIOI_Fns_struct define */
#include "adioi.h"

struct ADIOI_Fns_struct ADIO_LOGFS_operations = {
    ADIOI_LOGFS_Open,   /* Open */
    ADIOI_GEN_OpenColl, /* OpenColl */
    ADIOI_LOGFS_ReadContig,     /* ReadContig */
    ADIOI_LOGFS_WriteContig,    /* WriteContig */
    ADIOI_LOGFS_ReadStridedColl,        /* ReadStridedColl */
    ADIOI_LOGFS_WriteStridedColl,       /* WriteStridedColl */
    ADIOI_LOGFS_SeekIndividual, /* SeekIndividual */
    ADIOI_LOGFS_Fcntl,  /* Fcntl */
    ADIOI_LOGFS_SetInfo,        /* SetInfo */
    ADIOI_LOGFS_ReadStrided,    /* ReadStrided */
    ADIOI_LOGFS_WriteStrided,   /* WriteStrided */
    ADIOI_LOGFS_Close,  /* Close */
    ADIOI_LOGFS_IreadContig,    /* IreadContig */
    ADIOI_LOGFS_IwriteContig,   /* IwriteContig */
    ADIOI_GEN_IODone,   /* ReadDone */
    ADIOI_GEN_IODone,   /* WriteDone */
    ADIOI_GEN_IOComplete,       /* ReadComplete */
    ADIOI_GEN_IOComplete,       /* WriteComplete */
    ADIOI_LOGFS_IreadStrided,   /* IreadStrided */
    ADIOI_LOGFS_IwriteStrided,  /* IwriteStrided */
    ADIOI_LOGFS_Flush,  /* Flush */
    ADIOI_LOGFS_Resize, /* Resize */
    ADIOI_LOGFS_Delete, /* Delete */
    ADIOI_LOGFS_Feature,        /* Features */
    "LOGFS: logging layer for ROMIO drivers",
    ADIOI_GEN_IreadStridedColl, /* IreadStridedColl */
    ADIOI_GEN_IwriteStridedColl /* IwriteStridedColl */
};
