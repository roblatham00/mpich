/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *
 *   Copyright (C) 1997 University of Chicago.
 *   Copyright (C) 2017 DataDirect Networks.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_mochio.h"

/* adioi.h has the ADIOI_Fns_struct define */
#include "adioi.h"

struct ADIOI_Fns_struct ADIO_MOCHIO_operations = {
    ADIOI_MOCHIO_Open,  /* Open */
    ADIOI_SCALEABLE_OpenColl, /* OpenColl */ /*XXX*/
        ADIOI_MOCHIO_ReadContig,        /* ReadContig */
    ADIOI_MOCHIO_WriteContig,   /* WriteContig */
    /* no collectives for MOCHIO: same as indep */
    ADIOI_MOCHIO_ReadStrided,   /* ReadStridedColl */
    ADIOI_MOCHIO_WriteStrided,  /* WriteStridedColl */
    ADIOI_GEN_SeekIndividual,   /* SeekIndividual */
    ADIOI_GEN_Fcntl,    /* Fcntl */
    ADIOI_GEN_SetInfo,  /* SetInfo */
    ADIOI_MOCHIO_ReadStrided,   /* ReadStrided */
    ADIOI_MOCHIO_WriteStrided,  /* WriteStrided */
    ADIOI_MOCHIO_Close, /* Close */
    ADIOI_FAKE_IreadContig,     /* IreadContig */
    ADIOI_FAKE_IwriteContig,    /* IwriteContig */
    ADIOI_FAKE_IODone,  /* ReadDone */
    ADIOI_FAKE_IODone,  /* WriteDone */
    ADIOI_FAKE_IOComplete,      /* ReadComplete */
    ADIOI_FAKE_IOComplete,      /* WriteComplete */
    ADIOI_FAKE_IreadStrided,    /* IreadStrided */
    ADIOI_FAKE_IwriteStrided,   /* IwriteStrided */
    ADIOI_GEN_Flush,    /* Flush */
    ADIOI_GEN_Resize,   /* Resize */
    ADIOI_GEN_Delete,   /* Delete */
    ADIOI_MOCHIO_Feature,
    "MOCHIO: ROMIO + Mochi",
    ADIOI_GEN_IreadStridedColl,
    ADIOI_GEN_IwriteStridedColl
};
