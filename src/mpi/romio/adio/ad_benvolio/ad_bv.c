/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *
 *   Copyright (C) 1997 University of Chicago.
 *   Copyright (C) 2017 DataDirect Networks.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_bv.h"

/* adioi.h has the ADIOI_Fns_struct define */
#include "adioi.h"

struct ADIOI_Fns_struct ADIO_BV_operations = {
    ADIOI_BV_Open,      /* Open */
    ADIOI_SCALEABLE_OpenColl, /* OpenColl */ /*XXX*/
        ADIOI_BV_ReadContig,    /* ReadContig */
    ADIOI_BV_WriteContig,       /* WriteContig */
    /* no collectives for BV: same as indep */
    ADIOI_BV_ReadStrided,       /* ReadStridedColl */
    ADIOI_BV_WriteStrided,      /* WriteStridedColl */
    ADIOI_GEN_SeekIndividual,   /* SeekIndividual */
    /* custom fcntl to obtain file size */
    ADIOI_BV_Fcntl,     /* Fcntl */
    ADIOI_GEN_SetInfo,  /* SetInfo */
    ADIOI_BV_ReadStrided,       /* ReadStrided */
    ADIOI_BV_WriteStrided,      /* WriteStrided */
    ADIOI_BV_Close,     /* Close */
    ADIOI_FAKE_IreadContig,     /* IreadContig */
    ADIOI_FAKE_IwriteContig,    /* IwriteContig */
    ADIOI_FAKE_IODone,  /* ReadDone */
    ADIOI_FAKE_IODone,  /* WriteDone */
    ADIOI_FAKE_IOComplete,      /* ReadComplete */
    ADIOI_FAKE_IOComplete,      /* WriteComplete */
    ADIOI_FAKE_IreadStrided,    /* IreadStrided */
    ADIOI_FAKE_IwriteStrided,   /* IwriteStrided */
    ADIOI_BV_Flush,    /* Flush */
    ADIOI_BV_Resize,   /* Resize */
    ADIOI_BV_Delete,    /* Delete */
    ADIOI_BV_Feature,
    "BENVOLIO: ROMIO + Mochi",
    ADIOI_GEN_IreadStridedColl,
    ADIOI_GEN_IwriteStridedColl
};
