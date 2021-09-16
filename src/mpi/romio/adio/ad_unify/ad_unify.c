/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "ad_unify.h"

#include "adio.h"

struct ADIOI_Fns_struct ADIO_UNIFY_operations = {
    ADIOI_UNIFY_Open,   /* Open */
    /* cannot use SCALEABLE version until Unify implements "gfid lookup" */
    ADIOI_GEN_OpenColl, /* OpenColl */
    ADIOI_UNIFY_ReadContig,     /* ReadContig */
    ADIOI_UNIFY_WriteContig,    /* WriteContig */
    ADIOI_GEN_ReadStridedColl,  /* ReadStridedColl */
    ADIOI_GEN_WriteStridedColl, /* WriteStridedColl */
    ADIOI_GEN_SeekIndividual,   /* SeekIndividual */
    ADIOI_UNIFY_Fcntl,  /* Fcntl */
    ADIOI_GEN_SetInfo,  /* SetInfo */
    ADIOI_GEN_ReadStrided,      /* ReadStrided */
    ADIOI_GEN_WriteStrided_naive,       /* WriteStrided */
    ADIOI_UNIFY_Close,  /* Close */
#ifdef HAVE_MPI_GREQUEST_EXTENSIONS
    ADIOI_GEN_IReadContig,      /* IreadContig */
    ADIOI_GEN_IWriteContig,     /* IwriteContig */
#else
    ADIOI_FAKE_IreadContig,     /* IreadContig */
    ADIOI_FAKE_IwriteContig,    /* IwriteContig */
#endif
    ADIOI_FAKE_IODone,  /* ReadDone */
    ADIOI_FAKE_IODone,  /* WriteDone */
    ADIOI_FAKE_IOComplete,      /* ReadComplete */
    ADIOI_FAKE_IOComplete,      /* WriteComplete */
    ADIOI_FAKE_IreadStrided,    /* IreadStrided */
    ADIOI_FAKE_IwriteStrided,   /* IwriteStrided */
    ADIOI_UNIFY_Flush,  /* Flush */
    ADIOI_UNIFY_Resize, /* Resize */
    ADIOI_UNIFY_Delete, /* Delete */
    ADIOI_UNIFY_Feature,
    "UNIFY: the burst-buffer file system",
    ADIOI_GEN_IreadStridedColl, /* IreadStridedColl */
    ADIOI_GEN_IwriteStridedColl,        /* IwriteStridedColl */
#if defined(F_SETLKW64)
    ADIOI_GEN_SetLock   /* SetLock */
#else
    ADIOI_GEN_SetLock64 /* SetLock */
#endif
};
