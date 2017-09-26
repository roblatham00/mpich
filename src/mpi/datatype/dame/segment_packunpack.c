/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "dame.h"
#include "veccpy.h"

/* #define MPICH_DEBUG_SEGMENT_MOVE */
/* TODO: Consider integrating this with the general debug support. */
/* Note: This does not use the CVAR support for the environment variable
   because (a) this is a temporary code and (b) it is expert developer
   only */
#ifdef MPICH_DEBUG_SEGMENT_MOVE
static int printSegment = -1;
static void setPrint(void)
{
    char *s = getenv("MPICH_DATALOOP_PRINT");
    if (s && (strcmp(s, "yes") == 0 || strcmp(s, "YES") == 0)) {
        printSegment = 1;
    } else {
        printSegment = 0;
    }
}

#define DBG_SEGMENT(_a) do { if (printSegment < 0) setPrint();  \
        if (printSegment) { _a; } } while (0)
#else
#define DBG_SEGMENT(_a)
#endif


void MPIR_Segment_pack(DAME_Segment * segp, DAME_Offset first, DAME_Offset * lastp, void *streambuf)
{
    /* TODO: TLP: Use max(size_t) */
    size_t outSize = LONG_MAX;
    if (*lastp != SEGMENT_IGNORE_LAST)
        outSize = (*lastp - first);

    DAME_Offset size, extent;
    DAME_Handle_get_size_macro(segp->handle, size);
    DAME_Handle_get_extent_macro(segp->handle, extent);

    DAME_Offset totalcopied = 0;

    DAME_Dame *dl = NULL;
    DAME_Handle_get_loopptr_macro(segp->handle, dl);

    if (!dl) {
        /* TODO: If no dataloop is defined, it should be a basic type for which
         * size == extent. Need to confirm this though */
        MPIR_Assert(size == extent);
        DAME_Dame defDL[3];
        defDL[0].kind = DL_EXIT;
        defDL[1].kind = DL_CONTIGFINAL;
        defDL[1].count = segp->count;
        defDL[1].size = size * segp->count;
        defDL[1].extent = extent;
        defDL[1].s.c_t.basesize = size;
        defDL[1].s.c_t.baseextent = extent;
        defDL[2].kind = DL_BOTTOM;

        DL_pack(segp->ptr + segp->in_offset[SEGMENT_PACK_IDX],
                streambuf, outSize, defDL, 1, &segp->state[SEGMENT_PACK_IDX], &totalcopied);
        segp->curcount[SEGMENT_PACK_IDX] = totalcopied / size;
        segp->in_offset[SEGMENT_PACK_IDX] += (segp->curcount[SEGMENT_PACK_IDX] * size);

        goto exit;
    }

    DAME_Count i;
    for (i = segp->curcount[SEGMENT_PACK_IDX]; i < segp->count; i++) {
        int rv;
        DAME_Offset copied = 0;
        DAME_Offset outsz = outSize - totalcopied;

        rv = DL_pack(segp->ptr + segp->in_offset[SEGMENT_PACK_IDX],
                     streambuf + totalcopied,
                     outsz, dl, 0, &segp->state[SEGMENT_PACK_IDX], &copied);

        totalcopied += copied;
        /* TODO: TLP: Use the MPICH error handling */
        switch (rv) {
            case DAME_PACKUNPACK_COMPLETED:
                break;
            case DAME_PACKUNPACK_PARTIAL:
                segp->curcount[SEGMENT_PACK_IDX] = i;
                goto exit;
                break;
            default:
                fprintf(stderr, "*** ERROR in pack(%d)\n", rv);
                MPIR_Assert(0);
                break;
        }

        segp->in_offset[SEGMENT_PACK_IDX] += extent;
    }

  exit:
    *lastp = first + totalcopied;
    return;
}

void MPIR_Segment_unpack(DAME_Segment * segp,
                         DAME_Offset first, DAME_Offset * lastp, void *streambuf)
{
    /* TODO: TLP: Use max(size_t) */
    size_t inSize = LONG_MAX;
    if (*lastp != SEGMENT_IGNORE_LAST)
        inSize = (*lastp - first);

    DAME_Offset size, extent;
    DAME_Handle_get_size_macro(segp->handle, size);
    DAME_Handle_get_extent_macro(segp->handle, extent);

    DAME_Offset totalcopied = 0;

    DAME_Dame *dl;
    DAME_Handle_get_loopptr_macro(segp->handle, dl);

    if (!dl) {
        /* If no dataloop is defined, it should be a basic type for which
         * size == extent. Need to confirm this though */
        MPIR_Assert(size == extent);
        DAME_Dame defDL[3];
        defDL[0].kind = DL_EXIT;
        defDL[1].kind = DL_CONTIGFINAL;
        defDL[1].count = segp->count;
        defDL[1].size = size * segp->count;
        defDL[1].extent = extent;
        defDL[1].s.c_t.basesize = size;
        defDL[1].s.c_t.baseextent = extent;
        defDL[2].kind = DL_BOTTOM;

        DL_unpack(streambuf,
                  segp->ptr + segp->out_offset[SEGMENT_UNPACK_IDX],
                  inSize, defDL, 1, &segp->state[SEGMENT_UNPACK_IDX], &totalcopied);
        segp->curcount[SEGMENT_UNPACK_IDX] = totalcopied / size;
        segp->out_offset[SEGMENT_UNPACK_IDX] += (segp->curcount[SEGMENT_UNPACK_IDX] * extent);

        goto exit;
    }

    DAME_Count i;
    for (i = segp->curcount[SEGMENT_UNPACK_IDX]; i < segp->count; i++) {
        int rv;
        DAME_Offset copied = 0;
        DAME_Offset insz = inSize - totalcopied;
        rv = DL_unpack(streambuf + totalcopied,
                       segp->ptr + segp->out_offset[SEGMENT_UNPACK_IDX],
                       insz, dl, 0, &segp->state[SEGMENT_UNPACK_IDX], &copied);

        totalcopied += copied;
        /* TODO: TLP: Use the MPICH error handling */
        switch (rv) {
            case DAME_PACKUNPACK_COMPLETED:
                break;
            case DAME_PACKUNPACK_PARTIAL:
                segp->curcount[SEGMENT_UNPACK_IDX] = i;
                goto exit;
                break;
            default:
                fprintf(stderr, "*** ERROR in unpack(%d)\n", rv);
                MPIR_Assert(0);
                break;
        }

        segp->out_offset[SEGMENT_UNPACK_IDX] += extent;
    }
  exit:
    *lastp = first + totalcopied;
    return;
}
