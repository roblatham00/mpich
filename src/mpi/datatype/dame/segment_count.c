/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "mpiimpl.h"


/* DAME_Segment_count_contig_blocks()
 *
 * Count number of contiguous regions in segment between first and last.
 */
void MPIR_Segment_count_contig_blocks(DAME_Segment * segp,
                                      DAME_Offset first, DAME_Offset * lastp, DAME_Count * countp)
{
    fprintf(stderr, "*** Unimplemented: Segment_count_contig_blocks\n");
    MPIR_Assert(0);
    return;
}

/*
 * Local variables:
 * c-indent-tabs-mode: nil
 * End:
 */
