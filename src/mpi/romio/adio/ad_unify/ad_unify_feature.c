/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "adio.h"
#include "ad_unify.h"

int ADIOI_UNIFY_Feature(ADIO_File fd, int flag)
{
    switch (flag) {
        case ADIO_SCALABLE_RESIZE:
        case ADIO_TWO_PHASE:
            return 1;
        case ADIO_LOCKS:
        case ADIO_SHARED_FP:
        case ADIO_ATOMIC_MODE:
        case ADIO_DATA_SIEVING_WRITES:
        case ADIO_SCALABLE_OPEN:
        case ADIO_SEQUENTIAL:
        default:
            return 0;
    }
}
