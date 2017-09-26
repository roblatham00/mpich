/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
#ifndef TYPESIZE_SUPPORT_H
#define TYPESIZE_SUPPORT_H

#include "./dame.h"

#define DAME_Type_footprint MPIR_Type_footprint

typedef struct MPIR_Type_footprint_s {
    DAME_Offset size, extent;

    /* these are only needed for calculating footprint of types
     * built using this type. no reason to expose these.
     */
    DAME_Offset lb, ub, alignsz;
    DAME_Offset true_lb, true_ub;
    int has_sticky_lb;
    int has_sticky_ub;
} DAME_Type_footprint;

void MPIR_Type_calc_footprint(MPI_Datatype type, DAME_Type_footprint * tfp);

#endif
