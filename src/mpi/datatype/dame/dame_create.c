/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h>
#include <limits.h>

#include "mpiimpl.h"

static void Dame_create_named(MPI_Datatype type, DAME_Dame ** dl, MPI_Aint * size, int *depth);

void MPIR_Dame_create(MPI_Datatype type, DAME_Dame ** dl, MPI_Aint * size, int *depth)
{
    int i;
    int err ATTRIBUTE((unused));

    int nr_ints, nr_aints, nr_types, combiner;
    MPI_Datatype *types;
    int *ints;
    MPI_Aint *aints;

    DAME_Dame *old_dlp;
    MPI_Aint old_dlsize;
    int old_dldepth;

    int dummy1, dummy2, dummy3, type0_combiner, ndims;
    MPI_Datatype tmptype;

    MPI_Aint stride;
    MPI_Aint *disps;
    DAME_Size *blklen;

    MPIR_Type_get_envelope(type, &nr_ints, &nr_aints, &nr_types, &combiner);

    /* some named types do need dataloops; handle separately. */
    if (combiner == MPI_COMBINER_NAMED) {
        Dame_create_named(type, dl, size, depth);
        return;
    } else if (combiner == MPI_COMBINER_F90_REAL ||
               combiner == MPI_COMBINER_F90_COMPLEX || combiner == MPI_COMBINER_F90_INTEGER) {
        MPI_Datatype f90basetype;
        DAME_Handle_get_basic_type_macro(type, f90basetype);
        MPIR_Dame_create_contiguous(1, f90basetype, dl, size, depth);
        return;
    }

    DAME_Handle_get_loopptr_macro(type, old_dlp);
    DAME_Handle_get_loopsize_macro(type, old_dlsize);
    DAME_Handle_get_loopdepth_macro(type, old_dldepth);
    if (old_dlp != NULL) {
        /* dataloop already created; just return it. */
        *dl = old_dlp;
        *depth = old_dldepth;
        return;
    }

    if (!*dl) {
        MPIR_Dame_alloc(dl);
        if (*depth == -1)
            *depth = 1;
    }

    MPIR_Type_access_contents(type, &ints, &aints, &types);

    /* first check for zero count on types where that makes sense */
    switch (combiner) {
        case MPI_COMBINER_CONTIGUOUS:
        case MPI_COMBINER_VECTOR:
        case MPI_COMBINER_HVECTOR_INTEGER:
        case MPI_COMBINER_HVECTOR:
        case MPI_COMBINER_INDEXED_BLOCK:
        case MPI_COMBINER_HINDEXED_BLOCK:
        case MPI_COMBINER_INDEXED:
        case MPI_COMBINER_HINDEXED_INTEGER:
        case MPI_COMBINER_HINDEXED:
        case MPI_COMBINER_STRUCT_INTEGER:
        case MPI_COMBINER_STRUCT:
            if (ints[0] == 0) {
                MPIR_Dame_create_contiguous(0, MPI_INT, dl, size, depth);
                goto clean_exit;
            }
            break;
        default:
            break;
    }

    MPIR_Type_get_envelope(types[0], &dummy1, &dummy2, &dummy3, &type0_combiner);
    if (type0_combiner != MPI_COMBINER_NAMED) {
        DAME_Handle_get_loopptr_macro(types[0], old_dlp);
        if (old_dlp == NULL) {
            /* no dataloop already present; create and store one */
            MPIR_Dame_create(types[0], &old_dlp, &old_dlsize, &old_dldepth);

            DAME_Handle_set_loopptr_macro(types[0], old_dlp);
            DAME_Handle_set_loopdepth_macro(types[0], old_dlsize);
            DAME_Handle_set_loopdepth_macro(types[0], old_dldepth);
        }
    }


    switch (combiner) {
        case MPI_COMBINER_DUP:
            if (type0_combiner != MPI_COMBINER_NAMED) {
                MPIR_Dame_dup(old_dlp, 0, dl);
                *depth = old_dldepth;
            } else {
                MPIR_Dame_create_contiguous(1, types[0], dl, size, depth);
            }
            break;
        case MPI_COMBINER_RESIZED:
            if (type0_combiner != MPI_COMBINER_NAMED) {
                MPIR_Dame_dup(old_dlp, 0, dl);
                *depth = old_dldepth - 1;
            } else {
                MPIR_Dame_create_contiguous(1, types[0], dl, size, depth);
            }
            /* TODO: TLP: FIXME: I think this is correct, the old
             * code seemed to do something of this sort as well. Plus, there
             * doesn't seem to be another clean way of getting the correct
             * extent in place */
            (*dl)[1].extent = aints[1];
            break;
        case MPI_COMBINER_CONTIGUOUS:
            MPIR_Dame_create_contiguous(ints[0] /* count */ ,
                                        types[0] /* oldtype */ ,
                                        dl, size, depth);
            break;
        case MPI_COMBINER_VECTOR:
            MPIR_Dame_create_vector(ints[0] /* count */ ,
                                    ints[1] /* blklen */ ,
                                    ints[2] /* stride */ ,
                                    0 /* stride not bytes */ ,
                                    type, types[0] /* oldtype */ ,
                                    dl, size, depth);
            break;
        case MPI_COMBINER_HVECTOR_INTEGER:
        case MPI_COMBINER_HVECTOR:
            /* fortran hvector has integer stride in bytes */
            if (combiner == MPI_COMBINER_HVECTOR_INTEGER) {
                stride = (MPI_Aint) ints[2];
            } else {
                stride = aints[0];
            }
            MPIR_Dame_create_vector(ints[0] /* count */ ,
                                    ints[1] /* blklen */ ,
                                    stride, 1 /* stride in bytes */ ,
                                    type, types[0] /* oldtype */ ,
                                    dl, size, depth);
            break;
        case MPI_COMBINER_INDEXED_BLOCK:
            MPIR_Dame_create_blockindexed(ints[0] /* count */ ,
                                          ints[1] /* blklen */ ,
                                          &ints[2] /* disps */ ,
                                          0 /* disp not bytes */ ,
                                          type, types[0] /* oldtype */ ,
                                          dl, size, depth);
            break;
        case MPI_COMBINER_HINDEXED_BLOCK:
            disps = (MPI_Aint *) DAME_Malloc(ints[0] * sizeof(MPI_Aint));
            for (i = 0; i < ints[0]; i++)
                disps[i] = aints[i];
            MPIR_Dame_create_blockindexed(ints[0] /* count */ ,
                                          ints[1] /* blklen */ ,
                                          disps /* disps */ ,
                                          1 /* disp is bytes */ ,
                                          type, types[0] /* oldtype */ ,
                                          dl, size, depth);
            DAME_Free(disps);
            break;
        case MPI_COMBINER_INDEXED:
            blklen = (DAME_Size *) DAME_Malloc(ints[0] * sizeof(DAME_Size));
            for (i = 0; i < ints[0]; i++)
                blklen[i] = ints[1 + i];
            MPIR_Dame_create_indexed(ints[0] /* count */ ,
                                     blklen /* blklens */ ,
                                     &ints[ints[0] + 1] /* disp */ ,
                                     0 /* disp not in bytes */ ,
                                     type, types[0] /* oldtype */ ,
                                     dl, size, depth);
            DAME_Free(blklen);
            break;
        case MPI_COMBINER_HINDEXED_INTEGER:
        case MPI_COMBINER_HINDEXED:
            if (combiner == MPI_COMBINER_HINDEXED_INTEGER) {
                disps = (MPI_Aint *) DAME_Malloc(ints[0] * sizeof(MPI_Aint));

                for (i = 0; i < ints[0]; i++) {
                    disps[i] = (MPI_Aint) ints[ints[0] + 1 + i];
                }
            } else {
                disps = aints;
            }

            blklen = (DAME_Size *) DAME_Malloc(ints[0] * sizeof(DAME_Size));
            for (i = 0; i < ints[0]; i++)
                blklen[i] = (DAME_Size) ints[1 + i];
            MPIR_Dame_create_indexed(ints[0] /* count */ ,
                                     blklen /* blklens */ ,
                                     disps, 1 /* disp in bytes */ ,
                                     type, types[0] /* oldtype */ ,
                                     dl, size, depth);

            if (combiner == MPI_COMBINER_HINDEXED_INTEGER) {
                DAME_Free(disps);
            }
            DAME_Free(blklen);

            break;
        case MPI_COMBINER_STRUCT_INTEGER:
        case MPI_COMBINER_STRUCT:
            if (combiner == MPI_COMBINER_STRUCT_INTEGER) {
                disps = (MPI_Aint *) DAME_Malloc(ints[0] * sizeof(MPI_Aint));

                for (i = 0; i < ints[0]; i++) {
                    disps[i] = (MPI_Aint) ints[ints[0] + 1 + i];
                }
            } else {
                disps = aints;
            }

            err = MPIR_Dame_create_struct(ints[0] /* count */ ,
                                          &ints[1] /* blklens */ ,
                                          disps, type, types /* oldtype array */ ,
                                          dl, size, depth);
            /* TODO if/when this function returns error codes, propagate this failure instead */
            DAME_Assert(0 == err);
            /* if (err) return err; */

            if (combiner == MPI_COMBINER_STRUCT_INTEGER) {
                DAME_Free(disps);
            }
            break;
        case MPI_COMBINER_SUBARRAY:
            ndims = ints[0];
            MPIR_Type_convert_subarray(ndims, &ints[1] /* sizes */ ,
                                       &ints[1 + ndims] /* subsizes */ ,
                                       &ints[1 + 2 * ndims] /* starts */ ,
                                       ints[1 + 3 * ndims] /* order */ ,
                                       types[0], &tmptype);

            MPIR_Dame_create(tmptype, dl, size, depth);

            MPIR_Type_free_impl(&tmptype);
            break;
        case MPI_COMBINER_DARRAY:
            ndims = ints[2];
            MPIR_Type_convert_darray(ints[0] /* size */ ,
                                     ints[1] /* rank */ ,
                                     ndims, &ints[3] /* gsizes */ ,
                                     &ints[3 + ndims] /*distribs */ ,
                                     &ints[3 + 2 * ndims] /* dargs */ ,
                                     &ints[3 + 3 * ndims] /* psizes */ ,
                                     ints[3 + 4 * ndims] /* order */ ,
                                     types[0], &tmptype);

            MPIR_Dame_create(tmptype, dl, size, depth);

            MPIR_Type_free_impl(&tmptype);
            break;
        default:
            DAME_Assert(0);
            break;
    }

  clean_exit:

    *depth += 1;
    MPIR_Type_release_contents(type, &ints, &aints, &types);

    /* for now we just leave the intermediate dataloops in place.
     * could remove them to save space if we wanted.
     */

    return;
}

/*@
  Dame_create_named - create a dataloop for a "named" type
  if necessary.

  "named" types are ones for which MPI_Type_get_envelope() returns a
  combiner of MPI_COMBINER_NAMED. some types that fit this category,
  such as MPI_SHORT_INT, have multiple elements with potential gaps
  and padding. these types need dataloops for correct processing.
  @*/
static void Dame_create_named(MPI_Datatype type, DAME_Dame ** dl, MPI_Aint * size, int *depth)
{
    /* special case: pairtypes need dataloops too.
     *
     * note: not dealing with MPI_2INT because size == extent
     *       in all cases for that type.
     *
     * note: MPICH always precreates these, so we will never call
     *       Dame_create_pairtype() from here in the MPICH
     *       case.
     */
    if (type == MPI_FLOAT_INT || type == MPI_DOUBLE_INT ||
        type == MPI_LONG_INT || type == MPI_SHORT_INT || type == MPI_LONG_DOUBLE_INT) {
        DAME_Dame *dlp;
        int dldepth;
        DAME_Handle_get_loopptr_macro(type, dlp);
        DAME_Handle_get_loopdepth_macro(type, dldepth);
        if (dlp != NULL) {
            /* dataloop already created; just return it. */
            *dl = dlp;
            *depth = dldepth;
        } else {
            MPIR_Dame_create_pairtype(type, dl, size, depth);
        }
        return;
    }
    /* no other combiners need dataloops; exit. */
    else {
        *dl = NULL;
        return;
    }
}
