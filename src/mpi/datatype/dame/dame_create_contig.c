/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpiimpl.h"

/*@
  Dame_contiguous - create the dataloop representation for a
  contiguous datatype

  Arguments:
  +  int icount,
  .  MPI_Datatype oldtype,
  .  DAME_Dame **dlp_p,
  .  int *dlsz_p,
  .  int *dldepth_p,
  -  int flag

  .N Errors
  .N Returns 0 on success, -1 on failure.
  @*/
int MPIR_Dame_create_contiguous(DAME_Count count,
                                DAME_Type oldtype, DAME_Dame ** dl, MPI_Aint * size, int *depth)
{
    if (!*dl) {
        MPIR_Assert(*depth == -1);
        MPIR_Dame_alloc(dl);
        *depth = 1;
    }

    int is_builtin;
    is_builtin = (DAME_Handle_hasloop_macro(oldtype)) ? 0 : 1;

    DAME_Size oldsize = 0, oldextent = 0;
    DAME_Handle_get_size_macro(oldtype, oldsize);
    DAME_Handle_get_extent_macro(oldtype, oldextent);

    (*dl)[*depth].kind = DL_CONTIG;
    (*dl)[*depth].count = count;
    (*dl)[*depth].size = count * oldsize;
    (*dl)[*depth].extent = count * oldextent;
    (*dl)[*depth].s.c_t.basesize = oldsize;
    (*dl)[*depth].s.c_t.baseextent = oldextent;

    if (is_builtin) {
        (*dl)[*depth].kind = DL_CONTIGFINAL;
        (*dl)[*depth].size = count * oldsize;
        (*dl)[*depth].extent = count * oldextent;
        *depth += 1;
    } else {
        MPIR_Datatype *old_ptr;
        MPIR_Datatype_get_ptr(oldtype, old_ptr);

        MPIR_Dame_create(oldtype,
                         &old_ptr->dataloop, &old_ptr->dataloop_size, &old_ptr->dataloop_depth);

        if ((old_ptr->dataloop[1].kind == DL_CONTIGFINAL)
            && (oldsize == oldextent)) {
            (*dl)[*depth].kind = DL_CONTIGFINAL;
            (*dl)[*depth].size = count * oldsize;
            (*dl)[*depth].extent = count * oldextent;
            *depth += 1;
        } else {
            if (count > 1) {
                *depth += 1;
                DAME_Dame *pos = &(*dl)[*depth];
                MPIR_Dame_dup(&old_ptr->dataloop[1], 0, &pos);
                /* We have to ignore EXIT and BOTTOM */
                *depth = *depth + old_ptr->dataloop_depth - 2;
            } else {
                DAME_Dame *pos = &(*dl)[*depth];
                MPIR_Dame_dup(&old_ptr->dataloop[1], 0, &pos);
                *depth = *depth + old_ptr->dataloop_depth;
            }
        }
    }

    (*dl)[*depth].kind = DL_BOTTOM;

    return 0;
}
