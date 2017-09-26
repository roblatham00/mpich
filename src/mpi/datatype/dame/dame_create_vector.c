/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpiimpl.h"

/*@
  Dame_create_vector

  Arguments:
  +  int icount
  .  int iblocklength
  .  MPI_Aint astride
  .  int strideinbytes
  .  MPI_Datatype oldtype
  .  DAME_Dame **dlp_p
  .  int *dlsz_p
  .  int *dldepth_p
  -  int flag

  Returns 0 on success, -1 on failure.

  @*/
int MPIR_Dame_create_vector(DAME_Count count,
                            DAME_Size blocklength,
                            MPI_Aint stride,
                            int strideinbytes,
                            MPI_Datatype type,
                            DAME_Type oldtype, DAME_Dame ** dl, MPI_Aint * size, int *depth)
{
    if (!*dl) {
        MPIR_Assert(*depth == -1);
        MPIR_Dame_alloc(dl);
        *depth = 1;
    }

    /* if count or blocklength are zero, handle with contig code,
     * call it a int
     */
    if (count == 0 || blocklength == 0) {
        int err;
        err = MPIR_Dame_create_contiguous(0, MPI_INT, dl, size, depth);
        return err;
    }

    /* optimization:
     *
     * if count == 1, store as a contiguous rather than a vector dataloop.
     */
    if (count == 1) {
        int err;
        err = MPIR_Dame_create_contiguous(blocklength, oldtype, dl, size, depth);
        return err;
    }

    int is_builtin;
    DAME_Size oldsize = 0, oldextent = 0;
    DAME_Handle_get_size_macro(oldtype, oldsize);
    DAME_Handle_get_extent_macro(oldtype, oldextent);

    MPIR_Datatype *oldtype_ptr;
    MPIR_Datatype_get_ptr(oldtype, oldtype_ptr);

    MPI_Aint typesize, typeextent;
    DAME_Handle_get_size_macro(type, typesize);
    DAME_Handle_get_extent_macro(type, typeextent);

    is_builtin = (DAME_Handle_hasloop_macro(oldtype)) ? 0 : 1;

    (*dl)[*depth].kind = DL_VECTOR;
    (*dl)[*depth].count = count;
    (*dl)[*depth].size = typesize;
    (*dl)[*depth].extent = typeextent;
    (*dl)[*depth].s.v_t.stride = strideinbytes ? stride : stride * oldextent;
    (*dl)[*depth].s.v_t.blklen = blocklength;
    (*dl)[*depth].s.v_t.oldsize = oldsize;

    if (!is_builtin) {
        MPIR_Dame_create(oldtype,
                         &oldtype_ptr->dataloop,
                         &oldtype_ptr->dataloop_size, &oldtype_ptr->dataloop_depth);
    }

    if (is_builtin || Dame_is_contig(oldtype_ptr->dataloop)) {
        (*dl)[*depth].kind = DL_VECTORFINAL;
        if (stride % oldsize == 0)
            DAME_opt_set_aligned((*dl)[*depth]);
        if (blocklength < DAME_MEMCPY_THRESHOLD)
            DAME_opt_set_isshort((*dl)[*depth]);

        *depth += 1;
        (*dl)[*depth].kind = DL_CONTIGFINAL;
        (*dl)[*depth].size = blocklength * oldsize;
        (*dl)[*depth].extent = blocklength * oldextent;
        (*dl)[*depth].s.c_t.basesize = oldsize;
        (*dl)[*depth].s.c_t.baseextent = oldextent;
        *depth += 1;
    } else {
        if (blocklength > 1) {
            *depth += 1;
            (*dl)[*depth].kind = DL_CONTIGCHILD;
            (*dl)[*depth].size = oldsize;
            (*dl)[*depth].extent = typeextent;
        } else {
            (*dl)[*depth].kind = DL_VECTOR1;
        }

        *depth += 1;
        /* The first element of the innerdl will be DL_EXIT */
        DAME_Dame *pos = &(*dl)[*depth];
        MPIR_Dame_dup(&oldtype_ptr->dataloop[1], 0, &pos);
        /* The last element will be DL_BOTTOM */
        *depth = *depth + oldtype_ptr->dataloop_depth - 2;
    }

    return 0;
}
