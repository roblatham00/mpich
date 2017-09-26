/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>

#include "mpiimpl.h"

static void blockindexed_array_copy(DAME_Count count,
                                    const void *disp_array,
                                    DAME_Offset * out_disp_array,
                                    int dispinbytes, DAME_Offset old_extent, int *aligned);

/*@
  Dame_create_blockindexed - create blockindexed dataloop

  Arguments:
  +  DAME_Count count
  .  void *displacement_array (array of either MPI_Aints or ints)
  .  int displacement_in_bytes (boolean)
  .  MPI_Datatype old_type
  .  DAME_Dame **output_dataloop_ptr
  .  int output_dataloop_size
  .  int output_dataloop_depth
  -  int flag

  .N Errors
  .N Returns 0 on success, -1 on failure.
  @*/
int MPIR_Dame_create_blockindexed(DAME_Count count,
                                  DAME_Count blklen,
                                  const void *disp_array,
                                  int dispinbytes,
                                  MPI_Datatype type,
                                  DAME_Type oldtype,
                                  DAME_Dame ** dl, MPI_Aint * dlsize, int *dldepth)
{
    if (!*dl) {
        MPIR_Assert(*dldepth == -1);
        MPIR_Dame_alloc(dl);
        *dldepth = 1;
    }

    int err, is_builtin, is_vectorizable = 1;
    int i;

    DAME_Count contig_count;
    DAME_Offset old_size, old_extent, eff_disp0, eff_disp1, last_stride;
    int all_aligned = 1;

    /* if count or blklen are zero, handle with contig code, call it a int */
    if (count == 0 || blklen == 0) {
        err = MPIR_Dame_create_contiguous(0, MPI_INT, dl, dlsize, dldepth);
        return err;
    }

    is_builtin = (DAME_Handle_hasloop_macro(oldtype)) ? 0 : 1;

    DAME_Handle_get_size_macro(oldtype, old_size);
    DAME_Handle_get_extent_macro(oldtype, old_extent);


    contig_count = MPIR_Type_blockindexed_count_contig(count,
                                                       blklen, disp_array, dispinbytes, old_extent);

    /* optimization:
     *
     * if contig_count == 1 and block starts at displacement 0,
     * store it as a contiguous rather than a blockindexed dataloop.
     */
    if ((contig_count == 1) &&
        ((!dispinbytes && ((int *) disp_array)[0] == 0) ||
         (dispinbytes && ((MPI_Aint *) disp_array)[0] == 0))) {
        err = MPIR_Dame_create_contiguous(count * blklen, oldtype, dl, dlsize, dldepth);
        return err;
    }

    /* optimization:
     *
     * if contig_count == 1 store it as a blockindexed with one
     * element rather than as a lot of individual blocks.
     */
    if (contig_count == 1) {
        /* adjust count and blklen and drop through */
        blklen *= count;
        count = 1;
    }

    /* optimization:
     *
     * if displacements start at zero and result in a fixed stride,
     * store it as a vector rather than a blockindexed dataloop.
     */
    eff_disp0 = (dispinbytes) ? ((DAME_Offset) ((MPI_Aint *) disp_array)[0]) :
        (((DAME_Offset) ((int *) disp_array)[0]) * old_extent);

    if (count > 1 && eff_disp0 == (DAME_Offset) 0) {
        eff_disp1 = (dispinbytes) ?
            ((DAME_Offset) ((MPI_Aint *) disp_array)[1]) :
            (((DAME_Offset) ((int *) disp_array)[1]) * old_extent);
        last_stride = eff_disp1 - eff_disp0;

        for (i = 2; i < count; i++) {
            eff_disp0 = eff_disp1;
            eff_disp1 = (dispinbytes) ?
                ((DAME_Offset) ((MPI_Aint *) disp_array)[i]) :
                (((DAME_Offset) ((int *) disp_array)[i]) * old_extent);
            if (eff_disp1 - eff_disp0 != last_stride) {
                is_vectorizable = 0;
                break;
            }
        }
        if (is_vectorizable) {
            /* TODO: TLP: Write test to see if this optimization works
             * properly */
            err = MPIR_Dame_create_vector(count, blklen, last_stride, 1,        /* strideinbytes */
                                          type, oldtype, dl, dlsize, dldepth);
            return err;
        }
    }

    /* TODO: optimization:
     *
     * if displacements result in a fixed stride, but first displacement
     * is not zero, store it as a blockindexed (blklen == 1) of a vector.
     */

    /* TODO: optimization:
     *
     * if a blockindexed of a contig, absorb the contig into the blocklen
     * parameter and keep the same overall depth
     */

    /* TODO: optimization:
     *
     * Coalesce contig blocks
     */

    /* otherwise storing as a blockindexed dataloop */

    /* Q: HOW CAN WE TELL IF IT IS WORTH IT TO STORE AS AN
     * INDEXED WITH FEWER CONTIG BLOCKS (IF CONTIG_COUNT IS SMALL)?
     */
    MPI_Aint typesize, typeextent;
    DAME_Handle_get_size_macro(type, typesize);
    DAME_Handle_get_extent_macro(type, typeextent);

    MPIR_Datatype *oldtype_ptr, *type_ptr;
    MPIR_Datatype_get_ptr(oldtype, oldtype_ptr);
    MPIR_Datatype_get_ptr(type, type_ptr);

    (*dl)[*dldepth].kind = DL_BLOCKINDEX;

    (*dl)[*dldepth].count = count;
    (*dl)[*dldepth].size = typesize;
    (*dl)[*dldepth].extent = typeextent;

    (*dl)[*dldepth].s.bi_t.blklen = blklen;
    (*dl)[*dldepth].s.bi_t.oldsize = old_size;
    /* copy in displacement parameters
     *
     * regardless of dispinbytes, we store displacements in bytes in loop.
     */
    DAME_Offset *tmpdisps = (DAME_Offset *) DAME_Malloc(count * sizeof(DAME_Offset));
    blockindexed_array_copy(count, disp_array, tmpdisps, dispinbytes, old_extent, &all_aligned);
    (*dl)[*dldepth].s.bi_t.offsets = tmpdisps;

    if (!is_builtin) {
        MPIR_Dame_create(oldtype,
                         &oldtype_ptr->dataloop,
                         &oldtype_ptr->dataloop_size, &oldtype_ptr->dataloop_depth);
    }

    if (is_builtin || Dame_is_contig(oldtype_ptr->dataloop)) {
        (*dl)[*dldepth].kind = DL_BLOCKINDEXFINAL;

        if (all_aligned)
            DAME_opt_set_aligned((*dl)[*dldepth]);
        if (blklen < DAME_MEMCPY_THRESHOLD)
            DAME_opt_set_isshort((*dl)[*dldepth]);

        *dldepth += 1;
        (*dl)[*dldepth].kind = DL_CONTIGFINAL;
        (*dl)[*dldepth].size = blklen * old_size;
        (*dl)[*dldepth].extent = blklen * old_extent;
        (*dl)[*dldepth].s.c_t.basesize = old_size;
        (*dl)[*dldepth].s.c_t.baseextent = old_extent;
        *dldepth += 1;
    } else {
        if (blklen > 1) {
            *dldepth += 1;
            (*dl)[*dldepth].kind = DL_CONTIGCHILD;
            (*dl)[*dldepth].size = old_size;
            (*dl)[*dldepth].extent = typeextent;
            int i;
            /* TODO: TLP: Use MPI_Aint max or something */
            MPI_Aint min = LONG_MAX;
            for (i = 0; i < count; i++)
                if (tmpdisps[i] < min)
                    min = tmpdisps[i];
            (*dl)[*dldepth].size = min;
        } else {
            (*dl)[*dldepth].kind = DL_BLOCKINDEX1;
        }

        *dldepth += 1;
        /* The first element of the innerdl will be DL_EXIT */
        DAME_Dame *pos = &(*dl)[*dldepth];
        MPIR_Dame_dup(&oldtype_ptr->dataloop[1], 0, &pos);
        /* The last element will be DL_BOTTOM */
        *dldepth = *dldepth + oldtype_ptr->dataloop_depth - 2;
    }

    return 0;
}

/* DAME_Type_blockindexed_array_copy
 *
 * Unlike the indexed version, this one does not compact adjacent
 * blocks, because that would really mess up the blockindexed type!
 */
static void blockindexed_array_copy(DAME_Count count,
                                    const void *in_disp_array,
                                    DAME_Offset * out_disp_array,
                                    int dispinbytes, DAME_Offset old_extent, int *all_aligned)
{
    int i;
    *all_aligned = 1;
    if (!dispinbytes) {
        for (i = 0; i < count; i++) {
            DAME_Offset offset = ((int *) in_disp_array)[i] * old_extent;
            out_disp_array[i] = offset;
            switch (old_extent) {
                case 8:
                case 4:
                case 2:
                    if (offset % old_extent != 0)
                        *all_aligned = 0;
                    break;
                default:
                    *all_aligned = 0;
                    break;
            }
        }
    } else {
        for (i = 0; i < count; i++) {
            DAME_Offset offset = ((MPI_Aint *) in_disp_array)[i];
            out_disp_array[i] = offset;
            switch (old_extent) {
                case 8:
                case 4:
                case 2:
                    if (offset % old_extent != 0)
                        *all_aligned = 0;
                    break;
                default:
                    *all_aligned = 0;
                    break;
            }
        }
    }
    return;
}

DAME_Count MPIR_Type_blockindexed_count_contig(DAME_Count count,
                                               DAME_Count blklen,
                                               const void *disp_array,
                                               int dispinbytes, DAME_Offset old_extent)
{
    int i, contig_count = 1;

    if (!dispinbytes) {
        /* this is from the MPI type, is of type int */
        DAME_Offset cur_tdisp = (DAME_Offset) ((int *) disp_array)[0];

        for (i = 1; i < count; i++) {
            DAME_Offset next_tdisp = (DAME_Offset) ((int *) disp_array)[i];

            if (cur_tdisp + (DAME_Offset) blklen != next_tdisp) {
                contig_count++;
            }
            cur_tdisp = next_tdisp;
        }
    } else {
        /* this is from the MPI type, is of type MPI_Aint */
        DAME_Offset cur_bdisp = (DAME_Offset) ((MPI_Aint *) disp_array)[0];

        for (i = 1; i < count; i++) {
            DAME_Offset next_bdisp = (DAME_Offset) ((MPI_Aint *) disp_array)[i];

            if (cur_bdisp + (DAME_Offset) blklen * old_extent != next_bdisp) {
                contig_count++;
            }
            cur_bdisp = next_bdisp;
        }
    }
    return contig_count;
}
