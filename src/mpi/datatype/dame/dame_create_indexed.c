/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h>

#include "mpiimpl.h"

static void indexed_array_copy(DAME_Count count,
                               DAME_Count contig_count,
                               const DAME_Size * input_blocklength_array,
                               const void *input_displacement_array,
                               DAME_Size * output_blocklength_array,
                               DAME_Offset * out_disp_array,
                               int dispinbytes,
                               DAME_Offset old_extent, int *all_aligned, int *is_short);

/*@
  Dame_create_indexed

  Arguments:
  +  int icount
  .  DAME_Size *iblocklength_array
  .  void *displacement_array (either ints or MPI_Aints)
  .  int dispinbytes
  .  MPI_Datatype oldtype
  .  DAME_Dame **dlp_p
  .  int *dlsz_p
  .  int *dldepth_p
  -  int flag

  .N Errors
  .N Returns 0 on success, -1 on error.
  @*/

int MPIR_Dame_create_indexed(DAME_Count count,
                             const DAME_Size * blocklength_array,
                             const void *displacement_array,
                             int dispinbytes,
                             MPI_Datatype type,
                             MPI_Datatype oldtype, DAME_Dame ** dl, MPI_Aint * size, int *depth)
{
    if (!*dl) {
        MPIR_Assert(*depth == -1);
        MPIR_Dame_alloc(dl);
        *depth = 1;
    }

    int err, is_builtin;
    MPI_Aint i;
    DAME_Size blksz;
    DAME_Count first;
    int all_aligned = 1, is_short = 1;

    DAME_Count old_type_count = 0, contig_count;
    DAME_Offset old_size, old_extent;

    /* if count is zero, handle with contig code, call it an int */
    if (count == 0) {
        err = MPIR_Dame_create_contiguous(0, MPI_INT, dl, size, depth);
        return err;
    }

    /* Skip any initial zero-length blocks */
    for (first = 0; first < count; first++)
        if ((DAME_Count) blocklength_array[first])
            break;


    is_builtin = (DAME_Handle_hasloop_macro(oldtype)) ? 0 : 1;
    DAME_Handle_get_size_macro(oldtype, old_size);
    DAME_Handle_get_extent_macro(oldtype, old_extent);

    for (i = first; i < count; i++) {
        old_type_count += (DAME_Count) blocklength_array[i];
    }

    contig_count = MPIR_Type_indexed_count_contig(count,
                                                  blocklength_array,
                                                  displacement_array, dispinbytes, old_extent);

    /* if contig_count is zero (no data), handle with contig code */
    if (contig_count == 0) {
        err = MPIR_Dame_create_contiguous(0, MPI_INT, dl, size, depth);
        return err;
    }

    /* optimization:
     *
     * if contig_count == 1 and block starts at displacement 0,
     * store it as a contiguous rather than an indexed dataloop.
     */
    if ((contig_count == 1) &&
        ((!dispinbytes && ((int *) displacement_array)[first] == 0) ||
         (dispinbytes && ((MPI_Aint *) displacement_array)[first] == 0))) {
        err = MPIR_Dame_create_contiguous(old_type_count, oldtype, dl, size, depth);
        return err;
    }

    /* optimization:
     *
     * if contig_count == 1 (and displacement != 0), store this as
     * a single element blockindexed rather than a lot of individual
     * blocks.
     */
    if (contig_count == 1) {
        const void *disp_arr_tmp;       /* no ternary assignment to avoid clang warnings */
        if (dispinbytes)
            disp_arr_tmp = &(((const MPI_Aint *) displacement_array)[first]);
        else
            disp_arr_tmp = &(((const int *) displacement_array)[first]);
        err = MPIR_Dame_create_blockindexed(1,
                                            old_type_count,
                                            disp_arr_tmp,
                                            dispinbytes, type, oldtype, dl, size, depth);
        return err;
    }

    /* optimization:
     *
     * if block length is the same for all blocks, store it as a
     * blockindexed rather than an indexed dataloop.
     */
    blksz = blocklength_array[first];
    for (i = first + 1; i < count; i++) {
        if (blocklength_array[i] != blksz) {
            blksz--;
            break;
        }
    }
    if (blksz == blocklength_array[first]) {
        const void *disp_arr_tmp;       /* no ternary assignment to avoid clang warnings */
        if (dispinbytes)
            disp_arr_tmp = &(((const MPI_Aint *) displacement_array)[first]);
        else
            disp_arr_tmp = &(((const int *) displacement_array)[first]);
        err = MPIR_Dame_create_blockindexed(count - first,
                                            blksz,
                                            disp_arr_tmp,
                                            dispinbytes, type, oldtype, dl, size, depth);

        return err;
    }

    /* note: blockindexed looks for the vector optimization */

    /* TODO: optimization:
     *
     * if an indexed of a contig, absorb the contig into the blocklen array
     * and keep the same overall depth
     */

    /* otherwise storing as an indexed dataloop */

    MPI_Aint typesize, typeextent;
    DAME_Handle_get_size_macro(type, typesize);
    DAME_Handle_get_extent_macro(type, typeextent);

    MPIR_Datatype *oldtype_ptr;
    MPIR_Datatype_get_ptr(oldtype, oldtype_ptr);

    (*dl)[*depth].kind = DL_INDEX;
    (*dl)[*depth].count = contig_count;
    (*dl)[*depth].size = typesize;
    (*dl)[*depth].extent = typeextent;

    (*dl)[*depth].s.i_t.oldsize = old_size;
    /* copy in displacement parameters
     *
     * regardless of dispinbytes, we store displacements in bytes in loop.
     */
    DAME_Offset *tmpdisps = (DAME_Offset *) DAME_Malloc((contig_count) * sizeof(DAME_Offset));
    DAME_Size *tmpblks = (DAME_Size *) DAME_Malloc((contig_count) * sizeof(DAME_Size));

    indexed_array_copy(count,
                       contig_count,
                       blocklength_array,
                       displacement_array,
                       tmpblks, tmpdisps, dispinbytes, old_extent, &all_aligned, &is_short);
    (*dl)[*depth].s.i_t.offsets = tmpdisps;
    (*dl)[*depth].s.i_t.blklens = tmpblks;

    if (!is_builtin) {
        MPIR_Dame_create(oldtype,
                         &oldtype_ptr->dataloop,
                         &oldtype_ptr->dataloop_size, &oldtype_ptr->dataloop_depth);
    }

    if (is_builtin || Dame_is_contig(oldtype_ptr->dataloop)) {
        (*dl)[*depth].kind = DL_INDEXFINAL;
        if (all_aligned)
            DAME_opt_set_aligned((*dl)[*depth]);
        if (is_short)
            DAME_opt_set_isshort((*dl)[*depth]);

        *depth += 1;
        (*dl)[*depth].kind = DL_CONTIGFINAL;
        /* We don't use the size and extent for contigfinal in the interpreted
         * code, but we need this data in the JITted version */
        (*dl)[*depth].size = old_size;
        (*dl)[*depth].extent = old_extent;
        (*dl)[*depth].s.c_t.basesize = old_size;
        (*dl)[*depth].s.c_t.baseextent = old_extent;
        *depth += 1;
    } else {
        *depth += 1;
        (*dl)[*depth].kind = DL_CONTIGCHILD;
        (*dl)[*depth].size = old_size;
        (*dl)[*depth].extent = typeextent;
        int i;
        /* TODO: TLP: Use MPI_Aint max or something */
        MPI_Aint min = LONG_MAX;
        for (i = 0; i < contig_count; i++)
            if (tmpdisps[i] < min)
                min = tmpdisps[i];
        (*dl)[*depth].size = min;

        *depth += 1;
        /* The first element of the innerdl will be DL_EXIT */
        DAME_Dame *pos = &(*dl)[*depth];
        MPIR_Dame_dup(&oldtype_ptr->dataloop[1], 0, &pos);
        /* The last element will be DL_BOTTOM */
        *depth = *depth + oldtype_ptr->dataloop_depth - 2;
    }

    return MPI_SUCCESS;
}

/* indexed_array_copy()
 *
 * Copies arrays into place, combining adjacent contiguous regions and
 * dropping zero-length regions.
 *
 * Extent passed in is for the original type.
 *
 * Output displacements are always output in bytes, while block
 * lengths are always output in terms of the base type.
 */
static void indexed_array_copy(DAME_Count count,
                               DAME_Count contig_count,
                               const DAME_Size * in_blklen_array,
                               const void *in_disp_array,
                               DAME_Size * out_blklen_array,
                               DAME_Offset * out_disp_array,
                               int dispinbytes,
                               DAME_Offset old_extent, int *all_aligned, int *is_short)
{
    DAME_Count i, first, cur_idx = 0;
    *all_aligned = 1;
    *is_short = 1;

    /* Skip any initial zero-length blocks */
    for (first = 0; first < count; ++first)
        if ((DAME_Count) in_blklen_array[first])
            break;

    out_blklen_array[0] = (DAME_Count) in_blklen_array[first];

    if (!dispinbytes) {
        DAME_Offset offset0 = ((int *) in_disp_array)[first] * old_extent;
        out_disp_array[0] = offset0;
        switch (old_extent) {
            case 8:
            case 4:
            case 2:
                if (offset0 % old_extent != 0)
                    *all_aligned = 0;
                break;
            default:
                *all_aligned = 0;
                break;
        }

        for (i = first + 1; i < count; ++i) {
            if (in_blklen_array[i] == 0) {
                continue;
            } else if (out_disp_array[cur_idx] +
                       ((DAME_Offset) out_blklen_array[cur_idx]) * old_extent ==
                       ((DAME_Offset) ((int *) in_disp_array)[i]) * old_extent) {
                /* adjacent to current block; add to block */
                out_blklen_array[cur_idx] += (DAME_Count) in_blklen_array[i];
            } else {
                DAME_Offset offset = ((int *) in_disp_array)[i] * old_extent;
                cur_idx++;
                DAME_Assert(cur_idx < contig_count);
                out_disp_array[cur_idx] = offset;
                out_blklen_array[cur_idx] = in_blklen_array[i];
                if (in_blklen_array[i] > DAME_MEMCPY_THRESHOLD)
                    *is_short = 0;
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
    } else {    /* input displacements already in bytes */

        out_disp_array[0] = (DAME_Offset) ((MPI_Aint *) in_disp_array)[first];

        for (i = first + 1; i < count; ++i) {
            if (in_blklen_array[i] == 0) {
                continue;
            } else if (out_disp_array[cur_idx] +
                       ((DAME_Offset) out_blklen_array[cur_idx]) * old_extent ==
                       ((DAME_Offset) ((MPI_Aint *) in_disp_array)[i])) {
                /* adjacent to current block; add to block */
                out_blklen_array[cur_idx] += in_blklen_array[i];
            } else {
                DAME_Offset offset = ((MPI_Aint *) in_disp_array)[i];
                cur_idx++;
                DAME_Assert(cur_idx < contig_count);
                out_disp_array[cur_idx] = offset;
                out_blklen_array[cur_idx] = (DAME_Size) in_blklen_array[i];
                if (in_blklen_array[i] > DAME_MEMCPY_THRESHOLD)
                    *is_short = 0;
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
    }

    DAME_Assert(cur_idx == contig_count - 1);
    return;
}

/* indexed_count_contig()
 *
 * Determines the actual number of contiguous blocks represented by the
 * blocklength/displacement arrays.  This might be less than count (as
 * few as 1).
 *
 * Extent passed in is for the original type.
 */
DAME_Count MPIR_Type_indexed_count_contig(DAME_Count count,
                                          const DAME_Size * blocklength_array,
                                          const void *displacement_array,
                                          int dispinbytes, DAME_Offset old_extent)
{
    DAME_Count i, contig_count = 1;
    DAME_Count cur_blklen, first;

    if (count) {
        /* Skip any initial zero-length blocks */
        for (first = 0; first < count; ++first)
            if ((DAME_Count) blocklength_array[first])
                break;

        if (first == count) {   /* avoid invalid reads later on */
            contig_count = 0;
            return contig_count;
        }

        cur_blklen = (DAME_Count) blocklength_array[first];
        if (!dispinbytes) {
            DAME_Offset cur_tdisp = (DAME_Offset) ((int *) displacement_array)[first];

            for (i = first + 1; i < count; ++i) {
                if (blocklength_array[i] == 0) {
                    continue;
                } else if (cur_tdisp + (DAME_Offset) cur_blklen ==
                           (DAME_Offset) ((int *) displacement_array)[i]) {
                    /* adjacent to current block; add to block */
                    cur_blklen += (DAME_Count) blocklength_array[i];
                } else {
                    cur_tdisp = (DAME_Offset) ((int *) displacement_array)[i];
                    cur_blklen = (DAME_Count) blocklength_array[i];
                    contig_count++;
                }
            }
        } else {
            DAME_Offset cur_bdisp = (DAME_Offset) ((MPI_Aint *) displacement_array)[first];

            for (i = first + 1; i < count; ++i) {
                if (blocklength_array[i] == 0) {
                    continue;
                } else if (cur_bdisp + (DAME_Offset) cur_blklen * old_extent ==
                           (DAME_Offset) ((MPI_Aint *) displacement_array)[i]) {
                    /* adjacent to current block; add to block */
                    cur_blklen += (DAME_Count) blocklength_array[i];
                } else {
                    cur_bdisp = (DAME_Offset) ((MPI_Aint *) displacement_array)[i];
                    cur_blklen = (DAME_Count) blocklength_array[i];
                    contig_count++;
                }
            }
        }
    }
    return contig_count;
}
