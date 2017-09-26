/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpiimpl.h"

static int Dame_create_struct_memory_error(void);
static int Dame_create_unique_type_struct(DAME_Count count,
                                          const int *blklens,
                                          const MPI_Aint * disps,
                                          MPI_Datatype type,
                                          const DAME_Type * oldtypes,
                                          int type_pos,
                                          DAME_Dame ** dl, MPI_Aint * size, int *depth);
static int Dame_create_contig_all_bytes_struct(DAME_Count count,
                                               const int *blklens,
                                               const MPI_Aint * disps,
                                               MPI_Datatype type,
                                               const DAME_Type * oldtypes,
                                               DAME_Dame ** dl, MPI_Aint * dlsize, int *depth);

/*@
  Dame_create_struct - create the dataloop representation for a
  struct datatype

  Input Parameters:
  + count - number of blocks in vector
  . blklens - number of elements in each block
  . disps - offsets of blocks from start of type in bytes
  - oldtypes - types (using handle) of datatypes on which vector is based

  Output Parameters:
  + dlp_p - pointer to address in which to place pointer to new dataloop
  - dlsz_p - pointer to address in which to place size of new dataloop

  Return Value:
  0 on success, -1 on failure.

  Notes:
  This function relies on others, like Dame_create_indexed, to create
  types in some cases. This call (like all the rest) takes int blklens
  and MPI_Aint displacements, so it's possible to overflow when working
  with a particularly large struct type in some cases. This isn't detected
  or corrected in this code at this time.

  @*/
int MPIR_Dame_create_struct(DAME_Count count,
                            const int *blklens,
                            const MPI_Aint * disps,
                            MPI_Datatype type,
                            const DAME_Type * oldtypes,
                            DAME_Dame ** dl, MPI_Aint * dlsize, int *depth)
{
    if (!*dl) {
        MPIR_Assert(*depth == 0);
        MPIR_Dame_alloc(dl);
        *depth += 1;
    }

    int err, i, nr_basics = 0, nr_derived = 0, type_pos = 0, nr_contigs = 0;

    DAME_Type first_basic = MPI_DATATYPE_NULL, first_derived = MPI_DATATYPE_NULL;

    /* variables used in general case only */
    int loop_idx;

    /* if count is zero, handle with contig code, call it a int */
    if (count == 0) {
        err = MPIR_Dame_create_contiguous(0, MPI_INT, dl, dlsize, depth);
        return err;
    }

    /* browse the old types and characterize */
    for (i = 0; i < count; i++) {
        /* ignore type elements with a zero blklen */
        if (blklens[i] == 0)
            continue;

        if (oldtypes[i] != MPI_LB && oldtypes[i] != MPI_UB) {
            int is_builtin;

            is_builtin = (DAME_Handle_hasloop_macro(oldtypes[i])) ? 0 : 1;

            if (is_builtin) {
                if (nr_basics == 0) {
                    first_basic = oldtypes[i];
                    type_pos = i;
                } else if (oldtypes[i] != first_basic) {
                    first_basic = MPI_DATATYPE_NULL;
                }
                nr_basics++;
                nr_contigs++;
            } else {    /* derived type */

                if (nr_derived == 0) {
                    first_derived = oldtypes[i];
                    type_pos = i;
                } else if (oldtypes[i] != first_derived) {
                    first_derived = MPI_DATATYPE_NULL;
                }
                MPIR_Datatype *old_ptr;
                MPIR_Datatype_get_ptr(oldtypes[i], old_ptr);
                if (old_ptr->is_contig)
                    nr_contigs++;
                nr_derived++;
            }
        }
    }

    /* note on optimizations:
     *
     * because LB, UB, and extent calculations are handled as part of
     * the Datatype, we can safely ignore them in all our calculations
     * here.
     */

    /* optimization:
     *
     * if there were only MPI_LBs and MPI_UBs in the struct type,
     * treat it as a zero-element contiguous (just as count == 0).
     */
    if (nr_basics == 0 && nr_derived == 0) {
        err = MPIR_Dame_create_contiguous(0, MPI_INT, dl, dlsize, depth);
        return err;
    }

    /* optimization:
     *
     * If all the subtypes are contigs, create an hindexed type
     */
    if (nr_contigs == count) {
        err = Dame_create_contig_all_bytes_struct(count,
                                                  blklens,
                                                  disps, type, oldtypes, dl, dlsize, depth);
        return err;
    }

    MPI_Aint size = 0, extent = 0;
    DAME_Handle_get_size_macro(type, size);
    DAME_Handle_get_extent_macro(type, extent);

    /* optimization:
     *
     * if there is only one unique instance of a type in the struct, treat it
     * as a blockindexed type.
     *
     * notes:
     *
     * if the displacement happens to be zero, the blockindexed code will
     * optimize this into a contig.
     */
    if (nr_basics + nr_derived == 1) {
        /* type_pos is index to only real type in array */
        err = MPIR_Dame_create_blockindexed(1,  /* count */
                                            blklens[type_pos], &disps[type_pos], 1,     /* displacement in bytes */
                                            type, oldtypes[type_pos], dl, dlsize, depth);

        return err;
    }

    /* optimization:
     *
     * if there only one unique type (more than one instance) in the
     * struct, treat it as an indexed type.
     *
     * notes:
     *
     * this will apply to a single type with an LB/UB, as those
     * are handled elsewhere.
     *
     */
    if (((nr_derived == 0) && (first_basic != MPI_DATATYPE_NULL)) ||
        ((nr_basics == 0) && (first_derived != MPI_DATATYPE_NULL))) {
        return Dame_create_unique_type_struct(count,
                                              blklens,
                                              disps, type, oldtypes, type_pos, dl, dlsize, depth);
    }

    MPIR_Dame_struct_alloc((DAME_Count) nr_basics + nr_derived, &(*dl)[*depth]);

    MPI_Aint bufsize = (nr_basics + nr_derived) * sizeof(MPI_Aint);
    MPI_Aint *tmpoffsets = (MPI_Aint *) MPL_malloc(bufsize);
    MPI_Aint *tmpblklens = (MPI_Aint *) MPL_malloc(bufsize);
    MPI_Aint *tmpoldsizes = (MPI_Aint *) MPL_malloc(bufsize);

    int d = *depth;
    for (i = 0, loop_idx = 0; i < count; i++) {
        *depth = d;
        MPI_Aint inner_dlsize = 0;
        MPI_Aint oldsize = 0, oldextent = 0;
        DAME_Handle_get_size_macro(oldtypes[i], oldsize);
        DAME_Handle_get_extent_macro(oldtypes[i], oldsize);
        int is_builtin;

        /* ignore type elements with a zero blklen */
        if (blklens[i] == 0)
            continue;

        is_builtin = (DAME_Handle_hasloop_macro(oldtypes[i])) ? 0 : 1;

        if (is_builtin) {
            /* LBs and UBs already taken care of -- skip them */
            if (oldtypes[i] == MPI_LB || oldtypes[i] == MPI_UB) {
                continue;
            }

            /* build a contig dataloop for this basic and point to that
             *
             * optimization:
             *
             * push the count (blklen) from the struct down into the
             * contig so we can process more at the leaf.
             */
            *depth += 1;
            err = MPIR_Dame_create_contiguous(blklens[i],
                                              oldtypes[i],
                                              &(*dl)[d].s.s_t.dls[loop_idx], &inner_dlsize, depth);

            /* --BEGIN ERROR HANDLING-- */
            if (err) {
                /* TODO: FREE ALLOCATED RESOURCES */
                return -1;
            }
            /* --END ERROR HANDLING-- */
        } else {
            *depth += 1;
            MPIR_Dame_create_contiguous(blklens[i],
                                        oldtypes[i],
                                        &(*dl)[d].s.s_t.dls[loop_idx], &inner_dlsize, depth);
        }
        *depth = d;
        (*dl)[*depth].s.s_t.dls[loop_idx][*depth].kind = DL_RETURNTO;
        (*dl)[*depth].s.s_t.dls[loop_idx][*depth].returnto = *depth;

        tmpoldsizes[loop_idx] = oldsize;
        tmpoffsets[loop_idx] = disps[i];
        tmpblklens[loop_idx] = blklens[i];
        loop_idx++;
    }
    MPIR_Assert(loop_idx == (nr_basics + nr_derived));
    (*dl)[*depth].count = nr_basics + nr_derived;
    (*dl)[*depth].size = size;
    (*dl)[*depth].extent = extent;

    (*dl)[*depth].s.s_t.oldsizes = tmpoldsizes;
    (*dl)[*depth].s.s_t.offsets = tmpoffsets;
    (*dl)[*depth].s.s_t.blklens = tmpblklens;

    *depth += 1;

    return 0;
}


/* --BEGIN ERROR HANDLING-- */
static int Dame_create_struct_memory_error(void)
{
    return -1;
}

/* --END ERROR HANDLING-- */

static int Dame_create_unique_type_struct(DAME_Count count,
                                          const int *blklens,
                                          const MPI_Aint * disps,
                                          MPI_Datatype type,
                                          const DAME_Type * oldtypes,
                                          int type_pos,
                                          DAME_Dame ** dl, MPI_Aint * size, int *depth)
{
    /* the same type used more than once in the array; type_pos
     * indexes to the first of these.
     */
    int i, err, cur_pos = 0;
    DAME_Size *tmp_blklens;
    DAME_Offset *tmp_disps;

    /* count is an upper bound on number of type instances */
    tmp_blklens = (DAME_Size *) DAME_Malloc(count * sizeof(DAME_Size));
    /* --BEGIN ERROR HANDLING-- */
    if (!tmp_blklens) {
        /* TODO: ??? */
        return Dame_create_struct_memory_error();
    }
    /* --END ERROR HANDLING-- */

    tmp_disps = (DAME_Offset *)
        DAME_Malloc(count * sizeof(DAME_Offset));
    /* --BEGIN ERROR HANDLING-- */
    if (!tmp_disps) {
        DAME_Free(tmp_blklens);
        /* TODO: ??? */
        return Dame_create_struct_memory_error();
    }
    /* --END ERROR HANDLING-- */

    for (i = type_pos; i < count; i++) {
        if (oldtypes[i] == oldtypes[type_pos] && blklens != 0) {
            tmp_blklens[cur_pos] = blklens[i];
            tmp_disps[cur_pos] = disps[i];
            cur_pos++;
        }
    }

    err = MPIR_Dame_create_indexed(cur_pos, tmp_blklens, tmp_disps, 1,  /* disp in bytes */
                                   type, oldtypes[type_pos], dl, size, depth);

    DAME_Free(tmp_blklens);
    DAME_Free(tmp_disps);

    return err;
}

static int Dame_create_contig_all_bytes_struct(DAME_Count count,
                                               const int *blklens,
                                               const MPI_Aint * disps,
                                               MPI_Datatype type,
                                               const DAME_Type * oldtypes,
                                               DAME_Dame ** dl, MPI_Aint * size, int *depth)
{
    int i, err, cur_pos = 0;
    DAME_Size *tmp_blklens;
    MPI_Aint *tmp_disps;

    /* count is an upper bound on number of type instances */
    tmp_blklens = (DAME_Size *) DAME_Malloc(count * sizeof(DAME_Size));

    /* --BEGIN ERROR HANDLING-- */
    if (!tmp_blklens) {
        return Dame_create_struct_memory_error();
    }
    /* --END ERROR HANDLING-- */

    tmp_disps = (MPI_Aint *) DAME_Malloc(count * sizeof(MPI_Aint));

    /* --BEGIN ERROR HANDLING-- */
    if (!tmp_disps) {
        DAME_Free(tmp_blklens);
        return Dame_create_struct_memory_error();
    }
    /* --END ERROR HANDLING-- */

    for (i = 0; i < count; i++) {
        if (oldtypes[i] != MPI_LB && oldtypes[i] != MPI_UB && blklens[i] != 0) {
            DAME_Offset sz, extent;
            MPIR_Datatype *old_ptr;
            MPIR_Datatype_get_ptr(oldtypes[i], old_ptr);

            DAME_Handle_get_size_macro(oldtypes[i], sz);
            DAME_Handle_get_extent_macro(oldtypes[i], extent);
            tmp_blklens[cur_pos] = sz * blklens[i];
            tmp_disps[cur_pos] = disps[i] + old_ptr->lb;
            cur_pos++;
        }
    }

    err = MPIR_Dame_create_indexed(cur_pos, tmp_blklens, tmp_disps, 1,  /* disp in bytes */
                                   type, MPI_BYTE, dl, size, depth);

    DAME_Free(tmp_blklens);
    DAME_Free(tmp_disps);

    return err;
}
