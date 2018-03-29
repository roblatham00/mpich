/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpiimpl.h"

#undef DAME_DEBUG_MEMORY

typedef unsigned char byte;

#define GET_COMPACTED_BUF(s) ((void*)((byte*)(&(s)) + (long)(s)))
#define SET_COMPACTED_BUF(s, b) s = (void*)((byte*)(b) - (byte*)(&(s)))

/* Dataloops
 *
 * The functions here are used for the creation, copying, update, and display
 * of Dame structures and trees of these structures.
 *
 */


/*@
  Dame_free - deallocate the resources used to store a dataloop

  Input/output Parameters:
  . dataloop - pointer to dataloop structure
  @*/
void MPIR_Dame_free(DAME_Dame ** dl, int is_compact)
{
    if (*dl == NULL)
        return;

#ifdef DAME_DEBUG_MEMORY
    DAME_dbg_printf("Dame_free: freeing loop @ %x.\n", (int) *dl);
#endif
    unsigned i = 0, j = 0;
    /* why the memset? */
    for (i = 0; i < DAME_MAX_DEPTH; i++) {
        switch ((*dl)[i].kind) {
            case DL_BLOCKINDEX:
            case DL_BLOCKINDEX1:
            case DL_BLOCKINDEXFINAL:
                is_compact ?  DAME_Free(
				GET_COMPACTED_BUF((*dl)[i].s.bi_t.offsets)) :
			DAME_Free((void *) ((*dl)[i].s.bi_t.offsets)) ;
                (*dl)[i].s.bi_t.offsets = NULL;
                break;
            case DL_INDEX:
            case DL_INDEXFINAL:
                is_compact ?  DAME_Free(
				GET_COMPACTED_BUF((*dl)[i].s.i_t.offsets) ) :
			DAME_Free((void *) ((*dl)[i].s.i_t.offsets)) ;
                (*dl)[i].s.i_t.offsets = NULL;
                is_compact ? DAME_Free(
				GET_COMPACTED_BUF((*dl)[i].s.i_t.blklens) ) :
			DAME_Free((void *) ((*dl)[i].s.i_t.blklens));
                (*dl)[i].s.i_t.blklens = NULL;
                break;
            case DL_STRUCT:
                is_compact ? DAME_Free(
				GET_COMPACTED_BUF((*dl)[i].s.s_t.oldsizes)) :
			DAME_Free((void *) ((*dl)[i].s.s_t.oldsizes));
                is_compact ? DAME_Free(
				GET_COMPACTED_BUF((*dl)[i].s.s_t.blklens)) :
			DAME_Free((void *) ((*dl)[i].s.s_t.blklens));
                is_compact ?  DAME_Free(
				GET_COMPACTED_BUF((*dl)[i].s.s_t.offsets)) :
			DAME_Free((void *) ((*dl)[i].s.s_t.offsets));
                (*dl)[i].s.s_t.offsets = NULL;
                for (j = 0; j < (*dl)[i].count; j++) {
                    MPIR_Dame_free(&((*dl)[i].s.s_t.dls[j]), is_compact);
                    (*dl)[i].s.s_t.dls[j] = NULL;
                }
                is_compact ? DAME_Free(
				GET_COMPACTED_BUF((*dl)[i].s.s_t.dls)) :
			DAME_Free((void *) ((*dl)[i].s.s_t.dls));
            default:
                break;
        }
    }

    memset(*dl, 0, DAME_MAX_DEPTH * sizeof(DAME_Dame));
    DAME_Free(*dl);
    *dl = NULL;
    return;
}

/* Allocate space for a complete Dame program.
 * Returns 0 on success, -1 on error
 */
int MPIR_Dame_alloc(DAME_Dame ** new_dl)
{
    *new_dl = (DAME_Dame *) DAME_Malloc(DAME_MAX_DEPTH * sizeof(DAME_Dame));
    if (!*new_dl)
        return -1;
    unsigned i = 0;
    (*new_dl)[0].kind = DL_EXIT;
    for (i = 1; i < DAME_MAX_DEPTH; i++) {
        (*new_dl)[i].kind = DL_BOTTOM;
        (*new_dl)[i].count = 0;
        (*new_dl)[i].flags = 0x0;
    }
    return 0;
}

/* Allocate space for a complete Dame program.
 * Returns 0 on success, -1 on error
 */
int MPIR_Dame_alloc_compact(MPI_Aint dataloop_size, DAME_Dame ** new_dl)
{
    *new_dl = (DAME_Dame *) DAME_Malloc(dataloop_size);
    if (!*new_dl)
        return -1;
    return 0;
}

/* Allocate space for a struct dataloop.
 * Returns 0 on success, -1 on error
 */
int MPIR_Dame_struct_alloc(DAME_Count count, DAME_Dame * dl)
{
    dl->kind = DL_STRUCT;
    dl->s.s_t.dls = (DAME_Dame **) MPL_malloc(count * sizeof(DAME_Dame *), MPL_MEM_DATATYPE);
    if (!dl->s.s_t.dls)
        return -1;

    unsigned i;
    for (i = 0; i < count; i++)
        MPIR_Dame_alloc(&dl->s.s_t.dls[i]);

    return 0;
}

/*@
  Dame_dup - make a copy of a dataloop

  Returns 0 on success, -1 on failure.
  @*/
int MPIR_Dame_dup(DAME_Dame * old_dl, DAME_Size unused, DAME_Dame ** new_dl)
{
    if (*new_dl == 0) {
        if (MPIR_Dame_alloc(new_dl) == -1)
            return -1;
    }
    int i, j;
    for (i = 0; i < DAME_MAX_DEPTH && old_dl[i].kind != DL_BOTTOM; i++) {
        (*new_dl)[i].kind = old_dl[i].kind;
        (*new_dl)[i].count = old_dl[i].count;
        (*new_dl)[i].size = old_dl[i].size;
        (*new_dl)[i].extent = old_dl[i].extent;
        (*new_dl)[i].returnto = old_dl[i].returnto;
        switch (old_dl[i].kind) {
            case DL_CONTIG:
            case DL_CONTIGFINAL:
                (*new_dl)[i].s.c_t.basesize = old_dl[i].s.c_t.basesize;
                (*new_dl)[i].s.c_t.baseextent = old_dl[i].s.c_t.baseextent;
                break;
            case DL_VECTOR:
            case DL_VECTOR1:
            case DL_VECTORFINAL:
                (*new_dl)[i].s.v_t.oldsize = old_dl[i].s.v_t.oldsize;
                (*new_dl)[i].s.v_t.blklen = old_dl[i].s.v_t.blklen;
                (*new_dl)[i].s.v_t.stride = old_dl[i].s.v_t.stride;
                break;
            case DL_BLOCKINDEX:
            case DL_BLOCKINDEX1:
            case DL_BLOCKINDEXFINAL:
                (*new_dl)[i].s.bi_t.oldsize = old_dl[i].s.bi_t.oldsize;
                (*new_dl)[i].s.bi_t.blklen = old_dl[i].s.bi_t.blklen;
                (*new_dl)[i].s.bi_t.offsets =
                    (MPI_Aint *) MPL_malloc(old_dl[i].count * sizeof(MPI_Aint),
				    MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.bi_t.offsets)
                    goto error_exit;
                memcpy((MPI_Aint *) (*new_dl)[i].s.bi_t.offsets,
                       old_dl[i].s.bi_t.offsets, old_dl[i].count * sizeof(MPI_Aint));
                break;
            case DL_INDEX:
            case DL_INDEXFINAL:
                (*new_dl)[i].s.i_t.oldsize = old_dl[i].s.i_t.oldsize;
                (*new_dl)[i].s.i_t.blklens =
                    (MPI_Aint *) MPL_malloc(old_dl[i].count * sizeof(MPI_Aint),
				    MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.i_t.blklens)
                    goto error_exit;
                memcpy((MPI_Aint *) (*new_dl)[i].s.i_t.blklens,
                       old_dl[i].s.i_t.blklens, old_dl[i].count * sizeof(MPI_Aint));
                (*new_dl)[i].s.i_t.offsets =
                    (MPI_Aint *) MPL_malloc(old_dl[i].count * sizeof(MPI_Aint),
				    MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.i_t.offsets)
                    goto error_exit;
                memcpy((MPI_Aint *) (*new_dl)[i].s.i_t.offsets,
                       old_dl[i].s.i_t.offsets, old_dl[i].count * sizeof(MPI_Aint));
                break;
            case DL_STRUCT:
                (*new_dl)[i].s.s_t.oldsizes =
                    (MPI_Aint *) MPL_malloc(old_dl[i].count * sizeof(MPI_Aint),
				    MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.s_t.oldsizes)
                    goto error_exit;
                memcpy((MPI_Aint *) (*new_dl)[i].s.s_t.oldsizes,
                       old_dl[i].s.s_t.oldsizes, old_dl[i].count * sizeof(MPI_Aint));
                (*new_dl)[i].s.s_t.blklens =
                    (MPI_Aint *) MPL_malloc(old_dl[i].count * sizeof(MPI_Aint),
				    MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.s_t.blklens)
                    goto error_exit;
                memcpy((MPI_Aint *) (*new_dl)[i].s.s_t.blklens,
                       old_dl[i].s.s_t.blklens, old_dl[i].count * sizeof(MPI_Aint));
                (*new_dl)[i].s.s_t.offsets =
                    (MPI_Aint *) MPL_malloc(old_dl[i].count * sizeof(MPI_Aint), MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.s_t.offsets)
                    goto error_exit;
                memcpy((MPI_Aint *) (*new_dl)[i].s.s_t.offsets,
                       old_dl[i].s.s_t.offsets, old_dl[i].count * sizeof(MPI_Aint));

                MPIR_Dame_struct_alloc(old_dl[i].count, &(*new_dl)[i]);
                for (j = 0; j < old_dl[i].count; j++)
                    MPIR_Dame_dup(old_dl[i].s.s_t.dls[j], unused, &(*new_dl)[i].s.s_t.dls[j]);
                break;
            case DL_CONTIGCHILD:
                break;
            case DL_EXIT:
                break;
            case DL_BOTTOM:
                goto clean_exit;
                break;
            default:
                break;
        }
    }
  clean_exit:
    return 0;
  error_exit:
    return -1;
}

/*@
  Dame_dup_compact - make a copy of the compact dataloop

  Returns 0 on success, -1 on failure.
  @*/
int MPIR_Dame_dup_compact(DAME_Dame * old_compact,
                          MPI_Aint old_dataloop_size, DAME_Dame ** new_compact)
{
    if (*new_compact == 0) {
        if (MPIR_Dame_alloc_compact(old_dataloop_size, new_compact) != 0)
            return -1;
    }
    DAME_Memcpy(*new_compact, old_compact, old_dataloop_size);
    return 0;
}

/*@
  Dataloop_stream_size - return the size of the data described by the dataloop

Input Parameters:
+ dl_p   - pointer to dataloop for which we will return the size
- sizefn - function for determining size of types in the corresponding stream
           (passing NULL will instead result in el_size values being used)

@*/
DAME_Offset MPIR_Dame_stream_size(DAME_Dame * dl_p)
{
    fprintf(stderr, "Dame_stream_size unimplemented\n");
    DAME_Assert(0);
    return 0;
}

/* --BEGIN ERROR HANDLING-- */
/*@
  Dame_print - dump a dataloop tree to stdout for debugging
  purposes

  Input Parameters:
  + dataloop - root of tree to dump
  @*/
void MPIR_Dame_print(DAME_Dame * dl)
{
    int i, j;
    int sp = 0;
    do {
        switch (dl[sp].kind) {
            case DL_CONTIG:
                fprintf(stdout, "CONTIG(%d): %ld, %ld, %ld\n\t%ld, %ld\n",
                        sp,
                        (long) dl[sp].count, (long) dl[sp].size,
                        (long) dl[sp].extent,
                        (long) dl[sp].s.c_t.basesize, (long) dl[sp].s.c_t.baseextent);
                break;
            case DL_VECTOR:
                fprintf(stdout, "VEC(%d): %ld, %ld, %ld\n\t%ld, %ld\n\t%ld\n",
                        sp,
                        (long) dl[sp].count, (long) dl[sp].size,
                        (long) dl[sp].extent,
                        (long) dl[sp].s.v_t.blklen, (long) dl[sp].s.v_t.oldsize,
                        (long) dl[sp].s.v_t.stride);
                break;
            case DL_VECTOR1:
                fprintf(stdout, "VEC1(%d): %ld, %ld, %ld\n\t%ld\n\t%ld\n",
                        sp,
                        (long) dl[sp].count, (long) dl[sp].size,
                        (long) dl[sp].extent,
                        (long) dl[sp].s.v_t.oldsize, (long) dl[sp].s.v_t.stride);
                break;
            case DL_BLOCKINDEX:
                fprintf(stdout, "BLKINDEX(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld, %ld\n",
                        (long) dl[sp].s.bi_t.blklen, (long) dl[sp].s.bi_t.oldsize);
                fprintf(stdout, "\t%ld ", dl[sp].s.bi_t.offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", dl[sp].s.bi_t.offsets[i]);
                fprintf(stdout, "\n");
                break;
            case DL_BLOCKINDEX1:
                fprintf(stdout, "BLKINDEX1(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld\n", (long) dl[sp].s.bi_t.oldsize);
                fprintf(stdout, "\t%ld ", dl[sp].s.bi_t.offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", dl[sp].s.bi_t.offsets[i]);
                fprintf(stdout, "\n");
                break;
            case DL_INDEX:
                fprintf(stdout, "INDEX(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld\n", dl[sp].s.i_t.oldsize);
                fprintf(stdout, "\t%ld ", dl[sp].s.i_t.blklens[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", dl[sp].s.i_t.blklens[i]);
                fprintf(stdout, "\n");
                fprintf(stdout, "\t%ld ", dl[sp].s.i_t.offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", dl[sp].s.i_t.offsets[i]);
                fprintf(stdout, "\n");
                break;
            case DL_VECTORFINAL:
                fprintf(stdout, "VECFINAL(%d): %ld, %ld, %ld\n\t%ld, %ld\n\t%ld\n",
                        sp,
                        (long) dl[sp].count, (long) dl[sp].size,
                        (long) dl[sp].extent,
                        (long) dl[sp].s.v_t.blklen, (long) dl[sp].s.v_t.oldsize,
                        (long) dl[sp].s.v_t.stride);
                break;
            case DL_BLOCKINDEXFINAL:
                fprintf(stdout, "BLKINDEXFINAL(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld, %ld\n",
                        (long) dl[sp].s.bi_t.blklen, (long) dl[sp].s.bi_t.oldsize);
                fprintf(stdout, "\t%ld ", dl[sp].s.bi_t.offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", dl[sp].s.bi_t.offsets[i]);
                fprintf(stdout, "\n");
                break;
            case DL_CONTIGFINAL:
                fprintf(stdout, "CONTIGFINAL(%d): %ld, %ld, %ld\n\t%ld, %ld\n", sp,
                        (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent,
                        (long) dl[sp].s.c_t.basesize, (long) dl[sp].s.c_t.baseextent);
                break;
            case DL_INDEXFINAL:
                fprintf(stdout, "INDEXFINAL(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld", dl[sp].s.i_t.oldsize);
                fprintf(stdout, "\t%ld ", dl[sp].s.i_t.blklens[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", dl[sp].s.i_t.blklens[i]);
                fprintf(stdout, "\n");
                fprintf(stdout, "\t%ld ", dl[sp].s.i_t.offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", dl[sp].s.i_t.offsets[i]);
                fprintf(stdout, "\n");
                break;
            case DL_CONTIGCHILD:
                fprintf(stdout, "CONTIGCHILD(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                break;
            case DL_STRUCT:
                fprintf(stdout, "STRUCT(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld ", dl[sp].s.s_t.blklens[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", dl[sp].s.s_t.blklens[i]);
                fprintf(stdout, "\n");
                fprintf(stdout, "\t%ld ", dl[sp].s.s_t.offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", dl[sp].s.s_t.offsets[i]);
                fprintf(stdout, "\n");
                for (j = 0; j < dl[sp].count; j++) {
                    fprintf(stdout, " --- %d --- \n", j);
                    MPIR_Dame_print(dl[sp].s.s_t.dls[j]);
                    fprintf(stdout, "\n");
                }
                break;
            case DL_RETURNTO:
                fprintf(stdout, "RETURNTO(%d): %ld\n", sp, dl[sp].returnto);
                break;
            case DL_EXIT:
                fprintf(stdout, "EXIT(%d)\n", sp);
                break;
            case DL_BOTTOM:
                fprintf(stdout, "BOTTOM(%d)\n", sp);
                break;
            default:
                fprintf(stdout, "UNKNOWN(%d)\n", sp);
                break;
        }
        sp++;
    } while (dl[sp].kind != DL_BOTTOM && sp < DAME_MAX_DEPTH);
    fprintf(stdout, "BOTTOM(%d)\n", sp);

}

/* --END ERROR HANDLING-- */


/* --BEGIN ERROR HANDLING-- */
/*@
  Dame_print_compact - dump a serialized DAME program to stdout for
  debugging purposes

  Input Parameters:
  + dataloop - root of tree to dump
  @*/
void MPIR_Dame_print_compact(DAME_Dame * dl)
{
    int i, j;
    int sp = 0;
    MPI_Aint *offsets = NULL;
    MPI_Aint *blklens = NULL;
    DAME_Dame **dls = NULL;
    do {
        switch (dl[sp].kind) {
            case DL_CONTIG:
                fprintf(stdout, "CONTIG(%d): %ld, %ld, %ld\n\t%ld, %ld\n",
                        sp,
                        (long) dl[sp].count, (long) dl[sp].size,
                        (long) dl[sp].extent,
                        (long) dl[sp].s.c_t.basesize, (long) dl[sp].s.c_t.baseextent);
                break;
            case DL_VECTOR:
                fprintf(stdout, "VEC(%d): %ld, %ld, %ld\n\t%ld, %ld\n\t%ld\n",
                        sp,
                        (long) dl[sp].count, (long) dl[sp].size,
                        (long) dl[sp].extent,
                        (long) dl[sp].s.v_t.blklen, (long) dl[sp].s.v_t.oldsize,
                        (long) dl[sp].s.v_t.stride);
                break;
            case DL_VECTOR1:
                fprintf(stdout, "VEC1(%d): %ld, %ld, %ld\n\t%ld\n\t%ld\n",
                        sp,
                        (long) dl[sp].count, (long) dl[sp].size,
                        (long) dl[sp].extent,
                        (long) dl[sp].s.v_t.oldsize, (long) dl[sp].s.v_t.stride);
                break;
            case DL_BLOCKINDEX:
                fprintf(stdout, "BLKINDEX(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld, %ld\n",
                        (long) dl[sp].s.bi_t.blklen, (long) dl[sp].s.bi_t.oldsize);
                offsets = GET_COMPACTED_BUF(dl[sp].s.bi_t.offsets);
                fprintf(stdout, "\t%ld ", offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", offsets[i]);
                fprintf(stdout, "\n");
                break;
            case DL_BLOCKINDEX1:
                fprintf(stdout, "BLKINDEX1(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld, %ld\n",
                        (long) dl[sp].s.bi_t.blklen, (long) dl[sp].s.bi_t.oldsize);
                offsets = GET_COMPACTED_BUF(dl[sp].s.bi_t.offsets);
                fprintf(stdout, "\t%ld ", offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", offsets[i]);
                fprintf(stdout, "\n");
                break;
            case DL_INDEX:
                fprintf(stdout, "INDEX(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld", dl[sp].s.i_t.oldsize);
                blklens = GET_COMPACTED_BUF(dl[sp].s.i_t.blklens);
                fprintf(stdout, "\t%ld ", blklens[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", blklens[i]);
                fprintf(stdout, "\n");
                offsets = GET_COMPACTED_BUF(dl[sp].s.i_t.offsets);
                fprintf(stdout, "\t%ld ", offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", offsets[i]);
                fprintf(stdout, "\n");
                break;
            case DL_VECTORFINAL:
                fprintf(stdout, "VECFINAL(%d): %ld, %ld, %ld\n\t%ld, %ld\n\t%ld\n",
                        sp,
                        (long) dl[sp].count, (long) dl[sp].size,
                        (long) dl[sp].extent,
                        (long) dl[sp].s.v_t.blklen, (long) dl[sp].s.v_t.oldsize,
                        (long) dl[sp].s.v_t.stride);
                break;
            case DL_BLOCKINDEXFINAL:
                fprintf(stdout, "BLKINDEXFINAL(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld, %ld\n",
                        (long) dl[sp].s.bi_t.blklen, (long) dl[sp].s.bi_t.oldsize);
                offsets = GET_COMPACTED_BUF(dl[sp].s.bi_t.offsets);
                fprintf(stdout, "\t%ld ", offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", offsets[i]);
                fprintf(stdout, "\n");
                break;
            case DL_CONTIGFINAL:
                fprintf(stdout, "CONTIGFINAL(%d): %ld, %ld, %ld\n\t%ld, %ld\n", sp,
                        (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent,
                        (long) dl[sp].s.c_t.basesize, (long) dl[sp].s.c_t.baseextent);
                break;
            case DL_INDEXFINAL:
                fprintf(stdout, "INDEXFINAL(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                fprintf(stdout, "\t%ld", dl[sp].s.i_t.oldsize);
                blklens = GET_COMPACTED_BUF(dl[sp].s.i_t.blklens);
                fprintf(stdout, "\t%ld ", blklens[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", blklens[i]);
                fprintf(stdout, "\n");
                offsets = GET_COMPACTED_BUF(dl[sp].s.i_t.offsets);
                fprintf(stdout, "\t%ld ", offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", offsets[i]);
                fprintf(stdout, "\n");
                break;
            case DL_CONTIGCHILD:
                fprintf(stdout, "CONTIGCHILD(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                break;
            case DL_STRUCT:
                fprintf(stdout, "STRUCT(%d): %ld, %ld, %ld\n",
                        sp, (long) dl[sp].count, (long) dl[sp].size, (long) dl[sp].extent);
                blklens = GET_COMPACTED_BUF(dl[sp].s.s_t.blklens);
                fprintf(stdout, "\t%ld ", blklens[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", blklens[i]);
                fprintf(stdout, "\n");
                offsets = GET_COMPACTED_BUF(dl[sp].s.s_t.offsets);
                fprintf(stdout, "\t%ld ", offsets[0]);
                for (i = 1; i < dl[sp].count; i++)
                    fprintf(stdout, "%ld ", offsets[i]);
                fprintf(stdout, "\n");
                dls = GET_COMPACTED_BUF(dl[sp].s.s_t.dls);
                for (j = 0; j < dl[sp].count; j++) {
                    fprintf(stdout, " --- %d --- \n", j);
                    MPIR_Dame_print_compact(GET_COMPACTED_BUF(dls[j]));
                    fprintf(stdout, "\n");
                }
                break;
            case DL_RETURNTO:
                fprintf(stdout, "RETURNTO(%d): %ld\n", sp, dl[sp].returnto);
                break;
            case DL_EXIT:
                fprintf(stdout, "EXIT(%d)\n", sp);
                break;
            case DL_BOTTOM:
                fprintf(stdout, "BOTTOM(%d)\n", sp);
                break;
            default:
                fprintf(stdout, "UNKNOWN(%d)\n", sp);
                break;
        }
        sp++;
    } while (dl[sp].kind != DL_BOTTOM && sp < DAME_MAX_DEPTH);
    fprintf(stdout, "BOTTOM(%d)\n", sp);

}

/* --END ERROR HANDLING-- */


/* Moves a dataloop specification from old to new
 * Returns 0 on success, -1 on failure
 */
int MPIR_Dame_move(DAME_Dame * old, DAME_Dame * new)
{
    int j;
    for (j = 0; j < DAME_MAX_DEPTH && old[j].kind != DL_BOTTOM;)
        j++;
    j++;
    DAME_Dame *tmp;
    if (MPIR_Dame_alloc(&tmp) == -1)
        return -1;
    memcpy(tmp, old, j * sizeof(DAME_Dame));
    memcpy(new, tmp, j * sizeof(DAME_Dame));
    MPL_free(tmp);
    return 0;
}

/* A DAME program with STRUCT's will look like this:

   0: EXIT
   1: CONTIG
   2: STRUCT {
   --- 0 ---
   0: EXIT
   1: RETURN(1)
   2: VECFINAL
   3: CONTIGFINAL
   4: BOTTOM
   --- 1 ---
   0: EXIT
   1: RETURN(1)
   2: CONTIGFINAL
   3: BOTTOM
   }
   3: BOTTOM

   where the n: indicates the n'th slot in the stack which is used in
   DL_Pack and DL_Unpack to process this type. --- m --- indicates the m'th
   element of the struct.

   This can make things difficult when resuming from deeply nested struct
   types because you have to keep track of how deep in the stack you are and then
   have to save each return point.

   To avoid this, we adjust the dataloop so that the pointer to the type being
   processed increases or decreases monotonically. The resulting program is

   0: EXIT
   1: CONTIG
   2: STRUCT {
   --- 0 ---
   0: EXIT
   1: EXIT
   2: RETURNTO(2)
   3: VECFINAL
   4: CONTIGFINAL
   5: BOTTOM
   --- 1 ---
   0: EXIT
   1: EXIT
   2: RETURNTO(2)
   3: CONTIGFINAL
   4: BOTTOM
   }
   3: BOTTOM

   Return 0 on success, -1 on failure
*/
int MPIR_Dame_update(DAME_Dame * dl, MPI_Aint depth)
{
    int i, j;
    for (i = 0; i < DAME_MAX_DEPTH && dl[i].kind != DL_BOTTOM; i++) {
        switch (dl[i].kind) {
            case DL_RETURNTO:
                if (dl[i].returnto != depth) {
                    dl[i].returnto = depth;
                    if (MPIR_Dame_move(&dl[i], &dl[depth]) == -1)
                        return -1;
                    for (j = i; j < depth; j++)
                        dl[j].kind = DL_EXIT;
                }
                break;

            case DL_STRUCT:
                for (j = 0; j < dl[i].count; j++) {
                    if (MPIR_Dame_update(dl[i].s.s_t.dls[j], i) == -1)
                        return -1;
                }
                break;
            default:
                break;
        }
    }
    return 0;
}

void MPIR_Dame_calculate_size(DAME_Dame * dl, MPI_Aint * dataloop_size, int only_primitives)
{
    int i;
    int sp = 0;
    do {
        switch (dl[sp].kind) {
            case DL_CONTIG:
            case DL_CONTIGFINAL:
            case DL_VECTOR:
            case DL_VECTOR1:
            case DL_VECTORFINAL:
                break;
            case DL_BLOCKINDEX:
            case DL_BLOCKINDEX1:
            case DL_BLOCKINDEXFINAL:
                if (!only_primitives)
                    *dataloop_size += (dl[sp].count * sizeof(DAME_Offset));
                break;
            case DL_INDEX:
            case DL_INDEXFINAL:
                if (!only_primitives) {
                    *dataloop_size += (dl[sp].count * sizeof(DAME_Offset));
                    *dataloop_size += (dl[sp].count * sizeof(DAME_Size));
                }
                break;
            case DL_STRUCT:
                for (i = 0; i < dl[sp].count; i++)
                    MPIR_Dame_calculate_size(dl[sp].s.s_t.dls[i], dataloop_size, only_primitives);
                if (!only_primitives) {
                    *dataloop_size += (dl[sp].count * sizeof(MPI_Aint));
                    *dataloop_size += (dl[sp].count * sizeof(DAME_Offset));
                    *dataloop_size += (dl[sp].count * sizeof(DAME_Size));
                    *dataloop_size += (dl[sp].count * sizeof(DAME_Dame *));
                }
                break;
            case DL_CONTIGCHILD:
            case DL_RETURNTO:
            case DL_EXIT:
                break;
            default:
                DAME_Assert(0);
                break;
        }
        *dataloop_size += sizeof(DAME_Dame);
        sp += 1;
    } while (dl[sp].kind != DL_BOTTOM && sp < DAME_MAX_DEPTH);
    /* One last one for BOTTOM */
    *dataloop_size += sizeof(DAME_Dame);
}

/* Actual implementation of the serialization
 *
 * The format of the serialized DAME program is as follows:
 *
 * The serialized size of the program is dataloop_size bytes
 * The first m bytes are a code section which contains
 * all the DAME primitives that constitute the program (including all the
 * primitives within struct types.
 * The remaining (dataloop_size-m) bytes is a data section which contains
 * offsets, blocklengths and anything else that indexed and struct types need
 *
 *
 * 0                 m            dataloop_size
 * + ------------------------------------+
 * | DAME primitives |       data        |
 * + ------------------------------------+
 *
 *
 * The pointers in the indexed and struct primitives contain the offset from
 * that field to the start of the data in the data section. The offset is
 * guaranteed to be positive.
 *
 * Returns 0 on success, -1 on failure.
 */
static int MPIR_Dame_serialize_impl(DAME_Dame * dl,
                                    MPI_Aint dataloop_size, void **curr_buf, DAME_Dame ** curr_dl)
{
    int i;
    int sp = 0;
    DAME_Dame *saved_dl = NULL;
    DAME_Dame *saved_buf = NULL;
    DAME_Dame *allocated = NULL;

    /* If this is a recursive call from a struct that we are already processing,
     * all the pointers will have been set correctly and we can just use it.
     * If not, set everything up */
    if (!*curr_dl) {
        MPI_Aint buff_offset = 0;
        MPIR_Dame_calculate_size(dl, &buff_offset, 1);
        if (MPIR_Dame_alloc_compact(dataloop_size, &allocated) != 0)
            return -1;
        *curr_dl = allocated;
        *curr_buf = ((byte *) (*curr_dl) + buff_offset);
    }

    do {
        (*curr_dl)->kind = dl[sp].kind;
        (*curr_dl)->count = dl[sp].count;
        (*curr_dl)->size = dl[sp].size;
        (*curr_dl)->extent = dl[sp].extent;
        (*curr_dl)->returnto = dl[sp].returnto;
        switch (dl[sp].kind) {
            case DL_CONTIG:
            case DL_CONTIGFINAL:
                (*curr_dl)->s.c_t.basesize = dl[sp].s.c_t.basesize;
                (*curr_dl)->s.c_t.baseextent = dl[sp].s.c_t.baseextent;
                break;
            case DL_VECTOR:
            case DL_VECTOR1:
            case DL_VECTORFINAL:
                (*curr_dl)->s.v_t.oldsize = dl[sp].s.v_t.oldsize;
                (*curr_dl)->s.v_t.blklen = dl[sp].s.v_t.blklen;
                (*curr_dl)->s.v_t.stride = dl[sp].s.v_t.stride;
                break;
            case DL_BLOCKINDEX:
            case DL_BLOCKINDEX1:
            case DL_BLOCKINDEXFINAL:
                (*curr_dl)->s.bi_t.oldsize = dl[sp].s.bi_t.oldsize;
                (*curr_dl)->s.bi_t.blklen = dl[sp].s.bi_t.blklen;
                memcpy((MPI_Aint *) (*curr_buf),
                       dl[sp].s.bi_t.offsets, dl[sp].count * sizeof(MPI_Aint));
                SET_COMPACTED_BUF((*curr_dl)->s.bi_t.offsets, *curr_buf);
                (*curr_buf) = ((byte *) * curr_buf) + dl[sp].count * sizeof(MPI_Aint);
                break;
            case DL_INDEX:
            case DL_INDEXFINAL:
                (*curr_dl)->s.i_t.oldsize = dl[sp].s.i_t.oldsize;
                memcpy((MPI_Aint *) (*curr_buf),
                       dl[sp].s.i_t.blklens, dl[sp].count * sizeof(MPI_Aint));
                SET_COMPACTED_BUF((*curr_dl)->s.i_t.blklens, *curr_buf);
                (*curr_buf) = ((byte *) * curr_buf) + dl[sp].count * sizeof(MPI_Aint);
                memcpy((MPI_Aint *) (*curr_buf),
                       dl[sp].s.i_t.offsets, dl[sp].count * sizeof(MPI_Aint));
                SET_COMPACTED_BUF((*curr_dl)->s.i_t.offsets, *curr_buf);
                (*curr_buf) = ((byte *) * curr_buf) + dl[sp].count * sizeof(MPI_Aint);
                break;
            case DL_STRUCT:
                memcpy((MPI_Aint *) (*curr_buf),
                       dl[sp].s.s_t.oldsizes, dl[sp].count * sizeof(MPI_Aint));
                SET_COMPACTED_BUF((*curr_dl)->s.s_t.oldsizes, *curr_buf);
                (*curr_buf) = ((byte *) * curr_buf) + dl[sp].count * sizeof(MPI_Aint);
                memcpy((MPI_Aint *) (*curr_buf),
                       dl[sp].s.s_t.blklens, dl[sp].count * sizeof(MPI_Aint));
                SET_COMPACTED_BUF((*curr_dl)->s.s_t.blklens, *curr_buf);
                (*curr_buf) = ((byte *) * curr_buf) + dl[sp].count * sizeof(MPI_Aint);
                memcpy((MPI_Aint *) (*curr_buf),
                       dl[sp].s.s_t.offsets, dl[sp].count * sizeof(MPI_Aint));
                SET_COMPACTED_BUF((*curr_dl)->s.s_t.offsets, *curr_buf);
                (*curr_buf) = ((byte *) * curr_buf) + dl[sp].count * sizeof(MPI_Aint);
                saved_buf = *curr_buf;
                (*curr_dl)->s.s_t.dls = *curr_buf;
                (*curr_buf) = ((byte *) * curr_buf) + dl[sp].count * sizeof(DAME_Dame *);
                /* A struct is guaranteed to be the last non-bottom element in a
                 * DAME program by construction. The next available slot will
                 * have to be the one after the next BOTTOM primitive.
                 * It doesn't really matter in what order the subprograms appear
                 * as long as each subprogram is contiguous. So I'll just write out
                 * the BOTTOM here and exit */
                saved_dl = *curr_dl;
                (*curr_dl) += 1;
                (*curr_dl)->kind = DL_BOTTOM;
                (*curr_dl) += 1;
                for (i = 0; i < dl[sp].count; i++) {
                    SET_COMPACTED_BUF(saved_dl->s.s_t.dls[i], *curr_dl);
                    if (MPIR_Dame_serialize_impl(dl[sp].s.s_t.dls[i],
                                                 dataloop_size, curr_buf, curr_dl) == -1)
                        return -1;
                }
                SET_COMPACTED_BUF(saved_dl->s.s_t.dls, saved_buf);
                goto clean_exit;
                break;
            case DL_CONTIGCHILD:
            case DL_EXIT:
            case DL_RETURNTO:
                break;
            default:
                DAME_Assert(0);
        }
        (*curr_dl) += 1;
        sp += 1;
    } while (dl[sp].kind != DL_BOTTOM && sp < DAME_MAX_DEPTH);
    (*curr_dl)->kind = DL_BOTTOM;
    (*curr_dl) += 1;
  clean_exit:
    if (allocated) {
        *curr_dl = allocated;
    }
    return 0;
}

/* This takes a dataloop before it has been updated and serializes (compacts)
 * it into a contiguous buffer.
 * This will work even with an updated (structs adjusted for partial pack
 * efficiency) dataloop, but there will be wasted space in the dataloop.
 * The dataloop_size is the size of this compacted dataloop.
 * Returns 0 on success, -1 on failure
 */
int MPIR_Dame_serialize(DAME_Dame * dl, MPI_Aint dataloop_size, DAME_Dame ** compact_dataloop)
{
    void *tmp = NULL;
    return MPIR_Dame_serialize_impl(dl, dataloop_size, &tmp, compact_dataloop);
}

/* The reverse operation of serialize.
 *
 * Returns 0 on success, -1 on failure
 */
int MPIR_Dame_deserialize(DAME_Dame * compact_dl,
                          DAME_Offset dataloop_size,
                          DAME_Dame ** new_compact_dl, DAME_Dame ** new_dl)
{
    int i, j;
    MPI_Aint *compact_offsets = NULL;
    MPI_Aint *compact_blklens = NULL;
    MPI_Aint *compact_oldsizes = NULL;
    DAME_Dame **compact_dls = NULL;

    if (new_compact_dl) {
        if (MPIR_Dame_dup_compact(compact_dl, dataloop_size, new_compact_dl) != 0)
            return -1;
    }

    if (*new_dl == 0) {
        if (MPIR_Dame_alloc(new_dl) == -1)
            return -1;
    }

    for (i = 0; i < DAME_MAX_DEPTH && compact_dl[i].kind != DL_BOTTOM; i++) {
        (*new_dl)[i].kind = compact_dl[i].kind;
        (*new_dl)[i].count = compact_dl[i].count;
        (*new_dl)[i].size = compact_dl[i].size;
        (*new_dl)[i].extent = compact_dl[i].extent;
        (*new_dl)[i].returnto = compact_dl[i].returnto;
        switch (compact_dl[i].kind) {
            case DL_CONTIG:
            case DL_CONTIGFINAL:
                (*new_dl)[i].s.c_t.basesize = compact_dl[i].s.c_t.basesize;
                (*new_dl)[i].s.c_t.baseextent = compact_dl[i].s.c_t.baseextent;
                break;
            case DL_VECTOR:
            case DL_VECTOR1:
            case DL_VECTORFINAL:
                (*new_dl)[i].s.v_t.oldsize = compact_dl[i].s.v_t.oldsize;
                (*new_dl)[i].s.v_t.blklen = compact_dl[i].s.v_t.blklen;
                (*new_dl)[i].s.v_t.stride = compact_dl[i].s.v_t.stride;
                break;
            case DL_BLOCKINDEX:
            case DL_BLOCKINDEX1:
            case DL_BLOCKINDEXFINAL:
                (*new_dl)[i].s.bi_t.oldsize = compact_dl[i].s.bi_t.oldsize;
                (*new_dl)[i].s.bi_t.blklen = compact_dl[i].s.bi_t.blklen;
                (*new_dl)[i].s.bi_t.offsets
                    = (MPI_Aint *) MPL_malloc(compact_dl[i].count * sizeof(MPI_Aint), MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.bi_t.offsets)
                    goto error_exit;
                compact_offsets = GET_COMPACTED_BUF(compact_dl[i].s.bi_t.offsets);
                memcpy((MPI_Aint *) (*new_dl)[i].s.bi_t.offsets,
                       compact_offsets, compact_dl[i].count * sizeof(MPI_Aint));
                break;
            case DL_INDEX:
            case DL_INDEXFINAL:
                (*new_dl)[i].s.i_t.oldsize = compact_dl[i].s.i_t.oldsize;
                (*new_dl)[i].s.i_t.blklens
                    = (MPI_Aint *) MPL_malloc(compact_dl[i].count * sizeof(MPI_Aint), MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.i_t.blklens)
                    goto error_exit;
                compact_blklens = GET_COMPACTED_BUF(compact_dl[i].s.i_t.blklens);
                memcpy((MPI_Aint *) (*new_dl)[i].s.i_t.blklens,
                       compact_blklens, compact_dl[i].count * sizeof(MPI_Aint));
                (*new_dl)[i].s.i_t.offsets
                    = (MPI_Aint *) MPL_malloc(compact_dl[i].count * sizeof(MPI_Aint), MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.i_t.offsets)
                    goto error_exit;
                compact_offsets = GET_COMPACTED_BUF(compact_dl[i].s.i_t.offsets);
                memcpy((MPI_Aint *) (*new_dl)[i].s.i_t.offsets,
                       compact_offsets, compact_dl[i].count * sizeof(MPI_Aint));
                break;
            case DL_STRUCT:
                (*new_dl)[i].s.s_t.oldsizes
                   = (MPI_Aint *) MPL_malloc(compact_dl[i].count * sizeof(MPI_Aint), MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.s_t.oldsizes)
                    goto error_exit;
                compact_oldsizes = GET_COMPACTED_BUF(compact_dl[i].s.s_t.oldsizes);
                memcpy((MPI_Aint *) (*new_dl)[i].s.s_t.oldsizes,
                       compact_oldsizes, compact_dl[i].count * sizeof(MPI_Aint));
                (*new_dl)[i].s.s_t.blklens
                    = (MPI_Aint *) MPL_malloc(compact_dl[i].count * sizeof(MPI_Aint), MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.s_t.blklens)
                    goto error_exit;
                compact_blklens = GET_COMPACTED_BUF(compact_dl[i].s.s_t.blklens);
                memcpy((MPI_Aint *) (*new_dl)[i].s.s_t.blklens,
                       compact_blklens, compact_dl[i].count * sizeof(MPI_Aint));
                (*new_dl)[i].s.s_t.offsets
                    = (MPI_Aint *) MPL_malloc(compact_dl[i].count * sizeof(MPI_Aint), MPL_MEM_DATATYPE);
                if (!(*new_dl)[i].s.s_t.offsets)
                    goto error_exit;
                compact_offsets = GET_COMPACTED_BUF(compact_dl[i].s.s_t.offsets);
                memcpy((MPI_Aint *) (*new_dl)[i].s.s_t.offsets,
                       compact_offsets, compact_dl[i].count * sizeof(MPI_Aint));

                MPIR_Dame_struct_alloc(compact_dl[i].count, &(*new_dl)[i]);
                compact_dls = GET_COMPACTED_BUF(compact_dl[i].s.s_t.dls);
                for (j = 0; j < compact_dl[i].count; j++)
                    MPIR_Dame_deserialize
                        (GET_COMPACTED_BUF(compact_dls[j]), 0, NULL, &(*new_dl)[i].s.s_t.dls[j]);
                break;
            case DL_CONTIGCHILD:
                break;
            case DL_EXIT:
                break;
            case DL_BOTTOM:
                break;
            default:
                break;
        }
    }

    return MPIR_Dame_update(*new_dl, 0);

  error_exit:
    return -1;
}
