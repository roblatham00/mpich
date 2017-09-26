/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h>
#include "veccpy.h"
#include "mpiimpl.h"

#define DBG_PRINTF(s, ...)                      \
        /* #define DBG_PRINTF fprintf */

static MPI_Aint min(MPI_Aint l, MPI_Aint r)
{
    return l < r ? l : r;
}

static void DL_pack_contig(const void *in, void *out, const DAME_Dame * dl, MPI_Aint count)
{
    DAME_Offset basesize = dl->s.c_t.basesize;
    if (count > DAME_MEMCPY_THRESHOLD) {
        DAME_Memcpy(out, in, count);
    } else {
        DAME_Offset full = count / basesize;
        DAME_Offset part = count % basesize;
        switch (basesize) {
            case 8:
                MPIDI_COPY_BLOCK_UNALIGNED(in, out, int64_t, full);
                MPIDI_COPY_BLOCK_UNALIGNED(in, out, int8_t, part);
                break;
            case 4:
                MPIDI_COPY_BLOCK_UNALIGNED(in, out, int32_t, full);
                MPIDI_COPY_BLOCK_UNALIGNED(in, out, int8_t, part);
                break;
            case 2:
                MPIDI_COPY_BLOCK_UNALIGNED(in, out, int16_t, full);
                MPIDI_COPY_BLOCK_UNALIGNED(in, out, int8_t, part);
                break;
            default:
                MPIDI_COPY_BLOCK_UNALIGNED(in, out, int8_t, count);
                break;
        }
    }
}

static void DL_pack_iov_contig(void *in, DAME_VECTOR ** out, const DAME_Dame * dl, MPI_Aint count)
{
    char *restrict src = (char *) in;
    DAME_VECTOR *restrict dst = (DAME_VECTOR *) * out;

    dst->MPL_IOV_BUF = src;
    dst->MPL_IOV_LEN = count;
    dst += 1;
    *out = dst;
}

static void DL_pack_vector(const void *in, void **out, const DAME_Dame * dl, MPI_Aint count)
{
    const char *restrict src = (const char *) in;
    char *restrict dst = (char *) *out;
    MPI_Aint i, j;

    MPI_Aint stride = dl->s.v_t.stride;
    MPI_Aint blklen = dl->s.v_t.blklen;
    MPI_Aint oldsize = dl->s.v_t.oldsize;
    MPI_Aint blksize = blklen * oldsize;

    MPI_Aint oldextent = (dl + 1)->s.c_t.baseextent;

    int aligned = DAME_opt_get_aligned(*dl);
    int isshort = DAME_opt_get_isshort(*dl);

    if (oldsize == oldextent) {
        switch (oldsize) {
            case 8:
                MPIDI_COPY_FROM_VEC_ALIGNED(src, dst, stride / 8, int64_t, blklen, count, 8);
                break;
            case 4:
                MPIDI_COPY_FROM_VEC_ALIGNED(src, dst, stride / 4, int32_t, blklen, count, 4);
                break;
            case 2:
                MPIDI_COPY_FROM_VEC_ALIGNED(src, dst, stride / 2, int16_t, blklen, count, 2);
                break;
            default:
                /* We could simplify the code by having a separate loop over the
                 * blocks, but doing it this way widens the memcpy
                 * in the common case */
                if (isshort) {
                    for (i = 0; i < count; i++) {
                        MPIDI_COPY_BLOCK_UNALIGNED(src, dst, int8_t, blklen * oldsize);
                        dst += (blklen * oldsize);
                        src += stride;
                    }
                } else {
                    for (i = 0; i < count; i++) {
                        DAME_Memcpy(dst, src, blklen * oldsize);
                        dst += (blklen * oldsize);
                        src += stride;
                    }
                }
        }
    } else {
        if (isshort) {
            for (i = 0; i < count; i++) {
                const char *tmpsrc = src;
                for (j = 0; j < blklen; j++) {
                    MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int8_t, oldsize);
                    dst += oldsize;
                    tmpsrc += oldextent;
                }
                src += stride;
            }
        } else {
            for (i = 0; i < count; i++) {
                const char *tmpsrc = src;
                for (j = 0; j < blklen; j++) {
                    DAME_Memcpy(dst, tmpsrc, oldsize);
                    dst += oldsize;
                    tmpsrc += oldextent;
                }
                src += stride;
            }
        }
    }
    *out = (void *) ((char *) *out + count * blklen * oldsize);
}

static void DL_pack_iov_vector(const void *in,
                               DAME_VECTOR ** out, const DAME_Dame * dl, MPI_Aint count)
{
    char *restrict src = (char *) in;
    DAME_VECTOR *restrict dst = (DAME_VECTOR *) (*out);
    MPI_Aint i, j;

    MPI_Aint stride = dl->s.v_t.stride;
    MPI_Aint blklen = dl->s.v_t.blklen;
    MPI_Aint oldsize = dl->s.v_t.oldsize;
    MPI_Aint blksize = blklen * oldsize;

    MPI_Aint oldextent = (dl + 1)->s.c_t.baseextent;

    if (oldsize == oldextent) {
        for (i = 0; i < count; i++) {
            dst->MPL_IOV_BUF = src;
            dst->MPL_IOV_LEN = blklen * oldsize;
            dst++;
            src += stride;
        }
    } else {
        for (i = 0; i < count; i++) {
            char *tmpsrc = src;
            for (j = 0; j < blklen; j++) {
                dst->MPL_IOV_BUF = tmpsrc;
                dst->MPL_IOV_LEN = oldsize;
                dst++;
                tmpsrc += oldextent;
            }
            src += stride;
        }
    }
    *out += count;
}

static void DL_pack_blockindexed(const void *in,
                                 void **out, const DAME_Dame * dl, MPI_Aint begin, MPI_Aint end)
{
    const char *restrict src = (const char *) in;
    char *restrict dst = (char *) *out;
    MPI_Aint i, j;

    MPI_Aint blklen = dl->s.bi_t.blklen;
    MPI_Aint oldsize = dl->s.bi_t.oldsize;
    const MPI_Aint *offsets = dl->s.bi_t.offsets;
    MPI_Aint blksize = blklen * oldsize;

    MPI_Aint oldextent = (dl + 1)->s.c_t.baseextent;

    int aligned = DAME_opt_get_aligned(*dl);
    int isshort = DAME_opt_get_isshort(*dl);

    const char *restrict tmpsrc = 0;
    if (oldsize == oldextent) {
        switch (oldsize) {
            case 8:
                for (i = begin; i < end; i++) {
                    tmpsrc = src + offsets[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int64_t, blklen);
                    dst += blklen * oldsize;
                }
                break;
            case 4:
                for (i = begin; i < end; i++) {
                    tmpsrc = src + offsets[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int32_t, blklen);
                    dst += blklen * oldsize;
                }
                break;
            case 2:
                for (i = begin; i < end; i++) {
                    tmpsrc = src + offsets[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int16_t, blklen);
                    dst += blklen * oldsize;
                }
                break;
            default:
                /* Splitting this into two to make the memcpy as wider for
                 * the common case */
                if (isshort) {
                    for (i = begin; i < end; i++) {
                        tmpsrc = src + offsets[i];
                        MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int8_t, blklen * oldsize);
                        dst += (blklen * oldsize);
                    }
                } else {
                    for (i = begin; i < end; i++) {
                        tmpsrc = src + offsets[i];
                        DAME_Memcpy(dst, tmpsrc, blklen * oldsize);
                        dst += (blklen * oldsize);
                    }
                }
                break;
        }
    } else {
        if (isshort) {
            for (i = begin; i < end; i++) {
                tmpsrc = src + offsets[i];
                for (j = 0; j < blklen; j++) {
                    MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int8_t, oldsize);
                    dst += oldsize;
                    tmpsrc += oldextent;
                }
            }
        } else {
            for (i = begin; i < end; i++) {
                tmpsrc = src + offsets[i];
                for (j = 0; j < blklen; j++) {
                    DAME_Memcpy(dst, tmpsrc, oldsize);
                    dst += oldsize;
                    tmpsrc += oldextent;
                }
            }
        }
    }
    *out = (void *) dst;
}

static void DL_pack_iov_blockindexed(void *in,
                                     DAME_VECTOR ** out,
                                     const DAME_Dame * dl, MPI_Aint begin, MPI_Aint end)
{
    char *restrict src = (char *) in;
    DAME_VECTOR *restrict dst = (DAME_VECTOR *) * out;
    MPI_Aint i, j;

    MPI_Aint blklen = dl->s.bi_t.blklen;
    MPI_Aint oldsize = dl->s.bi_t.oldsize;
    const MPI_Aint *offsets = dl->s.bi_t.offsets;
    MPI_Aint blksize = blklen * oldsize;

    MPI_Aint oldextent = (dl + 1)->s.c_t.baseextent;

    char *restrict tmpsrc = 0;
    if (oldsize == oldextent) {
        for (i = begin; i < end; i++) {
            tmpsrc = src + offsets[i];
            dst->MPL_IOV_BUF = tmpsrc;
            dst->MPL_IOV_LEN = blklen * oldsize;
            dst += 1;
        }
    } else {
        for (i = begin; i < end; i++) {
            tmpsrc = src + offsets[i];
            for (j = 0; j < blklen; j++) {
                dst->MPL_IOV_BUF = tmpsrc;
                dst->MPL_IOV_LEN = oldsize;
                tmpsrc += oldextent;
                dst += 1;
            }
        }
    }

    *out = dst;
}

static void DL_pack_indexed(const void *in,
                            void **out, const DAME_Dame * dl, MPI_Aint begin, MPI_Aint end)
{
    const char *restrict src = (const char *) in;
    char *restrict dst = (char *) *out;
    MPI_Aint i, j;

    MPI_Aint blklen;
    MPI_Aint oldsize = dl->s.i_t.oldsize;
    const MPI_Aint *blklens = dl->s.i_t.blklens;
    const MPI_Aint *offsets = dl->s.i_t.offsets;

    MPI_Aint oldextent = (dl + 1)->s.c_t.baseextent;

    int aligned = DAME_opt_get_aligned(*dl);
    int isshort = DAME_opt_get_isshort(*dl);

    const char *restrict tmpsrc = 0;
    if (oldsize == oldextent) {
        switch (oldsize) {
            case 8:
                for (i = begin; i < end; i++) {
                    tmpsrc = src + offsets[i];
                    blklen = blklens[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int64_t, blklen);
                    dst += blklen * oldsize;
                }
                break;
            case 4:
                for (i = begin; i < end; i++) {
                    tmpsrc = src + offsets[i];
                    blklen = blklens[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int32_t, blklen);
                    dst += blklen * oldsize;
                }
                break;
            case 2:
                for (i = begin; i < end; i++) {
                    tmpsrc = src + offsets[i];
                    blklen = blklens[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int16_t, blklen);
                    dst += blklen * oldsize;
                }
                break;
            default:
                if (isshort) {
                    for (i = begin; i < end; i++) {
                        tmpsrc = src + offsets[i];
                        blklen = blklens[i];
                        MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int8_t, blklen * oldsize);
                        dst += (blklen * oldsize);
                    }
                } else {
                    for (i = begin; i < end; i++) {
                        tmpsrc = src + offsets[i];
                        blklen = blklens[i];
                        DAME_Memcpy(dst, tmpsrc, blklen * oldsize);
                        dst += (blklen * oldsize);
                    }
                }
                break;
        }
    } else {
        if (isshort) {
            for (i = begin; i < end; i++) {
                tmpsrc = src + offsets[i];
                for (j = 0; j < blklens[i]; j++) {
                    MPIDI_COPY_BLOCK_UNALIGNED(tmpsrc, dst, int8_t, oldsize);
                    dst += oldsize;
                    tmpsrc += oldextent;
                }
            }
        } else {
            for (i = begin; i < end; i++) {
                tmpsrc = src + offsets[i];
                for (j = 0; j < blklens[i]; j++) {
                    DAME_Memcpy(dst, tmpsrc, oldsize);
                    dst += oldsize;
                    tmpsrc += oldextent;
                }
            }
        }
    }
    *out = (void *) dst;
}

static void DL_pack_iov_indexed(void *in,
                                DAME_VECTOR ** out,
                                const DAME_Dame * dl, MPI_Aint begin, MPI_Aint end)
{
    char *restrict src = (char *) in;
    DAME_VECTOR *restrict dst = (DAME_VECTOR *) * out;
    MPI_Aint i, j;

    MPI_Aint blklen;
    MPI_Aint oldsize = dl->s.i_t.oldsize;
    const MPI_Aint *blklens = dl->s.i_t.blklens;
    const MPI_Aint *offsets = dl->s.i_t.offsets;

    MPI_Aint oldextent = (dl + 1)->s.c_t.baseextent;

    char *tmpsrc = 0;
    if (oldsize == oldextent) {
        for (i = begin; i < end; i++) {
            tmpsrc = src + offsets[i];
            dst->MPL_IOV_BUF = tmpsrc;
            dst->MPL_IOV_LEN = blklens[i] * oldsize;
            dst += 1;
        }
    } else {
        for (i = begin; i < end; i++) {
            tmpsrc = src + offsets[i];
            for (j = 0; j < blklens[i]; j++) {
                dst->MPL_IOV_BUF = tmpsrc;
                dst->MPL_IOV_LEN = oldsize;
                dst += 1;
                tmpsrc += oldextent;
            }
        }
    }

    *out = dst;
}

static void DL_unpack_vector(const void **in, void *out, const DAME_Dame * dl, MPI_Aint count)
{
    const char *restrict src = (const char *) *in;
    char *restrict dst = (char *) out;
    MPI_Aint i, j;

    MPI_Aint stride = dl->s.v_t.stride;
    MPI_Aint blklen = dl->s.v_t.blklen;
    MPI_Aint oldsize = dl->s.v_t.oldsize;

    MPI_Aint oldextent = (dl + 1)->s.c_t.baseextent;

    int aligned = DAME_opt_get_aligned(*dl);
    int isshort = DAME_opt_get_isshort(*dl);

    if (oldsize == oldextent) {
        switch (oldsize) {
            case 8:
                MPIDI_COPY_TO_VEC_ALIGNED(src, dst, stride / 8, int64_t, blklen, count, 8);
                break;
            case 4:
                MPIDI_COPY_TO_VEC_ALIGNED(src, dst, stride / 4, int32_t, blklen, count, 4);
                break;
            case 2:
                MPIDI_COPY_TO_VEC_ALIGNED(src, dst, stride / 2, int16_t, blklen, count, 2);
                break;
            default:
                if (isshort) {
                    for (i = 0; i < count; i++) {
                        MPIDI_COPY_BLOCK_UNALIGNED(src, dst, int8_t, blklen * oldsize);
                        src += (blklen * oldsize);
                        dst += stride;
                    }
                } else {
                    for (i = 0; i < count; i++) {
                        DAME_Memcpy(dst, src, blklen * oldsize);
                        src += (blklen * oldsize);
                        dst += stride;
                    }
                }
                break;
        }
    } else {
        if (isshort) {
            for (i = 0; i < count; i++) {
                char *tmpdst = dst;
                for (j = 0; j < blklen; j++) {
                    MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int8_t, oldsize);
                    src += oldsize;
                    tmpdst += oldextent;
                }
                dst += stride;
            }
        } else {
            for (i = 0; i < count; i++) {
                char *tmpdst = dst;
                for (j = 0; j < blklen; j++) {
                    DAME_Memcpy(tmpdst, src, oldsize);
                    src += oldsize;
                    tmpdst += oldextent;
                }
                dst += stride;
            }
        }
    }
    *in = (void *) ((char *) *in + count * blklen * oldsize);
}

static void DL_unpack_blockindexed(const void **in,
                                   void *out, const DAME_Dame * dl, MPI_Aint begin, MPI_Aint end)
{
    const char *restrict src = (const char *) *in;
    char *restrict dst = (char *) out;
    MPI_Aint i, j;

    MPI_Aint blklen = dl->s.bi_t.blklen;
    MPI_Aint oldsize = dl->s.bi_t.oldsize;
    const MPI_Aint *offsets = dl->s.bi_t.offsets;

    MPI_Aint oldextent = (dl + 1)->s.c_t.baseextent;

    int aligned = DAME_opt_get_aligned(*dl);
    int isshort = DAME_opt_get_isshort(*dl);

    char *restrict tmpdst = 0;
    if (oldsize == oldextent) {
        switch (oldsize) {
            case 8:
                for (i = begin; i < end; i++) {
                    tmpdst = dst + offsets[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int64_t, blklen);
                    src += blklen * oldsize;
                }
                break;
            case 4:
                for (i = begin; i < end; i++) {
                    tmpdst = dst + offsets[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int32_t, blklen);
                    src += blklen * oldsize;
                }
                break;
            case 2:
                for (i = begin; i < end; i++) {
                    tmpdst = dst + offsets[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int16_t, blklen);
                    src += blklen * oldsize;
                }
                break;
            default:
                if (isshort) {
                    for (i = begin; i < end; i++) {
                        tmpdst = dst + offsets[i];
                        MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int8_t, blklen * oldsize);
                        src += (blklen * oldsize);
                    }
                } else {
                    for (i = begin; i < end; i++) {
                        tmpdst = dst + offsets[i];
                        DAME_Memcpy(tmpdst, src, blklen * oldsize);
                        src += (blklen * oldsize);
                    }
                }
                break;
        }
    } else {
        if (isshort) {
            for (i = begin; i < end; i++) {
                tmpdst = dst + offsets[i];
                for (j = 0; j < blklen; j++) {
                    MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int8_t, oldsize);
                    src += oldsize;
                    tmpdst += oldextent;
                }
            }
        } else {
            for (i = begin; i < end; i++) {
                tmpdst = dst + offsets[i];
                for (j = 0; j < blklen; j++) {
                    DAME_Memcpy(tmpdst, src, oldsize);
                    src += oldsize;
                    tmpdst += oldextent;
                }
            }
        }
    }
    *in = (void *) src;
}

static void DL_unpack_indexed(const void **in,
                              void *out, const DAME_Dame * dl, MPI_Aint begin, MPI_Aint end)
{
    const char *restrict src = (const char *) *in;
    char *restrict dst = (char *) out;
    MPI_Aint i, j;

    MPI_Aint blklen;
    MPI_Aint oldsize = dl->s.i_t.oldsize;
    const MPI_Aint *blklens = dl->s.i_t.blklens;
    const MPI_Aint *offsets = dl->s.i_t.offsets;

    MPI_Aint oldextent = (dl + 1)->s.c_t.baseextent;

    int aligned = DAME_opt_get_aligned(*dl);
    int isshort = DAME_opt_get_isshort(*dl);

    char *restrict tmpdst = 0;
    if (oldsize == oldextent) {
        switch (oldsize) {
            case 8:
                for (i = begin; i < end; i++) {
                    tmpdst = dst + offsets[i];
                    blklen = blklens[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int64_t, blklen);
                    src += blklen * oldsize;
                }
                break;
            case 4:
                for (i = begin; i < end; i++) {
                    tmpdst = dst + offsets[i];
                    blklen = blklens[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int32_t, blklen);
                    src += blklen * oldsize;
                }
                break;
            case 2:
                for (i = begin; i < end; i++) {
                    tmpdst = dst + offsets[i];
                    blklen = blklens[i];
                    MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int16_t, blklen);
                    src += blklen * oldsize;
                }
                break;
            default:
                if (isshort) {
                    for (i = begin; i < end; i++) {
                        tmpdst = dst + offsets[i];
                        blklen = blklens[i];
                        MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int8_t, blklen * oldsize);
                        src += blklen * oldsize;
                    }
                } else {
                    for (i = begin; i < end; i++) {
                        tmpdst = dst + offsets[i];
                        blklen = blklens[i];
                        DAME_Memcpy(tmpdst, src, blklen * oldsize);
                        src += blklen * oldsize;
                    }
                }
                break;
        }
    } else {
        if (isshort) {
            for (i = begin; i < end; i++) {
                tmpdst = dst + offsets[i];
                for (j = 0; j < blklens[i]; j++) {
                    MPIDI_COPY_BLOCK_UNALIGNED(src, tmpdst, int8_t, oldsize);
                    src += oldsize;
                    tmpdst += oldextent;
                }
            }
        } else {
            for (i = begin; i < end; i++) {
                tmpdst = dst + offsets[i];
                for (j = 0; j < blklens[i]; j++) {
                    DAME_Memcpy(tmpdst, src, oldsize);
                    src += oldsize;
                    tmpdst += oldextent;
                }
            }
        }
    }
    *in = (void *) src;
}


static int shiftIOVs(MPI_Aint count, DAME_VECTOR * curr)
{
    int j;
    DAME_VECTOR *prev = curr - 1;
    if (prev->MPL_IOV_BUF + prev->MPL_IOV_LEN == curr->MPL_IOV_BUF) {
        prev->MPL_IOV_LEN += curr->MPL_IOV_LEN;
        for (j = 0; j < count - 1; j++) {
            curr[j].MPL_IOV_BUF = curr[j + 1].MPL_IOV_BUF;
            curr[j].MPL_IOV_LEN = curr[j + 1].MPL_IOV_LEN;
        }
        return 1;
    }
    return 0;
}


static void dlStackCopy(long sp, const DAME_Dame * dl, MPI_Aint partial, DAME_State * state)
{
    state->sp = sp;
    state->partial = partial;
    state->dl = dl;
}


/*
  Return values:
  -1: Error
  0: Done
  1: Needs to continue; state in state_p
*/
int DL_pack_iov(const void *inPtr,
                DAME_VECTOR * outPtr,
                DAME_Offset outSize,
                int outSlots,
                int merge,
                const DAME_Dame * dl,
                int isbuiltin, DAME_State * state, DAME_Offset * copy_p, int *iovs_p)
{
    MPI_Offset inbase = (MPI_Offset) inPtr;
    MPI_Offset outbase = (MPI_Offset) outPtr;
    long sp = 0;
    MPI_Aint sizeleft = outSize;
    int slotsleft = outSlots;

    /* Local variables */
    DAME_Stackelm *stack = state->stack;
    /* If it's a resume, this will get set correctly later */
    DAME_Offset partial = 0;

    *copy_p = 0;

    /* This is a resume */
    if (state->sp != 0) {
        sp = state->sp;
        inbase = state->stack[sp].base;
        partial = state->partial;

        /* NOTE: For builtin types, we create a dataloop as needed on the stack
         * so we don't save it in the state because when we resume, the same
         * type will be recreated in segment_packunpack. */
        if (!isbuiltin)
            dl = state->dl;

        /* By construction, the last state pushed into the stack will always
         * be a CONTIGFINAL. Need to decrement the stack pointer to ensure
         * that it is set to the right place by the immediate next call to
         * dlpush */
        sp--;
        goto dlpush;
    }

    do {
      dlpush:
        sp++;
        switch (dl[sp].kind) {
            case DL_CONTIG:{
                    DBG_PRINTF(stderr, "contig%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    goto dlpush;
                }
            case DL_VECTOR:{
                    DBG_PRINTF(stderr, "vec%d_push(%p, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, slotsleft);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp + 1].base = inbase;
                    stack[sp + 1].countLeft = dl[sp].s.v_t.blklen;
                    goto dlpush;
                    break;
                }
            case DL_VECTOR1:{
                    DBG_PRINTF(stderr, "vec%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX:{
                    DBG_PRINTF(stderr, "blkidx%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp + 1].base = inbase;
                    stack[sp + 1].countLeft = dl[sp].s.bi_t.blklen;
                    inbase += dl[sp].s.bi_t.offsets[0];
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX1:{
                    DBG_PRINTF(stderr, "blkidx%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    inbase += dl[sp].s.bi_t.offsets[0];
                    goto dlpush;
                    break;
                }
            case DL_INDEX:{
                    DBG_PRINTF(stderr, "idx%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp + 1].base = inbase;
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[0];
                    inbase += dl[sp].s.i_t.offsets[0];
                    goto dlpush;
                    break;
                }
            case DL_STRUCT:{
                    DBG_PRINTF(stderr, "struct%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp].prevdl = dl;
                    inbase += dl[sp].s.s_t.offsets[0];
                    dl = dl[sp].s.s_t.dls[0];
                    goto dlpush;
                    break;
                }
            case DL_CONTIGCHILD:{
                    DBG_PRINTF(stderr, "cc%d_push\n", sp);
                    /* countLeft set by parent */
                    /* inbase set by parent */
                    goto dlpush;
                }
            case DL_VECTORFINAL:{
                    DBG_PRINTF(stderr, "vf%d_push(%p, %ld, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft, slotsleft);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;

                    MPI_Aint blklen = dl[sp].s.v_t.blklen;
                    MPI_Aint oldsize = dl[sp].s.v_t.oldsize;
                    MPI_Aint stride = dl[sp].s.v_t.stride;

                    MPI_Aint oldextent = dl[sp + 1].s.c_t.baseextent;

                    DAME_VECTOR *curr = (DAME_VECTOR *) outbase;
                    MPI_Aint numblocks = dl[sp].count;
                    if (oldsize != oldextent)
                        numblocks *= blklen;
                    if (sizeleft >= dl[sp].size && slotsleft >= numblocks) {
                        if (outbase) {
                            DL_pack_iov_vector((void *) inbase,
                                               (DAME_VECTOR **) & outbase, &dl[sp], dl[sp].count);
                            if (merge || (outSlots != slotsleft)) {
                                if (shiftIOVs(dl[sp].count, curr)) {
                                    slotsleft += 1;
                                    outbase -= sizeof(DAME_VECTOR);
                                }
                            }
                        }
                        stack[sp].countLeft = 0;
                        sizeleft -= dl[sp].size;
                        slotsleft -= numblocks;
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                        break;
                    }
                    /* destination too small.  Determine how far we can go */
                    DAME_Offset count = min(sizeleft / (blklen * oldsize), slotsleft);
                    numblocks = count;
                    if (oldsize != oldextent) {
                        numblocks = min(count * blklen, slotsleft);
                        count = numblocks / blklen;
                    }
                    if (count > 0) {
                        if (outbase) {
                            DL_pack_iov_vector((void *) inbase,
                                               (DAME_VECTOR **) & outbase, &dl[sp], count);
                            if (merge || (outSlots != slotsleft)) {
                                if (shiftIOVs(count, curr)) {
                                    slotsleft += 1;
                                    outbase -= sizeof(DAME_VECTOR);
                                }
                            }
                        }
                        sizeleft -= count * blklen * oldsize;
                        slotsleft -= numblocks;
                        stack[sp].countLeft -= count;
                    }
                    DAME_Offset idx = count;
                    stack[sp + 1].countLeft = blklen;
                    inbase = stack[sp].base + idx * stride;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEXFINAL:{
                    DBG_PRINTF(stderr, "bf%d_push(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;

                    MPI_Aint blklen = dl[sp].s.bi_t.blklen;
                    MPI_Aint oldsize = dl[sp].s.bi_t.oldsize;

                    MPI_Aint oldextent = dl[sp + 1].s.c_t.baseextent;

                    DAME_VECTOR *curr = (DAME_VECTOR *) outbase;
                    MPI_Aint numblocks = dl[sp].count;
                    if (oldsize != oldextent)
                        numblocks *= blklen;

                    if (sizeleft >= dl[sp].size && slotsleft >= numblocks) {
                        if (outbase) {
                            DL_pack_iov_blockindexed((void *) stack[sp].base,
                                                     (DAME_VECTOR **) & outbase,
                                                     &dl[sp], 0, dl[sp].count);
                            if (merge || (outSlots != slotsleft)) {
                                if (shiftIOVs(dl[sp].count, curr)) {
                                    slotsleft += 1;
                                    outbase -= sizeof(DAME_VECTOR);
                                }
                            }
                        }
                        stack[sp].countLeft = 0;
                        sizeleft -= dl[sp].size;
                        slotsleft -= numblocks;
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                        break;
                    }
                    /* destination too small.  Determine how far we can go */
                    DAME_Offset count = min(sizeleft / (blklen * oldsize), slotsleft);
                    numblocks = count;
                    if (oldsize != oldextent) {
                        numblocks = min(count * blklen, slotsleft);
                        count = numblocks / blklen;
                    }
                    if (count > 0) {
                        if (outbase)
                            DL_pack_iov_blockindexed((void *) stack[sp].base,
                                                     (DAME_VECTOR **) & outbase, &dl[sp], 0, count);
                        if (merge || (outSlots != slotsleft)) {
                            if (shiftIOVs(count, curr)) {
                                slotsleft += 1;
                                outbase -= sizeof(DAME_VECTOR);
                            }
                        }
                        sizeleft -= count * blklen * oldsize;
                        slotsleft -= numblocks;
                        stack[sp].countLeft -= count;
                    }
                    DAME_Offset idx = count;
                    stack[sp + 1].countLeft = blklen;
                    inbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_INDEXFINAL:{
                    DBG_PRINTF(stderr, "if%d_push(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;

                    MPI_Aint oldsize = dl[sp].s.i_t.oldsize;

                    MPI_Aint oldextent = dl[sp + 1].s.c_t.baseextent;

                    DAME_VECTOR *curr = (DAME_VECTOR *) outbase;
                    MPI_Aint numblocks = dl[sp].count;
                    if (oldsize != oldextent) {
                        int i;
                        for (i = 0; i < dl[sp].count; i++)
                            numblocks += dl[sp].s.i_t.blklens[i];
                    }
                    if (sizeleft >= dl[sp].size && slotsleft >= numblocks) {
                        if (outbase) {
                            DL_pack_iov_indexed((void *) stack[sp].base,
                                                (DAME_VECTOR **) & outbase,
                                                &dl[sp], 0, dl[sp].count);
                            if (merge || (outSlots != slotsleft)) {
                                if (shiftIOVs(dl[sp].count, curr)) {
                                    slotsleft += 1;
                                    outbase -= sizeof(DAME_VECTOR);
                                }
                            }
                        }
                        stack[sp].countLeft = 0;
                        sizeleft -= dl[sp].size;
                        slotsleft -= numblocks;
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }

                    /* destination too small.  Determine how far we can go */
                    DAME_Offset count, copied;
                    if (oldsize != oldextent) {
                        for (count = 0, copied = 0, numblocks = 0;
                             count < dl[sp].count
                             && numblocks + dl[sp].s.i_t.blklens[count] < slotsleft
                             && copied + dl[sp].s.i_t.blklens[count] * oldsize <= sizeleft;
                             count++) {
                            copied += dl[sp].s.i_t.blklens[count] * oldsize;
                            numblocks += dl[sp].s.i_t.blklens[count];
                        }
                    } else {
                        for (count = 0, copied = 0, numblocks = 0;
                             count < dl[sp].count && count < slotsleft
                             && copied + dl[sp].s.i_t.blklens[count] * oldsize <= sizeleft;
                             count++, numblocks++) {
                            copied += dl[sp].s.i_t.blklens[count] * oldsize;
                        }
                    }

                    if (count > 0) {
                        if (outbase) {
                            DL_pack_iov_indexed((void *) stack[sp].base,
                                                (DAME_VECTOR **) & outbase, &dl[sp], 0, count);
                            if (merge || (outSlots != slotsleft)) {
                                if (shiftIOVs(count, curr)) {
                                    slotsleft += 1;
                                    outbase -= sizeof(DAME_VECTOR);
                                }
                            }
                        }
                        sizeleft -= copied;
                        slotsleft -= count;
                        stack[sp].countLeft -= count;
                    }
                    DAME_Offset idx = count;
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[idx];
                    inbase = stack[sp].base + dl[sp].s.i_t.offsets[idx];
                    goto dlpush;
                    break;
                }

            case DL_CONTIGFINAL:{
                    DBG_PRINTF(stderr, "cf%d_push(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);

                    stack[sp].base = inbase;
                    DAME_Offset basesize = dl[sp].s.c_t.basesize;

                    /* Quick check to see if there's enough buffer space to copy.
                     * sizeleft can be zero if the previous *FINAL type has blksize
                     * which is a multiple of in/outsize */
                    if (sizeleft == 0 || slotsleft == 0) {
                        dlStackCopy(sp, dl, 0, state);
                        *copy_p = outSize - sizeleft;
                        *iovs_p = outSlots - slotsleft;
                        return 1;
                    }

                    /* This case occurs when the contigfinal is directly below a
                     * contig which is part of a struct. There ought to be a better
                     * way of doing this, and eventually,
                     * I'll try and get rid of this check */
                    if (stack[sp].countLeft == 0 && partial == 0)
                        stack[sp].countLeft = dl[sp].count;

                    /* partial contains the remainder of an element that has been
                     * partially packed
                     * countLeft contains the number of complete blocks left */
                    DAME_VECTOR *curr = (DAME_VECTOR *) outbase;
                    DAME_Offset tocopy = stack[sp].countLeft * basesize + partial;
                    if (tocopy > sizeleft) {
                        if (outbase) {
                            DL_pack_iov_contig((void *) inbase,
                                               (DAME_VECTOR **) & outbase, &dl[sp], sizeleft);
                            if (merge || (outSlots != slotsleft)) {
                                if (shiftIOVs(1, curr)) {
                                    slotsleft += 1;
                                    outbase -= sizeof(DAME_VECTOR);
                                }
                            }
                        }
                        DAME_Offset newLeft = stack[sp].countLeft;
                        if (partial >= sizeleft) {
                            partial -= sizeleft;
                        } else {
                            DAME_Offset avail = sizeleft - partial;
                            newLeft -= (avail / basesize);
                            partial = 0;
                            if (avail % basesize) {
                                newLeft -= 1;
                                partial = basesize - avail % basesize;
                            }
                        }
                        stack[sp].base = inbase + sizeleft;
                        stack[sp].countLeft = newLeft;
                        dlStackCopy(sp, dl, partial, state);
                        slotsleft -= 1;
                        sizeleft -= sizeleft;
                        *copy_p = outSize - sizeleft;
                        *iovs_p = outSlots - slotsleft;
                        return 1;
                    } else {
                        if (outbase) {
                            DL_pack_iov_contig((void *) inbase,
                                               (DAME_VECTOR **) & outbase, &dl[sp], tocopy);
                            if (merge || (outSlots != slotsleft)) {
                                if (shiftIOVs(1, curr)) {
                                    slotsleft += 1;
                                    outbase -= sizeof(DAME_VECTOR);
                                }
                            }
                        }
                        partial = 0;
                        inbase += tocopy;
                        sizeleft -= tocopy;
                        slotsleft -= 1;
                        stack[sp].countLeft = 0;
                        goto dlpop;
                    }
                    break;
                }
            default:{
                    /* This can never happen */
                    DBG_PRINTF(stderr, "Unknown dataloop kind(%d)!\n", dl[sp].kind);
                    MPIR_Assert(0);
                    return -1;
                }
        }
      dlpop:
        sp--;

        switch (dl[sp].kind) {
            case DL_CONTIG:{
                    DBG_PRINTF(stderr, "contig%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + idx * dl[sp].s.c_t.baseextent;
                    goto dlpush;
                    break;
                }
            case DL_VECTOR:{
                    DBG_PRINTF(stderr, "vec%d_pop(%p, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, slotsleft);
                    /* For vectors, we don't need to compute offsets from the base.
                     * The stride is the distance from each element, so we can use
                     * the stack contents to update */
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + idx * dl[sp].s.v_t.stride;
                    stack[sp + 1].countLeft = dl[sp].s.v_t.blklen;
                    goto dlpush;
                    break;
                }
            case DL_VECTOR1:{
                    DBG_PRINTF(stderr, "vec%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + idx * dl[sp].s.v_t.stride;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX:{
                    DBG_PRINTF(stderr, "blkidx%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    stack[sp + 1].countLeft = dl[sp].s.bi_t.blklen;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX1:{
                    DBG_PRINTF(stderr, "blkidx%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_INDEX:{
                    DBG_PRINTF(stderr, "idx%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + dl[sp].s.i_t.offsets[idx];
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[idx];
                    goto dlpush;
                    break;
                }
            case DL_STRUCT:{
                    DBG_PRINTF(stderr, "struct%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0) {
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + dl[sp].s.s_t.offsets[idx];
                    dl = dl[sp].s.s_t.dls[idx];
                    goto dlpush;
                    break;
                }
            case DL_CONTIGCHILD:{
                    DBG_PRINTF(stderr, "cc%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0) {
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }
                    goto dlpush;
                    break;
                }
            case DL_VECTORFINAL:{
                    DBG_PRINTF(stderr, "vf%d_pop(%p, %ld, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft, slotsleft);

                    MPI_Aint blklen = dl[sp].s.v_t.blklen;
                    MPI_Aint oldsize = dl[sp].s.v_t.oldsize;
                    MPI_Aint stride = dl[sp].s.v_t.stride;

                    MPI_Aint oldextent = dl[sp + 1].s.c_t.baseextent;

                    DAME_VECTOR *curr = (DAME_VECTOR *) outbase;

                    if (--stack[sp].countLeft == 0) {
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    DAME_Offset count = min(sizeleft / (blklen * oldsize), slotsleft);
                    MPI_Aint numblocks = count;
                    if (oldsize != oldextent) {
                        numblocks = min(count * blklen, slotsleft);
                        count = numblocks / blklen;
                    }
                    if (count > 0) {
                        if (count > stack[sp].countLeft) {
                            if (oldsize != oldextent) {
                                count = min(stack[sp].countLeft * blklen, slotsleft);
                                numblocks = min(count * blklen, slotsleft);
                                count = numblocks / blklen;
                            } else {
                                count = min(stack[sp].countLeft, slotsleft);
                                numblocks = count;
                            }
                        }
                        if (outbase) {
                            DL_pack_iov_vector((void *) (stack[sp].base + idx * stride),
                                               (DAME_VECTOR **) & outbase, &dl[sp], count);
                            if (merge || (outSlots != slotsleft)) {
                                if (shiftIOVs(count, curr)) {
                                    slotsleft += 1;
                                    outbase -= sizeof(DAME_VECTOR);
                                }
                            }
                        }
                        sizeleft -= count * blklen * oldsize;
                        slotsleft -= numblocks;
                        stack[sp].countLeft -= count;
                        if (stack[sp].countLeft == 0) {
                            inbase = stack[sp].base + dl[sp].extent;
                            goto dlpop;
                        }
                    }
                    idx += count;
                    stack[sp + 1].countLeft = blklen;
                    inbase = stack[sp].base + idx * stride;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEXFINAL:{
                    DBG_PRINTF(stderr, "bf%d_pop(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);
                    if (--stack[sp].countLeft == 0) {
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }

                    MPI_Aint blklen = dl[sp].s.bi_t.blklen;
                    MPI_Aint oldsize = dl[sp].s.bi_t.oldsize;

                    MPI_Aint oldextent = dl[sp + 1].s.c_t.baseextent;

                    DAME_VECTOR *curr = (DAME_VECTOR *) outbase;

                    /* If there's still blocks left, try to do several of them */
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    DAME_Offset count = min(sizeleft / (blklen * oldsize), slotsleft);
                    MPI_Aint numblocks = count;
                    if (oldsize != oldextent) {
                        numblocks = min(count * blklen, slotsleft);
                        count = numblocks / blklen;
                    }
                    if (count > 0) {
                        if (count > stack[sp].countLeft) {
                            if (oldsize != oldextent) {
                                count = min(stack[sp].countLeft * blklen, slotsleft);
                                numblocks = min(count * blklen, slotsleft);
                                count = numblocks / blklen;
                            } else {
                                count = min(stack[sp].countLeft, slotsleft);
                                numblocks = count;
                            }
                        }
                        if (outbase) {
                            DL_pack_iov_blockindexed((void *) stack[sp].base,
                                                     (DAME_VECTOR **) & outbase,
                                                     &dl[sp], idx, idx + count);
                            if (merge || (outSlots != slotsleft)) {
                                if (shiftIOVs(count, curr)) {
                                    slotsleft += 1;
                                    outbase -= sizeof(DAME_VECTOR);
                                }
                            }
                        }
                        sizeleft -= count * blklen * oldsize;
                        slotsleft -= numblocks;
                        stack[sp].countLeft -= count;
                        if (stack[sp].countLeft == 0) {
                            inbase = stack[sp].base + dl[sp].extent;
                            goto dlpop;
                        }
                    }
                    idx += count;
                    stack[sp + 1].countLeft = blklen;
                    inbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_INDEXFINAL:{
                    DBG_PRINTF(stderr, "if%d_pop(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);
                    if (--stack[sp].countLeft == 0) {
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }

                    MPI_Aint oldsize = dl[sp].s.i_t.oldsize;

                    MPI_Aint oldextent = dl[sp + 1].s.c_t.baseextent;

                    DAME_VECTOR *curr = (DAME_VECTOR *) outbase;

                    /* If there are still blocks left, try to do several at once */
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    DAME_Offset count, copied, numblocks;
                    if (oldsize != oldextent) {
                        for (count = idx, copied = 0, numblocks = 0;
                             count < dl[sp].count
                             && numblocks + dl[sp].s.i_t.blklens[count] <= slotsleft
                             && copied + dl[sp].s.i_t.blklens[count] * oldsize <= sizeleft;
                             count++) {
                            copied += dl[sp].s.i_t.blklens[count] * oldsize;
                            numblocks += dl[sp].s.i_t.blklens[count];
                        }
                    } else {
                        for (count = idx, copied = 0, numblocks = 0;
                             count < dl[sp].count
                             && count < slotsleft
                             && copied + dl[sp].s.i_t.blklens[count] * oldsize <= sizeleft;
                             count++, numblocks++) {
                            copied += dl[sp].s.i_t.blklens[count] * oldsize;
                        }
                    }
                    if (count > idx) {
                        if (outbase) {
                            DL_pack_iov_indexed((void *) stack[sp].base,
                                                (DAME_VECTOR **) & outbase, &dl[sp], idx, count);
                            if (merge || (outSlots != slotsleft)) {
                                if (shiftIOVs(count - idx, curr)) {
                                    slotsleft += 1;
                                    outbase -= sizeof(DAME_VECTOR);
                                }
                            }
                        }
                        sizeleft -= copied;
                        slotsleft -= numblocks;
                        stack[sp].countLeft -= (count - idx);
                        if (stack[sp].countLeft == 0) {
                            inbase = stack[sp].base + dl[sp].extent;
                            goto dlpop;
                        }
                    }
                    idx = count;
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[idx];
                    inbase = stack[sp].base + dl[sp].s.i_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_CONTIGFINAL:{
                    DBG_PRINTF(stderr, "Invalid FINAL state in pop(%d)!\n", dl[sp].kind);
                    MPIR_Assert(0);
                    return -1;
                }
            case DL_RETURNTO:{
                    DBG_PRINTF(stderr, "POP: RETURNTO[%d](%returnto=d)\n", sp, dl[sp].returnto);
                    int returnto = dl[sp].returnto;
                    dl = stack[sp].prevdl;
                    // We add 1 to sp because when we pop it, it will be decremented
                    sp = returnto + 1;
                    goto dlpop;
                    break;
                }
                /* By using an EXIT command on pop, we never need to check the stack
                 * pointer */
            case DL_EXIT:{
                    DBG_PRINTF(stderr, "POP:EXIT[%d]\n", sp);
                    state->sp = 0;
                    state->partial = 0;
                    *copy_p = outSize - sizeleft;
                    *iovs_p = outSlots - slotsleft;
                    return 0;
                }
            case DL_BOTTOM:{
                    DBG_PRINTF(stderr, "POP:BOTTOM[%d] (Error!)\n", sp);
                    MPIR_Assert(0);
                    return -1;
                }
        }
    } while (1);
    /* No way to reach here? */
    return -1;
}

int DL_pack(const void *inPtr,
            void *outPtr,
            DAME_Offset outSize,
            const DAME_Dame * dl, int isbuiltin, DAME_State * state, DAME_Offset * copy_p)
{
    MPI_Offset inbase = (MPI_Offset) inPtr;
    MPI_Offset outbase = (MPI_Offset) outPtr;
    long sp = 0;
    MPI_Aint sizeleft = outSize;

    /* Local variables */
    DAME_Stackelm *stack = state->stack;
    /* If it's a resume, this will get set correctly later */
    DAME_Offset partial = 0;

    *copy_p = 0;

    /* This is a resume */
    if (state->sp != 0) {
        sp = state->sp;
        inbase = state->stack[sp].base;
        partial = state->partial;

        /* NOTE: For builtin types, we create a dataloop as needed on the stack
         * so we don't save it in the state because when we resume, the same
         * type will be recreated in segment_packunpack. */
        if (!isbuiltin)
            dl = state->dl;

        /* By construction, the last state pushed into the stack will always
         * be a CONTIGFINAL. Need to decrement the stack pointer to ensure
         * that it is set to the right place by the immediate next call to
         * dlpush */
        sp--;
        goto dlpush;
    }

    do {
      dlpush:
        sp++;
        switch (dl[sp].kind) {
            case DL_CONTIG:{
                    DBG_PRINTF(stderr, "contig%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    goto dlpush;
                }
            case DL_VECTOR:{
                    DBG_PRINTF(stderr, "vec%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp + 1].base = inbase;
                    stack[sp + 1].countLeft = dl[sp].s.v_t.blklen;
                    goto dlpush;
                    break;
                }
            case DL_VECTOR1:{
                    DBG_PRINTF(stderr, "vec%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX:{
                    DBG_PRINTF(stderr, "blkidx%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp + 1].base = inbase;
                    stack[sp + 1].countLeft = dl[sp].s.bi_t.blklen;
                    inbase += dl[sp].s.bi_t.offsets[0];
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX1:{
                    DBG_PRINTF(stderr, "blkidx%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    inbase += dl[sp].s.bi_t.offsets[0];
                    goto dlpush;
                    break;
                }
            case DL_INDEX:{
                    DBG_PRINTF(stderr, "idx%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp + 1].base = inbase;
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[0];
                    inbase += dl[sp].s.i_t.offsets[0];
                    goto dlpush;
                    break;
                }
            case DL_STRUCT:{
                    DBG_PRINTF(stderr, "struct%d_push(%p, %ld)\n", sp, inbase, dl[sp].count);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp].prevdl = dl;
                    inbase += dl[sp].s.s_t.offsets[0];
                    dl = dl[sp].s.s_t.dls[0];
                    goto dlpush;
                    break;
                }
            case DL_CONTIGCHILD:{
                    DBG_PRINTF(stderr, "cc%d_push\n", sp);
                    /* countLeft set by parent */
                    /* inbase set by parent */
                    goto dlpush;
                }
            case DL_VECTORFINAL:{
                    DBG_PRINTF(stderr, "vf%d_push(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;

                    MPI_Aint blklen = dl[sp].s.v_t.blklen;
                    MPI_Aint oldsize = dl[sp].s.v_t.oldsize;
                    MPI_Aint stride = dl[sp].s.v_t.stride;

                    if (sizeleft >= dl[sp].size) {
                        if (outbase)
                            DL_pack_vector((const void *) inbase,
                                           (void **) &outbase, &dl[sp], dl[sp].count);
                        stack[sp].countLeft = 0;
                        sizeleft -= dl[sp].size;
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                        break;
                    }
                    /* destination too small.  Determine how far we can go */
                    DAME_Offset count = sizeleft / (blklen * oldsize);
                    if (count > 0) {
                        if (outbase)
                            DL_pack_vector((const void *) inbase,
                                           (void **) &outbase, &dl[sp], count);
                        sizeleft -= count * blklen * oldsize;
                        stack[sp].countLeft -= count;
                    }
                    DAME_Offset idx = count;
                    stack[sp + 1].countLeft = blklen;
                    inbase = stack[sp].base + idx * stride;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEXFINAL:{
                    DBG_PRINTF(stderr, "bf%d_push(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;

                    MPI_Aint blklen = dl[sp].s.bi_t.blklen;
                    MPI_Aint oldsize = dl[sp].s.bi_t.oldsize;

                    if (sizeleft >= dl[sp].size) {
                        if (outbase)
                            DL_pack_blockindexed((const void *) stack[sp].base,
                                                 (void **) &outbase, &dl[sp], 0, dl[sp].count);
                        stack[sp].countLeft = 0;
                        sizeleft -= dl[sp].size;
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                        break;
                    }
                    /* destination too small.  Determine how far we can go */
                    DAME_Offset count = sizeleft / (blklen * oldsize);
                    if (count > 0) {
                        if (outbase)
                            DL_pack_blockindexed((const void *) stack[sp].base,
                                                 (void **) &outbase, &dl[sp], 0, count);
                        sizeleft -= count * blklen * oldsize;
                        stack[sp].countLeft -= count;
                    }
                    DAME_Offset idx = count;
                    stack[sp + 1].countLeft = blklen;
                    inbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_INDEXFINAL:{
                    DBG_PRINTF(stderr, "if%d_push(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);
                    stack[sp].base = inbase;
                    stack[sp].countLeft = dl[sp].count;

                    MPI_Aint oldsize = dl[sp].s.i_t.oldsize;

                    if (sizeleft >= dl[sp].size) {
                        if (outbase)
                            DL_pack_indexed((const void *) stack[sp].base,
                                            (void **) &outbase, &dl[sp], 0, dl[sp].count);
                        stack[sp].countLeft = 0;
                        sizeleft -= dl[sp].size;
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }

                    /* destination too small.  Determine how far we can go */
                    DAME_Offset count, copied;
                    for (count = 0, copied = 0;
                         count < dl[sp].count
                         && copied + dl[sp].s.i_t.blklens[count] * oldsize <= sizeleft; count++)
                        copied += dl[sp].s.i_t.blklens[count] * oldsize;

                    if (count > 0) {
                        if (outbase)
                            DL_pack_indexed((const void *) stack[sp].base,
                                            (void **) &outbase, &dl[sp], 0, count);
                        sizeleft -= copied;
                        stack[sp].countLeft -= count;
                    }
                    DAME_Offset idx = count;
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[idx];
                    inbase = stack[sp].base + dl[sp].s.i_t.offsets[idx];
                    goto dlpush;
                    break;
                }

            case DL_CONTIGFINAL:{
                    DBG_PRINTF(stderr, "cf%d_push(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);

                    stack[sp].base = inbase;
                    DAME_Offset basesize = dl[sp].s.c_t.basesize;

                    /* Quick check to see if there's enough buffer space to copy.
                     * sizeleft can be zero if the previous *FINAL type has blksize
                     * which is a multiple of in/outsize */
                    if (sizeleft == 0) {
                        dlStackCopy(sp, dl, 0, state);
                        *copy_p = outSize - sizeleft;
                        return 1;
                    }

                    /* This case occurs when the contigfinal is directly below a
                     * contig which is part of a struct. There ought to be a better
                     * way of doing this, and eventually,
                     * I'll try and get rid of this check */
                    if (stack[sp].countLeft == 0 && partial == 0)
                        stack[sp].countLeft = dl[sp].count;

                    /* partial contains the remainder of an element that has been
                     * partially packed
                     * countLeft contains the number of complete blocks left */
                    DAME_Offset tocopy = stack[sp].countLeft * basesize + partial;
                    if (tocopy > sizeleft) {
                        if (outbase)
                            DL_pack_contig((const void *) inbase,
                                           (void *) outbase, &dl[sp], sizeleft);
                        DAME_Offset newLeft = stack[sp].countLeft;
                        if (partial >= sizeleft) {
                            partial -= sizeleft;
                        } else {
                            DAME_Offset avail = sizeleft - partial;
                            newLeft -= (avail / basesize);
                            partial = 0;
                            if (avail % basesize) {
                                newLeft -= 1;
                                partial = basesize - avail % basesize;
                            }
                        }
                        stack[sp].base = inbase + sizeleft;
                        stack[sp].countLeft = newLeft;
                        dlStackCopy(sp, dl, partial, state);
                        sizeleft -= sizeleft;
                        *copy_p = outSize - sizeleft;
                        return 1;
                    } else {
                        if (outbase) {
                            DL_pack_contig((const void *) inbase,
                                           (void *) outbase, &dl[sp], tocopy);
                            outbase += tocopy;
                        }
                        partial = 0;
                        inbase += tocopy;
                        sizeleft -= tocopy;
                        stack[sp].countLeft = 0;
                        goto dlpop;
                    }
                    break;
                }
            default:{
                    /* This can never happen */
                    DBG_PRINTF(stderr, "Unknown dataloop kind(%d)!\n", dl[sp].kind);
                    MPIR_Assert(0);
                    return -1;
                }
        }
      dlpop:
        sp--;

        switch (dl[sp].kind) {
            case DL_CONTIG:{
                    DBG_PRINTF(stderr, "contig%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + idx * dl[sp].s.c_t.baseextent;
                    goto dlpush;
                    break;
                }
            case DL_VECTOR:{
                    DBG_PRINTF(stderr, "vec%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    /* For vectors, we don't need to compute offsets from the base.
                     * The stride is the distance from each element, so we can use
                     * the stack contents to update */
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + idx * dl[sp].s.v_t.stride;
                    stack[sp + 1].countLeft = dl[sp].s.v_t.blklen;
                    goto dlpush;
                    break;
                }
            case DL_VECTOR1:{
                    DBG_PRINTF(stderr, "vec%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + idx * dl[sp].s.v_t.stride;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX:{
                    DBG_PRINTF(stderr, "blkidx%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    stack[sp + 1].countLeft = dl[sp].s.bi_t.blklen;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX1:{
                    DBG_PRINTF(stderr, "blkidx%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_INDEX:{
                    DBG_PRINTF(stderr, "idx%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + dl[sp].s.i_t.offsets[idx];
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[idx];
                    goto dlpush;
                    break;
                }
            case DL_STRUCT:{
                    DBG_PRINTF(stderr, "struct%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0) {
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    inbase = stack[sp].base + dl[sp].s.s_t.offsets[idx];
                    dl = dl[sp].s.s_t.dls[idx];
                    goto dlpush;
                    break;
                }
            case DL_CONTIGCHILD:{
                    DBG_PRINTF(stderr, "cc%d_pop(%p, %ld)\n", sp, inbase, dl[sp].count);
                    if (--stack[sp].countLeft == 0) {
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }
                    goto dlpush;
                    break;
                }
            case DL_VECTORFINAL:{
                    DBG_PRINTF(stderr, "vf%d_pop(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);

                    MPI_Aint blklen = dl[sp].s.v_t.blklen;
                    MPI_Aint oldsize = dl[sp].s.v_t.oldsize;
                    MPI_Aint stride = dl[sp].s.v_t.stride;

                    if (--stack[sp].countLeft == 0) {
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    DAME_Offset count = sizeleft / (blklen * oldsize);
                    if (count > 0) {
                        if (count > stack[sp].countLeft)
                            count = stack[sp].countLeft;
                        if (outbase)
                            DL_pack_vector((const void *) (stack[sp].base + idx * stride),
                                           (void **) &outbase, &dl[sp], count);
                        sizeleft -= count * blklen * oldsize;
                        stack[sp].countLeft -= count;
                        if (stack[sp].countLeft == 0) {
                            inbase = stack[sp].base + dl[sp].extent;
                            goto dlpop;
                        }
                    }
                    idx += count;
                    stack[sp + 1].countLeft = blklen;
                    inbase = stack[sp].base + idx * stride;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEXFINAL:{
                    DBG_PRINTF(stderr, "bf%d_pop(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);
                    if (--stack[sp].countLeft == 0) {
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }

                    MPI_Aint blklen = dl[sp].s.bi_t.blklen;
                    MPI_Aint oldsize = dl[sp].s.bi_t.oldsize;

                    /* If there's still blocks left, try to do several of them */
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    DAME_Offset count = sizeleft / (blklen * oldsize);
                    if (count > 0) {
                        if (count > stack[sp].countLeft)
                            count = stack[sp].countLeft;
                        if (outbase)
                            DL_pack_blockindexed((const void *) stack[sp].base,
                                                 (void **) &outbase, &dl[sp], idx, idx + count);
                        sizeleft -= count * blklen * oldsize;
                        stack[sp].countLeft -= count;
                        if (stack[sp].countLeft == 0) {
                            inbase = stack[sp].base + dl[sp].extent;
                            goto dlpop;
                        }
                    }
                    idx += count;
                    stack[sp + 1].countLeft = blklen;
                    inbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_INDEXFINAL:{
                    DBG_PRINTF(stderr, "if%d_pop(%p, %ld, %ld, %ld)\n",
                               sp, inbase, dl[sp].count, dl[sp].size, sizeleft);
                    if (--stack[sp].countLeft == 0) {
                        inbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }

                    MPI_Aint oldsize = dl[sp].s.i_t.oldsize;

                    /* If there are still blocks left, try to do several at once */
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    DAME_Offset count, copied;
                    for (count = idx, copied = 0;
                         count < dl[sp].count
                         && copied + dl[sp].s.i_t.blklens[count] * oldsize <= sizeleft; count++)
                        copied += dl[sp].s.i_t.blklens[count] * oldsize;
                    if (count > idx) {
                        if (outbase)
                            DL_pack_indexed((const void *) stack[sp].base,
                                            (void **) &outbase, &dl[sp], idx, count);
                        sizeleft -= copied;
                        stack[sp].countLeft -= (count - idx);
                        if (stack[sp].countLeft == 0) {
                            inbase = stack[sp].base + dl[sp].extent;
                            goto dlpop;
                        }
                    }
                    idx = count;
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[idx];
                    inbase = stack[sp].base + dl[sp].s.i_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_CONTIGFINAL:{
                    DBG_PRINTF(stderr, "Invalid FINAL state in pop(%d)!\n", dl[sp].kind);
                    MPIR_Assert(0);
                    return -1;
                }
            case DL_RETURNTO:{
                    DBG_PRINTF(stderr, "POP: RETURNTO[%d](%returnto=d)\n", sp, dl[sp].returnto);
                    int returnto = dl[sp].returnto;
                    dl = stack[sp].prevdl;
                    // We add 1 to sp because when we pop it, it will be decremented
                    sp = returnto + 1;
                    goto dlpop;
                    break;
                }
                /* By using an EXIT command on pop, we never need to check the stack
                 * pointer */
            case DL_EXIT:{
                    DBG_PRINTF(stderr, "POP:EXIT[%d]\n", sp);
                    state->sp = 0;
                    state->partial = 0;
                    *copy_p = outSize - sizeleft;
                    return 0;
                }
            case DL_BOTTOM:{
                    DBG_PRINTF(stderr, "POP:BOTTOM[%d] (Error!)\n", sp);
                    MPIR_Assert(0);
                    return -1;
                }
        }
    } while (1);
    /* No way to reach here? */
    return -1;
}

/********** UNPACK **************/

/*
  Return values:
  -1: Error
  0: Done
  1: Needs to continue; state in state_p
*/
int DL_unpack(const void *inPtr,
              void *outPtr,
              DAME_Offset outSize,
              const DAME_Dame * dl, int isbuiltin, DAME_State * state, DAME_Offset * copy_p)
{
    MPI_Offset inbase = (MPI_Offset) inPtr;
    MPI_Offset outbase = (MPI_Offset) outPtr;
    long sp = 0;
    MPI_Aint sizeleft = outSize;

    /* Local variables */
    DAME_Stackelm *stack = state->stack;
    /* If it's a resume, this will get set correctly later */
    DAME_Offset partial = 0;

    *copy_p = 0;

    /* This is a resume */
    if (state->sp != 0) {
        sp = state->sp;
        outbase = state->stack[sp].base;
        partial = state->partial;

        /* NOTE: For builtin types, we create a dataloop as needed on the stack
         * so we don't save it in the state because when we resume, the same
         * type will be recreated in segment_packunpack. */
        if (!isbuiltin)
            dl = state->dl;

        /* By construction, the last state pushed into the stack will always
         * be a CONTIGFINAL. Need to decrement the stack pointer to ensure
         * that it is set to the right place by the immediate next call to
         * dlpush */
        sp--;
        goto dlpush;
    }

    do {
      dlpush:
        sp++;
        switch (dl[sp].kind) {
            case DL_CONTIG:{
                    DBG_PRINTF(stderr, "PUSH: CONTIG[%d]\n", sp);
                    stack[sp].base = outbase;
                    stack[sp].countLeft = dl[sp].count;
                    goto dlpush;
                }
            case DL_VECTOR:{
                    DBG_PRINTF(stderr, "PUSH: VECTOR[%d]\n", sp);
                    stack[sp].base = outbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp + 1].base = outbase;
                    stack[sp + 1].countLeft = dl[sp].s.v_t.blklen;
                    goto dlpush;
                    break;
                }
            case DL_VECTOR1:{
                    DBG_PRINTF(stderr, "PUSH: VECTOR1[%d]\n", sp);
                    stack[sp].base = outbase;
                    stack[sp].countLeft = dl[sp].count;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX:{
                    DBG_PRINTF(stderr, "PUSH: BLOCKINDEX[%d]\n", sp);
                    stack[sp].base = outbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp + 1].base = outbase;
                    stack[sp + 1].countLeft = dl[sp].s.bi_t.blklen;
                    outbase += dl[sp].s.bi_t.offsets[0];
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX1:{
                    DBG_PRINTF(stderr, "PUSH: BLOCKINDEX1[%d]\n", sp);
                    stack[sp].base = outbase;
                    stack[sp].countLeft = dl[sp].count;
                    outbase += dl[sp].s.bi_t.offsets[0];
                    goto dlpush;
                    break;
                }
            case DL_INDEX:{
                    DBG_PRINTF(stderr, "PUSH: INDEX[%d]\n", sp);
                    stack[sp].base = outbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp + 1].base = outbase;
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[0];
                    outbase += dl[sp].s.i_t.offsets[0];
                    goto dlpush;
                    break;
                }
            case DL_STRUCT:{
                    DBG_PRINTF(stderr, "PUSH: STRUCT[%d](absbase=%p)\n", sp, inbase);
                    stack[sp].base = outbase;
                    stack[sp].countLeft = dl[sp].count;
                    stack[sp].prevdl = dl;
                    outbase += dl[sp].s.s_t.offsets[0];
                    dl = dl[sp].s.s_t.dls[0];
                    goto dlpush;
                    break;
                }
            case DL_CONTIGCHILD:{
                    DBG_PRINTF(stderr, "PUSH: CONTIGCHILD[%d]\n", sp);
                    /* countLeft set by parent */
                    /* base set by parent */
                    goto dlpush;
                }
            case DL_VECTORFINAL:{
                    DBG_PRINTF(stderr, "PUSH: VECFINAL[%d] (sizeleft=%ld)\n", sp, sizeleft);
                    stack[sp].base = outbase;
                    stack[sp].countLeft = dl[sp].count;

                    if (sizeleft >= dl[sp].size) {
                        if (inbase)
                            DL_unpack_vector((const void **) &inbase,
                                             (void *) outbase, &dl[sp], dl[sp].count);
                        stack[sp].countLeft = 0;
                        sizeleft -= dl[sp].size;
                        outbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                        break;  /* to dlpop */
                    }
                    /* destination too small.  Determine how far we can go */
                    DAME_Offset count = sizeleft / (dl[sp].s.v_t.blklen * dl[sp].s.v_t.oldsize);
                    if (count > 0) {
                        if (inbase)
                            DL_unpack_vector((const void **) &inbase,
                                             (void *) outbase, &dl[sp], count);
                        sizeleft -= count * dl[sp].s.v_t.blklen * dl[sp].s.v_t.oldsize;
                        stack[sp].countLeft -= count;
                    }
                    DAME_Offset idx = count;
                    stack[sp + 1].countLeft = dl[sp].s.v_t.blklen;
                    outbase = stack[sp].base + idx * dl[sp].s.v_t.stride;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEXFINAL:{
                    DBG_PRINTF(stderr, "(u)PUSH: BLOCKINDEXFINAL[%d](sizeleft=%ld)\n",
                               sp, sizeleft);
                    stack[sp].base = outbase;
                    stack[sp].countLeft = dl[sp].count;

                    if (sizeleft >= dl[sp].size) {
                        if (inbase)
                            DL_unpack_blockindexed((const void **) &inbase,
                                                   (void *) stack[sp].base,
                                                   &dl[sp], 0, dl[sp].count);
                        stack[sp].countLeft = 0;
                        sizeleft -= dl[sp].size;
                        outbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                        break;  /* to dlpop */
                    }
                    /* destination too small.  Determine how far we can go */
                    DAME_Offset count = sizeleft / (dl[sp].s.bi_t.blklen * dl[sp].s.bi_t.oldsize);
                    if (count > 0) {
                        if (inbase)
                            DL_unpack_blockindexed((const void **) &inbase,
                                                   (void *) stack[sp].base, &dl[sp], 0, count);
                        sizeleft -= count * dl[sp].s.bi_t.blklen * dl[sp].s.bi_t.oldsize;
                        stack[sp].countLeft -= count;
                    }
                    DAME_Offset idx = count;
                    stack[sp + 1].countLeft = dl[sp].s.bi_t.blklen;
                    outbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_INDEXFINAL:{
                    DBG_PRINTF(stderr, "(u)PUSH: INDEXFINAL[%d](sizeleft=%ld)\n", sp, sizeleft);
                    stack[sp].base = outbase;
                    stack[sp].countLeft = dl[sp].count;

                    if (sizeleft >= dl[sp].size) {
                        if (inbase)
                            DL_unpack_indexed((const void **) &inbase,
                                              (void *) stack[sp].base, &dl[sp], 0, dl[sp].count);
                        stack[sp].countLeft = 0;
                        sizeleft -= dl[sp].size;
                        outbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                        break;  /* to dlpop */
                    }

                    /* destination too small.  Determine how far we can go */
                    DAME_Offset count, copied;
                    for (count = 0, copied = 0;
                         copied + dl[sp].s.i_t.blklens[count] * dl[sp].s.i_t.oldsize <= sizeleft;
                         count++)
                        copied += dl[sp].s.i_t.blklens[count] * dl[sp].s.i_t.oldsize;
                    if (count > 0) {
                        DL_unpack_indexed((const void **) &inbase,
                                          (void *) stack[sp].base, &dl[sp], 0, count);
                        sizeleft -= copied;
                        stack[sp].countLeft -= count;
                    }
                    DAME_Offset idx = count;
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[idx];
                    outbase = stack[sp].base + dl[sp].s.i_t.offsets[idx];
                    goto dlpush;
                    break;
                }

            case DL_CONTIGFINAL:{
                    DBG_PRINTF(stderr, "(u)PUSH: CONTIGFINAL[%d](sizeleft=%ld)\n", sp, sizeleft);

                    stack[sp].base = outbase;

                    /* Quick check to see if there's enough buffer space to copy.
                     * sizeleft can be zero if the previous *FINAL type has blksize
                     * which is a multiple of in/outsize */
                    if (sizeleft == 0) {
                        dlStackCopy(sp, dl, 0, state);
                        *copy_p = outSize - sizeleft;
                        return 1;
                    }

                    /* This case occurs when the contigfinal is directly below a contig
                     * which is part of a struct. There ought to be a better way of
                     * doing this, and eventually, I'll try and get rid of this check */
                    if (stack[sp].countLeft == 0 && partial == 0)
                        stack[sp].countLeft = dl[sp].count;

                    DAME_Offset basesize = dl[sp].s.c_t.basesize;
                    DAME_Offset tocopy = stack[sp].countLeft * basesize + partial;
                    if (tocopy > sizeleft) {
                        if (inbase)
                            DL_pack_contig((const void *) inbase,
                                           (void *) outbase, &dl[sp], sizeleft);
                        DAME_Offset newLeft = stack[sp].countLeft;
                        if (partial >= sizeleft) {
                            partial -= sizeleft;
                        } else {
                            DAME_Offset avail = sizeleft - partial;
                            newLeft -= (avail / basesize);
                            partial = 0;
                            if (avail % basesize) {
                                newLeft -= 1;
                                partial = basesize - avail % basesize;
                            }
                        }
                        stack[sp].base = outbase + sizeleft;
                        stack[sp].countLeft = newLeft;
                        dlStackCopy(sp, dl, partial, state);
                        sizeleft -= sizeleft;
                        *copy_p = outSize - sizeleft;
                        return 1;
                    } else {
                        if (inbase) {
                            DL_pack_contig((const void *) inbase,
                                           (void *) outbase, &dl[sp], tocopy);
                            inbase += tocopy;
                        }
                        partial = 0;
                        outbase += tocopy;
                        sizeleft -= tocopy;
                        stack[sp].countLeft = 0;
                        goto dlpop;
                    }
                    break;
                }

            default:{
                    /* This can never happen */
                    DBG_PRINTF(stderr, "The impossible happened(%d)!\n", dl[sp].kind);
                    MPIR_Assert(0);
                    return -1;
                }
        }

      dlpop:
        sp--;

        /* Question: can we ever pop INTO an xxFINAL? Or will those all
         * be handled in the continue-from-state? */

        switch (dl[sp].kind) {
            case DL_CONTIG:{
                    DBG_PRINTF(stderr, "POP: CONTIG[%d]\n", sp);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    outbase = stack[sp].base + idx * dl[sp].s.c_t.baseextent;
                    goto dlpush;
                    break;
                }
            case DL_VECTOR:{
                    DBG_PRINTF(stderr, "POP: VEC[%d]\n", sp);
                    /* For vectors, we don't need to compute offsets from the base.
                     * The stride is the distance from each element, so we can use
                     * the stack contents to update */
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    outbase = stack[sp].base + idx * dl[sp].s.v_t.stride;
                    stack[sp + 1].countLeft = dl[sp].s.v_t.blklen;
                    goto dlpush;
                    break;
                }
            case DL_VECTOR1:{
                    DBG_PRINTF(stderr, "POP: VECTOR1[%d]\n", sp);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    outbase = stack[sp].base + idx * dl[sp].s.v_t.stride;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX:{
                    DBG_PRINTF(stderr, "POP: BLOCKINDEX[%d] (countLeft=%ld)\n",
                               sp, stack[sp].countLeft);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    outbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    stack[sp + 1].countLeft = dl[sp].s.bi_t.blklen;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEX1:{
                    DBG_PRINTF(stderr, "POP: BLOCKINDEX1[%d]\n", sp);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    outbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_INDEX:{
                    DBG_PRINTF(stderr, "POP: INDEX[%d] (countLeft=%ld)\n", sp, stack[sp].countLeft);
                    if (--stack[sp].countLeft == 0)
                        goto dlpop;
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[idx];
                    outbase = stack[sp].base + dl[sp].s.i_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_STRUCT:{
                    DBG_PRINTF(stderr, "POP: STRUCT[%d] (countLeft=%ld)\n",
                               sp, stack[sp].countLeft - 1);
                    if (--stack[sp].countLeft == 0) {
                        outbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    outbase = stack[sp].base + dl[sp].s.s_t.offsets[idx];
                    dl = dl[sp].s.s_t.dls[idx];
                    goto dlpush;
                    break;
                }
            case DL_CONTIGCHILD:{
                    DBG_PRINTF(stderr, "POP: CONTIGCHILD[%d](countLeft=%d)\n", sp,
                               stack[sp].countLeft - 1);
                    if (--stack[sp].countLeft == 0) {
                        outbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }
                    goto dlpush;
                    break;
                }
            case DL_VECTORFINAL:{
                    DBG_PRINTF(stderr,
                               "POP: VECFINAL[%d](countLeft=%ld)\n", sp, stack[sp].countLeft);
                    if (--stack[sp].countLeft == 0) {
                        outbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    DAME_Offset count = sizeleft / (dl[sp].s.v_t.blklen * dl[sp].s.v_t.oldsize);
                    if (count > 0) {
                        if (count > stack[sp].countLeft)
                            count = dl[sp].count - idx;
                        if (inbase)
                            DL_unpack_vector((const void **) &inbase,
                                             (void *) (stack[sp].base + idx * dl[sp].s.v_t.stride),
                                             &dl[sp], count);
                        sizeleft -= count * dl[sp].s.v_t.blklen * dl[sp].s.v_t.oldsize;
                        stack[sp].countLeft -= count;
                        if (stack[sp].countLeft == 0) {
                            outbase = stack[sp].base + dl[sp].extent;
                            goto dlpop;
                        }
                    }
                    idx += count;
                    stack[sp + 1].countLeft = dl[sp].s.v_t.blklen;
                    outbase = stack[sp].base + idx * dl[sp].s.v_t.stride;
                    goto dlpush;
                    break;
                }
            case DL_BLOCKINDEXFINAL:{
                    DBG_PRINTF(stderr, "POP: BLOCKINDEXFINAL[%d]\n", sp);
                    if (--stack[sp].countLeft == 0) {
                        outbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }

                    /* If there's still blocks left, try to do several of them */
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    DAME_Offset count = sizeleft / (dl[sp].s.bi_t.blklen * dl[sp].s.bi_t.oldsize);
                    if (count > 0) {
                        if (count > stack[sp].countLeft)
                            count = dl[sp].count - idx;
                        if (inbase)
                            DL_unpack_blockindexed((const void **) &inbase,
                                                   (void *) stack[sp].base,
                                                   &dl[sp], idx, idx + count);
                        sizeleft -= count * dl[sp].s.bi_t.blklen * dl[sp].s.bi_t.oldsize;
                        stack[sp].countLeft -= count;
                        if (stack[sp].countLeft == 0) {
                            outbase = stack[sp].base + dl[sp].extent;
                            goto dlpop;
                        }
                    }
                    idx += count;
                    stack[sp + 1].countLeft = dl[sp].s.bi_t.blklen;
                    outbase = stack[sp].base + dl[sp].s.bi_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_INDEXFINAL:{
                    DBG_PRINTF(stderr, "POP: INDEXFINAL[%d]\n", sp);
                    if (--stack[sp].countLeft == 0) {
                        outbase = stack[sp].base + dl[sp].extent;
                        goto dlpop;
                    }

                    /* If there are still blocks left, try to do several at once */
                    DAME_Offset idx = dl[sp].count - stack[sp].countLeft;
                    DAME_Offset count = idx, copied = 0;
                    for (count = idx, copied = 0;
                         count < dl[sp].count
                         && copied + dl[sp].s.i_t.blklens[count] * dl[sp].s.i_t.oldsize <= sizeleft;
                         count++)
                        copied += dl[sp].s.i_t.blklens[count] * dl[sp].s.i_t.oldsize;
                    if (count > idx) {
                        if (inbase)
                            DL_unpack_indexed((const void **) &inbase,
                                              (void *) stack[sp].base, &dl[sp], idx, count);
                        sizeleft -= copied;
                        stack[sp].countLeft -= (count - idx);
                        if (stack[sp].countLeft == 0) {
                            outbase = stack[sp].base + dl[sp].extent;
                            goto dlpop;
                        }
                    }
                    idx = count;
                    stack[sp + 1].countLeft = dl[sp].s.i_t.blklens[idx];
                    outbase = stack[sp].base + dl[sp].s.i_t.offsets[idx];
                    goto dlpush;
                    break;
                }
            case DL_CONTIGFINAL:{
                    DBG_PRINTF(stderr, "Invalid FINAL state in pop(%d)!\n", dl[sp].kind);
                    MPIR_Assert(0);
                    return -1;
                }
            case DL_RETURNTO:{
                    DBG_PRINTF(stderr, "POP: RETURNTO[%d](%returnto=%d)\n", sp, dl[sp].returnto);
                    int returnto = dl[sp].returnto;
                    dl = stack[sp].prevdl;
                    // We add 1 to sp because at dlpop it will be decremented
                    sp = returnto + 1;
                    goto dlpop;
                    break;
                }
                /* By using an EXIT command on pop, we never need to check the stack
                 * pointer */
            case DL_EXIT:{
                    DBG_PRINTF(stderr, "POP:EXIT[%d]\n", sp);
                    state->sp = 0;
                    state->partial = 0;
                    *copy_p = outSize - sizeleft;
                    return 0;
                }
            case DL_BOTTOM:{
                    DBG_PRINTF(stderr, "POP:BOTTOM[%d] (Error!)\n", sp);
                    MPIR_Assert(0);
                    return -1;
                }
        }
    } while (1);
    /* No way to reach here? */
    return -1;
}

/* Segment_init
 *
 * buf    - datatype buffer location
 * count  - number of instances of the datatype in the buffer
 * handle - handle for datatype (could be derived or not)
 * segp   - pointer to previously allocated segment structure
 * flag   - flag is unused
 *
 * Notes:
 * - Assumes that the segment has been allocated.
 * - Older MPICH code may pass "0" to indicate HETEROGENEOUS or "1" to
 *   indicate HETEROGENEOUS.
 *
 */
int MPIR_Segment_init(const DAME_Buffer buf,
                      DAME_Count count, DAME_Handle handle, struct DAME_Segment *segp, int flag)
{
    int i;
    segp->ptr = (DAME_Buffer) buf;
    segp->handle = handle;
    segp->count = count;

    for (i = 0; i < 3; i++) {
        segp->in_offset[i] = 0;
        segp->out_offset[i] = 0;
        segp->curcount[i] = 0;

        segp->state[i].sp = 0;
        segp->state[i].dl = 0;
        segp->state[i].partial = 0;
        memset(segp->state[i].stack, 0, DAME_MAX_STACK_SIZE * sizeof(DAME_Stackelm));
    }

    return 0;
}

/* Segment_alloc
 *
 */
struct DAME_Segment *MPIR_Segment_alloc(void)
{
    int i;
    DAME_Segment *segp = (struct DAME_Segment *) DAME_Malloc(sizeof(struct DAME_Segment));
    MPIR_Segment_init(NULL, 0, MPI_DATATYPE_NULL, segp, 0);
    return segp;
}

/* Segment_free
 *
 * Input Parameters:
 * segp - pointer to segment
 */
void MPIR_Segment_free(struct DAME_Segment *segp)
{
    DAME_Free(segp);
    return;
}
