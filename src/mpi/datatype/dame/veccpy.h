/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef VECCPY_H
#define VECCPY_H

#ifdef HAVE_ANY_INT64_T_ALIGNEMENT
#define MPIR_ALIGN8_TEST(p1,p2) 0
#else
#define MPIR_ALIGN8_TEST(p1,p2) (((DAME_VOID_PTR_CAST_TO_OFFSET p1 | DAME_VOID_PTR_CAST_TO_OFFSET p2) & 0x7) == 0)
#endif

#ifdef HAVE_ANY_INT32_T_ALIGNEMENT
#define MPIR_ALIGN4_TEST(p1,p2) 0
#else
#define MPIR_ALIGN4_TEST(p1,p2) (((DAME_VOID_PTR_CAST_TO_OFFSET p1 | DAME_VOID_PTR_CAST_TO_OFFSET p2) & 0x3) == 0)
#endif


#define MPIR_SRC_DEF(type, src, aligned)  type *l_src = (type*)src;
#define MPIR_DST_DEF(type, dest, aligned)  type *l_dest = (type*)dest;
#define MPIR_TMPSRC_DEF(type, l_src, aligned) type *tmp_src = l_src;
#define MPIR_TMPDST_DEF(type, l_dest, aligned) type *tmp_dest = l_dest;


#define MPIDI_COPY_FROM_VEC_ALIGNED(src,dest,stride,type,nelms,count,aligned) \
  {                                                                     \
    MPIR_SRC_DEF(type, src, aligned);                                   \
    MPIR_DST_DEF(type, dest, aligned);                                  \
    MPIR_TMPSRC_DEF(type, l_src, aligned);                              \
                                                                        \
    register int k;                                                     \
    register unsigned long _i, j;                                       \
    unsigned long total_count = count * nelms;                          \
    const DAME_Offset l_stride = stride;                               \
                                                                        \
    DAME_Assert(stride <= INT_MAX);                                    \
    DAME_Assert(total_count <= INT_MAX);                               \
    DAME_Assert(nelms <= INT_MAX);                                     \
    if (nelms == 1) {                                                   \
      for (_i = (int)total_count; _i; _i--) {                           \
        *l_dest++ = *l_src;                                             \
        l_src += l_stride;                                              \
      }                                                                 \
    }                                                                   \
    else if (nelms == 2) {                                              \
      for (_i = (int)total_count; _i; _i -= 2) {                        \
        *l_dest++ = l_src[0];                                           \
        *l_dest++ = l_src[1];                                           \
        l_src += l_stride;                                              \
      }                                                                 \
    }                                                                   \
    else if (nelms == 3) {                                              \
      for (_i = (int)total_count; _i; _i -= 3) {                        \
        *l_dest++ = l_src[0];                                           \
        *l_dest++ = l_src[1];                                           \
        *l_dest++ = l_src[2];                                           \
        l_src += l_stride;                                              \
      }                                                                 \
    }                                                                   \
    else if (nelms == 4) {                                              \
      for (_i = (int)total_count; _i; _i -= 4) {                        \
        *l_dest++ = l_src[0];                                           \
        *l_dest++ = l_src[1];                                           \
        *l_dest++ = l_src[2];                                           \
        *l_dest++ = l_src[3];                                           \
        l_src += l_stride;                                              \
      }                                                                 \
    }                                                                   \
    else if (nelms == 5) {                                              \
      for (_i = (int)total_count; _i; _i -= 5) {                        \
        *l_dest++ = l_src[0];                                           \
        *l_dest++ = l_src[1];                                           \
        *l_dest++ = l_src[2];                                           \
        *l_dest++ = l_src[3];                                           \
        *l_dest++ = l_src[4];                                           \
        l_src += l_stride;                                              \
      }                                                                 \
    }                                                                   \
    else if (nelms == 6) {                                              \
      for (_i = (int)total_count; _i; _i -= 6) {                        \
        *l_dest++ = l_src[0];                                           \
        *l_dest++ = l_src[1];                                           \
        *l_dest++ = l_src[2];                                           \
        *l_dest++ = l_src[3];                                           \
        *l_dest++ = l_src[4];                                           \
        *l_dest++ = l_src[5];                                           \
        l_src += l_stride;                                              \
      }                                                                 \
    }                                                                   \
    else if (nelms == 7) {                                              \
      for (_i = (int)total_count; _i; _i -= 7) {                        \
        *l_dest++ = l_src[0];                                           \
        *l_dest++ = l_src[1];                                           \
        *l_dest++ = l_src[2];                                           \
        *l_dest++ = l_src[3];                                           \
        *l_dest++ = l_src[4];                                           \
        *l_dest++ = l_src[5];                                           \
        *l_dest++ = l_src[6];                                           \
        l_src += l_stride;                                              \
      }                                                                 \
    }                                                                   \
    else if (nelms == 8) {                                              \
      for (_i = (int)total_count; _i; _i -= 8) {                        \
        *l_dest++ = l_src[0];                                           \
        *l_dest++ = l_src[1];                                           \
        *l_dest++ = l_src[2];                                           \
        *l_dest++ = l_src[3];                                           \
        *l_dest++ = l_src[4];                                           \
        *l_dest++ = l_src[5];                                           \
        *l_dest++ = l_src[6];                                           \
        *l_dest++ = l_src[7];                                           \
        l_src += l_stride;                                              \
      }                                                                 \
    }                                                                   \
    else {                                                              \
      _i = (int)total_count;                                            \
      while (_i) {                                                      \
        tmp_src = l_src;                                                \
        j = (int)nelms;                                                 \
        while (j >= 8) {                                                \
          *l_dest++ = tmp_src[0];                                       \
          *l_dest++ = tmp_src[1];                                       \
          *l_dest++ = tmp_src[2];                                       \
          *l_dest++ = tmp_src[3];                                       \
          *l_dest++ = tmp_src[4];                                       \
          *l_dest++ = tmp_src[5];                                       \
          *l_dest++ = tmp_src[6];                                       \
          *l_dest++ = tmp_src[7];                                       \
          j -= 8;                                                       \
          tmp_src += 8;                                                 \
        }                                                               \
        for (k = 0; k < j; k++) {                                       \
          *l_dest++ = *tmp_src++;                                       \
        }                                                               \
        l_src += l_stride;                                              \
        _i -= nelms;                                                    \
      }                                                                 \
    }                                                                   \
  }

#define MPIDI_COPY_TO_VEC_ALIGNED(src,dest,stride,type,nelms,count,aligned) \
  {                                                                     \
    MPIR_SRC_DEF(type, src, aligned);                                   \
    MPIR_DST_DEF(type, dest, aligned);                                  \
    MPIR_TMPDST_DEF(type, l_dest, aligned);                             \
                                                                        \
    register int k;                                                     \
    register unsigned long _i, j;                                       \
    unsigned long total_count = count * nelms;                          \
    const DAME_Offset l_stride = stride;                               \
                                                                        \
    DAME_Assert(stride <= INT_MAX);                                    \
    DAME_Assert(total_count <= INT_MAX);                               \
    DAME_Assert(nelms <= INT_MAX);                                     \
    if (nelms == 1) {                                                   \
      for (_i = (int)total_count; _i; _i--) {                           \
        *l_dest = *l_src++;                                             \
        l_dest += l_stride;                                             \
      }                                                                 \
    }                                                                   \
    else if (nelms == 2) {                                              \
      for (_i = (int)total_count; _i; _i -= 2) {                        \
        l_dest[0] = *l_src++;                                           \
        l_dest[1] = *l_src++;                                           \
        l_dest += l_stride;                                             \
      }                                                                 \
    }                                                                   \
    else if (nelms == 3) {                                              \
      for (_i = (int)total_count; _i; _i -= 3) {                        \
        l_dest[0] = *l_src++;                                           \
        l_dest[1] = *l_src++;                                           \
        l_dest[2] = *l_src++;                                           \
        l_dest += l_stride;                                             \
      }                                                                 \
    }                                                                   \
    else if (nelms == 4) {                                              \
      for (_i = (int)total_count; _i; _i -= 4) {                        \
        l_dest[0] = *l_src++;                                           \
        l_dest[1] = *l_src++;                                           \
        l_dest[2] = *l_src++;                                           \
        l_dest[3] = *l_src++;                                           \
        l_dest += l_stride;                                             \
      }                                                                 \
    }                                                                   \
    else if (nelms == 5) {                                              \
      for (_i = (int)total_count; _i; _i -= 5) {                        \
        l_dest[0] = *l_src++;                                           \
        l_dest[1] = *l_src++;                                           \
        l_dest[2] = *l_src++;                                           \
        l_dest[3] = *l_src++;                                           \
        l_dest[4] = *l_src++;                                           \
        l_dest += l_stride;                                             \
      }                                                                 \
    }                                                                   \
    else if (nelms == 6) {                                              \
      for (_i = (int)total_count; _i; _i -= 6) {                        \
        l_dest[0] = *l_src++;                                           \
        l_dest[1] = *l_src++;                                           \
        l_dest[2] = *l_src++;                                           \
        l_dest[3] = *l_src++;                                           \
        l_dest[4] = *l_src++;                                           \
        l_dest[5] = *l_src++;                                           \
        l_dest += l_stride;                                             \
      }                                                                 \
    }                                                                   \
    else if (nelms == 7) {                                              \
      for (_i = (int)total_count; _i; _i -= 7) {                        \
        l_dest[0] = *l_src++;                                           \
        l_dest[1] = *l_src++;                                           \
        l_dest[2] = *l_src++;                                           \
        l_dest[3] = *l_src++;                                           \
        l_dest[4] = *l_src++;                                           \
        l_dest[5] = *l_src++;                                           \
        l_dest[6] = *l_src++;                                           \
        l_dest += l_stride;                                             \
      }                                                                 \
    }                                                                   \
    else if (nelms == 8) {                                              \
      for (_i = (int)total_count; _i; _i -= 8) {                        \
        l_dest[0] = *l_src++;                                           \
        l_dest[1] = *l_src++;                                           \
        l_dest[2] = *l_src++;                                           \
        l_dest[3] = *l_src++;                                           \
        l_dest[4] = *l_src++;                                           \
        l_dest[5] = *l_src++;                                           \
        l_dest[6] = *l_src++;                                           \
        l_dest[7] = *l_src++;                                           \
        l_dest += l_stride;                                             \
      }                                                                 \
    }                                                                   \
    else {                                                              \
      _i = (int)total_count;                                            \
      while (_i) {                                                      \
        tmp_dest = l_dest;                                              \
        j = (int)nelms;                                                 \
        while (j >= 8) {                                                \
          tmp_dest[0] = *l_src++;                                       \
          tmp_dest[1] = *l_src++;                                       \
          tmp_dest[2] = *l_src++;                                       \
          tmp_dest[3] = *l_src++;                                       \
          tmp_dest[4] = *l_src++;                                       \
          tmp_dest[5] = *l_src++;                                       \
          tmp_dest[6] = *l_src++;                                       \
          tmp_dest[7] = *l_src++;                                       \
          j -= 8;                                                       \
          tmp_dest += 8;                                                \
        }                                                               \
        for (k = 0; k < j; k++) {                                       \
          *tmp_dest++ = *l_src++;                                       \
        }                                                               \
        l_dest += l_stride;                                             \
        _i -= nelms;                                                    \
      }                                                                 \
    }                                                                   \
  }                                                                     \

#define MPIDI_COPY_BLOCK_UNALIGNED(src,dest,type,nelms) \
  {                                                     \
    type* l_src = (type*)src;                           \
    type* l_dest = (type*)dest;                         \
    type* tmp_src = l_src;                              \
                                                        \
    register long k;                                    \
    register long j;                                    \
                                                        \
    DAME_Assert(nelms <= INT_MAX);                     \
    if (nelms == 1) {                                    \
      *l_dest++ = *l_src;                               \
    } else if (nelms == 2) {                             \
      *l_dest++ = l_src[0];                             \
      *l_dest++ = l_src[1];                             \
    } else if (nelms == 3) {                             \
      *l_dest++ = l_src[0];                             \
      *l_dest++ = l_src[1];                             \
      *l_dest++ = l_src[2];                             \
    } else if (nelms == 4) {                             \
      *l_dest++ = l_src[0];                             \
      *l_dest++ = l_src[1];                             \
      *l_dest++ = l_src[2];                             \
      *l_dest++ = l_src[3];                             \
    } else if (nelms == 5) {                             \
      *l_dest++ = l_src[0];                             \
      *l_dest++ = l_src[1];                             \
      *l_dest++ = l_src[2];                             \
      *l_dest++ = l_src[3];                             \
      *l_dest++ = l_src[4];                             \
    } else if (nelms == 6) {                             \
      *l_dest++ = l_src[0];                             \
      *l_dest++ = l_src[1];                             \
      *l_dest++ = l_src[2];                             \
      *l_dest++ = l_src[3];                             \
      *l_dest++ = l_src[4];                             \
      *l_dest++ = l_src[5];                             \
    } else if (nelms == 7) {                             \
      *l_dest++ = l_src[0];                             \
      *l_dest++ = l_src[1];                             \
      *l_dest++ = l_src[2];                             \
      *l_dest++ = l_src[3];                             \
      *l_dest++ = l_src[4];                             \
      *l_dest++ = l_src[5];                             \
      *l_dest++ = l_src[6];                             \
    } else if (nelms == 8) {                             \
      *l_dest++ = l_src[0];                             \
      *l_dest++ = l_src[1];                             \
      *l_dest++ = l_src[2];                             \
      *l_dest++ = l_src[3];                             \
      *l_dest++ = l_src[4];                             \
      *l_dest++ = l_src[5];                             \
      *l_dest++ = l_src[6];                             \
      *l_dest++ = l_src[7];                             \
    } else {                                            \
      tmp_src = l_src;                                  \
      j = (long)nelms;                                  \
      while (j >= 8) {                                   \
        *l_dest++ = tmp_src[0];                         \
        *l_dest++ = tmp_src[1];                         \
        *l_dest++ = tmp_src[2];                         \
        *l_dest++ = tmp_src[3];                         \
        *l_dest++ = tmp_src[4];                         \
        *l_dest++ = tmp_src[5];                         \
        *l_dest++ = tmp_src[6];                         \
        *l_dest++ = tmp_src[7];                         \
        j -= 8;                                         \
        tmp_src += 8;                                   \
      }                                                 \
      for(k = 0; k < j; k++) {                          \
        *l_dest++ = *tmp_src++;                         \
      }                                                 \
    }                                                   \
  }                                                     \


#endif /* VECCPY_H */

/*
 * Local variables:
 * c-indent-tabs-mode: nil
 * End:
 */
