/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef MPIR_DAME_H
#define MPIR_DAME_H

#include <mpi.h>

#define MPIR_DATALOOP_GET_FIELD(hetero_,value_,fieldname_) do {         \
    value_ = ((MPIR_Datatype *)ptr)->dataloop##fieldname_;              \
  } while (0)

#define MPIR_DATALOOP_SET_FIELD(hetero_,value_,fieldname_) do {         \
    ((MPIR_Datatype *)ptr)->dataloop##fieldname_ = value_;              \
  } while (0)

#define DAME_Offset     MPI_Aint
#define DAME_Count      MPI_Aint
#define DAME_Handle     MPI_Datatype
#define DAME_Type       MPI_Datatype
#define DAME_Buffer     void *
#define DAME_VECTOR     MPL_IOV
#define DAME_VECTOR_LEN MPL_IOV_LEN
#define DAME_VECTOR_BUF MPL_IOV_BUF
#define DAME_Size       MPI_Aint

/* These flags are used at creation time to specify what types of
 * optimizations may be applied. They are also passed in at Segment_init
 * time to specify which dataloop to use.
 *
 * Note: The flag to MPIR_Segment_init() was originally simply "hetero"
 * and was a boolean value (0 meaning homogeneous). Some MPICH code
 * may still rely on HOMOGENEOUS being "0" and HETEROGENEOUS being "1".
 */
#define DAME_DATALOOP_HOMOGENEOUS   0
#define DAME_DATALOOP_HETEROGENEOUS 1
#define DAME_DATALOOP_ALL_BYTES     2

/* USE THE NOTATION THAT BILL USED IN MPIIMPL.H AND MAKE THESE MACROS */

/* NOTE: put get size into mpiimpl.h; the others go here until such time
 * as we see that we need them elsewhere.
 */
#define DAME_Handle_get_loopdepth_macro(handle_,depth_)                       \
  MPIR_Datatype_get_loopdepth_macro(handle_, depth_, DAME_DATALOOP_HOMOGENEOUS)

#define DAME_Handle_get_loopsize_macro(handle_,size_)                         \
  MPIR_Datatype_get_loopsize_macro(handle_, size_, DAME_DATALOOP_HOMOGENEOUS)

#define DAME_Handle_get_loopptr_macro(handle_,lptr_)                          \
  MPIR_Datatype_get_loopptr_macro(handle_, lptr_, DAME_DATALOOP_HOMOGENEOUS)

#define DAME_Handle_set_loopptr_macro(handle_,lptr_)                          \
  MPIR_Datatype_set_loopptr_macro(handle_, lptr_, DAME_DATALOOP_HOMOGENEOUS)

#define DAME_Handle_set_loopdepth_macro(handle_,depth_)                 \
  MPIR_Datatype_set_loopdepth_macro(handle_, depth_, DAME_DATALOOP_HOMOGENEOUS)

#define DAME_Handle_set_loopsize_macro(handle_,size_)                   \
  MPIR_Datatype_set_loopsize_macro(handle_, size_, DAME_DATALOOP_HOMOGENEOUS)

#define DAME_Handle_get_size_macro(handle_,size_) \
  MPIR_Datatype_get_size_macro(handle_,size_)

#define DAME_Handle_get_basic_type_macro(handle_,basic_type_) \
    MPIR_Datatype_get_basic_type(handle_, basic_type_)

#define DAME_Handle_get_extent_macro(handle_,extent_) \
    MPIR_Datatype_get_extent_macro(handle_,extent_)

#define DAME_Handle_hasloop_macro(handle_)                           \
    ((HANDLE_GET_KIND(handle_) == HANDLE_KIND_BUILTIN) ? 0 : 1)

#define DAME_Ensure_Offset_fits_in_pointer(value_) \
    MPIR_Ensure_Aint_fits_in_pointer(value_)

/* allocate and free functions must also be defined. */
#define DAME_Malloc(x) MPL_malloc((x), MPL_MEM_DATATYPE)
#define DAME_Free   MPL_free

/* assert function */
#define DAME_Assert MPIR_Assert

/* memory copy function */
#define DAME_Memcpy MPIR_Memcpy

/* casting macros */
#define DAME_OFFSET_CAST_TO_VOID_PTR MPIR_AINT_CAST_TO_VOID_PTR
#define DAME_VOID_PTR_CAST_TO_OFFSET MPIR_VOID_PTR_CAST_TO_MPI_AINT
#define DAME_PTR_DISP_CAST_TO_OFFSET MPIR_PTR_DISP_CAST_TO_MPI_AINT

/* Redefine all of the internal structures in terms of the prefix */
#define DAME_Dame              MPIR_Dataloop
#define DAME_Dame_contig       MPIR_Dataloop_contig
#define DAME_Dame_vector       MPIR_Dataloop_vector
#define DAME_Dame_blockindexed MPIR_Dataloop_blockindexed
#define DAME_Dame_indexed      MPIR_Dataloop_indexed
#define DAME_Dame_struct       MPIR_Dataloop_struct
#define DAME_Dame_common       MPIR_Dataloop_common
#define DAME_Segment           MPIR_Segment
#define DAME_Stackelm          MPIR_Dataloop_stackelm

#define DAME_PACKUNPACK_COMPLETED 0
#define DAME_PACKUNPACK_PARTIAL 1
#define DAME_PACKUNPACK_ERROR -1

// Forward declaration
DAME_VECTOR;

typedef enum {
    DL_BOTTOM = -1,
    DL_EXIT = 0,

    DL_CONTIG = 1,
    DL_CONTIGFINAL = 11,

    DL_VECTOR = 2,
    DL_VECTOR1 = 3,
    DL_VECTORFINAL = 21,
    DL_VECFINAL_1 = 22,
    DL_VECFINAL_2 = 23,
    DL_VECFINAL_4 = 24,
    DL_VECFINAL_8 = 25,

    DL_BLOCKINDEX = 4,
    DL_BLOCKINDEX1 = 5,
    DL_BLOCKINDEXFINAL = 41,
    DL_BLKINDEXFINAL_1 = 42,
    DL_BLKINDEXFINAL_2 = 43,
    DL_BLKINDEXFINAL_4 = 44,
    DL_BLKINDEXFINAL_8 = 45,

    DL_INDEX = 6,
    DL_INDEXFINAL = 61,

    DL_STRUCT = 7,

    DL_CONTIGCHILD = 8,

    DL_RETURNTO = 31
} Dame_kind;

// Forward declaration
struct DAME_Dame;

typedef struct DAME_Dame_contig {
    MPI_Aint basesize;
    MPI_Aint baseextent;
} DAME_Dame_contig;

typedef struct DAME_Dame_vector {
    MPI_Aint oldsize;
    MPI_Aint blklen;            /* This is the blocklength in #elements */
    MPI_Aint stride;
} DAME_Dame_vector;

typedef struct DAME_Dame_blockindexed {
    MPI_Aint oldsize;
    MPI_Aint blklen;
    const MPI_Aint *offsets;
} DAME_Dame_blockindexed;

typedef struct DAME_Dame_indexed {
    MPI_Aint oldsize;
    const MPI_Aint *blklens;
    const MPI_Aint *offsets;
} DAME_Dame_indexed;

typedef struct DAME_Dame_struct {
    const MPI_Aint *oldsizes;
    const MPI_Aint *blklens;
    const MPI_Aint *offsets;
    struct DAME_Dame **dls;
} DAME_Dame_struct;

typedef struct DAME_Dame {
    Dame_kind kind;
    MPI_Aint count;
    MPI_Aint size;
    MPI_Aint extent;
    /* This is used when operating on a struct type. It indicates where in the
     * stack the control should return */
    MPI_Aint returnto;
    /* Additional information about the dataloop, potentially to enable
     * optimizations */
    unsigned flags;
    union {
        DAME_Dame_contig c_t;
        DAME_Dame_vector v_t;
        DAME_Dame_blockindexed bi_t;
        DAME_Dame_indexed i_t;
        DAME_Dame_struct s_t;
    } s;
} DAME_Dame;

#define DLOOP_FINAL_MASK  0x00000008
#define DLOOP_KIND_MASK   0x00000007
#define DLOOP_KIND_CONTIG 0x1
#define DLOOP_KIND_VECTOR 0x2
#define DLOOP_KIND_BLOCKINDEXED 0x3
#define DLOOP_KIND_INDEXED 0x4
#define DLOOP_KIND_STRUCT 0x5

/* The max datatype depth is the maximum depth of the stack used to
   evaluate datatypes.  It represents the length of the chain of
   datatype dependencies.  Defining this and testing when a datatype
   is created removes a test in the datatype evaluation loop. */
#define DAME_MAX_DEPTH 64
#define DAME_MAX_STACK_SIZE DAME_MAX_DEPTH
typedef struct Dame_Stackelm {
    /* This is the base address used to keep track of
     * the src pointer when packing and the dst pointer when unpacking */
    MPI_Offset base;
    /* When we encounter a struct type, we push the current dataloop here
     * and use the dataloop of the struct element being processed */
    const DAME_Dame *prevdl;
    /* For the CONTIGFINAL type, this contains the number of complete elements
     * that are left to be copied (packed/unpacked).
     * bytes = stack[i].countLeft*c_t.basesize + state.partial;
     * is the number of bytes remaining in the type to be processed
     * For all other types, this is simply the number of elements left
     */
    MPI_Aint countLeft;
} DAME_Stackelm;

typedef struct DAME_State {
    /* Pointer in dataloop stack to which to resume */
    long sp;
    /* Actual dataloop of type to resume */
    const DAME_Dame *dl;
    /* For the CONTIGFINAL type, this gives the number of bytes in a partially
     * copied element that are yet to be copied */
    DAME_Offset partial;
    DAME_Stackelm stack[DAME_MAX_STACK_SIZE];
} DAME_State;


/* Optimization flags */
#define DAME_OPT_FLAG_ALIGNED 0x1
#define DAME_OPT_FLAG_ISSHORT 0x2

/* When copying blocks larger than this, use memcpy.
 * TODO: Provide a way of configuring this useing ./configure */
#define DAME_MEMCPY_THRESHOLD 1024

/* Helper macros */
#define Dame_is_contig(dl)                              \
    ((dl)[1].kind == DL_CONTIGFINAL)                    \

#define set_dlflag(fvec, f)                   \
    ((fvec) | ~(f))                           \

#define get_dlflag(fvec, f)                   \
    ((fvec) & (f))                            \

#define DAME_opt_get_aligned(dl) get_dlflag((dl).flags, DAME_OPT_FLAG_ALIGNED)
#define DAME_opt_set_aligned(dl) set_dlflag((dl).flags, DAME_OPT_FLAG_ALIGNED)

#define DAME_opt_get_isshort(dl) get_dlflag((dl).flags, DAME_OPT_FLAG_ISSHORT)
#define DAME_opt_set_isshort(dl) set_dlflag((dl).flags, DAME_OPT_FLAG_ISSHORT)

/* Redefine the Dame functions to have the dataloop names so that the
   functions in src/mpi/datatypes don't need to be changed */
#define MPIR_Dame_alloc           MPIR_Dataloop_alloc
#define MPIR_Dame_alloc_and_copy  MPIR_Dataloop_alloc_and_copy
#define MPIR_Dame_struct_alloc    MPIR_Dataloop_struct_alloc
#define MPIR_Dame_dup             MPIR_Dataloop_dup
#define MPIR_Dame_free            MPIR_Dataloop_free
#define MPIR_Dame_move            MPIR_Dataloop_move
#define MPIR_Dame_update          MPIR_Dataloop_update
#define MPIR_Dame_serialize       MPIR_Dataloop_serialize
#define MPIR_Dame_deserialize     MPIR_Dataloop_deserialize

/* Dame functions (dame.c) */
void MPIR_Dame_copy(void *dest, void *src, DAME_Size size);
DAME_Offset MPIR_Dame_stream_size(DAME_Dame * dl);
void MPIR_Dame_print(DAME_Dame * dataloop);
void MPIR_Dame_print_compact(DAME_Dame * compact_dataloop);

int MPIR_Dame_alloc(DAME_Dame ** new_loop_p);
void MPIR_Dame_alloc_and_copy(DAME_Dame * old_loop, DAME_Dame ** new_loop_p);
int MPIR_Dame_struct_alloc(DAME_Count count, DAME_Dame * dl);
int MPIR_Dame_dup(DAME_Dame * old_loop, DAME_Size unused, DAME_Dame ** new_loop_p);

void MPIR_Dame_free(DAME_Dame ** dataloop, int is_compact);
int MPIR_Dame_move(DAME_Dame * old, DAME_Dame * new);

int MPIR_Dame_update(DAME_Dame * dataloop, MPI_Aint depth);
void MPIR_Dame_calculate_size(DAME_Dame * dl, MPI_Aint * dl_size, int only_primitives);

int MPIR_Dame_serialize(DAME_Dame * dl, MPI_Aint dataloop_size, DAME_Dame ** compact_dataloop);
int MPIR_Dame_deserialize(DAME_Dame * compact_dl,
                          MPI_Aint dataloop_size, DAME_Dame ** new_compact_dl, DAME_Dame ** new_dl);


/* NOTE: ASSUMING LAST TYPE IS SIGNED */
#define SEGMENT_IGNORE_LAST ((DAME_Offset) -1)


#define SEGMENT_PACK_IDX 0
#define SEGMENT_PACK_IOV_IDX 1
#define SEGMENT_UNPACK_IDX 2

/*S
  DAME_Segment - Description of the Segment datatype

  Notes:
  This has no corresponding MPI datatype.

  Module:
  Segment

  Questions:
  Should this have an id for allocation and similarity purposes?
  S*/
typedef struct DAME_Segment {
    // The same segment object could be used to pack and IO vector and later
    // to unpack a type. One of them - could be a partial pack which would
    // then result in the resumption of a different operation. So we create
    // separate state objects for each possible use of the segment. Similarly
    // for some of the other fields
    void *ptr;                  /* pointer to datatype buffer */
    DAME_Handle handle;
    DAME_Count count;
    DAME_Offset in_offset[3];
    DAME_Offset out_offset[3];
    DAME_Count curcount[3];
    DAME_State state[3];
} DAME_Segment;

/* Segment functions (segment.c) */
DAME_Segment *MPIR_Segment_alloc(void);

void MPIR_Segment_free(DAME_Segment * segp);

int MPIR_Segment_init(const DAME_Buffer buf,
                      DAME_Count count, DAME_Handle handle, DAME_Segment * segp, int hetero);

int DL_pack(const void *inPtr, void *outPtr, DAME_Offset outSize,
            const DAME_Dame * dl, int isbuiltin, DAME_State * state_p, DAME_Offset * copy_p);

int DL_pack_iov(const void *inPtr, DAME_VECTOR * outPtr,
                DAME_Offset outsize, int outSlots, int merge,
                const DAME_Dame * dl,
                int isbuiltin, DAME_State * state_p, DAME_Offset * copy_p, int *iovs_p);

int DL_unpack(const void *inPtr, void *outPtr, DAME_Offset outSize,
              const DAME_Dame * dl, int isbuiltin, DAME_State * state_p, DAME_Offset * copy_p);

int DL_unpack_iov(const void *inPtr, DAME_VECTOR * outPtr,
                  DAME_Offset outsize, int outSlots,
                  const DAME_Dame * dl,
                  int isbuiltin, DAME_State * state_p, DAME_Offset * copy_p, int *iovs_p);

/* Common segment operations (segment_ops.c) */
void MPIR_Segment_count_contig_blocks(DAME_Segment * segp,
                                      DAME_Offset first, DAME_Offset * lastp, DAME_Count * countp);
void MPIR_Segment_pack(struct DAME_Segment *segp,
                       DAME_Offset first, DAME_Offset * lastp, void *streambuf);
void MPIR_Segment_unpack(struct DAME_Segment *segp,
                         DAME_Offset first, DAME_Offset * lastp, void *streambuf);


#define MPIR_DATALOOP_HOMOGENOUS    DAME_DATALOOP_HOMOGENEOUS
#define MPIR_DATALOOP_HETEROGENEOUS DAME_DATALOOP_HETEROGENEOUS
#define MPIR_DATALOOP_ALL_BYTES     DAME_DATALOOP_ALL_BYTES

void MPIR_Dame_create(MPI_Datatype type, DAME_Dame ** dlp_p, MPI_Aint * size, int *depth);

int MPIR_Dame_create_contiguous(DAME_Count count,
                                MPI_Datatype oldtype, DAME_Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_vector(DAME_Count count,
                            DAME_Size blocklength,
                            MPI_Aint stride,
                            int strideinbytes,
                            MPI_Datatype type,
                            MPI_Datatype oldtype, DAME_Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_blockindexed(DAME_Count count,
                                  DAME_Size blklen,
                                  const void *disp_array,
                                  int dispinbytes,
                                  MPI_Datatype type,
                                  MPI_Datatype oldtype,
                                  DAME_Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_indexed(DAME_Count count,
                             const DAME_Size * blklen,
                             const void *disp_array,
                             int dispinbytes,
                             MPI_Datatype type,
                             MPI_Datatype oldtype, DAME_Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_struct(DAME_Count count,
                            const int *blklen_array,
                            const MPI_Aint * disp_array,
                            MPI_Datatype type,
                            const MPI_Datatype * oldtype_array,
                            DAME_Dame ** dl, MPI_Aint * size, int *depth);

/* we bump up the size of the blocklength array because create_struct might use
 * create_indexed in an optimization, and in course of doing so, generate a
 * request of a large blocklength. */
int MPIR_Dame_create_pairtype(MPI_Datatype type, DAME_Dame ** dl, MPI_Aint * size, int *depth);

/* Helper functions for dataloop construction */
int MPIR_Type_convert_subarray(int ndims,
                               int *array_of_sizes,
                               int *array_of_subsizes,
                               int *array_of_starts,
                               int order, MPI_Datatype oldtype, MPI_Datatype * newtype);
int MPIR_Type_convert_darray(int size,
                             int rank,
                             int ndims,
                             int *array_of_gsizes,
                             int *array_of_distribs,
                             int *array_of_dargs,
                             int *array_of_psizes,
                             int order, MPI_Datatype oldtype, MPI_Datatype * newtype);

DAME_Count MPIR_Type_blockindexed_count_contig(DAME_Count count,
                                               DAME_Count blklen,
                                               const void *disp_array,
                                               int dispinbytes, DAME_Offset old_extent);

DAME_Count MPIR_Type_indexed_count_contig(DAME_Count count,
                                          const DAME_Size * blocklength_array,
                                          const void *displacement_array,
                                          int dispinbytes, DAME_Offset old_extent);

int MPIR_Type_blockindexed(int count,
                           int blocklength,
                           const void *displacement_array,
                           int dispinbytes, MPI_Datatype oldtype, MPI_Datatype * newtype);

int MPIR_Type_commit(MPI_Datatype * type);

/* Segment functions specific to MPICH */
void MPIR_Segment_pack_vector(struct DAME_Segment *segp,
                              DAME_Offset first,
                              DAME_Offset * lastp, DAME_VECTOR * vector, int *lengthp);

void MPIR_Segment_unpack_vector(struct DAME_Segment *segp,
                                DAME_Offset first,
                                DAME_Offset * lastp, DAME_VECTOR * vector, int *lengthp);

void MPIR_Segment_flatten(struct DAME_Segment *segp,
                          DAME_Offset first,
                          DAME_Offset * lastp,
                          DAME_Offset * offp, DAME_Size * sizep, DAME_Offset * lengthp);

void MPIR_Segment_pack_external32(struct DAME_Segment *segp,
                                  DAME_Offset first, DAME_Offset * lastp, void *pack_buffer);

void MPIR_Segment_unpack_external32(struct DAME_Segment *segp,
                                    DAME_Offset first,
                                    DAME_Offset * lastp, DAME_Buffer unpack_buffer);

/* The RMA code still refers to these with the dloop prefix */

#define DLOOP_Offset     DAME_Offset
#define DLOOP_Count      DAME_Count
#define DLOOP_Handle     DAME_Handle
#define DLOOP_Type       DAME_Type
#define DLOOP_Buffer     DAME_Buffer
#define DLOOP_VECTOR     DAME_VECTOR
#define DLOOP_VECTOR_LEN DAME_VECTOR_LEN
#define DLOOP_VECTOR_BUF DAME_VECTOR_BUF
#define DLOOP_Size       DAME_Size

#define DLOOP_Segment    DAME_Segment

#define DLOOP_Assert     DAME_Assert

#endif
