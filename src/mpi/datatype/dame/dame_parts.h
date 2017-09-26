/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef DAME_PARTS_H
#define DAME_PARTS_H

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
struct Dame;

typedef struct Dame_contig {
    MPI_Aint basesize;
    MPI_Aint baseextent;
} Dame_contig;

typedef struct Dame_vector {
    MPI_Aint oldsize;
    MPI_Aint blklen;            /* This is the blocklength in #elements */
    MPI_Aint stride;
} Dame_vector;

typedef struct Dame_blockindexed {
    MPI_Aint oldsize;
    MPI_Aint blklen;
    const MPI_Aint *offsets;
} Dame_blockindexed;

typedef struct Dame_indexed {
    MPI_Aint oldsize;
    const MPI_Aint *blklens;
    const MPI_Aint *offsets;
} Dame_indexed;

typedef struct Dame_struct {
    const MPI_Aint *oldsizes;
    const MPI_Aint *blklens;
    const MPI_Aint *offsets;
    struct Dame **dls;
} Dame_struct;

typedef struct Dame {
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
        Dame_contig c_t;
        Dame_vector v_t;
        Dame_blockindexed bi_t;
        Dame_indexed i_t;
        Dame_struct s_t;
    } s;
} Dame;

/* The max datatype depth is the maximum depth of the stack used to
   evaluate datatypes.  It represents the length of the chain of
   datatype dependencies.  Defining this and testing when a datatype
   is created removes a test in the datatype evaluation loop. */
#define DAME_MAX_DEPTH 64
#define DAME_MAX_STACK_SIZE DAME_MAX_DEPTH
typedef struct Dame_stackelm {
    /* This is the base address used to keep track of
     * the src pointer when packing and the dst pointer when unpacking */
    MPI_Offset base;
    /* When we encounter a struct type, we push the current dataloop here
     * and use the dataloop of the struct element being processed */
    const Dame *prevdl;
    /* For the CONTIGFINAL type, this contains the number of complete elements
     * that are left to be copied (packed/unpacked).
     * bytes = stack[i].countLeft*c_t.basesize + state.partial;
     * is the number of bytes remaining in the type to be processed
     * For all other types, this is simply the number of elements left
     */
    MPI_Aint countLeft;
} Dame_stackelm;

typedef struct Dame_state {
    /* Pointer in dataloop stack to which to resume */
    long sp;
    /* Actual dataloop of type to resume */
    const Dame *dl;
    /* For the CONTIGFINAL type, this gives the number of bytes in a partially
     * copied element that are yet to be copied */
    DAME_Offset partial;
    Dame_stackelm stack[DAME_MAX_STACK_SIZE];
} Dame_state;


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


/* Dataloop functions (dataloop.c) */
void MPIR_Dataloop_copy(void *dest, void *src, DAME_Size size);
DAME_Offset MPIR_Dataloop_stream_size(Dame * dl_p, DAME_Offset(*sizefn) (DAME_Type el_type));
void MPIR_Dataloop_print(Dame * dataloop);
void MPIR_Dataloop_print_compact(Dame * compact_dataloop);

int MPIR_Dataloop_alloc(Dame ** new_loop_p);
void MPIR_Dataloop_alloc_and_copy(Dame * old_loop, Dame ** new_loop_p);
int MPIR_Dataloop_struct_alloc(DAME_Count count, Dame * dl);
int MPIR_Dataloop_dup(Dame * old_loop, Dame ** new_loop_p);

void MPIR_Dataloop_free(Dame ** dataloop);
int MPIR_Dataloop_move(Dame * old, Dame * new);
int MPIR_Dataloop_update(Dame * dataloop, MPI_Aint depth);
void MPIR_Dataloop_calculate_size(struct Dame *dl, MPI_Aint * dl_size, int only_primitives);
int MPIR_Dataloop_serialize(Dame * dl, MPI_Aint dataloop_size, Dame ** compact_dataloop);
int MPIR_Dataloop_deserialize(Dame * compact_dl,
                              MPI_Aint dataloop_size, Dame ** new_compact_dl, Dame ** new_dl);



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
    Dame_state state[3];
} DAME_Segment;

/* Segment functions (segment.c) */
DAME_Segment *MPIR_Segment_alloc(void);

void MPIR_Segment_free(DAME_Segment * segp);

int MPIR_Segment_init(const DAME_Buffer buf,
                      DAME_Count count, DAME_Handle handle, DAME_Segment * segp, int hetero);

int DL_pack(const void *inPtr, void *outPtr, DAME_Offset outSize,
            const Dame * dl, int isbuiltin, Dame_state * state_p, DAME_Offset * copy_p);

int DL_pack_iov(const void *inPtr, DAME_VECTOR * outPtr,
                DAME_Offset outsize, int outSlots, int merge,
                const Dame * dl,
                int isbuiltin, Dame_state * state_p, DAME_Offset * copy_p, int *iovs_p);

int DL_unpack(const void *inPtr, void *outPtr, DAME_Offset outSize,
              const Dame * dl, int isbuiltin, Dame_state * state_p, DAME_Offset * copy_p);

int DL_unpack_iov(const void *inPtr, DAME_VECTOR * outPtr,
                  DAME_Offset outsize, int outSlots,
                  const Dame * dl,
                  int isbuiltin, Dame_state * state_p, DAME_Offset * copy_p, int *iovs_p);

/* Common segment operations (segment_ops.c) */
void MPIR_Segment_count_contig_blocks(DAME_Segment * segp,
                                      DAME_Offset first, DAME_Offset * lastp, DAME_Count * countp);
void MPIR_Segment_pack(struct DAME_Segment *segp,
                       DAME_Offset first, DAME_Offset * lastp, void *streambuf);
void MPIR_Segment_unpack(struct DAME_Segment *segp,
                         DAME_Offset first, DAME_Offset * lastp, void *streambuf);

void MPIR_Dame_create(MPI_Datatype type, Dame ** dlp_p, MPI_Aint * size, int *depth);

int MPIR_Dame_create_contiguous(DAME_Count count,
                                MPI_Datatype oldtype, Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_vector(DAME_Count count,
                            DAME_Size blocklength,
                            MPI_Aint stride,
                            int strideinbytes,
                            MPI_Datatype type,
                            MPI_Datatype oldtype, Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_blockindexed(DAME_Count count,
                                  DAME_Size blklen,
                                  const void *disp_array,
                                  int dispinbytes,
                                  MPI_Datatype type,
                                  MPI_Datatype oldtype, Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_indexed(DAME_Count count,
                             const DAME_Size * blklen,
                             const void *disp_array,
                             int dispinbytes,
                             MPI_Datatype type,
                             MPI_Datatype oldtype, Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_struct(DAME_Count count,
                            const int *blklen_array,
                            const MPI_Aint * disp_array,
                            MPI_Datatype type,
                            const MPI_Datatype * oldtype_array,
                            Dame ** dl, MPI_Aint * size, int *depth);

/* we bump up the size of the blocklength array because create_struct might use
 * create_indexed in an optimization, and in course of doing so, generate a
 * request of a large blocklength. */
int MPIR_Dame_create_pairtype(MPI_Datatype type, Dame ** dl, MPI_Aint * size, int *depth);

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

#endif
