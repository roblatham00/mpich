/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*-
 *   vim: ts=8 sts=4 sw=4 noexpandtab
 *
 * This is adapted from Avery Ching's "new List I/O" code -- one third the
 * number of lines as the "classic" list i/o
 */

#include <assert.h>
#include "adio.h"
#include "adio_extern.h"
#include "ad_mochio.h"
#include "ad_mochio_common.h"

#if 0
#define DEBUG_LIST 1
#define DEBUG_LIST2 1
#endif

#define COALESCE_REGIONS        /* TODO: would we ever want to *not* coalesce? */
#define MAX_OL_COUNT 1024

/* To do: Fix the code to coalesce the offset-length pairs for memory
 * and file. */

/* gen_listio_arr - fills in offset-length pairs for memory and file
 * for list I/O */
static int gen_listio_arr(ADIOI_Flatlist_node * flat_buf_p,     /* IN: memory description */
                          int *flat_buf_index_p,        /* IN: partial-progress cursor */
                          int64_t * cur_flat_buf_reg_off_p,     /* ??? */
                          int flat_buf_size,    /* IN: total bytes described in memory */
                          int flat_buf_extent,  /* IN: extent of memory description */
                          ADIOI_Flatlist_node * flat_file_p,    /* IN: file description */
                          int *flat_file_index_p,       /* IN: partial-progress cursor */
                          int64_t * cur_flat_file_reg_off_p,    /* ??? */
                          int flat_file_size,   /* IN: total bytes of file description */
                          int flat_file_extent, /* IN: extent of file description */
                          int max_ol_count,     /* IN: largest size offset-length arrays to make */
                          ADIO_Offset disp,     /* IN: displacement in file */
                          int bytes_into_filetype,      /* IN: partial progress: how far into file type to resume */
                          int64_t * bytes_completed,    /* OUT: bytes described in this collection of offset-length arrays */
                          int64_t total_io_size,        /* ??? */
                          void *buffer, /* IN: starting address for relative offsets */
                          char *buf_addr_arr[], /* OUT: array of memory offsets */
                          uint64_t buf_len_arr[],       /* OUT: array of memory lengths */
                          int64_t * buf_ol_count_p,     /* OUT: number of elements of memory offset/length arrays */
                          off_t file_off_arr[], uint64_t file_len_arr[], int64_t * file_ol_count_p)
{       /* OUT: same but for file arrays */
    int region_size = -1;

    /* parameters for flattened memory and file datatypes */
    int64_t cur_flat_buf_reg_left = 0;
    int64_t cur_flat_file_reg_left = 0;

#ifdef DEBUG_LIST2
    fprintf(stderr, "gen_list_arr:\n");
#endif

    if ((*buf_ol_count_p) != 0 || (*file_ol_count_p) != 0) {
        fprintf(stderr, "buf_ol_count != 0 || file_ol_count != 0\n");
        return -1;
    }

    /* Start on a non-zero memory and file region
     * Note this does not affect the bytes_completed
     * since no data is in these regions.  Initialize the
     * first memory and file offsets. */
    while (flat_buf_p->blocklens[(*flat_buf_index_p)] == 0) {
        (*flat_buf_index_p) = ((*flat_buf_index_p) + 1) % flat_buf_p->count;
    }
    buf_addr_arr[*buf_ol_count_p] = (char *) buffer +
        ((*bytes_completed / flat_buf_size) *
         flat_buf_extent + flat_buf_p->indices[*flat_buf_index_p] + *cur_flat_buf_reg_off_p);
    buf_len_arr[*buf_ol_count_p] = 0;

    while (flat_file_p->blocklens[(*flat_file_index_p)] == 0) {
        (*flat_file_index_p) = ((*flat_file_index_p) + 1) % flat_file_p->count;
    }
    file_off_arr[*file_ol_count_p] = disp +
        (((bytes_into_filetype + *bytes_completed) / flat_file_size) *
         flat_file_extent) + flat_file_p->indices[*flat_file_index_p] + *cur_flat_file_reg_off_p;
    file_len_arr[*file_ol_count_p] = 0;

#ifdef DEBUG_LIST2
    fprintf(stderr, "initial buf_addr_arr[%ld] = %p\n", *buf_ol_count_p,
            buf_addr_arr[*buf_ol_count_p]);
    fprintf(stderr, "initial file_off_arr[%ld] = %ld\n", *file_ol_count_p,
            file_off_arr[*file_ol_count_p]);
#endif

    while (*bytes_completed != total_io_size
           && (*buf_ol_count_p) < max_ol_count && (*file_ol_count_p) < max_ol_count) {
        /* How much data is left in the current piece in
         * the flattened datatypes */
        cur_flat_buf_reg_left = flat_buf_p->blocklens[*flat_buf_index_p]
            - *cur_flat_buf_reg_off_p;
        cur_flat_file_reg_left = flat_file_p->blocklens[*flat_file_index_p]
            - *cur_flat_file_reg_off_p;

#ifdef DEBUG_LIST2
        fprintf(stderr,
                "flat_buf_index=%d flat_buf->blocklens[%d]=%lld\n"
                "cur_flat_buf_reg_left=%ld "
                "*cur_flat_buf_reg_off_p=%ld\n"
                "flat_file_index=%d flat_file->blocklens[%d]=%lld\n"
                "cur_flat_file_reg_left=%ld "
                "*cur_flat_file_reg_off_p=%ld\n"
                "bytes_completed=%ld\n"
                "buf_ol_count=%ld file_ol_count=%ld\n"
                "buf_addr_arr[%ld]=%p file_len_arr[%ld]=%ld\n\n",
                *flat_buf_index_p, *flat_buf_index_p,
                flat_buf_p->blocklens[*flat_buf_index_p],
                cur_flat_buf_reg_left,
                *cur_flat_buf_reg_off_p,
                *flat_file_index_p, *flat_file_index_p,
                flat_file_p->blocklens[*flat_file_index_p],
                cur_flat_file_reg_left,
                *cur_flat_file_reg_off_p,
                *bytes_completed,
                *buf_ol_count_p, *file_ol_count_p,
                *buf_ol_count_p,
                buf_addr_arr[*buf_ol_count_p], *file_ol_count_p, file_len_arr[*file_ol_count_p]);
#endif

        /* What is the size of the next contiguous region agreed
         * upon by both memory and file regions that does not
         * surpass the file size */
        if (cur_flat_buf_reg_left > cur_flat_file_reg_left)
            region_size = cur_flat_file_reg_left;
        else
            region_size = cur_flat_buf_reg_left;

        if (region_size > total_io_size - *bytes_completed)
            region_size = total_io_size - *bytes_completed;

        /* Add this piece to both the mem and file arrays
         * coalescing offset-length pairs if possible and advance
         * the pointers through the flatten mem and file datatypes
         * as well Note: no more than a single piece can be done
         * since we take the smallest one possible */

        if (cur_flat_buf_reg_left == region_size) {
#ifdef DEBUG_LIST2
            fprintf(stderr, "reached end of memory block...\n");
#endif
            (*flat_buf_index_p) = ((*flat_buf_index_p) + 1) % flat_buf_p->count;
            while (flat_buf_p->blocklens[(*flat_buf_index_p)] == 0) {
                (*flat_buf_index_p) = ((*flat_buf_index_p) + 1) % flat_buf_p->count;
            }
            *cur_flat_buf_reg_off_p = 0;

#ifdef COALESCE_REGIONS
            if (*buf_ol_count_p != 0) {
                if (buf_addr_arr[(*buf_ol_count_p) - 1] +
                    buf_len_arr[(*buf_ol_count_p) - 1] == buf_addr_arr[*buf_ol_count_p]) {
                    buf_len_arr[(*buf_ol_count_p) - 1] += region_size;
                } else {
                    buf_len_arr[*buf_ol_count_p] += region_size;
                    (*buf_ol_count_p)++;
                }
            } else {
#endif
                buf_len_arr[*buf_ol_count_p] += region_size;
                (*buf_ol_count_p)++;
#ifdef COALESCE_REGIONS
            }
#endif

            /* Don't prepare for the next piece if we have reached
             * the limit or else it will segment fault. */
            if ((*buf_ol_count_p) != max_ol_count) {
                buf_addr_arr[*buf_ol_count_p] = (char *) buffer +
                    (((*bytes_completed + region_size) / flat_buf_size) *
                     flat_buf_extent +
                     flat_buf_p->indices[*flat_buf_index_p] + (*cur_flat_buf_reg_off_p));
                buf_len_arr[*buf_ol_count_p] = 0;
            }
        } else if (cur_flat_buf_reg_left > region_size) {
#ifdef DEBUG_LIST2
            fprintf(stderr, "advanced %d in memory block...\n", region_size);
#endif
            (*cur_flat_buf_reg_off_p) += region_size;
            buf_len_arr[*buf_ol_count_p] += region_size;
        } else {
            fprintf(stderr, "gen_listio_arr: Error\n");
        }

        /* To calculate the absolute file offset we need to
         * add the disp, how many filetypes we have gone through,
         * the relative block offset in the filetype and how far
         * into the block we have gone. */
        if (cur_flat_file_reg_left == region_size) {
#ifdef DEBUG_LIST2
            fprintf(stderr, "reached end of file block...\n");
#endif
            (*flat_file_index_p) = ((*flat_file_index_p) + 1) % flat_file_p->count;
            while (flat_file_p->blocklens[(*flat_file_index_p)] == 0) {
                (*flat_file_index_p) = ((*flat_file_index_p) + 1) % flat_file_p->count;
            }
            (*cur_flat_file_reg_off_p) = 0;

#ifdef COALESCE_REGIONS
            if (*file_ol_count_p != 0) {
                if (file_off_arr[(*file_ol_count_p) - 1] +
                    file_len_arr[(*file_ol_count_p) - 1] == file_off_arr[*file_ol_count_p]) {
                    file_len_arr[(*file_ol_count_p) - 1] += region_size;
                } else {
                    file_len_arr[*file_ol_count_p] += region_size;
                    (*file_ol_count_p)++;
                }
            } else {
#endif
                file_len_arr[*file_ol_count_p] += region_size;
                (*file_ol_count_p)++;
#ifdef COALESCE_REGIONS
            }
#endif

            /* Don't prepare for the next piece if we have reached
             * the limit or else it will segment fault. */
            if ((*file_ol_count_p) != max_ol_count) {
                file_off_arr[*file_ol_count_p] = disp +
                    (((bytes_into_filetype + *bytes_completed + region_size)
                      / flat_file_size) *
                     flat_file_extent) +
                    flat_file_p->indices[*flat_file_index_p] + (*cur_flat_file_reg_off_p);
                file_len_arr[*file_ol_count_p] = 0;
            }
        } else if (cur_flat_file_reg_left > region_size) {
#ifdef DEBUG_LIST2
            fprintf(stderr, "advanced %d in file block...\n", region_size);
#endif
            (*cur_flat_file_reg_off_p) += region_size;
            file_len_arr[*file_ol_count_p] += region_size;
        } else {
            fprintf(stderr, "gen_listio_arr: Error\n");
        }
#ifdef DEBUG_LIST2
        fprintf(stderr, "------------------------------\n\n");
#endif
        *bytes_completed += region_size;
    }
    /* Increment the count if we stopped in the middle of a
     * memory or file region */
    if (*cur_flat_buf_reg_off_p != 0)
        (*buf_ol_count_p)++;
    if (*cur_flat_file_reg_off_p != 0)
        (*file_ol_count_p)++;

    return 0;
}

#ifdef DEBUG_LIST2
static void print_buf_file_ol_pairs(char *buf_off_arr[],
                                    uint64_t buf_len_arr[],
                                    uint64_t buf_ol_count,
                                    int64_t file_off_arr[],
                                    uint64_t file_len_arr[], int64_t file_ol_count, void *buf,
                                    int rw_type)
{
    int i = -1;

    fprintf(stderr, "buf_ol_pairs(offset,length) count = %ld\n", buf_ol_count);
    for (i = 0; i < buf_ol_count; i++) {
        fprintf(stderr, "(%lld, %ld) ", (long long) buf_off_arr[i], buf_len_arr[i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "file_ol_pairs(offset,length) count = %ld\n", file_ol_count);
    for (i = 0; i < file_ol_count; i++) {
        fprintf(stderr, "(%lld, %ld) ", (long long) file_off_arr[i], file_len_arr[i]);
    }
    fprintf(stderr, "\n\n");

}
#endif
int ADIOI_MOCHIO_StridedListIO(ADIO_File fd, void *buf, int count,
                               MPI_Datatype datatype, int file_ptr_type,
                               ADIO_Offset offset, ADIO_Status * status,
                               int *error_code, int rw_type)
{
    /* list I/O parameters */
    int i = -1;
    int tmp_filetype_size = -1;
    int64_t cur_io_size = 0, io_size = 0;
    int etype_size = -1;
    int num_etypes_in_filetype = -1, num_filetypes = -1;
    int etypes_in_filetype = -1, size_into_filetype = -1;
    int bytes_into_filetype = 0;
    MPI_Offset response = 0, total_bytes_accessed = 0;

    /* parameters for offset-length pairs arrays */
    char *buf_off_arr[MAX_OL_COUNT];
    uint64_t buf_len_arr[MAX_OL_COUNT];
    off_t file_off_arr[MAX_OL_COUNT];
    uint64_t file_len_arr[MAX_OL_COUNT];
    int64_t buf_ol_count = 0;
    int64_t file_ol_count = 0;

    /* parameters for flattened memory and file datatypes */
    int flat_buf_index = 0;
    int flat_file_index = 0;
    int64_t cur_flat_buf_reg_off = 0;
    int64_t cur_flat_file_reg_off = 0;
    ADIOI_Flatlist_node *flat_buf_p, *flat_file_p;
    MPI_Count buftype_size = -1, filetype_size = -1;
    MPI_Aint filetype_extent = -1, buftype_extent = -1;;

    static char myname[] = "ADIOI_MOCHIO_STRIDED_LISTIO";

    if (fd->atomicity) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           myname, __LINE__,
                                           MPI_ERR_ARG,
                                           "Atomic noncontiguous writes"
                                           " are not supported by MOCHIO", 0);
        return -1;
    }

    MPI_Type_size_x(fd->filetype, &filetype_size);
    if (filetype_size == 0) {
        *error_code = MPI_SUCCESS;
        return -1;
    }
    MPI_Type_extent(fd->filetype, &filetype_extent);
    MPI_Type_size_x(datatype, &buftype_size);
    MPI_Type_extent(datatype, &buftype_extent);
    io_size = buftype_size * count;


    flat_buf_p = ADIOI_Flatten_and_find(datatype);
    flat_file_p = ADIOI_Flatten_and_find(fd->filetype);

    /* Find out where we are in the flattened filetype (the block index,
     * how far into the block, and how many bytes_into_filetype)
     * If the file_ptr_type == ADIO_INDIVIDUAL we will use disp, fp_ind
     * to figure this out (offset should always be zero)
     * If file_ptr_type == ADIO_EXPLICIT, we will use disp and offset
     * to figure this out. */

    etype_size = fd->etype_size;
    num_etypes_in_filetype = filetype_size / etype_size;
    if (file_ptr_type == ADIO_INDIVIDUAL) {
        int flag = 0;
        /* Should have already been flattened in ADIO_Open */
        num_filetypes = -1;
        while (!flag) {
            num_filetypes++;
            for (i = 0; i < flat_file_p->count; i++) {
                /* Start on a non zero-length region */
                if (flat_file_p->blocklens[i]) {
                    if (fd->disp + flat_file_p->indices[i] +
                        (num_filetypes * filetype_extent) +
                        flat_file_p->blocklens[i] > fd->fp_ind &&
                        fd->disp + flat_file_p->indices[i] <= fd->fp_ind) {
                        flat_file_index = i;
                        cur_flat_file_reg_off = fd->fp_ind -
                            (fd->disp + flat_file_p->indices[i] +
                             (num_filetypes * filetype_extent));
                        flag = 1;
                        break;
                    } else
                        bytes_into_filetype += flat_file_p->blocklens[i];
                }
            }
        }
        /* Impossible that we don't find it in this datatype */
        assert(i != flat_file_p->count);
    } else {
        num_filetypes = (int) (offset / num_etypes_in_filetype);
        etypes_in_filetype = (int) (offset % num_etypes_in_filetype);
        size_into_filetype = etypes_in_filetype * etype_size;

        tmp_filetype_size = 0;
        for (i = 0; i < flat_file_p->count; i++) {
            tmp_filetype_size += flat_file_p->blocklens[i];
            if (tmp_filetype_size > size_into_filetype) {
                flat_file_index = i;
                cur_flat_file_reg_off = flat_file_p->blocklens[i] -
                    (tmp_filetype_size - size_into_filetype);
                //bytes_into_filetype = offset * filetype_size - flat_file_p->blocklens[i];
                bytes_into_filetype = offset * filetype_size;
                break;
            }
        }
    }
#ifdef DEBUG_LIST
    fprintf(stderr, "ADIOI_MOCHIO_StridedListIO: (fd->fp_ind=%lld,fd->disp=%lld,"
            " offset=%Ld)\n(flat_file_index=%d,cur_flat_file_reg_off=%ld,"
            "bytes_into_filetype=%d)\n",
            fd->fp_ind, fd->disp, offset, flat_file_index,
            cur_flat_file_reg_off, bytes_into_filetype);
#endif
#ifdef DEBUG_LIST2
    fprintf(stderr, "flat_buf:\n");
    for (i = 0; i < flat_buf_p->count; i++)
        fprintf(stderr, "(offset, length) = (%Ld, %lld)\n",
                flat_buf_p->indices[i], flat_buf_p->blocklens[i]);
    fprintf(stderr, "flat_file:\n");
    for (i = 0; i < flat_file_p->count; i++)
        fprintf(stderr, "(offset, length) = (%Ld, %lld)\n",
                flat_file_p->indices[i], flat_file_p->blocklens[i]);
#endif

    /* total data written */
    cur_io_size = 0;
    while (cur_io_size != io_size) {
        /* Initialize the temporarily unrolling lists and
         * and associated variables */
        buf_ol_count = 0;
        file_ol_count = 0;
        for (i = 0; i < MAX_OL_COUNT; i++) {
            buf_off_arr[i] = 0;
            buf_len_arr[i] = 0;
            file_off_arr[i] = 0;
            file_len_arr[i] = 0;
        }

        /* Generate the offset-length pairs for a
         * list I/O operation */
        gen_listio_arr(flat_buf_p,
                       &flat_buf_index,
                       &cur_flat_buf_reg_off,
                       buftype_size,
                       buftype_extent,
                       flat_file_p,
                       &flat_file_index,
                       &cur_flat_file_reg_off,
                       filetype_size,
                       filetype_extent,
                       MAX_OL_COUNT,
                       fd->disp,
                       bytes_into_filetype,
                       &cur_io_size,
                       io_size,
                       buf,
                       buf_off_arr,
                       buf_len_arr, &buf_ol_count, file_off_arr, file_len_arr, &file_ol_count);

        assert(buf_ol_count <= MAX_OL_COUNT);
        assert(file_ol_count <= MAX_OL_COUNT);
#ifdef DEBUG_LIST2
        print_buf_file_ol_pairs(buf_off_arr,
                                buf_len_arr,
                                buf_ol_count,
                                file_off_arr, file_len_arr, file_ol_count, buf, rw_type);
#endif
#ifdef DEBUG_LIST2
#if 0
        do {
            int y, z;
            fprintf(stderr, "ad_pvfs2_io_list.c::\n");
            for (y = 0; y < buf_ol_count; y++) {
                for (z = 0; z < buf_len_arr[y]; z++) {
                    fprintf(stderr, "buf[%d][%d]=%c\n", y, z, ((char *) buf + buf_off_arr[y])[z]);
                }
            }
        } while (0);
#endif
#endif

        /* Run list I/O operation */
        if (rw_type == READ_OP) {
            response =
                mochio_read(fd->fs_ptr, fd->filename, buf_ol_count, (const char **) buf_off_arr,
                            buf_len_arr, file_ol_count, file_off_arr, file_len_arr);
        } else {
            response =
                mochio_write(fd->fs_ptr, fd->filename, buf_ol_count, (const char **) buf_off_arr,
                             buf_len_arr, file_ol_count, file_off_arr, file_len_arr);
        }
        if (response == -1) {
            fprintf(stderr, "ADIOI_MOCHIO_StridedListIO: Warning - mochio_"
                    "read/write completed %lld bytes.\n", (long long) response);
            *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                               MPIR_ERR_RECOVERABLE,
                                               myname, __LINE__, -1, "Error in mochio i/o\n", 0);
            goto error_state;
        }
        total_bytes_accessed += response;

    }

#ifdef DEBUG_LIST
    fprintf(stderr, "ADIOI_MOCHIO_StridedListIO: "
            "total_bytes_accessed=%Ld,ret=%lld\n", total_bytes_accessed, response);
#endif

    if (file_ptr_type == ADIO_INDIVIDUAL)
        fd->fp_ind += total_bytes_accessed;
    *error_code = MPI_SUCCESS;

  error_state:
#ifdef HAVE_STATUS_SET_BYTES
    /* TODO: why the cast? */
    MPIR_Status_set_bytes(status, datatype, total_bytes_accessed);
/* This is a temporary way of filling in status. The right way is to
   keep track of how much data was actually written by ADIOI_BUFFERED_WRITE. */
#endif

    return 0;
}
