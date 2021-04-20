/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "adio.h"
#include "adio_extern.h"
#include "ad_unify.h"

#include "ad_unify_common.h"

#include <stdint.h>
#include <string.h>


#define UNIFY_READ 0
#define UNIFY_WRITE 1

static void UNIFY_IOContig(ADIO_File fd,
                           void *buf,
                           int count,
                           MPI_Datatype datatype,
                           int file_ptr_type,
                           ADIO_Offset offset, ADIO_Status * status, int io_flag, int *error_code)
{
    MPI_Count datatype_size;
    size_t mem_len;
    uint64_t file_offset = offset;
    static char myname[] = "ADIOI_UNIFY_IOCONTIG";
    unifyfs_rc ret = UNIFYFS_SUCCESS;
    ADIOI_UNIFY_fs *unifyfs_blob = fd->fs_ptr;

    if (count == 0) {
        ret = 0;
        goto fn_exit;
    }

    MPI_Type_size_x(datatype, &datatype_size);
    mem_len = datatype_size * count;

    if (file_ptr_type == ADIO_INDIVIDUAL)
        file_offset = fd->fp_ind;


    unifyfs_io_request req = {
        .user_buf = buf,
        .nbytes = mem_len,
        .offset = offset,
        .gfid = unifyfs_blob->gfid,
    };
    switch (io_flag) {
        case UNIFY_READ:
            req.op = UNIFYFS_IOREQ_OP_READ;
            break;
        case UNIFY_WRITE:
            req.op = UNIFYFS_IOREQ_OP_WRITE;
            break;
        default:
            /* 'unknown flag' handeled in unifyfs_dispatch_io */
            break;
    };

    ret = unifyfs_dispatch_io(unifyfs_blob->fshdl, 1, &req);

    if (ret != UNIFYFS_SUCCESS) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           myname, __LINE__, MPI_ERR_IO, "Error in Dispatch:", "%s",
                                           strerror(ret));
        goto exit;
    }
    ret = unifyfs_wait_io(unifyfs_blob->fshdl, 1, &req, 1);
    if (ret != UNIFYFS_SUCCESS) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           myname, __LINE__, MPI_ERR_IO, "Error in I/O Wait:", "%s",
                                           strerror(ret));
        goto exit;
    }
    if (req.result.error != UNIFYFS_SUCCESS) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           myname, __LINE__, MPI_ERR_IO, "Error in I/O operation:",
                                           "%s", strerror(req.result.error));
        goto exit;
    }


    if (file_ptr_type == ADIO_INDIVIDUAL)
        fd->fp_ind += ret;
    fd->fp_sys_posn = file_offset + ret;

  fn_exit:
#ifdef HAVE_STATUS_SET_BYTES
    if (status)
        MPIR_Status_set_bytes(status, datatype, ret);
#endif

    *error_code = MPI_SUCCESS;

  exit:
    return;
}

void ADIOI_UNIFY_ReadContig(ADIO_File fd,
                            void *buf,
                            int count,
                            MPI_Datatype datatype,
                            int file_ptr_type,
                            ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    UNIFY_IOContig(fd, buf, count, datatype, file_ptr_type, offset, status, UNIFY_READ, error_code);
}

void ADIOI_UNIFY_WriteContig(ADIO_File fd,
                             const void *buf,
                             int count,
                             MPI_Datatype datatype,
                             int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    UNIFY_IOContig(fd,
                   (void *) buf,
                   count, datatype, file_ptr_type, offset, status, UNIFY_WRITE, error_code);
}
