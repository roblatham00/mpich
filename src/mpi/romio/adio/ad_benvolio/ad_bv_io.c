/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*-
 *
 *
 *   Copyright (C) 1997 University of Chicago.
 *   Copyright (C) 2017 DataDirect Networks.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "adio.h"
#include "adio_extern.h"
#include "ad_bv.h"
#include "ad_bv_common.h"

#include <stdint.h>

#include <bv.h>

#define BV_READ 0
#define BV_WRITE 1

static void BV_IOContig(ADIO_File fd,
                            void *buf,
                            int count,
                            MPI_Datatype datatype,
                            int file_ptr_type,
                            ADIO_Offset offset, ADIO_Status * status, int io_flag, int *error_code)
{
    ssize_t ret;
    MPI_Count datatype_size;
    size_t mem_len;
    off_t file_offset = offset;
    static char myname[] = "ADIOI_BV_IOCONTIG";
    const char *mem_addr;
    uint64_t file_size;

    MPI_Type_size_x(datatype, &datatype_size);
    mem_len = datatype_size * count;

    if (file_ptr_type == ADIO_INDIVIDUAL)
        file_offset = fd->fp_ind;

    mem_addr = buf;
    file_size = mem_len;

    switch (io_flag) {
        case BV_READ:
            ret =
                bv_read(fd->fs_ptr, fd->filename, 1, &mem_addr, &mem_len, 1, &file_offset,
                            &file_size);
            break;
        case BV_WRITE:
            ret =
                bv_write(fd->fs_ptr, fd->filename, 1, &mem_addr, &mem_len, 1, &file_offset,
                             &file_size);
            break;
        default:
            *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                               MPIR_ERR_RECOVERABLE,
                                               myname, __LINE__, MPI_ERR_IO, "Unknown flag", 0);
            goto exit;

            break;
    };

    /* Let the application decide how to fail */
    if (ret < 0) {
        *error_code = MPI_SUCCESS;
        goto exit;
    }

    if (file_ptr_type == ADIO_INDIVIDUAL)
        fd->fp_ind += ret;
    fd->fp_sys_posn = file_offset + ret;

#ifdef HAVE_STATUS_SET_BYTES
    MPIR_Status_set_bytes(status, datatype, ret);
#endif

    *error_code = MPI_SUCCESS;

  exit:
    return;
}

void ADIOI_BV_ReadContig(ADIO_File fd,
                             void *buf,
                             int count,
                             MPI_Datatype datatype,
                             int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    BV_IOContig(fd, buf, count, datatype, file_ptr_type, offset, status, BV_READ,
                    error_code);
}

void ADIOI_BV_WriteContig(ADIO_File fd,
                              const void *buf,
                              int count,
                              MPI_Datatype datatype,
                              int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    BV_IOContig(fd,
                    (void *) buf,
                    count, datatype, file_ptr_type, offset, status, BV_WRITE, error_code);
}

void ADIOI_BV_ReadStrided(ADIO_File fd,
                              void *buf,
                              int count,
                              MPI_Datatype datatype,
                              int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    ADIOI_BV_OldStridedListIO(fd, buf, count, datatype, file_ptr_type, offset, status,
                                  error_code, READ_OP);
}

void ADIOI_BV_WriteStrided(ADIO_File fd,
                               const void *buf,
                               int count,
                               MPI_Datatype datatype,
                               int file_ptr_type,
                               ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    ADIOI_BV_OldStridedListIO(fd, (void *) buf, count, datatype, file_ptr_type, offset, status,
                                  error_code, WRITE_OP);
}
