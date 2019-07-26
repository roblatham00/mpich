/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*-
 *
 *
 *   Copyright (C) 1997 University of Chicago.
 *   Copyright (C) 2017 DataDirect Networks.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "adio.h"
#include "adio_extern.h"
#include "ad_mochio.h"
#include "ad_mochio_common.h"

#include <stdint.h>

#include <mochio.h>

#define MOCHIO_READ 0
#define MOCHIO_WRITE 1

static void MOCHIO_IOContig(ADIO_File fd,
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
    static char myname[] = "ADIOI_MOCHIO_IOCONTIG";
    const char *mem_addr;
    uint64_t file_size;

    MPI_Type_size_x(datatype, &datatype_size);
    mem_len = datatype_size * count;

    if (file_ptr_type == ADIO_INDIVIDUAL)
        file_offset = fd->fp_ind;

    mem_addr = buf;
    file_size = mem_len;

    switch (io_flag) {
        case MOCHIO_READ:
            ret =
                mochio_read(fd->fs_ptr, fd->filename, 1, &mem_addr, &mem_len, 1, &file_offset,
                            &file_size);
            break;
        case MOCHIO_WRITE:
            ret =
                mochio_write(fd->fs_ptr, fd->filename, 1, &mem_addr, &mem_len, 1, &file_offset,
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

void ADIOI_MOCHIO_ReadContig(ADIO_File fd,
                             void *buf,
                             int count,
                             MPI_Datatype datatype,
                             int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    MOCHIO_IOContig(fd, buf, count, datatype, file_ptr_type, offset, status, MOCHIO_READ,
                    error_code);
}

void ADIOI_MOCHIO_WriteContig(ADIO_File fd,
                              const void *buf,
                              int count,
                              MPI_Datatype datatype,
                              int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    MOCHIO_IOContig(fd,
                    (void *) buf,
                    count, datatype, file_ptr_type, offset, status, MOCHIO_WRITE, error_code);
}

void ADIOI_MOCHIO_ReadStrided(ADIO_File fd,
                              void *buf,
                              int count,
                              MPI_Datatype datatype,
                              int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    ADIOI_MOCHIO_StridedListIO(fd, buf, count, datatype, file_ptr_type, offset, status, error_code,
                               READ_OP);
}

void ADIOI_MOCHIO_WriteStrided(ADIO_File fd,
                               const void *buf,
                               int count,
                               MPI_Datatype datatype,
                               int file_ptr_type,
                               ADIO_Offset offset, ADIO_Status * status, int *error_code)
{
    ADIOI_MOCHIO_StridedListIO(fd, (void *) buf, count, datatype, file_ptr_type, offset, status,
                               error_code, WRITE_OP);
}
