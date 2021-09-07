/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#ifndef AD_UNIFY_H_INCLUDED
#define AD_UNIFY_H_INCLUDED

#include "adio.h"

void ADIOI_UNIFY_Open(ADIO_File fd, int *error_code);
void ADIOI_UNIFY_Close(ADIO_File fd, int *error_code);
void ADIOI_UNIFY_ReadContig(ADIO_File fd, void *buf, int count,
                            MPI_Datatype datatype, int file_ptr_type,
                            ADIO_Offset offset, ADIO_Status * status, int
                            *error_code);
void ADIOI_UNIFY_WriteContig(ADIO_File fd, const void *buf, int count,
                             MPI_Datatype datatype, int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status * status, int
                             *error_code);
void ADIOI_GEN_Fcntl(ADIO_File fd, int flag, ADIO_Fcntl_t * fcntl_struct, int
                     *error_code);
void ADIOI_GEN_WriteStrided(ADIO_File fd, const void *buf, int count,
                            MPI_Datatype datatype, int file_ptr_type,
                            ADIO_Offset offset, ADIO_Status * status, int
                            *error_code);
void ADIOI_GEN_ReadStrided(ADIO_File fd, void *buf, int count,
                           MPI_Datatype datatype, int file_ptr_type,
                           ADIO_Offset offset, ADIO_Status * status, int
                           *error_code);
void ADIOI_UNIFY_Flush(ADIO_File fd, int *error_code);
void ADIOI_GEN_Delete(const char *filename, int *error_code);
void ADIOI_GEN_Resize(ADIO_File fd, ADIO_Offset size, int *error_code);
void ADIOI_GEN_SetInfo(ADIO_File fd, MPI_Info users_info, int *error_code);
int ADIOI_GEN_Feature(ADIO_File fd, int flag);

void ADIOI_GEN_IReadContig(ADIO_File fd, void *buf, int count,
                           MPI_Datatype datatype, int file_ptr_type,
                           ADIO_Offset offset, MPI_Request * request, int *error_code);
void ADIOI_GEN_IWriteContig(ADIO_File fd, const void *buf, int count,
                            MPI_Datatype datatype, int file_ptr_type,
                            ADIO_Offset offset, MPI_Request * request, int *error_code);
void ADIOI_GEN_AIO_contig(ADIO_File fd, void *buf, int count,
                          MPI_Datatype datatype, int file_ptr_type,
                          ADIO_Offset offset, MPI_Request * request, int flag, int *error_code);
#endif /* AD_UNIFY_H_INCLUDED */
