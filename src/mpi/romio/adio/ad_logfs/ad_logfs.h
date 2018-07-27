/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs.h,v 1.2 2002/10/24 17:01:03 gropp Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#ifndef AD_LOGFS_INCLUDE
#define AD_LOGFS_INCLUDE

#include <unistd.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <fcntl.h>
#include "adio.h"

void ADIOI_LOGFS_Open(ADIO_File fd, int *error_code);
void ADIOI_LOGFS_Close(ADIO_File fd, int *error_code);
void ADIOI_LOGFS_ReadContig(ADIO_File fd, void *buf, int count,
                            MPI_Datatype datatype, int file_ptr_type,
                            ADIO_Offset offset, ADIO_Status * status, int
                            *error_code);
void ADIOI_LOGFS_WriteContig(ADIO_File fd, const void *buf, int count,
                             MPI_Datatype datatype, int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status * status, int
                             *error_code);
void ADIOI_LOGFS_IwriteContig(ADIO_File fd, const void *buf, int count,
                              MPI_Datatype datatype, int file_ptr_type,
                              ADIO_Offset offset, ADIO_Request * request, int
                              *error_code);
void ADIOI_LOGFS_IreadContig(ADIO_File fd, void *buf, int count,
                             MPI_Datatype datatype, int file_ptr_type,
                             ADIO_Offset offset, ADIO_Request * request, int
                             *error_code);
int ADIOI_LOGFS_ReadDone(ADIO_Request * request, ADIO_Status * status, int
                         *error_code);
int ADIOI_LOGFS_WriteDone(ADIO_Request * request, ADIO_Status * status, int
                          *error_code);
void ADIOI_LOGFS_ReadComplete(ADIO_Request * request, ADIO_Status * status, int
                              *error_code);
void ADIOI_LOGFS_WriteComplete(ADIO_Request * request, ADIO_Status * status, int *error_code);
void ADIOI_LOGFS_Fcntl(ADIO_File fd, int flag, ADIO_Fcntl_t * fcntl_struct, int *error_code);
void ADIOI_LOGFS_WriteStrided(ADIO_File fd, const void *buf, int count,
                              MPI_Datatype datatype, int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status * status, int *error_code);
void ADIOI_LOGFS_ReadStrided(ADIO_File fd, void *buf, int count,
                             MPI_Datatype datatype, int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status * status, int
                             *error_code);
void ADIOI_LOGFS_WriteStridedColl(ADIO_File fd, const void *buf, int count,
                                  MPI_Datatype datatype, int file_ptr_type,
                                  ADIO_Offset offset, ADIO_Status * status, int
                                  *error_code);
void ADIOI_LOGFS_ReadStridedColl(ADIO_File fd, void *buf, int count,
                                 MPI_Datatype datatype, int file_ptr_type,
                                 ADIO_Offset offset, ADIO_Status * status, int
                                 *error_code);
void ADIOI_LOGFS_IreadStrided(ADIO_File fd, void *buf, int count,
                              MPI_Datatype datatype, int file_ptr_type,
                              ADIO_Offset offset, ADIO_Request * request, int
                              *error_code);
void ADIOI_LOGFS_IwriteStrided(ADIO_File fd, const void *buf, int count,
                               MPI_Datatype datatype, int file_ptr_type,
                               ADIO_Offset offset, ADIO_Request * request, int
                               *error_code);
void ADIOI_LOGFS_Flush(ADIO_File fd, int *error_code);
void ADIOI_LOGFS_Resize(ADIO_File fd, ADIO_Offset size, int *error_code);
ADIO_Offset ADIOI_LOGFS_SeekIndividual(ADIO_File fd, ADIO_Offset offset,
                                       int whence, int *error_code);
void ADIOI_LOGFS_SetInfo(ADIO_File fd, MPI_Info users_info, int *error_code);
void ADIOI_LOGFS_Get_shared_fp(ADIO_File fd, int size, ADIO_Offset * shared_fp, int *error_code);
void ADIOI_LOGFS_Set_shared_fp(ADIO_File fd, ADIO_Offset offset, int *error_code);
void ADIOI_LOGFS_Delete(const char *filename, int *error_code);


void ADIOI_LOGFS_Set_slave(ADIO_File fd, ADIOI_Fns * slaveops);

int ADIOI_LOGFS_Feature(ADIO_File fd, int flag);


#endif
