#include "adio.h"

void ADIOI_BV_Open(ADIO_File fd, int *error_code);

void ADIOI_BV_Close(ADIO_File fd, int *error_code);

void ADIOI_BV_ReadContig(ADIO_File fd,
                             void *buf,
                             int count,
                             MPI_Datatype datatype,
                             int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status * status, int *error_code);

void ADIOI_BV_WriteContig(ADIO_File fd,
                              const void *buf,
                              int count,
                              MPI_Datatype datatype,
                              int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status * status, int *error_code);
void ADIOI_BV_ReadStrided(ADIO_File fd,
                              void *buf,
                              int count,
                              MPI_Datatype datatype,
                              int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status * status, int *error_code);

void ADIOI_BV_WriteStrided(ADIO_File fd,
                               const void *buf,
                               int count,
                               MPI_Datatype datatype,
                               int file_ptr_type,
                               ADIO_Offset offset, ADIO_Status * status, int *error_code);

int ADIOI_BV_StridedListIO(ADIO_File fd, void *buf, int count,
                               MPI_Datatype datatype, int file_ptr_type,
                               ADIO_Offset offset, ADIO_Status * status,
                               int *error_code, int rw_type);
void ADIOI_BV_OldStridedListIO(ADIO_File fd, void *buf, int count,
                                   MPI_Datatype datatype, int file_ptr_type,
                                   ADIO_Offset offset, ADIO_Status * status,
                                   int *error_code, int rw_type);

void ADIOI_BV_Fcntl(ADIO_File fd, int flag, ADIO_Fcntl_t * fcntl_struct, int *error_code);

int ADIOI_BV_Feature(ADIO_File fd, int flag);
