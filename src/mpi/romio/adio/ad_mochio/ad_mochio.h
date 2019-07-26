#include "adio.h"

void ADIOI_MOCHIO_Open(ADIO_File fd, int *error_code);

void ADIOI_MOCHIO_Close(ADIO_File fd, int *error_code);

void ADIOI_MOCHIO_ReadContig(ADIO_File fd,
                             void *buf,
                             int count,
                             MPI_Datatype datatype,
                             int file_ptr_type,
                             ADIO_Offset offset, ADIO_Status * status, int *error_code);

void ADIOI_MOCHIO_WriteContig(ADIO_File fd,
                              const void *buf,
                              int count,
                              MPI_Datatype datatype,
                              int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status * status, int *error_code);
void ADIOI_MOCHIO_ReadStrided(ADIO_File fd,
                              void *buf,
                              int count,
                              MPI_Datatype datatype,
                              int file_ptr_type,
                              ADIO_Offset offset, ADIO_Status * status, int *error_code);

void ADIOI_MOCHIO_WriteStrided(ADIO_File fd,
                               const void *buf,
                               int count,
                               MPI_Datatype datatype,
                               int file_ptr_type,
                               ADIO_Offset offset, ADIO_Status * status, int *error_code);

int ADIOI_MOCHIO_StridedListIO(ADIO_File fd, void *buf, int count,
                               MPI_Datatype datatype, int file_ptr_type,
                               ADIO_Offset offset, ADIO_Status * status,
                               int *error_code, int rw_type);

int ADIOI_MOCHIO_Feature(ADIO_File fd, int flag);
