#include <mochio.h>

enum {
    READ_OP,
    WRITE_OP
};
mochio_client_t ADIOI_MOCHIO_Init(MPI_Comm comm, int *error_code);
