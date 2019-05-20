#include <bv.h>

enum {
    READ_OP,
    WRITE_OP
};
bv_client_t ADIOI_BV_Init(MPI_Comm comm, int *error_code);
