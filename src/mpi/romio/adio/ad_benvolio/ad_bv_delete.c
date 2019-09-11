
#include "ad_bv.h"
#include "adio.h"

#include "ad_bv_common.h"

#include <bv.h>

void ADIOI_BV_Delete(const char *filename, int *error_code)
{
    int ret;
    bv_client_t client = ADIOI_BV_Init(MPI_COMM_SELF, error_code);

    ret = bv_delete(client, filename);
    /* --BEGIN ERROR HANDLING-- */
    if (ret == -1) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           "ADIOI_BV_Delete", __LINE__, MPI_ERR_IO,
                                           "Error in bv_delete", " ", ret);
    }
    /* --END ERROR HANDLING-- */
    else
        *error_code = MPI_SUCCESS;


    return;
}
