
#include "ad_unify.h"
#include "adio.h"
#include "ad_unify_common.h"

void ADIOI_UNIFY_Delete(const char *filename, int *error_code)
{
    int ret;
    ADIOI_UNIFY_fs *unifyfs_blob = (ADIOI_UNIFY_fs *) ADIOI_Malloc(sizeof(ADIOI_UNIFY_fs));
    unifyfs_blob->fshdl = ADIOI_UNIFY_Init(error_code);
    if (*error_code != MPI_SUCCESS)
        return;

    ret = unifyfs_remove(unifyfs_blob->fshdl, filename);
    if (ret != UNIFYFS_SUCCESS) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           __func__, __LINE__,
                                           MPI_ERR_FILE, "Error in unifyfs_remove", 0);
    } else {
        *error_code = MPI_SUCCESS;
    }
    return;
}
