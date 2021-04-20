
#include "ad_unify.h"
#include "ad_unify_common.h"

#include <unifyfs/unifyfs_api.h>
int ADIOI_UNIFY_Initialized = MPI_KEYVAL_INVALID;

void ADIOI_UNIFY_End(unifyfs_handle fshdl, int *error_code)
{
    int ret;
    static char myname[] = "ADIOI_UNIFY_END";
    ret = unifyfs_finalize(fshdl);

    if (ret != 0) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           myname, __LINE__,
                                           MPI_ERR_FILE, "Error in unifyfs_finalize", "%s",
                                           strerror(ret));
        return;
    }
    *error_code = MPI_SUCCESS;
}

int ADIOI_UNIFY_End_call(MPI_Comm com, int keyval, void *attribute_val, void *extra_state)
{
    int error_code;
    MPI_Comm_free_keyval(&keyval);
    unifyfs_handle fshdl = (unifyfs_handle) extra_state;
    ADIOI_UNIFY_End(fshdl, &error_code);
    return error_code;
}

unifyfs_handle ADIOI_UNIFY_Init(int *error_code)
{
    unifyfs_handle fshdl = NULL;
    const char *env, *path;
    int flag;
    unifyfs_rc ret;

    /* Using MPI attribute to cache unify state object */
    if (ADIOI_UNIFY_Initialized != MPI_KEYVAL_INVALID) {
        MPI_Comm_get_attr(MPI_COMM_SELF, ADIOI_UNIFY_Initialized, &fshdl, &flag);
        *error_code = MPI_SUCCESS;
        return fshdl;
    }

    env = getenv("UNIFYFS_MOUNTPOINT");
    if (env == NULL)
        path = "/tmp";
    else
        path = env;

    /* here would be a good place to do some tuning, but we might not have
     * access to users hints -- the init calls happen pretty early in open and
     * delete
     *
     * One option we do need to crank down a bit is the size of each log
     * extent.  ROMIO can easily generate 4 byte requests (or smaller) but each
     * request becomes a 4 MiB unify extent by default.
     *
     * What are the performance implications of smaller chunk sizes? */
    unifyfs_cfg_option configs[] = { {.opt_name = "logio.chunk_size",
                                      .opt_value = "32768"}
    ,
    };

    ret = unifyfs_initialize(path, configs, sizeof(configs) / sizeof(configs[0]), &fshdl);
    if (ret != UNIFYFS_SUCCESS) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS, MPIR_ERR_RECOVERABLE,
                                           "Unify_init", __LINE__, MPI_ERR_IO,
                                           "Error in unify_initialize", "%s", strerror(ret));
        return UNIFYFS_INVALID_HANDLE;

    }

    MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, ADIOI_UNIFY_End_call, &ADIOI_UNIFY_Initialized,
                           fshdl);

    MPI_Comm_set_attr(MPI_COMM_SELF, ADIOI_UNIFY_Initialized, fshdl);

    *error_code = MPI_SUCCESS;
    return fshdl;
}
