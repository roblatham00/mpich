#include "ad_bv.h"
#include "ad_bv_common.h"

int ADIOI_BV_Initialized = MPI_KEYVAL_INVALID;

static void ADIOI_BV_End(int *error_code, bv_client_t client_info)
{
    bv_finalize(client_info);
    *error_code = MPI_SUCCESS;
}

static int ADIOI_BV_End_call(MPI_Comm comm, int keyval, void *attribute_val, void *extra_state)
{
    int error_code;
    bv_client_t client_info = (bv_client_t) extra_state;
    ADIOI_BV_End(&error_code, client_info);
    MPI_Keyval_free(&keyval);
    return error_code;
}

bv_client_t ADIOI_BV_Init(MPI_Comm comm, int *error_code)
{
    bv_client_t client_info = NULL;
    bv_config_t cfg = NULL;
    ssize_t cfg_size = 0;
    int can_skip, is_initialized = 0;
    int flag;
    /* hate to put another collective call here, but agreeing we all need to
     * initialize is cheaper than always re-attaching to the mochi ssg group.
     * you can see this pattern in ADIO_Open when we try to figure out if we
     * need to read the system-wide hints file */
    if (ADIOI_BV_Initialized != MPI_KEYVAL_INVALID)
        is_initialized = 1;

    MPI_Allreduce(&is_initialized, &can_skip, 1, MPI_INT, MPI_MIN, comm);

    if (!can_skip) {
        if (is_initialized) {
            /* need to basically un-initialize any prior bv so we can start with a clean slate */
            MPI_Comm_delete_attr(MPI_COMM_SELF, ADIOI_BV_Initialized);
            /*
             * bv_client_t old_info;
             * MPI_Comm_get_attr(MPI_COMM_SELF, ADIOI_BV_Initialized, &old_info, &flag);
             * if (old_info != NULL) bv_finalize(old_info);
             */
        }
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) {
            cfg = bvutil_cfg_get(getenv("BV_STATEFILE"));
            cfg_size = bvutil_cfg_getsize(cfg);
        }
        MPI_Bcast(&cfg_size, 8, MPI_BYTE, 0, comm);
        if (cfg == NULL) cfg = ADIOI_Malloc(cfg_size);
        MPI_Bcast(cfg, cfg_size, MPI_BYTE, 0, comm);

        client_info = bv_init(cfg);
        if (rank == 0)
            bvutil_cfg_free(cfg);
        else
            ADIOI_Free(cfg);
        MPI_Comm_create_keyval(MPI_NULL_COPY_FN, ADIOI_BV_End_call, &ADIOI_BV_Initialized,
                               (void *) client_info);
        MPI_Comm_set_attr(MPI_COMM_SELF, ADIOI_BV_Initialized, client_info);
    } else {
        /* everyone has already initialized bv in an earlier open call */
        MPI_Comm_get_attr(MPI_COMM_SELF, ADIOI_BV_Initialized, &client_info, &flag);
    }

    *error_code = MPI_SUCCESS;
    return client_info;
}
