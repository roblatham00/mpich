#include "ad_mochio.h"
#include "ad_mochio_common.h"

int ADIOI_MOCHIO_Initialized = MPI_KEYVAL_INVALID;

static void ADIOI_MOCHIO_End(int *error_code, mochio_client_t client_info)
{
    mochio_finalize(client_info);
    *error_code = MPI_SUCCESS;
}

static int ADIOI_MOCHIO_End_call(MPI_Comm comm, int keyval, void *attribute_val, void *extra_state)
{
    int error_code;
    mochio_client_t client_info = (mochio_client_t) extra_state;
    ADIOI_MOCHIO_End(&error_code, client_info);
    MPI_Keyval_free(&keyval);
    return error_code;
}

mochio_client_t ADIOI_MOCHIO_Init(MPI_Comm comm, int *error_code)
{
    mochio_client_t client_info = NULL;
    int can_skip, is_initialized = 0;
    int flag;
    /* hate to put another collective call here, but agreeing we all need to
     * initialize is cheaper than always re-attaching to the mochi ssg group.
     * you can see this pattern in ADIO_Open when we try to figure out if we
     * need to read the system-wide hints file */
    if (ADIOI_MOCHIO_Initialized != MPI_KEYVAL_INVALID)
        is_initialized = 1;

    MPI_Allreduce(&is_initialized, &can_skip, 1, MPI_INT, MPI_MIN, comm);

    if (!can_skip) {
        if (is_initialized) {
            /* need to basically un-initialize any prior mochio so we can start with a clean slate */
            MPI_Comm_delete_attr(MPI_COMM_SELF, ADIOI_MOCHIO_Initialized);
            /*
             * mochio_client_t old_info;
             * MPI_Comm_get_attr(MPI_COMM_SELF, ADIOI_MOCHIO_Initialized, &old_info, &flag);
             * if (old_info != NULL) mochio_finalize(old_info);
             */
        }
        client_info = mochio_init(comm, getenv("MOCHIO_STATEFILE"));
        MPI_Comm_create_keyval(MPI_NULL_COPY_FN, ADIOI_MOCHIO_End_call, &ADIOI_MOCHIO_Initialized,
                               (void *) client_info);
        MPI_Comm_set_attr(MPI_COMM_SELF, ADIOI_MOCHIO_Initialized, client_info);
    } else {
        /* everyone has already initialized mochio in an earlier open call */
        MPI_Comm_get_attr(MPI_COMM_SELF, ADIOI_MOCHIO_Initialized, &client_info, &flag);
    }

    *error_code = MPI_SUCCESS;
    return client_info;
}
