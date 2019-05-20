#include "ad_mochio.h"

#include <mochio.h>
/* open and close are no-ops as romio-mochi operates on file names only
 */
void ADIOI_MOCHIO_Close(ADIO_File fd, int *error_code)
{
    int rank;
    mochio_client_t client_info = fd->fs_ptr;
    MPI_Comm_rank(fd->comm, &rank);
    if (rank == 0 && getenv("MOCHIO_SHOW_STATS"))
        mochio_statistics(client_info);
    *error_code = MPI_SUCCESS;
}
