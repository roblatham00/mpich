#include "ad_bv.h"
#include "ad_bv_common.h"

/* collectively called among all processes, which means we can use the "open +
 * bcast" optimization from pvfs2 */
void ADIOI_BV_Open(ADIO_File fd, int *error_code)
{
    int perm, old_mask, amode, ret;
    int rank;

    if (fd->perm == ADIO_PERM_NULL) {
        old_mask = umask(022);
        umask(old_mask);
        perm = old_mask ^ 0666;
    } else
        perm = fd->perm;

    amode = 0;
    if (fd->access_mode & ADIO_CREATE)
        amode = amode | O_CREAT;
    if (fd->access_mode & ADIO_RDONLY)
        amode = amode | O_RDONLY;
    if (fd->access_mode & ADIO_WRONLY)
        amode = amode | O_WRONLY;
    if (fd->access_mode & ADIO_RDWR)
        amode = amode | O_RDWR;
    if (fd->access_mode & ADIO_EXCL)
        amode = amode | O_EXCL;



    struct bv_stats file_stats;

    /* Note that ADIOI_BV_Init is collective (maybe I should give it an _all
     * suffix).  In some drivers we try to initialize the underlying library
     * just once.  We wrote ADIOI_BV_Init, however, so that we can call it lots
     * of times (don't do that, though).  If benvolio has already been
     * initialized, ADIOI_BV_Init won't do any work and will instead pull data
     * from the attribute we placed on COMM_SELF */
    fd->fs_ptr = ADIOI_BV_Init(fd->comm, error_code);

    MPI_Comm_rank(fd->comm, &rank);
    if (rank == fd->hints->ranklist[0] )
        ret = bv_declare(fd->fs_ptr, fd->filename, amode, perm);
    MPI_Bcast(&ret, 1, MPI_INT, fd->hints->ranklist[0], fd->comm);

    if (ret != 0) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           "bv_declare", __LINE__, MPI_ERR_FILE, "Benvolio error", 0);
        return;
    }
    /* because of some internal state, we have to have all clients call stat */
    ret = bv_stat(fd->fs_ptr, fd->filename, &file_stats);

    if (ret != 0) {
        *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                           MPIR_ERR_RECOVERABLE,
                                           "bv_stat", __LINE__, MPI_ERR_FILE, "Benvolio error", 0);
        return;
    }

    *error_code = MPI_SUCCESS;
}
