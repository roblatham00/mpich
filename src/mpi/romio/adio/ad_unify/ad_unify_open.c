
#include "ad_unify.h"
#include "ad_unify_common.h"

#include <libgen.h>


void ADIOI_UNIFY_Open(ADIO_File fd, int *error_code)
{
    int ret, rank = 0;

    if (!fd) {
        *error_code = MPI_ERR_FILE;
        return;
    }

    if (!error_code) {
        *error_code = MPI_ERR_FILE;
        return;
    }

    /* permissions and modes are not part of Unify's api */

    MPI_Comm_rank(fd->comm, &rank);
    ADIOI_UNIFY_fs *unifyfs_blob = (ADIOI_UNIFY_fs *) ADIOI_Malloc(sizeof(ADIOI_UNIFY_fs));

    unifyfs_blob->fshdl = ADIOI_UNIFY_Init(error_code);
    unifyfs_blob->gfid = UNIFYFS_INVALID_GFID;
    if (*error_code != MPI_SUCCESS)
        return;

    /* Unify's Global FS ID could let us use the PVFS2-era "open on one,
     * broadcast to all" optimization, but there is a small bit of client-side
     * library state which gets initialized in the create/open call */
    if (rank == fd->hints->ranklist[0] && fd->fs_ptr == NULL) {
        if (fd->access_mode & ADIO_CREATE) {
            /* only one process needs to call create */
            ret = unifyfs_create(unifyfs_blob->fshdl, 0, fd->filename, &(unifyfs_blob->gfid));
            /* Unify create is always exclusive, but MPI-IO create is exclusive
             * only if flag set  */
            if (ret == EEXIST && fd->access_mode & ADIO_EXCL) {
                *error_code = MPI_ERR_FILE_EXISTS;
                return;
            }
        }
    }
    /* gfid valid only if CREAT flag set and file did not exist.
     * everyone is going to open an existing file now, except rank 0 if it
     * created a brand new file */
    if (unifyfs_blob->gfid == UNIFYFS_INVALID_GFID)
        ret = unifyfs_open(unifyfs_blob->fshdl, 0, fd->filename, &(unifyfs_blob->gfid));

    if (ret != UNIFYFS_SUCCESS) {
        ADIOI_Free(unifyfs_blob);
        fd->fs_ptr = NULL;
        *error_code = ADIOI_Err_create_code("ADIOI_Unify_Open", fd->filename, ret);
        return;
    }

    fd->fs_ptr = unifyfs_blob;
    *error_code = MPI_SUCCESS;

    return;
}
