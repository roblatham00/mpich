
#include "ad_unify.h"
#include "ad_unify_common.h"

#include <libgen.h>


void ADIOI_UNIFY_Open(ADIO_File fd, int *error_code)
{
    static char myname[] = "ADIO_UNIFY_OPEN";
    int perm, amode = 0, ret, rank = 0;
    mode_t old_mask;

    if (!fd) {
        *error_code = MPI_ERR_FILE;
        return;
    }

    if (!error_code) {
        *error_code = MPI_ERR_FILE;
        return;
    }

    if (fd->perm == ADIO_PERM_NULL) {
        old_mask = umask(022);
        umask(old_mask);
        perm = old_mask ^ 0666;
    } else
        perm = fd->perm;

    /* setup the file access mode */
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

    MPI_Comm_rank(fd->comm, &rank);
    ADIOI_UNIFY_fs *unifyfs_blob = (ADIOI_UNIFY_fs *) ADIOI_Malloc(sizeof(ADIOI_UNIFY_fs));

    unifyfs_blob->fshdl = ADIOI_UNIFY_Init(error_code);
    if (*error_code != MPI_SUCCESS)
        return;

    /* Unify's Global FS ID allows us to use the PVFS2-era "open on one,
     * broadcast to all" optimization */
    if (rank == fd->hints->ranklist[0] && fd->fs_ptr == NULL) {
        if (amode & O_CREAT) {
            ret = unifyfs_create(unifyfs_blob->fshdl, amode, fd->filename, &(unifyfs_blob->gfid));
        } else {
            ret = unifyfs_open(unifyfs_blob->fshdl, amode, fd->filename, &(unifyfs_blob->gfid));
        }

        /* OK to read from laminated file but cannot write to it.  We'll assume
         * the caller was just being lazy -- few if any users specify precisely
         * the open flags they requre. most open all files RDWR even if they
         * will only read */
        if (ret == EROFS && (amode & O_RDWR)) {
            amode ^= O_RDWR;
            amode |= O_RDONLY;
        }
        if (amode & O_CREAT) {
            ret = unifyfs_create(unifyfs_blob->fshdl, amode, fd->filename, &(unifyfs_blob->gfid));
        } else {
            ret = unifyfs_open(unifyfs_blob->fshdl, amode, fd->filename, &(unifyfs_blob->gfid));
        }
    }
    MPI_Bcast(&ret, 1, MPI_INT, fd->hints->ranklist[0], fd->comm);
    if (ret != UNIFYFS_SUCCESS) {
        ADIOI_Free(unifyfs_blob);
        fd->fs_ptr = NULL;
        *error_code = ADIOI_Err_create_code("ADIOI_Unify_Open", fd->filename, ret);
        return;
    }
    MPI_Bcast(&(unifyfs_blob->gfid), 1, MPI_INT32_T, fd->hints->ranklist[0], fd->comm);

    fd->fs_ptr = unifyfs_blob;
    *error_code = MPI_SUCCESS;

    return;
}
