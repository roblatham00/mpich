#include "ad_bv.h"
#include "ad_bv_common.h"

/* collectively called among all processes, which means we can collectively execute romio_init() */
void ADIOI_BV_Open(ADIO_File fd, int *error_code)
{
    int perm, old_mask, amode;

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
    fd->fs_ptr = ADIOI_BV_Init(fd->comm, error_code);
    bv_declare(fd->fs_ptr, fd->filename, amode, perm);
    bv_stat(fd->fs_ptr, fd->filename, &file_stats);
    *error_code = MPI_SUCCESS;
}
