#include "ad_bv.h"
#include "ad_bv_common.h"

/* collectively called among all processes, which means we can collectively execute romio_init() */
void ADIOI_BV_Open(ADIO_File fd, int *error_code)
{
    fd->fs_ptr = ADIOI_BV_Init(fd->comm, error_code);
    *error_code = MPI_SUCCESS;
}
