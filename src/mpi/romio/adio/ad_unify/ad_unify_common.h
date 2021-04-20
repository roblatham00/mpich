
#ifndef AD_UNIFY_COMMON_H
#define AD_UNIFY_COMMON_H

#include <unifyfs/unifyfs_api.h>


extern struct ADIOI_UNIFY_fs_s {
    unifyfs_handle fshdl;
    unifyfs_gfid gfid;
} ADIOI_UNIFY_fs_s;

typedef struct ADIOI_UNIFY_fs_s ADIOI_UNIFY_fs;

unifyfs_handle ADIOI_UNIFY_Init(int *error_code);
void ADIOI_UNIFY_End(unifyfs_handle fshdl, int *error_code);
int ADIOI_UNIFY_End_call(MPI_Comm comm, int keyval, void *attribute_val, void *extra_state);

#endif
