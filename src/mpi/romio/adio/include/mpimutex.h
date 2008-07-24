/*
 * $COPYRIGHT$
 */

/* observation:  this datatype can be extended to create a single memory region
 * containing 'fp' followed by the N bytes of 'waitlist' */

typedef struct fp_data {
	MPI_Offset fp;
	unsigned char waitlist[1];
} fp_data_t;

typedef struct mpimutex {
	int nprocs, myrank, homerank;
	MPI_Comm comm;
	MPI_Win waitlistwin;
	MPI_Datatype waitlisttype, fptype;
	fp_data_t *data;
} *mpimutex_t;



int ADIOI_MPIMUTEX_Create(int homerank, MPI_Comm comm, mpimutex_t *mutex_p);
int ADIOI_MPIMUTEX_Lock(mpimutex_t mutex);
int ADIOI_MPIMUTEX_Unlock(mpimutex_t mutex);
int ADIOI_MPIMUTEX_Free(mpimutex_t *mutex_p);

int ADIOI_MPIMUTEX_Fetch_and_increment(mpimutex_t mutex,
		MPI_Offset *current, MPI_Offset increment);
int ADIOI_MPIMUTEX_Set(mpimutex_t mutex, MPI_Offset value);

int ADIOI_MPIMUTEX_Get(mpimutex_t mutex, MPI_Offset *value);
