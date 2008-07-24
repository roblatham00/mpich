/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 * $COPYRIGHT$
 */

/* hybrid rma one-sided with point-to-point approach as outlined in
 * "Implementing MPI-IO Atomic Mode Without File System Support". Robert Ross,
 * Robert Latham, William Gropp, Rajeev Thakur, Brian Toonen; CCGrid 2005 and
 * "Implementing MPI-IO Shared File Pointers without File System Support"
 * Robert Latham, Robert Ross, Rajeev Thakur, Brian Toonen; EuroPVM/MPI 2005
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include "mpimutex.h"

#define WAKEUP 13

/* the create, lock, unlock, free, set and get routines are the same for both
 * Atomic mode and shared file ponter.  For atomic mode we will ignore the file
 * offset.  */
int ADIOI_MPIMUTEX_Create(int rank, MPI_Comm comm, mpimutex_t *mutex_p)
{
    int nprocs, myrank, mpi_err;
    mpimutex_t mutex = NULL;

    int blklens[2];
    MPI_Aint disps[2];

    mutex = malloc(sizeof(struct mpimutex));
    if (!mutex) goto err_return;

    memset(mutex, 0, sizeof(struct mpimutex));

    mutex->homerank = rank;
    mutex->waitlistwin = MPI_WIN_NULL;
    mutex->waitlisttype = MPI_DATATYPE_NULL;
    mutex->comm = MPI_COMM_NULL;

    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);
    mutex->myrank = myrank;
    mutex->nprocs = nprocs;

    blklens[0] = mutex->myrank;
    disps[0]   = 0;
    blklens[1] = mutex->nprocs - mutex->myrank - 1;
    disps[1]   = mutex->myrank + 1;

    mpi_err = MPI_Type_indexed(2, blklens, disps, MPI_BYTE,
			       &mutex->waitlisttype);
    if (mpi_err != MPI_SUCCESS) goto err_return;
    mpi_err = MPI_Type_commit(&mutex->waitlisttype);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    /* no pre-defined type for MPI_Offset. what a pain */
    mpi_err = MPI_Type_contiguous(sizeof(MPI_Offset), MPI_BYTE,
		    &mutex->fptype);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    mpi_err = MPI_Type_commit(&mutex->fptype);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    MPI_Comm_dup(comm, &mutex->comm);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    if (myrank == rank) {
	/* allocate file pointer and waitlist in one contiguous chunk. waitlist
	 * uses the "nprocs-1" additional bytes.  why '- 1'? because
	 * 'sizeof(struct fp_data) is one MPI_Offset + one single-element
	 * waitlist array */
	mpi_err = MPI_Alloc_mem(sizeof(struct fp_data) + (nprocs-1),
			MPI_INFO_NULL, &mutex->data);
	if (mpi_err != MPI_SUCCESS) goto err_return;
	memset(mutex->data, 0, sizeof(struct fp_data) + (nprocs -1));

	mpi_err = MPI_Win_create(mutex->data, sizeof(struct fp_data)+(nprocs -1),
			1, MPI_INFO_NULL, mutex->comm, &mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

    }
    else {
	mpi_err = MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL,
				 mutex->comm, &mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;
    }

    *mutex_p = mutex;
    return MPI_SUCCESS;

 err_return:
    printf("error!\n");
    if (mutex) {
	if (mutex->waitlistwin != MPI_WIN_NULL) {
	    MPI_Win_free(&mutex->waitlistwin);
	}

	if (mutex->data) {
	    MPI_Free_mem(mutex->data);
	}

	if (mutex->waitlisttype != MPI_DATATYPE_NULL) {
	    MPI_Type_free(&mutex->waitlisttype);
	}

	if (mutex->comm != MPI_COMM_NULL) {
	    MPI_Comm_free(&mutex->comm);
	}

	free(mutex);
    }
    return MPI_ERR_UNKNOWN;
}

/* obtain lock by writing a '1' into desginated array element.  Read everyone
 * elses (with that datatype created at init time).  If no contention (a '0' in
 * all others), we have the lock.  If contention, ANY_SRC recieve and whomever
 * holds the lock will tell us when
 * they are done
 * note that shared file pointer access done in fetch_and_increment below.
 * Lock and Unlock do not touch that value */

int ADIOI_MPIMUTEX_Lock(mpimutex_t mutex)
{
    int mpi_err, i;
    unsigned char val = 1;
    unsigned char *waitlistcopy = NULL;

    waitlistcopy = malloc(mutex->nprocs-1);
    if (!waitlistcopy) goto err_return;

    /* add self to waitlist */
    mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, mutex->homerank, 0,
			   mutex->waitlistwin);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    mpi_err = MPI_Get(waitlistcopy,
	    mutex->nprocs-1, MPI_BYTE, mutex->homerank,
	    sizeof(MPI_Offset), 1, mutex->waitlisttype,
	    mutex->waitlistwin);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    mpi_err = MPI_Put(&val,
	    1, MPI_BYTE, mutex->homerank,
	    sizeof(MPI_Offset) + mutex->myrank, 1, MPI_BYTE,
	    mutex->waitlistwin);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    mpi_err = MPI_Win_unlock(mutex->homerank, mutex->waitlistwin);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    /* check to see if lock is already held */
    for (i=0; i < (mutex->nprocs - 1) && waitlistcopy[i] == 0; i++);

    if (i < mutex->nprocs - 1) {
	/* wait for notification from some other process */
	mpi_err = MPI_Recv(NULL, 0, MPI_BYTE, MPI_ANY_SOURCE, WAKEUP,
			   mutex->comm, MPI_STATUS_IGNORE);
	if (mpi_err != MPI_SUCCESS) goto err_return;
    }

    free(waitlistcopy);

    return MPI_SUCCESS;

 err_return:
    printf("error!\n");
    if (waitlistcopy) free(waitlistcopy);

    return MPI_ERR_UNKNOWN;
}

/* similar to Lock:  read all bytes of waitlist except ours.  Put a 0 in ours.
 * If any procs are waiting for a lock, send a wakeup message to that waiting
 * process */
int ADIOI_MPIMUTEX_Unlock(mpimutex_t mutex)
{
    int mpi_err, i;
    unsigned char val = 0;
    unsigned char *waitlistcopy;

    /* TODO: MOVE INTO INITIALIZE */
    waitlistcopy = malloc(mutex->nprocs-1);
    if (!waitlistcopy) goto err_return;

    /* remove self from waitlist */
    mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, mutex->homerank, 0,
			   mutex->waitlistwin);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    mpi_err = MPI_Get(waitlistcopy, mutex->nprocs-1, MPI_BYTE,
		      mutex->homerank, sizeof(MPI_Offset), 1, mutex->waitlisttype,
		      mutex->waitlistwin);
    if (mpi_err != MPI_SUCCESS) goto err_return;
    mpi_err = MPI_Put(&val, 1, MPI_BYTE,
		      mutex->homerank, sizeof(MPI_Offset) + mutex->myrank, 1, MPI_BYTE,
		      mutex->waitlistwin);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    mpi_err = MPI_Win_unlock(mutex->homerank, mutex->waitlistwin);
    if (mpi_err != MPI_SUCCESS) goto err_return;

    /* check to see if lock is already held */
    for (i=0; i < (mutex->nprocs - 1) && waitlistcopy[i] == 0; i++);

    if (i < mutex->nprocs - 1) {
	int nextrank;

	/* find the next rank waiting for the lock.  we start with the
	 * rank after ours and look in order to ensure fairness.
	 */
	nextrank = mutex->myrank;
	while (nextrank < (mutex->nprocs - 1) && waitlistcopy[nextrank] == 0)
	{
	    nextrank++;
	}
	if (nextrank < mutex->nprocs - 1) {
	    /* we got a valid rank, but "nextrank" is off by one */
	    nextrank++;
	}
	else {
	    nextrank = 0;
	    while (nextrank < mutex->myrank && waitlistcopy[nextrank] == 0)
	    {
		nextrank++;
	    }

	    assert(nextrank != mutex->myrank);
	}


	/* notify next rank */
	mpi_err = MPI_Send(NULL, 0, MPI_BYTE, nextrank, WAKEUP,
			   mutex->comm);
	if (mpi_err != MPI_SUCCESS) goto err_return;
    }

    free(waitlistcopy);

    return MPI_SUCCESS;

 err_return:
    printf("error!\n");
    if (waitlistcopy) free(waitlistcopy);
    return MPI_ERR_UNKNOWN;
}

int ADIOI_MPIMUTEX_Free(mpimutex_t *mutex_p)
{
    mpimutex_t mutex = *mutex_p;

    if (mutex) {
	if (mutex->waitlistwin != MPI_WIN_NULL) {
	    MPI_Win_free(&mutex->waitlistwin);
	}

	if (mutex->data) {
	    MPI_Free_mem(mutex->data);
	}

	if (mutex->fptype != MPI_DATATYPE_NULL) {
	    MPI_Type_free(&mutex->fptype);
	}

	if (mutex->waitlisttype != MPI_DATATYPE_NULL) {
	    MPI_Type_free(&mutex->waitlisttype);
	}

	if (mutex->comm != MPI_COMM_NULL) {
	    MPI_Comm_free(&mutex->comm);
	}

	free(mutex);
    }

    *mutex_p = NULL;

    return MPI_SUCCESS;
}

/* increment the shared resource by 'increment'.  return present value (before
 * incrementing) in 'current'.  we do end up copying a lot of code from
 * lock/unlock, but we need to get the shared fp value in between the
 * MPI_Win_Lock and MPI_Win_Unlock calls */
int ADIOI_MPIMUTEX_Fetch_and_increment(mpimutex_t mutex,
		MPI_Offset *current, MPI_Offset increment)
{
	int mpi_err, i;
	unsigned char add = 1, remove = 0;

	unsigned char *waitlistcopy = NULL;
	MPI_Offset fpcopy = 0;

	waitlistcopy = malloc(mutex->nprocs-1);
	if (!waitlistcopy) goto err_return;

	/* add self to waitlist */
	mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, mutex->homerank, 0,
			mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	/* non-zero target displacement: contiguous region for waitlist
	 * after fp */
	mpi_err = MPI_Get(waitlistcopy, mutex->nprocs-1, MPI_BYTE,
			mutex->homerank, sizeof(MPI_Offset), 1,
			mutex->waitlisttype, mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	mpi_err = MPI_Put(&add, 1, MPI_BYTE,
			mutex->homerank, sizeof(MPI_Offset) + mutex->myrank, 1, MPI_BYTE,
			mutex->waitlistwin);

	/* shared fp value is in the same window */
	mpi_err = MPI_Get(&fpcopy, 1, mutex->fptype,
			mutex->homerank, 0, 1,
			mutex->fptype, mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	mpi_err = MPI_Win_unlock(mutex->homerank, mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	/* scan waitlist to see if lock already held */
	for (i=0; i< (mutex->nprocs -1) && waitlistcopy[i] == 0; i++);

	if (i < mutex->nprocs - 1) {
		/* someone else has the lock, so wait for notification (the new
		 * value of the file pointer)  */

		mpi_err = MPI_Recv(&fpcopy, 1, mutex->fptype, MPI_ANY_SOURCE,
				WAKEUP, mutex->comm, MPI_STATUS_IGNORE);
		if (mpi_err != MPI_SUCCESS) goto err_return;
	}
	/* now in either case we've got the latest value of the shared fp */
	*current = fpcopy;
	fpcopy += increment;

	/* remove self from waitlist */
	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, mutex->homerank, 0,
				mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	mpi_err = MPI_Get(waitlistcopy, mutex->nprocs-1, MPI_BYTE,
			mutex->homerank, sizeof(MPI_Offset), 1, mutex->waitlisttype,
			mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	mpi_err = MPI_Put(&remove, 1, MPI_BYTE,
			mutex->homerank, sizeof(MPI_Offset) + mutex->myrank, 1, MPI_BYTE,
			mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	mpi_err = MPI_Put(&fpcopy, 1, mutex->fptype,
			mutex->homerank, 0, 1,
			mutex->fptype, mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	mpi_err = MPI_Win_unlock(mutex->homerank, mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	/* check to see if anyone waiting for lock */
	for (i=0; i < (mutex->nprocs - 1) && waitlistcopy[i] == 0; i++);

	if(i < mutex->nprocs - 1) {
		int nextrank;
		/* find the next rank waiting for the lock.  we start with the
		 * rank after ours and look in order to ensure fairness.
		 */
		nextrank = mutex->myrank;
		while (nextrank < (mutex->nprocs - 1) &&
				waitlistcopy[nextrank] == 0)
		{
			nextrank++;
		}
		if (nextrank < mutex->nprocs - 1) {
			/* got a valid rnk, but 'nextrank' is off by one */
			nextrank++;
		} else {
			nextrank = 0;
			while (nextrank < mutex->myrank &&
					waitlistcopy[nextrank] == 0)
			{
				nextrank++;
			}
			assert(nextrank != mutex->myrank);
		}

		/* notify next rank */
		mpi_err = MPI_Send(&fpcopy, 1,  mutex->fptype, nextrank,
				WAKEUP, mutex->comm);
		if (mpi_err != MPI_SUCCESS) goto err_return;
	}

	free(waitlistcopy);
	return MPI_SUCCESS;

 err_return:
	printf("error!\n");
	if (waitlistcopy) free(waitlistcopy);
	return MPI_ERR_UNKNOWN;
}

/* Get/Set directly manipulate the memory window without coordinating among
 * other processors.  It is intened that the caller will make the necessary
 * effort to ensure no conflict with these calls */

int ADIOI_MPIMUTEX_Set(mpimutex_t mutex, MPI_Offset value)
{
	int mpi_err;

	mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, mutex->homerank, 0,
			mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	mpi_err = MPI_Put(&value, 1, mutex->fptype,
			mutex->homerank, 0, 1, mutex->fptype,
			mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	mpi_err = MPI_Win_unlock(mutex->homerank, mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	return MPI_SUCCESS;

err_return:
	printf("error!\n");
	return MPI_ERR_UNKNOWN;
}

int ADIOI_MPIMUTEX_Get(mpimutex_t mutex, MPI_Offset *value)
{
	int mpi_err;

	mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, mutex->homerank, 0,
			mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	mpi_err = MPI_Get(value, 1, mutex->fptype,
			mutex->homerank, 0, 1, mutex->fptype,
			mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	mpi_err = MPI_Win_unlock(mutex->homerank, mutex->waitlistwin);
	if (mpi_err != MPI_SUCCESS) goto err_return;

	return MPI_SUCCESS;

err_return:
	printf("error!\n");
	return MPI_ERR_UNKNOWN;
}

/*
 *  * vim: ts=8 sts=4 sw=4 noexpandtab
 */
