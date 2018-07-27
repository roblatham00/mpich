/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *   $Id: ad_testfs_seek.c,v 1.6 2004/10/25 18:46:00 robl Exp $
 *
 *   Copyright (C) 2001 University of Chicago.
 *   See COPYRIGHT notice in top-level directory.
 */

#include "ad_logfs.h"
#include "adioi.h"
#include "adio_extern.h"

/* ADIOI_LOGFS_SeekIndividual()
 *
 * Implements SEEK_SET only (and doesn't test for whence type); all
 * other types of whence must be converted before calling this.
 *
 * Returns an absolute offset in bytes.  The offset passed into the call is in
 * terms of the etype relative to the filetype, so some calculations are
 * necessary.
 */
ADIO_Offset ADIOI_LOGFS_SeekIndividual(ADIO_File fd, ADIO_Offset offset,
                                       int whence, int *error_code)
{
    int myrank, nprocs;

    *error_code = MPI_SUCCESS;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);
    /*   FPRINTF(stdout, "[%d/%d] ADIOI_LOGFS_SeekIndividual called on %s\n",
     * myrank, nprocs, fd->filename);
     */
    /* don't need to pass seek operations to the slave;
     * The file positions are reset on a view switch, so we just need
     * to keep track of the seek position in memory.
     *
     * ADIOI_Gen_SeekIndividual does this for us
     */


    /* have ADIO_GEN_SeekIndividual update the individual filepointer
     * (fd->fp_ind) */
    return ADIOI_GEN_SeekIndividual(fd, offset, whence, error_code);

}
