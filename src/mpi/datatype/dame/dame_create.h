/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */

/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef DAME_CREATE_H
#define DAME_CREATE_H

void MPIR_Dame_create(MPI_Datatype type, Dame ** dlp_p, MPI_Aint * size, int *depth);

int MPIR_Dame_create_contiguous(DAME_Count count,
                                MPI_Datatype oldtype, Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_vector(DAME_Count count,
                            DAME_Size blocklength,
                            MPI_Aint stride,
                            int strideinbytes,
                            MPI_Datatype type,
                            MPI_Datatype oldtype, Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_blockindexed(DAME_Count count,
                                  DAME_Size blklen,
                                  const void *disp_array,
                                  int dispinbytes,
                                  MPI_Datatype type,
                                  MPI_Datatype oldtype, Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_indexed(DAME_Count count,
                             const DAME_Size * blklen,
                             const void *disp_array,
                             int dispinbytes,
                             MPI_Datatype type,
                             MPI_Datatype oldtype, Dame ** dl, MPI_Aint * size, int *depth);

int MPIR_Dame_create_struct(DAME_Count count,
                            const int *blklen_array,
                            const MPI_Aint * disp_array,
                            MPI_Datatype type,
                            const MPI_Datatype * oldtype_array,
                            Dame ** dl, MPI_Aint * size, int *depth);

/* we bump up the size of the blocklength array because create_struct might use
 * create_indexed in an optimization, and in course of doing so, generate a
 * request of a large blocklength. */
int MPIR_Dame_create_pairtype(MPI_Datatype type, Dame ** dl, MPI_Aint * size, int *depth);

/* Helper functions for dataloop construction */
int MPIR_Type_convert_subarray(int ndims,
                               int *array_of_sizes,
                               int *array_of_subsizes,
                               int *array_of_starts,
                               int order, MPI_Datatype oldtype, MPI_Datatype * newtype);
int MPIR_Type_convert_darray(int size,
                             int rank,
                             int ndims,
                             int *array_of_gsizes,
                             int *array_of_distribs,
                             int *array_of_dargs,
                             int *array_of_psizes,
                             int order, MPI_Datatype oldtype, MPI_Datatype * newtype);

DAME_Count MPIR_Type_blockindexed_count_contig(DAME_Count count,
                                               DAME_Count blklen,
                                               const void *disp_array,
                                               int dispinbytes, DAME_Offset old_extent);

DAME_Count MPIR_Type_indexed_count_contig(DAME_Count count,
                                          const DAME_Size * blocklength_array,
                                          const void *displacement_array,
                                          int dispinbytes, DAME_Offset old_extent);


#endif
