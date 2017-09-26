## -*- Mode: Makefile; -*-
## vim: set ft=automake :
##
## (C) 2011 by Argonne National Laboratory.
##     See COPYRIGHT in top-level directory.
##

## this file is already guarded by an "if BUILD_MPID_COMMON_DATATYPE"

mpi_core_sources +=                                    \
src/mpi/datatype/dame/darray_support.c               \
src/mpi/datatype/dame/dame.c                     \
src/mpi/datatype/dame/dame_create.c              \
src/mpi/datatype/dame/dame_create_blockindexed.c \
src/mpi/datatype/dame/dame_create_contig.c       \
src/mpi/datatype/dame/dame_create_indexed.c      \
src/mpi/datatype/dame/dame_create_pairtype.c     \
src/mpi/datatype/dame/dame_create_struct.c       \
src/mpi/datatype/dame/dame_create_vector.c       \
src/mpi/datatype/dame/segment.c                      \
src/mpi/datatype/dame/segment_count.c                \
src/mpi/datatype/dame/segment_packunpack.c           \
src/mpi/datatype/dame/subarray_support.c

# several headers are included by the rest of MPICH
AM_CPPFLAGS += -I$(top_srcdir)/src/mpi/datatype

noinst_HEADERS +=                                        \
src/mpi/datatype/dame/dame.h         \
src/mpi/datatype/dame/typesize_support.h \
src/mpi/datatype/dame/veccpy.h

