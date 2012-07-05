## -*- Mode: Makefile; -*-
## vim: set ft=automake :
##
## (C) 2011 by Argonne National Laboratory.
##     See COPYRIGHT in top-level directory.
##

if BUILD_NEMESIS_NETMOD_SCIF

lib_lib@MPILIBNAME@_la_SOURCES +=                                   \
    src/mpid/ch3/channels/nemesis/nemesis/netmod/scif/scif_finalize.c \
    src/mpid/ch3/channels/nemesis/nemesis/netmod/scif/scif_init.c     \
    src/mpid/ch3/channels/nemesis/nemesis/netmod/scif/scif_send.c     \
    src/mpid/ch3/channels/nemesis/nemesis/netmod/scif/scifrw.c  \
    src/mpid/ch3/channels/nemesis/nemesis/netmod/scif/scifsm.c

noinst_HEADERS +=                                                \
    src/mpid/ch3/channels/nemesis/nemesis/netmod/scif/scifrw.h    \
    src/mpid/ch3/channels/nemesis/nemesis/netmod/scif/scif_impl.h

endif BUILD_NEMESIS_NETMOD_SCIF
