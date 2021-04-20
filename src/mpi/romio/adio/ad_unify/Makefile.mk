##
## Copyright (C) by Argonne National Laboratory
##     See COPYRIGHT in top-level directory
##

if BUILD_AD_UNIFY

noinst_HEADERS +=                       \
    adio/ad_unify/ad_unify.h

romio_other_sources +=                          \
    adio/ad_unify/ad_unify.c        \
    adio/ad_unify/ad_unify_common.c \
    adio/ad_unify/ad_unify_open.c \
    adio/ad_unify/ad_unify_close.c \
    adio/ad_unify/ad_unify_io.c

endif BUILD_AD_UNIFY
