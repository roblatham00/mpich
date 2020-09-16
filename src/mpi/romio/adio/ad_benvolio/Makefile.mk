## -*- Mode: Makefile; -*-
## vim: set ft=automake :
##

if BUILD_AD_BENVOLIO

noinst_HEADERS += adio/ad_benvolio/ad_bv.h \
                  adio/ad_benvolio/ad_bv_common.h

romio_other_sources +=                 \
    adio/ad_benvolio/ad_bv_features.c   \
    adio/ad_benvolio/ad_bv.c           \
    adio/ad_benvolio/ad_bv_common.c   \
    adio/ad_benvolio/ad_bv_open.c           \
    adio/ad_benvolio/ad_bv_close.c           \
    adio/ad_benvolio/ad_bv_io.c \
    adio/ad_benvolio/ad_bv_io_list.c \
    adio/ad_benvolio/ad_bv_io_list_classic.c \
    adio/ad_benvolio/ad_bv_fcntl.c \
    adio/ad_benvolio/ad_bv_delete.c \
    adio/ad_benvolio/ad_bv_flush.c \
    adio/ad_benvolio/ad_bv_resize.c

endif BUILD_AD_BENVOLIO
