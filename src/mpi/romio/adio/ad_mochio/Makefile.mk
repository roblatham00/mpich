## -*- Mode: Makefile; -*-
## vim: set ft=automake :
##

if BUILD_AD_MOCHIO

noinst_HEADERS += adio/ad_mochio/ad_mochio.h \
                  adio/ad_mochio/ad_mochio_common.h

romio_other_sources +=                 \
    adio/ad_mochio/ad_mochio_features.c   \
    adio/ad_mochio/ad_mochio.c           \
    adio/ad_mochio/ad_mochio_common.c   \
    adio/ad_mochio/ad_mochio_open.c           \
    adio/ad_mochio/ad_mochio_close.c           \
    adio/ad_mochio/ad_mochio_io.c \
    adio/ad_mochio/ad_mochio_io_list.c \
    adio/ad_mochio/ad_mochio_io_list_classic.c

endif BUILD_AD_MOCHIO
