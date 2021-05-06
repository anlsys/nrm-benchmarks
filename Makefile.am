###############################################################################
# BOILERPLATE
###############################################################################

ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = autogen.sh README.markdown

AM_CPPFLAGS = -Isrc $(OPENMP_CFLAGS) $(LIBNRM_CFLAGS)
AM_LDFLAGS = $(OPENMP_CFLAGS) $(LIBNRM_LIBS)

include_HEADERS=src/nrm-benchmarks.h
UTILS_SOURCES = src/utils.c

###############################################################################
# BENCHMARKS
###############################################################################

progress_ones_copy_SOURCES = $(UTILS_SOURCES) src/progress/ones/copy.c

bin_PROGRAMS = progress-ones-copy
