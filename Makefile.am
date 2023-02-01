###############################################################################
# BOILERPLATE
###############################################################################

ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = autogen.sh README.markdown

AM_CPPFLAGS = -I$(top_srcdir)/src $(OPENMP_CFLAGS) $(LIBNRM_CFLAGS)
AM_LDFLAGS = $(OPENMP_CFLAGS) $(LIBNRM_LIBS)

include_HEADERS=src/nrm-benchmarks.h
UTILS_SOURCES = src/utils.c

###############################################################################
# BENCHMARKS
###############################################################################

ones_stream_copy_SOURCES = $(UTILS_SOURCES) src/progress/ones/stream/copy.c
ones_stream_scale_SOURCES = $(UTILS_SOURCES) src/progress/ones/stream/scale.c
ones_stream_add_SOURCES = $(UTILS_SOURCES) src/progress/ones/stream/add.c
ones_stream_triad_SOURCES = $(UTILS_SOURCES) src/progress/ones/stream/triad.c
ones_stream_full_SOURCES = $(UTILS_SOURCES) src/progress/ones/stream/full.c

NPB_UTILS_SOURCES = src/progress/ones/npb/randdp.c
ones_npb_ep_SOURCES = $(UTILS_SOURCES) $(NPB_UTILS_SOURCES) src/progress/ones/npb/ep.c

bin_PROGRAMS = ones-stream-copy \
	       ones-stream-scale \
	       ones-stream-add \
	       ones-stream-triad \
	       ones-stream-full \
	       ones-npb-ep
