###############################################################################
# BOILERPLATE
###############################################################################

ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = autogen.sh README.markdown

AM_CPPFLAGS = -I$(top_srcdir)/src $(OPENMP_CFLAGS) $(LIBNRM_CFLAGS) $(BLAS_CFLAGS)
AM_LDFLAGS = $(OPENMP_CFLAGS) $(LIBNRM_LIBS) $(BLAS_LIBS)

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
noprogress_stream_full_SOURCES = $(UTILS_SOURCES) src/noprogress/stream/full.c

ones_solvers_cg_SOURCES = $(UTILS_SOURCES) \
			  src/progress/ones/iterative_solvers/common.h \
			  src/progress/ones/iterative_solvers/cg.c
ones_solvers_bicgstab_SOURCES = $(UTILS_SOURCES) \
			  src/progress/ones/iterative_solvers/common.h \
			  src/progress/ones/iterative_solvers/bicgstab.c

NPB_UTILS_SOURCES = src/progress/ones/npb/randdp.c
ones_npb_ep_SOURCES = $(UTILS_SOURCES) $(NPB_UTILS_SOURCES) src/progress/ones/npb/ep.c

phases_stream_full_SOURCES = $(UTILS_SOURCES) src/progress/phases/stream/full.c

phases_solvers_cg_SOURCES = $(UTILS_SOURCES) \
			  src/progress/phases/iterative_solvers/common.h \
			  src/progress/phases/iterative_solvers/cg.c

phases_solvers_bicgstab_SOURCES = $(UTILS_SOURCES) \
			  src/progress/phases/iterative_solvers/common.h \
			  src/progress/phases/iterative_solvers/bicgstab.c

bin_PROGRAMS = ones-solvers-cg \
	       ones-solvers-bicgstab \
	       ones-stream-copy \
	       ones-stream-scale \
	       ones-stream-add \
	       ones-stream-triad \
	       ones-stream-full \
	       ones-npb-ep \
	       phases-stream-full \
	       phases-solvers-cg \
	       phases-solvers-bicgstab \
	       noprogress-stream-full
