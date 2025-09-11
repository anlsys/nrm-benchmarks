AM_CPPFLAGS = \
	      -I $(top_builddir)/src $(OPENMP_CFLAGS) $(LIBNRM_CFLAGS) $(BLAS_CFLAGS)
AM_LDFLAGS = @MPILIBS@ $(OPENMP_CFLAGS) $(LIBNRM_LIBS) $(BLAS_LIBS)
CC=@MPICC@

mpi_noprogress_stream_triad_SOURCES = $(top_srcdir)/src/utils.c \
				      ./noprogress/stream/triad.c

mpi_noprogress_stream_scale_SOURCES = $(top_srcdir)/src/utils.c \
				      ./noprogress/stream/scale.c

mpi_noprogress_stream_add_SOURCES = $(top_srcdir)/src/utils.c \
				      ./noprogress/stream/add.c

mpi_noprogress_stream_copy_SOURCES = $(top_srcdir)/src/utils.c \
				      ./noprogress/stream/copy.c

mpi_noprogress_stream_full_SOURCES = $(top_srcdir)/src/utils.c \
				      ./noprogress/stream/full.c
bin_PROGRAMS = \
	       mpi-noprogress-stream-triad \
	       mpi-noprogress-stream-scale \
	       mpi-noprogress-stream-add \
	       mpi-noprogress-stream-copy \
	       mpi-noprogress-stream-full

