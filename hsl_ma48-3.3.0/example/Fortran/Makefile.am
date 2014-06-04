AM_FCFLAGS = -I $(top_builddir)/src

TESTS_ENVIRONMENT = $(srcdir)/runtest.sh
EXTRA_DIST = runtest.sh
TESTS = \
	hsl_ma48ss \
	hsl_ma48ds

check_PROGRAMS = \
	hsl_ma48ss \
	hsl_ma48ds

hsl_ma48ss_SOURCES = hsl_ma48ss.f90
hsl_ma48ss_LDADD = $(top_builddir)/src/libhsl_ma48.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma48ss.output hsl_ma48ss.data
hsl_ma48ds_SOURCES = hsl_ma48ds.f90
hsl_ma48ds_LDADD = $(top_builddir)/src/libhsl_ma48.a $(LAPACK_LIBS) $(BLAS_LIBS) $(METIS_LIBS) $(FCLIBS)
EXTRA_DIST += hsl_ma48ds.output hsl_ma48ds.data
