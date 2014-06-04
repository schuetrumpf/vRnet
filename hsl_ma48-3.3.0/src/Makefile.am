lib_LIBRARIES = libhsl_ma48.a

modulesdir = modules
nodist_modules_HEADERS = \
	hsl_ma48_single.mod \
	hsl_ma48_double.mod

libhsl_ma48_a_SOURCES = \
	hsl_ma48s.f90 \
	hsl_ma48d.f90 \
	common90.f90 \
	sdeps.f \
	sdeps90.f90 \
	ddeps.f \
	ddeps90.f90 \
	hsl_ma48s_ciface.f90 \
	hsl_ma48d_ciface.f90

hsl_ma48s.$(OBJEXT): sdeps90.$(OBJEXT)
hsl_ma48s.$(OBJEXT): common90.$(OBJEXT)
hsl_ma48d.$(OBJEXT): ddeps90.$(OBJEXT)
hsl_ma48d.$(OBJEXT): common90.$(OBJEXT)
sdeps90.$(OBJEXT): common90.$(OBJEXT)
ddeps90.$(OBJEXT): common90.$(OBJEXT)
hsl_ma48s_ciface.$(OBJEXT): hsl_ma48s.$(OBJEXT)
hsl_ma48d_ciface.$(OBJEXT): hsl_ma48d.$(OBJEXT)
