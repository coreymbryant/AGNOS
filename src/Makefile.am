#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS = agnos 

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(top_srcdir)/include \
							-I$(top_srcdir)/contrib

#-------------------------------------
# External library include directories
#-------------------------------------
AM_CPPFLAGS += $(BOOST_CPPFLAGS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
						 $(BOOST_PROGRAM_OPTIONS_LIBS)

agnos_SOURCES =$(top_srcdir)/src/agnos.C

#------------------------------------
# Problem sources
#------------------------------------
#oedH2O_SOURCES = $(top_srcdir)/src/computeH2O.C \
#								 $(top_srcdir)/src/main.C \
#								 $(top_srcdir)/extra/inputs.C \
#								 $(top_srcdir)/extra/example_likelihood.C \
#								 $(top_srcdir)/extra/daxpy.c \
#								 $(top_srcdir)/extra/ddott.c \
#								 $(top_srcdir)/extra/dgefa.c \
#								 $(top_srcdir)/extra/dgesl.c \
#								 $(top_srcdir)/extra/dscal.c \
#								 $(top_srcdir)/extra/idamax.c \
#								 $(top_srcdir)/extra/lsoda.c 

#------------------------------------
# Problem sources
#------------------------------------
convectiveDispersion_SOURCES = $(top_srcdir)/src/computeConvectiveDispersion.C \
                     $(top_srcdir)/src/main.C 

convectiveDispersion_CPPFLAGS = $(AM_CPPFLAGS) \
																-I$(HOME)/bin/eigen/include/eigen3/
 
#------------------------------------
# Problem sources
#------------------------------------
cdGP_SOURCES = $(top_srcdir)/src/computeConvectiveDispersion_gp.C \
                     $(top_srcdir)/src/main.C 

cdGP_CPPFLAGS = $(AM_CPPFLAGS) \
																-I$(HOME)/bin/eigen/include/eigen3/

#------------------------------------
# Problem sources
#------------------------------------
cdBiasGP_SOURCES = $(top_srcdir)/src/computeConvectiveDispersion_biasGP.C \
                     $(top_srcdir)/src/main.C 

cdBiasGP_CPPFLAGS = $(AM_CPPFLAGS) \
																-I$(HOME)/bin/eigen/include/eigen3/
#
# Required for AX_AM_MACROS
@INC_AMINCLUDE@
