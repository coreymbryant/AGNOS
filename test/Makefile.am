#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

check_PROGRAMS = polynomialSurrogates 

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(top_srcdir)/include \
							-I$(top_srcdir)/contrib/sandia_rules

#-------------------------------------
# External library include directories
#-------------------------------------
AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
							 $(LIBMESH_CPPFLAGS)
# AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
						 $(BOOST_PROGRAM_OPTIONS_LIBS) \
						 -lboost_unit_test_framework

polynomialSurrogates_SOURCES =$(top_srcdir)/test/PolynomialTest.cpp

TESTS = $(top_builddir)/test/polynomialSurrogates
