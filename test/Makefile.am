#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

check_PROGRAMS = polynomialTest catenaryTest evaluationTest

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
						 $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS) \
						 $(BOOST_UNIT_TEST_FRAMEWORK_LIBS) 

polynomialTest_SOURCES 	=$(top_srcdir)/test/PolynomialTest.cpp
catenaryTest_SOURCES 		=$(top_srcdir)/test/CatenaryTest.cpp
evaluationTest_SOURCES 	=$(top_srcdir)/test/EvaluationTest.cpp

TESTS  = $(top_builddir)/test/polynomialTest.sh
TESTS += $(top_builddir)/test/catenaryTest.sh
TESTS += $(top_builddir)/test/evaluationTest.sh