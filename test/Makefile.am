#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

check_PROGRAMS = polynomialTest \
								 catenaryTest \
								 evaluationTest \
								 parallelCatenary \
								 catenaryLibmeshTest

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(top_srcdir)/include \
							-I$(top_srcdir)/contrib/sandia_rules

#-------------------------------------
# External library include directories
#-------------------------------------
AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
							 $(LIBMESH_CPPFLAGS) \
							 $(CPPUNIT_CFLAGS)
# AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)
LIBS = $(CPPUNIT_LIBS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
						 $(BOOST_PROGRAM_OPTIONS_LIBS) \
						 $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS) \
						 $(BOOST_UNIT_TEST_FRAMEWORK_LIBS) \
						 $(LIBMESH_LDFLAGS)

polynomialTest_SOURCES 			= libmeshTests.C PolynomialTest.cpp
catenaryTest_SOURCES 				= libmeshTests.C CatenaryTest.cpp
evaluationTest_SOURCES 			= libmeshTests.C EvaluationTest.cpp
parallelCatenary_SOURCES 		= libmeshTests.C ParallelCatenary.cpp
catenaryLibmeshTest_SOURCES = libmeshTests.C LibmeshCatenaryTest.cpp


TESTS  = polynomialTest.sh
TESTS += catenaryTest.sh
TESTS += evaluationTest.sh
TESTS += parallelCatenary.sh
TESTS += catenaryLibmeshTest

shellfiles_src = polynomialTest.sh \
								 catenaryTest.sh \
								 evaluationTest.sh \
								 parallelCatenary.sh

# Want these put with the distro so we can run make check
EXTRA_DIST = $(shellfiles_src) input_files test_data 


