#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

check_PROGRAMS = polynomialTest \
								 catenaryTest \
								 evaluationTest \
								 parallelCatenary \
								 catenaryLibmeshTest \
								 viscousBurgersTest \
								 elementTest

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(top_srcdir)/include \
							-I$(top_srcdir)/include/physics \
							-I$(top_srcdir)/include/surrogate \
							-I$(top_srcdir)/contrib/sandia_rules

#-------------------------------------
# External library include directories
#-------------------------------------
AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
							 $(GSL_CFLAGS) \
							 $(LIBMESH_CPPFLAGS) \
							 $(HDF5_CFLAGS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
             $(BOOST_PROGRAM_OPTIONS_LIBS) \
						 $(GSL_LIBS) \
             $(LIBMESH_LIBS) \
             $(LIBMESH_LDFLAGS) \
						 ${HDF5_CXXLIBS} \
						 $(CPPUNIT_LIBS) \
						 -L$(top_srcdir)/src/ -lagnos

TESTS =
shellfiles_src = 

#------------------------------------
# libmesh test
#------------------------------------
check_PROGRAMS       += LibmeshTest
LibmeshTest_SOURCES   = LibmeshTest.C
TESTS += LibmeshTest.sh
shellfiles_src += LibmeshTest.sh

#------------------------------------
# HDF5 test
#------------------------------------
check_PROGRAMS       += hdf5Test
hdf5Test_SOURCES   = hdf5Test.C
TESTS += hdf5Test.sh
shellfiles_src += hdf5Test.sh

#------------------------------------
# Evaluator test
#------------------------------------
check_PROGRAMS       += EvaluatorTest
EvaluatorTest_SOURCES   = EvaluatorTest.C
TESTS += EvaluatorTest.sh
shellfiles_src += EvaluatorTest.sh

#------------------------------------
# Driver test
#------------------------------------
check_PROGRAMS       += DriverTest
DriverTest_SOURCES   = DriverTest.C
TESTS += DriverTest.sh
shellfiles_src += DriverTest.sh



polynomialTest_LDADD = $(CPPUNIT_LIBS)
polynomialTest_SOURCES 			= testDriver.C PolynomialTest.cpp
catenaryTest_SOURCES 				= testDriver.C CatenaryTest.cpp
evaluationTest_SOURCES 			= testDriver.C EvaluationTest.cpp
parallelCatenary_SOURCES 		= testDriver.C ParallelCatenary.cpp
catenaryLibmeshTest_SOURCES = testDriver.C LibmeshCatenaryTest.cpp
viscousBurgersTest_SOURCES  = testDriver.C ViscousBurgersTest.cpp
elementTest_SOURCES  				= testDriver.C ElementTest.cpp


TESTS += polynomialTest.sh
TESTS += catenaryTest.sh
TESTS += evaluationTest.sh
TESTS += parallelCatenary.sh
TESTS += catenaryLibmeshTest.sh
TESTS += viscousBurgersTest.sh
TESTS += elementTest.sh

shellfiles_src += polynomialTest.sh \
								 catenaryTest.sh \
								 evaluationTest.sh \
								 parallelCatenary.sh \
								 catenaryLibmeshTest.sh \
								 viscousBurgersTest.sh \
								 elementTest.sh 

#------------------------------------
# channel flow test
#------------------------------------
if ENABLE_CHANNELFLOW
check_PROGRAMS       += ChannelTest
ChannelTest_SOURCES   = ChannelTest.C
ChannelTest_CXXFLAGS  = -I$(CHANNELFLOW_DIR)/include
ChannelTest_LDADD     = -L$(CHANNELFLOW_DIR)/lib/ -lchannelflow
TESTS += ChannelTest.sh
shellfiles_src += ChannelTest.sh
endif

#------------------------------------
# Grins test
#------------------------------------
if ENABLE_GRINS
check_PROGRAMS       += GrinsTest
GrinsTest_SOURCES   = GrinsTest.C
GrinsTest_CXXFLAGS  = -I$(GRINS_DIR)/include
GrinsTest_LDADD     = -L$(GRINS_DIR)/lib/ -lgrins
TESTS += GrinsTest.sh
shellfiles_src += GrinsTest.sh
endif

# Want these put with the distro so we can run make check
EXTRA_DIST = $(shellfiles_src) input_files test_data 

CLEANFILES = *.out* test.log *.h5

