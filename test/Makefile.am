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
							-I$(top_srcdir)/contrib/sandia_rules

#-------------------------------------
# External library include directories
#-------------------------------------
AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
							 $(LIBMESH_CPPFLAGS)

# AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)
# LIBS = $(CPPUNIT_LIBS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
             $(BOOST_PROGRAM_OPTIONS_LIBS) \
             $(LIBMESH_LIBS) \
             $(LIBMESH_LDFLAGS) \
						 $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS) \
						 $(BOOST_UNIT_TEST_FRAMEWORK_LIBS) \
						 $(CPPUNIT_LIBS) \
						 $(BOOST_LDFLAGS) \
						 $(BOOST_LIBS)

TESTS =

#------------------------------------
# libmesh test
#------------------------------------
check_PROGRAMS       += LibmeshTest
LibmeshTest_SOURCES   = LibmeshTest.C
# LibmeshTest_CXXFLAGS  = -I$(CHANNELFLOW_DIR)/include
# LibmeshTest_LDADD     = -L$(CHANNELFLOW_DIR)/lib/ -lchannelflow
TESTS += LibmeshTest

#------------------------------------
# channel flow test
#------------------------------------
if ENABLE_CHANNELFLOW
check_PROGRAMS       += ChannelTest
ChannelTest_SOURCES   = ChannelTest.C
ChannelTest_CXXFLAGS  = -I$(CHANNELFLOW_DIR)/include
ChannelTest_LDADD     = -L$(CHANNELFLOW_DIR)/lib/ -lchannelflow
TESTS += ChannelTest
endif



polynomialTest_LDADD = $(CPPUNIT_LIBS)
polynomialTest_SOURCES 			= driver.C PolynomialTest.cpp
catenaryTest_SOURCES 				= driver.C CatenaryTest.cpp
evaluationTest_SOURCES 			= driver.C EvaluationTest.cpp
parallelCatenary_SOURCES 		= driver.C ParallelCatenary.cpp
catenaryLibmeshTest_SOURCES = driver.C LibmeshCatenaryTest.cpp
viscousBurgersTest_SOURCES  = driver.C ViscousBurgersTest.cpp
elementTest_SOURCES  				= driver.C ElementTest.cpp


TESTS += polynomialTest.sh
TESTS += catenaryTest.sh
TESTS += evaluationTest.sh
TESTS += parallelCatenary.sh
TESTS += catenaryLibmeshTest.sh
TESTS += viscousBurgersTest.sh
TESTS += elementTest.sh

shellfiles_src = polynomialTest.sh \
								 catenaryTest.sh \
								 evaluationTest.sh \
								 parallelCatenary.sh \
								 catenaryLibmeshTest.sh \
								 viscousBurgersTest.sh \
								 elementTest.sh 

# Want these put with the distro so we can run make check
EXTRA_DIST = $(shellfiles_src) input_files test_data 


