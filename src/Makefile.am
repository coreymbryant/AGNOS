
#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(top_srcdir) \
              -I$(top_builddir)/include \
							-I$(top_srcdir)/include/physics \
							-I$(top_srcdir)/include/surrogate \
							-I$(top_srcdir)/contrib/sandia_rules \
							-I$(top_srcdir)/contrib/getpot 

#-------------------------------------
# External library include directories
#-------------------------------------
AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
							 $(LIBMESH_CPPFLAGS)
# AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
# LIBS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
#              $(BOOST_PROGRAM_OPTIONS_LIBS) \
#              $(LIBMESH_LIBS) \
#              $(LIBMESH_LDFLAGS)
LIBS = $(LIBMESH_LIBS)

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------
lib_LTLIBRARIES = libagnos.la
libagnos_la_SOURCES = $(top_srcdir)/contrib/sandia_rules/sandia_rules.cpp \
											surrogate/SurrogateModel.C \
											surrogate/Parameter.C \
											surrogate/QuadratureRule.C \
											surrogate/QuadratureTensorProduct.C \
											surrogate/SurrogatePseudoSpectral.C \
											surrogate/PseudoSpectralTensorProduct.C \
											physics/PhysicsUser.C \
											physics/PhysicsCatenary.C \
											physics/PhysicsLibmesh.C \
											physics/PhysicsViscousBurgers.C \
											physics/PhysicsCatenaryLibmesh.C \
											physics/PhysicsDiffusion.C 

											
						 
if ENABLE_CHANNELFLOW
AM_CPPFLAGS += -I$(CHANNELFLOW_DIR)/include
LIBS += -L$(CHANNELFLOW_DIR)/lib/ -lchannelflow
libagnos_la_SOURCES += physics/PhysicsChannelFlow.C 
endif

libagnos_la_SOURCES += Element.C Driver.C 

#-------------------------------------
# Standalone executable
#-------------------------------------
bin_PROGRAMS = agnos 
agnos_SOURCES = agnos.C 
agnos_LDADD = libagnos.la

METASOURCES   = AUTO
#
# Required for AX_AM_MACROS
##@INC_AMINCLUDE@
MAINTAINERCLEANFILES = Makefile.in
