
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
							 $(GSL_CFLAGS) \
							 $(LIBMESH_CPPFLAGS) \
							 $(HDF5_CFLAGS)
# AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
# LIBS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
#              $(BOOST_PROGRAM_OPTIONS_LIBS) \
#              $(LIBMESH_LIBS) \
#              $(LIBMESH_LDFLAGS)
LIBS = $(LIBMESH_LIBS) \
			 $(GSL_LIBS) \
			 $(HDF5_CXXLIBS) 

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------
lib_LTLIBRARIES = libagnos.la
libagnos_la_SOURCES = $(top_srcdir)/contrib/sandia_rules/sandia_rules.cpp \
											surrogate/SurrogateModelBase.C \
											surrogate/SurrogateModel.C \
											surrogate/SurrogateEvaluator.C \
											surrogate/Parameter.C \
											surrogate/QuadratureRule.C \
											surrogate/QuadratureTensorProduct.C \
											surrogate/SurrogatePseudoSpectral.C \
											surrogate/PseudoSpectralTensorProduct.C \
											surrogate/EvaluatorPseudoSpectral.C \
											physics/PhysicsUser.C \
											physics/PhysicsCatenary.C \
											physics/PhysicsLibmesh.C \
											physics/PhysicsViscousBurgers.C \
											physics/PhysicsCatenaryLibmesh.C \
											physics/PhysicsDiffusion.C 

											
						 
if ENABLE_GRINS
AM_CPPFLAGS += -I$(GRINS_DIR)/include -I$(GRVY_DIR)/include
LIBS += -L$(GRINS_DIR)/lib/ -lgrins
libagnos_la_SOURCES += physics/PhysicsGrins.C 
endif

libagnos_la_SOURCES += ioHandler.C Element.C Driver.C 

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
