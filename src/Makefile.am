#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS = agnos 

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(top_srcdir)/include \
							-I$(top_srcdir)/contrib/sandia_rules \
							-I$(top_srcdir)/contrib/getpot \
							-I$(LIBMESH_CPPFLAGS)

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
             $(LIBMESH_LIBS) \
             $(LIBMESH_LDFLAGS)
						 
if ENABLE_CHANNELFLOW
AM_CPPFLAGS += -I$(CHANNELFLOW_DIR)/include
AM_LDFLAGS += -L$(CHANNELFLOW_DIR)/lib/ -lchannelflow
endif
						 
						 

# run_SOURCES =$(top_srcdir)/src/run.C \
# 							 $(top_srcdir)/src/agnos.C
# agnos_LDADD = $(LIBMESH_LIBS)
agnos_SOURCES = agnos.C

#
# Required for AX_AM_MACROS
##@INC_AMINCLUDE@
MAINTAINERCLEANFILES = Makefile.in
