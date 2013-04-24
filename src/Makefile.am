#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS = agnos 

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(top_srcdir)/include \
							-I$(top_srcdir)/contrib/sandia_rules

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

#
# Required for AX_AM_MACROS
##@INC_AMINCLUDE@
