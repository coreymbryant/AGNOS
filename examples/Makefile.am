# AUTOMAKE_OPTIONS=subdir-objects
#--------------------------------------
SUBDIRS          = catenary \
									 catenaryLibmesh \
									 viscousBurgers
#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = -I$(top_srcdir)/include \
							-I$(top_srcdir)/contrib/sandia_rules \
							-I$(LIBMESH_CPPFLAGS)

#-------------------------------------
# External library include directories
#-------------------------------------
AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
							 $(LIBMESH_CPPFLAGS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
						 $(BOOST_PROGRAM_OPTIONS_LIBS) \
						 $(LIBMESH_LDFLAGS)

# #-------------------------------------
# # Catenary example
# #-------------------------------------
# catenarydir = $(prefix)/examples/catenary

# catenary_DATA = $(top_srcdir)/examples/catenary/catenary.in

# catenary_PROGRAMS = catenary

# catenary_SOURCES = $(top_srcdir)/examples/catenary/runCatenary.C \
# 									 $(top_srcdir)/src/agnos.C

MAINTAINERCLEANFILES = Makefile.in

