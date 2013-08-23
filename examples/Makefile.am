EXTRA_DIST = 
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

#-------------------------------------
# Catenary example
#-------------------------------------
catenarydir = $(prefix)/examples/catenary
catenary_SCRIPTS = $(top_builddir)/examples/catenary/run.sh
catenary_DATA = $(top_srcdir)/examples/catenary/catenary.in
EXTRA_DIST += $(catenary_DATA)

#-------------------------------------
# Catenary libmesh example
#-------------------------------------
catenarylibmeshdir = $(prefix)/examples/catenaryLibmesh
catenarylibmesh_DATA = $(top_srcdir)/examples/catenaryLibmesh/catenaryLibmesh.in
catenarylibmesh_SCRIPTS = $(top_builddir)/examples/catenaryLibmesh/run.sh
EXTRA_DIST += $(catenarylibmesh_DATA)

#-------------------------------------
# Viscous Burgers example
#-------------------------------------
viscousburgersdir = $(prefix)/examples/viscousBurgers
viscousburgers_DATA = $(top_srcdir)/examples/viscousBurgers/viscousBurgers.in
viscousburgers_SCRIPTS = $(top_builddir)/examples/viscousBurgers/run.sh
EXTRA_DIST += $(viscousburgers_DATA)


MAINTAINERCLEANFILES = Makefile.in \
											 catenary/testFile.dat \
											 catenaryLibmesh/testFile.dat \
											 viscousBurgers/testFile.dat 

