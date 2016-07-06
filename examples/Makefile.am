
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
							 $(GSL_CFLAGS) \
							 $(LIBMESH_CPPFLAGS)

#------------------------------------
# Flags we need to pass to the linker
#------------------------------------
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
						 $(BOOST_PROGRAM_OPTIONS_LIBS) \
						 $(GSL_LIBS) \
						 $(LIBMESH_LDFLAGS)

#-------------------------------------
# Catenary example
#-------------------------------------
catenarydir = $(prefix)/examples/catenary
catenary_SCRIPTS = $(top_builddir)/examples/catenary/run.sh
catenary_DATA = $(top_srcdir)/examples/catenary/catenary.in
EXTRA_DIST += $(catenary_DATA)

#-------------------------------------
# Diffusion
#-------------------------------------
diffusiondir = $(prefix)/examples/diffusion
diffusion_SCRIPTS = $(top_builddir)/examples/diffusion/run.sh
diffusion_DATA = $(top_srcdir)/examples/diffusion/diffusion.in
EXTRA_DIST += $(diffusion_DATA)

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


#-------------------------------------
# grins
#-------------------------------------
if ENABLE_GRINS
grinsdir = $(prefix)/examples/grins
grins_SCRIPTS = $(top_builddir)/examples/grins/run.sh
grins_DATA = $(top_srcdir)/examples/grins/grinsExample.in $(top_srcdir)/examples/grins/grins.in $(top_srcdir)/examples/grins/fat_box_with_hole.e
EXTRA_DIST += $(grins_DATA)
endif

MAINTAINERCLEANFILES = Makefile.in \
											 catenary/testFile.dat \
											 diffusion/testFile.dat \
											 catenaryLibmesh/testFile.dat \
											 catenaryLibmesh/*.txt \
											 catenaryLibmesh/*.gp.* \
											 viscousBurgers/testFile.dat \
											 viscousBurgers/*.txt \
											 viscousBurgers/*.gp.* \
											 grins/*.dat  \
											 grins/output/* \
											 grins/*.exo \
											 grins/*.pdf \
											 grins/*.txt \
											 grins/*.out*

