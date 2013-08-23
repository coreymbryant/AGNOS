include $(top_srcdir)/doxygen/aminclude.am

AUTOMAKE_OPTIONS = foreign 
ACLOCAL_AMFLAGS  = -I m4 -I m4/common

SUBDIRS          = src test examples doxygen 
#SUBDIRS 				+= contrib
# EXTRA_DIST       = CHANGES LICENSE COPYING docs contrib
EXTRA_DIST       = docs contrib

includedir = $(prefix)/include/
include_HEADERS = $(top_builddir)/agnos_config.h
#
# Eliminate .git directories in dist tarball
dist-hook:
	rm -rf `find $(distdir)/ -name .git`

# Tools in the auxiliary directory

AUX_DIST         = build-aux/install-sh
AUX_DIST        += build-aux/missing
AUX_DIST        += build-aux/config.guess
AUX_DIST        += build-aux/config.sub
AUX_DIST        += build-aux/depcomp
AUX_DIST        += build-aux/ltmain.sh


# Additional files to be deleted by 'make distclean'
DISTCLEANFILES = _configs.sed

# Files to be deleted by 'make maintainer-clean'
MAINTAINERCLEANFILES = aclocal.m4                               \
                       aminclude.am                             \
                       autom4te.cache/*                         \
                       $(AUX_DIST)                              \
                       config.log                               \
                       config.status                            \
                       config.sub                               \
                       configure                                \
                       agnos_config.h                           \
                       agnos_config.h.tmp.in                    \
                       agnos_config.h.tmp.in~                   \
                       Makefile.in                              \
                       src/Makefile.in                          \
											 test/Makefile.in 												\
											 examples/Makefile.in 										\
											 examples/catenary/Makefile.in 						\
											 examples/catenaryLibmesh/Makefile.in 		\
											 examples/viscousBurgers/Makefile.in 			\
											 doxygen/Makefile.in


MOSTLYCLEANFILES     = $(DX_CLEANFILES) \
                       docs/lcov
