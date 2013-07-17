include $(top_srcdir)/doxygen/aminclude.am

AUTOMAKE_OPTIONS = foreign
#ACLOCAL_AMFLAGS  = 

SUBDIRS          = src test examples doxygen 
#SUBDIRS 				+= contrib
#EXTRA_DIST       =

# Tools in the auxiliary directory

AUX_DIST         = build-aux/install-sh
AUX_DIST        += build-aux/missing
AUX_DIST        += build-aux/config.guess
AUX_DIST        += build-aux/config.sub
AUX_DIST        += build-aux/depcomp
AUX_DIST        += build-aux/ltmain.sh

# Files to be deleted by 'make maintainer-clean'
MAINTAINERCLEANFILES = aclocal.m4                               \
                       aminclude.am                             \
                       autom4te.cache/*                         \
                       $(AUX_DIST)                              \
                       config.log                               \
                       config.status                            \
                       config.sub                               \
                       configure                                \
                       config.h                                 \
                       config.h.in                              \
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
