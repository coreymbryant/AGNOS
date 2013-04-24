include $(top_srcdir)doxygen/aminclude.am

AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS  = -I m4

SUBDIRS          = src src/contrib test doxygen 
#EXTRA_DIST       =

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
											 test/Makefile.in 


