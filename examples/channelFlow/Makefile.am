bin_PROGRAMS  = channelFlow

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
						 $(LIBMESH_LDFLAGS) \
						 -L$(CHANNELFLOW_DIR)/lib/ -lchannelflow



channelFlow_SOURCES = channelFlow.C


CLEANFILES = *.out
CLEANFILES += output/*

