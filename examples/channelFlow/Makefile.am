bin_PROGRAMS  = 

EXTRA_DIST = 

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




channelFlowdir = $(prefix)/examples/channelFlow
channelFlow_SCRIPTS = $(top_builddir)/examples/channelFlow/run.sh
channelFlow_DATA = $(top_srcdir)/examples/channelFlow/channelFlow.in $(top_srcdir)/examples/channelFlow/flow.in $(top_srcdir)/examples/channelFlow/channel_Re944_grid.xda
EXTRA_DIST += $(channelFlow_DATA)


CLEANFILES = *.out
CLEANFILES += output/*

