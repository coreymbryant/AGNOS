# The lib headers
libagnos_includedir  = $(prefix)/include

libagnos_include_HEADERS = \
											$(top_srcdir)/contrib/sandia_rules/sandia_rules.hpp \
											surrogate/Parameter.h \
											surrogate/QuadratureRule.h \
											surrogate/QuadratureTensorProduct.h \
											surrogate/SurrogateModel.h \
											surrogate/SurrogatePseudoSpectral.h \
											surrogate/PseudoSpectralTensorProduct.h \
											physics/PhysicsModel.h \
											physics/PhysicsUser.h \
											physics/PhysicsCatenary.h \
											physics/PhysicsLibmesh.h \
											physics/BurgersSystem.h \
											physics/CatenarySystem.h \
											physics/DiffusionSystem.h \
											physics/PhysicsViscousBurgers.h \
											physics/PhysicsCatenaryLibmesh.h \
											physics/PhysicsDiffusion.h 

if ENABLE_CHANNELFLOW
libagnos_include_HEADERS += physics/PhysicsChannelFlow.h
endif

libagnos_include_HEADERS += Element.h Driver.h 

# Headers required from channel exec that are not in libchannelflow
agnos_includedir = $(prefix)/include
agnos_include_HEADERS = 


