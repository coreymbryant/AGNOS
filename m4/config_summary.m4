# SYNOPSIS
#
#   Summarizes configuration settings.
#
#   AX_SUMMARIZE_CONFIG([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Outputs a summary of relevant configuration settings.
#
# LAST MODIFICATION
#
#   2010-03-24
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version............... : $PACKAGE-$VERSION
echo
echo C++ compiler.................. : $CXX
echo C++ compiler flags............ : $CXXFLAGS
echo Install dir................... : $prefix 
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo GIT revision.................. : $BUILD_VERSION
echo
echo Library Dependencies:
echo libMesh....................... : $LIBMESH_PREFIX
echo libMesh CXXFLAGS.............. : $LIBMESH_CXXFLAGS
echo libMesh INCLUDE............... : $LIBMESH_INCLUDE
echo Boost......................... : $BOOST_ROOT
echo GSL........................... : $GSL_PREFIX
#echo GRINS......................... : $GRINS_PREFIX
#echo QUESO......................... : $QUESO_PREFIX
#echo Trilinos...................... : $TRILINOS_PREFIX
#echo HDF5.......................... : $HDF5_PREFIX
#echo GLPK.......................... : $GLPK_PREFIX
echo
echo Optional Packages:
echo '  'cppunit.......................... : $enablecppunit
  if (test "x$enablecppunit" = "xyes"); then
  echo '     'CPPUNIT_CFLAGS................ : $CPPUNIT_CFLAGS
  echo '     'CPPUNIT_LIBS.................. : $CPPUNIT_LIBS
  fi
echo '  'channelflow...................... : $enablechannelflow
  #TODO
  if (test "x$enablechannelflow" = "xyes"); then
  echo '     'CHANNELFLOW_DIR............... : $CHANNELFLOW_DIR
  fi
echo '  'grins............................ : $enablegrins
  #TODO
  if (test "x$enablegrins" = "xyes"); then
  echo '     'GRINS_DIR..................... : $GRINS_DIR
  fi

echo
echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
