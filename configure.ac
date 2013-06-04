#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(agnos, 0.1.0, cbryant@ices.utexas.edu)
AC_CONFIG_HEADER(config.h)
#AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE(agnos,0.1.0)

AC_PROG_LIBTOOL

#------------------------------
# Checks for compilers and MPI
#------------------------------
AC_PROG_CC
AC_LANG([C])
ACX_MPI([CC="$MPICC"], [AC_MSG_ERROR([Could not find MPI.])])

AC_PROG_CXX
AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], [AC_MSG_ERROR([Could not find MPI.])])

#---------------------------------------------------------
# Checks for library prerequisites for other libraries...
#---------------------------------------------------------

#----------------
# Check for boost
#----------------
BOOST_REQUIRE([1.36])
BOOST_MATH()
BOOST_PROGRAM_OPTIONS()
BOOST_TEST()

dnl----------------
dnl Libmesh Config
dnl----------------
AX_PATH_LIBMESH_NEW(0.9.0,yes)

AC_ARG_ENABLE([libmesh-flags],
               AC_HELP_STRING([--enable-libmesh-flags],[enable use of libmesh cxxflags for CXXFLAGS]),
               USE_LIBMESH_CXX_FLAGS=1
               AC_DEFINE(USE_LIBMESH_CXX_FLAGS,1,[Define if libmesh cxxflags enabled]),[])


if test "$USE_LIBMESH_CXX_FLAGS" = "1"; then
   CXXFLAGS=`$LIBMESH_CONFIG --cppflags --cxxflags`
else
 dnl TODO: WE NEED TO DISTINGUISH BETWEEN COMPILERS IN A BETTER WAY
   CXXFLAGS="-O0 -g -fno-unsafe-math-optimizations"
fi

dnl -Wall warnings, -Wall the time.
dnl AX_CXXFLAGS_WARN_ALL

#--------------------------------
# Query configuration environment
#------------------------------
AX_SUMMARIZE_ENV

dnl------------------------------------------
dnl            Doxygen support
dnl------------------------------------------

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN(AGNOS, doxygen/agnos.dox, doc/doxygen)


dnl-----------------------------------------------
dnl Generate Makefiles
dnl-----------------------------------------------
AC_CONFIG_FILES([
                 Makefile 
                 src/Makefile 
                 test/Makefile 
                 examples/Makefile 
                 examples/catenary/Makefile
                 doxygen/Makefile])


#-------------------
# Generate makefiles
#-------------------
AC_OUTPUT()


#------------------------------
# Final summary
#------------------------------
AX_SUMMARIZE_CONFIG

