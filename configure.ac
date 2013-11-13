#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(agnos, 0.1.0, cbryant@ices.utexas.edu)
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADER(agnos_config.h.tmp)
AX_PREFIX_CONFIG_H(agnos_config.h,AGNOS,agnos_config.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([subdir-objects dist-xz dist-bzip2 tar-ustar color-tests 1.11])
# AM_INIT_AUTOMAKE([dist-xz dist-bzip2 tar-ustar color-tests 1.11])
AM_SILENT_RULES(yes)  # use silent rules where available - automake 1.11

LT_INIT

dnl Release versioning
AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=$AX_MAJOR_VERSION
GENERIC_MINOR_VERSION=$AX_MINOR_VERSION
GENERIC_MICRO_VERSION=$AX_POINT_VERSION
GENERIC_LIBRARY_NAME=agnos

dnl shared library versioning

PACKAGE=$GENERIC_LIBRARY_NAME
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

dnl AC_SUBST(GENERIC_LIBRARY_VERSION)
dnl AC_SUBST(GENERIC_LIBRARY_NAME)

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

AC_DEFINE_UNQUOTED([MAJOR_VERSION],$GENERIC_MAJOR_VERSION,[Major version])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$GENERIC_MINOR_VERSION,[Minor version])
AC_DEFINE_UNQUOTED([MICRO_VERSION],$GENERIC_MICRO_VERSION,[Micro version])


#------------------------------
# Checks for compilers and MPI
#------------------------------
AC_PROG_CC
AC_LANG([C])
ACX_MPI([CC="$MPICC"], [AC_MSG_ERROR([Could not find MPI.])])

AC_PROG_CXX
AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], [AC_MSG_ERROR([Could not find MPI.])])

AC_CHECK_LIB(mpich, MPI_Init, [USING_MPICH=1
                             AC_DEFINE(USING_MPICH,1,[Define if we are using
                                                      mpich])])
#---------------------------------------------------------
# Checks for library prerequisites for other libraries...
#---------------------------------------------------------

#----------------
# Check for boost
#----------------
BOOST_REQUIRE([1.36])
BOOST_MATH()
# BOOST_PROGRAM_OPTIONS()
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

CXXFLAGS+=" -std=c++0x"

dnl -Wall warnings, -Wall the time.
dnl AX_CXXFLAGS_WARN_ALL


# -------------------------------------------------------------
# cppunit C++ unit testing -- enabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(cppunit,
             AC_HELP_STRING([--enable-cppunit],
                            [Build with cppunit C++ unit testing support]),
		[case "${enableval}" in
		  yes)  enablecppunit=yes ;;
		   no)  enablecppunit=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-cppunit) ;;
		 esac],
		 [enablecppunit=$enableoptional])
if (test "$enablecppunit" = yes) ; then
  AM_PATH_CPPUNIT([1.10.0],[enablecppunit=yes],[enablecppunit=no])
fi
AM_CONDITIONAL(LIBMESH_ENABLE_CPPUNIT, test x$enablecppunit = xyes)
# -------------------------------------------------------------

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

DX_INIT_DOXYGEN(AGNOS, doxygen/agnos.dox, doc/)
# AC_CONFIG_FILES(doxygen/txt_common/about.page)

dnl-----------------------------------------------
dnl Generate files for unit and regression testing
dnl-----------------------------------------------
AC_CONFIG_FILES(test/polynomialTest.sh,       [chmod +x test/polynomialTest.sh])
AC_CONFIG_FILES(test/evaluationTest.sh,       [chmod +x test/evaluationTest.sh])
AC_CONFIG_FILES(test/catenaryTest.sh,	        [chmod +x test/catenaryTest.sh])
AC_CONFIG_FILES(test/catenaryLibmeshTest.sh,	[chmod +x test/catenaryLibmeshTest.sh])
AC_CONFIG_FILES(test/parallelCatenary.sh,	    [chmod +x test/parallelCatenary.sh])
AC_CONFIG_FILES(test/viscousBurgersTest.sh,	  [chmod +x test/viscousBurgersTest.sh])
AC_CONFIG_FILES(test/elementTest.sh,	        [chmod +x test/elementTest.sh])

dnl-----------------------------------------------
dnl Generate run scripts for examples
dnl-----------------------------------------------
AC_CONFIG_FILES(examples/catenary/run.sh,	      [chmod +x examples/catenary/run.sh])
AC_CONFIG_FILES(examples/catenaryLibmesh/run.sh,[chmod +x examples/catenaryLibmesh/run.sh])
AC_CONFIG_FILES(examples/higherOrderDiffusion/run.sh,[chmod +x examples/higherOrderDiffusion/run.sh])
AC_CONFIG_FILES(examples/viscousBurgers/run.sh,	[chmod +x examples/viscousBurgers/run.sh])

dnl-----------------------------------------------
dnl Generate Makefiles
dnl-----------------------------------------------
AC_CONFIG_FILES([
                 Makefile 
                 src/Makefile 
                 test/Makefile 
                 examples/Makefile 
                 doxygen/Makefile])


#-------------------
# Generate makefiles
#-------------------
AC_OUTPUT()


#------------------------------
# Final summary
#------------------------------
AX_SUMMARIZE_CONFIG

