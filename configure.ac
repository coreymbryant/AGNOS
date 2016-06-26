#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([agnos], [0.1.0], [cbryant@ices.utexas.edu])

# AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADER(agnos_config.h.tmp)
AX_PREFIX_CONFIG_H(agnos_config.h,AGNOS,agnos_config.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
# AM_INIT_AUTOMAKE
AM_INIT_AUTOMAKE([color-tests])
# AM_INIT_AUTOMAKE([dist-xz dist-bzip2 tar-ustar color-tests 1.11])
AM_SILENT_RULES(yes)  # use silent rules where available - automake 1.11

# LT_INIT

dnl Release versioning
AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=1
GENERIC_MICRO_VERSION=0
GENERIC_LIBRARY_NAME=agnos


dnl shared library versioning
# PACKAGE=$GENERIC_LIBRARY_NAME
# GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

dnl AC_SUBST(GENERIC_LIBRARY_VERSION)
dnl AC_SUBST(GENERIC_LIBRARY_NAME)

# AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

# AC_DEFINE_UNQUOTED([MAJOR_VERSION],$GENERIC_MAJOR_VERSION,[Major version])
# AC_DEFINE_UNQUOTED([MINOR_VERSION],$GENERIC_MINOR_VERSION,[Minor version])
# AC_DEFINE_UNQUOTED([MICRO_VERSION],$GENERIC_MICRO_VERSION,[Micro version])

# Required for some macros
AX_AM_MACROS

#------------------------------
# Checks for compilers and MPI
#------------------------------
AC_PROG_CC
AC_PROG_CXX

AC_LANG([C])
ACX_MPI([CC="$MPICC"], [AC_MSG_ERROR([Could not find MPI.])])
AC_LANG([C++])
ACX_MPI([CXX="$MPICXX"], [AC_MSG_ERROR([Could not find MPI.])])

AC_CHECK_LIB(mpich, MPI_Init, [USING_MPICH=1
                             AC_DEFINE(USING_MPICH,1,[Define if we are using
                                                      mpich])])
#---------------------------------------------------------
# Checks for library prerequisites for other libraries...
#---------------------------------------------------------


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

# CXXFLAGS+=" -std=c++0x"
dnl--------------------------
dnl Even if no dependency needed it,
dnl C++11 makes our asserts nicer,
dnl so let's try and enable it.
dnl--------------------------
AC_ARG_ENABLE(cxx11,
              AC_HELP_STRING([--enable-cxx11],
                             [build with C++11 support]),
              [case "${enableval}" in
                yes)  enablecxx11=yes ;;
                 no)  enablecxx11=no ;;
                  *)  AC_MSG_ERROR(bad value ${enableval} for --enable-cxx11) ;;
               esac],
               [enablecxx11=optional])

if (test x$enablecxx11 = xyes); then
  AX_CXX_COMPILE_STDCXX_11(noext, mandatory)
elif (test x$enablecxx11 = xoptional); then
  AX_CXX_COMPILE_STDCXX_11(noext, optional)
fi

dnl -Wall warnings, -Wall the time.
dnl AX_CXXFLAGS_WARN_ALL

# Check for available programs.
AC_PROG_RANLIB
AC_PROG_CXX
AC_DISABLE_STATIC
AC_PROG_LIBTOOL
AC_LANG([C++])



dnl --------------------------------------------------------------
dnl Get the source tree directory
dnl --------------------------------------------------------------
PWD=`pwd`
AC_SUBST(PWD)

#-------------------------
# External Library Checks
#-------------------------

#----------------
# Check for boost
#----------------
AC_LANG([C++])
BOOST_REQUIRE([1.36])
BOOST_MATH()
BOOST_PROGRAM_OPTIONS()
BOOST_TEST()
AC_CACHE_SAVE

# -------------------------------------------------------------
# Check for gsl
# -------------------------------------------------------------
AM_PATH_GSL([1.15],
            [],
            [AC_MSG_ERROR(GSL version is not new enough)]
            )

# -------------------------------------------------------------
# Check for hdf5
# -------------------------------------------------------------
AX_PATH_HDF5_NEW([1.8.0],[yes])

# -------------------------------------------------------------
# cppunit C++ unit testing -- disabled by default
# -------------------------------------------------------------
AC_ARG_ENABLE(cppunit,
             AC_HELP_STRING([--enable-cppunit],
                            [Build with cppunit C++ unit testing support]),
		[case "${enableval}" in
		  yes)  enablecppunit=yes ;;
		   no)  enablecppunit=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-cppunit) ;;
		 esac],
		 [enablecppunit=no])
if (test "$enablecppunit" = yes) ; then
  AM_PATH_CPPUNIT([1.10.0],[enablecppunit=yes],[enablecppunit=no])
fi
AM_CONDITIONAL(ENABLE_CPPUNIT, test x$enablecppunit = xyes)
# -------------------------------------------------------------


#--------------------------------
# Check for channelFlow library
#--------------------------------
AC_ARG_ENABLE(channelflow,
             AC_HELP_STRING([--enable-channelflow],
                            [Build with channelflow library for RANS simulations]),
		[case "${enableval}" in
		  yes)  enablechannelflow=yes ;;
		   no)  enablechannelflow=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-channelflow) ;;
		 esac],
		 [enablechannelflow=no])
if (test "$enablechannelflow" = yes) ; then
  AC_DEFINE_UNQUOTED([CHANNELFLOW_DIR],$CHANNELFLOW_DIR,[channelflow dir
                      location])
  AC_SUBST(CHANNELFLOW_DIR)
  AC_DEFINE(ENABLE_CHANNELFLOW,1,[Define if we are using channelflow])
fi
AM_CONDITIONAL(ENABLE_CHANNELFLOW, test x$enablechannelflow = xyes)

#--------------------------------
# Check for grins library
#--------------------------------
AC_ARG_ENABLE(grins,
             AC_HELP_STRING([--enable-grins],
                            [Build with grins library for NS simulations]),
		[case "${enableval}" in
		  yes)  enablegrins=yes ;;
		   no)  enablegrins=no ;;
 		    *)  AC_MSG_ERROR(bad value ${enableval} for --enable-grins) ;;
		 esac],
		 [enablegrins=no])
if (test "$enablegrins" = yes) ; then
  AC_DEFINE_UNQUOTED([GRINS_DIR],$GRINS_DIR,[grins dir location])
  AC_SUBST(GRINS_DIR)
  AC_DEFINE(ENABLE_GRINS,1,[Define if we are using grins])
fi
AM_CONDITIONAL(ENABLE_GRINS, test x$enablegrins = xyes)


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
AC_CONFIG_FILES(test/ChannelTest.sh,	        [chmod +x test/ChannelTest.sh])
AC_CONFIG_FILES(test/GrinsTest.sh,	          [chmod +x test/GrinsTest.sh])
AC_CONFIG_FILES(test/LibmeshTest.sh,	        [chmod +x test/LibmeshTest.sh])
AC_CONFIG_FILES(test/hdf5Test.sh,	            [chmod +x test/hdf5Test.sh])
AC_CONFIG_FILES(test/EvaluatorTest.sh,	      [chmod +x test/EvaluatorTest.sh])
AC_CONFIG_FILES(test/DriverTest.sh,	          [chmod +x test/DriverTest.sh])

dnl-----------------------------------------------
dnl Generate run scripts for examples
dnl-----------------------------------------------
AC_CONFIG_FILES(examples/catenary/run.sh,	      [chmod +x examples/catenary/run.sh])
AC_CONFIG_FILES(examples/catenaryLibmesh/run.sh,[chmod +x examples/catenaryLibmesh/run.sh])
AC_CONFIG_FILES(examples/diffusion/run.sh,[chmod +x examples/diffusion/run.sh])
AC_CONFIG_FILES(examples/viscousBurgers/run.sh,	[chmod +x examples/viscousBurgers/run.sh])
AC_CONFIG_FILES(examples/channelFlow/run.sh,	[chmod +x examples/channelFlow/run.sh])
AC_CONFIG_FILES(examples/grins/run.sh,	[chmod +x examples/grins/run.sh])

dnl-----------------------------------------------
dnl Generate Makefiles
dnl-----------------------------------------------
AC_CONFIG_FILES([
                 Makefile 
                 include/Makefile 
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

