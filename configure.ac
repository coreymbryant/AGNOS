#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(agnos, 0.1.0, cbryant@ices.utexas.edu)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_AUX_DIR([build-aux])
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

DX_INIT_DOXYGEN(AGNOS, doxygen/adnos.dox, doc/doxygen)


dnl-----------------------------------------------
dnl Generate Makefiles
dnl-----------------------------------------------
AC_CONFIG_FILES([Makefile src/Makefile test/Makefile doxygen/Makefile])


#-------------------
# Generate makefiles
#-------------------
AC_OUTPUT()


#------------------------------
# Final summary
#------------------------------
AX_SUMMARIZE_CONFIG

