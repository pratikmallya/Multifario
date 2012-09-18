#*
#*  PROGRAM NAME:  multifario
#*
#*  (c) COPYRIGHT INTERNATIONAL BUSINESS MACHINES
#*  CORPORATION 12/1/2001.  ALL RIGHTS RESERVED.
#*
#*  Please refer to the LICENSE file in the top directory
#*
#*      author: Mike Henderson mhender@watson.ibm.com
#*
echo "@(#) $Id: configure.ac,v 1.5 2010/12/07 16:55:15 mhender Exp $"

AC_INIT(multifario,v0.9i,[mhender@us.ibm.com],multifario)
AC_PREREQ([2.0])
AC_PREFIX_DEFAULT(.)
AC_CONFIG_HEADER(include/multifarioConfig.h)
AC_CONFIG_FILES([Makefile])
#AC_CONFIG_SUBDIRS([objs Doc])
AC_CONFIG_FILES([objs/Makefile])

AC_PROG_LIBTOOL
AM_INIT_AUTOMAKE
LT_INIT
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_LIBOBJ_DIR([objs])

AC_PROG_CC
AC_SUBST(INSTDIR)
AC_PREFIX_DEFAULT($(PWD))

AC_PROG_F77
AC_PROG_FC
AC_PROG_INSTALL
AC_CHECK_HEADERS(getopt.h)

AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([gettimeofday])
AC_CHECK_FUNCS([pow])
AC_CHECK_FUNCS([sqrt])
AC_CHECK_FUNCS([strerror])
AC_CHECK_FUNCS([strrchr])
AC_CHECK_FUNCS([strstr])
AC_CHECK_HEADERS([float.h])
AC_CHECK_HEADERS([sys/time.h])
AC_C_CONST
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_VPRINTF
AC_HEADER_TIME
AC_STRUCT_TM
AC_TYPE_SIZE_T
#--------------------------------------------------------------------------------------

# Fortran libraries.

AC_FC_LIBRARY_LDFLAGS()
AC_FC_DUMMY_MAIN()
AC_FC_WRAPPERS

AC_F77_LIBRARY_LDFLAGS()
AC_F77_WRAPPERS

#--------------------------------------------------------------------------------------

AC_ARG_WITH(libs,
      [AS_HELP_STRING([--with-libs=<"-Ldirpath1 -Ldirpath2 ..."],[Directories to search for libraries])],
      [LIBS="$LIBS ${with_libs}"],[])

#--------------------------------------------------------------------------------------

AC_ARG_WITH(incs,
      [AS_HELP_STRING([--with-incs=<"-Idirpath1 -Idirpath2 ..."],[Directories to search for includes])],
      [INCS="$INCS ${with_incs}"],[])

#--------------------------------------------------------------------------------------

# libTiff.

AC_LANG(C)

echo "--------------------------------------------------------------------------------------"
echo "     testing for libtiff    "

AC_ARG_WITH(libtiff,
      [AS_HELP_STRING([--with-libtiff=<path>],[The directory in which lib/libtiff.a can be found])],
      [LIBS="$LIBS -L${with_libtiff}/lib -lm";INCS="$INCS -I${with_libtiff}/include"],[])

AC_CHECK_LIB(tiff,TIFFSetField,[HAVE_LIBTIFF=yes;LIBS="$LIBS -ltiff"],[HAVE_LIBTIFF=no])
echo "     HAVE_LIBTIFF="$HAVE_LIBTIFF
echo "     LIBS        ="$LIBS

AH_TEMPLATE([HAVE_LIBTIFF],[Set if the library libtiff is available])
AS_IF([ test x$HAVE_LIBTIFF = xyes ],[AC_DEFINE(HAVE_LIBTIFF)])
#--------------------------------------------------------------------------------------

# fitpack

echo "--------------------------------------------------------------------------------------"
echo "     testing for fitpack    "

AC_ARG_WITH(fitpack,
      [AS_HELP_STRING([--with-fitpack=<path>],[The directory in which lib/libFITPACK.a can be found])],
      [LIBS="$LIBS -L${with_fitpack}/lib -lm";INCS="$INCS -I${with_fitpack}/include"],
      [])

AC_LANG(C)

AC_CHECK_LIB(FITPACK,curv1,[HAVE_FITPACK=yes;LIBS="$LIBS -lFITPACK"],[HAVE_FITPACK=no])

AS_IF([ test x$HAVE_FITPACK = xno ], 
  [
    AC_CHECK_LIB(FITPACK,curv1_, [HAVE_FITPACK=yes;LIBS="$LIBS -lFITPACK"] , [HAVE_FITPACK=no] )
  ])

AS_IF(
  [ test x$HAVE_FITPACK = xno ], 
  [
    AC_CHECK_LIB(FITPACK,CURV1, [HAVE_FITPACK=yes;LIBS="$LIBS -lFITPACK"] , [HAVE_FITPACK=no] )
  ])

AS_IF(
  [ test x$HAVE_FITPACK = xno ], 
  [
     AC_CHECK_LIB(FITPACK,CURV1_,[HAVE_FITPACK=yes;LIBS="$LIBS -lFITPACK"] , [HAVE_FITPACK=no] )
  ])

echo "     HAVE_FITPACK="$HAVE_FITPACK
echo "     LIBS        ="$LIBS

AM_CONDITIONAL([HAVEFITPACK],[test x$HAVE_FITPACK == xyes ])
AH_TEMPLATE([HAVE_FITPACK],[Set if the library libFITPACK is available])
AS_IF([ test x$HAVE_FITPACK = xyes ],[AC_DEFINE(HAVE_FITPACK)])

AC_LANG(Fortran)
fitpack_flavor="unset"
AS_IF(
      [ test x$HAVE_FITPACK = xyes ] , 
      [echo "Look for curv1";AC_CHECK_FUNCS([curv1], [echo "Found"; fitpack_flavor="fc"], [echo "Not found"; fitpack_flavor="f77"])],
      [echo "Don't have Fitpack";fitpack_flavor="undetermined"])
echo "     fitpack_flavor="$fitpack_flavor

#--------------------------------------------------------------------------------------

# dierckx

echo "--------------------------------------------------------------------------------------"
echo "     testing for dierckx    "

AC_ARG_WITH(dierckx,
      [AS_HELP_STRING([--with-dierckx=<path>],[The directory in which lib/libDierckx.a can be found])],
      [LIBS="$LIBS -L${with_dierckx}/lib -lm";INCS="$INCS -I${with_dierckx}/include"],
      [])

AC_LANG(C)

AC_CHECK_LIB(Dierckx,percur,[HAVE_DIECKXx=yes;LIBS="$LIBS -lDierckx"],[HAVE_DIERCKX=no])

AS_IF([ test x$HAVE_DIERCKX = xno ], 
  [
    AC_CHECK_LIB(Dierckx,percur_, [HAVE_DIERCKX=yes;LIBS="$LIBS -lDierckx"] , [HAVE_DIERCKX=no] )
  ])

AS_IF(
  [ test x$HAVE_DIERCKX = xno ], 
  [
    AC_CHECK_LIB(Dierckx,PERCUR, [HAVE_DIERCKX=yes;LIBS="$LIBS -lDierckx"] , [HAVE_DIERCKX=no] )
  ])

AS_IF(
  [ test x$HAVE_DIERCKX = xno ], 
  [
     AC_CHECK_LIB(Dierckx,PERCUR_,[HAVE_DIERCKX=yes;LIBS="$LIBS -lDierckx"] , [HAVE_DIERCKX=no] )
  ])

echo "     HAVE_DIERCKX="$HAVE_DIERCKX
echo "     LIBS        ="$LIBS

AM_CONDITIONAL([HAVEDIERCKX],[test x$HAVE_DIERCKX = xyes ])
AH_TEMPLATE([HAVE_DIERCKX],[Set if the library libDierckx is available])
AS_IF([ test x$HAVE_DIERCKX = xyes ],[AC_DEFINE(HAVE_DIERCKX)])

AC_LANG(Fortran)
dierckx_flavor="unset"
AS_IF(
      [ test x$HAVE_DIERCKX = xyes ] , 
      [echo "Look for percur";AC_CHECK_FUNCS([percur], [echo "Found";dierckx_flavor="fc"], [echo "Not found";dierckx_flavor="f77"])],
      [echo "Don't have Dierckx";dierckx_flavor="undetermined"])

echo "     dierckx_flavor="$dierckx_flavor

#--------------------------------------------------------------------------------------

# Blas.

echo "--------------------------------------------------------------------------------------"
echo "     testing for blas    "

AC_LANG(C)
AX_BLAS([HAVE_BLAS=yes],[HAVE_BLAS=no])
echo "     HAVE_BLAS="$HAVE_BLAS

AH_TEMPLATE([HAVE_BLAS],[Set if the library libblas is available])
AS_IF([ test x$HAVE_BLAS = xyes ],[AC_DEFINE(HAVE_BLAS)])

AC_LANG(Fortran)
blas_flavor="none"
AS_IF([ test x$HAVE_BLAS = xyes ], [AC_CHECK_FUNCS([daxpy],[blas_flavor="fc"],[blas_flavor="f77"])],[])
echo "     blas_flavor="$blas_flavor

AH_TEMPLATE([NEED_BLAS_WRAPPERS],[Set if the fortran compilers for automake and blas are different, in
                                  this case the link will fail unless wrappers with the automake fortran's
                                  mangled names are added to call the blas mangled names])

AS_IF([ test x$blas_flavor = xfc ],[ echo "Setting NEED_BLAS_WRAPPERS";AC_DEFINE(NEED_BLAS_WRAPPERS)],[])

#--------------------------------------------------------------------------------------

# Lapack.

echo "--------------------------------------------------------------------------------------"
echo "     testing for lapack    "

AC_LANG(C)
AX_LAPACK([HAVE_LAPACK=yes],[HAVE_LAPACK=no])
AM_CONDITIONAL([HAVELAPACK],[test x$HAVE_LAPACK = xyes])
echo "     HAVE_LAPACK="$HAVE_LAPACK

AH_TEMPLATE([HAVE_LAPACK],[Set if the library liblapack is available])
AS_IF([ test x$HAVE_LAPACK = xyes ],[AC_DEFINE(HAVE_LAPACK)])

AC_LANG(Fortran)
AS_IF(
      [ test x$HAVE_LAPACK = xyes ] , 
      [AC_CHECK_FUNCS([dgeev],
            [lapack_flavor="fc"],
            [lapack_flavor="f77"])
      ],[])
echo "     lapack_flavor="$lapack_flavor


#--------------------------------------------------------------------------------------

# Doxygen

echo "--------------------------------------------------------------------------------------"
echo "     testing for doxygen    "

AC_ARG_ENABLE([doxygen],
            [  --enable-doxygen              enable DOXYGEN (default is true)],
            [ case "${enableval}" in
              yes) doxygen=true ;;
              no)  doxygen=false ;;
              *) AC_MSG_ERROR([bad value ${enableval} for --enable-doxygen]) ;;
      esac], [doxygen=true])
AM_CONDITIONAL([DOXYGEN], [test x$doxygen = xtrue])
if test x$doxygen = xtrue
then
m4_include([m4/ax_prog_doxygen.m4])
DX_INIT_DOXYGEN(multifario, multifario.dox, localDocs)
DX_DOXYGEN_FEATURE(ON)
#DX_DOT_FEATURE(ON)
DX_HTML_FEATURE(ON)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
fi

#--------------------------------------------------------------------------------------

echo "--------------------------------------------------------------------------------------"
echo "     setting up lapack names for calls from C    "

AH_TEMPLATE([CALLDGEEV ],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLDGESVD],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLDGETRF],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLDGETRS],[Set to the correct name to use to call.])

AC_DEFINE(CALLDGEEV,[dggeev_not_found],[unmangled name])
AC_DEFINE(CALLDGESVD,[dggesvd_not_found],[unmangled name])
AC_DEFINE(CALLDGETRF,[dggetrf_not_found],[unmangled name])
AC_DEFINE(CALLDGETRS,[dggetrs_not_found],[unmangled name])

AS_IF([ test  x$HAVE_LAPACK = xyes -a x$lapack_flavor = xf77 ] ,
  [

   echo "defining Lapack calls as F77."

   AC_F77_FUNC(dgeev, CALLDGEEV)
   AC_F77_FUNC(dgesvd, CALLDGESVD)
   AC_F77_FUNC(dgetrf, CALLDGETRF)
   AC_F77_FUNC(dgetrs, CALLDGETRS)

   echo "CALLDGEEV="$CALLDGEEV
   echo "CALLDGESVD="$CALLDGESVD
   echo "CALLDGETRF="$CALLDGETRF
   echo "CALLDGETRS="$CALLDGETRS

   AC_DEFINE_UNQUOTED(CALLDGEEV,$CALLDGEEV,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLDGESVD,$CALLDGESVD,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLDGETRF,$CALLDGETRF,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLDGETRS,$CALLDGETRS,[unmangled name])
  ],
  [])

AS_IF([ test x$HAVE_LAPACK = xyes -a x$lapack_flavor = xfc ] ,
  [

   echo "defining Lapack calls as FC."

   AC_FC_FUNC(dgeev, CALLDGEEV)
   AC_FC_FUNC(dgesvd, CALLDGESVD)
   AC_FC_FUNC(dgetrf, CALLDGETRF)
   AC_FC_FUNC(dgetrs, CALLDGETRS)

   echo "CALLDGEEV="$CALLDGEEV
   echo "CALLDGESVD="$CALLDGESVD
   echo "CALLDGETRF="$CALLDGETRF
   echo "CALLDGETRS="$CALLDGETRS
 
   AC_DEFINE_UNQUOTED(CALLDGEEV,$CALLDGEEV,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLDGESVD,$CALLDGESVD,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLDGETRF,$CALLDGETRF,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLDGETRS,$CALLDGETRS,[unmangled name])
  ],
  [])

echo "--------------------------------------------------------------------------------------"
echo "     setting up Bordered names for calls from C    "

AH_TEMPLATE([CALLDBOSL ],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLDBONV ],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLDBOFA ],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLDBOSVD],[Set to the correct name to use to call.])

   echo "defining Bordered calls as F77."

AC_F77_FUNC(dbosl,CALLDBOSL)
AC_F77_FUNC(dbonv,CALLDBONV)
AC_F77_FUNC(dbofa,CALLDBOFA)
AC_F77_FUNC(dbosvd,CALLDBOSVD)

echo "CALLDBOSL="$CALLDBOSL
echo "CALLDBONV="$CALLDBONV
echo "CALLDBOFA="$CALLDBOFA
echo "CALLDBOSVD="$CALLDBOSVD

AC_DEFINE_UNQUOTED(CALLDBOSL,$CALLDBOSL,[unmangled name])
AC_DEFINE_UNQUOTED(CALLDBONV,$CALLDBONV,[unmangled name])
AC_DEFINE_UNQUOTED(CALLDBOFA,$CALLDBOFA,[unmangled name])
AC_DEFINE_UNQUOTED(CALLDBOSVD,$CALLDBOSVD,[unmangled name])

echo "--------------------------------------------------------------------------------------"
echo "     setting up FITPACK names for calls from C    "

AH_TEMPLATE([CALLCURV1 ],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLCURVP1],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLFCURV2 ],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLFCURVP2],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLFDCURVP2],[Set to the correct name to use to call.])

AC_DEFINE(CALLCURV1,"curv1_not_found",[unmangled name])
AC_DEFINE(CALLCURVP1,"curvp1_not_found",[unmangled name])
AC_DEFINE(CALLFCURV2,"fcurv2_not_found",[unmangled name])
AC_DEFINE(CALLFCURVP2,"fcurvp2_not_found",[unmangled name])
AC_DEFINE(CALLFDCURVP2,"fdcurvp2_not_found",[unmangled name])

   echo "HAVE_FITPACK="$HAVE_FITPACK
   echo "fitpack_flavor="$fitpack_flavor

AS_IF([ test x$HAVE_FITPACK = xyes -a x$fitpack_flavor = xf77 ] ,
  [

   echo "defining Fitpack calls as F77."

   AC_F77_FUNC(curv1,CALLCURV1)
   AC_F77_FUNC(curvp1,CALLCURVP1)
   AC_F77_FUNC(fcurv2,CALLFCURV2)
   AC_F77_FUNC(fcurvp2,CALLFCURVP2)
   AC_F77_FUNC(fdcurvp2,CALLFDCURVP2)

   echo "CALLCURV1   ="$CALLCURV1
   echo "CALLCURVP1  ="$CALLCURVP1
   echo "CALLFCURVP2 ="$CALLFCURVP2
   echo "CALLFCURVP2 ="$CALLFCURVP2
   echo "CALLFDCURVP2="$CALLFDCURVP2

   AC_DEFINE_UNQUOTED(CALLCURV1,$CALLCURV1,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLCURVP1,$CALLCURVP1,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLFCURV2,$CALLFCURV2,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLFCURVP2,$CALLFCURVP2,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLFDCURVP2,$CALLFDCURVP2,[unmangled name])
  ],[])

   echo "HAVE_FITPACK="$HAVE_FITPACK
   echo "fitpack_flavor="$fitpack_flavor

AS_IF([ test x$HAVE_FITPACK = xyes -a x$fitpack_flavor = xfc ] ,
  [

   echo "defining Fitpack calls as FC."

   AC_FC_FUNC(curv1,CALLCURV1)
   AC_FC_FUNC(curvp1,CALLCURVP1)
   AC_FC_FUNC(fcurv2,CALLFCURV2)
   AC_FC_FUNC(fcurvp2,CALLFCURVP2)
   AC_FC_FUNC(fdcurvp2,CALLFDCURVP2)

   echo "CALLCURV1   ="$CALLCURV1
   echo "CALLCURVP1  ="$CALLCURV1
   echo "CALLFCURV2="$CALLFCURV2
   echo "CALLFCURVP2="$CALLFCURVP2
   echo "CALLFDCURVP2="$CALLFDCURVP2

   AC_DEFINE_UNQUOTED(CALLCURV1,$CALLCURV1,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLCURVP1,$CALLCURVP1,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLFCURV2,$CALLFCURV2,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLFCURVP2,$CALLFCURVP2,[unmangled name])
   AC_DEFINE_UNQUOTED(CALLFDCURVP2,$CALLFDCURVP2,[unmangled name])
  ],
  [])

echo "--------------------------------------------------------------------------------------"
echo "     setting up Dierckx names for calls from C    "

AH_TEMPLATE([CALLPERCUR ],[Set to the correct name to use to call.])
AH_TEMPLATE([CALLSPLEV  ],[Set to the correct name to use to call.])

AC_DEFINE(CALLPERCUR,"percur_not_found",[unmangled name])
AC_DEFINE(CALLSPLEV,"splev_not_found",[unmangled name])

AS_IF([ test x$HAVE_DIERCKX = xyes -a x$dierckx_flavor = xf77 ] ,
  [

   echo "defining dierckx calls as F77."

   AC_F77_FUNC(percur,CALLPERCUR)
   AC_F77_FUNC(splev,CALLSPLEV)

   echo "CALLPERCUR="$CALLPERCUR
   echo "CALLSPLEV="$CALLSPLEV

   AC_DEFINE_UNQUOTED(CALLPERCUR,$CALLPERCUR,[])
   AC_DEFINE_UNQUOTED(CALLSPLEV,$CALLSPLEV,[])
  ],
  [])

AS_IF([ test x$HAVE_DIERCKX = xyes -a x$dierckx_flavor = xfc ] ,
  [

   echo "defining dierckx calls as FC."

   AC_FC_FUNC(percur,CALLPERCUR)
   AC_FC_FUNC(splev,CALLSPLEV)

   echo "CALLPERCUR="$CALLPERCUR
   echo "CALLSPLEV="$CALLSPLEV

   AC_DEFINE_UNQUOTED(CALLPERCUR,$CALLPERCUR,[])
   AC_DEFINE_UNQUOTED(CALLSPLEV,$CALLSPLEV,[])
  ],
  [])

echo "--------------------------------------------------------------------------------------"

AC_OUTPUT
