/* include/multifarioConfig.h.  Generated from multifarioConfig.h.in by configure.  */
/* include/multifarioConfig.h.in.  Generated from configure.ac by autoheader.  */

/* unmangled name */
#define CALLCURV1 curv1_

/* unmangled name */
#define CALLCURVP1 curvp1_

/* unmangled name */
#define CALLDBOFA dbofa_

/* unmangled name */
#define CALLDBONV dbonv_

/* unmangled name */
#define CALLDBOSL dbosl_

/* unmangled name */
#define CALLDBOSVD dbosvd_

/* unmangled name */
#define CALLDGEEV dgeev_

/* unmangled name */
#define CALLDGESVD dgesvd_

/* unmangled name */
#define CALLDGETRF dgetrf_

/* unmangled name */
#define CALLDGETRS dgetrs_

/* unmangled name */
#define CALLFCURV2 fcurv2_

/* unmangled name */
#define CALLFCURVP2 fcurvp2_

/* unmangled name */
#define CALLFDCURVP2 fdcurvp2_

/* unmangled name */
#define CALLPERCUR percur_

/* unmangled name */
#define CALLSPLEV splev_

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## __

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## __

/* Set if the library libblas is available */
#define HAVE_BLAS 1

/* Define to 1 if you have the `curv1' function. */
#define HAVE_CURV1 1

/* Define to 1 if you have the `daxpy' function. */
/* #undef HAVE_DAXPY */

/* Define to 1 if you have the `dgeev' function. */
/* #undef HAVE_DGEEV */

/* Set if the library libDierckx is available */
#define HAVE_DIERCKX 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you don't have `vprintf' but do have `_doprnt.' */
/* #undef HAVE_DOPRNT */

/* Set if the library libFITPACK is available */
#define HAVE_FITPACK 1

/* Define to 1 if you have the <float.h> header file. */
#define HAVE_FLOAT_H 1

/* Define to 1 if you have the `floor' function. */
#define HAVE_FLOOR 1

/* Define to 1 if you have the <getopt.h> header file. */
#define HAVE_GETOPT_H 1

/* Define to 1 if you have the `gettimeofday' function. */
#define HAVE_GETTIMEOFDAY 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Set if the library liblapack is available */
#define HAVE_LAPACK 1

/* Set if the library libtiff is available */
#define HAVE_LIBTIFF 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `percur' function. */
#define HAVE_PERCUR 1

/* Define to 1 if you have the `pow' function. */
#define HAVE_POW 1

/* Define to 1 if your system has a GNU libc compatible `realloc' function,
   and to 0 otherwise. */
#define HAVE_REALLOC 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strerror' function. */
#define HAVE_STRERROR 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strrchr' function. */
#define HAVE_STRRCHR 1

/* Define to 1 if you have the `strstr' function. */
#define HAVE_STRSTR 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vprintf' function. */
#define HAVE_VPRINTF 1

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#define LT_OBJDIR ".libs/"

/* Set if the fortran compilers for automake and blas are different, in this
   case the link will fail unless wrappers with the automake fortran's mangled
   names are added to call the blas mangled names */
/* #undef NEED_BLAS_WRAPPERS */

/* Name of package */
#define PACKAGE "multifario"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "mhender@us.ibm.com"

/* Define to the full name of this package. */
#define PACKAGE_NAME "multifario"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "multifario v0.9i"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "multifario"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "v0.9i"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Define to 1 if your <sys/time.h> declares `struct tm'. */
/* #undef TM_IN_SYS_TIME */

/* Version number of package */
#define VERSION "v0.9i"

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to rpl_realloc if the replacement function should be used. */
/* #undef realloc */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
