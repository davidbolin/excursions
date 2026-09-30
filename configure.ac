AC_INIT
AC_CONFIG_SRCDIR([src/gsl_fix.h])

# From https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Configure-and-cleanup :
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

CXX=`"${R_HOME}/bin/R" CMD config CXX`
if test -z "$CXX"; then
  AC_MSG_ERROR([No C++ compiler is available])
fi
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_LANG(C++)

AC_PROG_CC()
AC_PROG_CXX()

# From https://raw.githubusercontent.com/RobinHankin/gsl/master/configure.ac
# Following lines kindly supplied by Dirk Eddelbuettel
## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
   # Use gsl-config for header and linker arguments
   GSL_CFLAGS=`${GSL_CONFIG} --cflags`
   GSL_LIBS=`${GSL_CONFIG} --libs`
   CFLAGS="$CFLAGS $GSL_CFLAGS"  # kindly supplied by Ray Brownrigg.
   CXXFLAGS="$CFLAGS $GSL_CFLAGS"
else
   AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

AC_MSG_CHECKING([if GSL version >= 2.1])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <string.h>
#include <gsl/gsl_version.h>
#include <stdio.h>
#include <stdlib.h>
int main() {
#ifdef GSL_VERSION
   int major, minor;
   char *gslv = GSL_VERSION;
   if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
     exit (1);
   }
   exit(major < 2);
   if (major == 2) {
     exit(minor < 1);
   } else {
     exit(0);
   }
#else
  exit(1);
#endif
}
]])],
[gsl_version_ok=yes],
[gsl_version_ok=no],
[gsl_version_ok=yes])
if test "${gsl_version_ok}" = no; then
  AC_MSG_ERROR([Need GSL version >= 2.1])
else
  AC_MSG_RESULT([yes])
fi

## OpenMP
##
## On Linux, R provides the OpenMP flags in SHLIB_OPENMP_CXXFLAGS. On macOS
## these are empty, since Apple's clang has OpenMP disabled, but the compiler
## supports it with -Xclang -fopenmp, and the CRAN build of R ships the OpenMP
## runtime libomp.dylib in R_HOME/lib, which packages are encouraged to use
## (https://mac.r-project.org/openmp/). Test if a small OpenMP program can be
## built and run with R's runtime, and otherwise build without OpenMP. The
## runtime must be R's own, and not for example one from Homebrew, so that
## all packages use the same runtime and binary packages do not depend on it.
## Set the environment variable EXCURSIONS_OPENMP=no to build without OpenMP.
## (Square brackets are avoided below, since autoconf uses them for quoting.)
OPENMP_CXXFLAGS='$(SHLIB_OPENMP_CXXFLAGS)'
OPENMP_LIBS='$(SHLIB_OPENMP_CXXFLAGS)'
if test "${EXCURSIONS_OPENMP}" = "no"; then
  OPENMP_CXXFLAGS=""
  OPENMP_LIBS=""
  AC_MSG_NOTICE([OpenMP disabled by EXCURSIONS_OPENMP=no])
elif test "`uname -s`" = "Darwin"; then
  AC_MSG_CHECKING([for OpenMP with the runtime of R])
  omp_dir=`mktemp -d 2>/dev/null || echo "./omp-test"`
  mkdir -p "${omp_dir}"
  cat > "${omp_dir}/test_omp.cpp" <<'OMPTEST'
#include <omp.h>
extern "C" void test_omp(int *n) {
  int count = 0, sum = 0;
#pragma omp parallel num_threads(2) reduction(+:count)
  count += 1;
#pragma omp parallel for num_threads(2) schedule(dynamic) reduction(+:sum)
  for (int i = 1; i <= 100; i++) sum += i;
  *n = (sum == 5050) ? count : -1;
}
OMPTEST
  test_openmp() {
    rm -f "${omp_dir}/test_omp.o" "${omp_dir}/test_omp.so"
    (cd "${omp_dir}" && PKG_CXXFLAGS="$1" PKG_LIBS="$2" \
      "${R_HOME}/bin/R" CMD SHLIB test_omp.cpp) >> config.log 2>&1 || return 1
    "${R_HOME}/bin/Rscript" --vanilla -e \
      "dyn.load('${omp_dir}/test_omp.so'); n <- .C('test_omp', n = 0L)\$n; quit(status = as.integer(n < 1))" \
      >> config.log 2>&1 || return 1
    ## The runtime must be the one in R_HOME/lib
    if command -v otool > /dev/null 2>&1; then
      omp_lib=`otool -L "${omp_dir}/test_omp.so" | awk '/libomp/ {print $1}'`
      r_lib=`cd "${R_HOME}/lib" && pwd -P`
      case "${omp_lib}" in
        "${R_HOME}/lib/libomp.dylib" | "${r_lib}/libomp.dylib" | "@rpath/libomp.dylib") ;;
        *) echo "OpenMP runtime ${omp_lib} is not the one of R" >> config.log; return 1 ;;
      esac
    fi
    return 0
  }
  brew_prefix=`brew --prefix 2>/dev/null`
  if test ! -e "${R_HOME}/lib/libomp.dylib"; then
    OPENMP_CXXFLAGS=""
    OPENMP_LIBS=""
    AC_MSG_RESULT([no, R does not provide libomp.dylib])
  elif test_openmp "-Xclang -fopenmp" "-lomp"; then
    OPENMP_CXXFLAGS="-Xclang -fopenmp"
    OPENMP_LIBS="-lomp"
    AC_MSG_RESULT([yes])
  elif test -n "${brew_prefix}" && test -e "${brew_prefix}/opt/libomp/include/omp.h" && \
       test_openmp "-I${brew_prefix}/opt/libomp/include -Xclang -fopenmp" "-lomp"; then
    OPENMP_CXXFLAGS="-I${brew_prefix}/opt/libomp/include -Xclang -fopenmp"
    OPENMP_LIBS="-lomp"
    AC_MSG_RESULT([yes, with the header from Homebrew])
  else
    OPENMP_CXXFLAGS=""
    OPENMP_LIBS=""
    AC_MSG_RESULT([no, building without OpenMP])
  fi
  rm -rf "${omp_dir}"
fi

# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
AC_SUBST(OPENMP_CXXFLAGS)
AC_SUBST(OPENMP_LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

