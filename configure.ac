# See semver.org for version info
#################################

# Semver for benchmarks is a bit silly, but let's define it as:
#  - major: break in behavior, not comparable in perf
#  - minor: added benchmarks, comparable
#  - patch: fixes

m4_define([VERSION_MAJOR], [0])
m4_define([VERSION_MINOR], [1])
m4_define([VERSION_PATCH], [0])
m4_define([VERSION_STRING], VERSION_MAJOR.VERSION_MINOR.VERSION_PATCH)

# Init build tools
##################

AC_INIT([nrm-benchmarks],[VERSION_STRING],[swann@anl.gov])
AC_CONFIG_SRCDIR([src/nrm-benchmarks.h])
AC_CONFIG_AUX_DIR([m4])
AC_CONFIG_MACRO_DIR([m4])
# automake should fail on any error
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects 1.12])

# Detect features
#################

AC_LANG([C])
AC_USE_SYSTEM_EXTENSIONS
AC_PROG_CC
AC_PROG_CC_C99
AM_PROG_CC_C_O
AC_PROG_CPP
AC_TYPE_SIZE_T
AC_TYPE_INTPTR_T
AM_PROG_AR
LT_INIT

# Extra dependencies, configuration
###################################

AC_OPENMP

PKG_CHECK_MODULES([LIBNRM], [libnrm])
PKG_CHECK_MODULES([BLAS], [cblas])

AC_SEARCH_LIBS([log], [m], [], [
  AC_MSG_ERROR([unable to find the log() function])
])

# Feature flags
###############

# Select whether to validate benchmark computation after its execution.
# Default is to perform validation.
AC_ARG_ENABLE([post-validation],
  [AS_HELP_STRING([--disable-post-validation], [skip post validation])],
  [case "${enableval}" in
     yes | no ) enable_post_validation="${enableval}" ;;
     *) AC_MSG_ERROR([bad value ${enableval} for --disable-post-validation]) ;;
   esac],
  [enable_post_validation=yes]
)
AM_CONDITIONAL([ENABLE_POST_VALIDATION], [test "x${enable_post_validation}" = "xyes"])
AS_IF([test "x${enable_post_validation}" = "xyes"],
  [
    AC_DEFINE([ENABLE_POST_VALIDATION], [1], [Define to 1 to perform a validation of the benchmark computation])
  ],
  [
    AC_MSG_NOTICE([post validation disabled])
  ]
)

# Output
########

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

cat <<EOF

-------------------------------------------------------------------------------
NRM BENCHMARKS

Version: $PACKAGE_VERSION
OpenMP FLAGS: $OPENMP_CFLAGS
Libnrm CFLAGS: $LIBNRM_CFLAGS
Libnrm LDFLAGS: $LIBNRM_LIBS
BLAS CFLAGS: $BLAS_CFLAGS
BLAS LDFLAGS: $BLAS_LIBS

-------------------------------------------------------------------------------
EOF
