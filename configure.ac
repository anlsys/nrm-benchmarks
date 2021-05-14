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

AC_SEARCH_LIBS([log], [m], [], [
  AC_MSG_ERROR([unable to find the log() function])
])

# Output
########

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

cat <<EOF

-------------------------------------------------------------------------------
NRM BENCHMARKS

Version: $PACKAGE_VERSION
OpenMP Flags: $OPENMP_CFLAGS
Libnrm Cflags: $LIBNRM_CFLAGS
Libnrm LDflags: $LIBNRM_LIBS

-------------------------------------------------------------------------------
EOF
