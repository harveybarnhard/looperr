AC_PREREQ(2.61)

AC_INIT([looperr], m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))


## Determine Install Location of R
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi

## Setup RBin
RBIN="${R_HOME}/bin/R"
CXX=`"${RBIN}" CMD config CXX`
CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`
CXXFLAGS=`"${RBIN}" CMD config CXXFLAGS`

## We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP

## Check the C++ compiler using the CXX value set
AC_PROG_CXX

## If the compiler is g++, then GXX is set so version information can be exaimed
if test "${GXX}" = yes; then
    AC_MSG_CHECKING([whether g++ version is sufficient])
    gxx_version=$(${CXX} -v 2>&1 | awk '/^.*g.. version/ {print $3}')
    case ${gxx_version} in
        4.9.*)
             AC_MSG_RESULT([(${gxx_version}) yes])
        ;;

 1.*|2.*|3.*|4.0.*|4.1.*|4.2.*|4.3.*|4.4.*|4.5.*|4.6.*|4.7.*|4.8.*|5.*|6.*|7.*)
             AC_MSG_RESULT([no])
             AC_MSG_WARN([Only g++ version 4.9.* be used to compile this package.])
             AC_MSG_ERROR([Please set the default compiler to gcc++-4.9.])
    ;;
    esac
else
 AC_MSG_WARN([The package uses an external shared library that only compiles with gcc++-4.9])
 AC_MSG_ERROR([Please set the default compiler to gcc++-4.9.])
fi

AC_OUTPUT
