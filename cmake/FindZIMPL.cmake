# - Try to find ZIMPL LP solver library
# See http://zimpl.zib.de/ for more information on ZIMPL
#
# Variables used by this module:
#
# ZIMPL_ROOT_DIR
#   Set this variable to the ZIMPL source or install path.
#   Otherwise, default paths are searched, e.g. /usr/local/
# ZIMPL_BUILD
#   Set this variable to "opt", "dbg" or "prf" for corresponding builds.
#   The default is "opt".
# ZIMPL_OSTYPE
#   Set this variable to any of ZIMPL's OS types, e.g. "linux", "win", etc.
#   The default is determined from `uname -s`.
# ZIMPL_ARCH
#   Set this variable to any of ZIMPL's architecture types, e.g. "x86", "x86_64", etc.
#   The default is determined from `uname -m`.
# ZIMPL_COMP
#   Set this variable to any of ZIMPL's compiler types, e.g. "gnu", "intel", etc.
#   The default is "gnu".
#
# Once done, this will define
#
#  ZIMPL_LIBRARIES      - list of libraries when using ZIMPL.
#  ZIMPL_FOUND          - true if ZIMPL was found.
#
# Author:
# Matthias Walter <matthias.walter@ovgu.de>
#
# Distributed under the Boost Software License, Version 1.0.
# (See http://www.boost.org/LICENSE_1_0.txt)

# Hints and paths for the search
set(_ZIMPL_ROOT_HINTS $ENV{ZIMPL_ROOT_DIR} ${ZIMPL_ROOT_DIR})
set(_ZIMPL_ROOT_PATHS $ENV{ZIMPL_ROOT_DIR} ${ZIMPL_ROOT_DIR})


# Read ZIMPL_BUILD from (environment) variable.
if (${ZIMPL_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_ZIMPL_BUILD ${ZIMPL_BUILD})
elseif ($ENV{ZIMPL_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_ZIMPL_BUILD $ENV{ZIMPL_BUILD})
else()
  set(_ZIMPL_BUILD "opt")
endif()

# Note: To see how ZIMPL determines OSTYPE and ARCH, look at zimpl/Makefile.

# Read ZIMPL_OSTYPE from (environment) variable or from `uname -s`
if (${ZIMPL_OSTYPE} MATCHES "^(aix|cygwin|darwin|freebsd|hp-ux|irix|linux|mingw|osf1|sunos|win)$")
  set(_ZIMPL_OSTYPE ${ZIMPL_OSTYPE})
elseif ($ENV{ZIMPL_OSTYPE} MATCHES "^(aix|cygwin|darwin|freebsd|hp-ux|irix|linux|mingw|osf1|sunos|win)$")
  set(_ZIMPL_OSTYPE $ENV{ZIMPL_OSTYPE})
else()
  execute_process(COMMAND uname -s OUTPUT_VARIABLE _ZIMPL_OSTYPE OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(TOLOWER ${_ZIMPL_OSTYPE} _ZIMPL_OSTYPE)
  string(REGEX REPLACE "cygwin.*" "cygwin" _ZIMPL_OSTYPE ${_ZIMPL_OSTYPE})
  string(REGEX REPLACE "irix.." "irix" _ZIMPL_OSTYPE ${_ZIMPL_OSTYPE})
  string(REGEX REPLACE "windows.*" "windows" _ZIMPL_OSTYPE ${_ZIMPL_OSTYPE})
  string(REGEX REPLACE "mingw.*" "mingw" _ZIMPL_OSTYPE ${_ZIMPL_OSTYPE})
endif()

# Read ZIMPL_ARCH from (environment) variable or from `uname -m`
if (${ZIMPL_ARCH} MATCHES "^(alpha|arm|clang|gnu|hppa|intel|mips|ppc|pwr4|sparc|x86|x86_64)$")
  set(_ZIMPL_ARCH ${ZIMPL_ARCH})
elseif ($ENV{ZIMPL_ARCH} MATCHES "^(alpha|arm|clang|gnu|hppa|intel|mips|ppc|pwr4|sparc|x86|x86_64)$")
  set(_ZIMPL_ARCH $ENV{ZIMPL_ARCH})
else()
  execute_process(COMMAND uname -m OUTPUT_VARIABLE _ZIMPL_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "sun.." "sparc" _ZIMPL_ARCH ${_ZIMPL_ARCH})
  string(REGEX REPLACE "i.86" "x86" _ZIMPL_ARCH ${_ZIMPL_ARCH})
  string(REGEX REPLACE "i86pc" "x86" _ZIMPL_ARCH ${_ZIMPL_ARCH})
  string(REGEX REPLACE "[0-9]86" "x86" _ZIMPL_ARCH ${_ZIMPL_ARCH})
  string(REGEX REPLACE "amd64" "x86_64" _ZIMPL_ARCH ${_ZIMPL_ARCH})
  string(REGEX REPLACE "IP.." "mips" _ZIMPL_ARCH ${_ZIMPL_ARCH})
  string(REGEX REPLACE "9000...." "hppa" _ZIMPL_ARCH ${_ZIMPL_ARCH})
  string(REGEX REPLACE "Power\ Macintosh" "ppc" _ZIMPL_ARCH ${_ZIMPL_ARCH})
  string(REGEX REPLACE "00.........." "pwr4" _ZIMPL_ARCH ${_ZIMPL_ARCH})
  string(REGEX REPLACE "arm.*" "arm" _ZIMPL_ARCH ${_ZIMPL_ARCH})
endif()

# Read ZIMPL_COMP from (environment) variable.
if (${ZIMPL_COMP} MATCHES "^(clang|compaq|gnu|hp|ibm|insure|intel|msv|purify|sgi|sun)$")
  set(_ZIMPL_COMP ${ZIMPL_COMP})
elseif ($ENV{ZIMPL_COMP} MATCHES "^(clang|compaq|gnu|hp|ibm|insure|intel|msv|purify|sgi|sun)$")
  set(_ZIMPL_COMP $ENV{ZIMPL_COMP})
else()
  set(_ZIMPL_COMP "gnu")
endif()

# Run `zimpl -V` to get the version.
execute_process(COMMAND zimpl -V OUTPUT_VARIABLE _ZIMPL_VERSION_STR OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REGEX REPLACE "^([0-9]).*$" "\\1" ZIMPL_VERSION_MAJOR "${_ZIMPL_VERSION_STR}")
string(REGEX REPLACE "^[0-9]\\.([0-9]).*$" "\\1" ZIMPL_VERSION_MINOR "${_ZIMPL_VERSION_STR}")
string(REGEX REPLACE "^[0-9]\\.[0-9]\\.([0-9]).*$" "\\1" ZIMPL_VERSION_PATCH "${_ZIMPL_VERSION_STR}")
set(ZIMPL_VERSION_STRING "${ZIMPL_VERSION_MAJOR}.${ZIMPL_VERSION_MINOR}.${ZIMPL_VERSION_PATCH}.${_ZIMPL_OSTYPE}.${_ZIMPL_ARCH}.${_ZIMPL_COMP}.${_ZIMPL_BUILD}")
 
# Search for library corresponding to version.
find_library(ZIMPL_LIBRARIES NAMES "zimpl-${ZIMPL_VERSION_STRING}" PATHS ${_ZIMPL_ROOT_PATHS} PATH_SUFFIXES lib)

# Let cmake process everything.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ZIMPL REQUIRED_VARS ZIMPL_LIBRARIES VERSION_VAR ZIMPL_VERSION_STRING)

