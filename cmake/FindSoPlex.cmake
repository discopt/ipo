# Tries to find the SoPlex library.
#
# Parameter Variables:
#
# SoPlex_ROOT_DIR
#   Set this variable to the SoPlex source or install path.
#   Otherwise, default paths are searched, e.g. /usr/local/
# SoPlex_BUILD
#   Set this variable to "opt", "dbg" or "prf" for corresponding builds.
#   The default is "opt".
# SoPlex_OSTYPE
#   Set this variable to any of SoPlex's OS types, e.g. "linux", "win", etc.
#   The default is determined from `uname -s`.
# SoPlex_ARCH
#   Set this variable to any of SoPlex's architecture types, e.g. "x86", "x86_64", etc.
#   The default is determined from `uname -m`.
# SoPlex_COMP
#   Set this variable to any of SoPlex's compiler types, e.g. "gnu", "intel", etc.
#   The default is "gnu".
# SoPlex_USE_STATIC_LIBS
#   Set to TRUE for linking with static library.
#
# Defines Variables:
#
# SoPlex_FOUND
#   True if SoPlex was found.
# SoPlex_INCLUDE_DIRS
#   Include directories.
# SoPlex_LIBRARIES
#   Path of libraries.
# SoPlex_VERSION
#   Version found.
#
# Author:
# 
# Matthias Walter <matthias@matthiaswalter.org>
#
# Distributed under the Boost Software License, Version 1.0.
# (See http://www.boost.org/LICENSE_1_0.txt)

# Dependencies.
find_package(UBSan)
if (HAVE_UNDEFINED_BEHAVIOR_SANITIZER)
  set(UBSAN_LIBRARIES "-lubsan")
else()
  set(UBSAN_LIBRARIES "")
endif()

if(SoPlex_FIND_REQUIRED)
  find_package(GMP REQUIRED)
  find_package(ZLIB REQUIRED)
else()
  find_package(GMP)
  find_package(ZLIB)
endif()

# Handle SoPlex_ROOT_DIR.
set(_SoPlex_ROOT_HINTS ${SoPlex_ROOT_DIR} ENV ${SoPlex_ROOT_DIR})

# Handle SoPlex_USE_STATIC_LIBS.
if(SoPlex_USE_STATIC_LIBS)
  set(_SoPlex_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a )
  endif()
endif()

# Read SoPlex_BUILD from (environment) variable.
if (${SoPlex_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SoPlex_BUILD ${SoPlex_BUILD})
elseif ($ENV{SoPlex_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SoPlex_BUILD $ENV{SoPlex_BUILD})
else()
  set(_SoPlex_BUILD "opt")
endif()

# Note: To see how SoPlex determines OSTYPE and ARCH, look at soplex/make/make.detecthost.

# Read SoPlex_OSTYPE from (environment) variable or from `uname -s`
if (${SoPlex_OSTYPE} MATCHES "^(aix|cygwin|darwin|freebsd|hp-ux|irix|linux|mingw|osf1|sunos|win)$")
  set(_SoPlex_OSTYPE ${SoPlex_OSTYPE})
elseif ($ENV{SoPlex_OSTYPE} MATCHES "^(aix|cygwin|darwin|freebsd|hp-ux|irix|linux|mingw|osf1|sunos|win)$")
  set(_SoPlex_OSTYPE $ENV{SoPlex_OSTYPE})
else()
  execute_process(COMMAND uname -s OUTPUT_VARIABLE _SoPlex_OSTYPE OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(TOLOWER ${_SoPlex_OSTYPE} _SoPlex_OSTYPE)
  string(REGEX REPLACE "cygwin.*" "cygwin" _SoPlex_OSTYPE ${_SoPlex_OSTYPE})
  string(REGEX REPLACE "irix.." "irix" _SoPlex_OSTYPE ${_SoPlex_OSTYPE})
  string(REGEX REPLACE "windows.*" "windows" _SoPlex_OSTYPE ${_SoPlex_OSTYPE})
  string(REGEX REPLACE "mingw.*" "mingw" _SoPlex_OSTYPE ${_SoPlex_OSTYPE})
endif()

# Read SoPlex_ARCH from (environment) variable or from `uname -m`
if (${SoPlex_ARCH} MATCHES "^(alpha|arm|clang|gnu|hppa|intel|mips|ppc|pwr4|sparc|x86|x86_64)$")
  set(_SoPlex_ARCH ${SoPlex_ARCH})
elseif ($ENV{SoPlex_ARCH} MATCHES "^(alpha|arm|clang|gnu|hppa|intel|mips|ppc|pwr4|sparc|x86|x86_64)$")
  set(_SoPlex_ARCH $ENV{SoPlex_ARCH})
else()
  execute_process(COMMAND uname -m OUTPUT_VARIABLE _SoPlex_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "sun.." "sparc" _SoPlex_ARCH ${_SoPlex_ARCH})
  string(REGEX REPLACE "i.86" "x86" _SoPlex_ARCH ${_SoPlex_ARCH})
  string(REGEX REPLACE "i86pc" "x86" _SoPlex_ARCH ${_SoPlex_ARCH})
  string(REGEX REPLACE "[0-9]86" "x86" _SoPlex_ARCH ${_SoPlex_ARCH})
  string(REGEX REPLACE "amd64" "x86_64" _SoPlex_ARCH ${_SoPlex_ARCH})
  string(REGEX REPLACE "IP.." "mips" _SoPlex_ARCH ${_SoPlex_ARCH})
  string(REGEX REPLACE "9000...." "hppa" _SoPlex_ARCH ${_SoPlex_ARCH})
  string(REGEX REPLACE "Power\ Macintosh" "ppc" _SoPlex_ARCH ${_SoPlex_ARCH})
  string(REGEX REPLACE "00.........." "pwr4" _SoPlex_ARCH ${_SoPlex_ARCH})
  string(REGEX REPLACE "arm.*" "arm" _SoPlex_ARCH ${_SoPlex_ARCH})
endif()

# Read SoPlex_COMP from (environment) variable.
if (${SoPlex_COMP} MATCHES "^(clang|compaq|gnu|hp|ibm|insure|intel|msv|purify|sgi|sun)$")
  set(_SoPlex_COMP ${SoPlex_COMP})
elseif ($ENV{SoPlex_COMP} MATCHES "^(clang|compaq|gnu|hp|ibm|insure|intel|msv|purify|sgi|sun)$")
  set(_SoPlex_COMP $ENV{SoPlex_COMP})
else()
  set(_SoPlex_COMP "gnu")
endif()

# Root path.
find_path(SoPlex_ROOT_DIR NAMES include/soplex.h soplex.h HINTS ${_SoPlex_ROOT_HINTS}
  NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH)
find_path(SoPlex_ROOT_DIR NAMES include/soplex.h soplex.h HINTS ${_SoPlex_ROOT_HINTS})

# Includes
find_path(_SoPlex_INCLUDE_DIR NAMES soplex.h PATHS ${SoPlex_ROOT_DIR} PATH_SUFFIXES include src
  NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH)
find_path(_SoPlex_INCLUDE_DIR NAMES soplex.h PATHS ${SoPlex_ROOT_DIR} PATH_SUFFIXES include src)

if (_SoPlex_INCLUDE_DIR)
  set(SoPlex_INCLUDE_DIRS ${_SoPlex_INCLUDE_DIR})
  
  # Extract version from spxdefines.h
  file(STRINGS "${_SoPlex_INCLUDE_DIR}/spxdefines.h" _SoPlex_VERSION_STR REGEX "^#define[\t ]+SOPLEX_VERSION[\t ]+[0-9][0-9][0-9].*")
  string(REGEX REPLACE "^.*SOPLEX_VERSION[\t ]+([0-9]).*$" "\\1" SoPlex_VERSION_MAJOR "${_SoPlex_VERSION_STR}")
  string(REGEX REPLACE "^.*SOPLEX_VERSION[\t ]+[0-9]([0-9]).*$" "\\1" SoPlex_VERSION_MINOR "${_SoPlex_VERSION_STR}")
  string(REGEX REPLACE "^.*SOPLEX_VERSION[\t ]+[0-9][0-9]([0-9]).*$" "\\1" SoPlex_VERSION_PATCH "${_SoPlex_VERSION_STR}")
  file(STRINGS "${_SoPlex_INCLUDE_DIR}/spxdefines.h" _SoPlex_SUBVERSION_STR REGEX "^#define[\t ]+SOPLEX_SUBVERSION[\t ]+[0-9].*")
  string(REGEX REPLACE "^.*SOPLEX_SUBVERSION[\t ]+([0-9]).*$" "\\1" SoPlex_VERSION_SUBVERSION "${_SoPlex_SUBVERSION_STR}")

  # Search for library with version (release) or version.subversion (development).

  set(SoPlex_VERSION "${SoPlex_VERSION_MAJOR}.${SoPlex_VERSION_MINOR}.${SoPlex_VERSION_PATCH}")
  set(_SoPlex_VERSION_SUFFIX "${_SoPlex_OSTYPE}.${_SoPlex_ARCH}.${_SoPlex_COMP}.${_SoPlex_BUILD}")
  find_library(SoPlex_LIBRARY NAMES "soplex-${SoPlex_VERSION}.${_SoPlex_VERSION_SUFFIX}" "soplex-${SoPlex_VERSION}.${SoPlex_VERSION_SUBVERSION}.${_SoPlex_VERSION_SUFFIX}" PATHS ${SoPlex_ROOT_DIR} PATH_SUFFIXES lib NO_DEFAULT_PATH)
  if(NOT SoPlex_LIBRARY)
    find_library(SoPlex_LIBRARY NAMES "soplex-${SoPlex_VERSION}.${_SoPlex_VERSION_SUFFIX}" "soplex-${SoPlex_VERSION}.${SoPlex_VERSION_SUBVERSION}.${_SoPlex_VERSION_SUFFIX}" PATHS ${SoPlex_ROOT_DIR} PATH_SUFFIXES lib)
  endif()
  if(SoPlex_LIBRARY)
    string(REGEX REPLACE "^.*libsoplex-(.*)\\..*$" "\\1" SoPlex_VERSION_STRING "${SoPlex_LIBRARY}")
    set(SoPlex_LIBRARIES ${SoPlex_LIBRARY} ${GMP_LIBRARIES} ${ZLIB_LIBRARIES} ${UBSAN_LIBRARIES})
  endif()
endif()

# Let cmake process everything.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SoPlex REQUIRED_VARS SoPlex_LIBRARY SoPlex_ROOT_DIR SoPlex_INCLUDE_DIRS SoPlex_LIBRARIES VERSION_VAR SoPlex_VERSION)

# Restore the original find_library ordering.
if(SoPlex_USE_STATIC_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_SoPlex_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

