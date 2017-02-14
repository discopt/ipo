# Tries to find the SCIP library.
#
# Parameter Variables:
#
# SCIP_ROOT_DIR
#   Set this variable to the SCIP source or install path.
#   Otherwise, default paths are searched, e.g. /usr/local/
# SCIP_LPS
#   Set this variable to "spx" (="spx2") or "spx1" (others not implemented, yet).
#   The default is "spx".
# SCIP_TPI
#   Set this variable to "none", "openmp" or "tnycthrd".
#   The default is "none"
# SCIP_BUILD
#   Set this variable to "opt", "dbg" or "prf" for corresponding builds.
#   The default is "opt".
# SCIP_LPS_BUILD
#   Set this variable to "opt", "dbg" or "prf" for corresponding builds of the LP solver.
#   The default is "opt".
#   Note that this value is ignored if the LP solver was already detected,
#   e.g., if find_package(SoPlex) was called before calling find_package(SCIP).
# SCIP_OSTYPE
#   Set this variable to any of SCIP's OS types, e.g. "linux", "win", etc.
#   The default is determined from `uname -s`.
# SCIP_ARCH
#   Set this variable to any of SCIP's architecture types, e.g. "x86", "x86_64", etc.
#   The default is determined from `uname -m`.
# SCIP_COMP
#   Set this variable to any of SCIP's compiler types, e.g. "gnu", "intel", etc.
#   The default is "gnu".
# SCIP_USE_STATIC_LIBS
#   Set to TRUE for linking with static library.
#
# Defines Variables:
#
# SCIP_FOUND
#   True if SCIP was found.
# SCIP_INCLUDE_DIRS
#   Include directories.
# SCIP_LIBRARIES
#   Path of libraries.
# SCIP_VERSION
#   Version found.
# SCIP_LPS
#   LP solver used by SCIP.
#
# Author:
# 
# Matthias Walter <matthias@matthiaswalter.org>
#
# Distributed under the Boost Software License, Version 1.0.
# (See http://www.boost.org/LICENSE_1_0.txt)

# Dependencies.
if(SCIP_FIND_REQUIRED)
  find_package(GMP REQUIRED)
  find_package(ZLIB REQUIRED)
  find_package(Readline REQUIRED)
else()
  find_package(GMP)
  find_package(ZLIB)
  find_package(Readline)
endif()

find_package(ZIMPL)

# Handle ROOT_DIR.
set(_SCIP_ROOT_HINTS ${SCIP_ROOT_DIR} ENV ${SCIP_ROOT_DIR})

# Handle SCIP_USE_STATIC_LIBS.
if(SCIP_USE_STATIC_LIBS)
  set(_SCIP_LIB_PATHS lib lib/static)
  set(_SCIP_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a )
  endif()
else()
  set(_SCIP_LIB_PATHS lib lib/shared)
endif()

# Read SCIP_BUILD from (environment) variable.
if (${SCIP_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SCIP_BUILD ${SCIP_BUILD})
elseif ($ENV{SCIP_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SCIP_BUILD $ENV{SCIP_BUILD})
else()
  set(_SCIP_BUILD "opt")
endif()

# Read SCIP_LPS from (environment) variable.
if (${SCIP_LPS} MATCHES "^(spx|spx2|spx1)$")
elseif ($ENV{SCIP_LPS} MATCHES "^(spx|spx2|spx1)$")
  set(SCIP_LPS $ENV{SCIP_LPS})
else()
  set(SCIP_LPS "spx")
endif()

# Read SCIP_TPI from (environment) variable.
if (${SCIP_TPI} MATCHES "^(none|openmp|tnycthrd)$")
elseif ($ENV{SCIP_TPI} MATCHES "^(none|openmp|tnycthrd)$")
  set(SCIP_TPI $ENV{SCIP_TPI})
else()
  set(SCIP_TPI "none")
endif()

# Read SCIP_LPS_BUILD from (environment) variable.
if (${SCIP_LPS_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SCIP_LPS_BUILD ${SCIP_LPS_BUILD})
elseif ($ENV{SCIP_LPS_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SCIP_LPS_BUILD $ENV{SCIP_LPS_BUILD})
else()
  set(_SCIP_LPS_BUILD "opt")
endif()

# Note: To see how SCIP determines OSTYPE and ARCH, look at scip/make/make.detecthost.

# Read SCIP_OSTYPE from (environment) variable or from `uname -s`
if (${SCIP_OSTYPE} MATCHES "^(aix|cygwin|darwin|freebsd|hp-ux|irix|linux|mingw|osf1|sunos|win)$")
  set(_SCIP_OSTYPE ${SCIP_OSTYPE})
elseif ($ENV{SCIP_OSTYPE} MATCHES "^(aix|cygwin|darwin|freebsd|hp-ux|irix|linux|mingw|osf1|sunos|win)$")
  set(_SCIP_OSTYPE $ENV{SCIP_OSTYPE})
else()
  execute_process(COMMAND uname -s OUTPUT_VARIABLE _SCIP_OSTYPE OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(TOLOWER ${_SCIP_OSTYPE} _SCIP_OSTYPE)
  string(REGEX REPLACE "cygwin.*" "cygwin" _SCIP_OSTYPE ${_SCIP_OSTYPE})
  string(REGEX REPLACE "irix.." "irix" _SCIP_OSTYPE ${_SCIP_OSTYPE})
  string(REGEX REPLACE "windows.*" "windows" _SCIP_OSTYPE ${_SCIP_OSTYPE})
  string(REGEX REPLACE "mingw.*" "mingw" _SCIP_OSTYPE ${_SCIP_OSTYPE})
endif()

# Read SCIP_ARCH from (environment) variable or from `uname -m`
if (${SCIP_ARCH} MATCHES "^(alpha|arm|clang|gnu|hppa|intel|mips|ppc|pwr4|sparc|x86|x86_64)$")
  set(_SCIP_ARCH ${SCIP_ARCH})
elseif ($ENV{SCIP_ARCH} MATCHES "^(alpha|arm|clang|gnu|hppa|intel|mips|ppc|pwr4|sparc|x86|x86_64)$")
  set(_SCIP_ARCH $ENV{SCIP_ARCH})
else()
  execute_process(COMMAND uname -m OUTPUT_VARIABLE _SCIP_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "sun.." "sparc" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "i.86" "x86" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "i86pc" "x86" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "[0-9]86" "x86" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "amd64" "x86_64" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "IP.." "mips" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "9000...." "hppa" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "Power\ Macintosh" "ppc" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "00.........." "pwr4" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "arm.*" "arm" _SCIP_ARCH ${_SCIP_ARCH})
endif()

# Read SCIP_COMP from (environment) variable.
if (${SCIP_COMP} MATCHES "^(clang|compaq|gnu|hp|ibm|insure|intel|msv|purify|sgi|sun)$")
  set(_SCIP_COMP ${SCIP_COMP})
elseif ($ENV{SCIP_COMP} MATCHES "^(clang|compaq|gnu|hp|ibm|insure|intel|msv|purify|sgi|sun)$")
  set(_SCIP_COMP $ENV{SCIP_COMP})
else()
  set(_SCIP_COMP "gnu")
endif()

# Root path.
find_path(SCIP_ROOT_DIR NAMES include/scip/scip.h src/scip/scip.h HINTS ${_SCIP_ROOT_HINTS})

# Include directory.
find_path(SCIP_INCLUDE_DIR NAMES scip/scip.h HINTS ${SCIP_ROOT_DIR} PATH_SUFFIXES include src
  NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH)
find_path(SCIP_INCLUDE_DIR NAMES scip/scip.h HINTS ${SCIP_ROOT_DIR} PATH_SUFFIXES include src)

if (SCIP_INCLUDE_DIR)
  set(SCIP_INCLUDE_DIRS ${SCIP_INCLUDE_DIR})

  # Extract version from scip/def.h
  file(STRINGS "${SCIP_INCLUDE_DIR}/scip/def.h" _SCIP_VERSION_STR REGEX "^#define[\t ]+SCIP_VERSION[\t ]+[0-9][0-9][0-9].*")
  string(REGEX REPLACE "^.*SCIP_VERSION[\t ]+([0-9]).*$" "\\1" SCIP_VERSION_MAJOR "${_SCIP_VERSION_STR}")
  string(REGEX REPLACE "^.*SCIP_VERSION[\t ]+[0-9]([0-9]).*$" "\\1" SCIP_VERSION_MINOR "${_SCIP_VERSION_STR}")
  string(REGEX REPLACE "^.*SCIP_VERSION[\t ]+[0-9][0-9]([0-9]).*$" "\\1" SCIP_VERSION_PATCH "${_SCIP_VERSION_STR}")
  file(STRINGS "${SCIP_INCLUDE_DIR}/scip/def.h" _SCIP_SUBVERSION_STR REGEX "^#define[\t ]+SCIP_SUBVERSION[\t ]+[0-9].*")
  string(REGEX REPLACE "^.*SCIP_SUBVERSION[\t ]+([0-9]).*$" "\\1" SCIP_VERSION_SUBVERSION "${_SCIP_SUBVERSION_STR}")

#  set(SCIP_VERSION_STRING "${SCIP_VERSION_MAJOR}.${SCIP_VERSION_MINOR}.${SCIP_VERSION_PATCH}.${SCIP_VERSION_SUBVERSION}.${_SCIP_OSTYPE}.${_SCIP_ARCH}.${_SCIP_COMP}.${_SCIP_BUILD}")

  set(_SCIP_FOUND_ALL TRUE)

  # Search for libscip corresponding to version.
  set(SCIP_VERSION "${SCIP_VERSION_MAJOR}.${SCIP_VERSION_MINOR}.${SCIP_VERSION_PATCH}")
  set(_SCIP_VERSION_SUFFIX "${_SCIP_OSTYPE}.${_SCIP_ARCH}.${_SCIP_COMP}.${_SCIP_BUILD}")

  find_library(SCIP_SCIP_LIBRARY NAMES "scip-${SCIP_VERSION}.${_SCIP_VERSION_SUFFIX}" "scip-${SCIP_VERSION}.${SCIP_VERSION_SUBVERSION}.${_SCIP_VERSION_SUFFIX}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS} NO_DEFAULT_PATH)

  if(NOT SCIP_SCIP_LIBRARY)
    find_library(SCIP_SCIP_LIBRARY NAMES "scip-${SCIP_VERSION}.${_SCIP_VERSION_SUFFIX}" "scip-${SCIP_VERSION}.${SCIP_VERSION_SUBVERSION}.${_SCIP_VERSION_SUFFIX}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS})
  endif()

  if (SCIP_SCIP_LIBRARY)
    string(REGEX REPLACE "^.*libscip-(.*)\\..*$" "\\1" SCIP_VERSION_STRING "${SCIP_SCIP_LIBRARY}")
  else()
    set(_SCIP_FOUND_ALL FALSE)
    if(NOT SCIP_FIND_QUIETLY)
      message(STATUS "SCIP library libscip-${SCIP_VERSION}[.${SCIP_VERSION_SUBVERSION}].${_SCIP_VERSION_SUFFIX} was not found. Search paths: ${_SCIP_ROOT_HINTS}")

      # Check if some other version was found and report to user.
      find_library(SCIP_SCIP_LIBRARY_TEST NAMES "scip" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS})
      if (SCIP_SCIP_LIBRARY_TEST)
        message(STATUS "Found SCIP library ${SCIP_SCIP_LIBRARY_TEST} different from the one promised by ${SCIP_INCLUDE_DIRS}/scip/def.h")
      endif()
    endif()
  endif()
  
  # Search for libobjscip
  find_library(SCIP_OBJSCIP_LIBRARY NAMES "objscip-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS} NO_DEFAULT_PATH)
  if(NOT ${SCIP_OBJSCIP_LIBRARY})
    find_library(SCIP_OBJSCIP_LIBRARY NAMES "objscip-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS})
  endif()
  if(NOT ${SCIP_OBJSCIP_LIBRARY} MATCHES "objscip")
    set(_SCIP_FOUND_ALL FALSE)
    if(NOT SCIP_FIND_QUIETLY)
      message(STATUS "SCIP library libobjscip-${SCIP_VERSION_STRING} was not found.")
    endif()
  endif()
  
  # Search for nlpi. TODO: cppad is currently hard-coded, while ipopt is not recognized.
  find_library(SCIP_NLPI_LIBRARY NAMES "nlpi.cppad-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS} NO_DEFAULT_PATH)
  if(NOT ${SCIP_NLPI_LIBRARY})
    find_library(SCIP_NLPI_LIBRARY NAMES "nlpi.cppad-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS})
  endif()
  if(NOT ${SCIP_NLPI_LIBRARY} MATCHES "nlpi")
    set(_SCIP_FOUND_ALL FALSE)
    if(NOT SCIP_FIND_QUIETLY)
      message(STATUS "SCIP library libnlpi.cppad-${SCIP_VERSION_STRING} was not found.")
    endif()
  endif()

  if(${SCIP_VERSION_MAJOR} GREATER 3)
    # Search for libtpi${SCIP_TPI}
    find_library(SCIP_TPI_LIBRARY NAMES "tpi${SCIP_TPI}-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS} NO_DEFAULT_PATH)
    if(NOT ${SCIP_TPI_LIBRARY})
      find_library(SCIP_TPI_LIBRARY NAMES "lpi${SCIP_TPI}-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS})
    endif()
    if(NOT ${SCIP_TPI_LIBRARY} MATCHES "${SCIP_TPI}")
      set(_SCIP_FOUND_ALL FALSE)
      if(NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP library libtpi${SCIP_TPI}-${SCIP_VERSION_STRING} was not found.")
      endif()
    endif()
  else()
    set(SCIP_TPI_LIBRARY)
  endif()

  # Search for the LP solver: spx
  if(${SCIP_LPS} STREQUAL "spx")
    # Search for liblpispx
    find_library(SCIP_LPI_LIBRARY NAMES "lpispx-${SCIP_VERSION_STRING}" "lpispx2-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS} NO_DEFAULT_PATH)
    if(NOT ${SCIP_LPI_LIBRARY})
      find_library(SCIP_LPI_LIBRARY NAMES "lpispx-${SCIP_VERSION_STRING}" "lpispx2-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS})
    endif()
    if(NOT ${SCIP_LPI_LIBRARY} MATCHES "spx")
      set(_SCIP_FOUND_ALL FALSE)
      if(NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP library liblpispx[2]-${SCIP_VERSION_STRING} was not found.")
      endif()
    endif()
  
    # Search for SoPlex.
    if (NOT SoPlex_FOUND AND NOT SOPLEX_FOUND)
      set(SoPlex_BUILD ${SCIP_LPS_BUILD})
      find_package(SoPlex REQUIRED)
    endif()
    if (SoPlex_FOUND OR SOPLEX_FOUND)
      set(_SCIP_LIB_LPSOLVER ${SoPlex_LIBRARIES})
    else()
      set(_SCIP_FOUND_ALL FALSE)
      if (NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP dependency SoPlex was not found.")
      endif()
    endif()
  endif()
  
  # Search for the LP solver: spx2
  if (${SCIP_LPS} STREQUAL "spx2")
    # Search for liblpispx2
    find_library(SCIP_LPI_LIBRARY NAMES "lpispx2-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS} NO_DEFAULT_PATH)
    if(NOT ${SCIP_LPI_LIBRARY})
      find_library(SCIP_LPI_LIBRARY NAMES "lpispx2-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS})
    endif()
    if(NOT ${SCIP_LPI_LIBRARY} MATCHES "spx2")
      set(_SCIP_FOUND_ALL FALSE)
      if (NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP library liblpispx2-${SCIP_VERSION_STRING} was not found.")
      endif()
    endif()
  
    # Search for SoPlex.
    if (NOT SoPlex_FOUND AND NOT SOPLEX_FOUND)
      set(SoPlex_BUILD ${SCIP_LPS_BUILD})
      find_package(SoPlex REQUIRED)
    endif()
    if (SoPlex_FOUND OR SOPLEX_FOUND)
      set(_SCIP_LIB_LPSOLVER ${SoPlex_LIBRARIES})
    else()
      set(_SCIP_FOUND_ALL FALSE)
      if (NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP dependency SoPlex was not found.")
      endif()
    endif()
  endif()
  
  # Search for the LP solver: spx1
  if (${SCIP_LPS} STREQUAL "spx1")
    # Search for liblpispx1
    find_library(SCIP_LPI_LIBRARY NAMES "lpispx1-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS} NO_DEFAULT_PATH)
    if(NOT ${SCIP_LPI_LIBRARY})
      find_library(SCIP_LPI_LIBRARY NAMES "lpispx1-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${_SCIP_LIB_PATHS})
    endif()
    if(NOT ${SCIP_LPI_LIBRARY} MATCHES "liblpispx1")
      set(_SCIP_FOUND_ALL FALSE)
      if(NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP library liblpispx1-${SCIP_VERSION_STRING} was not found.")
      endif()
    endif()
  
    # Search for SoPlex.
    if (NOT SoPlex_FOUND AND NOT SOPLEX_FOUND)
      set(SoPlex_BUILD ${SCIP_LPS_BUILD})
      find_package(SoPlex REQUIRED)
    endif()
    if (SoPlex_FOUND OR SOPLEX_FOUND)
      set(_SCIP_LIB_LPSOLVER ${SoPlex_LIBRARIES})
    else()
      set(_SCIP_FOUND_ALL FALSE)
      if (NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP dependency SoPlex was not found.")
      endif()
    endif()
  endif()
  
  if (_SCIP_FOUND_ALL)
    set(SCIP_LIBRARIES
      ${SCIP_OBJSCIP_LIBRARY} ${SCIP_SCIP_LIBRARY} ${SCIP_TPI_LIBRARY} ${SCIP_LPI_LIBRARY} ${SCIP_NLPI_LIBRARY} ${_SCIP_LIB_LPSOLVER}
      ${ZLIB_LIBRARIES} ${Readline_LIBRARY} ${GMP_LIBRARIES}
    )
    if (ZIMPL_FOUND)
      set(SCIP_LIBRARIES ${SCIP_LIBRARIES} ${ZIMPL_LIBRARIES})
    endif()
  else()
    set(SCIP_VERSION_STRING)
  endif()
endif()

# Let cmake process everything.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCIP REQUIRED_VARS SCIP_SCIP_LIBRARY SCIP_ROOT_DIR SCIP_INCLUDE_DIRS SCIP_LIBRARIES SCIP_LPS VERSION_VAR SCIP_VERSION)

# Restore the original find_library ordering.
if(SCIP_USE_STATIC_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_SCIP_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

