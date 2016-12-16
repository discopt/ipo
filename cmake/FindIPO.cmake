# Tries to find the IPO library.
#
# Parameter Variables:
#
# IPO_ROOT_DIR
#   Root directory of IPO installation.
# IPO_USE_STATIC_LIBS
#   Set to TRUE for linking with static library.
#
# Defines Variables:
#
# IPO_FOUND
#   True if IPO was found.
# IPO_INCLUDE_DIRS
#   Directory where gmp.h resides.
# IPO_LIBRARIES
#   Libraries including dependencies.
# IPO_VERSION
#   Version found.
#
# Author:
# 
# Matthias Walter <matthias@matthiaswalter.org>
#
# Distributed under the Boost Software License, Version 1.0.
# (See http://www.boost.org/LICENSE_1_0.txt)

# Dependencies.
if(IPO_FIND_REQUIRED)
  find_package(SoPlex REQUIRED)
else()
  find_package(SoPlex)
endif()
find_package(SCIP)

# Handle ROOT_DIR.
set(_IPO_ROOT_HINTS ${IPO_ROOT_DIR} ENV ${IPO_ROOT_DIR})

# Root
find_path(IPO_ROOT_DIR NAMES include/ipo/oracles.h ipo/oracles.h HINTS ${_IPO_ROOT_HINTS})

# Handle IPO_USE_STATIC_LIBS.
if(IPO_USE_STATIC_LIBS)
  set(_IPO_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a )
  endif()
endif()

# Includes
find_path(_IPO_INCLUDE NAMES ipo/oracles.h PATHS ${IPO_ROOT_DIR} PATH_SUFFIXES include)
if(_IPO_INCLUDE)
  set(IPO_INCLUDE_DIRS ${IPO_INCLUDE_DIRS} ${SOPLEX_INCLUDE_DIRS})
endif()

# Libraries
find_library(_IPO_LIB NAMES ipo PATHS ${IPO_ROOT_DIR} PATH_SUFFIXES lib)

if(_IPO_LIB)
  set(IPO_LIBRARIES ${_IPO_LIB} ${SOPLEX_LIBRARIES})
  if(SCIP_FOUND)
    set(IPO_INCLUDE_DIRS ${IPO_INCLUDE_DIRS} ${SCIP_INCLUDE_DIRS})
    set(IPO_LIBRARIES ${IPO_LIBRARIES} ${SCIP_LIBRARIES})
  endif()
endif()

# Let cmake process everything.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPO REQUIRED_VARS IPO_ROOT_DIR IPO_INCLUDE_DIRS IPO_LIBRARIES)

# Restore the original find_library ordering.
if(IPO_USE_STATIC_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_IPO_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()
