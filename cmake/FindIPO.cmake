# Hints and paths for the search
set(_IPO_ROOT_HINTS $ENV{IPO_ROOT_DIR} ${IPO_ROOT_DIR})
set(_IPO_ROOT_PATHS $ENV{IPO_ROOT_DIR} ${IPO_ROOT_DIR})

find_package(SoPlex REQUIRED)

# Root

find_path(IPO_ROOT_DIR NAMES include/ipo/oracles.h ipo/oracles.h HINTS ${_IPO_ROOT_HINTS} PATHS ${_IPO_ROOT_PATHS})

# Includes

find_path(_IPO_INCLUDE NAMES ipo/oracles.h PATHS ${IPO_ROOT_DIR} PATH_SUFFIXES include)
if (_IPO_INCLUDE)
  set(IPO_INCLUDE_DIRS ${IPO_INCLUDE_DIRS} ${SoPlex_INCLUDE_DIRS})
endif()

# Libraries

find_library(_IPO_LIB NAMES ipo PATHS ${IPO_ROOT_DIR} PATH_SUFFIXES lib)

if (_IPO_LIB)
  set(IPO_LIBRARIES ${_IPO_LIB} ${SoPlex_LIBRARIES})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPO REQUIRED_VARS IPO_ROOT_DIR IPO_INCLUDE_DIRS IPO_LIBRARIES)
