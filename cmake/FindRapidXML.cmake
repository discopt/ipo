# - Try to find RapidXML lib
#
# Once done this will define
#
#  RapidXML_FOUND - system has eigen lib with correct version
#  RapidXML_INCLUDE_DIR - the eigen include directory

# Breannan Smith (smith@cs.columbia.edu)

find_path(RapidXML_INCLUDE_DIR NAMES rapidxml.hpp
    PATHS
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_INSTALL_PREFIX}/include
    PATH_SUFFIXES rapidxml
  )

file(READ ${RapidXML_INCLUDE_DIR}/rapidxml.hpp HEADER)
string(REGEX MATCH "// Version *([0-9]*)\.([0-9]*) *" _ ${HEADER})
set(RapidXML_MAJOR ${CMAKE_MATCH_1})
set(RapidXML_MINOR ${CMAKE_MATCH_2})
set(RapidXML_VERSION "${RapidXML_MAJOR}.${RapidXML_MINOR}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(RapidXML DEFAULT_MSG RapidXML_VERSION)

if(RapidXML_FOUND AND NOT TARGET RapidXML::RapidXML)
  add_library(RapidXML INTERFACE)
  set_target_properties(RapidXML PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${RapidXML_INCLUDE_DIR}
  )
  add_library(RapidXML::RapidXML ALIAS RapidXML)
endif()

mark_as_advanced(RapidXML_INCLUDE_DIR)
