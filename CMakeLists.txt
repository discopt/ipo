cmake_minimum_required(VERSION 3.5)

project(IPO VERSION 2.0 LANGUAGES CXX)
   
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)

set(RATIONAL AUTO CACHE STRING "Support rational arithmetic algorithms; requires at least GMP and Boost 1.65.")
set_property(CACHE RATIONAL PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER ${RATIONAL} RATIONAL_UPPER) # Also allow "auto"

set(SCIP AUTO CACHE STRING "Support SCIP as a MIP solver.")
set_property(CACHE SCIP PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER ${SCIP} SCIP_UPPER) # Also allow "auto"

set(SOPLEX AUTO CACHE STRING "Support SoPlex as an LP solver.")
set_property(CACHE SOPLEX PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER ${SOPLEX} SOPLEX_UPPER) # Also allow "auto"

set(GUROBI AUTO CACHE STRING "Support Gurobi as a MIP and LP solver.")
set_property(CACHE GUROBI PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER ${GUROBI} GUROBI_UPPER) # Also allow "auto"


option(SOPLEX "Support SoPlex as an LP solver." ON)
option(GUROBI "Support Gurobi as MIP and LP solver." ON)
set(LP AUTO CACHE STRING "Default floating-point LP solver")
set_property(CACHE LP PROPERTY STRINGS AUTO GUROBI SOPLEX OFF)

option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})

# Add cmake/ to CMAKE_MODULE_PATH.
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

# Set default build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release"
    CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Find dependencies.

if(RATIONAL_UPPER)
  find_package(GMP)
  if(GMP_FOUND)
    find_package(Boost 1.65.0)
    if (Boost_FOUND)
      set(IPO_WITH_RATIONAL ON)
    elseif(RATIONAL_UPPER STREQUAL "AUTO")
      message(STATUS "Boost >= 1.65 not found. Disabling rational arithmetic.")
      set(IPO_WITH_RATIONAL OFF)
    else()
      message(FATAL_ERROR "Boost >= 1.65 not found, but required for rational arithmetic.")
    endif()
  elseif(RATIONAL_UPPER STREQUAL "AUTO")
    message(STATUS "GMP not found. Disabling rational arithmetic.")
  else()
    message(FATAL_ERROR "GMP not found, but required for rational arithmetic.")
    set(IPO_WITH_RATIONAL OFF)
  endif()
else()
  set(IPO_WITH_RATIONAL OFF)
endif()

if(SCIP_UPPER)
  find_package(SCIP CONFIG HINTS ${SCIP_DIR})
  if(SCIP_FOUND)
    message(STATUS "Found SCIP: ${SCIP_INCLUDE_DIRS} / ${SCIP_LIBRARIES}.")
  elseif(SCIP_UPPER STREQUAL "AUTO")
    message(STATUS "SCIP not found.")
  else()
    message(FATAL_ERROR "SCIP not found.")
  endif()
endif()

if(SOPLEX_UPPER)
  if(NOT SOPLEX_FOUND)
    find_package(SOPLEX HINTS ${SOPLEX_DIR})
  endif()
  if(SOPLEX_FOUND)
    message(STATUS "Found SoPlex: ${SOPLEX_INCLUDE_DIRS} / ${SOPLEX_LIBRARIES}.")
  elseif(SOPLEX_UPPER STREQUAL "AUTO")
    message(STATUS "SoPlex not found.")
  else()
    message(FATAL_ERROR "SoPlex not found.")
  endif()
endif()

if(GUROBI_UPPER)
  find_package(GUROBI)
  if(GUROBI_FOUND)
    message(STATUS "Found Gurobi: ${GUROBI_INCLUDE_DIRS} / ${GUROBI_LIBRARIES}.")
  elseif(GUROBI_UPPER STREQUAL "AUTO")
    message(STATUS "Gurobi not found.")
  else()
    message(FATAL_ERROR "Gurobi not found.")
  endif()
endif()

if(GUROBI_UPPER AND GUROBI_FOUND)
  set(IPO_WITH_GUROBI ON)
else()
  set(IPO_WITH_GUROBI OFF)
endif()

if(SCIP_UPPER AND SCIP_FOUND)
  set(IPO_WITH_SCIP ON)
else()
  set(IPO_WITH_SCIP OFF)
endif()

if(SOPLEX_UPPER AND SOPLEX_FOUND)
  set(IPO_WITH_SOPLEX ON)
else()
  set(IPO_WITH_SOPLEX OFF)
endif()

set(IPO_WITH_DOUBLE_LP OFF)
set(IPO_WITH_DOUBLE_LP_GUROBI OFF)
set(IPO_WITH_DOUBLE_LP_SOPLEX OFF)
if(IPO_WITH_SOPLEX)
  set(IPO_WITH_DOUBLE_LP_SOPLEX ON)
  set(IPO_WITH_DOUBLE_LP ON)
elseif(IPO_WITH_GUROBI)
  set(IPO_WITH_DOUBLE_LP_GUROBI ON)
  set(IPO_WITH_DOUBLE_LP ON)
endif()

set(IPO_WITH_RATIONAL_LP OFF)
set(IPO_WITH_RATIONAL_LP_SOPLEX OFF)
if(IPO_WITH_SOPLEX)
  set(IPO_WITH_RATIONAL_LP_SOPLEX ON)
  set(IPO_WITH_RATIONAL_LP ON)
endif()

# RapidXML
find_package(RapidXML)
if(RapidXML_FOUND)
  message(STATUS "Found RapidXML ${RapidXML_VERSION}")
endif()

# Boost iostreams
find_package(Boost COMPONENTS iostreams HINTS ${BOOST_DIR})
if(Boost_IOSTREAMS_FOUND)
  message(STATUS "Found Boost iostreams ${Boost_VERSION}")
endif()

# Target for the IPO library.
add_library(ipo
  src/ipo/lp.cpp
  src/ipo/affine_hull.cpp
  src/ipo/arithmetic.cpp
  src/ipo/constraint.cpp
  src/ipo/dominant.cpp
  src/ipo/lu.cpp
  src/ipo/mip.cpp
  src/ipo/oracles.cpp
  src/ipo/oracles_forest.cpp
  src/ipo/oracles_polar.cpp
  src/ipo/oracles_radialcone.cpp
  src/ipo/polyhedron.cpp
  src/ipo/projection.cpp
  src/ipo/space.cpp
  src/ipo/submissive.cpp
  src/ipo/trust_region.cpp
)

# Link to/include SoPlex.
if(IPO_WITH_SOPLEX)
  target_link_libraries(ipo PRIVATE libsoplex-pic)
  # TODO: target_include_directories(ipo ... did not work here.
  include_directories(${SOPLEX_INCLUDE_DIRS})
endif()

# Whether to build the SCIP oracles.
if (IPO_WITH_SCIP)
  target_sources(ipo
    PRIVATE src/ipo/oracles_scip.cpp
  )
  target_link_libraries(ipo
    PRIVATE libscip
  )
  # TODO: target_include_directories(ipo ... did not work here.
  include_directories(
    ${SCIP_INCLUDE_DIRS}
  )
endif()

# Whether to build the Gurobi oracles.
if (IPO_WITH_GUROBI)
  target_sources(ipo
    PRIVATE src/ipo/oracles_gurobi.cpp
  )
  target_link_libraries(ipo
    PRIVATE ${GUROBI_LIBRARIES}
  )
  # TODO: target_include_directories(ipo ... did not work here.
  include_directories(
    ${GUROBI_INCLUDE_DIRS}
  )
endif()

# The ipo-scip binary.
if(IPO_WITH_SCIP)
  add_executable(ipo_scip src/mip/ipo_scip.cpp)
  target_link_libraries(ipo_scip PRIVATE IPO::IPO)
  set_target_properties(ipo_scip PROPERTIES OUTPUT_NAME ipo-scip)
endif()

# The ipo-gurobi binary.
if(IPO_WITH_GUROBI)
  add_executable(ipo_gurobi src/mip/ipo_gurobi.cpp)
  target_link_libraries(ipo_gurobi PRIVATE IPO::IPO)
  set_target_properties(ipo_gurobi PROPERTIES OUTPUT_NAME ipo-gurobi)
endif()

if(RapidXML_FOUND AND Boost_IOSTREAMS_FOUND)
  # The inverse-scip binary.
  if(IPO_WITH_SCIP)
    add_executable(inverse_scip src/mip/inverse_scip.cpp)
    target_link_libraries(inverse_scip PRIVATE RapidXML::RapidXML)
    target_link_libraries(inverse_scip PRIVATE Boost::iostreams)
    target_link_libraries(inverse_scip PRIVATE IPO::IPO)
    set_target_properties(inverse_scip PROPERTIES OUTPUT_NAME inverse-scip)
  endif()

  # The inverse-gurobi binary.
  if(IPO_WITH_GUROBI)
    add_executable(inverse_gurobi src/mip/inverse_gurobi.cpp)
    target_link_libraries(inverse_gurobi PRIVATE RapidXML::RapidXML)
    target_link_libraries(inverse_gurobi PRIVATE Boost::iostreams)
    target_link_libraries(inverse_gurobi PRIVATE IPO::IPO)
    set_target_properties(inverse_gurobi PROPERTIES OUTPUT_NAME inverse-gurobi)
  endif()
endif()

# Write compilation settings to ipo/config.hpp.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/ipo/config.hpp.in ${CMAKE_BINARY_DIR}/ipo/config.hpp @ONLY)

# Write export settings to ipo/export.hpp.
include(GenerateExportHeader)
generate_export_header(ipo EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/ipo/export.hpp)

# TODO: Why both things below?
# Add an alias so that library can be used inside the build tree.
add_library(IPO::IPO ALIAS ipo)

# This is required so that the exported target has the name IPO and not just ipo.
set_target_properties(ipo PROPERTIES EXPORT_NAME IPO)

# Set target properties.
target_include_directories(ipo
  PUBLIC
    $<INSTALL_INTERFACE:include> # <PREFIX>/include/ contains all installed headers.
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> # contains all regular headers.
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # contains all configured headers such as ipo/config.hpp and ipo/export.hpp.
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ipo/
)

target_compile_features(ipo PRIVATE cxx_auto_type)
target_compile_options(ipo PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall>)

if(IPO_WITH_RATIONAL)
  target_link_libraries(ipo
    INTERFACE
      GMP::GMP
  )
endif()

# Installation
include(GNUInstallDirs)

if(IPO_WITH_SCIP)
  install(TARGETS ipo_scip RUNTIME DESTINATION bin)
  if(TARGET inverse_scip)
    install(TARGETS inverse_scip RUNTIME DESTINATION bin)
  endif()
endif()
if(IPO_WITH_GUROBI)
  install(TARGETS ipo_gurobi RUNTIME DESTINATION bin)
  if(TARGET inverse_gurobi)
    install(TARGETS inverse_gurobi RUNTIME DESTINATION bin)
  endif()
endif()

install(TARGETS ipo
  EXPORT ipo-targets
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export the targets to a script.
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ipo/)
install(EXPORT ipo-targets
  FILE
    IPOTargets.cmake
    NAMESPACE IPO::
    DESTINATION ${INSTALL_CONFIGDIR}
)

include(CMakePackageConfigHelpers)

# Create IPOConfig.cmake
configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/IPOConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/IPOConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Create IPOConfigVersion.cmake.
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/IPOConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# Install IPOConfig.cmake and IPOConfigVersion.cmake.
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/IPOConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/IPOConfigVersion.cmake
  DESTINATION ${INSTALL_CONFIGDIR}
)

# Write exported targets.
export(EXPORT ipo-targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/IPOTargets.cmake NAMESPACE IPO::)

# Register package in user's package registry.
export(PACKAGE IPO)

add_subdirectory(doc)

enable_testing()
add_subdirectory(test)

