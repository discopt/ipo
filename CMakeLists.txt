cmake_minimum_required(VERSION 3.5)

project(IPO VERSION 2.0 LANGUAGES CXX)
   
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)

option(DOUBLE "Compile algorithms with double arithmetic." ON)
option(RATIONAL "Compile algorithms with rational arithmetic; requires at least GMP and Boost 1.65" ON)
option(SCIP "Support SCIP as a MIP solver." ON)
option(GUROBI "Support Gurobi as a MIP solver." ON)
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

if(DOUBLE)
  set(IPO_DOUBLE TRUE)
else()
  set(IPO_DOUBLE FALSE)
endif()

set(IPO_DOUBLE_LP FALSE)
set(IPO_DOUBLE_LP_SOPLEX FALSE)
set(IPO_DOUBLE_LP_GUROBI FALSE)
set(IPO_DOUBLE_MIP FALSE)
set(IPO_DOUBLE_MIP_SCIP FALSE)
set(IPO_DOUBLE_MIP_GUROBI FALSE)
set(IPO_RATIONAL FALSE)
set(IPO_RATIONAL_LP FALSE)
set(IPO_RATIONAL_LP_SOPLEX FALSE)
set(IPO_RATIONAL_MIP FALSE)
set(IPO_RATIONAL_MIP_SCIP FALSE)
set(IPO_RATIONAL_MIP_GUROBI FALSE)

if(RATIONAL)
  find_package(GMP)
  if(${GMP_FOUND})
    find_package(Boost 1.65.0)
    if (${Boost_FOUND})
      set(IPO_RATIONAL TRUE)
    else()
      message(SEND_ERROR "ERROR: Option RATIONAL requires Boost library version 1.65.0 or greater.")
      set(IPO_RATIONAL FALSE)
    endif()
  else()
    message(SEND_ERROR "ERROR: Option RATIONAL requires GMP library.")
    set(IPO_RATIONAL FALSE)
  endif()
else()
  set(IPO_RATIONAL FALSE)
endif()

if(IPO_DOUBLE OR IPO_RATIONAL) # Search for SoPlex if requested.
  find_package(SOPLEX QUIET)
  if(SOPLEX_FOUND)
    message(STATUS "Found SoPlex version ${SOPLEX_VERSION} with include dir ${SOPLEX_INCLUDE_DIRS}.")
    if(IPO_DOUBLE)
      set(IPO_DOUBLE_LP TRUE)
      set(IPO_DOUBLE_LP_SOPLEX TRUE)
    endif()
    if(IPO_RATIONAL)
      set(IPO_RATIONAL_LP TRUE)
      set(IPO_RATIONAL_LP_SOPLEX TRUE)
    endif()
  else()
    message(STATUS "Could NOT find SoPlex. Please set SOPLEX_DIR to point to build or installation directory.")
  endif()
endif()

if(SCIP AND (IPO_DOUBLE OR IPO_RATIONAL_LP)) # Search for SCIP if requested.
  find_package(SCIP QUIET)
  if(SCIP_FOUND)
    message(STATUS "Found SCIP version ${SCIP_VERSION} with include dir ${SCIP_INCLUDE_DIRS}.")
    if(IPO_DOUBLE)
      set(IPO_DOUBLE_MIP TRUE)
      set(IPO_DOUBLE_MIP_SCIP TRUE)
    endif()
    if(IPO_RATIONAL_LP)
      set(IPO_RATIONAL_MIP TRUE)
      set(IPO_RATIONAL_MIP_SCIP TRUE)
    endif()
  else()
    message(STATUS "Could NOT find SCIP. Please set SCIP_DIR to point to build or installation directory.")
  endif()
endif()

if(GUROBI AND (IPO_DOUBLE OR IPO_RATIONAL_LP)) # Search for Gurobi if requested.
  find_package(GUROBI QUIET)
  if(GUROBI_FOUND)
    message(STATUS "Found Gurobi version ${GUROBI_VERSION} with include dir ${GUROBI_INCLUDE_DIRS} and libraries ${GUROBI_LIBRARIES}.")
    if(IPO_DOUBLE)
      set(IPO_DOUBLE_MIP TRUE)
      set(IPO_DOUBLE_MIP_GUROBI TRUE)
    endif()
    if(IPO_RATIONAL_LP)
      set(IPO_RATIONAL_MIP TRUE)
      set(IPO_RATIONAL_MIP_GUROBI TRUE)
    endif()
  else()
    message(STATUS "Could NOT find Gurobi. Please set GUROBI_DIR to point to build or installation directory.")
  endif()
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
  src/ipo/polyhedron.cpp
  src/ipo/projection.cpp
  src/ipo/space.cpp
  src/ipo/submissive.cpp
)

# Link to/include SoPlex.
if(IPO_DOUBLE_LP_SOPLEX OR IPO_RATIONAL_LP_SOPLEX)
  target_link_libraries(ipo PRIVATE libsoplex-pic)
  # TODO: target_include_directories(ipo ... did not work here.
  include_directories(${SOPLEX_INCLUDE_DIRS})
endif()

# Whether to build the SCIP oracles.
if (SCIP AND (IPO_DOUBLE_MIP_SCIP OR IPO_RATIONAL_MIP_SCIP))
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
if (GUROBI AND (IPO_DOUBLE_MIP_GUROBI OR IPO_RATIONAL_MIP_GUROBI))
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
if(IPO_DOUBLE_MIP_SCIP OR IPO_RATIONAL_MIP_SCIP)
  add_executable(ipo_scip src/mip/scip.cpp)
  target_link_libraries(ipo_scip PRIVATE IPO::IPO)
  set_target_properties(ipo_scip PROPERTIES OUTPUT_NAME ipo-scip)
endif()

# The ipo-gurobi binary.
if(IPO_DOUBLE_MIP_GUROBI OR IPO_RATIONAL_MIP_GUROBI)
  add_executable(ipo_gurobi src/mip/gurobi.cpp)
  target_link_libraries(ipo_gurobi PRIVATE IPO::IPO)
  set_target_properties(ipo_gurobi PROPERTIES OUTPUT_NAME ipo-gurobi)
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

if(IPO_RATIONAL)
  target_link_libraries(ipo
    INTERFACE
      GMP::GMP
  )
endif()

### Installation ###
include(GNUInstallDirs)

if(IPO_DOUBLE_MIP_SCIP OR IPO_RATIONAL_MIP_SCIP)
  install(TARGETS ipo_scip RUNTIME DESTINATION bin)
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

