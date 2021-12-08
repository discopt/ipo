cmake_minimum_required(VERSION 3.5)

project(IPO
   VERSION 2.0
   LANGUAGES CXX)
   
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)

option(GMP "Use GMP" ON)
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
if(GMP)
  find_package(GMPXX)
  set(IPO_WITH_GMP ${GMP_FOUND} AND ${GMPXX_FOUND})
else()
  set(IPO_WITH_GMP FALSE)
endif()

# Target for the IPO library.
add_library(ipo
  src/ipo/affine_hull.cpp
  src/ipo/arithmetic.cpp
  src/ipo/constraint.cpp
  src/ipo/dominant.cpp
  src/ipo/facet.cpp
  src/ipo/lp.cpp
  src/ipo/lu.cpp
  src/ipo/oracles.cpp
  src/ipo/mip.cpp
  src/ipo/oracles_forest.cpp
  src/ipo/polyhedron.cpp
  src/ipo/projection.cpp
  src/ipo/space.cpp
  src/ipo/submissive.cpp
)

# SoPlex
find_package(SOPLEX QUIET)
if(SOPLEX_FOUND)
  message(STATUS "Found SoPlex version ${SOPLEX_VERSION}.")
else()
  message(STATUS
    "Could NOT find SoPlex. Please set SOPLEX_DIR to point to build or installation directory.")
endif()
set(IPO_WITH_SOPLEX ${SOPLEX_FOUND})
if(SOPLEX_FOUND)
  target_link_libraries(ipo
    PRIVATE libsoplex-pic
  )
  # TODO: target_include_directories(ipo ... did not work here.
  include_directories(
    ${SOPLEX_INCLUDE_DIRS}
  )
endif()

# LP solver.
if (SOPLEX_FOUND)
  set(LP_DOUBLE TRUE)
  set(LP_DOUBLE_SOPLEX TRUE)
  set(LP_GMP TRUE)
  set(LP_GMP_SOPLEX TRUE)
endif()

# SCIP.
find_package(SCIP QUIET)
if(SCIP_FOUND)
  message(STATUS "Found SCIP version ${SCIP_VERSION}.")
else()
  message(STATUS
    "Could NOT find SCIP. Please set SCIP_DIR to point to build or installation directory.")
endif()
set(IPO_WITH_SCIP ${SCIP_FOUND})
if(SCIP_FOUND)
  message(STATUS "SCIP found for include dir ${SCIP_INCLUDE_DIRS}")
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

# Target for the ipo-mip binary.
add_executable(ipo_mip
  src/mip/main.cpp)
target_link_libraries(ipo_mip
  PRIVATE
    IPO::IPO
)
set_target_properties(ipo_mip PROPERTIES OUTPUT_NAME ipo-mip)

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

if(GMPXX_FOUND)
  target_link_libraries(ipo
    INTERFACE
      GMP::GMP
      GMP::GMPXX
  )
endif()

### Installation ###
include(GNUInstallDirs)

install(TARGETS ipo_mip
  RUNTIME
    DESTINATION bin
)

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

