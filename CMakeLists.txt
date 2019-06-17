cmake_minimum_required(VERSION 3.5)

project(IPO
   VERSION 2.0
   LANGUAGES CXX)

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

# Find dependencies.
find_package(GMPXX)
set(IPO_WITH_GMP ${GMPXX_FOUND})

# Target for the IPO library.
add_library(ipo
   src/ipo/oracles.cpp
   src/ipo/space.cpp
)

# SCIP.
find_package(SCIP QUIET)
if(SCIP_FOUND)
   message(STATUS "Found SCIP version ${SCIP_VERSION}.")
else()
   message(STATUS
      "Could not find SCIP. Please set SCIP_DIR to point to build or installation directory.")
endif()
set(IPO_WITH_SCIP ${SCIP_FOUND})
if(SCIP_FOUND)
   target_sources(ipo
      PRIVATE src/ipo/oracles_scip.cpp
   )
   target_link_libraries(ipo
      PRIVATE libscip
   )
endif()

# Target for the ipo-mip binary.
add_executable(ipo_mip
   src/mip/main.cpp)
target_link_libraries(ipo_mip
   PRIVATE
      IPO::ipo
)
set_target_properties(ipo_mip PROPERTIES OUTPUT_NAME ipo-mip)

# Write compilation settings to ipo/config.hpp.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/ipo/config.hpp.in ${CMAKE_BINARY_DIR}/ipo/config.hpp @ONLY)

# Write export settings to ipo/export.hpp.
include(GenerateExportHeader)
generate_export_header(ipo EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/ipo/export.hpp)

# TODO: Why both things below?
# Add an alias so that library can be used inside the build tree.
add_library(IPO::ipo ALIAS ipo)

# This is required so that the exported target has the name IPO and not just ipo.
set_target_properties(ipo PROPERTIES EXPORT_NAME IPO)

# Hide non-exported symbols in shared library.
set_target_properties(ipo PROPERTIES CXX_VISIBILITY_PRESET hidden)
set_target_properties(ipo PROPERTIES VISIBILITY_INLINES_HIDDEN 1)

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

