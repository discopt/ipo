CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

## Here comes the name of your project:

project(IPO)
# The version number.
set (IPO_VERSION_MAJOR 1)
set (IPO_VERSION_MINOR 0)

add_definitions("-std=gnu++11")


# Set default BUILD_TYPE:
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
  "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
  FORCE)
endif()

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

IF(UNIX)
  SET(INSTALL_BIN_DIR "bin" CACHE STRING
    "Subdir for installing the binaries")
  SET(INSTALL_DOC_DIR "share/doc/${PROJECT_NAME}"
    CACHE STRING
    "Subdir for installing the doc")
  SET(INSTALL_DATA_DIR "share/data/${PROJECT_NAME}"
    CACHE STRING
    "Subdir for installing the map data")
  SET(INSTALL_EXAMPLES_DIR "share/doc/${PROJECT_NAME}/examples"
    CACHE STRING
    "Subdir for installing the examples")
ELSE(UNIX)
  SET(INSTALL_BIN_DIR "." CACHE STRING
                  "Subdir for installing the binaries")
  SET(INSTALL_DOC_DIR "doc" CACHE STRING
                  "Subdir for installing the doc")
  SET(INSTALL_DATA_DIR "data" CACHE STRING
                  "Subdir for installing the map data")
  SET(INSTALL_EXAMPLES_DIR "examples" CACHE STRING
                       "Subdir for installing the examples")
ENDIF(UNIX)

## ZIMPL, SoPlex and SCIP

find_package(SoPlex REQUIRED)
find_package(SCIP)

## IPO

FIND_PATH(IPO_SOURCE_ROOT_DIR CMakeLists.txt
  PATHS ${CMAKE_SOURCE_DIR}/ipo ${CMAKE_SOURCE_DIR}/deps/ipo
  NO_DEFAULT_PATH
  DOC "Location of IPO source as a CMAKE subproject")

IF(EXISTS ${IPO_SOURCE_ROOT_DIR})
  ADD_SUBDIRECTORY(${IPO_SOURCE_ROOT_DIR} deps/ipo)
  SET(IPO_INCLUDE_DIRS
    ${IPO_SOURCE_ROOT_DIR}
    ${CMAKE_BINARY_DIR}/deps/ipo
  )
  SET(IPO_LIBRARIES ipo)
  UNSET(IPO_ROOT_DIR CACHE)
  UNSET(IPO_DIR CACHE)
  UNSET(IPO_INCLUDE_DIR CACHE)
  UNSET(IPO_LIBRARY CACHE)
ELSE()
  FIND_PACKAGE(IPO QUIET NO_MODULE)
  FIND_PACKAGE(IPO REQUIRED)
ENDIF()




## This line finds doxygen (for document creation)

find_package(Doxygen)
find_package(Ghostscript)
include(FindPythonInterp)

## These are the include directories used by the compiler.

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

ADD_SUBDIRECTORY(tools)

ENABLE_TESTING()

## The auxiliary doxygen files (.dox) should be placed in the 'doc'
## subdirectory. The next line includes the CMAKE config of that directory.

ADD_SUBDIRECTORY(doc)

