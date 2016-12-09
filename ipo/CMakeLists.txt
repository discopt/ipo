cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME "Investigating Polyhedra by Oracles")
project(${PROJECT_NAME})

set(PROJECT_VERSION 1.0)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Doxygen)

set(IPO_SOURCES
  parameters.cpp
  polyhedron.cpp
  linear_constraint.cpp
  unique_vectors.cpp
  vectors.cpp
  rational.cpp
  spx_gmp.cpp
  space.cpp
  oracles.cpp
  projection.cpp
  cache_oracle.cpp
  vector_space_generators.cpp
  reconstruct.cpp
  vector_2d.cpp
  min_norm_2d.cpp
  mip.cpp
  facets.cpp
  affine_hull.cpp
  timer.cpp
  polar_lp.cpp
  smallest_face.cpp
  statistics_oracle.cpp
  external_oracle.cpp
  exactscip_oracle.cpp
  console_app.cpp
  parser.cpp
)

set(IPO_INSTALL_HEADERS
  affine_hull.h
  polyhedron.h
  cache_oracle.h
  console_app.h
  facets.h
  min_norm_2d.h
  mip.h
  oracles.h
  parameters.h
  parser.h
  polar_lp.h
  projection.h
  smallest_face.h
  space.h
  spx_gmp.h
  statistics_oracle.h
  timer.h
  unique_vectors.h
  exactscip_oracle.h
  rational.h
  reconstruct.h
  vector_2d.h
  vector_space_generators.h
)

set(IPO_OTHER_HEADERS
  common-devel.h
  common-release.h
  common-debug.h
  common-relwithdebinfo.h
  common-minsizerel.h
)

if(SCIP_FOUND)
  set(IPO_SOURCES ${IPO_SOURCES} scip_oracle.cpp)
  set(IPO_INSTALL_HEADERS ${IPO_INSTALL_HEADERS} scip_oracle.h scip_exception.hpp)
  string(TOUPPER ${SCIP_LPS} SCIP_LPI)
  add_definitions(-DWITH_SCIP -DWITH_SCIP_LPI_${SCIP_LPI})
endif()

if(ExactSCIP_FOUND)
  add_definitions(-DWITH_ExactSCIP)
  configure_file(exactscip_config.h.in ${CMAKE_BINARY_DIR}/exactscip_config.h)
  include_directories(${CMAKE_BINARY_DIR})
endif()

add_library(ipo-static STATIC ${IPO_SOURCES})
include_directories(ipo-static ${SoPlex_INCLUDE_DIRS})
target_link_libraries(ipo-static
  ${SoPlex_LIBRARIES}
  ${GMP_LIBRARIES}
)
set_target_properties(ipo-static PROPERTIES OUTPUT_NAME ipo)

add_library(ipo-shared SHARED ${IPO_SOURCES})
include_directories(ipo-shared ${SoPlex_INCLUDE_DIRS})
target_link_libraries(ipo-shared
  ${SoPlex_LIBRARIES}
  ${GMP_LIBRARIES}
)
set_target_properties(ipo-shared PROPERTIES OUTPUT_NAME ipo)

if(SCIP_FOUND)
  include_directories(${SCIP_INCLUDE_DIRS})
  target_link_libraries(ipo-static ${SCIP_LIBRARIES})
  target_link_libraries(ipo-shared ${SCIP_LIBRARIES})
endif()

install(TARGETS ipo-static   ARCHIVE DESTINATION lib   LIBRARY DESTINATION lib   COMPONENT library)
install(TARGETS ipo-shared   ARCHIVE DESTINATION lib   LIBRARY DESTINATION lib   COMPONENT library)

install(FILES ${IPO_INSTALL_HEADERS}                    DESTINATION include/ipo  COMPONENT headers)

install(FILES vectors-pub.h            RENAME vectors.h           DESTINATION include/ipo  COMPONENT headers)
install(FILES linear_constraint-pub.h  RENAME linear_constraint.h DESTINATION include/ipo  COMPONENT headers)

install(FILES common-debug.h           RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS Debug)
install(FILES common-release.h         RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS Release)
install(FILES common-relwithdebinfo.h  RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS RelWithDebInfo)
install(FILES common-minsizerel.h      RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS MinSizeRel)

