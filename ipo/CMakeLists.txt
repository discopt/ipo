cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME "Investigating Polyhedra by Oracles")
project(${PROJECT_NAME})

set(PROJECT_VERSION 1.0)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Doxygen)

set(IPO_SOURCES
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
  cpu_timer.cpp
  polar_lp.cpp
  smallest_face.cpp
  external_oracle.cpp
  console_app.cpp
  parser.cpp
)

set(IPO_INSTALL_HEADERS
  unique_vectors.h
  affine_hull.h
  facets.h
  mixed_integer_program.h
  polar_lp.h
  rows.h
  smallest_face.h
  unique_rational_vectors.h
  vector_space_generators.h
  cpu_timer.h
  min_norm_2d.h
  oracles.h
  projection.h
  cache_oracle.h
  reconstruct.h
  spx_gmp.h
  vector_2d.h
  console_app.h
  parser.h
  space.h
)

set(IPO_OTHER_HEADERS
  common-devel.h
  common-release.h
  common-debug.h
  common-relwithdebinfo.h
  common-minsizerel.h
)

if (SCIP_FOUND)
  set(IPO_SOURCES ${IPO_SOURCES} scip_oracles.cpp)
  set(IPO_INSTALL_HEADERS ${IPO_INSTALL_HEADERS} scip_oracles.h scip_exception.hpp)
  string(TOUPPER ${SCIP_LPS} SCIP_LPI)
  add_definitions(-DWITH_SCIP -DWITH_SCIP_LPI_${SCIP_LPI})
endif()

add_library(ipo ${IPO_SOURCES})
include_directories(ipo ${SOPLEX_INCLUDE_DIRS})
target_link_libraries(ipo
  ${SOPLEX_LIBRARIES}
  ${GMP_LIBRARIES}
)

if (SCIP_FOUND)
  include_directories(${SCIP_INCLUDE_DIRS})
  target_link_libraries(ipo ${SCIP_LIBRARIES})
endif()

install(TARGETS ipo   ARCHIVE DESTINATION lib   LIBRARY DESTINATION lib   COMPONENT library)

install(FILES ${IPO_INSTALL_HEADERS}                    DESTINATION include/ipo  COMPONENT headers)

install(FILES vectors-pub.h            RENAME vectors.h           DESTINATION include/ipo  COMPONENT headers)
install(FILES linear_constraint-pub.h  RENAME linear_constraint.h DESTINATION include/ipo  COMPONENT headers)

install(FILES common-debug.h           RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS Debug)
install(FILES common-release.h         RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS Release)
install(FILES common-relwithdebinfo.h  RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS RelWithDebInfo)
install(FILES common-minsizerel.h      RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS MinSizeRel)

