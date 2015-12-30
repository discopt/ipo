CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(PROJECT_NAME "Investigating Polyhedra by Oracle")
PROJECT(${PROJECT_NAME})

INCLUDE(FindGMP)

SET(PROJECT_VERSION 1.0)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

FIND_PACKAGE(Doxygen)

SET(IPO_SOURCES
  spx_gmp.cpp
  vector_space_generators.cpp
  reconstruct.cpp
  vector_2d.cpp
  min_norm_2d.cpp
  unique_rational_vectors.cpp
  rows.cpp
  oracles.cpp
  scip_oracles.cpp
  mixed_integer_program.cpp
  facets.cpp
  affine_hull.cpp
  cpu_timer.cpp
  polar_lp.cpp
  smallest_face.cpp
  wrapper_oracle.cpp
)

SET(IPO_INSTALL_HEADERS
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
  reconstruct.h
  scip_oracles.h
  spx_gmp.h
  vector_2d.h
  wrapper_oracle.h
)

ADD_LIBRARY(ipo ${IPO_SOURCES})

TARGET_LINK_LIBRARIES(
  ipo
  gmpxx
)

INSTALL(
  TARGETS ipo
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  COMPONENT library
)

INSTALL(
  FILES ${IPO_INSTALL_HEADERS}
  DESTINATION include/ipo
  COMPONENT headers
)
