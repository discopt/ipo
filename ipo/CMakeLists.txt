cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME "Investigating Polyhedra by Oracles")
project(${PROJECT_NAME})

set(PROJECT_VERSION 1.0)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Doxygen)

set(IPO_SOURCES
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

set(IPO_INSTALL_HEADERS
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

set(IPO_OTHER_HEADERS
  ipo-devel.h
  ipo-release.h
  ipo-debug.h
  ipo-relwithdebinfo.h
  ipo-minsizerel.h
)

include_directories(ipo
  ${SCIP_INCLUDE_DIRS}
  ${SOPLEX_INCLUDE_DIRS}
)

add_library(ipo ${IPO_SOURCES})

target_link_libraries(ipo
  ${SCIP_LIBRARIES}
  ${SOPLEX_LIBRARIES}
  ${GMP_LIBRARIES}
)

install(TARGETS ipo   ARCHIVE DESTINATION lib   LIBRARY DESTINATION lib   COMPONENT library)

install(FILES ${IPO_INSTALL_HEADERS}              DESTINATION include/ipo  COMPONENT headers)
install(FILES ipo-debug.h           RENAME ipo.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS Debug)
install(FILES ipo-release.h         RENAME ipo.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS Release)
install(FILES ipo-relwithdebinfo.h  RENAME ipo.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS RelWithDebInfo)
install(FILES ipo-minsizerel.h      RENAME ipo.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS MinSizeRel)

