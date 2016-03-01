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
  mixed_integer_program.cpp
  facets.cpp
  affine_hull.cpp
  cpu_timer.cpp
  polar_lp.cpp
  smallest_face.cpp
  wrapper_oracle.cpp
  console_app.cpp
  parser.cpp
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
  spx_gmp.h
  vector_2d.h
  wrapper_oracle.h
  console_app.h
  parser.h
)

set(IPO_OTHER_HEADERS
  ipo-devel.h
  ipo-release.h
  ipo-debug.h
  ipo-relwithdebinfo.h
  ipo-minsizerel.h
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

install(FILES ${IPO_INSTALL_HEADERS}              DESTINATION include/ipo  COMPONENT headers)
install(FILES ipo-debug.h           RENAME ipo.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS Debug)
install(FILES ipo-release.h         RENAME ipo.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS Release)
install(FILES ipo-relwithdebinfo.h  RENAME ipo.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS RelWithDebInfo)
install(FILES ipo-minsizerel.h      RENAME ipo.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS MinSizeRel)

