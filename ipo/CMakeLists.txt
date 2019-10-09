set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(IPO_SOURCES
  affine_hull.cpp
  cache_oracle.cpp
  console_app.cpp
  dominant.cpp
  exactscip_oracle.cpp
  external_oracle.cpp
  facets.cpp
  linear_constraint.cpp
  lp.cpp
  min_norm_2d.cpp
  mip_oracle.cpp
  oracles.cpp
  oracle_wrapper.cpp
  parameters.cpp
  parser.cpp
  polar_lp.cpp
  polyhedron.cpp
  projection.cpp
  rational.cpp
  reconstruct.cpp
  scip_oracle.cpp
  smallest_face.cpp
  soplex_reproduce.cpp
  space.cpp
  spx_gmp.cpp
  statistics_oracle.cpp
  submissive.cpp
  test.cpp
  timer.cpp
  unique_vectors.cpp
  vectors.cpp
  vector_2d.cpp
  vector_space_generators.cpp
)

set(IPO_INSTALL_HEADERS
  affine_hull.h
  polyhedron.h
  cache_oracle.h
  console_app.h
  facets.h
  min_norm_2d.h
  mip_oracle.h
  oracles.h
  parameters.h
  parser.h
  polar_lp.h
  projection.h
  smallest_face.h
  space.h
  spx_gmp.h
  submissive.h
  dominant.h
  statistics_oracle.h
  timer.h
  unique_vectors.h
  scip_oracle.h
  scip_exception.h
  exactscip_oracle.h
  rational.h
  reconstruct.h
  vector_2d.h
  vector_space_generators.h
  test.h
  lp.h
)

if(SCIP_FOUND)
#  string(TOUPPER ${SCIP_LPS} SCIP_LPI)
  add_definitions(-DWITH_SCIP)# -DWITH_SCIP_LPI_${SCIP_LPI})
  set(DEFINE_IPO_SCIP "#define IPO_WITH_SCIP\n#define IPO_WITH_SCIP_${SCIP_LPI}")
else()
  set(DEFINE_IPO_SCIP "#undef IPO_WITH_SCIP")
endif()

if(ExactSCIP_FOUND OR EXACTSCIP_FOUND)
  set(DEFINE_IPO_EXACT_SCIP "#define IPO_WITH_EXACT_SCIP\n#define IPO_EXACT_SCIP_PATH \"${ExactSCIP_EXECUTABLE}\"\n")
else()
  set(DEFINE_IPO_EXACT_SCIP "#undef IPO_WITH_EXACT_SCIP")
endif()

include_directories(${CMAKE_BINARY_DIR}/deps/ipo/) # for ipo/common.h
configure_file(common.h.in ${CMAKE_BINARY_DIR}/deps/ipo/common.h)
#configure_file(${CMAKE_BINARY_DIR}/common.h ${CMAKE_BINARY_DIR}/ipo/common.h COPYONLY)

add_library(ipo-static STATIC ${IPO_SOURCES})
include_directories(ipo-static ${SOPLEX_INCLUDE_DIRS})
target_link_libraries(ipo-static
  ${SOPLEX_LIBRARIES}
  ${GMP_LIBRARIES}
  ${GMPXX_LIBRARIES}
)
set_target_properties(ipo-static PROPERTIES OUTPUT_NAME ipo)

add_library(ipo-shared SHARED ${IPO_SOURCES})
include_directories(ipo-shared ${SOPLEX_INCLUDE_DIRS})


target_link_libraries(ipo-shared
  ${SOPLEX_PIC_LIBRARIES}
  ${GMP_LIBRARIES}
  ${GMPXX_LIBRARIES}
)
set_target_properties(ipo-shared PROPERTIES OUTPUT_NAME ipo)

if(SCIP_FOUND)
  include_directories(${SCIP_INCLUDE_DIRS})
#  target_link_libraries(ipo-static ${SCIP_LIBRARIES})
#target_link_libraries(ipo-shared ${SCIP_LIBRARIES})
endif()

install(TARGETS ipo-static   ARCHIVE DESTINATION lib   LIBRARY DESTINATION lib   COMPONENT library)
install(TARGETS ipo-shared   ARCHIVE DESTINATION lib   LIBRARY DESTINATION lib   COMPONENT library)

install(FILES ${IPO_INSTALL_HEADERS}                    DESTINATION include/ipo  COMPONENT headers)

install(FILES vectors-pub.h            RENAME vectors.h           DESTINATION include/ipo  COMPONENT headers)
install(FILES linear_constraint-pub.h  RENAME linear_constraint.h DESTINATION include/ipo  COMPONENT headers)

install(FILES ${CMAKE_BINARY_DIR}/deps/ipo/common.h DESTINATION include/ipo  COMPONENT headers)
