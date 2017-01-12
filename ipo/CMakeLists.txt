set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

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
  oracle_wrapper.cpp
  timer.cpp
  polar_lp.cpp
  smallest_face.cpp
  statistics_oracle.cpp
  external_oracle.cpp
  exactscip_oracle.cpp
  scip_oracle.cpp
  console_app.cpp
  parser.cpp
  python_wrapper.cpp
  test.cpp
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
  scip_oracle.h
  scip_exception.h
  exactscip_oracle.h
  rational.h
  reconstruct.h
  vector_2d.h
  vector_space_generators.h
  python_wrapper.h
  test.h
)

if(SCIP_FOUND)
  string(TOUPPER ${SCIP_LPS} SCIP_LPI)
  add_definitions(-DWITH_SCIP -DWITH_SCIP_LPI_${SCIP_LPI})
  set(DEFINE_IPO_SCIP "#define IPO_WITH_SCIP\n#define IPO_WITH_SCIP_${SCIP_LPI}")
else()
  set(DEFINE_IPO_SCIP "#undef IPO_WITH_SCIP")
endif()

if(ExactSCIP_FOUND OR EXACTSCIP_FOUND)
  set(DEFINE_IPO_EXACT_SCIP "#define IPO_WITH_EXACT_SCIP\n#define IPO_EXACT_SCIP_PATH \"${ExactSCIP_EXECUTABLE}\"\n")
else()
  set(DEFINE_IPO_EXACT_SCIP "#undef IPO_WITH_EXACT_SCIP")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(DEFINE_IPO_DEBUG "#define IPO_DEBUG")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(DEFINE_IPO_DEBUG "#define IPO_DEBUG")
elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
  set(DEFINE_IPO_DEBUG "#undef IPO_DEBUG")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(DEFINE_IPO_DEBUG "#undef IPO_DEBUG")
endif()

include_directories(${CMAKE_BINARY_DIR}/deps/ipo/) # for ipo/common.h
configure_file(common.h.in ${CMAKE_BINARY_DIR}/deps/ipo/common.h)
#configure_file(${CMAKE_BINARY_DIR}/common.h ${CMAKE_BINARY_DIR}/ipo/common.h COPYONLY)

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

install(FILES ${CMAKE_BINARY_DIR}/common.h DESTINATION include/ipo  COMPONENT headers)

# install(FILES common-debug.h           RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS Debug)
# install(FILES common-release.h         RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS Release)
# install(FILES common-relwithdebinfo.h  RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS RelWithDebInfo)
# install(FILES common-minsizerel.h      RENAME common.h  DESTINATION include/ipo  COMPONENT headers CONFIGURATIONS MinSizeRel)

