include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src
  ${SOPLEX_INCLUDE_DIRS}
)

if (SCIP_FOUND)
  string(TOUPPER ${SCIP_LPS} SCIP_LPI)
  add_definitions(-DWITH_SCIP -DWITH_SCIP_LPI_${SCIP_LPI})
  include_directories(${SCIP_INCLUDE_DIRS})
endif()

# paper-affine-hull

add_executable(ipo-paper-affine-hull
  paper_affine_hull.cpp
)

target_link_libraries(ipo-paper-affine-hull
  ${IPO_LIBRARIES}
)

if (SCIP_FOUND)
  target_link_libraries(ipo-paper-affine-hull ${SCIP_LIBRARIES})
endif()


# paper-facets

add_executable(ipo-paper-facets
  paper_facets.cpp
)

target_link_libraries(ipo-paper-facets
  ${IPO_LIBRARIES}
)

if (SCIP_FOUND)
  target_link_libraries(ipo-paper-facets ${SCIP_LIBRARIES})
endif()

# sage-playground

add_executable(sage-playground
  sage_playground.cpp
)

target_link_libraries(sage-playground
  ${IPO_LIBRARIES}
)

if (SCIP_FOUND)
  target_link_libraries(sage-playground ${SCIP_LIBRARIES})
endif()
