include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/deps/ipo/ # for ipo/common.h
  ${SoPlex_INCLUDE_DIRS}
)

if (SCIP_FOUND)
  string(TOUPPER ${SCIP_LPS} SCIP_LPI)
  add_definitions(-DWITH_SCIP -DWITH_SCIP_LPI_${SCIP_LPI})
  include_directories(${SCIP_INCLUDE_DIRS})
endif()

# paper-affine-hull

if(SCIP_FOUND)
  add_executable(ipo-paper-affine-hull
    paper_affine_hull.cpp
  )

  target_link_libraries(ipo-paper-affine-hull
    ${IPO_LIBRARIES}
    ${SCIP_LIBRARIES}
  )

  install(TARGETS ipo-paper-affine-hull DESTINATION bin COMPONENT runtime)
endif()

# paper-facets

if(SCIP_FOUND)
  add_executable(ipo-paper-facets
    paper_facets.cpp
  )

  target_link_libraries(ipo-paper-facets
    ${IPO_LIBRARIES}
    ${SCIP_LIBRARIES}
  )

  install(TARGETS ipo-paper-facets DESTINATION bin COMPONENT runtime)
endif()
