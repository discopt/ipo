include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src
  ${SoPlex_INCLUDE_DIRS}
)

if (SCIP_FOUND)
  string(TOUPPER ${SCIP_LPS} SCIP_LPI)
  add_definitions(-DWITH_SCIP -DWITH_SCIP_LPI_${SCIP_LPI})
  include_directories(${SCIP_INCLUDE_DIRS})
endif()

# ipo

add_executable(ipo-main
  ipo.cpp
)
set_target_properties(ipo-main PROPERTIES OUTPUT_NAME ipo)


target_link_libraries(ipo-main
  ${IPO_LIBRARIES}
)

if (SCIP_FOUND)
  target_link_libraries(ipo-main ${SCIP_LIBRARIES})
endif()

install(
  TARGETS ipo-main
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)

# paper-affine-hull

if(SCIP_FOUND)
  add_executable(ipo-paper-affine-hull
    paper_affine_hull.cpp
  )

  target_link_libraries(ipo-paper-affine-hull
    ${IPO_LIBRARIES}
  )

  if (SCIP_FOUND)
    target_link_libraries(ipo-paper-affine-hull ${SCIP_LIBRARIES})
  endif()
endif()

# paper-facets

if(SCIP_FOUND)
  add_executable(ipo-paper-facets
    paper_facets.cpp
  )

  target_link_libraries(ipo-paper-facets
    ${IPO_LIBRARIES}
  )

  if (SCIP_FOUND)
    target_link_libraries(ipo-paper-facets ${SCIP_LIBRARIES})
  endif()
endif()