include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/deps/ipo/ # for ipo/common.h
  ${SOPLEX_INCLUDE_DIRS}
)

if (SCIP_FOUND)
#  string(TOUPPER ${SCIP_LPS} SCIP_LPI)
#  add_definitions(-DWITH_SCIP -DWITH_SCIP_LPI_${SCIP_LPI})
  include_directories(${SCIP_INCLUDE_DIRS})
endif()

# paper-affine-hull

if(SCIP_FOUND)
  add_executable(paper-affine-hull
    paper_affine_hull.cpp
  )

  target_link_libraries(paper-affine-hull
    ${IPO_LIBRARIES}
    ${SCIP_LIBRARIES}
#    ${SOPLEX_LIBRARIES}
  )
endif()

# paper-facets

if(SCIP_FOUND)
  add_executable(paper-facets
    paper_facets.cpp
  )

  target_link_libraries(paper-facets
    ${IPO_LIBRARIES}
    ${SCIP_LIBRARIES}
#    ${SOPLEX_LIBRARIES}
  )
endif()

# paper-lp-relax

if(SCIP_FOUND)
  add_executable(paper-lp-relax
    paper_lp_relax.cpp
  )

  target_link_libraries(paper-lp-relax
    ${SCIP_LIBRARIES}
#    ${SOPLEX_LIBRARIES}
  )
endif()
