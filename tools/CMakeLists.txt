include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/deps/ipo/ # for ipo/common.h
  ${SOPLEX_INCLUDE_DIRS}
)

if (SCIP_FOUND)
  include_directories(${SCIP_INCLUDE_DIRS})
endif()

# ipo-affine-hull

if(SCIP_FOUND)
  add_executable(ipo-affine-hull
    paper_affine_hull.cpp
  )

  target_link_libraries(ipo-affine-hull
    ${IPO_LIBRARIES}
  )
endif()

# ipo-facets

if(SCIP_FOUND)
  add_executable(ipo-facets
    paper_facets.cpp
  )

  target_link_libraries(ipo-facets
    ${IPO_LIBRARIES}
  )
endif()

# ipo-lp-relax

if(SCIP_FOUND)
  add_executable(ipo-lp-relax
    paper_lp_relax.cpp
  )

  target_link_libraries(ipo-lp-relax
    ${SCIP_LIBRARIES}
  )
endif()
