include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src
  ${SOPLEX_INCLUDE_DIRS}
)

if (SCIP_FOUND)
  add_definitions(-DWITH_SCIP)
  include_directories(${SCIP_INCLUDE_DIRS})
endif()

# ipo-dimensions

add_executable(ipo-dimensions
  ipo_dimensions.cpp
  common.cpp
)

target_link_libraries(ipo-dimensions
  ${IPO_LIBRARIES}
)

if (SCIP_FOUND)
  target_link_libraries(ipo-dimensions ${SCIP_LIBRARIES})
endif()

install(
  TARGETS ipo-dimensions 
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)

# ipo-facets

add_executable(ipo-facets
  ipo_facets.cpp
  common.cpp
)

target_link_libraries(ipo-facets
  ${IPO_LIBRARIES}
)

install(
  TARGETS ipo-facets 
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)

# ipo-smallest-face

add_executable(ipo-smallest-face
  ipo_smallest_face.cpp
  common.cpp
)

target_link_libraries(ipo-smallest-face
  ${IPO_LIBRARIES}
)

install(
  TARGETS ipo-smallest-face
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)

if (SCIP_FOUND)
  target_link_libraries(ipo-dimensions ${SCIP_LIBRARIES})
  target_link_libraries(ipo-facets ${SCIP_LIBRARIES})
  target_link_libraries(ipo-smallest-face ${SCIP_LIBRARIES})
endif()

