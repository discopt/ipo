include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src
  ${SCIP_INCLUDE_DIRS}
  ${SOPLEX_INCLUDE_DIRS}
)

add_executable(ipo-dimensions
  ipo_dimensions.cpp
  common.cpp
)

target_link_libraries(ipo-dimensions
  ${IPO_LIBRARIES}
  ${SCIP_LIBRARIES}
)

install(
  TARGETS ipo-dimensions 
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)

add_executable(ipo-facets
  ipo_facets.cpp
  common.cpp
)

target_link_libraries(ipo-facets
  ${IPO_LIBRARIES}
  ${SCIP_LIBRARIES}
  ${SOPLEX_LIBRARIES}
)

install(
  TARGETS ipo-facets 
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)


add_executable(ipo-smallest-face
  ipo_smallest_face.cpp
  common.cpp
)

target_link_libraries(ipo-smallest-face
  ${IPO_LIBRARIES}
  ${SCIP_LIBRARIES}
  ${SOPLEX_LIBRARIES}
)

install(
  TARGETS ipo-smallest-face
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)

