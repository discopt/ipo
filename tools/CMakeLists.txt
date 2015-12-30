INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_BINARY_DIR}/src
  ${SCIP_INCLUDE_DIRS}
)

LINK_DIRECTORIES(
  ${SCIP_LIBRARY_DIRS}
)

ADD_EXECUTABLE(ipo-dimensions
  ipo_dimensions.cpp
  common.cpp
)

TARGET_LINK_LIBRARIES(ipo-dimensions
  ${IPO_LIBRARIES}
  ${SCIP_LIBRARIES}
)

INSTALL(
  TARGETS ipo-dimensions 
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)

ADD_EXECUTABLE(ipo-facets
  ipo_facets.cpp
  common.cpp
)

TARGET_LINK_LIBRARIES(ipo-facets
  ${IPO_LIBRARIES}
  ${SCIP_LIBRARIES}
)

INSTALL(
  TARGETS ipo-facets 
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)


ADD_EXECUTABLE(ipo-smallest-face
  ipo_smallest_face.cpp
  common.cpp
)

TARGET_LINK_LIBRARIES(ipo-smallest-face
  ${IPO_LIBRARIES}
  ${SCIP_LIBRARIES}
)

INSTALL(
  TARGETS ipo-smallest-face
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  COMPONENT bin
)

