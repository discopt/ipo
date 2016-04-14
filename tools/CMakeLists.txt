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

