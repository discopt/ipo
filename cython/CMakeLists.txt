if(CYTHON_FOUND AND PYTHONINTERP_FOUND)

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CYTHON_IPO_DEBUG "True")
  else()
    set(CYTHON_IPO_DEBUG "False")
  endif()

  if(SCIP_FOUND)
    set(CYTHON_IPO_WITH_SCIP "True")
  else()
    set(CYTHON_IPO_WITH_SCIP "False")
  endif()

  if(ExactSCIP_FOUND OR EXACTSCIP_FOUND)
    set(CYTHON_IPO_WITH_EXACT_SCIP "True")
  else()
    set(CYTHON_IPO_WITH_EXACT_SCIP "False")
  endif()

  configure_file(setup.py.in setup.py)
  configure_file(IPO.pyx IPO.pyx COPYONLY)
  configure_file(config.pxi config.pxi)

  set(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build")
  set(TIMESTAMP "${OUTPUT}/ipo-cython.timestamp")

  add_custom_command(OUTPUT ${TIMESTAMP}
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py build
    COMMAND ${CMAKE_COMMAND} -E touch ${TIMESTAMP}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in ${CMAKE_CURRENT_SOURCE_DIR}/IPO.pyx
  )

  add_custom_target(ipo-cython ALL DEPENDS ${TIMESTAMP})

  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/setup.py install --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/)")

endif()
