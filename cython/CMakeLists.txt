if(Cython_FOUND)

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

endif()
