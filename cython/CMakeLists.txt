if(Cython_FOUND)
  
  configure_file(setup.py.in setup.py)
  configure_file(IPO.pyx IPO.pyx COPYONLY)

endif()



  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
      message(STATUS "${_variableName}=${${_variableName}}")
      endforeach()

