if(Cython_FOUND)
  
  configure_file(setup.py.in setup.py)
  configure_file(IPO.pyx IPO.pyx COPYONLY)

endif()
