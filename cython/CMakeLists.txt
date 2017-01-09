if(Cython_FOUND)
  set_source_files_properties(IPO.pyx PROPERTIES CYTHON_IS_CXX 1)
  cython_add_module(IPO IPO.pyx)
endif()

