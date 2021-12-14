# Find Doxygen package.
find_package(Doxygen)

if(DOXYGEN_FOUND)
   set(DOXYGEN_PREDEFINED)
   if(IPO_RATIONAL)
      list(APPEND DOXYGEN_PREDEFINED IPO_RATIONAL)
   endif()
   if(IPO_RATIONAL_LP)
      list(APPEND DOXYGEN_PREDEFINED IPO_RATIONAL_LP)
   endif()
   set(DOXYGEN_EXTRA_PACKAGES amsmath amssymb mathtools)
   doxygen_add_docs(doc
      ${CMAKE_CURRENT_SOURCE_DIR}/README.md ${CMAKE_SOURCE_DIR}/include/ipo ${CMAKE_SOURCE_DIR}/src/ipo
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating doxygen documentation in ${CMAKE_CURRENT_BINARY_DIR}.")
else()
   message(WARNING "Could not find Doxygen. Generating a documentation is not supported.")
endif()
