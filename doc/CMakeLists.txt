# Find Doxygen package.
find_package(Doxygen)

if(DOXYGEN_FOUND)
   set(DOXYGEN_PREDEFINED)
   if(GMPXX_FOUND)
      list(APPEND DOXYGEN_PREDEFINED IPO_WITH_GMP)
   endif()
   if(SCIP_FOUND)
      list(APPEND DOXYGEN_PREDEFINED IPO_WITH_SCIP)
   endif()
   set(DOXYGEN_EXTRA_PACKAGES amsmath amssymb mathtools)
   doxygen_add_docs(doc
      ${CMAKE_CURRENT_SOURCE_DIR}/README.md ${CMAKE_SOURCE_DIR}/include/ipo ${CMAKE_SOURCE_DIR}/src/ipo
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating doxygen documentation in ${CMAKE_CURRENT_BINARY_DIR}.")
else()
   message(WARNING "Could not find Doxygen. Generating a documentation is not supported.")
endif()
