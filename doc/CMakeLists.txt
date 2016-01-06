SET(PACKAGE_NAME ${PROJECT_NAME})
SET(PACKAGE_VERSION ${PROJECT_VERSION})
SET(abs_top_srcdir ${PROJECT_SOURCE_DIR})
SET(abs_top_builddir ${PROJECT_BINARY_DIR})

SET(IPO_DOC_SOURCE_BROWSER "NO" CACHE STRING "Include source into the doc (YES/NO).")

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in
  ${PROJECT_BINARY_DIR}/doc/Doxyfile
  @ONLY
)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/doc/mainpage.dox.in
  ${PROJECT_BINARY_DIR}/doc/mainpage.dox
  @ONLY
)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/doc/footer.html.in
  ${PROJECT_BINARY_DIR}/doc/footer.html
  @ONLY
)

# Copy doc from source (if exists)
IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/html AND 
    NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
  MESSAGE(STATUS "Copy doc from source tree")
  EXECUTE_PROCESS(
    COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/html ${CMAKE_CURRENT_BINARY_DIR}/html
    )
ENDIF()

IF(DOXYGEN_EXECUTABLE)
  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/)
  SET(GHOSTSCRIPT_OPTIONS -dNOPAUSE -dBATCH -q -dEPSCrop -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sDEVICE=pngalpha)
  ADD_CUSTOM_TARGET(html
    COMMAND ${CMAKE_COMMAND} -E remove_directory gen-images
    COMMAND ${CMAKE_COMMAND} -E make_directory gen-images
    COMMAND ${CMAKE_COMMAND} -E remove_directory html
#    COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/scripts/bib2dox.py ${CMAKE_CURRENT_SOURCE_DIR}/references.bib >references.dox
    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  SET_TARGET_PROPERTIES(html PROPERTIES PROJECT_LABEL BUILD_DOC)

  IF(UNIX)
    INSTALL(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
      DESTINATION share/doc/ipo/html
      COMPONENT html_documentation
    )
  ELSEIF(WIN32)
    INSTALL(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
      DESTINATION doc
      COMPONENT html_documentation
    )
  ENDIF()

ENDIF()

IF(WGET_FOUND)
ADD_CUSTOM_TARGET(update-external-tags
  COMMAND ${CMAKE_COMMAND} -E make_directory dl
  # COMMAND ${CMAKE_COMMAND} -E copy libstdc++.tag dl
  COMMAND ${WGET_EXECUTABLE} wget -P dl -N libstdc++.tag.tmp http://gcc.gnu.org/onlinedocs/libstdc++/latest-doxygen/libstdc++.tag
  COMMAND ${CMAKE_COMMAND} -E rename dl/libstdc++.tag libstdc++.tag
  COMMAND ${CMAKE_COMMAND} -E remove dl/libstdc++.tag
  COMMAND ${CMAKE_COMMAND} -E remove_directory dl
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
ENDIF()
