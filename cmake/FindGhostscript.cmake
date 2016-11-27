# Tries to find Ghostscript.
#
# Parameter Variables:
#
# Defines Variables:
#
# Ghostscript_FOUND
#   True if Ghostscript was found.
# Ghostscript_EXECUTABLE
#   Path to the gs executable.
# Ghostscript_VERSION
#   Version found.
#
# Author:
# 
# Matthias Walter <matthias@matthiaswalter.org>
#
# Distributed under the Boost Software License, Version 1.0.
# (See http://www.boost.org/LICENSE_1_0.txt)
# Tries to find Ghostscript.
#
# Parameters:
#
# Defines Variables:
#
#  Ghostscript_FOUND        - true if Ghostscript was found.
#  Ghostscript_EXECUTABLE   - path to the gs executable.
#
# Author:
# Matthias Walter <matthias@matthiaswalter.org>
#
# Distributed under the Boost Software License, Version 1.0.
# (See http://www.boost.org/LICENSE_1_0.txt)

find_program(Ghostscript_EXECUTABLE
  NAMES gs gswin32c
  PATHS "$ENV{ProgramFiles}/gs"
  PATH_SUFFIXES gs8.61/bin gs8.62/bin gs8.63/bin gs8.64/bin gs8.65/bin
  DOC "Ghostscript: PostScript and PDF language interpreter and previewer."
)
if(Ghostscript_EXECUTABLE)
  execute_process(COMMAND ${Ghostscript_EXECUTABLE} -version OUTPUT_VARIABLE _Ghostscript_VERSION_STR OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "^GPL Ghostscript ([0-9][0-9]*)\\.[0-9][0-9]* .*$" "\\1" Ghostscript_VERSION_MAJOR "${_Ghostscript_VERSION_STR}")
  string(REGEX REPLACE "^GPL Ghostscript [0-9][0-9]*\\.([0-9][0-9]*) .*$" "\\1" Ghostscript_VERSION_MINOR "${_Ghostscript_VERSION_STR}")
  set(Ghostscript_VERSION "${Ghostscript_VERSION_MAJOR}.${Ghostscript_VERSION_MINOR}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Ghostscript REQUIRED_VARS Ghostscript_EXECUTABLE VERSION_VAR Ghostscript_VERSION)
