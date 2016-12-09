# Tries to find ExactSCIP.
#
# Parameter Variables:
#
# ExactSCIP_ROOT_DIR
#   Set this variable to the SCIP source or install path
#   of the exact scip solver.
#
# Defines Variables:
#
# ExactSCIP_FOUND
#   True if ExactSCIP was found.
# ExactSCIP_EXECUTABLE
#   Path to the gs executable.
# ExactSCIP_VERSION
#   Version found.
#
# Author:
# 
# Matthias Walter <matthias@matthiaswalter.org>
#
# Distributed under the Boost Software License, Version 1.0.
# (See http://www.boost.org/LICENSE_1_0.txt)
# Tries to find ExactSCIP.
#
# Parameters:
#
# Defines Variables:
#
#  ExactSCIP_FOUND        - true if ExactSCIP was found.
#  ExactSCIP_EXECUTABLE   - path to the scip executable.
#
# Author:
# Matthias Walter <matthias@matthiaswalter.org>
#
# Distributed under the Boost Software License, Version 1.0.
# (See http://www.boost.org/LICENSE_1_0.txt)

find_program(ExactSCIP_EXECUTABLE NAMES scip PATHS ${ExactSCIP_ROOT_DIR}/bin NO_DEFAULT_PATH)
if(ExactSCIP_EXECUTABLE)
  execute_process(COMMAND ${ExactSCIP_EXECUTABLE} -c quit OUTPUT_VARIABLE _ExactSCIP_VERSION_STR OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE ".*SCIP version ([0-9.][0-9.]+) \\[.*$" "\\1" ExactSCIP_VERSION "${_ExactSCIP_VERSION_STR}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ExactSCIP REQUIRED_VARS ExactSCIP_EXECUTABLE VERSION_VAR ExactSCIP_VERSION)
