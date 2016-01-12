# ${GMP_INCLUDE_DIRS} contains the paths to gmp.h (and gmpxx.h) if GMP is found.
# ${GMP_LIBRARIES} contains libgmp and libgmpxx if GMP is found.

FIND_PATH(GMP_INCLUDE_DIRS NAMES gmp.h gmpxx.h)

FIND_LIBRARY(GMP_LIBRARY_GMP gmp)
FIND_LIBRARY(GMP_LIBRARY_GMPXX gmpxx)

SET(GMP_LIBRARIES ${GMP_LIBRARY_GMP} ${GMP_LIBRARY_GMPXX})

IF (GMP_INCLUDE_DIRS AND GMP_LIBRARY)
   SET(GMP_FOUND TRUE)
ENDIF (GMP_INCLUDE_DIRS AND GMP_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDE_DIRS GMP_LIBRARIES)
