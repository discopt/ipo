# ${GMP_INCLUDE_DIRS} contains the paths to gmp.h (and gmpxx.h) if GMP is found.
# ${GMP_LIBRARIES} contains libgmp and libgmpxx if GMP is found.

FIND_PATH(GMP_INCLUDE_DIRS NAMES gmp.h gmpxx.h)

FIND_LIBRARY(GMP_LIBRARY gmp)
FIND_LIBRARY(GMPXX_LIBRARY gmpxx)

SET(GMP_LIBRARIES ${GMP_LIBRARY} ${GMPXX_LIBRARY})

IF (GMP_INCLUDE_DIRS AND GMP_LIBRARY)
  SET(GMP_FOUND TRUE)
ENDIF (GMP_INCLUDE_DIRS AND GMP_LIBRARY)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_INCLUDE_DIRS GMP_LIBRARIES)
