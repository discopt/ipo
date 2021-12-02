find_package(GTest)
if(NOT GTEST_FOUND)
  find_package(GTestSources)
  if(NOT GTESTSOURCES_FOUND)
    find_package(GTestGit)
  endif()
endif()

# Target for the googletest.
add_executable(ipo_gtest
  test_affine_complement.cpp
  test_lu.cpp
  test_main.cpp
  test_oracles_forest.cpp
  test_redundancy.cpp
  test_space.cpp)

if(SCIP_FOUND)
  target_sources(ipo_gtest
    PRIVATE test_oracles_scip.cpp
  )
endif()

# Set target properties.
target_include_directories(ipo_gtest
  PUBLIC
    $<INSTALL_INTERFACE:include> # <PREFIX>/include/ contains all installed headers.
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> # contains all regular headers.
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # contains all configured headers such as ipo/config.hpp and ipo/export.hpp.
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ipo/
)

# Add tests for non-exported functions only for static library.
#if(NOT SHARED)
#    target_sources(tu_gtest
#       PRIVATE
#          private.cpp)
#endif()

# Configure tu_gtest target.
target_compile_features(ipo_gtest PRIVATE cxx_auto_type)
target_link_libraries(ipo_gtest gtest_main IPO::IPO)
if(SCIP_FOUND)
  target_link_libraries(ipo_gtest libscip)
endif()

include(GoogleTest)
gtest_discover_tests(ipo_gtest)
