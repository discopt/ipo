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
  test_redundancy.cpp
  test_space.cpp)

if(SCIP_FOUND)
  target_sources(ipo_gtest
    PRIVATE test_oracles_scip.cpp
  )
endif()

# Add tests for non-exported functions only for static library.
#if(NOT SHARED)
#    target_sources(tu_gtest
#       PRIVATE
#          private.cpp)
#endif()

# Configure tu_gtest target.
target_compile_features(ipo_gtest PRIVATE cxx_auto_type)
target_link_libraries(ipo_gtest gtest_main IPO::ipo)
if(SCIP_FOUND)
  target_link_libraries(ipo_gtest libscip)
endif()

include(GoogleTest)
gtest_discover_tests(ipo_gtest)
