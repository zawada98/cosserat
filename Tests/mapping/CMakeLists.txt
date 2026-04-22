# ── BeamContactMappingTest ────────────────────────────────────────────────────

# Use SOFA's helper macro if available (handles gtest_discover_tests,
# SOFA_BUILD_TESTS guard, and install rules automatically).
sofa_add_gtest(
    NAME    BeamContactMappingTest
    SOURCES BeamContactMappingTest.cpp
)

# If sofa_add_gtest is not available in your SOFA version, fall back to:
#
find_package(GTest REQUIRED)
add_executable(BeamContactMappingTest BeamContactMappingTest.cpp)
target_link_libraries(BeamContactMappingTest
     PRIVATE
         GTest::GTest
         GTest::Main
         Cosserat
         Sofa.Core
         Sofa.DefaultType
         Sofa.Helper
 )
include(GoogleTest)
gtest_discover_tests(BeamContactMappingTest)