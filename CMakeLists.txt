cmake_minimum_required(VERSION 3.10)
project(UnitTest)

# Add include directories for Google Test headers and your source
include_directories(
    ./gtest/googletest/include
    ./gtest/googletest  # Needed for internal gtest source includes
    ./
)

# Add the test source
add_executable(unit_test ./src/matrix_math_test.cpp)

# Add the Google Test sources directly (if not building as a separate library)
add_subdirectory(gtest)

# Link with Google Test and pthread
target_link_libraries(unit_test
    gtest
    gtest_main
    pthread
)
