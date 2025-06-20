cmake_minimum_required(VERSION 3.10)
find_package(Eigen3 CONFIG REQUIRED)
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
find_package(Eigen3 CONFIG REQUIRED)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(unit_test
    gtest
    gtest_main
    pthread
    Eigen3::Eigen
)
