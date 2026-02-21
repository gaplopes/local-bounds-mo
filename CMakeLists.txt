cmake_minimum_required(VERSION 3.14)

project(LocalBounds
    VERSION 1.0.0
    DESCRIPTION "Local bounds library for multiobjective optimization"
    LANGUAGES CXX
)

# C++17 is required for std::optional and if constexpr
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Compiler warnings
add_compile_options(-Wall -Wextra -Wpedantic)

# Options
option(BUILD_TESTS "Build test programs" ON)
option(BUILD_EXAMPLES "Build example programs" ON)

# Header-only library
add_library(local_bounds INTERFACE)
target_include_directories(local_bounds INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(local_bounds INTERFACE cxx_std_17)

# Tests
if(BUILD_TESTS)
    enable_testing()

    add_executable(test_dominance tests/test_dominance.cpp)
    target_link_libraries(test_dominance PRIVATE local_bounds)
    add_test(NAME test_dominance COMMAND test_dominance)

    add_executable(test_bound_set tests/test_bound_set.cpp)
    target_link_libraries(test_bound_set PRIVATE local_bounds)
    add_test(NAME test_bound_set COMMAND test_bound_set)
endif()

# Examples
if(BUILD_EXAMPLES)
    add_executable(basic_usage examples/basic_usage.cpp)
    target_link_libraries(basic_usage PRIVATE local_bounds)
endif()

# Benchmark
option(BUILD_BENCHMARK "Build benchmark program" ON)
if(BUILD_BENCHMARK)
    add_executable(benchmark benchmark/benchmark.cpp)
    target_link_libraries(benchmark PRIVATE local_bounds)
endif()

# Installation
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS local_bounds
    EXPORT local_bounds-targets
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT local_bounds-targets
    FILE local_bounds-config.cmake
    NAMESPACE local_bounds::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/local_bounds
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/local_bounds-config-version.cmake"
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/local_bounds-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/local_bounds
)
