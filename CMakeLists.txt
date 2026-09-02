cmake_minimum_required(VERSION 3.15)
project(WelfordOnlineVariance LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler Warnings
if(MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

# Core Library
add_library(welford_lib
    src/accumulator.cpp
)

target_include_directories(welford_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# CLI Demo Target
add_executable(welford_demo apps/main.cpp)
target_link_libraries(welford_demo PRIVATE welford_lib)

# Unit Test Target
enable_testing()
add_executable(welford_tests tests/test_welford.cpp)
target_link_libraries(welford_tests PRIVATE welford_lib)
add_test(NAME WelfordUnitTest COMMAND welford_tests)

# Benchmark Target
add_executable(welford_benchmark benchmarks/bench_welford.cpp)
target_link_libraries(welford_benchmark PRIVATE welford_lib)
