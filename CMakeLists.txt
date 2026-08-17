cmake_minimum_required(VERSION 3.20)
project(butter VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(butter INTERFACE)
target_include_directories(butter INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_compile_features(butter INTERFACE cxx_std_20)

if(MSVC)
    target_compile_options(butter INTERFACE /utf-8)
endif()

option(BUTTER_BUILD_TESTS "Build tests" ON)
option(BUTTER_BUILD_EXAMPLES "Build examples" ON)
option(BUTTER_BUILD_BENCHMARKS "Build benchmarks" OFF)

if(BUTTER_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(BUTTER_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(BUTTER_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
