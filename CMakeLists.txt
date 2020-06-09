# Top-level CMakeLists.txt
cmake_minimum_required(VERSION 3.8)
project(compass,
    LANGUAGES CXX
    VERSION 0.0.0)

set(CMAKE_CXX_COMPILER clang++-10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ -v")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libstdc++")

add_subdirectory("src")

# Temporary; enables inclusion of Python.h
target_include_directories(compass PRIVATE /usr/include/python3.8)

# Numpy directory
target_include_directories(compass PRIVATE /usr/local/lib/python3.8/dist-packages/numpy/core/include)

# Enable matplotlib bindings
target_include_directories(compass PRIVATE ${CMAKE_SOURCE_DIR}/lib/matplotlib-cpp/)
