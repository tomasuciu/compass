# Top-level CMakeLists.txt
cmake_minimum_required(VERSION 3.8)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(compass,
    LANGUAGES CXX
    VERSION 0.0.0)

set(CMAKE_CXX_COMPILER clang++-10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ -v -lpython3.8")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libstdc++")

SET(CMAKE_FIND_PACKAGE_SORT_ORDER NATURAL)
SET(CMAKE_FIND_PACKAGE_SORT_DIRECTION DEC)


# Uses latest Eigen version, different from stable release installed locally
set(EIGEN3_DIR ${CMAKE_SOURCE_DIR}/thirdparty/eigen-master/build/)

find_package(Eigen3 3.3 REQUIRED)

message(STATUS "${EIGEN3_DIR}")
message(STATUS "${EIGEN3_INCLUDE_DIRS}")

add_subdirectory("src")

# Temporary; enables inclusion of Python.h
target_include_directories(compass PRIVATE /usr/include/python3.8)

# Numpy directory
target_include_directories(compass PRIVATE /usr/local/lib/python3.8/dist-packages/numpy/core/include)

# Enable matplotlib bindings
target_include_directories(compass PRIVATE ${CMAKE_SOURCE_DIR}/lib/matplotlib-cpp/)
