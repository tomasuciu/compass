# Top-level CMakeLists.txt
cmake_minimum_required(VERSION 3.8)

project(compass
    LANGUAGES CXX
    VERSION 0.0.0)

string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_CAPS)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

set(CMAKE_CXX_COMPILER clang++-10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++") #-v -lpython3.8")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libstdc++")

#SET(CMAKE_FIND_PACKAGE_SORT_ORDER NATURAL)
#SET(CMAKE_FIND_PACKAGE_SORT_DIRECTION DEC)


# Uses latest Eigen version, different from stable release installed locally
set(EIGEN3_DIR ${CMAKE_SOURCE_DIR}/thirdparty/eigen-master/Eigen)

find_package(Eigen3 3.3 REQUIRED)

message(STATUS "${EIGEN3_DIR}")
message(STATUS "${EIGEN3_INCLUDE_DIRS}")
message(STATUS "${EIGEN3_FOUND}")

if(NOT EIGEN3_FOUND)
  LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
  find_package(Eigen3 REQUIRED)
endif()

# Note that eigen 3.2 (on Ubuntu Wily) only provides EIGEN3_INCLUDE_DIR,
# not EIGEN3_INCLUDE_DIRS, so we have to set the latter from the former.
if(NOT EIGEN3_INCLUDE_DIRS)
  set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
endif()

# Necessary for Ubuntu 16.04's god-awful FindEigen3.cmake.
if((NOT (DEFINED EIGEN3_VERSION)) AND (DEFINED EIGEN3_VERSION_STRING))
  set(EIGEN3_VERSION ${EIGEN3_VERSION_STRING})
endif()

if(${EIGEN3_VERSION} VERSION_EQUAL "3.3.6")
  message(WARNING "Eigen3 version ${EIGEN3_VERSION} found in ${EIGEN3_INCLUDE_DIRS},"
                  "but this version has a [bug](http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1643)")
endif(${EIGEN3_VERSION} VERSION_EQUAL "3.3.6")


add_library(${PROJECT_NAME} INTERFACE)

message(STATUS "${PROJECT_SOURCE_DIR}")

target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)

# Add Eigen interface dependency, depending on available cmake info
if(TARGET Eigen3::Eigen)
  target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
  set(Eigen3_DEPENDENCY "find_dependency(Eigen3 ${Eigen3_VERSION})")
else(TARGET Eigen3::Eigen)
  target_include_directories(${PROJECT_NAME} SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIRS})
endif(TARGET Eigen3::Eigen)

set(config_install_dir "share/${PROJECT_NAME}/cmake")
set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")

# Targets:
install(
  TARGETS  ${PROJECT_NAME}
  EXPORT   ${PROJECT_NAME}Targets
  INCLUDES DESTINATION include
)

install(
  EXPORT      ${PROJECT_NAME}Targets
  NAMESPACE   ${PROJECT_NAME_CAPS}::
  DESTINATION "${config_install_dir}"
)

export(
  TARGETS ${PROJECT_NAME}
  NAMESPACE ${PROJECT_NAME_CAPS}::
  FILE ${PROJECT_NAME}Targets.cmake
)

export(PACKAGE ${PROJECT_NAME})

## Configuration

#
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${generated_dir}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${config_install_dir}"
)

# Configure '<PROJECT-NAME>ConfigVersion.cmake'
write_basic_package_version_file(
  "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION 0.0.0
  COMPATIBILITY AnyNewerVersion
)

# Config
install(
  FILES       "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake"
              "${generated_dir}/${PROJECT_NAME}Config.cmake"
  DESTINATION "${config_install_dir}"
)

install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/Compass"
    DESTINATION include
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

