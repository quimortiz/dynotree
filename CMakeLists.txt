# CMakeTemplate A template for CMake projects

# Build system is based on https://dominikberner.ch/cmake-interface-lib/
# https://github.com/bernedom/SI

cmake_minimum_required(VERSION 3.12)

project(
  dynotree
  VERSION 0.0.4
  DESCRIPTION "A dynamic nearest neighbor search library"
  LANGUAGES CXX
  HOMEPAGE_URL "https://github.com/quimortiz/dyn_kdtree/")

include(GNUInstallDirs)

# Is it required to find eigen here? find_package(Eigen3 REQUIRED)

add_library(dynotree INTERFACE)
add_library(dynotree::dynotree ALIAS dynotree)

target_include_directories(
  dynotree INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                     $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(dynotree INTERFACE cxx_std_17)

# add_library(dynobench::dynobench ALIAS dynobench)
option(BUILD_PYTHON_BINDINGS OFF)
option(BUILD_TESTING OFF)
option(BUILD_EXAMPLES OFF)

#
message(STATUS "BUILD_TESTING: ${BUILD_TESTING}")
message(STATUS "BUILD_PYTHON_BINDINGS ${BUILD_PYTHON_BINDINGS}")
message(STATUS "BUILD_EXAMPLES: ${BUILD_EXAMPLES}")
#

#

# endif(BUILD_TESTING)

install(
  TARGETS dynotree
  EXPORT dynotree_Targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "dynotreeConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/dynotreeConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/dynotreeConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/dynotree/cmake)

install(
  EXPORT dynotree_Targets
  FILE dynotreeTargets.cmake
  NAMESPACE dynotree::
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/dynotree/cmake)

install(FILES "${PROJECT_BINARY_DIR}/dynotreeConfig.cmake"
              "${PROJECT_BINARY_DIR}/dynotreeConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/dynotree/cmake)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/dynotree DESTINATION include)

option(ENABLE_TEST_COVERAGE OFF)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(test/cpp)
endif()

if(BUILD_PYTHON_BINDINGS)
  add_subdirectory(bindings/python)
endif()
