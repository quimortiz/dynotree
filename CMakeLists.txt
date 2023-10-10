# CMakeTemplate
# A template for CMake projects
project(kdtree)




# set c++17
set(CMAKE_CXX_STANDARD 17)

cmake_minimum_required(VERSION 3.0)
find_package(Eigen3 REQUIRED)


option(BUILD_PYTHON ON)
option(BUILD_TESTING OFF)

message(STATUS "BUILD_TESTING: ${BUILD_TESTING}")
message(STATUS "BUILD_PYTHON: ${BUILD_PYTHON}")


IF (BUILD_PYTHON)
  add_subdirectory(deps/pybind11)
  pybind11_add_module(dynotree src/pykdtree.cpp)
  target_include_directories(dynotree PRIVATE include )
  target_link_libraries(dynotree PRIVATE Eigen3::Eigen)
  set_property(TARGET dynotree PROPERTY POSITION_INDEPENDENT_CODE ON)
ENDIF (BUILD_PYTHON)


if (BUILD_TESTING)
  find_package(ompl 1.6 REQUIRED)
  message(STATUS "OMPL_FOUND: ${OMPL_FOUND}")
  message(STATUS "OMPL_INCLUDE_DIRS: ${OMPL_INCLUDE_DIRS}")
  message(STATUS "OMPL_LIBRARIES: ${OMPL_LIBRARIES}")
  message(STATUS "OMPL_VERSION: ${OMPL_VERSION}")



  find_package(Boost REQUIRED COMPONENTS unit_test_framework system filesystem serialization)



  add_executable(main src/main.cpp)


  target_include_directories(main PRIVATE include deps/nigh/src ${OMPL_INCLUDE_DIRS})
  target_link_libraries(main PRIVATE Eigen3::Eigen Boost::unit_test_framework Boost::filesystem ${OMPL_LIBRARIES})
ENDIF (BUILD_TESTING)
