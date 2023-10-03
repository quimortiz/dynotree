# CMakeTemplate
# A template for CMake projects
project(kdtree)




# set c++17
set(CMAKE_CXX_STANDARD 17)



find_package(Eigen3 REQUIRED)



find_package(Boost REQUIRED COMPONENTS unit_test_framework)
add_subdirectory(deps/pybind11)



cmake_minimum_required(VERSION 3.0)
add_executable(main src/main.cpp)


target_include_directories(main PUBLIC include deps/nigh/src)
target_link_libraries(main PUBLIC Eigen3::Eigen Boost::unit_test_framework)


pybind11_add_module(pykdtree src/pykdtree.cpp)
target_include_directories(pykdtree PUBLIC include)
target_link_libraries(pykdtree PRIVATE Eigen3::Eigen)

set_property(TARGET pykdtree PROPERTY POSITION_INDEPENDENT_CODE ON)
# set_property(TARGET pykdtree PROPERTY POSITION_INDEPENDENT_CODE ON)



