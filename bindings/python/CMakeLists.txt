# TODO

find_package(Eigen3 REQUIRED)

include(FetchContent)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11/
  GIT_TAG v2.11.1)

FetchContent_MakeAvailable(pybind11)

pybind11_add_module(pydynotree dynotree_pybindings.cpp)
target_link_libraries(pydynotree PRIVATE dynotree::dynotree Eigen3::Eigen)
set_property(TARGET pydynotree PROPERTY POSITION_INDEPENDENT_CODE ON)
