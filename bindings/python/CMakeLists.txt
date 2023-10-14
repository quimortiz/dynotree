# TODO

find_package(Eigen3 REQUIRED)

include(FetchContent)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11/
  GIT_TAG v2.11.1)
# it's much better to use a specific Git revision or Git tag for reproducibility

FetchContent_MakeAvailable(pybind11)

# add_subdirectory(../../deps/pybind11)
pybind11_add_module(pydynotree dynotree_pybindings.cpp)
# target_include_directories(dynotree PRIVATE include)
target_link_libraries(pydynotree PRIVATE dynotree::dynotree Eigen3::Eigen)
set_property(TARGET pydynotree PROPERTY POSITION_INDEPENDENT_CODE ON)
