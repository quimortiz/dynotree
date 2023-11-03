find_package(ompl 1.6 REQUIRED)
find_package(Eigen3 REQUIRED)

include(FetchContent)
FetchContent_Declare(nigh GIT_REPOSITORY https://github.com/quimortiz/nigh/)

FetchContent_MakeAvailable(nigh)

message(STATUS "OMPL_FOUND: ${OMPL_FOUND}")
message(STATUS "OMPL_INCLUDE_DIRS: ${OMPL_INCLUDE_DIRS}")
message(STATUS "OMPL_LIBRARIES: ${OMPL_LIBRARIES}")
message(STATUS "OMPL_VERSION: ${OMPL_VERSION}")

find_package(Boost REQUIRED COMPONENTS unit_test_framework system filesystem
                                       serialization)

add_executable(main main.cpp)

target_include_directories(main PRIVATE ${OMPL_INCLUDE_DIRS})

target_link_libraries(
  main PRIVATE dynotree::dynotree Eigen3::Eigen Boost::unit_test_framework
               nigh::nigh Boost::filesystem ${OMPL_LIBRARIES})

add_test(NAME test COMMAND main)
