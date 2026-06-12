cmake_minimum_required(VERSION 3.20)
project(safs_corefine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(corefine_pair corefine_pair.cpp)
add_executable(corefine_set  corefine_set.cpp)
add_executable(check_self_intersect check_self_intersect.cpp)
add_executable(mesh_volume   mesh_volume.cpp)
add_executable(autorefine_merged autorefine_merged.cpp)
add_executable(remesh_graded remesh_graded.cpp)

foreach(t corefine_pair corefine_set check_self_intersect mesh_volume autorefine_merged remesh_graded)
  target_link_libraries(${t} PRIVATE CGAL::CGAL Eigen3::Eigen)
  target_compile_options(${t} PRIVATE -O2 -Wall -Wno-deprecated-copy)
endforeach()
