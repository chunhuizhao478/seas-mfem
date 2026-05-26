cmake_minimum_required(VERSION 3.12)
project(remesh_fault_cgal CXX)

# Build inside the cgal-61 conda env:
#   CG=/Users/chunhuizhao/miniforge/envs/cgal-61
#   $CG/bin/cmake -S . -B build -DCMAKE_PREFIX_PATH=$CG \
#       -DCMAKE_CXX_COMPILER=$CG/bin/clang++
#   $CG/bin/cmake --build build -j
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL REQUIRED)

add_executable(remesh_fault_cgal remesh_fault_cgal.cpp)
target_link_libraries(remesh_fault_cgal PRIVATE CGAL::CGAL)
