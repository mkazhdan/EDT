project(MarchingCubes)
cmake_minimum_required(VERSION 3.9)

add_executable(
  MarchingCubes
  MarchingCubes.cpp)

if(OpenMP_CXX_FOUND)
  target_link_libraries(MarchingCubes PUBLIC OpenMP::OpenMP_CXX)
endif()
