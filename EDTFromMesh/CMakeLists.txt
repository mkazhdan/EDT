project(EDTFromMesh)
cmake_minimum_required(VERSION 3.9)

add_executable(
  EDTFromMesh
  EDTFromMesh.cpp)

if(OpenMP_CXX_FOUND)
  target_link_libraries(EDTFromMesh PUBLIC OpenMP::OpenMP_CXX)
endif()
