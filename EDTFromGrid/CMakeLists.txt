project(EDTFromGrid)
cmake_minimum_required(VERSION 3.9)

find_package(JPEG REQUIRED)
find_package(ZLIB REQUIRED)

include_directories(
  ${JPEG_INCLUDE_DIR}
  ${ZLIB_INCLUDE_DIRS})

add_executable(
  EDTFromGrid
  EDTFromGrid.cpp)

target_link_libraries(
  EDTFromGrid
  PUBLIC
  ${JPEG_LIBRARIES}
  mypng
  ${ZLIB_LIBRARIES})

if(OpenMP_CXX_FOUND)
  target_link_libraries(EDTFromGrid PUBLIC OpenMP::OpenMP_CXX)
endif()
