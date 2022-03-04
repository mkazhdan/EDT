project(EDT)
cmake_minimum_required(VERSION 3.0)

# C++ code requires CXX 14
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

# code used everywhere
find_package(OpenMP)


# to have more warnings (mostly for GCC and CLANG)
if (NOT MSVC)
  include (CheckCXXCompilerFlag)
  set (WARNINGS "no-deprecated" "no-unused-result" "no-format")
  foreach(warning ${WARNINGS})
    check_cxx_compiler_flag("-W${warning}" CXX_SUPPORTS_WARNING)
    if(CXX_SUPPORTS_WARNING)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W${warning}")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -W${warning}")
    endif()
  endforeach()

  set(CFLAGS_RELEASE "unroll-loops" "fast-math")
  foreach(cflag ${CFLAGS_RELEASE})
    check_cxx_compiler_flag("-f${cflag}" CXX_SUPPORTS_FLAG)
    if(CXX_SUPPORTS_FLAG)
      set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -f${cflag}")
      set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -f${cflag}")
    endif()
  endforeach()

endif(NOT MSVC)

# to include header files related to this directory
include_directories(${EDT_SOURCE_DIR})

# make sure all libraries go in the same place
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${EDT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${EDT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EDT_BINARY_DIR}/bin)

# now compile libraries and executables
add_subdirectory(PNG)
add_subdirectory(EDTFromMesh)
add_subdirectory(EDTFromGrid)
add_subdirectory(MarchingCubes)
