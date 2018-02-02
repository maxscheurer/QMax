cmake_minimum_required(VERSION 3.6)
project(QMax)
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_LAZYTEN_DIR "path/to/lazyten" CACHE STRING "lazyten path")
set(CMAKE_LIBINT_DIR "path/to/libint" CACHE STRING "libint path")
set(CMAKE_PELIB_IFC_DIR "path/to/cppe" CACHE STRING "cppe path")

include_directories(
  ${CMAKE_LAZYTEN_DIR}/build/external/krims/src
  ${CMAKE_LAZYTEN_DIR}/external/krims/src/
  ${CMAKE_LAZYTEN_DIR}/external/krims/src/
  ${CMAKE_LAZYTEN_DIR}/build/src
  ${CMAKE_LAZYTEN_DIR}/src
  ${CMAKE_LIBINT_DIR}/include
  ${CMAKE_LIBINT_DIR}/libint2
  ${CMAKE_PELIB_IFC_DIR}
)

link_directories(
  ${CMAKE_LAZYTEN_DIR}/build/external/krims/src/krims
  ${CMAKE_LAZYTEN_DIR}/build/src/lazyten
  ${CMAKE_LIBINT_DIR}/lib
  ${CMAKE_PELIB_IFC_DIR}
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/include/eigen3/")
set(SOURCE_FILES main.cpp)
add_executable(QMax ${SOURCE_FILES})
target_link_libraries(QMax armadillo krims lazyten int2 lapack arpack cppe)
