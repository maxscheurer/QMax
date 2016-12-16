cmake_minimum_required(VERSION 3.6)
project(QMax)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -L/home/max/Libraries/linalgwrap/build/src -I/home/max/Libraries/linalgwrap/src -I/usr/local/libint/2.3.0-beta.1/include -I/usr/local/libint/2.3.0-beta.1/include/libint2 -I/usr/include/eigen3 -L/usr/local/libint/2.3.0-beta.1/lib ")
set(SOURCE_FILES main.cpp)
add_executable(QMax ${SOURCE_FILES})
target_link_libraries(QMax int2)