cmake_minimum_required(VERSION 3.6)
project(QMax)
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -L/home/max/Libraries/lazyten/build/external/krims/src/krims \
# -I/home/max/Libraries/lazyten/external/krims/src -L/home/max/Libraries/lazyten/build/src/lazyten -I/home/max/Libraries/lazyten/src \
# -I/usr/local/libint/2.3.0-beta.1/include -I/usr/local/libint/2.3.0-beta.1/include/libint2 -I/usr/include/eigen3 -L/usr/local/libint/2.3.0-beta.1/lib ")
# find_package(krims REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -I /Users/maxscheurer/Software/lazyten/build/external/krims/src -I/Users/maxscheurer/Software/lazyten/external/krims/src/ -I/usr/local/include/eigen3/ -I/Users/maxscheurer/Software/lazyten/build/src -I/Users/maxscheurer/Software/lazyten/src -I/Users/maxscheurer/Software/libint_final/include -I/Users/maxscheurer/Software/libint_final/libint2 -L/Users/maxscheurer/Software/libint_final/lib -L/Users/maxscheurer/Software/lazyten/build/external/krims/src/krims -L/Users/maxscheurer/Software/lazyten/build/src/lazyten")
set(SOURCE_FILES main.cpp)
add_executable(QMax ${SOURCE_FILES})
target_link_libraries(QMax armadillo krims lazyten int2 lapack arpack)
