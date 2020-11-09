cmake_minimum_required(VERSION 3.16)
project(fast_sine)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Ofast -ftree-vectorize -march=native -fno-associative-math ")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Ofast -march=native ")

find_package(PkgConfig REQUIRED)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

find_package(Boost REQUIRED COMPONENTS Timer)

add_executable(fast_sine main.cpp)
#target_link_libraries(fast_sine Boost::timer)
target_link_libraries(fast_sine Boost::timer PkgConfig::mpfr)