cmake_minimum_required(VERSION 3.24)

project(RCS)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_PREFIX_PATH /home/user/lib/cmake)

find_package(MKL CONFIG REQUIRED)

add_executable(RCS main.cpp)

target_link_libraries(RCS PUBLIC MKL::MKL stdc++)

target_compile_options(RCS PUBLIC -march=native -DMKL_ILP64 -qmkl=parallel -fopenmp)