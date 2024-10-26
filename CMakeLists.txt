cmake_minimum_required(VERSION 3.22)
project(cmatrix)

set(CMAKE_CXX_STANDARD 17)

add_executable(cmatrix test.cpp
        matrix.cpp
)
