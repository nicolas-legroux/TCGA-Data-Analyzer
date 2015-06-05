cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)

project(TCGA-Clusterer CXX)

set(Eigen_INCLUDE_DIR /home/nicolas/opt/Eigen/include/eigen3)

file(GLOB TCGA-Clusterer_SRC
    "src/*.hpp"
    "src/*.cpp"
    "src/tests/*.hpp"
    "src/tests/*.cpp"
    "src/lodePNG/lodepng.cpp"
    "src/lodePNG/lodepng.h"
)

include_directories(${Eigen_INCLUDE_DIR})

add_executable(TCGA-Clusterer ${TCGA-Clusterer_SRC})

SET(CMAKE_CXX_FLAGS "-O3 -std=c++14")