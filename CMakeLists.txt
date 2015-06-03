cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)

project(TCGA-Clusterer CXX)

file(GLOB TCGA-Clusterer_SRC
    "src/*.hpp"
    "src/*.cpp"
    "src/tests/*.hpp"
    "src/tests/*.cpp"
    "src/lodePNG/lodepng.cpp"
    "src/lodePNG/lodepng.h"
)

add_executable(TCGA-Clusterer ${TCGA-Clusterer_SRC})

SET(CMAKE_CXX_FLAGS "-std=c++17")