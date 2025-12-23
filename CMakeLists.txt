cmake_minimum_required(VERSION 3.16)
project(algeng-matrix-multiplication LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(bench
        main.cpp
        strassen.cpp
        alphaevolve.cpp
        naive.cpp
)

target_compile_options(bench PRIVATE -O2 -march=native -DNDEBUG -Wall -Wextra)
