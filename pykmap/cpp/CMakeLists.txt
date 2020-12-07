cmake_minimum_required(VERSION 3.1)
project(pykmapcore)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall -Wextra -pedantic " )

add_subdirectory(pybind11)
include_directories("${CMAKE_SOURCE_DIR}/../../")
include_directories("${CMAKE_SOURCE_DIR}/../../external/")
pybind11_add_module(pykmapcore main.cpp)
