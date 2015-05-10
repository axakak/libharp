cmake_minimum_required(VERSION 2.6)
project(LIBHARP CXX)

set(HARP_VERSION_MAJOR "1")
set(HARP_VERSION_MINOR "0")
set(HARP_VERSION_PATCH "0")

#TODO: expand debug build specs
#TODO: add release build specs
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 11)
add_compile_options(-Wall)

find_package(PkgConfig)
pkg_check_modules(YAMLCPP REQUIRED yaml-cpp>=0.5.1)

add_subdirectory(libharp)
add_subdirectory(util)
