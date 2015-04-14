# CMakeLists files in this project can
# refer to the root source directory of the project as ${HARP_SOURCE_DIR} and
# to the root binary directory of the project as ${HARP_BINARY_DIR}.
cmake_minimum_required(VERSION 2.6)
project(HARP CXX)

set(HARP_VERSION_MAJOR "0")
set(HARP_VERSION_MINOR "1")
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
