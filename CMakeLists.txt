cmake_minimum_required(VERSION 3.5.1)

project(mandelbrot)

include(ExternalProject)

set(CMAKE_CXX_STANDARD 14)

################################################################################
# mandelbrot target

set(SOURCE_FILES main.cpp)
add_executable(mandelbrot ${SOURCE_FILES})

set_target_properties(
        mandelbrot
        PROPERTIES COMPILE_FLAGS "-msse4.2 -DSG14_USE_INT128 -Wall -Wextra"
)

################################################################################
# fixed_point

include(ExternalProject)

ExternalProject_Add(
        fixed_point
        GIT_REPOSITORY "https://github.com/johnmcfarlane/fixed_point.git"
        GIT_TAG "376e6d5aeb896c19bdf759b1661fb7f31fb4631c"
        UPDATE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
)

ExternalProject_Get_Property(fixed_point source_dir)
target_include_directories(
        mandelbrot
        PRIVATE
        SYSTEM
        ${source_dir}/include
)

######################################################################
# Boost.SIMD

ExternalProject_Add(
        boost.simd
        URL "https://github.com/johnmcfarlane/boost.simd/archive/5d67a72d77411dc610126151a182ebbae1bfbac6.zip"
        URL_MD5 "97d1aad967965ee79a751d099cd0ee6f"
        UPDATE_COMMAND ""
        INSTALL_COMMAND ""
)

ExternalProject_Get_Property(boost.simd source_dir)
target_include_directories(
        mandelbrot
        PRIVATE
        SYSTEM
        "${source_dir}/include"
)

add_dependencies(
        mandelbrot
        fixed_point
        boost.simd
)
