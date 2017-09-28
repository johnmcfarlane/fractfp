cmake_minimum_required(VERSION 3.5.1)

project(mandelbrot)

include(ExternalProject)

set(CMAKE_CXX_STANDARD 14)

################################################################################
# mandelbrot target

set(SOURCE_FILES main.cpp mandelbrot.h)
add_executable(mandelbrot ${SOURCE_FILES})

set(CXX_FLAGS "-mavx -DSG14_USE_INT128 -Wall -Wextra")

set_target_properties(
        mandelbrot
        PROPERTIES COMPILE_FLAGS ${CXX_FLAGS}
)

################################################################################
# cnl

include(ExternalProject)

ExternalProject_Add(
        cnl
        GIT_REPOSITORY "https://github.com/johnmcfarlane/cnl.git"
        GIT_TAG "754e4c65567d843081fb4cf45552764c35e30130"
        UPDATE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
)

ExternalProject_Get_Property(cnl source_dir)
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
        URL "https://github.com/johnmcfarlane/boost.simd/archive/b1e046e3ffa720226e45512c76c505f34cde7be2.zip"
        URL_MD5 "74fc0bc279a18f931e4d2ec026af8dd4"
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

######################################################################
# Boost

find_package(Boost REQUIRED)

if(Boost_VERSION LESS 106200)
    message(FATAL_ERROR "Boost version 1.62 is required by Boost.SIMD")
endif(Boost_VERSION LESS 106200)

target_link_libraries(mandelbrot Boost::boost)

######################################################################

add_dependencies(
        mandelbrot
        cnl
        boost.simd
        Boost::boost
)
