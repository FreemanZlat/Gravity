cmake_minimum_required(VERSION 2.8)

set(ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

set(INSTALL_DIR ${ROOT_DIR}/release)

set(INSTALL_DIR_EXTLIBS ${ROOT_DIR}/ext_libs_bin)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(INSTALL_DIR_EXTLIBS "${INSTALL_DIR_EXTLIBS}_debug")
endif()

include(${ROOT_DIR}/cmake/macro.cmake)

add_subdirectory(ext_libs)

add_subdirectory(app)
