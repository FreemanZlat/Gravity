cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME Gravity)
project(${PROJECT_NAME})

include_directories(${GLM_INCLUDE_DIR})
include_directories(${SFML_INCLUDE_DIR})
link_directories(${SFML_LINK_DIR})

set(CMAKE_CXX_FLAGS "-std=c++11")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG")

set(SOURCES
    src/Gravity.cpp
    src/GravityTree.cpp
    src/main.cpp
)

set(LIBS
    ${SFGUI}
    ${SFML_AUDIO}
    ${SFML_NETWORK}
    ${SFML_GRAPHICS}
    ${SFML_WINDOW}
    ${SFML_SYSTEM}
)

add_executable(${PROJECT_NAME} ${SOURCES})
sfml_link(${PROJECT_NAME} ${LIBS})

add_dependencies(${PROJECT_NAME} ext_SFML)
add_dependencies(${PROJECT_NAME} ext_GLM)

install(TARGETS ${PROJECT_NAME} DESTINATION ${INSTALL_DIR}/bin)
