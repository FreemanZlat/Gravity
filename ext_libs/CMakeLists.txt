cmake_minimum_required(VERSION 2.8)

include(ExternalProject)

ExternalProject_Add(
    ext_SFML
    GIT_REPOSITORY https://github.com/LaurentGomila/SFML.git
    GIT_TAG 2.3.2
    UPDATE_COMMAND ""
    CMAKE_ARGS -DBUILD_SHARED_LIBS=FALSE
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR_EXTLIBS}/SFML
    BUILD_COMMAND ${BUILD_COMAND_PARAMS}
)

ExternalProject_Add(
    ext_GLM
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG 0.9.7.2
    UPDATE_COMMAND ""
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR_EXTLIBS}/glm
    BUILD_COMMAND ${BUILD_COMAND_PARAMS}
)
