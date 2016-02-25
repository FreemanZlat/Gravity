include(CheckTypeSize)
check_type_size(void* SIZEOF_VOID_PTR)
if("${SIZEOF_VOID_PTR}" STREQUAL "4")
    set(ARCH_32BITS 1)
    set(ARCH_STR "x86")
elseif("${SIZEOF_VOID_PTR}" STREQUAL "8")
    set(ARCH_64BITS 1)
    set(ARCH_STR "x64")
else()
    message(FATAL_ERROR "Unsupported architecture")
    return()
endif()

if(MINGW)
    math(EXPR NUMBER_OF_CPU $ENV{NUMBER_OF_PROCESSORS}*2)
    set(BUILD_COMAND_PARAMS mingw32-make -j ${NUMBER_OF_CPU})
elseif(UNIX)
    EXEC_PROGRAM(nproc OUTPUT_VARIABLE NUM_CPU)
    math(EXPR NUMBER_OF_CPU ${NUM_CPU}*2)
    set(BUILD_COMAND_PARAMS make -j ${NUMBER_OF_CPU})
else()
    set(NUMBER_OF_CPU 1)
endif()

set(GLM_INCLUDE_DIR ${INSTALL_DIR_EXTLIBS}/glm/include)

set(SFML_INCLUDE_DIR ${INSTALL_DIR_EXTLIBS}/SFML/include)
set(SFML_LINK_DIR ${INSTALL_DIR_EXTLIBS}/SFML/lib)

set(SFML_DEBUG "")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(SFML_DEBUG "-d")
endif()

set(SFML_VER "-s${SFML_DEBUG}")

set(SFML_AUDIO "sfml-audio${SFML_VER}")
set(SFML_NETWORK "sfml-network${SFML_VER}")
set(SFML_GRAPHICS "sfml-graphics${SFML_VER}")
set(SFML_WINDOW "sfml-window${SFML_VER}")
set(SFML_SYSTEM "sfml-system${SFML_VER}")

macro(sfml_link name)
    set(libs)
    foreach(arg ${ARGN})
        set(libs ${libs} ${arg})
    endforeach()

    if(UNIX)
        add_definitions(-DSFML_STATIC)
        target_link_libraries(${name}
            ${libs}

            X11
            X11-xcb
            xcb
            xcb-randr
            xcb-image

            GL

            openal
            vorbisfile
            vorbis
            vorbisenc
            ogg
            FLAC

            pthread

            udev
            freetype
            jpeg
        )
    elseif(WIN32)
        add_definitions(-DSFML_STATIC)
        if(MINGW)
            set_target_properties(${name} PROPERTIES LINK_FLAGS "-static -static-libgcc -static-libstdc++")
        endif()
        target_link_libraries(${name}
            sfml-main${SFML_DEBUG}
            ${libs}

            opengl32

            winmm

            ws2_32
            openal32
            vorbisfile
            vorbis
            vorbisenc
            ogg
            FLAC

            freetype
            jpeg
        )
    endif()
endmacro()

