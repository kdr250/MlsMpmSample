cmake_minimum_required(VERSION 3.19)

if (APPLE AND DEFINED ENV{RPATH})
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_INSTALL_RPATH $ENV{RPATH})
endif()

if (DEFINED CMAKE_TOOLCHAIN_FILE)
    set(VCPKG_CHAINLOAD_TOOLCHAIN_FILE ${CMAKE_TOOLCHAIN_FILE})
endif()

if (DEFINED ENV{VCPKG_ROOT})
    cmake_path(SET VCPKG_PATH NORMALIZE $ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake)
    set(CMAKE_TOOLCHAIN_FILE ${VCPKG_PATH})
else()
    message(FATAL_ERROR "Please install vcpkg and set 'VCPKG_ROOT' environmental variable")    
endif()

project(main VERSION 1.0)

if (EMSCRIPTEN)
    set(VCPKG_TARGET_TRIPLET "wasm32-emscripten")
else()
    find_package(SDL2 CONFIG REQUIRED)
    find_package(sdl2-gfx CONFIG REQUIRED)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/Zc:__cplusplus>")

file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(
    main
    ${SOURCES}
)

if (EMSCRIPTEN)

    set(USE_FLAGS "-s USE_SDL=2 -s USE_SDL_GFX=2 -s ALLOW_MEMORY_GROWTH --preload-file resources/")
    set(CMAKE_CXX_FLAGS ${USE_FLAGS})

    set_target_properties(main PROPERTIES SUFFIX ".html")

    file(COPY ${PROJECT_SOURCE_DIR}/resources DESTINATION ${CMAKE_BINARY_DIR})

else()

    target_link_libraries(
        main PRIVATE
        $<TARGET_NAME_IF_EXISTS:SDL2::SDL2main>
        $<IF:$<TARGET_EXISTS:SDL2::SDL2>,SDL2::SDL2,SDL2::SDL2-static>
        SDL2::SDL2_gfx
    )

    add_custom_command(
        TARGET main POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/resources $<TARGET_FILE_DIR:main>/resources
    )

endif()
