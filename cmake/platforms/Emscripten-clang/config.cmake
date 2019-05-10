set(VORPALINE_ARCH_64 true)
include(${GEOGRAM_SOURCE_DIR}/cmake/platforms/Emscripten-clang.cmake)
add_flags(CMAKE_CXX_FLAGS -m64 -s FORCE_FILESYSTEM=1)
add_flags(CMAKE_C_FLAGS -m64 -s FORCE_FILESYSTEM=1)

