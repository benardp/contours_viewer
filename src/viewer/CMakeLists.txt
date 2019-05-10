include(${PROJECT_SOURCE_DIR}/cmake/opengl.cmake)

include_directories(.)
aux_source_directories(SOURCES "" .)
add_executable(contours_viewer ${SOURCES})
target_link_libraries(contours_viewer geogram_gfx geogram ${GLFW_LIBRARIES})
install_runtime_targets(contours_viewer)

set_target_properties(contours_viewer PROPERTIES FOLDER "GEOGRAM")
