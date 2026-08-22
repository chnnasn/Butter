foreach(example_name hello_butter falling_boxes chain playground physics2d_playground)
    add_executable(${example_name} ${example_name}.cpp)
    target_link_libraries(${example_name} PRIVATE butter)
endforeach()

# Optional GLFW example. GLFW is treated as an external dependency and is never
# vendored into Butter. Point BUTTER_GLFW_DIR at a GLFW source tree (for example
# E:/Github/glfw) or rely on the Windows default below.
set(BUTTER_GLFW_DIR "" CACHE PATH "Path to a GLFW source directory (e.g. E:/Github/glfw)")
if(NOT BUTTER_GLFW_DIR)
    if(WIN32 AND EXISTS "E:/Github/glfw/CMakeLists.txt")
        set(BUTTER_GLFW_DIR "E:/Github/glfw")
    endif()
endif()

if(BUTTER_GLFW_DIR AND EXISTS "${BUTTER_GLFW_DIR}/CMakeLists.txt")
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

    add_subdirectory("${BUTTER_GLFW_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/glfw" EXCLUDE_FROM_ALL)

    find_package(OpenGL REQUIRED)
    add_executable(exploding_crates exploding_crates.cpp)
    target_link_libraries(exploding_crates PRIVATE butter glfw OpenGL::GL)
    target_compile_definitions(exploding_crates PRIVATE BUTTER_EXPLOSION_USE_GLFW=1)
else()
    message(STATUS "GLFW not found - skipping exploding_crates (set BUTTER_GLFW_DIR)")
endif()
