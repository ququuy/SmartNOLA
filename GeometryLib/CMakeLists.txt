# CMakeList.txt : CMake project for GeometryLib, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

add_definitions( -DNOMINMAX -D_USE_MATH_DEFINES  )


include(D:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake)
# CGAL
find_package(CGAL CONFIG REQUIRED)

# CGAL_OSQP
set(OSQP_DIR "D:/dev/osqp/build/install/lib/cmake/osqp")
find_package(OSQP CONFIG REQUIRED)
include(CGAL_OSQP_support)

# CGAL_pointmatcher
set(libpointmatcher_DIR "C:/Program Files (x86)/libpointmatcher/share/libpointmatcher/cmake")
find_package(libpointmatcher REQUIRED)
include(CGAL_pointmatcher_support)

# GLM
find_package(glm CONFIG REQUIRED)

# gl
set(DEV_DIR "D://dev")
set(GL_DIR ${DEV_DIR}/opengl/opengl)
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glad REQUIRED)

# Add source to this project's executable.
#add_executable (GeometryLib "GeometryLib.cpp" "GeometryLib.h" "main.cpp" "IO.h" "IO.cpp" "vis/RenderObject.h" "vis/RenderObject.cpp" "vis/window_util.h")
add_executable (GeometryLib "GeometryLib.cpp" "GeometryLib.h" "IO.h" "IO.cpp" "interact.cpp" "vis/RenderObject.h" "vis/RenderObject.cpp" "vis/window_util.h" "vis/camera.h" "vis/PlaneNode.h" "vis/PlaneNode.cpp"  "vis/status_util.h" "vis/status_util.cpp" "vis/camera.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET GeometryLib PROPERTY CXX_STANDARD 20)
endif()


target_compile_options(
    #${target_name} 
    GeometryLib
    PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/MP>
    #PRIVATE $<$<CONFIG:Debug>:/MDd>
    #PRIVATE $<$<CONFIG:Release>:/MD>
    )

target_include_directories(GeometryLib PRIVATE ${OSQP_INCLUDE_DIR})
target_include_directories(GeometryLib PRIVATE ${libpointmatcher_INCLUDE_DIR})
target_link_directories(GeometryLib PRIVATE ${libpointmatcher_LIBRARIES})
target_link_libraries(GeometryLib PRIVATE CGAL::CGAL)
target_link_libraries(GeometryLib PRIVATE CGAL::OSQP_support)
target_link_libraries(GeometryLib PRIVATE CGAL::pointmatcher_support)
target_link_libraries(GeometryLib PRIVATE glm::glm)

target_include_directories(GeometryLib PRIVATE ${GL_DIR}/include)
target_link_libraries(GeometryLib PRIVATE OpenGL::GL)
target_link_libraries(GeometryLib PRIVATE glad::glad)
target_link_libraries(GeometryLib PRIVATE glfw)

