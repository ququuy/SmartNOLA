# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project ("BuildingPoints")

message(aaaa ${PROJECT_SOURCE_DIR})
set(ROOT_PATH ${PROJECT_SOURCE_DIR})
configure_file(
	"${PROJECT_SOURCE_DIR}/GeometryLib/cmake_definition.h.in"
	"${PROJECT_BINARY_DIR}/GeometryLib/cmake_definition.h"
)

include_directories(${PROJECT_BINARY_DIR}/GeometryLib)

# Include sub-projects.
add_subdirectory ("GeometryLib")
add_subdirectory ("imgui")
