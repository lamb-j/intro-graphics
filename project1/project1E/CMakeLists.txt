cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project1E)
SET(VTK_DIR /path/to/VTK6.0.0)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(project1E project1E)

if(VTK_LIBRARIES)
  target_link_libraries(project1E ${VTK_LIBRARIES})
else()
  target_link_libraries(project1E vtkHybrid)
endif()

