# CMakeLists.txt for hubomz
cmake_minimum_required(VERSION 2.6)

project(ccdwrapper)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

if(APPLE)
  include_directories(/opt/local/include) # MacPorts
  link_directories(/opt/local/lib)
  find_library(OPENGL_LIBRARY OpenGL)
else()
  find_library(OPENGL_LIBRARY GL)
  find_library(GLU_LIBRARY GLU)
  set(OPENGL_LIBRARY ${OPENGL_LIBRARY} ${GLU_LIBRARY})
endif()

find_library(GLUT_LIBRARY glut)

include(FindPkgConfig)

pkg_search_module(EIGEN3 REQUIRED eigen3>=3)
include_directories(${EIGEN3_INCLUDE_DIRS})

pkg_search_module(CCD REQUIRED libccd ccd)

set(CMAKE_C_FLAGS "-Wall -g")
set(CMAKE_CXX_FLAGS "-Wall -g")

add_executable(ccddemo ccddemo.cpp CCDWrapper.cpp GlCamera.cpp)
target_link_libraries(ccddemo ${GLUT_LIBRARY} ${CCD_LIBRARIES} ${OPENGL_LIBRARY})