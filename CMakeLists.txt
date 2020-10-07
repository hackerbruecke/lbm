cmake_minimum_required(VERSION 3.10)
project(lbm)

set(CMAKE_CXX_STANDARD 17)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(KEEP_RPATHS NO_OUTPUT_DIRS)

include_directories(include)
include_directories(include/io)

#include(${VTK_USE_FILE})
set(HEADERS
  include/io/configuration.h
  include/io/scenario.h
  include/io/vtk.h
  include/io/vtk.hpp
  include/io/vtk_legacy.hpp
  include/boundary.h
  include/boundary.hpp
  include/cell.h
  include/cell.hpp
  include/collision.h
  include/collision.hpp
  include/domain.h
  include/domain.hpp
  include/helper.h
  include/lbmdefinitions.h
  include/model.h
  include/parallel.h
)

set(SOURCES
  src/main.cpp
)

add_executable(lbm ${HEADERS} ${SOURCES})

target_compile_options(lbm PUBLIC -fopenmp)
target_link_libraries(lbm PRIVATE ${CONAN_LIBS})
target_link_libraries(lbm PRIVATE -fopenmp ${CMAKE_DL_LIBS})
