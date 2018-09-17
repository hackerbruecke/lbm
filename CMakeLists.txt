cmake_minimum_required(VERSION 3.12)
project(lbm)


set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /usr/lib64/cmake/vtk)

find_package(VTK REQUIRED)
find_package(pugixml REQUIRED)

include_directories(include)
include_directories(include/io)

include(${VTK_USE_FILE})

add_executable(lbm
        include/io/configuration.h
        include/io/scenario.h
        include/io/vtk.h
        include/io/vtk.hpp
        include/io/vtk_legacy.hpp
        #include/pugixml/pugiconfig.hpp
        #include/pugixml/pugixml.hpp
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
        src/main.cpp)


target_compile_options(lbm PUBLIC -fopenmp)
target_link_libraries(lbm pugixml boost_program_options boost_filesystem boost_system ${VTK_LIBRARIES})
target_link_libraries(lbm -fopenmp)