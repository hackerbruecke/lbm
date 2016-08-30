# Program name
NAME=lbm
# Files
CPPFILES=main.cpp
# Tools
CC=g++
LD=g++
# We will use gcc (instead of ld) to link our object files together because it automatically discovers appropriate linking flags
# Flags/Options
CFLAGS=-std=c++11 -O3 -funroll-loops -ftree-vectorize
LDFLAGS=-fopenmp -Wl,-R -Wl,build
LIBS=-lpugixml -lboost_mpi -lboost_program_options -lboost_filesystem -lboost_system -lvtkIO -lvtkFiltering -lvtkCommon -lvtkCharts  -lvtkCommon  -lvtkDICOMParser  -lvtkexoIIc  -lvtkexpat  -lvtkFiltering  -lvtkfreetype  -lvtkftgl  -lvtkGenericFiltering  -lvtkGeovis  -lvtkGraphics  -lvtkhdf5  -lvtkhdf5_hl  -lvtkHybrid  -lvtkImaging  -lvtkInfovis  -lvtkIO  -lvtkjpeg  -lvtkmetaio  -lvtkNetCDF  -lvtkpng  -lvtkproj4  -lvtkRendering  -lvtksqlite  -lvtksys  -lvtktiff  -lvtkverdict  -lvtkViews  -lvtkVolumeRendering  -lvtkWidgets -lvtkzlib

BOOST_INC_LOCATION=$(TACC_BOOST_MPI_INCLUDE) /home1/03822/malchera/boost/boost_1_55_0/libs/include
BOOST_LIB_LOCATION=$(TACC_BOOST_MPI_LIB) /home1/03822/malchera/boost/boost_1_55_0/libs/lib
VTK_INC_LOCATION=./include/vtk-5.10
VTK_LIB_LOCATION=./lib/vtk-5.10

#######################################################

INCLUDE_DIR=./include $(BOOST_INC_LOCATION) $(VTK_INC_LOCATION)
LIBRARY_DIR=./lib $(VTK_LIB_LOCATION) $(BOOST_LIB_LOCATION) 

# Tasks
all: clean cc ld
	@echo Done!

cc:
	@$(foreach file, $(CPPFILES),\
		echo Compiling $(file)...;\
		$(CC) $(CFLAGS) $(foreach d, $(INCLUDE_DIR), -I$d) -c $(CPPFILES)\
		;)

ld:
	@echo Linking...
	$(LD) $(foreach d, $(LIBRARY_DIR), -L$d) $(LDFLAGS) -o $(NAME) *.o $(LIBS)

clean:
	@echo Deleting all object files...
	@rm -f *.o $(NAME)

