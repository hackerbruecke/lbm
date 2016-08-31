# Program name
NAME=lbm
BUILDDIR=build
SOURCEDIR=src
# Files
OBJFILES=main.o
# Tools
CC=g++
LD=g++
# We will use gcc (instead of ld) to link our object files together because it automatically discovers appropriate linking flags
# Flags/Options
CFLAGS=-std=c++11 -O3 -funroll-loops -ftree-vectorize -Wno-deprecated
LDFLAGS=-fopenmp #-Wl,-R -Wl,`pwd`
LIBS=-lpugixml -lboost_program_options -lboost_filesystem -lboost_system -lvtkIO -lvtkFiltering -lvtkCommon -lvtksys -ldl -lvtkzlib

BOOST_INC_LOCATION=$(TACC_BOOST_INC)
BOOST_LIB_LOCATION=$(TACC_BOOST_LIB)
VTK_INC_LOCATION=./include/vtk-5.10
VTK_LIB_LOCATION=./lib/vtk-5.10 

#######################################################

INCLUDE_DIR=./include $(BOOST_INC_LOCATION) $(VTK_INC_LOCATION)
LIBRARY_DIR=./lib $(VTK_LIB_LOCATION) $(BOOST_LIB_LOCATION) 

# Tasks
all: cc ld
	@echo Done!

$(BUILDDIR)/%.o : $(SOURCEDIR)/%.cpp
	$(CC) $(CFLAGS) $(foreach d, $(INCLUDE_DIR), -I$d) -c $< -o $@

cc: $(foreach file, $(OBJFILES), $(BUILDDIR)/$(file))

ld:
	@echo Linking...
	$(LD) $(foreach d, $(LIBRARY_DIR), -L$d) $(LDFLAGS) -o $(BUILDDIR)/$(NAME) $(foreach obj, $(OBJFILES), $(BUILDDIR)/$(obj)) $(LIBS)

clean:
	@echo Deleting all object files...
	@cd $(BUILDDIR) && \
	rm -f $(OBJFILES) $(NAME) && \
	cd -

