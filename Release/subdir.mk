################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../main.cpp 

OBJS += \
./main.o 

CPP_DEPS += \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++11 -I/home1/03822/malchera/vtkinclude -I/usr/include/vtk-5.8 -I/usr/include/mpich -I../include -I/home1/03822/malchera/boost/boost_1_55_0/libs/include -O3 -funroll-loops -ftree-vectorize -pedantic -Wall -c -fmessage-length=0 -Wno-deprecated -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<" -DLEGACY_WRITER
	@echo 'Finished building: $<'
	@echo ' '


