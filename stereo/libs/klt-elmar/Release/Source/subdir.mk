################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Source/convolve.c \
../Source/error.c \
../Source/klt.c \
../Source/klt_util.c \
../Source/pnmio.c \
../Source/pyramid.c \
../Source/selectGoodFeatures.c \
../Source/stereo_util.c \
../Source/storeFeatures.c \
../Source/trackFeatures.c \
../Source/writeFeatures.c 

OBJS += \
./Source/convolve.o \
./Source/error.o \
./Source/klt.o \
./Source/klt_util.o \
./Source/pnmio.o \
./Source/pyramid.o \
./Source/selectGoodFeatures.o \
./Source/stereo_util.o \
./Source/storeFeatures.o \
./Source/trackFeatures.o \
./Source/writeFeatures.o 

C_DEPS += \
./Source/convolve.d \
./Source/error.d \
./Source/klt.d \
./Source/klt_util.d \
./Source/pnmio.d \
./Source/pyramid.d \
./Source/selectGoodFeatures.d \
./Source/stereo_util.d \
./Source/storeFeatures.d \
./Source/trackFeatures.d \
./Source/writeFeatures.d 


# Each subdirectory must supply rules for building sources it contributes
Source/%.o: ../Source/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -DHAVE_IPP -I$(IPPROOT)/include -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


