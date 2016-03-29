cmake_minimum_required(VERSION 2.6)
project(stuprotm1516)

#========================================================================
#               Executables         
#========================================================================

# Globale cxx-Library
FILE(GLOB SOURCES "*.cu" "*.cpp")# "*.c" "*.h" "*.hpp") # adding all source files
ADD_LIBRARY(SRCLIB ${SOURCES}) 
FILE(GLOB DEMO "*_demo.cpp")
IF(NOT DEMO STREQUAL "")
  CUDA_ADD_EXECUTABLE(demo ${DEMO})
  TARGET_LINK_LIBRARIES(demo SRCLIB ${MPI_LIBRARIES})
ENDIF(NOT DEMO STREQUAL "")

#Globale cuda-Kernels
CUDA_ADD_LIBRARY(linalg "linalg.cu")

#========================================================================
#               Compiler            
#========================================================================

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math -Wall")

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)

