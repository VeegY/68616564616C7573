cmake_minimum_required(VERSION 2.6)
project(stuprotm1516)

FIND_PACKAGE(CUDA REQUIRED) # Lade CUDA
FIND_PACKAGE(MPI REQUIRED)  # Lade MPI

INCLUDE(FindCUDA)
INCLUDE_DIRECTORIES(/usr/local/cuda/include ${MPI_INCLUDE_PATH})

#========================================================================
#               Executables         
#========================================================================

FILE(GLOB SOURCES "*.cu" "*.cpp")# "*.c" "*.h" "*.hpp") # adding all source files
SET(SRCFILES ${SOURCES})        # there must be a
ADD_LIBRARY(SRCLIB ${SRCFILES}) # more compact way

CUDA_ADD_EXECUTABLE(icarus ${SOURCES})
TARGET_LINK_LIBRARIES(icarus SRCLIB ${MPI_LIBRARIES})


#========================================================================
#               Custom_Target/Unit_Test
#========================================================================

#enables make code
add_custom_target(
     code ALL
     DEPENDS icarus
     )


#========================================================================
#               Compiler            
#========================================================================

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math -Wall")
#-I/usr/local/include -pthread -L/usr/local/lib -lmpi")


LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)


#========================================================================
#               MODE        
#========================================================================

if (DEBUG_MODE)
LIST(APPEND CMAKE_NVCC_FLAGS "") #TODO insert nvcc debug flags
LIST(APPEND CMAKE_CXX_FLAGS "") #TODO insert gcc debug flags
endif (DEBUG_MODE)

if (RELEASE_MODE)
LIST(APPEND CMAKE_NVCC_FLAGS "") #TODO insert nvcc release flags
LIST(APPEND CMAKE_CXX_FLAGS "-DNODEBUG") #TODO insert gcc release flags
endif (RELEASE_MODE) 
