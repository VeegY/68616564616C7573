#========================================================================
# 			CMake Quell-Verzeichnis				
#========================================================================
# Aufsetzen des Projekts
cmake_minimum_required(VERSION 2.6)
project(stuprotm1516)


#========================================================================
# 				Executables			
#========================================================================

FILE(GLOB SOURCES "*.cu" "*.cpp" "*.c" "*.h")
CUDA_ADD_EXECUTABLE(benchmark_ax ${SOURCES})


#========================================================================
# 				Custom_Target/Unit_Test
#========================================================================

#enables make code
add_custom_target(
     code ALL
     DEPENDS benchmark_ax
     )


#========================================================================
# 				Compiler			
#========================================================================

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++0x -O3 -ffast-math -Wall")
#-I/usr/local/include -pthread -L/usr/local/lib -lmpi")

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)

TARGET_LINK_LIBRARIES(icarus ${MPI_LIBRARIES})


#========================================================================
# 				MODE		
#========================================================================

if (DEBUG_MODE)
LIST(APPEND CMAKE_NVCC_FLAGS "") #TODO insert nvcc debug flags
LIST(APPEND CMAKE_CXX_FLAGS "") #TODO insert gcc debug flags
endif (DEBUG_MODE)

if (RELEASE_MODE)
LIST(APPEND CMAKE_NVCC_FLAGS "") #TODO insert nvcc release flags
LIST(APPEND CMAKE_CXX_FLAGS "-DNODEBUG") #TODO insert gcc release flags
endif (RELEASE_MODE) 


