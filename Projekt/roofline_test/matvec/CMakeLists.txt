cmake_minimum_required(VERSION 2.6)
project(stuprotm1516)


#========================================================================
# 				Executables			
#========================================================================

FILE(GLOB SOURCES "*.cu" "*.cpp" "*.c" "*.h")
CUDA_ADD_EXECUTABLE(matvec ${SOURCES})
TARGET_LINK_LIBRARIES(matvec ${MPI_LIBRARIES})


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


