#========================================================================
# 			CMake Quell-Verzeichnis				
#========================================================================
# Aufsetzen des Projekts
cmake_minimum_required(VERSION 2.6)
project(stuprotm1516)
#------------------------------------------------------------------------
# 			Lade benötigte Pakete				
#------------------------------------------------------------------------

FIND_PACKAGE(Doxygen)		# Lade Doxygen
FIND_PACKAGE(CUDA REQUIRED)	# Lade CUDA
FIND_PACKAGE(MPI REQUIRED)    	# Lade MPI

INCLUDE(FindCUDA)
INCLUDE_DIRECTORIES(/usr/local/cuda/include ${MPI_INCLUDE_PATH})


#------------------------------------------------------------------------
# 			Setze globale Flags			
#------------------------------------------------------------------------

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++0x -O3 -ffast-math -Wall")

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -Xptxas -dlcm=cg -O3 -use_fast_math)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
#LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_32,code=sm_32)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)


ADD_SUBDIRECTORY(axpy)
ADD_SUBDIRECTORY(copy)
ADD_SUBDIRECTORY(dotproduct)
ADD_SUBDIRECTORY(l2norm)
ADD_SUBDIRECTORY(matvec)
ADD_SUBDIRECTORY(maxnorm)


