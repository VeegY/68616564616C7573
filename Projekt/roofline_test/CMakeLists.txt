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

ADD_SUBDIRECTORY(axpy)
ADD_SUBDIRECTORY(copy)
ADD_SUBDIRECTORY(dotproduct)
ADD_SUBDIRECTORY(l2norm)
ADD_SUBDIRECTORY(l2norm_test)
ADD_SUBDIRECTORY(matvec)
ADD_SUBDIRECTORY(maxnorm)
ADD_SUBDIRECTORY(vmvr)


