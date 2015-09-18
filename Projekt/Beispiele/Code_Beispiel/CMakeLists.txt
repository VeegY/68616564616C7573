#========================================================================
# 			CMake Quell-Verzeichnis				
#========================================================================
# Aufsetzen des Projekts
cmake_minimum_required(VERSION 2.8)
project(Doxygen_Beispiel)
#------------------------------------------------------------------------
# 			Lade benötigte Pakete				
#------------------------------------------------------------------------
find_package(Doxygen)		# Lade Doxygen
find_package(CUDA REQUIRED)	# Lade CUDA
#find_package(MPI REQUIRED)	# Lade MPI

#========================================================================
# 				 MACROS			
#========================================================================
#------------------------------------------------------------------------
# 		  Macro fuer die Doxygen Dokumentation			
#------------------------------------------------------------------------
# Beispiel: make doxy
# Bemerke: CMAKE_CURRENT_SOURCE_DIR = arg wo ([dir]/cmake [arg]) hier: ../
# Bemerke: CMAKE_CURRENT_BINARY_DIR = dir wo ([dir]/cmake [arg]) hier: ./ bzw pwd
if(DOXYGEN_FOUND)
configure_file(			# Lade Doxygenkonfiguration
	${CMAKE_CURRENT_SOURCE_DIR}/Dokumentation/doxygen_config
	${CMAKE_CURRENT_BINARY_DIR}/../Dokumentation/doxygen_config
	@ONLY)
add_custom_target(		# Konfiguriere built macro
	doxy
	${DOXYGEN_EXECUTABLE} 
	${CMAKE_CURRENT_BINARY_DIR}/../Dokumentation/doxygen_config
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Erzeuge API Dokumentation mit Doxygen" VERBATIM
	)
endif(DOXYGEN_FOUND)

#------------------------------------------------------------------------
# 		            Macro fuer CTest			
#------------------------------------------------------------------------

# ctest anschalten
enable_testing()

macro(create_test target)
	add_test(
	  NAME ${target} 
	  COMMAND $<TARGET_FILE:${target}>)
endmacro(create_test)

#========================================================================
# 				Compiler			
#========================================================================
#------------------------------------------------------------------------
# 		  	     CUDA Compiler			
#------------------------------------------------------------------------
set(CUDA_NVCC_FLAGS "-g -lineinfo -O3 --restrict --use_fast_math --compiler-options -Wall; -g")

#if (DEBUG_MODE) #WIRHABENHIEREXTRADEBUGFLAGS
#  set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-v -Xptxas=-v -G --compiler-options=-rdynamic")
#endif(DEBUG_MODE)

# CMake 2.8.10
#if (CUDA_HOST_COMPILER) #WASSINDDASJETZTFUERFLAGS?!
#  set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -ccbin " ${FEAST_CUDA_HOST_COMPILER})
#endif (CUDA_HOST_COMPILER)

# Host Compiler Flags NICHT uebernehmen
set (CUDA_PROPAGATE_HOST_FLAGS OFF)

#------------------------------------------------------------------------
# 		  	      GCC Compiler			
#------------------------------------------------------------------------
#Standard Flags:
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fopenmp") 

#Standard Flags mit CMAKE_BUILD_TYPE=Debug: -g
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -DCONSTRUCTOR_TEST")

#Standard Flags mit CMAKE_BUILD_TYPE=Debug: -O3 -DNDEBUG
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse2 -mfpmath=sse")


#------------------------------------------------------------------------
# 	     - - - WIP - - -  MPI Compiler  - - - WIP - - -
#------------------------------------------------------------------------
#include_directories(${MPI_INCLUDE_PATH})
#set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

#========================================================================
# 			     Verzeichnis-Linking			
#========================================================================
#------------------------------------------------------------------------
# 		  	    Verzeichnisarchitektur strukturieren			
#------------------------------------------------------------------------
include_directories(			# Verzeichnisse linken
	${PROJECT_SOURCE_DIR} Code
	)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY 	# Executables
	${CMAKE_BINARY_DIR}/EXE
	)
#------------------------------------------------------------------------
# 		  	    CMakeLists Linking			
#------------------------------------------------------------------------
add_subdirectory(Code)
