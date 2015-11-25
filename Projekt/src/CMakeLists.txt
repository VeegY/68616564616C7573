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
FIND_PACKAGE(MPI REQUIRED)    # Lade MPI

INCLUDE(FindCUDA)
INCLUDE_DIRECTORIES(/usr/local/cuda/include ${MPI_INCLUDE_PATH})

#========================================================================
# 				Doxy		
#========================================================================

if(DOXYGEN_FOUND)
configure_file(			# Lade Doxygenkonfiguration
	${CMAKE_CURRENT_SOURCE_DIR}/../Dokumentation/doxygen_config
	${CMAKE_CURRENT_BINARY_DIR}/../Dokumentation/doxygen_config
	@ONLY)
add_custom_target(		# Konfiguriere built macro
	doxy ALL
	${DOXYGEN_EXECUTABLE} 
	${CMAKE_CURRENT_BINARY_DIR}/../Dokumentation/doxygen_config
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Erzeuge API Dokumentation mit Doxygen" VERBATIM
	)
endif(DOXYGEN_FOUND)


#========================================================================
# 				Executables			
#========================================================================

FILE(GLOB SOURCES "*.cu" "*.cpp" "*.c" "*.h")
CUDA_ADD_EXECUTABLE(icarus ${SOURCES})


#========================================================================
# 				Custom_Target/Unit_Test
#========================================================================

#enables make code
add_custom_target(
     code ALL
     DEPENDS icarus
     )


# ctest anschalten
enable_testing()

# Beispiel: create_test(logger_test)
macro(create_test target)
	add_test(
	  NAME ${target} 
	  COMMAND $<TARGET_FILE:${target}>)
endmacro(create_test)



#========================================================================
# 				Compiler			
#========================================================================

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++0x -O3 -ffast-math -Wall ")

#-L/sfw/openmpi/gcc4.9.x/1.8.5/lib/ -I/sfw/openmpi/gcc4.9.x/1.8.5/include -L/usr/local/cuda-7.5/targets/x86_64-linux/lib/ -I/usr/local/cuda-7.5/targets/x86_64-linux/include/ -lmpi -lcudart

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)

#TARGET_LINK_LIBRARIES(icarus /usr/local/cuda/lib64/libcudart.so ${MPI_LIBRARIES})


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


