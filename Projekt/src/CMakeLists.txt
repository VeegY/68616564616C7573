#========================================================================
#           CMake Quell-Verzeichnis             
#========================================================================
# Aufsetzen des Projekts
cmake_minimum_required(VERSION 2.6)
project(stuprotm1516)
#------------------------------------------------------------------------
#           Lade benoetigte Pakete               
#------------------------------------------------------------------------

FIND_PACKAGE(Doxygen)       # Lade Doxygen
FIND_PACKAGE(CUDA REQUIRED) # Lade CUDA
FIND_PACKAGE(MPI REQUIRED)  # Lade MPI

INCLUDE(FindCUDA)
INCLUDE_DIRECTORIES(/usr/local/cuda/include ${MPI_INCLUDE_PATH})

#========================================================================
#               Doxy        
#========================================================================

if(DOXYGEN_FOUND)
configure_file(         # Lade Doxygenkonfiguration
    ${CMAKE_CURRENT_SOURCE_DIR}/../Dokumentation/doxygen_config
    ${CMAKE_CURRENT_BINARY_DIR}/../Dokumentation/doxygen_config
    @ONLY)
add_custom_target(      # Konfiguriere built macro
    doxy ALL
    ${DOXYGEN_EXECUTABLE} 
    ${CMAKE_CURRENT_BINARY_DIR}/../Dokumentation/doxygen_config
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Erzeuge API Dokumentation mit Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)


#========================================================================
#               Executables         
#========================================================================

FILE(GLOB SOURCES "*.cu" "*.cpp" "*.c" "*.h" "*.hpp")
CUDA_ADD_EXECUTABLE(icarus ${SOURCES})


#========================================================================
#               Custom_Target/Unit_Test
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
#               Compiler            
#========================================================================

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math -Wall")
#-I/usr/local/include -pthread -L/usr/local/lib -lmpi")



LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)

TARGET_LINK_LIBRARIES(icarus ${MPI_LIBRARIES})


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




