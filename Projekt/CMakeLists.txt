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
    ${CMAKE_CURRENT_SOURCE_DIR}/documentation/doxygen_config
    ${CMAKE_CURRENT_BINARY_DIR}/documentation/doxygen_config
    @ONLY)
add_custom_target(      # Konfiguriere built macro
    doxy ALL
    ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/documentation/doxygen_config
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Erzeuge API Dokumentation mit Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

INCLUDE(CTest) # beinhaltet enable_testing()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests/out)
