# CMake Code-Verzeichnis

# Aufsetzen
cmake_minimum_required(VERSION 2.6)

# Preprocessor-Variable DIR_CODE
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDIR_CODE=\"${CMAKE_CURRENT_SOURCE_DIR}\"")

# Include-Verzeichnis
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR} vector GPU_Code
	)

# CMakeLists.txt der Unterverzeichnisse referenzieren
add_subdirectory(GPU_Code)
#add_subdirectory(MPI_Code)
add_subdirectory(vector)

