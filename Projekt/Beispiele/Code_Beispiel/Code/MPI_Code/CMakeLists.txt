# CMake Code-Verzeichnis

# Aufsetzen
cmake_minimum_required(VERSION 2.6)

# Preprocessor-Variable DIR_VECTOR
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDIR_MPI=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
# debug-Ausgabe
message(${CMAKE_CXX_FLAGS})

# Erzeugen der Executables
add_executable(MPI_TEST mpi_test.cpp)

# Linken der Bibliotheken

target_link_libraries(MPI_TEST ${MPI_LIBRARIES})
