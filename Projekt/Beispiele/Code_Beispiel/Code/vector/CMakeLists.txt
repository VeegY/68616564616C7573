# CMake Code-Verzeichnis

# Aufsetzen
cmake_minimum_required(VERSION 2.6)

# Preprocessor-Variable DIR_VECTOR
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDIR_VECTOR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
# debug-Ausgabe
message(${CMAKE_CXX_FLAGS})

# Erzeugen der Executables
add_executable(VECTOR_CPU_Test vector_CPU_TEST.cpp)

# Aufsetzen der Tests
add_test(VECTOR_CPU_Test VECTOR_CPU_Test)
