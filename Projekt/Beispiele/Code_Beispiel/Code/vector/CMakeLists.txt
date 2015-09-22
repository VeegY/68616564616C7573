# CMake Code-Verzeichnis

# Aufsetzen
cmake_minimum_required(VERSION 2.8)

# Preprocessor-Variable DIR_VECTOR
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDIR_VECTOR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")

# Erzeugen der Executables
add_executable(VECTOR_CPU_TEST vector_CPU_TEST.cpp)

# Aufsetzen der Tests
create_test(VECTOR_CPU_TEST)
