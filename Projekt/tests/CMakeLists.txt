CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(stuprotm1516)
#========================================================================
#               Executables         
#========================================================================

FILE(GLOB UNIT_TESTS "*_test.cpp")
FILE(GLOB UNIT_MPITESTS "*_mpitest.cpp")

# alle Testfunktionen werden ctest "bereitgestellt"
FOREACH(testpath ${UNIT_TESTS})
    string(REGEX REPLACE "([^ ]*)/tests/" "" filename ${testpath})
    string(REPLACE ".cpp" "" testname ${filename})
    CUDA_ADD_EXECUTABLE(${testname} ${filename})
    TARGET_LINK_LIBRARIES(${testname} SRCLIB linalg ${MPI_LIBRARIES} ${CUDA_LIBRARIES})
    add_test(${testname} ${testname})
ENDFOREACH(testpath)

# TODO: Ist es sinnvoll, jeden Test pauschal auf 4 Knoten (ohne dezidierte Hostfile) zu starten?
FOREACH(testpath ${UNIT_MPITESTS})
    string(REGEX REPLACE "([^ ]*)/tests/" "" filename ${testpath})
    string(REPLACE ".cpp" "" testname ${filename})
    CUDA_ADD_EXECUTABLE(${testname} ${filename})
    TARGET_LINK_LIBRARIES(${testname} SRCLIB linalg ${MPI_LIBRARIES} ${CUDA_LIBRARIES})
    add_test(${testname} mpirun -np 4 ${testname})
ENDFOREACH(testpath)

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math -Wall")

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)
