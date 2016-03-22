CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(stuprotm1516)
#========================================================================
#               Executables         
#========================================================================

FILE(GLOB UNIT_TESTS "*_test.cpp")
FILE(GLOB UNIT_MPITESTS "*_mpitest.cpp")


#========================================================================
#               Custom_Target/Unit_Test
#========================================================================
CUDA_ADD_LIBRARY(cuda "../src/benchmark_ax_cuda.cu")

# alle Testfunktionen werden ctest "bereitgestellt"
FOREACH(testpath ${UNIT_TESTS})
    # nicht unbedingt optimal, aber es funktioniert fuers erste
    string(REGEX REPLACE "([^ ]*)/tests/" "" filename ${testpath})
    string(REPLACE ".cpp" "" testname ${filename})
    #LIST(APPEND filename "../src/benchmark_ax_cuda.cu") got replaced by cuda_add_library, see above
    CUDA_ADD_EXECUTABLE(${testname} ${filename})
    TARGET_LINK_LIBRARIES(${testname} SRCLIB ${MPI_LIBRARIES})
    TARGET_LINK_LIBRARIES(${testname} cuda ${CUDA_LIBARIES})
    add_test(${testname} ${testname})
ENDFOREACH(testpath)

FOREACH(testpath ${UNIT_MPITESTS})
    # nicht unbedingt optimal, aber es funktioniert fuers erste
    string(REGEX REPLACE "([^ ]*)/tests/" "" filename ${testpath})
    string(REPLACE ".cpp" "" testname ${filename})
    #LIST(APPEND filename "../src/benchmark_ax_cuda.cu") got replaced by cuda_add_library
    CUDA_ADD_EXECUTABLE(${testname} ${filename})
    TARGET_LINK_LIBRARIES(${testname} SRCLIB ${MPI_LIBRARIES})
    TARGET_LINK_LIBRARIES(${testname} cuda ${CUDA_LIBRARIES})
    add_test(${testname} mpirun -np 4 ${testname})
ENDFOREACH(testpath)




#========================================================================
#               Compiler            
#========================================================================

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math -Wall")
#-I/usr/local/include -pthread -L/usr/local/lib -lmpi")

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)


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
