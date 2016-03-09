cmake_minimum_required(VERSION 2.6)
project(stuprotm1516)

CUDA_ADD_EXECUTABLE(solverbench "solverbench.cpp")
TARGET_LINK_LIBRARIES(solverbench SRCLIB ${MPI_LIBRARIES})

target_include_directories(solverbench 
     PUBLIC "include" 
     PUBLIC "../src/include")

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math -Wall")

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)
