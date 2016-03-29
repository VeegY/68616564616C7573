cmake_minimum_required(VERSION 2.6)
project(stuprotm1516)

CUDA_ADD_EXECUTABLE(mv_bench_gpu "mv_bench_gpu.cpp")
ADD_EXECUTABLE(mv_bench "mv_bench.cpp")
TARGET_LINK_LIBRARIES(mv_bench_gpu SRCLIB ${MPI_LIBRARIES} linalg ${CUDA_LIBRARIES})
TARGET_LINK_LIBRARIES(mv_bench SRCLIB ${MPI_LIBRARIES})

target_include_directories(mv_bench 
     PUBLIC "../src/include")
target_include_directories(mv_bench_gpu 
     PUBLIC "../src/include")

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -ffast-math -Wall")

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)
