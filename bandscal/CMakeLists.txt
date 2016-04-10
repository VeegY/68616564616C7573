cmake_minimum_required(VERSION 2.6)
project(bandscal)

FIND_PACKAGE(CUDA REQUIRED)
FIND_PACKAGE(MPI REQUIRED)

INCLUDE(FindCUDA)
INCLUDE_DIRECTORIES(/usr/local/cuda/include ${MPI_INCLUDE_PATH})

CUDA_ADD_EXECUTABLE(bandscal "bandscal.cpp")
TARGET_LINK_LIBRARIES(bandscal ${MPI_LIBRARIES} ${CUDA_cusparse_LIBRARY})
CUDA_ADD_CUBLAS_TO_TARGET(bandscal)

LIST(APPEND CMAKE_CXX_FLAGS "-fopenmp -std=c++11 -O3 -ffast-math -Wall")

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)
