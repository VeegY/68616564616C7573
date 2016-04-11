cmake_minimum_required(VERSION 2.6)
project(bandscal)

FIND_PACKAGE(CUDA REQUIRED)
FIND_PACKAGE(MPI REQUIRED)

INCLUDE(FindCUDA)
INCLUDE_DIRECTORIES(/usr/local/cuda/include ${MPI_INCLUDE_PATH})

CUDA_ADD_EXECUTABLE(slvb "slvb.cpp")
TARGET_LINK_LIBRARIES(slvb ${MPI_LIBRARIES} ${CUDA_cusparse_LIBRARY})
CUDA_ADD_CUBLAS_TO_TARGET(slvb)

CUDA_ADD_EXECUTABLE(mvb "mvb.cpp")
TARGET_LINK_LIBRARIES(mvb ${MPI_LIBRARIES} ${CUDA_cusparse_LIBRARY})
CUDA_ADD_CUBLAS_TO_TARGET(mvb)

LIST(APPEND CMAKE_CXX_FLAGS "-fopenmp -std=c++11 -O3 -ffast-math -Wall")

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)
