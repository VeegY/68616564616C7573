/*
* Autor       : David Schneider
* Datum       : 9.4.2016
* Beschreibung: Demo zur Skalierbarkeit der MV-Multiplikation
*               von Bandmatrizen mit NN-Kommunikation
*/

#ifndef __CUDAHELPER_H_
#define __CUDAHELPER_H_

#include <cusparse.h>
#include <cublas.h>
#include <cuda_runtime.h>

template<class T>
void cublas_axpy(cublasHandle_t handle, int n,
    const T *alpha,
    const T *x, int incx,
    T *y, int incy);

template<>
inline void cublas_axpy<float>(cublasHandle_t handle, int n,
    const float *alpha,
    const float *x, int incx,
    float *y, int incy)
{
    if (cublasSaxpy(handle, n, alpha, x, incx, y, incy) != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("cuBLAS axpy failed.");
}

template<>
inline void cublas_axpy<double>(cublasHandle_t handle, int n,
    const double *alpha,
    const double *x, int incx,
    double *y, int incy)
{
    if (cublasDaxpy(handle, n, alpha, x, incx, y, incy) != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("cuBLAS axpy failed.");
}

template<class T>
void cublas_copy(cublasHandle_t handle, int n,
    const T *x, int incx,
    T *y, int incy);

template<>
inline void cublas_copy<float>(cublasHandle_t handle, int n,
    const float *x, int incx,
    float *y, int incy)
{
    if (cublasScopy(handle, n, x, incx, y, incy) != CUBLAS_STATUS_SUCCESS)
            throw std::runtime_error("cuBLAS copy failed.");
}

template<>
inline void cublas_copy<double>(cublasHandle_t handle, int n,
    const double *x, int incx,
    double *y, int incy)
{
    if (cublasDcopy(handle, n, x, incx, y, incy) != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("cuBLAS copy failed.");
}

template<class T>
void cublas_dot(cublasHandle_t handle, int n,
    const T *x, int incx,
    const T *y, int incy,
    T *result);

template<>
inline void cublas_dot<float>(cublasHandle_t handle, int n,
    const float *x, int incx,
    const float *y, int incy,
    float *result)
{
    if (cublasSdot(handle, n, x, incx, y, incy, result) != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("cuBLAS dot failed.");
}

template<>
inline void cublas_dot<double>(cublasHandle_t handle, int n,
    const double *x, int incx,
    const double *y, int incy,
    double *result)
{
    if (cublasDdot(handle, n, x, incx, y, incy, result) != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("cuBLAS dot failed.");
}

template<class T>
void cublas_nrm2(cublasHandle_t handle, int n,
    const T *x, int incx, T *result);

template<>
inline void cublas_nrm2<float>(cublasHandle_t handle, int n,
    const float *x, int incx, float *result)
{
    if (cublasSnrm2(handle, n, x, incx, result) != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("cuBLAS nrm2 failed.");
}

template<>
inline void cublas_nrm2<double>(cublasHandle_t handle, int n,
    const double *x, int incx, double *result)
{
    if (cublasDnrm2(handle, n, x, incx, result) != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("cuBLAS nrm2 failed.");
}

template<class T>
void cublas_scal(cublasHandle_t handle, int n,
    const T *alpha,
    T *x, int incx);

template<>
inline void cublas_scal<float>(cublasHandle_t handle, int n,
    const float *alpha,
    float *x, int incx)
{
    if (cublasSscal(handle, n, alpha, x, incx) != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("cuBLAS scal failed.");
}

template<>
inline void cublas_scal<double>(cublasHandle_t handle, int n,
    const double *alpha,
    double *x, int incx)
{
    if (cublasDscal(handle, n, alpha, x, incx) != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("cuBLAS scal failed.");
}


template <class T>
void cusparse_csrmv(cusparseHandle_t handle, cusparseOperation_t transA,
    int m, int n, int nnz, const T *alpha,
    const cusparseMatDescr_t descrA,
    const T *csrValA,
    const int *csrRowPtrA, const int *csrColIndA,
    const T *x, const T *beta,
    T *y);

template<>
inline void cusparse_csrmv<float>(cusparseHandle_t handle, cusparseOperation_t transA,
    int m, int n, int nnz, const float *alpha,
    const cusparseMatDescr_t descrA,
    const float *csrValA,
    const int *csrRowPtrA, const int *csrColIndA,
    const float *x, const float *beta,
    float *y)
{
    if (cusparseScsrmv(handle, transA, m, n, nnz, alpha,
        descrA, csrValA, csrRowPtrA, csrColIndA, x, beta, y)
        != CUSPARSE_STATUS_SUCCESS)
        throw std::runtime_error("cuSPARSE csrmv failed.");
}

template<>
inline void cusparse_csrmv<double>(cusparseHandle_t handle, cusparseOperation_t transA,
    int m, int n, int nnz, const double *alpha,
    const cusparseMatDescr_t descrA,
    const double *csrValA,
    const int *csrRowPtrA, const int *csrColIndA,
    const double *x, const double *beta,
    double *y)
{
    if (cusparseDcsrmv(handle, transA, m, n, nnz, alpha,
        descrA, csrValA, csrRowPtrA, csrColIndA, x, beta, y)
        != CUSPARSE_STATUS_SUCCESS)
        throw std::runtime_error("cuSPARSE csrmv failed.");
}


#endif // __CUDAHELPER_H_
