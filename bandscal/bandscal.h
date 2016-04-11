/*
 * Autor       : David Schneider
 * Datum       : 9.4.2016
 * Beschreibung: Schnelle Demo zur Skalierbarkeit der MV-Multiplikation 
 *               von Bandmatrizen mit NN-Kommunikation
 */

#ifndef __BANDSCAL_H_
#define __BANDSCAL_H_

#include <map>
#include <thread>
#include <algorithm>
#include <exception>
#include <cassert>
#include <vector>
#include <iostream>
#include <chrono>
#include <mpi.h>

#include "cudahelper.h"

enum arch_t { ARCH_CPU, ARCH_GPU };

template<class Scalar>
struct ScalarTraits;

template<>
struct ScalarTraits<float>
{
    static constexpr MPI_Datatype MPI_Type = MPI_FLOAT;
};

template<>
struct ScalarTraits<double>
{
    static constexpr MPI_Datatype MPI_Type = MPI_DOUBLE;
};

void distribute_among_nodes(int n, int& l, int& p)
{
    int myrank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (n % nprocs == 0)
        l = p = n / nprocs;
    else
    {
        l = n / nprocs + 1;
        assert((nprocs - 1)*l < n);
        p = n - (nprocs - 1)*l;
    }
    assert(l > 0);
}

template <class Scalar>
class BVector
{
    Scalar *_data;

    int _n, _m, _B, _l, _p;
    int _myrank, _nprocs;
    int _offset;
    bool _iam_first, _iam_last;

    arch_t _arch;
    int _nloc, _length;

    MPI_Group _mygroup;
    MPI_Win _win_prev, _win_next;
    Scalar *_buf_prev, *_buf_next;

    int _local_offset;

    cublasHandle_t _cublas_handle;

public:
    BVector(int n, int B, arch_t arch, cublasHandle_t cublas_handle) :
        _data(nullptr),
        _n(n), _B(B),
        _arch(arch),
        _cublas_handle(cublas_handle)
    {
        assert(_B > 0);
        MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &_myrank);
        distribute_among_nodes(_n, _l, _p);
        assert(_l >= _B);
        if (_myrank == 0)
            _offset = 0;
        else _offset = _myrank*_l - _B;
        _iam_first = _myrank == 0;
        _iam_last = _myrank == _nprocs - 1;
        _nloc = _iam_last ? _p : _l;
        if (_nprocs == 1)
            _length = _l;
        else if (_iam_first)
            _length = _l + _B;
        else if (_iam_last)
            _length = _p + _B;
        else
            _length = _l + 2 * _B;

        _local_offset = (!_iam_first) * _B;

        // Speicher allokieren
        switch (arch)
        {
        case ARCH_CPU:
            _data = new Scalar[_length];
            break;
        case ARCH_GPU:
            gpuMalloc(&_data,sizeof(Scalar)*_length);
            break;
        }
        memset(_data, 0, sizeof(Scalar)*_length);

        // Mpi windows erstellen
        int win_prev_len = _iam_first ? 0 : _B;
        int win_next_len = _iam_last ? 0 : _B;
        _buf_prev = _iam_first ? nullptr : _data;
        _buf_next = _iam_last ? nullptr : (_data + win_prev_len + _nloc);
        MPI_Comm_group(MPI_COMM_WORLD, &_mygroup);
        int mpie1 = MPI_Win_create(
            _buf_prev, win_prev_len*sizeof(Scalar), sizeof(Scalar),
            MPI_INFO_NULL, MPI_COMM_WORLD, &_win_prev);
        int mpie2 = MPI_Win_create(
            _buf_next, win_next_len*sizeof(Scalar), sizeof(Scalar),
            MPI_INFO_NULL, MPI_COMM_WORLD, &_win_next);
        if (mpie1 || mpie2) 
            throw std::runtime_error("Creation of MPI_Windows failed.");
    }

    ~BVector()
    {
        MPI_Win_free(&_win_prev);
        MPI_Win_free(&_win_next);

        switch (_arch)
        {
        case ARCH_CPU:
            if (_data) delete[] _data;
            break;
        case ARCH_GPU:
            if(_data) cudaFree(_data);
            break;
        }
    }

    cublasHandle_t get_cublas_handle() const { return _cublas_handle; }

    void start_comm()
    {
        cudaDeviceSynchronize();
        MPI_Win_post(_mygroup, 0, _win_prev);
        MPI_Win_post(_mygroup, 0, _win_next);
        MPI_Win_start(_mygroup, 0, _win_prev);
        MPI_Win_start(_mygroup, 0, _win_next);
        // zum vorgänger
        if (!_iam_first)
            MPI_Put(_data+_B, _B, ScalarTraits<Scalar>::MPI_Type,
                _myrank - 1, 0, _B, ScalarTraits<Scalar>::MPI_Type, _win_next);
        // zum nachfolger
        if (!_iam_last)
            MPI_Put(_data+_length-2*_B, _B, ScalarTraits<Scalar>::MPI_Type,
                _myrank + 1, 0, _B, ScalarTraits<Scalar>::MPI_Type, _win_prev);
    }

    void finish_comm()
    {
        cudaDeviceSynchronize();
        MPI_Win_complete(_win_prev);
        MPI_Win_complete(_win_next);
        MPI_Win_wait(_win_prev);
        MPI_Win_wait(_win_next);
    }

    Scalar* data()
    {
        return _data;
    }

    const Scalar* data() const
    {
        return _data;
    }

    Scalar* local_data()
    {
        return _data + _local_offset;
    }

    const Scalar* local_data() const
    {
        return _data + _local_offset;
    }

    // Zugriff auf den Datenteil (mit Buffer)
    Scalar& operator[](int idx)
    {
        assert(idx < _length);
        return _data[idx];
    }

    // Zugriff auf den Datenteil (mit Buffer)
    const Scalar& operator[](int idx) const
    {
        assert(idx < _length);
        return _data[idx];
    }

    // Zugriff auf den Datenteil (ohne Buffer)
    Scalar& operator()(int idx)
    {
        assert(idx < _nloc);
        return _data[idx + _local_offset];
    }

    // Zugriff auf den Datenteil (ohne Buffer)
    const Scalar& operator()(int idx) const
    {
        assert(idx < _nloc);
        return _data[idx + _local_offset];
    }

    void fill_with(const Scalar& val)
    {
        for (int i = _local_offset; i < _local_offset + _nloc; i++)
            _data[i] = val;
    }
    
    int get_nloc() const { return _nloc; }

    void copy(BVector& dst) const
    {
        assert(_nloc == dst._nloc);
        switch (_arch)
        {
        case ARCH_CPU:
#           pragma omp parallel for
            for (int i = _local_offset; i < _nloc + _local_offset; i++)
                dst[i] = _data[i];
            break;
        case ARCH_GPU:
            cublas_copy(_cublas_handle, _nloc, _data + _local_offset, 1, dst._data + dst._local_offset, 1);
            break;
        }
    }

    Scalar dot(const BVector& other) const
    {
        assert(_nloc == other.get_nloc());
        Scalar res = 0, resloc = 0;
        switch (_arch)
        {
        case ARCH_CPU:
#           pragma omp parallel for
            for (int i = _local_offset; i < _nloc + _local_offset; i++)
                resloc += _data[i] * other[i];
            break;
            
        case ARCH_GPU:
            cublas_dot(_cublas_handle, _nloc, _data + _local_offset, 1, 
                other._data + other._local_offset, 1, &resloc);
            cudaDeviceSynchronize();
            break;
        }
        MPI_Allreduce(&resloc, &res, 1, ScalarTraits<Scalar>::MPI_Type, MPI_SUM, MPI_COMM_WORLD);
        return res;
    }

    int get_length() { return _length; }

    void scal(const Scalar& alpha)
    {
        switch (_arch)
        {
        case ARCH_CPU:
#           pragma omp parallel for
            for (int i = _local_offset; i < _nloc + _local_offset; i++)
                _data[i] *= alpha;
            break;
        case ARCH_GPU:
            cublas_scal(_cublas_handle, _nloc, &alpha, _data + _local_offset, 1);
            break;
        }
    }

    Scalar nrm2() const
    {
        Scalar res = 0.0, resloc = 0.0;
        switch (_arch)
        {
        case ARCH_CPU:
#           pragma omp parallel for
            for (int i = _local_offset; i < _nloc + _local_offset; i++)
                resloc += _data[i] * _data[i];
            break;
        case ARCH_GPU:
            cublas_nrm2(_cublas_handle, _nloc, _data + _local_offset, 1, &resloc);
            resloc *= resloc;
            cudaDeviceSynchronize();
            break;
        }
        MPI_Allreduce(&resloc, &res, 1, ScalarTraits<Scalar>::MPI_Type, MPI_SUM, MPI_COMM_WORLD);
        return sqrt(res);
    }

    void axpy(const Scalar& alpha, const BVector& x)
    {
        assert(_nloc == x._nloc);
        switch (_arch)
        {
        case ARCH_CPU:
#           pragma omp parallel for
            for (int i = _local_offset; i < _nloc + _local_offset; i++)
                _data[i] += alpha * x[i];
            break;
        case ARCH_GPU:
            cublas_axpy(_cublas_handle, _nloc, &alpha,
                x._data + x._local_offset, 1, _data + _local_offset, 1);
            break;
        }
    }

    void print(std::ostream& out) const
    {
        for (int rank = 0; rank < _nprocs; rank++)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            if (_myrank != rank)
            {
                // schnelle loesung zum geordneten ausgeben
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                continue;
            }
            out << "PROCESS " << rank << ":" << std::endl;
            out << "Prev-Buffer: " << std::endl;
            if (_iam_first) out << "EMPTY" << std::endl;
            else
            {
                for (int i = 0; i < _B; i++)
                    out << "(" << i << "):\t" << _data[i] << std::endl;
            }
            out << std::endl << "Main-Buffer: " << std::endl;
            for (int i = (_iam_first ? 0 : _B); i < (_iam_last ? _length : (_length - _B)); i++)
                out << "(" << i << "):\t" << _data[i] << std::endl;
            out << std::endl << "Next-Buffer: " << std::endl;
            if (_iam_last) out << "EMPTY" << std::endl;
            else
            {
                for (int i = (_iam_first ? _nloc : (_nloc + _B)); i < _length; i++)
                    out << "(" << i << "):\t" << _data[i] << std::endl;
            }
            out << std::endl;
        }
    }

    void swap(BVector<Scalar>& other)
    {
        using std::swap;
        swap(_data, other._data);
        swap(_n, other._n);
        swap(_m, other._m);
        swap(_B, other._B);
        swap(_l, other._l);
        swap(_p, other._p);
        swap(_myrank, other._myrank);
        swap(_offset, other._offset);
        swap(_iam_first, other._iam_first);
        swap(_iam_last, other._iam_last);
        swap(_arch, other._arch);
        swap(_nloc, other._nloc);
        swap(_length, other._length);
        swap(_mygroup, other._mygroup);
        swap(_win_prev, other._win_prev);
        swap(_win_next, other._win_next);
        swap(_buf_prev, other._buf_prev);
        swap(_buf_next, other._buf_next);
        swap(_mygroup, other._mygroup);
        swap(_local_offset, other._local_offset);
    }

    BVector(const BVector& other) = delete;
    BVector& operator=(const BVector& other) = delete;
};

template <class Scalar>
void swap(BVector<Scalar>& v1, BVector<Scalar>& v2)
{
    v1.swap(v2);
}

struct CsrIndexPair
{
    int i, j;

    // Zeilenordnung
    bool operator<(const CsrIndexPair& other) const
    {
        if (i != other.i)
            return i < other.i;
        return j < other.j;
    }
};

template<class Scalar>
class BCsrMatrix;

template<class Scalar>
class RawSparseMatrix
{
    friend class BCsrMatrix<Scalar>;

    std::map<CsrIndexPair, Scalar> _local_data;
    int _n, _B, _l, _p;
    int _myrank, _nprocs;
    int _offset;

public:
    RawSparseMatrix(int n, int B) : _n(n), _B(B)
    {
        //assert(_B > 0);
        MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &_myrank);
        distribute_among_nodes(_n, _l, _p);

        // prüfe, ob lokal überhaupt arbeit getan werden kann
        // (sonst besteht der algorithmus nur aus kommunikation)
        if (_l < 2 * _B || _p < _B)
            throw std::runtime_error("Matrix too small for the specified #processes.");
        assert(_l >= 2 * _B);
        assert(_p >= _B);

        if (_myrank == 0)
            _offset = 0;
        else _offset = _myrank*_l - _B;
    }

    void insert(int i, int j, const Scalar& v)
    {
        assert(i < _n);
        assert(j < _n);
        if (i >= _l*_myrank && i < _l*(_myrank + 1))
            _local_data.insert({ { i - _l*_myrank,j - _offset}, v });
    }

    int erase(int i, int j)
    {
        return _local_data.erase({ i - _l*_myrank, j - _offset });
    }

    void clear()
    {
        _local_data.clear();
    }
};

template<class Scalar>
class BCsrMatrix
{
    Scalar *_val;
    int *_row_ptr, *_col_ind;

    int _n, _m, _B, _l, _p;
    int _myrank, _nprocs;
    int _offset;
    bool _iam_first, _iam_last;

    arch_t _arch;
    int _nloc;
    int _localoffset;
    int _length;

    cusparseHandle_t _cusp_handle;
    cusparseMatDescr_t _desc;

public:
    int get_n() const { return _n; }
    arch_t get_arch() const { return _arch; }
    int get_B() const { return _B; }
    
    BCsrMatrix(RawSparseMatrix<Scalar>& mat, arch_t arch, cusparseHandle_t cusp_handle) :
        _val(nullptr), _row_ptr(nullptr), _col_ind(nullptr),
        _n(mat._n), _B(mat._B), _l(mat._l), _p(mat._p),
        _myrank(mat._myrank),
        _nprocs(mat._nprocs),
        _offset(mat._offset),
        _iam_first(_myrank == 0),
        _iam_last(_myrank == _nprocs - 1),
        _arch(arch),
        _nloc(_iam_last ? _p : _l),
        _localoffset(_B*(!_iam_first)),
        _cusp_handle(cusp_handle)
    {
        if (_nprocs == 1)
            _length = _l;
        else if (_iam_first)
            _length = _l + _B;
        else if (_iam_last)
            _length = _p + _B;
        else
            _length = _l + 2 * _B;
        
        std::vector<int> epl(_nloc); // Anzahl Einträge pro Zeile
        for (auto& p : mat._local_data)
              epl[p.first.i]++;
        _m = *std::max_element(epl.begin(), epl.end());

        switch (_arch)
        {
        case ARCH_CPU:
            _val = new Scalar[_nloc * _m];
            _row_ptr = new int[_nloc + 1];
            _col_ind = new int[_nloc * _m];
            break;
        case ARCH_GPU:
            gpuMalloc(&_val, sizeof(Scalar)*_nloc*_m);
            gpuMalloc(&_row_ptr, sizeof(Scalar)*(_nloc+1));
            gpuMalloc(&_col_ind, sizeof(Scalar)*_nloc*_m);

            cusparseCreateMatDescr(&_desc);
            break;
        }
        memset(_val, 0, sizeof(Scalar)*_nloc*_m);
        // Zugriffe auf Fuelleinträge haben Wert Null und
        // sollten lesend in gültigen Speicher (einfach erster Eintrag) 
        // (und nicht irgendwohin) gehen
        memset(_col_ind, 0, sizeof(int)*_nloc*_m);

        // Fülle row_ptr
        for (int i = 0; i <= _nloc; i++)
            _row_ptr[i] = i*_m;

        // Fülle val und colctr
        int colctr = 0, iprev = 0;
        for (auto p : mat._local_data)
        {
            if (p.first.i != iprev)
            {
                iprev = p.first.i;
                colctr = 0;
            }
            assert(colctr < _m);
            _val[idx2(p.first.i, colctr)] = p.second;
            _col_ind[idx2(p.first.i, colctr)] = p.first.j;
            colctr++;
        }
    }

    ~BCsrMatrix()
    {
        switch (_arch)
        {
        case ARCH_CPU:
            if (_val) delete[] _val;
            if (_row_ptr) delete[] _row_ptr;
            if (_col_ind) delete[] _col_ind;
            break;
        case ARCH_GPU:
            if(_val) cudaFree(_val);
            if(_row_ptr) cudaFree(_row_ptr);
            if(_col_ind) cudaFree(_col_ind);

            cusparseDestroyMatDescr(_desc);
            break;
        }
    }

    int get_myrank() const { return _myrank; }

    int idx2(int i, int j) const
    {
        assert(i < _nloc);
        assert(j < _m);
        return _m*i + j; 
    };

    void print(std::ostream& out) const
    {
        for (int rank = 0; rank < _nprocs; rank++)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            if (_myrank != rank)
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                continue;
            }
            out << "PROCESS " << rank << ":" << std::endl;
            out << "rowptr: " << std::endl;
            for (int i = 0; i < _nloc+1; i++)
                out << "(" << i << "):\t" << _row_ptr[i] << std::endl;
          
            out << std::endl << "col_ind: " << std::endl;
            for (int i = 0; i < _nloc; i++)
                for (int j = 0; j < _m; j++)
                    out << "(" << i << "," << j <<"):\t" << _col_ind[idx2(i,j)] << std::endl;

            out << std::endl << "vals: " << std::endl;
            for (int i = 0; i < _nloc; i++)
                for (int j = 0; j < _m; j++)
                    out << "(" << i << "," << j << "):\t" << _val[idx2(i, j)] << std::endl;

            out << std::endl;
            out.flush();
        }
    }

    void spmv(BVector<Scalar>& src, BVector<Scalar>& dst) const
    {
        //std::cout << "startmv: " << std::chrono::system_clock::now().time_since_epoch().count() << std::endl;
        src.start_comm();
        
        // lokale elemente
        switch (_arch)
        {
        case ARCH_CPU:
        {
            if (_nprocs == 1)
                compute_dst(0, _nloc, src, dst);
          
            else if (_iam_first)
                compute_dst(0, _nloc - _B, src, dst);
            
            else if (_iam_last)
                compute_dst(_B, _nloc, src, dst);
           
            else // ich bin weder erster noch letzter
                compute_dst(_B, _nloc - _B, src, dst);
            
            break;
        }
        case ARCH_GPU:
        {
            // ACHTUNG: Der rowptr bleibt bzgl. der verschobenen
            // vals, colindcs nullbasiert, d.h. dieser Trick
            // funktioniert nur fuer festes _m
            const Scalar one = 1.0, zero = 0.0;
            if (_nprocs == 1)
            {
                cusparse_csrmv(_cusp_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                _nloc, _length, _nloc * _m, &one,
                _desc,
                _val,
                _row_ptr, _col_ind,
                src.data(), &zero,
                dst.local_data());
            }
            else if (_iam_first)
            {
                cusparse_csrmv(_cusp_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                _nloc - _B, _length, (_nloc - _B) * _m, &one,
                _desc,
                _val,
                _row_ptr, _col_ind,
                src.data(), &zero,
                dst.local_data());
            }
            else if (_iam_last)
            {
                cusparse_csrmv(_cusp_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                _nloc - _B, _length, (_nloc - _B) * _m, &one,
                _desc,
                _val + _B*_m,
                _row_ptr, _col_ind + _B*_m,
                src.data(), &zero,
                dst.data() + 2*_B);
            }
            else
            {
                cusparse_csrmv(_cusp_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                _nloc - 2 * _B, _length, (_nloc - 2 * _B) * _m, &one,
                _desc,
                _val + _B*_m,
                _row_ptr, _col_ind + _B*_m,
                src.data(), &zero,
                dst.local_data() + _B);
            }
            break;
        }
        }
        src.finish_comm();
        
        // restliche, neu eingetroffene elemente (immer mit CPU)  
        if (!_iam_first)
            compute_dst(0, _B, src, dst);
        
        if (!_iam_last)
            compute_dst(_nloc - _B, _nloc, src, dst);
        
        //std::cout << "endmv: " << std::chrono::system_clock::now().time_since_epoch().count() << std::endl;
    }

    void compute_dst(const int begin, const int end, const BVector<Scalar>& src, BVector<Scalar>& dst) const
    {
        Scalar* dstptr = dst.local_data();
#pragma omp parallel for
        for (int i = begin; i < end; i++)
        {
            dstptr[i] = 0;
            for (int j = 0; j < _m; j++)
                dstptr[i] += _val[idx2(i, j)] * src[_col_ind[idx2(i, j)]];
        }
    }

    int get_nloc() const { return _nloc; }
};

template <class Scalar>
BCsrMatrix<Scalar> construct_model_matrix(int m, arch_t arch, cusparseHandle_t cusp)
{
    RawSparseMatrix<Scalar> raw(m*m, m);

    raw.insert(0, 0, 4);
    raw.insert(0, 1, -1);
    raw.insert(0, m, -1);
    // 4, -1 [...] -1
    for (int i = 1; i < m; i++)
    {
        raw.insert(i, i - 1, -1);
        raw.insert(i, i, 4);
        raw.insert(i, i + 1, -1);
        raw.insert(i, i + m, -1);
        // -1, 4,-1 [...] - 1
    }
    for (int i = m*(m - 1); i < m*m - 1; i++)
    {
        raw.insert(i, i - m, -1);
        raw.insert(i, i - 1, -1);
        raw.insert(i, i, 4);
        raw.insert(i, i + 1, -1);
        // -1 [...] -1 4 -1 
    }

    for (int i = m; i < m*(m - 1); i++)
    {
        raw.insert(i, i - m, -1);
        raw.insert(i, i - 1, -1);
        raw.insert(i, i, 4);
        raw.insert(i, i + 1, -1);
        raw.insert(i, i + m, -1);
        // -1 [...] -1 4 -1 [...] -1
    }

    raw.insert(m*m - 1, m*m - 1 - m, -1);
    raw.insert(m*m - 1, m*m - 1 - 1, -1);
    raw.insert(m*m - 1, m*m - 1, 4);
    // -1 [...] -1 4

    // Entferne linke, rechte Nachbarn am Rand
    for (int i = 0; i < m; i++)
    {
        raw.erase(i*m, i*m - 1);
        raw.erase((i + 1)*m - 1, (i + 1)*m);
    }
    return BCsrMatrix<Scalar>(raw, arch, cusp);
}

template <class Scalar>
void cg_solve(const BCsrMatrix<Scalar>& mat, 
    const BVector<Scalar>& b, BVector<Scalar>& x0, 
    Scalar TOL = 1e-9, unsigned MAX_ITER = 1000000)
{
    assert(mat.get_nloc() == b.get_nloc());
    assert(mat.get_nloc() == x0.get_nloc());

    int n = mat.get_n();
    int B = mat.get_B();
    arch_t arch = mat.get_arch();

    BVector<Scalar> x1(n,B,arch, b.get_cublas_handle()), z(n, B, arch, b.get_cublas_handle()), 
        d0(n, B, arch, b.get_cublas_handle()), d1(n, B, arch, b.get_cublas_handle()),
        r0(n, B, arch, b.get_cublas_handle()), r1(n, B, arch, b.get_cublas_handle());
    Scalar alpha, beta;
    
    mat.spmv(x0, z);
    b.copy(r0);
    r0.axpy(-1.0, z);

    r0.copy(d0);

    for (unsigned i = 0; i < MAX_ITER; i++)
    {
       mat.spmv(d0, z);
        
        alpha = r0.dot(r0) / d0.dot(z);
       
        x0.copy(x1);
        x1.axpy(alpha, d0);

        r0.copy(r1);
        r1.axpy(-alpha, z);

        Scalar res = r1.nrm2();
#       ifndef NDEBUG
        if(mat.get_myrank() == 0)
            std::cout << "Iteration " << i + 1 << ": res = " << res << std::endl;
#       endif
        if (res < TOL)
        {
            x0.copy(x1);
            if(mat.get_myrank() == 0)
                std::cout << "CG converged after " << i + 1 << " iterations." << std::endl;
            break;
        }

        beta = r1.dot(r1) / r0.dot(r0);

        r1.copy(d1);
        d1.axpy(beta, d0);

        // k->k+1
        swap(x0, x1);
        swap(r0, r1);
        swap(d0, d1);
    }
}

#endif
