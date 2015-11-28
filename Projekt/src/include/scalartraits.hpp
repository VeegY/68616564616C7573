/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                scalartraits.hpp
* Erstellt:                 23.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Traits für die skalaren Datentypen.
*/

#ifndef __SCALARTRAITS_HPP_
#define __SCALARTRAITS_HPP_

#include <complex>
#include "mpi.h"

namespace Icarus
{
    template <typename Scalar>
    struct ScalarTraits;

    template<>
    struct ScalarTraits<double>
    {
        static constexpr MPI_Datatype mpi_type = MPI_DOUBLE;
        typedef double RealType;
        static double abs2(double d) { return d*d; }
        static double abs(double d) { return fabs(d); }
        static double smult(double d1, double d2) { return d1*d2; }
    };

    template<>
    struct ScalarTraits<float>
    {
        static constexpr MPI_Datatype mpi_type = MPI_FLOAT;
        typedef float RealType;
        static float abs2(float f) { return f*f; }
        static float abs(float f) { return fabs(f); }
        static float smult(float f1, float f2) { return f1*f2; }
    };

    template<>
    struct ScalarTraits<std::complex<float> >
    {
        static constexpr MPI_Datatype mpi_type = MPI_COMPLEX;
        typedef float RealType;
        static float abs2(const std::complex<float>& c) { return c.real()*c.real() + c.imag()*c.imag(); }
        static double abs(const std::complex<float>& c) { return sqrt(abs2(c)); }
        static std::complex<float> smult(
                const std::complex<float>& c1,
                const std::complex<float>& c2)
        { return c1 * std::conj(c2); }
    };

    template<>
    struct ScalarTraits<std::complex<double> >
    {
        static constexpr MPI_Datatype mpi_type = MPI_DOUBLE_COMPLEX;
        typedef double RealType;
        static double abs2(const std::complex<double>& c) { return c.real()*c.real() + c.imag()*c.imag(); }
        static double abs(const std::complex<float>& c) { return sqrt(abs2(c)); }
        static std::complex<double> smult(
                const std::complex<double>& c1,
                const std::complex<double>& c2)
        { return c1 * std::conj(c2); }

    };
}

#endif // __SCALARTRAITS_HPP_
