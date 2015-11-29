/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                vector.hpp
* Erstellt:                 24.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Definiert das Interface Vektor, das alle Operationen vorstellt,
*   die ein Vektor implementieren muss.
* - Die "ist ein" Beziehung wird durch das CRTP, also als statischer
*   Polymorphismus realisiert.
*/

#ifndef __VECTOR_HPP_
#define __VECTOR_HPP_

namespace Icarus
{

template<typename T>
struct VectorTraits;

// das interface für alle vektoren
template<class Child>
class Vector
{
public:

    typedef typename VectorTraits<Child>::RealType RealType;
    typedef typename VectorTraits<Child>::ScalarType ScalarType;

protected:
    Vector() { }

public:

    Child& leaf() { return static_cast<Child&>(*this); }
    const Child& leaf() const { return static_cast<const Child&>(*this); }

    RealType l2norm2() const { return leaf().l2norm2_impl(); }

    RealType l2norm() const { return sqrt(leaf().l2norm2_impl()); }

    RealType maxnorm() const { return leaf().maxnorm_impl(); }

    void clear() { leaf().clear_impl(); }

    void fill_const(const ScalarType& s) { leaf().fill_const_impl(s); }

    ScalarType scal_prod(const Vector<Child>& other) { return leaf().scal_prod_impl(other.leaf()); }

    void axpy(const ScalarType& alpha, const Vector<Child>& y) { leaf().axpy_impl(alpha, y.leaf()); }

    void scal(const ScalarType& alpha) { leaf().scal_impl(alpha); }

    void copy(const Vector<Child>& other) {leaf().copy_impl(other.leaf());}

    size_t get_dim() const { leaf().get_dim_impl(); }

    void swap(Vector<Child>& other) { leaf().swap_impl(other.leaf()); }

};

}

#endif // __VECTOR_HPP_
