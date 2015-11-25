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
    Child& _leaf;
public:

    typedef typename VectorTraits<Child>::RealType RealType;
    typedef typename VectorTraits<Child>::ScalarType ScalarType;

    Vector() : _leaf(*static_cast<Child*>(this)) { }

    RealType l2norm2() const
    {
        return _leaf.l2norm2_impl();
    }

    RealType l2norm() const
    {
        return sqrt(_leaf.l2norm2_impl());
    }

    RealType maxnorm() const
    {
        return _leaf.maxnorm_impl();
    }

    void clear()
    {
        _leaf.clear_impl();
    }

    void fill_const(const ScalarType& s)
    {
        _leaf.fill_const_impl(s);
    }

};

}

#endif // __VECTOR_HPP_
