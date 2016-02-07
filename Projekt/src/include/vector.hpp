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

/**
  * \brief Basisklasse (Interface) für alle Vektortypen
  *
  * Diese Klasse definiert das Interface Vector, das alle Operationen vorstellt,
  * die ein Vektor implementieren muss. Die entsprechenden Implementierungen
  * werden als private Methoden der abgeleiteten Klassen mit dem Suffix "_impl"
  * bereitgestellt.
  *
  * \tparam   Child   Der Typ der abgeleiteten Klasse (siehe CRTP).   
  */
template<class Child>
class Vector
{
public:

    typedef typename VectorTraits<Child>::RealType RealType;
    typedef typename VectorTraits<Child>::ScalarType ScalarType;

protected:
    Vector() { }

public:

    /**
      * \brief   Zugriff auf den abgeleiteten Typ.
      */ 
    Child& leaf() { return static_cast<Child&>(*this); }
    
    /**
      * \brief   Zugriff auf den abgeleiteten Typ, konstante Variante.
      */ 
    const Child& leaf() const { return static_cast<const Child&>(*this); }

    /**
      * \brief   Berechnet das Quadrat der L2-Norm des Vektors.
      */ 
    RealType l2norm2() const { return leaf().l2norm2_impl(); }

    /**
      * \brief   Berechnet die L2-Norm des Vektors.
      */ 
    RealType l2norm() const { return sqrt(leaf().l2norm2_impl()); }

    /**
      * \brief   Berechnet die Maximum-Norm des Vektors.
      */ 
    RealType maxnorm() const { return leaf().maxnorm_impl(); }

    /**
      * \brief   Füllt den gesamten Vektor mit Scalar(0).
      */ 
    void clear() { leaf().clear_impl(); }

    /**
      * \brief   Füllt den gesamten Vektor mit Kopien von s.
      *
      * \param  s   Objekt, das an alle Positionen des Vektors geschrieben werden soll.
      */ 
    void fill_const(const ScalarType& s) { leaf().fill_const_impl(s); }

    /**
      * \brief   Berechne das Skalarprodukt mit einem zweiten Vektor desselben Typs.
      *
      * Bei komplexen Datentypen wird der zweite Faktor automatisch konjugiert.
      *
      * \param  other   Der zweite Vektor in dem Skalarprodukt.
      */ 
    ScalarType scal_prod(const Vector<Child>& other) { return leaf().scal_prod_impl(other.leaf()); }

    /**
      * \brief   Berechne die BLAS-Operation x <- x + alpha y.
      *
      * \param  alpha   Der Skalar, mit dem der zweite Vektor multipliziert wird.
      * \param  y       Der zweite Vektor.
      */ 
    void axpy(const ScalarType& alpha, const Vector<Child>& y) { leaf().axpy_impl(alpha, y.leaf()); }

    /**
      * \brief   Skaliere den Vektor um einen konstanten Skalar.
      *
      * Multipliziere alle Komponenten des Vektors mit alpha.
      *
      * \param  alpha   Der Skalar, mit dem der Vektor multipliziert werden.
      */ 
    void scal(const ScalarType& alpha) { leaf().scal_impl(alpha); }

    /**
      * \brief   Gibt die Dimension des Vektors zurück.
      */ 
    size_t get_dim() const { return leaf().get_dim_impl(); }

    /**
      * \brief   Explizite Kopieroperation.
      *
      * Kopiert den Inhalt eines zweiten Vektors in den Vektor.
      *
      * \param   other   Vektor, dessen Inhalt in diesen Vektor kopiert wird.   
      */
    void copy(const Vector<Child>& other) {leaf().copy_impl(other.leaf());}

    /**
      * \brief   Explizite Tauschoperation.
      *
      * Tauscht den Inhalt eines zweiten Vektors mit diesem Vektor.
      *
      * \param   other   Vektor, dessen Inhalt mit diesem Vektor vertauscht wird.   
      */
    void swap(Vector<Child>& other) { leaf().swap_impl(other.leaf()); }

};

}

#endif // __VECTOR_HPP_
