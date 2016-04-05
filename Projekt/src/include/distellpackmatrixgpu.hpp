#ifndef __DISTELLPACKMATRIXGPU_HPP_
#define __DISTELLPACKMATRIXGPU_HPP_

#include "matrix.hpp"
#include "utility.hpp"
#include "mpihandler.hpp"
#include "fullvectorgpu.hpp"
#include <fstream>
#include <iterator>
#include <algorithm>


namespace Icarus
{

/**
  * \brief  Dünnbesetzte, quadratische Matrix, deren Zeilen gleichverteilt auf einer Menge von Nodes liegen.
  *
  * Die Zeilen dieser Matrix liegen (annähernd) gleichverteilt auf einer Menge
  * von Nodes der zugeordneten Prozessgruppe. Lokal auf der Node werden die Zeilen
  * im Ellpack-Format gespeichert für maximale Effizienz der MV-Multiplikation
  * in CUDA.
  *
  * Dieser Matrixtyp kann nur sequentiell zeilenweise gefüllt werden, wobei
  * die maximale Zeilenlänge (node-weise) bekannt sein muss. Siehe dazu auch
  * die Dokumtation der Funktion sequential_fill.
  *
  * \tparam   Scalar  Skalarer Typ der Einträge.
  */
template <typename Scalar>
class DistEllpackMatrixGpu : public Matrix<DistEllpackMatrixGpu<Scalar>>
{
    friend class Matrix<DistEllpackMatrixGpu<Scalar>>;

    // MPI Eigenschaften
    MPI_Comm _my_comm;
    int _my_rank, _num_nodes;

    // Mit PAD wird das padding durchgeführt
    static const int PAD = 0;

    size_t _dim_global, _dim_local, _dim_local_nopad, _max_row_length;

    size_t * _indices;

    Scalar* _data;

    // hilfsvariablen zum sequentiellen füllen der matrix
    size_t _row_ptr, _col_ptr;
    bool _filled;

public:

    /// Zugeordneter (d.h. bezüglich der Operatoren vertäglicher) Vektor-Typ.
    typedef typename MatrixTraits<DistEllpackMatrixGpu<Scalar>>::VectorType VectorType;

    /**
     * \brief   Standardkonstruktor.
     *
     * Erzeugt einen Vektor der Dimension dim, der komplett auf jeder Node der
     * Prozessgruppe my_comm liegt.
     *
     * \param   dim_globaö     Dimension der Matrix.
     * \param   my_comm         Kommunikator in die Prozessgruppe der Matrix.
     */
    DistEllpackMatrixGpu(size_t dim_global, MPI_Comm my_comm = MPI_COMM_WORLD);

    ~DistEllpackMatrixGpu();

    DistEllpackMatrixGpu(DistEllpackMatrixGpu&& other);

    DistEllpackMatrixGpu(const DistEllpackMatrixGpu& other);

    DistEllpackMatrixGpu& operator=(DistEllpackMatrixGpu&& other);

    DistEllpackMatrixGpu& operator=(const DistEllpackMatrixGpu& other);

     /**
     * \brief   Gibt den Kommunikator in die Prozessgruppe der Matrix zurück.
     */
    MPI_Comm get_comm() const { return _my_comm; }

    /**
      * \brief   Gibt die lokale Dimension der Matrix, d.h. die Anzahl
      *          der auf der aufrufenden Node gespeicherten Zeilen zurück.
      */
    size_t get_dim_local() const { return _dim_local; }

    /**
      * \brief   Gibt die Anzahl der Zeilen, wie sie auf einer der ersten N-1 Nodes liegen, zurück.
      */
    size_t get_dim_local_nopad() const { return _dim_local_nopad; }

    /**
      * \brief   Gibt die globale Dimension, d.h. die Anzahl der Zeilen der Matrix, zurück.
      */
    size_t get_dim_global() const { return _dim_global; }

    //TODO: get, set local/global

    /**
      * \brief   Bereitet den zeilenweisen Füllvorgang der Matrix vor.
      *
      * Diese Funktion muss von jeder Node genau einmal zu Beginn des Füllvorgangs aufgerufen werden.
      * Anschließend können die Zeilen mit sequential_fill und end_of_row gefüllt werden, beginnend
      * bei der lokal ersten Zeile.
      * Die maximale auf dieser Node auftretende Zeilenlänge muss vorher bekannt sein.
      *
      * \param max_row_length Maximal auf dieser Node auftretende Zeilenlänge.
      */
    void prepare_sequential_fill(size_t max_row_length);

    /**
      * \brief   Fülle die Matrix zeilenweise.
      *
      * Diese Funktion fügt der aktuellen Zeile den Wert val mit Spaltenindex col hinzu.
      * Vor dem erste Aufruf dieser Funktion muss mit prepare_sequential_fill die maximal
      * auf der füllenden Node auftretende Zeilenlänge gesetzt werden.
      * Nachdem der letzte Eintrag einer Zeile gesetzt wurde, wird mit end_of_row
      * die Zeile beendet.
      *
      * \param colind Spaltenindex des einzutregenden Werts. Es muss colind < dim_global gelten.
      * \param val    Wert, der an die Position colind geschrieben werden soll.
      */
    void sequential_fill(size_t colind, const Scalar& val);

    /**
      * \brief   Beende eine Zeile beim zeilenweisen Füllen der Matrix.
      *
      * Diese Funktion beendet die aktuelle Zeile beim Füllvorgang und setzt den
      * Füllcursor auf die nächste Zeile, oder beendet den Füllvorgang, falls die
      * letzte auf der Node vorhandene Zeile beendet wurde.
      *
      * Vor dem erste Aufruf dieser Funktion muss mit prepare_sequential_fill die maximal
      * auf der füllenden Node auftretende Zeilenlänge gesetzt werden.
      *
      * Diese Funktion muss während eines Füllvorgangs auf der Node genau dim_local mal
      * aufgerufen werden.
      */
    void end_of_row();

    /**
      * \brief   Prüft, ob die Matrix korrekt gefüllt wurde.
      *
      * Ein positiver Rückgabewert dieser Funktion ist einerseits ein Indikator für einen erfolgreich
      * abgeschlossenen Füllvorgang und andererseits die Voraussetzung für sämtliche algebraische
      * Operationen mit der Matrix.
      *
      * \return Gibt zurück, ob die Matrix korrekt gefüllt wurde.
      */
    bool is_filled() const { return _filled; }

    /**
      * \brief   Gibt den globalen Index der ersten auf der Node liegenden Zeile zurück.
      */
    size_t first_row_on_node() const { return _my_rank * _dim_local_nopad; }

    /**
      * \brief   Erstellt einen zu der Matrix passenden Äquilibrierungsvorkonditionierer.
      *
      * \return  Der Vorkonditionierer hat denselben Typ wie das Objekt, auf das die Funktion
      *          aufgerufen wird.
      */
    DistEllpackMatrixGpu precond_equi() const;

    /**
      * \brief   Erstellt einen zu der Matrix passenden Äquilibrierungsvorkonditionierer.
      *
      * Wenn in einer Zeile eine Null auf der Diagonalen steht, wird diese Zeile
      * durch die Vorkonditionierung nicht verändert.
      *
      * \return  Der Vorkonditionierer hat denselben Typ wie das Objekt, auf das die Funktion
      *          aufgerufen wird.
      */
    DistEllpackMatrixGpu precond_jacobi() const;

    /**
      * \brief   Schreibe den lokalen Inhalt des Block in den Stream out.
      *
      * Für die Verwendung dieser Funktion muss eine entsprechende Überladung des
      * Operators std::ostream::operator<<(Scalar) existieren.
      *
      * \param out  Stream, in den die Ausgabe geschrieben werden soll.
      */
    void print_local_data(std::ostream &os) const;

     /**
      * \brief   Lese eine DistEllpackMatrixGpu aus einem CSR-artigen Dateiformat ein.
      *
      * Das benötigte Dateiformat wird von dem MATLAB-Skript /util/csrwrite.m
      * erzeugt. Die Informationen werden in drei Dateien gespeichert, deren Namen
      * aus einem gemeinsamen Präfix und verschiedenen Endungen bestehen.
      * Dieser Funktion wird (wie auch csrwrite.m) dieser Präfix übergeben.
      *
      * \param filename  Präfix des Dateinamens des Dateitripels, das eingelesen werden soll.
      * \param new_comm  Kommunikator in die Prozessgruppe, der die neu erzeugte Matrix
      *                  gehören soll.
      *
      * \return Gibt die aus dem Dateitripel erzeugte DistEllpackMatrixGpu zurück.
      */
    static DistEllpackMatrixGpu import_csr_file(const std::string& filename, MPI_Comm new_comm = MPI_COMM_WORLD);

private:

    size_t get_dim_impl() const {return _dim_global;}

    void mult_vec_impl(const VectorType& vec, VectorType& result) const;
};

}//namespace Icarus

#include "distellpackmatrixgpu.tpp"

#endif // __DISTELLPACKMATRIXGPU_HPP_
