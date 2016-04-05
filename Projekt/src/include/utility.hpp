#ifndef __UTILITY_HPP_
#define __UTILITY_HPP_

#include <fstream>
#include <limits>
#include <algorithm>
#include <iterator>
#include <ostream>

#include "mpihandler.hpp"

namespace Icarus
{

/**
* \brief Geordnete Ausgabe eines Objekts, das die print_local_data
*        Methode unterstützt.
* \tparam T Typ des Objekts, das ausgegeben werden soll
* \param  obj   Objekt, das ausgegeben werden soll.
* \param  out   Stream, auf dem das Ergebnis auegegeben werden soll
* \param  comm  Communicator in die Gruppe, auf der das Objekt verteilt liegt
*/
template <typename T>
void print_sliced_object(const T& obj, std::ostream& out = std::cout, MPI_Comm comm = MPI_COMM_WORLD)
{
    int nprocs, myrank;
    MPI_SCALL(MPI_Comm_rank(comm, &myrank));
    MPI_SCALL(MPI_Comm_size(comm, &nprocs));

    MPI_SCALL(MPI_Barrier(comm));

    for (int i = 0; i < nprocs; i++)
    {
        if (myrank == i)
        {
            out << "node " << i << ":" << std::endl;
            obj.print_local_data(out);
        }
        MPI_SCALL(MPI_Barrier(comm));
    }
}

/**
 * \brief Wandle linearen Index in (x,y,z) Tripel um.
 */
inline void deflatten_3d(size_t index,
    size_t nx, size_t ny,
    size_t& ax, size_t& ay, size_t& az)
{
    az = (index / nx) / ny;
    ay = index / nx - az*ny;
    ax = index - ay*nx - az*nx*ny;
}

/// \brief Setze Lesepointer auf Zeile line. Die Zeilennummern
/// werden nullbasiert gezählt.
/// \param f    Dateistream, auf dem operiert wird
/// \param line Zeilennummer, zu der gesprungen wird
inline std::istream& go_to_line(std::istream& f, size_t line)
{
    f.seekg(0);
    for (size_t i = 0; i < line; i++)
        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return f;
}

/// \brief Zählt die Zeilen (definiert durch '\n') in der Datei f
/// \param f    Dateistream, auf dem operiert wird
inline size_t get_num_lines(std::ifstream& f)
{
    f.seekg(0);
    f.unsetf(std::ios_base::skipws);
    size_t res = std::count(
        std::istream_iterator<char>(f),
        std::istream_iterator<char>(),
        '\n');
    f.setf(std::ios_base::skipws);
    f.clear(std::ios_base::eofbit);
    f.seekg(0);
    return res;
}

/// \brief Streamt ein Objekt und setzt den Streampointer
///        anschließend in den ursprünglichen Zustand zurück.
/// \tparam T    Typ des Objekts, in das gestreamt wird
/// \param  f    Dateistream, auf dem operiert wird
/// \param  dst  Zielobjekt, in das gestreamt wird
template<typename T>
void peek_obj(std::istream& f, T& dst)
{
    const std::streampos pos = f.tellg();
    f >> dst;
    f.clear(std::ios_base::eofbit);
    f.seekg(pos);
}

/// \brief Spult den Stream eine Anzahl von Zeilen vor
/// \param f         Dateistream, auf dem operiert wird
/// \param num_lines Anzahl der Zeilen, um die vorgespult wird
inline std::istream& skip_lines(std::istream& f, size_t num_lines)
{
    for (size_t i = 0; i < num_lines; i++)
        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return f;
}

}//namespace Icarus

#endif // __UTILITY_HPP_
