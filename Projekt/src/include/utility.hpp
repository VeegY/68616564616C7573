#ifndef __UTILITY_HPP_
#define __UTILITY_HPP_

#include <fstream>
#include <limits>
#include <algorithm>
#include <iterator>

namespace Icarus
{

/// \brief Setze Lesepointer auf Zeile line. Die Zeilennummern
/// werden nullbasiert gezählt.
/// \param f    Dateistream, auf dem operiert wird
/// \param line Zeilennummer, zu der gesprungen wird
std::istream& go_to_line(std::istream& f, size_t line)
{
    f.seekg(0);
    for(size_t i=0; i<line; i++)
        f.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    return f;
}

/// \brief Zählt die Zeilen (definiert durch '\n') in der Datei f
/// \param f    Dateistream, auf dem operiert wird
size_t get_num_lines(std::ifstream& f)
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
std::istream& skip_lines(std::istream& f, size_t num_lines)
{
    for(size_t i=0; i<num_lines; i++)
        f.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    return f;
}

}

#endif // __UTILITY_HPP_
