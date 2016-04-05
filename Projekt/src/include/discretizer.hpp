#ifndef __DISCRETIZER_HPP_
#define __DISCRETIZER_HPP_

#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include <cassert>

/// \file   discretizer.hpp

namespace Icarus
{

/// \brief  Simple Struktur um einen Raumpunkt kompakt zu speichern.
///         Kann ebenso dazu benutzt werden einen Vektor zu speichern.
struct Vertex {
    float x, y, z;
};

/// \brief  Eine Flaeche im Raum, aufgespannt von beliebig vielen Punkten.
///         Der Benutzer muss selbst darauf achten, dass die Eckpunkte alle in einer Ebene liegen.
class Face
{
public:
    /// \brief  Konstruktor
    /// \param  num_vertices    Anzahl Eckpunkte, die die Flaeche aufspannen
    /// \param  id_vertices     ID's geben an welche Punkte von 'vertices' zu der Flaeche gehoeren.
    /// \param  vertices        Liste von Punkten, die mindestens die Eckpunkte der Flaeche enthaelt.
    /// \param  normal          Normalenvektor zur Flaeche.
    Face(int num_vertices, std::vector<int>& id_vertices,
        std::vector<Vertex>& vertices, Vertex normal);

    /// \brief  Ueberprueft, ob ein Punkt in der Projektion der Flaeche auf die Y-Z-Ebene liegt.
    /// \param  point   Raumpunkt der ueberprueft werden soll.
    /// \return Gibt zurueck, ob der Punkt in der Projektion der Flaeche auf die Y-Z-Ebene liegt oder nicht.
    bool pointInsideYz(Vertex point);

    /// \brief  Gibt einen Eckpunkt der Fleache zureuck.
    /// \param  local_id    ID des gewuenschten Punktes bzgl der Nummerierung der Punkte der Flaeche.
    /// \return ID-ter Eckpunkt der Flaeche.
    Vertex get_vertex(int local_id) { assert(local_id < _num_vertices); return _vertices[local_id]; }

    /// \brief  Gibt Normalenvektor zurueck.
    /// \return Normalenvektor der Flaeche.
    Vertex get_normal() { return _normal; }

private:
    int _num_vertices;
    std::vector<Vertex> _vertices;
    Vertex _normal;
    Vertex _last_point;
    bool _last_check;
};

/// \brief  Ein Objekt im Raum.
///         Wird durch Flaechen aufgespannt. Ermoeglicht die Ueberpruefung,
///         ob ein Raumpunkt im Inneren, auf dem Rand oder ausserhalb des Objektes liegt.
class Object
{
public:
    /// \brief  Konstruktor
    /// \param  name    Optionaler Name des Objektes um daraus zB Material o.ae. abzuleiten.
    Object(std::string name = "");

    // TODO TODISCUSS: set_vertex oder add_vertex? Hier ist ein Hinzufuegen gemeint.
    // ebenso weiter unten
    /// \brief  Fuegt einen Raumpunkt hinzu.
    /// \param  x   x-Koordinate
    /// \param  y   y-Koordinate
    /// \param  z   z-Koordinate
    void set_vertex(float x, float y, float z);

    /// \brief  Fuegt einen Normalenvektor hinzu, der spaeter einer Flaeche zugewiesen werden kann.
    /// \param  x   x-Koordinate
    /// \param  y   y-Koordinate
    /// \param  z   z-Koordinate
    void set_normal(float x, float y, float z);

    /// \brief  Fuegt eine Flaeche hinzu.
    ///         Benutzt werden dafuer bereits gespeicherte Punkte und Normalenvektoren.
    ///         Es wird nicht ueberprueft, ob alle Punkte in einer Ebene liegen.
    /// \param  num_vertices    Anzahl an Eckpunkten der Flaeche.
    /// \param  id_vertices     Liste von ID's bzgl aller Punkte des Objektes, die die Flaeche aufspannen.
    /// \param  local_id_normal ID bzgl aller Normalenvektoren des Objektes, der der Flaeche zugeordnet werden soll.
    void set_face(int num_vertices, std::vector<int>& id_vertices, int local_id_normal);

    /// \brief  Ueberprueft, ob Punkt in, auf oder ausserhalb vom Objekt liegt.
    /// \param  point   Zu ueberpruefender Punkt.
    /// \return Gibt 'o' zurück, wenn sich der Punkt im Objekt befindet, sonst 'a'.
    char pointInside(Vertex point);

private:
    std::string _name;
    int _num_vertices;
    int _num_faces;
    std::vector<Vertex> _vertices;
    std::vector<Vertex> _normals;
    std::vector<Face> _faces;
};

/// \brief  "Diskretisiert" Raum anhand einer obj-Datei.
///         Die Diskretisierung beginnt im Ursprung und geht in positve Koordinatenrichtungen.
/// \param  inputFile   Relativer Pfad zur Input-Datei.
/// \param  h           Aequidistante Schrittweite der Diskretisierung (h=h_x=h_y=h_z).
/// \param  nx          Anzahl an Diskretisierungen in x-Richtung.
/// \param  ny          Anzahl an Diskretisierungen in y-Richtung.
/// \param  nz          Anzahl an Diskretisierungen in z-Richtung.
/// \return Vektor der jeden Gitterpunkt mit der Indizierung i=x+nx*y+nx*ny*z enthaelt.
///         'a' -> frei, 'o' -> im Innern eines Objektes, 'b' -> Am Rand eines Objektes oder des gesamten Gebietes.
std::vector<char> discretizer(std::string inputFile,
    float h, int nx, int ny, int nz);

// eventuelles TODO: Moeglichkeit einbauen nicht vom Ursprung zu diskretisieren. 
// zB erkennen, von wo bis wo Objekte liegen und ein Padding einbauen.
// Oder einfach weitere Eingabeparameter, die Eckpunkte angeben.

/// \brief  "Diskretiesiert" mit 'discretizer' und schreibt den Rueckgabevektor in eine Datei.
/// \param  outputFile  Relativer Pfad, in den die "Diskretisierung" geschrieben wird.
/* \ s a     #discretizer */ //TODO TOCHECK klappt sowas?
void save_discretizer(std::vector<char> discretized_points,
                      std::string outputFile,
                      int nx, int ny, int nz);

}//namespace Icarus

#endif//__DISCRETIZER_HPP_
